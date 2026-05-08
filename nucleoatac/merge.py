"""
Script to merge nuc positions

@author: Alicia Schep
"""
##### IMPORT MODULES #####
#

from pyatac.utils import shell_command, save_params_json
from pyatac.chunk import Chunk, ChunkList
from pyatac.fragmentsizes import FragmentSizes
import gzip
import numpy as np
import pysam
import os

class MergedNuc(Chunk):
    def __init__(self, chrom, start, end, occ, occ_lower, occ_upper, reads, source):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.occ = occ
        self.occ_lower = occ_lower
        self.occ_upper = occ_upper
        self.reads = reads
        self.source = source
    def asBed(self):
        out = "\t".join(map(str,[self.chrom,self.start,self.end,self.occ,self.occ_lower,self.occ_upper,self.reads, self.source])) 
        return out
    def write(self, handle):
        """write bed line for peak"""
        handle.write(self.asBed() + "\n")

class NucList(ChunkList):
    def __init__(self, *args):
        list.__init__(self,args)
    @staticmethod
    def read(bedfile, source, min_occ = 0):
        """Make a list of chunks from a tab-delimited bedfile"""
        opener = gzip.open if bedfile[-3:] == '.gz' else open
        out = NucList()
        with opener(bedfile,"rt") as infile:
            for line in infile:
                in_line = line.rstrip('\n').split("\t")
                start = int(in_line[1])
                end = int(in_line[2])
                if source == "occ":
                    occ = float(in_line[3])
                    occ_lower = float(in_line[4])
                    occ_upper = float(in_line[5])
                    reads = float(in_line[6])
                elif source == "nuc":
                    occ = float(in_line[4])
                    occ_lower = float(in_line[5])
                    occ_upper = float(in_line[6])
                    reads = float(in_line[10]) + float(in_line[11])
                else:
                    raise Exception("source must be 'occ' or 'nuc'")
                if occ_lower >= min_occ:
                    out.append(MergedNuc(in_line[0],start, end, occ, occ_lower, occ_upper, reads, source))
        return out




def merge(occ_peaks, nuc_calls, sep = 120):
    keep = NucList()
    i = 0 #index for occ peaks
    j = 0 #index for nuc calls
    while i < len(occ_peaks) and j < len(nuc_calls):
        if occ_peaks[i].chrom < nuc_calls[j].chrom:
            keep.append(occ_peaks[i])
            i += 1
        elif occ_peaks[i].chrom > nuc_calls[j].chrom:
            keep.append(nuc_calls[j])
            j += 1
        elif occ_peaks[i].start < (nuc_calls[j].start - sep):
            keep.append(occ_peaks[i])
            i += 1
        elif occ_peaks[i].start > (nuc_calls[j].start + sep):
            keep.append(nuc_calls[j])
            j += 1
        else:
            i += 1
    while j < len(nuc_calls):
        keep.append(nuc_calls[j])
        j += 1
    while i < len(occ_peaks):
        keep.append(occ_peaks[i])
        i += 1
    return keep

def _merge_fragmentsizes(chrom_files, out_file):
    """Merge per-chromosome fragmentsizes.txt files.

    Two valid workflows are auto-detected:

    1. All per-chrom occ runs computed independent fragment-size distributions
       (no --sizes). Each file carries a #raw_counts vector; we sum the raw
       counts and renormalise. Result equals a combined-occ run on the same
       chroms only if peak calls happen to line up; in general nuc_dist may
       drift from a combined-occ run because each per-chrom occ used a
       chrom-specific NFR fit.

    2. All per-chrom occ runs shared a global --sizes input (recommended
       parallel workflow). No file has #raw_counts, and all files are
       byte-identical copies of the input. We copy the first one through.
       This is the only configuration where merged nuc_dist matches a
       combined-occ run exactly.

    Any other combination (mixed presence, or files without #raw_counts that
    differ) is pathological and raises.
    """
    objs = [FragmentSizes.open(f) for f in chrom_files]
    lower, upper = objs[0].lower, objs[0].upper
    for f, o in zip(chrom_files, objs):
        assert (o.lower, o.upper) == (lower, upper), (
            f"FragmentSizes bound mismatch in {f}: "
            f"{o.lower}-{o.upper} vs {lower}-{upper}")
    have_raw = [o.raw is not None for o in objs]
    if all(have_raw):
        total_raw = np.sum([o.raw for o in objs], axis=0).astype(np.int64)
        total = total_raw.sum()
        merged = FragmentSizes(lower, upper)
        merged.raw = total_raw
        merged.vals = total_raw / (total + (total == 0))
        merged.save(out_file)
    elif not any(have_raw) and all(np.allclose(o.vals, objs[0].vals, atol=0, rtol=0)
                                   for o in objs[1:]):
        import shutil
        shutil.copy2(chrom_files[0], out_file)
    else:
        partial = [f for f, h in zip(chrom_files, have_raw) if not h]
        raise RuntimeError(
            "Cannot merge fragmentsizes.txt: per-chrom files are inconsistent. "
            "Either all files must carry the #raw_counts section "
            "(independent per-chrom occ, no --sizes), or all must lack it AND "
            "be byte-identical (per-chrom occ runs with a shared --sizes "
            "input). Files lacking #raw_counts:\n  "
            + "\n  ".join(partial) +
            "\nThis typically indicates a mix of workflows or files generated "
            "by an older NucleoATAC version. Re-run per-chrom occ uniformly "
            "(see the parallelization section of the README).")


def _merge_nuc_dist(chrom_files, out_file):
    """Merge per-chromosome nuc_dist.txt files by summing element-wise.
    Each per-chrom file is itself an unnormalised sum of per-peak normalised
    distributions, so chromosome-level summation reproduces the genome-wide
    file exactly."""
    objs = [FragmentSizes.open(f) for f in chrom_files]
    lower, upper = objs[0].lower, objs[0].upper
    for f, o in zip(chrom_files, objs):
        assert (o.lower, o.upper) == (lower, upper), (
            f"FragmentSizes bound mismatch in {f}: "
            f"{o.lower}-{o.upper} vs {lower}-{upper}")
    summed = np.sum([o.vals for o in objs], axis=0)
    out = FragmentSizes(lower, upper, vals=summed)
    out.save(out_file)


def run_merge_chroms(args):
    """Merge per-chromosome NucleoATAC output files.

    Expects per-chromosome outputs named as <prefix>.<chrom>.<suffix>
    and produces merged outputs named as <out>.<suffix>.
    """
    chroms = args.chroms.split(',')
    prefix = args.prefix
    out = args.out

    # Determine which file suffixes to merge based on what exists
    occ_suffixes = [
        'occ.bedgraph.gz',
        'occ.lower_bound.bedgraph.gz',
        'occ.upper_bound.bedgraph.gz',
        'occpeaks.bed.gz',
        'fragmentsizes.txt',
        'nuc_dist.txt',
    ]
    nuc_suffixes = [
        'nucpos.bed.gz',
        'nucpos.redundant.bed.gz',
        'nucleoatac_signal.bedgraph.gz',
        'nucleoatac_signal.smooth.bedgraph.gz',
        'nucleoatac_background.bedgraph.gz',
        'nucleoatac_raw.bedgraph.gz',
    ]
    all_suffixes = occ_suffixes + nuc_suffixes

    for suffix in all_suffixes:
        # Check if per-chrom files exist for this suffix
        chrom_files = []
        for chrom in chroms:
            f = f"{prefix}.{chrom}.{suffix}"
            if os.path.exists(f):
                chrom_files.append(f)

        if not chrom_files:
            continue

        out_file = f"{out}.{suffix}"

        if suffix.endswith('.gz'):
            # Concatenate compressed bedgraph/bed files:
            # decompress, concatenate, recompress, reindex
            uncompressed = out_file[:-3]  # remove .gz
            with open(uncompressed, 'w') as out_handle:
                for f in chrom_files:
                    with gzip.open(f, 'rt') as in_handle:
                        for line in in_handle:
                            out_handle.write(line)
            pysam.tabix_compress(uncompressed, out_file, force=True)
            os.remove(uncompressed)
            pysam.tabix_index(out_file, preset="bed", force=True)
            print(f"  merged {len(chrom_files)} files -> {out_file}")
        elif suffix == 'fragmentsizes.txt':
            _merge_fragmentsizes(chrom_files, out_file)
            print(f"  merged {len(chrom_files)} files -> {out_file}")
        elif suffix == 'nuc_dist.txt':
            _merge_nuc_dist(chrom_files, out_file)
            print(f"  merged {len(chrom_files)} files -> {out_file}")

    # Copy PDF files if they exist (from first chrom)
    for suffix in ['fragmentsizes.pdf', 'occ_fit.pdf', 'nuc_dist.pdf']:
        first = f"{prefix}.{chroms[0]}.{suffix}"
        if os.path.exists(first):
            import shutil
            shutil.copy2(first, f"{out}.{suffix}")

    print(f"merge-chroms complete: {len(chroms)} chromosomes merged to {out}.*")


def run_merge(args):
    if not args.out:
        args.out = '.'.join(os.path.basename(args.nucpos).split('.')[0:-3])
    occ = NucList.read(args.occpeaks, "occ", args.min_occ)
    nuc = NucList.read(args.nucpos, "nuc", args.min_occ)
    new = merge(occ, nuc, args.sep)
    save_params_json(args.out + '.merge.params.json', 'merge', {
        "occpeaks": args.occpeaks,
        "nucpos": args.nucpos,
        "out": args.out,
        "sep": args.sep,
        "min_occ": args.min_occ,
    })
    bed_path = args.out + '.nucmap_combined.bed'
    with open(bed_path, 'w') as out:
        out.write(new.asBed())
    # Sort by (chrom, start) so tabix can index — required when occpeaks/nucpos
    # come from `merge_chroms` of per-chromosome shards, which can produce
    # non-contiguous chromosome blocks.
    with open(bed_path) as f:
        lines = [ln for ln in f if ln.strip()]
    lines.sort(key=lambda ln: (ln.split('\t', 2)[0], int(ln.split('\t', 2)[1])))
    with open(bed_path, 'w') as f:
        f.writelines(lines)
    pysam.tabix_compress(bed_path, bed_path + '.gz', force=True)
    os.remove(bed_path)
    pysam.tabix_index(bed_path + '.gz', preset="bed", force=True)
 



