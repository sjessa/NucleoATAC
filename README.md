# NucleoATAC

## Update 2024/08/26

This is a fork of the NucleoATAC package authored by Alicia Schep in the Greenleaf Lab,
a package for calling nucleosome occupancy from ATAC-seq data.
We modify the package to take fragments files as inputs instead of BAM files.
This fork is not being officially maintained. Please refer to the original package
for the last stable version.


### Installation

1. Create a conda environment

```bash
# loading the module first seems to be necessary for properly loading cython
module load python/2.7.13
conda create -n "nucleoatac" python=2.7.13
conda activate nucleoatac
```

2. Install the package

```bash
cd NucleoATAC
pip install --editable .
```

**_NOTE_**: to edit the package and install the new version, update the version
number in `setup.py` for documentation, then re-run the above pip command.
Test the install with `pyatac --version`.


### Updates to functionality

**_NOTE_**: `nucleoatac run` doesn't work with the modifications for fragment
sizes yet. Run commands individually as per [the docs](https://nucleoatac.readthedocs.io/en/latest/nucleoatac/).

#### Allowing fragments files as input
- Goal: allow nucleosome calling to be performed with tabix-indexed compressed fragments files in BED format
- The `pyatac.fragmentsizes` module has been updated so that the `FragmentSizes` class
  can now accept a fragments file instead of a BAM file for calculation of the fragment
  size distribution. This makes use of new functions in the `nucleoatac.fragments_handling` module.
- The following commands now can now accept a fragments file instead of BAM file:
  - `pyatac sizes`
  - `nucleoatac occ` 
  - `nucleoatac nuc`
- The `pyatac.chunkmat2d.FragmentMat2D` module has been updated so that the `FragmentMat2D` class
  can now accept a fragments file instead of a BAM file for calculation of the fragment
  size distribution, which will fetch the chunk region from the fragments file
  using `pysam`'s tabix interface, and then populate the 2D matrix
- Add chunk handling for the fragemnts file as well
- When commands required BAM files for calculation of fragment sizes, the fragment sizes output
  from previous commands is now required instead


#### Restrict analysis to chromosomes of interest
- Goal: calculate the fragment size distribution on the whole dataset, but restrict
  the occupancy calculation to a subset of chromosomes to allow users to debug
  or further parallelize the analysis
- The `nucleotac occ` command now accepts a `--chroms_keep` argument that allows users
  to specify a list of chromosomes to restrict the occupancy calculation to, and a chunk
  list will only be generated for these chromosomes.
- This can be used to manually parallelize the runs, or do a quick run for debugging


#### Misc
- `nucleoatac occ` and `pyatac sizes` changed to produce .pdf instead of .eps files
- More verbose messaging added in `nucleoatac occ` and `nucleoatac nuc`



### Timing

For reference, running only on chr1 (with 7,267,633 fragments) on 4 cores and 64G of memory took about 14 minutes.

```
Command run:  /path/to/.local/bin/nucleoatac occ --fragments ../data/Eye_c11__sorted.tsv.gz --bed ../data/Eye_c11__peaks_overlap_filtered.narrowPeak --fasta
 /path/to/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --out ../out/04-Eye_c11__chr1 --cores 4 --chroms_keep chr
1
nucleoatac version 0.3.4
start run at: 2024-08-27 14:48
---------Computing Occupancy and Nucleosomal Insert Distribution----------------
@ NOTE: restricting analysis to chromosomes: chr1
@ calculating fragment sizes...
@ fitting fragment size distribution...
@ calculating occupancy...
@ compressing and indexing output files...
@ making figure
@ done.
end run at: 2024-08-27 15:02
```


## Previously

**This package is no longer being actively maintained; feel free to post issues that others in the community may respond to, but this package will likely not be updated further. Additionally, if anyone wants to maintain a fork of the package or has developed an alternative package for similar purposes, would be happy to link to that repo here.**

Python package for calling nucleosomes using ATAC-seq data.
Also includes general scripts for working with paired-end ATAC-seq data (or potentially other paired-end data).

Please cite our paper at [Genome Research](http://genome.cshlp.org/content/25/11/1757) if you use this tool in your research.

Please use GitHub Issues to bring up any errors that occur with software rather than emailing authors.

Note on Versions:  

* version 0 represents code used for biorxiv manuscript
* version 0.2.1 was used for Genome Research manuscript (See Supplemental Information as well)

Documentation  can be found at http://nucleoatac.readthedocs.org/en/latest/

If you want to easily read in NucleoATAC outputs into R for further processing or exploration, check out [NucleoATACR](https://github.com/GreenleafLab/NucleoATACR/)

Currently NucleoATAC only supports Python 2.7 (No Python 3).
