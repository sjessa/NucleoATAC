#!/usr/bin/env bash
# Run nucleoatac nuc per-chromosome using shared global occ+vprocess outputs,
# then merge into a single set of files.
#
# This demonstrates the recommended HPC parallelization strategy:
#   - occ + vprocess run once genome-wide (Phase 1)
#   - nuc run per-chromosome in parallel with shared --sizes, --vmat, --occ_track (Phase 2)
#   - merge_chroms combines the per-chromosome nuc outputs (Phase 3)
#
# Usage:
#   bash scripts/run_parallel_nuc.sh <global_occ_dir> <out_dir> [cores]
#
# Example (using the regression test outputs as the global occ):
#   bash scripts/run_parallel_nuc.sh \
#       tests/py3_regression_output \
#       profiling/py3_parallel \
#       2

set -euo pipefail

GLOBAL_DIR="${1:-tests/py3_regression_output}"
OUT_DIR="${2:-profiling/py3_parallel}"
CORES="${3:-2}"

PYTHON=/oak/stanford/groups/akundaje/sjessa/software/miniconda3/envs/nucleoatac2/bin/python
CHROMS="chrII chrIV chrIX chrV chrVI chrVII chrXI chrXII chrXIII chrXIV"

GLOBAL_SIZES="$GLOBAL_DIR/example.fragmentsizes.txt"
GLOBAL_VMAT="$GLOBAL_DIR/example.VMat"
GLOBAL_OCC="$GLOBAL_DIR/example.occ.bedgraph.gz"

echo "Global occ dir: $GLOBAL_DIR"
echo "Output dir:     $OUT_DIR"
echo "Chromosomes:    $CHROMS"
echo ""

mkdir -p "$OUT_DIR"

# Phase 2: per-chromosome nuc
for chrom in $CHROMS; do
    echo "--- nuc: $chrom ---"
    $PYTHON -c "
from nucleoatac.cli import nucleoatac_parser, nucleoatac_main
parser = nucleoatac_parser()
args = parser.parse_args([
    'nuc',
    '--bed',       'example/example.bed',
    '--bam',       'example/example.bam',
    '--fasta',     'example/sacCer3.fa',
    '--vmat',      '$GLOBAL_VMAT',
    '--sizes',     '$GLOBAL_SIZES',
    '--occ_track', '$GLOBAL_OCC',
    '--chroms_keep', '$chrom',
    '--out',       '$OUT_DIR/example.$chrom',
    '--cores',     '$CORES',
])
nucleoatac_main(args)
"
done

# Phase 3: merge per-chromosome nuc outputs
echo ""
echo "--- merge_chroms ---"
CHROM_LIST=$(echo $CHROMS | tr ' ' ',')
$PYTHON -c "
from nucleoatac.cli import nucleoatac_parser, nucleoatac_main
parser = nucleoatac_parser()
args = parser.parse_args([
    'merge_chroms',
    '--prefix', '$OUT_DIR/example',
    '--chroms', '$CHROM_LIST',
    '--out',    '$OUT_DIR/example',
])
nucleoatac_main(args)
"

echo ""
echo "Done. Merged outputs in $OUT_DIR/example.*"
