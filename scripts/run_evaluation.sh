#!/usr/bin/env bash
# Generate pipeline comparison figures.
#
# Produces a 3-way comparison:
#   - Python 2 reference
#   - Python 3 whole-genome run
#   - Python 3 per-chromosome parallel run (shared global FLD)
#
# Run the parallel nuc pipeline first:
#   bash scripts/run_parallel_nuc.sh tests/py3_regression_output profiling/py3_parallel 2
#
# Then run this script:
#   bash scripts/run_evaluation.sh

set -euo pipefail

PYTHON=/oak/stanford/groups/akundaje/sjessa/software/miniconda3/envs/nucleoatac2/bin/python

$PYTHON scripts/evaluate_pipeline.py \
    --ref   tests/py2_reference \
    --test  tests/py3_regression_output \
    --test2 profiling/py3_parallel \
    --prefix example \
    --label-ref   "Python 2" \
    --label-test  "Python 3 (whole-genome)" \
    --label-test2 "Python 3 (per-chrom)" \
    --out profiling/py2_vs_py3/
