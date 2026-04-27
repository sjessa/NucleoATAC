#!/usr/bin/env python
"""Compare NucleoATAC pipeline outputs and generate summary figures.

Supports 2-way or 3-way comparison (--test2 / --label-test2 for a third dataset).

Typical usage
-------------
# 2-way: Python 2 reference vs Python 3 output
python scripts/evaluate_pipeline.py \
    --ref tests/py2_reference \
    --test tests/py3_regression_output \
    --label-ref "Python 2" --label-test "Python 3 (whole-genome)" \
    --out profiling/py2_vs_py3/

# 3-way: add per-chromosome parallel as a third comparison
python scripts/evaluate_pipeline.py \
    --ref tests/py2_reference \
    --test tests/py3_regression_output \
    --test2 profiling/py3_parallel \
    --label-ref "Python 2" \
    --label-test "Python 3 (whole-genome)" \
    --label-test2 "Python 3 (per-chrom)" \
    --out profiling/three_way/
"""

import argparse
import gzip
import os

import numpy as np
import seaborn
import matplotlib.pyplot as plt
from scipy import stats

seaborn.set_style("whitegrid")

COLORS = ["#4C72B0", "#DD8452", "#55A868"]  # blue, orange, green


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_bedgraph_gz(path, chroms=None):
    data = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if chroms and parts[0] not in chroms:
                continue
            key = (parts[0], int(parts[1]))
            data[key] = float(parts[3])
    return data


def load_bed_gz(path, n_float_cols, has_str_col=False, chroms=None):
    data = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if chroms and parts[0] not in chroms:
                continue
            key = (parts[0], int(parts[1]))
            floats = tuple(float(parts[i]) for i in range(3, 3 + n_float_cols))
            label = parts[3 + n_float_cols] if has_str_col else None
            data[key] = (floats, label)
    return data


def load_distribution_txt(path):
    with open(path) as f:
        lines = f.read().strip().split("\n")
    lower = int(lines[1])
    upper = int(lines[3])
    values = np.array([float(x) for x in lines[5].split("\t")])
    return lower, upper, values


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _save(fig, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


def _scatter_with_identity(ax, x, y, xlabel, ylabel, subtitle, color):
    ax.scatter(x, y, s=4, alpha=0.4, color=color, rasterized=True)
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, "k--", lw=1, zorder=5)
    r, _ = stats.pearsonr(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(subtitle, fontsize=10)
    ax.text(0.05, 0.92, f"r = {r:.6f}", transform=ax.transAxes, fontsize=9)
    return r


# ---------------------------------------------------------------------------
# Comparison functions
# ---------------------------------------------------------------------------

def compare_occupancy(ref_dir, test_dirs, out_dir, prefix, labels, chroms=None):
    n = len(test_dirs)
    ref_path = os.path.join(ref_dir, f"{prefix}.occ.bedgraph.gz")
    if not os.path.exists(ref_path):
        print("  [skip] occ.bedgraph.gz not found in ref directory")
        return {}

    ref = load_bedgraph_gz(ref_path, chroms=chroms)
    test_data, avail_labels = [], []
    for td, tlab in zip(test_dirs, labels[1:]):
        p = os.path.join(td, f"{prefix}.occ.bedgraph.gz")
        if not os.path.exists(p):
            print(f"  [skip occ] {td} has no occ.bedgraph.gz")
            continue
        test_data.append(load_bedgraph_gz(p, chroms=chroms))
        avail_labels.append(tlab)

    if not test_data:
        print("  [skip] occ.bedgraph.gz not found in any test directory")
        return {}
    n = len(test_data)

    # --- scatter panels (ref vs each available test) ---
    fig, axes = plt.subplots(1, n + 1, figsize=(5 * (n + 1), 4))
    if n + 1 == 2:
        axes = list(axes)
    rs = []
    for i, (td, tlab) in enumerate(zip(test_data, avail_labels)):
        shared = sorted(set(ref) & set(td))
        ref_vals = np.array([ref[k] for k in shared])
        test_vals = np.array([td[k] for k in shared])
        r = _scatter_with_identity(
            axes[i], ref_vals, test_vals,
            labels[0], tlab,
            f"Occupancy: {labels[0]} vs {tlab}\n(n={len(shared):,})",
            color=COLORS[i + 1])
        rs.append(r)

    # --- KDE panel ---
    ax_kde = axes[n]
    seaborn.kdeplot(list(ref.values()), ax=ax_kde, label=labels[0],
                    fill=True, alpha=0.35, color=COLORS[0])
    for i, (td, tlab) in enumerate(zip(test_data, avail_labels)):
        seaborn.kdeplot(list(td.values()), ax=ax_kde, label=tlab,
                        fill=True, alpha=0.35, color=COLORS[i + 1])
    ax_kde.set_xlabel("Occupancy")
    ax_kde.set_title("Occupancy distribution")
    ax_kde.legend(fontsize=8)

    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "occupancy.png"))

    out = {}
    for i, tlab in enumerate(avail_labels):
        key = tlab.replace(" ", "_").replace("(", "").replace(")", "")
        out[f"occ_pearson_r_{key}"] = rs[i]
    return out


def compare_signal(ref_dir, test_dirs, out_dir, prefix, labels, chroms=None):
    n = len(test_dirs)
    ref_path = os.path.join(ref_dir, f"{prefix}.nucleoatac_signal.bedgraph.gz")
    if not os.path.exists(ref_path):
        print("  [skip] nucleoatac_signal.bedgraph.gz not found in ref directory")
        return {}

    ref = load_bedgraph_gz(ref_path, chroms=chroms)
    test_data = []
    for td in test_dirs:
        p = os.path.join(td, f"{prefix}.nucleoatac_signal.bedgraph.gz")
        if not os.path.exists(p):
            print(f"  [skip] nucleoatac_signal.bedgraph.gz not found in {td}")
            return {}
        test_data.append(load_bedgraph_gz(p, chroms=chroms))

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]
    rs = []
    for i, (td, tlab) in enumerate(zip(test_data, labels[1:])):
        shared = sorted(set(ref) & set(td))
        ref_vals = np.array([ref[k] for k in shared])
        test_vals = np.array([td[k] for k in shared])
        r = _scatter_with_identity(
            axes[i], ref_vals, test_vals,
            labels[0], tlab,
            f"NucleoATAC signal: {labels[0]} vs {tlab}\n(n={len(shared):,})",
            color=COLORS[i + 1])
        rs.append(r)

    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "signal.png"))

    out = {}
    for i, tlab in enumerate(labels[1:]):
        key = tlab.replace(" ", "_").replace("(", "").replace(")", "")
        out[f"signal_pearson_r_{key}"] = rs[i]
    return out


def compare_nucpos(ref_dir, test_dirs, out_dir, prefix, labels, chroms=None):
    n = len(test_dirs)
    ref_path = os.path.join(ref_dir, f"{prefix}.nucpos.bed.gz")
    if not os.path.exists(ref_path):
        print("  [skip] nucpos.bed.gz not found in ref directory")
        return {}

    ref = load_bed_gz(ref_path, n_float_cols=10, chroms=chroms)
    test_data = []
    for td in test_dirs:
        p = os.path.join(td, f"{prefix}.nucpos.bed.gz")
        if not os.path.exists(p):
            print(f"  [skip] nucpos.bed.gz not found in {td}")
            return {}
        test_data.append(load_bed_gz(p, n_float_cols=10, chroms=chroms))

    fig, axes = plt.subplots(1, n + 1, figsize=(5 * (n + 1), 4))

    # --- overlap bar chart ---
    ax_bar = axes[0]
    categories = []
    for tlab in labels[1:]:
        categories.append(f"{tlab}\nonly")
    categories = [f"{labels[0]}\nonly"] + categories + ["Shared\n(all)"]

    ref_keys = set(ref.keys())
    all_test_keys = [set(td.keys()) for td in test_data]
    shared_all = ref_keys.copy()
    for tk in all_test_keys:
        shared_all &= tk

    counts = [len(ref_keys - set.union(*all_test_keys))]
    for i, tk in enumerate(all_test_keys):
        others = set.union(*(all_test_keys[:i] + all_test_keys[i+1:] + [ref_keys])) if len(all_test_keys) > 1 else ref_keys
        counts.append(len(tk - ref_keys - set.union(*[atk for j, atk in enumerate(all_test_keys) if j != i])))
    counts.append(len(shared_all))

    bar_colors = [COLORS[0]] + COLORS[1:n+1] + ["grey"]
    ax_bar.bar(categories, counts, color=bar_colors)
    ax_bar.set_ylabel("# nucleosome calls")
    ax_bar.set_title("Nucleosome call overlap")
    ax_bar.tick_params(axis='x', rotation=45)
    for label in ax_bar.get_xticklabels():
        label.set_ha('right')

    # --- occupancy scatter panels (ref vs each test) ---
    rs = []
    for i, (td, tlab) in enumerate(zip(test_data, labels[1:])):
        shared = sorted(set(ref) & set(td))
        ref_occ = np.array([ref[k][0][4] for k in shared])
        test_occ = np.array([td[k][0][4] for k in shared])
        r = _scatter_with_identity(
            axes[i + 1], ref_occ, test_occ,
            f"{labels[0]} occupancy", f"{tlab} occupancy",
            f"Nucpos occupancy: {labels[0]} vs {tlab}\n(n={len(shared):,})",
            color=COLORS[i + 1])
        rs.append(r)

    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "nucpos.png"))

    out = {"nucpos_n_ref": len(ref)}
    for i, (td, tlab) in enumerate(zip(test_data, labels[1:])):
        key = tlab.replace(" ", "_").replace("(", "").replace(")", "")
        out[f"nucpos_n_{key}"] = len(td)
        out[f"nucpos_n_shared_ref_{key}"] = len(set(ref) & set(td))
        out[f"nucpos_occ_pearson_r_{key}"] = rs[i]
    return out


def compare_fragmentsizes(ref_dir, test_dirs, out_dir, prefix, labels):
    stats_out = {}
    for suffix, title in [("fragmentsizes.txt", "Fragment size distribution"),
                           ("nuc_dist.txt", "Nucleosomal insert distribution")]:
        ref_path = os.path.join(ref_dir, f"{prefix}.{suffix}")
        if not os.path.exists(ref_path):
            print(f"  [skip] {suffix} not found in ref directory")
            continue

        ref_lower, ref_upper, ref_vals = load_distribution_txt(ref_path)
        x = np.arange(ref_lower, ref_upper)

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(x, ref_vals, label=labels[0], alpha=0.85, color=COLORS[0])

        for i, (td, tlab) in enumerate(zip(test_dirs, labels[1:])):
            tp = os.path.join(td, f"{prefix}.{suffix}")
            if not os.path.exists(tp):
                print(f"  [skip] {suffix} not found in {td}")
                continue
            _, _, test_vals = load_distribution_txt(tp)
            ax.plot(x, test_vals, label=tlab, alpha=0.85,
                    linestyle="--", color=COLORS[i + 1])
            key = f"{suffix.replace('.txt','').replace('.','_')}_{tlab.replace(' ','_').replace('(','').replace(')','')}_max_diff"
            stats_out[key] = float(np.max(np.abs(ref_vals - test_vals)))

        ax.set_xlabel("Insert size (bp)")
        ax.set_ylabel("Frequency")
        ax.set_title(title)
        ax.legend(fontsize=8)
        fig.tight_layout()
        _save(fig, os.path.join(out_dir, suffix.replace(".txt", ".png")))

    return stats_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare NucleoATAC pipeline outputs and generate summary figures.")
    parser.add_argument("--ref", required=True, help="Reference output directory")
    parser.add_argument("--test", required=True, help="Test output directory")
    parser.add_argument("--test2", default=None, help="Optional second test directory")
    parser.add_argument("--prefix", default="example", help="File prefix (default: example)")
    parser.add_argument("--out", default="scripts/figures/", help="Output directory for figures")
    parser.add_argument("--label-ref", default="Reference")
    parser.add_argument("--label-test", default="Test")
    parser.add_argument("--label-test2", default="Test 2")
    parser.add_argument("--chroms", default=None,
                        help="Comma-separated list of chromosomes to include (default: all)")
    args = parser.parse_args()

    chroms = set(args.chroms.split(",")) if args.chroms else None
    test_dirs = [args.test]
    labels = [args.label_ref, args.label_test]
    if args.test2:
        test_dirs.append(args.test2)
        labels.append(args.label_test2)

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    print(f"Comparing:\n  ref:   {args.ref}")
    for i, td in enumerate(test_dirs):
        print(f"  test{i+1}: {td}")
    if chroms:
        print(f"  chromosomes: {', '.join(sorted(chroms))}")
    print(f"Output: {out_dir}\n")

    all_stats = {}
    print("Occupancy...")
    all_stats.update(compare_occupancy(args.ref, test_dirs, out_dir, args.prefix, labels, chroms))
    print("NucleoATAC signal...")
    all_stats.update(compare_signal(args.ref, test_dirs, out_dir, args.prefix, labels, chroms))
    print("Nucleosome calls...")
    all_stats.update(compare_nucpos(args.ref, test_dirs, out_dir, args.prefix, labels, chroms))
    print("Fragment size distributions...")
    all_stats.update(compare_fragmentsizes(args.ref, test_dirs, out_dir, args.prefix, labels))

    print("\n--- Summary ---")
    for k, v in all_stats.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6g}")
        else:
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
