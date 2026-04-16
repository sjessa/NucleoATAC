"""
Generate before/after performance comparison plots for NucleoATAC optimizations.

Usage:
    python profiling/plot_comparison.py
"""

import re
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_mprof(path):
    """Load an mprof .dat file. Returns (timestamps, memory_mb)."""
    times, mems = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('MEM'):
                parts = line.split()
                mems.append(float(parts[1]))
                times.append(float(parts[2]))
    times = np.array(times)
    mems = np.array(mems)
    times -= times[0]  # relative time starting at 0
    return times, mems


def plot_mprof_comparison(before_dat, after_dat, outfile):
    """Plot memory-over-time for before and after, side by side."""
    t_before, m_before = load_mprof(before_dat)
    t_after, m_after = load_mprof(after_dat)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('NucleoATAC Runtime & Memory: Before vs After Optimization',
                 fontsize=13, fontweight='bold')

    for ax, times, mems, label, color in [
        (axes[0], t_before, m_before, 'before optimization', '#333333'),
        (axes[1], t_after,  m_after,  'after optimization',  '#2271b3'),
    ]:
        ax.plot(times, mems, '.', ms=2, color=color)
        peak = mems.max()
        ax.axhline(peak, color='red', linestyle='--', linewidth=0.8, label=f'peak: {peak:.0f} MB')
        ax.set_xlabel('time (in seconds)', fontsize=11)
        ax.set_ylabel('memory used (in MB)', fontsize=11)
        total_time = times[-1]
        ax.set_title(f'NucleoATAC v0.5.0\n{label}\n(total: {total_time:.0f}s, peak: {peak:.0f} MB)',
                     fontsize=10, fontfamily='monospace')
        ax.legend(fontsize=9)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outfile}")


def plot_function_speedups(outfile):
    """Bar chart of per-function speedups from benchmark results."""
    functions = [
        'calculateOccupancy()\n(vectorized)',
        'calculateCov()\n(O(n²) → O(n))',
        'makeFragmentMat()\n(vectorized)',
        'makeBiasMat()\n(vectorized)',
    ]
    before_ms = [1.5597, 233.5, 6.47, 9.857]
    after_ms  = [0.1911, 0.0262, 2.297, 1.038]
    speedups = [b/a for b, a in zip(before_ms, after_ms)]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: absolute times
    x = np.arange(len(functions))
    w = 0.35
    bars_b = axes[0].bar(x - w/2, before_ms, w, label='before', color='#999999')
    bars_a = axes[0].bar(x + w/2, after_ms, w, label='after', color='#2271b3')
    axes[0].set_yscale('log')
    axes[0].set_ylabel('time per call (ms, log scale)', fontsize=11)
    axes[0].set_title('Function-level timing: before vs after', fontsize=11)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(functions, fontsize=8)
    axes[0].legend()
    axes[0].grid(True, axis='y', alpha=0.3)
    for bar, val in zip(bars_b, before_ms):
        axes[0].text(bar.get_x() + bar.get_width()/2, val * 1.1,
                     f'{val:.3g}ms', ha='center', va='bottom', fontsize=7, color='#555555')
    for bar, val in zip(bars_a, after_ms):
        axes[0].text(bar.get_x() + bar.get_width()/2, val * 1.1,
                     f'{val:.3g}ms', ha='center', va='bottom', fontsize=7, color='#2271b3')

    # Right: speedup factors
    colors = ['#e69f00' if s < 100 else '#cc79a7' for s in speedups]
    bars = axes[1].bar(x, speedups, color=colors)
    axes[1].set_ylabel('speedup factor (×)', fontsize=11)
    axes[1].set_title('Speedup per function', fontsize=11)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(functions, fontsize=8)
    axes[1].grid(True, axis='y', alpha=0.3)
    for bar, s in zip(bars, speedups):
        label = f'{s:.0f}×' if s >= 10 else f'{s:.1f}×'
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                     label, ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outfile}")


if __name__ == '__main__':
    before_dat = os.path.join(REPO, 'profiling/before/mprof.dat')
    after_dat  = os.path.join(REPO, 'profiling/after/mprof.dat')

    if os.path.exists(before_dat) and os.path.exists(after_dat):
        plot_mprof_comparison(before_dat, after_dat,
                              os.path.join(REPO, 'profiling/memory_comparison.png'))
    else:
        print("Warning: mprof .dat files not found, skipping memory plot")

    plot_function_speedups(os.path.join(REPO, 'profiling/speedup_comparison.png'))
    print("Done.")
