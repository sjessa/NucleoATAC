"""
Micro-benchmark for NucleoATAC hot functions.
Run before and after optimization to measure speedup.

Usage:
    python profiling/benchmark_hotfunctions.py
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

NREPS = 1000  # repetitions for calculateOccupancy (called thousands of times per chunk)
NREPS_COV = 100   # fewer reps since calculateCov is slower


def time_it(fn, nreps=100, label=""):
    """Run fn nreps times, return mean time in milliseconds."""
    # warm-up
    for _ in range(min(5, nreps)):
        fn()
    t0 = time.perf_counter()
    for _ in range(nreps):
        fn()
    elapsed = (time.perf_counter() - t0) / nreps * 1000  # ms per call
    print(f"  {label:50s}: {elapsed:10.4f} ms/call  ({nreps} reps)")
    return elapsed


def bench_calculateOccupancy():
    """Benchmark calculateOccupancy() — called per genomic position."""
    print("\n=== calculateOccupancy() ===")
    print(f"  (called ~2000x per 10kb chunk, once per 5bp step)")

    from nucleoatac.Occupancy import calculateOccupancy, OccupancyCalcParams

    # Simulate realistic inputs: upper=251 (yeast example uses 251)
    upper = 251
    rng = np.random.default_rng(42)

    # Build a fake insert distribution
    class FakeInsertDist:
        class nuc_fit:
            @staticmethod
            def get(lower, upper):
                v = rng.random(upper - lower) + 0.01
                return v / v.sum()
        class nfr_fit:
            @staticmethod
            def get(lower, upper):
                v = rng.random(upper - lower) + 0.01
                return v / v.sum()

    params = OccupancyCalcParams(0, upper, FakeInsertDist())
    inserts = rng.integers(0, 10, size=upper).astype(float)
    bias = rng.random(upper) + 0.5

    t = time_it(lambda: calculateOccupancy(inserts, bias, params),
                nreps=NREPS, label="calculateOccupancy(upper=251)")
    return t


def bench_calculateCov():
    """Benchmark calculateCov() — called per nucleosome candidate."""
    print("\n=== calculateCov() ===")
    print(f"  (called ~50x per chunk for nucleosome z-score computation)")

    import pyximport; pyximport.install(setup_args={"include_dirs": np.get_include()})
    from nucleoatac.multinomial_cov import calculateCov

    # Realistic size: vmat trimmed to lower=105, upper=251, w=60 → shape (146, 121)
    # flattened: 146*121 = 17666 elements
    n = 146 * 121  # ~17666
    rng = np.random.default_rng(42)
    p = rng.random(n)
    p /= p.sum()
    v = rng.random(n)
    r = 50

    print(f"  n = {n} (flattened V-plot elements)")
    t = time_it(lambda: calculateCov(p, v, r),
                nreps=NREPS_COV, label=f"calculateCov(n={n})")
    return t


def bench_makeFragmentMat():
    """Benchmark makeFragmentMatFromFragments() — per chunk."""
    print("\n=== makeFragmentMatFromFragments() ===")
    print(f"  (called once per chunk, processes all fragments)")

    # We can't easily benchmark the tabix I/O, but we can benchmark
    # the matrix filling loop vs the vectorized version.
    # Simulate 5000 fragments per chunk (typical ATAC-seq coverage)
    n_fragments = 5000
    lower, upper = 0, 251
    start, end = 100000, 110000
    nrow = upper - lower
    ncol = end - start

    rng = np.random.default_rng(42)
    frag_starts = rng.integers(start - upper, end, size=n_fragments)
    frag_ilens = rng.integers(lower, upper, size=n_fragments)
    frag_ends = frag_starts + frag_ilens

    def loop_version():
        mat = np.zeros((nrow, ncol), dtype=np.float64)
        for i in range(n_fragments):
            fs = int(frag_starts[i])
            fe = int(frag_ends[i])
            ilen = fe - fs
            if lower <= ilen < upper:
                r = ilen - lower
                c = (ilen - 1) // 2 + fs - start
                if 0 <= c < ncol and 0 <= r < nrow:
                    mat[r, c] += 1
        return mat

    def vectorized_version():
        mat = np.zeros((nrow, ncol), dtype=np.float64)
        fs = frag_starts
        ilens = frag_ends - fs
        mask = (ilens >= lower) & (ilens < upper)
        fs_m, ilens_m = fs[mask], ilens[mask]
        rows = ilens_m - lower
        cols = (ilens_m - 1) // 2 + fs_m - start
        valid = (cols >= 0) & (cols < ncol) & (rows >= 0) & (rows < nrow)
        np.add.at(mat, (rows[valid], cols[valid]), 1)
        return mat

    # Verify they produce the same result
    r1 = loop_version()
    r2 = vectorized_version()
    assert np.allclose(r1, r2), "Results differ!"
    print(f"  n_fragments = {n_fragments}, matrix shape = ({nrow}, {ncol})")

    t_loop = time_it(loop_version, nreps=20, label="loop version (current)")
    t_vec = time_it(vectorized_version, nreps=20, label="vectorized version (proposed)")
    print(f"  Speedup: {t_loop/t_vec:.1f}x")
    return t_loop, t_vec


def bench_makeBiasMat():
    """Benchmark BiasMat2D.makeBiasMat() row-by-row convolution."""
    print("\n=== BiasMat2D.makeBiasMat() ===")
    print(f"  (called once per chunk, ~2000 convolve calls)")

    from scipy import signal as scsignal
    import numpy as np

    # Realistic parameters: upper=251, lower=0, bias_track length varies
    lower, upper = 0, 251
    nrow = upper - lower  # 251 rows
    ncol_out = 500  # typical chunk width after convolution

    rng = np.random.default_rng(42)
    mid = upper // 2  # 125
    pattern = np.zeros((nrow, upper + (upper - 1) % 2))
    for i in range(lower, upper):
        pattern[i - lower, mid + (i - 1) // 2] = 1
        pattern[i - lower, mid - (i // 2)] = 1

    bias_len = ncol_out + pattern.shape[1] - 1
    bias = rng.random(bias_len)

    def loop_version():
        mat = np.zeros((nrow, ncol_out))
        for i in range(nrow):
            mat[i] = np.exp(np.convolve(bias, pattern[i, :], mode='valid'))
        return mat

    def correlate2d_version():
        return np.exp(scsignal.correlate2d(bias[np.newaxis, :], np.flipud(pattern), mode='valid'))

    # Verify equivalence
    r1 = loop_version()
    r2 = correlate2d_version()
    if not np.allclose(r1, r2, atol=1e-8):
        max_diff = np.max(np.abs(r1 - r2))
        print(f"  WARNING: results differ by up to {max_diff:.2e}")
    else:
        print(f"  Results match within 1e-8")

    print(f"  nrow = {nrow}, pattern shape = {pattern.shape}, bias len = {bias_len}")
    t_loop = time_it(loop_version, nreps=5, label="loop version (current)")
    t_corr = time_it(correlate2d_version, nreps=5, label="correlate2d version (proposed)")
    print(f"  Speedup: {t_loop/t_corr:.1f}x")
    return t_loop, t_corr


def bench_calculateCov_numpy():
    """Compare calculateCov Cython O(n^2) vs numpy O(n) formula."""
    print("\n=== calculateCov(): Cython O(n\u00b2) vs numpy O(n) ===")

    import pyximport; pyximport.install(setup_args={"include_dirs": np.get_include()})
    from nucleoatac.multinomial_cov import calculateCov as calculateCov_cython

    def calculateCov_numpy(p, v, r):
        return r * (np.dot(p, v**2) - np.dot(p, v)**2)

    n = 146 * 121
    rng = np.random.default_rng(42)
    p = rng.random(n)
    p /= p.sum()
    v = rng.random(n)
    r = 50

    # Verify equivalence
    ref = calculateCov_cython(p, v, r)
    new = calculateCov_numpy(p, v, r)
    diff = abs(ref - new)
    print(f"  n = {n}")
    print(f"  Cython result: {ref:.6f}, numpy result: {new:.6f}, diff: {diff:.2e}")

    t_cython = time_it(lambda: calculateCov_cython(p, v, r),
                       nreps=NREPS_COV, label="Cython O(n\u00b2)")
    t_numpy = time_it(lambda: calculateCov_numpy(p, v, r),
                      nreps=NREPS_COV, label="numpy O(n)")
    print(f"  Speedup: {t_cython/t_numpy:.1f}x")
    return t_cython, t_numpy


if __name__ == "__main__":
    print("=" * 70)
    print("NucleoATAC Hot-Function Benchmark")
    print("=" * 70)

    results = {}
    results['calculateOccupancy'] = bench_calculateOccupancy()
    results['calculateCov_cython'], results['calculateCov_numpy'] = bench_calculateCov_numpy()
    results['fragmat_loop'], results['fragmat_vec'] = bench_makeFragmentMat()
    results['biasmat_loop'], results['biasmat_corr'] = bench_makeBiasMat()

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  calculateOccupancy:   {results['calculateOccupancy']:.4f} ms/call (to be improved)")
    print(f"  calculateCov Cython:  {results['calculateCov_cython']:.4f} ms/call  →  numpy: {results['calculateCov_numpy']:.4f} ms/call  ({results['calculateCov_cython']/results['calculateCov_numpy']:.0f}x)")
    print(f"  fragmat loop:         {results['fragmat_loop']:.4f} ms/call  →  vector: {results['fragmat_vec']:.4f} ms/call  ({results['fragmat_loop']/results['fragmat_vec']:.0f}x)")
    print(f"  biasmat loop:         {results['biasmat_loop']:.4f} ms/call  →  corr2d: {results['biasmat_corr']:.4f} ms/call  ({results['biasmat_loop']/results['biasmat_corr']:.0f}x)")
