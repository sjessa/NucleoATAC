"""
Closed-form O(n) covariance calculation for multinomial distributions.

Replaces nucleoatac/multinomial_cov.pyx (O(n^2) Cython loop).

Mathematical identity:
    Var[sum_i v_i * X_i] for X ~ Multinomial(r, p)
    = r * (sum_i p_i * v_i^2 - (sum_i p_i * v_i)^2)

This is equivalent to the pairwise loop:
    sum_i p_i(1-p_i)*v_i^2  +  2 * sum_{i<j} (-p_i*p_j*v_i*v_j)
    = sum_i p_i*v_i^2 - (sum_i p_i*v_i)^2
"""

import numpy as np


def calculateCov(p, v, r):
    """
    Calculate variance of sum(v_i * X_i) for X ~ Multinomial(r, p).

    Args:
        p: probability vector (1D numpy array, must sum to 1)
        v: value vector (1D numpy array, same shape as p)
        r: number of trials (int)

    Returns:
        Variance as a float.
    """
    return r * (np.dot(p, v ** 2) - np.dot(p, v) ** 2)
