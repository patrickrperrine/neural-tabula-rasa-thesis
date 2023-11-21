import numpy as np
from scipy.stats import binom

def B(r, p, k):
    return binom.sf(k-1, r, p)

def T(r, p, j):
    return binom.pmf(j, r, p)

def compare_first_x_digits(a, b, x):
    return str(a)[:x] == str(b)[:x]
