from src.binom_fns import *
from src.two_step_disjoint import check_eq_2_4, check_eq_2_5

def check_eq_2_1_dash(n, p, k, r, c_1=1.0999):
    r_dash = round(r**2 / n)
    total = 0.0
    for i in range(0, k):
        total += T(r_dash, p, i)*(B(r-r_dash, p, k-i))**2
    p_dash = B(r_dash, p, k) + total
    expect = round(n * p_dash, 6)
    dev = expect * (c_1 - 1)
    return r >= (expect - dev) and r <= (expect + dev)

def check_eq_2_2_dash(n, p, k, r):
    r_dash = round(r**2 / n)
    total = 0.0
    for i in range(0, k):
        total += T(r_dash, p, i)*(B(r-r_dash, p, k-i))*B((r/2)-r_dash, p, k-i)
    p_dash = B(r_dash, p, k) + total
    expect = B(n, p_dash, r/10)
    return round(expect, 6) == 0.0

def check_eq_2_3_dash(n, p, k, r):
    r_dash = round(r**2 / n)
    r_ddash = round(r**3 / n**2)
    r_caret = r_dash - r_ddash
    p_dash = 0.0
    for i in range(0, k+1):
        total_j = 0.0
        for j in range(0, k-i+1):
            total_l = 0.0
            for l in range(0, k-i-j+1):
                total_m = 0.0
                for m in range(0, k-i-j-l+1):
                    total_m += T(r_ddash, p, m) * (B(r-2*r_dash+r_ddash, p, k-i-j-m) * B(r-2*r_dash+r_ddash, p, k-j-l-m) * (1-B(r-2*r_dash + r_ddash, p, k-i-l-m)))
                total_l += T(r_caret, p, l) * total_m
            total_j += T(r_caret, p, j) * total_l
        p_dash += T(r_caret, p, i) * total_j
    expect = B(n, p_dash, 2*r/3)
    return round(expect, 5) == 1.0

def check_eq_2_4_dash(n, p, k, r):
    return check_eq_2_4(n, p, k, r)

def check_eq_2_5_dash(n, p, k, r):
    return check_eq_2_5(n, p, k, r)

def check_eq_2_6_dash(n, p, k, r, t=1):
    r_dash = round(r**2 / n)
    total = 0.0
    for i in range(0, k):
        total += T(r_dash, p, i)*(B(r-r_dash, p, k-i))**2
    p_dash = B(r_dash, p, k) + total
    p_ddash = p_dash * B(n, p_dash, k)
    expect = B(r, p_ddash, r/2)
    return round(expect, 6) == 0.0

def run_twostep_shared_checks(n, p, k, r, c_1=1.0999):
    print("For a Shared Representation with Two-Step Mechanisms,")
    header = ' set n={0}, d={1}, k={2}, c_1={3},\n and test r={4}:\n'
    print(header.format(n, round(n*p), k, c_1, r))
    check_1 = check_eq_2_1_dash(n, p, k, r, c_1)
    print("Equation 2.1':", check_1)
    check_2 = check_eq_2_2_dash(n, p, k, r)
    print("Equation 2.2':", check_2)
    check_3 = check_eq_2_3_dash(n, p, k, r)
    print("Equation 2.3':", check_3)
    check_4 = check_eq_2_4_dash(n, p, k, r)
    print("Equation 2.4':", check_4)
    check_5 = check_eq_2_5_dash(n, p, k, r)
    print("Equation 2.5':", check_5)
    check_6 = check_eq_2_6_dash(n, p, k, r)
    print("Equation 2.6':", check_6, "\n")
    if check_1 and check_2 and check_3 and check_4 and check_5 and check_6:
        return "Passed"
    return "Failed"