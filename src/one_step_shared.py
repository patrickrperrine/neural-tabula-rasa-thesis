from src.binom_fns import *

def check_eq_2_1_ddash(n, p, k_m, r):
    r_dash = round(r**2 / n)
    p_dash = B((2*r) - r_dash, p, k_m)
    expect = round(n * p_dash, 2)
    if n < 10**6:
        stdev = 2 * np.sqrt(expect)
        return r >= (expect - stdev) and r <= (expect + stdev)
    return compare_first_x_digits(expect, r, 2)

def check_eq_2_2_ddash(n, p, k_m, r):
    p_dash = B(3*r/2, p, k_m) 
    expect = B(n, p_dash, r/10)
    return round(expect, 6) == 0.0

def eq_2_3_ddash_p_dash_disjoint(n, p, k_m, r):
    p_dash = 0.0
    for s in range(0, k_m):
        p_dash += T(r, p, s) * (B(r, p, k_m - s)) * (1 - B(r, p, k_m - s))
    return p_dash

def eq_2_3_ddash_p_dash_shared(n, p, k_m, r):
    r_dash = round(r**2 / n)
    r_ddash = round(r**3 / n**2)
    r_caret = r_dash - r_ddash
    r_sharp = r - 2*r_dash + r_ddash # double minus in paper?
    p_dash = 0.0
    for s in range(0, k_m):
        total_i = 0.0
        for i in range(0, k_m-s):
            total_j = 0.0
            for j in range(0, k_m-i-s):
                total_m = 0.0
                for m in range(0, k_m-i-j-s):
                    total_l = 0.0
                    for l in range(0, k_m-i-j-m-s):
                        total_l += T(r_caret, p, l) * ((B(r_sharp, p, k_m-s-i-j-l-m)) * (1 - B(r_sharp, p, k_m-s-i-j-l-m)))
                    total_m += T(r_ddash, p, m) * total_l
                total_j += T(r_caret, p, j) * total_m
            total_i += T(r_caret, p, i) * total_j
        p_dash += T(r_sharp, p, s) * total_i
    return p_dash

def check_eq_2_3_ddash(n, p, k_m, r, complete_share=True):
    if not complete_share:
        p_dash = eq_2_3_ddash_p_dash_disjoint(n, p, k_m, r)
    else:
        p_dash = eq_2_3_ddash_p_dash_shared(n, p, k_m, r)
    expect = B(n, p_dash, 2*r/3)
    return round(expect, 6) == 1.0

def check_eq_2_4_ddash(n, p, k_a, r):
    p_dash = p * B(r, p, k_a)
    Y = B(n, p_dash, k_a)
    return round(Y, 6) == 1.0

def check_eq_2_5_ddash(n, p, k_a, r):
    p_dash = p * B(r/2, p, k_a)
    p_ddash = B(n, p_dash, k_a)
    expect = B(r, p_ddash, r/2)
    return round(expect, 6) == 0.0

def check_eq_2_6_ddash(n, p, k_a, r):
    r_dash = round(r**2 / n)
    total = 0.0
    for i in range(0, k_a):
        total += T(r_dash, p, i)*(B(r-r_dash, p, k_a-i))**2
    p_dash = p * (B(r_dash, p, k_a) + total)
    p_ddash = B(n, p_dash, k_a)
    expect = B(r, p_ddash, r/2)
    return round(expect, 6) == 0.0

def run_onestep_shared_checks(n, p, k, k_adj, r, complete_share=True):
    k_a = k
    k_m = round(k_adj * k_a)
    if not complete_share:
        print("For a Partially-Shared Representation with One-Step Mechanisms,")
    else:
        print("For a Shared Representation with One-Step Mechanisms,")
    header = ' set n={0}, d={1}, k_a={2}, k_m={3},\n and test r={4}:\n'
    print(header.format(n, round(n*p), k_a, k_m, r))
    check_1 = check_eq_2_1_ddash(n, p, k_m, r)
    print("Equation 2.1'':", check_1)
    check_2 = check_eq_2_2_ddash(n, p, k_m, r)
    print("Equation 2.2'':", check_2)
    print("* Calculations of Equation 2.3'' (Shared) take a while... *", end='\r')
    check_3 = check_eq_2_3_ddash(n, p, k_m, r, complete_share)
    if not complete_share:
        print("Equation 2.3'' (Disjoint):", check_3, ' '*45)
    else:
        print("Equation 2.3'':", check_3, ' '*45)
    print("Equation 2.3'':", check_3)
    check_4 = check_eq_2_4_ddash(n, p, k_a, r)
    print("Equation 2.4'':", check_4)
    check_5 = check_eq_2_5_ddash(n, p, k_a, r)
    print("Equation 2.5'':", check_5)
    check_6 = check_eq_2_6_ddash(n, p, k_a, r)
    print("Equation 2.6'':", check_6, "\n")
    if check_1 and check_2 and check_3 and check_4 and check_5 and check_6:
        return "Passed"
    return "Failed"