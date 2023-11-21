from src.binom_fns import *

def check_eq_2_1(n, p, k, r):
    expect = round(n * B(r, p, k)**2, 2)
    if n < 10**6:
        stdev = 2 * np.sqrt(expect)
        return r >= (expect - stdev) and r <= (expect + stdev)
    return compare_first_x_digits(expect, r, 2)

def check_eq_2_2(n, p, k, r):
    p_dash = B(r/2, p, k) * B(r, p, k)
    expect = B(n, p_dash, r/10)
    return round(expect, 6) == 0.0

def check_eq_2_3(n, p, k, r):
    p_dash = (1 - B(r, p, k)) * (B(r, p, k))**2
    expect = B(n, p_dash, 2*r/3)
    return round(expect, 6) == 1.0

def check_eq_2_4(n, p, k, r):
    p_dash = p * B(r, p, k)
    Y = B(n, p_dash, k)
    return round(Y, 6) == 1.0

def check_eq_2_5(n, p, k, r):
    p_dash = B(n, p * B(r/2, p, k), k)
    expect = B(r, p_dash, r/2)
    return round(expect, 6) == 0.0

def check_eq_2_6(n, p, k, r, t):
    p_dash = p * (1 - (1 - B(r, p, k))**t)
    p_ddash = B(n, p_dash * B(r, p, k), k)
    approx = B(r, p_ddash, r/2)
    return round(approx, 6) == 0.0

def run_twostep_disjoint_checks(n, p, k, r, t=1):
    print("For a Disjoint Representation with Two-Step Mechanisms,")
    header = ' set n={0}, d={1}, k={2}, t={3},\n and test r={4}:\n'
    print(header.format(n, round(n*p), k, t, r))
    check_1 = check_eq_2_1(n, p, k, r)
    print("Equation 2.1:", check_1)
    check_2 = check_eq_2_2(n, p, k, r)
    print("Equation 2.2:", check_2)
    check_3 = check_eq_2_3(n, p, k, r)
    print("Equation 2.3:", check_3)
    check_4 = check_eq_2_4(n, p, k, r)
    print("Equation 2.4:", check_4)
    check_5 = check_eq_2_5(n, p, k, r)
    print("Equation 2.5:", check_5)
    check_6 = check_eq_2_6(n, p, k, r, t)
    print("Equation 2.6:", check_6, "\n")
    if check_1 and check_2 and check_3 and check_4 and check_5 and check_6:
        return "Passed"
    return "Failed"