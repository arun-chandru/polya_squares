#!/usr/bin/env python3
"""Exact finite certificates for the two Pólya theorems in the manuscript.

Python 3.8+; standard library only; no input files, network, or floating point.
Run: python verify_polya.py [--output certificate_results.json]

The program verifies the finite statements in manuscript.tex, Sections
3--5 and Appendix B. The mathematical interval reductions and infinite
tails are proved there; this program checks their explicit numerical
constants, not the variational or analytic arguments. Internal comments
label the arithmetic blocks C0--C7 for ease of source review.
The only optional file write is the explicitly named --output destination.
Assertions are deliberately avoided: verification also runs under python -O.
"""

import argparse
from array import array
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json


CAP = 4_000_000
PI_LO, PI_HI = F(333, 106), F(355, 113)
Q_LO, Q_HI = F(97, 200), F(11, 20)
R = F(243, 500)
EXCEPTIONS = (2, 5, 6, 10, 13, 20, 26, 29, 41, 53, 74, 85, 90, 130)

# C5: static candidate certificate. Columns: n, q_lo, q_hi, e_lo, e_hi,
# comparison (S=subtraction, M=mixed), claimed integer upper count.
# Every count, Weyl margin, and union-coverage claim is recomputed below.
RECTANGLES = (
    (2,   '97/200', '11/20',   '2',     '3',     'S', 1),
    (5,   '97/200', '11/20',   '5',     '6',     'M', 2),
    (6,   '97/200', '11/20',   '6',     '7',     'S', 3),
    (10,  '97/200', '11/20',   '10',    '11',    'S', 5),
    (13,  '97/200', '11/20',   '13',    '14',    'S', 7),
    (26,  '97/200', '11/20',   '26',    '27',    'S', 14),
    (29,  '97/200', '207/400', '29',    '30',    'S', 16),
    (29,  '207/400','11/20',   '29',    '30',    'M', 14),
    (41,  '97/200', '11/20',   '41',    '42',    'M', 22),
    (53,  '97/200', '207/400', '53',    '54',    'M', 30),
    (53,  '207/400','11/20',   '53',    '54',    'S', 29),
    (74,  '97/200', '207/400', '74',    '75',    'S', 42),
    (74,  '207/400','11/20',   '74',    '149/2', 'M', 40),
    (74,  '207/400','427/800', '149/2', '75',    'S', 41),
    (74,  '427/800','11/20',   '149/2', '75',    'S', 39),
    (85,  '97/200', '207/400', '85',    '86',    'M', 48),
    (85,  '207/400','11/20',   '85',    '86',    'M', 46),
    (90,  '97/200', '207/400', '90',    '91',    'S', 51),
    (90,  '207/400','11/20',   '90',    '91',    'M', 48),
    (130, '97/200', '207/400', '130',   '131',   'M', 74),
    (130, '207/400','11/20',   '130',   '131',   'M', 70),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def atan_enclosure(q, terms):
    """Alternating arctangent series at 1/q, including the next term."""
    require(q > 1 and terms > 0, 'Invalid arctangent parameters')
    s = sum((F((-1)**k, (2*k+1)*q**(2*k+1))
             for k in range(terms)), F(0))
    t = s + F((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(s, t), max(s, t)


def certify_constants():
    a_lo, a_hi = atan_enclosure(5, 12)
    b_lo, b_hi = atan_enclosure(239, 4)
    machin_lo, machin_hi = 16*a_lo-4*b_hi, 16*a_hi-4*b_lo
    require(PI_LO < machin_lo < machin_hi < PI_HI, 'C0: pi enclosure')
    # tan(4 atan(1/5) - atan(1/239)) = 1 exactly.  Angle selection
    # in (0, pi/2), hence equality to pi/4, is justified in the paper.
    tangent = F(1, 5)
    for _ in range(2):
        tangent = 2*tangent/(1-tangent*tangent)
    require((tangent-F(1, 239))/(1+tangent/F(239)) == 1,
            'C0: Machin tangent identity')
    require(4*a_lo-b_hi > 0 and 4*a_hi-b_lo < 1,
            'C0: Machin angle lies between zero and one')

    s2_lo, s2_hi = 1-F(12, 5)/PI_LO, 1-F(12, 5)/PI_HI
    require(Q_LO**2 < s2_lo < s2_hi < R**2, 'C1: side bounds')
    transitions = []
    for k, left, right in ((2, 8, 9), (5, 21, 22), (18, 76, 77)):
        require(left < F(k)/s2_hi <= F(k)/s2_lo < right,
                'C1: envelope transition')
        transitions.append({'numerator': k, 'strict_lower': left,
                            'strict_upper': right})
    require(PI_HI**2*F(2, 16) < F(1309, 1000)**2 < F(4, 3)**2,
            'C1: increasing envelope slopes')
    require(0 < s2_hi < F(2, 5), 'C1: sharpness witness interval')

    # C7: constants used in the analytic discrepancy proof.
    require(F(16, 45) < F(9, 25) and F(18, 25) < 1,
            'C7: elementary Bessel bounds on [1,2]')
    require(F(32, 13)/PI_LO < 1, 'C7: Bessel energy constant')
    require(F(157, 50) < PI_LO and PI_HI < F(22, 7),
            'C7: coarser pi bounds')
    require(F(1, 16) < F(2, 5)**3, 'C7: 2^(-4/3) < 2/5')
    require(F(19, 12)**3 < 4, 'C7: 2^(2/3) > 19/12')
    require(F(12, 5) < F(31, 20)**2, 'C7: square-root majorant')
    discrepancy_constant = F(121, 35)+F(22320, 2983)
    require(discrepancy_constant == F(1142143, 104405) < 11,
            'C7: discrepancy coefficient')
    lower1000 = 1-F(3, 4000)-F(11, 40)
    upper1000 = 1+F(1, 4000)+F(11, 40)
    require(10**3 == 1000 and lower1000 == F(2897, 4000)
            and lower1000 > F(81, 125) and upper1000 == F(5101, 4000)
            and upper1000 < F(4, 3), 'C7: tail at radius 1000')
    lower2000 = 1-F(3, 8000)-F(11, 48)
    upper2000 = 1+F(1, 8000)+F(11, 48)
    require(12**3 < 2000 and lower2000 == F(18491, 24000)
            and lower2000 > F(11, 15) and upper2000 == F(29503, 24000)
            and upper2000 < F(4, 3), 'C7: tail at radius 2000')
    require(F(8, 27) < F(9, 25), 'C7: 2 sqrt(2/27) < 3/5')
    require(F(3, 5)+F(8, 3)*(1-Q_HI) == F(9, 5)
            and 1-2*Q_HI == -F(1, 10) and 18**2 == 324 < 400,
            'C7: thin-frame tail')
    return {
        'pi_lower': str(PI_LO), 'pi_upper': str(PI_HI),
        'machin_lower': str(machin_lo), 'machin_upper': str(machin_hi),
        'side_threshold': 'sqrt(1-12/(5*pi))',
        'side_lower': str(Q_LO), 'side_upper': str(R),
        'squared_transition_intervals': transitions,
        'discrepancy_constant_upper': str(discrepancy_constant),
        'radius1000_lower_slope': str(lower1000),
        'radius1000_upper_slope': str(upper1000),
        'radius2000_lower_slope': str(lower2000),
        'radius2000_upper_slope': str(upper2000),
        'thin_tail_start': 324,
    }


def build_square_counts():
    """C_n = number of ordered positive pairs a^2+b^2 <= n.

    Unordered pairs are weighted by 2 off the diagonal and 1 on it.
    The single array occupies 4*(CAP+1) bytes on normal Python platforms;
    all subsequent products are unbounded Python integers.
    """
    counts = array('I', [0])*(CAP+1)
    require(counts.itemsize >= 4, 'Need an unsigned array type of >=32 bits')
    for a in range(1, isqrt(CAP//2)+1):
        a2 = a*a
        for b in range(a, isqrt(CAP-a2)+1):
            counts[a2+b*b] += 1 if a == b else 2
    for n in range(1, CAP+1):
        counts[n] += counts[n-1]
    require(counts[-1] == 3_139_586 and counts[1_000_000] == 784_387,
            'C2: count checkpoints')
    return counts


def sqrt_floor(r):
    require(r >= 0, 'Square-root radicand must be nonnegative')
    return isqrt(r.numerator//r.denominator)


def square_count(r):
    """Independent row count at rational squared radius r."""
    return sum(sqrt_floor(r-a*a) for a in range(1, sqrt_floor(r)+1))


def certify_deficits(counts):
    """C2: inclusive lower bounds and left-limit upper bounds."""
    exceptions, sharp_middle_exceptions = [], []
    for n in range(1, CAP+1):
        upper = 355*n-452*counts[n-1]
        require(upper <= 0 or 9*upper*upper <= 1808**2*n,
                'C2: global upper deficit at n=%s' % n)
        if n <= 18:
            require(upper <= 0 or (1000*upper)**2 <= (1309*452)**2*n,
                    'C2: prefix upper deficit at n=%s' % n)
        lower = 333*n-424*counts[n]
        if lower < 0 or (15*lower)**2 < (11*424)**2*n:
            exceptions.append(n)
        if 21 <= n <= 76:
            if lower < 0 or (500000*lower)**2 < (318087*424)**2*n:
                sharp_middle_exceptions.append(n)
        if 76 <= n <= 1_000_000:
            require(lower >= 0 and (125*lower)**2 >= (81*424)**2*n,
                    'C2: sharp high deficit at n=%s' % n)
    require(tuple(exceptions) == EXCEPTIONS, 'C2: exception bands')
    require(sharp_middle_exceptions == [53], 'C2: sharp middle exception')
    for n in range(9):
        require(3*n >= 5*counts[n], 'C2: first quadratic branch')
    for n in range(8, 22):
        lower = 333*n-424*counts[n]
        require(3*n >= 5*(counts[n]-1), 'C2: second quadratic branch')
        require(lower >= 0 and (500*lower)**2 >= 2*n*(243*333)**2,
                'C2: intermediate linear branch')
    # An independent, differently organized enumeration checks the entire
    # low-energy table and checkpoints near both analytic tail joins.
    for n in list(range(132))+[999999, 1000000, 3999999, CAP]:
        require(counts[n] == square_count(F(n)), 'C2: independent row count')
    return {
        'cap_squared_radius': CAP, 'count_at_cap': counts[-1],
        'count_at_1000000': counts[1_000_000],
        'global_upper_left_limit_failures': [],
        'prefix_upper_left_limit_failures': [],
        'middle_frame_exception_bands': exceptions,
        'sharp_envelope_exception_bands': sharp_middle_exceptions,
        'sharp_high_lower_failures': [], 'low_branch_failures': [],
        'independent_row_count_checks': 136,
        'low_counts_n_0_through_131': list(counts[:132]),
    }


def certify_partition(counts):
    """C3: the finite endpoint constants in the aggregate-hole lemma."""
    budget, ceiling = F(18, 5), F(9, 2)
    require(R*R*54 < budget*budget < 13, 'C3: exceptional-band budget')
    require(F(3, 2)**2 > 2, 'C3: radical upper bound')
    small = PI_HI*F(3, 2)*budget/4
    require(small == F(1917, 452) < ceiling, 'C3: all-small partition')
    endpoints = []
    # The first endpoint is inclusive sqrt(2); the others are left limits.
    for k, root_lo, count_before, expected in (
        (2, F(7, 5), 1, F(2488, 565)),
        (5, F(223, 100), 1, F(49973, 11300)),
        (8, F(14, 5), 3, F(442, 113)),
        (10, F(79, 25), 4, F(11349, 2825)),
    ):
        require(0 < root_lo and root_lo**2 < k, 'C3: radical lower bound')
        actual_count = counts[k] if k == 2 else counts[k-1]
        require(actual_count == count_before, 'C3: endpoint count')
        bound = PI_HI*(2*k+budget**2-2*budget*root_lo)/4-count_before
        require(bound == expected < ceiling, 'C3: partition endpoint')
        endpoints.append({'squared_radius': k, 'radical_lower': str(root_lo),
                          'count': count_before, 'upper': str(bound)})
    require(counts[12] == 6 and counts[53] == 37, 'C3: final counts')
    last = PI_HI*budget**2/4-6
    outer = 53*PI_LO/4-counts[53]
    require(last == F(2361, 565) < ceiling < outer == F(37, 8),
            'C3: partition closing inequalities')
    return {'budget': str(budget), 'deficit_ceiling': str(ceiling),
            'all_small_upper': str(small), 'endpoint_bounds': endpoints,
            'last_endpoint_upper': str(last), 'outer_53_lower': str(outer)}


def mixed_upper(qlo, qhi, ehi):
    """Inclusive count of modes with separately minimized coefficients.

    Summing the b ranges by isqrt avoids the original double mode loop.
    In the vertical term +1 retains b=0, including a zero radicand.
    """
    horizontal = vertical = 0
    for a in range(1, sqrt_floor(ehi)+1):
        r = (1-qlo)**2*(ehi-a*a)
        if r >= 1:
            horizontal += (sqrt_floor(r)+1)//2
        r = qhi*qhi*(ehi-4*a*a/(1-qlo)**2)
        if r >= 0:
            vertical += sqrt_floor(r)+1
    return 2*(horizontal+vertical)


def certify_operator_transfer(counts):
    """C4: numerical input to the analytic min--max transfer argument."""
    half_count = max(counts[n]-counts[n//4] for n in range(25))
    coefficient = Q_LO/(4*(1-Q_LO)**3)
    weyl = PI_LO*(1-Q_HI**2)*20/4
    require(half_count == 10, 'C4: half-frame count below energy 25')
    require(coefficient == F(970000, 1092727) and 25*coefficient > 22,
            'C4: transfer coefficient')
    require(weyl == F(92907, 8480) > 10, 'C4: energy-20 Weyl lower bound')
    q, e = F(49, 100), F(20)
    subtraction = square_count(e)-square_count(q*q*e)
    mixed = mixed_upper(q, q, e)
    weyl_upper = PI_HI*(1-q*q)*e/4
    require(subtraction == mixed == 12 and weyl_upper < 12,
            'C4: raw comparison obstruction')
    return {'half_frame_count_below_25_max': half_count,
            'minimum_transfer_coefficient': str(coefficient),
            'lambda11_over_pi_squared_lower': str(25*coefficient),
            'energy20_weyl_lower': str(weyl),
            'raw_comparison_obstruction': {'q': str(q), 'squared_energy': 20,
                'subtraction': subtraction, 'mixed': mixed,
                'weyl_upper': str(weyl_upper)}}


def atomic_representatives(edges, include_last):
    """Each axis endpoint plus one representative of each intervening cell."""
    return (edges if include_last else edges[:-1]) + [
        (a+b)/2 for a, b in zip(edges, edges[1:])]


def certify_middle_cover(counts):
    """C5: validate the static rectangles as an exact continuum cover."""
    bands, result = {}, []
    min_slack = None
    for n, ql, qh, el, eh, method, claimed in RECTANGLES:
        qlo, qhi, elo, ehi = map(F, (ql, qh, el, eh))
        require(n in EXCEPTIONS and n != 20, 'C5: extraneous band')
        require(Q_LO <= qlo < qhi <= Q_HI and n <= elo < ehi <= n+1,
                'C5: rectangle endpoints')
        if method == 'S':
            inner = qlo*qlo*elo
            upper = counts[n]-counts[inner.numerator//inner.denominator]
        else:
            require(method == 'M', 'C5: unknown comparison')
            upper = mixed_upper(qlo, qhi, ehi)
        weyl = PI_LO*(1-qhi*qhi)*elo/4
        require(upper == claimed and upper <= weyl, 'C5: rectangle inequality')
        slack = weyl-upper
        min_slack = slack if min_slack is None else min(min_slack, slack)
        bands.setdefault(n, []).append((qlo, qhi, elo, ehi))
        result.append({'band': n, 'qlo': ql, 'qhi': qh, 'energy_lo': el,
                       'energy_hi': eh, 'method': method,
                       'count_upper': upper, 'weyl_lower': str(weyl)})
    require(set(bands) == set(EXCEPTIONS)-{20}, 'C5: missing band')
    probes = 0
    for n, rects in sorted(bands.items()):
        qs = sorted({Q_LO, Q_HI} | {v for r in rects for v in r[:2]})
        es = sorted({F(n), F(n+1)} | {v for r in rects for v in r[2:]})
        for q in atomic_representatives(qs, True):
            for e in atomic_representatives(es, False):
                require(any(ql <= q <= qh and el <= e < eh
                            for ql, qh, el, eh in rects), 'C5: uncovered cell')
                probes += 1
    require(len(result) == 21 and probes == 120 and min_slack == F(107, 3200),
            'C5: cover checkpoints')
    return {'rectangles': result, 'rectangle_count': len(result),
            'atomic_cell_and_endpoint_checks': probes,
            'minimum_weyl_minus_count': str(min_slack)}


def certify_thin_frames():
    """C6: integer interval arithmetic on 1 <= x < 400 at q=11/20.

    This deterministic bisection is a proof-producing interval check; no
    sampled parameter values or floating-point optimization occur. Every
    initial half-open interval is partitioned into verified covering leaves.
    """
    grid, root_den = 4096, 1 << 20
    require(root_den*root_den//grid == 1 << 28, 'C6: root scale')
    accepted = tested = max_depth = 0
    min_margin = None
    total_width = 0
    for j in range(1, 400):
        initial_lo, initial_hi = j*j*grid, (j+1)*(j+1)*grid
        stack = [(initial_lo, initial_hi, 0)]
        squares = [k*k*grid for k in range(1, j+1)]
        previous_right = initial_lo
        while stack:
            lo, hi, depth = stack.pop()
            tested += 1
            max_depth = max(max_depth, depth)
            odd = even = 0
            for k, square in enumerate(squares, 1):
                upper = isqrt((hi-square) << 28)+1
                if k & 1:
                    odd += upper
                else:
                    even += upper
            margin = (103230*lo*root_den
                      - 8480*grid*(20*odd+11*even+9*(j//2)*root_den))
            if margin >= 0:
                # Leaves are visited left to right; these checks explicitly
                # verify no gaps, overlap, or lost endpoint in the partition.
                require(lo == previous_right and lo < hi, 'C6: leaf coverage')
                previous_right = hi
                total_width += hi-lo
                accepted += 1
                min_margin = margin if min_margin is None else min(min_margin, margin)
            else:
                require(hi-lo > 1, 'C6: uncertified unit-width box')
                mid = (lo+hi)//2
                stack.append((mid, hi, depth+1))
                stack.append((lo, mid, depth+1))
        require(previous_right == initial_hi, 'C6: terminal endpoint coverage')
    require(total_width == (400**2-1)*grid, 'C6: total covering width')
    require((tested, accepted, max_depth, min_margin)
            == (25837, 13118, 7, 169637314560), 'C6: scalar checkpoints')
    return {'q': str(Q_HI), 'x_interval': '[1,400)',
            'grid_denominator': grid, 'root_upper_denominator': root_den,
            'initial_intervals': 399, 'boxes_tested': tested,
            'covering_boxes': accepted, 'maximum_depth': max_depth,
            'failed_boxes': 0, 'least_integer_margin': min_margin,
            'covered_squared_energy_width_numerator': total_width,
            'analytic_tail_start': 324,
            'all_q_extension': 'Affine interpolation between q=11/20 and q=1'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Explicit JSON destination; omit to write no file')
    args = parser.parse_args()
    constants = certify_constants()
    counts = build_square_counts()
    result = {
        'certificate_version': 1,
        'arithmetic': 'Python integers, integer isqrt, and fractions.Fraction',
        'C0_C1_C7_constants': constants,
        'C2_deficits': certify_deficits(counts),
        'C3_partition': certify_partition(counts),
        'C4_operator_transfer_constants': certify_operator_transfer(counts),
        'C5_middle_rectangle_cover': certify_middle_cover(counts),
        'C6_thin_frame_intervals': certify_thin_frames(),
        'all_checks_passed': True,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output is not None:
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded, end='')


if __name__ == '__main__':
    main()
