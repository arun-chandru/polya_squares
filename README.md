# Pólya inequalities for square frames and perforated squares

**Author:** R. Arun Chandru  
**Date:** 27 September 2026

## Abstract

We prove Pólya's inequality for the Dirichlet eigenvalues of every concentric square frame whose inner and outer sides are parallel, for every hole-to-outer-side ratio strictly between zero and one and every eigenvalue. We also prove a position- and orientation-independent result for finitely many square holes: the sum of their side lengths, divided by the outer side length, may be as large as `sqrt(1 - 12/(5*pi))`. This constant is optimal for the universal square-count subtraction criterion, not for the geometric inequality. A further mixed-boundary comparison treats a single arbitrarily displaced, strictly interior parallel rectangular hole when both its side lengths are at least eleven twelfths of the outer square's side length. Explicit analytic estimates cover all unbounded spectral ranges. Integer and rational certificates verify the remaining compact comparisons without numerically approximating eigenvalues of a perforated domain.

## Scope

The operator throughout is the ordinary Euclidean Dirichlet Laplacian. Each result gives `area(Omega) * lambda_k(Omega) >= 4*pi*k` for every positive integer `k`, with multiplicities.

1. **Centered square frame:** every side ratio `0 < q < 1`, parallel inner and outer squares.
2. **Several small square holes:** finitely many closed, pairwise disjoint squares strictly inside an outer square of side `L`, with arbitrary positions and orientations, provided `sum(l_i) <= sqrt(1 - 12/(5*pi))*L`.
3. **One large rectangular hole:** a closed rectangle strictly inside an outer square, with parallel sides and both relative side lengths in `[11/12, 1)`. Its center is arbitrary subject to containment. The two gap totals need not be comparable.

The sharpness statement concerns the scalar comparison in item 2 only. The paper does **not** prove Pólya's conjecture for arbitrary planar domains, all displaced intermediate-sized holes, large rotated holes, general outer rectangles, or Neumann eigenvalues. The constant `11/12` is sufficient and is not asserted optimal. The geometric assertions include arbitrarily small positive boundary gaps, not actual boundary contact.

## Package contents

| File | Purpose |
|---|---|
| `manuscript.pdf` | Complete paper, including the analytic and exact-certificate appendices |
| `manuscript.tex` | Single standalone LaTeX source; bibliography included |
| `verify_polya.py` | Exact certificates for the centered-frame and finite-square-hole results |
| `verify_rectangular.py` | Exact certificates for the large rectangular-hole result |
| `certificate_results.json` | Deterministic output of the first verifier |
| `rectangular_results.json` | Deterministic output of the second verifier |
| `README.md` | This description and reproduction instructions |

The JSON files are output records, not inputs or trusted premises. Neither verifier reads them.

## Reproduce the exact computations

Use Python 3.8 or later and its standard library. From the folder containing the scripts:

```text
python verify_polya.py
python verify_rectangular.py
```

To save independent output records:

```text
python verify_polya.py --output checked_polya.json
python verify_rectangular.py --output checked_rectangular.json
```

The programs use arbitrary-precision integers, exact rational arithmetic, and integer square roots. No floating-point decision enters a certificate. Their explicit failure checks remain active with `python -O`. The first program uses an unsigned-integer array of 4,000,001 entries; on the tested platform that array occupies about 16 MB. The second program uses bounded memoization. Neither program demands network access or a GPU.

The finite computations are parts of the proof, not substitutes for the analytic tail estimates and variational comparisons supplied in the manuscript. Each script terminates with an `all_checks_passed` result on successful verification and raises an exception if a required inequality or coverage check fails.

## Compile the manuscript

The source uses standard LaTeX packages and the `amsart` class. It has no external figures, bibliography database, or shell-escape requirement. With a current TeX distribution:

```text
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

Alternatively, `latexmk -pdf manuscript.tex` or `tectonic manuscript.tex` manages the required reruns.

