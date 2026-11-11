# Python demo fits

`run_fits.py --observations PATH --output RUNS` reads the 360 effective measured
observations from the demo data archive. `--only pod` or `--only ml` selects
one workflow for a separate dependency environment. Output paths and scientific
fields match those used by the optional reference comparison.

- POD groups each endpoint and sampling day by concentration, then runs the
  original continuous analyzer with 1 SD BMR and AIC-weighted model averaging.
  `pod_continuous.py` restricts the analyzer to the case settings from
  `qAOPweaver/app/pod/continuous.py` at revision
  `0143dd56dd21a5268b17ecbe8383647b1c6323c9` (upstream file SHA-256
  `91941d70cfb962905981d0913b9c4758a14cddae221576667b6a54d1e6b1315e`).
  Grouping and fit invocation in `run_fits.py` follow the continuous branch of
  `app/pod/batch.py` at the same revision.
- ML forms separate means within each endpoint, concentration, and day and
  fits the three adjacent KERs on 18 pooled conditions. `ml_case.py` adapts
  the calculations in `app/aiml/gp.py`, the SciPy candidate-formula branch in
  `app/aiml/symbolic.py`, and the held-concentration-out diagnostics in
  `app/aiml/pipeline.py` and `app/aiml/data/validation.py` from that revision.
  It emits both prediction bands, fitted data metrics, the selected SciPy
  formula per link, and the MIE-to-AO holdout diagnostics. The near-constant
  control responses make the control-fold R² unstable.

The Python 3.12.7 Windows lock reproduces the local ML environment.
The Python 3.12.16 Ubuntu x86-64 lock for POD uses the numerical package
versions recorded for the Linux POD replay. Install either lock into a
fresh virtual environment. The fit code runs directly with these public
dependencies and does not import the application package.

The Ubuntu POD test rebuilt all 72 concentration summaries from the original
measured CSV, matched all 72 hosted groups exactly, and matched the 36 reference values
at displayed precision. The largest absolute difference from the
hosted BMD, BMDL, or BMDU was 2.3e-13 nM. The Windows test of the same analyzer
matched 17 of 36 at printed precision, consistent with the archived
Windows comparison. Both environments emitted the native optimizer warning
about the disabled NLopt Luksan code. Its presence alone did not distinguish
the matching and differing runs.

The Windows ML test using the archived effective observations matched all 21
reference values. All 600 grid points and prediction limits matched the archived
Windows ML refits exactly. Maximum absolute curve differences from the hosted
exports were 2.2e-6 percentage points for GP and 1.8e-5 for SciPy.
