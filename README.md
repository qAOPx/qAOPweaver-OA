# qAOPweaver-OpenAccess

R and Python reproduction code for the qAOPweaver demo case.

## Install the dependencies

The verified setup uses 64-bit Windows Python 3.12.7 and R 4.6.0 for preparation, ML and the four R workflows. POD uses Ubuntu x86-64 and Python 3.12.16 through WSL. In a fresh Windows virtual environment, install the pinned public packages.

~~~powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-windows.txt
~~~

Use a clean R 4.6.0 installation. To restore the packages into a separate library, create the library directory and install `renv`. From the repository root, run `renv::restore(project = ".", library = "public-r-library", prompt = FALSE)`. Pass `--r-library public-r-library` to the runner when using that library.

In Ubuntu, create a Python 3.12.16 virtual environment and install `python/fits/requirements-linux-pod.lock`. Use the absolute Linux path to that environment's Python executable.

## Run the demo

Use the `demo-replay-1.0` tag with the [Zenodo data archive](https://doi.org/10.5281/zenodo.24107875). The runner checks the ZIP and its member checksums, then rebuilds the inputs and reruns all six workflows.

~~~powershell
.\.venv\Scripts\python.exe replay.py --archive "C:\path\demo-data.zip" --fit-only --out "demo-run" --rscript "C:\path\Rscript.exe" --pod-wsl-python "/home/you/.venv-pod/bin/python"
~~~

For the optional 487-value comparison, replace `--fit-only` with `--si "C:\path\reference-tables.docx"` and supply the matching reference tables as a separate Word document.

Choose a new output directory for each run. On a system that can run both interpreters natively, `--pod-python` selects a separate Python executable for POD. Use the Ubuntu environment for the reference POD estimates. The optional comparison reports Windows numerical differences and exits with code 2 if a reference value differs at its displayed precision.

The `runs` folder contains the fitted results. `validation.json` records the ZIP checksum and whether the optional comparison ran. The optional comparison also records the reference document checksum, counts and individual checks in `comparison`. `qker-ode-archive-comparison.json` compares the qKER candidates and the synthetic and measured ODE results with the archived runs.

## Code scope

| Path | Case calculations |
| --- | --- |
| python/preparation | Read the 360 measured rows and the AOP 263 path, preserve the recorded effective observations and construct the three seeded pairings and their source provenance. Redraw disjoint training and test cases for each of the five treatment-held-out folds. |
| python/fits | Fit continuous POD at a one-SD benchmark response with AIC-weighted model averaging. Fit the case's Gaussian process and SciPy candidate formula analyses and calculate the associated prediction curves and holdout diagnostics. |
| r/qker_ode | Rebuild qKER inputs, fit and rank the selected model forms, and rerun the archived candidate set. Recalculate the 1000 nM synthetic and exploratory measured ODE fits using the recorded RK4 and bounded optimization settings. |
| r/prepared | Fit the three weighted tPSEM realizations and score their held-out conditions. Fit the DBN CPTs and method comparisons, calculate posterior coefficient summaries, and query the fitted networks. |
| replay.py and compare_results.py | Verify the ZIP, run the fits, and check the results against the archived runs and optional reference tables. |

The module documentation records the source revisions and adaptations. See [validation results](VALIDATION.md) for the checks and observed platform differences.

The code uses the [MIT License](LICENSE). For data and model terms, see `LICENSES.md` in the Zenodo archive.
