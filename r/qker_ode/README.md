# Demo qKER and ODE calculations

Run from the repository root with R 4.6.0.

    Rscript --vanilla r/qker_ode/run.R EXTRACTED_DEMO_ARCHIVE RESULTS_DIRECTORY all

The last argument can be set to `qker` or `ode`. The input directory is the extracted Zenodo archive. The launcher recomputes concentration-by-time condition means from the measured 360-row CSV and the synthetic 324-row CSV. qKER uses the means at days 2, 4, and 7 and the three AOP263 edges. ODE uses the 1,000 nM series and the saved model, time grid, bounds, and fitting configuration for each dataset.

The only non-base packages used directly are jsonlite 2.0.0 for data exchange and mgcv 1.9-4 for the v0320 GAM candidate. The nine GAM attempts reproduce the archived fitting failures because each pair has six distinct covariate points. The other 27 candidates are fitted and exported. The two ODE fits use the original fixed-step RK4 loop and bounded L-BFGS-B objective. Their parameter starts follow the archived bound-midpoint rule.

`upstream.R` contains selected function definitions from qAOPweaver 0.29.0 at `0143dd56dd21a5268b17ecbe8383647b1c6323c9`. `origins.tsv` records the source paths and line ranges. `run.R` contains the local case adapters for input preparation, per-time qKER ranking, the fixed ODE optimizer objective, and output assembly. The integrator preserves the original RK4 branch. The original qKER input dispatcher, R fit route, and ODE inverse fitter provide the behavior used in these adapters. `settings.json` holds the analysis settings from the archived requests with session and job fields omitted.

`verify.py` is an additional fixture comparison. It checks all 12 ranked qKER combinations, all 36 v0320 candidate rows, and both ODE fits against the archived Windows run. Run it with the full extracted archive and the new results directory.

    python r/qker_ode/verify.py EXTRACTED_DEMO_ARCHIVE RESULTS_DIRECTORY

The main runner can also compare the new qKER and synthetic ODE values with separate reference tables.
