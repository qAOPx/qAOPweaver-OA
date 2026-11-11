# Demo validation results

The runner verified the Zenodo ZIP with SHA-256 `5ee65e07c019e62889edeb355374c7f82896988b9280abd084b7eacfe67e4f66` and all 278 manifest entries. The optional comparison used a separate Word reference document and recorded its SHA-256 in `validation.json`. The checks ran with Windows Python 3.12.7, an isolated R 4.6.0 installation, and Ubuntu WSL Python 3.12.16 for POD. The environments used the recorded public dependencies.

| Check | Matches at reference precision |
| --- | ---: |
| Input summary | 216/216 |
| POD | 36/36 |
| Empirical qKER | 54/54 |
| tPSEM | 48/48 |
| DBN | 93/93 |
| ML | 21/21 |
| Synthetic ODE | 19/19 |
| Total | 487/487 |

The runner reproduced the 360 effective measured observations, the 90 weighted paired cases, the three temporal edges, and the 360 source tuples and ledger records. It regenerated 450 cases across the five disjoint treatment-held-out folds. The largest absolute discrepancy in the tPSEM fold RMSE and MAE was 1.10e-13 percentage points. DBN direct inference agreed with independent enumeration within 6.39e-16 in probability.

All 72 POD concentration summaries matched the hosted exports. Ubuntu reproduced all 36 BMD, BMDL and BMDU reference values at displayed precision. An independent Windows POD check matched 17/36 at the same precision, consistent with the archived Windows comparison. The two operating systems use different dependency sets and both emitted the native optimizer warning. Use the Ubuntu lock for the reference comparison.

All 12 ranked qKER combinations and all 36 candidate records matched the archived local run. Nine GAM candidate records reported the original fit failure for too few unique covariate combinations. Both ODE fits converged. Fit metrics agreed with the archived local runs at their stored precision, with maximum differences of about 0.0001 in the archived fitted parameters, predictions or residuals. Both ODE fits were compared with their archived results. The synthetic fit was also compared with the reference values at displayed precision.

The comparison writes individual checks and numerical diagnostics to the `comparison` folder. It measures agreement for this demo under the recorded settings.
