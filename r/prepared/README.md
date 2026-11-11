# Demo tPSEM and DBN

```sh
Rscript --vanilla r/prepared/run.R prepared_cases.csv temporal_edges.csv fold_cases.csv output
```

The case table contains the three 30-case realizations at 2, 2, 4, and 7 days. The edge table identifies the three selected relationships in path order. The fold table adds `held_out_condition_id` and `partition` to the case table fields. Generate training and test cases separately from disjoint source observations for each treated concentration. The 0 nM control remains in the training set for every assessment.

The runner fits three weighted tPSEM models with condition-cluster inference and scores the five condition-held-out folds for each realization. It fits DBN conditional probability tables by posterior predictive simulation, observed counts, and fixed-coefficient simulation with residual variation. The first realization sets the shared state cutpoints for each configuration. The baseline query and low-state query use a mixture of complete fitted networks. The demo query conditions on low day-2 KE1446 and returns marginals for all four nodes.

`prepared_fit_core.R` contains selected functions from `qAOPweaver/qAOPweaverR/R/prepared_analysis_v1.R`. `bn_query_core.R` contains selected functions from `qAOPweaver/qAOPweaverR/R/bn_static_cases.R`. Both source files are from qAOPweaver 0.29.0, revision `0143dd56dd21a5268b17ecbe8383647b1c6323c9`. The case runner sources these files directly. Prediction, state cutpoints, and root priors use the settings specified for this case. Its public R dependencies are `jsonlite`, `gRain`, and `gRbase` (plus their dependencies). The project dependency manifest records their versions.
