# Extracted from qAOPweaver/app/upstream/replicate_handling.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import pandas as pd

from app.upstream.profile import ReplicateConfig

_HIGH_CV_PCT = 30.0

_MAX_PAIRWISE_ROWS = 100_000

@dataclass
class AggregatedOutput:
    """
    Format A — one row per dose, columns: KE_mean, KE_sd, KE_n.
    Used by: POD workflow (pybmds), Empirical qKER R service, PSEM, BN, ODE.
    """
    df: pd.DataFrame
    # columns: dose | KE1_mean | KE1_sd | KE1_n | KE2_mean | ...
    quality_summary: dict

@dataclass
class PairedRowsOutput:
    """
    Format B — one row per paired observation.
    Rows within a single pairing are independent.
    Used by: ML single-run KER fitting, Empirical qKER bootstrap.
    """
    df: pd.DataFrame
    # columns: dose | KE1 | KE2 | ... | KEn
    pairing_info: dict

@dataclass
class BootstrapEnsembleOutput:
    """
    Format C — B paired DataFrames plus aggregated CI.
    Used by: ML pipeline (all_runs), uncertainty decomposition.
    """
    runs: list[pd.DataFrame]          # B DataFrames, each = one random pairing
    mean_ci: pd.DataFrame             # dose | KE_mean | KE_lower | KE_upper | KE_sd
    n_bootstrap: int
    seed: int

@dataclass
class ReplicateHandlingResult:
    """
    Complete output of _step2_replicate_handling().
    Always contains format_a; format_b/format_c depend on mode.
    """
    mode: str
    format_a: AggregatedOutput | None          # always populated
    format_b: PairedRowsOutput | None          # for modes 1, 2, 3, 4
    format_c: BootstrapEnsembleOutput | None   # for mode 3 only
    inflation_warning: str | None              # non-None for mode 2
    quality_flags: list[str] = field(default_factory=list)
    # One source-index matrix per stored run, captured during the actual pairing.
    source_runs: list[pd.DataFrame] = field(default_factory=list)

class ReplicateHandler:
    """Repeat the one-to-one pairing used by the demo."""

    def handle(
        self,
        df: pd.DataFrame,
        ke_cols: list[str],
        dose_col: str,
        rep_col: str,
        config: ReplicateConfig,
    ) -> ReplicateHandlingResult:
        if config.mode != "repeated_one_to_one":
            raise ValueError("The case replay uses repeated one-to-one pairing.")
        result = self._mode3_bootstrap(
            df, ke_cols, dose_col, rep_col, config.n_bootstrap,
            config.bootstrap_seed, "all_runs",
        )
        result.mode = "repeated_one_to_one"
        return result

    def _mode3_bootstrap(
        self,
        df: pd.DataFrame,
        ke_cols: list[str],
        dose_col: str,
        rep_col: str,
        n_bootstrap: int,
        seed: int,
        aggregation: str,  # "mean_ci" or "all_runs"
    ) -> ReplicateHandlingResult:
        """Repeat random one-to-one pairing B times; each iteration shuffles
        each KE's replicates within a dose group and pairs by position.
        Rows within one iteration ARE independent. Aggregate to mean + 95% CI.

        Leave-one-dose-out CV is applied OUTSIDE this module (ML/POD pipeline)."""
        rng = np.random.default_rng(seed)
        doses = sorted(df[dose_col].unique())
        dose_groups = {dose: df[df[dose_col] == dose] for dose in doses}
        source_indices = {
            dose: {ke: grp.index[grp[ke].notna()].to_numpy() for ke in ke_cols}
            for dose, grp in dose_groups.items()
        }
        positional_index = isinstance(df.index, pd.RangeIndex) and df.index.start == 0 and df.index.step == 1
        value_arrays = {ke: df[ke].to_numpy() for ke in ke_cols} if positional_index else {}
        all_runs: list[pd.DataFrame] = []
        source_runs = []
        run_means = {(dose, ke): [] for dose in doses for ke in ke_cols}

        for iteration in range(n_bootstrap):
            store_run = aggregation == "all_runs" or iteration == 0
            rows_b = []
            sources_b = []
            for dose in doses:
                ke_vals = {}
                ke_sources = {}
                for ke in ke_cols:
                    indices = source_indices[dose][ke].copy()
                    rng.shuffle(indices)
                    ke_sources[ke] = indices
                    ke_vals[ke] = value_arrays[ke][indices] if positional_index else df.loc[indices, ke].to_numpy()

                n_min = min(len(v) for v in ke_vals.values())
                if n_min:
                    for ke in ke_cols:
                        run_means[(dose, ke)].append(float(np.mean(ke_vals[ke][:n_min])))
                if store_run:
                    for i in range(n_min):
                        row = {dose_col: dose}
                        for ke in ke_cols:
                            row[ke] = ke_vals[ke][i]
                        rows_b.append(row)
                        sources_b.append({dose_col: dose, **{ke: (ke_sources[ke][i],) for ke in ke_cols}})

            if store_run:
                all_runs.append(pd.DataFrame(rows_b, columns=[dose_col] + ke_cols))
                source_runs.append(pd.DataFrame(sources_b, columns=[dose_col] + ke_cols))

        paired_doses = sorted(all_runs[0][dose_col].unique())

        # Aggregate across bootstrap runs → mean + 95% CI per (dose, KE).
        mean_ci_rows = []
        for dose in paired_doses:
            row = {dose_col: dose}
            for ke in ke_cols:
                vals_across_runs = run_means[(dose, ke)]
                row[f"{ke}_mean"] = float(np.mean(vals_across_runs))
                row[f"{ke}_lower"] = float(np.percentile(vals_across_runs, 2.5))
                row[f"{ke}_upper"] = float(np.percentile(vals_across_runs, 97.5))
                row[f"{ke}_sd"] = (
                    float(np.std(vals_across_runs, ddof=1))
                    if len(vals_across_runs) > 1 else 0.0
                )
            mean_ci_rows.append(row)
        mean_ci_columns = [dose_col] + [f"{ke}_{stat}" for ke in ke_cols for stat in ("mean", "lower", "upper", "sd")]
        mean_ci_df = pd.DataFrame(mean_ci_rows, columns=mean_ci_columns)

        # Format A: mean and SD from the bootstrap distribution.
        agg_rows = []
        quality: dict = {}
        for dose in paired_doses:
            row = {dose_col: dose}
            quality[dose] = {}
            dose_grp = dose_groups[dose]
            for ke in ke_cols:
                n_min = int(dose_grp[ke].dropna().shape[0])
                ke_mean = float(mean_ci_df.loc[mean_ci_df[dose_col] == dose, f"{ke}_mean"].values[0])
                ke_sd = float(mean_ci_df.loc[mean_ci_df[dose_col] == dose, f"{ke}_sd"].values[0])
                row[f"{ke}_mean"] = ke_mean
                row[f"{ke}_sd"] = ke_sd
                row[f"{ke}_n"] = n_min
                cv_pct = (ke_sd / abs(ke_mean) * 100) if ke_mean != 0 else np.nan
                quality[dose][ke] = {
                    "n": n_min, "mean": ke_mean, "sd": ke_sd, "cv_pct": cv_pct,
                    "flag": "high_cv" if (np.isfinite(cv_pct) and cv_pct > _HIGH_CV_PCT) else "ok",
                    "n_bootstrap": n_bootstrap,
                }
            agg_rows.append(row)
        agg_columns = [dose_col] + [f"{ke}_{stat}" for ke in ke_cols for stat in ("mean", "sd", "n")]
        agg_df = pd.DataFrame(agg_rows, columns=agg_columns).sort_values(dose_col).reset_index(drop=True)

        format_c = BootstrapEnsembleOutput(
            runs=all_runs if aggregation == "all_runs" else [],
            mean_ci=mean_ci_df,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )

        return ReplicateHandlingResult(
            mode="random_pairing_bootstrap",
            format_a=AggregatedOutput(df=agg_df, quality_summary=quality),
            format_b=PairedRowsOutput(
                df=all_runs[0].copy(),
                pairing_info={"method": "bootstrap_run_0", "seed": seed},
            ),
            format_c=format_c,
            inflation_warning=None,
            quality_flags=self._collect_quality_flags(quality),
            source_runs=source_runs,
        )


    def _collect_quality_flags(self, quality: dict) -> list[str]:
        """Scan quality_summary and collect any warning flags."""
        flags = []
        for dose, ke_dict in quality.items():
            for ke, metrics in ke_dict.items():
                flag = metrics.get("flag", "ok")
                if flag == "high_cv":
                    flags.append(
                        f"HIGH_CV: dose={dose} ke={ke} cv={metrics.get('cv_pct', 0):.1f}%"
                    )
                elif flag == "no_data":
                    flags.append(f"NO_DATA: dose={dose} ke={ke}")
        return flags
