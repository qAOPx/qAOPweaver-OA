# Extracted from qAOPweaver/app/upstream/preprocess.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

from typing import Any, Literal

import numpy as np

import pandas as pd

from app.core.contracts import (
    ActionLogEntry,
    PreprocessSummary,
)

MissingStrategy = Literal["retain_na"]
ReplicateStrategy = Literal["retain_replicates"]
NormStrategy = Literal["none"]
HarmonizeStrategy = Literal["keep_original"]


_DOSE_FACTORS_TO_M: dict[str, float] = {
    "m": 1.0, "mol/l": 1.0, "moll-1": 1.0, "moll^-1": 1.0,
    "mm": 1e-3, "mmol/l": 1e-3, "mmoll-1": 1e-3, "mmoll^-1": 1e-3,
    "um": 1e-6, "μm": 1e-6, "umol/l": 1e-6, "μmol/l": 1e-6,
    "nm": 1e-9, "nmol/l": 1e-9, "nmoll-1": 1e-9,
    "pm": 1e-12, "pmol/l": 1e-12, "pmoll-1": 1e-12,
}

_TIME_FACTORS_TO_H: dict[str, float] = {
    "h": 1.0, "hr": 1.0, "hrs": 1.0, "hour": 1.0, "hours": 1.0,
    "d": 24.0, "day": 24.0, "days": 24.0,
    "min": 1 / 60, "mins": 1 / 60, "minute": 1 / 60, "minutes": 1 / 60,
    "s": 1 / 3600, "sec": 1 / 3600, "second": 1 / 3600, "seconds": 1 / 3600,
    "week": 168.0, "weeks": 168.0,
}

def _normalize_unit(u: str) -> str:
    u = str(u).lower().strip()
    u = u.replace("µ", "u").replace("μ", "u").replace("·", " ")
    u = " ".join(u.split())
    return u.replace(" ", "")

def _dose_factor(unit: str | None) -> float | None:
    if not unit:
        return None
    return _DOSE_FACTORS_TO_M.get(_normalize_unit(str(unit)))

def _time_factor(unit: str | None) -> float | None:
    if not unit:
        return None
    return _TIME_FACTORS_TO_H.get(_normalize_unit(str(unit)))

def _standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure canonical analysis columns exist; create them from aliases if needed."""
    out = df.copy()

    # endpoint_id ← ke_id (preferred) or endpoint
    if "endpoint_id" not in out.columns:
        for src in ("ke_id", "endpoint"):
            if src in out.columns:
                out["endpoint_id"] = out[src].astype(str)
                break

    # time ← time_point
    if "time" not in out.columns:
        for src in ("time_point",):
            if src in out.columns:
                out["time"] = pd.to_numeric(out[src], errors="coerce")
                break

    # response ← response_value (numeric)
    if "response" not in out.columns:
        for src in ("response_value", "value", "effect"):
            if src in out.columns:
                out["response"] = pd.to_numeric(out[src], errors="coerce")
                break

    # guarantee columns exist
    for col, default in [
        ("dose", np.nan),
        ("dose_unit", pd.NA),
        ("time_unit", pd.NA),
        ("response_unit", pd.NA),
        ("replicate_id", pd.NA),
        ("normalized_response", np.nan),
        ("transformed_dose", np.nan),
        ("transformed_response", np.nan),
        ("preprocessing_flags", ""),
    ]:
        if col not in out.columns:
            out[col] = default

    # Ensure dose and response are numeric
    out["dose"] = pd.to_numeric(out["dose"], errors="coerce")
    if "response" in out.columns:
        out["response"] = pd.to_numeric(out["response"], errors="coerce")

    return out

def _step_filter(
    df: pd.DataFrame,
    *,
    source_filters: dict[str, list[str | float | bool | None]] | None = None,
    ke_ids: list[str] | None = None,
    time_points: list[float] | None = None,
    doses: list[float] | None = None,
    dose_min: float | None = None,
    dose_max: float | None = None,
    keep_control: bool = True,
    control_dose: float = 0.0,
) -> pd.DataFrame:
    selectors = (source_filters, ke_ids, time_points, doses, dose_min, dose_max)
    if any(value is not None for value in selectors) or not keep_control or control_dose != 0:
        raise ValueError("The demo retains all source rows and the 0 nM control.")
    return df.copy().reset_index(drop=True)


def _step_harmonize(df: pd.DataFrame, dose_strategy: str, time_strategy: str) -> pd.DataFrame:
    if dose_strategy != "keep_original" or time_strategy != "keep_original":
        raise ValueError("The demo retains the original dose and time units.")
    out = df.copy()
    out["dose_original"] = pd.to_numeric(out["dose"], errors="coerce")
    out["dose_unit_original"] = out["dose_unit"].astype(str)
    out["time_original"] = pd.to_numeric(out.get("time", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["time_unit_original"] = out["time_unit"].astype(str) if "time_unit" in out.columns else ""
    out["dose_harmonized"] = np.nan
    out["dose_unit_harmonized"] = ""
    out["time_harmonized"] = np.nan
    out["time_unit_harmonized"] = ""
    out["unit_harmonization_status"] = "original_units_retained"
    out["unit_harmonization_warnings"] = ""
    out["dose_conversion_factor_to_M"] = np.nan
    out["time_conversion_factor_to_h"] = np.nan

    out["dose_harmonized"] = out["dose_original"]
    out["dose_unit_harmonized"] = out["dose_unit_original"]
    out["time_harmonized"] = out["time_original"]
    out["time_unit_harmonized"] = out["time_unit_original"]

    out["harmonized_dose_unit"] = out["dose_unit_harmonized"]
    out["harmonized_time_unit"] = out["time_unit_harmonized"]
    out["harmonized_response_unit"] = out.get("response_unit", "")
    return out

def _is_control(df: pd.DataFrame, control_dose: float = 0.0) -> pd.Series:
    ctrl = pd.Series(False, index=df.index)
    if "is_control" in df.columns:
        ctrl = ctrl | df["is_control"].isin([True, "TRUE", "true", "yes", "control"])
    for col in ("group_id", "treatment", "treatment_id"):
        if col in df.columns:
            ctrl = ctrl | df[col].astype(str).str.lower().str.contains("control", na=False)
    ctrl = ctrl | (pd.to_numeric(df["dose"], errors="coerce") == control_dose)
    return ctrl

def _log(process: str, action: str, before: int, after: int, details: str = "") -> ActionLogEntry:
    return ActionLogEntry(
        process=process,
        action=action,
        rows_before=before,
        rows_after=after,
        rows_affected=abs(after - before),
        details=details,
    )

def _long_to_wide_reps(
    df: pd.DataFrame, *, source_map: dict | None = None,
) -> tuple[pd.DataFrame, list[str], str, str]:
    """Pivot long KE data to wide-with-replicates for ReplicateHandler.

    Returns (wide_df, ke_cols, dose_col, rep_col).
    wide_df columns: dose_col | rep_col | ke1 | ke2 | ...
    The caller supplies one sampling-time cell at a time.
    """
    ke_col = next((c for c in ("ke_id", "ke_display_label", "endpoint_id") if c in df.columns), None)
    response_col = next((c for c in ("response_value", "response", "value") if c in df.columns), None)
    rep_col = next((c for c in ("replicate_id", "replicate", "rep") if c in df.columns), None)
    dose_col = "dose"

    if not ke_col or not response_col or dose_col not in df.columns:
        return df.copy(), [], dose_col, rep_col or "replicate_id"
    work = df.copy()
    if rep_col is None:
        rep_col = "replicate_id"
        work[rep_col] = pd.NA
    original_reps = work["replicate_id"].copy() if "replicate_id" in work else work[rep_col].copy()
    # Missing IDs are independent observations. A shared placeholder would
    # collapse distinct measurements before sample counts and SD are computed.
    missing_reps = work[rep_col].isna() | work[rep_col].astype(str).str.strip().eq("")
    if missing_reps.any():
        work.loc[missing_reps, rep_col] = [
            f"__unlabeled_observation_{position}"
            for position in np.flatnonzero(missing_reps.to_numpy())
        ]

    work[response_col] = pd.to_numeric(work[response_col], errors="coerce")
    # Deduplicate repeated measurements of one replicate within this time cell.
    group_cols = [dose_col, rep_col, ke_col]
    agg = work.groupby(group_cols, dropna=False)[response_col].mean().reset_index()

    wide = agg.pivot_table(
        index=[dose_col, rep_col], columns=ke_col, values=response_col, aggfunc="mean",
    ).reset_index()
    wide.columns.name = None
    ke_cols_out = [c for c in wide.columns if c not in (dose_col, rep_col)]
    if source_map is not None:
        for index, row in wide.iterrows():
            for ke in ke_cols_out:
                mask = (work[dose_col] == row[dose_col]) & (work[rep_col] == row[rep_col]) & (work[ke_col] == ke) & work[response_col].notna()
                observations = work.loc[mask, "_bn_observation_id"].tolist()
                reps = original_reps.loc[mask].dropna().astype(str)
                source_map[(index, ke)] = (observations, sorted(set(reps[reps.str.strip().ne("")].tolist())))
    return wide, ke_cols_out, dose_col, rep_col

def _wide_agg_to_long(
    agg_df: pd.DataFrame,
    dose_col: str,
    ke_cols: list[str],
    meta_df: pd.DataFrame,
    ke_src_col: str = "ke_id",
    *,
    time_value: float | None = None,
    independent_wide: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Melt Format A (wide aggregated) back to long format for pipeline Steps 4–6.

    Input : dose | KE1_mean | KE1_sd | KE1_n | KE2_mean | ...
    Output: ke_src_col | dose | response_value | response |
            replicate_sd | replicate_n | time | time_point | ...
    When independent_wide is provided, paired/bootstrap rows never determine n or SD.
    """
    rows = []
    for _, row in agg_df.iterrows():
        for ke in ke_cols:
            if independent_wide is not None:
                observed = pd.to_numeric(
                    independent_wide.loc[independent_wide[dose_col] == row[dose_col], ke],
                    errors="coerce",
                ).dropna()
                n = len(observed)
                response = observed.mean()
                spread = observed.std(ddof=1) if n > 1 else 0.0
            else:
                n_raw = row.get(f"{ke}_n", 1)
                n = int(n_raw) if np.isfinite(float(n_raw if n_raw is not None else 1)) else 1
                response = row.get(f"{ke}_mean", np.nan)
                spread = row.get(f"{ke}_sd", np.nan)
            rows.append({
                ke_src_col: ke,
                dose_col: row[dose_col],
                "response_value": response,
                "response": response,
                "replicate_sd": spread,
                "replicate_n": n,
                "time": time_value,
                "time_point": time_value,
            })

    long_df = pd.DataFrame(rows)

    # Carry over per-KE metadata from original long-format DataFrame
    meta_candidates = [c for c in (
        "ke_display_label", "ke_name", "node_type", "ke_type",
        "response_unit", "dose_unit", "time_unit", "endpoint_id", "dataset_id",
    ) if c in meta_df.columns]
    if meta_candidates and ke_src_col in meta_df.columns:
        meta_keys = [ke_src_col, dose_col]
        meta_cols = [c for c in meta_candidates if c not in meta_keys]
        meta = meta_df.groupby(meta_keys, dropna=False)[meta_cols].first().reset_index()
        long_df = long_df.merge(meta, on=meta_keys, how="left")

    if "endpoint_id" not in long_df.columns:
        long_df["endpoint_id"] = long_df[ke_src_col].astype(str)

    return long_df

def preprocess_dataframe(
    df: pd.DataFrame,
    *,
    # Step 1: selection
    source_filters: dict[str, list[str | float | bool | None]] | None = None,
    ke_ids: list[str] | None = None,
    time_points: list[float] | None = None,
    doses: list[float] | None = None,
    dose_min: float | None = None,
    dose_max: float | None = None,
    keep_control: bool = True,
    control_dose: float = 0.0,
    # Step 2: missing
    missing: MissingStrategy = "retain_na",
    # Step 3: replicates (legacy 3-mode OR new 4-mode via replicate_config)
    replicates: ReplicateStrategy = "retain_replicates",
    replicate_config: Any | None = None,
    # Step 4: harmonization
    harmonization: HarmonizeStrategy = "keep_original",
    # Step 5: normalization
    normalization: NormStrategy = "none",
    # Step 6: transformation
    log_dose: bool = False,
    log_response: bool = False,
    reverse_direction: bool = False,
    bn_artifacts: dict[str, Any] | None = None,
    strict_source_units: bool = False,
    preserve_source_ids: bool = False,
) -> tuple[pd.DataFrame, PreprocessSummary, list[ActionLogEntry]]:
    if (
        missing != "retain_na"
        or replicates != "retain_replicates"
        or harmonization != "keep_original"
        or normalization != "none"
        or log_dose
        or log_response
        or reverse_direction
    ):
        raise ValueError("The demo uses the recorded preparation options.")
    out = _standardize_columns(df)
    has_sampling_time = "time" in out.columns
    capture = {} if (bn_artifacts is not None) and replicate_config is not None else None
    if capture is not None and (not preserve_source_ids or "_bn_observation_id" not in out):
        out["_bn_observation_id"] = [f"observation_{i}" for i in range(len(out))]
    actions: list[ActionLogEntry] = []

    before = len(out)
    out = _step_filter(
        out, source_filters=source_filters, ke_ids=ke_ids, time_points=time_points, doses=doses,
        dose_min=dose_min, dose_max=dose_max,
        keep_control=keep_control, control_dose=control_dose,
    )
    selection_details = f"source_filters={source_filters};ke_ids={ke_ids};time_points={time_points};dose_min={dose_min};dose_max={dose_max}"
    if doses is not None:
        selection_details += f";doses={doses}"
    actions.append(_log("data_selection_filtering", "selection_applied", before, len(out), selection_details))

    combines_rows = replicate_config is not None
    if combines_rows and "dataset_id" in out.columns:
        datasets = out["dataset_id"].dropna().astype(str).str.strip()
        if datasets[datasets.ne("")].nunique() > 1:
            raise ValueError(
                "Step 2 cannot combine multiple dataset_id values. Select or upload one "
                "compatible analysis series before replicate handling or normalization."
            )
    if combines_rows:
        for column in ("stressor", "chemical_name", "stressor_name"):
            if column in out.columns:
                identities = out[column].fillna("<missing>").astype(str).str.strip().replace("", "<missing>")
                if identities.nunique() > 1:
                    raise ValueError(
                        f"Step 2 cannot combine multiple {column} values. Select one compatible "
                        "stressor in shared data selection before replicate handling or normalization."
                    )
        for column in ("experiment_id", "study_id"):
            if column in out.columns:
                group_keys = ["endpoint_id", "dose"] + (["time"] if "time" in out.columns else [])
                for _, group in out.groupby(group_keys, dropna=False):
                    identities = group[column].fillna("<missing>").astype(str).str.strip().replace("", "<missing>")
                    if identities.nunique() > 1:
                        raise ValueError(
                            f"Step 2 cannot combine multiple {column} values within one KE/dose/time cell. "
                            "Select a compatible experiment in shared data selection."
                        )
    if strict_source_units:
        for column, factor_for_unit in (
            ("dose_unit", _dose_factor),
            ("time_unit", _time_factor),
        ):
            if column == "time_unit" and not has_sampling_time:
                continue
            identities = set()
            for value in out[column].unique():
                if pd.isna(value) or not str(value).strip():
                    identities.add(("missing", ""))
                else:
                    factor = factor_for_unit(str(value))
                    identities.add(("factor", factor) if factor is not None
                                   else ("literal", _normalize_unit(str(value))))
            if len(identities) > 1:
                raise ValueError(
                    f"Step 2 has mixed {column} values that cannot be combined before "
                    "replicate handling or normalization. Use one source unit per run, or "
                    "retain individual rows with compatible unit harmonization."
                )
    for endpoint, rows in (out.groupby("endpoint_id", dropna=False) if "endpoint_id" in out else []):
        response_units = {
            "<missing>" if pd.isna(value) or not str(value).strip() else str(value).strip().casefold()
            for value in rows["response_unit"].unique()
        }
        if len(response_units) > 1:
            raise ValueError(
                f"Step 2 response_unit differs within KE '{endpoint}'. "
                "Provide a consistent response scale for each KE before preprocessing."
            )

    before = len(out)
    _has_response = "response" in out.columns
    missing_before = int(out["response"].isna().sum()) if _has_response else 0
    missing_after = missing_before
    actions.append(_log("missing_data_handling", missing, before, len(out),
                        f"missing_before={missing_before};missing_after={missing_after}"))

    before = len(out)
    no_prepared_cases = False
    if replicate_config is not None:
        # Four-mode path: wide-format pivot → ReplicateHandler → back to long
        from app.upstream.replicate_handling import ReplicateHandler
        groups = out.groupby("time", dropna=False) if "time" in out.columns else [(None, out)]
        result_parts = []
        quality_flags = []
        inflation_warning = None
        pairing_attempted = False
        for time_value, time_rows in groups:
            metadata_rows = time_rows
            source_map = {} if capture is not None else None
            wide, ke_cols_rh, dose_col_rh, rep_col_rh = _long_to_wide_reps(
                time_rows, source_map=source_map,
            )
            if not ke_cols_rh:
                continue
            pairing_attempted = True
            ke_src = next(
                (c for c in ("ke_id", "ke_display_label", "endpoint_id") if c in time_rows.columns),
                "ke_id",
            )
            rh_result = ReplicateHandler().handle(
                wide, ke_cols_rh, dose_col_rh, rep_col_rh, replicate_config
            )
            if capture is not None:
                from app.upstream.bn_cases import capture_pairing
                capture_pairing(capture, rh_result, wide, source_map, metadata_rows, time_value,
                                ke_cols_rh, dose_col_rh, ke_src)
            if not rh_result.format_a.df.empty:
                result_parts.append(_wide_agg_to_long(
                    rh_result.format_a.df, dose_col_rh, ke_cols_rh, metadata_rows,
                    ke_src_col=ke_src, time_value=time_value, independent_wide=wide,
                ))
            quality_flags.extend(rh_result.quality_flags)
            inflation_warning = rh_result.inflation_warning or inflation_warning
        if result_parts:
            out = _standardize_columns(pd.concat(result_parts, ignore_index=True))
            out.attrs["_replicate_quality_flags"] = quality_flags
            out.attrs["_replicate_inflation_warning"] = inflation_warning
            actions.append(_log(
                "replicate_handling", replicate_config.mode, before, len(out),
                f"mode={replicate_config.mode};n_flags={len(quality_flags)}",
            ))
        elif pairing_attempted:
            # A selected pairing with no complete dose cells is not raw-row input.
            no_prepared_cases = True
            out = _standardize_columns(pd.DataFrame(columns=[
                ke_src, "endpoint_id", "dose", "time", "time_point", "response", "response_value",
                "replicate_sd", "replicate_n", "dose_unit", "time_unit", "response_unit",
            ]))
            out.attrs["_replicate_quality_flags"] = quality_flags
            out.attrs["_replicate_inflation_warning"] = inflation_warning
            actions.append(_log("replicate_handling", replicate_config.mode, before, 0, "no_prepared_cases"))
        else:
            out.attrs["_replicate_quality_flags"] = []
            out.attrs["_replicate_inflation_warning"] = None
            actions.append(_log("replicate_handling", replicate_config.mode, before, len(out),
                                "no_ke_cols_detected"))
    else:
        actions.append(_log("replicate_handling", replicates, before, len(out)))

    before = len(out)
    out = _step_harmonize(out, "keep_original", "keep_original")
    actions.append(_log("dose_time_unit_harmonization", "original_units_retained", before, len(out)))

    out["normalization_method"] = "none"
    actions.append(_log("scaling_normalization", "none", before, len(out),
                        f"normalized_rows={int(out['normalized_response'].notna().sum())}"))
    actions.append(_log("transformation", "log_dose=False;log_response=False;reverse=False",
                        before, len(out)))

    def require_finite(values: pd.Series, operation: str) -> pd.Series:
        numeric = pd.to_numeric(values, errors="coerce")
        invalid = ~np.isfinite(numeric)
        if invalid.any():
            raise ValueError(f"Step 2 {operation} is undefined for {int(invalid.sum())} row(s).")
        return numeric

    if out.empty and not no_prepared_cases:
        raise ValueError("Step 2 selected no rows for model input.")
    out["fit_dose"] = require_finite(out["dose"], "dose")
    if not has_sampling_time:
        out = out.drop(columns=["time", "time_point"], errors="ignore")
    if "time" in out.columns:
        out["fit_time"] = require_finite(out["time"], "time")

    response_source = pd.to_numeric(out["response"], errors="coerce")
    if np.isinf(response_source).any():
        raise ValueError("Step 2 response contains an infinite effective value.")
    out["fit_response"] = response_source

    numeric_cols = out.select_dtypes(include="number").columns.tolist()
    # Significant digits retain the magnitude of small measured responses;
    # fixed decimal places can misleadingly display a nonzero value as zero.
    column_stats: dict[str, Any] = {
        col: {k: (None if not np.isfinite(v) else int(v) if k == "count"
                  else float(f"{float(v):.8g}"))
              for k, v in out[col].describe().items()}
        for col in numeric_cols
        if col in out.columns
    }
    summary = PreprocessSummary(
        n_rows_before=len(df),
        n_rows_after=len(out),
        n_cols=len(out.columns),
        missing_handled=int(out["preprocessing_flags"].str.contains("missing", na=False).sum()),
        missing_response_rows=int(out["fit_response"].isna().sum()),
        outliers_removed=0,
        normalization=normalization,
        column_stats=column_stats,
    )
    if capture is not None:
        from app.upstream.bn_cases import finalize_artifacts
        prepared_artifacts = finalize_artifacts(capture, {
            "harmonization": harmonization, "normalization": normalization,
            "control_dose": control_dose, "log_dose": log_dose,
            "log_response": log_response, "reverse_direction": reverse_direction,
            "replicate_config": replicate_config.model_dump(),
        })
        if bn_artifacts is not None:
            bn_artifacts.update(prepared_artifacts)
    return out, summary, actions
