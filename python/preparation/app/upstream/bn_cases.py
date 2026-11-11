# Extracted from qAOPweaver/app/upstream/bn_cases.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

import hashlib

import json

from typing import Any

import numpy as np

import pandas as pd

from app.upstream.preprocessed_input import MissingPreprocessedInput

def _source_conditions(metadata, dose, time):
    return list({(float(dose), None if pd.isna(time) else float(time),
                  str(row.get("dose_unit", "")), str(row.get("time_unit", "")))
                 for _, row in metadata.iterrows()})

def capture_pairing(capture, result, wide, source_map, time_rows, time_value,
                    nodes, dose_col, ke_col):
    """Capture values and indices from the selected pairing, without running it again."""
    case_rows = capture.setdefault("case_rows", [])
    lineage = capture.setdefault("lineage", [])
    counts = capture.setdefault("counts", [])
    capture.setdefault("nodes", set()).update(str(node) for node in nodes)
    capture.setdefault("warnings", []).extend([result.inflation_warning] if result.inflation_warning else [])
    capture["mode"] = result.mode
    for dose, group in wide.groupby(dose_col):
        for node in nodes:
            for source_index in group.index[group[node].notna()]:
                observations, _ = source_map[(source_index, node)]
                metadata = time_rows.loc[time_rows["_bn_observation_id"].isin(observations)]
                first = metadata.iloc[0]
                counts.append({"node_id": str(node), "time": time_value, "dose": dose,
                               "dose_unit": first.get("dose_unit", ""),
                               "time_unit": first.get("time_unit", ""),
                               "response_unit": first.get("response_unit", ""), "n": 1,
                               "_source_conditions": _source_conditions(metadata, dose, time_value)})

    if result.mode == "aggregate_replicates":
        values = result.format_a.df.rename(columns={f"{node}_mean": node for node in nodes})
        runs = [values]
    elif result.mode in ("random_pairing_bootstrap", "repeated_one_to_one", "with_replacement_network") and result.format_c.runs:
        runs = result.format_c.runs
    else:
        runs = [result.format_b.df]
    capture.setdefault("available_runs", set()).update(str(run_id) for run_id in range(len(runs)))
    for run_id, (values, sources) in enumerate(zip(runs, result.source_runs)):
        for position, row in values.iterrows():
            case_id = f"case_{len(case_rows)}"
            for node in nodes:
                dose = row[dose_col]
                observations, reps = [], []
                for source_index in sources.iloc[position][node]:
                    source_observations, source_reps = source_map[(source_index, node)]
                    observations.extend(source_observations)
                    reps.extend(source_reps)
                observations = sorted(set(observations))
                metadata = time_rows.loc[time_rows["_bn_observation_id"].isin(observations)]
                source_conditions = _source_conditions(metadata, dose, time_value)
                if metadata.empty:
                    metadata = time_rows.loc[(time_rows[ke_col] == node) & (time_rows[dose_col] == dose)]
                first = metadata.iloc[0] if not metadata.empty else {}
                case_rows.append({
                    "case_id": case_id, "run_id": str(run_id), "node_id": str(node),
                    "local_edge_id": row.get("local_edge_id") if result.mode == "ker_local_all_pairs" else None,
                    "dose": dose, "time": time_value, "response": row[node],
                    "dose_unit": first.get("dose_unit", ""),
                    "time_unit": first.get("time_unit", ""),
                    "response_unit": first.get("response_unit", ""),
                    "_source_conditions": source_conditions,
                })
                source_key = (
                    f"condition:{case_id}:{node}" if result.mode == "aggregate_replicates"
                    else "source:" + hashlib.sha256(json.dumps(observations).encode()).hexdigest()
                ) if observations else ""
                lineage.append({
                    "case_id": case_id, "run_id": str(run_id), "node_id": str(node),
                    "source_key": source_key,
                    "source_observation_ids": json.dumps(observations),
                    "source_replicate_ids": json.dumps(sorted(set(reps))),
                })

def _effective_conditions(frame, operations):
    from app.upstream.preprocess import _step_harmonize
    if operations["harmonization"] != "keep_original" or operations["log_dose"]:
        raise ValueError("The demo retains the original dose and time scales.")
    frame = _step_harmonize(frame, "keep_original", "keep_original")
    frame["fit_dose"] = frame.dose_harmonized
    frame["fit_time"] = frame.time_harmonized
    return frame

def _verify_source_conditions(capture, operations):
    """A value already averaged across incompatible conditions cannot be relocated."""
    from app.upstream.preprocess import _dose_factor, _normalize_unit, _time_factor

    def unit_identity(unit, kind):
        factor = (_dose_factor if kind == "dose" else _time_factor)(unit)
        return (kind, factor) if factor is not None else (kind, _normalize_unit(str(unit)))

    records = capture["case_rows"] + capture["counts"]
    conditions = sorted({condition for record in records for condition in record["_source_conditions"]},
                        key=str)
    original = pd.DataFrame(conditions, columns=["dose", "time", "dose_unit", "time_unit"])
    effective = _effective_conditions(original, operations)
    mapped = {
        condition: (None if pd.isna(row.fit_time) else float(row.fit_time),
                    None if pd.isna(row.fit_dose) else float(row.fit_dose),
                    unit_identity(row.dose_unit_harmonized, "dose"),
                    unit_identity(row.time_unit_harmonized, "time"))
        for condition, row in zip(conditions, effective.itertuples())
    }
    for record in records:
        if len({mapped[condition] for condition in record["_source_conditions"]}) > 1:
            raise ValueError(f"Step 2 prepared value for node '{record['node_id']}' spans multiple effective conditions; "
                             "make source conditions compatible before preparing BN cases.")

def _case_tables_by_effective_condition(frame, capture, nodes):
    """Split an existing pairing when its nodes resolve to different conditions."""
    source_lineage = {(row["case_id"], row["node_id"]): row for row in capture["lineage"]}
    cases, lineage = [], []
    split = False
    for original_id, prepared in frame.groupby("case_id", sort=False):
        conditions = list(prepared.groupby(["fit_time", "fit_dose"], sort=False, dropna=False))
        split = split or len(conditions) > 1
        for number, ((time, dose), values) in enumerate(conditions):
            case_id = f"{original_id}_condition_{number}" if len(conditions) > 1 else original_id
            run_id = values.run_id.iloc[0]
            case = {"case_id": case_id, "run_id": run_id, "fit_time": time, "fit_dose": dose,
                    "local_edge_id": values.local_edge_id.iloc[0] if "local_edge_id" in values else None,
                    **{node: np.nan for node in nodes}}
            condition_nodes = set(values.node_id)
            case.update(dict(zip(values.node_id, values.fit_response)))
            cases.append(case)
            for node in nodes:
                original = source_lineage.get((original_id, node)) if node in condition_nodes else None
                if original is not None:
                    record = {**original, "case_id": case_id}
                    if capture["mode"] == "aggregate_replicates" and record["source_key"]:
                        record["source_key"] = f"condition:{case_id}:{node}"
                else:
                    record = {"case_id": case_id, "run_id": run_id, "node_id": node,
                              "source_key": "", "source_observation_ids": "[]", "source_replicate_ids": "[]"}
                lineage.append(record)
    case_columns = (["case_id", "run_id", "fit_time", "fit_dose"] +
                    (["local_edge_id"] if capture["mode"] == "ker_local_all_pairs" else []) + nodes)
    return (pd.DataFrame(cases, columns=case_columns),
            pd.DataFrame(lineage, columns=["case_id", "run_id", "node_id", "source_key",
                                           "source_observation_ids", "source_replicate_ids"]), split)

def _independent_condition_coverage(capture, operations, cases, nodes):
    """Keep original conditions even when the selected pairing produced no case."""
    original = _effective_conditions(pd.DataFrame(capture["counts"]), operations)
    counts, coverage = [], []
    for (time, dose), condition in original.groupby(["fit_time", "fit_dose"], sort=False, dropna=False):
        n_by_node = condition.groupby("node_id")["n"].sum().to_dict()
        effective_time = None if pd.isna(time) else float(time)
        effective_dose = None if pd.isna(dose) else float(dose)
        matching = (cases.fit_time.isna() if pd.isna(time) else cases.fit_time.eq(time)) & cases.fit_dose.eq(dose)
        coverage.append({"fit_time": effective_time, "fit_dose": effective_dose,
                         "observed_node_ids": [node for node in nodes if n_by_node.get(node, 0) > 0],
                         "missing_node_ids": [node for node in nodes if n_by_node.get(node, 0) == 0],
                         "prepared_case_count": int(matching.sum())})
        for node in nodes:
            counts.append({"node_id": node, "fit_time": effective_time,
                           "fit_dose": effective_dose, "n": int(n_by_node.get(node, 0))})
    return counts, coverage

def finalize_artifacts(capture: dict, operations: dict) -> dict:
    """Record paired values on the original response scale."""
    if operations["normalization"] != "none" or operations["log_response"] or operations["reverse_direction"]:
        raise ValueError("The demo retains the original response scale.")
    if not capture.get("case_rows") and not capture.get("counts"):
        return {}
    _verify_source_conditions(capture, operations)
    original_cases = pd.DataFrame(capture["case_rows"]) if capture["case_rows"] else pd.DataFrame(columns=[
        "case_id", "run_id", "node_id", "dose", "time", "response", "dose_unit", "time_unit", "response_unit",
    ])
    frame = _effective_conditions(original_cases, operations)
    frame["fit_response"] = pd.to_numeric(frame.response, errors="coerce")
    nodes = sorted(capture["nodes"])
    cases, lineage, split = _case_tables_by_effective_condition(frame, capture, nodes)
    counts, coverage = _independent_condition_coverage(capture, operations, cases, nodes)
    units = {}
    unit_frame = frame if not frame.empty else _effective_conditions(pd.DataFrame(capture["counts"]), operations)
    for kind, col in (("dose", "dose_unit_harmonized"), ("time", "time_unit_harmonized"), ("response", "response_unit")):
        units[kind] = sorted(set(str(value) for value in unit_frame[col].dropna()))
    warnings = list(dict.fromkeys(capture["warnings"]))
    if split:
        warnings.append("Prepared cases with different per-node effective conditions were split into partial cases; no cross-condition pairing was performed.")
    if any(condition["prepared_case_count"] == 0 for condition in coverage):
        warnings.append("Some original effective conditions have zero prepared cases; independent counts and missing-node coverage are retained.")
    if not np.isfinite(cases[nodes].to_numpy(dtype=float)).all():
        warnings.append("Some prepared cases have missing or undefined node values; evaluate complete-case coverage per time.")
    metadata = {
        "schema_version": "bn_step2_cases_v1", "node_ids": nodes,
        "pairing_mode": capture["mode"], "available_runs": sorted(capture["available_runs"], key=int),
        "effective_units": units, "effective_operations": operations,
        "original_independent_counts": counts, "dependency_warnings": warnings,
        "condition_coverage": coverage,
    }
    return {"bn_cases": cases, "bn_lineage": lineage, "bn_case_metadata": metadata}
