# Extracted from qAOPweaver/app/upstream/prepared_view.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

import hashlib

import json

from typing import Any

import numpy as np

import pandas as pd

from app.upstream.preprocessed_input import MissingPreprocessedInput

SCHEMA_VERSION = "step2_prepared_view_v1"

_CASE_COLUMNS = [
    "case_id", "condition_id", "run_id", "realization_id", "fit_dose",
    "fit_time", "time_unit", "true_unit_id", "local_edge_id", "case_weight", "case_complete",
]

_SOURCE_TUPLE_COLUMNS = [
    "case_id", "condition_id", "run_id", "node_id", "source_key",
    "source_observation_ids", "source_replicate_ids",
]

_SOURCE_LEDGER_COLUMNS = [
    "source_observation_id", "source_row_index", "uploaded_observation_id",
    "true_unit_id", "node_id", "dataset_id", "stressor", "experiment_id",
    "dose", "time", "response", "dose_unit", "time_unit", "response_unit",
    "fit_dose", "fit_time", "condition_id", "selected_for_preparation",
    "used_in_cases", "exclusion_reason",
]

_ROLE_NAMES = {"MIE", "KE", "AO", "NAME NOT PROVIDED"}

def _digest_parts(*parts: Any) -> str:
    digest = hashlib.sha256()
    for part in parts:
        if isinstance(part, pd.DataFrame):
            payload = part.to_json(orient="records", double_precision=15).encode()
        else:
            payload = json.dumps(part, sort_keys=True, allow_nan=False, separators=(",", ":")).encode()
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()

def _clean(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text if text and text.casefold() not in ("nan", "none", "null", "<na>") else None

def _first_nonempty(row: pd.Series, columns: tuple[str, ...]) -> str | None:
    for column in columns:
        if column in row:
            value = _clean(row[column])
            if value:
                return value
    return None

def _source_ledger(raw: pd.DataFrame, options: dict[str, Any]) -> pd.DataFrame:
    from app.upstream.preprocess import _standardize_columns, _step_filter, _step_harmonize

    source = _standardize_columns(raw).reset_index(drop=True)
    source["_prepared_source_id"] = [f"observation_{index}" for index in range(len(source))]
    selected_common = _step_filter(
        source,
        source_filters=options.get("source_filters"),
        ke_ids=options.get("ke_ids"), time_points=options.get("time_points"),
        doses=options.get("doses"), dose_min=options.get("dose_min"),
        dose_max=options.get("dose_max"),
        keep_control=options.get("keep_control", True),
        control_dose=options.get("control_dose", 0.0),
    )
    common_ids = set(selected_common["_prepared_source_id"])
    selected = selected_common
    selected_ids = common_ids
    effective = _step_harmonize(source, "keep_original", "keep_original")
    dose_values = pd.to_numeric(effective["dose_harmonized"], errors="coerce")
    time_values = pd.to_numeric(effective["time_harmonized"], errors="coerce")
    rows = []
    for index, row in source.iterrows():
        source_id = row["_prepared_source_id"]
        kept = source_id in selected_ids
        rows.append({
            "source_observation_id": source_id,
            "source_row_index": index,
            "uploaded_observation_id": _first_nonempty(row, ("observation_id", "obs_id")),
            # A replicate label is not proof of a repeated experimental unit.
            "true_unit_id": None,
            "node_id": _first_nonempty(row, ("endpoint_id", "ke_id")),
            "dataset_id": _first_nonempty(row, ("dataset_id",)),
            "stressor": _first_nonempty(row, ("stressor", "chemical_name", "stressor_name")),
            "experiment_id": _first_nonempty(row, ("experiment_id", "study_id")),
            "dose": row.get("dose"), "time": row.get("time"),
            "response": row.get("response"),
            "dose_unit": row.get("dose_unit"), "time_unit": row.get("time_unit"),
            "response_unit": row.get("response_unit"),
            "fit_dose": dose_values.iloc[index], "fit_time": time_values.iloc[index],
            "condition_id": None,
            "selected_for_preparation": kept,
            "used_in_cases": False,
            "exclusion_reason": "" if kept else "source_filter",
        })
    return pd.DataFrame(rows, columns=_SOURCE_LEDGER_COLUMNS)

def _resolved_labels(raw: pd.DataFrame, path: pd.DataFrame | None, node_ids: list[str]) -> dict[str, str]:
    from app.upstream.preprocess import _standardize_columns

    source = _standardize_columns(raw)
    found: dict[str, set[str]] = {node: set() for node in node_ids}

    def add(node: Any, value: Any) -> None:
        identifier, label = _clean(node), _clean(value)
        if not identifier or identifier not in found or not label:
            return
        if label == identifier or label.upper() in _ROLE_NAMES:
            return
        found[identifier].add(label)

    if "endpoint_id" in source:
        for row in source.to_dict("records"):
            for column in ("ke_display_label", "ke_name"):
                add(row["endpoint_id"], row.get(column))
    if path is not None and not path.empty:
        for row in path.to_dict("records"):
            for prefix in ("from", "to"):
                node = row.get(f"{prefix}_ke_id", row.get(f"{prefix}_node_id"))
                for column in (f"{prefix}_display_label", f"{prefix}_ke_name", f"{prefix}_ke_label"):
                    add(node, row.get(column))
    conflicts = {node: sorted(labels) for node, labels in found.items() if len(labels) > 1}
    if conflicts:
        detail = "; ".join(f"{node}: {', '.join(labels)}" for node, labels in conflicts.items())
        raise MissingPreprocessedInput(
            f"Step 2 has conflicting KE names for the same stable ID ({detail}). "
            "Resolve the source or path names before fitting."
        )
    return {node: next(iter(found[node]), "Name not provided") for node in node_ids}

def _scope(ledger: pd.DataFrame) -> dict[str, list[str]]:
    selected = ledger.loc[ledger["selected_for_preparation"]]
    return {
        column: sorted({str(value) for value in selected[column].dropna() if str(value).strip()})
        for column in ("dataset_id", "stressor", "experiment_id")
    }

def _mixed_source_fields(ledger: pd.DataFrame) -> list[str]:
    """Detect source identities that cannot share one prepared condition."""
    selected = ledger.loc[ledger["selected_for_preparation"]]
    return [column for column in ("dataset_id", "stressor", "experiment_id")
            if len({_clean(value) for value in selected[column]}) > 1]

def _path_in_scope(path: pd.DataFrame | None, scope: dict[str, list[str]]) -> pd.DataFrame | None:
    if path is None:
        return None
    selected = path.copy()
    for column, values in scope.items():
        if column in selected and values:
            selected = selected.loc[selected[column].astype(str).isin(values)]
    return selected

def _selected_path_edge_count(path: pd.DataFrame | None, node_ids: list[str]) -> int:
    if path is None or path.empty:
        return 0
    selected_nodes = set(node_ids)
    edges = set()
    for row in path.to_dict("records"):
        upstream = _clean(row.get("from_ke_id", row.get("from_node_id", row.get("from"))))
        downstream = _clean(row.get("to_ke_id", row.get("to_node_id", row.get("to"))))
        if upstream in selected_nodes and downstream in selected_nodes and upstream != downstream:
            edges.add((upstream, downstream))
    return len(edges)

def _unit_identity(values: pd.Series) -> str | None:
    observed = {_clean(value) for value in values}
    observed.discard(None)
    return next(iter(observed)) if len(observed) == 1 else None

def _effective_units(processed: pd.DataFrame, ledger: pd.DataFrame,
                     node_ids: list[str], options: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Expose source and fit scales, leaving unknown or mixed units explicit."""
    dose_col = "dose_unit_harmonized" if "dose_unit_harmonized" in processed else "dose_unit"
    time_col = "time_unit_harmonized" if "time_unit_harmonized" in processed else "time_unit"
    dose_unit = _unit_identity(processed[dose_col]) if dose_col in processed else None
    time_unit = _unit_identity(processed[time_col]) if time_col in processed else None
    selected = ledger.loc[ledger["selected_for_preparation"]]
    source_units = {
        "dose": _unit_identity(selected["dose_unit"]),
        "time": _unit_identity(selected["time_unit"]),
    }
    node_units = {}
    for node in node_ids:
        node_source = selected.loc[selected["node_id"].eq(node)]
        raw_response_unit = _unit_identity(node_source["response_unit"])
        node_units[node] = {
            "source_response_unit": raw_response_unit,
            "effective_response_unit": raw_response_unit,
        }
    return ({
        "source_dose_unit": source_units["dose"],
        "source_time_unit": source_units["time"],
        "dose_unit": dose_unit,
        "time_unit": time_unit,
        "dose_transform": "none",
        "response_transform": "none",
    }, node_units)

def _cases_and_sources(artifacts: dict[str, Any], node_ids: list[str],
                       scope: dict[str, list[str]], pairing_mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_cases = artifacts.get("bn_cases")
    source_tuples = artifacts.get("bn_lineage")
    if not isinstance(source_cases, pd.DataFrame) or not isinstance(source_tuples, pd.DataFrame):
        return (pd.DataFrame(columns=_CASE_COLUMNS + node_ids),
                pd.DataFrame(columns=_SOURCE_TUPLE_COLUMNS))
    required = {"case_id", "run_id", "fit_dose", "fit_time", *node_ids}
    if not required.issubset(source_cases) or not {"case_id", "node_id", "source_key", "source_observation_ids"}.issubset(source_tuples):
        raise MissingPreprocessedInput("Step 2 prepared cases or source lineage are incomplete; run Step 2 again.")
    cases = source_cases.copy()
    if cases["case_id"].duplicated().any():
        raise MissingPreprocessedInput("Step 2 prepared case IDs are not unique; run Step 2 again.")
    metadata = artifacts.get("bn_case_metadata", {})
    time_units = metadata.get("effective_units", {}).get("time", []) if isinstance(metadata, dict) else []
    time_unit = time_units[0] if len(time_units) == 1 else None
    namespace = {column: scope[column] for column in ("dataset_id", "stressor", "experiment_id")}
    cases["condition_id"] = ["condition_" + _digest_parts(namespace, float(dose))[:24]
                             for dose in cases["fit_dose"]]
    cases["run_id"] = cases["run_id"].astype(str)
    cases["realization_id"] = cases["run_id"]
    cases["time_unit"] = time_unit
    cases["true_unit_id"] = None
    cases["local_edge_id"] = source_cases["local_edge_id"] if "local_edge_id" in source_cases else None
    numeric_nodes = cases[node_ids].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    cases["case_complete"] = np.isfinite(numeric_nodes).all(axis=1)
    counts = cases.groupby(["condition_id", "fit_time", "run_id", "local_edge_id"], dropna=False)["case_id"].transform("count")
    cases["case_weight"] = 1.0 / counts
    cases = cases[_CASE_COLUMNS + node_ids].copy()
    tuples = source_tuples.copy()
    tuples["condition_id"] = tuples["case_id"].map(cases.set_index("case_id")["condition_id"])
    tuples["run_id"] = tuples["run_id"].astype(str)
    for column in _SOURCE_TUPLE_COLUMNS:
        if column not in tuples:
            tuples[column] = "[]" if column.endswith("_ids") else ""
    tuples = tuples[_SOURCE_TUPLE_COLUMNS].copy()
    if tuples[["case_id", "node_id"]].duplicated().any() or tuples["condition_id"].isna().any():
        raise MissingPreprocessedInput("Step 2 source tuples do not match unique prepared cases.")
    return cases, tuples

def _fold_replay_reason(view_kind: str, pairing_mode: str, options: dict[str, Any],
                        source_compatibility: str, conditions: list[dict[str, Any]],
                        reference_conditions: set[str]) -> str | None:
    deterministic_constructed = {"user_selected_pairing", "all_pairwise_replicate_combinations",
                                  "whole_network_cartesian"}
    stochastic_constructed = {"random_pairing_bootstrap", "repeated_one_to_one",
                              "with_replacement_network"}
    if view_kind not in ("condition_mean", "verified_matched_cases") and not (
        view_kind == "constructed_cases" and pairing_mode in deterministic_constructed | stochastic_constructed
    ):
        return "This representation has no verified training-only fold replay."
    if options.get("bn_temporal") and (options.get("temporal_connection") or {}).get("mode") in (None, "same_time_only"):
        return "Legacy temporal provenance has no explicit expanded Step 2 network to replay."
    normalization = options.get("normalization", "none")
    if normalization not in ("none", "percent_of_control", "fold_change_from_control", "delta_from_control"):
        return "Normalization parameters must be learned from training sources within each fold."
    if normalization != "none" and not reference_conditions:
        return "Control-relative normalization needs saved reference-control observations."
    if options.get("missing", "retain_na") not in ("retain_na", "drop_missing_response", "flag_missing"):
        return "Missing-value imputation must be learned from training sources within each fold."
    if source_compatibility != "single_source_scope":
        return "Choose one compatible source scope before grouped condition assessment."
    if len(conditions) - (len(reference_conditions) if normalization != "none" else 0) < 2:
        return "Grouped assessment needs at least two selected treatment conditions."
    return None

def _expanded_selected_contexts(
    cases: pd.DataFrame, tuples: pd.DataFrame, base_nodes: list[str],
    labels: dict[str, str], path: pd.DataFrame, connection: dict[str, Any],
    time_unit: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Materialize each chosen time graph separately from exact saved cases."""
    mode = connection["mode"]
    contexts = connection["selected_contexts"]
    if not isinstance(contexts, list) or not contexts:
        raise MissingPreprocessedInput("Select at least one coherent Step 2 time context.")
    if connection.get("memory_nodes"):
        raise MissingPreprocessedInput("Selected KER contexts do not include optional self-memory edges.")
    path_edges = {(_clean(row.get("from_ke_id")), _clean(row.get("to_ke_id")))
                  for row in path.to_dict("records")}
    node_times: set[tuple[str, float]] = set()
    edges = []
    seen_edges = set()
    for context in contexts:
        for edge in context["edges"]:
            upstream, downstream = edge["from_node_id"], edge["to_node_id"]
            earlier, later = float(edge["from_time"]), float(edge["to_time"])
            if ((upstream, downstream) not in path_edges or
                    upstream not in base_nodes or downstream not in base_nodes or earlier > later):
                raise MissingPreprocessedInput("A selected Step 2 KER time pair is not an eligible uploaded edge.")
            node_times.update(((upstream, earlier), (downstream, later)))
            signature = upstream, earlier, downstream, later
            if signature not in seen_edges:
                seen_edges.add(signature)
                edges.append({"from": f"{upstream}@{earlier:g}",
                              "to": f"{downstream}@{later:g}", "kind": "aop",
                              "delta_time": later - earlier})
    graph_nodes = [{"id": f"{node}@{time:g}", "node_id": node, "fit_time": time,
                    "label": labels[node], "role": "KE"}
                   for node, time in sorted(node_times, key=lambda item: (item[1], item[0]))]
    expanded_ids = [node["id"] for node in graph_nodes]
    source_by_case_node = {(row.case_id, row.node_id): row for row in tuples.itertuples(index=False)}
    run_count = connection.get("n_realizations", 1)
    seed = connection.get("seed", 42)
    if mode != "hypothetical_cross_time" or run_count != 1 or seed != 42:
        raise MissingPreprocessedInput("The demo uses one hypothetical cross-time realization with seed 42.")
    rows, source_rows = [], []
    for context in contexts:
        required: dict[float, set[str]] = {}
        for node, time in context["node_times"].items():
            required.setdefault(float(time), set()).add(node)
        allowed_doses = set(context["matched_doses"])
        scoped = cases.loc[cases["fit_dose"].isin(allowed_doses)]
        for (condition_id, run_id), group in scoped.groupby(["condition_id", "run_id"], sort=True):
            by_time = {}
            for time, nodes in required.items():
                at_time = group.loc[group["fit_time"].eq(time)].copy()
                if at_time.empty:
                    break
                finite = np.isfinite(at_time[list(nodes)].apply(pd.to_numeric, errors="coerce")).all(axis=1)
                at_time = at_time.loc[finite].sort_values("case_id")
                if at_time.empty:
                    break
                by_time[time] = at_time
            if len(by_time) != len(required):
                continue
            count = min(map(len, by_time.values()))
            rng = np.random.default_rng(int(_digest_parts(seed, context["context_id"], condition_id, run_id)[:16], 16))
            runs = []
            for _ in range(run_count):
                sampled = {time: frame.iloc[rng.permutation(len(frame))[:count]]
                           for time, frame in by_time.items()}
                runs.append([{time: frame.iloc[index] for time, frame in sampled.items()}
                             for index in range(count)])
            for realization, sequences in enumerate(runs):
                for sequence in sequences:
                    case_id = f"temporal_case_{len(rows)}"
                    first = next(iter(sequence.values()))
                    row = {"case_id": case_id, "condition_id": condition_id,
                           "run_id": str(run_id), "realization_id": f"{run_id}:{realization}",
                           "time_context_id": context["context_id"],
                           "fit_dose": float(first["fit_dose"]), "fit_time": np.nan,
                           "time_unit": time_unit,
                           "true_unit_id": None,
                           "case_weight": 1.0 / len(sequences),
                           "case_complete": False}
                    for time, saved in sequence.items():
                        for node in required[time]:
                            expanded = f"{node}@{time:g}"
                            row[expanded] = saved[node]
                            source = source_by_case_node.get((saved["case_id"], node))
                            if source is None:
                                raise MissingPreprocessedInput("Selected temporal case lacks saved source lineage.")
                            source_rows.append({"case_id": case_id, "condition_id": condition_id,
                                                "run_id": str(run_id), "node_id": expanded,
                                                "source_key": source.source_key,
                                                "source_observation_ids": source.source_observation_ids,
                                                "source_replicate_ids": source.source_replicate_ids})
                    row["case_complete"] = all(pd.notna(row.get(node)) for node in expanded_ids)
                    rows.append(row)
    if not rows:
        raise MissingPreprocessedInput("Selected KER time contexts have no jointly prepared cases.")
    network = {"status": "prepared", "mode": mode, "nodes": graph_nodes, "edges": edges,
               "time_grid": sorted({time for _, time in node_times}), "seed": seed,
               "n_realizations": run_count, "selected_context_count": len(contexts)}
    return (pd.DataFrame(rows, columns=_CASE_COLUMNS + ["time_context_id"] + expanded_ids),
            pd.DataFrame(source_rows, columns=_SOURCE_TUPLE_COLUMNS), network,
            pd.DataFrame(edges))

def _expanded_temporal_view(cases: pd.DataFrame, tuples: pd.DataFrame,
                            base_nodes: list[str], labels: dict[str, str],
                            path: pd.DataFrame | None, connection: dict[str, Any],
                            time_unit: str | None) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Join saved within-time tuples only under an explicit Step 2 time rule."""
    mode = connection.get("mode", "same_time_only")
    if mode != "hypothetical_cross_time":
        raise MissingPreprocessedInput(f"Unsupported Step 2 temporal connection: {mode}.")
    if path is None or path.empty:
        raise MissingPreprocessedInput("Step 2 temporal connection needs the uploaded AOP path structure.")
    if cases.empty or cases["fit_time"].isna().any():
        raise MissingPreprocessedInput("Step 2 temporal connection needs prepared cases at recorded effective times.")
    times = sorted(float(time) for time in cases["fit_time"].unique())
    if len(times) < 2:
        raise MissingPreprocessedInput("Step 2 temporal connection needs at least two selected effective times.")
    selected_grid = connection.get("time_grid")
    if selected_grid is not None and sorted(float(time) for time in selected_grid) != times:
        raise MissingPreprocessedInput("Step 2 temporal time grid must equal the selected Step 2 times.")
    if not connection.get("selected_contexts"):
        raise MissingPreprocessedInput("The demo needs the selected time context.")
    return _expanded_selected_contexts(cases, tuples, base_nodes, labels, path,
                                       connection, time_unit)


def build_prepared_view_artifacts(
    raw: pd.DataFrame, processed: pd.DataFrame, options: dict[str, Any],
    existing_artifacts: dict[str, Any], *, generation: str,
    path: pd.DataFrame | None = None, source_session_id: str | None = None,
) -> dict[str, Any]:
    """Describe already prepared values and cases without changing either one."""
    for field, label in (("fit_time", "effective sampling time"),
                         ("fit_dose", "effective dose")):
        if field in processed and pd.to_numeric(processed[field], errors="coerce").isin(
                [float("inf"), float("-inf")]).any():
            raise MissingPreprocessedInput(
                f"Step 2 requires a finite {label}; correct the selected source data before publishing."
            )
    metadata = existing_artifacts.get("bn_case_metadata", {})
    if isinstance(metadata, dict) and metadata.get("node_ids"):
        node_ids = [str(node) for node in metadata["node_ids"]]
    else:
        node_key = "endpoint_id" if "endpoint_id" in processed else "ke_id"
        node_ids = sorted({str(value) for value in processed[node_key].dropna()}) if node_key in processed else []
    ledger = _source_ledger(raw, options)
    scope = _scope(ledger)
    selected_indices = ledger.loc[ledger["selected_for_preparation"], "source_row_index"].astype(int).tolist()
    selected_path = _path_in_scope(path, scope)
    scoped_path = selected_path.copy() if selected_path is not None else None
    chosen = options.get("coverage_selection")
    selection_result = None
    if chosen is not None:
        from app.upstream.coverage_selection import (
            build_selection_candidates, edge_ids_for_rows, validate_source_issue_resolutions,
        )

        try:
            validate_source_issue_resolutions(raw, path, chosen)
            selection_result = build_selection_candidates(processed, selected_path, chosen)
        except ValueError as exc:
            raise MissingPreprocessedInput(str(exc)) from exc
        selected_ids = set(selection_result["selected_ker_ids"])
        selected_path = selected_path.iloc[
            [index for index, edge_id in enumerate(edge_ids_for_rows(selected_path))
             if edge_id in selected_ids]
        ].copy()
        selected_nodes = set(selected_path["from_ke_id"]) | set(selected_path["to_ke_id"])
        node_ids = [node for node in node_ids if node in selected_nodes] if selected_ids else node_ids
        if selected_ids and not node_ids:
            raise MissingPreprocessedInput("Selected KERs have no prepared KE nodes.")
    labels = _resolved_labels(raw.iloc[selected_indices], selected_path, node_ids)
    namespace = {column: scope[column] for column in ("dataset_id", "stressor", "experiment_id")}
    ledger["condition_id"] = [
        "condition_" + _digest_parts(namespace, float(dose))[:24]
        if selected and pd.notna(dose) and np.isfinite(float(dose)) else None
        for selected, dose in zip(ledger["selected_for_preparation"], ledger["fit_dose"])
    ]
    pairing_mode = metadata.get("pairing_mode", options.get("replicates", "retain_replicates")) if isinstance(metadata, dict) else options.get("replicates", "retain_replicates")
    cases, tuples = _cases_and_sources(existing_artifacts, node_ids, scope, pairing_mode)
    used = set()
    for encoded in tuples["source_observation_ids"]:
        try:
            source_ids = json.loads(encoded)
        except (TypeError, ValueError) as exc:
            raise MissingPreprocessedInput("Step 2 source tuple has invalid observation IDs.") from exc
        if not isinstance(source_ids, list):
            raise MissingPreprocessedInput("Step 2 source tuple observation IDs must be a list.")
        used.update(source_ids)
    if not used.issubset(set(ledger.loc[ledger["selected_for_preparation"], "source_observation_id"])):
        raise MissingPreprocessedInput("Step 2 case lineage references a source row outside the selected upload.")
    ledger["used_in_cases"] = ledger["source_observation_id"].isin(used)
    if pairing_mode != "repeated_one_to_one":
        raise MissingPreprocessedInput("The demo uses repeated one-to-one pairing.")
    view_kind = "constructed_cases"
    mixed_source_fields = _mixed_source_fields(ledger)
    if not cases.empty and view_kind != "observed_rows" and mixed_source_fields:
        raise MissingPreprocessedInput(
            "Step 2 cannot combine multiple " + ", ".join(mixed_source_fields) +
            " source identities in one prepared case view. Select one compatible "
            "analysis series using source filters in shared Step 2 before fitting."
        )
    connection = options.get("temporal_connection") or {}
    if selection_result is None or not selection_result["contexts"]:
        raise MissingPreprocessedInput("The demo needs the selected KER time context.")
    connection = {**connection, "selected_contexts": selection_result["contexts"]}
    base_node_ids = node_ids.copy()
    cases, tuples, temporal_network, temporal_edges = _expanded_temporal_view(
        cases, tuples, base_node_ids, labels, selected_path, connection,
        _effective_units(processed, ledger, base_node_ids, options)[0]["time_unit"],
    )
    node_ids = [node["id"] for node in temporal_network["nodes"]]
    labels = {node["id"]: node["label"] for node in temporal_network["nodes"]}
    temporal_policy = connection["mode"]
    used = {source_id for encoded in tuples["source_observation_ids"]
            for source_id in json.loads(encoded)}
    if not used.issubset(set(ledger.loc[ledger["selected_for_preparation"], "source_observation_id"])):
        raise MissingPreprocessedInput("Step 2 temporal lineage references a source outside the selected upload.")
    ledger["used_in_cases"] = ledger["source_observation_id"].isin(used)
    recipe_hash = _digest_parts(options)
    selected_doses = sorted(float(value) for value in pd.to_numeric(
        processed.get("fit_dose", pd.Series(dtype=float)), errors="coerce").dropna().unique())
    selected_times = sorted(float(value) for value in pd.to_numeric(
        processed.get("fit_time", pd.Series(dtype=float)), errors="coerce").dropna().unique())
    selected_conditions = [{"condition_id": "condition_" + _digest_parts(namespace, dose)[:24], "fit_dose": dose}
                           for dose in selected_doses]
    source_compatibility = "single_source_scope" if not mixed_source_fields else "mixed_source_scope"
    effective_units, base_node_units = _effective_units(processed, ledger, base_node_ids, options)
    node_units = ({node["id"]: base_node_units[node["node_id"]] for node in temporal_network["nodes"]}
                  if temporal_network["status"] == "prepared" else base_node_units)
    from app.upstream.preprocess import _is_control, _standardize_columns
    selected_source = _standardize_columns(raw).reset_index(drop=True).iloc[selected_indices]
    control_mask = _is_control(selected_source, options.get("control_dose", 0.0))
    control_indices = selected_source.index[control_mask]
    source_control_conditions = set(ledger.loc[control_indices, "condition_id"].dropna())
    measured_control_indices = [index for index in control_indices if ledger.loc[index, "used_in_cases"]]
    reference_control_source_ids = ledger.loc[measured_control_indices, "source_observation_id"].tolist()
    reference_conditions = set(ledger.loc[measured_control_indices, "condition_id"].dropna())
    assessment_reason = _fold_replay_reason(
        view_kind, pairing_mode, options, source_compatibility, selected_conditions, reference_conditions,
    )
    if selection_result is not None:
        if not selection_result["contexts"]:
            assessment_reason = "No KER graph was selected in Step 2."
    if assessment_reason is None and ledger.loc[ledger["selected_for_preparation"], "condition_id"].isna().any():
        assessment_reason = "Selected source rows lack a finite effective treatment condition."
    from app.upstream.coverage_selection import audit_effective_matching, edge_ids_for_rows
    matching_audit = (audit_effective_matching(processed, scoped_path,
                      (options.get("temporal_connection") or {}).get("mode", "same_time_only"))
                      if selection_result is not None else {})
    scoped_ker_ids = set()
    if scoped_path is not None and not scoped_path.empty:
        scoped_ker_ids = set(edge_ids_for_rows(scoped_path))
    core = {
        "schema_version": SCHEMA_VERSION,
        "source_session_id": source_session_id,
        "view_kind": view_kind,
        "pairing_mode": pairing_mode,
        "temporal_policy": temporal_policy,
        "temporal_network": temporal_network,
        "candidate_contexts": selection_result["contexts"] if selection_result is not None else [],
        "graph_selection_reviewed": selection_result is not None,
        "source_issue_resolutions": chosen.get("issue_resolutions", {}) if selection_result is not None else {},
        "matching_audit": matching_audit,
        "selected_ker_ids": selection_result["selected_ker_ids"] if selection_result is not None else [],
        "excluded_ker_ids": sorted(scoped_ker_ids - set(selection_result["selected_ker_ids"]))
        if selection_result is not None else [],
        "base_node_ids": base_node_ids,
        "node_ids": node_ids,
        "node_labels": labels,
        "path_edge_count": len({(edge["from"], edge["to"])
                                for edge in temporal_network["edges"]})
        if temporal_network["status"] == "prepared" else
        _selected_path_edge_count(selected_path, base_node_ids),
        "source_scope": scope,
        "source_compatibility": source_compatibility,
        "effective_units": effective_units,
        "node_units": node_units,
        "selected_conditions": selected_conditions,
        "assessment_condition_ids": [item["condition_id"] for item in selected_conditions
                                     if item["condition_id"] in set(cases["condition_id"])
                                     and item["condition_id"] not in source_control_conditions],
        "reference_control_source_ids": reference_control_source_ids,
        "reference_control_condition_ids": sorted(reference_conditions),
        "selected_grid": {
            "doses": selected_doses,
            "times": selected_times,
            "time_unit": metadata.get("effective_units", {}).get("time", [None])[0]
            if isinstance(metadata, dict) and len(metadata.get("effective_units", {}).get("time", [])) == 1 else None,
        },
        "recipe_hash": recipe_hash,
        "recipe": json.loads(json.dumps(options, allow_nan=False)),
        "assessment_ready": assessment_reason is None,
        "assessment_unavailable_reason": assessment_reason,
        "case_count": len(cases),
        "complete_case_count": int(cases["case_complete"].sum()) if not cases.empty else 0,
        "source_count": int(ledger["selected_for_preparation"].sum()),
    }
    if options.get("source_replicate_handling"):
        core["source_replicate_handling"] = options["source_replicate_handling"]
        core["derived_fit_view"] = options["derived_fit_view"]
    if options.get("derived_fit_unavailable_reason"):
        core["derived_fit_unavailable_reason"] = options["derived_fit_unavailable_reason"]
    core["view_digest"] = _digest_parts(core, cases, tuples, ledger, temporal_edges)
    manifest = {**core, "generation": generation}
    return {
        "prepared_view_manifest": manifest,
        "prepared_view_cases": cases,
        "prepared_view_source_tuples": tuples,
        "prepared_view_source_ledger": ledger,
        "prepared_view_temporal_edges": temporal_edges,
    }

def prepare_fold(raw: pd.DataFrame, path: pd.DataFrame, view: dict[str, Any], held_out_condition_id: str,
                 excluded_training_condition_ids: tuple[str, ...] = ()) -> dict[str, Any]:
    """Replay a safe Step 2 recipe separately on training and held-out sources."""
    manifest = view["manifest"]
    valid = {item["condition_id"] for item in manifest["selected_conditions"]}
    if held_out_condition_id not in valid:
        raise MissingPreprocessedInput("Unknown Step 2 held-out condition ID; reload the selected preparation.")
    excluded = tuple(sorted(excluded_training_condition_ids))
    if (len(excluded) != len(set(excluded)) or any(item not in valid for item in excluded) or
            held_out_condition_id in excluded):
        raise MissingPreprocessedInput("Nested fold exclusions must be distinct other selected condition IDs.")
    if not manifest["assessment_ready"]:
        raise MissingPreprocessedInput(
            "Step 2 fold replay is unavailable: " + manifest["assessment_unavailable_reason"]
        )
    from app.upstream.preprocess import preprocess_dataframe
    from app.upstream.profile import ReplicateConfig

    selected = view["source_ledger"].loc[
        lambda frame: frame["selected_for_preparation"] & ~frame["condition_id"].isin(excluded)
    ]
    raw = raw.reset_index(drop=True)
    if len(raw) != len(view["source_ledger"]):
        raise MissingPreprocessedInput("Step 2 source ledger does not match the current raw upload.")
    raw["_bn_observation_id"] = [f"observation_{i}" for i in range(len(raw))]
    options = manifest["recipe"]
    scope = manifest["source_scope"]
    normalization = options.get("normalization", "none")
    if normalization != "none" or manifest["pairing_mode"] != "repeated_one_to_one":
        raise MissingPreprocessedInput("The demo folds redraw the recorded pairing on the original response scale.")
    reference_control_ids = []
    fold_id = "loco_" + _digest_parts(manifest["generation"], held_out_condition_id, excluded)[:24]
    original_seed = (options.get("replicate_config") or {}).get("bootstrap_seed", 42)
    pairing_seeds = {
        side: int(_digest_parts(
            "step2_fold_pairing_v1", original_seed, held_out_condition_id, excluded, side,
        )[:16], 16) % (2 ** 32)
        for side in ("training", "test")
    }
    if pairing_seeds["training"] == pairing_seeds["test"]:
        pairing_seeds["test"] = (pairing_seeds["test"] + 1) % (2 ** 32)
    partitions = {}
    partition_sources = {}
    for side, condition_mask in (
        ("training", selected["condition_id"].ne(held_out_condition_id)),
        ("test", selected["condition_id"].eq(held_out_condition_id)),
    ):
        source_ids = selected.loc[condition_mask, "source_observation_id"].tolist()
        if not source_ids:
            raise MissingPreprocessedInput(f"Step 2 {side} fold has no selected source observations.")
        indices = selected.loc[condition_mask, "source_row_index"].astype(int).tolist()
        side_raw = raw.iloc[indices].copy()
        artifacts: dict[str, Any] = {}
        config = options.get("replicate_config")
        if not isinstance(config, dict):
            raise MissingPreprocessedInput("Step 2 stochastic fold needs its saved pairing configuration.")
        config = {**config, "bootstrap_seed": pairing_seeds[side]}
        try:
            replayed, _, _ = preprocess_dataframe(
                side_raw,
                source_filters=options.get("source_filters"),
                ke_ids=options.get("ke_ids"), time_points=options.get("time_points"),
                doses=options.get("doses"), dose_min=options.get("dose_min"),
                dose_max=options.get("dose_max"),
                keep_control=options.get("keep_control", True),
                control_dose=options.get("control_dose", 0.0),
                missing=options.get("missing", "retain_na"),
                replicates=options.get("replicates", "retain_replicates"),
                replicate_config=ReplicateConfig.model_validate(config) if config else None,
                harmonization=options.get("harmonization", "keep_original"),
                normalization=options.get("normalization", "none"),
                log_dose=options.get("log_dose", False),
                log_response=options.get("log_response", False),
                reverse_direction=options.get("reverse_direction", False),
                bn_artifacts=artifacts, preserve_source_ids=True,
                strict_source_units=options.get("source_module") not in (None, "BN"),
            )
        except (ValueError, TypeError) as exc:
            raise MissingPreprocessedInput(f"Step 2 {side} fold replay failed: {exc}") from exc
        if replayed.empty:
            raise MissingPreprocessedInput(f"Step 2 {side} fold replay produced no model input.")
        base_nodes = manifest.get("base_node_ids", manifest["node_ids"])
        cases, tuples = _cases_and_sources(artifacts, base_nodes, scope, manifest["pairing_mode"])
        connection = options.get("temporal_connection") or {}
        if connection.get("mode") not in (None, "same_time_only"):
            if manifest.get("graph_selection_reviewed"):
                connection = {**connection, "selected_contexts": manifest["candidate_contexts"]}
            base_labels = {node["node_id"]: node["label"]
                           for node in manifest["temporal_network"]["nodes"]}
            cases, tuples, _, _ = _expanded_temporal_view(
                cases, tuples, base_nodes, base_labels,
                _path_in_scope(path, scope), connection,
                manifest["effective_units"]["time_unit"],
            )
        if cases.empty or (not manifest.get("candidate_contexts") and not cases["case_complete"].any()):
            raise MissingPreprocessedInput(f"Step 2 {side} fold has no complete prepared cases.")
        mapped = {case_id: f"{fold_id}:{side}:{case_id}" for case_id in cases["case_id"]}
        cases["case_id"] = cases["case_id"].map(mapped)
        tuples["case_id"] = tuples["case_id"].map(mapped)
        used = {source_id for encoded in tuples["source_observation_ids"]
                for source_id in json.loads(encoded)}
        if not used.issubset(set(source_ids)):
            raise MissingPreprocessedInput("Step 2 fold replay crossed the held-out source boundary.")
        original = view["cases"].loc[view["cases"]["condition_id"].isin(set(cases["condition_id"]))]
        # The same selected algorithm is redrawn within each source-disjoint
        # fold. Its case combinations need not equal the full-data draw.
        expected_conditions = set(original["condition_id"])
        if set(cases["condition_id"]) != expected_conditions or not set(
            cases["fit_dose"].dropna()
        ).issubset(set(original["fit_dose"].dropna())) or not set(
            cases["fit_time"].dropna()
        ).issubset(set(original["fit_time"].dropna())):
            raise MissingPreprocessedInput("Step 2 stochastic fold changed the selected condition grid.")
        partitions[side] = cases
        partition_sources[side] = sorted(used)
    if set(partition_sources["training"]) & set(partition_sources["test"]):
        raise MissingPreprocessedInput("Step 2 fold replay has training and test source overlap.")
    return {
        "manifest": manifest,
        "training_cases": partitions["training"],
        "test_cases": partitions["test"],
        "fold_id": fold_id,
        "held_out_condition_id": held_out_condition_id,
        "excluded_training_condition_ids": list(excluded),
        "training_source_ids": partition_sources["training"],
        "test_source_ids": partition_sources["test"],
        "reference_control_ids": reference_control_ids,
        "learned_parameters": {"normalization": normalization, "imputation": "none",
                               "reference_control_source_ids": reference_control_ids,
                               "pairing_seed_by_partition": pairing_seeds},
    }
