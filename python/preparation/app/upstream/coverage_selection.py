# Extracted from qAOPweaver/app/upstream/coverage_selection.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

from typing import Any

import numpy as np

import pandas as pd

from app.upstream.preprocess import _normalize_unit, _time_factor

from app.upstream.coverage import audit_source_coverage

def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    value = str(value).strip()
    return "" if value.casefold() in ("none", "nan", "null", "<na>") else value

def _finite(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Selected KER time pairs need finite numeric times.") from exc
    if not np.isfinite(number):
        raise ValueError("Selected KER time pairs need finite numeric times.")
    return number

def _canonical_time_unit(value: Any) -> str:
    factor = _time_factor(_text(value))
    return {1.0: "h", 24.0: "day", 1 / 60: "min", 1 / 3600: "s",
            168.0: "week"}.get(factor, "")

def _selected_source(processed: pd.DataFrame) -> pd.DataFrame:
    required = {"ke_id", "fit_dose", "fit_response"}
    if not required.issubset(processed.columns):
        raise ValueError("Step 2 effective KE, dose and response fields are required for matching.")
    frame = processed.copy()
    frame["fit_dose"] = pd.to_numeric(frame["fit_dose"], errors="coerce")
    frame["fit_response"] = pd.to_numeric(frame["fit_response"], errors="coerce")
    frame = frame.loc[np.isfinite(frame["fit_dose"]) & np.isfinite(frame["fit_response"])].copy()
    if "fit_time" not in frame:
        frame["fit_time"] = np.nan
    else:
        frame["fit_time"] = pd.to_numeric(frame["fit_time"], errors="coerce")
    for field in ("dataset_id", "stressor"):
        if field not in frame:
            frame[field] = ""
        frame[field] = frame[field].map(_text)
    frame["ke_id"] = frame["ke_id"].map(_text)
    for name in ("dose_unit", "time_unit"):
        effective = f"{name}_harmonized"
        if effective in frame:
            frame[name] = frame[effective].map(_text)
        elif name in frame:
            frame[name] = frame[name].map(_text)
        else:
            frame[name] = ""
    frame["dose_unit"] = frame["dose_unit"].map(lambda value: _normalize_unit(value) if value else "")
    frame["time_unit"] = frame["time_unit"].map(_canonical_time_unit)
    return frame

def _condition_keys(rows: pd.DataFrame) -> set[tuple[str, str, str, float]]:
    return {(row.dataset_id, row.stressor, row.dose_unit, float(row.fit_dose))
            for row in rows.itertuples(index=False) if row.dose_unit}

def _condition_records(keys: set[tuple[str, str, str, float]]) -> list[dict[str, Any]]:
    return [{"dataset_id": dataset, "stressor": stressor,
             "dose_unit": unit, "dose": dose}
            for dataset, stressor, unit, dose in sorted(keys)]

def edge_ids_for_rows(path: pd.DataFrame) -> list[str]:
    """Use one KER identifier policy for matching and saved path selection."""
    return [_text(row.get("edge_id", row.get("ker_id", f"KER_{index + 1}")))
            for index, row in enumerate(path.to_dict("records"))]

def _edge_rows(path: pd.DataFrame) -> dict[str, tuple[str, str, list[str]]]:
    result: dict[str, tuple[str, str, list[str]]] = {}
    for edge_id, row in zip(edge_ids_for_rows(path), path.to_dict("records")):
        endpoints = (_text(row["from_ke_id"]), _text(row["to_ke_id"]))
        if edge_id in result and result[edge_id][:2] != endpoints:
            raise ValueError(f"KER {edge_id} has conflicting uploaded endpoints.")
        paths = result[edge_id][2] if edge_id in result else []
        path_id = _text(row.get("path_id"))
        if path_id and path_id not in paths:
            paths.append(path_id)
        result[edge_id] = (*endpoints, paths)
    return result

def validate_source_issue_resolutions(raw: pd.DataFrame, path: pd.DataFrame | None,
                                      selection: dict[str, Any]) -> dict[str, Any]:
    """Require a recorded user decision for each source finding, including exclusions."""
    if path is None or path.empty or not {"from_ke_id", "to_ke_id"}.issubset(path.columns):
        raise ValueError("Step 2 matching selection needs the uploaded KER path structure.")
    audit = audit_source_coverage(raw, path)
    decisions = selection.get("issue_resolutions")
    expected = {issue["id"] for issue in audit["issues"]}
    if not isinstance(decisions, dict) or set(decisions) != expected:
        raise ValueError("Review all Step 1 source-data findings in Step 2 before saving the selected graph.")
    selected = set(selection.get("selected_ker_ids") or [])
    edges = _edge_rows(path)
    paths = {row["path_id"]: set(row["ker_ids"]) for row in audit["path_rows"]}
    for issue in audit["issues"]:
        action = decisions[issue["id"]]
        if action == "accept_recommendation":
            continue
        if issue["kind"] == "missing_time_unit":
            raise ValueError("Missing sampling-time units need a Step 2 declaration or corrected source data.")
        if action != "exclude_affected_item":
            raise ValueError(f"Unsupported Step 2 resolution for {issue['id']}.")
        subject = issue["subject"]
        if issue["kind"] in ("dose_mismatch", "replicate_count_mismatch"):
            subject = subject.split(":", 1)[0]
        if subject in edges and subject in selected:
            raise ValueError(f"KER {subject} is selected despite its Step 2 exclusion decision.")
        if subject in paths and selected & paths[subject]:
            raise ValueError(f"Path {subject} has selected KERs despite its Step 2 exclusion decision.")
        if subject not in edges and subject not in paths:
            connected = {identifier for identifier, (upstream, downstream, _) in edges.items()
                         if subject in (upstream, downstream)}
            if selected & connected:
                raise ValueError(f"KE {subject} has selected KERs despite its Step 2 exclusion decision.")
    return audit

def audit_effective_matching(processed: pd.DataFrame, path: pd.DataFrame | None,
                             temporal_mode: str = "same_time_only") -> dict[str, Any]:
    """List eligible effective KER pairs without choosing a graph or mutating rows."""
    if path is None or path.empty or not {"from_ke_id", "to_ke_id"}.issubset(path.columns):
        return {"kers": [], "paths": [], "eligible_pair_count": 0}
    frame = _selected_source(processed)
    edges = _edge_rows(path)
    kers = []
    for edge_id, (upstream, downstream, path_ids) in edges.items():
        before = frame.loc[frame["ke_id"].eq(upstream)]
        after = frame.loc[frame["ke_id"].eq(downstream)]
        pairs = []
        reason = "No ordered measured time pair under this Step 2 time policy."
        if before.empty or after.empty:
            reason = "One or both KER endpoints have no selected measured response."
        for earlier in sorted(before["fit_time"].dropna().unique()):
            for later in sorted(after["fit_time"].dropna().unique()):
                if earlier > later or (temporal_mode == "same_time_only" and earlier != later):
                    continue
                left = before.loc[before["fit_time"].eq(earlier)]
                right = after.loc[after["fit_time"].eq(later)]
                units = set(left["time_unit"]) | set(right["time_unit"])
                if len(units) != 1 or not next(iter(units)):
                    reason = "Sampling time units are missing or incompatible; declare missing units or revise the source selection."
                    continue
                common = _condition_keys(left) & _condition_keys(right)
                if not common:
                    reason = "No shared measured dose and unit under the selected source conditions."
                    continue
                excluded = (_condition_keys(left) | _condition_keys(right)) - common
                matched_rows = sum(
                    (row.dataset_id, row.stressor, row.dose_unit, float(row.fit_dose)) in common
                    for side in (left, right) for row in side.itertuples(index=False)
                )
                pairs.append({"ker_id": edge_id, "upstream_time": float(earlier),
                              "downstream_time": float(later), "time_unit": next(iter(units)),
                              "matched_doses": sorted({key[3] for key in common}),
                              "excluded_doses": sorted({key[3] for key in excluded}),
                              "matched_conditions": _condition_records(common),
                              "excluded_conditions": _condition_records(excluded),
                              "matched_condition_count": len(common),
                              "effective_row_count": matched_rows,
                              "link_kind": "same_time" if earlier == later else temporal_mode})
        kers.append({"ker_id": edge_id, "from_ke_id": upstream,
                     "to_ke_id": downstream, "path_ids": path_ids,
                     "pairs": pairs, "status": "eligible" if pairs else "unavailable",
                     "reason": "" if pairs else reason})
    paths = []
    for path_id in sorted({identifier for row in kers for identifier in row["path_ids"]}):
        members = [row for row in kers if path_id in row["path_ids"]]
        paths.append({"path_id": path_id, "ker_ids": [row["ker_id"] for row in members],
                      "status": "candidate" if all(row["pairs"] for row in members) else "partial",
                      "unavailable_ker_ids": [row["ker_id"] for row in members if not row["pairs"]]})
    return {"kers": kers, "paths": paths,
            "eligible_pair_count": sum(len(row["pairs"]) for row in kers)}

def build_selection_candidates(
    processed: pd.DataFrame, path: pd.DataFrame, selection: dict[str, Any],
) -> dict[str, Any]:
    """Enumerate coherent selected graphs; never pool repeated source observations."""
    if not isinstance(selection, dict) or not selection.get("reviewed"):
        raise ValueError("Review and confirm the Step 1 coverage issues in Step 2 before Execute.")
    if path is None or path.empty or not {"from_ke_id", "to_ke_id"}.issubset(path.columns):
        raise ValueError("Step 2 needs the uploaded KER path to select matching contexts.")
    frame = _selected_source(processed)
    selected_ids = selection.get("selected_ker_ids")
    if not isinstance(selected_ids, list) or len(selected_ids) != len(set(selected_ids)):
        raise ValueError("Step 2 selected KER IDs must be a unique list.")
    if not selected_ids:
        if selection.get("time_pairs"):
            raise ValueError("No KER may be selected with time pairs when the graph is excluded.")
        return {"contexts": [], "selected_ker_ids": [], "candidate_count": 0}
    by_edge = {key: value[:2] for key, value in _edge_rows(path).items()}
    if not set(selected_ids).issubset(by_edge):
        raise ValueError("Selected KER ID is absent from the uploaded path.")

    requested = selection.get("time_pairs")
    if not isinstance(requested, list) or not requested:
        raise ValueError("Select one or more eligible KER time pairs in Step 2.")
    pair_map: dict[str, list[dict[str, Any]]] = {edge_id: [] for edge_id in selected_ids}
    seen_pairs = set()
    all_conditions = _condition_keys(frame)
    for item in requested:
        if not isinstance(item, dict) or item.get("ker_id") not in pair_map:
            raise ValueError("Every selected time pair must name a selected uploaded KER.")
        edge_id = item["ker_id"]
        earlier, later = _finite(item.get("upstream_time")), _finite(item.get("downstream_time"))
        if earlier > later:
            raise ValueError("A selected KER upstream time cannot be later than its downstream time.")
        signature = edge_id, earlier, later
        if signature in seen_pairs:
            raise ValueError("A selected KER time pair appears more than once.")
        seen_pairs.add(signature)
        upstream, downstream = by_edge[edge_id]
        before = frame.loc[frame["ke_id"].eq(upstream) & frame["fit_time"].eq(earlier)]
        after = frame.loc[frame["ke_id"].eq(downstream) & frame["fit_time"].eq(later)]
        if before.empty or after.empty:
            raise ValueError(f"Selected KER {edge_id} has no measured endpoints at that time pair.")
        units = set(before["time_unit"]) | set(after["time_unit"])
        if len(units) != 1 or not next(iter(units)):
            raise ValueError(f"Selected KER {edge_id} needs one comparable declared time unit.")
        matched = _condition_keys(before) & _condition_keys(after)
        if not matched:
            raise ValueError(f"Selected KER {edge_id} has no shared measured dose/condition at this time pair.")
        pair_map[edge_id].append({"ker_id": edge_id, "from_node_id": upstream,
                                  "to_node_id": downstream, "from_time": earlier,
                                  "to_time": later, "_matched_conditions": matched})
    if any(not pairs for pairs in pair_map.values()):
        raise ValueError("Every selected KER needs at least one selected eligible time pair.")

    contexts: list[dict[str, Any]] = []
    order = sorted(selected_ids, key=lambda edge_id: len(pair_map[edge_id]))

    def visit(index: int, assignment: dict[str, float], chosen: list[dict[str, Any]],
              conditions: set[tuple[Any, ...]] | None) -> None:
        if len(contexts) > 500:
            raise ValueError("Too many candidate time contexts; narrow the Step 2 time-pair selection.")
        if index == len(order):
            if not conditions:
                return
            matched = sorted({float(condition[3]) for condition in conditions})
            all_doses = sorted({float(condition[3]) for condition in all_conditions})
            contexts.append({"context_id": f"time_context_{len(contexts) + 1}",
                             "edges": [{key: value for key, value in edge.items()
                                        if key != "_matched_conditions"} for edge in chosen],
                             "node_times": dict(sorted(assignment.items())),
                             "matched_conditions": [
                                 {"dataset_id": dataset, "stressor": stressor,
                                  "dose_unit": unit, "dose": float(dose)}
                                 for dataset, stressor, unit, dose in sorted(conditions)
                             ], "matched_doses": matched,
                             "excluded_doses": [dose for dose in all_doses if dose not in matched]})
            return
        for pair in pair_map[order[index]]:
            upstream, downstream = pair["from_node_id"], pair["to_node_id"]
            earlier, later = pair["from_time"], pair["to_time"]
            if upstream in assignment and assignment[upstream] != earlier:
                continue
            if downstream in assignment and assignment[downstream] != later:
                continue
            next_conditions = pair["_matched_conditions"] if conditions is None else (
                conditions & pair["_matched_conditions"])
            if next_conditions:
                visit(index + 1, {**assignment, upstream: earlier, downstream: later},
                      [*chosen, pair], next_conditions)

    visit(0, {}, [], None)
    if not contexts:
        raise ValueError("Selected KER time pairs have no coherent path with a shared measured condition.")
    used_pairs = {(edge["ker_id"], edge["from_time"], edge["to_time"])
                  for context in contexts for edge in context["edges"]}
    unused = seen_pairs - used_pairs
    if unused:
        details = ", ".join(f"{ker} {earlier:g}→{later:g}"
                            for ker, earlier, later in sorted(unused))
        raise ValueError(f"Selected time pair(s) have no coherent context with the other chosen KERs: {details}. Review the selection.")
    return {"contexts": contexts, "selected_ker_ids": sorted(selected_ids),
            "candidate_count": len(contexts)}
