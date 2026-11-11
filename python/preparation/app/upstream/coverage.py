# Extracted from qAOPweaver/app/upstream/coverage.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

from collections import defaultdict

from typing import Any

import numpy as np

import pandas as pd

from app.upstream.preprocess import _dose_factor, _normalize_unit

def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in ("nan", "none", "null", "<na>") else text

def _numeric(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None

def _dose_key(value: Any, unit: Any) -> tuple[str, float] | None:
    dose = _numeric(value)
    if dose is None:
        return None
    unit_text = _normalize_unit(_text(unit))
    factor = _dose_factor(unit_text)
    if factor is not None:
        return "M", round(dose * factor, 15)
    if not unit_text:
        return None
    return unit_text, dose

def _issue(kind: str, subject: str, message: str, recommendation: str) -> dict[str, str]:
    return {"id": f"{kind}:{subject}", "kind": kind, "subject": subject,
            "message": message, "recommendation": recommendation}

def _condition_key(row: dict[str, Any]) -> tuple[str, str, tuple[str, float]] | None:
    dose = _dose_key(row.get("dose"), row.get("dose_unit"))
    if dose is None:
        return None
    return _text(row.get("dataset_id")), _text(row.get("stressor")), dose

def _pair_time(upstream: float | None, downstream: float | None) -> bool:
    return (upstream is None and downstream is None) or (
        upstream is not None and downstream is not None and upstream <= downstream
    )

def _path_has_joint_assignment(
    edges: list[tuple[str, str, list[dict[str, Any]]]],
    condition: tuple[str, str, tuple[str, float]],
) -> bool:
    """Find one coherent node-time assignment without treating pair rows as samples."""
    options: list[tuple[str, str, list[tuple[float | None, float | None]]]] = []
    for upstream, downstream, pairs in edges:
        times = [(pair["upstream_time"], pair["downstream_time"])
                 for pair in pairs if condition in pair["_condition_keys"]]
        if not times:
            return False
        options.append((upstream, downstream, times))

    def visit(index: int, assigned: dict[str, float | None]) -> bool:
        if index == len(options):
            return True
        upstream, downstream, pairs = options[index]
        for earlier, later in pairs:
            if upstream in assigned and assigned[upstream] != earlier:
                continue
            if downstream in assigned and assigned[downstream] != later:
                continue
            next_assignment = {**assigned, upstream: earlier, downstream: later}
            if visit(index + 1, next_assignment):
                return True
        return False

    return visit(0, {})

def audit_source_coverage(ke: pd.DataFrame, path: pd.DataFrame) -> dict[str, Any]:
    """Describe observed source coverage without selecting or changing any data."""
    empty = {"summary": {"ke_count": 0, "ker_count": 0, "path_count": 0,
                         "issue_count": 0, "cross_time_ker_count": 0,
                         "unmatched_ker_count": 0},
             "ke_rows": [], "ker_rows": [], "path_rows": [], "issues": []}
    if ke.empty or path.empty or not {"ke_id", "dose", "response_value"}.issubset(ke.columns) or not {
        "from_ke_id", "to_ke_id"
    }.issubset(path.columns):
        return empty

    records = ke.to_dict("records")
    observations: dict[str, dict[float | None, dict[tuple[str, str, tuple[str, float]], list[dict[str, Any]]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    ke_rows: list[dict[str, Any]] = []
    issues: list[dict[str, str]] = []
    for row in records:
        node = _text(row.get("ke_id"))
        response = _numeric(row.get("response_value"))
        time = _numeric(row.get("time_point"))
        condition = _condition_key(row)
        if node and response is not None and condition is not None:
            observations[node][time][condition].append(row)
    for node in sorted({_text(row.get("ke_id")) for row in records} - {""}):
        source_rows = [row for row in records if _text(row.get("ke_id")) == node]
        missing = sum(_numeric(row.get("response_value")) is None for row in source_rows)
        measured = sum(_numeric(row.get("response_value")) is not None for row in source_rows)
        names = sorted({_text(row.get("ke_display_label")) for row in source_rows} - {""})
        times = sorted({_numeric(row.get("time_point")) for row in source_rows} - {None})
        doses = sorted({float(row["dose"]) for row in source_rows
                        if _numeric(row.get("dose")) is not None})
        ke_rows.append({"ke_id": node, "name": names[0] if names else node,
                        "source_rows": len(source_rows), "measured_rows": measured,
                        "missing_response_rows": missing, "times": times,
                        "doses": doses,
                        "dose_units": sorted({_text(row.get("dose_unit")) for row in source_rows}),
                        "time_units": sorted({_text(row.get("time_unit")) for row in source_rows}),
                        "replicate_counts": [
                            {"time": time, "dose_unit": condition[2][0],
                             "effective_dose": condition[2][1], "n": len(values)}
                            for time, conditions in observations[node].items()
                            for condition, values in conditions.items()
                        ]})
        if missing:
            issues.append(_issue("missing_response", node, f"{missing} source response(s) are missing or nonfinite.",
                                 "Review missing-data handling in Step 2; do not count these as measured endpoints."))
        if any(_numeric(row.get("dose")) is not None and not _text(row.get("dose_unit"))
               for row in source_rows):
            issues.append(_issue("missing_dose_unit", node,
                                 "One or more dose values have no declared unit; numeric equality cannot verify a KER match.",
                                 "Provide a declared source dose unit or exclude the uncertain condition in Step 2."))
    missing_time_units = [row for row in records
                          if _numeric(row.get("time_point")) is not None
                          and not _text(row.get("time_unit"))]
    if missing_time_units:
        affected = len({_text(row.get("ke_id")) for row in missing_time_units})
        issues.append(_issue("missing_time_unit", "Sampling times",
                             f"Sampling-time unit is not declared for {len(missing_time_units)} source row(s) across {affected} KE(s).",
                             "Declare one missing time unit in Step 2, or supply source units if groups use different units."))

    all_cells: dict[tuple[str, str, float | None], set[tuple[str, float]]] = defaultdict(set)
    for by_time in observations.values():
        for time, by_condition in by_time.items():
            for dataset, stressor, dose in by_condition:
                all_cells[(dataset, stressor, time)].add(dose)
    for ke_row in ke_rows:
        node = ke_row["ke_id"]
        scopes = {(_text(row.get("dataset_id")), _text(row.get("stressor")))
                  for row in records if _text(row.get("ke_id")) == node}
        missing_cells = []
        for dataset, stressor in sorted(scopes):
            sampled_times = {_numeric(row.get("time_point")) for row in records
                             if _text(row.get("ke_id")) == node and
                             (_text(row.get("dataset_id")), _text(row.get("stressor"))) == (dataset, stressor)}
            for time in sorted(sampled_times, key=lambda value: (value is None, value or 0)):
                actual = {condition[2] for condition in observations[node].get(time, {})
                          if condition[:2] == (dataset, stressor)}
                for dose in sorted(all_cells[(dataset, stressor, time)] - actual):
                    missing_cells.append({"dataset_id": dataset, "stressor": stressor,
                                          "time": time, "dose_unit": dose[0], "effective_dose": dose[1]})
        ke_row["missing_condition_cells"] = missing_cells
        if missing_cells:
            issues.append(_issue("missing_ke_condition_cells", node,
                                 f"{len(missing_cells)} measured dose condition cell(s) are absent at this KE's sampled time(s).",
                                 "Review which matched KER conditions to retain in Step 2; absent responses are not imputed."))

    edge_map: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for index, row in enumerate(path.to_dict("records")):
        upstream, downstream = _text(row.get("from_ke_id")), _text(row.get("to_ke_id"))
        edge_id = _text(row.get("edge_id", row.get("ker_id", f"KER_{index + 1}")))
        edge_map[(edge_id, upstream, downstream)].add(_text(row.get("path_id")) or "<unspecified>")

    ker_rows: list[dict[str, Any]] = []
    pair_lookup: dict[str, tuple[str, str, list[dict[str, Any]]]] = {}
    for (edge_id, upstream, downstream), memberships in sorted(edge_map.items()):
        pairs: list[dict[str, Any]] = []
        upstream_times = sorted(observations[upstream], key=lambda value: (value is not None, value or 0))
        downstream_times = sorted(observations[downstream], key=lambda value: (value is not None, value or 0))
        for earlier in upstream_times:
            for later in downstream_times:
                if not _pair_time(earlier, later):
                    continue
                before = observations[upstream][earlier]
                after = observations[downstream][later]
                common = sorted(set(before) & set(after))
                common_doses = sorted({condition[2][1] for condition in common})
                before_only = sorted({condition[2][1] for condition in set(before) - set(after)})
                after_only = sorted({condition[2][1] for condition in set(after) - set(before)})
                counts = {f"{condition[2][1]:g}": {"upstream": len(before[condition]),
                                                  "downstream": len(after[condition])}
                          for condition in common}
                pairs.append({"upstream_time": earlier, "downstream_time": later,
                              "common_doses": common_doses,
                              "upstream_only_doses": before_only,
                              "downstream_only_doses": after_only,
                              "replicate_counts": counts,
                              "matched_condition_count": len(common),
                              "_condition_keys": common})
                if before_only or after_only:
                    issues.append(_issue("dose_mismatch", f"{edge_id}:{earlier}:{later}",
                                         "The endpoints do not share every dose/condition cell.",
                                         "Keep common conditions or choose another matching action in Step 2."))
                if any(value["upstream"] != value["downstream"] for value in counts.values()):
                    issues.append(_issue("replicate_count_mismatch", f"{edge_id}:{earlier}:{later}",
                                         "Matched conditions have unequal source replicate counts.",
                                         "Review the Step 2 replicate policy and matched-case coverage."))
        visible_pairs = [{key: value for key, value in pair.items() if key != "_condition_keys"}
                         for pair in pairs]
        status = "matched" if any(pair["matched_condition_count"] for pair in pairs) else "unmatched"
        if status == "unmatched":
            before_units = {_dose_key(row.get("dose"), row.get("dose_unit"))[0]
                            for time in observations[upstream].values() for rows in time.values()
                            for row in rows if _dose_key(row.get("dose"), row.get("dose_unit"))}
            after_units = {_dose_key(row.get("dose"), row.get("dose_unit"))[0]
                           for time in observations[downstream].values() for rows in time.values()
                           for row in rows if _dose_key(row.get("dose"), row.get("dose_unit"))}
            kind = "dose_unit_mismatch" if before_units and after_units and not (before_units & after_units) else "unmatched_ker"
            issues.append(_issue(kind, edge_id, "No measured condition matches both KER endpoints.",
                                 "Review dose units, dose values and time connections in Step 2."))
        ker_rows.append({"ker_id": edge_id, "from_ke_id": upstream, "to_ke_id": downstream,
                         "path_ids": sorted(memberships), "status": status,
                         "candidate_time_pairs": visible_pairs})
        pair_lookup[edge_id] = upstream, downstream, pairs

    path_rows: list[dict[str, Any]] = []
    for path_id in sorted({membership for members in edge_map.values() for membership in members}):
        ids = [edge_id for (edge_id, _, _), members in edge_map.items() if path_id in members]
        edges = [pair_lookup[edge_id] for edge_id in ids]
        conditions = set.intersection(*[
            {condition for pair in pairs for condition in pair["_condition_keys"]}
            for _, _, pairs in edges
        ]) if edges else set()
        supported = any(_path_has_joint_assignment(edges, condition) for condition in conditions)
        path_rows.append({"path_id": path_id, "ker_ids": sorted(set(ids)),
                          "joint_condition_count": sum(_path_has_joint_assignment(edges, condition)
                                                       for condition in conditions),
                          "status": "candidate" if supported else "unmatched"})
        if not supported:
            issues.append(_issue("unmatched_path", path_id,
                                 "No one condition has a coherent measured time assignment for this path.",
                                 "Review its KERs and conditions in Step 2."))

    return {"summary": {"ke_count": len(ke_rows), "ker_count": len(ker_rows),
                        "path_count": len(path_rows), "issue_count": len(issues),
                        "cross_time_ker_count": sum(any(
                            pair["upstream_time"] != pair["downstream_time"] and pair["matched_condition_count"]
                            for pair in row["candidate_time_pairs"]) for row in ker_rows),
                        "unmatched_ker_count": sum(row["status"] == "unmatched" for row in ker_rows)},
            "ke_rows": ke_rows, "ker_rows": ker_rows, "path_rows": path_rows,
            "issues": issues}
