"""Run the demo's continuous POD and response-response ML calculations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

LABELS = {
    "KE1446": "Decrease, Coupling of OXPHOS",
    "KE1771": "Decrease, ATP pool",
    "KE1821": "Decrease, Cell proliferation",
    "KE1521": "Decrease, Growth",
}
DAYS = (2, 4, 7)
DOSES = (0, 125, 250, 500, 1000, 2000)
SOURCE_REVISION = "0143dd56dd21a5268b17ecbe8383647b1c6323c9"


def measured_rows(path):
    df = pd.read_csv(path, float_precision="round_trip")
    required = {"ke_id", "ke_display_label", "fit_dose", "fit_time", "fit_response"}
    if required - set(df.columns):
        raise ValueError(f"Missing effective-response columns: {sorted(required - set(df.columns))}")
    for column in ("fit_dose", "fit_time", "fit_response"):
        df[column] = pd.to_numeric(df[column], errors="raise")
        if not np.isfinite(df[column]).all():
            raise ValueError(f"{column} contains a nonfinite observation.")
    if len(df) != 360 or set(df["ke_id"]) != set(LABELS) or set(df["fit_time"]) != set(DAYS):
        raise ValueError("Expected the 360 measured demo observations and the four endpoints.")
    if not df["ke_id"].map(LABELS).eq(df["ke_display_label"]).all():
        raise ValueError("Unexpected demo endpoint labels.")
    groups = df.groupby(["ke_id", "fit_time", "fit_dose"], sort=True)
    if len(groups) != 72 or not groups.size().eq(5).all() or set(df["fit_dose"]) != set(DOSES):
        raise ValueError("Expected five observations at each of 72 endpoint/time/concentration combinations.")
    return df


def safe_json(value):
    if isinstance(value, dict):
        return {key: safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, np.generic):
        return safe_json(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(safe_json(payload), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def provenance(path, packages):
    return {
        "source_revision": SOURCE_REVISION,
        "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {package: version(package) for package in packages},
    }


def fit_pod(df, path, output):
    from pod_continuous import ContinuousPODAnalyzer

    analyzer = ContinuousPODAnalyzer()
    provenance_info = provenance(path, ("numpy", "pandas", "scipy", "pybmds"))
    for day in DAYS:
        day_rows = df[df["fit_time"] == day]
        results = {}
        for endpoint, label in LABELS.items():
            responses = day_rows.loc[day_rows["ke_id"] == endpoint]
            if "replicate_id" in responses:
                responses = responses.sort_values(["fit_dose", "replicate_id"])
            groups = responses.groupby("fit_dose", sort=True)["fit_response"]
            doses = sorted(groups.groups)
            means = [float(groups.get_group(dose).mean()) for dose in doses]
            stdevs = [max(float(groups.get_group(dose).std(ddof=1)), 1e-8) for dose in doses]
            counts = [int(groups.get_group(dose).count()) for dose in doses]
            results[label] = analyzer.analyze(
                doses=doses,
                responses=means,
                ns=counts,
                stdevs=stdevs,
            )
            if not results[label]["model_averaged"]:
                raise RuntimeError(f"No model average for {endpoint} on day {day}.")
        write_json(output / "POD" / f"day-{day}" / "results.json", {
            "workflow": "pod",
            "data_type": "continuous",
            "bmr_type": "std_dev",
            "bmr_value": 1.0,
            "fit_input_unit": "observation",
            "time_value": day,
            "ke_results": results,
            "provenance": {"engine": "pybmds", **provenance_info},
        })


def condition_matrix(df):
    # Apply a Series mean to each group, as in as_fit_condition_means().
    grouped = df.groupby(["ke_id", "fit_dose", "fit_time"], sort=True)["fit_response"]
    rows = [(*key, float(group.mean())) for key, group in grouped]
    means = pd.DataFrame(rows, columns=["ke_id", "fit_dose", "fit_time", "mean"])
    means["label"] = means["ke_id"].map(LABELS)
    matrix = means.pivot(index=["fit_dose", "fit_time"], columns="label", values="mean").sort_index()
    if matrix.shape != (18, 4) or matrix.isna().any().any():
        raise ValueError("Expected 18 complete concentration-time condition means per endpoint.")
    return matrix


def fit_ml(df, path, output):
    from ml_case import anchor_holdout, gaussian_process, scipy_fallback

    matrix = condition_matrix(df)
    labels = list(LABELS.values())
    provenance_info = provenance(path, ("numpy", "pandas", "scipy", "scikit-learn"))
    methods = (
        ("Gaussian_process", "gaussian_process", "scikit-learn GaussianProcessRegressor", gaussian_process),
        ("Symbolic_fallback", "symbolic_regression", "scipy_fallback", scipy_fallback),
    )
    for folder, method, library, calculate in methods:
        per_edge_metrics = {}
        per_edge_formulas = {}
        predictions_map = {}
        scatter_data = {}
        for source, target in zip(labels, labels[1:]):
            x, y = matrix[source].to_numpy(), matrix[target].to_numpy()
            fitted = calculate(x, y, source)
            key = f"{source}→{target}"
            per_edge_metrics[key] = fitted["metrics"]
            per_edge_formulas[key] = fitted["formula"]
            predictions_map[key] = [
                {**row, "x_label": source, "y_label": target} for row in fitted["predictions"]
            ]
            scatter_data[key] = [{source: float(a), target: float(b)} for a, b in zip(x, y)]

        source, target = labels[0], labels[-1]
        anchor = calculate(matrix[source].to_numpy(), matrix[target].to_numpy(), source)
        write_json(output / "ML" / folder / "results.json", {
            "method": method,
            "target_scope": "AOP",
            "fit": {
                "params": anchor["params"],
                "predictions": anchor["predictions"],
                "per_edge_metrics": per_edge_metrics,
                "per_edge_formulas": per_edge_formulas,
                "predictions_map": predictions_map,
                "scatter_data": scatter_data,
            },
            "diagnostics": {**anchor["metrics"], **anchor_holdout(matrix, calculate, source, target)},
            "provenance": {"library": library, **provenance_info},
            "input_provenance": {
                "analysis_unit": "condition_mean",
                "time_scope": "pooled_descriptive",
                "fitted_times": list(DAYS),
                "fitted_rows": len(matrix),
            },
        })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--only", choices=("pod", "ml"))
    args = parser.parse_args()
    observations = measured_rows(args.observations)
    if args.only != "ml":
        fit_pod(observations, args.observations, args.output)
    if args.only != "pod":
        fit_ml(observations, args.observations, args.output)


if __name__ == "__main__":
    main()
