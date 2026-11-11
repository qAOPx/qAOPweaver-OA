"""Check the archived qKER candidate and ODE outputs against a new run."""

import argparse
import csv
import json
import math
from pathlib import Path


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def scalar(value):
    return value[0] if isinstance(value, list) and len(value) == 1 else value


def compare_means(actual, reference):
    def keyed(path):
        with path.open(encoding="utf-8-sig", newline="") as stream:
            return {
                (row["ke_id"], float(row["fit_dose"]), float(row["fit_time"])):
                float(row["fit_response"])
                for row in csv.DictReader(stream)
            }

    computed = keyed(actual)
    archived = keyed(reference)
    assert computed.keys() == archived.keys(), f"Different condition keys in {actual}"
    delta = max(abs(computed[key] - archived[key]) for key in computed)
    assert delta < 1e-10, f"{actual}: condition means differ by up to {delta}"
    return len(computed), delta


def compare_qker(actual, reference):
    def chosen(path):
        return read_json(path)["combinations"]

    # The reported table uses these rounded values and the ranked combination IDs.
    assert chosen(actual / "results.json") == chosen(reference / "results.json")
    fresh = read_json(actual / "all-candidates" / "results.json")["candidates"]
    archived = read_json(reference / "all-candidates" / "results.json")["candidates"]
    assert fresh == archived, "qKER all-candidate rows differ from the archived run"
    return len(chosen(actual / "results.json")), len(fresh)


def compare_ode(actual, reference):
    fresh = read_json(actual / "results.json")["fitting_result"]
    archived = read_json(reference / "results.json")["fitting_result"]
    assert scalar(fresh["status"]) == scalar(archived["status"]) == "fitted"
    assert fresh["optimizer_summary"][0]["convergence_code"] == archived["optimizer_summary"][0]["convergence_code"]

    def keyed(data, section, key):
        return {tuple(row[col] for col in key): row for row in data[section]}

    deltas = {}
    for section, keys, fields in (
        ("fitted_parameters", ("parameter_id",), ("estimate", "initial_value", "lower", "upper")),
        ("fitted_predictions", ("time", "state_id"), ("observed_value", "predicted_value")),
        ("residuals", ("time", "state_id"), ("residual",)),
    ):
        new = keyed(fresh, section, keys)
        old = keyed(archived, section, keys)
        assert new.keys() == old.keys()
        deltas[section] = max(abs(new[k][f] - old[k][f]) for k in new for f in fields)
    for field in ("RMSE", "MAE", "bias", "n_observations"):
        new = fresh["fit_metrics"][0][field]
        old = archived["fit_metrics"][0][field]
        deltas[field] = abs(new - old)
    assert all(math.isfinite(value) and value < 0.001 for value in deltas.values()), deltas
    return deltas


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Extracted demo archive")
    parser.add_argument("results", type=Path, help="Fresh standalone result directory")
    args = parser.parse_args()
    reference = args.archive / "rerun" / "windows_api"
    result = {"qKER": compare_qker(args.results / "qKER", reference / "qKER")}
    for label in ("synthetic", "measured"):
        actual = args.results / "ODE" / label
        expected = reference / "ODE" / label
        result[f"ODE {label}"] = compare_ode(actual, expected)
        result[f"ODE {label} means"] = compare_means(
            actual / "Step2" / "preprocessed.csv",
            expected / "Step2" / "preprocessed.csv",
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
