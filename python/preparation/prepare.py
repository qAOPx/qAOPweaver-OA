"""Rebuild the demo observations, paired cases and held-out folds."""

import argparse
import json
from pathlib import Path

import pandas as pd

from app.upstream.preprocess import preprocess_dataframe
from app.upstream.prepared_view import build_prepared_view_artifacts, prepare_fold
from app.upstream.profile import ReplicateConfig


EDGES = ("KER_KE1446_KE1771", "KER_KE1771_KE1821", "KER_KE1821_KE1521")
OPTIONS = {
    "source_module": "PSEM",
    "keep_control": True,
    "control_dose": 0.0,
    "missing": "retain_na",
    "replicates": "retain_replicates",
    "replicate_config": {"mode": "repeated_one_to_one", "n_bootstrap": 3, "bootstrap_seed": 42},
    "harmonization": "keep_original",
    "normalization": "none",
    "log_dose": False,
    "log_response": False,
    "reverse_direction": False,
    "temporal_connection": {"mode": "hypothetical_cross_time", "seed": 42, "n_realizations": 1},
    "coverage_selection": {
        "reviewed": True,
        "selected_ker_ids": list(EDGES),
        "issue_resolutions": {},
        "time_pairs": [
            {"ker_id": edge, "upstream_time": source, "downstream_time": target}
            for edge, (source, target) in zip(EDGES, ((2, 2), (2, 4), (4, 7)), strict=True)
        ],
    },
}


def rebuild(archive: Path, output: Path) -> None:
    measured = archive / "inputs" / "measured"
    raw = pd.read_csv(measured / "lemna_cccp_ke_response_data_intake_ready.csv").rename(
        columns={"ke_name": "ke_display_label", "time": "time_point", "response": "response_value"}
    )
    path = pd.read_csv(measured / "lemna_cccp_path_structure_data_intake_ready.csv").rename(
        columns={"ker_id": "edge_id"}
    )
    if len(raw) != 360 or list(path["edge_id"]) != list(EDGES):
        raise ValueError("The inputs do not match the qAOPweaver demo design.")

    output.mkdir(parents=True, exist_ok=True)
    observed, _, _ = preprocess_dataframe(raw, strict_source_units=True)
    observed.to_csv(output / "measured_observations.csv", index=False)

    artifacts = {}
    processed, _, _ = preprocess_dataframe(
        raw, replicate_config=ReplicateConfig.model_validate(OPTIONS["replicate_config"]),
        bn_artifacts=artifacts, strict_source_units=True,
    )
    prepared = build_prepared_view_artifacts(
        raw, processed, OPTIONS, artifacts, generation="lemna_cccp_v16", path=path,
    )
    view = {
        key.removeprefix("prepared_view_"): value
        for key, value in prepared.items()
    }
    for key in ("cases", "source_tuples", "source_ledger", "temporal_edges"):
        name = "prepared_cases" if key == "cases" else key
        view[key].to_csv(output / f"{name}.csv", index=False)
    (output / "prepared_view_manifest.json").write_text(
        json.dumps(view["manifest"], indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    fold_cases = []
    fold_provenance = []
    for condition_id in view["manifest"]["assessment_condition_ids"]:
        fold = prepare_fold(raw, path, view, condition_id)
        fold_provenance.append({
            "held_out_condition_id": condition_id,
            "training_source_ids": fold["training_source_ids"],
            "test_source_ids": fold["test_source_ids"],
            "pairing_seeds": fold["learned_parameters"]["pairing_seed_by_partition"],
        })
        for side in ("training", "test"):
            fold_cases.append(fold[f"{side}_cases"].assign(
                held_out_condition_id=condition_id, partition=side
            ))
    pd.concat(fold_cases, ignore_index=True).to_csv(output / "fold_cases.csv", index=False)
    (output / "fold_provenance.json").write_text(
        json.dumps(fold_provenance, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    rebuild(arguments.archive, arguments.output)
