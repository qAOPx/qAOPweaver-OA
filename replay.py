"""Rerun the qAOPweaver demo from its data ZIP and optionally compare reference tables."""

import argparse
from collections import Counter
import csv
from hashlib import sha256
from importlib import metadata
from importlib.util import find_spec
import io
import json
import os
from pathlib import Path, PurePosixPath
import platform
import subprocess
import sys

import pandas as pd
from zipfile import BadZipFile, ZIP_DEFLATED, ZipFile

from compare_results import Audit


ROOT = Path(__file__).resolve().parent
ZIP_HASH = "5ee65e07c019e62889edeb355374c7f82896988b9280abd084b7eacfe67e4f66"
MEASURED = "inputs/measured/lemna_cccp_ke_response_data_intake_ready.csv"
ARCHIVED_RUN = PurePosixPath("rerun/windows_api")
PREPARATION_CHECKS = {
    "measured_observations.csv": ARCHIVED_RUN / "POD/Step2/preprocessed.csv",
    **{
        name: ARCHIVED_RUN / "Prepared/Step2" / name
        for name in ("prepared_cases.csv", "temporal_edges.csv", "source_ledger.csv", "source_tuples.csv")
    },
}
ADDITIONAL_REFERENCES = {
    *PREPARATION_CHECKS.values(),
    ARCHIVED_RUN / "qKER/results.json",
    ARCHIVED_RUN / "qKER/all-candidates/results.json",
    *{
        ARCHIVED_RUN / "ODE" / scenario / file
        for scenario in ("synthetic", "measured")
        for file in ("results.json", "Step2/preprocessed.csv")
    },
}
EXPECTED = {
    ("Input reconstruction", "S1"): 216,
    ("POD", "S2"): 36,
    ("qKER", "S3"): 54,
    ("tPSEM", "S4"): 36,
    ("tPSEM", "S5"): 12,
    ("DBN", "S6"): 32,
    ("DBN", "S6A"): 16,
    ("DBN", "S6B"): 45,
    ("ML", "S7"): 21,
    ("ODE", "S8"): 15,
    ("ODE", "Fig 11/S1"): 4,
}


def extract_archive(source: Path, destination: Path) -> str:
    archive_hash = sha256(source.read_bytes()).hexdigest()
    if archive_hash != ZIP_HASH:
        raise ValueError("Use the demo data ZIP linked in the README.")
    with ZipFile(source) as archive:
        manifest = dict(
            (name, digest) for digest, name in
            (line.split("  ", 1) for line in archive.read("SHA256SUMS.txt").decode().splitlines())
        )
        if set(manifest) != set(archive.namelist()) - {"SHA256SUMS.txt"}:
            raise ValueError("Archive contents and checksum manifest differ.")
        for name, digest in manifest.items():
            member = PurePosixPath(name)
            if member.is_absolute() or ".." in member.parts or "\\" in name:
                raise ValueError(f"Unsafe archive path: {name}")
            if member.parts[0].lower() in {"manuscript", "supplement"} or member.suffix.lower() in {".doc", ".docx"}:
                raise ValueError(f"Supply reference documents separately from the demo archive. Found {name}")
            content = archive.read(name)
            if sha256(content).hexdigest() != digest:
                raise ValueError(f"Archive checksum mismatch: {name}")
            # Extract only the study inputs and archived comparison targets.
            if member.parts[0] in {"inputs", "reference"} or member in ADDITIONAL_REFERENCES:
                target = destination.joinpath(*member.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
    return archive_hash


def run(command: list, log: Path, *, r_library: Path | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "R_LIBS", "R_LIBS_SITE", "R_LIBS_USER"):
        env.pop(key, None)
    env.update({"PYTHONUTF8": "1", "PYTHONNOUSERSITE": "1", "R_LIBS_USER": "NULL", "R_LIBS_SITE": "NULL"})
    if r_library is not None:
        env["R_LIBS"] = str(r_library.resolve(strict=True))
    with log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}). See {log}.")


def linux_path(path: Path) -> str:
    path = path.resolve()
    if len(path.drive) != 2 or path.drive[1] != ":":
        raise ValueError(f"WSL POD requires a local drive path: {path}")
    return "/mnt/" + path.drive[0].lower() + "/" + Path(*path.parts[1:]).as_posix()


def pack_model_details(results: Path, prep: Path) -> None:
    # File shapes match the comparison routine for the archived exports.
    for family in ("tPSEM", "DBN"):
        for model in sorted((results / family).glob("realization-*")):
            cases = model / "prepared_cases.csv"
            frozen = model / "frozen_model.json"
            if not cases.is_file() or not frozen.is_file():
                raise FileNotFoundError(f"Missing regenerated model details in {model}")
            with (prep / "prepared_cases.csv").open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                case_fields = reader.fieldnames
                realization = model.name.removeprefix("realization-").replace("-", ":")
                scoped_cases = [row for row in reader if row["realization_id"] == realization]
            case_ids = {row["case_id"] for row in scoped_cases}
            if len(case_ids) != 30:
                raise ValueError(f"Expected 30 regenerated cases for {realization}.")
            def csv_bytes(fields, rows):
                buffer = io.StringIO(newline="")
                writer = csv.DictWriter(buffer, fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                return buffer.getvalue().encode()
            with (prep / "source_tuples.csv").open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                tuple_fields = reader.fieldnames
                tuples = [row for row in reader if row["case_id"] in case_ids]
            if len(tuples) != 120:
                raise ValueError(f"Expected 120 regenerated source tuples for {realization}.")
            with ZipFile(model / "model-details.zip", "w", ZIP_DEFLATED) as archive:
                archive.write(frozen, "frozen_model.json")
                archive.writestr("prepared_cases.csv", csv_bytes(case_fields, scoped_cases))
                archive.writestr("source_tuples.csv", csv_bytes(tuple_fields, tuples))
                archive.write(prep / "source_ledger.csv", "source_ledger.csv")


def verify_preparation(archive: Path, prep: Path) -> dict:
    result = {}
    for name, archived in PREPARATION_CHECKS.items():
        rebuilt = pd.read_csv(prep / name)
        reference = pd.read_csv(archive / archived)
        if name == "source_tuples.csv":
            keys = ["case_id", "node_id"]
            rebuilt = rebuilt.sort_values(keys).reset_index(drop=True)
            reference = reference.sort_values(keys).reset_index(drop=True)
        pd.testing.assert_frame_equal(rebuilt, reference, check_dtype=False, rtol=1e-12, atol=1e-12)
        result[name] = len(rebuilt)
    fold_cases = pd.read_csv(prep / "fold_cases.csv")
    if len(fold_cases) != 450 or fold_cases["held_out_condition_id"].nunique() != 5:
        raise ValueError("The five treatment-held-out folds were not rebuilt.")
    result["source_disjoint_fold_cases"] = len(fold_cases)
    return result


def verify_paper(archive: Path, results: Path, output: Path, si: bytes) -> dict:
    audit = Audit(argparse.Namespace(
        si=io.BytesIO(si),
        hosted_reference=archive / "reference" / "hosted",
        runs=results,
        out=output / "comparison",
        linux_pod=None,
        measured_input=archive / MEASURED,
    ))
    for operation in (audit.check_input, audit.pod, audit.qker, audit.psem, audit.dbn, audit.ml, audit.ode):
        operation()
    audit.write()
    counts = Counter((row["module"], row["reference"]) for row in audit.rows)
    if counts != EXPECTED:
        raise ValueError(f"Unexpected reference comparison counts {dict(counts)}")
    matches = Counter((row["module"], row["reference"]) for row in audit.rows
                      if row["verdict"] == "matches_display")
    summary = [
        {"module": module, "reference": reference, "checks": total,
         "at_printed_precision": matches[module, reference]}
        for (module, reference), total in sorted(counts.items())
    ]
    return {"checks": sum(counts.values()), "at_printed_precision": sum(matches.values()), "by_table": summary}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Local copy of the Zenodo demo ZIP")
    comparison = parser.add_mutually_exclusive_group(required=True)
    comparison.add_argument("--si", type=Path, help="Separate Word document with the demo reference tables")
    comparison.add_argument("--fit-only", action="store_true", help="Run all six workflows and archive checks using only the demo ZIP")
    parser.add_argument("--out", type=Path, required=True, help="New directory for results and logs")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable with public dependencies only")
    parser.add_argument("--r-library", type=Path, help="Optional isolated library restored from renv.lock")
    parser.add_argument("--pod-python", help="Optional Python executable for POD with the Linux dependency lock")
    parser.add_argument("--pod-wsl-python", help="Optional Python executable inside a WSL distribution for Linux POD")
    parser.add_argument("--wsl-distribution", default="Ubuntu")
    args = parser.parse_args()
    if args.pod_python and args.pod_wsl_python:
        parser.error("Select at most one separate POD Python executable.")

    if find_spec("app") is not None or any(
        dist.metadata["Name"].casefold() in {"qaopweaver", "qaopx"}
        for dist in metadata.distributions()
    ):
        raise RuntimeError("Use a Python environment containing only the public dependencies.")
    si = None
    if args.si is not None:
        if args.si.suffix.lower() != ".docx":
            parser.error("--si requires a separate Word document (.docx) with the demo reference tables.")
        try:
            si = args.si.read_bytes()
            with ZipFile(io.BytesIO(si)) as word_file:
                if "word/document.xml" not in word_file.namelist():
                    parser.error("--si must point to a valid Word document.")
        except (OSError, BadZipFile):
            parser.error(f"Cannot read the reference document at {args.si}.")

    output = args.out.resolve()
    source = args.archive.resolve()
    if output.exists():
        raise FileExistsError(f"Use a new output directory: {output}")
    output.mkdir(parents=True)
    archive = output / "archive"
    archive_hash = extract_archive(source, archive)
    results = output / "runs"
    prep = results / "Prepared" / "Step2"
    logs = output / "logs"
    r = [args.rscript, "--vanilla"]
    run([*r, "-e", "stopifnot(!requireNamespace('qAOPweaver', quietly=TRUE)); print(R.version); print(.libPaths())"],
        logs / "runtime-r.log", r_library=args.r_library)
    run([sys.executable, ROOT / "python" / "preparation" / "prepare.py", archive, prep],
        logs / "preparation.log")
    preparation = verify_preparation(archive, prep)
    observations = prep / "measured_observations.csv"
    fits = ROOT / "python" / "fits" / "run_fits.py"
    if args.pod_wsl_python:
        command = ["wsl.exe", "-d", args.wsl_distribution, "--", args.pod_wsl_python,
                   linux_path(fits), "--observations", linux_path(observations),
                   "--output", linux_path(results), "--only", "pod"]
    else:
        command = [args.pod_python or sys.executable, fits, "--observations", observations,
                   "--output", results, "--only", "pod"]
    run(command, logs / "pod.log")
    run([sys.executable, fits, "--observations", observations, "--output", results, "--only", "ml"],
        logs / "ml.log")
    run([*r, ROOT / "r" / "qker_ode" / "run.R", archive, results], logs / "qker-ode.log", r_library=args.r_library)
    run([*r, ROOT / "r" / "prepared" / "run.R", prep / "prepared_cases.csv",
         prep / "temporal_edges.csv", prep / "fold_cases.csv", results], logs / "tpsem-dbn.log", r_library=args.r_library)
    additional_report = output / "qker-ode-archive-comparison.json"
    run([sys.executable, ROOT / "r" / "qker_ode" / "verify.py", archive, results], additional_report)
    pack_model_details(results, prep)
    if si is None:
        summary = {
            "paper_comparison": "skipped",
            "paper_comparison_reason": "The 487 reference table checks require the separate Word document.",
        }
    else:
        summary = verify_paper(archive, results, output, si)
        summary["paper_comparison"] = "completed"
        summary["si_sha256"] = sha256(si).hexdigest()
    summary["preparation_rows_matching"] = preparation
    summary["archive_sha256"] = archive_hash
    summary["pod_platform"] = "Linux via WSL" if args.pod_wsl_python else platform.system()
    summary["python"] = platform.python_version()
    summary["python_platform"] = platform.platform()
    (output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if si is None:
        print("All six workflows and archive checks completed.")
        print(summary["paper_comparison_reason"])
    else:
        for row in summary["by_table"]:
            print(f"{row['module']} {row['reference']}: {row['at_printed_precision']}/{row['checks']}")
    print(f"Report: {output / 'validation.json'}")
    return 0 if si is None or summary["at_printed_precision"] == summary["checks"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
