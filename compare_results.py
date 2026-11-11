"""Compare fresh qAOPweaver demo results with the reference tables."""

import argparse
import csv
import io
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from zipfile import ZipFile

from docx import Document


KE_LABELS = {
    "KE1446": "Decrease, Coupling of OXPHOS",
    "KE1771": "Decrease, ATP pool",
    "KE1821": "Decrease, Cell proliferation",
    "KE1521": "Decrease, Growth",
}
SOURCE_NODE = "KE1446@2"
TARGET_NODE = "KE1521@7"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def zip_json(path, member="results/analysis.json"):
    with ZipFile(path) as archive:
        return json.loads(archive.read(member))


def zip_csv(path, member):
    with ZipFile(path) as archive:
        return list(csv.DictReader(io.StringIO(archive.read(member).decode("utf-8-sig"))))


def compare_case_exports(observed_path, reference_path, reference_prefix="data/"):
    comparison = {}
    for name, fields in (
        ("prepared_cases.csv", ("case_id",)),
        ("source_tuples.csv", ("case_id", "node_id")),
        ("source_ledger.csv", ("source_observation_id",)),
    ):
        observed = zip_csv(observed_path, name)
        expected = zip_csv(reference_path, reference_prefix + name)
        def identifier(row):
            return tuple(row[field] for field in fields)
        lookup = {identifier(row): row for row in expected}
        assert len(lookup) == len(expected)
        comparison[name] = {
            "compared_rows": len(observed),
            "identical_rows": sum(row == lookup.get(identifier(row)) for row in observed),
        }
    return comparison


def single(value):
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def equation_terms(text):
    text = text.translate(str.maketrans({"−": "-", "²": "^2"}))
    return [float(term.replace(" ", ""))
            for term in re.findall(r"[-+]?\s*\d+\.\d+", text)]


def read_tables(path):
    tables = Document(path).tables
    output = {}
    for label, index in (("S1", 1), ("S2", 2), ("S3", 3), ("S4", 4),
                         ("S5", 5), ("S6", 6), ("S6A", 7), ("S6B", 8),
                         ("S7", 9), ("S8", 10)):
        output[label] = [[cell.text.strip() for cell in row.cells]
                         for row in tables[index].rows[1:]]
    return output


def baseline_pod(reference, day):
    hosted = reference / "POD/hosted-repeat"
    base = zip_json(reference / f"POD/day-{day}/research-results.zip")
    if day == 4:
        base["ke_results"].update(zip_json(hosted / "POD_day4_research_results.zip")["ke_results"])
    if day == 7:
        base["ke_results"].update(zip_json(hosted / "POD_day7_Growth_research_results.zip")["ke_results"])
    return base["ke_results"]

def dbn_query(contexts, evidence=None):
    """Independent exact enumeration of a mixture of complete frozen networks."""
    evidence = evidence or {}
    component_marginals = []
    component_likelihoods = []
    for context in contexts:
        model = context["bn"]
        state_lists = model["node_states"]
        nodes = list(state_lists)
        totals = {node: dict.fromkeys(state_lists[node], 0.0) for node in nodes}
        likelihood = 0.0
        cpt_lookup = {}
        for node in nodes:
            cpt = model["cpts"][node]
            parents = cpt["parent_ids"]
            rows = {
                tuple(row["parent_assignment"][p] for p in parents):
                    row["probabilities"] for row in cpt["rows"]
            }
            cpt_lookup[node] = (parents, rows)
        for states in itertools.product(*(state_lists[node] for node in nodes)):
            assignment = dict(zip(nodes, states))
            if any(assignment[node] != state for node, state in evidence.items()):
                continue
            probability = math.prod(
                rows[tuple(assignment[p] for p in parents)][assignment[node]]
                for node, (parents, rows) in cpt_lookup.items()
            )
            likelihood += probability
            for node, state in assignment.items():
                totals[node][state] += probability
        component_marginals.append({
            node: {state: mass / likelihood for state, mass in by_state.items()}
            for node, by_state in totals.items()
        })
        component_likelihoods.append(likelihood)
    mixture = {
        node: {
            state: sum(weight * marginal[node][state] for weight, marginal in
                       zip(component_likelihoods, component_marginals)) / sum(component_likelihoods)
            for state in by_state
        } for node, by_state in component_marginals[0].items()
    }
    return mixture, component_marginals, component_likelihoods


class Audit:
    def __init__(self, args):
        self.runs = args.runs
        self.reference = args.hosted_reference
        self.tables = read_tables(args.si)
        self.out = args.out
        self.rows = []
        self.notes = {}
        self.linux_pod_file = args.linux_pod
        self.measured_input = args.measured_input

    def check_input(self):
        if self.measured_input is None:
            return
        grouped = defaultdict(list)
        row_count = 0
        with self.measured_input.open(encoding="utf-8-sig", newline="") as stream:
            for row in csv.DictReader(stream):
                grouped[(row["ke_id"], float(row["time"]), float(row["dose"]))].append(
                    float(row["response"])
                )
                row_count += 1
        assert len(grouped) == len(self.tables["S1"]) == 72
        for ke, day, hour, dose, paper_mean, paper_sd, n in self.tables["S1"]:
            group = grouped[(ke, float(day), float(dose))]
            item = f"{ke} day {day}, {dose} nM"
            for metric, printed, observed in (
                ("mean_pct_control", paper_mean, mean(group)),
                ("sample_sd_pct_control", paper_sd, stdev(group)),
                ("n", n, len(group)),
            ):
                self.measure("Input reconstruction", "S1", item, metric, printed, observed)
        self.notes["Input reconstruction"] = {"rows": row_count, "condition_groups": len(grouped)}

    def measure(self, module, table, item, metric, paper, fresh, hosted=None, format_spec=None):
        baseline = paper if hosted is None else hosted
        if isinstance(fresh, (float, int)):
            paper = str(paper)
            if paper.startswith(">"):
                verdict = "matches_display" if fresh > float(paper[1:]) else "differs_from_display"
            elif paper.startswith("<"):
                verdict = "matches_display" if fresh < float(paper[1:]) else "differs_from_display"
            else:
                if format_spec is None:
                    decimals = len(paper.partition(".")[2])
                    format_spec = f".{decimals}f"
                verdict = ("matches_display" if float(format(fresh, format_spec)) == float(paper)
                           else "differs_from_display")
            delta = fresh - float(baseline) if not str(baseline).startswith(("<", ">")) else None
        else:
            verdict = "matches_display" if paper == fresh else "differs_from_display"
            delta = None
        self.rows.append({
            "module": module, "reference": table, "item": item, "metric": metric,
            "printed_v16": paper, "hosted_reference": hosted, "fresh_local": fresh,
            "delta_from_hosted_or_displayed": delta, "verdict": verdict,
        })

    def pod(self):
        fresh = {day: read_json(self.runs / f"POD/day-{day}/results.json")
                 for day in (2, 4, 7)}
        reference = {day: baseline_pod(self.reference, day) for day in (2, 4, 7)}
        for ke, day_text, *paper in self.tables["S2"]:
            day = int(day_text)
            name = KE_LABELS[ke]
            observed = fresh[day]["ke_results"][name]
            hosted = reference[day][name]
            for metric, printed in zip(("bmd", "bmdl", "bmdu"), paper):
                self.measure("POD", "S2", f"{ke} day {day}", metric, printed,
                             observed["model_averaged"][metric], hosted["model_averaged"][metric])
        notes = {}
        for day in (2, 4, 7):
            local_order = sorted(fresh[day]["ke_results"], key=lambda k: fresh[day]["ke_results"][k]["model_averaged"]["bmd"])
            hosted_order = sorted(reference[day], key=lambda k: reference[day][k]["model_averaged"]["bmd"])
            notes[f"day_{day}_bmd_rank_matches"] = local_order == hosted_order
            matched = 0
            fields = ("dose", "mean", "sd", "n")
            for name in hosted_order:
                original = reference[day][name]["observed_groups"]
                replay = fresh[day]["ke_results"][name]["observed_groups"]
                assert len(original) == len(replay)
                matched += sum(all(a[field] == b[field] for field in fields)
                               for a, b in zip(original, replay))
            notes[f"day_{day}_observed_groups_exact"] = matched
        notes["day_4_growth"] = {}
        for label, data in (("hosted", reference[4]), ("local", fresh[4]["ke_results"])):
            growth = data[KE_LABELS["KE1521"]]
            hill = next(model for model in growth["all_models"] if model["name"] == "Hill")
            notes["day_4_growth"][label] = {
                "model_average": growth["model_averaged"],
                "hill_converged": hill["converged"], "hill_flagged": hill["flagged"],
                "hill_bmd": hill["bmd"],
            }
        notes["provenance"] = fresh[4]["provenance"]
        self.notes["POD"] = notes

    def linux_pod(self):
        if self.linux_pod_file is None:
            return
        replay = read_json(self.linux_pod_file)
        reference = {day: baseline_pod(self.reference, day) for day in (2, 4, 7)}
        for ke, day_text, *paper in self.tables["S2"]:
            day = int(day_text)
            name = KE_LABELS[ke]
            observed = replay["fits"][f"{day}:{name}"]["model_averaged"]
            hosted = reference[day][name]["model_averaged"]
            for metric, printed in zip(("bmd", "bmdl", "bmdu"), paper):
                self.measure("POD Linux fit replay", "S2", f"{ke} day {day}", metric, printed,
                             observed[metric], hosted[metric])
        self.notes["POD Linux fit replay"] = {
            "scope": "Direct replay of archived endpoint groups (dose, mean, SD, n)",
            "runtime": replay["runtime"],
        }

    def qker(self):
        result = read_json(self.runs / "qKER/results.json")
        result = result.get("engine_result", result)
        combos = {single(entry["id"]): entry for entry in result["combinations"]}
        hosted = zip_json(self.reference / "qKER/research-results.zip")
        hosted = {single(entry["id"]): entry for entry in hosted["engine_result"]["combinations"]}
        for combo_day, ker, r2, aic, pval, n, equation in self.tables["S3"]:
            combo, day = re.fullmatch(r"(C\d+) \((\d+)\)", combo_day).groups()
            upstream, downstream = ker.split("→")
            def edge_record(combinations):
                return next(row for row in combinations[combo]["kers"]
                            if (single(row["from"]), single(row["to"])) == (upstream, downstream))
            value = edge_record(combos)
            ref = edge_record(hosted)
            for field, printed in (("r2", r2), ("aic", aic), ("pval", pval), ("n", n)):
                self.measure("qKER", "S3", f"{combo_day} {ker}", field, printed,
                             float(single(value[field])), float(single(ref[field])))
            self.measure("qKER", "S3", f"{combo_day} {ker}", "equation_terms",
                         equation_terms(equation), equation_terms(single(value["eq"])),
                         equation_terms(single(ref["eq"])))
            self.measure("qKER", "S3", f"{combo_day} {ker}", "time_pair",
                         f"{day} → {day}", single(value["time_pair"]))
        day4 = [combo for combo in combos.values() if str(single(combo["temporal"])).startswith("4")]
        self.notes["qKER"] = {
            "combination_count": len(combos),
            "best_ids": list(combos)[:3],
            "all_native_combinations_identical": combos == hosted,
            "mean_r2_by_top_combination": {key: single(combos[key]["r2"]) for key in ("C1", "C2", "C3")},
            "day_4_summed_aic": {single(combo["complexity"]):
                                 sum(float(single(edge["aic"])) for edge in combo["kers"]) for combo in day4},
        }

    def psem(self):
        fresh = read_json(self.runs / "tPSEM/results.json")
        hosted = zip_json(self.reference / "tPSEM/research-results.zip")
        contexts = {str(item["context"]["realization_id"]): item for item in fresh["contexts"]}
        refs = {str(item["context"]["realization_id"]): item for item in hosted["contexts"]}
        for run, link, slope, interval, p in self.tables["S4"]:
            parent, target = link.split(" → ")
            limits = interval.split(" to ")
            def edge(context):
                return next(edge for edge in context[run]["network"]["edges"]
                            if (edge["from"], edge["to"]) == (parent, target))
            value, ref = edge(contexts), edge(refs)
            v_inf, h_inf = value["path_inference"], ref["path_inference"]
            for metric, printed, observed, reference in (
                ("slope", slope, value["coefficient"], ref["coefficient"]),
                ("lower", limits[0], v_inf["confidence_interval_95"][0], h_inf["confidence_interval_95"][0]),
                ("upper", limits[1], v_inf["confidence_interval_95"][1], h_inf["confidence_interval_95"][1]),
                ("p", p, v_inf["p_value"], h_inf["p_value"]),
            ):
                self.measure("tPSEM", "S4", f"{run} {link}", metric, printed, observed, reference)
        self.notes["tPSEM"] = {
            "contexts": [{"realization": run, "cases": context["case_count"], "weight": context["weight_sum"]}
                         for run, context in contexts.items()],
        }
        for run in contexts:
            details = self.runs / f"tPSEM/realization-{run.replace(':', '-')}/model-details.zip"
            if details.exists():
                reference = self.reference / "tPSEM/research-results.zip"
                self.notes["tPSEM"].setdefault("case_exports", {})[run] = compare_case_exports(details, reference)
        for path in sorted((self.runs / "tPSEM").rglob("*.json")):
            result = read_json(path)
            if not isinstance(result, dict) or not result.get("folds") or not result.get("summary"):
                continue
            run = str(result["realization_id"])
            _, folds, n, rmse, mae = next(row for row in self.tables["S5"] if row[0] == run)
            summary = result["summary"]
            host_file = self.reference / f"tPSEM/realization-{run.replace(':', '-')}-held-out-folds.csv"
            with host_file.open(encoding="utf-8-sig", newline="") as stream:
                old_folds = list(csv.DictReader(stream))
            host_rmse = math.sqrt(sum(float(fold["rmse"]) ** 2 for fold in old_folds) / len(old_folds))
            host_mae = sum(float(fold["mae"]) for fold in old_folds) / len(old_folds)
            for field, printed, reference in (("fold_count", folds, len(old_folds)), ("scored_cases", n, 25),
                                               ("rmse", rmse, host_rmse), ("mae", mae, host_mae)):
                self.measure("tPSEM", "S5", run, field, printed, summary[field], reference)
            by_condition = {fold["held_out_condition_id"]: fold["score"] for fold in result["folds"]}
            deltas = {
                field: max(abs(float(fold[field]) -
                               by_condition[fold["held_out_condition_id"]][field])
                           for fold in old_folds)
                for field in ("rmse", "mae")
            }
            self.notes["tPSEM"].setdefault("max_fold_metric_delta", {})[run] = deltas

    def dbn(self):
        result = read_json(self.runs / "DBN/results.json")
        nodes = result["contexts"][0]["discretization"]["nodes"]
        cutpoints = {node["node_id"]: node["cutpoints"] for node in nodes}
        before, baseline_pairings, _ = dbn_query(result["contexts"])
        after, conditioned_pairings, weights = dbn_query(result["contexts"], {SOURCE_NODE: "state_1"})
        reference = read_json(self.reference / "DBN/posterior/growth-query.json")
        hosted_cutpoints = {row["node_id"]: row["cutpoints"] for row in reference["result"]["discretization"]["nodes"]}
        hosted_graph = {node["id"]: node["states"] for node in reference["graph_view"]["nodes"]}
        for node, cut1, cut2, printed_base, printed_after in self.tables["S6"]:
            for i, printed in enumerate((cut1, cut2)):
                self.measure("DBN", "S6", node, f"cutpoint_{i+1}", printed, cutpoints[node][i], hosted_cutpoints[node][i])
            for label, printed, distribution, field in (("baseline", printed_base, before, "baseline"),
                                                        ("conditioned", printed_after, after, "current")):
                for i, paper in enumerate(printed.split(" / "), start=1):
                    self.measure("DBN", "S6", node, f"{label}_state_{i}_pct", paper,
                                 100 * distribution[node][f"state_{i}"], 100 * hosted_graph[node][i - 1][field])
        notes = {
            "day_7_growth_low_baseline_pct": before[TARGET_NODE]["state_1"] * 100,
            "day_7_growth_low_conditioned_pct": after[TARGET_NODE]["state_1"] * 100,
            "pairing_low_baseline_pct": [m[TARGET_NODE]["state_1"] * 100 for m in baseline_pairings],
            "pairing_low_conditioned_pct": [m[TARGET_NODE]["state_1"] * 100 for m in conditioned_pairings],
            "evidence_likelihoods": weights,
        }
        query_file = self.runs / "DBN/query.json"
        if query_file.exists():
            query = read_json(query_file)
            query = query.get("conditioned", query.get("query", query))
            if "marginals" in query:
                notes["max_abs_native_query_minus_independent_enumeration"] = max(
                    abs(value - query["marginals"][node][state])
                    for node, states in after.items() if node in query["marginals"]
                    for state, value in states.items()
                )
        for mode, relative_path, archive in (
            ("Counts / 3 quantile", "comparison/counts-3", "alternatives/DBN_count_frontend.zip"),
            ("Parameter posterior / 3 quantile", None, "posterior/research-results.zip"),
            ("Residual only / 3 quantile", "comparison/residual-3", "alternatives/DBN_regression_frontend.zip"),
            ("Regression / 5 equal width", "comparison/residual-5", "alternatives/DBN_regression_equalwidth5_frontend.zip"),
        ):
            variant = self.runs / "DBN" / relative_path / "results.json" if relative_path else None
            if variant and not variant.exists():
                continue
            current = read_json(variant) if variant else result
            v_base, _, _ = dbn_query(current["contexts"])
            v_after, _, _ = dbn_query(current["contexts"], {SOURCE_NODE: "state_1"})
            hosted_mode = zip_json(self.reference / "DBN" / archive)
            h_base, _, _ = dbn_query(hosted_mode["contexts"])
            h_after, _, _ = dbn_query(hosted_mode["contexts"], {SOURCE_NODE: "state_1"})
            _, source, target, b, a = next(row for row in self.tables["S6A"] if row[0] == mode)
            bounds = {n["node_id"]: n["cutpoints"][0] for n in current["contexts"][0]["discretization"]["nodes"]}
            ref_bounds = {n["node_id"]: n["cutpoints"][0]
                          for n in hosted_mode["contexts"][0]["discretization"]["nodes"]}
            self.measure("DBN", "S6A", mode, "source_cutpoint", source.lstrip("<"), bounds[SOURCE_NODE], ref_bounds[SOURCE_NODE])
            self.measure("DBN", "S6A", mode, "target_cutpoint", target.lstrip("<"), bounds[TARGET_NODE], ref_bounds[TARGET_NODE])
            self.measure("DBN", "S6A", mode, "baseline_low_pct", b,
                         v_base[TARGET_NODE]["state_1"] * 100, h_base[TARGET_NODE]["state_1"] * 100)
            self.measure("DBN", "S6A", mode, "conditioned_low_pct", a,
                         v_after[TARGET_NODE]["state_1"] * 100, h_after[TARGET_NODE]["state_1"] * 100)
        with (self.reference / "DBN/posterior/edge_parameters_frontend.csv").open(
            encoding="utf-8-sig", newline=""
        ) as stream:
            ref_edges = {(row["realization"], row["from"], row["to"]): row
                         for row in csv.DictReader(stream)}
        models = {}
        for run, link, slope, interval, counts in self.tables["S6B"]:
            model_zip = self.runs / f"DBN/realization-{run}-0/model-details.zip"
            if not model_zip.exists():
                continue
            if run not in models:
                models[run] = zip_json(model_zip, "frozen_model.json")
                path = self.runs / f"tPSEM/realization-{run}-0/model-details.zip"
                if path.exists():
                    notes.setdefault("case_exports_vs_fresh_tPSEM", {})[run] = compare_case_exports(model_zip, path, "")
                path = self.reference / "DBN/posterior/research-results.zip"
                notes.setdefault("case_exports_vs_hosted_DBN", {})[run] = compare_case_exports(model_zip, path)
            parent, target = link.split(" → ")
            edge = next(row for row in models[run]["regression"]["edges"]
                        if row["from"] == parent and row["to"] == target)
            ref = ref_edges[(run, parent, target)]
            limits = interval.strip("[]").split(", ")
            posterior = edge["posterior"]
            condition, weight = counts.split(" / ")
            for metric, printed, observed, key in (
                ("slope", slope, edge["coefficients"]["parent"], "slope"),
                ("lower", limits[0], posterior["coefficient_intervals"]["parent"][0], "posterior_slope_2.5pct"),
                ("upper", limits[1], posterior["coefficient_intervals"]["parent"][1], "posterior_slope_97.5pct"),
                ("conditions", condition, posterior["independent_conditions"], "independent_conditions"),
                ("weight", weight, posterior["effective_weight_sum"], "effective_weight_sum"),
            ):
                self.measure("DBN", "S6B", f"{run} {link}", metric, printed, observed, float(ref[key]))
        self.notes["DBN"] = notes

    def ml(self):
        for method, directory, hosted_directory in (
            ("Gaussian process", "Gaussian_process", "GP"),
            ("SciPy fallback", "Symbolic_fallback", "SciPy_candidate"),
        ):
            fresh = read_json(self.runs / f"ML/{directory}/results.json")
            hosted = zip_json(self.reference / f"ML/{hosted_directory}/research-results.zip")
            for label, link, n, r2, rmse, equation in self.tables["S7"]:
                if label != method:
                    continue
                key = "→".join(KE_LABELS[ke] for ke in link.split("→"))
                values = fresh["fit"]["per_edge_metrics"][key]
                ref = hosted["fit"]["per_edge_metrics"][key]
                for metric, printed in (("n", n), ("r2", r2), ("rmse", rmse)):
                    self.measure("ML", "S7", f"{method} {link}", metric, printed, values[metric], ref[metric])
                if method == "SciPy fallback":
                    fit_formula = fresh["fit"]["per_edge_formulas"][key]
                    ref_formula = hosted["fit"]["per_edge_formulas"][key]
                    self.measure("ML", "S7", f"{method} {link}", "formula_terms", equation_terms(equation),
                                 equation_terms(fit_formula.replace(KE_LABELS[link.split("→")[0]], "x")), equation_terms(ref_formula))
            self.notes.setdefault("ML", {})[method] = {"provenance": fresh["provenance"],
                                                       "input_provenance": fresh["input_provenance"],
                                                       "control_holdout_r2": fresh["diagnostics"]["cv_r2_per_dose"]["0"],
                                                       "curves": {}}
            curves = self.notes["ML"][method]["curves"]
            for key, curve in fresh["fit"]["predictions_map"].items():
                ref = hosted["fit"]["predictions_map"][key]
                assert len(curve) == len(ref)
                curves[key] = {
                    "grid_points": len(curve),
                    "x_grid_identical": all(a["x"] == b["x"] for a, b in zip(curve, ref)),
                    "scatter_identical": fresh["fit"]["scatter_data"][key] == hosted["fit"]["scatter_data"][key],
                    "max_abs_difference": {
                        field: max(abs(a[field] - b[field]) for a, b in zip(curve, ref))
                        for field in ("y_mean", "y_lower", "y_upper")
                    },
                }

    def ode(self):
        fresh = read_json(self.runs / "ODE/synthetic/results.json")
        reference = self.reference / "ODE/synthetic/research-results.zip"
        hosted = zip_json(reference)
        fit = fresh["fitting_result"]
        parameters = {row["parameter_id"]: row for row in fit["fitted_parameters"]}
        host_parameters = {row["parameter_id"]: row for row in hosted["fitting_result"]["fitted_parameters"]}
        for label, estimate, lower, upper, unit in self.tables["S8"]:
            key = "param_" + label.split(" / ")[0]
            for metric, printed in (("estimate", estimate), ("lower", lower), ("upper", upper)):
                self.measure("ODE", "S8", label, metric, printed,
                             parameters[key][metric], host_parameters[key][metric])
        for metric, paper in (("RMSE", "2.0177"), ("MAE", "1.6405"), ("bias", "0.1955"), ("n_observations", "36")):
            self.measure("ODE", "Fig 11/S1", "synthetic 1000 nM", metric, paper,
                         single(fit["fit_metrics"])[metric], single(hosted["fitting_result"]["fit_metrics"])[metric])
        def values(rows):
            return {(row["state_id"], row["time"]): row["predicted_value"] for row in rows}
        current_predictions = values(fit["fitted_predictions"])
        hosted_predictions = values(hosted["fitting_result"]["fitted_predictions"])
        assert current_predictions.keys() == hosted_predictions.keys()
        self.notes["ODE"] = {"synthetic_status": single(fit["status"]),
            "max_abs_synthetic_prediction_delta": max(abs(v - hosted_predictions[k]) for k, v in current_predictions.items()),
            "synthetic_boundary_hits": [name for name, row in parameters.items() if row["boundary_hit"]],
            "effective_start": {name: row["initial_value"] for name, row in parameters.items()}}
        step2 = self.runs / "ODE/synthetic/Step2/preprocessed.csv"
        if step2.exists():
            def fit_values(rows):
                return {(row["ke_id"], float(row["fit_dose"]), float(row["fit_time"])):
                        float(row["fit_response"]) for row in rows}
            with step2.open(encoding="utf-8-sig", newline="") as stream:
                means = fit_values(csv.DictReader(stream))
            old_means = fit_values(zip_csv(reference, "data/preprocessed.csv"))
            assert means.keys() == old_means.keys()
            self.notes["ODE"]["synthetic_step2"] = {
                "condition_means": len(means),
                "identical_condition_means": sum(value == old_means[key] for key, value in means.items()),
                "selected_1000_nm_means": sum(dose == 1000 for _, dose, _ in means),
            }
        measured = self.runs / "ODE/measured/results.json"
        if measured.exists():
            result = read_json(measured)
            self.notes["ODE"]["measured"] = {
                "status": single(result.get("status")), "input_provenance": result.get("input_provenance"),
                "fit_metrics": single((result.get("fitting_result") or {}).get("fit_metrics")),
                "fitted_parameters": (result.get("fitting_result") or {}).get("fitted_parameters"),
                "errors": result.get("errors"),
            }

    def write(self):
        self.out.mkdir(parents=True, exist_ok=True)
        with (self.out / "v16_reproducibility.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
            writer.writeheader()
            writer.writerows(self.rows)
        (self.out / "module_comparison_details.json").write_text(
            json.dumps(self.notes, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    for arg in ("si", "hosted-reference", "runs", "out"):
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--linux-pod", type=Path)
    parser.add_argument("--measured-input", type=Path)
    args = parser.parse_args()
    audit = Audit(args)
    for operation in (audit.check_input, audit.pod, audit.linux_pod, audit.qker, audit.psem,
                      audit.dbn, audit.ml, audit.ode):
        try:
            operation()
        except FileNotFoundError as exc:
            print(operation.__name__, "pending:", exc.filename)
    if audit.rows:
        audit.write()
    groups = sorted({(row["module"], row["reference"]) for row in audit.rows})
    for module, table in groups:
        rows = [row for row in audit.rows if row["module"] == module and row["reference"] == table]
        print(module, table, sum(row["verdict"] == "matches_display" for row in rows), "/", len(rows),
              "match displayed precision")
        for row in rows:
            if row["verdict"] == "differs_from_display":
                print(" ", row["item"], row["metric"], row["printed_v16"], row["fresh_local"])


if __name__ == "__main__":
    main()
