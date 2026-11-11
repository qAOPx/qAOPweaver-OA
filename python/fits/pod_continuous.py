"""
Continuous dose-response analysis using EPA pybmds 25.x.
Follows EPA BMD Technical Guidance (EPA/100/R-12/001).

pybmds 25.x API notes (verified against installed package):
  - session.models[i].name()          → str  (method, not property)
  - result.bmd / .bmdl / .bmdu        → float | None
  - result.fit.aic                     → float
  - result.tests.names / .p_values    → list  (Test 4 = GOF p-value)
  - result.parameters.names / .values → arrays
  - result.plotting.dr_x / .dr_y      → arrays (predicted curve)
  - result.has_completed              → bool
  - session.recommender.results.recommended_model_index → int | None
  - ContinuousRiskType.StandardDeviation
"""
from __future__ import annotations

import warnings as _warnings

import numpy as np
import pandas as pd
import pybmds
from scipy.stats import dunnett as _scipy_dunnett

_warnings.filterwarnings("ignore", category=RuntimeWarning)


class ContinuousPODAnalyzer:
    """
    Run pybmds continuous dose-response analysis for one KE column.

    Model suite: Hill, Power, Exponential M3/M5, Linear, Polynomial-2.
    Benchmark response is one standard deviation.
    """

    def analyze(
        self,
        doses: list[float],
        responses: list[float],
        ns: list[int],
        stdevs: list[float],
        ecx_values: list[int] | None = None,
        raw_df: pd.DataFrame | None = None,
        dose_col: str | None = None,
        ke_col: str | None = None,
    ) -> dict:
        if ecx_values is None:
            ecx_values = [10, 20, 50]

        dataset = pybmds.ContinuousDataset(
            doses=doses,
            means=responses,
            ns=ns,
            stdevs=stdevs,
        )

        session = pybmds.Session(dataset=dataset)
        session.add_default_models()

        for model in session.models:
            model.settings.bmr_type = pybmds.ContinuousRiskType.StandardDeviation
            model.settings.bmr = 1.0

        session.execute()
        session.recommend()

        all_models = self._extract_all_models(session)
        best_fit = self._get_best_fit(session, all_models)
        ma_result = self._model_average(all_models)

        noael, loael = self._dunnett_noael_loael(doses, responses, ns, stdevs)

        # Build ECx from the recommended model when available.
        best_fit_model_idx = session.recommender.results.recommended_model_index
        ecx = self._compute_ecx(
            ecx_values,
            session.models[best_fit_model_idx] if best_fit_model_idx is not None else None,
            responses[0],
            max(doses),
        )

        # The recommender can decline to select a model while the AIC fallback
        # still supplies the reported best fit. Plot that reported fit instead
        # of leaving its curve blank; ECx estimates retain their prior policy.
        curve_model_idx = next(
            (idx for idx, model in enumerate(all_models) if model is best_fit), None
        )
        curve = self._build_curve(session, curve_model_idx, doses)
        observed_groups = self._build_observed_groups(doses, responses, stdevs, ns)

        return {
            "best_fit_model": best_fit,
            "model_averaged": ma_result,
            "noael": noael,
            "loael": loael,
            "ecx": ecx,
            "all_models": all_models,
            "dose_response_curve": curve,
            "observed_groups": observed_groups,
            "flags": self._collect_flags(best_fit, all_models),
            "n_doses": len(doses),
            "n_replicates": int(round(sum(ns) / len(ns))),
        }

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _extract_all_models(self, session) -> list[dict]:
        results = []
        for model in session.models:
            r = model.results
            name = model.name()
            if r is None or not r.has_completed:
                results.append({
                    "name": name,
                    "bmd": None, "bmdl": None, "bmdu": None,
                    "aic": None, "goodness_of_fit_p": None,
                    "converged": False, "flagged": False,
                    "parameters": {},
                })
                continue

            gof_p = None
            if hasattr(r, "tests") and r.tests is not None:
                pvs = r.tests.p_values
                gof_p = float(pvs[3]) if len(pvs) > 3 else None

            params = {}
            if r.parameters is not None:
                params = dict(zip(r.parameters.names, r.parameters.values.tolist()))

            bmd = r.bmd if (r.bmd is not None and r.bmd > 0) else None
            results.append({
                "name": name,
                "bmd": bmd,
                "bmdl": r.bmdl if (r.bmdl and r.bmdl > 0) else None,
                "bmdu": r.bmdu if (r.bmdu and r.bmdu > 0) else None,
                "aic": float(r.fit.aic) if r.fit is not None else None,
                "goodness_of_fit_p": gof_p,
                "converged": bmd is not None,
                "flagged": self._is_hill_flagged(name, params),
                "parameters": params,
            })
        return results

    def _is_hill_flagged(self, name: str, params: dict) -> bool:
        """Flag Hill if n (power parameter) is at its lower bound (≤1), indicating poor fit."""
        if "Hill" in name:
            n = params.get("n")
            if n is not None and n <= 1.001:
                return True
        return False

    def _get_best_fit(self, session, all_models: list[dict]) -> dict | None:
        idx = session.recommender.results.recommended_model_index
        if idx is None or idx >= len(all_models):
            # Fallback: lowest AIC among converged, unflagged models
            candidates = [m for m in all_models if m["converged"] and not m["flagged"] and m["aic"] is not None]
            if not candidates:
                return None
            return min(candidates, key=lambda m: m["aic"])
        return all_models[idx]

    def _model_average(self, all_models: list[dict]) -> dict | None:
        """AIC-weighted model average over all converged, unflagged models."""
        converged = [
            m for m in all_models
            if m["converged"] and m["aic"] is not None and not m["flagged"]
            and m["bmd"] is not None
        ]
        if not converged:
            return None

        aics = np.array([m["aic"] for m in converged])
        delta = aics - aics.min()
        weights = np.exp(-0.5 * delta)
        weights /= weights.sum()

        def _wa(key: str) -> float | None:
            vals = [m[key] for m in converged]
            if any(v is None for v in vals):
                return None
            return float(np.dot(weights, vals))

        return {
            "bmd": _wa("bmd"),
            "bmdl": _wa("bmdl"),
            "bmdu": _wa("bmdu"),
            "model_weights": {
                m["name"]: float(w) for m, w in zip(converged, weights)
            },
        }

    def _compute_ecx(
        self,
        x_pcts: list[int],
        best_model,
        control_response: float,
        max_dose: float,
    ) -> dict:
        ecx: dict = {f"ec{x}": None for x in x_pcts}
        ecx["custom"] = {}

        if best_model is None:
            return ecx
        r = best_model.results
        if r is None or not r.has_completed:
            return ecx

        dr_x = np.array(r.plotting.dr_x)
        dr_y = np.array(r.plotting.dr_y)

        if len(dr_x) < 2 or len(dr_y) < 2:
            return ecx

        increasing = dr_y[-1] > dr_y[0]

        for x in x_pcts:
            target = control_response * (1 - x / 100) if increasing else control_response * (1 + x / 100)
            try:
                # Find crossing point via interpolation
                diff = dr_y - target
                sign_changes = np.where(np.diff(np.sign(diff)))[0]
                if len(sign_changes) == 0:
                    ecx[f"ec{x}"] = None
                    continue
                i = sign_changes[0]
                # Linear interpolation between i and i+1
                x0, x1 = dr_x[i], dr_x[i + 1]
                y0, y1 = dr_y[i] - target, dr_y[i + 1] - target
                dose_ec = float(x0 - y0 * (x1 - x0) / (y1 - y0))
                ecx[f"ec{x}"] = dose_ec if 0 < dose_ec <= max_dose * 10 else None
            except Exception:
                ecx[f"ec{x}"] = None

        return ecx

    def _dunnett_noael_loael(
        self,
        doses: list[float],
        responses: list[float],
        ns: list[int],
        stdevs: list[float],
    ) -> tuple[float | None, float | None]:
        """
        Dunnett's test (scipy >= 1.9): compare each treatment group vs control.
        NOAEL = highest non-significant dose; LOAEL = lowest significant dose.
        """
        if len(doses) < 2:
            return None, None

        control_mean = responses[0]
        control_n = ns[0]
        control_sd = stdevs[0]

        # Reconstruct per-replicate arrays (approximate from group stats)
        control_samples = np.random.default_rng(42).normal(
            loc=control_mean, scale=max(control_sd, 1e-8), size=control_n
        )

        noael: float | None = None
        loael: float | None = None

        for i in range(1, len(doses)):
            treatment = np.random.default_rng(i).normal(
                loc=responses[i], scale=max(stdevs[i], 1e-8), size=ns[i]
            )
            try:
                res = _scipy_dunnett(treatment, control=control_samples)
                significant = bool(res.pvalue[0] < 0.05)
            except Exception:
                # Fallback: pooled-t approximation
                se = np.sqrt(control_sd ** 2 / control_n + stdevs[i] ** 2 / ns[i])
                t_stat = abs(responses[i] - control_mean) / max(se, 1e-10)
                significant = t_stat > 2.5

            if not significant:
                noael = doses[i]
            elif loael is None:
                loael = doses[i]

        return noael, loael

    def _build_curve(self, session, best_idx: int | None, doses: list[float]) -> list[dict]:
        if best_idx is None:
            return []
        m = session.models[best_idx]
        r = m.results
        if r is None or not r.has_completed:
            return []
        dr_x = r.plotting.dr_x
        dr_y = r.plotting.dr_y
        # NOTE on ci_lower/ci_upper: pybmds only exposes bmd_y/bmdl_y/bmdu_y —
        # single scalar y-values at the BMD/BMDL/BMDU dose points, not a
        # per-dose confidence band matching the full dr_x/dr_y curve array.
        # Hand-rolling a per-point CI would be new (unreviewed) statistics,
        # so this stays None until pybmds exposes a real per-dose CI.
        return [
            {"dose": float(x), "predicted": float(y), "ci_lower": None, "ci_upper": None}
            for x, y in zip(dr_x, dr_y)
        ]

    def _build_observed_groups(
        self,
        doses: list[float],
        responses: list[float],
        stdevs: list[float],
        ns: list[int],
    ) -> list[dict]:
        """Echo the observed per-dose-group summary stats (the plot's data points)."""
        return [
            {"dose": float(d), "mean": float(m), "sd": float(s), "n": int(n)}
            for d, m, s, n in zip(doses, responses, stdevs, ns)
        ]

    def _collect_flags(self, best_fit: dict | None, all_models: list[dict]) -> list[str]:
        flags: list[str] = []
        converged = [m for m in all_models if m["converged"]]
        if not converged:
            flags.append("NO_MODELS_CONVERGED")
        elif len(converged) == 1:
            flags.append("ONLY_ONE_MODEL_CONVERGED")
        if best_fit and best_fit.get("flagged"):
            flags.append("BEST_FIT_MODEL_FLAGGED")
        return flags
