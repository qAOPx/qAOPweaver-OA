"""Demo fits adapted from qAOPweaver 0.29.0, revision 0143dd56."""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler


def gaussian_process(x, y, label):
    """The GP fitting and prediction path from app/aiml/gp.py."""
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x.reshape(-1, 1))
    kernel = (
        ConstantKernel(1.0, (1e-3, 1e3))
        * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 10.0))
        + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-5, 1.0))
    )
    model = GaussianProcessRegressor(
        kernel=kernel, n_restarts_optimizer=5, normalize_y=True, random_state=42,
    )
    model.fit(x_scaled, y)
    grid = np.linspace(float(x.min()), float(x.max()), 100)
    mean, spread = model.predict(scaler.transform(grid.reshape(-1, 1)), return_std=True)
    fitted, _ = model.predict(x_scaled, return_std=True)
    return {
        "formula": None,
        "params": {
            "kernel": str(model.kernel_),
            "log_marginal_likelihood": float(model.log_marginal_likelihood_value_),
        },
        "predictions": prediction_rows(grid, mean, 1.96 * spread),
        "metrics": fit_metrics(y, fitted),
    }


def candidate_formulae():
    """The ordered SciPy candidates and starting values from app/aiml/symbolic.py."""
    def linear(x, a, b):
        return a * x + b

    def power(x, a, b, c):
        return a * np.power(np.abs(x) + 1e-9, b) + c

    def exponential(x, a, b, c):
        return a * np.exp(b * x) + c

    def hill(x, emax, ec50, n):
        xn = np.power(np.abs(x) + 1e-9, np.abs(n) + 1e-9)
        ec50n = np.power(np.abs(ec50) + 1e-9, np.abs(n) + 1e-9)
        return emax * xn / (ec50n + xn)

    def logarithmic(x, a, b):
        return a * np.log(np.abs(x) + 1e-9) + b

    return [
        (linear, lambda a, b, name: f"{a:.4g} * {name} + {b:.4g}", [1.0, 0.0]),
        (power, lambda a, b, c, name: f"{a:.4g} * {name}^{b:.4g} + {c:.4g}", [1.0, 0.5, 0.0]),
        (
            exponential,
            lambda a, b, c, name: f"{a:.4g} * exp({b:.4g} * {name}) + {c:.4g}",
            [1.0, 0.05, 0.0],
        ),
        (
            hill,
            lambda emax, ec50, n, name: (
                f"{emax:.4g} * {name}^{n:.4g} / ({ec50:.4g}^{n:.4g} + {name}^{n:.4g})"
            ),
            [95.0, 8.0, 1.5],
        ),
        (logarithmic, lambda a, b, name: f"{a:.4g} * log({name}) + {b:.4g}", [20.0, 0.0]),
    ]


def scipy_fallback(x, y, label):
    """Select and evaluate the fixed formula families used for this case."""
    best_r2 = -np.inf
    best_pred = np.zeros_like(y)
    best_formula = f"0 * {label}"
    best_params = {}
    for function, formula, initial in candidate_formulae():
        try:
            parameters, _ = curve_fit(function, x, y, p0=initial, maxfev=10_000)
            fitted = function(x, *parameters)
            r2 = float(r2_score(y, fitted))
        except (ValueError, RuntimeError, FloatingPointError, OverflowError):
            continue
        if r2 > best_r2:
            best_r2 = r2
            best_pred = fitted
            best_formula = formula(*parameters, label)
            best_params = {f"p{i}": float(value) for i, value in enumerate(parameters)}

    grid = np.linspace(float(x.min()), float(x.max()), 100)
    order = np.argsort(x)
    mean = np.interp(grid, x[order], best_pred[order])
    half_width = np.full(grid.shape, 1.96 * np.std(y - best_pred))
    return {
        "formula": best_formula,
        "params": best_params,
        "predictions": prediction_rows(grid, mean, half_width),
        "metrics": fit_metrics(y, best_pred),
    }


def prediction_rows(grid, mean, half_width):
    return [
        {"x": float(x), "y_mean": float(y), "y_lower": float(lo), "y_upper": float(hi)}
        for x, y, lo, hi in zip(grid, mean, mean - half_width, mean + half_width)
    ]


def fit_metrics(observed, fitted):
    return {
        "r2": float(r2_score(observed, fitted)),
        "rmse": float(np.sqrt(np.mean((observed - fitted) ** 2))),
        "n": len(observed),
    }


def anchor_holdout(matrix, fit, source, target):
    """MIE-to-AO held-concentration-out diagnostics used in the saved results."""
    per_dose = {}
    covered = 0
    observations = 0
    for dose in sorted(matrix.index.get_level_values("fit_dose").unique()):
        held_out = matrix.index.get_level_values("fit_dose") == dose
        train, test = matrix.loc[~held_out], matrix.loc[held_out]
        curve = fit(train[source].to_numpy(), train[target].to_numpy(), source)["predictions"]
        test_x, test_y = test[source].to_numpy(), test[target].to_numpy()
        grid = [row["x"] for row in curve]
        mean = np.interp(test_x, grid, [row["y_mean"] for row in curve])
        lower = np.interp(test_x, grid, [row["y_lower"] for row in curve])
        upper = np.interp(test_x, grid, [row["y_upper"] for row in curve])
        residual = float(np.sum((test_y - mean) ** 2))
        spread = float(np.sum((test_y - test_y.mean()) ** 2))
        r2 = 1.0 - residual / spread if spread else float(residual == 0)
        per_dose[f"{dose:g}"] = round(r2, 6)
        covered += int(np.sum((test_y >= lower) & (test_y <= upper)))
        observations += len(test_y)
    values = list(per_dose.values())
    return {
        "cv_r2_per_dose": per_dose,
        "cv_r2_mean": round(float(np.mean(values)), 6),
        "cv_r2_sd": round(float(np.std(values)), 6),
        "ci_coverage": round(covered / observations, 4),
        "problematic_doses": [dose for dose, r2 in per_dose.items() if r2 < 0.5],
    }
