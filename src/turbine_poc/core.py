"""Numerical core. All public-data targets retain the paper's units and signs."""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np

INPUTS = ("x3", "x4", "x5", "x6")
TARGETS = ("delta_eta_pp", "delta_torque_nm")
PAIRS = tuple(combinations(range(4), 2))
TERMS = ("intercept", *INPUTS, *(f"{x}^2" for x in INPUTS),
         *(f"{INPUTS[i]}*{INPUTS[j]}" for i, j in PAIRS))
DATA_DIR = Path(__file__).parent / "data"
CONTEXT = "K1_2025_TABLE2_ORIGINAL"
REFERENCE = "K1_NO_LEAKAGE_REFERENCE"


def finite_array(values, columns: int) -> np.ndarray:
    a = np.asarray(values, dtype=np.float64)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.ndim != 2 or a.shape[1] != columns or len(a) == 0 or not np.isfinite(a).all():
        raise ValueError(f"Expected a finite, non-empty (n, {columns}) array.")
    return a


def features(x) -> np.ndarray:
    x = finite_array(x, 4)
    return np.column_stack((np.ones(len(x)), x, x*x,
                            *(x[:, i]*x[:, j] for i, j in PAIRS)))


def in_support(x) -> np.ndarray:
    """Exact convex hull of the K1 four-factor Box-Behnken design.

    Vertices are all permutations of (+/-1, +/-1, 0, 0).
    Being inside this hull is not a prediction-accuracy certificate.
    """
    x = finite_array(x, 4)
    return (np.max(np.abs(x), axis=1) <= 1+1e-10) & (np.abs(x).sum(axis=1) <= 2+1e-10)


def load_k1():
    path = DATA_DIR / "k1_table2.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    ids = [r["case_id"] for r in rows]
    x = finite_array([[float(r[k]) for k in INPUTS] for r in rows], 4)
    y = finite_array([[float(r[k]) for k in TARGETS] for r in rows], 2)
    if len(rows) != 25 or ids != [f"{i:02d}" for i in range(25)]:
        raise ValueError("K1 table must contain exactly cases 00-24 in order.")
    expected = {(0., 0., 0., 0.)}
    for i, j in PAIRS:
        for a in (-1., 1.):
            for b in (-1., 1.):
                point = [0.]*4
                point[i], point[j] = a, b
                expected.add(tuple(point))
    if set(map(tuple, x)) != expected:
        raise ValueError("K1 design coordinates do not match the published Box-Behnken design.")
    return ids, x, y, hashlib.sha256(path.read_bytes()).hexdigest()


def metrics(observed, predicted) -> dict:
    a, b = np.asarray(observed, float), np.asarray(predicted, float)
    if a.shape != b.shape or a.size == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Metric inputs must have identical, finite, non-empty shapes.")
    residual = b-a
    variance = np.sum((a-a.mean())**2)
    return {"rmse": float(np.sqrt(np.mean(residual**2))),
            "mae": float(np.mean(np.abs(residual))),
            "max_abs_error": float(np.max(np.abs(residual))),
            "r2": float(1-np.sum(residual**2)/variance) if variance > 0 else None}


@dataclass
class QuadraticRSM:
    coefficients: np.ndarray
    rank: int
    singular_values: np.ndarray

    @classmethod
    def fit(cls, x, y):
        a, y = features(x), finite_array(y, 2)
        if len(a) != len(y):
            raise ValueError("Inputs and targets have different row counts.")
        coef, _, rank, singular = np.linalg.lstsq(a, y, rcond=None)
        if rank < len(TERMS):
            raise ValueError(f"Quadratic is not identifiable: rank={rank}, required={len(TERMS)}.")
        return cls(coef, int(rank), singular)

    def predict(self, x) -> np.ndarray:
        x = finite_array(x, 4)
        if not in_support(x).all():
            raise ValueError("Prediction outside the K1 design convex hull is disabled.")
        return features(x) @ self.coefficients


def loo_diagnostics(ids, x, y) -> list[dict]:
    """Leave-one-out with rank checks; never silently use an arbitrary pseudoinverse.

    Removing the only center point aliases intercept with the sum of squares.
    That fold is explicitly undefined for this unconstrained quadratic model.
    """
    a = features(x)
    leverage = np.sum(a * np.linalg.pinv(a).T, axis=1)
    out = []
    for i, case in enumerate(ids):
        keep = np.arange(len(x)) != i
        rank = int(np.linalg.matrix_rank(a[keep]))
        row = {"case_id": case, "training_rank": rank, "leverage": float(leverage[i]),
               "status": "ok" if rank == len(TERMS) else "unidentifiable_fold"}
        predicted = QuadraticRSM.fit(x[keep], y[keep]).predict(x[i])[0] if rank == len(TERMS) else None
        for j, target in enumerate(TARGETS):
            row[f"observed_{target}"] = float(y[i, j])
            row[f"predicted_{target}"] = float(predicted[j]) if predicted is not None else None
            row[f"error_{target}"] = float(predicted[j]-y[i, j]) if predicted is not None else None
        out.append(row)
    return out


def endpoint_contrasts(x, y) -> list[dict]:
    """Matched endpoint contrasts. These do NOT measure curvature or per-kg/s sensitivity."""
    rows = []
    for k, name in enumerate(INPUTS):
        differences = []
        for i in np.flatnonzero(x[:, k] == -1):
            other = x[i].copy()
            other[k] = 1
            match = np.flatnonzero(np.all(np.isclose(x, other, atol=1e-12, rtol=0), axis=1))
            if len(match) == 1:
                differences.append(y[match[0]]-y[i])
        d = np.asarray(differences)
        for j, target in enumerate(TARGETS):
            rows.append({"variable": name, "target": target, "matched_pairs": len(d),
                         "mean_plus_minus": float(d[:, j].mean()),
                         "min_plus_minus": float(d[:, j].min()),
                         "max_plus_minus": float(d[:, j].max())})
    return rows


@dataclass
class Baseline:
    source_kind: str
    provenance: str
    flows: dict[str, float | None]

    @classmethod
    def load(cls, path):
        obj = json.loads(Path(path).read_text(encoding="utf-8"))
        kind, source = obj["source_kind"], obj["provenance"]
        if kind not in ("illustrative", "lab_provided") or not isinstance(source, str) or not source.strip():
            raise ValueError("Explicit source_kind and provenance are required.")
        flows = {}
        for i in range(1, 7):
            value = obj["flow_kg_s"].get(str(i))
            if value is not None:
                if isinstance(value, bool) or not np.isfinite(float(value)) or float(value) < 0:
                    raise ValueError("Flows must be finite non-negative kg/s values or null.")
                value = float(value)
            flows[str(i)] = value
        if any(flows[k] is None or flows[k] <= 0 for k in ("3", "4")):
            raise ValueError("Positive design flows for passages 3 and 4 are required; do not guess lab values.")
        if kind == "lab_provided" and source.startswith("REPLACE"):
            raise ValueError("Replace the template provenance with the authorized source.")
        return cls(kind, source, flows)

    @property
    def m34(self):
        return self.flows["3"] + self.flows["4"]

    @property
    def alpha0(self):
        return self.flows["3"] / self.m34

    @property
    def total6(self):
        return sum(self.flows.values()) if all(v is not None for v in self.flows.values()) else None

    def alpha_bounds(self, budget_scale=1.):
        if not np.isfinite(budget_scale) or budget_scale <= 0:
            raise ValueError("budget_scale must be positive and finite.")
        total = self.m34 * budget_scale
        m3, m4 = self.flows["3"], self.flows["4"]
        low = max(0., 0.5*m3/total, 1-1.5*m4/total)
        high = min(1., 1.5*m3/total, 1-0.5*m4/total)
        if high < low-1e-12:
            raise ValueError("No allocations satisfy both passage bounds at this total budget.")
        if high < low:
            high = low
        return float(low), float(high)

    def allocation(self, alpha, budget_scale=1.):
        a = np.atleast_1d(np.asarray(alpha, float))
        lo, hi = self.alpha_bounds(budget_scale)
        if a.ndim != 1 or not np.isfinite(a).all() or np.any(a < lo-1e-12) or np.any(a > hi+1e-12):
            raise ValueError(f"alpha must lie within [{lo:.8g}, {hi:.8g}].")
        total = self.m34*budget_scale
        m3, m4 = a*total, (1-a)*total
        x = np.column_stack((2*(m3/self.flows["3"]-1),
                             2*(m4/self.flows["4"]-1), np.zeros(len(a)), np.zeros(len(a))))
        return m3, m4, x


def best_on_slice(model, baseline, budget_scale=1., target_index=0):
    """Exact maximization of the fitted quadratic on a bounded one-dimensional slice."""
    if target_index not in (0, 1):
        raise ValueError("Unknown target index.")
    lo, hi = baseline.alpha_bounds(budget_scale)
    if hi-lo < 1e-12:
        return lo, float(model.predict(baseline.allocation(lo, budget_scale)[2])[0, target_index])
    mid, half = (lo+hi)/2, (hi-lo)/2
    a = np.array([lo, mid, hi])
    values = model.predict(baseline.allocation(a, budget_scale)[2])[:, target_index]
    c0, c1, c2 = np.polynomial.polynomial.polyfit([-1, 0, 1], values, 2)
    candidates = [lo, hi]
    if c2 < -1e-14:
        t = -c1/(2*c2)
        if -1 <= t <= 1:
            candidates.append(mid+half*t)
    scores = model.predict(baseline.allocation(candidates, budget_scale)[2])[:, target_index]
    j = int(np.argmax(scores))  # delta eta is negative; maximize, do NOT minimize!
    return float(candidates[j]), float(scores[j])
