"""Explicit provenance and fixed-reference swirl features; no hidden data joins."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .physics import effectiveness

META = ["condition_id", "configuration", "flow_coefficient", "beta_a", "sigma_beta_a",
        "beta_c_ref", "sigma_beta_c_ref", "phi_ref", "r_epsilon_over_b", "r_swirl_over_b",
        "density_ratio", "rig_id", "seal_id", "normalization", "source_kind", "source"]
NUMERIC = ["phi0", "epsilon", "sigma_epsilon", "sigma_phi0", "flow_coefficient", "beta_a",
           "sigma_beta_a", "beta_c_ref", "sigma_beta_c_ref", "phi_ref", "r_epsilon_over_b",
           "r_swirl_over_b", "density_ratio"]
COLUMNS = ["curve_id", "condition_id", "configuration", "flow_coefficient", "phi0", "epsilon",
           "sigma_epsilon", "sigma_phi0", "beta_a", "sigma_beta_a", "beta_c_ref", "sigma_beta_c_ref",
           "phi_ref", "r_epsilon_over_b", "r_swirl_over_b", "density_ratio", "rig_id", "seal_id",
           "normalization", "source_kind", "source"]


def validate(frame: pd.DataFrame) -> pd.DataFrame:
    """One dataset = one physical setup and one measurement-definition protocol.

    Each curve is a purge sweep at one condition/configuration. beta_c_ref is
    measured at a predeclared common reference purge, NOT the instantaneous
    beta_c at every point of the effectiveness sweep. Never silently pool radii.
    """
    d = frame.copy()
    missing = sorted(set(COLUMNS)-set(d))
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if d.empty:
        raise ValueError("Empty dataset; fill the template with authorized measurements.")
    for c in NUMERIC:
        d[c] = pd.to_numeric(d[c], errors='raise')
    if not np.isfinite(d[NUMERIC].to_numpy(float)).all():
        raise ValueError("Numeric fields require finite values; missing uncertainty is not zero.")
    for c in set(COLUMNS)-set(NUMERIC):
        if d[c].isna().any() or d[c].astype(str).str.strip().eq('').any():
            raise ValueError(f"Missing provenance/identifier in {c}.")
        d[c] = d[c].astype(str)
    if (d.sigma_epsilon <= 0).any() or (d[["sigma_phi0","sigma_beta_a","sigma_beta_c_ref"]] < 0).any().any():
        raise ValueError("sigma_epsilon must be >0; other standard uncertainties >=0.")
    if (d.phi0 < 0).any() or (d.phi_ref <= 0).any() or (d.flow_coefficient <= 0).any():
        raise ValueError("Require Phi0>=0, phi_ref>0 and flow_coefficient>0.")
    if not np.allclose(d.density_ratio, 1, rtol=0, atol=1e-8):
        raise ValueError("This release implements DR=1 only. Density scaling is not implemented.")
    if not d.configuration.isin(["baseline","rotor_swirler"]).all():
        raise ValueError("Supported configurations: baseline and rotor_swirler.")
    if not d.source_kind.isin(["synthetic","published_digitized","lab_measured"]).all():
        raise ValueError("Declare source_kind as synthetic, published_digitized or lab_measured.")
    for c in ["rig_id","seal_id","normalization","source_kind","r_epsilon_over_b","r_swirl_over_b","phi_ref"]:
        if d[c].nunique() != 1:
            raise ValueError(f"Do not pool different {c}: analyze separately.")
    if (d[["r_epsilon_over_b","r_swirl_over_b"]] <= 0).any().any():
        raise ValueError("Measurement radii must be positive.")
    for key, g in d.groupby('curve_id', sort=False):
        for c in META:
            if g[c].nunique() != 1:
                raise ValueError(f"Curve {key}: {c} varies within a purge sweep.")
        if g.phi0.nunique() < 3:
            raise ValueError(f"Curve {key}: at least three distinct purge levels required.")
    for key, g in d.groupby('condition_id', sort=False):
        # beta_a uses a shared condition-level measurement error in Monte Carlo.
        for c in ["flow_coefficient", "beta_a", "sigma_beta_a"]:
            if g[c].nunique() != 1:
                raise ValueError(f"Condition {key}: inconsistent {c}; condition-level pairing is invalid.")
    if d.groupby('flow_coefficient').condition_id.nunique().max() > 1:
        raise ValueError("Use the same condition_id for every curve at the same flow coefficient.")
    return d.reset_index(drop=True)


def synthetic_data(scenario: str = "swirl", seed: int = 20261005) -> pd.DataFrame:
    """NOT Song/Vella data. Separate null and misspecification negative controls.

    The swirl scenario is deliberately generated from the candidate correction;
    recovering it demonstrates software operation, not scientific evidence.
    """
    if scenario not in ("swirl", "null", "purge_dependence"):
        raise ValueError("Unknown synthetic scenario.")
    rng = np.random.default_rng(seed)
    rows = []
    for j, cf in enumerate([0.30,0.35,0.40,0.45,0.50]):
        ba = 0.18 + 2.0*cf
        bc_base = 0.44 + 0.045*(ba-1)
        for cfg in ["baseline","rotor_swirler"]:
            delta = 0 if cfg == "baseline" else 0.085+0.07*(cf-0.30)
            bc = bc_base + delta
            a = (0.015+0.115*ba) * np.exp((-2.5 if scenario != "null" else 0)*delta)
            # Same imposed purge grid across both configurations; no truth-derived grid leak.
            for p in np.linspace(0.012,0.205,17):
                actual_a = a*(1+1.8*(p-0.075)) if scenario == "purge_dependence" else a
                eps = float(effectiveness(p, actual_a)) + rng.normal(0,0.006)
                rows.append(dict(curve_id=f"{j}_{cfg}",condition_id=f"CF{cf:.2f}",configuration=cfg,
                    flow_coefficient=cf,phi0=p,epsilon=eps,sigma_epsilon=0.006,sigma_phi0=0.00035,
                    beta_a=ba,sigma_beta_a=0.002,beta_c_ref=bc,sigma_beta_c_ref=0.002,
                    phi_ref=0.075,r_epsilon_over_b=0.941,r_swirl_over_b=0.923,density_ratio=1.,
                    rig_id="SYNTHETIC_NOT_SNU",seal_id="SYNTHETIC",normalization="synthetic_consistent_phi",
                    source_kind="synthetic",source=f"synthetic:{scenario};seed={seed}"))
    return pd.DataFrame(rows)[COLUMNS]


def perturb(d: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Measurement-error propagation, NOT an iid bootstrap or a posterior.

    epsilon and Phi0: independent point errors; beta_a: shared per condition;
    beta_c_ref: shared per curve. No clipping of noisy epsilon values. Positive
    Phi0 errors use moment-matched lognormal draws; Phi0=0 requires sigma=0.
    """
    q = d.copy()
    q.epsilon = d.epsilon + rng.normal(size=len(d))*d.sigma_epsilon
    p, s = d.phi0.to_numpy(), d.sigma_phi0.to_numpy()
    if np.any((p == 0) & (s > 0)):
        raise ValueError("For zero Phi0, supply sigma_phi0=0 or a dedicated censoring model.")
    mask = p > 0
    v = np.log1p((s[mask]/p[mask])**2)
    q.loc[mask,'phi0'] = rng.lognormal(np.log(p[mask])-v/2, np.sqrt(v))
    for _, idx in d.groupby('condition_id').groups.items():
        q.loc[idx,'beta_a'] += rng.normal()*d.loc[idx,'sigma_beta_a'].iloc[0]
    for _, idx in d.groupby('curve_id').groups.items():
        q.loc[idx,'beta_c_ref'] += rng.normal()*d.loc[idx,'sigma_beta_c_ref'].iloc[0]
    return q
