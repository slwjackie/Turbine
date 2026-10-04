"""Conditional WLS calibration; proposed swirl correction is not a published law."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, least_squares
from .physics import effectiveness
from .data import META


class NotIdentifiable(ValueError):
    """Insufficient information; do not replace an unidentified fit by a number."""


def fit_curve(curve: pd.DataFrame) -> dict:
    """Fit one constant Phi_A to a purge sweep using absolute sigma_epsilon.

    Phi0 is held fixed in this fit; x-error is propagated by the outer Monte Carlo.
    se_phi_a is a local, conditional linearization, not a total model uncertainty.
    No clipping of observations and no artificial residual-variance inflation.
    """
    p, y, s = (curve[c].to_numpy(float) for c in ('phi0','epsilon','sigma_epsilon'))
    if not (np.isfinite([p,y,s]).all() and (p >= 0).all() and (s > 0).all()):
        raise ValueError('Invalid curve values.')
    informative = (p > 0) & (y > 2*s) & (y < 1-2*s)
    if len(np.unique(p[informative])) < 3:
        raise NotIdentifiable('Need >=3 distinct unsaturated, informative purge levels.')
    scale = float(max(p.max(), 1e-6))
    bounds = (np.log(scale*1e-4), np.log(scale*1e3))
    def objective(z):
        return float(np.sum(((effectiveness(p, np.exp(z))-y)/s)**2))
    # A broad direct golden search can get stuck on the low-A saturation plateau.
    # Locate the best log-grid cell first, then refine it (one-dimensional search).
    grid = np.linspace(*bounds, 129)
    costs = np.sum(((effectiveness(p[:,None], np.exp(grid)[None,:])-y[:,None])/s[:,None])**2,axis=0)
    j = int(np.argmin(costs))
    if j in (0,len(grid)-1):
        raise NotIdentifiable('Amplitude optimum lies outside the search grid.')
    fit = minimize_scalar(objective, bounds=(grid[j-1],grid[j+1]), method='bounded',
                          options={'xatol':1e-11})
    if not fit.success or min(fit.x-bounds[0], bounds[1]-fit.x) < 1e-4:
        raise NotIdentifiable('Amplitude fit failed or reached search boundary.')
    a = float(np.exp(fit.x))
    h = 1e-5
    derivative = (effectiveness(p,a*np.exp(h))-effectiveness(p,a*np.exp(-h)))/(2*h*a)
    information = float(np.sum((derivative/s)**2))
    if information <= 1e-12:
        raise NotIdentifiable('Effectiveness curve is insensitive to amplitude.')
    return dict(phi_a=a, se_phi_a=information**-0.5, chi2=objective(fit.x),
                dof=len(p)-1, reduced_chi2=objective(fit.x)/(len(p)-1), n_points=len(p),
                n_informative=int(informative.sum()))


def fit_curves(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, g in data.groupby('curve_id', sort=True):
        try:
            fit = fit_curve(g)
        except NotIdentifiable as exc:
            raise NotIdentifiable(f'{key}: {exc}') from exc
        rows.append(dict(curve_id=key, **g.iloc[0][META].to_dict(), **fit))
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class ScalingFit:
    model: str
    beta_low: float
    beta_high: float
    a_low: float
    a_high: float
    k: float
    cavity_intercept: float
    cavity_slope: float
    rank: int
    jacobian_condition: float
    weighted_cost: float

    def predict(self, beta_a, beta_c_ref):
        ba, bc = np.broadcast_arrays(np.asarray(beta_a,float), np.asarray(beta_c_ref,float))
        if not np.isfinite([ba,bc]).all():
            raise ValueError('Swirl features must be finite.')
        w = (ba-self.beta_low)/(self.beta_high-self.beta_low)
        base = self.a_low+(self.a_high-self.a_low)*w
        if (base <= 0).any():
            raise NotIdentifiable('Linear annulus scaling becomes nonpositive on extrapolation.')
        delta = bc-(self.cavity_intercept+self.cavity_slope*ba)
        result = base*np.exp(self.k*delta)
        if not np.isfinite(result).all() or (result <= 0).any():
            raise NotIdentifiable('Invalid amplitude prediction.')
        return result

    def extrapolation(self, beta_a):
        ba = np.asarray(beta_a,float)
        return (ba < self.beta_low) | (ba > self.beta_high)

    def record(self):
        r = asdict(self)
        r['annulus_slope'] = (self.a_high-self.a_low)/(self.beta_high-self.beta_low)
        r['annulus_intercept'] = self.a_low-r['annulus_slope']*self.beta_low
        return r


def fit_scaling(curves: pd.DataFrame, model: str) -> ScalingFit:
    """Train-only feature centering and amplitude regression.

    M0: A=a+b*beta_a (Vella's linear hypothesis; coefficients re-fit to this rig).
    M1: A=(a+b*beta_a)*exp(k*(beta_c_ref-baseline_beta_c(beta_a))).
    M1 is our one-extra-parameter candidate, not a formula attributed to Vella.
    Both models see the SAME training curves. No test epsilon is read here.
    """
    if model not in ('annulus','swirl'):
        raise ValueError('Model must be annulus or swirl.')
    ba, bc, y, s = (curves[c].to_numpy(float) for c in ('beta_a','beta_c_ref','phi_a','se_phi_a'))
    if not np.isfinite([ba,bc,y,s]).all() or (y <= 0).any() or (s <= 0).any():
        raise ValueError('Invalid amplitude training data.')
    lo, hi = float(ba.min()), float(ba.max())
    if hi-lo < 1e-6 or len(y) < 3:
        raise NotIdentifiable('Annulus scaling requires distinct annulus swirl levels.')
    intercept, slope = 0., 0.
    if model == 'swirl':
        baseline = curves[curves.configuration == 'baseline']
        if len(baseline) < 2 or not (curves.configuration == 'rotor_swirler').any():
            raise NotIdentifiable('Swirl correction needs training baseline AND swirler curves.')
        xb = np.column_stack([np.ones(len(baseline)), baseline.beta_a])
        if np.linalg.matrix_rank(xb) < 2:
            raise NotIdentifiable('Baseline reference-swirl trend is not identifiable.')
        intercept, slope = map(float, np.linalg.lstsq(xb,baseline.beta_c_ref,rcond=None)[0])
    delta = bc-(intercept+slope*ba)
    w = (ba-lo)/(hi-lo)
    design = np.column_stack([np.ones(len(ba)), w] + ([delta] if model=='swirl' else []))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise NotIdentifiable('Annulus and cavity predictors are confounded.')
    def prediction(t):
        base = np.exp(t[0])*(1-w)+np.exp(t[1])*w
        return base*np.exp(t[2]*delta) if model=='swirl' else base
    n = 3 if model=='swirl' else 2
    t0 = [np.log(np.median(y))]*2 + ([0.] if n==3 else [])
    lower = [-18.,-18.] + ([-30.] if n==3 else [])
    upper = [5.,5.] + ([30.] if n==3 else [])
    fit = least_squares(lambda t:(prediction(t)-y)/s, t0, bounds=(lower,upper), max_nfev=2000,
                        xtol=1e-11, ftol=1e-11, gtol=1e-11)
    rank = int(np.linalg.matrix_rank(fit.jac))
    cond = float(np.linalg.cond(fit.jac))
    if not fit.success or rank < n or cond > 1e9 or np.any(fit.active_mask):
        raise NotIdentifiable('Scaling fit failed, is ill-conditioned, or reached a bound.')
    return ScalingFit(model,lo,hi,float(np.exp(fit.x[0])),float(np.exp(fit.x[1])),
                      float(fit.x[2]) if n==3 else 0.,intercept,slope,rank,cond,float(2*fit.cost))
