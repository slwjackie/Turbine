"""Tang et al. (2024), doi:10.1016/j.ijheatfluidflow.2024.109300.

Equations 18, 30-35; equal-density sealing-effectiveness component ONLY.
No prediction of Phi_A from the Navier-Stokes equations is implied.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq


def _pair(phi0, phi_a):
    p, a = np.broadcast_arrays(np.asarray(phi0, float), np.asarray(phi_a, float))
    if not (np.isfinite(p).all() and np.isfinite(a).all()):
        raise ValueError("Phi0 and Phi_A must be finite.")
    if (p < 0).any() or (a <= 0).any():
        raise ValueError("Require Phi0 >= 0 and Phi_A > 0.")
    return p, a


def fluxes(phi0, phi_a):
    """Return nonnegative (ingress, egress); egress-ingress == Phi0."""
    p, a = _pair(phi0, phi_a)
    x = np.minimum(p / a, 1.0)
    xi = np.arcsin(x)
    egress = a / np.pi * np.sqrt(np.maximum(1 - x*x, 0)) + p / np.pi * (xi + np.pi/2)
    ingress = np.maximum(egress - p, 0)
    return ingress, egress


def effectiveness(phi0, phi_a):
    """Eq. 33, evaluated without the Phi_A/Phi0 singularity at Phi0=0."""
    p, a = _pair(phi0, phi_a)
    _, egress = fluxes(p, a)
    return np.clip(p / egress, 0, 1)


def required_flow(phi_a, target: float = 0.95):
    """For constant Phi_A only. Phi_95 != Phi_min; Phi_min=Phi_A."""
    _, a = _pair(0, phi_a)
    if not np.isfinite(target) or not 0 < target <= 1:
        raise ValueError("Require 0 < target <= 1.")
    ratio = 1.0 if target == 1 else brentq(lambda x: float(effectiveness(x, 1))-target, 0, 1, xtol=1e-14)
    return a * ratio


def amplitude_from_point(phi0: float, epsilon: float) -> float:
    """Interior-point inverse, not a unique estimate at epsilon=1 or Phi0=0."""
    if not np.isfinite(phi0) or phi0 <= 0 or not np.isfinite(epsilon) or not 0 < epsilon < 1:
        raise ValueError("A unique point inversion requires Phi0>0 and 0<epsilon<1.")
    return float(phi0 / required_flow(1.0, epsilon))


def mass_to_phi(mass_flow, density: float, omega: float, radius: float, clearance: float):
    """Eq. 18. kg/s, kg/m^3, rad/s, m, m; explicitly selected reference radius/gap."""
    if not all(np.isfinite(v) and v > 0 for v in (density, omega, radius, clearance)):
        raise ValueError("Density, omega, radius and clearance must be positive SI values.")
    m = np.asarray(mass_flow, float)
    if not np.isfinite(m).all() or (m < 0).any():
        raise ValueError("Mass flow must be finite and nonnegative.")
    return m / (2*np.pi*density*omega*radius**2*clearance)


def renormalize_phi(phi, *, old_density, old_omega, old_radius, old_gap,
                    new_density, new_omega, new_radius, new_gap):
    """Preserve dimensional mass flow, rather than mixing two Phi conventions."""
    old_factor = 1 / mass_to_phi(1, old_density, old_omega, old_radius, old_gap)
    return mass_to_phi(np.asarray(phi)*old_factor, new_density, new_omega, new_radius, new_gap)


def tracer_effectiveness(c_sample: float, c_annulus: float, c_purge: float, covariance=None):
    """(cs-ca)/(c0-ca), with optional 3x3 covariance in order (cs,ca,c0).

    Output is NOT clipped: noisy measured concentration ratios may exceed [0,1].
    Uncertainty is a first-order delta-method standard uncertainty, not a 95% bound.
    """
    cs, ca, c0 = map(float, (c_sample, c_annulus, c_purge))
    if not np.isfinite([cs, ca, c0]).all() or c0 <= ca:
        raise ValueError("Require finite concentrations in common units and c_purge>c_annulus.")
    d = c0-ca
    eps = (cs-ca)/d
    if covariance is None:
        return eps, None
    cov = np.asarray(covariance, float)
    if cov.shape != (3,3) or not np.isfinite(cov).all() or not np.allclose(cov, cov.T):
        raise ValueError("Covariance must be finite, symmetric, and 3x3.")
    if np.linalg.eigvalsh(cov).min() < -1e-14:
        raise ValueError("Covariance must be positive semidefinite.")
    grad = np.array([1/d, (eps-1)/d, -eps/d])
    return eps, float(np.sqrt(max(grad @ cov @ grad, 0)))


def radial_reference(radius, beta, sigma, target: float):
    """Interpolate one measured radial profile; never extrapolate.

    Endpoint errors assumed independent. Interpolation-model uncertainty and
    shared calibration errors are not included. Do not mix profiles at different purge.
    """
    r, b, s = (np.asarray(v, float) for v in (radius, beta, sigma))
    if r.ndim != 1 or r.shape != b.shape or b.shape != s.shape or len(r) < 2:
        raise ValueError("Need >=2 equally sized 1-D profile arrays.")
    if not all(np.isfinite(v).all() for v in (r,b,s)) or (s < 0).any() or not np.isfinite(target):
        raise ValueError("Profile and target must be finite; sigma >= 0.")
    order = np.argsort(r); r, b, s = r[order], b[order], s[order]
    if (np.diff(r) <= 0).any() or not r[0] <= target <= r[-1]:
        raise ValueError("Duplicate radii or requested extrapolation.")
    j = min(max(np.searchsorted(r, target, side='right')-1, 0), len(r)-2)
    w = (target-r[j])/(r[j+1]-r[j])
    return float((1-w)*b[j]+w*b[j+1]), float(np.hypot((1-w)*s[j], w*s[j+1]))
