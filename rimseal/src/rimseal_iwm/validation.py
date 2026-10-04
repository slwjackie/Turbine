"""Whole-condition cross-validation and explicitly conditional error propagation."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .data import validate, perturb
from .calibration import fit_curves, fit_scaling, NotIdentifiable
from .physics import effectiveness, required_flow, amplitude_from_point

MODELS = ('annulus','swirl')


def observed_threshold(curve: pd.DataFrame, target: float=.95):
    """Measured curve interpolation ONLY, never a fitted IWM 'observation'.

    Repeats at identical purge are averaged. Ambiguous/non-bracketed crossings
    are withheld. This is a descriptive estimate, not a noise-free ground truth.
    """
    z = curve.groupby('phi0',sort=True).epsilon.mean()
    p, y = z.index.to_numpy(float), z.to_numpy(float)
    if len(p) < 2:
        return np.nan, 'insufficient_levels'
    up = np.flatnonzero((y[:-1] < target) & (y[1:] >= target))
    down = np.flatnonzero((y[:-1] >= target) & (y[1:] < target))
    if len(up) != 1 or len(down) > 0:
        return np.nan, 'ambiguous_crossing' if len(up)+len(down)>1 else 'not_bracketed'
    j = up[0]
    return float(p[j]+(target-y[j])*(p[j+1]-p[j])/(y[j+1]-y[j])), 'bracketed_linear_interpolation'


def cross_validate(data: pd.DataFrame, *, split: str='condition_id', target: float=.95,
                   evaluation_phi0=None) -> dict[str,pd.DataFrame]:
    """No row-wise random split. Test epsilon is used ONLY for scoring.

    evaluation_phi0 may fix prediction abscissae during uncertainty propagation;
    training data always use their own (possibly perturbed) purge measurements.
    """
    if split not in ('condition_id','configuration'):
        raise ValueError('Split whole condition_id or configuration, never individual rows.')
    if not np.isfinite(target) or not 0 < target < 1:
        raise ValueError('Evaluation target must be strictly between 0 and 1.')
    if data[split].nunique() < 2:
        raise ValueError('Cross-validation needs at least two held-out groups.')
    predictions, thresholds, coefficients, metrics, membership = [], [], [], [], []
    for heldout in sorted(data[split].unique()):
        train, test = data[data[split]!=heldout], data[data[split]==heldout]
        membership.append(dict(heldout=str(heldout), train_curves=';'.join(sorted(train.curve_id.unique())),
                               test_curves=';'.join(sorted(test.curve_id.unique()))))
        try:
            fitted = fit_curves(train)  # every fit and baseline reference is TRAIN-ONLY
            train_error = None
        except NotIdentifiable as exc:
            train_error = str(exc)
        for model in MODELS:
            try:
                if train_error:
                    raise NotIdentifiable(train_error)
                fit = fit_scaling(fitted,model)
                a = fit.predict(test.beta_a,test.beta_c_ref)
                p = test.phi0.to_numpy() if evaluation_phi0 is None else np.asarray(evaluation_phi0)[test.index]
                pred = effectiveness(p,a)
                err = pred-test.epsilon.to_numpy()
                metrics.append(dict(heldout=str(heldout),model=model,status='ok',reason='',n=len(test),
                    rmse_epsilon=float(np.sqrt(np.mean(err**2))),
                    standardized_rmse=float(np.sqrt(np.mean((err/test.sigma_epsilon)**2)))))
                coefficients.append(dict(heldout=str(heldout),**fit.record()))
                for pos, (idx,row) in enumerate(test.iterrows()):
                    predictions.append(dict(row_id=int(idx),heldout=str(heldout),model=model,
                        curve_id=row.curve_id,configuration=row.configuration,flow_coefficient=row.flow_coefficient,
                        phi0=float(p[pos]),observed_epsilon=row.epsilon,predicted_epsilon=float(pred[pos]),
                        phi_a_prediction=float(a[pos]),extrapolation=bool(fit.extrapolation(row.beta_a))))
                for key,g in test.groupby('curve_id',sort=True):
                    row = g.iloc[0]
                    aa = float(fit.predict(row.beta_a,row.beta_c_ref))
                    measured,status = observed_threshold(g,target)
                    phi_target = float(required_flow(aa,target))
                    thresholds.append(dict(heldout=str(heldout),model=model,curve_id=key,
                        configuration=row.configuration,flow_coefficient=row.flow_coefficient,target=target,
                        predicted_phi_target=phi_target,model_phi_min=aa,
                        observed_phi_target=measured,observation_status=status,
                        phi_target_error=phi_target-measured,
                        extrapolation=bool(fit.extrapolation(row.beta_a))))
            except NotIdentifiable as exc:
                metrics.append(dict(heldout=str(heldout),model=model,status='unidentified',reason=str(exc),n=len(test),
                                    rmse_epsilon=np.nan,standardized_rmse=np.nan))
    return {name:pd.DataFrame(rows) for name,rows in dict(predictions=predictions,thresholds=thresholds,
            coefficients=coefficients,fold_metrics=metrics,fold_membership=membership).items()}


def _intervals(nominal, samples, keys, value, draws):
    result = nominal.copy()
    if result.empty:
        return result
    source = pd.concat(samples,ignore_index=True) if samples else pd.DataFrame()
    lookup = {} if source.empty else {k:g[value].to_numpy(float) for k,g in source.groupby(keys,sort=False)}
    lo,hi,count = [],[],[]
    for _,row in result.iterrows():
        key = tuple(row[k] for k in keys)
        values = lookup.get(key,np.array([]))
        values = values[np.isfinite(values)]
        count.append(len(values))
        # Do not present percentiles of an overwhelmingly failed calculation.
        if draws >= 20 and len(values) >= max(20,int(np.ceil(.8*draws))):
            low,high = np.quantile(values,[.025,.975])
        else:
            low,high = np.nan,np.nan
        lo.append(low); hi.append(high)
    result[value+'_mc_low'] = lo
    result[value+'_mc_high'] = hi
    result['mc_valid_draws'] = count
    result['mc_requested_draws'] = draws
    return result


def evaluate(data: pd.DataFrame, *, draws: int=0, seed: int=20261005,
             split: str='condition_id', target: float=.95):
    """Monte Carlo SENSITIVITY intervals, not confidence/validated prediction bands.

    Refit training curves, baseline reference trend and scaling on every draw.
    Annulus errors are shared across baseline/swirler of the same condition.
    Held-out measured swirl is an auxiliary input (not a zero-measurement forecast).
    Predictions use nominal test Phi0. Observed target-flow interpolation is not
    given a confidence interval in this release; its uncertainty is not ignored
    by claiming predicted-minus-observed errors are statistically significant.
    """
    if not isinstance(draws,int) or draws < 0:
        raise ValueError('draws must be a nonnegative integer.')
    d = validate(data)
    result = cross_validate(d,split=split,target=target)
    nominal_p = d.phi0.to_numpy()
    samples_p,samples_t,failures = [],[],[]
    rng = np.random.default_rng(seed)
    for draw in range(draws):
        q = perturb(d,rng)
        sample = cross_validate(q,split=split,target=target,evaluation_phi0=nominal_p)
        for name, destination in [('predictions',samples_p),('thresholds',samples_t)]:
            frame = sample[name].copy()
            if not frame.empty:
                frame['draw'] = draw
                destination.append(frame)
        failed = sample['fold_metrics'][sample['fold_metrics'].status!='ok']
        if not failed.empty:
            failed = failed.copy(); failed['draw'] = draw
            failures.append(failed)
    result['predictions'] = _intervals(result['predictions'],samples_p,['row_id','model'],
                                      'predicted_epsilon',draws)
    result['thresholds'] = _intervals(result['thresholds'],samples_t,['curve_id','model'],
                                     'predicted_phi_target',draws)
    result['mc_failures'] = pd.concat(failures,ignore_index=True) if failures else pd.DataFrame(
        columns=['heldout','model','status','reason','n','rmse_epsilon','standardized_rmse','draw'])
    # Full-data fits are labeled diagnostics and NEVER fed to cross-validation.
    diagnostic = []
    for key,g in d.groupby('curve_id',sort=True):
        try:
            row = fit_curves(g).iloc[0].to_dict(); row['status']='ok'; row['reason']=''
        except NotIdentifiable as exc:
            row = dict(curve_id=key,status='unidentified',reason=str(exc))
        diagnostic.append(row)
    result['curve_fits_DIAGNOSTIC_ONLY'] = pd.DataFrame(diagnostic)
    pointwise = []
    for idx,row in d.iterrows():
        informative = row.phi0 > 0 and 2*row.sigma_epsilon < row.epsilon < 1-2*row.sigma_epsilon
        pointwise.append(dict(row_id=idx,curve_id=row.curve_id,phi0=row.phi0,epsilon=row.epsilon,
            beta_c_ref=row.beta_c_ref,
            phi_a_point=amplitude_from_point(row.phi0,row.epsilon) if informative else np.nan,
            status='interior_point_inverse_not_validation' if informative else 'withheld_near_endpoint'))
    result['pointwise_amplitudes_DIAGNOSTIC_ONLY'] = pd.DataFrame(pointwise)
    return result
