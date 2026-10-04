"""Small reproducible command-line workflows; existing outputs are never overwritten."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from . import __version__
from .data import COLUMNS, synthetic_data, validate
from .physics import amplitude_from_point, required_flow
from .validation import evaluate


REFERENCES = {
    'Tang2024': {'doi':'10.1016/j.ijheatfluidflow.2024.109300',
                'implemented':'Eqs. 18 and 30-35; sealing effectiveness, not the pressure model'},
    'Choi2024': {'doi':'10.1115/1.4064230',
                'used':'Table 4, Fig. 9, measurement-location and uncertainty cautions'},
    'Vella2026': {'doi':'10.1115/1.4070053',
                 'used':'Linear Phi_A versus annulus swirl hypothesis; coefficients NOT transferred'},
}


def new_directory(path):
    p = Path(path)
    p.mkdir(parents=True,exist_ok=False)
    return p


def published_check(out):
    """Published scalar checks only; no reconstructed curves or invented errors."""
    out = new_directory(out)
    rows = []
    for cfg,eps,bc in [('baseline',.811,.387),('rotor_swirler',.908,.500)]:
        rows.append(dict(configuration=cfg,phi0=.074,epsilon=eps,beta_c=bc,r_over_b=.923,
                         pointwise_effective_phi_a=amplitude_from_point(.074,eps),
                         source='Choi2024 Table 4; doi:10.1115/1.4064230'))
    pd.DataFrame(rows).to_csv(out/'table4_pointwise_inversion.csv',index=False)
    pd.DataFrame([
        dict(configuration='baseline',published_phi_min=.151,r_over_b=.941),
        dict(configuration='rotor_swirler',published_phi_min=.117,r_over_b=.941),
    ]).to_csv(out/'fig9_thresholds_SEPARATE_RADIUS.csv',index=False)
    facts = dict(status='PUBLISHED_SCALAR_SANITY_CHECK_NOT_VALIDATION',
                 phi95_over_constant_phi_a=float(required_flow(1,.95)),
                 published_minimum_flow_reduction_fraction=1-.117/.151,
                 beta_annulus_approximately=1.15,beta_annulus_method='Estimated, not directly measured',
                 correction_coefficient_fitted=False,references=REFERENCES)
    (out/'summary.json').write_text(json.dumps(facts,indent=2),encoding='utf8')
    (out/'README.md').write_text('''# 공개 수치 sanity check — 새 실험 검증 아님

Table 4의 두 점을 IWM에 역대입한 유효 진폭입니다. 곡선 전체의 최적 적합값이 아닙니다.
같은 점을 역대입해서 되찾는 것은 항등적 계산이지 모델 검증이 아닙니다.

Table 4: Phi0=0.074, r/b=0.923. Fig.9의 minimum sealing flow: r/b=0.941.
서로 다른 위치/유량의 값을 짝지어 보정계수 k를 맞추지 않았습니다.
Fig.9 원곡선을 디지털화한 자료도 아니며 불확실성은 제공되지 않아 만들지 않았습니다.
공개 임계유량과 `Phi_A=Phi_min`의 모델상 정의를 자동으로 동일시하지 않습니다.

출처: Choi et al. (2024), doi:10.1115/1.4064230, Table 4 및 Fig.9.
역산식: Tang et al. (2024), doi:10.1016/j.ijheatfluidflow.2024.109300, Eqs.30–35.
''',encoding='utf8')
    return facts


def make_plots(result, out, label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    predictions = result['predictions']
    if predictions.empty:
        return
    conditions = sorted(predictions.heldout.unique())
    heldout = conditions[len(conditions)//2]
    p = predictions[predictions.heldout==heldout]
    fig,ax = plt.subplots(figsize=(8,5))
    for cfg,g in p.groupby('configuration',sort=True):
        observed = g.drop_duplicates('row_id').sort_values('phi0')
        ax.scatter(observed.phi0,observed.observed_epsilon,label=f'{cfg}: observations',s=22)
        for model,z in g.groupby('model',sort=True):
            z=z.sort_values('phi0')
            ax.plot(z.phi0,z.predicted_epsilon,label=f'{cfg}: {model}',
                    linestyle='--' if model=='annulus' else '-')
    ax.set(xlabel='Sealing flow parameter Phi0',ylabel='Sealing effectiveness',
           title=f'{label}\nWhole-condition holdout: {heldout}',ylim=(-.03,1.07))
    ax.legend(fontsize=8); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(out/'heldout_effectiveness.png',dpi=160); plt.close(fig)
    thresholds = result['thresholds']
    fig,ax=plt.subplots(figsize=(8,5))
    for (model,cfg),g in thresholds.groupby(['model','configuration'],sort=True):
        g=g.sort_values('flow_coefficient')
        ax.plot(g.flow_coefficient,g.predicted_phi_target,marker='o',label=f'{cfg}: {model}')
    obs=thresholds.drop_duplicates('curve_id')
    obs=obs[np.isfinite(obs.observed_phi_target)]
    for cfg,g in obs.groupby('configuration',sort=True):
        ax.scatter(g.flow_coefficient,g.observed_phi_target,marker='x',s=55,label=f'{cfg}: bracketed observation')
    target=thresholds.target.iloc[0]
    ax.set(xlabel='Flow coefficient',ylabel=f'Flow at effectiveness {target:g}',
           title=f'{label}\nHeld-out target-flow predictions (end folds extrapolate)')
    ax.legend(fontsize=8); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(out/'heldout_target_flow.png',dpi=160); plt.close(fig)
    fitted=result['curve_fits_DIAGNOSTIC_ONLY']
    fitted=fitted[fitted.status=='ok']
    if not fitted.empty:
        fig,ax=plt.subplots(figsize=(8,5))
        for cfg,g in fitted.groupby('configuration',sort=True):
            ax.errorbar(g.beta_c_ref,g.phi_a,xerr=g.sigma_beta_c_ref,yerr=g.se_phi_a,
                        fmt='o',capsize=3,label=cfg)
            for _,r in g.iterrows():
                ax.annotate(f'CF={r.flow_coefficient:.2f}',(r.beta_c_ref,r.phi_a),
                            xytext=(4,4),textcoords='offset points',fontsize=8)
        ax.set(xlabel='Reference wheel-space swirl beta_c_ref',ylabel='Curve-fitted effective Phi_A',
               title=f'{label}\nAll-data diagnostic; error bars are conditional 1-sigma only')
        ax.legend(fontsize=8); ax.grid(alpha=.25); fig.tight_layout()
        fig.savefig(out/'amplitude_vs_swirl_DIAGNOSTIC_ONLY.png',dpi=160); plt.close(fig)


def analyze(data, out, *, draws=0,seed=20261005,split='condition_id',target=.95,no_plots=False):
    d=validate(data)
    out=new_directory(out)
    # This snapshot may contain private measurements: outputs/ must stay untracked.
    snapshot=d.to_csv(index=False)
    (out/'input_snapshot.csv').write_text(snapshot,encoding='utf8')
    result=evaluate(d,draws=draws,seed=seed,split=split,target=target)
    for name,frame in result.items():
        frame.to_csv(out/f'{name}.csv',index=False)
    metrics=result['fold_metrics']
    pivot=metrics.pivot(index='heldout',columns='model',values='rmse_epsilon')
    complete=pivot.dropna() if set(('annulus','swirl')).issubset(pivot.columns) else pd.DataFrame()
    scores={} if complete.empty else {m:float(complete[m].mean()) for m in ('annulus','swirl')}
    diagnostics=result['curve_fits_DIAGNOSTIC_ONLY']
    ok=diagnostics[diagnostics.status=='ok']
    poor=int((ok.reduced_chi2 > 3).sum()) if not ok.empty else 0
    synthetic=d.source_kind.iloc[0]=='synthetic'
    label='SYNTHETIC SOFTWARE DEMO - NOT EXPERIMENTAL VALIDATION' if synthetic else 'RESEARCH ANALYSIS - NOT ENGINE CERTIFICATION'
    metadata=dict(status=label,created_utc=datetime.now(timezone.utc).isoformat(),
        package_version=__version__,python=platform.python_version(),
        versions={name:version(name) for name in ('numpy','scipy','pandas','matplotlib')},
        input_sha256=hashlib.sha256(snapshot.encode()).hexdigest(),seed=seed,draws=draws,
        split=split,target=target,source_kind=d.source_kind.iloc[0],
        successful_matched_folds=len(complete),total_folds=int(d[split].nunique()),
        mean_fold_rmse_on_matched_folds=scores,mc_failed_model_folds=len(result['mc_failures']),
        curves_with_reduced_chi2_above_3=poor,
        interval_type='conditional_measurement_perturbation_percentiles_NOT_confidence_or_prediction_intervals',
        uncertainty_exclusions=['shared sensor bias','radius/gap uncertainty','model discrepancy',
            'observed target-flow interpolation uncertainty','full measurement covariance'],references=REFERENCES)
    (out/'manifest.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding='utf8')
    lines=[f'# {label}', '', f'Source: `{d.source.iloc[0]}`',
           f'Whole-group split: `{split}`; matched successful folds: {len(complete)}/{d[split].nunique()}.',
           '', '## Held-out effectiveness error',
           '| Model | Mean fold RMSE, matched successful folds only |','|---|---:|']
    for model in ('annulus','swirl'):
        lines.append(f'| {model} | {scores[model]:.8f} |' if model in scores else f'| {model} | withheld |')
    lines += ['', 'These numbers are not evidence of a real swirler benefit when inputs are synthetic.',
        'A failed fold is listed in fold_metrics.csv; never report a partial score as complete cross-validation.',
        '', '## Interpretation and limits',
        f'- {poor} full-data diagnostic curves have reduced chi-square >3. This is a heuristic lack-of-fit flag, '
        'not a calibrated hypothesis test; x-error and shared sensor errors affect it.',
        '- M0 is a rig-specific fitted linear annulus-swirl law, not the published numerical coefficients of another rig.',
        '- M1 adds response coefficient k, plus a train-only baseline-swirl reference regression with two nuisance coefficients. Its sign is unconstrained. This is a candidate, not derived physics.',
        '- Both models use the same training curves. All centering and curve/scaling fits are train-only.',
        '- Held-out cavity swirl at a predeclared reference purge is measured auxiliary information. '
        'This is NOT a forecast requiring no measurements in the new condition.',
        '- First/last flow-coefficient holdouts extrapolate beyond training annulus swirl. See the extrapolation column.',
        '- Phi_target is the constant-amplitude model root. Model Phi_min=Phi_A is NOT Phi95 or a measured saturation threshold.',
        '- Observed target flow is only reported for an unambiguous measured bracket, using linear interpolation.',
        '- Full-data amplitude fits are diagnostics, never independent predictions or inputs to held-out fits.',
        '', '## Uncertainty',
        f'{draws} measurement-perturbation draws; {len(result["mc_failures"])} failed model-fold refits.',
        'Reported 2.5/97.5 percentiles are sensitivity intervals for predicted mean effectiveness at nominal test Phi0 '
        'and predicted target flow. They are NOT calibrated confidence intervals or future-observation prediction intervals.',
        'At least 20 valid draws and 80% success are required; otherwise bounds are withheld. '
        'Use several hundred draws and check quantile stability for a real analysis.',
        'The draw model uses independent epsilon/Phi0 point errors, shared annulus error per condition, '
        'and shared reference-cavity-swirl error per curve. Positive Phi0 uses moment-matched lognormal draws.',
        'Shared instrument bias, geometry uncertainty, model discrepancy and observed-threshold uncertainty are not propagated. '
        'Therefore prediction improvements cannot yet be claimed statistically significant.',
        '', '## Sources', *[f'- {k}: doi:{v["doi"]}' for k,v in REFERENCES.items()],
        '', '## Outputs',
        '`predictions.csv`, `thresholds.csv`, `fold_metrics.csv`, `fold_membership.csv`, `coefficients.csv`, '
        '`curve_fits_DIAGNOSTIC_ONLY.csv`, `pointwise_amplitudes_DIAGNOSTIC_ONLY.csv`, `mc_failures.csv`, and the reproducibility manifest.',
        'The input snapshot can contain confidential measurements. Do not commit outputs or private_data to a public repository.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    if not no_plots:
        make_plots(result,out,label)
    return metadata


def main(argv=None):
    parser=argparse.ArgumentParser(description='IWM pre-meeting research scaffold; not a validated engine model.')
    sub=parser.add_subparsers(dest='command',required=True)
    for cmd in ('demo','analyze'):
        p=sub.add_parser(cmd)
        p.add_argument('--out',required=True)
        p.add_argument('--draws',type=int,default=0)
        p.add_argument('--seed',type=int,default=20261005)
        p.add_argument('--target',type=float,default=.95)
        p.add_argument('--split',choices=['condition_id','configuration'],default='condition_id')
        p.add_argument('--no-plots',action='store_true')
        if cmd=='demo':
            p.add_argument('--scenario',choices=['swirl','null','purge_dependence'],default='swirl')
        else:
            p.add_argument('--data',required=True)
    sub.add_parser('published-check').add_argument('--out',required=True)
    sub.add_parser('template').add_argument('--out',required=True)
    args=parser.parse_args(argv)
    try:
        if args.command=='template':
            path=Path(args.out); path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x',encoding='utf8') as f:
                pd.DataFrame(columns=COLUMNS).to_csv(f,index=False)
            print(f'Empty template created: {path}'); return
        if args.command=='published-check':
            record=published_check(args.out)
        else:
            data=synthetic_data(args.scenario,args.seed) if args.command=='demo' else pd.read_csv(args.data)
            record=analyze(data,args.out,draws=args.draws,seed=args.seed,split=args.split,
                           target=args.target,no_plots=args.no_plots)
        print(json.dumps(record,indent=2,ensure_ascii=False))
    except (ValueError,OSError) as exc:
        parser.exit(2,f'Error: {exc}\n')


if __name__=='__main__':
    main()
