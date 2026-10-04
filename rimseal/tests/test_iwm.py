import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad
from rimseal_iwm.physics import (effectiveness,fluxes,required_flow,amplitude_from_point,
    mass_to_phi,renormalize_phi,tracer_effectiveness,radial_reference)
from rimseal_iwm.data import synthetic_data,validate,perturb
from rimseal_iwm.calibration import fit_curve,fit_curves,fit_scaling,NotIdentifiable
from rimseal_iwm.validation import cross_validate,evaluate,observed_threshold
from rimseal_iwm.cli import main,published_check


@pytest.mark.parametrize('p,a',[(0,.1),(.01,.1),(.099999,.1),(.1,.1),(.2,.1)])
def test_fluxes_against_independent_angular_quadrature(p,a):
    # Integrate the positive/negative radial wave, independent of the analytic expression.
    i,e=fluxes(p,a)
    eq=quad(lambda t:max(p+a*np.sin(t),0),0,2*np.pi,epsabs=1e-10,limit=200)[0]/(2*np.pi)
    iq=quad(lambda t:max(-(p+a*np.sin(t)),0),0,2*np.pi,epsabs=1e-10,limit=200)[0]/(2*np.pi)
    assert e==pytest.approx(eq,abs=1e-8)
    assert i==pytest.approx(iq,abs=1e-8)
    assert e-i==pytest.approx(p,abs=1e-14)


def test_limits_and_broadcasting():
    p=np.linspace(0,.3,201)
    e=effectiveness(p,.1)
    assert e[0]==0 and e[-1]==1
    assert np.all(np.diff(e)>=-1e-14)
    assert effectiveness(np.array([[0],[.2]]),[.1,.3]).shape==(2,2)


@pytest.mark.parametrize('p,a',[(-1,.1),(0,0),(1,-1),(np.nan,1),(1,np.inf)])
def test_invalid_physics(p,a):
    with pytest.raises(ValueError): effectiveness(p,a)


@pytest.mark.parametrize('target',[.01,.5,.95,.99,1.])
def test_target_roots(target):
    p=required_flow(.173,target)
    assert effectiveness(p,.173)==pytest.approx(target,abs=2e-12)
    assert required_flow(.173,1)==pytest.approx(.173)
    if target<1:
        assert amplitude_from_point(p,target)==pytest.approx(.173,rel=1e-10)


@pytest.mark.parametrize('p,e',[(0,.5),(.1,0),(.1,1),(.1,1.1)])
def test_inverse_rejects_nonunique_endpoints(p,e):
    with pytest.raises(ValueError): amplitude_from_point(p,e)


def test_units_and_tracer_covariance():
    p=mass_to_phi(.1,1.2,100,.3,.001)
    q=renormalize_phi(p,old_density=1.2,old_omega=100,old_radius=.3,old_gap=.001,
                       new_density=1.2,new_omega=100,new_radius=.6,new_gap=.001)
    assert q==pytest.approx(p/4)
    cov=np.diag([.01,.02,.03])**2
    e,s=tracer_effectiveness(.6,.1,1.1,cov)
    assert e==pytest.approx(.5)
    assert s==pytest.approx(np.sqrt(.01**2+(.5*.02)**2+(.5*.03)**2))
    assert tracer_effectiveness(1.2,.1,1.1)[0]>1  # don't clip observations
    with pytest.raises(ValueError): tracer_effectiveness(.6,.1,1.1,-cov)


def test_radial_interpolation_not_extrapolation():
    b,s=radial_reference([.9,.95],[.4,.6],[.01,.02],.925)
    assert b==pytest.approx(.5)
    assert s==pytest.approx(np.hypot(.005,.01))
    with pytest.raises(ValueError): radial_reference([.9,.95],[.4,.6],[.01,.02],.96)


def noiseless(scenario='swirl'):
    d=synthetic_data(scenario)
    for _,g in d.groupby('curve_id'):
        r=g.iloc[0]
        delta=r.beta_c_ref-(.44+.045*(r.beta_a-1))
        a=(.015+.115*r.beta_a)*np.exp((-2.5 if scenario=='swirl' else 0)*delta)
        d.loc[g.index,'epsilon']=effectiveness(g.phi0,a)
    return validate(d)


def test_amplitude_recovery_and_unidentified_plateau():
    d=noiseless(); g=d[d.curve_id==d.curve_id.iloc[0]]
    a=.015+.115*g.beta_a.iloc[0]
    fit=fit_curve(g)
    assert fit['phi_a']==pytest.approx(a,rel=1e-6)
    assert fit['se_phi_a']>0
    g=g.copy(); g.epsilon=1
    with pytest.raises(NotIdentifiable): fit_curve(g)


def test_candidate_recovers_synthetic_coefficient_not_scientific_evidence():
    f=fit_scaling(fit_curves(noiseless()),'swirl')
    assert f.k==pytest.approx(-2.5,abs=1e-5)
    assert f.record()['annulus_slope']==pytest.approx(.115,abs=1e-6)


def test_null_scenario_has_no_swirl_effect():
    d=noiseless('null'); f=fit_scaling(fit_curves(d),'swirl')
    assert f.k==pytest.approx(0,abs=1e-5)


def test_configuration_holdout_cannot_learn_swirl_correction():
    d=noiseless(); fitted=fit_curves(d[d.configuration=='baseline'])
    with pytest.raises(NotIdentifiable): fit_scaling(fitted,'swirl')
    cv=cross_validate(d,split='configuration')
    assert (cv['fold_metrics'].query("model=='swirl'").status=='unidentified').all()


def test_whole_condition_membership_and_no_test_outcome_leakage():
    d=noiseless(); held='CF0.40'
    base=cross_validate(d)
    changed=d.copy(); changed.loc[d.condition_id==held,'epsilon']+=.2
    altered=cross_validate(changed)
    for cv in (base,altered):
        row=cv['fold_membership'].query('heldout==@held').iloc[0]
        assert not set(row.train_curves.split(';'))&set(row.test_curves.split(';'))
        assert len(row.test_curves.split(';'))==2
    p=base['predictions'].query('heldout==@held').predicted_epsilon.to_numpy()
    q=altered['predictions'].query('heldout==@held').predicted_epsilon.to_numpy()
    np.testing.assert_array_equal(p,q)
    a=base['coefficients'].query('heldout==@held').reset_index(drop=True)
    b=altered['coefficients'].query('heldout==@held').reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b)


@pytest.mark.parametrize('column,value',[('density_ratio',1.5),('sigma_epsilon',0),
                                          ('r_swirl_over_b',.90),('phi_ref',.08)])
def test_input_definition_guards(column,value):
    d=synthetic_data(); d.loc[0,column]=value
    with pytest.raises(ValueError): validate(d)


def test_duplicate_condition_split_is_rejected():
    d=synthetic_data(); mask=(d.configuration=='rotor_swirler')&(d.condition_id=='CF0.30')
    d.loc[mask,'condition_id']='secretly_split_same_condition'
    with pytest.raises(ValueError): validate(d)


def test_correlated_feature_draws_and_reproducibility():
    d=validate(synthetic_data())
    a=perturb(d,np.random.default_rng(9)); b=perturb(d,np.random.default_rng(9))
    pd.testing.assert_frame_equal(a,b)
    assert a.groupby('condition_id').beta_a.nunique().max()==1
    assert a.groupby('curve_id').beta_c_ref.nunique().max()==1
    assert (a.phi0>0).all()


def test_observed_target_is_bracketed_not_fitted():
    d=pd.DataFrame({'phi0':[.01,.1,.2],'epsilon':[.1,.9,1.]})
    p,status=observed_threshold(d,.95)
    assert p==pytest.approx(.15) and status.startswith('bracketed')
    d.epsilon=[.1,.3,.5]
    assert np.isnan(observed_threshold(d,.95)[0])
    d=pd.DataFrame({'phi0':[.01,.1,.15,.2],'epsilon':[.1,.99,.9,1.]})
    assert observed_threshold(d,.95)[1]=='ambiguous_crossing'


def test_small_mc_does_not_issue_spurious_intervals():
    d=noiseless()
    a=evaluate(d,draws=2,seed=11); b=evaluate(d,draws=2,seed=11)
    pd.testing.assert_frame_equal(a['predictions'],b['predictions'])
    assert a['predictions'].predicted_epsilon_mc_low.isna().all()
    assert (a['predictions'].mc_valid_draws==2).all()


def test_cli_output_and_no_overwrite(tmp_path):
    out=tmp_path/'demo'
    main(['demo','--out',str(out),'--no-plots'])
    assert (out/'report.md').exists() and (out/'manifest.json').exists()
    with pytest.raises(SystemExit): main(['demo','--out',str(out),'--no-plots'])
    template=tmp_path/'private_data'/'measurements.csv'
    main(['template','--out',str(template)])
    assert pd.read_csv(template).empty


def test_published_scalars_are_not_paired_across_radii(tmp_path):
    out=tmp_path/'anchors'; facts=published_check(out)
    t=pd.read_csv(out/'table4_pointwise_inversion.csv')
    f=pd.read_csv(out/'fig9_thresholds_SEPARATE_RADIUS.csv')
    assert t.r_over_b.iloc[0]!=f.r_over_b.iloc[0]
    assert facts['correction_coefficient_fitted'] is False


def test_interval_success_threshold():
    from rimseal_iwm.validation import _intervals
    nominal=pd.DataFrame([dict(row_id=1,model='swirl',prediction=.8)])
    draws=[pd.DataFrame([dict(row_id=1,model='swirl',prediction=.8+i*.001)]) for i in range(20)]
    good=_intervals(nominal,draws,['row_id','model'],'prediction',20)
    assert .8 < good.prediction_mc_low.iloc[0] < good.prediction_mc_high.iloc[0] < .82
    bad=_intervals(nominal,draws,['row_id','model'],'prediction',30)
    assert bad.prediction_mc_low.isna().all() and bad.mc_valid_draws.iloc[0]==20


def test_pointwise_diagnostic_withholds_saturated_points():
    out=evaluate(noiseless(),draws=0)['pointwise_amplitudes_DIAGNOSTIC_ONLY']
    assert out.phi_a_point.isna().any()
    assert out.loc[out.epsilon==1,'phi_a_point'].isna().all()
