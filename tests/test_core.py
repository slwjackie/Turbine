import json

import numpy as np
import pytest

from turbine_poc.core import (Baseline, INPUTS, TARGETS, TERMS, QuadraticRSM,
                              best_on_slice, endpoint_contrasts, features,
                              in_support, load_k1, loo_diagnostics, metrics)


@pytest.fixture
def dataset():
    return load_k1()


@pytest.fixture
def model(dataset):
    return QuadraticRSM.fit(dataset[1],dataset[2])


@pytest.fixture
def baseline():
    return Baseline("illustrative","test only",{"1":None,"2":None,"3":.4,"4":.6,"5":None,"6":None})


def test_published_rows_and_units(dataset):
    ids,x,y,digest=dataset
    assert ids==[f"{i:02d}" for i in range(25)]
    np.testing.assert_array_equal(x[7],[0,0,-1,1])
    np.testing.assert_allclose(y[0],[-.3091,-.4692],atol=0,rtol=0)
    np.testing.assert_allclose(y[24],[-.4034,-.6399],atol=0,rtol=0)
    assert len(digest)==64


def test_feature_basis_and_rank(dataset):
    a=features(dataset[1])
    assert a.shape==(25,15)
    assert len(TERMS)==15
    assert np.linalg.matrix_rank(a)==15
    assert np.linalg.matrix_rank(a[1:])==14
    np.testing.assert_allclose((dataset[1][1:]**2).sum(axis=1),2)


def test_reproduction_numbers(dataset,model):
    m=metrics(dataset[2][:,0],model.predict(dataset[1])[:,0])
    assert m['r2']==pytest.approx(.9999170810809803,abs=1e-12)
    assert m['rmse']==pytest.approx(.0006476624635307034,abs=1e-12)


def test_center_loo_explicitly_undefined(dataset):
    rows=loo_diagnostics(*dataset[:3])
    assert rows[0]['status']=='unidentifiable_fold'
    assert rows[0]['training_rank']==14
    assert rows[0]['leverage']==pytest.approx(1)
    assert rows[0]['predicted_delta_eta_pp'] is None
    assert sum(r['status']=='ok' for r in rows)==24


def test_no_silent_minimum_norm_fit(dataset):
    with pytest.raises(ValueError,match='not identifiable'):
        QuadraticRSM.fit(dataset[1][1:],dataset[2][1:])


@pytest.mark.parametrize('x',[[0,0,0,0],[1,1,0,0],[.5,.5,.5,.5],[-1,0,0,1]])
def test_support_includes_design_hull(x):
    assert in_support(x)[0]


@pytest.mark.parametrize('x',[[1,1,1,1],[1.1,0,0,0],[.8,.8,.8,0]])
def test_box_alone_is_not_enough(x,model):
    assert not in_support(x)[0]
    with pytest.raises(ValueError,match='convex hull'):
        model.predict(x)


@pytest.mark.parametrize('x',[[np.nan,0,0,0],[np.inf,0,0,0],[0,0,0],[]])
def test_nonfinite_and_shape_guard(x):
    with pytest.raises(ValueError):
        features(x)


def test_pairwise_contrasts(dataset):
    rows=endpoint_contrasts(dataset[1],dataset[2])
    assert len(rows)==8
    assert all(r['matched_pairs']==6 for r in rows)
    eta={r['variable']:r['mean_plus_minus'] for r in rows if r['target']==TARGETS[0]}
    assert eta['x3']<0 and eta['x4']<0 and eta['x5']>0 and eta['x6']>0
    assert abs(eta['x6']/eta['x3'])==pytest.approx(.333,abs=.006)


def test_unequal_baselines_mass_conservation(baseline):
    lo,hi=baseline.alpha_bounds()
    assert (lo,hi)==pytest.approx((.2,.6))
    m3,m4,x=baseline.allocation(np.linspace(lo,hi,57))
    np.testing.assert_allclose(m3+m4,1,atol=1e-14)
    np.testing.assert_allclose(.4*x[:,0]+.6*x[:,1],0,atol=1e-14)
    assert not np.allclose(x[:,0]+x[:,1],0)
    np.testing.assert_allclose(x[:,2:],0)
    assert in_support(x).all()
    assert baseline.total6 is None


def test_total_six_is_not_the_partial_budget(baseline):
    baseline.flows.update({'1':.1,'2':.1,'5':.2,'6':.3})
    assert baseline.total6==pytest.approx(1.7)
    assert baseline.m34==pytest.approx(1)


@pytest.mark.parametrize('scale',[.5,.8,1,1.2,1.5])
def test_variable_budget(scale,baseline):
    lo,hi=baseline.alpha_bounds(scale)
    m3,m4,x=baseline.allocation(np.linspace(lo,hi,13),scale)
    np.testing.assert_allclose(m3+m4,scale,atol=1e-14)
    assert in_support(x).all()


@pytest.mark.parametrize('scale',[0,-1,.49,1.51,float('nan')])
def test_impossible_budget_rejected(scale,baseline):
    with pytest.raises(ValueError):
        baseline.alpha_bounds(scale)


def test_exact_optimum_maximizes_less_negative_efficiency(baseline,model):
    alpha,value=best_on_slice(model,baseline)
    a=np.linspace(*baseline.alpha_bounds(),10001)
    values=model.predict(baseline.allocation(a)[2])[:,0]
    assert value>=values.max()-1e-12
    assert alpha==pytest.approx(.2917022032693679,abs=1e-10)
    assert value>model.predict([0,0,0,0])[0,0]


def test_template_rejects_unknown_design_flows(tmp_path):
    path=tmp_path/'bad.json'
    path.write_text(json.dumps({'source_kind':'lab_provided','provenance':'lab',
                               'flow_kg_s':{str(i):None for i in range(1,7)}}))
    with pytest.raises(ValueError,match='design flows'):
        Baseline.load(path)


def test_constant_target_r2_is_not_invented():
    assert metrics([1,1],[1,1])['r2'] is None
