import csv
import json
from pathlib import Path

import numpy as np
import pytest

from turbine_poc.cli import main
from turbine_poc.core import CONTEXT, REFERENCE
from turbine_poc.workflows import (OBS_FIELDS, paired_comparison, read_observations,
                                   reproduce, sweep, validate_cfd, write_csv)

ROOT=Path(__file__).parents[1]
EXAMPLE=ROOT/'configs/illustrative.json'


def observation(case='test',allocation='B',**overrides):
    row=dict(zip(OBS_FIELDS,['']*len(OBS_FIELDS)))
    row.update(case_id=case,allocation_id=allocation,context_id=CONTEXT,
               phase_id='paper_zero_offset',mesh_id='paper_medium',reference_id=REFERENCE,
               x3=.2,x4=-.1,x5=0,x6=0,delta_eta_pp=-.32,delta_torque_nm=-.48,
               sealing_effectiveness=.85,ingress_kg_s=.01,
               status='completed',origin='cfd',role='holdout')
    row.update(overrides)
    return row


def test_demo_cli_and_artifacts(tmp_path):
    out=tmp_path/'demo'
    assert main(['demo','--baseline',str(EXAMPLE),'--out',str(out),'--no-plots'])==0
    metrics=json.loads((out/'reproduction/metrics.json').read_text())
    assert metrics['independent_cfd_validation']=='not_performed'
    assert not metrics['loo']['complete_loo_score_available']
    summary=json.loads((out/'allocation/allocation_summary.json').read_text())
    assert summary['source_kind']=='illustrative'
    assert summary['sealing_constraint']=='NOT_EVALUATED_NO_DATA'
    assert summary['k1_six_passage_total_kg_s'] is None
    rows=list(csv.DictReader((out/'allocation/cfd_results_template.csv').open()))
    assert 5<=len(rows)<=7
    assert all(r['delta_eta_pp']=='' and r['status']=='pending' for r in rows)
    assert all(r['sealing_effectiveness']=='' for r in rows)


def test_output_overwrite_refused(tmp_path):
    out=tmp_path/'run'
    reproduce(out,plots=False)
    assert main(['reproduce','--out',str(out),'--no-plots'])==2


def test_pending_template_cannot_be_validated(tmp_path):
    out=tmp_path/'run'
    sweep(EXAMPLE,out,plots=False)
    with pytest.raises(ValueError,match='No completed'):
        read_observations(out/'cfd_results_template.csv')


def test_new_cfd_validation_path_with_synthetic_test_fixture(tmp_path):
    # Fixture only, to test the importer. NEVER included among scientific results.
    csvpath=tmp_path/'fixture.csv'; write_csv(csvpath,[observation()],OBS_FIELDS)
    report=validate_cfd(csvpath,tmp_path/'validation')
    assert report['scores']['holdout']['n']==1
    assert report['scores']['verification']['n']==0


@pytest.mark.parametrize('change',[
    {'context_id':'H2'}, {'phase_id':'other_phase'}, {'mesh_id':'finer'},
    {'reference_id':'OTHER_REFERENCE'}, {'x3':0,'x4':0},
    {'origin':'rsm'}, {'sealing_effectiveness':1.2}, {'ingress_kg_s':-.1},
    {'delta_eta_pp':float('nan')}
])
def test_invalid_validation_or_physics_rejected(tmp_path,change):
    path=tmp_path/'invalid.csv'; write_csv(path,[observation(**change)],OBS_FIELDS)
    with pytest.raises(ValueError):
        validate_cfd(path,tmp_path/'out')


def test_paired_subtraction_preserves_phase_matching(tmp_path):
    # Arbitrary synthetic data for software testing, not K1 or new CFD measurements.
    rows=[]
    for phase,base,gain in [('p0',-.3,.02),('p1',-.5,.025)]:
        rows.extend([observation(case=f'{phase}_base',allocation='A',phase_id=phase,
                                 delta_eta_pp=base,role='verification'),
                     observation(case=f'{phase}_candidate',allocation='B',phase_id=phase,
                                 delta_eta_pp=base+gain)])
    path=tmp_path/'synthetic.csv'; write_csv(path,rows,OBS_FIELDS)
    report=paired_comparison(path,'A',tmp_path/'paired')
    assert report['n_pairs']==2
    out=list(csv.DictReader((tmp_path/'paired/paired_differences.csv').open()))
    np.testing.assert_allclose([float(r['difference_delta_eta_pp']) for r in out],[.02,.025])


def test_paired_requires_every_matched_reference(tmp_path):
    path=tmp_path/'fixture.csv'; write_csv(path,[observation()],OBS_FIELDS)
    with pytest.raises(ValueError,match='Missing matched baseline'):
        paired_comparison(path,'A',tmp_path/'out')


def test_duplicate_observation_refused(tmp_path):
    path=tmp_path/'fixture.csv'; write_csv(path,[observation(),observation()],OBS_FIELDS)
    with pytest.raises(ValueError,match='Duplicate'):
        read_observations(path)
