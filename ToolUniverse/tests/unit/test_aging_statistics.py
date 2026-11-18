"""Scientific contract tests, executed in the cluster worker runtime."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from tooluniverse.exceptions import ToolValidationError
from tooluniverse.aging_statistics_tool import (
    AgingJointSummaryTool, AgingAssociationTool, AgingCompositionTool,
    AgingValidationTool, _fit, _bh,
)


def table(tmp_path, frame, name='input.csv'):
    path = tmp_path / name
    pd.DataFrame(frame).to_csv(path, index=False)
    return str(path)


def call(cls, op, tmp_path, **kwargs):
    r = {'operation':op,'output_dir':str(tmp_path),**kwargs}
    return getattr(object.__new__(cls), op)(r)


def result(response):
    return pd.read_csv(response['table']['path'], dtype={'donor_id':str})


@pytest.fixture
def cohort(tmp_path):
    rng = np.random.default_rng(912)
    n = 16
    group = np.array(['control']*8+['treated']*8)
    age = rng.uniform(40,80,n)
    proportion = rng.uniform(.1,.8,n)
    score = 4*(group=='treated')+.12*age+2*proportion+rng.normal(0,.6,n)
    frame = pd.DataFrame({'donor_id':[f'{i:03}' for i in range(n)],'group':group,'age':age,
                          'score':score,'predicted_age':age+rng.normal(1,2,n),
                          'proportion':proportion,'batch':['A','B']*8})
    return frame, table(tmp_path,frame)


def test_ols_matches_independent_statsmodels_hc3(cohort):
    sm = pytest.importorskip('statsmodels.api')
    df,_ = cohort
    df['treatment']=(df.group=='treated').astype(int)
    actual=_fit(df.copy(),'score','treatment',['age','proportion'])
    expected=sm.OLS(df.score,sm.add_constant(df[['treatment','age','proportion']])).fit(cov_type='HC3',use_t=True)
    assert actual['estimate']==pytest.approx(expected.params['treatment'])
    assert actual['standard_error']==pytest.approx(expected.bse['treatment'])
    assert actual['p_value']==pytest.approx(expected.pvalues['treatment'])
    assert actual['ci_lower']==pytest.approx(expected.conf_int().loc['treatment',0])


def test_bh_known_values():
    np.testing.assert_allclose(_bh([.01,.04,.03,.002]),[.02,.04,.04,.008])


def test_group_age_and_model_associations(tmp_path,cohort):
    _,path=cohort
    group=result(call(AgingAssociationTool,'compare_groups',tmp_path,input_path=path,value_columns=['score'],group_key='group',contrast=['control','treated'],covariates=['age','batch'],categorical_covariates=['batch']))
    assert group.estimate.iloc[0]>3
    assert group.n_donors.iloc[0]==16
    age=result(call(AgingAssociationTool,'analyze_age_association',tmp_path,input_path=path,value_columns=['predicted_age'],predictor_key='age'))
    assert .8 < age.estimate.iloc[0]<1.2
    corr=result(call(AgingAssociationTool,'compare_model_readouts',tmp_path,input_path=path,readouts=['age','predicted_age']))
    assert corr.spearman_rho.iloc[0]>.9


def test_rejects_pseudoreplication_and_rank_deficiency(tmp_path,cohort):
    df,path=cohort
    dup=table(tmp_path,pd.concat([df,df.iloc[:1]]),'duplicate.csv')
    with pytest.raises(ToolValidationError,match='independent'):
        call(AgingAssociationTool,'compare_groups',tmp_path,input_path=dup,value_columns=['score'],group_key='group',contrast=['control','treated'])
    df['same_age']=df.age
    with pytest.raises(ToolValidationError,match='rank deficient'):
        _fit(df,'score','age',['same_age'])


def test_paired_change_and_incomplete_pairs(tmp_path):
    from scipy.stats import ttest_rel
    before=np.array([9.,11,12,15,17])
    after=before+np.array([1.,3,2,4,2])
    df=pd.DataFrame({'donor_id':['a','b','c','d','e']*2+['unpaired'],
                     'time':['before']*5+['after']*5+['before'],'score':[*before,*after,12]})
    response=call(AgingAssociationTool,'analyze_paired_change',tmp_path,input_path=table(tmp_path,df),time_key='time',contrast=['before','after'],value_columns=['score'])
    out=result(response)
    assert out.estimate.iloc[0]==pytest.approx(2.4)
    assert out.p_value.iloc[0]==pytest.approx(ttest_rel(after,before).pvalue)
    assert response['excluded_pairs'][0]['excluded_donors']==['unpaired']


def test_summary_retains_metadata_and_ids(tmp_path):
    df={'donor_id':['001','001','002','002'],'cell_id':['a','b','c','d'],
        'cell_type':['T','T','T','T'],'group':['old','old','young','young'],'score':[1,3,5,7]}
    path=table(tmp_path,df)
    for op in ['summarize_by_donor','summarize_by_cell_type']:
        out=result(call(AgingJointSummaryTool,op,tmp_path,input_path=path,value_columns=['score'],metadata_columns=['group'],observation_key='cell_id'))
        assert list(out.donor_id)==['001','002']
        assert list(out.score)==[2,6]
        assert list(out.n_observations)==[2,2]
    df['group']=['old','young','young','young']
    with pytest.raises(ToolValidationError,match='varies'):
        call(AgingJointSummaryTool,'summarize_by_donor',tmp_path,input_path=table(tmp_path,df,'bad.csv'),value_columns=['score'],metadata_columns=['group'])


def test_join_reports_missing_and_rejects_many_to_many(tmp_path):
    a=table(tmp_path,{'cell_id':['a','b'],'score':[1,2]},'a.csv')
    b=table(tmp_path,{'cell_id':['b','c'],'score':[4,5]},'b.csv')
    inputs=[{'name':'clock','path':a,'value_columns':['score']},{'name':'sen','path':b,'value_columns':['score']}]
    out=call(AgingJointSummaryTool,'join_model_results',tmp_path,join_keys=['cell_id'],inputs=inputs)
    assert out['joined_rows']==3
    assert out['coverage'][0]['missing_joined_rows']==1
    assert list(result(out).columns)==['cell_id','clock__score','sen__score']
    inputs[1]['path']=table(tmp_path,{'cell_id':['b','b'],'score':[4,5]},'dup.csv')
    with pytest.raises(ToolValidationError,match='independent'):
        call(AgingJointSummaryTool,'join_model_results',tmp_path,join_keys=['cell_id'],inputs=inputs)


def test_composition_zero_fills_absent_types(tmp_path):
    path=table(tmp_path,{'donor_id':['a','a','b'],'cell_type':['T','B','T'],'n':[3,1,4],'group':['old','old','young']})
    out=result(call(AgingCompositionTool,'summarize_proportions',tmp_path,input_path=path,count_key='n',metadata_columns=['group']))
    assert len(out)==4
    assert out.loc[(out.donor_id=='b')&(out.cell_type=='B'),'proportion'].iloc[0]==0
    assert np.allclose(out.groupby('donor_id').proportion.sum(),1)


def test_composition_group_and_sensitivity(tmp_path,cohort):
    df,path=cohort
    comp=df[['donor_id','group','proportion']].copy()
    comp['cell_type']='T'
    cp=table(tmp_path,comp,'composition.csv')
    out=result(call(AgingCompositionTool,'test_composition_changes',tmp_path,input_path=cp,group_key='group',contrast=['control','treated']))
    assert out.n_donors.iloc[0]==16
    out=result(call(AgingCompositionTool,'assess_composition_sensitivity',tmp_path,input_path=path,value_columns=['score'],group_key='group',contrast=['control','treated'],composition_columns=['proportion'],covariates=['age']))
    assert set(out.adjustment)=={'baseline','composition_adjusted'}


def test_prediction_metrics_and_validation(tmp_path,cohort):
    df,path=cohort
    out=result(call(AgingValidationTool,'evaluate_age_prediction',tmp_path,input_path=path,age_key='age',prediction_key='predicted_age'))
    assert out.mae.iloc[0]==pytest.approx(abs(df.predicted_age-df.age).mean())
    disc=table(tmp_path,{'readout':['a','b'],'estimate':[2,3],'p_value':[.01,.02]},'discovery.csv')
    ext=table(tmp_path,{'readout':['a','b'],'estimate':[1,-3],'p_value':[.001,.001]},'external.csv')
    out=result(call(AgingValidationTool,'validate_external_cohort',tmp_path,discovery_results_path=disc,external_results_path=ext))
    assert list(out.direction_and_fdr_supported)==[True,False]
    args=dict(input_path=path,value_columns=['score'],group_key='group',contrast=['control','treated'],covariates=['age'])
    robust=call(AgingValidationTool,'assess_robustness',tmp_path,**args)
    assert len(result(robust))==17
    assert not robust['failed_refits']
    perturb=result(call(AgingValidationTool,'evaluate_perturbation',tmp_path,**args))
    assert perturb.estimate.iloc[0]>3


def test_schema_rejects_unknown_arguments(tmp_path):
    root=Path(__file__).resolve().parents[2]
    configs=json.loads((root/'examples/aging/schema_statistics.json').read_text())
    cfg=next(c for c in configs if c['name']=='AgingAssociation')
    tool=AgingAssociationTool(cfg)
    assert tool.validate_parameters({'request':{'operation':'compare_groups','input_path':'x.csv','output_dir':str(tmp_path),'value_columns':['score'],'group_key':'group','contrast':['a','b'],'pretend':True}}) is not None


def test_join_metadata_then_serial_two_model_summary(tmp_path):
    a=table(tmp_path,{'cell_id':['a','b','c','d'],'donor_id':['001','001','002','002'],'cell_type':['T']*4,'group':['old','old','young','young'],'predicted_age':[60.,62.,30.,32.]},'clock.csv')
    b=table(tmp_path,{'cell_id':['a','b','c','d'],'donor_id':['001','001','002','002'],'cell_type':['T']*4,'SID1_score':[.8,.6,.2,.1]},'sen.csv')
    inputs=[{'name':'clock','path':a,'value_columns':['predicted_age'],'metadata_columns':['donor_id','cell_type','group']},
            {'name':'sen','path':b,'value_columns':['SID1_score'],'metadata_columns':['donor_id','cell_type']}]
    joined=call(AgingJointSummaryTool,'join_model_results',tmp_path,join_keys=['cell_id'],inputs=inputs)
    summary=call(AgingJointSummaryTool,'summarize_by_cell_type',tmp_path,input_path=joined['table']['path'],observation_key='cell_id',value_columns=['clock__predicted_age','sen__SID1_score'],metadata_columns=['group'])
    out=result(summary)
    assert list(out.donor_id)==['001','002']
    np.testing.assert_allclose(out.clock__predicted_age,[61.,31.])
    np.testing.assert_allclose(out.sen__SID1_score,[.7,.15])
    bad=pd.read_csv(b,dtype=str)
    bad.loc[0,'donor_id']='999'
    inputs[1]['path']=table(tmp_path,bad,'conflict.csv')
    with pytest.raises(ToolValidationError,match='Conflicting metadata'):
        call(AgingJointSummaryTool,'join_model_results',tmp_path,join_keys=['cell_id'],inputs=inputs)


def test_external_validation_rejects_reversed_contrast(tmp_path):
    a=table(tmp_path,{'readout':['age'],'estimate':[2.],'p_value':[.01],'contrast':['treated - control']},'discovery.csv')
    b=table(tmp_path,{'readout':['age'],'estimate':[2.],'p_value':[.01],'contrast':['control - treated']},'external.csv')
    with pytest.raises(ToolValidationError,match='contrast differs'):
        call(AgingValidationTool,'validate_external_cohort',tmp_path,discovery_results_path=a,external_results_path=b)


def test_statistics_return_schemas_reject_wrong_operation_payload(tmp_path,cohort):
    from jsonschema import Draft202012Validator, ValidationError
    root=Path(__file__).resolve().parents[2]
    configs=json.loads((root/'examples/aging/schema_statistics.json').read_text())
    cfg=next(c for c in configs if c['name']=='AgingAssociation')
    _,path=cohort
    out=AgingAssociationTool(cfg).run({'request':{'operation':'compare_groups','input_path':path,'output_dir':str(tmp_path),'value_columns':['score'],'group_key':'group','contrast':['control','treated']}})
    Draft202012Validator(cfg['return_schema']).validate(out)
    out['operation']='analyze_paired_change'
    with pytest.raises(ValidationError):
        Draft202012Validator(cfg['return_schema']).validate(out)
