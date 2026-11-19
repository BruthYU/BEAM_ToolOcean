"""Real small-data mechanism tests; execute on the cluster analysis worker."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from jsonschema import Draft202012Validator

from tooluniverse.aging_mechanism_tool import (
    AgingDifferentialExpressionTool, AgingPathwayTool, _design,
)
from tooluniverse.exceptions import ToolValidationError


@pytest.fixture
def tools():
    root = Path(__file__).resolve().parents[2]
    configs = json.loads((root / 'examples/aging/schema_mechanism.json').read_text())
    return {c['name']: cls(c) for c, cls in zip(configs, [AgingDifferentialExpressionTool, AgingPathwayTool])}


def call(tools, name, request):
    tool = tools[name]
    error = tool.validate_parameters({'request': request})
    assert error is None, str(error)
    result = tool.run({'request': request})
    Draft202012Validator(tool.tool_config['return_schema']).validate(result)
    assert Path(result['results']['path']).exists()
    return result


@pytest.fixture
def counts_data(tmp_path):
    rng = np.random.default_rng(91)
    n, g = 12, 120
    mean = rng.uniform(30, 120, (1, g)) * np.ones((n, 1))
    mean[6:, :12] *= 5
    values = rng.negative_binomial(30, 30 / (mean + 30))
    counts = pd.DataFrame(values, index=[f's{i:02d}' for i in range(n)], columns=[f'g{i}' for i in range(g)])
    counts['zero'] = 0
    counts.to_csv(tmp_path / 'counts.csv', index_label='sample_id')
    meta = pd.DataFrame({'sample_id':counts.index, 'donor_id':[f'd{i}' for i in range(n)],
                         'group':['control']*6+['treated']*6, 'score':np.arange(n)})
    meta.to_csv(tmp_path / 'metadata.csv', index=False)
    return {'input_path':str(tmp_path/'counts.csv'),'metadata_path':str(tmp_path/'metadata.csv'),
            'output_dir':str(tmp_path/'results'),'expression_scale':'counts'}, meta


def test_deseq_groups_and_gsea_serial(tools, counts_data, tmp_path):
    request, _ = counts_data
    result = call(tools,'AgingDifferentialExpression',dict(request,operation='compare_groups',compare_groups=['control','treated']))
    de = pd.read_csv(result['results']['path'])
    assert de.shape[0] == 121
    assert de.loc[de.gene_id == 'zero','status'].item() == 'low_count_filtered'
    assert de.loc[de.gene_id.isin([f'g{i}' for i in range(12)]),'log2_fold_change'].median() > 1
    assert de.padj.dropna().between(0,1).all()
    enrichment = call(tools,'AgingPathway',{'operation':'run_preranked_gsea','input_path':result['results']['path'],
        'gene_sets':{'induced':[f'g{i}' for i in range(12)],'other':[f'g{i}' for i in range(30,50)]},
        'gene_set_source':'synthetic known-signal fixture','n_permutations':50,'output_dir':str(tmp_path/'gsea')})
    gsea = pd.read_csv(enrichment['results']['path'])
    assert gsea.loc[gsea.gene_set == 'induced','es'].item() > 0
    assert enrichment['n_excluded_nonfinite'] >= 1


def test_deseq_continuous_and_paired(tools, counts_data, tmp_path):
    request, meta = counts_data
    result = call(tools,'AgingDifferentialExpression',dict(request,operation='associate_with_score',score_key='score'))
    assert result['n_genes_tested'] > 0
    meta['donor_id'] = [f'd{i}' for i in range(6)] * 2
    meta.to_csv(request['metadata_path'],index=False)
    result = call(tools,'AgingDifferentialExpression',dict(request,operation='analyze_paired_change',compare_groups=['control','treated']))
    assert result['n_donors'] == 6
    assert result['n_genes_tested'] > 0


def test_donor_and_design_validation(counts_data):
    request, meta = counts_data
    meta.index = meta.sample_id
    meta['donor_id'] = ['same'] * len(meta)
    with pytest.raises(ToolValidationError, match='one sample per donor'):
        _design(dict(request,operation='compare_groups',compare_groups=['control','treated']),meta)
    meta['donor_id'] = meta.sample_id
    meta['confounded'] = (meta.group == 'treated').astype(float)
    with pytest.raises(ToolValidationError,match='rank deficient'):
        _design(dict(request,operation='compare_groups',compare_groups=['control','treated'],covariates=['confounded']),meta)


def test_overrepresentation_known_probability(tools,tmp_path):
    req={'operation':'run_overrepresentation','selected_genes':['a','b'],'background_genes':['a','b','c','d'],
         'gene_sets':{'hit':['a','b'],'miss':['c','d']},'gene_set_source':'fixture','output_dir':str(tmp_path)}
    result=call(tools,'AgingPathway',req)
    table=pd.read_csv(result['results']['path']).set_index('gene_set')
    assert table.loc['hit','pvalue'] == pytest.approx(1/6)
    assert table.loc['hit','padj'] == pytest.approx(1/3)
    assert table.loc['miss','pvalue'] == 1
    req['selected_genes']=['absent']
    with pytest.raises(ToolValidationError,match='subset'):
        tools['AgingPathway'].run({'request':req})


def test_score_gene_sets_coverage_and_ids(tools,tmp_path):
    pd.DataFrame({'a':[1.,2.],'b':[3.,4.]},index=['001','002']).to_csv(tmp_path/'expr.csv')
    req={'operation':'score_gene_sets','input_path':str(tmp_path/'expr.csv'),'expression_scale':'log_normalized',
         'gene_sets':{'covered':['a','b'],'uncovered':['c']},'gene_set_source':'fixture','output_dir':str(tmp_path)}
    result=call(tools,'AgingPathway',req)
    table=pd.read_csv(result['results']['path'],dtype={'observation_id':str})
    assert table.observation_id.tolist()==['001','002']
    assert table.covered.tolist()==[2.,3.]
    assert table.uncovered.isna().all()
    req['expression_scale']='counts'
    assert tools['AgingPathway'].validate_parameters({'request':req}) is not None


def test_schema_rejects_missing_background_and_conflicting_gene_sets(tools,tmp_path):
    req={'operation':'run_overrepresentation','selected_genes':['a'],'gene_sets':{'set':['a']},
         'gene_set_source':'fixture','output_dir':str(tmp_path)}
    assert tools['AgingPathway'].validate_parameters({'request':req}) is not None
    req.update(background_genes=['a','b'],gmt_path='unused.gmt')
    assert tools['AgingPathway'].validate_parameters({'request':req}) is not None
