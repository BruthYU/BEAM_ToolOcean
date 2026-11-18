"""Meaningful data-integrity tests; run in the cluster analysis environment."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy.sparse import csr_matrix
from tooluniverse.aging_data_tool import AgingIOTool,AgingPreprocessTool,AgingPseudobulkTool
from tooluniverse.exceptions import ToolValidationError


def tool(cls):
    configs=json.loads((Path(__file__).parents[2]/'examples/aging/schema_data.json').read_text())
    return cls(next(x for x in configs if x['type']==cls.__name__))


@pytest.fixture
def source(tmp_path):
    x=np.array([[1,2,3],[2,3,4],[4,5,6],[7,8,9],[3,4,5],[6,7,8]])
    obs=pd.DataFrame({'donor_id':['001']*3+['002']*3,'cell_type':['T','T','B','T','T','B'],
                      'group':['control']*3+['treated']*3},index=[f'c{i}' for i in range(6)])
    a=ad.AnnData(csr_matrix(x),obs=obs,var=pd.DataFrame(index=['A','B','C']))
    path=tmp_path/'counts.h5ad';a.write_h5ad(path)
    return a,path


def test_pseudobulk_preserves_sums_and_ids(source,tmp_path):
    a,path=source
    result=tool(AgingPseudobulkTool).run({'request':{'operation':'aggregate_counts','input_path':str(path),'output_dir':str(tmp_path),'metadata_columns':['group']}})
    p=ad.read_h5ad(result['expression']['path'])
    np.testing.assert_array_equal(np.asarray(p.X.sum(axis=0)),np.asarray(a.X.sum(axis=0)))
    assert set(p.obs.donor_id)=={'001','002'}
    for row in range(p.n_obs):
        mask=(a.obs.donor_id==p.obs.iloc[row].donor_id)&(a.obs.cell_type==p.obs.iloc[row].cell_type)
        np.testing.assert_array_equal(p.X[row].toarray(),np.asarray(a.X[mask.to_numpy()].sum(axis=0)))


def test_conflicting_pseudobulk_metadata_rejected(source,tmp_path):
    a,path=source;a.obs.loc['c1','group']='treated';a.write_h5ad(path)
    with pytest.raises(ToolValidationError,match='constant'):
        tool(AgingPseudobulkTool).run({'request':{'operation':'aggregate_counts','input_path':str(path),'output_dir':str(tmp_path),'metadata_columns':['group']}})


def test_import_preserves_leading_zero_ids(tmp_path):
    csv=tmp_path/'bulk.csv';csv.write_text('gene,001,002\nA,1,2\nB,3,4\n')
    result=tool(AgingIOTool).run({'request':{'operation':'import_expression','input_path':str(csv),'output_dir':str(tmp_path),'orientation':'genes_by_observations','species':'human','tissue':'blood','expression_scale':'counts','gene_id_type':'SYMBOL','data_kind':'bulk'}})
    a=ad.read_h5ad(result['expression']['path']);assert list(a.obs_names)==['001','002']
    np.testing.assert_array_equal(a.X,[[1,3],[2,4]])


def test_qc_and_export_chain(source,tmp_path):
    a,path=source
    result=tool(AgingPreprocessTool).run({'request':{'operation':'filter_qc','input_path':str(path),'output_dir':str(tmp_path),'min_counts':10}})
    final=tool(AgingIOTool).run({'request':{'operation':'export_results','input_path':result['expression']['path'],'output_dir':str(tmp_path),'genes':['A']}})
    frame=pd.read_csv(final['table']['path']);assert len(frame)==4
    assert list(frame['expression:A'])==[4,7,3,6]


def test_no_silent_duplicate_gene_mapping(source,tmp_path):
    _,path=source;mapping=tmp_path/'map.csv';mapping.write_text('source_id,target_id\nA,C\n')
    with pytest.raises(ToolValidationError,match='duplicate'):
        tool(AgingPreprocessTool).run({'request':{'operation':'align_gene_ids','input_path':str(path),'output_dir':str(tmp_path),'mapping_path':str(mapping)}})


def test_pseudobulk_not_accepted_as_bulk_model_input(source,tmp_path):
    _,path=source
    with pytest.raises(ToolValidationError,match='bulk'):
        tool(AgingPreprocessTool).run({'request':{'operation':'prepare_model_input','input_path':str(path),'output_dir':str(tmp_path),'model':'RNAAgeCalc','species':'human','expression_scale':'counts','data_kind':'pseudobulk'}})


def test_normalization_preserves_raw_counts_and_metadata(source,tmp_path):
    a,path=source
    result=tool(AgingPreprocessTool).run({'request':{'operation':'normalize_expression','input_path':str(path),'output_dir':str(tmp_path),'target_sum':100}})
    prepared=ad.read_h5ad(result['expression']['path'])
    np.testing.assert_array_equal(prepared.X.toarray(),a.X.toarray())
    expected=np.log1p(a.X.toarray()/np.asarray(a.X.sum(axis=1))*100)
    np.testing.assert_allclose(prepared.layers[result['layer']].toarray(),expected)
    with pytest.raises(ToolValidationError,match='already exists'):
        tool(AgingPreprocessTool).run({'request':{'operation':'normalize_expression','input_path':result['expression']['path'],'output_dir':str(tmp_path)}})


def test_metadata_inspection_and_gene_coverage(source,tmp_path):
    a,path=source
    meta=tmp_path/'metadata.csv';meta.write_text('cell_id,age\nc0,40\nc1,40\nc2,40\nc3,70\nc4,70\nc5,70\n')
    attached=tool(AgingIOTool).run({'request':{'operation':'attach_metadata','input_path':str(path),'metadata_path':str(meta),'id_key':'cell_id','output_dir':str(tmp_path)}})
    inspected=tool(AgingIOTool).run({'request':{'operation':'inspect_dataset','input_path':attached['expression']['path']}})
    assert 'age' in inspected['obs_columns']
    features=tmp_path/'features.csv';features.write_text('gene\nA\nB\nD\nE\n')
    with pytest.raises(ToolValidationError,match='coverage'):
        tool(AgingPreprocessTool).run({'request':{'operation':'prepare_model_input','model':'SenCID','input_path':str(path),'output_dir':str(tmp_path),'species':'human','expression_scale':'counts','data_kind':'single_cell','feature_path':str(features)}})


def test_duplicate_csv_genes_rejected(tmp_path):
    path=tmp_path/'duplicate.csv';path.write_text('id,A,A\n001,1,2\n')
    with pytest.raises(ToolValidationError,match='unique column'):
        tool(AgingIOTool).run({'request':{'operation':'inspect_dataset','input_path':str(path)}})


def test_normalization_rejects_primary_matrix_layer_name(source,tmp_path):
    _,path=source
    with pytest.raises(ToolValidationError,match='reserved'):
        tool(AgingPreprocessTool).run({'request':{'operation':'normalize_expression',
            'input_path':str(path),'output_dir':str(tmp_path),'output_layer':'X'}})


@pytest.mark.parametrize('overrides',[
    {'donor_key':'n_cells'}, {'sample_key':'pseudobulk_id'},
    {'cell_type_key':'retained'}, {'metadata_columns':['n_cells']},
    {'metadata_columns':['pseudobulk_id']}, {'metadata_columns':['retained']},
])
def test_pseudobulk_rejects_reserved_output_collisions(source,tmp_path,overrides):
    _,path=source
    with pytest.raises(ToolValidationError,match='reserved output names'):
        tool(AgingPseudobulkTool).run({'request':{'operation':'aggregate_counts',
            'input_path':str(path),'output_dir':str(tmp_path),**overrides}})


@pytest.mark.parametrize('overrides',[
    {'cell_type_key':'donor_id'}, {'sample_key':'cell_type'},
])
def test_pseudobulk_rejects_aliased_biological_roles(source,tmp_path,overrides):
    _,path=source
    with pytest.raises(ToolValidationError,match='must differ'):
        tool(AgingPseudobulkTool).run({'request':{'operation':'aggregate_counts',
            'input_path':str(path),'output_dir':str(tmp_path),**overrides}})


def test_pseudobulk_inspection_is_strict_json(source,tmp_path):
    _,path=source
    aggregated=tool(AgingPseudobulkTool).run({'request':{'operation':'aggregate_counts',
        'input_path':str(path),'output_dir':str(tmp_path)}})
    inspected=tool(AgingIOTool).run({'request':{'operation':'inspect_dataset',
        'input_path':aggregated['expression']['path']}})
    serialized=json.dumps(inspected,allow_nan=False)
    restored=json.loads(serialized)
    assert restored['recorded_metadata']['aggregation_keys']==['donor_id','cell_type']
