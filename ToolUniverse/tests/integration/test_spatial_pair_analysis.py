"""Small numerical post-analysis checks; run on the cluster worker, not developer host."""
import json
from pathlib import Path
import pytest

pytestmark = pytest.mark.integration


@pytest.fixture
def tools():
    from tooluniverse.spatial_pair_analysis_tool import SpatialPairValidationTool, SpatialPairSummaryTool
    return SpatialPairValidationTool({}), SpatialPairSummaryTool({})


@pytest.fixture
def cells(tmp_path):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    rows, xy = [], []
    for sample in ['s1','s2','s3','s4']:
        for frame in ['f1','f2']:
            for j, x in enumerate([0, 1, 10, 11]):
                rows.append({'id':f'{sample}_{frame}_{j}', 'sample_id':sample, 'frame_id':frame, 'cell_type':['T','B','T','B'][j]})
                xy.append([x, 0])
    obs = pd.DataFrame(rows).set_index('id')
    data = ad.AnnData(csr_matrix(np.ones((len(obs), 2))), obs=obs, var=pd.DataFrame(index=['G1','G2']), obsm={'spatial':np.array(xy,dtype=float)})
    path = tmp_path/'cells.h5ad'
    data.write_h5ad(path)
    return path


def benchmark(tool, cells, tmp_path):
    return tool.run({'request':dict(operation='build_pair_benchmark', input_path=str(cells), sample_key='sample_id', frame_key='frame_id', coordinate_unit='micrometer', radius=2, pairs_per_class=10, output_dir=str(tmp_path))})


def test_serial_postanalysis(tools, cells, tmp_path):
    import pandas as pd
    validation, summary_tool = tools
    result = benchmark(validation, cells, tmp_path)
    truth = pd.read_csv(result['benchmark_path'])
    assert len(truth) == 48  # 8 frames, 6 distinct unordered pairs in each.
    assert truth.is_neighbor.sum() == 16
    assert all(a.split('_')[:2] == b.split('_')[:2] for a,b in zip(truth.cell_id_a, truth.cell_id_b))
    pred = truth[['sample_id','cell_id_a','cell_id_b']].copy()
    # Synthetic oracle is exclusively a validation fixture, not model performance evidence.
    pred['proximity_score'] = truth.is_neighbor * .8 + .1
    prediction = tmp_path/'pred.csv'
    pred.to_csv(prediction,index=False)
    evaluated = validation.run({'request':dict(operation='evaluate_pairs', prediction_path=str(prediction), benchmark_path=result['benchmark_path'], threshold=.5, output_dir=str(tmp_path))})
    assert evaluated['metrics']['auroc'] == 1
    assert evaluated['metrics']['auprc'] == 1
    assert Path(evaluated['figure_path']).is_file()
    summary = summary_tool.run({'request':dict(operation='summarize_type_pairs', input_path=str(cells), prediction_path=str(prediction), benchmark_path=result['benchmark_path'], cell_type_key='cell_type', sample_key='sample_id', frame_key='frame_id', sampling_description='Enumerated all within-frame pairs in synthetic fixture', output_dir=str(tmp_path))})
    metadata=tmp_path/'metadata.csv'
    pd.DataFrame(dict(sample_id=['s1','s2','s3','s4'],patient_id=['p1','p2','p3','p4'],group=['a','a','b','b'])).to_csv(metadata,index=False)
    comparison=summary_tool.run({'request':dict(operation='compare_groups', summary_path=summary['summary_path'],metadata_path=str(metadata),compare_groups=['a','b'],n_permutations=99,output_dir=str(tmp_path))})
    report=json.load(open(comparison['report_path']))
    assert report['number_testable'] == 3
    assert all(r['n_patients_by_group'] == {'a':2,'b':2} for r in report['results'])
    assert all(r['p_value'] == 1 for r in report['results'])


def test_exact_pair_join_rejects_missing_and_frame_mismatch(tools, cells, tmp_path):
    import pandas as pd
    from tooluniverse.exceptions import ToolValidationError
    validation, _ = tools
    result=benchmark(validation,cells,tmp_path)
    truth=pd.read_csv(result['benchmark_path'])
    pred=truth[['sample_id','frame_id','cell_id_a','cell_id_b']].copy()
    pred['proximity_score']=.5
    path=tmp_path/'pred.csv'
    pred.iloc[:-1].to_csv(path,index=False)
    request=dict(operation='evaluate_pairs',prediction_path=str(path),benchmark_path=result['benchmark_path'],output_dir=str(tmp_path))
    with pytest.raises(ToolValidationError,match='match exactly'):
        validation.run({'request':request})
    pred['frame_id']='wrong'
    pred.to_csv(path,index=False)
    with pytest.raises(ToolValidationError,match='frames disagree'):
        validation.run({'request':request})


def test_multiple_samples_do_not_create_fake_patients(tools,tmp_path):
    import pandas as pd
    _,summary=tools
    data=tmp_path/'summary.csv'
    meta=tmp_path/'meta.csv'
    pd.DataFrame(dict(sample_id=['s1','s2','s3','s4'],type_a=['B']*4,type_b=['T']*4,mean_score=[.9,.8,.1,.2])).to_csv(data,index=False)
    pd.DataFrame(dict(sample_id=['s1','s2','s3','s4'],patient_id=['p1','p1','p2','p2'],group=['a','a','b','b'])).to_csv(meta,index=False)
    result=summary.run({'request':dict(operation='compare_groups',summary_path=str(data),metadata_path=str(meta),compare_groups=['a','b'],output_dir=str(tmp_path))})
    report=json.load(open(result['report_path']))
    assert report['number_testable']==0
    assert report['results'][0]['status']=='insufficient_patients'


def test_local_frames_and_units_are_not_silently_mixed(tools,cells,tmp_path):
    import anndata as ad
    from tooluniverse.exceptions import ToolValidationError
    validation,_=tools
    data=ad.read_h5ad(cells)
    data.obs['fov_id']=data.obs['frame_id']
    data.uns['spatial_single_cell']={'coordinate_unit':'pixel','coordinate_frame':'local'}
    data.write_h5ad(cells)
    with pytest.raises(ToolValidationError,match='coordinate_unit'):
        benchmark(validation,cells,tmp_path)
    request=dict(operation='build_pair_benchmark',input_path=str(cells),sample_key='sample_id',frame_key='sample_id',coordinate_unit='pixel',radius=2,pairs_per_class=10,output_dir=str(tmp_path))
    with pytest.raises(ToolValidationError,match='field-of-view'):
        validation.run({'request':request})
