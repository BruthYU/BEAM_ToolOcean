"""Small serial conversion/reference tests; execute on the cluster worker."""
import json
import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse
from tooluniverse import ToolUniverse

pytestmark = pytest.mark.integration


@pytest.fixture
def universe(tmp_path):
    tu = ToolUniverse(workspace=str(tmp_path / 'workspace'), load_workspace=False)
    tu.load_tools(categories=['spatial_transcriptomics'], include_tools=['SpatialSingleCellIO'])
    yield tu
    tu.close()


def call(tu, operation, **kwargs):
    result = tu.run_one_function({'name': 'SpatialSingleCellIO', 'arguments': {'request': {'operation': operation, **kwargs}}})
    assert 'error' not in result, result
    assert result['operation'] == operation
    return result


@pytest.mark.parametrize('platform', ['Xenium', 'MERFISH', 'CosMx'])
def test_serial_conversion_annotation_reference(universe, tmp_path, platform):
    # Leading zero IDs and repeated local cell IDs must remain distinguishable.
    raw = pd.DataFrame({'fov': ['01', '01', '02', '02'], 'cell_id': ['001', '002', '001', '002'],
                        'x': [0, 3, 0, 3], 'y': [1, 1, 1, 1], 'G1': [2, 1, 5, 1], 'G2': [1, 5, 1, 6],
                        'NegativeControl': [100, 100, 100, 100]})
    path = tmp_path / 'counts.csv'
    raw.to_csv(path, index=False)
    imported = call(universe, 'import_cells', input_path=str(path), platform=platform,
                    sample_id='sample', cell_id_columns=['fov', 'cell_id'], gene_columns=['G1', 'G2'],
                    coordinate_columns=['x', 'y'], coordinate_unit='micrometer', coordinate_frame='local',
                    fov_key='fov', output_dir=str(tmp_path / 'results'))
    data = ad.read_h5ad(imported['output_path'])
    assert data.shape == (4, 2)
    np.testing.assert_array_equal(data.X.toarray(), raw[['G1', 'G2']].to_numpy())
    assert data.obs_names.is_unique
    assert set(data.obs['fov_id']) == {'01', '02'}
    assert 'cell_type' not in data.obs
    assert data.uns['spatial_single_cell']['coordinate_frame'] == 'local'
    assert data.var['gene_name'].tolist() == ['G1', 'G2']
    # Reverse annotation order verifies that alignment is by ID, never row order.
    annotations = raw[['fov', 'cell_id']].copy()
    annotations['label'] = ['A', 'B', 'A', 'B']
    annotations.iloc[::-1].to_csv(tmp_path / 'annotations.csv', index=False)
    annotated = call(universe, 'attach_annotations', input_path=imported['output_path'],
                     annotation_path=str(tmp_path / 'annotations.csv'), annotation_id_columns=['fov', 'cell_id'],
                     input_id_key='source_cell_id', label_column='label', output_dir=str(tmp_path / 'results'))
    target = ad.AnnData(sparse.csr_matrix([[3, 4]]), obs=pd.DataFrame(index=['spot']), var=pd.DataFrame(index=['G1', 'G3']))
    target.write_h5ad(tmp_path / 'spatial.h5ad')
    reference = call(universe, 'prepare_reference', input_path=annotated['output_path'], cell_type_key='cell_type',
                     spatial_path=str(tmp_path / 'spatial.h5ad'), subset_shared_genes=True, output_dir=str(tmp_path / 'results'))
    prepared = ad.read_h5ad(reference['reference_path'])
    assert prepared.shape == (4, 1)
    assert prepared.obs['cell_type'].tolist() == ['A', 'B', 'A', 'B']
    assert reference['shared_gene_count'] == 1
    assert 'spatial' in prepared.obsm


def test_reject_ambiguous_ids_and_missing_labels(universe, tmp_path):
    pd.DataFrame({'cell_id': ['1', '1'], 'G1': [1, 2]}).to_csv(tmp_path / 'bad.csv', index=False)
    result = universe.run_one_function({'name': 'SpatialSingleCellIO', 'arguments': {'request': {
        'operation': 'import_cells', 'platform': 'CosMx', 'sample_id': 's', 'input_path': str(tmp_path / 'bad.csv'),
        'gene_columns': ['G1'], 'output_dir': str(tmp_path / 'out')}}})
    assert 'error' in result
    data = ad.AnnData(sparse.csr_matrix([[1], [2]]), obs=pd.DataFrame(index=['1', '2']), var=pd.DataFrame(index=['G1']))
    data.layers['counts'] = data.X.copy()
    data.write_h5ad(tmp_path / 'cells.h5ad')
    pd.DataFrame({'cell_id': ['1'], 'label': ['T']}).to_csv(tmp_path / 'ann.csv', index=False)
    result = universe.run_one_function({'name': 'SpatialSingleCellIO', 'arguments': {'request': {
        'operation': 'attach_annotations', 'input_path': str(tmp_path / 'cells.h5ad'), 'annotation_path': str(tmp_path / 'ann.csv'),
        'label_column': 'label', 'output_dir': str(tmp_path / 'out')}}})
    assert 'error' in result


def test_xenium_native_h5_control_filter(universe, tmp_path):
    import h5py
    # Native 10x genes-by-cells CSC layout, including a negative-control feature.
    matrix = sparse.csc_matrix([[2, 3], [4, 5], [99, 99]])
    with h5py.File(tmp_path / 'cell_feature_matrix.h5', 'w') as f:
        m = f.create_group('matrix')
        for key, value in [('data', matrix.data), ('indices', matrix.indices), ('indptr', matrix.indptr), ('shape', matrix.shape)]:
            m.create_dataset(key, data=value)
        m.create_dataset('barcodes', data=np.array(['001', '002'], dtype='S'))
        features = m.create_group('features')
        for key, values in {'id': ['g1', 'g2', 'ctrl'], 'name': ['G1', 'G2', 'NegativeControl'],
                            'feature_type': ['Gene Expression', 'Gene Expression', 'Negative Control Probe'],
                            'genome': ['human', 'human', 'human']}.items():
            features.create_dataset(key, data=np.array(values, dtype='S'))
    pd.DataFrame({'cell_id': ['002', '001'], 'x_centroid': [2, 1], 'y_centroid': [4, 3]}).to_csv(tmp_path / 'cells.csv', index=False)
    result = call(universe, 'import_cells', input_path=str(tmp_path), platform='Xenium', sample_id='s',
                  coordinate_unit='micrometer', coordinate_frame='global', output_dir=str(tmp_path / 'out'))
    imported = ad.read_h5ad(result['output_path'])
    assert result['controls_removed'] == 1
    np.testing.assert_array_equal(imported.X.toarray(), [[2, 4], [3, 5]])
    np.testing.assert_array_equal(imported.obsm['spatial'], [[1, 3], [2, 4]])
