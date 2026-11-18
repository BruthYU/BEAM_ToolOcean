"""Input and artifact contracts; real checkpoint inference is exercised on the GPU worker."""
import hashlib
import json
from pathlib import Path

import pytest

from tooluniverse.exceptions import ToolValidationError
from tooluniverse.spatialformer_tool import SpatialFormerTool, _bundle


def test_bundle_rejects_corruption_and_wrong_mode(tmp_path):
    files = {}
    for key in ('checkpoint', 'config', 'tokenizer', 'gene_medians'):
        path = tmp_path / key
        path.write_bytes(b'original')
        files[key] = {'path': key, 'sha256': hashlib.sha256(b'original').hexdigest()}
    manifest = {'mode': 'single', 'model_id': 'test', 'source_commit': 'test',
                'package_version': 'test', 'files': files}
    (tmp_path / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ToolValidationError, match='pair'):
        _bundle(tmp_path, 'pair')
    (tmp_path / 'checkpoint').write_bytes(b'corrupt')
    with pytest.raises(ToolValidationError, match='checksum'):
        _bundle(tmp_path, 'single')


def test_pair_sample_mismatch_rejected_before_model_call(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import tooluniverse.spatialformer_tool as module

    data = ad.AnnData(np.ones((2, 3)), obs=pd.DataFrame({'sample_id': ['A', 'B']}, index=['cell_a', 'cell_b']))
    monkeypatch.setattr(module, '_prepare', lambda *_: (data, data, {}, {}, {}))
    monkeypatch.setattr(module, '_infer', lambda *_a, **_k: pytest.fail('must not call inference'))
    path = tmp_path / 'pairs.csv'
    pd.DataFrame({'sample_id': ['A'], 'cell_id_a': ['cell_a'], 'cell_id_b': ['cell_b']}).to_csv(path, index=False)
    tool = object.__new__(SpatialFormerTool)
    with pytest.raises(ToolValidationError, match='cross-sample'):
        tool.predict_pairs({'pairs_path': str(path)})


def test_embedding_output_preserves_expression_and_ids(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import tooluniverse.spatialformer_tool as module

    data = ad.AnnData(np.array([[1., 2., 3.], [4., 5., 6.]]), obs=pd.DataFrame(index=['a', 'b']))
    work = data[:, :2].copy()
    output = work.copy()
    output.obsm['X_SpaF'] = np.ones((2, 8))
    monkeypatch.setattr(module, '_prepare', lambda *_: (data.copy(), work, {}, {}, {}))
    monkeypatch.setattr(module, '_infer', lambda *_a, **_k: (output, 'cpu'))
    result = object.__new__(SpatialFormerTool).extract_embeddings({'output_dir': str(tmp_path)})
    saved = ad.read_h5ad(result['output_path'])
    np.testing.assert_array_equal(saved.X, data.X)
    assert list(saved.obs_names) == ['a', 'b']
    assert saved.obsm['X_SpaF'].shape == (2, 8)


def test_single_collate_accepts_exactly_two_genes():
    from tooluniverse.spatialformer_tool import _single_collate
    batch = [([27, 28], [1, 25, 8, 3, 6], [1949]),
             ([29], [1, 25, 8, 3, 6], [1949])]
    result = _single_collate(batch)
    assert result["indices"].shape == (2, 8)
    assert result["sequence_length"].tolist() == [[8, 0], [7, 0]]
    assert result["attention_mask"].sum(dim=1).tolist() == [8, 7]
