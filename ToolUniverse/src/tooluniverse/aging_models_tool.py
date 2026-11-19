"""Inference-only adapters for published transcriptomic aging models.

RNAAgeCalc uses an exported official elastic-net coefficient bundle. SenCID's
non-denoising backend implements the published CPM/log2/z-score preprocessing
and calls the original sklearn classifiers. Neither adapter retrains a model.
"""
from pathlib import Path
import hashlib
import importlib.util
import json

from .aging_common import AgingTool, read_expression, result_directory, write_table, write_h5ad
from .exceptions import ToolValidationError, ToolUnavailableError
from .tool_registry import register_tool


def _sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _bundle(path, expected):
    root = Path(path).expanduser().resolve()
    manifest = json.loads((root / 'manifest.json').read_text())
    if manifest.get('model') != expected or not manifest.get('source'):
        raise ToolValidationError(f'Expected a sourced {expected} model bundle.')
    for entry in manifest.get('files', {}).values():
        local = (root / entry['path']).resolve()
        if not local.is_relative_to(root) or not local.is_file() or _sha(local) != entry['sha256']:
            raise ToolValidationError('Model bundle path or checksum mismatch.')
    return root, manifest


def _file(root, manifest, key):
    if key not in manifest.get('files', {}):
        raise ToolValidationError(f'Model bundle missing {key}.')
    return root / manifest['files'][key]['path']


def _human(request):
    if request.get('species') != 'human':
        raise ToolValidationError('These bundled models support species=human only; no implicit ortholog conversion.')


def _dense(matrix):
    import numpy as np
    from scipy.sparse import issparse
    return matrix.toarray() if issparse(matrix) else np.asarray(matrix)


def _coverage(genes, required, threshold):
    found = set(map(str, genes))
    missing = [g for g in required if g not in found]
    fraction = 1 - len(missing) / len(required) if required else 0
    if fraction < threshold:
        raise ToolValidationError(f'Model gene coverage {fraction:.3f} is below {threshold:.3f}; missing {len(missing)} genes.')
    return {'required_genes': len(required), 'matched_genes': len(required) - len(missing),
            'coverage_fraction': fraction, 'missing_genes': missing}


def _metadata(data, id_name, request):
    import pandas as pd
    out = pd.DataFrame({id_name: data.obs_names.astype(str)})
    for key in request.get('metadata_columns', []):
        if key not in data.obs or key in {id_name, 'predicted_age', 'recommended_sid', *[f'SID{i}_{suffix}' for i in range(1, 7) for suffix in ('score', 'decision')]}:
            raise ToolValidationError(f'Invalid metadata column: {key}')
        out[key] = data.obs[key].to_numpy()
    return out


@register_tool('RNAAgeCalcTool')
class RNAAgeCalcTool(AgingTool):
    operations = ('predict_age',)

    def predict_age(self, request):
        import numpy as np
        import pandas as pd
        _human(request)
        root, manifest = _bundle(request['model_dir'], 'RNAAgeCalc')
        if request.get('tissue') != manifest.get('tissue'):
            raise ToolValidationError('Requested tissue must match the fixed model bundle tissue.')
        if request.get('data_kind') != 'bulk':
            raise ToolValidationError('RNAAgeCalc requires data_kind=bulk; pseudobulk is not implicitly validated.')
        expression_type = request['expression_type']
        if expression_type not in ('FPKM', 'counts'):
            raise ToolValidationError('RNAAgeCalc requires FPKM or counts plus gene lengths.')
        if request.get('gene_id_type') != 'SYMBOL':
            raise ToolValidationError('Align genes to human SYMBOL IDs before RNAAgeCalc.')
        data, matrix = read_expression(request['input_path'], request.get('layer'), counts=expression_type == 'counts')
        coefficients = pd.read_csv(_file(root, manifest, 'coefficients'), dtype={'gene': str})
        if not {'gene', 'coefficient'}.issubset(coefficients) or coefficients.gene.duplicated().any():
            raise ToolValidationError('Coefficients must contain unique gene and coefficient columns.')
        coefficients['coefficient'] = pd.to_numeric(coefficients.coefficient, errors='raise')
        if not np.isfinite(coefficients.coefficient).all() or not np.isfinite(manifest.get('intercept', float('nan'))):
            raise ToolValidationError('Coefficient bundle has non-finite weights/intercept.')
        coverage = _coverage(data.var_names, coefficients.gene.tolist(), 1.0)
        # Strict full coverage avoids silently replacing upstream KNN imputation.
        positions = data.var_names.get_indexer(coefficients.gene)
        values = _dense(matrix[:, positions]).astype(float)
        if expression_type == 'counts':
            key = request.get('gene_length_key')
            if not key or key not in data.var:
                raise ToolValidationError('Counts require a positive gene_length_key in h5ad var (base pairs).')
            lengths = pd.to_numeric(data.var[key], errors='raise').to_numpy()
            totals = np.asarray(matrix.sum(axis=1)).ravel()
            if not np.isfinite(lengths).all() or (lengths <= 0).any() or (totals <= 0).any():
                raise ToolValidationError('Gene lengths and library sizes must be positive and finite.')
            values = values / totals[:, None] / lengths[positions][None, :] * 1e9
        predicted = np.log2(values + 1) @ coefficients.coefficient.to_numpy() + float(manifest['intercept'])
        if not np.isfinite(predicted).all():
            raise ToolValidationError('Non-finite predicted ages.')
        table = _metadata(data, 'sample_id', request)
        table['predicted_age'] = predicted
        with result_directory(request['output_dir'], 'rnaagecalc') as out:
            return {'predictions': write_table(table, out, 'predictions.csv'), 'coverage': coverage,
                    'backend': 'official_coefficients_numpy_no_imputation',
                    'model': manifest, 'device': 'cpu', 'age_unit': 'years'}


@register_tool('ScAgeClockTool')
class ScAgeClockTool(AgingTool):
    operations = ('prepare_input', 'predict_age')

    def prepare_input(self, request):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy.sparse import csr_matrix, hstack
        _human(request)
        if request.get('expression_type') != 'counts':
            raise ToolValidationError('scAgeClock preparation requires raw counts.')
        root, manifest = _bundle(request['model_dir'], 'scAgeClock')
        data, matrix = read_expression(request['input_path'], request.get('layer'), counts=True)
        features = json.loads(_file(root, manifest, 'features').read_text())
        cats = manifest['architecture']['cat_cardinalities']
        genes = features[len(cats):]
        coverage = _coverage(data.var_names, genes, request.get('min_gene_coverage', 0.8))
        # Mapping values are explicit; no unknown label is silently encoded as zero.
        names = features[:len(cats)]
        mappings = json.loads(_file(root, manifest, 'categories').read_text())
        metadata_keys = request.get('category_columns', {name: name for name in names})
        use_indices = request.get('categories_are_indices', False)
        cat_values = []
        for name, size in zip(names, cats):
            column = metadata_keys.get(name)
            if not column or column not in data.obs:
                raise ToolValidationError(f'Missing categorical metadata for {name}.')
            raw = data.obs[column]
            if raw.isna().any():
                raise ToolValidationError(f'Missing values in categorical metadata {column}.')
            values = pd.to_numeric(raw, errors='coerce') if use_indices else raw.astype(str).map(mappings[name])
            values = values.to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values != np.floor(values)).any() or (values < 0).any() or (values >= size).any():
                raise ToolValidationError(f'Unmapped or invalid categories in {column}.')
            cat_values.append(values)
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        if (totals <= 0).any():
            raise ToolValidationError('Cannot normalize zero-count observations.')
        # Published default: median library size across all input genes, then natural log1p.
        target = float(request.get('target_sum', np.median(totals)))
        if not np.isfinite(target) or target <= 0:
            raise ToolValidationError('target_sum must be positive and finite.')
        normalized = csr_matrix(matrix, dtype=float).multiply((target / totals)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        positions = data.var_names.get_indexer(genes)
        blocks = [normalized[:, pos] if pos >= 0 else csr_matrix((data.n_obs, 1)) for pos in positions]
        formatted = ad.AnnData(hstack([csr_matrix(np.column_stack(cat_values)), *blocks], format='csr'),
                              obs=data.obs.copy(), var=pd.DataFrame(index=features))
        formatted.uns['aging_scageclock_preparation'] = {'model_manifest_sha256': _sha(root / 'manifest.json'),
                                                        'normalization_target_sum': target,
                                                        'missing_gene_policy': 'zero_fill', 'expression_type': 'scageclock_formatted'}
        with result_directory(request['output_dir'], 'prepare_scageclock') as out:
            return {'output': write_h5ad(formatted, out), 'coverage': coverage,
                    'expression_type': 'scageclock_formatted', 'normalization_target_sum': target}


    def predict_age(self, request):
        import numpy as np
        import torch
        _human(request)
        root, manifest = _bundle(request['model_dir'], 'scAgeClock')
        data, matrix = read_expression(request['input_path'], request.get('layer'))
        if request.get('expression_type') != 'scageclock_formatted':
            raise ToolValidationError('scAgeClock requires its exact formatted features and categorical indices.')
        features = json.loads(_file(root, manifest, 'features').read_text())
        prepared = data.uns.get('aging_scageclock_preparation', {})
        if prepared and prepared.get('model_manifest_sha256') != _sha(root / 'manifest.json'):
            raise ToolValidationError('Input was prepared for a different scAgeClock bundle.')
        if list(data.var_names) != features:
            raise ToolValidationError('scAgeClock feature IDs/order do not match the checkpoint bundle.')
        config = manifest['architecture']
        cats = config['cat_cardinalities']
        if len(features) != len(cats) + config['num_numeric_features']:
            raise ToolValidationError('Feature manifest disagrees with the model architecture.')
        category_values = _dense(matrix[:, :len(cats)])
        for i, size in enumerate(cats):
            v = category_values[:, i]
            if (v != np.floor(v)).any() or (v < 0).any() or (v >= size).any():
                raise ToolValidationError(f'Categorical feature {features[i]} contains invalid indices.')
        device = request.get('device', 'cpu')
        if device not in ('cpu', 'cuda'):
            raise ToolValidationError('device must be cpu or cuda.')
        if device == 'cuda' and not torch.cuda.is_available():
            raise ToolUnavailableError('CUDA requested but unavailable.')
        # Import only the published inference architecture, not training-only dependencies.
        architecture_path = _file(root, manifest, 'architecture_source')
        spec = importlib.util.spec_from_file_location('_aging_scage_gma', architecture_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        model = module.GatedMultiheadAttentionFCNet(**config)
        state = torch.load(_file(root, manifest, 'weights'), map_location='cpu', weights_only=True)
        model.load_state_dict(state)
        model.to(device).eval()
        batch_size = request.get('batch_size', 128)
        if not isinstance(batch_size, int) or batch_size < 1:
            raise ToolValidationError('batch_size must be positive.')
        batches = []
        with torch.inference_mode():
            for start in range(0, data.n_obs, batch_size):
                x = torch.as_tensor(_dense(matrix[start:start + batch_size]), dtype=torch.float32, device=device)
                batches.append(model(x).flatten().cpu().numpy())
        predicted = np.concatenate(batches)
        if len(predicted) != data.n_obs or not np.isfinite(predicted).all():
            raise ToolValidationError('scAgeClock produced invalid predictions.')
        table = _metadata(data, 'cell_id', request)
        table['predicted_age'] = predicted
        with result_directory(request['output_dir'], 'scageclock') as out:
            return {'predictions': write_table(table, out, 'predictions.csv'),
                    'backend': 'official_gma_batched', 'model': manifest,
                    'device': device, 'age_unit': 'years', 'n_cells': data.n_obs}


@register_tool('SenCIDTool')
class SenCIDTool(AgingTool):
    operations = ('score_senescence',)

    def score_senescence(self, request):
        import numpy as np
        import pandas as pd
        import joblib
        from sklearn.preprocessing import StandardScaler
        _human(request)
        if request.get('expression_type') != 'counts' or request.get('gene_id_type') != 'SYMBOL':
            raise ToolValidationError('SenCID requires raw counts with human SYMBOL gene IDs.')
        root, manifest = _bundle(request['model_dir'], 'SenCID')
        data, matrix = read_expression(request['input_path'], request.get('layer'), counts=True)
        features = json.loads(_file(root, manifest, 'features').read_text())
        coverage = _coverage(data.var_names, features, request.get('min_gene_coverage', 0.8))
        sidnums = request.get('sidnums', [1, 2, 3, 4, 5, 6])
        if not sidnums or len(set(sidnums)) != len(sidnums) or any(s not in range(1, 7) for s in sidnums):
            raise ToolValidationError('sidnums must be distinct integers from 1 through 6.')
        # Match upstream ScaleData(denoising=False): retain expressed seneset genes,
        # add absent/all-zero genes with 0.01, CPM over the seneset, then z-score
        # log2(CPM+1) across genes independently for each observation.
        shared = [g for g in features if g in data.var_names]
        selected = _dense(matrix[:, data.var_names.get_indexer(shared)]).astype(float)
        if (selected.sum(axis=1) <= 0).any():
            raise ToolValidationError('SenCID requires measured seneset counts in every observation.')
        expressed = selected.sum(axis=0) >= 1
        if len(shared) <= 10:
            raise ToolValidationError('SenCID requires more than 10 overlapping seneset genes.')
        values = np.full((data.n_obs, len(features)), 0.01, dtype=float)
        feature_index = {g: i for i, g in enumerate(features)}
        for j, gene in enumerate(shared):
            if expressed[j]:
                values[:, feature_index[gene]] = selected[:, j]
        totals = values.sum(axis=1)
        if (totals <= 0).any():
            raise ToolValidationError('SenCID seneset totals must be positive.')
        scaled = StandardScaler().fit_transform(np.log2(values / totals[:, None] * 1e6 + 1).T).T
        all_features = pd.DataFrame(scaled, index=data.obs_names, columns=features)
        table = _metadata(data, 'cell_id' if request.get('data_kind', 'single_cell') == 'single_cell' else 'sample_id', request)
        sid_features = json.loads(_file(root, manifest, 'sid_features').read_text())
        for sid in sidnums:
            x = all_features.loc[:, sid_features[str(sid)]]
            classifier = joblib.load(_file(root, manifest, f'SID{sid}'))
            probability = joblib.load(_file(root, manifest, f'SID{sid}_L'))
            table[f'SID{sid}_score'] = probability.predict_proba(x)[:, 1]
            table[f'SID{sid}_decision'] = classifier.decision_function(x)
        recommender = joblib.load(_file(root, manifest, 'recommend_model'))
        table['recommended_sid'] = recommender.predict(all_features)
        scores = table[[f'SID{s}_score' for s in sidnums]].to_numpy()
        if not np.isfinite(scores).all() or (scores < 0).any() or (scores > 1).any():
            raise ToolValidationError('SenCID produced invalid scores.')
        with result_directory(request['output_dir'], 'sencid') as out:
            return {'scores': write_table(table, out, 'scores.csv'), 'coverage': coverage,
                    'backend': 'official_classifiers_no_dca', 'denoising': False,
                    'model': manifest, 'device': 'cpu',
                    'interpretation': 'SID scores are senescence readouts, not ages or universal binary labels.'}
