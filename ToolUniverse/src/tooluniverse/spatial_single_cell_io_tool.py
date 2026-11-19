"""Import already segmented imaging RNA data as single-cell AnnData references."""
import json
from pathlib import Path
from .exceptions import ToolValidationError
from .spatial_transcriptomics_tool import _SpatialTool, _expression, _labels, _read_h5ad, _result_directory
from .tool_registry import register_tool


def _table(path):
    import pandas as pd
    path = Path(path).expanduser()
    if path.suffix == '.parquet':
        return pd.read_parquet(path)
    return pd.read_csv(path, sep='\t' if '.tsv' in path.name else ',', dtype=str)


def _ids(table, keys):
    if not keys or any(k not in table for k in keys):
        raise ToolValidationError(f'Missing identifier columns: {keys}')
    parts = table[keys]
    if parts.isna().any().any() or parts.astype(str).apply(lambda c: c.str.strip().eq('')).any().any():
        raise ToolValidationError('Cell identifiers contain empty values')
    values = parts.iloc[:, 0].astype(str) if len(keys) == 1 else parts.astype(str).apply(lambda r: json.dumps(r.tolist(), separators=(',', ':')), axis=1)
    if values.duplicated().any():
        raise ToolValidationError('Duplicate cell IDs; specify composite identifiers including FOV')
    return values.to_numpy()


def _save(adata, directory, report, name='cells.h5ad'):
    path = directory / name
    adata.write_h5ad(path, compression='gzip')
    report_path = directory / 'report.json'
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
    return {'output_path': str(path), 'report_path': str(report_path), **report}


@register_tool('SpatialSingleCellIOTool')
class SpatialSingleCellIOTool(_SpatialTool):
    operations = ('import_cells', 'attach_annotations', 'prepare_reference')

    def import_cells(self, request):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy import sparse
        platform = request['platform']
        path = Path(request['input_path']).expanduser()
        fmt = request.get('format', 'auto')
        metadata_path = request.get('metadata_path')
        keys = request.get('cell_id_columns', ['cell_id'])
        if fmt == 'auto':
            if path.is_dir() and platform == 'Xenium':
                directory = path
                path = directory / 'cell_feature_matrix.h5'
                if metadata_path is None:
                    metadata_path = next((str(directory / p) for p in ('cells.csv.gz', 'cells.csv', 'cells.parquet') if (directory / p).exists()), None)
                fmt = '10x_h5'
            else:
                fmt = 'h5ad' if path.suffix == '.h5ad' else '10x_h5' if path.suffix == '.h5' else 'csv'
        controls_removed = 0
        if fmt == 'h5ad':
            adata = _read_h5ad(path)
            layer = request.get('count_layer', 'X')
            adata.X = sparse.csr_matrix(_expression(adata, None if layer == 'X' else layer, counts=True))
        elif fmt == '10x_h5':
            import scanpy as sc
            adata = sc.read_10x_h5(path, gex_only=False)
            adata.var['gene_name'] = adata.var_names.astype(str)
            if not adata.var_names.is_unique:
                if 'gene_ids' not in adata.var or not adata.var['gene_ids'].is_unique:
                    raise ToolValidationError('Input has ambiguous duplicate gene identifiers')
                adata.var_names = adata.var['gene_ids'].astype(str)
        elif fmt == 'csv':
            table = _table(path)
            ids = _ids(table, keys)
            genes = request.get('gene_columns')
            if not genes:
                raise ToolValidationError('CSV imports require explicit gene_columns excluding metadata and control probes')
            if len(genes) != len(set(genes)) or any(g not in table for g in genes):
                raise ToolValidationError('gene_columns must be distinct existing columns')
            try:
                values = table[genes].apply(pd.to_numeric, errors='raise').to_numpy()
            except (ValueError, TypeError) as exc:
                raise ToolValidationError('Gene columns must contain numeric counts') from exc
            metadata_columns = [c for c in table.columns if c not in genes]
            obs = table[metadata_columns].fillna('').copy()
            obs.index = ids
            adata = ad.AnnData(sparse.csr_matrix(values), obs=obs, var=pd.DataFrame(index=genes))
        else:
            raise ToolValidationError('Supported formats are auto, h5ad, 10x_h5, csv')
        if 'feature_types' in adata.var:
            keep = adata.var['feature_types'].astype(str).eq('Gene Expression').to_numpy()
            controls_removed = int((~keep).sum())
            adata = adata[:, keep].copy()
        if not adata.n_obs or not adata.n_vars or not adata.obs_names.is_unique or not adata.var_names.is_unique:
            raise ToolValidationError('Import requires nonempty cells and genes with unique identifiers')
        adata.X = sparse.csr_matrix(_expression(adata, counts=True))
        adata.layers['counts'] = adata.X.copy()
        if metadata_path:
            metadata = _table(metadata_path)
            metadata.index = _ids(metadata, request.get('metadata_id_columns', keys))
            missing = adata.obs_names.difference(metadata.index)
            if len(missing):
                raise ToolValidationError(f'Metadata missing {len(missing)} expression cells; exact IDs required')
            for col in metadata:
                adata.obs[col] = metadata.loc[adata.obs_names, col].fillna('').to_numpy()
        gene_key = request.get('gene_name_key')
        if gene_key:
            if gene_key not in adata.var or adata.var[gene_key].isna().any():
                raise ToolValidationError('Missing gene name column or values')
            adata.var['gene_name'] = adata.var[gene_key].astype(str)
        elif 'gene_name' not in adata.var:
            adata.var['gene_name'] = adata.var_names.astype(str)
        if adata.var['gene_name'].isna().any() or adata.var['gene_name'].astype(str).str.strip().eq('').any():
            raise ToolValidationError('Gene names cannot be empty')
        original = adata.obs_names.astype(str).to_numpy()
        sample_id = request['sample_id']
        adata.obs['source_cell_id'] = original
        adata.obs['sample_id'] = sample_id
        adata.obs_names = [json.dumps([sample_id, cell], separators=(',', ':')) for cell in original]
        xy = request.get('coordinate_columns')
        if xy is None and platform == 'Xenium' and all(k in adata.obs for k in ('x_centroid', 'y_centroid')):
            xy = ['x_centroid', 'y_centroid']
        if xy:
            if len(xy) != 2 or any(k not in adata.obs for k in xy):
                raise ToolValidationError('coordinate_columns must name two metadata columns')
            try:
                adata.obsm['spatial'] = adata.obs[xy].astype(float).to_numpy()
            except (ValueError, TypeError) as exc:
                raise ToolValidationError('Coordinates must be numeric') from exc
        has_coordinates = 'spatial' in adata.obsm
        if has_coordinates:
            coords = np.asarray(adata.obsm['spatial'])
            if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
                raise ToolValidationError('Coordinates must be finite cell-by-2 values')
            if not request.get('coordinate_unit') or not request.get('coordinate_frame'):
                raise ToolValidationError('Coordinates require explicit coordinate_unit and coordinate_frame')
            if request['coordinate_frame'] == 'local':
                if not request.get('fov_key'):
                    raise ToolValidationError('Local coordinates require fov_key; do not join independent fields of view')
                adata.obs['fov_id'] = _labels(adata, request['fov_key']).to_numpy()
        adata.uns['spatial_single_cell'] = {'platform': platform, 'source_path': str(path.resolve()), 'sample_id': sample_id, 'observation_unit': 'cell', 'measurement': 'RNA', 'coordinate_unit': request.get('coordinate_unit', 'unknown'), 'coordinate_frame': request.get('coordinate_frame', 'unknown')}
        warnings = ['Conversion preserves targeted assay coverage; unmeasured genes are not zero-expression genes.']
        if adata.var['gene_name'].duplicated().any():
            warnings.append('Duplicate gene symbols need explicit resolution before tokenization.')
        report = {'n_cells': adata.n_obs, 'n_genes': adata.n_vars, 'platform': platform, 'controls_removed': controls_removed, 'has_coordinates': has_coordinates, 'zero_count_cells': int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum()), 'warnings': warnings}
        with _result_directory(request['output_dir'], 'import_cells') as directory:
            return _save(adata, directory, report)

    def attach_annotations(self, request):
        import pandas as pd
        adata = _read_h5ad(request['input_path'])
        table = _table(request['annotation_path'])
        table.index = _ids(table, request.get('annotation_id_columns', ['cell_id']))
        label_key = request['label_column']
        output_key = request.get('output_key', 'cell_type')
        if label_key not in table or output_key in adata.obs:
            raise ToolValidationError('Missing label column or output label column already exists')
        id_key = request.get('input_id_key', 'obs_names')
        ids = pd.Index((adata.obs_names if id_key == 'obs_names' else _labels(adata, id_key)).astype(str))
        if not ids.is_unique:
            raise ToolValidationError('Chosen input IDs are ambiguous; use standardized cell IDs')
        aligned = table[label_key].reindex(ids)
        missing = aligned.isna() | aligned.fillna('').astype(str).str.strip().eq('')
        if missing.any() and request.get('require_all', True):
            raise ToolValidationError(f'Annotations missing or empty for {int(missing.sum())} cells')
        adata.obs[output_key] = aligned.fillna('unassigned').astype(str).to_numpy()
        adata.obs.loc[missing.to_numpy(), output_key] = 'unassigned'
        report = {'n_cells': adata.n_obs, 'matched_cells': int((~missing).sum()), 'unassigned_cells': int(missing.sum()), 'extra_annotation_ids': int(len(table.index.difference(ids))), 'cell_type_key': output_key}
        with _result_directory(request['output_dir'], 'attach_annotations') as directory:
            pd.DataFrame({'cell_id': adata.obs_names, 'lookup_id': ids, 'matched': (~missing).to_numpy(), 'label': adata.obs[output_key].to_numpy()}).to_csv(directory / 'annotation_audit.csv', index=False)
            result = _save(adata, directory, report)
            result['audit_path'] = str(directory / 'annotation_audit.csv')
            return result

    def prepare_reference(self, request):
        import numpy as np
        from scipy import sparse
        adata = _read_h5ad(request['input_path'])
        layer = request.get('count_layer', 'counts')
        adata.X = sparse.csr_matrix(_expression(adata, None if layer == 'X' else layer, counts=True))
        adata.layers['counts'] = adata.X.copy()
        key = request['cell_type_key']
        labels = _labels(adata, key)
        if labels.str.lower().isin(['unassigned', 'unknown', 'unlabeled']).any():
            raise ToolValidationError('Reference contains unassigned labels; annotate or subset cells first')
        counts = labels.value_counts()
        if (counts < request.get('min_cells_per_type', 1)).any():
            raise ToolValidationError('At least one type has fewer than min_cells_per_type cells')
        if (np.asarray(adata.X.sum(axis=1)).ravel() == 0).any():
            raise ToolValidationError('Reference contains zero-count cells; filter before fitting')
        if request.get('batch_key'):
            _labels(adata, request['batch_key'])
        report = {'n_cells': adata.n_obs, 'n_genes': adata.n_vars, 'cell_type_key': key, 'count_layer': 'counts', 'cells_per_type': {str(k): int(v) for k, v in counts.items()}, 'warnings': ['Format compatibility is not evidence of valid cross-platform reference accuracy.']}
        shared = None
        if request.get('spatial_path'):
            spatial = _read_h5ad(request['spatial_path'])
            shared = adata.var_names.intersection(spatial.var_names, sort=False)
            report.update(shared_gene_count=len(shared), reference_gene_coverage=len(shared) / adata.n_vars)
            if not len(shared):
                raise ToolValidationError('No exact shared gene IDs; run SpatialPreprocess.align_gene_ids first')
            if request.get('subset_shared_genes', False):
                adata = adata[:, shared].copy()
                if (np.asarray(adata.X.sum(axis=1)).ravel() == 0).any():
                    raise ToolValidationError('Shared-gene filtering creates zero-count reference cells')
                report['n_genes'] = adata.n_vars
        elif request.get('subset_shared_genes', False):
            raise ToolValidationError('subset_shared_genes requires spatial_path')
        with _result_directory(request['output_dir'], 'prepare_reference') as directory:
            result = _save(adata, directory, report, 'reference.h5ad')
            result['reference_path'] = result['output_path']
            if shared is not None:
                import pandas as pd
                pd.DataFrame({'gene_id': shared}).to_csv(directory / 'shared_genes.csv', index=False)
                result['shared_genes_path'] = str(directory / 'shared_genes.csv')
            return result
