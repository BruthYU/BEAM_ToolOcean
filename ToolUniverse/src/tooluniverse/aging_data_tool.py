"""Expression preparation and donor-aware pseudobulk for aging research."""
from pathlib import Path
import json
from .aging_common import (AgingTool, read_table, write_table, read_expression,
                           write_h5ad, result_directory, validate_ids, require_columns)
from .exceptions import ToolValidationError
from .tool_registry import register_tool


def _metadata(data, path, id_key):
    meta = read_table(path)
    validate_ids(meta, id_key)
    meta = meta.set_index(id_key)
    missing = data.obs_names.difference(meta.index)
    if len(missing):
        raise ToolValidationError(f'Metadata missing {len(missing)} observations, e.g. {missing[:5].tolist()}')
    overlap = data.obs.columns.intersection(meta.columns)
    for col in overlap:
        if not data.obs[col].astype(str).equals(meta.loc[data.obs_names, col].astype(str).rename_axis(data.obs.index.name)):
            # Compare values independent of index/series metadata.
            if (data.obs[col].astype(str).to_numpy() != meta.loc[data.obs_names, col].astype(str).to_numpy()).any():
                raise ToolValidationError(f'Conflicting existing metadata column {col}')
    for col in meta.columns:
        data.obs[col] = meta.loc[data.obs_names, col].to_numpy()
    return len(meta.index.difference(data.obs_names))


def _matrix_values(matrix):
    from scipy.sparse import issparse
    import numpy as np
    return matrix.toarray() if issparse(matrix) else np.asarray(matrix)


@register_tool('AgingIOTool')
class AgingIOTool(AgingTool):
    operations = ('inspect_dataset', 'import_expression', 'attach_metadata', 'export_results')

    def inspect_dataset(self, r):
        import numpy as np
        from scipy.sparse import issparse
        data, x = read_expression(r['input_path'], r.get('layer'))
        values = x.data if issparse(x) else np.asarray(x)
        return {'shape': list(data.shape), 'obs_columns': list(data.obs.columns),
                'var_columns': list(data.var.columns), 'layers': list(data.layers),
                'integer_valued': bool((values == np.floor(values)).all()),
                'recorded_metadata': dict(data.uns.get('aging', {})),
                'zero_total_observations': int((np.asarray(x.sum(axis=1)).ravel() == 0).sum())}

    def import_expression(self, r):
        import anndata as ad
        import pandas as pd
        import numpy as np
        if Path(r['input_path']).suffix == '.h5ad':
            data, _ = read_expression(r['input_path'])
        else:
            frame = read_table(r['input_path'])
            validate_ids(frame, frame.columns[0])
            values = frame.iloc[:, 1:].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float)
            row_ids, col_ids = frame.iloc[:, 0].astype(str).tolist(), frame.columns[1:].tolist()
            if r['orientation'] == 'genes_by_observations':
                values, row_ids, col_ids = values.T, col_ids, row_ids
            if not np.isfinite(values).all() or (values < 0).any():
                raise ToolValidationError('Expression must be finite and nonnegative.')
            data = ad.AnnData(values, obs=pd.DataFrame(index=row_ids), var=pd.DataFrame(index=col_ids))
        if not data.obs_names.is_unique or not data.var_names.is_unique:
            raise ToolValidationError('Observation and gene IDs must be unique.')
        if r['expression_scale'] == 'counts':
            from scipy.sparse import issparse
            values = data.X.data if issparse(data.X) else data.X
            if (values != np.floor(values)).any():
                raise ToolValidationError('counts must be integer valued.')
        if r.get('metadata_path'):
            _metadata(data, r['metadata_path'], r.get('id_key', 'sample_id'))
        data.uns['aging'] = {k:r[k] for k in ('species', 'tissue', 'expression_scale', 'gene_id_type', 'data_kind')}
        with result_directory(r['output_dir'], 'import_expression') as out:
            return {'expression': write_h5ad(data, out)}

    def attach_metadata(self, r):
        data, _ = read_expression(r['input_path'])
        extra = _metadata(data, r['metadata_path'], r['id_key'])
        with result_directory(r['output_dir'], 'attach_metadata') as out:
            return {'expression': write_h5ad(data, out), 'extra_metadata_rows': extra}

    def export_results(self, r):
        import pandas as pd
        data, x = read_expression(r['input_path'], r.get('layer'))
        columns = r.get('obs_columns', list(data.obs.columns))
        require_columns(data.obs, columns)
        frame = data.obs[columns].copy()
        frame.insert(0, 'observation_id', data.obs_names)
        if r.get('genes'):
            missing = set(r['genes']) - set(data.var_names)
            if missing:
                raise ToolValidationError(f'Missing genes: {sorted(missing)}')
            for gene in r['genes']:
                frame['expression:' + gene] = _matrix_values(x[:, data.var_names.get_indexer([gene])]).ravel()
        with result_directory(r['output_dir'], 'export_results') as out:
            return {'table': write_table(frame, out, 'results.csv')}


@register_tool('AgingPreprocessTool')
class AgingPreprocessTool(AgingTool):
    operations = ('filter_qc', 'align_gene_ids', 'normalize_expression', 'prepare_model_input')

    def filter_qc(self, r):
        import numpy as np
        import pandas as pd
        data, x = read_expression(r['input_path'], r.get('layer'), counts=True)
        total = np.asarray(x.sum(axis=1)).ravel()
        genes = np.asarray((x > 0).sum(axis=1)).ravel()
        keep = (total >= r.get('min_counts', 1)) & (genes >= r.get('min_genes', 1))
        if 'max_counts' in r:
            keep &= total <= r['max_counts']
        audit = pd.DataFrame({'observation_id': data.obs_names, 'total_counts':total,
                              'detected_genes':genes, 'retained':keep})
        if not keep.any():
            raise ToolValidationError('QC would remove every observation.')
        filtered = data[keep].copy()
        with result_directory(r['output_dir'], 'filter_qc') as out:
            return {'expression':write_h5ad(filtered,out), 'audit':write_table(audit,out,'qc.csv'),
                    'retained':int(keep.sum()),'excluded':int((~keep).sum())}

    def align_gene_ids(self, r):
        import pandas as pd
        data, _ = read_expression(r['input_path'])
        source = data.var_names.astype(str)
        if r.get('strip_ensembl_version', False):
            source = source.str.replace(r'\.\d+$', '', regex=True)
        if not source.is_unique:
            raise ToolValidationError('Stripping gene versions creates duplicate IDs.')
        mapped = source.copy()
        if r.get('mapping_path'):
            mapping = read_table(r['mapping_path'])
            validate_ids(mapping, 'source_id'); validate_ids(mapping, 'target_id')
            lookup = mapping.set_index('source_id')['target_id']
            mapped = pd.Index([lookup.get(g, g) for g in source])
        if not mapped.is_unique:
            raise ToolValidationError('Gene mapping creates duplicate target IDs.')
        audit = pd.DataFrame({'original_id':data.var_names, 'aligned_id':mapped})
        data.var['original_gene_id'] = data.var_names
        data.var_names = mapped
        if r.get('mapping_path'):
            data.uns.setdefault('aging', {})['gene_id_type'] = r.get('target_gene_id_type', 'unknown')
        with result_directory(r['output_dir'], 'align_gene_ids') as out:
            return {'expression':write_h5ad(data,out), 'audit':write_table(audit,out,'gene_mapping.csv')}

    def normalize_expression(self, r):
        import numpy as np
        from scipy.sparse import issparse, diags
        data, x = read_expression(r['input_path'], r.get('layer'), counts=True)
        totals = np.asarray(x.sum(axis=1)).ravel()
        if (totals <= 0).any():
            raise ToolValidationError('Remove zero-count observations before normalization.')
        target = r.get('target_sum', 10000.0)
        if not np.isfinite(target) or target <= 0:
            raise ToolValidationError('target_sum must be positive and finite.')
        name = r.get('output_layer', 'aging_lognorm')
        if name == 'X':
            raise ToolValidationError('Output layer X is reserved for the primary expression matrix.')
        if name in data.layers:
            raise ToolValidationError('Output layer already exists; choose a new layer name.')
        normalized = diags(target / totals) @ x if issparse(x) else x * (target / totals)[:, None]
        if issparse(normalized):
            normalized = normalized.tocsr()
            normalized.data = np.log1p(normalized.data)
        else:
            normalized = np.log1p(normalized)
        data.layers[name] = normalized
        data.uns.setdefault('aging_layer_scales', {})[name] = 'log_normalized'
        with result_directory(r['output_dir'], 'normalize_expression') as out:
            return {'expression':write_h5ad(data,out), 'layer':name, 'target_sum':float(target),
                    'method':'library_size_scale_then_natural_log1p'}

    def prepare_model_input(self, r):
        import numpy as np
        import pandas as pd
        data, x = read_expression(r['input_path'], r.get('layer'), counts=r['expression_scale']=='counts')
        model = r['model']
        # This operation selects/validates source expression; adapters perform model-specific transforms.
        if r['species'] != 'human':
            raise ToolValidationError('Initial model adapters support human inputs only; no implicit ortholog mapping.')
        if model == 'RNAAgeCalc' and r['data_kind'] != 'bulk':
            raise ToolValidationError('RNAAgeCalc requires bulk input; pseudobulk extrapolation is not enabled.')
        if model == 'scAgeClock' and r['data_kind'] != 'single_cell':
            raise ToolValidationError('scAgeClock requires single_cell input.')
        allowed = {'RNAAgeCalc':('counts','fpkm'), 'scAgeClock':('counts',), 'SenCID':('counts',)}
        if r['expression_scale'] not in allowed[model]:
            raise ToolValidationError(f'{model} accepts source scales {allowed[model]}')
        data = data.copy(); data.X = x.copy()
        # Preserve every input gene. Model-specific feature filtering must occur after normalization.
        data.uns['aging'] = {**dict(data.uns.get('aging', {})), **{k:r[k] for k in ('species','expression_scale','data_kind')}}
        data.uns['aging']['prepared_for'] = model
        coverage = None
        if r.get('feature_path'):
            features = read_table(r['feature_path'])
            validate_ids(features, features.columns[0])
            genes = features.iloc[:,0].astype(str)
            if genes.duplicated().any():
                raise ToolValidationError('Feature genes must be unique.')
            present = genes.isin(data.var_names)
            coverage = {'requested':len(genes), 'present':int(present.sum()), 'fraction':float(present.mean())}
            if coverage['fraction'] < r.get('min_gene_coverage', 0.8):
                raise ToolValidationError(f'Insufficient gene coverage: {coverage}')
        if model == 'scAgeClock':
            fields = r.get('categorical_keys', ['assay','cell_type','tissue_general','sex'])
            require_columns(data.obs, fields)
            validate_ids(data.obs, fields, unique=False)
        if model == 'RNAAgeCalc' and r['expression_scale']=='counts':
            key = r.get('gene_length_key', 'gene_length')
            require_columns(data.var,[key])
            length = pd.to_numeric(data.var[key],errors='coerce').to_numpy()
            if not np.isfinite(length).all() or (length<=0).any():
                raise ToolValidationError('Counts to FPKM requires positive finite gene lengths in bases.')
        with result_directory(r['output_dir'], 'prepare_model_input') as out:
            return {'expression':write_h5ad(data,out), 'coverage':coverage,
                    'preprocessing':'Source selected and validated; normalization belongs to model adapter.'}


@register_tool('AgingPseudobulkTool')
class AgingPseudobulkTool(AgingTool):
    operations = ('aggregate_counts','summarize_cell_counts')

    def _aggregate(self, r, write_counts):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy.sparse import csr_matrix
        data, x = read_expression(r['input_path'], r.get('layer'), counts=True)
        donor = r.get('donor_key','donor_id'); celltype = r.get('cell_type_key','cell_type')
        sample = r.get('sample_key', donor)
        retain = r.get('metadata_columns', [])
        reserved = {'pseudobulk_id', 'n_cells', 'retained'}
        collisions = reserved.intersection([donor, sample, celltype, *retain])
        if collisions:
            raise ToolValidationError(f'Aggregation fields use reserved output names: {sorted(collisions)}')
        if celltype in {donor, sample}:
            raise ToolValidationError('Cell type key must differ from donor and sample keys.')
        group = list(dict.fromkeys([donor,sample,celltype]))
        validate_ids(data.obs, group, unique=False)
        require_columns(data.obs, retain)
        # A biological sample cannot belong to multiple donors.
        if data.obs.groupby(sample,observed=True)[donor].nunique().gt(1).any():
            raise ToolValidationError('Each sample must map to exactly one donor.')
        grouped = data.obs.groupby(group, sort=True, observed=True, dropna=False)
        records=[]; rows=[]; positions=[]
        for i,(key,frame) in enumerate(grouped):
            key = key if isinstance(key,tuple) else (key,)
            record = dict(zip(group,map(str,key)))
            for col in retain:
                if frame[col].isna().any() or frame[col].nunique()!=1:
                    raise ToolValidationError(f'Metadata {col} must be constant within each aggregation group.')
                record[col]=str(frame[col].iloc[0])
            record.update(pseudobulk_id=f'pb_{i:06d}', n_cells=len(frame), retained=len(frame)>=r.get('min_cells',1))
            records.append(record)
            if record['retained']:
                rows.extend([len(positions)]*len(frame)); positions.append(data.obs_names.get_indexer(frame.index))
        table=pd.DataFrame(records)
        with result_directory(r['output_dir'], 'aggregate_counts' if write_counts else 'summarize_cell_counts') as out:
            result={'table':write_table(table,out,'cell_counts.csv')}
            if write_counts:
                if not positions:
                    raise ToolValidationError('No aggregation group passes min_cells.')
                columns=np.concatenate(positions)
                membership=csr_matrix((np.ones(len(columns)),(rows,columns)),shape=(len(positions),data.n_obs))
                matrix=csr_matrix(membership @ x)
                obs=table.loc[table.retained].drop(columns='retained').set_index('pseudobulk_id')
                aggregated=ad.AnnData(matrix,obs=obs,var=data.var.copy())
                aggregated.uns['aging']={'data_kind':'pseudobulk','expression_scale':'counts',
                                          'aggregation_keys':group,'source_path':str(Path(r['input_path']).resolve())}
                result['expression']=write_h5ad(aggregated,out,'pseudobulk.h5ad')
                result['metadata']=write_table(obs.reset_index(),out,'metadata.csv')
            return result

    def aggregate_counts(self,r):
        return self._aggregate(r,True)

    def summarize_cell_counts(self,r):
        return self._aggregate(r,False)
