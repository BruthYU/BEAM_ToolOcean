"""Coordinate benchmarks and patient-level summaries of predicted spatial proximity."""
import json
from pathlib import Path
from .exceptions import ToolValidationError
from .spatial_transcriptomics_tool import _SpatialTool, _read_h5ad, _labels, _result_directory
from .tool_registry import register_tool

KEYS = ['sample_id', 'cell_id_a', 'cell_id_b']


def _table(path):
    import pandas as pd
    path = Path(path).expanduser()
    try:
        return pd.read_parquet(path) if path.suffix == '.parquet' else pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, ValueError) as exc:
        raise ToolValidationError(f'Cannot read {path}: {exc}') from exc


def _numeric(frame, key, binary=False):
    import numpy as np
    import pandas as pd
    frame[key] = pd.to_numeric(frame[key], errors='coerce')
    if not np.isfinite(frame[key]).all() or not frame[key].between(0, 1).all():
        raise ToolValidationError(f'{key} must be finite in [0, 1]')
    if binary and not frame[key].isin([0, 1]).all():
        raise ToolValidationError(f'{key} must be binary')


def _identifiers(frame, keys):
    if frame.empty or not set(keys).issubset(frame):
        raise ToolValidationError(f'Expected a nonempty table with columns {keys}')
    for key in keys:
        if frame[key].isna().any() or frame[key].astype(str).str.strip().eq('').any():
            raise ToolValidationError(f'Missing identifier in {key}')
        frame[key] = frame[key].astype(str)


def _pairs(path, value):
    frame = _table(path)
    _identifiers(frame, KEYS)
    if value not in frame:
        raise ToolValidationError(f'Missing {value}')
    if frame.duplicated(KEYS).any() or (frame.cell_id_a == frame.cell_id_b).any():
        raise ToolValidationError('Duplicate ordered pairs or self pairs are not supported')
    _numeric(frame, value, value == 'is_neighbor')
    return frame


def _join(prediction, benchmark):
    pred, truth = _pairs(prediction, 'proximity_score'), _pairs(benchmark, 'is_neighbor')
    extra = [c for c in ('distance', 'frame_id') if c in truth]
    if 'distance' in truth:
        import numpy as np
        import pandas as pd
        truth['distance'] = pd.to_numeric(truth['distance'], errors='coerce')
        if not np.isfinite(truth.distance).all() or truth.distance.lt(0).any():
            raise ToolValidationError('Benchmark distance must be finite and nonnegative')
    if 'frame_id' in pred and 'frame_id' in truth:
        frames = pred[KEYS + ['frame_id']].merge(truth[KEYS + ['frame_id']], on=KEYS, suffixes=('_pred', '_truth'))
        if frames.frame_id_pred.astype(str).ne(frames.frame_id_truth.astype(str)).any():
            raise ToolValidationError('Prediction and benchmark coordinate frames disagree')
    result = pred[KEYS + ['proximity_score']].merge(truth[KEYS + ['is_neighbor'] + extra], on=KEYS, how='outer', validate='one_to_one', indicator=True)
    if result._merge.ne('both').any():
        raise ToolValidationError('Prediction and benchmark ordered cell-pair IDs must match exactly')
    return result.drop(columns='_merge')


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def _metrics(frame, threshold=None):
    from sklearn.metrics import roc_auc_score, average_precision_score, confusion_matrix, precision_score, recall_score
    y, score = frame.is_neighbor.astype(int), frame.proximity_score.astype(float)
    both = y.nunique() == 2
    result = dict(n_pairs=len(frame), positive_fraction=float(y.mean()),
                  auroc=float(roc_auc_score(y, score)) if both else None,
                  auprc=float(average_precision_score(y, score)) if both else None,
                  status='ok' if both else 'single_class')
    if threshold is not None:
        pred = score >= threshold
        result.update(threshold=threshold, confusion_matrix=confusion_matrix(y, pred, labels=[0, 1]).tolist(),
                      precision=float(precision_score(y, pred, zero_division=0)), recall=float(recall_score(y, pred, zero_division=0)))
    return result


@register_tool('SpatialPairValidationTool')
class SpatialPairValidationTool(_SpatialTool):
    operations = ('build_pair_benchmark', 'evaluate_pairs')

    def build_pair_benchmark(self, request):
        import numpy as np
        import pandas as pd
        from scipy.spatial import cKDTree
        adata = _read_h5ad(request['input_path'])
        sample, frame = (_labels(adata, request[k]) for k in ('sample_key', 'frame_key'))
        metadata = adata.uns.get('spatial_single_cell', {})
        unit = metadata.get('coordinate_unit')
        if unit and unit != 'unknown' and unit != request['coordinate_unit']:
            raise ToolValidationError('coordinate_unit disagrees with recorded coordinate units')
        if metadata.get('coordinate_frame') == 'local':
            if request['frame_key'] == request['sample_key']:
                raise ToolValidationError('Local coordinates require a distinct field-of-view frame key')
            if 'fov_id' not in adata.obs:
                raise ToolValidationError('Local coordinates require recorded fov_id metadata')
            checks = pd.DataFrame({'sample': sample.to_numpy(), 'frame': frame.to_numpy(), 'fov': _labels(adata, 'fov_id').to_numpy()})
            if checks.groupby(['sample', 'frame']).fov.nunique().gt(1).any():
                raise ToolValidationError('A local coordinate frame cannot combine distinct fov_id values')
        key = request.get('coordinate_key', 'spatial')
        if key not in adata.obsm:
            raise ToolValidationError(f'Missing coordinates {key}')
        xy = np.asarray(adata.obsm[key], dtype=float)
        if xy.ndim != 2 or xy.shape[1] not in (2, 3) or not np.isfinite(xy).all():
            raise ToolValidationError('Expected finite 2D or 3D coordinates')
        radius, target = request['radius'], request.get('pairs_per_class', 1000)
        rng = np.random.default_rng(request.get('random_seed', 0))
        rows, audits = [], []
        groups = pd.DataFrame({'sample': sample.to_numpy(), 'frame': frame.to_numpy()})
        for (sid, fid), idx in groups.groupby(['sample', 'frame'], sort=True).indices.items():
            coords, ids = xy[idx], adata.obs_names[idx].astype(str).to_numpy()
            n = len(ids)
            if n < 2:
                audits.append(dict(sample_id=sid, frame_id=fid, n_cells=n, status='too_few_cells'))
                continue
            tree = cKDTree(coords)
            degree = np.asarray(tree.query_ball_point(coords, radius, return_length=True)) - 1
            positive, total = int(degree.sum() // 2), n * (n - 1) // 2
            probabilities = degree / degree.sum() if degree.sum() else None
            selected = {0: set(), 1: set()}
            for label, available in ((1, positive), (0, total - positive)):
                goal = min(target, available)
                for _ in range(max(1000, goal * 100)):
                    if len(selected[label]) >= goal:
                        break
                    if label:
                        a = int(rng.choice(n, p=probabilities))
                        neighbors = [b for b in tree.query_ball_point(coords[a], radius) if b != a]
                        b = int(rng.choice(neighbors))
                    else:
                        a, b = rng.choice(n, 2, replace=False).tolist()
                        if np.linalg.norm(coords[a] - coords[b]) <= radius:
                            continue
                    selected[label].add(tuple(sorted((a, b))))
                for a, b in sorted(selected[label]):
                    if ids[a] > ids[b]:
                        a, b = b, a
                    rows.append(dict(sample_id=sid, frame_id=fid, cell_id_a=ids[a], cell_id_b=ids[b],
                                     distance=float(np.linalg.norm(coords[a] - coords[b])), is_neighbor=label))
            complete = len(selected[1]) == min(target, positive) and len(selected[0]) == min(target, total-positive)
            audits.append(dict(sample_id=sid, frame_id=fid, n_cells=n, population_pairs=total,
                               population_positive_pairs=positive, sampled_positive=len(selected[1]), sampled_negative=len(selected[0]),
                               status='ok' if complete else 'sampling_budget_exhausted'))
        if not rows:
            raise ToolValidationError('No valid within-frame pairs')
        result = pd.DataFrame(rows)
        with _result_directory(request['output_dir'], 'pair_benchmark') as out:
            result.to_csv(out / 'benchmark.csv', index=False)
            result[KEYS].to_csv(out / 'pairs.csv', index=False)
            _json(out / 'report.json', dict(radius=radius, coordinate_unit=request['coordinate_unit'],
                  sample_key=request['sample_key'], frame_key=request['frame_key'], random_seed=request.get('random_seed', 0),
                  sampling='uniform_without_replacement_within_class_and_frame_bounded_attempts', frames=audits,
                  limitation='Balanced sampling changes prevalence. Pretraining overlap is not assessed.'))
            return dict(benchmark_path=str(out/'benchmark.csv'), pairs_path=str(out/'pairs.csv'), report_path=str(out/'report.json'), n_pairs=len(result))

    def evaluate_pairs(self, request):
        import pandas as pd
        frame = _join(request['prediction_path'], request['benchmark_path'])
        threshold = request.get('threshold')
        overall = _metrics(frame, threshold)
        samples = [dict(sample_id=sid, **_metrics(part, threshold)) for sid, part in frame.groupby('sample_id', sort=True)]
        with _result_directory(request['output_dir'], 'pair_evaluation') as out:
            frame.to_csv(out/'comparison.csv', index=False)
            pd.DataFrame(samples).to_csv(out/'sample_metrics.csv', index=False)
            distance_path = None
            if 'distance' in frame:
                frame['distance_interval'] = pd.cut(frame.distance.astype(float), bins=5, include_lowest=True)
                distance_summary = frame.groupby('distance_interval', observed=True).agg(n_pairs=('proximity_score','size'), mean_score=('proximity_score','mean'), observed_neighbor_fraction=('is_neighbor','mean')).reset_index()
                distance_path = str(out/'distance_summary.csv')
                distance_summary.to_csv(distance_path, index=False)
            _json(out/'metrics.json', dict(overall=overall, samples=samples, pretraining_overlap=request.get('pretraining_overlap', 'unknown'),
                  limitation='Shared cells induce dependent pairs. Metrics describe supplied sampling, not natural prevalence.'))
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            fig = Figure(figsize=(6, 4), layout='constrained')
            FigureCanvasAgg(fig)
            ax = fig.subplots()
            for label, title in ((0, 'Distant'), (1, 'Neighbor')):
                values = frame.loc[frame.is_neighbor == label, 'proximity_score']
                if len(values):
                    ax.hist(values, bins=20, range=(0, 1), alpha=.5, label=title)
            ax.set(xlabel='Predicted proximity score', ylabel='Sampled pair count')
            ax.legend()
            fig.savefig(out/'score_distribution.png', dpi=120)
            return dict(metrics=overall, metrics_path=str(out/'metrics.json'), comparison_path=str(out/'comparison.csv'),
                        sample_metrics_path=str(out/'sample_metrics.csv'), distance_summary_path=distance_path, figure_path=str(out/'score_distribution.png'))


@register_tool('SpatialPairSummaryTool')
class SpatialPairSummaryTool(_SpatialTool):
    operations = ('summarize_type_pairs', 'compare_groups')

    def summarize_type_pairs(self, request):
        import pandas as pd
        frame = _join(request['prediction_path'], request['benchmark_path']) if request.get('benchmark_path') else _pairs(request['prediction_path'], 'proximity_score')
        adata = _read_h5ad(request['input_path'])
        labels, samples = (_labels(adata, request[k]) for k in ('cell_type_key', 'sample_key'))
        frames = _labels(adata, request['frame_key'])
        if adata.uns.get('spatial_single_cell', {}).get('coordinate_frame') == 'local':
            if request['frame_key'] == request['sample_key'] or 'fov_id' not in adata.obs:
                raise ToolValidationError('Local coordinates require independent field-of-view frames')
            checks = pd.DataFrame({'sample': samples.to_numpy(), 'frame': frames.to_numpy(), 'fov': _labels(adata, 'fov_id').to_numpy()})
            if checks.groupby(['sample', 'frame']).fov.nunique().gt(1).any():
                raise ToolValidationError('Local frames cannot combine distinct fov_id values')
        lookup = pd.DataFrame(dict(sample_id=samples.to_numpy(), cell_id=adata.obs_names.astype(str), cell_type=labels.to_numpy(), annotation_frame=frames.to_numpy()))
        for side in ('a', 'b'):
            frame = frame.merge(lookup.rename(columns={'cell_id': f'cell_id_{side}', 'cell_type': f'cell_type_{side}', 'annotation_frame': f'annotation_frame_{side}'}),
                                on=['sample_id', f'cell_id_{side}'], how='left', validate='many_to_one')
            if frame[f'cell_type_{side}'].isna().any():
                raise ToolValidationError('Pair IDs and samples must exactly match annotated cells')
        if frame.annotation_frame_a.ne(frame.annotation_frame_b).any():
            raise ToolValidationError('Pairs must stay within a common coordinate frame')
        if 'frame_id' in frame and frame.frame_id.astype(str).ne(frame.annotation_frame_a).any():
            raise ToolValidationError('Benchmark and annotation frames disagree')
        pairs = [tuple(sorted(p)) for p in zip(frame.cell_id_a, frame.cell_id_b)]
        types = [tuple(sorted(p)) for p in zip(frame.cell_type_a, frame.cell_type_b)]
        for i, side in enumerate(('a', 'b')):
            frame[f'canonical_{side}'] = [p[i] for p in pairs]
            frame[f'type_{side}'] = [p[i] for p in types]
        agg = {'proximity_score': 'mean'}
        if 'is_neighbor' in frame:
            agg['is_neighbor'] = 'mean'
        if 'is_neighbor' in frame and frame.groupby(['sample_id', 'canonical_a', 'canonical_b']).is_neighbor.nunique().gt(1).any():
            raise ToolValidationError('Reciprocal pairs have inconsistent benchmark labels')
        unique = frame.groupby(['sample_id', 'canonical_a', 'canonical_b', 'type_a', 'type_b'], as_index=False).agg(agg)
        grouping = ['sample_id', 'type_a', 'type_b']
        summary = unique.groupby(grouping, as_index=False).agg(n_pairs=('proximity_score', 'size'), mean_score=('proximity_score', 'mean'))
        if 'is_neighbor' in unique:
            observed = unique.groupby(grouping, as_index=False).agg(observed_neighbor_fraction=('is_neighbor', 'mean'))
            summary = summary.merge(observed, on=grouping, validate='one_to_one')
        with _result_directory(request['output_dir'], 'pair_summary') as out:
            summary.to_csv(out/'type_pairs.csv', index=False)
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            import numpy as np
            averaged = summary.groupby(['type_a','type_b'], as_index=False).mean(numeric_only=True)
            type_names = sorted(set(averaged.type_a) | set(averaged.type_b))
            # Limit only visualization size; complete statistics are retained in CSV.
            shown = type_names[:30]
            index = {name:i for i,name in enumerate(shown)}
            matrix = np.full((len(shown),len(shown)),np.nan)
            for row in averaged.itertuples():
                if row.type_a in index and row.type_b in index:
                    i,j=index[row.type_a],index[row.type_b]
                    matrix[i,j]=matrix[j,i]=row.mean_score
            figure = Figure(figsize=(max(5,len(shown)*.3),max(4,len(shown)*.3)),layout='constrained')
            FigureCanvasAgg(figure)
            axis=figure.subplots()
            im=axis.imshow(np.ma.masked_invalid(matrix),vmin=0,vmax=1,cmap='viridis')
            axis.set(xticks=range(len(shown)),yticks=range(len(shown)),xticklabels=shown,yticklabels=shown,title='Mean proximity score (equal sample weights)')
            axis.tick_params(axis='x',labelrotation=90)
            figure.colorbar(im,ax=axis)
            figure.savefig(out/'type_pair_heatmap.png',dpi=120)
            network=Figure(figsize=(7,7),layout='constrained')
            FigureCanvasAgg(network)
            axis=network.subplots()
            angles=np.linspace(0,2*np.pi,len(shown),endpoint=False)
            positions=np.column_stack([np.cos(angles),np.sin(angles)])
            for i in range(len(shown)):
                for j in range(i+1,len(shown)):
                    if np.isfinite(matrix[i,j]):
                        axis.plot(positions[[i,j],0],positions[[i,j],1],color='#34678c',alpha=.15+.7*matrix[i,j],linewidth=.2+3*matrix[i,j])
            axis.scatter(positions[:,0],positions[:,1],s=90,color='#34678c')
            for name,(x,y) in zip(shown,positions):
                axis.text(x*1.13,y*1.13,name,ha='center',va='center',fontsize=8)
            axis.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5),title='Predicted spatial proximity, not communication')
            axis.set_aspect('equal');axis.axis('off')
            network.savefig(out/'type_pair_network.png',dpi=120)
            _json(out/'report.json', dict(sampling_description=request['sampling_description'], plotted_types=shown, omitted_plot_types=type_names[30:], reciprocal_policy='mean_once_per_unordered_pair',
                  limitation='Statistics describe sampled pairs, not population enrichment, abundance, or causal communication.'))
            return dict(summary_path=str(out/'type_pairs.csv'), report_path=str(out/'report.json'), heatmap_path=str(out/'type_pair_heatmap.png'), network_path=str(out/'type_pair_network.png'), n_type_pairs=len(summary))

    def compare_groups(self, request):
        import numpy as np
        import pandas as pd
        from .spatial_association_tool import _patient_comparison
        summary, meta = _table(request['summary_path']), _table(request['metadata_path'])
        metric = request.get('metric', 'mean_score')
        _identifiers(summary, ['sample_id', 'type_a', 'type_b'])
        _identifiers(meta, ['sample_id', 'patient_id', 'group'])
        if metric not in summary:
            raise ToolValidationError(f'Missing summary metric {metric}')
        _numeric(summary, metric)
        if meta.sample_id.duplicated().any() or summary.duplicated(['sample_id', 'type_a', 'type_b']).any():
            raise ToolValidationError('Duplicate sample metadata or sample/type-pair summary')
        if meta.groupby('patient_id')['group'].nunique().gt(1).any():
            raise ToolValidationError('Each patient must belong to exactly one group')
        merged = summary.merge(meta[['sample_id', 'patient_id', 'group']], on='sample_id', how='left', validate='many_to_one')
        if merged.patient_id.isna().any():
            raise ToolValidationError('Every sample requires patient metadata')
        groups = request['compare_groups']
        merged = merged[merged.group.isin(groups)]
        if merged.empty:
            raise ToolValidationError('No samples in requested groups')
        patients = merged.groupby(['patient_id', 'group', 'type_a', 'type_b'], as_index=False)[metric].mean()
        rng = np.random.default_rng(request.get('random_seed', 0))
        rows = []
        for (a, b), part in patients.groupby(['type_a', 'type_b'], sort=True):
            values = [dict(group=row.group, score=float(getattr(row, metric))) for row in part.itertuples()]
            result = _patient_comparison(values, groups, request.get('n_permutations', 9999), rng)
            rows.append(dict(type_a=a, type_b=b, q_value=None, **result))
        valid = [i for i, row in enumerate(rows) if row['p_value'] is not None]
        if valid:
            p = np.array([rows[i]['p_value'] for i in valid])
            order = np.argsort(p)
            factor = sum(1/i for i in range(1, len(p)+1)) if request.get('fdr_method', 'bh') == 'by' else 1
            q = np.minimum.accumulate((p[order]*len(p)*factor/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
            for rank, idx in enumerate(order):
                rows[valid[idx]]['q_value'] = float(q[rank])
        with _result_directory(request['output_dir'], 'pair_group_comparison') as out:
            patients.to_csv(out/'patient_summary.csv', index=False)
            pd.DataFrame(rows).to_csv(out/'group_comparison.csv', index=False)
            _json(out/'report.json', dict(unit='patient', metric=metric, compare_groups=groups,
                  sample_aggregation='equal_weight_mean_within_patient', fdr_method=request.get('fdr_method','bh'),
                  number_tests=len(rows), number_testable=len(valid), results=rows,
                  limitation='Missing type pairs are not zero. Comparable sampling is required; scores are not abundance or communication.'))
            return dict(comparison_path=str(out/'group_comparison.csv'), patient_summary_path=str(out/'patient_summary.csv'),
                        report_path=str(out/'report.json'), number_tests=len(rows), number_testable=len(valid))
