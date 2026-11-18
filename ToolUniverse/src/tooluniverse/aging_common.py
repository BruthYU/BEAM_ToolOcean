"""Shared local I/O, validation and provenance for the Aging category."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import importlib.metadata
import json
import shutil
import tempfile

from .base_tool import BaseTool
from .exceptions import ToolValidationError, ToolUnavailableError


class AgingTool(BaseTool):
    operations = ()

    def run(self, arguments):
        request = arguments.get('request', {})
        operation = request.get('operation')
        if operation not in self.operations:
            raise ToolValidationError(f'Unsupported operation {operation!r}; choose {self.operations}')
        # Validate direct calls too, not only calls made through the SDK.
        if self.tool_config.get('parameter'):
            from jsonschema import Draft202012Validator
            from jsonschema.exceptions import best_match
            problem = best_match(Draft202012Validator(self.tool_config['parameter']).iter_errors(arguments))
            if problem:
                raise ToolValidationError(f'Invalid {operation} request: {problem.message}')
        try:
            result = {'operation': operation, **getattr(self, operation)(request)}
            result['provenance'] = {
                'tool': self.__class__.__name__, 'parameters': request,
                'versions': versions(),
            }
            # Do not emit non-standard JSON NaN or Infinity.
            return clean_json(result)
        except ImportError as exc:
            raise ToolUnavailableError(f'{operation} dependency unavailable: {exc}', retriable=False) from exc


def clean_json(value):
    import numpy as np
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, np.generic):
        return clean_json(value.item())
    if isinstance(value, float) and not __import__('math').isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def versions():
    result = {}
    for name in ('tooluniverse', 'anndata', 'numpy', 'pandas', 'scipy', 'pydeseq2', 'gseapy', 'scikit-learn', 'joblib', 'torch', 'rdata'):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return result


@contextmanager
def result_directory(output_dir, operation):
    parent = Path(output_dir).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix=operation + '_', dir=parent))
    try:
        yield directory
    except BaseException:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def read_table(path):
    import pandas as pd
    path = Path(path).expanduser()
    if path.suffix.lower() not in ('.csv', '.tsv', '.txt'):
        raise ToolValidationError('Use a CSV or TSV table.')
    import csv
    try:
        with path.open(newline='') as handle:
            header = next(csv.reader(handle, delimiter='\t' if path.suffix.lower() in ('.tsv', '.txt') else ','), [])
        if not header or len(set(header)) != len(header) or any(not h.strip() for h in header[1:]):
            raise ToolValidationError('Table requires nonempty, unique column headers.')
        # Read strings first so IDs like 001 are not silently converted to integers.
        frame = pd.read_csv(path, sep='\t' if path.suffix.lower() in ('.tsv', '.txt') else ',', dtype=str, keep_default_na=False)
    except (OSError, ValueError) as exc:
        raise ToolValidationError(f'Cannot read {path}: {exc}') from exc
    if frame.empty or not frame.columns.is_unique:
        raise ToolValidationError('Table must have rows and unique column names.')
    return frame


def write_table(frame, outdir, filename):
    path = Path(outdir) / filename
    frame.to_csv(path, index=False)
    return {'path': str(path.resolve()), 'shape': list(frame.shape), 'columns': list(frame.columns)}


def require_columns(frame, columns):
    missing = [c for c in columns if c not in frame]
    if missing:
        raise ToolValidationError(f'Missing columns: {missing}')


def validate_ids(frame, columns, unique=True):
    if isinstance(columns, str):
        columns = [columns]
    require_columns(frame, columns)
    if frame[columns].isna().any().any() or frame[columns].astype(str).apply(lambda x:x.str.strip().eq('')).any().any():
        raise ToolValidationError(f'Missing identifiers in {columns}')
    if unique and frame.duplicated(columns).any():
        raise ToolValidationError(f'Duplicate identifiers in {columns}')


def read_expression(path, layer=None, counts=False):
    import anndata as ad
    import numpy as np
    from scipy.sparse import issparse
    path = Path(path)
    if path.suffix.lower() == '.h5ad':
        data = ad.read_h5ad(path)
    else:
        import pandas as pd
        frame = read_table(path)
        validate_ids(frame, frame.columns[0])
        data = ad.AnnData(frame.iloc[:, 1:].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float),
                          obs=pd.DataFrame(index=frame.iloc[:, 0].astype(str)),
                          var=pd.DataFrame(index=frame.columns[1:].astype(str)))
    if not data.n_obs or not data.n_vars or not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ToolValidationError('Expression needs nonempty unique observation and gene IDs.')
    if layer not in (None, 'X') and layer not in data.layers:
        raise ToolValidationError(f'Unknown expression layer: {layer}')
    matrix = data.X if layer in (None, 'X') else data.layers[layer]
    if matrix is None:
        raise ToolValidationError('Selected matrix is empty.')
    values = matrix.data if issparse(matrix) else np.asarray(matrix)
    if values.dtype.kind not in 'biuf' or not np.isfinite(values).all() or (values < 0).any():
        raise ToolValidationError('Expression must be finite and nonnegative.')
    if counts and (values != np.floor(values)).any():
        raise ToolValidationError('Raw integer counts are required.')
    return data, matrix


def write_h5ad(data, directory, filename='expression.h5ad'):
    path = Path(directory) / filename
    data.write_h5ad(path)
    return {'path': str(path.resolve()), 'shape': list(data.shape), 'axes': ['observations', 'genes']}
