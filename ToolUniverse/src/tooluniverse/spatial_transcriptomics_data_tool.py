"""File I/O and preprocessing for local spatial transcriptomics workflows."""

from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import tempfile

from .exceptions import ToolValidationError
from .spatial_transcriptomics_tool import (
    _SpatialTool,
    _expression,
    _read_h5ad,
    _result_directory,
)
from .tool_registry import register_tool


def _path(value):
    return Path(value).expanduser().resolve()


@contextmanager
def _open_h5ad(path):
    from anndata import read_h5ad

    try:
        adata = read_h5ad(_path(path), backed="r")
    except (OSError, ValueError, KeyError) as exc:
        raise ToolValidationError(f"Cannot open AnnData file {path!r}: {exc}") from exc
    try:
        yield adata
    finally:
        adata.file.close()


@contextmanager
def _new_file(path):
    """Publish a complete result without overwriting an existing file."""
    path = _path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            raise FileExistsError(f"Output already exists at {path}")
        with tempfile.TemporaryDirectory(dir=path.parent) as directory:
            temporary = Path(directory) / path.name
            yield temporary
            os.link(temporary, path)
    except OSError as exc:
        raise ToolValidationError(f"Cannot write {path}: {exc}") from exc


def _write_h5ad(adata, path):
    with _new_file(path) as temporary:
        adata.write_h5ad(temporary, compression="gzip")
    return str(_path(path))


def _count_matrix(adata, request):
    layer = request.get("count_layer", "X")
    return _expression(adata, None if layer == "X" else layer, counts=True)


def _dense(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else matrix


def _coordinates(adata, key):
    import numpy as np

    values = adata.obsm.get(key)
    if values is None:
        raise ToolValidationError(f"Coordinate field obsm[{key!r}] was not found.")
    try:
        coordinates = np.asarray(_dense(values), dtype=float)
    except (TypeError, ValueError) as exc:
        raise ToolValidationError(f"Coordinate field obsm[{key!r}] must be numeric.") from exc
    if coordinates.shape != (adata.n_obs, 2):
        raise ToolValidationError(f"Coordinate field obsm[{key!r}] must have two columns.")
    return coordinates


def _valid_positions(adata):
    import numpy as np

    coordinates = _coordinates(adata, "spatial")
    if not np.isfinite(coordinates).all():
        raise ToolValidationError("Spatial coordinates must be finite.")
    return coordinates


def _matrix_check(matrix, mode):
    import numpy as np

    if mode == "none":
        return {"status": "not_checked", "rows_checked": 0}
    if matrix is None:
        return {"status": "missing", "rows_checked": 0}
    rows, columns = matrix.shape
    if not rows or not columns:
        return {"status": "empty", "rows_checked": 0}
    if mode == "sample":
        indices = np.linspace(0, rows - 1, min(rows, 128), dtype=int)
        chunks = (matrix[indices, :],)
        checked = len(indices)
    else:
        chunks = (matrix[start:start + 256, :] for start in range(0, rows, 256))
        checked = rows
    from scipy.sparse import issparse

    invalid = 0
    for chunk in chunks:
        values = chunk.data if issparse(chunk) else np.asarray(chunk)
        if values.dtype.kind not in "buif":
            return {"status": "non_numeric", "rows_checked": checked}
        invalid += int((~np.isfinite(values) | (values < 0) | (values != np.floor(values))).sum())
    return {
        "status": "invalid" if invalid else "compatible",
        "rows_checked": checked,
        "invalid_values": invalid,
    }


def _read_visium_positions(directory):
    import pandas as pd

    columns = [
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres",
    ]
    modern = directory / "tissue_positions.csv"
    legacy = directory / "tissue_positions_list.csv"
    if modern.is_file():
        table = pd.read_csv(modern, dtype={"barcode": str})
    elif legacy.is_file():
        table = pd.read_csv(legacy, header=None, names=columns, dtype={"barcode": str})
    else:
        raise ToolValidationError("Expected tissue_positions.csv or tissue_positions_list.csv.")
    if any(column not in table for column in columns):
        raise ToolValidationError("The Visium position table is missing required columns.")
    if table["barcode"].isna().any() or not table["barcode"].is_unique:
        raise ToolValidationError("The Visium position table must contain unique, nonempty barcodes.")
    table = table.loc[:, columns].set_index("barcode")
    for column in columns[1:]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if not table["in_tissue"].isin([0, 1]).all():
        raise ToolValidationError("The Visium in_tissue column must contain 0 or 1.")
    return table


@register_tool("SpatialIOTool")
class SpatialIOTool(_SpatialTool):
    operations = ("import_visium", "inspect_spatial_h5ad", "export_spatial_table")

    def import_visium(self, request):
        import numpy as np
        import scanpy as sc

        source = _path(request["input_dir"])
        kind = request.get("matrix_kind", "filtered")
        prefix = f"{kind}_feature_bc_matrix"
        h5 = source / f"{prefix}.h5"
        mtx = source / prefix
        try:
            if h5.is_file():
                adata = sc.read_10x_h5(h5)
                matrix_path = h5
            elif mtx.is_dir():
                adata = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=False)
                matrix_path = mtx
            else:
                raise ToolValidationError(f"No {kind} 10x count matrix was found in {source}.")
            positions = _read_visium_positions(source / "spatial")
        except (OSError, ValueError, KeyError) as exc:
            raise ToolValidationError(f"Cannot read Space Ranger output: {exc}") from exc

        if not adata.obs_names.is_unique:
            raise ToolValidationError("The count matrix must contain unique barcodes.")
        adata.var["gene_symbol"] = adata.var_names.to_numpy()
        ids = adata.var.get("gene_ids")
        if ids is None or ids.isna().any() or not ids.is_unique or ids.astype(str).str.strip().eq("").any():
            raise ToolValidationError("The count matrix must contain unique, nonempty gene IDs.")
        adata.var_names = ids.astype(str).to_numpy()

        matched = adata.obs_names.intersection(positions.index)
        unmatched = adata.obs_names.difference(positions.index)
        if len(unmatched):
            raise ToolValidationError(
                f"{len(unmatched)} matrix barcodes have no spatial position. Examples: {unmatched[:5].tolist()}."
            )
        tissue_only = request.get("tissue_only", True)
        in_tissue = positions.loc[matched, "in_tissue"].eq(1)
        selected = matched[in_tissue.to_numpy()] if tissue_only else matched
        if not len(selected):
            raise ToolValidationError("No Visium spots remain after applying tissue_only.")
        adata = adata[selected].copy()
        table = positions.loc[selected]
        if not np.isfinite(table[["array_row", "array_col"]]).to_numpy().all():
            raise ToolValidationError("Visium array row and column values must be finite.")
        for key in ("in_tissue", "array_row", "array_col"):
            adata.obs[key] = table[key]
        adata.obsm["spatial"] = table[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        _valid_positions(adata)
        _expression(adata, counts=True)

        default_library = source.parent.name if source.name == "outs" else source.name
        library = request.get("library_id", default_library)
        if not library.strip() or "/" in library or library == ".":
            raise ToolValidationError("library_id must be nonempty and cannot contain a slash or be '.'.")
        adata.obs["library_id"] = library
        levels = request.get("image_levels", ["lowres"])
        if levels:
            from matplotlib.image import imread

            spatial = source / "spatial"
            try:
                scales = json.loads((spatial / "scalefactors_json.json").read_text())
                images = {level: imread(spatial / f"tissue_{level}_image.png") for level in levels}
                factors = [float(scales[f"tissue_{level}_scalef"]) for level in levels]
                if any(not math.isfinite(value) or value <= 0 for value in factors):
                    raise ValueError("Image scale factors must be finite and positive.")
            except (OSError, ValueError, KeyError, TypeError) as exc:
                raise ToolValidationError(f"Cannot read the requested Visium images and scales: {exc}") from exc
            adata.uns["spatial"] = {library: {"images": images, "scalefactors": scales}}

        output = _write_h5ad(adata, request["output_path"])
        return {
            "output_path": output, "library_id": library, "matrix_path": str(matrix_path),
            "n_spots": adata.n_obs, "n_genes": adata.n_vars,
            "matched_barcodes": len(matched), "n_outside_tissue": int((~in_tissue).sum()),
            "n_outside_tissue_retained": int(adata.obs["in_tissue"].eq(0).sum()),
            "coordinate_system": "full_resolution_image_pixels_xy", "image_levels": levels,
        }

    def inspect_spatial_h5ad(self, request):
        import numpy as np

        path = request["input_path"]
        layer = request.get("count_layer", "X")
        mode = request.get("count_check", "sample")
        role = request.get("role", "spatial")
        key = request.get("coordinate_key", "spatial")
        label_key = request.get("cell_type_key")
        with _open_h5ad(path) as adata:
            issues = []
            if not adata.n_obs or not adata.n_vars:
                issues.append("The file must contain observations and genes.")
            duplicate_obs = int(adata.obs_names.duplicated().sum())
            duplicate_genes = int(adata.var_names.duplicated().sum())
            if duplicate_obs or duplicate_genes:
                issues.append("Observation and gene IDs must be unique.")

            matrix = adata.X if layer == "X" else adata.layers.get(layer)
            count = _matrix_check(matrix, mode)
            if matrix is None:
                issues.append(f"Count matrix {layer!r} was not found.")
            elif count["status"] in {"invalid", "non_numeric", "empty"}:
                issues.append(f"Matrix {layer!r} is incompatible with nonnegative integer counts.")

            spatial = {"key": key, "present": key in adata.obsm}
            if spatial["present"] or role == "spatial":
                try:
                    coords = _coordinates(adata, key)
                    invalid = int((~np.isfinite(coords).all(axis=1)).sum())
                    spatial.update(invalid_rows=invalid, valid=not invalid)
                    if invalid:
                        issues.append(f"Coordinate field {key!r} has {invalid} non-finite rows.")
                except ToolValidationError as exc:
                    spatial["valid"] = False
                    issues.append(str(exc))

            labels = {"key": label_key, "checked": label_key is not None}
            if label_key is None and role == "reference":
                issues.append("Set cell_type_key to check reference cell type labels.")
            if label_key is not None:
                if label_key not in adata.obs:
                    labels["present"] = False
                    issues.append(f"Cell type column {label_key!r} was not found.")
                else:
                    series = adata.obs[label_key].astype("string").str.strip()
                    series = series.mask(series.eq(""))
                    missing = int(series.isna().sum())
                    labels.update(
                        present=True, missing=missing,
                        counts={str(k): int(v) for k, v in series.value_counts().items()},
                    )
                    if missing:
                        issues.append(f"Cell type column {label_key!r} has {missing} missing values.")

            return {
                "input_path": str(_path(path)), "role": role, "n_obs": adata.n_obs,
                "n_genes": adata.n_vars, "count_layer": layer, "count_check": {"mode": mode, **count},
                "coordinate_check": spatial, "cell_type_check": labels, "obs_columns": list(adata.obs),
                "layers": list(adata.layers), "obsm_keys": list(adata.obsm),
                "duplicate_obs_ids": duplicate_obs, "duplicate_gene_ids": duplicate_genes,
                "issues": issues,
            }

    def export_spatial_table(self, request):
        import pandas as pd

        adata = _read_h5ad(request["input_path"])
        output = request["output_path"]
        spot_ids = request.get("spot_ids")
        selected = adata.obs_names if spot_ids is None else pd.Index(spot_ids)
        row_indices = adata.obs_names.get_indexer(selected)
        if (row_indices < 0).any():
            raise ToolValidationError(f"Unknown spot IDs: {selected[row_indices < 0].tolist()}.")

        sources = {"spot_id": "obs_names"}
        columns = {}
        if request.get("include_coordinates", True):
            coords = _valid_positions(adata)
            columns.update(x=(coords, 0), y=(coords, 1))
            sources.update(x="obsm/spatial/0", y="obsm/spatial/1")
        for key in request.get("obs_columns", []):
            if key not in adata.obs:
                raise ToolValidationError(f"Observation field obs[{key!r}] was not found.")
            name = f"obs:{key}"
            columns[name] = (adata.obs[key].to_numpy(), None)
            sources[name] = f"obs/{key}"

        for selection in request.get("matrix_columns", []):
            key = selection["key"]
            matrix = adata.obsm.get(key)
            if matrix is None:
                raise ToolValidationError(f"Matrix obsm[{key!r}] was not found.")
            if getattr(matrix, "ndim", 0) != 2:
                raise ToolValidationError(f"Matrix obsm[{key!r}] must be two-dimensional.")
            names = matrix.columns if isinstance(matrix, pd.DataFrame) else None
            if names is not None and not names.is_unique:
                raise ToolValidationError(f"Matrix obsm[{key!r}] has duplicate column labels.")
            for column in selection["columns"]:
                if isinstance(column, int) and 0 <= column < matrix.shape[1]:
                    index = column
                elif names is not None and column in names:
                    index = names.get_loc(column)
                else:
                    raise ToolValidationError(f"Column {column!r} was not found in obsm[{key!r}].")
                label = names[index] if names is not None else index
                name = f"obsm:{key}:{label}"
                if name in columns:
                    raise ToolValidationError(f"Duplicate export column {name!r}.")
                columns[name] = (matrix, index)
                sources[name] = f"obsm/{key}/{label}"

        genes = request.get("genes", [])
        if genes:
            layer = request.get("expression_layer", "X")
            matrix = adata.X if layer == "X" else adata.layers.get(layer)
            if matrix is None:
                raise ToolValidationError(f"Expression matrix {layer!r} was not found.")
            indices = adata.var_names.get_indexer(genes)
            if (indices < 0).any():
                missing = [gene for gene, index in zip(genes, indices) if index < 0]
                raise ToolValidationError(f"Unknown gene IDs: {missing}.")
            for gene, index in zip(genes, indices):
                name = f"gene:{gene}"
                columns[name] = (matrix, index)
                sources[name] = f"X/{gene}" if layer == "X" else f"layers/{layer}/{gene}"

        separator = "\t" if request.get("format", "csv") == "tsv" else ","
        with _new_file(output) as temporary:
            with temporary.open("w", encoding="utf-8", newline="") as stream:
                for start in range(0, len(selected), 2048):
                    rows = row_indices[start:start + 2048]
                    frame = {"spot_id": selected[start:start + len(rows)]}
                    for name, (matrix, index) in columns.items():
                        if index is None:
                            values = matrix[rows]
                        elif isinstance(matrix, pd.DataFrame):
                            values = matrix.iloc[rows, index].to_numpy()
                        else:
                            values = _dense(matrix[rows, index]).reshape(-1)
                        frame[name] = values
                    pd.DataFrame(frame).to_csv(stream, index=False, header=start == 0, sep=separator)
        return {"output_path": str(_path(output)), "n_rows": len(selected), "column_sources": sources}


def _mitochondrial_mask(adata, request):
    import pandas as pd

    flag = request.get("mito_flag_key")
    symbol_key = request.get("gene_symbol_key")
    if flag is not None:
        if flag not in adata.var:
            raise ToolValidationError(f"Mitochondrial flag column var[{flag!r}] was not found.")
        values = adata.var[flag]
        if not pd.api.types.is_bool_dtype(values) or values.isna().any():
            raise ToolValidationError(f"Mitochondrial flag column var[{flag!r}] must contain booleans.")
        return values.to_numpy(dtype=bool)
    if symbol_key is not None and symbol_key not in adata.var:
        raise ToolValidationError(f"Gene symbol column var[{symbol_key!r}] was not found.")
    if symbol_key is None:
        symbol_key = next(
            (key for key in ("gene_symbol", "gene_symbols", "SYMBOL") if key in adata.var), None
        )
    symbols = adata.var_names if symbol_key is None else adata.var[symbol_key]
    mask = pd.Index(symbols, dtype="string").str.upper().str.startswith("MT-", na=False)
    if mask.any():
        return mask.to_numpy(dtype=bool) if hasattr(mask, "to_numpy") else mask
    return None


def _qc_metrics(adata, request):
    import numpy as np
    import pandas as pd
    from scipy.sparse import issparse

    matrix = _count_matrix(adata, request)
    if issparse(matrix):
        matrix = matrix.tocsr(copy=True)
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        detected = np.diff(matrix.indptr)
    else:
        detected = np.count_nonzero(matrix, axis=1)
    totals = np.asarray(matrix.sum(axis=1, dtype=np.float64)).ravel()
    frame = pd.DataFrame({"total_counts": totals, "n_detected_genes": detected}, index=adata.obs_names)
    mito = _mitochondrial_mask(adata, request)
    if mito is not None:
        numerator = np.asarray(matrix[:, mito].sum(axis=1, dtype=np.float64)).ravel()
        frame["pct_mito"] = np.divide(
            100 * numerator, totals, out=np.full(totals.shape, np.nan), where=totals != 0,
        )
    return frame


def _gene_ids(adata, key, strip_version):
    import pandas as pd

    if key == "var_names":
        values = adata.var_names
    elif key in adata.var:
        values = adata.var[key]
    else:
        raise ToolValidationError(f"Gene ID column var[{key!r}] was not found.")
    identifiers = pd.Index(values, dtype="string").str.strip()
    if identifiers.isna().any() or (identifiers == "").any():
        raise ToolValidationError(f"Gene ID field {key!r} contains missing or empty values.")
    if strip_version:
        identifiers = identifiers.str.replace(r"^(ENS[A-Z]*G[0-9]+)\.[0-9]+$", r"\1", regex=True)
    return pd.Index(identifiers.astype(str))


def _read_gene_mapping(path, strip_version):
    import pandas as pd

    path = _path(path)
    try:
        table = pd.read_csv(
            path, sep="\t" if path.suffix.lower() == ".tsv" else ",", dtype=str, keep_default_na=False,
        )
        table = table[["source_id", "target_id"]].apply(lambda column: column.str.strip())
    except (OSError, ValueError, KeyError) as exc:
        raise ToolValidationError("The mapping file must contain source_id and target_id columns.") from exc
    if table.empty or table.eq("").any().any():
        raise ToolValidationError("The mapping file must contain nonempty source_id and target_id values.")
    if strip_version:
        table = table.replace(r"^(ENS[A-Z]*G[0-9]+)\.[0-9]+$", r"\1", regex=True)
    table = table.drop_duplicates()
    duplicates = table.loc[table["source_id"].duplicated(), "source_id"]
    if len(duplicates):
        raise ToolValidationError(f"Ambiguous source IDs in mapping: {duplicates.head(5).tolist()}.")
    return dict(zip(table["source_id"], table["target_id"]))


def _align_side(adata, key, strip_version, mapping, side):
    import pandas as pd

    source = _gene_ids(adata, key, strip_version)
    final = source
    matched = None
    if mapping is not None:
        matched = source.isin(mapping)
        final = pd.Index([mapping.get(identifier, identifier) for identifier in source])
    duplicates = final[final.duplicated()].unique()
    if len(duplicates):
        raise ToolValidationError(f"Duplicate target gene IDs on the {side} side: {duplicates[:5].tolist()}.")
    audit = pd.DataFrame({
        "side": side, "original_var_name": adata.var_names.to_numpy(), "source_id": source,
        "aligned_id": final, "mapping_status": (
            ["not_requested"] * len(source) if matched is None
            else ["mapped" if item else "unmapped" for item in matched]
        ),
    })
    adata.var_names = final
    return audit


@register_tool("SpatialPreprocessTool")
class SpatialPreprocessTool(_SpatialTool):
    operations = ("filter_spatial_qc", "align_gene_ids", "normalize_spatial_expression")

    def filter_spatial_qc(self, request):
        adata = _read_h5ad(request["input_path"])
        frame = _qc_metrics(adata, request)
        group = request.get("group_by")
        if group is not None:
            if group not in adata.obs:
                raise ToolValidationError(f"Grouping column obs[{group!r}] was not found.")
            frame["group"] = adata.obs[group].astype("string")
        rules = {
            "min_counts": ("total_counts", "ge"),
            "max_counts": ("total_counts", "le"),
            "min_detected_genes": ("n_detected_genes", "ge"),
            "max_detected_genes": ("n_detected_genes", "le"),
            "max_pct_mito": ("pct_mito", "le"),
        }
        thresholds = {key: request[key] for key in rules if key in request}
        for lower, upper in (("min_counts", "max_counts"), ("min_detected_genes", "max_detected_genes")):
            if lower in thresholds and upper in thresholds and thresholds[lower] > thresholds[upper]:
                raise ToolValidationError(f"{lower} must be less than or equal to {upper}.")
        if "max_pct_mito" in thresholds and "pct_mito" not in frame:
            raise ToolValidationError("Mitochondrial genes could not be identified. Set mito_flag_key or gene_symbol_key.")

        failures = {}
        for name, threshold in thresholds.items():
            column, method = rules[name]
            failures[name] = ~getattr(frame[column], method)(threshold)
        frame["failed_rules"] = [
            ";".join(name for name, flag in zip(failures, row) if flag)
            for row in zip(*(series.to_numpy() for series in failures.values()))
        ] if failures else ""
        frame["passes_qc"] = frame["failed_rules"].eq("")
        n_pass = int(frame["passes_qc"].sum())
        report = {
            "mode": request.get("mode", "report"), "n_spots": len(frame),
            "n_pass": n_pass, "n_fail": len(frame) - n_pass,
            "mitochondrial_metric_available": "pct_mito" in frame,
            "failed_by_rule": {key: int(value.sum()) for key, value in failures.items()},
            "groups": [],
        }
        if group is not None:
            for name, subset in frame.groupby("group", dropna=False, sort=False):
                n = int(subset["passes_qc"].sum())
                report["groups"].append({
                    "group": None if subset["group"].isna().all() else str(name),
                    "n_spots": len(subset), "n_pass": n, "pass_fraction": n / len(subset),
                })

        output = request.get("output_path")
        if report["mode"] == "filter":
            if not thresholds:
                raise ToolValidationError("Filter mode requires at least one QC threshold.")
            if n_pass == 0:
                raise ToolValidationError("No spots pass the requested QC thresholds.")
            filtered = adata[frame["passes_qc"].to_numpy()].copy()
            for key in ("total_counts", "n_detected_genes", "pct_mito"):
                if key in frame:
                    filtered.obs[key] = frame.loc[filtered.obs_names, key]
            with _result_directory(_path(output).parent, "spatial_qc") as directory:
                audit = directory / "spot_qc.csv"
                frame.to_csv(audit, index_label="spot_id")
                report["output_path"] = _write_h5ad(filtered, output)
        elif output:
            audit = _path(output)
            with _new_file(audit) as temporary:
                frame.to_csv(temporary, index_label="spot_id")
        else:
            with _result_directory(tempfile.gettempdir(), "spatial_qc") as directory:
                audit = directory / "spot_qc.csv"
                frame.to_csv(audit, index_label="spot_id")
        return {**report, "report_path": str(audit)}

    def align_gene_ids(self, request):
        import pandas as pd

        # Existing duplicate var_names can be repaired with a unique var column.
        inputs = {}
        for side in ("spatial", "reference"):
            with _open_h5ad(request[f"{side}_path"]) as source:
                if not source.n_obs or not source.n_vars or not source.obs_names.is_unique:
                    raise ToolValidationError(f"The {side} input requires unique observations and nonempty axes.")
                inputs[side] = source.to_memory()
        strip = request.get("strip_ensembl_version", False)
        mapping_side = request.get("mapping_side")
        mapping_path = request.get("mapping_path")
        mapping = _read_gene_mapping(mapping_path, strip) if mapping_path else None
        audit = []
        for side, adata in inputs.items():
            audit.append(_align_side(
                adata, request.get(f"{side}_id_key", "var_names"), strip,
                mapping if side == mapping_side else None, side,
            ))

        spatial, reference = inputs["spatial"], inputs["reference"]
        shared = spatial.var_names.intersection(reference.var_names, sort=False)
        if not len(shared):
            raise ToolValidationError("No shared gene IDs remain after alignment.")
        for frame in audit:
            frame["shared"] = frame["aligned_id"].isin(shared)
        mode = request.get("mode", "preserve")
        if mode == "intersection":
            spatial = spatial[:, shared].copy()
            reference = reference[:, shared].copy()

        with _result_directory(request["output_dir"], "align_gene_ids") as directory:
            spatial_path = directory / "spatial_aligned.h5ad"
            reference_path = directory / "reference_aligned.h5ad"
            report_path = directory / "gene_alignment.csv"
            spatial.write_h5ad(spatial_path, compression="gzip")
            reference.write_h5ad(reference_path, compression="gzip")
            details = pd.concat(audit, ignore_index=True)
            details.to_csv(report_path, index=False)
        return {
            "spatial_path": str(spatial_path), "reference_path": str(reference_path),
            "report_path": str(report_path), "mode": mode, "n_shared": len(shared),
            "n_spatial_only": int((~audit[0]["shared"]).sum()),
            "n_reference_only": int((~audit[1]["shared"]).sum()),
            "n_unmapped": int(details["mapping_status"].eq("unmapped").sum()),
            "n_spatial_genes_written": spatial.n_vars, "n_reference_genes_written": reference.n_vars,
        }

    def normalize_spatial_expression(self, request):
        import anndata as ad
        import numpy as np
        import scanpy as sc

        adata = _read_h5ad(request["input_path"])
        matrix = _count_matrix(adata, request)
        target = request.get("target_sum", 10000)
        if not np.isfinite(target) or target <= 0:
            raise ToolValidationError("target_sum must be finite and positive.")
        log1p = request.get("log1p", True)
        layer = request.get("output_layer", "lognorm" if log1p else "normalized")
        factor = request.get("size_factor_key")
        if layer == "X" or layer in adata.layers:
            raise ToolValidationError(f"Output layer {layer!r} already exists or is reserved.")
        if factor is not None and factor in adata.obs:
            raise ToolValidationError(f"Size factor column obs[{factor!r}] already exists.")

        totals = np.asarray(matrix.sum(axis=1, dtype=np.float64)).ravel()
        working = ad.AnnData(X=matrix.astype(np.float64, copy=True))
        if np.any(totals):
            sc.pp.normalize_total(working, target_sum=target)
        if log1p:
            sc.pp.log1p(working)
        adata.layers[layer] = working.X
        if factor is not None:
            adata.obs[factor] = totals / target
        return {
            "output_path": _write_h5ad(adata, request["output_path"]), "output_layer": layer,
            "count_layer": request.get("count_layer", "X"), "target_sum": target,
            "log1p": log1p, "n_zero_spots": int((totals == 0).sum()),
        }
