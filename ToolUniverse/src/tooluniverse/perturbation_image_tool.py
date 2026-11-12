"""Index local microscopy images without retaining pixel arrays in tool responses."""

import hashlib
import json
import re
import uuid
from pathlib import Path

from .base_tool import BaseTool
from .perturbation_dataset_tool import DatasetLoadError, _fingerprint
from .tool_registry import register_tool

_SUFFIXES = {
    ".tif": "tiff",
    ".tiff": "tiff",
    ".tf2": "tiff",
    ".tf8": "tiff",
    ".btf": "tiff",
    ".png": "png",
    ".jpg": "jpeg",
    ".jpeg": "jpeg",
    ".bmp": "bmp",
    ".czi": "czi",
    ".nd2": "nd2",
    ".lif": "lif",
}


def _is_store(path):
    return path.is_dir() and any(
        (path / name).is_file() for name in (".zgroup", "zarr.json")
    )


def _image_fingerprint(path):
    if not path.is_dir():
        return _fingerprint(path)
    digest = hashlib.sha256()
    members = sorted(member for member in path.rglob("*") if member.is_file())
    if not members:
        raise DatasetLoadError("INVALID_DATA", "Empty local image store.")
    total = 0
    for member in members:
        stamp = _fingerprint(member)
        total += stamp["bytes"]
        stamp["path"] = str(member.relative_to(path))
        digest.update(json.dumps(stamp, sort_keys=True).encode())
    return {
        "path": str(path),
        "kind": "directory",
        "n_files": len(members),
        "bytes": total,
        "sha256": digest.hexdigest(),
    }


def _directory_files(source, pattern=None):
    candidates = sorted(source.glob(pattern or "*"))
    files = []
    for path in candidates:
        if any(store in path.parents for store in files if store.is_dir()):
            continue
        if _is_store(path) or (path.is_file() and path.suffix.lower() in _SUFFIXES):
            files.append(path.resolve())
    if not files:
        raise DatasetLoadError("INVALID_DATA", "No matching images found.")
    return files


def _image_header(path, arguments):
    from .microscopy_image_readers import read_native, read_tiff, read_zarr

    fmt = (
        "ome_zarr"
        if isinstance(path, str) or _is_store(path)
        else _SUFFIXES.get(path.suffix.lower())
    )
    if fmt is None:
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT", f"Unsupported image format: {path}"
        )
    if fmt in ("tiff", "ome_zarr", "czi", "nd2", "lif"):
        if fmt == "ome_zarr":
            header = read_zarr(path, arguments)
        elif fmt == "tiff":
            header = read_tiff(path, arguments)
        else:
            header = read_native(path, arguments)
        shape, axes, dtype = header["shape"], header["axes"], header["dtype"]
        fmt = header["format"]
    else:
        if any(
            key in arguments
            for key in (
                "series_index",
                "resolution_level",
                "group_path",
                "position_index",
            )
        ):
            raise DatasetLoadError(
                "INVALID_ARGUMENT",
                "Series, resolution and group selectors do not apply to raster files.",
            )
        from PIL import Image

        with Image.open(path) as image:
            if image.format != {"png": "PNG", "jpeg": "JPEG", "bmp": "BMP"}[fmt]:
                raise DatasetLoadError(
                    "INVALID_DATA", "File content does not match its extension."
                )
            dtypes = {
                "1": "bool",
                "L": "uint8",
                "LA": "uint8",
                "RGB": "uint8",
                "RGBA": "uint8",
                "I": "int32",
                "I;16": "uint16",
                "I;16B": "uint16",
            }
            if image.mode not in dtypes:
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT", f"Unsupported raster mode: {image.mode}"
                )
            bands = len(image.getbands())
            shape = [image.height, image.width] + ([bands] if bands > 1 else [])
            axes, dtype = ("YXC" if bands > 1 else "YX"), dtypes[image.mode]
            image.load() if not arguments.get("lazy", True) else image.verify()
        header = {"native_channels": None, "image_metadata": {}}
    requested = arguments.get("format", "auto")
    if requested not in ("auto", fmt) and not (
        requested == "tiff" and fmt == "ome_tiff"
    ):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT", f"Detected {fmt}, requested {requested}."
        )
    if fmt != "ome_zarr" and "group_path" in arguments:
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "group_path applies only to OME-Zarr."
        )
    if (
        fmt in ("ome_tiff", "ome_zarr", "czi", "nd2", "lif")
        and arguments.get("axes", axes) != axes
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "axes cannot override authoritative OME axis metadata."
        )
    if fmt in ("czi", "nd2", "lif") and not set(axes) <= set("TCZYX"):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT", f"Unsupported vendor acquisition axes: {axes}"
        )
    axes = arguments.get("axes", axes)
    if (
        len(axes) != len(shape)
        or len(set(axes)) != len(axes)
        or not set(axes) <= set("TCZYX")
        or not {"Y", "X"} <= set(axes)
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            f"Ambiguous axes {axes!r} for shape {shape}; provide explicit axes using TCZYX.",
        )
    if min(shape) <= 0:
        raise DatasetLoadError("INVALID_DATA", "Image dimensions must be positive.")
    channel_count = shape[axes.index("C")] if "C" in axes else 1
    mapping = arguments.get("channel_map", {})
    if mapping and set(mapping) != {str(i) for i in range(channel_count)}:
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "channel_map must name every zero-based channel exactly once.",
        )
    native = header.get("native_channels")
    channels = [
        mapping.get(
            str(i), native[i] if native and len(native) == channel_count else str(i)
        )
        for i in range(channel_count)
    ]
    if len(set(channels)) != len(channels):
        raise DatasetLoadError("INVALID_ARGUMENT", "Channel names must be unique.")
    return {
        "format": fmt,
        "shape": shape,
        "axes": axes,
        "dtype": dtype,
        "channels": channels,
        "spatial_shape": [shape[axes.index("Y")], shape[axes.index("X")]],
        "image_metadata": header["image_metadata"],
        "related_files": header.get("related_files", []),
    }


def _manifest_rows(source, arguments):
    import csv

    import pandas as pd

    if source.suffix.lower() not in (".csv", ".tsv"):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT", "Image manifest must be CSV or TSV."
        )
    separator = "\t" if source.suffix.lower() == ".tsv" else ","
    with source.open(newline="", encoding="utf-8-sig") as stream:
        header = next(csv.reader(stream, delimiter=separator), [])
    if (
        not header
        or any(not value for value in header)
        or len(header) != len(set(header))
    ):
        raise DatasetLoadError(
            "METADATA_MISMATCH", "Manifest column names must be unique and non-empty."
        )
    table = pd.read_csv(
        source,
        sep="\t" if source.suffix.lower() == ".tsv" else ",",
        dtype=str,
        keep_default_na=False,
    )
    paths = arguments.get("image_path_columns", {})
    mapping = arguments.get("metadata_column_map", {})
    if not paths:
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "Manifest input requires image_path_columns."
        )
    if arguments.get("channel_map"):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "Use manifest channel names or channel_map, not both."
        )
    missing = (set(paths.values()) | set(mapping.values())) - set(table.columns)
    if missing:
        raise DatasetLoadError(
            "METADATA_MISMATCH", f"Missing manifest columns: {sorted(missing)}"
        )
    keys = arguments.get("field_keys", [])
    if set(keys) - set(mapping):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "field_keys must refer to mapped metadata fields."
        )
    root = (
        Path(arguments["image_root"]).expanduser().resolve()
        if arguments.get("image_root")
        else source.parent
    )
    # Flat manifests with many missing downloads otherwise issue thousands of
    # negative GPFS lookups. Snapshot names only for the explicit skip mode.
    root_names = (
        {entry.name for entry in root.iterdir()}
        if arguments.get("missing_images") == "skip_field" and root.is_dir()
        else None
    )
    rows, skipped, identities = [], 0, set()
    for number, row in table.iterrows():
        metadata = {standard: row[column] for standard, column in mapping.items()}
        if any(not metadata[key] for key in keys):
            raise DatasetLoadError(
                "METADATA_MISMATCH", "Field key values must be non-empty."
            )
        identity = tuple(metadata[key] for key in keys) if keys else (str(number),)
        if identity in identities:
            raise DatasetLoadError(
                "METADATA_MISMATCH", "Duplicate field identity in manifest."
            )
        identities.add(identity)
        members = []
        for channel, column in paths.items():
            value = row[column]
            path = Path(value).expanduser() if value else None
            if path is not None and not path.is_absolute():
                path = root / path
            members.append((channel, path))

        def available(path):
            if path is None:
                return False
            if (
                root_names is not None
                and path.parent == root
                and path.name not in root_names
            ):
                return False
            return path.is_file() or _is_store(path)

        if any(not available(path) for _, path in members):
            if arguments.get("missing_images", "error") == "skip_field":
                skipped += 1
                continue
            raise DatasetLoadError(
                "SOURCE_NOT_FOUND", f"Missing image in manifest row {number}."
            )
        rows.append(
            (
                str(number),
                [(channel, path.resolve()) for channel, path in members],
                metadata,
            )
        )
    if not rows:
        raise DatasetLoadError("INVALID_DATA", "No complete image fields remain.")
    return rows, skipped


def _inspect_images(arguments):
    if arguments["source_type"] == "remote":
        from .microscopy_remote import remote_metadata_stamp

        source = arguments["source"]
        if any(
            key in arguments
            for key in (
                "image_path_columns",
                "metadata_column_map",
                "image_root",
                "field_keys",
                "missing_images",
                "file_pattern",
                "position_index",
            )
        ):
            raise DatasetLoadError(
                "INVALID_ARGUMENT",
                "Remote mode accepts one anonymous OME-Zarr image, not local manifest options.",
            )
        before = remote_metadata_stamp(source, arguments)
        header = _image_header(source, arguments)
        if remote_metadata_stamp(source, arguments) != before:
            raise DatasetLoadError(
                "INVALID_DATA", "Remote metadata changed during loading."
            )
        summary = {
            "modality": "image",
            "n_files": 1,
            "n_fields": None,
            "n_input_files": 1,
            "channels": header["channels"],
            "metadata_fields": [],
            "integrity_scope": "metadata_only",
            "pixel_validation": "header" if arguments.get("lazy", True) else "decoded",
        }
        record = {"path": source, "field_id": None, "metadata": {}, **header}
        return (
            summary,
            [record],
            [before],
            [
                {
                    "code": "REMOTE_METADATA_ONLY",
                    "message": "Remote pixels are not content-fingerprinted; source chunk changes without metadata changes are not detected.",
                }
            ],
        )
    if "remote_timeout" in arguments or "remote_endpoint" in arguments:
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "remote_timeout and remote_endpoint apply only to remote mode.",
        )
    source = Path(arguments["source"]).expanduser().resolve()
    if not source.exists():
        raise DatasetLoadError("SOURCE_NOT_FOUND", f"Source does not exist: {source}")
    mode = arguments["source_type"]
    if mode == "directory":
        if not source.is_dir():
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "directory source_type requires a directory."
            )
        rows = [
            (None, [(None, p)], {})
            for p in _directory_files(source, arguments.get("file_pattern"))
        ]
        inputs = []
    elif mode == "file":
        if not (source.is_file() or _is_store(source)):
            raise DatasetLoadError(
                "INVALID_ARGUMENT",
                "file source_type requires an image file or local OME-Zarr store.",
            )
        rows, inputs = [(None, [(None, source)], {})], []
    elif mode == "manifest":
        if not source.is_file():
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "manifest source_type requires a file."
            )
        inputs = [_fingerprint(source)]
        rows, skipped = _manifest_rows(source, arguments)
    else:
        raise DatasetLoadError("INVALID_ARGUMENT", "Invalid source_type.")
    if mode != "manifest" and any(
        k in arguments
        for k in (
            "image_path_columns",
            "metadata_column_map",
            "image_root",
            "field_keys",
            "missing_images",
        )
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "Manifest options require source_type=manifest."
        )
    if mode != "directory" and "file_pattern" in arguments:
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "file_pattern requires source_type=directory."
        )
    records, seen, channels = [], set(), set()
    for field_id, members, metadata in rows:
        shapes = set()
        for channel, path in members:
            if path in seen and mode != "directory":
                raise DatasetLoadError(
                    "INVALID_DATA", f"Image appears more than once: {path}"
                )
            from .microscopy_image_readers import ome_tiff_members

            related = (
                ome_tiff_members(path)
                if _SUFFIXES.get(path.suffix.lower()) == "tiff"
                else [path]
            )
            related_set = set(related)
            if mode == "directory" and related_set.issubset(seen):
                continue
            if related_set & seen:
                raise DatasetLoadError(
                    "INVALID_DATA",
                    "Overlapping image filesets were declared more than once.",
                )
            seen.update(related)
            stamps = [_image_fingerprint(member) for member in related]
            header = _image_header(path, arguments)
            if channel is not None:
                if len(header["channels"]) != 1:
                    raise DatasetLoadError(
                        "INVALID_ARGUMENT",
                        "Wide manifest channels require single-channel image files.",
                    )
                header["channels"] = [channel]
            shapes.add(
                tuple(
                    (axis, size)
                    for axis, size in zip(header["axes"], header["shape"], strict=True)
                    if axis != "C"
                )
            )
            channels.update(header["channels"])
            inputs.extend(stamps)
            records.append(
                {
                    "path": str(path),
                    "field_id": field_id,
                    "metadata": metadata,
                    **header,
                }
            )
        if not shapes and mode == "directory":
            continue
        if len(shapes) != 1:
            raise DatasetLoadError(
                "INVALID_DATA", f"Channel image dimensions differ in field {field_id}."
            )
    for stamp in inputs:
        if _image_fingerprint(Path(stamp["path"])) != stamp:
            raise DatasetLoadError("INVALID_DATA", "An input changed while loading.")
    warnings = []
    if mode == "manifest" and skipped:
        warnings.append({"code": "SKIPPED_FIELDS", "count": skipped})
    if mode != "manifest":
        warnings.append(
            {
                "code": "FIELD_GROUPING_UNKNOWN",
                "message": "File names are not interpreted as experimental field identities.",
            }
        )
    summary = {
        "modality": "image",
        "n_files": len(records),
        "n_input_files": len(seen),
        "integrity_scope": "full_content",
        "n_fields": len(rows) if mode == "manifest" else None,
        "channels": sorted(channels),
        "metadata_fields": list(arguments.get("metadata_column_map", {})),
        "pixel_validation": "header" if arguments.get("lazy", True) else "decoded",
    }
    return summary, records, inputs, warnings


def resolve_image_dataset_ref(dataset_ref):
    """Validate originals and the persistent JSONL index, returning the descriptor."""
    descriptor = json.loads(Path(dataset_ref).read_text())
    if descriptor.get("schema_version") != 1 or descriptor.get("modality") != "image":
        raise DatasetLoadError("INVALID_DATA", "Unsupported image reference.")
    for stamp in descriptor["inputs"] + [descriptor["index"]]:
        if stamp.get("kind") == "remote":
            from .microscopy_remote import remote_metadata_stamp

            if remote_metadata_stamp(stamp["path"], stamp) != stamp:
                raise DatasetLoadError(
                    "INVALID_DATA", "Remote metadata changed; load again."
                )
            continue
        path = Path(stamp["path"])
        if not path.exists():
            raise DatasetLoadError(
                "SOURCE_NOT_FOUND", f"Referenced file missing: {path}"
            )
        if _image_fingerprint(path) != stamp:
            raise DatasetLoadError(
                "INVALID_DATA", "Referenced input or index changed; load again."
            )
    options = descriptor["read_options"]
    if options["source_type"] == "directory":
        files = {
            str(p)
            for p in _directory_files(
                Path(options["source"]), options.get("file_pattern")
            )
        }
        if files != {s["path"] for s in descriptor["inputs"]}:
            raise DatasetLoadError(
                "INVALID_DATA", "Directory image membership changed; load again."
            )
    return descriptor


@register_tool("PerturbationImageDatasetTool")
class PerturbationImageDatasetTool(BaseTool):
    def run(self, arguments):
        directory = None
        try:
            dataset_id = arguments.get("dataset_id", uuid.uuid4().hex)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,127}", dataset_id):
                raise DatasetLoadError("INVALID_ARGUMENT", "Invalid dataset_id.")
            summary, records, inputs, warnings = _inspect_images(arguments)
            destination = (
                Path(arguments["output_dir"]).expanduser().resolve() / dataset_id
            )
            if any(
                stamp.get("kind") == "directory"
                and destination.is_relative_to(Path(stamp["path"]))
                for stamp in inputs
            ):
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "Output must be outside input OME-Zarr stores."
                )
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                destination.mkdir()
            except FileExistsError:
                raise DatasetLoadError(
                    "DATASET_EXISTS", f"Dataset already exists: {destination}"
                ) from None
            directory = destination
            index = directory / "images.jsonl"
            with index.open("w") as stream:
                for record in records:
                    stream.write(json.dumps(record) + "\n")
            options = dict(arguments)
            for key in ("source", "image_root", "output_dir"):
                if key in options and not (
                    key == "source" and options.get("source_type") == "remote"
                ):
                    options[key] = str(Path(options[key]).expanduser().resolve())
            descriptor = {
                "schema_version": 1,
                "dataset_id": dataset_id,
                **summary,
                "inputs": inputs,
                "index": _fingerprint(index),
                "read_options": options,
                "warnings": warnings,
            }
            temporary = directory / "manifest.json.part"
            temporary.write_text(json.dumps(descriptor, indent=2) + "\n")
            manifest = directory / "manifest.json"
            temporary.replace(manifest)
            return {
                "status": "success",
                "data": {
                    "dataset_ref": str(manifest),
                    "index_path": str(index),
                    "dataset_id": dataset_id,
                    **summary,
                },
                "warnings": warnings,
            }
        except DatasetLoadError as exc:
            return {"status": "error", "error_code": exc.code, "error": str(exc)}
        except ImportError as exc:
            return {
                "status": "error",
                "error_code": "DEPENDENCY_MISSING",
                "error": str(exc),
            }
        except OSError as exc:
            return {"status": "error", "error_code": "IO_ERROR", "error": str(exc)}
        except (ValueError, KeyError, TypeError) as exc:
            return {"status": "error", "error_code": "INVALID_DATA", "error": str(exc)}
        finally:
            if directory is not None and not (directory / "manifest.json").exists():
                for name in ("images.jsonl", "manifest.json.part"):
                    (directory / name).unlink(missing_ok=True)
                directory.rmdir()
