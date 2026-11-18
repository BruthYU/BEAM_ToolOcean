"""Local Cellpose segmentation with legacy and persistent batch interfaces.

Complete mode exports real selected pixels, uint32 instance masks, object tables,
previews and optional native flows/styles. Batches process logical images in
order, reuse a device/model cache, and publish resumable per-image results.
The image-reference bridge requires explicit non-singleton Z/T selections.
Cellpose 4.x rejects explicit legacy model_type/channels; physical channel
selection is handled by channel_names/channel_indices before inference. No volumetric
inference, training, or background job service is exposed.
"""

import math
import os
from threading import RLock
from typing import ClassVar

from .base_tool import BaseTool
from .tool_registry import register_tool

# Optional dependency import at module load so a missing package becomes a clean
# error rather than an exception (framework optional-dependency pattern).
CELLPOSE_AVAILABLE = False
_IMPORT_ERROR = None
try:
    import numpy as np
    from cellpose import models as _cp_models

    CELLPOSE_AVAILABLE = True
except Exception as exc:  # noqa: BLE001 - optional third-party imports can fail downstream
    np = None
    _cp_models = None
    _IMPORT_ERROR = str(exc)

# Legacy cellpose-3.x model types accepted from the caller. In 4.x these are
# ignored by the unified model but still accepted so callers/tests are stable.
_LEGACY_MODEL_TYPES = ("cyto3", "cyto2", "cyto", "nuclei")

# Image extensions we know how to load.
_TIFF_EXTS = (".tif", ".tiff")
_PIL_EXTS = (".png", ".jpg", ".jpeg", ".bmp")


@register_tool("CellposeTool")
class CellposeTool(BaseTool):
    """Local 2D segmentation with legacy single-image and complete result modes.

    Provide exactly one image_path, image_paths, or image_ref. output_dir enables
    persistent masks, object tables and optional native outputs, with per-image
    failure handling and validated resume. See data/cellpose_tools.json and
    examples/cellpose/TOOL_SPEC.md for parameters and result contracts.
    """

    # Cache loaded models per (class, key) — model construction (and the 4.x
    # weight load) is expensive, so reuse within a process.
    _model_cache: ClassVar[dict] = {}
    _model_weights: ClassVar[dict] = {}
    _model_lock: ClassVar[RLock] = RLock()

    @staticmethod
    def _load_image(image_path):
        """Load an image into a 2D/3D numpy array. Returns (array, error_or_None)."""
        ext = os.path.splitext(image_path)[1].lower()
        try:
            if ext in _TIFF_EXTS:
                try:
                    import tifffile

                    with tifffile.TiffFile(image_path) as handle:
                        series = handle.series[0]
                        if any(
                            axis in "TZIQ" and size > 1
                            for axis, size in zip(series.axes, series.shape)
                        ):
                            return None, (
                                "Only a single 2D image is supported; select one "
                                "Z/T plane before calling Cellpose_segment_image."
                            )
                        arr = series.asarray()
                except ImportError:
                    from PIL import Image

                    arr = np.array(Image.open(image_path))
            elif ext in _PIL_EXTS:
                from PIL import Image

                arr = np.array(Image.open(image_path))
            else:
                # Last-ditch attempt via PIL; report unsupported extension if it fails.
                try:
                    from PIL import Image

                    arr = np.array(Image.open(image_path))
                except (OSError, ValueError):
                    return None, (
                        f"Unsupported image extension '{ext}'. Supported: "
                        + ", ".join(_TIFF_EXTS + _PIL_EXTS)
                    )
        except Exception as exc:  # noqa: BLE001 - return third-party failures as tool errors
            return None, f"Failed to read image '{image_path}': {exc}"

        if arr is None or getattr(arr, "size", 0) == 0:
            return None, f"Image '{image_path}' is empty or could not be decoded."
        if arr.ndim not in (2, 3) or (
            arr.ndim == 3 and arr.shape[0] > 4 and arr.shape[-1] > 4
        ):
            return None, "Only a single 2D grayscale or channel image is supported."
        if not np.issubdtype(arr.dtype, np.number) or not np.isfinite(arr).all():
            return None, "Image pixels must be finite numeric values."
        return arr, None

    @staticmethod
    def _object_stats(masks):
        """Compute count and per-object area/centroid from an integer label mask."""
        labels, inverse = np.unique(masks, return_inverse=True)
        inverse = inverse.ravel()
        positions = np.arange(masks.size)
        counts = np.bincount(inverse)
        sums_y = np.bincount(inverse, weights=positions // masks.shape[1])
        sums_x = np.bincount(inverse, weights=positions % masks.shape[1])
        return [
            {
                "label": int(label),
                "area": int(counts[index]),
                "centroid_y": float(sums_y[index] / counts[index]),
                "centroid_x": float(sums_x[index] / counts[index]),
            }
            for index, label in enumerate(labels)
            if label != 0
        ]

    def _save_mask(self, masks, image_path, mask_output_path):
        """Write the label mask to disk. Returns (path_or_None, error_or_None)."""
        out_path = mask_output_path or (
            os.path.splitext(image_path)[0] + "_cp_masks.png"
        )
        try:
            from PIL import Image

            if os.path.realpath(out_path) == os.path.realpath(image_path):
                return None, "Mask output must not overwrite the source image."
            if int(masks.max(initial=0)) > np.iinfo(np.uint16).max:
                return (
                    None,
                    "Mask labels exceed the 16-bit limit (65535); mask was not saved.",
                )
            # 16-bit grayscale preserves up to 65535 distinct object labels.
            mask16 = masks.astype(np.uint16)
            Image.fromarray(mask16).save(out_path)
            return out_path, None
        except Exception as exc:  # noqa: BLE001 - return third-party failures as tool errors
            return None, f"Failed to save mask to '{out_path}': {exc}"

    @staticmethod
    def _resolve_device(value):
        import torch

        if value == "auto":
            return "cuda:0" if torch.cuda.is_available() else "cpu"
        if value == "cpu":
            return value
        if value == "mps":
            if not torch.backends.mps.is_available():
                raise ValueError("Requested MPS device is unavailable.")
            return value
        if value == "cuda" or value.startswith("cuda:"):
            index = int(value.split(":")[1]) if ":" in value else 0
            if not torch.cuda.is_available() or index >= torch.cuda.device_count():
                raise ValueError("Requested CUDA device is unavailable.")
            return f"cuda:{index}"
        raise ValueError("device must be auto, cpu, cuda, cuda:N or mps.")

    @classmethod
    def _get_model(cls, model_type, device="cpu", pretrained_model=None):
        import inspect
        from pathlib import Path

        import torch

        from .perturbation_dataset_tool import _fingerprint

        legacy = hasattr(_cp_models, "Cellpose")
        if legacy and pretrained_model:
            raise ValueError(
                "pretrained_model selection requires Cellpose 4.x in this wrapper."
            )
        if pretrained_model:
            path = Path(pretrained_model).expanduser()
            if path.is_file():
                pretrained_model = str(path.resolve())
            elif pretrained_model not in getattr(_cp_models, "MODEL_NAMES", []):
                raise ValueError(
                    "Unknown pretrained_model; use a supported built-in name or existing local file."
                )
        key = (
            "Cellpose" if legacy else "CellposeModel",
            model_type if legacy else pretrained_model,
            device,
        )
        with cls._model_lock:
            cached = cls._model_cache.get(key)
            if cached:
                if (
                    cached[2] is None
                    or _fingerprint(Path(cached[2]["path"])) == cached[2]
                ):
                    return cached[:2]
                del cls._model_cache[key]
            constructor = _cp_models.Cellpose if legacy else _cp_models.CellposeModel
            options = {"gpu": device != "cpu"}
            if "device" in inspect.signature(constructor).parameters:
                options["device"] = torch.device(device)
            elif device not in ("cpu", "cuda:0"):
                raise ValueError(
                    "Installed Cellpose cannot select this device explicitly."
                )
            if legacy:
                options["model_type"] = model_type
            elif pretrained_model:
                options["pretrained_model"] = pretrained_model
            model = constructor(**options)
            actual = getattr(model, "pretrained_model", None)
            if legacy:
                actual = getattr(model.cp, "pretrained_model", actual)
            if isinstance(actual, (list, tuple)):
                actual = actual[0] if actual else None
            weight = (
                _fingerprint(Path(actual))
                if actual and Path(actual).is_file()
                else None
            )
            used = (
                model_type
                if legacy
                else os.path.basename(str(actual or "cellpose-default"))
            )
            cls._model_cache[key] = (model, used, weight)
            cls._model_weights[id(model)] = weight
            return model, used

    @classmethod
    def _infer(cls, image, arguments, backend):
        import inspect

        model, _used = backend
        options = {
            "diameter": arguments.get("diameter") or None,
            "batch_size": arguments.get("batch_size", 8),
            "flow_threshold": arguments.get("flow_threshold", 0.4),
            "cellprob_threshold": arguments.get("cellprob_threshold", 0.0),
        }
        if hasattr(_cp_models, "Cellpose"):
            options["channels"] = arguments.get("channels", [0, 0])
        if "channel_axis" in inspect.signature(model.eval).parameters:
            options["channel_axis"] = -1 if image.ndim == 3 else None
        with cls._model_lock:
            native = model.eval(image.copy(), **options)
        masks = np.asarray(native[0])
        if (
            masks.ndim != 2
            or not np.issubdtype(masks.dtype, np.integer)
            or (masks < 0).any()
        ):
            raise ValueError("Cellpose returned an invalid 2D integer label mask.")
        if list(masks.shape) != list(image.shape[:2]):
            raise ValueError("Cellpose mask does not match the selected image plane.")
        return native

    @classmethod
    def _model_info(cls, backend, arguments, device):
        model, used = backend
        return cls._model_details(
            used, cls._model_weights.get(id(model)), arguments, device
        )

    @staticmethod
    def _model_details(used, weights, arguments, device):
        import importlib.metadata

        return {
            "model_used": used,
            "cellpose_version": importlib.metadata.version("cellpose"),
            "device": device,
            "weights": weights,
            "model_type_requested": arguments.get("model_type", "cyto3"),
            "pretrained_model_requested": arguments.get("pretrained_model"),
            "ignored_parameters": (
                ["model_type"] + (["channels"] if "channels" in arguments else [])
            )
            if not hasattr(_cp_models, "Cellpose")
            else [],
            "diameter": arguments.get("diameter") or None,
            "flow_threshold": arguments.get("flow_threshold", 0.4),
            "cellprob_threshold": arguments.get("cellprob_threshold", 0.0),
            "batch_size": arguments.get("batch_size", 8),
        }

    @classmethod
    def _resume_model_info(cls, arguments, device):
        """Identify installed 4.x weights without constructing/loading a network.

        Be conservative on unknown backends/defaults: the normal model-loading
        path remains the fallback. Never download weights just to validate a hit.
        """
        import inspect
        from pathlib import Path

        from .perturbation_dataset_tool import _fingerprint

        if hasattr(_cp_models, "Cellpose"):
            return None
        selected = arguments.get("pretrained_model")
        if selected is None:
            parameter = inspect.signature(_cp_models.CellposeModel).parameters.get(
                "pretrained_model"
            )
            selected = parameter.default if parameter else None
        if not isinstance(selected, str):
            return None
        path = Path(selected).expanduser()
        if not path.is_file():
            if selected not in getattr(_cp_models, "MODEL_NAMES", []):
                return None
            directory = getattr(_cp_models, "MODEL_DIR", None)
            if directory is None:
                return None
            path = Path(directory) / selected
        if not path.is_file():
            return None
        path = path.resolve()
        return cls._model_details(path.name, _fingerprint(path), arguments, device)

    @staticmethod
    def _canonical_request(arguments, device):
        """Compare effective options, including old manifests lacking defaults."""
        import inspect
        from pathlib import Path

        ignored = {
            "resume",
            "preview_rows",
            "on_error",
            "output_dir",
            "run_id",
            "device",
            "max_run_seconds",
            "cancel_file",
            "resolved_device",
            "image_ref_stamp",
        }
        result = {key: value for key, value in arguments.items() if key not in ignored}
        defaults = {
            "model_type": "cyto3",
            "diameter": 0,
            "batch_size": 8,
            "flow_threshold": 0.4,
            "cellprob_threshold": 0.0,
            "max_pixels": 16777216,
            "mask_format": "tiff",
            "save_input": True,
            "save_preview": True,
            "save_flows": False,
            "save_styles": False,
            "z_index": 0,
            "t_index": 0,
        }
        for key, value in defaults.items():
            result.setdefault(key, value)
        result["diameter"] = result["diameter"] or 0
        for key in ("image_path", "image_ref"):
            if key in result:
                result[key] = str(Path(result[key]).expanduser().resolve())
        if "image_paths" in result:
            result["image_paths"] = [
                str(Path(p).expanduser().resolve()) for p in result["image_paths"]
            ]
        if not hasattr(_cp_models, "Cellpose"):
            # Legacy manifests could explicitly supply these ignored options.
            result.pop("model_type", None)
            result.pop("channels", None)
            selected = result.get("pretrained_model")
            if selected is None:
                parameter = inspect.signature(_cp_models.CellposeModel).parameters.get(
                    "pretrained_model"
                )
                selected = parameter.default if parameter else None
            if isinstance(selected, str):
                path = Path(selected).expanduser()
                if not path.is_file() and selected in getattr(
                    _cp_models, "MODEL_NAMES", []
                ):
                    directory = getattr(_cp_models, "MODEL_DIR", None)
                    if directory is not None:
                        path = Path(directory) / selected
                result["pretrained_model"] = (
                    str(path.resolve()) if path.is_file() else selected
                )
        else:
            result.setdefault("channels", [0, 0])
        result["resolved_device"] = device
        stamp = arguments.get("image_ref_stamp")
        result["image_ref_stamp"] = (
            {**stamp, "path": str(Path(stamp["path"]).expanduser().resolve())}
            if stamp
            else None
        )
        return result

    @staticmethod
    def _same_model_info(old, current):
        if not isinstance(old, dict):
            return False
        ignored = {"pretrained_model_requested", "ignored_parameters"}
        if not hasattr(_cp_models, "Cellpose"):
            ignored.add("model_type_requested")
        return {key: value for key, value in old.items() if key not in ignored} == {
            key: value for key, value in current.items() if key not in ignored
        }

    def run(self, arguments=None, stream_callback=None):
        from .cellpose_image_io import image_items, read_item

        try:
            import jsonschema

            arguments = {} if arguments is None else arguments
            jsonschema.validate(arguments, self.tool_config["parameter"])
            for key in (
                "diameter",
                "flow_threshold",
                "cellprob_threshold",
                "max_run_seconds",
            ):
                if key in arguments and (
                    isinstance(arguments[key], bool)
                    or not math.isfinite(arguments[key])
                ):
                    raise ValueError(f"{key} must be a finite number.")
            if not CELLPOSE_AVAILABLE:
                return {
                    "status": "error",
                    "error": f"Cellpose dependency unavailable: {_IMPORT_ERROR}",
                }
            if not hasattr(_cp_models, "Cellpose"):
                unsupported = [
                    key for key in ("model_type", "channels") if key in arguments
                ]
                if unsupported:
                    raise ValueError(
                        "Cellpose 4.x does not apply "
                        + ", ".join(unsupported)
                        + "; omit these legacy parameters. Use pretrained_model for weights "
                        "and channel_names/channel_indices for physical channels."
                    )
            device = self._resolve_device(arguments.get("device", "auto"))
            items, descriptor = image_items(arguments)
            if not 1 <= len(items) <= 512:
                raise ValueError(
                    "Select 1-512 logical images; narrow image_indices for larger references."
                )
            artifact_mode = "output_dir" in arguments
            if not artifact_mode and (
                "image_path" not in arguments
                or any(
                    key in arguments
                    for key in (
                        "save_flows",
                        "save_styles",
                        "mask_format",
                        "preview_rows",
                        "on_error",
                        "save_input",
                        "save_preview",
                        "run_id",
                        "resume",
                        "max_run_seconds",
                        "cancel_file",
                    )
                )
            ):
                raise ValueError(
                    "Batch/reference inputs and complete outputs require output_dir."
                )
            if artifact_mode and any(
                key in arguments for key in ("save_mask", "mask_output_path")
            ):
                raise ValueError(
                    "output_dir always saves masks; omit legacy save_mask/mask_output_path."
                )
            if arguments.get("resume") and not arguments.get("run_id"):
                raise ValueError("resume requires an explicit run_id.")
            if not artifact_mode:
                image, _source = read_item(items[0], arguments)
                backend = self._get_model(
                    arguments.get("model_type", "cyto3"),
                    device,
                    arguments.get("pretrained_model"),
                )
                native = self._infer(image, arguments, backend)
                objects = self._object_stats(native[0])
                data = {
                    "image_path": arguments["image_path"],
                    "model_type_requested": arguments.get("model_type", "cyto3"),
                    "model_used": backend[1],
                    "device": device,
                    "num_objects": len(objects),
                    "image_shape": list(native[0].shape),
                    "objects": objects,
                    "mask_path": None,
                }
                if not hasattr(_cp_models, "Cellpose"):
                    data["note"] = (
                        "Cellpose 4.x uses the selected/default pretrained model; legacy model_type/channels are not applied. Omitted diameter uses default scaling, not size estimation."
                    )
                if arguments.get("save_mask", False):
                    data["mask_path"], error = self._save_mask(
                        native[0],
                        arguments["image_path"],
                        arguments.get("mask_output_path"),
                    )
                    if error:
                        data["mask_save_error"] = error
                return {"status": "success", "data": data}
            return self._run_collection(
                arguments, items, descriptor, device, stream_callback
            )
        except jsonschema.ValidationError as exc:
            return {"status": "error", "error": exc.message}
        except Exception as exc:  # noqa: BLE001 - boundary for optional codecs/models
            return {"status": "error", "error": str(exc)}

    def _run_collection(
        self, arguments, items, descriptor, device, stream_callback=None
    ):
        import uuid
        from pathlib import Path

        from filelock import FileLock

        arguments = dict(arguments)
        arguments.setdefault("run_id", uuid.uuid4().hex)
        parent = Path(arguments["output_dir"]).expanduser().resolve()
        destination = parent / arguments["run_id"]
        for item in items:
            for path in item["paths"]:
                if (
                    "://" not in path
                    and Path(path).is_dir()
                    and destination.is_relative_to(Path(path).resolve())
                ):
                    raise ValueError("output_dir must be outside input image stores.")
        if (
            descriptor
            and descriptor["read_options"]["source_type"] == "directory"
            and destination.is_relative_to(Path(descriptor["read_options"]["source"]))
        ):
            raise ValueError(
                "output_dir must be outside the image_ref source directory."
            )
        parent.mkdir(parents=True, exist_ok=True)
        with FileLock(str(parent / ("." + arguments["run_id"] + ".lock")), timeout=0):
            return self._run_collection_locked(
                arguments, items, descriptor, device, stream_callback
            )

    def _run_collection_locked(
        self, arguments, items, descriptor, device, stream_callback=None
    ):
        import json
        import time
        import uuid
        from pathlib import Path

        from .cellpose_image_io import read_item
        from .cellpose_results import (
            atomic_json,
            atomic_records,
            export_result,
            pixel_hash,
            valid_result,
        )
        from .perturbation_dataset_tool import _fingerprint

        started = time.monotonic()
        identifier = arguments.get("run_id", uuid.uuid4().hex)
        directory = Path(arguments["output_dir"]).expanduser().resolve() / identifier
        manifest_path = directory / "run.json"
        records_path = directory / "results.jsonl"
        # Prevent outputs from becoming new members of an input image store/index.
        source_paths = [
            Path(p).resolve() for item in items for p in item["paths"] if "://" not in p
        ]
        if any(
            path.is_dir() and directory.is_relative_to(path) for path in source_paths
        ):
            raise ValueError("output_dir must be outside input image stores.")
        if (
            descriptor
            and descriptor["read_options"]["source_type"] == "directory"
            and directory.is_relative_to(Path(descriptor["read_options"]["source"]))
        ):
            raise ValueError(
                "output_dir must be outside the image_ref source directory."
            )
        request = {
            key: value
            for key, value in arguments.items()
            if key not in ("resume", "preview_rows", "on_error", "output_dir", "run_id")
        }
        request["resolved_device"] = device
        request["image_ref_stamp"] = (
            _fingerprint(Path(arguments["image_ref"]).expanduser().resolve())
            if descriptor
            else None
        )
        previous = {}
        if arguments.get("resume") and not directory.exists():
            raise ValueError("Cannot resume a missing run_id.")
        if directory.exists():
            if not arguments.get("resume"):
                raise ValueError(
                    "Run already exists; choose another run_id or resume=true."
                )
            manifest = json.loads(manifest_path.read_text())
            if (
                self._canonical_request(
                    manifest.get("request", {}),
                    manifest.get("request", {}).get("resolved_device", device),
                )
                != self._canonical_request(request, device)
                or manifest.get("schema_version") != 1
            ):
                raise ValueError("Resume request differs from the original run.")
            if records_path.exists():
                with records_path.open() as stream:
                    for line in stream:
                        record = json.loads(line)
                        old = (
                            record
                            if record.get("status") == "success"
                            else record.get("previous_success")
                        )
                        if old is not None:
                            index = old["item_index"]
                            if not isinstance(index, int) or not 0 <= index < len(
                                items
                            ):
                                raise ValueError(
                                    "Invalid item index in resume results."
                                )
                            previous[index] = old
        else:
            directory.mkdir(parents=True)
        # Keep every prior success in the new snapshot BEFORE processing starts.
        # Unvisited successes remain available even if this resume is interrupted
        # repeatedly, stops early, or temporarily cannot read its source.
        records = []
        for index, item in enumerate(items):
            row = {
                "item_index": index,
                "paths": item["paths"],
                "field_id": item.get("field_id"),
                "status": "not_processed",
                "reason": "Not yet processed in this invocation.",
            }
            if index in previous:
                row["previous_success"] = previous[index]
            records.append(row)
        backend, model_error, resume_info = None, None, None
        resume_info_checked = False
        manifest = {
            "schema_version": 1,
            "tool": "Cellpose_segment_image",
            "run_id": identifier,
            "request": request,
            "reference_integrity_scope": descriptor.get("integrity_scope")
            if descriptor
            else None,
            "status": "running",
            "n_requested": len(items),
            "n_processed": 0,
            "n_succeeded": 0,
            "n_failed": 0,
            "results_path": str(records_path),
            "resumed_successes": 0,
            "n_unprocessed": len(items),
            "phase": "ready",
        }
        atomic_records(records_path, records)
        atomic_json(manifest_path, manifest)

        def publish():
            manifest["elapsed_seconds"] = time.monotonic() - started
            atomic_json(manifest_path, manifest)
            if stream_callback is not None:
                try:
                    stream_callback(
                        json.dumps(
                            {
                                "event": "cellpose_progress",
                                "run_ref": str(manifest_path),
                                "results_path": str(records_path),
                                **manifest,
                            }
                        )
                    )
                except Exception as exc:  # noqa: BLE001 - observer failure must not abort inference
                    manifest["progress_warning"] = str(exc) or type(exc).__name__
                    atomic_json(manifest_path, manifest)

        def stop_reason():
            if (
                arguments.get("cancel_file")
                and Path(arguments["cancel_file"]).expanduser().exists()
            ):
                return "cancel_requested"
            budget = arguments.get("max_run_seconds")
            if budget is not None and time.monotonic() - started >= budget:
                return "time_budget_exhausted"
            return None

        def ensure_model():
            nonlocal backend, model_error
            if backend is None:
                manifest["phase"] = "loading_model"
                publish()
                try:
                    backend = self._get_model(
                        arguments.get("model_type", "cyto3"),
                        device,
                        arguments.get("pretrained_model"),
                    )
                except Exception as exc:
                    model_error = str(exc) or type(exc).__name__
                    raise
            return self._model_info(backend, arguments, device)

        publish()
        fatal = None
        try:
            for index, item in enumerate(items):
                reason = stop_reason()
                if reason:
                    manifest["stop_reason"] = reason
                    break
                start = time.monotonic()
                record = {
                    "item_index": index,
                    "paths": item["paths"],
                    "field_id": item.get("field_id"),
                }
                try:
                    manifest.update(phase="reading", current_item_index=index)
                    publish()
                    if model_error is not None:
                        raise ValueError(model_error)
                    image, source = read_item(item, arguments)
                    source["pixel_sha256"] = pixel_hash(image)
                    old = previous.get(index)
                    if old and not resume_info_checked:
                        resume_info = self._resume_model_info(arguments, device)
                        resume_info_checked = True
                    info = (
                        self._model_info(backend, arguments, device)
                        if backend
                        else resume_info
                    )
                    if old and info is None:
                        reason = stop_reason()
                        if reason:
                            manifest["stop_reason"] = reason
                            break
                        info = ensure_model()
                    if (
                        old
                        and info is not None
                        and old.get("status") == "success"
                        and old.get("input_pixel_sha256") == source["pixel_sha256"]
                        and self._same_model_info(old.get("model"), info)
                        and old.get("paths") == item["paths"]
                        and old.get("field_id") == item.get("field_id")
                        and valid_result(old, source)
                    ):
                        record = {**old, "resumed": True}
                        manifest["resumed_successes"] += 1
                    else:
                        reason = stop_reason()
                        if reason:
                            manifest["stop_reason"] = reason
                            break
                        info = ensure_model()
                        reason = stop_reason()
                        if reason:
                            manifest["stop_reason"] = reason
                            break
                        manifest["phase"] = "inference"
                        publish()
                        inference_start = time.monotonic()
                        native = self._infer(image, arguments, backend)
                        record["inference_seconds"] = time.monotonic() - inference_start
                        objects = self._object_stats(native[0])
                        attempt = (
                            directory / "items" / f"{index:05d}" / uuid.uuid4().hex
                        )
                        record["attempt_dir"] = str(attempt)
                        summary = export_result(
                            attempt, image, native, objects, source, info, arguments
                        )
                        record.update(status="success", model=info, **summary)
                except Exception as exc:  # noqa: BLE001 - one failed image does not hide others
                    record.update(status="error", error=str(exc) or type(exc).__name__)
                    if index in previous:
                        record["previous_success"] = previous[index]
                record["elapsed_seconds"] = time.monotonic() - start
                records[index] = record
                atomic_records(records_path, records)
                manifest.update(
                    n_processed=index + 1,
                    n_unprocessed=len(items) - index - 1,
                    n_succeeded=sum(r["status"] == "success" for r in records),
                    n_failed=sum(r["status"] == "error" for r in records),
                )
                manifest["phase"] = "image_complete"
                publish()
                if (
                    record["status"] == "error"
                    and arguments.get("on_error", "continue") == "stop"
                ):
                    manifest["stop_reason"] = "image_error"
                    break
            if manifest.get("stop_reason"):
                for row in records:
                    if row["status"] == "not_processed":
                        row["reason"] = manifest["stop_reason"]
                atomic_records(records_path, records)
            if descriptor:
                from .perturbation_image_tool import resolve_image_dataset_ref

                resolve_image_dataset_ref(arguments["image_ref"])
        except BaseException as exc:
            fatal = exc
            manifest["last_error"] = str(exc)
            if not isinstance(exc, Exception):
                raise
        finally:
            if fatal:
                manifest["status"] = (
                    "interrupted" if not isinstance(fatal, Exception) else "error"
                )
            elif manifest["n_failed"] or manifest["n_processed"] < len(items):
                manifest["status"] = (
                    "partial_success" if manifest["n_succeeded"] else "error"
                )
            else:
                manifest["status"] = "success"
            manifest["phase"] = "finished"
            publish()
        preview_limit = arguments.get("preview_rows", 20)
        response_records = []
        for row in records[:preview_limit]:
            preview = dict(row)
            if "objects_preview" in preview:
                preview["objects_preview"] = preview["objects_preview"][:preview_limit]
            if "previous_success" in preview:
                preview["previous_success"] = {
                    key: value
                    for key, value in preview["previous_success"].items()
                    if key in ("item_index", "status", "prediction_ref", "num_objects")
                }
            response_records.append(preview)
        data = {
            "run_ref": str(manifest_path),
            "results_path": str(records_path),
            "run_id": identifier,
            "n_requested": len(items),
            "n_processed": manifest["n_processed"],
            "n_succeeded": manifest["n_succeeded"],
            "n_failed": manifest["n_failed"],
            "resumed_successes": manifest["resumed_successes"],
            "n_unprocessed": len(items) - manifest["n_processed"],
            "results": response_records,
            "device": device,
        }
        if manifest.get("stop_reason"):
            data["stop_reason"] = manifest["stop_reason"]
        result = {"status": manifest["status"], "data": data}
        if result["status"] == "error":
            result["error"] = (
                manifest.get("last_error")
                or manifest.get("stop_reason")
                or (
                    "No image completed successfully; inspect results_path for per-image errors."
                )
            )
        return result
