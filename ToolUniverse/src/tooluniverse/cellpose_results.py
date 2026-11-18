"""Persistent, lossless Cellpose arrays and bounded API summaries."""

import hashlib
import json
from pathlib import Path

from .perturbation_dataset_tool import _fingerprint


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def atomic_records(path, rows):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w") as stream:
        for row in rows:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
    temporary.replace(path)


def pixel_hash(array):
    import numpy as np

    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256(json.dumps([list(array.shape), str(array.dtype)]).encode())
    view = memoryview(contiguous).cast("B")
    for start in range(0, len(view), 4 * 1024**2):
        digest.update(view[start : start + 4 * 1024**2])
    return digest.hexdigest()


def _save_array(path, array):
    import numpy as np

    array = np.asarray(array)
    if not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all():
        raise ValueError("Backend output arrays must be finite and numeric.")
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("wb") as stream:
        np.save(stream, array, allow_pickle=False)
    temporary.replace(path)
    return {**_fingerprint(path), "shape": list(array.shape), "dtype": str(array.dtype)}


def _preview(path, image, masks):
    import numpy as np
    from PIL import Image

    # A display-only percentile stretch; input.npy retains the original values.
    stride = max(1, (max(masks.shape) + 1023) // 1024)
    samples = image[::stride, ::stride]
    labels = masks[::stride, ::stride]
    if samples.ndim == 2:
        samples = samples[..., None]
    rgb = np.zeros((*samples.shape[:2], 3), dtype=np.uint8)
    for index in range(samples.shape[-1]):
        plane = samples[..., index].astype(np.float32)
        low, high = np.percentile(plane, [1, 99])
        scaled = np.clip(
            (plane - low) / max(float(high - low), 1e-12) * 255, 0, 255
        ).astype(np.uint8)
        if samples.shape[-1] == 1:
            rgb[:] = scaled[..., None]
        else:
            rgb[..., index] = scaled
    boundary = np.zeros(labels.shape, dtype=bool)
    boundary[1:] |= labels[1:] != labels[:-1]
    boundary[:, 1:] |= labels[:, 1:] != labels[:, :-1]
    boundary &= labels > 0
    rgb[boundary] = [255, 60, 0]
    Image.fromarray(rgb).save(path)
    return {
        **_fingerprint(path),
        "display_only": True,
        "downsample_stride": stride,
        "contrast_percentiles": [1, 99],
    }


def export_result(directory, image, native, objects, source, model, arguments):
    import numpy as np
    import tifffile
    from PIL import Image

    masks = np.asarray(native[0])
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    fmt = arguments.get("mask_format", "tiff")
    if int(masks.max(initial=0)) > np.iinfo(np.uint32).max:
        raise ValueError("Mask labels exceed uint32 storage capacity.")
    labels = masks.astype(np.uint32)
    mask_path = (
        directory / {"tiff": "mask.tif", "png": "mask.png", "npy": "mask.npy"}[fmt]
    )
    if fmt == "tiff":
        tifffile.imwrite(mask_path, labels, metadata={"axes": "YX"})
    elif fmt == "npy":
        with mask_path.open("wb") as stream:
            np.save(stream, labels, allow_pickle=False)
    else:
        if int(labels.max(initial=0)) > 65535:
            raise ValueError("PNG labels exceed 65535; choose tiff or npy.")
        Image.fromarray(labels.astype(np.uint16)).save(mask_path)
    object_path = directory / "objects.jsonl"
    atomic_records(object_path, objects)
    artifacts = {
        "mask": {
            **_fingerprint(mask_path),
            "shape": list(labels.shape),
            "dtype": "uint16" if fmt == "png" else "uint32",
            "axes": "YX",
            "background_label": 0,
        },
        "objects": _fingerprint(object_path),
    }
    if arguments.get("save_input", True):
        artifacts["input"] = _save_array(directory / "input.npy", image)
    if arguments.get("save_preview", True):
        artifacts["preview"] = _preview(directory / "preview.png", image, labels)
    if arguments.get("save_flows", False):
        names = ["flow_rgb", "flow_vector", "cell_probability", "pixel_positions"]
        artifacts["flows"] = []
        for index, array in enumerate(native[1]):
            name = names[index] if index < len(names) else f"native_flow_{index}"
            stamp = _save_array(directory / (name + ".npy"), array)
            artifacts["flows"].append({"native_index": index, "name": name, **stamp})
    if arguments.get("save_styles", False):
        if len(native) < 3 or native[2] is None:
            raise ValueError("This backend did not return styles.")
        artifacts["styles"] = _save_array(directory / "native_styles.npy", native[2])
    data = {
        "schema_version": 1,
        "source": source,
        "model": model,
        "image_shape": list(labels.shape),
        "num_objects": len(objects),
        "objects_preview": objects[: arguments.get("preview_rows", 20)],
        "object_coordinate_frame": "selected ROI; add source.origin_yx for source-plane coordinates",
        "artifacts": artifacts,
    }
    if len(native) > 3 and native[3] is not None:
        artifacts["native_diameters"] = _save_array(
            directory / "native_diameters.npy", native[3]
        )
    result_path = directory / "result.json"
    atomic_json(result_path, data)
    return {
        "prediction_ref": str(result_path),
        "num_objects": len(objects),
        "image_shape": list(labels.shape),
        "mask_path": str(mask_path),
        "objects_path": str(object_path),
        "objects_preview": data["objects_preview"],
        "artifacts": artifacts,
        "model_used": model["model_used"],
        "device": model["device"],
        "input_pixel_sha256": source["pixel_sha256"],
        "result_stamp": _fingerprint(result_path),
    }


def valid_result(result, source=None):
    """Resume only immutable successes, including every optional output file."""
    try:
        if _fingerprint(Path(result["result_stamp"]["path"])) != result["result_stamp"]:
            return False
        saved = json.loads(Path(result["result_stamp"]["path"]).read_text())
        if source is not None:

            def comparable(value):
                value = dict(value)
                selection = dict(value.get("selection", {}))
                # Explicit zero selects the same plane as an omitted singleton
                # axis. Non-singleton axes are already required by the reader.
                for key in ("t_index", "z_index"):
                    if selection.get(key) == 0:
                        selection.pop(key)
                if "selection" in value:
                    value["selection"] = selection
                return value

            if comparable(saved["source"]) != comparable(source):
                return False
        if (
            saved["num_objects"] != result["num_objects"]
            or saved["image_shape"] != result["image_shape"]
            or saved["artifacts"] != result["artifacts"]
            or saved["model"] != result["model"]
            or saved["source"]["pixel_sha256"] != result["input_pixel_sha256"]
        ):
            return False
        stamps = []
        for value in result["artifacts"].values():
            stamps.extend(value if isinstance(value, list) else [value])
        for stamp in stamps:
            current = _fingerprint(Path(stamp["path"]))
            if any(current[key] != stamp[key] for key in current):
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False
