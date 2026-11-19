"""Read explicit 2D planes from the microscopy loader's persistent references."""

import errno
import json
import math
import shutil
import tempfile
from pathlib import Path

from .perturbation_dataset_tool import DatasetLoadError


def image_items(arguments):
    """Return logical images; wide-manifest channels sharing field_id stay together."""
    supplied = [
        key for key in ("image_path", "image_paths", "image_ref") if key in arguments
    ]
    if len(supplied) != 1:
        raise ValueError("Provide exactly one of image_path, image_paths or image_ref.")
    if "image_ref" not in arguments:
        if "image_indices" in arguments:
            raise ValueError("image_indices applies only to image_ref.")
        paths = arguments.get("image_paths", [arguments.get("image_path")])
        return [
            {"paths": [str(Path(path).expanduser().resolve())], "field_id": None}
            for path in paths
        ], None
    from .perturbation_image_tool import resolve_image_dataset_ref

    descriptor = resolve_image_dataset_ref(arguments["image_ref"])
    if Path(descriptor["index"]["path"]).stat().st_size > 128 * 1024**2:
        raise ValueError(
            "Image index exceeds 128 MiB; supply a smaller image reference."
        )
    groups, ordered = {}, []
    with Path(descriptor["index"]["path"]).open() as stream:
        for index, line in enumerate(stream):
            if index >= 100000:
                raise ValueError("Image index exceeds 100000 records.")
            record = json.loads(line)
            field = record["field_id"]
            key = ("field", field) if field is not None else ("record", index)
            if key not in groups:
                item = {
                    "paths": [],
                    "records": [],
                    "field_id": field,
                    "logical_index": len(ordered),
                    "read_options": descriptor["read_options"],
                }
                groups[key] = item
                ordered.append(item)
            groups[key]["paths"].append(record["path"])
            groups[key]["records"].append(record)
    selection = arguments.get("image_indices", list(range(len(ordered))))
    if any(index >= len(ordered) for index in selection):
        raise ValueError("image_indices contains an out-of-range logical image index.")
    return [ordered[index] for index in selection], descriptor


def _plane_selection(record, arguments):
    axes, shape = record["axes"], record["shape"]
    if not set(axes) <= set("TCZYX") or len(axes) != len(shape):
        raise ValueError(f"Unsupported source axes: {axes}")
    selection, remaining = [], []
    roi = arguments.get("roi")
    for axis, size in zip(axes, shape, strict=True):
        if axis in "TZ":
            key = "t_index" if axis == "T" else "z_index"
            if size > 1 and key not in arguments:
                raise ValueError(
                    f"2D inference requires a plane selection: source has {size} {axis} planes; specify {key}."
                )
            value = arguments.get(key, 0)
            if value >= size:
                raise ValueError(f"{key} is out of range.")
            selection.append(value)
        elif axis == "C":
            selection.append(0)  # replaced separately for each selected channel
        else:
            start = roi[axis.lower()] if roi else 0
            length = roi["height" if axis == "Y" else "width"] if roi else size
            if start + length > size:
                raise ValueError("ROI exceeds source image bounds.")
            selection.append(slice(start, start + length))
            remaining.append(axis)
    for axis, key in [("T", "t_index"), ("Z", "z_index")]:
        if axis not in axes and arguments.get(key, 0) != 0:
            raise ValueError(f"{key} cannot select an absent axis.")
    return selection, remaining


def _slice(array, record, arguments, channels):
    import numpy as np

    selection, remaining = _plane_selection(record, arguments)
    planes = []
    for channel in channels:
        if "C" in record["axes"]:
            selection[record["axes"].index("C")] = channel
        plane = np.asarray(array[tuple(selection)])
        if plane.ndim != 2:
            raise ValueError("Plane reader did not produce two spatial axes.")
        planes.append(plane.transpose([remaining.index("Y"), remaining.index("X")]))
    return planes


def _read_record(record, options, arguments, channels):
    import numpy as np

    fmt, source = record["format"], record["path"]
    meta = record["image_metadata"]
    if fmt in ("tiff", "ome_tiff"):
        import tifffile

        with tifffile.TiffFile(source) as handle:
            series = handle.series[meta["series_index"]].levels[
                meta["resolution_level"]
            ]
            try:
                import zarr
            except ImportError:
                if (
                    math.prod(record["shape"]) * np.dtype(record["dtype"]).itemsize
                    > 256 * 1024**2
                ):
                    raise ImportError(
                        "Large TIFF plane slicing requires zarr."
                    ) from None
                return _slice(series.asarray(), record, arguments, channels)
            store = series.aszarr()
            try:
                array = zarr.open(store, mode="r")
                return _slice(array, record, arguments, channels)
            finally:
                store.close()
    if fmt in ("png", "jpeg", "bmp"):
        from PIL import Image

        with Image.open(source) as image:
            return _slice(np.asarray(image), record, arguments, channels)
    if fmt == "ome_zarr":
        import zarr

        from .microscopy_remote import close_remote_fs, remote_group_url, remote_options

        remote = "://" in source
        if remote:
            url = remote_group_url(source, meta["group_path"])
            store = zarr.storage.FsspecStore.from_url(
                url,
                storage_options={
                    **remote_options(
                        url,
                        options.get("remote_timeout", 30),
                        options.get("remote_endpoint"),
                    ),
                    "skip_instance_cache": True,
                },
                read_only=True,
            )
        else:
            store = zarr.storage.LocalStore(
                str(Path(source) / meta["group_path"]), read_only=True
            )
        try:
            if remote and hasattr(store.fs, "set_session"):
                from zarr.core.sync import sync

                sync(store.fs.set_session())
                if source.startswith("s3://"):
                    try:
                        sync(
                            store.fs._cat_file(
                                store.path.rstrip("/") + "/.zgroup", start=0, end=1
                            )
                        )
                    except FileNotFoundError:
                        pass
            group = zarr.open_group(store, mode="r", use_consolidated=False)
            return _slice(group[meta["array_path"]], record, arguments, channels)
        finally:
            if remote:
                close_remote_fs(store.fs)
            store.close()
    if fmt == "lif":
        from readlif.reader import LifFile

        image = LifFile(source).get_image(meta["series_index"])
        selection, _remaining = _plane_selection(record, arguments)
        # Plane decoding in readlif does not load other time points or channels.
        spatial = {
            axis: value
            for axis, value in zip(record["axes"], selection, strict=True)
            if axis in "YX"
        }
        return [
            np.asarray(
                image.get_frame(
                    t=arguments.get("t_index", 0),
                    z=arguments.get("z_index", 0),
                    c=channel,
                    m=meta.get("position_index", 0),
                )
            )[spatial["Y"], spatial["X"]]
            for channel in channels
        ]
    # Vendor backends lacking a plane-slicing API are deliberately bounded.
    if math.prod(record["shape"]) * np.dtype(record["dtype"]).itemsize > 256 * 1024**2:
        raise ValueError(
            "Vendor acquisition exceeds 256 MiB decoder limit; export the selected plane to TIFF."
        )
    if fmt == "czi":
        import czifile

        with czifile.CziFile(source) as handle:
            scene = list(handle.scenes.values())[meta["series_index"]]
            return _slice(
                scene.levels[meta["resolution_level"]].asarray(),
                record,
                arguments,
                channels,
            )
    if fmt == "nd2":
        import nd2

        temporary = None
        try:
            try:
                handle = nd2.ND2File(source)
            except OSError as exc:
                if exc.errno not in (errno.ENODEV, errno.ENOSYS, errno.ENOTSUP):
                    raise
                temporary = tempfile.TemporaryDirectory(prefix="beam-cellpose-nd2-")
                staged = Path(temporary.name) / Path(source).name
                shutil.copyfile(source, staged)
                handle = nd2.ND2File(staged)
            with handle:
                array = handle.asarray(
                    position=meta["series_index"] if "P" in handle.sizes else None
                )
                return _slice(array, record, arguments, channels)
        finally:
            if temporary:
                temporary.cleanup()
    raise ValueError(f"Unsupported image reference format: {fmt}")


def read_item(item, arguments):
    """Return HWC/YX real pixels and explicit selection metadata, never inferred Z/T."""
    import numpy as np

    from .perturbation_image_tool import _image_header

    records = item.get("records")
    if records is None:
        path = Path(item["paths"][0])
        if not path.exists():
            raise FileNotFoundError(f"Image not found: {path}")
        header = _image_header(path, {"lazy": True})
        records = [{"path": str(path), **header}]
    channels = [
        (record, channel, name)
        for record in records
        for channel, name in enumerate(record["channels"])
    ]
    if "channel_indices" in arguments and "channel_names" in arguments:
        raise ValueError("Choose channel_indices or channel_names, not both.")
    if "channel_names" in arguments:
        names = [channel[2] for channel in channels]
        if any(name not in names for name in arguments["channel_names"]):
            raise ValueError(
                "A requested channel name is absent from the logical image."
            )
        chosen = [names.index(name) for name in arguments["channel_names"]]
    else:
        chosen = arguments.get("channel_indices", list(range(len(channels))))
    if not 1 <= len(chosen) <= 3 or any(i >= len(channels) for i in chosen):
        raise ValueError(
            "Select 1-3 valid channel_indices; channels are never silently discarded."
        )
    shape = records[0]["spatial_shape"]
    roi = arguments.get("roi")
    yx = [roi["height"], roi["width"]] if roi else shape
    if yx[0] * yx[1] * len(chosen) > arguments.get("max_pixels", 16777216):
        raise ValueError(
            "Selected pixel samples exceed max_pixels; choose a smaller ROI/level."
        )
    planes = []
    for index in chosen:
        record, channel, _name = channels[index]
        planes.extend(
            _read_record(record, item.get("read_options", {}), arguments, [channel])
        )
    if len({plane.shape for plane in planes}) != 1:
        raise ValueError("Channel planes have different spatial shapes.")
    image = planes[0] if len(planes) == 1 else np.stack(planes, axis=-1)
    if not np.issubdtype(image.dtype, np.number) or not np.isfinite(image).all():
        raise ValueError("Selected pixels must be finite numeric values.")
    selection = {
        key: arguments[key] for key in ("t_index", "z_index", "roi") if key in arguments
    }
    return image, {
        "paths": item["paths"],
        "field_id": item.get("field_id"),
        "logical_index": item.get("logical_index"),
        "channel_indices": chosen,
        "channels": [channels[i][2] for i in chosen],
        "selection": selection,
        "origin_yx": [roi["y"], roi["x"]] if roi else [0, 0],
        "records": records,
    }
