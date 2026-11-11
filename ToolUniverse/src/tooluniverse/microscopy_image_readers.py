"""Metadata readers for single-image TIFF/OME-TIFF and local OME-Zarr groups."""

import xml.etree.ElementTree as ET
from pathlib import Path

from .perturbation_dataset_tool import DatasetLoadError


def _select(count, index, label, explicit=False):
    if count > 1 and not explicit and label == "series":
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "Multiple image series found; provide series_index explicitly.",
        )
    if not 0 <= index < count:
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            f"{label} index {index} out of range for {count} entries.",
        )
    return index


def ome_tiff_members(path):
    """Validate the entire OME UUID/file graph before tifffile assembles pixels."""
    import tifffile

    root_dir = path.parent.resolve()
    pending, checked = [(path.resolve(), None)], {}
    while pending:
        member, expected_uuid = pending.pop()
        if member in checked:
            if expected_uuid and checked[member] != expected_uuid:
                raise DatasetLoadError("INVALID_DATA", "OME UUID references disagree.")
            continue
        if not member.is_file():
            raise DatasetLoadError(
                "SOURCE_NOT_FOUND", f"Missing OME fileset member: {member}"
            )
        with tifffile.TiffFile(member) as image:
            xml = image.ome_metadata
        if not xml:
            if expected_uuid:
                raise DatasetLoadError(
                    "INVALID_DATA", "Referenced OME member has no OME metadata."
                )
            checked[member] = None
            continue
        root = ET.fromstring(xml)
        actual_uuid = root.get("UUID")
        if expected_uuid and actual_uuid != expected_uuid:
            raise DatasetLoadError("INVALID_DATA", f"OME UUID mismatch: {member.name}")
        checked[member] = actual_uuid
        for node in root.iter():
            if node.tag.rsplit("}", 1)[-1] != "UUID":
                continue
            reference_uuid = (node.text or "").strip() or None
            if reference_uuid and reference_uuid == actual_uuid:
                continue
            filename = node.get("FileName")
            if not filename:
                raise DatasetLoadError(
                    "INVALID_DATA", "External OME UUID requires an explicit FileName."
                )
            target = (member.parent / filename).resolve()
            if not target.is_relative_to(root_dir):
                raise DatasetLoadError(
                    "INVALID_DATA",
                    "OME fileset reference escapes its source directory.",
                )
            pending.append((target, reference_uuid))
    return sorted(checked)


def read_tiff(path, arguments):
    import tifffile

    members = ome_tiff_members(path)
    with tifffile.TiffFile(path) as image:
        root = ET.fromstring(image.ome_metadata) if image.is_ome else None
        index = _select(
            len(image.series),
            arguments.get("series_index", 0),
            "series",
            "series_index" in arguments,
        )
        base = image.series[index]
        level = _select(
            len(base.levels), arguments.get("resolution_level", 0), "resolution"
        )
        selected = base.levels[level]
        shape, axes = list(selected.shape), selected.axes
        if "S" in axes and "C" not in axes:
            axes = axes.replace("S", "C")
        metadata = {
            "series_index": index,
            "n_series": len(image.series),
            "resolution_level": level,
            "n_levels": len(base.levels),
            "is_bigtiff": image.is_bigtiff,
            "physical_sizes": {},
            "channel_metadata": [],
        }
        channels = None
        if root is not None:
            images = root.findall("{*}Image")
            if index >= len(images):
                raise DatasetLoadError(
                    "INVALID_DATA", "TIFF series does not map to an OME Image."
                )
            node = images[index]
            pixels = node.find("{*}Pixels")
            if pixels is None:
                raise DatasetLoadError(
                    "INVALID_DATA", "OME Image has no Pixels metadata."
                )
            metadata["image_name"] = node.get("Name")
            metadata["dimension_order"] = pixels.get("DimensionOrder")
            for axis in "XYZ":
                key = "PhysicalSize" + axis
                if pixels.get(key) is not None:
                    value = float(pixels.get(key))
                    if axis in axes and axis in base.axes:
                        value *= (
                            base.shape[base.axes.index(axis)] / shape[axes.index(axis)]
                        )
                    metadata["physical_sizes"][axis] = {
                        "value": value,
                        "unit": pixels.get(key + "Unit"),
                    }
            if pixels.get("TimeIncrement") is not None:
                metadata["time_increment"] = {
                    "value": float(pixels.get("TimeIncrement")),
                    "unit": pixels.get("TimeIncrementUnit"),
                }
            channel_nodes = pixels.findall("{*}Channel")
            metadata["channel_metadata"] = [
                dict(channel.attrib) for channel in channel_nodes
            ]
            count = shape[axes.index("C")] if "C" in axes else 1
            if len(channel_nodes) == count:
                channels = [
                    channel.get("Name") or channel.get("ID") or str(i)
                    for i, channel in enumerate(channel_nodes)
                ]
        if not arguments.get("lazy", True):
            selected.asarray()
        return {
            "format": "ome_tiff" if root is not None else "tiff",
            "shape": shape,
            "axes": axes,
            "dtype": str(selected.dtype),
            "native_channels": channels,
            "image_metadata": metadata,
            "related_files": [str(member) for member in members],
        }


def read_zarr(path, arguments):
    import zarr

    group_path = arguments.get("group_path", "")
    if isinstance(path, str):
        from .microscopy_remote import close_remote_fs, remote_group_url, remote_options

        url = remote_group_url(path, group_path)
        options = remote_options(
            url, arguments.get("remote_timeout", 30), arguments.get("remote_endpoint")
        )
        store = zarr.storage.FsspecStore.from_url(
            url,
            storage_options={**options, "skip_instance_cache": True},
            read_only=True,
        )
        try:
            # Initialize once before Zarr concurrently probes metadata keys.
            # Otherwise S3 backends can race and orphan several HTTP clients.
            from zarr.core.sync import sync

            if hasattr(store.fs, "set_session"):
                sync(store.fs.set_session())
            if url.startswith("s3://"):
                # Some s3fs versions treat a not-yet-used client's empty session
                # pool as closed. Complete one request before concurrent probes.
                try:
                    sync(
                        store.fs._cat_file(
                            store.path.rstrip("/") + "/.zgroup", start=0, end=1
                        )
                    )
                except FileNotFoundError:
                    pass  # Zarr 3 groups use zarr.json instead.
            root = zarr.open_group(store, mode="r", use_consolidated=False)
            return _read_zarr_group(root, arguments)
        finally:
            close_remote_fs(store.fs)
            store.close()
    selected_path = (path / group_path).resolve()
    if not selected_path.is_relative_to(path.resolve()):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "group_path must remain inside the local store."
        )
    root = zarr.open_group(str(selected_path), mode="r", use_consolidated=False)
    try:
        return _read_zarr_group(root, arguments, selected_path)
    finally:
        root.store.close()


def _read_zarr_group(root, arguments, selected_path=None):
    import numpy as np

    group_path = arguments.get("group_path", "")
    remote = selected_path is None
    attributes = dict(root.attrs)
    metadata = attributes.get("ome", attributes)
    if not isinstance(metadata, dict):
        raise DatasetLoadError("INVALID_DATA", "OME-Zarr attributes must be objects.")
    version = metadata.get("version")
    multiscales = metadata.get("multiscales")
    if not isinstance(multiscales, list) or not multiscales:
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT",
            "Select an OME-Zarr image group with multiscales, not a plate or well container.",
        )
    index = _select(
        len(multiscales),
        arguments.get("series_index", 0),
        "series",
        "series_index" in arguments,
    )
    multiscale = multiscales[index]
    version = version or multiscale.get("version")
    if version not in ("0.4", "0.5"):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT",
            f"Unsupported OME-Zarr metadata version: {version}; supported: 0.4, 0.5.",
        )
    if getattr(root.metadata, "zarr_format", None) != (2 if version == "0.4" else 3):
        raise DatasetLoadError(
            "INVALID_DATA", "NGFF version does not match its Zarr storage version."
        )
    if multiscale.get("version", version) != version:
        raise DatasetLoadError("INVALID_DATA", "Conflicting NGFF version declarations.")
    axes_metadata = multiscale.get("axes", [])
    if not axes_metadata or any(not isinstance(axis, dict) for axis in axes_metadata):
        raise DatasetLoadError(
            "INVALID_DATA", "OME-Zarr requires explicit axis objects."
        )
    datasets = multiscale.get("datasets", [])
    level = _select(len(datasets), arguments.get("resolution_level", 0), "resolution")
    dataset = datasets[level]
    from pathlib import PurePosixPath

    dataset_path = PurePosixPath(dataset["path"])
    if dataset_path.is_absolute() or ".." in dataset_path.parts:
        raise DatasetLoadError(
            "INVALID_DATA", "Multiscale dataset path escapes its image group."
        )
    if not remote and not (selected_path / dataset["path"]).resolve().is_relative_to(
        selected_path
    ):
        raise DatasetLoadError(
            "INVALID_DATA", "Multiscale dataset path escapes its image group."
        )
    array = root[dataset["path"]]
    axes = "".join(axis["name"].upper() for axis in axes_metadata)
    if len(axes_metadata) != array.ndim:
        raise DatasetLoadError(
            "INVALID_DATA", "OME-Zarr axes do not match array dimensionality."
        )
    if version == "0.5" and tuple(
        getattr(array.metadata, "dimension_names", ()) or ()
    ) != tuple(axis["name"] for axis in axes_metadata):
        raise DatasetLoadError(
            "INVALID_DATA", "Zarr 3 dimension_names must match NGFF axes."
        )
    transforms = dataset.get("coordinateTransformations", [])
    if sum(t.get("type") == "scale" for t in transforms) != 1:
        raise DatasetLoadError(
            "INVALID_DATA", "Each NGFF level must declare exactly one scale transform."
        )
    global_transforms = multiscale.get("coordinateTransformations", [])
    scale, translation = np.ones(array.ndim), np.zeros(array.ndim)
    for transform in transforms + global_transforms:
        kind = transform.get("type")
        if kind not in ("scale", "translation"):
            raise DatasetLoadError(
                "UNSUPPORTED_FORMAT", f"Unsupported coordinate transformation: {kind}"
            )
        values = np.asarray(transform[kind], dtype=float)
        if values.shape != (array.ndim,) or not np.isfinite(values).all():
            raise DatasetLoadError(
                "INVALID_DATA", "Invalid coordinate transform dimensions or values."
            )
        if kind == "scale":
            if (values <= 0).any():
                raise DatasetLoadError("INVALID_DATA", "Scale values must be positive.")
            scale *= values
            translation *= values
        else:
            translation += values
    omero = metadata.get("omero", {})
    channel_metadata = omero.get("channels", [])
    channel_count = array.shape[axes.index("C")] if "C" in axes else 1
    channels = None
    if channel_metadata:
        if len(channel_metadata) != channel_count:
            raise DatasetLoadError(
                "INVALID_DATA", "OME-Zarr channel metadata count does not match array."
            )
        channels = [
            channel.get("label") or str(i) for i, channel in enumerate(channel_metadata)
        ]
    if not arguments.get("lazy", True):
        # Decode a chunk at a time rather than materializing a whole image pyramid.
        import itertools

        chunks = array.chunks
        for starts in itertools.product(
            *(
                range(0, size, chunk)
                for size, chunk in zip(array.shape, chunks, strict=True)
            )
        ):
            array[
                tuple(
                    slice(start, min(start + chunk, size))
                    for start, chunk, size in zip(
                        starts, chunks, array.shape, strict=True
                    )
                )
            ]
    return {
        "format": "ome_zarr",
        "shape": list(array.shape),
        "axes": axes,
        "dtype": str(array.dtype),
        "native_channels": channels,
        "image_metadata": {
            "ngff_version": version,
            "group_path": group_path,
            "array_path": dataset["path"],
            "series_index": index,
            "n_series": len(multiscales),
            "resolution_level": level,
            "n_levels": len(datasets),
            "axes_metadata": axes_metadata,
            "scale": scale.tolist(),
            "translation": translation.tolist(),
            "coordinate_transformations": transforms,
            "global_coordinate_transformations": global_transforms,
            "channel_metadata": channel_metadata,
        },
    }


def read_native(path, arguments):
    """Read vendor metadata with optional backends, retaining explicit selections."""
    import errno
    import shutil
    import tempfile

    if arguments.get("group_path"):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "group_path applies only to OME-Zarr."
        )
    fmt = path.suffix.lower()[1:]
    staged = False
    temporary = None
    image = None
    try:
        if fmt == "nd2":
            import nd2

            try:
                image = nd2.ND2File(path)
            except OSError as exc:
                if exc.errno not in (errno.ENODEV, errno.ENOSYS, errno.ENOTSUP):
                    raise
                temporary = tempfile.TemporaryDirectory(prefix="beam-nd2-")
                local = Path(temporary.name) / path.name
                shutil.copyfile(path, local)
                image = nd2.ND2File(local)
                staged = True
            sizes = dict(image.sizes)
            count = sizes.pop("P", 1)
            index = _select(
                count,
                arguments.get("series_index", 0),
                "series",
                "series_index" in arguments,
            )
            if arguments.get("resolution_level", 0) != 0:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "ND2 has no selectable pyramid level."
                )
            if "position_index" in arguments:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "Use series_index for ND2 positions."
                )
            axes = "".join(sizes)
            if "S" in axes and "C" not in axes:
                axes = axes.replace("S", "C")
            shape = list(sizes.values())
            channel_nodes = image.metadata.channels or []
            names = [channel.channel.name for channel in channel_nodes]
            voxel = image.voxel_size()
            calibrated = (
                channel_nodes[0].volume.axesCalibrated
                if channel_nodes
                else (False, False, False)
            )
            metadata = {
                "physical_sizes": {
                    axis: {"value": float(value), "unit": "micrometer"}
                    for axis, value, known in zip("XYZ", voxel, calibrated, strict=True)
                    if axis in axes and known
                },
                "series_index": index,
                "n_series": count,
                "resolution_level": 0,
                "n_levels": 1,
                "staged_local_copy": staged,
                "channel_names": names,
                "is_rgb": bool(image.is_rgb),
            }
            dtype = str(image.dtype)
            if not arguments.get("lazy", True):
                image.asarray(position=index if "P" in image.sizes else None)
        elif fmt == "czi":
            import czifile

            if not hasattr(czifile.CziFile, "scenes"):
                raise ImportError(
                    "CZI loading requires czifile>=2026.11.8 and Python>=3.12."
                )
            image = czifile.CziFile(path)
            scenes = list(image.scenes.items())
            index = _select(
                len(scenes),
                arguments.get("series_index", 0),
                "series",
                "series_index" in arguments,
            )
            scene_key, scene = scenes[index]
            level = _select(
                len(scene.levels), arguments.get("resolution_level", 0), "resolution"
            )
            selected = scene.levels[level]
            if "position_index" in arguments:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "Use series_index for CZI scenes."
                )
            axes, shape = selected.axes, list(selected.shape)
            if "S" in axes and "C" not in axes:
                axes = axes.replace("S", "C")
            names = list(selected.channels)
            dtype = str(selected.dtype)
            metadata = {
                "series_index": index,
                "scene_coordinate": int(scene_key),
                "n_series": len(scenes),
                "resolution_level": level,
                "n_levels": len(scene.levels),
                "physical_sizes": {
                    axis: {
                        "value": float(value),
                        "unit": selected.coord_units.get(axis),
                    }
                    for axis, value in selected.coord_scales.items()
                },
                "coordinate_offsets": {
                    axis: float(value) for axis, value in selected.coord_offsets.items()
                },
                "channel_metadata": selected.channels,
            }
            if not arguments.get("lazy", True):
                selected.asarray()
        elif fmt == "lif":
            from readlif.reader import LifFile

            # readlif closes its metadata handle and opens/closes each decoded plane.
            container = LifFile(str(path))
            index = _select(
                container.num_images,
                arguments.get("series_index", 0),
                "series",
                "series_index" in arguments,
            )
            selected = container.get_image(index)
            if arguments.get("resolution_level", 0) != 0:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT",
                    "LIF has no selectable pyramid level in this backend.",
                )
            if selected.display_dims != (1, 2) or any(
                k not in (1, 2, 3, 4, 10) and v > 1 for k, v in selected.dims_n.items()
            ):
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT",
                    "Unsupported non-XY or extra-dimensional LIF acquisition.",
                )
            position = _select(
                selected.n_mosaic, arguments.get("position_index", 0), "position"
            )
            if selected.n_mosaic > 1 and "position_index" not in arguments:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "Provide position_index for a LIF mosaic."
                )
            if any(b < 8 or b > 16 for b in selected.bit_depth):
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT",
                    "LIF backend supports 8-to-16-bit unsigned samples only.",
                )
            depths = {8 if b == 8 else 16 for b in selected.bit_depth}
            if len(depths) != 1 or max(depths) > 16:
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT",
                    "LIF backend supports uniform 8/16-bit channel storage.",
                )
            if selected.offsets[1] == 0 or sum(selected.offsets) > path.stat().st_size:
                raise DatasetLoadError(
                    "INVALID_DATA", "LIF image block is truncated or absent."
                )
            axes = "TCZYX"
            shape = [
                selected.nt,
                selected.channels,
                selected.nz,
                selected.dims.y,
                selected.dims.x,
            ]
            dtype = "uint8" if next(iter(depths)) == 8 else "uint16"
            names = []
            metadata = {
                "series_index": index,
                "n_series": container.num_images,
                "resolution_level": 0,
                "n_levels": 1,
                "position_index": position,
                "n_positions": selected.n_mosaic,
                "image_name": selected.name,
                "channel_bit_depths": list(selected.bit_depth),
                "physical_sizes": {
                    axis: {"value": 1.0 / value, "unit": "micrometer"}
                    for axis, value in zip("XYZ", selected.scale[:3], strict=True)
                    if value > 0
                },
            }
            if selected.scale[3] > 0:
                metadata["time_increment"] = {
                    "value": 1.0 / selected.scale[3],
                    "unit": "second",
                }
            if not arguments.get("lazy", True):
                for t in range(selected.nt):
                    for c in range(selected.channels):
                        for z in range(selected.nz):
                            selected.get_frame(t=t, c=c, z=z, m=position)
        else:
            raise DatasetLoadError(
                "UNSUPPORTED_FORMAT", f"Unsupported vendor format: {fmt}"
            )
        count = shape[axes.index("C")] if "C" in axes else 1
        return {
            "format": fmt,
            "shape": shape,
            "axes": axes,
            "dtype": dtype,
            "native_channels": names if len(names) == count else None,
            "image_metadata": metadata,
        }
    finally:
        if image is not None:
            image.close()
        if temporary is not None:
            temporary.cleanup()
