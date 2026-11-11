"""Anonymous HTTP(S)/S3 OME-Zarr metadata snapshots; never persist credentials."""

import hashlib
import json
from pathlib import PurePosixPath
from urllib.parse import urlsplit

from .perturbation_dataset_tool import DatasetLoadError


def remote_group_url(source, group_path=""):
    url = urlsplit(source)
    if url.scheme not in ("http", "https", "s3") or not url.netloc:
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "Remote sources require an http, https or s3 URL."
        )
    if url.username or url.password or url.query or url.fragment:
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "Use anonymous URLs without embedded credentials, query or fragment.",
        )
    if (
        PurePosixPath(group_path).is_absolute()
        or ".." in PurePosixPath(group_path).parts
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "group_path must stay within its remote store."
        )
    return source.rstrip("/") + ("/" + group_path.strip("/") if group_path else "")


def remote_options(source, timeout=30, endpoint=None):
    if endpoint:
        checked = urlsplit(endpoint)
        if (
            checked.scheme not in ("http", "https")
            or not checked.netloc
            or checked.username
            or checked.password
            or checked.query
            or checked.fragment
        ):
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "Invalid anonymous S3 endpoint URL."
            )
        if urlsplit(source).scheme != "s3":
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "remote_endpoint applies only to s3 URLs."
            )
    if urlsplit(source).scheme == "s3":
        options = {
            "anon": True,
            "config_kwargs": {"connect_timeout": timeout, "read_timeout": timeout},
        }
        if endpoint:
            options["client_kwargs"] = {"endpoint_url": endpoint}
        return options
    import aiohttp

    return {
        "client_kwargs": {
            "timeout": aiohttp.ClientTimeout(total=timeout),
            "trust_env": True,
        }
    }


def remote_metadata_stamp(source, arguments):
    import fsspec

    group = remote_group_url(source, arguments.get("group_path", ""))
    options = remote_options(
        group, arguments.get("remote_timeout", 30), arguments.get("remote_endpoint")
    )
    if urlsplit(group).scheme == "s3":
        fs, prefix = fsspec.core.url_to_fs(group, **options, skip_instance_cache=True)
    else:
        fs, prefix = None, group
    try:
        if fs is not None and hasattr(fs, "set_session"):
            from fsspec.asyn import sync

            sync(fs.loop, fs.set_session)
        return _collect_remote_metadata(source, arguments, fs, prefix)
    finally:
        if fs is not None:
            close_remote_fs(fs)


def _collect_remote_metadata(source, arguments, fs, prefix):
    documents = {}

    def load(relative, required=True):
        target = prefix.rstrip("/") + "/" + relative
        try:
            if fs is None:
                import requests

                timeout = arguments.get("remote_timeout", 30)
                with requests.get(
                    target, stream=True, timeout=(timeout, timeout)
                ) as response:
                    if response.status_code == 404:
                        raise FileNotFoundError(target)
                    response.raise_for_status()
                    pieces, length = [], 0
                    for block in response.iter_content(64 * 1024):
                        length += len(block)
                        if length > 4 * 1024**2:
                            raise DatasetLoadError(
                                "INVALID_DATA", "Remote metadata object exceeds 4 MiB."
                            )
                        pieces.append(block)
                    raw = b"".join(pieces)
            else:
                with fs.open(target, "rb") as stream:
                    raw = stream.read(4 * 1024**2 + 1)
        except FileNotFoundError:
            if required:
                raise DatasetLoadError(
                    "SOURCE_NOT_FOUND", f"Remote metadata missing: {relative}"
                ) from None
            return None
        if len(raw) > 4 * 1024**2:
            raise DatasetLoadError(
                "INVALID_DATA", "Remote metadata object exceeds 4 MiB."
            )
        documents[relative] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    group_v3 = load("zarr.json", False)
    if group_v3 is not None:
        metadata = group_v3.get("attributes", {}).get("ome", {})
    else:
        load(".zgroup")
        metadata = load(".zattrs")
    scales = metadata.get("multiscales", [])
    index = arguments.get("series_index", 0)
    if not 0 <= index < len(scales):
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "Remote multiscales index is out of range or image group is absent.",
        )
    levels = scales[index].get("datasets", [])
    level = arguments.get("resolution_level", 0)
    if not 0 <= level < len(levels):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "Remote resolution index out of range."
        )
    path = levels[level]["path"]
    if PurePosixPath(path).is_absolute() or ".." in PurePosixPath(path).parts:
        raise DatasetLoadError("INVALID_DATA", "Remote dataset path escapes its group.")
    if group_v3 is not None:
        load(path.rstrip("/") + "/zarr.json")
    else:
        load(path.rstrip("/") + "/.zarray")
        load(path.rstrip("/") + "/.zattrs", False)
    digest = hashlib.sha256(json.dumps(documents, sort_keys=True).encode()).hexdigest()
    return {
        "path": source,
        "kind": "remote",
        "integrity_scope": "metadata_only",
        "group_path": arguments.get("group_path", ""),
        "series_index": index,
        "resolution_level": level,
        "remote_timeout": arguments.get("remote_timeout", 30),
        "remote_endpoint": arguments.get("remote_endpoint"),
        "metadata_objects": documents,
        "sha256": digest,
    }


def close_remote_fs(fs):
    """Close backend sessions on the event loop that created them."""

    async def close():
        creator = getattr(fs, "_s3creator", None)
        if creator is not None:
            await creator.__aexit__(None, None, None)
            fs._s3creator = None
            fs._s3 = None
        session = getattr(fs, "_session", None)
        if session is not None and not session.closed:
            await session.close()

    if getattr(fs, "asynchronous", False):
        from zarr.core.sync import sync

        sync(close())
    else:
        from fsspec.asyn import sync

        sync(fs.loop, close, timeout=30)
