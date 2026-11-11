"""Exercise the local spatial transcriptomics tools through the SDK and MCP."""

import asyncio
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from tooluniverse import ToolUniverse

CATEGORY = "spatial_transcriptomics"
TOOL_NAMES = {"cell2location", "Tangram", "GraphST"}
pytestmark = pytest.mark.integration


@pytest.fixture
def universe(tmp_path):
    tu = ToolUniverse(workspace=str(tmp_path / "workspace"), load_workspace=False)
    try:
        tu.load_tools(categories=[CATEGORY], clear_existing=True)
        yield tu
    finally:
        tu.close()


@pytest.fixture
def expression_data(tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(10)
    genes = [f"gene{i}" for i in range(90)]
    profiles = rng.uniform(1, 8, size=(3, len(genes)))
    for i in range(3):
        profiles[i, i * 20 : (i + 1) * 20] *= 3

    cells = pd.DataFrame(
        {"label": np.repeat(["type_a", "type_b", "type_c"], 30)},
        index=[f"cell{i}" for i in range(90)],
    )
    coordinates = np.array([(x, y) for x in range(9) for y in range(8)])
    reference = ad.AnnData(
        csr_matrix(rng.poisson(np.repeat(profiles, 30, axis=0)).astype("float32")),
        obs=cells,
        var=pd.DataFrame(index=genes),
    )
    spatial = ad.AnnData(
        csr_matrix(
            rng.poisson(profiles[coordinates[:, 0] // 3] * 4 + profiles.sum(axis=0))
            .astype("float32")
        ),
        obs=pd.DataFrame(index=[f"spot{i}" for i in range(72)]),
        var=pd.DataFrame(index=genes),
        obsm={"spatial": coordinates},
    )
    for name, adata in {"reference": reference, "spatial": spatial}.items():
        adata.write_h5ad(tmp_path / f"{name}.h5ad")
    return reference, spatial


def _call(tu, name, operation, output_dir, **kwargs):
    request = {"operation": operation, "output_dir": str(output_dir), "max_epochs": 3, **kwargs}
    result = tu.run_one_function({"name": name, "arguments": {"request": request}})
    assert result.get("operation") == operation, result
    Draft202012Validator(tu.all_tool_dict[name]["return_schema"]).validate(result)
    return result


def test_category_and_invalid_request(universe, tmp_path):
    """Expose exactly three named tools and reject invalid operation arguments."""
    assert set(universe.all_tool_dict) == TOOL_NAMES
    assert set(universe.tool_category_dicts) >= {CATEGORY}

    for name in TOOL_NAMES:
        schema = universe.all_tool_dict[name]
        Draft202012Validator.check_schema(schema["parameter"])
        Draft202012Validator.check_schema(schema["return_schema"])

    invalid = {
        "operation": "deconvolve",
        "spatial_path": "unused.h5ad",
        "output_dir": str(tmp_path / "invalid"),
    }
    result = universe.run_one_function({"name": "GraphST", "arguments": {"request": invalid}})
    assert "reference_path" in result["error"]
    assert "deconvolve" in result["error"]
    assert not Path(invalid["output_dir"]).exists()


def test_cell2location(universe, expression_data, tmp_path):
    """Feed the fitted reference CSV back into spatial abundance estimation."""
    pytest.importorskip("cell2location")
    import anndata as ad
    import numpy as np
    import pandas as pd

    result = _call(
        universe, "cell2location", "fit_reference", tmp_path / "results",
        reference_path=str(tmp_path / "reference.h5ad"), cell_type_key="label",
    )
    artifact = result["reference_signatures"]
    signatures = pd.read_csv(artifact["path"], index_col=0)
    assert list(signatures.index) == list(expression_data[0].var_names)
    assert list(signatures.columns) == ["type_a", "type_b", "type_c"]

    result = _call(
        universe, "cell2location", "map_abundance", tmp_path / "results",
        spatial_path=str(tmp_path / "spatial.h5ad"),
        reference_signatures_path=artifact["path"], n_cells_per_location=6,
    )
    result_data = ad.read_h5ad(result["spatial_abundance"]["path"])
    assert result_data.obs_names.equals(expression_data[1].obs_names)
    assert result["spatial_abundance"]["shape"] == [72, 3]
    for key in ("means_cell_abundance_w_sf", "q05_cell_abundance_w_sf"):
        assert result_data.obsm[key].shape == (72, 3)
        assert list(result_data.obsm[key].columns) == list(signatures.columns)
        assert np.isfinite(result_data.obsm[key].to_numpy()).all()


@pytest.mark.parametrize("mode,n_rows", [("cells", 90), ("clusters", 3)])
def test_tangram(universe, expression_data, tmp_path, mode, n_rows):
    """Map both cells and groups and read back the optional projections."""
    pytest.importorskip("tangram")
    import anndata as ad
    import numpy as np

    result = _call(
        universe, "Tangram", "map_to_space", tmp_path / "results",
        reference_path=str(tmp_path / "reference.h5ad"),
        spatial_path=str(tmp_path / "spatial.h5ad"),
        mode=mode, cell_type_key="label", project_cell_types=True, project_genes=True,
    )
    mapping = ad.read_h5ad(result["mapping"]["path"])
    assert mapping.shape == (n_rows, 72)
    assert mapping.var_names.equals(expression_data[1].obs_names)
    np.testing.assert_allclose(mapping.X.sum(axis=1), 1, rtol=1e-6)
    types = ad.read_h5ad(result["cell_type_projection"]["path"])
    assert types.obsm["tangram_ct_pred"].shape == (72, 3)
    assert ad.read_h5ad(result["gene_projection"]["path"]).shape == (72, 90)


def test_graphst(universe, expression_data, tmp_path):
    """Check domains, integrated slice IDs, and deconvolution column ordering."""
    pytest.importorskip("GraphST")
    pytest.importorskip("igraph")
    import anndata as ad
    import numpy as np

    reference, spatial = expression_data
    spatial_path = str(tmp_path / "spatial.h5ad")
    output = tmp_path / "results"
    result = _call(universe, "GraphST", "cluster_domains", output, spatial_path=spatial_path)
    domains = ad.read_h5ad(result["domains"]["path"])
    assert domains.obs_names.equals(spatial.obs_names)
    assert domains.obsm["emb"].shape[0] == 72
    assert np.isfinite(domains.obsm["emb"]).all()
    assert domains.obs["domain"].notna().all()

    slices = []
    for slice_id, section in {"a": spatial[:36], "b": spatial[36:]}.items():
        section = section.copy()
        section.obs_names = spatial.obs_names[:36]  # Both slices reuse spot IDs.
        path = tmp_path / f"{slice_id}.h5ad"
        section.write_h5ad(path)
        slices.append({"slice_id": slice_id, "spatial_path": str(path)})
    result = _call(universe, "GraphST", "integrate_slices", output, slices=slices)
    integrated = ad.read_h5ad(result["domains"]["path"])
    assert integrated.n_obs == 72
    assert integrated.obs_names.is_unique
    assert set(integrated.obs["slice_id"]) == {"a", "b"}
    assert integrated.obs["original_spot_id"].duplicated().sum() == 36
    assert integrated.obs["domain"].notna().all()

    result = _call(
        universe, "GraphST", "deconvolve", output, spatial_path=spatial_path,
        reference_path=str(tmp_path / "reference.h5ad"), cell_type_key="label",
    )
    deconvolved = ad.read_h5ad(result["mapping"]["path"])
    assert deconvolved.obsm["map_matrix"].shape == (72, 90)
    assert list(deconvolved.uns["map_reference_cell_ids"]) == list(reference.obs_names)
    assert deconvolved.obsm["cell_type_proportions"].shape == (72, 3)
    np.testing.assert_allclose(deconvolved.obsm["cell_type_proportions"].sum(axis=1), 1)


def test_mcp_schemas(universe):
    """Preserve the full operation schemas on the MCP discovery surface."""
    from fastmcp import Client

    from tooluniverse.smcp import SMCP

    async def check():
        server = SMCP(
            tooluniverse_config=universe, tool_categories=[CATEGORY],
            search_enabled=False, auto_expose_tools=True,
        )
        try:
            async with Client(server) as client:
                tools = {tool.name: tool for tool in await client.list_tools()}
                assert set(tools) >= TOOL_NAMES
                for name in TOOL_NAMES:
                    assert tools[name].inputSchema == universe.all_tool_dict[name]["parameter"]
        finally:
            await server.close()

    asyncio.run(check())
