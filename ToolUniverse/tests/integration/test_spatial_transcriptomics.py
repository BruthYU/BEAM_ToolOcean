"""Exercise the local spatial transcriptomics tools through the SDK and MCP."""

import asyncio
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from tooluniverse import ToolUniverse

CATEGORY = "spatial_transcriptomics"
TOOL_NAMES = {
    "cell2location", "Tangram", "GraphST", "SpatialIO", "SpatialPreprocess", "SpatialAssociation",
}
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
    """Expose exactly the expected tools and reject invalid operation arguments."""
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


def _data_call(tu, name, operation, **request):
    args = {"name": name, "arguments": {"request": {"operation": operation, **request}}}
    result = tu.run_one_function(args)
    if "error" not in result:
        Draft202012Validator(tu.all_tool_dict[name]["return_schema"]).validate(result)
    return result


def test_import_visium_matrix_market(universe, tmp_path):
    """Read a raw matrix and match a headered position file by barcode."""
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    from scipy.io import mmwrite
    from scipy.sparse import coo_matrix

    source = tmp_path / "sample" / "outs"
    matrix = source / "raw_feature_bc_matrix"
    spatial = source / "spatial"
    matrix.mkdir(parents=True)
    spatial.mkdir()
    mmwrite(matrix / "matrix.mtx", coo_matrix([[1, 0], [3, 5]]))
    (matrix / "genes.tsv").write_text("ENSG0001\tA\nENSG0002\tA\n")
    (matrix / "barcodes.tsv").write_text("spot-a\nspot-b\n")
    (spatial / "tissue_positions.csv").write_text(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
        "spot-b,0,1,2,100,200\nspot-a,1,0,1,70,90\n"
    )
    result = _data_call(
        universe, "SpatialIO", "import_visium", input_dir=str(source),
        output_path=str(tmp_path / "imported.h5ad"), matrix_kind="raw",
        image_levels=[], tissue_only=False,
    )
    imported = ad.read_h5ad(result["output_path"])
    assert imported.obs_names.tolist() == ["spot-a", "spot-b"]
    assert imported.var_names.tolist() == ["ENSG0001", "ENSG0002"]
    assert imported.obsm["spatial"].tolist() == [[90, 70], [200, 100]]
    assert result["library_id"] == "sample"
    assert result["n_outside_tissue_retained"] == 1


def test_spatial_qc_normalization_and_export(universe, tmp_path):
    """Preserve counts and coordinates while filtering, normalizing, and exporting."""
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    spots = pd.DataFrame({"slice_id": ["a", "a", "b"]}, index=["s1", "s2", "s3"])
    genes = pd.DataFrame({"gene_symbol": ["MT-A", "B"]}, index=["gene1", "gene2"])
    adata = ad.AnnData(
        X=csr_matrix([[1, 3], [0, 0], [0, 9]]),
        obs=spots, var=genes,
        obsm={"spatial": np.array([[0, 0], [1, 1], [2, 2]])},
    )
    source = str(tmp_path / "source.h5ad")
    adata.write_h5ad(source)
    inspection = _data_call(universe, "SpatialIO", "inspect_spatial_h5ad", input_path=source)
    assert inspection["issues"] == []
    assert inspection["count_check"]["rows_checked"] == 3

    preview = _data_call(
        universe, "SpatialPreprocess", "filter_spatial_qc",
        input_path=source, output_path=str(tmp_path / "preview.csv"),
        max_pct_mito=20, group_by="slice_id",
    )
    assert preview["n_pass"] == 1
    assert preview["failed_by_rule"] == {"max_pct_mito": 2}
    assert [group["n_pass"] for group in preview["groups"]] == [0, 1]
    filtered = _data_call(
        universe, "SpatialPreprocess", "filter_spatial_qc", mode="filter",
        input_path=source, output_path=str(tmp_path / "filtered.h5ad"), min_counts=1,
    )
    assert ad.read_h5ad(filtered["output_path"]).obs_names.tolist() == ["s1", "s3"]

    normalized = _data_call(
        universe, "SpatialPreprocess", "normalize_spatial_expression",
        input_path=source, output_path=str(tmp_path / "norm.h5ad"),
        log1p=False, target_sum=100, size_factor_key="size_factor",
    )
    data = ad.read_h5ad(normalized["output_path"])
    np.testing.assert_array_equal(data.X.toarray(), adata.X.toarray())
    np.testing.assert_allclose(data.layers["normalized"].sum(axis=1).A1, [100, 0, 100])
    np.testing.assert_allclose(data.obs["size_factor"], [0.04, 0, 0.09])
    assert normalized["n_zero_spots"] == 1

    request = {
        "input_path": normalized["output_path"], "output_path": str(tmp_path / "spots.tsv"),
        "spot_ids": ["s3", "s1"], "format": "tsv", "obs_columns": ["slice_id"],
        "genes": ["gene2"], "expression_layer": "normalized",
    }
    exported = _data_call(universe, "SpatialIO", "export_spatial_table", **request)
    table = pd.read_csv(exported["output_path"], sep="\t")
    assert table["spot_id"].tolist() == ["s3", "s1"]
    assert exported["column_sources"]["gene:gene2"] == "layers/normalized/gene2"
    np.testing.assert_allclose(table["gene:gene2"], [100, 75])
    rejected = _data_call(universe, "SpatialIO", "export_spatial_table", **request)
    assert "Output already exists" in rejected["error"]


def test_align_gene_ids_with_mapping(universe, tmp_path):
    """Keep reference-only genes and reject ambiguous mappings."""
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd

    spatial = ad.AnnData(X=np.array([[1, 2]]), var=pd.DataFrame(index=["A", "B"]))
    reference = ad.AnnData(
        X=np.array([[3, 4, 5]]),
        var=pd.DataFrame(index=["ENSG0001.1", "ENSG0002.2", "ENSG0003.1"]),
    )
    paths = {side: str(tmp_path / f"{side}.h5ad") for side in ("spatial", "reference")}
    spatial.write_h5ad(paths["spatial"])
    reference.write_h5ad(paths["reference"])
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("source_id,target_id\nA,ENSG0002\nB,ENSG0001\n")
    request = {
        "spatial_path": paths["spatial"], "reference_path": paths["reference"],
        "output_dir": str(tmp_path), "mapping_path": str(mapping), "mapping_side": "spatial",
        "strip_ensembl_version": True,
    }
    preserved = _data_call(universe, "SpatialPreprocess", "align_gene_ids", **request)
    assert preserved["n_shared"] == 2
    assert preserved["n_reference_only"] == 1
    assert ad.read_h5ad(preserved["reference_path"]).var_names.tolist() == [
        "ENSG0001", "ENSG0002", "ENSG0003",
    ]
    aligned = _data_call(universe, "SpatialPreprocess", "align_gene_ids", mode="intersection", **request)
    result = ad.read_h5ad(aligned["reference_path"])
    assert result.var_names.tolist() == ["ENSG0002", "ENSG0001"]
    np.testing.assert_array_equal(result.X, [[4, 3]])
    mapping.write_text("source_id,target_id\nA,ENSG0001\nB,ENSG0001\n")
    rejected = _data_call(universe, "SpatialPreprocess", "align_gene_ids", **request)
    assert "Duplicate target gene IDs" in rejected["error"]


def test_spatial_association_uses_patient_as_unit(universe, tmp_path):
    """Weight slices within patients and permute patient labels across groups."""
    pytest.importorskip("matplotlib")
    import pandas as pd

    rows = []
    sections = [
        ("P01", "response", "large", 30, 1),
        ("P01", "response", "small", 5, -1),
        ("P02", "response", "large", 7, 1),
        ("P03", "no_response", "large", 5, -1),
        ("P04", "no_response", "large", 5, -1),
    ]
    for patient, group, slice_id, n, direction in sections:
        for i in range(n):
            rows.append({
                "patient": patient, "response": group, "slice_id": slice_id,
                "domain": "d", "x": i, "y": 0,
                "T": i + 1, "B": i + 1 if direction == 1 else n - i,
            })
    path = tmp_path / "abundance.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    request = {
        "input_path": str(path), "output_dir": str(tmp_path), "cell_types": ["T", "B"],
        "patient_key": "patient", "group_key": "response", "domain_key": "domain",
        "method": "colocalization", "compare_groups": ["response", "no_response"],
        "n_permutations": 99,
    }
    result = _data_call(universe, "SpatialAssociation", "analyze_spatial_association", **request)
    comparison = result["group_comparison"]
    assert comparison["unit"] == "patient"
    assert comparison["n_patients_by_group"] == {"response": 2, "no_response": 2}
    assert comparison["permutation_method"] == "exact"
    assert comparison["n_permutations"] == 6
    assert comparison["difference"] == pytest.approx(1.5)
    assert comparison["p_value"] == pytest.approx(1 / 3)
    patients = pd.read_csv(result["patient_summary_path"]).set_index("patient_id")
    assert patients.loc["P01", "n_spots_used"] == 35
    assert patients.loc["P01", "score"] == pytest.approx(0)
    assert patients["spatial_p_value"].notna().all()
    assert Path(result["figure_path"]).is_file()

    rows[0]["response"] = "no_response"
    pd.DataFrame(rows).to_csv(path, index=False)
    rejected = _data_call(universe, "SpatialAssociation", "analyze_spatial_association", **request)
    assert "Each patient must belong to exactly one group" in rejected["error"]


def test_spatial_screen_global_family_and_slice_weighting(universe, tmp_path):
    """Retain untestable pairs in global FDR and weight patients independently of spot counts."""
    import pandas as pd
    import numpy as np
    import anndata as ad

    rows = []
    for patient, group, section, n, direction in [
        ('p1', 'R', 's1', 30, 1), ('p1', 'R', 's2', 5, -1),
        ('p2', 'R', 's1', 7, 1), ('p3', 'NR', 's1', 5, -1), ('p4', 'NR', 's1', 5, -1),
    ]:
        for i in range(n):
            rows.append(dict(patient=patient, group=group, slice_id=section, domain='d',
                             x=i, y=0, A=i + 1, B=i + 1 if direction == 1 else n - i, C=1))
    frame = pd.DataFrame(rows)
    csv = tmp_path / 'screen.csv'
    frame.to_csv(csv, index=False)
    h5ad = tmp_path / 'screen.h5ad'
    obs = frame.drop(columns=['x', 'y']).copy()
    obs.index = obs.index.astype(str)
    ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs, obsm={'spatial': frame[['x', 'y']].to_numpy()}).write_h5ad(h5ad)
    request = dict(datasets=[dict(dataset_id='csv', input_path=str(csv)), dict(dataset_id='h5ad', input_path=str(h5ad))],
                   output_dir=str(tmp_path), cell_types=['A', 'B', 'C'], patient_key='patient',
                   group_key='group', domain_key='domain', compare_groups=['R', 'NR'], top_n=1)
    result = _data_call(universe, 'SpatialAssociation', 'screen_spatial_associations', **request)
    assert result['number_tests'] == 6
    assert result['number_testable'] == 2
    assert len(result['top_pairs']) == 1
    table = pd.read_csv(result['screen_table_path'])
    valid = table[table.status == 'ok']
    assert valid.difference.tolist() == pytest.approx([1.5, 1.5])
    assert valid.p_value.tolist() == pytest.approx([1/3, 1/3])
    assert valid.q_value.tolist() == pytest.approx([1, 1])
    assert table[table.status != 'ok'].q_value.isna().all()
    patients = pd.read_csv(result['patient_summary_path'])
    assert patients.query("patient_id == 'p1' and cell_a == 'A' and cell_b == 'B'").score.tolist() == pytest.approx([0, 0])
    request['datasets'][1]['dataset_id'] = 'csv'
    rejected = _data_call(universe, 'SpatialAssociation', 'screen_spatial_associations', **request)
    assert 'dataset_id values must be unique' in rejected['error']
