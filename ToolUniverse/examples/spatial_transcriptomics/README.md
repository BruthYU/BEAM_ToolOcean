# Spatial transcriptomics

The `spatial_transcriptomics` category contains `SpatialIO`, `SpatialPreprocess`, `SpatialAssociation`, `cell2location`, `Tangram`, and `GraphST`. Each tool accepts a `request` object whose `operation` selects the workflow.

## Installation

Use the Python environment that runs ToolUniverse. Install the dependencies for the tools you intend to call. The three Git revisions below are the upstream versions used in the integration tests.

```bash
python -m pip install -e ./ToolUniverse
```

`SpatialIO` and `SpatialPreprocess` use AnnData, NumPy, pandas, and SciPy. Visium import and normalization also use Scanpy. Importing tissue images uses Matplotlib. Install these dependencies in the same environment.

```bash
python -m pip install anndata 'scanpy>=1.12,<1.13'
```

```bash
# cell2location
python -m pip install 'pyro-ppl==1.9.1' 'scvi-tools==1.5.2' \
  'cell2location @ git+https://github.com/BayraktarLab/cell2location.git@afb07fda7f89458e44a4afdf21c26aeb7251ebfa'

# Tangram
python -m pip install \
  'tangram-sc @ git+https://github.com/broadinstitute/Tangram.git@4c68995a418f41dc8caef567598c4d9b47781a13'

# GraphST
python -m pip install 'scanpy[skmisc]>=1.12,<1.13' torch POT igraph \
  'GraphST @ git+https://github.com/JinmiaoChenLab/GraphST.git@d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74'
```

The local wrappers run on CPU. GraphST clustering uses Scanpy's igraph Leiden implementation. The reference and mapping stages of cell2location use Pyro 1.9.1 because this upstream revision imports a helper removed in Pyro 1.9.2.

## Call a tool

Run this example from Python with ToolUniverse installed.

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse(load_workspace=False)
try:
    tu.load_tools(categories=["spatial_transcriptomics"], clear_existing=True)
    result = tu.run_one_function({
        "name": "GraphST",
        "arguments": {
            "request": {
                "operation": "cluster_domains",
                "spatial_path": "/data/section.h5ad",
                "resolution": 0.8,
                "output_dir": "/data/results",
            }
        },
    })
    print(result)
finally:
    tu.close()
```

Paths refer to the filesystem of the process running ToolUniverse, including when it serves MCP requests. Model operations require `output_dir` and write their results into a new subdirectory. Each model matrix descriptor returns the absolute path, shape, and axis names. The tool configuration defines the complete input and result schemas.

## Data I/O and preprocessing

Use `SpatialIO` to import or inspect a slice, then use `SpatialPreprocess` to prepare it. Files created by these operations can be passed directly to model operations. Operations with an explicit `output_path` require a new file path.

| Tool and operation | Required inputs | Options and outputs |
| --- | --- | --- |
| `SpatialIO.import_visium` | `input_dir`, `output_path` | Reads a classic Visium Space Ranger directory with a 10x HDF5 or Matrix Market count matrix. `matrix_kind` defaults to `filtered`, and `tissue_only` defaults to true. Set `matrix_kind="raw"` and `tissue_only=false` to retain off-tissue spots. `image_levels` defaults to `["lowres"]` and also accepts `"hires"` or an empty list. `library_id` defaults to the directory name, or the parent directory name for an `outs` directory. The output preserves original counts in `X`, Ensembl or source gene IDs in `var_names`, gene symbols in `var["gene_symbol"]`, and full-resolution pixel coordinates in x, y order in `obsm["spatial"]`. |
| `SpatialIO.inspect_spatial_h5ad` | `input_path` | Set `role="spatial"` to check coordinates or `role="reference"` and `cell_type_key` to check cell labels. `coordinate_key` defaults to `spatial`. `count_layer` defaults to `X`. `count_check` accepts `sample` (up to 128 evenly spaced rows), `full`, or `none`. The report lists issues and the scope of the count check. Count compatibility describes the values in the checked rows. |
| `SpatialIO.export_spatial_table` | `input_path`, `output_path` | Set `obs_columns`, `matrix_columns`, or `genes` to choose values. For example, set `matrix_columns=[{"key":"means_cell_abundance_w_sf","columns":["B_naive"]}]` to export one cell type. Columns may also be selected by zero-based index. `expression_layer` defaults to `X`, `format` defaults to `csv`, and `include_coordinates` defaults to true. `spot_ids` selects rows in the requested order. The table includes `spot_id` and the response reports the source of each column. |
| `SpatialPreprocess.filter_spatial_qc` | `input_path`, plus `output_path` in filter mode | `mode="report"` previews QC and `mode="filter"` writes a new h5ad. Thresholds are inclusive and optional. Set `min_counts`, `max_counts`, `min_detected_genes`, `max_detected_genes`, or `max_pct_mito`. Select counts with `count_layer`, mitochondrial genes with `mito_flag_key` or `gene_symbol_key`, and group summaries with `group_by`. Each mode returns an audit CSV with the rule failures per spot. |
| `SpatialPreprocess.align_gene_ids` | `spatial_path`, `reference_path`, `output_dir` | Select IDs with `spatial_id_key` and `reference_id_key`, both defaulting to `var_names`. Use `mapping_path` with `mapping_side` for an explicit CSV or TSV containing `source_id` and `target_id`. Set `strip_ensembl_version=true` to remove Ensembl suffixes. The default `preserve` mode keeps side-specific genes. `intersection` uses the same shared gene order in both outputs. |
| `SpatialPreprocess.normalize_spatial_expression` | `input_path`, `output_path` | Choose the source with `count_layer` (default `X`) and the scaling total with `target_sum` (default 10,000). `log1p` defaults to true. `output_layer` defaults to `lognorm` or to `normalized` with `log1p=false`. Optionally store original totals divided by the target sum in `obs` using `size_factor_key`. Zero-count spots stay zero. Original counts stay in place. |

QC report mode accepts an optional `output_path` for the audit CSV. When omitted, the response points to a CSV in a new local temporary directory. Filter mode requires at least one threshold and writes the audit CSV next to the new h5ad in a run-specific subdirectory. The default mitochondrial check uses the `MT-` prefix when gene symbols are available. Use a boolean flag column for other naming conventions. The percentage is undefined for zero-count spots, so those spots fail a requested mitochondrial percentage threshold.

Gene alignment rejects ambiguous mappings and duplicate target IDs. Unmapped genes keep their original IDs in `preserve` mode. The audit table records original and aligned IDs, mapping status, and shared status. Alignment changes the active gene axis and its layers. Existing `.raw` snapshots retain their own gene IDs. Keep `preserve` mode when Tangram needs to predict genes present only in the reference.

Run these CLI calls to import a local Visium slice and preview its quality metrics.

```bash
tu run SpatialIO '{"request":{"operation":"import_visium","input_dir":"/data/sample/outs","output_path":"/results/sample.h5ad"}}'
tu run SpatialIO '{"request":{"operation":"inspect_spatial_h5ad","input_path":"/results/sample.h5ad","count_check":"full"}}'
tu run SpatialPreprocess '{"request":{"operation":"filter_spatial_qc","input_path":"/results/sample.h5ad","min_counts":1000,"output_path":"/results/qc_preview.csv"}}'
```

## Model operations

| Tool | Operation | Input | Results |
| --- | --- | --- | --- |
| `cell2location` | `fit_reference` | `reference_path` and `cell_type_key` | A CSV with genes in rows and cell types in columns. Pass the returned `reference_signatures.path` to `map_abundance`. |
| `cell2location` | `map_abundance` | `spatial_path`, `reference_signatures_path`, and `n_cells_per_location` | A spatial AnnData with spot-by-cell-type estimates in `obsm["means_cell_abundance_w_sf"]` and `obsm["q05_cell_abundance_w_sf"]`. |
| `Tangram` | `map_to_space` | `reference_path` and `spatial_path` | An AnnData with cell-by-spot or group-by-spot mapping probabilities in `X`. Optional projections are written as separate AnnData files. |
| `GraphST` | `cluster_domains` | `spatial_path` | A spatial AnnData with the spot representation in `obsm["emb"]` and the domain label in `obs["domain"]`. |
| `GraphST` | `integrate_slices` | `slices`, a list of at least two objects with `slice_id` and `spatial_path` | A combined spatial AnnData with a joint representation, domain labels, `obs["slice_id"]`, and `obs["original_spot_id"]`. |
| `GraphST` | `deconvolve` | `spatial_path`, `reference_path`, and `cell_type_key` | A spatial AnnData with the spot-by-reference-cell matrix in `obsm["map_matrix"]`, its column IDs in `uns["map_reference_cell_ids"]`, and spot-by-cell-type fractions in `obsm["cell_type_proportions"]`. |

Use unique observation and gene IDs in each input. For cell2location, supply raw integer counts in `X` or select `count_layer`. Both operations accept `batch_key`. Set `detection_alpha` for the spatial assay and tissue. The wrapper defaults to 20.

For Tangram, `mode` defaults to `cells`. Set `cell_type_key` for `clusters` mode or when `project_cell_types` is true. Set `project_genes` to export predicted expression. Tangram matches gene IDs case-insensitively and intersects `training_genes`, when provided, with shared expressed genes.

For GraphST, provide raw counts in `X` and two spatial coordinates per spot in `obsm["spatial"]`. The slice integration operation expects coordinates already registered in the same coordinate system. The `resolution` argument for clustering and integration controls Leiden granularity. Cell-type projection retains the highest-weighted 10% of reference cells per spot (at least one) before normalizing the resulting cell-type weights.

All six model operations accept `max_epochs`. The JSON schema lists the default for each workflow.

## Spatial association and patient groups

Use `SpatialAssociation.analyze_spatial_association` to compare two estimated cell abundance profiles in tissue and summarize the comparison by patient group. A CSV must have one row per spot and named abundance, patient, group, and domain columns. Coordinate columns default to `x` and `y` and can be set with `x_key` and `y_key`. For h5ad input, metadata comes from `obs` and coordinates default to `obsm["spatial"]`. Read abundances from named `obs` columns or set `abundance_key` to a named `obsm` DataFrame such as `means_cell_abundance_w_sf`. Keep the abundance scale consistent across sections and groups.

```bash
tu run SpatialAssociation '{
  "request": {
    "operation": "analyze_spatial_association",
    "input_path": "/data/cohort.h5ad",
    "output_dir": "/results",
    "cell_types": ["B_naive", "T_CD4+"],
    "abundance_key": "means_cell_abundance_w_sf",
    "patient_key": "patient_id",
    "slice_key": "slice_id",
    "group_key": "response",
    "domain_key": "domain",
    "method": "neighborhood",
    "n_neighbors": 6,
    "compare_groups": ["responder", "nonresponder"],
    "n_permutations": 999,
    "random_seed": 0
  }
}'
```

- **Slice score.** `colocalization` uses Pearson correlation between the two abundances at the same spots. `neighborhood` uses a symmetric cross-abundance score based on centered, unit-length abundance vectors and a degree-normalized nearest-neighbor graph. Both scores range from −1 to 1. Positive values indicate that higher values of one profile tend to occur with higher values of the other at the chosen scale. The graph is built separately for every patient and slice. The optional `max_neighbor_distance` caps edges in the input coordinate units. Supply `slice_key` for data with multiple sections per patient. When omitted, the tool uses `slice_id` or `library_id` if available, then falls back to one slice per patient. Use comparable labels across slices when restricting the analysis with `domains`.
- **Within-patient spatial permutation.** The graph stays fixed and the second abundance is shuffled within each slice and domain. Slice scores are averaged equally within the patient. The spatial p-value is a two-sided Monte Carlo result under this shuffle model. It assumes spot values are exchangeable within each block and can be sensitive to spatial autocorrelation. Treat it as exploratory for an autocorrelated slice.
- **Between-group permutation.** Each patient supplies one score. The effect is the mean score in the first requested group minus the mean in the second. The two-sided test shuffles patient labels while preserving group sizes, enumerating all assignments when they fit the permutation budget. It requires at least two usable patients in each group and assumes patient labels are exchangeable under the null. The comparison is unadjusted and expects one group per patient.

Outputs include `slices.csv` with slice scores and any reason a score or permutation could not be computed, `patients.csv` with one row per patient and the spatial p-value, `group_comparison.json` with the patient-level comparison, and `patients.png` with one dot per scored patient. With one group or insufficient patient replication, the tool still writes descriptive results and returns a null group p-value.
