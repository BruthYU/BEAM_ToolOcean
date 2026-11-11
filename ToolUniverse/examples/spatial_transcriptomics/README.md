# Spatial transcriptomics

The `spatial_transcriptomics` category contains three local tools named `cell2location`, `Tangram`, and `GraphST`. Each tool accepts a `request` object whose `operation` selects the workflow.

## Installation

Use the Python environment that runs ToolUniverse. Install the dependencies for the tools you intend to call. The three Git revisions below are the upstream versions used in the integration tests.

```bash
python -m pip install -e ./ToolUniverse
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

Paths refer to the filesystem of the process running ToolUniverse, including when it serves MCP requests. Every call requires `output_dir` and writes its results into a new subdirectory. Each matrix descriptor returns the absolute path, shape, and axis names. The tool configuration defines the complete input and result schemas.

## Operations

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

All six operations accept `max_epochs`. The JSON schema lists the default for each workflow.
