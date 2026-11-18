# SpatialFormer and imaging single-cell data

The `spatial_transcriptomics` category adds four tools: `SpatialSingleCellIO`, `SpatialFormer`, `SpatialPairValidation`, and `SpatialPairSummary`. All operations write to a new directory under `output_dir`; reuse returned paths rather than predicting generated directory names.

`SpatialSingleCellIO` standardizes **already segmented RNA counts**. It does not segment images, assign unassigned molecules, infer cell types, convert targeted assays into whole-transcriptome scRNA-seq, or fill unmeasured genes with zeros. Importing CosMx protein measurements as RNA is unsupported. Labels are optional for SpatialFormer inference but mandatory for reference fitting and cell-type summaries.

## Shared call helper

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse(load_workspace=False)
tu.load_tools(categories=["spatial_transcriptomics"], clear_existing=True)

def call(name, operation, **parameters):
    result = tu.run_one_function({
        "name": name,
        "arguments": {"request": {"operation": operation, **parameters}},
    })
    if "error" in result:
        raise RuntimeError(result)
    return result

# Call tu.close() when the workflow finishes, preferably in a finally block.
```

## Import segmented RNA counts

Common outputs are `output_path` (AnnData) and `report_path` (JSON). Counts are stored in `X` and `layers["counts"]`; `var["gene_name"]` supplies model gene symbols. Standardized cell IDs encode the sample and original ID without delimiter ambiguity, and `obs["source_cell_id"]` preserves the original ID. Metadata includes `sample_id`, optional `fov_id`, and provenance in `uns["spatial_single_cell"]`.

Coordinates, when present, require explicit `coordinate_unit` (`micrometer` or `pixel`) and `coordinate_frame` (`global` or `local`). Local coordinates require `fov_key`; independent fields of view must not be treated as one continuous slide. Import does not register fields of view or convert pixels to micrometers.

### Xenium directory

The native directory adapter reads `cell_feature_matrix.h5` and, if available, `cells.csv.gz`, `cells.csv`, or `cells.parquet`. It recognizes `x_centroid`/`y_centroid`. Feature types other than `Gene Expression` are excluded and counted in the report.

```python
xenium = call(
    "SpatialSingleCellIO", "import_cells",
    platform="Xenium", input_path="/data/xenium/outs",
    sample_id="xenium_01", coordinate_unit="micrometer",
    coordinate_frame="global", output_dir="/results/imports",
)
```

### MERFISH matrix and metadata tables

MERFISH exports vary across processing pipelines. Use explicit identifier and coordinate columns. This example assumes a cell-by-gene CSV with a `cell_id` column, and a separate metadata CSV with matching IDs and global coordinates. Replace the example gene list with the complete intended RNA panel; do not include metadata/control columns.

```python
merfish = call(
    "SpatialSingleCellIO", "import_cells",
    platform="MERFISH", format="csv",
    input_path="/data/merfish/cell_by_gene.csv",
    metadata_path="/data/merfish/cell_metadata.csv",
    sample_id="merfish_01", cell_id_columns=["cell_id"],
    metadata_id_columns=["cell_id"], gene_columns=["CD3D", "MS4A1", "COL1A1"],
    coordinate_columns=["center_x", "center_y"],
    coordinate_unit="micrometer", coordinate_frame="global",
    output_dir="/results/imports",
)
```

### CosMx matrix with composite FOV/cell IDs

If cell IDs repeat between fields of view, use both FOV and cell ID as a composite key, in the same order in expression and metadata tables. This example intentionally uses local pixel coordinates. Such coordinates support within-FOV validation in pixels only; use calibrated global coordinates if you need a physical radius in micrometers.

```python
cosmx = call(
    "SpatialSingleCellIO", "import_cells",
    platform="CosMx", format="csv",
    input_path="/data/cosmx/expression.csv",
    metadata_path="/data/cosmx/metadata.csv",
    sample_id="cosmx_01", cell_id_columns=["fov", "cell_ID"],
    metadata_id_columns=["fov", "cell_ID"],
    gene_columns=["CD3D", "MS4A1", "COL1A1"],
    coordinate_columns=["CenterX_local_px", "CenterY_local_px"],
    coordinate_unit="pixel", coordinate_frame="local", fov_key="fov",
    output_dir="/results/imports",
)
```

The tabular reader accepts CSV, TSV, and Parquet; `format="csv"` selects this tabular adapter. Non-gene columns in an expression table are retained as cell metadata. Missing/duplicate cell IDs, unmatched metadata rows, noninteger counts, and nonfinite coordinates are rejected rather than silently repaired.

### Existing AnnData

```python
converted = call(
    "SpatialSingleCellIO", "import_cells",
    platform="Xenium", format="h5ad", input_path="/data/segmented.h5ad",
    sample_id="sample_01", count_layer="counts", gene_name_key="gene_symbol",
    coordinate_unit="micrometer", coordinate_frame="global",
    output_dir="/results/imports",
)
```

The selected layer must contain original integer RNA counts. Do not select log-normalized expression. If coordinates are absent, coordinate arguments can be omitted; embedding inference still works, but coordinate validation does not.

## Attach externally verified annotations

Labels are matched by exact ID, never by row position. Use standardized `obs_names` by default, or `source_cell_id` for source identifiers. For composite source identifiers, provide the original identifier columns in the annotation CSV.

```python
annotated = call(
    "SpatialSingleCellIO", "attach_annotations",
    input_path=cosmx["output_path"],
    annotation_path="/data/cosmx/cell_types.csv",
    annotation_id_columns=["fov", "cell_ID"],
    input_id_key="source_cell_id", label_column="author_cell_type",
    output_key="cell_type", require_all=True, output_dir="/results/annotations",
)
```

Outputs include an `audit_path` describing each match. `require_all=False` records unmatched/blank labels as `unassigned`; these labels must be resolved or explicitly filtered before reference preparation. Existing output label columns are never overwritten.

## Prepare a cell2location/Tangram reference

```python
reference = call(
    "SpatialSingleCellIO", "prepare_reference",
    input_path=annotated["output_path"], cell_type_key="cell_type",
    count_layer="counts", batch_key="sample_id", min_cells_per_type=5,
    spatial_path="/data/visium.h5ad", subset_shared_genes=False,
    output_dir="/results/references",
)

signatures = call(
    "cell2location", "fit_reference",
    reference_path=reference["reference_path"], cell_type_key="cell_type",
    count_layer="counts", batch_key="sample_id", output_dir="/results/cell2location",
)
```

The reference report lists cells per type and, when `spatial_path` is supplied, exact gene-ID overlap. Use `SpatialPreprocess.align_gene_ids` first if one dataset uses symbols and the other uses Ensembl IDs. Reference conversion does not synthesize missing genes. Keep `subset_shared_genes=False` when retaining reference-only genes for Tangram projection. A format-compatible imaging reference still requires validation for panel coverage, segmentation contamination, tissue/type coverage, and assay differences.

## SpatialFormer inference only

Prepare single-input and paired-input model bundles **during deployment** on the network-connected host. Make them available to the offline GPU worker. A bundle contains pinned checkpoint/configuration/tokenizer resources and its manifest. The deployment script and validated bundle paths are documented with the deployment artifacts; downloading weights is not a tool operation.

Only `extract_embeddings` and `predict_pairs` are exposed. No training or fine-tuning is invoked. `model_bundle_path` must match the operation's single/pair mode. The wrapper validates metadata against the actual model vocabulary. Import uses platform name `MERFISH`, while upstream model token spelling is **`Merfish`**.

```python
embedded = call(
    "SpatialFormer", "extract_embeddings",
    input_path=xenium["output_path"], count_layer="counts",
    model_bundle_path="/models/spatialformer/single",
    assay="Xenium", tissue="Lung", condition="Disease",
    device="cuda", batch_size=4, num_workers=0,
    output_dir="/results/spatialformer",
)
```

Output `output_path` contains the original data plus `obsm["X_SpaF"]`. The vector dimension is read from the model output, not assumed to be 512. `summary_path` records gene coverage, device, model provenance, and platform warnings. MERFISH/CosMx inference is a cross-platform application of these pretrained weights; accepted tokens alone do not establish accuracy. Ordinary dissociated scRNA-seq is not given an invented assay token.

## Build a coordinate benchmark and predict explicit pairs

```python
benchmark = call(
    "SpatialPairValidation", "build_pair_benchmark",
    input_path=xenium["output_path"], sample_key="sample_id",
    frame_key="sample_id", coordinate_unit="micrometer", radius=30,
    pairs_per_class=100, random_seed=0, output_dir="/results/pairs",
)

predicted = call(
    "SpatialFormer", "predict_pairs",
    input_path=xenium["output_path"], count_layer="counts",
    model_bundle_path="/models/spatialformer/pair",
    assay="Xenium", tissue="Lung", condition="Disease",
    pairs_path=benchmark["pairs_path"], sample_key="sample_id",
    device="cuda", batch_size=4, num_workers=0, max_pairs=10000,
    output_dir="/results/spatialformer",
)

validation = call(
    "SpatialPairValidation", "evaluate_pairs",
    prediction_path=predicted["predictions_path"],
    benchmark_path=benchmark["benchmark_path"], pretraining_overlap="unknown",
    output_dir="/results/validation",
)
```

Use `frame_key="fov_id"` for independent local fields of view. Set radius in the declared input coordinate units. `frame_key="sample_id"` is appropriate only when coordinates in each sample form a shared registered frame. The benchmark includes spatially near and far pairs, rather than only known neighbors. It writes a pair list for inference and a separate table with distances/labels for evaluation.

Custom pair CSVs require `sample_id`, `cell_id_a`, and `cell_id_b`, matching standardized AnnData observation IDs exactly. Cross-sample pairs, self-pairs, duplicate ordered pairs, and unrecognized IDs are rejected. Pair inference never creates an unbounded all-by-all matrix. `max_pairs` and `max_dense_bytes` constrain the upstream implementation's memory use.

Predictions include `proximity_score` (positive-class output), model ID, and status. A score is not a measured distance, independently calibrated contact probability, or evidence of cell communication. Pair ordering is preserved; reverse-pair averaging is not silently performed. Evaluation joins the same ordered pairs and reports AUROC/AUPRC and group metrics. An optional fixed `threshold` enables threshold metrics; it is not tuned on the evaluated results. Report `pretraining_overlap` honestly: unknown overlap does not establish a fully independent test. Precision on a balanced sampled benchmark does not estimate prevalence in intact tissue.

## Summarize cell-type pairs and compare patient groups

Use the same underlying cells with externally attached type labels. Predictions alone do not create these labels.

```python
summary = call(
    "SpatialPairSummary", "summarize_type_pairs",
    prediction_path=predicted["predictions_path"],
    input_path="/results/xenium_with_cell_types.h5ad",
    cell_type_key="cell_type", sample_key="sample_id", frame_key="sample_id",
    benchmark_path=benchmark["benchmark_path"],
    sampling_description="Seeded near/far benchmark; equal requested counts per class",
    output_dir="/results/type_pairs",
)
```

The output `summary_path` records denominators and average scores per type pair. If a benchmark is supplied it also records observed near-pair fractions. These describe the sampled pairs; biased sampling is not automatically reweighted into population contact frequencies.

For multiple patients, concatenate compatible summary tables from the same model and sampling protocol. Supply metadata CSV columns `sample_id,patient_id,group`; each sample has one metadata row and each patient belongs to one group.

```python
comparison = call(
    "SpatialPairSummary", "compare_groups",
    summary_path="/results/cohort_type_pairs.csv",
    metadata_path="/data/cohort_metadata.csv",
    compare_groups=["treated", "control"], metric="mean_score",
    n_permutations=9999, random_seed=0, fdr_method="bh",
    output_dir="/results/group_comparison",
)
```

Samples contribute equally within a patient, and inference uses patients rather than individual cell pairs as replicates. Insufficient replication is reported rather than turned into a significant result. Group differences and FDR apply to the declared analysis, not previous searches or model/threshold selection.

`SpatialPairSummary` complements `SpatialAssociation`: the former analyzes predicted cell-pair proximity; the latter analyzes spatial abundance profiles from cell2location/Tangram/GraphST. Compare biological trends at matched sample/type/region levels, not numerical equality of these different scores. Cell-level and spot-level data need an explicit spatial correspondence before any direct comparison.

## Worker smoke checks

Run scientific tests on the cluster worker, with dependencies and weights prepared in advance. `tests/integration/test_spatial_single_cell_io.py` exercises three platform imports, exact annotation alignment, reference export, native Xenium control filtering, and rejection of ambiguous IDs. The full smoke workflow should additionally run both real pretrained inference modes and the coordinate-validation/type-summary chain. Do not describe cached predictions or mock inference as a fresh model run.

## Provisioning official weights (deployment only)

The tools never download or train models. On the networked development host, run:

```bash
python ToolUniverse/examples/spatialformer/prepare_runtime.py --output-dir /shared/spatialformer-runtime
python -m pip download --no-deps --dest /shared/spatialformer-wheels communities==3.0.0 transformers==4.57.3 huggingface-hub==0.36.0
```

The preparation script pins the official Git commit, verifies Git blob hashes and published checkpoint MD5 values, and records SHA-256 hashes in each bundle manifest. Bundles are `/shared/spatialformer-runtime/single` and `/shared/spatialformer-runtime/pair`. They include references to the adjacent `source` resources: keep the entire runtime directory together. By default, both inference modes share the public pair-trained checkpoint, without duplicating its storage. The single-input manifest retains `training_mode="pair"` and the same model ID; it does not represent a separately single-trained model. Pass `--include-single-checkpoint` to optionally prepare the distinct single-trained checkpoint. Interrupted checkpoint downloads retain a `.partial` file and resume before final checksum verification. The source is MIT-licensed; the published checkpoint metadata specifies CC BY 4.0. License metadata and attribution URLs are retained in the manifests.

For the validated slime base image, install offline on the worker:

```bash
bash ToolUniverse/examples/spatialformer/install_offline.sh /shared/spatialformer-runtime /shared/spatialformer-wheels
/opt/spatialformer-venv/bin/python ToolUniverse/examples/spatialformer/smoke_pipeline.py --runtime /shared/spatialformer-runtime --output /shared/new-smoke-run
```

The smoke deliberately uses small synthetic counts with real model weights, calls the registered SDK tools serially, and also passes the prepared reference to cell2location. It is an execution test, not a biological performance claim. The output directory must be new, avoiding accidental reuse of cached results.

The inference adapter also handles two upstream compatibility cases without editing the checkpoint or source package: single-cell batches whose first cell has exactly two expressed genes, and trusted, checksum-verified Lightning checkpoints on recent PyTorch. Tissue context is restricted to actual tissue tokens rather than arbitrary vocabulary entries.

The pinned upstream tokenizer always adds the `Human` species token. This integration therefore supports human samples only; platform compatibility does not imply support for mouse or other species.

## Validated cluster deployment — 2026-11-18

Project root: `/mnt/shared-storage-gpfs2/gpfs2-agent-math/huangwenxuan/BEAM_ToolOcean`.

Published image: `registry.h.pjlab.org.cn/ailab-deepdivegzy-deepdivegzy_gpu/slime:v0`, verified digest `sha256:2747b6f7517de852ad6121438179567c63ec9f68ed5beba1a086d3c4eaf29242`. The previous 117 layers are unchanged; one 22 MB layer adds the SpatialFormer environment and four tool definitions. Use `/opt/spatialformer-venv/bin/python` for workflows containing these tools. The original `/opt/toolocean-venv` remains available.

The shared runtime is `<project root>/cluster-runtime/spatialformer/runtime`. Set `model_bundle_path` to its `single` directory for `extract_embeddings` or `pair` directory for `predict_pairs`. Both refer to the same official pair-trained checkpoint, whose SHA-256 is `e01526ba2683051503070f0f1eedbd98ba024992986569423b1db5ee4d7b0b16`. Weights remain on shared storage; they are not embedded in the image. Mount `gpfs://gpfs2/gpfs2-agent-math` at `/mnt/shared-storage-gpfs2/gpfs2-agent-math` on workers.

Validation completed on NVIDIA H200 with CUDA enabled:

- 14 data, post-analysis, model-adapter and category tests passed; the separate MCP schema-discovery test also passed.
- A fresh candidate-image run used installed `/opt` code with `PYTHONPATH` unset and passed all 16 serial SDK calls: four imports and annotation joins, reference export, cell2location reference fitting, both SpatialFormer inference operations, coordinate benchmark/evaluation, type-pair summary and patient-group comparison.
- The 16-cell synthetic input produced finite `16 × 512` embeddings and predictions for 16 supplied pairs. Only cell2location's existing reference stage used a one-epoch fit; SpatialFormer performed inference only.
- GPU tasks exited successfully and their resources were released. This small execution check does not estimate biological accuracy.

Evidence under `<project root>/cluster-runtime/spatialformer/`: `jobs/12-tests.sh.log`, `jobs/13-mcp.sh.log`, `candidate-validation.log`, `candidate-smoke-01/events.json`, `candidate-smoke-01/SUCCESS.json`, and `image-validation.json`. The candidate tag `spatialformer-candidate-20261118` retains the same published digest.
