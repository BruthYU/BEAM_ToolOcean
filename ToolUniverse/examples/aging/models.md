# Aging model adapters

These inference adapters use published weights. The category does not train or
fine-tune models. Deploy assets with `download_models.py` on the networked host,
then run `prepare_model_bundles.py` on a compute worker. Every deployed bundle
contains a source URL, license identifier, model configuration and checksums.
Do not load bundles from an untrusted source: SenCID uses upstream joblib files.

## RNAAgeCalc

`predict_age` evaluates the original published elastic-net coefficients as
`intercept + sum(coefficient * log2(FPKM + 1))`. The included deployment script
exports the `all_tissues / GTExAge / all` predictor from the authors' RData, using
`rdata` only during deployment. This is an explicitly documented NumPy backend,
not a call to the R package.

The wrapper requires `species="human"`, `data_kind="bulk"`,
`tissue="all_tissues"`, `gene_id_type="SYMBOL"`, and either FPKM or integer counts.
Counts require a positive gene-length column (base pairs) in AnnData `var`;
conversion uses the total counts across the supplied complete gene universe.
Use `AgingPreprocess` for explicit gene-ID mapping. All signature genes must be
present. The adapter does not reproduce the R package's KNN missing-gene
imputation, automatic ID mapping, automatic gene-length database, or tissue
fallback. A different tissue/signature needs an explicitly exported matching
bundle. Pseudobulk input is not automatically validated as compatible bulk data.

The smoke test uses the author-provided FPKM example and an independent scalar
application of the same coefficient formula. This is algebraic numerical parity,
not an executed R-runtime comparison, and the report labels it accordingly.

## scAgeClock

`prepare_input` converts raw counts to the published feature layout. It maps
explicit assay/cell-type/tissue/sex labels using the bundled categorical tables,
or accepts explicitly declared integer category-index columns. Unknown categories
fail. Counts are normalized across all supplied genes to their median library
size by default and natural-log transformed, matching the upstream formatter's
default. Absent model genes are zero-filled and coverage is reported; insufficient
coverage fails. `target_sum` can explicitly fix the normalization reference.

`predict_age` requires the exact feature ordering and valid category-index ranges.
It uses the original GMA inference classes and checkpoint, loaded with
`torch.load(weights_only=True)`, in bounded batches on CPU or CUDA. The deployment
script extracts the three unchanged inference class bodies from upstream source
to avoid importing unrelated training algorithms. The original source is kept
and checksummed in the bundle. Normalization is performed before batching.

Official formatted h5ad input is supported without requiring true ages.
Preparation records its bundle hash; predictions reject mismatched prepared
bundles. Smoke compares full-batch CPU to batched CUDA predictions on the
32-cell subset of the official test fixture. No donor or biological group labels
are invented for this fixture.

## SenCID

`score_senescence` requires human gene-symbol raw counts and uses the original
published SID classifiers and recommendation model. The P0 adapter implements
exactly the upstream **non-denoising** path: select measured seneset genes, drop
all-zero genes, fill absent genes with 0.01, normalize seneset counts to CPM,
apply log2(CPM+1), and z-score across genes within each observation. The smoke
test compares all six scores against the original `ScaleData(denoising=False)`
and `Pred` functions. Only the unused eager import of DCA is omitted to run the
upstream reference without installing TensorFlow.

DCA denoising is not implemented in P0, and the tool reports `denoising=false`.
This choice is different from upstream's default denoising setting for single
cells. Zero measured seneset observations are rejected rather than scored from
imputed values alone. Missing-feature coverage is reported and checked.
Continuous scores are exported separately for each SID; the adapter does not
apply a universal binary senescence cutoff or convert scores to ages.

## Common contract

Input h5ad rows are observations and columns genes. CSV input has observations
in rows, a first ID column, and named gene columns. Model result tables retain
exact observation IDs and optional named metadata columns. Each request creates
a new result directory; no earlier results are overwritten. The output records
the backend and actual execution device. RNA age and SID scores must not be
combined into a single undifferentiated score or treated as independent evidence.

The published SenCID pickles record scikit-learn 0.24.2. The current worker smoke
loads them with the deployed scikit-learn runtime and compares to the original
upstream Python functions under that same runtime. This verifies wrapper parity
in the tested environment, not numerical equivalence to a separately restored
0.24.2 environment; scikit-learn's version warning is retained in the run log.
