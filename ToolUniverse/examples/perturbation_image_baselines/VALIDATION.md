# IMPA integration checks

The SDK ran the released pretrained weights on original BBBC021 DMSO cells.
A CPU run generated two 96×96×3 uint8 PNGs, and both independent evaluation
tools consumed the returned image index. The 26 training condition names were
verified by exact matches against the frozen fingerprint rows.

Calling the canonical and legacy SDK names with the same input and seed
produced identical PNG hashes. Two requests with caching enabled produced
separate runs. These checks used a small integration subset.

Local evidence is available under `perturb_tools/test_runs/image_baselines/`,
`perturb_tools/test_runs/tool_batch_quality/` and
`perturb_tools/test_runs/morphology_name_validation/`.
