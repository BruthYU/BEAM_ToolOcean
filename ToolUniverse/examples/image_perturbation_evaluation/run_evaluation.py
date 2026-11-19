"""Call independent image evaluation tools for any named generator."""

import argparse
import json
from pathlib import Path

from tooluniverse import ToolUniverse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metric", choices=["moa", "fid_kid"], required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--asset-dir", type=Path, required=True)
    parser.add_argument("--reference-metadata", type=Path)
    parser.add_argument("--reference-images", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument(
        "--profile",
        choices=[
            "bbbc021_standard",
            "bbbc021_unseen_compounds",
            "bbbc021_smoke",
            "bbbc021_native_inception",
        ],
        default="bbbc021_standard",
    )
    parser.add_argument("--split", default="test")
    parser.add_argument("--variant", choices=["standard", "ood"], default="standard")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--kid-subset-size", type=int, default=100)
    parser.add_argument("--kid-subsets", type=int, default=100)
    args = parser.parse_args()
    arguments = {
        "predictions": str(args.predictions.resolve()),
        "model_name": args.model_name,
        "evaluation_profile": args.profile,
        "split": args.split,
        "output_dir": str(args.output_dir.resolve()),
        "python_executable": args.python,
        "device": args.device,
    }
    assets = args.asset_dir.resolve()
    if args.metric == "moa":
        if not args.reference_metadata:
            parser.error("--reference-metadata is required for MoA.")
        arguments.update(
            {
                "reference_metadata": str(args.reference_metadata.resolve()),
                "evaluator_checkpoint": str(assets / f"moa_{args.variant}.pth"),
                "label_map": str(assets / "bbbc021_labels.json"),
            }
        )
    else:
        if not args.reference_images:
            parser.error("--reference-images is required for FID/KID.")
        arguments.update(
            {
                "reference_images": str(args.reference_images.resolve()),
                "feature_extractor_weights": str(
                    assets / f"inception_from_moa_{args.variant}.pth"
                ),
                "kid_subset_size": args.kid_subset_size,
                "kid_subsets": args.kid_subsets,
            }
        )
    tu = ToolUniverse()
    tu.load_tools(categories=["image_perturbation_evaluation"])
    result = tu.run_one_function(
        {
            "name": {
                "moa": "ImagePerturbation_evaluate_moa",
                "fid_kid": "ImagePerturbation_calculate_fid_kid",
            }[args.metric],
            "arguments": arguments,
        },
        use_cache=False,
    )
    print(json.dumps(result, indent=2))
    if not isinstance(result, dict) or result.get("status") != "success":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
