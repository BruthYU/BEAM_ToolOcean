"""Invoke the IMPA image generation tool with explicit pretrained assets."""

import argparse
import json
from pathlib import Path

from tooluniverse import ToolUniverse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path, required=True)
    parser.add_argument("--control-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument("--compound", action="append", default=None)
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--allow-unseen-compounds", action="store_true")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = parser.parse_args()
    root = args.asset_dir.resolve()
    arguments = {
        "control_manifest": str(args.control_manifest.resolve()),
        "compounds": args.compound or ["floxuridine"],
        "model_dir": str(root / "IMPA"),
        "output_dir": str(args.output_dir.resolve()),
        "python_executable": args.python,
        "device": args.device,
        "samples_per_condition": args.samples,
        "batch_size": 1,
    }
    arguments.update(
        {
            "checkpoint_path": str(root / "weights/impa/bbbc021/nets.ckpt"),
            "embedding_checkpoint": str(root / "weights/impa/bbbc021/embeddings.ckpt"),
            "embeddings_path": str(root / "emb_fp.csv"),
            "allow_unseen_compounds": args.allow_unseen_compounds,
        }
    )
    tu = ToolUniverse()
    tu.load_tools(categories=["image_perturbation_models"])
    result = tu.run_one_function(
        {
            "name": "IMPA_predict_morphology",
            "arguments": arguments,
        },
        use_cache=False,
    )
    print(json.dumps(result, indent=2))
    if not isinstance(result, dict) or result.get("status") != "success":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
