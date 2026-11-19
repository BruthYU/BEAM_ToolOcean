"""Call CellFlux through ToolUniverse with explicitly prepared local assets."""

import argparse
import json
from pathlib import Path

from tooluniverse import ToolUniverse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--python", required=True, help="Model environment Python")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--compound", action="append", default=None)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--samples", type=int, default=1)
    args = parser.parse_args()
    tu = ToolUniverse()
    tu.load_tools(categories=["cellflux"])
    result = tu.run_one_function(
        {
            "name": "CellFlux_predict_morphology",
            "arguments": {
                "control_manifest": str((args.data_dir / "controls.csv").resolve()),
                "model_dir": str((args.asset_dir / "CellFlux").resolve()),
                "checkpoint_path": str((args.asset_dir / "checkpoint.pth").resolve()),
                "embeddings_path": str((args.asset_dir / "emb_fp.csv").resolve()),
                "python_executable": args.python,
                "output_dir": str(args.output_dir.resolve()),
                "compounds": args.compound or ["floxuridine"],
                "samples_per_condition": args.samples,
                "batch_size": 1,
                "sampling_steps": args.sampling_steps,
                "device": args.device,
            },
        },
        use_cache=False,
    )
    print(json.dumps(result, indent=2))
    if not isinstance(result, dict) or result.get("status") != "success":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
