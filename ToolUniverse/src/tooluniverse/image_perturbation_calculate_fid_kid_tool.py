"""Calculate FID and KID for generated and reference images."""

from .image_perturbation_metrics_tool import ImagePerturbationEvaluationTool
from .tool_registry import register_tool


@register_tool("ImagePerturbationCalculateFIDKIDTool")
class ImagePerturbationCalculateFIDKIDTool(ImagePerturbationEvaluationTool):
    """Calculate FID and KID for generated and reference images."""

    operation = "fid_kid"
