"""Evaluate perturbation images with a frozen MoA classifier."""

from .image_perturbation_metrics_tool import ImagePerturbationEvaluationTool
from .tool_registry import register_tool


@register_tool("ImagePerturbationEvaluateMoATool")
class ImagePerturbationEvaluateMoATool(ImagePerturbationEvaluationTool):
    """Evaluate perturbation images with a frozen MoA classifier."""

    operation = "moa"
