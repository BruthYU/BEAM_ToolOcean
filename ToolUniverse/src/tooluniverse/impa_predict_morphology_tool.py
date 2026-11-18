"""IMPA image perturbation prediction with explicit local pretrained assets."""

from .perturbation_image_generation import BaselineImageGenerationTool
from .tool_registry import register_tool


@register_tool("IMPATool")
class IMPATool(BaselineImageGenerationTool):
    backend = "impa"
