from .base_tool import BaseTool
from .tool_registry import register_tool


@register_tool("H5ADReaderTool")
class H5ADReaderTool(BaseTool):
    def run(self, arguments):
        from anndata.io import read_h5ad

        adata = read_h5ad(arguments["file_path"], backed="r")
        try:
            return {
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
            }
        finally:
            adata.file.close()
