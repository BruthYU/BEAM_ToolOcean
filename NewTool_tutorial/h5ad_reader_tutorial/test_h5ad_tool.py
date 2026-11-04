"""完成教程前三步后运行，验证默认配置能够找到工具。"""

import anndata as ad
import numpy as np
from tooluniverse import ToolUniverse


def test_read_h5ad(tmp_path):
    # 准备一个 2 个观测、3 个特征的测试文件。
    file_path = tmp_path / "sample.h5ad"
    ad.AnnData(X=np.zeros((2, 3))).write_h5ad(file_path)

    tu = ToolUniverse(workspace=str(tmp_path / "workspace"), load_workspace=False)
    try:
        # h5ad 对应 default_config.py 中新增的键。
        tu.load_tools(categories=["h5ad"], clear_existing=True)

        # read_h5ad 对应 JSON 中的 name。
        result = tu.run_one_function({
            "name": "read_h5ad",
            "arguments": {"file_path": str(file_path)},
        })
        assert result == {"n_obs": 2, "n_vars": 3}
    finally:
        tu.close()
