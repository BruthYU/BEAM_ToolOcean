"""测试：模拟 LLM 发起请求，真实工具执行，再根据结果回答。"""

import anndata as ad
import numpy as np

from llm_call_demo import run_demo


def test_llm_tool_dispatch(tmp_path):
    file_path = tmp_path / "sample.h5ad"
    ad.AnnData(X=np.zeros((2, 3))).write_h5ad(file_path)

    result = run_demo(file_path, workspace=tmp_path / "workspace")

    assert result["tool_call"] == {
        "name": "read_h5ad",
        "arguments": {"file_path": str(file_path)},
    }
    assert result["tool_result"] == {"n_obs": 2, "n_vars": 3}
    assert result["answer"] == "这个文件有 2 个观测、3 个特征。"
