"""模拟 LLM 的两次回复，使用真实 ToolUniverse 执行 h5ad 工具。"""

import json
import tempfile
from pathlib import Path


def mock_llm(request):
    """固定场景的模拟函数：第一次要求调用工具，第二次根据结果回答。"""
    if "tool_result" not in request:
        # 本例只有一个工具，固定选它；这里没有真实的语言理解。
        return {
            "tool_call": {
                "name": request["tools"][0]["name"],
                "arguments": {"file_path": request["file_path"]},
            }
        }

    result = request["tool_result"]
    if result.get("status") == "error":
        return {"answer": f"读取失败：{result['error']}"}
    return {"answer": f"这个文件有 {result['n_obs']} 个观测、{result['n_vars']} 个特征。"}


def run_demo(file_path, workspace):
    from tooluniverse import ToolUniverse

    tu = ToolUniverse(workspace=str(workspace), load_workspace=False)
    try:
        # 第一步：从 default_config.py 加载 h5ad 类别。
        tu.load_tools(categories=["h5ad"], clear_existing=True)
        if "read_h5ad" not in tu.all_tool_dict:
            raise RuntimeError("请先按教程将 h5ad 的 JSON、类别映射和实现接入源码。")

        # 第二步：取出工具名称、描述、参数，作为模拟 LLM 的输入。
        config = tu.all_tool_dict["read_h5ad"]
        tool_spec = {key: config[key] for key in ("name", "description", "parameter")}
        request = {
            "question": "请读取这个 h5ad，告诉我观测数和特征数。",
            "file_path": str(file_path),
            "tools": [tool_spec],
        }
        print("[1] 用户问题：", request["question"])
        print("[2] 给模拟 LLM 的工具描述：", json.dumps(tool_spec, ensure_ascii=False))

        # 第三步：模拟 LLM 生成 name + arguments。
        response = mock_llm(request)
        tool_call = response["tool_call"]
        print("[3] 模拟 LLM 请求调用：", json.dumps(tool_call, ensure_ascii=False))

        # 第四步：由 ToolUniverse 查找工具、校验参数并执行 run()。
        tool_result = tu.run_one_function(tool_call, use_cache=False, validate=True)
        print("[4] ToolUniverse 执行结果：", json.dumps(tool_result, ensure_ascii=False))

        # 第五步：把执行结果交回模拟 LLM，得到最终回答。
        request["tool_result"] = tool_result
        final_response = mock_llm(request)
        print("[5] 模拟 LLM 最终回答：", final_response["answer"])
        return {
            "tool_call": tool_call,
            "tool_result": tool_result,
            "answer": final_response["answer"],
        }
    finally:
        tu.close()


if __name__ == "__main__":
    import anndata as ad
    import numpy as np

    # 自动生成小文件，演示结束后清理，不需要准备真实数据。
    with tempfile.TemporaryDirectory(prefix="tu-llm-demo-") as directory:
        root = Path(directory)
        file_path = root / "sample.h5ad"
        ad.AnnData(X=np.zeros((2, 3))).write_h5ad(file_path)
        run_demo(file_path, workspace=root / "workspace")
