# 分步骤模拟 LLM 调用 ToolUniverse 工具

继续使用前面的 `read_h5ad` 工具，模拟这个过程：

```text
用户提出问题
  → 给 LLM 看工具描述
  → LLM 返回工具名和参数
  → 主程序交给 ToolUniverse 执行
  → 把执行结果交回 LLM
  → LLM 给用户最终回答
```

本例用 `mock_llm()` 模拟两次模型回复，不需要 API key；ToolUniverse 的调度和 h5ad 读取都使用真实代码。模拟函数针对一个固定问题选择工具，使用简化字典表示回复。

完整代码：[llm_call_demo.py](llm_call_demo.py)。测试：[test_llm_call.py](test_llm_call.py)。

**准备：先接入上一例的 h5ad 工具**

下面所有命令都在本篇教程目录执行：

```bash
cd /mnt/shared-storage-gpfs2/beam-gpfs02/yulang/master/CELL/BEAM_ToolOcean/NewTool_tutorial/llm_call_tutorial
```

若尚未完成前面的注册步骤，先复制两个文件：

```bash
cp ../h5ad_reader_tutorial/h5ad_tools.json ../../ToolUniverse/src/tooluniverse/data/h5ad_tools.json
cp ../h5ad_reader_tutorial/h5ad_tool.py ../../ToolUniverse/src/tooluniverse/h5ad_tool.py
```

然后在 [default_config.py](../../ToolUniverse/src/tooluniverse/default_config.py) 的 `default_tool_files` 字典中添加这一项，已有则跳过：

```python
"h5ad": os.path.join(current_dir, "data", "h5ad_tools.json"),
```

安装依赖：

```bash
python -m pip install -e ../../ToolUniverse 'anndata>=0.11' pytest
```

**第一步：主程序加载工具**

```python
tu.load_tools(categories=["h5ad"], clear_existing=True)
```

框架根据 `default_config.py` 找到 `h5ad_tools.json`，把其中的工具配置放入当前实例。此时可以用 `tu.all_tool_dict["read_h5ad"]` 取得配置。

**第二步：把工具描述交给模拟 LLM**

```python
config = tu.all_tool_dict["read_h5ad"]
tool_spec = {key: config[key] for key in ("name", "description", "parameter")}
request = {
    "question": "请读取这个 h5ad，告诉我观测数和特征数。",
    "file_path": str(file_path),
    "tools": [tool_spec],
}
```

工具描述告诉模型：有什么工具、什么时候使用、需要什么参数。这里把文件路径也作为用户提供的上下文交给模型。

**第三步：模拟 LLM 生成调用请求**

```python
response = mock_llm(request)
tool_call = response["tool_call"]
```

第一次回复的内容类似：

```json
{
  "tool_call": {
    "name": "read_h5ad",
    "arguments": {"file_path": "/tmp/example/sample.h5ad"}
  }
}
```

这里还没有读取文件，只生成了工具名称和参数。

**第四步：主程序将请求交给 ToolUniverse 执行**

```python
tool_result = tu.run_one_function(tool_call, use_cache=False, validate=True)
```

这一行触发的核心过程是：

```text
name = read_h5ad
  → 查到 JSON 配置中的 type = H5ADReaderTool
  → 从注册表找到类，创建或复用实例
  → 根据 parameter 校验 file_path
  → 调用实例.run(arguments)
  → 读取 h5ad，返回 {"n_obs": 2, "n_vars": 3}
```

对应源码集中在 [execute_function.py](../../ToolUniverse/src/tooluniverse/execute_function.py)：入口是 `run_one_function()`，实例化在 `init_tool()`，最终调用发生在 `_execute_tool_with_stream()`。

**第五步：把工具结果交回模拟 LLM**

```python
request["tool_result"] = tool_result
final_response = mock_llm(request)
print(final_response["answer"])
```

模拟函数第二次被调用时，根据工具实际返回的数字生成回答：

```text
这个文件有 2 个观测、3 个特征。
```

这两次模型调用由我们的主程序组织。ToolUniverse 负责中间的工具执行；读取完成后，主程序把结果继续交给模型。

**第六步：运行演示和 test**

```bash
python llm_call_demo.py
python -m pytest test_llm_call.py -q
```

演示脚本自动生成一个 2 行、3 列的小 h5ad，并依次打印用户问题、工具描述、调用请求、工具结果和最终回答。

test 检查三件事：调用的工具和参数正确、真实读取结果正确、最终回答使用了返回结果。预期输出为 `1 passed`。

本次已检查 Python 语法，并独立验证模拟 LLM 的请求生成和结果回答；真实 SDK 与 h5ad 的端到端执行，需先完成上面的工具接入和依赖安装。
