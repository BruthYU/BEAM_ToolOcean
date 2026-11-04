# 按官方方式新增一个 h5ad 工具

顺序就是：**写 JSON → 修改 default_config.py → 写工具类并加装饰器 → 运行 test**。

下面提供的是操作教程和示例文件。按步骤将示例放进仓库根目录下的 ToolUniverse 源码目录，即可完成接入。
所有 shell 命令均从本教程目录 `NewTool_tutorial/h5ad_reader_tutorial` 执行。在仓库根目录运行：

```bash
cd NewTool_tutorial/h5ad_reader_tutorial
```

**第一步：在 data 文件夹下写 JSON**

目标文件：`../../ToolUniverse/src/tooluniverse/data/h5ad_tools.json`。
这个文件声明工具的名字、用途和输入参数。内容如下，注意最外层是列表：

```json
[
  {
    "name": "read_h5ad",
    "type": "H5ADReaderTool",
    "description": "读取本地 h5ad 文件，返回观测数 n_obs 和特征数 n_vars。",
    "parameter": {
      "type": "object",
      "properties": {
        "file_path": {
          "type": "string",
          "description": "本地 h5ad 文件路径"
        }
      },
      "required": [
        "file_path"
      ],
      "additionalProperties": false
    },
    "cacheable": false
  }
]
```

示例已写在 [h5ad_tools.json](./h5ad_tools.json)，直接复制到目标位置：

```bash
cp ./h5ad_tools.json ../../ToolUniverse/src/tooluniverse/data/h5ad_tools.json
```

**第二步：修改 default_config.py，让框架知道这个 JSON 在哪里**

打开 [default_config.py](../../ToolUniverse/src/tooluniverse/default_config.py)，在 `default_tool_files = { ... }` 字典内部添加一项：

```python
"h5ad": os.path.join(current_dir, "data", "h5ad_tools.json"),
```

例如放在已有的 `"bar"` 那一行后面。这里的 `h5ad` 是类别名，后面加载工具时会用到。

**第三步：实现工具类，加上装饰器**

目标文件：`../../ToolUniverse/src/tooluniverse/h5ad_tool.py`。
内容如下：

```python
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
```

示例已写在 [h5ad_tool.py](./h5ad_tool.py)，复制到目标位置：

```bash
cp ./h5ad_tool.py ../../ToolUniverse/src/tooluniverse/h5ad_tool.py
```

装饰器里的 `H5ADReaderTool` 必须和 JSON 的 `type` 一致。
`run()` 接收调用参数字典，在这里读取文件并返回结果。配置已经写在 JSON 中，所以装饰器只需要写类型名。

**第四步：运行一个 test**

先在当前 Python 环境安装源码包和测试所需依赖：

```bash
python -m pip install -e ../../ToolUniverse 'anndata>=0.11' pytest
```

测试文件是 [test_h5ad_tool.py](./test_h5ad_tool.py)。它生成一个小 h5ad，通过刚才添加的类别加载工具，再检查结果：

```python
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
```

运行：

```bash
python -m pytest ./test_h5ad_tool.py -q
```

预期输出：`1 passed`。测试使用 `load_tools(categories=["h5ad"])`，会经过第二步添加的默认配置。

三个名字分别用于：

| 名字 | 写在哪里 | 用在哪里 |
| --- | --- | --- |
| `h5ad` | `default_config.py` 的字典键 | `load_tools(categories=["h5ad"])` |
| `read_h5ad` | JSON 的 `name` | 调用请求的 `name` |
| `H5ADReaderTool` | JSON 的 `type`、装饰器参数 | 框架查找实现类 |

本次已检查示例的语法和 JSON 配置。端到端 test 需要先完成前三步并安装依赖后执行。
