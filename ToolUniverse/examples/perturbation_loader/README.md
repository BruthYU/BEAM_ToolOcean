# 扰动表达数据加载工具

`perturbation_load_dataset` 读取本地表达数据，检查矩阵和标识，按显式键关联实验元数据，并返回 JSON 摘要及持久化引用。兼容别名为 `perturbation.load_dataset`。它不执行标准化、过滤或预测，也不把完整矩阵返回给调用方。

## 支持格式

| format | 输入 | 约束 |
| --- | --- | --- |
| `h5ad` | H5AD 文件 | 支持 X、raw/X、layers/名称，默认 backed 读取 |
| `10x_h5` | 现代 10x H5 | 必须有 /matrix 布局 |
| `10x_mtx` | 10x MTX 目录 | 直接包含 matrix、features 或 genes、barcodes，支持 gzip |
| `csv`、`tsv` | 数值表达表 | 支持 gzip，必须指定索引列和矩阵方向 |

省略 format 时自动识别。旧版多基因组 H5、混合模态 10x、图像、表型特征及 Zarr 不属于当前支持范围。10x 保持稀疏矩阵并统一为观测 × 特征；CSV/TSV 为内存读取。

## 安装与调用

在 ToolUniverse 项目目录使用现有 Python 环境安装源码及可选读取依赖：

```bash
python -m pip install -e . anndata scipy h5py
```

替换示例中的源路径后执行：

```python
import json
from tooluniverse import ToolUniverse

engine = ToolUniverse(load_workspace=False)
try:
    engine.load_tools(categories=["perturbation"], clear_existing=True)
    result = engine.run_one_function({
        "name": "perturbation_load_dataset",
        "arguments": {
            "source": "/data/expression.h5ad",
            "output_dir": "/tmp/beam-expression-datasets"
        }
    }, validate=True, use_cache=False)
    print(json.dumps(result, ensure_ascii=False, indent=2))
finally:
    engine.close()
```

也可用 `from tooluniverse.tools import perturbation_load_dataset` 调用 SDK 包装入口。

## 参数

必填参数为 source 与 output_dir。dataset_id 可选，省略时自动生成；同名输出不覆盖。format 默认为 auto。H5AD 的 matrix_key 默认为 X，可选 raw/X 或 layers/名称；backed 默认 true，其他格式不得显式传入 backed。

CSV/TSV 附加参数示例：

```json
{
  "source": "/data/expression.csv",
  "output_dir": "/tmp/beam-expression-datasets",
  "index_column": "cell_id",
  "matrix_orientation": "observations_by_features"
}
```

矩阵也可指定为 features_by_observations。除索引列外必须全部为有限数值，观测和基因标识必须唯一且非空，不自动猜测方向。

obs_column_map 是标准字段到原始观测列名的映射，如 `{"perturbation": "condition"}`。可关联外部 CSV/TSV/Parquet 元数据，Parquet 另需 pyarrow。关联需提供 metadata_path、metadata_join_keys，例如 `{"data": ["__obs_index__"], "metadata": ["cell_id"]}`；默认严格匹配，禁止按行号拼接或一对多展开。

## 返回与数据引用

成功返回 status、data、warnings。data 包含 dataset_ref、dataset_id、格式、shape、dtype、矩阵位置及 metadata_fields；metadata_fields 仅为字段名，不包含具体数据记录。

dataset_ref 指向 manifest.json，记录源文件绝对路径、读取规则、字段映射及各输入文件的指纹。后续进程可验证并解析：

```python
from tooluniverse.perturbation_dataset_tool import resolve_dataset_ref

descriptor = resolve_dataset_ref(result["data"]["dataset_ref"])
print(descriptor["inputs"], descriptor["read_options"])
```

解析返回描述，不返回打开的 AnnData。外部元数据不会写回原始文件；下游按描述中的键和映射重建关联。源被移动或修改后需要重新加载。

## 资源与验证范围

H5AD backed 模式可将主矩阵保留在磁盘，obs 和部分 layers 仍可能进入内存。当前加载前后各计算一次完整 SHA256，解析引用时再次计算；没有内部固定超时或后台任务接口，客户端可能另有超时限制。

本地验证覆盖真实 H5AD、10x H5、MTX、CSV/TSV 转换副本、gzip、矩阵数值、元数据关联和错误输入。一次约 12.27 GiB、4705402 × 62710 的 H5AD backed 实测成功，工具约 42.15 秒、进程峰值 RSS 约 1.70 GiB；结果依赖运行环境，不代表后续完整矩阵分析的资源需求。测试目录、数据文件和运行产物不随本工具提交。

完整构造、参数、错误和验收要求见 [工具规范](TOOL_SPEC.md)。

## 两个 dataloader 的组合调用

本 PR 提供表达数据和图像数据两个通用入口，通过 `format="auto"` 自动识别支持的格式，也可通过 `format` 显式选择；不按文件后缀拆分工具名。

| 工具名 | 实现类（JSON type） | 分类 |
|---|---|---|
| `perturbation_load_dataset` | `PerturbationDatasetTool` | `perturbation` |
| `image_perturbation_load_images` | `PerturbationImageDatasetTool` | `image_perturbation` |

```python
from tooluniverse import ToolUniverse

engine = ToolUniverse(load_workspace=False)
try:
    engine.load_tools(categories=["perturbation", "image_perturbation"])
    expression = engine.run_one_function({
        "name": "perturbation_load_dataset",
        "arguments": {
            "source": "/data/expression.h5ad",
            "output_dir": "/tmp/expression-datasets",
            "format": "auto",
        },
    }, validate=True, use_cache=False)
    images = engine.run_one_function({
        "name": "image_perturbation_load_images",
        "arguments": {
            "source": "/data/images",
            "source_type": "directory",
            "file_pattern": "*.png",
            "output_dir": "/tmp/image-datasets",
            "format": "auto",
        },
    }, validate=True, use_cache=False)
finally:
    engine.close()
```

完整图像参数、可选读取依赖及持久化引用说明见 [图像加载说明](../image_perturbation_loader/README.md)。这两个入口独立加载各自数据，不自动匹配细胞与图像。
