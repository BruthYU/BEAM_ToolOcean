# 扰动表达数据加载工具

`perturbation.load_dataset` 读取本地表达数据、检查矩阵与标识、关联实验注释，并保存可供后续工具使用的数据引用。它返回 JSON 摘要，不返回完整 AnnData 或表达矩阵，不执行标准化、质控过滤或预测。

## 支持格式

| 输入 | format | 要求 |
| --- | --- | --- |
| H5AD | `h5ad` | 可选 `X`、`raw/X`、`layers/<name>`，默认 backed 读取 |
| 现代 10x H5 | `10x_h5` | 使用 `/matrix` 布局 |
| 10x MTX 目录 | `10x_mtx` | 直接包含 matrix、features 或 genes、barcodes 文件，支持 gzip |
| CSV 或 TSV | `csv` 或 `tsv` | 支持 gzip；必须指定索引列与矩阵方向 |

`format` 默认 `auto`。10x 保持稀疏矩阵并统一为观测 × 特征方向。旧版多基因组 10x H5、混合模态 10x、图像、表型特征、Zarr、H5MU 和 GCTX 不在当前支持范围。

接口定义见 [工具规范](TOOL_SPEC.md)。

## 环境

从仓库根目录执行。使用项目 Python 环境安装源码及表达读取依赖：

```bash
python -m pip install -e ./ToolUniverse anndata scipy h5py pandas requests pytest
```

如需外部 Parquet 元数据，再安装 `pyarrow`。下面的 `PYTHONPATH` 确保使用本仓库源码，避免误用旧的已安装版本。

## 调用

最少参数为 `source` 和 `output_dir`。文件不存在时先按下一节下载真实测试数据。

```bash
PYTHONPATH=ToolUniverse/src python - <<'PY'
import json
from tooluniverse import ToolUniverse

engine = ToolUniverse(load_workspace=False)
try:
    engine.load_tools(categories=["perturbation"], clear_existing=True)
    result = engine.run_one_function({
        "name": "perturbation.load_dataset",
        "arguments": {
            "source": "perturb_tools/test_data/real/molecular/adamson_2016.h5ad",
            "output_dir": "/tmp/beam-expression-datasets",
            "obs_column_map": {"perturbation": "perturbation"}
        }
    }, validate=True, use_cache=False)
    print(json.dumps(result, ensure_ascii=False, indent=2))
finally:
    engine.close()
PY
```

也可使用 Python SDK 包装入口 `from tooluniverse.tools import perturbation_load_dataset`；它映射到同一个带点号的工具名。

10x MTX 的 `source` 应指向 `perturb_tools/test_data/real/format_reference/pbmc3k_mtx/filtered_gene_bc_matrices/hg19`。CSV/TSV 另外需要：

```json
{
  "source": "/data/expression.csv",
  "output_dir": "/tmp/beam-expression-datasets",
  "index_column": "cell_id",
  "matrix_orientation": "observations_by_features"
}
```

表格除索引列外必须全部为有限数值。输入为基因 × 细胞时使用 `features_by_observations`。大规模稀疏数据优先使用 H5AD 或 MTX；CSV/TSV 为内存读取。

## 数据引用

成功结果包含 `status`、`data`、`warnings`。`data` 包含 `dataset_ref`、数据集 ID、格式、维度、dtype、矩阵位置及映射字段名。字段名列表不等于实际元数据记录。

`dataset_ref` 指向 `manifest.json`，保存源文件绝对路径、读取选项和各文件指纹。`resolve_dataset_ref` 验证源未改变并返回描述，不打开矩阵，也不将外部注释写回原文件：

```python
from tooluniverse.perturbation_dataset_tool import resolve_dataset_ref

descriptor = resolve_dataset_ref(result["data"]["dataset_ref"])
print(descriptor["inputs"], descriptor["read_options"])
```

输出目录中同名数据集不能覆盖；省略 `dataset_id` 时自动生成名称。原数据移动或修改后，引用需要重新生成。

## 真实数据准备与验证

下载脚本获取三个原始公开文件，不生成合成观测，也不会在 pytest 运行时自动联网：

```bash
python perturb_tools/download_expression_data.py
```

数据保存在 `test_data/real/`，不提交 Git。`expression_sources.json` 记录原始 URL、大小和本地 SHA256；Adamson 文件另外核对发布方 MD5。下载内容包括：

| 数据 | 用途 |
| --- | --- |
| Adamson 2016 CRISPRi，5768 × 35635 | H5AD 结构、矩阵选择和真实表达值 |
| 10x PBMC 1k，1222 × 33538 | 现代 10x H5，作为格式测试 |
| 10x PBMC3k，2700 × 32738 | MTX 与标签，作为格式测试 |

来源：[scPerturb 发布记录](https://zenodo.org/records/13350497)。两个 10x 原始下载 URL 保存在下载脚本和生成的来源清单中。PBMC 不作为扰动预测评价数据。

CSV/TSV 测试使用真实 Adamson 子集的格式转换副本；gzip MTX 测试压缩真实 10x 文件。临时关联标签仅用于校验关联行为，不作为新的生物实验标签。这些测试仅在本地开发环境执行，不随本工具 PR 提交。

本次 PR 不包含测试目录或本地数据。下载完成后，可用上面的工具调用示例检查 `status == "success"` 和 `shape == [5768, 35635]`。本地开发测试另外覆盖格式转换、元数据关联、异常输入和 SDK 映射。

## 大文件行为

H5AD 默认 `backed=True`，主矩阵可保留在磁盘上，但 `obs`、部分 layers 和其他注释仍可能进入内存。不要将 backed 理解为恒定内存读取。

当前加载前后各计算一次完整 SHA256；解析引用时再次计算。它没有内部固定超时，也没有进度或后台任务接口，客户端可能有独立超时。一次本地测试用 12.27 GiB、4705402 × 62710 的 H5AD 成功完成：工具约 42.15 秒、进程峰值 RSS 约 1.70 GiB，其中两次哈希约 30.62 秒。该结果依赖当时的共享存储和运行环境，不保证其他机器相同，也不代表后续完整矩阵分析的资源需求。
