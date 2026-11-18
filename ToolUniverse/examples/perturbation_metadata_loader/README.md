# 扰动实验元数据加载工具

`perturbation_load_metadata` 专门处理表达及图像扰动实验的条件注释，兼容 `perturbation.load_metadata`。它将原始列映射为扰动、剂量、时间、对照、重复和批次等标准字段，验证实验记录的关联键，生成标准化文件与引用，并返回少量实际记录预览。

它不绑定某个 dataset，也不写回或自动关联表达矩阵、图像索引。下游工具通过 association_keys 关联实验数据。它不推断药物机制、控制标签或剂量单位，也不执行完整实验设计评估。通用表格读取在内部 metadata_io 模块复用；通用 data.load_metadata 尚未注册。

## 格式与依赖

支持 CSV、TSV（含 gzip）、XLSX、Parquet，以及 H5AD 的 obs。CSV/TSV 保留字符串标识和字面值 NA，不自动把它们当作缺失；标准字段中的空白按缺失处理。XLSX 读取第一个 sheet，也可指定 sheet_name。H5AD 直接读取 obs，不加载表达矩阵或 layers，但源文件完整性校验仍会扫描整个文件。

使用项目 Python 环境，在 ToolUniverse 目录安装源码：

```bash
python -m pip install -e .
```

Parquet 另需 pyarrow，XLSX 需要 openpyxl，H5AD obs 需要 anndata 和 h5py。依赖按格式导入，普通 CSV 不需要安装 anndata。

## 必填参数

| 参数 | 含义 |
| --- | --- |
| source | 本地实验元数据路径 |
| output_dir | 保存标准化文件及引用的目录 |
| experiment_level | cell、sample、well 或 field |
| field_map | 标准字段到原始列名的映射 |
| association_keys | 能唯一标识记录的标准字段组合 |

每个层级必须包含下列 ID 字段，且这些字段必须出现在 association_keys 中：

| 层级 | 必需 ID |
| --- | --- |
| cell | cell_id |
| sample | sample_id |
| well | plate_id、well_id |
| field | plate_id、well_id、site_id |

每条保留记录必须有 perturbation。关联键可额外包含 batch、time 等字段，但不得缺失、重复或指向未提供字段。H5AD 的观测索引以保留列 __obs_index__ 暴露。

## 图像扰动实验调用

下面使用 [BBBC021 官方图像注释](https://bbbc.broadinstitute.org/BBBC021) 的列名。替换 source 为本地文件路径：

```python
from tooluniverse import ToolUniverse

engine = ToolUniverse(load_workspace=False)
try:
    engine.load_tools(categories=["perturbation"], clear_existing=True)
    result = engine.run_one_function({
        "name": "perturbation.load_metadata",
        "arguments": {
            "source": "/data/BBBC021_v1_image.csv",
            "output_dir": "/tmp/beam-experiment-metadata",
            "experiment_level": "well",
            "field_map": {
                "plate_id": "Image_Metadata_Plate_DAPI",
                "well_id": "Image_Metadata_Well_DAPI",
                "perturbation": "Image_Metadata_Compound",
                "dose": "Image_Metadata_Concentration"
            },
            "association_keys": ["plate_id", "well_id"],
            "control_values": ["DMSO"],
            "duplicate_keys": "collapse_identical"
        }
    }, validate=True, use_cache=False)
    print(result)
finally:
    engine.close()
```

该表含每孔多个视野。显式 collapse_identical 只合并标准化后完全一致的记录；同一个键下扰动、剂量或其他选定注释有冲突时仍失败。未映射的原始列不参与一致性判断，因此应映射所有影响实验条件的字段。原始表保持不变。

本地真实验证得到 13200 条源记录、3300 条孔级记录，其中 DMSO 对照为 330 条。剂量单位未声明，因此返回 UNIT_UNSPECIFIED，不自动猜测单位。

## 表达扰动实验调用参数

```json
{
  "source": "/data/experiment.h5ad",
  "output_dir": "/tmp/beam-experiment-metadata",
  "experiment_level": "cell",
  "field_map": {
    "cell_id": "__obs_index__",
    "perturbation": "perturbation",
    "cell_type": "celltype"
  },
  "association_keys": ["cell_id"]
}
```

perturbation 只保留调用方声明的标签，不将空标签或任意特殊符号推断为对照。缺失标签默认报错；确认可排除未标注记录时显式设置 missing_perturbation=drop_rows，并检查 DROPPED_UNANNOTATED_ROWS 的数量。若全部被排除则失败。

## 可选规则

- constants：显式常量，如 `{"dose_unit": "uM"}`，不能覆盖 field_map 中的字段。
- control_values：精确匹配扰动标签。若源表也映射了 is_control，两种声明必须一致。
- is_control：可由源列或常量提供，接受 true/false、1/0、yes/no，不接受缺失或任意文本。
- missing_optional：默认 preserve，保留可选字段空值并报告数量；error 则拒绝。
- duplicate_keys：默认 error，只有显式 collapse_identical 才合并一致的重复记录。
- preview_rows：默认返回 5 条实际标准化记录，范围 0 到 20。
- metadata_id：可选输出名称，省略自动生成，同名不覆盖。

dose 和 time 必须为有限非负数值或明确缺失值；当前 time 表示非负实验时间，不支持负的相对时间偏移。单位以源值或常量保留，不做跨单位换算。没有对照声明时不生成 is_control，n_controls 为 null，并返回 CONTROL_UNSPECIFIED。

## 返回与引用

成功返回 status、data、warnings。data 包含 metadata_ref、table_path、层级、源记录数、输出记录数、标准字段、关联键、扰动类别数、对照数和 preview。preview 是实际标准化记录，不只是统计值。

metadata.jsonl 保存标准化记录，manifest.json 保存源文件与输出表的指纹、读取和字段映射规则、摘要及警告。未映射的原始字段留在源文件中，不全部复制进标准化输出。

```python
from tooluniverse.perturbation_metadata_tool import resolve_perturbation_metadata_ref

descriptor = resolve_perturbation_metadata_ref(result["data"]["metadata_ref"])
print(descriptor["table"]["path"], descriptor["summary"]["association_keys"])
```

引用解析检查源与标准化表是否改变，不自动打开或修改已有 dataset。下游应按显式键关联，禁止按行号拼接。大表读取和完整哈希校验仍有内存、磁盘与等待开销，不承诺固定时延。

构造与验收约定见 [工具规范](TOOL_SPEC.md)。测试目录、真实数据及运行产物仅保留本地。
