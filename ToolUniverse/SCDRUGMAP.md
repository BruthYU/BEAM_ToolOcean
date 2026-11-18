# scDrugMap 队列和单细胞表达读取

BEAM_ToolOcean 的 `scdrugmap` 分类包含六个只读接口，共享一个 ToolUniverse 后端。每个 `dataset_id` 对应作者的一个分组，例如 `primary:Data1` 与 `external:Data1` 是两组不同数据。同一 GEO 可以按疗法或取样时间出现多次，先按 `accession` 查目录，再选分组。

| 接口 | 输入 | 返回 |
|---|---|---|
| `scDrugMap_list_datasets` | 可选目录过滤、分页 | 分组 ID、研究、癌种、组织、方案、平台、细胞与基因数、已核实样本与患者数 |
| `scDrugMap_get_dataset_info` | `dataset_id` | 标签分布和单位、临床映射覆盖、原始矩阵方向和尺度、源文件名与校验值 |
| `scDrugMap_query_samples` | `dataset_id`、可选临床过滤、分页 | 每份已回连活检一行，含患者、时间点、细胞数及原始结局；缺少映射时返回 `sample_mapping_unavailable` |
| `scDrugMap_query_cells` | `dataset_id`、可选细胞／临床过滤、分页 | 作者的细胞 ID、原标签，以及可用的患者、活检和临床信息 |
| `scDrugMap_list_genes` | `dataset_id`、可选基因 ID 前缀、分页 | 源文件中的确切基因 ID 和列／行位置，前缀按字面匹配 |
| `scDrugMap_get_expression` | `dataset_id`、1–64 个 `genes`、1–128 个 `cell_ids` | 按请求顺序返回基因 × 细胞数值、轴名、表达尺度及缺失值个数 |

分页大小为 1–200，下一页由 `next_offset` 指定。过滤字段接受精确值列表，字段内取并集，字段间取交集。成功响应为 `status=success,data=…`；查询失败时返回 `status=error,code,error`。未知基因或条码会得到错误。源表达中的缺失值返回 JSON null，数值零按原值返回。完整矩阵留在服务器，路径不作为 Agent 调用参数。

## 从非 coding Agent 调用

先通过以下两次调用选定队列并检查独立样本数。将第一步得到的具体 ID 传给后续工具。

```json
{"name":"scDrugMap_list_datasets","arguments":{"filters":{"accession":["GSE120575"]}}}
```

```json
{"name":"scDrugMap_query_samples","arguments":{"dataset_id":"external:Data20","filters":{"timepoint":["pre"]},"limit":50}}
```

该分组为抗 PD-1；治疗前有 12 位患者各一份活检，对应 3,795 个细胞。4 份病灶应答，8 份未应答。需要读取表达时，先用 `query_cells` 获取条码，再传入基因和细胞 ID。下面这个独立读取示例可用于核对轴顺序。

```json
{"name":"scDrugMap_get_expression","arguments":{"dataset_id":"external:Data20","genes":["NKG7","TCF7"],"cell_ids":["A10_P4_M15","A10_P1_M15"]}}
```

预期数值为 `[[12.46,10.38],[0,2.67]]`。这些数值是作者发布的表达量，按原精度读取，不作 UMI、TPM 或统一归一化的假定。

## 数据范围与临床口径

- 目录含主集 36 组、外部集 24 组，归档共有 665,602 条组内细胞记录。跨组可能重用同一份表达矩阵。这个总数不能作为独立患者数或论文抽样后的 495,237 个细胞使用。
- 36 张主集矩阵按细胞 × 基因排列。24 张外部矩阵按基因 × 细胞排列，细胞标签来自另一个文件。返回值统一为基因 × 细胞；公开条码以注释为准，连接符差异只在逐项唯一匹配后用于关联。
- 临床联结覆盖两个 GEO 的四个分组，共 54,965 条细胞记录。GSE120575 三个治疗组有患者重叠，合计 32 位不同患者、48 份样本；原始 R/NR 定义于取样病灶，`post` 在这个研究中表示治疗中取样。GSE123813 的归档分组有 8 位有前后配对的患者、18 份活检；`response_original` 保留 Yes/No，`sample_regimen` 保留回连的实际方案。其他 56 组的患者、样本数以 null 表示未核实，仍可读取细胞、基因与表达。
- GSE247684 的 `external:Data6` 包含 673 条 `persistent`。GSE243649 的 `external:Data9` 与 `external:Data10` 在归档中重用同一表达矩阵，进行留出评测时需按来源绑定。

数据出处：[作者归档](https://zenodo.org/records/17561695) · [scDrugMap 论文](https://www.nature.com/articles/s41467-025-67481-2) · [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) · [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813) 。GSE123813 的临床响应来源记录在分组详情中。

## 注册、索引与复跑

服务端程序位于 `src/tooluniverse/scdrugmap_tool.py`，注册配置位于 `src/tooluniverse/data/scdrugmap_tools.json`，默认分类在 `default_config.py` 中登记。`ToolUniverse.load_tools(categories=["scdrugmap"])` 可加载这组工具，MCP 服务可向 Agent 开放相同分类。

在支持标准 MCP stdio 的客户端添加以下启动项，并将 `command` 替换为已安装该仓库 ToolUniverse 的 Python 解释器。

```json
{
  "mcpServers": {
    "scdrugmap": {
      "command": "/absolute/path/to/venv/bin/python",
      "args": ["-c", "from tooluniverse.smcp_server import run_stdio_server; run_stdio_server()", "--categories", "scdrugmap"]
    }
  }
}
```

索引的 `registry.sqlite` 保存分组、基因、细胞、临床映射及表达源行的字节位置。运行时以只读方式打开索引和矩阵，并核对源矩阵的大小、修改时间以及请求行的 ID 和宽度。构建时记录矩阵 SHA-256。移动或修改表达文件后需用新输出目录重建索引并切换配置。索引不会复制整库表达值。

从 BEAM 仓库源码运行时，默认索引位于仓库上级目录下的 `dataset/scDrugMap/query_data`。其他部署可通过 `TOOLUNIVERSE_SCDRUGMAP_DATA` 指向索引目录。构建程序只需 Python 标准库。目录描述随代码提供，两个 GEO 文件是按固定 SHA-256 校验的可选临床输入。

```sh
python ToolUniverse/scripts/prepare_scdrugmap.py \
  --source /path/to/dataset/scDrugMap \
  --output /path/to/new/query_data \
  --geo /path/to/GSE120575_patient_ID_single_cells.txt.gz \
  --bcc /path/to/GSE123813_bcc_all_metadata.txt.gz
TEST_SCDRUGMAP_DATA=/path/to/new/query_data python ToolUniverse/tests/integration/test_scdrugmap_reader.py
```

## 还需补哪些分析能力

以下为拟开发接口。公开读取工具面向数据发现与核对；批量分析接口应在服务端接收分组、过滤条件或上一步返回的队列 ID，直接运算，避免让 Agent 分页传递几十万个条码。

| 顺序 | 分析能力 | 输入与输出 | 首次验收 |
|---|---|---|---|
| 1 | 构建分析队列和配对表 | 选择研究、方案、时间窗、终点与单位，返回版本化 `cohort_id`、纳入／排除病例、配对表、重用数据和缺失信息 | GSE120575 抗 PD-1 基线 12 位；GSE123813 归档内 8 位有前后配对。临床终点不清或仅有细胞系时说明可用任务 |
| 2 | 计算细胞状态并汇总到样本 | 对指定细胞群计算标志物检出比例、患者内表达程序均值或已有注释模型的结果，返回逐活检特征、分子、分母、基因覆盖和测量版本 | 先把单 GSE 原型的比例与程序分数接到通用队列 ID；固定门控下复核 TCF7。调用 CellTypist 等模型前检查输入尺度、全基因覆盖 |
| 3 | 检验基线关联与前后变化 | 基线病例间比较和同患者配对比较使用各自的统计单位；输出效应量、不确定性、留一敏感性和多重校正；保留重复活检处理规则 | 基线在 GSE120575 复算；纵向在 GSE123813 按患者算变化。治疗前后配对指患者，不能由两次活检认定同一个细胞 |
| 4 | 统一疗效预测和盲评 | 冻结模型接收基线 `cohort_id` 与 `model_id`，返回逐病例分数。划分与测试标签由评测端托管，按患者或整研究留出，比较简单基线、单细胞聚合模型和预定义标志物 | 同一病例全部细胞在同一折；GSE243649 的重复矩阵绑定；报告 AUROC、AUPRC、校准和不确定性，注明病灶或患者终点 |
| 5 | 评测虚拟细胞的状态变化预测 | 模型接收基线状态与受支持的治疗条件，预测细胞群状态分布；与留出的真实治疗中／后样本比较，再测其对临床预测的增量 | 先确认模型支持相应药物、细胞背景与时点，记录取材差异；已有扰动文件读取器可复用，影像形态分数和病灶应答分开评测 |

工具负责读取、测量、预测和统计。上层 Agent 或研究者给出科学问题，并据输出解释证据。
