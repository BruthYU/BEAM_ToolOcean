# 扰动表达数据加载工具规范

本规范定义 `perturbation.load_dataset` 的已实现接口及验收要求，接口版本为 0.1。图像加载、表型加载、数据预览和自动分析不属于该工具。

## 注册与配置

实现继承 `BaseTool`，以 `@register_tool("PerturbationDatasetLoaderTool")` 注册，实现 `run(arguments)`。构造阶段不读取用户数据。

- category：`perturbation`，通过 `default_config.py` 映射配置。
- name：`perturbation.load_dataset`，用于运行时请求。
- type：`PerturbationDatasetLoaderTool`，与注册名一致。
- SDK callable：`perturbation_load_dataset`，保留点号形式的运行时调用名。

配置文件为 `ToolUniverse/src/tooluniverse/data/perturbation_tools.json`，顶层为列表。参数 Schema 拒绝未声明字段；配置声明 `requires_local_input: true`、`cacheable: false` 和 anndata 依赖。完整输入与返回 Schema 以该 JSON 为准。

## 参数

| 参数 | 类型 | 约束 |
| --- | --- | --- |
| `source` | string | 必填；本地文件或直接包含 10x 文件的目录 |
| `output_dir` | string | 必填；描述文件输出目录 |
| `dataset_id` | string | 可选；1 至 128 个字母、数字、下划线或连字符，以字母或数字开始 |
| `format` | string | 默认 auto；h5ad、10x_h5、10x_mtx、csv、tsv |
| `matrix_key` | string | 默认 X；H5AD 可选 raw/X、layers/名称 |
| `backed` | boolean | H5AD 默认 true；其他格式禁止显式传入 |
| `gene_id_column` | string | 基因标识列；省略时使用变量索引，必须唯一且非空 |
| `obs_column_map` | object | 标准字段到原始注释列名的映射，不修改源列名 |
| `control_values` | array[string] | 需显式映射 perturbation，且至少一个观测匹配 |
| `expression_scale` | string | counts、normalized、log1p、unknown；默认 unknown，不自动推断 |
| `matrix_orientation` | string | CSV/TSV 必填，observations_by_features 或 features_by_observations |
| `index_column` | string | CSV/TSV 必填，行 ID 所在列的确切表头名称 |
| `metadata_path` | string | 可选，外部 CSV/TSV/Parquet 注释表 |
| `metadata_join_keys` | object | 外部元数据必填，data 和 metadata 是等长非空列名数组 |
| `metadata_match` | string | strict 或 allow_missing，默认 strict |

相对路径以调用进程当前目录解析，描述中保存绝对源路径。CSV/TSV 的方向和索引列参数不得用于其他格式，非 X 矩阵选择只适用于 H5AD。

## 格式与数据校验

输出方向统一为观测 × 特征。10x 保持稀疏存储，要求特征数与条码数匹配且值有限、非负，不支持非 Gene Expression 混合特征。现代 H5 必须有 `/matrix` 布局；旧版多基因组 H5 及普通 HDF5 返回不支持错误。

MTX 目录必须各有一个 `matrix.mtx[.gz]`、`features.tsv[.gz]` 或 `genes.tsv[.gz]`、`barcodes.tsv[.gz]`。同一角色有多个候选时必须拒绝，不任意选取。特征文件为两列或三列，条码文件为一列。

CSV/TSV 支持 gzip，要求唯一非空的表头及观测、基因标识，除索引列外全部为有限数值。重复标识必须失败，不能静默追加后缀。输入不接受夹杂的实验注释列，注释应通过外部表关联。

H5AD 使用 backed 读取时，layers 仍可能被 AnnData 载入内存；选中 layer 时返回相应提示。不得为了摘要把完整稀疏表达转换为稠密矩阵。

## 元数据关联

`metadata_join_keys.data` 引用表达观测列，其中保留键 `__obs_index__` 表示观测索引。`metadata_join_keys.metadata` 引用外部列，支持复合键。保留键与原始列冲突必须失败。

关联按键进行并保留观测顺序，禁止按行号拼接。关联键不得为空，元数据键必须唯一，禁止一对多展开。默认每条数据记录必须匹配，元数据允许有未使用额外行；allow_missing 保留未匹配观测并在 warnings 返回数量。

外部数据列与原始观测列重名时必须报错，不自动覆盖。显式声明的映射列不存在时必须失败。关联规则保存在描述文件中；源数据及外部表保持原样，下游需按规则重建关联。

## 返回与持久化

成功返回 JSON 字典，含 status=success、data 和 warnings。data 包含 dataset_ref、dataset_id、modality=expression、format、shape、dtype、matrix_key、orientation、expression_scale 和 metadata_fields。metadata_fields 仅为字段名；提供 control_values 时另返回 n_controls。

不得返回完整矩阵、AnnData、NumPy 标量或打开的文件句柄。dataset_ref 为 JSON 描述文件路径，描述至少包含 schema_version=1、数据集标识、摘要、inputs、read_options 和 warnings。

inputs 为每个原始文件的绝对路径、字节数、mtime_ns、SHA256；MTX 的三个成员和外部注释均需登记。read_options 保存重建读取所需的矩阵选择、字段映射和关联规则。

加载前后检查输入指纹；验证成功后才创建唯一数据集目录，通过临时文件原子发布 manifest.json。已有数据集不静默覆盖。失败不留下成功引用，所有 backed 句柄均关闭。

`resolve_dataset_ref` 检查版本、各输入文件和 MTX 候选成员，返回描述，不返回打开的数据。源文件改变或移动后必须重新加载；临时进程内 ID 不能作为数据引用。

## 错误

工具检测到的错误包含 status=error、error_code、error。框架入参校验可在进入 run 前返回框架自身的错误结构。

| error_code | 条件 |
| --- | --- |
| INVALID_ARGUMENT | 参数或格式选项冲突 |
| SOURCE_NOT_FOUND | 输入不存在 |
| UNSUPPORTED_FORMAT | 格式或布局不支持 |
| DEPENDENCY_MISSING | 缺少读取依赖 |
| INVALID_DATA | 文件内容、维度或标识不合法 |
| METADATA_MISMATCH | 关联键、映射列或匹配结果不合法 |
| DATASET_EXISTS | 目标数据集已存在 |
| IO_ERROR | 底层读取或写入失败 |

空字典、None 或空矩阵不能作为成功结果。损坏文件的底层解码失败可能归入 INVALID_DATA 或 IO_ERROR。

## 验收与变更

通过类别加载和 `run_one_function(validate=True)` 验证，不只直接调用 Python 类。成功路径使用真实公开数据，转换副本明确标记；失败场景可修改临时副本，禁止覆盖源文件。

必须核对矩阵数值、方向与标识，覆盖 gzip、raw 和 layer、元数据打乱顺序及重复缺失键、源变化、重复输出、参数校验和跨进程引用。pytest 不自动下载；缺少真实数据时显式跳过集成测试。SDK 生成与调用名映射测试不依赖大文件。

大规模 backed 读取仍可能读取全部 obs，完整哈希会顺序扫描文件；不得宣称恒定内存或固定时延。新增格式必须补真实数据与验收用例，新增必填字段或不兼容语义变化必须同步版本、Schema、文档及测试。
