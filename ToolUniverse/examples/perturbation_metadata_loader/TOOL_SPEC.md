# 扰动实验元数据工具规范

接口版本 0.1。专用工具解释扰动实验字段，通用 metadata_io 只负责读取表格；不增加按药物、基因、孔板拆分的独立工具，也不实现通用 data.load_metadata 入口。

## 接入与职责

category 为 perturbation，name 为 perturbation_load_metadata，别名为 perturbation.load_metadata，type 为 PerturbationMetadataTool。实现继承 BaseTool，构造阶段不读取用户文件。配置追加在 data/perturbation_tools.json，不覆盖已有表达 loader。SDK 与静态注册表同时生成。

输入 Schema 拒绝未声明字段，配置 requires_local_input=true、cacheable=false。参数包括 source、output_dir、experiment_level、field_map、association_keys 及可选格式、sheet、常量、对照、缺失与去重规则。完整 Schema 以配置为准。

本工具只标准化和验证实验注释，不连接远端数据库，不推断机制，不修改 source 或已有数据集，不以加载成功代替完整实验设计验证。

## 数据规则

必须支持 CSV/TSV gzip、XLSX、Parquet 和 H5AD obs，按格式延迟导入依赖。H5AD 只读取 obs，不实例化完整 AnnData。列名必须唯一非空，输入不得为空。表格字符串标识保留，不根据内容自动识别生物学含义。

标准字段为 perturbation、perturbation_type、cell_id、sample_id、plate_id、well_id、site_id、cell_type、dose、dose_unit、time、time_unit、batch、replicate、is_control。field_map 从标准字段指向原始列；constants 显式提供标准字段常量，两者不得重叠。原始未映射列保留在源，不进入标准化输出。

cell 必须有 cell_id；sample 有 sample_id；well 有 plate_id 和 well_id；field 再有 site_id。相应 ID 必须属于 association_keys，调用方可增加 batch、time 等键。关联键不得缺失，输出中必须唯一。不将本工具绑定到某个 dataset_ref；下游按键关联。

perturbation 必须非空。默认错误；只有显式 missing_perturbation=drop_rows 才丢弃缺失标签记录并报告数量，不能填充未知处理或将缺失视为对照。全部记录被丢弃时失败。

可选字段的空值按 missing_optional 处理。dose/time 必须有限非负，不把非数值文本当成合法缺失。单位不推断、不自动换算；有值而无单位时警告。非标量列不得转换成伪标签字符串。空白标识不能作为关联键。

is_control 只来自显式字段、常量或 control_values。后者精确匹配扰动标签，并必须至少命中一条；同时声明时必须一致。未声明时不生成全 false 的控制列，返回 n_controls=null 和 CONTROL_UNSPECIFIED。

默认重复键失败。collapse_identical 仅折叠所有标准化字段完全一致的重复记录，冲突始终失败。原始不同但未映射的列不参与比较。丢弃未标注记录数和折叠数分别报告；可选缺失及单位警告以实际持久化行数计数。

## 输出与完整性

响应必须为 JSON，成功包含 status=success、data、warnings。data 提供 metadata_ref、table_path、metadata_id、experiment_level、n_source_rows、n_rows、fields、association_keys、n_perturbations、n_controls 和至多 20 条 preview。preview 来自实际标准化表，不是生成的示例。

JSONL 标准表使用 null 表示缺失，使用数值表示 dose/time，使用布尔表示已声明的 is_control；其余标准字段作为字符串。不得将完整表或 AnnData 放入响应。默认 preview 为 5，设 0 时返回空列表。

描述文件 schema_version=1、modality=perturbation_metadata，保存源与标准表的绝对路径、大小、修改时间与 SHA256，以及 read_options、summary、warnings。source 在处理前后校验，同名输出不覆盖，全部验证成功后原子发布 manifest。失败清理未发布产物。

resolve_perturbation_metadata_ref 检查源与标准表。文件移动或改变后引用失效，不能将短期进程内对象作为持久化引用。下游可按 descriptor.table.path 打开标准表，按 summary.association_keys 关联。

## 错误与验收

工具错误统一返回 status=error、error_code、error。使用 INVALID_ARGUMENT、SOURCE_NOT_FOUND、UNSUPPORTED_FORMAT、DEPENDENCY_MISSING、INVALID_DATA、METADATA_MISMATCH、DATASET_EXISTS、IO_ERROR。框架参数校验可在 run 前返回框架错误结构。

真实 BBBC021 验证视野到孔的显式一致记录折叠、对照数、剂量及重复冲突；真实 Adamson obs 验证细胞级标注、缺失扰动标签和不加载表达矩阵。不同表格格式用真实注释转换副本测试，不生成合成观测。

必须验证类别加载、别名、SDK、控制声明冲突、缺失、非有限或负剂量、重复键、输出预览、已有目录、源与标准表变化、跨进程引用。测试离线读取本地数据，不在 pytest 中自动下载；测试目录与数据不随工具提交。
