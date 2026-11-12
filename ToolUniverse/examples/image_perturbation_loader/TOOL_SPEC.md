# 图像数据加载工具构造规范

接口版本为 0.1。本工具登记 TIFF、OME-TIFF、BigTIFF、PNG、JPEG、BMP 、CZI、ND2、LIF 与本地／匿名远端 OME-Zarr 0.4／0.5，输出索引和数据描述。表达数据、表型特征、分割及预测由其他工具处理。

## 框架接入

category 为 image_perturbation；规范 name 为 image_perturbation_load_images；兼容别名为 image_perturbation.load_images；注册 type 为 PerturbationImageDatasetTool。继承 BaseTool，实现 run(arguments)，构造阶段不读取用户数据。

配置位于 src/tooluniverse/data/image_perturbation_tools.json，由 default_config.py 映射；配置拒绝额外参数，声明 requires_local_input=true、cacheable=false。依赖按格式导入，TIFF 需要 tifffile，PNG／JPEG／BMP 需要 Pillow，OME-Zarr 需要 zarr，远端还需 fsspec 与 aiohttp 或 s3fs；CZI／ND2／LIF 分别需要 czifile>=2026.11.8、nd2 和 readlif。CZI 后端要求 Python>=3.12。SDK 须随配置同步生成；实现类通过 @register_tool 注册，并由源码 AST 扫描自动发现，本 PR 不修改静态注册表。

## 参数约束

| 参数 | 类型 | 规则 |
| --- | --- | --- |
| source | string | 必填，本地文件、目录或清单 |
| source_type | string | 必填，file、directory、manifest、remote |
| output_dir | string | 必填，索引与引用输出位置 |
| dataset_id | string | 可选，唯一标识，禁止路径分隔符及目录跳转 |
| format | string | auto、tiff、ome_tiff、ome_zarr、png、jpeg、bmp、czi、nd2、lif，默认 auto |
| series_index | integer | 零起始 series，多个 series 时必须显式提供 |
| resolution_level | integer | 零起始层级，默认 0，适用于 TIFF 和 OME-Zarr |
| group_path | string | OME-Zarr 内部图像组，默认根组；不得跳出存储 |
| position_index | integer | LIF mosaic 位置，零起始；多个位置时必须显式提供 |
| remote_timeout | number | 仅 remote，每次请求超时，默认 30 秒，最大 300 秒 |
| remote_endpoint | string | 可选匿名 S3-compatible endpoint，仅适用于 s3 URL |
| file_pattern | string | 仅目录模式，默认第一层支持格式；递归须显式声明 |
| image_path_columns | object | 清单模式必填，通道名到路径列的映射 |
| image_root | string | 仅清单模式，可覆盖相对图像路径的基准目录 |
| metadata_column_map | object | 仅清单模式，标准字段到原始清单列 |
| field_keys | array[string] | 清单视野唯一键，使用映射后字段；省略用行号 |
| missing_images | string | error 或 skip_field，默认 error |
| axes | string | 可选，使用 TCZYX 的非重复轴，须含 YX 且长度等于 ndim |
| channel_map | object | 内部通道索引到名称，须完整且唯一；不与清单通道混用 |
| lazy | boolean | 默认 true，false 逐文件解码检查 |

相对 source、output_dir 和 image_root 以调用进程目录解析；清单内相对路径以 image_root 或清单目录解析。清单只接受 CSV/TSV，列名须唯一非空，映射列须存在，元数据值以原始文本保存。

## 图像与实验分组

TIFF 与 OME-TIFF 支持显式选择 series 和分辨率层级，多 series 不得静默选择第一幅。BigTIFF 由容器标记识别。多文件 OME-TIFF 必须在拼装前验证全部 UUID 和 FileName，登记每个成员并关闭全部句柄；不能漏登记关联文件或用零填补缺失像素。目录模式合并同一文件集入口，清单重复声明时报错。PNG／JPEG／BMP 根据实际文件签名和模式检查。无法确定轴含义时返回 INVALID_ARGUMENT，禁止自动解释任意轴。OME 的通道名称、轴和尺度元数据须保留；禁止不一致的轴覆盖。OME-TIFF 物理像素尺寸须按所选金字塔层级换算。

清单每行是一个逻辑视野，宽表各通道必须是单通道文件，并具有一致的非通道维度。声明 field_keys 时值非空且组合唯一；同一图像路径重复登记必须失败。目录或单文件模式不推断实验分组，n_fields 为 null。

缺失任一通道默认失败；显式 skip_field 时跳过整个视野并记录数量。不得静默丢弃缺失图像、任意选择重复记录或保留残缺通道组。所有图像缺失时不能返回成功。

## 持久化与返回

成功响应为 status=success、data、warnings。data 包含 dataset_ref、index_path、dataset_id、modality=image、n_files、n_fields、channels、metadata_fields、pixel_validation、n_input_files、integrity_scope。n_files 表示图像索引条目数；多文件 OME 图像可对应多个物理成员。像素检查程度为 header 或 decoded，不能混淆两者。

索引为 JSONL，每行含图像绝对路径、field_id、shape、axes、dtype、channels、spatial_shape 和 metadata。完整索引不放进响应，避免大清单占用模型上下文。

描述 schema_version=1，保存读取参数、原图和清单指纹、索引指纹及摘要。OME-Zarr 保存聚合目录指纹，覆盖所有文件名、内容和修改时间；增删块或元数据均需重新加载，输出不得位于输入存储内。发布 manifest 前必须完成输入校验；失败清理未发布产物；已有 dataset_id 不覆盖。源文件保持只读，所有图像句柄在成功和失败路径关闭。

resolve_image_dataset_ref 验证引用、源与索引完整性，目录模式额外检查匹配成员。返回描述，不持有解码数组或图像句柄。引用依赖原始数据位置，源移动或改变后重新登记。

## 错误与性能

工具错误统一返回 status=error、error_code、error。代码包括 INVALID_ARGUMENT、SOURCE_NOT_FOUND、UNSUPPORTED_FORMAT、DEPENDENCY_MISSING、INVALID_DATA、METADATA_MISMATCH、DATASET_EXISTS、IO_ERROR。框架输入 Schema 校验可在 run 之前返回框架错误。

lazy 模式不解码 TIFF 像素，但文件仍有两次完整 SHA256 读取；解析引用再次校验。清单和索引构造占用与记录数量相关的内存，TIFF 解码模式还需容纳所选 series／层数组，OME-Zarr 解码按块进行，不承诺恒定内存或固定耗时。

## 本地验收

使用真实 BBBC021 三通道处理组与对照组，以及真实 CellProfiler PNG。核对索引维度、dtype、实际通道路径和处理标签，不只检查 status。覆盖文件与目录、CSV/TSV、内部多通道、lazy 不解码、显式解码、缺失整组跳过、重复键、损坏输入、源及索引变化、跨进程引用和 SDK 路由。

由真实像素重新封装的多通道或 OME 文件仅作格式测试，不标为新的原始实验。测试目录、源图像和运行产物不提交。新增格式或参数须同步 Schema、SDK、文档和真实输入验证。


## OME Zarr 元数据约束

当前仅接受 NGFF 0.4／0.5 图像组。0.4 从组属性读取 multiscales，0.5 从 ome 命名空间读取。plate／well 容器必须通过 group_path 明确选择其图像组，不自动遍历。多个 multiscales 通过 series_index 选择。

每层必须有唯一 scale 变换，可附加 translation，并按序组合数据集级和图像级变换；索引同时保留原始变换与合成 scale／translation。轴数量须匹配数组 ndim，支持 TCZYX 及原始单位，通道数须匹配 omero.channels。路径、尺度和层级越界必须报错。

图像索引的 image_metadata 记录各后端可用的通道、物理尺度、series、层级及 NGFF 元数据。不存在的信息保持缺失，不伪造单位或实验标识。任意仿射变换及其他 NGFF 版本尚未接入；远端和仪器文件按下列约束处理。


## 仪器与远端约束

CZI 场景、ND2 位置、LIF 内部图像统一由 series_index 选择。LIF mosaic 由 position_index 单独选择。不得把无法表达的特殊维度静默丢弃；必须返回不支持错误。CZI 保留后端通道与带单位坐标；ND2 保留轴和已知物理尺寸；LIF 将 pixels-per-micrometer 校准正确取倒数，位深仅支持统一 8/16 位存储。

ND2 仅在 mmap 明确不受文件系统支持时回退到本地临时副本，记录 staged_local_copy，任何成功或失败路径均清理临时文件并关闭 reader，输入指纹仍针对原始文件。

remote 只接受匿名 http／https／s3 OME-Zarr URL，显式限定元数据版本与图像组选取。remote_endpoint 只用于匿名 S3-compatible 服务。拒绝 URL 内嵌凭据和查询参数，不把鉴权对象或客户端配置写入引用。

本地 integrity_scope=full_content；远端 integrity_scope=metadata_only，并返回 REMOTE_METADATA_ONLY 警告。远端引用校验所选组和数组的元数据内容摘要，不能声称发现仅像素块变化。lazy=false 只是解码检查，不将远端内容升级为完整指纹。每个远端元数据对象最大 4 MiB，请求超时可配置，整个工具无固定超时。

索引、Schema、SDK 和文档须同时体现这些边界。真实 native 文件、本地 HTTP 服务承载的真实像素以及公开远端元数据用于验证；不提交测试目录、样例数据或运行产物。
