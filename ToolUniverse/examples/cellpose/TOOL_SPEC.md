# Cellpose_segment_image 接口规范

注册类型：`CellposeTool`；分类：`cellpose`；默认禁用 ToolUniverse 结果缓存。
JSON 参数契约以 `src/tooluniverse/data/cellpose_tools.json` 为准；运行时同时校验参数关联约束。

## 输入与参数

| 参数 | 类型 / 默认 | 语义与限制 |
| --- | --- | --- |
| `image_path` | string | 单个本地文件；与另外两种输入互斥 |
| `image_paths` | string[] | 1–512 个本地文件，保持请求顺序；需 output_dir |
| `image_ref` | string | loader 返回的 dataset_ref；需 output_dir |
| `image_indices` | int[] | 仅 image_ref 可用；去重、非负，宽表按 field_id 合并后的逻辑图索引 |
| `output_dir` | string | 启用完整输出；运行目录为 output_dir/run_id；应在输入图像目录/图像存储外 |
| `run_id` | string / 自动 UUID | 1–64 位字母数字、下划线、短横线，首位字母数字 |
| `resume` | bool / false | 续跑已有 run_id；需显式 run_id，缺失运行拒绝 |
| `device` | string / auto | auto 使用可用 CUDA:0，否则 CPU；cpu、cuda、cuda:N、mps 显式不可用时报错 |
| `pretrained_model` | string | 4.x 支持安装版本内置名称或可信的已有 checkpoint 路径；未知名称拒绝；3.x 不接受此选项 |
| `model_type` | string / cyto3 | cyto3、cyto2、cyto、nuclei；仅 3.x 生效；4.x 显式填写时报错，省略则使用安装版本默认权重 |
| `diameter` | number / 后端默认 | 非负有限数；0 同未指定。3.x 可估计尺寸，4.x 默认缩放，不自动估计尺寸 |
| `channels` | int[2] / [0,0] | 旧 3.x 细胞/细胞核通道码，取值 0–3；4.x 显式填写时报错 |
| `channel_names` | string[] | 1–3 个唯一的 loader 声明通道名，按请求顺序输入模型；与 channel_indices 互斥 |
| `channel_indices` | int[] | 1–3 个唯一的非负物理通道索引；宽表按 image_path_columns 顺序；所有后端均生效 |
| `z_index`, `t_index` | int / 单例为 0 | 多平面轴必须明确选择；不存在的轴只接受 0；不做隐式投影 |
| `roi` | object | y、x、height、width 均必填整数；起点非负，尺寸正数，不能超出平面 |
| `max_pixels` | int / 16777216 | 选择后的像素样本数，包含通道，最大 67108864 |
| `batch_size` | int / 8 | 单图 tile 微批大小，1–64；外层图像串行处理 |
| `flow_threshold` | number / 0.4 | 非负有限数，流误差阈值 |
| `cellprob_threshold` | number / 0.0 | 有限数，细胞概率阈值 |
| `mask_format` | string / tiff | tiff、npy 均 uint32；png 为 uint16，标签超过 65535 报错 |
| `save_input` | bool / true | 保存模型调用前选择后的原始像素 input.npy |
| `save_preview` | bool / true | 保存展示用叠加图，最长边最多 1024，独立对比度拉伸 |
| `save_flows` | bool / false | 保存后端返回的全部原生 flow 数组 |
| `save_styles` | bool / false | 保存后端原生 styles；某些版本返回占位数组，不保证有用的 embedding |
| `preview_rows` | int / 20 | 0–100，同时限制响应内图像记录数和对象预览数；文件保留完整记录 |
| `on_error` | string / continue | continue 继续其他图；stop 将剩余图标为 not_processed |
| `max_run_seconds` | number / 不限 | 正有限秒数，仅完整模式；在安全边界检查，不硬中断当前操作 |
| `cancel_file` | string | 仅完整模式；调用方创建此本地标记后在安全边界停止，工具不会删除标记 |
| `save_mask` | bool / false | 仅旧单图模式：保存 uint16 PNG |
| `mask_output_path` | string | 仅旧单图模式：默认去掉 image_path 扩展名后加 _cp_masks.png，禁止覆盖源文件 |

`image_path`、`image_paths`、`image_ref` 必须且只能提供一种。
省略通道选择时，源通道超过 3 个会报错；不能静默丢弃通道。
`image_indices` 仅适用于 image_ref。当前最多一次处理 512 张逻辑图；引用索引读取上限
100000 条记录、128 MiB。大引用需先缩小数据范围或选择索引。

以下选项需完整输出模式：批量/引用输入、run_id、resume、mask_format、save_input、
save_preview、save_flows、save_styles、preview_rows、on_error、max_run_seconds、cancel_file。
设置 output_dir 后必须省略 save_mask 和 mask_output_path（即使 save_mask=false）。

## 图像读取

本地路径使用现有 loader 的头信息解析。image_ref 沿用 loader 已选的 series、
resolution level、NGFF group、stage position，不覆盖这些选择。

| 格式 | 读取策略 |
| --- | --- |
| TIFF / OME-TIFF | tifffile + Zarr 选取平面/ROI；无 Zarr 时仅允许整体解码不超过 256 MiB 的 acquisition |
| PNG / JPEG / BMP | PIL 解码，再选通道与 ROI |
| OME-NGFF Zarr | 本地或 loader 支持的匿名 HTTPS/S3，按选中 array/level 切片；结束关闭连接 |
| LIF | 按指定 T/Z/channel 读取平面 |
| CZI / ND2 | 解码 loader 选择的 acquisition 后切片；整体大小上限 256 MiB |

CZI、ND2、LIF 和 NGFF 需对应可选依赖。ND2 在 GPFS 不支持 mmap 时临时暂存本地，
退出读取后清理。远端超时、endpoint 等连接参数沿用 loader 引用；不新增认证接口。
输出为 YX 或 YXC，拒绝非有限数、非数值像素、不支持的轴组合、越界选择。
这不是三维模型或时间序列自动展开接口。

## 返回契约

无 output_dir 的旧单图成功返回：

```json
{"status":"success","data":{"image_path":"...","model_type_requested":"cyto3","model_used":"cpsam_v2","device":"cpu","num_objects":15,"image_shape":[256,256],"objects":[{"label":1,"area":30,"centroid_y":0.83,"centroid_x":185.03}],"mask_path":null}}
```

objects 为完整对象列表；4.x 另有 note。旧 PNG 保存失败时可能出现 mask_save_error，
推理统计仍可成功返回；完整输出模式导出失败则该项标为 error。

设置 output_dir 后，单图和批量统一返回：

```json
{"status":"success","data":{"run_ref":".../run.json","results_path":".../results.jsonl","run_id":"experiment_001","n_requested":2,"n_processed":2,"n_succeeded":2,"n_failed":0,"n_unprocessed":0,"resumed_successes":0,"device":"cpu","results":[]}}
```

示例 results=[] 对应 preview_rows=0；默认返回最多 20 条。
每条成功记录含 item_index、paths、field_id、status、prediction_ref、mask_path、
objects_path、num_objects、image_shape、objects_preview、artifacts、model、device、
input_pixel_sha256、result_stamp 及耗时。重新生成的成功项含 attempt_dir；跳过项含 resumed=true。
错误记录含 item_index、paths、field_id、status=error、error、耗时；失败时可能保留 attempt_dir。
未处理记录 status=not_processed，包含停止原因。每次调用的计数只计算本轮已处理项。
未处理项或临时失败项若曾成功，完整索引记录的 previous_success 保存上次成功记录，
供下次重新校验后复用；响应预览中的 previous_success 只展示简短引用。
主动停止时 data.stop_reason 为 cancel_requested、time_budget_exhausted 或 image_error。

整体 status：全部成功为 success；存在成功和失败/未处理为 partial_success；无成功为 error。
参数、设备、引用开始前全局校验失败返回 status=error、error，可能没有 data；运行中可捕获的全局异常返回 error 及运行路径 data，已写文件保留。
应先检查 status，再访问 data。分割零实例属于成功，mask 全背景，对象表为空。

## 文件契约

```text
output_dir/
  .<run_id>.lock
  <run_id>/
    run.json
    results.jsonl
    items/00000/<attempt_uuid>/
      result.json
      mask.tif             # 或 mask.npy、mask.png
      objects.jsonl
      input.npy            # 默认保存，可关闭
      preview.png          # 默认保存，可关闭
      flow_*.npy           # 可选，实际命名见 artifacts.flows
      cell_probability.npy # 后端有此 flow 时
      native_styles.npy    # 可选
      native_diameters.npy # 后端有第四个返回值时
```

- run.json：schema_version=1、请求、状态、处理/成功/失败计数、续跑跳过计数、结果索引。
  每张图提交后更新；运行中 status=running；含 phase、current_item_index、elapsed_seconds、
  n_unprocessed；可通过普通文件读取或 stream_callback 观察。不是单独的查询工具。
- results.jsonl：从运行开始即包含全部图像位置的最新原子快照；未处理位置保留 previous_success，避免续跑中断丢失历史成果。
- result.json：schema_version=1、source、model、num_objects、对象预览、坐标说明、artifacts。
- mask：二维实例标签，0 为背景；TIFF/NPY 为 uint32，PNG 为 uint16。
- objects.jsonl：每个实例 label、area、centroid_y、centroid_x；面积为像素数，质心单位为像素。
  坐标相对 ROI；加 source.origin_yx 转回来源平面。
- input.npy：确切选择的像素，包括通道顺序；不含模型内部归一化或补齐操作。
- flows：逐个保存原生数组，标记 native_index、name、shape、dtype。常见索引 0 为颜色流图，
  1 为向量场，2 为 cell probability，3 为后端返回的位置数组；实际数量由后端决定。
- styles：直接保存原生数值数组，不能保证其生物学含义或可用于下游表征任务。
- artifacts：文件路径、字节数、mtime_ns、SHA256；数值文件附带 shape/dtype。
  result_stamp 为 result.json 指纹；input_pixel_sha256 同时覆盖选中像素的 dtype、shape 与内容。

## 续跑、并发与执行边界

run_id 已存在时需 resume=true；不存在时不能 resume。
续跑按规范化后的有效参数比较；已知数值/保存默认值、直达路径和单例 Z/T=0 统一处理。
允许调整 preview_rows、on_error、max_run_seconds、cancel_file；output_dir/run_id 用于定位运行。
解析后的设备与 image_ref manifest 指纹也属于匹配条件。仍建议重用原参数字典；通道选择、
图像索引等依赖具体来源的选项不能自动推断等价性。旧 schema_version=1 快照仍可读取。

跳过成功项之前重新读取输入像素，并验证模型版本/权重/推理选项、来源元数据、
result.json 与全部导出文件指纹。失败项、像素变化项、损坏项重新推理，每次写入独立 attempt 目录。
输入引用在开始和结束均校验；如果 loader 自身发现本地源/索引变化，则需重新 load。
远端引用的整体校验范围沿用 loader 元数据快照；选中像素哈希及 input.npy 记录本次实际输入，
不保证整个远端数据集不可变。
4.x 可通过已安装的默认/指定权重路径及完整指纹识别模型，不构造网络即可跳过成功项。
未知默认、缺失权重或旧后端退回正常模型加载验证；不以牺牲权重校验换取跳过速度。

相同 run_id 使用非阻塞文件锁，拒绝并发写入。模型按设备/选择缓存在进程内，模型初始化及
推理使用锁；外层图像串行执行。run.json/results.jsonl 原子更新，但不提供完整历史日志。
中断后可复用已提交的成功记录；当前未提交项会重算。进程被强杀可能保留 running 状态，
不妨碍锁释放后的续跑。

工具的 run 是同步计算，可用 ToolUniverse.run_one_function_async 在线程中执行并异步等待。
没有独立后台 job/poll 服务。stream_callback 接收 JSON 进度事件（event=cellpose_progress），
包含运行路径、阶段与计数；观察者普通异常记录为 progress_warning，计算继续。
max_run_seconds 是本次调用的协作预算，cancel_file 是调用方控制的停止标记。
二者在每张图之前及开始推理之前检查；读图、模型加载、推理本身不可被硬中断。
因此超时可能在当前操作结束后生效；调用方等待超时/取消等待不表示计算已停止。
建议 shield 同一个异步任务、定期检查 run.json，取消时创建标记并等待安全退出。
完整模式若无成功项而主动停止会返回 error + data.stop_reason；有成功项则 partial_success。
调用方和远端存储可有各自超时。大图需控制 ROI、max_pixels 和 tile batch_size，
这些限制不保证 GPU 显存足够。训练、体分割和自动展开 T/Z 均不在当前接口范围。
