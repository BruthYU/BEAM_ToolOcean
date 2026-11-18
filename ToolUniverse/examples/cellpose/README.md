# Cellpose 图像分割工具

`Cellpose_segment_image` 在本地运行 Cellpose。支持单图、路径列表，以及
`image_perturbation_load_images` 返回的 `dataset_ref`。推荐设置 `output_dir`，
获得完整、可追溯的结果文件。参数与返回契约见 [TOOL_SPEC.md](TOOL_SPEC.md)。

## 环境

本机已配置独立环境和官方 CPSAMv2 权重：

```bash
source /root/.tooluniverse/cellpose/activate.sh
```

本机命令行入口也支持完整参数 JSON：

```bash
python /root/.tooluniverse/cellpose/call.py --arguments-json /absolute/path/arguments.json
```

JSON 文件内容就是下方示例中的 `arguments` 对象，不包含工具名称或外层封装。
脚本与 activate.sh 是本机配置文件，不随仓库分发。

其他机器可通过 `python -m pip install -e './ToolUniverse[cellpose]'` 安装已验证的 4.x 后端及基本输出依赖；其他图像格式仍需对应可选依赖。
TIFF 使用 tifffile；大 TIFF 切片和 NGFF 使用 Zarr；CZI、ND2、LIF 分别需要
对应读取库。权重使用 Cellpose 的官方下载机制或指定可信的本地 checkpoint。
首次调用可能包含权重下载、加载开销。

## 单图与批量调用

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse(load_workspace=False)
try:
    tu.load_tools(categories=["cellpose"], clear_existing=True)
    arguments = {
        "image_paths": ["/absolute/path/cells_a.tif", "/absolute/path/cells_b.png"],
        "output_dir": "/absolute/path/cellpose_results",
        "run_id": "experiment_001",
        "device": "cpu",
        "save_flows": True,
        "save_styles": True,
        "on_error": "continue",
    }
    result = tu.run_one_function(
        {"name": "Cellpose_segment_image", "arguments": arguments},
        validate=True, use_cache=False,
    )
    print(result)
finally:
    tu.close()
```

单图完整输出把 `image_paths` 替换为 `image_path`。三种输入只能选择一种。
默认生成 uint32 TIFF mask、完整对象 JSONL、选中输入像素 NPY 和分割预览。
`save_flows` / `save_styles` 默认关闭，开启时保存后端的原生数值数组。

已有运行不能直接覆盖。续跑使用原参数，并增加 `resume=True`；必须提供同一个
`run_id`。成功项会校验选中像素、模型、来源信息及结果文件后跳过；失败项和损坏项重算。
`preview_rows`、`on_error`、`max_run_seconds`、`cancel_file` 可以调整。
默认参数会统一比较，例如首次省略 batch_size、续跑填写 8 等价。
未处理项的历史成功结果保留在索引的 `previous_success` 中；即使续跑再次中断也能恢复。
4.x 可从已安装权重校验成功项，不必先实例化模型；仍需读取输入及校验文件哈希。

也可以使用 SDK：

```python
from tooluniverse.tools import Cellpose_segment_image

result = Cellpose_segment_image(
    image_path="/absolute/path/cells.tif",
    output_dir="/absolute/path/cellpose_results",
    device="cpu",
    use_cache=False,
)
```

## 接入图像 loader

取得 loader 返回的 `dataset_ref` 后，将其作为 `image_ref`：

```python
arguments = {
    "image_ref": loaded["data"]["dataset_ref"],
    "image_indices": [0],
    "channel_names": ["Actin", "DNA"],
    "roi": {"y": 256, "x": 256, "height": 256, "width": 256},
    "output_dir": "/absolute/path/cellpose_results",
    "device": "cpu",
}
```

通道名称必须与 loader 的记录一致；也可用 `channel_indices` 按物理通道顺序选择。
宽表中相同 `field_id` 的通道合并为一张逻辑图，`image_indices` 在合并后编号。
Z/T 轴有多个平面时必须显式提供 `z_index` / `t_index`。此接口执行二维分割。

## 读取结果

```python
import json
import numpy as np
import tifffile

assert result["status"] in ("success", "partial_success")
with open(result["data"]["results_path"]) as f:
    records = [json.loads(line) for line in f]
for record in records:
    if record["status"] != "success":
        print(record)
        continue
    mask = tifffile.imread(record["mask_path"])  # 默认 TIFF；NPY 用 np.load
    with open(record["prediction_ref"]) as f:
        prediction = json.load(f)
    with open(record["objects_path"]) as f:
        objects = [json.loads(line) for line in f]
    if "input" in record["artifacts"]:
        pixels = np.load(record["artifacts"]["input"]["path"], allow_pickle=False)
```

mask 为二维实例标签，0 表示背景。对象面积单位是像素，质心坐标以选中 ROI 为原点；
加上 `prediction["source"]["origin_yx"]` 即得到来源平面的坐标。
`input.npy` 保留选择、排序、裁剪后的原始像素；模型仍会执行其默认内部归一化。
预览图只用于展示，不能替代原始输入或 mask。

不设置 `output_dir` 的旧单图调用继续返回完整对象统计，默认不保存 mask。
旧接口的 `save_mask`、`mask_output_path` 不能与完整输出模式同时使用。
Cellpose 4.x 下显式填写旧 `model_type`、`channels` 会报错；应省略它们，使用 `pretrained_model`
选择权重，用 `channel_names` / `channel_indices` 选择实际输入通道。

## 已验证范围

本机使用 Cellpose 4.2.1.1、官方 CPSAMv2、CPU，在真实 BBBC021 图像上完成：
单图及批量推理、完整文件导出、缺失输入的逐项失败、续跑跳过、损坏输出重算、
补齐缺失文件后重试、多通道 loader 引用、ROI 原始像素逐元素核对、SDK 调用和旧接口兼容。
256×256 DAPI 裁剪得到 15 个实例；这是执行与输出一致性检查，不是分割精度评估。
续跑耗时包含像素及文件指纹校验，首次导入依赖也需要时间；推理耗时随图像尺寸、设备和首次模型加载而变化。

相关验证覆盖参数、真实图像读取、完整结果及续跑恢复；新增稳定性测试覆盖重复中断、临时读取失败、默认参数等价、取消、时间预算和观察者断连。读取层另外验证了真实 OME-TIFF、CZI、ND2、LIF、
NGFF 转换数据及公共 HTTPS NGFF 的选中像素；不表示所有格式均做过模型精度评估。
GPU/MPS 本机无设备，未做硬件推理验证。3.x 兼容分支保留，但本机未完成 3.x
运行验证（该版本 NumPy 约束与当前 ToolUniverse 环境冲突）。

## 长任务进度、取消与异步等待

同步调用可传 `stream_callback` 接收 JSON 进度消息。开始、读取、模型加载、推理、
单图完成和最终状态都会通知，并附带 run_ref/results_path。观察者抛出普通异常会记录
progress_warning，不会中止或重复推理。run.json 同时记录 phase、current_item_index、
elapsed_seconds 和当前调用的计数；不需要等待整批完成即可检查进度。

```python
import asyncio
import json
from pathlib import Path

async def run_batch(tu, arguments):
    arguments = {
        **arguments,
        "max_run_seconds": 3600,
        "cancel_file": "/absolute/path/cancel_experiment_001",
    }
    # 创建一个 task，继续等待同一 task；不要因等待超时重新提交同一个 run_id。
    task = asyncio.create_task(tu.run_one_function_async(
        {"name": "Cellpose_segment_image", "arguments": arguments},
        stream_callback=lambda text: print(json.loads(text)["phase"]),
        use_cache=False, validate=True,
    ))
    while not task.done():
        try:
            return await asyncio.wait_for(asyncio.shield(task), timeout=30)
        except asyncio.TimeoutError:
            print("仍在运行，可检查 run.json；需要取消时创建 cancel_file。")
    return await task

# 在另一进程/终端创建标记请求停止：
# Path("/absolute/path/cancel_experiment_001").touch()
```

时间预算与取消标记在图像之间及开始推理前检查，**不能硬中断正在执行的读图、加载或推理**。
当前图会安全提交，后续项保留可续跑状态；因此实际时间可能超过 max_run_seconds。
返回 data.stop_reason 区分 cancel_requested、time_budget_exhausted、image_error。
取消后重新运行需先移除标记（工具不会自动删除），再以 resume=true 续跑。
只取消 Python 的等待 task 也不会杀死正在运行的计算线程；应使用取消标记并等它安全退出。

外层图像依次执行并复用模型；`batch_size` 指每张图的模型 tile 微批大小。
没有独立后台任务服务、训练或三维体分割接口。测试数据、结果和权重保留本地。
