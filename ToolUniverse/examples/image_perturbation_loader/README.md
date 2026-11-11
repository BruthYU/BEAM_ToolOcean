# 图像扰动数据加载工具

`image_perturbation_load_images` 将本地显微图像登记为可追踪的数据集，返回摘要、JSONL 图像索引和持久化引用。兼容别名为 `image_perturbation.load_images`，类别为 `image_perturbation`。

## 支持范围

支持普通 TIFF、OME-TIFF、BigTIFF、PNG、JPEG、BMP、本地及匿名远端 OME-Zarr 0.4／0.5、图像目录及 CSV/TSV 宽表清单。清单每行代表一个视野，各通道的路径分别存放在不同列。单文件也可包含多个通道，须有明确的轴定义。

当前支持多文件 OME-TIFF 和 CZI／ND2／LIF 常见成像布局；不自动展开整板 OME-Zarr，也不支持其他 NGFF 版本、需要登录的远端存储、分割掩膜专用语义或独立元数据表关联。它不做光照校正、分割、形态特征提取或预测。清单中的实验注释保留在索引中。

## 安装与基本调用

在 ToolUniverse 项目目录安装源码和读取依赖：

```bash
python -m pip install -e . tifffile pillow zarr fsspec aiohttp nd2 readlif "czifile>=2026.11.8" imagecodecs
```

替换为实际图像目录后调用：

```python
import json
from tooluniverse import ToolUniverse

engine = ToolUniverse(load_workspace=False)
try:
    engine.load_tools(categories=["image_perturbation"], clear_existing=True)
    result = engine.run_one_function({
        "name": "image_perturbation_load_images",
        "arguments": {
            "source": "/data/images",
            "source_type": "directory",
            "output_dir": "/tmp/beam-image-datasets",
            "file_pattern": "**/*.tif",
            "lazy": True
        }
    }, validate=True, use_cache=False)
    print(json.dumps(result, ensure_ascii=False, indent=2))
finally:
    engine.close()
```

最少参数为 source、source_type 和 output_dir。source_type 可为 file、directory、manifest 或 remote。CZI 后端使用 czifile 2026.11.8 及之后的场景 API，需要 Python 3.12 或更新版本；其他依赖按所选格式导入，未安装时会返回依赖缺失错误。S3 访问另需安装 s3fs。SDK 入口为 `from tooluniverse.tools import image_perturbation_load_images`。

目录模式默认选择第一层支持格式的图像和 OME-Zarr 存储，不递归。file_pattern 可显式指定 glob；非图像 sidecar 文件会忽略，Zarr 内部元数据与块不会被当成独立图像。目录和单文件模式不根据文件名猜测视野、通道生物学含义或处理组，n_fields 返回 null，并给出分组未知提示。

## 清单与通道

清单可包含 dna_path、actin_path、plate、well、site、compound 等列，其中图像路径可为绝对路径或相对路径。默认相对清单目录解析，image_root 可覆盖路径基准。必须使用有效本地路径，历史服务器路径不会被猜测替换。

调用参数示例：

```json
{
  "source": "/data/image_manifest.csv",
  "source_type": "manifest",
  "output_dir": "/tmp/beam-image-datasets",
  "image_root": "/data/images",
  "image_path_columns": {
    "DNA": "dna_path",
    "Actin": "actin_path"
  },
  "metadata_column_map": {
    "plate_id": "plate",
    "well_id": "well",
    "site_id": "site",
    "perturbation": "compound"
  },
  "field_keys": ["plate_id", "well_id", "site_id"]
}
```

field_keys 指定映射后的视野唯一键，要求值非空且组合唯一；省略时以清单行号标识视野。同一孔多个视野必须包含 site_id 等区分信息。多通道文件不能直接作为宽表中的单通道文件；清单通道与单文件内部 channel_map 不混用。

缺失图像默认报错。明确设置 missing_images=skip_field 时，整个不完整视野被跳过，并在 warnings 报告数量，不保留残缺通道组。若没有完整视野则失败。真实实验清单中含未下载图像时，应选择完整本地视野或显式接受跳过。

## 轴与像素读取

支持明确的 T、C、Z、Y、X 轴，必须包含 Y 和 X。普通 TIFF 自带轴不明确时必须传 axes，不能猜测 Q 或 I 轴含义。OME-TIFF 和 OME-Zarr 使用权威元数据，禁止用不一致的 axes 覆盖；OME 原生通道名称默认保留。单文件 channel_map 按零起始内部通道索引命名，例如 `{"0": "DNA", "1": "Actin"}`，须覆盖所有通道。

lazy 默认为 true，TIFF 只检查头部和布局，不解码像素；PNG 做结构校验。lazy=false 对 TIFF 解码所选 series 和层级后丢弃像素，对 OME-Zarr 逐块解码，不把整个存储累计到内存。解码单个大文件仍可能需要较多内存，lazy 检查通过也不代表全部压缩像素载荷已验证。

## 返回与后续使用

成功返回 status、data、warnings。data 包含 dataset_ref、index_path、dataset_id、n_files、n_input_files、n_fields、channels、metadata_fields、pixel_validation 和 integrity_scope。n_files 是索引条目数，一个多文件 OME 图像算一条；n_input_files 给出关联源成员数。响应不包含像素数组。

images.jsonl 每行保存绝对图像路径、视野标识、shape、axes、dtype、通道和已映射元数据。manifest.json 保存源文件指纹、索引指纹和读取规则。可在另一进程验证：

```python
from tooluniverse.perturbation_image_tool import resolve_image_dataset_ref

descriptor = resolve_image_dataset_ref(result["data"]["dataset_ref"])
with open(descriptor["index"]["path"]) as stream:
    first_image = json.loads(next(stream))
print(first_image)
```

后续工具按索引路径读取像素。源图像、清单或索引改变后引用失效；OME-Zarr 将全部存储成员聚合为目录指纹，块或元数据增删、内容改变均使引用失效，输出目录必须放在输入存储之外；目录模式还检查成员变化。本地加载时会做完整 SHA256 校验，不是完全不读取文件字节，只是不默认解码像素。远端引用只校验选中组和数组元数据，不保证像素块内容不变，返回 integrity_scope=metadata_only 及相应警告。原图和清单保持原样，同名数据集不得覆盖。

接口、错误和验收约定见 [构造规范](TOOL_SPEC.md)。真实图像和测试目录仅保留在本地。


## OME 与金字塔选择

OME-TIFF 和 BigTIFF 使用文件输入。存在多个 series 时必须显式传入 series_index；resolution_level 默认 0，可选择金字塔下采样层。索引中记录 series 总数、层级总数、BigTIFF 标记、OME 通道及物理像素尺寸。选择下采样层后，物理像素尺寸按该层与基础层的大小比例换算。

```json
{
  "source": "/data/image.ome.tiff",
  "source_type": "file",
  "output_dir": "/tmp/beam-image-datasets",
  "series_index": 0,
  "resolution_level": 1
}
```

如果只有一个层级，不要传 resolution_level=1。series 和层级均为零起始；超出范围时报参数错误，不自动退回第一幅。

本地 OME-Zarr 也用 source_type=file，source 指向整个存储目录。group_path 可选择整板或容器内部的一个图像组；当前不自动遍历整板。

```json
{
  "source": "/data/plate.ome.zarr",
  "source_type": "file",
  "output_dir": "/tmp/beam-image-datasets",
  "format": "ome_zarr",
  "group_path": "A/1/0",
  "resolution_level": 0
}
```

图像组必须有 multiscales、明确的 TCZYX 轴和层级坐标变换。保存的 image_metadata 包含 NGFF 版本、轴单位、所选 array_path、scale、translation 和通道元数据。支持 NGFF 0.4（Zarr 2）与 0.5（Zarr 3）；未声明支持的版本明确拒绝。group_path 不得跳出存储。远端 URL 使用下一节的 remote 模式。

当前用真实 BBBC 像素封装的金字塔 BigTIFF 核验容器和层级读取，并不代表完成超过 4 GiB 文件的性能压测。官方 OME-TIFF 样例来源为 [OME 图像资源](https://downloads.openmicroscopy.org/images/OME-TIFF/2016-06/MitoCheck/)。坐标元数据遵循 [NGFF 0.4](https://ngff.openmicroscopy.org/0.4/) 和 [NGFF 0.5](https://ngff.openmicroscopy.org/0.5/)。


## 仪器原生文件

CZI、ND2、LIF 使用 source_type=file，也可放在目录或清单中。format 可省略或设为 czi、nd2、lif。series_index 对应 CZI 场景、ND2 位置或 LIF 内部图像；多个场景／位置／图像时须明确选择，不能静默取第一幅。

CZI 支持其后端提供的金字塔 resolution_level，并保留通道名称、坐标偏移和带单位的物理尺寸。ND2 和 LIF 当前只接受 resolution_level=0。LIF 的多 mosaic 位置另需 position_index，当前支持 XY 平面、TCZ 维度及统一的 8/16 位通道存储；特殊扫描维度和更高位深会报不支持。LIF 后端未提供可靠通道名称时保留数字索引，不将颜色表当成生物学通道名。

ND2 库需要 mmap。如果共享文件系统返回 ENODEV／ENOSYS／ENOTSUP，工具会复制到安全的本地临时目录，读取后清理；image_metadata.staged_local_copy 标记该行为，引用仍登记原始文件。此时需要与源文件大小相当的临时磁盘空间；默认 lazy 不等于不复制文件。其他读取失败不会被不加区分地转为复制重试。

```json
{
  "source": "/data/acquisition.czi",
  "source_type": "file",
  "output_dir": "/tmp/beam-image-datasets",
  "series_index": 0,
  "lazy": true
}
```

## 多文件 OME TIFF

以文件集任一入口文件调用。工具在像素拼装前遍历 OME UUID／FileName 引用，要求关联成员存在、UUID 匹配且位于入口文件目录之内。随后登记所有关联文件指纹，索引 related_files 保存成员路径。缺失成员立即失败，不能由底层库用零补齐。

目录模式自动合并同一完整文件集的重复入口，只返回一条图像索引；清单模式重复声明同一文件集则报错。关联成员修改或丢失后引用失效。原始文件名应保持一致；若移动或重命名文件集，需同步合法的 OME 文件引用。无 FileName 的外部 UUID 自动搜索及目录外关联尚不支持。

## 远端 OME Zarr

远端模式读取一个匿名 HTTP(S) 或 S3 图像组，默认只取元数据，不先下载全部像素。remote_timeout 默认为每次请求 30 秒，可设至 300 秒；它不是整个工具的执行时限。lazy=false 按块解码所选层级，可能产生大量网络流量。

```json
{
  "source": "https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0062A/6001240.zarr",
  "source_type": "remote",
  "output_dir": "/tmp/beam-image-datasets",
  "resolution_level": 2,
  "remote_timeout": 30,
  "lazy": true
}
```

S3-compatible 服务可使用 source=s3://bucket/path 和 remote_endpoint=https://service.example.org；仅匿名读取，安装 s3fs 后启用。URL 和 endpoint 不接受内嵌用户名、密码、查询参数或 fragment，凭据不会写入引用；需要登录或签名链接的数据应先下载到本地。

远端 manifest 保存被读取的组和所选层级数组元数据摘要。源像素块在元数据不变时被替换，resolve_image_dataset_ref 不会发现，因此不能把 metadata_only 当成完整内容验证。需要完整内容指纹时先准备本地存储，再用本地 file 模式。

本地真实验证样例包括 OME 官方两文件 tubhiswt、IDR CZI、Nikon ND2 和 readlif 分发的 LIF；样例数据、下载来源记录和运行产物均保留本地，不随工具提交。
