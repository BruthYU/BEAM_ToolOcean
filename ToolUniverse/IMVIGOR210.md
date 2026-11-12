# IMvigor210：先下载，再查询

四个工具读取本地官方 IMvigor210CoreBiologies 1.0.1 快照。首次使用须执行以下准备流程；查询本身不访问网络、不下载数据。数据无需上传 GitHub。

## 1. 安装工具及准备脚本依赖

在 BEAM_ToolOcean 仓库根目录、同一个 Python 环境中执行：

```bash
python -m pip install -e ./ToolUniverse
python -m pip install 'rdata==0.11.2'
```

rdata 仅供一次性转换使用，不需要安装 R。numpy 和 pandas 已在 ToolUniverse 的依赖中。

## 2. 下载官方数据并转换

```bash
python ToolUniverse/scripts/prepare_imvigor210.py \
  --output /absolute/path/dataset/IMvigor210/query_data
```

脚本下载官方约74 MB的数据包，检查包和 cds.RData 的固定 SHA256，再导出 SQLite、NumPy计数矩阵及 manifest。原站当前使用 HTTP，内容由固定 SHA256核验；无需API key。

如果服务器无法访问原站，可在本地下载官方包后传入服务器：

```bash
curl -fL --retry 3 \
  http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.1.tar.gz \
  -o IMvigor210CoreBiologies_1.0.1.tar.gz
python ToolUniverse/scripts/prepare_imvigor210.py \
  --archive /absolute/path/IMvigor210CoreBiologies_1.0.1.tar.gz \
  --output /absolute/path/dataset/IMvigor210/query_data
```

输出目录必须不存在，脚本不覆盖现有数据。转换失败不会留下可被查询工具误用的半成品目录；可修复网络或依赖后重试。脚本不下载原始FASTQ。

```text
query_data/
├── imvigor210.sqlite
├── counts.npy
└── manifest.json
```

实际快照包含348份样本、347个不同患者ID、31286个基因。manifest记录来源、版本、单位、字段映射和文件校验值。

## 3. 配置数据路径

```bash
export TOOLUNIVERSE_IMVIGOR210_DATA=/absolute/path/dataset/IMvigor210/query_data
```

通过 MCP 或其他后台服务启动 ToolUniverse 时，也要在该进程的环境中配置此变量。

## 4. 查询

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse(load_workspace=False)
tu.load_tools()

def call(name, arguments):
    result = tu.run_one_function({"name": name, "arguments": arguments}, use_cache=False)
    if result.get("status") != "success":
        raise RuntimeError(result)
    return result["data"]

info = call("IMvigor210_get_dataset_info", {})
print(info["samples"], info["patients"], info["genes"])

page = call("IMvigor210_query_samples", {
    "filters": {"response": ["PD"], "immune_phenotype": ["excluded"]},
    "fields": ["sample_id", "patient_id", "response"],
    "limit": 50,
    "offset": 0,
})
print(page)
# 返回 next_offset 时可继续分页。
if page["rows"]:
    expression = call("IMvigor210_get_expression", {
        "genes": ["CD274", "CXCL9"],
        "sample_ids": [row["sample_id"] for row in page["rows"]],
        "expression_type": "raw_counts",
    })
    print(expression)

print(call("IMvigor210_summarize_cohort", {
    "filters": {"immune_phenotype": ["excluded"]},
    "group_by": "response",
}))
tu.close()
```

不同筛选字段取AND，同字段列表取OR；使用`None`匹配缺失值。字段名区分大小写，从dataset_info读取。应答标签包括CR、PR、SD、PD、NE。原始计数不是TPM/FPKM；缺失值不填零。表达接口最多100个基因、500个样本；基因优先匹配Entrez ID，再精确匹配symbol。同名基因返回全部匹配项，缺失基因列入missing_genes。

工具仅提供数据查询与描述性统计，不自动进行差异表达、生存分析或临床预测。共享存储不支持mmap时会回退到普通读取。
