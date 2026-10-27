# 虚拟细胞模型综述表（100 篇，以 2025–2026 前沿工作为主）

> 更新日期：2026-10-27
>
> 年份结构：2025–2026 年 75 篇（75%）；2024 年及以前 25 篇（25%）。

## 口径说明

- 本表采用广义“虚拟细胞”口径：模型能够表示细胞状态，并至少支持状态生成、扰动响应、时空演化、跨模态重建、形态预测或机制模拟中的一项。
- “数据尺度”按要求只使用“单细胞”或“Bulk”。空间 spot、组织切片和细胞群模型归入“单细胞”，并在“所需数据”中注明实际分辨率。
- “主要作者单位｜年份”列优先列第一作者或通讯作者的 1–2 个核心单位；不代表完整作者单位清单。
- 年份优先采用正式版本年份；尚未正式发表的工作按首次公开版本列出并标注“预印本”。因此同一模型的新闻发布时间、预印本年份和正式出版年份可能不同。
- 表中 2025–2026 年的快速迭代工作较多。预印本表示“公开可核查”，不等同于已经过同行评议或独立复现。

## A. 前沿扰动模型、世界模型与可执行虚拟细胞（30）

| 论文题目 | 数据尺度 | 组学模态 | 支持任务 | 所需数据 | 主要作者单位｜年份 |
|---|---|---|---|---|---|
| 001. **Predicting Cellular Responses to Perturbation across Diverse Contexts with State**（State） | 单细胞 | 转录组、遗传/化学/细胞因子扰动 | 跨细胞背景扰动预测、细胞状态嵌入、分布生成 | 大规模观察性 scRNA-seq、Perturb-seq/药物扰动、细胞背景与处理标签 | **Arc Institute**｜2026 |
| 002. **Stack: In-Context Learning of Single-Cell Biology**（Stack） | 单细胞 | 转录组、扰动组 | 以细胞集合为提示的跨条件预测、未见背景迁移 | 对照与扰动 scRNA-seq、上下文细胞集合、扰动标签 | **Arc Institute**｜2026（预印本） |
| 003. **PIE: Generalizing Perturbation Effects across Unseen Perturbations, Contexts and Datasets**（PIE） | 单细胞 | 转录组、生物知识库 | 未见扰动×未见背景预测、差异表达和效应量预测 | Perturb-seq、基线表达、通路/PPI/文献知识和背景元数据 | **Arc Institute**｜加州大学旧金山分校｜2026（预印本） |
| 004. **Tahoe-x1: Scaling Perturbation-Trained Single-Cell Foundation Models to 3 Billion Parameters**（Tahoe-x1） | 单细胞 | 转录组、药物扰动 | 药物响应生成、未见细胞系迁移、基因—细胞—化合物表征 | Tahoe-100M 等扰动 scRNA-seq、药物与剂量、细胞系标签 | **Tahoe Therapeutics**｜2025（预印本） |
| 005. **X-Cell: Scaling Causal Perturbation Prediction across Diverse Cellular Contexts via Diffusion Language Models**（X-Cell） | 单细胞 | 转录组、CRISPRi、多模态先验 | 跨背景零样本基因扰动预测、靶点发现、机制优先级排序 | X-Atlas/Pisces Perturb-seq、对照细胞、靶基因及蛋白/PPI等先验 | **Xaira Therapeutics**｜2026（预印本） |
| 006. **Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling** | 单细胞 | 转录组、药物扰动 | 药物—细胞状态建模、机制分群、扰动基础模型训练 | 50 个癌细胞系、约 1,100 个小分子条件和约一亿个单细胞转录组 | **Tahoe Therapeutics**｜Parse Biosciences｜2025（预印本） |
| 007. **SCALE: Scalable Conditional Atlas-Level Endpoint Transport for Virtual Cell Perturbation Prediction**（SCALE） | 单细胞 | 转录组、扰动组 | 条件流匹配、细胞群分布迁移、图谱级训练与推理 | 对照/扰动 scRNA-seq、处理与背景标签；可使用 Tahoe-100M | 上海人工智能实验室｜NVIDIA｜2026（预印本） |
| 008. **Towards Building a World Model to Simulate Perturbation-Induced Cellular Dynamics by AlphaCell**（AlphaCell） | 单细胞 | 全转录组、扰动组 | 全基因状态重建、连续状态转移、未见背景零样本预测 | 全转录组 scRNA-seq、对照/扰动群体、细胞背景与处理标签 | 同济大学｜2026（预印本） |
| 009. **CellOS: Learning a World Model of Cellular State through Joint Embedding Prediction**（CellOS） | 单细胞 | 转录组 | 多视图细胞表示、批次整合、扰动响应迁移 | 大规模 scRNA-seq、表达排序与群体信息排序视图、任务标签 | Vitaura｜中国科学院动物研究所｜2026（预印本） |
| 010. **VCWorld: A Biological World Model for Virtual Cell Simulation**（VCWorld） | 单细胞 | 转录组、通路/PPI/文本知识 | 白盒扰动推理、差异表达方向预测、机制链生成 | 扰动表达、药物/基因描述、通路与相互作用知识图谱 | 上海交通大学｜NeoLife AI｜2026 |
| 011. **VCLMU: Mechanism-Centric Virtual Cell World Modeling for Perturbation Response**（VCLMU） | 单细胞 | 转录组、多模态知识先验 | 潜在机制单元状态转移、未见遗传扰动预测、机制解释 | pseudo-Bulk 扰动谱、单细胞扰动数据、PPI/通路/功能证据 | UTHealth Houston｜2026（预印本） |
| 012. **Lingshu-Cell: A Generative Cellular World Model for Transcriptome Modeling toward Virtual Cells**（Lingshu-Cell） | 单细胞 | 全转录组、遗传/细胞因子扰动 | 离散扩散生成、全转录组模拟、跨组织/物种与组合条件预测 | 多物种 scRNA-seq、细胞身份/供体、扰动条件及全基因计数 | 阿里巴巴达摩院｜2026（预印本） |
| 013. **HarmonyCell: Automating Single-Cell Perturbation Modeling under Semantic and Distribution Shifts**（HarmonyCell） | 单细胞 | 转录组、扰动元数据 | 自动元数据统一、模型搜索、跨数据集分布迁移 | 异构扰动 scRNA-seq、数据字典、背景/扰动/剂量字段 | 上海人工智能实验室｜2026（预印本） |
| 014. **OCOO-T: A Simple and Scalable Virtual Cell Model for Transcriptional Perturbation Response Prediction**（OCOO-T） | 单细胞 | 转录组、遗传/化学/细胞因子扰动 | 连续时间去噪、跨细胞类型与剂量响应预测 | 对照/处理 scRNA-seq、扰动、剂量和细胞背景标签 | INFevo｜晶泰科技｜2026（预印本） |
| 015. **scDFM: Distributional Flow Matching Model for Robust Single-Cell Perturbation Prediction**（scDFM） | 单细胞 | 转录组、基因关系图 | 细胞群分布生成、未见/组合扰动预测 | 对照与扰动 scRNA-seq、靶基因、基因相互作用图 | 西湖大学｜2026 |
| 016. **SAVE: A Generalizable Framework for Multi-Condition Single-Cell Generation with Gene Block Attention**（SAVE） | 单细胞 | 转录组、多条件元数据 | 条件生成、批次校正、扰动预测、合成细胞 | 多条件 scRNA-seq、批次/细胞类型/处理标签 | 浙江大学｜2026 |
| 017. **TxPert: Using Multiple Knowledge Graphs for Prediction of Transcriptomic Perturbation Effects**（TxPert） | 单细胞 | 转录组、多知识图谱 | 未见单/双基因扰动、跨新细胞系预测 | Perturb-seq、对照表达、PPI/通路/共表达图 | Helmholtz Munich｜慕尼黑工业大学｜2026 |
| 018. **UniPert-G2CP Bridges Genetic and Chemical Screens from Molecular Representation to Phenotype Modeling**（UniPert-G2CP） | 单细胞 | 转录组、遗传扰动、化学结构 | 遗传筛选到化合物筛选迁移、药物表型预测 | CRISPR 与药物扰动表达、靶点、化合物结构及表型 | 浙江大学｜2026 |
| 019. **Reinforcement Learning-Based Design of Sequential Drug Treatment Targeting the Evolving Tumour Landscape with SequenTx**（SequenTx） | Bulk | 肿瘤转录组、药物、时间序列 | 序贯用药设计、肿瘤状态演化与耐药控制 | 时序肿瘤表达、药物作用与剂量、体内外验证读出 | 同济大学｜2026 |
| 020. **scDRP: Disentangled Representation Learning for Predicting Single-Cell Responses to Perturbations and Estimating Individual Treatment Effects**（scDRP） | 单细胞 | 转录组 | 个体化处理效应、反事实配对、剂量与组合外推 | 未配对对照/处理 scRNA-seq、协变量和处理标签 | 哈尔滨工业大学｜2026 |
| 021. **Unbalanced Perturbation Dynamics for Cell Fate Design**（U-Pert） | 单细胞 | 转录组、遗传/药物/细胞因子扰动 | 同时预测表达和细胞数量、细胞命运逆向设计 | 对照/处理 scRNA-seq、细胞丰度、剂量、背景与目标状态 | 清华大学｜2026（预印本） |
| 022. **Modeling Cell-Cycle-Aware Single-Cell Drug Perturbation Responses**（scCycleMol） | 单细胞 | 转录组、化学结构、细胞周期 | 药物响应、细胞周期变化、未见药物/剂量预测 | sci-Plex3、SMILES、剂量、细胞系和细胞周期标签 | 密歇根大学｜2026（预印本） |
| 023. **Dynamic Graph Prior-Guided Latent Adaptation for Generalizable Single-Cell Perturbation Response Prediction**（DGP-LA） | 单细胞 | 转录组、药物结构、动态图先验 | 未见细胞类型药物响应、知识引导状态迁移 | Tahoe-100M/免疫扰动数据、药物、背景与动态基因图 | 华南理工大学｜2026 |
| 024. **scPPDM: A Diffusion Model for Single-Cell Drug-Response Prediction**（scPPDM） | 单细胞 | 转录组、药物与剂量 | 未见药物、未见协变量组合的扩散式响应生成 | Tahoe-100M、对照细胞、药物、剂量和细胞系 | 南京大学｜2025（预印本） |
| 025. **CFM-GP: Unified Conditional Flow Matching to Learn Gene Perturbation across Cell Types**（CFM-GP） | 单细胞 | 转录组、遗传扰动 | 跨细胞类型条件流、未见靶基因响应预测 | 多细胞类型 Perturb-seq、对照、靶基因和类型标签 | 加州大学圣地亚哥分校｜2025（预印本） |
| 026. **GPerturb: Gaussian Process Modelling of Single-Cell Perturbation Data**（GPerturb） | 单细胞 | 转录组、扰动组 | 稀疏扰动效应估计、响应基因发现、不确定性量化 | 对照/处理 scRNA-seq、扰动和协变量标签 | 赫尔辛基大学｜2025 |
| 027. **Morphologically Intelligent Perturbation Prediction with FORM**（FORM） | 单细胞 | 3D 显微图像、形态组、信号读出 | 3D 形态变化、组合扰动、未见扰动与信号活性预测 | 扰动前后 3D 细胞图像、处理标签、信号通路读出 | Allen Institute for Cell Science｜2025（预印本） |
| 028. **PertAdapt: Unlocking Single-Cell Foundation Models for Genetic Perturbation Prediction via Condition-Sensitive Adaptation**（PertAdapt） | 单细胞 | 转录组、遗传扰动 | 轻量适配基础模型、单/双基因未见扰动预测 | 预训练细胞模型、Perturb-seq、靶基因与条件标签 | 北京大学｜2025（预印本） |
| 029. **Unified Multimodal Learning Enables Generalized Cellular Response Prediction to Diverse Perturbations**（X-Pert） | 单细胞 | 转录组、遗传/化学扰动、剂量 | 统一多扰动类型、组合与剂量响应预测 | 基因和药物扰动 scRNA-seq、SMILES、剂量与细胞背景 | 中国科学院｜2025（预印本） |
| 030. **Virtual-Cell Verification Enables Self-Auditing AI Discovery for Immune Rejuvenation**（CellQ） | 单细胞 | 转录组、文本、遗传扰动 | 少样本扰动预测、模型自审计、免疫年轻化靶点发现 | 大规模 scRNA-seq 预训练、扰动示例、自然语言问题与验证实验 | Shift Bioscience｜2026（预印本） |

## B. 单细胞与多组学基础模型（20）

| 论文题目 | 数据尺度 | 组学模态 | 支持任务 | 所需数据 | 主要作者单位｜年份 |
|---|---|---|---|---|---|
| 031. **Universal Cell Embeddings: A Foundation Model for Cell Biology**（UCE） | 单细胞 | 转录组、蛋白序列先验 | 跨物种/组织细胞表征、细胞类型标注、零样本迁移 | 多物种大规模 scRNA-seq、基因同源关系和蛋白序列嵌入 | Stanford University｜2026 |
| 032. **TranscriptFormer: A Generative Foundation Model for Transcriptomes**（TranscriptFormer） | 单细胞 | 转录组、多物种图谱 | 跨物种细胞表示、表达生成、细胞类型与扰动任务迁移 | 多物种单细胞转录组、基因同源映射、组织与细胞注释 | Chan Zuckerberg Initiative｜CZ Biohub｜2026 |
| 033. **RegFormer: A Single-Cell Foundation Model Powered by Gene Regulatory Hierarchies**（RegFormer） | 单细胞 | 转录组、基因调控网络 | 细胞表征、类型注释、批次整合、调控关系和扰动预测 | scRNA-seq、转录因子—靶基因层级及细胞标签 | BGI Research｜Beijing Graphene Institute｜2026 |
| 034. **CAPTAIN: A Multimodal Foundation Model Pretrained on Co-Assayed Single-Cell RNA and Protein**（CAPTAIN） | 单细胞 | 转录组、蛋白质组 | RNA/蛋白联合表征、缺失模态补全、细胞注释和跨数据迁移 | CITE-seq 等同细胞 RNA+表面蛋白数据及细胞标签 | 北京大学｜中国科学院｜2026 |
| 035. **HoloCell: A Generative Foundation Model for Holistic Cellular Modeling**（HoloCell） | 单细胞 | 转录组、多组学先验 | 细胞状态生成、条件模拟、跨任务迁移与虚拟细胞表征 | 大规模 scRNA-seq、细胞/组织/条件元数据和任务标签 | 哈尔滨工业大学｜2026（预印本） |
| 036. **CLM-X: A Generalist Cellular Language Model for Single-Cell Biology**（CLM-X） | 单细胞 | 转录组、文本/知识先验 | 通用细胞表示、注释、检索、扰动与跨数据集推理 | 大规模 scRNA-seq、基因词表、细胞标签和生物医学知识 | 中国科学院杭州医学研究所｜2026（预印本） |
| 037. **xVERSE: A Transcriptomics-Native Foundation Model for Universal Cell Representation and Virtual Cell Synthesis**（xVERSE） | 单细胞 | 转录组 | 通用细胞嵌入、条件细胞生成、跨组织/疾病迁移 | 多来源 scRNA-seq、组织/疾病/细胞类型元数据 | Duke University｜2026（预印本） |
| 038. **scLong: A Billion-Parameter Foundation Model for Capturing Long-Range Gene Context in Single-Cell Transcriptomics**（scLong） | 单细胞 | 转录组 | 长程基因依赖建模、细胞注释、批次整合与扰动预测 | 大规模 scRNA-seq、基因排序/表达分箱和下游标签 | 上海交通大学｜2026 |
| 039. **BulkFormer: A Large-Scale Foundation Model for Bulk Transcriptomes**（BulkFormer） | Bulk | 转录组 | Bulk 样本表征、疾病分类、预后和药物响应迁移 | 大规模 Bulk RNA-seq/芯片表达、组织/疾病/临床标签 | Korea Advanced Institute of Science and Technology｜2026（预印本） |
| 040. **Scalable Single-Cell Gene Expression Generation with Latent Diffusion Models**（scLDM） | 单细胞 | 转录组 | 图谱级细胞生成、稀有状态扩增、条件表达模拟 | 大规模 scRNA-seq、细胞类型与组织/疾病条件标签 | Chan Zuckerberg Initiative｜2025（预印本） |
| 041. **EpiAgent: A Generative Foundation Model for Single-Cell Epigenomics**（EpiAgent） | 单细胞 | 表观组（scATAC-seq） | 染色质可及性生成、细胞表征、类型标注和跨批次迁移 | 单细胞 ATAC-seq peaks/fragments、基因组位置和细胞标签 | 清华大学｜2025 |
| 042. **EpiFoundation: A Foundation Model for Single-Cell Epigenomics**（EpiFoundation） | 单细胞 | 表观组（scATAC-seq） | 可及性表征、细胞类型注释、调控元件与跨组织迁移 | 大规模 scATAC-seq、peak/基因组序列及细胞元数据 | University of Texas at Dallas｜2025（预印本） |
| 043. **Language May Be All Omics Needs: Harmonizing Multimodal Data for Omics Understanding with CellHermes**（CellHermes） | 单细胞 | 转录组、表观组、蛋白质组、多组学 | 异构组学统一编码、跨模态检索/补全、细胞注释 | scRNA-seq、scATAC-seq、CITE-seq 等及模态/细胞标签 | 哈尔滨工业大学｜2025（预印本） |
| 044. **SCARF: Single-Cell ATAC-seq and RNA-seq Foundation Model**（SCARF） | 单细胞 | 转录组、表观组 | RNA/ATAC 联合嵌入、跨模态预测、细胞类型与调控分析 | 配对或非配对 scRNA-seq/scATAC-seq、基因—peak 映射 | 中国科学院｜2025（预印本） |
| 045. **CellFM: A Large-Scale Foundation Model Pre-Trained on Transcriptomics of 100 Million Human Cells**（CellFM） | 单细胞 | 转录组 | 细胞注释、聚类、批次校正、基因网络与扰动任务 | 约一亿细胞的 scRNA-seq 预训练库及任务标签 | 重庆大学｜中山大学｜2025 |
| 046. **scPRINT: A Large Cell Model for Gene Network Inference from Single-Cell Data**（scPRINT） | 单细胞 | 转录组、基因网络 | 细胞嵌入、细胞类型注释、基因调控网络推断 | 大规模 scRNA-seq、基因 ID、批次与细胞元数据 | Institut Pasteur｜Université Paris Cité｜2025 |
| 047. **SCimilarity: A Foundation Model for Single-Cell Transcriptomics**（SCimilarity） | 单细胞 | 转录组 | 细胞相似性检索、类型注释、异常/疾病状态发现 | 大规模带注释 scRNA-seq 参考图谱和查询细胞 | Genentech｜Broad Institute｜2025 |
| 048. **scNET: Learning Context-Specific Gene and Cell Embeddings by Single-Cell Foundation Models**（scNET） | 单细胞 | 转录组、基因网络 | 基因/细胞嵌入、注释、基因模块与网络发现 | 大规模 scRNA-seq、基因共表达上下文及下游标签 | Tel Aviv University｜2025 |
| 049. **scKGBERT: A Knowledge-Graph-Enhanced Language Model for Single-Cell Transcriptomics**（scKGBERT） | 单细胞 | 转录组、知识图谱 | 细胞类型注释、基因关系建模、可解释细胞表示 | scRNA-seq、基因本体/PPI/通路知识图谱和标签 | Northeast Forestry University｜KAUST｜2025 |
| 050. **A Foundation Model of Transcription across Human Cell Types**（GET） | 单细胞 | 表观组、转录组调控读出 | 顺式调控活性、基因表达和未见细胞类型转录预测 | 单细胞染色质可及性、DNA 序列、基因表达与细胞类型 | Columbia University｜MBZUAI｜2025 |

## C. 空间组织、细胞生态与多细胞虚拟环境（15）

| 论文题目 | 数据尺度 | 组学模态 | 支持任务 | 所需数据 | 主要作者单位｜年份 |
|---|---|---|---|---|---|
| 051. **Nicheformer: A Foundation Model for Single-Cell and Spatial Omics**（Nicheformer） | 单细胞 | 转录组、空间转录组 | 空间域/生态位识别、细胞类型注释、空间迁移学习 | 大规模 scRNA-seq 与空间转录组、坐标、组织和细胞标签 | Helmholtz Munich｜TUM｜2025 |
| 052. **SToFM: A Multi-Scale Foundation Model for Spatial Transcriptomics**（SToFM） | 单细胞 | 空间转录组、组织学图像 | 跨平台空间表征、域识别、基因补全与组织迁移 | spot/单细胞空间表达、坐标、可选 H&E 图像和区域标签 | 清华大学｜PharMolix｜2025（预印本） |
| 053. **Quantitative Characterization of Cell Niches in Spatially Resolved Omics Data**（NicheCompass） | 单细胞 | 空间转录组、通路知识 | 生态位识别、空间基因程序、细胞通讯与条件比较 | 单细胞分辨空间表达、坐标、邻接图和基因程序先验 | Helmholtz Munich｜TUM｜2025 |
| 054. **Novae: A Graph-Based Foundation Model for Spatial Transcriptomics Data**（Novae） | 单细胞 | 空间转录组 | 跨样本组织域划分、空间表示、参考图谱迁移 | 多切片空间表达、空间邻接图和可选区域注释 | CentraleSupélec｜Gustave Roussy｜2025 |
| 055. **spaGFM Is a Scalable Graph Foundation Model for Spatial Transcriptomics Analyses**（spaGFM） | 单细胞 | 空间转录组 | 可扩展空间嵌入、域检测、基因插补与跨技术迁移 | 多平台空间表达、坐标/图结构和下游标签 | Indiana University Indianapolis｜Houston Methodist｜2026 |
| 056. **STARNet: A Foundation Model for Spatial Transcriptomics**（STARNet） | 单细胞 | 空间转录组、组织学图像 | 空间表征、组织区域分割、基因表达/功能预测 | 空间转录组、坐标、H&E 图像与区域/细胞类型标签 | 西湖大学｜2025（预印本） |
| 057. **SpatialProp: Tissue Perturbation Modeling with Spatially Resolved Single-Cell Transcriptomics**（SpatialProp） | 单细胞 | 空间转录组、扰动组 | 组织内扰动传播、空间响应预测、邻域效应分解 | 对照/扰动空间单细胞表达、坐标、组织和处理标签 | Stanford University｜2025（预印本） |
| 058. **Inferring Virtual Cell Environments Using Multi-Agent Reinforcement Learning**（CellTRIP） | 单细胞 | 转录组、空间转录组 | 从表达反推虚拟空间、细胞邻域重建、组织结构模拟 | scRNA-seq 或空间表达、细胞类型及可选坐标/邻接约束 | University of Wisconsin–Madison｜2025 |
| 059. **The Virtual Tissues Foundation Model Resolves Spatial Proteomics across Scales**（VirTues） | 单细胞 | 空间蛋白质组、显微图像 | 多尺度组织表示、细胞/区域表型、跨队列预后与迁移 | 高多重成像蛋白数据、细胞坐标、组织图像和临床标签 | EPFL｜2026 |
| 060. **STPath: A Generative Foundation Model for Integrating Spatial Transcriptomics and Whole-Slide Images**（STPath） | 单细胞 | 空间转录组、病理全切片图像 | 图像—表达对齐、空间表达生成、区域表征和病理迁移 | 配对空间转录组与 H&E 全切片图、坐标和组织标签 | 香港中文大学｜2025 |
| 061. **FmH2ST: Foundation Model-Based Spatial Transcriptomics Generation from Histological Images**（FmH2ST） | 单细胞 | 空间转录组、组织学图像 | 从 H&E 预测空间表达、超分辨生成、组织区域分析 | 配对 H&E 图像与 spot 表达、坐标及基因面板 | 电子科技大学｜2025 |
| 062. **SPATIA: Multimodal Generation and Prediction of Spatial Cell Phenotypes**（SPATIA） | 单细胞 | 空间转录组、空间蛋白质组、显微图像 | 空间细胞表型生成、缺失模态预测、组织状态分类 | 配对或对齐的空间表达/蛋白/图像、坐标和细胞标签 | University of Cambridge｜2025（预印本） |
| 063. **HEIST: A Graph Foundation Model for Spatial Transcriptomics and Proteomics Data**（HEIST） | 单细胞 | 空间转录组、空间蛋白质组 | 跨模态空间图表示、域识别、缺失模态补全 | 空间 RNA/蛋白表达、细胞坐标、邻接图和区域标签 | Harvard Medical School｜2025（预印本） |
| 064. **Multi-Scale Modeling of Human Tissues from Spatial Transcriptomics with TERRA**（TERRA） | 单细胞 | 空间转录组 | 细胞—生态位—组织多尺度表示、跨组织预测和检索 | 多组织空间表达、细胞坐标、样本/疾病与区域元数据 | Broad Institute｜2026（预印本） |
| 065. **Bridging the Dimensional Gap from Planar Spatial Transcriptomics to 3D Cell Atlases**（SpatialZ） | 单细胞 | 空间转录组、3D 组织结构 | 连续切片配准、三维细胞图谱重建、空间域分析 | 多张连续空间转录组切片、坐标、形态图像和配准标记 | ETH Zürich｜2025 |

## D. 显微图像、超微结构与机制型整细胞模型（10）

| 论文题目 | 数据尺度 | 组学模态 | 支持任务 | 所需数据 | 主要作者单位｜年份 |
|---|---|---|---|---|---|
| 066. **An Operational Perturbation Proteomics-Based Virtual Cell Model**（ProteinTalks） | 单细胞 | 蛋白质组、扰动组 | 蛋白状态表征、扰动响应预测、通路和靶点机制解释 | 扰动蛋白质组、蛋白互作网络、处理条件及表型读出 | 清华大学｜2026 |
| 067. **Bringing the Genetically Minimal Cell to Life on a Computer in 4D**（4D Syn3A） | 单细胞 | 基因组、转录组、蛋白质组、代谢组、3D 超微结构 | 最小细胞全生命周期模拟、分裂、代谢与空间动力学联动 | JCVI-syn3A 全基因组、多组学参数、cryo-ET/显微结构及生理数据 | University of Illinois Urbana-Champaign｜J. Craig Venter Institute｜2026 |
| 068. **Abundant Positively Charged Proteins Underlie JCVI-syn3A's Expanded Nucleoid and Ribosome Distribution** | 单细胞 | 超微结构、蛋白质组、基因组 | 核区与核糖体空间组织模拟、分子机制检验 | 最小细胞蛋白组成、分子结构、cryo-ET 约束和粗粒化力场 | University of Illinois Urbana-Champaign｜2026 |
| 069. **CellFluxV2: An Image Generative Foundation Model for Virtual Cell Modeling**（CellFluxV2） | 单细胞 | 显微图像、形态组、扰动元数据 | 条件细胞图像生成、形态扰动预测、时序状态模拟 | 大规模细胞显微图像、处理/细胞系/时间标签及可选分子条件 | University of California, San Francisco｜2026（预印本） |
| 070. **Virtual Experiments Bridge Sequence and Microscopy with Generative Models**（CELL-FM） | 单细胞 | 显微图像、蛋白序列、形态组 | 序列到定位/形态预测、虚拟成像、突变效应模拟 | 蛋白序列、带通道标注的细胞图像、定位与扰动标签 | University of California, San Francisco｜2026（预印本） |
| 071. **A Foundation AI Model Enhances Electron Microscopy Image Analysis**（DF5T） | 单细胞 | 电子显微图像、超微结构 | EM 分割、检测、表征迁移和少样本分析 | 大规模 2D/3D 电镜图像、像素/实例标注或自监督图块 | Howard Hughes Medical Institute｜2026（预印本） |
| 072. **Unifying the Electron Microscopy Multiverse through a Large-Scale Foundation Model**（EM-DINO / OmniEM） | 单细胞 | 电子显微图像、超微结构 | 跨器官/尺度 EM 表征、分割、检索与域迁移 | 多来源大规模电镜图像、图块、可选结构分割与数据域标签 | European Molecular Biology Laboratory｜2025（预印本） |
| 073. **EMCF Ecosystem: Towards Pretrained Foundation Model for Electron Microscopy Image Analysis**（EMCF） | 单细胞 | 电子显微图像、超微结构 | 通用 EM 编码、分割、分类、检测与模型基准 | 跨模态/跨尺度 EM 图像、预训练图块和任务标注 | 中国科学院自动化研究所｜2025（预印本） |
| 074. **CellSAM: A Foundation Model for Cell Segmentation**（CellSAM） | 单细胞 | 显微图像 | 跨成像模态细胞检测与分割、零/少样本迁移 | 荧光、明场、相差或组织显微图像及细胞实例掩膜 | Brigham and Women's Hospital｜Harvard Medical School｜2025 |
| 075. **MONET—Virtual Cell Painting of Brightfield Images and Time Lapses Using Reference-Consistent Diffusion**（MONET） | 单细胞 | 明场/荧光显微图像、时序影像 | 虚拟染色、时间序列荧光生成、形态动态追踪 | 配对明场—荧光图像、时间序列和细胞/通道标识 | Helmholtz Munich｜2025（预印本） |

## E. 经典基础与 2024 年以前代表性模型（25）

| 论文题目 | 数据尺度 | 组学模态 | 支持任务 | 所需数据 | 主要作者单位｜年份 |
|---|---|---|---|---|---|
| 076. **E-CELL: Software Environment for Whole-Cell Simulation**（E-CELL） | 单细胞 | 基因组、代谢组、生化反应网络 | 整细胞动力学模拟、代谢/基因表达过程联动 | 基因组注释、反应网络、酶动力学参数和初始浓度 | Keio University｜1999 |
| 077. **Virtual Cell: A Software Environment for Computational Cell Biology**（VCell） | 单细胞 | 显微图像、生化反应网络、空间几何 | 反应—扩散、膜转运、电生理与空间细胞模拟 | 细胞几何、反应式、扩散/动力学参数和边界条件 | University of Connecticut Health Center｜2001 |
| 078. **A Whole-Cell Computational Model Predicts Phenotype from Genotype**（M. genitalium 全细胞模型） | 单细胞 | 基因组、转录组、蛋白质组、代谢组 | 基因型到表型模拟、基因敲除、细胞周期和生长预测 | 最小细菌多组学、反应/动力学参数和单细胞验证数据 | Stanford University｜2012 |
| 079. **Simultaneous Cross-Evaluation of Heterogeneous E. coli Datasets via Mechanistic Simulation**（E. coli 全细胞模型） | 单细胞 | 基因组、转录组、蛋白质组、代谢组 | 细菌生长、环境响应、跨数据集机制一致性检验 | E. coli 多组学、代谢/表达参数、培养条件与表型数据 | Stanford University｜SRI International｜2020 |
| 080. **Fundamental Behaviors Emerge from Simulations of a Living Minimal Cell**（Syn3A） | 单细胞 | 基因组、代谢组、空间反应网络 | 最小细胞三维动力学、生长、能量代谢与分子扩散 | JCVI-syn3A 基因组、动力学/结构参数和细胞物理数据 | University of Illinois Urbana-Champaign｜J. Craig Venter Institute｜2022 |
| 081. **PhysiCell: An Open Source Physics-Based Cell Simulator for 3-D Multicellular Systems** | 单细胞 | 细胞表型、空间/力学、微环境 | 三维多细胞生长、细胞周期、死亡、迁移和药物扩散 | 细胞规则、力学参数、微环境场和初始空间分布 | Indiana University｜2018 |
| 082. **CompuCell3D: Modeling of Multicellular Systems with Complex Cell Behaviors** | 单细胞 | 细胞表型、空间/力学、信号网络 | 多细胞组织形态发生、黏附、迁移与反应—扩散 | Potts 模型参数、细胞类型/规则、场方程和组织几何 | Indiana University｜2012 |
| 083. **Deep Generative Modeling for Single-Cell Transcriptomics**（scVI） | 单细胞 | 转录组 | 批次整合、降维、差异表达、去噪与生成 | UMI count 矩阵、批次和可选细胞标签 | University of California, Berkeley｜2018 |
| 084. **scGen Predicts Single-Cell Perturbation Responses**（scGen） | 单细胞 | 转录组、扰动组 | 跨细胞类型扰动响应预测、潜空间向量算术 | 对照/处理 scRNA-seq、细胞类型和处理标签 | Helmholtz Munich｜TUM｜2019 |
| 085. **Predicting Cellular Responses to Complex Perturbations in High-Throughput Screens**（CPA） | 单细胞 | 转录组、化学结构、遗传扰动 | 药物/基因、剂量、时间和组合响应外推 | 对照/扰动 scRNA-seq、SMILES/靶基因、剂量与协变量 | Helmholtz Munich｜TUM｜2023 |
| 086. **Learning Single-Cell Perturbation Responses Using Neural Optimal Transport**（CellOT） | 单细胞 | 转录组、扰动组 | 未配对对照到处理分布映射、反事实响应预测 | 对照/处理 scRNA-seq、处理条件和细胞协变量 | ETH Zürich｜2023 |
| 087. **Predicting Transcriptional Outcomes of Novel Multigene Perturbations with GEARS**（GEARS） | 单细胞 | 转录组、遗传扰动、基因知识图 | 未见单/多基因扰动与组合效应预测 | Perturb-seq、基因相互作用图、靶基因和对照细胞 | Stanford University｜2024 |
| 088. **scGPT: Toward Building a Foundation Model for Single-Cell Multi-Omics Using Generative AI**（scGPT） | 单细胞 | 转录组、表观组、蛋白质组、多组学 | 细胞注释、批次整合、多组学整合和扰动预测 | 大规模单细胞表达预训练数据及任务特定组学/标签 | University of Toronto｜Vector Institute｜2024 |
| 089. **Transfer Learning Enables Predictions in Network Biology**（Geneformer） | 单细胞 | 转录组、基因网络 | 细胞状态分类、基因剂量敏感性与网络扰动预测 | 大规模 scRNA-seq、表达排序和疾病/功能标签 | Harvard Medical School｜Broad Institute｜2023 |
| 090. **Large-Scale Foundation Model on Single-Cell Transcriptomics**（scFoundation） | 单细胞 | 转录组 | 细胞和基因嵌入、注释、药物响应与扰动任务 | 大规模人类 scRNA-seq、基因表达和下游任务标签 | Tsinghua University｜BioMap｜2024 |
| 091. **Joint Probabilistic Modeling of Single-Cell Multi-Omic Data with totalVI**（totalVI） | 单细胞 | 转录组、蛋白质组 | CITE-seq 联合潜变量、去噪、差异分析和缺失值估计 | 同细胞 RNA count 与抗体衍生标签、批次信息 | University of California, Berkeley｜2021 |
| 092. **MultiVI: Deep Generative Model for the Integration of Multimodal Data**（MultiVI） | 单细胞 | 转录组、表观组 | RNA/ATAC 联合表示、配对/非配对整合和模态补全 | paired 或 unpaired scRNA-seq/scATAC-seq、批次标签 | University of California, Berkeley｜2023 |
| 093. **Multi-Omics Single-Cell Data Integration and Regulatory Inference with Graph-Linked Embedding**（scGLUE） | 单细胞 | 转录组、表观组、多组学 | 非配对多组学整合、细胞对齐和调控网络推断 | scRNA-seq、scATAC-seq、先验基因—peak 图和批次信息 | Peking University｜2022 |
| 094. **MOFA+: A Statistical Framework for Comprehensive Integration of Multi-Modal Single-Cell Data**（MOFA+） | 单细胞 | 转录组、表观组、蛋白质组、代谢组 | 多视图潜因子、变异来源分解、亚群与特征发现 | 匹配样本/细胞的多组学矩阵及批次/分组元数据 | EMBL-EBI｜EMBL｜2020 |
| 095. **Deep Learning and Alignment of Spatially Resolved Single-Cell Transcriptomes with Tangram**（Tangram） | 单细胞 | 转录组、空间转录组、组织学图像 | 单细胞到空间映射、基因表达补全和空间细胞组成 | scRNA-seq 参考、空间表达/图像和共同基因集合 | Broad Institute｜MIT｜2021 |
| 096. **Cell2location Maps Fine-Grained Cell Types in Spatial Transcriptomics**（cell2location） | 单细胞 | 转录组、空间转录组 | 空间去卷积、细胞类型丰度和组织生态位定位 | 带注释 scRNA-seq 参考、spot 表达与空间坐标 | Wellcome Sanger Institute｜2022 |
| 097. **High-Resolution Alignment of Single-Cell and Spatial Transcriptomes with CytoSPACE**（CytoSPACE） | 单细胞 | 转录组、空间转录组 | 单细胞到空间高分辨率定位、细胞组成重建 | scRNA-seq、空间 spot 表达、坐标与共同基因 | Stanford University｜2023 |
| 098. **CellOrganizer: Learning and Using Cell Geometries for Spatial Cell Simulations**（CellOrganizer） | 单细胞 | 3D 显微图像、亚细胞结构 | 概率细胞形态建模、亚细胞结构合成与空间仿真 | 分割的 2D/3D 荧光图像、细胞/核/细胞器通道 | Carnegie Mellon University｜2019 |
| 099. **Integrated Intracellular Organization and Its Variations in Human iPS Cells** | 单细胞 | 3D 显微图像、蛋白定位、形态组 | 细胞形态空间坐标系、细胞器位置预测和状态比较 | 活细胞 3D 荧光图像、标记蛋白、细胞/核形态 | Allen Institute for Cell Science｜2023 |
| 100. **A Deep Generative Model of 3D Single-Cell Organization** | 单细胞 | 3D 显微图像、亚细胞结构 | 三维细胞与细胞器形态生成、结构关系和变异分析 | 分割的活细胞 3D 图像、细胞膜/核/细胞器通道 | Allen Institute for Cell Science｜2022 |

## 核查来源与使用建议

以下入口用于核对题名、年份、机构、版本状态与模型定位；其中综述负责交叉索引，正式论文或机构主页优先用于确认单项信息。

- 总体框架与 2025–2026 模型索引：[Virtual cells: progress and challenges from mechanistic modeling to data-driven integration](https://pmc.ncbi.nlm.nih.gov/articles/PMC13677070/)
- Arc Institute：[State](https://arcinstitute.org/news/virtual-cell-model-state)、[PIE](https://arcinstitute.org/news/pie)；Stack 及相关模型亦可由上述综述的原始文献索引核查。
- Tahoe Therapeutics：[Tahoe-100M](https://www.biorxiv.org/content/10.1101/2025.02.20.639398v3)；Tahoe-x1 使用该数据图谱并在 2025 年公开。
- Xaira Therapeutics：[X-Cell 官方说明](https://www.xaira.com/news/announcing-x-cell-our-virtual-cell-model-trained-on-billions-of-genomic-data-points-2)
- CZI 虚拟细胞模型入口：[Virtual Cell Models](https://virtualcellmodels.cziscience.com/models)
- 空间基础模型：[Novae](https://www.nature.com/articles/s41592-025-02980-0)、[NicheCompass](https://www.nature.com/articles/s41592-025-02814-z)
- 多尺度虚拟组织：[VirTues](https://www.nature.com/articles/s41586-026-10884-y)
- 机制型整细胞模型：[Bringing the Genetically Minimal Cell to Life on a Computer in 4D](https://doi.org/10.1016/j.cell.2026.02.009)

> 使用建议：如用于正式综述或项目立项，请对预印本在引用当日再次检查版本号、正式发表状态和作者单位；对模型能力的比较应同时核对训练数据是否公开、评测是否跨数据集，以及是否真正覆盖未见扰动与未见细胞背景。
