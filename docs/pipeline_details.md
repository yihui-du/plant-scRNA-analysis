# 完整分析流程详解

## 概述

本文档详细说明了植物单细胞RNA-seq数据的完整分析流程，包括从原始FASTQ文件到功能注释的所有步骤。

## 分析流程

### 第1阶段：数据准备和FASTQ处理

#### 步骤1.1：FASTQ质量控制
```bash
cd scripts/01_fastq_conversion
python fastq_conversion.py \
    --input ../../data/raw \
    --output ../../data/processed \
    --threads 8
```

**输入**：原始FASTQ文件
**输出**：质量检查报告、过滤后的FASTQ文件

#### 步骤1.2：FASTQ格式转换
- 如需要，转换为其他格式（如FASTA）
- 进行适配器移除和质量修剪

### 第2阶段：下游分析

#### 步骤2.1：数据加载和质控
```bash
cd scripts/02_downstream_analysis
python main_analysis.py \
    --config ../../configs/analysis_config.yaml \
    --input ../../data/processed/expression_matrix.csv
```

**关键步骤**：
1. 加载表达矩阵
2. 计算QC指标（基因数、UMI数、线粒体百分比）
3. 过滤低质量细胞

**输出**：
- `qc_metrics.csv` - 质量控制指标
- `qc_plots.pdf` - QC可视化图表

#### 步骤2.2：预处理和正规化
- 对数转换（log1p）
- 中心化和标准化
- 高变基因选择（默认2000个）

#### 步骤2.3：降维
- **PCA** (50成分)
- **UMAP** (用于可视化)
- **t-SNE** (可选，用于验证)

#### 步骤2.4：聚类
- 使用Leiden聚类算法
- 推荐分辨率：0.4-1.2（根据细胞类型复杂度调整）

**输出**：
- `clusters.csv` - 聚类标签
- `umap_plot.pdf` - UMAP可视化

#### 步骤2.5：细胞类型注释
- 使用marker基因进行手动注释
- 或使用数据库进行自动注释

**输出**：
- `cell_annotations.csv` - 细胞类型注释

### 第3阶段：HDWGCNA共表达分析

#### 步骤3.1：表达矩阵准备
- 选择top变异基因
- 按聚类分组（如果有多个条件）

#### 步骤3.2：软功率选择
```bash
cd scripts/03_hdwgcna
python hdwgcna_analysis.py \
    --input ../../data/processed/expression_matrix.csv \
    --config ../../configs/analysis_config.yaml
```

**软功率选择**：
- 使用无标度拓扑适应度（R² > 0.8）选择软功率
- 植物数据通常使用软功率 6-12

#### 步骤3.3：网络构建和模块检测
- 构建加权相关性网络
- 层级聚类和动态模块检测
- 模块合并（相似性> 0.25）

#### 步骤3.4：模块特征基因识别
- 计算模块特征值（module eigengene）
- 识别每个模块的中心基因（hub genes）

**输出**：
- `gene_modules.csv` - 基因模块分配
- `module_eigengenes.csv` - 模块特征值
- `module_plots.pdf` - 模块可视化

### 第4阶段：功能富集分析

#### 步骤4.1：GO富集分析
```bash
cd scripts/04_enrichment
python enrichment_analysis.py \
    --input ../../results/gene_modules.csv \
    --config ../../configs/analysis_config.yaml \
    --gene-id-type gene_symbol
```

**本体**：
- BP (Biological Process)
- CC (Cellular Component)
- MF (Molecular Function)

**输出**：
- `go_enrichment_bp.csv` - BP富集结果
- `go_enrichment_cc.csv` - CC富集结果
- `go_enrichment_mf.csv` - MF富集结果

#### 步骤4.2：KEGG富集分析

**物种代码**：
- `ath` - 拟南芥 (Arabidopsis thaliana)
- `osa` - 水稻 (Oryza sativa)
- `zma` - 玉米 (Zea mays)

**输出**：
- `kegg_enrichment.csv` - KEGG通路富集结果

#### 步骤4.3：结果整合和可视化
- 气泡图展示富集结果
- 网络图展示通路间的关系
- 热力图展示模块-功能的对应关系

## 关键参数说明

### 质控参数
| 参数 | 默认值 | 说明 |
|------|-------|------|
| min_genes | 200 | 每个细胞的最小基因数 |
| max_genes | 5000 | 每个细胞的最大基因数 |
| mt_percent | 10 | 线粒体基因比例阈值(%) |

### 聚类参数
| 参数 | 默认值 | 说明 |
|------|-------|------|
| resolution | 0.8 | Leiden聚类分辨率(越高越多clusters) |
| n_neighbors | 15 | KNN算法的邻居数 |

### WGCNA参数
| 参数 | 默认值 | 说明 |
|------|-------|------|
| soft_power | 6 | 软功率（根据无标度适应度调整） |
| min_module_size | 30 | 最小模块大小 |

## 常见问题

### Q1: 如何选择合适的软功率？
**A**: 使用sft.fit函数进行软功率测试，选择使R²>0.8且平均连接度不过低的最小软功率。

### Q2: 聚类分辨率应该设为多少？
**A**: 这取决于你期望的细胞类型细节程度。推荐尝试0.4-1.2的范围，观察聚类结果的生物学意义。

### Q3: 如何处理批次效应？
**A**: 可在预处理阶段使用harmony或Seurat的ScaleData + vars.to.regress参数。

## 输出文件说明

| 文件 | 说明 |
|------|------|
| qc_metrics.csv | 质控指标 |
| clusters.csv | 聚类结果 |
| cell_annotations.csv | 细胞类型 |
| gene_modules.csv | WGCNA模块分配 |
| enrichment_results.csv | 富集分析结果 |
| plots/ | 所有可视化图表 |

## 推荐的参考文献

1. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome biology, 19(1), 15.

2. Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 559.

3. Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284-287.

---

**最后更新**: 2026年1月14日
