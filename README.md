# 植物单细胞测序下游分析流程

## 项目简介

这个仓库包含植物单细胞RNA-seq（scRNA-seq）数据的完整下游分析流程代码，包括FASTQ文件转换、质量控制、数据处理、共表达分析和功能富集分析。

## 项目结构

```
plant-scRNA-analysis/
├── scripts/                          # 分析脚本目录
│   ├── 01_fastq_conversion/         # FASTQ格式转换脚本
│   ├── 02_downstream_analysis/      # 主要下游分析流程
│   ├── 03_hdwgcna/                  # HDWGCNA共表达分析
│   └── 04_enrichment/               # GO/KEGG富集分析
├── data/                            # 数据存储目录
│   ├── raw/                         # 原始数据
│   └── processed/                   # 处理后的数据
├── results/                         # 分析结果输出
├── configs/                         # 配置文件
├── docs/                            # 详细文档
├── issues/                          # 问题记录和解决方案
├── requirements.txt                 # Python依赖包
├── environment.yml                  # Conda环境配置
└── LICENSE                          # 许可证
```

## 各模块说明

### 1. FASTQ转换 (`01_fastq_conversion/`)
- **功能**: FASTQ格式转换、预处理和质量检查
- **主要脚本**:
  - `fastq_conversion.py` - 格式转换主程序
  - `quality_control.py` - 质量检查脚本
  
### 2. 下游分析 (`02_downstream_analysis/`)
- **功能**: 包括比对、定量、质控、聚类、注释等
- **主要脚本**:
  - `read_mapping.py` - 测序数据比对
  - `quantification.py` - 基因定量
  - `quality_metrics.py` - 质量指标计算
  - `clustering.py` - 细胞聚类分析
  - `cell_annotation.py` - 细胞类型注释
  
### 3. HDWGCNA分析 (`03_hdwgcna/`)
- **功能**: 加权基因共表达网络分析
- **主要脚本**:
  - `wgcna_preprocessing.py` - 数据预处理
  - `hdwgcna_analysis.py` - HDWGCNA主分析
  - `module_analysis.py` - 模块分析和可视化
  - `trait_correlation.py` - 性状关联分析

### 4. 功能富集分析 (`04_enrichment/`)
- **功能**: GO和KEGG通路富集分析
- **主要脚本**:
  - `go_enrichment.py` - Gene Ontology富集
  - `kegg_enrichment.py` - KEGG通路富集
  - `enrichment_visualization.py` - 富集结果可视化

## 环境要求

### Python版本
- Python 3.8+

### 主要依赖包
```
scanpy>=1.9.0          # 单细胞分析
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
scikit-learn>=1.0.0
anndata>=0.8.0
```

### R依赖（部分分析需要）
- R 4.0+
- hdwgcna
- ggplot2
- clusterProfiler

## 快速开始

### 1. 克隆仓库
```bash
git clone https://github.com/your-username/plant-scRNA-analysis.git
cd plant-scRNA-analysis
```

### 2. 创建环境
使用Conda创建环境：
```bash
conda env create -f environment.yml
conda activate plant-scRNA
```

或使用pip：
```bash
pip install -r requirements.txt
```

### 3. 准备数据
将原始FASTQ文件放在 `data/raw/` 目录下

### 4. 运行分析

**第一步：FASTQ转换**
```bash
cd scripts/01_fastq_conversion
python fastq_conversion.py --input ../../data/raw --output ../../data/processed
```

**第二步：下游分析**
```bash
cd ../02_downstream_analysis
python main_analysis.py --config ../../configs/analysis_config.yaml
```

**第三步：HDWGCNA分析**
```bash
cd ../03_hdwgcna
python hdwgcna_analysis.py --input ../../data/processed/seurat_object.h5ad
```

**第四步：富集分析**
```bash
cd ../04_enrichment
python enrichment_analysis.py --input ../../results/gene_modules.csv
```

## 输出文件

分析完成后，结果将保存在 `results/` 目录下：
- `qc_metrics.csv` - 质量控制指标
- `seurat_object.h5ad` - Seurat对象（AnnData格式）
- `clusters.csv` - 细胞聚类结果
- `gene_modules.csv` - HDWGCNA基因模块
- `enrichment_results.csv` - 富集分析结果
- `plots/` - 各类可视化图表

## 配置文件

`configs/analysis_config.yaml` 包含所有分析参数：
```yaml
# 质量控制
qc:
  min_genes: 200
  max_genes: 5000
  mt_percent: 10

# 聚类参数
clustering:
  resolution: 0.8
  n_neighbors: 15
  
# HDWGCNA参数
wgcna:
  soft_power: 6
  min_module_size: 30
```

## 常见问题

详见 `issues/` 目录中的问题记录和解决方案。

## 文档

更多详细说明请查看 `docs/` 目录：
- `pipeline_details.md` - 完整流程说明
- `parameter_explanation.md` - 参数说明
- `troubleshooting.md` - 故障排除

## 最新更新

- **2026-01-14**: 项目初始化，创建基础框架

## 引用

如果使用了本分析流程，请引用以下工具：
- Scanpy: [Wolf et al., 2018](https://doi.org/10.1186/s13059-017-1382-0)
- HDWGCNA: [Mostafavi et al., 2020](https://doi.org/10.1038/s41467-020-20241-9)

## 许可证

MIT License - 详见 LICENSE 文件

## 作者

[Your Name]

## 联系方式

如有问题或建议，请提交Issue或联系作者。

---

**最后更新**: 2026年1月14日
