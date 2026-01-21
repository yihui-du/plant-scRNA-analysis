#!/bin/bash

# 第二部分：定量分析 - cellranger count
# 用于 10x 单细胞测序数据的质量控制、比对和定量

# ============================================================
# 参数配置
# ============================================================

# 任务ID（输出目录名）
ID="SRR23825002_new_counts"

# FASTQ 文件所在目录
FASTQS="/data/duyihui/20_PRJNA941486_2025/redownload/SRR23825002"

# 参考转录组路径
TRANSCRIPTOME="/data/duyihui/reference/cellranger_arabidopsis/arabidopsis_tair10"

# 样本名称（与FASTQ文件名中的样本标签对应）
SAMPLE="SRR23825002"

# 本地计算核心数
LOCALCORES=16

# 强制细胞数阈值
FORCE_CELLS=10000

# 本地内存（GB）
LOCALMEM=64

# 是否生成 BAM 文件（true/false）
CREATE_BAM="false"

# ============================================================
# 运行 cellranger count
# ============================================================

cellranger count \
    --id="${ID}" \
    --create-bam="${CREATE_BAM}" \
    --fastqs="${FASTQS}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --sample="${SAMPLE}" \
    --localcores="${LOCALCORES}" \
    --force-cells="${FORCE_CELLS}" \
    --localmem="${LOCALMEM}"

# ============================================================
# 输出说明
# ============================================================
# 
# 运行完成后，在 ${ID} 目录中会生成以下文件：
# 
# 主要输出：
#   - outs/filtered_feature_bc_matrix/     细胞×基因表达矩阵（HTSeq 格式）
#   - outs/raw_feature_bc_matrix/          原始计数矩阵
#   - outs/molecule_info.h5                分子级别统计信息
#   - outs/metrics_summary.csv             定量统计汇总
#   - outs/web_summary.html                交互式质控报告
# 
# 可选输出（如果 CREATE_BAM=true）：
#   - outs/possorted_genome_bam.bam        排序后的 BAM 文件
#   - outs/possorted_genome_bam.bam.bai    BAM 索引文件
# 
# 后续分析通常使用 filtered_feature_bc_matrix 或 molecule_info.h5
#
