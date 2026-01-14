# 本脚本用于记录已上传的代码文件

## 已上传代码

### 第1阶段：FASTQ转换
- [ ] `scripts/01_fastq_conversion/fastq_conversion.py` - 格式转换主程序
- [ ] `scripts/01_fastq_conversion/quality_control.py` - 质量检查
- [ ] `scripts/01_fastq_conversion/adapter_removal.py` - 适配器移除
- [ ] `scripts/01_fastq_conversion/utils.py` - 辅助函数

### 第2阶段：下游分析
- [ ] `scripts/02_downstream_analysis/main_analysis.py` - 主程序
- [ ] `scripts/02_downstream_analysis/quality_metrics.py` - 质控指标
- [ ] `scripts/02_downstream_analysis/preprocessing.py` - 预处理
- [ ] `scripts/02_downstream_analysis/clustering.py` - 聚类
- [ ] `scripts/02_downstream_analysis/visualization.py` - 可视化

### 第3阶段：HDWGCNA分析
- [ ] `scripts/03_hdwgcna/hdwgcna_analysis.py` - 主分析程序
- [ ] `scripts/03_hdwgcna/module_analysis.py` - 模块分析
- [ ] `scripts/03_hdwgcna/trait_correlation.py` - 性状关联
- [ ] `scripts/03_hdwgcna/visualization.py` - 可视化

### 第4阶段：富集分析
- [ ] `scripts/04_enrichment/go_enrichment.py` - GO富集
- [ ] `scripts/04_enrichment/kegg_enrichment.py` - KEGG富集
- [ ] `scripts/04_enrichment/enrichment_visualization.py` - 可视化

### 辅助脚本
- [ ] `scripts/utils/data_loader.py` - 数据加载器
- [ ] `scripts/utils/config_manager.py` - 配置管理
- [ ] `scripts/utils/logger_setup.py` - 日志设置

---

## 上传流程

1. 复制代码文件到对应目录
2. 检查代码是否包含敏感信息
3. 更新此文件中的复选框
4. 提交并推送到GitHub

---

## 文件上传检查清单

每个文件上传前请检查：
- [ ] 代码已注释和文档化
- [ ] 移除了绝对路径（使用相对路径）
- [ ] 移除了调试代码
- [ ] 包含了适当的错误处理
- [ ] 使用了日志记录
- [ ] 代码风格符合PEP 8

---

**最后更新**: 2026年1月14日
