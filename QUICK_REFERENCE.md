# 🎯 快速参考卡片

## 项目目录位置
```
D:\test\plant-scRNA-analysis\
```

---

## 📚 文档导航（按使用频率排序）

### 🔴 必读（第一次使用）
1. **项目总结** → `PROJECT_SUMMARY.md`
   - 了解项目全貌和下一步操作

2. **README** → `README.md`
   - 项目介绍、快速开始、环境要求

### 🟠 常用（日常工作）
3. **配置文件** → `configs/analysis_config.yaml`
   - 调整分析参数

4. **参数说明** → `docs/parameter_explanation.md`
   - 理解每个参数的含义

5. **完整流程** → `docs/pipeline_details.md`
   - 了解分析流程的每一步

### 🟡 备用（遇到问题）
6. **故障排除** → `docs/troubleshooting.md`
   - 常见问题解决方案

7. **更新日志** → `issues/CHANGELOG.md`
   - 记录问题和改进

### 🟢 高级（上传代码）
8. **Git上传指南** → `GIT_UPLOAD_GUIDE.md`
   - 完整的Git/GitHub使用说明

9. **上传检查** → `issues/CODE_UPLOAD_CHECKLIST.md`
   - 上传前的检查清单

---

## 🚀 快速命令

### 环境设置
```bash
# 使用Conda
conda env create -f environment.yml
conda activate plant-scRNA

# 或使用pip
pip install -r requirements.txt
```

### 运行分析
```bash
# FASTQ转换
cd scripts/01_fastq_conversion
python fastq_conversion.py --input ../../data/raw --output ../../data/processed

# 下游分析
cd ../02_downstream_analysis
python main_analysis.py --config ../../configs/analysis_config.yaml

# HDWGCNA
cd ../03_hdwgcna
python hdwgcna_analysis.py --input ../../data/processed/expression_matrix.csv

# 富集分析
cd ../04_enrichment
python enrichment_analysis.py --input ../../results/gene_modules.csv
```

### Git操作
```bash
# 初始化
git init
git add .
git commit -m "Initial commit"
git remote add origin https://github.com/YOUR-USERNAME/plant-scRNA-analysis.git

# 推送
git branch -M main
git push -u origin main

# 更新
git status
git add .
git commit -m "Update: [描述]"
git push
```

---

## 📂 文件结构速查

```
plant-scRNA-analysis/
├── 📘 README.md                       ← 项目说明
├── 📋 PROJECT_SUMMARY.md              ← 项目总结 ⭐
├── 📚 GIT_UPLOAD_GUIDE.md            ← Git上传指南
├── ⚙️ configs/analysis_config.yaml   ← 参数配置文件 ⭐
│
├── 📁 scripts/                        ← 分析脚本
│   ├── 01_fastq_conversion/          ← FASTQ处理
│   ├── 02_downstream_analysis/       ← 主分析
│   ├── 03_hdwgcna/                   ← 共表达
│   └── 04_enrichment/                ← 富集分析
│
├── 📁 data/                           ← 数据
│   ├── raw/                          ← 原始数据
│   └── processed/                    ← 处理数据
│
├── 📁 results/                        ← 结果输出
│
├── 📁 docs/                           ← 文档库 ⭐
│   ├── pipeline_details.md           ← 流程说明
│   ├── parameter_explanation.md      ← 参数详解
│   └── troubleshooting.md            ← 常见问题
│
└── 📁 issues/                         ← 问题记录
    └── CHANGELOG.md                  ← 更新日志
```

---

## 🎓 学习路径

### 第1天：了解项目
- [ ] 读 `PROJECT_SUMMARY.md`
- [ ] 读 `README.md`
- [ ] 查看项目文件结构

### 第2天：准备环境
- [ ] 复制代码文件到 `scripts/` 目录
- [ ] 运行 `pip install -r requirements.txt`
- [ ] 测试导入必要的库

### 第3天：理解流程
- [ ] 读 `docs/pipeline_details.md`
- [ ] 查看 `configs/analysis_config.yaml`
- [ ] 阅读 `docs/parameter_explanation.md`

### 第4天+：开始工作
- [ ] 准备数据到 `data/raw/`
- [ ] 根据需要修改 `configs/analysis_config.yaml`
- [ ] 运行各个分析脚本
- [ ] 保存结果到 `results/`

### 遇到问题时
- [ ] 查 `docs/troubleshooting.md`
- [ ] 在 `issues/CHANGELOG.md` 中记录
- [ ] 更新问题和解决方案文档

---

## 💾 关键操作

### 将代码复制到项目
```bash
# 假设你的代码在其他地方
# 复制到对应的scripts目录

# 例如：FASTQ转换脚本
copy your_fastq_script.py scripts\01_fastq_conversion\

# 更新检查清单
编辑 issues\CODE_UPLOAD_CHECKLIST.md
```

### 记录问题和解决方案
```bash
# 编辑更新日志
编辑 issues\CHANGELOG.md

# 编辑故障排除指南
编辑 docs\troubleshooting.md

# 提交更新
git add .
git commit -m "Update: Record issue and solution"
git push
```

### 调整分析参数
```bash
# 编辑配置文件
编辑 configs\analysis_config.yaml

# 运行分析时自动使用新参数
python scripts/02_downstream_analysis/main_analysis.py \
    --config configs/analysis_config.yaml
```

---

## ✨ 特色功能

| 功能 | 位置 | 说明 |
|------|------|------|
| **统一配置管理** | `configs/analysis_config.yaml` | 所有参数集中管理 |
| **详细参数说明** | `docs/parameter_explanation.md` | 每个参数的含义和建议 |
| **完整流程文档** | `docs/pipeline_details.md` | 从FASTQ到富集的全过程 |
| **故障排除指南** | `docs/troubleshooting.md` | 10+常见问题的解决方案 |
| **问题记录系统** | `issues/CHANGELOG.md` | 跟踪项目更新和问题 |
| **Git上传指南** | `GIT_UPLOAD_GUIDE.md` | 完整的GitHub上传说明 |
| **代码模板** | `scripts/*/` | 四个分析阶段的代码模板 |
| **环境管理** | `requirements.txt` / `environment.yml` | 轻松创建分析环境 |

---

## 🎯 下一步行动

### 立即做
1. ✅ 读 `PROJECT_SUMMARY.md`（5分钟）
2. ✅ 查看项目文件结构（2分钟）
3. ✅ 复制你的代码到 `scripts/` 目录（10分钟）

### 稍后做
4. 设置Git仓库（参考 `GIT_UPLOAD_GUIDE.md`）
5. 创建GitHub仓库
6. 推送到GitHub

### 持续做
7. 更新 `CHANGELOG.md` 记录进展
8. 在 `troubleshooting.md` 中记录问题
9. 定期 `git push` 保存工作

---

## 💡 小提示

- **快速查找**：使用Ctrl+F搜索关键词
- **参数调试**：先在配置文件中尝试不同参数
- **备份工作**：定期 `git push` 到GitHub
- **寻求帮助**：参考故障排除文档，记录在CHANGELOG中

---

## 📞 紧急快速命令

```bash
# 需要帮助？
cat README.md                              # 项目总览
cat docs/troubleshooting.md               # 常见问题
cat GIT_UPLOAD_GUIDE.md                   # Git帮助

# 快速测试
python -c "import scanpy; print(scanpy.__version__)"

# 查看环境
pip list | grep -E "scanpy|numpy|pandas"

# 列出文件
dir /b /s scripts\                        # 查看脚本
dir /b data\                              # 查看数据
```

---

**⏰ 使用此卡片节省时间！**

**📅 最后更新**: 2026年1月14日

