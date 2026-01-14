# 🚀 开始使用指南 - 三步快速上手

## ⏱️ 总耗时：15分钟

---

## 第1步：了解项目（5分钟）

### 打开这三个文件（按顺序）：

#### 1️⃣ `QUICK_REFERENCE.md` ⭐
```
快速参考卡片，包含所有常用命令和文档导航
位置: D:\test\plant-scRNA-analysis\QUICK_REFERENCE.md
```

#### 2️⃣ `README.md`
```
项目总览、环境要求、快速开始
位置: D:\test\plant-scRNA-analysis\README.md
```

#### 3️⃣ `configs/analysis_config.yaml`
```
所有参数配置，包含详细的中文注释
位置: D:\test\plant-scRNA-analysis\configs\analysis_config.yaml
```

---

## 第2步：复制你的代码（5分钟）

### 将你的代码文件按以下结构复制：

```
你的代码位置                              →  目标位置
───────────────────────────────────────────────────────────

fastq_conversion.py                      →  scripts/01_fastq_conversion/
your_fastq_code.py                       →  scripts/01_fastq_conversion/

downstream_analysis.py                   →  scripts/02_downstream_analysis/
clustering.py                            →  scripts/02_downstream_analysis/

hdwgcna.py                               →  scripts/03_hdwgcna/
module_analysis.py                       →  scripts/03_hdwgcna/

enrichment.py                            →  scripts/04_enrichment/
go_kegg.py                               →  scripts/04_enrichment/
```

### 检查清单：
- [ ] FASTQ转换代码已复制
- [ ] 下游分析代码已复制
- [ ] HDWGCNA代码已复制
- [ ] 富集分析代码已复制
- [ ] 更新了 `issues/CODE_UPLOAD_CHECKLIST.md`

---

## 第3步：安装环境（5分钟）

### 选择方式1：使用Conda（推荐）
```bash
conda env create -f environment.yml
conda activate plant-scRNA
```

### 或选择方式2：使用pip
```bash
pip install -r requirements.txt
```

### 验证安装
```bash
python -c "import scanpy; print(f'Scanpy版本: {scanpy.__version__}')"
```

---

## 现在你可以做什么

### ✅ 运行分析
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

### ✅ 调整参数
编辑 `configs/analysis_config.yaml` 中的参数，无需修改代码

### ✅ 记录问题
遇到问题时，在 `docs/troubleshooting.md` 中查找解决方案或记录新问题

### ✅ 上传到GitHub
参考 `GIT_UPLOAD_GUIDE.md` 中的完整说明

---

## 📚 根据需要查看的文档

| 需求 | 查看文件 | 耗时 |
|------|--------|------|
| 快速参考所有命令 | `QUICK_REFERENCE.md` | 5分钟 |
| 理解项目结构 | `README.md` | 10分钟 |
| 调整分析参数 | `configs/analysis_config.yaml` | 根据需要 |
| 理解每个参数 | `docs/parameter_explanation.md` | 查询时间 |
| 了解完整流程 | `docs/pipeline_details.md` | 15分钟 |
| 遇到问题 | `docs/troubleshooting.md` | 查询时间 |
| 上传到GitHub | `GIT_UPLOAD_GUIDE.md` | 30分钟 |

---

## 🎯 常见场景应对

### 场景1：我想快速测试分析流程
```bash
1. 放一些测试数据到 data/raw/
2. 编辑 configs/analysis_config.yaml（如需要）
3. 运行各个脚本
4. 查看 results/ 下的输出
```

### 场景2：我想调整分析参数
```bash
1. 打开 configs/analysis_config.yaml
2. 找到需要调整的参数
3. 参考 docs/parameter_explanation.md 理解含义
4. 修改参数值
5. 重新运行脚本（代码会自动读取新参数）
```

### 场景3：我遇到了问题
```bash
1. 查看 docs/troubleshooting.md
2. 找到类似问题和解决方案
3. 如果没有找到，记录在 issues/CHANGELOG.md
4. 解决后更新 docs/troubleshooting.md
```

### 场景4：我想上传代码到GitHub
```bash
1. 阅读 GIT_UPLOAD_GUIDE.md
2. 安装Git（如未安装）
3. 创建GitHub账户和仓库
4. 按步骤推送代码
5. 后续更新时定期 git push
```

---

## 💡 小贴士

- **修改参数最简单的方式**: 编辑 `configs/analysis_config.yaml`
- **快速查找命令**: 查阅 `QUICK_REFERENCE.md`
- **理解流程**: 按顺序读 `README.md` → `docs/pipeline_details.md`
- **记录工作**: 使用 `issues/CHANGELOG.md` 跟踪进度
- **备份代码**: 定期 `git push` 到GitHub

---

## ❓ 常见问题

### Q: 如何修改分析参数？
A: 编辑 `configs/analysis_config.yaml` 文件，所有参数都在这里

### Q: 代码该放在哪里？
A: 按照四个阶段放在 `scripts/01_*` 到 `scripts/04_*` 目录下

### Q: 数据该放在哪里？
A: 原始数据放 `data/raw/`，处理后的放 `data/processed/`

### Q: 结果会存到哪里？
A: 所有结果都会保存到 `results/` 目录

### Q: 如何查看详细的参数说明？
A: 查看 `docs/parameter_explanation.md` 或 `configs/analysis_config.yaml` 中的注释

### Q: 如何上传到GitHub？
A: 参考 `GIT_UPLOAD_GUIDE.md` 中的完整说明

---

## 🔗 快速链接

| 文件 | 路径 |
|-----|------|
| 快速参考 | `QUICK_REFERENCE.md` |
| 项目总结 | `PROJECT_SUMMARY.md` |
| 这个文件 | `START_HERE.md` |
| README | `README.md` |
| 配置文件 | `configs/analysis_config.yaml` |
| 完整流程 | `docs/pipeline_details.md` |
| 参数说明 | `docs/parameter_explanation.md` |
| 故障排除 | `docs/troubleshooting.md` |
| Git上传 | `GIT_UPLOAD_GUIDE.md` |

---

## ✨ 下一步

现在你已经准备好了！选择一个开始：

1. **立即开始**: 复制代码到 `scripts/` 目录
2. **深入学习**: 阅读 `docs/pipeline_details.md`
3. **快速查找**: 使用 `QUICK_REFERENCE.md`
4. **准备上传**: 按照 `GIT_UPLOAD_GUIDE.md`

---

**🎉 祝你使用愉快！**

有任何问题，可以查阅相关文档。

**最后更新**: 2026年1月14日

