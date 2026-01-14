# README模板 - 更新说明

## 如何使用本仓库

### 第一次使用

1. **克隆仓库**
   ```bash
   git clone https://github.com/your-username/plant-scRNA-analysis.git
   cd plant-scRNA-analysis
   ```

2. **安装依赖**
   ```bash
   # 使用Conda
   conda env create -f environment.yml
   conda activate plant-scRNA
   
   # 或使用pip
   pip install -r requirements.txt
   ```

3. **准备数据**
   - 将原始FASTQ文件放在 `data/raw/` 目录下
   - 或将处理后的表达矩阵放在 `data/processed/` 目录下

### 运行分析

按照 `docs/pipeline_details.md` 中的步骤运行分析。

### 查看结果

分析完成后，结果将保存在 `results/` 目录下：
- `plots/` - 所有可视化图表
- `*.csv` - 分析结果表格

### 遇到问题

查阅 `docs/troubleshooting.md` 或提交Issue。

---

## 贡献指南

如果你有新的分析脚本或改进建议，欢迎提交Pull Request！

### 文件结构约定

- 脚本文件放在 `scripts/` 下的相应目录
- 每个脚本应包含详细的文档字符串
- 使用 `loguru` 库记录日志

### 代码规范

- 使用PEP 8风格
- 变量和函数使用有意义的英文名称
- 添加类型提示
- 包含测试代码

---

## 许可证

MIT License - 详见 LICENSE 文件

---

**最后更新**: 2026年1月14日
