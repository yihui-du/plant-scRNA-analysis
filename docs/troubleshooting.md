# 常见问题和故障排除指南

## 数据加载问题

### 问题1：无法读取FASTQ文件

**错误信息**：`Error reading FASTQ file`

**可能原因**：
1. 文件格式不正确
2. 文件编码问题
3. 文件被损坏

**解决方案**：
```bash
# 检查文件格式
file your_file.fastq

# 查看文件头部
head -n 20 your_file.fastq

# 使用FastQC进行质量检查
fastqc your_file.fastq
```

### 问题2：表达矩阵加载失败

**错误信息**：`ValueError: Could not load data`

**可能原因**：
1. 行列标注不正确
2. 数据格式不一致
3. 文件编码问题

**解决方案**：
```python
import pandas as pd

# 检查文件格式
df = pd.read_csv('expression_matrix.csv', index_col=0, nrows=5)
print(df.shape)
print(df.index)
print(df.columns)
```

---

## 质控问题

### 问题3：过滤后细胞数量过多

**现象**：设定的质控参数很严格，但仍然过滤不掉足够的细胞

**可能原因**：
1. 参数设置不合理
2. 数据本身质量较好

**解决方案**：
```python
# 查看QC指标的分布
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.hist(adata.obs['n_genes_by_counts'], bins=50)
plt.xlabel('Number of genes')

plt.subplot(1, 3, 2)
plt.hist(adata.obs['total_counts'], bins=50)
plt.xlabel('Total UMI')

plt.subplot(1, 3, 3)
plt.hist(adata.obs['pct_counts_mt'], bins=50)
plt.xlabel('Mitochondrial %')
plt.tight_layout()
plt.savefig('qc_distributions.pdf')
```

### 问题4：过滤后细胞数量过少

**现象**：过滤后剩余细胞太少（<100个）

**可能原因**：
1. 原始数据质量差
2. 参数过于严格

**解决方案**：
```yaml
# 调整配置文件中的参数
qc:
  min_genes: 100       # 降低最小基因数
  max_genes: 8000      # 提高最大基因数
  mt_percent_threshold: 15  # 放宽线粒体百分比阈值
```

---

## 聚类问题

### 问题5：聚类结果看起来不合理

**现象**：Cluster数量过多或过少，或者某个cluster太大

**可能原因**：
1. 分辨率参数不合适
2. 降维不充分
3. 数据质量问题

**解决方案**：
```python
# 尝试不同的分辨率
import scanpy as sc

for res in [0.4, 0.6, 0.8, 1.0, 1.2]:
    sc.tl.leiden(adata, resolution=res)
    print(f"Resolution {res}: {adata.obs['leiden'].nunique()} clusters")

# 观察每个分辨率下的聚类质量
```

### 问题6：UMAP图上的cluster分布不均匀

**现象**：某些cluster聚集在一起，难以区分

**可能原因**：
1. 批次效应存在
2. 细胞类型间差异不明显

**解决方案**：
```python
# 检查是否存在批次效应
sc.pl.umap(adata, color='batch', show=False)
plt.savefig('batch_effect.pdf')

# 使用harmony进行批次校正
import harmony
adata_corrected = harmony.integrate(adata, 'batch')
```

---

## HDWGCNA问题

### 问题7：软功率选择困难

**现象**：无法找到合适的软功率（R² > 0.8）

**可能原因**：
1. 基因选择不当
2. 样本量太少

**解决方案**：
```python
# 减少基因数量
expr_matrix = expr_matrix.iloc[:, :2000]  # 选择top 2000基因

# 或尝试更宽泛的软功率范围
soft_powers = range(1, 31)
```

### 问题8：模块数量过多或过少

**现象**：检测到的模块数量不符合预期

**可能原因**：
1. 最小模块大小设置不当
2. 模块合并阈值不合理

**解决方案**：
```yaml
wgcna:
  min_module_size: 20    # 降低可检测更多小模块
  merge_threshold: 0.3   # 降低可合并更少模块
```

---

## 富集分析问题

### 问题9：GO/KEGG富集结果为空

**现象**：没有显著的富集项

**可能原因**：
1. 基因ID格式不正确
2. 基因数量太少
3. 物种设置错误

**解决方案**：
```python
# 检查基因ID格式
print(genes[:10])

# 确保使用正确的ID类型
# 常见格式：ENSEMBL ID、基因名称、Entrez ID

# 检查基因数量
print(f"Total genes: {len(genes)}")

# 确保使用正确的物种代码
# 拟南芥: ath
# 水稻: osa
# 玉米: zma
```

### 问题10：富集分析运行缓慢

**现象**：富集分析耗时过长

**可能原因**：
1. 基因数量太多
2. 注释数据库太大

**解决方案**：
```python
# 在进行富集分析前选择高质量的基因
# 例如，只分析显著差异的基因

genes_filtered = genes[genes['padj'] < 0.05]

# 或增加p值阈值以减少输入基因数
```

---

## 内存和计算问题

### 问题11：内存溢出错误

**错误信息**：`MemoryError`

**可能原因**：
1. 数据集太大
2. 参数设置导致中间文件过大

**解决方案**：
```python
# 限制基因数量
adata = adata[:, :adata.n_vars // 2]

# 使用稀疏矩阵
import scipy.sparse as sp
adata.X = sp.csr_matrix(adata.X)

# 分批处理
# 可将大数据集分成多个小批次处理
```

### 问题12：计算超时

**现象**：某个步骤执行时间过长

**可能原因**：
1. 数据规模太大
2. 参数不合理

**解决方案**：
```python
# 使用多线程处理
import multiprocessing
n_jobs = -1  # 使用所有CPU

# 降低邻居数
sc.pp.neighbors(adata, n_neighbors=10)

# 减少计算精度
# 某些步骤可使用更粗糙的参数加快计算
```

---

## 文件导出问题

### 问题13：无法保存大型h5ad文件

**错误信息**：`Error writing file`

**可能原因**：
1. 磁盘空间不足
2. 权限问题
3. 文件系统限制

**解决方案**：
```bash
# 检查磁盘空间
df -h

# 使用分块保存
python -c "
import scanpy as sc
adata = sc.read_h5ad('large_file.h5ad', backed='r')
# 分块处理和保存
"
```

---

## 可视化问题

### 问题14：图表无法正常显示或保存

**可能原因**：
1. 后端问题
2. 字体问题
3. 参数不合理

**解决方案**：
```python
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端

# 检查可用的绘图函数
import scanpy as sc
sc.settings.figdir = './results/plots/'
sc.set_figure_params(dpi=100, vector_friendly=True)
```

---

## 报告bug

如果遇到未列出的问题，请提交Issue：

1. 提供完整的错误信息
2. 提供重现问题的最小代码例子
3. 提供环境信息（Python版本、包版本等）

```bash
# 获取环境信息
python -c "
import sys, scanpy, numpy, pandas
print(f'Python: {sys.version}')
print(f'Scanpy: {scanpy.__version__}')
print(f'NumPy: {numpy.__version__}')
print(f'Pandas: {pandas.__version__}')
"
```

---

**最后更新**: 2026年1月14日
