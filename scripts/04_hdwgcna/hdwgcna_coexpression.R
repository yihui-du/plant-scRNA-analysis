## ================== 基本参数手动设置 ==================
file_prex <- "Atha_leaf_PRJNA648028"  
work_wd   <- "/data/duyihui/scRNA"  # 你的 rds 所在目录

## ================== 加载包 ==================
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 20)

## ================== 读入数据 ==================
setwd(work_wd)
seurat_obj <- readRDS(paste0(file_prex, "_allgene_rename.rds"))

## ================== 【新增 】输出细胞类型/数量 + 基因数量 ==================
cat("总细胞数:", ncol(seurat_obj), "\n")
# 按 RNA assay 统计总基因数（feature 数）
cat("总基因数(RNA features):", nrow(seurat_obj[["RNA"]]), "\n")
# 这里定义 genes_all
genes_all <- rownames(seurat_obj[["RNA"]])
cat("总基因数(RNA features):", length(genes_all), "\n\n")
cat("细胞类型与数量：\n")
print(sort(table(seurat_obj$cell_type), decreasing = TRUE))

is_chloro_all <- grepl("^ATCG", genes_all)           # 叶绿体
is_mito_all   <- grepl("^ATMG", genes_all)           # 线粒体
is_nuclear_all <- grepl("^AT[1-5]G", genes_all)      # 核基因 AT1G-AT5G
cat("\n===== 全样本基因分类（按基因名） =====\n")
cat("叶绿体基因数 (ATCG*):", sum(is_chloro_all), "\n")
cat("线粒体基因数 (ATMG*):", sum(is_mito_all), "\n")
cat("核基因数 (AT1G–AT5G):", sum(is_nuclear_all), "\n")

types    <- "Vasculature"  # 想分析的细胞类型：比如 "Mesophyll" / "Epidermis" / "Vasculature"
## ================== 设置初始输出目录（占位，稍后会按匹配到的真实类型覆盖） ==================
file_types <- gsub(" ", "_", types)
out_dir <- file.path(work_wd, paste0(file_prex, "_", file_types, "_hdwgcna"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
cat("工作目录：", getwd(), "\n")

## ================== 1. SetupForWGCNA（按原脚本写法） ==================
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.01,  # 基因至少在 1% 细胞中有表达
  wgcna_name = paste0(file_prex, "_", file_types)
)

## ================== 2. 构建 Metacells ==================
if (length(unique(seurat_obj@meta.data$orig.ident)) > 1 ){
  seurat_obj <- MetacellsByGroups(
    seurat_obj   = seurat_obj,
    group.by     = c("cell_type","orig.ident"),
    reduction    = "harmony",
    k            = 25,
    max_shared   = 10,
    ident.group  = "cell_type"
  )
} else {
  seurat_obj <- MetacellsByGroups(
    seurat_obj   = seurat_obj,
    group.by     = c("cell_type","orig.ident"),
    reduction    = "pca",
    k            = 25,
    max_shared   = 10,
    ident.group  = "cell_type"
  )
}

# 归一化 metacell 表达矩阵
seurat_obj <- NormalizeMetacells(seurat_obj)

## 3.1 看看 metacell 后的 cell_type 列表
celltypes_chr <- as.character(seurat_obj$cell_type)
ct_levels     <- sort(unique(celltypes_chr))

cat("\n===== Metacell 之后的 cell_type levels =====\n")
print(ct_levels)
cat("============================================\n\n")

## 3.2 逻辑上想选的名称关键字（你手动改这里）
ct_keyword <- "Vasculature"   # 这里你写 “Epidermis” / “Atrichoblast” 等关键字

## 用 grep 在 level 里模糊匹配
match_ct <- grep(ct_keyword, ct_levels, value = TRUE)

if (length(match_ct) == 0) {
  stop("在 cell_type levels 中找不到包含 '", ct_keyword, "' 的条目，请检查关键字。")
}
if (length(match_ct) > 1) {
  cat("警告：关键字 '", ct_keyword, "' 命中了多个 cell_type：\n",
      paste(match_ct, collapse = " | "), "\n",
      "先使用第一个：", match_ct[1], "\n\n", sep = "")
}

## 这个才是真正存在于对象里的 cell_type 名称（可能带奇怪空格）
types <- match_ct[1]
file_types <- gsub(" ", "_", types)

## 根据最终 types 重新设置输出目录，并切换工作目录
out_dir <- file.path(work_wd, paste0(file_prex, "_", file_types, "_hdwgcna"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
cat("最终工作目录：", getwd(), "\n")

cat("选中的真实 types：'", types, "'\n", sep = "")
## 3 为指定细胞类型设置 WGCNA 输入矩阵（只调用一次）
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name    = types,
  group.by      = "cell_type",
  assay         = "RNA",
  slot          = "data",
  return_seurat = TRUE
)
## 从对象里正规地取出 datExpr
datExpr <- GetDatExpr(seurat_obj)

cat("\n===== ", types, " 进入 WGCNA 的矩阵信息 =====\n", sep = "")
cat("metacell 数:", nrow(datExpr), "\n")
cat("用于 WGCNA 的基因数:", ncol(datExpr), "\n")

# 选中细胞类型中参与 WGCNA 的基因分类
genes_ct <- colnames(datExpr)
is_chloro_ct  <- grepl("^ATCG", genes_ct)
is_mito_ct    <- grepl("^ATMG", genes_ct)
is_nuclear_ct <- grepl("^AT[1-5]G", genes_ct)

cat("\n--- ", types, " 中参与 WGCNA 的基因分类 ---\n", sep = "")
cat("叶绿体基因数 (ATCG*):", sum(is_chloro_ct), "\n")
cat("线粒体基因数 (ATMG*):", sum(is_mito_ct), "\n")
cat("核基因数 (AT1G–AT5G):", sum(is_nuclear_ct), "\n")

# 保存 datExpr（可选）
saveRDS(datExpr,
        file = file.path(out_dir, paste0(file_prex, "_", file_types, "_datExpr.rds")))

# 计算并保存基因-基因相关矩阵（后期建子网用）
cor_mat <- cor(datExpr, use = "pairwise.complete.obs", method = "pearson")
saveRDS(cor_mat,
        file = file.path(out_dir, paste0(file_prex, "_", file_types, "_cor_mat.rds")))

cat("相关系数矩阵已保存：",
    file.path(out_dir, paste0(file_prex, "_", file_types, "_cor_mat.rds")), "\n")

## ================== 4. 软阈值测试 ==================
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "unsigned"
)

# 画 SoftPower 图保存
plot_list <- PlotSoftPowers(seurat_obj)
pdf(file.path(out_dir, paste0(file_prex, "_", file_types, "_01_SoftPower.pdf")),
    width = 10, height = 8)
wrap_plots(plot_list, ncol = 2)
dev.off()

# 也导出 power table
power_table <- GetPowerTable(seurat_obj)
write.csv(power_table,
          file = file.path(out_dir, paste0(file_prex, "_", file_types, "_SoftPower_table.csv")),
          row.names = FALSE)

cat("SoftPower 结果已输出，记得人工选一个合适的 power。\n")

## ================== 【新增 C】人工填写软阈值的位置 ==================
# 你打开 01_SoftPower.pdf 后，把这里改成你要用的阈值
soft_power <- 6   # <<< 你看图后改这里

## ================== 5. 构建共表达网络 ==================
# 为兼容不同 hdWGCNA 版本（有的参数叫 soft_power，有的叫 power）
construct_formals <- names(formals(hdWGCNA::ConstructNetwork))
construct_args <- list(
  seurat_obj = seurat_obj,
  setDatExpr = FALSE,
  tom_name   = paste0(file_types, "_tom")
)
if ("soft_power" %in% construct_formals) construct_args$soft_power <- soft_power
if ("power" %in% construct_formals) construct_args$power <- soft_power

seurat_obj <- do.call(hdWGCNA::ConstructNetwork, construct_args)

# 保存模块树状图
pdf(file.path(out_dir, paste0(file_prex, "_", file_types, "_02_Dendrogram.pdf")),
    width = 10, height = 6)
PlotDendrogram(seurat_obj, main = paste0(file_types, " hdWGCNA Dendrogram"))
dev.off()

## ================== 6. 导出 TOM（可选，文件会很大） ==================
TOM <- GetTOM(seurat_obj)
write.csv(TOM,
          file = file.path(out_dir, paste0(file_prex, "_", file_types, "_TOM.csv")),
          row.names = FALSE)

## ================== 7. 模块特征值（ME） & 批次矫正 ==================
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

if (length(unique(seurat_obj@meta.data$orig.ident)) > 1 ) {
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = "orig.ident"
  )
} else {
  seurat_obj <- ModuleEigengenes(seurat_obj)
}

hMEs <- GetMEs(seurat_obj, harmonized = TRUE)

## ================== 8. 计算 kME（模块连接度） ==================
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by   = "cell_type",
  group_name = types
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(file_types, "-M")
)

## ================== 9. 导出模块表 & 保存对象 ==================
modules <- GetModules(seurat_obj)
write.csv(modules,
          file = file.path(out_dir, paste0(file_prex, "_", file_types, "_modules.csv")),
          row.names = FALSE)

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

saveRDS(seurat_obj,
        file = file.path(out_dir, paste0(file_prex, "_", file_types, "_hdwgcna.rds")))

## ================== 10. 画 DotPlot ==================
p <- DotPlot(seurat_obj, features = mods, group.by = "cell_type") +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue") +
  ggtitle(types) +
  theme(
    plot.title  = element_text(hjust = 0.5, size = rel(1.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )

pdf(file.path(out_dir, paste0(file_prex, "_", file_types, "_03_Module_DotPlot.pdf")),
    width = 18, height = 10)
print(p)
dev.off()

cat("==== hdWGCNA 完成，所有结果已写入：", out_dir, " ====\n")
