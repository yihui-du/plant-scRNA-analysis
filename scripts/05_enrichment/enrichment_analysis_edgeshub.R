############################################################
## Scheme 2: Custom hub genes from Cytoscape edge table
## Input: hdWGCNA rds + per-module edges csv (Gene1/Gene2)
## Output: GO(BP/CC/MF) + KEGG dotplot PDFs (no overwrite Scheme1)
############################################################

############################################################
## 0) 环境准备：加载包（假设你已装好）
############################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.At.tair.db)

############################################################
## 1) 读入 hdWGCNA rds：抽出 active_wgcna，构建 universe
############################################################

# >>>>>> 你只需要改这里 <<<<<<
rds_file <- "Atha_leaf_PRJNA648028_Mesophyll_hdwgcna.rds"
# >>>>>> 改完就行 <<<<<<

obj <- readRDS(rds_file)

stopifnot("active_wgcna" %in% names(obj@misc))
active_name <- obj@misc$active_wgcna
w <- obj@misc[[active_name]]

mods <- w$wgcna_modules
deg  <- w$wgcna_degrees


mods$gene_name <- as.character(mods$gene_name)
mods$module    <- as.character(mods$module)

mods_ng  <- subset(mods, module != "grey")
universe <- unique(mods_ng$gene_name)

cat("Universe size:", length(universe), "\n")
cat("Modules:\n"); print(table(mods_ng$module))

############################################################
## 2) 富集函数：GO + KEGG（TAIR ID 直接跑）
############################################################

run_go <- function(gene, universe, ont = "BP",
                   pAdjustMethod="BH", qvalueCutoff=0.05, min_n=10) {
  gene <- unique(gene)
  gene <- gene[gene %in% universe]
  if (length(gene) < min_n) return(NULL)

  eg <- enrichGO(
    gene          = gene,
    universe      = universe,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = ont,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff  = qvalueCutoff,
    readable      = TRUE
  )
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
  eg
}

run_kegg_tair <- function(gene, universe, organism="ath",
                          pAdjustMethod="BH", qvalueCutoff=0.05, min_n=10) {
  gene <- unique(gene)
  gene <- gene[gene %in% universe]
  if (length(gene) < min_n) return(NULL)

  ek <- enrichKEGG(
    gene          = gene,
    universe      = universe,
    organism      = organism,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff  = qvalueCutoff
  )
  if (is.null(ek) || nrow(as.data.frame(ek)) == 0) return(NULL)
  ek
}

############################################################
## 3) 读入“每个模块的 edges 文件” -> 得到自选基因集
############################################################
# 你说：M2 是一个文件；M1/M3 也各有自己的文件
# 其他细胞类型模块数量不同：你就按实际模块增删这一段即可

# >>>>>> 你只需要改这里：模块名 -> 对应 edges 文件 <<<<<<
module_edges <- list(
  "Mesophyll-M2" = "Mesophyll_M2_edges_r_gt0.7_noNucNuc.csv",
  "Mesophyll-M1" = "Mesophyll_M1_edges_noNucNuc_absR_gt_0.5.csv"
)
# >>>>>> 改完就行 <<<<<<

# 安全检查：文件存在
for (m in names(module_edges)) {
  f <- module_edges[[m]]
  if (!file.exists(f)) stop("Edges file not found for ", m, ": ", f)
}

############################################################
## 4) 出图函数：dotplot（自动换行避免Y轴重叠）
############################################################

safe_dotplot_pdf <- function(x, pdf_file,
                             showCategory = 20,
                             title = "",
                             label_format = 35,
                             font.size = 9,
                             width = 10,
                             height = 8) {
  pdf(pdf_file, width = width, height = height)
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo enriched terms passed cutoff."), cex = 1)
  } else {
    p <- dotplot(
      x,
      showCategory = showCategory,
      label_format = label_format,
      font.size = font.size
    ) + ggtitle(title)
    print(p)
  }
  dev.off()
}

############################################################
## 5) 方案二输出目录（不会覆盖方案一）
############################################################

plot_dir <- "Scheme2_CustomHub_ByModule_EnrichPlots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

############################################################
## 6) 逐模块：从 edges 得到基因集 -> GO/KEGG -> 画图
############################################################

showN <- 20
wrapN <- 35
fs    <- 9
W     <- 10
H     <- 8

# 存结果（可选，方便你后面保存rds）
res_go_bp <- list()
res_go_cc <- list()
res_go_mf <- list()
res_kegg  <- list()
res_genes <- list()

for (m in names(module_edges)) {

  cat("\n=== Processing module:", m, "===\n")

  edges <- fread(module_edges[[m]])

  # 兼容你当前 edges 表的列名 Gene1/Gene2
  stopifnot(all(c("Gene1","Gene2") %in% names(edges)))

  gene_list_raw <- unique(c(edges$Gene1, edges$Gene2))

  # 关键：与背景集取交集（保证富集统计一致、可解释）
  gene_list_in_universe <- intersect(gene_list_raw, universe)

  cat("Raw genes:", length(gene_list_raw),
      " | In universe:", length(gene_list_in_universe), "\n")

  res_genes[[m]] <- gene_list_in_universe

  # 富集
  ego_bp <- run_go(gene_list_in_universe, universe, "BP")
  ego_cc <- run_go(gene_list_in_universe, universe, "CC")
  ego_mf <- run_go(gene_list_in_universe, universe, "MF")
  ekegg  <- run_kegg_tair(gene_list_in_universe, universe, organism="ath")

  res_go_bp[[m]] <- ego_bp
  res_go_cc[[m]] <- ego_cc
  res_go_mf[[m]] <- ego_mf
  res_kegg[[m]]  <- ekegg

  # 画图（文件名不改也没事，因为目录已换成Scheme2，不会覆盖方案一）
  safe_dotplot_pdf(ego_bp,
                   file.path(plot_dir, paste0(m, "_CustomHub_GO_BP_dotplot.pdf")),
                   showN, paste0(m, " CustomHub GO BP"), wrapN, fs, W, H)

  safe_dotplot_pdf(ego_cc,
                   file.path(plot_dir, paste0(m, "_CustomHub_GO_CC_dotplot.pdf")),
                   showN, paste0(m, " CustomHub GO CC"), wrapN, fs, W, H)

  safe_dotplot_pdf(ego_mf,
                   file.path(plot_dir, paste0(m, "_CustomHub_GO_MF_dotplot.pdf")),
                   showN, paste0(m, " CustomHub GO MF"), wrapN, fs, W, H)

  safe_dotplot_pdf(ekegg,
                   file.path(plot_dir, paste0(m, "_CustomHub_KEGG_dotplot.pdf")),
                   showN, paste0(m, " CustomHub KEGG"), wrapN, fs, W, H)
}

cat("\nDONE.\nPlots in: ", normalizePath(plot_dir), "\n", sep="")

############################################################
## 7) （可选）把方案二结果保存成 rds（方便你后续复用/对比）
############################################################
scheme2_out <- list(
  rds_file   = rds_file,
  universe   = universe,
  module_edges = module_edges,
  genes_by_module = res_genes,
  go_bp = res_go_bp,
  go_cc = res_go_cc,
  go_mf = res_go_mf,
  kegg  = res_kegg
)

saveRDS(scheme2_out, file = "Scheme2_CustomHub_Enrichment_Result.rds")
cat("Saved: Scheme2_CustomHub_Enrichment_Result.rds\n")
