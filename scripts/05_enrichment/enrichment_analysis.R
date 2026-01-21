############################################################
## 0) 环境准备：包检查 + 安装（缺啥装啥）
############################################################

# CRAN 包
cran_pkgs <- c("data.table", "dplyr", "ggplot2")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos="https://cloud.r-project.org")
}

# Bioconductor 包：clusterProfiler + org.At.tair.db
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
bioc_pkgs <- c("clusterProfiler", "org.At.tair.db")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

library(data.table)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.At.tair.db)

############################################################
## 1) 读入 hdWGCNA rds：抽出 active_wgcna
############################################################

# >>>>>> 你只需要改这里 <<<<<<
rds_file     <- "Atha_leaf_PRJNA648028_Mesophyll_hdwgcna.rds"  # 改成你的rds
celltype_tag <- "Mesophyll"                                   # 只是输出文件名标签
# >>>>>> 改完就行 <<<<<<

obj <- readRDS(rds_file)

# hdWGCNA 通常放在 obj@misc$active_wgcna
stopifnot("active_wgcna" %in% names(obj@misc))
w <- obj@misc$active_wgcna

# 模块/连边/中心性主要用到这俩表
mods <- w$wgcna_modules   # gene -> module + kME columns
deg  <- w$wgcna_degrees   # gene -> module + degree/weighted_degree

# 基础清洗
mods$gene_name <- as.character(mods$gene_name)
mods$module    <- as.character(mods$module)
deg$gene_name  <- as.character(deg$gene_name)
deg$module     <- as.character(deg$module)

# 去掉 grey（未分配模块的基因）
mods_ng  <- subset(mods, module != "grey")

# 背景基因（universe）：本细胞类型进入 WGCNA 且不在 grey 的基因
universe <- unique(mods_ng$gene_name)

# 每个模块的基因列表
module_genes <- split(mods_ng$gene_name, mods_ng$module)
modules <- names(module_genes)

cat("Modules:\n"); print(table(mods_ng$module))
cat("Universe size:", length(universe), "\n")

############################################################
## 2) 富集函数：GO（BP/CC/MF） + KEGG（用 TAIR ID 直接跑）
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

# 说明：你这套环境下 enrichKEGG(organism="ath") 能直接吃 TAIR (ATxGxxxxx)
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
## 3) 模块富集：GO BP/CC/MF + KEGG
############################################################

mod_go_bp <- list(); mod_go_cc <- list(); mod_go_mf <- list(); mod_kegg <- list()

for (m in modules) {
  mod_go_bp[[m]] <- run_go(module_genes[[m]], universe, "BP")
  mod_go_cc[[m]] <- run_go(module_genes[[m]], universe, "CC")
  mod_go_mf[[m]] <- run_go(module_genes[[m]], universe, "MF")
  mod_kegg[[m]]  <- run_kegg_tair(module_genes[[m]], universe, organism="ath")
}

# （可选）快速查看 top5（避免 geneID 超长输出导致断线）
peek_enrich <- function(x, n=5) {
  if (is.null(x)) return(NULL)
  head(as.data.frame(x)[, c("ID","Description","p.adjust","Count")], n)
}
cat("\n== Module GO BP top5 ==\n")
print(lapply(mod_go_bp, peek_enrich, n=5))

############################################################
## 4) Hub 基因：kME / weighted_degree / 交集 / hub_score
############################################################

# 取“模块内kME列名”，例如 kME_Mesophyll-M1
get_kme_col <- function(module_name) paste0("kME_", module_name)

# 把 mods 与 deg 合并（按 gene+module）
md <- mods_ng %>%
  select(gene_name, module, starts_with("kME_")) %>%
  left_join(deg %>% select(gene_name, module, weighted_degree),
            by=c("gene_name","module"))

# 计算每个基因在所属模块的 kME_inModule
md$kME_inModule <- NA_real_
for (m in modules) {
  kcol <- get_kme_col(m)
  if (!kcol %in% colnames(md)) stop("Missing kME column: ", kcol)
  idx <- md$module == m
  md$kME_inModule[idx] <- md[idx, kcol]
}

# hub_score：z(kME_inModule) + z(weighted_degree)（在模块内做z）
md <- md %>%
  group_by(module) %>%
  mutate(
    z_kME = as.numeric(scale(kME_inModule)),
    z_wD  = as.numeric(scale(weighted_degree)),
    hub_score = z_kME + z_wD
  ) %>%
  ungroup()

# 输出每个模块的 topN hub（你可以改 topN）
topN <- 100

hub_kME <- list()
hub_deg <- list()
hub_intersect <- list()
hub_score <- list()

for (m in modules) {
  sub <- md %>% filter(module == m)

  hub_kME[[m]] <- sub %>% arrange(desc(kME_inModule)) %>% head(topN) %>%
    select(gene_name, kME_inModule, weighted_degree, hub_score)

  hub_deg[[m]] <- sub %>% arrange(desc(weighted_degree)) %>% head(topN) %>%
    select(gene_name, kME_inModule, weighted_degree, hub_score)

  hub_score[[m]] <- sub %>% arrange(desc(hub_score)) %>% head(topN) %>%
    select(gene_name, kME_inModule, weighted_degree, hub_score)

  hub_intersect[[m]] <- intersect(hub_kME[[m]]$gene_name, hub_deg[[m]]$gene_name)
}

# （可选）看看每个模块 top10 by hub_score
for (m in modules) {
  cat("\n==== ", m, " Top10 by hub_score ====\n", sep="")
  print(hub_score[[m]] %>% head(10))
}

# 导出 hub 表
out_tab_dir <- "HubTables"
if (!dir.exists(out_tab_dir)) dir.create(out_tab_dir)

for (m in modules) {
  fwrite(hub_kME[[m]],  file = file.path(out_tab_dir, paste0(m, "_HubTop", topN, "_by_kME.tsv")), sep="\t")
  fwrite(hub_deg[[m]],  file = file.path(out_tab_dir, paste0(m, "_HubTop", topN, "_by_wDegree.tsv")), sep="\t")
  fwrite(data.frame(gene_name = hub_intersect[[m]]),
         file = file.path(out_tab_dir, paste0(m, "_HubIntersectTop", topN, ".tsv")), sep="\t")
  fwrite(hub_score[[m]], file = file.path(out_tab_dir, paste0(m, "_HubTop", topN, "_by_Score.tsv")), sep="\t")
}

############################################################
## 5) Hub 富集：先取 topN 的 hub_score 基因，再做 GO/KEGG
############################################################

# hub 基因集合：每模块取 hub_score topN
hub_genes <- lapply(modules, function(m) hub_score[[m]]$gene_name[1:topN])
names(hub_genes) <- modules

# 去掉奇怪 ID（lncRNA / STRG）
is_weird <- function(x) grepl("^ATHA-LNC|^STRG\\.", x)
hub_genes_clean <- lapply(hub_genes, function(v) v[!is_weird(v)])
names(hub_genes_clean) <- modules

hub_go_bp <- list(); hub_go_cc <- list(); hub_go_mf <- list(); hub_kegg <- list()

for (m in modules) {
  hub_go_bp[[m]] <- run_go(hub_genes_clean[[m]], universe, "BP")
  hub_go_cc[[m]] <- run_go(hub_genes_clean[[m]], universe, "CC")
  hub_go_mf[[m]] <- run_go(hub_genes_clean[[m]], universe, "MF")
  hub_kegg[[m]]  <- run_kegg_tair(hub_genes_clean[[m]], universe, organism="ath")
}

# （可选）快速查看 hub GO BP top10（避免 geneID）
peek_enrich2 <- function(x, n=10) {
  if (is.null(x)) return(NULL)
  head(as.data.frame(x)[, c("Description","p.adjust","Count")], n)
}
cat("\n== Hub GO BP top10 ==\n")
print(lapply(hub_go_bp, peek_enrich2, n=10))

############################################################
## 6) 出图：统一风格（label_format 自动换行，避免纵轴重叠）
############################################################

plot_dir <- "EnrichPlots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

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

# 参数可调：showN 越小越不挤；wrapN 越小换行越频繁；H 越大越不挤
showN <- 20
wrapN <- 35
fs    <- 9
W     <- 10
H     <- 8

for (m in modules) {

  ## --- Module GO ---
  safe_dotplot_pdf(mod_go_bp[[m]],
                   file.path(plot_dir, paste0(m, "_Module_GO_BP_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Module GO BP"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  safe_dotplot_pdf(mod_go_cc[[m]],
                   file.path(plot_dir, paste0(m, "_Module_GO_CC_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Module GO CC"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  safe_dotplot_pdf(mod_go_mf[[m]],
                   file.path(plot_dir, paste0(m, "_Module_GO_MF_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Module GO MF"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  ## --- Module KEGG ---
  safe_dotplot_pdf(mod_kegg[[m]],
                   file.path(plot_dir, paste0(m, "_Module_KEGG_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Module KEGG"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  ## --- Hub GO ---
  safe_dotplot_pdf(hub_go_bp[[m]],
                   file.path(plot_dir, paste0(m, "_Hub_GO_BP_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Hub GO BP"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  safe_dotplot_pdf(hub_go_cc[[m]],
                   file.path(plot_dir, paste0(m, "_Hub_GO_CC_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Hub GO CC"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  safe_dotplot_pdf(hub_go_mf[[m]],
                   file.path(plot_dir, paste0(m, "_Hub_GO_MF_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Hub GO MF"),
                   label_format=wrapN, font.size=fs, width=W, height=H)

  ## --- Hub KEGG ---
  safe_dotplot_pdf(hub_kegg[[m]],
                   file.path(plot_dir, paste0(m, "_Hub_KEGG_dotplot.pdf")),
                   showCategory=showN, title=paste0(m, " Hub KEGG"),
                   label_format=wrapN, font.size=fs, width=W, height=H)
}

cat("\nDONE.\n",
    "Hub tables in: ", normalizePath(out_tab_dir), "\n",
    "Plots in:      ", normalizePath(plot_dir), "\n", sep="")
