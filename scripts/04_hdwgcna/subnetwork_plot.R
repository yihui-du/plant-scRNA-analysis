library(SeuratObject)
library(data.table)
library(dplyr)
library(tidyr)

## 0) 参数
rds_file <- "Atha_root_PRJNA497883_Atrichoblast_hdwgcna.rds"
cutoff   <- 0.3   # |r| 阈值

## 1) 读取对象 & 取 WGCNA 实验

obj <- readRDS(rds_file)
# 如果有 active_wgcna，就用它
wgc_name <- obj@misc$active_wgcna
cat("Active WGCNA experiment:", wgc_name, "\n")
exp     <- obj@misc[[wgc_name]]
mods    <- exp$wgcna_modules
datExpr <- exp$datExpr

## 2) 先看所有模块的基因信息

type_of <- function(g){
  ifelse(grepl("^ATMG", g), "Mitochondria",
         ifelse(grepl("^ATCG", g), "Chloroplast",
                ifelse(grepl("^AT[1-5]G", g), "Nuclear", NA_character_)))
}
mods2 <- mods %>%
  mutate(
    module   = as.character(module),
    gene_type = type_of(gene_name)
  )
modules_all <- sort(unique(mods2$module))
cat("共有模块数:", length(modules_all), "\n")
cat("模块列表:\n")
print(modules_all)

gene_summary <- mods2 %>%
  filter(!is.na(gene_type)) %>%
  count(module, gene_type) %>%
  pivot_wider(
    names_from  = gene_type,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total_genes = Chloroplast + Mitochondria + Nuclear
  ) %>%
  arrange(module) %>%
  select(module, Chloroplast, Mitochondria, Nuclear, total_genes)
cat("\n每个模块的基因构成（叶绿体/线粒体/核）：\n")
print(gene_summary)

## 3) 在看到信息之后，选择一个模块做子网

cat("\n请输入要做子网分析的模块名（例如 Epidermis-M2），直接回车确认：\n")
module_name <- readline(prompt = "module = ")
module_name <- trimws(module_name)

out_file <- paste0(module_name, "_subnetwork_absR_gt_0.3_merged.csv")

## 4) 取这个模块的基因

m_genes <- mods2 %>%
  filter(module == module_name) %>%
  pull(gene_name) %>%
  unique()
cat("该模块基因数（wgcna_modules 中）:", length(m_genes), "\n")
m_genes <- intersect(m_genes, colnames(datExpr))
cat("该模块基因数（在 datExpr 中）:", length(m_genes), "\n")
mat <- as.matrix(datExpr[, m_genes, drop = FALSE])

## 5) 模块内相关矩阵 & 转边表（只取上三角 + |r|>cutoff）

cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
idx <- which(abs(cor_mat) > cutoff & upper.tri(cor_mat), arr.ind = TRUE)
edges <- data.table(
  Gene1 = colnames(cor_mat)[idx[, 1]],
  Gene2 = colnames(cor_mat)[idx[, 2]],
  r     = cor_mat[idx]
)
cat("Edges before type filter:", nrow(edges), "\n")

## 6) 严格基因类型 & EdgeClass

edges[, Type1 := type_of(Gene1)]
edges[, Type2 := type_of(Gene2)]
edges <- edges[!is.na(Type1) & !is.na(Type2)]
cat("Edges after type filter:", nrow(edges), "\n")
edges[, EdgeClass := fifelse(
  Type1 == "Nuclear" & Type2 == "Nuclear", "Nuc-Nuc",
  fifelse(Type1 != "Nuclear" & Type2 != "Nuclear", "Org-Org", "Org-Nuc")
)]
cat("\nEdgeClass counts:\n")
print(edges[, .N, by = EdgeClass][order(-N)])

## 7) 按 |r| 从大到小排序并输出

edges[, abs_r := abs(r)]
setorder(edges, -abs_r)
edges[, abs_r := NULL]
fwrite(edges, out_file)
cat("\n子网已保存到:", out_file, "\n")

## 8) 在此基础上继续筛一份：去掉 Nuc-Nuc，且 |r| > 0.7

edges_strict <- edges[EdgeClass != "Nuc-Nuc" & abs(r) > 0.7]
out_file2 <- paste0(module_name, "_subnetwork_noNucNuc_absR_gt_0.7.csv")
fwrite(edges_strict, out_file2)
cat("去掉 Nuc-Nuc 且 |r|>0.7 的子网已保存到:", out_file2, "\n")

edges_strict <- edges[EdgeClass != "Nuc-Nuc" & abs(r) > 0.5]
out_file3 <- paste0(module_name, "_subnetwork_noNucNuc_absR_gt_0.5.csv")
fwrite(edges_strict, out_file3)
cat("去掉 Nuc-Nuc 且 |r|>0.5 的子网已保存到:", out_file3, "\n")
