suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# =========================
# Parameters
# =========================
rds_path <- "/data/duyihui/scRNA/Atha_leaf_PRJNA648028_allgene_rename.rds"
outdir   <- "/data/duyihui/scRNA/Marker_results/Atha_leaf_PRJNA648028_allgene_rename"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

celltype_col <- "cell_type"

# Marker defaults
min_pct <- 0.10
logfc_threshold <- 0.25
padj_cutoff <- 0.05
test_use <- "wilcox"
only_pos <- TRUE
top_n <- 50

# How to handle NA cell types
# Option A: drop NA (recommended)
drop_na_celltype <- TRUE
# Option B: keep NA as "unknown"
# drop_na_celltype <- FALSE

# =========================
# Load object
# =========================
obj <- readRDS(rds_path)

stopifnot("RNA" %in% Assays(obj))
DefaultAssay(obj) <- "RNA"

# =========================
# Handle celltype NA
# =========================
ct <- obj@meta.data[[celltype_col]]

if (drop_na_celltype) {
  keep <- !is.na(ct) & ct != "" & ct != "NA"
  obj <- subset(obj, cells = colnames(obj)[keep])
} else {
  obj@meta.data[[celltype_col]] <- as.character(ct)
  obj@meta.data[[celltype_col]][is.na(obj@meta.data[[celltype_col]])] <- "unknown"
  obj@meta.data[[celltype_col]] <- factor(obj@meta.data[[celltype_col]])
}

# Save celltype counts (sanity check)
ct_counts <- table(obj@meta.data[[celltype_col]])
write.csv(as.data.frame(ct_counts),
          file.path(outdir, "celltype_counts_after_filter.csv"),
          row.names = FALSE)

# Set identities to cell type
Idents(obj) <- obj@meta.data[[celltype_col]]

# =========================
# FindAllMarkers across all cell types
# =========================
markers_all <- FindAllMarkers(
  object = obj,
  only.pos = only_pos,
  test.use = test_use,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# Harmonize column names across Seurat versions
if ("avg_logFC" %in% colnames(markers_all) && !"avg_log2FC" %in% colnames(markers_all)) {
  markers_all$avg_log2FC <- markers_all$avg_logFC
}

if ("p_val_adj" %in% colnames(markers_all)) {
  markers_all$signif <- markers_all$p_val_adj < padj_cutoff
}

# Write full marker table
write.csv(markers_all,
          file.path(outdir, "markers_FindAllMarkers.csv"),
          row.names = FALSE)

# =========================
# Top N markers per cell type
# =========================
if ("p_val_adj" %in% colnames(markers_all)) {
  top_tbl <- markers_all %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = top_n) %>%
    ungroup()
} else {
  top_tbl <- markers_all %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

write.csv(top_tbl,
          file.path(outdir, paste0("markers_Top", top_n, "_per_celltype.csv")),
          row.names = FALSE)

# =========================
# Summary counts per cell type
# =========================
if ("p_val_adj" %in% colnames(markers_all)) {
  summary_tbl <- markers_all %>%
    group_by(cluster) %>%
    summarise(
      n_total = n(),
      n_sig = sum(p_val_adj < padj_cutoff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig))

  write.csv(summary_tbl,
            file.path(outdir, "markers_summary_by_celltype.csv"),
            row.names = FALSE)
}

# =========================
# Optional: one-vs-rest example
# =========================
target <- "Guard cells"
if (target %in% levels(Idents(obj))) {
  mk_one <- FindMarkers(
    object = obj,
    ident.1 = target,
    ident.2 = NULL,
    test.use = test_use,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = only_pos
  )
  mk_one$gene <- rownames(mk_one)
  write.csv(mk_one,
            file.path(outdir, paste0("markers_", gsub(" ", "_", target), "_vs_rest.csv")),
            row.names = FALSE)
}

message("Done. Results in: ", outdir)
