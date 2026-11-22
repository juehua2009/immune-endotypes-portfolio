# Module 1 (R twin) — 02: Harmony integration (the R counterpart of notebook 02's scVI)
# batch = Sample, as in the scVI setup. Then UMAP + clustering + purity vs author labels.
# Starts from the counts-only so_qc.rds saved by 01; PCA runs on the HVG subset.
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(patchwork)
})
datadir <- path.expand("~/Documents/immune-endotypes/data/smillie_uc")
figdir <- path.expand("~/Documents/github/immune-endotypes-portfolio/modules/01_ibd/figures")
so <- readRDS(file.path(datadir, "so_qc.rds"))

so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, selection.method = "mvp", nfeatures = 3000, verbose = FALSE)
so <- subset(so, features = VariableFeatures(so))
gc()
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 50, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "Sample", verbose = FALSE)
so <- FindNeighbors(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
cat("harmony clusters:", length(unique(so$seurat_clusters)), "(scVI gave 28)\n")

p1 <- DimPlot(so, group.by = "condition", raster = TRUE) + ggtitle("condition")
p2 <- DimPlot(so, group.by = "Major_l2", raster = TRUE) + ggtitle("Major_l2")
p3 <- DimPlot(so, group.by = "author_label", raster = TRUE, label = TRUE, repel = TRUE) +
  ggtitle("author_label") + NoLegend()
p4 <- DimPlot(so, group.by = "seurat_clusters", raster = TRUE, label = TRUE) +
  ggtitle("Harmony clusters (res 0.5)") + NoLegend()
ggsave(file.path(figdir, "umap_r_harmony_overview.png"),
       (p1 | p2) / (p3 | p4), width = 16, height = 14, dpi = 120)

# Purity: dominant author label per cluster (mirrors notebook 02's crosstab)
ct <- table(so$seurat_clusters, so$author_label)
purity <- data.frame(
  cluster = rownames(ct),
  top_label = colnames(ct)[max.col(ct)],
  share = round(apply(ct, 1, max) / rowSums(ct), 2)
)
print(purity[order(-purity$share), ], row.names = FALSE)

saveRDS(so, file.path(datadir, "so_harmony.rds"))
cat("saved so_harmony.rds\n")
