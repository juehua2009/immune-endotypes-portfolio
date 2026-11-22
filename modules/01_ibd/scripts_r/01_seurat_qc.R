# Module 1 (R twin) — 01: Seurat load, QC, PCA, UMAP, clustering
# Mirrors Python notebook 01 (scanpy). Input: r_export/ from 00_export_from_python.py
# Memory notes (352k cells on a 32 Gb machine): the counts-only object is saved FIRST
# (so_qc.rds, the input for scripts 02 and 03); PCA/UMAP then run on an HVG-subset
# object so ScaleData's dense matrix stays small.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork)
})
base <- path.expand("~/Documents/immune-endotypes/data/smillie_uc/r_export")
datadir <- path.expand("~/Documents/immune-endotypes/data/smillie_uc")
figdir <- path.expand("~/Documents/github/immune-endotypes-portfolio/modules/01_ibd/figures")

counts <- readMM(file.path(base, "counts_genes_by_cells.mtx.gz"))
obs <- read.csv(file.path(base, "obs.csv"), row.names = 1, check.names = FALSE)
genes <- read.csv(file.path(base, "genes.csv"))$gene
barcodes <- read.csv(file.path(base, "barcodes.csv"))$barcode
rownames(counts) <- genes
colnames(counts) <- barcodes
stopifnot(all(colnames(counts) == rownames(obs)))

so <- CreateSeuratObject(counts = counts, meta.data = obs)
rm(counts, obs, genes, barcodes); gc()
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
# Same gates as the Python notebook: >200 genes per cell, <50% mitochondrial
so <- subset(so, subset = nFeature_RNA > 200 & percent.mt < 50)
cat("cells after QC:", ncol(so), "(Python notebook 01 kept 351,952)\n")
saveRDS(so, file.path(datadir, "so_qc.rds"))
cat("saved counts-only so_qc.rds (input for 02 and 03)\n")

# Continue on an HVG-subset object (PCA/UMAP never use all genes)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, selection.method = "mvp", nfeatures = 3000, verbose = FALSE)
so <- subset(so, features = VariableFeatures(so))
gc()
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 50, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:40, verbose = FALSE)
so <- RunUMAP(so, dims = 1:40, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
cat("clusters:", length(unique(so$seurat_clusters)), "\n")

p1 <- DimPlot(so, group.by = "condition", raster = TRUE) + ggtitle("condition")
p2 <- DimPlot(so, group.by = "Major_l2", raster = TRUE) + ggtitle("Major_l2")
p3 <- DimPlot(so, group.by = "author_label", raster = TRUE, label = TRUE, repel = TRUE) +
  ggtitle("author_label") + NoLegend()
p4 <- DimPlot(so, group.by = "seurat_clusters", raster = TRUE, label = TRUE) +
  ggtitle("Seurat clusters (res 0.5)") + NoLegend()
ggsave(file.path(figdir, "umap_r_seurat_overview.png"),
       (p1 | p2) / (p3 | p4), width = 16, height = 14, dpi = 120)
saveRDS(so, file.path(datadir, "so_umap.rds"))
cat("saved so_umap.rds\n")
