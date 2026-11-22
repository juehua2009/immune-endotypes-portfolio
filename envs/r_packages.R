# R environment for Module 1 (Seurat version) and pseudobulk DE
install.packages(c("Seurat", "harmony", "tidyverse", "patchwork"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "muscat", "SingleR", "celldex"))
# CellChat: remotes::install_github("sqjin/CellChat")
