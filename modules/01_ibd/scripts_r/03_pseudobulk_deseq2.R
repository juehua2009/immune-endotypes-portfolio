# Module 1 (R twin) — 03: Pseudobulk DESeq2 + anti-TNF module
# Mirrors Python notebook 03. Pseudobulk unit = Subject x condition x author_label
# (the file's "Sample" column is a processing-batch code, not a sample ID).
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(DESeq2); library(ggplot2); library(patchwork)
})
repo <- path.expand("~/Documents/github/immune-endotypes-portfolio")
resdir <- file.path(repo, "modules/01_ibd/results")
figdir <- file.path(repo, "modules/01_ibd/figures")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
so <- readRDS(path.expand("~/Documents/immune-endotypes/data/smillie_uc/so_qc.rds"))
# so_qc.rds is counts-only; AddModuleScore below needs the normalized layer
so <- NormalizeData(so, verbose = FALSE)

## ---- Pseudobulk aggregation (sparse indicator multiply) ----
md <- so[[]]
key <- paste(md$Subject, md$condition, md$author_label, sep = "|")
groups <- unique(key)
gi <- match(key, groups)
ind <- sparseMatrix(i = gi, j = seq_len(ncol(so)), dims = c(length(groups), ncol(so)))
pb <- ind %*% t(GetAssayData(so, layer = "counts"))  # groups x genes
pb <- as(pb, "dgCMatrix")
rownames(pb) <- groups  # sparseMatrix dropped the dimnames; the DE loop subsets by name
meta <- do.call(rbind, strsplit(groups, "|", fixed = TRUE))
meta <- data.frame(Subject = meta[, 1], condition = meta[, 2], label = meta[, 3],
                   n_cells = as.integer(table(factor(gi, levels = seq_along(groups)))),
                   row.names = groups, stringsAsFactors = FALSE)
keep <- meta$n_cells >= 20
pb <- pb[keep, ]; meta <- meta[keep, ]
cat("pseudobulk groups kept:", nrow(meta), "\n")

## ---- DESeq2 per label: Inflamed vs Healthy ----
results <- list()
for (lab in sort(unique(meta$label))) {
  m <- meta[meta$label == lab, ]
  if (sum(m$condition == "Healthy") < 3 || sum(m$condition == "Inflamed") < 3) next
  cm <- t(as.matrix(pb[rownames(m), ]))  # genes x samples
  cm <- cm[rowSums(cm >= 10) >= 3, , drop = FALSE]
  if (nrow(cm) < 100) next
  m$condition <- factor(m$condition, levels = c("Healthy", "Non-inflamed", "Inflamed"))
  dds <- DESeqDataSetFromMatrix(countData = round(cm), colData = m, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "Inflamed", "Healthy")))
  res <- res[order(res$padj), ]
  results[[lab]] <- res
  write.csv(res, file.path(resdir, paste0("deseq2_r_inflamed_vs_healthy_",
              gsub(" ", "_", gsub("\\+", "pos", lab)), ".csv")))
  cat(sprintf("%-35s DE genes (padj<0.05): %d\n", lab, sum(res$padj < 0.05, na.rm = TRUE)))
}
cat("labels tested:", length(results), "\n")

## ---- Volcano: Plasma (the label that lit up in Python) ----
if ("Plasma" %in% names(results)) {
  r <- results[["Plasma"]]
  r$gene <- rownames(r)
  r$sig <- !is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) > 1
  top <- head(r[r$sig, ][order(r[r$sig, ]$padj), ], 6)
  p <- ggplot(r, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_text(data = top, aes(label = gene), size = 3, vjust = -0.6, colour = "black") +
    scale_colour_manual(values = c("grey", "firebrick"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Plasma — Inflamed vs Healthy (R/DESeq2)",
         x = "log2 fold change", y = "-log10 padj") + theme_bw()
  ggsave(file.path(figdir, "volcano_r_plasma.png"), p, width = 6, height = 5, dpi = 150)
}

## ---- Anti-TNF module: abundance + program scores per Subject x condition ----
fib_label <- grep("inflammatory fibroblast", unique(md$author_label),
                  ignore.case = TRUE, value = TRUE)
mono_label <- grep("inflammatory monocyte", unique(md$author_label),
                   ignore.case = TRUE, value = TRUE)
abundance <- function(target, compartment) {
  comp <- md[md$Major_l2 == compartment, ]
  tot <- as.data.frame(table(comp$Subject, comp$condition))
  tgt <- as.data.frame(table(comp$Subject[comp$author_label == target],
                             comp$condition[comp$author_label == target]))
  out <- merge(tot, tgt, by = c("Var1", "Var2"), all.x = TRUE, suffixes = c("_tot", "_tgt"))
  out$Freq_tgt[is.na(out$Freq_tgt)] <- 0
  data.frame(Subject = out$Var1, condition = out$Var2,
             fraction = out$Freq_tgt / out$Freq_tot)
}
ab_fib <- abundance(fib_label, "Stromal")
ab_mono <- abundance(mono_label, "Myeloid")

fib_genes <- list(c("IL11", "IL13RA2", "IL24", "CHI3L1", "MMP3", "MMP10",
                    "IL6", "CXCL1", "CXCL5"))
mono_genes <- list(c("OSM", "OSMR", "IL1B", "IL6", "TNF", "S100A8", "S100A9"))
so <- AddModuleScore(so, features = fib_genes, name = "fib_program")
so <- AddModuleScore(so, features = mono_genes, name = "mono_program")
md2 <- so[[]]
sample_scores <- function(label, col) {
  sub <- md2[md2$author_label == label, ]
  agg <- aggregate(sub[[col]], by = list(sub$Subject, sub$condition), mean)
  data.frame(Subject = agg$Group.1, condition = agg$Group.2, score = agg$x)
}
sc_fib <- sample_scores(fib_label, "fib_program1")
sc_mono <- sample_scores(mono_label, "mono_program1")

strip <- function(df, y, title) {
  df$condition <- factor(df$condition, levels = c("Healthy", "Non-inflamed", "Inflamed"))
  ggplot(df, aes(condition, .data[[y]])) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.7) +
    ggtitle(title) + theme_bw() + xlab(NULL)
}
p <- strip(ab_fib, "fraction", "Inflammatory fibroblasts\n(share of stromal)") |
     strip(sc_fib, "score", "Fibroblast program score\n(IL11 / IL13RA2 axis)") |
     strip(ab_mono, "fraction", "Inflammatory monocytes\n(share of myeloid)") |
     strip(sc_mono, "score", "Monocyte program score\n(OSM / IL1B axis)")
ggsave(file.path(figdir, "antitnf_module_summary_r.png"), p, width = 16, height = 4.2, dpi = 150)
cat("module figure saved\n")
