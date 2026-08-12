# ==============================================
# COMPLETE TB TRANSCRIPTOMICS ANALYSIS
# GSE37250: Active TB vs LTBI
# ==============================================

# 1. SETUP & LIBRARIES
setwd("D:/OneDrive/Desktop/TB-Host-Transcriptomics-main/TB-Host-Transcriptomics-main")

# install.packages("tidyverse")
# install.packages("pheatmap")
# BiocManager::install("limma")
# BiocManager::install("GEOquery")

library(limma)
library(GEOquery)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(dplyr)

# 2. LOAD RAW DATA (Pre-normalized, but linear)
cat("STEP 1: Loading raw expression matrix...\n")
expr <- read.csv("expression_matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("metadata.csv", check.names = FALSE)
all(colnames(expr) == meta$gsm_id)

# 3. LOG2 TRANSFORMATION WITH OFFSET (Handles negatives)
cat("STEP 2: Applying log2 transformation with offset...\n")
offset <- -min(expr) + 1
cat("Offset applied:", offset, "\n")
expr_log <- log2(expr + offset)
cat("Log2 range:", round(min(expr_log), 2), "to", round(max(expr_log), 2), "\n")

# 4. PROBE FILTERING (IQR > 0.1)
cat("STEP 3: Filtering low-variance probes...\n")
probe_iqr <- apply(expr_log, 1, IQR, na.rm = TRUE)
keep <- probe_iqr > 0.1
expr_final <- expr_log[keep, ]
cat("Probes kept:", nrow(expr_final), "out of", nrow(expr_log), "\n")

# 5. ANNOTATION (Fetch from GEO)
cat("STEP 4: Downloading annotation from GPL10558...\n")
gpl <- getGEO("GPL10558", destdir = ".")
annot_table <- Table(gpl)
cat("Annotation rows:", nrow(annot_table), "\n")
if (!"Symbol" %in% colnames(annot_table)) stop("Could not find 'Symbol' column!")
annot_lookup <- annot_table[, c("ID", "Symbol")]
colnames(annot_lookup) <- c("ProbeID", "GeneSymbol")
annot_lookup <- annot_lookup[!is.na(annot_lookup$GeneSymbol) & annot_lookup$GeneSymbol != "", ]

probe_ids <- rownames(expr_final)
gene_annotation <- data.frame(ProbeID = probe_ids, stringsAsFactors = FALSE)
gene_annotation$GeneSymbol <- annot_lookup$GeneSymbol[match(probe_ids, annot_lookup$ProbeID)]

keep_mapped <- !is.na(gene_annotation$GeneSymbol)
expr_mapped <- expr_final[keep_mapped, ]
gene_annotation_mapped <- gene_annotation[keep_mapped, ]
cat("Mapped probes:", nrow(expr_mapped), "\n")
cat("Unique genes:", length(unique(gene_annotation_mapped$GeneSymbol)), "\n")

# 6. LIMMA DIFFERENTIAL EXPRESSION (WITH REGION COVARIATE)
cat("STEP 5: Running Limma with Region as covariate...\n")

design <- model.matrix(~ 0 + disease + region, data = meta)
colnames(design) <- c("ActiveTB", "LatentTB", "RegionSA")
cat("Design columns:\n")
print(colnames(design))

fit <- lmFit(expr_mapped, design)

contrast_matrix <- makeContrasts(Active_vs_LTBI = ActiveTB - LatentTB, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = "Active_vs_LTBI", number = Inf, adjust.method = "BH")
results$GeneSymbol <- gene_annotation_mapped$GeneSymbol[match(rownames(results), rownames(expr_mapped))]

sig_deg <- results[results$adj.P.Val < 0.05, ]
up <- sig_deg[sig_deg$logFC > 0.5, ]
down <- sig_deg[sig_deg$logFC < -0.5, ]
cat("\n--- RESULTS ---\n")
cat("Total probes tested:", nrow(results), "\n")
cat("DEGs (FDR < 0.05):", nrow(sig_deg), "\n")
cat("Upregulated (logFC > 0.5):", nrow(up), "\n")
cat("Downregulated (logFC < -0.5):", nrow(down), "\n")

# 7. SAVE OUTPUTS
write.csv(results, "DE_results_limma_region.csv", row.names = TRUE)
write.csv(up$GeneSymbol, "upregulated_genes.csv", row.names = FALSE)
write.csv(down$GeneSymbol, "downregulated_genes.csv", row.names = FALSE)
cat("Saved results, up/down gene lists.\n")

# 8. VISUALIZATIONS
cat("STEP 6: Generating plots...\n")

# Volcano
results$Sig <- "Not Significant"
results$Sig[results$adj.P.Val < 0.05 & results$logFC > 0.5] <- "Upregulated"
results$Sig[results$adj.P.Val < 0.05 & results$logFC < -0.5] <- "Downregulated"

volcano <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "steelblue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: Active TB vs LTBI (Region-adjusted)", x = "Log2 FC", y = "-Log10 FDR") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave("volcano_plot.png", volcano, width = 10, height = 7, dpi = 300)

# Heatmap (Top 50)
top_genes <- rownames(results[order(results$adj.P.Val), ])[1:50]
top_expr <- expr_mapped[top_genes, ]
top_scaled <- t(scale(t(top_expr)))
row_labels <- gene_annotation_mapped$GeneSymbol[match(rownames(top_expr), rownames(expr_mapped))]
annotation_col <- data.frame(Disease = meta$disease, Region = meta$region)
rownames(annotation_col) <- colnames(top_expr)

pheatmap(top_scaled, show_rownames = TRUE, show_colnames = FALSE,
         annotation_col = annotation_col, labels_row = row_labels,
         main = "Top 50 DEGs (Active TB vs LTBI, Region-adjusted)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "heatmap.png", width = 10, height = 8)

cat("\n--- ANALYSIS COMPLETE ---\n")