## 03_GSE70947_QC_DE.R - STEPS 4-6: discovery QC, paired limma DE, RCD DE
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(pheatmap); library(matrixStats)
  library(ggrepel); library(patchwork); library(ggpubr)
})

run_wrapped("03", {
  hlog("=== 03 GSE70947 QC & DIFFERENTIAL EXPRESSION ===")
  emat <- readRDS(file.path(LOG, "cache", "GSE70947_std_expr.rds"))
  meta <- readRDS(file.path(LOG, "cache", "GSE70947_meta.rds"))
  sett <- readRDS(file.path(LOG, "cache", "rcd_sets.rds"))
  mypal <- c(normal = "#1b9e77", tumor = "#d95f02")

  ## ---------- missing data handling ----------
  miss_frac <- rowMeans(is.na(emat))
  keep_gene <- miss_frac <= 0.20
  hlog("genes with >20% missing excluded:", sum(!keep_gene))
  keynum("GSE70947_genes_excl_missing20", sum(!keep_gene))
  emat <- emat[keep_gene, , drop = FALSE]
  ## per-sample missingness
  samp_miss <- colMeans(is.na(emat))
  hlog("sample-wise missingness after gene filter: max", round(max(samp_miss), 4))
  ## median imputation (per gene)
  med <- apply(emat, 1, median, na.rm = TRUE)
  for (j in which(rowSums(is.na(emat)) > 0)) emat[j, is.na(emat[j, ])] <- med[j]
  hlog("imputation done; remaining NA:", sum(is.na(emat)))

  ## ---------- expression distribution ----------
  set.seed(SEED)
  ds <- data.frame(sample = rep(colnames(emat), each = 500),
                   condition = rep(meta$condition, each = 500))
  idx <- unlist(lapply(seq_len(ncol(emat)), function(j) sample(seq_len(nrow(emat)), 500)))
  ds$expr <- as.vector(emat[cbind(idx, rep(seq_len(ncol(emat)), each = 500))])
  p <- ggplot(ds, aes(x = sample, y = expr, fill = condition)) +
    geom_boxplot(outlier.size = 0.2, lwd = 0.25) +
    scale_fill_manual(values = mypal) +
    labs(x = "Sample", y = "Expression (log2)", title = "GSE70947 sample expression distributions") +
    theme_bw(base_size = 10) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                     legend.position = "top")
  save_gg(p, "GSE70947_expression_boxplot.png", w = 11, h = 5)

  set.seed(SEED)
  dens <- rbind(data.frame(condition = "normal",
                           expr = as.vector(emat[, meta$condition == "normal"])),
                data.frame(condition = "tumor",
                           expr = as.vector(emat[, meta$condition == "tumor"])))
  dens <- dens[!is.na(dens$expr), ]
  dens <- dens[sample(seq_len(nrow(dens)), min(1e5, nrow(dens))), ]
  p2 <- ggplot(dens, aes(x = expr, fill = condition)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = mypal) +
    labs(title = "GSE70947 expression density (sampled genes)", x = "Expression (log2)") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  save_gg(p2, "GSE70947_expression_density.png", w = 7, h = 5)

  ## ---------- sample correlation / clustering ----------
  cmat <- cor(emat, use = "pairwise.complete.obs")
  p3 <- pheatmap(cmat, show_rownames = FALSE, show_colnames = FALSE,
                 annotation_col = data.frame(condition = meta$condition, row.names = colnames(emat)),
                 annotation_colors = list(condition = mypal),
                 color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                 clustering_distance_cols = "euclidean", clustering_method = "average",
                 treeheight_col = 40, treeheight_row = 40, main = "GSE70947 sample correlation")
  save_base("GSE70947_sample_correlation_heatmap.png", w = 8, h = 7, print(p3))

  hcl <- hclust(as.dist(1 - cmat), method = "average")
  condv <- meta$condition[hcl$order]
  save_base("GSE70947_sample_clustering.png", w = 12, h = 6, {
    par(mar = c(4, 4, 2, 8), cex.axis = 0.7)
    plot(hcl, labels = FALSE, main = "GSE70947 hierarchical sample clustering (1 - correlation)", xlab = "", sub = "")
    cc <- ifelse(condv == "normal", "#1b9e77", "#d95f02")
    rect.hclust(hcl, k = 2, border = 3)
    legend("topright", legend = c("normal", "tumor"), col = c("#1b9e77", "#d95f02"), pch = 15, inset = -0.04, xpd = TRUE)
  })
## ---------- PCA ----------
  mv <- rowVars(emat, na.rm = TRUE)
  top <- order(mv, decreasing = TRUE)[seq_len(min(5000, nrow(emat)))]
  pca <- prcomp(t(emat[top, ]), center = TRUE, scale. = FALSE)
  pv <- summary(pca)$importance[2, 1:5] * 100
  pcd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    condition = meta$condition, sample = colnames(emat))
  p5 <- ggplot(pcd, aes(PC1, PC2, color = condition)) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_color_manual(values = mypal) +
    labs(title = "GSE70947 PCA (top 5,000 variable genes)",
         x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
    theme_bw(base_size = 12) + theme(legend.position = "top") +
    stat_ellipse(aes(group = condition), level = 0.9)
  save_gg(p5, "GSE70947_PCA_tumor_normal.png", w = 7, h = 6)
  saveRDS(pca, file.path(LOG, "cache", "GSE70947_pca.rds"))

  ## ---------- objective sample outlier detection ----------
  hv <- cor(emat[top, ], use = "pairwise.complete.obs")   ## sample cor from top variable genes
  mc <- rowMeans(hv, na.rm = TRUE)
  zk <- (mc - median(mc)) / (mad(mc) + 1e-12)
  out_samples <- names(zk)[zk < -5]
  rem <- data.frame(sample = out_samples,
                    reason = if (length(out_samples)) "inter-sample mean correlation z < -5 (low connectivity outlier)" else character(0),
                    criterion = if (length(out_samples)) "sample-network connectivity z-score < -5" else character(0))
  if (!nrow(rem)) rem <- data.frame(sample = NA, reason = "no samples removed", criterion = "sample-network connectivity z-score >= -5 for all samples")
  write_tab(rem, "GSE70947_removed_samples.csv")
  hlog("GSE70947 outlier samples (z<-5):", length(out_samples), if (length(out_samples)) paste(out_samples, collapse = ",") else "NONE")
  if (length(out_samples)) {
    keep_s <- !(colnames(emat) %in% out_samples)
    meta  <- meta[keep_s, ]
    emat  <- emat[, keep_s]
    hlog("removed", length(out_samples), "samples from downstream analysis")
  }

  ## retain only complete matched pairs for paired DE
  pat_tab <- table(meta$patient_id[meta$condition == "tumor"])
  pat_n   <- table(meta$patient_id[meta$condition == "normal"])
  complete <- intersect(names(pat_tab)[pat_tab == 1], names(pat_n)[pat_n == 1])
  use_s <- meta$sample[meta$patient_id %in% complete]
  hlog("complete matched pairs for paired DE:", length(complete), "pairs /",
       sum(meta$sample %in% use_s & meta$condition == "tumor"), "tumors")
  keynum("GSE70947_pairs_complete", length(complete))
## ==================== STEP 5: PAIRED DIFFERENTIAL EXPRESSION ====================
  E <- emat[, use_s]
  M <- meta[meta$sample %in% use_s, ]
  M$condition <- factor(M$condition, levels = c("normal", "tumor"))
  M$patient_id <- factor(M$patient_id)

  dupcor <- duplicateCorrelation(E, design = model.matrix(~ patient_id + condition, data = M),
                                 block = M$patient_id)
  hlog("consensus correlation (paired):", round(dupcor$consensus.correlation, 4))
  keynum("GSE70947_dupcor_consensus", round(dupcor$consensus.correlation, 4))

  design <- model.matrix(~ patient_id + condition, data = M)
  fit <- lmFit(E, design, block = M$patient_id, correlation = dupcor$consensus.correlation)
  fit <- eBayes(fit)
  coef_name <- tail(colnames(design), 1)   ## conditiontumor
  de_all <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  de_all <- de_all[order(de_all$adj.P.Val, -abs(de_all$logFC)), ]
  de_all$Gene <- rownames(de_all)
  de_all <- de_all[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  keynum("GSE70947_n_genes_DE", nrow(de_all))

  degs05 <- de_all[de_all$adj.P.Val < 0.05 & abs(de_all$logFC) >= 0.5, ]
  degs1  <- de_all[de_all$adj.P.Val < 0.05 & abs(de_all$logFC) >= 1, ]
  up <- degs05[degs05$logFC > 0, ]; dn <- degs05[degs05$logFC < 0, ]
  hlog("DEGs (FDR<0.05, |logFC|>=0.5):", nrow(degs05), " up:", nrow(up), " down:", nrow(dn))
  keynum("GSE70947_DEG_total", nrow(degs05)); keynum("GSE70947_DEG_up", nrow(up)); keynum("GSE70947_DEG_down", nrow(dn))

  write_tab(de_all, "GSE70947_DE_all_genes.csv")
  write_tab(degs05, "GSE70947_DEGs_FDR005_logFC05.csv")
  write_tab(degs1,  "GSE70947_DEGs_FDR005_logFC1.csv")
  write_tab(up,     "GSE70947_upregulated_DEGs.csv")
  write_tab(dn,     "GSE70947_downregulated_DEGs.csv")
  de_sum <- data.frame(
    contrast = "tumor vs normal (paired, duplicateCorrelation)",
    ntotal = length(complete), n_pairs = length(complete),
    n_total_genes = nrow(de_all), n_DEG_FDR05_logFC05 = nrow(degs05),
    n_up = nrow(up), n_down = nrow(dn),
    n_DEG_FDR05_logFC1 = nrow(degs1),
    consensus_correlation = dupcor$consensus.correlation,
    stringsAsFactors = FALSE)
  write_tab(de_sum, "GSE70947_DE_summary.csv")

  ## ---------- volcano ----------
  vd <- de_all
  vd$status <- ifelse(vd$adj.P.Val < 0.05 & vd$logFC > 0.5, "Upregulated",
               ifelse(vd$adj.P.Val < 0.05 & vd$logFC < -0.5, "Downregulated", "Not significant"))
  vd$neglog10p <- -log10(vd$P.Value)
  lab <- vd[vd$status != "Not significant" & vd$neglog10p > 10, ]
  if (nrow(lab) > 20) lab <- lab[order(lab$neglog10p, decreasing = TRUE), ][1:20, ]
  pv <- ggplot(vd, aes(logFC, neglog10p, color = status)) +
    geom_point(size = 0.55, alpha = 0.7) +
    scale_color_manual(values = c(Upregulated = "#d95f02", Downregulated = "#1b9e77", `Not significant` = "grey70")) +
    geom_vline(xintercept = c(-0.5, 0.5), lty = 2, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey40") +
    geom_text_repel(data = lab, aes(label = Gene), size = 2.4, max.overlaps = 30,
                    segment.size = 0.2) +
    labs(title = "GSE70947 tumor vs normal (paired limma)",
         x = "log2 fold change", y = "-log10(raw P)") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  save_gg(pv, "GSE70947_volcano_plot.png", w = 8, h = 7)

  ## ---------- DEG heatmaps ----------
  ## z-score scaled across samples
  zmat <- t(scale(t(E)))

  ann <- data.frame(condition = M$condition, row.names = colnames(E))
  anncol <- list(condition = mypal)

  if (nrow(degs05) >= 2 && nrow(degs05) <= 3000) {
    genes_hm <- rownames(degs05)
    hm_mat <- zmat[genes_hm, ]
    hm_mat[is.na(hm_mat)] <- 0
    png(file.path(FIG, "GSE70947_DEG_heatmap.png"), width = 9, height = 8, units = "in", res = 600)
    pheatmap(hm_mat, annotation_col = ann, annotation_colors = anncol,
             show_rownames = FALSE, show_colnames = FALSE,
             color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
             main = paste("GSE70947 DEGs (FDR<0.05, |logFC|>=0.5):", nrow(hm_mat), "genes"),
             treeheight_col = 25, treeheight_row = 5)
    dev.off(); hlog("saved figure GSE70947_DEG_heatmap.png")
  } else hlog("DEG heatmap skipped: too many/too few genes", nrow(degs05))

  top30 <- head(degs05, 30)
  hm2 <- zmat[rownames(top30), ]
  hm2[is.na(hm2)] <- 0
  png(file.path(FIG, "GSE70947_top_DEGs_heatmap.png"), width = 9, height = 7, units = "in", res = 600)
  pheatmap(hm2, annotation_col = ann, annotation_colors = anncol,
           show_rownames = TRUE, show_colnames = FALSE,
           fontsize_row = 6.5,
           color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           main = "GSE70947 top 30 DEGs")
  dev.off(); hlog("saved figure GSE70947_top_DEGs_heatmap.png")

  saveRDS(list(E = E, M = M, de_all = de_all, degs05 = degs05, degs1 = degs1,
               zmat = zmat, emat_all = emat, meta_all = meta),
          file.path(LOG, "cache", "GSE70947_DE.rds"))
## ==================== STEP 6: RCD DIFFERENTIAL EXPRESSION ====================
  pathway_list <- list(Ferroptosis = sett$Ferroptosis, Necroptosis = sett$Necroptosis,
                       Cuproptosis = sett$Cuproptosis)
  de_rc_tabs <- list(); sig_rc <- list()
  for (nm in names(pathway_list)) {
    genes <- intersect(pathway_list[[nm]], de_all$Gene)
    sub <- de_all[de_all$Gene %in% genes, ]
    sub$Pathway <- nm
    sub$DE_significance <- ifelse(sub$adj.P.Val < 0.05, "significant", "not significant")
    sub$DE_direction <- ifelse(sub$logFC > 0, "up",
                        ifelse(sub$logFC < 0, "down", "no change"))
    sub <- sub[, c("Gene", "Pathway", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
                   "DE_direction", "DE_significance")]
    de_rc_tabs[[nm]] <- sub
    sig_rc[[nm]] <- sub[sub$DE_significance == "significant" & abs(sub$logFC) >= 0.5, ]
    hlog(nm, ": measured", nrow(sub), " significant DEGs:", nrow(sig_rc[[nm]]))
    keynum(paste0("GSE70947_", tolower(nm), "_DEGs"), nrow(sig_rc[[nm]]))
    write_tab(sub, paste0("GSE70947_", tolower(nm), "_DE_results.csv"))
  }
  rcd_all <- do.call(rbind, de_rc_tabs)
  rcd_sig <- do.call(rbind, sig_rc)
  write_tab(rcd_sig, "GSE70947_significant_RCD_DEGs.csv")
  keynum("GSE70947_significant_RCD_DEGs", nrow(rcd_sig))

  ## combined RCD DEG matrix for plotting (z-scores)
  if (nrow(rcd_sig) >= 2) {
    genes_rc <- unique(rcd_sig$Gene)
    hm_mat <- zmat[genes_rc, ]            ## zmat computed in STEP 5 block
    hm_mat[is.na(hm_mat)] <- 0
    membership <- ifelse(genes_rc %in% sett$Ferroptosis, "Ferroptosis", "")
    membership[genes_rc %in% sett$Necroptosis] <- paste(membership[genes_rc %in% sett$Necroptosis], "Necroptosis")
    membership[genes_rc %in% sett$Cuproptosis] <- paste(membership[genes_rc %in% sett$Cuproptosis], "Cuproptosis")
    membership <- trimws(gsub("^\\s+|\\s+$", "", membership))
    ann2 <- data.frame(row.names = rownames(hm_mat), RCD_pathway = membership)
    pal_base <- c(Ferroptosis = "#66c2a5", Necroptosis = "#fc8d62", Cuproptosis = "#e78ac3",
                 "Ferroptosis Necroptosis" = "#a6d854", "Ferroptosis Cuproptosis" = "#ffd92f",
                 "Necroptosis Cuproptosis" = "#8da0cb", "Ferroptosis Necroptosis Cuproptosis" = "#984ea3")
    anncols2 <- list(RCD_pathway = pal_base[unique(as.character(ann2$RCD_pathway))])
    anncols2$RCD_pathway <- anncols2$RCD_pathway[!is.na(names(anncols2$RCD_pathway))]
    png(file.path(FIG, "GSE70947_RCD_DEG_heatmap.png"), width = 9, height = 8, units = "in", res = 600)
    pheatmap(hm_mat, annotation_col = ann, annotation_row = ann2,
             annotation_colors = c(list(condition = mypal), anncols2),
             show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 5,
             color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
             main = paste("GSE70947 significant RCD DEGs (", nrow(hm_mat), "genes )"),
             treeheight_col = 20, treeheight_row = 10)
    dev.off(); hlog("saved figure GSE70947_RCD_DEG_heatmap.png")

    ## logFC barplot by pathway
    bp <- rcd_sig[order(rcd_sig$Pathway, rcd_sig$logFC), ]
    ## a multi-pathway gene would otherwise create duplicated factor levels; keep
    ## one row per gene using its most extreme entry
    bp <- bp[order(bp$Pathway, -abs(bp$logFC)), ]
    bp <- bp[!duplicated(bp$Gene), ]
    bp <- bp[order(bp$Pathway, bp$logFC), ]
    bp$Gene <- factor(bp$Gene, levels = bp$Gene)
    pb <- ggplot(bp, aes(Gene, logFC, fill = Pathway)) +
      geom_col(width = 0.8) +
      facet_wrap(~ Pathway, scales = "free_x", ncol = 1) +
      coord_flip() +
      scale_fill_manual(values = c(Ferroptosis = "#e41a1c", Necroptosis = "#377eb8", Cuproptosis = "#4daf4a")) +
      labs(title = "GSE70947 RCD DEG log2 fold changes", y = "log2 FC (tumor vs normal)") +
      theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 0),
                                       strip.text = element_text(face = "bold"), legend.position = "none")
    save_gg(pb, "GSE70947_RCD_logFC_barplot.png", w = 8, h = 9)
  } else hlog("no significant RCD DEGs to plot")

  saveRDS(rcd_all, file.path(LOG, "cache", "GSE70947_RCD_DE_results.rds"))
  hlog("STEP 6 complete")
})