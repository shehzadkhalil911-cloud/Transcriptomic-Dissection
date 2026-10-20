## 07_validation_datasets.R - STEPS 16-18,20: process & validate GSE42568, GSE15852, GSE29431
## DE per dataset (paired for GSE15852), candidate validation, pathway score validation
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(matrixStats); library(pheatmap); library(GSVA)
})

run_wrapped("07", {
  hlog("=== 07 VALIDATION DATASETS ===")
  frozen <- readRDS(file.path(LOG, "cache", "frozen_candidates.rds"))
  de_disc <- readRDS(file.path(LOG, "cache", "GSE70947_DE.rds"))$de_all
  sett <- readRDS(file.path(LOG, "cache", "rcd_sets.rds"))
  gsets <- list(Ferroptosis = sett$Ferroptosis, Necroptosis = sett$Necroptosis, Cuproptosis = sett$Cuproptosis)

  valid_gses <- c("GSE42568", "GSE15852", "GSE29431")
  cand_valid <- list(); path_summary <- list(); val_objs <- list()

  for (g in valid_gses) {
    hlog("================ validation dataset:", g, "================")
    emat <- readRDS(file.path(LOG, "cache", paste0(g, "_std_expr.rds")))
    meta <- readRDS(file.path(LOG, "cache", paste0(g, "_meta.rds")))
    paired <- if (g == "GSE15852") TRUE else FALSE

    ## ---------- QC: PCA ----------
    miss <- rowMeans(is.na(emat))
    e0 <- emat[miss <= 0.20, ]
    for (j in which(rowSums(is.na(e0)) > 0)) e0[j, is.na(e0[j, ])] <- median(e0[j, ], na.rm = TRUE)
    hlog(g, "genes retained for analysis (missing <=20%):", nrow(e0))
    mv <- rowVars(e0)
    top <- order(mv, decreasing = TRUE)[seq_len(min(2000, nrow(e0)))]
    pca <- tryCatch(prcomp(t(e0[top, ]), center = TRUE), error = function(e) NULL)
    if (!is.null(pca)) {
      pv <- summary(pca)$importance[2, 1:2] * 100
      pcd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], condition = meta$condition)
      pp <- ggplot(pcd, aes(PC1, PC2, color = condition)) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = c(normal = "#1b9e77", tumor = "#d95f02")) +
        labs(title = paste(g, "PCA"), x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
        theme_bw(base_size = 11) + theme(legend.position = "top")
      save_gg(pp, paste0(g, "_PCA.png"), w = 6.5, h = 5.5)
    }

    ## ---------- DE ----------
    E <- e0
    M <- meta
    if (paired) {
      M$condition <- factor(M$condition, levels = c("normal", "tumor"))
      M$patient_id <- factor(M$patient_id)
      dc <- duplicateCorrelation(E, design = model.matrix(~ patient_id + condition, data = M),
                                 block = M$patient_id)
      design <- model.matrix(~ patient_id + condition, data = M)
      fit <- lmFit(E, design, block = M$patient_id, correlation = dc$consensus.correlation)
      fit <- eBayes(fit)
      coef_name <- tail(colnames(design), 1)
      hlog(g, "paired DE with duplicateCorrelation; consensus cor =", round(dc$consensus.correlation, 4))
    } else {
      M$condition <- factor(M$condition, levels = c("normal", "tumor"))
      design <- model.matrix(~ condition, data = M)
      fit <- lmFit(E, design)
      fit <- eBayes(fit)
      coef_name <- "conditiontumor"
    }
    de_all <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    de_all$Gene <- rownames(de_all)
    de_all <- de_all[order(de_all$adj.P.Val, -abs(de_all$logFC)), c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    write_tab(de_all, paste0(g, "_DE_all_genes.csv"))
    hlog(g, "DE done; genes:", nrow(de_all), " significant (FDR<0.05):", sum(de_all$adj.P.Val < 0.05))

    ## ---------- candidate validation ----------
    fc <- de_all[de_all$Gene %in% frozen$Gene, ]
    fc$Pathway <- frozen$`RCD_pathway(s)`[match(fc$Gene, frozen$Gene)]
    fc$Discovery_logFC <- de_disc$logFC[match(fc$Gene, de_disc$Gene)]
    fc$Discovery_FDR   <- de_disc$adj.P.Val[match(fc$Gene, de_disc$Gene)]
    fc$Direction <- ifelse(fc$logFC > 0, "up", ifelse(fc$logFC < 0, "down", "nc"))
    fc$Direction_concordant_with_discovery <- sign(fc$logFC) == sign(fc$Discovery_logFC) & sign(fc$logFC) != 0
    fc$P_sig <- fc$P.Value < 0.05
    fc$FDR_sig <- fc$adj.P.Val < 0.05
    fc <- fc[, c("Gene", "Pathway", "Discovery_logFC", "Discovery_FDR", "logFC", "t",
                 "P.Value", "adj.P.Val", "Direction", "Direction_concordant_with_discovery")]
    names(fc)[names(fc) == "logFC"] <- "Validation_logFC"
    names(fc)[names(fc) == "adj.P.Val"] <- "Validation_FDR"
    cand_valid[[g]] <- fc
    write_tab(fc, paste0(g, "_candidate_validation.csv"))
    hlog(g, "candidate validation rows:", nrow(fc))
## ---------- pathway scores (ssgsea) ----------
    eimat <- e0   ## already NA-filtered and imputed
    gs_use <- lapply(gsets, function(gs) intersect(gs, rownames(eimat)))
    hlog(g, "ssgsea gene coverage:", sapply(gs_use, length))
    scores_raw <- gsva(ssgseaParam(exprData = as.matrix(eimat), geneSets = gs_use,
                              normalize = TRUE), verbose = FALSE)
    scores <- t(scores_raw)   ## samples x pathways
    sc_df <- data.frame(sample = rownames(scores), scores, stringsAsFactors = FALSE)
    sc_df <- merge(sc_df, meta[, c("sample", "condition", "patient_id")], by = "sample", sort = FALSE)
    write_tab(sc_df, paste0(g, "_RCD_pathway_scores.csv"))

    ## compare tumor vs normal
    comp_rows <- list()
    for (nm in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
      y <- scores[, nm]
      dft <- data.frame(y = y, condition = meta$condition, patient_id = meta$patient_id)
      tT <- dft$y[dft$condition == "tumor"]; tN <- dft$y[dft$condition == "normal"]
      if (paired) {
        mm <- merge(subset(dft, condition == "tumor"), subset(dft, condition == "normal"),
                    by = "patient_id", suffixes = c("_T", "_N"))
        pw <- suppressWarnings(wilcox.test(mm$y_T, mm$y_N, paired = TRUE)$p.value)
        eff <- mean(mm$y_T - mm$y_N) / (sd(mm$y_T - mm$y_N) + 1e-12)
      } else {
        pw <- suppressWarnings(wilcox.test(tT, tN)$p.value)
        eff <- (mean(tT) - mean(tN)) / sqrt((var(tT) + var(tN)) / 2 + 1e-12)
      }
      comp_rows[[nm]] <- data.frame(Dataset = g, Pathway = nm, n_tumor = length(tT), n_normal = length(tN),
                                    Mean_tumor = mean(tT), Mean_normal = mean(tN),
                                    Effect_size = eff, Wilcox_P = pw, stringsAsFactors = FALSE)
    }
    path_summary[[g]] <- do.call(rbind, comp_rows)
    val_objs[[g]] <- list(de_all = de_all, scores = scores, sc_df = sc_df)
  }

  path_comb <- do.call(rbind, path_summary)
  path_comb$Wilcox_FDR <- p.adjust(path_comb$Wilcox_P, "BH")
  write_tab(path_comb, "RCD_pathway_score_validation_summary.csv")
  hlog("pathway score validation summary written")

  ## ---------- STEP 18: cross-dataset validation matrix ----------
  allg <- unique(frozen$Gene)
  cvm <- data.frame(Gene = allg, stringsAsFactors = FALSE)
  cvm$RCD_Pathway <- frozen$`RCD_pathway(s)`[match(allg, frozen$Gene)]
  cvm$Discovery_logFC <- de_disc$logFC[match(allg, de_disc$Gene)]
  cvm$Discovery_FDR   <- de_disc$adj.P.Val[match(allg, de_disc$Gene)]
  cvm$N_datasets_measured <- 0L; cvm$N_same_direction <- 0L; cvm$N_FDR_significant <- 0L; cvm$N_nominal_sig <- 0L
  for (g in valid_gses) {
    f <- cand_valid[[g]]
    cvm[[paste0(g, "_logFC")]] <- f$Validation_logFC[match(allg, f$Gene)]
    cvm[[paste0(g, "_FDR")]]    <- f$Validation_FDR[match(allg, f$Gene)]
    cvm[[paste0(g, "_P")]]      <- f$P.Value[match(allg, f$Gene)]
    cvm[[paste0(g, "_direction_concordant")]] <- f$Direction_concordant_with_discovery[match(allg, f$Gene)]
    meas <- !is.na(cvm[[paste0(g, "_logFC")]])
    cvm$N_datasets_measured <- cvm$N_datasets_measured + meas
    same <- meas & !is.na(cvm[[paste0(g, "_direction_concordant")]]) & cvm[[paste0(g, "_direction_concordant")]]
    cvm$N_same_direction <- cvm$N_same_direction + same
    cvm$N_FDR_significant <- cvm$N_FDR_significant + (!is.na(cvm[[paste0(g, "_FDR")]]) & cvm[[paste0(g, "_FDR")]] < 0.05)
    cvm$N_nominal_sig <- cvm$N_nominal_sig + (!is.na(cvm[[paste0(g, "_P")]]) & cvm[[paste0(g, "_P")]] < 0.05)
  }
  ## classification
  status <- rep("NO_EVIDENCE", nrow(cvm))
  m <- cvm$N_datasets_measured
  sd <- cvm$N_same_direction
  fdr <- cvm$N_FDR_significant
  nom <- cvm$N_nominal_sig
  robust <- m == 3 & sd == 3 & fdr >= 2
  ## support: dominant same direction (>=2 of 3, or 2 of 2) and at least one cohort
  ## with statistical evidence (nominal or FDR), no more than one cohort opposing
  supported <- m >= 2 & sd >= min(2, m) & (fdr > 0 | nom > 0) & (m - sd) <= 1 & !robust
  inconsistent <- m >= 2 & (m - sd) >= 2
  status[robust] <- "ROBUST"
  status[supported] <- "SUPPORTED"
  status[inconsistent] <- "INCONSISTENT"
  status[status == "NO_EVIDENCE" & m == 0] <- "NOT_MEASURED"
  cvm$Overall_validation_status <- status
  write_tab(cvm, "RCD_candidate_cross_dataset_validation.csv")
  hlog("validation status counts:", paste(names(table(status)), table(status), collapse = "; "))
  for (s in c("ROBUST", "SUPPORTED", "INCONSISTENT")) keynum(paste0("validation_", s), sum(status == s))
  saveRDS(list(cand_valid = cand_valid, val_objs = val_objs, path_comb = path_comb, cvm = cvm),
          file.path(LOG, "cache", "validation_results.rds"))
  hlog("saved validation_results.rds (early checkpoint)")
## ---------- STEP 18 figures ----------
  ## validation heatmap (signed -log10 FDR* direction) for candidates with data
  hm_genes <- cvm$Gene[(!is.na(cvm$GSE42568_logFC) | !is.na(cvm$GSE15852_logFC) | !is.na(cvm$GSE29431_logFC))]
  if (length(hm_genes)) {
    hmmat <- sapply(valid_gses, function(g) {
      lfc <- cvm[[paste0(g, "_logFC")]][match(hm_genes, cvm$Gene)]
      fdr <- cvm[[paste0(g, "_FDR")]][match(hm_genes, cvm$Gene)]
      v <- -log10(fdr + 1e-320) * sign(lfc)
      v[is.na(v)] <- 0
      v
    })
    rownames(hmmat) <- hm_genes; colnames(hmmat) <- valid_gses
    disc_v <- -log10(cvm$Discovery_FDR[match(hm_genes, cvm$Gene)] + 1e-320) * sign(cvm$Discovery_logFC[match(hm_genes, cvm$Gene)])
    hmmat <- cbind(Discovery = disc_v, hmmat)
    pw_raw <- cvm$RCD_Pathway[match(hm_genes, cvm$Gene)]
    pw_levels <- sort(unique(pw_raw))
    annrow <- data.frame(Pathway = factor(pw_raw, levels = pw_levels), row.names = hm_genes)
    pw_cols <- c(Ferroptosis = "#66c2a5", Necroptosis = "#fc8d62", Cuproptosis = "#e78ac3",
                 `Ferroptosis;Necroptosis` = "#a6d854", `Ferroptosis;Cuproptosis` = "#ffd92f",
                 `Necroptosis;Cuproptosis` = "#e5c494",
                 `Ferroptosis;Necroptosis;Cuproptosis` = "#b3b3b3")
    anncolr <- list(Pathway = pw_cols[pw_levels])
    hm_gap <- max(3, nrow(hmmat))
    hmmat2 <- hmmat
    hmmat2[hmmat2 == 0] <- NA
    png(file.path(FIG, "RCD_candidate_validation_heatmap.png"), width = 7, height = max(6, hm_gap * 0.14), units = "in", res = 600)
    pheatmap(hmmat2, annotation_row = annrow, annotation_colors = anncolr,
             cluster_cols = FALSE,
             color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
             na_col = "grey92",
             show_rownames = TRUE, show_colnames = TRUE, fontsize_row = 6,
             main = "RCD candidate validation: signed -log10(FDR) by cohort")
    dev.off(); hlog("saved figure RCD_candidate_validation_heatmap.png")
  }

  ## logFC across datasets (lines per gene)
  lfc_long <- do.call(rbind, lapply(valid_gses, function(g) data.frame(
    Gene = cvm$Gene, Dataset = g,
    logFC = cvm[[paste0(g, "_logFC")]],
    Pathway = cvm$RCD_Pathway, stringsAsFactors = FALSE)))
  lfc_long <- lfc_long[!is.na(lfc_long$logFC), ]
  lfc_long <- lfc_long[lfc_long$Gene %in% names(which(table(lfc_long$Gene) >= 2)), ]
  lfc_long$Pathway <- factor(lfc_long$Pathway, levels = sort(unique(lfc_long$Pathway)))
  pw_cols2 <- c(Ferroptosis = "#e41a1c", Necroptosis = "#377eb8", Cuproptosis = "#4daf4a",
                `Ferroptosis;Necroptosis` = "#a6d854", `Ferroptosis;Cuproptosis` = "#ffd92f",
                `Necroptosis;Cuproptosis` = "#e5c494",
                `Ferroptosis;Necroptosis;Cuproptosis` = "#b3b3b3")
  pl <- ggplot(lfc_long, aes(Dataset, logFC, group = Gene, color = Pathway)) +
    geom_line(alpha = 0.3) + geom_point(size = 0.8) +
    geom_hline(yintercept = 0, lty = 2, color = "grey30") +
    facet_wrap(~ Pathway, ncol = 1) +
    scale_color_manual(values = pw_cols2[levels(lfc_long$Pathway)]) +
    labs(title = "RCD candidate log2FC across validation cohorts", y = "log2FC (tumor vs normal)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")
  save_gg(pl, "RCD_candidate_logFC_across_datasets.png", w = 7, h = 11)

  ## ---------- STEP 20 figure: pathway score validation ----------
  pvl <- do.call(rbind, lapply(valid_gses, function(g) {
    s <- val_objs[[g]]$sc_df
    reshape2::melt(s, id.vars = c("sample", "condition"),
                   measure.vars = c("Ferroptosis", "Necroptosis", "Cuproptosis"),
                   variable.name = "Pathway", value.name = "score")
  }))
  pvl <- pvl[!is.na(pvl$condition), ]
  pv <- ggplot(pvl, aes(condition, score, fill = condition)) +
    geom_boxplot(alpha = 0.7, width = 0.55, outlier.size = 0.5) +
    facet_grid(Pathway ~ Dataset, scales = "free_y") +
    scale_fill_manual(values = c(normal = "#1b9e77", tumor = "#d95f02")) +
    labs(title = "RCD pathway score validation across cohorts", x = "", y = "ssGSEA score") +
    theme_bw(base_size = 10) + theme(legend.position = "top", strip.text = element_text(face = "bold"))
  save_gg(pv, "RCD_pathway_score_validation.png", w = 10, h = 9)

  saveRDS(list(cand_valid = cand_valid, val_objs = val_objs, path_comb = path_comb, cvm = cvm),
          file.path(LOG, "cache", "validation_results.rds"))
  hlog("STEP 07 complete")
})