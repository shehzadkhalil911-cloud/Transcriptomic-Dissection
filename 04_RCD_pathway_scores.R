## 04_RCD_pathway_scores.R - STEP 7: GSVA/ssGSEA pathway activity in GSE70947
Sys.setenv(BRCA_NO_WGCNA_THREADS = "1")   ## avoid BLAS-thread conflict with GSVA
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({ library(GSVA); library(ggplot2); library(ggpubr) })

run_wrapped("04", {
  hlog("=== 04 RCD PATHWAY SCORES (GSE70947) ===")
  emat <- readRDS(file.path(LOG, "cache", "GSE70947_std_expr.rds"))
  meta <- readRDS(file.path(LOG, "cache", "GSE70947_meta.rds"))
  sett <- readRDS(file.path(LOG, "cache", "rcd_sets.rds"))

  gsets <- list(Ferroptosis = intersect(sett$Ferroptosis, rownames(emat)),
                Necroptosis = intersect(sett$Necroptosis, rownames(emat)),
                Cuproptosis = intersect(sett$Cuproptosis, rownames(emat)))

  ## ---------- gene coverage across all four datasets ----------
  cov_rows <- list()
  for (g in c("GSE70947", "GSE42568", "GSE15852", "GSE29431")) {
    e <- readRDS(file.path(LOG, "cache", paste0(g, "_std_expr.rds")))
    for (nm in names(gsets)) {
      cov_rows[[paste(g, nm)]] <- data.frame(
        Dataset = g, Pathway = nm,
        Pathway_size_original = length(gsets[[nm]]),   # placeholder corrected below
        Genes_measured_in_dataset = sum(gsets[[nm]] %in% rownames(e)),
        Fraction = round(sum(gsets[[nm]] %in% rownames(e)) / length(gsets[[nm]]), 3),
        stringsAsFactors = FALSE)
    }
  }
  orig <- c(Ferroptosis = length(sett$Ferroptosis), Necroptosis = length(sett$Necroptosis),
            Cuproptosis = length(sett$Cuproptosis))
  cov_rows2 <- do.call(rbind, cov_rows)
  cov_rows2$Pathway_size_original <- unname(orig[cov_rows2$Pathway])
  write_tab(cov_rows2, "RCD_pathway_gene_coverage_by_dataset.csv")
  hlog("pathway coverage table written")

  ## ---------- GSVA ssgsea for GSE70947 ----------
  ## genes with <=20% missing are retained and remaining NAs median-imputed
  miss <- rowMeans(is.na(emat))
  eimat <- emat[miss <= 0.20, ]
  hlog("genes retained (missing <=20%):", nrow(eimat))
  if (any(rowSums(is.na(eimat)) > 0)) {
    for (j in which(rowSums(is.na(eimat)) > 0))
      eimat[j, is.na(eimat[j, ])] <- median(eimat[j, ], na.rm = TRUE)
    hlog("median-imputed remaining NAs for ssgsea input")
  }
  hlog("genes with complete data for GSVA:", sum(!apply(is.na(eimat), 1, any)))
  scores <- gsva(ssgseaParam(exprData = as.matrix(eimat), geneSets = gsets,
                             normalize = TRUE), verbose = FALSE)
  scores <- t(scores)   # samples x pathways
  sc_df <- data.frame(sample = rownames(scores), scores, stringsAsFactors = FALSE)
  sc_df <- merge(sc_df, meta[, c("sample", "condition", "patient_id", "er", "her2", "age")],
                 by = "sample", sort = FALSE)
  write_tab(sc_df, "GSE70947_RCD_pathway_scores.csv")
  keynum("GSE70947_ssgsea_ferro_genes", length(gsets$Ferroptosis))
  keynum("GSE70947_ssgsea_necro_genes", length(gsets$Necroptosis))
  keynum("GSE70947_ssgsea_cupro_genes", length(gsets$Cuproptosis))

  ## ---------- paired comparisons tumor vs normal ----------
  comp_rows <- list()
  for (nm in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
    y <- scores[, nm]
    df <- data.frame(y = y, condition = meta$condition, patient_id = meta$patient_id)
    dfT <- df[df$condition == "tumor", ]; dfN <- df[df$condition == "normal", ]
    m <- merge(dfT, dfN, by = "patient_id", suffixes = c("_T", "_N"))
    d <- m$y_T - m$y_N
    pt <- tryCatch(t.test(d)$p.value, error = function(e) NA)
    pw <- tryCatch(wilcox.test(d)$p.value, error = function(e) NA)
    eff <- mean(d) / (sd(d) + 1e-12)
    comp_rows[[nm]] <- data.frame(Pathway = nm,
                                  Mean_tumor = mean(dfT$y), Mean_normal = mean(dfN$y),
                                  Mean_diff = mean(d), SD_diff = sd(d),
                                  Paired_cohens_d = eff, Paired_t_P = pt, Paired_wilcox_P = pw,
                                  stringsAsFactors = FALSE)
  }
  comp <- do.call(rbind, comp_rows)
  comp$Paired_t_FDR <- p.adjust(comp$Paired_t_P, "BH")
  comp$Paired_wilcox_FDR <- p.adjust(comp$Paired_wilcox_P, "BH")
  write_tab(comp, "GSE70947_RCD_pathway_score_comparison.csv")
  hlog("pathway score comparison:", apply(comp, 1, function(r) paste(r["Pathway"], "tc=", round(as.numeric(r["Paired_t_FDR"]), 4))))
  for (i in seq_len(nrow(comp))) keynum(paste0("GSE70947_score_", comp$Pathway[i], "_paired_t_FDR"), comp$Paired_t_FDR[i])

  ## ---------- figure: tumor vs normal with paired lines ----------
  plotdf <- reshape2::melt(sc_df, id.vars = c("sample", "condition", "patient_id"),
                           measure.vars = c("Ferroptosis", "Necroptosis", "Cuproptosis"),
                           variable.name = "Pathway", value.name = "score")
  pm <- merge(subset(plotdf, condition == "tumor"),
              subset(plotdf, condition == "normal"), by = c("patient_id", "Pathway"),
              suffixes = c("_T", "_N"))
  p <- ggplot(plotdf, aes(condition, score, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_line(data = pm, aes(patient_id, score_T, group = patient_id),
              inherit.aes = FALSE, color = "grey75", alpha = 0.5, show.legend = FALSE) +
    geom_point(data = pm, aes(patient_id, score_N), inherit.aes = FALSE, color = "grey40", size = 0.6) +
    facet_wrap(~ Pathway, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(normal = "#1b9e77", tumor = "#d95f02")) +
    labs(title = "GSE70947 RCD pathway scores: paired tumor vs normal", x = "", y = "ssGSEA score") +
    theme_bw(base_size = 12) + theme(legend.position = "top", strip.text = element_text(face = "bold"))
  labdf <- comp[, c("Pathway", "Paired_wilcox_P", "Paired_wilcox_FDR")]
  ymaxs <- aggregate(score ~ Pathway, plotdf, FUN = max, na.rm = TRUE)
  labdf$y <- ymaxs$score[match(labdf$Pathway, ymaxs$Pathway)] * 1.15
  p <- p +
    geom_text(data = labdf, aes(x = "tumor", y = y,
             label = sprintf("paired wilcox p=%.3g (FDR %.3g)", Paired_wilcox_P, Paired_wilcox_FDR)),
              size = 3, inherit.aes = FALSE, check_overlap = FALSE)
  save_gg(p, "GSE70947_RCD_pathway_scores_tumor_vs_normal.png", w = 7.5, h = 9)

  saveRDS(list(scores = scores, sc_df = sc_df, comp = comp),
          file.path(LOG, "cache", "GSE70947_scores.rds"))
  ## tumor-only score matrix for WGCNA traits
  tumors <- meta$sample[meta$condition == "tumor"]
  tscores <- scores[tumors, , drop = FALSE]
  saveRDS(list(tscores = tscores, tumors = tumors), file.path(LOG, "cache", "GSE70947_tumor_scores.rds"))
  hlog("STEP 7 complete")
})