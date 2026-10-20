## 08_meta_analysis.R - STEP 19: gene-level meta-analysis across validation cohorts
## Effect measure: unpaired Hedges' g (SMD) computed identically per cohort;
## fixed-effect and random-effects models with metafor.
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({ library(metafor); library(ggplot2) })

run_wrapped("08", {
  hlog("=== 08 META-ANALYSIS (validation cohorts) ===")
  frozen <- readRDS(file.path(LOG, "cache", "frozen_candidates.rds"))
  vres <- readRDS(file.path(LOG, "cache", "validation_results.rds"))
  valid_gses <- c("GSE42568", "GSE15852", "GSE29431")

  ## per-dataset per-gene SMD
  smd_list <- list()
  for (g in valid_gses) {
    emat <- readRDS(file.path(LOG, "cache", paste0(g, "_std_expr.rds")))
    meta <- readRDS(file.path(LOG, "cache", paste0(g, "_meta.rds")))
    genes <- intersect(frozen$Gene, rownames(emat))
    mt <- meta$sample[meta$condition == "tumor"]; mn <- meta$sample[meta$condition == "normal"]
    rows <- lapply(genes, function(ge) {
      t <- emat[ge, mt]; n <- emat[ge, mn]
      t <- t[!is.na(t)]; n <- n[!is.na(n)]
      if (length(t) < 3 || length(n) < 3) return(NULL)
      es <- tryCatch(escalc(measure = "SMD", m1i = mean(t), m2i = mean(n),
                            sd1i = sd(t), sd2i = sd(n), n1i = length(t), n2i = length(n)),
                     error = function(e) NULL)
      if (is.null(es)) return(NULL)
      data.frame(Gene = ge, Dataset = g, SMD = es$yi, SMD_var = es$vi,
                 nT = length(t), nN = length(n), stringsAsFactors = FALSE)
    })
    smd_list[[g]] <- do.call(rbind, rows)
  }
  smd_all <- do.call(rbind, smd_list)

  meta_rows <- list()
  genes_ok <- unique(smd_all$Gene)
  for (ge in genes_ok) {
    sub <- smd_all[smd_all$Gene == ge, ]
    if (nrow(sub) < 2) next
    for (mth in c("FE", "REML")) {
      fit <- tryCatch(rma(yi = SMD, vi = SMD_var, data = sub, method = if (mth == "FE") "FE" else "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) next
      meta_rows[[paste(ge, mth)]] <- data.frame(
        Gene = ge, Model = mth, Pooled_effect = fit$b[1, 1],
        CI_lower = fit$ci.lb, CI_upper = fit$ci.ub,
        Meta_P_value = fit$pval, Q = fit$QE, Q_P_value = fit$QEp,
        I2 = if (mth == "FE") NA else tryCatch(fit$I2, error = function(e) NA_real_),
        Tau2 = if (mth == "FE") 0 else fit$tau2,
        n_datasets = nrow(sub), stringsAsFactors = FALSE)
    }
  }
  meta_all <- do.call(rbind, meta_rows)
  meta_all$Meta_FDR <- p.adjust(meta_all$Meta_P_value, "BH")
  meta_all <- meta_all[order(meta_all$Gene, meta_all$Model), ]
  write_tab(meta_all, "RCD_candidate_validation_meta_analysis.csv")
  hlog("meta-analysis complete;", nrow(meta_all), "rows for", length(unique(meta_all$Gene)), "genes")
  keynum("meta_genes_analyzed", length(unique(meta_all$Gene)))

  ## ---------- forest plots for strongest genes ----------
  ## rank genes by REML pooled |effect| among those with FDR < 0.05 (FE or REML);
  ## fall back to top by |effect|
  rem <- meta_all[meta_all$Model == "REML", ]
  sig <- rem[rem$Meta_FDR < 0.05 & !is.na(rem$Meta_FDR), ]
  top_genes <- if (nrow(sig)) sig$Gene[order(-abs(sig$Pooled_effect))] else
    rem$Gene[order(-abs(rem$Pooled_effect))]
  top_genes <- head(unique(top_genes), 12)

  for (ge in top_genes) {
    sub <- smd_all[smd_all$Gene == ge, ]
    fitF <- tryCatch(rma(yi = SMD, vi = SMD_var, data = sub, method = "FE"), error = function(e) NULL)
    fitR <- tryCatch(rma(yi = SMD, vi = SMD_var, data = sub, method = "REML"), error = function(e) NULL)
    fitP <- if (!is.null(fitR)) fitR else fitF
    if (is.null(fitP)) next
    pool_lab <- if (!is.null(fitR)) "Pooled (random)" else "Pooled (fixed)"
    fp <- data.frame(
      label = c(as.character(sub$Dataset), pool_lab),
      est = c(sub$SMD, fitP$b[1, 1]),
      lo  = c(sub$SMD - 1.96 * sqrt(sub$SMD_var), fitP$ci.lb),
      hi  = c(sub$SMD + 1.96 * sqrt(sub$SMD_var), fitP$ci.ub),
      type = c(rep("study", nrow(sub)), "pooled"),
      stringsAsFactors = FALSE)
    fp$label <- factor(fp$label, levels = rev(fp$label))
    p <- ggplot(fp, aes(est, label)) +
      geom_vline(xintercept = 0, color = "grey40", lty = 2) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) +
      geom_point(aes(color = type, size = type)) +
      scale_color_manual(values = c(study = "#2c7fb8", pooled = "#d95f02")) +
      scale_size_manual(values = c(study = 2.5, pooled = 4)) +
      labs(title = paste("Meta-analysis forest plot:", ge),
           subtitle = sprintf("FE %s | REML I2=%s  P=%s (FDR %s)",
                              if (!is.null(fitF)) sprintf("est %.2f [%.2f, %.2f]", fitF$b[1,1], fitF$ci.lb, fitF$ci.ub) else "n/a",
                              if (!is.null(fitR)) round(tryCatch(fitR$I2, error = function(e) NA_real_), 0) else "n/a",
                              if (!is.null(fitP)) signif(fitP$pval, 3) else "n/a",
                              if (!is.null(fitP)) signif(p.adjust(fitP$pval, "BH"), 3) else "n/a"),
           x = "Standardized mean difference (Hedges' g)", y = "") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
    save_gg(p, paste0("meta_forest_", ge, ".png"), w = 7, h = 4.5)
  }
  hlog("forest plots for", length(top_genes), "genes")

  saveRDS(list(smd_all = smd_all, meta_all = meta_all), file.path(LOG, "cache", "meta_results.rds"))
  hlog("STEP 19 complete")
})