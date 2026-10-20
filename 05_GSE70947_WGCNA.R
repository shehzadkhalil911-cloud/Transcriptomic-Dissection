## 05_GSE70947_WGCNA.R - STEPS 8-10: signed WGCNA on GSE70947 tumors only
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({
  library(WGCNA); library(limma); library(ggplot2); library(matrixStats)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

run_wrapped("05", {
  hlog("=== 05 WGCNA (GSE70947 tumors) ===")
  emat <- readRDS(file.path(LOG, "cache", "GSE70947_std_expr.rds"))
  meta <- readRDS(file.path(LOG, "cache", "GSE70947_meta.rds"))
  tscores_obj <- readRDS(file.path(LOG, "cache", "GSE70947_tumor_scores.rds"))
  tscores <- tscores_obj$tscores; tumors <- tscores_obj$tumors

  ## ---------- tumor-only expression, NA filter, impute ----------
  E <- emat[, tumors]
  miss <- rowMeans(is.na(E))
  keep_gene <- miss <= 0.20
  E <- E[keep_gene, , drop = FALSE]
  for (j in which(rowSums(is.na(E)) > 0)) E[j, is.na(E[j, ])] <- median(E[j, ], na.rm = TRUE)
  hlog("tumor genes retained after NA filter:", nrow(E))

  ## ---------- informative gene selection: top 75% by MAD ----------
  madv <- rowMads(E)
  thr <- quantile(madv, probs = 0.25)
  keep <- madv >= thr
  datExpr <- t(E[keep, ])   # samples x genes
  hlog("MAD filter: kept", sum(keep), "of", nrow(E),
       "genes (top 75% by MAD, threshold", round(thr, 4), ")")
  keynum("GSE70947_WGCNA_genes_before_MAD", nrow(E))
  keynum("GSE70947_WGCNA_genes_used", sum(keep))

  ## ---------- gene filtering summary ----------
  filt <- data.frame(
    stage = c("standardized", "after_NA20", "after_MAD_top75"),
    n_genes = c(nrow(emat), nrow(E), sum(keep)),
    stringsAsFactors = FALSE)
  write_tab(filt, "GSE70947_WGCNA_gene_filtering_summary.csv")
  write_tab(data.frame(Gene = colnames(datExpr)), "GSE70947_WGCNA_genes_used.csv")

  ## ---------- goodSamplesGenes ----------
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    hlog("goodSamplesGenes: removing", sum(!gsg$goodGenes), "genes and",
         sum(!gsg$goodSamples), "samples")
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  ## ---------- tumor sample clustering / outliers ----------
  sampleTree <- hclust(dist(datExpr), method = "average")
  save_base("GSE70947_WGCNA_sample_clustering.png", w = 12, h = 6, {
    par(cex = 0.7, mar = c(5, 4, 2, 2))
    plot(sampleTree, main = "GSE70947 tumor samples, sample clustering (distance)",
         sub = "", xlab = "", cex.lab = 1.2, cex.axis = 0.8)
    abline(h = quantile(sampleTree$height, 0.999), col = "red", lty = 2)
  })
  ## objective outlier: connectivity z-score
  mcor <- cor(datExpr, use = "p"); km <- rowMeans(mcor)
  zk <- (km - median(km)) / (mad(km) + 1e-12)
  outw <- names(zk)[zk < -5]
  if (length(outw)) {
    write_tab(data.frame(sample = outw, reason = "WGCNA sample connectivity z < -5",
                         criterion = "network connectivity z-score"), "GSE70947_WGCNA_removed_samples.csv")
    datExpr <- datExpr[!(rownames(datExpr) %in% outw), ]
    hlog("WGCNA removed outlier samples:", paste(outw, collapse = ","))
  } else {
    write_tab(data.frame(sample = NA, reason = "no WGCNA samples removed",
                         criterion = "connectivity z >= -5 for all"), "GSE70947_WGCNA_removed_samples.csv")
    hlog("WGCNA no sample outliers removed")
  }
  tumors_final <- rownames(datExpr)
  keynum("GSE70947_WGCNA_n_tumors", length(tumors_final))
## ---------- soft threshold selection ----------
  powers <- c(seq(1, 20), 22, 24, 26, 28, 30)
  st <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 3,
                          networkType = "signed", corFnc = "bicor")
  save_base("GSE70947_WGCNA_soft_threshold.png", w = 10, h = 5, {
    par(mfrow = c(1, 2), cex = 0.8, mar = c(5, 5, 3, 2))
    fit <- -sign(st$fitIndices[, 3]) * st$fitIndices[, 2]
    plot(st$fitIndices[, 1], fit, xlab = "Soft threshold power", ylab = "Scale-free fit index (R^2)", type = "n",
         main = "Scale-free topology fit")
    text(st$fitIndices[, 1], fit, labels = powers, col = "red", cex = 0.9)
    abline(h = 0.85, col = "blue", lty = 2)
    plot(st$fitIndices[, 1], st$fitIndices[, 5], xlab = "Soft threshold power", ylab = "Mean connectivity",
         type = "n", main = "Mean connectivity")
    text(st$fitIndices[, 1], st$fitIndices[, 5], labels = powers, col = "red", cex = 0.9)
  })
  sft_stats <- data.frame(power = st$fitIndices[, 1],
                          SFT_Rsq = -sign(st$fitIndices[, 3]) * st$fitIndices[, 2],
                          mean_k = st$fitIndices[, 5],
                          stringsAsFactors = FALSE)
  write_tab(sft_stats, "GSE70947_WGCNA_soft_threshold_stats.csv")

  ## choose the smallest power >= 6 achieving SFT R^2 >= 0.85;
  ## fallback to max R^2 if none reaches the cut-off
  cand <- sft_stats$power[sft_stats$SFT_Rsq >= 0.85 & sft_stats$power >= 6]
  power <- if (length(cand)) min(cand) else sft_stats$power[which.max(sft_stats$SFT_Rsq)]
  hlog("chosen soft threshold power:", power)
  keynum("GSE70947_WGCNA_soft_power", power)

  ## ---------- network construction (blockwise, signed) ----------
  net <- blockwiseModules(datExpr, power = power, networkType = "signed",
                         TOMType = "signed", corType = "bicor",
                         minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25,
                         detectCutHeight = 0.995, pamRespectsDendro = FALSE,
                         numericLabels = TRUE, saveTOMFiles = FALSE,
                         verbose = 3, maxBlockSize = 5000,
                         randomSeed = SEED)
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  MEs0 <- net$MEs
  hlog("modules:", length(unique(moduleLabels)) - 1, "; grey(unassigned):", sum(moduleLabels == 0))
  keynum("GSE70947_WGCNA_n_modules", length(unique(moduleLabels)) - 1)
  keynum("GSE70947_WGCNA_n_genes_network", length(moduleLabels))
## ---------- module dendrogram ----------
  save_base("GSE70947_WGCNA_gene_dendrogram_modules.png", w = 12, h = 7, {
    par(cex = 0.6)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = paste0("GSE70947 tumor gene dendrogram and module colors (power=", power, ")"))
  })

  ## ---------- module assignments and sizes ----------
  genes_net <- colnames(datExpr)
  mdf <- data.frame(Gene = genes_net,
                    Module = moduleColors,
                    ModuleColor = moduleColors,
                    stringsAsFactors = FALSE)
  write_tab(mdf, "GSE70947_WGCNA_module_assignments.csv")
  msizes <- sort(table(moduleColors), decreasing = TRUE)
  msize_df <- data.frame(ModuleColor = names(msizes), ModuleSize = as.integer(msizes),
                         stringsAsFactors = FALSE)
  write_tab(msize_df, "GSE70947_WGCNA_module_sizes.csv")
  hlog("module sizes:", paste(names(msizes)[1:min(8, length(msizes))], msizes[1:min(8, length(msizes))], collapse = "; "), "...")

  ## ---------- module eigengenes ----------
  MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  ## rename ME columns by color
  colnames(MEs) <- gsub("^ME", "", colnames(MEs))
  save_base("GSE70947_WGCNA_module_eigengene_dendrogram.png", w = 10, h = 6, {
    par(cex = 0.7)
    plotEigengeneNetworks(MEs, "GSE70947 module eigengene clustering",
                          marDendro = c(3, 3, 2, 6), marHeatmap = c(3, 3, 2, 4),
                          plotDendrograms = TRUE, plotHeatmaps = FALSE)
  })

  ## ---------- module-trait analysis ----------
  ## traits: pathway scores + available clinical vars
  tt2 <- tscores[tumors_final, , drop = FALSE]
  colnames(tt2) <- paste0("score_", colnames(tt2))
  met_t <- meta[meta$sample %in% tumors_final, ]
  rownames(met_t) <- met_t$sample
  tit <- data.frame(ER = as.numeric(met_t[tumors_final, "er"]),
                    HER2 = ifelse(met_t[tumors_final, "her2"] == "pos", 1, ifelse(met_t[tumors_final, "her2"] == "neg", 0, NA)),
                    Age = as.numeric(met_t[tumors_final, "age"]),
                    BMI = as.numeric(met_t[tumors_final, "bmi"]),
                    Tumor_size_cm = as.numeric(met_t[tumors_final, "tumor_size_cm"]),
                    row.names = tumors_final)
  traits <- cbind(tt2, tit)
  traits <- as.data.frame(lapply(traits, as.numeric))
  ## trait availability
  hlog("trait NA counts:", paste(names(traits), colSums(is.na(traits)), collapse = "; "))

  nSamples <- nrow(datExpr)
  mtCor <- apply(MEs, 2, function(me) apply(traits, 2, function(tr) {
    ok <- !is.na(tr)
    if (sum(ok) < 10) return(NA_real_)
    suppressWarnings(cor(me[ok], tr[ok], method = "pearson"))
  }))
  mtP <- apply(MEs, 2, function(me) apply(traits, 2, function(tr) {
    ok <- !is.na(tr)
    if (sum(ok) < 10) return(NA_real_)
    ct <- suppressWarnings(cor.test(me[ok], tr[ok], method = "pearson"))
    ct$p.value
  }))
  mtCor <- t(mtCor); mtP <- t(mtP)
  mtFDR <- mtP
  for (jj in seq_len(ncol(mtFDR))) mtFDR[, jj] <- p.adjust(mtP[, jj], "BH")
  mtCor_df <- data.frame(Trait = rep(colnames(traits), each = nrow(mtCor)),
                         Module = rep(rownames(mtCor), ncol(traits)),
                         Correlation = as.vector(mtCor),
                         P_value = as.vector(mtP),
                         FDR = as.vector(mtFDR),
                         stringsAsFactors = FALSE)
  write_tab(mtCor_df, "GSE70947_WGCNA_module_trait_correlations.csv")
  write_tab(mtCor_df[, c("Trait", "Module", "FDR", "P_value", "Correlation")],
            "GSE70947_WGCNA_module_trait_FDR.csv")

  ## heatmap
  textM <- paste0(signif(mtCor, 2), "\n(", signif(mtFDR, 1), ")")
  dim(textM) <- dim(mtCor)
  save_base("GSE70947_WGCNA_module_trait_heatmap.png", w = 9, h = 7, {
    par(cex = 0.8, mar = c(9, 9, 2, 2))
    labeledHeatmap(Matrix = mtCor, xLabels = colnames(traits), yLabels = rownames(mtCor),
                   ySymbols = rownames(mtCor), colorLabels = FALSE, colors = blueWhiteRed(50),
                   textMatrix = textM, setStdMargins = FALSE, cex.text = 0.7,
                   main = "GSE70947 module-trait correlations (FDR in parentheses)")
  })

  ## ---------- save all WGCNA objects ----------
  params_df <- data.frame(
    parameter = c("networkType", "TOMType", "corType", "power", "minModuleSize",
                  "deepSplit", "mergeCutHeight", "detectCutHeight", "pamRespectsDendro",
                  "maxBlockSize", "MAD_top_quantile", "MAD_threshold", "n_genes_used",
                  "n_tumors_used", "n_modules", "n_grey"),
    value = c("signed", "signed", "bicor", power, 30, 2, 0.25, 0.995, FALSE, 5000,
              0.75, round(thr, 4), ncol(datExpr), nrow(datExpr),
              length(unique(moduleLabels)) - 1, sum(moduleLabels == 0)),
    stringsAsFactors = FALSE)
  write_tab(params_df, "GSE70947_WGCNA_parameters.csv")

  saveRDS(list(net = net, datExpr = datExpr, moduleLabels = moduleLabels,
               moduleColors = moduleColors, MEs = MEs,
               power = power, traits = traits, mtCor = mtCor, mtP = mtP, mtFDR = mtFDR,
               tumors_final = tumors_final, mdf = mdf, thrMAD = thr),
          file.path(LOG, "cache", "GSE70947_WGCNA.rds"))
  hlog("STEP 8-10 complete")
})