## 12_final_integration.R - STEP 24: final evidence table & classifications
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(reshape2) })

run_wrapped("12", {
  hlog("=== 12 FINAL INTEGRATION ===")
  frozen <- readRDS(file.path(LOG, "cache", "frozen_candidates.rds"))
  sett <- readRDS(file.path(LOG, "cache", "rcd_sets.rds"))
  de_disc <- readRDS(file.path(LOG, "cache", "GSE70947_DE.rds"))$de_all
  w <- readRDS(file.path(LOG, "cache", "GSE70947_WGCNA.rds"))
  metric <- read.csv(file.path(TAB, "GSE70947_WGCNA_gene_level_metrics.csv"), stringsAsFactors = FALSE)
  sa <- read.csv(file.path(TAB, "GSE70947_RCD_associated_modules.csv"), stringsAsFactors = FALSE)
  vres <- readRDS(file.path(LOG, "cache", "validation_results.rds"))
  cvm <- vres$cvm
  meta <- tryCatch(readRDS(file.path(LOG, "cache", "meta_results.rds")), error = function(e) NULL)
  pres <- tryCatch(readRDS(file.path(LOG, "cache", "preservation_results.rds")), error = function(e) NULL)
  clin <- tryCatch(read.csv(file.path(TAB, "GSE42568_clinical_validation.csv"), stringsAsFactors = FALSE), error = function(e) NULL)

  valid_gses <- c("GSE42568", "GSE15852", "GSE29431")
  genes <- frozen$Gene
  n <- length(genes)

  fev <- data.frame(Gene = genes, stringsAsFactors = FALSE)
  fev$Ferroptosis_membership <- as.integer(fev$Gene %in% sett$Ferroptosis)
  fev$Necroptosis_membership <- as.integer(fev$Gene %in% sett$Necroptosis)
  fev$Cuproptosis_membership <- as.integer(fev$Gene %in% sett$Cuproptosis)
  fev$Number_RCD_pathways <- fev$Ferroptosis_membership + fev$Necroptosis_membership + fev$Cuproptosis_membership

  m <- match(fev$Gene, frozen$Gene)
  fev$Discovery_logFC <- frozen$logFC[m]
  fev$Discovery_FDR   <- frozen$DE_FDR[m]
  fev$WGCNA_module    <- frozen$Module[m]
  fev$ModuleMembership <- frozen$ModuleMembership[m]
  fev$GeneSignificance <- frozen$GeneSignificance[m]
  fev$Connectivity     <- frozen$IntramodularConnectivity[m]
  fev$RCD_module_trait_evidence <- ifelse(frozen$Module_trait_FDR[m] < 0.05, "yes", "no")
  fev$RCD_module_enrichment_evidence <- ifelse(frozen$Module_RCD_enrichment_FDR[m] < 0.05, "yes", "no")

  ## validation logFC/FDR
  for (g in valid_gses) {
    fev[[paste0(g, "_logFC")]] <- cvm[[paste0(g, "_logFC")]][match(fev$Gene, cvm$Gene)]
    fev[[paste0(g, "_FDR")]]   <- cvm[[paste0(g, "_FDR")]][match(fev$Gene, cvm$Gene)]
  }
  fev$Number_concordant_validation_cohorts <- cvm$N_same_direction[match(fev$Gene, cvm$Gene)]
  fev$Number_significant_validation_cohorts <- cvm$N_FDR_significant[match(fev$Gene, cvm$Gene)]
  fev$Overall_validation_status <- cvm$Overall_validation_status[match(fev$Gene, cvm$Gene)]

  ## meta-analysis (REML if available else FE)
  if (!is.null(meta)) {
    mm <- meta$meta_all[meta$meta_all$Model == "REML", ]
    if (!nrow(mm)) mm <- meta$meta_all[meta$meta_all$Model == "FE", ]
    rownames(mm) <- NULL
    if (nrow(mm)) {
      fev$Meta_effect <- mm$Pooled_effect[match(fev$Gene, mm$Gene)]
      fev$Meta_FDR    <- mm$Meta_FDR[match(fev$Gene, mm$Gene)]
      fev$I2          <- mm$I2[match(fev$Gene, mm$Gene)]
    }
  }
  if (!"Meta_effect" %in% colnames(fev)) fev$Meta_effect <- fev$Meta_FDR <- fev$I2 <- NA_real_

  ## module preservation evidence
  if (!is.null(pres) && length(pres)) {
    ptab <- do.call(rbind, pres)
    if (!is.null(ptab) && nrow(ptab)) {
      p_genes <- sapply(fev$Gene, function(ge) {
        mod <- fev$WGCNA_module[fev$Gene == ge]
        sub <- ptab[ptab$Module == mod, ]
        if (!nrow(sub)) return(NA_character_)
        paste(paste0(sub$Dataset, "=", round(sub$Zsummary, 1)),
              collapse = ";")
      })
      fev$Module_preservation_evidence <- p_genes
    }
  }
  if (!"Module_preservation_evidence" %in% colnames(fev)) fev$Module_preservation_evidence <- NA_character_

  ## clinical validation evidence (GSE42568): min FDR among ER/grade/LN/age/size
  if (!is.null(clin) && nrow(clin)) {
    clin_g <- clin[clin$Gene %in% fev$Gene, ]
    if (nrow(clin_g)) {
      best <- aggregate(clin_g$FDR, by = list(Gene = clin_g$Gene), FUN = min, na.rm = TRUE)
      fev$Clinical_validation_minFDR <- best$x[match(fev$Gene, best$Gene)]
    }
  }
  if (!"Clinical_validation_minFDR" %in% colnames(fev)) fev$Clinical_validation_minFDR <- NA_real_

  ## ---------- classification ----------
  classification <- vector("character", n)
  for (i in seq_len(n)) {
    st <- fev$Overall_validation_status[i]
    meta_sig <- !is.na(fev$Meta_FDR[i]) && fev$Meta_FDR[i] < 0.05
    disc_sig <- fev$Discovery_FDR[i] < 0.05
    hub <- !is.na(fev$ModuleMembership[i]) && abs(fev$ModuleMembership[i]) >= 0.8
    class_i <- "DISCOVERY_ONLY"
    if (st == "ROBUST" && meta_sig && hub) class_i <- "CORE_ROBUST_RCD_HUB"
    else if (st == "ROBUST" || (st %in% c("ROBUST", "SUPPORTED") && meta_sig)) class_i <- "VALIDATED_RCD_HUB"
    else if (st == "SUPPORTED") class_i <- "PARTIALLY_VALIDATED"
    else if (st == "INCONSISTENT") class_i <- "INCONSISTENT"
    else if (st %in% c("NOT_MEASURED", "NO_EVIDENCE")) class_i <- "DISCOVERY_ONLY"
    classification[i] <- class_i
  }
  fev$Final_classification <- classification

  write_tab(fev, "FINAL_RCD_candidate_evidence_table.csv")
  strong <- fev[fev$Final_classification %in% c("CORE_ROBUST_RCD_HUB", "VALIDATED_RCD_HUB"), ]
  write_tab(strong, "FINAL_core_robust_RCD_hub_genes.csv")
  keynum("final_core_robust_hub_genes", nrow(strong))
  keynum("final_validated_hub_genes", sum(fev$Final_classification == "VALIDATED_RCD_HUB"))
  hlog("final classification counts:", paste(names(table(classification)), table(classification), collapse = "; "))
  saveRDS(fev, file.path(LOG, "cache", "final_candidates.rds"))

  hlog("STEP 24 complete; core robust genes:", nrow(strong))
  if (nrow(strong)) print(strong[, c("Gene", "Final_classification", "Meta_FDR", "Number_concordant_validation_cohorts")])
})