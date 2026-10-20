## 06_RCD_WGCNA_integration.R - STEPS 11-15: module-RCD enrichment,
## RCD-associated modules, gene-level metrics, hub candidates, freeze list
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2); library(matrixStats); library(ggrepel)
})
options(stringsAsFactors = FALSE)

run_wrapped("06", {
  hlog("=== 06 RCD-WGCNA INTEGRATION ===")
  w   <- readRDS(file.path(LOG, "cache", "GSE70947_WGCNA.rds"))
  de  <- readRDS(file.path(LOG, "cache", "GSE70947_DE.rds"))
  rc  <- readRDS(file.path(LOG, "cache", "GSE70947_RCD_DE_results.rds"))
  sett <- readRDS(file.path(LOG, "cache", "rcd_sets.rds"))
  sc  <- readRDS(file.path(LOG, "cache", "GSE70947_scores.rds"))

  datExpr <- w$datExpr; mods <- w$moduleColors; MEs <- w$MEs; power <- w$power
  de_all <- de$de_all
  genes_net <- colnames(datExpr)
  tumors  <- w$tumors_final
  all_mods <- setdiff(unique(mods), "grey")

  ## ================= STEP 11: module enrichment for RCD gene sets =================
  pathways <- list(Ferroptosis = sett$Ferroptosis, Necroptosis = sett$Necroptosis,
                   Cuproptosis = sett$Cuproptosis, RCD_union = sett$RCD_union)
  enrich_rows <- list()
  for (m in c(all_mods, "grey")) {
    m_genes  <- genes_net[mods == m]
    nm  <- length(m_genes)
    for (pm in names(pathways)) {
      pg <- intersect(pathways[[pm]], genes_net)
      np <- length(pg)
      ov <- intersect(m_genes, pg)
      no <- length(ov)
      a <- no; b <- nm - no; c <- np - no; d <- length(genes_net) - np - nm + no
      ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
                     error = function(e) NULL)
      enrich_rows[[paste(m, pm)]] <- data.frame(
        Module = m, Pathway = pm, Module_gene_count = nm, RCD_gene_count = np,
        Overlap_count = no, Overlap_genes = paste(ov, collapse = ";"),
        Odds_ratio = if (is.null(ft)) NA else unname(ft$estimate),
        P_value = if (is.null(ft)) NA else ft$p.value,
        stringsAsFactors = FALSE)
    }
  }
  enrich <- do.call(rbind, enrich_rows)
  enrich$FDR <- p.adjust(enrich$P_value, "BH")
  hlog("significant module enrichments (FDR<0.05):", sum(enrich$FDR < 0.05, na.rm = TRUE))
  enrich_sig <- enrich[enrich$FDR < 0.05 & !enrich$Module %in% "grey", ]
  if (nrow(enrich_sig)) { print(enrich_sig[order(enrich_sig$FDR), c("Module", "Pathway", "Overlap_count", "P_value", "FDR")]) }
  write_tab(enrich, "GSE70947_WGCNA_RCD_module_enrichment.csv")

  ## enrichment dot plot
  ed <- enrich[enrich$Module != "grey", ]
  ed$neglogFDR <- -log10(ed$FDR + 1e-320)
  p1 <- ggplot(ed, aes(Pathway, Module, size = Overlap_count, color = neglogFDR)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
    scale_size(range = c(2, 8), name = "Overlap") +
    labs(title = "GSE70947 WGCNA module enrichment for RCD gene sets", x = "", y = "Module") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_gg(p1, "GSE70947_WGCNA_RCD_module_enrichment.png", w = 8, h = 7)

  enrich_red <- enrich[enrich$Pathway %in% c("Ferroptosis", "Necroptosis", "Cuproptosis"), ]
## ============ STEP 12: RCD-associated modules ============
  mtFDR <- w$mtFDR; mtCor <- w$mtCor
  score_cols <- grep("^score_", colnames(mtFDR), value = TRUE)
  ## per module, best (min) FDR across the three RCD score traits
  mod_score_fdr <- apply(mtFDR[, score_cols, drop = FALSE], 1, min, na.rm = TRUE)
  mod_score_cor <- apply(mtCor[, score_cols, drop = FALSE], 1, function(r) r[which.min(abs(r))])
  enrich_fdr <- tapply(enrich_red$FDR, enrich_red$Module, min, na.rm = TRUE)
  enrich_ov   <- tapply(enrich_red$Overlap_count, enrich_red$Module, sum, na.rm = TRUE)

  assoc_trait <- names(mod_score_fdr)[mod_score_fdr < 0.05]
  assoc_enrich <- names(enrich_fdr)[enrich_fdr < 0.05 & enrich_fdr != Inf]
  associated <- unique(c(assoc_trait, assoc_enrich))
  hlog("RCD-associated modules (trait FDR<0.05):", paste(assoc_trait, collapse = ","))
  hlog("RCD-associated modules (enrichment FDR<0.05):", paste(assoc_enrich, collapse = ","))
  hlog("combined RCD-associated modules:", paste(associated, collapse = ","))
  keynum("GSE70947_n_RCD_associated_modules", length(associated))

  sa <- data.frame(
    Module = names(mod_score_fdr),
    Module_size = sapply(names(mod_score_fdr), function(m) sum(mods == m)),
    score_trait_minFDR = unname(mod_score_fdr),
    score_trait_cor = unname(mod_score_cor),
    enrichment_minFDR = unname(enrich_fdr[names(mod_score_fdr)]),
    enrichment_overlap = unname(enrich_ov[names(mod_score_fdr)]),
    RCD_associated = names(mod_score_fdr) %in% associated,
    stringsAsFactors = FALSE)
  sa$Which_evidence <- ifelse(sa$Module %in% assoc_trait & sa$Module %in% assoc_enrich, "both",
                       ifelse(sa$Module %in% assoc_trait, "trait", "enrichment"))
  sa$Pathway_association <- sapply(sa$Module, function(m) {
    sub <- enrich_red[enrich_red$Module == m & enrich_red$FDR < 0.05, ]
    if (nrow(sub) == 0) {
      ## fall back to which score trait has min FDR
      s <- score_cols[which.min(mtFDR[m, score_cols])]
      return(gsub("score_", "", s))
    }
    paste(sub$Pathway[sub$FDR == min(sub$FDR)], collapse = "/")
  })
  sa$Overlapping_RCD_genes <- sapply(sa$Module, function(m) {
    sub <- enrich_red[enrich_red$Module == m, ]
    paste(unique(unlist(strsplit(sub$Overlap_genes[sub$Overlap_count > 0], ";"))), collapse = ";")
  })
  write_tab(sa, "GSE70947_RCD_associated_modules.csv")
  hlog("step 12 done; associated modules:", paste(associated, collapse = ","))

  ## ============ STEP 13: gene-level WGCNA metrics ============
  ## module membership (kME) for every gene: cor(gene expr, its module ME)
  ## (MEs columns are bare colors here; direct cor avoids naming pitfalls)
  ## gene significance to RCD scores
  gs_list <- list()
  for (s in score_cols) {
    tv <- w$traits[, s]
    gs <- apply(datExpr, 2, function(g) {
      ok <- !is.na(tv)
      if (sum(ok) < 10) return(NA_real_)
      suppressWarnings(cor(g[ok], tv[ok]))
    })
    gs_list[[s]] <- gs
  }
  gsmat <- do.call(cbind, gs_list)   # genes x scores
  colnames(gsmat) <- gsub("score_", "GS_", score_cols)

  ## per-gene: MM to own module, GS to relevant score, kWithin
  metric_rows <- list()
  for (m in c(all_mods, "grey")) {
    mg <- genes_net[mods == m]
    if (length(mg) < 2) next
    ME_m <- MEs[, m]                     ## columns of MEs are bare colors
    mm <- apply(datExpr[, mg, drop = FALSE], 2, function(g) cor(g, ME_m, use = "p"))
    kw <- tryCatch(softConnectivity(datExpr[, mg], power = power, type = "signed",
                                    corFnc = "bicor"), error = function(e) rep(NA, length(mg)))
    metric_rows[[m]] <- data.frame(
      Gene = mg, Module = m, kME = mm,
      GS_ferroptosis = gsmat[mg, "GS_Ferroptosis"],
      GS_necroptosis = gsmat[mg, "GS_Necroptosis"],
      GS_cuproptosis = gsmat[mg, "GS_Cuproptosis"],
      kWithin = kw, stringsAsFactors = FALSE)
  }
  metrics <- do.call(rbind, metric_rows)
  rownames(metrics) <- NULL
  write_tab(metrics, "GSE70947_WGCNA_gene_level_metrics.csv")
  hlog("gene-level metrics written for", nrow(metrics), "genes")

  ## MM vs GS plots for top associated modules per pathway
  for (pm in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
    am <- sa$Module[sa$RCD_associated & grepl(pm, sa$Pathway_association)]
    if (!length(am)) am <- sa$Module[sa$RCD_associated][1]
    am <- am[1]
    g <- metrics[metrics$Module == am & !is.na(metrics$kME), ]
    col <- gsub("GS_", "", paste0("GS_", pm))
    if (col %in% colnames(g)) {
      pp <- ggplot(g, aes(kME, .data[[col]])) +
        geom_point(color = "steelblue", alpha = 0.6, size = 1.2) +
        geom_smooth(method = "lm", se = TRUE, color = "darkred") +
        labs(title = paste("GSE70947 module", am, "MM vs GS", pm),
             x = "Module membership (kME)", y = paste("Gene significance (", pm, ")")) +
        theme_bw(base_size = 11)
      save_gg(pp, paste0("GSE70947_MM_vs_GS_", tolower(pm), "_module_", am, ".png"), w = 6.5, h = 5.5)
    }
  }
  hlog("step 13 done")
## ============ STEP 14: discovery RCD hub candidates ============
  ## candidate logic: RCD gene AND significant DEG AND in RCD-associated module
  deg_sig <- de_all[de_all$adj.P.Val < 0.05, ]
  rcd_genes <- sett$RCD_union
  assoc_genes <- metrics[metrics$Module %in% associated, ]
  cand <- merge(assoc_genes, deg_sig[, c("Gene", "logFC", "adj.P.Val")], by = "Gene")
  cand <- cand[cand$Gene %in% rcd_genes, ]
  cand <- merge(cand, enrich_red[, c("Module", "Pathway", "FDR")],
                by = "Module", all.x = TRUE)
  ## module trait FDR for the module
  cand$Module_trait_FDR <- unname(mod_score_fdr[cand$Module])
  cand$DE_direction <- ifelse(cand$logFC > 0, "up", "down")
  cand$RCD_pathways <- sapply(cand$Gene, function(g)
    paste(names(pathways[1:3])[sapply(pathways[1:3], function(s) g %in% s)], collapse = ";"))

  ## integrated hub score: mean of rank-normalized |kME|, |GS|, |kWithin|
  ## GS = the relevant pathway score for this gene (max abs across three)
  cand$GS_maxAbs <- pmax(abs(cand$GS_ferroptosis), abs(cand$GS_necroptosis), abs(cand$GS_cuproptosis), na.rm = TRUE)
  rk <- function(v) rank(-v, na.last = "keep", ties.method = "average")
  cand$Hub_score <- rowMeans(cbind(rk(abs(cand$kME)), rk(cand$GS_maxAbs), rk(cand$kWithin)), na.rm = TRUE)
  ## one row per gene: merge with enrichment table can duplicate genes across pathways
  cand <- cand[order(cand$Hub_score, decreasing = TRUE,
                     cand$Module_trait_FDR), ]
  cand <- cand[!duplicated(cand$Gene), ]
  cand$Hub_rank <- seq_len(nrow(cand))
  cand$Strong_hub <- abs(cand$kME) >= 0.8

  hub_tab <- cand[, c("Gene", "RCD_pathways", "logFC", "adj.P.Val", "DE_direction",
                      "Module", "kME", "GS_maxAbs", "kWithin", "Module_trait_FDR",
                      "Hub_score", "Hub_rank", "Strong_hub")]
  names(hub_tab) <- c("Gene", "RCD_pathway(s)", "logFC", "DE_FDR", "DE_direction",
                      "Module", "ModuleMembership", "GeneSignificance", "IntramodularConnectivity",
                      "Module_trait_FDR", "Hub_score", "Hub_rank", "Strong_hub")
  ## module RCD enrichment FDR (min across pathways for this module)
  hub_tab$Module_RCD_enrichment_FDR <- unname(enrich_fdr[hub_tab$Module])
  hub_tab$Module_RCD_association <- hub_tab$Module %in% associated
  hub_tab <- hub_tab[, c("Gene", "RCD_pathway(s)", "logFC", "DE_FDR", "DE_direction",
                         "Module", "ModuleMembership", "GeneSignificance", "IntramodularConnectivity",
                         "Module_RCD_association", "Module_trait_FDR", "Module_RCD_enrichment_FDR",
                         "Hub_rank", "Hub_score", "Strong_hub")]

  write_tab(hub_tab, "GSE70947_discovery_RCD_hub_candidates.csv")
  keynum("GSE70947_n_hub_candidates", nrow(hub_tab))
  hlog("hub candidates:", nrow(hub_tab))

  ## pathway-specific candidate lists
  for (pm in names(pathways[1:3])) {
    pgenes <- pathways[[pm]]
    sub <- hub_tab[hub_tab$Gene %in% pgenes, ]
    write_tab(sub, paste0("GSE70947_", tolower(pm), "_hub_candidates.csv"))
    keynum(paste0("GSE70947_", tolower(pm), "_hub_candidates"), nrow(sub))
  }

  ## ---------- integrated candidate evidence figure ----------
  plotC <- hub_tab[1:min(25, nrow(hub_tab)), ]
  plotC$Gene <- factor(plotC$Gene, levels = rev(plotC$Gene))
  plotC$absMM <- abs(plotC$ModuleMembership)
  plotC$absGS <- abs(plotC$GeneSignificance)
  plotC$neglogFDR <- -log10(plotC$DE_FDR + 1e-320)
  long <- reshape2::melt(plotC[, c("Gene", "absMM", "absGS", "IntramodularConnectivity",
                                   "neglogFDR", "RCD_pathway(s)", "Module")],
                         id.vars = c("Gene", "RCD_pathway(s)", "Module"),
                         variable.name = "Evidence", value.name = "Value")
  p2 <- ggplot(long, aes(Evidence, Gene, fill = Value)) +
    geom_tile(color = "white") +
    facet_grid(`RCD_pathway(s)` ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradientn(colours = c("#2166ac", "white", "#b2182b")) +
    labs(title = "GSE70947 integrated RCD hub candidate evidence (top 25)",
         x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_gg(p2, "GSE70947_integrated_RCD_candidate_evidence.png", w = 8, h = 8)

  ## ============ STEP 15: freeze discovery candidate list ============
  frozen <- hub_tab
  ## ensure frozen list has a stable Gene + pathway columns for validation joins
  frozen$Source_pipeline <- "GSE70947 discovery"
  write_tab(frozen, "FROZEN_discovery_RCD_candidate_list.csv")
  hlog("FROZEN discovery candidate list:", nrow(frozen), "genes")
  keynum("FROZEN_candidates", nrow(frozen))

  saveRDS(frozen, file.path(LOG, "cache", "frozen_candidates.rds"))
  hlog("STEP 14-15 complete")
})