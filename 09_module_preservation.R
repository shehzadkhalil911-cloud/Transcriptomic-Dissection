## 09_module_preservation.R - STEP 21: preserve GSE70947 modules in validation tumor cohorts
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({ library(WGCNA); library(ggplot2) })
options(stringsAsFactors = FALSE)

## helper: recursively find a data.frame containing any of want_names as columns
find_named_df <- function(x, want_names, depth = 0) {
  if (depth > 5) return(NULL)
  if (is.data.frame(x)) {
    if (any(want_names %in% colnames(x))) return(x)
    return(NULL)
  }
  if (is.list(x)) for (i in seq_along(x)) {
    r <- find_named_df(x[[i]], want_names, depth + 1)
    if (!is.null(r)) return(r)
  }
  NULL
}

empty_pres_tab <- function(g) {
  data.frame(Dataset = g, Module = character(0), Module_size = integer(0),
             Zsummary = numeric(0), medianRank = numeric(0), Zsummary_p = numeric(0),
             Interpretation = character(0), stringsAsFactors = FALSE)
}

run_wrapped("09", {
  hlog("=== 09 MODULE PRESERVATION ===")
  w <- readRDS(file.path(LOG, "cache", "GSE70947_WGCNA.rds"))
  datExpr <- w$datExpr
  modCols <- w$moduleColors
  genes_all <- colnames(datExpr)

  ## reference gene set: top 100 genes/module (non-grey) by |kME| (permutation tractability)
  mdf <- read.csv(file.path(TAB, "GSE70947_WGCNA_gene_level_metrics.csv"), stringsAsFactors = FALSE)
  sel <- integer(0)
  for (m in setdiff(unique(modCols), "grey")) {
    mm <- mdf[mdf$Module == m, ]
    mm <- mm[order(-abs(mm$kME)), ]
    sel <- c(sel, which(genes_all %in% head(mm$Gene, 100)))
  }
  sel <- unique(sel)
  refGenes <- genes_all[sel]
  hlog("preservation reference gene set:", length(refGenes), "genes (top ~100/module)")

  nPerm <- suppressWarnings(as.integer(Sys.getenv("PRES_NPERM", "100")))
  if (is.na(nPerm) || nPerm < 2) nPerm <- 100
  hlog("nPermutations:", nPerm)

  setNameRef <- "DISCOVERY"
  valid_gses <- c("GSE42568", "GSE15852", "GSE29431")
  pres_rows <- list()

  ## per-cohort execution mode: PRES_GSE=<cohort> runs one cohort and caches it;
  ## PRES_GSE=COMBINE assembles cached results; unset runs all cohorts + combine
  oneGse <- Sys.getenv("PRES_GSE", "")
  do_combine <- identical(oneGse, "COMBINE")
  run_gses <- if (oneGse %in% valid_gses) oneGse else valid_gses
  if (oneGse %in% valid_gses) hlog("single-cohort mode:", oneGse)
  if (do_combine) hlog("combine mode: assembling cached preservation results")
  if (!do_combine) for (g in run_gses) {
    hlog("preservation in", g)
    emat <- readRDS(file.path(LOG, "cache", paste0(g, "_std_expr.rds")))
    meta <- readRDS(file.path(LOG, "cache", paste0(g, "_meta.rds")))
    tumors <- meta$sample[meta$condition == "tumor"]
    genes_use <- intersect(refGenes, rownames(emat))
    if (length(genes_use) < 100) { hlog(g, "insufficient overlapping genes; skipping preservation"); next }
    nTumors <- length(tumors)
    hlog(g, "tumor samples for preservation:", nTumors, "; candidate genes:", length(genes_use))
    if (nTumors < 10) {
      hlog(g, "tumor cohort too small for stable preservation; recording and skipping")
      write_tab(empty_pres_tab(g), paste0(g, "_module_preservation.csv"))
      next
    }

    testExpr <- t(emat[genes_use, tumors, drop = FALSE])
    testExpr <- testExpr[, match(genes_use, colnames(testExpr)), drop = FALSE]
    gsgT <- goodSamplesGenes(testExpr)
    if (!all(gsgT$goodGenes)) {
      hlog(g, "removing", sum(!gsgT$goodGenes), "zero-variance/missing genes from test set")
      testExpr <- testExpr[, gsgT$goodGenes, drop = FALSE]
      genes_use <- genes_use[gsgT$goodGenes]
    }
    if (ncol(testExpr) < 50) { hlog(g, "too few genes after QC (", ncol(testExpr), "); skipping"); next }

    refExpr <- datExpr[, genes_use, drop = FALSE]
    gsgR <- goodSamplesGenes(refExpr)
    if (!all(gsgR$goodGenes)) {
      hlog(g, "removing", sum(!gsgR$goodGenes), "zero-variance/missing genes from reference set")
      refExpr <- refExpr[, gsgR$goodGenes, drop = FALSE]
      testExpr <- testExpr[, gsgR$goodGenes, drop = FALSE]
      genes_use <- genes_use[gsgR$goodGenes]
    }
    if (ncol(refExpr) < 50) { hlog(g, "too few genes after QC (", ncol(refExpr), "); skipping"); next }

    refCols_use <- modCols[match(genes_use, genes_all)]
    multiData <- list()
    multiData[[setNameRef]] <- list(data = as.data.frame(refExpr))
    multiData[[g]]          <- list(data = as.data.frame(testExpr))
    multiColor <- list()
    multiColor[[setNameRef]] <- refCols_use
    multiColor[[g]]          <- refCols_use

    mp <- tryCatch(
      modulePreservation(multiData, multiColor,
                         referenceNetworks = 1, nPermutations = nPerm,
                         networkType = "signed", corFnc = "bicor",
                         randomSeed = SEED, verbose = 0, indent = 0),
      error = function(e) { hlog(g, "modulePreservation error:", conditionMessage(e)); NULL })
    if (is.null(mp)) next
    refTag <- paste0("ref.", setNameRef)
    Zs <- tryCatch(mp$preservation$Z[[refTag]]$Zsummary.pres, error = function(e) NULL)
    if (is.null(Zs)) Zs <- find_named_df(mp, c("Zsummary.pres", "Zsummary"))
    MRs <- tryCatch(mp$preservation$medianRank[[refTag]]$medianRank.pres, error = function(e) NULL)
    if (is.null(MRs)) MRs <- find_named_df(mp, c("medianRank.pres", "medianRank"))
    Zp <- tryCatch(mp$preservation$Z[[refTag]]$Zsummary.p, error = function(e) NULL)
    if (is.null(Zp)) Zp <- find_named_df(mp, c("Zsummary.p", "Zsummary.pval"))

    if (is.null(Zs) || !nrow(Zs)) {
      hlog(g, "no Zsummary found; top-level preservation names:",
           paste(names(mp$preservation), collapse = ","))
      write_tab(empty_pres_tab(g), paste0(g, "_module_preservation.csv"))
      next
    }

    zs <- data.frame(Module = rownames(Zs), stringsAsFactors = FALSE)
    zcol <- intersect(c("Zsummary.pres", "Zsummary"), colnames(Zs))[1]
    zs$Zsummary <- Zs[[zcol]]
    if ("module.size" %in% colnames(Zs)) zs$Module_size <- Zs[["module.size"]]
    if (!is.null(MRs) && nrow(MRs)) {
      mrcol <- intersect(c("medianRank.pres", "medianRank"), colnames(MRs))[1]
      if (!is.na(mrcol)) zs$medianRank <- MRs[match(zs$Module, rownames(MRs)), mrcol]
      else zs$medianRank <- NA_real_
    } else zs$medianRank <- NA_real_
    if (!is.null(Zp) && nrow(Zp)) {
      zpcol <- intersect(c("Zsummary.p", "Zsummary.pval"), colnames(Zp))[1]
      if (!is.na(zpcol)) zs$Zsummary_p <- Zp[match(zs$Module, rownames(Zp)), zpcol]
      else zs$Zsummary_p <- NA_real_
    } else zs$Zsummary_p <- NA_real_
    zs$Dataset <- g
    zs <- zs[, c("Dataset", "Module", setdiff(names(zs), c("Dataset", "Module"))), drop = FALSE]
    zs <- zs[zs$Module != "grey", , drop = FALSE]
    zs$Interpretation <- ifelse(is.na(zs$Zsummary), NA_character_,
                          ifelse(zs$Zsummary > 10, "strong",
                          ifelse(zs$Zsummary > 2, "moderate", "weak/none")))
    write_tab(zs, paste0(g, "_module_preservation.csv"))
    saveRDS(list(tab = zs), file.path(LOG, "cache", paste0("pres_", g, ".rds")))
    pres_rows[[g]] <- zs
    hlog(g, "preservation summary written; modules:", nrow(zs),
         "; Zsummary range:", paste(sprintf("%.2f", range(zs$Zsummary, na.rm = TRUE)), collapse = "-"))
  }
  ## combined summary figure (assembled from per-cohort caches)
  cache_files <- file.path(LOG, "cache", paste0("pres_", valid_gses, ".rds"))
  cached <- lapply(cache_files, function(f) tryCatch(readRDS(f)$tab, error = function(e) NULL))
  names(cached) <- valid_gses
  ## fallback: read per-cohort CSV tables when rds caches are missing (single-process runs)
  for (gse in valid_gses) {
    if (is.null(cached[[gse]]) || !nrow(cached[[gse]])) {
      tt <- tryCatch(read.csv(file.path(TAB, paste0(gse, "_module_preservation.csv")),
                              stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(tt) && nrow(tt)) cached[[gse]] <- tt
    }
  }
  cached <- cached[sapply(cached, function(x) !is.null(x) && nrow(x) > 0)]
  if (length(cached)) allp <- do.call(rbind, cached) else allp <- NULL
  if (oneGse %in% valid_gses) {
    hlog("single-cohort result cached; figure deferred to COMBINE pass")
    allp <- NULL
  }
  if (!is.null(allp) && nrow(allp)) {
    allp$Module <- factor(allp$Module, levels = sort(unique(allp$Module)))
    p <- ggplot(allp, aes(Module, Zsummary, fill = Dataset)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 2,  lty = 2, color = "orange") +
      geom_hline(yintercept = 10, lty = 2, color = "red") +
      scale_fill_manual(values = c(GSE42568 = "#66c2a5", GSE15852 = "#fc8d62", GSE29431 = "#8da0cb")) +
      labs(title = "Preservation of GSE70947 WGCNA modules in validation tumor cohorts",
           subtitle = "Dashed lines: Zsummary = 2 (moderate) and 10 (strong)",
           y = "Zsummary", x = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_gg(p, "WGCNA_module_preservation_summary.png", w = 10, h = 6)
  } else {
    hlog("no preservation results to plot")
  }
  saveRDS(cached, file.path(LOG, "cache", "preservation_results.rds"))
  hlog("STEP 21 complete")
})