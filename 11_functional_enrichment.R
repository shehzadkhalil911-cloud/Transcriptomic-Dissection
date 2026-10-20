## 11_functional_enrichment.R - STEP 23: GO/KEGG ORA & GSEA on full ranked list
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(enrichplot)
})

run_wrapped("11", {
  hlog("=== 11 FUNCTIONAL ENRICHMENT ===")
  de_all <- readRDS(file.path(LOG, "cache", "GSE70947_DE.rds"))$de_all
  frozen <- readRDS(file.path(LOG, "cache", "frozen_candidates.rds"))
  w <- readRDS(file.path(LOG, "cache", "GSE70947_WGCNA.rds"))

  ## ---------- ranked list for GSEA (FULL transcriptome, no FDR filter) ----------
  ranked <- de_all$t
  names(ranked) <- de_all$Gene
  ranked <- sort(ranked, decreasing = TRUE)
  ## map to Entrez for clusterProfiler-based GSEA
  suppressMessages({
    eg <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(ranked),
                                column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  })
  keep <- !is.na(eg) & !duplicated(eg)
  ranked_eg <- ranked[keep]
  names(ranked_eg) <- eg[keep]
  rank_desc <- sort(ranked_eg, decreasing = TRUE)
  hlog("GSEA ranked list:", length(rank_desc), "genes (entrez)")

  gsea_go <- tryCatch(
    gseGO(geneList = rank_desc, OrgDb = org.Hs.eg.db, ont = "BP",
          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
          verbose = FALSE, seed = TRUE),
    error = function(e) NULL)
  if (is.null(gsea_go)) hlog("gseGO failed / no results") else {
    go_df <- as.data.frame(gsea_go)
    if (nrow(go_df)) {
      go_df <- go_df[order(go_df$p.adjust), ]
      go_df$Gene <- NULL
      write_tab(go_df, "GSEA_GSE70947_results.csv")
      hlog("GSEA GO significant terms:", nrow(go_df))
      ## selected pathways figure: top 12 up- and down-regulated enriched terms
      go_top <- go_df[order(go_df$p.adjust), ][1:min(12, nrow(go_df)), ]
      go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))
      gp <- ggplot(go_top, aes(NES, Description)) +
        geom_col(aes(fill = p.adjust), width = 0.7) +
        scale_fill_gradient(low = "red", high = "blue") +
        labs(title = "GSEA: top GO BP terms (GSE70947 tumor vs normal full ranked)",
             x = "Normalized Enrichment Score", y = "") +
        theme_bw(base_size = 10)
      save_gg(gp, "GSEA_selected_pathways.png", w = 9, h = 7)
    }
  }
  keynum("GSEA_ranked_genes", length(rank_desc))

  ## ---------- ORA on final robust candidates ----------
  final <- tryCatch(readRDS(file.path(LOG, "cache", "final_candidates.rds")),
                    error = function(e) {
                      hlog("final_candidates.rds not found; falling back to frozen list")
                      frozen
                    })
  cand_genes <- unique(final$Gene)
  hlog("ORA candidate gene count:", length(cand_genes))

  enrich_go <- enrich_kegg <- NULL
  if (length(cand_genes) >= 5) {
    suppressMessages({
      eg <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = cand_genes,
                                  column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    })
    eg <- unique(unname(eg[!is.na(eg)]))
    hlog("candidates mapped to Entrez:", length(eg))
    if (length(eg) >= 5) {
      enrich_go <- tryCatch(enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP",
                                     pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.2, readable = TRUE),
                            error = function(e) NULL)
      enrich_kegg <- tryCatch(enrichKEGG(gene = eg, organism = "hsa",
                                         pvalueCutoff = 0.05),
                              error = function(e) NULL)
      if (!is.null(enrich_go) && nrow(as.data.frame(enrich_go))) {
        go_ora <- as.data.frame(enrich_go)
        write_tab(go_ora[order(go_ora$p.adjust), ], "GO_enrichment_final_RCD_candidates.csv")
        p1 <- dotplot(enrich_go, showCategory = 12) +
          labs(title = "GO BP enrichment: final robust RCD candidates")
        save_gg(p1, "GO_enrichment_final_candidates.png", w = 9, h = 7)
        hlog("GO ORA terms:", nrow(go_ora))
      }
      if (!is.null(enrich_kegg) && nrow(as.data.frame(enrich_kegg))) {
        keg_df <- as.data.frame(enrich_kegg)
        write_tab(keg_df[order(keg_df$p.adjust), ], "KEGG_enrichment_final_RCD_candidates.csv")
        p2 <- dotplot(enrich_kegg, showCategory = 12) +
          labs(title = "KEGG enrichment: final robust RCD candidates")
        save_gg(p2, "KEGG_enrichment_final_candidates.png", w = 9, h = 6)
        hlog("KEGG ORA terms:", nrow(keg_df))
      }
    }
  } else {
    hlog("ORA skipped: candidate gene count < 5 (module-level and GSEA emphasized instead)")
  }

  saveRDS(list(gsea_go = if (is.null(gsea_go)) NULL else as.data.frame(gsea_go),
               enrich_go = if (is.null(enrich_go)) NULL else as.data.frame(enrich_go),
               enrich_kegg = if (is.null(enrich_kegg)) NULL else as.data.frame(enrich_kegg)),
          file.path(LOG, "cache", "enrichment_results.rds"))
  hlog("STEP 23 complete")
})