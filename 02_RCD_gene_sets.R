## 02_RCD_gene_sets.R - STEP 3: prepare RCD gene sets & overlap membership
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({library(readxl); library(ComplexUpset); library(patchwork)})

## Robust single-column gene list reader (handles idiosyncratic headers)
read_gene_list <- function(f) {
  if (grepl("\\.xlsx$", f, ignore.case = TRUE)) {
    d <- read_excel(f, sheet = 1, col_names = FALSE)
  } else {
    d <- read.csv(f, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  g <- as.character(d[[1]])
  g <- trimws(g)
  lab <- grepl("^(symbol|gene|genes|genesymbol|gene_name|hgnc_symbol|name|names)$", g, ignore.case = TRUE)
  g <- g[!lab]
  g <- sub("^([A-Za-z0-9._-]+).*$", "\\1", g)   ## keep first token if messy entries
  g <- g[!is.na(g) & nzchar(g)]
  g <- unique(g)
  g
}

run_wrapped("02", {
  hlog("=== 02 RCD GENE SETS ===")
  ferr  <- read_gene_list(file.path(DD, "Ferroptosis Genes.xlsx"))
  necr  <- read_gene_list(file.path(DD, "Necroptosis genes.xlsx"))
  cupr  <- read_gene_list(file.path(DD, "Cuproptosis genes.csv"))

  hlog("raw counts: ferroptosis", length(ferr), "| necroptosis", length(necr), "| cuproptosis", length(cupr))
  keynum("RCD_raw_ferroptosis", length(ferr)); keynum("RCD_raw_necroptosis", length(necr)); keynum("RCD_raw_cuproptosis", length(cupr))

  Ferroptosis_gene_set <- ferr; Necroptosis_gene_set <- necr; Cuproptosis_gene_set <- cupr
  RCD_union <- unique(c(ferr, necr, cupr))

  overlap3 <- Reduce(intersect, list(ferr, necr, cupr))
  pair_fer_nec <- intersect(ferr, necr)
  pair_fer_cup <- intersect(ferr, cupr)
  pair_nec_cup <- intersect(necr, cupr)

  ## one row per unique gene with binary membership
  genes <- sort(RCD_union)
  membership <- data.frame(
    Gene = genes,
    Ferroptosis  = as.integer(genes %in% ferr),
    Necroptosis  = as.integer(genes %in% necr),
    Cuproptosis  = as.integer(genes %in% cupr),
    stringsAsFactors = FALSE)
  membership$Number_of_RCD_pathways <- rowSums(membership[, 2:4, drop = FALSE])

  ## clean CSV outputs
  write_tab(data.frame(Gene = ferr), "ferroptosis_gene_list_clean.csv")
  write_tab(data.frame(Gene = necr), "necroptosis_gene_list_clean.csv")
  write_tab(data.frame(Gene = cupr), "cuproptosis_gene_list_clean.csv")
  write_tab(data.frame(Gene = RCD_union, Pathway_overlap = ifelse(RCD_union %in% overlap3, "Ferroptosis+Necroptosis+Cuproptosis",
                        ifelse(RCD_union %in% pair_fer_nec, "Ferroptosis+Necroptosis",
                        ifelse(RCD_union %in% pair_fer_cup, "Ferroptosis+Cuproptosis",
                        ifelse(RCD_union %in% pair_nec_cup, "Necroptosis+Cuproptosis", "single"))))), "RCD_union_gene_list.csv")
  write_tab(membership, "RCD_gene_overlap_membership.csv")

  ## UpSet plot
  sdata <- data.frame(row.names = genes,
                      Ferroptosis = as.integer(genes %in% ferr),
                      Necroptosis = as.integer(genes %in% necr),
                      Cuproptosis = as.integer(genes %in% cupr))
  ## NOTE: UpSetR::upset silently renders blank PNGs in this R 4.6/ggplot2 environment
  ## (verified in isolation) -> use ComplexUpset (ggplot2-based) instead.
  p_up <- ComplexUpset::upset(sdata, intersect = c("Ferroptosis", "Necroptosis", "Cuproptosis"),
                              base_annotations = list(
                                "Intersection size" = ComplexUpset::intersection_size(counts = TRUE)),
                              width_ratio = 0.15) +
    patchwork::plot_annotation(
      title = "Overlap of ferroptosis, necroptosis and cuproptosis gene sets",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14)))
  save_gg(p_up, "RCD_gene_set_overlap_upset.png", w = 8, h = 6)

  ## save sets for downstream use
  settlist <- list(Ferroptosis = ferr, Necroptosis = necr, Cuproptosis = cupr, RCD_union = RCD_union,
                   overlap3 = overlap3, fer_nec = pair_fer_nec, fer_cup = pair_fer_cup, nec_cup = pair_nec_cup)
  saveRDS(settlist, file.path(LOG, "cache", "rcd_sets.rds"))

  hlog("RCD overlaps: FnN=", length(pair_fer_nec), " FnC=", length(pair_fer_cup),
       " NnC=", length(pair_nec_cup), " triple=", length(overlap3))
  keynum("RCD_union_count", length(RCD_union))
  keynum("RCD_overlap_fer_nec", length(pair_fer_nec)); keynum("RCD_overlap_fer_cup", length(pair_fer_cup))
  keynum("RCD_overlap_nec_cup", length(pair_nec_cup)); keynum("RCD_overlap_triple", length(overlap3))
})