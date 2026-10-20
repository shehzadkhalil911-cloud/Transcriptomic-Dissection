## 00_setup_and_inventory.R - STEP 1: data inventory & structure inspection
## Creates dirs, tables/dataset_inventory.csv, logs/dataset_structure_summary.txt
source(file.path("scripts", "common.R"))

run_wrapped("00", {
  hlog("=== 00 SETUP & INVENTORY ===")
  write_session_info()
  dir.create(TAB, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)
  dir.create(SCR, showWarnings = FALSE); dir.create(LOG, showWarnings = FALSE)

  ## ---- parameters record (STEP 27 partial, completed by 13) ----
  params <- c(
    "random_seed=12345",
    "stringAsFactors=FALSE",
    "DE_FDR_threshold=0.05",
    "DE_logFC_primary_threshold=0.5",
    "DE_logFC_sensitivity_threshold=1",
    "WGCNA_networkType=signed",
    "WGCNA_TOMType=signed",
    "WGCNA_gene_filter_method=top75_percent_by_MAD_overall",
    "WGCNA_minModuleSize=30",
    "WGCNA_deepSplit=2",
    "WGCNA_mergeCutHeight=0.25",
    "WGCNA_corType=bicor",
    "WGCNA_maxBlockSize=8000",
    "hub_strong_MM_threshold=0.8",
    "hub_corFDR_threshold=0.05",
    "pathway_score_method=GSVA_ssgsea",
    "module_RCD_enrichment_test=Fisher_exact",
    "module_RCD_enrichment_FDR=0.05",
    "validation_robust=same_direction_in_all_3_and_FDR_significant_in_>=2",
    "validation_supported=same_direction_in_>=2_with_nominal_psupport",
    "meta_analysis=metafor_SMD_fixed_and_random",
    "module_preservation_permutations=200"
  )
  writeLines(params, file.path(LOG, "analysis_parameters.txt"))
  hlog("parameters skeleton written")

  ## ---- inventory ----
  gses <- c("GSE70947", "GSE42568", "GSE15852", "GSE29431")
  inv <- data.frame(
    GEO_accession = gses,
    filename_RDS  = paste0(gses, "_processed.rds"),
    platform      = NA_character_,
    object_class  = NA_character_,
    n_genes       = NA_integer_,
    n_samples     = NA_integer_,
    n_tumor       = NA_integer_,
    n_normal      = NA_integer_,
    paired        = NA_character_,
    clinical_vars = NA_character_,
    missing_info  = NA_character_,
    stringsAsFactors = FALSE
  )

  structure_con <- file(file.path(LOG, "dataset_structure_summary.txt"), open = "wt")
for (i in seq_along(gses)) {
    g <- gses[i]
    hlog("inspecting", g)
    x <- load_gse(g)
    cls <- paste(class(x), collapse = ",")
    inv$object_class[i] <- cls
    inv$platform[i] <- if (length(x$platform)) paste(x$platform, collapse = ",") else NA_character_
    emat <- x$expression
    pd   <- x$phenotype

    # ---- label inference (documented rules) ----
    rule <- NA_character_
    if (g == "GSE70947") {
      tissue <- trimws(sub("tissue: *", "", as.character(pd[["tissue:ch1"]])))
      if (all(tissue %in% c("normal", "tumor"))) { cond <- tissue; rule <- "tissue:ch1 [normal/tumor] + sample id CM016.N/T"
        sid <- as.character(pd[["sample id:ch1"]])
        patient <- sub("\\.c[0-9].*$", "", sid); patient <- sub("\\..*$", "", patient)
        inv$paired[i] <- "yes (148 tumor-normal matches by sample id)"
      } else cond <- NA
    } else if (g == "GSE42568") {
      src <- as.character(pd$source_name_ch1)
      if (all(src %in% c("Breast tissue, cancer", "Breast tissue, normal"))) {
        cond <- ifelse(grepl("cancer", src), "tumor", "normal"); rule <- "source_name_ch1 ['Breast tissue, cancer/normal']"
        inv$paired[i] <- "no (104T + 17N unpaired)"
      } else cond <- NA
    } else if (g == "GSE15852") {
      src <- as.character(pd$source_name_ch1)
      if (all(src %in% c("breast tumor tissue", "normal breast tissue"))) {
        cond <- ifelse(grepl("tumor", src), "tumor", "normal")
        rule <- "source_name_ch1 + title 'Normal/Cancer BCxxxxN/T' (43 pairs)"
        inv$paired[i] <- "yes (43 tumor-normal matches)"
      } else cond <- NA
    } else if (g == "GSE29431") {
      ds <- as.character(pd[["disease state:ch1"]])
      if (all(grepl("^Tumor|^none$", ds))) { cond <- ifelse(grepl("^none$", ds), "normal", "tumor")
        rule <- "disease state:ch1 [Tumor... / none]; 12 normal + 54 tumor"
        inv$paired[i] <- "no (54T + 12N unpaired)"
      } else cond <- NA
    } else cond <- NA

    nT <- sum(cond == "tumor", na.rm = TRUE); nN <- sum(cond == "normal", na.rm = TRUE)
    inv$n_tumor[i] <- nT; inv$n_normal[i] <- nN

    clin <- grep(":ch1", colnames(pd), value = TRUE)
    clin <- setdiff(clin, c("tissue:ch1", "tumor_size_cm:ch1", "age_in_years:ch1",
                            "bmi:ch1", "er:ch1", "her2:ch1", "cohort:ch1"))
## ---- text structure summary ----
    writeLines(sprintf("\n================ %s ================", g), structure_con)
    writeLines(sprintf("RDS class: %s", cls), structure_con)
    writeLines(sprintf("platform: %s", inv$platform[i]), structure_con)
    writeLines(sprintf("expression (gene-level) dim: %d x %d", nrow(emat), ncol(emat)), structure_con)
    writeLines(sprintf("expression_raw dim: %d x %d", nrow(x$expression_raw), ncol(x$expression_raw)), structure_con)
    writeLines(sprintf("expression_norm dim: %d x %d", nrow(x$expression_norm), ncol(x$expression_norm)), structure_con)
    writeLines(sprintf("expression rownames[1:8]: %s", paste(head(rownames(emat), 8), collapse = " | ")), structure_con)
    writeLines(sprintf("expression colnames[1:5]: %s", paste(head(colnames(emat), 5), collapse = " | ")), structure_con)
    writeLines(sprintf("expression NA:%d (%.3f%%)  | duplicate rownames:%d", sum(is.na(emat)),
                       100 * mean(is.na(emat)), sum(duplicated(rownames(emat)))), structure_con)
    writeLines(sprintf("expression value range: %.2f .. %.2f", min(emat, na.rm = TRUE), max(emat, na.rm = TRUE)), structure_con)
    writeLines(sprintf("scale_info: %s", paste(unlist(lapply(x$scale_info, function(z) paste(head(z, 1)))), collapse = " | ")), structure_con)
    writeLines(sprintf("phenotype dim: %d x %d", nrow(pd), ncol(pd)), structure_con)
    writeLines(sprintf("phenotype columns: %s", paste(colnames(pd), collapse = "; ")), structure_con)
    writeLines(sprintf("pData title[1:4]: %s", paste(head(pd$title, 4), collapse = " | ")), structure_con)
    writeLines(sprintf("pData title[last:4]: %s", paste(tail(pd$title, 4), collapse = " | ")), structure_con)
    writeLines(sprintf("label inference rule: %s", rule), structure_con)
    writeLines(sprintf("inferred condition counts: tumor=%d normal=%d", nT, nN), structure_con)
    if (g == "GSE70947") {
      sid <- sub("\\..*$", "", as.character(pd[["sample id:ch1"]]))
      writeLines(sprintf("GSE70947 pairing: %d distinct patient ids; paired N/T = %d",
                         length(unique(sid)), sum(table(sid) == 2)), structure_con)
      writeLines(sprintf("GSE70947 er pos:%d neg:%d | her2 pos:%d neg:%d NA:%d | age missing:%d | bmi missing:%d | size missing:%d",
                         sum(pd[["er:ch1"]] == "pos"), sum(pd[["er:ch1"]] == "neg"),
                         sum(pd[["her2:ch1"]] == "pos"), sum(pd[["her2:ch1"]] == "neg"), sum(is.na(pd[["her2:ch1"]])),
                         sum(is.na(pd[["age_in_years:ch1"]])), sum(is.na(pd[["bmi:ch1"]])), sum(is.na(pd[["tumor_size_cm:ch1"]]))), structure_con)
    }
    if (g == "GSE42568") {
      writeLines(sprintf("GSE42568 er_status 1:%d 0:%d NA:%d | grade 1:%d 2:%d 3:%d",
                         sum(pd[["er_status:ch1"]] == 1, na.rm = TRUE), sum(pd[["er_status:ch1"]] == 0, na.rm = TRUE), sum(is.na(pd[["er_status:ch1"]])),
                         sum(pd[["grade:ch1"]] == 1, na.rm = TRUE), sum(pd[["grade:ch1"]] == 2, na.rm = TRUE), sum(pd[["grade:ch1"]] == 3, na.rm = TRUE)), structure_con)
      writeLines(sprintf("GSE42568 OS events:%d of %d non-NA | RFS events:%d | age missing:%d",
                         sum(pd[["overall survival event:ch1"]] == 1, na.rm = TRUE), sum(!is.na(pd[["overall survival event:ch1"]])),
                         sum(pd[["relapse free survival event:ch1"]] == 1, na.rm = TRUE),
                         sum(is.na(pd[["age:ch1"]]))), structure_con)
    }
    writeLines(paste(rep("=", 60), collapse = ""), structure_con)
    inv$n_genes[i] <- nrow(emat); inv$n_samples[i] <- ncol(emat)
  }
  close(structure_con)
  write_tab(inv, "dataset_inventory.csv")
  hlog("inventory saved")
})
    inv$clinical_vars[i] <- paste(clin, collapse = ";")
    inv$missing_info[i] <- ifelse(sum(is.na(pd)) > 0,
                                  paste("NA cells in phenotype:", sum(is.na(pd))), "none observed")
