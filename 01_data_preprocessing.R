## 01_data_preprocessing.R - STEP 2: standardize gene identifiers, build metadata
## Outputs: logs/cache/<GSE>_std_expr.rds, tables/<GSE>_sample_metadata.csv,
##          tables/<GSE>_gene_annotation_summary.csv
source(file.path("scripts", "common.R"))
dir.create(file.path(LOG, "cache"), showWarnings = FALSE)

meta_columns <- c("sample", "GSE", "title", "source_name", "condition", "patient_id",
                  "er", "her2", "age", "tumor_size_cm", "bmi",
                  "grade", "lymph_node_status", "os_event", "os_time_days",
                  "rfs_event", "rfs_time_days", "extra")

make_meta <- function(g, pd, cond, extra = NULL) {
  m <- data.frame(sample = rownames(pd), GSE = g, title = pd$title,
                  source_name = pd$source_name_ch1, condition = cond,
                  stringsAsFactors = FALSE)
  m$patient_id <- NA_character_
  m$er <- m$her2 <- m$age <- m$tumor_size_cm <- m$bmi <- NA_real_
  m$grade <- m$lymph_node_status <- NA_character_
  m$os_event <- m$os_time_days <- m$rfs_event <- m$rfs_time_days <- NA_real_
  m$extra <- ""
  if (!is.null(extra)) for (nm in names(extra)) {
    if (!nm %in% colnames(m)) m[[nm]] <- NA
    m[[nm]] <- extra[[nm]]
  }
  m[, intersect(meta_columns, colnames(m))]
}

run_wrapped("01", {
  hlog("=== 01 DATA PREPROCESSING ===")
  gses <- c("GSE70947", "GSE42568", "GSE15852", "GSE29431")

  for (g in gses) {
    hlog("processing", g)
    x <- load_gse(g)
    emat <- x$expression
    pd   <- x$phenotype

    ## --- condition + patient rules (validated in 00) ---
    if (g == "GSE70947") {
      cond <- trimws(sub("tissue: *", "", as.character(pd[["tissue:ch1"]])))
      sid  <- as.character(pd[["sample id:ch1"]])
      patient <- sub("\\.c[0-9].*$", "", sid); patient <- sub("\\..*$", "", patient)
      extra <- list(er = ifelse(pd[["er:ch1"]] == "pos", 1, ifelse(pd[["er:ch1"]] == "neg", 0, NA)),
                    her2 = ifelse(pd[["her2:ch1"]] == "pos", "pos", ifelse(pd[["her2:ch1"]] == "neg", "neg", NA)),
                    age = suppressWarnings(as.numeric(pd[["age_in_years:ch1"]])),
                    tumor_size_cm = suppressWarnings(as.numeric(pd[["tumor_size_cm:ch1"]])),
                    bmi = suppressWarnings(as.numeric(pd[["bmi:ch1"]])),
                    patient_id = patient)
      extra$cohort <- as.character(pd[["cohort:ch1"]])
    } else if (g == "GSE42568") {
      cond <- ifelse(grepl("cancer", pd$source_name_ch1), "tumor", "normal")
      extra <- list(er = suppressWarnings(as.numeric(pd[["er_status:ch1"]])),
                    age = suppressWarnings(as.numeric(pd[["age:ch1"]])),
                    grade = as.character(pd[["grade:ch1"]]),
                    lymph_node_status = suppressWarnings(as.numeric(pd[["lymph node status:ch1"]])),
                    os_event = suppressWarnings(as.numeric(pd[["overall survival event:ch1"]])),
                    os_time_days = suppressWarnings(as.numeric(pd[["overall survival time_days:ch1"]])),
                    rfs_event = suppressWarnings(as.numeric(pd[["relapse free survival event:ch1"]])),
                    rfs_time_days = suppressWarnings(as.numeric(pd[["relapse free survival time_days:ch1"]])),
                    tumor_size_cm = suppressWarnings(as.numeric(pd[["size:ch1"]])))
    } else if (g == "GSE15852") {
      cond <- ifelse(grepl("tumor", pd$source_name_ch1), "tumor", "normal")
      patient <- sub("^[A-Za-z]+ *", "", pd$title); patient <- sub("[NT]$", "", patient)
      extra <- list(patient_id = patient,
                    age = suppressWarnings(as.numeric(pd[["age in years:ch1"]])),
                    grade = gsub(".*grade ", "", as.character(pd[["grade:ch1"]])),
                    race = as.character(pd[["race:ch1"]]))
    } else if (g == "GSE29431") {
      ds <- as.character(pd[["disease state:ch1"]])
      cond <- ifelse(grepl("^none$", ds), "normal", "tumor")
      ## 'Stage:ch1' actually records grade I/II/III in this cohort
      stage <- as.character(pd[["Stage:ch1"]])
      grade <- gsub("grade ", "", stage)
      grade[grade %in% c("", "NA", "I", "II", "III", "1", "2", "3")] <- NA
      re <- as.character(pd[["re:ch1"]])
      er <- ifelse(grepl("^\\+", re), 1, ifelse(grepl("^\\-", re), 0, NA))
      extra <- list(age = suppressWarnings(as.numeric(pd[["age:ch1"]])),
                    er = er,
                    grade = stage,
                    her2 = ifelse(grepl("IHC3|FISH\\+", pd$title), "pos", "neg"),
                    ki67 = suppressWarnings(as.numeric(pd[["ki67:ch1"]])),
                    p53 = as.character(pd[["p53:ch1"]]),
                    metastasis = as.character(pd[["metastasis:ch1"]]))
    } else stop("unknown gse")

    cond <- ifelse(cond %in% c("tumor", "normal"), cond, NA_character_)
    m <- make_meta(g, pd, cond, if (exists("extra")) extra else NULL)
    rownames(m) <- NULL

    ## --- standardize expression ---
    emat_std <- standardize_matrix(emat)
    gsum <- data.frame(GSE = g,
                       n_expression_rows   = nrow(emat),
                       n_expression_cols   = ncol(emat),
                       n_controls_removed  = sum(grepl(control_pattern, clean_symbols(rownames(emat)))),
                       n_dup_symbols_collapsed = sum(duplicated(clean_symbols(rownames(emat)))),
                       n_genes_final = nrow(emat_std),
                       mapping_note = "gene symbols taken from curated rds 'expression' matrix (probe->gene mapping pre-applied by data generator)",
                       stringsAsFactors = FALSE)

    saveRDS(emat_std, file.path(LOG, "cache", paste0(g, "_std_expr.rds")))
    saveRDS(m,     file.path(LOG, "cache", paste0(g, "_meta.rds")))
    write_tab(gsum, paste0(g, "_gene_annotation_summary.csv"))
    write_tab(m,    paste0(g, "_sample_metadata.csv"))
    keynum(paste0(g, "_n_tumor"), sum(m$condition == "tumor", na.rm = TRUE))
    keynum(paste0(g, "_n_normal"), sum(m$condition == "normal", na.rm = TRUE))
    keynum(paste0(g, "_n_genes_std"), nrow(emat_std))
  }
  hlog("preprocessing complete")
})