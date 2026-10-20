## 10_clinical_validation.R - STEP 22: clinical validation in GSE42568 tumors
## Secondary validation only; does NOT alter the frozen candidate list.
source(file.path("scripts", "common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })

run_wrapped("10", {
  hlog("=== 10 CLINICAL VALIDATION (GSE42568) ===")
  frozen <- readRDS(file.path(LOG, "cache", "frozen_candidates.rds"))
  emat <- readRDS(file.path(LOG, "cache", "GSE42568_std_expr.rds"))
  meta <- readRDS(file.path(LOG, "cache", "GSE42568_meta.rds"))

  tumors <- meta$sample[meta$condition == "tumor"]
  tm <- meta[meta$sample %in% tumors, ]
  rownames(tm) <- tm$sample
  E <- emat[, tumors]

  genes <- intersect(frozen$Gene, rownames(E))
  nER <- sum(!is.na(tm$er)); nGR <- sum(!is.na(tm$grade)); nLN <- sum(!is.na(tm$lymph_node_status))
  nAge <- sum(!is.na(tm$age)); nSz <- sum(!is.na(tm$tumor_size_cm))
  nOS <- sum(!is.na(tm$os_event) & !is.na(tm$os_time_days)); nOSev <- sum(tm$os_event == 1, na.rm = TRUE)
  nRFS <- sum(!is.na(tm$rfs_event) & !is.na(tm$rfs_time_days)); nRFSev <- sum(tm$rfs_event == 1, na.rm = TRUE)
  hlog("GSE42568 tumors:", length(tumors),
       "; ER:", nER, "grade:", nGR, "LN:", nLN, "age:", nAge, "size:", nSz,
       "; OS events:", nOSev, "/", nOS, "; RFS events:", nRFSev, "/", nRFS)

  ## ---------- figures: ER & grade association for top candidates ----------
  topg <- frozen$Gene[seq_len(min(12, nrow(frozen)))]
  topg <- intersect(topg, rownames(E))
  if (nER >= 10 && length(topg)) {
    erd <- do.call(rbind, lapply(topg, function(ge) {
      data.frame(Gene = ge, ER = factor(ifelse(tm$er == 1, "positive", "negative"),
                                        levels = c("negative", "positive")), expr = E[ge, ])
    }))
    erd <- erd[!is.na(erd$ER), ]
    if (nrow(erd)) {
      p1 <- ggplot(erd, aes(ER, expr)) +
        geom_boxplot(aes(fill = ER), alpha = 0.7, width = 0.5) +
        facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c(negative = "#1b9e77", positive = "#d95f02")) +
        labs(title = "GSE42568: candidate expression by ER status", x = "", y = "Expression (log2)") +
        theme_bw(base_size = 9) + theme(legend.position = "none", strip.text = element_text(size = 7))
      save_gg(p1, "GSE42568_candidate_ER_association.png", w = 10, h = 7)
    }
  } else {
    hlog("GSE42568 ER association figure skipped: nER =", nER)
  }

  if (nGR >= 10 && length(topg)) {
    grd <- do.call(rbind, lapply(topg, function(ge) {
      data.frame(Gene = ge, Grade = factor(tm$grade), expr = E[ge, ])
    }))
    grd <- grd[!is.na(grd$Grade), ]
    if (nrow(grd)) {
      p2 <- ggplot(grd, aes(Grade, expr)) +
        geom_boxplot(aes(fill = Grade), alpha = 0.7, width = 0.6) +
        facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
        labs(title = "GSE42568: candidate expression by tumor grade", x = "Grade", y = "Expression (log2)") +
        theme_bw(base_size = 9) + theme(legend.position = "none", strip.text = element_text(size = 7))
      save_gg(p2, "GSE42568_candidate_grade_association.png", w = 10, h = 7)
    }
  } else {
    hlog("GSE42568 grade association figure skipped: nGR =", nGR)
  }



  rows <- list()
  for (ge in genes) {
    y <- E[ge, ]
    ## ER (1 vs 0)
    if (nER >= 10) {
      w <- suppressWarnings(wilcox.test(y[tm$er == 1], y[tm$er == 0]))
      rows[[paste(ge, "ER")]] <- data.frame(Gene = ge, Clinical_variable = "ER_status",
        Test = "Wilcoxon", Group1 = "ER+", Group2 = "ER-",
        Statistic = w$statistic, P_value = w$p.value, stringsAsFactors = FALSE)
    }
if (nGR >= 10) {
      kk <- suppressWarnings(kruskal.test(y ~ factor(tm$grade)))
      rows[[paste(ge, "GRADE")]] <- data.frame(Gene = ge, Clinical_variable = "Grade",
        Test = "Kruskal-Wallis", Group1 = "grade 1/2/3", Group2 = "",
        Statistic = kk$statistic, P_value = kk$p.value, stringsAsFactors = FALSE)
    }
    ## lymph node (1 vs 0)
    if (nLN >= 10) {
      w <- suppressWarnings(wilcox.test(y[tm$lymph_node_status == 1], y[tm$lymph_node_status == 0]))
      rows[[paste(ge, "LN")]] <- data.frame(Gene = ge, Clinical_variable = "Lymph_node_status",
        Test = "Wilcoxon", Group1 = "LN+", Group2 = "LN-",
        Statistic = w$statistic, P_value = w$p.value, stringsAsFactors = FALSE)
    }
    ## age (Spearman)
    if (nAge >= 10) {
      ct <- suppressWarnings(cor.test(y, tm$age, method = "spearman"))
      rows[[paste(ge, "AGE")]] <- data.frame(Gene = ge, Clinical_variable = "Age",
        Test ="Spearman correlation", Group1 = "", Group2 = "",
        Statistic = ct$estimate, P_value = ct$p.value, stringsAsFactors = FALSE)
    }
    ## size (Spearman)
    if (nSz >= 10) {
      ct <- suppressWarnings(cor.test(y, tm$tumor_size_cm, method = "spearman"))
      rows[[paste(ge, "SIZE")]] <- data.frame(Gene = ge, Clinical_variable = "Tumor_size",
        Test ="Spearman correlation", Group1 = "", Group2 = "",
        Statistic = ct$estimate, P_value = ct$p.value, stringsAsFactors = FALSE)
    }
  }
  clin <- do.call(rbind, rows)
  clin$FDR <- p.adjust(clin$P_value, "BH")
  write_tab(clin, "GSE42568_clinical_validation.csv")
  hlog("clinical validation rows:", nrow(clin), " significant (FDR<0.05):", sum(clin$FDR < 0.05))
## ---------- figures: ER & grade association for top candidates ----------
  topg <- frozen$Gene[seq_len(min(12,nrow(frozen)))]
  topg <- intersect(topg, rownames(E))
  erd <- do.call(rbind, lapply(topg, function(ge) {
    data.frame(Gene = ge, ER = factor(ifelse(tm$er == 1, "positive", "negative"),
                                             levels = c("negative", "positive")), expr = E[ge, ])
  }))
  erd <- erd[!is.na(erd$ER), ]
  p1 <- ggplot(erd, aes(ER, expr)) +
    geom_boxplot(aes(fill = ER), alpha = 0.7, width = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol =  4) +
    scale_fill_manual(values = c(negative = "#1b9e77", positive = "#d95f02")) +
    labs(title = "GSE42568: candidate expression by ER status", x = "", y = "Expression (log2)") +
    theme_bw(base_size =  9) + theme(legend.position = "none", strip.text = element_text(size =  7))
  save_gg(p1, "GSE42568_candidate_ER_association.png", w =  10, h =  7)

  grd <- do.call(rbind, lapply(topg, function(ge) {
    data.frame(Gene = ge, Grade = factor(tm$grade), expr = E[ge, ])
  }))
  grd <- grd[!is.na(grd$Grade), ]
  p2 <- ggplot(grd, aes(Grade, expr)) +
    geom_boxplot(aes(fill = Grade), alpha =0.7, width =0.6) +
    facet_wrap(~ Gene, scales = "free_y", ncol =  4) +
    labs(title = "GSE42568: candidate expression by tumor grade", x = "Grade", y = "Expression (log2)") +
    theme_bw(base_size =  9) + theme(legend.position = "none", strip.text = element_text(size =  7))
  save_gg(p2, "GSE42568_candidate_grade_association.png", w =  10, h =  7)
## ---------- survival (Cox, OS & RFS) ----------
  surv_rows <- list()
  if (nOS >=  30 && nOSev >=  5) {
    if (requireNamespace("survival", quietly = TRUE)) {
      suppressPackageStartupMessages(library(survival))
      osdf <- data.frame(time = tm$os_time_days, event = tm$os_event)
      rfsdf <- data.frame(time = tm$rfs_time_days, event = tm$rfs_event)
      for (ge in genes) {
        for (endpoint in c("OS", "RFS")) {
          d <- if (endpoint == "OS") osdf else rfsdf
          ok <- !is.na(d$time) & !is.na(d$event)
          if (sum(ok) < 30 || sum(d$event[ok]) < 5) next
          dd <- data.frame(time = d$time[ok], event = d$event[ok], expr = E[ge, ok])
          fit <- tryCatch(coxph(Surv(time, event) ~ expr, data = dd), error = function(e) NULL)
          if (is.null(fit)) next
          sf <- summary(fit)
          surv_rows[[paste(ge, endpoint)]] <- data.frame(
            Gene = ge, Endpoint = endpoint, n = sum(ok), Events = sum(d$event[ok]),
            HR = sf$conf.int[1, 1], HR_lower = sf$conf.int[1, 3], HR_upper = sf$conf.int[1,  4],
            P_value = sf$coefficients[1, 5], stringsAsFactors = FALSE)
        }
      }
      surv <- do.call(rbind, surv_rows)
      surv$FDR <- p.adjust(surv$P_value, "BH")
      write_tab(surv, "GSE42568_survival_analysis.csv")
      hlog("survival Cox rows:", nrow(surv), " FDR<0.05:", sum(surv$FDR < 0.05))
      ## KM curve for the most significant gene if any
      if (nrow(surv) && min(surv$P_value) < 0.05) {

        bg <- surv$Gene[which.min(surv$P_value)]
        ep <- surv$Endpoint[which.min(surv$P_value)]
        d <- if (ep == "OS") osdf else rfsdf
        dd0 <- data.frame(time = d$time, event = d$event, expr = E[bg, ])
        ok <- !is.na(dd0$time) & !is.na(dd0$event)
        dd <- dd0[ok, ]
        dd$group <- ifelse(dd$expr > median(dd$expr), "High", "Low")
        km_fit <- survfit(Surv(time, event) ~ group, data = dd)
        km_df <- data.frame(time = km_fit$time, surv = km_fit$surv,
                          group = rep(names(km_fit$strata), km_fit$strata))
        km_df$group <- gsub("group=", "", km_df$group)
        p3 <- ggplot(km_df, aes(time, surv, group = group, color = group)) +
          geom_step(lwd = 1.1) +
          scale_color_manual(values = c(Low = "#1b9e77", High = "#d95f02")) +
          labs(title = paste("GSE42568", ep, "by", bg, "expression"),
               x = "Days", y = "Survival probability") +
          theme_bw(base_size =  11)
        save_gg(p3, paste0("GSE42568_KM_", ep, "_", bg, ".png"), w =  6.5, h =  5)
      }
    } else {
      hlog("GSE42568 survival skipped: survival package unavailable")
    }
  } else {
    hlog("GSE42568 survival skipped: OS n=", nOS, " events=", nOSev, " < required (30 /  5)")
  }
  saveRDS(clin, file.path(LOG, "cache", "GSE42568_clinical.rds"))
  hlog("STEP  22 complete")
})
