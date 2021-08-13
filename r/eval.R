
#' Evaluate DOSE scores
#'
#' Provided one or several possible DOSE scores and one or several datasets
#' (e.g. data for different time windows) provided in indicator encoded
#' fashion, `dose_eval` computes the DOSE scores across all datasets.
#'
#' @param scores A single named numeric vector or a list of named numeric
#' vectors (the names of which are used to name the score columns)
#' @param dat A single `id_tbl` or a list of `id_tbl` datasets
#' @param label `NULL` or the name of the label column
#'
#' @return An `id_tbl` with a column per combination of score and dataset and
#' if provided a column containing the labels
#'
dose_eval <- function(scores, dat, label = "death") {

  eval_score <- function(score, name, dat, label) {

    eval_one <- function(dat, name) {

      assert_that(is_id_tbl(dat), has_name(dat, c(label, names(score))))

      score <- as.vector(
        as.matrix(dat[, names(score), with = FALSE]) %*% score
      )

      res <- dat[, c(id_var(dat), label), with = FALSE]
      res <- res[, c(name) := score]

      res
    }

    if (is_id_tbl(dat)) {
      eval_one(dat, name)
    } else {
      Reduce(merge_all, Map(eval_one, dat, paste(name, names(dat))))
    }
  }

  if (!is.list(scores)) {
    scores <- list(scores)
  }

  if (is.null(names(scores))) {

    if (length(scores) == 1L) {
      names(scores) <- "DOSE"
    } else {
      names(scores) <- paste0("DOSE_", seq_along(scores))
    }
  }

  Reduce(merge_all,
    Map(eval_score, scores, names(scores), MoreArgs = list(dat, label))
  )
}

#' Evaluate sofa scores
#'
#' For a given data source and a single or multiple time points, the SOFA
#' score is evaluated. As `cohort` argument, either a numeric vector of
#' patient ids or an `id_tbl` object as returned from [dose_eval()] may be
#' provided. In the latter case, this object is concatenated with the table
#' of SOFA scores and provides the label column which in the former case,
#' the `death` outcome is retrieved and added as label column.
#'
#' @param src String-valued data source name
#' @param upr Time-point(s) for which SOFA is evaluated
#' @param cohort Either a numeric vector of patient IDs or an `id_tbl`
#'
#' @return An `id_tbl` with columns per SOFA time-point, potentially merged
#' with the object passed as `cohort` argument.
#'
sofa_eval <- function(src, upr, cohort = si_cohort(src)) {

  extract_score <- function(x, name) {
    x <- x[, c(id_var(x), "sofa"), with = FALSE]
    x <- rename_cols(x, name, "sofa")
    x
  }

  merge_two <- function(x, y) merge(x, y, all = TRUE)

  if (is_id_tbl(cohort)) {
    out    <- cohort
    cohort <- out[[id_var(out)]]
  } else {
    out <- load_concepts("death", src, patient_ids = cohort)
  }

  res <- load_concepts("sofa", src, patient_ids = cohort, explicit_wins = upr)

  unt <- time_unit(res)
  res <- split(res, by = index_var(res), keep.by = FALSE)

  unt <- as.difftime(as.numeric(names(res)), units = unt)
  res <- Map(extract_score, res,
             paste0("SOFA [", format(unt - 24L), ", ", format(unt), "]"))

  Reduce(merge_all, c(res, list(out)))
}

#' Calculate ROC/PR
#'
#' For an `id_tbl` containing several scores (such as DOSE or SOFA),
#' potentially over multiple time windows, for each column vs the column
#' designated as `label`, resampled RO and PR curves are calculated.
#'
#' @param dat `id_tbl` holding scores and label as columns
#' @param lable String-valued name of the label column
#' @param n_rep Number of resampling steps
#' @param frac Fraction of data included in each sampling step
#'
#' @return See [precrec::evalmod()]
#'
eval_score <- function(dat, label = "death", n_rep = 10L, frac = 0.75) {

  extract_sco <- function(i, x, sco) c(x[i, sco, with = FALSE])
  extract_lab <- function(i, x, lab, n) rep(list(x[[lab]][i]), n)

  assert_that(is_id_tbl(dat), has_name(dat, label))

  scores <- setdiff(names(dat), c(id_var(dat), label))

  folds  <- replicate(n_rep, sample(nrow(dat), frac * nrow(dat)),
                      simplify = FALSE)

  splits <- lapply(folds, extract_sco, dat, scores)
  labels <- lapply(folds, extract_lab, dat, label, length(scores))

  precrec::evalmod(scores = splits, labels = labels,
                   modnames = rep(scores, n_rep),
                   dsids = rep(seq_len(n_rep), each = length(scores)))
}

merge_all <- function(x, y) merge_cols(x, y, all = TRUE)

merge_cols <- function(x, y, ...) {

  non_id_cols <- function(z) setdiff(colnames(z), id_var(z))

  common <- intersect(non_id_cols(x), non_id_cols(y))
  suffix <- c(".x", ".y")

  res <- merge(x, y, ..., suffixes = suffix)

  if (length(common) > 0L) {

    for (col in (common)) {

      cols <- paste0(col, suffix)

      if (!identical(res[[cols[1L]]], res[[cols[2L]]])) browser()

      assert_that(identical(res[[cols[1L]]], res[[cols[2L]]]))

      res[, c(col, cols) := list(get(cols[1L]), NULL, NULL)]
    }
  }

  res
}

eval_at <- function(eval_dat, type, nboot = 1000L, object = FALSE) {

  if (is.list(eval_dat) & !data.table::is.data.table(eval_dat)) {
    
    if (object) {
      return(lapply(eval_dat, eval_at, type = type, nboot = nboot, 
                    object = object))
    }
    
    return(cowplot::plot_grid(
        plotlist = lapply(eval_dat, eval_at, type = type, nboot = nboot, 
                          object = object),
        ncol = length(eval_dat)
      )
    )

  }

  methods <- setdiff(names(eval_dat), "death")

  trans <- list(four_c = "ISARIC 4C", age = "Age",
                random_forest = "Random Forest", naps = "RECOILS",
                sofa2 = "SOFA", saps_3 = "SAPS-III")
  methods_full <- lapply(1:length(methods), function(i) trans[[methods[i]]])

  if (nboot > 0L) {

    scores <- list()
    labels <- list()

    for(bt in 1:nboot) {

      ids <- sample(1:nrow(eval_dat), nrow(eval_dat), replace = T)
      scores <- c(scores,
                  lapply(methods, function(x) as.numeric(eval_dat[[x]][ids])))
      labels <- c(labels,
                  lapply(1:length(methods), function(x) eval_dat[["death"]][ids]))

    }

    dsids <- rep(1:nboot, each = length(methods))
    modnames <- rep(methods, nboot)

  } else {

    scores <- lapply(methods, function(x) as.numeric(eval_dat[[x]]))
    labels <- lapply(1:length(scores), function(i) eval_dat[["death"]])
    dsids <- 1:length(scores)
    modnames = methods[1:length(scores)]

  }
  
  if (object) return(list(scores = scores, labels = labels, dsids = dsids,
                          modnames = modnames))
  eval <- evalmod(scores = scores, labels = labels, dsids = dsids,
                  modnames = modnames)

  aurocs <- auc(eval)$aucs[c(T, F)]
  auprcs <- auc(eval)$aucs[c(F, T)]

  baseline <- which(methods_full == "RECOILS")

  aurocs_full <- tapply(aurocs, factor(modnames, levels = unique(modnames)),
                        function(x) x, simplify = F)
  names(aurocs_full) <- methods_full

  # more powerful, paired hypothesis testing
  pair.pvals <- rep(0, length(aurocs_full))
  if (length(baseline) == 1L) {

    base_aucs <- aurocs_full[[baseline]]

    for (i in 1:length(aurocs_full)) {

      cmp_aucs <- aurocs_full[[i]]
      pair.pvals[i] <- pnorm(mean(base_aucs - cmp_aucs) /
                             sd(base_aucs - cmp_aucs),
                             lower.tail = FALSE)

    }

  }

  aurocs <- tapply(aurocs, factor(modnames, levels = unique(modnames)),
                 function(x) {

                   c(mean(x), mean(x) - 1.96 * sd(x), mean(x) + 1.96 * sd(x),
                     sd(x))

                 }, simplify = F)

  auprcs <- tapply(auprcs, factor(modnames, levels = unique(modnames)),
                   function(x) {

                     c(mean(x), mean(x) - 1.96 * sd(x), mean(x) + 1.96 * sd(x))

                   }, simplify = F)
  
  if (length(methods_full) == 1L) {
    aurocs <- as.data.frame(t(Reduce(rbind, aurocs)), row.names = methods_full)
  } else {
    aurocs <- as.data.frame(Reduce(rbind, aurocs), row.names = methods_full)
  }
  names(aurocs) <- c("auc", "auc.min", "auc.max", "sd")

  if (length(baseline) == 1L) {

    sdb <- aurocs$sd[baseline]
    mb <- aurocs$auc[baseline]

    sds <- sqrt(sdb^2 + aurocs$sd^2)
    mdiff <- mb - aurocs$auc

    p.values <- pnorm(mdiff / sds, lower.tail = FALSE)

    aurocs <- cbind(aurocs, p.values = p.values)
    aurocs <- cbind(aurocs, pair.pvals = pair.pvals)

  }

  print(aurocs)
  print(auprcs)

  roc_labels <- paste0(methods_full, " (", round(aurocs$auc, 2), ")")
  
  if (length(methods_full) == 1L) {
    auprcs <- as.data.frame(t(Reduce(rbind, auprcs)), row.names = methods_full)
  } else {
    auprcs <- as.data.frame(Reduce(rbind, auprcs), row.names = methods_full)
  }
  names(auprcs) <- c("auc", "auc.min", "auc.max")

  pr_labels <- paste0(methods_full, " (", round(auprcs$auc, 2), ")")

  auprcs$aci <- paste0(round(auprcs$auc, 3), " (",
                       round(auprcs$auc.min, 3), "-",
                       round(auprcs$auc.max, 3), ")")
  aurocs$aci <- paste0(round(aurocs$auc, 3), " (",
                       round(aurocs$auc.min, 3), "-",
                       round(aurocs$auc.max, 3), ")")

  aurocs <- aurocs[, c("aci"), drop = FALSE]
  auprcs <- auprcs[, c("aci"), drop = FALSE]

  names(aurocs) <- "AUROC (95% CI)"
  names(auprcs) <- "AUPRC (95% CI)"

  aurocs$Score <- methods_full

  my_doc <- read_docx()

  my_doc <- my_doc %>%
    body_add_table(cbind(aurocs, auprcs), style = "table_template")

  print(my_doc, target = "~/aucs.docx")

  total_roc <- autoplot(eval, "ROC", show_cb = TRUE) +
    ggtitle("ROC curves at 24 hours into ICU") +
    theme(legend.position = c(0.7, 0.2),
          legend.text = element_text(size=10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")) +
    scale_color_discrete(labels = roc_labels)

  total_prc <- autoplot(eval, "PRC", show_cb = TRUE) +
    ggtitle("PR curves at 24 hours into ICU") +
    theme(legend.position = c(0.7, 0.8),
          legend.text = element_text(size=10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")) +
    scale_color_discrete(labels = pr_labels)
  
  cowplot::plot_grid(total_roc, total_prc, ncol = 1L) +
    theme(plot.background = element_rect(color = "black"))

}

eval_strat <- function(eval_dat, type, nboot = 1000L, lab_1 = "1st wave",
                       lab_2 = "2nd wave") {
  
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  res <- eval_at(eval_dat, nboot = nboot, object = TRUE)
  
  res[[1]][["dsids"]] <- res[[1]][["dsids"]] + nboot
  res[[1]][["modnames"]] <- gsub("naps", "RECOILS", paste(res[[1]][["modnames"]], 
                                                          lab_1))
  res[[2]][["modnames"]] <- gsub("naps", "RECOILS", paste(res[[2]][["modnames"]], 
                                                          lab_2))
  
  res <- Map(c, res[[1]], res[[2]])
  eval <- evalmod(scores = res[["scores"]], labels = res[["labels"]], 
                  dsids = res[["dsids"]], modnames = res[["modnames"]],
                  calc_avg = TRUE, cb_alpha = 0.05)
  aucc <- data.table(sc = res[["modnames"]], roc = auc(eval)$aucs[c(T, F)])
  aucc <- aucc[, list(mean(roc), 1.96 * sd(roc)), 
               by = "sc"]
  labels <- paste0(
    aucc$sc, " (", specify_decimal(aucc$V1, 3), " ± ", specify_decimal(aucc$V2, 3), ")"
  )
  
  autoplot(eval, "ROC", show_cb = TRUE) +
    ggtitle("ROC curves at 24 hours into ICU") +
    theme(legend.position = c(0.625, 0.15),
          legend.text = element_text(size=10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")) +
    scale_color_discrete(labels = labels)
  
}
