
#' Learn a DOSE score from data
#'
#' Using preprocessed data that is indicator encoded as returned by
#' [load_data()], a dose score is created from data. `train_dose()` serves to
#' set up data for, as well as check results of the function passed as `method`
#' argument, which carries out the actual score learning. The `method` function
#' is expected to take at least two arguments matrix `x` and vector `y` and
#' should return a named vector of length `ncol(x)` with the same names are
#' `colnames(x)`.
#'
#' @param dat An `id_tbl` object that holds indicator encoded data alongside
#' a response vector
#' @param response String valued name of the response
#' @param seed Integer valued seed used to fix e.g. data splitting in cross
#' validation
#' @param method Function that actually carries out learning of the score
#' @param ... Further arguments are passed to `emthod`
#'
#' @return A named numeric vector with the same attributes as the columns of
#' `dat` combined.
#'
train_dose <- function(dat, response = "death", seed = 2020L,
                       method = dose_glm, ...) {

  get_attrs <- function(x) {

    att <- attributes(x[[1L]])
    ats <- Map(vapply, list(x), list(attr), att, names(att))
    names(ats) <- names(att)

    lapply(ats, unname)
  }

  set_attrs <- function(x, attribs) {
    attributes(x) <- c(attributes(x), attribs)
    x
  }

  assert_that(is_id_tbl(dat), has_name(dat, response))

  x <- dat[, setdiff(colnames(dat), c(id_vars(dat), response)), with = FALSE]

  attribs <- get_attrs(x)

  set.seed(seed)

  res <- method(set_attrs(as.matrix(x), attribs), dat[[response]], ...)

  #assert_that(is.numeric(res), identical(names(res), colnames(x)))

  #set_attrs(res, attribs)

  res
}

dose_glm <- function(x, y, model_ind = 22L, model_ind2 = NULL, coef_frac = 0.1, n_folds = 5L,
                     n_cores = get_cores(), zoom_lambda = F, ...) {

  as_score <- function(x) round(x / min(abs(x)[abs(x) > 0]))

  coefs <- function(x) coef(x[["glmnet.fit"]])[-1L, ]

  plot_glmnet <- function(fit, counts, model = NULL) {

    dat <- data.frame(index = seq_along(counts), counts = counts,
                      auc = fit$cvm, df = fit$glmnet.fit$df)

    p1 <- ggplot(dat, aes(counts, auc)) + geom_point() +
      labs(x = "Complexity", y = "Cross-validated AUC",
           title = "Complexity-AUC path") + theme_bw()
    p2 <- ggplot(dat, aes(index, counts)) + geom_point() +
      labs(x = "Index", y = "Counts", title = "Features used") + theme_bw()
    p3 <- ggplot(dat, aes(index, df)) + geom_point() +
      labs(x = "Index", y = "DF", title = "Degrees of freedom") + theme_bw()
    p4 <- ggplot(dat, aes(index, auc)) + geom_point() +
      labs(x = "Index", y = "AUC", title = "Cross-validated AUC") + theme_bw()

    if (!is.null(model)) {
      hit <- dat[dat[["index"]] == model, ]
      p1 <- p1 + geom_vline(xintercept = hit[["counts"]]) +
                 geom_hline(yintercept = hit[["auc"]])
      p2 <- p2 + geom_vline(xintercept = model)
      p3 <- p3 + geom_vline(xintercept = model)
      p4 <- p4 + geom_vline(xintercept = model)
    }

    print(cowplot::plot_grid(p1, p2, p3, p4))
  }

  auto_mod <- function(fit, counts) {
    opts <- counts < 20
    counts[fit$cvm < quantile(fit$cvm[opts], probs = 0.8) | !opts] <- Inf
    which.min(counts)
  }

  if (n_cores > 1L) {
    doMC::registerDoMC(cores = n_cores)
  }

  fit <- glmnet::cv.glmnet(x, y, family = "binomial", type.measure = "auc",
                           parallel = n_cores > 1L, nfolds = n_folds, ...)

  counts <- lengths(
    apply(coefs(fit), 2L, which_feats, coef_frac, attr(x, "concept"))
  )

  assert_that(is.count(model_ind))

  res <- coefs(fit)[, model_ind]
  res[!gt_thresh(res, coef_frac)] <- 0

  list(
    fit = fit, counts = counts, plot = plot_glmnet(fit, counts, model_ind),
    score = as_score(res)
  )
}

gt_thresh <- function(col, frac) abs(col) >= (max(abs(col)) * frac)

which_feats <- function(col, frac, feat) unique(feat[gt_thresh(col, frac)])
