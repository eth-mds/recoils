
#' Load preprocessed data
#'
#' For a given data source (alongside a config file, and a cohort) and for a
#' single time window or multiple time windows returns preprocessed data. The
#' Returned data first is preprocessed by [preproc()] and subsequently
#' indicator encoded, using [indicator_encoding()]. Finally the `death`
#' outcome is merged in.
#'
#' @param src String valued data source
#' @param cfg A named list containing a slot per data column, each with
#' entries `direction`, `upper`, `lower` and `step`
#' @param lwr,upr Lower and upper bound of data windows; the length of one has
#' to be a multiple of the other
#' @param cohort Vector of patient ids
#'
#' @return Either a single or a list of `id_tbl` objects, depending on whether
#' a single or multiple time windows were specified.
#'
load_data <- function(src, cfg, lwr, upr, cohort = NULL, impute_vals = T, cts = F) {

  load_win <- function(lwr, upr, cfg, dat, out, impute_vals, cts) {

    dat <- preproc(dat, cfg, lwr, upr, impute_vals)
    if (!cts) dat <- indicator_encoding(dat, cfg)

    ret <- merge(dat, out, all.x = T)
    ret[is.na(death), "death"] <- F

    ret
  }

  load_wins <- function(lwr, upr, cfg, dat, out, impute_vals, cts) {
    Map(function(a, b) {
      if (is_id_tbl(a) & is_id_tbl(b)) {
        return(load_win(a, b, cfg, dat, out, impute_vals, cts))
      } else {
        return(load_win(as.difftime(a, units = units(lwr)),
                 as.difftime(b, units = units(upr)), cfg, dat, out, impute_vals, cts))
      }
    }, lwr, upr)
  }

  load_from <- FALSE

  if (file.exists("/home/drago/dat.RData")) {

    load("/home/drago/dat.RData")
    if (!all(names(cfg) %in% names(dat))) {

      dat <- load_concepts(names(cfg), src, aggregate = aggreg_fun(cfg))
      save(dat, file = "/home/drago/dat.RData")

    }

  } else {

    dat <- load_concepts(names(cfg), src, aggregate = aggreg_fun(cfg))
    if (src == "covid19") save(dat, file = "/home/drago/dat.RData")

  }

  dat <- dat[get(id_vars(dat)) %in% cohort]
  dat <- dat[, c(meta_vars(dat), names(cfg)), with=F]

  out <- load_concepts("death", src, patient_ids = cohort)
  out[, c(index_var(out)) := NULL]

  res <- load_wins(lwr, upr, cfg, dat, out, impute_vals, cts)

  if (!is.list(lwr)) {
    names(res) <- paste0(ifelse(is.finite(lwr), "[", "("), format(lwr), ", ",
                         format(upr), ifelse(is.finite(upr), "]", ")"))
  } else {
    names(res) <- "Mech_vent_explicit_wins"
  }

  ids <- Map(`[[`, res, lapply(res, id_var))
  ids <- Reduce(intersect, ids)
  res <- lapply(res, function(x) x[get(id_vars(x)[1L]) %in% ids])
  ids <- Map(`[[`, res, lapply(res, id_var))

  assert_that(all(vapply(ids, identical, logical(1L), ids[[1L]])))

  res
}

#' Data preprocessing
#'
#' Fist, rows are filtered out where the index is strictly smaller than the
#' `win_lwr` or strictly larger than the `win_upr` argument. From this reduced
#' data, feature-wise medians are calculated for imputing `NA` data later on.
#' Per id group, the worst value is calculated according to the `cfg` argument
#' which for the value `increasing` means a maximum and for `decreasing` a
#' minimum. If not a single value is available for a given feature and id
#' group, the median value for this feature calculated earlier (over all rows)
#' is imputed.
#'
#' @param dat A `ts_tbl` object holding patient data
#' @param cfg Named list (names have to correspond to feature columns in
#' `dat`), each slot containing a `direction` string, indicating the direction
#' of the given feature
#' @param win_lwr,win_upr Filter criteria for rows
#'
#' @return An `id_tbl` object
#'
preproc <- function(dat, cfg, win_lwr = hours(-Inf),
                    win_upr = hours(Inf), impute_vals = TRUE) {

  do_call <- function(fun, x) do.call(fun, list(x))

  repl_na <- function(x, val) replace(x, is.na(x), val)

  assert_that(is_ts_tbl(dat), is.list(cfg), all(data_vars(dat) %in% names(cfg)))

  cfg <- cfg[data_vars(dat)]

  agg <- aggreg_fun(cfg, list(max_or_na), list(min_or_na))

  med <- lapply(cfg, function(x) {

      if(is.null(x[["impute_val"]]))
        return(
          ifelse(x[["direction"]] == "increasing", x[["lower"]]*0.99, x[["upper"]]*1.01)
        )

      x[["impute_val"]]

  })

  if(is_id_tbl(win_lwr) & is_id_tbl(win_upr)) {

    assert_that("min_time" %in% names(win_lwr), "max_time" %in% names(win_upr))

    win_lwr <- merge(dat[, id_vars(dat), with=F],
      win_lwr, by = id_vars(dat), all.x = T)[["min_time"]]
    win_upr <- merge(dat[, id_vars(dat), with=F],
      win_upr, by = id_vars(dat), all.x = T)[["max_time"]]

  }

  res <- dat[(get(index_var(dat)) >= win_lwr) & (get(index_var(dat)) <= win_upr), ]
  res <- res[, Map(do_call, agg, .SD), .SDcols = names(agg), by = c(id_vars(dat))]
  if (impute_vals)
    res <- res[, c(names(med)) := Map(repl_na, .SD, med), .SDcols = names(med)]

  as_id_tbl(res, id_vars(dat))
}

#' Convert to indicator encoded data
#'
#' Turn an `id_tbl` object containing numeric data columns into an indicator
#' encoded `id_tbl` object consisting of logical columns according to intervals
#' as specified by the configuration list passed as `cfg`. Each slot in `cfg`
#' is expected to contain an entry for `direction`, alongside a discretization
#' specification as `upper`, `lower` and `step`.
#'
#' @param dat An `id_tbl` object containing numeric data
#' @param cfg A named list containing a slot per data column, each with
#' entries `direction`, `upper`, `lower` and `step`.
#'
#' @return An `id_tbl` object
#'
indicator_encoding <- function(dat, cfg) {

  encode <- function(x, sequ, is_inc, name, len) {
    if (is.null(x)) {
      ival <- cut(sequ, sequ, right = !is_inc)
    } else {
      ival <- cut(x, sequ, right = !is_inc)
    }

    if (is_inc) {
      col_names <- sub(",\\d+(\\.\\d+)?\\)$", ",Inf)", levels(ival))
    } else {
      #browser()
      col_names <- sub("^\\(\\d+(\\.\\d+)?,", "(-Inf,", levels(ival))
    }


    res <- matrix(FALSE, nrow = len, ncol = nlevels(ival),
                  dimnames = list(NULL, col_names))

    if (!is.null(x)) {

      if (is_inc) {

        seq_fun <- function(i) seq.int(1L, i)
        res[!is.na(x) & x >= cfg[[name]][["upper"]], ] <- TRUE

      } else {
        seq_fun <- function(i) seq.int(i, ncol(res))
        res[!is.na(x) & x <= cfg[[name]][["lower"]], ] <- TRUE
      }

      lvls <- as.integer(ival)

      for (i in seq_len(nlevels(ival))) {
        res[!is.na(lvls) & lvls == i, seq_fun(i)] <- TRUE
      }
    }

    res <- asplit(res, 2L)
    res <- lapply(res, `dimnames<-`, NULL)
    res <- lapply(res, `attr<-`, "concept", name)
    res <- Map(`attr<-`, res, "threshold",
               if (is_inc) sequ[-length(sequ)] else sequ[-1])
    res <- lapply(res, `attr<-`, "right", is_inc)

    res
  }

  assert_that(is_id_tbl(dat), is.list(cfg),
              all(data_vars(dat) %in% names(cfg)))

  seqs <- Map(seq,
    vapply(cfg, `[[`, numeric(1L), "lower"),
    vapply(cfg, `[[`, numeric(1L), "upper"),
    vapply(cfg, `[[`, numeric(1L), "step")
  )

  res <- Map(encode, c(dat)[names(cfg)], seqs, aggreg_fun(cfg, TRUE, FALSE),
             names(cfg), nrow(dat))
  names(res) <- names(cfg)

  res <- data.table::setDT(unlist(res, recursive = FALSE))
  res <- cbind(dat[[id_vars(dat)]], res)
  res <- data.table::setnames(res, "V1", id_vars(dat))

  as_id_tbl(res, id_vars(dat))
}

set_design <- function(src = "covid19", type = "icu24", output = "mat",
                       pids = config(paste0(type, "-clean"))[[src]],
                       add_cmp = TRUE, cfg = get_config("features4", config_dir()),
                       ...) {

  if (type == "mechvent") {

    upr <- load_concepts("vent_ind2", src, patient_ids = pids)
    upr <- upr[, head(.SD, n = 1L), by = c(id_vars(upr))]
    upr[, max_time := get(index_var(upr))]
    upr <- upr[, c(id_vars(upr), "max_time"), with=F]

    lwr <- data.table::copy(upr)
    lwr[, min_time := max_time - hours(48L)]
    lwr[, max_time := NULL]

    cohort <- unique(id_col(upr))
    lwr <- list(lwr)
    upr <- list(upr)

  } else if (type == "icu24"){

    upr <- hours(c(24L, 48L))
    lwr <- upr - hours(48L)
    cohort <- config("8-of-10")[[src]]

  } else if (type == "conclusive") {

    upr <- hours(c(24L, 48L))
    lwr <- upr - hours(48L)
    cohort <- config("conclusive")[[src]]

  } else return(cat(type, " does not exist"))

  dat <- load_data(src, cfg, lwr, upr, cohort = cohort, cts = (output == "cts"),
                   ...)

  ### add all baselines
  if (add_cmp) {

    # adding Random Forest
    rf_data <- load_data(src, get_config("features2", config_dir()), lwr, upr,
                         cohort = id_col(dat[[1L]]), cts = TRUE)
    form <- as.formula(paste("death ~ . -", id_var(dat[[1L]])))

    rf_pred <- lapply(rf_data, function(x) {

      rf <- ranger(form, data = x, probability = TRUE, importance = "impurity")
      cat("\n\n\n")
      print(sort(importance(rf)))

      rf$predictions[, "TRUE"]

    })

    dat <- Map(function(x, y) {

      x[, random_forest := y]
      x
    }, dat, rf_pred)

    # add sofa
    load(paste0("~/sofa_", type, ".RData"))
    dat <- Map(function(x, y) merge(x, y, all.x = TRUE), dat, sofa_baseline)

    # add ISARIC4C
    load(paste0("~/four_c_", type, ".RData"))
    dat <- Map(function(x, y) merge(x, y, all.x = TRUE), dat, four_c_baseline)

    # add SAPS-III
    # add ISARIC4C
    load(paste0("~/saps_3_", type, ".RData"))
    saps_baseline <- lapply(saps_baseline, function(x) x[, c(meta_vars(x), "saps_3"),
                                                         with=FALSE])
    dat <- Map(function(x, y) merge(x, y, all.x = TRUE), dat, saps_baseline)

    # add age if needed
    if (output == "mat") {

      dat <- lapply(dat, function(dt) {

          res <- merge(dt, load_concepts("age", src), all.x = T)
          res <- replace_na(res, 65L)

          res

       })

    }

  }

  if (output == "mat") {

    X <- as.matrix(dat[[1L]][, setdiff(names(dat[[1L]]),
                                       c(meta_vars(dat[[1L]]), "death", "age",
                                         "min_time", "max_time", "saps_3",
                                         "random_forest", "sofa2", "four_c")),
                             with=FALSE])
    y <- dat[[1L]][["death"]]
    ids <- dat[[1L]][, c(id_vars(dat[[1L]])), with=FALSE]
    ids <- merge(ids, load_concepts(c("hosp", "wave"), src), all.x = T)

    return(list(X = X, y = y, ids = ids, design = dat))

  }

  dat

}
