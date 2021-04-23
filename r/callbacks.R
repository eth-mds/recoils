covid19_death <- function(x, ...) {
  x[, gender := T]
}

covid19_conclusive <- function(x, val_var, extra_var, ...) {
  x[, c(val_var) := (!get(val_var) & !get(extra_var))]
  x[, c(extra_var) := NULL]

  x
}

span_window <- function (x, val_var, end_var, ...)
{
  hr <- `units<-`(hours(1L), time_unit(x))
  idx <- index_var(x)
  x <- x[get(end_var) - get(idx) >= 0, ]
  x <- x[get(end_var) - get(idx) == 0, `:=`(c(end_var), get(idx) +
      hr)]
  x <- x[, `:=`(c(val_var), T)]
  res <- expand(x, idx, end_var, keep_vars = c(id_vars(x),
    val_var))
  res
}

safi <- function (..., match_win = hours(2L), mode = c("match_vals",
  "extreme_vals", "fill_gaps"), fix_na_fio2 = TRUE, interval = NULL)
{
  mode <- match.arg(mode)
  cnc <- c("o2sat", "fio2")
  res <- ricu:::collect_dots(cnc, interval, ...)
  assert_that(is_interval(match_win), match_win > ricu:::check_interval(res),
    is.flag(fix_na_fio2))
  if (identical(mode, "match_vals")) {
    on12 <- paste(meta_vars(res[[1L]]), "==", meta_vars(res[[2L]]))
    on21 <- paste(meta_vars(res[[2L]]), "==", meta_vars(res[[1L]]))
    res <- rbind(res[[1L]][res[[2L]], on = on12, roll = match_win],
      res[[2L]][res[[1L]], on = on21, roll = match_win])
    res <- unique(res)
  }
  else {
    res <- reduce(merge, res, all = TRUE)
    if (identical(mode, "fill_gaps")) {
      res <- fill_gaps(res)
    }
    win_expr <- substitute(list(po2 = min_fun(get(cnc[1L])),
      fio2 = max_fun(get(cnc[2L]))), list(min_fun = min_or_na,
        max_fun = max_or_na))
    res <- slide(res, !!win_expr, before = match_win, full_window = FALSE)
  }
  if (fix_na_fio2) {
    res <- res[is.na(get(cnc[2L])), `:=`(c(cnc[2L]), 21)]
  }
  res <- res[!is.na(get(cnc[1L])) & !is.na(get(cnc[2L])) &
      get(cnc[2L]) != 0, ]
  res <- res[, `:=`(c("safi"), 100 * get(cnc[1L])/get(cnc[2L]))]
  res <- rm_cols(res, cnc)
  res
}

nlr_callback <- function (..., match_win = hours(6L), interval = NULL) {

  cnc <- c("neut", "lymph")
  res <- ricu:::collect_dots(cnc, interval, ...)
  assert_that(is_interval(match_win), match_win > ricu:::check_interval(res))

  on12 <- paste(meta_vars(res[[1L]]), "==", meta_vars(res[[2L]]))
  on21 <- paste(meta_vars(res[[2L]]), "==", meta_vars(res[[1L]]))

  res <- rbind(res[[1L]][res[[2L]], on = on12, roll = match_win],
    res[[2L]][res[[1L]], on = on21, roll = match_win])
  res <- unique(res)

  res <- res[!is.na(get(cnc[1L])) & !is.na(get(cnc[2L])), ]
  res <- res[, `:=`(c("nlr"), get(cnc[1L])/get(cnc[2L]))]
  res <- rm_cols(res, cnc)
  res
}

los_callback <- function(x, id_type, interval) {

  as_day <- function(x) as.double(x, units = "days")

  win <- x[["win_type"]]
  cfg <- as_id_cfg(x)

  if (identical(win, id_type)) {

    res <- id_map(x, id_vars(cfg[id_type]), id_vars(cfg[win]), NULL, "end")

    res <- res[, c("val_var", "end") := list(as_day(get("end")), NULL)]

  } else {

    res <- id_map(x, id_vars(cfg[id_type]), id_vars(cfg[win]), "start", "end")

    res <- res[, c("val_var", "start", "end") := list(
      as_day(get("end") - get("start")), NULL, NULL
    )]

    res <- rm_cols(res, id_vars(cfg[win]), by_ref = TRUE)

    if (cfg[win] > cfg[id_type]) {
      res <- unique(res)
    }
  }

  res
}

los_hosp <- function(x, val_var, end_var, ...) {

  as_day <- function(x) as.double(x, units = "days")

  x <- x[, c(val_var) := (as_day(get(end_var)) - as_day(get(val_var)))]

  x

}

los_icu <- function(x, val_var, end_var, ...) {

  as_day <- function(x) as.double(x, units = "days")

  x <- x[, c(val_var) := (as_day(get(end_var)) - as_day(get(val_var)))]

  x

}

four_c_score_callback <- function (..., win_length = hours(48L),
                                   interval = NULL,
                                   keep_components = FALSE,
                                   explicit_wins = FALSE)  {

  score_calc <- function(age, sex, n_comb, resp, o2sat, gcs, bun, crp) {

    na_to_zero <- function(x) {
      x[is.na(x)] <- 0
      x
    }

    eval_vec <- function(x, vals, breaks) na_to_zero(vals[.bincode(x, breaks)])

    pts <- list(
      age = list(c(0, 2, 4, 6, 7), c(0, 50, 60, 70, 80, Inf)),
      n_comb = list(c(0, 1, 2), c(-0.01, 0.5, 1.5, Inf)),
      resp = list(c(0, 1, 2), c(0, 20, 30, Inf)),
      o2sat = list(c(2, 0), c(0, 92, Inf)),
      gcs = list(c(2, 0), c(0, 14.5, Inf)),
      bun = list(c(0, 1, 3), c(0, 19.6, 39.2, Inf)),
      crp = list(c(0, 1, 2), c(0, 50, 100, Inf))
    )

    score <- na_to_zero(sex == "Male")
    comps <- c("age", "n_comb", "resp", "o2sat", "gcs", "bun", "crp")
    for (cmp in comps) {

      score <- score + eval_vec(get(cmp), pts[[cmp]][[1]], pts[[cmp]][[2]])

    }

    score

  }

  cnc <- c("age", "sex", "n_comb", "resp", "o2sat", "gcs", "bun", "crp")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)

  cnc_lb <- c("resp", "o2sat", "gcs", "bun", "crp")
  cnc_const <- c("age", "sex", "n_comb")
  worst_val_funs <- list(max_or_na, min_or_na, min_or_na, max_or_na, max_or_na)
  app_fun <- function(f, x) f(x)
  expr <- substitute(Map(app_fun, fun, .SD), list(fun = worst_val_funs))

  if (isFALSE(explicit_wins)) {
    res <- fill_gaps(dat)
    res <- slide(res, !!expr, before = win_length, full_window = FALSE,
                 .SDcols = cnc_lb)
  } else {

    if (is_id_tbl(explicit_wins)) {

      res <- hop(dat, !!expr, windows = explicit_wins, .SDcols = cnc_lb)

    } else {

      res <- slide_index(dat, !!expr, explicit_wins, before = win_length,
                         full_window = FALSE, .SDcols = cnc_lb)

    }

  }

  res <- rename_cols(res, cnc_lb, paste0("V", 1:length(cnc_lb)))

  dat <- dat[, c(id_vars(dat), cnc_const), with = F]
  dat[, n_comb := max_or_na(n_comb), by = c(id_vars(dat))]
  dat[, n_comb := data.table::fifelse(is.na(n_comb), 0L, n_comb)]
  dat <- dat[complete.cases(dat)]
  dat <- unique(dat)
  dat <- merge(res, dat)

  if("max_time" %in% names(dat)) dat <- as_ts_tbl(dat, index_var = "max_time")

  dat <- dat[, `:=`(c("four_c"), score_calc(get("age"),
    get("sex"), get("n_comb"), get("resp"), get("o2sat"),
    get("gcs"), get("bun"), get("crp")))]
  if (isFALSE(keep_components)) {
    extra_cols <- setdiff(names(dat), c(meta_vars(dat), "four_c"))
    dat <- rm_cols(dat, extra_cols, by_ref = TRUE)
  }
  dat
}

saps_3_callback <- function (..., win_length = hours(48L),
                                   interval = NULL,
                                   keep_components = FALSE,
                                   explicit_wins = FALSE)  {

  score_calc <- function(age, los_bef_icu, norepi_rate, gcs, bili, temp, crea,
                         hr, wbc, ph,
                         plt, sbp, pafi, po2, vent_ind2) {

    na_to_zero <- function(x) {
      x[is.na(x)] <- 0
      x
    }

    pafi[is.na(pafi)] <- 476.19

    eval_vec <- function(x, vals, breaks) na_to_zero(vals[.bincode(x, breaks)])

    pts <- list(
      age = list(c(0, 5, 9, 13, 15, 18), c(0, 40, 60, 70, 75, 80, Inf)),
      los_bef_icu = list(c(7, 6, 0), c(-Inf, -672, -336, Inf)),
      norepi_rate = list(c(0, 3), c(-0.01, 0.000001, Inf)),
      gcs = list(c(15, 10, 7, 2, 0), c(0, 4.5, 5.5, 6.5, 12.5, Inf)),
      bili = list(c(0, 4, 5), c(0, 2, 6, Inf)),
      temp = list(c(7, 0), c(0, 35, Inf)),
      crea = list(c(0, 2, 7, 8), c(0, 1.2, 2, 3.5, Inf)),
      hr = list(c(0, 5, 7), c(0, 120, 160, Inf)),
      wbc = list(c(0, 2), c(0, 15, Inf)),
      ph = list(c(3, 0), c(0, 7.25, Inf)),
      plt = list(c(13, 8, 5, 0), c(0, 20, 50, 100, Inf)),
      sbp = list(c(11, 8, 3, 0), c(0, 40, 70, 120, Inf)),
      pafi = list(c(11, 7), c(0, 100, Inf)),
      po2 = list(c(5, 0), c(0, 60, Inf))
    )

    score <- 0L
    comps <- c("age", "los_bef_icu", "norepi_rate",
               "gcs", "bili", "temp", "crea", "hr", "wbc", "ph", "plt", "sbp",
               "oxygenation")

    for (cmp in comps) {

      if (cmp == "oxygenation") {

        pafi_eval <- eval_vec(pafi, pts[["pafi"]][[1]], pts[["pafi"]][[2]])
        po2_eval <- eval_vec(po2, pts[["po2"]][[1]], pts[["po2"]][[2]])

        score <- score + na_to_zero(ifelse(ricu:::is_true(vent_ind2), pafi_eval,
                                            po2_eval))

      } else score <- score + eval_vec(get(cmp), pts[[cmp]][[1]],
                                       pts[[cmp]][[2]])

    }

    score

  }

  cnc_lb <- c("norepi_rate",
              "gcs", "bili", "temp", "crea", "hr", "wbc", "ph", "plt", "sbp",
              "pafi", "po2", "vent_ind2")
  cnc_const <- c("age", "los_bef_icu")
  cnc <- c(cnc_const, cnc_lb)
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)

  worst_val_funs <- list(max_or_na, min_or_na, max_or_na, max_or_na, max_or_na,
                         max_or_na, max_or_na, min_or_na, min_or_na, min_or_na,
                         min_or_na, min_or_na, any)
  app_fun <- function(f, x) f(x)
  expr <- substitute(Map(app_fun, fun, .SD), list(fun = worst_val_funs))

  if (isFALSE(explicit_wins)) {
    res <- fill_gaps(dat)
    res <- slide(res, !!expr, before = win_length, full_window = FALSE,
                 .SDcols = cnc_lb)

  } else {

    if (is_id_tbl(explicit_wins)) {

      res <- hop(dat, !!expr, windows = explicit_wins, .SDcols = cnc_lb)

    } else {

      res <- slide_index(dat, !!expr, explicit_wins, before = win_length,
                         full_window = FALSE, .SDcols = cnc_lb)

    }

  }

  res <- rename_cols(res, cnc_lb, paste0("V", 1:length(cnc_lb)))

  dat <- dat[, c(id_vars(dat), cnc_const), with = F]
  #dat[, n_comb := max_or_na(n_comb), by = c(id_vars(dat))]
  #dat[, n_comb := data.table::fifelse(is.na(n_comb), 0L, n_comb)]
  dat[is.na(los_bef_icu), "los_bef_icu"] <- hours(0L)
  dat <- dat[complete.cases(dat)]
  dat <- unique(dat)
  dat <- merge(res, dat)

  if("max_time" %in% names(dat)) dat <- as_ts_tbl(dat, index_var = "max_time")
  dat <- dat[, `:=`(c("saps_3"), score_calc(get("age"), get("los_bef_icu"),
                                            get("norepi_rate"),
                                            get("gcs"), get("bili"),
                                            get("temp"), get("crea"),
                                            get("hr"), get("wbc"),
                                            get("ph"), get("plt"),
                                            get("sbp"), get("pafi"),
                                            get("po2"), get("vent_ind2")))]
  if (isFALSE(keep_components)) {
    extra_cols <- setdiff(names(dat), c(meta_vars(dat), "saps_3"))
    dat <- rm_cols(dat, extra_cols, by_ref = TRUE)
  }
  dat
}

grs_callback <- function (..., win_length = hours(48L), interval = NULL,
  keep_components = FALSE, explicit_wins = FALSE)  {

  score_calc <- function(age, sex, o2sat, neut_total, crp, alb, crea, DM, COPD) {

    na_to_zero <- function(x) {
      x[is.na(x)] <- 0
      x
    }

    score <- na_to_zero(age > 0) + na_to_zero(sex == "Male") +
      na_to_zero(o2sat < 93) + na_to_zero(neut_total > 8) +
      na_to_zero(crp > 40) + na_to_zero(alb < 3.4) + na_to_zero(crea > 1.13) +
      na_to_zero(DM) + na_to_zero(COPD)

    score

  }

  cnc <- c("age", "sex", "o2sat", "neut_total", "crp", "alb", "crea", "DM",
           "COPD")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)

  cnc_lb <- c("o2sat", "neut_total", "crp", "alb", "crea")
  cnc_const <- c("age", "sex", "DM", "COPD")
  worst_val_funs <- list(min_or_na, max_or_na, max_or_na, min_or_na, max_or_na)
  app_fun <- function(f, x) f(x)
  expr <- substitute(Map(app_fun, fun, .SD), list(fun = worst_val_funs))

  if (isFALSE(explicit_wins)) {
    res <- fill_gaps(dat)
    res <- slide(res, !!expr, before = win_length, full_window = FALSE,
                 .SDcols = cnc_lb)
  } else {

    if (is_id_tbl(explicit_wins)) {

      res <- hop(dat, !!expr, windows = explicit_wins, .SDcols = cnc_lb)

    } else {

      res <- slide_index(dat, !!expr, explicit_wins, before = win_length,
                         full_window = FALSE, .SDcols = cnc_lb)

    }

  }

  res <- rename_cols(res, cnc_lb, paste0("V", 1:length(cnc_lb)))

  dat <- dat[, c(id_vars(dat), cnc_const), with = F]
  dat <- dat[complete.cases(dat)]
  dat <- unique(dat)
  dat <- merge(res, dat)

  if("max_time" %in% names(dat)) dat <- as_ts_tbl(dat, index_var = "max_time")


  dat <- dat[, `:=`(c("grs"), score_calc(get("age"),
    get("sex"), get("o2sat"), get("neut_total"), get("crp"),
    get("alb"), get("crea"), get("DM"), get("COPD")))]
  if (isFALSE(keep_components)) {
    extra_cols <- setdiff(names(dat), c(meta_vars(dat), "grs"))
    dat <- rm_cols(dat, extra_cols, by_ref = TRUE)
  }
  dat
}

blood_cell_callback <- function (x, val_var, env, ...)
{
  x <- add_concept(x, env, "wbc")
  x <- x[, `:=`(c(val_var, "wbc"), list(100 * get(val_var)/get("wbc"),
    NULL))]
  x
}

wave_callback <- function(x, val_var, env, ...) {

  real_time <- load_src("episodes", cols = c("episode_id",
                                             "admission_timestamp"),
    src = env)

  x <- merge(x, real_time, by = "episode_id")
  x[, c(val_var) := ifelse(admission_timestamp <
                             as.POSIXlt("2020-08-01 00:00:00"), "1st", "2nd")]

}

cmbd_callback <- function(x, val_var, ...) {

  x[, n_comb := rowSums(.SD, na.rm=T), by = c(id_vars(x))]

  x <- x[, c(id_vars(x), "n_comb"), with=F]

  rename_cols(x, val_var, "n_comb")
}

nice_callback <- function(x, val_var, ...) {
  x[, c(val_var) := ricu:::is_true(get(val_var))]
}

vent_total_callback <- function(vent_ind2, ...) {

  res <- vent_ind2[, list(vent_total = sum(vent_ind2)),
                   by = c(id_var(vent_ind2))]

  res
}

is_vent_callback <- function(vent_ind2, interval, ...) {

  vent_ind2[, c(index_var(vent_ind2)) := NULL]
  vent_ind2 <- unique(vent_ind2)

  vent_ind2[, is_vent := vent_ind2]
  vent_ind2[, vent_ind2 := NULL]

  vent_ind2

}

DM_callback <- function(x, ...) {

  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {

    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]

  }

  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(icd::comorbid_charlson(intm)[, c("DM", "DMcx")]) > 0

  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )

  names(res) <- names(x)

  res
}

CPD_callback <- function(x, ...) {

  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {

    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]

  }

  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(icd::comorbid_charlson(intm)[, c("Pulmonary"), drop = F]) > 0

  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )

  names(res) <- names(x)
  res
}

norepi_rate_kg <- function (x, val_var, stop_var, env, ...)
{

  x <- add_concept(x, env, "sex")
  shift <- as.double(hours(1L), units = units(interval(x)))
  x[, c(stop_var) := data.table::fifelse(is.na(get(stop_var)),
                                         get(index_var(x)) + shift,
                                         get(stop_var))]

  x <- unique(x)

  x[is.na(get(val_var)), c(val_var) := 0.0001]

  x[, w_norm := data.table::fifelse(sex == "Male", 80, 60)]

  x[, c(val_var) := 1000 * get(val_var) /
                    (w_norm * as.double(get(stop_var) - get(index_var(x)),
                                             units = "mins"))]

  expand(x, index_var(x), stop_var, keep_vars = c(id_vars(x), val_var))
}

norepi_dur <- function (x, val_var, stop_var, env, ...)
{

  shift <- as.double(hours(1L), units = units(interval(x)))
  x[, c(stop_var) := data.table::fifelse(is.na(get(stop_var)),
                                         get(index_var(x)) + shift,
                                         get(stop_var))]
  x <- unique(x)

  x[, c(val_var) := get(stop_var) - get(index_var(x))]

  x

}

multiply_wbc <- function(..., interval = NULL) {

  cnc <- c("neut", "wbc")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)

  dat[, neut_total := neut*wbc/100]
  dat[, c("wbc", "neut") := NULL]
  dat
}

sofa_resp2 <- function (..., interval = NULL)
{
  score_calc <- function(x) {
    data.table::fifelse(is_true(x < 100), 4L,
      data.table::fifelse(is_true(x < 200), 3L,
        data.table::fifelse(is_true(x < 300), 2L,
          data.table::fifelse(is_true(x < 400), 1L, 0L))))
  }
  cnc <- c("pafi", "vent_ind2")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)
  dat <- dat[is_true(get("pafi") < 200) & !is_true(get("vent_ind2")),
             `:=`(c("pafi"), 200)]
  dat <- dat[, `:=`(c("sofa_resp2"), score_calc(get("pafi")))]
  dat <- rm_cols(dat, cnc, by_ref = TRUE)
  dat
}

sofa_score2 <- function (..., worst_val_fun = max_or_na, explicit_wins = FALSE,
          win_length = hours(24L), keep_components = FALSE, interval = NULL)
{

  cnc <- c("sofa_resp2", "sofa_coag", "sofa_liver", "sofa_cardio",
           "sofa_cns", "sofa_renal")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)
  expr <- substitute(lapply(.SD, fun), list(fun = worst_val_fun))
  if (isFALSE(explicit_wins)) {
    res <- fill_gaps(dat)
    res <- slide(res, !!expr, before = win_length, full_window = FALSE,
                 .SDcols = cnc)
  }
  else {
    if (is_id_tbl(explicit_wins)) {

      res <- hop(dat, !!expr, windows = explicit_wins, .SDcols = cnc)
      if("max_time" %in% names(res))
        res <- as_ts_tbl(res, index_var = "max_time")

    } else {
      res <- slide_index(dat, !!expr, explicit_wins, before = win_length,
                         full_window = FALSE, .SDcols = cnc)
    }
  }
  res <- res[, `:=`(c("sofa2"), rowSums(.SD, na.rm = TRUE)),
             .SDcols = cnc]
  if (isTRUE(keep_components)) {
    res <- rename_cols(res, paste0(cnc, "_comp"), cnc, by_ref = TRUE)
  }
  else {
    res <- rm_cols(res, cnc, by_ref = TRUE)
  }
  res
}
