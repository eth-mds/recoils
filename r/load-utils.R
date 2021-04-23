
load_difftime.covid19_tbl <- function(x, rows, cols = colnames(x),
  id_hint = id_vars(x),
  time_vars = ricu::time_vars(x), ...) {
  
  ricu:::warn_dots(...)
  
  ricu:::load_mihi(x, {{ rows }}, cols, id_hint, time_vars)
}

id_win_helper.covid19_env <- function(x) {
  
  merge_inter <- function(x, y) {
    merge(x, y, by = intersect(colnames(x), colnames(y)))
  }
  
  get_id_tbl <- function(tbl, id, start, end, aux) {
    as_src_tbl(x, tbl)[, c(id, start, end, aux)]
  }
  
  cfg <- sort(as_id_cfg(x), decreasing = TRUE)
  
  ids  <- vctrs::field(cfg, "id")
  sta <- vctrs::field(cfg, "start")
  end  <- vctrs::field(cfg, "end")
  
  res <- Map(get_id_tbl, vctrs::field(cfg, "table"), ids, sta,
    end, c(as.list(ids[-1L]), list(NULL)))
  res <- Reduce(merge_inter, res)
  
  res <- res[, c(sta, end) := lapply(.SD, ricu:::as_dt_min, get(sta[1L])),
    .SDcols = c(sta, end)]
  
  ricu:::order_rename(res, ids, sta, end)
}
