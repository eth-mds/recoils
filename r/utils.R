
proj_root <- function() rprojroot::find_root(".git/index")

config_dir <- function() file.path(proj_root(), "config")

aggreg_fun <- function(cfg, inc = "max", dec = "min") {

  dir <- vapply(cfg, `[[`, character(1L), "direction")

  #assert_that(setequal(dir, c("increasing", "decreasing")))

  ifelse(dir == "increasing", inc, dec)
}

si_cohort <- function(source, age_threshold = 18L, ...) {

  if (grepl("eicu", source)) {

    patient_ids <- get_config("eicu-cohort", config_dir())[["eicu"]][["all"]]

    susp_infec <- load_concepts("susp_inf", source, abx_min_count = 2L,
      si_mode = "or", patient_ids = patient_ids, ...)

  } else if (identical(source, "hirid")) {

    susp_infec <- load_concepts("susp_inf", source, abx_min_count = 2L,
                     si_mode = "or", ...)

  } else {

    susp_infec <- load_concepts("susp_inf", source, id_type = "icustay", ...)

  }

  si_lwr <- hours(48L)
  si_upr <- hours(24L)
  susp_infec <- susp_infec[is_true(get("susp_inf")), ]
  susp_infec <- susp_infec[, c("susp_inf") := NULL]
  susp_infec <- susp_infec[, c("si_lwr", "si_upr") := list(
    get(index_var(susp_infec)) - si_lwr,
    get(index_var(susp_infec)) + si_upr
  )]

  susp_infec <- susp_infec[(si_lwr <= 0L) & (si_upr >= 0L)]

  above_age <- load_concepts("age", source)[age > age_threshold]

  unique(intersect(id_col(susp_infec), id_col(above_age)))
}

get_cohort <- function(src, type = "config", hospitals = NULL) {

  if (type == "hospital" & src == "covid19")
    return(id_col(load_concepts("hosp", src)[hosp %in% hospitals]))


  config("cohort")[[src]]

}

strat_samp <- function(dat, by, frac = 0.75, replace = FALSE) {

  samp <- function(x, repl) sample(x, length(x) * frac, repl)

  assert_that(inherits(dat, "data.frame"), has_name(dat, by))

  rows <- split(seq_len(nrow(dat)), dat[, by, with = FALSE])

  unlist(lapply(rows, samp, replace))
}

score2table <- function(score, cfg) {
  flip_sum <- function(x) {
    x <- if (all(x[["Dir"]])) x else x[rev(seq_len(nrow(x))), ]
    x <- x[x[["Points"]] > 0, ]
    x[["Points"]] <- cumsum(x[["Points"]])
    x
  }


  full_names <- sapply(attr(score, "concept"),
    function(x) {
      name <- stringr:::str_to_title(cfg[[x]][["full_name"]])

      name <- paste0(name, " (", cfg[[x]][["unit"]][1L], ")")

      name

    })
  tbl <- data.frame(
    Feature = full_names, Dir = attr(score, "right"),
    Thresh = attr(score, "threshold"), Points = score, stringsAsFactors = FALSE
  )
  
  tbl <- do.call(rbind, lapply(split(tbl, tbl[["Feature"]]), flip_sum))
  tbl[["Thresh"]] <- paste(ifelse(tbl[["Dir"]], "&ge;", "&lt;"),
                           tbl[["Thresh"]])
  res <- kableExtra::kable(tbl[ , c("Thresh", "Points")], row.names = FALSE,
                           col.names = c("", "Points"), escape = FALSE)
  ind <- rle(tbl[ , "Feature"])
  ind <- stats::setNames(ind[["lengths"]], gsub("_", " ", ind[["values"]]))
  res <- kableExtra::pack_rows(res, index = ind)
  res
}

beta2table <- function(score, cfg) {
  flip_order <- function(x) {
    x <- if (all(x[["Dir"]])) x else x[rev(seq_len(nrow(x))), ]
    x <- x[x[["Beta"]] > 0, ]
    x
  }


  full_names <- sapply(attr(score, "concept"),
    function(x) {
      name <- stringr:::str_to_sentence(cfg[[x]][["full_name"]])

      #name <- paste0(name, " (", cfg[[x]][["unit"]][1L], ")")

      name

    })
  tbl <- data.frame(
    Feature = full_names, Dir = attr(score, "right"),
    Thresh = attr(score, "threshold"), Beta = score, stringsAsFactors = FALSE
  )

  tbl[["unit"]] <- vapply(attr(beta, "concept"),
                          function(x) cfg[[x]][["unit"]][[1L]], FUN.VALUE = "")

  tbl <- do.call(rbind, lapply(split(tbl, tbl[["Feature"]]), flip_order))
  tbl[["Thresh"]] <- paste(ifelse(tbl[["Dir"]], "≥", "<"),
                           tbl[["Thresh"]])


  tbl[["Score threshold"]] <- paste0(tbl[["Feature"]], " ",
                                     tbl[["Thresh"]], " ",
                                     tbl[["unit"]])
  tbl[["Beta x4 (rounded)"]] <- round(tbl[["Beta"]] * 4)
  tbl[["Beta"]] <- round(tbl[["Beta"]], 3)
  tbl[, c("Score threshold", "Beta", "Beta x4 (rounded)")]

}

get_cores <- function(n_cores = NULL) {

  max_cores <- as.integer(
    Sys.getenv("LSB_DJOB_NUMPROC", unset = parallel::detectCores())
  )

  if (is.null(n_cores)) {
    n_cores <- max_cores
  }

  if (n_cores > max_cores) {
    warning("asked for ", n_cores, " cores but only ", max_cores,
            " are available")
    n_cores <- max_cores
  }

  n_cores
}

roc_pr_plot <- function(x, ..., nrow_legend = 1) {

  roc <- autoplot(x, "ROC", ...)
  prc <- autoplot(x, "PR", ...)

  legend <- cowplot::get_legend(
    roc + guides(color = guide_legend(nrow = nrow_legend)) +
          theme(legend.position = "bottom")
  )

  plot <- cowplot::plot_grid(roc + theme(legend.position = "none"),
                             prc + theme(legend.position = "none"), ncol = 2)

  cowplot::plot_grid(plot, legend, nrow = 2, rel_heights = c(1, .1))
}

auc <- function(x) data.table::setDT(precrec::auc(x))

is_difftime <- function(x) inherits(x, "difftime")

is_interval <- function(x) {
  assert_that(is_difftime(x), length(x) > 0L) && all(x >= 0)
}

noreq <- function(source, press = "map", beta = 10, patient_ids = NULL,
  upr = hours(24L), lwr = hours(0L)) {

  na_z <- function(x) ifelse(is.na(x), 0, x)
  col_name <- paste0(press, beta)
  imp_val <- cfg[[press]][["upper"]]*1.01

  tbl <- load_concepts(c("norepi_rate", "epi_rate", "dopa_rate", "dobu_rate", press),
    source, patient_ids = patient_ids)
  tbl[, noreq := (na_z(norepi_rate) + na_z(epi_rate) + na_z(dopa_rate)/150 +
                  na_z(dobu_rate > 0)*0.015)]

  tbl[, cf_bp := nafill(get(press), "locf")]

  tbl[, c(col_name) := (cf_bp - beta*noreq)]

  tbl[, c(meta_vars(tbl), col_name), with = FALSE]

  tbl <- tbl[get(index_var(tbl)) <= upr & get(index_var(tbl)) >= lwr,
      -min(get(col_name), na.rm = T), by = eval(id_var(tbl))]

  setnames(tbl, "V1", col_name)
  tbl <- replace_na(tbl, imp_val)

  tbl
}

add_titles <- function(p) {
  
  title1 <- ggdraw() + 
    draw_label(
      "Development",
      fontface = 'bold'
    ) + theme(plot.background = element_rect(color = "black"))
  
  title2 <- ggdraw() + 
    draw_label(
      "Validation",
      fontface = 'bold'
    ) + theme(plot.background = element_rect(color = "black"))
  
  titles <- plot_grid(title1, title2, ncol = 2L)
  
  plot_grid(titles, p, ncol = 1L, rel_heights = c(1, 20))
  
}