library(ricu)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

type <- "icu24"

if (type == "icu24") {
  
  tbl <- load_concepts("sofa2", "covid19", explicit_wins = hours(c(24L, 48L)))
  
  sofa_baseline <- list(
    tbl[get(index_var(tbl)) == hours(24L)],
    tbl[get(index_var(tbl)) == hours(48L)]
  )
  
  save(sofa_baseline, file = "~/sofa_icu24.RData")
  
  tbl <- load_concepts("four_c", "covid19", explicit_wins = hours(c(24L, 48L)))
  
  four_c_baseline <- list(
    tbl[get(index_var(tbl)) == hours(24L)],
    tbl[get(index_var(tbl)) == hours(48L)]
  )
  
  save(four_c_baseline, file = "~/four_c_icu24.RData")
  
  saps <- load_concepts("saps_3", "covid19", explicit_wins = hours(c(24L, 48L)),
                        keep_components = TRUE)
  saps_baseline <- list(
    saps[get(index_var(saps)) == hours(24L)],
    saps[get(index_var(saps)) == hours(48L)]
  )
  
  save(saps_baseline, file = "~/saps_3_icu24.RData")
  
} else {
  
  upr <- load_concepts("vent_ind2", "covid19")
  upr <- upr[, head(.SD, n = 1L), by = c(id_vars(upr))]
  upr[, max_time := get(index_var(upr))]
  upr <- upr[, c(id_vars(upr), "max_time"), with=F]
  
  lwr <- data.table::copy(upr)
  lwr[, min_time := max_time - hours(48L)]
  lwr[, max_time := NULL]
  
  lwr <- list(lwr)
  upr <- list(upr)
  
  tbl <- load_concepts("sofa2", "covid19",
                       explicit_wins = merge(lwr[[1L]], upr[[1L]]))
  
  sofa_baseline <- list(tbl)
  
  save(sofa_baseline, file = "~/sofa_mechvent.RData")
  
  tbl <- load_concepts(c("four_c"), "covid19",
                explicit_wins = merge(lwr[[1L]], upr[[1L]]))
  
  four_c_baseline <- list(tbl)
  
  save(four_c_baseline, file = "~/four_c_mechvent.RData")
  
}


