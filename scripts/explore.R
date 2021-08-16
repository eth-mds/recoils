library(ricu)
library(ggplot2)
library(precrec)
library(ranger)
library(magrittr)
library(officer)
library(assertthat)
library(data.table)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- "covid19"
type <- "config"
res_type <- "tabular"
split_waves <- FALSE
mech_vent <- FALSE
# cohort <- get_cohort(src, type = type,
#                      hospitals = c("stantonius", "martini", "maasstad",
#                                    "geldersevallei", "etz", "umcu", "vumc",
#                                    "franciscus", "albertschweitzer", "amphia",
#                                    "laurentius"))
cohort <- config("conclusive")[[src]]

if (src == "covid19") {
  cfg <- get_config("features2", config_dir())
  vent_cnc <- "vent_ind2"
} else {
  cfg <- get_config("features", config_dir())
  vent_cnc <- "vent_ind"
}

if (mech_vent) {
  
  upr <- load_concepts(vent_cnc, src, patient_ids = cohort)
  upr <- upr[, head(.SD, n = 1L), by = c(id_vars(upr))]
  upr[, max_time := get(index_var(upr))]
  upr <- upr[, c(id_vars(upr), "max_time"), with=F]
  
  lwr <- data.table::copy(upr)
  lwr[, min_time := max_time - hours(24L)]
  lwr[, max_time := NULL]
  
  cohort <- id_col(upr)
  lwr <- list(lwr)
  upr <- list(upr)
  fig_ext <- "_MV"
  
  
} else {
  
  upr <- hours(24L)
  lwr <- upr - hours(48L)
  cohort <- config("conclusive")[[src]]
  fig_ext <- ""
  
}

dat <- load_data(src, cfg, lwr, upr, impute_vals = T, cts = T, 
                 cohort = cohort)[[1L]]
plt <- component_plots(dat, src, split_waves, res_type = res_type)

sys_select <- c("age", "metabolic", "coag", "resp", "inflam", "bone_marrow", 
                "liver", "renal")
final_plots <- lapply(sys_select, function(s) plt[[src]][[s]])
names(final_plots) <- sys_select

if (res_type != "tabular") {
  
  if (!file.exists("res")) dir.create("res")
  if (!file.exists(file.path("res", type))) dir.create(file.path("res", type))
  for (sys in sys_select) {
    
    ggsave(file.path(".", "res", type,
                     paste0(sys, fig_ext, ".png")), plot = final_plots[[sys]], 
           height = 8, width = 8+8*split_waves)
    
  }
  
  
} else {
  
  tbl <- Reduce(rbind, final_plots)
  
  RF <- ranger(death ~ . - episode_id, dat, importance = "impurity")
  imp <- importance(RF)
  feat_names <- sapply(names(imp), function(x) {
    if(!is.null(cfg[[x]][["full_name"]])) return(cfg[[x]][["full_name"]])
    return(x)
  })
  rf_imp <- data.table::data.table(importance = imp, feature = feat_names)
  
  obj <- merge(tbl, rf_imp, by = "feature", all = T)
  
  # add range, increase/decrease and imputation value
  rng <- paste0(
    sapply(cfg, `[[`, "lower"), "-",
    sapply(cfg, `[[`, "upper")
  )
  inc <- paste(
    (-1)^(sapply(cfg, `[[`, "direction") == "decreasing") * 
      sapply(cfg, `[[`, "step")
  )
  ridi_names <- sapply(names(cfg), function(x) {
    if(!is.null(cfg[[x]][["full_name"]])) return(cfg[[x]][["full_name"]])
    return(x)
  })
  imp <- sapply(cfg, `[[`, "impute_val")
  unt <- sapply(cfg, function(x) x[["unit"]][[1L]])
  ridi <- data.table::data.table(feature = ridi_names, rng = rng, inc = inc,
                                 imp = imp, unt = unt)
  
  obj <- merge(obj, ridi, by = "feature")
  
  obj[["AUROC"]] <- round(obj[["AUROC"]], 3)
  obj[["importance"]] <- round(obj[["importance"]], 2)
  obj[["feature"]] <- paste0(
    stringr::str_to_sentence(obj[["feature"]]), " (",
    obj[["unt"]], ")"
  )
  setorderv(obj, c("AUROC", "importance"), -1L)
  
  my_doc <- read_docx()
  
  my_doc <- my_doc %>%
    body_add_table(obj, style = "table_template")
  
  print(my_doc, target = "~/Supplementary_Table2.docx")
  
}

# cohort <- list(
#   covid19 = id_col(load_concepts("conclusive", "covid19")[conclusive == 1L]),
#   mimic_demo = unique(mimic_demo$icustays$icustay_id)
# )
# 
# config("conclusive")
