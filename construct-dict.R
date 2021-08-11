library(ricu)

ricu_dict <- get_config("concept-dict", ricu:::default_config_path())
# unit_list <- subset(covid19$single_timestamp,
#   !duplicated(data.table::setDT(list(pacmed_subname, unit_name))),
#   select = c("pacmed_subname", "unit_name"), part_safe = T)
#save(unit_list, file = "/home/drago/unit_list.RData")
load("/home/drago/unit_list.RData")

cncpt <- c("alkaline", "^basophils_", "^bicar", "^bilirubin_[tu]",
  "^calcium_[bt]", "^chloride", "^creatinine_[bu]", "bp_diastolic", "^eosinophils_",
  "^fio2", "heart_rate_[mu]", "hemato", "^hemoglobin_", "_inr", "^lactate_[abmuw]",
  "^lymphocytes_", "^magnes", "bp_mean", "^neutrophils_",
  "^so2_", "^pco2_[au]", "po2_arterial|po2_unspecified$", "^ph_[abu]",
  "^ph_[cmv]", "^phosphate_", "^thrombocytes_[bu]",
  "^potass", "^prothromb", "^erythrocytes", "^respiratory_rate_",
  "^sodium", "bp_systolic", "^ureum_[bu]", "^leuko", "temp", "^albumin_",
  "^creatine_kin", "^end_tidal", "^crp_", "^gcs_eye", "^gcs_motor", "^gcs_verbal",
  "^gcs_total", "^d_dimer", "^ventilatory_ratio", "^lung_compliance",
  "^pco2_minus_end_tidal_co2", "^aptt_", "^fibrinogen_", "^procalcitonin_",
  "^lactate_dehydrogenase_", "^adjusted_sofa_total$")

dict_names <- c("alp", "basos", "bicar", "bili", "ca", "cl", "crea", "dbp",
  "eos", "fio2", "hr", "hct", "hgb", "inr_pt", "lact", "lymph",
  "mg", "map", "neut", "o2sat", "pco2",
  "po2", "ph", "ph_ven", "phos", "plt", "k", "pt", "rbc", "resp",
  "na", "sbp", "bun", "wbc", "temp", "alb", "ck", "etco2", "crp",
  "egcs", "mgcs", "vgcs", "tgcs", "d_dim", "vent_ratio", "lung_comp",
  "dead_space", "ptt", "fgn", "prcl", "ldh", "sofa_total")

L <- lapply(cncpt, function(cnc)
  grep(cnc, covid19$parameters$pacmed_subname, value = T, ignore.case = T))

names(L) <- dict_names

make_item <- function(c_names) {
  
  target_table <- unique(
    subset(covid19$parameters, pacmed_subname %in% c_names)[["table_source"]]
  )
  
  if (length(target_table) != 1L) browser()
  assertthat::assert_that(length(target_table) == 1L)
  
  if (!grepl("eosinophils|lymphocytes|basophils|eosinophils", c_names[[1L]])) {
    
    itm <- list(
      ids = c_names,
      table = target_table,
      sub_var = "pacmed_subname", class = "sel_itm"
    )

    
    if (grepl("so2_", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 100))"
    if (grepl("^fibrinogen_", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 100))"
    if (grepl("^bilirubin_", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 0.058467))"
    if (grepl("^ureum_[bu]", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 2.8))"
    if (grepl("^magnes", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 2.431))"
    if (grepl("^phosphate_", c_names[[1L]])) itm[["callback"]] <- 
        "transform_fun(binary_op(`*`, 3.097521))"
    if (grepl("^creatinine_[bu]", c_names[[1L]])) itm[["callback"]] <-
        "transform_fun(binary_op(`*`, 0.011312))"
    if (grepl("^albumin_", c_names[[1L]])) itm[["callback"]] <-
        "transform_fun(binary_op(`*`, 0.1))"
    if (grepl("^calcium_[bt]", c_names[[1L]])) itm[["callback"]] <-
        "transform_fun(binary_op(`*`, 4.008))"
    if (grepl("^hemoglobin_", c_names[[1L]])) itm[["callback"]] <-
        "transform_fun(binary_op(`*`, 1.611344))"
    
    itm <- list(itm)
    
  } else {
    
    perc_item <- grep("_perc_", c_names, value = T)
    itm <- list(
      list(
        ids = perc_item,
        table = target_table,
        sub_var = "pacmed_subname", class = "sel_itm"
      ),
      list(
        ids = setdiff(c_names, perc_item),
        table = target_table, callback = "blood_cell_callback",
        sub_var = "pacmed_subname", class = "sel_itm"
      )
    )
    
  }
  
  itm

}

dict <- Map(function(x, y) {
  
  if (length(y) == 0L) browser()
  res <- list(
    sources = list(
      covid19 = make_item(y)
    )
  )

  # check units
  covid19_units <- unique(unit_list$unit_name[unit_list$pacmed_subname %in% y])
  covid19_units <- covid19_units[!is.na(covid19_units)]
  target <- ricu_dict[[x]][["unit"]]

  solved <- c(
    "o2sat", "bun", "bili", "fgn", "lact", "crp", "alb",
    "plt", "na", "k", "rbc", "mg", "wbc", "phos", "crea",
    "ca", "hgb", "basos", "eos", "lymph", "neut"
  )

  if(!all(covid19_units %in% target) & !(x %in% solved)) {

    cat("Possible unit mismatch for", x, ":",
      "target is", paste(target, collapse = ","),
      ", but current:", paste(covid19_units, collapse = ","), "\n")

  }

  res

}, dict_names, L)


dict[["death"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = c("gender"),
           class = "col_itm", callback = "covid19_death",
           index_var = "death_timestamp")
    )
  )
)

dict[["conclusive"]] <- list(
  target = "id_tbl",
  sources = list(
    covid19 = list(
      list(table = "outcomes", val_var = c("still_admitted"),
           extra_var = c("transfer"),
           class = "col_itm", callback = "covid19_conclusive")
    )
  )
)


dict[["age"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = "age", class = "col_itm")
    )
  )
)

dict[["nlr"]] <- list(
  concepts = c("neut", "lymph"), callback = "nlr_callback", class = "rec_cncpt"
)

dict[["vent_ind2"]] <- list(
  class = "lgl_cncpt",
  sources = list(
    covid19 = list(
      list(table = "intubations", val_var = "successful_extubation",
        end_var = "end_intubation", callback = "span_window", class = "col_itm")
    )
  )
)

dict[["los_bef_icu"]] <- list(
  target = "id_tbl",
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = c("hospital_admission_timestamp"),
           class = "col_itm")
    )
  )
)

dict[["los_hosp"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = c("hospital_admission_timestamp"),
        end_var = "hospital_discharge_timestamp",
        class = "col_itm",
        callback = "los_hosp")
    )
  )
)

dict[["hosp"]] <- list(
  target = "id_tbl",
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = c("hospital_first"),
        class = "col_itm")
    )
  )
)

dict[["wave"]] <- list(
  target = "id_tbl",
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = c("gender"),
        class = "col_itm", callback = "wave_callback")
    )
  )
)


dict[["n_comb"]] <- list(
  target = "id_tbl",
  description = "No. comorbidities",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", 
           val_var = "chronic_dialysis", cmb_2 = "chronic_renal_insufficiency", 
           cmb_3 = "cirrhosis", cmb_4 = "copd", cmb_5 = "diabetes", 
           cmb_6 = "hematologic_malignancy", cmb_7 = "immunodeficiency", 
           cmb_8 = "neoplasm", cmb_9 = "respiratory_insufficiency", 
           cmb_10 = "cardiovascular_insufficiency", cmb_11 = "acute_kidney_injury",
           class = "col_itm",
           callback = "cmbd_callback"
      )
    )
  )
)


dict[["AKI"]] <- list(
  target = "id_tbl",
  description = "Acute kidney injury",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("acute_kidney_injury"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["DM"]] <- list(
  target = "id_tbl",
  description = "Diabetes Mellitus",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("diabetes"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["COPD"]] <- list(
  target = "id_tbl",
  description = "Chronic obstructive pulmonary disease",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("copd"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["CD"]] <- list(
  target = "id_tbl",
  description = "Chronic dialysis",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("chronic_dialysis"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["CRI"]] <- list(
  target = "id_tbl",
  description = "Chronic renal insufficiency",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("chronic_renal_insufficiency"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["RI"]] <- list(
  target = "id_tbl",
  description = "Chronic respiratory insufficiency",
  sources = list(
    covid19 = list(
      list(table = "comorbidities", val_var = c("respiratory_insufficiency"),
        class = "col_itm", callback = "nice_callback"
      )
    )
  )
)

dict[["sex"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "episodes", val_var = "gender", class = "col_itm", 
           callback = "apply_map(c(M = 'Male', V = 'Female'))")
    )
  )
)

dict[["los_icu"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "admissions", val_var = c("admission_timestamp"),
        end_var = "discharge_timestamp",
        class = "col_itm",
        callback = "los_icu")
    )
  )
)

dict[["norepi_rate"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "medications", ids = c("med_norepinephrine"),
           stop_var = "end_timestamp", sub_var = "pacmed_subname", 
           class = "sel_itm", callback = "norepi_rate_kg")
    )
  )
)

dict[["norepi_dur"]] <- list(
  sources = list(
    covid19 = list(
      list(table = "medications", ids = c("med_norepinephrine"),
           stop_var = "end_timestamp", sub_var = "pacmed_subname", 
           class = "sel_itm", callback = "norepi_dur")
    )
  )
)

dict[["vent_total"]] <- list(
  concepts = c("vent_ind2"), callback = "vent_total_callback", 
  class = "rec_cncpt", target = "id_tbl"
)

dict[["is_vent"]] <- list(
  concepts = c("vent_ind2"), callback = "is_vent_callback", class = "rec_cncpt",
  target = "id_tbl"
)

dict[["grs"]] <- list(
  description = "Galloway Risk Score",
  concepts = c("age", "sex", "o2sat", "neut_total", "crp", "alb", "crea", "DM",
               "COPD"), 
  callback = "grs_callback", class = "rec_cncpt"
)

dict[["four_c"]] <- list(
  description = "ISARIC 4C Mortality score",
  concepts = c("age", "sex", "n_comb", "resp", "o2sat", "gcs", "bun", "crp"), 
  callback = "four_c_score_callback", class = "rec_cncpt"
)

dict[["saps_3"]] <- list(
  description = "SAPS-III Mortality score",
  concepts = c("age", "los_bef_icu", "norepi_rate",  
               "gcs", "bili", "temp", "crea", "hr", "wbc", "ph", "plt", "sbp",
               "pafi", "po2", "vent_ind2"), 
  callback = "saps_3_callback", class = "rec_cncpt"
)

dict[["neut_total"]] <- list(
  description = "Neutrophils total count",
  concepts = c("neut", "wbc"), 
  callback = "multiply_wbc", class = "rec_cncpt"
)

dict[["sofa_resp2"]] <- list(
  concepts = c("pafi", "vent_ind2"),
  description = "SOFA respiratory component",
  category = "outcome",
  callback = "sofa_resp2",
  class = "rec_cncpt"
)

dict[["sofa2"]] <- list(
  concepts = c("sofa_resp2", "sofa_coag", "sofa_liver", "sofa_cardio",
               "sofa_cns", "sofa_renal"),
  description = "sequential organ failure assessment score 2",
  category = "outcome",
  callback = "sofa_score2",
  class = "rec_cncpt"
)

jsonlite::write_json(dict, file.path("~/config/concept-dict.json"), pretty = T)
