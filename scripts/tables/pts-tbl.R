library(ricu)
library(jsonlite)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(plyr)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- c("covid19")
cohorts <- ifelse(src == "covid19", 
                  list(config("8-of-10")[["covid19"]]), 
                  list(mimic_demo$icustays$icustay_id))
names(cohorts) <- src

vars <- list(
  age = list(
    concept = "age",
    callback = med_iqr
  ),
  death = list(
    concept = "death",
    callback = percent_fun
  ),
  los_icu = list(
    concept = "los_icu",
    callback = med_iqr
  ),
  los_hosp = list(
    concept = "los_hosp",
    callback = med_iqr
  ),
  gender = list(
    concept = "sex",
    callback = tab_design
  ),
  is_vent = list(
    concept = "is_vent",
    callback = count_percent
  ),
  AKI = list(
    concept = "AKI",
    callback = count_percent
  ),
  RI = list(
    concept = "RI",
    callback = count_percent
  ),
  DM = list(
    concept = "DM",
    callback = count_percent
  ),
  CD = list(
    concept = "CD",
    callback = count_percent
  ),
  CRI = list(
    concept = "CRI",
    callback = count_percent
  ),
  COPD = list(
    concept = "COPD",
    callback = count_percent
  ) #,
  # sofa = list(
  #   concept = "sofa",
  #   callback = multi_med_iqr
  # )
)

pts_source_sum <- function(source, patient_ids) {

  tbl_list <- list()
  for (var in names(vars)) {
    x <- vars[[var]]
    
    skip_vars <- c("DM", "COPD", "CRI", "AKI", "CD", "RI", "is_vent", "vent_total")
    if (var %in% skip_vars & source %in% c("mimic_demo")) next

    data <- load_concepts(x[["concept"]], source, patient_ids = patient_ids, 
                          keep_components = T)

    if (var == "sofa") dat <- data.table::copy(data)

    tbl_list[[var]] <- x[["callback"]](data, patient_ids)
  }

  pts_tbl <- Reduce(rbind,
    lapply(
      tbl_list,
      function(x) data.frame(Reduce(cbind, x))
    )
  )

  cohort_info <- as.data.frame(cbind("Cohort size", "n", length(patient_ids)))

  names(cohort_info) <- names(pts_tbl)

  pts_tbl <- rbind(
    cohort_info,
    pts_tbl
  )

  names(pts_tbl) <- c("Variable", "Reported", srcwrap(source))

  cfg <- get_config("features2", config_dir())
  vars <- load_data(source, cfg, hours(-24L), hours(24L), cohort = patient_ids, 
            cts = T)[[1L]]
  
  var_part <- Reduce(rbind,
         lapply(c("pafi", "bun", "crp", "gcs", "resp", "pco2", "temp", "plt",
                  "ph"), 
                function(cnc) {
         Reduce(cbind, med_iqr(vars[, c(id_vars(vars), cnc), with=F], 1:10))
           
        })
  )
  var_part <- as.data.frame(var_part)
  names(var_part) <- names(pts_tbl)
  
  pts_tbl <- rbind(pts_tbl, var_part)
  
  pts_tbl$Variable <- mapvalues(pts_tbl$Variable,
    from = names(concept_translator),
    to = sapply(names(concept_translator), function(x) concept_translator[[x]])
  )

  pts_tbl

}

res <- Reduce(
  function(x, y) merge(x, y, by = c("Variable", "Reported"), sort = F, all = T),
  Map(pts_source_sum, src, cohorts)
)

my_doc <- read_docx()

my_doc <- my_doc %>%
  body_add_table(res, style = "table_template")

print(my_doc, target = "~/Table1.docx")
