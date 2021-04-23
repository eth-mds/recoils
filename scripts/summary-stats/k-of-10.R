library(ricu)
library(ggplot2)
library(precrec)
library(assertthat)
library(data.table)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- "covid19"
regen <- TRUE
type <- "conclusive"

# generate the adult conclusive cohort
if (regen) {
  
  config(
    "conclusive",
    list(
      covid19 = id_col(load_concepts(c("conclusive", "age"), 
                                     "covid19")[conclusive == 1 & 
                                                (age >= 18L | is.na(age))]
                       ),
      mimic_demo = mimic_demo$icustays$icustay_id
    )
  )
  
}

dat <- set_design(type = type, output = "cts", impute_vals = FALSE, 
                  pids = config("conclusive")[["covid19"]],
                  cfg = get_config("features3", config_dir()),
                  add_cmp = FALSE)[[1]]

k_miss <- rowSums(is.na(dat))

dat[, k_miss := k_miss]

k_thresh <- 2

coh <- id_col(dat[k_miss <= k_thresh])

overlap <- mean(config("icu24-clean")[["covid19"]] %in% coh)

k_of_10 <- list(
  covid19 = coh,
  mimic_demo = mimic_demo$icustays$icustay_id
)

config("8-of-10", value = k_of_10)
