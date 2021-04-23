library(ricu)
library(ggplot2)
library(precrec)
library(assertthat)
library(data.table)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- "covid19"
type <- "conclusive"

dat <- set_design(type = type, output = "cts", impute_vals = FALSE, 
                  pids = config("conclusive")[["covid19"]],
                  add_cmp = FALSE)[[1]]

res <- miss_summ(dat, feats = c("bun", "crp", "gcs", "o2sat", "resp", "pafi", 
                                "safi", "crea", "gcs", "temp", "ph", "pco2"))

save(res, file = file.path("~", paste0("miss_", Sys.Date(), ".RData")))

### missingness intersection
MSI <- function(dat, cnc, type = "intersect") {
  
  cpc <- lapply(cnc, function(cc) id_col(dat[is.na(get(cc))]))
  cat("Total number appearing", length(unique(Reduce(union, cpc))), "\n")
  cat("Size of intersection", length(unique(Reduce(intersect, cpc))), "\n")
  
  if (type == "union") return(unique(Reduce(union, cpc)))
  unique(Reduce(intersect, cpc))
}

rm_cohort <- MSI(dat, c("bun", "resp", "crp", "pafi",  "pco2", "ph"), 
                 type = "union")

length()

### generate the cohort

cohort <- list(
  covid19 = intersect(config("conclusive")[["covid19"]], 
                      setdiff(id_col(dat), rm_cohort)),
  mimic_demo = mimic_demo$icustays$icustay_id
)

#config("mechvent-clean", cohort)
#config("icu24-clean", config("cohort-clean"))
