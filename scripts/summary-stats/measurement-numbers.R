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
cohort <- config("8-of-10")[[src]]
cfg <- config("features3")

load("~/dat1.RData")

res <- dat[get(id_var(dat)) %in% cohort &
           get(index_var(dat)) > hours(-24L) &
           get(index_var(dat)) <= hours(24L), setdiff(names(cfg), "age"), 
           with = FALSE]

cat("We evaluated a total of",
    paste(
      colSums(!is.na(res)),
      vapply(names(res), function(x) cfg[[x]][["full_name"]], ""),
      collapse = ", "
    ), 
    "measurements collected in the 48 hour period prior to 24 hours into ICU admission")




