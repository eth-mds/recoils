library(ricu)
library(ggplot2)
library(precrec)
library(assertthat)
library(data.table)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

cfg <- get_config("features2", config_dir())
set.seed(2020)

spital <- load_concepts("hosp", "covid19")[, list(total = .N), by = "hosp"]
res <- NULL

master_tbl <- load_concepts(c(names(cfg), "hosp"), "covid19")
master_tbl <- master_tbl[!is.na(get(index_var(master_tbl)))]

for (feat in names(cfg)) {
  
  tbl <- master_tbl[, c(meta_vars(master_tbl), feat, "hosp"), with = F]
  
  if (feat == "age") {
    tbl[, c(index_var(tbl)) := hours(0L)]
    tbl <- unique(tbl)
  }
  
  tbl <- tbl[!is.na(get(feat))]
  
  for (mark in hours(c(24, 48, Inf))) {
    
    bot <- merge(
      spital,
      tbl[get(index_var(tbl)) < mark, 
          list(count = length(unique(get(id_var(tbl))))), by = "hosp"], 
      all.x = T
    )
    
    bot <- cbind(bot, feature = feat, mark = mark)
    res <- rbind(res, bot)
    
  }
  
}

res[is.na(count), "count"] <- 0
res[, proportion := 100*(count/total)]

plots <- lapply(hours(c(24, 48, Inf)), function(upr_win) {
  
  ggplot(res[mark == upr_win & feature != "fgn"], aes(x = feature, y = hosp, fill = proportion)) +
    geom_tile() + theme_minimal() + scale_fill_viridis_c() +
    geom_text(aes(label = round(proportion))) +
    ggtitle(paste("Feature missingnes by hospital at", upr_win, "hours"))
  
})

save(res, file = "missingness_26Feb.RData")

