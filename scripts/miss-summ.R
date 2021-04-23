
miss_summ <- function(dat, cohort = config("conclusive")[["covid19"]], feats) {
  
  w_dat <- merge(dat, load_concepts(c("wave", "hosp"), src), all.x = T)
  
  df <- data.table(
    id = id_col(dat),
    count = rowSums(!is.na(dat[, -c(id_vars(dat), "death"), with = FALSE])),
    wave = w_dat[["wave"]], hospital = w_dat[["hosp"]]
  )
  
  # missingness by hospital
  hosp_view <- ggplot(df[, list(No_meas = mean(count), size = .N), 
                         by = "hospital"], aes(y = No_meas, x = hospital)) +
                  geom_col() + theme_minimal()
  
  # missingness by feature
  if (!is.null(cohort)) w_dat <- w_dat[get(id_vars(w_dat)) %in% cohort]
  cat("Looking at ", nrow(dat), "patients\n")
  feat_proj <- sapply(w_dat, function(x) mean(is.na(x)))
  
  ids <- which(names(feat_proj) %in% feats)
  feat_proj <- feat_proj[ids]
  
  feat_df <- data.frame(names(feat_proj), feat_proj)
  names(feat_df) <- c("Feature", "Proportion")
  
  feat_view <- ggplot(feat_df, aes(y = Proportion, x = Feature)) +
    geom_col() + theme_minimal() + geom_hline(yintercept = 0.33, color = "red")
  
  
  list(hosp = hosp_view, feat = feat_view)
  
}