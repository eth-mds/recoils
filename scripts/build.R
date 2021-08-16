library(ricu)
library(ggplot2)
library(assertthat)
library(precrec)
library(cowplot)
library(ranger)
library(magrittr)
library(officer)
library(glmnet)
library(data.table)

r_dir <- file.path(rprojroot::find_root(".git/index"), "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

#set.seed(7)
#f_name <- paste0("build", fig_ext, "_", Sys.Date(), ".RData")

# Choose the type
type <- "icu24"
feature_set <- "features3"
dev_val <- TRUE

# Get the design
design <- set_design(type = type, output = ifelse(dev_val, "mat", "cts"),
                     impute_vals = FALSE, 
                     cfg = get_config(feature_set, config_dir()))
save(design, file = paste0("~/design_", type, ".RData"))

# Run locally from here
cat("scp -r drago@10.60.0.18:~/design_", type, ".RData  ~/covid-19/naps-res/",
    sep = "")
load(paste0("~/covid-19/naps-res/design_", type, ".RData"))

# Train the model
if (dev_val) {

  X <- design[["X"]]
  y <- design[["y"]]
  ids <- design[["ids"]]
  des <- design[["design"]]

  dev_info <- scaling_factor(X, y, ids, by_hosp = TRUE, seed = 0.7) # 24L
  
  score <- dev_info[["score"]] 
  score[score < 0] <- 0
  attr(score, "concept") <- gsub("\\..*", "", names(score))
  attr(score, "threshold") <- 
    Reduce(c, lapply(config(feature_set),
                     function(x) {
                       
                       sq <- seq(x[["lower"]], x[["upper"]], x[["step"]])
                       if (x[["direction"]] == "increasing") 
                         sub_idx <- 1:(length(sq)-1) else sub_idx <- -1
                       
                         sq[sub_idx]
                       
                     }
                     ))
  attr(score, "right") <- 
    Reduce(c, lapply(config(feature_set),
                     function(x) {
                       sq <- seq(x[["lower"]], x[["upper"]], x[["step"]])
                       rep(x[["direction"]] == "increasing", length(sq) - 1L)
                     }))
  
  cfg <- config(feature_set)
  score2table(score, cfg)
  
  naps_s <- as.matrix(des[[1L]][, names(score),with=FALSE]) %*% score

  des[[1L]] <- des[[1L]][, naps := naps_s]

}

# Evaluate the model
models <- c("four_c", "sofa2", "saps_3", "age", "death")

if (dev_val) {

  models <- c(models, "naps")
  set.seed(123)
  eval_dat <- list(
    des[[1L]][get(id_vars(des[[1L]])) %in% dev_info[["dev"]], models, with=FALSE],
    des[[1L]][get(id_vars(des[[1L]])) %in% dev_info[["val"]], models, with=FALSE]
  )

} else {
  
  des <- design[["design"]]
  eval_dat <- des[[1L]][, models, with = FALSE]

}

deval <- eval_at(eval_dat, type = type)
add_titles(deval)
ggsave(filename = file.path(r_dir, "..", "paper", "Figure1.tiff"),
       width = 10, height = 10, bg = "white")

### do the calibration plot
{
  res <- ids
  res <- res[, naps := naps_s]
  res[, death := y]
  
  res[, mort := mean(death), by = "naps"]
  res[, bin := .bincode(naps, c(-Inf, seq(2.1, 13.1), Inf))]
  res[, bin_mort := mean(death), by = "bin"]
  
  brier <- mean((res$mort - res$death)^2)
  cat("Brier score is", brier, "\n")

  expit <- function(x) exp(x) / (1 + exp(x))
  res[, fit := expit(-2.886 + naps/4)]
  res[, development := get(id_var(res)) %in% dev_info[["dev"]]]
  
  nbt <- 100
  set.seed(123)
  bplot <- NULL
  bdat <- res[, c("naps", "death"), with=FALSE]
  calib_fit(bdat) + ggtitle("Mortality rate per score level")
  ggsave(filename = file.path(r_dir, "..", "paper", "Figure2.tiff"),
         width = 10, height = 6.5)
  
  cb_fit <- calib_fit(bdat, plot = FALSE)
  cb <- read_docx()
  cb <- cb %>%
    body_add_table(cb_fit, style = "table_template")
  print(cb, target = file.path(r_dir, "..", "paper", 
                                   "Supplementary_Table4.docx"))
}

### Beta coefficient
{
  
  beta <- score
  beta[seq_along(beta)] <- 0
  beta[seq_along(beta)] <- dev_info[["beta"]]
  
  coef <- beta2table(beta, cfg)
  
  my_doc <- read_docx()
  
  my_doc <- my_doc %>%
    body_add_table(coef, style = "table_template")
  
  print(my_doc, target = file.path(r_dir, "..", "paper", 
                                   "Supplementary_Table3.docx"))
  
}

### Stratify on wave and sex

# sex
# config("female",
#        value = list(id_col(load_concepts("sex", "covid19")[sex == "Female"])))
fem <- config("female")[[1L]]
sx_split <- list(
  des[[1L]][get(id_vars(des[[1L]])) %in% fem, c("naps", "death"), with=FALSE],
  des[[1L]][!(get(id_vars(des[[1L]])) %in% fem), c("naps", "death"), with=FALSE]
)

p_sex <- eval_strat(sx_split, lab_1 = "Female", lab_2 = "Male")

# wave
# config("first-wave", 
#        value = list(id_col(load_concepts("wave", "covid19")[wave == "1st"])))
first <- config("first-wave")[[1L]]
wav_split <- list(
  des[[1L]][get(id_vars(des[[1L]])) %in% first, c("naps", "death"), with=FALSE],
  des[[1L]][!(get(id_vars(des[[1L]])) %in% first), c("naps", "death"), with=FALSE]
)

p_wav <- eval_strat(wav_split)

aux <- cowplot::plot_grid(p_sex, p_wav, ncol = 2L)
ggsave(filename = file.path(r_dir, "..", "paper", "eFigure1.tiff"),
       width = 10, height = 6)
