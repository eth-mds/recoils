eval_otp <- function(test, dose, data_src) {

  within <- dose_eval(dose, test)
  scor0  <- sofa_eval(data_src, times, cohort = within)


  aucs <- auc(eval_score(scor0))
  aucs <- aucs[, list(mean = mean(aucs), sd = sd(aucs)),
    by = c("modnames", "curvetypes")]
  aucs <- aucs[, c("lwr", "upr") := list(mean - sd, mean + sd)]
  aucs <- aucs[, c("type", "time") := list(sub(" .+", "", modnames),
    as.integer(sub(".+, ", "", sub(" hours]$", "", modnames)))
  )]

  ggplot(aucs, aes(x = time, y = mean, color = type, fill = type)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, linetype = 0) +
    facet_grid(rows = vars(curvetypes), scales = "free_y") +
    theme_bw() +
    xlab("Hours in ICU") + ylab("Mean AUC") +
    theme(legend.position = "bottom") +
    scale_colour_brewer(type = "qual", palette = 6, direction = 1)

}

eval_fxtp <- function(test_24, score) {

  fx_plot <- list()
  fx_roc <- list()

  for (cp in names(score)) {

    # if (length(grep(paste0("^sofa.", cp), names(test_24))) == 0) {
    #   sofa.cp <- rep(0, nrow(test_24))
    #   sofa.name <- "Prev. no"
    # } else {
    #   sofa.cp <- test_24[[grep(paste0("^sofa.", cp), names(test_24))]]
    #   sofa.name <- "SOFA"
    # }

    eval <- evalmod( # again need to beware here
      scores = list(rowSums(test_24[, score[[cp]], with=F])),#, replace_na(sofa.cp, 0)),
      labels = list(test_24[["death"]]), #, test_24[["death"]]),
      #dsids = 1:2, modnames = paste(c("DOSE", sofa.name), cp), cb_alpha = 0.05
      modnames = paste("ECOS", cp)
    )

    fx_roc[[cp]] <-  auc(eval)$aucs[c(T, F)]

    fx_plot[[cp]] <- autoplot(eval, "ROC") + geom_line(size = 1) +
    ggtitle("ROC curve at 24 hours") +
      theme(legend.position = "bottom",
        legend.text = element_text(size=10))

  }

  list(fx_plot = fx_plot, fx_roc = fx_roc)

}

component_plots <- function(dat, data_src, split_waves = FALSE, res_type = "plot") {

  if (split_waves) dat <- merge(dat, load_concepts("wave", data_src), all.x = T)

  direc <- function(x) ifelse(is.null(cfg[[x]]), 1L,
    1 - 2*(cfg[[x]][["direction"]] == "decreasing"))

    sys_perf <- function(dat, split_waves, res_type = "plot") {

      if (split_waves) {

        cat("Hit the splitting \n")

        return(
          cowplot::plot_grid(
            sys_perf(dat[wave == "1st"], FALSE),
            sys_perf(dat[wave == "2nd"], FALSE), ncol = 2L
          )
        )

      }

      scores <- lapply(grid[[i]], function(x) {
        if(is.null(dat[[x]])) return(rep(0, nrow(dat)))
        dat[[x]]*direc(x)
      })
      #scores <- append(scores, list(sofa.cp))
      labels <- lapply(1:length(scores), function(i) dat[["death"]])
      modnames <- sapply(grid[[i]], function(x) {
          if(!is.null(cfg[[x]][["full_name"]])) return(cfg[[x]][["full_name"]])
          return(x)
      })

      eval <- evalmod(
        scores = scores,
        labels = labels,
        dsids = 1:length(scores),
        modnames = modnames
      )

      plot.title <- stringr::str_to_upper(data_src)
      auc.str <- paste(round(auc(eval)$aucs[c(TRUE, FALSE)], 3), collapse = ",")
      plot.title <- paste0(plot.title, " ", sys, ": ((", auc.str, "))")

      if (!split_waves & res_type == "tabular") {

        ret <- auc(eval)
        ret <- ret[curvetypes == "ROC"]
        ret <- ret[, c("modnames", "aucs"), with = FALSE]
        ret <- data.table::setnames(ret, c("modnames", "aucs"), c("feature", "AUROC"))

        return(ret)

      } else return(autoplot(eval, "ROC") + geom_line(size = 2) + ggtitle(plot.title))

    }


  grid <- list(
    bone_marrow = c("lymph", "neut", "eos", "basos", "wbc", "nlr"), #
    renal = c("bun", "crea", "na", "k", "cl", "mg"),
    liver = c("bili", "alb", "phos", "alp"),
    age = c("age", "n_comb", "gcs"), #
    metabolic = c("lact", "bicar", "ph", "ldh"), #
    coag = c("d_dim", "ptt", "inr_pt", "plt", "fgn"), #
    resp = c("vent_ratio", "lung_comp", "dead_space", "pafi", "safi", "po2",
             "pco2", "o2sat", "resp"), #
    inflam = c("crp", "prcl", "temp") #
  )

  plots <- list()

  for (i in 1:length(grid)) {

    sys <- names(grid)[i]

    plots[[data_src]][[sys]] <- sys_perf(dat, split_waves = split_waves, res_type = res_type)

  }

  plots

}
