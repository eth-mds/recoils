auc_optimizer <- function(train_data, cfg, systems, ...) {

  best <- lapply(systems, function(x) list(auc = 0.5))
  names(best) <- systems
  comp_points <- list(
    age = 6L, renal = 4L, coag = 4L, resp = 4L, metabolic = 4L
  )

  for (sys in systems) {

    components <- NULL
    for (i in 1:length(cfg))
      if(cfg[[i]][["category"]] == sys) components <- append(components,
                                                             names(cfg)[i])

      for (cp in components) {

        print(cp)
        cp.idx <- grep(paste0("^", cp, "[.]"), names(train_data))
        mat <- as.matrix(train_data[, cp.idx, with=FALSE])
        cpt <- replicate(500, sample(1:ncol(mat), comp_points[[sys]], F))

        for (k in 1:ncol(cpt)) {

          auc <- PRROC::roc.curve(weights.class0 = train_data[["death"]],
            scores.class0 = rowSums(mat[, cpt[, k]]))$auc

          if (auc > best[[sys]][["auc"]]) {

            best[[sys]][["auc"]] <- auc
            best[[sys]][["feature"]] <- cp
            best[[sys]][["cols"]] <- colnames(mat)[cpt[, k]] # careful here

            eval <- evalmod(
              scores = list(rowSums(mat[, cpt[, k]])),#, replace_na(sofa.cp, 0)),
              labels = list(train_data[["death"]]),#, train_data[["death"]]),
              #dsids = 1:2,
              #modnames = paste(c("DOSE", sofa.name), sys)
              modnames = paste("ECOS", sys)
            )

            best[[sys]][["plot"]]<- autoplot(eval, "ROC") + geom_line(size = 2) +
              theme(legend.position = "bottom",
                legend.text = element_text(size=15),
                plot.title = element_blank())

          }

        }

      }

  }

  best

}

running_decorr <- function(train_data, cfg, score, max_epoch = 50, ...) {

  sys_components <- list()
  for (i in 1:length(cfg)) {

    catg <- cfg[[i]][["category"]]
    sys_components[[catg]] <- c(sys_components[[catg]], names(cfg)[i])

  }

  run_auc <- PRROC::roc.curve(weights.class0 = train_data[["death"]],
    scores.class0 = rowSums(train_data[, unlist(score), with = F]))$auc

  for (epoch in 1:max_epoch) {

    change <- F

    for (sys in names(score)) {

      smsys <- rowSums(train_data[, unlist(score[-which(names(score) == sys)]),
                                  with=F])

      for (cp in sys_components[[sys]]) {

        print(cp)
        cp.idx <- grep(paste0("^", cp, "[.]"), names(train_data))
        mat <- as.matrix(train_data[, cp.idx, with=FALSE])
        cpt <- replicate(500, sample(1:ncol(mat), 4, F))

        for (k in 1:ncol(cpt)) {

          auc <- PRROC::roc.curve(
            weights.class0 = train_data[["death"]],
            scores.class0 = smsys + rowSums(mat[, cpt[, k]]))$auc

          if (auc > run_auc) {

            cat("AUC improves from", round(run_auc, 4), "to", round(auc, 4),
                "\n")
            run_auc <- auc
            score[[sys]] <- colnames(mat)[cpt[, k]]
            change <- T

          }

        }

        nuke_auc <- PRROC::roc.curve(
          weights.class0 = train_data[["death"]],
          scores.class0 = smsys)$auc

        if (nuke_auc > run_auc) {

          cat("Nuke the", sys, "component\n")
          score[[sys]] <- character(0)

        }

      }

    }


    if (!change) break

  }

  score

}

scaling_factor <- function(X, y, ids, by_hosp = FALSE, seed = 7L, nfolds = 5L) {

  set.seed(seed)
  n_target <- round(nrow(X) * 0.6)

  if (by_hosp) {

    hosp_select <- config("hospital-split")[["dev"]]
    dev <- id_col(ids[hosp %in% hosp_select])

  } else dev <- sample(id_col(ids), n_target)

  dev_idx <- (id_col(ids) %in% dev)
  X <- X[dev_idx, ]
  y <- y[dev_idx]

  folds <- sample(1:nrow(X))

  folds <- rep(1:nfolds, ceiling(nrow(X) / nfolds))[1:nrow(X)][folds]

  res <- cv.glmnet(X, y, family = "binomial", foldid = folds, type = "auc")

  df <- NULL
  for (i in 1:length(res$lambda)) {

    lam <- res$lambda[i]

    for (fold in 1:nfolds) {

      ind <- (folds != fold)
      model <- glmnet(X[ind, ], y[ind], family = "binomial", lambda = lam)
      beta <- as.vector(model$beta)

      for (m in c(3, 4, 5, 10^6)) {

        auc <- PRROC::roc.curve(
          scores.class0 = X[!ind, ] %*% round(beta * m),
          weights.class0 = y[!ind]
        )$auc

        df <- rbind(df, c(i, fold, m, auc))

      }

    }

  }


  df <- data.table::as.data.table(df)
  data.table::setnames(df, names(df), c("iteration", "fold", "m", "auc"))

  df[, mean_auc := mean(auc), by = c("m", "iteration")]

  plot_df <- unique(df[, c("iteration", "m", "mean_auc")])

  pp <- ggplot(plot_df, aes(x = iteration, y = mean_auc, color = factor(m))) +
    geom_line() + theme_bw() + ggtitle("Performance paths by scaling factor") +
    scale_color_discrete(name = "Scaling factor") +
    theme(
      legend.position = c(0.8, 0.8),
      legend.box.background = element_rect()
    ) + ylim(c(0.6, 0.8)) + geom_vline(xintercept = 30, linetype = "dashed",
                                       color = "red")

  print(pp)

  # interactive
  model_ind <- as.integer(readline("choose model: "))
  beta <- res$glmnet.fit$beta[, model_ind]
  score <- round(beta * 4)

  return(list(score = score, beta = beta, dev = dev,
              val = setdiff(id_col(ids), dev)))

}
