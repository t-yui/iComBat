source("icombat.R")
source("utils.R")

suppressPackageStartupMessages({
  library(sva)
  library(matrixStats)
  library(BiocParallel)
  library(jsonlite)
})

get_pvalues_age <- function(data, group, age) {
  pvals <- apply(data, 1, function(x) {
    fit <- summary(lm(x ~ group + age))
    pf <- pf((fit$coefficients[2, "t value"])^2, 1, fit$df[2], lower.tail = FALSE)
    return(pf)
  })
  return(pvals)
}

# data generation
generateDataFlex <- function(
  G = 500,
  batch_sizes,
  batch_means,
  batch_sds,
  case_prob,
  global_mean = 5,
  global_sd = 1,
  Delta_value = 0.5,
  age_mu = 50,
  age_sd = 10,
  age_beta = 0.02,
  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  total_samples <- sum(batch_sizes)
  dat <- matrix(NA, nrow = G, ncol = total_samples)
  batch <- numeric(total_samples)
  group <- numeric(total_samples)
  age <- rnorm(total_samples, mean = age_mu, sd = age_sd)
  
  Delta <- c(rep(Delta_value, min(50, G)), rep(0, G - min(50, G)))
  Age_beta <- c(rep(age_beta, G))
  
  start_idx <- 1
  for (k in seq_along(batch_sizes)) {
    nk <- batch_sizes[k]
    idx <- start_idx:(start_idx + nk - 1)
    grp <- rbinom(nk, 1, prob = case_prob[k])
    base <- matrix(rnorm(G * nk, mean = global_mean, sd = global_sd), nrow = G)
    dat[, idx] <- (base + batch_means[k] + Delta %*% t(as.matrix(grp)) +
                     Age_beta %*% t(as.matrix(age[idx]))) * batch_sds[k]
    batch[idx] <- k
    group[idx] <- grp
    start_idx <- start_idx + nk
  }
  list(dat = dat,
       batch = factor(batch),
       group = factor(group, levels = c(0, 1)),
       age = age)
}

# scenario settings
G <- 500
Delta_value <- 0.5
global_mean <- 5
global_sd <- 1
age_mu <- 50
age_sd <- 10

systematic_scenarios <- list(
  # Basic scenario
  S01 = list(
    name = "Baseline",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 3
  ),
  
  # Case-Control imbalance scenarios
  S02 = list(
    name = "Mild imbalance",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.3, 0.7, 0.4, 0.6),
    age_beta = 0.01,
    n_old_batches = 3
  ),
  
  S03 = list(
    name = "Extreme imbalance",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.1, 0.9, 0.2, 0.8),
    age_beta = 0.01,
    n_old_batches = 3
  ),
  
  # Sample size variations
  S04 = list(
    name = "Small samples",
    batch_sizes = c(25, 15, 20, 15),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 3
  ),
  
  S05 = list(
    name = "Large samples",
    batch_sizes = c(250, 150, 200, 150),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 3
  ),

  S06 = list(
    name = "Mixed sizes",
    batch_sizes = c(200, 10, 50, 100),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 3
  ),
  
  # Batch number variations
  S07 = list(
    name = "Few batches (2+1)",
    batch_sizes = c(50, 30, 40),
    batch_means = c(0, 2.5, -2),
    batch_sds = c(1, 1.5, 1.2),
    case_prob = c(0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 2
  ),
  
  S08 = list(
    name = "Many batches (6+3)",
    batch_sizes = c(20, 20, 20, 20, 20, 20, 15, 15, 15),
    batch_means = c(0, 0.5, 1, -0.5, -1, 1.5, -2, 2.5, -3),
    batch_sds = c(1, 1.1, 1.2, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6),
    case_prob = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    age_beta = 0.01,
    n_old_batches = 6
  ),
  
  # Strong covariate effect
  S09 = list(
    name = "Strong covariate",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 0.1,
    n_old_batches = 3
  ),

  S10 = list(
    name = "Extremely Strong covariate",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.5, 0.5, 0.5, 0.5),
    age_beta = 1.0,
    n_old_batches = 3
  ),
  
  # Complex scenarios
  S11 = list(
    name = "Imbalance + Strong covariate",
    batch_sizes = c(50, 30, 40, 30),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.3, 0.7, 0.4, 0.6),
    age_beta = 0.10,
    n_old_batches = 3
  ),
  
  S12 = list(
    name = "Imbalance + Strong covariate + Small samples",
    batch_sizes = c(25, 15, 20, 15),
    batch_means = c(0, 2.5, -2, -4),
    batch_sds = c(1, 1.5, 1.2, 1.8),
    case_prob = c(0.3, 0.7, 0.4, 0.6),
    age_beta = 0.10,
    n_old_batches = 3
  ),
  
  S13 = list(
    name = "Imbalance + Strong covariate + Small samples + Few batches",
    batch_sizes = c(25, 15, 20),
    batch_means = c(0, 2.5, -2),
    batch_sds = c(1, 1.5, 1.2),
    case_prob = c(0.3, 0.7, 0.4),
    age_beta = 0.10,
    n_old_batches = 2
  )
)

run_scenario_simulation <- function(scenario, n_rep = 1000, alpha = 0.05) {
  cat(sprintf("\nRunning scenario: %s\n", scenario$name))
  
  TPR_raw <- FPR_raw <- numeric(n_rep)
  TPR_full <- FPR_full <- numeric(n_rep)
  TPR_incr <- FPR_incr <- numeric(n_rep)
  lambda_raw <- lambda_full <- lambda_incr <- numeric(n_rep)
  nsv_raw <- nsv_full <- nsv_incr <- numeric(n_rep)
  correlation_vals <- numeric(n_rep)
  
  n_old <- scenario$n_old_batches
  batch_sizes_old <- scenario$batch_sizes[1:n_old]
  batch_means_old <- scenario$batch_means[1:n_old]
  batch_sds_old <- scenario$batch_sds[1:n_old]
  case_prob_old <- scenario$case_prob[1:n_old]
  
  batch_sizes_new <- scenario$batch_sizes[(n_old + 1):length(scenario$batch_sizes)]
  batch_means_new <- scenario$batch_means[(n_old + 1):length(scenario$batch_means)]
  batch_sds_new <- scenario$batch_sds[(n_old + 1):length(scenario$batch_sds)]
  case_prob_new <- scenario$case_prob[(n_old + 1):length(scenario$case_prob)]
  
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3)
  for (rep in seq_len(n_rep)) {
    setTxtProgressBar(pb, rep)

    old <- generateDataFlex(
      G = G,
      batch_sizes = batch_sizes_old,
      batch_means = batch_means_old,
      batch_sds = batch_sds_old,
      case_prob = case_prob_old,
      global_mean = global_mean,
      global_sd = global_sd,
      Delta_value = Delta_value,
      age_mu = age_mu,
      age_sd = age_sd,
      age_beta = scenario$age_beta,
      seed = 10000 * rep
    )

    newb <- generateDataFlex(
      G = G,
      batch_sizes = batch_sizes_new,
      batch_means = batch_means_new,
      batch_sds = batch_sds_new,
      case_prob = case_prob_new,
      global_mean = global_mean,
      global_sd = global_sd,
      Delta_value = Delta_value,
      age_mu = age_mu,
      age_sd = age_sd,
      age_beta = scenario$age_beta,
      seed = 20000 * rep
    )

    dat.raw <- cbind(old$dat, newb$dat)
    
    new_batch_labels <- character()
    for (i in seq_along(batch_sizes_new)) {
      new_batch_labels <- c(new_batch_labels, 
                            rep(as.character(n_old + i), batch_sizes_new[i]))
    }
    batch <- factor(c(as.character(old$batch), new_batch_labels))
    
    group <- factor(c(as.character(old$group), as.character(newb$group)), levels = c(0, 1))
    age <- c(old$age, newb$age)
    
    # (a) test for raw data
    p_raw <- get_pvalues_age(dat.raw, group, age)
    TPR_raw[rep] <- mean(p_raw[1:50] < alpha)
    FPR_raw[rep] <- mean(p_raw[-(1:50)] < alpha)
    lambda_raw[rep] <- gc_lambda(p_raw)
    
    # ComBat
    mod_full <- model.matrix(~ group + age)
    tryCatch({
      dat_full <- ComBat(dat = dat.raw, batch = batch, mod = mod_full, 
                         par.prior = TRUE, mean.only = FALSE)
      p_full <- get_pvalues_age(dat_full, group, age)
      TPR_full[rep] <- mean(p_full[1:50] < alpha)
      FPR_full[rep] <- mean(p_full[-(1:50)] < alpha)
      lambda_full[rep] <- gc_lambda(p_full)
      nsv_full[rep] <- get_n_sva(dat_full, mod_full)
    }, error = function(e) {
      cat("\nComBat error in rep", rep, ":", e$message, "\n")
      TPR_full[rep] <- NA
      FPR_full[rep] <- NA
      lambda_full[rep] <- NA
      nsv_full[rep] <- NA
    })
    
    # iComBat
    n_old_samples <- ncol(old$dat)
    group_old <- group[1:n_old_samples]
    age_old <- age[1:n_old_samples]
    mod_old <- model.matrix(~ group_old + age_old)
    
    tryCatch({
      fit_old <- ComBatInitialize(dat = old$dat, batch = old$batch, 
                                  mod = mod_old, mean.only = FALSE)
      dat_old_corr <- ComBat(dat = old$dat, batch = old$batch, 
                             mod = mod_old, par.prior = TRUE, mean.only = FALSE)

      group_new <- group[(n_old_samples + 1):length(group)]
      age_new <- age[(n_old_samples + 1):length(age)]
      mod_new <- model.matrix(~ group_new + age_new)

      dat_new_corr_list <- list()
      start_idx <- 1
      for (i in seq_along(batch_sizes_new)) {
        end_idx <- start_idx + batch_sizes_new[i] - 1
        batch_data <- newb$dat[, start_idx:end_idx, drop = FALSE]
        batch_label <- factor(rep(n_old + i, batch_sizes_new[i]))
        
        dat_new_corr_list[[i]] <- ComBatIncremental(
          dat.new = batch_data,
          new.batch = batch_label,
          fit = fit_old,
          par.prior = TRUE,
          mod.new = mod_new[start_idx:end_idx, , drop = FALSE]
        )
        start_idx <- end_idx + 1
      }
      dat_new_corr <- do.call(cbind, dat_new_corr_list)
      
      dat_incr <- cbind(dat_old_corr, dat_new_corr)
      p_incr <- get_pvalues_age(dat_incr, group, age)
      TPR_incr[rep] <- mean(p_incr[1:50] < alpha)
      FPR_incr[rep] <- mean(p_incr[-(1:50)] < alpha)
      lambda_incr[rep] <- gc_lambda(p_incr)
      nsv_incr[rep] <- get_n_sva(dat_incr, mod_full)

      if (!is.na(TPR_full[rep])) {
        correlation_vals[rep] <- cor(dat_full[, 1], dat_incr[, 1])
      }
      
    }, error = function(e) {
      cat("\niComBat error in rep", rep, ":", e$message, "\n")
      TPR_incr[rep] <- NA
      FPR_incr[rep] <- NA
      lambda_incr[rep] <- NA
      nsv_incr[rep] <- NA
      correlation_vals[rep] <- NA
    })

    nsv_raw[rep] <- get_n_sva(dat.raw, mod_full)
  }
  close(pb)
  
  # aggregate results
  results <- list(
    scenario_name = scenario$name,
    n_samples_total = sum(scenario$batch_sizes),
    n_batches = length(scenario$batch_sizes),
    n_old_batches = scenario$n_old_batches,
    n_new_batches = length(scenario$batch_sizes) - scenario$n_old_batches,

    TPR_raw = mean(TPR_raw, na.rm = TRUE),
    TPR_full = mean(TPR_full, na.rm = TRUE),
    TPR_incr = mean(TPR_incr, na.rm = TRUE),

    FPR_raw = mean(FPR_raw, na.rm = TRUE),
    FPR_full = mean(FPR_full, na.rm = TRUE),
    FPR_incr = mean(FPR_incr, na.rm = TRUE),

    lambda_raw = mean(lambda_raw, na.rm = TRUE),
    lambda_full = mean(lambda_full, na.rm = TRUE),
    lambda_incr = mean(lambda_incr, na.rm = TRUE),

    nsv_raw = mean(nsv_raw, na.rm = TRUE),
    nsv_full = mean(nsv_full, na.rm = TRUE),
    nsv_incr = mean(nsv_incr, na.rm = TRUE),

    combat_icombat_cor = mean(correlation_vals, na.rm = TRUE),

    success_rate_combat = mean(!is.na(TPR_full)),
    success_rate_icombat = mean(!is.na(TPR_incr)),

    raw_data = list(
      TPR_raw = TPR_raw,
      TPR_full = TPR_full,
      TPR_incr = TPR_incr,
      FPR_raw = FPR_raw,
      FPR_full = FPR_full,
      FPR_incr = FPR_incr,
      lambda_raw = lambda_raw,
      lambda_full = lambda_full,
      lambda_incr = lambda_incr,
      nsv_raw = nsv_raw,
      nsv_full = nsv_full,
      nsv_incr = nsv_incr,
      correlation_vals = correlation_vals
    )
  )
  return(results)
}

run_all_scenarios <- function(scenarios = systematic_scenarios, n_rep = 1000, alpha = 0.05) {
  all_results <- list()
  cat("=== Running All Scenarios ===\n")
  cat(sprintf("Number of scenarios: %d\n", length(scenarios)))
  cat(sprintf("Repetitions per scenario: %d\n", n_rep))
  cat(sprintf("Significance level: %.3f\n\n", alpha))
  
  start_time <- Sys.time()
  for (i in seq_along(scenarios)) {
    scenario_name <- names(scenarios)[i]
    scenario <- scenarios[[scenario_name]]
    cat(sprintf("\n[%d/%d] ", i, length(scenarios)))
    result <- run_scenario_simulation(scenario, n_rep = n_rep, alpha = alpha)
    all_results[[scenario_name]] <- result
  }
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")

  cat(sprintf("\n\nTotal execution time: %.1f minutes\n", as.numeric(total_time)))
  return(all_results)
}

compare_scenarios <- function(all_results) {
  comparison_df <- do.call(rbind, lapply(names(all_results), function(name) {
    res <- all_results[[name]]
    data.frame(
      Scenario = res$scenario_name,
      N_Total = res$n_samples_total,
      N_Batches = res$n_batches,
      TPR_Raw = res$TPR_raw,
      TPR_ComBat = res$TPR_full,
      TPR_iComBat = res$TPR_incr,
      FPR_Raw = res$FPR_raw,
      FPR_ComBat = res$FPR_full,
      FPR_iComBat = res$FPR_incr,
      Lambda_Raw = res$lambda_raw,
      Lambda_ComBat = res$lambda_full,
      Lambda_iComBat = res$lambda_incr,
      nSV_Raw = res$nsv_raw,
      nSV_ComBat = res$nsv_full,
      nSV_iComBat = res$nsv_incr,
      Correlation = res$combat_icombat_cor,
      stringsAsFactors = FALSE
    )
  }))
  return(comparison_df)
}


all_results <- run_all_scenarios(n_rep = 1000)
save(all_results, file = "all_results.RData")
write_json(
  all_results,
  path       = "all_results.json",
  pretty     = TRUE,
  auto_unbox = TRUE
)
