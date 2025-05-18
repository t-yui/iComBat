source("/path_to_directory/iComBat/scripts/icombat.R")

# data generation
generateTestDataExperiment <- function(
    G = 500,
    batch_sizes,
    batch_means,
    batch_sds,
    global_mean = 5,
    global_sd = 1,
    Delta_value = 0.5,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  total_samples <- sum(batch_sizes)
  dat <- matrix(NA, nrow = G, ncol = total_samples)
  batch <- numeric(total_samples)
  group <- numeric(total_samples)
  start_idx <- 1
  Delta <- c(rep(Delta_value, min(50, G)), rep(0, G - min(50, G))) # 群効果: 最初の50遺伝子のみ
  for (k in seq_along(batch_sizes)) {
    nk <- batch_sizes[k]
    end_idx <- start_idx + nk - 1
    grp <- sample(rep(0:1, length.out = nk))
    base_expr <- matrix(rnorm(G * nk, mean = global_mean, sd = global_sd), nrow = G, ncol = nk)
    dat[, start_idx:end_idx] <- (base_expr + batch_means[k] + Delta %*% t(as.matrix(grp))) * batch_sds[k]
    batch[start_idx:end_idx] <- k
    group[start_idx:end_idx] <- grp
    start_idx <- end_idx + 1
  }
  batch <- factor(batch, levels = seq_along(batch_sizes))
  group <- factor(group, levels = c(0, 1))
  list(dat = dat, batch = batch, group = group, true_G = G, true_N = total_samples)
}

# get p values by welch's t-test
get_pvalues <- function(data, group) {
  pvals <- apply(data, 1, function(x) {
    t.test(x ~ group)$p.value
  })
  return(pvals)
}

# simulation study
n_rep <- 1000

TPR_raw_vec <- numeric(n_rep)
FPR_raw_vec <- numeric(n_rep)
TPR_full_vec <- numeric(n_rep)
FPR_full_vec <- numeric(n_rep)
TPR_incr_vec <- numeric(n_rep)
FPR_incr_vec <- numeric(n_rep)

batch_sizes_old <- c(50, 30, 40)
batch_means_old <- c(0, 2.5, -2)
batch_sds_old <- c(1, 1.5, 1.2)
batch_sizes_new <- c(30)
batch_means_new <- c(-4)
batch_sds_new <- c(1.8)

alpha <- 0.05 # significance level

for (i in 1:n_rep) {
  # data generation for existing batches
  synthetic_old <- generateTestDataExperiment(
    G = 500,
    batch_sizes = batch_sizes_old,
    batch_means = batch_means_old,
    batch_sds = batch_sds_old,
    global_mean = 5,
    global_sd = 1,
    seed = 1000 + i
  )
  dat.old <- synthetic_old$dat
  batch.old <- synthetic_old$batch
  group.old <- synthetic_old$group

  # data generation for new batches
  synthetic_new <- generateTestDataExperiment(
    G = 500,
    batch_sizes = batch_sizes_new,
    batch_means = batch_means_new,
    batch_sds = batch_sds_new,
    global_mean = 5,
    global_sd = 1,
    seed = 2000 + i
  )
  dat.new <- synthetic_new$dat
  batch.new <- factor(rep(4, ncol(dat.new)), levels = 1:4)
  group.new <- synthetic_new$group

  # non-corrected combined data
  dat.raw.combined <- cbind(dat.old, dat.new)
  batch.combined <- factor(c(as.character(batch.old), as.character(batch.new)), levels = 1:4)
  group.combined <- factor(c(as.character(group.old), as.character(group.new)), levels = c(0, 1))

  # (a) results: no correction
  p_raw <- get_pvalues(dat.raw.combined, group.combined)

  # (b) results: combat
  mod_full <- model.matrix(~group, data = data.frame(group = group.combined))
  corrected_all <- ComBat(dat = dat.raw.combined, batch = batch.combined, mod = mod_full)
  p_full <- get_pvalues(corrected_all, group.combined)

  # (c) results: icombat
  mod_old <- model.matrix(~group, data = data.frame(group = group.old))
  fit.old <- ComBatInitialize(
    dat = dat.old, batch = batch.old, mod = mod_old,
    par.prior = TRUE, mean.only = FALSE, ref.batch = NULL
  )
  corrected_old <- ComBat(dat = dat.old, batch = batch.old, mod = mod_old)
  mod_new <- model.matrix(~group, data = data.frame(group = group.new))
  corrected_new <- ComBatIncremental(
    dat.new = dat.new, new.batch = batch.new,
    fit = fit.old, mod.new = mod_new
  )
  corrected_incr <- cbind(corrected_old, corrected_new)
  group_incr <- factor(c(as.character(group.old), as.character(group.new)), levels = c(0, 1))
  p_incr <- get_pvalues(corrected_incr, group_incr)

  TPR_raw_vec[i] <- sum(p_raw[1:50] < alpha) / 50
  FPR_raw_vec[i] <- sum(p_raw[51:500] < alpha) / 450

  TPR_full_vec[i] <- sum(p_full[1:50] < alpha) / 50
  FPR_full_vec[i] <- sum(p_full[51:500] < alpha) / 450

  TPR_incr_vec[i] <- sum(p_incr[1:50] < alpha) / 50
  FPR_incr_vec[i] <- sum(p_incr[51:500] < alpha) / 450
}

cat("No Correction: Avg.TPR =", mean(TPR_raw_vec), " Avg.FPR =", mean(FPR_raw_vec), "\n")
cat("Corrected by ComBat: Avg.TPR =", mean(TPR_full_vec), " Avg.FPR =", mean(FPR_full_vec), "\n")
cat("Corrected by iComBat: Avg.TPR =", mean(TPR_incr_vec), " Avg.FPR =", mean(FPR_incr_vec), "\n")

# PCA plot for last generated and corrected data

# (a) no correction
pca_raw <- prcomp(t(dat.raw.combined), scale. = FALSE)
par(mfrow = c(1, 3))
plot(pca_raw$x[, 1], pca_raw$x[, 2],
  col = as.factor(batch.combined),
  pch = 19, xlab = "PC1", ylab = "PC2",
  main = "No Correction"
)
legend("topright", legend = levels(batch.combined), pch = 19, col = 1:length(levels(batch.combined)))

# (b) combat
pca_full <- prcomp(t(corrected_all), scale. = FALSE)
plot(pca_full$x[, 1], pca_full$x[, 2],
  col = as.factor(batch.combined),
  pch = 19, xlab = "PC1", ylab = "PC2",
  main = "Corrected by ComBat"
)
legend("topright", legend = levels(batch.combined), pch = 19, col = 1:length(levels(batch.combined)))

# (c) icombat
pca_incr <- prcomp(t(corrected_incr), scale. = FALSE)
plot(pca_incr$x[, 1], pca_incr$x[, 2],
  col = as.factor(c(as.character(batch.old), rep(4, ncol(corrected_new)))),
  pch = 19, xlab = "PC1", ylab = "PC2",
  main = "Corrected by iComBat"
)
legend("topright",
  legend = unique(c(as.character(batch.old), 4)),
  pch = 19, col = 1:length(unique(c(as.character(batch.old), 4)))
)
par(mfrow = c(1, 1))
