source("icombat.R")
source("utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
  library(limma)
  library(ggplot2)
  library(sva)
  library(viridis)
  library(patchwork)
  library(hexbin)
  library(dplyr)
  library(reshape2)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# load data
meth_data <- fread("GSE42861_processed_methylation_matrix.txt")
probe_names <- meth_data$ID_REF
meth_matrix <- as.matrix(meth_data[, -1])
rownames(meth_matrix) <- probe_names
colnames(meth_matrix) <- colnames(meth_data)[-1]

# exclude chrX/Y probe
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
xy_probes <- rownames(anno_450k)[anno_450k$chr %in% c("chrX", "chrY")]
to_drop <- intersect(rownames(meth_matrix), xy_probes)
meth_matrix <- meth_matrix[setdiff(rownames(meth_matrix), to_drop), , drop = FALSE]
probe_names <- rownames(meth_matrix)
cat(sprintf("Removed %d X/Y probes; %d probes remain.\n",
            length(to_drop), nrow(meth_matrix)))

# convert Beta to M-value
epsilon <- 1e-6
meth_clipped <- pmin(pmax(meth_matrix, epsilon), 1 - epsilon)
mval_matrix <- log2(meth_clipped / (1 - meth_clipped))
rm(meth_matrix, meth_clipped); gc()

# get sample information from soft file
soft_lines <- readLines("GSE42861_family.soft")
sample_indices <- grep("^\\^SAMPLE = ", soft_lines)
sample_info <- data.frame(
  sample_id = character(),
  disease_state = character(),
  smoking_status = character(),
  age = character(),
  gender = character(),
  stringsAsFactors = FALSE
)
gsm_sentrix_map <- data.frame(gsm = character(), sentrix = character(), stringsAsFactors = FALSE)

for (i in seq_along(sample_indices)) {
  start <- sample_indices[i]
  end <- if (i < length(sample_indices)) sample_indices[i + 1] - 1 else length(soft_lines)
  block <- soft_lines[start:end]

  sample_id <- sub("^\\^SAMPLE = ", "", block[1])

  disease_line <- grep("!Sample_characteristics_ch1 = disease state:", block, value = TRUE)
  disease_state <- if (length(disease_line) > 0) sub("!Sample_characteristics_ch1 = disease state: ", "", disease_line) else NA

  smoking_line <- grep("!Sample_characteristics_ch1 = smoking status:", block, value = TRUE)
  smoking_status <- if (length(smoking_line) > 0) sub("!Sample_characteristics_ch1 = smoking status: ", "", smoking_line) else NA

  age_line <- grep("!Sample_characteristics_ch1 = age:", block, value = TRUE)
  age <- if (length(age_line) > 0) sub("!Sample_characteristics_ch1 = age: ", "", age_line) else NA
  
  gender_line <- grep("!Sample_characteristics_ch1 = gender:", block, value = TRUE)
  gender <- if (length(gender_line) > 0) sub("!Sample_characteristics_ch1 = gender: ", "", gender_line) else NA
  
  sample_info <- rbind(
    sample_info,
    data.frame(
      sample_id = sample_id,
      disease_state = disease_state,
      smoking_status = smoking_status,
      age = age,
      gender = gender,
      stringsAsFactors = FALSE
    )
  )
  sentrix_line <- grep("_Grn.idat.gz", block, value = TRUE)
  sentrix_id <- if (length(sentrix_line) > 0) sub(".*_(\\d{10}_R\\d{2}C\\d{2})_Grn\\.idat\\.gz", "\\1", sentrix_line[1]) else NA
  gsm_sentrix_map <- rbind(gsm_sentrix_map, data.frame(gsm = sample_id, sentrix = sentrix_id, stringsAsFactors = FALSE))
}
sample_info$smoking_history <- ifelse(sample_info$smoking_status %in% c("current", "ex"), "yes", "no")
sample_info <- merge(sample_info, gsm_sentrix_map, by.x = "sample_id", by.y = "gsm", all.x = TRUE)
pheno_soft <- unique(sample_info[, c("sample_id","sentrix","disease_state","smoking_status","age","gender","smoking_history")])

# create batch
sample_ids <- colnames(mval_matrix)
batches <- sub("_.*", "", sample_ids)
pheno_data <- data.frame(sample_id = sample_ids, batch = batches, stringsAsFactors = FALSE)
pheno_data <- merge(pheno_data, pheno_soft[, c("sentrix","disease_state","smoking_status","age","gender","smoking_history")],
                    by.x = "sample_id", by.y = "sentrix", all.x = TRUE)
pheno_data <- subset(pheno_data, sample_id != "Pval")
pheno_data <- pheno_data[match(colnames(mval_matrix), pheno_data$sample_id), ]

# remove probes with zero variance
batches_fac <- factor(pheno_data$batch)
keep_probe <- apply(mval_matrix, 1, function(vals) {
  vars <- tapply(vals, batches_fac, var)
  all(vars > 0 & !is.na(vars))
})
saveRDS(keep_probe, file = "GSE42861_keep_probe.rds", compress = "xz")
mval_matrix <- mval_matrix[keep_probe, ]
cat(sprintf("Probes kept after per-batch variance filter: %d\n", sum(keep_probe)))

# create model of with covariates
pheno_cb <- pheno_data
pheno_cb$disease_state   <- factor(pheno_cb$disease_state, levels = c("Normal", "rheumatoid arthritis"))
pheno_cb$smoking_history <- factor(pheno_cb$smoking_history, levels = c("no","yes"))
pheno_cb$gender          <- factor(pheno_cb$gender, levels = c("f","m"))
pheno_cb$age             <- suppressWarnings(as.numeric(pheno_cb$age))
keep_cb <- complete.cases(pheno_cb[, c("disease_state","age","smoking_history","gender","batch")])
mval_cb  <- mval_matrix[, keep_cb]
batch_cb <- droplevels(factor(pheno_cb$batch[keep_cb]))
pheno_cb <- pheno_cb[keep_cb, ]
mod_cb <- model.matrix(~ disease_state + age + smoking_history + gender, data = pheno_cb)

# perform ComBat
set.seed(1)
corrected_full <- ComBat(
  dat = mval_cb,
  batch = batch_cb,
  mod = mod_cb,
  par.prior = TRUE,
  mean.only = FALSE
)

# perform iComBat
b_levels <- levels(batch_cb)
existing_batches <- grep("^5", b_levels, value = TRUE)
new_batches <- grep("^7", b_levels, value = TRUE)
n_b <- length(b_levels)

existing_samples <- pheno_cb$sample_id[batch_cb %in% existing_batches]
new_samples      <- pheno_cb$sample_id[batch_cb %in% new_batches]

dat.old   <- mval_cb[, pheno_cb$sample_id %in% existing_samples]
dat.new   <- mval_cb[, pheno_cb$sample_id %in% new_samples]
batch.old <- droplevels(batch_cb[pheno_cb$sample_id %in% existing_samples])
batch.new <- droplevels(batch_cb[pheno_cb$sample_id %in% new_samples])

mod.old <- model.matrix(~ disease_state + age + smoking_history + gender,
                        data = pheno_cb[pheno_cb$sample_id %in% existing_samples, ])
mod.new <- model.matrix(~ disease_state + age + smoking_history + gender,
                        data = pheno_cb[pheno_cb$sample_id %in% new_samples, ])

corrected_old <- ComBat(dat = dat.old, batch = batch.old, mod = mod.old,
                        par.prior = TRUE, mean.only = FALSE)

## initializeation
fit.old <- ComBatInitialize(
  dat = dat.old, batch = batch.old, mod = mod.old,
  mean.only = FALSE
)

# incremantally add new batches
new_batch_levels <- levels(batch.new)
corrected_new_list <- list()
if (length(new_batch_levels) > 0) {
  pb <- txtProgressBar(min = 0, max = length(new_batch_levels), style = 3)
  for (i in seq_along(new_batch_levels)) {
    current_batch <- new_batch_levels[i]
    idx <- which(batch.new == current_batch)
    batch_data <- dat.new[, idx, drop = FALSE]
    batch_mod  <- mod.new[idx, , drop = FALSE]
    
    corrected_batch <- ComBatIncremental(
      dat.new = batch_data,
      new.batch = factor(rep(current_batch, ncol(batch_data))),
      fit = fit.old,
      mod.new = batch_mod
    )
    corrected_new_list[[i]] <- corrected_batch
    setTxtProgressBar(pb, i)
  }
  close(pb)
  corrected_new <- do.call(cbind, corrected_new_list)
  colnames(corrected_new) <- colnames(dat.new)
  corrected_incr <- cbind(corrected_old, corrected_new)
  corrected_incr <- corrected_incr[, match(pheno_cb$sample_id, colnames(corrected_incr))]
} else {
  corrected_incr <- corrected_old
  corrected_incr <- corrected_incr[, match(pheno_cb$sample_id, colnames(corrected_incr))]
}

# PCA
n_top_cpgs <- 10000
select_top_variable_cpgs <- function(dat, n_cpgs = 10000) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    cpg_var <- matrixStats::rowVars(dat, na.rm = TRUE)
  } else {
    cpg_var <- apply(dat, 1, var, na.rm = TRUE)
  }
  cpg_var <- cpg_var[!is.na(cpg_var) & cpg_var > 0]
  n_select <- min(n_cpgs, length(cpg_var))
  names(sort(cpg_var, decreasing = TRUE))[seq_len(n_select)]
}
top_cpgs <- select_top_variable_cpgs(mval_cb, n_cpgs = n_top_cpgs)

mval_cb_pca       <- mval_cb[top_cpgs, ]
corrected_full_pca <- corrected_full[top_cpgs, ]
corrected_incr_pca <- corrected_incr[top_cpgs, ]

pca_raw     <- prcomp(t(mval_cb_pca), scale. = FALSE, center = TRUE)
pca_combat  <- prcomp(t(corrected_full_pca), scale. = FALSE, center = TRUE)
pca_icombat <- prcomp(t(corrected_incr_pca), scale. = FALSE, center = TRUE)

var_raw <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2)) * 100
var_combat <- (pca_combat$sdev^2 / sum(pca_combat$sdev^2)) * 100
var_icombat <- (pca_icombat$sdev^2 / sum(pca_icombat$sdev^2)) * 100

n_batch <- length(levels(batch_cb))
batch_colors <- rainbow(n_batch)
batch_cb <- as.factor(sample_info$sentrix)

# data with PCs
pca_data_raw <- data.frame(PC1 = pca_raw$x[,1], PC2 = pca_raw$x[,2], PC3 = pca_raw$x[,3],
                           Batch = droplevels(as.factor(pheno_cb$batch)),
                           Sex = pheno_cb$gender,
                           Age = pheno_cb$age,
                           IsNew = pheno_cb$sample_id %in% new_samples)
pca_data_combat <- data.frame(PC1 = pca_combat$x[,1], PC2 = pca_combat$x[,2], PC3 = pca_combat$x[,3],
                              Batch = droplevels(as.factor(pheno_cb$batch)),
                              Sex = pheno_cb$gender,
                              Age = pheno_cb$age,
                              IsNew = pheno_cb$sample_id %in% new_samples)
pca_data_icombat <- data.frame(PC1 = pca_icombat$x[,1], PC2 = pca_icombat$x[,2], PC3 = pca_icombat$x[,3],
                               Batch = droplevels(as.factor(pheno_cb$batch)),
                               Sex = pheno_cb$gender,
                               Age = pheno_cb$age,
                               IsNew = pheno_cb$sample_id %in% new_samples)

all_batch_levels <- levels(droplevels(as.factor(pheno_cb$batch)))
existing_levels  <- intersect(all_batch_levels, existing_batches)
new_levels       <- intersect(all_batch_levels, new_batches)

# グラデーション（依存パッケージ不要・baseの色で作成）
reds  <- grDevices::colorRampPalette(
  c("#fee5d9","#fcae91","#fb6a4a","#de2d26","#a50f15")
)(max(1, length(existing_levels)))

blues <- grDevices::colorRampPalette(
  c("#deebf7","#9ecae1","#6baed6","#3182bd","#08519c")
)(max(1, length(new_levels)))

batch_palette <- c(stats::setNames(reds,  existing_levels),
                   stats::setNames(blues, new_levels))

# 凡例の並び（既存→新規）
batch_breaks <- c(existing_levels, new_levels)

p1 <- ggplot(pca_data_raw, aes(PC1, PC2, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), y = sprintf("PC2 (%.1f%%)", var_raw[2]), title = "Raw") +
  theme_light_plot(12) + theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p2 <- ggplot(pca_data_combat, aes(PC1, PC2, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat[1]), y = sprintf("PC2 (%.1f%%)", var_combat[2]), title = "ComBat") +
  theme_light_plot(12) + theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p3 <- ggplot(pca_data_icombat, aes(PC1, PC2, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat[1]), y = sprintf("PC2 (%.1f%%)", var_icombat[2]), title = "iComBat") +
  theme_light_plot(12) + theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p1_13 <- ggplot(pca_data_raw, aes(PC1, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), y = sprintf("PC3 (%.1f%%)", var_raw[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p2_13 <- ggplot(pca_data_combat, aes(PC1, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat[1]), y = sprintf("PC3 (%.1f%%)", var_combat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p3_13 <- ggplot(pca_data_icombat, aes(PC1, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat[1]), y = sprintf("PC3 (%.1f%%)", var_icombat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p1_23 <- ggplot(pca_data_raw, aes(PC2, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC2 (%.1f%%)", var_raw[2]), y = sprintf("PC3 (%.1f%%)", var_raw[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p2_23 <- ggplot(pca_data_combat, aes(PC2, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC2 (%.1f%%)", var_combat[2]), y = sprintf("PC3 (%.1f%%)", var_combat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p3_23 <- ggplot(pca_data_icombat, aes(PC2, PC3, color = Batch))+
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC2 (%.1f%%)", var_icombat[2]), y = sprintf("PC3 (%.1f%%)", var_icombat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

pca_combined <- (p1 | p2 | p3) / (p1_13 | p2_13 | p3_13) / (p1_23 | p2_23 | p3_23)
print(pca_combined)

legend_df <- data.frame(
  Batch = batch_breaks,
  Color = unname(batch_palette[batch_breaks]),
  Group = ifelse(batch_breaks %in% existing_batches, "Existing batches", "New batches"),
  stringsAsFactors = FALSE
)
ncol_per_group <- 5
split_list <- split(legend_df, legend_df$Group, drop = TRUE)
for (g in names(split_list)) {
  dfg <- split_list[[g]]
  n <- nrow(dfg)
  nrow_per_col <- ceiling(n / ncol_per_group)
  dfg$col <- ((seq_len(n) - 1) %/% nrow_per_col) + 1
  dfg$row <- ((seq_len(n) - 1) %%  nrow_per_col) + 1
  split_list[[g]] <- dfg
}
legend_df_layout <- do.call(rbind, split_list)
p_legend_vertical <-
  ggplot(legend_df_layout, aes(x = col, y = -row)) +
  geom_point(aes(color = Color), size = 3) +
  geom_text(aes(label = Batch), hjust = 0, nudge_x = 0.25, size = 2.8) +
  scale_color_identity() +
  facet_wrap(~ Group, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = NULL,
                     expand = expansion(mult = c(0.02, 0.50)))+
  scale_y_continuous(breaks = NULL) +
  labs(title = "Batch (Sentrix IDs)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.margin = margin(t = 6, r = 24, b = 6, l = 8),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_cartesian(clip = "off")
print(p_legend_vertical)

# test association with PCs
res_raw <- check_PC_covariate_association(pca_data_raw)
res_combat <- check_PC_covariate_association(pca_data_combat)
res_icombat <- check_PC_covariate_association(pca_data_icombat)
res_raw$PC1
res_combat$PC1
res_icombat$PC1

# sample to sample density plot
combat_values  <- as.vector(corrected_full)
icombat_values <- as.vector(corrected_incr)
cor_value <- cor(combat_values, icombat_values, use = "complete.obs")
cat(sprintf("\nCorrelation between ComBat and iComBat: r = %.4f\n", cor_value))

n_points <- length(combat_values)
if (n_points > 100000) {
  set.seed(1)
  idx <- sample.int(n_points, 100000)
  plot_combat  <- combat_values[idx]
  plot_icombat <- icombat_values[idx]
} else {
  plot_combat  <- combat_values
  plot_icombat <- icombat_values
}

density_plot <- ggplot(data.frame(ComBat = plot_combat, iComBat = plot_icombat), aes(ComBat, iComBat)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
  scale_fill_viridis(name = "Density", option = "plasma") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "ComBat corrected values", y = "iComBat corrected values",  title = sprintf("r = %.3f", cor_value)) +
  theme_light_plot() + theme(aspect.ratio = 1) + coord_fixed()
print(density_plot)

# EWAS analysis
pheno_full_subset <- pheno_cb
corrected_full_subset <- corrected_full
corrected_incr_subset <- corrected_incr
mval_subset <- mval_cb

design_full_disease <- model.matrix(~ disease_state + age + smoking_history + gender, data = pheno_full_subset)
design_full_smoke   <- model.matrix(~ smoking_history + age + disease_state + gender, data = pheno_full_subset)

## disease
fit_raw_disease <- eBayes(lmFit(mval_subset, design_full_disease))
fit_full_disease <- eBayes(lmFit(corrected_full_subset, design_full_disease))
fit_incr_disease <- eBayes(lmFit(corrected_incr_subset, design_full_disease))

results_raw_disease  <- topTable(fit_raw_disease,  coef = "disease_staterheumatoid arthritis", number = Inf, adjust.method = "BH")
results_full_disease <- topTable(fit_full_disease, coef = "disease_staterheumatoid arthritis", number = Inf, adjust.method = "BH")
results_incr_disease <- topTable(fit_incr_disease, coef = "disease_staterheumatoid arthritis", number = Inf, adjust.method = "BH")

gc_raw_disease    <- gc_lambda(results_raw_disease$P.Value)
gc_combat_disease <- gc_lambda(results_full_disease$P.Value)
gc_icombat_disease<- gc_lambda(results_incr_disease$P.Value)

## smoking
fit_raw_smoke <- eBayes(lmFit(mval_subset, design_full_smoke))
fit_full_smoke <- eBayes(lmFit(corrected_full_subset, design_full_smoke))
fit_incr_smoke <- eBayes(lmFit(corrected_incr_subset, design_full_smoke))

results_raw_smoke  <- topTable(fit_raw_smoke,  coef = "smoking_historyyes", number = Inf, adjust.method = "BH")
results_full_smoke <- topTable(fit_full_smoke, coef = "smoking_historyyes", number = Inf, adjust.method = "BH")
results_incr_smoke <- topTable(fit_incr_smoke, coef = "smoking_historyyes", number = Inf, adjust.method = "BH")

gc_raw_smoke    <- gc_lambda(results_raw_smoke$P.Value)
gc_combat_smoke <- gc_lambda(results_full_smoke$P.Value)
gc_icombat_smoke<- gc_lambda(results_incr_smoke$P.Value)

# nSV
mod_sva <- model.matrix(~ disease_state + age + smoking_history + gender, data = pheno_full_subset)
n_sva_raw     <- get_n_sva(mval_subset, mod_sva)
n_sva_combat  <- get_n_sva(corrected_full_subset, mod_sva)
n_sva_icombat <- get_n_sva(corrected_incr_subset, mod_sva)
