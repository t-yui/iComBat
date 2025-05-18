source("/path_to_directory/iComBat/scripts/icombat.R")
library(data.table)
library(irlba)

# 1. data loading
meth_data <- fread("/path_to_data/GSE42861_processed_methylation_matrix.txt")
probe_names <- meth_data$ID_REF
meth_matrix <- as.matrix(meth_data[, -1])
rownames(meth_matrix) <- probe_names
colnames(meth_matrix) <- colnames(meth_data)[-1]

# 2. convert beta value into M value
epsilon <- 1e-6
meth_clipped <- pmin(pmax(meth_matrix, epsilon), 1 - epsilon)
mval_matrix <- log2(meth_clipped / (1 - meth_clipped))

# 3. obtain case/control information of samples
soft_lines <- readLines("/path_to_data/GSE42861_family.soft")
sample_indices <- grep("^\\^SAMPLE = ", soft_lines)
sample_info <- data.frame(sample_id = character(), disease_state = character())

for (i in seq_along(sample_indices)) {
  start <- sample_indices[i]
  end <- if (i < length(sample_indices)) sample_indices[i + 1] - 1 else length(soft_lines)
  block <- soft_lines[start:end]
  sample_id <- sub("^\\^SAMPLE = ", "", block[1])
  disease_line <- grep("!Sample_characteristics_ch1 = disease state:", block, value = TRUE)
  disease_state <- if (length(disease_line) > 0) sub("!Sample_characteristics_ch1 = disease state: ", "", disease_line) else NA
  sample_info <- rbind(sample_info, data.frame(sample_id, disease_state, stringsAsFactors = FALSE))
}

gsm_sentrix_map <- data.frame(gsm = character(), sentrix = character(), stringsAsFactors = FALSE)
for (i in seq_along(sample_indices)) {
  start <- sample_indices[i]
  end <- if (i < length(sample_indices)) sample_indices[i + 1] - 1 else length(soft_lines)
  block <- soft_lines[start:end]

  gsm_id <- sub("^\\^SAMPLE = ", "", block[1])
  sentrix_line <- grep("_Grn.idat.gz", block, value = TRUE)
  sentrix_id <- if (length(sentrix_line) > 0) {
    sub(".*_(\\d{10}_R\\d{2}C\\d{2})_Grn\\.idat\\.gz", "\\1", sentrix_line[1])
  } else {
    NA
  }
  gsm_sentrix_map <- rbind(
    gsm_sentrix_map,
    data.frame(gsm = gsm_id, sentrix = sentrix_id, stringsAsFactors = FALSE)
  )
}
gsm_sentrix_clean <- unique(sample_info[, c("disease_state", "sentrix")])

# 4. data processing
sample_ids <- colnames(mval_matrix)
batches <- sub("_.*", "", sample_ids)
pheno_data <- data.frame(sample_id = sample_ids, batch = batches, stringsAsFactors = FALSE)
sample_info <- merge(sample_info, gsm_sentrix_map, by.x = "sample_id", by.y = "gsm")

pheno_data <- merge(pheno_data, gsm_sentrix_clean,
  by.x = "sample_id", by.y = "sentrix",
  all.x = TRUE
)
pheno_data <- subset(pheno_data, sample_id != "Pval")
pheno_data <- subset(pheno_data, disease_state == "Normal")

top5_batches <- names(table(pheno_data$batch))[1:20]
pheno_top5 <- pheno_data[pheno_data$batch %in% top5_batches, ]
pheno_top5$batch <- factor(pheno_top5$batch, levels = unique(pheno_top5$batch))

existing_batches <- levels(pheno_top5$batch)[1:10]
new_batch <- levels(pheno_top5$batch)[11:20]
existing_samples <- pheno_top5$sample_id[pheno_top5$batch %in% existing_batches]
new_samples <- pheno_top5$sample_id[pheno_top5$batch == new_batch]

dat.old <- mval_matrix[, existing_samples]
dat.new <- mval_matrix[, new_samples]
mval_top5 <- mval_matrix[, pheno_top5$sample_id]

common_cpgs <- intersect(rownames(dat.old), rownames(dat.new))
dat.old <- dat.old[common_cpgs, ]
dat.new <- dat.new[common_cpgs, ]
mval_top5 <- mval_top5[common_cpgs, ]

# 5. perform combat and icombat

# ComBat

## check if variance > 0 for all batches and filtering
batches <- unique(pheno_top5$batch)
keep_probe <- apply(mval_top5, 1, function(vals) {
  vars <- sapply(batches, function(b) {
    ix <- which(pheno_top5$batch == b)
    var(vals[ix])
  })
  all(vars > 0)
})
mval_top5 <- mval_top5[keep_probe, ]
dat.old <- dat.old[keep_probe, ]
dat.new <- dat.new[keep_probe, ]

message(sprintf(
  "number of probes after removing batches with zero variance: %d to %d",
  length(keep_probe), sum(keep_probe)
))

mod.full <- model.matrix(~1, data = data.frame(sample = pheno_top5$sample_id))
corrected_full <- ComBat(
  dat = mval_top5,
  batch = pheno_top5$batch,
  mod = NULL
)

# iComBat
batch.old <- droplevels(pheno_top5$batch[pheno_top5$sample_id %in% existing_samples])
batch.new <- droplevels(pheno_top5$batch[pheno_top5$sample_id %in% new_samples])

mod.old <- model.matrix(~1, data = data.frame(sample = existing_samples))
mod.new <- model.matrix(~1, data = data.frame(sample = new_samples))

corrected_old <- ComBat(
  dat   = dat.old,
  batch = batch.old,
  mod   = mod.old
)

fit.old <- ComBatInitialize(
  dat = dat.old,
  batch = batch.old,
  mod = mod.old,
  par.prior = TRUE, mean.only = FALSE
)

corrected_new <- ComBatIncremental(
  dat.new   = dat.new,
  new.batch = batch.new,
  fit       = fit.old,
  mod.new   = mod.new
)

corrected_incr <- cbind(corrected_old, corrected_new)

# 6. visualization
top_cpgs <- order(rowVars(mval_top5), decreasing = TRUE)[1:10000]
sample_order <- c(existing_samples, new_samples)
batch_colors <- pheno_top5$batch[match(sample_order, pheno_top5$sample_id)]

pca_raw <- prcomp_irlba(t(mval_top5[top_cpgs, sample_order]), n = 2)
pca_combat <- prcomp_irlba(t(corrected_full[top_cpgs, sample_order]), n = 2)
pca_icombat <- prcomp_irlba(t(corrected_incr[top_cpgs, sample_order]), n = 2)

par(mfrow = c(1, 3))
col_vector <- as.numeric(batch_colors)
batch_levels <- levels(batch_colors)

# (a) No Correction
plot(pca_raw$x,
  col = col_vector, pch = 19,
  xlab = "PC1", ylab = "PC2", main = "No Correction"
)
legend("topright",
  legend = paste0("Batch ", seq_along(batch_levels), ": ", batch_levels),
  col = seq_along(batch_levels),
  pch = 19, cex = 0.8
)
# (2) ComBat
plot(pca_combat$x,
  col = col_vector, pch = 19,
  xlab = "PC1", ylab = "PC2", main = "Corrected by ComBat"
)
legend("topright",
  legend = paste0("Batch ", seq_along(batch_levels), ": ", batch_levels),
  col = seq_along(batch_levels),
  pch = 19, cex = 0.8
)
# (3) iComBat
plot(pca_icombat$x,
  col = col_vector, pch = 19,
  xlab = "PC1", ylab = "PC2", main = "Corrected by iComBat"
)
legend("topright",
  legend = paste0("Batch ", seq_along(batch_levels), ": ", batch_levels),
  col = seq_along(batch_levels),
  pch = 19, cex = 0.8
)
par(mfrow = c(1, 1))
