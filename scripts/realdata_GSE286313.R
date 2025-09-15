source("icombat.R")
source("utils.R")

suppressPackageStartupMessages({
  library(sesame)
  library(sesameData)
  library(BiocParallel)
  library(data.table)
  library(GenomicRanges)
  library(R.utils)
  library(ggplot2)
  library(patchwork)
  library(sva)
  library(reshape2)
  library(RColorBrewer)
  library(methylclock)
})


# ===== DATA QC =====

pval_cutoff <- 0.01  # pOOBAH cutoff
probe_detect_min  <- 0.90  # probe cutoff
drop_sex_chr <- TRUE
genome_build <- "hg19"

decompress_idat_gz(getwd())

# load metadata from family.soft
read_family_soft_metadata <- function(soft_file) {
  if (!file.exists(soft_file)) {
    warning("family.soft file not found")
    return(NULL)
  }
  lines <- readLines(soft_file, warn = FALSE)
  sample_starts <- which(grepl("^\\^SAMPLE = ", lines))
  if (length(sample_starts) == 0) return(NULL)
  metadata_list <- list()
  for (i in seq_along(sample_starts)) {
    start <- sample_starts[i]
    end <- if (i < length(sample_starts)) sample_starts[i+1] - 1 else length(lines)
    sample_lines <- lines[start:end]

    gsm_line <- grep("^!Sample_geo_accession = ", sample_lines, value = TRUE)
    if (length(gsm_line) == 0) next
    gsm <- gsub("^!Sample_geo_accession = ", "", gsm_line[1])

    title_line <- grep("^!Sample_title = ", sample_lines, value = TRUE)
    title <- ifelse(length(title_line) > 0, 
                    gsub("^!Sample_title = ", "", title_line[1]), NA)

    desc_line <- grep("^!Sample_description = ", sample_lines, value = TRUE)
    description <- ifelse(length(desc_line) > 0, 
                          gsub("^!Sample_description = ", "", desc_line[1]), NA)

    char_lines <- grep("^!Sample_characteristics_ch1 = ", sample_lines, value = TRUE)

    sex <- NA
    age <- NA
    tissue <- NA
    
    for (char_line in char_lines) {
      content <- gsub("^!Sample_characteristics_ch1 = ", "", char_line)
      if (grepl("^Sex:", content)) {
        sex <- gsub("^Sex:\\s*", "", content)
      } else if (grepl("^age:", content)) {
        age <- as.numeric(gsub("^age:\\s*", "", content))
      } else if (grepl("^tissue:", content)) {
        tissue <- gsub("^tissue:\\s*", "", content)
      }
    }
    metadata_list[[i]] <- data.frame(
      GSM_ID = gsm,
      title = title,
      description = description,
      sex = sex,
      age = age,
      tissue = tissue,
      stringsAsFactors = FALSE
    )
  }
  metadata_df <- do.call(rbind, metadata_list)
  return(metadata_df)
}
soft_file <- "../GSE286313_family.soft"
soft_metadata <- read_family_soft_metadata(soft_file)
epicv2_samples <- soft_metadata[grepl("EPICv2", soft_metadata$description), ]
epicv2_samples$study <- gsub(",.*", "", epicv2_samples$description)
study_counts <- table(epicv2_samples$study)
epicv2_gsm_ids <- epicv2_samples$GSM_ID

find_idat_files_gse286313 <- function(base_dir = getwd(), gsm_filter = NULL) {
  idat_files <- list.files(base_dir, pattern = "(?i)\\.idat$", 
                           recursive = TRUE, full.names = TRUE)
  if (length(idat_files) == 0) stop("No *.idat files found after decompression.")
  prefixes <- unique(gsub("_(Grn|Red|Green|RED|GRN)\\.idat$", "", basename(idat_files), ignore.case = TRUE))
  if (!is.null(gsm_filter)) {
    gsm_in_file <- gsub("_.*", "", prefixes)
    prefixes <- prefixes[gsm_in_file %in% gsm_filter]
  }
  full_prefixes <- file.path(dirname(idat_files[1]), prefixes)
  has_pair <- function(pref) {
    base <- basename(pref)
    grn_exists <- any(grepl(paste0(base, "_(Grn|Green|GRN)\\.idat$"), basename(idat_files), ignore.case = TRUE))
    red_exists <- any(grepl(paste0(base, "_(Red|RED)\\.idat$"), basename(idat_files), ignore.case = TRUE))
    grn_exists && red_exists
  }
  valid_prefixes <- full_prefixes[vapply(full_prefixes, has_pair, logical(1))]
  if (length(valid_prefixes) == 0) stop("Found *.idat files but no valid Grn/Red pairs.")
  return(valid_prefixes)
}
decompress_idat_gz(getwd())
idat_prefixes <- find_idat_files_gse286313(getwd(), epicv2_gsm_ids)
gsm_ids <- gsub("_.*", "", basename(idat_prefixes))
sample_names  <- gsub("^GSM[0-9]+_", "", basename(idat_prefixes))
full_names <- basename(idat_prefixes) 
message(sprintf("Detected %d EPICv2 IDAT pairs for GSE286313.", length(idat_prefixes)))
try(sesameDataCache("EPICv2"), silent = TRUE)

# per-sample QC
qc_list <- vector("list", length(idat_prefixes))
betas_list <- vector("list", length(idat_prefixes))
names(betas_list) <- sample_names
for (i in seq_along(idat_prefixes)) {
  pref <- idat_prefixes[i]
  sname <- sample_names[i]
  gsm <- gsm_ids[i]
  message(sprintf("[%d/%d] Processing: %s", i, length(idat_prefixes), basename(pref)))
  tryCatch({
    sdf <- readIDATpair(pref, platform = "EPICv2")
    sdf <- noob(sdf)
    sdf <- dyeBiasCorrTypeINorm(sdf)
    
    # detection p-value
    sdf <- pOOBAH(sdf, pval.threshold = pval_cutoff)
    
    # QC stats
    qcd <- as.data.frame(sesameQC_calcStats(sdf))
    qcd$Sample_Name <- sname  # メタデータとマッチする名前
    qcd$GSM_ID <- gsm  # GSM番号も保存
    qcd$median_total_intensity <- medianTotalIntensity(sdf)
    
    # masked betas
    b <- getBetas(sdf, mask = TRUE)
    qcd$detection_rate <- mean(!is.na(b))
    qcd$n_probes <- length(b)
    
    qc_list[[i]] <- qcd
    betas_list[[i]] <- b
    
  }, error = function(e) {
    message(sprintf("Error processing sample %s: %s", basename(pref), e$message))
    qc_list[[i]] <- data.frame(Sample_Name = sname, 
                               GSM_ID = gsm, 
                               Error = e$message,
                               stringsAsFactors = FALSE)
    betas_list[[i]] <- NULL
  })
}

# remove error samples
valid_idx <- !sapply(betas_list, is.null)
if (sum(!valid_idx) > 0) {
  message(sprintf("Warning: %d samples failed to process.", sum(!valid_idx)))
  error_samples <- basename(idat_prefixes[!valid_idx])
  error_info <- qc_list[!valid_idx]
  error_df <- do.call(rbind, error_info)
  qc_list <- qc_list[valid_idx]
  betas_list <- betas_list[valid_idx]
  sample_names <- sample_names[valid_idx]
  gsm_ids <- gsm_ids[valid_idx]
}

## QC results
qc_df <- rbindlist(qc_list, fill = TRUE)

# data preparation
extract_sentrix_info <- function(sample_names) {
  sentrix_id <- gsub("_[Rr][0-9]{2}C[0-9]{2}$", "", sample_names)
  sentrix_pos <- gsub("^.*_([Rr][0-9]{2}C[0-9]{2})$", "\\1", sample_names)
  sentrix_pos <- toupper(sentrix_pos)
  
  data.frame(
    Sample_Name = sample_names,
    Sentrix_ID = sentrix_id,
    Sentrix_Pos = sentrix_pos,
    stringsAsFactors = FALSE
  )
}
sentrix_info <- extract_sentrix_info(qc_df$Sample_Name)
qc_df <- merge(qc_df, sentrix_info, by = "Sample_Name", all.x = TRUE)
qc_df <- merge(qc_df, epicv2_samples, by = "GSM_ID", all.x = TRUE)
fwrite(qc_df, "GSE286313_EPICv2_sesame_QC_metrics.csv")

all_probes <- unique(unlist(lapply(betas_list, names)))

# create Beta value matrix
beta_mat <- do.call(cbind, lapply(betas_list, function(v) v[all_probes]))
colnames(beta_mat) <- qc_df$Sample_Name
rownames(beta_mat) <- all_probes

# probe filtering by detection across samples 
probe_detect_rate <- rowMeans(!is.na(beta_mat))
keep_probes <- probe_detect_rate >= probe_detect_min
beta_mat <- beta_mat[keep_probes, , drop = FALSE]
message(sprintf("Probes retained after detection filtering: %d", nrow(beta_mat)))

# drop sex-chromosome probes
if (drop_sex_chr) {
  tryCatch({
    man <- sesameData_getManifestGRanges(platform = "EPICv2", genome = genome_build)
    chr <- as.character(seqnames(man))
    autosomal <- names(man)[!(chr %in% c("chrX","chrY"))]
    keep <- rownames(beta_mat) %in% autosomal
    beta_mat <- beta_mat[keep, , drop = FALSE]
    message(sprintf("Probes retained after removing sex chromosomes: %d", nrow(beta_mat)))
  }, error = function(e) {
    warning(sprintf("Error loading EPICv2 manifest: %s", e$message))
  })
}

# convert Beta to M-values
m_mat <- log2(beta_mat/(1 - beta_mat))

# saving
fwrite(as.data.table(beta_mat, keep.rownames = "Probe_ID"), "GSE286313_EPICv2_beta_values_sesame_QCed.csv")
fwrite(as.data.table(m_mat,    keep.rownames = "Probe_ID"), "GSE286313_EPICv2_m_values_sesame_QCed.csv")
fwrite(qc_df, "GSE286313_EPICv2_full_metadata.csv")
saveRDS(list(beta = beta_mat,
             m    = m_mat,
             qc   = qc_df,
             genome = genome_build),
        file = "GSE286313_EPICv2_sesame_QCed_data.rds")


# ===== ANALYSIS =====

# load data
qc_data <- readRDS("GSE286313_EPICv2_sesame_QCed_data.rds")
beta_mat_full <- qc_data$beta
m_mat_full <- qc_data$m
metadata_full <- qc_data$qc

extract_cpg_name <- function(cpg_with_suffix) {
  sub("_.*", "", cpg_with_suffix)
}
original_cpg_names <- rownames(beta_mat_full)
n_cpgs <- length(original_cpg_names)
pb <- txtProgressBar(min = 0, max = n_cpgs, style = 3)
clean_cpg_names <- character(n_cpgs)
for (i in seq_len(n_cpgs)) {
  clean_cpg_names[i] <- extract_cpg_name(original_cpg_names[i])
  if (i %% 1000 == 0) {
    setTxtProgressBar(pb, i)
  }
}
setTxtProgressBar(pb, n_cpgs)
close(pb)
if (any(dup_status)) {
  keep_indices <- !dup_status
  beta_mat_full <- beta_mat_full[keep_indices, ]
  m_mat_full <- m_mat_full[keep_indices, ]
  rownames(beta_mat_full) <- clean_cpg_names[keep_indices]
  rownames(m_mat_full) <- clean_cpg_names[keep_indices]
} else {
  rownames(beta_mat_full) <- clean_cpg_names
  rownames(m_mat_full) <- clean_cpg_names
}

# select EPIC v2 samples
epicv2_idx <- grepl("EPICv2", metadata_full$description)
beta_mat <- beta_mat_full[, epicv2_idx]
m_mat <- m_mat_full[, epicv2_idx]
metadata <- metadata_full[epicv2_idx, ]
common_samples <- intersect(colnames(beta_mat), metadata$Sample_Name)
beta_mat <- beta_mat[, common_samples]
m_mat <- m_mat[, common_samples]
metadata <- metadata[match(common_samples, metadata$Sample_Name), ]

# batch order
batch_order <- c("CALERIE", "BeCOME", "CLHNS", "VHAS")
metadata$study <- factor(metadata$study, levels = batch_order)


## ---- 2) Horvath CpGsリストの定義と照合 ----
cat("\n=== 2. Horvath CpGsの定義と照合 ===\n")

# Horvath CpGs list (without intercept)
methylclock::load_DNAm_Clocks_data()
horvath_cpgs <- as.vector(coefHorvath$CpGmarker)[2:length(coefHorvath$CpGmarker)]

horvath_cpgs_present <- intersect(horvath_cpgs, rownames(beta_mat))
cat(sprintf("Horvath CpGs total: %d\n", length(horvath_cpgs)))
cat(sprintf("Horvath CpGs present in data: %d (%.1f%%)\n", 
            length(horvath_cpgs_present), 
            100 * length(horvath_cpgs_present) / length(horvath_cpgs)))

# covariates
metadata$sex_numeric <- ifelse(metadata$sex == "F", 0, 1)

m_mat_work <- m_mat
beta_mat_work <- beta_mat
batch_data_m <- list()
batch_data_beta <- list()
batch_mod <- list()

for (batch in batch_order) {
  idx <- which(metadata$study == batch)
  batch_data_m[[batch]] <- m_mat_work[, idx, drop = FALSE]
  batch_data_beta[[batch]] <- beta_mat_work[, idx, drop = FALSE]
  mod_batch <- model.matrix(~ metadata$age[idx] + metadata$sex_numeric[idx])
  colnames(mod_batch) <- c("Intercept", "Age", "Sex")
  batch_mod[[batch]] <- mod_batch
}

anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}

calculate_horvath_age <- function(beta_values, impute_missing = TRUE) {
  horvath_cpgs_only <- coefHorvath$CpGmarker[-1]
  horvath_cpgs_present <- intersect(horvath_cpgs_only, rownames(beta_values))
  missing_cpgs <- setdiff(horvath_cpgs_only, rownames(beta_values))
  
  if (length(missing_cpgs) > 0 && impute_missing) {
    mean_beta <- mean(beta_values, na.rm = TRUE)
    missing_matrix <- matrix(mean_beta, 
                             nrow = length(missing_cpgs), 
                             ncol = ncol(beta_values))
    rownames(missing_matrix) <- missing_cpgs
    colnames(missing_matrix) <- colnames(beta_values)
    
    beta_horvath <- rbind(beta_values[horvath_cpgs_present, , drop = FALSE], 
                          missing_matrix)
    beta_horvath <- beta_horvath[horvath_cpgs_only, , drop = FALSE]
    used_cpgs <- horvath_cpgs_only
  } else {
    beta_horvath <- beta_values[horvath_cpgs_present, , drop = FALSE]
    used_cpgs <- horvath_cpgs_present
  }
  beta_df <- as.data.frame(t(beta_horvath))

  coef_indices <- which(coefHorvath$CpGmarker %in% c("(Intercept)", used_cpgs))
  coef_horvath_subset <- coefHorvath[coef_indices, ]
  
  intercept <- coef_horvath_subset$CoefficientTraining[coef_horvath_subset$CpGmarker == "(Intercept)"]
  cpg_coefs <- coef_horvath_subset$CoefficientTraining[coef_horvath_subset$CpGmarker != "(Intercept)"]
  cpg_names <- coef_horvath_subset$CpGmarker[coef_horvath_subset$CpGmarker != "(Intercept)"]
  
  n_samples <- nrow(beta_df)
  horvath_age <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    sample_values <- beta_df[i, cpg_names]
    non_na_idx <- !is.na(sample_values)
    if (sum(non_na_idx) == 0) {
      horvath_age[i] <- NA
    } else {
      valid_values <- as.numeric(sample_values[non_na_idx])
      valid_coefs <- cpg_coefs[non_na_idx]
      weighted_sum <- sum(valid_values * valid_coefs)
      linear_predictor <- weighted_sum + intercept
      horvath_age[i] <- anti.trafo(linear_predictor)
    }
  }
  na_count <- sum(is.na(horvath_age))
  if (na_count > 0) {
    cat(sprintf("%d samples failed to calculate epigenetic age\n", na_count))
  }
  return(horvath_age)
}

# ComBat
mod_all <- model.matrix(~ metadata$age + metadata$sex_numeric)
colnames(mod_all) <- c("Intercept", "Age", "Sex")

m_combat_all <- ComBat(
  dat = m_mat_work,
  batch = metadata$study,
  mod = mod_all,
  par.prior = TRUE,
  mean.only = TRUE
)
beta_combat_all <- 2^m_combat_all / (1 + 2^m_combat_all)
horvath_age_all <- calculate_horvath_age(beta_combat_all)

# Standard ComBat with incremantal addition
incremental_standard_results_m <- list()
incremental_standard_results_beta <- list()
horvath_age_incremental_standard <- matrix(NA, nrow = nrow(metadata), 
                                           ncol = length(batch_order))
colnames(horvath_age_incremental_standard) <- paste0("After_", batch_order)

for (i in 1:length(batch_order)) {
  cat(sprintf("\n  Step %d: Adding %s\n", i, batch_order[i]))
  current_batches <- batch_order[1:i]
  current_idx <- which(metadata$study %in% current_batches)
  current_data_m <- m_mat_work[, current_idx, drop = FALSE]
  current_metadata <- metadata[current_idx, ]
  current_metadata$study <- factor(as.character(current_metadata$study), 
                                   levels = current_batches)
  current_mod <- model.matrix(~ current_metadata$age + current_metadata$sex_numeric)
  colnames(current_mod) <- c("Intercept", "Age", "Sex")
  if (i == 1) {
    m_combat_current <- current_data_m
  } else {
    m_combat_current <- ComBat(
      dat = current_data_m,
      batch = current_metadata$study,
      mod = current_mod,
      par.prior = TRUE,
      mean.only = TRUE
    )
  }
  incremental_standard_results_m[[i]] <- m_combat_current
  beta_combat_current <- 2^m_combat_current / (1 + 2^m_combat_current)
  incremental_standard_results_beta[[i]] <- beta_combat_current
  horvath_age_current <- calculate_horvath_age(beta_combat_current)
  horvath_age_incremental_standard[current_idx, i] <- horvath_age_current
}

# iComBat with incremantal addition
incremental_combat_results_m <- list()
incremental_combat_results_beta <- list()
horvath_age_incremental_combat <- matrix(NA, nrow = nrow(metadata), 
                                         ncol = length(batch_order))
colnames(horvath_age_incremental_combat) <- paste0("After_", batch_order)
combat_fit <- NULL
for (i in 1:length(batch_order)) {
  cat(sprintf("\n  Step %d: Adding %s\n", i, batch_order[i]))
  current_batch <- batch_order[i]
  current_idx <- which(metadata$study == current_batch)
  if (i == 1) {
    m_combat_current <- batch_data_m[[current_batch]]
    incremental_combat_results_m[[i]] <- m_combat_current
    beta_combat_current <- 2^m_combat_current / (1 + 2^m_combat_current)
    incremental_combat_results_beta[[i]] <- beta_combat_current
    horvath_age_current <- calculate_horvath_age(beta_combat_current)
    horvath_age_incremental_combat[current_idx, i] <- horvath_age_current
    first_metadata <- metadata[current_idx, ]
    first_metadata$study <- factor(as.character(first_metadata$study), 
                                   levels = current_batch)
    first_mod <- model.matrix(~ first_metadata$age + first_metadata$sex_numeric)
    colnames(first_mod) <- c("Intercept", "Age", "Sex")
    
    combat_fit <- ComBatInitialize(
      dat = m_combat_current,
      batch = first_metadata$study,
      mod = first_mod,
      mean.only = TRUE
    )
    
  } else {
    prev_batches <- batch_order[1:(i-1)]
    prev_idx <- which(metadata$study %in% prev_batches)
    new_data <- batch_data_m[[current_batch]]
    new_mod <- batch_mod[[current_batch]]
    
    m_combat_new <- ComBatIncremental(
      dat.new = new_data,
      new.batch = current_batch,
      fit = combat_fit,
      mod.new = new_mod
    )
    rownames(m_combat_new) <- rownames(new_data)

    all_idx <- c(prev_idx, current_idx)
    if (i == 2) {
      m_combat_all_current <- cbind(incremental_combat_results_m[[1]], m_combat_new)
    } else {
      m_combat_all_current <- cbind(incremental_combat_results_m[[i-1]], m_combat_new)
    }
    incremental_combat_results_m[[i]] <- m_combat_all_current
    beta_combat_new <- 2^m_combat_new / (1 + 2^m_combat_new)
    if (i == 2) {
      beta_combat_all_current <- cbind(incremental_combat_results_beta[[1]], beta_combat_new)
    } else {
      beta_combat_all_current <- cbind(incremental_combat_results_beta[[i-1]], beta_combat_new)
    }
    incremental_combat_results_beta[[i]] <- beta_combat_all_current
    horvath_age_incremental_combat[prev_idx, i] <- horvath_age_incremental_combat[prev_idx, i-1]
    horvath_age_new <- calculate_horvath_age(beta_combat_new)
    horvath_age_incremental_combat[current_idx, i] <- horvath_age_new
  }
}

# mean change of epigenetic age in existing batches
mean_changes_existing <- data.frame(
  Step = character(),
  Method = character(),
  Mean_Change_Existing = numeric(),
  SD_Change_Existing = numeric(),
  Max_Change_Existing = numeric(),
  stringsAsFactors = FALSE
)

for (i in 2:length(batch_order)) {
  existing_batches <- batch_order[1:(i-1)]
  existing_idx <- which(metadata$study %in% existing_batches)
  
  if (length(existing_idx) > 0) {
    change_standard <- horvath_age_incremental_standard[existing_idx, i] - 
      horvath_age_incremental_standard[existing_idx, i-1]
    
    mean_changes_existing <- rbind(mean_changes_existing,
                                   data.frame(
                                     Step = paste0("Add_", batch_order[i]),
                                     Method = "Incremental_Standard",
                                     Mean_Change_Existing = mean(change_standard, na.rm = TRUE),
                                     SD_Change_Existing = sd(change_standard, na.rm = TRUE),
                                     Max_Change_Existing = max(abs(change_standard), na.rm = TRUE)
                                   )
    )
    
    change_combat <- horvath_age_incremental_combat[existing_idx, i] - 
      horvath_age_incremental_combat[existing_idx, i-1]
    
    mean_changes_existing <- rbind(mean_changes_existing,
                                   data.frame(
                                     Step = paste0("Add_", batch_order[i]),
                                     Method = "Incremental_ComBat",
                                     Mean_Change_Existing = mean(change_combat, na.rm = TRUE),
                                     SD_Change_Existing = sd(change_combat, na.rm = TRUE),
                                     Max_Change_Existing = max(abs(change_combat), na.rm = TRUE)
                                   )
    )
  }
}

change_data <- data.frame()
for (i in 2:length(batch_order)) {
  existing_batches <- batch_order[1:(i-1)]
  existing_idx <- which(metadata$study %in% existing_batches)
  
  if (length(existing_idx) > 0) {
    change_standard <- horvath_age_incremental_standard[existing_idx, i] - 
      horvath_age_incremental_standard[existing_idx, i-1]
    
    change_data <- rbind(change_data, data.frame(
      Step = paste0("Add ", batch_order[i]),
      Method = "Incremental Correction by ComBat",
      Change = change_standard,
      Sample_Study = metadata$study[existing_idx]
    ))
    
    change_combat <- horvath_age_incremental_combat[existing_idx, i] - 
      horvath_age_incremental_combat[existing_idx, i-1]
    
    change_data <- rbind(change_data, data.frame(
      Step = paste0("Add ", batch_order[i]),
      Method = "Incremental Correction by iComBat",
      Change = change_combat,
      Sample_Study = metadata$study[existing_idx]
    ))
  }
}
p_change <- ggplot(change_data, aes(x = Step, y = Change, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Incremental Correction by ComBat" = "#FF6B6B", 
                               "Incremental Correction by iComBat" = "#4ECDC4")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Step",
    y = "Change in Epigenetic Age (years)",
    fill = "Method"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) 
print(p_change)

# compare with all-batch ComBat
final_comparison <- data.frame(
  Sample_Name = metadata$Sample_Name,
  Study = metadata$study,
  Horvath_All_Batch = horvath_age_all,
  Horvath_Inc_Standard = horvath_age_incremental_standard[, length(batch_order)],
  Horvath_Inc_Combat = horvath_age_incremental_combat[, length(batch_order)],
  Chronological_Age = metadata$age
)

# PCA
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
n_top_cpgs <- 10000
top_cpgs <- select_top_variable_cpgs(m_mat, n_cpgs = n_top_cpgs)
m_mat_pca <- m_mat[top_cpgs, ]

impute_median <- function(mat) {
  mat_imputed <- mat
  for (i in 1:nrow(mat)) {
    row_data <- mat[i, ]
    na_idx <- is.na(row_data)
    if (any(na_idx)) {
      row_median <- median(row_data, na.rm = TRUE)
      if (!is.na(row_median)) {
        mat_imputed[i, na_idx] <- row_median
      } else {
        mat_imputed[i, na_idx] <- median(mat, na.rm = TRUE)
      }
    }
  }
  return(mat_imputed)
}
m_mat_pca <- impute_median(m_mat_pca)

m_combat_all_pca <- m_combat_all[top_cpgs, ]
m_combat_all_pca <- impute_median(m_combat_all_pca)

m_inc_combat_final <- incremental_combat_results_m[[length(batch_order)]]
m_inc_combat_pca <- m_inc_combat_final[top_cpgs, ]
m_inc_combat_pca <- impute_median(m_inc_combat_pca)

pca_raw <- prcomp(t(m_mat_pca), scale. = FALSE, center = TRUE)
pca_combat_all <- prcomp(t(m_combat_all_pca), scale. = FALSE, center = TRUE)
pca_inc_combat <- prcomp(t(m_inc_combat_pca), scale. = FALSE, center = TRUE)

# 分散説明率
var_raw <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2)) * 100
var_combat_all <- (pca_combat_all$sdev^2 / sum(pca_combat_all$sdev^2)) * 100
var_inc_combat <- (pca_inc_combat$sdev^2 / sum(pca_inc_combat$sdev^2)) * 100

get_color_palette <- function(n) {
  base_colors <- c(
    "#377EB8",  # Blue
    "#E41A1C",  # Red
    "#4CAF50",  # Green
    "#FF9800"   # Orange
  )
  return(base_colors[1:n])
}
n_batches <- length(batch_order)
batch_colors <- get_color_palette(n_batches)
batch_palette <- stats::setNames(batch_colors, batch_order)
batch_breaks <- batch_order
new_batches <- c("BeCOME","CLHNS","VHAS")

pca_data_raw <- data.frame(
  PC1 = pca_raw$x[,1], 
  PC2 = pca_raw$x[,2], 
  PC3 = pca_raw$x[,3],
  Sex = metadata$sex,
  Age = metadata$age,
  Batch = metadata$study,
  IsNew = metadata$study %in% new_batches
)
pca_data_combat_all <- data.frame(
  PC1 = pca_combat_all$x[,1], 
  PC2 = pca_combat_all$x[,2], 
  PC3 = pca_combat_all$x[,3],
  Sex = metadata$sex,
  Age = metadata$age,
  Batch = metadata$study,
  IsNew = metadata$study %in% new_batches
)
pca_data_inc_combat <- data.frame(
  PC1 = pca_inc_combat$x[,1], 
  PC2 = pca_inc_combat$x[,2], 
  PC3 = pca_inc_combat$x[,3],
  Sex = metadata$sex,
  Age = metadata$age,
  Batch = metadata$study,
  IsNew = metadata$study %in% new_batches
)

p1_raw <- ggplot(pca_data_raw, aes(PC1, PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), 
       y = sprintf("PC2 (%.1f%%)", var_raw[2]), 
       title = "Raw") +
  theme_light_plot(12) + 
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p1_combat_all <- ggplot(pca_data_combat_all, aes(PC1, PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat_all[1]), 
       y = sprintf("PC2 (%.1f%%)", var_combat_all[2]), 
       title = "All-Batch ComBat") +
  theme_light_plot(12) + 
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p1_inc_combat <- ggplot(pca_data_inc_combat, aes(PC1, PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_inc_combat[1]), 
       y = sprintf("PC2 (%.1f%%)", var_inc_combat[2]), 
       title = "iComBat") +
  theme_light_plot(12) + 
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p13_raw <- ggplot(pca_data_raw, aes(PC1, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), 
       y = sprintf("PC3 (%.1f%%)", var_raw[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p13_combat_all <- ggplot(pca_data_combat_all, aes(PC1, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat_all[1]), 
       y = sprintf("PC3 (%.1f%%)", var_combat_all[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p13_inc_combat <- ggplot(pca_data_inc_combat, aes(PC1, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", var_inc_combat[1]), 
       y = sprintf("PC3 (%.1f%%)", var_inc_combat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p23_raw <- ggplot(pca_data_raw, aes(PC2, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC2 (%.1f%%)", var_raw[2]), 
       y = sprintf("PC3 (%.1f%%)", var_raw[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

p23_combat_all <- ggplot(pca_data_combat_all, aes(PC2, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = "Batch (Cohort)",
                     labels = c(
                       "CALERIE" = "CALERIE (Existing)",
                       "BeCOME"  = "BeCOME (New)",
                       "CLHNS"   = "CLHNS (New)",
                       "VHAS"    = "VHAS (New)"
                     )) +
  labs(x = sprintf("PC2 (%.1f%%)", var_combat_all[2]), 
       y = sprintf("PC3 (%.1f%%)", var_combat_all[3])) +
  theme_light_plot(12) + theme(legend.position = "bottom")

p23_inc_combat <- ggplot(pca_data_inc_combat, aes(PC2, PC3, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = NULL) +
  labs(x = sprintf("PC2 (%.1f%%)", var_inc_combat[2]), 
       y = sprintf("PC3 (%.1f%%)", var_inc_combat[3])) +
  theme_light_plot(12) + theme(legend.position = "none")

pca_combined <- (p1_raw | p1_combat_all | p1_inc_combat) / 
  (p13_raw | p13_combat_all | p13_inc_combat) / 
  (p23_raw | p23_combat_all | p23_inc_combat)
print(pca_combined)

res_raw <- check_PC_covariate(pca_data_raw)
res_combat <- check_PC_covariate(pca_data_combat_all)
res_icombat <- check_PC_covariate(pca_data_inc_combat)

res_raw$PC1
res_combat$PC1
res_icombat$PC1

# epigenetic age comparison: All batch combat vs incremental combat
cor_value <- cor(final_comparison$Horvath_All_Batch, 
                 final_comparison$Horvath_Inc_Combat, 
                 use = "complete.obs")
p_comparison <- ggplot(final_comparison, aes(x = Horvath_All_Batch, y = Horvath_Inc_Combat, 
                                   color = Study)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_manual(values = batch_palette, breaks = batch_breaks,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = "Batch") +
  annotate("text", 
           x = min(final_comparison$Horvath_All_Batch, na.rm = TRUE) + 
             0.05 * diff(range(final_comparison$Horvath_All_Batch, na.rm = TRUE)),
           y = max(final_comparison$Horvath_Inc_Combat, na.rm = TRUE) - 
             0.05 * diff(range(final_comparison$Horvath_Inc_Combat, na.rm = TRUE)),
           label = paste0("r = ", sprintf("%.3f", cor_value)),
           size = 5,
           hjust = 0,
           fontface = "bold") +
  labs(
    title = NULL,
    x = "Horvath Epigenetic Age (All-Batch ComBat)",
    y = "Horvath Epigenetic Age (iComBat)",
    color = "Study"
  ) +
  theme_bw() + theme(panel.grid = element_blank()) +  
  coord_fixed()
print(p_ea_comparison)

# nSV
m_mat_sva <- impute_median(m_mat)
m_combat_all_sva <- impute_median(m_combat_all)
m_inc_combat_final_sva <- impute_median(incremental_combat_results_m[[length(batch_order)]])

mod_sva_full <- model.matrix(~ metadata$age + metadata$sex_numeric)
colnames(mod_sva_full) <- c("Intercept", "Age", "Sex")

n_sva_raw_original <- get_n_sva(m_mat_sva, mod_sva_full)
n_sva_combat_all_original <- get_n_sva(m_combat_all_sva, mod_sva_full)
n_sva_inc_combat_original <- get_n_sva(m_inc_combat_final_sva, mod_sva_full)

# sample to sample density plot
combat_all_values <- as.vector(m_combat_all)
icombat_values <- as.vector(m_inc_combat_final)
cor_value <- cor(combat_all_values, icombat_values, use = "complete.obs")
cat(sprintf("\nCorrelation between All-Batch ComBat and iComBat: r = %.4f\n", cor_value))

n_points <- length(combat_all_values)
if (n_points > 100000) {
  set.seed(1)
  idx <- sample.int(n_points, 100000)
  plot_combat_all <- combat_all_values[idx]
  plot_icombat <- icombat_values[idx]
} else {
  plot_combat_all <- combat_all_values
  plot_icombat <- icombat_values
}

density_plot <- ggplot(data.frame(ComBat = plot_combat_all, iComBat = plot_icombat), 
                       aes(ComBat, iComBat)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
  scale_fill_viridis(name = "Density", option = "plasma") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "All-Batch ComBat corrected values", 
       y = "iComBat corrected values",
       title = sprintf("r = %.3f", cor_value)) +
  theme_light_plot() + theme(aspect.ratio = 1) + coord_fixed()
print(density_plot)
