source("icombat.R")
source("utils.R")

suppressPackageStartupMessages({
  library(sesame)
  library(sesameData)
  library(BiocParallel)
  library(data.table)
  library(GenomicRanges)
  library(R.utils)
  library(data.table)
  library(irlba)
  library(limma)
  library(ggplot2)
  library(sva)
  library(viridis)
  library(patchwork)
  library(hexbin)
})


# ===== DATA QC =====

pval_cutoff <- 0.01  # pOOBAH cutoff
probe_detect_min  <- 0.90  # probe cutoff
drop_sex_chr <- TRUE
genome_build <- "hg19"

platform_450K_samples <- c(
  "100973330057_R03C01", "100973330105_R03C02", "100973330105_R06C01", 
  "3999002065_R01C01", "3999002065_R05C01", "3999078138_R06C02",
  "3999142124_R02C01", "9406922007_R01C02", "9406922007_R03C01",
  "9406922007_R03C02", "9406922009_R02C01", "9611518029_R01C02",
  "9611518029_R03C02", "9611518029_R04C01", "9611518029_R04C02",
  "9611518029_R06C01", "9611519086_R01C02", "9701417056_R01C02",
  "9701417056_R02C02", "9701417056_R03C02", "9701417056_R04C02", "9701417098_R01C02")

# find IDAT files and create prefix mapping
find_idat_files_gse224218 <- function(base_dir = getwd()) {
  idat_files <- list.files(base_dir, pattern = "(?i)\\.idat$", 
                           recursive = TRUE, full.names = TRUE)
  if (length(idat_files) == 0) stop("No *.idat files found after decompression.")
  prefixes <- unique(gsub("_(Grn|Red|Green|RED|GRN)\\.idat$", "", basename(idat_files), ignore.case = TRUE))
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
idat_prefixes <- find_idat_files_gse224218(getwd())
gsm_ids <- gsub("_.*", "", basename(idat_prefixes))
sample_names  <- gsub("^GSM[0-9]+_", "", basename(idat_prefixes))
full_names <- basename(idat_prefixes)
message(sprintf("Detected %d IDAT pairs for GSE224218.", length(idat_prefixes)))
try(sesameDataCache("HM450"), silent = TRUE)
try(sesameDataCache("EPIC"), silent = TRUE)

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
    sdf <- readIDATpair(pref)
    if (sname %in% platform_450K_samples) {
      platform <- "450K"
    } else {
      platform <- "EPIC"
    }
    message(sprintf("  Platform: %s", platform))
    
    ## normalize + dye-bias correction
    sdf <- noob(sdf)
    sdf <- dyeBiasCorrTypeINorm(sdf)
    
    ## detection p-value
    sdf <- pOOBAH(sdf, pval.threshold = pval_cutoff)
    
    ## QC stats
    qcd <- as.data.frame(sesameQC_calcStats(sdf))
    qcd$Sample_Name <- sname 
    qcd$GSM_ID <- gsm
    qcd$Platform <- platform
    qcd$median_total_intensity <- medianTotalIntensity(sdf)
    
    ## masked betas
    b <- getBetas(sdf, mask = TRUE)
    qcd$detection_rate <- mean(!is.na(b))
    qcd$n_probes <- length(b)
    
    qc_list[[i]] <- qcd
    betas_list[[i]] <- b
    
  }, error = function(e) {
    message(sprintf("Error processing sample %s: %s", basename(pref), e$message))
    qc_list[[i]] <- data.frame(Sample_Name = sname, 
                               GSM_ID = gsm, 
                               Platform = NA, 
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
}

# QC results
qc_df <- rbindlist(qc_list, fill = TRUE)

# make Beta matrix
probes_by_platform <- list()
for (plat in unique(qc_df$Platform[!is.na(qc_df$Platform)])) {
  plat_samples <- qc_df$Sample_Name[qc_df$Platform == plat]
  plat_betas <- betas_list[plat_samples]
  plat_probes <- unique(unlist(lapply(plat_betas, names)))
  probes_by_platform[[plat]] <- plat_probes
  message(sprintf("%s platform: %d probes detected", plat, length(plat_probes)))
}

# find common probes among platforms
if (length(probes_by_platform) > 1) {
  common_probes <- Reduce(intersect, probes_by_platform)
} else {
  common_probes <- probes_by_platform[[1]]
}
message(sprintf("Number of common probes across platforms: %d", length(common_probes)))

# create Beta value matrix
beta_mat <- do.call(cbind, lapply(betas_list, function(v) v[common_probes]))
colnames(beta_mat) <- qc_df$Sample_Name
rownames(beta_mat) <- common_probes

# probe filtering by detection across samples
probe_detect_rate <- rowMeans(!is.na(beta_mat))
keep_probes <- probe_detect_rate >= probe_detect_min
beta_mat <- beta_mat[keep_probes, , drop = FALSE]
message(sprintf("Probes retained after detection filtering: %d", nrow(beta_mat)))

# drop sex-chromosome probes
if (drop_sex_chr) {
  platforms_used <- unique(qc_df$Platform[!is.na(qc_df$Platform)])
  message(sprintf("Platforms in dataset: %s", paste(platforms_used, collapse = ", ")))
  autosomal_probes_all <- c()
  for (plat in platforms_used) {
    tryCatch({
      if (plat == "450K") {
        man <- sesameData_getManifestGRanges(platform = "HM450", genome = genome_build)
      } else if (plat == "EPIC") {
        man <- sesameData_getManifestGRanges(platform = "EPIC", genome = genome_build)
      } else {
        warning(sprintf("Unknown platform: %s. Skipping.", plat))
        next
      }
      chr <- as.character(seqnames(man))
      autosomal <- names(man)[!(chr %in% c("chrX","chrY"))]
      if (length(autosomal_probes_all) == 0) {
        autosomal_probes_all <- autosomal
      } else {
        autosomal_probes_all <- intersect(autosomal_probes_all, autosomal)
      }
    }, error = function(e) {
      warning(sprintf("Error loading manifest for %s: %s", plat, e$message))
    })
  }
  if (length(autosomal_probes_all) > 0) {
    keep <- rownames(beta_mat) %in% autosomal_probes_all
    beta_mat <- beta_mat[keep, , drop = FALSE]
    message(sprintf("Probes retained after removing sex chromosomes: %d", nrow(beta_mat)))
  }
}

# convert Beta to M-values
m_mat <- log2(beta_mat/(1 - beta_mat))

# saving
fwrite(as.data.table(beta_mat, keep.rownames = "Probe_ID"), "GSE224218_beta_values_sesame_QCed.csv")
fwrite(as.data.table(m_mat,    keep.rownames = "Probe_ID"), "GSE224218_m_values_sesame_QCed.csv")
saveRDS(list(beta = beta_mat,
             m    = m_mat,
             qc   = qc_df,
             sample_mapping = sample_mapping,
             common_probes = common_probes,
             platforms_used = unique(qc_df$Platform[!is.na(qc_df$Platform)]),
             genome = genome_build),
        file = "GSE224218_sesame_QCed_data.rds")

# load metadata from family.soft
read_family_soft_metadata <- function(path) {
  if (!file.exists(path)) stop("Soft file is not found: ", path)
  lines <- readLines(path, warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  meta_list <- list()
  
  for (i in seq_along(sample_starts)) {
    start_line <- sample_starts[i]
    end_line <- if (i < length(sample_starts)) sample_starts[i + 1] - 1 else length(lines)
    sample_lines <- lines[start_line:end_line]
    si <- list(
      Sample = NA_character_, GSM = NA_character_, platform = NA_character_,
      institute = NA_character_, sex = NA_character_, age = NA_integer_,
      incoming_diagnosis = NA_character_, tumor_location = NA_character_,
      tert_mutation_status = NA_character_, chromosome_6_status = NA_character_,
      progression_free_survival_months = NA_real_, progression_0_no_1_yes = NA_integer_,
      overall_survival_months = NA_real_, death_0_no_1_yes = NA_integer_,
      material_type = NA_character_, Sentrix_ID = NA_character_
    )
    
    si$GSM <- sub("^\\^SAMPLE = ", "", sample_lines[1])
    title_line <- grep("^!Sample_title = ", sample_lines, value = TRUE)[1]
    if (!is.na(title_line)) si$Sample <- sub("^!Sample_title = ", "", title_line)
    
    platform_line <- grep("^!Sample_platform_id = ", sample_lines, value = TRUE)[1]
    if (!is.na(platform_line)) {
      platform_id <- sub("^!Sample_platform_id = ", "", platform_line)
      si$platform <- ifelse(platform_id == "GPL21145", "EPIC", "450K")
    }
    
    institute_line <- grep("^!Sample_contact_institute = ", sample_lines, value = TRUE)[1]
    if (!is.na(institute_line)) si$institute <- sub("^!Sample_contact_institute = ", "", institute_line)
    
    char_lines <- grep("^!Sample_characteristics_ch1 = ", sample_lines, value = TRUE)
    for (cl in char_lines) {
      value <- sub("^!Sample_characteristics_ch1 = ", "", cl)
      if (grepl("^Sex:", value)) {
        si$sex <- sub("^Sex: ", "", value)
      } else if (grepl("^age:", value)) {
        si$age <- as.integer(sub("^age: ", "", value))
      } else if (grepl("^incoming diagnosis:", value)) {
        si$incoming_diagnosis <- sub("^incoming diagnosis: ", "", value)
      } else if (grepl("^tumor location:", value)) {
        si$tumor_location <- sub("^tumor location: ", "", value)
      } else if (grepl("^tert mutation status:", value)) {
        si$tert_mutation_status <- sub("^tert mutation status: ", "", value)
      } else if (grepl("^chromosome 6 status:", value)) {
        si$chromosome_6_status <- sub("^chromosome 6 status: ", "", value)
      } else if (grepl("^progression-free survival \\(months\\):", value)) {
        si$progression_free_survival_months <- as.numeric(sub("^progression-free survival \\(months\\): ", "", value))
      } else if (grepl("^progression \\(0=no, 1=yes\\):", value)) {
        si$progression_0_no_1_yes <- as.integer(sub("^progression \\(0=no, 1=yes\\): ", "", value))
      } else if (grepl("^overall survival \\(months\\):", value)) {
        si$overall_survival_months <- as.numeric(sub("^overall survival \\(months\\): ", "", value))
      } else if (grepl("^death \\(0=no, 1=yes\\):", value)) {
        si$death_0_no_1_yes <- as.integer(sub("^death \\(0=no, 1=yes\\): ", "", value))
      } else if (grepl("^material type:", value)) {
        si$material_type <- sub("^material type: ", "", value)
      }
    }
    if (!is.na(si$Sample)) si$Sentrix_ID <- sub("_.*", "", si$Sample)
    meta_list[[i]] <- si
    
    if (i %% 10 == 0) cat("Processed:", i, "/", length(sample_starts), "(", si$GSM, ":", si$platform, ")\n")
  }
  do.call(rbind, lapply(meta_list, function(x) data.frame(x, stringsAsFactors = FALSE)))
}


# ===== ANALYSIS =====

# load data
qc_data <- readRDS("GSE224218_sesame_QCed_data.rds")
m_mat   <- qc_data$m
qc_df   <- qc_data$qc
meta <- read_family_soft_metadata("GSE224218_family.soft")

# merge samples and metadata
common_samples <- intersect(meta$Sample, colnames(m_mat))
meta_clean <- meta[match(common_samples, meta$Sample), ]
m_mat      <- m_mat[, common_samples]

# outcome
complete_outcome_idx <- which(!is.na(meta_clean$progression_0_no_1_yes) &
                                !is.na(meta_clean$death_0_no_1_yes) &
                                !is.na(meta_clean$tumor_location))
meta_complete <- meta_clean[complete_outcome_idx, ]
m_mat_complete <- m_mat[, meta_complete$Sample]

# covariates
meta_complete$age_numeric <- as.numeric(meta_complete$age)
age_median <- median(meta_complete$age_numeric, na.rm = TRUE)
meta_complete$age_numeric[is.na(meta_complete$age_numeric)] <- age_median
meta_complete$sex_numeric <- ifelse(meta_complete$sex == "Male", 1,
                                    ifelse(meta_complete$sex == "Female", 0, NA_integer_))
meta_complete$tumor_location_numeric <- ifelse(meta_complete$tumor_location == "infratentorial", 1, 0)
covar_cols <- c("progression_0_no_1_yes", "death_0_no_1_yes",
                "tumor_location_numeric", "age_numeric", "sex_numeric")
keep_cb <- complete.cases(meta_complete[, covar_cols])
if (!all(keep_cb)) {
  drop_n <- sum(!keep_cb)
}
meta_complete <- meta_complete[keep_cb, ]
m_mat_complete <- m_mat_complete[, meta_complete$Sample]

# preparation for batch effect correction
epic_idx <- which(meta_complete$platform == "EPIC")
k450_idx <- which(meta_complete$platform == "450K")
m_mat_epic <- m_mat_complete[, epic_idx, drop = FALSE]
m_mat_450k <- m_mat_complete[, k450_idx, drop = FALSE]
meta_epic  <- meta_complete[epic_idx, ]
meta_450k  <- meta_complete[k450_idx, ]

mod_all  <- model.matrix(~ progression_0_no_1_yes + death_0_no_1_yes +
                           tumor_location_numeric + age_numeric + sex_numeric,
                         data = meta_complete)
mod_epic <- model.matrix(~ progression_0_no_1_yes + death_0_no_1_yes +
                           tumor_location_numeric + age_numeric + sex_numeric,
                         data = meta_epic)
mod_450k <- model.matrix(~ progression_0_no_1_yes + death_0_no_1_yes +
                           tumor_location_numeric + age_numeric + sex_numeric,
                         data = meta_450k)
mod_450k_before_epic <- model.matrix(~ progression_0_no_1_yes + 
                                       tumor_location_numeric + age_numeric + sex_numeric,
                                     data = meta_450k)
mod_epic_after_450k <- model.matrix(~ progression_0_no_1_yes + 
                                      tumor_location_numeric + age_numeric + sex_numeric,
                                    data = meta_epic)

# ComBat correction
batch_all <- meta_complete$platform
set.seed(1)
corrected_combat <- ComBat(
  dat = m_mat_complete,
  batch = batch_all,
  mod = mod_all,
  par.prior = TRUE,
  mean.only = FALSE
)

# iComBat: EPIC to 450K
batch_epic <- rep("EPIC", ncol(m_mat_epic))
fit_epic <- ComBatInitialize(
  dat = m_mat_epic,
  batch = batch_epic,
  mod   = mod_epic,
  mean.only = FALSE
)
m_mat_450k_corrected_v1 <- ComBatIncremental(
  dat.new   = m_mat_450k,
  new.batch = "450K",
  fit       = fit_epic,
  mod.new   = mod_450k
)
m_mat_icombat_v1 <- cbind(m_mat_epic, m_mat_450k_corrected_v1)
m_mat_icombat_v1 <- m_mat_icombat_v1[, meta_complete$Sample]

# iComBat: 450K to EPIC
batch_450k <- rep("450K", ncol(m_mat_450k))
fit_450k <- ComBatInitialize(
  dat = m_mat_450k,
  batch = batch_450k,
  mod   = mod_450k_before_epic,
  mean.only = TRUE
)
m_mat_epic_corrected_v2 <- ComBatIncremental(
  dat.new   = m_mat_epic,
  new.batch = "EPIC",
  fit       = fit_450k,
  mod.new   = mod_epic_after_450k
)
m_mat_icombat_v2 <- cbind(m_mat_450k, m_mat_epic_corrected_v2)
m_mat_icombat_v2 <- m_mat_icombat_v2[, meta_complete$Sample]

# PCA
n_top_cpgs <- min(10000, nrow(m_mat_complete))
cpg_var <- apply(m_mat_complete, 1, var)
cpg_var <- cpg_var[!is.na(cpg_var) & cpg_var > 0]
top_cpgs <- names(sort(cpg_var, decreasing = TRUE)[1:n_top_cpgs])

pca_raw        <- prcomp(t(m_mat_complete[top_cpgs, ]),        scale. = TRUE, center = TRUE)
pca_combat     <- prcomp(t(corrected_combat[top_cpgs, ]),      scale. = TRUE, center = TRUE)
pca_icombat_v1 <- prcomp(t(m_mat_icombat_v1[top_cpgs, ]),      scale. = TRUE, center = TRUE)
pca_icombat_v2 <- prcomp(t(m_mat_icombat_v2[top_cpgs, ]),      scale. = TRUE, center = TRUE)

var_raw        <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2)) * 100
var_combat     <- (pca_combat$sdev^2 / sum(pca_combat$sdev^2)) * 100
var_icombat_v1 <- (pca_icombat_v1$sdev^2 / sum(pca_icombat_v1$sdev^2)) * 100
var_icombat_v2 <- (pca_icombat_v2$sdev^2 / sum(pca_icombat_v2$sdev^2)) * 100

create_pca_df <- function(pca_result, meta_data) {
  data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], PC3 = pca_result$x[,3],
             Platform = meta_data$platform,
             Batch = meta_data$platform,
             Sex = meta_data$sex,
             Age = meta_data$age)
}
pca_df_raw        <- create_pca_df(pca_raw, meta_complete)
pca_df_combat     <- create_pca_df(pca_combat, meta_complete)
pca_df_icombat_v1 <- create_pca_df(pca_icombat_v1, meta_complete)
pca_df_icombat_v2 <- create_pca_df(pca_icombat_v2, meta_complete)

platform_colors <- c("450K" = "#E41A1C", "EPIC" = "#377EB8")

p1_12 <- ggplot(pca_df_raw, aes(PC1, PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), y = sprintf("PC2 (%.1f%%)", var_raw[2]),
       title = "Raw") + 
  theme_light_plot() +
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p2_12 <- ggplot(pca_df_combat, aes(PC1, PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat[1]), y = sprintf("PC2 (%.1f%%)", var_combat[2]),
       title = "ComBat") + 
  theme_light_plot() +
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p3_12 <- ggplot(pca_df_icombat_v1, aes(PC1, PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat_v1[1]), y = sprintf("PC2 (%.1f%%)", var_icombat_v1[2]),
       title = "iComBat: EPIC→450K") + 
  theme_light_plot() +
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p4_12 <- ggplot(pca_df_icombat_v2, aes(PC1, PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat_v2[1]), y = sprintf("PC2 (%.1f%%)", var_icombat_v2[2]),
       title = "iComBat: 450K→EPIC") + 
  theme_light_plot() +
  theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "bold"))

p1_13 <- ggplot(pca_df_raw, aes(PC1, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_raw[1]), y = sprintf("PC3 (%.1f%%)", var_raw[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p2_13 <- ggplot(pca_df_combat, aes(PC1, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_combat[1]), y = sprintf("PC3 (%.1f%%)", var_combat[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p3_13 <- ggplot(pca_df_icombat_v1, aes(PC1, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat_v1[1]), y = sprintf("PC3 (%.1f%%)", var_icombat_v1[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p4_13 <- ggplot(pca_df_icombat_v2, aes(PC1, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_icombat_v2[1]), y = sprintf("PC3 (%.1f%%)", var_icombat_v2[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p1_23 <- ggplot(pca_df_raw, aes(PC2, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC2 (%.1f%%)", var_raw[2]), y = sprintf("PC3 (%.1f%%)", var_raw[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p2_23 <- ggplot(pca_df_combat, aes(PC2, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values =platform_colors,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
                     name = "Batch (Platform)") +
  labs(x = sprintf("PC2 (%.1f%%)", var_combat[2]), y = sprintf("PC3 (%.1f%%)", var_combat[3])) + 
  theme_light_plot() +
  theme(legend.position = "bottom")

p3_23 <- ggplot(pca_df_icombat_v1, aes(PC2, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC2 (%.1f%%)", var_icombat_v1[2]), y = sprintf("PC3 (%.1f%%)", var_icombat_v1[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

p4_23 <- ggplot(pca_df_icombat_v2, aes(PC2, PC3, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = platform_colors) +
  labs(x = sprintf("PC2 (%.1f%%)", var_icombat_v2[2]), y = sprintf("PC3 (%.1f%%)", var_icombat_v2[3])) + 
  theme_light_plot() +
  theme(legend.position = "none")

pca_combined <- (p1_12 | p2_12 | p3_12 | p4_12) / 
  (p1_13 | p2_13 | p3_13 | p4_13) / 
  (p1_23 | p2_23 | p3_23 | p4_23) 
print(pca_combined)

res_raw <- check_PC_covariate_association(pca_df_raw)
res_combat <- check_PC_covariate_association(pca_df_combat)
res_icombat_v1 <- check_PC_covariate_association(pca_df_icombat_v1)
res_icombat_v2 <- check_PC_covariate_association(pca_df_icombat_v2)

res_raw$PC1
res_combat$PC1
res_icombat_v1$PC1
res_icombat_v2$PC1

# EWAS analysis
## Progression
design_prog <- model.matrix(~ progression_0_no_1_yes + age_numeric + sex_numeric, data = meta_complete)
fit_raw_prog <- eBayes(lmFit(m_mat_complete, design_prog))
fit_combat_prog <- eBayes(lmFit(corrected_combat, design_prog))
fit_icombat_v1_prog <- eBayes(lmFit(m_mat_icombat_v1, design_prog))
fit_icombat_v2_prog <- eBayes(lmFit(m_mat_icombat_v2, design_prog))

results_raw_prog        <- topTable(fit_raw_prog, coef = "progression_0_no_1_yes", number = Inf, adjust.method = "BH")
results_combat_prog     <- topTable(fit_combat_prog, coef = "progression_0_no_1_yes", number = Inf, adjust.method = "BH")
results_icombat_v1_prog <- topTable(fit_icombat_v1_prog, coef = "progression_0_no_1_yes", number = Inf, adjust.method = "BH")
results_icombat_v2_prog <- topTable(fit_icombat_v2_prog, coef = "progression_0_no_1_yes", number = Inf, adjust.method = "BH")

gc_raw_prog        <- gc_lambda(results_raw_prog$P.Value)
gc_combat_prog     <- gc_lambda(results_combat_prog$P.Value)
gc_icombat_v1_prog <- gc_lambda(results_icombat_v1_prog$P.Value)
gc_icombat_v2_prog <- gc_lambda(results_icombat_v2_prog$P.Value)

## Death
design_death <- model.matrix(~ death_0_no_1_yes + age_numeric + sex_numeric, data = meta_complete)
fit_raw_death <- eBayes(lmFit(m_mat_complete, design_death))
fit_combat_death <- eBayes(lmFit(corrected_combat, design_death))
fit_icombat_v1_death <- eBayes(lmFit(m_mat_icombat_v1, design_death))
fit_icombat_v2_death <- eBayes(lmFit(m_mat_icombat_v2, design_death))

results_raw_death        <- topTable(fit_raw_death, coef = "death_0_no_1_yes", number = Inf, adjust.method = "BH")
results_combat_death     <- topTable(fit_combat_death, coef = "death_0_no_1_yes", number = Inf, adjust.method = "BH")
results_icombat_v1_death <- topTable(fit_icombat_v1_death, coef = "death_0_no_1_yes", number = Inf, adjust.method = "BH")
results_icombat_v2_death <- topTable(fit_icombat_v2_death, coef = "death_0_no_1_yes", number = Inf, adjust.method = "BH")

gc_raw_death        <- gc_lambda(results_raw_death$P.Value)
gc_combat_death     <- gc_lambda(results_combat_death$P.Value)
gc_icombat_v1_death <- gc_lambda(results_icombat_v1_death$P.Value)
gc_icombat_v2_death <- gc_lambda(results_icombat_v2_death$P.Value)

## Tumor Location
design_tumor <- model.matrix(~ tumor_location_numeric + age_numeric + sex_numeric, data = meta_complete)
fit_raw_tumor <- eBayes(lmFit(m_mat_complete, design_tumor))
fit_combat_tumor <- eBayes(lmFit(corrected_combat, design_tumor))
fit_icombat_v1_tumor <- eBayes(lmFit(m_mat_icombat_v1, design_tumor))
fit_icombat_v2_tumor <- eBayes(lmFit(m_mat_icombat_v2, design_tumor))

results_raw_tumor        <- topTable(fit_raw_tumor, coef = "tumor_location_numeric", number = Inf, adjust.method = "BH")
results_combat_tumor     <- topTable(fit_combat_tumor, coef = "tumor_location_numeric", number = Inf, adjust.method = "BH")
results_icombat_v1_tumor <- topTable(fit_icombat_v1_tumor, coef = "tumor_location_numeric", number = Inf, adjust.method = "BH")
results_icombat_v2_tumor <- topTable(fit_icombat_v2_tumor, coef = "tumor_location_numeric", number = Inf, adjust.method = "BH")

gc_raw_tumor        <- gc_lambda(results_raw_tumor$P.Value)
gc_combat_tumor     <- gc_lambda(results_combat_tumor$P.Value)
gc_icombat_v1_tumor <- gc_lambda(results_icombat_v1_tumor$P.Value)
gc_icombat_v2_tumor <- gc_lambda(results_icombat_v2_tumor$P.Value)

# nSV
impute_with_median <- function(mat, name) {
  n_na_before <- sum(is.na(mat))
  if (n_na_before == 0) { cat("No missing values\n"); return(mat) }
  for (i in 1:nrow(mat)) {
    na_idx <- is.na(mat[i, ])
    if (any(na_idx)) mat[i, na_idx] <- median(mat[i, !na_idx], na.rm = TRUE)
  }
  n_na_after <- sum(is.na(mat))
  cat(sprintf("Imputed missing values: %d → %d NA values\n", n_na_before, n_na_after))
  mat
}
m_mat_complete_imputed   <- impute_with_median(m_mat_complete, "Raw data")
corrected_combat_imputed <- impute_with_median(corrected_combat, "ComBat")
m_mat_icombat_v1_imputed <- impute_with_median(m_mat_icombat_v1, "iComBat v1")
m_mat_icombat_v2_imputed <- impute_with_median(m_mat_icombat_v2, "iComBat v2")
mod_sva <- model.matrix(~ progression_0_no_1_yes + death_0_no_1_yes +
                          tumor_location_numeric + age_numeric + sex_numeric,
                        data = meta_complete)
n_sva_raw <- get_n_sva(m_mat_complete_imputed, mod_sva)
n_sva_combat <- get_n_sva(corrected_combat_imputed, mod_sva)
n_sva_icombat_v1 <- get_n_sva(m_mat_icombat_v1_imputed, mod_sva)
n_sva_icombat_v2 <- get_n_sva(m_mat_icombat_v2_imputed, mod_sva)

# sample to sample density
combat_vals    <- as.vector(corrected_combat)
icombat_v1_vals<- as.vector(m_mat_icombat_v1)
icombat_v2_vals<- as.vector(m_mat_icombat_v2)

cor_v1   <- cor(combat_vals,    icombat_v1_vals, use = "complete.obs")
cor_v2   <- cor(combat_vals,    icombat_v2_vals, use = "complete.obs")
cor_v1v2 <- cor(icombat_v1_vals, icombat_v2_vals, use = "complete.obs")

n_points <- length(combat_vals)
if (n_points > 100000) {
  set.seed(123); sample_idx <- sample.int(n_points, 100000)
  plot_combat <- combat_vals[sample_idx]
  plot_icombat_v1 <- icombat_v1_vals[sample_idx]
  plot_icombat_v2 <- icombat_v2_vals[sample_idx]
} else {
  plot_combat <- combat_vals
  plot_icombat_v1 <- icombat_v1_vals
  plot_icombat_v2 <- icombat_v2_vals
}

p_scatter1 <- ggplot(data.frame(ComBat = plot_combat, iComBat_v1 = plot_icombat_v1),
                     aes(ComBat, iComBat_v1)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
  scale_fill_viridis(name = "Density", option = "plasma") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "ComBat corrected values", y = "iComBat corrected values (EPIC→450K)", title = ifelse(cor_v1 > 0.999, "r > 0.999", sprintf("r = %.3f", cor_v1))) +
  theme_light_plot() + theme(aspect.ratio = 1) + coord_fixed()

p_scatter2 <- ggplot(data.frame(ComBat = plot_combat, iComBat_v2 = plot_icombat_v2),
                     aes(ComBat, iComBat_v2)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
  scale_fill_viridis(name = "Density", option = "plasma") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "ComBat corrected values", y = "iComBat corrected values (450K→EPIC)", title = ifelse(cor_v2 > 0.999, "r > 0.999", sprintf("r = %.3f", cor_v2))) +
  theme_light_plot() + theme(aspect.ratio = 1) + coord_fixed()

print(p_scatter1)
print(p_scatter2)
