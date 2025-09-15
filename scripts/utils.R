suppressPackageStartupMessages({
  library(ggplot2)
  library(sva)
})

gc_lambda <- function(p) {
  chisq <- qchisq(1 - p, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  return(lambda)
}

get_n_sva <- function(dat, mod) {
  tryCatch(num.sv(dat, mod = mod, method = "leek"), error = function(e) NA)
}

theme_light_plot <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "grey90", color = NA)
    )
}

decompress_idat_gz <- function(base_dir = getwd()) {
  gz_files <- list.files(base_dir, pattern = "(?i)\\.idat\\.gz$", 
                         recursive = TRUE, full.names = TRUE)
  message(sprintf("Found %d gzipped IDAT files", length(gz_files)))
  if (length(gz_files) == 0) {
    message("No gzipped IDAT files found. Skipping decompression.")
    return(invisible(TRUE))
  }
  for (f in gz_files) {
    out <- sub("\\.gz$", "", f)
    if (!file.exists(out)) {
      try(R.utils::gunzip(f, destname = out, overwrite = TRUE, remove = FALSE), silent = TRUE)
    }
  }
  invisible(TRUE)
}

check_PC_covariate_association <- function(pca_data, pcs = c("PC1","PC2","PC3"), covariates = c("Age","Sex","Batch")) {
  results <- list()
  for(pc in pcs) {
    res_pc <- data.frame(Covariate = character(), Type = character(), Statistic = numeric(), p.value = numeric(), stringsAsFactors = FALSE)
    for(cov in covariates) {
      x <- pca_data[[pc]]
      y <- pca_data[[cov]]
      if(is.numeric(y)) {
        test <- cor.test(x, y)
        res_pc <- rbind(res_pc, data.frame(Covariate = cov,
                                           Type = "Correlation",
                                           Statistic = test$estimate,
                                           p.value = test$p.value))
      } 
      else if(length(unique(y)) == 2) {
        group1 <- x[y == unique(y)[1]]
        group2 <- x[y == unique(y)[2]]
        test <- t.test(group1, group2)
        res_pc <- rbind(res_pc, data.frame(Covariate = cov,
                                           Type = "t-test",
                                           Statistic = test$statistic,
                                           p.value = test$p.value))
      } 
      else {
        test <- aov(x ~ y)
        test_sum <- summary(test)[[1]]
        fstat <- test_sum[1, "F value"]
        pval <- test_sum[1, "Pr(>F)"]
        res_pc <- rbind(res_pc, data.frame(Covariate = cov,
                                           Type = "ANOVA",
                                           Statistic = fstat,
                                           p.value = pval))
      }
    }
    results[[pc]] <- res_pc
  }
  return(results)
}