# This file is a modified version of the original ComBat function from the sva package.
# The original code is licensed under the Artistic License 2.0.
# This modified version is also distributed under the same license.
# See: https://opensource.org/licenses/Artistic-2.0
#
# Original citation:
# Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2025). sva: Surrogate Variable Analysis. doi:10.18129/B9.bioc.sva, R package version 3.56.0, https://bioconductor.org/packages/sva.
# Citation of this modification:
# Tomo Y. & Nakaki R. (2025). iComBat: An Incremental Framework for Batch Effect Correction in DNA Methylation Array Data. doi:10.1101/2025.05.06.652337, bioRxiv. 


library(matrixStats)
library(BiocParallel)
library(sva)

# Utility functions

aprior <- function(gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2 * s2 + m^2) / s2
}

bprior <- function(gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m * s2 + m^3) / s2
}

postmean <- function(g.hat, g.bar, n, d.star, t2) {
  (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
}

postvar <- function(sum2, n, a, b) {
  (0.5 * sum2 + b) / (n / 2 + a - 1)
}

it.sol <- function(sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) {
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  while (change > conv) {
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, na.rm = TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old) / pmax(abs(g.old), 1e-8),
                  abs(d.new - d.old) / pmax(abs(d.old), 1e-8))
    g.old <- g.new
    d.old <- d.new
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}

int.eprior <- function(sdat, g.hat, d.hat) {
  g.star <- d.star <- numeric(nrow(sdat))
  for (i in seq_len(nrow(sdat))) {
    g <- g.hat[-i]
    d <- d.hat[-i]
    x <- sdat[i, !is.na(sdat[i, ])]
    if (length(x) == 0) next
    dat <- matrix(as.numeric(x), length(g), length(x), byrow = TRUE)
    resid2 <- (dat - g)^2
    LH <- 1 / (2 * pi * d)^(length(x) / 2) * exp(-rowSums(resid2) / (2 * d))
    g.star[i] <- sum(g * LH) / sum(LH)
    d.star[i] <- sum(d * LH) / sum(LH)
  }
  adjust <- rbind(g.star, d.star)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}

# iComBat: Initialize and Incremental function

ComBatInitialize <- function(dat, 
                             batch, 
                             mod = NULL, 
                             par.prior = TRUE, 
                             mean.only = FALSE,
                             ref.batch = NULL,
                             BPPARAM = bpparam("SerialParam")) {
  
  if (length(dim(batch)) > 1) {
    stop("This incremental version of ComBat only allows one batch variable")
  }
  
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
  
  zero.rows.lst <- lapply(levels(batch), function(b_level){
    if(sum(batch == b_level) > 1){
      return(which(apply(dat[, batch == b_level], 1, var) == 0))
    } else{
      return(integer(0))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(seq_len(nrow(dat)), zero.rows)
  if (length(zero.rows) > 0) {
    message(sprintf("Found %d genes with zero variance in at least one batch; ignoring those for fit.",
                    length(zero.rows)))
  }
  dat <- dat[keep.rows, ]
  
  if(any(table(batch) == 1)) {
    mean.only <- TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  
  if (mean.only) {
    message("Using the 'mean only' version of ComBat")
  }
  
  batchmod <- model.matrix(~ -1 + batch)
  n.batch <- nlevels(batch)
  batches <- lapply(seq_len(n.batch), function(i) which(batch == levels(batch)[i]))
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  ref <- NULL
  if(!is.null(ref.batch)) {
    if(!(ref.batch %in% levels(batch))) {
      stop("reference batch is not one of the levels of the batch variable")
    }
    ref <- which(levels(batch) == ref.batch)
    message("Using batch =", ref.batch, "as reference.")
    batchmod[, ref] <- 1
  }
  
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)) {
    check[ref] <- FALSE
  }
  design <- as.matrix(design[, !check, drop = FALSE])
  
  if (qr(design)$rank < ncol(design)) {
    stop("Design matrix is confounded! Remove or adjust some covariates/batches.")
  }
  
  B.hat <- solve(crossprod(design), t(design) %*% t(dat))
  
  if(!is.null(ref)) {
    grand.mean <- B.hat[ref, ]
  } else {
    grand.mean <- as.numeric((n.batches / n.array) %*% B.hat[1:n.batch, ])
  }
  
  fitted.values <- t(design %*% B.hat)
  resid <- dat - fitted.values
  var.pooled <- rowMeans(resid^2, na.rm = TRUE)
  
  fit <- list(
    mean.only = mean.only,
    par.prior = par.prior,
    ref.batch = ref.batch,
    B.hat = B.hat,
    grand.mean = grand.mean,
    var.pooled = var.pooled,
    keep.rows <- keep.rows
  )
  
  return(fit)
}

ComBatIncremental <- function(dat.new, new.batch, fit, mod.new = NULL) {
  
  dat.new <- as.matrix(dat.new)
  if (nrow(dat.new) == 0) {
    stop("No features in dat.new.")
  }
  mean.only  <- fit$mean.only
  par.prior  <- fit$par.prior
  grand.mean <- fit$grand.mean
  var.pooled <- fit$var.pooled
  
  stand.mean <- matrix(grand.mean, nrow = length(grand.mean), ncol = ncol(dat.new))
  s.data <- (dat.new - stand.mean) / sqrt(var.pooled) 
  
  gamma.hat <- t(rowMeans(s.data))
  delta.hat <- t(rowVars(s.data))
  
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  
  gamma.hat.new <- rowMeans(s.data, na.rm = TRUE)
  delta.hat.new <- apply(s.data, 1, var, na.rm = TRUE)
  
  if (par.prior) {
    gamma.star.new <- numeric(length(gamma.hat.new))
    delta.star.new <- numeric(length(delta.hat.new))
    if (mean.only) {
      gamma.star.new <- postmean(gamma.hat.new, gamma.bar, 1, 1, t2)
      delta.star.new <- 1
    } else {
      temp <- it.sol(
        sdat = s.data,
        g.hat = gamma.hat.new,
        d.hat = delta.hat.new,
        g.bar = gamma.bar,
        t2 = t2,
        a = a.prior,
        b = b.prior
      )
      gamma.star.new <- temp[1, ]
      delta.star.new <- temp[2, ]
    }
  } else {
    gamma.star.new <- numeric(length(gamma.hat.new))
    delta.star.new <- numeric(length(delta.hat.new))
    if (mean.only) {
      gamma.star.new <- gamma.hat.new
      delta.star.new <- 1
    } else {
      temp <- int.eprior(as.matrix(s.data),
                         gamma.hat.new,
                         delta.hat.new)
      gamma.star.new <- temp[1, ]
      delta.star.new <- temp[2, ]
    }
  }
  
  corrected.new <- s.data
  for (g in seq_len(nrow(s.data))) {
    corrected.new[g, ] <- (s.data[g, ] - gamma.star.new[g]) / sqrt(delta.star.new[g])
    corrected.new[g, ] <- corrected.new[g, ] * sqrt(var.pooled[g]) + grand.mean[g]
  }
  return(corrected.new)
}

