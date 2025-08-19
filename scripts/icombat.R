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

# Utility functions

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2*s2 + m^2) / s2
}

bprior <- function(gamma.hat){
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m*s2 + m^3) / s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2 + b) / (n/2 + a - 1)
}

# Inverse gamma distribution density function. (Note: does not do any bounds checking on arguments)
dinvgamma <- function (x, shape, rate = 1/scale, scale = 1) {
  # PDF taken from https://en.wikipedia.org/wiki/Inverse-gamma_distribution
  # Note: alpha = shape, beta = rate
  stopifnot(shape > 0)
  stopifnot(rate > 0)
  ifelse(x <= 0, 0, ((rate ^ shape) / gamma(shape)) * x ^ (-shape - 1) * exp(-rate/x))
}

# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  ## cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

## likelihood function used below
L <- function(x,g.hat,d.hat){
  prod(dnorm(x, g.hat, sqrt(d.hat)))
}

## Monte Carlo integration functions
int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
    ## if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust	
} 

## fits the L/S model in the presence of missing data values

Beta.NA <- function(y,X){
  des <- X[!is.na(y),]
  y1 <- y[!is.na(y)]
  B <- solve(crossprod(des), crossprod(des, y1))
  B
}

# iComBat: Initialize and Incremental function

ComBatInitialize <- function(dat, batch, mod = NULL, mean.only = FALSE) {
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  } 
  
  ## coerce dat into a matrix
  dat <- as.matrix(dat)
  
  ## find genes with zero variance in any of the batches
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level){
    if(sum(batch==batch_level)>1){
      return(which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0})))
    }else{
      return(which(rep(1,3)==2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", length(zero.rows)))
    # keep a copy of the original data matrix and remove zero var rows
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
  
  ## make batch a factor and make a set of indicators for batch
  if(any(table(batch)==1)){mean.only=TRUE}
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  
  n.batch <- nlevels(batch)
  if (n.batch == 1) {
    batchmod <- matrix(1, nrow = length(batch), ncol = 1)
    colnames(batchmod) <- paste0("batch", levels(batch))
  } else {
    batchmod <- model.matrix(~ -1 + batch)
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if (ncol(design) >= n.batch) check[1:n.batch] <- FALSE
  design <- as.matrix(design[, !check, drop = FALSE])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  if (!is.null(mod)) {
    n.covariates.original <- ncol(mod)
    n.covariates <- ncol(design) - ncol(batchmod)
  } else {
    design <- batchmod
    n.covariates.original <- 0
    n.covariates <- 0
  }
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  message('Standardizing Data across genes')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  if (!NAs){
    var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
  } else {
    var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  # Regression coefficients for save
  cov_idx_start <- ncol(design) - n.covariates + 1
  cov_idx_end   <- ncol(design)
  mod.coef <- if (n.covariates > 0) {
    B.hat[cov_idx_start:cov_idx_end, , drop = FALSE]
  } else NULL
  
  list(
    n.batch = n.batch,
    n.genes.full = length(keep.rows) + length(zero.rows),
    keep.rows = keep.rows,
    zero.rows = zero.rows,
    grand.mean = grand.mean,
    var.pooled = var.pooled,
    n.covariates.original = n.covariates.original,
    n.covariates = n.covariates,
    mod.coef = mod.coef,
    mean.only = mean.only
  )
}

ComBatIncremental <- function(dat.new, new.batch, fit, mod.new = NULL, par.prior = TRUE, prior.plots=FALSE, intercept.col = 1) {
  colnames.dat.new <- colnames(dat.new)
  rownames.dat.new <- rownames(dat.new)
  
  # Align genes to the fitted model (respect removed zeros)
  if (nrow(dat.new) == fit$n.genes.full && !is.null(fit$keep.rows)) {
    dat.new.orig <- dat.new
    dat.new <- dat.new[fit$keep.rows, , drop = FALSE]
  } else if (nrow(dat.new) != length(fit$grand.mean)) {
    stop(sprintf("Number of genes (%d) doesn't match the fitted model (%d)",
                 nrow(dat.new), length(fit$grand.mean)))
  }
  
  n.samples.new <- ncol(dat.new)
  
  # Build stand.mean for new samples = grand.mean + covariate part from fitted mod.coef
  if (!is.null(mod.new) && !is.null(fit$mod.coef)) {
    if (nrow(mod.new) != n.samples.new) {
      stop("Number of rows in mod.new must equal number of new samples")
    }
    if (intercept.col >= 1) {
      mod.new.no.intercept <- mod.new[, -intercept.col, drop = FALSE]
    } else {
      mod.new.no.intercept <- mod.new
    }
    if (ncol(mod.new.no.intercept) != fit$n.covariates) {
      stop(sprintf("Number of non-intercept covariates in mod.new (%d) doesn't match fitted model (%d)",
                   ncol(mod.new.no.intercept), fit$n.covariates))
    }
    cov.effects.new <- t(mod.new.no.intercept %*% fit$mod.coef)
  } else {
    cov.effects.new <- matrix(0, nrow = nrow(dat.new), ncol = n.samples.new)
  }
  
  stand.mean <- matrix(fit$grand.mean, nrow = nrow(dat.new), ncol = n.samples.new) + cov.effects.new
  s.data <- as.matrix((dat.new - stand.mean) / sqrt(as.vector(fit$var.pooled)))
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  gamma.hat <- rowMeans(s.data, na.rm = TRUE)
  if (fit$mean.only==TRUE) {
    delta.hat <- rep(1, nrow(s.data))
  } else {
    delta.hat <- rowVars(s.data, na.rm = TRUE)
  }
  
  ##Find Priors
  gamma.bar <- mean(gamma.hat)
  t2 <- var(gamma.hat)
  a.prior <- aprior(delta.hat)
  b.prior <- bprior(delta.hat)
  
  ## Plot empirical and parametric priors
  
  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow=c(2,2))
    
    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    
    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)
    
    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)
    
    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }
  
  ## Find EB batch adjustments
  gamma.star <- delta.star <- matrix(NA, nrow=1, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    if (fit$mean.only) {
      gamma.star <- postmean(gamma.hat, gamma.bar, 1, 1, t2)
      delta.star <- rep(1, nrow(s.data))
    }
    else {
      temp <- it.sol(s.data, gamma.hat,
                     delta.hat, gamma.bar, t2, a.prior,
                     b.prior)
      gamma.star <- temp[1, ]
      delta.star <- temp[2, ]
    }
  } else {
    message("Finding nonparametric adjustments")
    if (fit$mean.only) {
      delta.hat = 1
    }
    temp <- int.eprior(as.matrix(s.data),
                       gamma.hat, delta.hat)
    gamma.star=temp[1,]
    delta.star=temp[2,]
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- as.matrix(s.data)
  bayesdata <- (bayesdata-gamma.star) / sqrt(pmax(delta.star, 1e-8))
  bayesdata <- (bayesdata*(sqrt(as.vector(fit$var.pooled)))) + stand.mean
  
  colnames(bayesdata) <- colnames.dat.new
  rownames(bayesdata) <- rownames.dat.new
  
  ## put genes with 0 variance in any batch back in data
  if (length(fit$zero.rows) > 0) {
    dat.new.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.new.orig
  }
  
  return(bayesdata)
}

