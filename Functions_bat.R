# --- Required packages ---
library(extraDistr)
library(SCI)
library(mev)
library(dplyr)
library(ggplot2)
library(zoo)
library(readr)
library(stats)
library(tidyverse)
library(nortest)   # for AD test
library(furrr)     # for parallel processing

# =============================================================
#   Helper Functions
# =============================================================

.mon.vect <- function(first.mon, n) {
  if (!(first.mon %in% (1:12))) stop("'first.mon' should be in 1:12")
  mm <- (1:12) + (first.mon - 1)
  mm[mm > 12] <- mm[mm > 12] - 12
  rep(mm, length.out = n)
}

# ------------------------------
# BAT Distribution Functions
# ------------------------------
Upsilon <- function(y) log1p(exp(y))
Upsilon_prime <- function(y) 1 / (1 + exp(-y))
Upsilon_inv <- function(x) { if (any(x <= 0)) stop("Upsilon_inv: x must be > 0"); log(expm1(x)) }

H_theta <- function(y, theta) {
  alpha1 <- theta[1]; beta1 <- theta[2]; gamma1 <- theta[3]
  alpha2 <- theta[4]; beta2 <- theta[5]; gamma2 <- theta[6]
  (1 + gamma2 * Upsilon((y - alpha2) / beta2))^(1 / gamma2) -
    (1 + gamma1 * Upsilon((alpha1 - y) / beta1))^(1 / gamma1)
}

H_theta_prime <- function(y, theta) {
  alpha1 <- theta[1]; beta1 <- theta[2]; gamma1 <- theta[3]
  alpha2 <- theta[4]; beta2 <- theta[5]; gamma2 <- theta[6]
  u1 <- (y - alpha2) / beta2
  u2 <- (alpha1 - y) / beta1
  (1 / beta2) * (1 + gamma2 * Upsilon(u1))^(1 / gamma2 - 1) * Upsilon_prime(u1) +
    (1 / beta1) * (1 + gamma1 * Upsilon(u2))^(1 / gamma1 - 1) * Upsilon_prime(u2)
}

dbat <- function(y, theta, nu) {
  Hval <- H_theta(y, theta)
  Hprime <- H_theta_prime(y, theta)
  dens <- dt(Hval, df = nu) * Hprime
  dens[!is.finite(dens) | dens <= 0] <- NA
  return(dens)
}

pbat <- function(y, theta, nu) pt(H_theta(y, theta), df = nu)

loglik_bat <- function(theta, y, nu) {
  dens <- dbat(y, theta, nu)
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)
  -sum(log(dens))
}

fit_bats <- function(y, init_theta = NULL, init_nu = 3) {
  if (is.null(init_theta))
    init_theta <- c(median(y) - 1, 1, 0.1, median(y) + 1, 1, 0.1)
  init <- c(init_theta, init_nu)
  fit <- optim(
    par = init,
    fn = function(par, y) {
      theta <- par[1:6]; nu <- par[7]
      if (nu <= 1 || any(par[c(2,5)] <= 0)) return(1e10)
      loglik_bat(theta, y, nu)
    },
    y = y,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6, -5, -Inf, 1e-6, -5, 1.01),
    upper = c(Inf, Inf, 5, Inf, Inf, 5, 1000),
    control = list(maxit = 1000)
  )
  list(theta = fit$par[1:6], nu = fit$par[7])
}

fitSCI_bat <- function(x, first.mon, time.scale = 1, scaling = c("no","max","sd")) {
  x <- as.numeric(x)
  scaling <- match.arg(scaling)
  scale.val <- switch(scaling,
                      no = 1, max = max(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
  x <- x / scale.val
  if (time.scale > 1) x <- as.numeric(stats::filter(x, rep(1, time.scale)/time.scale, sides=1))
  mmon <- .mon.vect(first.mon, length(x))
  param_list <- vector("list", 12)
  for (mm in 1:12) {
    xm <- x[mmon == mm]; xm <- xm[is.finite(xm)]
    if (length(xm) < 3 || all(xm == xm[1])) { param_list[[mm]] <- rep(NA,7); next }
    fit <- try(fit_bats(xm), silent=TRUE)
    if (inherits(fit, "try-error")) param_list[[mm]] <- rep(NA,7)
    else param_list[[mm]] <- c(fit$theta, fit$nu)
  }
  list(dist.para = do.call(cbind, param_list), time.scale = time.scale, scaling = scaling)
}

transformSCI_bat <- function(x, first.mon, obj, scaling=c("no","max","sd")) {
  x <- as.numeric(x)
  scaling <- match.arg(scaling)
  scale.val <- switch(scaling, no=1, max=max(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE))
  x <- x/scale.val
  if (obj$time.scale > 1) x <- as.numeric(stats::filter(x, rep(1, obj$time.scale)/obj$time.scale, sides=1))
  mmon <- .mon.vect(first.mon, length(x))
  pars <- obj$dist.para
  for (mm in 1:12) {
    idx <- which(mmon==mm)
    if (anyNA(pars[,mm])) {x[idx] <- NA; next}
    theta <- pars[1:6,mm]; nu <- pars[7,mm]
    x[idx] <- tryCatch(pbat(x[idx], theta, nu), error=function(e) rep(NA,length(idx)))
  }
  qnorm(x)
}
