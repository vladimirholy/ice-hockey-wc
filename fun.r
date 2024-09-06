
model_rank <- function(data, time_col = "time", team_col = "team", rank_col = "rank", predictor_cols = character(), dynamics = "mean_reverting", lambda = 0, bounded = TRUE, coef_start = NULL) {
  
  time_vec <- sort(unique(data[[time_col]]))
  team_vec <- sort(unique(data[[team_col]]))
  
  t <- length(time_vec)
  n <- length(team_vec)
  m <- length(predictor_cols)
  
  data_idx <- expand.grid(time_vec, team_vec)
  colnames(data_idx) <- c(time_col , team_col)
  
  data_full <- merge(x = data_idx, y = data, by = c(time_col , team_col), all.x = TRUE)
  
  y <- matrix(data_full[[rank_col]], nrow = t, ncol = n, byrow = TRUE)
  rownames(y) <- time_vec
  colnames(y) <- team_vec
  
  x <- list()
  if (m > 0) {
    for (i in 1:m) {
      x[[i]] <- matrix(data_full[[predictor_cols[i]]], nrow = t, ncol = n, byrow = TRUE)
      rownames(x[[i]]) <- time_vec
      colnames(x[[i]]) <- team_vec
    }
    names(x) <- predictor_cols
  }
  
  data <- list(y = y, x = x)
  
  if (dynamics == "static") {
    theta_start <- rep(0, n + m - 1)
    theta_lb <- rep(-Inf, n + m - 1)
    theta_ub <- rep(Inf, n + m - 1)
    coef_name <- c(team_vec, predictor_cols)
  } else if (dynamics == "mean_reverting") {
    theta_start <- c(rep(0, n + m - 1), 0.50, 0.50)
    theta_lb <- c(rep(-Inf, n + m - 1), 1e-12, 1e-12)
    theta_ub <- c(rep(Inf, n + m - 1), 1 - 1e-12, Inf)
    coef_name <- c(team_vec, predictor_cols, "phi", "alpha")
  } else if (dynamics == "persistent") {
    theta_start <- c(rep(0, n + m - 1), 0.50)  
    theta_lb <- c(rep(-Inf, n + m - 1), 1e-12)
    theta_ub <- c(rep(Inf, n + m - 1), Inf)
    coef_name <- c(team_vec, predictor_cols, "kappa")
  }
  
  if (!is.null(coef_start)) {
    theta_start <- coef_start[-n]
  }
  
  if (bounded) {
    optim_res <- nloptr::nloptr(x0 = theta_start, lb = theta_lb, ub = theta_ub, eval_f = likelihood_objective, y = y, x = x, dynamics = dynamics, time_idx = time_vec, lambda = lambda, opts = list(algorithm = 'NLOPT_LN_SBPLX', ftol_rel = 0, ftol_abs = 0, xtol_rel = 0, xtol_rel = 0, maxeval = 1e6))
  } else {
    optim_res <- nloptr::nloptr(x0 = theta_start, eval_f = likelihood_objective, y = y, x = x, dynamics = dynamics, time_idx = time_vec, lambda = lambda, opts = list(algorithm = 'NLOPT_LN_SBPLX', ftol_rel = 0, ftol_abs = 0, xtol_rel = 0, xtol_rel = 0, maxeval = 1e6))
  }
  
  theta_est <- optim_res$solution
  
  if (length(theta_est) > n - 1) {
    coef_est <- c(theta_est[1:(n - 1)], -sum(theta_est[1:(n - 1)]), theta_est[n:length(theta_est)])
  } else {
    coef_est <- c(theta_est[1:(n - 1)], -sum(theta_est[1:(n - 1)]))
  }
  names(coef_est) <- coef_name
  
  hess_res <- try(stats::optimHess(par = theta_est, fn = likelihood_objective, y = y, x = x, dynamics = dynamics, time_idx = time_vec, lambda = lambda))
  
  if (inherits(hess_res, 'try-error')) {
    theta_vcov <- matrix(NA, nrow = length(theta_est), ncol = length(theta_est))
  } else {
    theta_vcov <- pracma::pinv(hess_res)
  }
  
  amat <- rbind(cbind(diag(n - 1), -1, matrix(0, nrow = n - 1, ncol = length(coef_est) - n)), cbind(matrix(0, nrow = length(coef_est) - n, ncol = n), diag(length(coef_est) - n)))
  coef_vcov <- t(amat) %*% theta_vcov %*% amat
  rownames(coef_vcov) <- coef_name
  colnames(coef_vcov) <- coef_name
  
  coef_sd <- sqrt(diag(coef_vcov))
  coef_zstat <- coef_est / coef_sd
  coef_pval <- 2 * stats::pnorm(-abs(coef_zstat))
  
  coef <- list(est = coef_est, sd = coef_sd, zstat = coef_zstat, pval = coef_pval, vcov = coef_vcov)
  
  tv <- likelihood_evaluate(theta = theta_est, y = y, x = x, dynamics = dynamics, time_idx = time_vec)
  
  loglik_sum <- sum(tv$l, na.rm = TRUE)
  aic <- 2 * length(theta_est) - 2 * loglik_sum
  bic <- log(t) * length(theta_est) - 2 * loglik_sum
  
  fit <- list(loglik_sum = loglik_sum, aic = aic, bic = bic)
  
  return(list(data = data, coef = coef, tv = tv, fit = fit))
  
}

likelihood_evaluate <- function(theta, y, x, dynamics, time_idx) {
  
  t <- nrow(y)
  n <- ncol(y)
  m <- length(x)
  
  play <- !is.na(y)
  
  f <- matrix(c(theta[1:(n - 1)], -sum(theta[1:(n - 1)])), nrow = t, ncol = n, byrow = TRUE)
  
  if (m > 0) {
    f[play] <- f[play] + Reduce("+", lapply(1:m, function(i) { theta[n + i - 1] * x[[i]] }))[play]
  }
  
  u <- matrix(0, nrow = t, ncol = n)
  s <- matrix(0, nrow = t, ncol = n)
  l <- rep(0, times = t)
  
  s[1, play[1, ]] <- 1 - exp(f[1, play[1, ]]) * cumsum(1 / rev(cumsum(exp(rev(f[1, play[1, ]][Matrix::invPerm(y[1, play[1, ]])])))))[y[1, play[1, ]]]
  l[1] <- sum(f[1, play[1, ]]) - sum(log(cumsum(exp(rev(f[1, play[1, ]][Matrix::invPerm(y[1, play[1, ]])])))))
  
  if (dynamics == "static") {
    
    for (i in 2:t) {
      
      s[i, play[i, ]] <- 1 - exp(f[i, play[i, ]]) * cumsum(1 / rev(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))[y[i, play[i, ]]]
      l[i] <- sum(f[i, play[i, ]]) - sum(log(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))
      
    }
    
  } else if (dynamics == "mean_reverting") {
    
    for (i in 2:t) {
      
      u[i, ] <- theta[n + m]^(time_idx[i] - time_idx[i - 1] - 1L) * (theta[n + m] * u[i - 1, ] + theta[n + m + 1] * s[i - 1, ])
      f[i, ] <- f[i, ] + u[i, ]
      s[i, play[i, ]] <- 1 - exp(f[i, play[i, ]]) * cumsum(1 / rev(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))[y[i, play[i, ]]]
      l[i] <- sum(f[i, play[i, ]]) - sum(log(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))
      
    }
    
  } else if (dynamics == "persistent") {
    
    for (i in 2:t) {
      
      u[i, ] <- u[i - 1, ] + theta[n + m] * s[i - 1, ]
      f[i, ] <- f[i, ] + u[i, ]
      s[i, play[i, ]] <- 1 - exp(f[i, play[i, ]]) * cumsum(1 / rev(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))[y[i, play[i, ]]]
      l[i] <- sum(f[i, play[i, ]]) - sum(log(cumsum(exp(rev(f[i, play[i, ]][Matrix::invPerm(y[i, play[i, ]])])))))
      
    }
    
  } 
  
  return(list(loglik = l, strength = f, dyn = u, score = s))
  
}

likelihood_objective <- function(theta, y, x, dynamics, time_idx, lambda) {
  
  tv <- likelihood_evaluate(theta = theta, y = y, x = x, dynamics = dynamics, time_idx = time_idx)
  
  obj <- min(-sum(tv$loglik) + lambda * sum(tv$strength[!is.na(y)]^2), 1e100, na.rm = TRUE)
  
  return(obj)
  
}

prob_first_k <- function(y, f, k) {
  
  ok <- !is.na(y) & is.finite(f)
  
  y <- y[ok]
  f <- f[ok]
  
  if(k == 1) {
    
    prob <- exp(f[y == 1]) / sum(exp(f))
    
  } else {
    
    f_exp <- matrix(exp(f), nrow = 1)
    
    y[y > k] <- Inf
    
    y_all <- matrix(y, nrow = factorial(k), ncol = length(f), byrow = TRUE)
    
    y_all[is.finite(y_all)] <- arrangements::permutations(k)
    
    prob <- 0
    
    for (i in 1:nrow(y_all)) {
      
      prob <- prob + exp(distr_pluce_worth_loglik(y = y_all[i, , drop = FALSE], f = f_exp))
      
    }
    
  }
  
  as.vector(prob)
  
}

distr_pluce_worth_loglik <- function(y, f) {
  
  t <- nrow(f)
  n <- ncol(f)
  
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  y_inv <- matrix(NA_real_, nrow = t, ncol = n)
  
  for (i in 1:t) {
    
    k <- sum(is.finite(y[i, ]))
    
    if (k < n) {
      y[i, !is.finite(y[i, ])] <- (k + 1):n
    }
    
    y_inv[i, ] <- Matrix::invPerm(y[i, ])
    
    o1 <- t(matrix(f[i, y_inv[i, ]], nrow = n, ncol = n))
    o1[lower.tri(o1)] <- 0
    o1 <- o1[1:k, 1:n]
    
    o2 <- log(rowSums(o1))
    
    res_loglik[i, ] <- sum(log(f[i, y_inv[i, 1:k]])) - sum(o2)
    
  }
  
  return(res_loglik)
  
}
