library(MASS)
library(Matrix)
library(matrixcalc)
library(matrixStats)
library(mvtnorm)

blm <- function(formula, data, subset,
                prior.mean=0, prior.precision=0.0001,
                prior.df=0.0001, prior.scale=1,
                cov.structure=1,
                prior.a=prior.df / 2,
                prior.b=(prior.df * prior.scale ^ 2) / 2,
                design=FALSE) {
  
  # add names to data frame if needed
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("X", 1:ncol(data))
  }
  
  # build call to stats::model.frame from passed arguments
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  # extract terms from the model frame
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")   # response vector
  x <- model.matrix(mt, mf, contrasts) # design matrix
  x.names <- colnames(x)
  if (design) {
    return(x)
  }
  
  # make sure arguments are passed from at most one pair
  # (prior.df, prior.scale), (prior.a, prior.b)
  if ((!missing(prior.df) || !missing(prior.scale)) &&
        (!missing(prior.a)  || !missing(prior.b))) {
    stop(paste("Can only set arguments from one pair:",
               "(prior.df, prior.scale) or (prior.a, prior.b)."))
  }
  
  # expand prior.mean, prior.precision, cov.structure
  n <- nrow(x)
  p <- ncol(x)
  
  if (length(prior.mean) == 0) {
    stop("prior.mean cannot have length zero.")
  } else if (length(prior.mean) == 1) {
    prior.mean <- rep(prior.mean, p)
  } else if (length(prior.mean) != p) {
    stop(paste("Length of prior.mean is not 1 and does not",
               "match number of columns in design matrix."))
  }
  names(prior.mean) <- x.names
  
  if (is.null(ncol(prior.precision))) {
    if (length(prior.precision) == 0) {
      stop("prior.precision cannot have length zero.")
    } else if (length(prior.precision) == 1) {
      prior.precision <- rep(prior.precision, p)
    } else if (length(prior.precision) != p) {
      stop(paste("Length of prior.precision is not 1 and does not",
                 "match the number of columns in design matrix."))
    }
    prior.precision <- diag(prior.precision, p)
  } else if (ncol(prior.precision) != p ||
               nrow(prior.precision) != p) {
    stop(paste("Number of rows and columns of prior.precision must",
               "match the number of columns in design matrix."))
  }
  colnames(prior.precision) <- rownames(prior.precision) <- x.names
  
  if (is.null(ncol(cov.structure))) {
    if (length(cov.structure) == 0) {
      stop("cov.structure cannot have length zero.")
    } else if (length(cov.structure) == 1) {
      cov.structure <- rep(cov.structure, n)
    } else if (length(cov.structure) != n) {
      stop(paste("Length of cov.structure is not 1 and does not",
                 "match the number of rows in design matrix."))
    }
    cov.structure <- diag(cov.structure, n)
  } else if (ncol(cov.structure) != n ||
               nrow(cov.structure) != n) {
    stop(paste("Number of rows and columns of cov.structure must",
               "match the number of rows in design matrix."))
  }
  
  # fit model
  fit <- blm.fit(x, y,
                 prior.mean, prior.precision,
                 prior.a, prior.b,
                 cov.structure)
  
  # prepare return value
  prior.params <- list(mean=prior.mean,
                       precision=prior.precision,
                       df=prior.df,
                       scale=prior.scale,
                       a=prior.a,
                       b=prior.b)
  fit$prior <- prior.params
  fit$cov.structure <- cov.structure
  fit$call <- match.call()
  fit$model <- mf
  fit$terms <- mt
  fit$xlevels <- .getXlevels(mt, mf)
  class(fit) <- c("blm")
  fit
}

blm.fit <- function(x, y,
                    prior.mean, prior.precision,
                    prior.a, prior.b,
                    cov.structure) {
  n <- nrow(x)
  p <- ncol(x)
  rank <- as.numeric(rankMatrix(x))
  
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", 1:p)
  }
  x.names <- colnames(x)
  
  # error checking
  if (length(y) != n) {
    stop("Length of y must match number of rows of x.")
  }
  
  #   if (rank < p) {
  #     stop("x is rank deficient.")
  #   }
  
  if (!is.positive.semi.definite(cov.structure)) {
    stop("cov.structure must be positive semidefinite.")
  }
  
  if (!is.positive.definite(prior.precision)) {
    stop("prior.precision must be positive definite.")
  }
  
  # fit model
  Sigma.inv <- svd.inverse(cov.structure)
  
  H.inv <- t(x) %*% Sigma.inv %*% x + prior.precision
  H <- svd.inverse(H.inv)
  h <- t(x) %*% Sigma.inv %*% y + prior.precision %*% prior.mean
  a <- drop(prior.a + n / 2)
  b <- drop(prior.b + (t(y) %*% Sigma.inv %*% y +
                         t(prior.mean) %*% prior.precision %*% prior.mean -
                         t(h) %*% H %*% h
  ) / 2)
  
  # prepare return value
  mean <- as.vector(H %*% h)
  names(mean) <- x.names
  precision <- H.inv
  colnames(precision) <- rownames(precision) <- x.names
  df <- 2 * a
  scale <- sqrt(b / a)
  
  t.mean <- mean
  t.cov  <- scale ^ 2 * H
  t.df   <- 2 * a
  colnames(t.cov) <- rownames(t.cov) <- x.names
  
  coefficients <- mean
  fitted.values <- as.vector(x %*% coefficients)
  residuals <- y - fitted.values
  
  posterior.params <- list(mean=mean,
                           precision=precision,
                           df=df,
                           scale=scale,
                           H=H,
                           h=h,
                           a=a,
                           b=b,
                           t.mean=t.mean,
                           t.cov=t.cov,
                           t.df=t.df)
  
  list(coefficients=coefficients,
       residuals=residuals,
       fitted.values=fitted.values,
       rank=rank,
       df.residual=df,
       posterior=posterior.params)
}

print.blm <- function(blm) {
  cat("\nCall:\n")
  print(blm$call)
  cat("\nCoefficients:\n")
  print(blm$posterior$mean)
}

summary.blm <- function(blm) {
  cat("\nCall:\n")
  print(blm$call)
  
  se <- sqrt(diag(blm$posterior$t.cov) *
               (blm$posterior$a / (blm$posterior$a - 1)))
  coeff <- cbind(blm$posterior$t.mean, se)
  colnames(coeff) <- c("Post. Mean", "Mar. Post. SE")
  cat("\nCoefficients:\n")
  print(coeff)
  
  cat("\nEstimated Variance:\n")
  cat(blm$posterior$b / blm$posterior$a)
  
  cat("\n\nResidual degrees of freedom:\n")
  cat(blm$posterior$t.df - length(blm$posterior$mean))
}

rposterior.blm <- function(blm, n, regression.only=TRUE) {
  # TODO: sample for variance if regression.only == FALSE
  post <- blm$posterior
  sample <- rmvt(n=n, sigma=post$t.cov, df=post$t.df, delta=post$t.mean)
  colnames(sample) <- names(blm$coefficients)
  return(sample)
}

qf.min <- function(z, delta, mean, cov) {
  # Minimizes the quadratic form used in criteria for
  # membership in credible subgroups
  #
  # Args:
  #   z: A predictive covariate vector (include the intercept)
  #   delta: The treatment effect threshold
  #   mean: The mean vector of the posterior predictive regression
  #         parameter distribution
  #   cov: The covariance matrix of the posterior predictive regression
  #        parameter distribution
  #
  # Returns:
  #   The minimum of the quadratic form
  
  nz2 <- sum(z ^ 2)
  dz <- (delta / nz2) * z
  v <- z / sqrt(nz2)
  P <- diag(length(z)) - v %o% v
  cov.inv <- solve(cov)
  beta <- ginv(t(P) %*% cov.inv %*% P) %*% t(P) %*% cov.inv %*% (mean - dz)
  
  t(P %*% beta + dz - mean) %*% cov.inv %*% (P %*% beta + dz - mean)
}

find.credible.subgroups.pb <- function(Z, post.sample, post.mean, post.cov,
                                    delta=0, cred=0.95,
                                    tolerance=0.001, max.iters=20,
                                    verbose=FALSE) {
  # Finds exclusive and inclusive credible subgroups via
  # "pure Bayesian (PB)" method
  #
  # Args:
  #   Z: Matrix of predictive covariate points (include intercept)
  #   post.sample: A sample from the posterior distribution
  #                of predictive regression parameters
  #   post.mean: The posterior mean vector of predictive regression parameters
  #   post.cov: The covariance matrix of predictive regression parameters
  #   delta: Treatment effect threshold
  #   cred: credibility level of the credible subgroups
  #   tolerance: Maximum difference between target and final credibility
  #   max.iters: Maximum number of algorithm iterations
  #   verbose: (logical) if true, print progress
  #
  # Returns: Logical vectors indicating inclusion of points of Z
  #          in the exclusive and inclusive credible subgroups
  M <- nrow(post.sample)
  n <- nrow(Z)
  q <- ncol(Z)
  
  score <- apply(Z, 1, qf.min, delta=delta, mean=post.mean, cov=post.cov)
  est.qual <- Z %*% post.mean
  if (is.null(cred)) {
    print("cred")
  }
  if (is.null(q)) {
    print("q")
    print(Z)
  }
  if (is.null(n)) {
    print("n")
    print(Z)
  }
  m.u <- q * qf(cred, q, n)
  m.l <- 0
  
  iter <- 0

  repeat {
    iter <- iter + 1
    m <- (m.u + m.l) / 2
    
    # Test credible subgroups
    in.w.D <- score > m & est.qual > delta
    in.w.Sc <- score >= m & est.qual <= delta
    working.D <- Z[in.w.D, ]
    working.Sc <- Z[in.w.Sc, ]
    
    csp.check <- function(post.draw, working.D, working.Sc, delta) {
      return(all(working.D %*% post.draw > delta) &&
               all(working.Sc %*% post.draw <= delta))
    }
    
    win <- apply(post.sample, 1, csp.check,
                 working.D=working.D, working.Sc=working.Sc, delta=delta)
    
    p <- mean(win)
    if (verbose) {
      print(p)
    }
    
    # Tighten bounds
    if (p > cred) {
      m.u <- m
    } else if (p < cred) {
      m.l <- m
    }
    
    if (cred <= p && p < cred + tolerance) {
      break;
    }
    
    if (iter == max.iters) {
      warning(paste("Credible subgroup algorithm did not converge within maximum number of iterations.",
                    "Final credible level:", p))
      break;
    }
  }
  
  return(list(D=in.w.D,
              S=!in.w.Sc))
}

find.credible.subgroups.hpd <- function(Z, post.mean, post.cov, df,
                                           delta=0, cred=0.95) {
  # Finds exclusive and inclusive credible subgroups via
  # "highest posterior density (HPD)" method
  #
  # Args:
  #   Z: Matrix of predictive covariate points (include intercept)
  #   post.mean: The posterior mean vector of predictive regression parameters
  #   post.cov: The covariance matrix of predictive regression parameters
  #   df: The posterior degrees of freedom of the
  #       predictive regression parameter distribution
  #   delta: Treatment effect threshold
  #   cred: credibility level of the credible subgroups
  #
  # Returns: Logical vectors indicating inclusion of points of Z
  #          in the exclusive and inclusive credible subgroups
  q <- ncol(Z)
  n <- nrow(Z)
  
  m <- q * qf(cred, q, df)
  
  score <- apply(Z, 1, qf.min, delta=delta, mean=post.mean, cov=post.cov)
  est.qual <- Z %*% post.mean
  
  in.D <- score > m & est.qual > delta
  in.Sc <- score >= m & est.qual <= delta
  
  return(list(D=in.D,
              S=!in.Sc))
}

find.credible.subgroups.rcs <- function(Z, post.sample, post.mean, post.cov,
                                         delta=0, cred=0.95) {
  # Finds exclusive and inclusive credible subgroups via
  # "restricted covariate space (RCS)" method
  #
  # Args:
  #   Z: Matrix of predictive covariate points (include intercept)
  #   post.sample: A sample from the posterior distribution
  #                of predictive regression parameters
  #   post.mean: The posterior mean vector of predictive regression parameters
  #   post.cov: The covariance matrix of predictive regression parameters
  #   delta: Treatment effect threshold
  #   cred: credibility level of the credible subgroups
  #
  # Returns: Logical vectors indicating inclusion of points of Z
  #          in the exclusive and inclusive credible subgroups
  se <- sqrt(rowSums((Z %*% post.cov) * Z))
  dist <- t(abs(Z %*% (t(post.sample) - post.mean)) / se)
  
  sups <- rowMaxs(dist)
  mdq <- quantile(sups, cred)
  
  exclusive <- as.vector(Z %*% post.mean - mdq * se > delta)
  inclusive <- as.vector(Z %*% post.mean + mdq * se >= delta)
  
  return(list(D=exclusive, S=inclusive))
}

find.credible.subgroups.pointwise <- function(Z, post.mean, post.cov, df,
                                              delta=0, cred=0.95) {
  # Finds exclusive and inclusive credible subgroups via
  # a competing "pointwise" method---not accounting for multiplicity
  #
  # Args:
  #   Z: Matrix of predictive covariate points (include intercept)
  #   post.mean: The posterior mean vector of predictive regression parameters
  #   post.cov: The covariance matrix of predictive regression parameters
  #   df: The posterior degrees of freedom of the
  #       predictive regression parameter distribution
  #   delta: Treatment effect threshold
  #   cred: credibility level of the credible subgroups
  #
  # Returns: Logical vectors indicating inclusion of points of Z
  #          in the exclusive and inclusive credible subgroups
  se <- sqrt(rowSums((Z %*% post.cov) * Z))
  mdq <- qt(cred, df)
  
  exclusive <- as.vector(Z %*% post.mean - mdq * se > delta)
  inclusive <- as.vector(Z %*% post.mean + mdq * se >= delta)
  
  return(list(D=exclusive, S=inclusive))
}

find.credible.subgroups <- function(predictors, fit, treatment,
                                    delta=0, cred=0.95,
                                    post.sample.size=10000,
                                    method="rcs", df=NA,
                                    tolerance=0.001, max.iters=20,
                                    verbose=FALSE) {
  # (Wrapper) Finds exclusive and inclusive credible subgroups
  # via specified method
  #
  # Args:
  #   predictors: Data frame containing predictive covariate points
  #   fit: A "blm" object
  #   delta: Treatment effect threshold
  #   cred: credibility level of the credible subgroups
  #   method: String "PB", "RCS", or "HPD" specifying
  #           which method to call. Can also be "pointwise" but
  #           this method does not give true credible subgroups
  #   df: ("HPD" and "pointwise methods only, then required)
  #       The posterior degrees of freedom of the predictive regression
  #       parameter distribution
  #   tolerance: ("PB" method only) Maximum difference between target
  #              and final credibility level
  #   max.iters: ("PB" method only) Maximum number of algorithm iterations
  #   verbose: ("PB" method only) If true, print progress
  #
  # Returns: Logical vectors indicating inclusion of points of Z
  #          in the exclusive and inclusive credible subgroups
  
  if (method == "hpd" && is.na(df)) {
    stop("Method \"hpd\" requires \"df\" argument.")
  } else if (method == "pointwise" && is.na(df)) {
    stop("Method \"pointwise\" requires \"df\" argument.")
  }
  
  # Extract list of predictive covariates
  factors <- attr(delete.response(terms(fit)), "factors")
  pred.term <- factors[treatment, ]
  pred.term <- as.logical(pred.term)
  pred.vars <- names(which(rowSums(factors[, pred.term]) > 0))
  pred.vars <- pred.vars[which(pred.vars != treatment)]
  
  # Remove columns that are not predictive covariates
  predictors <- predictors[, pred.vars]
  
  # Extract list of purely prognostic covariates
  prog.vars <- rownames(factors)[which(!(rownames(factors) %in% pred.vars))]
  prog.vars <- prog.vars[which(prog.vars != treatment)]
  
  # Make sure all predictive covariates are included
  incl.vars <- colnames(predictors)
  if (!all(pred.vars %in% incl.vars)) {
    stop("All predictive covariates must be included in new data frame.")
  }
  
  # Append column of treatments
  if (treatment %in% names(fit$xlevels)) {
    predictors <- cbind(fit$xlevels[[treatment]][2], predictors)
  } else {
    predictors <- cbind(1, predictors)
  }
  
  colnames(predictors)[1] <- treatment
  
  # Fill out prognostic covariates with 0's
  for (prog in rev(prog.vars)) {
    if (prog %in% names(fit$xlevels)) {
      predictors <- cbind(fit$xlevels[[prog]][1], predictors)
    } else {
      predictors <- cbind(0, predictors)
    }
    colnames(predictors)[1] <- prog
  }
  
  # Create treatment and control model matrices
  trt.mm <- model.matrix(delete.response(terms(fit)), predictors, xlev=fit$xlevels)
  if (treatment %in% names(fit$xlevels)) {
    predictors[, treatment] <- fit$xlevels[[treatment]][1]
  } else {
    predictors[, treatment] <- 0
  }
  ctl.mm <- model.matrix(delete.response(terms(fit)), predictors, xlev=fit$xlevels)
  
  # Create predictor matrix
  predictors <- trt.mm - ctl.mm
  
  # Get posterior parameters
  post.mean <- fit$posterior$t.mean
  post.cov  <- fit$posterior$t.cov
  post.df   <- fit$posterior$t.df
  
  # Trim columns of zero predictors
  rel.cols <- which(colSums(abs(predictors)) != 0)
  post.mean <- post.mean[rel.cols]
  post.cov  <- post.cov[rel.cols, rel.cols]
  predictors <- predictors[, rel.cols]
  
  # Generate posterior sample
  if (method %in% c("rcs", "pb")) {
    post.sample <- rposterior.blm(fit, post.sample.size)
    post.sample <- post.sample[, rel.cols]
  }
  
  # Dispatch to specified method
  switch(method,
         pb = find.credible.subgroups.pb(predictors, post.sample,
                                               post.mean, post.cov,
                                               delta, cred,
                                               tolerance, max.iters, verbose),
         hpd = find.credible.subgroups.hpd(predictors,
                                                 post.mean, post.cov, df,
                                                 delta, cred),
         rcs = find.credible.subgroups.rcs(predictors, post.sample,
                                             post.mean, post.cov,
                                             delta, cred),
         pointwise = find.credible.subgroups.pointwise(predictors,
                                                       post.mean, post.cov,
                                                       df, delta, cred),
         stop("Invalid method name.")
  )
}

credible.subgroup.pair.size <- function(D, S) {
  # Estimates the size (measure) of the credible subgroup pair.
  #
  # Args:
  #   D: Logical vector indicating membership of sample
  #      from covariate space in exclusive credible subgroup
  #   S: Logical vector indicating membership of sample
  #      from covariate space in inclusive credible subgroup
  #
  # Returns: Size of credible subgroup pair
  return(mean(S & !D))
}

test.summary <- function(subset, true.subset) {
  # Estimates properties of using a credible subgroup as a diagnostic
  # test for positive treatment effect.
  #
  # Args:
  #   subset: Logical vector indicating membership of rows of
  #           Z in credible subgroup
  #   Z: Matrix of sample drawn from predictive covariate distribution
  #   gamma: True value of predictive parameters
  #   delta: Efficacy threshold
  #
  # Returns: Vector of four estimates:
  #          Positive predictive value,
  #          Negative predictive value,
  #          Sensitivity,
  #          Specificity.
  true.positive <- true.subset
  true.negative <- !true.positive
  test.positive <- subset
  test.negative <- !subset
  
  positive.predictive.power <- sum(test.positive & true.positive) / sum(test.positive)
  negative.predictive.power <- sum(test.negative & true.negative) / sum(test.negative)
  sensitivity <- sum(test.positive & true.positive) / sum(true.positive)
  specificity <- sum(test.negative & true.negative) / sum(true.negative)
  
  return(c(pppw=positive.predictive.power,
           nppw=negative.predictive.power,
           sens=sensitivity,
           spec=specificity))
}

average.overestimate <- function(subset, real.effect, delta=0) {
  # Estimates average overestimate from using exclusive credible subgroup
  # membership to bound treatment effect.
  #
  # Args:
  #   subset: Logical vector indicating membership of rows of
  #           Z in credible subgroup
  #   Z: Matrix of sample drawn from predictive covariate distribution
  #   gamma: True value of predictive parameters
  #   delta: Efficacy threshold
  #
  # Returns: Average overestimate
  true.negative <- real.effect <= delta
  test.positive <- subset
  overestimated <- true.negative & test.positive
  overestimate  <- delta - real.effect[overestimated]
  return(mean(overestimate))
}

maximum.overestimate <- function(subset, real.effect, delta=0) {
  # Estimates maximum overestimate from using exclusive credible subgroup
  # membership to bound treatment effect.
  #
  # Args:
  #   subset: Logical vector indicating membership of rows of
  #           Z in credible subgroup
  #   Z: Matrix of sample drawn from predictive covariate distribution
  #   gamma: True value of predictive parameters
  #   delta: Efficacy threshold
  #
  # Returns: Maximum overestimate
  true.negative <- real.effect <= delta
  test.positive <- subset
  overestimated <- true.negative & test.positive
  overestimate  <- delta - real.effect[overestimated]
  return(max(overestimate))
}
