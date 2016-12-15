library(methods) # allows running via Rscript
library(BayesTree)

source("bayes-linear.R")

generate.y.linear <- function(sex, age, treat, phi, radius) {
  W <- cbind(1, sex, age, treat, treat * sex, treat * age)
  rnorm(nrow(W), W %*% phi)
}

generate.y.sqrt <- function(sex, age, treat, phi, radius) {
  W <- cbind(1, sex, age, treat, treat * sex,
             treat * (sqrt(age + radius) - sqrt(radius)))
  rnorm(nrow(W), W %*% phi)
}

generate.y.sigmoid <- function(sex, age, treat, phi, radius) {
  W <- cbind(1, sex, age, treat, treat * sex,
             treat * (sign(age) * abs(age) ^ (1/3)))
  rnorm(nrow(W), W %*% phi)
}

generate.y.invu <- function(sex, age, treat, phi, radius) {
  W <- cbind(1, sex, age, treat, treat * sex,
             treat * (0.5 - (age / radius) ^ 2))
  rnorm(nrow(W), W %*% phi)
}

effect.linear <- function(sex, age, gamma, radius) {
  Z <- cbind(1, sex, age)
  Z %*% gamma
}

effect.sqrt <- function(sex, age, gamma, radius) {
  Z <- cbind(1, sex, sqrt(age + radius) - sqrt(radius))
  Z %*% gamma
}

effect.sigmoid <- function(sex, age, gamma, radius) {
  Z <- cbind(1, sex, sign(age) * abs(age) ^ (1/3))
  Z %*% gamma
}

effect.invu <- function(sex, age, gamma, radius) {
  Z <- cbind(1, sex, 0.5 - (age / radius) ^ 2)
  Z %*% gamma
}

fit.linear <- function(data, test.set) {
  # Conservative prior for interactions, vague prior for others
  R  <- diag(c(10000, 10000, 10000, 10000, 1, 1))
  
  blm(y ~ trt*sex + trt*age, data=data, prior.precision=solve(R))
}

fit.bart <- function(data, test.set) {
  test.mat <- as.matrix(rbind(cbind(trt=0, test.set),
                    cbind(trt=1, test.set)))
  bart(cbind(data$trt, data$sex, data$age), data$y, x.test=test.mat,
        ndpost=500, verbose=FALSE)
}

est.effect.linear <- function(test.set, fit) {
  p      <- 3
  q      <- 3
  gamma.hat <- fit$posterior$t.mean[(p+1):(p+q)]
  Z <- cbind(1, test.set$sex, test.set$age)
  Z %*% gamma.hat
}

est.effect.bart <- function(test.set, fit) {
  est.effect.f <- function(i, matrix) {
    mean(matrix[, ncol(matrix) / 2 + i] - matrix[, i])
  }
  sapply(1:nrow(test.set), est.effect.f,
          matrix=fit$yhat.test)
}

het.test.linear <- function(fit, cred) {
  p <- 3
  q <- 3
  n <- length(fit$fitted.values)
  g.cov  <- fit$posterior$t.cov[(p+1):(p+q), (p+1):(p+q)]
  g.mean <- fit$posterior$t.mean[(p+1):(p+q)]
  het.stat <- t(g.mean[-1]) %*% solve(g.cov[-1, -1]) %*% g.mean[-1]
  het.test <- pf((n-q+1)/((q-1)*(n-1)) * het.stat, q-1, n-(q-1)) > cred
}

het.test.bart <- function(fit, cred) {
  NA
}

credsubs.linear <- function(space, fit, cred, post.size, method, delta) {
  find.credible.subgroups(space, fit, "trt",
                          delta=delta, cred=cred, df=fit$df.residual,
                          post.sample.size=post.size,
                          method=method, tolerance=0.01)
}

credsubs.bart <- function(space, fit, cred, post.size, method, delta) {
  is.in.D <- function(i, matrix, delta, cred) {
    return(mean(matrix[, ncol(matrix) / 2 + i] - matrix[, i] > delta) >= cred)
  }
  is.in.S <- function(i, matrix, delta, cred) {
    return(mean(matrix[, ncol(matrix) / 2 + i] - matrix[, i] > delta) >= 1 - cred)
  }
  in.D <- sapply(1:nrow(space), is.in.D,
                 matrix=fit$yhat.test,
                 delta=delta, cred=cred)
  in.S <- sapply(1:nrow(space), is.in.S,
                 matrix=fit$yhat.test,
                 delta=delta, cred=cred)
  credibles <- list(D=in.D, S=in.S)
}

simulate <- function(seed, n=100, delta=0, cred=0.95, method="hpd",
                     generate.y.func, effect.func, gamma=c(0, 0, 0), radius=3,
                     fit.func, est.effect.func, credsubs.func, het.test.func,
                     post.size=10000, verbose=TRUE, plots=FALSE) {
  # Simulation configuration
  set.seed(seed)
  if (verbose) {
    print(seed)
  }
  
  # Approximation of covariate space
  test.set <- expand.grid(sex=c(0, 1),
                          age=seq(-3, 3, length.out=61))
  test.mat <- as.matrix(cbind(1, test.set))
  
  # True parameters
  beta  <- c(0, 0, 0)
  phi   <- c(beta, gamma)
  
  real.effect <- effect.func(test.set$sex, test.set$age, gamma, radius)
  real.subgroup <- real.effect > delta
  
  # Generate data
  sex <- rbinom(n, 1, 0.5)
  age <- runif(n, -3, 3)
  trt <- rbinom(n, 1, 0.5)
  Z   <- cbind(1, sex, age)
  y   <- generate.y.func(sex, age, trt, phi, radius)
  
  data <- data.frame(cbind(y, trt, sex, age))
  
  # Model fit
  fit <- fit.func(data, test.set)
  
  # Model diagnostics
  est.effect <- est.effect.func(test.set, fit)
  mse.effect <- mean((est.effect - real.effect) ^ 2)
  het.test <- het.test.func(fit, cred)
  
  # Find credible subgroups
  credibles <- credsubs.func(test.set, fit, cred, post.size, method, delta)
  
  # Credible subgroup diagnostics
  csp.size     <- credible.subgroup.pair.size(credibles$D, credibles$S)
  summ.stats.D <- test.summary(credibles$D, real.subgroup)
  summ.stats.S <- test.summary(credibles$S, real.subgroup)
  avg.over     <- average.overestimate(credibles$D, real.effect, delta)
  max.over     <- max(maximum.overestimate(credibles$D, real.effect, delta), 0)
  
  # Plotting
  if (plots) {
    plot(test.set[, 2] ~ test.set[, 3], pch=ifelse(credibles$D,   19, 4), main="D")
    plot(test.set[, 2] ~ test.set[, 3], pch=ifelse(real.subgroup, 19, 4), main="Truth")
    plot(test.set[, 2] ~ test.set[, 3], pch=ifelse(credibles$S,   19, 4), main="S")
  }
  
  # Check credible subgroup coverage
  success.D <- all(real.subgroup[credibles$D]) # D is exclusive
  success.S <- all(credibles$S[real.subgroup]) # S is inclusive
  success <- success.D && success.S
  
  return(c(success.tot=success, success.D=success.D, success.S=success.S, csp.size=csp.size, #error,
           summ.stats.D, summ.stats.S, avg.over=avg.over, max.over=max.over,
           mse.effect=mse.effect, het.test=het.test))
}

##########################
### Perform Simulation ###
##########################

# Configuration
n <- 40
M <- 1000
radius <- 3
methods <- c("pb", "rcs", "hpd", "pointwise", "bart")
generators <- c(rep("linear", 6), "sqrt", "sigmoid", "invu")
models <- c(rep("linear", 4), "bart")
gammas <- matrix(c(0, 0, 0,
                   0, 0, 1,
                   0, 1, 0,
                   0, 1, 1,
                   1, 0, 0,
                   1, 1, 1,
                   0, 0, 1,
                   0, 0, 1,
                   0, 0, 1),
                 ncol=3, byrow=TRUE)
#levels <- c(0.20, 0.50, 0.80, 0.95)
levels <- c(0.80)

# Results array
tables <- array(dim=c(16,
                      length(generators),
                      length(methods),
                      length(levels),
                      M),
                dimnames=list(c("Total Coverage", "Exclusive Coverage", "Inclusive Coverage", "Credible Pair Size",
                                "Pos Pred Pow of D", "Neg Pred Pow of D", "Sensitivity of D", "Specificity of D",
                                "Pos Pred Pow of S", "Neg Pred Pow of S", "Sensitivity of S", "Specificity of S",
                                "Avg Overestimate", "Max Overestimate", "Effect MSE", "Het Test"),
                              c("Linear (0, 0, 0)", "Linear (0, 0, 1)", "Linear (0, 1, 0)",
                                "Linear (0, 1, 1)", "Linear (1, 0, 0)", "Linear (1, 1, 1)",
                                "Square Root", "Sigmoid", "Inverted U"),
                              c("PB", "RCS", "HPD", "Pointwise", "BART"),
                              #c("20%", "50%", "80%", "95%")))
                              c("80%")))

# Simulate
for (gen in 1:length(generators)) {
  print(paste0("Generator :", generators[gen], " (", gen, ")"))
  gamma <- gammas[gen, ]
  print(gamma)
  for(m in 1:length(methods)) {
    print(paste0("Method: ", methods[m], " (", m, ")"))
    method <- methods[m]
    for(l in 1:length(levels)) {
      level <- levels[l]
      print(paste("Level:", l))
      tables[, gen, m, l, ] <- sapply(1:M, simulate, n=n, cred=level,
                                    post.size=1000, delta=0, gamma=gamma,
                                    generate.y.func=get(paste0("generate.y.", generators[[gen]])),
                                    effect.func=get(paste0("effect.", generators[[gen]])),
                                    fit.func=get(paste0("fit.", models[[m]])),
                                    est.effect.func=get(paste0("est.effect.", models[[m]])),
                                    het.test.func=get(paste0("het.test.", models[[m]])),
                                    credsubs.func=get(paste0("credsubs.", models[[m]])),
                                    method=method, verbose=FALSE)
    }
  }
}

# Save results to file
save(tables, file=paste0("tables-", n, ".RData"))

# Display table used in paper
for (truth in 1:9) {
  print(dimnames(tables)[[2]][truth])
  print(format(round(t(cbind(rowMeans(tables[c(1, 4, 7, 8, 15, 16), truth, , "80%", ], na.rm=TRUE, dims=2))), 2), nsmall=2))
}

