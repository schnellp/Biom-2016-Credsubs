library(plotrix)
library(vcd)

cs.decompose <- function(predictors, fit, treatment,
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
  
  ests <- predictors %*% post.mean
  ses <- sqrt(rowSums((predictors %*% post.cov) * predictors))
  
  return(list(ests=ests, ses=ses))
}

decomp <- cs.decompose(grid, fit.b, "TREATMENT", cred=0.80)

par(mar=c(4.1, 4.1, 4.1, 2.1))

### Estimated treatment effects ###

layout(matrix(c(1, 2,
                3, 4,
                5, 5), 
              3, 2, byrow=TRUE),
       height=c(3, 3, 1))
max.mag <- max(abs(decomp$ests))
for (i in c("F", "M")) {
  for (j in c("NON-CARRIER", "CARRIER")) {
    w <- which(
      grid[, 3] == i &
        grid[, 4] == j)
    image(t(matrix(decomp$ests[w], nrow=36, ncol=41, byrow=TRUE)),
          col=rev(c(color.scale(1:50, 0, c(1, 0.5), 0, color.spec="rgb"),
                    color.scale(1:50, c(0.5, 1), 0, 0, color.spec="rgb"))),
          zlim=c(-15, 15),
          xlab="Severity", ylab="Age",
          xaxt='n', yaxt='n')
    
    title(paste(titles.sex[i], titles.car[j]))
    axis(1, at=(seq(5, 45, by=5) - 5) / 40, labels=seq(5, 45, by=5))
    axis(2, at=(seq(55, 90, by=5) - 55) / 35, labels=seq(55, 90, by=5))
  }
}

par(mgp = c(0.2, 1, 0))
legend_image <- as.raster(matrix(testcol<-c(color.scale(1:50, c(1, 0.5), 0, 0, color.spec="rgb"),
                                            color.scale(1:50, 0, c(0.5, 1), 0, color.spec="rgb")), nrow=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '', cex.lab=1.5)
mtext('Estimated Treatment Effect', side=3, line=1)
axis(3, at=c(0, 2), labels=c(-15, 15))
rasterImage(legend_image, 0, -30, 2, 30)
par(mgp = c(3, 1, 0))


### Treatment effect standard errors ###

layout(matrix(c(1, 2,
                3, 4,
                5, 5), 
              3, 2, byrow=TRUE),
       height=c(3, 3, 1))
max.se <- max(decomp$ses)
for (i in c("F", "M")) {
  for (j in c("NON-CARRIER", "CARRIER")) {
    w <- which(
      grid[, 3] == i &
        grid[, 4] == j)
    image(t(matrix(decomp$ses[w], nrow=36, ncol=41, byrow=TRUE)),
          col=rev(color.scale(1:100,c(0, 0.75), 1,1,color.spec="hsv")),
          zlim=c(1.9, 6.7),
          xlab="Severity", ylab="Age",
          xaxt='n', yaxt='n')
    
    title(paste(titles.sex[i], titles.car[j]))
    axis(1, at=(seq(5, 45, by=5) - 5) / 40, labels=seq(5, 45, by=5))
    axis(2, at=(seq(55, 90, by=5) - 55) / 35, labels=seq(55, 90, by=5))
  }
}

par(mgp = c(0.2, 1, 0))
legend_image <- as.raster(matrix(testcol<-color.scale(100:1,c(0, 0.75), 1,1,color.spec="hsv"), nrow=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '', cex.lab=1.5)
mtext('SE(Treatment Effect)', side=3, line=1)
# rasterImage(legend_image, 0, 0.25, 2, 0.75)
rasterImage(legend_image, 0, -30, 2, 30)
axis(3, at=c(0, 2), labels=c(2, 7))
par(mgp = c(3, 1, 0))


### Observation plot ###

layout(matrix(c(1, 2,
                3, 4,
                5, 5),
              3, 2, byrow=TRUE),
       height=c(3, 3, 1))

for (i in c("F", "M")) {
  for (j in c("NON-CARRIER", "CARRIER")) {
    w <- which(
      data$SEX == i &
        data$CARRIER == j)
    
    
        plot(data$AGE[w] ~ data$SEVERITY[w],
             xlim=scale(c(5, 45), center=attr(data$SEVERITY, "scaled:center"),
                        scale=attr(data$SEVERITY, "scaled:scale")),
             ylim=scale(c(55, 90), center=attr(data$AGE, "scaled:center"),
                        scale=attr(data$AGE, "scaled:scale")),
             xaxt='n', yaxt='n',
             xlab="Severity", ylab="Age",
             cex=2,
             pch=ifelse(data$TREAT[w] == "low", 3, 1))
        
        title(paste(titles.sex[i], titles.car[j]))
        axis(1, at=scale(seq(5, 45, by=5), center=attr(data$SEVERITY, "scaled:center"),
                         scale=attr(data$SEVERITY, "scaled:scale")), labels=seq(5, 45, by=5))
        axis(2, at=scale(seq(55, 90, by=5), center=attr(data$AGE, "scaled:center"),
                         scale=attr(data$AGE, "scaled:scale")), labels=seq(55, 90, by=5))
      
    }
  }

scale(seq(5, 45, by=5), center=attr(data$SEVERITY, "scaled:center"),
      scale=attr(data$SEVERITY, "scaled:scale"))

plot.new()
g = grid_legend(0, 0,
                pch=c(1),
                labels=c("Control"), 
                frame=FALSE, draw=FALSE)
grid.draw(grobTree(g, vp=viewport(x=0.3, y=0.1)))

g = grid_legend(0, 0,
                pch=c(3),
                labels=c("Treatment"), 
                frame=FALSE, draw=FALSE)
grid.draw(grobTree(g, vp=viewport(x=0.7, y=0.1)))

par(mar=c(5.1, 4.1, 4.1, 2.1))
