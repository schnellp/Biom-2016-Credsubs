library(vcd)
source("bayes-linear.R")

set.seed(1)

data <- read.csv("alzheimers.csv", header=TRUE)

# OPTIONAL: set proper reference levels for factors
# Makes output more intuitive
data$TREATMENT <- relevel(data$TREATMENT, ref="placebo")
data$CARRIER   <- relevel(data$CARRIER,   ref="NON-CARRIER")

# OPTIONAL: standardize continuous covariates
# Better numerical accuracy, and more importantly,
# Easier to specify balanced priors
data$SEVERITY <- scale(data$SEVERITY)
data$AGE <- scale(data$AGE)

# Get design matrix to build conservative prior
design <- blm(CHANGE ~ 
                TREATMENT*SEVERITY + TREATMENT*AGE +
                TREATMENT*SEX*CARRIER,
              data=data, design=TRUE)

# Get column names that start with "TREATMENTlow:"
# These are treatment-covariate interactions
interaction.cols <- grep("TREATMENTlow:", colnames(design))

# Construct conservative prior for regression parameter precision
precision <- rep(1E-4, ncol(design))
names(precision) <- colnames(design)
precision[interaction.cols] <- 1

# Fit model
fit.b <- blm(CHANGE ~
               TREATMENT*SEVERITY + TREATMENT*AGE +
               TREATMENT*SEX*CARRIER,
             prior.precision=precision,
             data=data)
summary(fit.b)

# Construct grid on predictive covariate region of interest
grid <- expand.grid(SEVERITY=5:45,
                    AGE=55:90,
                    SEX=factor(c("F", "M")),
                    CARRIER=factor(c("NON-CARRIER", "CARRIER"))
)
# Transform to match original standardization
grid$SEVERITY <- scale(grid$SEVERITY,
                       center=attr(data$SEVERITY, "scaled:center"),
                       scale=attr(data$SEVERITY, "scaled:scale"))
grid$AGE <- scale(grid$AGE,
                  center=attr(data$AGE, "scaled:center"),
                  scale=attr(data$AGE, "scaled:scale"))

# Find credible subgroups
credibles <- find.credible.subgroups(grid, fit.b, "TREATMENT", cred=0.80)
credibles2 <- find.credible.subgroups(grid, fit.b, "TREATMENT", cred=0.50, delta=2)

################
### Plotting ###
################


par(mar=c(4.1, 4.1, 4.1, 2.1))

offset.plot <- TRUE
in.color <- TRUE

titles.sex <- c("M"="Male", "F"="Female")
titles.car <- c("NON-CARRIER"="Non-Carriers", "CARRIER"="Carriers")

in.d <- credibles$D
in.s <- credibles$S

if (offset.plot) {
  layout(matrix(c(1, 2,
                  3, 4,
                  5, 5),
                3, 2, byrow=TRUE),
         height=c(3, 3, 1))
} else {
  layout(matrix(c(1, 2,
                  3, 4), 
                2, 2, byrow=TRUE),
         width=c(1, 1))
}

for (i in c("F", "M")) {
  for (j in c("NON-CARRIER", "CARRIER")) {
    w <- which(
      grid[, 3] == i &
        grid[, 4] == j)
    image(t(matrix(in.d[w] + in.s[w], nrow=36, ncol=41, byrow=TRUE)),
          col=c("red", "yellow", "green"),
          zlim=c(0, 2),
          xlab="Severity", ylab="Age",
          xaxt='n', yaxt='n')
    
    title(paste(titles.sex[i], titles.car[j]))
    axis(1, at=(seq(5, 45, by=5) - 5) / 40, labels=seq(5, 45, by=5))
    axis(2, at=(seq(55, 90, by=5) - 55) / 35, labels=seq(55, 90, by=5))
    box()
  }
}

if (offset.plot) {
  plot.new()
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("D"),
                  col=c("black"),
                  gp=gpar(fill = c("green")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.225, y=0.1)))
  
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("S remove D"), 
                  col=c("black"),
                  gp=gpar(fill = c("yellow")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.425, y=0.1)))
  
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("S-complement"), 
                  col=c("black"),
                  gp=gpar(fill = c("red")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.725, y=0.1)))
}

in.d <- credibles2$D
in.s <- credibles2$S

if (offset.plot) {
  layout(matrix(c(1, 2,
                  3, 4,
                  5, 5),
                3, 2, byrow=TRUE),
         height=c(3, 3, 1))
} else {
  layout(matrix(c(1, 2,
                  3, 4), 
                2, 2, byrow=TRUE),
         width=c(1, 1))
}

for (i in c("F", "M")) {
  for (j in c("NON-CARRIER", "CARRIER")) {
    w <- which(
      grid[, 3] == i &
        grid[, 4] == j)
    image(t(matrix(in.d[w] + in.s[w], nrow=36, ncol=41, byrow=TRUE)),
          col=c("red", "yellow", "green"),
          zlim=c(0, 2),
          xlab="Severity", ylab="Age",
          xaxt='n', yaxt='n')
    
    title(paste(titles.sex[i], titles.car[j]))
    axis(1, at=(seq(5, 45, by=5) - 5) / 40, labels=seq(5, 45, by=5))
    axis(2, at=(seq(55, 90, by=5) - 55) / 35, labels=seq(55, 90, by=5))
    box()
  }
}

if (offset.plot) {
  plot.new()
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("D"),
                  col=c("black"),
                  gp=gpar(fill = c("green")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.225, y=0.1)))
  
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("S remove D"), 
                  col=c("black"),
                  gp=gpar(fill = c("yellow")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.425, y=0.1)))
  
  g = grid_legend(0, 0,
                  pch=c(22),
                  labels=c("S-complement"), 
                  col=c("black"),
                  gp=gpar(fill = c("red")),
                  frame=FALSE, draw=FALSE)
  grid.draw(grobTree(g, vp=viewport(x=0.725, y=0.1)))
}

par(mar=c(5.1, 4.1, 4.1, 2.1))
