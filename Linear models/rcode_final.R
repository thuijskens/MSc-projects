# We will need the following packages
library(reshape2)
library(dplyr)
library(ggplot2)
library(GGally)
library(corrplot)
library(corpcor)
library(GGplots)
library(leaps)
library(penalized)
library(glmnet)
library(ggpairs)
library(boot)

# Load data in R
boco <- read.table("C:/Users/Thomas/Documents/MSc in Applied Statistics/Statistical Methods/Practical 4/boco.txt", header = TRUE)

# Inspect data
head(boco)
summary(boco)
# No missing values, no factors. That's good.

# Subset the data so we only get LL, A, W, H and C
boco.vars <- boco[,c("LL", "A", "W" , "C", "H")]
names(boco.vars) <- c("Legs.LM", "Age", "Weight", "Circumference", "Height")

# Scatter plot matrix
pdf('scatterplot.pdf', width = 12, height = 10, paper = "special")
ggpairs(boco.vars, upper = list(continuous = "points"), axisLabels = "show", diag = list(continuous = "density"), title = "Matrix of scatter plots")
dev.off()

# Correlation matrix plots
p.cor <- -solve(cor(boco.vars))
diag(p.cor) <- -diag(p.cor)
p.cor <- cov2cor(p.cor)

# Seperate plots

par(mfrow=c(1,1))
pdf("corrplot.pdf")
corrplot(cor(boco.vars), tl.col = "black", mar=c(0,0,1,0))
title("Correlation matrix plot", line = -2)
dev.off()
# Since all the variables are heavily correlated with each other, this analysis does not help us a lot. 
# Therefore, we need to use partial correlations, where we remove the effect of the other variables

pdf("pcorrplot.pdf")
corrplot(p.cor, tl.col = "black", mar=c(0,0,1,0))
title("Partial correlation matrix plot", line = -2)
dev.off()
# A lot is going on in this plot. We see:
# 1. W and C are correlated with LL, when controlled for each other
# 2. W and C are correlated.
# 4. There is a small correlation between A and LL

# Reset the names, for easy reference
names(boco.vars) <- c("LL", "A", "W", "C", "H")

# Fit a first linear model to see what is going on
boco.fit <- lm(LL ~ A + W + C + H, data = boco)

# Get a summary of the model and adjust p-values
summary(boco.fit)
p.adjust(summary(boco.fit)$coefficients[,4], method = "fdr")

# Get diagnostic plots. Autoplot.lm is implemented by GGplots, custom made package, implementing plot.lm for ggplot2
pdf("diag_plots.pdf", height = 10, width = 10, paper = "special")
autoplot(boco.fit)
dev.off()

# part of the ANOVA Analysis
anova(lm(LL ~ A + W + H + C, data = boco))
anova(lm(LL ~ A + W + C + H, data = boco))
anova(lm(LL ~ W + C + A + H, data = boco))
anova(lm(LL ~ W + H + A + C, data = boco))
anova(lm(LL ~ H + A + W + C, data = boco))

# Calculate VIF
vif(boco.fit)

# We will use BIC as a model selection tool, since we want to explain LL with a small model, so we should penalize "large" models harder. 
bic.fit <- step(lm(LL ~ 1, data = boco), trace = TRUE, scope = ~ A + W + C + H, k = log(nrow(boco)))
# The model with the lowest BIC is LL ~ W + C. Get diagnostic plots
pdf("bic_diag_plots.pdf", height = 10, width = 10, paper = "special")
autoplot(bic.fit)
dev.off()

# Since some predictors are highly correlated, we might want to use Ridge regression as an alternative (Ridge is robust against multicollinearity)
resp.var <- boco[, "LL"]
pred.var <- boco[, c("A", "W", "C", "H")]

# With penalized package

lambda.opt <- optL2(response = resp.var, penalized = pred.var, standardize = T)
prof <- profL2(response = resp.var, penalized = pred.var, standardize = T, minlambda2 = 0.00001, maxlambda2 = 1, plot = T)

# Plot cross-validated likelihood
df <- data.frame(log.lambda = log(prof$lambda), cvl = prof$cvl)
opt.coords <- c(df[which(df[,2] == max(df[,2])), 1], max(df[,2]))

p1 <- ggplot(df, aes(x = log.lambda, y = cvl)) + 
  geom_line(col = "blue") +
  geom_vline(xintercept = opt.coords[1], linetype = 2) +
  scale_x_continuous(limits = c(-5, -2.5)) + 
  scale_y_continuous(limits = c(-164.3615, -164.356)) +
  labs(x = expression(log(lambda)), y = "Cross-validated likelihood") 

p2 <- ggplot(df, aes(x = log.lambda, y = cvl)) + 
  geom_line(col = "blue") +  
  labs(x = expression(log(lambda)), y = "Cross-validated likelihood")

pdf("ridge_plot2.pdf", width = 7, height = 5)
do.call(grid.arrange, c(list(p2,p1), main = paste("Cross-validated likelihood vs.", expression(log(lambda))), ncol = 2))
dev.off()

pen.fit <- penalized(response = resp.var, penalized = pred.var,lambda2 = lambda.opt$lambda, standardize = T)

# Ridge regression with smaller model

lambda.opt.small <- optL2(response = resp.var, penalized = pred.var[,c("C", "W")], standardize = T)
pen.fit.small <- penalized(response = resp.var, penalized = pred.var,lambda2 = lambda.opt.small$lambda, standardize = T)

# Leave-one-out cross-validation
library(boot)
boco.fit.big <- glm(LL ~ A + W + C + H, data = boco)
boco.fit.small <- glm(LL ~ W + C, data = boco)
cv.glm(boco.vars, boco.fit.big)$delta[2]
cv.glm(boco.vars, boco.fit.small)$delta[2]

set.seed(12345)
nr <- 50
mse <- matrix(0, nr, 2)
for(i in 1:nr) {
  mse[i,1] <- cv.glm(boco.vars, boco.fit.big, K =10)$delta[2]
  mse[i,2] <- cv.glm(boco.vars, boco.fit.small, K = 10)$delta[2]
}

apply(mse, 2, mean)