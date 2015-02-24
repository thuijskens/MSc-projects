# R code for the assessed practical for Logistic regression
library(ggplot2) # for nice graphics
library(RColorBrewer) # for nice graphics
library(gridExtra) # for nice graphics
library(corrplot) # for (partial) corelation matrix plot
library(penalized) # for Ridge and LASSO regression
library(caret) # for confusion matrix functions
library(logistf) # for penalized logistic regression model
library(dplyr) # for easy reshaping of data
library(xtable) # convert tables to LaTeX

# Load helper functions for computation of ROC curve
source('helperROCfunctions.R')

# Read in the data
asthma <- read.table('asthma.txt', header = TRUE)

# Make sure 'control' is the reference value for the factor
asthma$CASE <- relevel(asthma$CASE, ref = 'control')

# Inspect the data
head(asthma[,1:5])
summary(asthma[,1:5])
sapply(asthma[,1:5], unique)
summary(asthma$CASE) # 177 control vs 23 case. Unbalanced sample

# Split data in training and test set
# Set up parameters that define the split
n <- nrow(asthma)
ntraining <- 160
ntest <- n - ntraining
ntest_control <- 35
ntest_case <- ntest - ntest_control

# training set: 142 control vs 18 case. 18/142 ~= 0.14
# test set: 35 control vs 5 case. 5/35 ~=0.13
# So the two sets are balanced the same approximately

#Split the data
set.seed(123)
test_control <- sample(which(asthma$CASE == "control"), ntest_control)
test_case <- sample(which(asthma$CASE == 'case'), ntest_case)
test <- c(test_control, test_case)

asthma_test <- asthma[test,]
asthma_training <- asthma[-test,]

#####
# For the identification of predictors, we use only the training set
#####

# For each gene, we do a two sample t-test to see if they have a 
# statistical significant effect on the response

ngenes <- ncol(asthma_training) - 1
pvals <- numeric(ngenes)

for(i in seq_len(ngenes)) {
  # Form the formula for the t.test command  
  strformula <- paste0("SNP", i, " ~ CASE")
  
  # Compute the p-value for the t.test
  pvals[i] <- t.test(formula = as.formula(strformula), data = asthma_training[, c(1,i+1)])$p.value
}

# We need to do a multiplicity adjustment
pvals <- p.adjust(pvals, method = "fdr")

# Plot these p-values
pval_data <- data.frame(gene = 1:ngenes, 
                        genestr = paste("SNP", 1:ngenes), 
                        logpvals = -log10(pvals), 
                        significant = -log10(pvals) > 8)
pval_data_annotate <- filter(pval_data, significant == 1)

p.pvals <- pval_data %>% top_n(n = 10, wt = logpvals) %>%
  ggplot(aes(x = reorder(genestr, logpvals), y = logpvals)) +
  geom_point(aes(col = significant)) + 
  geom_hline(yintercept = 8, linetype = 2) +
  scale_color_manual(values = c("#377eb8", "#e41a1c")) +
  scale_y_continuous(breaks = seq(0, 20, 2), limits = c(0,20)) + 
  #scale_color_brewer(palette = "Set1"[2:1]) +
  labs(x = "", y = expression(-log[10](p)), title = "P-values of the two sample t-test") +
  theme(legend.position="none") +
  coord_flip() +
  theme(axis.text.y = element_text(hjust=0)) 

# We identify SNP85, SNP255 and SNP42 as significant genes

# We now try to identify significant genes based on another method
linmod_pvals <- numeric(ngenes)
dev <- numeric(ngenes)

for(i in seq_len(ngenes)) {
  # Form the formula for the glm command  
  strformula <- paste0("CASE ~ ", "SNP", i)
  
  # Fit logistic regression
  mod <- glm(formula = as.formula(strformula), 
             data = asthma_training[, c(1,i + 1)], 
             family = binomial(link = "logit"))
  
  # Do the deviance test
  linmod_pvals[i] <- 1 - pchisq(mod$null.deviance - mod$deviance, df = 1)
  dev[i] <- mod$null.deviance - mod$deviance
}

# We need to do a multiplicity adjustment
linmod_pvals <- p.adjust(linmod_pvals, method = "fdr")

# Plot these p-values
pval_data_mod <- data.frame(gene = 1:ngenes, 
                            genestr = paste("SNP", 1:ngenes), 
                            logpvals = -log10(linmod_pvals),
                            significant = -log10(linmod_pvals) > 3) #!!!!
pval_data_mod_annotate <- filter(pval_data_mod, significant == 1)

p.pvals_mod <- pval_data_mod %>% filter(rank(desc(logpvals), ties.method = "first") <= 10) %>%
  ggplot(aes(x = reorder(genestr, logpvals), y = logpvals)) +
  geom_point(aes(col = significant)) + 
  geom_hline(yintercept = 8, linetype = 2) +
  scale_color_manual(values = c("#377eb8", "#e41a1c")) +
  scale_y_continuous(breaks = seq(0, 22, 2), limits = c(0,20)) + 
  #scale_color_brewer(palette = "Set1"[2:1]) +
  labs(x = "", y = expression(-log[10](p)), title = "P-values of the deviance test") +
  theme(legend.position="none") +
  coord_flip() +
  theme(axis.text.y = element_text(hjust=0)) 

# Save the plots
pdf('pvals_t.pdf', height = 4, width = 4)
print(p.pvals)
dev.off()

pdf('pvals_dev.pdf', height = 4, width = 4)
print(p.pvals_mod)
dev.off()

# We find the same set of explanatory variables has high values of -log10(p)

# We fit a penalized logistic regression model to the full training dataset to see if there are significant predictors

# Find optimal lambda 
set.seed(12345)
opt.lambda2 <- optL2(response = asthma_training$CASE,
                     penalized = asthma_training[, -1],
                     model = "logistic",
                     fold = 10)
# 10-fold CV
# lambda = 31.65318
# cvl = -53.20332

# Fit ridge regression with optimal lambda
asthma_full_ridge <- penalized(response = asthma_training$CASE,
                               penalized = asthma_training[, -1],
                               model = "logistic",
                               lambda2 = opt.lambda2$lambda)

# Find optimal lambda for LASSO model
# The seed must be equal so that LASSO uses the same folds for CV as ridge
# Then the two can be compared correctly
set.seed(12345)
opt.lambda1 <- optL1(response = asthma_training$CASE,
                     penalized = asthma_training[, -1],
                     model = "logistic",
                     fold = 10)
# 10-fold CV
# lambda = 3.120683
# cvl = -39.75176

# Fit a LASSO model with optimal lambda
asthma_full_lasso <- penalized(response = asthma_training$CASE,
                               penalized = asthma_training[, -1],
                               model = "logistic",
                               lambda1 = opt.lambda1$lambda)

# LASSO might actually be better than Ridge regression in this case
# since we identify only three significant predictors in the preliminary analysis

# Get LASSO model output in LaTeX table form
xtable(as.data.frame(coef(asthma_full_lasso)))

# Plot coefficients of Ridge regression
ridge_coefs <- coef(asthma_full_ridge)
names(ridge_coefs) <- NULL
significant <- rep(FALSE, 500)
significant[c(42,85,255)] <- TRUE

ridge_coefs_df <- data.frame(SNP = 1:500,
                             Estimate = ridge_coefs[-1],
                             Significant = significant)

pdf('ridge_coefs.pdf', width = 7, height = 5)
ggplot(ridge_coefs_df, aes(x = SNP, y = Estimate)) +
  geom_point(aes(col = Significant)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("#377eb8", "#e41a1c")) +
  theme(legend.position="none")  +
  labs(x = "SNP", y = "Estimate", title = "Coefficient size for each SNP")
dev.off()

# We now focus on the dataset with the significant predictors 
asthma_training_subset <- asthma_training[, c(1, 86, 256, 43)]

# Correlation analysis
c <- cor(x = asthma_training_subset[, -1])
p.cor <- -solve(c)
diag(p.cor) <- -diag(p.cor)
p.cor <- cov2cor(p.cor)

corrplot(c, tl.col = "black", mar=c(0,0,1,0))
title("Correlation matrix plot", line = -2)

corrplot(p.cor, tl.col = "black", mar=c(0,0,1,0))
title("Partial correlation matrix plot", line = -2)

# These predictors are not correlated to eachother at all. A simple logistic regression might work well enough
asthma_logit <- step(glm(CASE ~ 1, 
                         data = asthma_training_subset,
                         family = binomial(link = "logit")),
                     trace = TRUE,
                     scope = ~ SNP85 + SNP255 + SNP42,
                     k = log(nrow(asthma_training_subset)))

summary(asthma_logit)
# Near perfect seperation 

asthma_subset %>% filter(CASE == 'case') %>% sapply(unique)

# From this we see that all the case observations have -1 for SNP85, SNP255 and SNP42
# Therefore, the logit model with just these three explanatory variables can not be used

# We use an alternative penalization method, called Firth's correction. The method is implemented in the R package 'logistf'
firth_data$CASE <- with(firth_data, ifelse(CASE == 'case', 1, 0))
asthma_logit_firth <- logistf(formula = CASE ~ SNP85 + SNP255 + SNP42,
                              data = asthma_training_subset)

print(asthma_logit_firth)

# Computation of log odds and confidence interval
cov_firth <- vcov(asthma_logit_firth)
odds_data <- model.matrix(object = CASE ~ SNP85 + SNP255 + SNP42, 
                           data = asthma_training_subset)
odds_fitted <- (odds_data %*% coef(asthma_logit_firth))[,,drop = TRUE]
odds_se <- diag(sqrt(odds_data %*% cov_firth %*% t(odds_data)))
odds_lower <- odds_fitted + qnorm(0.025) * odds_se
odds_upper <- odds_fitted + qnorm(0.975) * odds_se
odds_df <- data.frame(n = 1:dim(asthma_training_subset)[1],
                      odds = odds_fitted,
                      lower = odds_lower,
                      upper = odds_upper,
                      cover = odds_fitted >= 0)

# Plot the log odds
pdf('firth_logodds.pdf', width = 9, height = 4)
ggplot(odds_df, aes(x = reorder(n, odds), y = odds)) +
  geom_point(aes(col = cover)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(breaks = seq(from = -30, to = 6, by = 2)) +
  scale_color_manual(values = c("#2ca25f", "red"), labels = c("Control", "Case"), name = "Classification") +
  labs(y = expression(paste(log, " ", frac(pi, 1 - pi))), title = "Log-odds and confidence intervals for fitted values") +
  theme_classic() +
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        axis.title.y = element_text(angle = 0),
        legend.position = c(0.85, 0.1)) 
dev.off()

# Compare models on the training set:
fitData(asthma_full_ridge)$stats
fitData(asthma_full_lasso)$stats
fitData(asthma_logit_firth)$stats

thr <- seq(0, 1, 0.01)

# On the test set
roc_ridge_test <- constructROC(asthma_full_ridge, thr = thr)
roc_lasso_test <- constructROC(asthma_full_lasso, thr = thr)
roc_firth_test <- constructROC(asthma_logit_firth, thr = thr)

rocdata <- rbind(data.frame(roc_ridge_test, mod = "Ridge"), 
                  data.frame(roc_lasso_test, mod = "LASSO"), 
                  data.frame(roc_firth_test, mod = "Firth"))
rocdata <- arrange(rocdata2, mod, x, y)

pdf('roc_cv.pdf', width = 8, height = 7)
ggplot(rocdata, aes(x = x, y = y)) + 
  geom_line(aes(col = mod)) + 
  #geom_point(aes(col = mod)) + 
  geom_abline(slope = 1, linetype = 2) + 
  guides(col = guide_legend(title = "Model")) +
  theme(legend.position = "top") +
  #scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  #geom_line(data = data.frame(a = c(0, 1), b = c(0, 1)), aes(a, b, col = "black"), linetype = 2) +
  labs(x = "False positive rate", y = "True positive rate", title = "ROC curve of each model based on the test set")
dev.off()

# Sensitivity = # people predicted having asthma / # people actually having asthma (=TPR)
# Specificity = # people predicted to not have asthma / # people actually not having asthma (= 1 - FPR)

# So these models have a high sensitivity, but low specificity
