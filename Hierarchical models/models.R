# Code for model fitting
library(data.table)
library(dplyr)
library(nlme)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd('C:/Users/Thomas/Documents/MSc in Applied Statistics/Statistical Methods/Hierarchical Models')

# Load in the data
schools <- data.table::fread('schools.txt', header = TRUE, sep = ",")
schools <- as_data_frame(schools)

# Make sure the data is formatted correctly
schools$schoolid <- as.factor(schools$schoolid)

schools$cohort90 <- factor(schools$cohort90,
                           levels = c(-6, -4, -2, 0, 6, 8),
                           labels = c(1984, 1986, 1988, 1990, 1996, 1998))

schools$female <- factor(schools$female, 
                         levels = c(1, 0), 
                         labels = c("Female", "Male"))

schools$sclass <- factor(schools$sclass, 
                         levels = c(1, 2, 3, 4), 
                         labels = c("Managerial and professional",
                                    "Intermediate",
                                    "Working",
                                    "Unclassified"))

schools$schtype <- factor(schools$schtype,
                          levels = c(1, 0),
                          labels = c("Independent",
                                     "State-funded"))

schools$schurban <- factor(schools$schurban,
                           levels = c(1, 0),
                           labels = c("Urban",
                                      "Town or rural"))

schools$schdenom <- factor(schools$schdenom,
                           levels = c(1, 0),
                           labels = c("Roman Catholic",
                                      "Non-denominational"))


# First we fit an intercept-only mixed effects model, to show that
# there is indeed a need to use a multilevel model (confirming the suspicion
# that score depends on the school)

schools_base <- lme(fixed = score ~ 1,
                    random = score ~ 1 | schoolid,
                    data = schools)

# Need to do residual analysis

# Now find better models:
# 1. Interaction effect between schtype and schurban
# 2. Interaction effect between sclass and schtype (?)
# 3. cohort90 has an effect on the score

schools_type <- lme(fixed = score ~ 1 + schtype*schurban,
                    random = score ~ 1 | schoolid,
                    data = schools)

# This is already a lot better
