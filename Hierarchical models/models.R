# Code for model fitting
library(data.table)
library(dplyr)
library(nlme)
library(ggplot2)
library(gridExtra)

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

schools <- schools %>% filter(score != 0)

# First we fit a very basic mixed effects model, to show that
# there is indeed a need to use a multilevel model (confirming the suspicion
# that score depends on the school)

schools_base <- lme(fixed = score ~ 1 + cohort90,
                    random = score ~ 1 | schoolid,
                    data = schools)

schools_base_lm <- lm(score ~ 1 + cohort90, data = schools)

anova(schools_base, schools_base_lm)
# very significant

# Need to do residual analysis:
plot(schools_base, schoolid ~ resid(.), abline = 0, ylab = "School")
# Very bad plot, only thing that is clear from this is that the variability is still rather large.

qqnorm(schools_base, ~ resid(., type = "p"))
# Residuals of the within group-errors seem normally distributed.





# Now find better models:
# 1. Interaction effect between schtype and schurban
# 2. Interaction effect between sclass and schtype (?)
# 3. cohort90 has an effect on the score

schools_type <- lme(fixed = score ~ 1 + schtype*schurban + cohort90 + sclass*schtype,
                    random = score ~ 1 | schoolid,
                    data = schools)

schools_type_re <- lme(fixed = score ~ 1 + schtype*schurban + cohort90 + sclass*schtype,
                       random = score ~ 1 + schtype | schoolid,
                       data = schools)

anova(schools_type, schools_type_re)
# We reject the null hypothesis that schools_type is more appropriate) The random effects are significant

# Adding cohort90 to the random effects part is not appropriate. The model then becomes a three-level model
# since students are nested in cohorts and cohorts are nested in schools.

# Adding sclass to the random effects part is also not appropriate as the random effects model between-school
# variability and social class is a student-level variable.
schools_cohort <- lme(fixed =  score ~ 1 + cohort90 + sclass,
                      random = score ~ 1 | schoolid/cohort90,
                      data = schools)

schools_cohort2 <- lme(fixed =  score ~ 1 + cohort90 + sclass*schtype,
                      random = score ~ 1 | schoolid/cohort90,
                      data = schools)

# Difference in log likelihood is small

schools_cohort3 <- lme(fixed =  score ~ 1 + cohort90 + sclass + sclass*schurban,
                       random = score ~ 1 | schoolid/cohort90,
                       data = schools)

# sclass*schurban is not significant at all

schools_cohort4 <- lme(fixed =  score ~ 1 + cohort90 + sclass + schtype,
                       random = score ~ 1 | schoolid/cohort90,
                       data = schools)

# almost everything is highly significant!

