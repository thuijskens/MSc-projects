# Code for preliminary analysis
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

# Get some first statistics
head(schools)
tail(schools)
summary(schools)

# From the summary:
# 1. More urban than non-urban schools (24k vs 10k)
# 2. approx 50% male-female
# 3. 3k unclasscified working classes
# 4. far more state-funded schools (33k vs 1.5k)
# 5. Far more non-denominational schools (29k vs 5k)

##
# Top level analysis
##

# The data is structured as students nested in schools
# How many different schools are there?
schools %>% summarise(n_distinct(schoolid)) # 508

# How are the schools structured?
school_level <- schools %>% 
  group_by(schoolid, schtype, schdenom, schurban) %>% 
  summarise(med_score = median(as.numeric(score))) %>%
  summary()

# What is the distribution of students across schools?
schools %>% group_by(schoolid) %>% summarise(nrstudents = n()) %>%
  ggplot(aes(x = nrstudents)) + geom_histogram(fill = "white", colour = "black", binwidth = 10)

# Is there a diference in score for some of the explanatory variables
p_bp_cohort <- ggplot(schools, aes(x = cohort90, y = score)) + geom_boxplot() + labs(x = "Cohort", y = "Score")
p_bp_sex <- ggplot(schools, aes(x = female, y = score)) + geom_boxplot() + labs(x = "Sex", y = "Score")
p_bp_class <- ggplot(schools, aes(x = sclass, y = score)) + geom_boxplot() + labs(x = "Social class", y = "Score")
p_bp_type <- ggplot(schools, aes(x = schtype, y = score)) + geom_boxplot() + labs(x = "School type", y = "Score")
p_bp_urban <- ggplot(schools, aes(x = schurban, y = score)) + geom_boxplot() + labs(x = "School area", y = "Score")
p_bp_rel <- ggplot(schools, aes(x = schdenom, y = score)) + geom_boxplot() + labs(x = "School religion", y = "Score")

print(p_bp_cohort)
grid.arrange(p_bp_sex, p_bp_class, p_bp_type, p_bp_type, p_bp_urban, p_bp_rel, ncol = 3)

# The following variables seem to have an effect on score
# 1. Social class. (Higher class = higher score)
# 2. School type. (Independent schools score higher)
# 3. Cohort (Later cohorts score better)


# NOTE:
# There are students with a score of zero. What is going on with these students?
schools %>% filter(score == 0) %>% summary()

# There are 3158 students with a score of zero. Most of them in 1984 (1302/3158)
# Does a score of 0 mean that a student failed completely or that (s)he did not take the exam?
# If the latter: remove from the dataset


# What is the distribution of scores across schools?
# Note: plotting this for ALL schools is unreadable. So we just pick a couple to show 
# that there is a change in score distributon across different schools
# Or just plot the median scores

median_scores <- schools %>% group_by(schoolid) %>% summarise(median_score = median(as.numeric(score))) %>% arrange(desc(median_score))

# Option 1:
plotschools <- as.integer(unlist(median_scores[round(seq(from = 1, to = 508, length.out = 20)), 1]))

schools %>% filter(schoolid %in% plotschools) %>%
ggplot(aes(x = reorder(schoolid, score, FUN=median), y = score)) +
  geom_boxplot() +
  labs(x = "School number", y = "Score", title = "Boxplot of scores of students for different schools")

# Option 2:
ggplot(median_scores, aes(x = as.numeric(reorder(schoolid, median_score)), median_score)) + 
  geom_line() + 
  labs(x = "", y = "Median score", title = "Median scores for each schools arranged in ascending order")



# Is size (# students per cohort per school) of a school an expl. var?
school_sizes <- schools %>% group_by(schoolid, cohort90) %>% summarise(size = n())

median_scores %>% left_join(school_sizes) %>%
ggplot(aes(x = size, y = median_score)) + facet_wrap(~ cohort90) + geom_point()

# If we define size as the total number of all students per cohort than the plot becomes
school_sizes2 <- schools %>% group_by(schoolid) %>% summarise(size = n())
median_scores %>% left_join(school_sizes2) %>%
  ggplot(aes(x = size, y = median_score)) + geom_point()

# There is not a clear relationship visible. Large schools seem to score in a more narrow band than small schools though. 

# Cutoff point: 100.
schools_size <- school_sizes2 %>% mutate(schsize = ifelse(size >= 100, "Large", "Small")) %>% select(schoolid, schsize)
schools_size$schsize <- as.factor(schools_size$schsize)

# Join it to the original data set
schools <- schools %>% left_join(schools_size)

# Now see if there is a difference on top level
ggplot(schools, aes(x = schsize, y = score)) + geom_boxplot()

# Not a significant difference..






# Conclusion: Social class, school type and cohort seem to be good predictors of the score


# What happens if we subset the data, to look only at independent schools for example
schools %>%
  ggplot(aes(x = sclass, y = score)) +
  geom_boxplot() +
  #coord_flip() + 
  facet_wrap(~ schtype) +
  labs(x = "Class", y = "Score", title = "Scores against social class sliced by school type")

# We see different behaviour here (suggests interaction effect?)! What happens if we slice at the cohort level?
schools %>%
  ggplot(aes(x = sclass, y = score)) +
  geom_boxplot() +
  facet_wrap(~ cohort90) + 
  labs(x = "Class", y = "Score", title = "Score against social class for each cohort")

# Here we see consistent behaviour accross cohorts. Is it the same for independent schools vs state-funded schools?
schools %>%
  ggplot(aes(x = schtype, y = score)) + 
  geom_boxplot() + 
  facet_wrap(~ cohort90) +
  labs(x = "Class", y = "Score", title = "Score against school type for each cohort")

# Same kind of structure, variability gets smaller when the cohort increases
  
# There are only three types of school-level variables. School type, school area and schole denomination.
schools %>%
 ggplot(aes(x = schtype, y = score)) + 
  geom_boxplot() +
  #geom_point(aes(colour = schdenom)) + 
  #geom_jitter(aes(colour = schdenom)) +
  facet_wrap(~ schurban)

## !!
# Urban independent schools score very well, wheres town/rural independent school score very variable
# This suggests that there is an interaction effect going between schtype and schurban, since there is an overall effect of school type on score
## !!

# Updated conclusion:
# 1. Interaction effect between schtype and schurban
# 2. Interaction effect between sclass and schtype (?)
# 3. cohort90 has an effect on the score
