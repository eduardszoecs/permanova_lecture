###############################################################################
### Load data & packages
require(vegan)
data(pyrifos)
dim(pyrifos)
# 132 observations, 178 species
# data are already log-transformed
head(pyrifos[ , 1:10])


# environmental variables
ditch <- gl(12, 1, length = 132)
week <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12,15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44,0.9, 0, 6), 11))
pyrifos_env <- data.frame(ditch, week, dose)
str(pyrifos_env)


################################################################################
### PRC
mod_prc <- prc(pyrifos, treatment = dose, time = week)

## Numerical output
mod_prc
# 33.5% explained by treatment + treatment x time
# 21.9% explained by time



# The first RDA axis has an eigenvalue of 25.3.
mod_prc$CCA$eig[1]
#If we divide this
# eigenvalue by the sum of all eigenvalues, we get the proportion
# of explained variance which is displayed on the first axis1:
mod_prc$CCA$eig[1] / sum(mod_prc$CCA$eig) * 100

sum_mod_prc <- summary(mod_prc, scaling = 1)
sum_mod_prc
# species & site scores

## Plot
plot(mod_prc, scaling = 1)
# cluttered species scores
# display only species with |scores| > 1
plot(mod_prc, select = abs(sum_mod_prc$sp) > 1, scaling = 1)



## Significance tests
# Observations from a experimental ditch are not
# independent, since the same ditch was measured repeatedly during
# the experiment. We have to take this into account:
#   each ditch represents a time-series.
# We will permute the whole series of one ditch, keeping the temporal order.
# To setup such a permutation scheme we use the permute package,
# which is automatically loaded with vegan :
control <- how(plots = Plots(strata = ditch, type = "free"),
                within = Within(type = "none"), nperm = 99)
# This sets up our permutation scheme:
# plots   : We will permute ditches, without any restrictions.
# within  : But within one ditch there will be no permutations.
# nperm   : We request 99 permutations.

# Permutation tests using this scheme can be performed using anova:
set.seed(1234)
# Test of first eigenvalue
anova(mod_prc, permutations = control, first = TRUE)
# same as testing significance of axes
# but this take smuch longer, as all axes are tested
# anova(mod_prc, permutations = control, by = "axis")

# Test of terms
anova(mod_prc, permutations = control, by = "terms")

# Overall test
anova(mod_prc, permutations = control)



################################################################################
# We may be interested in recovery and reponses per week
# Idea is to run a RDA for every week and perform a permutation test
# for-loop computes for every week a RDA and a permutation test.
rdas <- NULL
for (i in levels(week)) {
  # subset data
  take_spec <- pyrifos[week == i, ]
  take_dose <- dose[week == i]
  # RDA
  rdas[[i]]$rda <- rda(take_spec ~ take_dose)
  # perm
  rdas[[i]]$anova <- anova(rdas[[i]]$rda, by = 'terms', permutations = 99) #99 permutations (to safe computation time)
}

str(rdas)
length(rdas)
# rdas is a big list

# we can extract the information we need by walking through the list
# here we extract the p-values from the permutation test
rdas[[1]]$
rdas[[1]]$anova
sapply(rdas, function(x) x$anova[1, 4])

# Note these results differ from
# Van den Brink, P. J. & Braak, C. J. F. . Principal response curves: Analysis of time-dependent multivariate responses of biological community to stress. Environmental Toxicology and Chemistry 18, 138–148 (1999).
# They used a regession like analysis (dose a continuous variable),
# we used a anova like analysis (dose a categorical )

# A package with some convenience function is currently in preparation!


#################################################################################
# NOEC
# Besides the overall significance of treatment, we can also look which treatments
# differ from control in order to get a no-observed-effect Concentration (NOEC) [=the concentration below the lowest significant concentration].
# Testing by permutations fails here, because there are not enough unique permutations
# (we have only 2 treated and 4 controls per sampling date). T

# One sultion is to break it down to a univariate problem.
# A PCA is performed and scores on first axis are compared using univariate techniques
# Therefore they applied a Williams Test [1] on the first principle component of a PCA on each sampling date.

df <- data.frame(dose = dose, week = week)
# package for multiple comparisons
require(multcomp)
out <- NULL

# loop through time, compute PCA, extract scores and do Dunnett-Test
for (i in levels(week)) {
  # PCA
  take_spec <- pyrifos[week == i, ]
  pca <- rda(take_spec)
  # scores of first principal component
  pca_scores <- scores(pca, display = "sites", choices = 1, scaling = 1)
  mod_aov <- aov(pca_scores ~ dose, data = df[week == i, ])
  out[[i]] <- summary(glht(mod_aov,
                                 alternative = "t",
                                 linfct = mcp(dose = "Dunnett")))
}
# extract p-values
result <- lapply(out, function(x) data.frame(pval = x$test$pvalues, sig = x$test$pvalues < 0.05))
# shows the results of Dunnett-Test on PCA-scores for week 1:
result[["1"]]
result[["0.1"]]

## Usually a Trend-test is performed (Williams test) - change mcp(dose = "Dunnett")
# to mcp(dose = "Williams ") for Williams type contrasts.
