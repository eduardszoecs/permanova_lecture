### ------------ Load data and package -----------------------------------------
### Load data and package
require(vegan)

# This downloads the data files directly from the web
werra_sp <- read.table('https://raw.githubusercontent.com/EDiLD/permanova_lecture/master/data/werra_sp.csv', 
                       sep = ',', header = TRUE, row.names = 1)
werra_env <- read.table('https://raw.githubusercontent.com/EDiLD/permanova_lecture/master/data/werra_env.csv', 
                        sep = ',', header = TRUE)

#! Note: If this does not work (probably because of https...), open the urls in your browser,
#! download them ('Save page as') and read in the downloaded csv

str(werra_sp)
str(werra_env)


### ------------ Distance matrix ----------------------------------------------
### Compute distance matrix using Bray-Curtis on double-root transformed abundances
range(werra_sp)
range(werra_sp^0.5)
range(werra_sp^0.25)

# We use a double root transformation
dist_werra <- vegdist(werra_sp^0.25, method = 'bray')



### ------------ NMDS ----------------------------------------------------------
### Run NMDS
nmds <- metaMDS(dist_werra)
### Plot NMDS
op <- ordiplot(nmds, type = 'n')

# points
cols = c('red', 'green')
points(nmds, cex = 2, pch = 16, col = cols[werra_env$position])

# decoration
ordispider(nmds, groups = werra_env$position, label = TRUE)
ordihull(nmds, groups = werra_env$position, lty = 'dotted')
legend("bottomleft", pch = 16, col = cols, legend = levels(werra_env$position))
# upstream and downstream communities are well separated (indicates an effect!)
# spread may be different (lower for upstream section)



### ------------ PERMANOVA -----------------------------------------------------
pmv <- adonis(werra_sp^0.25 ~ position, data = werra_env, 
              permutations = 999, 
              method = 'bray')
pmv
# communities differ statsitically significant between upstream/downstream
# the position explains 30.7% of variance

# plot permuted F-values
densityplot(permustats(pmv))

## But are the assumptions met?
# check assumptions visually
# on previouse plot
# graphically we would say that upstream has lower spread then downstream


### ------------ Distance based dispersion test -------------------------------
bd <- betadisper(dist_werra, werra_env$position)
bd
# also a eigenvalue based method

boxplot(bd)
# boxplot of Average distance to median shows also that there might be lower dispersion 
# for upstream sites

# F-Test
anova(bd)

# permutaion test
permutest(bd)
# We cannot find a statistically different dispersion 
# -> assumption of homogeneity is met



### ------------ SIMPER --------------------------------------------------------
sim <- simper(werra_sp^0.25, group = werra_env$position)
sim
summary(sim)
# contr :   contribution to dissimilarity between upstream and downstream
# sd    :   standard deviation of contribution (is the species response consitent?)
# ratio :   ratio between contr and sd (high ratio = high, consisten contribution)
# av.   : average abundance per groups
# cumsum:   cumulative contribution (rule of thumb : species till 70% are investigated)

# many species contribute to the difference between upstream and downstream sections
# Lype sp. decreases (from 6.29 to 2.49) and has the highest contribution
# Lasiocephala.basalis increases (from 3.24 to 7.09)
# Gammarus sp. descreases to zero (from 2.78 to 0) and also shows this consitently (high contribution to sd ratio)




### ------------ Exercise ------------------------------------------------------

# In this exercise we will use benthic community data around offshore installations. 
# The data was retrieved from the UK Benthos Database, Oil & Gas UK, 
# http://www.oilandgasuk.co.uk/knowledgecentre/uk_benthos_database.cfm.

# Benthic communities were sampled at sites with varying distances from an oil platform 
# along transects. 
# Moreover, several environmental parameters have been measured at the sampling sites.
# However for this exercise we only need the distance column of the table.

# We are interested if communities near the oil platform differ from communities apart.


# Task 1 - Read the data into R
oil_abu <- read.table('https://raw.githubusercontent.com/EDiLD/permanova_lecture/master/data/oil_abu.csv', 
                      sep = ',', header = TRUE)
oil_env <- read.table('https://raw.githubusercontent.com/EDiLD/permanova_lecture/master/data/oil_env.csv', 
                      sep = ',', header = TRUE)
# remove not needed variables from oil_env
oil_env <- oil_env[ , c('distance', 'site')]


# Task 2 - Group the sampling sites into three distance classes (0-250m, 251-750m, >750m). 
# Use the variable 'distance' in the oil_env table and 
# the cut() function to classify the distance classes.
# Add this new variable to the oil_env data.frame.

# Task 3 - The abundance data (oil_abu) consists of counts. 
# What distance measure may be appropriate for this kind of data?

# Task 4 - Run a NMDS on the community data using an appropriate distance measure. 
# Plot the results of the NMDS. 
# Display only sites and use different colors for the distance classes.

# Task 5 - Perform a PERMANOVA to test whether communities near the oil platform 
# differ from communities further apart.
# How much of the variance in the community data could be explained by the distance class?

# Task - 6 The PERMOVA might show that there is statistically significant difference 
# between different distance classes.
#  However, can we trust these results? 
# Perform a distance-based test for homogeneity of dispersion using betadisper.
