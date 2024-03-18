
## Load required packages 
library(phytools)
library(ggrepel) # Labelling plots with species names
library(ggplot2) # Plotting
library(dplyr) # Data wrangling 
library(nlme) # GLS analysis
library(caper) # PGLS 
library(ggpubr) # Combining plots

#################################################################################
#################################################################################

## Preparing the dataframes and trees

## Import the dataset and tree 
bop.tree <- read.nexus("Ligon.et.al._UltrametricTree") # Import the tree as .nexus
plot(bop.tree) # Quick look at the tree
summary(bop.tree)

data <- read.csv("MacGillavry.et.al._Brain.Data.csv") # Import the dataframe as .csv

## Log-transform continuous variables 
data$logECV <- log10(data$ECV)
data$logBehav.Richness <- log10(data$behavioral.richness.ligon)
data$logBehav.Diversity <- log10(data$behavioral.diversity.ligon)
data$logAcoust.Richness <- log10(data$acoustic.richness.ligon)
data$logAcoust.Diversity <- log10(data$acoustic.diversity.ligon)
data$logComplexity.Fux <- log10(data$complexity.fuxjager)
data$logMass <- log10(data$mass)

head(data, 23) # Quick look at the data
max(data$ECV)
min(data$ECV)
mean(data$ECV)
sd(data$ECV)

## Make a new dataframe including only 'Core' BoPs
data.poly <- filter(data, mating.system == "poly") 

## Create a tree including only species in our dataset 

## First, create an object containing all the species we have data for: 
data.species <- c("Astrapia_mayeri",
                  "Astrapia_rothschildi",
                  "Cicinnurus_magnificus",
                  "Cicinnurus_regius",
                  "Cicinnurus_respublica",
                  "Lophorina_superba",
                  "Paradisea_apoda",
                  "Paradisea_minor",
                  "Paradisea_rubra",
                  "Paradisea_raggiana",
                  "Parotia_lawesii",
                  "Pteridophora_alberti",
                  "Ptiloris_magnificus",
                  "Ptiloris_paradiseus",
                  "Ptiloris_victoriae",
                  "Seleucidis_melanoleucus",
                  "Semioptera_wallacii",
                  "Lycocorax_pyrrhopterus",
                  "Manucodia_ater",
                  "Manucodia_chalybatus",
                  "Manucodia_comrii",
                  "Manucodia_jobiensis",
                  "Phonygammus_keraudrenii")

## Drop everything except the birds-of-paradise. 
bop.tree <- keep.tip(bop.tree, data.species)
plot(bop.tree)

## Then drop every species except the core birds of paradise. 
core.bop.tree <- drop.tip(bop.tree, c(
  "Lycocorax_pyrrhopterus", 
  "Manucodia_ater", 
  "Manucodia_chalybatus", 
  "Manucodia_comrii", 
  "Manucodia_jobiensis", 
  "Phonygammus_keraudrenii"))

plot(core.bop.tree) # We will use this tree for the bivariate regressions with complexity metrics

## Create comparative datasets for caper

## Specify row names so the data matches the tree
rownames(data) <- data$species

## Create comparative dataset for caper 
comparative.data <- comparative.data(phy = bop.tree, 
                                     data = data, 
                                     names.col = species, 
                                     vcv = TRUE, 
                                     na.omit = FALSE, 
                                     warn.dropped = TRUE) 

## Again for the core bops
rownames(data.poly) <- data.poly$species
## Create comparative dataset for caper 
comparative.data.core <- comparative.data(phy = core.bop.tree, 
                                          data = data.poly, 
                                          names.col = species, 
                                          vcv = TRUE, 
                                          na.omit = FALSE, 
                                          warn.dropped = TRUE) 

#################################################################################
#################################################################################

##########################################################
### Analysis -- Correlations between complexity scores ###
##########################################################

## First, drop NA and non-real values 

data.1 <- filter(data, !(species %in% c("Lycocorax_pyrrhopterus", "Manucodia_jobiensis")))

data.1$species

# Create comparative dataset for caper 
comparative.data.1 <- comparative.data(phy = bop.tree, 
                                       data = data.1, 
                                       names.col = species, 
                                       vcv = TRUE, 
                                       na.omit = FALSE, 
                                       warn.dropped = TRUE) 

### Behavioural diversity versus complexity ###

## Model using OLS (no phylogeny)
m.ols.i<-gls(logBehav.Richness ~ logComplexity.Fux, data=data.1, method="ML")
summary(m.ols.i)

# Model using PGLS (phylogeny corrected) 
m.pgls.i <- pgls(logBehav.Richness ~ logComplexity.Fux, data = comparative.data.1, 
                 lambda = "ML")

summary(m.pgls.i)
nobs(m.pgls.i)
hist(residuals(m.pgls.i)) # Normality of residuals 

## Behavioural richness and behavioural complexity are not correlated!

### Behavioural richness versus complexity ###

## Model using OLS (no phylogeny)
m.ols.ii<-gls(logBehav.Diversity ~ logComplexity.Fux, data=data.1, method="ML")
summary(m.ols.ii)

# Model using PGLS (phylogeny corrected) 
m.pgls.ii <- pgls(logBehav.Diversity ~ logComplexity.Fux, data = comparative.data.1, 
                  lambda = "ML")
summary(m.pgls.ii)
hist(residuals(m.pgls.ii)) # Normality of residuals 

## Behavioural richness and behavioural diversity are also not correlated!

## What about the two metrics from the same publication (Ligon et al., 2018)? 

### Behavioural diversity versus richness ###

## Model using OLS (no phylogeny)
m.ols.iii<-gls(logBehav.Diversity ~ logBehav.Richness, data=data, method="ML")
summary(m.ols.iii)

# Model using PGLS (phylogeny corrected) 
m.pgls.iii <- pgls(logBehav.Diversity ~ logBehav.Richness, data = comparative.data.1, 
                   lambda = "ML")
summary(m.pgls.iii) 
hist(residuals(m.pgls.iii)) # Normality of residuals 

#################################################################################
#################################################################################

###################################################################
### Analysis -- OLS and PGLS for behavioural complexity metrics ###
###################################################################
### ABSOLUTE BRAIN VOLUME (NOT body size corrected) ### 
#######################################################

### Behavioural richness ### 

## Model using OLS (no phylogeny)
m.ols.1<-gls(logECV ~ logBehav.Richness, data=data.poly, method="ML")
summary(m.ols.1)

# Model using PGLS (phylogeny corrected) 
m.pgls.1 <- pgls(logECV ~ logBehav.Richness, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.1)
nobs(m.pgls.1)
hist(residuals(m.pgls.1)) # Normality of residuals; clearly not normally distributed. 

### Behavioural diversity ### 

## Model using OLS (no phylogeny)
m.ols.2<-gls(logECV ~ logBehav.Diversity, data=data.poly, method="ML")
summary(m.ols.2)

# Model using PGLS (phylogeny corrected) 
m.pgls.2 <- pgls(logECV ~ logBehav.Diversity, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.2)
hist(residuals(m.pgls.2)) # Normality of residuals; clearly not normally distributed. 

### Acoustic richness ### 

## Model using OLS (no phylogeny)
m.ols.3<-gls(logECV ~ logAcoust.Richness, data=data.poly, method="ML")
summary(m.ols.3)

# Model using PGLS (phylogeny corrected) 
m.pgls.3 <- pgls(logECV ~ logAcoust.Richness, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.3)
hist(residuals(m.pgls.3)) # Normality of residuals; clearly not normally distributed. 

### Acoustic diversity ### 

## Model using OLS (no phylogeny)
m.ols.4<-gls(logECV ~ logAcoust.Diversity, data=data.poly, method="ML")
summary(m.ols.4)

# Model using PGLS (phylogeny corrected) 
m.pgls.4 <- pgls(logECV ~ logAcoust.Diversity, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.4)
hist(residuals(m.pgls.4)) # Normality of residuals; clearly not normally distributed. 

### Behavioural complexity ### 

## Model using OLS (no phylogeny)
m.ols.5<-gls(logECV ~ logComplexity.Fux, data=data.poly, method="ML")
summary(m.ols.5)

# Model using PGLS (phylogeny corrected) 
m.pgls.5 <- pgls(logECV ~ logComplexity.Fux, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.5)
hist(residuals(m.pgls.5)) # Normality of residuals. 

#################################################################################
#################################################################################

###################################################################
### Analysis -- OLS and PGLS for behavioural complexity metrics ###
###################################################################
### RESIDUAL BRAIN VOLUME (BODY SIZE corrected) ### 
###################################################

### First we need to calculate relative brain size ###

## run a log-log pgls of ECV to body size and extract the residuals. 

## Quick inspection of the scaling relationship 
plot(data$logECV, data$logMass) 

# Run the pgls model
bm.pgls <- pgls(logECV ~ logMass, data = comparative.data, lambda = "ML")

## Check the output 
summary(bm.pgls) 

## Get the residuals
residuals(bm.pgls) # I added these to the data frame manually... 

## Did we successfully remove the effects of body size? 
summary(pgls(ResidECV ~ logECV, data = comparative.data, lambda = "ML"))
plot(data$logECV, data$ResidECV)
## It doesn't look like we *completely* removed the allometric effect of body size. 
## However, we can see that the adjusted R-squared value is rather low (0.1975), 
## so we accounted for body size quite a bit. 

###############

### Behavioural richness ### 

## Model using OLS (no phylogeny)
m.ols.6<-gls(ResidECV ~ logComplexity.Fux, data=data.poly, method="ML")
summary(m.ols.6)

# Model using PGLS (phylogeny corrected) 
m.pgls.6 <- pgls(ResidECV ~ logBehav.Richness, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.6)
hist(residuals(m.pgls.6)) # Normality of residuals. 

### Behavioural diversity ### 

## Model using OLS (no phylogeny)
m.ols.7<-gls(ResidECV ~ logBehav.Diversity, data=data.poly, method="ML")
summary(m.ols.7)

# Model using PGLS (phylogeny corrected) 
m.pgls.7 <- pgls(ResidECV ~ logBehav.Diversity, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.7)
hist(residuals(m.pgls.7)) # Normality of residuals. 

### Acoustic richness ### 

## Model using OLS (no phylogeny)
m.ols.8<-gls(ResidECV ~ logAcoust.Richness, data=data.poly, method="ML")
summary(m.ols.8)

# Model using PGLS (phylogeny corrected) 
m.pgls.8 <- pgls(ResidECV ~ logAcoust.Richness, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.8)
hist(residuals(m.pgls.8)) # Normality of residuals. 

### Acoustic diversity ### 

## Model using OLS (no phylogeny)
m.ols.9<-gls(ResidECV ~ logAcoust.Diversity, data=data.poly, method="ML")
summary(m.ols.9)

# Model using PGLS (phylogeny corrected) 
m.pgls.9 <- pgls(ResidECV ~ logAcoust.Diversity, data = comparative.data.core, 
                 lambda = "ML")
summary(m.pgls.9)
hist(residuals(m.pgls.9)) # Normality of residuals. 

### Behavioural complexity ### 

## Model using OLS (no phylogeny)
m.ols.10<-gls(ResidECV ~ logComplexity.Fux, data=data.poly, method="ML")
summary(m.ols.10)

# Model using PGLS (phylogeny corrected) 
m.pgls.10 <- pgls(ResidECV ~ logComplexity.Fux, data = comparative.data.core, 
                  lambda = "ML")
summary(m.pgls.10)
hist(residuals(m.pgls.10)) # Normality of residuals. 

#################################################################################
#################################################################################

### Clade (represented by nating system) and brain volume ### 

## We set lambda = 1 in both models (all manucodes are monogamous and form a single clade). 

m.pgls.11 <- pgls(logECV ~ as.factor(mating.system), data = comparative.data, 
                  lambda = 1) # We set lambda to '1' here as we don't need to compute it using ML, since all manucodes are one clade.
summary(m.pgls.11)
hist(residuals(m.pgls.11)) # Normality of residuals; clearly not normally distributed. 


# Model using PGLS (phylogeny corrected) 
m.pgls.12 <- pgls(ResidECV ~ as.factor(mating.system), data = comparative.data, 
                  lambda = 1) # We set lambda to '1' here as we don't need to compute it using ML, since all manucodes are one clade.
summary(m.pgls.12)
hist(residuals(m.pgls.12)) # Normality of residuals; more or less. 

#######################################
### Ancestral trait reconstructions ###
#######################################

## Load required packages for plotting 
library(viridis)
library(ggtree)
library(scico)

## First, let's plot the tree 
plotTree(bop.tree,fsize=0.9,ftype="i",lwd=1)

## Convert to a vector 
xx<-setNames(data$ECV,rownames(data))
xx

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(bop.tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

# Fit an ancestral state character reconstruction
fit <- phytools::fastAnc(bop.tree, xx, vars = TRUE, CI = TRUE)

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(bop.tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(bop.tree, d, by = 'node')

# Adjust plot margins and xlim
ggtree(
  tree, 
  aes(color = trait), 
  layout = "rectangular", 
  ladderize = FALSE, continuous = "color", size = 1.5) +
  scale_color_viridis(option = "plasma") +  
  #  scale_color_scico(palette = "buda") + 
  theme(
    legend.position = c(0.2, .74),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text.x = element_text(size = 10)) + # adjust x-axis text size
  labs(color = "ECV mm³") + 
  geom_point2(aes(subset=(node==30)), shape=16, size=3.5, colour='black') + 
  geom_point2(aes(subset=(node==25)), shape=16, size=3.5, colour='black') + 
  geom_tiplab(aes(label = gsub("_", " ", label), fontface = "italic"), size = 3, hjust = -0.05) + # Replace underscores with spaces and set font to italic
  xlim(-10, 50) 

### Residual ECV ###

## Convert to a vector 
yy<-setNames(data$ResidECV,rownames(data))
yy

## Change data frame to a vector
yy<-as.matrix(yy)[,1]
yy

## Estimate ancestral states 
fit <- fastAnc(bop.tree, yy, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(yy)

# Fit an ancestral state character reconstruction
fit <- phytools::fastAnc(bop.tree, yy, vars = TRUE, CI = TRUE)

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(bop.tree, names(yy)),
  trait = yy)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(bop.tree, d, by = 'node')

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(bop.tree, d, by = 'node')

# Adjust plot margins and xlim
ggtree(
  tree,
  aes(color = trait), 
  layout = "rectangular", 
  ladderize = FALSE, continuous = "color", size = 1.5) +
  scale_color_viridis(option = "viridis") +  
  #  scale_color_scico(palette = "buda") + 
  theme(
    legend.position = c(0.2, .74),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text.x = element_text(size = 10)) + # adjust x-axis text size
  labs(color = "Residual ECV") + 
  geom_point2(aes(subset=(node==30)), shape=16, size=4, colour='black') + 
  geom_point2(aes(subset=(node==25)), shape=16, size=4, colour='black')  + 
  geom_tiplab(aes(label = gsub("_", " ", label), fontface = "italic"), size = 3, hjust = -0.05) + # Replace underscores with spaces and set font to italic
  xlim(-10, 50) 

#################################################################################
#################################################################################

