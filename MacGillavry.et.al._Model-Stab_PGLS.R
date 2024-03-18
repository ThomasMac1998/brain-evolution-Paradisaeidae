## Packages
library(ggplot2) # Plotting
library(dplyr) # Data wrangling 
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis

#####################################################
### Create a new dataset for each missing species ###
#####################################################

# Create a list to store new datasets
dropped_datasets <- list()

# Loop through each observation and create a new dataset without that observation
for (i in 1:nrow(data.poly)) {
  dropped_datasets[[i]] <- data.poly[-i, ]
}

# Accessing any of the dropped datasets, for example, the first one:
dropped_datasets[[1]]

#######################################################
### Now turn all of these into comparative datasets ###
#######################################################

# Create a list to store comparative dataframes
comparative_dataframes <- list()

# Load the phylogenetic tree (assuming core.bop.tree is available)
core_bop_tree <- core.bop.tree

# Loop through each dropped dataset and create a comparative dataframe
for (i in 1:length(dropped_datasets)) {
  # Create comparative dataframe for the current dropped dataset
  comp_data <- comparative.data(phy = core_bop_tree,
                                data = dropped_datasets[[i]],
                                names.col = "species",
                                vcv = TRUE,
                                na.omit = FALSE,
                                warn.dropped = TRUE)
  
  # Add the comparative dataframe to the list
  comparative_dataframes[[i]] <- comp_data
}

# Accessing any of the comparative dataframes, for example, the first one:
comparative_dataframes[[1]]

######################################################
### Now run a model for each one of these datasets ###
######################################################

# Create an empty list to store the models and their summaries
model_summary_list <- list()

# Loop through each dropped dataset
for (i in 1:length(comparative_dataframes)) {
  # Fit the PGLS model on the current dataset
  m_pgls <- pgls(ResidECV ~ logComplexity.Fux, data = comparative_dataframes[[i]], lambda = "ML")
  
  # Obtain the summary of the model
  m_summary <- summary(m_pgls)
  
  # Store the model and its summary in the list
  model_summary_list[[i]] <- list(model = m_pgls, summary = m_summary)
}

# Now you have a list of PGLS models and their summaries stored in model_summary_list
model_summary_list[1]

##################################################################
### Extracting the estimate and p-val for behavioural richness ###
##################################################################

# Create an empty dataframe to store estimates and p-values
estimate_pvalue_dataframe <- data.frame(source_model = character(),
                                        estimate = numeric(),
                                        p_value = numeric(),
                                        stringsAsFactors = FALSE)

# Loop through each model summary
for (i in 1:length(model_summary_list)) {
  # Extract estimate for logBehav.Richness from the model summary
  estimate <- coef(model_summary_list[[i]]$model)["logComplexity.Fux"]
  
  # Extract p-value for logBehav.Richness from the model summary
  p_value <- model_summary_list[[i]]$summary$coefficients["logComplexity.Fux", "Pr(>|t|)"]
  
  # Add the values to the dataframe along with the source model
  estimate_pvalue_dataframe <- rbind(estimate_pvalue_dataframe, 
                                     data.frame(source_model = paste("Model", i),
                                                estimate = estimate,
                                                p_value = p_value,
                                                stringsAsFactors = FALSE))
}

################
### Plotting ###
################

library(ggplot2)

ggplot(estimate_pvalue_dataframe, aes(x = "", y = estimate, color = p_value)) +
  geom_jitter(size = 3, width = 0.15, shape = ifelse(estimate_pvalue_dataframe$p_value > 0.05, 1, 19)) + 
  geom_boxplot(alpha = 0, size = 0.5, width = 0.07) +
  geom_violin(alpha = 0, size = 0.5) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "",
       y = "Estimate",
       color = "p-value", 
       title = expression("Stability of PGLS model: Residual ECV ~ log"[10]*" Behavioural Complexity")) + # Adding a subscript '10' using expression
  theme_minimal() + 
  ylim(-0.05, 0.35) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_text(size = 20)) +
  guides(color = guide_colorbar(title = "p-val")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
  coord_flip() 

## Which species can be dropped to remove this effect? 

print(filter(estimate_pvalue_dataframe, p_value >= 0.05))

comparative_dataframes[[3]]$data # Ciccinurus magnificus
comparative_dataframes[[5]]$data # Ciccinurus respublica
comparative_dataframes[[8]]$data # Paradisaea minor
comparative_dataframes[[9]]$data # Paradisaea rubra
comparative_dataframes[[11]]$data # Parotia lawesii


