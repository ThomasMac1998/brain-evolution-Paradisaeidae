
## Load required packages 
library(grid)
library(cowplot)

## Note: first run the code in the analysis script to import relevant dataframes and phylogenetic trees. 

#################################################################################
#################################################################################

### First, we want to see how the brains of the birds of paradise compare to other taxa: 

## Import the full dataset by Ksepka et al.: 
all.birds <- read.csv("MacGillavry.et.al._Brain.Data.AllBirds.csv")
head(all.birds, 52) # Quick inspection 

## Figure 1b: 
## Plot showing relationship of body mass to ECV across all songbirds 
## in our extended µCT-based ECV dataset (see methods)

ggplot(data = filter(all.birds, Order == "Passeriformes"), 
       aes(x = log10(mass/1000), y = log10(ECV), color = Family == "Paradisaeidae")) +
  geom_point(size = 3, stroke = 0.5, shape = 1, colour = "darkgrey") + 
  geom_point(size = 3, shape = 16) + 
  scale_color_manual(values = c("grey90", "#49A4B9")) + 
  ylab(expression(paste("log"[10], " ECV (mm"^"3", ")"))) + 
  xlab(expression(paste("log"[10], " MASS (g)"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=17),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(color = FALSE) 

## Figure 1c: 
## Plot showing frequency histogram of absolute ECV across songbirds. 

p1 <- ggplot(data = filter(all.birds, Order == "Passeriformes"), 
             aes(x = (ECV/1000), fill = Family == "Paradisaeidae")) +
  geom_histogram(alpha = 1, bins = 35, colour = "black") + 
  scale_fill_manual(values = c("grey90", "#49A4B9")) + 
  xlab(expression(paste("ECV (cm"^"3", ")"))) + 
  ylab("Count") +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +  # Adjust y-axis ticks
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),  # Adjust size of axis tick labels
        axis.title = element_text(size = 15)) +  # Adjust size of axis titles
  guides(fill = FALSE)

## Figure 1d: 
## Plot showing frequency histogram of residual ECV across songbirds. 

p2 <- ggplot(data = filter(all.birds, Order == "Passeriformes"), 
             aes(x = ResidECV, fill = Family == "Paradisaeidae")) +
  geom_histogram(alpha = 1, bins = 35, colour = "black") + 
  scale_fill_manual(values = c("grey90", "#49A4B9")) + 
  xlab(expression(paste("Residual ECV"))) + 
  ylab("Count") +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),  # Adjust size of x-axis tick labels
        axis.title = element_text(size = 15)) +  # Adjust size of y-axis tick labels
  guides(fill = FALSE) + 
  ylab("Count")

# Arrange plots vertically
ggarrange(p1, p2,
          ncol = 1, nrow = 2, 
          align = "v")

#################################################################################
#################################################################################

## Figures 2a-f and 3a-d: 
## Bivariate plots showing correlations between brain size metrics and motor complexity scores 
## as well as acoustic complexity scores. 

## Behavioural richness (Ligon et al.)
br_1 <- ggplot(data = data.poly, 
               aes(x = log10(behavioral.richness.ligon), y = log10(ECV))) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#1B9E77") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") +  
  ylab(expression(paste("log"[10], " ECV (mm"^"3", ")"))) + 
  xlab(" ") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Behavioural diversity (Ligon et al.)
bd_1 <- ggplot(data = data.poly, 
               aes(x = log10(behavioral.diversity.ligon), y = log10(ECV))) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#D95F02") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") +  
  ylab(NULL) + 
  xlab(" ") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Behavioural complexity (Miles & Fuxjager)
bc_1 <- ggplot(data = data.poly,  
               aes(x = log10(complexity.fuxjager), y = log10(ECV))) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#7570B3") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") + 
  ylab(NULL) + 
  xlab(" ") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Acoustic richness (Ligon et al.)
ar_1 <- ggplot(data = data.poly, 
               aes(x = log10(acoustic.richness.ligon), y = log10(ECV))) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#E7298A") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") +  
  ylab(expression(paste("log"[10], " ECV (mm"^"3", ")"))) + 
  xlab(" ") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Acoustic diversity (Ligon et al.)
ad_1 <- ggplot(data = data.poly, 
               aes(x = log10(acoustic.diversity.ligon), y = log10(ECV))) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#66A61E") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") + 
  ylab(NULL) + 
  xlab(" ") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Behavioural richness (Ligon et al.)
br_2 <- ggplot(data = data.poly, 
               aes(x = log10(behavioral.richness.ligon), y = ResidECV)) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#1B9E77") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") +  
  ylab(expression(paste("Residual ECV"))) + 
  xlab(expression(paste("log"[10], " BEHAVIOURAL RICHNESS"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) +
  #  annotate("text", x = Inf, y = Inf, label = "n.s.", hjust = 1.1, vjust = 1.1, size = 6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Behavioural diversity (Ligon et al.)
bd_2 <- ggplot(data = data.poly, 
               aes(x = log10(behavioral.diversity.ligon), y = ResidECV)) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#D95F02") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") + 
  ylab(NULL) + 
  xlab(expression(paste("log"[10], " BEHAVIOURAL DIVERSITY"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Behavioural complexity (Miles & Fuxjager)
bc_2 <- ggplot(data = data.poly,  
               aes(x = log10(complexity.fuxjager), y = ResidECV)) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#7570B3") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") + 
  ylab("NULL") + 
  xlab(expression(paste("log"[10]," BEHAVIOURAL COMPLEXITY"))) +
  geom_abline(intercept = -0.27892, slope = 0.25481, color="black", size = 0.7, 
              linetype="dashed") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Acoustic richness (Ligon et al.)
ar_2 <- ggplot(data = data.poly, 
               aes(x = log10(acoustic.richness.ligon), y = ResidECV)) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#E7298A") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") +  
  ylab(expression(paste("Residual ECV"))) + 
  xlab(expression(paste("log"[10], " ACOUSTIC RICHNESS"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Acoustic diversity (Ligon et al.)
ad_2 <- ggplot(data = data.poly, 
               aes(x = log10(acoustic.diversity.ligon), y = ResidECV)) + 
  geom_point(size = 5, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 5, shape = 16, colour = "#66A61E") +
  geom_text(aes(label = species_number), vjust = 0.5, hjust = 0.5, size = 3, colour = "white") + 
  ylab("") + 
  xlab(expression(paste("log"[10], " ACOUSTIC DIVERSITY"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Combining all bivariate brain-complexity plots
ggarrange(br_1, bd_1, bc_1, 
          br_2, bd_2, bc_2, 
          ncol = 3, nrow = 2, 
          align = "v")

ggarrange(ar_1, ad_1, 
          ar_2, ad_2,
          ncol = 2, nrow = 2, 
          align = "v")

#################################################################################
#################################################################################

## Figures 4c, d: 
## Violin + boxplots showing ECV and residual ECV for the two major bop clades (represented by mating system). 

## Clade; absolute ECV
clade.plot.b <- ggplot(data = data, 
                        aes(x = mating.system, y = ResidECV, fill = mating.system)) + # Convert mass to g and ECV to cm^3
  geom_violin() +  
  geom_boxplot(width = 0.07) + 
  ylab(expression(paste("Residual ECV"))) + 
  xlab("CLADE") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=17),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#BCE4D8", "#49A4B9")) + 
  theme(legend.position = "NULL") +
  annotate("text", x = Inf, y = Inf, label = "n.s.", hjust = 1.1, vjust = 1.1, size = 6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Clade; residual ECV
clade.plot.a <- ggplot(data = data, 
                        aes(x = mating.system, y = log10(ECV), fill = mating.system)) +
  geom_violin() +  
  geom_boxplot(width = 0.07) + 
  ylab(expression(paste("log"[10], " ECV (mm"^"3", ")"))) + 
  xlab(NULL) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=17),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#BCE4D8", "#49A4B9")) + 
  theme(legend.position = "NULL") +
  annotate("text", x = Inf, y = Inf, label = "n.s.", hjust = 1.1, vjust = 1.1, size = 6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Combining mating system plots 
ggarrange(clade.plot.a, clade.plot.b,
          ncol = 1, nrow = 2)

#################################################################################
#################################################################################

## Figures S1: 
## Plots showing ECV values for all samples per species. 

sample.size <- read.csv("MacGillavry.et.al._Samples.csv")

## Calculate the count of observations for each group
observation_count <- sample.size %>%
  group_by(SPECIES) %>%
  summarize(count = n())

## Sort the data by the count of observations in increasing order
sorted_data <- observation_count %>%
  arrange(count)

## Plot the boxplot with sorted data and rotated x-axis labels
ggplot(data = sample.size, aes(x = factor(SPECIES, levels = sorted_data$SPECIES), y = ECV, colour = mating.system)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(xintercept = 1:length(levels(factor(sample.size$SPECIES))), color = "grey", linetype = "dashed") +  # Add vertical lines
  geom_point(size = 3, stroke = 1, shape = 1, colour = "black") + 
  geom_point(size = 3, shape = 16) +
  scale_colour_manual(name = "Clade", values = c("#BCE4D8", "#49A4B9")) +
  xlab(NULL) + 
  ylab(expression(paste("ECV (mm"^"3", ")"))) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) 

## Figures S2a-c: 
## Plots showing correlations between different behavioural complexity scores. 

## Richness vs complexity 
compl.plot.a <- ggplot(data = data.1, 
                       aes(x = logComplexity.Fux, y = logBehav.Richness, colour = mating.system)) + # Convert mass to g and ECV to cm^3
  geom_point(size = 4, shape = 1, colour = "black") + 
  geom_point(size = 4, shape = 16) +
  scale_colour_manual(values = c("#BCE4D8", "#49A4B9")) + 
  ylab(expression(paste("log"[10], " BEHAVIOURAL RICHNESS"))) + 
  xlab(expression(paste("log"[10], " BEHAVIOURAL COMPLEXITY"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank()) +
  annotate("text", x = Inf, y = Inf, label = "n.s.", hjust = 1.1, vjust = 1.1, size = 6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Diversity vs complexity 
compl.plot.b <- ggplot(data = data.1, 
                       aes(x = logComplexity.Fux, y = logBehav.Diversity, colour = mating.system)) + # Convert mass to g and ECV to cm^3
  geom_point(size = 4, shape = 1, colour = "black") + 
  geom_point(size = 4, shape = 16) +
  scale_colour_manual(values = c("#BCE4D8", "#49A4B9")) + 
  ylab(expression(paste("log"[10], " BEHAVIOURAL DIVERSITY"))) + 
  xlab(expression(paste("log"[10], " BEHAVIOURAL COMPLEXITY"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank()) +
  annotate("text", x = Inf, y = Inf, label = "n.s.", hjust = 1.1, vjust = 1.1, size = 6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Diversity versus richness
compl.plot.c <- ggplot(data = data, 
                       aes(x = logBehav.Richness, y = logBehav.Diversity, colour = mating.system)) + # Convert mass to g and ECV to cm^3
  geom_abline(intercept = -0.437271, slope = 1.236718, color="black", size = 0.7, 
              linetype="solid") + 
  geom_point(size = 4, shape = 1, colour = "black") + 
  geom_point(size = 4, shape = 16) + 
  scale_colour_manual(values = c("#BCE4D8", "#49A4B9")) + 
  ylab(expression(paste("log"[10], " BEHAVIOURAL DIVERSITY"))) + 
  xlab(expression(paste("log"[10], " BEHAVIOURAL RICHNESS"))) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

## Combining complexity score comparisons
ggarrange(compl.plot.a, compl.plot.b, compl.plot.c,
          ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "bottom", 
          labels = c("A", "B", "C"))


