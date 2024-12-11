# File which aggregates the dataset results into a single file and comparable plot

# Load the required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(xtable)

# Create a custom theme
custom_theme <- theme_minimal(base_size = 15) 

# Load and aggregate percentile plots

ausprivauto <- read.csv("results/ausprivate/percentile_data.csv")
freMTPLfreq <- read.csv("results/fretmtpl/percentile_data.csv")
swmotor <- read.csv("results/swmotor/percentile_data.csv")
brvehins1b <- read.csv("results/brveh/percentile_data.csv")

# Add labels
ausprivauto$dataset <- "ausPrivateAuto"
freMTPLfreq$dataset <- "freMTPLfreq"
swmotor$dataset <- "swmotor"
brvehins1b$dataset <- "brvehins1b"

# Aggregate the data
aggdata <- rbind(ausprivauto, freMTPLfreq, swmotor, brvehins1b)

aggdata %>% head(10)

# Create ggplot of the different datasets and add labels
p <- ggplot(aggdata, aes(x = percentile, y = W_avg, colour = dataset)) + 
    geom_line(aes(group = dataset), linewidth = 1, alpha = 0.7) + 
    labs(title = "Exposure means per percentile", x = "Percentile", y = "duration")+custom_theme
ggsave("results/plots/aggdurationplot.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)
p

# Create aggregated prediction to percentile with claims as scatter plot around
p <- ggplot(aggdata, aes(x = percentile, y = mu_sorted_avg, colour = dataset)) + 
    geom_line(aes(group = dataset), linewidth = 1, alpha = 0.7) + 
    geom_point(aes(x = percentile, y = ClaimNb_avg, colour = dataset), size = 0.5, alpha = 0.5) + 
    labs(title = "Aggregated Prediction to Percentile with Claims", x = "Percentile", y = "duration")+custom_theme
ggsave("results/plots/ModelFit_agg.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)

# Chec Dispersion Differences
p <- ggplot(aggdata, aes(x = percentile, y = rho_hat_p_avg, colour = dataset)) + 
    geom_point(aes(group = dataset), size = 0.5, alpha = 0.5) + 
    labs(title = "Dispersion Differences", x = "Percentile", y = "Deviance Estimates")+custom_theme+
    coord_cartesian(ylim = c(0, 4))
ggsave("results/plots/DispersionDifferences.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)


xtable(aggdata %>% group_by(dataset) %>% 
    summarise( mean(rho_hat_p_avg),mean(rho_bar_avg),mean(rho_tilde_avg)))