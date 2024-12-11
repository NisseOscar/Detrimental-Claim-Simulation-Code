## Code for creating plots of the simulated data

# Load the required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

################ Plot Results from Estimation of Parameters ################
df <- read.csv("./results/probability_estimates_results.csv") %>% select(-X)
df %>% head(10)

# Create a line plot with confidence interval for estimate p_hat and real value p
df_agg <- df %>% group_by(p) %>% summarize(mean_p_hat = mean(p_hat), lower = mean(p_hat) - 1.96 * sd(p_hat), upper = mean(p_hat) + 1.96 * sd(p_hat))
p <- ggplot(df_agg, aes(x = p, y = mean_p_hat)) +
    geom_ribbon(data = df_agg, aes(ymin = lower, ymax = upper), alpha = 0.2,fill="blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "Real p", y = "Estimated p") +
    coord_cartesian(ylim = c(0.01,0.99),xlim=c(0.05,0.95)) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom")
ggsave("./results/plots/estimated_vs_real_p.png", p, width = 18, height = 18, units = "cm")
p

# Plot histogram for one p group
p <-ggplot(df %>% filter(p==0.1), aes(x = p_hat)) +
    geom_histogram(bins=30, position = "dodge") +
    labs(title = "Histogram of estimated p values", x = "p", y = "Frequency") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom")
ggsave("./results/plots/histogram_p.png", p, width = 24, height = 18, units = "cm")

############## Plots Results from Simulation Estimates ################
df <- read.csv("./results/result_sim_estimates.csv") %>% select(-X)
df %>% head(10)

n_select = 100000
########## Plot duration
# Plot dispersion for the two groups
df_agg <- df %>%
    group_by(Duration_category,p) %>% 
    filter(Duration_category == "Without Detrimental Claims" | Duration_category== "With Detrimental Claims") %>%
    summarize(
        mean_dispersion = median(Dispersion), 
        lower = median(Dispersion)*qchisq(0.025, df = n_select-2)/(n_select-2),
        upper = median(Dispersion) * qchisq(0.975, df = n_select - 2) / (n_select - 2)
    )
df_agg

p <- ggplot(df_agg, aes(x = p, y = mean_dispersion, color = Duration_category)) +
    geom_line( size=0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper,fill = Duration_category,color=Duration_category), alpha = 0.2,colour = NA) +
    labs(title = paste("Dispersion estimate for different p values, n=100000"), x = "p", y = "Dispersion") +
    theme_minimal(base_size = 15)+ coord_cartesian(ylim = c(0, 10))+theme(legend.position = "bottom")
ggsave("./final_results/plots/dispersion_with_detrimental.png", p, width = 24, height = 18, units = "cm")
p

#### Coefficents for detrimental claims
df_agg <- df %>%
    filter(Duration_category == "With Detrimental Claims") %>%
    pivot_longer(
        cols = starts_with("B"),
        names_to = c("B", ".value"),
        names_pattern = "B(\\d+)_(.*)"
    ) %>%
    group_by(p, B) %>%
    summarize(est = mean(est), std = mean(std), true = mean(true))
df_agg$lower <- exp(df_agg$est - 1.96 * df_agg$std)
df_agg$upper <- exp(df_agg$est + 1.96 * df_agg$std)
df_agg$est <- exp(df_agg$est)
df_agg$true <- exp(df_agg$true)

p <- ggplot(df_agg, aes(x = p, y = est, color = B, fill = B)) +
    geom_line() +
    geom_line(aes(y = true), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = B, color = B), alpha = 0.2, colour = NA) +
    labs(title = "Estimated intensities for data with detrimental claims, n=100 000", x = "p", y = "G Coef") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 0.12))
ggsave("./final_results/plots/estimates_with_detrimental.png", p, width = 24, height = 18, units = "cm")
p

######## Results for adjusted duration

## Dispersion
df_agg <- df %>%
    group_by(Duration_category, p) %>%
    filter(Duration_category == "Detrimental Claim adjusted Duration" | Duration_category == "Without Detrimental Claims") %>%
    summarize(
        mean_dispersion = median(Dispersion),
        lower = median(Dispersion) * qchisq(0.025, df = n_select - 2) / (n_select - 2),
        upper = median(Dispersion) * qchisq(0.975, df = n_select - 2) / (n_select - 2)
    )
df_agg

p <- ggplot(df_agg, aes(x = p, y = mean_dispersion, color = Duration_category)) +
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Duration_category, color = Duration_category), alpha = 0.2, colour = NA) +
    labs(title = "Dispersion estimate after duration adjustment for different p values, n=100000", x = "p", y = "Dispersion") +
    theme_minimal(base_size = 15) +
    coord_cartesian(ylim = c(0.5, 1.5)) +
    theme(legend.position = "bottom")
ggsave("./final_results/plots/dispersion_adj_duration.png", p, width = 24, height = 18, units = "cm")
p


df %>% distinct(Duration_category)
#### Coefficents for detrimental claims
df_agg <- df %>%
    filter(Duration_category == "Detrimental Claim adjusted Duration") %>%
    pivot_longer(
        cols = starts_with("B"),
        names_to = c("B", ".value"),
        names_pattern = "B(\\d+)_(.*)"
    ) %>%
    group_by(p, B) %>%
    summarize(est = mean(est), std = mean(std), true = mean(true))
df_agg$lower <- exp(df_agg$est - 1.96 * df_agg$std)
df_agg$upper <- exp(df_agg$est + 1.96 * df_agg$std)
df_agg$est <- exp(df_agg$est)
df_agg$true <- exp(df_agg$true)

p <- ggplot(df_agg, aes(x = p, y = est, color = B, fill = B)) +
    geom_line() +
    geom_line(aes(y = true), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = B, color = B), alpha = 0.2, colour = NA) +
    labs(title = "Estimated intensities after adjusting duration to detrimental claims, n=100 000", x = "p", y = "G Coef") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 0.12))
ggsave("./final_results/plots/estimates_adj_duration.png", p, width = 24, height = 18, units = "cm")
p
