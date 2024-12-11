#### Replication notebook of paper ON DURATION EFFECTS IN NON-LIFE INSURANCE PRICING
# A notebook to replicate the results of the paper to be used for initial work of the project
# https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4474908
##########################################################################################

# packages
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)
library(gbm)
library(CASdatasets)

# Load the data
data(ausprivauto0405)
df_aus <- ausprivauto0405
df_aus$duration <- df_aus$Exposure
df_aus$z <- df_aus$ClaimNb
ausprivate_formula <- z ~ offset(log(duration)) + VehValue + C(VehAge) + C(VehBody) + C(Gender) + C(DrivAge)

# Load the data
data(brvehins1b)
brveh_df <- brvehins1b %>% filter(ExposTotal <=1) %>% filter(ExposTotal >0)
brveh_df$duration <- brveh_df$ExposTotal
brveh_df$z <- brveh_df$ClaimNbTotColl
brveh_formula <- z ~ offset(log(duration)) + C(Gender) + C(DrivAge) + C(VehYear) + C(Area)

# Load the data
fretmtpl_df <- data.frame(read.csv("Datasets/freMTPLfreq.csv"))
fretmtpl_df$duration <- fretmtpl_df$Exposure
fretmtpl_df$z <- fretmtpl_df$ClaimNb
fretmtpl_formula <- z ~ offset(log(duration)) + C(Power) + CarAge + DriverAge + C(Brand) + C(Gas) + C(Region) + Density

# Load the data
data(swmotorcycle)
swmotor_df <- swmotorcycle %>% filter(Exposure > 0) %>% filter(Exposure <=1)
swmotor_df$duration <- swmotor_df$Exposure
swmotor_df$z <- swmotor_df$ClaimNb
swmotor_formula <- z ~ offset(log(duration)) + OwnerAge + C(Gender) + C(Area) + C(RiskClass) + VehAge + C(BonusClass)


datasets <- c("ausprivate", "brveh", "fretmtpl", "swmotor")
formulas <- c(ausprivate_formula, brveh_formula, fretmtpl_formula, swmotor_formula)
df_list <- list(df_aus, brveh_df, fretmtpl_df, swmotor_df)

# Define helper functions
detrimental_adjustment <- function(W, Z) {
    # Calculate Deteriminetal adjusted duration
    p_n_hat <- mean(W[Z == 0] < 1)
    p_z_hat <- mean(W[Z > 0] < 1)
    p_hat <- (p_z_hat - p_n_hat / 2) / (1 - p_n_hat / 2)

    # Calculate adjusted duration
    q_hat <- sapply(W, function(x) ifelse(x > 0, (1 - p_hat), 0))
    q_hat[W >= 1] <- 0
    p_c_hat <- sapply(W, function(x) ifelse(x < 1, p_n_hat * x, 1))
    prob_natural <- q_hat * p_c_hat

    # # Expected Full Duration
    p_n_remain <- (1 - W) * p_n_hat
    expected_remaining_duration <- (1 - p_n_remain) + p_n_remain * (W + 1) / 2

    # Calculate duration adjusted for detrimental claims
    W_adj <- (prob_natural * W + (1 - prob_natural) * expected_remaining_duration)

    # Adjust durations for to follow expected values
    W_adj <- W_adj * sum(W) / sum(W_adj)

    return(W_adj)
}


for (i in 1:length(df_list)){
    set.seed(1234)
    df <- df_list[[i]]
    formula <- formulas[[i]]
    dataset <- datasets[i]
    n <- nrow(df)
 
    # Create gbm Model
    model <- gbm(formula, data = df,train.fraction = 0.8, distribution = "poisson", n.trees = 300, interaction.depth = 2, shrinkage = 0.1, verbose = TRUE)
    
    # Save optimal Model
    png(paste("./results/",dataset,"/plots/variable_importance.png",sep=""), width = 800, height = 600)
    summary(model)
    dev.off()

    # Save optimal number of trees
    png(paste("./results/",dataset,"/plots/optimal_para.png",sep=""), width = 800, height = 600)
    best_iter <- gbm.perf(model, method = "test")
    dev.off()
    best_iter
    df %>% head(19)
    # Study distribution of predictions
    mu <- predict(model, df, n.trees = best_iter, type = "response")
    test <- glm(df$z ~ log(df$duration), offset = log(mu), family = poisson)
    summary(test)

    ########## Plot Replication #####################################
    # Create parameters for moving average and percentile
    k <- floor(n / 200) + 1
    percentile <- 1:(floor(n / k) + 1) / (floor(n / k) + 1)

    # Create index of smallest to largest mu
    mu_sorted_idx <- order(mu, decreasing = FALSE)
    mu_sorted <- mu[mu_sorted_idx]
    mu_sorted_avg <- sapply(percentile, function(i) mean(mu_sorted[floor(i * n - k + 1):floor(i * n)]))
    
    # Create moving average  of exposure based on mu_sorted_idx
    W_sorted <- df$duration[mu_sorted_idx]
    W_mvavg <- sapply(k:(n - k), function(i) mean(W_sorted[(i - k + 1):(i + k)]))
    W_mvavg <- c(rep(W_mvavg[1], k), W_mvavg, rep(W_mvavg[length(W_mvavg)], k - 1))
    W_avg <- sapply(percentile, function(i) mean(W_sorted[floor(i * n - k + 1):floor(i * n)]))
    
    # Create Average Claim Frequency
    ClaimNb_sorted <- df$z[mu_sorted_idx]
    ClaimNb_avg <- sapply(percentile, function(i) mean(df$z[floor(i * n - k + 1):floor(i * n)]))

    # Calculate W adjusted
    W_detadj_sorted <- detrimental_adjustment(W_sorted,ClaimNb_sorted)
    
    # Plot W_avg and compare percentile
    p <- ggplot() +
        geom_line(aes(x = percentile, y = W_avg, colour = "Moving Average"), color = "#000000", linewidth = 1, alpha = 0.7) +
        labs(x = "Percentile", y = "duration") +
        ylim(0, max(W_avg)*1.1) #+ legend("topright", c("Moving Average", "Observations"), fill = c("#000000", "#686868"))
    ggsave(paste("./results/", dataset, "/plots/duration_means.png", sep = ""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

    # Plot Prediction and claim frequency towards percentile
    p <- ggplot() +
        geom_point(aes(x = percentile, y = mu_sorted_avg), color = "#000000", alpha = 0.7) +
        geom_point(aes(x = percentile, y = ClaimNb_avg), color = "#FF0000", size = 1, alpha = 0.7) +
        labs(x = "Percentile", y = "Average Claim Frequency")
    ggsave(paste("./results/",dataset,"/plots/model_fit.png",sep=""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

    # Calculate Dispersion Parameter
    sig_hat_p <- (ClaimNb_sorted - W_sorted * mu_sorted)^2
    sig_hat_p_avg <- sapply(percentile, function(i) mean(sig_hat_p[floor(i * n - k + 1):floor(i * n)]))
    rho_hat_p <- sig_hat_p / (W_sorted * mu_sorted)
    rho_hat_p_mvavg <- sapply(k:(n - k), function(i) mean(rho_hat_p[(i - k + 1):(i + k)]))
    rho_hat_p_avg <- sapply(percentile, function(i) mean(rho_hat_p[(i * n - k + 1):(i * n)]))
    rho_hat_p_mean <- mean(rho_hat_p_avg)
    rho_hat_p_mean

    # Dismissing the duration effect
    sig_bar <- (ClaimNb_sorted - W_mvavg * mu_sorted)^2
    sig_bar_mvavg <- sapply(percentile, function(i) mean(sig_bar[(i * n - k + 1):(i * n)]))
    sig_bar_avg <- sapply(percentile, function(i) mean(sig_bar[floor(i * n - k + 1):floor(i * n)]))
    rho_bar <- sig_bar / (W_mvavg * mu_sorted)
    rho_bar_avg <- sapply(percentile, function(i) mean(rho_bar[(i * n - k + 1):(i * n)]))
    rho_bar_mean <- mean(rho_bar_avg)
    rho_bar_mean

    # Adjusted duration effect
    sig_tilde <- (ClaimNb_sorted - W_detadj_sorted * mu_sorted)^2
    sig_tilde_avg <- sapply(percentile, function(i) mean(sig_tilde[(i * n - k + 1):(i * n)]))
    rho_tilde <- sig_tilde / (W_detadj_sorted * mu_sorted)
    rho_tilde_avg <- sapply(percentile, function(i) mean(rho_tilde[(i * n - k + 1):(i * n)]))
    rho_tilde_mean <- mean(rho_tilde_avg)
    rho_tilde_mean

    # Plot anualized amount of claims with standard deviations
    p <- ggplot() +
        geom_line(aes(x = percentile, y = mu_sorted_avg), color = "#000000", alpha = 0.7) +
        geom_point(aes(x = percentile, y = ClaimNb_avg), color = "#696969", size = 1, alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg + sqrt(sig_hat_p_avg / k)), color = "#FF7575", alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg - sqrt(sig_hat_p_avg / k)), color = "#FF7575", alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg - sqrt(sig_bar_avg / k)), color = "#ffee05", alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg + sqrt(sig_bar_avg / k)), color = "#ffee05", alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg - sqrt(sig_tilde_avg / k)), color = "#2966ff", alpha = 0.7) +
        geom_line(aes(x = percentile, y = mu_sorted_avg + sqrt(sig_tilde_avg / k)), color = "#2966ff", alpha = 0.7) +
        labs(x = "Percentile", y = "Average Claim Frequency")+
        theme_minimal(base_size = 15)
    ggsave(paste("./results/",dataset,"/plots/Anualized_Claim_Amount.png",sep=""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

    # Plot Deviance Estimates with mean
    p <- ggplot() +
        geom_point(aes(x = percentile, y = rho_bar_avg), color = "#FF7575", alpha = 0.7) +
        geom_point(aes(x = percentile, y = rho_hat_p_avg), color = "#000000", alpha = 0.7) +
        geom_point(aes(x = percentile, y = rho_tilde_avg), color = "#2966ff", alpha = 0.7) +
        geom_hline(yintercept = rho_bar_mean, color = "#FF7575", linetype = "dashed") +
        geom_hline(yintercept = rho_hat_p_mean, color = "#000000", linetype = "dashed") +
        geom_hline(yintercept = rho_tilde_mean, color = "#2966ff", linetype = "dashed") +
        labs(x = "Percentile", y = "Deviance Estimates")+ 
        theme_minimal(base_size = 15) + coord_cartesian(ylim = c(0, 1.5 * rho_hat_p_mean))
    ggsave(paste("./results/",dataset,"/plots/Deviance_Estimates.png",sep=""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

    # Create dataframe and save to csv file in the DataSets folder
    q_df <- data.frame(percentile, mu_sorted_avg, ClaimNb_avg, sig_hat_p_avg, sig_bar_avg,sig_tilde_avg, rho_hat_p_avg, rho_bar_avg,rho_tilde_avg, W_avg)
    write.csv(q_df, paste("./results/",dataset,"/percentile_data.csv",sep=""), row.names = FALSE)

}
