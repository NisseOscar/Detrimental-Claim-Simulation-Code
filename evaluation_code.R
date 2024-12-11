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
swmotor_df <- swmotorcycle %>% filter(Exposure > 0) %>% filter(Exposure<=1)
swmotor_df$duration <- swmotor_df$Exposure
swmotor_df$z <- swmotor_df$ClaimNb
swmotor_formula <- z ~ offset(log(duration)) + OwnerAge + C(Gender) + C(Area) + C(RiskClass) + VehAge + C(BonusClass)


datasets <- c("fretmtpl","ausprivate","swmotor", "brveh")
formulas <- c(fretmtpl_formula,ausprivate_formula,swmotor_formula, brveh_formula)
df_list <- list(fretmtpl_df,ausprivate_formula,swmotor_df, brveh_df)


# Define helper functions
detrimental_adjustment <- function(W,Z){
    # Calculate Deteriminetal adjusted duration
    p_n_hat <- mean(W[Z == 0] < 1)
    p_z_hat <- mean(W[Z > 0] < 1)
    p_hat <- (p_z_hat - p_n_hat / 2) / (1 - p_n_hat / 2)

    # Calculate adjusted duration
    q_hat <- sapply(W, function(x) ifelse(x > 0, (1 - p_hat), 0))
    q_hat[W >= 1] <- 0
    p_c_hat = sapply(W, function(x) ifelse(x < 1, p_n_hat * x, 1))
    prob_natural <- q_hat * p_c_hat

    # # Expected Full Duration
    p_n_remain = (1 - W) * p_n_hat
    expected_remaining_duration <- (1 - p_n_remain) + p_n_remain * (W + 1) / 2

    # Calculate duration adjusted for detrimental claims
    W_adj <- (prob_natural * W + (1 - prob_natural) * expected_remaining_duration)

    # Adjust durations for to follow expected values
    W_adj <- W_adj * sum(W) / sum(W_adj)

    return(W_adj)
}

lindskogNazar <- function(mu,W){
    # Take out estimates
    n <- length(mu)
    k <- floor(n / 200) + 1

    # Lindskog & Nazar Adjustment
    mu_sorted_idx <- order(mu, decreasing = FALSE)
    mu_sorted <- mu[mu_sorted_idx]
    W_sorted <- W[mu_sorted_idx]
    W_mvavg <- sapply(k:(n - k), function(i) mean(W_sorted[(i - k + 1):(i + k)]))
    W_ln <- c(rep(W_mvavg[1], k), W_mvavg, rep(W_mvavg[length(W_mvavg)], k - 1))

    return(W_ln)
}

for (i in 4:length(df_list)){
    set.seed(1234)
    df <- df_list[[i]]
    formula <- formulas[[i]]
    dataset <- datasets[i]

    folds <- cut(seq(1, nrow(df)), breaks = 10, labels = FALSE)
    cv_errors <- rep(0, 10)
    tot_res <- data.frame()
    for (j in 1:10) {
        testIndexes <- which(folds == j, arr.ind = TRUE)
        train <- df[-testIndexes,]
        test <- df[testIndexes,]

        # Fit base model
        W_obs <- train$duration
        model <- gbm(formula, data = train, train.fraction = 0.8, distribution = "poisson", n.trees = 200, interaction.depth = 2, shrinkage = 0.1)
        best_iter <- gbm.perf(model, method = "test", plot.it = FALSE)
        mu_train <- predict(model, train, n.trees = best_iter, type = "response")
        mu <- predict(model, test, n.trees = best_iter, type = "response")

        # Lindskog & Nazar Adjustment
        W_ln <- lindskogNazar(mu_train, W_obs)
        train$duration <- W_ln
        model_ln <- gbm(formula, data = train, train.fraction = 0.8, distribution = "poisson", n.trees = 200, interaction.depth = 2, shrinkage = 0.1)
        best_iter_ln <- gbm.perf(model_ln, method = "test", plot.it = FALSE)
        mu_ln <- predict(model_ln, test, n.trees = best_iter_ln, type = "response")

        # Calculate Deteriminetal adjusted duration
        W_adj <- detrimental_adjustment(train$duration, train$z)
        train$duration <- W_adj
        model_adj <- gbm(formula, data = train, train.fraction = 0.8, distribution = "poisson", n.trees = 200, interaction.depth = 2, shrinkage = 0.1)
        best_iter_adj <- gbm.perf(model_adj, method = "test", plot.it = FALSE)
        mu_adj <- predict(model_adj, test, n.trees = best_iter_adj, type = "response")

        # Join into a single dataframe
        res <- data.frame(mu=mu, mu_ln=mu_ln, mu_adj=mu_adj, z =test$z, duration=test$duration)
        tot_res <- rbind(tot_res, res)
    }

    ########## Plot Replication #####################################

    models <- c('base','LindskogNazar','Detrimental')
    mus <- c('mu','mu_ln','mu_adj')

    for (j in 1:3){
        mu <- tot_res[,mus[j]]
        W_obs <- tot_res$duration
        model <- models[j]

        n <- nrow(tot_res)
        k <- floor(n / 200) + 1
        percentile <- 1:(floor(n / k) + 1) / (floor(n / k) + 1)

        # Create index of smallest to largest mu
        mu_sorted_idx <- order(mu, decreasing = FALSE)
        mu_sorted <- mu[mu_sorted_idx]
        mu_sorted_avg <- sapply(percentile, function(i) mean(mu_sorted[floor(i * n - k + 1):floor(i * n)]))

        # Create moving average  of exposure based on mu_sorted_idx
        W_sorted <- df$duration[mu_sorted_idx]
        W_avg <- sapply(percentile, function(i) mean(W_sorted[floor(i * n - k + 1):floor(i * n)]))
    
        # Create Average Claim Frequency
        ClaimNb_sorted <- tot_res$z[mu_sorted_idx]
        ClaimNb_avg <- sapply(percentile, function(i) mean(ClaimNb_sorted[floor(i * n - k + 1):floor(i * n)]))

        # Plot W_avg and compare percentile
        p <- ggplot() +
            geom_line(aes(x = percentile, y = W_avg, colour = "Moving Average"), color = "#000000", linewidth = 1, alpha = 0.7) +
            labs(title = "Moving Average of Exposure", x = "Percentile", y = "duration") +
            ylim(0, max(W_avg)*1.1) #+ legend("topright", c("Moving Average", "Observations"), fill = c("#000000", "#686868"))
        ggsave(paste("./results/", dataset,"/",model,"_plots/duration_means.png", sep = ""), plot = p, width = 10, height = 6, units = "in", dpi = 300)
        
        # Plot Prediction and claim frequency towards percentile
        p <- ggplot() +
            geom_point(aes(x = percentile, y = mu_sorted_avg), color = "#000000", alpha = 0.7) +
            geom_point(aes(x = percentile, y = ClaimNb_avg), color = "#FF0000", size = 1, alpha = 0.7) +
            labs(title = "Prediction and Claim Frequency towards Percentile", x = "Percentile", y = "Average Claim Frequency")
        ggsave(paste("./results/",dataset,"/",model,"_plots/model_fit.png",sep=""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

        # Calculate Dispersion Parameter
        sig_hat_p <- (ClaimNb_sorted - W_sorted * mu_sorted)^2
        sig_hat_p_avg <- sapply(percentile, function(i) mean(sig_hat_p[floor(i * n - k + 1):floor(i * n)]))
        rho_hat_p <- sig_hat_p / (W_sorted * mu_sorted)
        rho_hat_p_avg <- sapply(percentile, function(i) mean(rho_hat_p[(i * n - k + 1):(i * n)]))
        rho_hat_p_mean <- mean(rho_hat_p_avg)
        
        # Plot anualized amount of claims with standard deviations
        p <- ggplot() +
            geom_line(aes(x = percentile, y = mu_sorted_avg), color = "#000000", alpha = 0.7) +
            geom_point(aes(x = percentile, y = ClaimNb_avg), color = "#696969", size = 1, alpha = 0.7) +
            geom_line(aes(x = percentile, y = mu_sorted_avg + sqrt(sig_hat_p_avg / k)), color = "#5429ff", alpha = 0.7) +
            geom_line(aes(x = percentile, y = mu_sorted_avg - sqrt(sig_hat_p_avg / k)), color = "#5429ff", alpha = 0.7) +
            labs(title = "Prediction and Claim Frequency towards Percentile", x = "Percentile", y = "Average Claim Frequency")
        ggsave(paste("./results/", dataset, "/", model, "_plots/Anualized_Claim_Amount.png", sep = ""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

        # Plot Deviance Estimates with mean
        p <- ggplot() +
            geom_point(aes(x = percentile, y = rho_hat_p_avg), color = "#FF0000", alpha = 0.7) +
            geom_hline(yintercept = rho_hat_p_mean, color = "#FF0000", linetype = "dashed") +
            labs(title = "Deviance Estimates with Mean", x = "Percentile", y = "Deviance Estimates")
        ggsave(paste("./results/", dataset, "/", model, "_plots/Deviance_Estimates.png", sep = ""), plot = p, width = 10, height = 6, units = "in", dpi = 300)

        # Create dataframe and save to csv file in the DataSets folder
        q_df <- data.frame(percentile, mu_sorted_avg, ClaimNb_avg, sig_hat_p_avg, rho_hat_p_avg, W_avg)
        write.csv(q_df, paste("./results/", dataset, "/", model, "_percentile_data.csv", sep = ""), row.names = FALSE)


        # Create Model Performance metrics
        earned <- mu_sorted*W_sorted
        # Add other metrics
        mse <- mean((earned - ClaimNb_sorted)^2)
        mae <- mean(abs(earned - ClaimNb_sorted))
        log_fix <- function(x) ifelse(x == 0, 0, log(x))
        deviance <- mean(2 * (ClaimNb_sorted * log_fix(ClaimNb_sorted / earned) - ClaimNb_sorted + earned))
        q75 <- quantile(abs(earned - ClaimNb_sorted), 0.75)
        q90 <- quantile(abs(earned - ClaimNb_sorted), 0.90)
        q95 <- quantile(abs(earned - ClaimNb_sorted), 0.95)
        q99 <- quantile(abs(earned - ClaimNb_sorted), 0.95)
        q100<- quantile(abs(earned - ClaimNb_sorted), 1)
        q50 <- quantile(abs(earned - ClaimNb_sorted), 0.50)
        q25 <- quantile(abs(earned - ClaimNb_sorted), 0.25)
        q10 <- quantile(abs(earned - ClaimNb_sorted), 0.10)
        q05 <- quantile(abs(earned - ClaimNb_sorted), 0.05)

        # Save metrics
        metrics <- data.frame(model = model, n=n,deviance = deviance, mse = mse, mae = mae, q05=q05,q10=q10, q25 = q25, q75 = q75, q90 = q90, q95 = q95, q99 = q99, q100 = q100)
        write.csv(metrics, paste("./results/", dataset, "/", model, "_metrics.csv", sep = ""), row.names = FALSE)

    }

    
}
