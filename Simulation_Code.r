################ Simulate and Measure Effect of Duration Distortion ################

# Load libraries
library(dplyr)
library(ggplot2)

# Parameters
max_intensity <- 0.1
min_intensity <- 0.02
n_groups <- 5
p_n <- 0.2
n <- 100000

# Create Functions for Simulation
data_generator <- function(n,n_groups,p,p_n){
    # Generalte policies
    policies = data.frame(
        policy_id = 1:n,
        risk_group = sample(1:n_groups, n, replace = TRUE)
    )
    policies$intensity = min_intensity + (max_intensity - min_intensity) * (policies$risk_group - 1) / (n_groups - 1)
    policies$risk_group <- as.factor(policies$risk_group)

    # Generate Natural Duration and Claims
    policies$duration = 1 - (runif(n) < p_n) * runif(n, 0, 1)
    policies$n_claims <- rpois(n, policies$intensity * policies$duration)

    # Introduce detrimental claims
    policies$duration_cens = policies$duration
    policies$n_claims_cens <- policies$n_claims
    for (i in 1:n) {
        if (policies$n_claims[i] == 0) {
            next
        }
        policy_id = policies$policy_id[i]
        claim_durations <- runif(policies$n_claims[i], 0, policies$duration[i])
        for (j in 1:policies$n_claims[i]) {
            claim_id = i + j
            claim_duration = claim_durations[j]
            cause_cancellation = runif(1) < p
            if (cause_cancellation) {
                policies$duration_cens[i] = claim_duration
                policies$n_claims_cens[i] = j
                break
            }
        }
    }
    return(policies)
}

eval_model <- function(model,df){
    # Calculate Estimates
    dispersion <- sum(residuals(model, type = "pearson")^2)/model$df.residual

    # Add coefficent estimates
    res <- data.frame(Dispersion = dispersion)
    for (i in 1:max(as.numeric(df$risk_group))) {
        res$B.est <- summary(model)$coefficients[i, 1]
        res$B.std <- summary(model)$coefficients[i, 2]
        res$B.true <-log(mean(df$intensity[df$risk_group == i]))
        colnames(res) <- c(colnames(res)[1:(length(res)-3)],paste("B", i, c("_est", "_std", "_true"),sep=""))
    }

    return(res)
}

########## Run Experiment #############
run_experiment <- function(n_experiments,n,n_groups,p,p_n){
    results = data.frame()
    for (i in 1:n_experiments) {
        policies = data_generator(n, n_groups, p, p_n)

        # Estimate parameters
        p_n_hat <- mean(policies$duration_cens[policies$n_claims_cens == 0] < 1)
        p_z_hat <- mean(policies$duration_cens[policies$n_claims_cens > 0] < 1)
        p_hat <- (p_z_hat - p_n_hat / 2) / (1 - p_n_hat / 2)

        # Calculate probability of being natural cancelled
        q_hat <- sapply(policies$n_claims_cens, function(x) ifelse(x > 0, (1 - p_hat), 0))
        q_hat[policies$duration_cens == 1] <- 0
        p_c_hat = sapply(policies$duration_cens, function(x) ifelse(x < 1, p_n_hat * x, 1))
        policies$prob_natural <- q_hat * p_c_hat

        # # Expected Full Duration
        p_n_remain = (1 - policies$duration_cens) * p_n_hat
        expected_remaining_duration <- (1 - p_n_remain) + p_n_remain * (policies$duration_cens + 1) / 2

        # Calculate duration adjusted for detrimental claims
        policies$dur_adj <- (policies$prob_natural * policies$duration_cens + (1 - policies$prob_natural) * expected_remaining_duration)

        # Adjust durations for to follow expected values
        policies$dur_adj <- policies$dur_adj * sum(policies$duration_cens) / sum(policies$dur_adj)

        # Estimate poisson model without Detrimental Claims
        model <- glm(n_claims ~ risk_group - 1, data = policies, family = poisson, offset = log(duration))
        res <- eval_model(model, policies)
        res$Duration_category <- "Without Detrimental Claims"
        results <- rbind(results, res)
        
        # Estimate poisson model with Detrimental Claims
        model_cens <- glm(n_claims_cens ~ risk_group - 1, data = policies, family = poisson, offset = log(duration_cens))
        res_cens <- eval_model(model_cens, policies)
        res_cens$Duration_category <- "With Detrimental Claims"
        results <- rbind(results, res_cens)
        
        # Estimate poisson model with Detrimental Claims and adjusted duration
        model_adj1 <- glm(n_claims_cens ~ risk_group - 1, data = policies, family = poisson, offset = log(dur_adj))
        res_adj1 <- eval_model(model_adj1, policies)
        res_adj1$Duration_category <- "Detrimental Claim adjusted Duration"

        results <- rbind(results, res_adj1)
    }
    return(results)
}

# Evaluate probability estimates
set.seed(1234)
tot_res <- data.frame()
for (p in 0:10/10) {
    for (i in 1:1000){
        results <- data_generator(n, n_groups, p,p_n)
        mean_z = mean(results$n_claims_cens[results$n_claims_cens > 0])
        p_n_hat <- mean(results$duration_cens[results$n_claims_cens == 0] < 1) 
        p_z_hat <- mean(results$duration_cens[results$n_claims_cens > 0] < 1) 
        p_hat <- (p_z_hat - p_n_hat / 2) / (1 - p_n_hat / 2)
        # Save phat results
        tot_res <- rbind(tot_res, data.frame(p_n_hat = p_n_hat, p_hat = p_hat, p = p, mean_z = mean_z))

    }
}
write.csv(tot_res, "./results/probability_estimates_results.csv")

# Run Experiment
set.seed(1234)
tot_res <- data.frame()
for (p in 0:100/100){
    results <- run_experiment(100, n, n_groups, p,p_n)
    results$p <- p
    tot_res <- rbind(tot_res, results)
    write.csv(tot_res, "./results/result_sim_estimates.csv")
}
