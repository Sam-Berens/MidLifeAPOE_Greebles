# Define the inverse link function (the logistic function) ...
# ... used to convert linear estimates into probability estimates
logistic <- function(l) {
  p <- 1 / (1 + exp(-l))
  return(p)
}

# Define a function that will return pCorrect estimates
get_pc_est <- function(x, bayesian) {
  if (!bayesian) {
    # --- Frequentest model ---

    # Get the parameter estimates
    b <- matrix(as.numeric(fixef(freq05f)))
    # Get the number of parameters
    q <- length(b)
    # Get the parameter covariances
    c <- matrix(as.numeric(vcov(freq05f)), q, q)

    # Compute the linear estimates for each condition
    lin_estim <- x %*% b

    # Compute the linear standard errors for each condition
    lin_stder <- sqrt(diag(x %*% c %*% t(x)))

    # Get the degrees of freedom
    df <- as.numeric(summary(freq05f)$AICtab[5])

    # Compute the confidence intervals for each condition
    lin_ciwid <- lin_stder * qt(1 - 0.025, df)

    # Compute the RT estimates for each condition
    est <- logistic(lin_estim)
    # Compute  95% CIs on the probability estimates for each condition
    low <- logistic(lin_estim - lin_ciwid)
    hig <- logistic(lin_estim + lin_ciwid)

    # Return the results
    return(list("Frequ.Est" = est, "Frequ.Low" = low, "Frequ.Hig" = hig))

  } else {
    # --- Bayesian model ---

    # Get the posterior draws from the full (saturated) Bayesian model
    post_samps <- as_draws_matrix(bayes05f)

    # Remove all the stuff we don't need to leave only fixed-effects estimates
    post_samps <- post_samps[, 1:12]

    # Take the posterior mean as parameter estimates
    b <- colMeans(post_samps)
    # Compute the parameter covariances
    c <- cov(post_samps)

    # Compute the linear estimates for each condition
    lin_estim <- x %*% b

    # Compute the linear standard errors for each condition
    lin_stder <- sqrt(diag(x %*% c %*% t(x)))

    # Compute the (linear) posterior samples for each condition
    lin_samps <- x %*% t(post_samps)

    # Compute the 95% linear credible intervals for each condition
    lin_low <- apply(lin_samps, 1, function(xx) quantile(xx, .025))
    lin_hig <- apply(lin_samps, 1, function(xx) quantile(xx, .975))

    # Compute the probability estimates for each condition
    est <- logistic(lin_estim)

    # Compute  95% CIs on the probability estimates for each condition
    low <- logistic(lin_low)
    hig <- logistic(lin_hig)

    # Compute the posterior probability estimates for each condition
    samps <- logistic(lin_samps)

    # Return the results
    return(list("Bayes.Est" = est, "Bayes.Low" = low, "Bayes.Hig" = hig,
                "Bayes.Samps" = samps))
  }
}

# --- Get the design matrix ---
x <- model.matrix(freq05l)

# --- Get the condition codes for zAge ---
x_zage <- sort(unique(x[, 3]))

# --- Specify the condition parameter weights for each condition ---
# The 12 numbers are the 12 model coefficients.
# Have a look at 'x' to see how they are derived.
# We only specify the first 4 terms and then...
# ... work out the interactions by multiplying main contrasts.
# Check the order of model coefficients (Intercept, hAmbi, zAge, fDose1, fDose2)
# e33
x_lo_age1_e33 <- t(c(1, 0, x_zage[1], 0, 0))
x_hi_age1_e33 <- t(c(1, 1, x_zage[1], 0, 0))
x_lo_age2_e33 <- t(c(1, 0, x_zage[2], 0, 0))
x_hi_age2_e33 <- t(c(1, 1, x_zage[2], 0, 0))
x_lo_age3_e33 <- t(c(1, 0, x_zage[3], 0, 0))
x_hi_age3_e33 <- t(c(1, 1, x_zage[3], 0, 0))
x_lo_age4_e33 <- t(c(1, 0, x_zage[4], 0, 0))
x_hi_age4_e33 <- t(c(1, 1, x_zage[4], 0, 0))
# e34
x_lo_age1_e34 <- t(c(1, 0, x_zage[1], 1, 0))
x_hi_age1_e34 <- t(c(1, 1, x_zage[1], 1, 0))
x_lo_age2_e34 <- t(c(1, 0, x_zage[2], 1, 0))
x_hi_age2_e34 <- t(c(1, 1, x_zage[2], 1, 0))
x_lo_age3_e34 <- t(c(1, 0, x_zage[3], 1, 0))
x_hi_age3_e34 <- t(c(1, 1, x_zage[3], 1, 0))
x_lo_age4_e34 <- t(c(1, 0, x_zage[4], 1, 0))
x_hi_age4_e34 <- t(c(1, 1, x_zage[4], 1, 0))
# e44
x_lo_age1_e44 <- t(c(1, 0, x_zage[1], 0, 1))
x_hi_age1_e44 <- t(c(1, 1, x_zage[1], 0, 1))
x_lo_age2_e44 <- t(c(1, 0, x_zage[2], 0, 1))
x_hi_age2_e44 <- t(c(1, 1, x_zage[2], 0, 1))
x_lo_age3_e44 <- t(c(1, 0, x_zage[3], 0, 1))
x_hi_age3_e44 <- t(c(1, 1, x_zage[3], 0, 1))
x_lo_age4_e44 <- t(c(1, 0, x_zage[4], 0, 1))
x_hi_age4_e44 <- t(c(1, 1, x_zage[4], 0, 1))

# --- Put all the condition parameter weights into a single matrix ---
x <- matrix(c(
              x_lo_age1_e33,
              x_lo_age2_e33,
              x_lo_age3_e33,
              x_lo_age4_e33,
              x_lo_age1_e34,
              x_lo_age2_e34,
              x_lo_age3_e34,
              x_lo_age4_e34,
              x_lo_age1_e44,
              x_lo_age2_e44,
              x_lo_age3_e44,
              x_lo_age4_e44,
              x_hi_age1_e33,
              x_hi_age2_e33,
              x_hi_age3_e33,
              x_hi_age4_e33,
              x_hi_age1_e34,
              x_hi_age2_e34,
              x_hi_age3_e34,
              x_hi_age4_e34,
              x_hi_age1_e44,
              x_hi_age2_e44,
              x_hi_age3_e44,
              x_hi_age4_e44), 24, 5, TRUE)
x <- cbind(x, matrix(nrow = 24, ncol = 7))

# --- Calcuate the interaction terms ---
x[, 06] <- x[, 2] * x[, 3]
x[, 07] <- x[, 2] * x[, 4]
x[, 08] <- x[, 2] * x[, 5]
x[, 09] <- x[, 3] * x[, 4]
x[, 10] <- x[, 3] * x[, 5]
x[, 11] <- x[, 2] * x[, 3] * x[, 4]
x[, 12] <- x[, 2] * x[, 3] * x[, 5]

# --- Produce a data frame listing all the condition estimates ---
conditions <- data.frame(
                         Condition = c(rep("Low", 12), rep("High", 12)),
                         Genotype = c(
                                      rep("e33", 4),
                                      rep("e34", 4),
                                      rep("e44", 4),
                                      rep("e33", 4),
                                      rep("e34", 4),
                                      rep("e44", 4)),
                         Age = rep(c("45-49", "50-54", "55-59", "60-65"), 6))
conditions <- cbind(conditions, as.data.frame(get_pc_est(x, FALSE)))
bayes_condest <- get_pc_est(x, TRUE)
conditions <- cbind(conditions,
                    as.data.frame(within(bayes_condest, rm(Bayes.Samps))))
print(conditions)

# --- Make a quick plot ---
conditions$Trial <- c(
                      "Lo_e33_45-49",
                      "Lo_e33_50-54",
                      "Lo_e33_55-59",
                      "Lo_e33_60-65",
                      "Lo_e34_45-49",
                      "Lo_e34_50-54",
                      "Lo_e34_55-59",
                      "Lo_e34_60-65",
                      "Lo_e44_45-49",
                      "Lo_e44_50-54",
                      "Lo_e44_55-59",
                      "Lo_e44_60-65",
                      "Hi_e33_45-49",
                      "Hi_e33_50-54",
                      "Hi_e33_55-59",
                      "Hi_e33_60-65",
                      "Hi_e34_45-49",
                      "Hi_e34_50-54",
                      "Hi_e34_55-59",
                      "Hi_e34_60-65",
                      "Hi_e44_45-49",
                      "Hi_e44_50-54",
                      "Hi_e44_55-59",
                      "Hi_e44_60-65")
plot01a <- ggplot(conditions, aes(x = Trial, y = Frequ.Est, fill = Genotype)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = Frequ.Low, ymax = Frequ.Hig),
                width = .2, position = position_dodge(.9)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Model Estimates - Pr(Correct)")
print(plot01a)