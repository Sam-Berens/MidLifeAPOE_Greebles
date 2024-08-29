# Define the inverse link function (the logistic function) ...
# ... used to convert linear estimates into probability estimates
logistic <- function(l) {
  p <- 1 / (1 + exp(-l))
  return(p)
}

# Define a function that will return RT estimates
get_pc_est <- function(x, bayesian) {
  if (!bayesian) {
    # --- Frequentest model ---

    # Get the parameter estimates
    b <- matrix(as.numeric(fixef(freq06c)))
    # Get the number of parameters
    q <- length(b)
    # Get the parameter covariances
    c <- matrix(as.numeric(vcov(freq06c)), q, q)

    # Compute the linear estimates for each condition
    lin_estim <- x %*% b

    # Compute the linear standard errors for each condition
    lin_stder <- sqrt(diag(x %*% c %*% t(x)))

    # Get the degrees of freedom
    df <- as.numeric(summary(freq06c)$AICtab[5])

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
    post_samps <- as_draws_matrix(bayes06c)

    # Remove all the stuff we don't need to leave only fixed-effects estimates
    post_samps <- post_samps[, 1:16]

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
x <- model.matrix(freq06c)

# --- Get the condition codes for zAge, e4, sex ---
x_zage <- sort(unique(x[, 3]))
x_e4 <- sort(unique(x[, 4]))
x_sex <- sort(unique(x[, 5]))

# --- Specify the condition parameter weights for each condition ---
# The 8 numbers are the 8 model coefficients.
# Have a look at 'x' to see how they are derived.
# We only specify the first 4 terms and then...
# ... work out the interactions by multiplying main contrasts.
# Check the order of model coefficients (Intercept, hAmbi, zAge, E4 Female)
# e33
x_lo_age1_e33_male <- t(c(1, 0, x_zage[1], x_e4[1], x_sex[1]))
x_hi_age1_e33_male <- t(c(1, 1, x_zage[1], x_e4[1], x_sex[1]))
x_lo_age2_e33_male <- t(c(1, 0, x_zage[2], x_e4[1], x_sex[1]))
x_hi_age2_e33_male <- t(c(1, 1, x_zage[2], x_e4[1], x_sex[1]))
x_lo_age3_e33_male <- t(c(1, 0, x_zage[3], x_e4[1], x_sex[1]))
x_hi_age3_e33_male <- t(c(1, 1, x_zage[3], x_e4[1], x_sex[1]))
x_lo_age4_e33_male <- t(c(1, 0, x_zage[4], x_e4[1], x_sex[1]))
x_hi_age4_e33_male <- t(c(1, 1, x_zage[4], x_e4[1], x_sex[1]))

x_lo_age1_e33_fema <- t(c(1, 0, x_zage[1], x_e4[1], x_sex[2]))
x_hi_age1_e33_fema <- t(c(1, 1, x_zage[1], x_e4[1], x_sex[2]))
x_lo_age2_e33_fema <- t(c(1, 0, x_zage[2], x_e4[1], x_sex[2]))
x_hi_age2_e33_fema <- t(c(1, 1, x_zage[2], x_e4[1], x_sex[2]))
x_lo_age3_e33_fema <- t(c(1, 0, x_zage[3], x_e4[1], x_sex[2]))
x_hi_age3_e33_fema <- t(c(1, 1, x_zage[3], x_e4[1], x_sex[2]))
x_lo_age4_e33_fema <- t(c(1, 0, x_zage[4], x_e4[1], x_sex[2]))
x_hi_age4_e33_fema <- t(c(1, 1, x_zage[4], x_e4[1], x_sex[2]))
# e4+
x_lo_age1_e4p_male <- t(c(1, 0, x_zage[1], x_e4[2], x_sex[1]))
x_hi_age1_e4p_male <- t(c(1, 1, x_zage[1], x_e4[2], x_sex[1]))
x_lo_age2_e4p_male <- t(c(1, 0, x_zage[2], x_e4[2], x_sex[1]))
x_hi_age2_e4p_male <- t(c(1, 1, x_zage[2], x_e4[2], x_sex[1]))
x_lo_age3_e4p_male <- t(c(1, 0, x_zage[3], x_e4[2], x_sex[1]))
x_hi_age3_e4p_male <- t(c(1, 1, x_zage[3], x_e4[2], x_sex[1]))
x_lo_age4_e4p_male <- t(c(1, 0, x_zage[4], x_e4[2], x_sex[1]))
x_hi_age4_e4p_male <- t(c(1, 1, x_zage[4], x_e4[2], x_sex[1]))

x_lo_age1_e4p_fema <- t(c(1, 0, x_zage[1], x_e4[2], x_sex[2]))
x_hi_age1_e4p_fema <- t(c(1, 1, x_zage[1], x_e4[2], x_sex[2]))
x_lo_age2_e4p_fema <- t(c(1, 0, x_zage[2], x_e4[2], x_sex[2]))
x_hi_age2_e4p_fema <- t(c(1, 1, x_zage[2], x_e4[2], x_sex[2]))
x_lo_age3_e4p_fema <- t(c(1, 0, x_zage[3], x_e4[2], x_sex[2]))
x_hi_age3_e4p_fema <- t(c(1, 1, x_zage[3], x_e4[2], x_sex[2]))
x_lo_age4_e4p_fema <- t(c(1, 0, x_zage[4], x_e4[2], x_sex[2]))
x_hi_age4_e4p_fema <- t(c(1, 1, x_zage[4], x_e4[2], x_sex[2]))
# all
x_lo_age1_male <- t(c(1, 0, x_zage[1], 0.5, x_sex[1]))
x_hi_age1_male <- t(c(1, 1, x_zage[1], 0.5, x_sex[1]))
x_lo_age2_male <- t(c(1, 0, x_zage[2], 0.5, x_sex[1]))
x_hi_age2_male <- t(c(1, 1, x_zage[2], 0.5, x_sex[1]))
x_lo_age3_male <- t(c(1, 0, x_zage[3], 0.5, x_sex[1]))
x_hi_age3_male <- t(c(1, 1, x_zage[3], 0.5, x_sex[1]))
x_lo_age4_male <- t(c(1, 0, x_zage[4], 0.5, x_sex[1]))
x_hi_age4_male <- t(c(1, 1, x_zage[4], 0.5, x_sex[1]))

x_lo_age1_fema <- t(c(1, 0, x_zage[1], 0.5, x_sex[2]))
x_hi_age1_fema <- t(c(1, 1, x_zage[1], 0.5, x_sex[2]))
x_lo_age2_fema <- t(c(1, 0, x_zage[2], 0.5, x_sex[2]))
x_hi_age2_fema <- t(c(1, 1, x_zage[2], 0.5, x_sex[2]))
x_lo_age3_fema <- t(c(1, 0, x_zage[3], 0.5, x_sex[2]))
x_hi_age3_fema <- t(c(1, 1, x_zage[3], 0.5, x_sex[2]))
x_lo_age4_fema <- t(c(1, 0, x_zage[4], 0.5, x_sex[2]))
x_hi_age4_fema <- t(c(1, 1, x_zage[4], 0.5, x_sex[2]))

x_lo_male <- t(c(1, 0, 0, 0.5, x_sex[1]))
x_hi_male <- t(c(1, 1, 0, 0.5, x_sex[1]))
x_lo_fema <- t(c(1, 0, 0, 0.5, x_sex[2]))
x_hi_fema <- t(c(1, 1, 0, 0.5, x_sex[2]))

# --- Put all the condition parameter weights into a single matrix ---
x <- matrix(c(
              x_lo_age1_e33_male,
              x_lo_age2_e33_male,
              x_lo_age3_e33_male,
              x_lo_age4_e33_male,
              x_lo_age1_e33_fema,
              x_lo_age2_e33_fema,
              x_lo_age3_e33_fema,
              x_lo_age4_e33_fema,
              x_lo_age1_e4p_male,
              x_lo_age2_e4p_male,
              x_lo_age3_e4p_male,
              x_lo_age4_e4p_male,
              x_lo_age1_e4p_fema,
              x_lo_age2_e4p_fema,
              x_lo_age3_e4p_fema,
              x_lo_age4_e4p_fema,
              x_hi_age1_e33_male,
              x_hi_age2_e33_male,
              x_hi_age3_e33_male,
              x_hi_age4_e33_male,
              x_hi_age1_e33_fema,
              x_hi_age2_e33_fema,
              x_hi_age3_e33_fema,
              x_hi_age4_e33_fema,
              x_hi_age1_e4p_male,
              x_hi_age2_e4p_male,
              x_hi_age3_e4p_male,
              x_hi_age4_e4p_male,
              x_hi_age1_e4p_fema,
              x_hi_age2_e4p_fema,
              x_hi_age3_e4p_fema,
              x_hi_age4_e4p_fema,
              x_lo_age1_male,
              x_lo_age2_male,
              x_lo_age3_male,
              x_lo_age4_male,
              x_lo_age1_fema,
              x_lo_age2_fema,
              x_lo_age3_fema,
              x_lo_age4_fema,
              x_hi_age1_male,
              x_hi_age2_male,
              x_hi_age3_male,
              x_hi_age4_male,
              x_hi_age1_fema,
              x_hi_age2_fema,
              x_hi_age3_fema,
              x_hi_age4_fema,
              x_lo_male,
              x_lo_fema,
              x_hi_male,
              x_hi_fema), 52, 5, TRUE)
x <- cbind(x, matrix(nrow = 52, ncol = 11))

# --- Calculate the interaction terms ---
x[, 06] <- x[, 2] * x[, 3]
x[, 07] <- x[, 2] * x[, 4]
x[, 08] <- x[, 3] * x[, 4]
x[, 09] <- x[, 2] * x[, 5]
x[, 10] <- x[, 3] * x[, 5]
x[, 11] <- x[, 4] * x[, 5]
x[, 12] <- x[, 2] * x[, 3] * x[, 4]
x[, 13] <- x[, 2] * x[, 3] * x[, 5]
x[, 14] <- x[, 2] * x[, 4] * x[, 5]
x[, 15] <- x[, 3] * x[, 4] * x[, 5]
x[, 16] <- x[, 2] * x[, 3] * x[, 4] * x[, 5]

# --- Produce a data frame listing all the condition estimates ---
conditions <- data.frame(
                         Condition = c(rep("Low", 16), rep("High", 16), rep("Low", 8), rep("High", 8), rep("Low", 2), rep("High", 2)),
                         Genotype = c(
                                      rep("e33", 8),
                                      rep("e4p", 8),
                                      rep("e33", 8),
                                      rep("e4p", 8),
                                      rep("ave", 20)),
                         Sex = c(
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 4),
                           rep("Female", 4),
                           rep("Male", 1),
                           rep("Female", 1),
                           rep("Male", 1),
                           rep("Female", 1)),
                         Age = c(rep(c("45-49", "50-54", "55-59", "60-65"), 12), rep(c("mean"),4)))
conditions <- cbind(conditions, as.data.frame(get_pc_est(x, FALSE)))
# bayes_condest <- get_pc_est(x, TRUE)
# conditions <- cbind(conditions,
#                     as.data.frame(within(bayes_condest, rm(Bayes.Samps))))
print(conditions)

# # --- Make a quick plot ---
# conditions$Trial <- c(
#                       "Lo_e33_45-49",
#                       "Lo_e33_50-54",
#                       "Lo_e33_55-59",
#                       "Lo_e33_60-65",
#                       "Lo_e4p_45-49",
#                       "Lo_e4p_50-54",
#                       "Lo_e4p_55-59",
#                       "Lo_e4p_60-65",
#                       "Hi_e33_45-49",
#                       "Hi_e33_50-54",
#                       "Hi_e33_55-59",
#                       "Hi_e33_60-65",
#                       "Hi_e4p_45-49",
#                       "Hi_e4p_50-54",
#                       "Hi_e4p_55-59",
#                       "Hi_e4p_60-65",
#                       "Lo_all_45-49",
#                       "Lo_all_50-54",
#                       "Lo_all_55-59",
#                       "Lo_all_60-65",
#                       "Hi_all_45-49",
#                       "Hi_all_50-54",
#                       "Hi_all_55-59",
#                       "Hi_all_60-65")
# plot01a <- ggplot(conditions, aes(x = Trial, y = Frequ.Est, fill = Genotype)) +
#   geom_bar(stat = "identity", color = "black", position = position_dodge()) +
#   geom_errorbar(aes(ymin = Frequ.Low, ymax = Frequ.Hig),
#                 width = .2, position = position_dodge(.9)) +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   ylab("Model Estimates - Pr(Correct)")
# print(plot01a)