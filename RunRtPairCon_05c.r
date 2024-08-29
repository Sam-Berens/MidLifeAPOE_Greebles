# --- Frequentest contrast function ---
comp_freq_cons <- function(h) {

  # Get the parameter estimates
  b <- matrix(as.numeric(fixef(freq05c)))
  # Get the number of parameters
  q <- length(b)
  # Get the parameter covariances
  c <- matrix(as.numeric(vcov(freq05c)), q, q)
  # Get the degrees of freedom
  df <- as.numeric(summary(freq05c)$AICtab[5])

  # Pre-allocate a data frame for the stats
  contrasts <- data.frame()

  # Loop trough the contrasts to estimate the stats
  for (con_name in names(h)) {
    con_mat <- h[[con_name]]
    name <- con_name
    if (dim(con_mat)[1] == 1) {
      test_type <- sprintf("t(%i)", df)
      test_stat <- (con_mat %*% b) / sqrt(con_mat %*% c %*% t(con_mat))
      pvalue <- 2 * pt(-abs(test_stat), df = df)

    } else {
      df1 <- rankMatrix(con_mat)[1]
      test_type <- sprintf("F(%i,%i)", df1, df)
      test_stat <- (t(con_mat %*% b) %*% (solve(con_mat %*% c %*% t(con_mat))) %*% (con_mat %*% b)) / df1 # nolint: line_length_linter.
      pvalue <- 1 -  pf(test_stat, df1, df)
    }
    contrasts <- rbind(contrasts, data.frame(name,
                                            test_type,
                                            test_stat,
                                            pvalue))
  }

  # Return the result
  return(contrasts)
}

# --- Bayesian contrast function ---
comp_baye_cons <- function(h) {

  # Get the posterior draws from the full (saturated) Bayesian model
  post_samps <- as_draws_matrix(bayes05c)

  # Remove all the stuff we don't need to leave only fixed-effects estimates
  post_samps <- post_samps[, 1:8]

  # Pre-allocate a data frame for the stats
  contrasts <- data.frame()

  # Loop trough the contrasts to estimate the stats
  for (con_name in names(h)) {
    con_mat <- h[[con_name]]
    name <- con_name
    if (dim(con_mat)[1] == 1) {

      # Compute the posterior samples
      lin_samps <- t(con_mat %*% t(post_samps))

      # Compute the posterior mean
      post_mean <- mean(lin_samps)

      # Compute the posterior median
      post_medi <- quantile(lin_samps, .5)

      # Compute the 95% credible intervals
      post_025 <- quantile(lin_samps, .025)
      post_975 <- quantile(lin_samps, .975)

    } else {
      # Compute the posterior samples
      lin_samps <- t(con_mat %*% t(post_samps))
      l2 <- as.matrix(sqrt(rowSums(lin_samps^2)))

      post_mean <- mean(l2)
      post_medi <- quantile(l2, .5)
      post_025 <- quantile(l2, .025)
      post_975 <- quantile(l2, .975)
    }

    contrasts <- rbind(contrasts,
                       data.frame(name,
                                  post_mean,
                                  post_medi,
                                  post_025,
                                  post_975))
  }

  # Return the result
  return(contrasts)
}

# --- Get the design matrix ---
x <- model.matrix(freq05c)

# --- Get the condition codes for zAge and lDose ---
x_zage <- sort(unique(x[, 3]))
x_e4 <- sort(unique(x[, 4]))

# --- Specify the condition parameter weights for each condition ---
# The 8 numbers are the 8 model coefficients.
# Have a look at 'x' to see how they are derived.
# We only specify the first 4 terms and then...
# ... work out the interactions by multiplying main contrasts.
# Check the order of model coefficients (Intercept, hAmbi, zAge, lDose)
# e33
x_lo_age1_e33 <- t(c(1, 0, x_zage[1], x_e4[1]))
x_hi_age1_e33 <- t(c(1, 1, x_zage[1], x_e4[1]))
x_lo_age2_e33 <- t(c(1, 0, x_zage[2], x_e4[1]))
x_hi_age2_e33 <- t(c(1, 1, x_zage[2], x_e4[1]))
x_lo_age3_e33 <- t(c(1, 0, x_zage[3], x_e4[1]))
x_hi_age3_e33 <- t(c(1, 1, x_zage[3], x_e4[1]))
x_lo_age4_e33 <- t(c(1, 0, x_zage[4], x_e4[1]))
x_hi_age4_e33 <- t(c(1, 1, x_zage[4], x_e4[1]))
# e4p
x_lo_age1_e4p <- t(c(1, 0, x_zage[1], x_e4[2]))
x_hi_age1_e4p <- t(c(1, 1, x_zage[1], x_e4[2]))
x_lo_age2_e4p <- t(c(1, 0, x_zage[2], x_e4[2]))
x_hi_age2_e4p <- t(c(1, 1, x_zage[2], x_e4[2]))
x_lo_age3_e4p <- t(c(1, 0, x_zage[3], x_e4[2]))
x_hi_age3_e4p <- t(c(1, 1, x_zage[3], x_e4[2]))
x_lo_age4_e4p <- t(c(1, 0, x_zage[4], x_e4[2]))
x_hi_age4_e4p <- t(c(1, 1, x_zage[4], x_e4[2]))

# --- Put all the condition parameter weights into a single matrix ---
x <- matrix(c(
              x_lo_age1_e33,
              x_lo_age2_e33,
              x_lo_age3_e33,
              x_lo_age4_e33,
              x_lo_age1_e4p,
              x_lo_age2_e4p,
              x_lo_age3_e4p,
              x_lo_age4_e4p,
              x_hi_age1_e33,
              x_hi_age2_e33,
              x_hi_age3_e33,
              x_hi_age4_e33,
              x_hi_age1_e4p,
              x_hi_age2_e4p,
              x_hi_age3_e4p,
              x_hi_age4_e4p), 16, 4, TRUE)
x <- cbind(x, matrix(nrow = 16, ncol = 4))

# --- Calcuate the interaction terms ---
x[, 5] <- x[, 2] * x[, 3]
x[, 6] <- x[, 2] * x[, 4]
x[, 7] <- x[, 3] * x[, 4]
x[, 8] <- x[, 2] * x[, 3] * x[, 4]

# --- Define the lDose contrasts by age-group
h_e4_4age1 <- c(-1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0) %*% x # nolint: line_length_linter.
h_e4_4age2 <- c(0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0) %*% x # nolint: line_length_linter.
h_e4_4age3 <- c(0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0) %*% x # nolint: line_length_linter.
h_e4_4age4 <- c(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1) %*% x # nolint: line_length_linter.

# --- Package the contrasts into a list ---
h <- list(
  E4_4Age1 = h_e4_4age1,
  E4_4Age2 = h_e4_4age2,
  E4_4Age3 = h_e4_4age3,
  E4_4Age4 = h_e4_4age4)

# --- Run the Frequentist contrasts ---
freq_cons <- comp_freq_cons(h)
print(format(freq_cons, scientific = F))

# --- Run the Bayesian contrasts ---
baye_cons <- comp_baye_cons(h)
print(baye_cons)