comp_baye_cons <- function(h, bayes) {

  # Get the posterior draws from the full (saturated) Bayesian model
  post_samps <- as_draws_matrix(bayes)

  # Remove all the stuff we don't need to leave only fixed-effects estimates
  nfixed = nrow(summary(bayes)$fixed)
  post_samps <- post_samps[, 1:nfixed]

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
      mean_post <- as.matrix(colMeans(lin_samps))
      unit <- mean_post / sqrt(sum(mean_post^2))
      dot <- lin_samps %*% unit

      post_mean <- mean(dot)
      post_medi <- quantile(dot, .5)
      post_025 <- quantile(dot, .025)
      post_975 <- quantile(dot, .975)
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