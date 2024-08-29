comp_freq_cons <- function(h, frequ) {

  # Get the parameter estimates
  b <- matrix(as.numeric(fixef(frequ)))

  # Get the number of parameters
  q <- length(b)

  # Get the parameter covariances
  cov <- matrix(as.numeric(vcov(frequ)), q, q)

  # Get the degrees of freedom
  df2 <- as.numeric(summary(frequ)$AICtab[5])

  # Pre-allocate a data frame for the stats
  contrasts <- data.frame()

  # Loop trough the contrasts to estimate the stats
  for (con_name in names(h)) {
    name <- con_name
    con_mat <- h[[con_name]]
    if (dim(con_mat)[1] == 1) {
      test_type <- sprintf("t(%i)", df2)
      test_stat <- (con_mat %*% b) / sqrt(con_mat %*% cov %*% t(con_mat))
      p_value <- 2 * pt(-abs(test_stat), df = df2)
    } else {
      df1 <- rankMatrix(con_mat)[1]
      test_type <- sprintf("F(%i,%i)", df1, df2)
      test_stat <- (t(con_mat %*% b) %*% 
                      (solve(con_mat %*% cov %*% t(con_mat))) %*%
                      (con_mat %*% b)) / df1
      p_value <- 1 - pf(test_stat, df1, df2)
    }
    contrasts <- rbind(contrasts,
                       data.frame(name, test_type, test_stat, p_value))
  }

  # Return the result
  return(contrasts)
}