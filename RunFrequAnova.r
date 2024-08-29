run_freq_anova <- function(frequ) {

  # --- Run the ANOVA ---
  anova_tab <- anova(frequ)

  # --- Get the p values ---
  p <- 1 - pf(anova_tab$`F value`,
              anova_tab$npar,
              summary(frequ)$AICtab[5])

  # --- Add the p values ---
  anova_tab$pValue <- p

  # --- Return the result ---
  return(anova_tab)
}