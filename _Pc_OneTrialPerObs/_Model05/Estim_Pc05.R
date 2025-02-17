# --- Clear all variables ---
rm(list = ls())

# --- Set whether we do a sloppy job of estimating the Bayesian models ---
is_dry_run <- FALSE

# --- Import libraries ---
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(brms)
library(BayesFactor)
library(emmeans)
library(ggplot2)
library(gridExtra)

# --- Load helper functions ---
source("../../RunFrequAnova.r")
source("../../RunFrequCons.r")
source("../../RunBayesCons.r")

# --- Description of dataframes ---
# data00 : Raw data
# data01 : Variables have been recoded, renamed and trimmed
# data02 : Long format
# data03 : One row per trial

# --- Read in the data ---
data00 <- read.csv("../../Greebles_SummaryData.csv")

# --- Select out key variables of interest  ---
cols2sel <- c(2, 7, 14:18, 20:21)
data01 <- data00[, cols2sel]

# --- Recoding variables  ---
data01$PpantId <- as.factor(data01$PpantId)
data01$Female <- recode(data01$gender, "Male" = 0, "Female" = 1)
data01$Age <- recode(data01$age_band,
                     "45-49" = 0, "50-54" = 1, "55-59" = 2, "60-65" = 3)
data01$E4 <- recode(data01$e4, "e3" = 0, "e4" = 1)
data01$lDose <- as.numeric(data01$dose)
data01$fDose <- as.factor(data01$dose)
data01$n_l <- data01$LowAmbiguity_n
data01$n_h <- data01$HighAmbiguity_n
data01$k_l <- data01$LowAmbiguity_accuracy
data01$k_h <- data01$HighAmbiguity_accuracy

# --- Delete unused varibales ---
data01 <- subset(data01, select = -c(gender,
                                     age_band,
                                     e4,
                                     dose,
                                     LowAmbiguity_n,
                                     HighAmbiguity_n,
                                     LowAmbiguity_accuracy,
                                     HighAmbiguity_accuracy))

# --- Z-score age (etc) ---
data01$Age <- (data01$Age - mean(data01$Age)) / sd(data01$Age)
names(data01)[names(data01) == "Age"] <- "zAge"

# --- Put the data table into long format ---
cols2sel <- c(1:6)
data02 <- rbind.data.frame(data01[, cols2sel], data01[, cols2sel])
data02$hAmbi <- c(rep(0, length(data01$PpantId)),
                  rep(1, length(data01$PpantId)))
data02$n <- c(data01$n_l, data01$n_h)
data02$k <- c(data01$k_l, data01$k_h)

# --- Make data03 (same as data02 but with one row per trial) ---
data03 <- data.frame()
num_rows <- nrow(data02)
iIn <- 0
for (ii in seq_len(num_rows)) {
  ppantId <- data02[ii, "PpantId"]
  female <- data02[ii, "Female"]
  zAge <- data02[ii, "zAge"]
  e4 <- data02[ii, "E4"]
  lDose <- data02[ii, "lDose"]
  fDose <- data02[ii, "fDose"]
  hAmbi <- data02[ii, "hAmbi"]
  n <- data02[ii, "n"]
  k <- data02[ii, "k"]
  for (jj in seq_len(n)) {
    iIn <- iIn + 1
    data03[iIn, "PpantId"] <- ppantId
    data03[iIn, "Female"] <- female
    data03[iIn, "zAge"] <- zAge
    data03[iIn, "E4"] <- e4
    data03[iIn, "lDose"] <- lDose
    data03[iIn, "fDose"] <- fDose
    data03[iIn, "hAmbi"] <- hAmbi
    data03[iIn, "n"] <- 1
    if (k >= jj) {
      data03[iIn, "k"] <- 1
    } else {
      data03[iIn, "k"] <- 0
    }
  }
}

# --- Specify the models of interest ---
formula05c <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * E4) + (1 | PpantId)
formula05l <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * lDose) + (1 | PpantId)
formula05f <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * fDose) + (1 | PpantId)

# --- Set the frequentist estimation options ---
control_opts <-
  glmerControl(
               # Generic stuff
               optimizer = c("bobyqa"),
               # "bobyqa", "Nelder_Mead", "optimx"
               restart_edge = FALSE,
               boundary.tol = 1e-5,
               calc.derivs = TRUE,
               use.last.params = FALSE,
               sparseX = FALSE,
               standardize.X = FALSE,

               # Convergence checking options
               check.conv.grad =
               .makeCC("warning", tol = 1e-3, relTol = NULL),
               check.conv.singular =
               .makeCC(action = "message", tol = formals(isSingular)$tol),
               check.conv.hess = .makeCC(action = "warning", tol = 1e-6),

               # Optimizer args
               optCtrl = list(maxfun = 200000),
               # optCtrl = list(method='nlminb'), # nolint
               # maxfun = 10000000, Added method, "nlminb"
               mod.type = "glmer",
               tolPwrss = 1e-7,
               compDev = TRUE,
               nAGQ0initStep = TRUE,
               check.response.not.const = "stop")

# --- Estimate the Frequentist models ---
freq05c <- glmer(formula = formula05c,
                 data = data03, family = binomial(link = "logit"),
                 control = control_opts)
freq05l <- glmer(formula = formula05l,
                 data = data03, family = binomial(link = "logit"),
                 control = control_opts)
freq05f <- glmer(formula = formula05f,
                 data = data03, family = binomial(link = "logit"),
                 control = control_opts)

# --- Frequentist ANOVAs ---
anova_f05c <- run_freq_anova(freq05c)
anova_f05l <- run_freq_anova(freq05l)
anova_f05f <- run_freq_anova(freq05f)