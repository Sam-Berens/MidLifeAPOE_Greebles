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
source("RunFrequAnova.r")
source("RunFrequCons.r")
source("RunBayesCons.r")

# --- Description of dataframes ---
# data00 : Raw data
# data01 : Variables have been recoded, renamed and trimmed
# data02 : Long format

# --- Read in the data ---
data00 <- read.csv("Greebles_SummaryData.csv")

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

# --- Specify the models of interest ---
formula05c <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * E4) + (1 | PpantId)
formula05l <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * lDose) + (1 | PpantId)
formula05f <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * fDose) + (1 | PpantId)

bormula05c <- k | trials(n) ~ 1 + (hAmbi * zAge * E4) + (1 | PpantId)
bormula05l <- k | trials(n) ~ 1 + (hAmbi * zAge * lDose) + (1 | PpantId)
bormula05f <- k | trials(n) ~ 1 + (hAmbi * zAge * fDose) + (1 | PpantId)

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
                 data = data02, family = binomial(link = "logit"),
                 control = control_opts)
freq05l <- glmer(formula = formula05l,
                 data = data02, family = binomial(link = "logit"),
                 control = control_opts)
freq05f <- glmer(formula = formula05f,
                 data = data02, family = binomial(link = "logit"),
                 control = control_opts)

# --- Frequentist ANOVAs ---
anova_f05c <- run_freq_anova(freq05c)
anova_f05l <- run_freq_anova(freq05l)
anova_f05f <- run_freq_anova(freq05f)

# --- Specify the within-model priors for the Bayesian models ---
priors <- c(
  set_prior("logistic(0, 1)", class = "Intercept"),
  set_prior("cauchy(0, sqrt(.5))", class = "b")
)

# --- Estimate the Bayesian models ---
if (is_dry_run) {
  n_warmup <- 100
  n_iter <- 1000
} else {
  n_warmup <- 1000
  n_iter <- 20000
}
bayes05c <- brm(bormula05c,
                data = data02, prior = priors,
                family = binomial(link = "logit"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes05l <- brm(bormula05l,
                data = data02, prior = priors,
                family = binomial(link = "logit"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes05f <- brm(bormula05f,
                data = data02, prior = priors,
                family = binomial(link = "logit"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)

# --- Contrast vectors for main effects and interactions ---
# For models 05c and 05l
h_ambi      <- t(c(0, 2, 0, 0, 0, 1, 0, 0))
h_e4        <- t(c(0, 0, 0, 2, 0, 1, 0, 0))
h_ambie4    <- t(c(0, 0, 0, 0, 0, 1, 0, 0))
h_age       <- t(c(0, 0, 4, 0, 2, 0, 2, 1))
h_ambiage   <- t(c(0, 0, 0, 0, 2, 0, 0, 1))
h_e4age     <- t(c(0, 0, 0, 0, 0, 0, 2, 1))
h_ambie4age <- t(c(0, 0, 0, 0, 0, 0, 0, 1))

# Package
h05cl <- list(
  Ambi = h_ambi,
  E4 = h_e4,
  AmbiE4 = h_ambie4,
  Age = h_age,
  AmbiAge = h_ambiage,
  E4Age = h_e4age,
  AmbiE4Age = h_ambie4age
)

# --- Contrast matrices for main effects and interactions ---
# For model 05f

############  1  2  3  4  5  6  7  8  9 10 11 12
h_ambi <- t(c(0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0))

#################  1  2  3  4  5  6  7  8  9 10 11 12  1  2  3  4  5  6  7  8  9 10 11 12 # nolint
h_e4 <- t(matrix(c(0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0), 12, 2)) # nolint

#####################  1  2  3  4  5  6  7  8  9 10 11 12  1  2  3  4  5  6  7  8  9 10 11 12 # nolint
h_ambie4 <- t(matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0), 12, 2)) # nolint

###########  1  2  3  4  5  6  7  8  9 10 11 12
h_age <- t(c(0, 0, 6, 0, 0, 3, 0, 0, 0, 0, 1, 1))

###############  1  2  3  4  5  6  7  8  9 10 11 12
h_ambiage <- t(c(0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1))

####################  1  2  3  4  5  6  7  8  9 10 11 12  1  2  3  4  5  6  7  8  9 10 11 12 # nolint
h_e4age <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1), 12, 2)) # nolint

########################  1  2  3  4  5  6  7  8  9 10 11 12  1  2  3  4  5  6  7  8  9 10 11 12 # nolint
h_ambie4age <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1), 12, 2)) # nolint

# Package
h05f <- list(
  Ambi = h_ambi,
  E4 = h_e4,
  AmbiE4 = h_ambie4,
  Age = h_age,
  AmbiAge = h_ambiage,
  E4Age = h_e4age,
  AmbiE4Age = h_ambie4age
)

# --- Bayesian ANOVAs ---
anova_b05c <- comp_baye_cons(h05cl, bayes05c)
anova_b05l <- comp_baye_cons(h05cl, bayes05l)
anova_b05f <- comp_baye_cons(h05f, bayes05f)

# --- Contrast matrices for pairwise contrasts ---
# For model 05f
c_e34me33 <- t(c(0, 0, 0,  2, 0, 0,  2, 0, 0, 0, 0, 0))
c_e44me34 <- t(c(0, 0, 0, -2, 2, 0, -2, 2, 0, 0, 0, 0))
c_e44me33 <- t(c(0, 0, 0,  0, 2, 0,  0, 2, 0, 0, 0, 0))

# Package
c05f <- list(
  E34_vs_E33 = c_e34me33,
  E44_vs_E34 = c_e44me34,
  E44_vs_E33 = c_e44me33
)

# --- Run the contrasts ---
paircon_f05f <- comp_freq_cons(c05f, freq05f)
paircon_b05f <- comp_baye_cons(c05f, bayes05f)