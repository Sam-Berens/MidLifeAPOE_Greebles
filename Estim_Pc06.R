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

# --- descriptives by gender
gender <- data01 %>% group_by(Female) %>% summarise(k_lM = mean(k_l), k_hM = mean(k_h), k_lP = (k_lM/24)*100, k_hP = (k_hM/36)*100)
genderDose <- data01 %>% group_by(Female, fDose) %>% summarise(k_lM = mean(k_l), k_hM = mean(k_h), k_lP = (k_lM/24)*100, k_hP = (k_hM/36)*100)

# --- Put the data table into long format ---
cols2sel <- c(1:6)
data02 <- rbind.data.frame(data01[, cols2sel], data01[, cols2sel])
data02$hAmbi <- c(rep(0, length(data01$PpantId)),
                  rep(1, length(data01$PpantId)))
data02$n <- c(data01$n_l, data01$n_h)
data02$k <- c(data01$k_l, data01$k_h)

# --- Specify the models of interest ---

  # Including gender
formula06c <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * E4 * Female) + (1 | PpantId)
bormula06c <- k | trials(n) ~ 1 + (hAmbi * zAge * E4 * Female) + (1 | PpantId)
formula06f <- cbind(k, n - k) ~ 1 + (hAmbi * zAge * fDose * Female) + (1 | PpantId)
bormula06f <- k | trials(n) ~ 1 + (hAmbi * zAge * fDose * Female) + (1 | PpantId)

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
freq06c <- glmer(formula = formula06c,
                 data = data02, family = binomial(link = "logit"),
                 control = control_opts)
freq06f <- glmer(formula = formula06f,
                 data = data02, family = binomial(link = "logit"),
                 control = control_opts)

# --- Frequentist ANOVAs ---
anova_f06c <- run_freq_anova(freq06c)
anova_f06f <- run_freq_anova(freq06f)

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

bayes06c <- brm(bormula06c,
                data = data02, prior = priors,
                family = binomial(link = "logit"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes06f <- brm(bormula06f,
                data = data02, prior = priors,
                family = binomial(link = "logit"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)

# --- Contrast matrices for main effects and interactions (06c) --
h_Ambi <- t(c(0,4,0,0,0,0,2,0,2,0,0,0,0,1,0,0))
h_E4 <- t(c(0,0,0,4,0,0,2,0,0,0,2,0,0,1,0,0))
h_AmbiE4 <- t(c(0,0,0,0,0,0,2,0,0,0,0,0,0,1,0,0))
h_Sex <- t(c(0,0,0,0,4,0,0,0,2,0,2,0,0,1,0,0))
h_AmbiSex <- t(c(0,0,0,0,0,0,0,0,2,0,0,0,0,1,0,0))
h_E4Sex <- t(c(0,0,0,0,0,0,0,0,0,0,2,0,0,1,0,0))
h_AmbiE4Sex <- t(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0))

h_Age <- t(c(0,0,8,0,0,4,0,4,0,4,0,2,2,0,2,1))
h_AmbiAge <- t(c(0,0,0,0,0,4,0,0,0,0,0,2,2,0,0,1))
h_AgeE4 <- t(c(0,0,0,0,0,0,0,4,0,0,0,2,0,0,2,1))
h_AgeSex <- t(c(0,0,0,0,0,0,0,0,0,4,0,0,2,0,2,1))
h_AmbiAgeE4 <- t(c(0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,1))
h_AgeE4Sex <- t(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1))
h_AmbiAgeSex <- t(c(0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,1))
h_AmbiAgeE4Sex <- t(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))

h_Sex_Low <- t(c(0,0,0,0,2,0,0,0,0,0,1,0,0,0,0,0))
h_Sex_High <- t(c(0,0,0,0,2,0,0,0,2,0,1,0,0,1,0,0))
h_Ambi_Female <- t(c(0,2,0,0,0,0,1,0,2,0,0,0,0,1,0,0))
h_Ambi_Male <- t(c(0,2,0,0,0,0,1,0,0,0,0,0,0,0,0,0))

# Package
h06cl <- list(
  Ambi = h_Ambi,
  E4 = h_E4,
  AmbiE4 = h_AmbiE4,
  Sex = h_Sex,
  AmbiSex = h_AmbiSex,
  E4Sex = h_E4Sex,
  AmbiE4Sex = h_AmbiE4Sex,
  Age = h_Age,
  AmbiAge = h_AmbiAge,
  AgeE4 = h_AgeE4,
  AgeSex = h_AgeSex,
  AmbiAgeE4 = h_AmbiAgeE4,
  AgeE4Sex = h_AgeE4Sex,
  AmbiAgeSex = h_AmbiAgeSex,
  AmbiAgeE4Sex = h_AmbiAgeE4Sex,
  Sex_Low = h_Sex_Low,
  Sex_High = h_Sex_High,
  Ambi_Female = h_Ambi_Female,
  Ambi_Male = h_Ambi_Male
)

# --- Bayesian ANOVAs ---
anova_b06c <- comp_baye_cons(h06cl, bayes06c)