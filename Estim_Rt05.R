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
# data02 : Rows have been trimmed to include correct responses only
# data03 : Rows have been trimmed to remove long RTs

# --- Read in the data ---
data00 <- read.csv("Greebles_TrialLevelData.csv")

# --- Select out key variables of interest  ---
cols2sel <- c(2, 10, 15, 24:27)
data01 <- data00[, cols2sel]

# --- Recoding variables  ---
data01$PpantId <- as.factor(data01$PpantId)
data01$Age <- recode(data01$age_band,
                     "45-49" = 0, "50-54" = 1, "55-59" = 2, "60-65" = 3)
data01$E4 <- recode(data01$e4, "e3" = 0, "e4" = 1)
data01$lDose <- as.numeric(data01$dose)
data01$fDose <- as.factor(data01$dose)
data01$hAmbi <- recode(data01$Ambiguity, "H" = 1, "L" = 0)
data01$Correct <- data01$accuracy
data01$Rt <- data01$RT / 1000

# --- Delete unused varibales ---
data01 <- subset(data01, select = -c(Ambiguity,
                                     RT,
                                     age_band,
                                     e4,
                                     dose,
                                     accuracy))

# --- Select only the RTs of interest (Correct trials only) ---
data02 <- data01 %>% filter(Correct == 1)
data02 <- subset(data02, select = -c(Correct))

# --- Cleaning: Remove RTs more than 3sd above the mean (...) ---
# (... per individual, per condition)
sum_stats_rt <- data02 %>%
  group_by(PpantId, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRt = sd(Rt), upBound = muRt + 3 * sdRt)
data03 <- left_join(data02, sum_stats_rt)
filter_rt <- data03$Rt <= data03$upBound
sprintf("%f%% of all RT datapoints above mu+3sd cutoff.",
        (1 - mean(as.numeric(filter_rt))) * 100)
data03 <- data03 %>% filter(filter_rt)
data03 <- subset(data03, select = -c(muRt, sdRt, upBound))

# --- Visualise the raw distribution of RTs ---
raw_rt_plot02 <- ggplot(data02, aes(x = Rt)) +
  geom_histogram(color = "black", fill = "coral2") + xlim(0, 40) +
  theme_classic() + labs(title = "RTs (pre-clean)", x = "RT(ms)", y = "Count")
raw_rt_plot03 <- ggplot(data03, aes(x = Rt)) +
  geom_histogram(color = "black", fill = "coral2") +  xlim(0, 40) +
  theme_classic() + labs(title = "RTs (post-clean)", x = "RT(ms)", y = "Count")
grid.arrange(raw_rt_plot02, raw_rt_plot03, ncol = 1, nrow = 2)

# --- Z-score age (etc) ---
data03$Age <- (data03$Age - mean(data03$Age)) / sd(data03$Age)
names(data03)[names(data03) == "Age"] <- "zAge"

# --- Descriptive stats ---
preclean_descstats_e4 <- data02 %>%
  group_by(E4, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))
preclean_descstats_dose <- data02 %>%
  group_by(fDose, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))
preclean_descstats_doseage <- data02 %>%
  group_by(fDose, Age, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))

descstats_e4 <- data03 %>%
  group_by(E4, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))
descstats_dose <- data03 %>%
  group_by(fDose, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))
descstats_doseage <- data03 %>%
  group_by(fDose, zAge, hAmbi) %>%
  summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))

# --- Specify the models of interest ---
formula05c <- Rt ~ 1 + (hAmbi * zAge * E4) + (1 | PpantId)
formula05l <- Rt ~ 1 + (hAmbi * zAge * lDose) + (1 | PpantId)
formula05f <- Rt ~ 1 + (hAmbi * zAge * fDose) + (1 | PpantId)

# --- Set the frequentist estimation options ---
control_opts <-
  glmerControl(
               # Generic stuff
               optimizer = c("Nelder_Mead"),
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
               optCtrl = list(maxfun = 100000),
               # optCtrl = list(method='nlminb'), # nolint
               # maxfun = 10000000, Added method, "nlminb"
               mod.type = "glmer",
               tolPwrss = 1e-7,
               compDev = TRUE,
               nAGQ0initStep = TRUE,
               check.response.not.const = "stop")

# --- Estimate the frequentist models ---
freq05c <- glmer(formula = formula05c,
                 data = data03, family = Gamma(link = "log"),
                 control = control_opts)
freq05l <- glmer(formula = formula05l,
                 data = data03, family = Gamma(link = "log"),
                 control = control_opts)
freq05f <- glmer(formula = formula05f,
                 data = data03, family = Gamma(link = "log"),
                 control = control_opts)

# --- Frequentist ANOVAs ---
anova_f05c <- run_freq_anova(freq05c)
anova_f05l <- run_freq_anova(freq05l)
anova_f05f <- run_freq_anova(freq05f)

# --- Specify the within-model priors for the Bayesian models ---
priors <- c(
  set_prior("lognormal(0, 1)", class = "Intercept"),
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
bayes05c <- brm(formula05c,
                data = data03, prior = priors,
                family = Gamma(link = "log"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes05l <- brm(formula05l,
                data = data03, prior = priors,
                family = Gamma(link = "log"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes05f <- brm(formula05f,
                data = data03, prior = priors,
                family = Gamma(link = "log"),
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
c_e34me33 <- t(c(0, 0, 0,  2, 0, 0,  1, 0, 0, 0, 0, 0))
c_e44me34 <- t(c(0, 0, 0, -2, 2, 0, -1, 1, 0, 0, 0, 0))
c_e44me33 <- t(c(0, 0, 0,  0, 2, 0,  0, 1, 0, 0, 0, 0))

c_e33_hml <- t(c(0, 1, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0))
c_e34_hml <- t(c(0, 1, 0,  0, 0, 0,  1, 0, 0, 0, 0, 0))
c_e44_hml <- t(c(0, 1, 0,  0, 0, 0,  0, 1, 0, 0, 0, 0))

c_e34hml_e33hml <- t(c(0, 0, 0, 0, 0, 0,  1, 0, 0, 0, 0, 0))
c_e44hml_e34hml <- t(c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0))
c_e44hml_e33hml <- t(c(0, 0, 0, 0, 0, 0,  0, 1, 0, 0, 0, 0))

# Package
c05f <- list(
  E34_vs_E33 = c_e34me33,
  E44_vs_E34 = c_e44me34,
  E44_vs_E33 = c_e44me33,
  E33_HmL = c_e33_hml,
  E34_HmL = c_e34_hml,
  E44_HmL = c_e44_hml,
  E34HmL_vs_E33HmL = c_e34hml_e33hml,
  E44HmL_vs_E34HmL = c_e44hml_e34hml,
  E44HmL_vs_E33HmL = c_e44hml_e33hml
)

# --- Run the contrasts ---
paircon_f05f <- comp_freq_cons(c05f, freq05f)
paircon_b05f <- comp_baye_cons(c05f, bayes05f)