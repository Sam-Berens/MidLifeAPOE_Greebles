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
cols2sel <- c(2, 10, 15, 17, 24:27)
data01 <- data00[, cols2sel]

# --- Recoding variables  ---
data01$PpantId <- as.factor(data01$PpantId)
data01$Female <- recode(data01$gender, "Male" = 0, "Female" = 1)
data01$Age <- recode(data01$age_band,
                     "45-49" = 0, "50-54" = 1, "55-59" = 2, "60-65" = 3)
data01$E4 <- recode(data01$e4, "e3" = 0, "e4" = 1)
data01$lDose <- as.numeric(data01$dose)
data01$fDose <- as.factor(data01$dose)
data01$hAmbi <- recode(data01$Ambiguity, "H" = 1, "L" = 0)
data01$Correct <- data01$accuracy
data01$Rt <- data01$RT / 1000

# --- Delete unused variables ---
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

#gender 
descstats_gender <- data03 %>% group_by(Female, hAmbi) %>% summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))
descstats_genderDose <- data03 %>% group_by(Female, lDose, hAmbi) %>% summarise(muRt = mean(Rt), sdRT = sd(Rt), medianRt = median(Rt))

# --- Specify the models of interest ---
formula06c <- Rt ~ 1 + (hAmbi * zAge * E4 * Female) + (1 | PpantId)
formula06l <- Rt ~ 1 + (hAmbi * zAge * lDose * Female) + (1 | PpantId)

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
freq06c <- glmer(formula = formula06c,
                 data = data03, family = Gamma(link = "log"),
                 control = control_opts)
freq06l <- glmer(formula = formula06l,
                 data = data03, family = Gamma(link = "log"),
                 control = control_opts)

# --- Frequentist ANOVAs ---
anova_f06c <- run_freq_anova(freq06c)
anova_f06l <- run_freq_anova(freq06l)

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

bayes06c <- brm(formula06c,
                data = data03, prior = priors,
                family = Gamma(link = "log"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)
bayes06l <- brm(formula06l,
                data = data03, prior = priors,
                family = Gamma(link = "log"),
                warmup = n_warmup, iter = n_iter,
                chains = 5, cores = 5,
                save_pars = save_pars(all = TRUE), sample_prior = TRUE,
                control = list(adapt_delta = 0.95),
                seed = 1729)


# --- Contrast vectors for main effects and interactions ---
# For models 06c and 06l
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
anova_b05c <- comp_baye_cons(h05cl, bayes05c)
anova_b05l <- comp_baye_cons(h05cl, bayes05l)