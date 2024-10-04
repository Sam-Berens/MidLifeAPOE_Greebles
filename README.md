The following sections provide a guide to identifying and understanding different analyses within this repository. All analyses involve estimating mixed-effects models on behavioural observations. The structure of each statistical model can be inferred from how it is named (see sections 1-3 below). Section 4 details how the behavioural data are stored and provides an index of other relevant files that are needed for the analysis.

# 0: Model naming convention

Each analysis script will only ever include analyses of one outcome (dependent) variable. The identity of this variable can be determined from the filename of the analysis script which will contain one of the two following strings:

- “Pc” for the probability of a correct response across all perceptual discrimination trials (i.e., accuracy);
- “Rt” for the response time to perceptual discrimination trials that were responded to correctly (i.e., response latency);

Within all analysis scripts, specific statistical models are identified by a string of characters composed of three substrings: “[type][number][coding]”, where:

- “[type]” indicates whether a model is estimated via frequentist means ([type]=”freq”, sometimes abbreviated to “f”) or via Bayesian methods ([type]=”bayes”, sometimes abbreviated to “b”);
- “[number]” refers to the model number and indicates the overall structure of fixed- and random-effects predators in the model (see section 2);
- “[coding]” indicates how genotype information about different participants is encoded within the model’s fixed-effects predators (i.e., either as E4 ‘carrier status’, or E4 ‘haplotype’, see section 3);

# 1: Model type

## 1.1: Frequentist models (”freq”)

Frequentist models and inferential statistics are estimated via the lme4 package in R.

## 1.2: Bayesian models (”bayes”)

Bayesian models and inferential statistics are estimated via the BRMS package in R (via RStan), and custom-written functions (included in this repository, see section 4).

# 2: Model number

All analyses may be identified by one of two model numbers: “05” and “06”. Models identified by these numbers may involve different outcome (dependent) variables, yet models with the same number will all have the same fixed- and random-effects predators (independent variables):

- “05”: y ~ 1 + (hAmbi * zAge * E4) + (1 | PpantId)
- “06”: y ~  1 + (hAmbi * zAge * E4 * Female) + (1 | PpantId)

Where:

- “y” is the outcome variable;
- “1” represents the intercept term;
- “hAmbi” is a binary predictor indicating whether an observation relates to the high ambiguity condition (hAmbi =1) or not (hAmbi=0);
- “zAge” is a mean-centred predictor encoding the age category of different participants;
- “E4” is either a binary, linear, or multinomial (categorical) predictor encoding E4 genotype information about different participants (see section 3).
- “PpantId” is a multinomial (categorical) predictor denoting the identity of each participant;

# 3: Model coding

All statistical models within this repository include predictor variables encoding the APOE genotype of participants. However, this information is encoded in different ways across different models. The model coding identifier indicates how this genetic information is represented in the model. Specifically, the model coding identifier can take on three different values indicating the following:

- “c” (for ‘carrier’) refers to models that denote whether participants have either zero, or at least one APOE e4 allele. As such, the “E4” predictor within these models is a binary variable;
- “l” (for ‘linear’) refers to models that encode APOE e4 genotype in a linear, dose-dependent manner reflecting the number of e4 alleles carried by each participant;
- “f” (for ‘factorial’) refers to models that encode all APOE genotypes in the sample (”E3/E3”, “E3/E4”, “E4/E4”) as separate categories via two dummy-coded predictors;

# 4: Index of files and directories

## 4.2: Data files

- “Greebles_SummaryData.csv” : A CSV file containing summary performance statistics for each participant in a wide data format (one row per participant). This is used to model the probability of a correct response across conditions/predictors of interest;
- “Greebles_TrialLevelData.csv” : A CSV file containing trial-by-trial performance statistics for each participant in a long data format (one row per trial). This is used to model response times across conditions/predictors of interest;

## 4.2: Executable script files

- “Estim_Pc05.R” : An R script that estimates all models predicting the probability of a correct response with model number “05”;
- “Estim_Pc06.R” : An R script that estimates all models predicting the probability of a correct response with model number “06”;
- “Estim_Rt05.R” : An R script that estimates all models predicting response time with model number “05”;
- “Estim_Rt06.R” : An R script that estimates all models predicting response time with model number “06”;

## 4.3: Helper script files

These script files contain functions that are invoked by other analysis scripts to perform specific computations:

- GetPcEstims_05c.r;
- GetPcEstims_05f.r;
- GetPcEstims_05l.r;
- GetPcEstims_06c.r;
- GetRtEstims_05c.r;
- GetRtEstims_05f.r;
- GetRtEstims_05l.r;
- GetRtEstims_06c.r;
- RunBayesCons.r;
- RunFrequAnova.r;
- RunFrequCons.r;
- RunRtPairCon_05c.r;
- RunRtPairCon_05l.r;

## 4.4: Directories containing analysis outputs

These directories contain various output files derived from the analyses including plots and CSV files of model parameter estimates and inferential statistics. The plots are numbered in a way that is consistent with the published manuscript and supplementary information.

- /_Model05/;
- /_Model06/;
- Figure2;
- Figure3;
- Figure4;
- SupFigure1;
- SupFigure2;
- SupFigure3;
