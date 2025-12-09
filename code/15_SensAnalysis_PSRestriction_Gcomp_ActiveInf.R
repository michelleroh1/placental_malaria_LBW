
##################################################################################
##	  Quantifying the contributions of active and past placental malaria infection
##    to low birthweight in a high transmission area of Uganda
##    15_Sensitivity Analyses Restricting Analyses of Active Infection to 
##       Individuals with Common Support
##
##    This script is the same as 7_Gcomp_ActiveInf, but excludes those 
##    that fell outside the regions of good overlap of propensity scores 
##################################################################################

rm(list=ls())

# ----------------------------
#   User inputs: 
# ----------------------------
# Edit line below to set path directory user's own path directory: 
setwd('[insert path directory]') 

# Create "Results" folder if it doesn't already exist
if (!dir.exists("Results")) {
  dir.create("Results")
}


# ----------------------------
#   Upload libraries
# ----------------------------
library(readstata13)   
library(boot)   
library(tibble)
library(dplyr)
library(openxlsx)
library(ggplot2)


# ----------------------------
#   Upload databases
# ----------------------------

#  Upload dataset with indicator variable of the individuals to exclude based on PS
pscore <- read.xlsx('Databases/Created Databases/Propensity Scores_Indicator for Restricted Analyses.xlsx') 

# Merge with analytic dataset  
db <- read.dta13('Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta', nonint.factors = TRUE) %>%
  mutate(propfibrincat_sim_cat = case_when(
    propfibrincat_sim == 0 ~ "None", 
    propfibrincat_sim == 1 ~ "<10%", 
    propfibrincat_sim == 2 ~ "10 to <30%", 
    propfibrincat_sim == 3 ~ ">=30%")) %>%
  left_join(., pscore, by='id')            # merge in pscore database
  

# -------------------------------
#   Define outcomes for loop
# -------------------------------
outcomes <- c("preterm", "SGA", "LBWdich")


# ==================================================================
#   PART 1: OVERALL G-COMPUTATION (NO GRAVIDITY STRATIFICATION)
# ==================================================================

# ---------------------------------------------
#   G-computation function 
# ---------------------------------------------
# This function:
#   - Fits an adjusted logistic regression model for outcome_var
#   - If model with interaction activeHP*graviddich term does not converge,
#     falls back to a main-effects model
#   - Constructs counterfactual datasets to:
#       * estimate risks and effect estimates among the exposed (AFE)
#       * estimate population-level risks and effect estimates (PAF)
#   - Returns a named vector of:
#       margin0, margin1, margin2, margin3,
#       afe_rr, afe_rd, afe,
#       pae_rr, pae_rd, paf

overall_compute_metrics <- function(data, outcome_var, indices) {
  # Resample rows according to indices (needed for when bootstrapping)
  d <- data[indices, ]
  
  # -----------------------
  # 1. Fit outcome model
  # -----------------------
  # Full model with interaction between active placental infection and gravidity
  formula_full <- as.formula(
    paste(outcome_var, "~ activeHP*graviddich + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
  )
  
  model <- try(
    glm(formula_full, data = d, family = binomial(link = "logit")),
    silent = TRUE
  )
  
  # If the interaction model fails (non-convergence),
  # refit a simpler main-effects model without the interaction.
  if (inherits(model, "try-error")) {
    formula_maineffects <- as.formula(
      paste(
        outcome_var,
        "~ activeHP + graviddich + treatmentarm + age_enroll + age2",
        "+ wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2",
        "+ GAenroll + gaenroll2 + educ_cat + male"
      )
    )
    
    model <- glm(formula_maineffects, data = d, family = binomial(link = "logit"))
  }
  
  # ------------------------------------------------------------
  # 2. Define counterfactual scenarios for g-computation
  # ------------------------------------------------------------

  # ---- AFE: subset of women with activeHP_sim == 1 ----
  data0 <- d %>%
    filter(activeHP_sim==1) %>%
    mutate(activeHP=0,              # set exposure to 0 (no active infection) for model prediction
           treatmentarm="SP")       # set treatment arm to SP for all
  data1 <- d %>%
    filter(activeHP_sim==1) %>%
    mutate(activeHP=1,              # set exposure to 1 (active infection) for model prediction
           treatmentarm="SP")       # set treatment arm to SP for all
  
  # ---- PAE: full population under natural course exposure assignments ----
  data2 <- d %>%
    mutate(activeHP = 0,            # set exposure to 0 (no active infection) for model prediction 
           treatmentarm = "SP")
  data3 <- d %>%
    mutate(activeHP = activeHP_sim, # set exposure to simulated natural course exposure distribution
           treatmentarm = "SP")
  
  # -------------------------------------------------
  # 3. Compute marginal risks under each scenario
  # -------------------------------------------------
  # Predict outcome probabilities for each participant, then estimate risk
  margin0 <- mean(predict(model, newdata=data0, type="response"))
  margin1 <- mean(predict(model, newdata=data1, type="response"))
  margin2 <- mean(predict(model, newdata=data2, type="response"))
  margin3 <- mean(predict(model, newdata=data3, type="response"))
  
  # -------------------------------------------------
  # 4. Attributable effects among the exposed (AFE)
  # -------------------------------------------------
  # Risk ratio among exposed (RRexp)
  afe_rr = margin1 / margin0
  # Risk difference among exposed (RDexp)
  afe_rd = margin1 - margin0
  # Attributable fraction among the exposed: (AFE)
  afe <- (margin1 - margin0) / margin1
  
  # -------------------------------------------------
  # 5. Population attributable effects (PAF)
  # -------------------------------------------------
  # Population risk ratio (RRpop)
  pae_rr <- margin3 / margin2
  # Population risk difference (RDpop)
  pae_rd <- margin3 - margin2
  # Population attributable fraction (PAF)
  paf <- (margin3 - margin2) / margin3
  
  # Return all metrics as a named numeric vector
  return(c(margin0=margin0, margin1=margin1, margin2=margin2, margin3=margin3,
           afe_rr=afe_rr, afe_rd=afe_rd, afe=afe, 
           pae_rr=pae_rr, pae_rd=pae_rd, paf=paf))
}


# -----------------------------------------------------------
#   Loop g-computation function over all birth outcomes
# -----------------------------------------------------------
# For each outcome:
#   - Compute point estimates using the full dataset
#   - Bootstrap the metrics to obtain BCa 95% CIs
#   - Store estimates + CIs in results vector

# Restrict dataset to only those with propensity scores that overlapped
activeinf_overall <- db %>%
  filter(overall_actinf_pscore==1)

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Generate point estimates using full dataset  -------
  full_est <- overall_compute_metrics(activeinf_overall, y, 1:nrow(activeinf_overall))
  
  # --- 2. Bootstrap for BCa 95% CIs --------
  boot_obj <- boot(
    data = activeinf_overall,
    statistic = function(d, i) overall_compute_metrics(d, y, i),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_margin0 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_margin1 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  ci_margin2 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_margin3 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_afe_rr <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_afe_rd <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  ci_afe  <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  
  ci_pae_rr  <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pae_rd  <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_paf  <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  # --- 3. Combine point estimates and 95% CIs -------
  res <- c(
    
    margin0   = full_est["margin0"],
    margin0_lb = ci_margin0[1],
    margin0_ub = ci_margin0[2],
    
    margin1   = full_est["margin1"],
    margin1_lb = ci_margin1[1],
    margin1_ub = ci_margin1[2],
    
    margin2   = full_est["margin2"],
    margin2_lb = ci_margin2[1],
    margin2_ub = ci_margin2[2],
    
    margin3   = full_est["margin3"],
    margin3_lb = ci_margin3[1],
    margin3_ub = ci_margin3[2],
    
    afe_rr        = full_est["afe_rr"],
    afe_rr_lb     = ci_afe_rr[1],
    afe_rr_ub     = ci_afe_rr[2],
    
    afe_rd        = full_est["afe_rd"],
    afe_rd_lb     = ci_afe_rd[1],
    afe_rd_ub     = ci_afe_rd[2],
    
    afe        = full_est["afe"],
    afe_lb     = ci_afe[1],
    afe_ub     = ci_afe[2],
    
    pae_rr        = full_est["pae_rr"],
    pae_rr_lb     = ci_pae_rr[1],
    pae_rr_ub     = ci_pae_rr[2],
    
    pae_rd        = full_est["pae_rd"],
    pae_rd_lb     = ci_pae_rd[1],
    pae_rd_ub     = ci_pae_rd[2],
    
    paf        = full_est["paf"],
    paf_lb     = ci_paf[1],
    paf_ub     = ci_paf[2]
  )
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# ---------------------------------------------
#     Assemble results into a data frame
# ---------------------------------------------
# Combine results into a dataframe
overall_df <- do.call(rbind, lapply(names(results_list), function(outcome) {
  vec <- results_list[[outcome]]        # the named vector
  df <- as.data.frame(t(unname(vec)))   # remove names, transpose, make 1-row df
  df$outcome <- outcome                 # add outcome column
  df })) %>%
  mutate(subgroup = 'Overall')          # label that these are overall estimates

# Assign column names
colnames(overall_df) <- c(
  "margin0", "margin0_lb", "margin0_ub",
  "margin1", "margin1_lb", "margin1_ub",
  "margin2", "margin2_lb", "margin2_ub",
  "margin3", "margin3_lb", "margin3_ub",
  "afe_rr", "afe_rr_lb", "afe_rr_ub",
  "afe_rd", "afe_rd_lb", "afe_rd_ub",
  "afe", "afe_lb", "afe_ub",
  "pae_rr", "pae_rr_lb", "pae_rr_ub",
  "pae_rd", "pae_rd_lb", "pae_rd_ub",
  "paf", "paf_lb", "paf_ub",
  "outcome", "subgroup"
)

# Reorder so outcome and subgroup appear first
overall_df <- overall_df[, c("outcome", "subgroup", setdiff(colnames(overall_df), c("outcome", "subgroup")))]

rownames(overall_df) <- NULL
print(overall_df)


# ==================================================================
#   PART 2: OVERALL G-COMPUTATION (GRAVIDITY-STRATIFICATION)
# ==================================================================

# --------------------------------------------------------------
#   Gravidity-specific g-computation analyses
# --------------------------------------------------------------
# This section repeats the g-computation procedure separately for:
#   - Primigravidae
#   - Multigravidae
#
# For each stratum and each outcome, we estimate:
#   (1) Marginal risks under counterfactual exposure scenarios
#   (2) Attributable effects among the exposed (RRexp, RDexp, AFE)
#   (3) Population attributable effects (RRpop, RDpop, PAF)
#
#   Function for gravidity-specific g-computation: 

gravid_compute_metrics <- function(data, outcome_var, indices, grav_subgroup = NULL) {
  # Resample rows according to indices (needed for when bootstrapping)
  d <- data[indices, ]
  
  # Filter by gravidity subgroup 
  if (!is.null(grav_subgroup)) {
    d <- d[d$graviddich == grav_subgroup, ]
  }
  
  # -----------------------
  # 1. Fit outcome model
  # -----------------------
  formula_full <- as.formula(
    paste(outcome_var, "~ activeHP + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
  )
  
  model <- try(glm(formula_full, data=d, family=binomial(link="logit")), silent=TRUE)
  
  # ------------------------------------------------------------
  # 2. Define counterfactual scenarios for g-computation
  # ------------------------------------------------------------
  # ---- AFE: subset of women with activeHP_sim == 1 ----
  data0 <- d %>%
    filter(activeHP_sim==1) %>%
    mutate(activeHP=0,               # set exposure to 0 (no active infection) for model prediction
           treatmentarm="SP")        # set treatment arm to SP for all
  # All exposed among exposed
  data1 <- d %>%
    filter(activeHP_sim==1) %>%
    mutate(activeHP=1,               # set exposure to 1 (active infection) for model prediction
           treatmentarm="SP")        # set treatment arm to SP for all
  
  # ---- PAE: full population under natural course exposure assignments ----
  data2 <- d %>%
    mutate(activeHP = 0,             # set exposure to 0 (no active infection) for model prediction 
           treatmentarm = "SP")
  data3 <- d %>%
    mutate(activeHP = activeHP_sim,  # set exposure to simulated natural course exposure distribution
           treatmentarm = "SP")
  
  # -------------------------------------------------
  # 3. Compute marginal risks under each scenario
  # -------------------------------------------------
  # Predict outcome probabilities for each participant, then estimate risk
  margin0 <- mean(predict(model, newdata=data0, type="response")) # none exposed among exposed
  margin1 <- mean(predict(model, newdata=data1, type="response")) # all exposed among exposed
  margin2 <- mean(predict(model, newdata=data2, type="response")) # none exposed in pop
  margin3 <- mean(predict(model, newdata=data3, type="response")) # NC exposure
  
  # -------------------------------------------------
  # 4. Attributable effects among the exposed (AFE)
  # -------------------------------------------------
  # Risk ratio among exposed (RRexp)
  afe_rr = margin1 / margin0
  # Risk difference among exposed (RDexp)
  afe_rd = margin1 - margin0
  # Attributable fraction among the exposed: (AFE)
  afe <- (margin1 - margin0) / margin1
  
  # -------------------------------------------------
  # 5. Population attributable effects (PAF)
  # -------------------------------------------------
  # Population risk ratio (RRpop)
  pae_rr <- margin3 / margin2
  # Population risk difference (RDpop)
  pae_rd <- margin3 - margin2
  # Population attributable fraction (PAF)
  paf <- (margin3 - margin2) / margin3
  
  # Return all metrics as a named numeric vector
  return(c(margin0=margin0, margin1=margin1, margin2=margin2, margin3=margin3,
           afe_rr=afe_rr, afe_rd=afe_rd, afe=afe, 
           pae_rr=pae_rr, pae_rd=pae_rd, paf=paf))
}


# --------------------------------------------------
#   Loop across outcomes for primigravidae 
# --------------------------------------------------

# Restrict dataset to only those with propensity scores that overlapped
activeinf_primi <- db %>%
  filter(primi_actinf_pscore==1)

# Create empty list to store results 
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Generate point estimates using full dataset restricted to Primigravidae  -------
  full_est <- gravid_compute_metrics(activeinf_primi, y, 1:nrow(activeinf_primi), grav_subgroup = "Primigravida")
  
  # --- 2. Bootstrap for BCa 95% CIs --------
  boot_obj <- boot(
    data = activeinf_primi,
    statistic = function(d, i) gravid_compute_metrics(d, y, i, grav_subgroup = "Primigravida"),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_margin0 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_margin1 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  ci_margin2 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_margin3 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_afe_rr <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_afe_rd <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  ci_afe  <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  
  ci_pae_rr  <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pae_rd  <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_paf  <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  # --- 3. Combine point estimates and 95% CIs -------
  res <- c(
    
    margin0   = full_est["margin0"],
    margin0_lb = ci_margin0[1],
    margin0_ub = ci_margin0[2],
    
    margin1   = full_est["margin1"],
    margin1_lb = ci_margin1[1],
    margin1_ub = ci_margin1[2],
    
    margin2   = full_est["margin2"],
    margin2_lb = ci_margin2[1],
    margin2_ub = ci_margin2[2],
    
    margin3   = full_est["margin3"],
    margin3_lb = ci_margin3[1],
    margin3_ub = ci_margin3[2],
    
    afe_rr        = full_est["afe_rr"],
    afe_rr_lb     = ci_afe_rr[1],
    afe_rr_ub     = ci_afe_rr[2],
    
    afe_rd        = full_est["afe_rd"],
    afe_rd_lb     = ci_afe_rd[1],
    afe_rd_ub     = ci_afe_rd[2],
    
    afe        = full_est["afe"],
    afe_lb     = ci_afe[1],
    afe_ub     = ci_afe[2],
    
    pae_rr        = full_est["pae_rr"],
    pae_rr_lb     = ci_pae_rr[1],
    pae_rr_ub     = ci_pae_rr[2],
    
    pae_rd        = full_est["pae_rd"],
    pae_rd_lb     = ci_pae_rd[1],
    pae_rd_ub     = ci_pae_rd[2],
    
    paf        = full_est["paf"],
    paf_lb     = ci_paf[1],
    paf_ub     = ci_paf[2]
  ) 
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# -----------------------------------------------------------
#     Assemble primigravidae results into a data frame
# -----------------------------------------------------------
# Combine results into a dataframe
primi_df <- do.call(rbind, lapply(names(results_list), function(outcome) {
  vec <- results_list[[outcome]]        # the named vector
  df <- as.data.frame(t(unname(vec)))   # remove names, transpose, make 1-row df
  df$outcome <- outcome                 # add outcome column
  df })) %>%
  mutate(subgroup = 'Primigravidae')    # label that these are primigravide estimates

# Assign column names
colnames(primi_df) <- c(
  "margin0", "margin0_lb", "margin0_ub",
  "margin1", "margin1_lb", "margin1_ub",
  "margin2", "margin2_lb", "margin2_ub",
  "margin3", "margin3_lb", "margin3_ub",
  "afe_rr", "afe_rr_lb", "afe_rr_ub",
  "afe_rd", "afe_rd_lb", "afe_rd_ub",
  "afe", "afe_lb", "afe_ub",
  "pae_rr", "pae_rr_lb", "pae_rr_ub",
  "pae_rd", "pae_rd_lb", "pae_rd_ub",
  "paf", "paf_lb", "paf_ub",
  "outcome", "subgroup"
)

# Reorder outcome to be first
primi_df <- primi_df[, c("outcome", "subgroup", setdiff(colnames(primi_df), c("outcome", "subgroup")))]

rownames(primi_df) <- NULL
print(primi_df)


# --------------------------------------------------
#   Loop across outcomes for multigravidae 
# --------------------------------------------------

# Restrict dataset to only those with propensity scores that overlapped
activeinf_multi <- db %>%
  filter(multi_actinf_pscore==1)

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Full dataset model for point estimates ------
  full_est <- gravid_compute_metrics(activeinf_multi, y, 1:nrow(activeinf_multi), grav_subgroup = "Multigravida")
  
  # --- 2. Bootstrap for BCa CIs ------
  boot_obj <- boot(
    data = activeinf_multi,
    statistic = function(d, i) gravid_compute_metrics(d, y, i, grav_subgroup = "Multigravida"),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_margin0 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_margin1 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  ci_margin2 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_margin3 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_afe_rr <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_afe_rd <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  ci_afe  <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  
  ci_pae_rr  <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pae_rd  <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_paf  <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  # --- 3. Combine point estimates (from full data) + BCa CIs ------
  res <- c(
    
    margin0   = full_est["margin0"],
    margin0_lb = ci_margin0[1],
    margin0_ub = ci_margin0[2],
    
    margin1   = full_est["margin1"],
    margin1_lb = ci_margin1[1],
    margin1_ub = ci_margin1[2],
    
    margin2   = full_est["margin2"],
    margin2_lb = ci_margin2[1],
    margin2_ub = ci_margin2[2],
    
    margin3   = full_est["margin3"],
    margin3_lb = ci_margin3[1],
    margin3_ub = ci_margin3[2],
    
    afe_rr        = full_est["afe_rr"],
    afe_rr_lb     = ci_afe_rr[1],
    afe_rr_ub     = ci_afe_rr[2],
    
    afe_rd        = full_est["afe_rd"],
    afe_rd_lb     = ci_afe_rd[1],
    afe_rd_ub     = ci_afe_rd[2],
    
    afe        = full_est["afe"],
    afe_lb     = ci_afe[1],
    afe_ub     = ci_afe[2],
    
    pae_rr        = full_est["pae_rr"],
    pae_rr_lb     = ci_pae_rr[1],
    pae_rr_ub     = ci_pae_rr[2],
    
    pae_rd        = full_est["pae_rd"],
    pae_rd_lb     = ci_pae_rd[1],
    pae_rd_ub     = ci_pae_rd[2],
    
    paf        = full_est["paf"],
    paf_lb     = ci_paf[1],
    paf_ub     = ci_paf[2]
  ) 
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# ----------------------------------------------------------
#     Assemble mutigravidae results into a data frame
# ----------------------------------------------------------
multi_df <- do.call(rbind, lapply(names(results_list), function(outcome) {
  vec <- results_list[[outcome]]        # the named vector
  df <- as.data.frame(t(unname(vec)))   # remove names, transpose, make 1-row df
  df$outcome <- outcome                 # add outcome column
  df })) %>%
  mutate(subgroup = 'Multigravidae')

# Assign proper column names
colnames(multi_df) <- c(
  "margin0", "margin0_lb", "margin0_ub",
  "margin1", "margin1_lb", "margin1_ub",
  "margin2", "margin2_lb", "margin2_ub",
  "margin3", "margin3_lb", "margin3_ub",
  "afe_rr", "afe_rr_lb", "afe_rr_ub",
  "afe_rd", "afe_rd_lb", "afe_rd_ub",
  "afe", "afe_lb", "afe_ub",
  "pae_rr", "pae_rr_lb", "pae_rr_ub",
  "pae_rd", "pae_rd_lb", "pae_rd_ub",
  "paf", "paf_lb", "paf_ub",
  "outcome", "subgroup"
)

# Reorder outcome to be first
multi_df <- multi_df[, c("outcome", "subgroup", setdiff(colnames(multi_df), c("outcome", "subgroup")))]

rownames(multi_df) <- NULL
print(multi_df)

# ---------------------------------------
#   Combine datasets and export 
# ---------------------------------------
finalresults <- rbind(overall_df, primi_df, multi_df ) %>%
  mutate(runexp_amongexp_descrip = paste0(round(100*margin0, 1), "% [", round(100*margin0_lb, 1), "%, ", round(100*margin0_ub, 1), "%]"),
         rexp_amongexp_descrip = paste0(round(100*margin1, 1), "% [", round(100*margin1_lb, 1), "%, ", round(100*margin1_ub, 1), "%]"),
         runexp_pop_descrip = paste0(round(100*margin2, 1), "% [", round(100*margin2_lb, 1), "%, ", round(100*margin2_ub, 1), "%]"),
         rnc_pop_descrip = paste0(round(100*margin3, 1), "% [", round(100*margin3_lb, 1), "%, ", round(100*margin3_ub, 1), "%]"),
         
         afe_rr_descrip = paste0(round(afe_rr,2), " [", round(afe_rr_lb, 2), ", ", round(afe_rr_ub, 2), "]"), 
         afe_rd_descrip = paste0(round(100*afe_rd,1), "% [", round(100*afe_rd_lb, 1), "%, ", round(100*afe_rd_ub, 1), "%]"),
         afe_descrip = paste0(round(100*afe, 1), "% [", round(100*afe_lb, 1), "%, ", round(100*afe_ub, 1), "%]"),
         
         pae_rr_descrip = paste0(round(pae_rr,2), " [", round(pae_rr_lb, 2), ", ", round(pae_rr_ub, 2), "]"), 
         pae_rd_descrip = paste0(round(100*pae_rd,1), "% [", round(100*pae_rd_lb, 1), "%, ", round(100*pae_rd_ub, 1), "%]"), 
         paf_descrip = paste0(round(100*paf, 1), "% [", round(100*paf_lb, 1), "%, ", round(100*paf_ub, 1), "%]"),
         
         pigment_code = 4,
         
         outcome_code = case_when(
           outcome == "preterm"  ~ 1,
           outcome == "SGA"      ~ 2,
           outcome == "LBWdich"  ~ 3), 
         subgroup_code = case_when(
           subgroup == "Overall" ~ 1, 
           subgroup == "Primigravidae" ~ 2, 
           subgroup == "Multigravidae" ~ 3))


# -----------------------
#   Export results 
# -----------------------
write.xlsx(finalresults, file="Results/DPSP_Positivity_ActInf_Gcomp_Boot.xlsx")








