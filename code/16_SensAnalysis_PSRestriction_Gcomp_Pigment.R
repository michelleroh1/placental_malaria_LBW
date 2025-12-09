
##################################################################################
##	  Quantifying the contributions of active and past placental malaria infection
##    to low birthweight in a high transmission area of Uganda
##    16_Sensitivity Analyses Restricting Analyses of Severe Pigment to 
##       Individuals with Common Support
##
##    This script is the same as 7_Gcomp_ActiveInf, but restricts analysis to 
#     those with severe or no pigment and excludes those 
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
library(logistf)


# ----------------------------
#   Upload databases
# ----------------------------

#  Upload dataset with indicator variable of the individuals to exclude based on PS
pscore <- read.xlsx('Databases/Created Databases/Propensity Scores_Indicator for Restricted Analyses.xlsx') 

# Full dataset 
db <- read.dta13('Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta', nonint.factors = TRUE) %>%
  mutate(propfibrincat_sim_cat = case_when(
    propfibrincat_sim == 0 ~ "None", 
    propfibrincat_sim == 1 ~ "<10%", 
    propfibrincat_sim == 2 ~ "10 to <30%", 
    propfibrincat_sim == 3 ~ ">=30%")) %>%
  left_join(., pscore, by='id') %>%           # merge in pscore database
  mutate(severe = ifelse(propfibrincat==">=30%", 1, 0),
         severe_sim = ifelse(propfibrincat_sim_cat==">=30%", 1, 0)) %>%
  filter(propfibrincat_sim_cat==">=30%" | propfibrincat_sim_cat=="None") %>%
  mutate(preterm = case_when(
    preterm == "Yes" ~ 0,
    preterm == "No" ~ 1), 
    LBWdich = case_when(
      LBWdich == "No" ~ 0, 
      LBWdich == "Yes" ~ 1)) 

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
#   - If logistic regression does not converge, try Firth's exact logistic regression
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
  
  # Define a safe wrapper to return NA if anything fails
  safe_return <- function() {
    setNames(rep(NA, 10),
             c("margin0", "margin1", "margin2", "margin3",
               "afe_rr", "afe_rd", "afe",
               "pae_rr", "pae_rd", "paf"))
  }
  
  # -----------------------
  # 1. Fit outcome model
  # -----------------------
  out <- tryCatch({
    # Full model with interaction between severe infection infection and gravidity
    formula_full <- as.formula(
      paste(outcome_var, "~ severe*graviddich + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
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
          "~ severe + graviddich + treatmentarm + age_enroll + age2",
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
      filter(severe_sim==1) %>%
      mutate(severe=0,              # set exposure to 0 (no severe pigment) for model prediction
             treatmentarm="SP")       # set treatment arm to SP for all
    data1 <- d %>%
      filter(severe_sim==1) %>%
      mutate(severe=1,              # set exposure to 1 (severe pigment) for model prediction
             treatmentarm="SP")       # set treatment arm to SP for all
    
    # ---- PAE: full population under natural course exposure assignments ----
    data2 <- d %>%
      mutate(severe = 0,            # set exposure to 0 (no severe pigment) for model prediction 
             treatmentarm = "SP")
    data3 <- d %>%
      mutate(severe = severe_sim, # set exposure to simulated natural course exposure distribution
             treatmentarm = "SP")
    
    # -------------------------------------------------
    # 3. Compute marginal risks under each scenario
    # -------------------------------------------------
    # Predict outcome probabilities for each participant, then estimate risk
    margin0 <- mean(predict(model, newdata=data0, type="response"))
    margin1 <- mean(predict(model, newdata=data1, type="response"))
    margin2 <- mean(predict(model, newdata=data2, type="response"))
    margin3 <- mean(predict(model, newdata=data3, type="response"))
    
    # Check for NA/NaN predictions
    if (any(!is.finite(c(margin0, margin1, margin2, margin3)))) return(safe_return())
    
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
    
    return(c(margin0 = margin0, margin1 = margin1, margin2 = margin2, margin3 = margin3,
             afe_rr = afe_rr, afe_rd = afe_rd, afe = afe,
             pae_rr = pae_rr, pae_rd = pae_rd, paf = paf))
  
  }, error = function(e) safe_return())
  
  return(out)
}

# -----------------------------------------------------------
#   Loop g-computation function over all birth outcomes
# -----------------------------------------------------------
# For each outcome:
#   - Compute point estimates using the full dataset
#   - Bootstrap the metrics to obtain BCa 95% CIs
#   - Store estimates + CIs in results vector

# Restrict dataset to only those with propensity scores that overlapped
severe_overall <- db %>%
  filter(overall_severe_pscore==1)

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Generate point estimates using full dataset  -------
  full_est <- overall_compute_metrics(severe_overall, y, 1:nrow(severe_overall))
  
  # --- 2. Bootstrap for BCa 95% CIs --------
  boot_obj <- boot(
    data = db,
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
gravid_compute_metrics <- function(data, outcome_var, indices, grav_subgroup = NULL) {
  # Resample rows according to indices (needed for when bootstrapping)
  d <- data[indices, ]
  
  # Filter by gravidity subgroup 
  if (!is.null(grav_subgroup)) {
    d <- d[d$graviddich == grav_subgroup, ]
  }
  
  # --- Guard clauses to avoid crashes ---
  if (nrow(d) < 10) {
    return(rep(NA_real_, 10))  # not enough data to fit a model
  }
  
  # --- If no variation in exposure or outcome
  if (length(unique(d$severe)) < 2 || length(unique(d[[outcome_var]])) < 2) {
    return(rep(NA_real_, 10))  
  }
  
  # -----------------------
  # 1. Fit outcome model
  # -----------------------
  formula_full <- as.formula(
    paste(
      outcome_var,
      "~ severe + treatmentarm + age_enroll + age2 + wealthcat_hd +",
      "bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male"
    )
  )
  
  # If glm failed or clearly did not converge, try Firth logistic regression
  model <- tryCatch(
    glm(formula_full, data=d, family=binomial(link="logit")),
    error = function(e) {
      message("Switching to Firth logistic regression due to separation")
      logistf(formula_full, data=d)
    }
  )
  
  # If Firth also failed, give up and return all NAs
  if (is.null(model)) {
    return(rep(NA_real_, 10))
  }
  
  # ------------------------------------------------------------
  # 2. Define counterfactual scenarios for g-computation
  # ------------------------------------------------------------
  # ---- AFE: subset of women with severe_sim == 1 ----
  data0 <- d %>%
    filter(severe_sim==1) %>%
    mutate(severe=0,               # set exposure to 0 (no severe pigment) for model prediction
           treatmentarm="SP")        # set treatment arm to SP for all
  # All exposed among exposed
  data1 <- d %>%
    filter(severe_sim==1) %>%
    mutate(severe=1,               # set exposure to 1 (severe pigment) for model prediction
           treatmentarm="SP")        # set treatment arm to SP for all
  
  # ---- PAE: full population under natural course exposure assignments ----
  data2 <- d %>%
    mutate(severe = 0,             # set exposure to 0 (severe pigment) for model prediction 
           treatmentarm = "SP")
  data3 <- d %>%
    mutate(severe = severe_sim,  # set exposure to simulated natural course exposure distribution
           treatmentarm = "SP")
  
  # If any of the counterfactual datasets are empty, bail out
  if (nrow(data0) == 0L || nrow(data1) == 0L || nrow(data2) == 0L || nrow(data3) == 0L) {
    return(rep(NA_real_, 10))
  }
  
  # -------------------------------------------------
  # 3. Compute marginal risks under each scenario
  # -------------------------------------------------
  safe_mean_pred <- function(newdata) {
    tryCatch(
      mean(predict(model, newdata = newdata, type = "response")),
      error = function(e) NA_real_
    )
  }
  
  # Predict outcome probabilities for each participant, then estimate risk
  margin0 <- safe_mean_pred(data0) # none exposed among exposed
  margin1 <- safe_mean_pred(data1) # all exposed among exposed
  margin2 <- safe_mean_pred(data2) # none exposed in pop
  margin3 <- safe_mean_pred(data3) # NC exposure
  
  # If any margins missing, just keep as NA
  if (any(!is.finite(c(margin0, margin1, margin2, margin3)))) {
    return(rep(NA_real_, 10))
  }
  
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
  
  # If any return as infinite, then keep all as NA
  if (any(!is.finite(c(afe_rr, afe_rd, afe, pae_rr, pae_rd, paf)))) {
    return(rep(NA_real_, 10))
  }
  
  return(c(
    margin0 = margin0, margin1 = margin1, margin2 = margin2, margin3 = margin3,
    afe_rr = afe_rr, afe_rd = afe_rd, afe = afe,
    pae_rr = pae_rr, pae_rd = pae_rd, paf = paf
  ))
}
  
# --------------------------------------------------
#   Loop outcomes for primigravidae
# --------------------------------------------------

# Restrict dataset to only those with propensity scores that overlapped
severe_primi <- db %>%
  filter(primi_severe_pscore==1)

# Create empty list to store results 
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Generate point estimates using full dataset restricted to Primigravidae  -------
  full_est <- gravid_compute_metrics(severe_primi, y, 1:nrow(severe_primi), grav_subgroup = "Primigravida")
  
  # --- 2. Bootstrap for BCa 95% CIs --------
  boot_obj <- boot(
    data = severe_primi,
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
severe_multi <- db %>%
  filter(multi_severe_pscore==1)

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Full dataset model for point estimates ------
  full_est <- gravid_compute_metrics(severe_multi, y, 1:nrow(severe_multi), grav_subgroup = "Multigravida")
  
  # --- 2. Bootstrap for BCa CIs ------
  boot_obj <- boot(
    data = severe_multi,
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
write.xlsx(finalresults, file="Results/DPSP_Positivity_SeverePigment_Gcomp_Boot.xlsx")


