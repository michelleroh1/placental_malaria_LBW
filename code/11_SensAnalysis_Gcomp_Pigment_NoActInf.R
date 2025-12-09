
##################################################################################
##	  Quantifying the contributions of active and past placental malaria infection
##    to low birthweight in a high transmission area of Uganda
##	  11_Sensitivity Analyses - Repeat g-computation analyses of pigment, 
##                             but exclude women with both active + severe pigment
##
##  This script is the same as 8_Gcomputation_Pigment, but exclude those 
##   both active infection + severe pigmentation 
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
library(purrr)

# ----------------------------
#   Upload databases
# ----------------------------
db <- read.dta13("Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta",
                 nonint.factors = TRUE) %>%
  filter(!(propfibrincat==">=30%" & activeHP==1)) %>%
  mutate(propfibrincat_sim_cat = case_when(
    propfibrincat_sim == 0 ~ "None", 
    propfibrincat_sim == 1 ~ "<10%", 
    propfibrincat_sim == 2 ~ "10 to <30%", 
    propfibrincat_sim == 3 ~ ">=30%"))


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
#  For a given outcome:
#    1. Fits a logistic model with propfibrincat 
#         - Includes propfibrincat*graviddich interaction term if model 
#           model converges, otherwise use main effects model only 
#    2. Constructs counterfactual datasets:
#         - data0 : everyone set to "None" pigment (no pigment)
#         - data1 : pigment set to simulated "natural course" under IPTp-SP
#         - data01: like data1, but mild pigment (<10%) removed
#         - data02: like data1, but moderate pigment (10–<30%) removed
#         - data03: like data1, but severe pigment (>=30%) removed
#    3. Computes:
#         - AFE-type quantities among women who would naturally
#           have mild/moderate/severe pigment under IPTp-SP
#           (exp_margin_0k, exp_margin_kk, RR, RD, AF)
#         - PAF-type quantities at the population level for
#           eliminating each pigment category in turn

overall_compute_metrics <- function(data, outcome_var, indices) {
  # Resample rows according to indices (needed for when bootstrapping)
  d <- data[indices, ]
  
  # Full model with interaction between pigment category and gravidity
  formula_full <- as.formula(
    paste(outcome_var, "~ propfibrincat*graviddich + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
  )
  
  model <- try(glm(formula_full, data=d, family=binomial(link="logit")), silent=TRUE)
  
  # If model with interaction term fails, remove interaction term
  if (inherits(model, "try-error") || all(abs(coef(model)[grep("propfibrincat:graviddich", names(coef(model)))] < 1e-12))) {
    formula_maineffects <- as.formula(
      paste(outcome_var, "~ propfibrincat + graviddich + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
    )
    
    model <- glm(formula_maineffects, data=d, family=binomial(link="logit"))
  }
  
  # ------------------------------------------------------------
  # 2. Define counterfactual scenarios for g-computation
  # ------------------------------------------------------------
  # Key variables:
  #   - propfibrincat: Observed pigment category at delivery:
  #                    "None", "<10%", "10 to <30%", ">=30%"
  #
  #   - propfibrincat_sim: Simulated "natural course" exposure 
  #                        (used as the counterfactual scenario had all women 
  #                        received standard-of-care IPTp-SP) 
  #                1 = mild    (<10%)
  #                2 = moderate (10 to <30%)
  #                3 = severe  (>=30%)
  #
  # We construct:
  #   data0 : All women were are assigned to "None" pigment category 
  #           (used to compute risks if pigment were eliminated)
  #
  #   data1 : Everyone assigned to their natural exposure (propfibrincat_sim_cat)
  #
  #   data01, data02, data03 :  Start from data1, then remove the pigment category
  #                             of interest to "None" (hypothetical intervention
  #                             to eliminate the specific pigment category)
  #           data01: Change only mild (<10%) pigment category to "None"
  #           data02: Change only moderate (10 to <30%) category to "None"
  #           data03: Change only severe (>=30%) category to "None"
  # 
  # From these, we define:
  #
  #   Among exposed (AFE-type), for each k = 1,2,3:
  #     exp_margin_0k: mean predicted risk if pigment set to "None"
  #                    among women with propfibrincat_sim == k
  #     exp_margin_kk: mean predicted risk under natural-course pigment k
  #
  #   Population (PAF-type), for each k = 1,2,3:
  #     pop_margin_0: mean risk in data1 (natural-course pigment distribution)
  #     pop_margin_k: mean risk in data0k (data01, data02, data03),
  #                   where category k is eliminated.
  #
  # ------------------------------------------------------------
  
  # Create datasets where treatmentarm = SP and exposure is set at varying levels
  # 1) All set to no pigment
  data0 <- d %>%
    mutate(propfibrincat = "None", 
           treatmentarm = "SP")
  # 2) Natural-course pigment under IPTp-SP
  data1 <- d %>%
    mutate(propfibrincat = propfibrincat_sim_cat, 
           treatmentarm = "SP")
  
  # 3) Datasets removing mild (01), moderate (02), or severe (03) pigment
  data01 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat=="<10%", "None", propfibrincat))
  data02 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat=="10 to <30%", "None", propfibrincat))
  data03 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat==">=30%", "None", propfibrincat))
  
  # -------------------------------------------------
  # 3. Compute marginal risks under each scenario
  # -------------------------------------------------
  # Predict outcome probabilities for each participant
  data0$propfibrincat_cf <- predict(model, newdata=data0, type="response")
  data1$propfibrincat_cf <- predict(model, newdata=data1, type="response")
  
  # # --- AFE-type: marginal risks among subset of women in each pigment category ---
  # Mild
  exp_margin_01 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==1])
  exp_margin_11 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==1])
  
  # Moderate
  exp_margin_02 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==2])
  exp_margin_22 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==2])
  
  # Severe
  exp_margin_03 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==3])
  exp_margin_33 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==3])
  
  # --- PAF-type: population risks under elimination of each pigment category ---
  # Predict outcomes in CF databases removing mild, moderate, or severe pigment
  data01$propfibrincat_cf <- predict(model, newdata=data01, type="response")
  data02$propfibrincat_cf <- predict(model, newdata=data02, type="response")
  data03$propfibrincat_cf <- predict(model, newdata=data03, type="response") 
  
  # Population mean risks
  pop_margin_0 <- mean(data1$propfibrincat_cf) # Natural course exposure distribution
  pop_margin_1 <- mean(data01$propfibrincat_cf) # Remove Mild only 
  pop_margin_2 <- mean(data02$propfibrincat_cf) # Remove Moderate only 
  pop_margin_3 <- mean(data03$propfibrincat_cf) # Remove Severe only 
  
  # -------------------------------------------------
  # 4. Attributable effects among the exposed (AFE)
  # -------------------------------------------------
  # Risk ratios among exposed in category k
  afe_rr_1 = exp_margin_11 / exp_margin_01
  afe_rr_2 = exp_margin_22 / exp_margin_02
  afe_rr_3 = exp_margin_33 / exp_margin_03
  
  # Risk differences among exposed in category k
  afe_rd_1 = exp_margin_11 - exp_margin_01
  afe_rd_2 = exp_margin_22 - exp_margin_02
  afe_rd_3 = exp_margin_33 - exp_margin_03
  
  # Attributable fraction among the exposed in category k
  afe_1  = afe_rd_1 / exp_margin_11
  afe_2  = afe_rd_2 / exp_margin_22
  afe_3  = afe_rd_3 / exp_margin_33
  
  # -------------------------------------------------
  # 5. Population attributable effects (PAF-type)
  # -------------------------------------------------
  # Population risk ratios when eliminating pigment category k
  pae_rr_1 = pop_margin_0 / pop_margin_1
  pae_rr_2 = pop_margin_0 / pop_margin_2
  pae_rr_3 = pop_margin_0 / pop_margin_3
  
  # Population risk differences when eliminating pigment category k
  pae_rd_1 = pop_margin_0 - pop_margin_1
  pae_rd_2 = pop_margin_0 - pop_margin_2
  pae_rd_3 = pop_margin_0 - pop_margin_3
  
  # Population attributable fractions for category k
  paf_1  = pae_rd_1 / pop_margin_0
  paf_2  = pae_rd_2 / pop_margin_0
  paf_3  = pae_rd_3 / pop_margin_0
  
  # Return all metrics as a named vector
  return(c(exp_margin_01=exp_margin_01, exp_margin_11=exp_margin_11,
           exp_margin_02=exp_margin_02, exp_margin_22=exp_margin_22,
           exp_margin_03=exp_margin_03, exp_margin_33=exp_margin_33,
           pop_margin_0=pop_margin_0, pop_margin_1=pop_margin_1, 
           pop_margin_2=pop_margin_2, pop_margin_3=pop_margin_3, 
           afe_rr_1=afe_rr_1, afe_rr_2=afe_rr_2, afe_rr_3=afe_rr_3,
           afe_rd_1=afe_rd_1, afe_rd_2=afe_rd_2, afe_rd_3=afe_rd_3, 
           afe_1=afe_1, afe_2=afe_2, afe_3=afe_3, 
           pae_rr_1=pae_rr_1, pae_rr_2=pae_rr_2, pae_rr_3=pae_rr_3,
           pae_rd_1=pae_rd_1, pae_rd_2=pae_rd_2, pae_rd_3=pae_rd_3,
           paf_1=paf_1, paf_2=paf_2, paf_3=paf_3))
}

# ------------------------------------------------------------------
#   Loop over outcomes: overall estimates + bootstrap CIs
# ------------------------------------------------------------------

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Full dataset model for point estimates -----
  full_est <- overall_compute_metrics(db, y, 1:nrow(db))
  
  # --- 2. Bootstrap for BCa CIs -----
  boot_obj <- boot(
    data = db,
    statistic = function(d, i) overall_compute_metrics(d, y, i),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_exp_margin_01 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_exp_margin_11 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  
  ci_exp_margin_02 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_exp_margin_22 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_exp_margin_03 <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_exp_margin_33 <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  
  ci_pop_margin_0 <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  ci_pop_margin_1 <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pop_margin_2 <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_pop_margin_3 <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  ci_afe_rr_1 <- boot.ci(boot_obj, type="bca", index=11)$bca[4:5]
  ci_afe_rr_2 <- boot.ci(boot_obj, type="bca", index=12)$bca[4:5]
  ci_afe_rr_3 <- boot.ci(boot_obj, type="bca", index=13)$bca[4:5]
  
  ci_afe_rd_1 <- boot.ci(boot_obj, type="bca", index=14)$bca[4:5]
  ci_afe_rd_2 <- boot.ci(boot_obj, type="bca", index=15)$bca[4:5]
  ci_afe_rd_3 <- boot.ci(boot_obj, type="bca", index=16)$bca[4:5]
  
  ci_afe_1 <- boot.ci(boot_obj, type="bca", index=17)$bca[4:5]
  ci_afe_2 <- boot.ci(boot_obj, type="bca", index=18)$bca[4:5]
  ci_afe_3 <- boot.ci(boot_obj, type="bca", index=19)$bca[4:5]
  
  ci_pae_rr_1 <- boot.ci(boot_obj, type="bca", index=20)$bca[4:5]
  ci_pae_rr_2 <- boot.ci(boot_obj, type="bca", index=21)$bca[4:5]
  ci_pae_rr_3 <- boot.ci(boot_obj, type="bca", index=22)$bca[4:5]
  
  ci_pae_rd_1 <- boot.ci(boot_obj, type="bca", index=23)$bca[4:5]
  ci_pae_rd_2 <- boot.ci(boot_obj, type="bca", index=24)$bca[4:5]
  ci_pae_rd_3 <- boot.ci(boot_obj, type="bca", index=25)$bca[4:5]
  
  ci_paf_1 <- boot.ci(boot_obj, type="bca", index=26)$bca[4:5]
  ci_paf_2 <- boot.ci(boot_obj, type="bca", index=27)$bca[4:5]
  ci_paf_3 <- boot.ci(boot_obj, type="bca", index=28)$bca[4:5]
  
  # --- 3. Combine point estimates + 95% CIs
  res <- c(
    
    exp_margin_01    = full_est["exp_margin_01"],
    exp_margin_01_lb   = ci_exp_margin_01[1],
    exp_margin_01_ub   = ci_exp_margin_01[2],
    
    exp_margin_11    = full_est["exp_margin_11"],
    exp_margin_11_lb   = ci_exp_margin_11[1],
    exp_margin_11_ub   = ci_exp_margin_11[2],
    
    exp_margin_02    = full_est["exp_margin_02"],
    exp_margin_02_lb   = ci_exp_margin_02[1],
    exp_margin_02_ub   = ci_exp_margin_02[2],
    
    exp_margin_22    = full_est["exp_margin_22"],
    exp_margin_22_lb   = ci_exp_margin_22[1],
    exp_margin_22_ub   = ci_exp_margin_22[2],
    
    exp_margin_03    = full_est["exp_margin_03"],
    exp_margin_03_lb   = ci_exp_margin_03[1],
    exp_margin_03_ub   = ci_exp_margin_03[2],
    
    exp_margin_33    = full_est["exp_margin_33"],
    exp_margin_33_lb   = ci_exp_margin_33[1],
    exp_margin_33_ub   = ci_exp_margin_33[2],
    
    pop_margin_0    = full_est["pop_margin_0"],
    pop_margin_0_lb   = ci_pop_margin_0[1],
    pop_margin_0_ub   = ci_pop_margin_0[2],
    
    pop_margin_1    = full_est["pop_margin_1"],
    pop_margin_1_lb   = ci_pop_margin_1[1],
    pop_margin_1_ub   = ci_pop_margin_1[2],
    
    pop_margin_2    = full_est["pop_margin_2"],
    pop_margin_2_lb   = ci_pop_margin_2[1],
    pop_margin_2_ub   = ci_pop_margin_2[2],
    
    pop_margin_3    = full_est["pop_margin_3"],
    pop_margin_3_lb   = ci_pop_margin_3[1],
    pop_margin_3_ub   = ci_pop_margin_3[2],
    
    afe_rr_1      = full_est["afe_rr_1"],
    afe_rr_1_lb   = ci_afe_rr_1[1],
    afe_rr_1_ub   = ci_afe_rr_1[2],
    
    afe_rr_2      = full_est["afe_rr_2"],
    afe_rr_2_lb   = ci_afe_rr_2[1],
    afe_rr_2_ub   = ci_afe_rr_2[2],
    
    afe_rr_3      = full_est["afe_rr_3"],
    afe_rr_3_lb   = ci_afe_rr_3[1],
    afe_rr_3_ub   = ci_afe_rr_3[2],
    
    afe_rd_1      = full_est["afe_rd_1"],
    afe_rd_1_lb   = ci_afe_rd_1[1],
    afe_rd_1_ub   = ci_afe_rd_1[2],
    
    afe_rd_2      = full_est["afe_rd_2"],
    afe_rd_2_lb   = ci_afe_rd_2[1],
    afe_rd_2_ub   = ci_afe_rd_2[2],
    
    afe_rd_3      = full_est["afe_rd_3"],
    afe_rd_3_lb   = ci_afe_rd_3[1],
    afe_rd_3_ub   = ci_afe_rd_3[2],
    
    afe_1      = full_est["afe_1"],
    afe_1_lb   = ci_afe_1[1],
    afe_1_ub   = ci_afe_1[2],
    
    afe_2      = full_est["afe_2"],
    afe_2_lb   = ci_afe_2[1],
    afe_2_ub   = ci_afe_2[2],
    
    afe_3      = full_est["afe_3"],
    afe_3_lb   = ci_afe_3[1],
    afe_3_ub   = ci_afe_3[2],
    
    pae_rr_1      = full_est["pae_rr_1"],
    pae_rr_1_lb   = ci_pae_rr_1[1],
    pae_rr_1_ub   = ci_pae_rr_1[2],
    
    pae_rr_2      = full_est["pae_rr_2"],
    pae_rr_2_lb   = ci_pae_rr_2[1],
    pae_rr_2_ub   = ci_pae_rr_2[2],
    
    pae_rr_3      = full_est["pae_rr_3"],
    pae_rr_3_lb   = ci_pae_rr_3[1],
    pae_rr_3_ub   = ci_pae_rr_3[2],
    
    pae_rd_1      = full_est["pae_rd_1"],
    pae_rd_1_lb   = ci_pae_rd_1[1],
    pae_rd_1_ub   = ci_pae_rd_1[2],
    
    pae_rd_2      = full_est["pae_rd_2"],
    pae_rd_2_lb   = ci_pae_rd_2[1],
    pae_rd_2_ub   = ci_pae_rd_2[2],
    
    pae_rd_3      = full_est["pae_rd_3"],
    pae_rd_3_lb   = ci_pae_rd_3[1],
    pae_rd_3_ub   = ci_pae_rd_3[2],
    
    paf_1      = full_est["paf_1"],
    paf_1_lb   = ci_paf_1[1],
    paf_1_ub   = ci_paf_1[2],
    
    paf_2      = full_est["paf_2"],
    paf_2_lb   = ci_paf_2[1],
    paf_2_ub   = ci_paf_2[2],
    
    paf_3      = full_est["paf_3"],
    paf_3_lb   = ci_paf_3[1],
    paf_3_ub   = ci_paf_3[2]
    
  )
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# ---------------------------------------------
#     Assemble results into a data frame
# ---------------------------------------------
# Combine results into a dataframe
overall_df <- imap_dfr(results_list, ~{
  as_tibble_row(.x) %>%
    mutate(outcome = .y, 
           subgroup = "Overall")
})

# Remove repeated column names
colnames(overall_df) <- sub("^([^.]+)\\.\\1$", "\\1", colnames(overall_df)) 

# Reorder so outcome and subgroup appear first
overall_df <- overall_df[, c("outcome", "subgroup", setdiff(colnames(overall_df), c("outcome", "subgroup")))]

# Create long database (a little too difficult to reshape all pigment categories; sort out each pigment category and rbind) -----
mild <- overall_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_01, exp_margin_01_lb, exp_margin_01_ub, 
                exp_margin_11,  exp_margin_11_lb,  exp_margin_11_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_1, pop_margin_1_lb, pop_margin_1_ub, 
                afe_rr_1, afe_rr_1_lb, afe_rr_1_ub, 
                afe_rd_1, afe_rd_1_lb, afe_rd_1_ub, 
                afe_1, afe_1_lb, afe_1_ub, 
                pae_rr_1, pae_rr_1_lb, pae_rr_1_ub, 
                pae_rd_1, pae_rd_1_lb, pae_rd_1_ub, 
                paf_1, paf_1_lb, paf_1_ub
  ) %>%
  rename_with(~ gsub("_01", "_0", .x)) %>%                # _01 → _0
  rename_with(~ gsub("_11", "", .x)) %>%                  # _11 → removed
  rename_with(~ gsub("_1(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _1 only if standalone
  mutate(pigment_code = 1)

moderate <- overall_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_02, exp_margin_02_lb, exp_margin_02_ub, 
                exp_margin_22,  exp_margin_22_lb,  exp_margin_22_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_2, pop_margin_2_lb, pop_margin_2_ub, 
                afe_rr_2, afe_rr_2_lb, afe_rr_2_ub, 
                afe_rd_2, afe_rd_2_lb, afe_rd_2_ub, 
                afe_2, afe_2_lb, afe_2_ub, 
                pae_rr_2, pae_rr_2_lb, pae_rr_2_ub, 
                pae_rd_2, pae_rd_2_lb, pae_rd_2_ub, 
                paf_2, paf_2_lb, paf_2_ub
  ) %>%
  rename_with(~ gsub("_02", "_0", .x)) %>%                # _02 → _0
  rename_with(~ gsub("_22", "", .x)) %>%                  # _22 → removed
  rename_with(~ gsub("_2(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _2 only if standalone
  mutate(pigment_code = 2)

severe <- overall_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_03, exp_margin_03_lb, exp_margin_03_ub, 
                exp_margin_33,  exp_margin_33_lb,  exp_margin_33_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_3, pop_margin_3_lb, pop_margin_3_ub, 
                afe_rr_3, afe_rr_3_lb, afe_rr_3_ub, 
                afe_rd_3, afe_rd_3_lb, afe_rd_3_ub, 
                afe_3, afe_3_lb, afe_3_ub, 
                pae_rr_3, pae_rr_3_lb, pae_rr_3_ub, 
                pae_rd_3, pae_rd_3_lb, pae_rd_3_ub, 
                paf_3,  paf_3_lb,  paf_3_ub
  ) %>%
  rename_with(~ gsub("_03", "_0", .x)) %>%                # _03 → _0
  rename_with(~ gsub("_33", "", .x)) %>%                  # _33 → removed
  rename_with(~ gsub("_3(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _3 only if standalone
  mutate(pigment_code = 3)

overall_df_final <- rbind(mild, moderate, severe)
print(overall_df_final)

# ==================================================================
#   PART 2: GRAVIDITY-SPECIFIC G-COMPUTATION
# ==================================================================

# ------------------------------------------------------------------
#  gravid_compute_metrics()
#
#  Same logic as overall_compute_metrics(), but applied within
#  given gravidity subgroups: "Primigravida" or "Multigravida"
#
#  Returns the same set of metrics (exp_margin_0k, exp_margin_kk,
#  pop_margin_k, afe_rr_k, afe_rd_k, afe_k, pae_rr_k, pae_rd_k, paf_k).
# ------------------------------------------------------------------

gravid_compute_metrics <- function(data, outcome_var, indices, grav_subgroup = NULL) {
  
  # Resample rows according to indices (needed for when bootstrapping)
  d <- data[indices, ]
  
  # Filter by gravidity subgroup
  if (!is.null(grav_subgroup)) {
    d <- d[d$graviddich == grav_subgroup, ]
  }
  
  # -------------------------------------------------------
  # 1. Fit outcome model (no interaction term here)
  # -------------------------------------------------------
  formula_full <- as.formula(
    paste(outcome_var, "~ propfibrincat + treatmentarm + age_enroll + age2 + wealthcat_hd + bmi_enroll_hd + bmi_enroll_hd2 + GAenroll + gaenroll2 + educ_cat + male")
  )
  
  model <- try(glm(formula_full, data=d, family=binomial(link="logit")), silent=TRUE)
  
  # ------------------------------------------------------------
  # 2. Define counterfactual scenarios for g-computation
  # ------------------------------------------------------------
  # 1) All set to no pigment
  data0 <- d %>%
    mutate(propfibrincat = "None", 
           treatmentarm = "SP")
  
  # 2) Natural-course pigment under IPTp-SP
  data1 <- d %>%
    mutate(propfibrincat = propfibrincat_sim_cat, 
           treatmentarm = "SP")
  
  # 3) Datasets removing mild (01), moderate (02), or severe (03) pigment
  data01 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat=="<10%", "None", propfibrincat))
  data02 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat=="10 to <30%", "None", propfibrincat))
  data03 <- data1 %>%
    mutate(propfibrincat = ifelse(propfibrincat==">=30%", "None", propfibrincat))
  
  # -------------------------------------------------
  # 3. Compute marginal risks under each scenario
  # -------------------------------------------------
  # Predict outcome probabilities for each participant
  data0$propfibrincat_cf <- predict(model, newdata=data0, type="response")
  data1$propfibrincat_cf <- predict(model, newdata=data1, type="response")
  
  # --- AFE-type: marginal risks among subset of women in each pigment category ---
  # Mild
  exp_margin_01 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==1])
  exp_margin_11 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==1])
  
  # Moderate
  exp_margin_02 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==2])
  exp_margin_22 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==2])
  
  # Severe
  exp_margin_03 <- mean(data0$propfibrincat_cf[data0$propfibrincat_sim==3])
  exp_margin_33 <- mean(data1$propfibrincat_cf[data1$propfibrincat_sim==3])
  
  # --- PAF-type: population risks under elimination of each pigment category ---
  # Predict outcomes in CF databases removing mild, moderate, or severe pigment
  data01$propfibrincat_cf <- predict(model, newdata=data01, type="response")
  data02$propfibrincat_cf <- predict(model, newdata=data02, type="response")
  data03$propfibrincat_cf <- predict(model, newdata=data03, type="response") 
  
  # Population mean risks 
  pop_margin_0 <- mean(data1$propfibrincat_cf) # NC exposure
  pop_margin_1 <- mean(data01$propfibrincat_cf) # Remove Mild only 
  pop_margin_2 <- mean(data02$propfibrincat_cf) # Remove Moderate only 
  pop_margin_3 <- mean(data03$propfibrincat_cf) # Remove Severe only 
  
  # -------------------------------------------------
  # 4. Attributable effects among the exposed (AFE)
  # -------------------------------------------------
  # Risk ratios among exposed in category k
  afe_rr_1 = exp_margin_11 / exp_margin_01
  afe_rr_2 = exp_margin_22 / exp_margin_02
  afe_rr_3 = exp_margin_33 / exp_margin_03
  
  # Risk differences among exposed in category k
  afe_rd_1 = exp_margin_11 - exp_margin_01
  afe_rd_2 = exp_margin_22 - exp_margin_02
  afe_rd_3 = exp_margin_33 - exp_margin_03
  
  # Attributable fraction among the exposed in category k
  afe_1  = afe_rd_1 / exp_margin_11
  afe_2  = afe_rd_2 / exp_margin_22
  afe_3  = afe_rd_3 / exp_margin_33
  
  # -------------------------------------------------
  # 5. Population attributable effects (PAF-type)
  # -------------------------------------------------
  # Population risk ratios when eliminating pigment category k
  pae_rr_1 = pop_margin_0 / pop_margin_1
  pae_rr_2 = pop_margin_0 / pop_margin_2
  pae_rr_3 = pop_margin_0 / pop_margin_3
  
  # Population risk differences when eliminating pigment category k
  pae_rd_1 = pop_margin_0 - pop_margin_1
  pae_rd_2 = pop_margin_0 - pop_margin_2
  pae_rd_3 = pop_margin_0 - pop_margin_3
  
  # Population attributable fractions for category k
  paf_1  = pae_rd_1 / pop_margin_0
  paf_2  = pae_rd_2 / pop_margin_0
  paf_3  = pae_rd_3 / pop_margin_0
  
  # Return all metrics as a named vector
  return(c(exp_margin_01=exp_margin_01, exp_margin_11=exp_margin_11,
           exp_margin_02=exp_margin_02, exp_margin_22=exp_margin_22,
           exp_margin_03=exp_margin_03, exp_margin_33=exp_margin_33,
           pop_margin_0=pop_margin_0, pop_margin_1=pop_margin_1, 
           pop_margin_2=pop_margin_2, pop_margin_3=pop_margin_3, 
           afe_rr_1=afe_rr_1, afe_rr_2=afe_rr_2, afe_rr_3=afe_rr_3,
           afe_rd_1=afe_rd_1, afe_rd_2=afe_rd_2, afe_rd_3=afe_rd_3, 
           afe_1=afe_1, afe_2=afe_2, afe_3=afe_3, 
           pae_rr_1=pae_rr_1, pae_rr_2=pae_rr_2, pae_rr_3=pae_rr_3,
           pae_rd_1=pae_rd_1, pae_rd_2=pae_rd_2, pae_rd_3=pae_rd_3,
           paf_1=paf_1, paf_2=paf_2, paf_3=paf_3))
}


# ------------------------------------------------------------------
#   Primigravidae estimates: Loop over outcomes + bootstrap CIs
# ------------------------------------------------------------------

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Full dataset model for point estimates ----
  full_est <- gravid_compute_metrics(db, y, 1:nrow(db), grav_subgroup = "Primigravida")
  
  # --- 2. Bootstrap for BCa CIs ----
  boot_obj <- boot(
    data = db,
    statistic = function(d, i) gravid_compute_metrics(d, y, i, grav_subgroup = "Primigravida"),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_exp_margin_01 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_exp_margin_11 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  
  ci_exp_margin_02 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_exp_margin_22 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_exp_margin_03 <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_exp_margin_33 <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  
  ci_pop_margin_0 <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  ci_pop_margin_1 <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pop_margin_2 <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_pop_margin_3 <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  ci_afe_rr_1 <- boot.ci(boot_obj, type="bca", index=11)$bca[4:5]
  ci_afe_rr_2 <- boot.ci(boot_obj, type="bca", index=12)$bca[4:5]
  ci_afe_rr_3 <- boot.ci(boot_obj, type="bca", index=13)$bca[4:5]
  
  ci_afe_rd_1 <- boot.ci(boot_obj, type="bca", index=14)$bca[4:5]
  ci_afe_rd_2 <- boot.ci(boot_obj, type="bca", index=15)$bca[4:5]
  ci_afe_rd_3 <- boot.ci(boot_obj, type="bca", index=16)$bca[4:5]
  
  ci_afe_1 <- boot.ci(boot_obj, type="bca", index=17)$bca[4:5]
  ci_afe_2 <- boot.ci(boot_obj, type="bca", index=18)$bca[4:5]
  ci_afe_3 <- boot.ci(boot_obj, type="bca", index=19)$bca[4:5]
  
  ci_pae_rr_1 <- boot.ci(boot_obj, type="bca", index=20)$bca[4:5]
  ci_pae_rr_2 <- boot.ci(boot_obj, type="bca", index=21)$bca[4:5]
  ci_pae_rr_3 <- boot.ci(boot_obj, type="bca", index=22)$bca[4:5]
  
  ci_pae_rd_1 <- boot.ci(boot_obj, type="bca", index=23)$bca[4:5]
  ci_pae_rd_2 <- boot.ci(boot_obj, type="bca", index=24)$bca[4:5]
  ci_pae_rd_3 <- boot.ci(boot_obj, type="bca", index=25)$bca[4:5]
  
  ci_paf_1 <- boot.ci(boot_obj, type="bca", index=26)$bca[4:5]
  ci_paf_2 <- boot.ci(boot_obj, type="bca", index=27)$bca[4:5]
  ci_paf_3 <- boot.ci(boot_obj, type="bca", index=28)$bca[4:5]
  
  # --- 3. Combine point estimates + 95% CIs --------
  res <- c(
    
    exp_margin_01    = full_est["exp_margin_01"],
    exp_margin_01_lb   = ci_exp_margin_01[1],
    exp_margin_01_ub   = ci_exp_margin_01[2],
    
    exp_margin_11    = full_est["exp_margin_11"],
    exp_margin_11_lb   = ci_exp_margin_11[1],
    exp_margin_11_ub   = ci_exp_margin_11[2],
    
    exp_margin_02    = full_est["exp_margin_02"],
    exp_margin_02_lb   = ci_exp_margin_02[1],
    exp_margin_02_ub   = ci_exp_margin_02[2],
    
    exp_margin_22    = full_est["exp_margin_22"],
    exp_margin_22_lb   = ci_exp_margin_22[1],
    exp_margin_22_ub   = ci_exp_margin_22[2],
    
    exp_margin_03    = full_est["exp_margin_03"],
    exp_margin_03_lb   = ci_exp_margin_03[1],
    exp_margin_03_ub   = ci_exp_margin_03[2],
    
    exp_margin_33    = full_est["exp_margin_33"],
    exp_margin_33_lb   = ci_exp_margin_33[1],
    exp_margin_33_ub   = ci_exp_margin_33[2],
    
    pop_margin_0    = full_est["pop_margin_0"],
    pop_margin_0_lb   = ci_pop_margin_0[1],
    pop_margin_0_ub   = ci_pop_margin_0[2],
    
    pop_margin_1    = full_est["pop_margin_1"],
    pop_margin_1_lb   = ci_pop_margin_1[1],
    pop_margin_1_ub   = ci_pop_margin_1[2],
    
    pop_margin_2    = full_est["pop_margin_2"],
    pop_margin_2_lb   = ci_pop_margin_2[1],
    pop_margin_2_ub   = ci_pop_margin_2[2],
    
    pop_margin_3    = full_est["pop_margin_3"],
    pop_margin_3_lb   = ci_pop_margin_3[1],
    pop_margin_3_ub   = ci_pop_margin_3[2],
    
    afe_rr_1      = full_est["afe_rr_1"],
    afe_rr_1_lb   = ci_afe_rr_1[1],
    afe_rr_1_ub   = ci_afe_rr_1[2],
    
    afe_rr_2      = full_est["afe_rr_2"],
    afe_rr_2_lb   = ci_afe_rr_2[1],
    afe_rr_2_ub   = ci_afe_rr_2[2],
    
    afe_rr_3      = full_est["afe_rr_3"],
    afe_rr_3_lb   = ci_afe_rr_3[1],
    afe_rr_3_ub   = ci_afe_rr_3[2],
    
    afe_rd_1      = full_est["afe_rd_1"],
    afe_rd_1_lb   = ci_afe_rd_1[1],
    afe_rd_1_ub   = ci_afe_rd_1[2],
    
    afe_rd_2      = full_est["afe_rd_2"],
    afe_rd_2_lb   = ci_afe_rd_2[1],
    afe_rd_2_ub   = ci_afe_rd_2[2],
    
    afe_rd_3      = full_est["afe_rd_3"],
    afe_rd_3_lb   = ci_afe_rd_3[1],
    afe_rd_3_ub   = ci_afe_rd_3[2],
    
    afe_1      = full_est["afe_1"],
    afe_1_lb   = ci_afe_1[1],
    afe_1_ub   = ci_afe_1[2],
    
    afe_2      = full_est["afe_2"],
    afe_2_lb   = ci_afe_2[1],
    afe_2_ub   = ci_afe_2[2],
    
    afe_3      = full_est["afe_3"],
    afe_3_lb   = ci_afe_3[1],
    afe_3_ub   = ci_afe_3[2],
    
    pae_rr_1      = full_est["pae_rr_1"],
    pae_rr_1_lb   = ci_pae_rr_1[1],
    pae_rr_1_ub   = ci_pae_rr_1[2],
    
    pae_rr_2      = full_est["pae_rr_2"],
    pae_rr_2_lb   = ci_pae_rr_2[1],
    pae_rr_2_ub   = ci_pae_rr_2[2],
    
    pae_rr_3      = full_est["pae_rr_3"],
    pae_rr_3_lb   = ci_pae_rr_3[1],
    pae_rr_3_ub   = ci_pae_rr_3[2],
    
    pae_rd_1      = full_est["pae_rd_1"],
    pae_rd_1_lb   = ci_pae_rd_1[1],
    pae_rd_1_ub   = ci_pae_rd_1[2],
    
    pae_rd_2      = full_est["pae_rd_2"],
    pae_rd_2_lb   = ci_pae_rd_2[1],
    pae_rd_2_ub   = ci_pae_rd_2[2],
    
    pae_rd_3      = full_est["pae_rd_3"],
    pae_rd_3_lb   = ci_pae_rd_3[1],
    pae_rd_3_ub   = ci_pae_rd_3[2],
    
    paf_1      = full_est["paf_1"],
    paf_1_lb   = ci_paf_1[1],
    paf_1_ub   = ci_paf_1[2],
    
    paf_2      = full_est["paf_2"],
    paf_2_lb   = ci_paf_2[1],
    paf_2_ub   = ci_paf_2[2],
    
    paf_3      = full_est["paf_3"],
    paf_3_lb   = ci_paf_3[1],
    paf_3_ub   = ci_paf_3[2]
    
  )
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# --------------------------------------------------------------
#     Assemble primigravidee results into a data frame
# -------------------------------------------------------------
# Combine results into a dataframe
primi_df <- imap_dfr(results_list, ~{
  as_tibble_row(.x) %>%
    mutate(outcome = .y, 
           subgroup = "Primigravidae")
})

# Remove repeated column names
colnames(primi_df) <- sub("^([^.]+)\\.\\1$", "\\1", colnames(primi_df)) 

# Reorder so outcome and subgroup appear first
primi_df <- primi_df[, c("outcome", "subgroup", setdiff(colnames(primi_df), c("outcome", "subgroup")))]

# Create long database (a little too difficult to reshape all pigment categories; sort out each pigment category and rbind) -----
mild <- primi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_01, exp_margin_01_lb, exp_margin_01_ub, 
                exp_margin_11,  exp_margin_11_lb,  exp_margin_11_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_1, pop_margin_1_lb, pop_margin_1_ub, 
                afe_rr_1, afe_rr_1_lb, afe_rr_1_ub, 
                afe_rd_1, afe_rd_1_lb, afe_rd_1_ub, 
                afe_1, afe_1_lb, afe_1_ub, 
                pae_rr_1, pae_rr_1_lb, pae_rr_1_ub, 
                pae_rd_1, pae_rd_1_lb, pae_rd_1_ub, 
                paf_1, paf_1_lb, paf_1_ub
  ) %>%
  rename_with(~ gsub("_01", "_0", .x)) %>%                # _01 → _0
  rename_with(~ gsub("_11", "", .x)) %>%                  # _11 → removed
  rename_with(~ gsub("_1(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _1 only if standalone
  mutate(pigment_code = 1)


moderate <- primi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_02, exp_margin_02_lb, exp_margin_02_ub, 
                exp_margin_22, exp_margin_22_lb, exp_margin_22_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_2, pop_margin_2_lb, pop_margin_2_ub, 
                afe_rr_2, afe_rr_2_lb, afe_rr_2_ub, 
                afe_rd_2, afe_rd_2_lb, afe_rd_2_ub, 
                afe_2, afe_2_lb, afe_2_ub, 
                pae_rr_2, pae_rr_2_lb, pae_rr_2_ub, 
                pae_rd_2, pae_rd_2_lb, pae_rd_2_ub, 
                paf_2, paf_2_lb, paf_2_ub
  ) %>%
  rename_with(~ gsub("_02", "_0", .x)) %>%                # _02 → _0
  rename_with(~ gsub("_22", "", .x)) %>%                  # _22 → removed
  rename_with(~ gsub("_2(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _2 only if standalone
  mutate(pigment_code = 2)

severe <- primi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_03, exp_margin_03_lb, exp_margin_03_ub, 
                exp_margin_33, exp_margin_33_lb, exp_margin_33_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_3, pop_margin_3_lb, pop_margin_3_ub, 
                afe_rr_3, afe_rr_3_lb, afe_rr_3_ub, 
                afe_rd_3, afe_rd_3_lb, afe_rd_3_ub, 
                afe_3, afe_3_lb, afe_3_ub, 
                pae_rr_3, pae_rr_3_lb, pae_rr_3_ub, 
                pae_rd_3, pae_rd_3_lb, pae_rd_3_ub, 
                paf_3, paf_3_lb, paf_3_ub
  ) %>%
  rename_with(~ gsub("_03", "_0", .x)) %>%                # _03 → _0
  rename_with(~ gsub("_33", "", .x)) %>%                  # _33 → removed
  rename_with(~ gsub("_3(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _3 only if standalone
  mutate(pigment_code = 3)

primi_df_final <- rbind(mild, moderate, severe)
print(primi_df_final)



# ------------------------------------------------------------------
#   Multigravidae estimates: Loop over outcomes + bootstrap CIs
# ------------------------------------------------------------------

# Store results for each outcome
results_list <- list()

# Set seed for reproducibility
set.seed(123)

for (y in outcomes) {
  
  # --- 1. Full dataset model for point estimates -----
  full_est <- gravid_compute_metrics(db, y, 1:nrow(db), grav_subgroup = "Multigravida")
  
  # --- 2. Bootstrap for BCa CIs -----
  boot_obj <- boot(
    data = db,
    statistic = function(d, i) gravid_compute_metrics(d, y, i, grav_subgroup = "Multigravida"),
    R = 1000,
    parallel = 'multicore',
    ncpus = parallel::detectCores() - 2
  )
  
  ci_exp_margin_01 <- boot.ci(boot_obj, type="bca", index=1)$bca[4:5]
  ci_exp_margin_11 <- boot.ci(boot_obj, type="bca", index=2)$bca[4:5]
  
  ci_exp_margin_02 <- boot.ci(boot_obj, type="bca", index=3)$bca[4:5]
  ci_exp_margin_22 <- boot.ci(boot_obj, type="bca", index=4)$bca[4:5]
  
  ci_exp_margin_03 <- boot.ci(boot_obj, type="bca", index=5)$bca[4:5]
  ci_exp_margin_33 <- boot.ci(boot_obj, type="bca", index=6)$bca[4:5]
  
  ci_pop_margin_0 <- boot.ci(boot_obj, type="bca", index=7)$bca[4:5]
  ci_pop_margin_1 <- boot.ci(boot_obj, type="bca", index=8)$bca[4:5]
  ci_pop_margin_2 <- boot.ci(boot_obj, type="bca", index=9)$bca[4:5]
  ci_pop_margin_3 <- boot.ci(boot_obj, type="bca", index=10)$bca[4:5]
  
  ci_afe_rr_1 <- boot.ci(boot_obj, type="bca", index=11)$bca[4:5]
  ci_afe_rr_2 <- boot.ci(boot_obj, type="bca", index=12)$bca[4:5]
  ci_afe_rr_3 <- boot.ci(boot_obj, type="bca", index=13)$bca[4:5]
  
  ci_afe_rd_1 <- boot.ci(boot_obj, type="bca", index=14)$bca[4:5]
  ci_afe_rd_2 <- boot.ci(boot_obj, type="bca", index=15)$bca[4:5]
  ci_afe_rd_3 <- boot.ci(boot_obj, type="bca", index=16)$bca[4:5]
  
  ci_afe_1 <- boot.ci(boot_obj, type="bca", index=17)$bca[4:5]
  ci_afe_2 <- boot.ci(boot_obj, type="bca", index=18)$bca[4:5]
  ci_afe_3 <- boot.ci(boot_obj, type="bca", index=19)$bca[4:5]
  
  ci_pae_rr_1 <- boot.ci(boot_obj, type="bca", index=20)$bca[4:5]
  ci_pae_rr_2 <- boot.ci(boot_obj, type="bca", index=21)$bca[4:5]
  ci_pae_rr_3 <- boot.ci(boot_obj, type="bca", index=22)$bca[4:5]
  
  ci_pae_rd_1 <- boot.ci(boot_obj, type="bca", index=23)$bca[4:5]
  ci_pae_rd_2 <- boot.ci(boot_obj, type="bca", index=24)$bca[4:5]
  ci_pae_rd_3 <- boot.ci(boot_obj, type="bca", index=25)$bca[4:5]
  
  ci_paf_1 <- boot.ci(boot_obj, type="bca", index=26)$bca[4:5]
  ci_paf_2 <- boot.ci(boot_obj, type="bca", index=27)$bca[4:5]
  ci_paf_3 <- boot.ci(boot_obj, type="bca", index=28)$bca[4:5]
  
  # --- 3. Combine point estimates (from full data) + BCa CIs -------
  res <- c(
    
    exp_margin_01    = full_est["exp_margin_01"],
    exp_margin_01_lb   = ci_exp_margin_01[1],
    exp_margin_01_ub   = ci_exp_margin_01[2],
    
    exp_margin_11    = full_est["exp_margin_11"],
    exp_margin_11_lb   = ci_exp_margin_11[1],
    exp_margin_11_ub   = ci_exp_margin_11[2],
    
    exp_margin_02    = full_est["exp_margin_02"],
    exp_margin_02_lb   = ci_exp_margin_02[1],
    exp_margin_02_ub   = ci_exp_margin_02[2],
    
    exp_margin_22    = full_est["exp_margin_22"],
    exp_margin_22_lb   = ci_exp_margin_22[1],
    exp_margin_22_ub   = ci_exp_margin_22[2],
    
    exp_margin_03    = full_est["exp_margin_03"],
    exp_margin_03_lb   = ci_exp_margin_03[1],
    exp_margin_03_ub   = ci_exp_margin_03[2],
    
    exp_margin_33    = full_est["exp_margin_33"],
    exp_margin_33_lb   = ci_exp_margin_33[1],
    exp_margin_33_ub   = ci_exp_margin_33[2],
    
    pop_margin_0    = full_est["pop_margin_0"],
    pop_margin_0_lb   = ci_pop_margin_0[1],
    pop_margin_0_ub   = ci_pop_margin_0[2],
    
    pop_margin_1    = full_est["pop_margin_1"],
    pop_margin_1_lb   = ci_pop_margin_1[1],
    pop_margin_1_ub   = ci_pop_margin_1[2],
    
    pop_margin_2    = full_est["pop_margin_2"],
    pop_margin_2_lb   = ci_pop_margin_2[1],
    pop_margin_2_ub   = ci_pop_margin_2[2],
    
    pop_margin_3    = full_est["pop_margin_3"],
    pop_margin_3_lb   = ci_pop_margin_3[1],
    pop_margin_3_ub   = ci_pop_margin_3[2],
    
    afe_rr_1      = full_est["afe_rr_1"],
    afe_rr_1_lb   = ci_afe_rr_1[1],
    afe_rr_1_ub   = ci_afe_rr_1[2],
    
    afe_rr_2      = full_est["afe_rr_2"],
    afe_rr_2_lb   = ci_afe_rr_2[1],
    afe_rr_2_ub   = ci_afe_rr_2[2],
    
    afe_rr_3      = full_est["afe_rr_3"],
    afe_rr_3_lb   = ci_afe_rr_3[1],
    afe_rr_3_ub   = ci_afe_rr_3[2],
    
    afe_rd_1      = full_est["afe_rd_1"],
    afe_rd_1_lb   = ci_afe_rd_1[1],
    afe_rd_1_ub   = ci_afe_rd_1[2],
    
    afe_rd_2      = full_est["afe_rd_2"],
    afe_rd_2_lb   = ci_afe_rd_2[1],
    afe_rd_2_ub   = ci_afe_rd_2[2],
    
    afe_rd_3      = full_est["afe_rd_3"],
    afe_rd_3_lb   = ci_afe_rd_3[1],
    afe_rd_3_ub   = ci_afe_rd_3[2],
    
    afe_1      = full_est["afe_1"],
    afe_1_lb   = ci_afe_1[1],
    afe_1_ub   = ci_afe_1[2],
    
    afe_2      = full_est["afe_2"],
    afe_2_lb   = ci_afe_2[1],
    afe_2_ub   = ci_afe_2[2],
    
    afe_3      = full_est["afe_3"],
    afe_3_lb   = ci_afe_3[1],
    afe_3_ub   = ci_afe_3[2],
    
    pae_rr_1      = full_est["pae_rr_1"],
    pae_rr_1_lb   = ci_pae_rr_1[1],
    pae_rr_1_ub   = ci_pae_rr_1[2],
    
    pae_rr_2      = full_est["pae_rr_2"],
    pae_rr_2_lb   = ci_pae_rr_2[1],
    pae_rr_2_ub   = ci_pae_rr_2[2],
    
    pae_rr_3      = full_est["pae_rr_3"],
    pae_rr_3_lb   = ci_pae_rr_3[1],
    pae_rr_3_ub   = ci_pae_rr_3[2],
    
    pae_rd_1      = full_est["pae_rd_1"],
    pae_rd_1_lb   = ci_pae_rd_1[1],
    pae_rd_1_ub   = ci_pae_rd_1[2],
    
    pae_rd_2      = full_est["pae_rd_2"],
    pae_rd_2_lb   = ci_pae_rd_2[1],
    pae_rd_2_ub   = ci_pae_rd_2[2],
    
    pae_rd_3      = full_est["pae_rd_3"],
    pae_rd_3_lb   = ci_pae_rd_3[1],
    pae_rd_3_ub   = ci_pae_rd_3[2],
    
    paf_1      = full_est["paf_1"],
    paf_1_lb   = ci_paf_1[1],
    paf_1_ub   = ci_paf_1[2],
    
    paf_2      = full_est["paf_2"],
    paf_2_lb   = ci_paf_2[1],
    paf_2_ub   = ci_paf_2[2],
    
    paf_3      = full_est["paf_3"],
    paf_3_lb   = ci_paf_3[1],
    paf_3_ub   = ci_paf_3[2]
    
  )
  
  # --- 4. Store estimates ------
  results_list[[y]] <- res
  
}

# --------------------------------------------------------------
#     Assemble multigravidae results into a data frame
# -------------------------------------------------------------
# Combine results into a dataframe
multi_df <- imap_dfr(results_list, ~{
  as_tibble_row(.x) %>%
    mutate(outcome = .y, 
           subgroup = "Multigravidae")
})

# Remove repeated column names
colnames(multi_df) <- sub("^([^.]+)\\.\\1$", "\\1", colnames(multi_df)) 

# Reorder so outcome and subgroup appear first
multi_df <- multi_df[, c("outcome", "subgroup", setdiff(colnames(multi_df), c("outcome", "subgroup")))]

# Create long database (a little too difficult to reshape all pigment categories; sort out each pigment category and rbind) -----
mild <- multi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_01, exp_margin_01_lb, exp_margin_01_ub, 
                exp_margin_11,  exp_margin_11_lb,  exp_margin_11_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_1, pop_margin_1_lb, pop_margin_1_ub, 
                afe_rr_1, afe_rr_1_lb, afe_rr_1_ub, 
                afe_rd_1, afe_rd_1_lb, afe_rd_1_ub, 
                afe_1, afe_1_lb, afe_1_ub, 
                pae_rr_1, pae_rr_1_lb, pae_rr_1_ub, 
                pae_rd_1, pae_rd_1_lb, pae_rd_1_ub, 
                paf_1, paf_1_lb, paf_1_ub
  ) %>%
  rename_with(~ gsub("_01", "_0", .x)) %>%                # _01 → _0
  rename_with(~ gsub("_11", "", .x)) %>%                  # _11 → removed
  rename_with(~ gsub("_1(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _1 only if standalone
  mutate(pigment_code = 1)


moderate <- multi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_02, exp_margin_02_lb, exp_margin_02_ub, 
                exp_margin_22, exp_margin_22_lb, exp_margin_22_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_2, pop_margin_2_lb, pop_margin_2_ub, 
                afe_rr_2, afe_rr_2_lb, afe_rr_2_ub, 
                afe_rd_2, afe_rd_2_lb, afe_rd_2_ub, 
                afe_2, afe_2_lb, afe_2_ub, 
                pae_rr_2, pae_rr_2_lb, pae_rr_2_ub, 
                pae_rd_2, pae_rd_2_lb, pae_rd_2_ub, 
                paf_2, paf_2_lb, paf_2_ub
  ) %>%
  rename_with(~ gsub("_02", "_0", .x)) %>%                # _02 → _0
  rename_with(~ gsub("_22", "", .x)) %>%                  # _22 → removed
  rename_with(~ gsub("_2(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _2 only if standalone
  mutate(pigment_code = 2)

severe <- multi_df %>%
  dplyr::select(outcome, subgroup, 
                exp_margin_03, exp_margin_03_lb, exp_margin_03_ub, 
                exp_margin_33, exp_margin_33_lb, exp_margin_33_ub,
                pop_margin_0, pop_margin_0_lb, pop_margin_0_ub, 
                pop_margin_3, pop_margin_3_lb, pop_margin_3_ub, 
                afe_rr_3, afe_rr_3_lb, afe_rr_3_ub, 
                afe_rd_3, afe_rd_3_lb, afe_rd_3_ub, 
                afe_3, afe_3_lb, afe_3_ub, 
                pae_rr_3, pae_rr_3_lb, pae_rr_3_ub, 
                pae_rd_3, pae_rd_3_lb, pae_rd_3_ub, 
                paf_3, paf_3_lb, paf_3_ub
  ) %>%
  rename_with(~ gsub("_03", "_0", .x)) %>%                # _03 → _0
  rename_with(~ gsub("_33", "", .x)) %>%                  # _33 → removed
  rename_with(~ gsub("_3(?![0-9])", "", .x, perl = TRUE)) %>%  # drop _3 only if standalone
  mutate(pigment_code = 3)

multi_df_final <- rbind(mild, moderate, severe)
print(multi_df_final)


# =====================================================
#   PART 3: COMBINE ALL RESULTS AND EXPORT 
# =====================================================
finalresults <- rbind(overall_df_final, primi_df_final, multi_df_final) %>%
  mutate(
    rexp_amongexp_descrip = paste0(round(100*exp_margin, 1), "% [", round(100*exp_margin_lb, 1), "%, ", round(100*exp_margin_ub, 1), "%]"),
    runexp_descrip = paste0(round(100*exp_margin_0, 1), "% [", round(100*exp_margin_0_lb, 1), "%, ", round(100*exp_margin_0_ub, 1), "%]"),
    rnc_pop_descrip = paste0(round(100*pop_margin_0, 1), "% [", round(100*pop_margin_0_lb, 1), "%, ", round(100*pop_margin_0_ub, 1), "%]"),
    runexp_pop_descrip = paste0(round(100*pop_margin, 1), "% [", round(100*pop_margin_lb, 1), "%, ", round(100*pop_margin_ub, 1), "%]"),
    
    afe_rr_descrip = paste0(round(afe_rr,2), " [", round(afe_rr_lb, 2), ", ", round(afe_rr_ub, 2), "]"), 
    afe_rd_descrip = paste0(round(100*afe_rd,1), "% [", round(100*afe_rd_lb, 1), "%, ", round(100*afe_rd_ub, 1), "%]"),
    afe_descrip = paste0(round(100*afe, 1), "% [", round(100*afe_lb, 1), "%, ", round(100*afe_ub, 1), "%]"),
    
    pae_rr_descrip = paste0(round(pae_rr,2), " [", round(pae_rr_lb, 2), ", ", round(pae_rr_ub, 2), "]"), 
    pae_rd_descrip = paste0(round(100*pae_rd,1), "% [", round(100*pae_rd_lb, 1), "%, ", round(100*pae_rd_ub, 1), "%]"), 
    paf_descrip = paste0(round(100*paf, 1), "% [", round(100*paf_lb, 1), "%, ", round(100*paf_ub, 1), "%]"),
    
    outcome_code = case_when(
      outcome == "preterm"  ~ 1,
      outcome == "SGA"      ~ 2,
      outcome == "LBWdich"  ~ 3), 
    subgroup_code = case_when(
      subgroup == "Overall" ~ 1, 
      subgroup == "Primigravidae" ~ 2, 
      subgroup == "Multigravidae" ~ 3),
    pigment_type = case_when(
      pigment_code == 1 ~ "Mild", 
      pigment_code == 2 ~ "Moderate", 
      pigment_code == 3 ~ "Severe"
    ))

# ----------------------------
#   Export results 
# ----------------------------
write.xlsx(finalresults, file="Results/DPSP_SensAnalysis_Pigment_NoActiveInf_Gcomp_Boot.xlsx")

