
##################################################################################
##	  Quantifying the contributions of active and past placental malaria infection
##    to low birthweight in a high transmission area of Uganda
##	  14_Sensitivity Analyses - Assessing for potential positivity violations 
##
##  This script assessed for potential positivity violations using propensity scores 
##  to assess whether there is adequate overlap in covariate distribution between
##  exposed and unexposed groups. 
##      - Positivity checks were done for active infection and severe pigment deposition
##      
##################################################################################

rm(list=ls())

# ----------------------------
#   User inputs: 
# ----------------------------
# Edit line below to set path directory user's own path directory: 
setwd('[insert path directory]') 

# Create "Figures" folder if it doesn't already exist
if (!dir.exists("Figures")) {
  dir.create("Figures")
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
db <- read.dta13('Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta', nonint.factors = TRUE) %>%
  mutate(propfibrincat_sim_cat = case_when(
    propfibrincat_sim == 0 ~ "None", 
    propfibrincat_sim == 1 ~ "<10%", 
    propfibrincat_sim == 2 ~ "10 to <30%", 
    propfibrincat_sim == 3 ~ ">=30%"))

# -------------------------------------------------
#   Define covariates included in PS model 
# -------------------------------------------------
# Excludes graviddich, which will added in the function for Overall estimates
base_covars <- c(
  "treatmentarm",
  "age_enroll", "age2",
  "wealthcat_hd",
  "bmi_enroll_hd", "bmi_enroll_hd2",
  "GAenroll", "gaenroll2",
  "educ_cat",
  "male"
)

# -----------------------------------------------------------------------------------
#   Initialize final dataset where propensity scores will be added to sequentially
# -----------------------------------------------------------------------------------
finaldb <- db


# ===============================================================================
#   GENERATE FUNCTIONS FOR ASSESSING POSITIVITY VIOLATIONS/COMMON SUPPORT
# ===============================================================================

# -------------------------------------------------------------------------------
#   common_support_quant: Identify common support bounds based on quantiles 
#
#     This function returns the lower and upper bounds of the quantiles where 
#     there is sufficient PS overlap between exposed and unexposed individuals. 
# -----------------------------------------------------------------------------
common_support_quant <- function(ps, a, q = 0.01) {
  # q: fraction of each tail to trim within each exposure group
  #    (e.g., q = 0.01 trims 1% in each tail before defining common support)
  
  lower0 <- stats::quantile(ps[a == 0], q,      na.rm = TRUE)
  lower1 <- stats::quantile(ps[a == 1], q,      na.rm = TRUE)
  upper0 <- stats::quantile(ps[a == 0], 1 - q,  na.rm = TRUE)
  upper1 <- stats::quantile(ps[a == 1], 1 - q,  na.rm = TRUE)
  
  L <- max(lower0, lower1)  # lower end of common support
  U <- min(upper0, upper1)  # upper end of common support
  
  c(lower = L, upper = U)
}


# ------------------------------------------------------------------------------  
#  run_ps_positivity: Function to run PS model, make density plot, and flag common support
#
# Arguments: 
#   data          : data.frame
#   exposure      : name of exposure variable (string), e.g. "activeHP" or "severe"
#   gravidity     : "Overall", "Primigravidae", "Multigravidae"
#   base_covars   : character vector of baseline covariate names
#   ps_var        : name of PS variable to create (string), e.g. "actinf_prop"
#   pdf_file      : file path for PDF plot
#   title_prefix  : text prefix for plot title
#   two_sided     : if TRUE, uses both lower and upper bounds; if FALSE,
#                   uses only the lower bound (matches one-sided trimming
#                   you were doing before)
#   indicator_name: name of variable returning 1 if individual's propensity score
#                   fell within the region of common support; 
#                   Used as an indicator value to restrict sensitivity analyses 
#                   to those with good overlap (common support)
# ------------------------------------------------------------------------------

run_ps_positivity <- function(data,
                              exposure,
                              gravidity = c("Overall", "Primigravida", "Multigravida"),
                              base_covars,
                              ps_var,
                              jpg_file,
                              title_prefix,
                              q = 0.01, 
                              two_sided=FALSE,
                              indicator_name = NULL) {
  
  gravidity <- match.arg(gravidity)
  
  # 1. Subset by gravidity if stratified analyses -----
  if (gravidity != "Overall") {
    data <- data %>%
      dplyr::filter(graviddich == gravidity)
  }
  
  # 2. Build formula: include graviddich only in Overall models -----
  covars <- base_covars
  if (gravidity == "Overall") {
    covars <- c("graviddich", covars)
  }
  
  form <- as.formula(
    paste(exposure, "~", paste(covars, collapse = " + "))
  )
  
  # 3. Fit logistic PS model -----
  model <- glm(form, data = data, family = binomial(link = "logit"))
  
  # 4. Add PS variable to data -----
  data[[ps_var]] <- predict(model, type = "response")
  
  # 5. Compute common support bounds using quantiles in each exposure group ----
  ps_vec <- data[[ps_var]]
  a_vec  <- data[[exposure]]
  
  cs <- common_support_quant(ps = ps_vec, a = a_vec, q = q)
  cutoff_lower <- unname(cs["lower"])
  cutoff_upper <- unname(cs["upper"])
  
  # 6. Density plot of PS by observed exposure status -----
  plot_df <- data %>%
    dplyr::mutate(
      ps       = .data[[ps_var]],
      exposure = .data[[exposure]]
    )
  
  gravid_label <- if (gravidity == "Overall") "Overall" else gravidity
  
  p <- ggplot(plot_df, aes(x = ps, fill = factor(exposure))) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(
      values = c("0" = "steelblue", "1" = "tomato"),
      labels = c("Unexposed", "Exposed"),
      name   = "Observed exposure"
    ) +
    labs(
      x     = "Propensity score",
      y     = "Density",
      title = paste0(title_prefix, " (", gravid_label, ")"), 
      subtitle = paste0(
        "Cut-off: ",
        round(cutoff_upper, 2)
      )
    ) +
    theme_minimal() + 
    geom_vline(xintercept = cutoff_upper, linetype = "dashed", color = "red")
  
  if (two_sided) {
    p <- p + geom_vline(xintercept = cutoff_lower, linetype = "dashed", color = "red")
  }
  
  ggsave(file = jpg_file, p, height = 3, width = 7, dpi=300)
  
  # 7. Flag high-PS observations (one-sided trimming by default) ----
  if (two_sided) {
    flagged <- data %>%
      dplyr::filter(.data[[ps_var]] >= cutoff_lower | .data[[ps_var]] <= cutoff_upper) %>%
      dplyr::transmute(id, !!indicator_name := 1L)
  } else {
    flagged <- data %>%
      dplyr::filter(.data[[ps_var]] <= cutoff_upper) %>%
      dplyr::transmute(id, !!indicator_name := 1L)
  }
  
  list(
    data    = data,
    model   = model,
    flagged = flagged
  )
}

# ============================================================================
# ACTIVE PLACENTAL INFECTION (activeHP)
# ============================================================================

# ----------------------------
#  Overall 
# ----------------------------
res_act_overall <- run_ps_positivity(
  data          = db,
  exposure      = "activeHP",
  gravidity     = "Overall",
  base_covars   = base_covars,
  ps_var        = "actinf_prop",
  jpg_file      = "Figures/PositivityCheck_ActInf_Overall_Cutoff.jpg",
  title_prefix  = "Positivity check: Active Infection in Overall Population",
  q             = 0.01,                      
  two_sided     = FALSE,                     
  indicator_name = "overall_actinf_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_act_overall$flagged, by = "id")


# ----------------------------
#  Primigravidae
# ----------------------------
res_act_primi <- run_ps_positivity(
  data          = db,
  exposure      = "activeHP",
  gravidity     = "Primigravida",
  base_covars   = base_covars,
  ps_var        = "actinf_prop",
  jpg_file      = "Figures/PositivityCheck_ActInf_Primi_Cutoff.jpg",
  title_prefix  = "Positivity check: Active Infection in Primigravidae",
  q             = 0.01,
  two_sided     = FALSE,
  indicator_name = "primi_actinf_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_act_primi$flagged, by = "id")


# ----------------------------
#  Multigravidae
# ----------------------------
res_act_multi <- run_ps_positivity(
  data          = db,
  exposure      = "activeHP",
  gravidity     = "Multigravida",
  base_covars   = base_covars,
  ps_var        = "actinf_prop",
  jpg_file      = "Figures/PositivityCheck_ActInf_Multi_Cutoff.jpg",
  title_prefix  = "Positivity check: Active Infection in Multigravidae",
  q             = 0.01,
  two_sided     = FALSE,
  indicator_name = "multi_actinf_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_act_multi$flagged, by = "id")




# ============================================================================
#     SEVERE PIGMENT DEPOSITION (>=30% vs none)
# ============================================================================

# Create a new variable "severe" as severe pigment vs none -------
severe_db <- db %>%
  dplyr::filter(propfibrincat == "None" | propfibrincat == ">=30%") %>%
  dplyr::mutate(severe = ifelse(propfibrincat == ">=30%", 1L, 0L))

# ----------------------------
#  Overall 
# ----------------------------
res_severe_overall <- run_ps_positivity(
  data          = severe_db,
  exposure      = "severe",
  gravidity     = "Overall",
  base_covars   = base_covars,
  ps_var        = "severe_prop",
  jpg_file      = "Figures/PositivityCheck_SeverePigment_Overall_Cutoff.jpg",
  title_prefix  = "Positivity check: Severe Pigment in Overall Population",
  q             = 0.01,
  two_sided     = FALSE,
  indicator_name = "overall_severe_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_severe_overall$flagged, by = "id")


# ----------------------------
#  Primigravidae 
# ----------------------------
res_severe_primi <- run_ps_positivity(
  data          = severe_db,
  exposure      = "severe",
  gravidity     = "Primigravida",
  base_covars   = base_covars,
  ps_var        = "severe_prop",
  jpg_file      = "Figures/PositivityCheck_SeverePigment_Primi_Cutoff.jpg",
  title_prefix  = "Positivity check: Severe Pigment in Primigravidae",
  q             = 0.01,
  two_sided     = FALSE,
  indicator_name = "primi_severe_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_severe_primi$flagged, by = "id")


# ----------------------------
#  Multigravidae 
# ----------------------------
res_severe_multi <- run_ps_positivity(
  data          = severe_db,
  exposure      = "severe",
  gravidity     = "Multigravida",
  base_covars   = base_covars,
  ps_var        = "severe_prop",
  jpg_file      = "Figures/PositivityCheck_SeverePigment_Multi_Cutoff.jpg",
  title_prefix  = "Positivity check: Severe Pigment in Multigravidae",
  q             = 0.01,
  two_sided     = FALSE,
  indicator_name = "multi_severe_pscore"
)

finaldb <- finaldb %>%
  dplyr::left_join(res_severe_multi$flagged, by = "id")


# ============================================================================
# EXPORT: indicators for restricted analyses
# ============================================================================

export_pscore <- finaldb %>%
  dplyr::select(
    id,
    overall_actinf_pscore,
    primi_actinf_pscore,
    multi_actinf_pscore,
    overall_severe_pscore,
    primi_severe_pscore,
    multi_severe_pscore
  )

write.xlsx(
  export_pscore,
  file = "Databases/Created Databases/Propensity Scores_Indicator for Restricted Analyses.xlsx"
)


