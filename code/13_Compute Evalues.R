
##################################################################################
##	Quantifying the contributions of active and past placental malaria infection
##  to low birthweight in a high transmission area of Uganda
##  13_Compute Evalues for Exposure Effect among Exposed 
##
##  This script contains a function (add_evalues) that computes E-values for RRexp
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
library(EValue)
library(openxlsx)
library(purrr)


# ----------------------------
#   Upload databases
# ----------------------------
active <- read.xlsx('Results/DPSP_ActInf_Gcomp_Boot.xlsx') 
pigment <- read.xlsx('Results/DPSP_Pigment_Gcomp_Boot.xlsx')
active_sev <- read.xlsx('Results/DPSP_SensAnalysis_AnyActorSev_Gcomp_Boot.xlsx')
active_nosev <- read.xlsx('Results/DPSP_SensAnalysis_ActInf_NoSevere_Gcomp_Boot.xlsx')
pig_noactive <- read.xlsx('Results/DPSP_SensAnalysis_Pigment_NoActiveInf_Gcomp_Boot.xlsx')
anypm <- read.xlsx('Results/DPSP_SensAnalysis_AnyPM_Gcomp_Boot.xlsx')

# -------------------------------------------------------------------
#   add_evalues: Adds e-values to existing results databases and exports
#                them as an Excel file with subscript "_FINAL.xlsx"
# -------------------------------------------------------------------
add_evalues <- function(df, label) {
  # Initialize evalue vectors 
  evalue_point <- evalue_lo <- evalue_hi <- numeric(nrow(df))
  
  # Loop through rows of dataframe
  for (i in 1:nrow(df)) {
    if (any(is.na(c(df$afe_rr[i], df$afe_rr_lb[i], df$afe_rr_ub[i])))) next
    
    ev <- evalues.RR(est = df$afe_rr[i], lo = df$afe_rr_lb[i], hi = df$afe_rr_ub[i], true = 1)
    evalue_point[i] <- ev[2, 1]
    evalue_lo[i] <- ev[2, 2]
    evalue_hi[i] <- ev[2, 3]
  }
  
  # Add only E-value columns to dataframe
  df <- df %>%
    mutate(evalue_point = evalue_point,
           evalue_lo    = evalue_lo,
           evalue_hi    = evalue_hi)
  
  # Export excel file 
  write.xlsx(df, file = paste0("Results/DPSP_", label, "_Gcomp_Boot_FINAL.xlsx")) 
  
  return(df)
}


# -------------------------------------------------------------------
#   Compute evalues and export excel files
# -------------------------------------------------------------------
active  <- add_evalues(active, "ActInf")
pigment <- add_evalues(pigment, "Pigment")
active_sev <- add_evalues(active_sev, "SensAnalysis_AnyActorSev")
active_nosev <- add_evalues(active_nosev, "SensAnalysis_ActInf_NoSevere")
pig_noactive <- add_evalues(pig_noactive, "SensAnalysis_Pigment_NoActiveInf")
anypm   <- add_evalues(anypm, "SensAnalysis_AnyPM")




