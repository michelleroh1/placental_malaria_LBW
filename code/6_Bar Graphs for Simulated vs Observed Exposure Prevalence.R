
##################################################################################
##	 Quantifying the contributions of active and past placental malaria infection 
##   to low birthweight in a high transmission area of Uganda
##   6_Bar Graphs for Simulated vs Observed Exposure Prevalence 
##################################################################################

rm(list=ls())

# ----------------------------
#   User inputs: 
# ----------------------------
setwd('[insert path directory]')

# Create "Results" folder if it doesn't already exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Create "Figures" folder if it doesn't already exist
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# ----------------------------
#   Upload libraries
# ----------------------------
library(openxlsx)
library(readstata13)
library(dplyr)
library(tidyr)


# ----------------------------
#   Upload databases
# ----------------------------
db <- read.dta13("Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta")


# --------------------------------------------
#   Generate databases of prevalence 
# --------------------------------------------

# 1. OVERALL -----------

  # Simulated -----
  overall_sim <- db %>% 
    summarise(
      gravidity = "Overall",
      type = "Simulated", 
      `Any PM`= mean(anyplcmal_sim, na.rm=TRUE),
      Active = mean(activeHP_sim, na.rm=TRUE), 
      Mild = mean(propfibrincat_sim==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat_sim==2, na.rm=TRUE),
      Severe = mean(propfibrincat_sim==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )

  # Observed prevalence in IPTp-SP arm only ---
  overall_sp <- db %>%
    filter(treatmentarm==1) %>%
    summarise(
      gravidity = "Overall",
      type = "IPTp-SP arm only",
      `Any PM`= mean(anyplcmal, na.rm=TRUE),
      Active = mean(activeHP, na.rm=TRUE), 
      Mild = mean(propfibrincat==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat==2, na.rm=TRUE),
      Severe = mean(propfibrincat==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )

# 1. PRIMIGRAVIDAE -----------
  
  # Simulated -----
  primi_sim <- db %>% 
    filter(graviddich==1) %>%
    summarise(
      gravidity = "Primigravidae", 
      type = "Simulated", 
      `Any PM`= mean(anyplcmal_sim, na.rm=TRUE),
      Active = mean(activeHP_sim, na.rm=TRUE), 
      Mild = mean(propfibrincat_sim==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat_sim==2, na.rm=TRUE),
      Severe = mean(propfibrincat_sim==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )
  
  # Observed prevalence in IPTp-SP arm only ---
  primi_sp <- db %>%
    filter(treatmentarm==1 & graviddich==1) %>%
    summarise(
      gravidity = "Primigravidae", 
      type = "IPTp-SP arm only",
      `Any PM`= mean(anyplcmal, na.rm=TRUE),
      Active = mean(activeHP, na.rm=TRUE), 
      Mild = mean(propfibrincat==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat==2, na.rm=TRUE),
      Severe = mean(propfibrincat==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )
  
# 3. MULTIGRAVIDAE -----------
  
  # Simulated -----
  multi_sim <- db %>% 
    filter(graviddich==2) %>%
    summarise(
      gravidity = "Multigravidae", 
      type = "Simulated", 
      `Any PM`= mean(anyplcmal_sim, na.rm=TRUE),
      Active = mean(activeHP_sim, na.rm=TRUE), 
      Mild = mean(propfibrincat_sim==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat_sim==2, na.rm=TRUE),
      Severe = mean(propfibrincat_sim==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )
  
  # Observed prevalence in IPTp-SP arm only ---
  multi_sp <- db %>%
    filter(treatmentarm==1 & graviddich==2) %>%
    summarise(
      gravidity = "Multigravidae", 
      type = "IPTp-SP arm only",
      `Any PM`= mean(anyplcmal, na.rm=TRUE),
      Active = mean(activeHP, na.rm=TRUE), 
      Mild = mean(propfibrincat==1, na.rm=TRUE), 
      Moderate = mean(propfibrincat==2, na.rm=TRUE),
      Severe = mean(propfibrincat==3, na.rm=TRUE)
    ) %>%
    pivot_longer(
      cols = `Any PM`:Severe,
      names_to = "pm_type",
      values_to = "prev"
    )

# --------------------------------------------
#   Combine databases 
# --------------------------------------------
prev <- rbind(overall_sim, overall_sp, 
              primi_sim, primi_sp, 
              multi_sim, multi_sp) %>%
  mutate(alpha = ifelse(type=="Simulated", 1, 0.5), 
         type = factor(type, levels=c('Simulated', 'IPTp-SP arm only')), 
         pm_type = factor(pm_type, levels=c('Any PM', 
                                            'Active', 
                                            'Mild', 
                                            'Moderate', 
                                            'Severe')), 
         gravidity = factor(gravidity, levels=c('Overall', 'Primigravidae', 
                                                'Multigravidae')),
         prev = 100*prev)  
  
# ---------------------------------------------------------------
#   Plot bar graph of simulated vs observed prevalence 
# ---------------------------------------------------------------

p <- ggplot(prev, aes(x = pm_type, 
                      y = prev, 
                      fill = pm_type, 
                      group=type,
                      alpha = alpha)) +
       facet_wrap(~gravidity, ncol=1) +
       geom_bar(stat = "identity", position = position_dodge(width = 0.90)) +
       theme_minimal() +
       labs(x = "Placental Malaria Type", y = "Prevalence") +
       scale_alpha_continuous(range = c(0.5, 1)) +
       theme(legend.position = "none") +
       theme_minimal() +
       ylab('') +
       xlab('') +
       scale_fill_manual(values=c('#E5B84D', '#6EB5B5', '#D48C8C', 
                                  '#8C9A66', '#D94F3D')) +
       scale_y_continuous(limits=c(0, 100))
       theme(axis.text = element_text(size=11), 
             strip.text = element_text(size = 12))

ggsave("Figures/Simulated vs Observed Exp Prev.jpg", p, width = 8, height = 8, dpi = 300)

