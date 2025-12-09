
##################################################################################
##	  Quantifying the contributions of active and past placental malaria infection
##    to low birthweight in a high transmission area of Uganda
##    17_Generate Forest plots
## 
##  This script reads in g-comp results, reshapes them into long format, and 
##   generates forest plots by RRexp stratified by gravidity.
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
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)


# ----------------------------------------------------------------------------------
#   Prior to uploading databases, clean up estimates using below function: 
#   zero_afe_rows: Some of the estimates were highly unstable and generated implausible estimates
#   (e.g., RRexp = 0), if this is the case, this function makes all estimates in that row NA
# ----------------------------------------------------------------------------------
unstable_afe_rows <- function(df) {
  # If afe_rr isn't present, just return unchanged
  if (!"afe_rr" %in% names(df)) return(df)
  
  df %>%
    mutate(
      across(
        where(is.numeric),
        ~ ifelse(afe_rr < 0.01, NA_real_, .)
      )
    )
}

# ----------------------------
#   Upload databases
# ----------------------------
active <- read.xlsx('Results/DPSP_ActInf_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()
pigment <- read.xlsx('Results/DPSP_Pigment_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()
anypm <- read.xlsx('Results/DPSP_SensAnalysis_AnyPM_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()
sens_pigment <- read.xlsx('Results/DPSP_SensAnalysis_Pigment_NoActiveInf_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()
active_nopig <- read.xlsx('Results/DPSP_SensAnalysis_ActInf_NoSevere_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()
actsev <- read.xlsx('Results/DPSP_SensAnalysis_AnyActorSev_Gcomp_Boot_FINAL.xlsx') %>%
  unstable_afe_rows()


# ==================================
#   GENERATE FOREST PLOTS 
# ==================================

# --------------------------------------------------------------------------
#   Input parameters that will be repeated throughout forest plots 
# --------------------------------------------------------------------------
# Common factor levels 
outcome_levels  <- c("LBWdich", "SGA", "preterm")
subgroup_levels <- c("Overall", "Primigravidae", "Multigravidae")
est_type_levels <- c("pae_rr", "afe_rr")

# Custom palette 
pm_palette <- c("#729e9e", "#d08770", "#e9c46a")


# --------------------------------------------------------------------------
#   reshape_long: Function to convert wide RRexp results into long format
# --------------------------------------------------------------------------
reshape_long <- function(df) {
  
  df %>%
    select(
      outcome, subgroup,
      afe_rr, afe_rr_lb, afe_rr_ub,
      pae_rr, pae_rr_lb, pae_rr_ub
    ) %>%
    rename(
      afe_rr_est = afe_rr,
      pae_rr_est = pae_rr
    ) %>%
    pivot_longer(
      cols = c(afe_rr_est, afe_rr_lb, afe_rr_ub,
               pae_rr_est, pae_rr_lb, pae_rr_ub),
      names_to     = c("est_type", ".value"),
      names_pattern = "(.*_rr)_(est|lb|ub)"
    ) %>%
    mutate(
      outcome  = factor(outcome,  levels = outcome_levels),
      subgroup = factor(subgroup, levels = subgroup_levels),
      est_type = factor(est_type, levels = est_type_levels)
    )
}


# ---------------------------------------------------------------
#   make_forest_plot: Function to create forest plots 
# ---------------------------------------------------------------
make_forest_plot <- function(df_long,
                             file,
                             x_breaks,
                             x_limits,
                             height = 5.3,
                             width  = 6.6) {
  
  plot_df <- df_long %>%
    mutate(
      outcome  = factor(.data$outcome,  levels = outcome_levels),
      subgroup = factor(.data$subgroup, levels = subgroup_levels),
      est_type = factor(.data$est_type, levels = est_type_levels)
    )
  
  p <- ggplot(
    data = df_long,
    aes(
      x     = log(est),
      y     = outcome,
      col   = outcome,
      fill  = outcome,
      group = interaction(outcome, est_type),
      shape = est_type
    )
  ) +
    geom_point(
      size = 3.4,
      position = position_dodge(width = 0.7)
    ) +
    geom_pointrange(
      aes(xmin = log(lb), 
          xmax = log(ub)),
      position = position_dodge(width = 0.7),
      size = 1
    ) +
    facet_wrap(~ subgroup, ncol = 1) +
    theme_minimal() +
    geom_vline(xintercept = 0, lty = "dashed") +
    scale_color_manual(values = pm_palette) +
    scale_fill_manual(values  = c("white", "white", "white")) +
    scale_shape_manual(values = c(21, 16)) +  
    scale_x_continuous(
      breaks = log(x_breaks),
      labels = x_breaks
    ) +
    coord_cartesian(xlim = log(x_limits)) +
    theme(
      axis.line.x = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black")
    ) +
    ylab("") +
    xlab("Risk Ratio")
  
  ggsave(
    filename = file,
    plot = p,
    height = height,
    width = width,
    dpi = 300
  )
}


# ---------------------------------------
#   Primary analysis: Active infection 
# ---------------------------------------
active_long <- reshape_long(active)

make_forest_plot(
  df_long  = active_long,
  file     = file.path("Figures", "Active Inf RR_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8, 14),
  x_limits = c(0.15, 14)
)

# -----------------------------------------
#    Primary analysis: Mild pigment
# -----------------------------------------
mild_long <- pigment %>%
  filter(pigment_type == "Mild") %>%
  reshape_long()

make_forest_plot(
  df_long  = mild_long,
  file     = file.path("Figures", "Mild Pigment RR_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.20, 8)
)



# -----------------------------------------
#    Primary analysis: Moderate pigment
# -----------------------------------------
moderate_long <- pigment %>%
  filter(pigment_type == "Moderate") %>%
  reshape_long()

make_forest_plot(
  df_long  = moderate_long,
  file     = file.path("Figures", "Moderate Pigment RR_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.20, 8)
)


# -----------------------------------------
#    Primary analysis: Severe pigment
# -----------------------------------------
severe_long <- pigment %>%
  filter(pigment_type == "Severe") %>%
  reshape_long()

make_forest_plot(
  df_long  = severe_long,
  file     = file.path("Figures", "Severe Pigment RR_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.20, 8)
)


# -------------------------------------------------------
#    Sensitivity analysis: Any Parasites or Pigment
# -------------------------------------------------------
anypm_long <- reshape_long(anypm)

make_forest_plot(
  df_long  = anypm_long,
  file     = file.path("Figures", "Any PM RR_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.15, 14)
)


# -------------------------------------------------------------------
#    Sensitivity analysis: Active infection w/o severe pigment
# -------------------------------------------------------------------
sens_act_long <- reshape_long(active_nopig)

make_forest_plot(
  df_long  = sens_act_long,
  file     = file.path("Figures", "Sens_ActInf_NoSevere_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.15, 14)
)


# -------------------------------------------------------------------
#    Sensitivity analysis: Severe pigment w/o active infection
# -------------------------------------------------------------------
sens_sev_long <- sens_pigment %>%
  filter(pigment_code == 3) %>%
  reshape_long()

make_forest_plot(
  df_long  = sens_sev_long,
  file     = file.path("Figures", "Sens_Severe_NoActiveInf_By Gravidity.jpg"),
  x_breaks = c(0.25, 0.5, 1, 2, 4, 8),
  x_limits = c(0.15, 14)
)



# -------------------------------------------------------------------
#    Sensitivity analysis: Active or severe pigment on LBW risk
# -------------------------------------------------------------------
active_or_sev <- actsev %>%
  select(c(outcome, subgroup, 
           afe_rr, afe_rr_lb, afe_rr_ub, 
           pae_rr, pae_rr_lb, pae_rr_ub)) %>%
  rename(
    afe_rr_est = afe_rr,
    pae_rr_est = pae_rr
  ) %>%
  pivot_longer(
    cols = c(afe_rr_est, afe_rr_lb, afe_rr_ub, pae_rr_est, pae_rr_lb, pae_rr_ub),
    names_to = c("est_type", ".value"),
    names_pattern = "(.*_rr)_(est|lb|ub)"
  ) %>%  
  mutate(outcome = factor(outcome, levels=c('LBWdich', 'SGA', 'preterm')), 
         subgroup = factor(subgroup, levels=c( 'Multigravidae', 'Primigravidae', 'Overall')),
         est_type = factor(est_type, levels=c('afe_rr', 'pae_rr')),
         alpha = ifelse(est_type == "pae_rr", 0.2, 1))

actsev_LBW_p <- ggplot(data=active_or_sev, aes(x=log(est), y=outcome,
                                             col=subgroup, fill=subgroup,
                                             group=interaction(outcome, subgroup), 
                                             shape=est_type),) +
  geom_point(size=3.4,
             position=position_dodge(width=0.9)) +
  geom_pointrange(aes(xmin=log(lb), xmax=log(ub)), 
                  position=position_dodge(width=0.9), size=1) +
  facet_wrap(~est_type, ncol=1) +  
  theme_minimal() +
  geom_vline(xintercept=0, lty='dashed') +
  scale_color_manual(values=c('#729e9e', '#d08770', '#e9c46a')) +
  scale_fill_manual(values=c('white', 'white', 'white')) +
  scale_shape_manual(values=c(16, 21)) +   # pick two distinct shapes for est_type
  scale_x_continuous(breaks=log(c(0.25, 0.5, 1, 2.0, 4.0, 8)), 
                     labels=c(0.25, 0.5, 1.0, 2.0, 4.0, 8.0)) + 
  coord_cartesian(xlim=c(log(0.25), log(8.0))) +
  theme(axis.line.x = element_line(color='black'), 
        axis.ticks.x = element_line(color='black')) +
  ylab('') +
  xlab('Risk Ratio')

ggsave('Figures/Active or Severe on LBW.jpg', 
       actsev_LBW_p, 
       height=3, width=5, 
       dpi=300)



