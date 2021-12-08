library(tidyverse)
library(brms)
source("./process_tms_data.R")

#REMOVE OUTLIER
tms_processed <- tms_processed %>% filter(fraction < 5)

#Plot concentration against time
ggplot(data = tms_processed,
       aes(x=days, y=value, col=sample_id, group=sample_id)) +
  geom_line() + 
  theme_minimal()+
  theme(legend.position = "none") + 
  facet_wrap(vars(biomarker), scales="free")


#Plot concentration (as proportion of baseline) against time
ggplot(data = tms_processed,
       aes(x=days, y=100*(fraction-1))) +
  ylab("% Change")+
  geom_line(aes(col=sample_id, group=sample_id)) + 
  theme_minimal()+
  theme(legend.position = "none") + 
  facet_wrap(vars(biomarker), scales="free") + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), size = 1)

#Plot conditional effects with mean starting condition against time
#with boxplots
plot_list <- list()
for (bio in unique(tms_processed$biomarker)){
  df <- tms_processed %>% filter(biomarker == bio)
  
  baseline_median <- df %>%
    filter(timepoint == "T0") %>%
    select(value) %>%
    pull() %>%
    median()
  
  model <- readRDS(paste0("./models/", bio, "_SFO_model2.rds"))
  ce_object <- conditional_effects(
    model,
      effects = "years"
    )
  
  cond_effect_plot_data <- ce_object$years %>%
    mutate(across(effect1__:upper__, ~ .x * baseline_median))
  
  p <- ggplot(data = cond_effect_plot_data) +
    geom_line(aes(x = years,  y = estimate__)) + 
    geom_ribbon(
      aes(x = years  , ymin =  lower__, ymax = upper__), 
      alpha=0.3, fill="red") + 
    geom_point(data=df, 
      aes(x = years, y = value)) + 
    geom_boxplot(data = df,  
                 aes(y = value, group = years, x=years), alpha=0.2) + 
    theme_minimal() + 
    ylab("Concentration (pmol/L)") + xlab("Time (years)")
  plot_list[[bio]] <- p
} 

ggpubr::ggarrange(plot_list$Testosterone,
                  plot_list$Androstenedione,
                  plot_list$`11KT`,
                  plot_list$`11OHA4`,
                  plot_list$`17OHP`, labels = c("A", "B", "C", "D", "E"))


#################################################
#Plot conditional effects with mean starting condition against time
#with boxplots
#################
##BUT AS FRACTION OF BASELINE
###
plot_list <- list()
for (bio in unique(tms_processed$biomarker)){
  df <- tms_processed %>% filter(biomarker == bio)
  
  model <- readRDS(paste0("./models/", bio, "_SFO_model2.rds"))
  ce_object <- conditional_effects(
    model,
    effects = "years"
  )
  
  cond_effect_plot_data <- ce_object$years
  
  p <- ggplot(data = cond_effect_plot_data) +
    geom_line(aes(x = years,  y = estimate__)) + 
    geom_ribbon(
      aes(x = years  , ymin =  lower__, ymax = upper__), 
      alpha=0.3, fill="red") + 
    geom_point(data=df, 
               aes(x = years, y = fraction)) + 
    geom_boxplot(data = df,  
                 aes(y = fraction, group = years, x=years), alpha=0.2) + 
    theme_minimal() + 
    ylab("Concentration (pmol/L)") + xlab("Time (years)")
  plot_list[[bio]] <- p
} 

ggpubr::ggarrange(plot_list$Testosterone,
                  plot_list$Androstenedione,
                  plot_list$`11KT`,
                  plot_list$`11OHA4`,
                  plot_list$`17OHP`, labels = c("A", "B", "C", "D", "E"))
  
