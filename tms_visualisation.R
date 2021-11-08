library(tidyverse)
source("stability_analysis_saliva/process_tms_data.R")

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


  
