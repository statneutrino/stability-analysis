library(tidyverse)
library(lubridate)
library(readxl)
#library(brms)
affinity <- read_xlsx("./data/affinity_final_stability_data.xlsx", sheet = "Stability study results") %>%
  as.data.frame()

affinity.biomarker.names <- c("Creatinine", "FSH", "LH", "Cotinine", "Cortisol", "DHEA", "OE2")

colnames(affinity)[10:17] <- c("Creatinine", "FSH", "LH", "Cotinine", "Cortisol", "DHEA", "OE2", "Date.Time")

##Process Affinity Data, steroids only
affinity_steroids <- affinity %>%
  select(UID, SampleType, Timepoint, AliquotNumber, OE2, DHEA, Date.Time) %>%
  filter(SampleType == "saliva") %>%
  select(-SampleType) %>%
  mutate(Date.Time = ymd_hms(Date.Time)) %>%
  mutate(days = interval(ymd("2020-02-21"), Date.Time) %/% days()) %>%
  mutate(years = round(days / 365.25,3)) %>%
  pivot_longer(
    cols=DHEA:OE2,
    names_to = "biomarker"
  ) %>%
  rename(sample_id = UID)

steroid_names <- c("DHEA", "OE2")

affinity_steroids$fraction <- affinity_steroids$value
for (subject in unique(affinity_steroids$sample_id)){
  subject.df <- affinity_steroids[affinity_steroids$sample_id == subject,]
  for(biomarker in steroid_names){
    biomarker_by_subject <- subject.df[subject.df$biomarker == biomarker,]
    baseline <- biomarker_by_subject[biomarker_by_subject$Timepoint == "T0","value"] %>% pull()
    affinity_steroids$fraction[affinity_steroids$sample_id == subject & affinity_steroids$biomarker == biomarker] <- biomarker_by_subject$value / baseline
  }
}

affinity_steroids <- affinity_steroids %>%
  mutate(sex = ifelse(substr(sample_id,1,1) == "M", "male", "female"))

##Process Affinity Data, URINE only
affinity_urine <- affinity %>%
  select(UID, SampleType, Timepoint, AliquotNumber, Creatinine, FSH, LH, Date.Time) %>%
  .[1:192,] %>%
  select(-SampleType) %>%
  mutate(Date.Time = ymd_hms(Date.Time)) %>%
  mutate(days = interval(ymd_hms("2020-02-24 11:30:00"), Date.Time) %/% days()) %>%
  mutate(years = round(days / 365.25,3)) %>%
  pivot_longer(
    cols=Creatinine:LH,
    names_to = "biomarker"
  ) %>%
  rename(sample_id = UID)

urine_names <- c("Creatinine", "FSH", "LH")

affinity_urine$fraction <- affinity_urine$value
for (subject in unique(affinity_urine$sample_id)){
  subject.df <- affinity_urine[affinity_urine$sample_id == subject,]
  for(biomarker in urine_names){
    biomarker_by_subject <- subject.df[subject.df$biomarker == biomarker,]
    baseline <- biomarker_by_subject[biomarker_by_subject$Timepoint == "T0","value"] %>% pull()
    affinity_urine$fraction[affinity_urine$sample_id == subject & affinity_urine$biomarker == biomarker] <- biomarker_by_subject$value / baseline
  }
}
         
