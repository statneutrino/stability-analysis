library(tidyverse)
library(readxl)

#Import stability results
tms.results <- read_xlsx("./data/TMS_saliva_stability_final.xlsx", sheet = "Sheet1")

head(tms.results)

#rename columns
colnames(tms.results)[3:7] <- colnames(tms.results)[2]
colnames(tms.results)[11:15] <- colnames(tms.results)[10]
colnames(tms.results)[19:23] <- colnames(tms.results)[18]
colnames(tms.results)[27:31] <- colnames(tms.results)[26]
colnames(tms.results)[35:39] <- colnames(tms.results)[34]

results.df <- tibble(sample_id="a", value=2, biomarker="c", timepoint="d")

for (column in 2:dim(tms.results)[2]){
  if(!is.na(tms.results[2,column])){
    x <-tibble(
      sample_id = tms.results[2:17,1] %>% pull(),
      value = tms.results[2:17,column] %>% pull(),
      biomarker = rep(colnames(tms.results)[column],16),
      timepoint = rep(tms.results[1,column]%>% pull(),16) 
    )
    results.df <- rbind(results.df, x)
  }
}
results.df <- results.df[-1,] %>%
  #Split analysis dates
  mutate(analysis.date = str_split(timepoint, " ", simplify = TRUE)[,2]) %>%
  mutate(analysis.date = gsub("\\(|\\)", "", analysis.date)) %>%
  #Create days since baseline variable
  mutate(analysis.date = as.Date(analysis.date, tryFormats = "%d/%m/%y")) %>%
  mutate(analysis.date = strptime(analysis.date, format="%Y-%m-%d", tz="GMT")) %>%
  mutate(days = difftime(analysis.date, strptime("2020-02-21", format="%Y-%m-%d", tz="UTC"), units = "days")) %>%
  mutate(days = as.numeric(days)) %>%
  #years
  mutate(years = days / 365.25) %>%
  #create timepoint and biomarker columns
  mutate(timepoint = str_split(timepoint, " ", simplify = TRUE)[,1]) %>%
  mutate(units = str_split(biomarker, "\\(", simplify = TRUE)[,2]) %>%
  mutate(units = gsub("\\(|\\)", "", units)) %>%
  mutate(biomarker = str_split(biomarker, "\\(", simplify = TRUE)[,1]) %>%
  mutate(across(where(is.character), str_trim)) %>%
  #change insufficient and mucinous to missing
  mutate(value = ifelse(grepl("Ins", .$value), NA, .$value)) %>%
  mutate(value = ifelse(grepl("Mucinous", .$value), NA, .$value)) %>%
  #REMOVE NAs from value
  filter(!is.na(value)) %>%
  mutate(value = as.numeric(value))

biomarker_names <- unique(results.df$biomarker)

results.df$fraction <- results.df$value
for (subject in unique(results.df$sample_id)){
  subject.df <- results.df[results.df$sample_id == subject,]

  for(biomarker in biomarker_names){
    biomarker_by_subject <- subject.df[subject.df$biomarker == biomarker,]
    baseline <- biomarker_by_subject[biomarker_by_subject$timepoint == "T0","value"] %>% pull()

    results.df$fraction[results.df$sample_id == subject & results.df$biomarker == biomarker] <- biomarker_by_subject$value / baseline
  }
}

write.csv(results.df, "./data/tms_processed.csv", row.names = FALSE)
rm(list = ls())
tms_processed <- read.csv("./data/tms_processed.csv")