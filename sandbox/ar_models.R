library(tidyverse)
library(brms)

source("./process_tms_data.R")

#REMOVE OUTLIER
tms_processed <- tms_processed %>% filter(fraction < 5)
tms_processed$conc_nmol <- tms_processed$value / 1000

bio <- "Testosterone"

#For each biomarker, fit
#1 linear model fixed k
#2 linear model k varies between-subject
#3 SFO fixed k
#4 SFO k varies between-subject

biomarkers <- unique(tms_processed$biomarker)

SFO_k_fixed_formula <- brmsformula(
  fraction ~ C0 * exp(-k * years),
  C0~1+(1|sample_id), 
  k~1,
  nl=TRUE)

SFO_k_random_formula <- brmsformula(
  fraction ~ C0 * exp(-k * years),
  C0~1+(1|sample_id), 
  k~1+(1|sample_id),
  nl=TRUE)

df <- tms_processed %>% filter(biomarker == bio)
  
#SET PRIORS FOR linear models
location_b <- mean(df$fraction) %>% signif(., 2)
scale_b <- sd(df$fraction)
  
b_class_prior <- paste0("normal(0, 1)")
intercept_prior <- paste0("normal(1, 0.02)")
sigma_prior <- paste0("student_t(4, 0, 0.2)")
group_variance_prior <- paste0("student_t(4, 0, 0.1)")

prior_linear_mod2 <- c(
    prior_string(b_class_prior, class="b"),
    prior_string(intercept_prior, class="Intercept"),
    prior_string(group_variance_prior, class="sd"),
    prior_string(sigma_prior, class="sigma"),
    prior(lkj(2), class="cor")
  )

##Fit linear model with ar
# linear_ar <- brm(
#     formula = fraction ~ years + (years|sample_id) + ar(time = years, gr = sample_id, p = 1, cov = FALSE),
#     data = df,
#     prior = prior_linear_mod2,
#     family = gaussian(),
#     chains = 6,
#     cores = 12,
#     iter = 3000,
#     warmup = 1000,
#     #backend = "cmdstanr",
#     save_pars = save_pars(all = TRUE)
# )
# linear_ar <- add_criterion(linear_ar, "loo", moment_match=TRUE)
# saveRDS(linear_ar, paste0("./sandbox/", bio, "_linear_ar.rds"))

##Fit sfo_model with ar
SFO_k_random_formula <- brmsformula(
    fraction ~ C0 * exp(-k * years) + ar(time = years, gr = sample_id, p = 1, cov = FALSE),
    C0~1+(1|sample_id), 
    k~1+(1|sample_id),
    nl=TRUE)

prior_SFO2 <- c(
    prior_string(intercept_prior, class="b", nlpar="C0"),
    prior_string("normal(0, 1)", lb=0, nlpar="k"),
    prior_string(group_variance_prior, class="sd", nlpar="C0"),
    prior_string(group_variance_prior, class="sd", nlpar="k"),
    prior_string(sigma_prior, class="sigma")
  )

SFO_ar <- brm(
    formula = SFO_k_random_formula,
    #sample_prior = TRUE,
    data = df,
    family=gaussian(),
    prior = prior_SFO2,
    chains = 6,
    cores = 12,
    iter = 3000,
    warmup = 1000,
    #backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  SFO_ar <- add_criterion(SFO_ar, "loo", moment_match=TRUE)
  saveRDS(SFO_ar, paste0("./sandbox/", bio, "_SFO_ar.rds"))


