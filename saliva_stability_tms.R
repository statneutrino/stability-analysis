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


for (bio in biomarkers){
  
  print("Calculating stability for")
  print(bio)
  
  df <- tms_processed %>% filter(biomarker == bio)
  
  #SET PRIORS FOR linear models
  location_b <- mean(df$fraction) %>% signif(., 2)
  scale_b <- sd(df$fraction)
  
  b_class_prior <- paste0("normal(0, 1)")
  intercept_prior <- paste0("normal(1, 0.02)")
  sigma_prior <- paste0("student_t(4, 0, 0.2)")
  group_variance_prior <- paste0("student_t(4, 0, 0.1)")
  
  prior_linear_mod1 <- c(
    prior_string(b_class_prior, class="b"),
    prior_string(intercept_prior, class="Intercept"),
    prior_string(group_variance_prior, class="sd"),
    prior_string(sigma_prior, class="sigma")
  )
  
  print("Fitting linear model (fixed k)")
  linear_mod1 <- brm(
    formula = fraction ~ years + (1|sample_id),
    data = df,
    prior = prior_linear_mod1,
    sample_prior = TRUE,
    family = gaussian(),
    chains = 6,
    cores = 12,
    iter = 3000,
    warmup = 1000,
    #backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  linear_mod1 <- add_criterion(linear_mod1, "loo", moment_match = TRUE)
  saveRDS(linear_mod1, paste0("stability_analysis_saliva/", bio, "_linear_mod1.rds"))

  Sys.sleep(10)
  print(summary(linear_mod1))

  print("Fitting linear model (k varies between subject)")
  prior_linear_mod2 <- c(
    prior_string(b_class_prior, class="b"),
    prior_string(intercept_prior, class="Intercept"),
    prior_string(group_variance_prior, class="sd"),
    prior_string(sigma_prior, class="sigma"),
    prior(lkj(2), class="cor")
  )
  
  
  linear_mod2 <- brm(
    formula = fraction ~ years + (years|sample_id),
    data = df,
    prior = prior_linear_mod2,
    sample_prior = TRUE,
    family = gaussian(),
    chains = 6,
    cores = 12,
    iter = 3000,
    warmup = 1000,
    #backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  linear_mod2 <- add_criterion(linear_mod2, "loo", moment_match=TRUE)
  saveRDS(linear_mod2, paste0("stability_analysis_saliva/", bio, "_linear_mod2.rds"))
  
  Sys.sleep(10)
  
  print(summary(linear_mod2))
  #Set priors
  prior_SFO1 <- c(
    prior_string(intercept_prior, class="b", nlpar="C0"),
    prior_string("normal(0, 1)", lb=0, nlpar="k"),
    prior_string(group_variance_prior, class="sd", nlpar="C0"),
    prior_string(sigma_prior, class="sigma")
  )
  print("Fitting SFO model (fixed k)")
  SFO_model1 <- brm(
    formula = SFO_k_fixed_formula,
    sample_prior = TRUE,
    data = df,
    family=gaussian(),
    prior = prior_SFO1,
    chains = 6,
    cores = 12,
    iter = 3000,
    warmup = 1000,
    #backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  SFO_model1 <- add_criterion(SFO_model1, "loo", moment_match=TRUE)
  saveRDS(SFO_model1, paste0("stability_analysis_saliva/", bio, "_SFO_model1.rds"))
  
  Sys.sleep(10)
  
  print(summary(SFO_model1))
  
  prior_SFO2 <- c(
    prior_string(intercept_prior, class="b", nlpar="C0"),
    prior_string("normal(0, 1)", lb=0, nlpar="k"),
    prior_string(group_variance_prior, class="sd", nlpar="C0"),
    prior_string(group_variance_prior, class="sd", nlpar="k"),
    prior_string(sigma_prior, class="sigma")
  )
  

  print("Fitting SFO model (k varies between subject)")
  SFO_model2 <- brm(
    formula = SFO_k_random_formula,
    sample_prior = TRUE,
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
  SFO_model2 <- add_criterion(SFO_model2, "loo", moment_match=TRUE)
  saveRDS(
    SFO_model2,
    paste0("stability_analysis_saliva/", bio, "_SFO_model2.rds")
    )

  Sys.sleep(10)

  print(summary(SFO_model2))

  print("comparing the LOO for models for")
  print(bio)
  print(loo_compare(linear_mod1, linear_mod2, SFO_model1, SFO_model2))
  print("completed analysis for")
  print(bio)
}


###Double FIRST ORDER IN PARALLEL
DFOP_k_fixed_formula <- brmsformula(
  fraction ~ C0 * (g * exp(-k1 * years) + (1 - g) * exp(-k2 * years)),
  C0~1, 
  k1~1,
  k2~1,
  g~1,
  nl=TRUE)

prior_DFOP_fixed <- c(
  set_prior("normal(1, 0.2)", nlpar="C0"),
  set_prior("normal(0, 1)", lb=0, nlpar="k1"),
  set_prior("normal(0, 1)", lb=0, nlpar="k2"),
  set_prior("beta(0.75,2)", ub=0.5, nlpar="g"),
  prior_string(sigma_prior, class="sigma")
)

DFOP_fixed <- brm(
  formula = DFOP_k_fixed_formula,
  sample_prior = TRUE,
  data = df,
  family=gaussian(),
  prior = prior_DFOP_fixed,
  chains = 3,
  cores = 4,
  iter = 3000,
  warmup = 1000,
  #backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)


###POTENTIAL IORE model
IORE_k_fixed_formula <- brmsformula(
  fraction ~ (C0^(1-N) - (1-N) * k * years)^(1/(1-N)),
  C0~1,
  k~1,
  N~1,
  nl=TRUE)

prior_IORE1 <- c(
  set_prior("normal(1, 0.2)", nlpar="C0"),
  set_prior("normal(0, 1)", lb=0, nlpar="k"),
  set_prior("uniform(1.0001, 5)", nlpar="N")
)

prior_IORE1 <- c(
  prior_string(intercept_prior, class="b", nlpar="C0"),
  prior_string("normal(0, 1)", lb=0, nlpar="k"),
  set_prior("normal(1, 1)", lb=1, nlpar="N"),
  prior_string(sigma_prior, class="sigma")
)


IORE_fixed <- brm(
  formula = IORE_k_fixed_formula,
  sample_prior = TRUE,
  data = df,
  family=gaussian(),
  prior = prior_IORE1,
  chains = 3,
  cores = 4,
  iter = 3000,
  warmup = 1000,
  #backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

conditional_effects(IORE_fixed)

