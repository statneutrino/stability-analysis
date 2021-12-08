library(tidyverse)
library(brms)

source("./process_tms_data.R")

#REMOVE OUTLIER
tms_processed <- tms_processed %>% filter(fraction < 5)
tms_processed$conc_nmol <- tms_processed$value / 1000

bio <- "11KT"

df <- tms_processed %>% filter(biomarker == bio)


fm1 <- lme4::nlmer(fraction ~ SSbiexp(years, A0, k, 0, 0) ~ A0|sample_id,
            data = df,
            start = c(k=0.2, A0 = 1))

nform <- ~A0 * exp(-k * years)
nfun <- deriv(nform,namevec=c("A0","k"),
              function.arg=c("years","A0","k"))

mod2 <- lme4::nlmer(fraction ~ nfun(years, k, A0) ~ k|sample_id,
                   data = df,
                   start = c(k=0.2, A0 = 1))


library(rstanarm)

fm1 <- stan_nlmer(fraction ~ SSbiexp(years, A0, k1, A1, k2) ~ A0|sample_id,
                   data = df,
                   chains = 3)


rstan::get_stanmodel(fm1$stanfit)

yrep <- posterior_predict(fm1, draws = 50)
color_scheme_set("brightblue")
bayesplot::ppc_dens_overlay(df$fraction, yrep)
