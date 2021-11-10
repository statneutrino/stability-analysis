library(tidyverse)
library(brms)

source("./process_tms_data.R")

#REMOVE OUTLIER
tms_processed <- tms_processed %>% filter(fraction < 5)
tms_processed$conc_nmol <- tms_processed$value / 1000

bio <- "11KT"

df <- tms_processed %>% filter(biomarker == bio)

models <- paste0("./models/",list.files(
  path = "./models",
  pattern = bio
))

conditions <- data.frame(
    zAge = c(-1, 0, 1)
  )

linear_mod1 <- readRDS(models[1])
linear_mod2 <- readRDS(models[2])
sfo_mod1 <- readRDS(models[3])
sfo_mod2 <- readRDS(models[4])

x1 <- conditional_effects(linear_mod1)$years %>% data.frame()
x2 <- conditional_effects(linear_mod2)$years %>% data.frame()
x3 <- conditional_effects(sfo_mod1)$years %>% data.frame()
x4 <- conditional_effects(sfo_mod2)$years %>% data.frame()

cond_effect_plot_data <- rbind(
  data.frame(
    years = x1$years,
    fraction = x1$fraction,
    x = x1$effect1__,
    est = x1$estimate__,
    lower = x1$lower__,
    upper = x1$upper__,
    model = rep("linear_fixed", 100)
  ),
  data.frame(
    years = x2$years,
    fraction = x2$fraction,
    x = x2$effect1__,
    est = x2$estimate__,
    lower = x2$lower__,
    upper = x2$upper__,
    model = rep("linear_random", 100)
  ),
  data.frame(
    years = x3$years,
    fraction = x3$fraction,
    x = x3$effect1__,
    est = x3$estimate__,
    lower = x3$lower__,
    upper = x3$upper__,
    model = rep("sfo_fixed", 100)
  ),
  data.frame(
    years = x4$years,
    fraction = x4$fraction,
    x = x4$effect1__,
    est = x4$estimate__,
    lower = x4$lower__,
    upper = x4$upper__,
    model = rep("sfo_random", 100)
  )
)

ggplot(data = cond_effect_plot_data %>% filter(model %in% c("linear_random", "sfo_random"))) + 
  geom_line(aes(x = x, y = est, col=model, group=model)) + 
  geom_ribbon(
    aes(x = x, ymin = lower, ymax = upper, fill=model, group=model), 
    alpha=0.2) + 
  geom_point(data=df, aes(x = years, y = fraction)) + 
  theme_minimal() 


conditions <- data.frame(sample_id = unique(df$sample_id))
rownames(conditions) <- unique(df$sample_id)
me_loss <- conditional_effects(
  IORE_fixed, conditions = conditions, 
  re_formula = NULL, method = "predict"
)
plot(me_loss, ncol = 5, points = TRUE)