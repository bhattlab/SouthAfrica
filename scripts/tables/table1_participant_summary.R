## participant summary table
library(here)
library(reshape2)
library(tidyverse)

load(here("RData/metadata.RData"))

summary_tbl <- za_pheno %>%
  mutate(age_mb = age) %>%
  group_by(sample) %>%
  mutate(
    bp_sys = mean(c(bp_sys_1, bp_sys_2, bp_sys_3), na.rm = T),
    bp_dia = mean(c(bp_dia_1, bp_dia_2, bp_dia_3), na.rm = T),
  ) %>%
  group_by(site) %>%
  summarise(
    n_participants = n(),
    age_mean = mean(age_mb),
    age_sd = sd(age_mb),
    age_low = min(age_mb),
    age_high = max(age_mb),
    bmi_mean = mean(bmi, na.rm = T),
    bmi_sd = sd(bmi, na.rm = T),
    bmi_low = min(bmi, na.rm = T),
    bmi_high = max(bmi, na.rm = T),
    bp.sys_mean = mean(bp_sys, na.rm = T),
    bp.sys_sd = sd(bp_sys, na.rm = T),
    bp.sys_low = min(bp_sys, na.rm = T),
    bp.sys_high = max(bp_sys, na.rm = T),
    bp.dia_mean = mean(bp_dia, na.rm = T),
    bp.dia_sd = sd(bp_dia, na.rm = T),
    bp.dia_low = min(bp_dia, na.rm = T),
    bp.dia_high = max(bp_dia, na.rm = T)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

s <- melt(summary_tbl, id.vars = c("site", "n_participants"))
s <- s %>%
  tidyr::separate(variable, into = c("var", "stat"), sep = "_")

s <- dcast(s, site + n_participants + var ~ stat, value.var = "value")
s <- s %>%
  group_by(site, var) %>%
  mutate(range = paste0(low, " - ", high)) %>%
  select(-low, -high, -n_participants) %>%
  select(var, site, mean, sd, range) %>%
  arrange(var)

write.table(s, here("final_tables/table_1_participants.txt"), sep = "\t",
            quote = F, row.names = F)