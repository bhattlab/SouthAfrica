# save all manuscript figures

library(here)

# main
figs <- list.files(here("scripts/figures"), pattern = "fig_.+.R")

for (f in figs){
  source(here("scripts/figures", f))
}

# supplementary
sfigs <- list.files(here("scripts/figures/supplementary"), pattern = "fig_.+.R")
sfigs <- sfigs[!grepl("S13", sfigs)]

for (f in sfigs){
  source(here("scripts/figures/supplementary", f))
}
