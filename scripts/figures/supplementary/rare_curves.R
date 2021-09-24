library(here)
library(tidyverse)
library(vegan)

## load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))

## rarefaction curves ----

set.seed(1)

cutoff <- 100

cts <- za_S

cts[cts < cutoff] <- 0
cts <- cts[which(rowSums(cts) > 0), ]

cts_t_rare <- rrarefy(t(cts), 2e5)

rarecurve(cts_t_rare, step = 1e3)
