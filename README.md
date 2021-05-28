# Bushbuckridge and Soweto microbiome project

Analysis scripts for Wits-Stanford shotgun microbiome data.

To generate all figures from the manuscript, install all required packages and run each	R script in `scripts/figures` and `scripts/figures/supplementary`

e.g. from the command line:

     find scripts/figures -name "*.R" | xargs -I foo Rscript foo

Required packages:

* cowplot
* DESeq2
* dplyr
* genefilter
* ggplot2
* ggpubr
* ggrepel
* ggspatial
* gtools
* harrietr
* here
* maps
* MASS
* RColorBrewer
* reshape2
* rgeos
* rnaturalearth
* rnaturalearthdata
* scales
* sf
* tidyverse
* vegan
