# Bushbuckridge and Soweto microbiome project

Analysis scripts for Wits-Stanford shotgun microbiome data.

To generate all figures from the manuscript, install all required packages and run each	R script in `scripts/figures` and `scripts/figures/supplementary`

e.g. from the command line:

     find scripts/figures -name "*.R" | xargs -I foo Rscript foo

Required packages:

* DESeq2
* MASS
* Maaslin2
* RColorBrewer
* cowplot
* dplyr
* genefilter
* ggplot2
* ggpubr
* ggrepel
* ggspatial
* gplots
* gtools
* harrietr
* here
* maps
* metagenomeSeq
* reshape2
* rgeos
* rnaturalearth
* rnaturalearthdata
* scales
* sf
* tidyverse
* vegan
