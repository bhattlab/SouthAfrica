## Plot a map of South Africa including Bushbuckridge and Soweto

library(cowplot)
library(ggrepel)
library(ggspatial)
library(here)
library(maps)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyverse)

# tutorial
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-2.html

# South Africa coordinate files
# https://gadm.org/download_country_v3.html

# world/country coordinates
world <- ne_countries(scale = "medium", returnclass = "sf")
sa <- ne_countries(scale = "medium", country = "south africa", returnclass = "sf")

# read SA province sf file
sa_sf <- readRDS(here("rds/gadm36_ZAF_1_sf.rds"))

# show Mpumalanga and Gauteng
sa_province_sf <- sa_sf[grepl("Mpumalanga|Gauteng", sa_sf$NAME_1), ]

# add color
sa_sf$fill_color <- ifelse(sa_sf$NAME_1 %in%
                             c("Gauteng", "Mpumalanga"), "fill", "no_fill")
sa_sf$fill_color <- factor(sa_sf$fill_color, levels = c("fill", "no_fill"))

site_symbols <- data.frame(
  site = c("Bushbuckridge", "Soweto"),
  long = c(31.3, 27.9),
  lat = c(-24.8, -26.2)
)

provinces <- data.frame(
  province = c("Mpumalanga", "Gauteng"),
  long = c(29.5, 26),
  lat = c(-24, -26.1)
)

countries <- data.frame(
  country = c("Namibia", "Botswana", "Zimbabwe", "Mozambique",
              "Lesotho", "Eswatini"),
  long = c(17.5, 24, 30.25, 33, 28.25, 31.45),
  lat = c(-27, -24.5, -21.5, -23, -29.5, -26.5)
  )

lines <- data.frame(x1 = c(29.6970, 26.6), x2 = c(30.7259, 27.2067),
                    y1 = c(-24.2516, -26.15), y2 = c(-24.5372, -26.3570))

ggplot(data = world) +
  geom_sf(fill = "snow") +
  geom_sf(data = sa, fill = "#eaf7e9") +
  geom_sf(data = sa_province_sf, fill = "#C9DBC7") +
  geom_text(data = countries, aes(x = long, y = lat, label = country),
            size = 3, fontface = "italic", color = "grey30") +
  geom_text(data = provinces, aes(x = long, y = lat, label = province),
            size = 4) +
  geom_point(data = site_symbols, aes(x = long, y = lat), size = 2,
             shape = 23, fill = "black") +
  geom_label_repel(data = site_symbols, aes(x = long, y = lat, label = site),
                   fontface = "bold", color = "black",
                   point.padding = unit(0.1, "cm"), ylim = c(-25, -27),
                   xlim = c(27, 32), label.size = 0) +
  annotate(geom = "text", x = 23, y = -30, label = "South Africa", size = 6,
           fontface = "italic", color = "grey30") +
  coord_sf(xlim = c(14.5, 34.5), ylim = c(-36, -21), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  theme_cowplot() +
  panel_border(color = "black") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_orienteering) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "lightblue1"),
    axis.line = element_blank()
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

ggsave(here("final_plots/supplementary/supp_info_map.png"),
       width = 8.64, height = 7.224)

