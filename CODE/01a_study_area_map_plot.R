# ============================================================
# 01a_study_area_map.R
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Fig. 1A: study-area map of Papua New Guinea.
#   left  = overview (smaller), with a red box marking the detail extent
#   right = detail (larger) with the eight sampling localities
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(elevatr)
library(terra)
library(sf)
library(rnaturalearth)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(ggnewscale)
library(cowplot)
library(dplyr)
library(svglite)

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Area of interest: Papua New Guinea
# ------------------------------------------------------------
aoi <- ne_countries(country = "Papua New Guinea",
                    scale = "medium", returnclass = "sf") |>
  st_make_valid()

# ------------------------------------------------------------
# 2. Download DEM (elevation) and clip to the country
# ------------------------------------------------------------
dem <- get_elev_raster(locations = aoi, z = 7, clip = "locations") |>
  rast()
names(dem) <- "elevation"
dem[dem <= 0] <- NA

dem_clip <- mask(crop(dem, vect(aoi)), vect(aoi))

# ------------------------------------------------------------
# 3. Hillshade
# ------------------------------------------------------------
slope  <- terrain(dem_clip, "slope",  unit = "radians")
aspect <- terrain(dem_clip, "aspect", unit = "radians")
hsd    <- shade(slope, aspect, angle = 45, direction = 315)
names(hsd) <- "hillshade"

# ------------------------------------------------------------
# 4. Sampling localities
# ------------------------------------------------------------
sites <- data.frame(
  locality = c("Elem", "Morox", "Niksek", "Ohu", "Utai",
               "Wamangu", "Wanang", "Yapsiei"),
  Lat  = c(-4.8167, -4.0167, -4.7000, -5.2333, -3.3841,
           -3.7871, -5.2309, -4.6283),
  Long = c(143.9167, 144.1000, 142.5333, 145.6833, 141.5859,
           143.6521, 145.1818, 141.0973)
) |>
  st_as_sf(coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# ------------------------------------------------------------
# 5. Elevation palette (tropical greens -> highlands -> peaks)
# ------------------------------------------------------------
elev_palette <- c(
  "#2d6a4f", "#52b788", "#b7e4c7",
  "#e9c46a", "#d4a373",
  "#9c6644", "#7f5539",
  "#d6ccc2", "#ffffff"
)

# ------------------------------------------------------------
# 6. Detail map: extent from the localities + padding
# ------------------------------------------------------------
sites_bbox <- st_bbox(sites)
pad <- 0.6
xlim_detail <- c(as.numeric(sites_bbox["xmin"]) - pad,
                 as.numeric(sites_bbox["xmax"]) + pad)
ylim_detail <- c(as.numeric(sites_bbox["ymin"]) - pad,
                 as.numeric(sites_bbox["ymax"]) + pad)

detail_map <- ggplot() +
  geom_spatraster(data = hsd, alpha = 1, show.legend = FALSE) +
  scale_fill_gradient(low = "black", high = "white",
                      na.value = "transparent", guide = "none") +
  new_scale_fill() +
  geom_spatraster(data = dem_clip, alpha = 0.6) +
  scale_fill_gradientn(colours = elev_palette,
                       na.value = "transparent",
                       guide = "none") +
  geom_sf(data = aoi, fill = NA, color = "grey20", linewidth = 0.4) +
  geom_sf(data = sites, color = "black", fill = "red",
          shape = 21, size = 5, stroke = 0.6) +
  annotation_scale(location = "bl",
                   width_hint = 0.4,               # wider bar -> larger units
                   pad_x = unit(1, "cm"),
                   pad_y = unit(0.5, "cm"),
                   height = unit(0.25, "cm"),       # bar thickness
                   text_cex = 1.2) +                # label size
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = xlim_detail, ylim = ylim_detail, expand = FALSE) +
  labs(title = "Papua New Guinea - Sampling Localities") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_line(color = "white", linewidth = 0.2),
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# ------------------------------------------------------------
# 7. Overview map: red box marks the detail extent
# ------------------------------------------------------------
overview_countries <- ne_countries(
  country = c("Papua New Guinea", "Indonesia", "Australia"),
  scale = "medium", returnclass = "sf"
) |>
  st_make_valid()

# rectangle marking the detail area (same limits as the detail coord_sf)
detail_bbox <- st_as_sfc(st_bbox(c(xmin = xlim_detail[1],
                                   xmax = xlim_detail[2],
                                   ymin = ylim_detail[1],
                                   ymax = ylim_detail[2]),
                                 crs = st_crs(4326)))
overview_extent <- st_bbox(c(xmin = 129, xmax = 154,
                             ymin = -22, ymax = 1),
                           crs = st_crs(4326))

overview_map <- ggplot() +
  geom_sf(data = overview_countries,
          fill = "grey85", color = "grey40", linewidth = 0.2) +
  geom_sf(data = overview_countries[overview_countries$name %in%
                                      c("Papua New Guinea", "Indonesia"), ],
          fill = "grey70", color = "grey30", linewidth = 0.2) +
  geom_sf(data = detail_bbox, fill = NA, color = "red", linewidth = 0.8) +
  coord_sf(xlim = c(overview_extent["xmin"], overview_extent["xmax"]),
           ylim = c(overview_extent["ymin"], overview_extent["ymax"]),
           expand = FALSE) +
  labs(title = "") +
  theme_minimal(base_size = 10) +
  theme(panel.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.4),
        panel.grid = element_line(color = "grey92", linewidth = 0.2),
        plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 7))

# ------------------------------------------------------------
# 8. Combine: overview (left, smaller) + detail (right, larger)
# ------------------------------------------------------------
final_map <- plot_grid(
  overview_map, detail_map,
  ncol = 2,
  rel_widths = c(1, 2.2),     # right panel ~2x wider
  align = "h",
  axis = "tb"
)

# ------------------------------------------------------------
# 9. Export (combined + detail alone)
# ------------------------------------------------------------
ggsave("output/fig/Figure_1A_localities_map.svg", final_map,
       width = 14, height = 8, dpi = 300, device = svglite::svglite)
ggsave("output/fig/Figure_1A_localities_map.pdf", final_map,
       width = 14, height = 8, dpi = 300)
ggsave("output/fig/Figure_1A_localities_map.png", final_map,
       width = 14, height = 8, dpi = 300, bg = "white")

detail_width  <- 10
detail_height <- 8
ggsave("output/fig/Figure_1A_detail_map.svg", detail_map,
       width = detail_width, height = detail_height, dpi = 300, device = svglite::svglite)
ggsave("output/fig/Figure_1A_detail_map.pdf", detail_map,
       width = detail_width, height = detail_height, dpi = 300)
ggsave("output/fig/Figure_1A_detail_map.png", detail_map,
       width = detail_width, height = detail_height, dpi = 300, bg = "white")