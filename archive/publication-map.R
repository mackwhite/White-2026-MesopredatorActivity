### Project: Shark River Trophic Ecology
### Goal(s): build a map with google satellite imagery
### Author(s): Mack White
### Date(s): Summer 2024
### Notes:

# Housekeeping ------------------------------------------------------------
librarian::shelf(tidyverse, readr, janitor, zoo, scales, 
                 lubridate, openintro, maps, ggmap,
                 ggthemes, shapefiles, broom, sf, ggspatial, GISTools,
                 patchwork, magick)

### set theme ---
theme_set(theme_minimal())

### set up google map key ---
api_secret <- 'AIzaSyB_Q-Ow1klpX9jblm7M2614k5KCVYUXTZM'
register_google(key = api_secret)
has_google_key()

### read in necessary data ---
sites <- read_csv("local-data/vbf-accelerometer-snook.csv") |> 
      select(station, latitude, longitude) |> 
      distinct()

data <- read_rds('local-data/archive/accelerometer-model-data-012025.RDS') |> 
      group_by(station, latitude, longitude) |> 
      count()

dat1 <- sites |> left_join(data, by = c("station", "latitude", "longitude")) |> 
      mutate(scale = rescale(n, to = c(1,10)),
             data = "Acoustic") |> 
      select(-n) |> 
      rename(lat = latitude,
             long = longitude) |> 
      select(station, data, lat, long, scale)
glimpse(dat1)

monitoring <- read_csv('../snook-forage-quality/mapping/efishing_monitoringstation_coords.csv') |> 
      janitor::clean_names() |> 
      rename(station = site_name,
             data = type,
             lat = latitude,
             long = longitude) |> 
      mutate(scale = NA_real_) |> 
      select(station, data, lat, long, scale)
glimpse(monitoring)

dat <- rbind(dat1, monitoring) |> 
      filter(data != "Electrofishing" & data != "Nutrients")
glimpse(dat)

layer_three <- get_map(
      ### determine bounding box: https://www.openstreetmap.org/#map=5/25.304/-69.412
      c(left = -80.96, bottom = 25.39, right = -80.82, top = 25.51),
      maptype = 'satellite',
      source = 'google',
      api_key = api_secret
)

ggmap(layer_three) +
      geom_point(data = dat,
                 aes(x = long, y = lat),
                 color = "white",
                 size = 5.5) +
      geom_point(data = dat,
                 aes(x = long, y = lat, color = data),
                 size = 3.5) +
      scale_color_manual(values = c('Acoustic' = '#ffcc00',
                                    'Temperature' = '#e41a1c',
                                    'Marsh Stage' = '#377eb8')) +
      scale_y_continuous(limits =c(25.410,25.475)) +
      annotation_north_arrow(location = 'tl',
                             style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                    fill = 'white',
                                                                    line_col = 'white',
                                                                    text_face = "bold",
                                                                    text_size = 18)) +
      annotation_scale(location = 'br', width_hint = 0.5, text_cex = 1.25,
                       text_face = "bold", text_col = 'white') +
      coord_sf(crs = 4326) +
      theme_map() +
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(face = 'bold',
                                       color = 'white',
                                       size = 16))

### save for publication
ggsave('output/figs/figure-one-layer-one.png',
       dpi = 600, units= 'in', height = 6, width = 12)

layer_two <- get_map(
      ### determine bounding box: https://www.openstreetmap.org/#map=5/25.304/-69.412
      c(left = -81.5804, bottom = 24.9163, right = -80.4268, top = 25.8147),
      maptype = 'satellite',
      source = 'google',
      api_key = api_secret
)

ggmap(layer_two) +
      theme_map() +
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(face = 'bold',
                                       color = 'white',
                                       size = 16))

ggsave('output/figs/figure-one-layer-two.png',
       dpi = 600, units= 'in', height = 6, width = 12)

### attempt to combine everything here ---
# img1 <- image_read('output/figs/figure-one-layer-one.png')
# img2 <- image_read('output/figs/figure-one-layer-two.png')
# 
# combined_img <- image_composite(img1, img2, offset = "+500+500")
# image_write(combined_img, 'output/figs/combined-map.png')
