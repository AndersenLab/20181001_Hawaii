library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)
library(janitor)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# color palettes from pals package: https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html
species_palette <- c("C. elegans" = "#BE0032", #7
                     #"C. oiwi" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     #"Panagrolaimus sp." = "#C2B280", #8
                     #"Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1")  #6
#"multiple" = "#F3C300")

island_palette <- c("Kauai" = "#E69F00",
                    "Oahu" = "#56B4E9",
                    "Molokai" = "#009E73",
                    "Maui" = "#F0E442",
                    "Big Island" = "#D55E00")

####################################################
#  A: Define Functions                             #
####################################################
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

islands = list(
  "Oahu" = c(-158.323116,21.112767,-157.623081,21.814254),
  "Big Island" = c(-156.3651,18.8049,-154.765,20.4064)
)

gtmap <- function(loc) {
  get_map(location = loc,
          maptype = "terrain-background",
          source = "stamen",
          scale = "auto")
}

mget_map <- memoise(gtmap)

map_overview <- function(df, geoms, face = "plain") {
  island_set = lapply(names(islands), function(i) {
    
    if (i == "ALL") {
      l_position = c(0.65, 0.8)
      imap <- df
      island_size = 0.8
      rects = sapply(2:6, function(x) {
        map_label = LETTERS[x-1]
        l = islands[[x]]
        c(annotate("rect", xmin = l[1], xmax = l[3], ymin=l[2], ymax=l[4], alpha = 0.0, color = "black"),
          annotate("text", x = l[3] + 0.34, y = l[4] + 0.34, label = map_label),
          annotate("segment", x = l[3] + 0.18, xend = l[3], y = l[4] + 0.18, yend = l[4], color = "black"))
      })
      #No stroke on mini-map
    } else {
      l_position = "none"
      island_size = 2
      imap <- df %>% dplyr::filter(island == i)
      rects = element_blank()
    }
    
    map = mget_map(islands[[i]])
    # Calculate scalebar
    bb <- attr(map,"bb")
    sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                       lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                       lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                       lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))
    
    sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                      c(sbar$lon.end,sbar$lat.end))
    
    scalebar.length <- 20
    sbar$lon.end <- sbar$lon.start +
      ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000
    ptspermm <- 2.83464567
    
    base_map <- ggplot(imap) +
      ggmap::inset_ggmap(map) +
      #rects +
      geoms +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(c(0,0,0,0), "lines"),
            axis.line = element_blank(),
            plot.title = element_text(lineheight=.8, face="bold", vjust=1),
            plot.margin = unit(c(0,0,0,0), "lines"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = l_position,
            legend.background = element_rect(fill="white"),
            legend.text=element_text(size=12, color = "black", face = face)) +
      coord_equal(ratio=1) +
      scale_x_continuous(limits = islands[[i]][c(1,3)], expand = c(0, 0)) +
      scale_y_continuous(limits = islands[[i]][c(2,4)], expand = c(0, 0)) +
      geom_segment(data = sbar,
                   aes(x = lon.start,
                       xend = lon.end,
                       y = lat.start,
                       yend = lat.end),
                   arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                               ends = "both", type = "open")) +
      geom_text(data = sbar,
                aes(x = (lon.start + lon.end)/2,
                    y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
                    label = paste(format(scalebar.length),
                                  'km')),
                hjust = 0.5,
                vjust = 0,
                size = 8/ptspermm)  +
      coord_map(projection = "mercator",
                xlim=c(bb$ll.lon, bb$ur.lon),
                ylim=c(bb$ll.lat, bb$ur.lat)) +
      scale_radius(range = c(island_size, island_size), guide = "none") #+
    #scale_shape_manual(values = shape)
    
    base_map
    
  })
  
  
  plot_grid(plotlist = island_set,
            labels = c("A - Oahu",
                       "B - Big Island"
            ),
            label_y = 0.98,
            hjust = 0,
            label_x = 0.06,
            align = "vh")
}

####################################################
#           Make overview plot                    # 
####################################################
# setup overview plot groups
plot_spp_ids <- df %>%
  dplyr::filter(confirmed_species_id %in% c("C. elegans", "C. tropicalis", "C. briggsae")) %>%
  #dplyr::filter(!is.na(collection_longitude_photo)) %>%
  dplyr::mutate(collection_type = confirmed_species_id) %>%
  dplyr::select(c_label_corrected, collection_type, island, confirmed_species_id, worms_on_sample, collection_longitude_photo, collection_latitude_photo) %>%
  dplyr::distinct(c_label_corrected, collection_type, .keep_all=T) %>% 
  dplyr::group_by(c_label_corrected) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(collection_type = forcats::as_factor(collection_type),
                collection_type = forcats::fct_relevel(collection_type,
                                                       "C. elegans",
                                                       #"C. oiwi",
                                                       "C. tropicalis",
                                                       "C. briggsae")) %>%
  dplyr::arrange(collection_type) %>% # arrange sets order for c-labels with multiples so highest priority collection type is on top
  dplyr::distinct(c_label_corrected, .keep_all = T) %>% # selects highest priority collection type from a c-label with multiple collection types on it
  dplyr::arrange(desc(collection_type)) # reorders c-labels so highest priority collections are plotted on top


# Make map overview plot
plot_spp_ids_p <- map_overview(plot_spp_ids,
                               c(geom_point(aes(x=collection_longitude_photo,
                                                y=collection_latitude_photo,
                                                fill=collection_type,
                                                size = 1),
                                            color="black",
                                            shape=21,
                                            stroke = 0.5
                               ),
                               scale_fill_manual("species", values = species_palette)
                               ),
                               face="italic"
)

# make nice legned
legend_species <- get_legend(ggplot(plot_spp_ids %>% dplyr::filter(collection_type %in% names(species_palette)) %>%
                                      dplyr::mutate(collection_type = factor(collection_type, levels = names(species_palette)))) +
                               geom_bar(aes(x = c_label_corrected, fill = collection_type)) +
                               scale_fill_manual(values=c(species_palette)) +
                               labs(fill = "species") +
                               theme(legend.position = "bottom"))


collection_overview <- cowplot::plot_grid(plot_spp_ids_p, legend_species,
                                 rel_heights = c(10, 1), nrow=2)

ggsave('plots/Figure1_Collection_Overview_caenorhabditis.pdf', width = 12, height = 8.6)
ggsave('plots/Figure1_Collection_Overview_caenorhabditis.png', width = 12, height = 8.6)

# Export dataframe with gps cordinates of collections and one with totals
count_dist_collections_sp <- df %>%
  dplyr::filter(confirmed_species_id %in% c("C. elegans", "C. tropicalis", "C. briggsae")) %>%
  dplyr::distinct(c_label_corrected, confirmed_species_id, .keep_all =T) %>%
  dplyr::group_by(confirmed_species_id) %>%
  dplyr::mutate(species_distinct_collections = n()) %>%
  dplyr::distinct(confirmed_species_id, .keep_all = T) %>%
  dplyr::select(confirmed_species_id, species_distinct_collections)

count_number_of_strains_per_species <- df %>%
  dplyr::filter(confirmed_species_id %in% c("C. elegans", "C. tropicalis", "C. briggsae")) %>%
  dplyr::group_by(confirmed_species_id) %>%
  dplyr::mutate(number_of_strains_per_species = n()) %>%
  dplyr::distinct(confirmed_species_id, .keep_all = T) %>%
  dplyr::select(confirmed_species_id, number_of_strains_per_species)

output_summary1 <- full_join(plot_spp_ids, count_dist_collections_sp)
output_summary2 <- full_join(output_summary1, count_number_of_strains_per_species)

output_summary2 %>%
  dplyr::filter(confirmed_species_id %in% c("C. elegans", "C. tropicalis", "C. briggsae")) %>%
  readr::write_tsv("data/collection_summary_caenorhabditis.tsv")