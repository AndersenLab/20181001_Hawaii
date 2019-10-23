library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# count number of collections = 1001 
collection_count <- df %>%
  dplyr::distinct(fulcrum_id_collection)

# count number of caenorhabditis elegans collected 
isolated_c_worms <- df %>%
  dplyr::filter

# from genotyping spreadsheet 41 C. elegans from 1001 collections 2018 trip
hit_rate_2018 <- 41/1001

# from 2017 trip 38 c. elgans from 2263 collections
hit_rate_2017 <-  38/2263

# 2018 isolates = 2026
# 2017 isolates = 2532
# total isolates = 4558

# 
test <- df %>%
  dplyr::filter(confirmed_species_id == "C. elegans", males_observed == "yes")

