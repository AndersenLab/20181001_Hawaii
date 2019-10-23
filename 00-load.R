library(tidyverse)
library(lubridate)
library(geosphere)
library(googlesheets)
library(googledrive)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#====================================#
# Part 1: functions and variables    #
#====================================#

# set project code(s)
projects <- "HawaiiOctober2018"

LL <- c("longitude", "latitude")
LLA <- c("longitude", "latitude", "altitude")

filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

FtoC <- function(F) {
  (F - 32)*(5/9)
}

#=========================================================#
# Part 2: (Optional) Read in isolation correction list    #
#=========================================================#
#manually build isolation correction files based on outside script.
isolation_correction <- readr::read_csv("data/isolation_correction_pairs.csv")

#==============================================================#
# Part 3: Read in nematode field sampling data (collection)    #
#==============================================================#
collection <- readr::read_csv("data/fulcrum/nematode_field_sampling.csv") %>%
  dplyr::mutate(c_label = stringr::str_to_upper(c_label)) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(duplicated_c_label = ifelse(n()>1, 1, 0)) %>%
  dplyr::mutate(duplicated_c_label_corrected = ifelse(duplicated_c_label == 1 & c_label %in% isolation_correction$c_label, 1,
                                                      ifelse(duplicated_c_label == 1 & !(c_label %in% isolation_correction$c_label), 0, NA))) %>%
  dplyr::ungroup() %>%
  tidyr::separate(sample_photo, into = c("a", "b", "sample_photo"), sep = ",", extra = "merge", fill = "left", remove = TRUE) %>%
  dplyr::rename(collected_by = created_by,
                collection_time = time,
                collection_latitude = latitude,
                collection_longitude = longitude) %>%
  dplyr::select(-updated_at,
                -system_created_at,
                -system_updated_at,
                -date,
                -a,
                -b,
                -sample_photo_caption) %>%
  dplyr::mutate(collection_datetime = lubridate::ymd_hms(created_at)) %>%
  dplyr::mutate(collection_date = lubridate::date(created_at)) %>%
  dplyr::select(-created_at) %>%
  # Label substrate moisture as NA for substrates that did not have reading (we entered 100 for no reading) 
  dplyr::mutate(substrate_moisture = ifelse(substrate_moisture == 100, NA, substrate_moisture)) %>%
  # Fix ambient temp F to C
  dplyr::mutate(ambient_temperature_c = ifelse(ambient_temperature_c > 50,
                                             FtoC(ambient_temperature_c),
                                             ambient_temperature_c))

#========================================================================================================#
# Part 4: Add data from collection photos lat, long, elevation. Only need to run once to extract data    #
#========================================================================================================#
# Read in data from photos. Need to install using ‘brew install exiftool’ in terminal.
# comm <- paste0("exiftool -coordFormat '%+.6f' -csv -ext jpg ",
#                getwd(),
#                "/photos/*")
# 
# # Exif Data
# exif <- readr::read_csv(pipe(comm)) %>%
#    dplyr::mutate(SourceFile = stringr::str_replace(basename(SourceFile), ".jpg", "")) %>%
#    dplyr::select(sample_photo = SourceFile,
#                 altitude = GPSAltitude,
#                 latitude = GPSLatitude,
#                 longitude = GPSLongitude,
#                 ExposureTime,
#                Artist,
#                 Aperture,
#                 BrightnessValue,
#                 PhotoDate = DateCreated,
#                 FOV) %>%
#    dplyr::mutate(altitude =  as.numeric(stringr::str_replace(altitude, " m", ""))) %>%
#    dplyr::mutate(FOV =  as.numeric(stringr::str_replace(FOV, " deg", ""))) %>%
#    dplyr::group_by(sample_photo) %>%
#    # Only retain data from one sample photo.
#    dplyr::distinct(.keep_all=T)
#  save(file = "data/exif.Rda", exif)
load("data/exif.Rda")


collection <- dplyr::full_join(collection, exif) %>%
  dplyr::rename(collection_latitude_photo = latitude,
                collection_longitude_photo = longitude,
                collection_elevation = altitude) %>%
  dplyr::select(-Artist, -ExposureTime, -Aperture, -BrightnessValue,
                -PhotoDate, -FOV) %>%
  #remove extra fulcrum photos
  dplyr::filter(!is.na(c_label)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gps_err = geosphere::distHaversine(c(collection_longitude, collection_latitude),
                                                   c(collection_longitude_photo, collection_latitude_photo))) %>%
  dplyr::ungroup()

#========================================================#
# Part 4: Read in nematode isolation data (isolation)    #
#========================================================#
isolation <- readr::read_csv("data/fulcrum/nematode_isolation.csv") %>%
  dplyr::mutate(date = lubridate::mdy(date)) %>%
  # Need to filter isolation data from wrong time period
  dplyr::filter(project == "HawaiiOctober2018" & date > ymd("2018-11-1") | is.na(project) & date > ymd("2018-11-1")) %>%
  dplyr::select(fulcrum_id,
              c_label_id = c_label,
              isolated_at = system_created_at,
              isolated_by = created_by,
              worms_on_sample,
              approximate_number_of_worms,
              males_observed,
              dauers_on_sample,
              approximate_number_of_worms,
              isolation_date = date,
              isolation_time = time,
              isolation_latitude = latitude,
              isolation_longitude = longitude)

# Add isolations that were accidentically assigned to other projects and were C. elegans, C. briggsae, or C. tropicalis positive.
isolation_cross_project_corrections <- readr::read_csv("data/nematode_isolation_cross_project_corrections.csv") %>%
  dplyr::mutate(date = lubridate::mdy(date)) %>%
  dplyr::select(fulcrum_id,
                c_label_id = c_label,
                isolated_at = system_created_at,
                isolated_by = created_by,
                worms_on_sample,
                approximate_number_of_worms,
                males_observed,
                dauers_on_sample,
                approximate_number_of_worms,
                isolation_date = date,
                isolation_time = time,
                isolation_latitude = latitude,
                isolation_longitude = longitude)

# Join isolation dataframes together
isolation <- full_join(isolation, isolation_cross_project_corrections)
  

#================================================================================================#
# Part 5: Read in data from isolation correction file and correct isolation data if necessary    #
#================================================================================================#
# apply corrections to isolation data
isolation <- left_join(isolation, isolation_correction, by = c("fulcrum_id" = "fulcrum_id_isolation")) %>%
  dplyr::mutate(c_label_id = ifelse(!is.na(fulcrum_id_collection), fulcrum_id_collection, c_label_id)) %>%
  dplyr::select(-c_label, -fulcrum_id_collection)



#=================================================================#
# Part 6: Read in nematode isolation s_plate labels (s_plates)    #
#=================================================================#
s_plates <- readr::read_csv("data/fulcrum/nematode_isolation_s_labeled_plates.csv") %>%
  # Need to filter isolation data from wrong time period
  dplyr::filter(created_at > ymd("2018-11-1")) %>%
  dplyr::select(fulcrum_parent_id, s_label)

# read s_plate data that was accidentally assigned to another project and was C. elegans, C. briggsae, or C. tropicalis positive.
s_plates_cross_project_corrections <- readr::read_csv("data/nematode_isolation_s_labeled_plates_cross_project_corrections.csv") %>%
  dplyr::select(fulcrum_parent_id, s_label)

# Join s_plate data frames together
s_plates <- full_join(s_plates, s_plates_cross_project_corrections)
  
#=================================================#
# Part 7: Read in genotyping data (species_ids)   #
#=================================================#
species_ids <- readr::read_csv("data/20180109_species_id.csv", na = c("", "#N/A"))

#====================================================#
# Part 8: Join genotyping data with s_plate labels   #
#====================================================#
s_labels_and_species_ids <- dplyr::left_join(s_plates, species_ids) %>%
  dplyr::select(fulcrum_parent_id, s_label, confirmed_species_id, ECA_name)

#====================================================================#
# Part 9: Join isolation data with s_plate labels and species ids   #
#====================================================================#
isolation_full <- dplyr::full_join(isolation, s_labels_and_species_ids, by = c("fulcrum_id" = "fulcrum_parent_id"))

#================================================================================================================#
# Part 10: Join collection, isolation, s_plate labels and rename duplicated c_labels that have been corrected    #
#================================================================================================================#
df <- dplyr::left_join(collection, isolation_full, by = c("fulcrum_id" = "c_label_id")) %>%
  select(c_label,
         s_label,
         confirmed_species_id,
         duplicated_c_label,
         duplicated_c_label_corrected,
         fulcrum_id_collection = fulcrum_id,
         fulcrum_id_isolation = fulcrum_id.y,
         collection_datetime,
         isolated_at,
         collection_latitude,
         collection_longitude,
         ambient_temperature_c,
         ambient_humidity,
         substrate_temperature,
         substrate_moisture,
         everything())

# Create unique C-labels for those that were duplicated and corrected
isolation_correction_rename <- isolation_correction %>%
  #filter c_labels that were not duplicated within the HawaiiOctober2018 collection. These are isolations that were mistakenly assigned to collections from other projects.
  dplyr::filter(!(c_label %in% c("C-1063", "C-2571", "C-2727"))) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(correction_vector1 = as.character(seq(n()))) %>%
  dplyr::ungroup() %>%
  tidyr::unite(c_label_cor, c_label, correction_vector1, sep = "-", remove = F) %>%
  dplyr::select(-correction_vector1)

# Replace duplicated C-labels with unique corrected labels
df <- left_join(df, isolation_correction_rename) %>%
  dplyr::mutate(c_label_corrected = ifelse(!is.na(c_label_cor), c_label_cor, c_label)) %>%
  dplyr::select(c_label,
                c_label_corrected,
                s_label,
                confirmed_species_id,
                duplicated_c_label,
                duplicated_c_label_corrected,
                everything()) %>%
  dplyr::select(-c_label_cor)


  
#========================================#
# Part 11: Set the islands and trails    #
#========================================#
# Create Island Column
df$island <- "?"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-158.3617,21.1968,-157.5117,21.7931)), "island"] <- "Oahu"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-159.9362, 21.6523, -159.1782, 22.472)), "island"] <- "Kauai"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.327, 21.0328, -156.685, 21.2574)), "island"] <- "Molokai"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-156.7061, 20.4712, -155.9289, 21.0743)), "island"] <- "Maui"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-156.1346, 18.6619, -154.6985, 20.4492)), "island"] <- "Big Island"
df$location <- NA
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.72537,21.303309,-157.71919,21.32122)), "location"] <- "Kuliouou Ridge Trail"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-158.0192352613,21.5014265529,-158.0145925283,21.5041245046)), "location"] <- "Wahiawa Botanical Garden"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.8598800302,21.3149311581,-157.855797708,21.3182194587)), "location"] <- "Foster Community Garden"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.7829487403,21.3569863645,-157.7752268314,21.3655295525)), "location"] <- "Maunawili Demonstration Trail"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.8014534712,21.3322593,-157.798127532,21.3427719396)), "location"] <- "Manoa Falls Trail"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-157.8135502338,21.3779082884,-157.7915561199,21.3970691079)), "location"] <- "Ho'omaluhia Botanical Garden"
df[filter_box(df$collection_longitude, df$collection_latitude, c(-159.613624,22.167098,-159.575601,22.226422)), "location"] <- "Na Pali Coast State Wilderness Park"

#=============================================================#
# Part 12: attach urls for photos on googledrive              #
#=============================================================#
# Need to authenticate when you run code. As of 1/14/2019 files are on Tim Crombie's google drive
# include automatic authentification here

# Identify photo ids stored on google drive ()
hi_photo_ids <- googledrive::drive_ls("~/HI_photos/") %>%
  dplyr::mutate(sample_photo = str_remove(name, ".jpg")) %>%
  dplyr::select(sample_photo, id)

# Merge googledrive ids with df by sample photo
df <- left_join(df, hi_photo_ids)

#=============================================================#
# Part 13: Exporting dataframe for use with other scripts     #
#=============================================================#
save(file = "data/fulcrum/df.Rda", df)

#=============================================================#
# Part 14: Exporting wild isolate ECA names for labguru       #
#=============================================================#
df_lab_guru <- df %>%
  dplyr::filter(!is.na(confirmed_species_id)) %>%
  dplyr::select(name = ECA_name,
                Alternative_name = ECA_name,
                Organism = confirmed_species_id,
                Isolation = substrate,
                Isolation_date = collection_date,
                GPS_latitude = collection_latitude_photo,
                GPS_longitude = collection_longitude_photo,
                Species = confirmed_species_id,
                Isolation_location = island,
                sampled_in_the_field_by = collected_by
                ) %>%
  dplyr::arrange(Alternative_name)

readr::write_csv(df_lab_guru, "data/lab_guru_upload.csv")

