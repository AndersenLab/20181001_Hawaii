library(tidyverse)
library(googlesheets)

# note: this script should not be needed once positive controls have been added to lysis worksheet.

# load google sheet. This requires googlesheet to be published to the web. 
dat <- googlesheets::gs_key("1reCq9SaF20CMO8jGYUt-IuQRp7fgRAuGkrNcLj3OI8w") %>%
  googlesheets::gs_read("lysis_sheet")

  
# build positive controls
dat_controls <- tibble(s_label = rep("positive control", times = 117),
                       plate_growth = "control", strip_id = 1:117,
                       well_id = rep(12, times = 117),
                       pcr_plate_id = rep(1:15, each = 8, length.out = 117),
                       lysis_notes = NA,
                       pcr_plate_well = NA,
                       band = NA,
                       sanger_plate_id = NA,
                       sanger_well_id = NA,
                       species_id = NA,
                       notes = NA)

# join positive controls and s-plates
dat_lyse <- rbind(dat_controls, dat) %>%
  dplyr::arrange(strip_id, well_id)

# export .csv with positive controls built in
write.csv(dat_lyse, file = "~/Dropbox/AndersenLab/LabFolders/Tim/2018HawaiiSampling/data/lysis_sheets.csv", na = "", row.names = FALSE)
