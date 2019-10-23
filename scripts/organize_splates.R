library(tidyverse)
library(googlesheets)
library(lubridate)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data
dat <- read.csv("~/Dropbox/AndersenLab/LabFolders/Tim/2018HawaiiSampling/data/nematode_isolation_ALL_before_edit.csv", header = TRUE)

# filter data to correct s-plates
proc_dat <- dat %>%
  dplyr::filter(project == "HawaiiOctober2018" | project == "",
                date )
