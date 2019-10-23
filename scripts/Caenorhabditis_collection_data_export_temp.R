library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('/Users/tim/repos/HawaiiMS/data/fulcrum/df.Rda')

# filter to Caenorhabditis strains
caen_cso <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. briggsae", "C. sp. 53", "C. kamaaina")) %>%
  dplyr::mutate(source_lab = "ECA") %>%
  dplyr::mutate(spp_id = ifelse(spp_id == "C. sp. 53", "C. oiwi", spp_id)) %>%
  dplyr::select(s_label,
                species = spp_id,
                strain,
                source_lab,
                latitude,
                longitude,
                landscape,
                substrate,
                substrate_temp = substrate_temperature,
                ambient_temp = ambient_temperature,
                substrate_moisture,
                ambient_humidity,
                photo = photo_url,
                isolation_date = date)

forzen_coiwi_cso <- caen_cso %>%
  dplyr::filter(species == "C. oiwi", s_label != "S-05077")

# Write the wild isolate info sheet .csv for HW2017 c. oiwi strains
readr::write_csv(forzen_coiwi_cso, 'data/2017_HW_processed_wild_isolate_info_sheet_export.csv')

# join freezing log data with collection data
join_proc_strain_df <- full_join(df, proc_strain_df)

# import column names from wild isolate sheet
df_col_names <- names(data.table::fread('data/HW_2018_strain_info_form.csv'))

# filter to strains that made it through processing. Adding strains ECA1253, ECA1255, ECA1258 b/c they will be frozen before 10/01/2019.
proc_strain_df_final <- join_proc_strain_df %>%
  dplyr::filter(initial_spp_id != "") %>%
  dplyr::mutate(`worms frozen` = ifelse(ECA_name %in% c("ECA1253", "ECA1255", "ECA1258"), "TBA", `worms frozen`)) %>%
  dplyr::filter(`worms frozen` != "") %>%
  dplyr::mutate(reference_strain = NA,
                isotype = NA,
                previous_names = NA,
                sequenced = NA,
                warning_message = NA,
                release = NA,
                source_lab = "ECA",
                collection_group = "hawaii_2018",
                substrate_comment = NA,
                sample_photo_url = ifelse(duplicated_c_label == 1, NA, sample_photo_url), # setting these to blank b/c not sure if duplicated url is correct. Need to check
                isolation_date_comment = NA,
                notes = NA,
                state = "Hawaii",
                country = "United States") %>%
  dplyr::select(species = confirmed_species_id,
                strain = ECA_name,
                reference_strain,
                isotype,
                previous_names,
                sequenced,
                warning_message,
                release,
                source_lab,
                collection_group,
                c_label,
                s_label,
                latitude = collection_latitude_photo,
                longitude = collection_longitude_photo,
                landscape,
                substrate,
                substrate_comment,
                substrate_temp = substrate_temperature,
                ambient_temp = ambient_temperature_c,
                substrate_moisture,
                ambient_humidity,
                photo = sample_photo_url,
                isolated_by,
                sampled_by = collected_by,
                isolation_date,
                isolation_date_comment,
                notes,
                state,
                country)

# Write the wild isolate info sheet .csv
readr::write_csv(proc_strain_df_final, 'data/HW_processed_wild_isolate_info_sheet_export.csv')
