#' ----------------------------------------------------
#' Species Thermal Indices (based on code 
#' from Tom Webb and adaptations by Amelia Hesketh and
#' Jake Lawlor)
#' 
#' Code pulls in species sampled by Sebens et al.
#' and calculates their distribution based thermal
#' preferences using layers from OBIS selected for
#' relevance
#' 
#' Current code adapted by Jarrett Byrnes <jarrett.byrnes@umb.edu>
#' for use with Ken Sebens subtidal rock wall data
#' ----------------------------------------------------

# Load packages with pacman to make sure they load properly
pacman::p_load(tidyverse, lubridate, 
               worrms, robis, raster, 
               sdmpredictors, naniar, taxize)

# some packages are now deprecated - but available from 
# github
pacman::p_load_gh("ropensci/bold")

# set wd to right place
setwd(here::here())

source("scripts/sdmpredictors_fixed.R") #

## --------------------------------------------------------------------------------------------------------------------
## Functions for working with Biooracle layers


get_temp_summ_by_sp <- function(sp_id, bo_lc = c("BO_sstmean", "BO21_tempmax_bdmean", 
                                                 "BO21_tempmax_ss", "BO21_tempmean_bdmean",
                                                 "BO21_tempmean_bdmin",
                                                 "BO21_tempmax_bdmin", "BO21_templtmax_bdmin"),
                                save_all_recs = TRUE, use_defaults = TRUE){
  
  # function to get OBIS records for a given species and match to Bio-Oracle data
  # if save_all_recs == TRUE, this full matched dataset will be saved before summarising
  # function returns a summary of temperature affinity of the species
  
  # run the full set of functions #109620
  sp_temp <- get_obis_recs(species_id = sp_id) %>% # the [1] is for cases we have to force multiple IDs
    get_bio_oracle_t(layercodes = bo_lc) %>%
    save_full_recs(save_recs = save_all_recs) %>%
    t_summary(layercodes = bo_lc) %>%
    mutate(species_id = sp_id[1]) %>% # the [1] is for cases we have to force multiple IDs
    dplyr::select(species_id, everything())
  
  # return the temperature summary for the species
  sp_temp
}



## ---- eval = FALSE---------------------------------------------------------------------------------------------------
##get_temp_summ_by_sp(sp_id = my_sp)


## --------------------------------------------------------------------------------------------------------------------
## Functions for getting information from OBIS
get_obis_recs <- function(species_id, missing_check = FALSE,
                          fields = c("decimalLongitude", "decimalLatitude", "depth", "date_year", #"month",
                                     "scientificName", "aphiaID")
){
  # Fuction to get OBIS records for a given species_id, which must be a
  # recognised WoRMS Aphia ID
  
  # NB OBIS returns records from all taxa gathered under the same valid Aphia
  # ID; the aphia ID returned is that of the taxon as recorded, not necessarily
  # the valid ID, so in order that the final dataset is correctly named we add
  # back in the 'correct' ID here as valid_AphiaID
  
  if(missing_check == TRUE){
    # catch invalid / unrecognised AphiaIDs here - but recommend doing this prior to calling these functions
    if(length(checklist(taxonid = species_id)) > 1){
      # get OBIS records for a given species ID, add year and month, set negative and missing depth to 0
      obis_recs <- occurrence(taxonid = species_id) %>% 
        as_tibble() %>% 
        dplyr::select(all_of(fields)) %>% 
        mutate(depth = as.numeric(depth),
               year = formatC(date_year),
              # month = formatC(month, width = 2, flag = "0"),
               depth0 = case_when(
                 is.na(depth) ~ 0,
                 depth < 0 ~ 0,
                 TRUE ~ depth),
               valid_AphiaID = species_id[1]) # the [1] is for when we are looking at multiple species due to taxonomy not catching up yet
    } else {
      # at present just returns an empty tibble, which causes problems with
      # other functions further down the pipeline, hence recommend checking
      # AphiaIDs prior to calling
      obis_recs <- tibble()
    }
  } else {
    obis_recs <- occurrence(taxonid = species_id) %>% 
      as_tibble() %>% 
      dplyr::select(all_of(fields)) %>% 
      mutate(depth = as.numeric(depth),
             year = formatC(date_year),
            # month = formatC(month, width = 2, flag = "0"),
             depth0 = case_when(
               is.na(depth) ~ 0,
               depth < 0 ~ 0,
               TRUE ~ depth),
             valid_AphiaID = species_id[1]) # the [1] is for when we are looking at multiple species due to taxonomy not catching up yet
  }
  # return the OBIS records
  obis_recs
}



## --------------------------------------------------------------------------------------------------------------------
get_bio_oracle_t <- function(obis_recs, layercodes, use_defaults = TRUE){
  # Function to match a set of OBIS occurrence recorods to the specified layers
  # from Bio-ORACLE Set path for where these two temperature datasets will be
  # stored
  bo_path <- ifelse(use_defaults,
                    paste0(file.path(getwd()), "/biooracle"),
                    paste0(bespoke_path, "/biooracle"))
  if(!dir.exists(file.path(bo_path))){
    dir.create(path = file.path(bo_path),
               recursive = FALSE, showWarnings = FALSE)}
  # load the layers
  bo_t_dat <- load_layers(layercodes,
                          equalarea =TRUE,
                          datadir = "biooracle")
  # Turn the OBIS occurrence locations into spatial points
  points <- SpatialPoints(
    obis_recs[,c("decimalLongitude", "decimalLatitude")],
    lonlatproj)
  # Reproject (could avoid this by setting equalarea = FALSE)
  points <- spTransform(points, equalareaproj)
  # Extract values from each layer for each point
  i <-  1:length(names(bo_t_dat))
  bo_temp <- sapply(i, function(i){raster::extract(bo_t_dat[[i]], points)})
  colnames(bo_temp) <- names(bo_t_dat)
  bo_temp <- as_tibble(bo_temp)
  
  # add these temperatures back to the OBIS records and return
  bind_cols(obis_recs, bo_temp)
  
}


## --------------------------------------------------------------------------------------------------------------------
save_full_recs <- function(rec_df, save_recs = TRUE, use_defaults = TRUE, bespoke_path = NULL){
  # if save_recs == TRUE, save the full set of obis records + BO layer values
  # for a species
  
  if(save_recs == TRUE){
    out_path <- ifelse(use_defaults, paste0(file.path(getwd()),
                                            "/t_matched_obis_recs"),
                       paste0(bespoke_path, "/t_matched_obis_recs"))
    
    if(!dir.exists(file.path(out_path))){
      dir.create(path = file.path(out_path), recursive = FALSE, showWarnings = FALSE)}
    
    # paste together the filename
    sp_filename <- paste0("aphia", rec_df$valid_AphiaID[1],
                          "_obis_iap_bo_", Sys.Date(), ".csv")
    
    # write the file
    write_csv(x = rec_df, file = file.path(paste(out_path, sp_filename, sep = "/")))
  }
  
  # Return the (unchanged) data to pass to next function
  rec_df
}



## --------------------------------------------------------------------------------------------------------------------
t_summary <- function(t_matched_dat, layercodes){
  # Function to get a range of  summary stats from a matched obis-bio-oracle data frame
  counts <- summarise(t_matched_dat, n_obis_rec = n())
  missings <- miss_var_summary(dplyr::select(
    t_matched_dat, layercodes))
  missings_df <- t(missings[, "n_miss"])
  colnames(missings_df) <- paste0(pull(missings, variable), "_NA")
  missings_df <- as_tibble(missings_df)
  # define separate functions for 5% and 95% quantiles
  q5 <- function(x, na.rm = TRUE){stats::quantile(x, 0.05, na.rm = TRUE)}
  q95 <- function(x, na.rm = TRUE){stats::quantile(x, 0.95, na.rm = TRUE)}
  
  # get a range of summary stats over all variables in the dataset
  t_stats <- summarise_at(t_matched_dat,
                          vars(layercodes),
                          tibble::lst(mean, min, max, median, sd, mad, q5, q95), na.rm = TRUE)
  
  # Tidy up and return the species-level summary
  t_summ <- bind_cols(counts, missings_df, t_stats)
  
  t_summ
}


########
# apply functios to our data -------------------------------------------------------
########
data <- read.csv("data/Sebens_found_sp_list.csv")

# pull only entries with spaces to ensure Genus species
data <- data %>%
  arrange(Accepted.Name) %>%
  filter(str_detect(Accepted.Name, " "))

# isolate names variable
names <- data %>% pull(Accepted.Name)

# use taxize package to get cleaned verifies name
# gna_verifier pulling from WORMS, inat, EOL, Catalogue of Life
names_clean <- taxize::gna_verifier(names,
                                    data_sources = c(1,9,12,180))%>% 
  filter(submittedName == matchedCanonicalSimple) %>%
  distinct(matchedCanonicalSimple) %>% pull(matchedCanonicalSimple)

# view the entries we cut - should be null
names[!names %in% names_clean]

# Because "Phymatolithon scabriusculum"  is not
# valid, after consultation, we are using the
# older name which has not been updated in
# many later databases and they are synonymous
names_clean <- c(names_clean, "Phymatolithon rugulosum")


# Not in OBIS - but very rare, so no data and
# will cause script to fail - not used in analysis
# so removing
names_clean <- names_clean[!names_clean == "Halisarca nahantensis"]

# Waernia mirabilis has entries from only one site
# so is not reliable for this analysis. Exclude.
names_clean <- names_clean[!names_clean == "Waernia mirabilis"]

# "Alcyonium siderium" is likely incorrectly listed
# as Alcyonium digitatum - A.s. only has one record in Europe
# while digitatum is all over the Atlantic. Likely a taxonomy
# issue, so, using digitatum is appropriate
names_clean <- gsub("Alcyonium siderium", "Alcyonium digitatum", names_clean)



# get occurrence data -----------------------------------------------------

#Pulling STI values for all species

# the wm_name2id() function throws errors when it can't identify species,
# and these errors stop the entire operation, so here we use Jarrett's 
# workaround from an open issue on github: https://github.com/ropensci/worrms/issues/27
get_id <- function(.x){
  
  # if names2 returns an error.....
  rec <- try(wm_records_name(.x))
  # if(class(rec) == "try-error"){return(NA)} #original
   if("try-error" %in% class(rec)){return(NA)}
  
  # return the valid ID
  return(rec$valid_AphiaID[1])
  # 
  # #... then use another option to find alphiaID
  # acc <- subset(rec, rec$status=="accepted")
  # 
  # # if we don't get back anything accepted, use the 
  # # accepted field
  # if(nrow(acc) == 0){
  #   return(rec$valid_AphiaID[1])
  # }
  # 
  # acc$AphiaID[1]
}

# create blank dataframe
coefout <- data.frame()

# go through names list and get spp occurrence info
for(i in names_clean){
 
  # add a clock to print species name as we start
  print(i)

  #my_sp <- wm_name2id(i)
  my_sp <- get_id(i)

  # if boreolithothamnion, also get lithothamnion glaciale
  # as it is not yet updated in OBIS
  if(i == "Boreolithothamnion glaciale") {
    # c(boreo, litho)
    my_sp <- c(1736213, 145170)
  }
  
  db <- get_temp_summ_by_sp(my_sp)
  db$gen_spp <- i
  coefout <- rbind(coefout, db)
  Sys.sleep(5) #so we don't kill the API
  
  print("Completed")
}


#Botryllus schlosseri# save csv ----------------------------------------------------------------
# move name column to beginning
coefout <- coefout %>% relocate(gen_spp)
write.csv(coefout, 
          file = here::here(#"outputs","datasets",
                      "data",
                      "Occurrence_based_species_thermal_indicies_Photos_20250904.csv"),
          row.names = F)

# read in and go from here

coefout <- read.csv("data/Occurrence_based_species_thermal_indicies_Photos_20250605.csv")

# plot data ---------------------------------------------------------------
# maybe put in a separate script?

#color1 <- "#782391"
color2 <- "#782391"
 coefout %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmean_mean) ) %>%
  
  ggplot(aes(x=gen_spp)) +
  
  
  geom_point(aes(y=BO21_tempmax_bdmean_q5),  color = color2, alpha=.5,
             position = position_nudge(x = 0.25))+
  geom_point(aes(y=BO21_tempmax_bdmean_q95), color = color2, alpha=.5,
             position = position_nudge(x = 0.25)) +
  geom_point(aes(y=BO21_tempmax_bdmean_mean), color = color2, alpha=1, size=1.5,
             position = position_nudge(x = 0.25)) +
  geom_segment(aes(xend=gen_spp,
                   y=BO21_tempmax_bdmean_q5,
                   yend=BO21_tempmax_bdmean_q95), color = color2, alpha=.5,
               position = position_nudge(x = 0.25)) +
  annotate(geom="text",
           x=3, y=30, 
           hjust=-0.1, vjust=0.2,
           label = 
             "\n Max temp at mean depth",
           color = color2,
           fontface="bold",
           size=5) +
  
  
  labs(x=NULL, y= expression(paste("Water Temperature in ", "\u00b0C"))) +
   theme_bw() +
  theme(#plot.margin = margin(l=25,b=5,unit="pt"),
        axis.text.x = element_text(angle = -90, hjust = 0))

ggsave("figures/thermal_preference_ranges.jpg")

