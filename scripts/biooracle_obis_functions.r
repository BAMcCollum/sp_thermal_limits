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
    mutate(year = suppressWarnings(as.numeric(year))) |> # some weird stuff with NA years
    add_sites_this_study() |> # puts halfway rock in
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
  
  t_summ |> 
    mutate(scientificName = t_matched_dat$scientificName[1],
           hasSebensData = t_matched_dat$hasSebensData[1])
}

## --------------------------------------------------------------------------------------------------------------------
## Add data from species from Sebens surveys
# for Thermal preference change across all sites
#sp_over_time <- read_csv("data/sp_presence_over_time.csv")

# for shallow-deep analysis
sp_over_time <- read_csv("data/sp_presence_over_time100825.csv")

add_sites_this_study <- function(dat){
  # is this something we are adding information to?
  if(!(dat$scientificName[1] %in% unique(sp_over_time$scientificName))){
    return(dat |>
             mutate(hasSebensData = FALSE) )
  }
  
  # where was this species found
  sp_found_at <- sp_over_time |>
    filter(scientificName == dat$scientificName[1])
  
  sites <- bind_rows(
    
    # note - biooracle grid cells to large to do inner/outer
    # so binding into one site for ease of use
    
    # shag rocks 42.414692, -70.906408
    # inner is 7,8m
    # outer is 8,9,10m
    tibble(
      area = c("Shag Rocks Inner", "Shag Rocks Outer"),
      decimalLongitude = -70.906408,
      decimalLatitude = 42.414692,
    ),
    # halfway rock 42.50227063025665, -70.77500303312874
    # inner is 7,8,15,16m
    # outer is 9,10,12,15m
    tibble(
      area = c("Halfway Rock Inner", "Halfway Rock Outer"),
      decimalLongitude = -70.77500303312874,
      decimalLatitude = 42.50227063025665,
    ),
    # dive beach 42.420571, -70.904461, 7m
    tibble(
      area = "Dive Beach", 
      decimalLongitude = -70.904461,
      decimalLatitude = 42.420571,
    )
  )
  
  # bind occurances to site info
  sp_dat <- left_join(sp_found_at, sites) |>
    dplyr::select(decimalLongitude,
                  decimalLatitude,
                  depth_m) |>
    rename(depth = depth_m)
  
  # put it together in robis format
  newdat <- tibble(
    sp_dat,
    date_year = NA,
    scientificName = dat$scientificName[1],
    aphiaID = dat$aphiaID[1],
    year = 0,
    depth0 = NA,
    valid_AphiaID = dat$valid_AphiaID[1]
  ) |>
    mutate(depth0 = depth)
  
  # return new info
  bind_rows(dat, newdat) |>
    mutate(hasSebensData = TRUE)
}

