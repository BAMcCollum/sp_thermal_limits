## Species Thermal Indices (courtesy of Dr. Tom Webb)##
#ON SEPT 16TH 2022 I changed the temp metric from mean at mean depth to max at mean depth at line 37+
library(tidyverse) # v1.2.1
library(lubridate) # v1.7.4
library(worrms) # v0.2.8
library(robis) # v1.0.1
library(raster) # v2.6-7
library(sdmpredictors) # v0.2.8
library(naniar) # v0.3.1
library(taxize)

setwd(here::here())
#setwd("~/Dropbox (Byrnes Lab)/Breck_GOM/Data/R_Projects/sp_thermal_limits")

# JL notes: -------------------------------------------------------------
# we should think more about which layers we are using, because there are loads
# of options. read about them here:
list_layers() %>% filter(str_detect(name, "Sea surface|Sea water temperature|Sea bottom")) %>%
 dplyr::select(layer_code,name) %>% arrange(layer_code)
list_layers() %>% filter(str_detect(layer_code,"BO2_tempmean_bdmean")) %>% dplyr::select(name)

## ---- eval = FALSE---------------------------------------------------------------------------------------------------
## bespoke_path <- "/my bespoke file path name"


## --------------------------------------------------------------------------------------------------------------------
use_defaults <- TRUE


## --------------------------------------------------------------------------------------------------------------------
#Modiolus modiolus, Metridium senile, Anurida maritima, Fucus distichus, Mastocarpus stellatus, Ceramium virgatum, Hildenbrandia rubra, Strongylocentrotus droebachiensis, Prasiola stipitata, Hiatella arctica, Cryptosula pallasiana, Rhizoclonium tortuosum, Leathesia marina, Phymatolithon lenormandii, Asterias rubens, Elachista fucicola, Lacuna vincta, Anomia simplex, Clathromorphum circumscriptum, Saccharina latissima, Crepidula fornicata, Asterias forbesi, Cancer borealis, Diadumene lineata, Pagurus acadianus, Semibalanus balanoides, Mytilus edulis, Nucella lapillus, Littorina obtusata, Littorina saxatilis, Ascophyllum nodosum, Chondrus crispus, Littorina littorea, Corallina officinalis, Fucus spiralis, Fucus vesiculosus, Tectura testudinalis, Alaria esculenta, Ulva lactuca, Codium fragile, Ulva intestinalis, Chordaria flagelliformis, Ahnfeltia plicata, Plumaria plumosa
my_sp <- wm_name2id()
my_sp


## --------------------------------------------------------------------------------------------------------------------
get_temp_summ_by_sp <- function(sp_id, bo_lc = c("BO_sstmean", "BO21_tempmax_bdmean", ""),
                                save_all_recs = TRUE, use_defaults = use_defaults){
  
  # function to get OBIS records for a given species and match to Bio-Oracle data
  # if save_all_recs == TRUE, this full matched dataset will be saved before summarising
  # function returns a summary of temperature affinity of the species
  
  # run the full set of functions
  sp_temp <- get_obis_recs(species_id = sp_id) %>%
    get_bio_oracle_t(layercodes = bo_lc) %>%
    save_full_recs(save_recs = save_all_recs) %>%
    t_summary(layercodes = bo_lc) %>%
    mutate(species_id = sp_id) %>%
    dplyr::select(species_id, everything())
  
  # return the temperature summary for the species
  sp_temp
}



## ---- eval = FALSE---------------------------------------------------------------------------------------------------
## get_temp_summ_by_sp(sp_id = my_sp)


## --------------------------------------------------------------------------------------------------------------------
get_obis_recs <- function(species_id, missing_check = FALSE,
                          fields = c("decimalLongitude", "decimalLatitude", "depth", "year", "month",
                                     "scientificName", "aphiaID")
){
  # Fuction to get OBIS records for a given species_id, which must be a recognised WoRMS Aphia ID
  
  # NB OBIS returns records from all taxa gathered under the same valid Aphia ID; the aphia ID returned is that of the taxon as recorded, not necessarily the valid ID, so in order that the final dataset is correctly named we add back in the 'correct' ID here as valid_AphiaID
  
  if(missing_check == TRUE){
    # catch invalid / unrecognised AphiaIDs here - but recommend doing this prior to calling these functions
    if(length(checklist(taxonid = species_id)) > 1){
      # get OBIS records for a given species ID, add year and month, set negative and missing depth to 0
      obis_recs <- occurrence(taxonid = species_id) %>% 
        as_tibble() %>% 
        dplyr::select(fields) %>% 
        mutate(depth = as.numeric(depth),
               year = formatC(year),
               month = formatC(month, width = 2, flag = "0"),
               depth0 = case_when(
                 is.na(depth) ~ 0,
                 depth < 0 ~ 0,
                 TRUE ~ depth),
               valid_AphiaID = species_id)
    } else {
      # at present just returns an empty tibble, which causes problems with other functions further down the pipeline, hence recommend checking AphiaIDs prior to calling
      obis_recs <- tibble()
    }
  } else {
    obis_recs <- occurrence(taxonid = species_id) %>% 
      as_tibble() %>% 
      dplyr::select(fields) %>% 
      mutate(depth = as.numeric(depth),
             year = formatC(year),
             month = formatC(month, width = 2, flag = "0"),
             depth0 = case_when(
               is.na(depth) ~ 0,
               depth < 0 ~ 0,
               TRUE ~ depth),
             valid_AphiaID = species_id)
  }
  # return the OBIS records
  obis_recs
}



## --------------------------------------------------------------------------------------------------------------------
get_bio_oracle_t <- function(obis_recs, layercodes, use_defaults = TRUE){
  # Function to match a set of OBIS occurrence recrods to the specified layers from Bio-ORACLE
  # Set path for where these two temperature datasets will be stored
  bo_path <- ifelse(use_defaults,
                    paste0(file.path(getwd()), "/biooracle"),
                    paste0(bespoke_path, "/biooracle"))
  if(!dir.exists(file.path(bo_path))){
    dir.create(path = file.path(bo_path),
               recursive = FALSE, showWarnings = FALSE)}
  # load the layers
  bo_t_dat <- load_layers(layercodes,
                          equalarea = TRUE, datadir = "biooracle")
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
  # if save_recs == TRUE, save the full set of obis records + BO layer values for a species
  
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
    write_csv(x = rec_df, path = file.path(paste(out_path, sp_filename, sep = "/")))
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






# apply to our data -------------------------------------------------------
# Alexis's Edits with Jake's Additiona
#webbdata= read_csv(url("https://raw.githubusercontent.com/tomjwebb/occurrence-derived-thermal-affinity/master/data/t_matched_globtherm_dat_full.csv"))
#data <- read.csv("./McCollum_Sebens_sp_list.csv")
View(data)

data <- read.csv("./Sebens_found_sp_list.csv")

# pull only entries with spaces to ensure Genus species
data <- data %>%
  arrange(name) %>%
  filter(str_detect(name, " "))

# isolate names variable
names <- data %>% pull(name)

# use taxize package
names_clean <- taxize::gnr_resolve(names) %>% 
  filter(user_supplied_name == matched_name) %>%
  distinct(matched_name) %>% pull(matched_name)

# view the entries we cut
names[!names %in% names_clean]

# make some manual changes bc these names don't work as is
# (usually it's because the name is not accepted anymore in worms)
#names_clean <- recode(names_clean, 
       "Carcinus maenus" =  "Carcinus maenas",
       "Giffordia granulosa" = "Hincksia granulosa",
       "Hippothoa hyalina" = "Celleporella hyalina",
       "Idotea baltica" = "Idotea balthica",
       "Nemalion helminthoides" = "Nemalion elminthoides",
       "Tectura testudinalis" = "Testudinalia testudinalis",
       "Tonicella rubra" = "Boreochiton ruber",
       "Botryllus schlosseri" = "Botryllus schlosseri",
       "Bugula turrita" = "Crisularia turrita",
       "Cuthona gymnota" = "Catriona gymnota",
       "Gellius arcoferus" = "Hemigellius arcofer",
       "Halichondria panicea" = "Halichondria (Halichondria) panicea",
       "Haliclona oculata" = "Haliclona (Haliclona) oculata",
       "Microciona prolifera" = "Clathria (Clathria) prolifera",
       "Myxilla fimbriata" = "Myxilla (Myxilla) fimbriata",
       "Notoacmea testudinalis" = "Testudinalia testudinalis",
       "Pholis gunnelus" = "Pholis gunnellus",
       "Porania pulvillus" = "Porania (Porania) pulvillus",
       "Sagartia elegans" = "Cylista elegans",
       "Scrupocellaria scabra" = "Aquiloniella scabra",
       "Tritonia plebeia" = "Duvaucelia plebeia",
        "Halisarca nahtantensis" = "Halisarca nahantensis",
       "Bugula turrita" = "Crisularia turrita",
       "Clathria prolifera" = "Clathria (Clathria) prolifera",
       "Flabellina verrucosa" = "Coryphella verrucosa",
       "Hymedesmia paupertas" = "Hymedesmia (Hymedesmia) paupertas",
       "Leptasterias polaris" = "Leptasterias (Hexasterias) polaris") %>% unique()

names_clean <- recode(names_clean,
                      "Tubularia crocea" = "Ectopleura crocea",
                      "Spirorbis borealis" = "Spirorbis (Spirorbis) spirorbis",
                      "Haliclona oculata" = "Haliclona (Haliclona) oculata",
                      "Halichondria panicea" = "Halichondria (Halichondria) panicea",
                      "Hymedesmia paupertas" = "Hymedesmia (Hymedesmia) paupertas",
                      "Clathria prolifera" = "Clathria (Clathria) prolifera") %>% unique()

# manually remove some that didn't work - not sure why,
# errors said these entries didn't have depth values
names_clean <- names_clean[!names_clean == "Halisarca nahantensis"]
names_clean <- names_clean[!names_clean == "Phymatolithon rugulosum"]

#names_clean <- names_clean[!names_clean =="Micrura affinis"]
#names_clean <- names_clean[!names_clean == "Protectocarpus speciosus"]
#names_clean <- names_clean[!names_clean == "Balanus spp."]
#names_clean <- names_clean[!names_clean == "Libinia spp."]
#names_clean <- names_clean[!names_clean == "Obelia spp."]
#names_clean <- names_clean[!names_clean == "Pagurus spp."]
#names_clean <- names_clean[!names_clean == "Aplysilla longispina"]
#names_clean <- names_clean[!names_clean == "Fagesia lineata"]




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
  
  #... then use another option to find alphiaID
  acc <- subset(rec, rec$status=="accepted")
  acc$AphiaID[1]
}

# create blank dataframe
coefout <- data.frame()

# go through names list and get spp occurrence info
for(i in names_clean){
 
  # add a clock to print species name as we start
  print(i)
  
  #my_sp <- wm_name2id(i)
  my_sp <- get_id(i)
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
                      "Occurrence_based_species_thermal_indicies_Photos.csv"),
          row.names = F)
read.csv(coefout, 
         file = here::here(#"outputs","datasets",
           "data",
           "Occurrence_based_species_thermal_indicies_Photos.csv"),
         row.names = F)
View(coefout)

# plot data ---------------------------------------------------------------

color1 <- "#782391"
color2 <- "#f2a23f"
 coefout %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmean_mean) ) %>%
  
  ggplot(aes(x=gen_spp)) +
  
  geom_point(aes(y=BO_sstmean_q5), color = color1, alpha=.5)+
  geom_point(aes(y=BO_sstmean_q95),  color = color1, alpha=.5) +
  geom_point(aes(y=BO_sstmean_mean),  color = color1, alpha=1, size=1.5 ) +
  geom_segment(aes(xend=gen_spp,
                   y=BO_sstmean_q5,
                   yend=BO_sstmean_q95),  color = color1, alpha=.5) +
  annotate(geom="text",
           x=3, y=30, 
           hjust=0, vjust=1,
           label = "Mean SST",
           color = color1,
           fontface="bold",
           size=5) +
  
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
           hjust=0, vjust=1,
           label = 
             "\n Max temp at mean depth",
           color = color2,
           fontface="bold",
           size=5) +
  
  
  labs(x=NULL, y= "Water Temperature in C")+
  theme(#plot.margin = margin(l=25,b=5,unit="pt"),
        axis.text.x = element_text(angle = -90, hjust = 0))
r

