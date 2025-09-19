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

# load functions we will need
source("scripts/sdmpredictors_fixed.R") 
source("scripts/biooracle_obis_functions.r") 

########
# apply functions to our data -------------------------------------------------------
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
  
  # if Phymatolithon, get Phymatolithon rugulosum
  # and Phymatolithon scabriusculum as it is not 
  # updated in OBIS yes
  if(i %in% c("Phymatolithon rugulosum", "Phymatolithon scabriusculum")) {
    # c(boreo, litho)
    my_sp <- c(157333, 1832160)
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
                      "Occurrence_based_species_thermal_indicies_Photos_20250919.csv"),
          row.names = F)
