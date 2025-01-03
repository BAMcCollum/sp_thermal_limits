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

thermal_indicies <- read.csv("data/Occurrence_based_species_thermal_indicies.csv")
View(thermal_indicies)

# we should think more about which layers we are using, because there are loads
# of options. read about them here:
list_layers() %>% filter(str_detect(name, "Sea surface|Sea water temperature|Sea bottom")) %>%
  dplyr::select(layer_code,name) %>% arrange(layer_code)
list_layers() %>% filter(str_detect(layer_code,"BO2_tempmean_bdmean")) %>% dplyr::select(name)


# plot data ---------------------------------------------------------------

color1 <- "#782391"
color2 <- "#f2a23f"

thermal_indicies %>% 
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

