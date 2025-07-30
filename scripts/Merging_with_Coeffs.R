library(tidyverse) # v1.2.1
library(lubridate) # v1.7.4
library(worrms) # v0.2.8
library(robis) # v1.0.1
library(raster) # v2.6-7
library(sdmpredictors) # v0.2.8
library(naniar) # v0.3.1
library(taxize)
library(flextable) 


setwd("~/Dropbox (Byrnes Lab)/Breck_GOM/Data/R_Projects/")

thermal_indicies <- read.csv("sp_thermal_limits/data/Occurrence_based_species_thermal_indicies_Photos_20250605.csv")
View(thermal_indicies)

combined_data <- read_csv("rock_wall_change_thermal_preference/data/coefs_with_indices.csv")
View(combined_data)

thermal_indicies2 <- left_join(combined_data, thermal_indicies, by = "gen_spp")
View(thermal_indicies2)

thermal_indicies2 <- thermal_indicies2 |>
filter(!is.na(gen_spp))

# we should think more about which layers we are using, because there are loads
# of options. read about them here:
layers <- tibble(list_layers(simplify = FALSE) %>% filter(str_detect(name, "Sea water temperature")) %>%
  dplyr::select(layer_code,name) %>% arrange(layer_code))

View(layers)

layers_table <- flextable(layers)
autofit(layers_table)

# plot data ---------------------------------------------------------------

#MAX TEMP AT MEAN DEPTH + SST
color1 <- "#782391"
color2 <- "#f2a23f"

maxtemp_meandepthSST_plot <- thermal_indicies %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmean_mean) ) %>% # BO21... doesn't exist in thermal_indicies 
  
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

#MAX TEMP AT MEAN DEPTH

color2 <- "#f2a23f"

maxtemp_meandepth_Pershing_plot <- thermal_indicies %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmean_mean.x) ) %>% # BO21... doesn't exist in thermal_indicies 
  
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
           hjust=0, vjust=1,
           label = 
             "\n Max temp at mean depth",
           color = color2,
           fontface="bold",
           size=5) +
  geom_hline(yintercept = 13.21, col = "turquoise")+
  
  labs(x=NULL, y= "Water Temperature in C")+
  theme(#plot.margin = margin(l=25,b=5,unit="pt"),
    axis.text.x = element_text(angle = -90, hjust = 0))

#MAX TEMP AT MIN DEPTH (quantiles)

color2 <- "#f2a23f"

maxtemp_mindepth_Pershing_plotq <- thermal_indicies %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmin_mean) ) %>% # BO21... doesn't exist in thermal_indicies 
  
  ggplot(aes(x=gen_spp)) +
  
  geom_point(aes(y=BO21_tempmax_bdmin_q5),  color = color2, alpha=.5,
             position = position_nudge(x = 0.25))+
  geom_point(aes(y=BO21_tempmax_bdmin_q95), color = color2, alpha=.5,
             position = position_nudge(x = 0.25)) +
  geom_point(aes(y=BO21_tempmax_bdmin_mean), color = color2, alpha=1, size=1.5,
             position = position_nudge(x = 0.25)) +
  geom_segment(aes(xend=gen_spp,
                   y=BO21_tempmax_bdmin_q5,
                   yend=BO21_tempmax_bdmin_q95), color = color2, alpha=.5,
               position = position_nudge(x = 0.25)) +
  annotate(geom="text",
           x=3, y=30, 
           hjust=0, vjust=1,
           label = 
             "\n Max temp at min depth",
           color = color2,
           fontface="bold",
           size=5) +
  
  geom_hline(yintercept = 13.21, col = "turquoise")+
  
  labs(x=NULL, y= "Water Temperature in C")+
  theme(#plot.margin = margin(l=25,b=5,unit="pt"),
    axis.text.x = element_text(angle = -90, hjust = 0))

#MAX TEMP AT MIN DEPTH (max and min)
color1 <- "#782391"
color2 <- "#f2a23f"

maxtemp_mindepth_tipping_point_plotminmax <- thermal_indicies2 %>% 
  mutate(gen_spp = forcats::fct_reorder(gen_spp, BO21_tempmax_bdmin_mean.x) ) %>% # BO21... doesn't exist in thermal_indicies 
  
  ggplot(aes(x=gen_spp)) +
  
  geom_point(aes(y=BO21_tempmax_bdmin_min.x),  color = color2, alpha = 1)+
  geom_point(aes(y=BO21_tempmax_bdmin_max.x), color = color2, alpha = 1) +
  geom_segment(aes(xend=gen_spp,
                   y=BO21_tempmax_bdmin_min.x,
                   yend=BO21_tempmax_bdmin_max.x), color = color2, alpha=1) +
  geom_point(aes(y=BO21_tempmax_bdmin_mean.x), color = color1, alpha=1, size=3) +
  theme_classic(base_size = 16)+
  
  annotate(geom="text",
           x=2, y=32, 
           hjust=0, vjust=0.2,
           label = 
             "\n Range of max temp at min depth",
           color = color2,
           fontface="bold",
           size=5) +
  
  geom_hline(yintercept = 15, col = "turquoise")+
  
  labs(x=NULL, y= "Water Temperature in C")+
  theme(#plot.margin = margin(l=25,b=5,unit="pt"),
    axis.text.x = element_text(angle = -90, hjust = 0))

ggsave("sp_thermal_limits/figures/maxtemp_mindepth_tipping_point_plotminmax.jpg")


