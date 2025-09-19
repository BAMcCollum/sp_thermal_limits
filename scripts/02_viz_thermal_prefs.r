library(ggplot2)
library(readr)

# read in and go from here

coefout <- read.csv("data/Occurrence_based_species_thermal_indicies_Photos_20250919.csv")

# plot data ---------------------------------------------------------------
# maybe put in a separate script?

#color1 <- "#782391"
color2 <- "#782391"
coefout %>% 
  filter(hasSebensData) |>
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

