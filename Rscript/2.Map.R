############ map
rm(list = ls())
set.seed(10086)
library(tableone)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggsci)
library(ggrepel)
library(scales)
library(tableone)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
install.packages("maps")
install.packages("mapproj")
install.packages('sjmisc')
install.packages("ggspatial")
devtools::install_github("Rimagination/ggmapcn", force = TRUE)
library("maps")
library("mapproj")
library("hchinamap")
library(sf)
library(RColorBrewer)
library(tidyverse)
library(ggspatial)
library(sjmisc)

clinic <- read.xlsx("DATA/TEST.xlsx")
clinic <- clinic[!duplicated(clinic$Patient_ID), ]
colnames(clinic)
map <- clinic[,c("Gene_new","Province_new")]

province_count <- map %>%
  group_by(Province_new) %>%
  summarise(Count = n(), .groups = 'drop')

china_prov =  read_sf("DATA/China.geojson") %>% 
  left_join(province_count, join_by(name_en==Province_new)) %>% 
  mutate(prod_lev = Count) 

china_prov$Count[which(is.na(china_prov$Count))] <- 0

china_line = read_sf("DATA/China_line.geojson")

pal = c("#50b98f", "#92ce6e", "#edda0e", "#e5b107", "#d25c02", "#c01f33")

ggplot() +
  geom_sf(data = china_prov,
          aes(fill = prod_lev),
          color = "#011910", linewidth = 0.15) +
  geom_sf(data = china_line,
          color = "#011910", linewidth = 0.5) +
  geom_sf_text(data = china_prov,  
               aes(label = Count),  
               color = "black",  
               size = 3,  
               check_overlap = TRUE) +  
  coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=104") +
  
  scale_fill_gradient(low = "#FAEBD7", high = "#c01f33", na.value = "grey95", name = "PCD Cases") +
  
  annotation_scale(location = "bl") +  
  annotation_north_arrow(location = "tl", which_north = "false",
                         style = north_arrow_fancy_orienteering) +  
  labs(title = "N=") +  
  
  theme_minimal() +
  theme(legend.position = c(0.25, 0.1),
        legend.direction = 'horizontal',
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        plot.background = element_rect(fill = "white", color = "white"))


