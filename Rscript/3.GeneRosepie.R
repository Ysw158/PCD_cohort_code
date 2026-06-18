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
library(eoffice)
library(showtext)

clinic <- read.xlsx("DATA/TEST.xlsx")


gene <- clinic
gene_with_percent <- gene %>%
  group_by(Gene_new) %>% 
  mutate(
    count_in_group = n(),  
    percentage = count_in_group / nrow(gene) * 100,  
    percentage_label = sprintf("%.1f%%", percentage)  
  ) %>%
  ungroup()

gene_with_percent$label <- paste0(gene_with_percent$percentage_label,":",gene_with_percent$Gene_new)


default_colors <- c("#92a4c3","#bebcbb","#baa4c2","#c9ba9b","#cedbe9","#dea296","#adcda3","#e8c49f",
                    "#d484ae","#5b9771","#ebe5c6","#c9b2d1","#9ec4c5","#f6f5b3","#a9793d","#dccf57",
                    "#9cb6d3","#72a348","#a29bc7","#97c5c4","#ac90a0","#d384ae","#28458d","#e8a543",
                    "#9fa0a0","#c0cd84","#ad90a1","grey")


default_colors <- rev(default_colors)

library(ggpie)
ggrosepie(data = gene,
          group_key = c("Gene_new"),
          count_type = "full",
          donut_frac = 0.5,
          donut_label = T,
          #donut_label_size=0,
          donut_label_color="#458B74",
          label_sep=":",
          show_label = T,
          show_tick = FALSE,
          #label_gap = 0.1,
          label_size = 3.8,
          label_color="black",
          tick_break = c(0,10,20,30,40),
          #border_color="#FFF5EE",
          border_size = 0.5,
          label_info="ratio" )+
  scale_fill_manual(values = default_colors)+
  #guides(fill = guide_legend(title = "Variant Type"))+
  labs(fill = "Variant Type") +  
  theme(
    #panel.grid.major.x = element_blank(),  
    panel.grid.major.y = element_blank()  
  )


