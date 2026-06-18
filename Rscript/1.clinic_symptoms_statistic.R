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


#######clinic statistic

clinic <- read.xlsx("DATA/TEST.xlsx")
clinic <- clinic[!duplicated(clinic$Patient_ID), ]
colnames(clinic)
myVars <- c("TEM_new","Gene_new","Zygosity","VariantSlassification_new","Onset.age_new",
            "Diagnosis.age_new","Enrollment.age_new","Gender_new","Height(m)_new","Weight(kg)_new",
            "BMI_new","Province_new","FnNO(ppb)_new","nNO.(nL/min)","CBF(nasal.cell)_new",
            "CBP_new2","EFV1%","chronic_respiratory_infections","chronic_wet_cough","bronchiectasis",
            "Pneumonia","bronchitis","pulmonary.atelectasis/pulmonary.consolidation_new","nasosinusitis",
            "Otitis_media/mastoiditis","Hearing_impairment_new","Situs_inversus","Neonatal.pulmonary.symptoms_new",
            "CHD_new","infertility_new","Brain.disease_new","Fetal.lateral.ventricle.widening_new","COVID-19",
            "maternal.miscarriage_new","Kindey.disease_new","Skeletal.abnormalities")

Median <- list()
Q25 <- list()
Q75 <- list()
Mean <- list()
SD <- list()
Non_NA_Count <- list()

a <- c("Onset.age_new","Diagnosis.age_new","Enrollment.age_new","Height(m)_new","Weight(kg)_new",
       "BMI_new","FnNO(ppb)_new","nNO.(nL/min)","CBF(nasal.cell)_new","EFV1%")
for (i in a){
  Median[i] <- median(clinic[,i],na.rm = TRUE)
  Q25[i] <- quantile(clinic[,i], 0.25,na.rm = TRUE)
  Q75[i] <- quantile(clinic[,i], 0.75,na.rm = TRUE)
  Mean[i] <- mean(clinic[,i],na.rm = TRUE)
  SD[i] <- sd(clinic[,i],na.rm = TRUE)
}

summary_table <- data.frame(
  Variable = a,
  Median = unlist(Median),
  Q25 = unlist(Q25),
  Q75 = unlist(Q75),
  Mean = unlist(Mean),
  SD = unlist(SD),
  row.names = NULL
)



for (i in a){
  current_data <- clinic[[i]]
  Non_NA_Count[[i]] <- sum(!is.na(current_data))

  if(Non_NA_Count[[i]] > 0) {
    Median[[i]] <- median(current_data, na.rm = TRUE)
    Q25[[i]] <- quantile(current_data, 0.25, na.rm = TRUE)
    Q75[[i]] <- quantile(current_data, 0.75, na.rm = TRUE)
    Mean[[i]] <- mean(current_data, na.rm = TRUE)
    SD[[i]] <- sd(current_data, na.rm = TRUE)
  } else {
    Median[[i]] <- NA
    Q25[[i]] <- NA
    Q75[[i]] <- NA
    Mean[[i]] <- NA
    SD[[i]] <- NA
  }
}

summary_table <- data.frame(
  Variable = a,
  N = unlist(Non_NA_Count[a]),  
  Median = unlist(Median[a]),
  Q25 = unlist(Q25[a]),
  Q75 = unlist(Q75[a]),
  Mean = unlist(Mean[a]),
  SD = unlist(SD[a]),
  row.names = NULL
)
print(summary_table)


tab1 <-CreateTableOne(vars = myVars, data = clinic, includeNA = TRUE)

tab <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
print(tab)


########## clinic sympotom percent barplot

clinic <- read.xlsx("DATA/TEST.xlsx")
clinic <- clinic[!duplicated(clinic$Patient_ID), ]
colnames(clinic)

clinic2 <- clinic[,c("chronic_respiratory_infections","chronic_wet_cough","bronchiectasis",
                     "Pneumonia","bronchitis","pulmonary.atelectasis/pulmonary.consolidation_new","nasosinusitis",
                     "Otitis_media/mastoiditis","Hearing_impairment_new","Situs_inversus","Neonatal.pulmonary.symptoms_new",
                     "CHD_new","infertility_new","Brain.disease_new","Fetal.lateral.ventricle.widening_new","COVID-19",
                     "maternal.miscarriage_new","Kindey.disease_new","Skeletal.abnormalities")]

symptom_percentages <- clinic2 %>%
  summarise(across(everything(), ~ {
    yes_count <- sum(. == "yes", na.rm = TRUE)
    no_count <- sum(. == "no", na.rm = TRUE)
    if (yes_count + no_count == 0) {
      NA 
    } else {
      yes_count / (yes_count + no_count) * 100
    }
  }))
symptom_percentages <- symptom_percentages %>%
  pivot_longer(everything(), names_to = "Symptom", values_to = "Percentage")

symptom_Counts <- clinic2 %>%
  summarise(across(everything(), ~ {
    yes_count <- sum(. == "yes", na.rm = TRUE)
    no_count <- sum(. == "no", na.rm = TRUE)
    if (yes_count + no_count == 0) {
      NA  
    } else {
      yes_count + no_count
    }
  }))

symptom_Counts <- symptom_Counts %>%
  pivot_longer(everything(), names_to = "Symptom", values_to = "Counts")

symptom <- merge(symptom_percentages,symptom_Counts)

symptom$Symptom


symptom$Percentage <- round(symptom$Percentage,2)
symptom$Symptom_new2 <- paste0(symptom$Symptom," (n=",symptom$Counts,")")


ggplot(symptom, aes(x = reorder(Symptom_new2, Percentage), y = Percentage)) +
  geom_bar(stat = "identity", fill = "#786937", width = 0.8, alpha = 1) +  
  geom_text(aes(label = paste0(Percentage, "%")),  
            hjust = -0.2, size = 4.5, color = "black", fontface = "bold") +  
  coord_flip() +  
  ylim(0,110)+
  labs(title = "Clinical Phenotype",  
       x = "Symptom",  
       y = "Percentage (%)") +  
  scale_y_continuous(
    breaks = seq(0, max(symptom$Percentage) + 10, by = 10),  
    expand = expansion(mult = c(0, 0.1)) 
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),  
    axis.title.x = element_text(face = "bold", size = 14, color = "black"),  
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),  
    axis.text.y = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),  
    axis.line.x = element_line(color = "black", size = 0.5), 
    axis.line.y = element_line(color = "black", size = 0.5),  
    xis.ticks.x = element_line(color = "black", size = 0.5),  
    axis.ticks.y = element_line(color = "black", size = 0.5),  
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.y = element_blank(),  
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)  
  )





