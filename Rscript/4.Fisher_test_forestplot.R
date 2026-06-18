rm(list = ls())
set.seed(10086)
library(tableone)
library(openxlsx)
library(meta)
library(dplyr)
library(tidyr)
library(tidyverse)
library(forestplot)
library(epiR)
library(grid)
sym_clinic <- read.xlsx("DATA/TEST.xlsx")
colnames(sym_clinic)
clinic_symptom <- c(
  "chronic_respiratory_infections","chronic_wet_cough","bronchiectasis",
  "Pneumonia","bronchitis" ,"pulmonary.atelectasis/pulmonary.consolidation_new",
  "nasosinusitis","Otitis_media/mastoiditis",
  "Hearing_impairment_new","Situs_inversus","Neonatal.pulmonary.symptoms_new",
  "CHD_new","infertility_new","Brain.disease_new","Fetal.lateral.ventricle.widening_new",
  "COVID-19","maternal.miscarriage_new","Kindey.disease_new","Skeletal.abnormalities")

contingency_table_sym <- list()


for (a in 1:length(clinic_symptom)){
  sym_gene <- sym_clinic[,c("Gene_new",clinic_symptom[a])]
  sym_gene <- na.omit(sym_gene)
  sym_gene$Clinical_Sign<- NA
  sym_gene$Clinical_Sign[which(sym_gene[,clinic_symptom[a]]=="yes")] <- 1
  sym_gene$Clinical_Sign[which(sym_gene[,clinic_symptom[a]]=="no")] <- 0
  
  contingency_table <- sym_gene %>%
    group_by(Gene_new, Clinical_Sign) %>%
    summarise(count = n()) %>%
    spread(Clinical_Sign, count, fill = 0) %>%
    rename("NO" = `0`, "YES" = `1`)
  
  contingency_table_sym[[a]] <- contingency_table
  names(contingency_table_sym)[a] <- clinic_symptom[a]
}

##############  fishertest function ############

forest_data_f <- function(contingency_table){
  group_get <- as.character(unique(contingency_table$Gene_new)) 
  forest_data_combine <- data.frame(
    Study = character(0),
    mean = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    pval = numeric(0),
    "OR(95% CI)" = character(0),
    Method=character(0),
    YES = numeric(0),
    NO = numeric(0),
    Count = numeric(0)
  )
  
  or_result_list <- list()
  colnames(forest_data_combine) <- c("Study","mean","lower","upper","pval","OR(95% CI)","Method","YES","NO","Count")
  for (i in 1:length(group_get)){
    #数据测试
    # print(i)
    #产生一个临时的2X2列联表，用于计算OR值；
    contingency_table_temp<-contingency_table[which(contingency_table$Gene_new==group_get[i]),]
    
    contingency_table_temp[2, ] <- list("others", 
                                        sum(contingency_table[-which(contingency_table$Gene_new == group_get[i]), "NO"]),
                                        sum(contingency_table[-which(contingency_table$Gene_new == group_get[i]), "YES"]))
    
    contingency_table_temp<-contingency_table_temp[,-1]
    contingency_table_temp <- as.data.frame(lapply(contingency_table_temp, as.numeric))
    rownames(contingency_table_temp)<-c(group_get[i],"others")
    contingency_table_temp<-contingency_table_temp[,c(2,1)]
    
    or_result <- fisher.test(as.matrix(contingency_table_temp))
    forest_data <- data.frame(
      Study = group_get[i],       # 分组标签
      mean = round(or_result$estimate,3) ,             # OR值
      lower = round(or_result$conf.int[1],3)  ,            # 95% CI 下限
      upper =  round(or_result$conf.int[2],3)  ,           # 95% CI 上限
      pval = round(or_result$p.value, 3) 
    )
    forest_data$'OR(95% CI)' <- paste0(forest_data$mean, "(", forest_data$lower, "-", forest_data$upper, ")")
    forest_data$Method <- "Fisher's exact test"
    
    or_result_list[[i]] <- or_result
    names(or_result_list)[i] <- group_get[i]
    
    fre_OR<-cbind(forest_data,contingency_table_temp[1,],sum(contingency_table_temp[1,]))
    forest_data_combine<-rbind(forest_data_combine,fre_OR)
  }
  forest_data_combine <- arrange(forest_data_combine,pval)
  forest_data_combine <- rbind(colnames(forest_data_combine), forest_data_combine)
  return(forest_data_combine)
}


####### forest plot function#######
forest_f <- function(forest_data){
  library(grid)
  library(forestplot)
  forestplot(as.matrix(forest_data[,c("Study","YES","NO","OR(95% CI)","pval")]), 
             mean = forest_data$mean,          
             lower = forest_data$lower,        
             upper = forest_data$upper,        
             zero = 1,                         
             graph.pos=2,   
             graphwidth=unit(70,"mm"),
             lineheight=unit(8,"mm"),
             colgap = unit(8,"mm"),
             boxsize=0.5, 
             is.summary = c(TRUE,rep(FALSE,100)), 
             ci.vertices = TRUE,
             ci.vertices.height=0.1,
             xticks=c(0,0.5,1.0,1,2,3,4), 
             clip=c(0.05,2),
             lwd.ci =2,
             xlab="OR(95% CI)",
             col= fpColors(lines = "#0079AF",
                           box="#0079AF",
                           zero = "black"))
  
}

############### test

forest_data_ls<-lapply(1:length(contingency_table_sym),function(x){
  forest_data_f(contingency_table_sym[[x]])
})

sheetnames <- names(contingency_table_sym)

names(forest_data_ls)<-c(
  "chronic_respiratory_infections","chronic_wet_cough","bronchiectasis",
  "Pneumonia","bronchitis" ,"pulmonary symptom",
  "nasosinusitis","Otitis_media/mastoiditis",
  "Hearing_impairment_new","Situs_inversus","Neonatal pulmonary symptoms",
  "CHD_new","infertility_new","Brain.disease","Lateral.ventricle.widening",
  "COVID-19","maternal.miscarriage","Kindey.disease","Skeletal.abnormalities")
#添加表头

#names(forest_data_ls) <- names(contingency_table_sym)
#names(forest_data_ls)[6] <- "pulmonary disease"
write.xlsx(forest_data_ls,"output260118/forest_data_VariantSlassification_new_clinic2.xlsx")


#####画图

for (i in 1:length(forest_data_ls)) {
  p <- forest_f(forest_data_ls[[i]])
  print(p)
  grid.text(names(forest_data_ls)[i],x = 0.5,
            y = 2)
}


