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
install.packages("logistf")
library(logistf)  
library(tidyverse)
packageVersion("logistf")


sym_clinic <- read.xlsx("DATA/TEST.xlsx")
#sym_clinic <- sym_clinic[!duplicated(sym_clinic$Patient_ID), ]
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
colnames(sym_clinic)[4] <- "Gene"



################################################## multivariate_analysis_Firth_regression
multivariate_analysis <- function(data, symptom_col, 
                                  covariates = c("Gender_new","Enrollment.age_new"),
                                  debug = FALSE) {
  
  if(!require(logistf, quietly = TRUE)) {
    stop("Please install logistf package: install.packages('logistf')")
  }
  sym_data <- data[, c("Gene", symptom_col, covariates)]
  sym_data <- na.omit(sym_data)
  sym_data$Outcome <- ifelse(sym_data[, symptom_col] == "yes", 1, 0)
  
  if("Gender_new" %in% covariates) {
    sym_data$Gender_new <- as.factor(sym_data$Gender_new)
  }
  
  all_genes <- unique(sym_data$Gene)
  
  results_adjusted <- data.frame()
  results_unadjusted <- data.frame()
  
  firth_success <- list(unadj = 0, adj = 0)
  
  for(gene in all_genes) {
    sym_data$Target_gene <- ifelse(sym_data$Gene == gene, 1, 0)
    
    if(debug) cat("\nProcessing:", gene)
    
    # ----- -----
    tryCatch({
      model_unadj <- logistf(Outcome ~ Target_gene, data = sym_data)
      
      if(debug) cat(" - Unadjusted: SUCCESS")
      firth_success$unadj <- firth_success$unadj + 1
      beta <- coef(model_unadj)["Target_gene"]
      p_val <- model_unadj$prob["Target_gene"] 
      ci_lower <- model_unadj$ci.lower["Target_gene"]
      ci_upper <- model_unadj$ci.upper["Target_gene"]
      results_unadjusted <- rbind(results_unadjusted, data.frame(
        Gene = gene,
        Method_unadj = "Firth",
        OR_unadj = round(exp(beta), 3),
        CI_lower_unadj = round(exp(ci_lower), 3),
        CI_upper_unadj = round(exp(ci_upper), 3),
        P_unadj = p_val
      ))
      
    }, error = function(e) {
      if(debug) cat(" - Unadjusted: FAILED -", e$message)
      
      model_unadj <- glm(Outcome ~ Target_gene, 
                         data = sym_data, 
                         family = binomial(link = "logit"))
      coef_unadj <- summary(model_unadj)$coefficients["Target_gene", ]
      
      results_unadjusted <<- rbind(results_unadjusted, data.frame(
        Gene = gene,
        Method_unadj = "GLM",
        OR_unadj = round(exp(coef_unadj["Estimate"]), 3),
        CI_lower_unadj = round(exp(coef_unadj["Estimate"] - 1.96 * coef_unadj["Std. Error"]), 3),
        CI_upper_unadj = round(exp(coef_unadj["Estimate"] + 1.96 * coef_unadj["Std. Error"]), 3),
        P_unadj = coef_unadj["Pr(>|z|)"]
      ))
    })
    
    # ----- -----
    tryCatch({
      formula_adj <- as.formula(paste("Outcome ~ Target_gene +", 
                                      paste(covariates, collapse = " + ")))
      
      model_adj <- logistf(formula_adj, data = sym_data)
      
      if(debug) cat(" - Adjusted: SUCCESS\n")
      firth_success$adj <- firth_success$adj + 1
      beta <- coef(model_adj)["Target_gene"]
      p_val <- model_adj$prob["Target_gene"]
      ci_lower <- model_adj$ci.lower["Target_gene"]
      ci_upper <- model_adj$ci.upper["Target_gene"]
      results_adjusted <- rbind(results_adjusted, data.frame(
        Gene = gene,
        Method_adj = "Firth",
        OR_adj = round(exp(beta), 3),
        CI_lower_adj = round(exp(ci_lower), 3),
        CI_upper_adj = round(exp(ci_upper), 3),
        P_adj = p_val
      ))
      
    }, error = function(e) {
      if(debug) cat(" - Adjusted: FAILED -", e$message, "\n")
      
      formula_adj <- as.formula(paste("Outcome ~ Target_gene +", 
                                      paste(covariates, collapse = " + ")))
      model_adj <- glm(formula_adj, 
                       data = sym_data, 
                       family = binomial(link = "logit"))
      coef_adj <- summary(model_adj)$coefficients["Target_gene", ]
      
      results_adjusted <<- rbind(results_adjusted, data.frame(
        Gene = gene,
        Method_adj = "GLM",
        OR_adj = round(exp(coef_adj["Estimate"]), 3),
        CI_lower_adj = round(exp(coef_adj["Estimate"] - 1.96 * coef_adj["Std. Error"]), 3),
        CI_upper_adj = round(exp(coef_adj["Estimate"] + 1.96 * coef_adj["Std. Error"]), 3),
        P_adj = coef_adj["Pr(>|z|)"]
      ))
    })
  }
  
  cat("\n=== Firth Regression Summary ===\n")
  cat("Unadjusted models:", firth_success$unadj, "/", length(all_genes), "succeeded\n")
  cat("Adjusted models:", firth_success$adj, "/", length(all_genes), "succeeded\n")
  
  results_unadjusted$FDR_unadj <- p.adjust(results_unadjusted$P_unadj, method = "BH")
  results_adjusted$FDR_adj <- p.adjust(results_adjusted$P_adj, method = "BH")

  final_results <- merge(results_unadjusted, results_adjusted, by = "Gene")

  final_results$OR_change_pct <- abs(final_results$OR_adj - final_results$OR_unadj) / 
    final_results$OR_unadj * 100

  final_results$Conclusion_stable <- (final_results$FDR_unadj < 0.05) == 
    (final_results$FDR_adj < 0.05)

  final_results$OR_CI_unadj <- paste0(final_results$OR_unadj, " (", 
                                      final_results$CI_lower_unadj, "-", 
                                      final_results$CI_upper_unadj, ")")
  final_results$OR_CI_adj <- paste0(final_results$OR_adj, " (", 
                                    final_results$CI_lower_adj, "-", 
                                    final_results$CI_upper_adj, ")")
  

  final_results <- final_results[order(final_results$P_adj), ]
  

  attr(final_results, "covariates") <- covariates
  attr(final_results, "method_summary") <- table(final_results$Method_unadj, 
                                                 final_results$Method_adj)
  attr(final_results, "firth_success") <- firth_success
  
  return(final_results)
}

batch_analysis_with_covariates <- function(data, symptom_list, 
                                           covariates = c("Gender_new","Enrollment.age_new")) {
  
  all_results <- list()
  
  for(i in 1:length(symptom_list)) {
    symptom_name <- symptom_list[i]
    cat("\n========================================\n")
    cat("Processing symptom", i, ":", symptom_name, "\n")
    cat("========================================\n")
    
    cat("\n>>> Multivariate Firth regression (adjusting for", 
        paste(covariates, collapse = ", "), ")...\n")
    adj_results <- multivariate_analysis(data, symptom_name, covariates)

    all_results[[symptom_name]] <- list(
      adjusted = adj_results
    )

    cat("\n>>> Summary for", symptom_name, ":\n")
    cat("  - Total genes analyzed:", nrow(adj_results), "\n")
    cat("  - Significant (unadjusted, FDR<0.05):", 
        sum(adj_results$FDR_unadj < 0.05, na.rm = TRUE), "\n")
    cat("  - Significant (adjusted, FDR<0.05):", 
        sum(adj_results$FDR_adj < 0.05, na.rm = TRUE), "\n")
    cat("  - Genes with stable conclusion:", 
        sum(adj_results$Conclusion_stable, na.rm = TRUE), "/", nrow(adj_results), "\n")
  }
  
  return(all_results)
}


####################################################
results <- batch_analysis_with_covariates(
  data = sym_clinic,
  symptom_list = clinic_symptom,
  covariates = c("Gender_new", "Enrollment.age_new"))




