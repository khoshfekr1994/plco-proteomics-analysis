# ========================================================
# Project Title: Comprehensive PLCO Proteomics Analysis
# Author: Hamid Khoshfekr Rudsari, Ph.D.
# Contact: hkhoshfekr@mdanderson.org khoshfekr1994@gmail.com
# Date: November, 2024
# ========================================================
# Description: Unified pipeline for protein and peptide-level proteomics
# analysis of PLCO data with support for IgB and FT partitions
# ========================================================

# Load required libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(parallel)
  library(tidyr)
  library(writexl)
})

# Configuration ----
CONFIG <- list(
  # Analysis parameters
  cut_off_95 = 0.95,
  cut_off_86 = 0.86,
  protein_column = 15,
  protein_column_sum = 10,
  
  # Cancer types
  cancer_types = c("Bladder", "Lung", "Ovarian", "Breast", "Colorectal", "Pancreatic"),
  scenarios = 1:2,
  
  # File paths (modify these according to your setup)
  data_dir = "/Users/hkhoshfekr/Library/CloudStorage/OneDrive-InsideMDAnderson/Projects/proteomics/PLCO/main_data/",
  output_dir = "/Users/hkhoshfekr/OneDrive - Inside MD Anderson/Projects/proteomics/PLCO/cmd/results/"
)

# Display startup message
cat("========================================================\n",
    "    PLCO Comprehensive Proteomics Analysis Pipeline\n",
    "    Author: Hamid Khoshfekr Rudsari, Ph.D.\n",
    "    Contact: hkhoshfekr@mdanderson.org\n",
    "    Date: November, 2024\n",
    "========================================================\n",
    "    Analysis Types: Protein & Peptide Level\n",
    "    Data Types: IgB & FT Partitions\n",
    "    Missing Value Strategies: Zero & NA Assignment\n",
    "========================================================\n")

# Utility Functions ----

#' Create output directory with date stamp
#' @param base_dir Base directory path
#' @param analysis_type Type of analysis (e.g., "protein", "peptide")
#' @param missing_strategy Missing value strategy ("zero_assigned" or "NA_assigned")
create_output_dir <- function(base_dir, analysis_type, missing_strategy) {
  dir_path <- file.path(base_dir, analysis_type, missing_strategy, format(Sys.Date(), "%B_%d"))
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}

#' Load and preprocess data
#' @param file_path Path to CSV file
#' @param missing_strategy Strategy for handling missing values
#' @param protein_column Starting column for protein data
load_data <- function(file_path, missing_strategy = "zero_assigned", protein_column = 15) {
  cat("Loading data from:", basename(file_path), "\n")
  
  if (!file.exists(file_path)) {
    stop("Data file not found: ", file_path)
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Remove last two columns if they exist (cleanup)
  if (ncol(data) > 500) {
    last_cols <- (ncol(data)-1):ncol(data)
    if (all(is.na(data[, last_cols]) | data[, last_cols] == "")) {
      data <- data[, -last_cols]
    }
  }
  
  # Handle missing values
  if (missing_strategy == "zero_assigned") {
    data[, protein_column:ncol(data)] <- lapply(data[, protein_column:ncol(data)], 
                                                function(x) ifelse(is.na(x), 0, x))
  } else if (missing_strategy == "NA_to_zero_for_peptide") {
    # Special case for peptide analysis
    data$Quantity[data$Quantity == 0] <- NA
  }
  
  cat("Data loaded successfully. Dimensions:", nrow(data), "x", ncol(data), "\n")
  return(data)
}

# Core Analysis Functions ----

#' Sum protein values from two dataframes (for FT + IgB combination)
sum_protein_dfs <- function(df1, df2, protein_column = 15, use_parallel = TRUE) {
  cat("Combining protein data from two datasets...\n")
  
  # Get protein column names
  protein_cols1 <- names(df1[, -c(1:(protein_column - 1))])
  protein_cols2 <- names(df2[, -c(1:(protein_column - 1))])
  all_proteins <- unique(c(protein_cols1, protein_cols2))
  
  # Determine reference dataframe (larger one)
  df_large <- max(nrow(df1), nrow(df2))
  df_ref <- if (nrow(df1) == df_large) df1 else df2
  
  # Initialize result dataframe
  df3 <- data.frame(
    Subject_ID = df_ref$Subject_ID,
    age = df_ref$age,
    cig_stat = df_ref$cig_stat,
    j_breast_er_status = df_ref$j_breast_er_status,
    j_breast_her2summ = df_ref$j_breast_her2summ,
    j_breast_pr_status = df_ref$j_breast_pr_status,
    cancer_type = df_ref$cancer_type,
    Gender = df_ref$Gender,
    Stage = df_ref$Stage
  )
  
  # Process proteins
  cat("Processing", length(all_proteins), "proteins...\n")
  
  if (use_parallel && length(all_proteins) > 100) {
    # Parallel processing for large datasets
    cl <- makeCluster(parallel::detectCores() - 1)
    clusterExport(cl, c("df1", "df2", "all_proteins", "df3"))
    
    results <- parLapply(cl, all_proteins, function(protein) {
      values <- numeric(nrow(df3))
      for (i in 1:nrow(df3)) {
        subject_id <- df3$Subject_ID[i]
        val1 <- if (protein %in% names(df1)) {
          idx <- which(df1$Subject_ID == subject_id)
          if (length(idx) > 0) df1[idx[1], protein] else 0
        } else 0
        
        val2 <- if (protein %in% names(df2)) {
          idx <- which(df2$Subject_ID == subject_id)
          if (length(idx) > 0) df2[idx[1], protein] else 0
        } else 0
        
        values[i] <- sum(c(val1, val2), na.rm = TRUE)
        if (values[i] == 0) values[i] <- NA
      }
      return(values)
    })
    
    stopCluster(cl)
    
    # Add results to dataframe
    for (i in seq_along(all_proteins)) {
      df3[[all_proteins[i]]] <- results[[i]]
    }
  } else {
    # Sequential processing
    for (protein in all_proteins) {
      values <- numeric(nrow(df3))
      for (i in 1:nrow(df3)) {
        subject_id <- df3$Subject_ID[i]
        val1 <- if (protein %in% names(df1)) {
          idx <- which(df1$Subject_ID == subject_id)
          if (length(idx) > 0) df1[idx[1], protein] else 0
        } else 0
        
        val2 <- if (protein %in% names(df2)) {
          idx <- which(df2$Subject_ID == subject_id)
          if (length(idx) > 0) df2[idx[1], protein] else 0
        } else 0
        
        values[i] <- sum(c(val1, val2), na.rm = TRUE)
        if (values[i] == 0) values[i] <- NA
      }
      df3[[protein]] <- values
    }
  }
  
  cat("Dataset combination completed.\n")
  return(df3)
}

#' Calculate fold change
calc_fold_change <- function(case_values, control_values) {
  mean_case <- mean(case_values, na.rm = TRUE)
  mean_control <- mean(control_values, na.rm = TRUE)
  FC <- ifelse(mean_control == 0, NA, mean_case / mean_control)
  return(FC)
}

#' Calculate comprehensive ROC metrics
calculate_roc_metrics <- function(data, value_column = "Quantity") {
  # Handle both protein and peptide data
  if (value_column == "Quantity") {
    values <- data$Quantity
  } else {
    values <- data[[value_column]]
  }
  
  # Skip if insufficient data
  if (sum(!is.na(values)) < 10) {
    return(list(AUC = NA, AUC_CI_L = NA, AUC_CI_H = NA, p_value = NA))
  }
  
  roc_result <- suppressMessages(
    tryCatch({
      roc(data$Label, values, direction = "<", quiet = TRUE)
    }, error = function(e) return(NULL))
  )
  
  if (is.null(roc_result)) {
    return(list(AUC = NA, AUC_CI_L = NA, AUC_CI_H = NA, p_value = NA))
  }
  
  auc_value <- auc(roc_result)
  dfROC <- data.frame(
    Spec = roc_result$specificities,
    Sens = roc_result$sensitivities,
    Thresh = roc_result$thresholds
  )
  
  # Calculate metrics at different cutoffs
  calc_cutoff_metrics <- function(cutoff_pct) {
    spec_at_sens <- dfROC[which.min(abs(dfROC$Sens - cutoff_pct)), 1]
    sens_at_spec <- dfROC[which.min(abs(dfROC$Spec - cutoff_pct)), 2]
    thresh_at_sens <- dfROC[which.min(abs(dfROC$Sens - cutoff_pct)), 3]
    thresh_at_spec <- dfROC[which.min(abs(dfROC$Spec - cutoff_pct)), 3]
    
    return(list(
      spec_at_sens = spec_at_sens,
      sens_at_spec = sens_at_spec,
      thresh_at_sens = thresh_at_sens,
      thresh_at_spec = thresh_at_spec
    ))
  }
  
  metrics_95 <- calc_cutoff_metrics(CONFIG$cut_off_95)
  metrics_86 <- calc_cutoff_metrics(CONFIG$cut_off_86)
  
  # Calculate additional performance metrics
  calc_performance <- function(threshold) {
    if (is.na(threshold)) return(list(TP = NA, TN = NA, FP = NA, FN = NA, Accuracy = NA, Precision = NA))
    
    predictions <- ifelse(values > threshold, 1, 0)
    actuals <- ifelse(data$Label, 1, 0)
    
    TP <- sum(predictions == 1 & actuals == 1, na.rm = TRUE)
    TN <- sum(predictions == 0 & actuals == 0, na.rm = TRUE)
    FP <- sum(predictions == 1 & actuals == 0, na.rm = TRUE)
    FN <- sum(predictions == 0 & actuals == 1, na.rm = TRUE)
    
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    precision <- if (TP + FP > 0) TP / (TP + FP) else NA
    
    return(list(TP = TP, TN = TN, FP = FP, FN = FN, Accuracy = accuracy, Precision = precision))
  }
  
  perf_95_spec <- calc_performance(metrics_95$thresh_at_spec)
  perf_95_sens <- calc_performance(metrics_95$thresh_at_sens)
  perf_86_spec <- calc_performance(metrics_86$thresh_at_spec)
  perf_86_sens <- calc_performance(metrics_86$thresh_at_sens)
  
  # Calculate confidence interval and p-value
  ci <- tryCatch(ci.auc(roc_result), error = function(e) c(NA, NA, NA))
  p_value <- tryCatch({
    2 * pt(-abs((auc_value - 0.5) / sqrt(pROC::var(roc_result))), df = Inf)
  }, error = function(e) NA)
  
  return(list(
    AUC = auc_value,
    AUC_CI_L = ci[1],
    AUC_CI_H = ci[3],
    p_value = p_value,
    
    # 95% cutoff metrics
    Sensitivity_At_95_Specificity = metrics_95$sens_at_spec,
    Specificity_At_95_Sensitivity = metrics_95$spec_at_sens,
    Cutoff_At_95_Specificity = metrics_95$thresh_at_spec,
    Cutoff_At_95_Sensitivity = metrics_95$thresh_at_sens,
    
    # 86% cutoff metrics
    Sensitivity_At_86_Specificity = metrics_86$sens_at_spec,
    Specificity_At_86_Sensitivity = metrics_86$spec_at_sens,
    Cutoff_At_86_Specificity = metrics_86$thresh_at_spec,
    Cutoff_At_86_Sensitivity = metrics_86$thresh_at_sens,
    
    # Performance at 95% specificity
    TP_At_95_Specificity = perf_95_spec$TP,
    TN_At_95_Specificity = perf_95_spec$TN,
    FP_At_95_Specificity = perf_95_spec$FP,
    FN_At_95_Specificity = perf_95_spec$FN,
    Accuracy_At_95_Specificity = perf_95_spec$Accuracy,
    Precision_At_95_Specificity = perf_95_spec$Precision,
    
    # Performance at 95% sensitivity  
    TP_At_95_Sensitivity = perf_95_sens$TP,
    TN_At_95_Sensitivity = perf_95_sens$TN,
    FP_At_95_Sensitivity = perf_95_sens$FP,
    FN_At_95_Sensitivity = perf_95_sens$FN,
    Accuracy_At_95_Sensitivity = perf_95_sens$Accuracy,
    Precision_At_95_Sensitivity = perf_95_sens$Precision,
    
    # Performance at 86% specificity
    TP_At_86_Specificity = perf_86_spec$TP,
    TN_At_86_Specificity = perf_86_spec$TN,
    FP_At_86_Specificity = perf_86_spec$FP,
    FN_At_86_Specificity = perf_86_spec$FN,
    Accuracy_At_86_Specificity = perf_86_spec$Accuracy,
    Precision_At_86_Specificity = perf_86_spec$Precision,
    
    # Performance at 86% sensitivity
    TP_At_86_Sensitivity = perf_86_sens$TP,
    TN_At_86_Sensitivity = perf_86_sens$TN,
    FP_At_86_Sensitivity = perf_86_sens$FP,
    FN_At_86_Sensitivity = perf_86_sens$FN,
    Accuracy_At_86_Sensitivity = perf_86_sens$Accuracy,
    Precision_At_86_Sensitivity = perf_86_sens$Precision
  ))
}

#' Determine breast cancer subtype
determine_subtype <- function(er, pr, her2) {
  if (is.na(er) || er %in% c(".N", "7") || 
      is.na(pr) || pr %in% c(".N", "7", "5") ||
      is.na(her2) || her2 %in% c(".N", "5")) {
    return("Unknown")
  }
  
  her2_num <- as.numeric(as.character(her2))
  
  if (er == "3" && pr == "3" && her2_num == 1) {
    return("Luminal")
  } else if (her2_num == 3) {
    return("HER2")
  } else if (er == "1" && pr == "1" && her2_num == 1) {
    return("TNBC")
  } else {
    return("Other")
  }
}

# Main Analysis Functions ----

#' Process protein-level analysis
process_protein_combination <- function(dataset, params_row, missing_strategy = "zero_assigned") {
  cancer <- as.character(params_row$cancer)
  protein <- as.character(params_row$protein)
  scenario <- params_row$scenario
  stage <- params_row$stage
  
  # Gender filtering
  if (cancer %in% c("Breast", "Ovarian")) {
    data_filtered <- dataset %>% filter(Gender == "F")
  } else {
    data_filtered <- dataset
  }
  
  # Select relevant columns
  data_filtered <- data_filtered[, c("Subject_ID", "age", "cig_stat", "j_breast_er_status",
                                     "j_breast_her2summ", "j_breast_pr_status", "cancer_type", 
                                     "Gender", "Stage", protein)]
  
  # Define cases and controls based on scenario
  data_scenario <- data_filtered %>%
    mutate(Label = case_when(
      scenario == 1 ~ cancer_type == cancer,
      scenario == 2 ~ cancer_type == cancer,
      scenario == 3 ~ cancer_type == cancer & Stage == 1,
      scenario == 4 ~ cancer_type == cancer & Stage == 1,
      scenario == 5 ~ cancer_type == cancer & Stage == 2,
      scenario == 6 ~ cancer_type == cancer & Stage == 2,
      scenario == 7 ~ cancer_type == cancer & Stage == 3,
      scenario == 8 ~ cancer_type == cancer & Stage == 3,
      scenario == 9 ~ cancer_type == cancer & Stage == 4,
      scenario == 10 ~ cancer_type == cancer & Stage == 4,
      TRUE ~ FALSE
    ))
  
  # Filter based on scenario (odd = Non-Case controls, even = other cancer controls)
  if (scenario %% 2 == 1) {
    data_scenario <- data_scenario %>% filter(Label | cancer_type == "Non-Case")
  } else {
    data_scenario <- data_scenario %>% filter(Label | !(cancer_type %in% c("Non-Case", cancer)))
  }
  
  # Calculate basic statistics
  age_cases <- data_scenario %>% filter(Label) %>% pull(age)
  age_controls <- data_scenario %>% filter(!Label) %>% pull(age)
  
  case_values <- data_scenario %>% filter(Label) %>% pull(!!sym(protein))
  control_values <- data_scenario %>% filter(!Label) %>% pull(!!sym(protein))
  
  # Handle missing values based on strategy
  if (missing_strategy == "zero_assigned") {
    total_cases_valid <- sum(case_values != 0, na.rm = TRUE)
    total_controls_valid <- sum(control_values != 0, na.rm = TRUE)
    missing_cases <- sum(case_values == 0, na.rm = TRUE)
    missing_controls <- sum(control_values == 0, na.rm = TRUE)
  } else {
    total_cases_valid <- sum(!is.na(case_values))
    total_controls_valid <- sum(!is.na(control_values))
    missing_cases <- sum(is.na(case_values))
    missing_controls <- sum(is.na(control_values))
  }
  
  total_cases <- length(case_values)
  total_controls <- length(control_values)
  
  # Calculate demographics and clinical characteristics
  demo_stats <- calculate_demographics(data_scenario, cancer, protein, missing_strategy, 
                                       total_cases_valid, total_controls_valid)
  
  # Perform statistical analysis if sufficient data
  if (total_cases_valid > 1 && total_controls_valid > 1) {
    roc_metrics <- calculate_roc_metrics(data_scenario, protein)
    fold_change <- calc_fold_change(case_values, control_values)
    
    p_value_fold <- if (total_cases_valid >= 2 && total_controls_valid >= 2) {
      tryCatch(t.test(case_values, control_values)$p.value, error = function(e) NA)
    } else NA
    
    return(create_result_list(cancer, protein, scenario, stage, demo_stats, roc_metrics,
                              fold_change, p_value_fold, case_values, control_values,
                              total_cases_valid, total_controls_valid, 
                              missing_cases, missing_controls, total_cases, total_controls))
  } else {
    return(create_empty_result_list(cancer, protein, scenario, stage, demo_stats,
                                    case_values, control_values, total_cases_valid, 
                                    total_controls_valid, missing_cases, missing_controls, 
                                    total_cases, total_controls))
  }
}

#' Process peptide-level analysis
process_peptide_combination <- function(dataset, params_row, missing_strategy = "zero_assigned") {
  cancer <- as.character(params_row$cancer)
  protein <- as.character(params_row$protein)
  scenario <- params_row$scenario
  stage <- params_row$stage
  
  # Gender filtering
  if (cancer %in% c("Breast", "Ovarian")) {
    df <- dataset %>% filter(Gender == "F")
  } else {
    df <- dataset
  }
  
  df <- df[df$PG.ProteinNames == protein, ]
  peptides <- unique(df$EG.PrecursorId)
  results <- list()
  
  for (peptide in peptides) {
    data_filtered <- df[df$EG.PrecursorId == peptide, ]
    
    # Get most frequent PTM flag
    table_values <- table(data_filtered$PTM_Flag)
    PTM <- names(which.max(table_values))
    
    # Define cases and controls
    data_scenario <- data_filtered %>%
      rowwise() %>%
      mutate(Label = case_when(
        scenario == 1 ~ cancer_type == cancer,
        scenario == 2 ~ cancer_type == cancer,
        scenario == 3 ~ cancer_type == cancer & Stage == 1,
        scenario == 4 ~ cancer_type == cancer & Stage == 1,
        scenario == 5 ~ cancer_type == cancer & Stage == 2,
        scenario == 6 ~ cancer_type == cancer & Stage == 2,
        scenario == 7 ~ cancer_type == cancer & Stage == 3,
        scenario == 8 ~ cancer_type == cancer & Stage == 3,
        scenario == 9 ~ cancer_type == cancer & Stage == 4,
        scenario == 10 ~ cancer_type == cancer & Stage == 4,
        TRUE ~ FALSE
      ))
    
    # Filter based on scenario
    if (scenario %% 2 == 1) {
      data_scenario <- data_scenario %>% filter(Label | cancer_type == "Non-Case")
    } else {
      data_scenario <- data_scenario %>% filter(Label | !(cancer_type %in% c("Non-Case", cancer)))
    }
    
    # Calculate statistics similar to protein analysis
    Protein_group <- unique(data_scenario$PG.ProteinGroups)[1]
    
    age_cases <- data_scenario %>% filter(Label) %>% pull(age)
    age_controls <- data_scenario %>% filter(!Label) %>% pull(age)
    
    case_values <- data_scenario %>% filter(Label) %>% pull(Quantity)
    control_values <- data_scenario %>% filter(!Label) %>% pull(Quantity)
    
    # Handle missing values
    if (missing_strategy == "zero_assigned") {
      total_cases_valid <- sum(case_values != 0, na.rm = TRUE)
      total_controls_valid <- sum(control_values != 0, na.rm = TRUE)
      missing_cases <- sum(case_values == 0, na.rm = TRUE)
      missing_controls <- sum(control_values == 0, na.rm = TRUE)
    } else {
      total_cases_valid <- sum(!is.na(case_values))
      total_controls_valid <- sum(!is.na(control_values))
      missing_cases <- sum(is.na(case_values))
      missing_controls <- sum(is.na(control_values))
    }
    
    total_cases <- length(case_values)
    total_controls <- length(control_values)
    
    # Calculate demographics
    demo_stats <- calculate_demographics(data_scenario, cancer, "Quantity", missing_strategy,
                                         total_cases_valid, total_controls_valid)
    
    # Statistical analysis
    if (total_cases_valid > 1 && total_controls_valid > 1) {
      roc_metrics <- calculate_roc_metrics(data_scenario)
      fold_change <- calc_fold_change(case_values, control_values)
      
      p_value_fold <- if (total_cases_valid >= 1 && total_controls_valid >= 1) {
        tryCatch(t.test(case_values, control_values)$p.value, error = function(e) NA)
      } else NA
      
      result <- create_peptide_result_list(cancer, protein, Protein_group, scenario, stage, 
                                           peptide, PTM, demo_stats, roc_metrics, fold_change,
                                           p_value_fold, case_values, control_values,
                                           total_cases_valid, total_controls_valid,
                                           missing_cases, missing_controls, total_cases, total_controls)
    } else {
      result <- create_empty_peptide_result_list(cancer, protein, Protein_group, scenario, 
                                                 stage, peptide, PTM, demo_stats, case_values, 
                                                 control_values, total_cases_valid, total_controls_valid,
                                                 missing_cases, missing_controls, total_cases, total_controls)
    }
    
    results[[length(results) + 1]] <- result
  }
  
  return(results)
}

# Helper functions for result creation and demographics
calculate_demographics <- function(data_scenario, cancer, protein_col, missing_strategy, 
                                   total_cases_valid, total_controls_valid) {
  # Age statistics
  age_cases <- data_scenario %>% filter(Label) %>% pull(age)
  age_controls <- data_scenario %>% filter(!Label) %>% pull(age)
  
  # Cigarette statistics
  if (protein_col == "Quantity") {
    filter_condition <- !is.na(Quantity) if missing_strategy != "zero_assigned" else Quantity != 0
    cig_cases <- data_scenario %>% filter(Label & !!filter_condition) %>% pull(cig_stat)
    cig_controls <- data_scenario %>% filter(!Label & !!filter_condition) %>% pull(cig_stat)
  } else {
    filter_condition <- eval(parse(text = paste("!is.na(", protein_col, ") | ", protein_col, "!= 0")))
    cig_cases <- data_scenario %>% filter(Label & eval(filter_condition)) %>% pull(cig_stat)
    cig_controls <- data_scenario %>% filter(!Label & eval(filter_condition)) %>% pull(cig_stat)
  }
  
  # Breast cancer subtypes
  subtype_stats <- if (cancer == "Breast") {
    calculate_breast_subtypes(data_scenario, protein_col, missing_strategy, total_cases_valid)
  } else {
    list(Unkwon = NA, Unkwon_percentage = NA, Luminal = NA, Luminal_percentage = NA,
         HER2 = NA, HER2_percentage = NA, TNBC = NA, TNBC_percentage = NA,
         Other = NA, Other_percentage = NA)
  }
  
  return(list(
    mean_age_cases = mean(age_cases, na.rm = TRUE),
    mean_age_controls = mean(age_controls, na.rm = TRUE),
    cig_stats = calculate_smoking_stats(cig_cases, cig_controls, total_cases_valid, total_controls_valid),
    subtype_stats = subtype_stats
  ))
}

calculate_smoking_stats <- function(cig_cases, cig_controls, total_cases_valid, total_controls_valid) {
  list(
    never_smoked_cases_count = sum(cig_cases == 0),
    never_smoked_controls_count = sum(cig_controls == 0),
    never_smoked_cases_percentage = sum(cig_cases == 0) / total_cases_valid * 100,
    never_smoked_controls_percentage = sum(cig_controls == 0) / total_controls_valid * 100,
    
    current_smoker_cases_count = sum(cig_cases == 1),
    current_smoker_controls_count = sum(cig_controls == 1),
    current_smoker_cases_percentage = sum(cig_cases == 1) / total_cases_valid * 100,
    current_smoker_controls_percentage = sum(cig_controls == 1) / total_controls_valid * 100,
    
    former_smoker_cases_count = sum(cig_cases == 2),
    former_smoker_controls_count = sum(cig_controls == 2),
    former_smoker_cases_percentage = sum(cig_cases == 2) / total_cases_valid * 100,
    former_smoker_controls_percentage = sum(cig_controls == 2) / total_controls_valid * 100
  )
}

calculate_breast_subtypes <- function(data_scenario, protein_col, missing_strategy, total_cases_valid) {
  if (protein_col == "Quantity") {
    filter_condition <- if (missing_strategy != "zero_assigned") !is.na(Quantity) else Quantity != 0
  } else {
    filter_condition <- eval(parse(text = paste("!is.na(", protein_col, ") | ", protein_col, "!= 0")))
  }
  
  data_scenario <- data_scenario %>%
    mutate(
      Subtype = ifelse(Label & eval(filter_condition),
                       mapply(determine_subtype, j_breast_er_status, 
                              j_breast_pr_status, j_breast_her2summ),
                       NA)
    )
  
  list(
    Unkwon = sum(data_scenario$Subtype == "Unknown", na.rm = TRUE),
    Unkwon_percentage = sum(data_scenario$Subtype == "Unknown", na.rm = TRUE) / total_cases_valid * 100,
    Luminal = sum(data_scenario$Subtype == "Luminal", na.rm = TRUE),
    Luminal_percentage = sum(data_scenario$Subtype == "Luminal", na.rm = TRUE) / total_cases_valid * 100,
    HER2 = sum(data_scenario$Subtype == "HER2", na.rm = TRUE),
    HER2_percentage = sum(data_scenario$Subtype == "HER2", na.rm = TRUE) / total_cases_valid * 100,
    TNBC = sum(data_scenario$Subtype == "TNBC", na.rm = TRUE),
    TNBC_percentage = sum(data_scenario$Subtype == "TNBC", na.rm = TRUE) / total_cases_valid * 100,
    Other = sum(data_scenario$Subtype == "Other", na.rm = TRUE),
    Other_percentage = sum(data_scenario$Subtype == "Other", na.rm = TRUE) / total_cases_valid * 100
  )
}

# Result list creation functions
create_result_list <- function(cancer, protein, scenario, stage, demo_stats, roc_metrics,
                                fold_change, p_value_fold, case_values, control_values,
                                total_cases_valid, total_controls_valid, missing_cases, 
                                missing_controls, total_cases, total_controls) {
  list(
    Cancer = cancer, Protein = protein, Scenario = scenario, Stage = stage,
    Mean_Age_Cases = demo_stats$mean_age_cases,
    Mean_Age_Controls = demo_stats$mean_age_controls,
    Cases = total_cases_valid, Controls = total_controls_valid,
    Missing_Cases = missing_cases, Missing_Controls = missing_controls,
    Missing_Cases_Percentage = (missing_cases / total_cases) * 100,
    Missing_Controls_Percentage = (missing_controls / total_controls) * 100,
    
    # Breast cancer subtypes
    BrC_Unkwon = demo_stats$subtype_stats$Unkwon,
    BrC_Unkwon_percentage = demo_stats$subtype_stats$Unkwon_percentage,
    BrC_Luminal = demo_stats$subtype_stats$Luminal,
    BrC_Luminal_percentage = demo_stats$subtype_stats$Luminal_percentage,
    BrC_HER2 = demo_stats$subtype_stats$HER2,
    BrC_HER2_percentage = demo_stats$subtype_stats$HER2_percentage,
    BrC_TNBC = demo_stats$subtype_stats$TNBC,
    BrC_TNBC_percentage = demo_stats$subtype_stats$TNBC_percentage,
    BrC_Other = demo_stats$subtype_stats$Other,
    BrC_Other_percentage = demo_stats$subtype_stats$Other_percentage,
    
    # Smoking statistics
    Never_Smoked_cig_cases_count = demo_stats$cig_stats$never_smoked_cases_count,
    Never_Smoked_cig_controls_count = demo_stats$cig_stats$never_smoked_controls_count,
    Never_Smoked_cig_cases_percentage = demo_stats$cig_stats$never_smoked_cases_percentage,
    Never_Smoked_cig_controls_percentage = demo_stats$cig_stats$never_smoked_controls_percentage,
    Current_Smoker_cig_cases_count = demo_stats$cig_stats$current_smoker_cases_count,
    Current_Smoker_cig_controls_count = demo_stats$cig_stats$current_smoker_controls_count,
    Current_Smoker_cig_cases_percentage = demo_stats$cig_stats$current_smoker_cases_percentage,
    Current_Smoker_cig_controls_percentage = demo_stats$cig_stats$current_smoker_controls_percentage,
    Former_Smoker_cig_cases_count = demo_stats$cig_stats$former_smoker_cases_count,
    Former_Smoker_cig_controls_count = demo_stats$cig_stats$former_smoker_controls_count,
    Former_Smoker_cig_cases_percentage = demo_stats$cig_stats$former_smoker_cases_percentage,
    Former_Smoker_cig_controls_percentage = demo_stats$cig_stats$former_smoker_controls_percentage,
    
    # ROC and performance metrics
    AUC = roc_metrics$AUC,
    AUC_CI_L = roc_metrics$AUC_CI_L,
    AUC_CI_H = roc_metrics$AUC_CI_H,
    AUC_p_value = roc_metrics$p_value,
    
    # 95% cutoff metrics
    Sensitivity_At_95_Specificity = roc_metrics$Sensitivity_At_95_Specificity,
    Cutoff_At_95_Specificity = roc_metrics$Cutoff_At_95_Specificity,
    Specificity_At_95_Sensitivity = roc_metrics$Specificity_At_95_Sensitivity,
    Cutoff_At_95_Sensitivity = roc_metrics$Cutoff_At_95_Sensitivity,
    TP_At_95_Specificity = roc_metrics$TP_At_95_Specificity,
    TN_At_95_Specificity = roc_metrics$TN_At_95_Specificity,
    FP_At_95_Specificity = roc_metrics$FP_At_95_Specificity,
    FN_At_95_Specificity = roc_metrics$FN_At_95_Specificity,
    Accuracy_At_95_Specificity = roc_metrics$Accuracy_At_95_Specificity,
    Precision_At_95_Specificity = roc_metrics$Precision_At_95_Specificity,
    TP_At_95_Sensitivity = roc_metrics$TP_At_95_Sensitivity,
    TN_At_95_Sensitivity = roc_metrics$TN_At_95_Sensitivity,
    FP_At_95_Sensitivity = roc_metrics$FP_At_95_Sensitivity,
    FN_At_95_Sensitivity = roc_metrics$FN_At_95_Sensitivity,
    Accuracy_At_95_Sensitivity = roc_metrics$Accuracy_At_95_Sensitivity,
    Precision_At_95_Sensitivity = roc_metrics$Precision_At_95_Sensitivity,
    
    # 86% cutoff metrics
    Sensitivity_At_86_Specificity = roc_metrics$Sensitivity_At_86_Specificity,
    Cutoff_At_86_Specificity = roc_metrics$Cutoff_At_86_Specificity,
    Specificity_At_86_Sensitivity = roc_metrics$Specificity_At_86_Sensitivity,
    Cutoff_At_86_Sensitivity = roc_metrics$Cutoff_At_86_Sensitivity,
    TP_At_86_Specificity = roc_metrics$TP_At_86_Specificity,
    TN_At_86_Specificity = roc_metrics$TN_At_86_Specificity,
    FP_At_86_Specificity = roc_metrics$FP_At_86_Specificity,
    FN_At_86_Specificity = roc_metrics$FN_At_86_Specificity,
    Accuracy_At_86_Specificity = roc_metrics$Accuracy_At_86_Specificity,
    Precision_At_86_Specificity = roc_metrics$Precision_At_86_Specificity,
    TP_At_86_Sensitivity = roc_metrics$TP_At_86_Sensitivity,
    TN_At_86_Sensitivity = roc_metrics$TN_At_86_Sensitivity,
    FP_At_86_Sensitivity = roc_metrics$FP_At_86_Sensitivity,
    FN_At_86_Sensitivity = roc_metrics$FN_At_86_Sensitivity,
    Accuracy_At_86_Sensitivity = roc_metrics$Accuracy_At_86_Sensitivity,
    Precision_At_86_Sensitivity = roc_metrics$Precision_At_86_Sensitivity,
    
    # Fold change and means
    Fold_Change = fold_change,
    Fold_Change_p_value = p_value_fold,
    Average_case_intensities = mean(case_values, na.rm = TRUE),
    Average_control_intensities = mean(control_values, na.rm = TRUE)
  )
}

create_empty_result_list <- function(cancer, protein, scenario, stage, demo_stats,
                                     case_values, control_values, total_cases_valid, 
                                     total_controls_valid, missing_cases, missing_controls,
                                     total_cases, total_controls) {
  result <- create_result_list(cancer, protein, scenario, stage, demo_stats,
                               list(AUC = NA, AUC_CI_L = NA, AUC_CI_H = NA, p_value = NA,
                                    Sensitivity_At_95_Specificity = NA, Cutoff_At_95_Specificity = NA,
                                    Specificity_At_95_Sensitivity = NA, Cutoff_At_95_Sensitivity = NA,
                                    Sensitivity_At_86_Specificity = NA, Cutoff_At_86_Specificity = NA,
                                    Specificity_At_86_Sensitivity = NA, Cutoff_At_86_Sensitivity = NA,
                                    TP_At_95_Specificity = NA, TN_At_95_Specificity = NA,
                                    FP_At_95_Specificity = NA, FN_At_95_Specificity = NA,
                                    Accuracy_At_95_Specificity = NA, Precision_At_95_Specificity = NA,
                                    TP_At_95_Sensitivity = NA, TN_At_95_Sensitivity = NA,
                                    FP_At_95_Sensitivity = NA, FN_At_95_Sensitivity = NA,
                                    Accuracy_At_95_Sensitivity = NA, Precision_At_95_Sensitivity = NA,
                                    TP_At_86_Specificity = NA, TN_At_86_Specificity = NA,
                                    FP_At_86_Specificity = NA, FN_At_86_Specificity = NA,
                                    Accuracy_At_86_Specificity = NA, Precision_At_86_Specificity = NA,
                                    TP_At_86_Sensitivity = NA, TN_At_86_Sensitivity = NA,
                                    FP_At_86_Sensitivity = NA, FN_At_86_Sensitivity = NA,
                                    Accuracy_At_86_Sensitivity = NA, Precision_At_86_Sensitivity = NA),
                               NA, NA, case_values, control_values, total_cases_valid, 
                               total_controls_valid, missing_cases, missing_controls, 
                               total_cases, total_controls)
  return(result)
}

create_peptide_result_list <- function(cancer, protein, protein_group, scenario, stage, peptide, 
                                       ptm, demo_stats, roc_metrics, fold_change, p_value_fold,
                                       case_values, control_values, total_cases_valid, 
                                       total_controls_valid, missing_cases, missing_controls,
                                       total_cases, total_controls) {
  result <- create_result_list(cancer, protein, scenario, stage, demo_stats, roc_metrics,
                               fold_change, p_value_fold, case_values, control_values,
                               total_cases_valid, total_controls_valid, missing_cases,
                               missing_controls, total_cases, total_controls)
  
  # Add peptide-specific fields
  result$Protein_Group <- protein_group
  result$Peptide <- peptide
  result$PTM_Flag_most_freq <- ptm
  result$Average_case_values <- mean(case_values, na.rm = TRUE)
  result$Average_control_values <- mean(control_values, na.rm = TRUE)
  result$SD_case_values <- sd(case_values, na.rm = TRUE)
  result$SD_control_values <- sd(control_values, na.rm = TRUE)
  
  # Remove protein-specific field
  result$Average_case_intensities <- NULL
  result$Average_control_intensities <- NULL
  
  return(result)
}

create_empty_peptide_result_list <- function(cancer, protein, protein_group, scenario, stage, 
                                             peptide, ptm, demo_stats, case_values, control_values,
                                             total_cases_valid, total_controls_valid, 
                                             missing_cases, missing_controls, total_cases, total_controls) {
  result <- create_empty_result_list(cancer, protein, scenario, stage, demo_stats,
                                     case_values, control_values, total_cases_valid,
                                     total_controls_valid, missing_cases, missing_controls,
                                     total_cases, total_controls)
  
  # Add peptide-specific fields
  result$Protein_Group <- protein_group
  result$Peptide <- peptide
  result$PTM_Flag_most_freq <- ptm
  result$Average_case_values <- mean(case_values, na.rm = TRUE)
  result$Average_control_values <- mean(control_values, na.rm = TRUE)
  result$SD_case_values <- sd(case_values, na.rm = TRUE)
  result$SD_control_values <- sd(control_values, na.rm = TRUE)
  
  # Remove protein-specific field
  result$Average_case_intensities <- NULL
  result$Average_control_intensities <- NULL
  
  return(result)
}

# Main Execution Function ----

#' Run comprehensive PLCO proteomics analysis
#' @param analysis_level Either "protein" or "peptide"
#' @param data_types Vector of data types to analyze (e.g., c("IgB", "FT", "SUM"))
#' @param missing_strategies Vector of missing value strategies
#' @param cancer_type Specific cancer type to analyze (optional, defaults to all)
#' @param use_parallel Whether to use parallel processing
run_plco_analysis <- function(analysis_level = "protein", 
                              data_types = c("IgB", "FT", "SUM"),
                              missing_strategies = c("zero_assigned", "NA_assigned"),
                              cancer_type = NULL,
                              use_parallel = TRUE) {
  
  cat("Starting PLCO", analysis_level, "analysis...\n")
  cat("Data types:", paste(data_types, collapse = ", "), "\n")
  cat("Missing value strategies:", paste(missing_strategies, collapse = ", "), "\n")
  
  # Define file paths
  files <- list(
    IgB = if (analysis_level == "protein") {
      "protein_IgB_PLCO.csv"
    } else {
      "peptide_IgB_PLCO.csv"
    },
    FT = if (analysis_level == "protein") {
      "protein_FT_PLCO.csv"  
    } else {
      "peptide_FT_PLCO.csv"
    }
  )
  
  # Cancer types to analyze
  cancers_to_analyze <- if (is.null(cancer_type)) CONFIG$cancer_types else cancer_type
  
  # Loop through missing value strategies
  for (missing_strategy in missing_strategies) {
    cat("\n=== Processing with", missing_strategy, "strategy ===\n")
    
    # Create output directory
    output_dir <- create_output_dir(CONFIG$output_dir, analysis_level, missing_strategy)
    cat("Output directory:", output_dir, "\n")
    
    # Load and preprocess data
    datasets <- list()
    
    if ("IgB" %in% data_types) {
      datasets$IgB <- load_data(file.path(CONFIG$data_dir, files$IgB), missing_strategy)
      if (analysis_level == "peptide") {
        datasets$IgB$PTM_Flag <- "Not_provided"
      }
    }
    
    if ("FT" %in% data_types) {
      datasets$FT <- load_data(file.path(CONFIG$data_dir, files$FT), missing_strategy)
      if (analysis_level == "peptide") {
        datasets$FT$PTM_Flag <- "Not_provided"
      }
    }
    
    if ("SUM" %in% data_types && analysis_level == "protein" && 
        "IgB" %in% data_types && "FT" %in% data_types) {
      cat("Creating combined dataset (IgB + FT)...\n")
      datasets$SUM <- sum_protein_dfs(datasets$FT, datasets$IgB, CONFIG$protein_column, use_parallel)
    }
    
    # Process each dataset
    for (data_type in names(datasets)) {
      cat("\n--- Processing", data_type, "data ---\n")
      dataset <- datasets[[data_type]]
      
      if (analysis_level == "peptide") {
        dataset <- dataset[!is.na(dataset$EG.PrecursorId), ]
        if (missing_strategy == "NA_assigned") {
          dataset$Quantity[dataset$Quantity == 0] <- NA
        }
      }
      
      # Determine proteins/features to analyze
      if (analysis_level == "protein") {
        if (data_type == "SUM") {
          proteins <- colnames(dataset[, -c(1:(CONFIG$protein_column_sum - 1))])
        } else {
          proteins <- colnames(dataset[, -c(1:(CONFIG$protein_column - 1))])
        }
      } else {
        proteins <- unique(dataset$PG.ProteinNames)
      }
      
      cat("Analyzing", length(proteins), analysis_level, "features...\n")
      
      # Create parameter combinations
      params <- expand.grid(
        cancer = cancers_to_analyze,
        protein = proteins,
        scenario = CONFIG$scenarios,
        stringsAsFactors = FALSE
      )
      
      # Add stage information
      lookup_table <- data.frame(
        scenario = seq(1, 10),
        stage = c("All", "All", rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))
      )
      params$stage <- lookup_table$stage[match(params$scenario, lookup_table$scenario)]
      
      cat("Processing", nrow(params), "parameter combinations...\n")
      
      # Run analysis
      if (use_parallel && nrow(params) > 100) {
        cat("Using parallel processing...\n")
        cl <- makeCluster(parallel::detectCores() - 1)
        clusterEvalQ(cl, {
          library(dplyr)
          library(pROC)
        })
        clusterExport(cl, c("dataset", "missing_strategy", "analysis_level",
                            ls(pattern = "^process_|^calculate_|^create_|^determine_"),
                            "CONFIG"))
        
        if (analysis_level == "protein") {
          results <- parLapply(cl, 1:nrow(params), function(i) {
            process_protein_combination(dataset, params[i, ], missing_strategy)
          })
        } else {
          results <- parLapply(cl, 1:nrow(params), function(i) {
            process_peptide_combination(dataset, params[i, ], missing_strategy)
          })
        }
        
        stopCluster(cl)
      } else {
        if (analysis_level == "protein") {
          results <- lapply(1:nrow(params), function(i) {
            if (i %% 1000 == 0) cat("Processed", i, "of", nrow(params), "combinations\n")
            process_protein_combination(dataset, params[i, ], missing_strategy)
          })
        } else {
          results <- lapply(1:nrow(params), function(i) {
            if (i %% 100 == 0) cat("Processed", i, "of", nrow(params), "combinations\n")
            process_peptide_combination(dataset, params[i, ], missing_strategy)
          })
        }
      }
      
      # Process results
      if (analysis_level == "peptide") {
        # Flatten peptide results
        results_flat <- do.call(c, results)
        df_results_final <- do.call(rbind, lapply(results_flat, as.data.frame))
      } else {
        # Process protein results
        df_2 <- list()
        for (i in seq_along(results)) {
          df <- results[[i]]
          if (is.list(df$AUC)) df$AUC <- df$AUC[[1]]
          df_2[[i]] <- as.data.frame(df)
        }
        df_results_final <- do.call(rbind, df_2)
      }
      
      # Add derived metrics
      if (nrow(df_results_final) > 0) {
        df_results_final$neg_log10_pval_FC <- -log10(df_results_final$Fold_Change_p_value)
        df_results_final$log2_FC <- log2(df_results_final$Fold_Change)
        
        # Save results
        output_prefix <- if (is.null(cancer_type)) {
          paste0("results_", analysis_level, "_", data_type, "_", missing_strategy)
        } else {
          paste0(cancer_type, "_", analysis_level, "_results_", data_type, "_", missing_strategy)
        }
        
        # Save as RData
        save(df_results_final, file = file.path(output_dir, paste0(output_prefix, ".RData")))
        
        # Save as Excel (split if too large)
        if (nrow(df_results_final) > 1e6) {
          split_data <- split(df_results_final, (seq(nrow(df_results_final)) - 1) %/% 1e6)
          lapply(seq_along(split_data), function(i) {
            write_xlsx(split_data[[i]], file.path(output_dir, paste0(output_prefix, "_part_", i, ".xlsx")))
          })
        } else {
          write_xlsx(df_results_final, file.path(output_dir, paste0(output_prefix, ".xlsx")))
        }
        
        cat("Results saved:", nrow(df_results_final), "rows\n")
      } else {
        cat("Warning: No results generated for", data_type, "\n")
      }
    }
  }
  
  cat("\n=== Analysis completed successfully! ===\n")
}

# Example usage and main execution ----
if (interactive()) {
  cat("PLCO Proteomics Analysis Pipeline Loaded Successfully!\n")
  cat("Use run_plco_analysis() to start analysis.\n")
  cat("Example usage:\n")
  cat('run_plco_analysis(analysis_level = "protein", data_types = c("IgB", "FT", "SUM"))\n')
  cat('run_plco_analysis(analysis_level = "peptide", data_types = "IgB", cancer_type = "Breast")\n')
} else {
  # Command line execution
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    cancer_type <- args[1]
    analysis_level <- if (length(args) > 1) args[2] else "protein"
    data_type <- if (length(args) > 2) args[3] else "IgB"
    
    run_plco_analysis(
      analysis_level = analysis_level,
      data_types = data_type,
      cancer_type = cancer_type
    )
  } else {
    # Run full analysis by default
    run_plco_analysis()
  }
}
