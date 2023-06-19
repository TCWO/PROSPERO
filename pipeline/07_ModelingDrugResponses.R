# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'mixtools', 'doParallel', 'BiocManager', 'doSNOW', 'parallel', 'minerva', 'caret')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.packages <- c('destiny', 'slingshot')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

library(tidyverse)
library(mixtools)
library(slingshot)
library(destiny)
library(doParallel)
library(doSNOW)
library(parallel)
library(minerva)
library(caret)

# Define functional markers --------------------------------------------------------

functional_markers <- c('p53', 'p21', 'MDM2', 'pH2AX', 'BAX', 'CC3', 'pATM')

# Load data --------------------------------------------------------

df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'data_annotated_tumor_subtypes_batches2and3.csv'), stringsAsFactors = FALSE) %>% 
  bind_rows(read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'data_annotated_tumor_subtypes_batch1.csv'), stringsAsFactors = FALSE)) %>% 
  filter(time.point == 16)
drug_treatment_meta <- read.csv2(file.path('..', 'Data', 'CellLines', 'drug_treatment_metadata_amg.csv'), stringsAsFactors = FALSE)

# Gaussian Mixture Models (GMM) for p53 --------------------------------------------------------

df_gmm <- data.frame()
for(i in unique(df_data$sample.id)){
  print(i)
  tmp_data <- df_data %>% filter(sample.id == i, treatment == 'Control') %>% filter(p53 > min(p53))
  set.seed(1234)
  tmp_gmm <- normalmixEM(tmp_data$p53)
  df_gmm <- df_gmm %>% bind_rows(data.frame(sample.id = i, mu = c(tmp_gmm$mu[1], tmp_gmm$mu[2]), sigma = c(tmp_gmm$sigma[1], tmp_gmm$sigma[2]), G = c('gaussian01', 'gaussian02')))
}

df_gmm <- df_gmm %>% mutate(Q95 = qnorm(0.05, mu, sigma, lower.tail = FALSE))
df_gmm <- df_gmm %>% group_by(sample.id) %>% filter(mu == min(mu)) %>% ungroup()

# Dicotomize signals --------------------------------------------------------

df_data <- df_data %>% 
  gather(marker, value, BAX:Vimentin) %>% 
  filter(marker %in% functional_markers) %>% 
  group_by(sample.id, marker) %>% 
  mutate(M = mean(value[treatment == 'Control' & value > min(value)], na.rm = TRUE),
         S = sd(value[treatment == 'Control' & value > min(value)], na.rm = TRUE),
         p.high = pnorm(value, mean = mean(value[treatment == 'Control' & value > min(value)], na.rm = TRUE), sd = sd(value[treatment == 'Control' & value > min(value)], na.rm = TRUE), lower.tail = FALSE),
         p.low = pnorm(value, mean = mean(value[treatment == 'Control' & value > min(value)], na.rm = TRUE), sd = sd(value[treatment == 'Control' & value > min(value)], na.rm = TRUE), lower.tail = TRUE)) %>% 
  ungroup() %>% 
  mutate(induction = ifelse(p.high < 0.05, 1, 0)) %>% 
  left_join(df_gmm) %>% 
  mutate(induction_p53 = ifelse(value > Q95, 1, 0)) %>% 
  mutate(induction = ifelse(is.na(induction_p53), induction, induction_p53)) %>% 
  select(sample.id, treatment, OID, marker, induction)

write.csv2(df_data, file.path('..', 'Data', 'CellLines', 'df_data_dichotomous.csv'), row.names = FALSE)

# Optimize panel --------------------------------------------------------

## AMG
tmp_amg <- df_data %>% 
  filter(marker %in% c('p21', 'p53', 'MDM2', 'BAX', 'CC3')) %>% 
  filter(treatment != 'RT')

df_panel_data_amg <- data.frame()
df_stats_amg <- data.frame()
n_cores <- 10

for(i in c(1:n_distinct(tmp_amg$marker))){
  print(i)
  tmp_combs <- combn(unique(tmp_amg$marker), i) %>% t()
  
  cl <- parallel::makeCluster(n_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = nrow(tmp_combs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  tmp_list <- foreach(j = c(1:nrow(tmp_combs)), .packages = c('tidyverse', 'slingshot', 'destiny', 'minerva', 'caret'), .options.snow = opts) %dopar% {
    # for(j in c(1:nrow(tmp_combs))){
    tmp_panel <- tmp_combs[j,]
    # print(tmp_panel)
    tmp_panel_data <- tmp_amg %>% 
      filter(marker %in% tmp_panel) %>% 
      group_by(sample.id, treatment, OID) %>% 
      summarise(ID = paste(marker, induction, sep = '_', collapse = '')) %>% 
      ungroup() %>% 
      group_by(sample.id, treatment, ID) %>% 
      summarise(N = n()) %>% 
      ungroup() %>% 
      group_by(sample.id, treatment) %>% 
      mutate(M = sum(N)) %>% 
      ungroup() %>% 
      mutate(P = 100*N/M) %>% 
      select(-N, -M) %>% 
      spread(treatment, P) %>% 
      left_join(drug_treatment_meta) %>% 
      mutate(D = AMG232-Control) %>% 
      select(sample.id, ID, D, response:MaxInhibition) %>% 
      spread(ID, D, fill = 0)
    
    tmp_panel_data <- tmp_panel_data %>% group_by(sample.id) %>% sample_n(1) %>% ungroup()
    
    tmp_rownames <- tmp_panel_data$sample.id
    tmp_matrix <- tmp_panel_data %>% select(-sample.id, -response, -p53.mut, -AUC, -MaxInhibition) %>% as.matrix()
    rownames(tmp_matrix) <- tmp_rownames
    
    dm <- DiffusionMap(tmp_matrix)
    
    tmp_panel_data <- tmp_panel_data %>%
      mutate(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3)
    
    # getLineages crashes if sample size small. Increase virtually repeating samples.
    tmp_panel_data <- tmp_panel_data %>% mutate(QC = 1) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0)) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0)) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0))
    
    cl <- tmp_panel_data$response %>% as.factor()
    tmp_lineages <- getLineages(data = tmp_panel_data %>% dplyr::select(DC1,DC2) %>% as.matrix(), clusterLabels = cl)
    tmp_lineages
    
    tmp_curves <- as.SlingshotDataSet(getCurves(tmp_lineages))
    
    tmp_panel_data <- tmp_panel_data %>% 
      mutate(curve.s1 = tmp_curves@curves$Lineage1$s[,1], 
             curve.s2 = tmp_curves@curves$Lineage1$s[,2], 
             order = tmp_curves@curves$Lineage1$ord,
             pseudotime = slingshot::slingPseudotime(tmp_curves))
    
    tmp_panel_data <- tmp_panel_data %>% filter(QC == 1) %>% select(-QC)
    
    tmp_model_max_inhib <- lm(tmp_panel_data$MaxInhibition~tmp_panel_data$pseudotime)
    tmp_predictions_max_inhib <- tmp_model_max_inhib %>% predict(tmp_panel_data, type = 'response')
    
    tmp_model_auc <- lm(tmp_panel_data$AUC~tmp_panel_data$pseudotime)
    tmp_predictions_auc <- tmp_model_auc %>% predict(tmp_panel_data, type = 'response')
    
    tmp_panel_data <- tmp_panel_data %>% mutate(predicted.MaxInhibition = tmp_predictions_max_inhib, predicted.AUC = tmp_predictions_auc)
    
    tmp_cor <- cor.test(tmp_panel_data$AUC, tmp_panel_data$pseudotime)
    tmp_mic <- minerva::mine(tmp_panel_data$AUC, tmp_panel_data$pseudotime)$MIC
    tmp_stats_auc <- data.frame(panel = paste(tmp_panel, collapse = '', sep = '_'), 
                                Nmarkers = i, parameter = 'AUC', R = tmp_cor$estimate, pVal = tmp_cor$p.value, 
                                MIC = tmp_mic, RMSE = RMSE(tmp_predictions_auc, tmp_panel_data$AUC), R2 = R2(tmp_predictions_auc, tmp_panel_data$AUC))
    
    tmp_cor <- cor.test(tmp_panel_data$MaxInhibition, tmp_panel_data$pseudotime)
    tmp_mic <- minerva::mine(tmp_panel_data$MaxInhibition, tmp_panel_data$pseudotime)$MIC
    tmp_stats_maxinhib <- data.frame(panel = paste(tmp_panel, collapse = '', sep = '_'), Nmarkers = i, parameter = 'MaxInhibition', R = tmp_cor$estimate, pVal = tmp_cor$p.value, MIC = tmp_mic, RMSE = RMSE(tmp_predictions_max_inhib, tmp_panel_data$MaxInhibition), R2 = R2(tmp_predictions_max_inhib, tmp_panel_data$MaxInhibition))
    
    tmp_stats <- tmp_stats_auc %>% bind_rows(tmp_stats_maxinhib)
    tmp_panel_data <- tmp_panel_data %>% gather(feature, value, -sample.id, -response, -p53.mut, -AUC, -predicted.AUC, -MaxInhibition, -predicted.MaxInhibition, -DC1, -DC2, -DC3, -curve.s1, -curve.s2, -order, -pseudotime) %>% mutate(panel = paste(tmp_panel, collapse = ', '))
    return(list(tmp_stats, tmp_panel_data))
  }
  stopCluster(cl)
  close(pb)
  
  # Extract the first data frame from each nested list
  first_elements <- lapply(tmp_list, function(x) x[[1]])
  # Bind the rows of the first elements
  tmp_stats <- do.call(rbind, first_elements)
  
  # Extract the second data frame from each nested list
  second_elements <- lapply(tmp_list, function(x) x[[2]])
  # Bind the rows of the first elements
  tmp_panel_data <- do.call(rbind, second_elements)
  
  df_stats_amg <- df_stats_amg %>% bind_rows(tmp_stats)
  df_panel_data_amg <- df_panel_data_amg %>% bind_rows(tmp_panel_data)
}

## RT
tmp_rt <- df_data %>% 
  filter(marker %in% c('p21', 'p53', 'MDM2', 'BAX', 'CC3', 'pATM', 'pH2AX')) %>% 
  filter(treatment != 'AMG232')

df_panel_data_rt <- data.frame()
df_stats_rt <- data.frame()

for(i in c(1:n_distinct(tmp_rt$marker))){
  print(i)
  tmp_combs <- combn(unique(tmp_rt$marker), i) %>% t()
  
  cl <- parallel::makeCluster(n_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = nrow(tmp_combs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  tmp_list <- foreach(j = c(1:nrow(tmp_combs)), .packages = c('tidyverse', 'slingshot', 'destiny', 'minerva', 'caret'), .options.snow = opts) %dopar% {
    # for(j in c(1:nrow(tmp_combs))){
    tmp_panel <- tmp_combs[j,]
    # print(tmp_panel)
    tmp_panel_data <- tmp_rt %>% 
      filter(marker %in% tmp_panel) %>% 
      group_by(sample.id, treatment, OID) %>% 
      summarise(ID = paste(marker, induction, sep = '_', collapse = '')) %>% 
      ungroup() %>% 
      group_by(sample.id, treatment, ID) %>% 
      summarise(N = n()) %>% 
      ungroup() %>% 
      group_by(sample.id, treatment) %>% 
      mutate(M = sum(N)) %>% 
      ungroup() %>% 
      mutate(P = 100*N/M) %>% 
      select(-N, -M) %>% 
      spread(treatment, P) %>% 
      left_join(drug_treatment_meta) %>% 
      mutate(D = RT-Control) %>% 
      select(sample.id, ID, D, response:MaxInhibition) %>% 
      spread(ID, D, fill = 0)
    
    tmp_panel_data <- tmp_panel_data %>% group_by(sample.id) %>% sample_n(1) %>% ungroup()
    
    tmp_rownames <- tmp_panel_data$sample.id
    tmp_matrix <- tmp_panel_data %>% select(-sample.id, -response, -p53.mut, -AUC, -MaxInhibition) %>% as.matrix()
    rownames(tmp_matrix) <- tmp_rownames
    
    dm <- DiffusionMap(tmp_matrix)
    
    tmp_panel_data <- tmp_panel_data %>%
      mutate(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3)
    
    # getLineages crashes if sample size small. Increase virtually repeating samples.
    tmp_panel_data <- tmp_panel_data %>% mutate(QC = 1) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0)) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0)) %>% 
      bind_rows(tmp_panel_data %>% mutate(QC = 0))
    
    cl <- tmp_panel_data$response %>% as.factor()
    tmp_lineages <- getLineages(data = tmp_panel_data %>% dplyr::select(DC1,DC2) %>% as.matrix(), clusterLabels = cl)
    tmp_lineages
    
    tmp_curves <- as.SlingshotDataSet(getCurves(tmp_lineages))
    
    tmp_panel_data <- tmp_panel_data %>% 
      mutate(curve.s1 = tmp_curves@curves$Lineage1$s[,1], 
             curve.s2 = tmp_curves@curves$Lineage1$s[,2], 
             order = tmp_curves@curves$Lineage1$ord,
             pseudotime = slingshot::slingPseudotime(tmp_curves))
    
    tmp_panel_data <- tmp_panel_data %>% filter(QC == 1) %>% select(-QC)
    
    tmp_model_max_inhib <- lm(tmp_panel_data$MaxInhibition~tmp_panel_data$pseudotime)
    tmp_predictions_max_inhib <- tmp_model_max_inhib %>% predict(tmp_panel_data, type = 'response')
    
    tmp_model_auc <- lm(tmp_panel_data$AUC~tmp_panel_data$pseudotime)
    tmp_predictions_auc <- tmp_model_auc %>% predict(tmp_panel_data, type = 'response')
    
    tmp_panel_data <- tmp_panel_data %>% mutate(predicted.MaxInhibition = tmp_predictions_max_inhib, predicted.AUC = tmp_predictions_auc)
    
    tmp_cor <- cor.test(tmp_panel_data$AUC, tmp_panel_data$pseudotime)
    tmp_mic <- minerva::mine(tmp_panel_data$AUC, tmp_panel_data$pseudotime)$MIC
    tmp_stats_auc <- data.frame(panel = paste(tmp_panel, collapse = '', sep = '_'), 
                                Nmarkers = i, parameter = 'AUC', R = tmp_cor$estimate, pVal = tmp_cor$p.value, 
                                MIC = tmp_mic, RMSE = RMSE(tmp_predictions_auc, tmp_panel_data$AUC), R2 = R2(tmp_predictions_auc, tmp_panel_data$AUC))
    
    tmp_cor <- cor.test(tmp_panel_data$MaxInhibition, tmp_panel_data$pseudotime)
    tmp_mic <- minerva::mine(tmp_panel_data$MaxInhibition, tmp_panel_data$pseudotime)$MIC
    tmp_stats_maxinhib <- data.frame(panel = paste(tmp_panel, collapse = '', sep = '_'), Nmarkers = i, parameter = 'MaxInhibition', R = tmp_cor$estimate, pVal = tmp_cor$p.value, MIC = tmp_mic, RMSE = RMSE(tmp_predictions_max_inhib, tmp_panel_data$MaxInhibition), R2 = R2(tmp_predictions_max_inhib, tmp_panel_data$MaxInhibition))
    
    tmp_stats <- tmp_stats_auc %>% bind_rows(tmp_stats_maxinhib)
    tmp_panel_data <- tmp_panel_data %>% gather(feature, value, -sample.id, -response, -p53.mut, -AUC, -predicted.AUC, -MaxInhibition, -predicted.MaxInhibition, -DC1, -DC2, -DC3, -curve.s1, -curve.s2, -order, -pseudotime) %>% mutate(panel = paste(tmp_panel, collapse = ', '))
    return(list(tmp_stats, tmp_panel_data))
  }
  stopCluster(cl)
  close(pb)
  
  # Extract the first data frame from each nested list
  first_elements <- lapply(tmp_list, function(x) x[[1]])
  # Bind the rows of the first elements
  tmp_stats <- do.call(rbind, first_elements)
  
  # Extract the second data frame from each nested list
  second_elements <- lapply(tmp_list, function(x) x[[2]])
  # Bind the rows of the first elements
  tmp_panel_data <- do.call(rbind, second_elements)
  
  df_stats_rt <- df_stats_rt %>% bind_rows(tmp_stats)
  df_panel_rt <- df_panel_data_rt %>% bind_rows(tmp_panel_data)
}

# Write data --------------------------------------------------------

write.csv2(df_stats_amg, file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_stats_amg.csv'), row.names = FALSE)
write.csv2(df_panel_data_amg, file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_panel_data_amg.csv'), row.names = FALSE)
write.csv2(df_stats_rt, file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_stats_rt.csv'), row.names = FALSE)
write.csv2(df_panel_rt, file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_panel_data_rt.csv'), row.names = FALSE)

