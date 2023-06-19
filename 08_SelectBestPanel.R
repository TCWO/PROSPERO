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

# Load data --------------------------------------------------------

df_stats_amg <- read.csv2(file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_stats_amg.csv'), stringsAsFactors = FALSE)
df_panel_data_amg <- read.csv2(file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_panel_data_amg.csv'), stringsAsFactors = FALSE)
df_stats_rt <- read.csv2(file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_stats_rt.csv'), stringsAsFactors = FALSE)
df_panel_data_rt <- read.csv2(file.path('..', 'Data', 'CellLines', 'PanelOptimization', 'df_panel_data_rt.csv'), stringsAsFactors = FALSE)
df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'df_data_dichotomous.csv'), stringsAsFactors = FALSE)
drug_treatment_meta <- read.csv2(file.path('..', 'Data', 'CellLines', 'drug_treatment_metadata_amg.csv'), stringsAsFactors = FALSE)

# AMG best panel: classification based on MDM2 and then regression based on p21 --------------------------------------------------------

best_amg <- df_data %>% 
  filter(treatment != 'RT') %>%
  filter(marker %in% c('MDM2', 'p21')) %>%
  group_by(sample.id, treatment, marker, induction) %>% 
  summarise(N = n()) %>% 
  ungroup() %>% 
  group_by(sample.id, treatment, marker) %>% 
  mutate(M = sum(N)) %>% 
  ungroup() %>% 
  mutate(P = 100*N/M) %>% 
  select(-N, -M) %>% 
  filter(induction == 1) %>% 
  spread(treatment, P) %>% 
  left_join(drug_treatment_meta) %>% 
  mutate(D = AMG232-Control) %>% 
  select(sample.id, marker, D, response:MaxInhibition) %>% 
  spread(marker, D, fill = 0) %>% 
  group_by(sample.id) %>% 
  sample_n(1) %>% 
  ungroup()

mdm2_induced_samples <- best_amg %>% filter(MDM2 > 15) %>% 
  mutate(predicted.MaxInhibition = mean(MaxInhibition))

p21_model <- best_amg %>% filter(MDM2 <= 15)
tmp_model <- mgcv::gam(MaxInhibition ~ s(p21, bs = 'cs'), data = p21_model, method = 'REML')
p21_model <- p21_model %>% 
  mutate(predicted.MaxInhibition = predict(tmp_model))

best_amg <- mdm2_induced_samples %>% bind_rows(p21_model)

write.csv2(best_amg, file.path('..', 'Data', 'CellLines', 'OptimalPanels', 'amg232_pli_predicted_pli.csv'), row.names = FALSE)

# RT best panel: regression based on p21, pATM, and pH2AX --------------------------------------------------------

best_rt <- df_panel_data_rt %>% 
  filter(panel == "p21, pATM, pH2AX")

write.csv2(best_rt, file.path('..', 'Data', 'CellLines', 'OptimalPanels', 'rt_pli_predicted_pli.csv'), row.names = FALSE)