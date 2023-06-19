# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'parallel', 'doSNOW')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(parallel)
library(doSNOW)

# we have added options for parallelizing this step. However, for simplicity, they are not active.

# Load data --------------------------------------------------------

df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'df_data.csv'), stringsAsFactors = FALSE) %>% 
  mutate(batch = as.character(batch)) %>% 
  mutate(value = asinh(value)) # transform values to hyperbolic asin

# Filter batch control (BC) samples --------------------------------------------------------

df_bc <- df_data %>% 
  dplyr::filter(label == 'BC') %>% 
  mutate(id = paste(sample.id, treatment, time.point, sep = '_'), 
         id2 = paste(sample.id, treatment, time.point, batch, sep = '_'))

# Remove astrocytes --------------------------------------------------------

df_bc <- df_bc %>% dplyr::filter(sample.id != 'Astrocytes')

# Source function for overlapping coefficient --------------------------------------------------------

source('auxiliary_functions.R')

# Estimate overlapping coefficients  --------------------------------------------------------

n_cores <- 10
n_cores_max <- detectCores()
cl <- parallel::makeCluster(min(n_cores, (n_cores_max-1)))
registerDoSNOW(cl)

overlapping_stats <- data.frame()
ref_batch <- 3 # reference batch (largest number of cells)

for(i in unique(df_bc$marker)){
  print(paste('Analysing marker:', i))
  for(k in setdiff(unique(df_bc$batch), ref_batch)){
    tmp_data <- df_bc %>% 
      dplyr::filter(marker == i) 
    
    best_result <- overlappingFunction(tmp_data, k)
    
    overlapping_stats <- overlapping_stats %>% bind_rows(best_result %>% mutate(marker = i, batch = k))
  }
}

stopCluster(cl)

# Apply normalization  --------------------------------------------------------

df_data_norm <- df_data %>%
  left_join(overlapping_stats) %>% 
  mutate(value = ifelse(batch != 3, tmp_alpha+value, value)) %>% 
  select(-tmp_alpha, -OVL, -OVL.ref)

# Write results  --------------------------------------------------------

write.csv2(overlapping_stats, file.path('..', 'Data', 'CellLines', 'overlapping_stats.csv'), row.names = FALSE)
write.csv2(df_data_norm, file.path('..', 'Data', 'CellLines', 'df_data_norm.csv'), row.names = FALSE)
