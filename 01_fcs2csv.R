# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'BiocManager')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

bioconductor.packages <- c('flowCore')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install("flowCore")

library(flowCore)
library(tidyverse)

# Load, merge, and tabulate data files --------------------------------------------------------

input_directory <- file.path('..', 'Data', 'CellLines', 'Pathsetter_cleaned')
tmp_files <- list.files(input_directory, pattern = 'fcs')

df_data <- lapply(tmp_files, function(x){
  tmp_data <- read.FCS(file.path(input_directory, x), transformation = FALSE, emptyValue = FALSE)
  tmp_exp <- tmp_data@exprs %>% as.data.frame() %>% mutate(file = x)
}) %>% bind_rows() %>%
  # mutate(ofile = file) %>% 
  mutate(file = gsub('-', '', file)) %>% 
  mutate(file = sub('20_80', '20.80', file)) %>% 
  mutate(file = sub('50_50', '50.50', file)) %>% 
  mutate(file = sub('80_20', '80.20', file)) %>% 
  group_by(file) %>% 
  mutate(OID = 1:n()) %>% 
  ungroup() %>% 
  separate(file, c('sample.id', 'treatment', 'time.point', 'label', 'batch'), sep = '_') %>% # tabulation of sample information based on file names 
  mutate(batch = sub('.fcs', '', batch))

# add marker information based on metadata 
df_parameters <- lapply(tmp_files, function(x){
  tmp_data <- read.FCS(file.path(input_directory, x), transformation = FALSE, emptyValue = FALSE)
  tmp_exp <- tmp_data@parameters@data %>% as.data.frame() %>% mutate(file = x)
}) %>% bind_rows() %>% 
  separate(desc, c('clone', 'marker'), sep = '_') %>% 
  select(name, marker) %>% 
  dplyr::filter(!is.na(marker)) %>% 
  dplyr::filter(!(marker %in% c('separation', 'mahalanobis_dist', 'dist'))) %>% 
  unique() %>% 
  mutate(name = as.character(name))

df_data <- df_data %>% 
  # dplyr::select_('sample.id', 'treatment', 'time.point', 'label', 'batch', 'OID', df_parameters$name) %>% 
  dplyr::select(sample.id:OID, df_parameters$name) %>% 
  gather(name, value, -sample.id, -treatment, -time.point, -label, -batch, -OID) %>% 
  left_join(df_parameters) %>% 
  dplyr::filter(!(marker %in% c('DNA1', 'DNA2', 'length')))

df_data <- df_data %>% select(-name)

# Write data in csv format --------------------------------------------------------

write.csv2(df_data, file.path('..', 'Data', 'CellLines', 'df_data.csv'), row.names = FALSE)