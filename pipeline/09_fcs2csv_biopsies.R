# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'BiocManager')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

bioconductor.packages <- c('flowCore')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install("flowCore")

library(flowCore)
library(tidyverse)

tmp_files <- list.files(file.path(wd, 'PathSetter'), pattern = 'fcs')

df_data <- lapply(tmp_files, function(x){
  tmp_data <- read.FCS(file.path(wd, 'PathSetter', x), transformation = FALSE, emptyValue = FALSE)
  tmp_exp <- tmp_data@exprs %>% as.data.frame() %>% mutate(file = x)
}) %>% bind_rows() %>%
  mutate(file = gsub('-', '', file)) %>% 
  group_by(file) %>% 
  mutate(OID = 1:n()) %>% 
  ungroup() %>% 
  separate(file, c('sample.id','region', 'treatment' ,'label', 'batch'), sep = '_') %>% 
  mutate(batch = sub('.fcs', '', batch))

df_parameters <- lapply(tmp_files, function(x){
  tmp_data <- read.FCS(file.path(wd, 'PathSetter', x), transformation = FALSE, emptyValue = FALSE)
  tmp_exp <- tmp_data@parameters@data %>% as.data.frame() %>% mutate(file = x)
}) %>% bind_rows() %>% 
  separate(desc, c('clone', 'marker'), sep = '_') %>% 
  select(name, marker) %>% 
  dplyr::filter(!is.na(marker)) %>% 
  dplyr::filter(!(marker %in% c('viability', 'separation', 'mahalanobis_dist', 'dist'))) %>% 
  unique() %>% 
  mutate(name = as.character(name))

df_data <- df_data %>% 
  dplyr::select(sample.id:OID, df_parameters$name) %>% 
  gather(name, value, -sample.id, -treatment, -region, -label, -batch, -OID) %>% 
  left_join(df_parameters) %>% 
  dplyr::filter(!(marker %in% c('DNA1', 'DNA2', 'length')))

dir.create(file.path(wd, 'Data'))
write.csv2(df_data, file.path(wd, 'Data', 'df_data.csv'), row.names = FALSE)




