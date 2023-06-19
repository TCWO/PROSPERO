# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(caret)
library(e1071)
library(doParallel)

# Load data --------------------------------------------------------

df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'df_data_norm.csv'), stringsAsFactors = FALSE) %>% 
  filter(label != 'BC') %>% # remove batch control samples
  filter(batch == 1) %>% # remove batch 1
  filter(sample.id != 'Astrocytes') %>% # remove astrocytes
  mutate(time.point = ifelse(sample.id == 'BT360' & treatment == 'Control' & time.point == 48, 16, time.point)) %>% # some samples were incorrectly annotated
  mutate(time.point = ifelse(sample.id == 'LBT005' & treatment == 'Control' & time.point == 48, 16, time.point)) %>% 
  mutate(treatment = ifelse(sample.id == 'LBT062', plyr::mapvalues(treatment, from = c('Control', 'AMG232', 'RT'), to = c('Control', 'RT', 'AMG232')), treatment))

cell_cycle_markers <- c('IdU', 'CCNB1', 'pRB1', 'pH3', 'Ki67')

# Manual gating of cell cycle markers --------------------------------------------------------

df_data_sub <- df_data %>% 
  filter(marker %in% cell_cycle_markers) %>% 
  group_by(sample.id, treatment, time.point, label, batch, OID, marker) %>% 
  summarise(value = max(value)) %>% 
  ungroup() %>% 
  spread(marker, value) %>% 
  mutate(cell_cycle_phase = ifelse(pH3 > 6, 'M', 
                                   ifelse(IdU > 4.5, 'S',
                                          ifelse(CCNB1 > 6.25, 'G2', 'G0G1')))) %>% 
  select(-cell_cycle_markers)

# Scale complete data --------------------------------------------------------

write.csv2(df_data_sub, file.path('..', 'Data', 'CellLines', 'CellCycle', 'cell_cycle_annotations_batch1.csv'), row.names = FALSE)
