# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'caret', 'e1071', 'doParallel')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(caret)
library(e1071)
library(doParallel)

# Load data --------------------------------------------------------

df_data_all <- read.csv2(file.path('..', 'Data', 'CellLines', 'df_data_norm.csv'), stringsAsFactors = FALSE) %>% 
  filter(label != 'BC') %>% # remove batch control samples
  filter(batch == 1) %>% # filter batch 1
  filter(sample.id != 'Astrocytes') %>% # remove astrocytes
  mutate(time.point = ifelse(sample.id == 'BT360' & treatment == 'Control' & time.point == 48, 16, time.point)) %>% # some samples were incorrectly annotated
  mutate(time.point = ifelse(sample.id == 'LBT005' & treatment == 'Control' & time.point == 48, 16, time.point)) %>% 
  mutate(treatment = ifelse(sample.id == 'LBT062', plyr::mapvalues(treatment, from = c('Control', 'AMG232', 'RT'), to = c('Control', 'RT', 'AMG232')), treatment))

df_data_training <- read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'df_data_consensus_batch1.csv'), stringsAsFactors = FALSE)

phenotypic_markers <- c('GFAP', 'Olig2', 'Nestin', 'Vimentin', 'PDGFRa', 'CD44', 'CD24')

# Train SVM model --------------------------------------------------------

library(doParallel)
n_cores <- 20
n_cores_max <- detectCores()
cl <- makePSOCKcluster(min(n_cores, (n_cores_max-1)))
registerDoParallel(cl)

# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

tmp_svm <- train(cell_type ~., data = df_data_training %>% select(cell_type, phenotypic_markers), method = "svmRadial", trControl = train_control, tuneLength = 20)
saveRDS(tmp_svm, file.path('..', 'Data', 'CellLines', 'svm_models', 'svm_model_tumor_subtypes_batch1.RDS'))

stopCluster(cl)

# Scale complete data --------------------------------------------------------

df_data_all <- df_data_all %>% 
  group_by(sample.id, treatment, time.point, label, batch, OID, marker) %>%
  summarise(value = max(value)) %>% 
  ungroup() %>% 
  group_by(marker) %>% 
  dplyr::mutate(value = scale(value)) %>% 
  ungroup() %>% 
  mutate(value = ifelse(value > 3, 3, value)) %>% 
  mutate(value = ifelse(value < -3, -3, value)) %>% 
  spread(marker, value)

df_data_all <- df_data_all %>% 
  mutate(cell_type = predict(tmp_svm, .) %>% as.character())

# Scale complete data --------------------------------------------------------

write.csv2(df_data_all, file.path('..', 'Data', 'CellLines', 'Clustering', 'data_annotated_tumor_subtypes_batch1.csv'), row.names = FALSE)
