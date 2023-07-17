# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse", 'umap', 'Rtsne', 'devtools', 'BiocManager', 'fpc', 'cluster')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

bioconductor.packages <- c('FlowSOM')
new.packages <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

github.packages <- c('Rphenograph')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) devtools::install_github("JinmiaoChenLab/Rphenograph")

library(tidyverse)
library(umap)
library(Rphenograph)
library(FlowSOM)
library(fpc)
library(cluster)

# Load data --------------------------------------------------------

tmp_files <- list.files(file.path(wd, 'SOX2P'), pattern = 'csv')

df_data <- lapply(tmp_files, function(x){
  tmp_data <- read.csv(file.path(wd, 'SOX2P', x))
  tmp_exp <- tmp_data %>% as.data.frame() %>% rename(SOX2P = selected_) %>% select(-c(Sox2h,CD45h))
}) %>% bind_rows() 

write.csv2(df_data, file.path(wd, 'Data', 'df_data3.csv'), row.names = FALSE)


#Load experimental design -------------------------------------------------
meta_table <- read.csv2(file.path(wd,'Data','Metatable_Biopsies_forPouya_FDS.csv')) %>% rename(sample.id=Sample_ID) %>%
  rename(region = Region) %>%
  rename(treatment = Treatment) %>%
  rename(batch = Batch) %>%
  rename(sample = Sample.BC) %>% select(-X,-X.1,-Original.filename)


df_data <- df_data %>% left_join(meta_table) %>% 
  gather(marker,value, BAX:Vimentin) %>% 
  dplyr::filter(label!='BC')


df_data <- df_data %>%
  group_by(sample.id, region, treatment, label, batch, OID, marker) %>% 
  summarise(value = max(value)) %>% 
  ungroup()


df_data <- df_data %>%
  group_by(marker) %>%
  mutate(value = scale(value)) %>% 
  ungroup()


df_data <- df_data %>% 
  mutate(value = ifelse(value > 5, 5, value),
         value = ifelse(value < -5, -5, value))

# Defining the markers used to cluster Sox2 positive cells ----------------------------------
markers_clust <- c('CD44','CD24', 'Vimentin','Nestin', 'PDGFRa', 'GFAP','Olig2')


# Sample data following a stratified random sampling scheme --------------------------------------------------------
df_data <- df_data %>% select(c('sample.id','region','treatment','label','batch','OID'),markers_clust)
set.seed(1234)

tmp_reduced <- df_data %>% 
  group_by(sample.id, region, treatment) %>% 
  mutate(N = n()) %>% 
  ungroup()

tmp_cells <- tmp_reduced %>% 
  select(sample.id, region, treatment, label, batch, N) %>% 
  unique() %>% 
  mutate(M = sum(N)) %>% 
  mutate(P = round(25000*N/M))

tmp_reduced <- tmp_reduced %>% 
  left_join(tmp_cells) %>% 
  group_by(sample.id, region, treatment) %>% 
  sample_n(P) %>% ungroup() %>% distinct(CD44,CD24, Olig2,Vimentin,Nestin,GFAP,PDGFRa , .keep_all = TRUE)

tmp_data <- tmp_reduced %>%
  select(markers_clust) %>% 
  as.matrix() %>% unique() 


custom.config <- umap.defaults

custom.config$metric <- 'cosine'
custom.config$n_neighbors <- 10

tmp_umap <- umap(tmp_data) #, custom.config)
tmp_tsne <- Rtsne::Rtsne(tmp_data, perplexity = 20)
tmp_pca <- prcomp(tmp_data)


tmp_phenograph <- Rphenograph(tmp_data)
tmp_clusters_phenograph <- as.character(membership(tmp_phenograph[[2]]))



tmp_reduced <- tmp_reduced %>%
  mutate(uMap1 = tmp_umap$layout[,1], uMap2 = tmp_umap$layout[,2]) %>%
  mutate(tsne1 = tmp_tsne$Y[,1], tsne2 = tmp_tsne$Y[,2]) %>%
  mutate(PC1 = tmp_pca$x[,1], PC2 = tmp_pca$x[,2]) %>% 
  mutate(phenograph = tmp_clusters_phenograph)

# Write data --------------------------------------------------------
write.csv2(x = tmp_reduced, file = file.path(wd, 'Data', 'df_data_clusters_Sox2P.csv'), row.names = FALSE)















