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

df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'df_data_norm.csv'), stringsAsFactors = FALSE) 

# Filter batch 1, remove astrocytes, and trim negative values to 0 --------------------------------------------------------

df_data <- df_data %>% 
  dplyr::filter(batch != 1) %>% 
  dplyr::filter(sample.id != 'Astrocytes') %>% 
  mutate(value = ifelse(value < 0, 0, value))

# Define phenotypic markers and filter in data --------------------------------------------------------

phenotypic_markers <- c('GFAP', 'Olig2', 'Nestin', 'Vimentin', 'PDGFRa', 'CD44', 'CD24')

df_data <- df_data %>% 
  filter(marker %in% phenotypic_markers)

# Scale data (z-scores) and trim in the [-3, +3] range --------------------------------------------------------

df_data <- df_data %>% 
  group_by(marker) %>% 
  mutate(value = scale(value)) %>% 
  ungroup() %>% 
  mutate(value = ifelse(value > 3, 3, value),
         value = ifelse(value < -3, -3, value)) %>% 
  spread(marker, value)

# Sample data following a stratified random sampling scheme --------------------------------------------------------

# count number of cells per sample
tmp_reduced <- df_data %>% 
  group_by(sample.id, treatment, time.point, label, batch) %>% 
  mutate(N = n()) %>% 
  ungroup()

# define number of cells to be sampled
n_cells <- 10000

tmp_cells <- tmp_reduced %>% 
  select(sample.id, treatment, time.point, label, batch, N) %>% 
  unique() %>% 
  mutate(M = sum(N)) %>% 
  mutate(P = round(n_cells*N/M))

# sample cells with seed for reproducibility
tmp_seed <- 1234
set.seed(tmp_seed)

tmp_reduced <- tmp_reduced %>% 
  left_join(tmp_cells) %>% 
  group_by(sample.id, treatment, time.point, label, batch) %>% 
  sample_n(unique(P)) %>% 
  ungroup()

# Dimensionality reduction --------------------------------------------------------

tmp_matrix <- tmp_reduced %>%
  select(phenotypic_markers) %>% 
  as.matrix() %>% 
  unique()

tmp_umap <- umap(tmp_matrix)
set.seed(tmp_seed)
tmp_tsne <- Rtsne::Rtsne(tmp_matrix, perplexity = 20)
tmp_pca <- prcomp(tmp_matrix)

# Clustering --------------------------------------------------------

tmp_phenograph <- Rphenograph(tmp_matrix, k = 10)
tmp_clusters_phenograph <- as.character(membership(tmp_phenograph[[2]]))

tmp_fsom <- FlowSOM(tmp_matrix, colsToUse = c(1:ncol(tmp_matrix)), nClus = n_distinct(tmp_clusters_phenograph)) #, maxMeta = n_distinct(membership(tmp_clusters_phenograph[[2]])))
tmp_clusters_fsom <- tmp_fsom$metaclustering[GetClusters(tmp_fsom)]

set.seed(tmp_seed)
tmp_cluster_kmeans <- kmeans(tmp_matrix,  n_distinct(tmp_clusters_phenograph), iter.max = 100)$cluster

tmp_cluster_pam <- clara(tmp_matrix, n_distinct(tmp_clusters_phenograph), metric = 'manhattan', stand=FALSE, samples=10000, pamLike=TRUE)
tmp_cluster_clara <- tmp_cluster_pam$clustering

tmp_matrix <- tmp_matrix %>%
  as.data.frame() %>% 
  mutate(uMap1 = tmp_umap$layout[,1], uMap2 = tmp_umap$layout[,2]) %>%
  mutate(tsne1 = tmp_tsne$Y[,1], tsne2 = tmp_tsne$Y[,2]) %>%
  mutate(PC1 = tmp_pca$x[,1], PC2 = tmp_pca$x[,2]) %>% 
  mutate(phenograph = tmp_clusters_phenograph,
         flowsom = tmp_clusters_fsom,
         kmeans = tmp_cluster_kmeans,
         clara = tmp_cluster_clara)

tmp_reduced <- tmp_reduced %>% left_join(tmp_matrix)

# Write data --------------------------------------------------------

write.csv2(x = tmp_reduced, file = file.path('..', 'Data', 'CellLines', 'Clustering', 'df_data_clusters_batches2and3.csv'), row.names = FALSE)

# Generate fingerprints for each cluster --------------------------------------------------------

df_heatmap <- tmp_reduced %>%
  select(phenograph, flowsom, kmeans,clara, phenotypic_markers) %>%
  gather(marker, value, -phenograph, -flowsom, -kmeans, -clara) %>%
  gather(clustering_method, cluster, -marker, -value) %>% 
  group_by(clustering_method, cluster, marker) %>%
  summarise(average_expression = mean(value)) %>%
  ungroup()

df_heatmap <- df_heatmap %>% mutate(cluster = factor(cluster, levels = c(1:n_distinct(df_heatmap$cluster))))

hmcol <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(9))
df_sorted <- data.frame()

for(i in unique(df_heatmap$clustering_method)){
  tmp_heatmap <- df_heatmap %>%
    filter(clustering_method == i) %>%
    droplevels()
  
  tmp_plot <- tmp_heatmap %>%
    group_by(cluster) %>%
    arrange(-average_expression) %>%
    mutate(N = 1:n()) %>%
    ungroup() %>%
    filter(N <= 15)
  
  tmp_plot <- tmp_plot %>%
    filter(N %in% factor(c(1:max(tmp_plot$N), 'Percentages'))) %>%
    dplyr::select(-average_expression) %>%
    spread(N, marker) %>%
    mutate(Rank = 1:n()) %>%
    arrange(-Rank) %>%
    dplyr::select(-Rank)
  write.csv2(x = tmp_plot, file = file.path('..', 'Data', 'CellLines', 'Clustering', 'Clusters_per_method', paste0(i, '_clusters_sorted_batches2and3.csv')), row.names = FALSE)
}
