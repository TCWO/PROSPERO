# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)

# Load data --------------------------------------------------------

df_dictionary <- read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'Clusters_annotated', 'flowsom_clusters_sorted_batch1_annotated.csv'), stringsAsFactors = FALSE) %>% 
  bind_rows(read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'Clusters_annotated', 'kmeans_clusters_sorted_batch1_annotated.csv'), stringsAsFactors = FALSE)) %>% 
  bind_rows(read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'Clusters_annotated', 'phenograph_clusters_sorted_batch1_annotated.csv'), stringsAsFactors = FALSE)) %>% 
  bind_rows(read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'Clusters_annotated', 'clara_clusters_sorted_batch1_annotated.csv'), stringsAsFactors = FALSE)) %>% 
  rename(Phenotype = Annotation) %>% 
  select(clustering_method, cluster, Phenotype) %>%
  mutate(Phenotype = ifelse(is.na(Phenotype), 'other', Phenotype)) %>% 
  mutate(Phenotype = ifelse(Phenotype == 'mixed', 'mixed.high', Phenotype),
         Phenotype = ifelse(Phenotype == 'blank', 'mixed.low', Phenotype))

df_data <- read.csv2(file.path('..', 'Data', 'CellLines', 'Clustering', 'df_data_clusters_batch1.csv'), stringsAsFactors = FALSE) %>% 
  mutate(ID = paste(sample.id ,treatment, time.point, label, batch, sep = '_'))

phenotypic_markers <- c('GFAP', 'Olig2', 'Nestin', 'Vimentin', 'PDGFRa', 'CD44', 'CD24')

# Add annotations to clusters --------------------------------------------------------

tmp_data <- df_data %>% 
  select(ID, OID, phenograph, kmeans, flowsom, clara) %>% 
  gather(clustering_method, cluster, -ID, -OID) %>% 
  left_join(df_dictionary)

table(tmp_data[c('Phenotype', 'clustering_method')])

# Naming consensus --------------------------------------------------------

tmp_clusters <- tmp_data %>% select(ID, OID, clustering_method, Phenotype) %>% spread(clustering_method, Phenotype) %>% 
  group_by(phenograph, flowsom, kmeans, clara) %>% summarise(N = n()) %>% ungroup()

df_agreement4 <- data.frame()
df_agreement3 <- data.frame()
df_agreement2 <- data.frame()
df_noise <- data.frame()

tmp_counts <- 0

for(i in 1:nrow(tmp_clusters)){
  tmp_group <- tmp_clusters[i,]
  
  tmp_name <- tmp_group %>% 
    select(-N) %>% 
    gather(CM, Name)
  
  new_names <- tmp_name
  
  new_names <- new_names %>% 
    select(Name) %>% 
    table() %>% 
    as.data.frame()
  colnames(new_names) <- c('Name', 'Freq')
  
  if(max(new_names$Freq) == 4) {
    tmp_max <- new_names %>% dplyr::filter(Freq == 4)
    if(nrow(tmp_max) > 1){
      asdsad
      tmp_cell_type <- (tmp_max %>% sample_n(1))$Name
    } else{
      tmp_cell_type <- (new_names %>% dplyr::filter(Freq == 4))$Name
    }
    df_agreement4 <- bind_rows(df_agreement4, tmp_group %>% mutate(cell_type = tmp_cell_type))
  }else if(max(new_names$Freq) %in% c(2,3)){
    tmp_max <- new_names %>% dplyr::filter(Freq == max(Freq))
    if(nrow(tmp_max) > 1){
      tmp_counts <- tmp_counts + 1
      tmp_cell_type <- (tmp_name %>% dplyr::filter(CM == 'clara'))$Name
      # tmp_cell_type <- paste(tmp_max$Name, collapse = '_')
    } else{
      tmp_cell_type <- (new_names %>% dplyr::filter(Freq == max(Freq)))$Name
    }
    if(max(new_names$Freq) == 3){
      df_agreement3 <- bind_rows(df_agreement3, tmp_group %>% mutate(cell_type = tmp_cell_type))
    } else {
      df_agreement2 <- bind_rows(df_agreement2, tmp_group %>% mutate(cell_type = tmp_cell_type))
    }
  }else{
    df_noise <- bind_rows(df_noise, tmp_group %>% mutate(cell_type = 'NOISE'))
  }
}

'Agreement of 4'
sum(df_agreement4$N)
'Agreement of 3'
sum(df_agreement3$N)
'Agreement of 2'
sum(df_agreement2$N)
'Inconsistent'
sum(df_noise$N)

# Substitute cluster annotations for final cell types --------------------------------------------------------

df_data <- df_data %>% 
  mutate(phenograph = plyr::mapvalues(x = phenograph, from = (df_dictionary %>% dplyr::filter(clustering_method == 'phenograph'))$cluster,
                                      to = (df_dictionary %>% dplyr::filter(clustering_method == 'phenograph'))$Phenotype),
         flowsom = plyr::mapvalues(x = flowsom, from = (df_dictionary %>% dplyr::filter(clustering_method == 'flowsom'))$cluster,
                                   to = (df_dictionary %>% dplyr::filter(clustering_method == 'flowsom'))$Phenotype),
         kmeans = plyr::mapvalues(x = kmeans, from = (df_dictionary %>% dplyr::filter(clustering_method == 'kmeans'))$cluster,
                                  to = (df_dictionary %>% dplyr::filter(clustering_method == 'kmeans'))$Phenotype),
         clara = plyr::mapvalues(x = clara, from = (df_dictionary %>% dplyr::filter(clustering_method == 'clara'))$cluster,
                                 to = (df_dictionary %>% dplyr::filter(clustering_method == 'clara'))$Phenotype))

df_data <- df_data %>% 
  mutate(phenograph = ifelse(phenograph %in% unique(df_dictionary$Phenotype), phenograph, 'other'),
         flowsom = ifelse(flowsom %in% unique(df_dictionary$Phenotype), flowsom, 'other'),
         kmeans = ifelse(kmeans %in% unique(df_dictionary$Phenotype), kmeans, 'other'),
         clara = ifelse(clara %in% unique(df_dictionary$Phenotype), clara, 'other'))

df_data <- df_data %>% select(-N, -M, -P)

df_data_r <- df_data %>% 
  left_join(df_agreement4 %>% bind_rows(df_agreement3) %>% bind_rows(df_agreement2)) %>% 
  dplyr::filter(!is.na(cell_type)) %>% 
  select(-phenograph, -flowsom, -kmeans, -clara, -N)

df_data_r <- df_data_r %>% 
  mutate(cell_type = ifelse(cell_type == 'blank', 'mixed_low', cell_type),
         cell_type = ifelse(cell_type == 'mixed', 'mixed_high', cell_type))

# Write results --------------------------------------------------------

write.csv2(df_data_r, file.path('..', 'Data', 'CellLines', 'Clustering', 'df_data_consensus_batch1.csv'), row.names = FALSE)
