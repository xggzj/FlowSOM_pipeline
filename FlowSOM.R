# FlowSOM clustering

# 0. Load packages ----------

library(FlowSOM)
library(flowCore)
library(tidyverse)
library(magrittr)

# 0. Load data ----------

ff = read.FCS('cytof_94pat_32marker_outlier99_arcsinh5_combat_data.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)
label = read.csv('cytof_94pat_32marker_outlier99_arcsinh5_combat_label.csv', 
                 check.names = F, stringsAsFactors = F) 
sample_info = read.csv('baseline_sample_info.csv', stringsAsFactors = FALSE) 

# 1. FlowSOM 30 clusters ----------

common_marker = colnames(ff@exprs) %>% unname()

# Run FlowSOM
fSOM_30 = FlowSOM(ff,
                  compensate = F,
                  transform = F,
                  scale = F,
                  colsToUse = common_marker,
                  nClus = 10,
                  xdim = 5, ydim = 6,
                  seed = 824)

# FlowSOM results summary
FlowSOMmary(fSOM_30, plotFile = 'flowsom_summary_30clusters.pdf')

# Median marker expression of FlowSOM cell clusters
MFI_30 = GetClusterMFIs(fSOM_30, colsUsed = T)
clusters_30 = GetClusters(fSOM_30)
freqClusters_30 = data.frame(clusters = clusters_30) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_30 = cbind(MFI_30, freqClusters_30)

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_30), 0.99)
pheatmap::pheatmap(as.matrix(MFI_30),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_30),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "30 clusters: Median marker expression per cluster")

# 2. Separate Neutrophil clusters ----------

# Identify Neutrophil clusters according to heatmap
nrow_neutrophils = which(clusters_30 %in% c(12,26,28,29,18,24,27,17,23,25))

# Separate Neutrophil data
label_neutrophils = label[nrow_neutrophils,]
clusters_neutrophils = clusters_30[nrow_neutrophils]
label_neutrophils$lineage = 'Neutrophils'
label_neutrophils$subtype = 'Neutrophils'
label_neutrophils$clusters = clusters_neutrophils
data_neutrophils = ff@exprs[nrow_neutrophils,]
res_neutrophils = res_30 %>% dplyr::filter(clusters %in% c(12,26,28,29,18,24,27,17,23,25))
res_neutrophils$lineage = 'Neutrophils'
res_neutrophils$subtype = 'Neutrophils'

# Remaining cells after removing Neutrophils
data_remaining = ff@exprs[-nrow_neutrophils,]
label_remaining = label[-nrow_neutrophils,]

# 3. FlowSOM 100 clusters on remaining cells ----------

ff_remaining = ff
ff_remaining@exprs = data_remaining

# Run FlowSOM
fSOM_100 = FlowSOM(ff_remaining,
                   compensate = F,
                   transform = F,
                   scale = F,
                   colsToUse = common_marker,
                   nClus = 10,
                   xdim = 10, ydim = 10,
                   seed = 824)

# FlowSOM results summary
FlowSOMmary(fSOM_100, plotFile = 'flowsom_summary_100clusters.pdf')

# Median marker expression of FlowSOM cell clusters
MFI_100 = GetClusterMFIs(fSOM_100, colsUsed = T)
clusters_100 = GetClusters(fSOM_100)
freqClusters_100 = data.frame(clusters = clusters_100) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_100 = cbind(MFI_100, freqClusters_100)

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_100), 0.99)
pheatmap::pheatmap(as.matrix(MFI_100),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_100),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "100 clusters: Median marker expression per cluster")

# Manual annotation
annotation = read.csv2('ISAC_cluster_annotation.csv', check.names = F, stringsAsFactors = F)
res_100 %<>% left_join(annotation)
label_remaining$clusters = clusters_100
label_remaining %<>% left_join(annotation)

# 4. Bind remaining cells with neutrophils ----------

data_bind = rbind(data_neutrophils, data_remaining)
label_bind = rbind(label_neutrophils, label_remaining)
label_bind %<>% mutate(cell_cluster = paste(subtype, '_', clusters, sep = '')) %>% left_join(sample_info)
all_bind = cbind(data_bind, label_bind) 
res_bind = rbind(res_neutrophils, res_100)
res_bind %<>% mutate(frequency = (n/sum(res_bind$n))*100)
res_bind %<>% mutate(cell_cluster = paste(subtype, '_', clusters, sep = ''))

# Heatmap of all cell clusters
q99 = quantile(as.matrix(res_bind[,1:length(common_marker)]), 0.99)
pheatmap::pheatmap(as.matrix(res_bind[,1:length(common_marker)]),
                   scale = 'none',
                   labels_row = paste(res_bind$cell_cluster,' (', round(res_bind$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster")

# Calculate relative frequency of each celltype in each sample
freq_lineage = label_bind %>%
  group_by(study_id) %>%
  count(.data$lineage) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  spread(lineage, percentage)
freq_subtype = label_bind %>%
  group_by(study_id) %>%
  count(.data$subtype) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  spread(subtype, percentage)
freq_cellcluster = label_bind %>%
  group_by(study_id) %>%
  count(.data$cell_cluster) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  spread(cell_cluster, percentage)

# 5. Save results ----------

save(all_bind, file = 'flowsom_results.RData')
write.csv(freq_lineage, file = 'freq_lineage.csv', row.names = F)
write.csv(freq_subtype, file = 'freq_subtype.csv', row.names = F)
write.csv(freq_cellcluster, file = 'freq_cellcluster.csv', row.names = F)

# Save median marker expression data for Network plot
write.table(res_bind, file = 'flowsom_clustered.txt', 
            sep="\t", row.names=F, col.names = T, quote = F)


