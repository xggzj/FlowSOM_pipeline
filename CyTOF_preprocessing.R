# CyTOF data preprocessing

# 0. Load packages ----------

library(tidyverse)
library(magrittr)
library(flowCore)
'%nin%' <- Negate('%in%')

library(ggridges) 
library(umap)
library(sva)

# 0. Load data ----------

# Baseline sample info
sample_info = read.csv('baseline_sample_info.csv', stringsAsFactors = FALSE) 
sample_info$grid_id = as.character(sample_info$grid_id)
sample_info$olink_id = as.character(sample_info$olink_id)
str(sample_info)

# Check cytof file exists
all(file.exists(sample_info$path_cytof_fcs))
all(file.exists(sample_info$path_cytof_csv))

# Read cytof data
files_fcs = lapply(sample_info$path_cytof_fcs, function(x){read.FCS(x, truncate_max_range = FALSE)})
names(files_fcs) = sample_info$path_cytof_fcs

# Read grid annotation
files_csv = lapply(sample_info$path_cytof_csv, function(x){read.csv(x, row.names = 1)})
names(files_csv) = sample_info$path_cytof_csv

# 1. Channel filter ----------

# Remove empty channels
channelFilter = function(flowframe){
  exprs = exprs(flowframe)[,-1]
  colnames(exprs) = markernames(flowframe)
  nonEmptyChannels = pData(parameters(flowframe))[c("name", "desc")] %>% dplyr::filter(name != desc)
  exprs = exprs[,colnames(exprs) %in% nonEmptyChannels$desc]
  return(exprs)
}
dat = lapply(files_fcs, channelFilter)

# Check marker names 
# (Note: some markers were named differently in different batches (gdTCR/TCRgd, CD3/CD3e))
lapply(dat, colnames)
# Replace TCRgd with gdTCR, replace CD3e with CD3
dat = lapply(dat, function(x){
  colnames(x)[which(colnames(x) == 'TCRgd')] = 'gdTCR'
  colnames(x)[which(colnames(x) == 'CD3e')] = 'CD3'
  return(x)
})

# Common markers
common_marker = colnames(dat[[1]])
for (i in 1:(length(dat)-1)) {
  common_marker = intersect(common_marker, colnames(dat[[i+1]]))
} 
# (check if 'EQBeads','DNA-Ir191','DNA-Ir193' were excluded)
# (exclude markers that are not needed)
common_marker

# Keep only common markers
dat = lapply(dat, function(x){x = x[,common_marker]})

# 2. Remove non-cells and dead cells (according to Grid annotation) ----------

# Check cell numbers before removal
Reduce('+', lapply(dat, nrow))

# Remove non-cells (level0 != cell) and dead cells (without level1)
for (i in 1:length(files_csv)) {
  level1_null = which(files_csv[[i]]$level1 == ' ')
  files_csv[[i]] = files_csv[[i]][-level1_null,]
  dat[[i]] = dat[[i]][-level1_null,]
}

# Check cell numbers after removal
Reduce('+', lapply(dat, nrow))

# 3. Remove outliers ----------

# Concatenate all files together
# .fcs
dat = lapply(dat, as.data.frame)
dat_all = data.table::rbindlist(dat)
remove(files_fcs)
# .csv
for (i in 1:length(files_csv)) {
  files_csv[[i]] %<>%
    mutate(path_cytof_csv = names(files_csv)[[i]])
}
label_all = data.table::rbindlist(files_csv)
label_all %<>% select(-level0) %>%
  mutate(level2 = ifelse(level2 == ' ', yes = level1, no = level2))
remove(files_csv)

# Diagnosis: 
# Check distribution of each marker
dat_all %>% summary()
marker_quantile = apply(dat_all, 2, function(x){quantile(x, c(0, 0.25, 0.5, 0.75, 0.99, 1))}) %>% t()
# Check cell counts of each population
label_all %>% group_by(level2) %>% tally()

# For each marker, remove top 1% cells with highest expression 
marker_q99 = apply(dat_all, 2, function(x){quantile(x, 0.99)})
row_keep = which(dat_all[,1] <= marker_q99[1])
# Note: ignore CD123 and gdTCR to keep pDC and basophils
which(colnames(dat_all) %in% c('CD123', 'gdTCR'))
for (i in c(1:ncol(dat_all))[-c(1,which(colnames(dat_all) %in% c('CD123', 'gdTCR')))]) {
  row_keep = intersect(row_keep, which(dat_all[,..i] <= marker_q99[i]))
}
dat_all_keep_99 = dat_all[row_keep,]
label_all_keep_99 = label_all[row_keep,]

# Diagnosis after removal: 
# Check distribution of each marker
dat_all_keep_99 %>% summary()
marker_quantile_99 = apply(dat_all_keep_99, 2, function(x){quantile(x, c(0, 0.25, 0.5, 0.75, 0.99, 1))}) %>% t()
# Check cell counts in total after removal
nrow(label_all_keep_99)
# Check cell counts of each population
label_all_keep_99 %>% group_by(level2) %>% tally()

# 4. Arcsinh transformation ----------

dat_arcsinh = asinh(dat_all_keep_99/5)

# Check marker distribution in each cell pop
dat_df = cbind(dat_arcsinh, label_all_keep_99) %>% left_join(sample_info)
dat_df_gather = dat_df[,c(1:31,33)] %>% gather(marker, expression, -level2)
# Density plot
ggplot(dat_df_gather, aes(x = expression, y = level2, fill = level2)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~marker, scales = 'free') +
  theme_bw()

# 5. Batch correction ----------

# 5.1 Check batch effects before correction
# Subsampling 1500 cells in each patient
dat_df_sub = dat_df %>% group_by(study_id) %>% sample_n(1500) %>% ungroup() 
# UMAP
umap.dat = umap(scale(dat_df_sub[,1:length(common_marker)]))
umap.df = data.frame(UMAP1 = umap.dat$layout[,1],
                     UMAP2 = umap.dat$layout[,2],
                     label = dat_df_sub$level2,
                     batch = dat_df_sub$cytof_batch)
# Plot
# Cell type
ggplot(umap.df, aes(UMAP1, UMAP2, color = label)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  scale_color_manual(values = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                                '#9467bd', '#8c564b', '#e377c2', 'yellow', 
                                '#bcbd22', '#17becf', 'blue')) + 
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cell type'))
# Batch
ggplot(umap.df, aes(UMAP1, UMAP2, color = batch)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cytof batch'))

# 5.2 Batch correction
# ComBat
dat_combat = ComBat(t(dat_df[,1:length(common_marker)]), batch = dat_df$cytof_batch) %>% t()
dat_df_combat = cbind(dat_combat, label_all_keep_99) %>% left_join(sample_info)

# 5.3 Check batch effects after correction
# Subsampling 1500 cells in each patient
dat_df_sub = dat_df_combat %>% group_by(study_id) %>% sample_n(1500) %>% ungroup() 
# UMAP
umap.dat = umap(scale(dat_df_sub[,1:length(common_marker)]))
umap.df = data.frame(UMAP1 = umap.dat$layout[,1],
                     UMAP2 = umap.dat$layout[,2],
                     label = dat_df_sub$level2,
                     batch = dat_df_sub$cytof_batch)
# Plot
# Cell type
ggplot(umap.df, aes(UMAP1, UMAP2, color = label)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  scale_color_manual(values = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                                '#9467bd', '#8c564b', '#e377c2', 'yellow', 
                                '#bcbd22', '#17becf', 'blue')) + 
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cell type'))
# Batch
ggplot(umap.df, aes(UMAP1, UMAP2, color = batch)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cytof batch'))

# 6. Z-score transformation ----------

#dat_scale = scale(dat_combat)

# 7. Save pre-processed data ----------

data = dat_combat
label = label_all_keep_99

# Save as RData
save(data, label, file = 'cytof_94pat_32marker_outlier99_arcsinh5_combat.RData')

# Save as .fcs file
ff = flowFrame(data.matrix(data))
write.FCS(ff, filename = 'cytof_94pat_32marker_outlier99_arcsinh5_combat_data.fcs')
write.csv(label, 'cytof_94pat_32marker_outlier99_arcsinh5_combat_label.csv', row.names = F)
