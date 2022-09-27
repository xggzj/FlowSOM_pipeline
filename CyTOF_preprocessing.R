# CyTOF data preprocessing ----------

# Rscript CyTOF_preprocessing.R -n_outlier -markers_ignore -ComBat -batch_filepath
# Example: Rscript CyTOF_preprocessing.R 1 CD3,gdT TRUE ./input/batch.csv


# Parameters:
# - n_outlier: an integer of threshold of outlier removal. for each marker, remove top n_outlier% cells with highest expression
# - markers_ignore: a character vector of markers to be ignored during outlier removal
# - ComBat: a logical value of doing ComBat batch correction or not
# - batch: path to batch table (a csv file of CyTOF exp_id. The first column: "file" - filename of CyTOF data; the second column: 'batch" - exp_id)

# 0. Load packages ----------

library(tidyverse)
library(magrittr)
library(flowCore)
'%nin%' <- Negate('%in%')

library(ggridges) 
library(umap)
library(sva)
library(ggpubr)

# 0. Load data ----------

args = commandArgs(trailingOnly=TRUE) 
n_outlier = as.integer(args[1][1]) 
markers_ignore = unlist(strsplit(args[2][1], split = ','))  
ComBat = args[3][1] 
batch = read.csv(args[4][1])
batch$file = as.character(batch$file)

# Read cytof data
filepaths_fcs = list.files(path = '.', pattern = '*.fcs', full.names = FALSE) 
files_fcs = lapply(filepaths_fcs, function(x){read.FCS(x, truncate_max_range = FALSE)})
names(files_fcs) = sapply(filepaths_fcs, function(x){strsplit(x, split = '\\.')[[1]][1]})

# Read grid annotation
filepaths_csv = list.files(path = '.', pattern = '*.csv', full.names = FALSE) 
files_csv = lapply(filepaths_csv, function(x){read.csv(x, row.names = 1)})
names(files_csv) = sapply(filepaths_csv, function(x){strsplit(x, split = '\\.')[[1]][1]})

# Check if files_fcs and files_csv are in the same order
identical(names(files_fcs), names(files_csv))

# 1. Channel filter ----------

# Remove non-marker channels
channelFilter = function(flowframe){
  exprs = exprs(flowframe)
  colnames(exprs)[-1] = markernames(flowframe)
  non_marker_col = c(which(colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")), grep("^[[:digit:]]+", colnames(exprs)))
  exprs = exprs[,-non_marker_col]
  return(exprs)
}
dat = lapply(files_fcs, channelFilter)

# Check marker names 
# (Note: some markers were named differently in different batches (gdTCR/TCRgd, CD3/CD3e))
#lapply(dat, colnames)
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
ncell_raw = Reduce('+', lapply(dat, nrow))

# Remove non-cells (level0 != cell) and dead cells (without level1)
for (i in 1:length(files_csv)) {
  level1_null = which(files_csv[[i]]$level1 == ' ')
  files_csv[[i]] = files_csv[[i]][-level1_null,]
  dat[[i]] = dat[[i]][-level1_null,]
}

# Check cell numbers after removal
ncell_remove_beads = Reduce('+', lapply(dat, nrow))

# 3. Remove outliers ----------

# Concatenate all files together
# .fcs
dat = lapply(dat, as.data.frame)
dat_all = data.table::rbindlist(dat)
remove(files_fcs)
# .csv
for (i in 1:length(files_csv)) {
  files_csv[[i]] %<>%
    mutate(file = names(files_csv)[[i]])
}
label_all = data.table::rbindlist(files_csv)
label_all %<>% select(-level0) %>%
  mutate(level2 = ifelse(level2 == ' ', yes = level1, no = level2))
remove(files_csv)

# Diagnosis: 
# Check distribution of each marker
#dat_all %>% summary()
marker_quantile = apply(dat_all, 2, function(x){quantile(x, c(0, 0.25, 0.5, 0.75, 0.99, 1))}) %>% t()
# Check cell counts of each population
ncell_per_pop_pre = label_all %>% group_by(level2) %>% tally()
colnames(ncell_per_pop_pre)[2] = 'pre'

# For each marker, remove top n_outlier% cells with highest expression 
marker_q99 = apply(dat_all, 2, function(x){quantile(x, 1-0.01*n_outlier)})
row_keep = which(dat_all[,1] <= marker_q99[1])
# Note: ignore CD123 and gdTCR to keep pDC and basophils
for (i in c(1:ncol(dat_all))[-c(1,which(colnames(dat_all) %in% markers_ignore))]) {
  row_keep = intersect(row_keep, which(dat_all[,..i] <= marker_q99[i]))
}
dat_all_keep_99 = dat_all[row_keep,]
label_all_keep_99 = label_all[row_keep,]

# Diagnosis after removal: 
# Check distribution of each marker
#dat_all_keep_99 %>% summary()
#marker_quantile_99 = apply(dat_all_keep_99, 2, function(x){quantile(x, c(0, 0.25, 0.5, 0.75, 0.99, 1))}) %>% t()
# Check cell counts in total after removal
ncell_remove_outlier = nrow(label_all_keep_99)
# Check cell counts of each population
ncell_per_pop_post = label_all_keep_99 %>% group_by(level2) %>% tally()
colnames(ncell_per_pop_post)[2] = 'post'

# 4. Arcsinh transformation ----------

dat_arcsinh = asinh(dat_all_keep_99/5)

# Check marker distribution in each cell pop
dat_df = cbind(dat_arcsinh, label_all_keep_99) %>% left_join(batch)
dat_df_gather = dat_df %>% select(c(all_of(common_marker), 'level2')) %>% gather(marker, expression, -level2)
# Density plot
gp_marker_density = ggplot(dat_df_gather, aes(x = expression, y = level2, fill = level2)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~marker, scales = 'free') +
  theme_bw() +
  labs(title = 'Marker distribution in all samples after outlier removal (arcsinh-transformed)')

# 5. Batch correction ----------

# 5.1 Check batch effects before correction
# Subsampling 1500 cells in each patient
dat_df_sub = dat_df %>% group_by(file) %>% sample_n(1500) %>% ungroup() 
# UMAP
umap.dat = umap(scale(dat_df_sub[,1:length(common_marker)]))
umap.df = data.frame(UMAP1 = umap.dat$layout[,1],
                     UMAP2 = umap.dat$layout[,2],
                     label = dat_df_sub$level2,
                     batch = dat_df_sub$batch)
# Plot
# Cell type
gp_umap_cell_before = ggplot(umap.df, aes(UMAP1, UMAP2, color = label)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cell type')) +
  labs(title = 'UMAP of raw data (subsampling 1500 cells per sample)')
# Batch
gp_umap_batch_before = ggplot(umap.df, aes(UMAP1, UMAP2, color = batch)) + 
  geom_point(size = 0.12, alpha = 0.08) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cytof batch')) +
  labs(title = 'UMAP of raw data (subsampling 1500 cells per sample)')

# 5.2 Batch correction
# ComBat
if(ComBat == 'TRUE'){
  dat_combat = ComBat(t(dat_arcsinh), batch = dat_df$batch) %>% t()
  dat_df_combat = cbind(dat_combat, dat_df %>% select(level1, level2, file, batch))
  # Check batch effects after correction
  # Subsampling 1500 cells in each patient
  dat_df_sub = dat_df_combat %>% group_by(file) %>% sample_n(1500) %>% ungroup() 
  # UMAP
  umap.dat = umap(scale(dat_df_sub[,1:length(common_marker)]))
  umap.df = data.frame(UMAP1 = umap.dat$layout[,1],
                       UMAP2 = umap.dat$layout[,2],
                       label = dat_df_sub$level2,
                       batch = dat_df_sub$batch)
  # Plot
  # Cell type
  gp_umap_cell_after = ggplot(umap.df, aes(UMAP1, UMAP2, color = label)) + 
    geom_point(size = 0.12, alpha = 0.08) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cell type')) +
    labs(title = 'UMAP after batch correction (subsampling 1500 cells per sample)')
  # Batch
  gp_umap_batch_after = ggplot(umap.df, aes(UMAP1, UMAP2, color = batch)) + 
    geom_point(size = 0.12, alpha = 0.08) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 5), title = 'cytof batch')) +
    labs(title = 'UMAP after batch correction (subsampling 1500 cells per sample)')
}else{
  dat_df_combat = dat_df
}

# 6. Z-score transformation ----------

#dat_scale = scale(dat_combat)

# 7. Save pre-processed data ----------

data = dat_combat
label = label_all_keep_99

# Save as RData
#save(data, label, file = './output/cytof_preprocessed.RData')

# Save as .fcs file
ff = flowFrame(data.matrix(data))
write.FCS(ff, filename = './output/cytof_preprocessed.fcs')
#write.csv(label, './output/cytof_preprocessed.csv', row.names = F)

# 8. Return a report ----------

common_marker_tb = data.frame(marker_name = paste(common_marker, collapse = ", ") %>% strwrap(width = 35) %>% paste(collapse = "\n"))
rownames(common_marker_tb) = 'Common marker'
table1 = ggtexttable(common_marker_tb, theme = ttheme("light")) %>% 
  tab_add_title(text = 'Markers used for clustering', face = "bold", padding = unit(3, "line"))

ncell = data.frame(ncell = c(ncell_raw, ncell_remove_beads, ncell_remove_outlier))
rownames(ncell) = c('Raw data', 'After removing non-cells and dead cells', 'After outlier removal')
table2 = ggtexttable(ncell, theme = ttheme("light"))  %>% 
  tab_add_title(text = 'Total cell number', face = "bold", padding = unit(3, "line"))

ncell_per_pop = ncell_per_pop_pre %>% left_join(ncell_per_pop_post) %>% column_to_rownames(var = 'level2')
table3 = ggtexttable(ncell_per_pop, theme = ttheme(base_style = "light")) %>% 
  tab_add_title(text = "Cell counts before and after outlier removal" %>% strwrap(width = 35) %>% paste(collapse = "\n"), face = "bold", padding = unit(3, "line"))

marker_quantile_1 = marker_quantile[1:round(nrow(marker_quantile)/2),] %>% round(digits = 2)
marker_quantile_2 = marker_quantile[(round(nrow(marker_quantile)/2)+1):nrow(marker_quantile),] %>% round(digits = 2)
table4 = ggtexttable(marker_quantile_1, theme = ttheme("light")) %>% 
  tab_add_title(text = 'Marker distribution in all samples before outlier removal', face = "bold", size = 10, padding = unit(3, "line")) 
table5 = ggtexttable(marker_quantile_2, theme = ttheme("light"))

grDevices::pdf('./output/cytof_preprocessing_report.pdf', width = 20, height = 15)
print(ggarrange(table1, table2, table3, ncol = 3))
print(ggarrange(table4, table5))
print(gp_marker_density)
if(ComBat == 'TRUE'){
  print(ggarrange(gp_umap_cell_before, gp_umap_batch_before, gp_umap_cell_after, gp_umap_batch_after))
}else{
  print(ggarrange(gp_umap_cell_before, gp_umap_batch_before))
}
dev.off()
