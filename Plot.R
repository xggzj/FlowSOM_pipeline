# Making plot

# Rscript Plot.R

# 0. Load packages ----------

library(FlowSOM)
library(vite)
library(ggraph)
library(igraph)
library(tidyverse)
library(magrittr)

# 0. Load data ----------

load(file = './output/fSOM.RData')

# 1. Heatmap of FlowSOM clusters ----------

# Median marker expression of FlowSOM cell clusters
MFI = GetClusterMFIs(fSOM, colsUsed = T)
clusters = GetClusters(fSOM)
freqClusters = data.frame(clusters = clusters) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res = cbind(MFI, freqClusters)

# Heatmap
q99 = quantile(as.matrix(MFI), 0.99) 
grDevices::pdf(file = './output/heatmap.pdf', width = 16, height = 15)
heatmap = pheatmap::pheatmap(as.matrix(MFI),
                             scale = 'none',
                             labels_row = paste(rownames(MFI),' (', round(freqClusters$percentage,1), '%', ')',sep = ''),
                             display_numbers = TRUE,
                             angle_col = 45,
                             breaks = seq(0, q99, q99/90),
                             main = "Median marker expression per cluster")
print(heatmap)
dev.off()

# 2. Force-directed graph of FlowSOM cell clusters ----------

# Create unsupervised graph ----------
common_marker = colnames(res)[1:(ncol(res)-3)]
set.seed(824)
G <- vite::get_unsupervised_graph(res, 
                                  col.names = common_marker, 
                                  filtering.threshold = 5,
                                  method = 'forceatlas2',
                                  process.clusters.data = FALSE)
vite::write_graph(G, "./output/unsupervised.graphml")

# Make network plot ----------
l_g = create_layout(G, layout="fr")
l_g$x <- V(G)$x
l_g$y <- V(G)$y

ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage), fill = "#FEEDC3", shape = 21,  alpha = 0.9) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(aes(label = clusters), size = 3, repel = F, max.overlaps = 1000) +
  theme_graph(base_family = 'Helvetica')
ggsave(filename = './output/network_graph.pdf', width = 18, height = 16)

