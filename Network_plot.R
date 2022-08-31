# Force-directed graph of FlowSOM cell clusters

# 0. Load packages ----------

library(vite)
library(ggraph)
library(igraph)
library(tidyverse)
library(magrittr)

# 0. Load data ----------

res_df = read.csv('treekoR_res_Neuroblastoma - Metastasis vs. Non-metastasis.csv')
res_bind = read.table('flowsom_clustered.txt', header = TRUE, sep = "\t", check.names = FALSE)

# 1. Create unsupervised graph ----------

common_marker = colnames(res_bind)[1:32]

set.seed(824)
G <- vite::get_unsupervised_graph(res_bind, 
                                  col.names = common_marker, 
                                  filtering.threshold = 5,
                                  method = 'forceatlas2',
                                  process.clusters.data = FALSE)
vite::write_graph(G, "unsupervised.graphml")

# 2. Make network plot ----------

l_g = create_layout(G, layout="fr")
l_g$x <- V(G)$x
l_g$y <- V(G)$y
l_g %<>% left_join(res_df)
ann_col <- c("pos"= "#F8766D", "neg" = "#619CFF")
title = 'Neuroblastoma - Metastasis vs. Non-metastasis'

# Colored by lineage
ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage, fill = lineage), shape=21,  alpha = 0.9) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(aes(label = cell_cluster), size = 3, repel = T, max.overlaps = 1000) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = title)
ggsave(filename = 'Network graph_color by lineage.pdf', width = 16, height = 15)

# Colored by sig_total
ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage), fill = "#FEEDC3", shape=21, alpha = 0.9) +
  geom_node_point(data = subset(l_g, pval_total < 0.05), aes(size= percentage, fill = log2fc_total_direction), shape=21, alpha = 0.9) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(aes(label = cell_cluster), size = 3, repel = T, max.overlaps = 1000) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = title)
ggsave(filename = 'Network graph_color by sig_total.pdf', width = 15, height = 15)

# Colored by sig_parent
ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(data = l_g, aes(size= percentage), fill = "#FEEDC3", shape=21, alpha = 0.9) +
  geom_node_point(data = subset(l_g, pval_parent < 0.05), aes(size= percentage, fill = log2fc_parent_direction), shape=21, alpha = 0.9) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(aes(label = cell_cluster), size = 3, repel = T, max.overlaps = 1000) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = title)
ggsave(filename = 'Network graph_color by sig_parent.pdf', width = 15, height = 15)

