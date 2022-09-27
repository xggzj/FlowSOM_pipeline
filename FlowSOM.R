# FlowSOM clustering

# Rscript FlowSOM.R -ncluster

# Parameters:
# - ncluster: number of FlowSOM cell clusters

# 0. Load packages ----------

library(FlowSOM)
library(flowCore)

# 0. Load data ----------

args = commandArgs(trailingOnly=TRUE) 
ncluster = as.integer(args[1][1]) 

ff = read.FCS('./output/cytof_preprocessed.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)

# 1. FlowSOM ----------

common_marker = colnames(ff@exprs) %>% unname()

# Run FlowSOM
fSOM = FlowSOM(ff,
               compensate = F,
               transform = F,
               scale = F,
               colsToUse = common_marker,
               nClus = round(sqrt(ncluster)), # max number of meta cluster is 90 (k in cutree)
               xdim = round(sqrt(ncluster)),
               ydim = round(sqrt(ncluster)),
               seed = 824)

save(fSOM, file = './output/fSOM.RData')



