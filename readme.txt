Section 1: CyTOF preprocessing ----------

Input:

- A table of file paths: one sample per row, the first column is path to the FCS file, the second column is path to the Grid annotation file.
- Threshold of outlier removal - n: for each marker, remove top n% cells with highest expression.
- A vector of markers to be ignored during outlier removal (example: ignore CD123 and gdTCR to keep pDC and Basophils).
- A logical value: do ComBat batch correction or not.

Output: 

- A fcs file of integrated preprocessed CyTOF data. 
- A csv file of integrated Grid annotation (cells in the same order of the integrated fcs file).

Optional output: 
- A vector of markers to be used for FlowSOM clustering.
- Number of cells in total before and after preprocessing.
- A table of the distribution of expression level for each marker after preprocessing: one marker per row, 0, 0.25, 0.5, 0.75, 0.99 and 1 percentiles in each column respectively.
- A UMAP plot showing batch effects before and after preprocessing.




Section 2: FlowSOM clustering ----------

Input: 

- A fcs file of integrated preprocessed CyTOF data. 
- A csv file of integrated Grid annotation (cells in the same order of the integrated fcs file).
- Number of FlowSOM cell clusters.

Output: 
- A table of median marker expression of each cell cluster: each row is a FlowSOM cell cluster, columns are marker name (median expression value), cell_cluster (cluster number), n (number of cells in each cluster), percentage (the proportion of cells in each cluster in the total cells).
- A table of relative frequency of each cell cluster in each sample: one sample per row, cell clusters in columns.
- A pdf file of FlowSOM results summary.




Section 3: Making plots ----------

Input:

- A table of median marker expression of each cell cluster: each row is a FlowSOM cell cluster, columns are marker name (median expression value), cell_cluster (cluster number), n (number of cells in each cluster), percentage (the proportion of cells in each cluster in the total cells).
- Filtering threshold of edges in the graph (for each node only the edges with rank less than the threshold are retained).

Output:

- A GraphML file of the graph.
- A heatmap of FlowSOM cell cluster phenotypes.
- A force-directed graph.

