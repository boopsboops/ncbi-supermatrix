# ncbi-supermatrix


```bash
# download - need to add the -a to append another search
scripts/download-sequences.R -c Ancistrus -m 500 -x 1500 -b 100 -a false

# dereplicate and cluster
scripts/clean-and-cluster.R -m -p

# pick the clusters
scripts/pick-clusters.R -c 8,4,6 -g cox1,cytb,rag1

# annotate the ncbi data with fishbase
scripts/annotate-ncbi.R -t 1 -c fishbase

# clean up the data and filter the species to one indiv 
scripts/filter-species.R

# align trim and concatentate
scripts/align-trim-concatenate.R -p 0.1 -t 4

```