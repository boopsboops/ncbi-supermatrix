# ncbi-supermatrix


```bash
# download - need to add the -a to append another search
scripts/download-sequences.R -c Ancistrus -n 500 -x 1500 -b 100 -a false
scripts/download-sequences.R -c Lasiancistrus -n 500 -x 1500 -b 50 -a true

# dereplicate and cluster
scripts/clean-and-cluster.R -n 10 -c 0.6

# pick the clusters
scripts/pick-clusters.R -c 8,4,6,7,9 -g cox1,cytb,rag1,rag2,rhod

# annotate the ncbi data with fishbase
scripts/annotate-ncbi.R -t 1 -c fishbase

# clean up the data and filter the species to one indiv 
scripts/filter-species.R

# align trim and concatentate
scripts/align-trim-concatenate.R -p 0.1 -t 4

```