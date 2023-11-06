# ncbi-supermatrix

This is pipeline that will assemble phylogenetic 'supermatrix' datasets directly from NCBI. Currently working on Ubuntu Linux.

### INSTALLATION

First ensure that [vsearch](https://github.com/torognes/vsearch), [trimal](https://github.com/inab/trimal), [mafft](https://mafft.cbrc.jp/alignment/software/) and [raxml-ng](https://github.com/amkozlov/raxml-ng) are installed and on your `$PATH`. Both mafft and vsearch can be installed from Ubuntu repositories with `sudo apt install mafft vsearch`. The others will need to be compiled. Check they are all installed by running e.g. `which mafft` or `which vsearch` in your home directory.x

```
# clone the repository 
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix
# install the R packages (takes a LONG time)
Rscript -e "renv::restore()"
```

### RUNNING

Here we search for our clade of interest (Ancistrus)

```bash
# flag '-c' [name] is the clade of interest
#    this clade can be any taxonomic level present in the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy)
#    in species names the space must be replaced with an underscore
# flag '-n' [integer] is the minimum length of the sequence searched for (bp)
#    about 400-500 is best for Sanger sequenced markers
# flag '-x' [integer] is the maximum length of the sequence searched for (bp)
#    about 1500-2500 is best for Sanger sequenced markers
# flag '-b' [integer] is the batch size of sequences to download from NCBI Entrez.
#    a batch size of 100 will download the sequences in batches of 100 or divide them 
#    this number must be smaller than the number of hits, but no larger than 9999
# flag '-a' [true/false] toggles appending of data onto a previous file.
#    use this to perform multiple searches, e.g. for outgroups
#    the 'false' option will overwrite previous files, while 'true' will add data
# flag '-d' [true/false] is the dry run option that does not download any sequence data
#    use this to estimate the batch size
scripts/download-sequences.R -c Ancistrus -n 500 -x 2500 -b 1 -a false -d true
scripts/download-sequences.R -c Lasiancistrus_schomburgkii -n 500 -x 2500 -b 1 -a false -d true

```


```bash
# dereplicate and cluster
scripts/clean-and-cluster.R -n 10 -c 0.6

# pick the clusters
scripts/pick-clusters.R -c 8,4,6,7 -g cox1,cytb,rag1,rag2

# annotate the ncbi data with fishbase
scripts/annotate-ncbi.R -t 1 -c fishbase

# clean up the data and filter the species to one indiv 
scripts/filter-species.R

# align trim and concatentate
scripts/align-trim-concatenate.R -p 0.2 -t 4

# run raxml tree search
scripts/tree-search.R -m TN93+G -v false -e 0.1 -t 4

# plot tree pdfs
scripts/tree-plot.R -s 50

```