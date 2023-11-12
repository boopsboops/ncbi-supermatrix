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

Here we search for our clade of interest (Ancistrus), and an outgroup. Because the search is clade based on taxonomic names, to add an outgroup we perform a second search and append the results to the previous one with the '-a' flag. We also first run a dry run with the '-d' flag to calculate an effective size of download batch.

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

# perform the dry runs to estimate batch size
scripts/download-sequences.R -c Ancistrus -n 500 -x 2500 -b 1 -a false -d true
scripts/download-sequences.R -c Lasiancistrus_schomburgkii -n 500 -x 2500 -b 1 -a false -d true

# now download the sequence data, setting batch size and appending our second search (Lasiancistrus) 
scripts/download-sequences.R -c Ancistrus -n 500 -x 2500 -b 100 -a false -d false
scripts/download-sequences.R -c Lasiancistrus_schomburgkii -n 500 -x 2500 -b 10 -a true -d false
```

Here we dereplicated and clean up the resulting fasta sequences downloaded from NCBI. The clustering groups the sequences into putative homologous loci based on their similarity, and regardless of their annotation. The script returns a table with a list of 

```bash
# flag '-n' is the maximum number of allowed missing data characters (Ns) in the sequence
#    Ns might indicate poor quality sequence data
# flag '-c' is the clustering threshold used to group the sequences into homologs
#    a lower value may mean multiple loci in the same cluster, and a high value may result in one locus split over multiple clusters 
scripts/clean-and-cluster.R -n 10 -c 0.6
```

```bash
# flag '-c' allows the user to pick clusters based on the numbered output of the 'clean-and-cluster.R' script
#    be sure to check that you are choosing the right cluster
# flag '-g' allows the use to give these cluster arbitrary names for reference
#    don't use spaces, commas, or other punctuation characters in the names
scripts/pick-clusters.R -c 9,5,7,0,10,8,11 -g cox1,cytb,rag1,rnas,rhod,rag2,myh6
```

```bash
# flag '-t' allows multithreading for obtaining metadata data from NCBI and FishBase
#    only required for very large search with thousands of species
#    do not use a '-t' value of more threads than is available on your machine
# flag '-c' chooses the database to use for taxonomic annotation
#    only 'ncbi'  and 'fishbase' are available
scripts/annotate-ncbi.R -t 1 -c fishbase
```

```bash
# flag '-n' is the length of species name
#    make some examples
scripts/filter-species.R -n 2
```

```bash
# flag '-p' is the minimum proportion of missing data per site
#    here sites with > 20% missing data are removed
# flag '-t' is multithreading for the mafft alignment
#    do not use a '-t' value of more threads than is available on your machine
scripts/align-trim-concatenate.R -p 0.2 -t 4
```

```bash
# flag '-m' is the phylogenetic model used by RAxML
#    here we use a simple model to generate the tree
# flag '-v' is the verbose option for RAxML
#    set to 'true' to see full output if running into errors
# flag '-e' is epsilon flag used used by RAxML
#    a smaller epsilon value gives a more thorough tree search 
# flag '-t' is multithreading for RAxML
#    do not use a '-t' value of more threads than is available on your machine
scripts/tree-search.R -m TN93+G -v false -e 0.1 -t 4
```

```bash
# flag '-w' is additional width proportion
#    if your tip labels are being cut off increase this value
#    or if too much white space, reduce
# flag '-h' is the height-width ratio
#    a ratio of 1.5 means the height will be 1.5 times the width
#    increase this value to stop tip labels getting bunched up
#    or decrease to make more compact
# flag '-s' is tree scaling factor (multiplicative)
#    bigger scaling factors are required for bigger trees
#    experiment with this to get the tree plotted suitably on the page
scripts/tree-plot.R -w 0.6 -h 1.5 -s 1
```
