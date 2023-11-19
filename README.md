[![DOI](https://zenodo.org/badge/689690400.svg)](https://zenodo.org/doi/10.5281/zenodo.10157383)

# NCBI Supermatrix

This is pipeline that will assemble phylogenetic 'supermatrix' datasets directly from the [NCBI nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide) database. Currently working on Ubuntu Linux.


### INSTALLATION

First ensure that [git v2.x](https://git-scm.com/downloads), [R v4.3.x](https://cloud.r-project.org/), [vsearch v2.x](https://github.com/torognes/vsearch), [trimal v1.4](https://github.com/inab/trimal), [mafft v7.x](https://mafft.cbrc.jp/alignment/software/) and [raxml-ng v1.2](https://github.com/amkozlov/raxml-ng) are installed and available on your `$PATH`. Newer versions of these programs will likely work, but R needs to be v4.3. Both git, mafft and vsearch can be installed from Ubuntu Linux repositories with `sudo apt install git mafft vsearch` or via Homebrew for Mac with `brew install git mafft vsearch`. The others will need to be installed via instructions on their respective GitHub pages. Check they are all installed by running e.g. `which mafft` or `which vsearch` in your home directory. The [renv-installer](https://github.com/jcrodriguez1989/renv-installer) software is useful for managing multiple R versions across projects.

To install NCBI-Supermatrix:

```
# clone the repository 
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix
# install the R packages (... takes a LONG time ...)
Rscript -e "renv::restore()"
```


### GENERAL NOTES

The dated results output directory (e.g. `temp/Results_2023-11-19`) is created in directory `temp/` using the date that the `download-sequences.R` script was initiated. All subsequent scripts use this most recent directory in `temp/` irrespective of day. Be aware of this if deleting directories. Always make backups.

Unpredictable results may occur when running scripts repeatedly with different settings, or out of order. Clusters may be created that are now longer required, for example. The scripts will often scan directories for files of a certain type and use all of them. Be aware of this, and if in doubt, delete files manually and start at the beginning for a new analysis.


### RUNNING

#### download-sequences.R

With the `download-sequences.R` script we search for our clade of interest (catfishes family Akysidae), and an outgroup (_Amblyceps mangois_). Because the search is clade-based using taxonomic names, to add an outgroup we perform a second search and append the results to the previous one with the '-a' flag. To estimate an appropriate maximum download batch size we can make a dry run with the '-d' flag.

```bash
# flag '-c' [character] is the clade of interest
#    this clade can be any taxonomic level present in the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy)
#    in species names the space must be replaced with an underscore
# flag '-n' [integer] is the minimum length of the sequence searched for (bp)
#    about 400-500 is best for Sanger-sequenced markers
# flag '-x' [integer] is the maximum length of the sequence searched for (bp)
#    about 1500-2500 is best for Sanger-sequenced markers
# flag '-b' [integer] is the maximum batch size of sequences to download from NCBI Entrez
#    a batch size of 100 will divide the sequences in batches of < 100 
#    this number must be smaller than the number of hits, but no larger than 9999
# flag '-a' [logical] toggles appending of data onto a previous file
#    use this to perform multiple searches, e.g. for outgroups
#    the 'false' option will overwrite previous files, while 'true' will add data
# flag '-d' [logical] is the dry run option that does not download any sequence data
#    use this to estimate the batch size

# perform the dry runs with '-d' to estimate batch size
scripts/download-sequences.R -c Akysidae -n 500 -x 2500 -b 1 -a false -d true
scripts/download-sequences.R -c Amblyceps_mangois -n 500 -x 2500 -b 1 -a false -d true

# now download the sequence data, setting maximum batch size to 30 and appending our outgroup Amblyceps mangois 
scripts/download-sequences.R -c Akysidae -n 500 -x 2500 -b 30 -a false -d false
scripts/download-sequences.R -c Amblyceps_mangois -n 500 -x 2500 -b 30 -a true -d false
```


#### clean-and-cluster.R

The `clean-and-cluster.R` script dereplicates and clean up the resulting fasta sequences downloaded from NCBI. The clustering groups the sequences into putative homologous loci based on their similarity, and regardless of their annotation. The script returns a table with a list of clusters and their size. If this step results in error run `rm cluster*` in the `temp/Results_today` directory and try again. To save computation time, use the '-m' flag to remove clusters with low numbers of sequences.

```bash
# flag '-n' [integer] is the maximum number of allowed missing data characters (Ns) in the sequence
#    Ns might indicate poor quality sequence data
#    here is arbitrarily set to 10
# flag '-c' [numeric] is the clustering threshold used to group the sequences into homologs
#    a lower value may mean multiple loci in the same cluster, 
#    a high value may result in one locus split over multiple clusters
#    a value of 0.6 appears to work well
# flag '-m' [integer] is the minimum retained cluster size
#    all clusters with value of '-m' or greater are retained
#    to remove singleton clusters use a value of 2
#    to retain all clusters use a value of 1
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2
```


#### pick-clusters.R

The `pick-clusters.R` script allows the user to pick loci for the analysis from the clustering table generated from the previous script. The '-c' flag is the cluster number in the table, and the '-g' flag are arbitrary names the user assigns for reference. Here we pick clusters 6, 3, and 2, which correspond to the cox1, rag2, and cytb loci.

```bash
# flag '-c' [character] allows the user to pick clusters
#    see table output from the 'clean-and-cluster.R' script
#    be sure to check that you are choosing the right cluster
#    separate the numbers with commas (no spaces)
# flag '-g' [character] allows the use to give these cluster arbitrary names for reference
#    separate the names with commas (no spaces)
#    don't use spaces, commas, or other punctuation characters in the names
scripts/pick-clusters.R -c 6,3,2 -g cox1,rag2,cytb
```


#### annotate-ncbi.R

The `annotate-ncbi.R` script adds additional higher taxonomic information from either FishBase or NCBI.

```bash
# flag '-t' [integer] allows multithreading for obtaining metadata data from NCBI and FishBase
#    only required for very large search with thousands of species
#    do not use a '-t' value of more threads than is available on your machine
# flag '-c' [character] chooses the database to use for taxonomic annotation
#    only 'ncbi'  and 'fishbase' are available
scripts/annotate-ncbi.R -t 1 -c fishbase
```


##### filter-species.R

The `filter-species.R` cleans up the species names and picks the longest sequence for each of the loci. The script also allows you to choose how many elements are in the species names (based on spaces), which is important if dealing with subspecies or undescribed species. For example "Akysis sp. INHS 93579" would be cut down to "Akysis sp.INHS" with '-n 2' and "Akysis sp.INHS 93579" with '-n 3'. If you don't care about lumping subspecies or undescribed species, then use '-n 2'. Check the `ncbi-clean.csv` table to see how names have been edited. The script can also exclude sequences that you do not want using the `assets/exclusions.csv`. Just add any sequences to this file and they are automatically removed next time the script is run. Use the `ncbi-clean.csv` table to get the accession numbers. 

```bash
# flag '-n' [integer] is the length of the number of elements in the species name (delimited by spaces)
#    the spaces after 'sp.' or 'aff.' or 'cf.' are automatically removed
#    '-n 2' will turn 'Akysis sp. INHS 93579' into 'Akysis sp.INHS'
#    '-n 3' will turn 'Akysis sp. INHS 93579' into 'Akysis sp.INHS 93579'
scripts/filter-species.R -n 3
```


#### align-trim-concatenate.R
The `align-trim-concatenate.R` script aligns each locus seperately and trims missing data to a user specified proportion. Here a value of 0.2 means that sites (alignment columns) that have less than 20% data are removed (i.e. comprise > 80% missing data). This is to remove poorly aligned regions and increase the efficiency of the phylogenetic analysis.

```bash
# flag '-p' [proportion] is the minimum proportion of data per site
#    use any value between 0 and 1
#    '-p 0.2' means sites with < 20% data are removed
# flag '-t' [integer] is multithreading for the mafft alignment
#    do not use a '-t' value of more threads than is available on your machine
scripts/align-trim-concatenate.R -p 0.2 -t 4
```


#### tree-search.R

The `tree-search.R` script runs RAxML Next Generation to generate a tree for the individual genes and for the concatenated matrix. Please refer to the raxml-ng manual for help understanding models and epsilon values.

```bash
# flag '-m' [character] is the phylogenetic model used by RAxML
#    here we use a simple model to generate the tree
# flag '-v' [logical] is the verbose option for RAxML
#    set to 'true' to see full output if running into errors
# flag '-e' [numeric] is epsilon flag used used by RAxML
#    a smaller epsilon value such as 0.01 gives a more thorough tree search
#    a value of 10 is much faster 
# flag '-t' [integer] is multithreading for RAxML
#    do not use a '-t' value of more threads than is available on your machine
scripts/tree-search.R -m TN93+G -v false -e 10 -t 4
```


#### tree-plot.R

The `tree-plot.R` script uses the ggtree R package to plot the trees as PDF. It will require some experimentation to get a tree looking how you like it, and settings may need to be changed for different loci, depending on number of taxa and tree length. Once you are happy with how a tree looks, make a copy to prevent it being overwritten.

```bash
# flag '-w' [proportion] is additional width proportion
#    if your tip labels are being cut off increase this value
#    or if too much white space, reduce
# flag '-h' [proportion] is the height-width ratio
#    a ratio of 0.5 means the height will be 0.5 times (i.e. half) the width
#    increase this value to stop tip labels getting bunched up
#    or decrease to make more compact
# flag '-s' [proportion] is tree scaling factor (multiplicative)
#    bigger scaling factors are required for bigger trees
#    experiment with this to get the tree plotted suitably on the page
scripts/tree-plot.R -w 0.6 -h 0.5 -s 1
```


#### tidy-results-directory.R

The `tidy-results-directory.R` script tidies the results directory, as there can be a lot of intermediate files created. The complete directory is also backed up in `backups` so can be retrieved if needed, or deleted to save disk space.

```bash
# clean up and tidy the results directory
scripts/tidy-results-directory.R
```
