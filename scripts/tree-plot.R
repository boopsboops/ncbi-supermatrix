#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Plotting phylogenetic trees ...\n")

# get args
option_list <- list( 
    make_option(c("-s","--scalefactor"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- read_csv(here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE)

# list trees
tree.files <- list.files(today.dir,pattern="*.fasta.raxml.bestTree$",recursive=FALSE,full.names=TRUE)


##### AUTO PLOT #####

# make tip labels
ncbi.clean.tips <- ncbi.clean |> distinct(scientificName,genus,family) |> mutate(tiplabel=glue("{family} | {str_replace_all(scientificName,'_',' ')}"))

# make tree plotting fun
ggtree_autoplot <- function(path,tb,scale.factor) {
    tr <- treeio::read.newick(path)
    tr <- phangorn::midpoint(tr)
    ntips <- length(tr$tip.label)
    tree.length <- max(castor::get_all_distances_to_root(tr)) + (max(castor::get_all_distances_to_root(tr)) * 0.5)
    tree.scale <- ntips/scale.factor
        p <- ggtree(tr, ladderize=TRUE,right=TRUE) %<+% tb
        pp <- p + geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE) +
        geom_nodelab(hjust=1.6,size=3) +
        geom_tippoint(aes(color=genus),size=2.5) +
        theme(legend.position="none") +
        xlim(0,tree.length)
    #ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=(tree.length*(2000/tree.length)),height=(ntips*(4000/ntips)),units="px")#,scale=tree.scale
    ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=210,height=297,units="mm",scale=tree.scale)
    #print(c(ntips,tree.length,tree.scale))
}

# BASIC
#ggtree_autoplot <- function(path,tb,scale.factor) {
#    tr <- treeio::read.newick(path)
#    tr <- phangorn::midpoint(tr)
#        p <- ggtree(tr, ladderize=TRUE,right=TRUE) %<+% tb
#        pp <- p + geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE) +
#        geom_nodelab(hjust=1.6,size=3) +
#        geom_tippoint(aes(color=genus),size=2.5) +
#        theme(legend.position="none")
#    ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE)
#}

#ggtree_autoplot(path=tree.files[1],tb=ncbi.clean.tips,scale.factor=40) 

purrr::walk(tree.files, \(x) ggtree_autoplot(path=x,tb=ncbi.clean.tips,scale.factor=opt$scalefactor))

# print
writeLines(glue("\nTrees written to:\n",.trim=FALSE))
writeLines(glue("{basename(tree.files)}.pdf\n",.trim=FALSE))
