#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Generating phylogenetic trees (may take a while) ...\n")

# get args
option_list <- list( 
    make_option(c("-p","--prop"), type="numeric"),
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$prop <- 0.1
#opt$threads <- 4


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- read_csv(here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE)

# list fasta
fasta.files <- list.files(today.dir,pattern="*.fasta$",recursive=FALSE,full.names=TRUE)

##### #####


# run on one
#xtr <- raxml_ng(file=fasta.files[6],model="TN93+G",maxthreads=8,epsilon=0.1)

# run on all 
# write out
atr <- purrr::walk(fasta.files, \(x) raxml_ng(file=x,model="TN93+G",maxthreads=8,epsilon=0.1,verbose="false"))


tree.files <- glue("{fasta.files}.raxml.bestTree")



ncbi.clean.tips <- ncbi.clean |> distinct(scientificName,genus,family) |> mutate(tiplabel=glue("{family} | {str_replace_all(scientificName,'_',' ')}"))

# make tree plotting fun
ggtree_autoplot <- function(path,tb,scale.factor) {
    tr <- treeio::read.newick(path)
    tr <- phangorn::midpoint(tr)
    ntips <- length(tr$tip.label)
    tree.length <- max(castor::get_all_distances_to_root(tr)) + (max(castor::get_all_distances_to_root(tr)) * 0.1)
    tree.scale <- ntips/scale.factor
        p <- ggtree(tr, ladderize=TRUE,right=TRUE) %<+% tb
        pp <- p + geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE) +
        geom_nodelab(hjust=1.6,size=3) +
        geom_tippoint(aes(color=genus),size=2.5) +
        theme(legend.position="none") +
        xlim(0,tree.length)
    ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=(tree.length*(2000/tree.length)),height=(ntips*(4000/ntips)),units="px")#,scale=tree.scale
    #ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=1000,height=4000,units="px",scale=tree.scale)
    print(c(ntips,tree.length,tree.scale))
}

# BASIC
ggtree_autoplot <- function(path,tb,scale.factor) {
    tr <- treeio::read.newick(path)
    tr <- phangorn::midpoint(tr)
        p <- ggtree(tr, ladderize=TRUE,right=TRUE) %<+% tb
        pp <- p + geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE) +
        geom_nodelab(hjust=1.6,size=3) +
        geom_tippoint(aes(color=genus),size=2.5) +
        theme(legend.position="none")
    ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE)
}



#ggtree_autoplot(path=tree.files[1],tb=ncbi.clean.tips,scale.factor=40) 

purrr::walk(tree.files, \(x) ggtree_autoplot(path=x,tb=ncbi.clean.tips,scale.factor=40))

