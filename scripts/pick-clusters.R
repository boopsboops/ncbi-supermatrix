#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("\nPicking clusters ...\n")

# get args
option_list <- list( 
    make_option(c("-c","--clusters"), type="character"),
    make_option(c("-g","--genes"), type="character")
    #make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# split the input
clusters.chosen <- unlist(str_split(opt$clusters,","))
genes.chosen <- unlist(str_split(opt$genes,","))
#clusters.chosen <- c("8","4","6")
#genes.chosen <- c("cox1","cytb","rag1")
# scripts/pick-clusters.R -c 8,4,6 -g cox1,cytb,rag1

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]

# load csv
clusters <- read_csv(here(today.dir,"clusters.csv"),show_col_types=FALSE)

# make cluster names
cluster.gene <- tibble(cluster=paste("cluster",clusters.chosen,sep="."),gene=genes.chosen)

# subset
cluster.gene.path <- cluster.gene %>% left_join(clusters,by=join_by(cluster))

# rename fasta files
invisible(mapply(function(x,y) load_and_rename(gene=x,path=y,wd=today.dir), x=pull(cluster.gene.path,gene), y=pull(cluster.gene.path,path), SIMPLIFY=TRUE,USE.NAMES=FALSE))

# print 
writeLines(paste("\nTotal",length(genes.chosen),"gene clusters written to:\n"))
writeLines(here(today.dir,list.files(today.dir,"gene\\.")))
