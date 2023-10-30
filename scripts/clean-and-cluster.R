#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("\nCleaning and clustering sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-n","--maxns"), type="numeric"),
    make_option(c("-c","--clustprop"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
fas.in <- here(today.dir,"genbank-dump.fasta")

# dereplicate
dereplicate_fasta(infile=fas.in)

# filter
filter_fasta(infile=fas.in,maxns=opt$maxns)

# cluster
cluster_fasta(infile=fas.in,identity=opt$clustprop)

# load clusters
clust.files <- list.files(here(today.dir),pattern="cluster\\.",full.names=TRUE)

# table of clusters
clust.tab <- tibble::as_tibble(file.info(clust.files),rownames="path") %>% 
    mutate(cluster=basename(path),size=size/1000) %>% 
    select(cluster,size,path) %>% 
    arrange(desc(size))

# add length data
clust.fas <- mapply(ape::read.FASTA,clust.files,USE.NAMES=TRUE,SIMPLIFY=TRUE)
clust.lens <- mapply(length,clust.fas,USE.NAMES=FALSE,SIMPLIFY=TRUE)
clust.lens.tab <- tibble(cluster=basename(names(clust.fas)),nseqs=clust.lens)
all.clusters <- left_join(clust.tab,clust.lens.tab,by=join_by(cluster))

# add description
first.acc <- mapply(function(x) labels(x)[[1]],clust.fas,USE.NAMES=TRUE,SIMPLIFY=TRUE)
ape.gb <- ape::read.GenBank(first.acc)
description.tab <- tibble(acc=labels(ape.gb),description=attr(ape.gb,"description")) %>% 
    left_join(tibble(path=names(first.acc),acc=unname(first.acc)),by=join_by(acc)) %>%
    mutate(description=str_split_fixed(description," ",2)[,2]) %>%
    select(description,acc,path)

# final table
all.clusters.desc <- all.clusters %>% 
    left_join(description.tab,by=join_by(path)) %>% 
    select(cluster,size,nseqs,description,acc,path)
write_csv(all.clusters.desc,file=here(today.dir,"clusters.csv"))

# print and write out table
writeLines("\nAll clusters with number sequences, file size in kb, and description of first sequence.\n")
writeLines(paste("\nTable written out to",here(today.dir,"clusters.csv")))
all.clusters.desc %>% select(-path) %>% knitr::kable()
