#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("\nCleaning and clustering sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-n","--maxns"), type="numeric"),
    make_option(c("-c","--clustprop"), type="numeric"),
    make_option(c("-m","--minclust"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$minclust <- 2

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
clust.tab <- tibble::as_tibble(file.info(clust.files),rownames="path") |>
    mutate(cluster=basename(path),size=size/1000) |> 
    select(cluster,size,path) |> 
    arrange(desc(size))

# add data
clust.fas <- mapply(ape::read.FASTA,clust.files,USE.NAMES=TRUE,SIMPLIFY=TRUE)
# add description
first.acc <- mapply(function(x) labels(x)[[1]],clust.fas,USE.NAMES=TRUE,SIMPLIFY=TRUE)
first.acc.tab <- tibble(cluster=basename(names(clust.fas)),facc=unname(first.acc))
# add length
clust.lens <- mapply(length,clust.fas,USE.NAMES=FALSE,SIMPLIFY=TRUE)
clust.lens.tab <- tibble(cluster=basename(names(clust.fas)),nseqs=clust.lens)

# join
all.clusters <- clust.tab |> left_join(first.acc.tab,by=join_by(cluster)) |> left_join(clust.lens.tab,by=join_by(cluster))

# filter the clusters by minimum size and make vector
all.clusters.red <- all.clusters |> filter(nseqs>=opt$minclust)
facc <- all.clusters.red |> dplyr::pull(facc)

# chunk 200 should result in string of around 2200 chars
if(length(facc) < 200) {
    chunk <- length(facc) } else {
    chunk <- 200
}
chunk.frag <- unname(split(facc, ceiling(seq_along(facc)/chunk)))

# run chunk fun and combine
ape.gb.list <- mapply(FUN=read_genbank,chunk.frag,SIMPLIFY=FALSE,USE.NAMES=FALSE)
ape.gb <- bind_rows(ape.gb.list) |> mutate(description=str_split_fixed(description," ",2)[,2])

# final table
all.clusters.desc <- all.clusters.red |> 
    left_join(ape.gb,by=join_by(facc)) |> 
    select(cluster,size,nseqs,description,facc,path) |>
    arrange(desc(nseqs))
write_csv(all.clusters.desc,file=here(today.dir,"clusters.csv"))

# print and write out table
writeLines(glue("\nTop 20 clusters with >={opt$minclust} sequences, including file size in kb, n sequences, and description of first sequence.\n",.trim=FALSE))
writeLines(glue("\nTable written out to {here(today.dir,'clusters.csv')}"))
all.clusters.desc |> slice_head(n=20) |> select(-path) |> knitr::kable()
