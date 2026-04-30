#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# info
#writeLines("\nCleaning and clustering sequence data ...\n")
cli_report(txt="Running 'clean-and-cluster.R' ... Cleaning and clustering sequence data ...",rule=FALSE,alert="info")

# get args
option_list <- list( 
    make_option(c("-n","--maxns"), type="numeric"),
    make_option(c("-c","--clustprop"), type="numeric"),
    make_option(c("-m","--minclust"), type="numeric"),
    make_option(c("-d","--derep"), type="character"),
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$minclust <- 2

# get latest dir
today.dir <- sort(grep("/Results_",list.dirs(here::here("temp"),recursive=FALSE),value=TRUE),decreasing=TRUE)[1]
fas.in <- here::here(today.dir,"genbank-dump.fasta")

# dereplicate
dereplicate_fasta(infile=fas.in,dereplicate=opt$derep,threads=opt$threads)

# filter
filter_fasta(infile=fas.in,maxns=opt$maxns,threads=opt$threads)

# cluster
cluster_fasta(infile=fas.in,identity=opt$clustprop,threads=opt$threads)

# load clusters
clust.files <- list.files(here::here(today.dir),pattern="cluster\\.",full.names=TRUE)

# table of clusters
clust.tab <- tibble::as_tibble(file.info(clust.files),rownames="path") |>
    dplyr::mutate(cluster=basename(path),size=size/1000) |> 
    dplyr::select(cluster,size,path) |> 
    dplyr::arrange(desc(size))

# add data
clust.fas <- mapply(ape::read.FASTA,clust.files,USE.NAMES=TRUE,SIMPLIFY=TRUE)
# add description
first.acc <- mapply(function(x) labels(x)[[1]],clust.fas,USE.NAMES=TRUE,SIMPLIFY=TRUE)
first.acc.tab <- tibble::tibble(cluster=basename(names(clust.fas)),facc=unname(first.acc))
# add length
clust.lens <- mapply(length,clust.fas,USE.NAMES=FALSE,SIMPLIFY=TRUE)
clust.lens.tab <- tibble::tibble(cluster=basename(names(clust.fas)),nseqs=clust.lens)

# join
all.clusters <- clust.tab |> dplyr::left_join(first.acc.tab,by=join_by(cluster)) |> dplyr::left_join(clust.lens.tab,by=join_by(cluster))

# filter the clusters by minimum size and make vector
all.clusters.red <- all.clusters |> dplyr::filter(nseqs>=opt$minclust)
facc <- all.clusters.red |> dplyr::pull(facc)

# chunk 200 should result in string of around 2200 chars
if(length(facc) < 200) {
    chunk <- length(facc) } else {
    chunk <- 200
}
chunk.frag <- unname(split(facc, ceiling(seq_along(facc)/chunk)))

# run chunk fun and combine
ape.gb.list <- mapply(FUN=read_genbank,chunk.frag,SIMPLIFY=FALSE,USE.NAMES=FALSE)
ape.gb <- dplyr::bind_rows(ape.gb.list) |> dplyr::mutate(description=str_split_fixed(description," ",2)[,2])

# final table
all.clusters.desc <- all.clusters.red |> 
    dplyr::left_join(ape.gb,by=join_by(facc)) |> 
    dplyr::select(cluster,size,nseqs,description,facc,path) |>
    dplyr::arrange(desc(nseqs))
readr::write_csv(all.clusters.desc,file=here::here(today.dir,"clusters.csv"))

# print and write out table
cli_report(txt=glue::glue("Top 20 clusters with >={opt$minclust} sequences, including file size in kb, n sequences, and description of first sequence."),rule=FALSE,alert="info")
#writeLines(glue::glue("\nTop 20 clusters with >={opt$minclust} sequences, including file size in kb, n sequences, and description of first sequence.\n",.trim=FALSE))
cli_report(txt=glue::glue("Table written out to {here(today.dir,'clusters.csv')}"),rule=FALSE,alert="info")
#writeLines(glue::glue("\nTable written out to {here(today.dir,'clusters.csv')}"))
all.clusters.desc |> dplyr::slice_head(n=20) |> dplyr::select(-path) |> knitr::kable()

# report
cli_report(txt="Clustering step completed.",rule=TRUE,alert="success")
