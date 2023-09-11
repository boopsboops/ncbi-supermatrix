#!/usr/bin/env Rscript

# LOAD PACKAGES
suppressMessages({
    library("here")
    library("rentrez")
    library("ape")
    library("tidyverse")
    library("optparse")
    library("crul")
    library("traits")
    library("parallel")
})


# ENTREZ PARALLEL SEARCH
entrez_fetch_parallel <- function(webhist,chunks,file) {
    chunk.first <- dplyr::first(chunks)
    chunk.last <- dplyr::last(chunks)
    chunk.len <- length(chunks)
    writeLines(paste0("Downloading hits ",chunk.first," to ",chunk.last," from NCBI ..."))
    fetch.fas <- rentrez::entrez_fetch(db="nucleotide",rettype="fasta",web_history=webhist,retstart=chunk.first,retmax=chunk.len)
    write(fetch.fas,file=file,append=TRUE)
    Sys.sleep(3)
    #invisible()
}


# ENTREZ DOWNLOAD
entrez_download <- function(clade,minlen,maxlen,batchsize,fasout) {
    file.create(fasout,overwrite=TRUE)
    rtm <- as.integer(batchsize)
    if (rtm > 9999) {stop("\nError! Batch size can be no larger than 9,999.\n")}
    es.res <- rentrez::entrez_search(db="nucleotide",term=paste0(clade,"[ORGANISM] AND ",minlen,":",maxlen,"[SLEN]"),retmax=999999,use_history=TRUE)
    if (rtm > es.res$count) {stop("\nError! Batch size value is larger than number of hits. Please use a smaller value.\n")}
    writeLines(paste("\nNCBI search reports",length(es.res$ids),"hits for",clade,"...\n"))
    if (length(es.res$ids) > 999999) {stop("\nError! You have more hits than can be downloaded. Please break up the search into smaller batches.\n")}
    es.post <- entrez_post(db="nucleotide",web_history=es.res$web_history)
    query.split <- split(0:es.res$count,cut(0:es.res$count,ceiling(es.res$count/rtm),labels=FALSE))
    x <- mapply(function(x) entrez_fetch_parallel(webhist=es.post,chunks=x,file=fasout), x=query.split,SIMPLIFY=TRUE,USE.NAMES=FALSE)
    es.fas <- ape::read.FASTA(fasout)
    #if (length(es.fas) != es.res$count) {stop("Error! Not all sequences were downloaded. Please try again.")}
    writeLines(paste("\nDone!",length(es.fas),"sequences downloaded in FASTA format.\n"))
    return(es.fas)
}


# FUNCTION TO RUN PARALLEL NCBI_BYID WITH TIMEOUT AND REPEAT
ncbi_byid_parallel <- function(accs){
    start_time <- Sys.time()
    Sys.sleep(time=runif(n=1,min=0,max=3))
    crul::set_opts(http_version=2)
    ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    if(class(ncbi.tab)!="data.frame") {
        #writeLines("Error found! Repeating ...")
        Sys.sleep(time=3)
        crul::set_opts(http_version=2)
        ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    } else {
        ncbi.tab <- ncbi.tab
    }
    if(class(ncbi.tab)!="data.frame") {
        stop(writeLines("Searches failed ... aborted")) 
    } else {
        end_time <- Sys.time()
        writeLines(paste0("Metadata for ",length(accs)," accessions downloaded (starting ",accs[1],"). Download took ",round(as.numeric(end_time-start_time),digits=2)," seconds."))
        return(ncbi.tab)
    }
}


# DEREPLICATE FASTA
dereplicate_fasta <- function(infile) {
    outfile <- str_replace_all(infile,"\\.fasta",".derep.fasta")
    exe.string <- paste("vsearch --threads 0 --fastx_uniques",infile,"--fastaout",outfile)
    system(exe.string)
}


# FILTER FASTA
filter_fasta <- function(infile,maxns) {
    infile <- str_replace_all(infile,"\\.fasta",".derep.fasta")
    outfile <- str_replace_all(infile,"\\.fasta",".filtered.fasta")
    exe.string <- paste("vsearch --threads 0 -fastx_filter",infile,"--fastq_maxns",maxns,"--fastaout",outfile)
    system(exe.string)
}


# CLUSTER FASTA
cluster_fasta <- function(infile,identity) {
    infile <- str_replace_all(infile,"\\.fasta",".derep.filtered.fasta")
    outfile <- here(today.dir,"cluster.")
    str_replace_all(infile,"genbank-dump.derep.filtered.fasta","cluster.")
    exe.string <- paste("vsearch --threads 0 --cluster_fast",infile,"--id",identity,"--clusters",outfile)
    system(exe.string)
}


# LOAD AND RENAME FASTA
load_and_rename <- function(gene,path,wd) {
    fas <- ape::read.FASTA(path)
    path.gene <- here(today.dir,paste("gene",gene,"fasta",sep="."))
    ape::write.FASTA(fas,path.gene)
}


# TABULATE GENES
tabulate_genes <- function(path) {
    fas <- ape::read.FASTA(path)
    accs <- names(fas)
    gene <- str_split_fixed(basename(path),"\\.",3)[,2]
    accs.gene <- tibble(acc=accs,gene=gene)
    return(accs.gene)
}

# REPORT
writeLines("\nPackages and functions loaded.\n")
