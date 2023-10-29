#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# make date dir for results
today.dir <- here("temp",paste0("Results_",Sys.Date()))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}

# info
writeLines("\nSearching GenBank for sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-c","--clade"), type="character"),
    make_option(c("-m","--minlen"), type="numeric"),
    make_option(c("-x","--maxlen"), type="numeric"),
    make_option(c("-b","--batch"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# run entrez
seqs.fas <- entrez_download(clade=opt$clade,minlen=opt$minlen,maxlen=opt$maxlen,batchsize=opt$batch,fasout=here(today.dir,"genbank-dump.fasta"))
