#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# make date dir for results
today.dir <- here::here("temp",paste0("Results_",Sys.Date()))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}

# info
writeLines("\nSearching GenBank for sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-c","--clade"), type="character"),
    make_option(c("-n","--minlen"), type="numeric"),
    make_option(c("-x","--maxlen"), type="numeric"),
    make_option(c("-b","--batch"), type="numeric"),
    make_option(c("-a","--append"), type="character"),
    make_option(c("-d","--dry"), type="character")    
    )

# set args
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list,add_help_option=FALSE))

# if an underscore then write out as an outgroup
if(stringr::str_detect(opt$clade,"_")) {
    writeLines(opt$clade, here::here(today.dir,"outgroup.txt"))
}

# edit clade to remove underscore
opt$clade <- stringr::str_replace_all(opt$clade,"_"," ")

# run entrez
seqs.fas <- entrez_download(clade=opt$clade,minlen=opt$minlen,maxlen=opt$maxlen,batchsize=opt$batch,append=opt$append,dry=opt$dry,fasout=here::here(today.dir,"genbank-dump.fasta"))

# remove annoying TRUE file
if (file.exists(here::here("TRUE"))) {invisible(file.remove(here::here("TRUE")))}
