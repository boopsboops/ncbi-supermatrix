#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Generating phylogenetic trees (may take a while) ...\n")

# get args
option_list <- list( 
    make_option(c("-m","--model"), type="character"),
    make_option(c("-v","--verbose"), type="character"),
    make_option(c("-e","--epsilon"), type="numeric"),
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$threads <- 4
#opt$model <- "TN93+G"
#opt$epsilon <- 0.1 
#opt$verbose <- "false"


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- read_csv(here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE)

# list fasta
fasta.files <- list.files(today.dir,pattern="*.aligned.trimmed.fasta$",recursive=FALSE,full.names=TRUE)
#print(fasta.files)

##### CHECK NUMBERS OF TAXA #####

# read files in and get length
files.fas <- purrr::map(fasta.files,ape::read.FASTA)
files.fas.length <- purrr::map(files.fas,length)

# filter length > 4
fasta.files.tree <- fasta.files[files.fas.length > 4]


##### RUN RAXML #####

# run on one
#xtr <- raxml_ng(file=fasta.files[6],model="TN93+G",maxthreads=8,epsilon=0.1)

# run on all 
purrr::walk(fasta.files.tree, \(x) raxml_ng(file=x,model=opt$model,maxthreads=opt$threads,epsilon=opt$epsilon,verbose=opt$verbose))

# print warning
if(length(fasta.files.tree) < length(fasta.files)) {
    # get names of length <= 4 
    fasta.files.err <- fasta.files[files.fas.length <= 4] |> basename() |> str_replace_all("\\.aligned\\.trimmed\\.fasta","")
    fasta.files.err.genes <- paste(fasta.files.err,collapse=", ")
    writeLines(glue("\nWARNING! The following loci had < 5 individuals, so gene trees could not be computed for these: {fasta.files.err.genes}\n",.trim=FALSE))
}
