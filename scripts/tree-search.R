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
fasta.files <- list.files(today.dir,pattern="*.fasta$",recursive=FALSE,full.names=TRUE)


##### RUN RAXML #####

# run on one
#xtr <- raxml_ng(file=fasta.files[6],model="TN93+G",maxthreads=8,epsilon=0.1)

# run on all 
purrr::walk(fasta.files, \(x) raxml_ng(file=x,model=opt$model,maxthreads=opt$threads,epsilon=opt$epsilon,verbose=opt$verbose))


##### CLEAN UP #####

# list files
rax.files <- list.files(today.dir,pattern="*.raxml.",include.dirs=FALSE,recursive=TRUE,full.names=TRUE)

# grep ones we want
want <- "raxml.bestTree$"
files.from <- rax.files[!str_detect(rax.files,want)]

# add temp dir path
if(!dir.exists(here(today.dir,"tempfiles"))) {dir.create(here(today.dir,"tempfiles"),recursive=TRUE)}
files.to <- here(today.dir,"tempfiles",basename(files.from))

# mv
invisible(file.rename(from=files.from,to=files.to))
