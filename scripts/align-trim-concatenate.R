#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Aligning and trimming sequence data ...\n")

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


#### ALIGN AND TRIM FASTA #####

# get gene names and paths
genes <- ncbi.clean |> distinct(gene) |> pull()
genes.files <- here(today.dir,glue("{genes}.fasta"))

# write out
purrr::walk(genes, \(x) write_fasta(df=ncbi.clean,gene=x,dir=today.dir))

# align
purrr::walk(genes.files, \(x) align_fasta(infile=x,threads=opt$threads))

# trim
purrr::walk(genes.files, \(x) trim_fasta(infile=x,prop=opt$prop))


##### CONCATENATE AND WRITE OUT #####
 
# read in with ape
ali.files <- here(today.dir,glue("{genes}.aligned.trimmed.fasta"))
ali.all <- purrr::map(ali.files, \(x) read.FASTA(file=x))

# convert to matrix
ali.all.mat <- purrr::map(ali.all, as.matrix)

# concatenate
genes.concat <- as.list(do.call(cbind.DNAbin,args=c(ali.all.mat,fill.with.gaps=TRUE)))

# write out
genes.concat |> write.FASTA(file=here(today.dir,"concatenated-matrix.fasta"))
genes.concat |> write.dna(file=here(today.dir,"concatenated-matrix.phy"),format="sequential",nbcol=-1,colsep="")
genes.concat |> write.nexus.data(file=here(today.dir,"concatenated-matrix.nex"),interleaved=FALSE)

# make partion file and write out
partition_table(mat=ali.all.mat) |> write_tsv(here(today.dir,"concatenated-matrix.parts"),col_names=FALSE)

# file
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated-matrix.fasta'.",.trim=FALSE))
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated-matrix.nex'.",.trim=FALSE))
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated-matrix.phy'.",.trim=FALSE))
writeLines(glue("\nRAxML partitions file written to 'temp/{basename(today.dir)}/concatenated-matrix.parts'.\n",.trim=FALSE))


###### CLEAN UP #####

# make dir
if(!dir.exists(here(today.dir,"tempfiles"))) {dir.create(here(today.dir,"tempfiles"),recursive=TRUE)}

# list files
today.files <- list.files(today.dir,include.dirs=FALSE,recursive=TRUE,full.names=TRUE)

# grep ones we want
want <- "concatenated-matrix|trimmed.fasta|ncbi-clean"
files.from <- today.files[!str_detect(today.files,want)]

# add temp dir path
files.to <- here(today.dir,"tempfiles",basename(files.from))

invisible(file.rename(from=files.from,to=files.to))
