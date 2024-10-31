#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Aligning and trimming sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-p","--prop"), type="numeric"),
    make_option(c("-t","--threads"), type="numeric"),
    make_option(c("-i","--indiv"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$prop <- 0.2
#opt$threads <- 4


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- read_csv(here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE,col_types=cols(.default=col_character()))


#### ALIGN AND TRIM FASTA #####

# get gene names and paths
genes <- ncbi.clean |> distinct(gene) |> pull()
genes.files <- here(today.dir,glue("{genes}.fasta"))

# write out
purrr::walk(genes, \(x) write_fasta(df=ncbi.clean,genez=x,dir=today.dir,pop=as.logical(opt$indiv)))

# align
purrr::walk(genes.files, \(x) align_fasta(infile=x,threads=opt$threads))

# trim
purrr::walk(genes.files, \(x) trim_fasta(infile=x,prop=opt$prop))


#### STOP IF POP ####

# end the script if in pop mode
if(opt$indiv == "true") {
   cli::cli_alert_success("Finished running in pop mode.\f")
   quit(save="no")
}


##### CONCATENATE AND WRITE OUT #####
 
# read in with ape
ali.files <- here(today.dir,glue("{genes}.aligned.trimmed.fasta"))
ali.all <- purrr::map(ali.files, \(x) read.FASTA(file=x))

# convert to matrix
ali.all.mat <- purrr::map(ali.all, as.matrix)

# concatenate
genes.concat <- as.list(do.call(cbind.DNAbin,args=c(ali.all.mat,fill.with.gaps=TRUE)))

# write out
genes.concat |> write.FASTA(file=here(today.dir,"concatenated.aligned.trimmed.fasta"))
genes.concat |> write.dna(file=here(today.dir,"concatenated.aligned.trimmed.phy"),format="sequential",nbcol=-1,colsep="")
genes.concat |> write.nexus.data(file=here(today.dir,"concatenated.aligned.trimmed.nex"),interleaved=FALSE)

# make partion file and write out
partition_table(mat=ali.all.mat) |> write_tsv(here(today.dir,"concatenated.aligned.trimmed.parts"),col_names=FALSE)

# file
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated.aligned.trimmed.fasta'.",.trim=FALSE))
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated.aligned.trimmed.nex'.",.trim=FALSE))
writeLines(glue("\nConcatenated matrix written to 'temp/{basename(today.dir)}/concatenated.aligned.trimmed.phy'.",.trim=FALSE))
writeLines(glue("\nRAxML partitions file written to 'temp/{basename(today.dir)}/concatenated.aligned.trimmed.parts'.\n",.trim=FALSE))
