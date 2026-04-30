#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
#writeLines("Aligning and trimming sequence data ...\n")
cli_report(txt="Running 'align-trim-concatenate.R' ... Aligning and trimming sequence data ...",rule=FALSE,alert="info")

# get args
option_list <- list( 
    make_option(c("-p","--prop"), type="numeric"),
    make_option(c("-t","--threads"), type="numeric"),
    make_option(c("-i","--indiv"), type="character")
    )

# set args
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$prop <- 0.2
#opt$threads <- 4


##### LOAD DATA #####

# get latest dir
today.dir <- sort(grep("/Results_",list.dirs(here::here("temp"),recursive=FALSE),value=TRUE),decreasing=TRUE)[1]
cli_report(txt=glue::glue("Working in directory 'temp/{basename(today.dir)}'."),rule=FALSE,alert="info")
#writeLines(glue::glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- readr::read_csv(here::here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE,col_types=cols(.default=col_character()))

#### ALIGN AND TRIM FASTA #####

# get gene names and paths
genes <- ncbi.clean |> dplyr::distinct(gene) |> pull()
genes.files <- here::here(today.dir,glue::glue("{genes}.fasta"))

# write out
purrr::walk(genes, \(x) write_fasta(df=ncbi.clean,genez=x,dir=today.dir,pop=as.logical(opt$indiv)))

# align
purrr::walk(genes.files, \(x) align_fasta(infile=x,threads=opt$threads))

# trim
purrr::walk(genes.files, \(x) trim_fasta(infile=x,prop=opt$prop))


#### STOP IF POP ####

# end the script if in pop mode
if(opt$indiv == "true") {
   cli_report(txt="Finished running in pop mode.",rule=FALSE,alert="success")
   #cli::cli_alert_success("Finished running in pop mode.\f")
   quit(save="no")
}


##### CONCATENATE AND WRITE OUT #####
 
# read in with ape
ali.files <- here::here(today.dir,glue::glue("{genes}.aligned.trimmed.fasta"))
ali.all <- purrr::map(ali.files, \(x) ape::read.FASTA(file=x))

# convert to matrix
ali.all.mat <- purrr::map(ali.all, as.matrix)

# concatenate
genes.concat <- as.list(do.call(cbind.DNAbin,args=c(ali.all.mat,fill.with.gaps=TRUE)))

# write out
genes.concat |> ape::write.FASTA(file=here::here(today.dir,"concatenated.aligned.trimmed.fasta"))
genes.concat |> ape::write.dna(file=here::here(today.dir,"concatenated.aligned.trimmed.phy"),format="sequential",nbcol=-1,colsep="")
genes.concat |> ape::write.nexus.data(file=here::here(today.dir,"concatenated.aligned.trimmed.nex"),interleaved=FALSE)

# make partion file and write out
partition_table(mat=ali.all.mat) |> readr::write_tsv(here::here(today.dir,"concatenated.aligned.trimmed.parts"),col_names=FALSE)

# file
cli_report(txt="Trimmed alignments and partitions files written to:",rule=FALSE,alert="info")
writeLines(glue::glue("'temp/{basename(today.dir)}/concatenated.aligned.trimmed.fasta'."))
writeLines(glue::glue("'temp/{basename(today.dir)}/concatenated.aligned.trimmed.nex'."))
writeLines(glue::glue("'temp/{basename(today.dir)}/concatenated.aligned.trimmed.phy'."))
writeLines(glue::glue("'temp/{basename(today.dir)}/concatenated.aligned.trimmed.parts'."))

cli_report(txt="Alignment and trimming completed.",rule=TRUE,alert="success")
