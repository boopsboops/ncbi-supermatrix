#!/usr/bin/env Rscript

GOT UP TO HERE


##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("\nFiltering sequence data ...\n")

# get args
option_list <- list(
    make_option(c("-f","--fishbase"), type="character"),
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
cores <- opt$threads
#cores <- 1


##### FILTER #####




# filter 
frag.df.clean <- clean_ncbi(df=frag.df)


# join with genes 
frag.df.clean.genes <- frag.df.clean |> 
    left_join(rename(ids.all.tab,gbAccession=acc),by=join_by(gbAccession)) |> 
    filter(!is.na(gene)) |>
    select(scientificName,label,dbid,gbAccession,gene,length,organelle,catalogNumber,country,publishedAs,publishedIn,publishedBy,date,decimalLatitude,decimalLongitude,notesGenBank,taxonomy,nucleotides)

### WRITE OUT


# filter by species and sequence length
frag.df.clean.genes.singles <- frag.df.clean.genes |> group_by(gene,scientificName) |>
    slice_max(order_by=length,with_ties=FALSE,n=1) |>
    ungroup()



# print and write out table
#writeLines("\nAll clusters with number sequences, file size in kb, and description of first sequence.\n")
#writeLines(paste("\nTable written out to",here(today.dir,"clusters.csv"),"\n"))
#all.clusters.desc |> select(-path) |> print(n=Inf,width=Inf)
