#!/usr/bin/env Rscript

##### EXAMPLE #####

# scripts/annotate-ncbi.R -t 1 -c fishbase


##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# get latest dir
today.dir <- sort(grep("/Results_",list.dirs(here::here("temp"),recursive=FALSE),value=TRUE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
#today.dir <- here("temp",paste0("Results_",Sys.Date()))
#if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}

# info
writeLines("Annotating sequence data ...")

# get args
option_list <- list(
    make_option(c("-c","--classification"), type="character"),
    make_option(c("-t","--threads"), type="numeric")
)

# set args
#opt <- NULL
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list,add_help_option=FALSE))
#cores <- opt$threads
#opt <- NULL
#opt$threads <- 1
#opt$classification <- "fishbase"


##### RETRIEVE METADATA FROM NCBI #####

# load up fasta
all.fas <- here::here(today.dir,list.files(today.dir,"gene\\."))

# tabulate genes
ids.all.tab <- dplyr::bind_rows(mapply(tabulate_genes,all.fas,SIMPLIFY=FALSE,USE.NAMES=FALSE))
ncbi.accs <- dplyr::pull(ids.all.tab,acc_no)

# chunk 200 should result in string of around 2200 chars
if(length(ncbi.accs) < 200) {
    chunk <- length(ncbi.accs) } else {
    chunk <- 200
}

# randomise accessions
set.seed(42)
ids.all <- sample(ncbi.accs)
# chunk
chunk.frag <- unname(split(ids.all, ceiling(seq_along(ids.all)/chunk)))
writeLines("\nRetrieving metadata from NCBI ...\n")

#threads error
err.msg <- glue::glue("Error! You requested {opt$threads} thread(s), but there are only {length(chunk.frag)} chunk(s). Please use equal or fewer threads than chunks.")
if(opt$threads > length(chunk.frag)) {stop(err.msg)}

# parallel ncbi
ncbi.frag <- mcmapply(FUN=ncbi_byid_parallel, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=opt$threads)

# check for errors (should all be "data.frame")
if(length(sapply(ncbi.frag,class)) == length(which(sapply(ncbi.frag,class) == "data.frame"))) {
    writeLines("\nNCBI metadata sucessfully retrieved.")
    } else {writeLines("\nNCBI search failed, try again")}

# join
frag.df <- tibble::as_tibble(bind_rows(ncbi.frag))


##### ANNOTATE WITH FISHBASE OR NCBI TAXONOMY #####

# info
writeLines("\nNow retrieving taxonomy ...\n")

# get genera
frag.df.gens <- frag.df |> dplyr::mutate(genus=stringr::str_split_fixed(taxon," ",2)[,1])
list.gens <- frag.df.gens |> dplyr::distinct(genus) |> dplyr::pull(genus)

# run if else for classification 
if (opt$classification=="ncbi") {
    # NCBI #
    writeLines("\nUsing NCBI taxonomy.\n")
    fb.gens <- ncbi_annotate(genera=list.gens) 
    } else if (opt$classification=="fishbase") {
    # FISHBASE #
    writeLines(glue::glue("Using FishBase taxonomy version {sort(rfishbase::available_releases(),decreasing=TRUE)[1]}."))
    # get fb table
    fb.gens <- fishbase_annotate(genera=list.gens)
    # close connections
    #rfishbase::db_disconnect()
    # stop
} else stop(writeLines("Error! The '-c' flag must be either 'ncbi' or 'fishbase'."))

# join to ncbi table and to genes table
frag.df.gens.genes <- frag.df.gens |> dplyr::left_join(fb.gens,by=join_by(genus)) |> dplyr::left_join(ids.all.tab,by=join_by(acc_no))

# write out
frag.df.gens.genes |> readr::write_csv(here::here(today.dir,"ncbi-raw.csv"))

# info
writeLines(glue::glue("\nAnnotation of sequence data completed. Output written to 'temp/{basename(today.dir)}/ncbi-raw.csv'.\n",.trim=FALSE))
