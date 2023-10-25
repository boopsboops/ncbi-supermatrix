#!/usr/bin/env Rscript

##### EXAMPLE #####

# scripts/annotate-ncbi.R -t 1 -c fishbase


##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
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
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
cores <- opt$threads
#opt$threads <- 1
#opt$classification <- "fishbase"


##### RETRIEVE METADATA FROM NCBI #####

# load up fasta
all.fas <- here(today.dir,list.files(today.dir,"gene\\."))

# tabulate genes
ids.all.tab <- bind_rows(mapply(tabulate_genes,all.fas,SIMPLIFY=FALSE,USE.NAMES=FALSE))
ncbi.accs <- pull(ids.all.tab,acc_no)

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
err.msg <- glue("Error! You requested {opt$threads} thread(s), but there are only {length(chunk.frag)} chunk(s). Please use equal or fewer threads than chunks.")
if(opt$threads > length(chunk.frag)) {stop(err.msg)}

# parallel ncbi
ncbi.frag <- mcmapply(FUN=ncbi_byid_parallel, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=opt$threads)

# check for errors (should all be "data.frame")
if(length(sapply(ncbi.frag,class)) == length(which(sapply(ncbi.frag,class) == "data.frame"))) {
    writeLines("\nNCBI metadata sucessfully retrieved.")
    } else {writeLines("\nNCBI search failed, try again")}

# join
frag.df <- as_tibble(bind_rows(ncbi.frag))


##### ANNOTATE WITH FISHBASE OR NCBI TAXONOMY #####

# info
writeLines("\nNow retrieving taxonomy ...\n")

# get genera
frag.df.gens <- frag.df |> mutate(genus=str_split_fixed(taxon," ",2)[,1])
list.gens <- frag.df.gens |> distinct(genus) |> pull(genus)

# run if else for classification 
if (opt$classification=="ncbi") {
    # NCBI #
    writeLines("\nUsing NCBI taxonomy.\n")
    fb.gens <- ncbi_annotate(genera=list.gens) 
    } else if (opt$classification=="fishbase") {
    # FISHBASE #
    writeLines(glue("Using FishBase taxonomy version {rfishbase::available_releases()[1]}."))
    # get fb table
    fb.gens <- fishbase_annotate(genera=list.gens)
    # close connections
    rfishbase::db_disconnect()
    # stop
} else stop(writeLines("Error! The '-c' flag must be either 'ncbi' or 'fishbase'."))

# join to ncbi table and to genes table
frag.df.gens.genes <- frag.df.gens |> left_join(fb.gens,by=join_by(genus)) |> left_join(ids.all.tab,by=join_by(acc_no))

# write out
frag.df.gens.genes |> write_csv(here(today.dir,"ncbi-raw.csv"))

# info
writeLines(glue("\nAnnotation of sequence data completed. Output written to 'temp/{basename(today.dir)}/ncbi-raw.csv'.\n",.trim=FALSE))
