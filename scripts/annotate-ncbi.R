#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("\nAnnotating sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# make date dir for results
today.dir <- here("temp",paste0("Results_",Sys.Date()))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}


# load up fasta
all.fas <- here(today.dir,list.files(today.dir,"gene\\."))

# tabulate genes
ids.all.tab <- bind_rows(mapply(tabulate_genes,all.fas,SIMPLIFY=FALSE,USE.NAMES=FALSE))
ncbi.accs <- pull(ids.all.tab,acc)

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
# parallel ncbi
if(opt$threads > length(chunk.frag)) {stop("Error! You requested more threads than chunks. Please use fewer threads.")}
ncbi.frag <- mcmapply(FUN=ncbi_byid_parallel, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=1)#opt$threads

# check for errors (should all be "data.frame")
if(length(sapply(ncbi.frag,class)) == length(which(sapply(ncbi.frag,class) == "data.frame"))) {
    writeLines("\nNCBI metadata sucessfully retrieved")
    } else {writeLines("\nNCBI search failed, try again")}

# join
frag.df <- as_tibble(bind_rows(ncbi.frag))

# filter 
frag.df.clean <- clean_ncbi(df=frag.df)


# join with genes 
frag.df.clean.genes <- frag.df.clean %>% 
    left_join(rename(ids.all.tab,gbAccession=acc),by=join_by(gbAccession)) %>% 
    filter(!is.na(gene)) %>%
    select(scientificName,label,dbid,gbAccession,gene,length,organelle,catalogNumber,country,publishedAs,publishedIn,publishedBy,date,decimalLatitude,decimalLongitude,notesGenBank,taxonomy,nucleotides)

### WRITE OUT


# filter by species and sequence length
frag.df.clean.genes.singles <- frag.df.clean.genes %>% group_by(gene,scientificName) %>%
    slice_max(order_by=length,with_ties=FALSE,n=1) %>%
    ungroup()



# print and write out table
#writeLines("\nAll clusters with number sequences, file size in kb, and description of first sequence.\n")
#writeLines(paste("\nTable written out to",here(today.dir,"clusters.csv"),"\n"))
#all.clusters.desc %>% select(-path) %>% print(n=Inf,width=Inf)
