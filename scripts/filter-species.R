#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Filtering sequence data ...\n")

# get args
option_list <- list( 
    make_option(c("-n","--names"), type="numeric"),
    make_option(c("-i","--indiv"), type="character"),
    make_option(c("-o","--outgroup"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$names <- 2
#opt$indiv <- "true"
#opt$outgroup <- "true"

##### LOAD DATA #####

# get latest dir
today.dir <- sort(grep("/Results_",list.dirs(here::here("temp"),recursive=FALSE),value=TRUE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.raw <- read_csv(here(today.dir,"ncbi-raw.csv"),show_col_types=FALSE,col_types=cols(.default=col_character()))
excl <- read_csv(here("assets/exclusions.csv"),show_col_types=FALSE,col_types=cols(.default=col_character()))

# og
if(file.exists(here(today.dir,"outgroup.txt"))) {
    og <- readLines(here(today.dir,"outgroup.txt"), warn=FALSE)
}

##### CLEAN AND REMOVE UNWANTED SEQUENCES #####

# clean up the table
ncbi.clean <- ncbi.raw |> clean_ncbi() |> clean_names(n=opt$names)
#ncbi.clean |> distinct(scientificName) |> print(n=Inf)

# drop exclusions
ncbi.clean.excl <- ncbi.clean |> filter(!gbAccession %in% pull(excl,acc))

# count drops and report
ndrop <- nrow(ncbi.clean) - nrow(ncbi.clean.excl)
#writeLines(glue("\nA total of {ndrop} sequence(s) have been excluded after filtering with 'assets/exclusions.csv'.\n",.trim=FALSE))


##### FILTER FOR LONGEST SINGLE OUTGROUP #####

if(opt$outgroup == "true") {
    # filter by outgroup
    cli::cli_alert_info("Filtering outgroup to one indiv.")
    ncbi.clean.excl <- bind_rows(
        ncbi.clean.excl |> filter(scientificName == og) |> slice_max(length,n=1,with_ties=FALSE),
        ncbi.clean.excl |> filter(scientificName != og)
    )
    } else if(opt$outgroup=="false") {
    # or do nothing
    ncbi.clean.excl <- ncbi.clean.excl
    } else {stop(writeLines("Error! the '-o' flag must be 'true' or 'false'."))
}


##### FILTER FOR LONGEST ONE INDIV PER SP #####

if(opt$indiv == "false") {
    # filter by species and sequence length
    cli::cli_alert_info("Running in species mode.")
    ncbi.clean.excl.filt <- ncbi.clean.excl |> 
        group_by(gene,scientificName) |>
        slice_max(order_by=length,with_ties=FALSE,n=1) |>
        ungroup() |> 
        arrange(scientificName,gene)
    } else if(opt$indiv == "true") {
    # or do nothing if pop
    cli::cli_alert_info("Running in pop mode.\f")
    ncbi.clean.excl.filt <- ncbi.clean.excl
    } else {stop(writeLines("Error! the '-i' flag must be 'true' or 'false'."))
}


##### WRITE OUT #####

# write out
ncbi.clean.excl.filt |> write_csv(here(today.dir,"ncbi-clean.csv"))

# get numbers
nseq <- nrow(ncbi.clean.excl.filt)
nspp <- distinct(ncbi.clean.excl.filt,scientificName) |> nrow()
ngene <- distinct(ncbi.clean.excl.filt,gene) |> nrow()

# print
writeLines(glue("After filtering, {nseq} sequence(s) for {nspp} species have been retained for {ngene} genes as follows:",.trim=FALSE))
ncbi.clean.excl.filt |> group_by(gene) |> summarise(n=n()) |> knitr::kable()

# file
writeLines(glue("\nCleaning of sequence data completed. Output written to 'temp/{basename(today.dir)}/ncbi-clean.csv'.\n",.trim=FALSE))
