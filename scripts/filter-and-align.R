#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Filtering sequence data ...\n")


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.raw <- read_csv(here(today.dir,"ncbi-raw.csv"),show_col_types=FALSE)
excl <- read_csv(here("assets/exclusions.csv"),show_col_types=FALSE)


##### CLEAN AND REMOVE UNWANTED SEQUENCES #####

# clean up the table
ncbi.clean <- ncbi.raw |> clean_ncbi() |> clean_names()

# drop exclusions
ncbi.clean.excl <- ncbi.clean |> filter(!gbAccession %in% pull(excl,acc))

# count drops and report
ndrop <- nrow(ncbi.clean) - nrow(ncbi.clean.excl)
#writeLines(glue("\nA total of {ndrop} sequence(s) have been excluded after filtering with 'assets/exclusions.csv'.\n",.trim=FALSE))


##### FILTER FOR LONGEST ONE INDIV PER SP #####

# filter by species and sequence length
ncbi.clean.excl.filt <- ncbi.clean.excl |> 
    group_by(gene,scientificName) |>
    slice_max(order_by=length,with_ties=FALSE,n=1) |>
    ungroup() |> 
    arrange(scientificName,gene)


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
