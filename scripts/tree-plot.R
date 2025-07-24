#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Plotting phylogenetic trees ...\n")

# get args
option_list <- list( 
    make_option(c("-s","--scalefactor"), type="numeric"),
    make_option(c("-h","--hratio"), type="numeric"),
    make_option(c("-w","--width"), type="numeric"),
    make_option(c("-c","--colour"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt <- NULL
#opt$scalefactor <- 1
#opt$hratio <- 0.5
#opt$width <- 0.6
#opt$colour <- "species"

##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))
ncbi.clean <- read_csv(here(today.dir,"ncbi-clean.csv"),show_col_types=FALSE,col_types=cols(.default=col_character()))

# list trees
tree.files <- list.files(today.dir,pattern="*.fasta.raxml.bestTree$",recursive=FALSE,full.names=TRUE)


##### AUTO PLOT #####

# make tip labels for family, genus, species
if(opt$colour == "family") {
    ncbi.clean.tips <- ncbi.clean |> distinct(scientificName,genus,family) |> mutate(tiplabel=glue("{family} | {str_replace_all(scientificName,'_',' ')}"))  |> mutate(tip.colour=family)
}
if(opt$colour == "genus") {
    ncbi.clean.tips <- ncbi.clean |> distinct(scientificName,genus,family) |> mutate(tiplabel=glue("{family} | {str_replace_all(scientificName,'_',' ')}"))  |> mutate(tip.colour=genus)
}
if(opt$colour == "species") {
    ncbi.clean.tips <- ncbi.clean |> distinct(gbAccession,scientificName,genus,family) |> mutate(tiplabel=glue("{gbAccession} | {str_replace_all(scientificName,'_',' ')}")) |> mutate(tip.colour=scientificName)
}

# fun plot function over all trees
purrr::walk(tree.files, \(x) ggtree_autoplot(path=x,tb=ncbi.clean.tips,scale.factor=opt$scalefactor,width=opt$width,hratio=opt$hratio))

# print
writeLines(glue("\nTrees written to:\n",.trim=FALSE))
writeLines(glue("{basename(tree.files)}.pdf\n",.trim=FALSE))
