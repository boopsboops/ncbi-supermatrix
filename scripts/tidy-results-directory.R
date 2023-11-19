#!/usr/bin/env Rscript

##### LOAD LIBS FUNS ARGS #####

source(here::here("scripts/load-libs-funs.R"))

# info
writeLines("Cleaning and tidying directory ...\n")


##### LOAD DATA #####

# get latest dir
today.dir <- sort(list.dirs(here("temp"),recursive=FALSE),decreasing=TRUE)[1]
writeLines(glue("Working in directory 'temp/{basename(today.dir)}'.\n",.trim=FALSE))

# list files
all.files <- list.files(today.dir,recursive=FALSE,full.names=TRUE)
#print(fasta.files)


##### MAKE OUTDIRS #####

if(!dir.exists(here(today.dir,"trees"))) {dir.create(here(today.dir,"trees"),recursive=TRUE)}
if(!dir.exists(here(today.dir,"alignments"))) {dir.create(here(today.dir,"alignments"),recursive=TRUE)}
if(!dir.exists(here(today.dir,"pdf"))) {dir.create(here(today.dir,"pdf"),recursive=TRUE)}
if(!dir.exists(here(today.dir,"metadata"))) {dir.create(here(today.dir,"metadata"),recursive=TRUE)}
if(!dir.exists(here(today.dir,"backup"))) {dir.create(here(today.dir,"backup"),recursive=TRUE)}


##### CLEAN UP #####

# copy
invisible(file.copy(str_subset(all.files,"raxml\\.bestTree$"),here(today.dir,"trees"),copy.date=TRUE))
invisible(file.copy(str_subset(all.files,"\\.aligned\\.trimmed\\.fasta$|\\.phy$|\\.nex$|\\.parts$"),here(today.dir,"alignments"),copy.date=TRUE))
invisible(file.copy(str_subset(all.files,"raxml\\.bestTree.pdf$"),here(today.dir,"pdf"),copy.date=TRUE))
invisible(file.copy(str_subset(all.files,"\\.csv$"),here(today.dir,"metadata"),copy.date=TRUE))
invisible(file.copy(all.files,here(today.dir,"backup"),copy.date=TRUE))

# delete
invisible(file.remove(all.files))


##### REPORT #####

writeLines(glue("All Newick tree files have been moved to 'temp/{basename(today.dir)}/trees'.\n",.trim=FALSE))
writeLines(glue("All PDF tree files have been moved to 'temp/{basename(today.dir)}/pdf'.\n",.trim=FALSE))
writeLines(glue("All final alignment files have been moved to 'temp/{basename(today.dir)}/alignments'.\n",.trim=FALSE))
writeLines(glue("All metadata files have been moved to 'temp/{basename(today.dir)}/metadata'.\n",.trim=FALSE))
writeLines(glue("A backup of all intermediate and final files have been moved to 'temp/{basename(today.dir)}/backup'.\n",.trim=FALSE))
