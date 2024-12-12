#!/usr/bin/env Rscript

# LOAD PACKAGES
suppressMessages({
    library("here")
    library("rentrez")
    library("ape")
    library("tidyverse")
    library("optparse")
    library("crul")
    library("traits")
    library("parallel")
    library("glue")
    library("rfishbase")
    library("taxize")
    library("knitr")
    library("castor")
    library("treeio")
    library("ggtree")
    library("phangorn")
    library("cli")
})


# LOAD FROM GITHUB
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/main/RScripts/tab2fas.R")


# FUN TO GET DESCRIPTIONS FROM GENBANK
read_genbank <- function(accs) {
    Sys.sleep(1)
    gb.accs <- ape::read.GenBank(accs)
    gb.accs.tib <- tibble(facc=labels(gb.accs),description=attr(gb.accs,"description"))
    writeLines(glue("Cluster information downloaded for {length(accs)} hits."))
    return(gb.accs.tib)
}


# ENTREZ PARALLEL SEARCH
entrez_fetch_parallel <- function(webhist,chunks,file) {
    chunk.first <- dplyr::first(chunks)
    chunk.last <- dplyr::last(chunks)
    chunk.len <- length(chunks)
    writeLines(paste0("Downloading hits ",chunk.first," to ",chunk.last," from NCBI ..."))
    fetch.fas <- rentrez::entrez_fetch(db="nucleotide",rettype="fasta",web_history=webhist,retstart=chunk.first,retmax=chunk.len)
    write(fetch.fas,file=file,append=TRUE)
    Sys.sleep(3)
    #invisible()
}


# ENTREZ DOWNLOAD
entrez_download <- function(clade,minlen,maxlen,batchsize,fasout,append,dry) {
    if (append=="false") {
        file.create(fasout,overwrite=TRUE)
    } else if (append=="true") {
        Sys.sleep(1)
    } else stop("\nError! the append flag '-a' must be 'true' or 'false'.\n")
    rtm <- as.integer(batchsize)
    nots <- glue("NOT \"PREDICTED\"[All Fields] NOT \"UNVERIFIED\"[All Fields] NOT \"microsatellite\"[All Fields] NOT \"pseudogene\"[All Fields] NOT \"genomic survey sequence\"[All Fields]")# NOT \"transcribed\"[All Fields] NOT \"cDNA\"[All Fields]  NOT \"mRNA\"[All Fields] NOT \"whole genome shotgun sequence\"[All Fields]
    if (rtm > 9999) {stop("\nError! Batch size can be no larger than 9,999.\n")}
    es.res <- rentrez::entrez_search(db="nucleotide",term=glue("(\"{clade}\"[ORGANISM] AND {minlen}:{maxlen}[SLEN]) {nots}"),retmax=999999,use_history=TRUE)
    if (dry=="true") {stop(glue("\rInformation: The number of hits for this search is {es.res$count}. Ensure that your batch size '-b' is lower than this value (but no larger than 9,999).\n",.trim=FALSE),call.=FALSE)}
    if (rtm > es.res$count) {stop(glue("\nError! Batch size value ({rtm}) is larger than number of hits ({es.res$count}). Please use a smaller '-b' value.\n",.trim=FALSE))}
    writeLines(paste("\nNCBI search reports",length(es.res$ids),"hits for",clade,"...\n"))
    if (length(es.res$ids) > 999999) {stop("\nError! You have more hits than can be downloaded. Please break up the search into seperate clades.\n")}
    es.post <- entrez_post(db="nucleotide",web_history=es.res$web_history)
    query.split <- split(0:es.res$count,cut(0:es.res$count,ceiling(es.res$count/rtm),labels=FALSE))
    x <- mapply(function(x) entrez_fetch_parallel(webhist=es.post,chunks=x,file=fasout), x=query.split,SIMPLIFY=TRUE,USE.NAMES=FALSE)
    es.fas <- ape::read.FASTA(fasout)
    #if (length(es.fas) != es.res$count) {stop("Error! Not all sequences were downloaded. Please try again.")}
    writeLines(paste("\nDone!",length(es.fas),"sequences downloaded in FASTA format.\n"))
    return(es.fas)
}


# FUNCTION TO RUN PARALLEL NCBI_BYID WITH TIMEOUT AND REPEAT
ncbi_byid_parallel <- function(accs){
    start_time <- Sys.time()
    Sys.sleep(time=runif(n=1,min=0,max=3))
    crul::set_opts(http_version=2)
    ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    if(class(ncbi.tab)!="data.frame") {
        #writeLines("Error found! Repeating ...")
        Sys.sleep(time=3)
        crul::set_opts(http_version=2)
        ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    } else {
        ncbi.tab <- ncbi.tab
    }
    if(class(ncbi.tab)!="data.frame") {
        stop(writeLines("Searches failed ... aborted")) 
    } else {
        end_time <- Sys.time()
        writeLines(paste0("Metadata for ",length(accs)," accessions downloaded (starting ",accs[1],"). Download took ",round(as.numeric(end_time-start_time),digits=2)," seconds."))
        return(ncbi.tab)
    }
}


# DEREPLICATE FASTA
dereplicate_fasta <- function(infile,dereplicate) {
    outfile <- str_replace_all(infile,"\\.fasta",".derep.fasta")
        if (dereplicate=="true") {
            exe.string <- paste("vsearch --threads 0 --fastx_uniques",infile,"--fastaout",outfile)
            system(exe.string)
        } else if (dereplicate=="false") {
            file.copy(infile,outfile)
        } else stop("\nError! the dereplicate flag '-d' must be 'true' or 'false'.\n")
}


# FILTER FASTA
filter_fasta <- function(infile,maxns) {
    infile <- str_replace_all(infile,"\\.fasta",".derep.fasta")
    outfile <- str_replace_all(infile,"\\.fasta",".filtered.fasta")
    exe.string <- paste("vsearch --threads 0 -fastx_filter",infile,"--fastq_maxns",maxns,"--fastaout",outfile)
    system(exe.string)
}


# CLUSTER FASTA
cluster_fasta <- function(infile,identity) {
    infile <- str_replace_all(infile,"\\.fasta",".derep.filtered.fasta")
    outfile <- here(today.dir,"cluster.")
    str_replace_all(infile,"genbank-dump.derep.filtered.fasta","cluster.")
    exe.string <- paste("vsearch --threads 0 --cluster_fast",infile,"--id",identity,"--clusters",outfile)
    system(exe.string)
}


# LOAD AND RENAME FASTA
load_and_rename <- function(gene,path,wd) {
    fas <- ape::read.FASTA(path)
    path.gene <- here(today.dir,paste("gene",gene,"fasta",sep="."))
    ape::write.FASTA(fas,path.gene)
}


# TABULATE GENES
tabulate_genes <- function(path) {
    fas <- ape::read.FASTA(path)
    accs <- names(fas)
    gene <- str_split_fixed(basename(path),"\\.",3)[,2]
    accs.gene <- tibble(acc_no=accs,gene=gene)
    return(accs.gene)
}


# CLEAN NCBI DATA
clean_ncbi <- function(df) {
    df.clean <- df |>
        distinct(gi_no,.keep_all=TRUE) |> 
        # filter
        filter(gi_no!="NCBI_GENOMES") |> 
        filter(!is.na(sequence)) |> 
        filter(!grepl("UNVERIFIED:",gene_desc)) |>
        filter(!grepl("PREDICTED:",gene_desc)) |>
        filter(!grepl("similar to",gene_desc)) |> 
        filter(!grepl("mRNA",gene_desc)) |> 
        filter(!grepl("cDNA",gene_desc)) |> 
        filter(!grepl("transcribed",gene_desc)) |> 
        filter(!grepl("-like",gene_desc)) |>
        # fix lan lon
        mutate(lat=paste(str_split_fixed(lat_lon, " ", 4)[,1], str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(str_split_fixed(lat_lon, " ", 4)[,3], str_split_fixed(lat_lon, " ", 4)[,4])) |>
        mutate(lat=if_else(grepl(" N",lat), true=str_replace_all(lat," N",""), false=if_else(grepl(" S",lat), true=paste0("-",str_replace_all(lat," S","")), false=lat))) |>
        mutate(lon=if_else(grepl(" E",lon), true=str_replace_all(lon," E",""), false=if_else(grepl(" W",lon), true=paste0("-",str_replace_all(lon," W","")), false=lon))) |> 
        mutate(lat=str_replace_all(lat,"^ ", NA_character_), lon=str_replace_all(lon,"^ ", NA_character_)) |>
        mutate(lat=suppressWarnings(as.numeric(lat)), lon=suppressWarnings(as.numeric(lon))) |> 
        # tidy up
        mutate(length=as.numeric(length)) |>
        mutate(gi_no=as.character(gi_no)) |>
        mutate(scientificNameGenBank=taxon) |>
        rename(scientificName=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence) |>
        select(scientificName,scientificNameGenBank,class,order,family,genus,dbid,gbAccession,gene,length,organelle,catalogNumber,country,publishedAs,publishedIn,publishedBy,date,decimalLatitude,decimalLongitude,nucleotides)
    return(df.clean)
}


# SPLIT NAMES FUNCTION
splitter <- function(x,n) {
    splitz <- stringr::str_split_fixed(x,"_",5)
    splitz1 <- splitz[,1]
    splitz2 <- splitz[,2]
    splitz3 <- splitz[,3]
    splitz4 <- splitz[,4]
    splitz5 <- splitz[,5]
    gstring <- case_when(
        n==2 ~ glue("{splitz1}_{splitz2}"),
        n==3 ~ glue("{splitz1}_{splitz2}_{splitz3}"),
        n==4 ~ glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}"),
        n==5 ~ glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}_{splitz5}")
        )
    return(gstring)
}


# CLEAN NCBI DATA
clean_names <- function(df,n) {
    df.clean <- df |>
        mutate(label=scientificName) |>
        mutate(label=str_replace_all(label," ","_")) |>
        mutate(label=str_replace_all(label,"'","")) |>
        mutate(label=str_replace_all(label,",","")) |>
        mutate(label=str_replace_all(label,"\\(","")) |>
        mutate(label=str_replace_all(label,"\\)","")) |>
        mutate(label=str_replace_all(label,":","-")) |>
        mutate(label=str_replace_all(label,"cf\\._","cf.")) |>
        mutate(label=str_replace_all(label,"aff\\._","aff.")) |>
        mutate(label=str_replace_all(label,"sp\\._","sp.")) |>
        mutate(elements=str_count(label,"_")+1) |>
        mutate(label=case_when(
            n == 2 ~ splitter(label,n=2),
            n == 3 ~ splitter(label,n=3),
            n == 4 ~ splitter(label,n=4),
            n == 5 ~ splitter(label,n=5),
            #str_detect(label,"cf\\.") ~ splitter(label,n=3),
            #elements > 4 ~ splitter(label,n=5)
            )) |>
        mutate(label=as.character(label)) |>
        mutate(label=str_replace_all(label,"_+$","")) |>
        mutate(scientificName=label) |>
        select(!c(label,elements))
    return(df.clean)
}


# ANNOTATE WITH FISHBASE TAXONOMY
fishbase_annotate <- function(genera) {
    genera.taxonomy <- rfishbase::fb_tbl("classification_tree",server="fishbase",version="latest") |> 
        distinct(Class,Order,Family,Genus) |>
        rename(class=Class,order=Order,family=Family,genus=Genus)
    genera.taxonomy <- genera.taxonomy |> filter(genus %in% genera)
    return(genera.taxonomy)
}


# ANNOTATE WITH NCBI TAXONOMY
ncbi_annotate <- function(genera) {
    tax <- c("class","order","family")
    tax.tab <- taxize::tax_name(sci=genera,get=tax,db="ncbi")
    tax.tab.tib <- tax.tab |> as_tibble() |> select(-db) |> rename(genus=query) |> relocate(genus,.after="family")
    return(tax.tab.tib)
}


# WRITE A FASTA FILE FROM DF
write_fasta <- function(df,genez,dir,pop) {
    if(isFALSE(pop)) {
        df.sub <- df |> filter(gene == genez) |> tab2fas(seqcol="nucleotides",namecol="scientificName")
    } else if(isTRUE(pop)) {
        df.sub <- df |> filter(gene == genez) |> tab2fas(seqcol="nucleotides",namecol="gbAccession")
    }
    df.sub |> write.FASTA(file=here(dir,glue("{genez}.fasta")))
    #writeLines(glue("Unaligned fasta file written to 'temp/{basename(dir)}/{genez}.fasta'.",.trim=FALSE))
}


# ALIGN FASTA
align_fasta <- function(infile,threads) {
    outfile <- str_replace_all(infile,"\\.fasta",".aligned.fasta")
    exe.string <- glue("mafft --quiet --auto --thread {threads} {infile} > {outfile}")
    system(exe.string)
    #writeLines(glue("Aligned fasta file written to '{basename(outfile)}'.",.trim=FALSE))
}


# TRIM FASTA
trim_fasta <- function(infile,prop) {
    infile <- str_replace_all(infile,"\\.fasta",".aligned.fasta")
    outfile <- str_replace_all(infile,"\\.fasta",".trimmed.fasta")
    exe.string <- glue("trimal -in {infile} -out {outfile} -gt {prop}")
    system(exe.string)
    writeLines(glue("\nAligned and trimmed fasta file written to '{basename(outfile)}'.\n",.trim=FALSE))
}


# MAKE PARTS
partition_table <- function(mat) {
    dims <- purrr::map(mat, \(x) dim(x)[2]) |> unlist()
    parts.tib <- tibble(gene=genes,len=dims,cs=cumsum(dims)) |> 
        mutate(cs1=cs+1,start=lag(cs1)) |>
        replace_na(list(start=1)) |>
        rename(end=cs) %>%
        select(gene,start,end) %>%
        mutate(parts=as.character(glue("DNA, {gene} = {start}-{end}"))) %>%
        select(parts) 
    return(parts.tib)
}


# NEW RAXML-NG FUN
raxml_ng <- function(file,model,maxthreads,epsilon,verbose) {
    string.search <- glue("raxml-ng --search -threads auto{{{maxthreads}}} --workers auto --tree pars{{1}} --redo --seed 42 --model {model} --lh-epsilon {epsilon} --msa {file}")
    if (verbose == "true") {
        system(command=string.search,ignore.stdout=FALSE)
        } else if (verbose == "false") {
        system(command=string.search,ignore.stdout=TRUE)
        } else stop(writeLines("Error! the '-v' flag must be 'true' or 'false'."))
    rax.tr <- ape::read.tree(file=glue("{file}.raxml.bestTree"))
    writeLines(glue("\nTree search completed for '{basename(file)}'.",.trim=FALSE))
    return(rax.tr)
}


# MAKE TREE PLOTTING FUN
ggtree_autoplot <- function(path,tb,scale.factor,width,hratio) {
    tr <- treeio::read.newick(path)
    tr <- tr |> castor::root_in_edge(root_edge=which.max(tr$edge.length))
    #tr <- phangorn::midpoint(tr)
    tree.length <- max(castor::get_all_distances_to_root(tr)) + (max(castor::get_all_distances_to_root(tr)) * width)
        p <- ggtree(tr, ladderize=TRUE,right=TRUE,size=0.7) %<+% tb
        pp <- p + geom_tiplab(offset=0,aes(label=tiplabel),align=FALSE,size=4) +
        geom_tippoint(aes(color=tip.colour),size=2.5) +
        theme(legend.position="none") +
        xlim(0,tree.length)
    ggsave(filename=glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=200,height=200*hratio,units="mm",scale=scale.factor)
}


# REPORT
writeLines("\nPackages and functions loaded.\n")
