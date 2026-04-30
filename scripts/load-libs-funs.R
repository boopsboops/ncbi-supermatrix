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
    gb.accs.tib <- tibble::tibble(facc=labels(gb.accs),description=attr(gb.accs,"description"))
    #writeLines(glue("Cluster information downloaded for {length(accs)} hits."))
    cli_report(txt=glue::glue("Cluster information downloaded for {length(accs)} hits."),rule=FALSE,alert="info")
    return(gb.accs.tib)
}


# ENTREZ PARALLEL SEARCH
entrez_fetch_parallel <- function(webhist,chunks,file) {
    chunk.first <- dplyr::first(chunks)
    chunk.last <- dplyr::last(chunks)
    chunk.len <- length(chunks)
    cli_report(txt=glue::glue("Downloading hits {chunk.first} to {chunk.last} from NCBI ..."),rule=FALSE,alert="info")
    #writeLines(paste0("Downloading hits ",chunk.first," to ",chunk.last," from NCBI ..."))
    fetch.fas <- rentrez::entrez_fetch(db="nucleotide",rettype="fasta",web_history=webhist,retstart=chunk.first,retmax=chunk.len)
    write(fetch.fas,file=file,append=TRUE)
    Sys.sleep(3)
    #invisible()
}


# ENTREZ DOWNLOAD
entrez_download <- function(clade,minlen,maxlen,batchsize,fasout,append,dry,mito) {
    if (append=="false") {
        file.create(fasout,overwrite=TRUE)
    } else if (append=="true") {
        Sys.sleep(1)
    } else stop(cli::cli_alert_danger("ERROR The append flag '-a' must be 'true' or 'false'."))
    rtm <- as.integer(batchsize)
    nots <- glue::glue("NOT \"PREDICTED\"[All Fields] NOT \"UNVERIFIED\"[All Fields] NOT \"microsatellite\"[All Fields] NOT \"pseudogene\"[All Fields] NOT \"genomic survey sequence\"[All Fields]")# NOT \"transcribed\"[All Fields] NOT \"cDNA\"[All Fields]  NOT \"mRNA\"[All Fields] NOT \"whole genome shotgun sequence\"[All Fields]
    if (rtm > 9999) {stop(cli::cli_alert_danger("ERROR Batch size can be no larger than 9,999."))}
    if (mito=="true") {
        es.res <- rentrez::entrez_search(db="nucleotide",term=glue::glue("(\"{clade}\"[ORGANISM] AND {minlen}:{maxlen}[SLEN]) AND (mitochondrion[TITLE] OR mitochondrial[TITLE]) {nots}"),retmax=999999,use_history=TRUE)
    }
    if (mito=="false") {
        es.res <- rentrez::entrez_search(db="nucleotide",term=glue::glue("(\"{clade}\"[ORGANISM] AND {minlen}:{maxlen}[SLEN]) {nots}"),retmax=999999,use_history=TRUE)
    }
    if (dry=="true") {stop(cli_report(txt=glue::glue("Information: The number of hits for this search is {es.res$count}. Ensure that your batch size '-b' is lower than this value (but no larger than 9,999)."),rule=FALSE,alert="info"),call.=FALSE)}
    if (rtm > es.res$count) {stop(cli::cli_alert_danger(glue::glue("Error! Batch size value ({rtm}) is larger than number of hits ({es.res$count}). Please use a smaller '-b' value.")))}
    cli_report(txt=glue::glue("NCBI search reports {length(es.res$ids)} hits for {clade}."),rule=FALSE,alert="info")
    #writeLines(paste("\nNCBI search reports",length(es.res$ids),"hits for",clade,"...\n"))
    if (es.res$count==0) {stop(cli::cli_alert_danger(glue::glue("Information: The number of hits for this search is {es.res$count}. Nothing to download.")),call.=FALSE)}
    if (length(es.res$ids) > 999999) {stop(cli::cli_alert_danger("ERROR You have more hits than can be downloaded. Please break up the search into seperate clades."))}
    es.post <- rentrez::entrez_post(db="nucleotide",web_history=es.res$web_history)
    if (rtm==0 & es.res$count >= 100) {
        query.split <- split(0:es.res$count,cut(0:es.res$count,ceiling(es.res$count/(es.res$count-1)),labels=FALSE))
    }
    if (rtm==0 & es.res$count > 1 & es.res$count <100) {
        query.split <- list(0:es.res$count)
    }
    if (rtm==0 & es.res$count == 1) {
        query.split <- list(0)
    }
    if (rtm>0) {
        query.split <- split(0:es.res$count,cut(0:es.res$count,ceiling(es.res$count/rtm),labels=FALSE))
    }
    x <- mapply(function(x) entrez_fetch_parallel(webhist=es.post,chunks=x,file=fasout), x=query.split,SIMPLIFY=TRUE,USE.NAMES=FALSE)
    es.fas <- ape::read.FASTA(fasout)
    #writeLines(paste("\nDone!",length(es.fas),"sequences downloaded in FASTA format.\n"))
    cli_report(txt=glue::glue("n={length(es.fas)} sequences downloaded in FASTA format."),rule=FALSE,alert="info")
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
        stop(cli::cli_alert_danger("ERROR Searches failed ... aborted.")) 
    } else {
        end_time <- Sys.time()
        #writeLines(paste0("Metadata for ",length(accs)," accessions downloaded (starting ",accs[1],"). Download took ",round(as.numeric(end_time-start_time),digits=2)," seconds."))
        cli_report(txt=glue::glue("Metadata for {length(accs)} accessions downloaded (starting {accs[1]}. Download took {round(as.numeric(end_time-start_time),digits=2)} seconds."),rule=FALSE,alert="success")
        return(ncbi.tab)
    }
}


# DEREPLICATE FASTA
dereplicate_fasta <- function(infile,dereplicate,threads) {
    outfile <- stringr::str_replace_all(infile,"\\.fasta",".derep.fasta")
        if (dereplicate=="true") {
            exe.string <- glue::glue("vsearch --threads {threads} --fastx_uniques {infile} --fastaout {outfile}")
            system(exe.string)
        } else if (dereplicate=="false") {
            file.copy(infile,outfile)
        } else stop(cli::cli_alert_danger("ERROR the dereplicate flag '-d' must be 'true' or 'false'."))
}


# FILTER FASTA
filter_fasta <- function(infile,maxns,threads) {
    infile <- stringr::str_replace_all(infile,"\\.fasta",".derep.fasta")
    outfile <- stringr::str_replace_all(infile,"\\.fasta",".filtered.fasta")
    exe.string <- glue::glue("vsearch --threads {threads} -fastx_filter {infile} --fastq_maxns {maxns} --fastaout {outfile}")
    system(exe.string)
}


# CLUSTER FASTA
cluster_fasta <- function(infile,identity,threads) {
    infile <- stringr::str_replace_all(infile,"\\.fasta",".derep.filtered.fasta")
    outfile <- here::here(today.dir,"cluster.")
    stringr::str_replace_all(infile,"genbank-dump.derep.filtered.fasta","cluster.")
    exe.string <- glue::glue("vsearch --threads {threads} --cluster_fast {infile} --id {identity} --clusters {outfile} --strand both")
    system(exe.string)
}


# LOAD AND RENAME FASTA
load_and_rename <- function(gene,path,wd) {
    fas <- ape::read.FASTA(path)
    path.gene <- here::here(today.dir,paste("gene",gene,"fasta",sep="."))
    ape::write.FASTA(fas,path.gene)
}


# TABULATE GENES
tabulate_genes <- function(path) {
    fas <- ape::read.FASTA(path)
    accs <- names(fas)
    gene <- stringr::str_split_fixed(basename(path),"\\.",3)[,2]
    accs.gene <- tibble::tibble(acc_no=accs,gene=gene)
    return(accs.gene)
}


# CLEAN NCBI DATA
clean_ncbi <- function(df) {
    df.clean <- df |>
        dplyr::distinct(gi_no,.keep_all=TRUE) |> 
        # filter
        dplyr::filter(gi_no!="NCBI_GENOMES") |> 
        dplyr::filter(!is.na(sequence)) |> 
        dplyr::filter(!grepl("UNVERIFIED:",gene_desc)) |>
        dplyr::filter(!grepl("PREDICTED:",gene_desc)) |>
        dplyr::filter(!grepl("similar to",gene_desc)) |> 
        dplyr::filter(!grepl("mRNA",gene_desc)) |> 
        dplyr::filter(!grepl("cDNA",gene_desc)) |> 
        dplyr::filter(!grepl("transcribed",gene_desc)) |> 
        dplyr::filter(!grepl("-like",gene_desc)) |>
        # fix lan lon
        dplyr::mutate(lat=paste(stringr::str_split_fixed(lat_lon, " ", 4)[,1], stringr::str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(stringr::str_split_fixed(lat_lon, " ", 4)[,3], stringr::str_split_fixed(lat_lon, " ", 4)[,4])) |>
        dplyr::mutate(lat=dplyr::if_else(grepl(" N",lat), true=stringr::str_replace_all(lat," N",""), false=dplyr::if_else(grepl(" S",lat), true=paste0("-",stringr::str_replace_all(lat," S","")), false=lat))) |>
        dplyr::mutate(lon=dplyr::if_else(grepl(" E",lon), true=stringr::str_replace_all(lon," E",""), false=dplyr::if_else(grepl(" W",lon), true=paste0("-",stringr::str_replace_all(lon," W","")), false=lon))) |> 
        dplyr::mutate(lat=stringr::str_replace_all(lat,"^ ", NA_character_), lon=stringr::str_replace_all(lon,"^ ", NA_character_)) |>
        dplyr::mutate(lat=suppressWarnings(as.numeric(lat)), lon=suppressWarnings(as.numeric(lon))) |> 
        # tidy up
        dplyr::mutate(length=as.numeric(length)) |>
        dplyr::mutate(gi_no=as.character(gi_no)) |>
        dplyr::mutate(scientificNameGenBank=taxon) |>
        dplyr::rename(scientificName=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence) |>
        dplyr::select(scientificName,scientificNameGenBank,class,order,family,genus,dbid,gbAccession,gene,length,organelle,catalogNumber,country,publishedAs,publishedIn,publishedBy,date,decimalLatitude,decimalLongitude,nucleotides)
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
    gstring <- if (n == 2) {
            glue::glue("{splitz1}_{splitz2}")
        } else if (n == 3) {
            glue::glue("{splitz1}_{splitz2}_{splitz3}")
        } else if (n == 4) {
            glue::glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}")
        } else if (n == 5) {
            glue::glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}_{splitz5}")
        } else {stop(cli::cli_alert_danger("ERROR Names value must be between 2 and 5."))
        }
#    gstring <- case_when(
#        n==2 ~ glue("{splitz1}_{splitz2}"),
#        n==3 ~ glue("{splitz1}_{splitz2}_{splitz3}"),
#        n==4 ~ glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}"),
#        n==5 ~ glue("{splitz1}_{splitz2}_{splitz3}_{splitz4}_{splitz5}")
#        )
    return(gstring)
}


# CLEAN NCBI DATA
clean_names <- function(df,n) {
    df.clean <- df |>
        dplyr::mutate(label=scientificName) |>
        dplyr::mutate(label=stringr::str_replace_all(label," ","_")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"'","")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,",","")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"\\(","")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"\\)","")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,":","-")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"cf\\._","cf.")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"aff\\._","aff.")) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"sp\\._","sp.")) |>
        dplyr::mutate(elements=stringr::str_count(label,"_")+1) |>
#        mutate(label=case_when(
#            n == 2 ~ splitter(label,n=2),
#            n == 3 ~ splitter(label,n=3),
#            n == 4 ~ splitter(label,n=4),
#            n == 5 ~ splitter(label,n=5),
#            )) |>
mutate(label=
    if (n == 2) {splitter(label,n=2)} 
    else if (n == 3) {splitter(label,n=3)}
    else if (n == 4) {splitter(label,n=4)}
    else if (n == 5) {splitter(label,n=5)}
    else {stop(cli::cli_alert_danger("ERROR Names value must be between 2 and 5."))}
) |>
        dplyr::mutate(label=as.character(label)) |>
        dplyr::mutate(label=stringr::str_replace_all(label,"_+$","")) |>
        dplyr::mutate(scientificName=label) |>
        dplyr::select(!c(label,elements))
    return(df.clean)
}


# ANNOTATE WITH FISHBASE TAXONOMY
fishbase_annotate <- function(genera) {
    genera.taxonomy <- rfishbase::fb_tbl("classification_tree",server="fishbase",version="latest") |> 
        dplyr::distinct(Class,Order,Family,Genus) |>
        dplyr::rename(class=Class,order=Order,family=Family,genus=Genus)
    genera.taxonomy <- genera.taxonomy |> dplyr::filter(genus %in% genera)
    return(genera.taxonomy)
}


# ANNOTATE WITH NCBI TAXONOMY
ncbi_annotate <- function(genera) {
    tax <- c("class","order","family")
    tax.tab <- taxize::tax_name(sci=genera,get=tax,db="ncbi")
    tax.tab.tib <- tax.tab |> tibble::as_tibble() |> dplyr::select(-db) |> dplyr::rename(genus=query) |> dplyr::relocate(genus,.after="family")
    return(tax.tab.tib)
}


# WRITE A FASTA FILE FROM DF
write_fasta <- function(df,genez,dir,pop) {
    if(isFALSE(pop)) {
        df.sub <- df |> dplyr::filter(gene == genez) |> tab2fas(seqcol="nucleotides",namecol="scientificName")
    } else if(isTRUE(pop)) {
        df.sub <- df |> dplyr::filter(gene == genez) |> tab2fas(seqcol="nucleotides",namecol="gbAccession")
    }
    df.sub |> ape::write.FASTA(file=here::here(dir,glue::glue("{genez}.fasta")))
    #writeLines(glue("Unaligned fasta file written to 'temp/{basename(dir)}/{genez}.fasta'.",.trim=FALSE))
}


# ALIGN FASTA
align_fasta <- function(infile,threads) {
    outfile <- stringr::str_replace_all(infile,"\\.fasta",".aligned.fasta")
    exe.string <- glue::glue("mafft --adjustdirection --quiet --auto --thread {threads} {infile} > {outfile}")
    system(exe.string)
    #writeLines(glue("Aligned fasta file written to '{basename(outfile)}'.",.trim=FALSE))
}


# TRIM FASTA
trim_fasta <- function(infile,prop) {
    infile <- stringr::str_replace_all(infile,"\\.fasta",".aligned.fasta")
    outfile <- stringr::str_replace_all(infile,"\\.fasta",".trimmed.fasta")
    exe.string <- glue::glue("trimal -in {infile} -out {outfile} -gt {prop}")
    system(exe.string)
    cli::cli_alert_info(glue::glue("Aligned and trimmed fasta gene file written to '{basename(outfile)}'."))
    #cli_report(txt=glue::glue("Aligned and trimmed fasta file written to '{basename(outfile)}'."),rule=FALSE,alert="info")
    #writeLines(glue::glue("\nAligned and trimmed fasta file written to '{basename(outfile)}'.\n",.trim=FALSE))
}


# MAKE PARTS
partition_table <- function(mat) {
    dims <- purrr::map(mat, \(x) dim(x)[2]) |> unlist()
    parts.tib <- tibble::tibble(gene=genes,len=dims,cs=cumsum(dims)) |> 
        dplyr::mutate(cs1=cs+1,start=lag(cs1)) |>
        tidyr::replace_na(list(start=1)) |>
        dplyr::rename(end=cs) %>%
        dplyr::select(gene,start,end) %>%
        dplyr::mutate(parts=as.character(glue::glue("DNA, {gene} = {start}-{end}"))) %>%
        dplyr::select(parts) 
    return(parts.tib)
}


# NEW RAXML-NG FUN
raxml_ng <- function(file,model,maxthreads,epsilon,verbose) {
    string.search <- glue::glue("raxml-ng --search -threads auto{{{maxthreads}}} --workers auto --tree pars{{1}} --redo --seed 42 --model {model} --lh-epsilon {epsilon} --msa {file}")
    if (verbose == "true") {
        system(command=string.search,ignore.stdout=FALSE)
        } else if (verbose == "false") {
        system(command=string.search,ignore.stdout=TRUE)
        } else stop(cli::cli_alert_danger("ERROR the '-v' flag must be 'true' or 'false'."))
    rax.tr <- ape::read.tree(file=glue::glue("{file}.raxml.bestTree"))
    cli::cli_alert_info(glue::glue("Tree search completed for '{basename(file)}'."))
    #writeLines(glue::glue("\nTree search completed for '{basename(file)}'.",.trim=FALSE))
    return(rax.tr)
}


# MAKE TREE PLOTTING FUN
ggtree_autoplot <- function(path,tb,scale.factor,width,hratio) {
    tr <- treeio::read.newick(path)
    tr <- tr |> castor::root_in_edge(root_edge=which.max(tr$edge.length))
    #tr <- phangorn::midpoint(tr)
    tree.length <- max(castor::get_all_distances_to_root(tr)) + (max(castor::get_all_distances_to_root(tr)) * width)
        p <- ggtree::ggtree(tr, ladderize=TRUE,right=TRUE,linewidth=0.7) %<+% tb
        pp <- p +  ggtree::geom_tiplab(offset=0,aes(label=tiplabel),align=FALSE,size=4) +
        ggtree::geom_tippoint(aes(color=tip.colour),size=2.5) +
        theme(legend.position="none") +
        xlim(0,tree.length)
    ggsave(filename=glue::glue("{path}.pdf"),plot=pp,limitsize=FALSE,width=200,height=200*hratio,units="mm",scale=scale.factor)
}


# MAKE A REPORTING FUNCTION
cli_report <- function(txt,rule,alert) {
    if(isTRUE(rule)) {
        cli::cli_rule()
            if(alert=="info") {
                cli::cli_alert_info(txt) 
            } else if(alert=="success") {
                cli::cli_alert_success(txt)
            }
        cli::cli_rule()
    } else if (isFALSE(rule)) {
        cli::cli_text("")
            if(alert=="info") {
                cli::cli_alert_info(txt) 
            } else if(alert=="success") {
                cli::cli_alert_success(txt)
            }
        cli::cli_text("")
    }
}

# report
cli_report(txt="Packages and functions loaded.",rule=FALSE,alert="success")
#writeLines("\nPackages and functions loaded.\n")

