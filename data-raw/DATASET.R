#!/usr/bin/env Rscript

library(usethis)
library(dplyr)

# Helper function to input the data objects
generate_data <- function(
    Row
    ){

    # What is the name of the object?
    data_object <- paste0(
            Row[["dataset"]],
            "_",
            Row[["focus"]],
            "_v",
            Row[["version"]]
        )

    # Where is it located?
    data_file <- paste0(
            getwd(),
            "/inst/extdata/",
            Row[["dataset"]],
            "/",
            Row[["version"]],
            "/",
            Row[["focus"]],
            ".tsv"
        )

    # Use the name of the object to store the associated data
    assign(
        data_object,
        read.delim(
            data_file
        )
    )
    # Generate .rda object with the data in it
    do.call(
        "use_data",
        list(as.name(data_object),
        overwrite = TRUE)
    )
}

# Coordinate-based data objects
# All possible names of datasets
datasets <- c(
    "somatic_hypermutation_locations"
)

# All possible genome builds
genome_builds <- c(
    "GRCh37",
    "GRCh38"
)

# All possible versions
# This will need to be modified with the release of new data
versions <- c(
    "0.0",
    "0.1",
    "0.2",
    "0.3",
    "_latest"
)

# All possible combinations of data/genome build/version
all_coordinates <- expand.grid(
    dataset = datasets,
    focus = genome_builds,
    version = versions,
    stringsAsFactors = FALSE)


# Pathology-based data objects
# All possible names of datasets
datasets <- c(
    "lymphoma_genes"
)

# All possible genome builds
pathologies <- c(
    "bl",
    "dlbcl",
    "mcl"
)

# All possible versions
# This will need to be modified with the release of new data
versions <- c(
    "0.1",
    "_latest"
)

# All possible combinations of data/pathology/version
all_pathologies <- expand.grid(
    dataset = datasets,
    focus = pathologies,
    version = versions,
    stringsAsFactors = FALSE)

all <- bind_rows(
    all_coordinates,
    all_pathologies
)

# Manually include legacy datasets
all <- bind_rows(
    all,
    tibble(
        dataset = c("lymphoma_genes"),
        focus = c("lymphoma_genes"),
        version = c("0.0")
    )
)

# Generate the rdas for each object
apply(
    all,
    1,
    generate_data
)


grch37_ashm_regions <- read.delim(
    paste0(
                getwd(),
                "/inst/extdata/",
                "somatic_hypermutation_locations",
                "/",
                "_latest",
                "/",
                "GRCh37",
                ".tsv"
            )
)

hg38_ashm_regions <- read.delim(
    paste0(
                getwd(),
                "/inst/extdata/",
                "somatic_hypermutation_locations",
                "/",
                "_latest",
                "/",
                "GRCh38",
                ".tsv"
            )
)

#create alias for lymphoma genes
lymphoma_genes <- read.delim(
    paste0(
                getwd(),
                "/inst/extdata/",
                "lymphoma_genes",
                "/",
                "0.0",
                "/",
                "lymphoma_genes",
                ".tsv"
            )
)

use_data(grch37_ashm_regions, overwrite = TRUE)
use_data(hg38_ashm_regions, overwrite = TRUE)
use_data(lymphoma_genes, overwrite = TRUE)

#migrated resources
gene_blacklist = system.file("extdata", "gene_blacklist_with_IG.tsv", package = "GAMBLR.data") %>%
  read_tsv()

usethis::use_data(gene_blacklist, overwrite = TRUE)

grch37_all_gene_coordinates = system.file("extdata", "grch37_gene_coordinates.tsv", package = "GAMBLR.data") %>%
  read_tsv() %>% dplyr::filter(grepl("PATCH",chromosome))

usethis::use_data(grch37_all_gene_coordinates, overwrite = TRUE)


hg38_oncogene = system.file("extdata", "oncogene_regions.hg38.tsv", package = "GAMBLR.data") %>%
  read_tsv(col_types="ciici")

usethis::use_data(hg38_oncogene, overwrite = TRUE)

grch37_oncogene = system.file("extdata", "oncogene_regions.grch37.tsv", package = "GAMBLR.data") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_oncogene, overwrite = TRUE)

hg38_partners = system.file("extdata","superenhancer_regions.hg38.tsv",package="GAMBLR.data") %>%
  read_tsv(col_types="ciici")

usethis::use_data(hg38_partners, overwrite = TRUE)

grch37_partners = system.file("extdata", "superenhancer_regions.grch37.tsv", package = "GAMBLR.data") %>%
  read_tsv(col_types="ciici")


#THis needs to be syncronized with the list above (it currently is out of sync!)
grch37_lymphoma_genes_bed = system.file("extdata","lymphoma_genes.grch37.bed",package="GAMBLR.data") %>%
  read_tsv()

usethis::use_data(grch37_lymphoma_genes_bed, overwrite = TRUE)

hg38_lymphoma_genes_bed = system.file("extdata","lymphoma_genes.hg38.bed",package="GAMBLR.data") %>%
  read_tsv()

usethis::use_data(hg38_lymphoma_genes_bed, overwrite = TRUE)


wright_genes_with_weights = system.file("extdata","WrightGenesWithWeights.txt",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1) %>% rename(Ensembl_ID=EnsemblGeneID,Hugo_Symbol=GeneName)

usethis::use_data(wright_genes_with_weights, overwrite = TRUE)

dhitsig_genes_with_weights = system.file("extdata","DHITsigGenesWithWeights.txt",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1) %>% rename(Ensembl_ID=ensembl_gene_id,Hugo_Symbol=GeneName)
usethis::use_data(dhitsig_genes_with_weights, overwrite = TRUE)

target_regions_hg38 = system.file("extdata","target_regions_hg38.txt",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1)
usethis::use_data(target_regions_hg38, overwrite = TRUE)

target_regions_grch37 = system.file("extdata","target_regions_grch37.txt",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1)
usethis::use_data(target_regions_grch37, overwrite = TRUE)

hotspot_regions_grch37 = system.file("extdata","hotspot_regions.grch37.tsv",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1) %>%
  column_to_rownames("gene")
usethis::use_data(hotspot_regions_grch37, overwrite = TRUE)

hotspot_regions_hg38 = system.file("extdata","hotspot_regions.hg38.tsv",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1) %>%
  column_to_rownames("gene")
usethis::use_data(hotspot_regions_hg38, overwrite = TRUE)

chromosome_arms_grch37 = system.file("extdata","chromosome_arms_grch37.tsv",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1)
usethis::use_data(chromosome_arms_grch37, overwrite = TRUE)

chromosome_arms_hg38 = system.file("extdata","chromosome_arms_hg38.tsv",package="GAMBLR.data") %>%
  read.table(sep="\t",header=1)
usethis::use_data(chromosome_arms_hg38, overwrite = TRUE)

lymphgen_entrez = system.file("extdata","lymphgen_genes_entrez.txt",package="GAMBLR.data") %>%
  read_tsv()

entrez_map = system.file("extdata","hugo2entrez.tsv",package="GAMBLR.data") %>%
  read_tsv()

lymphgen_anno = left_join(lymphgen_entrez,entrez_map) %>% dplyr::rename("Hugo_Symbol"="Approved symbol") %>% dplyr::select(1:3)


library("biomaRt")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
#need to get entrezgene_id, hgnc_symbol using ensembl_gene_id
gene_detail = getBM(attributes=c( 'ensembl_gene_id','entrezgene_id','hgnc_symbol'),
      filters = 'ensembl_gene_id',
      values = lymphoma_genes$ensembl_gene_id,
      mart = ensembl,useCache = FALSE)

#fill in missing entrezgene_id

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000036448","entrezgene_id"] = 9172

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000065526","entrezgene_id"] = 23013

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000102096","entrezgene_id"] = 11040

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000172578","entrezgene_id"] = 89857

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000205542","entrezgene_id"] = 7114

lymphoma_genes = left_join(lymphoma_genes,gene_detail)
lymphoma_genes$LymphGen=FALSE
lymphoma_genes[lymphoma_genes$entrezgene_id %in% lymphgen_anno$NCBI_Gene_ID,"LymphGen"] = TRUE

#load the Reddy gene list
reddy_genes = system.file("extdata","reddy_genes.tsv",package="GAMBLR") %>%
  read_tsv() %>% dplyr::rename("hgnc_symbol"="Approved symbol")

usethis::use_data(reddy_genes, overwrite = TRUE)

lymphoma_genes$Reddy = FALSE
lymphoma_genes[lymphoma_genes$hgnc_symbol %in% reddy_genes$hgnc_symbol,"Reddy"]=TRUE

usethis::use_data(lymphoma_genes, overwrite = TRUE)

reddy_only = reddy_genes[which(!reddy_genes$hgnc_symbol %in% lymphoma_genes$hgnc_symbol),"hgnc_symbol"]

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

reddy_detail = getBM(attributes=c( 'ensembl_gene_id','entrezgene_id','hgnc_symbol'),
                    filters = 'hgnc_symbol',
                    values = reddy_only$hgnc_symbol,
                    mart = ensembl,useCache = FALSE)

#update any based on mapping to Ensembl ID
lymphoma_genes[lymphoma_genes$ensembl_gene_id %in% reddy_detail$ensembl_gene_id,"Reddy"]=TRUE

#reddy_only = dplyr::filter(reddy_detail,!ensembl_gene_id %in% lymphoma_genes$ensembl_gene_id ) %>%
#  group_by(ensembl_gene_id,hgnc_symbol) %>% slice_head() %>% ungroup() %>% dplyr::select(-entrezgene_id) %>% dplyr::rename("Gene"="hgnc_symbol") %>%
#  mutate(Reddy=TRUE)
#lymphoma_genes_comprehensive = bind_rows(reddy_only,lymphoma_genes) %>% dplyr::select(ensembl_gene_id,Gene,Reddy,LymphGen,Chapuy)



chapuy_genes = system.file("extdata","chapuy_genes.tsv",package="GAMBLR.data") %>%
  read_tsv() %>% dplyr::rename("hgnc_symbol"="gene")

lymphoma_genes$Chapuy = FALSE
lymphoma_genes[lymphoma_genes$hgnc_symbol %in% chapuy_genes$hgnc_symbol,"Chapuy"]=TRUE

#lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% chapuy_genes$hgnc_symbol,"Chapuy"]=TRUE
#lymphoma_genes_comprehensive[!lymphoma_genes_comprehensive$Gene %in% chapuy_genes$hgnc_symbol,"Chapuy"]=FALSE


lymphoma_genes_comprehensive = read_tsv("inst/extdata/lymphoma_genes_comprehensive.tsv")

lacy = read_tsv("inst/extdata/lacy_genes.tsv") %>% dplyr::filter(`Included in statistical analysis`=='Yes',Feature!="Amplification")
lacy_genes = pull(lacy,Gene)
lymphoma_genes$Lacy = FALSE
lymphoma_genes[lymphoma_genes$hgnc_symbol %in% lacy_genes,"Lacy"]=TRUE

lymphoma_genes_comprehensive$Lacy = FALSE
lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% lacy_genes,"Lacy"]=TRUE

lacy_ashm = dplyr::filter(lacy,!is.na(`Annotation as aSHM`)) %>% pull(Gene)
lymphoma_genes_comprehensive$aSHM = FALSE
lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% lacy_ashm,"aSHM"]=TRUE
lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% grch37_ashm_regions$gene,"aSHM"]=TRUE
lymphoma_genes_comprehensive = mutate(lymphoma_genes_comprehensive,aSHM=ifelse(grepl("HIST",Gene),TRUE,aSHM))
usethis::use_data(lymphoma_genes_comprehensive, overwrite = TRUE)



#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
#need to get entrezgene_id, hgnc_symbol using ensembl_gene_id
#gene_detail = getBM(attributes=c( 'ensembl_gene_id','hgnc_symbol','chromosome_name'),
#                    filters = 'hgnc_symbol',
#                    values = chapuy_only$hgnc_symbol,
#                    mart = ensembl,useCache = FALSE) %>% dplyr::filter(chromosome_name %in% c(c(1:22),"X","Y")) %>%
#  dplyr::select(-chromosome_name) %>% dplyr::rename("Gene"="hgnc_symbol") %>%
#  mutate(Chapuy=TRUE) %>% dplyr:: filter(!Gene %in% lymphoma_genes_comprehensive$Gene)

#lymphoma_genes_comprehensive = bind_rows(gene_detail,lymphoma_genes_comprehensive) %>%
#  mutate(Reddy=ifelse(is.na(Reddy),FALSE,Reddy))

#lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% lymphgen_anno$Hugo_Symbol,"LymphGen"] = TRUE
#lymphoma_genes_comprehensive[!lymphoma_genes_comprehensive$Gene %in% lymphgen_anno$Hugo_Symbol,"LymphGen"] = FALSE

#lymphoma_genes_comprehensive %>% dplyr::filter(Chapuy ==FALSE, Reddy==FALSE, LymphGen == FALSE)

lymphoma_genes = dplyr::select(lymphoma_genes,-entrezgene_id) %>% group_by(Gene,ensembl_gene_id) %>% slice_head() %>% ungroup()

#DLBCL_curated_genes = lymphoma_genes %>% dplyr::filter(DLBCL==TRUE) %>% pull(Gene)

#lymphoma_genes_comprehensive$curated = FALSE
#lymphoma_genes_comprehensive[lymphoma_genes_comprehensive$Gene %in% DLBCL_curated_genes,]$curated = TRUE
#lymphoma_genes_comprehensive = dplyr::filter(lymphoma_genes_comprehensive,Reddy==TRUE | Chapuy == TRUE | LymphGen == TRUE | curated == TRUE)
#lymphoma_genes_comprehensive$other_support = ""

#download gene annotations
#load library
library(curl)

#get data
curl::curl_download(url = "ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz",
                    destfile = "../gene_ann/Homo_sapiens.GRCh38.86.gtf.gz",
                    quiet = FALSE)

#read gtf into R
gtf_hg38 = rtracklayer::import('../gene_ann/Homo_sapiens.GRCh38.86.gtf.gz')

#tidy
gene_annotations_hg38 = as.data.frame(gtf_hg38) %>% #convert to data frame
  dplyr::filter(type == "gene") %>% #only keep genes
  dplyr::select(gene_id, seqnames, start, end, gene_name) %>% #select relevant columns
  dplyr::distinct(gene_id, .keep_all = TRUE) %>% #only keep unique elements (based on gene_id column)
  dplyr::mutate(seqnames = paste0("chr", seqnames)) %>% #add chromosome (seqnames) prefix
  dplyr::mutate(seqnames = as.factor(seqnames)) %>% #convert chromosome ( seqnames) variable to factor
  dplyr::rename(ensembl_gene_id = gene_id, chromosome = seqnames) %>% #rename columns to amtch expected format
  dplyr::mutate(hugo_symbol = gene_name) %>% #duplicate gene_name and rename it to hugo_symbol
  dplyr::arrange(chromosome, start) #sort data frame on chromosom, then start coordinates

#save gene coordinates to file (rda)
save(gene_annotations_hg38, file = "data/hg38_gene_coordinates.rda")
