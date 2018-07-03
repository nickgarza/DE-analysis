Functional Enrichment on DE data
================
Nick Garza
7/2/2018

Load required functions and set directories
===========================================

``` r
# Investigating functional gene expression at different time points
library(DESeq2)
library(lattice)
library(bigWig)
library(dplyr)
library(enrichR)


# load functions and data
function.dir = paste0(base.dir,"/functions")
fn.name = "deseq_functions4.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "data_functions3.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "bedtools.R"
source( paste0(function.dir,"/",fn.name) )
```

Load data sets to look for functional enrichment
================================================

``` r
# data
load("pdgf_veh_t1hgenes.RData")
load("pdgf_veh_t4hgenes.RData")
load("pdgf_veh_t24hgenes.RData")
```

Loading functional expression databases through enrichR
=======================================================

``` r
# functional expression analysis through enrichR
dbs <- listEnrichrDbs()
dbs_ind = c(93:95,99,102,115:117)
Db = dbs[dbs_ind,1] 
enriched_1hr <- enrichr(t1hgenes, Db) %>% as.list
```

    ## Uploading data to Enrichr... Done.
    ##   Querying Reactome_2016... Done.
    ##   Querying KEGG_2016... Done.
    ##   Querying WikiPathways_2016... Done.
    ##   Querying BioCarta_2016... Done.
    ##   Querying Panther_2016... Done.
    ##   Querying GO_Cellular_Component_2017b... Done.
    ##   Querying GO_Molecular_Function_2017b... Done.
    ##   Querying GO_Biological_Process_2017b... Done.
    ## Parsing results... Done.

``` r
enriched_4hr <- enrichr(t4hgenes, Db) %>% as.list 
```

    ## Uploading data to Enrichr... Done.
    ##   Querying Reactome_2016... Done.
    ##   Querying KEGG_2016... Done.
    ##   Querying WikiPathways_2016... Done.
    ##   Querying BioCarta_2016... Done.
    ##   Querying Panther_2016... Done.
    ##   Querying GO_Cellular_Component_2017b... Done.
    ##   Querying GO_Molecular_Function_2017b... Done.
    ##   Querying GO_Biological_Process_2017b... Done.
    ## Parsing results... Done.

``` r
enriched_24hr <- enrichr(t24hgenes, Db) %>% as.list 
```

    ## Uploading data to Enrichr... Done.
    ##   Querying Reactome_2016... Done.
    ##   Querying KEGG_2016... Done.
    ##   Querying WikiPathways_2016... Done.
    ##   Querying BioCarta_2016... Done.
    ##   Querying Panther_2016... Done.
    ##   Querying GO_Cellular_Component_2017b... Done.
    ##   Querying GO_Molecular_Function_2017b... Done.
    ##   Querying GO_Biological_Process_2017b... Done.
    ## Parsing results... Done.

``` r
head(enriched_1hr[[1]][,c(1,3)])
```

    ##                                                                  Term
    ## 1 Molecules associated with elastic fibres_Homo sapiens_R-HSA-2129379
    ## 2                  Elastic fibre formation_Homo sapiens_R-HSA-1566948
    ## 3                   Sialic acid metabolism_Homo sapiens_R-HSA-4085001
    ## 4                   Glutathione conjugation_Homo sapiens_R-HSA-156590
    ## 5        Extracellular matrix organization_Homo sapiens_R-HSA-1474244
    ## 6                              cGMP effects_Homo sapiens_R-HSA-418457
    ##        P.value
    ## 1 0.0004147903
    ## 2 0.0017965240
    ## 3 0.0053081981
    ## 4 0.0088021531
    ## 5 0.0087308242
    ## 6 0.0063744190

Filter databases by FDR cutoff
==============================

``` r
# filter by FDR using list method
fdr_filter <- function(enriched_list=NULL, pcut=NULL){
  filtered <- lapply(enriched_list,function(x){
    xnew = x %>% filter(P.value < pcut)
  })
  
  return(filtered)
}
enrich_sig_1hr = fdr_filter(enriched_list = enriched_1hr, pcut = 0.05)
enrich_sig_4hr = fdr_filter(enriched_list = enriched_4hr, pcut = 0.05)
enrich_sig_24hr = fdr_filter(enriched_list = enriched_24hr, pcut = 0.05)
```

Move into a data frame and reduced data frame
=============================================

``` r
# concatenate the dataframes
library(plyr)
enrich_sig_1hr_df <- ldply(enrich_sig_1hr, data.frame)
enrich_sig_4hr_df <- ldply(enrich_sig_4hr, data.frame)
enrich_sig_24hr_df <- ldply(enrich_sig_24hr, data.frame)

total_enrich_df <- rbind(enrich_sig_1hr_df,enrich_sig_4hr_df,enrich_sig_24hr_df)
reduced_df <- select(total_enrich_df, .id, Term, P.value, Adjusted.P.value, Genes)

gene_list = reduced_df[5]

# first 3 genes
for(ii in 1:nrow(gene_list)){
  separate = strsplit(gene_list[ii,1],";")[[1]][1:3]
  gene_list[ii,1] = paste(paste(separate[1],separate[2], sep = ","), separate[3], sep = ",")

}
# paste back into data frame
reduced_df[5] = gene_list
```
