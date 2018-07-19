Functional Enrichment on DE data
================
Nick Garza
7/2/2018

Load required functions and set directories
===========================================

Notable used include: DESeq2 and enrichR

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

Example data: platelet derived growth factor treatment versus vehicle treatment at different times

``` r
# data
load("pdgf_veh_t1hgenes.RData")
load("pdgf_veh_t4hgenes.RData")
load("pdgf_veh_t24hgenes.RData")
```

Loading functional expression databases through enrichR
=======================================================

Query databases to find genes of interest that are differentially expressed at 1hr, 4hr, 24hr in response to inflammatory stimulation

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

Here we use a p-value filter function with a parameter to change the level of stringency. For our example we use a value of p = 0.05

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

We generate two data frames summarizing data: total\_enrich\_df and reduced\_df. reduced\_df has the info total\_enrich\_df limited to database, description, p-value, adjusted p-value, and 3 genes of interest.

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

head(total_enrich_df)
```

    ##             .id
    ## 1 Reactome_2016
    ## 2 Reactome_2016
    ## 3 Reactome_2016
    ## 4 Reactome_2016
    ## 5 Reactome_2016
    ## 6 Reactome_2016
    ##                                                                  Term
    ## 1 Molecules associated with elastic fibres_Homo sapiens_R-HSA-2129379
    ## 2                  Elastic fibre formation_Homo sapiens_R-HSA-1566948
    ## 3                   Sialic acid metabolism_Homo sapiens_R-HSA-4085001
    ## 4                   Glutathione conjugation_Homo sapiens_R-HSA-156590
    ## 5        Extracellular matrix organization_Homo sapiens_R-HSA-1474244
    ## 6                              cGMP effects_Homo sapiens_R-HSA-418457
    ##   Overlap      P.value Adjusted.P.value  Old.P.value Old.Adjusted.P.value
    ## 1    5/30 0.0004147903        0.1700640 0.0001311745           0.05378155
    ## 2    5/41 0.0017965240        0.3682874 0.0004880547           0.06670080
    ## 3    4/33 0.0053081981        0.6014805 0.0018720533           0.19188546
    ## 4    4/38 0.0088021531        0.6014805 0.0030018541           0.21892063
    ## 5  13/283 0.0087308242        0.6014805 0.0003244101           0.06650406
    ## 6    3/18 0.0063744190        0.6014805 0.0032037165           0.21892063
    ##     Z.score Combined.Score
    ## 1 -1.973547      15.369464
    ## 2 -1.964402      12.418753
    ## 3 -2.081774      10.905381
    ## 4 -2.096794       9.923620
    ## 5 -2.079843       9.860317
    ## 6 -1.905378       9.632568
    ##                                                                                 Genes
    ## 1                                                     EFEMP2;EFEMP1;ITGB8;FBLN1;FBLN5
    ## 2                                                     EFEMP2;EFEMP1;ITGB8;FBLN1;FBLN5
    ## 3                                               ST6GALNAC1;ST8SIA1;ST8SIA5;ST6GALNAC4
    ## 4                                                             GSTM4;GSTM2;GSTM1;GSTM5
    ## 5 ICAM2;FBLN1;LAMB1;ICAM5;FBLN5;CAPN11;MMP25;EFEMP2;EFEMP1;COL21A1;TLL2;ITGB8;ADAMTS9
    ## 6                                                                  PDE10A;MRVI1;PDE5A

``` r
head(reduced_df)
```

    ##             .id
    ## 1 Reactome_2016
    ## 2 Reactome_2016
    ## 3 Reactome_2016
    ## 4 Reactome_2016
    ## 5 Reactome_2016
    ## 6 Reactome_2016
    ##                                                                  Term
    ## 1 Molecules associated with elastic fibres_Homo sapiens_R-HSA-2129379
    ## 2                  Elastic fibre formation_Homo sapiens_R-HSA-1566948
    ## 3                   Sialic acid metabolism_Homo sapiens_R-HSA-4085001
    ## 4                   Glutathione conjugation_Homo sapiens_R-HSA-156590
    ## 5        Extracellular matrix organization_Homo sapiens_R-HSA-1474244
    ## 6                              cGMP effects_Homo sapiens_R-HSA-418457
    ##        P.value Adjusted.P.value                      Genes
    ## 1 0.0004147903        0.1700640        EFEMP2,EFEMP1,ITGB8
    ## 2 0.0017965240        0.3682874        EFEMP2,EFEMP1,ITGB8
    ## 3 0.0053081981        0.6014805 ST6GALNAC1,ST8SIA1,ST8SIA5
    ## 4 0.0088021531        0.6014805          GSTM4,GSTM2,GSTM1
    ## 5 0.0087308242        0.6014805          ICAM2,FBLN1,LAMB1
    ## 6 0.0063744190        0.6014805         PDE10A,MRVI1,PDE5A

``` r
# can write total_enrich_df and recued_df into csv's, lists as we please
```
