Venn Diagrams of DE data
================
Nick Garza
7/2/2018

Loading our libraries
---------------------

Load required backages and bedtool/bigwig functions we need for computations

``` r
## load functions and diff expr data
# load packages
library(DESeq2)
library(lattice)
library(bigWig)
library(dplyr)
library(tibble)
library(VennDiagram)

# load functions
function.dir = paste0(base.dir,"/functions")
fn.name = "deseq_functions4.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "data_functions3.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "bedtools.R"
source( paste0(function.dir,"/",fn.name) )

#filter degdata with significant candidates
load("degdata.RData")

deg_list = lapply(degdata,function(x){
  degs = x[["degAll"]] %>% as.data.frame %>% filter(abs(log2FoldChange)>1, padj<0.001)
  return(degs)
})
```

Load Differential Comparisons of interest
=========================================

Read in differential comparisons of interest from differential expression matrix

``` r
## acquire significant DE genes
t1h_list = degdata[["t1h_Veh_&_t1h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t1h_list_degs = t1h_list %>% filter(abs(log2FoldChange)>1, padj<0.001) # sig DE genes

t4h_list = degdata[["t4h_Veh_&_t4h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t4h_list_degs = t4h_list %>% filter(abs(log2FoldChange) >1, padj<0.001) # sig DE genes

t24h_list = degdata[["t24h_Veh_&_t24h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t24h_list_degs = t24h_list %>% filter(abs(log2FoldChange) >1, padj<0.001) # sig DE genes
```

Find Significant Up-regulated or Down-regulated genes
=====================================================

Input parameters for the minimum log base 2 fold change to mark significant differential expression

``` r
## find upreg or downreg
t1h_upreg = t1h_list_degs %>% filter(log2FoldChange > 1)
t1h_downreg = t1h_list_degs %>% filter(log2FoldChange < 1)

t4h_upreg = t4h_list_degs %>% filter(log2FoldChange > 1)
t4h_downreg = t4h_list_degs %>% filter(log2FoldChange < 1)

t24h_upreg = t24h_list_degs %>% filter(log2FoldChange > 1)
t24h_downreg = t24h_list_degs %>% filter(log2FoldChange < 1)
```

Intersect sets to find overlaps
===============================

Intersect our data's three sets for upregulation, downregulation

``` r
## arrange sets
upreg_intersect_1hr_4hr <- intersect(t1h_upreg$gene, t4h_upreg$gene) # A&B 23
upreg_intersect_1hr_24hr <- intersect(t1h_upreg$gene, t24h_upreg$gene) # A&C 38
upreg_intersect_4hr_24hr <- intersect(t4h_upreg$gene, t24h_upreg$gene) # B&C 28

downreg_intersect_1hr_4hr <- intersect(t1h_downreg$gene, t4h_downreg$gene) # A&B 6
downreg_intersect_1hr_24hr <- intersect(t1h_downreg$gene, t24h_downreg$gene) # A&C 11
downreg_intersect_4hr_24hr <- intersect(t4h_downreg$gene, t24h_downreg$gene) # B&C 6

total_upreg_intersect <- intersect(t24h_upreg$gene, upreg_intersect_1hr_4hr) # A&B&C 15
total_downreg_intersect <- intersect(t24h_downreg$gene, downreg_intersect_1hr_4hr) # A&B&C 3
```

Plot the venn diagrams
======================

Plot up-regulated genes, down-regulated genes, and differentially expressed genes

``` r
## Plot venn diagrams
# upreg genes
grid.newpage()
draw.triple.venn(area1 = 192, area2 = 180, area3 = 96, n12 = 23, n23 = 28, n13 = 38, 
                 n123 = 15, category = c("1 hr Up", "4 hr Up", "24 hr Up"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))
```

![](DE_Venn_git_files/figure-markdown_github/unnamed-chunk-4-1.png)

    ## (polygon[GRID.polygon.11], polygon[GRID.polygon.12], polygon[GRID.polygon.13], polygon[GRID.polygon.14], polygon[GRID.polygon.15], polygon[GRID.polygon.16], text[GRID.text.17], text[GRID.text.18], text[GRID.text.19], text[GRID.text.20], text[GRID.text.21], text[GRID.text.22], text[GRID.text.23], text[GRID.text.24], text[GRID.text.25], text[GRID.text.26])

``` r
# downreg genes
grid.newpage()
draw.triple.venn(area1 = 240, area2 = 134, area3 = 39, n12 = 6, n23 = 6, n13 = 11, 
                 n123 = 3, category = c("1 hr Down", "4 hr Down", "24 hr Down"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))
```

![](DE_Venn_git_files/figure-markdown_github/unnamed-chunk-4-2.png)

    ## (polygon[GRID.polygon.27], polygon[GRID.polygon.28], polygon[GRID.polygon.29], polygon[GRID.polygon.30], polygon[GRID.polygon.31], polygon[GRID.polygon.32], text[GRID.text.33], text[GRID.text.34], text[GRID.text.35], text[GRID.text.36], text[GRID.text.37], text[GRID.text.38], text[GRID.text.39], text[GRID.text.40], text[GRID.text.41], text[GRID.text.42])

``` r
#overall genes
grid.newpage()
draw.triple.venn(area1 = 192+240, area2 = 180+134, area3 = 96+39, n12 = 23+6, n23 = 28+6, n13 = 38+11, 
                 n123 = 15+3, category = c("1 hr", "4 hr", "24 hr"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))
```

![](DE_Venn_git_files/figure-markdown_github/unnamed-chunk-4-3.png)

    ## (polygon[GRID.polygon.43], polygon[GRID.polygon.44], polygon[GRID.polygon.45], polygon[GRID.polygon.46], polygon[GRID.polygon.47], polygon[GRID.polygon.48], text[GRID.text.49], text[GRID.text.50], text[GRID.text.51], text[GRID.text.52], text[GRID.text.53], text[GRID.text.54], text[GRID.text.55], text[GRID.text.56], text[GRID.text.57], text[GRID.text.58])
