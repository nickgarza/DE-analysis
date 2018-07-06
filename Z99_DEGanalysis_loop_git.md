Differential Expression Analysis
================
Nick Garza
7/3/2018

Intialize: set working directory and load required packages
===========================================================

Notable packages: DESeq2. We use human annotated genes and DESeq2 functions to prepare for data generation

``` r
# set working directory to outputs
base.dir = "~/Documents/UVA/Z99"
outputs.dir = paste0(base.dir,"/outputs")
setwd(outputs.dir)

# load packages
library(DESeq2)
library(lattice)
library(bigWig)
library(dplyr)

# load functions
function.dir = paste0(base.dir,"/functions")
fn.name = "deseq_functions4.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "data_functions3.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "bedtools.R"
source( paste0(function.dir,"/",fn.name) )

# get gene annotations
gene.ann.dir = paste0(base.dir,"/annotation/human/")
fname = "Homo_sapiens.GRCh38.87.bed"
file = paste0(gene.ann.dir, fname)
segment = c("tss","end")
gene.file = get.gene.file(fname=file,genomic.interval=segment,removeM=T,removeX=F,removeY=F)

# import data
data.dir = paste0(base.dir,"/bigwigs/")
expr0 = get.raw.counts.interval(df=gene.file, path.to.bigWig=data.dir, file.prefix = 'Z99_')
size_factors = estimateSizeFactorsForMatrix(expr0)
```

Initialize DESeq object with size factors and conditions of interest
====================================================================

We initialize the deseq object with a count matrix, coldata, and design. We lastly generate a list of conditions for which we want to compare differential expression.

``` r
# create deseq object
conditions = sapply(names(size_factors),function(x){
  o1 = strsplit(x,"Z99_")[[1]][2]
  return(strsplit(o1,"_rep")[[1]][1])
}) %>% as.factor
colData = as.data.frame(conditions)
rownames(colData) = colnames( expr0 )
deseq_obj = DESeqDataSetFromMatrix(countData=expr0, colData=colData, design=~conditions)
sizeFactors(deseq_obj) = size_factors
check = counts(deseq_obj, normalized=TRUE) %>% as.data.frame
full_annotation = colData
rownames(full_annotation) = 1:nrow(full_annotation)

condits = list("t0", "t1h_Veh", "t1h_PDGF", "t4h_Veh", "t4h_PDGF", "t24h_Veh", "t24h_PDGF") %>%  as.character
full_annotation$conditions = full_annotation$conditions %>% as.character
```

Function for pairwise expression analysis
=========================================

This function generates a list of list with the first entry data\[\[1\]\] being the DESeq2 comparison matrices and the second entry data\[\[2\]\] a list of ma-plots to represent the data compared at the different time points.

``` r
# function to generate pairwise expression analysis and ma-plots of results
ma_deg <- function(condits, annotations){ 
  # iterator to keep track of organizing output data
  i <- 1 
  plots = list()
  degdata = list()
  
  # make non-redundant pairwise comparisons of DE genes based on user input
  for(ii in 1:(length(unique(annotations))-1)){ 
    for(jj in (ii+1):length(unique(annotations))){
      # filter redundant annotations
      indices1 = which(annotations == unique(annotations)[ii])
      indices2 = which(annotations == unique(annotations)[jj])
      
      indices = c(indices1,indices2)
      expr = expr0[,indices]
      
      conditions = factor(c(rep(condits[ii],3), rep(condits[jj],3)), levels = c(condits[ii],condits[jj]))
      colData = as.data.frame(conditions)
      rownames(colData) = colnames(expr)
      deseq_obj = DESeqDataSetFromMatrix(countData=expr, colData=colData, design=~conditions)
      
      sizeFactors(deseq_obj) = size_factors[indices]
      
      # perform DEG analysis and plot MA
      dds = DESeq(deseq_obj,quiet=T)
      res = results(dds)
      res = res[order(res$padj),] %>% as.data.frame
      iteration = paste(paste(condits[ii],"&", sep = "_"), condits[jj], sep = "_")
      degdata[[iteration]] = res
      plots[[i]] = plot.ma.lattice(ma.df=res, filename=NULL, p=0.001, title.main= iteration,log2fold=1)
      i = i + 1
    }
  }
  data = list(degdata,plots)
  return(data) # return list with degdata compared & plots of comparisons
}
```

Call function and save data + plots
===================================

An example of calling this function and putting its output into two separate lists

``` r
# call function with given params
deg_output = ma_deg(condits, full_annotation$conditions)
```

    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode
    ## converting counts to integer mode

``` r
degdata = deg_output[[1]]
plots = deg_output[[2]]
```
