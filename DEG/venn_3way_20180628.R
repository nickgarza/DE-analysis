# This script analysis the DE between t0 and t1, t4, t24 respectively

# set working directory to outputs
base.dir = "~/Documents/UVA/Z99"
outputs.dir = paste0(base.dir,"/outputs")
setwd(outputs.dir)

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

# acquire significant DE genes
t1h_list = degdata[["t1h_Veh_&_t1h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t1h_list_degs = t1h_list %>% filter(abs(log2FoldChange)>1, padj<0.001) # sig DE genes

t4h_list = degdata[["t4h_Veh_&_t4h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t4h_list_degs = t4h_list %>% filter(abs(log2FoldChange) >1, padj<0.001) # sig DE genes

t24h_list = degdata[["t24h_Veh_&_t24h_PDGF"]] %>% as.data.frame %>% rownames_to_column('gene')
t24h_list_degs = t24h_list %>% filter(abs(log2FoldChange) >1, padj<0.001) # sig DE genes

# find upreg or downreg
t1h_upreg = t1h_list_degs %>% filter(log2FoldChange > 1)
t1h_downreg = t1h_list_degs %>% filter(log2FoldChange < 1)

t4h_upreg = t4h_list_degs %>% filter(log2FoldChange > 1)
t4h_downreg = t4h_list_degs %>% filter(log2FoldChange < 1)

t24h_upreg = t24h_list_degs %>% filter(log2FoldChange > 1)
t24h_downreg = t24h_list_degs %>% filter(log2FoldChange < 1)

# arrange sets
upreg_intersect_1hr_4hr <- intersect(t1h_upreg$gene, t4h_upreg$gene) # A&B 23
upreg_intersect_1hr_24hr <- intersect(t1h_upreg$gene, t24h_upreg$gene) # A&C 38
upreg_intersect_4hr_24hr <- intersect(t4h_upreg$gene, t24h_upreg$gene) # B&C 28

downreg_intersect_1hr_4hr <- intersect(t1h_downreg$gene, t4h_downreg$gene) # A&B 6
downreg_intersect_1hr_24hr <- intersect(t1h_downreg$gene, t24h_downreg$gene) # A&C 11
downreg_intersect_4hr_24hr <- intersect(t4h_downreg$gene, t24h_downreg$gene) # B&C 6

total_upreg_intersect <- intersect(t24h_upreg$gene, upreg_intersect_1hr_4hr) # A&B&C 15
total_downreg_intersect <- intersect(t24h_downreg$gene, downreg_intersect_1hr_4hr) # A&B&C 3


# plot venn diagrams
nrow(t1h_upreg) # 192
nrow(t4h_upreg) # 180
nrow(t24h_upreg) # 96

nrow(t1h_downreg) # 240
nrow(t4h_downreg) # 134
nrow(t24h_downreg) # 39

nrow(total_upreg_intersect) # 15
nrow(total_downreg_intersect) # 3

num_t1_upreg = nrow(t1h_upreg) %>% as.numeric
num_t4_upreg = nrow(t4h_upreg) %>% as.numeric
num_t1_downreg = nrow(t1h_downreg) %>% as.numeric
num_t4_downreg = nrow(t4h_downreg) %>% as.numeric
num_upreg_intersect = nrow(upreg_intersect) %>% as.numeric
num_downreg_intersect = nrow(downreg_intersect) %>% as.numeric


# upreg genes
grid.newpage()
draw.triple.venn(area1 = 192, area2 = 180, area3 = 96, n12 = 23, n23 = 28, n13 = 38, 
                 n123 = 15, category = c("1 hr Up", "4 hr Up", "24 hr Up"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))


# downreg genes
grid.newpage()
draw.triple.venn(area1 = 240, area2 = 134, area3 = 39, n12 = 6, n23 = 6, n13 = 11, 
                 n123 = 3, category = c("1 hr Down", "4 hr Down", "24 hr Down"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

#overall genes
grid.newpage()
draw.triple.venn(area1 = 192+240, area2 = 180+134, area3 = 96+39, n12 = 23+6, n23 = 28+6, n13 = 38+11, 
                 n123 = 15+3, category = c("1 hr", "4 hr", "24 hr"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))


t1hgenes = sapply(t1h_list_degs$gene,function(x){strsplit(x,"_")[[1]][2]})
t4hgenes = sapply(t4h_list_degs$gene,function(x){strsplit(x,"_")[[1]][2]})
t24hgenes = sapply(t24h_list_degs$gene,function(x){strsplit(x,"_")[[1]][2]})
save.image("pdgf_veh_t1hgenes.RData")
save.image("pdgf_veh_t4hgenes.RData")
save.image("pdgf_veh_t24hgenes.RData")
