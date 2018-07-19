# Investigating functional gene expression at different time points
library(DESeq2)
library(lattice)
library(bigWig)
library(dplyr)
library(enrichR)

# set working directory to outputs
base.dir = "~/Documents/UVA/Z99"
outputs.dir = paste0(base.dir,"/outputs")
setwd(outputs.dir)

# load functions and data
function.dir = paste0(base.dir,"/functions")
fn.name = "deseq_functions4.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "data_functions3.R"
source( paste0(function.dir,"/",fn.name) )
fn.name = "bedtools.R"
source( paste0(function.dir,"/",fn.name) )

# data
load("pdgf_veh_t1hgenes.RData")
load("pdgf_veh_t4hgenes.RData")
load("pdgf_veh_t24hgenes.RData")

# functional expression analysis through enrichR
dbs <- listEnrichrDbs()
dbs_ind = c(93:95,99,102,115:117)
Db = dbs[dbs_ind,1] 
enriched_1hr <- enrichr(t1hgenes, Db) %>% as.list
enriched_4hr <- enrichr(t4hgenes, Db) %>% as.list 
enriched_24hr <- enrichr(t24hgenes, Db) %>% as.list 

head(enriched_1hr[[1]][,c(1,3)])
head(enriched_4hr[[1]][,c(1,3)])

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
#y <- fdr_filter(enriched_1hr[[1]])

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


# write.csv(reduced_df,'pdgf_veh_vehicle_func_expr.csv')
