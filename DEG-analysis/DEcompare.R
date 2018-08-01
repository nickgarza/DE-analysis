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