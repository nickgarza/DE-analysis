
# loop to compute all pairwise DEG analyses
deg.comparisons = function(data=NULL,conditions=NULL,hm.cols=NULL,ann.cols=NULL,color.map=NULL,
                           heatmap.file=NULL,ma.file=NULL,cellwidth=20,cellheight=20,
                           fc=1,pval=0.001,fcUn=0.25,pvalUn=0.5){
  
  # organize data and condits if col.map is specified
  if(is.null(color.map)==F){
    orgs = organize.data(data,conditions,color.map)
    data = orgs$data; conditions = orgs$conditions; condits = orgs$condits
  }
  
  # loop through all pairwise comparisons
  # pairwise_deg() will return a matrix of deg data (degAll), a matrix of data for significant genes (degSig)
  # and lists of up-, down-, and un-regulated genes (genesUp, genesDown, genesUnreg)
  compares = condits
  deg_counts = matrix(c(0),length(compares),length(compares)) %>% as.data.frame
  deg_genes = c()
  deg_lists = list()
  deg_data = list()
  names(deg_counts) = rownames(deg_counts) = compares
  for(ii in 1:(length(compares)-1)){
    for(jj in (ii+1):length(compares)){
      
      t0 = compares[ii]
      tt = compares[jj]
      entry = paste0(t0,"_(reference)_verses_",tt)
      #print(entry)
      
      degs0 = pairwise_deg(data, conditions$time, t0, tt, fc, pval, fcUn, pvalUn)
      deg_counts[ii,jj] = as.data.frame(degs0$degSig) %>% filter(log2FoldChange > 0) %>% nrow %>% +1 %>% log10
      deg_counts[jj,ii] = as.data.frame(degs0$degSig) %>% filter(log2FoldChange < 0) %>% nrow %>% +1 %>% log10
      genes0 = as.data.frame(degs0$degSig) %>% rownames
      genes0 = sapply(genes0,function(x){strsplit(x,"_")[[1]][2]})
      deg_genes = c(deg_genes,unname(genes0, force = FALSE)) %>% unique
      
      sublist = list(up=degs0$genesUp, down=degs0$genesDown, unchanged=degs0$genesUnreg)
      deg_lists[[entry]] = sublist
      deg_data[[entry]] = degs0$degAll
    }
  }
  
  # heatmap of DEG counts
  if(is.null(heatmap.file)==F){
    plot.heatmap(plt_matrix=deg_counts, cols = c("gray"), sample.ann=colnames(deg_counts), colormap=color.map, plot.file=heatmap.file,
                 cellwidth=cellwidth, cellheight=cellheight)
  }
  
  # generate MA plots
  if(is.null(ma.file)==F){
    ma_plots = ma.plot.list(plt.data=deg_data, plt.file=ma.file)
  }
  
  out = list(deg.lists=deg_lists, deg.data=deg_data, deg.count.matrix=deg_counts, ma.plots=ma_plots)
  return(out)
  
} # deg.comparisons

#############################################################################################
## function to implement DESeq for a pairwise comparison
#############################################################################################
pairwise_deg = function(mat=NULL, condits=NULL, ref=NULL, comp=NULL, fcCut=NULL, fdrCut=NULL, fcUn = NULL, pvalUn=NULL){
  if(nrow(mat) > length(unique(rownames(mat)))){
    rownames(mat) = get.multiples(rownames(mat))
  }
  if(is.integer(mat[,1])==F){
    integer.or.not = apply(mat,2,function(y){
      if (!isTRUE(all(y == floor(y)))) stop("count data must only contain integer values")
    })
    rownamen = rownames(mat)
    mat = apply(mat,2,as.integer)
    rownames(mat) = rownamen
  }
  ind_t0 = which(condits == ref)
  ind_tt = which(condits == comp)
  deseq_dat = mat[rowSums(mat)>0, c(ind_t0,ind_tt)]
  sample.conditions = factor(c(rep("ref",length(ind_t0)), rep("comp",length(ind_tt))),levels=c("ref","comp"))
  deseq.counts.table = DESeqDataSetFromMatrix(deseq_dat, DataFrame(sample.conditions), ~ sample.conditions)
  colData(deseq.counts.table)$condition <- factor(colData(deseq.counts.table)$sample.conditions)
  cat(paste0('comparing ',ref,' (reference) to ',comp, ' (comparison) using DESeq\n'))
  dds = DESeq(deseq.counts.table,quiet=T)
  res = results(dds)
  res = res[order(res$padj),]
  ind_keep = intersect(which(abs(res$log2FoldChange) > fcCut), which(res$padj < fdrCut))
  ind_up = intersect( which(res$log2FoldChange > fcCut), which(res$padj < fdrCut) )
  up_list = sapply(rownames(res)[ind_up],function(x){strsplit(x,"_")[[1]][2]})
  up_list = sapply(up_list,function(x){strsplit(x,"_")[[1]][1]})
  ind_down = intersect( which(res$log2FoldChange < -fcCut), which(res$padj < fdrCut) )
  down_list = sapply(rownames(res)[ind_down],function(x){strsplit(x,"_")[[1]][2]})
  down_list = sapply(down_list,function(x){strsplit(x,"_")[[1]][1]})
  ind_unchanged = intersect( which(abs(res$log2FoldChange) > fcUn), which(res$padj > pvalUn) )
  un_list = sapply(rownames(res)[ind_unchanged],function(x){strsplit(x,"_")[[1]][2]})
  un_list = sapply(un_list,function(x){strsplit(x,"_")[[1]][1]})
  out = list(degAll=res, degSig=res[ind_keep,], genesUp=unique(unname(up_list)), genesDown=unique(unname(down_list)), genesUnreg=unique(unname(un_list)))
  return( out )
}

#############################################################################################
## function for basic MA plot
#############################################################################################
plot.ma.lattice <- function(ma.df=NULL, filename=NULL, p=NULL, title.main="Differential PRO Expression",log2fold=NULL) {
  y = ma.df$log2FoldChange
  x = log(ma.df$baseMean, base=10)
  dat = as.data.frame(cbind(x,y))
  names(dat)=c("x","y")
  out = xyplot(y ~ x, data=dat,
               groups=(ma.df$padj < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)),
               col=c("black","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"PRO fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
               panel=function(...) {
                 panel.xyplot(...)
                 panel.abline(h=0,col="gray")
               },
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2)));
 
  if(is.null(filename)==F){
    pdf(filename, width=3.83, height=3.83);
    #png(filename, width=3.83, height=3.83);
    print(out)
    dev.off()
  }
  return(out)
}



#############################################################################################
## generate file for a list of MA plots
#############################################################################################
library(gridExtra)
ma.plot.list = function(plt.data=NULL, plt.file=NULL, pval=NULL, fc=NULL){
  
  # generate a list of plots
  deg_data = plt.data
  ma_plots = list()
  for(ii in 1:length(deg_data)){
    namen0 = unlist( strsplit(names(deg_data)[ii],"_") )
    namen = c()
    for(jj in 1:length(namen0)){namen = paste(namen,namen0[jj])}
    ma_plots[[ii]] = plot.ma.lattice(ma.df=deg_data[[ii]],p=pval,title.main=namen,log2fold=fc)
  }
  
  pdf(plt.file)
  plt.per.page = 4
  for(plt in 1:(ceiling(length(ma_plots)/plt.per.page))){
    inds = c( ((plt-1)*plt.per.page+1):(plt*plt.per.page) )
    newlist=list()
    for (hh in 1:length(inds)){
      if(inds[hh]<=length(ma_plots)){
        newlist[[hh]]=ma_plots[[inds[hh]]] #[[1]]
      }
    }
    do.call(grid.arrange, c(newlist, list(ncol=2,nrow=2)))
  }
  dev.off()
  
  return(ma_plots)
}


#############################################################################################
## function for DEG analysis based on timing
## this function adjusts g and polII speed
## raw gene annotation and bigwig files are input
#############################################################################################

# gene file - gene annotation for the entire gene body
# bigwig.directory - directory for bigwig files
# file.prefix - prefix for bigwig data files
# condition.order = vector specifying the temporal order of the conditions
# time0.condition - designator for the zero timepoint condition (e.g., "t0")
# timing - matrix: col1 - conditions indicated by file names, col2 - corresponding times in minutes
# polspeed - speed of RNApolII in kb/min (default 2 kb/min)
# genomic.interval - c("tss","end"), c("tss", int), or c(int, "end") where int is an integer corresponding to 
#   the number of basepairs in the genomic interval of interest; e.g., c("tss",100) for pausing or c(200,"end")
#   for the gene body
# cut - genes with less than this number of base pairs will be removed from analysis
# fcCut - fold change cutoff for designating differential expression (default = 1)
# fdrCut - FDR cutoff for designating differential expression (default = 0.001)
# fcUn - fold change cutoff for designating absence of differential expression (default = 0.25)
# pvalUn - FDR cutoff for designating the absence of differential expression (default = 0.5)
# ma.file - file name for MA plots; without specifying this filename, no MA plots will be produced
# removeM - TRUE is mitochondrial chromosome data should be removed, FALSE otherwise
# removeX - TRUE is X chromosome data should be removed, FALSE otherwise
# removeY - TRUE is Y chromosome data should be removed, FALSE otherwise
# removeAuto - TRUE is autosomal chromosomes should be removed, FALSE otherwise
# note. this is the assumed file naming convention: general_condition_repNumber_strand.bigWig
#   only the 'condition' and 'repNumber' information is used for analysis
#   the 'general' and 'condition' indicators can be arbitrary
#   'repNumber' must be in the form 'rep1', 'rep2', etc

pairwise.deg.timing = function(gene.file=NULL, bigwig.directory=NULL, file.prefix=NULL, file.suffix=NULL, time0.condition=NULL, timing=NULL, 
                               polspeed=2, genomic.interval=NULL, cut=200, removeM=T, removeX=F, removeY=F, removeAuto=F,
                               condition.order=NULL, fcCut=1, fdrCut=0.001, fcUn=0.25, pvalUn=0.5, ma.file=NULL){
  
  # bigwig.directory = "/media/wa3j/Seagate/Documents/PRO/adipogenesis/bigWig"
  # file.prefix = 'adi_'
  # time0.condition = 't0'
  # genomic.interval = c(200,"end")
  # polspeed = 2
  # cut=200
  # removeM=T
  # removeX=F
  # removeY=F
  # removeAuto=F
  # fcCut=1
  # fdrCut=0.001
  # fcUn=0.25
  # pvalUn=0.5
  # ma.file = "adipogen_ma.pdf"
  
  deg.data = list()
  ma.data = list()
  
  file.dir.names = Sys.glob(file.path(bigwig.directory, paste(file.prefix,file.suffix,sep ='')))
  prefix.condition.rep = unname( sapply(file.dir.names,function(x){
    strsplit(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])], '_plus')[[1]][1]
  })) # these identifiers should correspond to the first column of the 'timing' matrix
  
  prefix1 = strsplit(prefix.condition.rep[1],"_")[[1]][1]
  prefix2 = "rep"
  conditions = get.conditions(prefix.condition.rep,prefix1,prefix2)
  file.data = cbind(file.dir.names,conditions)
  unique.conditions = unique(conditions[,2])
  if(all(condition.order %in% unique.conditions)==F) {warning("the conditions supplied in condition.order do not match the conditions identified based on file names")}
  if(nrow(timing) != nrow(file.data)){warning("the 'timing' data do not match the identified files")}
  
  # loop through all pairwise conditions: get genomic intervals and perform DE analysis
  for(ii in 1:(length(condition.order)-1)){
    file.cond1 = as.character(file.data[file.data[,3] == condition.order[ii], 1])
    all.rep.cond1 = as.character( file.data[which(file.data[,3]==condition.order[ii]), 2] )
    rep.cond1 = as.character( file.data[which(file.data[,3]==condition.order[ii])[1], 2] )
    min.cond1 = timing[timing[,1]==rep.cond1, 2]
    for(jj in (ii+1):length(condition.order)){
      file.cond2 = as.character(file.data[file.data[,3] == condition.order[jj], 1])
      all.rep.cond2 = as.character( file.data[which(file.data[,3]==condition.order[jj]), 2] )
      rep.cond2 = as.character( file.data[which(file.data[,3]==condition.order[jj])[1], 2] )
      min.cond2 = timing[timing[,1]==rep.cond2, 2]
      if(condition.order[ii]==time0.condition){min.cond1 = min.cond2} # set the time0 condition to the comparison time
      
      # gene annotations, deal with annotations for differening gene lengths due to timing
      gene1 = get.gene.timing(data=gene.file, time.min=min.cond1, polspeed=polspeed, genomic.interval=genomic.interval,
                               cut=cut,removeM=removeM,removeX=removeX,removeY=removeY,removeAuto=removeAuto)
      gene2 = get.gene.timing(data=gene.file, time.min=min.cond2, polspeed=polspeed, genomic.interval=genomic.interval,
                              cut=cut,removeM=removeM,removeX=removeX,removeY=removeY,removeAuto=removeAuto)
      gene.chr.str.cond1 = paste0(gene1[,1], ':', gene1[,2], '-', gene1[,3], '_', gene1[,4]) # gene annotation: chr:start-end_gene (1)
      gene.chr.str.cond2 = paste0(gene2[,1], ':', gene2[,2], '-', gene2[,3], '_', gene2[,4]) # gene annotation: chr:start-end_gene (2)
      gene.ann.12 = gene.chr.str.cond1
      if(length(gene.chr.str.cond1) != length(gene.chr.str.cond2) || all(gene.chr.str.cond1 == gene.chr.str.cond2)==F){
        ind.pos.cond1 = which(gene1[,6]=="+")
        ind.neg.cond1 = which(gene1[,6]=="-")
        ind.pos.cond2 = which(gene2[,6]=="+")
        ind.neg.cond2 = which(gene2[,6]=="-")
        ann.cond1 = rep("a",nrow(gene1))
        ann.cond2 = rep("a",nrow(gene2))
        ann.cond1[ind.pos.cond1] = paste0(gene1[ind.pos.cond1,1], ':', gene1[ind.pos.cond1,2], '_', gene1[ind.pos.cond1,4]) # chr:start_gene (1+)
        ann.cond2[ind.pos.cond2] = paste0(gene2[ind.pos.cond2,1], ':', gene2[ind.pos.cond2,2], '_', gene2[ind.pos.cond2,4]) # chr:start_gene (2+)
        ann.cond1[ind.neg.cond1] = paste0(gene1[ind.neg.cond1,1], ':', gene1[ind.neg.cond1,3], '_', gene1[ind.neg.cond1,4]) # chr:start_gene (1-)
        ann.cond2[ind.neg.cond2] = paste0(gene2[ind.neg.cond2,1], ':', gene2[ind.neg.cond2,3], '_', gene2[ind.neg.cond2,4]) # chr:start_gene (2-)
        if(length(ann.cond1) == length(ann.cond2)){
          if(all(ann.cond1 == ann.cond2)){
            gene.ann.12 = ann.cond1
          }
        } 
      }
      
      # get read data for ii and jj conditions
      expr.matrix.cond1 = c()
      for(iii in 1:length(file.cond1)){
        loaded.bw.plus = load.bigWig(file.cond1[iii])
        loaded.bw.minus = load.bigWig(gsub("plus", "minus", file.cond1[iii]))
        expr.matrix = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, gene1)
        expr.matrix.cond1 = cbind(expr.matrix.cond1,expr.matrix)
      }
      expr.matrix.cond2 = c()
      for(iii in 1:length(file.cond2)){
        loaded.bw.plus = load.bigWig(file.cond2[iii])
        loaded.bw.minus = load.bigWig(gsub("plus", "minus", file.cond2[iii]))
        expr.matrix = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, gene2)
        expr.matrix.cond2 = cbind(expr.matrix.cond2,expr.matrix)
      }
      expr.matrix.cond12 = cbind(expr.matrix.cond1, expr.matrix.cond2)
      colnames(expr.matrix.cond12) = c(all.rep.cond1, all.rep.cond2)
      rownames(expr.matrix.cond12) = gene.ann.12
      
      # perform differential expression analysis
      condit.id = paste0(condition.order[ii],"_vs_",condition.order[jj])
      condits = c(rep(condition.order[ii],length(all.rep.cond1)), rep(condition.order[jj],length(all.rep.cond2)))
      degs = pairwise_deg(mat=expr.matrix.cond12, condits=condits, ref=condition.order[ii], comp=condition.order[jj], 
                   fcCut=fcCut, fdrCut=fdrCut, fcUn=fcUn, pvalUn=pvalUn)
      deg.data[[condit.id]] = degs
      ma.data[[condit.id]] = degs$degAll
      
    } ## jj
  } ## ii
  
  # generate MA plots
  if(is.null(ma.file)==F){
    ma_plots = ma.plot.list(plt.data=ma.data, plt.file=ma.file, pval=fdrCut, fc=fcCut) 
  }
return(deg.data)
} # pairwise.deg.timing

# get multiple of items in a vector
# for duplicated items, annotate with _1, _2, etc
get.multiples = function(vec=NULL){
  out = vec
  vecU = unique(vec)
  for(ii in 1:length(vecU)){
    ind = which(vec == vecU[ii])
    if(length(ind)>1){
      for(jj in 1:length(ind)){out[ind[jj]] = paste0(out[ind[jj]],"_",jj)}
    }
  }
  return(out)
}

# get a matrix of conditions (col1) and replicate numbers (col2)
get.conditions = function(prefix.condition.rep=NULL, prefix1=NULL, prefix2=NULL){
  pf1 = paste0(prefix1,"_")
  pf2 = paste0("_",prefix2)
  base_condits = sapply( prefix.condition.rep, function(x) { strsplit(x, pf1)[[1]][2] })
  time_condits = sapply( base_condits, function(x) { strsplit(x, pf2)[[1]][1] })
  rep_condits = sapply( base_condits, function(x) { strsplit(x, pf2)[[1]][2] })
  conditions = cbind(base_condits,time_condits,rep_condits) %>% as.data.frame
  names(conditions) = c("condit","time","rep"); rownames(conditions) = c(1:nrow(conditions))
  conditions$rep = as.numeric(conditions$rep)
  return(conditions)
}



#### get all DEGs from the output of pairwise.deg.timing
# note that genesUp and genesDown positions in list are based on pairwise_deg()
get.all.degs = function(deg=NULL){
  alldeg = c()
  allup = c()
  alldown = c()
  for(ii in 1:length(deg)){
    degdat = c(deg[[ii]][[3]],deg[[ii]][[4]]) # genesUp and genesDown
    degup = deg[[ii]][[3]] # genesUp
    degdown = deg[[ii]][[4]] # genesDown
    alldeg = unique(c(alldeg,degdat))
    allup = unique(c(allup,degup))
    alldown = unique(c(alldown,degdown))
  }
  out = list(alldeg=alldeg, allup=allup, alldown=alldown)
  return(out)
}
