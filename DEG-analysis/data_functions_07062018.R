#############################################################################################
## function to import bigWig data
## irrespective of time, all conditions treated identically
#############################################################################################

get.raw.counts.interval <- function(df, path.to.bigWig, file.prefix = 'H') {
  vec.names = c()
  inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
  
  for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*_plus_body_0-mer.bigWig", sep ='')))) {
    factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
    #print(factor.name)
    vec.names = c(vec.names, factor.name)
    loaded.bw.plus = load.bigWig(mod.bigWig)
    #print(mod.bigWig)
    minus.bigwig = paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep='')
    #print(minus.bigwig)
    loaded.bw.minus = load.bigWig(minus.bigwig)
    mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
    inten.df = cbind(inten.df, mod.inten)
  }
  colnames(inten.df) = vec.names
  r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
  row.names(inten.df) = r.names
  return(inten.df)
}


#############################################################################################
## function to generate a data matrix for PCA
## this function implements DESeq to generate a DESeqDataSet
#############################################################################################

run.deseq.table <- function(mat) {
  sample.conditions = factor(colnames(mat))
  if(ncol(mat) > length(unique(colnames(mat)))){colnames(mat) = get.multiples(colnames(mat))}
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions)
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=levels(sample.conditions))
  dds = DESeq(deseq.counts.table)
  return(dds)
}

#############################################################################################
## function to get gene annotations
#############################################################################################
get.gene.file = function(fname=NULL,data=NULL,genomic.interval=NULL,cut=200,removeM=FALSE,removeX=FALSE,removeY=FALSE,removeAuto=FALSE){
  
  if(genomic.interval[1]!="tss" && genomic.interval[2]!="end"){warning("set segment[1]=tss or segment[2]=end (or both for whole gene)")}
  if(is.null(fname)==F && is.null(data)==F){warning("enter fname or data, not both")}
  
  # import file and initial processing
  if(is.null(data)==F){gene.file0 = data}
  if(is.null(fname)==F){
    auto = paste0("chr",c(1:22))
    chr_keep = c( auto, "chrX", "chrY", "chrMT" )
    bed6 = read.table(fname,stringsAsFactors=F)
    bed6b = bed6[bed6[,1] %in% chr_keep,]
    gene.file = bed6b
    gene.file[bed6b$V1=="chrMT",1] = "chrM"
    gene.file0 = gene.file[!duplicated(paste(gene.file[,1],gene.file[,2],gene.file[,3])),]
  }
  
  # remove genes based on class
  if(removeM==T){gene.file0 = gene.file0[gene.file0[,1] != "chrM",]}
  if(removeX==T){gene.file0 = gene.file0[gene.file0[,1] != "chrX",]}
  if(removeY==T){gene.file0 = gene.file0[gene.file0[,1] != "chrY",]}
  if(removeAuto==T){gene.file0 = gene.file0[!(gene.file0[,1] %in% auto),]}

  # remove genes that are below a length threshold
  gene.file0 = gene.file0[gene.file0[,3] > gene.file0[,2], ]
  gene.file0 = gene.file0[gene.file0[,3] - gene.file0[,2] > cut, ]
  
  # TSS through END of the gene (whole gene)
  if(genomic.interval[1]=="tss" && genomic.interval[2]=="end"){output = gene.file0}
  
  # TSS through a set interval (< whole gene)
  # remove genes that do not contain the entire interval
  if(genomic.interval[1]=="tss" && genomic.interval[2]!="end"){
    
    terminal_base = as.numeric(genomic.interval[2])
    ind_plus = which(gene.file0[,6] == "+")
    ind_minus = which(gene.file0[,6] == "-")
    
    # remove genes without the relevant interval
    ind_short_plus = intersect( which((gene.file0[,3] <= (gene.file0[,2]+terminal_base))==T), ind_plus )
    ind_short_minus = intersect( which((gene.file0[,2] >= (gene.file0[,3]-terminal_base))==T), ind_minus )
    ind_short = c(ind_short_plus,ind_short_minus)
    if(length(ind_short)>0){gene.file0 = gene.file0[-ind_short,]}
    
    # keep only the interval of interest
    ind_plus = which(gene.file0[,6] == "+")
    ind_minus = which(gene.file0[,6] == "-")
    gene.file0[ind_plus,3] = gene.file0[ind_plus,2] + terminal_base
    gene.file0[ind_minus,2] = gene.file0[ind_minus,3] - terminal_base
    output = gene.file0
  }
  
  # a set interval from the gene END (< whole gene)
  # remove genes that do not contain the entire interval
  if(genomic.interval[1]!="tss" && genomic.interval[2]=="end"){
    
    initial_base = as.numeric(genomic.interval[1])
    ind_plus = which(gene.file0[,6] == "+")
    ind_minus = which(gene.file0[,6] == "-")
    
    # remove genes without the relevant interval
    ind_short_plus = intersect( which((gene.file0[,3] <= (gene.file0[,2]+initial_base))==T), ind_plus )
    ind_short_minus = intersect( which((gene.file0[,2] >= (gene.file0[,3]-initial_base))==T), ind_minus )
    ind_short = c(ind_short_plus,ind_short_minus)
    if(length(ind_short)>0){gene.file0 = gene.file0[-ind_short,]}
    
    # keep only the interval of interest
    ind_plus = which(gene.file0[,6] == "+")
    ind_minus = which(gene.file0[,6] == "-")
    gene.file0[ind_plus,2] = gene.file0[ind_plus,2] + initial_base
    gene.file0[ind_minus,3] = gene.file0[ind_minus,3] - initial_base
    output = gene.file0
  }
  
  output = output[!duplicated(paste(output[,1],output[,2],output[,3],output[,4])),]
  return(output)
} # get.gene.file


#############################################################################################
## function to get file names
#############################################################################################

get.file.names = function(gene.file=NULL, directory=NULL, file.prefix=NULL){
  files = Sys.glob(file.path(directory, paste(file.prefix, "*_plus.bigWig", sep ='')))
  factor_names = sapply(files,function(x){
    out = strsplit(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])], '_plus')[[1]][1]
  }) %>% unname
  return(factor_names)
}


#############################################################################################
## function to take in gene annotation, gene length of interest based on timing 
## and adjust the gene annotation accordingly
#############################################################################################
# data - the gene annotation file (whole gene annotation assumed)
# time.min - time (min) to offset gene length based on polII transcription speed
# polspeed - polII transcription speed (kb/min)
# genomic.interval - c("tss","end"), c("tss", int), or c(int, "end") where int is an integer 
# cut - genes with less than this number of base pairs will be removed from analysis
# removeM - TRUE is mitochondrial chromosome data should be removed, FALSE otherwise
# removeX - TRUE is X chromosome data should be removed, FALSE otherwise
# removeY - TRUE is Y chromosome data should be removed, FALSE otherwise
# removeAuto - TRUE is autosomal chromosomes should be removed, FALSE otherwise
get.gene.timing = function(data=NULL, time.min=NULL, polspeed=NULL, genomic.interval=NULL,
                           cut=NULL,removeM=FALSE,removeX=FALSE,removeY=FALSE,removeAuto=FALSE){

  dist.bp = polspeed*time.min*1000 # polspeed is in kb/min
  gene.file0 = data
  ind_plus = which(gene.file0[,6] == "+")
  ind_minus = which(gene.file0[,6] == "-")
  
  # adjust appropriate genes to the relevant interval
  # if dist.bp = 0, no adjustments are made
  if(dist.bp != 0){
    ind_plus_long = intersect( which((gene.file0[,3] >= (gene.file0[,2]+dist.bp))==T), ind_plus )
    ind_minus_long = intersect( which((gene.file0[,2] <= (gene.file0[,3]-dist.bp))==T), ind_minus )
    if(length(ind_plus_long)>0){
      gene.file0[ind_plus_long,3] = gene.file0[ind_plus_long,2] + dist.bp
    } 
    if(length(ind_minus_long)>0){
      gene.file0[ind_minus_long,2] = gene.file0[ind_minus_long,3] - dist.bp
    } 
  }
  
  # run the original filtering function to get appropriate genomic interval
  out = get.gene.file(fname=NULL,data=gene.file0,genomic.interval=genomic.interval,cut=cut,
                      removeM=removeM,removeX=removeX,removeY=removeY,removeAuto=removeAuto)
  if(nrow(out)==0){warning("there is no gene annotation output, check the cutoff parameter, etc...")}
  return(out)
} # get.gene.timing


#############################################################################################
## function to take in gene annotation and timing information and get an expression matrix
## only common genes are included
## this function can be used to generate data to be displayed in a heatmap
#############################################################################################
get.timing.matrix = function(timing=NULL, polspeed=NULL, genomic.interval=NULL, gene.file=NULL, file.prefix=NULL, bigwig.directory=NULL,
                             time0.condition=NULL, condition.order=NULL, cut=NULL, removeM=FALSE, removeX=FALSE, removeY=FALSE,removeAuto=FALSE){

  ##################################
  ## this entire section is being done by get.beds.timing()
  # time0.condition=c("t0",20)
  file.dir.names = Sys.glob(file.path(bigwig.directory, paste(file.prefix, "*_plus.bigWig", sep ='')))
  prefix.condition.rep = unname( sapply(file.dir.names,function(x){
    strsplit(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])], '_plus')[[1]][1]
  })) 
  
  # basic sample annotation
  prefix1 = strsplit(prefix.condition.rep[1],"_")[[1]][1]
  prefix2 = "rep"
  conditions = get.conditions(prefix.condition.rep,prefix1,prefix2)
  file.data = cbind(file.dir.names,conditions)
  unique.conditions = unique(conditions[,2])
  if(all(condition.order %in% unique.conditions)==F) {warning("the conditions supplied in condition.order do not match the conditions identified based on file names")}
  if(nrow(timing) != nrow(file.data)){warning("the 'timing' data do not match the identified files")}

  # get gene coordinate annotation for each condition
  gene.data = list()
  for(ii in 1:length(condition.order)){
    file.cond1 = as.character(file.data[file.data[,3] == condition.order[ii], 1])
    all.rep.cond1 = as.character( file.data[which(file.data[,3]==condition.order[ii]), 2] )
    rep.cond1 = as.character( file.data[which(file.data[,3]==condition.order[ii])[1], 2] )
    min.cond1 = timing[timing[,1]==rep.cond1, 2]
    cond1 = strsplit(rep.cond1,"_")[[1]][1]
    if(is.null(time0.condition)==F){
      if(cond1 == time0.condition[1]){min.cond1 = as.numeric(time0.condition[2])}
    }
    gene.data[[ condition.order[ii] ]] = get.gene.timing(data=gene.file, time.min=min.cond1, polspeed=polspeed, genomic.interval=genomic.interval,
                            cut=cut,removeM=removeM,removeX=removeX,removeY=removeY,removeAuto=removeAuto)
  } # ii
  ##################################
  ## this entire section is being done by get.beds.timing()
  
  # get common  transcripts for analysis
  tr = list()
  for(ii in 1:length(gene.data)){
    gene1 = gene.data[[ii]]
    ind.pos.cond1 = which(gene1[,6]=="+")
    ind.neg.cond1 = which(gene1[,6]=="-")
    ann.cond1 = rep("a",nrow(gene1))
    ann.cond1[ind.pos.cond1] = paste0(gene1[ind.pos.cond1,1], ':', gene1[ind.pos.cond1,2], '_', gene1[ind.pos.cond1,4]) # chr:start_gene (1+)
    ann.cond1[ind.neg.cond1] = paste0(gene1[ind.neg.cond1,1], ':', gene1[ind.neg.cond1,3], '_', gene1[ind.neg.cond1,4]) # chr:start_gene (1-)
    if(length(ann.cond1) > length(unique(ann.cond1))){ann.cond1 = get.multiples(ann.cond1)}
    tr[[names(gene.data)[ii]]] = ann.cond1
  }
  common.genes = Reduce(intersect, tr)
  
  # get gene data for common transcripts
  gene.data.common = list()
  for(ii in 1:length(gene.data)){
    gene1 = gene.data[[ii]]
    ind.pos.cond1 = which(gene1[,6]=="+")
    ind.neg.cond1 = which(gene1[,6]=="-")
    ann.cond1 = rep("a",nrow(gene1))
    ann.cond1[ind.pos.cond1] = paste0(gene1[ind.pos.cond1,1], ':', gene1[ind.pos.cond1,2], '_', gene1[ind.pos.cond1,4]) # chr:start_gene (1+)
    ann.cond1[ind.neg.cond1] = paste0(gene1[ind.neg.cond1,1], ':', gene1[ind.neg.cond1,3], '_', gene1[ind.neg.cond1,4]) # chr:start_gene (1-)
    if(length(ann.cond1) > length(unique(ann.cond1))){ann.cond1 = get.multiples(ann.cond1)}
    indices = unlist(sapply(common.genes,function(x){which(ann.cond1==x)}))
    gene.data.common[[ condition.order[ii] ]] = gene1[indices,]
  } # ii
  
  # build expression matrix
  data.out = matrix(c(0),length(common.genes),nrow(conditions))
  rownames(data.out) = common.genes
  colnames(data.out) = conditions[,2]
  for(ii in 1:length(file.dir.names)){
    factor.name = strsplit(strsplit(file.dir.names[ii], "/")[[1]][length(strsplit(file.dir.names[ii], "/")[[1]])], '_plus')[[1]][1]
    loaded.bw.plus = load.bigWig(file.dir.names[ii])
    loaded.bw.minus = load.bigWig(paste0(bigwig.directory,'/',factor.name, '_minus.bigWig'))
    expr_ii = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, gene.data.common[[conditions[ii,2]]])
    data.out[,ii] = expr_ii
  }
  return(data.out)
} # get.timing.matrix



### function to generate heatmap plots as a function of time
# center.replicates - specify whether to center replicates, NULL (default), mean, or median
timing.heatmap = function(gene.file=NULL, bigwig.directory=NULL, file.prefix=NULL, time0.condition=NULL, timing=NULL, 
                          polspeed=2, genomic.interval=NULL, cut=200, removeM=TRUE, removeX=FALSE, removeY=FALSE, removeAuto=FALSE,
                          condition.order=NULL, hm.file=NULL, gene.list=NULL,z.score=TRUE,row.col.ann=FALSE,
                          cols=NULL,sample.ann=NULL,colors=NULL,colormap=NULL,cellwidth=NULL,cellheight=NULL,
                          center.replicates=NULL,peak.sort=TRUE){
  
  # get expression matrix and gene list
  mat = get.timing.matrix(gene.file=gene.file, bigwig.directory=bigwig.directory, file.prefix=file.prefix, time0.condition=time0.condition, 
                          timing=timing, polspeed=polspeed, genomic.interval=genomic.interval, cut=cut, 
                          removeM=removeM, removeX=removeX, removeY=removeY, removeAuto=removeAuto,
                          condition.order=condition.order)
  mat.genes = unname( sapply(rownames(mat),function(x){strsplit(x,"_")[[1]][[2]]}) )
  mat.genes = unname( sapply(mat.genes,function(x){strsplit(x,"_")[[1]][[1]]}) )
  
  # get genes of interest and matrix for heatmap plotting
  if(is.null(gene.list)==FALSE){
    inds = unlist(sapply(gene.list,function(x){which(mat.genes==x)}))
    if(length(inds)>0){
      mat.out = mat[inds,]
      genes.out = unname( sapply(rownames(mat.out),function(x){strsplit(x,"_")[[1]][[2]]}) )
      genes.out = unname( sapply(genes.out,function(x){strsplit(x,"_")[[1]][[1]]}) )
    }
  }
  
  # normalize/scale/average/clip data for heatmap
  deseq_obj = run.deseq.table(mat.out)
  log2_exp = rlogTransformation(deseq_obj)
  hmdat0 = assay(log2_exp)
  if(z.score==TRUE){hmdat0 = t(scale(t(hmdat0),center=TRUE,scale=TRUE))}
  if(is.null(center.replicates)==FALSE){hmdat0 = center.reps(hmdat0,center.replicates)}
  if(peak.sort==TRUE){hmdat0 = sort.peaks(hmdat0)}
  
  # produce heatmap
  condit = unique(sample.ann[,1])
  plot.heatmap(plt_matrix=hmdat0,sample.ann=condit,colors=colors,colormap=colormap,plot.file=hm.file,
               cellwidth=cellwidth,cellheight=cellheight,row.col.ann=row.col.ann)
  
  return(hmdat0)
} # timing.heatmap



# function to center multiple replicates
center.reps = function(mat=NULL,center.replicates=NULL){
  reps = sapply(colnames(mat),function(x){strsplit(x,"_")[[1]][1]})
  reps.unique = unique(reps)
  out = matrix(c(0),nrow(mat),length(reps.unique))
  rownames(out) = rownames(mat)
  colnames(out) = reps.unique
  for(ii in 1:length(reps.unique)){
    inds = which(reps == reps.unique[ii])
    if(center.replicates == "median"){ out[,ii] = apply(mat[,inds],1,median) }
    if(center.replicates == "mean"){ out[,ii] = apply(mat[,inds],1,mean) }
  }  
  return(out)
}

# function to sort a matrix based on peak level
sort.peaks = function(mat){
  max.inds = apply(mat,1,function(x){which(x==max(x))})
  sorted.data = sort(max.inds, decreasing=TRUE, index.return=TRUE)
  sorted.vals = sorted.data$x
  sorted.inds = sorted.data$ix
  mat.sorted = mat[sorted.inds,]
  out = matrix(c(0),nrow(mat),ncol(mat))
  rownames(out) = rownames(mat.sorted)
  colnames(out) = colnames(mat.sorted)
  for(ii in 1:length(unique(sorted.vals))){
    inds = which(sorted.vals == unique(sorted.vals)[ii])
    max.inds = sort(mat.sorted[inds,unique(sorted.vals)[ii]], index.return=TRUE)$ix
    out[inds,] = mat.sorted[inds[max.inds],]
  }
  return(out)
}
