library(tidyverse)
library(caret)
library(umap)
library(MOVICS)

flist <- c('DWI.csv',
           'DCE0.csv',
           'DCE2.csv',
           'T2WI.csv')

clust.data <- list()
for(fname in flist)
{
  dt <- read_csv(fname) %>% select(-ID)
  
  s.pre <- preProcess(dt, method = 'medianImpute')
  dt <- predict(s.pre, dt)
  
  cat(fname)
  print(which(is.na(dt)))
  
  dt.umap <- umap(dt, random_state = 123, n_neighbors = 10, n_components = 10)
  d.mat <- dt.umap$layout %>% t()
  
  colnames(d.mat) <- paste0('subj', 1:ncol(d.mat))
  rownames(d.mat) <- paste0('feat', 1:nrow(d.mat))
  
  clust.data[[str_replace(fname, '.csv', '')]] <- d.mat
}



optk.cluster <- getClustNum(data = clust.data,
                            is.binary = c(F, F, F, F),
                            try.N.clust = 2:8,
                            fig.name = 'Number of Cluster')

moic.res.list <- getMOIC(data        = clust.data,
                         methodslist = list("iClusterBayes", "SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian",  "gaussian", "gaussian"))


cmoic.clust <- getConsensusMOIC(moic.res.list = moic.res.list,
                                fig.name      = "CONSENSUS HEATMAP",
                                distance      = "euclidean",
                                linkage       = "average")

save(cmoic.clust, file = 'clustering.RData')
grp.info <- cmoic.clust$clust.res
write_csv(grp.info, 'grp_info.csv')

plotdata <- getStdiz(data       = clust.data,
                     halfwidth  = c(2,2,2, 2), # no truncation for mutation
                     centerFlag = c(T,T,T, T), # no center for mutation
                     scaleFlag  = c(T,T,T, T)) # no scale for mutation

dwi.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
dce0.col <- c("#6699CC", "white"  , "#FF3C38")
dce2.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
t2wi.col    <- c("#6699CC", "white"  , "#FF3C38")
col.list   <- list(dwi.col, dce0.col, dce2.col, t2wi.col)


feat   <- cmoic.clust$feat.res
feat1  <- feat[which(feat$dataset == "DWI"),][1:10,"feature"] 
feat2 <- feat[which(feat$dataset == "DCE0"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "DCE2"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "T2WI"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c('DWI', 'DCE0', 'DCE2', 'T2WI'),
             is.binary     = c(F,F,F, F), # the 4th data is mutation which is binary
             legend.name   = c('DWI', 'DCE0', 'DCE2', 'T2WI'),
             clust.res     = cmoic.clust$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F, F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")
