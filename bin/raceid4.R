#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Matrix)
library(data.table)
library(Seurat)

date = Sys.Date()
#download counts file from url: https://stackoverflow.com/questions/28986150/downloading-and-extracting-gz-data-file-using-r
#aut_counts dataset

#url <- "https://cells.ucsc.edu/autism/rawMatrix.zip"

#counts <- Read10X("data/")
genes <- data.frame("ID" = fread("data/exprMatrix.tsv")[["gene"]])
genes <- tidyr::separate(genes, ID, into = c("ID", "Gene"), sep = "\\|")
counts <- readMM("data/matrix.mtx")

meta <- read_delim("data/meta.txt", delim = "\t")
micr <- which(meta$cell_type == "Microglia" & meta$diagnosis == "Control") # & meta$diagnosis == "Control"
micr_counts <- Matrix(as.matrix(counts[,micr]), sparse = T)
dimnames(micr_counts) <- list(genes$Gene, meta$cell[micr])


save(micr_counts, file = "data/schirmer_ms_nuc_seq-microglia.Robj")

load("data/schirmer_ms_nuc_seq-microglia.Robj")

prdata <- micr_counts

sc <- SCseq(prdata)
# filtering of expression data
a <- apply(prdata, 2, sum)
sc <- filterdata(sc, mintotal=quantile(a, 0.1)) # exlcude the lower quartile of the cells
sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('JUN',
                           'FOS',
                           'ZFP36',
                           'HSPA1A|HSPA1B',
                           'DUSP1',
                           'EGR1',
                           'MALAT1'))

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc) 

plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)

sc <- clustexp(sc,cln=11,sat=FALSE) 
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
ord_clust <- clustheatmap(sc)
save(ord_clust, file = 'data/ord_clust.Robj')

pdf(paste0('plots/heatmaps/clustheatmap.pdf'))
clustheatmap(sc, final = T)
dev.off()

sc <- comptsne(sc)
sc <- compfr(sc,knn=10)

plotmap(sc)
plotmap(sc,fr=TRUE)
dev.off()

name2id <- function(x,id) {
  ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
  n <- c()
  for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
  }
  n
}

plotexpmap(sc,name2id("MRC1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("LYVE1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD163", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("TMEM119", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CX3CR1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("PTPRC", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD3E", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("ITGAM", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD8A", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("CD4", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("P2RY12", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("SLC2A5", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("^EGR1", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("JUN", rownames(sc@ndata)),logsc=F,fr=F)
plotexpmap(sc,name2id("GPR34", rownames(sc@ndata)),logsc=F,fr=F)

dg <- clustdiffgenes(sc,4,pvalue=.01)
head(dg,25)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

#Save sc file
save(sc, file = 'data/sc.Robj')

#write_csv(as.data.frame(micr_ids), "data/microglia-cell-ids.csv")
