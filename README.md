library(GEOquery)
setwd('work_directory')

# Getting processed data from GEO into R (parse the Series Matrix files)
gse153922 <- getGEO('GSE153922')
gse153922 <- gse153922[[1]]   #  ExpressionSet object

################## Working with experimental metadata #################
library(tidyverse)
# gse153922$supplementary_file
pd153922 <- pData(gse153922)   # phenotype data
# pd153922[,c("geo_accession","title")]
pd153922['cel_file'] <- str_split(pd153922$supplementary_file,"/") %>% map_chr(tail,1)

################## Reading CEL data ##################
library(oligo)
library(oligoClasses)
## download data: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153922&format=file
gse153922_celdata <- read.celfiles(paste0('GSE153922/',pd153922$cel_file),phenoData=phenoData(gse153922)) # they are in the correct order as the experimental data we have from getGEO

################## Microarray Data processing with RMA ##################
gse153922_eset <- rma(gse153922_celdata)
head(exprs(gse153922_eset))   # gene expression matrix

################## Identifying differentially expressed genes using linear models ##################
library(limma)
library(dplyr)

# create a design representing the different groups
pd153922$group <- gsub("1|2|3| - | ","",pd153922$title)
pd153922$group <- as.factor(gsub("-",".",pd153922$group))
design153922 <- model.matrix(~ 0 + pd153922$group)
colnames(design153922) <- levels(pd153922$group)
design153922

# specify contrasts that might be interesting
contrasts_matrix153922 <- makeContrasts(ICMV.TRF_vs_iCMV.Vector=ICMV.TRF - iCMV.Vector,
                                        Young_vs_Old=Young - Old,
                                        levels=design153922)

contrasts_matrix153922

# Specify model for differential gene expression analysis
pData(gse153922_eset) <- pd153922
gse153922_fit <- lmFit(gse153922_eset,design153922)
gse153922_fit2 <- contrasts.fit(gse153922_fit,contrasts=contrasts_matrix153922)
gse153922_fit2 <- eBayes(gse153922_fit2)
summary(decideTests(gse153922_fit2))

# annotation
print(gse153922@annotation)
GPL16686 <-getGEO('GPL16686',destdir =".")
GPL16686_anno <- Table(GPL16686)
head(GPL16686_anno)

################## Basic downstream analysis of microarray data ##################
library(RColorBrewer)
library(hgu133plus2.db)
FC = 4
PValue = 0.05
for(contrast in c("ICMV.TRF_vs_iCMV.Vector","Young_vs_Old")){
  message(contrast)
  fitted.ebayes <- gse153922_fit2[,contrast]
  
  data = topTable(gse153922_fit2[,contrast],coef=1, p.value=0.05, lfc=0,n=Inf)
  data$sig <- NA
  data$sig[(data$adj.P.Val > PValue | data$adj.P.Val == "NA") | (data$logFC < log2(FC) & data$logFC > -log2(FC))] <- "NotSig"
  data$sig[(data$adj.P.Val <= PValue & data$logFC >= log2(FC))] <- "Up"
  data$sig[(data$adj.P.Val <= PValue & data$logFC <= -log2(FC))] <- "Down"
  table(data$sig)
  nrDEG = na.omit(data) 
  ps <- rownames(nrDEG)
  assign(contrast,ps)
  
  id_anno <- AnnotationDbi::select(hgu133plus2.db,GPL16686_anno$GB_ACC[match(ps,GPL16686_anno$ID)],c("SYMBOL","ENTREZID","GENENAME"),keytype="REFSEQ")
  write.csv(cbind(nrDEG,id_anno),paste0("limma_",contrast,"_p0.05_results.csv"),quote = F)
  
  interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=2)
  write.csv(interesting_genes,paste0("limma_",contrast,"_p0.05_lfc2_results.csv"),quote = F)
  eset_of_interest <- gse153922_eset[rownames(interesting_genes),which(pd153922$group==str_split(contrast,"_vs_")[[1]][1]|pd153922$group==str_split(contrast,"_vs_")[[1]][2])]
  
  # Volcano plots
  data$label = "" 
  library(ggrepel)  
  library(ggplot2)
  
  pdf(paste0(contrast,".volcanoplot.pdf"),width = 4,height = 4)
  print(ggplot(data,aes(data$logFC,-1*log10(data$adj.P.Val))) +    
          geom_point(aes(color = sig)) +                           
          labs( x="log[2](FC)", y="-log[10](PValue)",title = sprintf("%d features pass our cutoffs",nrow(interesting_genes))) +  
          xlim(-5,5) +
          scale_color_manual(values = c("blue","black","red")) + 
          geom_hline(yintercept=-log10(PValue),linetype=2)+        
          geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2)+ 
          geom_text_repel(aes(x = data$logFC,                   
                              y = -1*log10(data$adj.P.Val),          
                              label=label),                       
                          max.overlaps = 10000,                    
                          size=3,                                 
                          box.padding=unit(0.5,'lines'),           
                          point.padding=unit(0.1, 'lines'), 
                          segment.color='black',                   
                          show.legend=FALSE) +
          theme_bw())
  dev.off()
  
  pdf(paste0("limma_",contrast,"_p0.05_results.pdf"),width = 5)
  # Heatmaps
  print(heatmap(exprs(eset_of_interest)))
  print(heatmap(exprs(eset_of_interest),
                labCol=eset_of_interest[['title']] ,labRow=NA,
                col = rev(brewer.pal(10, "RdBu")),
                distfun   = function(x) as.dist(1-cor(t(x)))))
  dev.off()
}

length(intersect(Young_vs_Old,ICMV.TRF_vs_iCMV.Vector))

eset_of_interest <- gse153922_eset[intersect(Young_vs_Old,ICMV.TRF_vs_iCMV.Vector),]
pdf("heatmap.pdf",width = 4,height = 6)
print(heatmap(exprs(eset_of_interest),
              labCol=eset_of_interest[['title']] ,labRow=NA,
              col = rev(brewer.pal(10, "RdBu")),
              distfun   = function(x) as.dist(1-cor(t(x)))))
dev.off()
