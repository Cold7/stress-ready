library("sva") #Note this exercise requires sva (>= v3.36.0) which is only available for R (>= v4.x)
library("ggplot2")
library("gridExtra")
library("UpSetR")
library ( DESeq2 )
library(EnhancedVolcano)
library (NMF)

readcounts <- read.csv("counts_ortholog.csv", header = TRUE, row.names = 1)
head(readcounts,n=3)
orig_names <- names(readcounts)

#orig_names
condition = c("Sly.drought","Sly.drought","Sly.drought","Sly.control","Sly.control","Sly.control","Spenn.drought","Spenn.drought","Spenn.drought", "Spenn.control","Spenn.control","Spenn.control")

sample_info <- data.frame(condition=condition, row.names=orig_names)

DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = sample_info,
                                   design = ~ condition)
## remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 10 , ]
#DESeq.ds<- DESeq.ds[rowSums(counts(DESeq.ds)>=10)>0,]

DESeq.ds <- estimateSizeFactors (DESeq.ds)
sizeFactors(DESeq.ds)

## if you check colData () again , you see that this now contains the sizeFactors
#colData(DESeq.ds)

#counts() allows you to immediately retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE )
counts.sf_normalized
#transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

DESeq.rlog <- rlog(DESeq.ds, blind = FALSE ) 
#This function transforms the count data to the log2 scale in a way which minimizes 
# differences between samples for rows with small counts, and which normalizes with 
#respect to library size. The rlog transformation produces a similar variance 
#stabilizing effect 
rlog.norm.counts <- assay(DESeq.rlog)



###########################
##
## DESeq2
##
###########################

#Now, running the DGE analysis
DESeq.ds <- DESeq(DESeq.ds)


######################################
#
# Sly control vs drought 
#
######################################

DGE.results <- results(DESeq.ds, contrast=c("condition", "Sly.drought", "Sly.control"), independentFiltering=TRUE, alpha=0.01, pAdjustMethod="BH")
#saving data
res <- DGE.results[order(DGE.results$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
write.csv(resdata, file="Sly.drought_Sly.control.csv")

######################################
#
# Spenn control vs drought 
#
######################################

DGE.results <- results(DESeq.ds, contrast=c("condition", "Spenn.drought", "Spenn.control"), independentFiltering=TRUE, alpha=0.01, pAdjustMethod="BH")
#saving data
res <- DGE.results[order(DGE.results$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
write.csv(resdata, file="Spenn.drought_Spenn.control.csv")

######################################
#
# Sly control vs Spenn drought 
#
######################################

DGE.results <- results(DESeq.ds, contrast=c("condition", "Spenn.drought", "Sly.control"), independentFiltering=TRUE, alpha=0.01, pAdjustMethod="BH")
#saving data
res <- DGE.results[order(DGE.results$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
write.csv(resdata, file="Spenn.drought_Sly.control.csv")

######################################
#
# Sly drought vs Spenn control 
#
######################################

DGE.results <- results(DESeq.ds, contrast=c("condition", "Sly.drought", "Spenn.control"), independentFiltering=TRUE, alpha=0.01, pAdjustMethod="BH")
#saving data
res <- DGE.results[order(DGE.results$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
write.csv(resdata, file="Sly.drought_Spenn.control.csv")


######################################
#
# Spenn control vs Sly control
#
######################################

DGE.results <- results(DESeq.ds, contrast=c("condition", "Spenn.control", "Sly.control"), independentFiltering=TRUE, alpha=0.01, pAdjustMethod="BH")
#saving data
res <- DGE.results[order(DGE.results$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
write.csv(resdata, file="Spenn.control_Sly.control.csv")
