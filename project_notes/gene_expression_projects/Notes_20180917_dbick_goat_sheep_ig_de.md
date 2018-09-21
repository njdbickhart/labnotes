# Goat and Sheep IG differential expression analysis
---
*9/17/2018*

These are my notes on running the EBseq pipeline on John H and John S's differential expression data on Goat and Sheep.


## Preparing the data

So, I first need to format the data into a useable matrix. I transposed the data in excel and copied into text files separated by shared tissue types. I also trimmed the gene comparisons down to the shared set of genes between both goat and sheep.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/side_projects/hammond_goat_sheep_ig

```bash
for i in *.tab; do echo $i; dos2unix $i; done
```

## Running the analysis and correcting for housekeeping gene signal

> pwd: F:/SharedFolders/side_projects/hammond_goat_sheep_ig/

```R
library(RUVSeq)
library(dplyr)
library(tibble)
library(EDASeq)
library(RColorBrewer)
library(edgeR)
setwd("F:/SharedFolders/side_projects/hammond_goat_sheep_ig/"

# Let's test this out on one tissue first
gdata <- read.delim("goat_alv_data.tab", header=TRUE)
sdata <- read.delim("sheep_alv_data.tab", header=TRUE)

rownames(gdata) <- gdata$Run
rownames(sdata) <- sdata$Run
gdata <- gdata[,-1]
sdata <- sdata[,-1]

goatsamps <- colnames(gdata)
sheepsamps <- colnames(sdata)

cdata <- left_join(rownames_to_column(gdata), rownames_to_column(sdata), by="rowname")

ctrls <- c("SDHA", "PPIA", "GAPDH", "ACTB", "RPL13A", "YWHAZ")
genes <- cdata$rowname[!(cdata$rowname %in% ctrls)]

# Final sorting of data, now organizing the data into the appropriate object for analysis
rownames(cdata) <- cdata$rowname
cdata <- cdata[,-1]
set <- newSeqExpressionSet(as.matrix(cdata), phenoData = data.frame(as.factor(sapply(colnames(cdata), function(v){if(v %in% goatsamps){return("Ctl")}else{return("Trt")}})), row.names=colnames(cdata)))

# Normalize the data by the upper quartile
plotRLE(set, outline=FALSE, col= brewer.pal(3,"Set2"))
dev.copy2pdf(file="alv_prenorm.pdf")
set <- betweenLaneNormalization(set, which="upper")
dev.copy2pdf(file="alv_postuqnorm.pdf")

# Identify weighted factors for each sample and then regress read depth values against them
set1 <- RUVg(set, ctrls, k=1)
plotRLE(set1, outline=FALSE, col= brewer.pal(3,"Set2"))
dev.copy2pdf(file="alv_postuw_regression.pdf")

# Just checking to see how the samples span in a PCA plot
plotPCA(set1, col=brewer.pal(3,"Set2"), cex=1.2)
dev.copy2pdf(file="alv_postuw_regression_pca.pdf")

# Now, we can do an EdgeR analysis!
# I really should have changed the factor name for the newSeqExpression function call...
design <- model.matrix(~as.factor.sapply.colnames.cdata...function.v... + W_1, data=pData(set1))

y <- DGEList(counts=counts(set1), group=set1$as.factor.sapply.colnames.cdata...function.v...)
y <- calcNormFactors(y, method="upperquartile"
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
# This prints out the ordered list of genes that are likely differentially expressed
topTags(lrt)
# Writing the table to file
write.table(lrt$table, file="alv_norm_tags.tab", quote=FALSE, sep="\t")

heatmap(cpm(y, prior.count=2, log=TRUE)[seq(1,14),], scale="column", Rowv = NA)
dev.copy2pdf(file="alv_heatmap_postnorm.pdf")
```

## Streamlining the process

OK, so now that I know the workflow, let's turn this into a function.

```R
runRUVseqPipe <- function(outtag, ctrlgenes, goatfile, sheepfile){
gdata <- read.delim(goatfile, header=TRUE)
sdata <- read.delim(sheepfile, header=TRUE)

goatsamps <- colnames(gdata[,-1])

cdata <- left_join(gdata, sdata, by="Run")

rownames(cdata) <- cdata$Run
cdata <- cdata[,-1]

x <- as.factor(sapply(colnames(cdata), function(v){if(v %in% goatsamps){return("Ctl")}else{return("Trt")}}))
set <- newSeqExpressionSet(as.matrix(cdata), phenoData = data.frame(x, row.names=colnames(cdata)))

bpal <- brewer.pal(3,"Set2")
pdf(file=paste0(outtag, "_prenorm.pdf"))
plotRLE(set, outline=FALSE, col=bpal[x])
dev.off()

set <- betweenLaneNormalization(set, which="upper")
pdf(file=paste0(outtag, "_postuqnorm.pdf"))
plotRLE(set, outline=FALSE, col=bpal[x])
dev.off()

set1 <- RUVg(set, ctrlgenes, k=1)
pdf(file=paste0(outtag, "_postuw_regression.pdf"))
plotRLE(set1, outline=FALSE, col=bpal[x])
dev.off()

pdf(file=paste0(outtag, "_postuw_regression_pca.pdf"))
plotPCA(set1, col=bpal[x], cex=1.2)
dev.off()

design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

write.table(lrt$table, file=paste0(outtag, "_norm_tags.tab"), quote=FALSE, sep="\t")

pdf(file=paste0(outtag, "_heatmap_postnorm.pdf"))
heatmap(cpm(y, prior.count=2, log=TRUE)[seq(1,14),], scale="column", Rowv = NA)
dev.off()
}
```

Then it was trivial to queue up the rest of the runs.

```R
runRUVseqPipe("liver", ctrls, "goat_liver_data.tab", "sheep_liver_data.tab")
runRUVseqPipe("bone", ctrls, "goat_bone_data.tab", "sheep_bone_data.tab")
# The spleen exited with an error! I see why from the pre-norm data. One of the sheep is so off-base it screwed everything else up
runRUVseqPipe("spleen", ctrls, "goat_spleen_data.tab", "sheep_spleen_data.tab")
runRUVseqPipe("thymus", ctrls, "goat_thymus_data.tab", "sheep_thymus_data.tab")
```