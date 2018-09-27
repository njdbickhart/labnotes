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

## Generating a full dataset heatmap with relative abundance estimates

> pwd: F:/SharedFolders/side_projects/hammond_goat_sheep_ig/

```R
library(reshape2)
load(file=".RData")


# Goat first
gtables <- lapply(list.files(pattern="goat.*\\.tab"), read.delim, header=TRUE)
gtables <- lapply(gtables, function(x){rownames(x) <- x$Run; x["Run"] <- NULL; x})
gdataframe <- do.call(cbind, gtables)

gfull <- newSeqExpressionSet(as.matrix(gdataframe))
gfull1 <- RUVg(gfull, ctrls, k=1)
# There were some genes that had very, very little expression post normalization
drops <- c("TARM1", "OSCAR", "KIR group 2", "KIR group 6", "GP6")
gfull1.m <- melt(normCounts(gfull1)[!(rownames(normCounts(gfull1)) %in% c(ctrls, drops)), ])

# Now I need to change the goat ids into their tissue type so that I can format the image easier
gtables <- lapply(list.files(pattern="goat.*\\.tab"), read.delim, header=TRUE)
names(gtables) <- sapply(list.files(pattern="goat.*\\.tab"), function(x) strsplit(x, split="\\_", perl=TRUE)[[1]][2])
gtissues <- Map(function(x, i){x["Run"] <- NULL; data.frame(tissue=i, sample=colnames(x))}, gtables, names(gtables))
gtissues <- do.call(rbind, gtissues)

pdf(file="goat_draft_heatmap_norm.pdf", useDingbats = FALSE)
ggplot(gfull1.m, aes(Var2, Var1)) + geom_tile(aes(fill=value), colour = "white") + scale_fill_gradient(low="white", high="steelblue", limits=c(0,500)) + scale_x_discrete(labels=gtissues$tissue) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))

# Now Sheep
stables <- lapply(list.files(pattern="sheep.*\\.tab"), read.delim, header=TRUE)
names(stables) <- sapply(list.files(pattern="sheep.*\\.tab"), function(x) strsplit(x, split="\\_", perl=TRUE)[[1]][2])
stissues <- Map(function(x, i){x["Run"] <- NULL; data.frame(tissue=i, sample=colnames(x))}, stables, names(stables))
stissues <- do.call(rbind, stissues)
stables <- lapply(stables, function(x){rownames(x) <- x$Run; x["Run"] <- NULL; x})
sdataframe <- do.call(cbind, stables)

sfull <- newSeqExpressionSet(as.matrix(sdataframe))
sfull1 <- RUVg(sfull, ctrls, k=1)
drops <- c("TARM1", "OSCAR", "NCR1", "GP6")
sfull1.m <- melt(normCounts(sfull1)[!(rownames(normCounts(sfull1)) %in% c(ctrls, drops)), ])

pdf(file="sheep_draft_heatmap_norm.pdf", useDingbats=FALSE)
ggplot(sfull1.m, aes(Var2, Var1)) + geom_tile(aes(fill=value), colour = "white") + scale_fill_gradient(low="white", high="steelblue", limits=c(0,500)) + scale_x_discrete(labels=stissues$tissue) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
```

So, I've created a common scale for the heatmap that peaks at about 500+ normalized, expressed transcripts. I will need to edit the image pdfs in inkscape to adjust the colors of those tiles (only because ggplot apparently colors them dark grey by default!). I will then need to replace the x axis names with bars for tissues. Finally, the plots should be ready for John's review.

## Preparing the last data files

```R
# Writing the normalized counts
write.table(normCounts(sfull1), file="sheep_normcounts_full.tab", sep = "\t", quote=FALSE)
write.table(normCounts(gfull1), file="goat_normcounts_full.tab", sep="\t", quote=FALSE)

# Testing something for John
library(dplyr)

sfull1.fm <- melt(normCounts(sfull1))
sdataframe <- mutate(sdataframe, ids = rownames(sdataframe))
sraw.fm <- melt(sdataframe, id.vars = "ids")
sfull1.fm$Raw <- sraw.fm$value
sfull1.fm <- mutate(sfull1.fm, tissue=str_extract(Var2, "^[^\\.$]"))
ggplot(sfull1.fm, aes(x=Raw, y=value, col=tissue, label=Var1)) + geom_point() + ylab(label = "Normalized value") + xlab(label="Raw Read Count") + geom_text(aes(label=ifelse(Var1 %in% ctrls, as.character(Var1),'')), hjust=0, vjust=-1)

dev.copy2pdf(file="sheep_example_normalization_effects.pdf", useDingbats=FALSE)
```