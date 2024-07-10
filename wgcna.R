### Running WGCNA on transcriptomic data from SHANK3 haploinsufficiency patients

## Load packages
library(WGCNA)
library(doParallel)
library(preprocessCore)
library(parallel)
library(reshape2)
library(dplyr)
library(AnnotationDbi)
library(dynamicTreeCut)
library(fastcluster)
library(flashClust)
library(multtest)
library(igraph)
library(org.Hs.eg.db)
library(tximport)
library(DESeq2)
library(ggrepel)

## Loading files

metadata = read.delim(file = "~/RNAseq-Elisa/wgcna_samples.txt", check.names = F)

rownames(metadata) = metadata$individual
metadata = metadata[!metadata$individual %in% "C5",] # remove sample C5 since it was discarded

setwd("~/RNAseq-Elisa/diff_exp") # path containing all quantification files
samples <- data.frame(row.names=metadata$individual,
                      condition=c(rep("control",6),rep("patient",11)),
                      sample=metadata$individual)

files <- list.files(pattern = "genes.results")
names(files) <- gsub("(.+?)(\\_.*)", "\\1", files)
all(file.exists(files))

txi = tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] = 0.01

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)
smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)

vsd<- varianceStabilizingTransformation(dds) # necessary normalization transformation for WGCNA
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
ncvsd2 <- ncvsd

dup_genes <- gsub(".*_","",rownames(ncvsd2))
ncvsd2 <- ncvsd2[!duplicated(dup_genes),] # removing the few genes with duplicated gene symbol
ncvsd2 <- as.data.frame(ncvsd2)
ncvsd2$Symbol <- gsub(".*_","",rownames(ncvsd2))
rownames(ncvsd2) <- ncvsd2$Symbol
ncvsd2 <- ncvsd2[,-18]

# Check for variance over the mean and possibliy remove some genes

variance <- apply(ncvsd2,1,sd)/apply(ncvsd2,1,mean)
hist(variance)
summary(variance)

keep=which(apply(ncvsd2,1,sd)/apply(ncvsd2,1,mean)>= 0.025) # calculated the variance over the mean and defined a minimum value in order to have variance enough to detect co-expression
ncvsd2=ncvsd2[keep,]

# Removing genes with low expression

dat=t(ncvsd2)
infoData=rownames(ncvsd2)
dim(dat)

gsg<-goodSamplesGenes(dat,verbose=3) # check for no/low counts transcripts
gsg$allOK
datExpr0 = dat[gsg$goodSamples, gsg$goodGenes]
dim(datExpr0)

## Hierarchical clustering

sampleTree = flashClust(dist(datExpr0), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

pdf(file="~/RNAseq-Elisa/hierarchical_cluster_inicial.pdf")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 125, col = "red");
dev.off()

## Run WGCNA

# run soft pick threshold
# test a series of powers to which co-expression similarity is raised to calculate adjacency

allowWGCNAThreads()
powers=c(c(1:10), seq(from = 10, to = 30, by = 2)) 
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, blockSize = 30000, networkType = "signed") # Signed gives an indication of positive vs negative correlations
write.table(sft,'~/RNAseq-Elisa/powers_table.txt', sep="\t",quote=F)

png(
  "~/RNAseq-Elisa/WGCNA_SoftThreshold.png",
  width     = 9,
  height    = 4.5,
  units     = "in",
  res       = 1200
)
par(
  mfrow = c(1,2),
  mar = c(2,5,2,5),
  cex = 0.6
)
cex1 = 0.9
# Fit to scale-free topology 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence")); 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red"); 
abline(h=c(0.8, 0.5), col="red")
# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold (power)", 
     ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col="red")
dev.off()

## Construction of the network

# Chosen soft threshold in this case is 24; it sets the fit to Scale-free Topology > 0.8 while retaining as much connectivity as possible
net=blockwiseModules(datExpr0, power=24, numericLabels=TRUE, networkType = "signed", deepSplit = 1,
                     minModuleSize=50, mergeCutHeight=0.2, saveTOMs=FALSE, verbose=6, minKMEtoStay = 0.1,
                     nThreads=24, maxBlockSize=30000, checkMissingData=FALSE)


## Labelling the modules with a colour tag
modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N") 
modules$Label=paste("M", modules$Label, sep=""); 
modules$Color=c("grey",labels2colors(modules$Label[-1])) 
moduleLabel=paste("M",net$colors, sep=""); 
moduleColor=modules$Color[match(moduleLabel, modules$Label)] 

### Get kME table

## Calculating kMEs
KMEs<-signedKME(datExpr0, net$MEs,outputColumnName = "M") # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme=data.frame(infoData[match(colnames(datExpr0), infoData)], moduleColor,moduleLabel, KMEs)
colnames(kme)[1]="Symbol"

#Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedata=kme[,-kmeInfoCols]; 
pvalBH=kmedata; pvalBH[,]=NA 
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedata[,j], nSamples=17), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme$newModule="NA" 
for (j in c(1:nrow(kmedata)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m=which(kmedata[j,]==max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedata[j,m]>0.5)) kme$newModule[j]=as.character(colnames(kmedata)[m])
}

## Assign genes not associated to any module to M0
kme$newModule[which(kme$newModule%in%"NA")]="M0" 

## Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)] 
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor 
kme=kme[,-grep("newModule", colnames(kme))];kme=kme[,-grep("newColor", colnames(kme))] 

# Saving kMEs 
mod=modules$Label[-1] 
kmeTable=kme[,kmeInfoCols]; 
for(j in c(1:length(mod)))
{
  kmeTable=cbind(kmeTable, kmedata[,match(mod[j],colnames(kmedata))]);colnames(kmeTable)[ncol(kmeTable)]=paste("kME", mod[j], sep="_")                                                                             
  kmeTable=cbind(kmeTable, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTable)[ncol(kmeTable)]=paste("pvalBH", mod[j], sep="_")
}

# kME Sup Table
write.csv(kmeTable, "~/RNAseq-Elisa/kME_total.csv", row.names=FALSE)

# Saving Module Eigengenes 
me<-data.frame(rownames(datExpr0), net$MEs) # Sample names bound to module eigengenes
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me) <- modules$Color[match(colnames(me),modules$Label)]
colnames(me)[1]="Sample"
write.csv(me, "~/RNAseq-Elisa/ME_total.csv", row.names=FALSE)

#Plotting Module Eigengene Values
colors=rep("red", nrow(me) )
colors[grep("C", me$Sample)]="turquoise"
  
# Barplots with kME of samples for each module
for(m in colnames(me)[-1]) 
{ 
  png(
    paste("~/RNAseq-Elisa/moduleBarplots_total",m,".png"),
    width     = 10,
    height    = 2,
    units     = "in",
    res       = 1200
  )
  par(
    mfrow = c(1,2),
    mar = c(2,5,2,5),
    cex = 0.6
  )
  j=match(m, colnames(me))
  col=kme$moduleColor[match(m, kme$moduleLabel)]
  barplot(me[,j],  xlab="Samples", ylab="ME",col=colors, main=m, names=me[,1], cex.names=0.5, axisnames = FALSE)
  dev.off()
}

## Dendogram of modules after genes were reassigned due to kME
  
png(
  "~/RNAseq-Elisa/cluster_dendogram.png",
  width     = 10,
  height    = 7.5,
  units     = "in",
  res       = 1200
)
par(
  mfrow = c(1,2),
  mar = c(2,5,2,5),
  cex = 0.6
)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## Get biological variables correlations

# Define numbers of genes and samples

Samples = rownames(datExpr0);
traitRows = match(Samples, tmp_sampleInfo$individual);
SampleInfo = tmp_sampleInfo[traitRows, -1];
rownames(SampleInfo) = tmp_sampleInfo[traitRows, 1];
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
SampleInfo <- apply(SampleInfo,2,as.numeric)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColor)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, SampleInfo, use = "p");

x <- MEs
y <- cbind(x,SampleInfo)

for (i in 1:25){
  z = y[,c(i,26:31)]
  colnames(z)[1] <- "modulo"
  tmp_lm <- lm(formula = modulo ~ sex + age + mut_type, data = z)
  summ_tmp_lm = summary(tmp_lm)
  print(colnames(y)[i])
  print(summ_tmp_lm$coefficients)
}

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

png(
  "~/RNAseq-Elisa/Module-TRait_total.png",
  width     = 12,
  height    = 5,
  units     = "in",
  res       = 1200
)
par(
  mfrow = c(1,2),
  mar = c(2,5,2,5),
  cex = 0.6
)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(SampleInfo),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

## Get clinical variables correlations
library(dplyr)

Samples = rownames(datExpr0);

clinicalTraits = read.delim(file = "~/RNAseq-Elisa/dados_traços_clinicos.txt", check.names = F)
clinicalTraits$P5[clinicalTraits$P5 == "1,5"] <- 1.5
clinicalTraits$P5 <- as.numeric(clinicalTraits$P5)
colnames(clinicalTraits)
clinicalTraits = t(clinicalTraits)
colnames(clinicalTraits) = clinicalTraits[1,]
clinicalTraits = as.data.frame(clinicalTraits[-1,])
clinicalTraits <- mutate_all(clinicalTraits, function(x) as.numeric(as.character(x)))
clinicalTraits$individual = rownames(clinicalTraits)
clinicalTraits = clinicalTraits[,c(-1,-6,-12,-13,-14)]
Samples[Samples == "P8"] <- "P2.2"
traitRows = match(Samples, clinicalTraits$individual);
SampleInfo = clinicalTraits[traitRows,];
SampleInfo <- SampleInfo[grepl("P",rownames(SampleInfo)),]
Samples <- Samples[grepl("P",Samples)]
traitRows <- na.omit(traitRows)
rownames(SampleInfo) = clinicalTraits[traitRows, 22];
SampleInfo = SampleInfo[,-22];
datExprCorr <- datExpr0[grepl("P",rownames(datExpr0)),]
nGenes = ncol(datExprCorr);
nSamples = nrow(datExprCorr);
SampleInfo <- apply(SampleInfo,2,as.numeric)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprCorr, moduleColor)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, SampleInfo, use = "p");

x <- MEs
y <- cbind(x,SampleInfo)

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

moduleTraitCor <- moduleTraitCor[rownames(moduleTraitCor) %in% c("MEturquoise", "MElightgreen", "MEroyalblue","MEgrey"),]
moduleTraitPvalue <- moduleTraitPvalue[rownames(moduleTraitPvalue) %in% c("MEturquoise", "MElightgreen", "MEroyalblue","MEgrey"),]

png(
  "~/RNAseq-Elisa/Module-clinicalTrait_selected.png",
  width     = 10,
  height    = 3,
  units     = "in",
  res       = 1200
)
par(
  #mfrow = c(1,2),
  mar = c(2,5,2,5),
  cex = 0.6
)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

nMEs <- MEs[,colnames(MEs) %in% c("MEturquoise", "MElightgreen", "MEroyalblue","MEgrey")]

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(SampleInfo),
               yLabels = names(nMEs),
               ySymbols = names(nMEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

### Module preservation using BrainSpan data

brainspan <- read.table("~/RNAseq-Elisa/brainspan/expression_matrix.csv", sep = ",")
brainspan <- brainspan[,-1]
header <- read.table("~/RNAseq-Elisa/brainspan/columns_metadata.csv", sep = ",", header = TRUE)

header_filtered <- header[header$age %in% c("16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"),]
header_filtered <- header_filtered[header_filtered$structure_name %in% c("occipital neocortex","primary motor-sensory cortex (samples)","posterior (caudal) superior temporal cortex (area 22c)","anterior (rostral) cingulate (medial prefrontal) cortex","dorsolateral prefrontal cortex","orbital frontal cortex","inferolateral temporal cortex (area TEv, area 20)","ventrolateral prefrontal cortex","parietal neocortex","temporal neocortex","primary auditory cortex (core)","primary visual cortex (striate cortex, area V1/17)","primary motor cortex (area M1, area 4)","posteroventral (inferior) parietal cortex","primary somatosensory cortex (area S1, areas 3,1,2)","medial ganglionic eminence","caudal ganglionic eminence","lateral ganglionic eminence"),]
header_filtered$filter <- paste(header_filtered$donor_id,header_filtered$age,header_filtered$structure_name,sep = "_")

genes_brainspan <- read.table("~/RNAseq-Elisa/brainspan/rows_metadata.csv", sep = ",", header = TRUE)
dup <- duplicated(genes_brainspan$gene_symbol)
brainspan <- brainspan[!dup,]
genes_brainspan <- genes_brainspan[!dup,]
row.names(brainspan) <- genes_brainspan$gene_symbol
colnames(brainspan) <- paste(header$donor_id,header$age,header$structure_name,sep = "_")
brainspan_filtered <- brainspan[,na.omit(match(colnames(brainspan),header_filtered$filter))]

keep <- rowSums(brainspan_filtered >= 1) >= 2
brainspan_filtered2 <- brainspan_filtered[keep,]

datBS=t(brainspan_filtered2)
infoData=rownames(brainspan_filtered2)
dim(datBS)

gsg<-goodSamplesGenes(datBS,verbose=3) # check for no/low counts transcripts
gsg$allOK
datBS00 = datBS[gsg$goodSamples, gsg$goodGenes]
dim(datBS00)

datBS0 <- -log2(datBS00 + 1)

# Create an object with module assignment of original dataset (datExpr0)

colorsNetwork = kmeTable$moduleColor[match(colnames(datExpr0),rownames(kmeTable))]
table(colorsNetwork == kmeTable$moduleColor)

## Module preservation

setLabels = c("PMDS"," BrainSpan")
multiExpr = list(PMDS = list(data = datExpr0), BrainSpan = list(data = datBS0))
multiColor = list(PMDS = colorsNetwork)

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );

# Save the results
save(mp, file = "~/RNAseq-Elisa/modulePreservation.RData")

# Analysis and display of module preservation results

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

write.table(statsObs,'~/RNAseq-Elisa/statsObs_preservationAnalysis.txt', sep="\t",quote=F)
write.table(statsZ,'~/RNAseq-Elisa/statsZ_preservationAnalysis.txt', sep="\t",quote=F)

## Plot results

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
#colnames(plotData) <- c("observed","Zscore")
#plotData = plotData[plotMods,]
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
#sizeGrWindow(10, 5);

png(
  "~/RNAseq-Elisa/modulePreservation-Zsummary-medianRankFinal.png",
  width     = 6,
  height    = 11,
  units     = "in",
  res       = 1200,
  pointsize = 31
)
par(
  mar = c(4.5,4.5,2.5,1),
  cex = 0.6
)
p = 2
min = min(plotData[, 1], na.rm = TRUE);
max = max(plotData[, 1], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
  if (min > -max/10) min = -max/10
  #min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
  ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
pos1 = rep(3, length(text))
pos1[text %in% c("darkturquoise","purple","tan","brown","midgnightblue")] <- 1

plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
     main = mains[p],
     cex = 2.4,
     ylab = mains[p], xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
geom_text_repel(moduleSizes[plotMods], plotData[plotMods, p], labels = text, cex = 1.5)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 0.5, pos=pos1);, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
  abline(h=0)
  abline(h=5, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();
