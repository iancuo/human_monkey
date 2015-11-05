library(edgeR)

setwd("/Users/iancuo/workDir/macaque/CEA")

geneReads =as.matrix(read.csv("countData/geneCountsCEA.csv",header=T,row.names=1))
geneNames=rownames(geneReads)

totalSampleReads=colSums(geneReads)

geneReads=geneReads[,1:30]
###################################################################################
UQnormFactors=calcNormFactors(geneReads, method=c("upperquartile"))

effectiveLibrarySizes= UQnormFactors*colSums(geneReads)
meanEffLibSize=mean(effectiveLibrarySizes)
countNormFactor= meanEffLibSize/effectiveLibrarySizes

normalizedGeneCountsUQ=0* geneReads

for (sample in 1:dim(normalizedGeneCountsUQ)[2]){
	
	normalizedGeneCountsUQ[,sample]= geneReads[, sample]* countNormFactor[sample]	
}

# librarySizesNormalized=colSums(normalizedGeneCountsUQ)

# plot(colSums(geneReads), countNormFactor)
# sd(librarySizesNormalized)
# max(librarySizesNormalized)-min(librarySizesNormalized)

selectedIdxs=vector(mode="numeric",length=dim(normalizedGeneCountsUQ)[1])
for (gene in 1:dim(normalizedGeneCountsUQ)[1]){
	meanCounts=mean(unlist(normalizedGeneCountsUQ[gene,]))
	
	if (meanCounts >=100){		
		selectedIdxs[gene]=1
	}

}
sum(selectedIdxs)
selectedIdxs =as.logical(selectedIdxs)
normalizedGeneCountsUQSelected = round(normalizedGeneCountsUQ[selectedIdxs,])
colnames(normalizedGeneCountsUQSelected)=colnames(geneReads)
save(normalizedGeneCountsUQSelected, countNormFactor, file="data/normGenesAndFactors.RData")
###################################################################################

exonCounts =as.matrix(read.csv("countData/exonCountsCEA.csv",header=T))
exonInfo= read.csv("countData/exonInfo.csv",header=T)

exonCounts=exonCounts[,1:30]

exonGeneTranscriptName= as.character(unlist(exonInfo[,"Gene"]))
exon_start=exonInfo[,"start"]
splitIDs=mapply(strsplit, exonGeneTranscriptName, MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 1))
exon_unique_id=mapply("paste", exonGeneName, exon_start, MoreArgs=list(sep="_"))

rownames(exonCounts)= exon_unique_id
sampleNames=as.vector(colnames(exonCounts))

colnames(exonCounts)==colnames(geneReads) # sanity check


selectedExonCounts=0* exonCounts

for (sample in 1:dim(exonCounts)[2]){
	
	selectedExonCounts[,sample]= exonCounts[, sample]* countNormFactor[sample]	
}

selectedExonCounts =round(selectedExonCounts)
rownames(selectedExonCounts)=rownames(exonCounts)

############################################################################################
# filter out by setting to 0 exons that have fewer than 10 average number of reads
# exon counts with 0 do not contribute to Canberra metric

for(exon in 1:dim(selectedExonCounts)[1]){
	meanCounts=mean(unlist(selectedExonCounts[exon,]))
	
	if (meanCounts < 10){		
		selectedExonCounts[exon,]=0
	}
	
}

#selectedGenes=intersect(exonGeneName[selectedIdxsExons==1], rownames(normalizedGeneCountsUQSelected))

zeroExonCountGenes=c()

for (i in 1:length(geneNames)){
	geneName=geneNames[i]

	#geneCounts= selectedGeneCounts[geneName,]

	exonCounts= selectedExonCounts[which(exonGeneName ==geneName),]
	
	if (sum(exonCounts)==0){
		zeroExonCountGenes=c(zeroExonCountGenes, geneName)
		
	}
}

exonGeneNameSelected=setdiff(exonGeneName, zeroExonCountGenes)

selectedGeneCounts= normalizedGeneCountsUQSelected[intersect(rownames(normalizedGeneCountsUQSelected), exonGeneNameSelected),]


###################################################################################
###################################################################################

save(selectedExonCounts, exonGeneName, exon_start,  selectedGeneCounts,  file="data/selectedCountData.RData")
load("data/selectedCountData.RData")



