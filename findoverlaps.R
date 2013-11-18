library(GenomicRanges)
library(Biostrings)
library(getopt)

optspec <- matrix(c(
	'grangesFile', 'g', 1, "character", "GRanges.RData file (required)",
	'gtfAnnotationFile', 'e', 1, "character", "ensembl_gtfannotation.RData file (required)",
	'fileGeneList', 'a', 1, "character", "a file containing a list of genes to target (required)",
	'geneList', 'd', 1, "character", "a comma separated list of genes to target (required)",
	'outFilePrefix', 'p', 1, "character", "text to prefix the output files with (optional)",
	'help', 'h', 0, "logical", "this help"
),ncol=5,byrow=T)

opt = getopt(optspec)

if(!is.null(opt$help) || is.null(opt$grangesFile) || is.null(opt$gtfAnnotationFile) ||
  (is.null(opt$geneList) && is.null(opt$fileGeneList)))
{
	cat(paste(getopt(optspec, usage=T),"\n"))
	q()
}

if(!is.null(opt$geneList))
{
  targetGenes <- scan(textConnection(gsub("\\s*,\\s*|\\s+", " ", opt$geneList, perl=TRUE)), what="character", sep=" ")
  targetGenes = targetGenes[targetGenes != ""]
} else
{
  targetGenes <- readLines(opt$fileGeneList)
}

cat(sprintf("Processing %d genes...\n", length(targetGenes)))

# load GRanges of Ensembl human gtf annotations (based on hg19)
if(!file.exists(opt$gtfAnnotationFile))
{
	cat(sprintf("GTF annotation file '%s' does not exist or cannot be accessed.\n", opt$gtfAnnotationFile))
	q("no", 1)
}

cat(sprintf("Loading GTF annotation file '%s'...\n", opt$gtfAnnotationFile))
load(opt$gtfAnnotationFile)

# load tophat2 aligned Rdata
if(!file.exists(opt$grangesFile))
{
	cat(sprintf("GRanges file '%s' does not exist or cannot be accessed.\n", opt$grangesFile))
	q("no", 1)
}

cat(sprintf("Loading GRanges file '%s'...\n", opt$grangesFile))
load(opt$grangesFile)

# findoverlaps between aligned reads and genomic annotations
cat("Finding overlaps...\n")
overlaps <- findOverlaps(GRbam,GR)

cat("Collating data")
overs <- data.frame(NA,rownames=c(1:length(overlaps)))
cat(".")
overs$queryHits<-queryHits(overlaps)
cat(".")
overs$subjectHits<-subjectHits(overlaps)

cat(".")
overs$readname <- values(GRbam)["name"][overs$queryHits,]
cat(".")
overs$readSeq <-(values(GRbam)[["sequence"]][overs$queryHits])
cat(".")
overs$geneid <-(values(GR)[["geneid"]][overs$subjectHits])
cat(".")
overs$transcriptid <-(values(GR)[["transcriptid"]][overs$subjectHits])
cat(".")
overs$exonnumber <-(values(GR)[["exonnumber"]][overs$subjectHits])
cat(".")
overs$genename <-(values(GR)[["genename"]][overs$subjectHits])
cat(".")
overs$biotype <-(values(GR)[["biotype"]][overs$subjectHits])

cat(".")
overs <- overs[,-c(1,2)]
cat("\n")

if(is.null(opt$outFilePrefix))
{
	prefix <- ""
} else
{
	prefix <- opt$outFilePrefix
}

for(t in 1:length(targetGenes))
{
	outFile = paste(prefix, targetGenes[t], ".tab",sep="")
	cat(sprintf("Writing '%s'...\n", outFile))
	oversSubset <- overs[overs$genename==targetGenes[t],]
	write.table(oversSubset,outFile,row.names=F,col.names=T,sep="\t",quote=F)
}

