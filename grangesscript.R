library(GenomicRanges)
library(ShortRead)
library(getopt)

optspec <- matrix(c(
	'bamFile', 'b', 1, "character", "input BAM file (required)",
	'chrLengthFile', 't', 1, "character", "chromosome length table file (required)",
	'outFile', 'o', 1, "character", "output file name (optional)",
	'help', 'h', 0, "logical", "this help"
),ncol=5,byrow=T)

opt = getopt(optspec)

if(!is.null(opt$help) || is.null(opt$bamFile) || is.null(opt$chrLengthFile))
{
	cat(paste(getopt(optspec, usage=T),"\n"))
	q()
}

if(!is.null(opt$outFile))
{
	outFile <- opt$outFile
} else
{
	outFile <- "GRanges.RData"
}

chrLenFile <- opt$chrLengthFile

if(!file.exists(chrLenFile))
{
	cat(sprintf("Chromosome length table file '%s' does not exist or cannot be accessed.\n", chrLenFile))
	q("no", 1)
}

cat(sprintf("Loading chromosome length file '%s'...\n", chrLenFile))
chromLengthTab <- read.table(chrLenFile, sep="\t", as.is=TRUE, header = TRUE,colClasses = c("character", "integer"))
chromLengthTab <- chromLengthTab[order(chromLengthTab[,1]),]
chrLens <- chromLengthTab[,2]
names(chrLens) <- chromLengthTab[,1]

bamFile <- opt$bamFile

if(!file.exists(bamFile))
{
	cat(sprintf("BAM file '%s' does not exist or cannot be accessed.\n", bamFile))
	q("no", 1)
}

cat(sprintf("Scanning BAM file '%s'...\n", bamFile))
what <- c("qname","rname","strand","pos","qwidth","seq")
flag <- scanBamFlag(isUnmappedQuery=FALSE)
param <- ScanBamParam(what=what, flag=flag)
bamHITS <- scanBam(bamFile,param=param)[[1]]

cat("Collating data...\n")
sequences <- as.data.frame(bamHITS$seq)	
input <- data.frame(bamHITS$qname,bamHITS$rname,bamHITS$strand,bamHITS$pos,bamHITS$qwidth,stringsAsFactors=FALSE)
input$Sequence <- as.character(sequences[,1]) 
colnames(input) <- c("Read","Seq","Strand","Start","Width","Sequence")

input <- input[input$Seq %in% names(chrLens),]

cat("Running GRanges...\n")
GRbam  <-   GRanges(
         seqnames   = Rle(factor(input$Seq,levels=names(chrLens))),
         ranges     = IRanges(start=input$Start, width=input$Width),
         strand     = Rle(factor(input$Strand)),
         seqlengths = chrLens,
         name       = input$Read,
	 sequence = input$Sequence
      )

cat(sprintf("Writing output file '%s'.\n", outFile))
save(GRbam,file=outFile)
