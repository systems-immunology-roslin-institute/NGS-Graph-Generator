library(GenomicRanges)
library(ShortRead)
library(getopt)

optspec <- matrix(c(
	'gtfFile', 'g', 1, "character", "input GTF file (required)",
	'chrLengthFile', 't', 1, "character", "chromosome length table file (required)",
	'outFile', 'o', 1, "character", "output file name (optional)",
	'help', 'h', 0, "logical", "this help"
),ncol=5,byrow=T)

opt = getopt(optspec)

if(!is.null(opt$help) || is.null(opt$gtfFile) || is.null(opt$chrLengthFile))
{
	cat(paste(getopt(optspec, usage=T),"\n"))
	q()
}

if(!is.null(opt$outFile))
{
	outFile <- opt$outFile
} else
{
	outFile <- "ensembl_gtfannotation.RData"
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

# read gtf file
gtfFile <- opt$gtfFile

if(!file.exists(gtfFile))
{
	cat(sprintf("GTF file '%s' does not exist or cannot be accessed.\n", gtfFile))
	q("no", 1)
}

cat(sprintf("Loading GTF file '%s'...\n", gtfFile))
gtf <- read.table(gtfFile, header=F,sep="\t")

gtf <- gtf[gtf$V3=="exon",]
gtf$gene_id = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",3))
gtf$transcript_id = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",5))
gtf$exon_number = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",7))
gtf$gene_name = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",9))
gtf$transcript_name = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",13))
gtf$biotype = sub(";","",sapply(strsplit(as.character(gtf$V9)," "), "[",11))
gtf <- gtf[gtf$V1 %in% names(chrLens),]

cat("Running GRanges...\n")
GR  <-   GRanges(
         seqnames   = Rle(factor(gtf$V1,levels=names(chrLens))),
         ranges     = IRanges(start=gtf$V4, end=gtf$V5),
         strand     = Rle(factor(gtf$V7)),
         seqlengths = chrLens,
         geneid       = gtf$gene_id,
	 transcriptid = gtf$transcript_id,
	 exonnumber = gtf$exon_number,
	 genename   = gtf$gene_name,
	 biotype = gtf$biotype
      )

cat(sprintf("Writing output file '%s'.\n", outFile))
save(GR,file=outFile)
