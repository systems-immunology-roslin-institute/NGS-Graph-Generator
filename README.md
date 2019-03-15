NGS Graph Generator: graph-based visualisation of RNA-seq data
==========================================

NGS graph generator supports the graph-based visualisation of RNA-seq data and provides a complementary approach to understanding transcript diversity and issues with assembly.  Following the mapping of reads to the reference genome, read to read comparison is performed on portions of the data (e.g. all reads mapping to a given gene) using Mega BLAST to provide a matrix of weighted similarity scores between reads.  

This is a script that is used to take next generation gene sequencing data and perform read to read comparison in order to construct a graph based representation. The output is in the form of a .layout file which may be opened by the graph visualisation tool [Graphia](https://kajeka.com/).

Prerequisites
-------------
* Python version ≥ 2.7.3 (https://www.python.org/)
* SAMtools version ≥ 0.1.18 (http://samtools.sourceforge.net/)
* R package version ≥ 2.15.2 (http://www.r-project.org/)
* Mega Blast version ≥  2.2.19 (http://blast.ncbi.nlm.nih.gov/)

System requirements
-------------------
* Linux operating system
* 16GB RAM recommended

Each of prerequisites needs to be installed on the server before running NGS graph generator. 

Running the program
-------------------

Input format
* a) Mapping file from alignment output (BAM file) such as TopHat or Bowtie
* b)	GTF/GFF file (the same GTF/GFF file used when running alignment)
* c)	Chromosome length file

To run the program, type the command-line:

```
$ create-biolayout-file.sh –b <bam file> -t <The chromosome length file> -g <The GTF file> -o <The directory in which to place the output> -d <"GENE" A list of genes to examine> -p <The percentage similarity value (default 98)> -l <The percentage coverage value (default 31)> -u <discard redundant reads> 
```

Where:
  * –b is a bam file after mapping process
  * –t is a chromosome length file 
  * –g is a GTF/GFF file for annotation of the nodes
  * –o is the desired output directory of this process
  * –d is an interested gene or set of genes for the analysis
  * –p is a percentage sequence similarity (overlap) between two reads (default = 98)
  * –l is a percentage coverage value (default = 31)
  * –u is an optional to remove redundant reads

For graph visualisation
Output format

The output of this program is a layout (.layout) file which can visualise using BioLayout Express3D software. Layout file is a text file that contains pairwise comparison of all reads of selected gene with edge weight. It also contains classsets which defines the nodes. For instance there will be a classsets called nodeclass which defines the nodes to exon and transcripts which it belongs.

In this example below is a .layout format which is a text file contains all information of selected gene. 

|Read1 				      |Read2 			      |Edge weight |
|:--------------------|:--------------------|:------------|   
| 1626463_t10_w100_x1 |	22175909_t12_w100_x1 | 185        |
| 1626463_t10_w100_x1 |	21060673_t12_w100_x1 | 185        |
| 1626463_t10_w100_x1 |	19760852_t13_w100_x1 | 185        |
| 1626463_t10_w100_x1 | 17188395_t12_w100_x1 | 185        |
….

| //NODECLASS	| //GENEID            | //EXONNUMBER | //TRANSCRIPTID  |
|:-------------|:---------------------|:--------------|:-----------------|
| //NODECLASS |	5431875_t14_w100_x1 | Exon 1       | ENST00000350051 |
| //NODECLASS	| 9230189_t14_w100_x1 | Exon 1       | ENST00000301633 |
| //NODECLASS	| 9230189_t14_w100_x1 | Exon 1       | ENST00000350051 |
| //NODECLASS	| 9230189_t14_w100_x1 | Exon 1       | ENST00000590946 |
| //NODECLASS	| 9230189_t14_w100_x1 | Exon 1       | ENST00000374948 |
| //NODECLASS	| 9230189_t14_w100_x1 | Exon 1       | ENST00000590449 |
….


For graph visualisation of RNA-seq data, a user has to download:
•	Graphia software (https://kajeka.com/)
Please note Graphia uses high end 3D visualisation technology and the machine running the program will require a good graphics card, at least 4 Gb of RAM and preferably a 64 bit OS.  The better your machine the bigger the graphs you will be able to render and the smoother your interaction with them.






 
Authors
-------

The scripts provided are based on the work of others. In particular:

* Harpreet Kaur and Fahmi W. Nazarie
 * findoverlaps.R
 * grangesscript.R
 * grangesscript_gtf.R

* Fioravante de Sapio and Fahmi W. Nazarie
 * read2read.py
 * megablast2ncol.py

* Tim Angus and Fahmi W. Nazarie
 * create-biolayout-file.sh
 * tab-to-fasta.sh
 * tab-to-nodeclass.sh
 * Modifications to above scripts


