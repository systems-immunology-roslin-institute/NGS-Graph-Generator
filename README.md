Next Generation Sequencing Graph Generator
==========================================

This is a script that is used to take next generation gene sequencing data and perform read to read comparison in order to construct a graph based representation. The output is in the form of a .layout file which may be opened by the graph visualisation tool [BioLayout](http://www.biolayout.org/).

Prerequisites
=============

Python >= 2.7.3
samtools >= 0.1.18
R >= 2.15.2

Authors
=======

The scripts provide are based on the work of others. In particular:

* Harpreet Kaur
.* findoverlaps.R
.* grangesscript.R
.* grangesscript_gtf.R

* Fioravante de Sapio
.* read2read.py
.* megablast2ncol.py

* Tim Angus
.* create-biolayout-file.sh
.* tab-to-fasta.sh
.* tab-to-nodeclass.sh
.* Modifications to above scripts
