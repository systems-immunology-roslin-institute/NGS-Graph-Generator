#! /bin/bash

TAB_FILE=$1
BASE_NAME=$(basename "$TAB_FILE")
NO_EXT="${BASE_NAME%.*}"

awk '{print "//NODECLASS\t\"" $3 "\"\t\"Exon " $7 "\"\t\"" $6 "\""}' ${TAB_FILE} | \
  sed 's,"readname",\//GENEID,' | sed 's,"Exon exonnumber",\//EXONNUMBER,' | \
  sed 's,"transcriptid",\//TRANSCRIPTID,'| sed '2d'
