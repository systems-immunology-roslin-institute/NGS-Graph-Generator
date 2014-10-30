#! /bin/bash

while getopts uet: ARG
do
  case ${ARG} in
    (u) UNIQUIFY="1";;
    (e) EXONS="1";;
    (t) TAB_FILE="${OPTARG}";;
  esac
done

if [ -z "${TAB_FILE}" ];
then
    echo "tab file not supplied"
    exit 1
fi

BASE_NAME=$(basename "$TAB_FILE")
NO_EXT="${BASE_NAME%.*}"

if [ "${EXONS}" == "1" ];
then
    awk '{print "//NODECLASS\t\"" $3 "\"\t\"Exon " $7 "\"\t\"" $6 "\""}' ${TAB_FILE} | \
      sed 's,"readname",\//GENEID,' | sed 's,"Exon exonnumber",\//EXONNUMBER,' | \
      sed 's,"transcriptid",\//TRANSCRIPTID,'| sed '2d'
fi

if [ "${UNIQUIFY}" == "1" ];
then
    cat ${TAB_FILE} | \
        # Skip first line
        sed -e '1,2d' | \
        # Look at third and fourth columns
        awk {'print $3" "$4'} | \
        # sort on second column
        sort -k 2,2 | \
        # uniq on second column
        uniq -cf 1 | \
        # re-sort by counts (not really necessary)
        sort -rn | \
        # Use read counts as node size
        awk {'print "//NODESIZE\t\""$1"-depth-read-"$2"\"\t"$1'}
fi
