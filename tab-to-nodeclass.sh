#! /bin/bash

while getopts eut: ARG
do
  case ${ARG} in
    (e) EXONS="1";;
    (u) UNIQUIFY="1";;
    (t) TAB_FILE="${OPTARG}";;
  esac
done

if [ -z "${TAB_FILE}" ];
then
    echo "tab file not supplied"
    exit 1
fi

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
BASE_NAME=$(basename "$TAB_FILE")
NO_EXT="${BASE_NAME%.*}"
COLOUR="${SCRIPT_DIR}/color.txt"

if [ "${UNIQUIFY}" == "1" ];
then
    #create nodeclass
    cat ${TAB_FILE} | \
        sed -e '1,2d' | \
	cut -f3,4 | \
	sort -k1 -u | \
	sort -k2 | \
	uniq -cf1| \
	sort -rn > t1

    cat ${TAB_FILE} | \
        sed -e '1,2d' | \
        awk {'print $3"\t "$6"\t" $7"\t "$4'} | \
        sort -k4 > t2

    awk 'FNR==NR {C[$2]=$1;next}FNR==1 \
        {print "Count Read_ID Sequence Exon Transcript_ID"; next}$1 in C \
        {print C[$1], $1, $4, $3, $2}' t1 t2 > t3

    cat t3 | awk '{print "//NODECLASS\t\"" $2"_"$1  "\"\t\"Exon " $4 "\"\t\"" $5 "\""}'

    #create nodesize
    cat t1 | \
	awk {'print "//NODESIZE\t\""$2"_"$1"\"\t"$1'} 

    rm t1 t2 t3
    
    #create nodeclasscolor
    cat ${TAB_FILE} | \
         sed -e '1,2d' | \
         cut -f7,6 | \
         sort -u | \
         sort | \
         awk '{print "//NODECLASSCOLOR\t\t\"Exon " $2 "\"\t\"" $1 "\""}' | \
         sort > t1    

    awk 'FNR==NR{B=$NF;$NF="";gsub(/[[:space:]]+$/,X,$0);A[$0]=B;next} \
        ($2" "$3 in A){print $0 OFS A[$2" "$3]}' ${COLOR} \
        OFS="\t" t1
   
    rm t1

else
        awk '{print "//NODECLASS\t\"" $3 "\"\t\"Exon " $7 "\"\t\"" $6 "\""}' ${TAB_FILE} | \
        sed 's,"readname",\//GENEID,' | sed 's,"Exon exonnumber",\//EXONNUMBER,' | \
        sed 's,"transcriptid",\//TRANSCRIPTID,'| sed '2d'
	
	#create nodeclasscolor
        cat ${TAB_FILE} | \
         sed -e '1,2d' | \
         cut -f7,6 | \
         sort -u | \
         sort | \
         awk '{print "//NODECLASSCOLOR\t\t\"Exon " $2 "\"\t\"" $1 "\""}' | \
         sort > t1

	awk 'FNR==NR{B=$NF;$NF="";gsub(/[[:space:]]+$/,X,$0);A[$0]=B;next} \
        ($2" "$3 in A){print $0 OFS A[$2" "$3]}' ${COLOR} \
        OFS="\t" t1

    rm t1
fi
