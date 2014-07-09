#! /bin/bash

while getopts ut: ARG
do
  case ${ARG} in
    (u) UNIQUIFY="1";;
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
        # Change to fasta format
        awk {'print ">"$2" "$1"\n"$3'}
else
    cat ${TAB_FILE} | \
        # Skip first line
        sed -e '1,2d' | \
        # Change to fasta format
        awk {'print ">"$3"\n"$4'}
fi
