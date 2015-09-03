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
	# Cut column read and sequence
	cut -f3,4 | \
	# Sort first column and remove duplicate
	sort -k1 -u | \
	# Sort second column and remove duplicate read name
	sort -k2 | \
	# Count occurrence of duplicate sequence
	uniq -cf1 | \
	# Sort (not really neccessary)
	sort -rn| \
	awk {'print ">"$2"_"$1" "$1"\n"$3'}
else
    cat ${TAB_FILE} | \
        # Skip first line
        sed -e '1,2d' | \
	# Cut column read and sequence
	cut -f3,4 | \
	# Sort first column and remove duplicate
	sort -k1 -u | \
	# Change to fasta format
	awk {'print ">"$1"\n"$2'}
fi
