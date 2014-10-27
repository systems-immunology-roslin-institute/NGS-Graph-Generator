#! /bin/bash

SCRIPT_NAME=$(basename $0)
DIR_NAME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo_timestamp()
{
  date +"%R:%S $*"
}

while getopts b:t:g:o:c:d:n:p:l:uh ARG
do
  case ${ARG} in
    (b) UNSORTED_BAM_FILE=$(readlink -f "$OPTARG");;
    (t) CHROMOSOME_LENGTH_FILE=$(readlink -f "$OPTARG");;
    (g) GTF_FILE=$(readlink -f "$OPTARG");;
    (o) OUTPUT_DIRECTORY=$(readlink -f "$OPTARG");;
    (c) CACHE_DIRECTORY=$(readlink -f "$OPTARG");;
    (d) GENE_LIST="$OPTARG";;
    (n) GENE_LIST=$(cat "$OPTARG");;
    (p) PERCENTAGE="$OPTARG";;
    (l) COVERAGE="$OPTARG";;
    (u) UNIQUIFY="-u";;
    (h)
      echo "${SCRIPT_NAME}"
      echo "  One of -d or -n must be specified in addition to all the other options"
      echo "  -b <file> The unsorted BAM file"
      echo "  -t <file> The chromosome length file"
      echo "  -g <file> The GTF file"
      echo "  -o <directory> The directory in which to place the output"
      echo "  -c <directory> The directory in which to cache intermediate files"
      echo "  -d \"GENE1, GENE2, ..., GENEN\" A list of genes to examine"
      echo "  -n <file> A file containing a list of genes to examine"
      echo "  -p <value> The percentage similarity value (default 85)"
      echo "  -l <value> The percentage coverage value (default 55)"
      echo "  -u discard redundant reads"
      exit 0;;

    (*)
      exit 1;;
  esac
done

if [ ! "${UNSORTED_BAM_FILE}" ] || [ ! -e "${UNSORTED_BAM_FILE}" ];
then
  echo "BAM file not supplied or not found; use -b option"
  exit 1;
fi

if [ ! "${CHROMOSOME_LENGTH_FILE}" ] || [ ! -e "${CHROMOSOME_LENGTH_FILE}" ];
then
  echo "Chromosome length file not supplied or not found; use -t option"
  exit 1;
fi

if [ ! "${GTF_FILE}" ] || [ ! -e "${GTF_FILE}" ];
then
  echo "GTF file not supplied or not found; use -g option"
  exit 1;
fi

if [ ! "${OUTPUT_DIRECTORY}" ];
then
  OUTPUT_DIRECTORY=$(readlink -f output)
fi

if [ ! "${CACHE_DIRECTORY}" ];
then
  CACHE_DIRECTORY=$(readlink -f cache)
fi

if [ ! "${GENE_LIST}" ];
then
  echo "Gene list is empty; use -d or -n options"
  exit 1
fi

if [ ! "${PERCENTAGE}" ];
then
  PERCENTAGE="85"
fi

if [ ! "${COVERAGE}" ];
then
  COVERAGE="55"
fi

mkdir -p ${OUTPUT_DIRECTORY}
EXITCODE="$?"
if [ "$EXITCODE" != 0 ];
then
  echo_timestamp "Cannot create ${OUTPUT_DIRECTORY}"
  exit $EXITCODE
fi

echo_timestamp "Computing valid gene names..."
VALID_GENES=$(perl -pe 's/.*gene_name "([^"]+)".*/\1/' ${GTF_FILE} | \
  tr '[:lower:]' '[:upper:]' | sort | uniq)
GENE_LIST=$(echo "$GENE_LIST" | tr '[:lower:]' '[:upper:]' | perl -pe 's/\s*,\s*|\s+/ /g')

for GENE in ${GENE_LIST}
do
  echo ${VALID_GENES} | grep -q "\(^\|.*\s\+\)${GENE}\($\|.*\s\+\)"
  if [ "$?" != 0 ];
  then
    echo_timestamp "Gene ${GENE} is not present in ${GTF_FILE}."
    exit 1
  fi
done

NUM_CORES=$(nproc)
echo_timestamp "Using ${NUM_CORES} cores..."

echo_timestamp "Computing hash..."
INPUT_HASH=$(cat ${UNSORTED_BAM_FILE} ${CHROMOSOME_LENGTH_FILE} ${GTF_FILE} | md5sum)
INPUT_HASH=${INPUT_HASH%% *}

echo_timestamp "BAM file: ${UNSORTED_BAM_FILE}"
echo_timestamp "Chromosome length file: ${CHROMOSOME_LENGTH_FILE}"
echo_timestamp "GTF file: ${GTF_FILE}"
echo_timestamp "Output directory: ${OUTPUT_DIRECTORY}"
echo_timestamp "Gene list: ${GENE_LIST}"
echo_timestamp "Input hash: ${INPUT_HASH}"

SAMTOOLS=samtools
R_SCRIPT=Rscript

COMPUTE_DATA=0
HASH_DIRECTORY="${CACHE_DIRECTORY}/${INPUT_HASH}"
CACHE_LOCK="${HASH_DIRECTORY}/lock"

takeCacheLock()
{
  if ( set -o noclobber; echo_timestamp "$$" > "$CACHE_LOCK") 2> /dev/null;
  then
    trap 'rm -f "$CACHE_LOCK"; exit $?' INT TERM EXIT

    echo_timestamp "Aquired lock ${CACHE_LOCK}..."
    return 1
  else
    return 0
  fi
}

takeCacheLockBlocking()
{
  echo_timestamp "Waiting for ${CACHE_LOCK}..."
  takeCacheLock
  while (( "$?" == "0" ))
  do
    sleep 30
    takeCacheLock
  done
}

waitForCacheToComplete()
{
  takeCacheLockBlocking
  releaseCacheLock
}

releaseCacheLock()
{
  rm -f "$CACHE_LOCK"
  trap - INT TERM EXIT
  echo_timestamp "Released lock ${CACHE_LOCK}..."
}

cacheIsLocked()
{
  takeCacheLock
  if [ "$?" == "1" ];
  then
    return 1
    releaseCacheLock
  else
    return 0
  fi
}

echo_timestamp Checking if ${HASH_DIRECTORY} exists...
if [ ! -d "${HASH_DIRECTORY}" ];
then
  COMPUTE_DATA=1
  mkdir -p ${HASH_DIRECTORY}
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo_timestamp "Cannot create ${HASH_DIRECTORY}"
    exit $EXITCODE
  fi
  takeCacheLock
else
  if [ cacheIsLocked ];
  then
    waitForCacheToComplete
  fi

  if [ ! -e "${HASH_DIRECTORY}/valid" ]
  then
    # There is a possible race condition where two or more processes
    # manage to get here simultaneously, but the chances of it happening
    # are fairly small
    COMPUTE_DATA=1
    takeCacheLock
    if [ "$?" == "0" ];
    then
      echo_timestamp "!!! Probable race condition encountered, exiting."
      exit 1
    fi
  fi
fi

BASENAME_BAM_FILE=$(basename ${UNSORTED_BAM_FILE})
NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"
SAMTOOLS_SORTED_BAM_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-sorted"
SORTED_BAM_FILE="${SAMTOOLS_SORTED_BAM_FILE}.bam"
GRANGES_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-GRanges.RData"
GTF_ANNOTATION_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-ensembl_gtfannotation.RData"

if [ "${COMPUTE_DATA}" == "1" ];
then
  echo_timestamp "Sorting ${UNSORTED_BAM_FILE}..."
  ${SAMTOOLS} sort ${UNSORTED_BAM_FILE} ${SAMTOOLS_SORTED_BAM_FILE}
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo_timestamp "samtools sort step failed"
    exit $EXITCODE
  fi

  ${R_SCRIPT} ${DIR_NAME}/grangesscript.R -b "${SORTED_BAM_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GRANGES_FILE}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo_timestamp "grangesscript.R failed"
    exit $EXITCODE
  fi

  ${R_SCRIPT} ${DIR_NAME}/grangesscript_gtf.R -g "${GTF_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GTF_ANNOTATION_FILE}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo_timestamp "grangesscript_gtf.R failed"
    exit $EXITCODE
  fi

  touch "${HASH_DIRECTORY}/valid"
  releaseCacheLock
else
  echo_timestamp "Using cached data..."
fi

# Make a symlink in the job directory to the cached data
ln -s ${HASH_DIRECTORY} ${OUTPUT_DIRECTORY}/${INPUT_HASH}

R2R_OUTPUT_DIR="${OUTPUT_DIRECTORY}/r2r_output"

for GENE in ${GENE_LIST}
do
  if [ ! -e "${HASH_DIRECTORY}/${GENE}.tab" ];
  then
    takeCacheLockBlocking
    echo_timestamp "Writing ${GENE}.tab"
    ${R_SCRIPT} ${DIR_NAME}/findoverlaps.R -g "${GRANGES_FILE}" -e "${GTF_ANNOTATION_FILE}" -d "${GENE}" \
      -p "${HASH_DIRECTORY}/"
    releaseCacheLock
    EXITCODE="$?"
    if [ "$EXITCODE" != 0 ];
    then
      echo_timestamp "findoverlaps.R failed"
      exit $EXITCODE
    fi
  else
    echo_timestamp "Using cached ${GENE}.tab"
  fi

  echo_timestamp "Writing ${GENE}.fasta"
  ${DIR_NAME}/tab-to-fasta.sh ${UNIQUIFY} -t "${HASH_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.fasta"

  echo_timestamp "Writing ${GENE}.nodeclass"
  ${DIR_NAME}/tab-to-nodeclass.sh ${UNIQUIFY} -e -t "${HASH_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass"

  rm -rf "${R2R_OUTPUT_DIR}"
  ${DIR_NAME}/read2read.py -p ${PERCENTAGE} -l ${COVERAGE} -a ${NUM_CORES} \
    "${OUTPUT_DIRECTORY}/${GENE}.fasta" "${R2R_OUTPUT_DIR}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo_timestamp "read2read.py failed"
    exit $EXITCODE
  fi
  
  OUTPUT_FILE_NAME="${GENE}-s${PERCENTAGE}-c${COVERAGE}${UNIQUIFY}.layout"
  cat "${R2R_OUTPUT_DIR}/${GENE}_pairwise.txt" \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass" > \
    "${OUTPUT_DIRECTORY}/${OUTPUT_FILE_NAME}"
done

BASE_DIRECTORY_NAME=$(basename ${OUTPUT_DIRECTORY})
zip -j ${OUTPUT_DIRECTORY}/${BASE_DIRECTORY_NAME}.zip ${OUTPUT_DIRECTORY}/*.layout
