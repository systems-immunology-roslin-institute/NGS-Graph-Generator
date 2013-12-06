#! /bin/bash

SCRIPT_NAME=$(basename $0)
DIR_NAME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while getopts b:t:g:o:c:d:n:p:l:h ARG
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

GENE_LIST=$(echo "$GENE_LIST" | perl -pe 's/\s*,\s*|\s+/ /g')

NUM_CORES=$(nproc)
echo "Using ${NUM_CORES} cores..."

echo "Computing hash..."
INPUT_HASH=$(cat ${UNSORTED_BAM_FILE} ${CHROMOSOME_LENGTH_FILE} ${GTF_FILE} | md5sum)
INPUT_HASH=${INPUT_HASH%% *}

echo "BAM file: ${UNSORTED_BAM_FILE}"
echo "Chromosome length file: ${CHROMOSOME_LENGTH_FILE}"
echo "GTF file: ${GTF_FILE}"
echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Gene list: ${GENE_LIST}"
echo "Input hash: ${INPUT_HASH}"

SAMTOOLS=samtools
R_SCRIPT=Rscript

mkdir -p ${OUTPUT_DIRECTORY}
EXITCODE="$?"
if [ "$EXITCODE" != 0 ];
then
  echo "Cannot create ${OUTPUT_DIRECTORY}"
  exit $EXITCODE
fi

COMPUTE_DATA=0
HASH_DIRECTORY="${CACHE_DIRECTORY}/${INPUT_HASH}"
echo Checking if ${HASH_DIRECTORY} exists...
if [ ! -d "${HASH_DIRECTORY}" ];
then
  COMPUTE_DATA=1
  mkdir -p ${HASH_DIRECTORY}
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo "Cannot create ${HASH_DIRECTORY}"
    exit $EXITCODE
  fi
elif [ ! -e "${HASH_DIRECTORY}/valid" ]
then
  COMPUTE_DATA=1
fi

BASENAME_BAM_FILE=$(basename ${UNSORTED_BAM_FILE})
NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"
SAMTOOLS_SORTED_BAM_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-sorted"
SORTED_BAM_FILE="${SAMTOOLS_SORTED_BAM_FILE}.bam"
GRANGES_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-GRanges.RData"
GTF_ANNOTATION_FILE="${HASH_DIRECTORY}/${NO_EXT_BAM_FILE}-ensembl_gtfannotation.RData"

if [ "${COMPUTE_DATA}" == "1" ];
then
  echo "Sorting ${UNSORTED_BAM_FILE}..."
  ${SAMTOOLS} sort ${UNSORTED_BAM_FILE} ${SAMTOOLS_SORTED_BAM_FILE}
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo "samtools sort step failed"
    exit $EXITCODE
  fi

  ${R_SCRIPT} ${DIR_NAME}/grangesscript.R -b "${SORTED_BAM_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GRANGES_FILE}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo "grangesscript.R failed"
    exit $EXITCODE
  fi

  ${R_SCRIPT} ${DIR_NAME}/grangesscript_gtf.R -g "${GTF_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GTF_ANNOTATION_FILE}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo "grangesscript_gtf.R failed"
    exit $EXITCODE
  fi

  touch "${HASH_DIRECTORY}/valid"
else
  echo Using cached data...
fi

${R_SCRIPT} ${DIR_NAME}/findoverlaps.R -g "${GRANGES_FILE}" -e "${GTF_ANNOTATION_FILE}" -d "${GENE_LIST}" \
  -p "${OUTPUT_DIRECTORY}/"
EXITCODE="$?"
if [ "$EXITCODE" != 0 ];
then
  echo "findoverlaps.R failed"
  exit $EXITCODE
fi

R2R_OUTPUT_DIR="${OUTPUT_DIRECTORY}/r2r_output"

for GENE in ${GENE_LIST}
do
  echo "Writing ${GENE}.fasta"
  ${DIR_NAME}/tab-to-fasta.sh "${OUTPUT_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.fasta"

  echo "Writing ${GENE}.nodeclass"
  ${DIR_NAME}/tab-to-nodeclass.sh "${OUTPUT_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass"

  rm -rf "${R2R_OUTPUT_DIR}"
  ${DIR_NAME}/read2read.py -a ${NUM_CORES} -p ${PERCENTAGE} -l ${COVERAGE} \
    "${OUTPUT_DIRECTORY}/${GENE}.fasta" "${R2R_OUTPUT_DIR}"
  EXITCODE="$?"
  if [ "$EXITCODE" != 0 ];
  then
    echo "read2read.py failed"
    exit $EXITCODE
  fi
  
  cat "${R2R_OUTPUT_DIR}/${GENE}_pairwise.txt" \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass" > \
    "${OUTPUT_DIRECTORY}/${GENE}.layout"
done
