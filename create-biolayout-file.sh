#! /bin/bash

SCRIPT_NAME=$(basename $0)

while getopts b:t:g:o:d:n:h ARG
do
  case ${ARG} in
    (b) UNSORTED_BAM_FILE=$(readlink -f "$OPTARG");;
    (t) CHROMOSOME_LENGTH_FILE=$(readlink -f "$OPTARG");;
    (g) GTF_FILE=$(readlink -f "$OPTARG");;
    (o) OUTPUT_DIRECTORY=$(readlink -f "$OPTARG");;
    (d) GENE_LIST="$OPTARG";;
    (n) GENE_LIST=$(cat "$OPTARG");;
    (h)
      echo "${SCRIPT_NAME}"
      echo "  One of -d or -n must be specified in addition to all the other options"
      echo "  -b <file> The unsorted BAM file"
      echo "  -t <file> The chromosome length file"
      echo "  -g <file> The GTF file"
      echo "  -o <directory> The directory in which to place the output"
      echo "  -d \"GENE1, GENE2, ..., GENEN\" A list of genes to examine"
      echo "  -n <file> A file containing a list of genes to examine"
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

if [ ! "${GENE_LIST}" ];
then
  echo "Gene list is empty; use -d or -n options"
  exit 1
fi

GENE_LIST=$(echo "$GENE_LIST" | perl -pe 's/\s*,\s*|\s+/ /g')

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
if [ "$?" != 0 ];
then
  echo "Cannot create ${OUTPUT_DIRECTORY}"
  exit $?
fi

COMPUTE_DATA=0
HASH_FILE="${OUTPUT_DIRECTORY}/hash"
if [ ! -e "${HASH_FILE}" ];
then
  COMPUTE_DATA=1
else
  EXISTING_HASH=$(cat ${HASH_FILE})
  if [ "${EXISTING_HASH}" != "${INPUT_HASH}" ];
  then
    COMPUTE_DATA=1
  fi
fi

BASENAME_BAM_FILE=$(basename ${UNSORTED_BAM_FILE})
NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"
SAMTOOLS_SORTED_BAM_FILE="${OUTPUT_DIRECTORY}/${NO_EXT_BAM_FILE}-sorted"
SORTED_BAM_FILE="${SAMTOOLS_SORTED_BAM_FILE}.bam"
GRANGES_FILE="${OUTPUT_DIRECTORY}/${NO_EXT_BAM_FILE}-GRanges.RData"
GTF_ANNOTATION_FILE="${OUTPUT_DIRECTORY}/${NO_EXT_BAM_FILE}-ensembl_gtfannotation.RData"

if [ "${COMPUTE_DATA}" == "1" ];
then
  echo "Sorting ${UNSORTED_BAM_FILE}..."
  ${SAMTOOLS} sort ${UNSORTED_BAM_FILE} ${SAMTOOLS_SORTED_BAM_FILE}
  if [ "$?" != 0 ];
  then
    echo "samtools sort step failed"
    exit $?
  fi

  ${R_SCRIPT} grangesscript.R -b "${SORTED_BAM_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GRANGES_FILE}"
  if [ "$?" != 0 ];
  then
    echo "grangesscript.R failed"
    exit $?
  fi

  ${R_SCRIPT} grangesscript_gtf.R -g "${GTF_FILE}" -t "${CHROMOSOME_LENGTH_FILE}" \
    -o "${GTF_ANNOTATION_FILE}"
  if [ "$?" != 0 ];
  then
    echo "grangesscript_gtf.R failed"
    exit $?
  fi

  ${R_SCRIPT} findoverlaps.R -g "${GRANGES_FILE}" -e "${GTF_ANNOTATION_FILE}" -d "${GENE_LIST}" \
    -p "${OUTPUT_DIRECTORY}/"
  if [ "$?" != 0 ];
  then
    echo "findoverlaps.R failed"
    exit $?
  fi
fi

echo ${INPUT_HASH} > ${HASH_FILE}

R2R_OUTPUT_DIR="${OUTPUT_DIRECTORY}/r2r_output"

for GENE in ${GENE_LIST}
do
  echo "Writing ${GENE}.fasta"
  ./tab-to-fasta.sh "${OUTPUT_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.fasta"

  echo "Writing ${GENE}.nodeclass"
  ./tab-to-nodeclass.sh "${OUTPUT_DIRECTORY}/${GENE}.tab" > \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass"

  rm -rf "${R2R_OUTPUT_DIR}"
  ./read2read.py "${OUTPUT_DIRECTORY}/${GENE}.fasta" "${R2R_OUTPUT_DIR}"
  if [ "$?" != 0 ];
  then
    echo "read2read.py failed"
    exit $?
  fi
  
  cat "${R2R_OUTPUT_DIR}/${GENE}_pairwise.txt" \
    "${OUTPUT_DIRECTORY}/${GENE}.nodeclass" > \
    "${OUTPUT_DIRECTORY}/${GENE}.layout"
done
