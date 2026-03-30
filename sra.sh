#!/usr/bin/env bash

conda activate SRA 

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 PRJNAxxxxxx OUTPUT_DIR" >&2
    exit 1
fi

PRJ="$1"          #PRJNA identifier
OUTDIR="$2"       # output directory (choose PRJNA identifier or study number) 
SRA_LIST="$3"   # accession list available in accessions
JOBS="${JOBS:-4}" # parallel jobs for fasterq-dump 

mkdir -p "${OUTDIR}"

# Loop over SRR accessions and download FASTQ using fasterq-dump
while read -r SRR; do
   
   [[ -z "${SRR}" ]] && continue
   echo "Downloading ${SRR} ..."
   fasterq-dump --split-files --outdir "${OUTDIR}" "${SRR}"
done < "${SRA_LIST}"
       
