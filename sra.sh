#!/usr/bin/env bash

conda activate SRA 

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 PRJNAxxxxxx OUTPUT_DIR" >&2
    exit 1
fi

PRJ="$1"          # e.g. PRJNA257197
OUTDIR="$2"       # where FASTQs will go
SRA_LIST="$3"
JOBS="${JOBS:-4}" # parallel jobs for fasterq-dump

mkdir -p "${OUTDIR}"

# 2) Loop over SRR accessions and download FASTQ using fasterq-dump
while read -r SRR; do
   
   [[ -z "${SRR}" ]] && continue
   echo "Downloading ${SRR} ..."
   fasterq-dump --split-files --outdir "${OUTDIR}" "${SRR}"
done < "${SRA_LIST}"
       
