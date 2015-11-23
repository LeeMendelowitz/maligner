#!/usr/bin/env bash

# This example will align E. coli K-12 assembled sequence contigs (contigs.fasta)
# to an in silico digest of the E. coli K-12 reference (ecoli_k12_ref.fasta)

####################################################################
# CONFIGURE PATH & PYTHON PATH

# set to the cmake build directory
MALIGNER_INSTALL_PATH=~/workspace/maligner/build

# set PATH and PYTHONPATH so maligner binaries and utility 
# scripts are in the bath
export PATH=${MALIGNER_INSTALL_PATH}/bin:$PATH
export PYTHONPATH=${MALIGNER_INSTALL_PATH}/lib:$PYTHONPATH


###################################################################
# DEFINE VARIABLES CONTROLLING INPUT FILES, OUTPUT FILES,
# AND SETTINGS

# inputs
REF_FASTA=ecoli_k12_ref.fasta
CONTIGS_FASTA=contigs.fasta

# outputs
REF_OUT_PFX=ecoli_k12_ref.BamHI
MAPS_OUT_PFX=contigs.BamHI
OUT_PFX=contigs_to_ref

# setttings
REC_SEQ=GGATCC # recognition sequence for BamHI
MIN_FRAG_SIZE=1000 # bp units
QUERY_MISS_PENALTY=3.0
REF_MISS_PENALTY=3.0
QUERY_MAX_MISSES=5
REF_MAX_MISSES=5
SD_RATE=0.05
MIN_SD=750
MAX_SCORE_PER_INNER_CHUNK=1.0
MAX_ALIGNMENTS_PER_QUERY=5


###################################################################
# PREPARE MALIGNER_DP INPUTS

# convert the contigs.fasta file to the Maligner maps format
echo -e 1>&2 "\n\n*********************************************"
echo -e 1>&2 "Converting contigs to Maligner maps file.\n\n"
make_insilico_map -o $MAPS_OUT_PFX $CONTIGS_FASTA $REC_SEQ

# smooth the maps file by merging consecutive fragments that
# are less than 1kb
echo -e 1>&2 "\n\n*********************************************"
echo -e 1>&2 "Smoothing small fragments from contigs maps file.\n\n"
smooth_maps_file -m $MIN_FRAG_SIZE ${MAPS_OUT_PFX}.maps > ${MAPS_OUT_PFX}.smoothed.maps

# convert the reference fasta file to the Maligner maps format
echo -e 1>&2 "\n\n*********************************************"
echo -e 1>&2 "Converting E. coli K-12 reference to Maligner maps file.\n\n"
make_insilico_map -o $REF_OUT_PFX $REF_FASTA $REC_SEQ

echo -e 1>&2 "\n\n*********************************************"
echo -e 1>&2 "Smoothing small fragments from reference maps file.\n\n"
smooth_maps_file -m $MIN_FRAG_SIZE ${REF_OUT_PFX}.maps > ${REF_OUT_PFX}.smoothed.maps

###################################################################
# RUN MALIGNER_DP

# Align the smoothed query maps file to the smoothed
# reference maps file with maligner_dp.
# 
# The alignments are saved to $OUT_PFX.aln,
# Log file written to $OUT_PFX.

echo -e 1>&2 "\n\n*********************************************"
echo -e 1>&2 "Running maligner_dp.\n\n"
maligner_dp \
  -q $QUERY_MISS_PENALTY \
  -r $REF_MISS_PENALTY \
  --query-max-misses $QUERY_MAX_MISSES \
  --ref-max-misses $REF_MAX_MISSES \
  --max-score-per-inner-chunk $MAX_SCORE_PER_INNER_CHUNK \
  --sd-rate $SD_RATE \
  --min-sd $MIN_SD \
  --reference-is-circular \
  --no-query-rescaling \
  --max-alignments $MAX_ALIGNMENTS_PER_QUERY \
  ${MAPS_OUT_PFX}.smoothed.maps \
  ${REF_OUT_PFX}.smoothed.maps \
  2>&1 1> ${OUT_PFX}.aln| tee ${OUT_PFX}.log


