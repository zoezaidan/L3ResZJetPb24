#!/bin/bash
# Condor job wrapper for Z+Jet analysis
# Arguments: INPUT_KEY IS_MC OUTPUT_TAG BATCH_INDEX TOTAL_BATCHES

INPUT_KEY=$1
IS_MC=$2
OUTPUT_TAG=$3
BATCH_INDEX=${4:--1}
TOTAL_BATCHES=${5:-1}

# Redirect all output to EOS log (keeps EOS paths out of the condor submit file)
LOG_DIR="/eos/home-z/zzaidanc/L3ResZJetpp24/fillhistograms/logs"
exec >> "${LOG_DIR}/zjet_${INPUT_KEY}_batch${BATCH_INDEX}.log" 2>&1

echo "Starting job at $(date)"
echo "Input key    : $INPUT_KEY"
echo "isMC         : $IS_MC"
echo "Output tag   : $OUTPUT_TAG"
echo "Batch index  : $BATCH_INDEX / $TOTAL_BATCHES"

# ---- Set your paths here ----
CMSSW_DIR="/eos/home-z/zzaidanc/CMSSW_15_1_0_patch3/src"
ANALYSIS_DIR="/eos/home-z/zzaidanc/L3ResZJetpp24/fillhistograms"
# -----------------------------

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_DIR
eval `scramv1 runtime -sh`

cd $ANALYSIS_DIR

root -l -b -q "analyse_ZJet.cc(\"${INPUT_KEY}\", \"${OUTPUT_TAG}\", ${IS_MC}, false, \"era\", -1, -1, \"\", ${BATCH_INDEX}, ${TOTAL_BATCHES})"

echo "Job finished at $(date)"
