#!/bin/bash
set -e

TAG=$1

echo "Running job: $TAG"

SCRATCH=${_CONDOR_SCRATCH_DIR:-$(pwd)}

# Set up CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd /eos/home-z/zzaidanc/CMSSW_15_1_0_patch3/src
eval $(scram runtime -sh)
cd $SCRATCH

# Set X509 proxy for xrootd access
export X509_USER_PROXY=/eos/home-z/zzaidanc/x509up_zzaidanc

# Detect if MC
if [[ "$TAG" == *"MC"* ]]; then
    IS_MC="true"
else
    IS_MC="false"
fi
echo "isMC: $IS_MC"

# Run the analysis (era mode: file path comes from configurations.h)
root -l -b -q "analyse_ZJet.cc+(\"${TAG}\", \"AK4_zjet\", ${IS_MC}, false, \"era\", -1, -1, \".\")"

# Copy output to EOS
OUTFILE="${TAG}_AK4_zjet_ak4.root"
if [ ! -f "$OUTFILE" ]; then
    echo "ERROR: Output file ${OUTFILE} not found. Files present:"
    ls *.root 2>/dev/null || echo "  none"
    exit 1
fi

EOSDIR=root://eosuser.cern.ch//eos/home-z/zzaidanc/L3ResZJetpp24/Outputs/
echo "Copying ${OUTFILE} to EOS..."
xrdcp -f ${OUTFILE} ${EOSDIR}/${OUTFILE}

echo "Done."
