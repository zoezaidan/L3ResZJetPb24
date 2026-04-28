#!/bin/bash

# submit condor jobs from lxplus

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/l/lamartik/work/HIJEC/MCTruth2023PbPb/CMSSW_13_2_6/src

eval `scramv1 runtime -sh`

cd /afs/cern.ch/user/l/lamartik/work/HIJEC/L2resdiualanalysisstuff/redo_feb2026/residualanalysis/fillhistograms

cmsswpath=/afs/cern.ch/user/l/lamartik/work/HIJEC/MCTruth2023PbPb/CMSSW_13_2_6/src

cd ${cmsswpath}
cmsenv
cd -

# arguments             = $(era) $(outputfiletag) $(ismc) $(checkjetid) $(iszb) $(dol2res) $(dojer) $(fillforJER) $(jtptlimitforalpha) $(outputpath) $(jesvariation) $(jervariation)
# given in the submission file

era=$1
outputfiletag=$2
ismc=$3
checkjetid=$4
iszb=$5
dol2res=$6
dojer=$7
fillforjer=$8
jtptlimitforalpha=$9
outputpath=${10}
jesvar=${11}
jersfvar=${12}


echo "outputpath is" $outputpath

resultpath=updated_MCtruth

if [ -d $outputpath ]; then
  echo "Given output directory already exists"
  echo "Given output path was "$outputpath
else
  mkdir $outputpath  
  echo "Created output directory "$outputpath
fi


fullpath=$outputpath'/'$resultpath
if [ -d $fullpath ]; then
  echo "Given output directory already exists"
  echo "Full output path was "$fullpath
else
  mkdir $fullpath
  echo "Created output directory "$fullpath
fi


root -l -x  'analyse.cc('\"${era}\",\"${outputfiletag}\",$ismc,$checkjetid,$iszb,$dol2res,$dojer,$fillforjer,$jtptlimitforalpha,\"${fullpath}\",\"${jesvar}\",\"${jersfvar}\"')'


