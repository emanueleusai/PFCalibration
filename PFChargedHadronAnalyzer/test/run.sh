#!/bin/bash

echo $SHELL

INPUTFILE=$1
OUTPUTFILE=$2
OUTPUTDIR=$3
CARDFILE=$4

source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW_9_4_1
cd CMSSW_9_4_1/src
eval `scram runtime -sh`
cd ../../
tar -xzf delphes.tar.gz
cd delphes
cp ../${CARDFILE} cards/CMS_PhaseII/
xrdcp root://cmseos.fnal.gov//store/user/eusai/MinBias_100k.pileup .
xrdcp root://cmseos.fnal.gov//store/user/eusai/${INPUTFILE} .
./DelphesCMSFWLite cards/CMS_PhaseII/${CARDFILE} ${OUTPUTFILE} ${INPUTFILE}
xrdcp -f ${OUTPUTFILE} ${OUTPUTDIR}/${OUTPUTFILE}


# #
# # variables from arguments string in jdl
# #

# PILEUPFILE=$1
# INPUTFILE=$2
# OUTPUTFILE=$3
# OUTDIR=$4
# echo ""
# echo "PILEUPFILE $PILEUPFILE"
# echo "INPUTFILE:   $INPUTFILE"
# echo "OUTPUTFILE:  $OUTPUTFILE"
# echo "OUTDIR $OUTDIR"


# tar -xzf CMSSW_9_1_0_pre3.tar.gz
# cd CMSSW_9_1_0_pre3
# scram b ProjectRename
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# eval `scramv1 runtime -sh`
# cd -
# cd CMSSW_9_1_0_pre3/src/delphes
# echo $PWD
# xrdcp root://cmseos.fnal.gov//store/user/bmahakud/Delphes/MinBias_100k.pileup .


# #run macro

# #./DelphesCMSFWLite cards/CMS_PhaseII/CMS_PhaseII_0PU_v02.tcl   outtest1.root root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring17GS/RSGluonToTT_M-5000_TuneCUETP8M1_14TeV-pythia8/GEN-SIM/91X_upgrade2023_realistic_v3-v1/00000/00EC1E4F-C754-E711-B899-E0071B7A7870.root

# ./DelphesCMSFWLite cards/CMS_PhaseII/${PILEUPFILE}   ${OUTPUTFILE} ${INPUTFILE}


# CMSEXIT=$?

# if [[ $CMSEXIT -ne 0 ]]; then
#   echo "exit code $CMSEXIT, skipping xrdcp"
#   return $CMSEXIT
# else
# #echo "xrdcp output for condor"
#   for FILE in *.root
#     do
#       echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
#       xrdcp -f ${FILE} ${OUTDIR}/${FILE}
#       rm ${FILE}
#     done
# fi






