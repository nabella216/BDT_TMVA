#!/bin/bash

##
trainlabel=20240701init

# -- signal sample for training
inputs=/afs/cern.ch/work/w/wangj/public/summer2024/mcfilelist_OfficialMCPthat5And8_GammaPlusGoing__add_zdcanalyzer-zdcdigi_sumPlus_sumMinus.root

# -- background sample for training
bkgstrategy=sideband
inputb=/afs/cern.ch/work/w/wangj/public/summer2024/filelist_May07__add_Dfinder-ntDkpi_Dgen.root

# -- mva application sample
inputms=(
    /afs/cern.ch/work/w/wangj/public/summer2024/filelist_May07__add_Dfinder-ntDkpi_Dgen.root
)
outputmvadir=/home/data/public/wangj/UPC/mva_output_${trainlabel}/

##
# -- PV + MET
cut="pprimaryVertexFilter && vz > -15 && vz < 15 && nVtx <= 1 && cscTightHalo2015Filter > 0"
# -- track kinematics
cut=$cut" && Dtrk1Pt > 0.9 && Dtrk2Pt > 0.9 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4"
# -- track qualirty
cut=$cut" && Dtrk1highPurity && Dtrk2highPurity && TMath::Abs(Dtrk1PtErr/Dtrk1Pt) < 0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt) < 0.1"
# -- X prefilter
cut=$cut" && TMath::Abs(Dy) < 2.4 && Dchi2cl > 0.1"
# -- tricky selection
cut=$cut" && DsvpvDisErr>1.e-5 && DsvpvDisErr_2D>1.e-5"

cuts=$cut
cutb=$cut
# cutb=$cutb" && HLT_HIUPC_SingleJet8_ZDC1nXOR_MaxPixelCluster50000_v2"
# cutb=$cutb" && HLT_HIUPC_SingleJet8_NotMBHF2AND_MaxPixelCluster50000_v2"
# cutb=$cutb" && ((sumPlus > 1100 && sumMinus <= 1000 && rapidityGapMinus) || (sumPlus <= 1100 && sumMinus > 1000 && rapidityGapPlus))"

##
# algo="BDT,BDTG,LD,DNN_GPU"
#algo="BDT,BDTG,CutsGA,LD"
algo = "BDT, BDTG, BDTF, BDTD, BDTB"
#algo="BDT"

# stages="0,1,3,2,4,5,6,7,10,11,12,13" ; sequence=0 ; # see definition below #
stages="0,1,2,5,6,7,12,13" ; sequence=1 ; # see definition below #

## ===== do not change the lines below =====
varstrategy=("Single set" "Sequence")

cuts=$cuts" && (Dgen == 23333 || Dgen == 41022 || Dgen == 41044)"
# -- HLT
# cutb=$cutb" && HLT_HIUPC_SingleJet8_ZDC1nAsymXOR_MaxPixelCluster50000_v1"
[[ $bkgstrategy == "sideband" ]] && cutb=$cutb" && TMath::Abs(Dmass-1.8648) > 0.05 && TMath::Abs(Dmass-1.8648) < 0.12" # sideband
[[ $bkgstrategy == "samesign" ]] && cutb=$cutb" && TMath::Abs(Dmass-1.8648) < 0.03" # samesign

## ===== do not do not do not change the lines below =====
function catspace() { echo -e $(cat "$@" | sed  's/$/\\n/' | sed 's/ /\\a /g') ; }
IFS=','; allstages=($stages); unset IFS;
echo -e '
##########################
# Variables \e[32m(To be used)\e[0m #
##########################'
vv=0
catspace TMVAClassification.h | grep --color=no 'mytmva::tmvavar("' > tmpvarlist.list
while read -r line
do
    for ss in ${allstages[@]} ; do [[ $ss == $vv ]] && { echo -en "\e[32m" ; break ; } ; done ;
    echo -e $line | sed 's/mytmva::tmvavar//' | sed 's/\*\///' | sed 's/\a \a \/\*//' ; echo -ne "\e[0m" ;
    vv=$((vv+1))
done < tmpvarlist.list
rm tmpvarlist.list

##
echo -e "
###########################
# Training Configurations #
###########################

>>>>>> Variables training strategy
  >>>> \e[32m[${varstrategy[sequence]}]\e[0m

>>>>>> Background strategy
  >>>> \e[32m[$bkgstrategy]\e[0m

>>>>>> Algorithms
  >>>> \e[32m[$algo]\e[0m

>>>>>> Input files
  >>>> Signal:      \e[32m$inputs\e[0m
  >>>> Background:  \e[32m$inputb\e[0m

>>>>>> Selections
  >>>> Prefilters
    >> \e[32m\"$cut\"\e[0m
  >>>> Signal cut
    >> \e[32m\"${cuts##$cut}\"\e[0m
  >>>> Background cut
    >> \e[32m\"${cutb##$cut}\"\e[0m
"

##
output=rootfiles/TMVA_trainD0_${trainlabel}_${bkgstrategy}
[[ -d $output ]] && rm -r $output
[[ ${3:-0} -eq 1 ]] && mkdir -p $outputmvadir
tmp=$(date +%y%m%d%H%M%S)

##
[[ $# -eq 0 ]] && echo "usage: ./run_TMVAClassification.sh [train] [draw curves] [create BDT tree]"
echo "Compiling .cc macros..."

echo -e "\e[35m==> (1/5) building TMVAClassification.C\e[0m"
make TMVAClassification || exit 1
echo -e "\e[35m==> (2/5) building guivariables.C\e[0m"
make guivariables || exit 1
echo -e "\e[35m==> (3/5) building guiefficiencies.C\e[0m"
make guiefficiencies || exit 1
echo -e "\e[35m==> (4/5) building guieffvar.C\e[0m"
make guieffvar || exit 1
echo -e "\e[35m==> (5/5) building mvaprod.C\e[0m"
make mvaprod || exit 1

for i in TMVAClassification guivariables guiefficiencies guieffvar mvaprod ; do
    cp $i ${i}_${tmp}.exe
done

[[ ${1:-0} -eq 1 ]] && {
    conf=
    echo -e "\e[2m==> Do you really want to run\e[0m \e[1mTMVAClassification.C\e[0m \e[2m(it might be very slow)?\e[0m [y/n]"
    read conf
    while [[ $conf != 'y' && $conf != 'n' ]] ; do { echo "warning: input [y/n]" ; read conf ; } ; done ;
    [[ $conf == 'n' ]] && { rm *_${tmp}.exe ; exit ; }
}
# SEE HERE
# Run TMVAClassification.C within ROOT context for each NTrees value
#for ntrees in "${ntrees_values[@]}"; do
#    methodName="BDT_NTrees${ntrees}"
#    output_n="${output}_${methodName}"
# train
stage=$stages
while [[ $stage == *,* ]]
do
 [[ ${1:-0} -eq 1 ]] && { ./TMVAClassification_${tmp}.exe $inputs $inputb "$cuts" "$cutb" $output "$algo" "$stage"; } 
# [[ ${1:-0} -eq 1 ]] && { ./TMVAClassification_${tmp}.exe $inputs $inputb "$cuts" "$cutb" $output "$algo" "$stage" "$ntrees; } 
    [[ $sequence -eq 0 ]] && break;
    while [[ $stage != *, ]] ; do stage=${stage%%[0-9]} ; done ;
    stage=${stage%%,}
done

# draw curves
[[ ${2:-0} -eq 1 ]] && { 
    ./guivariables_${tmp}.exe $output "$algo" "$stages"
    ./guiefficiencies_${tmp}.exe $output "$algo" "$stages"
}
# draw curve vs. var
[[ ${2:-0} -eq 1 && $sequence -eq 1 ]] && ./guieffvar_${tmp}.exe $output "$algo" "$stages"

# produce mva values
for inputm in ${inputms[@]} ; do 
    [[ ${3:-0} -eq 1 ]] && ./mvaprod_${tmp}.exe $inputm "Dfinder/ntDkpi" $output $outputmvadir "$algo" "${stages}"
done

##
rm *_${tmp}.exe
