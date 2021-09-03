#root -b -q -l "MakeDoAllFitsScript.C(collId, whichModel, PHIBINS, INTBIN, PTBINS, YBINS, CBINS)" 

#pp bins
root -b -q -l "MakeDoAllFitsScript.C(kPPDATA, 0, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE)"
./doAllFits.sh

#PbPb bins
root -b -q -l "MakeDoAllFitsScript.C(kAADATA, 0, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE)"
./doAllFits.sh

#PbPb dphi bins
root -b -q -l "MakeDoAllFitsScript.C(kAADATA, 0, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE)"
./doAllFits.sh
