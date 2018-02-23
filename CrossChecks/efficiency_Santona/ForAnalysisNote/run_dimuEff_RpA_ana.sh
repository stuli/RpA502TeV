#!/bin/bash
#$1 = OniaMode (1,2,3), $2 = ispPb (1,0), $3 = pPb or pp, $4 = Name (eg. 29_5)
mkdir eff_$3$4
cp /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_RpA_ana.C /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/VVV/$1/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/WWW/$2/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/XXX/$3/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/RpA_ana/oniaMode$1_$3_$4/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/TAG/$4/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
cp /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/eff_$3$4/dimuEff_oniaMode$1_$3_$4.C
root -l -q -b /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
#This line is optional ^^^ You can run later. But if you do, don't remove it below!
rm /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/ForAnalysisNote/dimuEff_oniaMode$1_$3_$4.C
