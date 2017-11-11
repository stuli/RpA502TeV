#!/bin/bash
#$1 = OniaMode (1,2,3), $2 = ispPb (1,0), $3 = pPb or pp, $4 = Name (eg. 29_5)
mkdir eff_$3$4
cp /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_RpA_ana.C /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/VVV/$1/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/WWW/$2/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/XXX/$3/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/RpA_ana/oniaMode$1_$3_$4/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
sed -i "s/TAG/$4/g;" /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
cp /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/eff_$3$4/dimuEff_oniaMode$1_$3_$4.C
#cp /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/{effCommon.h, output_official_5eta_cutG_all_nominal_v3.root, tdrstyle.C, CMS_lumi.C, CMS_lumi.h, tnp_weight.h} /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/eff_$3$4/
root -l -q -b /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/dimuEff_oniaMode$1_$3_$4.C
#This line is optional ^^^ You can run later
mv /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/*.png /home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/CrossChecks/efficiency_Santona/New/eff_$3$4/


