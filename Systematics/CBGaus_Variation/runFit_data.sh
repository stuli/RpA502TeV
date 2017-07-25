# int collId = kPPDATA,                   float ptLow=5, float ptHigh=100,
#void doFitUpsilon_Data(
#       int collId = kPPDATA,
#       float ptLow=0, float ptHigh=1,
#       float yLow=0, float yHigh=1.2,
#       int cLow=0, int cHigh=200,
#       float muPtCut=4.0,
#       bool fixParameters=1
#                        )


muPtCut='4'
fixParam='1'


kPPDATA='0'
kAADATA='2'
kPPMCUps1S='8'
kPPMCUps2S='9'
kPPMCUps3S='10'

#outputDir='fitResults/dataFit_fixParam'${fixParam}'MuPt'${muPtCut}'_2016_08_30'
outputDir='MassPlot'

#outputDir='Test/fixParam'${fixParam}'MuPt'${muPtCut}'_fullRap_mar01_bkgtest'${fbkg_ch}''
mkdir -p $outputDir
cp *.h tdrstyle.C CMS_lumi.C doFitUpsilon_Data_CBGaus.C  $outputDir

root -l -b <<EOF
.L doFitUpsilon_Data_CBGaus.C++
.q
EOF

echo entering to the loop...

for i in ${kPPDATA}  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
    for pt in '0,5'
    #for pt in '0,2.5' '2.5,5' '5,8' '8,15' '15,30' '0,5' '5,15'
    do
	 #   for rap in '0,2.4' '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4' '0,1.2' '1.2,2.4'
	    for rap in '0,2.4'
        do
            echo  $outputDir 'doFitUpsilon_Data_CBGaus.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
            #./condor_root.sh $outputDir 'doFitUpsilon_Data.C('$i','$pt','$rap',0,200,'$muPtCut',0)'
           # root -l -q -b $outputDir 'doFitUpsilon_Data_CBGaus.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
        done
    done
done

for i in ${kAADATA}  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
   #  for cent in '20,40'
    for cent in '0,10' '10,20' '20,40' '40,60' '0,20' '20,60'
      do
	      for rap in '0,2.4'
        do
	      echo  $outputDir 'doFitUpsilon_Data_CBGaus.C('$i',0,30,'$rap','$cent','$muPtCut','$fixParam')'
#	      ./condor_root.sh $outputDir 'doFitUpsilon_Data.C('$i',0,30,'$rap','$cent','$muPtCut','$fixParam')'
	      root -l -q -b $outputDir 'doFitUpsilon_Data_CBGaus.C('$i',0,30,'$rap','$cent','$muPtCut','$fixParam')'
	      done
      done
    for cent1 in '60,80' '80,100' '100,120' '120,140' '140,200' '60,100' '100,200'
 #   for cent1 in '100,120'
    do
      for rap in '0,2.4'
      do
	    echo  $outputDir 'doFitUpsilon_Data_CBGaus.C(6,0,30,'$rap','$cent1','$muPtCut','$fixParam')'
#	    ./condor_root.sh $outputDir 'doFitUpsilon_Data.C(6,0,30,'$rap','$cent1','$muPtCut','$fixParam')'
	    root -l -q -b $outputDir 'doFitUpsilon_Data_CBGaus.C(6,0,30,'$rap','$cent1','$muPtCut','$fixParam')'
	    done
    done

done




