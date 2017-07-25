# int collId = kPPDATA,                   float ptLow=5, float ptHigh=100,
# float yLow=0, float yHigh=1.2,                   int cLow=0, int cHigh=20,
#void doFitUpsilon_MC(
#   int collId = kAAMCUps1S,
#       float ptLow=0, float ptHigh=6,
#       float yLow=0, float yHigh=1.2,
#       int cLow=0, int cHigh=50,
#       float muPtCut=4.0



muPtCut='4'

kPPMCUps1S='8'
kPPMCUps2S='9'
kPPMCUps3S='10'
kAAMCUps1S='11'
kAAMCUps2S='12'
kAAMCUps3S='13'

outputDir='fitResults/mcFit_MuPt'${muPtCut}'_2016_11_04'

#outputDir='Test/fixParam'${fixParam}'MuPt'${muPtCut}'_fullRap_mar01_bkgtest'${fbkg_ch}''
mkdir -p $outputDir
cp onia2ySkim.C doFitUpsilon_MC.C PsetCollection.h rootFitHeaders.h commonUtility.h cutsAndBin.h $outputDir

root -l -b <<EOF
.L doFitUpsilon_MC.C++
.q
EOF

echo entering to the loop...

for i in ${kPPMCUps1S}  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
#    for pt in '0,2' '2,4' '4,6' '6,8' '8,12' '12,16' '16,30'  '0,8' '8,30' '0,30'
#    for pt in '0,2.5' '0,5' '2.5,5' '5,8' '5,15' '8,15' '15,30'
    for pt in '0,30'
    do
#	for rap in '0,1.2' '1.2,2.4' '0,2.4'
  #for rap in '0,2.4' '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4' '0,1.2' '1.2,2.4'
  for rap in '0,2.4'
	do
	    echo  $outputDir 'doFitUpsilon_MC.C('$i','$pt','$rap',0,200,'$muPtCut')'
	    root -l -b -q $outputDir 'doFitUpsilon_MC.C++('$i','$pt','$rap',0,200,'$muPtCut')'
	done
    done
done

for i in ${kAAMCUps1S}  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
    #for cent in '0,10' '10,20' '20,40' '40,60' '60,80' '80,100' '100,140' '140,200' '0,20' '20,60' '60,100' '100,200'
    for cent in '100,120' '120,140' '140,200'
    do
	for rap in '0,1.2' '1.2,2.4' '0,2.4'
	do
	    echo  $outputDir 'doFitUpsilon_MC.C('$i',0,30,'$rap','$cent','$muPtCut')'
#	    ./condor_root.sh $outputDir 'doFitUpsilon_MC.C('$i',0,30,'$rap','$cent','$muPtCut')'
	done
    done
done



