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
outputDir='4thorder'

#outputDir='Test/fixParam'${fixParam}'MuPt'${muPtCut}'_fullRap_mar01_bkgtest'${fbkg_ch}''
mkdir -p $outputDir
cp doSys_BkgVar.C PsetCollection.h cutsAndBin.h $outputDir

root -l -b <<EOF
.L doSys_BkgVar.C++
.q
EOF

echo entering to the loop...

for i in ${kPPDATA} ${kAADATA} # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
    for pt in '0,30'
    do
	    for rap in '0.8,1.2'
	    #for rap in '0,0.8' '0.8,1.6' '1.6,2.4' '0,1.2' '1.2,2.4' '0,2.4'
        do
          echo  $outputDir 'doSys_BkgVar.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
            ./condor_root.sh $outputDir 'doSys_BkgVar.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
        done
    done
    for pt in '0,4' '4,9' '9,30'
    #for pt in '0,2' '2,4' '4,6' '6,9' '9,12' '12,30' '0,4' '4,9' '9,30' '0,6' '6,30'
    do
	    for rap in '0,2.4'
        do
          echo  $outputDir 'doSys_BkgVar.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
        #    ./condor_root.sh $outputDir 'doSys_BkgVar.C('$i','$pt','$rap',0,200,'$muPtCut','$fixParam')'
        done
    done
done

for i in ${kAADATA}  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
do
    for cent in '0,10' '10,20' '20,40' '40,60' '0,60'
    do
	      for rap in '0,2.4'
        do
          echo  $outputDir 'doSys_BkgVar.C('$i',0,30,'$rap','$cent','$muPtCut','$fixParam')'
         # ./condor_root.sh $outputDir 'doSys_BkgVar.C('$i',0,30,'$rap','$cent','$muPtCut','$fixParam')'
	      done
      done
    for cent1 in '60,80' '80,100' '100,120' '120,140' '140,200' '60,200'
    do
      for rap in '0,2.4'
      do
        echo  $outputDir 'doSys_BkgVar.C(6,0,30,'$rap','$cent1','$muPtCut','$fixParam')'
	    #./condor_root.sh $outputDir 'doSys_BkgVar.C(6,0,30,'$rap','$cent1','$muPtCut','$fixParam')'
	    done
    done

done




