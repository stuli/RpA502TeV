#      int collId = kPPDATA,
#           float4ptLow=0, float ptHigh=5,
#           float yLow=0, float yHigh=2.4,
#           int cLow=0, int cHigh=160,
#           float muPtCut=4.0,
#           int inputOption=kErrExp, //kChPol3,
#           int nGen = 100000
#           int useCentIntBkgShape = 1
root -l -b <<EOF
.L toyMC.C++
.q
EOF

muPtCut='4'
nGen='1000000'
useCentIntBkgShape=1
nToys='1'
for (( irnd=501 ; irnd<=547 ; irnd++))
do
  for bkgOption in 3 # 1=Pol3, 2=ErrExp, #4=ErrExpExp #3 = pol4
  do
    #outputDir='ffsys_BkgVar/TEST'
    #outputDir='ffsys_BkgVar/finalbkg_Gen100K'
    outputDir='TOY_Pre'
    #outputDir='TOYPLOT_TEST_SB'
    #outputDir='ffsys_BkgVar/finalbkg_GenNom'
    #outputDir='ffsys_BkgVar/Aug3_toyMC_MuPt'${muPtCut}'_nGen'${nGen}
    mkdir -p $outputDir
    cp toyMC.C *.h $outputDir
    
    echo entering to the loop...
    #Centrality dependence 
    for i in  2   # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
    do
	    for pt in '0,30'
	    do
	      for rap in '0,2.4'
	      do
    		for cent in '0,200'
    		#for cent in '0,10' '10,20' '20,40' '40,60' '60,80' '80,100' '100,120' '120,140' '140,200'
          do
 #         echo ..
  	      root -l -b -q  $outputDir 'toyMC.C('$i','$pt','$rap','$cent','$muPtCut','$bkgOption','$nGen','$useCentIntBkgShape','$nToys','$irnd')'
		      done
	      done
	    done
    done
    
    for i in 0 2  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3 8=kPPMCUps1S 9=kPPMCUps2S
    do
	    for cent in '0,200'
        do
	        for pt in '0,30'
	        do
		        for rap in '0,2.4' '0.0,1.2' '1.2,2.4' '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4'
#		        for rap in '0.0,1.2' 
		        do
              echo ..
#		          root -l -q -b $outputDir 'toyMC.C('$i','$pt','$rap','$cent','$muPtCut','$bkgOption','$nGen','$useCentIntBkgShape','$nToys','$irnd')'
		        done
	        done
          for rap in '0,2.4'
          do
		        for pt in '0,5' '5,15' '15,30' '0,2.5' '2.5,5' '5,8' '8,15'
#		        for pt in '0,30' 
		        do
 #   	        root -l -q -b $outputDir 'toyMC.C('$i','$pt','$rap','$cent','$muPtCut','$bkgOption','$nGen','$useCentIntBkgShape','$nToys','$irnd')'
		        echo ..
            done
          done
	      done
    done
 
  done
done

