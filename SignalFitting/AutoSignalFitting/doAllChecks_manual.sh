#!/bin/bash
changeParam=1

if [ "$1" != "" ]
then
  changeParam=$1
  shift
fi

# PP DATA
<<COMMENT
#regular pt bins
./DoCheck.sh kPPDATA 0.0 2.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 2.0 4.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 4.0 6.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 6.0 9.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 9.0 12.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 12.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPPDATA 0.0 4.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 4.0 9.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 9.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPPDATA 0.0 6.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 6.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#regular y bins
./DoCheck.sh kPPDATA 0.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 0.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 0.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 0.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPPDATA 0.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPPDATA 0.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPPDATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#y bins in low and high pt for 1S and 2S
#./DoCheck.sh kPPDATA 0.0 6.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 0.0 6.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 0.0 6.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 0.0 6.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPPDATA 0.0 6.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 0.0 6.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPPDATA 6.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 6.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 6.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 6.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPPDATA 6.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPPDATA 6.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam
COMMENT
# PA DATA

#integrated bin
./DoCheck.sh kPADATA 0.0 30.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#regular pt bins
./DoCheck.sh kPADATA 0.0 2.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 2.0 4.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 4.0 6.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 9.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 9.0 12.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 12.0 30.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 4.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 4.0 9.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 9.0 30.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 6.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 30.0 -1.93 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#regular y bins
./DoCheck.sh kPADATA 0.0 30.0 -1.93 -1.2 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.2 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -0.8 -0.4 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -0.4 0.0 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 30.0 -1.93 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -0.8 0.0 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#pt bins in backward and forward rapitidy
#./DoCheck.sh kPADATA 0.0 2.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 2.0 4.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 4.0 6.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 9.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 9.0 12.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 12.0 30.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 0.0 4.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 4.0 9.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 9.0 30.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 6.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 30.0 -1.93 0.0 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 0.0 2.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 2.0 4.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 4.0 6.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 9.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 9.0 12.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 12.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 0.0 4.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 4.0 9.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 9.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 6.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 30.0 0.0 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#y bins in low and high pt
#./DoCheck.sh kPADATA 0.0 6.0 -1.93 -1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 -1.2 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 -0.8 -0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 -0.4 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 0.0 6.0 -1.93 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 -0.8 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 0.0 6.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 6.0 30.0 -1.93 -1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 -1.2 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 -0.8 -0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 -0.4 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#./DoCheck.sh kPADATA 6.0 30.0 -1.93 -0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 -0.8 0.0 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $changeParam
#./DoCheck.sh kPADATA 6.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#integrated bin in -2.87<y<1.93 in pPb.
./DoCheck.sh kPADATA 0.0 30.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#pt bins in rapidity range [-2.87,1.93]
./DoCheck.sh kPADATA 0.0 2.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 2.0 4.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 4.0 6.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 9.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 9.0 12.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 12.0 30.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 4.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 4.0 9.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 9.0 30.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 6.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 6.0 30.0 -2.87 1.93 0 200 4.0 0 120 0 400 0 $changeParam

#-2.87<y<-1.93 in pPb.
./DoCheck.sh kPADATA 0.0 30.0 -2.87 -1.93 0 200 4.0 0 120 0 400 0 $changeParam

#HF bins in backward and forward y
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 12 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 12 19 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 19 27 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 27 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 12 120 0 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 12 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 12 19 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 19 27 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 27 120 0 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 12 120 0 400 0 $changeParam

#HF bins in backward and forward y
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 0 40 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 40 62 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 62 88 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 88 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 -1.93 0.0 0 200 4.0 0 120 40 400 0 $changeParam

./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 0 40 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 40 62 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 62 88 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 88 400 0 $changeParam
./DoCheck.sh kPADATA 0.0 30.0 0.0 1.93 0 200 4.0 0 120 40 400 0 $changeParam





