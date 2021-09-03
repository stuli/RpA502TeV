#!/bin/bash
collId=kPADATA
ptLow=0.0
ptHigh=30.0
yLow=0.0
yHigh=2.4
cLow=0
cHigh=200
muPtCut=4.0
hfLow=0
hfHigh=120
ntracksLow=0
ntracksHigh=400
whichModel=0
changeParam=0
nTries=1

if [ "$1" != "" ]
then
  collId=$1
  shift
fi
if [ "$1" != "" ]
then
  ptLow=$1
  shift
fi
if [ "$1" != "" ]
then
  ptHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  yLow=$1
  shift
fi
if [ "$1" != "" ]
then
  yHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  cLow=$1
  shift
fi
if [ "$1" != "" ]
then
  cHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  muPtCut=$1
  shift
fi
if [ "$1" != "" ]
then
  hfLow=$1
  shift
fi
if [ "$1" != "" ]
then
  hfHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  ntracksLow=$1
  shift
fi
if [ "$1" != "" ]
then
  ntracksHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  whichModel=$1
  shift
fi
if [ "$1" != "" ]
then
  changeParam=$1
  shift
fi
if [ "$1" != "" ]
then
  nTries=$1
  shift
fi

if [ $nTries -lt 2 ]
then
  echo root -b -q -l "CheckFitExists.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,$changeParam)"
  root -b -q -l "CheckFitExists.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,$changeParam)" >> junk$changeParam
  isIt=$(tail -c 2 junk$changeParam)
  rm junk$changeParam
else
  isIt=1
fi

if [ $isIt -gt 0 ]
then
  if [ $nTries -lt 2 ]
  then
    echo "THE FIT EXISTS ALREADY! :)"
  fi
  echo root -b -q -l "CheckFitQuality.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,$changeParam)"
  root -b -q -l "CheckFitQuality.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,$changeParam)" >> junk$changeParam
  ToF=$(tail -c 2 junk$changeParam)
  rm junk$changeParam
else
  echo "THE FIT DOES NOT EXIST YET! :O"
  echo
  ToF=0
fi

if [ $ToF -gt 0 ]
then
  echo "THE FIT PASSED THE QUALITY CHECK! :)"
  echo
else
  echo "THE FIT FAILED THE QUALITY CHECK! :("
  echo
  if [ $nTries -lt 21 ]
  then
    echo "STARTING TRY #$nTries"
    if [ $nTries -lt 2 ]
    then
      echo root -b -q -l "FitDataLikeAPro.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,kFALSE,$changeParam)"
      root -b -q -l "FitDataLikeAPro.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,kFALSE,$changeParam)"
    else
      echo root -b -q -l "FitDataLikeAPro.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,kTRUE,$changeParam)"
      root -b -q -l "FitDataLikeAPro.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh,$whichModel,kTRUE,$changeParam)"
    fi
    nTries=$(($nTries+1))
    ./DoFit.sh $collId $ptLow $ptHigh $yLow $yHigh $cLow $cHigh $muPtCut $hfLow $hfHigh $ntracksLow $ntracksHigh $whichModel $changeParam $nTries
  else
    nTries=$(($nTries-1))
    echo "GIVING UP AFTER $nTries TRIES :("
    root -b -q -l "WriteToLog.C()"
    echo
  fi
fi

