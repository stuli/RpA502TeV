#!/bin/bash
collId=kPADATA
ptLow=0.0
ptHigh=30.0
yLow=0.0
yHigh=2.4
hfLow=0
hfHigh=120
ntracksLow=0
ntracksHigh=400

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
cd /afs/cern.ch/work/j/jjay/public/RpA502TeV/PseudoExperiments
root -b -q "FitPseudoData_Constrained.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$hfLow,$hfHigh,$ntracksLow,$ntracksHigh)";
