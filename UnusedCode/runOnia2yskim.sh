#!/bin/bash
'''
int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;
int kPPMCUps1S = 8 ;
int kPPMCUps2S = 9 ;
int kPPMCUps3S = 10 ;
int kAAMCUps1S = 11 ;
int kAAMCUps2S = 12 ;
int kAAMCUps3S = 13 ;

int kNoTrigSel       =0;
int kL1DoubleMu0      =1;
int kL3JpsiCentral    =2;
int kL3UpsilonCentral =3;
int kL1DoubleMu0Peripheral=4;
int kL1DoubleMu10     =5;

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;

void onia2ySkim( int nevt = -1,
                 int fileID = kAAMCUps1S,
                 int trigId=kL1DoubleMu0,
                 int epSelection = kEPOppositeHF,
                 bool saveTracks=false,
                 TString skimVersion="unIdentified",
                 bool DiMuSign = false
                 ) {
'''

#### create outputDir if it does not exist
if [ ! -d "$(pwd)/skimmedFiles" ]; then
	mkdir $(pwd)/skimmedFiles
fi

#### options
DiMuSign='0' ## (1:same sign, 0:opposite sign)
nevt=-1 ## -1: all
gitVer=$(git show | head -1 | awk '{print $2}')

#### official PPMC
root -l -q -b 'onia2ySkim.C+('$nevt', 8, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 1S
#### official PPMC
root -l -q -b 'onia2ySkim.C+('$nevt', 9, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 2S
#### official PPMC
root -l -q -b 'onia2ySkim.C+('$nevt', 10, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 3S

#### official AAMC
root -l -q -b 'onia2ySkim.C+('$nevt', 11, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 1S
#### official AAMC
root -l -q -b 'onia2ySkim.C+('$nevt', 12, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 2S
#### official AAMC
root -l -q -b 'onia2ySkim.C+('$nevt', 13, 1, 1, 0, "'$gitVer'", '$DiMuSign')'   # pp mc 3S
