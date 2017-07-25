#ifndef TriggerManipulation_h
#define TriggerManipulation_h

#include "cutsAndBin.h"

int kNoTrigSel       =0;
int kL1DoubleMu0      =1;
int kL3JpsiCentral    =2 ;
int kL3UpsilonCentral =3;
int kL1DoubleMu0Peripheral=4;
int kL1DoubleMu10     =5;
int kL1DoubleMu0NoHFBit     =11;
int kL1DoubleMuOpen2016     =12;



struct hltIndex { int doTrigSel, ind1, ind2,  ind3,  ind4; };

TString getTrig( int trigId ) {
  if ( trigId == kNoTrigSel )  return "noTrigSel" ;
  else if ( trigId == kL1DoubleMu0 ) return "L1DoubleMu0" ;
  else if ( trigId == kL3JpsiCentral ) return "L3JpsiCentral" ;
  else if ( trigId == kL3UpsilonCentral ) return "L3UpsilonCentral" ;
  else if ( trigId == kL1DoubleMu0Peripheral ) return "L1DoubleMu0Peripheral" ;
  else if ( trigId == kL1DoubleMu10 ) return "L1DoubleMu10" ;
  else if ( trigId == kL1DoubleMu0NoHFBit ) return "L1DoubleMu0NoHFRequired" ;
  else if ( trigId == kL1DoubleMuOpen2016 ) return "L1DoubleMuOpen2016" ;
  else return "none" ;
}

hltIndex getTrigIndex( int trigId, TString fname ) { 
  hltIndex ret = {1, -1, -1, -1, -1};
  TString trigName = getTrig(trigId);
  cout << endl << endl << "*==*==*==*==*==*==* Trigger selection  *==*==*==*==*==*==*" << endl;
  cout << " trigger selection : " << trigName << endl;

  TFile* f1 = new TFile(fname.Data());
  TH1F* hStats = (TH1F*)f1->Get("hionia/hStats");
  if ( hStats == 0 ){
    cout << " it's not the oniaAndFriends tree made by TFileServie.. hionia/hStats -> hStats myTree instead of hionia/myTree" << endl;
    hStats = (TH1F*)f1->Get("hStats");
  }
  cout << endl << "*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*" << endl;

  if ( trigId == kNoTrigSel ) {
    ret.doTrigSel = 0 ;
  }
  else { 
    ret.doTrigSel = 1 ;
    int nBins = hStats->GetNbinsX() ;
    for ( int ii=1 ; ii<=nBins/2 ; ii++) { 
      TString theName = hStats->GetXaxis()->GetBinLabel(ii); 
      
      if ( trigId == kL1DoubleMu0 ) { 
	if ( theName == "HLT_HIL1DoubleMu0_v1" )  {               ret.ind1 = ii-2;  
	  cout << "HLT_HIL1DoubleMu0_v1 : ind1 = 2^" << ret.ind1 << endl;     }
	if ( theName == "HLT_HIL1DoubleMu0_2HF_v1" )   {          ret.ind2 = ii-2;
	  cout << "HLT_HIL1DoubleMu0_2HF_v1 : ind2 = 2^" << ret.ind2 << endl;     }
	if ( theName == "HLT_HIL1DoubleMu0_2HF0_v1" )   { 	ret.ind3 = ii-2;
	  cout << "HLT_HIL1DoubleMu0_2HF0_v1 : ind3 = 2^" << ret.ind3 << endl;     }
	if ( theName == "HLT_HIL1DoubleMu0ForPPRef_v1" )   { 	ret.ind4 = ii-2;
	  cout << "HLT_HIL1DoubleMu0ForPPRef_v1 (pp Reference) : ind4 = 2^" << ret.ind4 << endl;     }
      }

      if ( trigId == kL1DoubleMuOpen2016 ) { 
	if ( theName == "HLT_PAL1DoubleMuOpen_v1" )  {               ret.ind1 = ii-2;
          cout << "kL1DoubleMuOpen2016 : ind1 = 2^" << ret.ind1 << endl;     }
      }

      if ( trigId == kL1DoubleMu0NoHFBit ) { 
	if ( theName == "HLT_HIL1DoubleMu0_v1" )  {               ret.ind1 = ii-2;  
	  cout << "HLT_HIL1DoubleMu0_v1 : ind1 = 2^" << ret.ind1 << endl;     }
      }
      
      else if ( trigId == kL3JpsiCentral ) { 
	if ( theName == "HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1" )  { 	ret.ind1 = ii-2;
	  cout << "HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 : ind1 = 2^" << ret.ind1 << endl;     }
      }
      
      else if ( trigId == kL3UpsilonCentral ) { 
	if ( theName == "HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1" ) {   	ret.ind1 = ii-2;
	  cout << "HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1 : ind1 = 2^"    << ret.ind1 << endl;     }
      }
      else if ( trigId == kL1DoubleMu0Peripheral ) {
	if ( theName == "HLT_HIL1DoubleMu0_2HF_Cent30100_v1" )  {         ret.ind1 = ii-2;
	  cout << "HLT_HIL1DoubleMu0_2HF_Cent30100_v1 : ind1 = 2^"    << ret.ind1 << endl;     }
	if ( theName == "HLT_HIL1DoubleMu0_2HF0_Cent30100_v1" )  {        ret.ind2 = ii-2;
	  cout << "HLT_HIL1DoubleMu0_2HF0_Cent30100_v1 : ind1 = 2^"   << ret.ind2 << endl;     }
      }
      
    }
  }
  
  cout << "Do trigger filter? = " << ret.doTrigSel << endl;
  cout << "index 1-4 : " << ret.ind1 << ", "<< ret.ind2 << ", "<<ret.ind3<<", "<<ret.ind4<<endl;
  return ret;
}




bool isTrigMatched(hltIndex hltind, ULong64_t hltInput) {
  bool ret=false;

  if ( hltind.doTrigSel == 0 ) 
    ret = true ; 
  else { 
    ULong64_t ind1_ = (ULong64_t)(pow(2,hltind.ind1));
    ULong64_t ind2_ = (ULong64_t)(pow(2,hltind.ind2));
    ULong64_t ind3_ = (ULong64_t)(pow(2,hltind.ind3));
    ULong64_t ind4_ = (ULong64_t)(pow(2,hltind.ind4));
    
    if (  (ind1_>0) && ( (hltInput&ind1_)==ind1_ ) )     
      ret = true;
    if (  (ind2_>0) && ( (hltInput&ind2_)==ind2_ ) )     
      ret = true;
    if (  (ind3_>0) && ( (hltInput&ind3_)==ind3_ ) )     
      ret = true;
    if (  (ind4_>0) && ( (hltInput&ind4_)==ind4_ ) )     
      ret = true;
  }
  return ret;
}


#endif


