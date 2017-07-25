#ifndef PsetCollection_h
#define PsetCollection_h

#include "cutsAndBin.h"


class  PSetUpsAndBkg { 
 public:
 PSetUpsAndBkg() :
  collId(0), muPtCut(0), ptLow(0), ptHigh(0), yLow(0), yHigh(0), cLow(0), cHigh(200),
    n1s_1(-1), n2s_1(0), n3s_1(0),
    n1s_2(0), n2s_2(0), n3s_2(0),
    alpha1s_1(0), alpha2s_1(0), alpha3s_1(0),
    alpha1s_2(0), alpha2s_2(0), alpha3s_2(0),
    sigma1s_1(0), sigma2s_1(0), sigma3s_1(0),
    sigma1s_2(0), sigma2s_2(0), sigma3s_2(0),
    mean1s(0), mean2s(0), mean3s(0),
    f1s(0), f2s(0), f3s(0),
    x1s(0),
    //    MCn(0), MCalpha(0), MCsigma1S(0), MCsigma2S(0), MCsigma3S(0), MCm0(0), MCf(0), MCx(0),
    
    bkg_mu(0), bkg_sigma(0), bkg_lambda(0),
    // Only For Systematics
    bkg_mu1(0), bkg_sigma1(0), bkg_lambda1(0), bkg_mu2(0), bkg_sigma2(0), bkg_lambda2(0), rBkg2over1(0), // double ErrFunction
    ch3_k1(0), ch3_k2(0), ch3_k3(0),
    ch4_k1(0), ch4_k2(0), ch4_k3(0), ch4_k4(0),
    nSignal1s(0), nSignal2s(0), nSignal3s(0), nBkg(0),
    bkg4_mu(0), bkg4_sigma(0), bkg4_lambda(0),  bkg4_lambda2(0), rBkg42over1(0) //bkg4 = exp + err*exp
    
    {}
  
  int collId;
  float muPtCut;
  float ptLow;
  float ptHigh;
  float yLow;
  float yHigh;
  int cLow;
  int cHigh;

  float mean1s;
  float mean2s;  float mean3s;

  float n1s_1;
  float n1s_2;
  float alpha1s_1;
  float alpha1s_2;
  float sigma1s_1;
  float sigma1s_2;
  float n2s_1;
  float n2s_2;
  float alpha2s_1;
  float alpha2s_2;
  float sigma2s_1;
  float sigma2s_2;
  float n3s_1;
  float n3s_2;
  float alpha3s_1;
  float alpha3s_2;
  float sigma3s_1;
  float sigma3s_2;
  float x1s;

  float f1s; 
  float f2s; 
  float f3s; 


  //  float MCn, MCalpha, MCsigma1S, MCsigma2S, MCsigma3S, MCm0, MCf, MCx;
  float bkg_mu, bkg_sigma, bkg_lambda;

  float bkg_mu1, bkg_sigma1, bkg_lambda1, bkg_mu2, bkg_sigma2, bkg_lambda2, rBkg2over1; // double ErrFunction
  float ch3_k1, ch3_k2, ch3_k3 ; 
  float ch4_k1, ch4_k2, ch4_k3, ch4_k4 ; 

  float nSignal1s, nSignal2s, nSignal3s, nBkg;

  float bkg4_mu, bkg4_sigma, bkg4_lambda, bkg4_lambda2, rBkg42over1;

  void setKine(int collId_, float muPtCut_, float ptLow_, float ptHigh_, float yLow_, float yHigh_, int cLow_, int cHigh_)
  {
    collId = collId_;    muPtCut = muPtCut_;
    ptLow  = ptLow_;     ptHigh  = ptHigh_;
    yLow   = yLow_;      yHigh   = yHigh_;
    cLow   = cLow_;      cHigh   = cHigh_;
    cout << " collId : " << getCollID( collId) << endl;
    cout << " pT     : " << ptLow << " - " << ptHigh << "GeV/c" << endl;
    cout << " y      : " << yLow << " - " << yHigh << "" << endl;
    cout << " cBin   : " << cLow << " - " << cHigh << "" << endl;
    cout << " Muon pT > " << muPtCut << endl;
  }

  void setParBkg(float bkg_mu_, float bkg_sigma_, float bkg_lambda_)
  {
    bkg_mu = bkg_mu_; bkg_sigma = bkg_sigma_; bkg_lambda = bkg_lambda_;
  }
  
  void setSignalParMC(float MCn_, float MCalpha_, float MCsigma1S_, float MCm0_, float MCf_, float MCx_)
  {
    n1s_1 = MCn_ ;    n2s_1 = MCn_ ;    n3s_1 = MCn_ ;
    n1s_2 = MCn_ ;    n2s_2 = MCn_ ;    n3s_2 = MCn_ ;
    alpha1s_1 = MCalpha_ ;     alpha2s_1 = MCalpha_ ;     alpha3s_1 = MCalpha_ ; 
    alpha1s_2 = MCalpha_ ;     alpha2s_2 = MCalpha_ ;     alpha3s_2 = MCalpha_ ; 
    sigma1s_1 = MCsigma1S_ ;         sigma2s_1 = MCsigma1S_ * (pdgMass.Y2S/pdgMass.Y1S) ;     sigma3s_1 = MCsigma1S_ * (pdgMass.Y3S/pdgMass.Y1S) ; 
    sigma1s_2 = sigma1s_1*MCx_ ;     sigma2s_2 = sigma2s_1*MCx_ ;  sigma3s_2 = sigma3s_1*MCx_  ;
    mean1s  = MCm0_ ;   mean2s  = MCm0_ * (pdgMass.Y2S/pdgMass.Y1S);  mean3s  = MCm0_ * (pdgMass.Y3S/pdgMass.Y1S); 
    f1s = MCf_;     f2s= MCf_;    f3s= MCf_;   x1s = MCx_;
  }

  void setSignalParPPDATA(float MCn_, float MCalpha_, float MCsigma1S_, float MCm0_, float MCf_, float MCx_)
  {
    n1s_1 = MCn_ ;
    alpha1s_1 = MCalpha_ ;  sigma1s_1 = MCsigma1S_ ;  sigma2s_1 = MCsigma1S_ * (pdgMass.Y2S/pdgMass.Y1S) ;  sigma3s_1 = MCsigma1S_ * (pdgMass.Y3S/pdgMass.Y1S) ; 
    sigma1s_2 = sigma1s_1*MCx_ ;     sigma2s_2 = sigma2s_1*MCx_ ;      sigma3s_2 = sigma3s_1*MCx_  ;
    mean1s  = MCm0_ ;   mean2s  = MCm0_ * (pdgMass.Y2S/pdgMass.Y1S);  mean3s  = MCm0_ * (pdgMass.Y3S/pdgMass.Y1S); 
    f1s = MCf_;     f2s= MCf_;    f3s= MCf_;
    x1s = MCx_;
  }

  void setParBkg2ErrExp(float bkg_mu1_, float bkg_sigma1_, float bkg_lambda1_, float bkg_mu2_, float bkg_sigma2_, float bkg_lambda2_, float rBkg2over1_)
  {
    bkg_mu1 = bkg_mu1_;  bkg_sigma1 = bkg_sigma1_;  bkg_lambda1 = bkg_lambda1_;
    bkg_mu2 = bkg_mu2_;  bkg_sigma2 = bkg_sigma2_;  bkg_lambda2 = bkg_lambda2_; rBkg2over1 = rBkg2over1_;
  }

  void setParBkgErrExpExp(float bkg4_mu_, float bkg4_sigma_, float bkg4_lambda_, float bkg4_lambda2_, float rBkg42over1_)
  {
    bkg4_mu = bkg4_mu_;  bkg4_sigma = bkg4_sigma_;  bkg4_lambda = bkg4_lambda_;
    bkg4_lambda2 = bkg4_lambda2_; rBkg42over1 = rBkg42over1_;
  }
  
  void setParBkgPol3(float k1_, float k2_, float k3_) 
  {
    ch3_k1 = k1_ ;      ch3_k2 = k2_ ;       ch3_k3 = k3_ ;
  }

  void setParBkgPol4(float k1_, float k2_, float k3_, float k4_) 
  {
    ch4_k1 = k1_ ;      ch4_k2 = k2_ ;       ch4_k3 = k3_ ;   ch4_k4 = k4_;
  }

  void setSig1sF21NBkg(float sig1s_, float f21_, float nBkg_)
  {
    nSignal1s = sig1s_;
    nSignal2s = sig1s_ * f21_;
    nBkg = nBkg_;
  }
  
  
  
  void reset() {
    n1s_1 = -1 ;  n1s_2 = 0 ;  n2s_1 = 0 ;  n2s_2 = 0 ;  n3s_1 = 0 ;  n3s_2 =0 ; 
    alpha1s_1 = 0 ;  alpha1s_2 = 0 ;  alpha2s_1 = 0 ;  alpha2s_2 = 0 ;  alpha3s_1 = 0 ;  alpha3s_2 = 0 ;
    sigma1s_1 = 0 ;  sigma1s_2 = 0 ;  sigma2s_1 = 0 ;  sigma2s_2 = 0 ;  sigma3s_1 = 0 ;  sigma3s_2 = 0 ;
    //    MCn = 0 ;  MCalpha = 0 ;  MCsigma1S = 0 ;  MCsigma2S = 0 ;  MCsigma3S = 0 ;  MCm0 = 0 ;  MCf = 0 ;  MCx = 0 ;
    bkg_mu = 0; bkg_sigma = 0; bkg_lambda = 0;
    bkg_mu1 = 0; bkg_sigma1 = 0; bkg_lambda1 = 0; bkg_mu2 = 0; bkg_sigma2 = 0; bkg_lambda2=0; rBkg2over1=0; // double ErrFunction
    ch3_k1 = 0; ch3_k2 = 0; ch3_k3=0 ;
    ch4_k1 = 0; ch4_k2 = 0; ch4_k3=0 ; ch4_k4=0;
    nSignal1s=0; nSignal2s=0; nSignal3s=0; nBkg=0;
    bkg4_mu = 0; bkg4_sigma = 0; bkg4_lambda = 0; bkg4_lambda2 = 0; rBkg42over1=0;
    
  }

  void SetMCSgl();
  void SetParPPDATASgl();
  void SetMCBkg();
  bool binMatched( float muonPtCut_, float ptLow_, float ptHigh_, float yLow_, float yHigh_, int cLow_=-1, int cHigh_=-1); 
  
};


PSetUpsAndBkg getUpsilonPsets( int collId = kPPDATA,
                               float ptLow=5, float ptHigh=100,
                               float yLow=1.2, float yHigh=2.4,
                               int cLow=0, int cHigh=200,
                               float muPtCut=4
                               )
{ 
  PSetUpsAndBkg ret;
  ret.setKine( collId, muPtCut, ptLow, ptHigh, yLow, yHigh, cLow, cHigh) ;
  return ret;
}


void PSetUpsAndBkg::SetMCBkg() {
  if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   {
    //Integrated Bin // No centrlaity dependence 
    if ( (muPtCut==4) && (ptLow == 0 ) && (ptHigh == 30 ) && (yLow == 0 ) && (yHigh == 2.4 ) )      
      {      setParBkg(8.52,1.28,7.95);    }

    else 
      {      setParBkg(7.86,1.02,6.08);}

  }
  else {  // pp
    if ( (muPtCut==4) && (ptLow == 0 ) && (ptHigh == 30 ) && (yLow == 0 ) && (yHigh == 2.4 ) )      
      {      setParBkg(7.86,1.02,6.08);}
    else 
      {      setParBkg(7.86,1.02,6.08);}
  }


}

void PSetUpsAndBkg::SetParPPDATASgl() {

    cout << " ///////////////////////////////////////////" << endl;
    cout << " Fixing the Parameters from PP Data..." << endl;
    cout << " ///////////////////////////////////////////" << endl;
    cout << " muPtCut = " << muPtCut << endl;
    cout << " ptLow = " << ptLow << endl;
    cout << " ptHigh = " << ptHigh << endl;
    cout << " yLow = " << yLow << endl;
    cout << " yHigh = " << yHigh << endl;

    if ( collId == kAADATA) {

if ( binMatched( 4, 0, 2, 0, 1.2) )   { setSignalParPPDATA( 1.62853, 1.88891, 0.113022, 9.4588, 0.276787, 0.534744 );} 
    }
}


void PSetUpsAndBkg::SetMCSgl() 
{
  
  /*    cout << " ///////////////////////////////////////////" << endl;
    cout << " MC signal PDFs are not defined for this bin" << endl; 
    cout << " ///////////////////////////////////////////" << endl;
  */
    cout << " ///////////////////////////////////////////" << endl;
    cout << " Fixing the Parameters..." << endl;
    cout << " ///////////////////////////////////////////" << endl;
    cout << " muPtCut = " << muPtCut << endl;
    cout << " ptLow = " << ptLow << endl;
    cout << " ptHigh = " << ptHigh << endl;
    cout << " yLow = " << yLow << endl;
    cout << " yHigh = " << yHigh << endl;
    
// mcFit_MuPt4_2016_11_30
if ( collId == kPPDATA) 
{
  //Bin for 1S 
  if ( binMatched( 4, 0, 2.5, 0, 2.4) )   { setSignalParMC( 3.4916, 1.564, 0.0656761, 9.45756, 0.54466, 1.97011);} 
  if ( binMatched( 4, 2.5, 5, 0, 2.4) )   { setSignalParMC( 3.35845, 1.63789, 0.0716633, 9.45802, 0.63576, 1.9029);} 
  if ( binMatched( 4, 5, 8, 0, 2.4) )   { setSignalParMC(  3.19046  ,  1.43575  ,  0.0730977  ,  9.45723  ,  0.639965  ,  1.78201);} 
  if ( binMatched( 4, 8, 15, 0, 2.4) )   { setSignalParMC(3.29824  ,  1.60905  ,  0.0678905  ,  9.45453  ,  0.537631  ,  1.87837 );} 
  if ( binMatched( 4, 15, 30, 0, 2.4) )   { setSignalParMC(3.25139  ,  1.60934  ,  0.0676854  ,  9.45591  ,  0.520767  ,  1.86066 );} 
  if ( binMatched( 4, 0, 30, 0, 0.4) )   { setSignalParMC(1.36109,  1.90884  ,  0.0889045  ,  9.4581  ,  0.189231  ,  0.608914 );} 
  if ( binMatched( 4, 0, 30, 0.4, 0.8) )   { setSignalParMC(2.30895  ,  1.75606  ,  0.0942427  ,  9.45866  ,  0.355442  ,  0.664363 );} 
  if ( binMatched( 4, 0, 30, 0.8, 1.2) )   { setSignalParMC(1.62317  ,  1.98761  ,  0.0777123  ,  9.45557  ,  0.497941  ,  1.45039 );} 
  if ( binMatched( 4, 0, 30, 1.2, 1.6) )   { setSignalParMC(1.38009  ,  2.12727  ,  0.0781824  ,  9.45002  ,  0.239506  ,  1.61113 );} 
  if ( binMatched( 4, 0, 30, 1.6, 2.0) )   { setSignalParMC(1.77962  ,  1.98774  ,  0.0726482  ,  9.45001  ,  0.0954821  ,  1.92647 );} 
  if ( binMatched( 4, 0, 30, 2.0, 2.4) )   { setSignalParMC(1.77998  ,  2.03591  ,  0.0864709  ,  9.45  ,  0.0852784  ,  1.99843 );} 
  if ( binMatched( 4, 0, 30, 0, 2.4) )   { setSignalParMC(3.30732  ,  1.60069  ,  0.0667393  ,  9.45626  ,  0.535528  ,  1.92417 );} 
  //Bin for 2S & 3S 
  if ( binMatched( 4, 0, 5, 0, 2.4) )   { setSignalParMC( 3.76, 1.59235, 0.0648477, 9.45739, 0.521473, 1.97614);} 
  if ( binMatched( 4, 5, 15, 0, 2.4) )   { setSignalParMC( 3.78201, 1.54632, 0.128471, 9.45551, 0.446142, 0.531346);} 
  if ( binMatched( 4, 0, 30, 0, 1.2) )   { setSignalParMC( 2.99101, 1.65623, 0.0981835, 9.45822, 0.419242, 0.605439);} 
  if ( binMatched( 4, 0, 30, 1.2, 2.4) )   { setSignalParMC(2.54509, 1.89829, 0.107104, 9.44942, 0.609942, 1.55899);} 
}

else if (collId == kAADATA ) 
{ 
 
//for pt bin
if ( binMatched( 4, 0, 2.5, 0, 2.4) )   { setSignalParMC( 3.45481, 1.59633, 0.0736743, 9.45525, 0.608694, 1.73631 );}
if ( binMatched( 4, 2.5, 5, 0, 2.4) )   { setSignalParMC( 3.40773, 1.62957, 0.0736023, 9.45613, 0.648529, 1.76715 );}
if ( binMatched( 4, 5, 8, 0, 2.4) )     { setSignalParMC( 3.72, 1.69586, 0.0676078, 9.45589, 0.553222, 1.90002 );}
if ( binMatched( 4, 8, 15, 0, 2.4) )    { setSignalParMC( 3.23, 1.51628, 0.0686302, 9.45512, 0.56218, 1.8871 );}
if ( binMatched( 4, 15, 30, 0, 2.4) )   { setSignalParMC( 3.3, 1.56701, 0.0734659, 9.45512, 0.612436, 1.85341 );}
if ( binMatched( 4, 0, 30, 0, 2.4) )    { setSignalParMC( 3.19451, 1.55708, 0.0684966, 9.45589, 0.560898, 1.90076 );}
if ( binMatched( 4, 0, 5, 0, 2.4) )     { setSignalParMC( 3.714, 1.70183, 0.0736797, 9.45614, 0.55781, 1.61859 );}
if ( binMatched( 4, 5, 15, 0, 2.4) )    { setSignalParMC( 3.72, 1.52551, 0.0683189, 9.45571, 0.563272, 1.89665 );}
//for eta bin
if ( binMatched( 4, 0, 30, 0, 0.4) )    { setSignalParMC( 1.53865, 1.88391, 0.0702923, 9.45862, 0.545296, 0.685079 );}
if ( binMatched( 4, 0, 30, 0.4, 0.8) )  { setSignalParMC( 1.92388, 1.84344, 0.0951642, 9.45858, 0.328146, 0.6797 );}
if ( binMatched( 4, 0, 30, 0.8, 1.2) )  { setSignalParMC( 3.49665, 1.75969, 0.071481, 9.45411, 0.335342, 1.49306 );}
if ( binMatched( 4, 0, 30, 1.2, 1.6) )  { setSignalParMC( 3.54659, 1.79847, 0.107393, 9.45, 0.87389, 1.59254 );}
if ( binMatched( 4, 0, 30, 1.6, 2.0) )  { setSignalParMC( 1.65287, 2.063, 0.0727294, 9.45, 0.0836645, 1.94962 );}
if ( binMatched( 4, 0, 30, 2.0, 2.4) )  { setSignalParMC( 1.22102, 2.05874, 0.0874083, 9.4581, 0.0987932, 1.99991 );}
if ( binMatched( 4, 0, 30, 0, 1.2) )    { setSignalParMC( 1.92198, 1.75, 0.0971168, 9.45798, 0.437838, 0.611932 );}
if ( binMatched( 4, 0, 30, 1.2, 2.4) )  { setSignalParMC( 3.73744, 1.45492, 0.0980838, 9.45048, 0.40192, 1.47978 );}

    }
}



bool PSetUpsAndBkg::binMatched( float muonPtCut_, float ptLow_, float ptHigh_, float yLow_, float yHigh_, int cLow_, int cHigh_) {
  if ( (float)muonPtCut_ != (float)muPtCut )
    return false;
  if  ( (float)ptLow_ != (float)ptLow )
    return false;
  if ( (float)ptHigh_ != (float)ptHigh ) 
    return false;
  if ( (float)yLow_ != (float)yLow )
    return false;
  if ( (float)yHigh_ != (float)yHigh ) 
    return false;
  if ( (cLow_>0) &&  ( (int)cLow_ != (int)cLow ) )
    return false;
  if ( (cHigh_>0) && ( (int)cHigh_ != (int)cHigh ) ) 
    return false;

  return true;

    
}

#endif
