/////Modified by Graham. Uniform distributions for mass, pT, and rapidity.

//
// Modeled after the PhaseSpace.C macro.
// This will use the TGenPhaseSpace Class
// to do a 2-body Upsilon Decay into mu+mu-
// Author: Manuel Calderon de la Barca Sanchez
void PhaseSpaceUpsilonMuMu_flat() {

  
  //if (!gROOT->GetClass("TGenPhaseSpace")) gSystem.Load("libPhysics");

   //(Momentum, Energy units are Gev/C, GeV)
  // double UpsilonMass = 9.460; // GeV/c^2
  // TLorentzVector UpsilonAtRest(0.0, 0.0, 0.0, UpsilonMass);
   // TLorentzVector beam(0.0, 0.0, .65, .65);
   // TLorentzVector W = beam + target;


  // Specify the massess of the decay products.
  // For this example, we will use a 2-body decay, so we will need
  // an array of size 2.  The mass of the muon is 0.105658 GeV/c^2
  double muMass = 0.105658; //
  const int nBodies = 2;
  Double_t masses[nBodies] = { muMass, muMass } ;

  // Set up the TGenPhaseSpace class. Constructor does not take any arguments.
  // Use the SetDecay class to control the decay kinematics.
  // First argument is 4-momentum of the decay particle, and should be a TLorentzVector
  // 2nd argument is the number of decay bodies (2-body, 3-body, etc), in this example it is 2.
  // 3d argument is the array of decay product masses
  
   TGenPhaseSpace event;

   TFile* outFile = new TFile("upsilonFlatDimuonMass2BodyNtuple.root","RECREATE");
   
   TNtuple *upsilonNtuple = new TNtuple("upsilonNtuple","Upsilon 2-body decay","UpsM:UpsPt:UpsRap:UpsP:UpsPhi:MuPPt:MuPPhi:MuPy:MuPE:MuMPt:MuMPhi:MuMy:MuME:cosTheta");

   // Upsilon pT distribution for a realistic case
   // Parameters are from a fit to the Upsilon CMS data at 7 TeV from pp collisions
   TF1* upsPtFit = new TF1("upsPtFit","[0]*x/(exp(x/[1])+1)",0,5);
   upsPtFit->SetParameters(100,2.71); // GeV/c^2
   upsPtFit->SetParNames("A","T");

   //TF1* upsRapidityFunc = new TF1("upsRapidityFunc","gaus(0)",-5,5);
   TF1* upsRapidityFunc = new TF1("upsRapidityFunc","gaus(0)+gaus(3)",-5,5);
   upsRapidityFunc->SetParameters(70,0.0,1.0,550,0.0,0.23);
   upsRapidityFunc->SetParNames("A1","mean1","sigma1","A2","mean2","sigma2");
   //TRandom3* randomGen = gRandom;
   gRandom->SetSeed(31415927);
   
   /*TF1* upsMassFunc = new TF1("upsMassFunc","[0]*exp(-x/[1])",6,20);
   upsMassFunc->SetParameters(100,3.4);
   upsMassFunc->SetParNames("A","#lambda");*/
   
   int nEvents = 20000000;
   int nInterval = nEvents/100;
   for (Int_t n=0;n<nEvents;n++) {

     // We can generate a pseudo-realistic distribution of Upsilons using
     // a Gaussian in rapidity and a power-law in pT
     // the phi distribution is uniform.

     double UpsilonMass = gRandom->Uniform(6,20);
	 //double UpsilonMass = upsMassFunc->GetRandom();

     //double UpsPt = upsPtFit->GetRandom();
	 double UpsPt = gRandom->Uniform(0,6);
     double UpsPhi = gRandom->Uniform(2.*TMath::Pi());
     //double UpsRap = upsRapidityFunc->GetRandom();
	 double UpsRap = gRandom->Uniform(-5,5);
     double UpsPx,UpsPy,UpsMt,UpsPz,UpsE, UpsP;
     UpsMt = sqrt(UpsilonMass*UpsilonMass + UpsPt * UpsPt);
     UpsPx = UpsPt*cos(UpsPhi);
     UpsPy = UpsPt*sin(UpsPhi);
     UpsPz = UpsMt * sinh(UpsRap);
     UpsE  = UpsMt * cosh(UpsRap);
     UpsP  = sqrt(UpsE * UpsE - UpsilonMass*UpsilonMass);
     TLorentzVector Upsilon4Vec(UpsPx,UpsPy,UpsPz,UpsE);
     
     //event.SetDecay(UpsilonAtRest, nBodies, masses);
     event.SetDecay(Upsilon4Vec, nBodies, masses);
     Double_t weight = event.Generate();
     
     TLorentzVector *pMuPlus  = event.GetDecay(0);
     TLorentzVector *pMuMinus = event.GetDecay(1);
     
     //double UpsPt = UpsilonAtRest.Pt();
     //double UpsRap = UpsilonAtRest.Rapidity();
     double MuPPt = pMuPlus->Pt();
     double MuPPhi = pMuPlus->Phi();
     double MuPy = pMuPlus->Rapidity();
     double MuPE = pMuPlus->E();
     double MuMPt = pMuMinus->Pt();
     double MuMPhi = pMuMinus->Phi();
     double MuMy = pMuMinus->Rapidity();
     double MuME = pMuMinus->E();
     double cosTheta = cos(pMuPlus->Angle(pMuMinus->Vect()));
     upsilonNtuple->Fill(UpsilonMass,UpsPt,UpsRap,UpsP,UpsPhi,MuPPt,MuPPhi,MuPy,MuPE,MuMPt,MuMPhi,MuMy,MuME,cosTheta);
     if (n%nInterval == 0) {

       cout << "event n= " << n << "  (" << n/nEvents*100 << "%)" << endl;
       cout << "   UpsPt " << UpsPt << endl;
       cout << "   UpsE  " << UpsE << endl;
       cout << "   UpsMt  " << UpsMt << endl;
       cout << "   UpsRap  " << UpsRap << endl;
       cout << "   cosTheta= " << cosTheta << endl;
       cout << "   MuPE= " << MuPE << endl;
       cout << "   MuME= " << MuME << endl;
       cout << "   weight  " << weight << endl;
     }
   }// event loop
   outFile->Write();
   outFile->Close();
   return;
   
}
