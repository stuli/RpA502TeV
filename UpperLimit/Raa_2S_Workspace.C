#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "TAxis.h"
#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TROOT.h>
#include <fstream>
#include <TGraph.h>
#include "TMath.h"
#include "TF1.h"
#include <RooMinuit.h>

#include "test_combine.C"

using namespace RooFit;
using namespace RooStats;

//void SingleRt_2S_1S_Workspace(const char* name_pbpb="fitResults/fitresults_upsilon_DoubleCB_AA_DATA_pt4.0-6.0_y0.0-1.2_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root", const char* name_pp="fitResults/fitresults_upsilon_DoubleCB_PP_DATA_pt4.0-6.0_y0.0-1.2_muPt4.0.root", const char* name_out="fitResults_combo_reducedDS_2S.root"){
void SingleRt_2S_1S_Workspace(const char* name_pbpb="fitResults/fitresults_upsilon_DoubleCB_AA_DATA_pt15.0-30.0_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI-1.root", const char* name_pp="fitResults/fitresults_upsilon_DoubleCB_PP_DATA_pt15.0-30.0_y0.0-2.4_muPt4.0-1.root", const char* name_out="fitResults_combo_reducedDS_2S_pt15_30.root"){

   RooWorkspace * ws = test_combine(name_pbpb, name_pp);

   RooAbsData * data = ws->data("data");

   RooRealVar* raa2 = new RooRealVar("raa2","Ratio(#Upsilon(2S)_{PbPb}/#Upsilon(2S)_{pp})",.5,-1,1);
   //RooRealVar* dubRat = new RooRealVar("dubRat","DoubleRatio(#Upsilon (3S))",0.5,0.,1.);
   RooRealVar* leftEdge = new RooRealVar("leftEdge","leftEdge",0);
   RooRealVar* rightEdge = new RooRealVar("rightEdge","rightEdge",1.);
   RooGenericPdf step("step", "step", "(@0 >= @1) && (@0 < @2)", RooArgList(*raa2, *leftEdge, *rightEdge));
   ws->import(step);
   ws->factory( "Uniform::flat(raa2)" );//ws->factory( "Uniform::flat(dubRat)" );

   RooRealVar* nsig1_hi = ws->var("nSig1s_hi");//RooRealVar* nsig1_hi = ws->var("N_{#varUpsilon(1S)}_hi");
   RooRealVar* nsig3_hi = ws->var("nSig3s_hi");
   //RooRealVar* nsig3_1_hi = ws->var("nSig3s_hi");
   RooRealVar* nsig2_pp = ws->var("nSig2s_pp");
   cout << "!!!!!!!# of N#Upsilon2S_pp!!!!!!! : " << nsig2_pp << endl;
   RooFormulaVar nsig2_hi_modified("nsig2_hi_modified", "@0*@1", RooArgList(*raa2, *nsig2_pp));
   ws->import(nsig2_hi_modified);
   // background yield with systematics
   ws->factory( "nbkg_hi_kappa[1.1]" );
   ws->factory( "expr::alpha_nbkg_hi('pow(nbkg_hi_kappa,beta_nbkg_hi)',nbkg_hi_kappa,beta_nbkg_hi[0,-5,5])" );
   ws->factory( "SUM::nbkg_hi_nom(alpha_nbkg_hi*bkgHighPt_hi)" );//ws->factory( "SUM::nbkg_hi_nom(alpha_nbkg_hi*bkgPdf_hi)" );
   ws->factory( "Gaussian::constr_nbkg_hi(beta_nbkg_hi,glob_nbkg_hi[0.,-5.,5.],1)" );
   RooAbsPdf* sig1S_hi = ws->pdf("cb1s_hi");
   RooAbsPdf* sig2S_hi = ws->pdf("cb2s_hi");
   RooAbsPdf* sig3S_hi = ws->pdf("cb3s_hi");
   RooAbsPdf* LSBackground_hi = ws->pdf("nbkg_hi_nom");
   //RooFormulaVar *nsig2_hi = (RooFormulaVar*)ws->function("N_{#varUpsilon(2S)}_hi");
   RooFormulaVar *nsig2_hi = (RooFormulaVar*)ws->function("nsig2_hi_modified"); 
   RooRealVar* norm_nbkg_hi = ws->var("nBkg_hi");
   //pp
   ws->factory( "nbkg_pp_kappa[1.03]" );
   ws->factory( "expr::alpha_nbkg_pp('pow(nbkg_pp_kappa,beta_nbkg_pp)',nbkg_pp_kappa,beta_nbkg_pp[0.,-5.,5.])" );
   ws->factory( "SUM::nbkg_pp_nom(alpha_nbkg_pp*bkgHighPt_pp)" );//ws->factory( "SUM::nbkg_pp_nom(alpha_nbkg_pp*bkgPdf_pp)" );
   ws->factory( "Gaussian::constr_nbkg_pp(beta_nbkg_pp,glob_nbkg_pp[0.,-5.,5.],1)" );
   RooAbsPdf* sig1S_pp = ws->pdf("cb1s_pp");//RooAbsPdf* sig1S_pp = ws->pdf("sig1S_pp");
   RooAbsPdf* sig2S_pp = ws->pdf("cb2s_pp");//RooAbsPdf* sig2S_pp = ws->pdf("sig2S_pp");
   RooAbsPdf* sig3S_pp = ws->pdf("cb3s_pp");//RooAbsPdf* sig3S_pp = ws->pdf("sig3S_pp");
   RooAbsPdf* LSBackground_pp = ws->pdf("nbkg_pp_nom");
   RooRealVar* nsig1_pp = ws->var("nSig1s_pp");//RooRealVar* nsig1_pp = ws->var("N_{#varUpsilon(1S)}_pp");
   RooRealVar* nsig3_pp = ws->var("nSig3s_pp");
   RooRealVar* norm_nbkg_pp = ws->var("nBkg_pp");//RooRealVar* norm_nbkg_pp = ws->var("n_{Bkgd}_pp");

   RooArgList pdfs_pp( *sig1S_pp,*sig2S_pp,*sig3S_pp, *LSBackground_pp);
   RooArgList norms_pp( *nsig1_pp,*nsig2_pp,*nsig3_pp,*norm_nbkg_pp);

   RooArgList pdfs_hi( *sig1S_hi,*sig2S_hi,*sig3S_hi, *LSBackground_hi);
   RooArgList norms_hi(*nsig1_hi,*nsig2_hi,*nsig3_hi, *norm_nbkg_hi);
   ////////////////////////////////////////////////////////////////////////////////
   RooAddPdf model_num("model_num", "model_num", pdfs_hi,norms_hi); 
   ws->import(model_num);
   ws->factory("PROD::models_hi(model_num, constr_nbkg_hi)");
   RooAddPdf model_den("model_den", "model_den", pdfs_pp,norms_pp); 
   ws->import(model_den);
   ws->factory("PROD::models_pp(model_den, constr_nbkg_pp)");
   ws->factory("SIMUL::joint(dataCat,hi=models_hi,pp=models_pp)");
   /////////////////////////////////////////////////////////////////////
   RooRealVar * pObs = ws->var("mass"); // get the pointer to the observable
   //RooRealVar * pObs = ws->var("invariantMass"); // get the pointer to the observable
   RooArgSet obs("observables");
   obs.add(*pObs);
   obs.add( *ws->cat("dataCat"));    
   /////////////////////////////////////////////////////////////////////
   //ws->var("glob_lumipp")->setConstant(true);
   //ws->var("glob_Taa")->setConstant(true);
   //ws->var("glob_effRat")->setConstant(true);
   ws->var("glob_nbkg_pp")->setConstant(true);
   ws->var("glob_nbkg_hi")->setConstant(true);
   RooArgSet globalObs("global_obs");
   //globalObs.add( *ws->var("glob_lumipp") );
   //globalObs.add( *ws->var("glob_Taa") );
   //globalObs.add( *ws->var("glob_effRat") );
   globalObs.add( *ws->var("glob_nbkg_pp") );
   globalObs.add( *ws->var("glob_nbkg_hi") );

   RooArgSet poi("poi");
   poi.add( *ws->var("raa2") );//poi.add( *ws->var("dubRat") );

   // create set of nuisance parameters
   RooArgSet nuis("nuis");
   //nuis.add( *ws->var("beta_lumipp") );
   nuis.add( *ws->var("beta_nbkg_hi") );
   nuis.add( *ws->var("beta_nbkg_pp") );
   //nuis.add( *ws->var("beta_Taa") );
   //nuis.add( *ws->var("beta_effRat") );

   ws->var("alpha1s_1_hi")->setConstant(true);
   ws->var("alpha1s_2_hi")->setConstant(true);
   ws->var("alpha1s_1_pp")->setConstant(true);
   ws->var("alpha1s_2_pp")->setConstant(true);
   ws->var("alpha2s_1_hi")->setConstant(true);
   ws->var("alpha2s_2_hi")->setConstant(true);
   ws->var("alpha2s_1_pp")->setConstant(true);
   ws->var("alpha2s_2_pp")->setConstant(true);
   ws->var("alpha3s_1_hi")->setConstant(true);
   ws->var("alpha3s_2_hi")->setConstant(true);
   ws->var("alpha3s_1_pp")->setConstant(true);
   ws->var("alpha3s_2_pp")->setConstant(true);
   ws->var("sigma1s_1_hi")->setConstant(true);
   ws->var("sigma1s_2_hi")->setConstant(true);
   ws->var("sigma1s_1_pp")->setConstant(true);
   ws->var("sigma1s_2_pp")->setConstant(true);
   ws->var("sigma2s_1_hi")->setConstant(true);
   ws->var("sigma2s_2_hi")->setConstant(true);
   ws->var("sigma2s_1_pp")->setConstant(true);
   ws->var("sigma2s_2_pp")->setConstant(true);
   ws->var("sigma3s_1_hi")->setConstant(true);
   ws->var("sigma3s_2_hi")->setConstant(true);
   ws->var("sigma3s_1_pp")->setConstant(true);
   ws->var("sigma3s_2_pp")->setConstant(true);
   //ws->var("Centrality")->setConstant(true);
   ws->var("nSig1s_hi")->setConstant(true);//ws->var("N_{#varUpsilon(1S)}_hi")->setConstant(true);
   ws->var("nSig1s_pp")->setConstant(true);//ws->var("N_{#varUpsilon(1S)}_pp")->setConstant(true);
   //ws->var("nSig2s_hi")->setConstant(true);//ws->var("N_{#varUpsilon(2S)}_hi")->setConstant(true);
   ws->var("nSig2s_pp")->setConstant(true);//ws->var("N_{#varUpsilon(2S)}_pp")->setConstant(true);
   ws->var("nSig3s_hi")->setConstant(true);//ws->var("N_{#varUpsilon(3S)}_hi")->setConstant(true);
   //ws->var("nsig3_hi_modified")->setConstant(true);//ws->var("N_{#varUpsilon(3S)}_hi")->setConstant(true);
   ws->var("nSig3s_pp")->setConstant(true);//ws->var("N_{#varUpsilon(3S)}_pp")->setConstant(true);
   //ws->var("Nmb_hi")->setConstant(true);
   //ws->var("Taa_hi")->setConstant(true);
   //ws->var("Taa_kappa")->setConstant(true);
   //ws->var("beta_Taa")->setConstant(true);
   //ws->var("beta_effRat")->setConstant(true);
   //ws->var("beta_lumipp")->setConstant(true);
   ws->var("beta_nbkg_hi")->setConstant(true);
   ws->var("beta_nbkg_pp")->setConstant(true);
   ws->var("#lambda_hi")->setConstant(true);//ws->var("decay_hi")->setConstant(true);
   ws->var("#lambda_pp")->setConstant(true);//ws->var("decay_pp")->setConstant(true);
   //ws->var("effRat_kappa")->setConstant(true);
   //ws->var("glob_Taa")->setConstant(true);
   //ws->var("glob_effRat")->setConstant(true);
   //ws->var("glob_lumipp")->setConstant(true);
   ws->var("glob_nbkg_hi")->setConstant(true);
   ws->var("glob_nbkg_pp")->setConstant(true);
   //ws->var("invariantMass")->setConstant(true);
   ws->var("mass")->setConstant(true);
   //ws->var("lumipp_hi")->setConstant(true);
   //ws->var("lumipp_kappa")->setConstant(true);
   ws->var("m_{#Upsilon(1S)}_hi")->setConstant(true);//ws->var("m_{ #varUpsilon(1S)}_hi")->setConstant(true);
   ws->var("m_{#Upsilon(1S)}_pp")->setConstant(true);//ws->var("m_{ #varUpsilon(1S)}_pp")->setConstant(true);
   //ws->var("muMinusPt")->setConstant(true);
   //ws->var("muPlusPt")->setConstant(true);
   ws->var("nBkg_hi")->setConstant(true);
   ws->var("nBkg_pp")->setConstant(true);
   ws->var("nbkg_hi_kappa")->setConstant(true);
   ws->var("nbkg_pp_kappa")->setConstant(true);
   //ws->var("n_{CB}_hi")->setConstant(true);
   //ws->var("n_{CB}_pp")->setConstant(true);
   ws->var("n1s_1_hi")->setConstant(true);
   ws->var("n1s_2_hi")->setConstant(true);
   ws->var("n1s_1_pp")->setConstant(true);
   ws->var("n1s_2_pp")->setConstant(true);
   ws->var("n2s_1_hi")->setConstant(true);
   ws->var("n2s_2_hi")->setConstant(true);
   ws->var("n2s_1_pp")->setConstant(true);
   ws->var("n2s_2_pp")->setConstant(true);
   ws->var("n3s_1_hi")->setConstant(true);
   ws->var("n3s_2_hi")->setConstant(true);
   ws->var("n3s_1_pp")->setConstant(true);
   ws->var("n3s_2_pp")->setConstant(true);
   //ws->var("npow")->setConstant(true);
   //ws->var("raa3")->setConstant(true);
   ws->var("leftEdge")->setConstant(true);
   ws->var("rightEdge")->setConstant(true);
   ws->var("f1s_hi")->setConstant(true);//ws->var("sigmaFraction_hi")->setConstant(true);
   ws->var("f1s_pp")->setConstant(true);//ws->var("sigmaFraction_pp")->setConstant(true);
   //ws->var("#mu_hi")->setConstant(true);//ws->var("turnOn_hi")->setConstant(true);
   //ws->var("#mu_pp")->setConstant(true);//ws->var("turnOn_pp")->setConstant(true);
   ws->var("pt")->setConstant(true);//ws->var("dimuPt")->setConstant(true); //ws->var("upsPt")->setConstant(true);
   //ws->("dimuRapidity")->setConstant(true); //ws->var("dimuRapidity")->setConstant(true); //ws->var("upsRapidity")->setConstant(true);
   //ws->var("vProb")->setConstant(true);
   //ws->var("#sigma_hi")->setConstant(true);//ws->var("width_hi")->setConstant(true);
   //ws->var("#sigma_pp")->setConstant(true);//ws->var("width_pp")->setConstant(true);
   //ws->var("x3raw")->setConstant(true);
   //RooArgSet fixed_again("fixed_again");
   //fixed_again.add( *ws->var("leftEdge") );
   //fixed_again.add( *ws->var("rightEdge") );
   //fixed_again.add( *ws->var("Taa_hi") );
   //fixed_again.add( *ws->var("Nmb_hi") );
   //fixed_again.add( *ws->var("lumipp_hi") );
   //fixed_again.add( *ws->var("effRat1_hi") );
   //fixed_again.add( *ws->var("effRat2_hi") );
   //fixed_again.add( *ws->var("effRat3_hi") );
   //fixed_again.add( *ws->var("nsig3_pp") );
   //fixed_again.add( *ws->var("nsig1_pp") );
   //fixed_again.add( *ws->var("nbkg_hi") );
   //fixed_again.add( *ws->var("alpha") );
   //fixed_again.add( *ws->var("nbkg_kappa") );
   //fixed_again.add( *ws->var("Taa_kappa") );
   //fixed_again.add( *ws->var("lumipp_kappa") );
   //fixed_again.add( *ws->var("mean_hi") );
   //fixed_again.add( *ws->var("mean_pp") );
   //fixed_again.add( *ws->var("width_hi") );
   //fixed_again.add( *ws->var("turnOn_hi") );
   //fixed_again.add( *ws->var("bkg_a1_pp") );
   //fixed_again.add( *ws->var("bkg_a2_pp") );
   //fixed_again.add( *ws->var("decay_hi") );
   //fixed_again.add( *ws->var("raa1") );
   //fixed_again.add( *ws->var("raa2") );
   //fixed_again.add( *ws->var("nsig2_pp") );
   //fixed_again.add( *ws->var("sigma1") );
   //fixed_again.add( *ws->var("nbkg_pp") );
   //fixed_again.add( *ws->var("npow") );
   //fixed_again.add( *ws->var("muPlusPt") );
   //fixed_again.add( *ws->var("muMinusPt") );
   //fixed_again.add( *ws->var("mscale_hi") );
   //fixed_again.add( *ws->var("mscale_pp") );

   //create signal+background Model Config
   RooStats::ModelConfig sbHypo("SbHypo");
   sbHypo.SetWorkspace( *ws );
   sbHypo.SetPdf( *ws->pdf("joint") );
   sbHypo.SetObservables( obs );
   sbHypo.SetGlobalObservables( globalObs );
   sbHypo.SetParametersOfInterest( poi );
   sbHypo.SetNuisanceParameters( nuis );
   sbHypo.SetPriorPdf( *ws->pdf("step") ); // this is optional

   // ws->Print();
   /////////////////////////////////////////////////////////////////////
   RooAbsReal * pNll = sbHypo.GetPdf()->createNLL( *data,NumCPU(4) );
   RooMinuit(*pNll).migrad(); // minimize likelihood wrt all parameters before making plots
   RooPlot *framepoi = ((RooRealVar *)poi.first())->frame(Bins(10),Range(0.,.3),Title("LL and profileLL in Yield(#Upsilon(2S))"));
   pNll->plotOn(framepoi,ShiftToZero());
   
   RooAbsReal * pProfile = pNll->createProfile( globalObs ); // do not profile global observables
   pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
   pProfile->plotOn(framepoi,LineColor(kRed));
   framepoi->SetMinimum(0);
   framepoi->SetMaximum(3);
   TCanvas *cpoi = new TCanvas();
   cpoi->cd(); framepoi->Draw();
   cpoi->SaveAs("cpoi_reducedDS_2S.pdf");

   ((RooRealVar *)poi.first())->setMin(0.);
   RooArgSet * pPoiAndNuisance = new RooArgSet("poiAndNuisance");
   //pPoiAndNuisance->add(*sbHypo.GetNuisanceParameters());
   //pPoiAndNuisance->add(*sbHypo.GetParametersOfInterest());
   pPoiAndNuisance->add( nuis );
   pPoiAndNuisance->add( poi );
   sbHypo.SetSnapshot(*pPoiAndNuisance);

   RooPlot* xframeSB = pObs->frame(Title("SBhypo"));
   data->plotOn(xframeSB,Cut("dataCat==dataCat::hi"));
   RooAbsPdf *pdfSB = sbHypo.GetPdf();
   RooCategory *dataCat = ws->cat("dataCat");
   pdfSB->plotOn(xframeSB,Slice(*dataCat,"hi"),ProjWData(*dataCat,*data));
   TCanvas *c1 = new TCanvas();
   c1->cd(); xframeSB->Draw();
   c1->SaveAs("c1_reducedDS_2S.pdf");

   delete pProfile;
   delete pNll;
   delete pPoiAndNuisance;
   ws->import( sbHypo );
   /////////////////////////////////////////////////////////////////////
   RooStats::ModelConfig bHypo = sbHypo;
   bHypo.SetName("BHypo");
   bHypo.SetWorkspace(*ws);
   pNll = bHypo.GetPdf()->createNLL( *data,NumCPU(4) );
   RooArgSet poiAndGlobalObs("poiAndGlobalObs");
   poiAndGlobalObs.add( poi );
   poiAndGlobalObs.add( globalObs );
   pProfile = pNll->createProfile( poiAndGlobalObs ); // do not profile POI and global observables
   ((RooRealVar *)poi.first())->setVal( 0 );  // set raa3=0 here
   pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
   pPoiAndNuisance = new RooArgSet( "poiAndNuisance" );
   pPoiAndNuisance->add( nuis );
   pPoiAndNuisance->add( poi );
   bHypo.SetSnapshot(*pPoiAndNuisance);

   RooPlot* xframeB = pObs->frame(Title("Bhypo"));
   data->plotOn(xframeB,Cut("dataCat==dataCat::hi"));
   RooAbsPdf *pdfB = bHypo.GetPdf();
   pdfB->plotOn(xframeB,Slice(*dataCat,"hi"),ProjWData(*dataCat,*data));
   TCanvas *c2 = new TCanvas();
   c2->cd(); xframeB->Draw();
   c2->SaveAs("c2_reducedDS_2S.pdf");

   delete pProfile;
   delete pNll;
   delete pPoiAndNuisance;

   // import model config into workspace
   bHypo.SetWorkspace(*ws);
   ws->import( bHypo );
   /////////////////////////////////////////////////////////////////////
   ws->Print();
   bHypo.Print();
   sbHypo.Print();

   // save workspace to file
   ws -> SaveAs(name_out);

   return;
}
