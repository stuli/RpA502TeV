#include "RooSimWSTool.h"
#include "RooSimPdfBuilder.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooKeysPdf.h"
#include "RooUniform.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "TTree.h"
#include "TFile.h"

float mmin = 7, mmax = 14.0;
TString dirname_ = "";
TString treeName = "UpsilonTree_allsign";
/*Only for track rotation*/
TString treeTrkRot = "UpsilonTree_trkRot";

bool buildPdf(RooWorkspace &ws, bool hi, int bkgdModel = 3, bool TrkRotBkgd = 0) {

	double const M1S(9.460);
	double const M2S(10.023);
	double const M3S(10.355);
	RooRealVar * mass = ws.var("invariantMass");
	if (!mass) {
		mass = new RooRealVar("invariantMass", "#mu#mu mass", mmin, mmax, 
				"GeV/c^{2}");
	}
	//   mass->setRange(mmin,mmax);
	RooRealVar mean("mean", "mean", M1S, M1S-0.3, M1S+0.3, "GeV/c^{2}");

	RooConstVar rat2("rat2", "rat2", M2S/M1S);
	RooConstVar rat3("rat3", "rat3", M3S/M1S);

	RooConstVar diff2("diff2", "diff2", M2S-M1S);
	RooConstVar diff3("diff3", "diff3", M3S-M1S);
	RooRealVar mscale("mscale", "mscale", 1.0);


	RooFormulaVar mean1S("mean1S", "@0", RooArgList(mean));
	RooFormulaVar mean2S("mean2S", "@0+@1*@2", RooArgList(mean,diff2,mscale));
	RooFormulaVar mean3S("mean3S", "@0+@1*@2", RooArgList(mean,diff3,mscale));

	RooRealVar sigma1("sigma1", "#sigma_{1}", 0.092, 0.01, 0.3);
	//sigma1.setConstant();
	RooFormulaVar sigma1S("sigma1S", "@0", RooArgList(sigma1));
	RooFormulaVar sigma2S("sigma2S", "@0*@1", RooArgList(sigma1,rat2));
	RooFormulaVar sigma3S("sigma3S", "@0*@1", RooArgList(sigma1,rat3));

	RooRealVar alpha("alpha", "#alpha", 1.6, 0.1, 10.);
	//alpha.setConstant();
	RooRealVar npow("npow", "n_{CB}", 2.3, 0.1, 10.);
	npow.setConstant();

	RooCBShape sig1S("sig1S", "sig1S", *mass, mean1S, sigma1S, alpha, npow);
	RooCBShape sig2S("sig2S", "sig2S", *mass, mean2S, sigma2S, alpha, npow);
	RooCBShape sig3S("sig3S", "sig3S", *mass, mean3S, sigma3S, alpha, npow);

	RooRealVar nsig1("nsig1", "N_{1S}", 100, -1000, 1e5);
	RooRealVar f23("f23", "(2S+3S)/1S", 0.5, 0, 1.3);
	f23.setConstant(false);
	//RooRealVar f2("f2", "2S/1S", 0.2, 0, 1);
	RooRealVar f3o2("f3o2", "3S/2S", 0.6, 0., 1.);
	f3o2.setConstant(false);
	RooRealVar x23("x23", "#chi_{23}", 0.2, -1, 1);
	x23.setConstant(false);
	RooRealVar x3o2("x3o2", "#chi_{3o2}", 0.3, -1, 1);
	x3o2.setConstant(false);
	//RooUniform prior("prior","",x23);
	//ws.import(prior, RooFit::Silence());
	
	RooFormulaVar *nsig2;
	RooFormulaVar *nsig3;
	if (!hi){
		nsig2 = new RooFormulaVar("nsig2", "@0*@1*(1/(1+@2))", RooArgList(nsig1, f23, f3o2));
		nsig3 = new RooFormulaVar("nsig3", "@0*@1*(@2/(1+@2))", RooArgList(nsig1, f23, f3o2));
	}
	else {
		nsig2 = new RooFormulaVar("nsig2", "@0*@1*@2*(1/(1+@3*@4))", RooArgList(nsig1,f23,x23,f3o2,x3o2)); 
    	nsig3 = new RooFormulaVar("nsig3", "@0*@1*@2*(@3*@4/(1+@3*@4))", RooArgList(nsig1,f23,x23,f3o2,x3o2));
	}
	
	//   RooFormulaVar nsig2("nsig2", "@0*@1", RooArgList(nsig1, f2));
	//   RooFormulaVar nsig3("nsig3", "@0*(@1-@2)", RooArgList(nsig1, f23, f2));
	RooRealVar turnOn("turnOn","turnOn", 6., 0., 20.);
	//turnOn.setConstant(false);
	RooRealVar width("width","width", 1., 0., 20.);
	RooRealVar decay("decay","decay", 7., 0., 10.);
	//decay.setConstant(false);
	RooGenericPdf bkgErfExp("bkgLikeSignPdf","bkg",
			"exp(-@0/@3)*(TMath::Erf((@0-@1)/@2)+1)",
			RooArgList(*mass, turnOn, width, decay));

	RooRealVar nbkg("nbkg", "N_{bkg}", 5000., -1000., 1e6);

	RooRealVar bkg_a1("bkg_a1", "a_{1,bkg}", 0., -2., 2.);
	RooRealVar bkg_a2("bkg_a2", "a_{2,bkg}", 0., -2., 2.);
	RooChebychev bkgPoly("bkgPoly", "bkg", *mass, 
			RooArgList(bkg_a1, bkg_a2));
	RooRealVar * nLikeSign = 0;
	RooKeysPdf * LikeSignPdf = 0;
	RooFormulaVar * nPoly = 0;

	TString likeSignCut("(QQsign == QQsign::PlusMinus) && (dataCat == dataCat::");
	if (bkgdModel == 1 || bkgdModel == 2) {
		if (TrkRotBkgd) likeSignCut += "TrkRot)";
		else 
		likeSignCut += (hi) ? "hi)" : "pp)";
		RooDataSet * likeSignData = 
			dynamic_cast<RooDataSet*>(ws.data("data")->reduce(likeSignCut));

		assert(likeSignData);
		LikeSignPdf = new RooKeysPdf("bkgLikeSignPdf", "bkgLikeSignPdf", *mass, 
				*likeSignData, RooKeysPdf::MirrorBoth, 1.4);
		nLikeSign = 
			new RooRealVar("nLikeSign", "nLikeSign", 
					likeSignData->sumEntries(TString::Format("(%s>%0.1f)&&"
							"(%s<%0.1f)",
							mass->GetName(),
							mmin, 
							mass->GetName(),
							mmax)));
			//nLikeSign->setVal(9046);
		delete likeSignData;
		std::cout << "like sign events: " << nLikeSign->getVal() << '\n';
		nPoly =  new RooFormulaVar ("nPoly", "nPoly", "@0-@1",
				RooArgList(nbkg,*nLikeSign));
		RooAddPdf bkg("bkg", "bkg", RooArgList(*LikeSignPdf, bkgPoly),
				RooArgList(*nLikeSign, *nPoly));
	}

	else if (bkgdModel == 6 || bkgdModel == 7) {
        if (TrkRotBkgd) likeSignCut += "TrkRot)";
        else 
		likeSignCut += (hi) ? "hi)" : "pp)";
		RooDataSet * likeSignData =
			dynamic_cast<RooDataSet*>(ws.data("data")->reduce(likeSignCut));

		assert(likeSignData);

		bkgErfExp.fitTo(*likeSignData);
        turnOn.setConstant(true);
        width.setConstant(true);
        decay.setConstant(true);

		nLikeSign =
			new RooRealVar("nLikeSign", "nLikeSign",
					likeSignData->sumEntries(TString::Format("(%s>%0.1f)&&"
							"(%s<%0.1f)",
							mass->GetName(),
							mmin,
							mass->GetName(),
							mmax)));
			//nLikeSign->setVal(9046);
		delete likeSignData;
		std::cout << "like sign events: " << nLikeSign->getVal() << '\n';
		nPoly =  new RooFormulaVar ("nPoly", "nPoly", "@0-@1",
				RooArgList(nbkg,*nLikeSign));
		RooAddPdf bkg("bkg", "bkg", RooArgList(bkgErfExp, bkgPoly),
				RooArgList(*nLikeSign, *nPoly));
	}

	//   RooExtendPdf sig1SN("sig1SN", "sig1SN", sig1SN, nsig1);
	//   RooExtendPdf sig2SN("sig2SN", "sig2SN", sig2SN, nsig2);
	//   RooExtendPdf sig3SN("sig3SN", "sig3SN", sig3SN, nsig3);
	//   RooExtendPdf bkgN("bkgN", "bkgN", bkgN, nbkg);

	// RooAddPdf pdf("pdf", "pdf", RooArgList(sig1SN, sig2SN, sig3SN, bkgN));
	RooArgList pdfs(sig1S, sig2S, sig3S);
	RooArgList norms(nsig1, *nsig2, *nsig3);
	RooAddPdf sig("sig", "sig", pdfs, norms);
	ws.import(sig, RooFit::RenameAllNodes((hi)?"hi":"pp"),
			RooFit::RenameAllVariablesExcept((hi)? "hi": "pp", 
				"npow,invariantMass,"
				//"prior,"
				//"mean,"
				//"turnOn,"
				"f23,f3o2,"
				"x23,x3o2,"
				"alpha,"
				"sigma1"
				), 
			RooFit::RecycleConflictNodes());

	switch (bkgdModel) {
		case 1 : //use RooKeysPdf to smooth the like-sign, then fit unlikesign with keys + pol1
			bkg_a2.setConstant(kTRUE);
			pdfs.add(RooArgList(*LikeSignPdf, bkgPoly));
			norms.add(RooArgList(*nLikeSign, *nPoly));
			break;

		case 2 : //use RooKeysPdf to smooth the like-sign, then fit unlikesign with keys + pol2
			bkg_a2.setConstant(kFALSE);
			pdfs.add(RooArgList(*LikeSignPdf, bkgPoly));
			norms.add(RooArgList(*nLikeSign, *nPoly));
			break;

		case 3 : //use erf * exp to fit the unlike-sign
			bkgErfExp.SetName("bkg");
			pdfs.add(bkgErfExp);
			norms.add(nbkg);
			break;

		case 4 : //use pol2 to fit the unlike-sign
			bkgPoly.SetName("bkg");
			bkg_a2.setConstant(kFALSE);
			pdfs.add(bkgPoly);
			norms.add(nbkg);
			break;

		case 5 : //use pol1 to fit the unlike-sign
			bkgPoly.SetName("bkg");
			bkg_a2.setConstant(kTRUE);
			pdfs.add(bkgPoly);
			norms.add(nbkg);
			break;

		case 6 : //use erf*exp to fit like-sign, then fit unlikesign with keys + pol1
			bkg_a2.setConstant(kTRUE);
			pdfs.add(RooArgList(bkgErfExp, bkgPoly));
			norms.add(RooArgList(*nLikeSign, *nPoly));
			break;

		case 7 : //use erf*exp to fit like-sign, then fit unlikesign with keys + pol2
			bkg_a2.setConstant(kFALSE);
			pdfs.add(RooArgList(bkgErfExp, bkgPoly));
			norms.add(RooArgList(*nLikeSign, *nPoly));
			break;

		default :
			break;
	}	
	//   RooAddPdf pdf(pdfName, pdfName, RooArgList(sig, bkg));

	RooAddPdf pdf("pdf", "pdf", pdfs, norms);

	//pdf.Print("v");
	ws.import(pdf, 
			RooFit::RenameAllNodes((hi)?"hi":"pp"),
			RooFit::RenameAllVariablesExcept((hi)? "hi": "pp", 
				"npow,invariantMass,"
				//"prior,"
				//"turnOn,"
				//"mean,"
				"f23,f3o2,"
				"x23,x3o2,"
				"alpha,"
				"sigma1"
				),
			RooFit::RecycleConflictNodes());

	if (bkgdModel == 1 || bkgdModel == 2 || bkgdModel == 6 || bkgdModel == 7) {
		delete nLikeSign;
		delete LikeSignPdf;
		delete nPoly;
	}
	return true;
}

RooSimultaneous* buildSimPdf(RooWorkspace &ws, RooCategory& dataCat) {

	if (ws.pdf("simPdf"))
		return dynamic_cast<RooSimultaneous *>(ws.pdf("simPdf"));

	RooSimultaneous simPdf("simPdf", "simPdf", dataCat);
	//RooAbsPdf * pdf_hi = ws.pdf("pdf_hi");
	//RooAbsPdf * pdf_pp = ws.pdf("pdf_pp");
	RooAbsPdf * pdf_hi = ws.pdf("model_hi");
	RooAbsPdf * pdf_pp = ws.pdf("model_pp");
	assert(pdf_hi);
	assert(pdf_pp);
	simPdf.addPdf(*pdf_hi, "hi");
	simPdf.addPdf(*pdf_pp, "pp");
	ws.import(simPdf, RooFit::Silence());
	return dynamic_cast<RooSimultaneous *>(ws.pdf("simPdf"));

}



// RooSimultaneous*  buildNullPdf(RooWorkspace &ws, RooCategory const &dataCat) {
//   if (!ws.pdf("pdf") && buildPdf(ws)) {
//     std::cout << "cannot get pdf.\n";
//     return 0;
//   }

//   RooSimPdfBuilder sb(RooArgSet(*(ws.pdf("pdf"))));
//   RooArgSet * config = sb.createProtoBuildConfig();
//   config->setStringValue("physModels", "pdf");
//   config->setStringValue("splitCats", dataCat.GetName());
//   config->setStringValue("pdf", TString(dataCat.GetName()) + " : f2,nsig1,nbkg,bkg_a1,bkg_a2");

//   RooArgSet deps(*(ws.var("invariantMass")),dataCat);

//   RooSimultaneous * sim = sb.buildPdf(*config, deps);
//   sim->SetName("simNullPdf");

//   ws.import(*sim, RooFit::RenameConflictNodes("null"),
// 	    RooFit::RenameAllVariablesExcept("null", "invariantMass,dataCat"));

//   return (RooSimultaneous *)ws.pdf("simNullPdf");
// }

bool readData(RooWorkspace &ws, TString HIfilename, TString ppFilename,
		TString extraCut = "") {

	RooRealVar * mass = ws.var("invariantMass");
	if (!mass) {
		mass = new RooRealVar("invariantMass", "#mu#mu mass", mmin, mmax, 
				"GeV/c^{2}");
		mass->setBins(70);
	}
	RooRealVar muppt("muPlusPt" ,"#mu+ pt",2,20,"GeV/c"); 
	RooRealVar mumpt("muMinusPt","#mu- pt",2,20,"GeV/c"); 
	RooRealVar upsPt("upsPt","p_{T}(#Upsilon)",0.,"GeV/c");
	RooRealVar vProb("vProb","vProb",0.05,1);
	//   RooRealVar upsEta("upsEta","#eta(#Upsilon)",0.,"");
	RooRealVar upsRapidity("upsRapidity", "upsRapidity", 0.);
	RooCategory QQsign("QQsign", "QQsign");
	QQsign.defineType("PlusMinus", 0);
	QQsign.defineType("PlusPlus", 1);
	QQsign.defineType("MinusMinus", 2);
	RooRealVar Centrality("Centrality", "Centrality", 0.);

	RooArgSet cols(*mass, muppt, mumpt, upsPt, vProb, upsRapidity, QQsign, Centrality);

	//import HI data
	TFile * hifile = TFile::Open(HIfilename);
	TTree * tree;
	TString dirTree = treeName;
	if (dirname_.Length() > 0)
		dirTree = dirname_ + "/" + treeName;
	hifile->GetObject(dirTree, tree);
	assert(tree);
	TFile temp("DeleteMe.root", "recreate");
	if (extraCut.Length() > 0) {
		TTree * tree2 = tree->CopyTree(extraCut);
		//TTree * tree2 = tree->CopyTree(extraCut + "&& Centrality>=2 && Centrality < 4");
		delete tree;
		tree = tree2;
	}

	RooDataSet hidata("hidata", "hidata", tree, cols);
	delete tree;
	delete hifile;
/*Only for track rotation*/
	//improt tack rotation data
    TFile * hiTrkRotfile = TFile::Open(HIfilename);
    tree = 0;
    TString dirTree_TrkRot = treeTrkRot;
    if (dirname_.Length() > 0) 
		dirTree_TrkRot = dirname_ + "/" + treeTrkRot;
    hiTrkRotfile->GetObject(dirTree_TrkRot, tree);
    assert(tree);
	temp.cd();
    if (extraCut.Length() > 0) {
        TTree * tree2 = tree->CopyTree(extraCut);
        delete tree;
        tree = tree2;
    }   

    RooDataSet hidata_TrkRot("hidata_TrkRot", "hidata_TrkRot", tree, cols);
    delete tree;
    delete hiTrkRotfile;

	//import pp data
	TFile * ppfile = TFile::Open(ppFilename);
	tree = 0;
	ppfile->GetObject(dirTree, tree);
	assert(tree);
	temp.cd();
	if (extraCut.Length() > 0) {
		TTree * tree2 = tree->CopyTree(extraCut);
		delete tree;
		tree = tree2;
	}

	RooDataSet ppdata("ppdata", "ppdata", tree, cols);
	delete tree;
	delete ppfile;

	RooCategory dataCat("dataCat", "dataCat");
	dataCat.defineType("hi");
	dataCat.defineType("pp");
	dataCat.defineType("TrkRot");

	RooDataSet data("data", "data", cols, RooFit::Index(dataCat),
/*Only for track rotation*/
		RooFit::Import("hi", hidata), RooFit::Import("pp", ppdata), RooFit::Import("TrkRot", hidata_TrkRot));
//          RooFit::Import("hi", hidata), RooFit::Import("pp", ppdata));
	  data.Print("v");
	return ws.import(data/*, RooFit::RecycleConflictNodes(), RooFit::Silence()*/);  
}

double computeRatio(RooRealVar& x, RooRealVar& y) {
	assert(y.getVal() != 0);
	return x.getVal()/y.getVal();
}

double computeRatioError(RooRealVar& x, RooRealVar& y, 
		double correlation = 0.) {
	double err2 = x.getError()*x.getError()/x.getVal()/x.getVal() +
		y.getError()*y.getError()/y.getVal()/y.getVal() - 
		2.*x.getError()*y.getError()/x.getVal()/y.getVal()*correlation;

	return fabs(computeRatio(x,y))*sqrt(err2);
}

/*
void modelconfig(RooWorkspace &ws){
  RooWorkspace * pWs = new RooWorkspace("myWS");
  pWs = &ws;
  // fix all other variables in model:
  // everything except observables, POI, and nuisance parameters
  // must be constant
  pWs->var("alpha")->setConstant(true);
  pWs->var("bkg_a1_pp")->setConstant(true);
  pWs->var("bkg_a2_pp")->setConstant(true);
  pWs->var("decay_hi")->setConstant(true);
  pWs->var("f2")->setConstant(true);
  pWs->var("f23")->setConstant(true);
  pWs->var("mean_hi")->setConstant(true);
  pWs->var("mean_pp")->setConstant(true);
  pWs->var("nbkg_hi")->setConstant(true);
  pWs->var("nbkg_pp")->setConstant(true);
  pWs->var("sigma1")->setConstant(true);
  pWs->var("turnOn_hi")->setConstant(true);
  pWs->var("width_hi")->setConstant(true);
  pWs->var("x2")->setConstant(true);
  pWs->var("x23")->setConstant(true);
  RooArgSet fixed("fixed");
  fixed.add( *(pWs->var("alpha")));
  fixed.add( *(pWs->var("bkg_a1_pp")));
  fixed.add( *(pWs->var("bkg_a2_pp")));
  fixed.add( *(pWs->var("decay_hi")));
  fixed.add( *(pWs->var("f2")));
  fixed.add( *(pWs->var("f23")));
  fixed.add( *(pWs->var("mean_hi")));
  fixed.add( *(pWs->var("mean_pp")));
  fixed.add( *(pWs->var("nbkg_hi")));
  fixed.add( *(pWs->var("nbkg_pp")));
  fixed.add( *(pWs->var("sigma1")));
  fixed.add( *(pWs->var("turnOn_hi")));
  fixed.add( *(pWs->var("width_hi")));
  fixed.add( *(pWs->var("x2")));
  fixed.add( *(pWs->var("x23")));
  // create signal+background Model Config
  RooStats::ModelConfig sbHypo("SbHypo");
  sbHypo.SetWorkspace( *pWs);
  sbHypo.SetPdf( *pWs->pdf("simPdf") );
  sbHypo.SetObservables( *pWs->var("invariantMass") );
  //sbHypo.SetGlobalObservables( globalObs );
  sbHypo.SetParametersOfInterest( *pWs->var("x23") );
  //sbHypo.SetNuisanceParameters( nuis );
  sbHypo.SetPriorPdf( *pWs->pdf("prior") ); // this is optional
  // import ModelConfig into workspace
  pWs->import( sbHypo );
  pWs->SaveAs("workspace.root");
  return;
}
*/
