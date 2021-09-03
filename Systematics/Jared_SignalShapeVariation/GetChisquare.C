#include "MakeCMSPlot.C"

void GetChisquare(int collId=kPADATA, int whichUpsilon=2, int whichn=2) {
  
  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[3] = {-1.93,0.0,1.93};
  }
  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int ptconst = numptbins+1;
  const int yconst = numybins+1;

  float ptshift, yshift;
  if (whichn==4) {
    ptshift = 0.3;
    yshift = 0.05;
  }
  else if (whichn==3) {
    ptshift = 0;
    yshift = 0;
  }
  else if (whichn==2) {
    ptshift = -0.3;
    yshift = -0.05;
  }
  float shiftedptbins[ptconst];
  float shiftedybins[yconst];
  cout << endl << "{";
  for (int ishift=0; ishift<ptconst; ishift++) {
    shiftedptbins[ishift] = ptbins[ishift]+ptshift;
    cout << shiftedptbins[ishift] << ",";
  }
  cout << "}" << endl;
  for (int ishifty=0; ishifty<yconst; ishifty++) {
    shiftedybins[ishifty] = ybins[ishifty]+yshift;
  }

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";
  TString outfilename = Form("ChisqHisto_%s_%is_fixedn%i.root",strId.Data(),whichUpsilon,whichn);
  TFile outFile(outfilename, "RECREATE");

  //declare histograms
  TH1D* hchisqpt = new TH1D("hchisqpt","chisq vs pt",numptbins,shiftedptbins);
  TH1D* hchisqy = new TH1D("hchisqy","chisq vs y",numybins,shiftedybins);

  float ptLow, ptHigh, yLow, yHigh;
  double chisqdof;

  //pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {
    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    if (collId==kPADATA) {
      yLow = -1.93;
      yHigh = 1.93;
    }
    else if (collId==kPPDATA) {
      yLow = 0.0;
      yHigh = 1.93;
    } 
    //extract chisq value 
    chisqdof = MakeCMSPlot(collId,ptLow,ptHigh,yLow,yHigh);
    cout << "chisq/dof = " << chisqdof << endl;
    hchisqpt->SetBinContent(ipt+1, chisqdof);
  }

  //y loop
  for (int iy = 0; iy<numybins; iy++) {
    ptLow = 0;
    ptHigh = 30;
    yLow = ybins[iy];
    yHigh = ybins[iy+1];
    if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(ybins[iy+1]);
      yHigh = TMath::Abs(ybins[iy]);
    }  
    //extract chisq value
    chisqdof = MakeCMSPlot(collId,ptLow,ptHigh,yLow,yHigh);
    cout << "chisq/dof = " << chisqdof << endl;
    hchisqy->SetBinContent(iy+1, chisqdof);
  }

  //save histograms
  outFile.cd();
  hchisqpt->Write();
  hchisqy->Write();
  outFile.Close();

  //delete hchisqy;
  //delete hchisqpt;
}
