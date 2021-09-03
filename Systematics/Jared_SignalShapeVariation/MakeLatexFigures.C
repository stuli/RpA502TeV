#include "../../HeaderFiles/cutsAndBin.h"

void MakeLatexFigures(int whichUpsilon=1) {

  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[11] = {-2.87,-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[7] = {-2.87,-2.4,-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[5] = {-2.87,-2.4,-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numtot = numptbins + numybins;

  float ptLow, ptHigh, yLow, yHigh, yLowPP, yHighPP;
  TString binvar, caption;

  for (int i = -1; i<numtot; i++) {

    if (i<numptbins){
      yLow = -1.93;
      yHigh = 1.93;
      yLowPP = 0.00;
      yHighPP = 1.93;
      if (i<0) {
        cout << "\%INTEGRATED BIN" << endl;
        ptLow = 0;
        ptHigh = 30;
      }
      else {
        ptLow = ptbins[i];
        ptHigh = ptbins[i+1];
      }
      caption = Form("	\\caption{Double Crystal Ball Fit to \\pp\\ (left) and pPb (right) for $\\pt~[\\GeVc]$ $\\in$ [%.1f-%.1f].",ptLow,ptHigh);
    }
    else {
      binvar = "$y$";
      ptLow = 0;
      ptHigh = 30;
      yLow = ybins[i-numptbins];
      yHigh = ybins[i-numptbins+1];
      if (yLow<0) {
        yLowPP = TMath::Abs(yHigh);
        yHighPP = TMath::Abs(yLow);
      }
      else {
        yLowPP = yLow;
        yHighPP = yHigh;
      }
      caption = Form("	\\caption{Double Crystal Ball Fit to \\pp\\ (left) and pPb (right) for $y$ $\\in$ [%.2f-%.2f].",yLow,yHigh);
    }

    cout << "\%" << endl;
    cout << "\\begin{figure}[hbtp]" << endl;
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PP_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLowPP, yHighPP) << endl;
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLow, yHigh) << endl;
    cout << caption.Data() << endl;
    cout << Form("	\\label{fig:Nominalpt%.0fto%.0fy%.3ito%.3i}}", ptLow, ptHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
    cout << "\\end{figure}" << endl;
    if (i<0) {
      cout << "\%INTEGRATED BIN" << endl;
    }
  }
  cout << "\%" << endl;
}
