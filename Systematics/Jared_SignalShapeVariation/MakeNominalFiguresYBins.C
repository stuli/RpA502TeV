#include "../../HeaderFiles/cutsAndBin.h"

void MakeNominalFiguresYBins(int whichUpsilon=1, int whichPtRange=0) {

  //choose a set of bins
  if (whichUpsilon==1) {
    float ybins[11] = {-2.87,-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[3] = {-1.93,0.0,1.93};
  }

  const int numybins = sizeof(ybins)/sizeof(float)-1;

  float ptLow, ptHigh, yLow, yHigh, yLowPP, yHighPP;
  TString caption;
  int istart = 0;

  if (whichPtRange==0) {
    ptLow = 0;
    ptHigh = 30;
  }
  else if (whichPtRange==1) {
    istart = 2;
    ptLow = 0;
    ptHigh = 6;
  }
  else if (whichPtRange==2) {
    istart = 2;
    ptLow = 6;
    ptHigh = 30;
  }

  for (int i = istart; i<numybins; i++) {

    yLow = ybins[i];
    yHigh = ybins[i+1];
    if (yLow<0) {
      yLowPP = TMath::Abs(yHigh);
      yHighPP = TMath::Abs(yLow);
    }
    else {
      yLowPP = yLow;
      yHighPP = yHigh;
    }
    caption = Form("	\\caption{Double Crystal Ball Fit to \\pp\\ (left) and pPb (right) for $\\pt~[\\GeVc]$ $\\in$ [%.1f,%.1f], y $\\in$ [%.2f,%.2f].}",ptLow,ptHigh,yLow,yHigh);

    cout << "\\begin{figure}[hbtp]" << endl;
    if (yLow>-2.8) {
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PP_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLowPP, yHighPP) << endl;
    }
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLow, yHigh) << endl;
    cout << caption.Data() << endl;
    //cout << Form("	\\label{fig:NomFitpt%.0fto%.0fy%.3ito%.3i}", ptLow, ptHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
    cout << "\\end{figure}" << endl;
    cout << "\%" << endl;
  }
}
