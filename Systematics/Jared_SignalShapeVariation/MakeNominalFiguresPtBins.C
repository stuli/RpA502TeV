#include "../../HeaderFiles/cutsAndBin.h"

void MakeNominalFiguresPtBins(int whichUpsilon=1, int whichYRange=0) {

  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;

  float ptLow, ptHigh, yLow, yHigh;
  float yLowPP = 0.0;
  float yHighPP = 1.93;
  TString caption;

  if (whichYRange==0) {
    yLow = -1.93;
    yHigh = 1.93;
  }
  else if (whichYRange==1) {
    yLow = -1.93;
    yHigh = 0.0;
  }
  else if (whichYRange==2) {
    yLow = 0.0;
    yHigh = 1.93;
  }

  for (int i = 0; i<numptbins; i++) {

    if (i<0) {
      ptLow = 0;
      ptHigh = 30;
    }
    else {
      ptLow = ptbins[i];
      ptHigh = ptbins[i+1];
    }
    caption = Form("	\\caption{Double Crystal Ball Fit to \\pp\\ (left) and pPb (right) for $\\pt~[\\GeVc]$ $\\in$ [%.1f,%.1f], y $\\in$ [%.2f,%.2f].}",ptLow,ptHigh,yLow,yHigh);

    cout << "\\begin{figure}[hbtp]" << endl;
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PP_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLowPP, yHighPP) << endl;
    cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", ptLow, ptHigh, yLow, yHigh) << endl;
    cout << caption.Data() << endl;
    //cout << Form("	\\label{fig:NomFitpt%.0fto%.0fy%.3ito%.3i}", ptLow, ptHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
    cout << "\\end{figure}" << endl;
    cout << "\%" << endl;
  }
}
