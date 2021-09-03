#include "../../HeaderFiles/cutsAndBin.h"

void MakePseudoExpFiguresPtBins(int whichUpsilon=1, int whichYRange=0) {

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
    caption = Form("	\\caption{Percent differences in $\\Upsilon$(%iS) yields between the nominal and alternative signal PDF fits to pseudodata generated in the bin: $\\pt$ $\\in$ [%.1f,%.1f], y $\\in$ [%.2f,%.2f].}",whichUpsilon,ptLow,ptHigh,yLow,yHigh);

    cout << "\\begin{figure}[hbtp]" << endl;
    cout << Form("	\\includegraphics[width=1.0\\textwidth]{{figures/systematics_SignalPDF/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0}.pdf}", whichUpsilon, ptLow, ptHigh, yLow, yHigh) << endl;
    cout << caption.Data() << endl;
    //cout << Form("	\\label{fig:PseudoExp%iSpt%.0fto%.0fy%.3ito%.3i}", whichUpsilon, ptLow, ptHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
    cout << "\\end{figure}" << endl;
    cout << "\%" << endl;
  }
}
