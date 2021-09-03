#include "MakeNominalFiguresPtBins.C"
#include "MakeNominalFiguresYBins.C"
#include "MakeNominalFiguresHFBins.C"
#include "MakeNominalFiguresNtracksBins.C"

void MakeNominalFiguresAll() {

  cout << "\\subsection{Nominal Fits to Data in $P_T$ Bins}" << endl << endl;
  cout << "%PT BINS" << endl;
  MakeNominalFiguresPtBins(1);
  MakeNominalFiguresPtBins(2);
  MakeNominalFiguresPtBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Nominal Fits to Data in Rapidity Bins}" << endl << endl;
  cout << "%RAPIDITY BINS" << endl;
  MakeNominalFiguresYBins(1);
  MakeNominalFiguresYBins(2);
  MakeNominalFiguresYBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Nominal Fits to Data in Differential Bins}" << endl << endl;
  cout << "%DIFFERENTIAL PT0-6-30 BINS" << endl;
  MakeNominalFiguresYBins(1,1);
  MakeNominalFiguresYBins(2,1);
  MakeNominalFiguresYBins(3,1);
  MakeNominalFiguresYBins(1,2);
  MakeNominalFiguresYBins(2,2);
  MakeNominalFiguresYBins(3,2);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "%DIFFERENTIAL Y-193-0-193 BINS" << endl;
  MakeNominalFiguresPtBins(1,1);
  MakeNominalFiguresPtBins(2,1);
  MakeNominalFiguresPtBins(3,1);
  MakeNominalFiguresPtBins(1,2);
  MakeNominalFiguresPtBins(2,2);
  MakeNominalFiguresPtBins(3,2);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Nominal Fits to Data in HF Bins}" << endl << endl;
  cout << "%HF BINS" << endl;
  MakeNominalFiguresHFBins(1);
  MakeNominalFiguresHFBins(2);
  MakeNominalFiguresHFBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Nominal Fits to Data in Ntracks Bins}" << endl << endl;
  cout << "%NTRACKS BINS" << endl;
  MakeNominalFiguresNtracksBins(1);
  MakeNominalFiguresNtracksBins(2);
  MakeNominalFiguresNtracksBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;
}
