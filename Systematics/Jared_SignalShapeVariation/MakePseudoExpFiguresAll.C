#include "MakePseudoExpFiguresPtBins.C"
#include "MakePseudoExpFiguresYBins.C"
#include "MakePseudoExpFiguresHFBins.C"
#include "MakePseudoExpFiguresNtracksBins.C"

void MakePseudoExpFiguresAll() {

  cout << "\\subsection{Results of Pseudoexperiments in $P_T$ Bins}" << endl << endl;
  cout << "%PT BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresPtBins(1);
  MakePseudoExpFiguresPtBins(2);
  MakePseudoExpFiguresPtBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Results of Pseudoexperiments in Rapidity Bins}" << endl << endl;
  cout << "%RAPIDITY BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresYBins(1);
  MakePseudoExpFiguresYBins(2);
  MakePseudoExpFiguresYBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Results of Pseudoexperiments in Differential Bins}" << endl << endl;
  cout << "%DIFFERENTIAL PT0-6-30 BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresYBins(1,1);
  MakePseudoExpFiguresYBins(2,1);
  MakePseudoExpFiguresYBins(3,1);
  MakePseudoExpFiguresYBins(1,2);
  MakePseudoExpFiguresYBins(2,2);
  MakePseudoExpFiguresYBins(3,2);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "%DIFFERENTIAL Y-193-0-193 BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresPtBins(1,1);
  MakePseudoExpFiguresPtBins(2,1);
  MakePseudoExpFiguresPtBins(3,1);
  MakePseudoExpFiguresPtBins(1,2);
  MakePseudoExpFiguresPtBins(2,2);
  MakePseudoExpFiguresPtBins(3,2);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Results of Pseudoexperiments in HF Bins}" << endl << endl;
  cout << "%HF BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresHFBins(1);
  MakePseudoExpFiguresHFBins(2);
  MakePseudoExpFiguresHFBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;

  cout << "\\subsection{Results of Pseudoexperiments in Ntracks Bins}" << endl << endl;
  cout << "%NTRACKS BINS PSEUDOEXPERIMENTS" << endl;
  MakePseudoExpFiguresNtracksBins(1);
  MakePseudoExpFiguresNtracksBins(2);
  MakePseudoExpFiguresNtracksBins(3);
  cout << endl;
  cout << "\\clearpage" << endl << endl;
}
