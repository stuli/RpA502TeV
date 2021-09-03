#include "../../HeaderFiles/cutsAndBin.h"

void PrintTheStuff(int collId=kPADATA, TString param="alpha") {

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

cout << "\%" << Form("\%s",param.Data()) << endl;
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_1sbins}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_2sbins}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_3sbins}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s versus \\pt in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_1sbins2}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_2sbins2}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_3sbins2}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s versus rapidity in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_1sbins_pt0to6to30}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_2sbins_pt0to6to30}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_3sbins_pt0to6to30}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s in low and high \\pt versus rapidity in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_1sbins_y193to0to193}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_2sbins_y193to0to193}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/%sfitted_%s_3sbins_y193to0to193}.png}",strId.Data(),param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s in backward and forward rapidity versus \\pt in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
if (collId==kPADATA){
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_1sbins_HF_y1.20-1.93}.png}",param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_2sbins_HF_y0.80-1.93}.png}",param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_3sbins_HF_y0.00-1.93}.png}",param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s in $1.2<|y|<1.93$ (left), $0.8<|y|<1.93$ (center), and $0.0<|y|<1.93$ (right) versus $E_T^{HF|\\eta|>4}$ in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
cout << Form("\\begin{figure}[htbp]") << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_1sbins_Ntracks_y1.20-1.93}.png}",param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_2sbins_Ntracks_y0.80-1.93}.png}",param.Data()) << endl;
cout << Form("	\\includegraphics[width=0.33\\textwidth]{{figures/NominalFits/PAfitted_%s_3sbins_Ntracks_y0.00-1.93}.png}",param.Data()) << endl;
cout << Form("	\\caption{Signal parameter %s in $1.2<|y|<1.93$ (left), $0.8<|y|<1.93$ (center), and $0.0<|y|<1.93$ (right) versus Ntracks in \\pPb\\.}",param.Data()) << endl;
cout << Form("\\end{figure}") << endl;
cout << "\%" << endl;
}

cout << Form("\\afterpage{\\clearpage}") << endl;

}

void MakeParameterLatexPlots(int SigOrBkg=0) {

if (SigOrBkg==0){
  TString params[9] = {"alpha","f1s","mass","n1s","nSig1s","nSig2s","nSig3s","sigma1s","x1s"};

  cout << "\\subsection{Signal parameters in pPb fits}" << endl << endl;
  for (int i=0; i<9; i++){
    PrintTheStuff(kPADATA,params[i].Data());
  }
  cout << "\\subsection{Signal parameters in pp fits}" << endl << endl;
  for (int i=0; i<9; i++){
    PrintTheStuff(kPPDATA,params[i].Data());
  }
}
else{
  TString bkgparams[4] = {"errmu","errsigma","lambda","nBkg"};
  cout << "\\subsection{Background parameters in pPb fits}" << endl << endl;
  for (int i=0; i<4; i++){
    PrintTheStuff(kPADATA,bkgparams[i].Data());
  }
  cout << "\\subsection{Background parameters in pp fits}" << endl << endl;
  for (int i=0; i<4; i++){
    PrintTheStuff(kPPDATA,bkgparams[i].Data());
  }
}
}
