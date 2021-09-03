#include "GetFitParamsHFNtracks.C"

void GetFitParamsInTableHFNtracks(int whichUpsilon=1) {

  //TNtuple* params = new TNtuple("params","params","m0:n:alpha:sigma0:f:x:mu:sigma:lambda",1);

  float paramsArray[13] = {0};
  float * paramsptr;
  paramsptr = &paramsArray[0];

  float errorsArray[13] = {0};
  float * errorsptr;
  errorsptr = &errorsArray[0];

  //choose a set of bins
  if (whichUpsilon==1) {
    float hfbins[7] = {0,12,19,27,120,12,120};
    int ntbins[7] = {0,40,62,88,400,40,400};
  }
  else if (whichUpsilon==2) {
    float hfbins[5] = {0,12,19,27,120};
    int ntbins[5] = {0,40,62,88,400};
  }
  else if (whichUpsilon==3) {
    float hfbins[3] = {0,12,120};
    int ntbins[3] = {0,40,400};
  }

  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntbins = sizeof(ntbins)/sizeof(float)-1;

  bool DoHF = kFALSE;
  bool DoNtracks = kTRUE;

  TString paramName[13] = {"$m_0$","$\\sigma_0$","$f$","$n$","$\\alpha$","$x$","$n_{1S}$","$n_{2S}$","$n_{3S}$","$\\mu$","$\\sigma$","$\\lambda$","$n_{Bkg}$"};

  TString strId = "PA";

  TString Options[4] = {"signal","fixed","yields","background"};

for (int l=0; l<4; l++) {
  TString whichParams=Options[l];//"signal","fixed","yields","background"
  //cout << whichParams.Data() << endl;
  int istart, iend;
  if (whichParams=="signal") {
    istart=0;
    iend=3;
  }
  else if (whichParams=="fixed") {
    istart=3;
    iend=6;
  }
  else if (whichParams=="yields") {
    istart=6;
    iend=9;
  }
  else if (whichParams=="background") {
    istart=9;
    iend=13;
  }

if (DoHF) {
  cout << "\\begin{table}[hbtp]" << endl;
  cout << "\\label{" << strId.Data() << "HFForwy" << whichParams.Data() <<   "params}" << endl;
  cout << "\\centering" << endl;
  cout << "\\begin{tabular}{|c|";
  for (int k=istart; k<iend; k++) {
    cout << "c|";
  }
  cout << "}" << endl;
  cout << "\\hline" << endl;
  cout << "Bin";
  for (int k=istart; k<iend; k++) {
    cout << " & " << paramName[k].Data();
  }
  cout << "\\\\" << endl;
  cout << "\\hline" << endl;
  for (int j=0; j<numhfbins; j++) {
    if (j==4) continue;
    GetFitParamsHFNtracks(kPADATA,0,30,0.0,1.93,hfbins[j],hfbins[j+1],0,400,paramsptr,errorsptr);
    cout << Form("$%.0f<E_T^{HF}<%.0f \\GeV$",hfbins[j],hfbins[j+1]);
    for (int i=istart; i<iend; i++) {
      cout << " & " << Form("%.2f$\\pm$%.2f",paramsArray[i],errorsArray[i]);
    }
    cout << "\\\\" << endl;
  }
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\caption{}" << endl;
  cout << "\\end{table}" << endl;
}
  //cout << "\\hline" << endl;
if (DoNtracks) {
  cout << "\\begin{table}[hbtp]" << endl;
  cout << "\\label{" << strId.Data() << "NtracksForwy" << whichParams.Data() <<   "params}" << endl;
  cout << "\\centering" << endl;
  cout << "\\begin{tabular}{|c|";
  for (int k=istart; k<iend; k++) {
    cout << "c|";
  }
  cout << "}" << endl;
  cout << "\\hline" << endl;
  cout << "Bin";
  for (int k=istart; k<iend; k++) {
    cout << " & " << paramName[k].Data();
  }
  cout << "\\\\" << endl;
  cout << "\\hline" << endl;
  for (int j=0; j<numntbins; j++) {
    if (j==4) continue;
    GetFitParamsHFNtracks(kPADATA,0,30,0.0,1.93,0,120,ntbins[j],ntbins[j+1],paramsptr,errorsptr);
    cout << Form("$%i<=\\mathrm{Ntracks}<%i$",ntbins[j],ntbins[j+1]);
    for (int i=istart; i<iend; i++) {
      cout << " & " << Form("%.2f$\\pm$%.2f",paramsArray[i],errorsArray[i]);
    }
    cout << "\\\\" << endl;
  }
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\caption{}" << endl;
  cout << "\\end{table}" << endl;
}
cout << endl;
}
}
