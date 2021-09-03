#include "GetFitParams.C"

void GetFitParamsInTable(int whichUpsilon=1, int collId=kPPDATA) {

  //TNtuple* params = new TNtuple("params","params","m0:n:alpha:sigma0:f:x:mu:sigma:lambda",1);

  float paramsArray[13] = {0};
  float * paramsptr;
  paramsptr = &paramsArray[0];

  float errorsArray[13] = {0};
  float * errorsptr;
  errorsptr = &errorsArray[0];

  int ycheck1, ycheck2;
  //choose a set of bins
  if (collId==kPADATA) {
    ycheck1=9;
    ycheck2=14;
    if (whichUpsilon==1) {
      float ptbins[14] = {0,2,4,6,9,12,30,0,4,9,30,0,6,30};
      float ybins[18] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93,-1.93,-0.8,0.0,0.8,1.93,-1.93,0.0,1.93};
    }
    else if (whichUpsilon==2) {
      float ptbins[4] = {0,4,9,30};
      float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
    }
    else if (whichUpsilon==3) {
      float ptbins[3] = {0,6,30};
      float ybins[3] = {-1.93,0.0,1.93};
    }
  }
  else if (collId==kPPDATA) {
    ycheck1=4;
    ycheck2=7;
    if (whichUpsilon==1) {
      float ptbins[14] = {0,2,4,6,9,12,30,0,4,9,30,0,6,30};
      float ybins[10] = {0.0,0.4,0.8,1.2,1.93,0.0,0.8,1.93,0.0,1.93};
    }
    else if (whichUpsilon==2) {
      float ptbins[4] = {0,4,9,30};
      float ybins[3] = {0.0,0.8,1.93};
    }
    else if (whichUpsilon==3) {
      float ptbins[3] = {0,6,30};
      float ybins[2] = {0.0,1.93};
    }
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;

  bool DoPt = kFALSE;
  bool DoY = kTRUE;

  TString paramName[13] = {"$m_0$","$\\sigma_0$","$f$","$n$","$\\alpha$","$x$","$n_{1S}$","$n_{2S}$","$n_{3S}$","$\\mu$","$\\sigma$","$\\lambda$","$n_{Bkg}$"};

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

  TString Options[4] = {"signal","fixed","yields","background"};

for (int l=3; l<4; l++) {
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

if (DoPt) {
  cout << "\\begin{table}[hbtp]" << endl;
  cout << "\\label{" << strId.Data() << "pt" << whichParams.Data() <<   "params}" << endl;
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
  for (int j=0; j<numptbins; j++) {
    if (j==6) {
      cout << "\\hline" << endl;
      continue;
    }
    if (j==10) {
      cout << "\\hline" << endl;
      continue;
    }
    GetFitParams(collId,ptbins[j],ptbins[j+1],0.0,1.93,paramsptr,errorsptr);
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$",ptbins[j],ptbins[j+1]);
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

if (DoY) {
  cout << "\\begin{table}[hbtp]" << endl;
  cout << "\\label{" << strId.Data() << "yHighPt" << whichParams.Data() <<   "params}" << endl;
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
  for (int j=0; j<numybins; j++) {
    if (j==ycheck1) {
      cout << "\\hline" << endl;
      continue;
    }
    if (j==ycheck2) {
      cout << "\\hline" << endl;
      continue;
    }
    GetFitParams(collId,6,30,ybins[j],ybins[j+1],paramsptr,errorsptr);
    cout << Form("$%.2f<y<%.2f$",ybins[j],ybins[j+1]);
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
}
}
