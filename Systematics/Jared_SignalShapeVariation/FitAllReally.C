#include "FitAll.C"
#include "FitAllHFNtracks.C"

void FitAllReally() {

  //Fit the integrated bin
  /*FitDataWithNominalSeeds((int)kPADATA, 0, 30, -1.93, 1.93);

  //Fit the regular pt and y bins
  FitAll((int)kPADATA, 0, 30, -1.93, 1.93, kTRUE, kTRUE);
*/
  //Fit the pt bins in backward and forward rapidity
 /* FitAll((int)kPADATA, 0, 30, -1.93, 0.0, kTRUE, kFALSE);
  FitAll((int)kPADATA, 0, 30, 0.0, 1.93, kTRUE, kFALSE);
 */
  //Fit the y bins in low and high pt
/*  FitAll((int)kPADATA, 0, 6, -1.93, 1.93, kFALSE, kTRUE, 2);
  FitAll((int)kPADATA, 6, 30, -1.93, 1.93, kFALSE, kTRUE, 2);

  //Fit the integrated bin in -2.87<y<1.93 in pPb.
  FitDataWithNominalSeeds((int)kPADATA, 0, 30, -2.87, 1.93);
*/
  //Fit the pt bins in rapidity range [-2.87,1.93]
 /* FitAll((int)kPADATA, 0, 30, -2.87, 1.93, kTRUE, kFALSE);

  //Fit the bin -2.87<y<-1.93 in pPb.
  FitDataWithNominalSeeds((int)kPADATA, 0, 30, -2.87, -1.93);
*/
  //Fit the regular pt and y bins in pp.
  /*FitAll((int)kPPDATA, 0, 30, 0.0, 1.93, kTRUE, kTRUE);//The 3S y bin is the integrated bin.

  //Fit all the y bins in low and high pt in pp for 1S and 2S.
  FitAll((int)kPPDATA, 0, 6, 0.0, 1.93, kFALSE, kTRUE, 2);
  FitAll((int)kPPDATA, 6, 30, 0.0, 1.93, kFALSE, kTRUE, 2);
*/
  //I fit all the HF and Ntracks bins in pPb.
  FitAllHFNtracks();

  gROOT->ProcessLine(".q");
}
