#include "MakeParamsHistos_pPbVSPbp.C"
#include "pPbVSPbp_PlotOneYield.C"

void pPbVSPbp_PlotOneYieldAll() {

  int collId=kPADATA;

  /*MakeParamsHistos_pPbVSPbp(collId,1,1);
  MakeParamsHistos_pPbVSPbp(collId,1,2);
  MakeParamsHistos_pPbVSPbp(collId,1,3);
  MakeParamsHistos_pPbVSPbp(collId,2,1);
  MakeParamsHistos_pPbVSPbp(collId,2,2);
  MakeParamsHistos_pPbVSPbp(collId,2,3);
  MakeParamsHistos_pPbVSPbp(collId,3,1);
  MakeParamsHistos_pPbVSPbp(collId,3,2);
  MakeParamsHistos_pPbVSPbp(collId,3,3);
*/
  for (int whichUpsilon=1; whichUpsilon<4; whichUpsilon++) {
    for (int whichParam=9; whichParam<13; whichParam++) {
      pPbVSPbp_PlotOneYield(collId,whichUpsilon,whichParam);
    }
  }

}
