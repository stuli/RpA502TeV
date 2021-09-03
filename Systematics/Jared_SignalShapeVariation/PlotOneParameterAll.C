#include "MakeParamsHistos_compare.C"
#include "PlotOneParameter.C"

void PlotOneParameterAll() {

  int collId=kPADATA;

  MakeParamsHistos_compare(collId,1,2);
  MakeParamsHistos_compare(collId,1,3);
  MakeParamsHistos_compare(collId,1,4);
  MakeParamsHistos_compare(collId,2,2);
  MakeParamsHistos_compare(collId,2,3);
  MakeParamsHistos_compare(collId,2,4);
  MakeParamsHistos_compare(collId,3,2);
  MakeParamsHistos_compare(collId,3,3);
  MakeParamsHistos_compare(collId,3,4);

  for (int whichUpsilon=1; whichUpsilon<4; whichUpsilon++) {
    for (int whichParam=0; whichParam<13; whichParam++) {
      PlotOneParameter(collId,whichUpsilon,whichParam);
    }
  }

}
