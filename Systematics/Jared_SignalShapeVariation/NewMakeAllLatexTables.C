#include "NewMakeLatexTable.C"
#include "NewMakeLatexTable_pt0to6to30.C"
#include "NewMakeLatexTable_y193to000to193.C"
#include "NewMakeLatexTable_HFNtracks.C"

void NewMakeAllLatexTables(TString dir="ErrorEstimates", TString whichSyst="") {

//dir = "ErrorEstimates", "ErrorEstimatesAsOfAug13", etc.
//whichSyst = "", "Combined", "_ParamFixingOnly"

  NewMakeLatexTable(1,dir,whichSyst);
  NewMakeLatexTable(2,dir,whichSyst);
  NewMakeLatexTable(3,dir,whichSyst);
  NewMakeLatexTable_pt0to6to30(1,dir,whichSyst);
  NewMakeLatexTable_pt0to6to30(2,dir,whichSyst);
  NewMakeLatexTable_pt0to6to30(3,dir,whichSyst);
  //NewMakeLatexTable_y193to000to193(1,dir,whichSyst);
  //NewMakeLatexTable_y193to000to193(2,dir,whichSyst);
  //NewMakeLatexTable_y193to000to193(3,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(1,0,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(2,0,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(3,0,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(1,1,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(2,1,dir,whichSyst);
  NewMakeLatexTable_HFNtracks(3,1,dir,whichSyst);

}
