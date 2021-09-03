#include "MakeLatexTable.C"
#include "MakeLatexTable_pt0to6to30.C"
#include "MakeLatexTable_y193to000to193.C"
#include "MakeLatexTable_HFNtracks.C"

void MakeAllLatexTables() {

  MakeLatexTable(1);
  MakeLatexTable(2);
  MakeLatexTable(3);
  MakeLatexTable_pt0to6to30(1);
  MakeLatexTable_pt0to6to30(2);
  MakeLatexTable_pt0to6to30(3);
  MakeLatexTable_y193to000to193(1);
  MakeLatexTable_y193to000to193(2);
  MakeLatexTable_y193to000to193(3);
  MakeLatexTable_HFNtracks(1,0);
  MakeLatexTable_HFNtracks(2,0);
  MakeLatexTable_HFNtracks(3,0);
  MakeLatexTable_HFNtracks(1,1);
  MakeLatexTable_HFNtracks(2,1);
  MakeLatexTable_HFNtracks(3,1);

}
