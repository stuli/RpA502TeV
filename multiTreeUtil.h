#ifndef MULTITREETUIL_H
#define MULTITREETUIL_H

// =========================================================================
// Multi Tree Utility
// Author: Yen-Jie Lee (MIT)
// Updated by Yongsun Kim (MIT)
//       A tool for multiTree handling
// =========================================================================
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TTree.h>
#include <TCut.h>
#include <TFile.h>
//#include <vector.h>
#include <TH1.h>


class multiTreeUtil
{
   public:
      multiTreeUtil() {
         nTrees=0;
	 for (int i=0;i<100;i++)
	 {
	    color[i]=i;
	 }
	 color[0]=kOrange-3;
	 color[1]=3;
	 color[2]=4;
	 color[3]=5;
	 color[4]=6;
	 color[5]=9;
	 color[6]=1;
         color[7]=kAzure+3;
         color[8]=kOrange+3;
	 color[9]=kSpring-7;
	 color[10]=kPink+4;

	 /*
	   handsomeTH1(tmp1[1],2);
	   handsomeTH1(tmp1[2],4);
	   handsomeTH1(tmp1[3],kGreen+2);
	   handsomeTH1(tmp1[4],6);
	   handsomeTH1(tmp1[5],7);
	   handsomeTH1(tmp1[6],kGreen+4);
	   handsomeTH1(tmp1[7],kBlue-2);
	   handsomeTH1(tmp1[8],kAzure+3);
	   handsomeTH1(tmp1[9],kOrange+3);
	   handsomeTH1(tmp1[10],kSpring-7);
	   handsomeTH1(tmp1[11],kPink+4);
	 */
	 
      }
      ~multiTreeUtil() {
	 
      }
      void addFile(TString filename, TString treeName, TCut cut, Float_t scaleFactor=1);
      void Draw(TH1D* h,TString expression,TCut cut = "" , TString indWeight="", TString opt="");
      void Draw2(TH1D* h,TString expression,TCut cut = "", TString indWeight="", TString opt="",bool hist=true);
      void Draw2D(TH2D* h,TString expression,TCut cut = "", TString indWeight="", TString opt="",bool hist=true);
      void Draw3(TH1D* h,TString expression,TCut cut = "", TString indWeight="", TString opt="",bool hist=true);
      void Print(Option_t* option = "");
      void AddFriend(TString friendname="");
      void SetAlias(TString nameAl="", TString varAl="");
      TTree *getTree(int i) { return trees_[i]; };
	
 private:
      int nTrees;
      vector<TFile*>  tfiles_;                // TFile 
      vector<TTree*>  trees_;                 // Trees
      vector<TString>   tfileNames_;            // File names
      vector<TCut>    tcuts_;
      vector<Float_t> scaleFactors_;
      int color[100];

      void cleanHist(TH1* h);

};

// =========================================================================
// Add file
// =========================================================================
void multiTreeUtil::addFile(TString filename, TString treeName, TCut cut, Float_t scaleFactor)
{
   TFile *f = new TFile(filename.Data());
   TTree *t = (TTree*) f->FindObjectAny(treeName.Data());
   
   if (t!=0) {
      nTrees++;
      tfileNames_.push_back(filename);
      tfiles_.push_back(f);
      trees_.push_back(t);
      tcuts_.push_back(cut);
      scaleFactors_.push_back(scaleFactor);
   }
   else
      cout << " there is no such tree in this file :" << treeName << endl;
}

// =========================================================================
// Draw
// =========================================================================
void multiTreeUtil::Draw(TH1D *h, TString expression, TCut cut,TString indWeight, TString opt)
{
  TH1::SetDefaultSumw2();

  h->Reset();
  //  h->Sumw2();
   for (int i=0;i<(int)trees_.size();i++)
   {
     TH1D *htmp = (TH1D*)h->Clone();
     htmp->SetName(Form("my_htmp_%d",i));
     htmp->Reset();
     //   htmp->Sumw2();
     TCut fCut = cut&&tcuts_[i];
     if ( indWeight=="")
       trees_[i]->Draw(Form("%s>>my_htmp_%d",expression.Data(),i),  Form("(%s)", fCut.GetTitle() ));
     else
       trees_[i]->Draw(Form("%s>>my_htmp_%d",expression.Data(),i),  Form("(%s) * (%s)", fCut.GetTitle(), indWeight.Data() ) );
     htmp->Scale(scaleFactors_[i]);
      h->Add(htmp);
      delete htmp;
   }
   
   h->Draw(opt.Data());
}

// =========================================================================
// Draw components
// =========================================================================
void multiTreeUtil::Draw2(TH1D *h, TString expression, TCut cut, TString indWeight, TString opt,bool hist)
{
  TH1::SetDefaultSumw2();

  //TH1D *hComponent[100];
   
   h->Reset();
   //   h->Sumw2();
   
   for (int i=0;i<(int)trees_.size();i++)
     {
       TString hName = Form("%s_Draw_%d",h->GetName(),i);
       TH1D *htmp= (TH1D*)h->Clone();
       htmp->SetName(hName.Data());
       htmp->Reset();
       //   htmp->Sumw2();
       
       TCut fCut = cut&&tcuts_[i];
       
       if ( fCut == "" ) fCut = "1==1";
       
       if ( indWeight=="" ) {
	 trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("(%s)", fCut.GetTitle() ));
       }
       else {
	 trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("( %s ) * (%s)", fCut.GetTitle(),  indWeight.Data()) );
       }
       cout << " Running draw2" << endl;
       cout << " Var             :  " << expression.Data() << endl;
       cout << "(cut) * (weight) :  " << endl << Form("( %s )     *     ( %s )", fCut.GetTitle(),  indWeight.Data()) << endl;
       cout << " and scale factor = " << double(scaleFactors_[i]) << endl;
       h->Add(htmp, double(scaleFactors_[i]) );
       delete htmp;
       //  hComponent[i]=htmp;
     }
   //   h->Sumw2();
   //  h->DrawCopy(opt.Data());
   
   /*
   for (int i=0;i<(int)trees_.size();i++)
     {
       hComponent[i]->SetLineColor(color[i]);
       hComponent[i]->SetMarkerColor(color[i]);
       hComponent[i]->SetFillColor(color[i]);
       if ( hist) opt ="hist";
       //   hComponent[i]->DrawCopy(Form("%s same",opt.Data()));
       delete hComponent[i];
     }
   */

}


void multiTreeUtil::Draw2D(TH2D *h, TString expression, TCut cut, TString indWeight, TString opt,bool hist)
{
  TH1::SetDefaultSumw2();

  TH2D *hComponent[20];

  h->Reset();


  for (int i=0;i<(int)trees_.size();i++)
    {
      TString hName = Form("%s_Draw_%d",h->GetName(),i);
      TH2D *htmp= (TH2D*)h->Clone();
      //      htmp->Sumw2();
      htmp->SetName(hName.Data());

      TCut fCut = cut&&tcuts_[i];

      if ( indWeight=="" ) {
        trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("(%s)", fCut.GetTitle() ));
      }
      else {
        trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("( %s ) * %s", fCut.GetTitle(),  indWeight.Data()) );
      }
      htmp->Scale(scaleFactors_[i]);
      h->Add(htmp);
      hComponent[i]=htmp;

    }
  //  h->Sumw2();
  //  h->DrawCopy(opt.Data());

  for (int i=0;i<(int)trees_.size();i++)
    {
      hComponent[i]->SetLineColor(color[i]);
      hComponent[i]->SetMarkerColor(color[i]);
      hComponent[i]->SetFillColor(color[i]);
      if ( hist) opt ="hist";
      //      hComponent[i]->DrawCopy(Form("%s same",opt.Data()));
      delete hComponent[i];
    }
}




void multiTreeUtil::Draw3(TH1D *h, TString expression, TCut cut, TString indWeight, TString opt,bool hist)
{
  TH1::SetDefaultSumw2();

   TH1D *hComponent[100];

   h->Reset();


   for (int i=0;i<(int)trees_.size();i++)
      {
	 TString hName = Form("%s_Draw_%d",h->GetName(),i);
	 TH1D *htmp= (TH1D*)h->Clone();
	 htmp->SetName(hName.Data());

	 TCut fCut = cut&&tcuts_[i];

	 if ( indWeight=="" )
	    trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("%s", fCut.GetTitle()));
	 else
	    trees_[i]->Draw(Form("%s>>%s",expression.Data(),hName.Data()), Form("(%s) * %s", fCut.GetTitle(),  indWeight.Data()) );
	 //	 htmp->Sumw2();
	 htmp->Scale(scaleFactors_[i]);
	 h->Add(htmp);
	 hComponent[i]=htmp;
	 hComponent[i]->SetLineColor(color[i]);
         hComponent[i]->SetMarkerColor(color[i]);
         hComponent[i]->SetFillColor(color[i]);
      }
   //   h->Sumw2();
   h->DrawCopy(opt.Data());

   
   for (int i=0;i<(int)trees_.size();i++)
      {
	 if ( i<(int)trees_.size() -1 ) 
	    hComponent[i+1]->Add( hComponent[i]);
      }
   
   for (int i=0;i<(int)trees_.size();i++)
      {
         if ( hist) opt ="hist";
         hComponent[trees_.size()-i-1]->DrawCopy(Form("%s same",opt.Data()));
      }
}


// =========================================================================
// Print variables from the first tree
// =========================================================================
void multiTreeUtil::AddFriend(TString friendname)
{
   for (int i=0;i<(int)trees_.size();i++)
      {
	 trees_[i]->AddFriend(friendname.Data(), tfileNames_[i]);
      }
}

void multiTreeUtil::SetAlias(TString nameAl, TString varAl) 
{
  for (int i=0;i<(int)trees_.size();i++)
    {
      trees_[i]->SetAlias( nameAl.Data(), varAl.Data());
    }
}



void multiTreeUtil::Print(Option_t* option)
{
   trees_[0]->Print(option);   
}

//--------------------------------PRIVATE-----------------------------------

// =========================================================================
// Clear
// =========================================================================
void multiTreeUtil::cleanHist(TH1* h)
{
   for (int i=0;i<=h->GetNbinsX()+1;i++) {
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
   }
}

#endif
