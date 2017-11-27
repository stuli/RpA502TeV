void mergeRuns()
{
	TString filename1 = "../../../RD2013_pa_1st_run_merged.root";
	TString filename2 = "../../../RD2013_pa_2nd_run_merged.root";
	TString outfilename = "../../../RD2013_pa_bothruns_merged.root";
	
	TFile* file1 = new TFile(filename1,"READ");
	TFile* file2 = new TFile(filename2,"READ");
	TFile* outfile = new TFile(outfilename,"RECREATE");
	
	TTree* tree1 = (TTree*)file1->Get("myTree");
	TTree* tree2 = (TTree*)file2->Get("myTree");
	
	//TCut kinCut = "Reco_QQ_mupl_4mom.Pt() > 4 && Reco_QQ_mumi_4mom.Pt() > 4 && TMath::Abs(Reco_QQ_mupl_4mom.Eta()) < 2.4 && TMath::Abs(Reco_QQ_mumi_4mom.Eta()) < 2.4";
	//TTree* cutTree1 = tree1->CopyTree(kinCut);
	//TTree* cutTree2 = tree2->CopyTree(kinCut);
	
	TList* treeList = new TList;
	treeList->Add(tree1);
	treeList->Add(tree2);
	TTree* mergedTree = TTree::MergeTrees(treeList);
	mergedTree->SetName("mergedTree");
	
	mergedTree->Write();
}