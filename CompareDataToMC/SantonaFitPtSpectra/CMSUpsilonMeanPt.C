
void CMSUpsilonMeanPt() {
  // The first CMS Paper has various rapidity ranges.
  // The closest one to the STAR acceptance is |y|<1.
  // The data for this rapidity range is the one plotted here.
  // In that paper, the data are given as sigma (cross section in
  // a given bin), so I divide by the bin width to get dsigma/dpT vs. pT.
  // The 2015 CMS Paper, the data are given as dsigma/dpT,
  //  but is in the range |y|<0.6, so assuming a flat rapidity, there is a correction of a factor 0.6.
  // To correct for the flat rapidity assumption, we note that the data in the 2015 CMS paper from |y|<1.2 compared to the
  // data from |y|<0.6 is lower by about 10%.  So we can apply a 10% correction to lower the data.
  // The overlap bin is 8-11. From pT=10 GeV/c, I use the 2015 data.  Since the 2011 paper has the 8-11 bin cross section as
  // 0.44 nb, I divide by the bin witdh of 3 GeV to get dsigma/dp, but I also apply an ad-hoc 10% correction because if one were
  // to integrate between 8-10 instead, there is also an ~11% correction, estimated using the fit function.
  const int nBins1s = 14;
  double xbins1s[] =     {0,           2.,      5.,     8.,     10.,       12., 14., 16., 18., 20., 22., 24., 26., 28., 30.};
  double dsigmadpt1s[] = {0.7/2., 1.54/3., 1.02/3., 0.44/3./0.89, .060936/0.6/1.1, .033828/0.6/1.1, .019670/0.6/1.1, .011504/0.6/1.1, .006914/0.6/1.1, .004235/0.6/1.1, .002721/0.6/1.1, .001848/0.6/1.1, .001117/0.6/1.1, .000845/0.6/1.1};

  double dsigmadpt1serr[] = {0.05*dsigmadpt1s[0],0.04*dsigmadpt1s[1], 0.05*dsigmadpt1s[2], 0.06*dsigmadpt1s[3], dsigmadpt1s[4]*0.005, dsigmadpt1s[5]*0.005, dsigmadpt1s[6]*0.006,  dsigmadpt1s[7]*0.008, dsigmadpt1s[8]*0.009, dsigmadpt1s[9]*0.011, dsigmadpt1s[10]*0.014, dsigmadpt1s[11]*0.017, dsigmadpt1s[12]*0.022, dsigmadpt1s[13]*0.025};
  
  TH1D* CMSUps1s = new TH1D("CMSUps1s","CMS #varUpsilon data;p_{T} (GeV/c);d#sigma/dp_{T} (nb/(GeV/c))",nBins1s,xbins1s);
  CMSUps1s->SetMarkerStyle(20);
  CMSUps1s->SetMarkerColor(4);
  CMSUps1s->GetYaxis()->SetTitleOffset(1.3);
  
  double dsigmady1s(0), dsigmady1serr(0);
  for (size_t i1s=0;i1s<nBins1s;++i1s) {
    CMSUps1s->SetBinContent(i1s+1,dsigmadpt1s[i1s]);
    CMSUps1s->SetBinError(i1s+1,dsigmadpt1serr[i1s]);
    //cout << "Pt bin " << i1s+1 << " sigma " << CMSUps1s->GetBinContent(i1s+1) << endl;
    dsigmady1s+=CMSUps1s->GetBinContent(i1s+1);
    dsigmady1serr+=pow(CMSUps1s->GetBinError(i1s+1),2);
  }
  TCanvas* CMSUpsData = new TCanvas("CMSUpsData","CMS Upsilon Data",500,500);
  CMSUpsData->SetLeftMargin(0.11);
  gPad->SetLogx(0);
  gPad->SetLogy(1);
  gPad->SetTicks(1,1);

  CMSUps1s->Draw("e");
  CMSUps1s->SetMaximum(4);
  CMSUps1s->SetMinimum(9e-4);

    // Tsallis function
  TF1* upsPtFit = new TF1("upsPtFit","[0]*([2]-1)*([2]-2)/([2]*[1])/([2]*[1]+9.46*([2]-2))*x/(pow(1.+((sqrt(x*x+9.46*9.46)-9.46)/([2]*[1])),[2]))",0,30); // Tsallis function
  upsPtFit->SetParameters(100,3.0,4.);
  upsPtFit->SetParNames("A","T","n");
  upsPtFit->FixParameter(2,6.0);
  upsPtFit->SetLineColor(kBlue);
  
  CMSUps1s->Fit("upsPtFit","","same",1,30); // fit for 1S
  cout << "X^2/dof = " << upsPtFit->GetChisquare()/upsPtFit->GetNDF() << endl;
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  // For calculating mean pT, need to multiply the fit function by pT
  // integrate it, then divide by the fit function integral, i.e.
  // <pT> = integral(pt*pdf) / integra(pdf)
  TF1* upsPtFitPt = new TF1("upsPtFitPt","[0]*([2]-1)*([2]-2)/([2]*[1])/([2]*[1]+9.46*([2]-2))*x*x/(pow(1.+((sqrt(x*x+9.46*9.46)-9.46)/([2]*[1])),[2]))",0,30);
  upsPtFitPt->SetParameter(0,upsPtFit->GetParameter(0));
  upsPtFitPt->SetParameter(1,upsPtFit->GetParameter(1));
  upsPtFitPt->SetParameter(2,upsPtFit->GetParameter(2));
  cout << "Fit function: <pT> = " << upsPtFitPt->Integral(0,30)/upsPtFit->Integral(0,30) << endl;
  cout << "Histogram   : <pT> = " << CMSUps1s->GetMean() << endl;
  return;
}
