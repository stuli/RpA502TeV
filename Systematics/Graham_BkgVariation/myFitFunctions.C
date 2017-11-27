void myFitFunctions()
{
	cout << "This does nothing.\n";
}

double fExpErf(double* xVar, double* par)
{
	double x = xVar[0];
	double norm = par[0];
	double mu = par[1];
	double sigma = par[2];
	double lambda = par[3];
	return norm*((TMath::Erf((x-mu)/(TMath::Sqrt(2)*sigma))+1)/2*TMath::Exp(-x/lambda));
}

double fSumErfExp(double* xVar, double* par)
{
	double x = xVar[0];
	double norm = par[0];
	double ratio = par[1];
	double muErf = par[2];
	double sigma = par[3];
	double muExp = par[4];
	double lambdaExp = par[5];
	double lambdaDecay = par[6];
	
	double erfpart = (TMath::Erf((x-muErf)/(TMath::Sqrt(2)*sigma))+1);
	double exppart = 0;
	if (x > muExp)
		exppart = 1-TMath::Exp(-(x-muExp)/lambdaExp);
	double decay = TMath::Exp(-x/lambdaDecay);
	
	return norm*(ratio*erfpart + exppart)*decay;
}

//Parameters from Heather
double fSumErfExpTotal(double* xVar, double* par)
{
	double x = xVar[0];
	double norm1 = par[1];
	double norm2 = par[2];
	double norm3 = par[3];
	double norm4 = par[4];
	double norm5 = par[5];
	double lambdaDecay = par[0];
	
	double decay = TMath::Exp(-x/lambdaDecay);
	
	//Individual params, fits
	double ratio1 = 43.566/346.165; //ratio
	double muErf1 = 8.65893; //erf mu
	double sigma1 = 0.31689; //erf sigma
	double muExp1 = 8.107; //exp mu
	double lambdaExp1 = 2.71228; // exp lambda
	
	double erfpart1 = (TMath::Erf((x-muErf1)/(TMath::Sqrt(2)*sigma1))+1);
	double exppart1 = 0;
	if (x > muExp1)
		exppart1 = 1-TMath::Exp(-(x-muExp1)/lambdaExp1);
	double fpart1 = norm1*(ratio1*erfpart1 + exppart1);
	
	double ratio2 = 238.397/563.232; //ratio
	double muErf2 = 9.05477; //erf mu
	double sigma2 = 0.530565; //erf sigma
	double muExp2 = 9.34286; //exp mu
	double lambdaExp2 = 2.91684; // exp lambda
	
	double erfpart2 = (TMath::Erf((x-muErf2)/(TMath::Sqrt(2)*sigma2))+1);
	double exppart2 = 0;
	if (x > muExp2)
		exppart2 = 1-TMath::Exp(-(x-muExp2)/lambdaExp2);
	double fpart2 = norm2*(ratio2*erfpart2 + exppart2);
	
	double ratio3 = 297.451/915.932; //ratio
	double muErf3 = 10.2115; //erf mu
	double sigma3 = 1.07635; //erf sigma
	double muExp3 = 7.99382; //exp mu
	double lambdaExp3 = 6.44623; // exp lambda
	
	double erfpart3 = (TMath::Erf((x-muErf3)/(TMath::Sqrt(2)*sigma3))+1);
	double exppart3 = 0;
	if (x > muExp3)
		exppart3 = 1-TMath::Exp(-(x-muExp3)/lambdaExp3);
	double fpart3 = norm3*(ratio3*erfpart3 + exppart3);
	
	double ratio4 = 356.276/1032.99; //ratio
	double muErf4 = 10.8803; //erf mu
	double sigma4 = 1.32576; //erf sigma
	double muExp4 = 7.37455; //exp mu
	double lambdaExp4 = 9.92836; // exp lambda
	
	double erfpart4 = (TMath::Erf((x-muErf4)/(TMath::Sqrt(2)*sigma4))+1);
	double exppart4 = 0;
	if (x > muExp4)
		exppart4 = 1-TMath::Exp(-(x-muExp4)/lambdaExp4);
	double fpart4 = norm4*(ratio4*erfpart4 + exppart4);
	
	double ratio5 = 299.226/1035.74; //ratio
	double muErf5 = 11.4665; //erf mu
	double sigma5 = 1.37887; //erf sigma
	double muExp5 = 6.94769; //exp mu
	double lambdaExp5 = 9.84612; // exp lambda
	
	double erfpart5 = (TMath::Erf((x-muErf5)/(TMath::Sqrt(2)*sigma5))+1);
	double exppart5 = 0;
	if (x > muExp5)
		exppart5 = 1-TMath::Exp(-(x-muExp5)/lambdaExp5);
	double fpart5 = norm5*(ratio5*erfpart5 + exppart5);
	
	
	return (fpart1+fpart2+fpart3+fpart4+fpart5)*decay;
	
}

double fSumErfExpTotal2(double* xVar, double* par)
{
	/*double a1 = par[0];
	double a2 = par[1]
	double a3 = par[2]
	double a4 = par[3]
	double a5 = par[4]*/
	double result;
	for (int i = 0; i<5; i++)
		result += par[i]*fSumErfExp(xVar,par+5+7*i);
	return result;
}

double fSumErfExpTotal3(double* xVar, double* par)
{
	double x = xVar[0];
	
	double norm1 = par[0];
	double ratio1 = par[1]; //ratio
	double muErf1 = par[2]; //erf mu
	double sigma1 = par[3]; //erf sigma
	double muExp1 = par[4]; //exp mu
	double lambdaExp1 = par[5]; // exp lambda
	double lambdaDecay1 = par[6];
	
	double erfpart1 = (TMath::Erf((x-muErf1)/(TMath::Sqrt(2)*sigma1))+1);
	double exppart1 = 0;
	if (x > muExp1)
		exppart1 = 1-TMath::Exp(-(x-muExp1)/lambdaExp1);
	double decay1 = TMath::Exp(-x/lambdaDecay1);
	double fpart1 = norm1*(ratio1*erfpart1 + exppart1)*decay1;
	
	double norm2 = par[7];
	double ratio2 = par[8]; //ratio
	double muErf2 = par[9]; //erf mu
	double sigma2 = par[10]; //erf sigma
	double muExp2 = par[11]; //exp mu
	double lambdaExp2 = par[12]; // exp lambda
	double lambdaDecay2 = par[13];
	
	double erfpart2 = (TMath::Erf((x-muErf2)/(TMath::Sqrt(2)*sigma2))+1);
	double exppart2 = 0;
	if (x > muExp2)
		exppart2 = 1-TMath::Exp(-(x-muExp2)/lambdaExp2);
	double decay2 = TMath::Exp(-x/lambdaDecay2);
	double fpart2 = norm2*(ratio2*erfpart2 + exppart2)*decay2;
	
	double norm3 = par[14];
	double ratio3 = par[15]; //ratio
	double muErf3 = par[16]; //erf mu
	double sigma3 = par[17]; //erf sigma
	double muExp3 = par[18]; //exp mu
	double lambdaExp3 = par[19]; // exp lambda
	double lambdaDecay3 = par[20];
	
	double erfpart3 = (TMath::Erf((x-muErf3)/(TMath::Sqrt(2)*sigma3))+1);
	double exppart3 = 0;
	if (x > muExp3)
		exppart3 = 1-TMath::Exp(-(x-muExp3)/lambdaExp3);
	double decay3 = TMath::Exp(-x/lambdaDecay3);
	double fpart3 = norm3*(ratio3*erfpart3 + exppart3)*decay3;
	
	double norm4 = par[21];
	double ratio4 = par[22]; //ratio
	double muErf4 = par[23]; //erf mu
	double sigma4 = par[24]; //erf sigma
	double muExp4 = par[25]; //exp mu
	double lambdaExp4 = par[26]; // exp lambda
	double lambdaDecay4 = par[27];
	
	double erfpart4 = (TMath::Erf((x-muErf4)/(TMath::Sqrt(2)*sigma4))+1);
	double exppart4 = 0;
	if (x > muExp4)
		exppart4 = 1-TMath::Exp(-(x-muExp4)/lambdaExp4);
	double decay4 = TMath::Exp(-x/lambdaDecay4);
	double fpart4 = norm4*(ratio4*erfpart4 + exppart4)*decay4;
	
	double norm5 = par[28];
	double ratio5 = par[29]; //ratio
	double muErf5 = par[30]; //erf mu
	double sigma5 = par[31]; //erf sigma
	double muExp5 = par[32]; //exp mu
	double lambdaExp5 = par[33]; // exp lambda
	double lambdaDecay5 = par[34];
	
	double erfpart5 = (TMath::Erf((x-muErf5)/(TMath::Sqrt(2)*sigma5))+1);
	double exppart5 = 0;
	if (x > muExp5)
		exppart5 = 1-TMath::Exp(-(x-muExp5)/lambdaExp5);
	double decay5 = TMath::Exp(-x/lambdaDecay5);
	double fpart5 = norm5*(ratio5*erfpart5 + exppart5)*decay5;
	
	return fpart1+fpart2+fpart3+fpart4+fpart5;
}

double fSumErfExpTotalFixed(double* xVar, double* par)
{
	double x = xVar[0];
	
	double norm[5];
	double ratio[5];
	double muErf[5];
	double sigma[5];
	double muExp[5];
	double lambdaExp[5];
	double lambdaDecay[5];
	
	for (int i = 0; i < 5; i++)
		norm[i] = par[i];
	
	//Fix all parameters except norms
	ratio[0] = 0.348998;
	muErf[0] = 8.67629;
	sigma[0] = 0.264901;
	muExp[0] = 8;
	lambdaExp[0] = 1.17301;
	lambdaDecay[0] = 3.5058;
	ratio[1] = 0.849509;
	muErf[1] = 9.28584;
	sigma[1] = 0.574242;
	muExp[1] = 7.96041;
	lambdaExp[1] = 1.30271;
	lambdaDecay[1] = 3.72635;
	ratio[2] = 2.61185;
	muErf[2] = 9.7654;
	sigma[2] = 0.864455;
	muExp[2] = 7.7708;
	lambdaExp[2] = 0.557547;
	lambdaDecay[2] = 3.99346;
	ratio[3] = 1.8412;
	muErf[3] = 10;
	sigma[3] = 1.19453;
	muExp[3] = 7.38372;
	lambdaExp[3] = 0.870146;
	lambdaDecay[3] = 4.45817;
	ratio[4] = 0.676893;
	muErf[4] = 9.83576;
	sigma[4] = 1.51859;
	muExp[4] = 6.81872;
	lambdaExp[4] = 1.52168;
	lambdaDecay[4] = 5.80538;
	
	//Calculate
	double result = 0;
	
	for (int i = 0; i < 5; i++)
	{
		double erfpart = (TMath::Erf((x-muErf[i])/(TMath::Sqrt(2)*sigma[i]))+1);
		double exppart = 0;
		if (x > muExp[i])
			exppart = 1-TMath::Exp(-(x-muExp[i])/lambdaExp[i]);
		double decay = TMath::Exp(-x/lambdaDecay[i]);
		
		result += norm[i]*(ratio[i]*erfpart + exppart)*decay;
	}
	
	return result;
}

double fRapidity(double* xVar, double* par)
{
	double x = xVar[0];
	double A = par[0];
	double pgmu = par[1];
	double pgsig = par[2];
	double ratio = par[3];
	double sgmu = par[4];
	double sgsig = par[5];
	double sgerf = par[6];
	
	double part1 = A*(TMath::Exp(-0.5*TMath::Power(x-pgmu,2)/TMath::Power(pgsig,2)) + ratio*TMath::Exp(-0.5*TMath::Power(x-sgmu,2)/TMath::Power(sgsig,2))*(1+TMath::Erf(sgerf*(x-sgmu)/TMath::Sqrt(2))));
	double part2 = A*(TMath::Exp(-0.5*TMath::Power(-x-pgmu,2)/TMath::Power(pgsig,2)) + ratio*TMath::Exp(-0.5*TMath::Power(-x-sgmu,2)/TMath::Power(sgsig,2))*(1+TMath::Erf(sgerf*(-x-sgmu)/TMath::Sqrt(2))));
	
	return part1+part2;
}

