 //gStyle->SetPalette(1);

void noNormHist5(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

	petAnalysis_xNoNorm_80->SetAxisRange(-50.,50.,"X");
	petAnalysis_yNoNorm_80->SetAxisRange(-50.,50.,"X");
	petAnalysis_zNoNorm_80->SetAxisRange(-50.,50.,"X");

//	petAnalysis_xNoNorm_80->Fit("gauF","","e",-30,30);
//	petAnalysis_yNoNorm_80->Fit("gauF","","e",-30,30);
//	petAnalysis_zNoNorm_80->Fit("gauF","","e",-30,30);

	c1->Divide(3,0);
	c1->cd(1);
	petAnalysis_xNoNorm_80->Draw();
	c1->cd(2);
	petAnalysis_yNoNorm_80->Draw();
	c1->cd(3);
	petAnalysis_zNoNorm_80->Draw();
}
