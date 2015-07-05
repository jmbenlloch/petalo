 //gStyle->SetPalette(1);

void normHist(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

	petAnalysis_xNorm_70->SetAxisRange(-50.,50.,"X");
	petAnalysis_yNorm_70->SetAxisRange(-50.,50.,"X");
	petAnalysis_zNorm_70->SetAxisRange(-50.,50.,"X");

//	petAnalysis_xNorm_70->Fit("gauF","","e",-30,30);
//	petAnalysis_yNorm_70->Fit("gauF","","e",-30,30);
//	petAnalysis_zNorm_70->Fit("gauF","","e",-30,30);

	c1->Divide(3,0);
	c1->cd(1);
	petAnalysis_xNorm_70->Draw();
	c1->cd(2);
	petAnalysis_yNorm_70->Draw();
	c1->cd(3);
	petAnalysis_zNorm_70->Draw();
}
