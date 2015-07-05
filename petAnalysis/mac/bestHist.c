 //gStyle->SetPalette(1);

void bestHist(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

	petAnalysis_xbestTrue_75->SetAxisRange(-50.,50.,"X");
	petAnalysis_ybestTrue_80->SetAxisRange(-50.,50.,"X");
	petAnalysis_zbestTrue_85->SetAxisRange(-50.,50.,"X");

	petAnalysis_xbestTrue_75->Fit("gauF","","e",-15,15);
	petAnalysis_ybestTrue_80->Fit("gauF","","e",-15,15);
	petAnalysis_zbestTrue_85->Fit("gauF","","e",-15,15);

	c1->Divide(3,0);
	c1->cd(1);
	petAnalysis_xbestTrue_75->Draw();
	c1->cd(2);
	petAnalysis_ybestTrue_80->Draw();
	c1->cd(3);
	petAnalysis_zbestTrue_85->Draw();
}
