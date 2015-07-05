 //gStyle->SetPalette(1);

void nearest(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);


	petAnalysis_x2Near_65->SetAxisRange(-50.,50.,"X");
	petAnalysis_y2Near_65->SetAxisRange(-50.,50.,"X");
	petAnalysis_z2Near_65->SetAxisRange(-50.,50.,"X");

	petAnalysis_x2Near_65->Fit("gauF","","e",-20,20);
	petAnalysis_y2Near_65->Fit("gauF","","e",-20,20);
	petAnalysis_z2Near_65->Fit("gauF","","e",-20,20);

	c1->Divide(3,0);
	c1->cd(1);
	petAnalysis_x2Near_65->Draw();
	c1->cd(2);
	petAnalysis_y2Near_65->Draw();
	c1->cd(3);
	petAnalysis_z2Near_65->Draw();
}
