 //gStyle->SetPalette(1);

void photEnergy(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	petAnalysis_PhotEnergy->Rebin(300);
	petAnalysis_XeGammaEnergy->Rebin(20);

	petAnalysis_PhotEnergy->SetAxisRange(0,1.2,"X");
	petAnalysis_XeGammaEnergy->SetAxisRange(0,0.04,"X");

	c1->Divide(2,0);
	c1->cd(1);
	petAnalysis_PhotEnergy->Draw();
	c1->cd(2);
	petAnalysis_XeGammaEnergy->Draw();
}
