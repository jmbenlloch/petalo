 //gStyle->SetPalette(1);

void bestHist(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("counts");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("counts");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("counts");

	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Fit("gauF","","e",-15,15);

	petAnalysis_xBest->Draw();
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_yBest->Draw();
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_zBest->Draw();

/*	TCanvas* c10 = new TCanvas("c10", "First canvas", 800,600);
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("counts");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();*/
}
