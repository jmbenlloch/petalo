void correctedBias(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

//	c1->SetFillColor(0);
//	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

/*	petAnalysis_xCoronna1->SetTitle("xRecons - xTrue");
	petAnalysis_yCoronna1->SetTitle("yRecons - yTrue");

	petAnalysis_xCoronna1->SetAxisRange(-20.,20.,"X");
	petAnalysis_xCoronna1->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xCoronna1->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xCoronna1->GetYaxis()->SetTitle("counts");

	petAnalysis_yCoronna1->SetAxisRange(-20.,20.,"X");
	petAnalysis_yCoronna1->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yCoronna1->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yCoronna1->GetYaxis()->SetTitle("counts");

	petAnalysis_xCoronna1->Fit("gauF","","e",-15,15);
	petAnalysis_yCoronna1->Fit("gauF","","e",-15,15);
*/
	//Best
	petAnalysis_xBest->SetTitle("xRecons - xTrue (-10 to -10)");
	petAnalysis_yBest->SetTitle("yRecons - yTrue (-10 to -10)");

	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");

	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");

	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Fit("gauF","","e",-15,15);

/*	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_xCoronna1->Draw();
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_yCoronna1->Draw();*/
	TCanvas* c4 = new TCanvas("c4", "First canvas", 800,600);
	petAnalysis_xBest->Draw();
	TCanvas* c5 = new TCanvas("c5", "First canvas", 800,600);
	petAnalysis_yBest->Draw();

	petAnalysis_xCorrected->SetTitle("xReconsCorrected - xTrue (-10 to 10)");
	petAnalysis_yCorrected->SetTitle("yReconsCorrected - yTrue (-10 to 10)");

	petAnalysis_xCorrected->SetAxisRange(-20.,20.,"X");
	petAnalysis_xCorrected->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xCorrected->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xCorrected->GetYaxis()->SetTitle("pes");

	petAnalysis_yCorrected->SetAxisRange(-20.,20.,"X");
	petAnalysis_yCorrected->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yCorrected->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yCorrected->GetYaxis()->SetTitle("pes");

	petAnalysis_xCorrected->Fit("gauF","","e",-15,15);
	petAnalysis_yCorrected->Fit("gauF","","e",-15,15);
	TCanvas* c7 = new TCanvas("c7", "First canvas", 800,600);
	petAnalysis_xCorrected->Draw();
	TCanvas* c8 = new TCanvas("c8", "First canvas", 800,600);
	petAnalysis_yCorrected->Draw();

}