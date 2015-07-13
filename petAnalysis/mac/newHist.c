void newHist(){

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_EPosX->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_EPosX->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EPosX->GetYaxis()->SetTitle("\deltaE (pe)");
	petAnalysis_EPosX->Draw("colz");
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_EPosZ->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_EPosZ->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EPosZ->GetYaxis()->SetTitle("\deltaE (pe)");
	petAnalysis_EPosZ->Draw("colz");

	TCanvas* c5 = new TCanvas("c5", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",0,10000);
	petAnalysis_Energy->Fit("gauF","","e",-25,25);

	petAnalysis_zRatio->Draw("colz");

//	gStyle->SetOptStat(kFALSE);

	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_xPosZ->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_xPosZ->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xPosZ->GetYaxis()->SetTitle("\deltax (mm)");
	petAnalysis_xPosZ->Rebin2D(2,2);
	petAnalysis_xPosZ->Draw("colz");

	TCanvas* c4 = new TCanvas("c4", "First canvas", 800,600);
	petAnalysis_yPosZ->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_yPosZ->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yPosZ->GetYaxis()->SetTitle("\deltay (mm)");
	petAnalysis_yPosZ->Rebin2D(2,2);
	petAnalysis_yPosZ->Draw("colz");

}
