void positionPlots(){
	gStyle->SetOptStat(kFALSE);
	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_xy->SetTitle("Event position (XY plane)");
	petAnalysis_xy->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_xy->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_xy->Draw("colz");

	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_xz->SetTitle("Event position (XZ plane)");
	petAnalysis_xz->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_xz->GetYaxis()->SetTitle("z (mm)");
	petAnalysis_xz->Draw("colz");

	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_yz->SetTitle("Event position (YZ plane)");
	petAnalysis_yz->GetXaxis()->SetTitle("y (mm)");
	petAnalysis_yz->GetYaxis()->SetTitle("z (mm)");
	petAnalysis_yz->Draw("colz");

	TCanvas* c4 = new TCanvas("c4", "First canvas", 800,600);
	petAnalysis_z->SetTitle("Proportion of events in Z position");
	petAnalysis_z->Rebin(100);
	petAnalysis_z->ComputeIntegral();
	petAnalysis_z->SetAxisRange(0.,50.,"X");
	Double_t *integral = petAnalysis_z->GetIntegral();
	petAnalysis_z->SetContent(integral);
	petAnalysis_z->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_z->GetYaxis()->SetTitle("Proportion");
	petAnalysis_z->Draw();

}
