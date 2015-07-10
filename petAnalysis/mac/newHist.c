void newHist(){

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_EPosX->Draw("colz");
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_EPosZ->Draw("colz");
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_xPosZ->Draw("colz");
	TCanvas* c4 = new TCanvas("c4", "First canvas", 800,600);
	petAnalysis_yPosZ->Draw("colz");

	TCanvas* c5 = new TCanvas("c5", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",0,10000);
	petAnalysis_Energy->Fit("gauF","","e",-25,25);

	petAnalysis_zRatio->Draw("colz");

}
