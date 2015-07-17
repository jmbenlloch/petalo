void lego(){

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_phot150->SetTitle("Charge distribution (Photoelectric event)");
	petAnalysis_phot150->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_phot150->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_phot150->GetXaxis()->SetTitleOffset(2);
	petAnalysis_phot150->GetYaxis()->SetTitleOffset(2);
	petAnalysis_phot150->Draw("LEGO");
	
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_noise->SetTitle("Charge distribution");
	petAnalysis_noise->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_noise->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_noise->GetXaxis()->SetTitleOffset(2);
	petAnalysis_noise->GetYaxis()->SetTitleOffset(2);
	petAnalysis_noise->Draw("LEGO");
	
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_compt1->SetTitle("Charge distribution (Compton)");
	petAnalysis_compt1->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_compt1->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_compt1->GetXaxis()->SetTitleOffset(2);
	petAnalysis_compt1->GetYaxis()->SetTitleOffset(2);
	petAnalysis_compt1->Draw("LEGO");
	

}
