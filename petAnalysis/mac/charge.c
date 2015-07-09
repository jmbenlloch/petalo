void charge(){
	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_maxCharge->SetTitle("Max charge read by a SiPM");
	petAnalysis_maxCharge->GetXaxis()->SetTitle("pe");
	petAnalysis_maxCharge->GetYaxis()->SetTitle("counts");
	petAnalysis_maxCharge->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_maxCharge->Draw();

	TCanvas* c2	= new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_avgCharge->SetTitle("Average charge");
	petAnalysis_avgCharge->GetXaxis()->SetTitle("pe");
	petAnalysis_avgCharge->GetYaxis()->SetTitle("counts");
	petAnalysis_avgCharge->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_avgCharge->Draw();

	TCanvas* c3	= new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_avgChargeNoMax->SetTitle("Average charge excluding max value");
	petAnalysis_avgChargeNoMax->GetXaxis()->SetTitle("pe");
	petAnalysis_avgChargeNoMax->GetYaxis()->SetTitle("counts");
	petAnalysis_avgChargeNoMax->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_avgChargeNoMax->Draw();
}
