void energy(){
//	gStyle->SetOptStat(kFALSE);

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_Energy->SetTitle("Total Charge");
	petAnalysis_Energy->GetXaxis()->SetTitle("pe");
	petAnalysis_Energy->GetYaxis()->SetTitle("counts");
	petAnalysis_Energy->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Energy->SetAxisRange(0., 7000.0,"X");

	TF1* gauF2 = new TF1("gauF2","gaus",0,10000);
	petAnalysis_Energy->Fit("gauF2","","e",5000,6000);

	petAnalysis_Energy->Draw();

	//petAnalysis_Energy->GetFunction("gauF")->SetBit(TF1::kNotDraw);
}
