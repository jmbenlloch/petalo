void energy(){
//	gStyle->SetOptStat(kFALSE);

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_Energy->SetTitle("Total Charge");
	petAnalysis_Energy->GetXaxis()->SetTitle("pe");
	petAnalysis_Energy->GetYaxis()->SetTitle("counts");
	petAnalysis_Energy->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Energy->SetAxisRange(0., 20000.0,"X");

	TF1* gauF2 = new TF1("gauF2","gaus",0,10000);
	petAnalysis_Energy->Fit("gauF2","","e",12000,14000);

	petAnalysis_Energy->Draw();


	TCanvas* c2	= new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_Energy->SetTitle("Total Charge Phot");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 20000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF2","","e",12000,14000);
	petAnalysis_EnergyPhot->Draw();

	TCanvas* c3	= new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_Energy->SetTitle("Total Charge Compt");
	petAnalysis_EnergyCompt->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyCompt->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyCompt->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyCompt->SetAxisRange(0., 20000.0,"X");
	petAnalysis_EnergyCompt->Fit("gauF2","","e",12000,14000);
	petAnalysis_EnergyCompt->Draw();


	//petAnalysis_Energy->GetFunction("gauF")->SetBit(TF1::kNotDraw);
}
