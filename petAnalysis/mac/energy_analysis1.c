void energy_analysis1(){

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",0,15000);

/*	TFile *file = TFile::Open("lxsc2_36.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",8000,10000);
	petAnalysis_EnergyPhot->Draw();
	file->Close();*/

	TFile *file = TFile::Open("lxsc2_49.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",8000,10000);
	c1->Update();
	TPaveStats *st = (TPaveStats*)petAnalysis_EnergyPhot->FindObject("stats");
	st->SetX1NDC(0.15);
	st->SetX2NDC(0.45);
	st->SetY1NDC(0.45);
	st->SetY2NDC(0.85);
	petAnalysis_EnergyPhot->Draw();
	c1->Print("deliverables1/energy_2_49.pdf");
	file->Close();

	TFile *file = TFile::Open("lxsc2_64.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",8000,10000);
	c1->Update();
	TPaveStats *st = (TPaveStats*)petAnalysis_EnergyPhot->FindObject("stats");
	st->SetX1NDC(0.15);
	st->SetX2NDC(0.45);
	st->SetY1NDC(0.45);
	st->SetY2NDC(0.85);
	petAnalysis_EnergyPhot->Draw();
	c1->Print("deliverables1/energy_2_64.pdf");
	file->Close();
	
/*	TFile *file = TFile::Open("lxsc4_36.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",9000,11000);
	petAnalysis_EnergyPhot->Draw();
	file->Close();*/

/*	TFile *file = TFile::Open("lxsc4_49.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",9000,11000);
	petAnalysis_EnergyPhot->Draw();
	file->Close();*/

	TFile *file = TFile::Open("lxsc4_64.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",9000,11000);
	c1->Update();
	TPaveStats *st = (TPaveStats*)petAnalysis_EnergyPhot->FindObject("stats");
	st->SetX1NDC(0.15);
	st->SetX2NDC(0.45);
	st->SetY1NDC(0.45);
	st->SetY2NDC(0.85);
	petAnalysis_EnergyPhot->Draw();
	c1->Print("deliverables1/energy_4_64.pdf");
	file->Close();
	
/*	TFile *file = TFile::Open("lxsc6_36.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",10000,11000);
	petAnalysis_EnergyPhot->Draw();
	file->Close();*/
	
	TFile *file = TFile::Open("lxsc6_49.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",9000,11000);
	c1->Update();
	TPaveStats *st = (TPaveStats*)petAnalysis_EnergyPhot->FindObject("stats");
	st->SetX1NDC(0.15);
	st->SetX2NDC(0.45);
	st->SetY1NDC(0.45);
	st->SetY2NDC(0.85);
	petAnalysis_EnergyPhot->Draw();
	c1->Print("deliverables1/energy_6_49.pdf");
	file->Close();

	TFile *file = TFile::Open("lxsc6_64.root");
	petAnalysis_EnergyPhot->GetXaxis()->SetTitle("pe");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitle("counts");
	petAnalysis_EnergyPhot->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_EnergyPhot->SetAxisRange(0., 15000.0,"X");
	petAnalysis_EnergyPhot->Fit("gauF","","e",10000,12000);
	c1->Update();
	TPaveStats *st = (TPaveStats*)petAnalysis_EnergyPhot->FindObject("stats");
	st->SetX1NDC(0.15);
	st->SetX2NDC(0.45);
	st->SetY1NDC(0.45);
	st->SetY2NDC(0.85);
	petAnalysis_EnergyPhot->Draw();
	c1->Print("deliverables1/energy_6_64.pdf");
	file->Close();

}
