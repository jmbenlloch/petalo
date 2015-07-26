void zReconsRatio_analysis2(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",-50,50);
	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TFile *file = TFile::Open("lxsc2_z2.root");
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("pes");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();
	c1->Print("deliverables2/zRatio_2_z2.pdf");
	file->Close();
	
	TFile *file = TFile::Open("lxsc2_z3.root");
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("pes");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();
	c1->Print("deliverables2/zRatio_2_z3.pdf");

	TFile *file = TFile::Open("lxsc2_z4.root");
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("pes");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();
	c1->Print("deliverables2/zRatio_2_z4.pdf");

	TFile *file = TFile::Open("lxsc2_z5.root");
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("pes");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();
	c1->Print("deliverables2/zRatio_2_z5.pdf");

	c1->Close();
}
