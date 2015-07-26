void bestCoord_analysis2(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",-50,50);
	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	//////////////
	// LXSC2_Z3 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z2.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables2/xBest_2_z2.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables2/yBest_2_z2.pdf");
	file->Close();

	//////////////
	// LXSC2_Z3 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z3.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables2/xBest_2_z3.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables2/yBest_2_z3.pdf");
	file->Close();

	//////////////
	// LXSC2_Z4 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z4.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables2/xBest_2_z4.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables2/yBest_2_z4.pdf");
	file->Close();
	
	//////////////
	// LXSC2_Z5 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z5.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables2/xBest_2_z5.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables2/yBest_2_z5.pdf");
	file->Close();

	c1->Close();
}
