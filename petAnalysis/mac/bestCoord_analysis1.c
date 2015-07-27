void bestCoord_analysis1(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",-50,50);
	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	//////////////
	// LXSC2_36 //
	//////////////
	TFile *file = TFile::Open("lxsc2_36.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_2_36.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_2_36.pdf");
	file->Close();

	//////////////
	// LXSC2_49 //
	//////////////
	TFile *file = TFile::Open("lxsc2_49.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_2_49.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_2_49.pdf");
	file->Close();
	
	//////////////
	// LXSC2_64 //
	//////////////
	TFile *file = TFile::Open("lxsc2_64.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_2_64.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_2_64.pdf");
	file->Close();

	//////////////
	// LXSC4_36 //
	//////////////
	TFile *file = TFile::Open("lxsc4_36.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_4_36.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_4_36.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_4_36.pdf");
	file->Close();

	//////////////
	// LXSC4_49 //
	//////////////
	TFile *file = TFile::Open("lxsc4_49.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_4_49.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_4_49.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_4_49.pdf");
	file->Close();
	
	//////////////
	// LXSC4_64 //
	//////////////
	TFile *file = TFile::Open("lxsc4_64.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_4_64.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_4_64.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_4_64.pdf");
	file->Close();

	//////////////
	// LXSC6_36 //
	//////////////
	TFile *file = TFile::Open("lxsc6_36.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_6_36.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_6_36.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_6_36.pdf");
	file->Close();
	
	//////////////
	// LXSC6_49 //
	//////////////
	TFile *file = TFile::Open("lxsc6_49.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_6_49.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_6_49.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_6_49.pdf");
	file->Close();

	//////////////
	// LXSC6_64 //
	//////////////
	TFile *file = TFile::Open("lxsc6_64.root");
	petAnalysis_xBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_xBest->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xBest->GetYaxis()->SetTitle("pes");
	petAnalysis_xBest->Fit("gauF","","e",-15,15);
	petAnalysis_xBest->Draw();
	c1->Print("deliverables1/xBest_6_64.pdf");
	petAnalysis_yBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_yBest->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_yBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_yBest->GetYaxis()->SetTitle("pes");
	petAnalysis_yBest->Fit("gauF","","e",-15,15);
	petAnalysis_yBest->Draw();
	c1->Print("deliverables1/yBest_6_64.pdf");
	petAnalysis_zBest->SetAxisRange(-20.,20.,"X");
	petAnalysis_zBest->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zBest->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zBest->GetYaxis()->SetTitle("pes");
	petAnalysis_zBest->Fit("gauF","","e",-15,15);
	petAnalysis_zBest->Draw();
	c1->Print("deliverables1/zBest_6_64.pdf");
	file->Close();

	c1->Close();
}
