void param_analysis2(){
	gStyle->SetOptStat(kFALSE);
//	gStyle->SetOptFit(1);
//	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	TF1* gauF = new TF1("gauF","gaus",-50,50);
	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	//////////////
	// LXSC2_Z2 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z2.root");
	petAnalysis_Param_SiPMMC_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane0->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane2->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p2_2_z2.pdf");
	petAnalysis_Param_SiPMMC_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p2_2_z2.pdf");

	petAnalysis_Param_SiPMMC_C1_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane0->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane2->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p2_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p2_2_z2.pdf");

	petAnalysis_Param_SiPMMC_C2_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane0->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p0_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane2->SetAxisRange(-10.,10.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p2_2_z2.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p2_2_z2.pdf");
	file->Close();

	//////////////
	// LXSC2_Z3 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z3.root");
	petAnalysis_Param_SiPMMC_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane0->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane2->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p2_2_z3.pdf");
	petAnalysis_Param_SiPMMC_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p2_2_z3.pdf");

	petAnalysis_Param_SiPMMC_C1_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane0->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane2->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p2_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p2_2_z3.pdf");

	petAnalysis_Param_SiPMMC_C2_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane0->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p0_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane2->SetAxisRange(-15.,15.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p2_2_z3.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p2_2_z3.pdf");
	file->Close();

	//////////////
	// LXSC2_Z4 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z4.root");
	petAnalysis_Param_SiPMMC_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane0->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane2->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p2_2_z4.pdf");
	petAnalysis_Param_SiPMMC_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p2_2_z4.pdf");

	petAnalysis_Param_SiPMMC_C1_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane0->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane2->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p2_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p2_2_z4.pdf");

	petAnalysis_Param_SiPMMC_C2_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane0->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p0_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane2->SetAxisRange(-20.,20.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p2_2_z4.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p2_2_z4.pdf");
	file->Close();
	
	//////////////
	// LXSC2_Z5 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z5.root");
	petAnalysis_Param_SiPMMC_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane0->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_Plane2->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_p2_2_z5.pdf");
	petAnalysis_Param_SiPMMC_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_p2_2_z5.pdf");

	petAnalysis_Param_SiPMMC_C1_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane0->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C1_Plane2->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C1_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C1_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c1_p2_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C1_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c1_p2_2_z5.pdf");

	petAnalysis_Param_SiPMMC_C2_Plane0->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane0->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane0->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane0->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane0->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p0_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->RebinY(8);
//	petAnalysis_Param_SiPMMC_C2_Plane2->SetAxisRange(-25.,25.,"X");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_Param_SiPMMC_C2_Plane2->GetYaxis()->SetTitle("pes");
	petAnalysis_Param_SiPMMC_C2_Plane2->Draw("colz");
	c1->Print("deliverables2/param_sipmmc_c2_p2_2_z5.pdf");
	petAnalysis_Param_SiPMMC_C2_Plane2->ProfileX()->Draw();
	c1->Print("deliverables2/prof_sipmmc_c2_p2_2_z5.pdf");
	file->Close();

	c1->Close();
}
