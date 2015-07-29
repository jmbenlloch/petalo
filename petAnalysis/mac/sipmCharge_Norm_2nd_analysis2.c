void sipmCharge_Norm_2nd_analysis2(){
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
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->SetAxisRange(-10.,10.,"X");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p0_norm_2nd_2_z2.pdf");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->SetAxisRange(-10.,10.,"X");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p2_norm_2nd_2_z2.pdf");

	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->SetAxisRange(-10.,10.,"X");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p0_norm_2nd_2_z2.pdf");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->SetAxisRange(-10.,10.,"X");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p2_norm_2nd_2_z2.pdf");
	file->Close();

	//////////////
	// LXSC2_Z3 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z3.root");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->SetAxisRange(-15.,15.,"X");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p0_norm_2nd_2_z3.pdf");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->SetAxisRange(-15.,15.,"X");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p2_norm_2nd_2_z3.pdf");

	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->SetAxisRange(-15.,15.,"X");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p0_norm_2nd_2_z3.pdf");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->SetAxisRange(-15.,15.,"X");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p2_norm_2nd_2_z3.pdf");
	file->Close();

	//////////////
	// LXSC2_Z4 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z4.root");

	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->SetAxisRange(-20.,20.,"X");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p0_norm_2nd_2_z4.pdf");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->SetAxisRange(-20.,20.,"X");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p2_norm_2nd_2_z4.pdf");

	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->SetAxisRange(-20.,20.,"X");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p0_norm_2nd_2_z4.pdf");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->SetAxisRange(-20.,20.,"X");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p2_norm_2nd_2_z4.pdf");
	file->Close();
	
	//////////////
	// LXSC2_Z5 //
	//////////////
	TFile *file = TFile::Open("lxsc2_z5.root");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->SetAxisRange(-25.,25.,"X");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p0_norm_2nd_2_z5.pdf");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->SetAxisRange(-25.,25.,"X");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C1_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c1_p2_norm_2nd_2_z5.pdf");

	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->SetAxisRange(-25.,25.,"X");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane0_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p0_norm_2nd_2_z5.pdf");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->SetAxisRange(-25.,25.,"X");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->GetYaxis()->SetTitle("pes");
	petAnalysis_SiPMMC_C2_Plane2_Norm_2nd->Draw("colz");
	c1->Print("deliverables2/sipmmc_c2_p2_norm_2nd_2_z5.pdf");
	file->Close();

	c1->Close();
}
