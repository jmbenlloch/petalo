 //gStyle->SetPalette(1);

void planesHist5(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 1000,500);

	petAnalysis_Plane0->SetAxisRange(0.,60.,"X");
	petAnalysis_Plane1->SetAxisRange(0.,60.,"X");
	petAnalysis_Plane2->SetAxisRange(0.,60.,"X");
	petAnalysis_Plane3->SetAxisRange(0.,60.,"X");
	petAnalysis_Plane4->SetAxisRange(0.,60.,"X");
	petAnalysis_Plane5->SetAxisRange(0.,60.,"X");

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	c1->Divide(3,2);
	c1->cd(1);
//	petAnalysis_Plane0->Draw();
	c1->cd(2);
	petAnalysis_Plane1->Draw();
	c1->cd(3);
	petAnalysis_Plane2->Draw();
	c1->cd(4);
	petAnalysis_Plane3->Draw();
	c1->cd(5);
	petAnalysis_Plane4->Draw();
	c1->cd(6);
	petAnalysis_Plane5->Draw();
}
