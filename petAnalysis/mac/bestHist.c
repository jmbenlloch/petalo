 //gStyle->SetPalette(1);

void bestHist(){
//	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1110);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);

	c1->SetFillColor(0);
	c1->SetFrameFillStyle(0);

	TF1* gauF = new TF1("gauF","gaus",-50,50);

/*	petAnalysis_xbestTrue_75->SetAxisRange(-50.,50.,"X");
	petAnalysis_xbestTrue_75->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_ybestTrue_80->SetAxisRange(-50.,50.,"X");
	petAnalysis_ybestTrue_80->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_zbestTrue_85->SetAxisRange(-50.,50.,"X");
	petAnalysis_zbestTrue_85->GetXaxis()->SetTitle("#Deltaz (mm)");
*/
	
	petAnalysis_xbestTrue_75->SetTitle("xRecons - xTrue");
	petAnalysis_ybestTrue_80->SetTitle("yRecons - yTrue");
	petAnalysis_zbestTrue_85->SetTitle("zRecons - zTrue");

	petAnalysis_xbestTrue_75->SetAxisRange(-20.,20.,"X");
	petAnalysis_xbestTrue_75->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xbestTrue_75->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xbestTrue_75->GetYaxis()->SetTitle("counts");
//	petAnalysis_xbestTrue_75->SetLineWidth(2);
	petAnalysis_ybestTrue_80->SetAxisRange(-20.,20.,"X");
	petAnalysis_ybestTrue_80->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_ybestTrue_80->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_ybestTrue_80->GetYaxis()->SetTitle("counts");
	petAnalysis_zbestTrue_85->SetAxisRange(-20.,20.,"X");
	petAnalysis_zbestTrue_85->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zbestTrue_85->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zbestTrue_85->GetYaxis()->SetTitle("counts");

	petAnalysis_xbestTrue_75->Fit("gauF","","e",-15,15);
	petAnalysis_ybestTrue_80->Fit("gauF","","e",-15,15);
	petAnalysis_zbestTrue_85->Fit("gauF","","e",-15,15);

//	c1->Divide(3,0);
//	c1->cd(1);
//	c1->cd(2);
//	c1->cd(3);

	petAnalysis_xbestTrue_75->Draw();
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_ybestTrue_80->Draw();
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_zbestTrue_85->Draw();
}
