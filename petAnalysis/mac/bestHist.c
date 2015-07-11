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
	
	petAnalysis_xbestTrue_60->SetTitle("xRecons - xTrue");
	petAnalysis_ybestTrue_60->SetTitle("yRecons - yTrue");
	petAnalysis_zbestTrue_60->SetTitle("zRecons - zTrue");

	petAnalysis_xbestTrue_60->SetAxisRange(-20.,20.,"X");
	petAnalysis_xbestTrue_60->GetXaxis()->SetTitle("#Deltax (mm)");
	petAnalysis_xbestTrue_60->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_xbestTrue_60->GetYaxis()->SetTitle("counts");
//	petAnalysis_xbestTrue_75->SetLineWidth(2);
	petAnalysis_ybestTrue_60->SetAxisRange(-20.,20.,"X");
	petAnalysis_ybestTrue_60->GetXaxis()->SetTitle("#Deltay (mm)");
	petAnalysis_ybestTrue_60->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_ybestTrue_60->GetYaxis()->SetTitle("counts");
	petAnalysis_zbestTrue_60->SetAxisRange(-20.,20.,"X");
	petAnalysis_zbestTrue_60->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zbestTrue_60->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zbestTrue_60->GetYaxis()->SetTitle("counts");

	petAnalysis_xbestTrue_60->Fit("gauF","","e",-15,15);
	petAnalysis_ybestTrue_60->Fit("gauF","","e",-15,15);
	petAnalysis_zbestTrue_60->Fit("gauF","","e",-15,15);

//	c1->Divide(3,0);
//	c1->cd(1);
//	c1->cd(2);
//	c1->cd(3);

	petAnalysis_xbestTrue_60->Draw();
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_ybestTrue_60->Draw();
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_zbestTrue_60->Draw();

/*
	TCanvas* c4 = new TCanvas("c4", "First canvas", 800,600);
	petAnalysis_xCoronna1->Fit("gauF","","e",-15,15);
	petAnalysis_yCoronna1->Fit("gauF","","e",-15,15);
	petAnalysis_zCoronna1->Fit("gauF","","e",-15,15);

	petAnalysis_xCoronna2->Fit("gauF","","e",-15,15);
	petAnalysis_yCoronna2->Fit("gauF","","e",-15,15);
	petAnalysis_zCoronna2->Fit("gauF","","e",-15,15);

	petAnalysis_xCoronna1->Draw();
	TCanvas* c5 = new TCanvas("c5", "First canvas", 800,600);
	petAnalysis_yCoronna1->Draw();
	TCanvas* c6 = new TCanvas("c6", "First canvas", 800,600);
	petAnalysis_zCoronna1->Draw();

	TCanvas* c7 = new TCanvas("c7", "First canvas", 800,600);
	petAnalysis_xCoronna2->Draw();
	TCanvas* c8 = new TCanvas("c8", "First canvas", 800,600);
	petAnalysis_yCoronna2->Draw();
	TCanvas* c9 = new TCanvas("c9", "First canvas", 800,600);
	petAnalysis_zCoronna2->Draw();*/

	TCanvas* c10 = new TCanvas("c10", "First canvas", 800,600);
	petAnalysis_zReconsRatio->SetAxisRange(-20.,20.,"X");
	petAnalysis_zReconsRatio->GetXaxis()->SetTitle("#Deltaz (mm)");
	petAnalysis_zReconsRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zReconsRatio->GetYaxis()->SetTitle("counts");
	petAnalysis_zReconsRatio->Fit("gauF","","e",-15,15);
	petAnalysis_zReconsRatio->Draw();
}
