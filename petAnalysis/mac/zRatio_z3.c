void zRatio_z3(){
	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_zRatio->GetXaxis()->SetTitle("z (mm)");
	petAnalysis_zRatio->GetYaxis()->SetTitle("ratio");
	petAnalysis_zRatio->GetYaxis()->SetTitleOffset(1.4);
	petAnalysis_zRatio->Draw("colz");

	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	TProfile *prof = petAnalysis_zRatio->ProfileX();
	prof->Fit("pol2");
	prof->Draw();;

	for(unsigned int i=20;i<80;i++){
		std::cout << ", " << prof->GetBinContent(i);
	}
	std::cout << std::endl;
}
