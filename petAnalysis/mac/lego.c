void lego(){

	gStyle->SetOptStat(kFALSE);

	TCanvas* c1 = new TCanvas("c1", "First canvas", 800,600);
	petAnalysis_phot150->SetTitle("");
	petAnalysis_phot150->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_phot150->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_phot150->GetZaxis()->SetTitle("Q (pe)");
	petAnalysis_phot150->GetXaxis()->SetTitleOffset(2);
	petAnalysis_phot150->GetYaxis()->SetTitleOffset(2);
	petAnalysis_phot150->GetZaxis()->SetTitleOffset(1.5);
	petAnalysis_phot150->Draw("LEGO");
	//c1->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_in_zfar.pdf");
	//c1->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_in_zmiddle.pdf");
	//c1->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_in_znear.pdf");
	
	TCanvas* c2 = new TCanvas("c2", "First canvas", 800,600);
	petAnalysis_phot150out->SetTitle("");
	petAnalysis_phot150out->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_phot150out->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_phot150out->GetZaxis()->SetTitle("Q (pe)");
	petAnalysis_phot150out->GetXaxis()->SetTitleOffset(2);
	petAnalysis_phot150out->GetYaxis()->SetTitleOffset(2);
	petAnalysis_phot150out->GetZaxis()->SetTitleOffset(1.5);
	petAnalysis_phot150out->Draw("LEGO");
	//c2->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_out_zfar.pdf");
	//c2->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_out_zmiddle.pdf");
	//c2->Print("/home/jmbenlloch/universidad/tfm/img/photoelectric_out_znear.pdf");
	
	TCanvas* c3 = new TCanvas("c3", "First canvas", 800,600);
	petAnalysis_compt1->SetTitle("");
	petAnalysis_compt1->GetXaxis()->SetTitle("x (mm)");
	petAnalysis_compt1->GetYaxis()->SetTitle("y (mm)");
	petAnalysis_compt1->GetZaxis()->SetTitle("Q (pe)");
	petAnalysis_compt1->GetXaxis()->SetTitleOffset(2);
	petAnalysis_compt1->GetYaxis()->SetTitleOffset(2);
	petAnalysis_compt1->GetZaxis()->SetTitleOffset(1.5);
	petAnalysis_compt1->Draw("LEGO");
	//c3->Print("/home/jmbenlloch/universidad/tfm/img/compton_2cluster.pdf");
	//c3->Print("/home/jmbenlloch/universidad/tfm/img/compton_1cluster.pdf");
	

}
