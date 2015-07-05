void generatePDF() {
   // example of script to loop on all the objects of a ROOT file directory
   // and print on Postscript all TH1 derived objects
   // This script uses the file generated by tutorial hsimple.C
   //Author: Rene Brun

	gStyle->SetOptStat(kFALSE);
   TFile *f1 = TFile::Open("algo_petAnalysis_histos.root");
   TIter next(f1->GetListOfKeys());
   TKey *key;
   TCanvas* c1 = new TCanvas("c1","c1",800,600); 
   //gStyle->SetPalette(53); 
   gStyle->SetPalette(54); 
//   c1->Print("pdf/hist.pdf")
   while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH2")) continue;
      TH2 *h = (TH2*)key->ReadObj();
      h->Draw("colz");
	  std::string filename = "pdf/";
	  filename += key->GetName();
	  filename += ".png";
	  cout << filename << endl;
      c1->Print(filename.c_str());
   }
}

