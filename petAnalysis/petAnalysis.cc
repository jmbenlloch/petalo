
#include "barycenterAlgorithm.hh"
#include "TF1.h"
#include<petAnalysis.h>

ClassImp(petAnalysis)

//==========================================================================
petAnalysis::petAnalysis(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"petAnalysis",0,label){
//==========================================================================


}

//==========================================================================
petAnalysis::petAnalysis(const gate::ParamStore& gs, 
			   gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"petAnalysis",0,label){
//==========================================================================


}

//==========================================================================
bool petAnalysis::initialize(){
//==========================================================================

  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  gate::Centella::instance()
    ->hman()->h1(this->alabel("Energy"),"Energy SiPM",30000,0,10000);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("Error"),"Distance from recons. to truth",30000,0,120);

  //gate::Run* run = &gate::Centella::instance()->getRun();
  //int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  //_m.message("Number of generated events in file:",nevt,gate::NORMAL);

  store("photoCount",0);
  store("photoCompton",0);

  return true;

}

//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  // Classify sensor hits per plane
  std::vector<gate::Hit*> plane0,plane1,plane2,plane3,plane4,plane5;

  for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
	  int id = evt.GetMCSensHits()[i]->GetSensorID();

	  if(id < 100){
		  plane0.push_back(evt.GetMCSensHits()[i]);
		  //std::cout << "Plane 0: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }else if(id < 2000){
		  plane1.push_back(evt.GetMCSensHits()[i]);
		  //std::cout << "Plane 1: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }else if(id < 3000){
		  plane2.push_back(evt.GetMCSensHits()[i]);
		  //std::cout << "Plane 2: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }else if(id < 4000){
		  plane3.push_back(evt.GetMCSensHits()[i]);
		  //std::cout << "Plane 3: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }else if(id < 5000){
		  plane4.push_back(evt.GetMCSensHits()[i]);
		  //std::cout << "Plane 4: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }else if(id < 6000){
		  plane5.push_back(evt.GetMCSensHits()[i]);
		//  std::cout << "Plane 5: " << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().z() << std::endl;
	  }
//	  std::cout << evt.GetMCSensHits()[i]->GetPosition().x() << ", \t" << evt.GetMCSensHits()[i]->GetPosition().y() << std::endl;
  }

  double x=0.,y=0.,z=0.,xNorm=0.,yNorm=0.,zNorm=0.;
  util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();
  // Plane 0
/*  barycenter->setPlane("yz");
  barycenter->computePosition(plane0);
  std::cout << "Plane 0: y = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
  y += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  yNorm += std::pow(barycenter->getX1Err(),-2);
  zNorm += std::pow(barycenter->getX2Err(),-2);*/
  // Plane 1
  barycenter->setPlane("yz");
  barycenter->computePosition(plane1);
 // std::cout << "Plane 1: y = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
  y += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  yNorm += std::pow(barycenter->getX1Err(),-2);
  zNorm += std::pow(barycenter->getX2Err(),-2);
  // Plane 2
  barycenter->setPlane("xy");
  barycenter->computePosition(plane2);
  //std::cout << "Plane 2: x = " << barycenter->getX1() << " ; y = " << barycenter->getX2() << std::endl;
  x += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  y += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  xNorm += std::pow(barycenter->getX1Err(),-2);
  yNorm += std::pow(barycenter->getX2Err(),-2);
  // Plane 3
  barycenter->setPlane("yz");
  barycenter->computePosition(plane3);
  //std::cout << "Plane 3: y = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
  y += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  yNorm += std::pow(barycenter->getX1Err(),-2);
  zNorm += std::pow(barycenter->getX2Err(),-2);
  // Plane 4
  barycenter->setPlane("xz");
  barycenter->computePosition(plane4);
  //std::cout << "Plane 4: x = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
  x += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  xNorm += std::pow(barycenter->getX1Err(),-2);
  zNorm += std::pow(barycenter->getX2Err(),-2);
  // Plane 5
  barycenter->setPlane("xz");
  barycenter->computePosition(plane4);
  //std::cout << "Plane 4: x = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
  x += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
  z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
  xNorm += std::pow(barycenter->getX1Err(),-2);
  zNorm += std::pow(barycenter->getX2Err(),-2);

  // Average
  x = x / xNorm;
  y = y / yNorm;
  z = z / zNorm;
  std::cout << "Avg: x = " << x << ";\t y = " << y << ";\t z = " << z << std::endl;

  //Fill energy histogram
  energyHist(evt);

  // Classify event as compton or photoelectric
  int countPhoto = fetch_istore("photoCount");
  int countCompton = fetch_istore("photoCompton");
  int comptFlag = 0;

  gate::MCParticle* primary;
  const gate::MCParticle* firstDaughter;
  double firstDaughterTime = 10000000.;

  //std::cout << "Particles: "<< evt.GetMCParticles().size() << "\n";
  for(unsigned int i=0;i<evt.GetMCParticles().size();i++){
	  gate::MCParticle* particle = evt.GetMCParticles()[i];
	  if(particle->IsPrimary()){

		  primary = particle;

//		  std::cout << "Primary daughters: " << particle->GetDaughters().size();
		  //std::vector<gate::MCParticle*> daughters = particle->GetDaughters();
		  for(unsigned int j=0; j<particle->GetDaughters().size();j++){
			  if(particle->GetDaughters()[j]->GetCreatorProc() == std::string("compt")){
				 comptFlag = 1;
//				 break;
			  }
	//		  std::cout << "Time: " << particle->GetDaughters()[j]->GetInitialVtx4D().GetT() <<  ";\t First Time: " << firstDaughterTime << std::endl;
			  //Search the first daughter
			  if(particle->GetDaughters()[j]->GetInitialVtx4D().GetT() < firstDaughterTime){
				  firstDaughter = particle->GetDaughters()[j];
				  firstDaughterTime = particle->GetDaughters()[j]->GetInitialVtx4D().GetT();
	//			  std::cout << "this\n";
			  }
		  }
	  }
  }

  //True position of the first interaction
//  std::cout << "True: x = " << firstDaughter->GetInitialVtx4D().GetX() << 
//	  ";\t y = " << firstDaughter->GetInitialVtx4D().GetY() << 
//	  ";\t z = " << firstDaughter->GetInitialVtx4D().GetZ() << std::endl;

  // %Error (at zero order...)
  double xErr = (x - firstDaughter->GetInitialVtx4D().GetX());
  double yErr = (y - firstDaughter->GetInitialVtx4D().GetY());
  double zErr = (z - firstDaughter->GetInitialVtx4D().GetZ());
  double error = std::sqrt(std::pow(xErr,2) + std::pow(yErr,2) + std::pow(zErr,2));
  //std::cout << "Error: x = " << error << std::endl;
  
  gate::Centella::instance()
    ->hman()->fill(this->alabel("Error"),error);

  if(comptFlag){
	  countCompton++;
	  //std::cout << "\t Compton\n";
	  //Need to find first interaction
	  //for(unsigned int j=0; j<primary->GetDaughters().size();j++){
		//  std::cout << "\t Vol: " <<primary->GetDaughters()[j]->GetInitialVol();
	  //}
	  //std::cout << std::endl;
  }else{
	  countPhoto++;
	 /* std::cout << "\t Photo ->";
	  for(unsigned int j=0; j<primary->GetDaughters().size();j++){
		  std::cout << "\t PDG/Proc/Energy: " << primary->GetDaughters()[j]->GetPDG() << "/" << primary->GetDaughters()[j]->GetCreatorProc() << "/" << primary->GetDaughters()[j]->GetInitialMom().GetE();
	   } 
	  std::cout << std::endl;
	  //If only one particle generated then should be energy deposition (due to geant4 IR cut)
	  if(primary->GetDaughters().size() == 1){
		  //Only one track
		  std::cout << "Tracks: " << primary->GetDaughters()[0]->GetTracks().size();
		  std::cout << "\t Hits: " << primary->GetDaughters()[0]->GetTracks()[0]->GetHits().size();
		  std::cout << "\t HitsEnergy: " << primary->GetDaughters()[0]->GetTracks()[0]->GetHitsEnergy() << std::endl;
	  }*/
  }

  fstore("photoCount",countPhoto);
  fstore("photoCompton",countCompton);

  return true;

}

//==========================================================================
bool petAnalysis::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);

  gate::Run* run = &gate::Centella::instance()->getRun();
  int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  _m.message("Number of generated events in file:",nevt,gate::NORMAL);

  std::cout << "Events with photoelectric first: " << fetch_istore("photoCount") << std::endl;
  std::cout << "Events with compton first: " << fetch_istore("photoCompton") << std::endl;

  // Actual hist name includes algorithm's name
  TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_Energy");
  hist->Rebin(125);
  //double maxV = hist->GetBinCenter( hist->GetMaximumBin() );
  TF1* gauF = new TF1("gauF","gaus",0,10000);
  hist->Fit("gauF","","e",5000,6000);
  std::cout << "FWHM res = " << 2.35*gauF->GetParameter(2)/gauF->GetParameter(1) << std::endl;
  
  return true;

}

void petAnalysis::energyHist(gate::Event& evt){
  // Energy histogram
  int energy = 0;
  for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
	  energy += evt.GetMCSensHits()[i]->GetAmplitude();
  }
  
  gate::Centella::instance()
    ->hman()->fill(this->alabel("Energy"),energy);
}
