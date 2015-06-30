
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
    ->hman()->h1(this->alabel("PhotEnergy"),"Photoelectron Energy",30000,0,2);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("XeGammaEnergy"),"Xe Gamma Energy",30000,0,1);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("Error"),"Distance from recons. to truth",30000,0,120);

  //gate::Run* run = &gate::Centella::instance()->getRun();
  //int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  //_m.message("Number of generated events in file:",nevt,gate::NORMAL);

  store("photoCount",0);
  store("photoCompton",0);
  store("photoEGamma",0);
  store("photoE",0);
  store("photoWall",0);

  return true;

}

//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  //Reconstruct point
//  gate::Point3D reconsPoint; 
//  reconstruction(evt,reconsPoint);

  //Fill energy histogram
  energyHist(evt);

  // Classify event as compton or photoelectric
  int countPhoto = fetch_istore("photoCount");
  int countCompton = fetch_istore("photoCompton");
  int photoEGamma = fetch_istore("photoEGamma");
  int photoE= fetch_istore("photoE");
  int photoWall = fetch_istore("photoWall");

  // Search primary particle and its first daughter
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  //True position of the first interaction
//  std::cout << "True: x = " << firstDaughter->GetInitialVtx4D().GetX() << 
//	  ";\t y = " << firstDaughter->GetInitialVtx4D().GetY() << 
//	  ";\t z = " << firstDaughter->GetInitialVtx4D().GetZ() << std::endl;

  // %Error (at zero order...)
/*  double xErr = (reconsPoint.x() - firstDaughter->GetInitialVtx4D().GetX());
  double yErr = (reconsPoint.y() - firstDaughter->GetInitialVtx4D().GetY());
  double zErr = (reconsPoint.z() - firstDaughter->GetInitialVtx4D().GetZ());
  double error = std::sqrt(std::pow(xErr,2) + std::pow(yErr,2) + std::pow(zErr,2));
  //std::cout << "Error: x = " << error << std::endl;
  
  gate::Centella::instance()
    ->hman()->fill(this->alabel("Error"),error);
*/

  if(firstDaughter.GetCreatorProc() == std::string("compt")){
	  countCompton++;
	  //std::cout << "\t Compton\n";
	  //Need to find first interaction
	  //for(unsigned int j=0; j<primary->GetDaughters().size();j++){
		//  std::cout << "\t Vol: " <<primary->GetDaughters()[j]->GetInitialVol();
	  //}
	  //std::cout << std::endl;
  }else{
	  countPhoto++;
	 // std::cout << "\t Photo ->";
	  for(unsigned int j=0; j<primary.GetDaughters().size();j++){
	//	  std::cout << "\t PDG/Proc/Energy: " << primary.GetDaughters()[j]->GetPDG() << "/" << primary.GetDaughters()[j]->GetCreatorProc() << "/" << primary.GetDaughters()[j]->GetInitialMom().GetE();

		  //Electron
		  if(primary.GetDaughters()[j]->GetPDG() == 11){
			  gate::Centella::instance()
				  ->hman()->fill(this->alabel("PhotEnergy"),primary.GetDaughters()[j]->GetInitialMom().GetE());
		  }else if(primary.GetDaughters()[j]->GetPDG() == 22){ //Photon
			  gate::Centella::instance()
				  ->hman()->fill(this->alabel("XeGammaEnergy"),primary.GetDaughters()[j]->GetInitialMom().GetE());
		  }
	   } 
	  if(primary.GetDaughters().size() == 1){
		  photoE++;
		  if(primary.GetFinalVol() == "WALL"){
			  photoWall++;
		  }
	  }else{
		  photoEGamma++;
	  }
/*	  std::cout << std::endl;
	  //If only one particle generated then should be energy deposition (due to geant4 IR cut)
	  if(primary.GetDaughters().size() == 1){
		  //Only one track
		  std::cout << "Tracks: " << primary.GetDaughters()[0]->GetTracks().size();
		  if(primary.GetDaughters()[0]->GetTracks().size() > 0){
			  std::cout << "\t Hits: " << primary.GetDaughters()[0]->GetTracks()[0]->GetHits().size();
			  std::cout << "\t HitsEnergy: " << primary.GetDaughters()[0]->GetTracks()[0]->GetHitsEnergy() << std::endl;
		  }else{
			  std::cout << " noTracks" << primary.GetDaughters()[0]->GetCreatorProc() << " - id:" << evt.GetEventID() << std::endl;
			  std::cout << "Initial x: " << primary.GetDaughters()[0]->GetInitialVtx4D().x() << 
			  " y: " << primary.GetDaughters()[0]->GetInitialVtx4D().y() << 
			  " z: " << primary.GetDaughters()[0]->GetInitialVtx4D().z() <<
			  " t: " << primary.GetDaughters()[0]->GetInitialVtx4D().GetT() << std::endl;
			  std::cout << "Final x: " << primary.GetDaughters()[0]->GetFinalVtx4D().x() << 
			  " y: " << primary.GetDaughters()[0]->GetFinalVtx4D().y() << 
			  " z: " << primary.GetDaughters()[0]->GetFinalVtx4D().z() << 
			  " t: " << primary.GetDaughters()[0]->GetFinalVtx4D().GetT() << std::endl;
			  std::cout << "Initial volume x: " << primary.GetDaughters()[0]->GetInitialVol() << std::endl;
			  std::cout << "Final volume x: " << primary.GetDaughters()[0]->GetFinalVol() << std::endl;
			  std::cout << "Path length: " << primary.GetDaughters()[0]->GetPathLength() << std::endl; 
		  }
	  }*/
  }

  fstore("photoCount",countPhoto);
  fstore("photoCompton",countCompton);

  fstore("photoEGamma",photoEGamma);
  fstore("photoE",photoE);
  fstore("photoWall",photoWall);

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
  std::cout << "Photoelectric with e- gamma: " << fetch_istore("photoEGamma") << std::endl;
  std::cout << "Photoelectric with e-: " << fetch_istore("photoE") << std::endl;
  std::cout << "Photoelectric in Wall: " << fetch_istore("photoWall") << std::endl;

  // Actual hist name includes algorithm's name
  TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_Energy");
  hist->Rebin(125);
  //double maxV = hist->GetBinCenter( hist->GetMaximumBin() );
  TF1* gauF = new TF1("gauF","gaus",0,10000);
  hist->Fit("gauF","","e",5000,6000);
  std::cout << "FWHM res = " << 2.35*gauF->GetParameter(2)/gauF->GetParameter(1) << std::endl;
  
  return true;

}

void petAnalysis::reconstruction(gate::Event& evt, gate::Point3D& pt){
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

  pt.x(x);
  pt.y(y);
  pt.z(z);

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

void petAnalysis::findFirstParticle(const std::vector<const gate::MCParticle*> particles, gate::MCParticle& first){
	double firstTime = 10000000.;
	for(unsigned int i=0; i<particles.size();i++){
		//Search the first daughter
		if(particles[i]->GetInitialVtx4D().GetT() < firstTime){
			first = *particles[i];
			firstTime = particles[i]->GetInitialVtx4D().GetT();
		}
	}
}

void petAnalysis::findFirstParticle(std::vector<gate::MCParticle*> particles, gate::MCParticle& first){
	double firstTime = 10000000.;
	for(unsigned int i=0; i<particles.size();i++){
		//Search the first daughter
		if(particles[i]->GetInitialVtx4D().GetT() < firstTime){
			first = *particles[i];
			firstTime = particles[i]->GetInitialVtx4D().GetT();
		}
	}
}
