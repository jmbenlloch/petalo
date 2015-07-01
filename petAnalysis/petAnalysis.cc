
#include <algorithm> 
#include "barycenterAlgorithm.hh"
#include "TF1.h"
#include<petAnalysis.h>

#define CUT 0.2

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
    ->hman()->h1(this->alabel("Compton"),"Number of Compton interactions",10,0,10);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("PhotEnergy"),"Photoelectron Energy",30000,0,2);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("XeGammaEnergy"),"Xe Gamma Energy",30000,0,1);

  //gate::Centella::instance()
  // ->hman()->h1(this->alabel("Error"),"Distance from recons. to truth",30000,0,120);

  for(unsigned int i=0;i<6;i++){
	  string histName = "Plane" + gate::to_string(i);
	  string histTitle = "Plane " + gate::to_string(i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(histName),histTitle,30000,0,120);
  }

	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y1"),"Plane 1 y-y0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z1"),"Plane 1 z-z0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x2"),"Plane 2 x-x0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y2"),"Plane 2 y-y0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y3"),"Plane 3 y-y0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z3"),"Plane 3 z-z0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x4"),"Plane 4 x-x0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z4"),"Plane 4 z-z0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x5"),"Plane 5 x-x0",30000,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z5"),"Plane 5 z-z0",30000,-120,120);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("x"),"x-x0",30000,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("y"),"y-y0",30000,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("z"),"z-z0",30000,-120,120);

  for(unsigned int i=0;i<6;i++){
	  string histName = "SiPM" + gate::to_string(i);
	  string histNameRel = "SiPM_Rel" + gate::to_string(i);
	  string histTitle = "SiPM Counts (Plane " + gate::to_string(i) + ")";
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histName),histTitle,100,0,100,100,0,100);
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histNameRel),histTitle,100,0,100,100,0.,1.);
  }
 
  //Hist2d events
/*  store("index",0);
  for(unsigned int i=0;i<100;i++){
	  string histName = "Event" + gate::to_string(i);
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histName),histName,40,0,40,20,0,20);
  }*/

  //gate::Run* run = &gate::Centella::instance()->getRun();
  //int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  //_m.message("Number of generated events in file:",nevt,gate::NORMAL);

  store("photoCount",0);
  store("photoCompton",0);
  store("photoEGamma",0);
  store("photoE",0);
  store("photoWall",0);
  store("comptWall",0);

  store("cut",CUT);

  return true;

}

//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  //Fill energy histogram
  energyHist(evt);

  // Search primary particle and its first daughter
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  //Fill compton hist
  fillComptonHist(primary);

  // Classify event as compton or photoelectric
 // classifyEvent(primary,firstDaughter);

 //Try only events with photoelectric and one vertex
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0){
	  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 
	  //Classify sensor hits per planes
	  std::vector<std::vector<gate::Hit*> > planes(6);
	  splitHitsPerPlane(evt,planes);

	  //Apply cut per plane
	  std::vector<std::vector<gate::Hit*> > planesCut(6);
	  for(unsigned int i=0; i<6;i++){
		  applyCut(planes[i],CUT,planesCut[i]);
	  }

	  //Point Reconstruction
	  gate::Point3D reconsPoint; 
	  reconsPerPlane(planesCut,trueVertex,reconsPoint);  

	  reconstruction(planesCut,reconsPoint);
	  gate::Centella::instance()
		->hman()->fill(this->alabel("x"), reconsPoint.x() - trueVertex.x());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("y"), reconsPoint.y() - trueVertex.y());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("z"), reconsPoint.z() - trueVertex.z());
  }

  //Hist2d to find the cut
  //hist2dHits(evt);

  //Hist2d event
  //hist2dEvent(evt);

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
  std::cout << "Compton in Wall: " << fetch_istore("comptWall") << std::endl;

  // Actual hist name includes algorithm's name
  TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_Energy");
  hist->Rebin(125);
  //double maxV = hist->GetBinCenter( hist->GetMaximumBin() );
  TF1* gauF = new TF1("gauF","gaus",0,10000);
  hist->Fit("gauF","","e",5000,6000);
  std::cout << "FWHM res = " << 2.35*gauF->GetParameter(2)/gauF->GetParameter(1) << std::endl;
  
  return true;

}

//only for 5 faces
void petAnalysis::bestPointRecons(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	double x=0.,y=0.,z=0.,error=0.;
	gate::Point3D auxPt;
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	double y1,z1,x2,y2,y3,z3,x4,z4,x5,z5;
	double y1Err,z1Err,x2Err,y2Err,y3Err,z3Err,x4Err,z4Err,x5Err,z5Err;
	//Plane 1
	barycenter->setPlane("yz");
	barycenter->computePosition(planes[1]);
	y1 = barycenter->getX1();
	z1 = barycenter->getX2();
	y1Err = barycenter->getX1Err();
	z1Err = barycenter->getX2Err();
	//Plane 2
	barycenter->setPlane("xy");
	barycenter->computePosition(planes[1]);
	x2 = barycenter->getX1();
	y2 = barycenter->getX2();
	x2Err = barycenter->getX1Err();
	y2Err = barycenter->getX2Err();
	//Plane 3
	barycenter->setPlane("yz");
	barycenter->computePosition(planes[1]);
	y3 = barycenter->getX1();
	z3 = barycenter->getX2();
	y3Err = barycenter->getX1Err();
	z3Err = barycenter->getX2Err();
	//Plane 4
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[1]);
	x4 = barycenter->getX1();
	z4 = barycenter->getX2();
	x4Err = barycenter->getX1Err();
	z4Err = barycenter->getX2Err();
	//Plane 5
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[1]);
	x5 = barycenter->getX1();
	z5 = barycenter->getX2();
	x5Err = barycenter->getX1Err();
	z5Err = barycenter->getX2Err();

	//Option 1: Planes 1-2
	x = x2;
	y = (y1/std::pow(y1Err,2) + y2/std::pow(y2Err,2)) / (std::pow(y1Err,-2) + std::pow(y1Err,-2));
	z = z1;
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	pt.x(x);
	pt.y(y);
	pt.z(z);
	error = distance(truePt, auxPt);
	
	//Option 2: Planes 1-4
	x = x4;
	y = y1;
	z = (z1/std::pow(z1Err,2) + z4/std::pow(z4Err,2)) / (std::pow(z1Err,-2) + std::pow(z4Err,-2));
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}

	//Option 3: Planes 1-5
	x = x5;
	y = y1;
	z = (z1/std::pow(z1Err,2) + z5/std::pow(z5Err,2)) / (std::pow(z1Err,-2) + std::pow(z5Err,-2));
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}

	//Option 4: Planes 2-3
	x = x2;
	y = (y2/std::pow(y2Err,2) + y3/std::pow(y3Err,2)) / (std::pow(y2Err,-2) + std::pow(y3Err,-2));
	z = z3;
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}

	//Option 5: Planes 2-4
	x = (x2/std::pow(x2Err,2) + x4/std::pow(x4Err,2)) / (std::pow(x2Err,-2) + std::pow(x4Err,-2));
	y = y2;
	z = z4;	
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}

	//Option 6: Planes 2-5
	x = (x2/std::pow(x2Err,2) + x5/std::pow(x5Err,2)) / (std::pow(x2Err,-2) + std::pow(x5Err,-2));
	y = y2;
	z = z5;	
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}

	//Option 7: Planes 3-4
	x = x4;
	y = y3;
	z = (z3/std::pow(z3Err,2) + z4/std::pow(z4Err,2)) / (std::pow(z3Err,-2) + std::pow(z4Err,-2));
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}
	
	//Option 8: Planes 3-5
	x = x5;
	y = y3;
	z = (z3/std::pow(z3Err,2) + z5/std::pow(z5Err,2)) / (std::pow(z3Err,-2) + std::pow(z5Err,-2));
	auxPt.x(x);
	auxPt.y(y);
	auxPt.z(z);
	if(distance(truePt, auxPt) < error){
		error = distance(truePt, auxPt);
		pt.x(x);
		pt.y(y);
		pt.z(z);
	}
}

void petAnalysis::reconsPerPlane(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	//Calculate barycenter
	double x=0.,y=0.,z=0.;
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	// Plane 1
	barycenter->setPlane("yz");
	barycenter->computePosition(planes[1]);
	y = barycenter->getX1();
	z = barycenter->getX2();
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane1"), std::sqrt(std::pow(y - truePt.y(),2)  + std::pow(z - truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("y1"), y - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z1"), z - truePt.z());

	// Plane 2
	barycenter->setPlane("xy");
	barycenter->computePosition(planes[2]);
	x = barycenter->getX1();
	y = barycenter->getX2();
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane2"), std::sqrt(std::pow(x - truePt.x(),2)  + std::pow(y - truePt.y(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("x2"), x - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("y2"), y - truePt.y());

	// Plane 3
	barycenter->setPlane("yz");
	barycenter->computePosition(planes[3]);
	y = barycenter->getX1();
	z = barycenter->getX2();
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane3"), std::sqrt(std::pow(y - truePt.y(),2)  + std::pow(z - truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("y3"), y - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z3"), z - truePt.z());
	
	// Plane 4
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[4]);
	x = barycenter->getX1();
	z = barycenter->getX2();
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane4"), std::sqrt(std::pow(x - truePt.x(),2)  + std::pow(z - truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("x4"), x - truePt.x(),2);
	gate::Centella::instance()
		->hman()->fill(this->alabel("z4"), z - truePt.z());

	// Plane 5
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[4]);
	x = barycenter->getX1();
	z = barycenter->getX2();
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane5"), std::sqrt(std::pow(x - truePt.x(),2)  + std::pow(z - truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("x5"), x - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z5"), z - truePt.z());

}

void petAnalysis::reconstruction(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
	//Calculate barycenter
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
	barycenter->computePosition(planes[1]);
	// std::cout << "Plane 1: y = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
	y += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
	z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
	yNorm += std::pow(barycenter->getX1Err(),-2);
	zNorm += std::pow(barycenter->getX2Err(),-2);
	xNorm += std::pow(barycenter->getX1Err(),-2);
	yNorm += std::pow(barycenter->getX2Err(),-2);
	// Plane 3
	barycenter->setPlane("yz");
	barycenter->computePosition(planes[3]);
	//std::cout << "Plane 3: y = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
	y += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
	z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
	yNorm += std::pow(barycenter->getX1Err(),-2);
	zNorm += std::pow(barycenter->getX2Err(),-2);
	// Plane 4
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[4]);
	//std::cout << "Plane 4: x = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
	x += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
	z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
	xNorm += std::pow(barycenter->getX1Err(),-2);
	zNorm += std::pow(barycenter->getX2Err(),-2);
	// Plane 5
	barycenter->setPlane("xz");
	barycenter->computePosition(planes[4]);
	//std::cout << "Plane 4: x = " << barycenter->getX1() << " ; z = " << barycenter->getX2() << std::endl;
	x += barycenter->getX1() / std::pow(barycenter->getX1Err(),2);
	z += barycenter->getX2() / std::pow(barycenter->getX2Err(),2);
	xNorm += std::pow(barycenter->getX1Err(),-2);
	zNorm += std::pow(barycenter->getX2Err(),-2);

	// Average
	x = x / xNorm;
	y = y / yNorm;
	z = z / zNorm;
	//  std::cout << "Avg: x = " << x << ";\t y = " << y << ";\t z = " << z << std::endl;

	pt.x(x);
	pt.y(y);
	pt.z(z);

	//  std::cout << "Avg: " << pt;
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

void petAnalysis::classifyEvent(gate::MCParticle& primary, gate::MCParticle& firstDaughter){
	int countPhoto = fetch_istore("photoCount");
	int countCompton = fetch_istore("photoCompton");
	int photoEGamma = fetch_istore("photoEGamma");
	int photoE= fetch_istore("photoE");
	int photoWall = fetch_istore("photoWall");
	int comptWall = fetch_istore("comptWall");

	if(firstDaughter.GetCreatorProc() == std::string("compt")){
		countCompton++;
//		std::cout << "\t Compton\n";
		if(primary.GetFinalVol() == "WALL"){
			comptWall++;
		}
/*		std::vector<const gate::MCParticle*> newPart(primary.GetDaughters());
		std::sort(newPart.begin(), newPart.end(), petAnalysis::timeOrderParticles);
		for(unsigned int j=0; j<newPart.size();j++){
		  //std::cout << "\t t: " <<newPart[j]->GetInitialVtx4D().GetT();
		  std::cout << "\t Vol: " <<newPart[j]->GetInitialVol();
		  if(newPart[j]->GetInitialVol() ==  std::string("LXE_DICE")){
			  std::cout << newPart[j]->GetInitialVtx();
		  }
		}
		std::cout << std::endl;
*/
	}else{
		countPhoto++;
		// std::cout << "\t Photo ->";
		for(unsigned int j=0; j<primary.GetDaughters().size();j++){
				  std::cout << "\t PDG/Proc/Energy: " << primary.GetDaughters()[j]->GetPDG() << "/" << primary.GetDaughters()[j]->GetCreatorProc() << "/" << primary.GetDaughters()[j]->GetInitialMom().GetE();

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
		std::cout << std::endl;
		//If only one particle generated then should be energy deposition (due to geant4 IR cut)
		if(primary.GetDaughters().size() == 1){
			//Only one track
			std::cout << "Tracks: " << primary.GetDaughters()[0]->GetTracks().size();
			if(primary.GetDaughters()[0]->GetTracks().size() > 0){
				std::cout << "\t Hits: " << primary.GetDaughters()[0]->GetTracks()[0]->GetHits().size();
				std::cout << "\t HitsEnergy: " << primary.GetDaughters()[0]->GetTracks()[0]->GetHitsEnergy() << std::endl;
			}else{
				std::cout << " noTracks" << primary.GetDaughters()[0]->GetCreatorProc() << std::endl;
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
		}
		//Look for bremsstrahlung
		std::cout << primary.GetDaughters()[0]->GetPDG() << "e- daughters: " << primary.GetDaughters()[0]->GetDaughters().size();
		for(unsigned int k=0; k < primary.GetDaughters()[0]->GetDaughters().size() ; k++){
			std::cout << "\t" << primary.GetDaughters()[0]->GetDaughters()[0]->GetPDG();
			std::cout << "/" << primary.GetDaughters()[0]->GetDaughters()[0]->GetCreatorProc();
			std::cout << "/" << primary.GetDaughters()[0]->GetDaughters()[0]->GetPathLength();
		}
		std::cout << std::endl;
	}

	fstore("photoCount",countPhoto);
	fstore("photoCompton",countCompton);
	fstore("photoEGamma",photoEGamma);
	fstore("photoE",photoE);
	fstore("photoWall",photoWall);
	fstore("comptWall",comptWall);
}

bool petAnalysis::timeOrderParticles(const gate::MCParticle* p1, const gate::MCParticle* p2){
	return (p1->GetInitialVtx4D().GetT() < p2->GetInitialVtx4D().GetT()) ;
}

bool petAnalysis::chargeOrderSensorsDesc(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() > s2->GetAmplitude()) ;
}

bool petAnalysis::chargeOrderSensorsAsc(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() < s2->GetAmplitude()) ;
}

double petAnalysis::distance(gate::Point3D& p1, gate::Point3D& p2){
	double x = (p1.x() - p2.x());
	double y = (p1.y() - p2.y());
	double z = (p1.z() - p2.z());
	return std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
}

void petAnalysis::hist2dEvent(gate::Event& evt){
	std::string histName = "Event" + gate::to_string(fetch_istore("index"));
	int counts[500];
	memset(counts, 0, 500*sizeof(int));
	for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		int id = evt.GetMCSensHits()[i]->GetSensorID();
		counts[(id/1000)*100+(id%100)-100] += evt.GetMCSensHits()[i]->GetAmplitude();
	}
	for(unsigned int i=0;i<100;i++){
		//Plane 1
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),i/10,10-i%10 -0.5,counts[i]);
		//Plane 5
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +10 -0.5, i%10,counts[i+400]);
		//Plane 3
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +20 -0.5, i%10,counts[i+200]);
		//Plane 4
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +30 -0.5, i%10,counts[i+300]);
		//Plane 2
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i%10 +10, 10 - i/10 + 10 -0.5,counts[i+100]);
	}
	fstore("index", fetch_istore("index")+1);
}

void petAnalysis::splitHitsPerPlane(gate::Event& evt, std::vector<std::vector<gate::Hit*> >& planes){
	for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		int id = evt.GetMCSensHits()[i]->GetSensorID();
		if(id < 100){
			planes[0].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 2000){
			planes[1].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 3000){
			planes[2].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 4000){
			planes[3].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 5000){
			planes[4].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 6000){
			planes[5].push_back(evt.GetMCSensHits()[i]);
		}
	}
}

void petAnalysis::hist2dHits(gate::Event& evt){
	std::vector<std::vector<gate::Hit*> >  planes(6);
	splitHitsPerPlane(evt,planes);
	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);

	for(unsigned int i=0; i<6;i++){
		std::sort(sortedPlanes[i].begin(), sortedPlanes[i].end(), petAnalysis::chargeOrderSensorsDesc);
		for(unsigned int j=0; j<sortedPlanes[i].size();j++){
			string histName = "SiPM" + gate::to_string(i);
			gate::Centella::instance()
				->hman()->fill2d(this->alabel(histName),j,sortedPlanes[i][j]->GetAmplitude());

			string histNameRel = "SiPM_Rel" + gate::to_string(i);
			gate::Centella::instance()
				->hman()->fill2d(this->alabel(histNameRel),j, sortedPlanes[i][j]->GetAmplitude() / sortedPlanes[i][0]->GetAmplitude());
		}
	}
}

// Recieve sensor hits and a percentage of the max count.
// Returns hits with amplitude >= cut*maxAmplitude in filtered
void petAnalysis::applyCut(const std::vector<gate::Hit*>& sensorHits, double cut, std::vector<gate::Hit*>& filtered){
	if(sensorHits.size() > 0){
		gate::Hit* max = *std::max_element(sensorHits.begin(),sensorHits.end(),chargeOrderSensorsAsc);

		for(unsigned int i=0; i<sensorHits.size(); i++){
			if(sensorHits[i]->GetAmplitude() >= cut*max->GetAmplitude()){
				filtered.push_back(sensorHits[i]);
			}
		}
	}
}

bool petAnalysis::nearPlane(gate::Point3D& pt, double distance){
	//Plane 0: z = -50.575
	//Plane 1: x = -50.575
	//Plane 2: z = 50.575
	//Plane 3: x = 50.575
	//Plane 4: y = 50.575
	//Plane 5: y = -50.575
	//Formula: D=(ax_0+by_0+cz_0+d)/(sqrt(a^2+b^2+c^2)), 
	/*
	//Planes 1-2
	bool p12 = ((std::abs(pt.x() + 50.575) <= distance) && (std::abs(pt.z() - 50.575) <= distance));
	//Planes 1-4
	bool p14 = ((std::abs(pt.x() + 50.575) <= distance) && (std::abs(pt.y() - 50.575) <= distance));
	//Planes 1-5
	bool p15 = ((std::abs(pt.x() + 50.575) <= distance) && (std::abs(pt.y() + 50.575) <= distance));
	//Planes 2-3
	bool p23 = ((std::abs(pt.z() - 50.575) <= distance) && (std::abs(pt.x() - 50.575) <= distance));
	//Planes 2-4
	bool p24 = ((std::abs(pt.z() - 50.575) <= distance) && (std::abs(pt.y() - 50.575) <= distance));
	//Planes 2-5
	bool p25 = ((std::abs(pt.z() - 50.575) <= distance) && (std::abs(pt.y() + 50.575) <= distance));
	//Planes 3-4
	bool p34 = ((std::abs(pt.x() - 50.575) <= distance) && (std::abs(pt.y() - 50.575) <= distance));
	//Planes 3-5
	bool p35 = ((std::abs(pt.x() - 50.575) <= distance) && (std::abs(pt.y() + 50.575) <= distance));
*/
//	std::cout << "p12: " << p12 << "; p14: " << p14 <<  "; p15: " << p15 
//		<< "; p23: "<< p23 << "; p24: "<< p24 << "; p25: "<< p25
//		<< "; p34: "<< p34 << "; p35: "<< p35 << std::endl;

	//return p12 || p14 || p15 || p23 || p24 || p25 || p34 || p35;

	return (std::abs(pt.x() - 50.575) <= distance) || (std::abs(pt.x() + 50.575) <= distance) ||
	(std::abs(pt.y() - 50.575) <= distance) || (std::abs(pt.y() + 50.575) <= distance) ||
	(std::abs(pt.z() - 50.575) <= distance);

}

void petAnalysis::fillComptonHist(gate::MCParticle& primary){
	int count = 0;
	for(unsigned int i=0; i<primary.GetDaughters().size();i++){
		if(primary.GetDaughters()[i]->GetCreatorProc() == std::string("compt")){
			count++;
		}
	}
	gate::Centella::instance()
		->hman()->fill(this->alabel("Compton"),count);
}

