#include <algorithm> 
#include "barycenterAlgorithm.hh"
#include "TF1.h"
#include<petAnalysis.h>

#define CUT 0.7

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
		  ->hman()->h1(this->alabel(histName),histTitle,100,0,120);
  }

	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x0"),"Plane 0 x-x0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y0"),"Plane 0 y-y0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y1"),"Plane 1 y-y0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z1"),"Plane 1 z-z0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x2"),"Plane 2 x-x0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y2"),"Plane 2 y-y0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("y3"),"Plane 3 y-y0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z3"),"Plane 3 z-z0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x4"),"Plane 4 x-x0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z4"),"Plane 4 z-z0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("x5"),"Plane 5 x-x0",100,-120,120);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel("z5"),"Plane 5 z-z0",100,-120,120);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xNorm"),"x-x0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yNorm"),"y-y0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zNorm"),"z-z0",100,-120,120);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("x"),"x-x0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("y"),"y-y0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("z"),"z-z0",100,-120,120);

  //best
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xbest"),"x-x0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("ybest"),"y-y0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zbest"),"z-z0",100,-120,120);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xbest2Near"),"x-x0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("ybest2Near"),"y-y0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zbest2Near"),"z-z0",100,-120,120);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xbest2NearBySiPM"),"x-x0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("ybest2NearBySiPM"),"y-y0",100,-120,120);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zbest2NearBySiPM"),"z-z0",100,-120,120);

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
  //Only counts comptons from the primary particle
  fillComptonHist(primary);

  // Classify event as compton or photoelectric
 // classifyEvent(primary,firstDaughter);

  //True Vertex
  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 

 //Try only events with photoelectric and one vertex
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0){

      std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

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
	  gate::Point3D reconsPoint2; 
	  gate::Point3D reconsPoint3; 
	  gate::Point3D reconsPoint4; 
	  gate::Point3D reconsPoint5; 
	  reconsPerPlane(planesCut,trueVertex,reconsPoint);  


	  reconstruc2NearestPlanes(planesCut, planes, reconsPoint4);
	  reconstruc2NearestPlanesByMaxSiPM(planesCut, planes, reconsPoint5);


	  reconstruction(planesCut,reconsPoint);
	  reconstructionNoNorm(planesCut,reconsPoint2);
	  bestPointReconsNoNorm(planesCut,trueVertex,reconsPoint3);
	  gate::Centella::instance()
		->hman()->fill(this->alabel("x"), reconsPoint2.x() - trueVertex.x());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("y"), reconsPoint2.y() - trueVertex.y());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("z"), reconsPoint2.z() - trueVertex.z());

	  std::cout << "Norm: x-x0 = " << reconsPoint.x() - trueVertex.x() << "\t y-y0 = " 
		  << reconsPoint.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint.z() - trueVertex.z() << std::endl;
	  std::cout << "NoNorm: x-x0 = " << reconsPoint2.x() - trueVertex.x() << "\t y-y0 = " 
		  << reconsPoint2.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint2.z() - trueVertex.z() << std::endl;
	  std::cout << "Best: x-x0 = " << reconsPoint3.x() - trueVertex.x() << "\t y-y0 = " 
		  << reconsPoint3.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint3.z() - trueVertex.z() << std::endl;
	  std::cout << "Best2Near: x-x0 = " << reconsPoint4.x() - trueVertex.x() << "\t y-y0 = " 
		  << reconsPoint4.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint4.z() - trueVertex.z() << std::endl;
	  std::cout << "Best2Near: x = " << reconsPoint4.x() << "\t y = " << reconsPoint4.y() << "\t z = " << reconsPoint4.z() << std::endl;
	  std::cout << "Best2NearBySiPM: x-x0 = " << reconsPoint5.x() - trueVertex.x() << "\t y-y0 = " 
		  << reconsPoint5.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint5.z() - trueVertex.z() << std::endl;
	  std::cout << "Best2NearBySiPM: x = " << reconsPoint5.x() << "\t y = " << reconsPoint5.y() << "\t z = " << reconsPoint5.z() << std::endl;

	  gate::Centella::instance()
		->hman()->fill(this->alabel("xNorm"), reconsPoint.x() - trueVertex.x());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("yNorm"), reconsPoint.y() - trueVertex.y());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("zNorm"), reconsPoint.z() - trueVertex.z());

	  gate::Centella::instance()
		->hman()->fill(this->alabel("xbest2Near"), reconsPoint4.x() - trueVertex.x());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("ybest2Near"), reconsPoint4.y() - trueVertex.y());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("zbest2Near"), reconsPoint4.z() - trueVertex.z());

	  gate::Centella::instance()
		->hman()->fill(this->alabel("xbest2NearBySiPM"), reconsPoint5.x() - trueVertex.x());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("ybest2NearBySiPM"), reconsPoint5.y() - trueVertex.y());
	  gate::Centella::instance()
		->hman()->fill(this->alabel("zbest2NearBySiPM"), reconsPoint5.z() - trueVertex.z());

//	  printSensors(planesCut);

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

void petAnalysis::reconsPerPlane(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

	std::cout << "x0: " << pointsRecons[0][0] << "\tx0-x = " 
		<< pointsRecons[0][0] - truePt.x() << "\t Var = " << errors[0][0] << std::endl;
	std::cout << "y0: " << pointsRecons[0][1] << "\ty0-y = "
	   	<< pointsRecons[0][1] - truePt.y() << "\t Var = " << errors[0][1] << std::endl;
	std::cout << "y1: " << pointsRecons[1][0] << "\ty1-y = "
	   	<< pointsRecons[1][0] - truePt.y() << "\t Var = " << errors[1][0] << std::endl;
	std::cout << "z1: " << pointsRecons[1][1] << "\tz1-z = "
	   	<< pointsRecons[1][1] - truePt.z() << "\t Var = " << errors[1][1] << std::endl;
	std::cout << "x2: " << pointsRecons[2][0] << "\tx2-x = "
	   	<< pointsRecons[2][0] - truePt.x() << "\t Var = " << errors[2][0] << std::endl;
	std::cout << "y2: " << pointsRecons[2][1] << "\ty2-y = "
	   	<< pointsRecons[2][1] - truePt.y() << "\t Var = " << errors[2][1] << std::endl;
	std::cout << "y3: " << pointsRecons[3][0] << "\ty3-y = "
	   	<< pointsRecons[3][0] - truePt.y() << "\t Var = " << errors[3][0] << std::endl;
	std::cout << "z3: " << pointsRecons[3][1] << "\tz3-z = " 
		<< pointsRecons[3][1] - truePt.z() << "\t Var = " << errors[3][1] << std::endl;
	std::cout << "x4: " << pointsRecons[4][0] << "\tx4-x = "
	   	<< pointsRecons[4][0] - truePt.x() << "\t Var = " << errors[4][0] << std::endl;
	std::cout << "z4: " << pointsRecons[4][1] << "\tz4-z = "
	   	<< pointsRecons[4][1] - truePt.z() << "\t Var = " << errors[4][1] << std::endl;
	std::cout << "x5: " << pointsRecons[5][0] << "\tx5-x = "
	   	<< pointsRecons[5][0] - truePt.x() << "\t Var = " << errors[5][0] << std::endl;
	std::cout << "z5: " << pointsRecons[5][1] << "\tz5-z = "
	   	<< pointsRecons[5][1] - truePt.z() << "\t Var = " << errors[5][1] << std::endl;

	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane0"), std::sqrt(std::pow(pointsRecons[0][0]-truePt.x(),2) + std::pow(pointsRecons[0][1]-truePt.y(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane1"), std::sqrt(std::pow(pointsRecons[1][0]-truePt.y(),2) + std::pow(pointsRecons[1][1]-truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane2"), std::sqrt(std::pow(pointsRecons[2][0]-truePt.x(),2) + std::pow(pointsRecons[2][1]-truePt.y(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane3"), std::sqrt(std::pow(pointsRecons[3][0]-truePt.y(),2) + std::pow(pointsRecons[3][1]-truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane4"), std::sqrt(std::pow(pointsRecons[4][0]-truePt.x(),2) + std::pow(pointsRecons[4][1]-truePt.z(),2)));
	gate::Centella::instance()
		->hman()->fill(this->alabel("Plane5"), std::sqrt(std::pow(pointsRecons[5][0]-truePt.x(),2) + std::pow(pointsRecons[5][1]-truePt.z(),2)));

	gate::Centella::instance()
		->hman()->fill(this->alabel("x0"), pointsRecons[0][0] - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("y0"), pointsRecons[0][1] - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("y1"), pointsRecons[1][0] - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z1"), pointsRecons[1][1] - truePt.z());
	gate::Centella::instance()
		->hman()->fill(this->alabel("x2"), pointsRecons[2][0] - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("y2"), pointsRecons[2][1] - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("y3"), pointsRecons[3][0] - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z3"), pointsRecons[3][1] - truePt.z());
	gate::Centella::instance()
		->hman()->fill(this->alabel("x4"), pointsRecons[4][0] - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z4"), pointsRecons[4][1] - truePt.z());
	gate::Centella::instance()
		->hman()->fill(this->alabel("x5"), pointsRecons[5][0] - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("z5"), pointsRecons[5][1] - truePt.z());
}

void petAnalysis::reconstruction(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
	//Calculate barycenter
	double x=0.,y=0.,z=0.,xNorm=0.,yNorm=0.,zNorm=0.;
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

	// Average
	x += pointsRecons[0][0] / std::pow(errors[0][0],2);
	x += pointsRecons[2][0] / std::pow(errors[2][0],2); 
	x += pointsRecons[4][0] / std::pow(errors[4][0],2);
	x += pointsRecons[5][0] / std::pow(errors[5][0],2);
	
	y += pointsRecons[0][1] / std::pow(errors[0][1],2); 
	y += pointsRecons[1][0] / std::pow(errors[1][0],2);
	y += pointsRecons[2][1] / std::pow(errors[2][1],2);
	y += pointsRecons[3][0] / std::pow(errors[3][0],2);

	z += pointsRecons[1][1] / std::pow(errors[1][1],2);
	z += pointsRecons[3][1] / std::pow(errors[3][1],2);
	z += pointsRecons[4][1] / std::pow(errors[4][1],2);
	z += pointsRecons[5][1] / std::pow(errors[5][1],2);

	xNorm = std::pow(errors[0][0],-2) + std::pow(errors[2][0],-2) + std::pow(errors[4][0],-2) + std::pow(errors[5][0],-2);
	yNorm = std::pow(errors[0][1],-2) + std::pow(errors[1][0],-2) + std::pow(errors[2][1],-2) + std::pow(errors[3][0],-2);
	zNorm = std::pow(errors[1][1],-2) + std::pow(errors[3][1],-2) + std::pow(errors[4][1],-2) + std::pow(errors[5][1],-2);

	pt.x(x/xNorm);
	pt.y(y/yNorm);
	pt.z(z/zNorm);
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
	//Plane 0: z = -50.575 && Plane 1: x = -50.575
	//Plane 2: z =  50.575 && Plane 3: x =  50.575
	//Plane 4: y =  50.575 && Plane 5: y = -50.575
	return (std::abs(pt.x() - 50.575) <= distance) || (std::abs(pt.x() + 50.575) <= distance) ||
	(std::abs(pt.y() - 50.575) <= distance) || (std::abs(pt.y() + 50.575) <= distance) ||
	(std::abs(pt.z() - 50.575) <= distance) || (std::abs(pt.z() + 50.575) <= distance);
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

void petAnalysis::printSensors(std::vector<std::vector<gate::Hit*> >& planes){
	std::cout << "---------- Sensors -----------\n";
	//for(unsigned int i=0;i<6;i++){
	//	for(unsigned int j=0;j<planes[i].size();j++){
	//		std::cout << "SensorID: " << planes[i][j]->GetSensorID() << "\tAmplitude: " << planes[i][j]->GetAmplitude() << std::endl;
	//	}
	//}
	int id;
	double count=0.,row=0.,colTotal=0.,total=0.,rowTotal=0.;
	double colSum[10];
	memset(colSum, 0, 10*sizeof(double));

	std::cout << "---------- Plane 1 (y,z) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = 1000 + i + (9-j)*10;
			count = findSensors(planes[1],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-45+j*10) * count;
			colSum[j] += (45-i*10)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<10;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 10*sizeof(double));

	std::cout << "---------- Plane 2 (x,y) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = 2000 + i*10 + (9-j);
			count = findSensors(planes[2],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-45+j*10) * count;
			colSum[j] += (45-i*10)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<10;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 10*sizeof(double));
	
	std::cout << "---------- Plane 3 (y,z) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = 3000 + (9-i) + (9-j)*10;
			count = findSensors(planes[3],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-45+j*10) * count;
			colSum[j] += (45-i*10)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<10;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 10*sizeof(double));

	std::cout << "---------- Plane 4 (x,z) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = 4000 + (9-i) + j*10;
			count = findSensors(planes[4],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-45+j*10) * count;
			colSum[j] += (45-i*10)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<10;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 10*sizeof(double));

	std::cout << "---------- Plane 5 (x,z) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45\t\tsum" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = 5000 + (9-i) + (9-j)*10;
			count = findSensors(planes[5],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-45+j*10) * count;
			colSum[j] += (45-i*10)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<10;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
}

double petAnalysis::findSensors(std::vector<gate::Hit*>& plane, int id){
	double amplitude = 0;
	for(unsigned int i=0;i<plane.size();i++){
		if(plane[i]->GetSensorID() == id){
			amplitude = plane[i]->GetAmplitude();
	//		std::cout << "(" << plane[i]->GetPosition().x() << "," << plane[i]->GetPosition().y() << "," << plane[i]->GetPosition().z() << ")";
			break;
		}
	}
	return amplitude;
}

void petAnalysis::reconstructionNoNorm(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
	//Calculate barycenter
	double x=0.,y=0.,z=0.;
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

	// Average
	// TODO: Adapt to the number of planes
	x = (pointsRecons[0][0] + pointsRecons[2][0] + pointsRecons[4][0] + pointsRecons[5][0]) / 4.0;
	y = (pointsRecons[0][1] + pointsRecons[1][0] + pointsRecons[2][1] + pointsRecons[3][0]) / 4.0;
	z = (pointsRecons[1][1] + pointsRecons[3][1] + pointsRecons[4][1] + pointsRecons[5][1]) / 4.0;

	pt.x(x);
	pt.y(y);
	pt.z(z);
}

void petAnalysis::bestPointReconsNoNorm(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	double error=0.;
	gate::Point3D auxPt;

	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

	//Select best x
	error = std::abs(pointsRecons[0][0] - truePt.x());
	pt.x(pointsRecons[0][0]);
	if(std::abs(pointsRecons[2][0] - truePt.x()) < error){
		pt.x(pointsRecons[2][0]);
		error = std::abs(pointsRecons[2][0] - truePt.x());
	}
	if(std::abs(pointsRecons[4][0] - truePt.x()) < error){
		pt.x(pointsRecons[4][0]);
		error = std::abs(pointsRecons[4][0] - truePt.x());
	}
	if(std::abs(pointsRecons[5][0] - truePt.x()) < error){
		pt.x(pointsRecons[5][0]);
		error = std::abs(pointsRecons[5][0] - truePt.x());
	}

	//Select best y
	error = std::abs(pointsRecons[0][1] - truePt.y());
	pt.y(pointsRecons[0][1]);
	if(std::abs(pointsRecons[1][0] - truePt.y()) < error){
		pt.y(pointsRecons[1][0]);
		error = std::abs(pointsRecons[1][0] - truePt.y());
	}
	if(std::abs(pointsRecons[2][1] - truePt.y()) < error){
		pt.y(pointsRecons[2][1]);
		error = std::abs(pointsRecons[2][1] - truePt.y());
	}
	if(std::abs(pointsRecons[3][0] - truePt.y()) < error){
		pt.y(pointsRecons[3][0]);
		error = std::abs(pointsRecons[3][0] - truePt.y());
	}

	//Select best z
	error = std::abs(pointsRecons[1][1] - truePt.z());
	pt.z(pointsRecons[1][1]);
	if(std::abs(pointsRecons[3][1] - truePt.z()) < error){
		pt.z(pointsRecons[3][1]);
		error = std::abs(pointsRecons[3][1] - truePt.z());
	}
	if(std::abs(pointsRecons[4][1] - truePt.z()) < error){
		pt.z(pointsRecons[4][1]);
		error = std::abs(pointsRecons[4][1] - truePt.z());
	}
	if(std::abs(pointsRecons[5][1] - truePt.z()) < error){
		pt.z(pointsRecons[5][1]);
		error = std::abs(pointsRecons[5][1] - truePt.z());
	}
	std::cout << std::endl;

	gate::Centella::instance()
		->hman()->fill(this->alabel("xbest"), pt.x() - truePt.x());
	gate::Centella::instance()
		->hman()->fill(this->alabel("ybest"), pt.y() - truePt.y());
	gate::Centella::instance()
		->hman()->fill(this->alabel("zbest"), pt.z() - truePt.z());
}

void petAnalysis::reconstruc2NearestPlanes(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<gate::Hit*> > planesNoCut, gate::Point3D& pt){
	int nonOrthogonal[6] = {2,3,0,1,5,4};
	int planesCoord[6][2] = {{0,1},{1,2},{0,1},{1,2},{0,2},{0,2}};
	double point[3] = {0.,0.,0.};

	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

//	for(unsigned int i=0;i<6;i++){
//		std::cout << "Plane " << i << "\t" << pointsRecons[i][0] << "\t" << pointsRecons[i][1] << std::endl;
//	}

	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);
	std::vector<std::pair<int, double> > planesOrder(6);
	for(unsigned int i=0; i<6; i++){
		planesOrder[i] = std::pair<int, double>(i,totalCharge(planesNoCut[i]));
	//	std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
	}
	std::sort(planesOrder.begin(), planesOrder.end(), petAnalysis::chargeOrderPlanesDesc);

	std::cout << "Ordering... " << std::endl;
	for(unsigned int i=0; i<6; i++){
		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
	}

	int fstPlane = planesOrder[0].first;
	int sndPlane = planesOrder[1].first;

	if(sndPlane == nonOrthogonal[planesOrder[0].first]){
		sndPlane = planesOrder[2].first;
	}
	point[planesCoord[fstPlane][0]] = pointsRecons[fstPlane][0];
	point[planesCoord[fstPlane][1]] = pointsRecons[fstPlane][1];

	if(planesCoord[sndPlane][0] != planesCoord[fstPlane][0] && planesCoord[sndPlane][0] != planesCoord[fstPlane][1]){
		point[planesCoord[sndPlane][0]] = pointsRecons[sndPlane][0];
	}else{
		point[planesCoord[sndPlane][1]] = pointsRecons[sndPlane][1];
	}

	pt.x(point[0]);
	pt.y(point[1]);
	pt.z(point[2]);

	std::cout << "Best planes: " << fstPlane << ", " << sndPlane << std::endl;
}

double petAnalysis::totalCharge(std::vector<gate::Hit*> plane){
	double charge = 0;
	for(unsigned int i=0; i<plane.size();i++){
		charge += plane[i]->GetAmplitude();
	}
	return charge;
}

bool petAnalysis::chargeOrderPlanesDesc(std::pair<int,double> s1, std::pair<int,double> s2){
	return (s1.second > s2.second);
}

void petAnalysis::reconstruc2NearestPlanesByMaxSiPM(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<gate::Hit*> > planesNoCut, gate::Point3D& pt){
	int nonOrthogonal[6] = {2,3,0,1,5,4};
	int planesCoord[6][2] = {{0,1},{1,2},{0,1},{1,2},{0,2},{0,2}};
	double point[3] = {0.,0.,0.};
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

//	for(unsigned int i=0;i<6;i++){
//		std::cout << "Plane " << i << "\t" << pointsRecons[i][0] << "\t" << pointsRecons[i][1] << std::endl;
//	}

	//max charge
	std::vector<std::vector<gate::Hit*> >  sortedSiPM(planes);
	for(unsigned int i=0; i<6; i++){
		std::sort(sortedSiPM[i].begin(), sortedSiPM[i].end(), petAnalysis::chargeOrderSensorsDesc);
	}

	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);
	std::vector<std::pair<int, double> > planesOrder(6);
	for(unsigned int i=0; i<6; i++){
		planesOrder[i] = std::pair<int, double>(i,sortedSiPM[i][0]->GetAmplitude());
		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << " - total: " << totalCharge(planesNoCut[i]) << std::endl;
	}
	std::sort(planesOrder.begin(), planesOrder.end(), petAnalysis::chargeOrderPlanesDesc);

	std::cout << "Ordering... " << std::endl;
	for(unsigned int i=0; i<6; i++){
		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
	}

	int fstPlane = planesOrder[0].first;
	int sndPlane = planesOrder[1].first;

	if(sndPlane == nonOrthogonal[planesOrder[0].first]){
		sndPlane = planesOrder[2].first;
	}
	point[planesCoord[fstPlane][0]] = pointsRecons[fstPlane][0];
	point[planesCoord[fstPlane][1]] = pointsRecons[fstPlane][1];

	if(planesCoord[sndPlane][0] != planesCoord[fstPlane][0] && planesCoord[sndPlane][0] != planesCoord[fstPlane][1]){
		point[planesCoord[sndPlane][0]] = pointsRecons[sndPlane][0];
	}else{
		point[planesCoord[sndPlane][1]] = pointsRecons[sndPlane][1];
	}

	pt.x(point[0]);
	pt.y(point[1]);
	pt.z(point[2]);

	std::cout << "Best planes: " << fstPlane << ", " << sndPlane << std::endl;
}

void petAnalysis::computeBarycenters(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<double> >& points, std::vector<std::vector<double> >& errors){
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	for(unsigned int i=0;i<6;i++){
		barycenter->setPlane(planesDirections[i]);
		barycenter->computePosition(planes[i]);
		points[i][0] = barycenter->getX1();
		points[i][1] = barycenter->getX2();
		errors[i][0] = barycenter->getX1Err();
		errors[i][1] = barycenter->getX2Err();
	}

}
