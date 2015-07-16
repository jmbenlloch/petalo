#include <algorithm> 
#include "barycenterAlgorithm.hh"
#include "findCluster.hh"
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
    ->hman()->h1(this->alabel("Energy"),"Energy SiPM",500,0,20000);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("EnergyPhot"),"Energy Phot SiPM",500,0,20000);
  gate::Centella::instance()
    ->hman()->h1(this->alabel("EnergyCompt"),"Energy Compton SiPM",500,0,20000);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("z"),"Event position",30000,-50,50);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("Compton"),"Number of Compton interactions",10,0,10);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("PhotEnergy"),"Photoelectron Energy",30000,0,2);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("XeGammaEnergy"),"Xe Gamma Energy",30000,0,1);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xBest"),"x-x0 best",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yBest"),"y-y0 best",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zBest"),"z-z0 best",100,-50,50);
  
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xCorrected"),"x-x0 best",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yCorrected"),"y-y0 best",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zCorrected"),"z-z0 best",100,-50,50);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xCoronna0"),"x-x0 coronna0",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yCoronna0"),"y-y0 coronna0",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zCoronna0"),"z-z0 coronna0",100,-50,50);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xCoronna1"),"x-x0 coronna1",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yCoronna1"),"y-y0 coronna1",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zCoronna1"),"z-z0 coronna1",100,-50,50);

  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xCoronna2"),"x-x0 coronna2",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yCoronna2"),"y-y0 coronna2",100,-50,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zCoronna2"),"z-z0 coronna2",100,-50,50);

  for(unsigned int i=0;i<6;i++){
	  string histName = "SiPM" + gate::to_string(i);
	  string histNameRel = "SiPM_Rel" + gate::to_string(i);
	  string histTitle = "SiPM Counts (Plane " + gate::to_string(i) + ")";
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histName),histTitle,100,0,100,100,0,100);
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histNameRel),histTitle,100,0,100,100,0.,1.);
  }

  //Charge histograms
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("maxCharge"),"Max charge",1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("avgCharge"),"Average charge",100,0,50);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("avgChargeNoMax"),"Average charge without max",100,0,50);
  
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("chargeCluster"),"Charge cluster",2000,0,2000);

  //2Planes Z recons
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("zRatio"),"Charge Ratio Plane 0 & Plane 2",100,-25,25,100,0,10);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zReconsRatio"),"zRecons-zTrue using ratio",100,-25,25);

  //Position studies
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xPosX"),"xRecons-xTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xPosY"),"xRecons-xTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xPosZ"),"xRecons-xTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yPosY"),"yRecons-yTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yPosZ"),"yRecons-yTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("zPosZ"),"yRecons-yTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("EPosX"),"E_Recons-E_True (x)",100,-25,25,100,-2000,2000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("EPosZ"),"E_Recons-E_True (z)",100,-25,25,100,-2000,2000);

  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xy"),"xy",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xz"),"xz",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yz"),"yz",100,-25,25,100,-25,25);

  //Energy
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xyEnergy"),"xy",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xzEnergy"),"xz",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yzEnergy"),"yz",100,-25,25,100,-25,25);

  store("photoCount",0);
  store("photoCompton",0);
  store("photoEGamma",0);
  store("photoE",0);
  store("photoWall",0);
  store("comptWall",0);

  store("reconsEvents",0);
  store("totalEvents",0);

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

  //Fill charge histograms
  gate::Centella::instance()
	  ->hman()->fill(this->alabel("maxCharge"),GetMaxCharge(evt));
  gate::Centella::instance()
	  ->hman()->fill(this->alabel("avgCharge"),GetAvgCharge(evt));
  gate::Centella::instance()
	  ->hman()->fill(this->alabel("avgChargeNoMax"),GetAvgChargeNoMax(evt));

  //Fill energy hists
  energyPhotCompt(evt);

  //Fill compton hist
  //Only counts comptons from the primary particle
  fillComptonHist(primary);

  // Classify event as compton or photoelectric
 // classifyEvent(primary,firstDaughter);

  //True Vertex
  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 

  //Fill position histograms
  gate::Centella::instance()
	  ->hman()->fill2d(this->alabel("xy"),trueVertex.x(),trueVertex.y());
  gate::Centella::instance()
	  ->hman()->fill2d(this->alabel("xz"),trueVertex.x(),trueVertex.z());
  gate::Centella::instance()
	  ->hman()->fill2d(this->alabel("yz"),trueVertex.y(),trueVertex.z());
  gate::Centella::instance()
	  ->hman()->fill(this->alabel("z"),trueVertex.z());

 //Try only events with photoelectric and one vertex
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0){

//	  std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

	  //Energy
	  int energy = 0;
	  for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		  energy += evt.GetMCSensHits()[i]->GetAmplitude();
	  }

	  //Classify sensor hits per planes
	  std::vector<std::vector<gate::Hit*> > planes(6);
	  splitHitsPerPlane(evt,planes);

	  //zRatio
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("zRatio"),trueVertex.z(),totalCharge(planes[0])/totalCharge(planes[2]));
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zReconsRatio"), zReconsRatio(totalCharge(planes[0])/totalCharge(planes[2])) - trueVertex.z());

	  //Energy
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xyEnergy"),trueVertex.x(),trueVertex.y(),energy);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xzEnergy"),trueVertex.x(),trueVertex.z(),energy);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yzEnergy"),trueVertex.y(),trueVertex.z(),energy);
	  
	  //findCoronna
	  util::findCluster* findCluster = new util::findCluster();
	  std::vector<std::vector<gate::Hit*> > clusters(6);
	  std::vector<std::vector<gate::Hit*> > clusters2(6);
	  std::vector<std::vector<gate::Hit*> > clusters0(6);
	  findCluster->findCoronnaAllPlanes(planes,clusters,1,80,80);
	  findCluster->findCoronnaAllPlanes(planes,clusters2,2,80,80);
	  //findCluster->findCoronnaAllPlanes(planes,clusters0,0,80,80);

	  gate::Point3D reconsPointCoronna; 
	  gate::Point3D reconsPointCoronna0; 
	  gate::Point3D reconsPointCoronna2; 
	  gate::Point3D reconsPointBest; 
	  gate::Point3D reconsPointCorrected; 
	  reconstruction(clusters,reconsPointCoronna);
	  reconstruction(clusters2,reconsPointCoronna2);
	//  reconstruction(clusters0,reconsPointCoronna0);
	  reconstructionCoronna0(planes,reconsPointCoronna0,80,80);
	  reconstructionCorrected(clusters,reconsPointCorrected);
	  bestPointRecons(clusters,trueVertex,reconsPointBest);

/*	  if(std::isnan(reconsPointCoronna.x())){
		  std::cout << "xNan" << std::endl;
	  }
	  if(std::isnan(reconsPointCoronna.y())){
		  std::cout << "yNan" << std::endl;
	  }
	  if(std::isnan(reconsPointCoronna.z())){
		  std::cout << "zNan" << std::endl;
	  }*/

	  fstore("totalEvents",fetch_istore("totalEvents")+1);
	  if(!(std::isnan(reconsPointCoronna.x()) || std::isnan(reconsPointCoronna.y()) | std::isnan(reconsPointCoronna.z()))){
		  fstore("reconsEvents",fetch_istore("reconsEvents")+1);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("xCoronna1"), reconsPointCoronna.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("yCoronna1"), reconsPointCoronna.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("zCoronna1"), reconsPointCoronna.z() - trueVertex.z());
	  }

	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCoronna2"), reconsPointCoronna2.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCoronna2"), reconsPointCoronna2.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCoronna2"), reconsPointCoronna2.z() - trueVertex.z());

	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCoronna0"), reconsPointCoronna0.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCoronna0"), reconsPointCoronna0.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCoronna0"), reconsPointCoronna0.z() - trueVertex.z());

	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xBest"), reconsPointBest.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yBest"), reconsPointBest.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zBest"), reconsPointBest.z() - trueVertex.z());

	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCorrected"), reconsPointCorrected.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCorrected"), reconsPointCorrected.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCorrected"), reconsPointCorrected.z() - trueVertex.z());

	  //chargeCluster
	  for(unsigned int i=0;i<6;i++){
		  std::vector<gate::Hit*> hits(evt.GetMCSensHits());
		  gate::Hit* max = *std::max_element(hits.begin(),hits.end(),chargeOrderSensorsAsc);
		  for(unsigned int j=0;j<clusters[i].size();j++){
			  if(max->GetSensorID() != clusters[i][j]->GetSensorID()){
				  gate::Centella::instance()
					  ->hman()->fill(this->alabel("chargeCluster"), clusters[i][j]->GetAmplitude());
			  }
		  }
	  }

	  //Bias
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosX"), trueVertex.x(), reconsPointBest.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosY"), trueVertex.y(), reconsPointBest.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosZ"), trueVertex.z(), reconsPointBest.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yPosY"), trueVertex.y(), reconsPointBest.y()-trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yPosZ"), trueVertex.z(), reconsPointBest.y()-trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("zPosZ"), trueVertex.z(), reconsPointBest.z()-trueVertex.z());
/*	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosX"), trueVertex.x(), reconsPointCoronna.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosY"), trueVertex.y(), reconsPointCoronna.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xPosZ"), trueVertex.z(), reconsPointCoronna.x()-trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yPosY"), trueVertex.y(), reconsPointCoronna.y()-trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yPosZ"), trueVertex.z(), reconsPointCoronna.y()-trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("zPosZ"), trueVertex.z(), reconsPointCoronna.z()-trueVertex.z());*/

	  //Energy in function of z
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosX"), trueVertex.x(), energy-10770);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosZ"), trueVertex.z(), energy-10770);

	//  printSensors(clusters);
  }

  //Hist2d to find the cut
  hist2dHits(evt);

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
  std::cout << "\n Reconstructed events: " << fetch_istore("reconsEvents") 
	  << "\t Total events: " << fetch_istore("totalEvents")
	  << "\t Efficiency: " << fetch_istore("reconsEvents")/(fetch_istore("totalEvents")*1.0) << std::endl;

  return true;

}

void petAnalysis::reconstruction(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
	//Calculate barycenter
	std::vector<std::vector<double> > points(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	std::vector<int> signalPlanes;
	double point[3] = {0.,0.,0.};
	double norm[3] = {0.,0.,0.};
	int planesCoord[6][2] = {{0,1},{1,2},{0,1},{1,2},{0,2},{0,2}};

	for(unsigned int i=0;i<6;i++){
		if(planes[i].size()>0){
			signalPlanes.push_back(i);

			barycenter->setPlane(planesDirections[i]);
			barycenter->computePosition(planes[i]);
			points[i][0] = barycenter->getX1();
			points[i][1] = barycenter->getX2();
			errors[i][0] = barycenter->getX1Err();
			errors[i][1] = barycenter->getX2Err();
		}
	}
	
	for(unsigned int i=0;i<signalPlanes.size();i++){
		point[planesCoord[signalPlanes[i]][0]] += points[signalPlanes[i]][0] / std::pow(errors[signalPlanes[i]][0],2);
		point[planesCoord[signalPlanes[i]][1]] += points[signalPlanes[i]][1] / std::pow(errors[signalPlanes[i]][1],2);

		norm[planesCoord[signalPlanes[i]][0]] += std::pow(errors[signalPlanes[i]][0],-2);
		norm[planesCoord[signalPlanes[i]][1]] += std::pow(errors[signalPlanes[i]][1],-2);
	}

	pt.x(point[0]/norm[0]);
	pt.y(point[1]/norm[1]);
	pt.z(point[2]/norm[2]);
/*	std::cout << "x: " << point[0] << "\t y: " << point[1] << "\t z: " << point[2] << std::endl;
	std::cout << "xVar: " << norm[0] << "\t yVar: " << norm[1] << "\t zVar: " << norm[2] << std::endl;
	std::cout << "x: " << pt.x() << "\t y: " << pt.y() << "\t z: " << pt.z() << std::endl;*/
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

void petAnalysis::splitHitsPerPlane(gate::Event& evt, std::vector<std::vector<gate::Hit*> >& planes){
	for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		int id = evt.GetMCSensHits()[i]->GetSensorID();
		if(id < 100){
		//	std::cout << "Plane 0, sensor id: " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
			planes[0].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 2000){
		//	std::cout << "Plane 1, sensor id: " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
			planes[1].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 3000){
		//	std::cout << "Plane 2, sensor id: " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
			planes[2].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 4000){
		//	std::cout << "Plane 3, sensor id: " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
			planes[3].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 5000){
		//	std::cout << "Plane 4, sensor id: " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
			planes[4].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 6000){
		//	std::cout << "Plane 5, sensor id " << evt.GetMCSensHits()[i]->GetSensorID() << "\t x: " << evt.GetMCSensHits()[i]->GetPosition().x() << "\t y: " << evt.GetMCSensHits()[i]->GetPosition().y() << "\t z: " << evt.GetMCSensHits()[i]->GetPosition().z()  << std::endl;
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
//	return (std::abs(pt.x() - 50.575) <= distance) || (std::abs(pt.x() + 50.575) <= distance) ||
//	(std::abs(pt.y() - 50.575) <= distance) || (std::abs(pt.y() + 50.575) <= distance) ||
//	(std::abs(pt.z() - 50.575) <= distance) || (std::abs(pt.z() + 50.575) <= distance);
	return (std::abs(pt.z() + 25.575) <= distance);
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

void petAnalysis::bestPointRecons(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	double error=0.;
	gate::Point3D auxPt;
	//Calculate barycenter
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();
	std::vector<double> x,y,z;

	for(unsigned int i=0;i<6;i++){
		if(planes[i].size()>0){
			barycenter->setPlane(planesDirections[i]);
			barycenter->computePosition(planes[i]);
			switch(i){
				case 0:
					x.push_back(barycenter->getX1());
					y.push_back(barycenter->getX2());
					break;
				case 1:
					y.push_back(barycenter->getX1());
					z.push_back(barycenter->getX2());
					break;
				case 2:
					x.push_back(barycenter->getX1());
					y.push_back(barycenter->getX2());
					break;
				case 3:
					y.push_back(barycenter->getX1());
					z.push_back(barycenter->getX2());
					break;
				case 4:
					x.push_back(barycenter->getX1());
					z.push_back(barycenter->getX2());
					break;
				case 5:
					x.push_back(barycenter->getX1());
					z.push_back(barycenter->getX2());
					break;
			}
		}
	}

	error = 100000.;
	for(unsigned int i=0;i<x.size();i++){
		if(std::abs(x[i] - truePt.x()) < error){
			pt.x(x[i]);
			error = std::abs(x[i] - truePt.x());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<y.size();i++){
		if(std::abs(y[i] - truePt.y()) < error){
			pt.y(y[i]);
			error = std::abs(y[i] - truePt.y());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<z.size();i++){
		if(std::abs(z[i] - truePt.z()) < error){
			pt.z(z[i]);
			error = std::abs(z[i] - truePt.z());
		}
	}
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

/*	std::cout << "Plane 0: x=" << points[0][0] << ", y=" << points[0][1] << std::endl;
	std::cout << "Plane 1: y=" << points[1][0] << ", z=" << points[1][1] << std::endl;
	std::cout << "Plane 2: x=" << points[2][0] << ", y=" << points[2][1] << std::endl;
	std::cout << "Plane 3: y=" << points[3][0] << ", z=" << points[3][1] << std::endl;
	std::cout << "Plane 4: x=" << points[4][0] << ", z=" << points[4][1] << std::endl;
	std::cout << "Plane 5: x=" << points[5][0] << ", z=" << points[5][1] << std::endl;
*/
}

double petAnalysis::GetMaxCharge(gate::Event& evt){
  std::vector<gate::Hit*> hits(evt.GetMCSensHits());
  gate::Hit* max = *std::max_element(hits.begin(),hits.end(),chargeOrderSensorsAsc);
  return max->GetAmplitude();
}

double petAnalysis::GetAvgCharge(gate::Event& evt){
	double avg = 0;
	for(unsigned int i=0;i<evt.GetMCSensHits().size();i++){
		avg += evt.GetMCSensHits()[i]->GetAmplitude();
	}
	return avg/evt.GetMCSensHits().size();
}

double petAnalysis::GetAvgChargeNoMax(gate::Event& evt){
	double avg = 0;
	for(unsigned int i=0;i<evt.GetMCSensHits().size();i++){
		avg += evt.GetMCSensHits()[i]->GetAmplitude();
	}
	avg -= GetMaxCharge(evt);
	return avg/(evt.GetMCSensHits().size()-1);
}


void petAnalysis::energyPhotCompt(gate::Event& evt){
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  int energy = 0;
  for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
	  energy += evt.GetMCSensHits()[i]->GetAmplitude();
  }

 // if(firstDaughter.GetCreatorProc() == std::string("phot")){
//	  gate::Centella::instance()
//		  ->hman()->fill(this->alabel("EnergyPhot"),energy);
  //}
  /*
  std::cout << "Event " << evt.GetEventID() <<  " - MCParticles: " << evt.GetMCParticles().size() << std::endl;
  if(firstDaughter.GetCreatorProc() == std::string("compt")){
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("EnergyCompt"),energy);
	  for(unsigned int j=0; j<primary.GetDaughters().size();j++){
		  std::cout << "\t PDG/Proc/Energy: " << primary.GetDaughters()[j]->GetPDG() << "/" << primary.GetDaughters()[j]->GetCreatorProc() << "/" << primary.GetDaughters()[j]->GetInitialMom().GetE() << " - daughters: " << primary.GetDaughters()[j]->GetDaughters().size() << " - t=" << primary.GetDaughters()[j]->GetInitialVtx4D().GetT() << std::endl;

		  const gate::MCParticle* p = primary.GetDaughters()[j];
		  for(unsigned int k=0; k< p->GetDaughters().size();k++){
			  std::cout << "\t\t PDG/Proc/Energy: " << p->GetDaughters()[k]->GetPDG() << 
				  "/" << p->GetDaughters()[k]->GetCreatorProc() << "/" << 
				  p->GetDaughters()[k]->GetInitialMom().GetE() << " - daughters: " << 
				  p->GetDaughters()[k]->GetDaughters().size() << " - t=" << 
				  p->GetDaughters()[k]->GetInitialVtx4D().GetT() << std::endl;
		  }
	  }
  }
*/
  int flagPhot=0;
  for(unsigned int i=0;i<evt.GetMCParticles().size();i++){
	  if(evt.GetMCParticles()[i]->GetCreatorProc() == std::string("phot")){
		  flagPhot = 1;
	//	  gate::Centella::instance()
	//		  ->hman()->fill(this->alabel("EnergyPhot"),energy);
		  break;
	  }
  }
  if(flagPhot==0){
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("EnergyCompt"),energy);
  }

  if(firstDaughter.GetCreatorProc() == std::string("phot")){
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("EnergyPhot"),energy);
  }
}


double petAnalysis::zReconsRatio(double ratio){
	double ratios[100] = {2.73205, 2.73906, 2.66909, 2.61167, 2.52151, 2.49822, 2.44918, 2.417, 2.35749, 2.32905, 2.26496, 2.23821, 2.19439, 2.12249, 2.10054, 2.06435, 2.01884, 1.9875, 1.94561, 1.90444, 1.87033, 1.83465, 1.79233, 1.76525, 1.7309, 1.68712, 1.65806, 1.6348, 1.59805, 1.5388, 1.52833, 1.4927, 1.45451, 1.43733, 1.40809, 1.36951, 1.34272, 1.32054, 1.29409, 1.25984, 1.25199, 1.21908, 1.19365, 1.15915, 1.13903, 1.12939, 1.1031, 1.059, 1.05392, 1.03404, 1.00917, 0.992857, 0.967647, 0.930952, 0.931102, 0.915556, 0.907547, 0.859375, 0.851111, 0.840588, 0.829592, 0.80619, 0.775581, 0.760185, 0.748851, 0.747826, 0.715487, 0.711039, 0.687398, 0.666667, 0.657407, 0.654545, 0.642391, 0.646629, 0.627632, 0.604348, 0.577869, 0.565663, 0.55396, 0.55, 0.55, 0.543333, 0.547297, 0.532716, 0.502747, 0.467978, 0.462766, 0.459211, 0.453704, 0.452667, 0.45, 0.45, 0.44726, 0.441803, 0.435, 0.43, 0.371154, 0.389726, 0.367143, 0.358};
	double z;
	for(unsigned int i=0;i<100;i++){
		if(ratio >= ratios[i]){
			z = -24.75 + 0.5*i;
			break;
		}
	}
	return z;
}

void petAnalysis::reconstructionCorrected(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
	reconstruction(planes,pt);
	double x[100] = {3.35526, 2.96429, 2.26071, 1.82045, 1.26136, 1.05429, 0.772785, 0.15375, 0.0522472, -0.252632, -0.406731, -0.746907, -1.22339, -1.39571, -1.12188, -1.57785, -1.72949, -1.83913, -2.32075, -1.95783, -2.21429, -1.90446, -1.98778, -1.65594, -1.35, -0.535714, -0.323864, -0.523043, -0.4025, -0.278205, -0.374211, -0.240845, -0.284483, -0.164423, -0.208065, -0.0845865, -0.193373, -0.128182, 0.0572727, 0.0466387, 0.0458678, -0.0175, -0.0201031, -0.102778, -0.0268212, -0.00867769, 0.054878, 0.0264706, -0.0516807, -0.0739726, 0.0430435, -0.0233333, 0.0615789, -0.0689189, 0.15, 0.114179, 0.00258621, 0.0293103, 0.075, 0.16978, 0.0209677, 0.0387931, 0.0456522, 0.0769565, 0.304054, 0.233193, 0.063913, 0.333673, 0.211682, 0.11129, 0.408696, 0.396154, 0.442035, 0.532418, 0.864286, 1.74737, 1.8075, 2.355, 1.7434, 2.74907, 2.32982, 1.42304, 2.31207, 2.23919, 1.77651, 1.52037, 1.39488, 0.833168, 0.692857, 0.233505, 0.472222, 0.00185185, -0.105, -0.380435, -0.854545, -1.58115, -2.01, -2.60357, -3.36, -3.4};
	double y[100] = {3.15, 2.72368, 2.53966, 1.93364, 1.51, 0.812069, 0.614286, 0.138235, 0.06, -0.132524, -0.429775, -0.839286, -0.853448, -1.1625, -1.824, -1.57683, -1.63983, -2.04391, -1.85388, -2.42328, -2.38659, -1.99128, -2.42405, -1.65928, -2.01375, -0.7325, -0.511856, -0.111538, -0.0463636, -0.260526, -0.168, -0.53087, -0.0517699, -0.287405, -0.276804, -0.248507, -0.157317, 0.214773, 0.111468, -0.168103, -0.0166667, 0.0163366, -0.107895, 0.127778, -0.0766187, -0.228151, 0.0679487, -0.198387, 0.0138462, -0.0371795, 0.0047619, 0.177027, 0.0525641, -0.114238, 0.146319, 0.0365546, 0.0735849, -0.0616438, -0.115289, 0.0508929, 0.054, 0.163725, 0.015, 0.0151376, 0.176866, -0.00157895, 0.204098, 0.205385, 0.36, 0.168919, 0.583696, 0.00535714, 0.294, 0.383945, 0.767857, 1.38333, 1.926, 2.5, 2.6775, 1.93919, 2.17443, 1.93816, 1.57174, 2.04179, 1.93, 1.90714, 1.44036, 0.838732, 0.644681, 0.321429, 0.1375, 0.066, 0.105172, -0.594304, -0.960638, -1.23197, -1.80455, -2.26154, -2.88214, -2.93571};
	double z[100] = {3.44826, 3.00467, 2.67919, 2.13037, 1.70339, 1.36382, 0.872513, 0.412857, 0.13871, 0.0216418, -0.272513, -0.584759, -1.03757, -1.33, -1.52165, -1.89146, -1.81909, -1.95, -2.78497, -2.17979, -2.25577, -2.15948, -2.07541, -2.41667, -1.97079, -0.747521, -0.567073, -0.586364, -0.378947, -0.256957, -0.0862832, -0.176923, -0.444, 0.0847826, 0.256542, 0.275, -0.0326087, -0.0924658, -0.145652, -0.0233945, -0.332474, -0.0182927, -0.176374, -0.132955, -0.118966, 0.119231, -0.109091, 0.389189, -0.0127119, -0.253125, 0.145714, -0.322727, -0.496667, -0.075, -0.0589286, 0.22619, 0.0663934, 0.334615, 0.291176, 0.173684, -0.00633803, 0.177632, 0.045, -0.17069, -0.124286, 0.337013, 0.144231, 0.0359155, 0.188961, -0.136667, 0.137143, 0.760714, 0.384545, 0.37381, 1.67609, 1.7325, 1.77128, 2.18462, 2.30132, 2.286, 2.0125, 2.49, 1.81849, 1.86892, 1.485, 1.47632, 0.807353, 0.893836, 0.718182, 0.0394737, -0.31, -0.183803, -0.551887, -0.697059, -1.16111, -1.51731, -2.34296, -2.58582, -2.97188, -3.59262};

	TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_xPosX");
	int xIndex = hist->GetXaxis()->FindBin(pt.x());
	int yIndex = hist->GetXaxis()->FindBin(pt.y());
	int zIndex = hist->GetXaxis()->FindBin(pt.z());

//	std::cout << "xPre: " << pt.x() << "\t yPre: " << pt.y() << "\t zPre: " << pt.z() << std::endl;
	pt.x(pt.x() + x[xIndex]);
	pt.y(pt.y() + y[yIndex]);
	pt.z(pt.z() + z[zIndex]);
//	std::cout << "xPost: " << pt.x() << "\t yPost: " << pt.y() << "\t zost: " << pt.z() << std::endl;
}

void petAnalysis::printSensors(std::vector<std::vector<gate::Hit*> >& planes){
	std::cout << "---------- Sensors -----------\n";
	int id;
	double count=0.,row=0.,colTotal=0.,total=0.,rowTotal=0.;
	double colSum[8];
	memset(colSum, 0, 8*sizeof(double));
	int numberSiPM=8;
	double pitch = 6.2;

	std::cout << "---------- Plane 0 (x,y) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<numberSiPM;i++){
		std::cout << (21.7-i*pitch);
		for(int j=0;j<numberSiPM;j++){
			id = i*numberSiPM + j;
			count = findSensors(planes[0],id);
			std::cout << "\tid00" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<8;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 8*sizeof(double));

	std::cout << "---------- Plane 1 (y,z) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<8;i++){
		std::cout << (21.7-i*pitch);
		for(int j=0;j<8;j++){
			id = 1000 + i + (7-j)*8;
			count = findSensors(planes[1],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<8;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 8*sizeof(double));

	std::cout << "---------- Plane 2 (x,y) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<8;i++){
		std::cout << (21.7-i*pitch);
		for(int j=0;j<8;j++){
			id = 2000 + i*8 + (7-j);
			count = findSensors(planes[2],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<8;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 8*sizeof(double));
	
	std::cout << "---------- Plane 3 (y,z) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<8;i++){
		std::cout << (21.7-i*pitch);
		for(int j=0;j<8;j++){
			id = 3000 + (7-i) + (7-j)*8;
			count = findSensors(planes[3],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<8;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 8*sizeof(double));

	std::cout << "---------- Plane 4 (x,z) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<8;i++){
		std::cout << (21.7-i*pitch);
		for(int j=0;j<8;j++){
			id = 4000 + (7-i) + j*8;
			count = findSensors(planes[4],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<8;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
	total=0;
	rowTotal=0.;
	colTotal=0.;
	memset(colSum, 0, 8*sizeof(double));

	std::cout << "---------- Plane 5 (x,z) -----------\n";
	std::cout << "\t-21.7\t\t-15.5\t\t-9.3\t\t-3.1\t\t+3.1\t\t+9.3\t\t+15.5\t\t+21.7" << std::endl;
	for(int i=0;i<8;i++){
		std::cout << (21.7-i*8);
		for(int j=0;j<8;j++){
			id = 5000 + (7-i) + (7-j)*8;
			count = findSensors(planes[5],id);
			std::cout << "\tid" << id << ": " << count;
			total += count;
			row += (-21.7+j*pitch) * count;
			colSum[j] += (21.7-i*pitch)*count;
		}
		std::cout << "\t " << row << std::endl;
		rowTotal += row;
		row=0.;
	}
	std::cout << "Col: ";
	for(int i=0;i<numberSiPM;i++){
		colTotal += colSum[i];
		std::cout << "\t\t" << colSum[i];
	}
	std::cout << std::endl;
	std::cout << "Total sum: " << total << "\t Total Row: " << rowTotal << "\t Total Col: " << colTotal << std::endl;
	std::cout << "Row/Sum: " << rowTotal/total << "\t Col/Sum: " << colTotal/total << std::endl;
}

void petAnalysis::reconstructionCoronna0(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt, double thresholdMax, double thresholdNeighbours){
	//Calculate barycenter
	std::vector<std::vector<double> > points(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	std::vector<int> signalPlanes;
	double point[3] = {0.,0.,0.};
	double norm[3] = {0.,0.,0.};
	int planesCoord[6][2] = {{0,1},{1,2},{0,1},{1,2},{0,2},{0,2}};

	for(unsigned int i=0;i<6;i++){
		barycenter->setPlane(planesDirections[i]);
		bool active = barycenter->computePositionCoronna0(planes[i],thresholdMax,thresholdNeighbours);
		if(active){
			signalPlanes.push_back(i);
			points[i][0] = barycenter->getX1();
			points[i][1] = barycenter->getX2();
			errors[i][0] = barycenter->getX1Err();
			errors[i][1] = barycenter->getX2Err();
		}
	}
	
	for(unsigned int i=0;i<signalPlanes.size();i++){
		point[planesCoord[signalPlanes[i]][0]] += points[signalPlanes[i]][0] / std::pow(errors[signalPlanes[i]][0],2);
		point[planesCoord[signalPlanes[i]][1]] += points[signalPlanes[i]][1] / std::pow(errors[signalPlanes[i]][1],2);

		norm[planesCoord[signalPlanes[i]][0]] += std::pow(errors[signalPlanes[i]][0],-2);
		norm[planesCoord[signalPlanes[i]][1]] += std::pow(errors[signalPlanes[i]][1],-2);
	}

	pt.x(point[0]/norm[0]);
	pt.y(point[1]/norm[1]);
	pt.z(point[2]/norm[2]);
/*	std::cout << "x: " << point[0] << "\t y: " << point[1] << "\t z: " << point[2] << std::endl;
	std::cout << "xVar: " << norm[0] << "\t yVar: " << norm[1] << "\t zVar: " << norm[2] << std::endl;
	std::cout << "x: " << pt.x() << "\t y: " << pt.y() << "\t z: " << pt.z() << std::endl;*/
}
