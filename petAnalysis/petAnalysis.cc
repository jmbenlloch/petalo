#include <algorithm> 
#include "barycenterAlgorithm.hh"
#include "findCluster.hh"
#include "TF1.h"
#include<petAnalysis.h>

#define CUT 0.7
#define STEP 5

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
	  ->hman()->h2(this->alabel("xPosZ"),"xRecons-xTrue",100,-25,25,100,-15,15);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yPosZ"),"yRecons-yTrue",100,-25,25,100,-15,15);
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

	//  std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

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
	  findCluster->findCoronnaAllPlanes(planes,clusters,1,80);
	  gate::Point3D reconsPointCoronna; 
	  reconstruction(clusters,reconsPointCoronna);
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCoronna1"), reconsPointCoronna.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCoronna1"), reconsPointCoronna.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCoronna1"), reconsPointCoronna.z() - trueVertex.z());

	  //chargeCluster
	  for(unsigned int i=0;i<6;i++){
		  std::vector<gate::Hit*> hits(evt.GetMCSensHits());
		  gate::Hit* max = *std::max_element(hits.begin(),hits.end(),chargeOrderSensorsAsc);
		  for(unsigned int j=0;j<clusters[i].size();j++){
			  if(max->GetSensorID() != clusters[i][j]->GetSensorID() &&
					  max->GetAmplitude() > 80){
				  gate::Centella::instance()
					  ->hman()->fill(this->alabel("chargeCluster"), clusters[i][j]->GetAmplitude());
			  }
		  }
	  }

	  //Energy in function of z
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosX"), trueVertex.x(), energy-10770);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosZ"), trueVertex.z(), energy-10770);
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
		point[planesCoord[i][0]] += points[i][0] / std::pow(errors[i][0],2);
		point[planesCoord[i][1]] += points[i][1] / std::pow(errors[i][1],2);

		norm[planesCoord[i][0]] += std::pow(errors[i][0],-2);
		norm[planesCoord[i][1]] += std::pow(errors[i][1],-2);
	}

	pt.x(point[0]/norm[0]);
	pt.y(point[1]/norm[1]);
	pt.z(point[2]/norm[2]);
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
