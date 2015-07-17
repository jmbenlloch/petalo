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

	  std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

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
	  //reconstruction(clusters,reconsPointCoronna);

	  reconstructionComplete(planes,clusters,reconsPointCoronna);

	  //reconstruction(clusters2,reconsPointCoronna2);
	//  reconstruction(clusters0,reconsPointCoronna0);
	  //reconstructionCoronna0(planes,reconsPointCoronna0,80,80);
	  //reconstructionCorrected(clusters,reconsPointCorrected);
	  //
	  ///////
	  std::vector<std::vector<gate::Hit*> > planes2(6);
	  for(unsigned int i=0;i<clusters.size();i++){
		  if(clusters[i].size() > 0){
			  planes2[i] = planes[i];
		  }
	  }
	  bestPointRecons(planes2,trueVertex,reconsPointBest);
	  reconstructionCorrected(clusters,reconsPointBest);
	  //////
	  //bestPointRecons(clusters,trueVertex,reconsPointBest);

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
	  if(!(std::isnan(reconsPointCoronna.x()) || std::isnan(reconsPointCoronna.y()) || std::isnan(reconsPointCoronna.z()))){
//	  if(!(std::isnan(reconsPointCoronna.x()) || std::isnan(reconsPointCoronna.y()))){
		  fstore("reconsEvents",fetch_istore("reconsEvents")+1);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("xCoronna1"), reconsPointCoronna.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("yCoronna1"), reconsPointCoronna.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("zCoronna1"), reconsPointCoronna.z() - trueVertex.z());
	  }
	  std::cout << "xC1: " << reconsPointCoronna.x() - trueVertex.x() << std::endl;
	  std::cout << "yC1: " << reconsPointCoronna.y() - trueVertex.y() << std::endl;
	  std::cout << "zC1: " << reconsPointCoronna.z() - trueVertex.z() << std::endl;

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

	  if(!(std::isnan(reconsPointCoronna.x()) || std::isnan(reconsPointCoronna.y()) || std::isnan(reconsPointCoronna.z()))){
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("xBest"), reconsPointBest.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("yBest"), reconsPointBest.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel("zBest"), reconsPointBest.z() - trueVertex.z());
	  }
	  std::cout << "xBest: " << reconsPointBest.x() - trueVertex.x() << std::endl;
	  std::cout << "yBest: " << reconsPointBest.y() - trueVertex.y() << std::endl;
	  std::cout << "zBest: " << reconsPointBest.z() - trueVertex.z() << std::endl;

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
	  if(!(std::isnan(reconsPointCoronna.x()) || std::isnan(reconsPointCoronna.y()) || std::isnan(reconsPointCoronna.z()))){
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
		  /*gate::Centella::instance()
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
	  }

	  //Energy in function of z
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosX"), trueVertex.x(), energy-10770);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosZ"), trueVertex.z(), energy-10770);

	  //printSensors(clusters);
	//  printSensors(planes);
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
//			std::cout << "Plane " << i << " active " << planesDirections[i] << std::endl;

			signalPlanes.push_back(i);

			barycenter->setPlane(planesDirections[i]);
			barycenter->computePosition(planes[i]);
			points[i][0] = barycenter->getX1();
			points[i][1] = barycenter->getX2();
			errors[i][0] = barycenter->getX1Err();
			errors[i][1] = barycenter->getX2Err();

//			std::cout << "x1: " << points[i][0] << " Var: " <<  errors[i][0] << std::endl;
//			std::cout << "x2: " << points[i][1] << " Var: " <<  errors[i][1] << std::endl;
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
			std::cout << "one more plane" << i << std::endl;
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
		std::cout << "\t x: " << x[i] << std::endl;
		if(std::abs(x[i] - truePt.x()) < error){
			pt.x(x[i]);
			error = std::abs(x[i] - truePt.x());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<y.size();i++){
		std::cout << "\t y: " << y[i] << std::endl;
		if(std::abs(y[i] - truePt.y()) < error){
			pt.y(y[i]);
			error = std::abs(y[i] - truePt.y());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<z.size();i++){
		std::cout << "\t z: " << z[i] << std::endl;
		if(std::abs(z[i] - truePt.z()) < error){
			pt.z(z[i]);
			error = std::abs(z[i] - truePt.z());
		}
	}

	std::cout << "Best x: " << pt.x() << "\t y: " << pt.y() << "\t z: " << pt.z() << std::endl;
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
	//reconstruction(planes,pt);
	double x[100] = {12, 12.15, 11.4879, 10.9, 12.1891, 11.5319, 10.0033, 9.63353, 8.475, 7.81579, 7.87358, 7.55286, 8.1619, 7.91803, 6.88302, 6.77022, 6.15, 6.06735, 5.83611, 5.41429, 5.36224, 4.98806, 4.64795, 4.37586, 4.35594, 3.978, 3.9934, 3.62209, 3.23731, 2.96405, 2.80385, 2.61761, 2.35714, 2.17358, 2.1063, 2.00693, 1.75357, 1.57973, 1.5698, 1.48486, 1.29911, 1.05, 0.973333, 0.771649, 0.591322, 0.55, 0.295455, 0.219, 0.0589888, 0.035, 0.0159574, -0.0857143, -0.207353, -0.311765, -0.511765, -0.702101, -0.891573, -0.906522, -1.16803, -1.26404, -1.43625, -1.53056, -1.72394, -1.84429, -1.67182, -2.04174, -2.31667, -2.56402, -2.81923, -3.08697, -3.12078, -3.20323, -3.3622, -4.01422, -4.07804, -4.11857, -4.48767, -5.06584, -4.82107, -5.6819, -5.62563, -5.7, -6.12025, -7.134, -6.96522, -7.51364, -7.45522, -7.90854, -7.98478, -8.38056, -9.31, -10.0995, -9.96045, -10.6576, -10.98, -10.2946, -10.5763, -10.1283, -11.8607, -11.01};
	double y[100] = { 0, 13.75, 11.3192, 10.6827, 10.4437, 10.0034, 10.081, 10.4089, 9.39615, 9.24231, 9.1067, 8.33842, 7.67256, 6.96176, 7.73919, 7.101, 6.4703, 6.40227, 5.81557, 5.90682, 5.2719, 4.73571, 4.92059, 4.5925, 4.31907, 3.87429, 3.62339, 3.3468, 3.33425, 3.06509, 3.03273, 2.55496, 2.46933, 2.25882, 1.79536, 1.98835, 1.78902, 1.72253, 1.54143, 1.37797, 1.32303, 1.11923, 0.992308, 0.758108, 0.587615, 0.309783, 0.295545, 0.242, 0.161811, -0.012069, 0.0472222, -0.0675824, -0.198913, -0.316957, -0.414179, -0.583696, -0.84375, -1.025, -1.18947, -1.32624, -1.5, -1.61165, -1.53367, -1.66121, -1.64586, -2.11634, -2.45149, -2.26957, -2.60693, -2.9388, -3.22941, -3.20122, -3.6325, -3.80081, -3.88404, -3.93443, -4.48404, -4.63767, -5.10506, -4.96463, -5.68061, -6.04371, -5.7197, -6.41096, -7.20116, -7.52857, -7.13313, -7.44474, -7.90922, -7.989, -7.91, -8.502, -10.066, -9.76111, -10.7541, -10.8392, -10.3263, -10.8346, -12.2022, -12.51};
	double z[100] = {10.8273, 11.6021, 11.1089, 10.9844, 10.8504, 10.493, 9.83265, 9.63533, 9.58622, 8.874, 8.4709, 8.01089, 8.1956, 7.67775, 7.13182, 7.34534, 6.67857, 6.60938, 6.41369, 6.09431, 5.65164, 5.62762, 5.02852, 4.77, 5.0325, 4.21887, 4.01918, 3.90469, 3.8125, 3.45, 3.25714, 3.09545, 2.79393, 2.5883, 2.4283, 2.25345, 2.10789, 1.99521, 1.68043, 1.6067, 1.45141, 1.35517, 1.14583, 1.01418, 0.866129, 0.597458, 0.45, 0.361475, 0.0954545, 0.125676, -0.0676471, 0.00849057, -0.132857, -0.454545, -0.482432, -0.772642, -0.897059, -1.08158, -1.19211, -1.32143, -1.50238, -1.74, -1.6, -1.82, -2.04565, -2.00921, -2.43103, -2.5911, -2.72931, -3.09286, -3.2981, -3.50217, -3.60692, -3.81462, -4.63235, -4.34492, -5.20574, -4.90851, -5.15333, -5.27561, -5.34623, -6.16286, -6.2537, -6.70465, -6.43621, -7.43427, -6.46644, -7.68462, -7.68529, -7.46579, -7.77923, -8.57466, -10.0254, -10.5868, -9.85, -9.71939, -9.82273, -10.5854, -11.8914, -11.406};

	TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_xPosX");
	int xIndex = hist->GetXaxis()->FindBin(pt.x());
	int yIndex = hist->GetXaxis()->FindBin(pt.y());
	int zIndex = hist->GetXaxis()->FindBin(pt.z());

	std::cout << "xPre: " << pt.x() << "\t yPre: " << pt.y() << "\t zPre: " << pt.z() << std::endl;
	pt.x(pt.x() - x[xIndex]);
	pt.y(pt.y() - y[yIndex]);
	pt.z(pt.z() - z[zIndex]);
	std::cout << "xPost: " << pt.x() << "\t yPost: " << pt.y() << "\t zPost: " << pt.z() << std::endl;
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

void petAnalysis::reconstructionComplete(std::vector<std::vector<gate::Hit*> > planesComplete, std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt){
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
			std::cout << "Plane " << i << " active " << planesDirections[i] << std::endl;

			signalPlanes.push_back(i);

			barycenter->setPlane(planesDirections[i]);
			barycenter->computePosition(planesComplete[i]);
			points[i][0] = barycenter->getX1();
			points[i][1] = barycenter->getX2();
			errors[i][0] = barycenter->getX1Err();
			errors[i][1] = barycenter->getX2Err();

			std::cout << "x1: " << points[i][0] << " Var: " <<  errors[i][0] << std::endl;
			std::cout << "x2: " << points[i][1] << " Var: " <<  errors[i][1] << std::endl;
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
	std::cout << "x: " << point[0] << "\t y: " << point[1] << "\t z: " << point[2] << std::endl;
	std::cout << "xVar: " << norm[0] << "\t yVar: " << norm[1] << "\t zVar: " << norm[2] << std::endl;
	std::cout << "x: " << pt.x() << "\t y: " << pt.y() << "\t z: " << pt.z() << std::endl;
}
