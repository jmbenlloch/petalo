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
    ->hman()->h1(this->alabel("correctSiPM"),"Max SiPM charge in true vertex projection",10,0,10);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("correctSiPMByMax"),"Max SiPM charge in true vertex projection (2 best planes)",10,0,10);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("correctSiPMByTotal"),"Max SiPM charge in true vertex projection (2 best planes)",10,0,10);

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


  for(unsigned int i=0;i<20;i++){
	  string namexNorm = "xNorm_" + gate::to_string(5*i);
	  string nameyNorm = "yNorm_" + gate::to_string(5*i);
	  string namezNorm = "zNorm_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namexNorm),"x-x0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameyNorm),"y-y0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namezNorm),"z-z0",100,-50,50);

	  string namexNoNorm = "xNoNorm_" + gate::to_string(5*i);
	  string nameyNoNorm = "yNoNorm_" + gate::to_string(5*i);
	  string namezNoNorm = "zNoNorm_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namexNoNorm),"x-x0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameyNoNorm),"y-y0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namezNoNorm),"z-z0",100,-50,50);

	  string namexbestTrue = "xbestTrue_" + gate::to_string(5*i);
	  string nameybestTrue = "ybestTrue_" + gate::to_string(5*i);
	  string namezbestTrue = "zbestTrue_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namexbestTrue),"x-x0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameybestTrue),"y-y0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namezbestTrue),"z-z0",100,-50,50);

	  string namex2Near = "x2Near_" + gate::to_string(5*i);
	  string namey2Near = "y2Near_" + gate::to_string(5*i);
	  string namez2Near = "z2Near_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namex2Near),"x-x0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namey2Near),"y-y0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namez2Near),"z-z0",100,-50,50);

	  string namex2NearBySiPM = "x2NearBySiPM_" + gate::to_string(5*i);
	  string namey2NearBySiPM = "y2NearBySiPM_" + gate::to_string(5*i);
	  string namez2NearBySiPM = "z2NearBySiPM_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namex2NearBySiPM),"x-x0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namey2NearBySiPM),"y-y0",100,-50,50);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(namez2NearBySiPM),"z-z0",100,-50,50);
  }

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

  //Hist2d events
/*  store("index",0);
  for(unsigned int i=0;i<100;i++){
	  string histName = "Event" + gate::to_string(i);
	  gate::Centella::instance()
		  ->hman()->h2(this->alabel(histName),histName,40,0,40,30,0,30);
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

/*  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0
		  && nearPlane(trueVertex,10)){*/

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
	  

	  // First centimeter
//	  if(nearPlane(trueVertex,10)){

	  //TODO test findCoronna
	  util::findCluster* findCluster = new util::findCluster();
	  std::vector<std::vector<gate::Hit*> > clusters(6);
	  findCluster->findCoronnaAllPlanes(planes,clusters,1);
	  gate::Point3D reconsPointCoronna; 
	  bestPointRecons(clusters,trueVertex,reconsPointCoronna);
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCoronna1"), reconsPointCoronna.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCoronna1"), reconsPointCoronna.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCoronna1"), reconsPointCoronna.z() - trueVertex.z());

	  //2ring
	  std::vector<std::vector<gate::Hit*> > clusters2(6);
	  findCluster->findCoronnaAllPlanes(planes,clusters2,2);
	  gate::Point3D reconsPointCoronna2; 
	  bestPointRecons(clusters2,trueVertex,reconsPointCoronna2);
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("xCoronna2"), reconsPointCoronna2.x() - trueVertex.x());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("yCoronna2"), reconsPointCoronna2.y() - trueVertex.y());
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zCoronna2"), reconsPointCoronna2.z() - trueVertex.z());

	  //Energy in function of z
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosX"), trueVertex.x(), energy-9641)
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("EPosZ"), trueVertex.z(), energy-9641);

	  //Find best cut
	  //TODO: Uncomment to find best value
	  for(unsigned int k=0;k<(100/STEP);k++){
//  for(unsigned int k=0;k<1;k++){

		  //Apply cut per plane
		  std::vector<std::vector<gate::Hit*> > planesCut(6);
		  for(unsigned int i=0; i<6;i++){
			  applyCut(planes[i], 0.01*STEP*k ,planesCut[i]);
		//	  applyCut(planes[i], 0.,planesCut[i]);
		  }
		  //std::cout << "cut: " << 0.01*STEP*k << std::endl;
		  

		  //Point Reconstruction
		  gate::Point3D reconsPoint; 
		  gate::Point3D reconsPoint2; 
		  gate::Point3D reconsPoint3; 
		  gate::Point3D reconsPoint4; 
		  gate::Point3D reconsPoint5; 

/*		  reconsPerPlane(clusters,trueVertex);  
		  reconstruction(clusters,reconsPoint);
		  reconstructionNoNorm(clusters,reconsPoint2);
		  bestPointRecons(clusters,trueVertex,reconsPoint3);
		  reconstruct2NearestPlanes(clusters, planes, reconsPoint4);
		  //reconstruct2NearestPlanesByMaxSiPM(clusters, planes, reconsPoint5);
*/
		  reconsPerPlane(planesCut,trueVertex);  
		  reconstruction(planesCut,reconsPoint);
		  reconstructionNoNorm(planesCut,reconsPoint2);
		  bestPointRecons(planesCut,trueVertex,reconsPoint3);
		  reconstruct2NearestPlanes(planesCut, planes, reconsPoint4);
//		  reconstruct2NearestPlanesByMaxSiPM(planesCut, planes, reconsPoint5);

		  //////////////
		  // Position //
		  //////////////
		  if(5*k == 60){
			  gate::Centella::instance()
				  ->hman()->fill2d(this->alabel("xPosZ"), trueVertex.z(), reconsPoint.x()-trueVertex.x());
			  gate::Centella::instance()
				  ->hman()->fill2d(this->alabel("yPosZ"), trueVertex.z(), reconsPoint.y()-trueVertex.y());
		  }
		  ///////

		  string namexNorm = "xNorm_" + gate::to_string(5*k);
		  string nameyNorm = "yNorm_" + gate::to_string(5*k);
		  string namezNorm = "zNorm_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namexNorm), reconsPoint.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameyNorm), reconsPoint.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namezNorm), reconsPoint.z() - trueVertex.z());

	/*	  std::cout << "Norm: x-x0 = " << reconsPoint.x() - trueVertex.x() << "\t y-y0 = " 
			  << reconsPoint.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint.z() - trueVertex.z() << std::endl;
		  std::cout << "NoNorm: x-x0 = " << reconsPoint2.x() - trueVertex.x() << "\t y-y0 = " 
			  << reconsPoint2.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint2.z() - trueVertex.z() << std::endl;
		  std::cout << "Best: x-x0 = " << reconsPoint3.x() - trueVertex.x() << "\t y-y0 = " 
			  << reconsPoint3.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint3.z() - trueVertex.z() << std::endl;
		  std::cout << "Best2Near: x-x0 = " << reconsPoint4.x() - trueVertex.x() << "\t y-y0 = " 
			  << reconsPoint4.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint4.z() - trueVertex.z() << std::endl;
		  std::cout << "Best2Near: x = " << reconsPoint4.x() << "\t y = " << reconsPoint4.y() << "\t z = " << reconsPoint4.z() << std::endl;*/
		  //std::cout << "Best2NearBySiPM: x-x0 = " << reconsPoint5.x() - trueVertex.x() << "\t y-y0 = " 
		//	  << reconsPoint5.y() - trueVertex.y() << "\t z-z0 = " << reconsPoint5.z() - trueVertex.z() << std::endl;
		  //std::cout << "Best2NearBySiPM: x = " << reconsPoint5.x() << "\t y = " << reconsPoint5.y() << "\t z = " << reconsPoint5.z() << std::endl;

		  string namexNoNorm = "xNoNorm_" + gate::to_string(5*k);
		  string nameyNoNorm = "yNoNorm_" + gate::to_string(5*k);
		  string namezNoNorm = "zNoNorm_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namexNoNorm), reconsPoint2.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameyNoNorm), reconsPoint2.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namezNoNorm), reconsPoint2.z() - trueVertex.z());

		  string namexbestTrue = "xbestTrue_" + gate::to_string(5*k);
		  string nameybestTrue = "ybestTrue_" + gate::to_string(5*k);
		  string namezbestTrue = "zbestTrue_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namexbestTrue), reconsPoint3.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameybestTrue), reconsPoint3.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namezbestTrue), reconsPoint3.z() - trueVertex.z());

		  string namex2Near = "x2Near_" + gate::to_string(5*k);
		  string namey2Near = "y2Near_" + gate::to_string(5*k);
		  string namez2Near = "z2Near_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namex2Near), reconsPoint4.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namey2Near), reconsPoint4.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namez2Near), reconsPoint4.z() - trueVertex.z());

		  string namex2NearBySiPM = "x2NearBySiPM_" + gate::to_string(5*k);
		  string namey2NearBySiPM = "y2NearBySiPM_" + gate::to_string(5*k);
		  string namez2NearBySiPM = "z2NearBySiPM_" + gate::to_string(5*k);
		  //std::cout << namez2NearBySiPM << std::endl;
		  //std::cout << reconsPoint5.z() - trueVertex.z() << std::endl;
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namex2NearBySiPM), reconsPoint5.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namey2NearBySiPM), reconsPoint5.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(namez2NearBySiPM), reconsPoint5.z() - trueVertex.z());

		  //printSensors(planesCut);
	//	  printSensors(clusters);
//		  checkMaxSiPMPosition(planesCut,planes,trueVertex);
	//	  projectPosition({,trueVertex);
	  }

//	  }
  }

  //Hist2d to find the cut
  hist2dHits(evt);

  //Hist2d event
 // hist2dEvent(evt);

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
//  TH1* hist = gate::Centella::instance()->hman()->operator[]("petAnalysis_Energy");
  //double maxV = hist->GetBinCenter( hist->GetMaximumBin() );
//  TF1* gauF = new TF1("gauF","gaus",0,10000);
//  hist->Fit("gauF","","e",5000,6000);
//  std::cout << "FWHM res = " << 2.35*gauF->GetParameter(2)/gauF->GetParameter(1) << std::endl;

  std::vector<double> xNormSD(100/STEP), xNoNormSD(100/STEP), xNearSD(100/STEP), xNearBySiPMSD(100/STEP),xBestSD(100/STEP);
  std::vector<double> yNormSD(100/STEP), yNoNormSD(100/STEP), yNearSD(100/STEP), yNearBySiPMSD(100/STEP),yBestSD(100/STEP);
  std::vector<double> zNormSD(100/STEP), zNoNormSD(100/STEP), zNearSD(100/STEP), zNearBySiPMSD(100/STEP),zBestSD(100/STEP);
  for(unsigned int i=0;i<(100/STEP);i++){
	  string namexNorm = "petAnalysis_xNorm_" + gate::to_string(5*i);
	  string nameyNorm = "petAnalysis_yNorm_" + gate::to_string(5*i);
	  string namezNorm = "petAnalysis_zNorm_" + gate::to_string(5*i);
	  string namexNoNorm = "petAnalysis_xNoNorm_" + gate::to_string(5*i);
	  string nameyNoNorm = "petAnalysis_yNoNorm_" + gate::to_string(5*i);
	  string namezNoNorm = "petAnalysis_zNoNorm_" + gate::to_string(5*i);
	  string namexbestTrue = "petAnalysis_xbestTrue_" + gate::to_string(5*i);
	  string nameybestTrue = "petAnalysis_ybestTrue_" + gate::to_string(5*i);
	  string namezbestTrue = "petAnalysis_zbestTrue_" + gate::to_string(5*i);
	  string namex2Near = "petAnalysis_x2Near_" + gate::to_string(5*i);
	  string namey2Near = "petAnalysis_y2Near_" + gate::to_string(5*i);
	  string namez2Near = "petAnalysis_z2Near_" + gate::to_string(5*i);
	  string namex2NearBySiPM = "petAnalysis_x2NearBySiPM_" + gate::to_string(5*i);
	  string namey2NearBySiPM = "petAnalysis_y2NearBySiPM_" + gate::to_string(5*i);
	  string namez2NearBySiPM = "petAnalysis_z2NearBySiPM_" + gate::to_string(5*i);

	  xNormSD[i] = gate::Centella::instance()->hman()->fetch(namexNorm)->GetStdDev();
	  yNormSD[i] = gate::Centella::instance()->hman()->fetch(nameyNorm)->GetStdDev();
	  zNormSD[i] = gate::Centella::instance()->hman()->fetch(namezNorm)->GetStdDev();
	  xNoNormSD[i] = gate::Centella::instance()->hman()->fetch(namexNoNorm)->GetStdDev();
	  yNoNormSD[i] = gate::Centella::instance()->hman()->fetch(nameyNoNorm)->GetStdDev();
	  zNoNormSD[i] = gate::Centella::instance()->hman()->fetch(namezNoNorm)->GetStdDev();
	  xBestSD[i] = gate::Centella::instance()->hman()->fetch(namexbestTrue)->GetStdDev();
	  yBestSD[i] = gate::Centella::instance()->hman()->fetch(nameybestTrue)->GetStdDev();
	  zBestSD[i] = gate::Centella::instance()->hman()->fetch(namezbestTrue)->GetStdDev();
	  xNearSD[i] = gate::Centella::instance()->hman()->fetch(namex2Near)->GetStdDev();
	  yNearSD[i] = gate::Centella::instance()->hman()->fetch(namey2Near)->GetStdDev();
	  zNearSD[i] = gate::Centella::instance()->hman()->fetch(namez2Near)->GetStdDev();
	  xNearBySiPMSD[i] = gate::Centella::instance()->hman()->fetch(namex2NearBySiPM)->GetStdDev();
	  yNearBySiPMSD[i] = gate::Centella::instance()->hman()->fetch(namey2NearBySiPM)->GetStdDev();
	  zNearBySiPMSD[i] = gate::Centella::instance()->hman()->fetch(namez2NearBySiPM)->GetStdDev();

/*	  std::cout << "xNorm stddev " << i << ": " << xNormSD[i] << std::endl;
	  std::cout << "yNorm stddev "<< i << ": " << yNoNormSD[i] << std::endl;
	  std::cout << "zNorm stddev "<< i << ": " << zNoNormSD[i] << std::endl;
	  std::cout << "xNoNorm stddev "<< i << ": " << xNoNormSD[i] << std::endl;
	  std::cout << "yNoNorm stddev "<< i << ": " << yNoNormSD[i]<< std::endl;
	  std::cout << "zNoNorm stddev "<< i << ": " << zNoNormSD[i]<< std::endl;
	  std::cout << "xBest stddev "<< i << ": " << xBestSD[i] << std::endl;
	  std::cout << "yBest stddev "<< i << ": " << yBestSD[i] << std::endl;
	  std::cout << "zBest stddev "<< i << ": " << zBestSD[i] << std::endl;
	  std::cout << "x2Near stddev "<< i << ": " << xNearSD[i] << std::endl;
	  std::cout << "y2Near stddev "<< i << ": " << yNearSD[i] << std::endl;
	  std::cout << "z2Near stddev " << i << ": "<< zNearSD[i] << std::endl;
	  std::cout << "x2NearBySiPM stddev "<< i << ": " << xNearBySiPMSD[i] << std::endl;
	  std::cout << "y2NearBySiPM stddev "<< i << ": " << yNearBySiPMSD[i] << std::endl;
	  std::cout << "z2NearBySiPM stddev "<< i << ": " << zNearBySiPMSD[i] << std::endl;*/
  }
  std::cout << "---------------------------------" << std::endl;
  std::cout << "xNormSD min value at " << std::min_element(xNormSD.begin(), xNormSD.end()) - xNormSD.begin() << 
	" val = " <<  xNormSD[std::min_element(xNormSD.begin(), xNormSD.end()) - xNormSD.begin()] << std::endl;
  std::cout << "yNormSD min value at " << std::min_element(yNormSD.begin(), yNormSD.end()) - yNormSD.begin() <<
	" val = " <<  yNormSD[std::min_element(yNormSD.begin(), yNormSD.end()) - yNormSD.begin()] << std::endl;
  std::cout << "zNormSD min value at " << std::min_element(zNormSD.begin(), zNormSD.end()) - zNormSD.begin() << 
	 " val = " << zNormSD[std::min_element(zNormSD.begin(), zNormSD.end()) - zNormSD.begin()] << std::endl;
  std::cout << std::endl;
  std::cout << "xNoNormSD min value at " << std::min_element(xNoNormSD.begin(), xNoNormSD.end()) - xNoNormSD.begin() << 
	 " val = " << xNoNormSD[std::min_element(xNoNormSD.begin(), xNoNormSD.end()) - xNoNormSD.begin()] << std::endl;
  std::cout << "yNoNormSD min value at " << std::min_element(yNoNormSD.begin(), yNoNormSD.end()) - yNoNormSD.begin() <<
    " val = " << yNoNormSD[std::min_element(yNoNormSD.begin(), yNoNormSD.end()) - yNoNormSD.begin()] << std::endl;
  std::cout << "zNoNormSD min value at " << std::min_element(zNoNormSD.begin(), zNoNormSD.end()) - zNoNormSD.begin() << 
	 " val = " << zNoNormSD[std::min_element(zNoNormSD.begin(), zNoNormSD.end()) - zNoNormSD.begin()] << std::endl;
  std::cout << std::endl;
  std::cout << "xBestSD min value at " << std::min_element(xBestSD.begin(), xBestSD.end()) - xBestSD.begin() << 
	  " val = " << xBestSD[std::min_element(xBestSD.begin(), xBestSD.end()) - xBestSD.begin()] << std::endl;
  std::cout << "yBestSD min value at " << std::min_element(yBestSD.begin(), yBestSD.end()) - yBestSD.begin() << 
	  " val = " << yBestSD[std::min_element(yBestSD.begin(), yBestSD.end()) - yBestSD.begin()] << std::endl;
  std::cout << "zBestSD min value at " << std::min_element(zBestSD.begin(), zBestSD.end()) - zBestSD.begin() << 
	  " val = " << zBestSD[std::min_element(zBestSD.begin(), zBestSD.end()) - zBestSD.begin()] << std::endl;
  std::cout << std::endl;
  std::cout << "xNearSD min value at " << std::min_element(xNearSD.begin(), xNearSD.end()) - xNearSD.begin() << 
	  " val = " << xNearSD[std::min_element(xNearSD.begin(), xNearSD.end()) - xNearSD.begin()] << std::endl;
  std::cout << "yNearSD min value at " << std::min_element(yNearSD.begin(), yNearSD.end()) - yNearSD.begin() << 
	  " val = " << yNearSD[std::min_element(yNearSD.begin(), yNearSD.end()) - yNearSD.begin()] <<std::endl;
  std::cout << "zNearSD min value at " << std::min_element(zNearSD.begin(), zNearSD.end()) - zNearSD.begin() << 
	  " val = " << zNearSD[std::min_element(zNearSD.begin(), zNearSD.end()) - zNearSD.begin()] << std::endl;
  std::cout << std::endl;
  std::cout << "xNearBySiPMSD min value at " << std::min_element(xNearBySiPMSD.begin(), xNearBySiPMSD.end()) - xNearBySiPMSD.begin() << 
	  " val = " << xNearBySiPMSD[std::min_element(xNearBySiPMSD.begin(), xNearBySiPMSD.end()) - xNearBySiPMSD.begin()] << std::endl;
  std::cout << "yNearBySiPMSD min value at " << std::min_element(yNearBySiPMSD.begin(), yNearBySiPMSD.end()) - yNearBySiPMSD.begin() <<
	 " val = " << yNearBySiPMSD[std::min_element(yNearBySiPMSD.begin(), yNearBySiPMSD.end()) - yNearBySiPMSD.begin()] << std::endl;
  std::cout << "zNearBySiPMSD min value at " << std::min_element(zNearBySiPMSD.begin(), zNearBySiPMSD.end()) - zNearBySiPMSD.begin() <<
	 " val = " << zNearBySiPMSD[std::min_element(zNearBySiPMSD.begin(), zNearBySiPMSD.end()) - zNearBySiPMSD.begin()] << std::endl;
  std::cout << "---------------------------------" << std::endl;

  return true;

}

void petAnalysis::reconsPerPlane(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt){
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

/*	std::cout << "x0: " << pointsRecons[0][0] << "\tx0-x = " 
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
*/
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
	int counts[600];
	memset(counts, 0, 600*sizeof(int));
	for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		int id = evt.GetMCSensHits()[i]->GetSensorID();
		counts[(id/1000)*100+(id%100)] += evt.GetMCSensHits()[i]->GetAmplitude();
	}
	for(unsigned int i=0;i<100;i++){
		//Plane 1
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),i/10,10 +10-i%10 -0.5,counts[100+i]);
		//Plane 5
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +10 -0.5,10+ i%10,counts[i+500]);
		//Plane 3
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +20 -0.5, 10+i%10,counts[i+300]);
		//Plane 4
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i/10 +30 -0.5, 10+i%10,counts[i+400]);
		//Plane 2
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),10-i%10 +10, 10 + 10 - i/10 + 10 -0.5,counts[i+200]);
		//Plane 0
		gate::Centella::instance()
			->hman()->fill2d(this->alabel(histName),i%10+10, i/10,counts[i]);
	}
	fstore("index", fetch_istore("index")+1);
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

	std::cout << "---------- Plane 0 (x,y) -----------\n";
	std::cout << "\t-45\t\t-35\t\t-25\t\t-15\t\t-05\t\t+05\t\t+15\t\t+25\t\t+35\t\t+45" << std::endl;
	for(int i=0;i<10;i++){
		std::cout << (45-i*10);
		for(int j=0;j<10;j++){
			id = i*10 + j;
			count = findSensors(planes[0],id);
			std::cout << "\tid00" << id << ": " << count;
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

void petAnalysis::reconstruct2NearestPlanes(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<gate::Hit*> > planesNoCut, gate::Point3D& pt){
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

//	std::cout << "Ordering... " << std::endl;
//	for(unsigned int i=0; i<6; i++){
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
//	}

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

//	std::cout << "Best planes: " << fstPlane << ", " << sndPlane << std::endl;
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

void petAnalysis::reconstruct2NearestPlanesByMaxSiPM(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<gate::Hit*> > planesNoCut, gate::Point3D& pt){
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

	//TODO: check position of max sipm
	//	std::cout << "Plane " <<i<< " max: x=" << sortedSiPM[i][0]->GetPosition().x() << ", y=" << sortedSiPM[i][0]->GetPosition().y() << ", z=" << sortedSiPM[i][0]->GetPosition().z() << "; id " << sortedSiPM[i][0]->GetSensorID() << "; Value: " << sortedSiPM[i][0]->GetAmplitude() << std::endl;
	}


	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);
	std::vector<std::pair<int, double> > planesOrder(6);
	for(unsigned int i=0; i<6; i++){
		planesOrder[i] = std::pair<int, double>(i,sortedSiPM[i][0]->GetAmplitude());
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << " - total: " << totalCharge(planesNoCut[i]) << std::endl;
	}
	std::sort(planesOrder.begin(), planesOrder.end(), petAnalysis::chargeOrderPlanesDesc);

//	std::cout << "Ordering... " << std::endl;
//	for(unsigned int i=0; i<6; i++){
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
//	}

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

//	std::cout << "Best planes: " << fstPlane << ", " << sndPlane << std::endl;
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


//TODO
void petAnalysis::checkMaxSiPMPosition(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<gate::Hit*> > planesNoCut, gate::Point3D& truePt){
	int nonOrthogonal[6] = {2,3,0,1,5,4};


	////////////////////////////////////////
	// SELECT 2 Planes by max SiPM charge //
	////////////////////////////////////////
	std::vector<std::vector<double> > pointsRecons(6,std::vector<double>(2));
	std::vector<std::vector<double> > errors(6,std::vector<double>(2));
	computeBarycenters(planes,pointsRecons,errors);

	//max charge
	std::vector<std::vector<gate::Hit*> >  sortedSiPM(planes);
	for(unsigned int i=0; i<6; i++){
		std::sort(sortedSiPM[i].begin(), sortedSiPM[i].end(), petAnalysis::chargeOrderSensorsDesc);

	//TODO: check position of max sipm
	//	std::cout << "Plane " <<i<< " max: x=" << sortedSiPM[i][0]->GetPosition().x() << ", y=" << sortedSiPM[i][0]->GetPosition().y() << ", z=" << sortedSiPM[i][0]->GetPosition().z() << "; id " << sortedSiPM[i][0]->GetSensorID() << "; Value: " << sortedSiPM[i][0]->GetAmplitude() << std::endl;
	}

	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);
	std::vector<std::pair<int, double> > planesOrder(6);
	for(unsigned int i=0; i<6; i++){
		planesOrder[i] = std::pair<int, double>(i,sortedSiPM[i][0]->GetAmplitude());
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << " - total: " << totalCharge(planesNoCut[i]) << std::endl;
	}
	std::sort(planesOrder.begin(), planesOrder.end(), petAnalysis::chargeOrderPlanesDesc);

//	std::cout << "Ordering... " << std::endl;
//	for(unsigned int i=0; i<6; i++){
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
//	}

	int fstPlane = planesOrder[0].first;
	int sndPlane = planesOrder[1].first;

	if(sndPlane == nonOrthogonal[planesOrder[0].first]){
		sndPlane = planesOrder[2].first;
	}

	//std::cout << "Best planes: " << fstPlane << " (max = id" << sortedSiPM[fstPlane][0]->GetSensorID()
	  // 	<< "), " << sndPlane << " (max = id" << sortedSiPM[sndPlane][0]->GetSensorID() << std::endl;
	std::vector<int> ids(6);
	projectPosition(ids,truePt);
	int countBySiPM=0;
	if(sortedSiPM[fstPlane][0]->GetSensorID() == ids[fstPlane]){
		countBySiPM++;
	}	
	if(sortedSiPM[sndPlane][0]->GetSensorID() == ids[fstPlane]){
		countBySiPM++;
	}	
//	std::cout << "countBySiPM: " << countBySiPM << std::endl;
	gate::Centella::instance()
		->hman()->fill(this->alabel("correctSiPMByMax"), countBySiPM);



	/////////////////////////////////////////
	// SELECT 2 Planes by max total charge //
	/////////////////////////////////////////
	std::vector<std::pair<int, double> > planesOrderTotal(6);
	for(unsigned int i=0; i<6; i++){
		planesOrder[i] = std::pair<int, double>(i,totalCharge(planesNoCut[i]));
	//	std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
	}
	std::sort(planesOrderTotal.begin(), planesOrderTotal.end(), petAnalysis::chargeOrderPlanesDesc);

//	std::cout << "Ordering... " << std::endl;
//	for(unsigned int i=0; i<6; i++){
//		std::cout << "Plane " << planesOrder[i].first << " - charge: " << planesOrder[i].second << std::endl;
//	}

	fstPlane = planesOrderTotal[0].first;
	sndPlane = planesOrderTotal[1].first;

	if(sndPlane == nonOrthogonal[planesOrderTotal[0].first]){
		sndPlane = planesOrderTotal[2].first;
	}

	int countByTotal=0;
	if(sortedSiPM[fstPlane][0]->GetSensorID() == ids[fstPlane]){
		countByTotal++;
	}	
	if(sortedSiPM[sndPlane][0]->GetSensorID() == ids[fstPlane]){
		countByTotal++;
	}	
//	std::cout << "countByTotal: " << countByTotal << std::endl;
	gate::Centella::instance()
		->hman()->fill(this->alabel("correctSiPMByTotal"), countByTotal);

	int correctSiPM=0;
	for(unsigned int i=0;i<6;i++){
		if(sortedSiPM[i][0]->GetSensorID() == ids[i]){
			correctSiPM++;
		}
	}
//	std::cout << "correctSiPM: " << correctSiPM << std::endl;
	gate::Centella::instance()
		->hman()->fill(this->alabel("correctSiPM"), correctSiPM);
}

void petAnalysis::projectPosition(std::vector<int>& ids, gate::Point3D& truePt){
	//int ids[6] = {-1,-1,-1,-1,-1,-1};
	int i,j;
	//std::cout << "x: " << truePt.x() << "\t y: " << truePt.y() << "\t z:" << truePt.z() << std::endl;
	//Plane 0
	j = floor((truePt.x()+50)/10);
	i = 9-floor((truePt.y()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[0] = i*10 + j;	
	//Plane 1
	j = floor((truePt.y()+50)/10);
	i = 9-floor((truePt.z()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[1] = 1000 + i + (9-j)*10;
	//Plane 2
	j = floor((truePt.x()+50)/10);
	i = 9-floor((truePt.y()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[2] = 2000 + i*10 + (9-j);
	//Plane 3
	j = floor((truePt.y()+50)/10);
	i = 9-floor((truePt.z()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[3] = 3000 + (9-i) + (9-j)*10;
	//Plane 4
	j = floor((truePt.x()+50)/10);
	i = 9-floor((truePt.z()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[4] = 4000 + (9-i) + j*10;
	//Plane 5
	j = floor((truePt.x()+50)/10);
	i = 9-floor((truePt.z()+50)/10);
	//std::cout << "i: " << i << "\t j: " << j << std::endl;
	ids[5] = 5000 + (9-i) + (9-j)*10;

/*	std::cout << "Plane 0: sensor " << ids[0] << std::endl;
	std::cout << "Plane 1: sensor " << ids[1] << std::endl;
	std::cout << "Plane 2: sensor " << ids[2] << std::endl;
	std::cout << "Plane 3: sensor " << ids[3] << std::endl;
	std::cout << "Plane 4: sensor " << ids[4] << std::endl;
	std::cout << "Plane 5: sensor " << ids[5] << std::endl;*/
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
	double ratios[100] = {2.93554, 2.81048, 2.73079, 2.69554, 2.6162, 2.59906, 2.55726, 2.48745, 2.41752, 2.41324, 2.35573, 2.31105, 2.24191, 2.19726, 2.17396, 2.13864, 2.08448, 2.06149, 2.01218, 1.95327, 1.91282, 1.86527, 1.84818, 1.80134, 1.76512, 1.72113, 1.7054, 1.666, 1.62286, 1.62326, 1.5602, 1.52193, 1.49063, 1.46494, 1.4121, 1.38277, 1.37409, 1.35741, 1.3119, 1.286, 1.26779, 1.22895, 1.20616, 1.17619, 1.14781, 1.11585, 1.10667, 1.06803, 1.05783, 1.03391, 1.00414, 0.979771, 0.962838, 0.947059, 0.931132, 0.898466, 0.879752, 0.863115, 0.853252, 0.830952, 0.819524, 0.788636, 0.762766, 0.751575, 0.746471, 0.73375, 0.708889, 0.708974, 0.676154, 0.661842, 0.652062, 0.644059, 0.63375, 0.62191, 0.605263, 0.567391, 0.554478, 0.55, 0.55, 0.54, 0.546774, 0.543617, 0.516667, 0.485294, 0.474561, 0.457576, 0.459009, 0.45, 0.45, 0.44505, 0.448592, 0.444231, 0.440625, 0.423171, 0.375926, 0.386842, 0.377692, 0.361628, 0.357353, 0.35};
	double z;
	for(unsigned int i=0;i<100;i++){
		if(ratio >= ratios[i]){
			z = -24.75 + 0.5*i;
			break;
		}
	}
	return z;
}
