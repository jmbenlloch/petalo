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

  // Energy
  gate::Centella::instance()
    ->hman()->h1(this->alabel("Energy"),"Energy SiPM " + fetch_sstore("CONF"),500,0,20000);
  gate::Centella::instance()
    ->hman()->h1(this->alabel("EnergyPhot"),"Energy Phot SiPM " + fetch_sstore("CONF"),500,0,20000);
  gate::Centella::instance()
    ->hman()->h1(this->alabel("EnergyCompt"),"Energy Compton SiPM " + fetch_sstore("CONF"),500,0,20000);

  // Best reconstruction
/*  gate::Centella::instance()
	  ->hman()->h1(this->alabel("xBest"),"xRecons-xTrue " + fetch_sstore("CONF"),100,-25,25);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("yBest"),"yRecons-yTrue " + fetch_sstore("CONF"),100,-25,25);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zBest"),"zRecons-zTrue " + fetch_sstore("CONF"),100,-25,25);*/

  //Find best cut
  for(unsigned int i=0;i<20;i++){
	  string nameX = "x_" + gate::to_string(5*i);
	  string nameY = "y_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameX),"xRecons-xTrue " + fetch_sstore("CONF"),100,-25,25);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameY),"yRecons-yTrue " + fetch_sstore("CONF"),100,-25,25);
  }

  // z Ratio
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("zRatio"),"Charge Ratio Planes 0-2 " + fetch_sstore("CONF"),100,-25,25,100,0,10);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zReconsRatio"),"zRecons-zTrue using ratio " + fetch_sstore("CONF"),100,-25,25);

  // SiPM Relative Charge
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane0"),"SiPM Relative Charge Plane 0 " + fetch_sstore("CONF"),64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane2"),"SiPM Relative Charge Plane 2 " + fetch_sstore("CONF"),64,0,64,100,0.,1.);

  //SiPMMC as function of z
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_Plane0"),"SiPM Max Charge (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_Plane2"),"SiPM Max Charge (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane0"),"SiPM Charge Corona 1 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane2"),"SiPM Charge Corona 1 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane0"),"SiPM Charge Corona 2 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane2"),"SiPM Charge Corona 2 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane0_Norm"),"SiPM Charge Corona 1 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,100,0,1000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane2_Norm"),"SiPM Charge Corona 1 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,100,0,1000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane0_Norm"),"SiPM Charge Corona 2 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,100,0,1000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane2_Norm"),"SiPM Charge Corona 2 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,100,0,1000);

  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC"),"Charge without SiPMMC (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC"),"Charge without SiPMMC (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C1"),"Charge without SiPMMC+Cluster1 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C1"),"Charge without SiPMMC+Cluster1 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C2"),"Charge without SiPMMC+Cluster2 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C2"),"Charge without SiPMMC+Cluster2 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C1_Norm"),"Charge without SiPMMC+Cluster1 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C1_Norm"),"Charge without SiPMMC+Cluster1 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C2_Norm"),"Charge without SiPMMC+Cluster2 (Plane 0) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C2_Norm"),"Charge without SiPMMC+Cluster2 (Plane 2) " + fetch_sstore("CONF"),100,-25,25,2000,0,6000);

  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0"),"Charge in Plane 0 " + fetch_sstore("CONF"),100,-25,25,2000,0,6500);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2"),"Charge in Plane 2 " + fetch_sstore("CONF"),100,-25,25,2000,0,6500);

  //Energy distribution
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xyEnergy"),"xy",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("xzEnergy"),"xz",100,-25,25,100,-25,25);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("yzEnergy"),"yz",100,-25,25,100,-25,25);

  //Parametrization
  int xmin,xmax,nBins;
  if(fetch_sstore("CONF") == "LXSC2_Z2"){
	  xmin = -10;
	  xmax = 10;
	  nBins = 10;
  }
  if(fetch_sstore("CONF") == "LXSC2_Z3"){
	  xmin = -15;
	  xmax = 15;
	  nBins = 15;
  }
  if(fetch_sstore("CONF") == "LXSC2_Z4"){
	  xmin = -20;
	  xmax = 20;
	  nBins = 20;
  }
  if(fetch_sstore("CONF") == "LXSC2_Z5"){
	  xmin = -25;
	  xmax = 25;
	  nBins = 25;
  }

  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_Plane0"),"SiPM Max Charge (Plane 0) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_Plane2"),"SiPM Max Charge (Plane 2) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_C1_Plane0"),"SiPM Charge Corona 1 (Plane 0) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_C1_Plane2"),"SiPM Charge Corona 1 (Plane 2) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_C2_Plane0"),"SiPM Charge Corona 2 (Plane 0) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Param_SiPMMC_C2_Plane2"),"SiPM Charge Corona 2 (Plane 2) " + fetch_sstore("CONF"),nBins,xmin,xmax,2000,0,6000);

  store("photoCount",0);
  store("comptCount",0);

  return true;
}


//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  //Fill energy hists
  energyPhotCompt(evt);

  // Search primary particle and its first daughter
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  //True Vertex
  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 

  //Compton count
  if(firstDaughter.GetCreatorProc() == std::string("compt") 
	//	  && firstDaughter.GetDaughters().size()==1
		  && firstDaughter.GetInitialVol() == "ACTIVE"
		  && totalCharge(evt.GetMCSensHits()) > 8000){
	  int comptCount = fetch_istore("comptCount");
	  fstore("comptCount",comptCount+1);
  }

 //Try only events with photoelectric and one vertex
//  if(firstDaughter.GetCreatorProc() == std::string("phot") 
//		  && firstDaughter.GetDaughters().size()==0){
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0 
		  && totalCharge(evt.GetMCSensHits()) > 8000 ){

	  //Store fraction of photoelectric
	  int photoCount = fetch_istore("photoCount");
	  fstore("photoCount",photoCount+1);

	  //std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

	  //Classify sensor hits per planes
	  std::vector<std::vector<gate::Hit*> > planes(6);
	  splitHitsPerPlane(evt,planes);

	  //Reconstruction
//	  reconstruction(planes,trueVertex);
	  reconstructionCorona1(planes,trueVertex);

	  //SiPMMC Charge Histograms
	  sipmmcHist(planes,trueVertex);

	  //Charge histograms
	  chargeHist2d(planes);

	  //Energy distribution
	  double energy = totalCharge(evt.GetMCSensHits());
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xyEnergy"),trueVertex.x(),trueVertex.y(),energy);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("xzEnergy"),trueVertex.x(),trueVertex.z(),energy);
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("yzEnergy"),trueVertex.y(),trueVertex.z(),energy);

	  //zRatio
	  gate::Centella::instance()
		  ->hman()->fill2d(this->alabel("zRatio"),trueVertex.z(),planeCharge(planes[0])/planeCharge(planes[2]));
	  gate::Centella::instance()
		  ->hman()->fill(this->alabel("zReconsRatio"), zReconsRatio(planeCharge(planes[0])/planeCharge(planes[2])) - trueVertex.z());

	  //Find best cut
	  for(unsigned int k=0;k<20;k++){
		  //Apply cut
		  std::vector<std::vector<gate::Hit*> > planesCut(6);
		  for(unsigned int i=0; i<6;i++){
			  applyCut(planes[i], 0.05*k,planesCut[i]);
		  }

		  gate::Point3D reconsPointBest; 
		  bestPointRecons(planesCut,trueVertex,reconsPointBest);

		  string nameX = "x_" + gate::to_string(5*k);
		  string nameY = "y_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameX), reconsPointBest.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameY), reconsPointBest.y() - trueVertex.y());
		  /*	  std::cout << "xBest: " << reconsPointBest.x() - trueVertex.x() << std::endl;
				  std::cout << "yBest: " << reconsPointBest.y() - trueVertex.y() << std::endl;*/
	  }

  }

  return true;
}

//==========================================================================
bool petAnalysis::finalize(){
	//==========================================================================

	_m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);

  gate::Run* run = &gate::Centella::instance()->getRun();
  int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  _m.message("Number of generated events in file:",nevt,gate::NORMAL);

  std::vector<double> sigmaX(20),sigmaY(20),sigmaAvg(20);
  for(unsigned int i=0;i<20;i++){
	  string nameX = "petAnalysis_x_" + gate::to_string(5*i);
	  string nameY = "petAnalysis_y_" + gate::to_string(5*i);
	  TF1* gauF = new TF1("gauF","gaus",-50,50);

	  TH1* histX = gate::Centella::instance()->hman()->operator[](nameX);
	  histX->Fit("gauF","","e",-15,15);
	  sigmaX[i] = gauF->GetParameter(2);

	  TH1* histY = gate::Centella::instance()->hman()->operator[](nameY);
	  histY->Fit("gauF","","e",-15,15);
	  sigmaY[i] = gauF->GetParameter(2);

	  sigmaAvg[i] = (sigmaX[i] + sigmaY[i])/2.0;
  }

/*  for(unsigned int i=0;i<20;i++){
	  std::cout << "x_" << i*0.05 << ": " << sigmaX[i] << std::endl;
  }
  for(unsigned int i=0;i<20;i++){
	  std::cout << "y_" << i*0.05 << ": " << sigmaY[i] << std::endl;
  }*/

  int indexBest = std::min_element(sigmaAvg.begin(), sigmaAvg.end()) - sigmaAvg.begin();
  string nameX = "petAnalysis_x_" + gate::to_string(5*indexBest);
  string nameY = "petAnalysis_y_" + gate::to_string(5*indexBest);

  gate::Centella::instance()->logman().addLog("USER","bestCut",5*indexBest);
  gate::Centella::instance()->logman().addLog("USER","photoCount",fetch_istore("photoCount"));
  gate::Centella::instance()->logman().addLog("USER","comptCount",fetch_istore("comptCount"));

/*  int nBins = gate::Centella::instance()->hman()->operator[]("petAnalysis_xBest")->GetNbinsX();
  //Copy best hist to xBest,yBest,zBest
  for(int i=0;i<nBins;i++){
	  gate::Centella::instance()->hman()->operator[]("petAnalysis_xBest")->SetBinContent(i,gate::Centella::instance()->hman()->operator[](nameX)->GetBinContent(i));
	  gate::Centella::instance()->hman()->operator[]("petAnalysis_yBest")->SetBinContent(i,gate::Centella::instance()->hman()->operator[](nameY)->GetBinContent(i));
	  gate::Centella::instance()->hman()->operator[]("petAnalysis_zBest")->SetBinContent(i,gate::Centella::instance()->hman()->operator[](nameZ)->GetBinContent(i));
  }*/

/*  std::cout << "Cut minValue (avg): " << 0.05*(std::min_element(sigmaAvg.begin(), sigmaAvg.end()) - sigmaAvg.begin()) << std::endl;
  std::cout << "Cut minValue (x): " << 0.05*(std::min_element(sigmaX.begin(), sigmaX.end()) - sigmaX.begin()) << std::endl;
  std::cout << "Cut minValue (y): " << 0.05*(std::min_element(sigmaY.begin(), sigmaY.end()) - sigmaY.begin()) << std::endl;*/

  //Try rename
  gate::Centella::instance()->hman()->operator[](nameX)->SetName("petAnalysis_xBest");
  gate::Centella::instance()->hman()->operator[](nameY)->SetName("petAnalysis_yBest");

  return true;

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

void petAnalysis::bestPointRecons(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt){
	double error=0.;
	gate::Point3D auxPt;
	//Calculate barycenter
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();
	std::vector<double> x,y,z;

	for(unsigned int i=0;i<6;i++){
		if(planes[i].size()>0){
//			std::cout << "one more plane" << i << std::endl;
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
//		std::cout << "\t x: " << x[i] << std::endl;
		if(std::abs(x[i] - truePt.x()) < error){
			pt.x(x[i]);
			error = std::abs(x[i] - truePt.x());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<y.size();i++){
//		std::cout << "\t y: " << y[i] << std::endl;
		if(std::abs(y[i] - truePt.y()) < error){
			pt.y(y[i]);
			error = std::abs(y[i] - truePt.y());
		}
	}
	error = 100000.;
	for(unsigned int i=0;i<z.size();i++){
//		std::cout << "\t z: " << z[i] << std::endl;
		if(std::abs(z[i] - truePt.z()) < error){
			pt.z(z[i]);
			error = std::abs(z[i] - truePt.z());
		}
	}

//	std::cout << "Best x: " << pt.x() << "\t y: " << pt.y() << "\t z: " << pt.z() << std::endl;
}

double petAnalysis::totalCharge(std::vector<gate::Hit*> hits){
	double total=0;
	for(unsigned int i=0;i<hits.size();i++){
		total += hits[i]->GetAmplitude();
	}
	return total;
}

double petAnalysis::planeCharge(std::vector<gate::Hit*> plane){
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

void petAnalysis::energyPhotCompt(gate::Event& evt){
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  int energy = 0;
  for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
	  energy += evt.GetMCSensHits()[i]->GetAmplitude();
  }

  gate::Centella::instance()
    ->hman()->fill(this->alabel("Energy"),energy);

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

//Must be update for each case
double petAnalysis::zReconsRatio(double ratio){
	double ratios2_z2[40] = {1.25, 1.23607, 1.2276, 1.20388, 1.20783, 1.19315, 1.18744, 1.17124, 1.15269, 1.13978, 1.11809, 1.1065, 1.09767, 1.06979, 1.06548, 1.05759, 1.05063, 1.0434, 1.04241, 1.02402, 1.00759, 0.996875, 0.967219, 0.959259, 0.952198, 0.950725, 0.947436, 0.937079, 0.924054, 0.918889};
	double ratios2_z3[60] = {1.51, 1.435, 1.426, 1.39255, 1.38953, 1.36842, 1.33793, 1.32885, 1.31818, 1.30128, 1.29118, 1.27571, 1.2475, 1.24111, 1.23571, 1.205, 1.19773, 1.18947, 1.16163, 1.15313, 1.14091, 1.13947, 1.12442, 1.10714, 1.075, 1.0575, 1.05323, 1.05, 1.05, 1.03421, 1.00313, 0.985484, 0.971739, 0.95, 0.953571, 0.95, 0.946667, 0.932143, 0.925, 0.885, 0.874138, 0.877586, 0.85, 0.85, 0.84697, 0.8375, 0.822727, 0.816667, 0.783333, 0.8, 0.782143, 0.790741, 0.75, 0.753846, 0.754167, 0.733333, 0.745238, 0.743333, 0.716667, 0.690909};
	double ratios2_z4[80] = {1.64211, 1.60769, 1.58801, 1.56278, 1.55812, 1.5251, 1.52488, 1.48744, 1.49459, 1.47769, 1.44801, 1.43564, 1.41747, 1.40481, 1.39052, 1.37195, 1.35944, 1.33555, 1.33009, 1.30196, 1.29475, 1.2742, 1.26105, 1.24628, 1.24241, 1.21957, 1.19497, 1.17833, 1.16769, 1.15207, 1.14789, 1.1372, 1.12735, 1.0956, 1.06859, 1.06318, 1.05129, 1.04478, 1.03976, 1.02669, 1.01277, 0.984711, 0.966456, 0.955797, 0.950769, 0.94916, 0.947436, 0.934043, 0.919828, 0.884146, 0.874, 0.8548, 0.851869, 0.85, 0.85, 0.844309, 0.834034, 0.814211, 0.800877, 0.791584, 0.768919, 0.752609, 0.753419, 0.75, 0.748889, 0.746429, 0.733607, 0.723469, 0.710177, 0.699425, 0.693836, 0.668391, 0.663793, 0.668033, 0.659184, 0.65, 0.655455, 0.651266, 0.641589, 0.619737};
	double ratios2_z5[100] = {1.86935, 1.85418, 1.81627, 1.78846, 1.78493, 1.73993, 1.75, 1.71221, 1.68512, 1.67117, 1.65769, 1.64313, 1.62394, 1.60115, 1.58216, 1.56194, 1.56026, 1.52786, 1.50802, 1.48767, 1.48294, 1.45776, 1.43457, 1.43051, 1.40761, 1.37722, 1.37313, 1.34653, 1.33877, 1.32772, 1.30695, 1.28504, 1.25714, 1.25221, 1.24211, 1.2244, 1.20361, 1.186, 1.17701, 1.15978, 1.15, 1.14773, 1.12761, 1.10662, 1.07695, 1.05957, 1.05, 1.05292, 1.04853, 1.03425, 1.00743, 0.988542, 0.969588, 0.954762, 0.95, 0.95, 0.94542, 0.939873, 0.892391, 0.890945, 0.870408, 0.857463, 0.854348, 0.85, 0.85, 0.842135, 0.823171, 0.807143, 0.780435, 0.766, 0.751563, 0.751408, 0.75, 0.75, 0.75, 0.746512, 0.733495, 0.736154, 0.695556, 0.680337, 0.658824, 0.674, 0.65, 0.653261, 0.65, 0.65, 0.65, 0.641358, 0.639474, 0.632278, 0.591111, 0.59625, 0.580488, 0.566129, 0.559259, 0.557692, 0.55, 0.559375, 0.553774, 0.548438};

	double*  ratios;
	int bins=0;
	double offset=0.;
	if(fetch_sstore("CONF") == "LXSC2_Z2"){
		ratios = ratios2_z2;
		bins = 40;
		offset = 15.;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z3"){
		ratios = ratios2_z3;
		bins = 60;
		offset = 10.;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z4"){
		ratios = ratios2_z4;
		bins = 80;
		offset = 5.;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z5"){
		ratios = ratios2_z5;
		bins = 100;
		offset = 0.;
	}

	double z;
	for(int i=0;i<bins;i++){
		if(ratio >= ratios[i]){
			z = -24.75 + offset + 0.5*i;
			break;
		}
	}
	return z;
}



void petAnalysis::chargeHist2d(std::vector<std::vector<gate::Hit*> > planes){
	std::vector<std::vector<gate::Hit*> >  sortedPlanes(planes);
	for(unsigned int i=0; i<6;i++){
		if(sortedPlanes[i].size() > 0){
			std::sort(sortedPlanes[i].begin(), sortedPlanes[i].end(), petAnalysis::chargeOrderSensorsDesc);
			for(unsigned int j=0; j<sortedPlanes[i].size();j++){
				string histNameRel = "SiPM_Plane" + gate::to_string(i);
				gate::Centella::instance()
					->hman()->fill2d(this->alabel(histNameRel),j, sortedPlanes[i][j]->GetAmplitude() / sortedPlanes[i][0]->GetAmplitude());
			}
		}
	}
}

void petAnalysis::sipmmcHist(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& trueVertex){
	gate::Hit* maxP0 = *std::max_element(planes[0].begin(),planes[0].end(),chargeOrderSensorsAsc);
	gate::Hit* maxP2 = *std::max_element(planes[2].begin(),planes[2].end(),chargeOrderSensorsAsc);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_Plane0"),trueVertex.z(),maxP0->GetAmplitude());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_Plane2"),trueVertex.z(),maxP2->GetAmplitude());

	util::findCluster* findCluster = new util::findCluster();
	std::vector<std::vector<gate::Hit*> > clusters1(6);
	std::vector<std::vector<gate::Hit*> > clusters2(6);

	findCluster->findCoronnaAllPlanes(planes,clusters1,1,0,0);
	findCluster->findCoronnaAllPlanes(planes,clusters2,2,0,0);

	double chargePlanes[6] = {0.,0.,0.,0.,0.,0.};
	double chargeC1[6] = {0.,0.,0.,0.,0.,0.};
	double chargeC2[6] = {0.,0.,0.,0.,0.,0.};

	for(unsigned int i=0;i<6;i++){
		chargePlanes[i] = planeCharge(planes[i]);
		chargeC1[i] = planeCharge(clusters1[i]);
		chargeC2[i] = planeCharge(clusters2[i]);
	}

	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C1_Plane0"),trueVertex.z(),chargeC1[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C1_Plane2"),trueVertex.z(),chargeC1[2]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C2_Plane0"),trueVertex.z(),chargeC2[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C2_Plane2"),trueVertex.z(),chargeC2[2]);

	//Normalize with number of SiPM
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C1_Plane0_Norm"),trueVertex.z(),chargeC1[0]/clusters1.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C1_Plane2_Norm"),trueVertex.z(),chargeC1[2]/clusters1.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C2_Plane0_Norm"),trueVertex.z(),chargeC2[0]/clusters2.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("SiPMMC_C2_Plane2_Norm"),trueVertex.z(),chargeC2[2]/clusters2.size());

	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0"),trueVertex.z(),chargePlanes[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2"),trueVertex.z(),chargePlanes[2]);

	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0_NoSiPMMC"),trueVertex.z(), chargePlanes[0] - maxP0->GetAmplitude());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2_NoSiPMMC"),trueVertex.z(), chargePlanes[2] - maxP2->GetAmplitude());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0_NoSiPMMC_C1"),trueVertex.z(), chargePlanes[0] - chargeC1[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2_NoSiPMMC_C1"),trueVertex.z(), chargePlanes[2] - chargeC1[2]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0_NoSiPMMC_C2"),trueVertex.z(), chargePlanes[0] - chargeC2[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2_NoSiPMMC_C2"),trueVertex.z(), chargePlanes[2] - chargeC2[2]);

	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0_NoSiPMMC_C1_Norm"),trueVertex.z(), chargePlanes[0] - chargeC1[0]/clusters1.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2_NoSiPMMC_C1_Norm"),trueVertex.z(), chargePlanes[2] - chargeC1[2]/clusters1.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane0_NoSiPMMC_C2_Norm"),trueVertex.z(), chargePlanes[0] - chargeC2[0]/clusters2.size());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Plane2_NoSiPMMC_C2_Norm"),trueVertex.z(), chargePlanes[2] - chargeC2[2]/clusters2.size());

	//Parametrization
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_Plane0"),trueVertex.z(),maxP0->GetAmplitude());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_Plane2"),trueVertex.z(),maxP2->GetAmplitude());
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane0"),trueVertex.z(),chargeC1[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane2"),trueVertex.z(),chargeC1[2]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_C2_Plane0"),trueVertex.z(),chargeC2[0]);
	gate::Centella::instance()
		->hman()->fill2d(this->alabel("Param_SiPMMC_C2_Plane2"),trueVertex.z(),chargeC2[2]);
}

void petAnalysis::reconstruction(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& trueVertex){
	double thresholdSiPMMC_P0_z2[2] = {723.032, 198.054};
	double thresholdSiPMMC_P2_z2[2] = {194.089, 688.701};
	double thresholdSiPMMC_C1_P0_z2[2] = {2045.35, 995.393};
	double thresholdSiPMMC_C1_P2_z2[2] = {985.271, 2008.15};
	double thresholdSiPMMC_C2_P0_z2[2] = {2370.08, 1546.52};
	double thresholdSiPMMC_C2_P2_z2[2] = {1541.09, 2347.96};

	double thresholdSiPMMC_P0_z3[3] = {683.887, 190.252, 115.318};
	double thresholdSiPMMC_P2_z3[3] = {113.575, 178.822, 645.578};
	double thresholdSiPMMC_C1_P0_z3[3] = {2017.25, 959.863, 650.798};
	double thresholdSiPMMC_C1_P2_z3[3] = {638.519, 924.236, 1938.94};
	double thresholdSiPMMC_C2_P0_z3[3] = {2362.58, 1548.82, 1172.13};
	double thresholdSiPMMC_C2_P2_z3[3] = {1168.1, 1511.04, 2325.42};

	double thresholdSiPMMC_P0_z4[4] = {624.569, 167.05, 108.444, 85.8442};
	double thresholdSiPMMC_P2_z4[4] = {85.6222, 107.816, 165.535, 610.492};
	double thresholdSiPMMC_C1_P0_z4[4] = {1981.54, 931.007, 639.871, 529.121};
	double thresholdSiPMMC_C1_P2_z4[4] = {533.26, 631.336, 923.163, 1939.23};
	double thresholdSiPMMC_C2_P0_z4[4] = {2426.35, 1582.51, 1221.23, 1084.15};
	double thresholdSiPMMC_C2_P2_z4[4] = {1090.43, 1235.6, 1568.96, 2373.5};

	double thresholdSiPMMC_P0_z5[5] = {649.707, 170.885, 111.769, 87.2472, 75.3164};
	double thresholdSiPMMC_P2_z5[5] = {75.2515, 86.6805, 111.023, 169.563, 604.085};
	double thresholdSiPMMC_C1_P0_z5[5] = {2017.07, 968.217, 669.531, 550.771, 473.303};
	double thresholdSiPMMC_C1_P2_z5[5] = {474.229, 543.411, 666.626, 956.245, 1982.23};
	double thresholdSiPMMC_C2_P0_z5[5] = {2481.14, 1621.77, 1295.75, 1138.08, 1011.12};
	double thresholdSiPMMC_C2_P2_z5[5] = {1019.72, 1124.34, 1274.16, 1627.46, 2443.23};

	double *thresholdSiPMMC_P0,*thresholdSiPMMC_P2,*thresholdSiPMMC_C1_P0,*thresholdSiPMMC_C1_P2,*thresholdSiPMMC_C2_P0,*thresholdSiPMMC_C2_P2;

	//Compute Z
	double z = zReconsRatio(planeCharge(planes[0])/planeCharge(planes[2]));
	double offset;
	if(fetch_sstore("CONF") == "LXSC2_Z2"){
		offset = 10;
		thresholdSiPMMC_P0 = thresholdSiPMMC_P0_z2;
		thresholdSiPMMC_P2 = thresholdSiPMMC_P2_z2;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z2;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z2;
		thresholdSiPMMC_C2_P0 = thresholdSiPMMC_C2_P0_z2;
		thresholdSiPMMC_C2_P2 = thresholdSiPMMC_C2_P2_z2;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z3"){
		offset = 15;
		thresholdSiPMMC_P0 = thresholdSiPMMC_P0_z3;
		thresholdSiPMMC_P2 = thresholdSiPMMC_P2_z3;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z3;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z3;
		thresholdSiPMMC_C2_P0 = thresholdSiPMMC_C2_P0_z3;
		thresholdSiPMMC_C2_P2 = thresholdSiPMMC_C2_P2_z3;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z4"){
		offset = 20;
		thresholdSiPMMC_P0 = thresholdSiPMMC_P0_z4;
		thresholdSiPMMC_P2 = thresholdSiPMMC_P2_z4;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z4;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z4;
		thresholdSiPMMC_C2_P0 = thresholdSiPMMC_C2_P0_z4;
		thresholdSiPMMC_C2_P2 = thresholdSiPMMC_C2_P2_z4;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z5"){
		offset = 25;
		thresholdSiPMMC_P0 = thresholdSiPMMC_P0_z5;
		thresholdSiPMMC_P2 = thresholdSiPMMC_P2_z5;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z5;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z5;
		thresholdSiPMMC_C2_P0 = thresholdSiPMMC_C2_P0_z5;
		thresholdSiPMMC_C2_P2 = thresholdSiPMMC_C2_P2_z5;
	}
	int index = floor((z+offset)/10);
	std::cout << "z: " << trueVertex.z() << "\t Index: " << index << std::endl;

	//Find Max
	gate::Hit* maxP0 = *std::max_element(planes[0].begin(),planes[0].end(),chargeOrderSensorsAsc);
	gate::Hit* maxP2 = *std::max_element(planes[2].begin(),planes[2].end(),chargeOrderSensorsAsc);

	util::findCluster* findCluster = new util::findCluster();
	std::vector<std::vector<gate::Hit*> > clusters1(6);
	std::vector<std::vector<gate::Hit*> > clusters2(6);

	findCluster->findCoronnaAllPlanes(planes,clusters1,1,0,0);
	findCluster->findCoronnaAllPlanes(planes,clusters2,2,0,0);

	//Planes 0 & 2
	double chargeC1[2] = {0.,0.};
	double chargeC2[2] = {0.,0.};

	for(unsigned int i=0;i<2;i++){
		chargeC1[i] = planeCharge(clusters1[i]);
		chargeC2[i] = planeCharge(clusters2[i]);
	}

	//2nd iteration
	std::vector<std::vector<gate::Hit*> > planesFiltered(6);

	std::vector<gate::Hit*> toFilter;

	int pe = 0;
	if(maxP0->GetAmplitude() >= thresholdSiPMMC_P0[index]){
		std::cout << "plane0_sipmmc\t";
		toFilter.push_back(maxP0);
		pe++;
	} else if(chargeC1[0] >= thresholdSiPMMC_C1_P0[index]){
		std::cout << "plane0_c1\t";
		toFilter = clusters1[0];
		pe++;
	} else if(chargeC2[0] >= thresholdSiPMMC_C2_P0[index]){
		std::cout << "plane0_c2\t";
		toFilter = clusters2[0];
		pe++;
	}
	filterHits(planes[0],toFilter,planesFiltered[0]);
//	std::cout << "toFilter: " << "\t size: " << toFilter.size() << std::endl;
	toFilter.clear();

	if(maxP2->GetAmplitude() >= thresholdSiPMMC_P2[index]){
		std::cout << "plane2_sipmmc\t";
		toFilter.push_back(maxP2);
		pe++;
	} else if(chargeC1[1] >= thresholdSiPMMC_C1_P2[index]){
		std::cout << "plane2_c1\t";
		toFilter = clusters1[2];
		pe++;
	} else if(chargeC2[1] >= thresholdSiPMMC_C2_P2[index]){
		std::cout << "plane2_c2\t";
		toFilter = clusters2[2];
		pe++;
	}
	filterHits(planes[2],toFilter,planesFiltered[2]);
//	std::cout << "toFilter: " << "\t size: " << toFilter.size() << std::endl;



//	for(unsigned int i=0;i<6;i++){
//		std::cout << "i: " << i << "\t size: " << planes[i].size() << std::endl;
//	}
//	std::cout << "clusters1[0]: " << "\t size: " << clusters1[0].size() << std::endl;
//	std::cout << "clusters2[0]: " << "\t size: " << clusters2[0].size() << std::endl;
//	std::cout << "clusters1[2]: " << "\t size: " << clusters1[2].size() << std::endl;
//	std::cout << "clusters2[2]: " << "\t size: " << clusters2[2].size() << std::endl;
//	for(unsigned int i=0;i<6;i++){
//		std::cout << "iFiltered: " << i << "\t size: " << planesFiltered[i].size() << std::endl;
//	}


	//2nd iteration
	gate::Hit* maxP0Filtered = *std::max_element(planesFiltered[0].begin(),planesFiltered[0].end(),chargeOrderSensorsAsc);
	gate::Hit* maxP2Filtered = *std::max_element(planesFiltered[2].begin(),planesFiltered[2].end(),chargeOrderSensorsAsc);
	std::vector<std::vector<gate::Hit*> > clusters1Filtered(6);
	std::vector<std::vector<gate::Hit*> > clusters2Filtered(6);
	findCluster->findCoronnaAllPlanes(planesFiltered,clusters1Filtered,1,0,0);
	findCluster->findCoronnaAllPlanes(planesFiltered,clusters2Filtered,2,0,0);
	
	//Planes 0 & 2
	double chargeC1Filtered[2] = {0.,0.};
	double chargeC2Filtered[2] = {0.,0.};

	for(unsigned int i=0;i<2;i++){
		chargeC1Filtered[i] = planeCharge(clusters1Filtered[i]);
		chargeC2Filtered[i] = planeCharge(clusters2Filtered[i]);
	}

	if(pe>0){
		int peFiltered = 0;
		if(maxP0Filtered->GetAmplitude() >= thresholdSiPMMC_P0[index]){
			std::cout << "filtered0_sipmmc\t";
			peFiltered++;
		} else if(chargeC1Filtered[0] >= thresholdSiPMMC_C1_P0[index]){
			std::cout << "filtered0_c1\t";
			peFiltered++;
		} else if(chargeC2Filtered[0] >= thresholdSiPMMC_C2_P0[index]){
			std::cout << "filtered0_c2\t";
			peFiltered++;
		}

		if(maxP2Filtered->GetAmplitude() >= thresholdSiPMMC_P2[index]){
			std::cout << "filtered2_sipmmc\t";
			peFiltered++;
		} else if(chargeC1Filtered[1] >= thresholdSiPMMC_C1_P2[index]){
			std::cout << "filtered2_c1\t";
			peFiltered++;
		} else if(chargeC2Filtered[1] >= thresholdSiPMMC_C2_P2[index]){
			std::cout << "filtered2_c2\t";
			peFiltered++;
		}
		std::cout << std::endl;
		std::cout << "PE-Filtered: " << peFiltered << std::endl;
	}
	std::cout << std::endl;
	std::cout << "PE: " << pe << std::endl;
}

void petAnalysis::filterHits(std::vector<gate::Hit*> plane, std::vector<gate::Hit*> hits, std::vector<gate::Hit*>& planeFiltered){
	bool found;
//	std::cout << "sizePlane: " << plane.size() << std::endl;
	for(unsigned int i=0;i<plane.size();i++){
		found = false;
		for(unsigned int j=0;j<hits.size();j++){
			if(plane[i]->GetSensorID() == hits[j]->GetSensorID()){
				found = true;
				break;
			}
		}
		if(!found){
			planeFiltered.push_back(plane[i]);
		}
	}
}


void petAnalysis::reconstructionCorona1(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& trueVertex){
/*	double thresholdSiPMMC_C1_P0_z2[2] = {2045.35, 995.393};
	double thresholdSiPMMC_C1_P2_z2[2] = {985.271, 2008.15};
	double thresholdSiPMMC_C1_P0_z3[3] = {2017.25, 959.863, 650.798};
	double thresholdSiPMMC_C1_P2_z3[3] = {638.519, 924.236, 1938.94};
	double thresholdSiPMMC_C1_P0_z4[4] = {1981.54, 931.007, 639.871, 529.121};
	double thresholdSiPMMC_C1_P2_z4[4] = {533.26, 631.336, 923.163, 1939.23};
	double thresholdSiPMMC_C1_P0_z5[5] = {2017.07, 968.217, 669.531, 550.771, 473.303};
	double thresholdSiPMMC_C1_P2_z5[5] = {474.229, 543.411, 666.626, 956.245, 1982.23};*/
	double thresholdSiPMMC_C1_P0_z2[10] = {2694.01, 2319.99, 1955.38, 1653.05, 1395.8, 1198.97, 1048.71, 959.848, 886.538, 855.157};
	double thresholdSiPMMC_C1_P2_z2[10] = {849.946, 887.385, 962.679, 1059.15, 1206.6, 1423.4, 1687.93, 1993.18, 2334.74, 2753.5};
	double thresholdSiPMMC_C1_P0_z3[15] = {2710, 2264.06, 1955.07, 1640.79, 1382.01, 1194.46, 1019.12, 908.468, 803.571, 760.052, 723.214, 659.118, 651.667, 619.382, 572.143};
	double thresholdSiPMMC_C1_P2_z3[15] = {589.167, 621.875, 636.486, 659.211, 690.287, 774.664, 832.721, 901.21, 1003.57, 1227.06, 1379.46, 1684.41, 1973.33, 2281.18, 2676.43};
	double thresholdSiPMMC_C1_P0_z4[20] = {2764.36, 2329.57, 1982.66, 1686.19, 1439.04, 1221.88, 1053.05, 946.933, 852.189, 781.672, 738.784, 691.697, 665.864, 626.528, 609.825, 595.18, 564.454, 561.752, 521.648, 504.653};
	double thresholdSiPMMC_C1_P2_z4[20] = {501.316, 520.428, 543.066, 575, 585.545, 613.824, 630, 669.709, 689.977, 746.706, 779.054, 840.619, 948.183, 1077.51, 1227.84, 1448.52, 1691.87, 1996.79, 2354.33, 2785.41};
	double thresholdSiPMMC_C1_P0_z5[25] = {476.713, 477.273, 494.618, 478.087, 481.038, 493.702, 515.164, 539.208, 563.86, 581.683, 616.713, 641.916, 656.591, 706.432, 724.551, 797.093, 824.204, 933.863, 1059.16, 1226.45, 1438.3, 1645.64, 1975.15, 2321.6, 2706.47};
	double thresholdSiPMMC_C1_P2_z5[25] = {2676.87, 2260.64, 1912.85, 1621.91, 1375.14, 1193.89, 1029.55, 909.55, 838.731, 790.594, 737.983, 700.234, 663.409, 620.539, 617.022, 591.279, 565.584, 542.057, 533.721, 497.688, 477.66, 481.831, 481.949, 468.689, 472.629};

	double *thresholdSiPMMC_C1_P0,*thresholdSiPMMC_C1_P2;

	//Compute Z
	double z = zReconsRatio(planeCharge(planes[0])/planeCharge(planes[2]));
	double offset;
	if(fetch_sstore("CONF") == "LXSC2_Z2"){
		offset = 10;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z2;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z2;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z3"){
		offset = 15;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z3;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z3;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z4"){
		offset = 20;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z4;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z4;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z5"){
		offset = 25;
		thresholdSiPMMC_C1_P0 = thresholdSiPMMC_C1_P0_z5;
		thresholdSiPMMC_C1_P2 = thresholdSiPMMC_C1_P2_z5;
	}
	//int index = floor((z+offset)/10);
	int index = floor((z+offset)/100);
	std::cout << "z: " << trueVertex.z() << "\t Index: " << index << "\tc1_p0: " << thresholdSiPMMC_C1_P0[index] << "\tc1_p2: " << thresholdSiPMMC_C1_P2[index] << std::endl;

	util::findCluster* findCluster = new util::findCluster();
	std::vector<std::vector<gate::Hit*> > clusters1(6);
	findCluster->findCoronnaAllPlanes(planes,clusters1,1,0,0);

	//Planes 0 & 2
	double chargeC1[2] = {0.,0.};
	for(unsigned int i=0;i<2;i++){
		chargeC1[i] = planeCharge(clusters1[i]);
	}

	std::cout << "c1_p0: " << chargeC1[0] << "\tc1_p2: " << chargeC1[1] << std::endl;

	//2nd iteration
	std::vector<std::vector<gate::Hit*> > planesFiltered(6);
	std::vector<gate::Hit*> toFilter;

	int pe = 0;
	if(chargeC1[0] >= thresholdSiPMMC_C1_P0[index]){
		std::cout << "plane0_c1\n";
		toFilter = clusters1[0];
		pe++;
	}
	filterHits(planes[0],toFilter,planesFiltered[0]);
//	std::cout << "toFilter: " << "\t size: " << toFilter.size() << std::endl;
	toFilter.clear();

	if(chargeC1[1] >= thresholdSiPMMC_C1_P2[index]){
		std::cout << "plane2_c1\n";
		toFilter = clusters1[2];
		pe++;
	}
	std::cout << "PE: " << pe << std::endl;
	filterHits(planes[2],toFilter,planesFiltered[2]);
	std::cout << "toFilter: " << "\t size: " << toFilter.size() << std::endl;

	for(unsigned int i=0;i<6;i++){
		std::cout << "i: " << i << "\t size: " << planes[i].size() << std::endl;
	}
	std::cout << "clusters1[0]: " << "\t size: " << clusters1[0].size() << std::endl;
	std::cout << "clusters1[2]: " << "\t size: " << clusters1[2].size() << std::endl;
	for(unsigned int i=0;i<6;i++){
		std::cout << "iFiltered: " << i << "\t size: " << planesFiltered[i].size() << std::endl;
	}


	//2nd iteration
	std::vector<std::vector<gate::Hit*> > clusters1Filtered(6);
	findCluster->findCoronnaAllPlanes(planesFiltered,clusters1Filtered,1,0,0);
	
	//Planes 0 & 2
	double chargeC1Filtered[2] = {0.,0.};
	for(unsigned int i=0;i<2;i++){
		chargeC1Filtered[i] = planeCharge(clusters1Filtered[i]);
	}
	std::cout << "c1_p0_Filtered: " << chargeC1Filtered[0] << "\tc1_p2_Filtered: " << chargeC1Filtered[1] << std::endl;

	if(pe>0){
		int peFiltered = 0;
		if(chargeC1Filtered[0] >= thresholdSiPMMC_C1_P0[index]){
			std::cout << "filtered0_c1\n";
			peFiltered++;
		}

		if(chargeC1Filtered[1] >= thresholdSiPMMC_C1_P2[index]){
			std::cout << "filtered2_c1\n";
			peFiltered++;
		}
		std::cout << "PE-Filtered: " << peFiltered << std::endl;
	}
	std::cout << std::endl;
}
