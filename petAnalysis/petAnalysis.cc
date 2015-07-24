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
	  ->hman()->h2(this->alabel("SiPMMC_Plane0"),"SiPM Max Charge (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_Plane2"),"SiPM Max Charge (Plane 2) " + fetch_sstore("CONF"),50,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane0"),"SiPM Max Charge (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C1_Plane2"),"SiPM Max Charge (Plane 2) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane0"),"SiPM Max Charge (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPMMC_C2_Plane2"),"SiPM Max Charge (Plane 2) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC"),"Charge without SiPMMC (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC"),"Charge without SiPMMC (Plane 2) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C1"),"Charge without SiPMMC+Cluster1 (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C1"),"Charge without SiPMMC+Cluster1 (Plane 2) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0_NoSiPMMC_C2"),"Charge without SiPMMC+Cluster2 (Plane 0) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2_NoSiPMMC_C2"),"Charge without SiPMMC+Cluster2 (Plane 2) " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane0"),"Charge in Plane 0 " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("Plane2"),"Charge in Plane 2 " + fetch_sstore("CONF"),100,25,-25,1000,0,3000);

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

 //Try only events with photoelectric and one vertex
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0){

	  //std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

	  //Classify sensor hits per planes
	  std::vector<std::vector<gate::Hit*> > planes(6);
	  splitHitsPerPlane(evt,planes);

	  //SiPMMC Charge Histograms
	  sipmmcHist(planes,trueVertex);

	  //Charge histograms
	  chargeHist2d(planes);

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
	double ratios2_z2[40];
	double ratios2_z3[60];
	double ratios2_z4[80];
	double ratios2_z5[100] = {1.86556, 1.85289, 1.81627, 1.78846, 1.78493, 1.73993, 1.75, 1.71221, 1.68512, 1.67117, 1.65769, 1.64313, 1.62394, 1.60115, 1.58216, 1.56194, 1.56026, 1.52786, 1.50802, 1.48767, 1.48294, 1.45776, 1.43457, 1.43051, 1.40761, 1.37722, 1.37313, 1.34653, 1.33877, 1.32772, 1.30695, 1.28504, 1.25714, 1.25221, 1.24211, 1.2244, 1.20361, 1.186, 1.17701, 1.15978, 1.15, 1.14773, 1.12761, 1.10662, 1.07695, 1.05957, 1.05, 1.05292, 1.04853, 1.03425, 1.00743, 0.988542, 0.969588, 0.954762, 0.95, 0.95, 0.94542, 0.939873, 0.892391, 0.890945, 0.870408, 0.857463, 0.854348, 0.85, 0.85, 0.842135, 0.823171, 0.807143, 0.780435, 0.766, 0.751563, 0.751408, 0.75, 0.75, 0.75, 0.746512, 0.733495, 0.736154, 0.695556, 0.680337, 0.658824, 0.674, 0.65, 0.653261, 0.65, 0.65, 0.65, 0.641358, 0.639474, 0.632278, 0.591111, 0.59625, 0.580488, 0.566129, 0.559259, 0.557692, 0.55, 0.559375, 0.553774, 0.548438};

	double*  ratios;
	if(fetch_sstore("CONF") == "LXSC2_Z2"){
		ratios = ratios2_z2;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z3"){
		ratios = ratios2_z3;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z4"){
		ratios = ratios2_z4;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z5"){
		ratios = ratios2_z5;
	}

	double z;
	for(unsigned int i=0;i<100;i++){
		if(ratio >= ratios[i]){
			z = -24.75 + 0.5*i;
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
}
