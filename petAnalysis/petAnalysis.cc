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
	  string nameZ = "z_" + gate::to_string(5*i);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameX),"xRecons-xTrue " + fetch_sstore("CONF"),100,-25,25);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameY),"yRecons-yTrue " + fetch_sstore("CONF"),100,-25,25);
	  gate::Centella::instance()
		  ->hman()->h1(this->alabel(nameZ),"zRecons-zTrue " + fetch_sstore("CONF"),100,-25,25);
  }

  // z Ratio
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("zRatio"),"Charge Ratio Planes 0-2 " + fetch_sstore("CONF"),100,-25,25,100,0,10);
  gate::Centella::instance()
	  ->hman()->h1(this->alabel("zReconsRatio"),"zRecons-zTrue using ratio " + fetch_sstore("CONF"),100,-25,25);

  // SiPM Relative Charge
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane0"),"SiPM Relative Charge Plane 0",64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane1"),"SiPM Relative Charge Plane 1",64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane2"),"SiPM Relative Charge Plane 2",64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane3"),"SiPM Relative Charge Plane 3",64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane4"),"SiPM Relative Charge Plane 4",64,0,64,100,0.,1.);
  gate::Centella::instance()
	  ->hman()->h2(this->alabel("SiPM_Plane5"),"SiPM Relative Charge Plane 5",64,0,64,100,0.,1.);

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
		  string nameZ = "z_" + gate::to_string(5*k);
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameX), reconsPointBest.x() - trueVertex.x());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameY), reconsPointBest.y() - trueVertex.y());
		  gate::Centella::instance()
			  ->hman()->fill(this->alabel(nameZ), reconsPointBest.z() - trueVertex.z());
		  /*	  std::cout << "xBest: " << reconsPointBest.x() - trueVertex.x() << std::endl;
				  std::cout << "yBest: " << reconsPointBest.y() - trueVertex.y() << std::endl;
				  std::cout << "zBest: " << reconsPointBest.z() - trueVertex.z() << std::endl;*/
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
  string nameZ = "petAnalysis_z_" + gate::to_string(5*indexBest);

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
  gate::Centella::instance()->hman()->operator[](nameZ)->SetName("petAnalysis_zBest");

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
	double ratios2_36[100];
	double ratios2_49[100] = {1.82903, 1.81367, 1.77798, 1.76424, 1.75, 1.72834, 1.70656, 1.68653, 1.66872, 1.65187, 1.64479, 1.61789, 1.6089, 1.57614, 1.57154, 1.55254, 1.52843, 1.52126, 1.47469, 1.49109, 1.44322, 1.4481, 1.43063, 1.41899, 1.40088, 1.37727, 1.34926, 1.34274, 1.33925, 1.31358, 1.28099, 1.27061, 1.25397, 1.24884, 1.23367, 1.22238, 1.20472, 1.18231, 1.16793, 1.15245, 1.14786, 1.13402, 1.11111, 1.09937, 1.07171, 1.0644, 1.05327, 1.04524, 1.04123, 1.04552, 1.01304, 0.996296, 0.972523, 0.959272, 0.951, 0.938095, 0.940991, 0.935065, 0.894578, 0.892391, 0.856542, 0.859709, 0.857778, 0.850862, 0.848667, 0.841089, 0.839189, 0.811446, 0.799451, 0.773611, 0.767857, 0.760924, 0.75, 0.75, 0.75, 0.75, 0.736842, 0.723611, 0.708163, 0.686585, 0.685294, 0.667778, 0.65, 0.656667, 0.65, 0.65, 0.65, 0.64703, 0.637234, 0.640244, 0.619231, 0.596296, 0.584043, 0.601562, 0.594, 0.578169, 0.564286, 0.556579, 0.548333, 0.539899};
	double ratios2_64[100] = {1.86556, 1.85289, 1.81627, 1.78846, 1.78493, 1.73993, 1.75, 1.71221, 1.68512, 1.67117, 1.65769, 1.64313, 1.62394, 1.60115, 1.58216, 1.56194, 1.56026, 1.52786, 1.50802, 1.48767, 1.48294, 1.45776, 1.43457, 1.43051, 1.40761, 1.37722, 1.37313, 1.34653, 1.33877, 1.32772, 1.30695, 1.28504, 1.25714, 1.25221, 1.24211, 1.2244, 1.20361, 1.186, 1.17701, 1.15978, 1.15, 1.14773, 1.12761, 1.10662, 1.07695, 1.05957, 1.05, 1.05292, 1.04853, 1.03425, 1.00743, 0.988542, 0.969588, 0.954762, 0.95, 0.95, 0.94542, 0.939873, 0.892391, 0.890945, 0.870408, 0.857463, 0.854348, 0.85, 0.85, 0.842135, 0.823171, 0.807143, 0.780435, 0.766, 0.751563, 0.751408, 0.75, 0.75, 0.75, 0.746512, 0.733495, 0.736154, 0.695556, 0.680337, 0.658824, 0.674, 0.65, 0.653261, 0.65, 0.65, 0.65, 0.641358, 0.639474, 0.632278, 0.591111, 0.59625, 0.580488, 0.566129, 0.559259, 0.557692, 0.55, 0.559375, 0.553774, 0.548438};
	double ratios6_36[100];
	double ratios6_49[100] = {2.93554, 2.81048, 2.73079, 2.69554, 2.6162, 2.59906, 2.55726, 2.48745, 2.41752, 2.41324, 2.35573, 2.31105, 2.24191, 2.19726, 2.17396, 2.13864, 2.08448, 2.06149, 2.01218, 1.95327, 1.91282, 1.86527, 1.84818, 1.80134, 1.76512, 1.72113, 1.7054, 1.666, 1.62286, 1.62326, 1.5602, 1.52193, 1.49063, 1.46494, 1.4121, 1.38277, 1.37409, 1.35741, 1.3119, 1.286, 1.26779, 1.22895, 1.20616, 1.17619, 1.14781, 1.11585, 1.10667, 1.06803, 1.05783, 1.03391, 1.00414, 0.979771, 0.962838, 0.947059, 0.931132, 0.898466, 0.879752, 0.863115, 0.853252, 0.830952, 0.819524, 0.788636, 0.762766, 0.751575, 0.746471, 0.73375, 0.708889, 0.708974, 0.676154, 0.661842, 0.652062, 0.644059, 0.63375, 0.62191, 0.605263, 0.567391, 0.554478, 0.55, 0.55, 0.54, 0.546774, 0.543617, 0.516667, 0.485294, 0.474561, 0.457576, 0.459009, 0.45, 0.45, 0.44505, 0.448592, 0.444231, 0.440625, 0.423171, 0.375926, 0.386842, 0.377692, 0.361628, 0.357353, 0.35};
	double ratios6_64[100] = {2.73205, 2.73906, 2.66909, 2.61167, 2.52151, 2.49822, 2.44918, 2.417, 2.35749, 2.32905, 2.26496, 2.23821, 2.19439, 2.12249, 2.10054, 2.06435, 2.01884, 1.9875, 1.94561, 1.90444, 1.87033, 1.83465, 1.79233, 1.76525, 1.7309, 1.68712, 1.65806, 1.6348, 1.59805, 1.5388, 1.52833, 1.4927, 1.45451, 1.43733, 1.40809, 1.36951, 1.34272, 1.32054, 1.29409, 1.25984, 1.25199, 1.21908, 1.19365, 1.15915, 1.13903, 1.12939, 1.1031, 1.059, 1.05392, 1.03404, 1.00917, 0.992857, 0.967647, 0.930952, 0.931102, 0.915556, 0.907547, 0.859375, 0.851111, 0.840588, 0.829592, 0.80619, 0.775581, 0.760185, 0.748851, 0.747826, 0.715487, 0.711039, 0.687398, 0.666667, 0.657407, 0.654545, 0.642391, 0.646629, 0.627632, 0.604348, 0.577869, 0.565663, 0.55396, 0.55, 0.55, 0.543333, 0.547297, 0.532716, 0.502747, 0.467978, 0.462766, 0.459211, 0.453704, 0.452667, 0.45, 0.45, 0.44726, 0.441803, 0.435, 0.43, 0.371154, 0.389726, 0.367143, 0.358};


	double*  ratios;
	ratios = ratios6_64; //to avoid seg fault on lxsc4
	if(fetch_sstore("CONF") == "LXSC2_36"){
		ratios = ratios2_36;
	}
	if(fetch_sstore("CONF") == "LXSC2_49"){
		ratios = ratios2_49;
	}
	if(fetch_sstore("CONF") == "LXSC2_64"){
		ratios = ratios2_64;
	}
	if(fetch_sstore("CONF") == "LXSC6_36"){
		ratios = ratios6_36;
	}
	if(fetch_sstore("CONF") == "LXSC6_49"){
		ratios = ratios6_49;
	}
	if(fetch_sstore("CONF") == "LXSC6_64"){
		ratios = ratios6_64;
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

