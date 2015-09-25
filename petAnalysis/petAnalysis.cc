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

  TFile *file = new TFile((fetch_sstore("CONF") + "_ntuple.root").c_str(), "RECREATE", "An Example ROOT file");
  setFile(file);
  TTree *tree = new TTree("petalo","petalo");
  setTree(tree);
  getTree()->Branch("entryPlane",getEntryPlane(),"entryPlane[64]/D");
 
  return true;
}


//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  //Fill energy hists
//  energyPhotCompt(evt);

  // Search primary particle and its first daughter
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  //True Vertex
  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 

  //Classify sensor hits per planes
  std::vector<std::vector<gate::Hit*> > planes(6);
  splitHitsPerPlane(evt,planes);

  float phot = 0;
  float compt = 0;
  if(firstDaughter.GetCreatorProc() == std::string("compt")){
	  compt = 1;
  }
  if(firstDaughter.GetCreatorProc() == std::string("phot")){
	  phot = 1;
  }

  /////////////////////////////////////
  // Preproc data for neural network //
  /////////////////////////////////////
  double entryPlane[64];
  for(int i=0; i<64; i++){
	  entryPlane[i] = findSensors(planes[0],i);
//	  std::cout << entryPlane[i] << ", " << std::endl;
  }
  std::vector<double> sipm_plane0(entryPlane,entryPlane+64);


  for(int i=0; i<64; i++){
	  getEntryPlane()[i] = findSensors(planes[0],i);
  }
  getTree()->Fill();

  this->fstore("plane0",sipm_plane0);
  this->fstore("x",trueVertex.x());
  this->fstore("y",trueVertex.y());
  this->fstore("z",trueVertex.z());
  this->fstore("eventID",evt.GetEventID());


  //Find second daughter
/*  std::vector<const gate::MCParticle*> daughters(primary.GetDaughters());
  std::sort(daughters.begin(), daughters.end(), petAnalysis::timeOrderParticles);
  std::cout << "Event number:" << evt.GetEventID() << std::endl; 
  for(unsigned int i=0;i<daughters.size();i++){
	  std::cout << "\t Vol: " <<daughters[i]->GetInitialVol();
	  std::cout << "\t t: " << daughters[i]->GetInitialVtx4D().GetT();
	  std::cout << "\t" << daughters[i]->GetCreatorProc() << "/" << daughters[i]->GetPDG();
	  std::cout << "\t\t daughters:" << daughters[i]->GetDaughters().size();
	  std::cout << std::endl;

	  for(unsigned int j=0;j<daughters[i]->GetDaughters().size();j++){
		  std::cout << "\t\t Vol: " << daughters[i]->GetDaughters()[j]->GetInitialVol();
		  std::cout << "\t\t t: " << daughters[i]->GetDaughters()[j]->GetInitialVtx4D().GetT();
		  std::cout << "\t\t" << daughters[i]->GetDaughters()[j]->GetCreatorProc() << "/" << daughters[i]->GetDaughters()[j]->GetPDG();
		  std::cout << "\t\t daughters:" << daughters[i]->GetDaughters()[j]->GetDaughters().size();
		  std::cout << std::endl;
	  }
  }*/
  
  //Energy MCParticles
  double energy = 0.;
  for(unsigned int i=0;i<evt.GetMCParticles().size();i++){
	  for(unsigned int j=0;j<evt.GetMCParticles()[i]->GetTracks().size();j++){
		  energy += evt.GetMCParticles()[i]->GetTracks()[j]->GetHitsEnergy();
	  }
  }
  //std::cout << "Energy all MCParticles: " << energy << std::endl;

  std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

  return true;
}

//==========================================================================
bool petAnalysis::finalize(){
	//==========================================================================

	_m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);

  gate::Run* run = &gate::Centella::instance()->getRun();
  int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  _m.message("Number of generated events in file:",nevt,gate::NORMAL);

  getFile()->cd();
  getTree()->Write();
  gDirectory->pwd();

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

//Must be update for each case
double petAnalysis::zReconsRatio(double ratio){
	double ratios2_z4[80] = {1.64211, 1.60769, 1.58801, 1.56278, 1.55812, 1.5251, 1.52488, 1.48744, 1.49459, 1.47769, 1.44801, 1.43564, 1.41747, 1.40481, 1.39052, 1.37195, 1.35944, 1.33555, 1.33009, 1.30196, 1.29475, 1.2742, 1.26105, 1.24628, 1.24241, 1.21957, 1.19497, 1.17833, 1.16769, 1.15207, 1.14789, 1.1372, 1.12735, 1.0956, 1.06859, 1.06318, 1.05129, 1.04478, 1.03976, 1.02669, 1.01277, 0.984711, 0.966456, 0.955797, 0.950769, 0.94916, 0.947436, 0.934043, 0.919828, 0.884146, 0.874, 0.8548, 0.851869, 0.85, 0.85, 0.844309, 0.834034, 0.814211, 0.800877, 0.791584, 0.768919, 0.752609, 0.753419, 0.75, 0.748889, 0.746429, 0.733607, 0.723469, 0.710177, 0.699425, 0.693836, 0.668391, 0.663793, 0.668033, 0.659184, 0.65, 0.655455, 0.651266, 0.641589, 0.619737};
	double ratios2_z5[100] = {1.86935, 1.85418, 1.81627, 1.78846, 1.78493, 1.73993, 1.75, 1.71221, 1.68512, 1.67117, 1.65769, 1.64313, 1.62394, 1.60115, 1.58216, 1.56194, 1.56026, 1.52786, 1.50802, 1.48767, 1.48294, 1.45776, 1.43457, 1.43051, 1.40761, 1.37722, 1.37313, 1.34653, 1.33877, 1.32772, 1.30695, 1.28504, 1.25714, 1.25221, 1.24211, 1.2244, 1.20361, 1.186, 1.17701, 1.15978, 1.15, 1.14773, 1.12761, 1.10662, 1.07695, 1.05957, 1.05, 1.05292, 1.04853, 1.03425, 1.00743, 0.988542, 0.969588, 0.954762, 0.95, 0.95, 0.94542, 0.939873, 0.892391, 0.890945, 0.870408, 0.857463, 0.854348, 0.85, 0.85, 0.842135, 0.823171, 0.807143, 0.780435, 0.766, 0.751563, 0.751408, 0.75, 0.75, 0.75, 0.746512, 0.733495, 0.736154, 0.695556, 0.680337, 0.658824, 0.674, 0.65, 0.653261, 0.65, 0.65, 0.65, 0.641358, 0.639474, 0.632278, 0.591111, 0.59625, 0.580488, 0.566129, 0.559259, 0.557692, 0.55, 0.559375, 0.553774, 0.548438};

	double*  ratios;
	int bins=0;
	double offset=0.;
	if(fetch_sstore("CONF") == "LXSC2_Z4_64"){
		ratios = ratios2_z4;
		bins = 80;
		offset = 5.;
	}
	if(fetch_sstore("CONF") == "LXSC2_Z5_64"){
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


void petAnalysis::sipmmcHist(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& trueVertex, int phot){
	util::findCluster* findCluster = new util::findCluster();
	std::vector<std::vector<gate::Hit*> > clusters1(6);
	std::vector<std::vector<gate::Hit*> > clusters2(6);

	findCluster->findCoronnaAllPlanes(planes,clusters1,1,0,0);
	findCluster->findCoronnaAllPlanes(planes,clusters2,2,0,0);

	double chargeC1[6] = {0.,0.,0.,0.,0.,0.};
	for(unsigned int i=0;i<6;i++){
		chargeC1[i] = planeCharge(clusters1[i]);
	}

	//Parametrization
	if(phot == 1){
		gate::Centella::instance()
			->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane0"),trueVertex.z(),chargeC1[0]);
		gate::Centella::instance()
			->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane2"),trueVertex.z(),chargeC1[2]);
	}else{
		gate::Centella::instance()
			->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane0_Compt"),trueVertex.z(),chargeC1[0]);
		gate::Centella::instance()
			->hman()->fill2d(this->alabel("Param_SiPMMC_C1_Plane2_Compt"),trueVertex.z(),chargeC1[2]);
	}
}

double petAnalysis::findSensors(std::vector<gate::Hit*>& plane, int id){
	double amplitude = 0;
	for(unsigned int i=0;i<plane.size();i++){
		if(plane[i]->GetSensorID() == id){
			amplitude = plane[i]->GetAmplitude();
	//		std::cout << plane[i]->GetSensorID() << ": " << amplitude << std::endl;
			break;
		}
	}
	return amplitude;
}
