#include <cmath>
#include "barycenterAlgorithm.hh"

util::barycenterAlgorithm::barycenterAlgorithm(){
  x1Pos_ = 0.;
  x2Pos_ = 0.;
  x1Err_ = 0.;
  x2Err_ = 0.;
  plane_ = 0;
}

bool
util::barycenterAlgorithm::computePosition(const vector<gate::Hit*>& sensors)
{
  clearData();

  // Calculate the sum of the energys
  double sumE = 0.;
  for(unsigned int i=0; i<sensors.size();i++){
	  sumE += sensors[i]->GetAmplitude();
  }

  //double sumE2 = 0.;
  double x1pos2=0.,x2pos2=0.;

  // loop over sensors and calculate position.
  for (unsigned int i=0; i<sensors.size();i++){
    
	double x1pos,x2pos;

    // Get the sensor position
	switch(plane_){
		case 0:
			x1pos = sensors[i]->GetPosition().x();
			x2pos = sensors[i]->GetPosition().y();
			break;
		case 1:
			x1pos = sensors[i]->GetPosition().x();
			x2pos = sensors[i]->GetPosition().z();
			break;
		case 2:
			x1pos = sensors[i]->GetPosition().y();
			x2pos = sensors[i]->GetPosition().z();
			break;
	}

	x1Pos_ += x1pos * sensors[i]->GetAmplitude() / sumE;
    x2Pos_ += x2pos * sensors[i]->GetAmplitude() / sumE;

    // sigma
	x1pos2 += std::pow(x1pos,2) * sensors[i]->GetAmplitude() / sumE;
	x2pos2 += std::pow(x2pos,2) * sensors[i]->GetAmplitude() / sumE;

  }
  // Sqrt to get error.
  x1Err_ = std::sqrt((std::pow(x1Pos_,2) + x1pos2)/sumE);
  x2Err_ = std::sqrt((std::pow(x2Pos_,2) + x2pos2)/sumE);

  return true;
}

void
util::barycenterAlgorithm::getX1withErr(std::pair<double,double>& xErr) const
{
  xErr.first = this->x1Pos_;
  xErr.second = this->x1Err_;
}

void
util::barycenterAlgorithm::getX2withErr(std::pair<double,double>& yErr) const
{
  yErr.first = this->x2Pos_;
  yErr.second = this->x2Err_;
}

void
util::barycenterAlgorithm::clearData()
{
  // Reset values.
  x1Pos_ = 0.;
  x2Pos_ = 0.;
  x1Err_ = 0.;
  x2Err_ = 0.;
}

void util::barycenterAlgorithm::setPlane(std::string planeDirection){
	if (planeDirection == "xy"){
		plane_ = 0;
	}else if(planeDirection == "xz"){
		plane_ = 1;
	}else if(planeDirection == "yz"){
		plane_ = 2;
	}
}

