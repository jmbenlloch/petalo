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
//  std::cout << "test-init\n";
  clearData();
 // std::cout << "cleared\n";

  // Calculate the sum of the energys
  double sumE = 0.;
  for(unsigned int i=0; i<sensors.size();i++){
	  sumE += sensors[i]->GetAmplitude();
  }

  //double sumE2 = 0.;

  // loop over sensors and calculate position.
  for (unsigned int i=0; i<sensors.size();i++){
	//  std::cout << "pos: i=" << i << std::endl;
    
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

/*    // Error on E taken as sqrt(E)
	xErr_ += std::pow(sens.second*deltaPos/sumE, 2) +
		std::pow((sumE-sens.second)*xpos/std::pow(sumE,2), 2)*sens.second;
    yErr_ += std::pow(sens.second*deltaPos/sumE, 2) +
      std::pow((sumE-sens.second)*ypos/std::pow(sumE,2), 2)*sens.second;

    // For variance calculation.
    sumE2 += std::pow(sens.second, 2);  */
  }
  // Sqrt to get error.
 // xErr_ = std::sqrt(xErr_);
 // yErr_ = std::sqrt(yErr_);

//  std::cout << "test\n";
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

 // std::cout << "test-clear\n";
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
