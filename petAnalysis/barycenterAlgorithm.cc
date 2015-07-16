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
  double x1pos2=0.,x2pos2=0.;

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

	//Paola
//	x1Err_ += std::pow((sumE-sensors[i]->GetAmplitude())*x1pos/std::pow(sumE,2),2) * sensors[i]->GetAmplitude();
//	x2Err_ += std::pow((sumE-sensors[i]->GetAmplitude())*x2pos/std::pow(sumE,2),2) * sensors[i]->GetAmplitude();


    // sigma
	x1pos2 += std::pow(x1pos,2) * sensors[i]->GetAmplitude() / sumE;
	x2pos2 += std::pow(x2pos,2) * sensors[i]->GetAmplitude() / sumE;

  }
  // Sqrt to get error.
  x1Err_ = std::sqrt((std::pow(x1Pos_,2) + x1pos2)/sumE);
  x2Err_ = std::sqrt((std::pow(x2Pos_,2) + x2pos2)/sumE);

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


bool util::barycenterAlgorithm::computePositionCoronna0(const vector<gate::Hit*>& sensors, double thresholdMax, double thresholdNeighbours)
{
	bool result = false;
	int x1,x2;
	std::vector<int> idsX1,idsX2;
	vector<gate::Hit*> sensorsX1,sensorsX2;

	gate::Hit* max = *std::max_element(sensors.begin(),sensors.end(),chargeOrderSensorsDec);
	int id = max->GetSensorID();
	int planeNumber = floor(id/1000);
	if(max->GetAmplitude() > thresholdMax){
		switch(planeNumber){
			case 0:
				x1 = id%8;
				x2 = floor((id%100)/8);
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + x2*8 + (x1-1));
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + x2*8 + (x1+1));
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + (x2-1)*8 + x1);
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + (x2+1)*8 + x1);
				}
				break;
			case 1:
				x1 = floor((id%100)/8);
				x2 = id%8;
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + (x1-1)*8 + x2);
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + (x1+1)*8 + x2);
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + x1*8 + (x2-1));
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + x1*8 + (x2+1));
				}
				break;
			case 2:
				x1 = id%8;
				x2 = floor((id%100)/8);
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + x2*8 + (x1-1));
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + x2*8 + (x1+1));
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + (x2-1)*8 + x1);
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + (x2+1)*8 + x1);
				}
				break;
			case 3:
				x1 = floor((id%100)/8);
				x2 = id%8;
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + (x1-1)*8 + x2);
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + (x1+1)*8 + x2);
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + x1*8 + (x2-1));
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + x1*8 + (x2+1));
				}
				break;
			case 4:
				x1 = floor((id%100)/8);
				x2 = id%8;
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + (x1-1)*8 + x2);
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + (x1+1)*8 + x2);
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + x1 + (x2-1));
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + x1 + (x2+1));
				}
				break;
			case 5:
				x1 = floor((id%100)/8);
				x2 = id%8;
				if((x1-1) >= 0){
					idsX1.push_back(planeNumber*1000 + (x1-1)*8 + x2);
				}
				if((x1+1) <= 8){
					idsX1.push_back(planeNumber*1000 + (x1+1)*8 + x2);
				}
				if((x2-1) >= 0){
					idsX2.push_back(planeNumber*1000 + x1 + (x2-1));
				}   
				if((x2+1) <= 8){
					idsX2.push_back(planeNumber*1000 + x1 + (x2+1));
				}
				break;
		}


		if(idsX1.size() > 0 && idsX2.size() > 0){
			result = true;
			//Rethink this loops. There ir probably a better way
			for(unsigned int j=0;j<idsX1.size();j++){
				for(unsigned int i=0;i<sensors.size();i++){
					if(sensors[i]->GetSensorID() == idsX1[j]){
						if(sensors[i]->GetAmplitude() > thresholdNeighbours){
							sensorsX1.push_back(sensors[i]);
						}
						break;
					}
				}
			}
			for(unsigned int j=0;j<idsX2.size();j++){
				for(unsigned int i=0;i<sensors.size();i++){
					if(sensors[i]->GetSensorID() == idsX2[j]){
						if(sensors[i]->GetAmplitude() > thresholdNeighbours){
							sensorsX2.push_back(sensors[i]);
						}
						break;
					}
				}
			}

			// Calculate the sum of the energys
			double sumE_x1 = 0., sumE_x2 = 0.;
			for(unsigned int i=0; i<sensorsX1.size();i++){
				sumE_x1 += sensorsX1[i]->GetAmplitude();
			}
			for(unsigned int i=0; i<sensorsX2.size();i++){
				sumE_x2 += sensorsX2[i]->GetAmplitude();
			}

			double x1pos2=0.,x2pos2=0.;
			// loop over sensors for X1.
			for (unsigned int i=0; i<sensorsX1.size();i++){
				double x1pos;
				switch(plane_){
					case 0:
						x1pos = sensors[i]->GetPosition().x();
						break;
					case 1:
						x1pos = sensors[i]->GetPosition().x();
						break;
					case 2:
						x1pos = sensors[i]->GetPosition().y();
						break;
				}
				x1pos = sensorsX1[i]->GetPosition().x();
				x1Pos_ += x1pos * sensorsX1[i]->GetAmplitude() / sumE_x1;
				// sigma
				x1pos2 += std::pow(x1pos,2) * sensorsX1[i]->GetAmplitude() / sumE_x1;
			}
			// Sqrt to get error.
			x1Err_ = std::sqrt((std::pow(x1Pos_,2) + x1pos2)/sumE_x1);

			// loop over sensors for X2
			for (unsigned int i=0; i<sensorsX2.size();i++){
				double x2pos;
				switch(plane_){
					case 0:
						x2pos = sensors[i]->GetPosition().y();
						break;
					case 1:
						x2pos = sensors[i]->GetPosition().z();
						break;
					case 2:
						x2pos = sensors[i]->GetPosition().z();
						break;
				}
				x2pos = sensorsX2[i]->GetPosition().y();
				x2Pos_ += x2pos * sensorsX2[i]->GetAmplitude() / sumE_x2;
				// sigma
				x2pos2 += std::pow(x2pos,2) * sensorsX2[i]->GetAmplitude() / sumE_x2;
			}
			// Sqrt to get error.
			x2Err_ = std::sqrt((std::pow(x2Pos_,2) + x2pos2)/sumE_x2);
		}
	}

	return result;
}

bool util::barycenterAlgorithm::chargeOrderSensorsDec(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() < s2->GetAmplitude());
}
