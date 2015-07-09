#include <cmath>
#include "findCluster.hh"

util::findCluster::findCluster(){
}

bool
util::findCluster::findCoronnaAllPlanes(const std::vector<std::vector<gate::Hit*> >& planes, std::vector<std::vector<gate::Hit*> >& clusters){
	for(unsigned int i=0;i<planes.size();i++){
		std::cout << "------ Plane " << i << "-------\n";
		findCoronna(planes[i],clusters[i]);
	}
	return true;
}

bool
util::findCluster::findCoronna(const std::vector<gate::Hit*>& plane, std::vector<gate::Hit*>& cluster){
	// Find max value in the plane
	// Select all SiPM around max (about a threshold)
	// Delete all selected SiPM from list and store first cluster
	// Loop until no more SiPM above threshold
	
	gate::Hit* max = *std::max_element(plane.begin(),plane.end(),chargeOrderSensorsDec);
	int id = max->GetSensorID();

	std::cout << "Max sensor: " << id << std::endl;

	std::vector<int> idsFirstRing;
	idsFirstRing.push_back(id);

	if((id%1000) == 0){
		idsFirstRing.push_back(id+1);
		idsFirstRing.push_back(id+11);
		idsFirstRing.push_back(id+10);
	}
	if((id%1000) == 9){
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id+9);
		idsFirstRing.push_back(id+10);
	}
	if((id%1000) == 90){
		idsFirstRing.push_back(id-10);
		idsFirstRing.push_back(id-9);
		idsFirstRing.push_back(id+1);
	}
	if((id%1000) == 99){
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id-11);
		idsFirstRing.push_back(id-10);
	}
	if((id%1000) > 0 && (id%1000) < 9){
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id+1);
		idsFirstRing.push_back(id+9);
		idsFirstRing.push_back(id+10);
		idsFirstRing.push_back(id+11);
	}
	if((id%1000) > 90 && (id%1000) < 99){
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id+1);
		idsFirstRing.push_back(id-9);
		idsFirstRing.push_back(id-10);
		idsFirstRing.push_back(id-11);
	}
	if((id%10) == 0 && (id%1000) > 0 && (id%1000) < 90){
		idsFirstRing.push_back(id-10);
		idsFirstRing.push_back(id+10);
		idsFirstRing.push_back(id-9);
		idsFirstRing.push_back(id+1);
		idsFirstRing.push_back(id+11);
	}
	if((id%10) == 9 && (id%1000) > 9 && (id%1000) < 99){
		idsFirstRing.push_back(id-10);
		idsFirstRing.push_back(id+10);
		idsFirstRing.push_back(id-11);
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id+9);
	}
	if((id%10) != 9 && (id%10) != 0 && floor((id%1000)/10) != 0 && floor((id%1000)/10) != 9){
		idsFirstRing.push_back(id-10);
		idsFirstRing.push_back(id+10);
		idsFirstRing.push_back(id-11);
		idsFirstRing.push_back(id+11);
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id-1);
		idsFirstRing.push_back(id-9);
		idsFirstRing.push_back(id+9);
	}

	std::cout << "First ring: ";
	for(unsigned int i=0;i<idsFirstRing.size();i++){
		std::cout << "\t" << idsFirstRing[i];
	}
	std::cout << std::endl;

	for(unsigned int j=0;j<idsFirstRing.size();j++){
		for(unsigned int i=0;i<plane.size();i++){
			if(plane[i]->GetSensorID() == idsFirstRing[j]){
				cluster.push_back(plane[i]);
				break;
			}
		}
	}

	std::cout << "Selected: ";
	for(unsigned int i=0;i<cluster.size();i++){
		std::cout << "\t" << cluster[i]->GetSensorID();
	}
	std::cout << std::endl;

	return true;
}

bool util::findCluster::chargeOrderSensorsDec(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() < s2->GetAmplitude());
}

