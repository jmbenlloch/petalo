#include <cmath>
#include "findCluster.hh"

util::findCluster::findCluster(){
}

bool
util::findCluster::findCoronnaAllPlanes(const std::vector<std::vector<gate::Hit*> >& planes, std::vector<std::vector<gate::Hit*> >& clusters, int rings){
	if(rings==1){
		for(unsigned int i=0;i<planes.size();i++){
	//		std::cout << "------ Plane " << i << "-------\n";
			findCoronna(planes[i],clusters[i]);
		}
	}
	if(rings==2){
		for(unsigned int i=0;i<planes.size();i++){
//			std::cout << "------ Plane " << i << "-------\n";
			findCoronna2Rings(planes[i],clusters[i]);
		}
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
	int planeNumber = floor(id/1000);
	int row = floor((id%100)/10);
	int col = id%10;

//	std::cout << "Max sensor: " << id << std::endl;
//	std::cout << "Plane: " << planeNumber << "\t row " << row << "\t col " << col << std::endl;

	std::vector<int> idsFirstRing;

	for(int i = (row-1); i <= (row+1); i++){
		if(i<0 || i>9){
			continue;
		}
		for(int j = (col-1); j <= (col+1); j++){
			if(j<0 || j>9){
				continue;
			}
			idsFirstRing.push_back(planeNumber*1000 + i*10 + j);
		}
	}

/*	std::cout << "First ring: ";
	for(unsigned int i=0;i<idsFirstRing.size();i++){
		std::cout << "\t" << idsFirstRing[i];
	}
	std::cout << std::endl;
*/
	for(unsigned int j=0;j<idsFirstRing.size();j++){
		for(unsigned int i=0;i<plane.size();i++){
			if(plane[i]->GetSensorID() == idsFirstRing[j]){
				cluster.push_back(plane[i]);
				break;
			}
		}
	}

/*	std::cout << "Selected: ";
	for(unsigned int i=0;i<cluster.size();i++){
		std::cout << "\t" << cluster[i]->GetSensorID();
	}
	std::cout << std::endl;
*/
	return true;
}

bool util::findCluster::chargeOrderSensorsDec(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() < s2->GetAmplitude());
}

bool
util::findCluster::findCoronna2Rings(const std::vector<gate::Hit*>& plane, std::vector<gate::Hit*>& cluster){
	// Find max value in the plane
	// Select all SiPM around max (about a threshold)
	// Delete all selected SiPM from list and store first cluster
	// Loop until no more SiPM above threshold
	
	gate::Hit* max = *std::max_element(plane.begin(),plane.end(),chargeOrderSensorsDec);
	int id = max->GetSensorID();
	int planeNumber = floor(id/1000);
	int row = floor((id%100)/10);
	int col = id%10;

//	std::cout << "Max sensor: " << id << std::endl;
//	std::cout << "Plane: " << planeNumber << "\t row " << row << "\t col " << col << std::endl;

	std::vector<int> idsFirstRing;
	//idsFirstRing.push_back(id);

	for(int i = (row-2); i <= (row+2); i++){
		if(i<0 || i>9){
			continue;
		}
		for(int j = (col-2); j <= (col+2); j++){
			if(j<0 || j>9){
				continue;
			}
			idsFirstRing.push_back(planeNumber*1000 + i*10 + j);
		}
	}

/*	std::cout << "First ring: ";
	for(unsigned int i=0;i<idsFirstRing.size();i++){
		std::cout << "\t" << idsFirstRing[i];
	}
	std::cout << std::endl;
*/
	for(unsigned int j=0;j<idsFirstRing.size();j++){
		for(unsigned int i=0;i<plane.size();i++){
			if(plane[i]->GetSensorID() == idsFirstRing[j]){
				cluster.push_back(plane[i]);
				break;
			}
		}
	}

/*	std::cout << "Selected: ";
	for(unsigned int i=0;i<cluster.size();i++){
		std::cout << "\t" << cluster[i]->GetSensorID();
	}
	std::cout << std::endl;
*/
	return true;
}

