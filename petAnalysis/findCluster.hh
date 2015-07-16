#include <string>
#include <GATE/Centella.h>

namespace util
{

  class findCluster
  {
  public:

    // Constructor
    findCluster();

    // Destructor
    ~findCluster() {}

    // Calculation function
    bool findCoronna(const std::vector<gate::Hit*>& plane, std::vector<gate::Hit*>& cluster, double thresholdMax, double thresholdNeighbours);
    bool findCoronna0(const std::vector<gate::Hit*>& plane, std::vector<gate::Hit*>& cluster, double thresholdMax, double thresholdNeighbours);
    bool findCoronna2Rings(const std::vector<gate::Hit*>& plane, std::vector<gate::Hit*>& cluster, double thresholdMax, double thresholdNeighbours);
	bool findCoronnaAllPlanes(const std::vector<std::vector<gate::Hit*> >& planes, std::vector<std::vector<gate::Hit*> >& clusters, int rings, double thresholdMax, double thresholdNeighbours);


	static bool chargeOrderSensorsDec(const gate::Hit* s1, const gate::Hit* s2);


 // private:

  };

} // namespace util
