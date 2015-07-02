#ifndef _petAnalysis__
#define _petAnalysis__

#include <GATE/Centella.h>

class petAnalysis : public gate::IAlgo {

 public:
  
  //! default contructor
  petAnalysis(gate::VLEVEL=gate::NORMAL,
	       std::string label="petAnalysisInstance");
  
  //! constructor with store with input parameters 
  petAnalysis(const gate::ParamStore& gs,
	       gate::VLEVEL=gate::NORMAL,
	       std::string label="petAnalysisInstance");
  
  //! destructor
  virtual ~petAnalysis(){};
  
  //! initialize algorithm
  bool initialize();        
  
  //! execute algorithm: process current event
  bool execute(gate::Event& evt);  
  
  //Position reconstruction using barycenter
  void reconstruction(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt);
  void reconstructionNoNorm(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& pt);
  void bestPointRecons(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt);
  void bestPointReconsNoNorm(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt);
  //Fill energy histogram
  void energyHist(gate::Event& evt);
  //Find first particle in a vector of particles by its creation time
  void findFirstParticle(const std::vector<const gate::MCParticle*> particles, gate::MCParticle& first);
  void findFirstParticle(std::vector<gate::MCParticle*> particles, gate::MCParticle& first);
  //Classify events as compton or photoelectric
  void classifyEvent(gate::MCParticle& primary, gate::MCParticle& firstDaughter);
  //Order particles
  static bool timeOrderParticles(const gate::MCParticle* p1, const gate::MCParticle* p2);
  //Order sensors by charge (descending)
  static bool chargeOrderSensorsDesc(const gate::Hit* s1, const gate::Hit* s2);
  static bool chargeOrderSensorsAsc(const gate::Hit* s1, const gate::Hit* s2);
  //Compute distance between two Point3D
  double distance(gate::Point3D& p1, gate::Point3D& p2);
  void hist2dEvent(gate::Event& evt);
  void hist2dHits(gate::Event& evt);
  void splitHitsPerPlane(gate::Event& evt, std::vector<std::vector<gate::Hit*> >& planes);
  void applyCut(const std::vector<gate::Hit*>& sensorHits, double cut, std::vector<gate::Hit*>& filtered);
  bool nearPlane(gate::Point3D& pt,double distance);
  void fillComptonHist(gate::MCParticle& primary);
  void reconsPerPlane(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt);
  void printSensors(std::vector<std::vector<gate::Hit*> >& planes);
  double findSensors(std::vector<gate::Hit*>& plane, int id);


  //! finalize algorithm
  bool finalize();          
  
 private:
  
  ClassDef(petAnalysis,0)
    
};

#endif
