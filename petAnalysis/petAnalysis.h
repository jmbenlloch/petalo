#ifndef _petAnalysis__
#define _petAnalysis__

#include <GATE/Centella.h>
#include <utility>
#include "TNtuple.h"
#include "TFile.h"

class petAnalysis : public gate::IAlgo {

	TTree *_tree;
	TFile *_file;
	double _entryPlane[64];
	double _x[1], _y[1], _z[1];
	double _maxEntry[1], _maxExit[1], _totalEntry[1], _totalExit[1], _ratio[1];

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
  void bestPointRecons(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& truePt, gate::Point3D& pt);
  //Find first particle in a vector of particles by its creation time
  void findFirstParticle(const std::vector<const gate::MCParticle*> particles, gate::MCParticle& first);
  void findFirstParticle(std::vector<gate::MCParticle*> particles, gate::MCParticle& first);
  //Order particles
  static bool timeOrderParticles(const gate::MCParticle* p1, const gate::MCParticle* p2);
  //Order sensors by charge (descending)
  static bool chargeOrderSensorsDesc(const gate::Hit* s1, const gate::Hit* s2);
  static bool chargeOrderSensorsAsc(const gate::Hit* s1, const gate::Hit* s2);
  static bool chargeOrderPlanesDesc(std::pair<int,double> s1, std::pair<int,double> s2);
  //Compute distance between two Point3D
  double distance(gate::Point3D& p1, gate::Point3D& p2);
  void splitHitsPerPlane(gate::Event& evt, std::vector<std::vector<gate::Hit*> >& planes);
  void applyCut(const std::vector<gate::Hit*>& sensorHits, double cut, std::vector<gate::Hit*>& filtered);
  double planeCharge(std::vector<gate::Hit*> plane);
  void computeBarycenters(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<double> >& points, std::vector<std::vector<double> >& errors);
  double zReconsRatio(double ratio);
  void chargeHist2d(std::vector<std::vector<gate::Hit*> > planes);
  void sipmmcHist(std::vector<std::vector<gate::Hit*> > planes, gate::Point3D& trueVertex, int phot);
  double totalCharge(std::vector<gate::Hit*> hits);
  void filterHits(std::vector<gate::Hit*> plane, std::vector<gate::Hit*> hits, std::vector<gate::Hit*>& planeFiltered);
  double findSensors(std::vector<gate::Hit*>& plane, int id);

  //! finalize algorithm
  bool finalize();          

  TTree * getTree() {return _tree;}
  void setTree(TTree* tree) {_tree = tree;}
  TFile * getFile() {return _file;}
  void setFile(TFile* file) {_file = file;}
  double * getEntryPlane() {return _entryPlane;}
  double * getX() {return _x;}
  double * getY() {return _y;}
  double * getZ() {return _z;}

  double * getMaxEntry() {return _maxEntry;}
  double * getMaxExit() {return _maxExit;}
  double * getTotalEntry() {return _totalEntry;}
  double * getTotalExit() {return _totalExit;}
  double * getRatio() {return _ratio;}

 private:
  
  ClassDef(petAnalysis,0)
    
};

#endif
