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
  
  //Fill energy histogram
  void energyHist(gate::Event& evt);


  //! finalize algorithm
  bool finalize();          
  
 private:
  
  ClassDef(petAnalysis,0)
    
};

#endif
