#include <string>
#include <GATE/Centella.h>

namespace util
{

  ///
  /// Barycentre algorithm for xy reconstruction from tracking
  /// plane information.
  ///

  class barycenterAlgorithm 
  {
  public:

    // Constructor
    barycenterAlgorithm();

    // Destructor
    ~barycenterAlgorithm() {}

    // Calculation function
    bool computePosition(const vector<gate::Hit*>& sensors);
    bool computePositionCoronna0(const vector<gate::Hit*>& sensors, double thresholdMax, double thresholdNeighbours);
	static bool chargeOrderSensorsDec(const gate::Hit* s1, const gate::Hit* s2);

    // Getters for the information calculated.
    double getX1() const;
    double getX2() const;
    double getX1Err() const;
    double getX2Err() const;
    void getX1withErr(std::pair<double,double>& xErr) const;
    void getX2withErr(std::pair<double,double>& yErr) const;
    //double getXSig() const override;
    //double getYSig() const override;
    //
	
	//Set plane type
	void setPlane(std::string planeDirection);

  private:

    // Clear the data members out.
    void clearData();

    // Calculate the sigmas.
    //void calculateSigmas(const std::map<int,double>& sensors,
	//		 double engSum, double sqEngSum);

    // The information that the algorithm will calculate and save.
    double x1Pos_;
    double x2Pos_;

    double x1Err_;
    double x2Err_;

	// Code: 0 -> XY, 1-> XZ, 2 -> YZ
	int plane_;

    //double xSig_;
    //double ySig_;
    //

  }; // class barycenterAlgorithm

  /// INLINE METHODS

  inline double barycenterAlgorithm::getX1() const { return x1Pos_; }
  inline double barycenterAlgorithm::getX2() const { return x2Pos_; }
  inline double barycenterAlgorithm::getX1Err() const { return x1Err_; }
  inline double barycenterAlgorithm::getX2Err() const { return x2Err_; }
  //inline double barycenterAlgorithm::getXSig() const { return xSig_; }
  //inline double barycenterAlgorithm::getYSig() const { return ySig_; }

} // namespace util
