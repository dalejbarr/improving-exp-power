#ifndef RNG_INCLUDED
#define RNG_INCLUDED

#include <vector>
#include <random>
#include <array>

class RNG_base {
public:
  // non-virtual member functions
  void normalize(std::vector<double>& v);
  
  double mean(const std::vector<double>& v) const;
  double stdDev(const std::vector<double>& v) const;

  virtual std::vector<double> rnorm2(const int& n,
				     const double& mean = 0,
				     const double& sd = 1) = 0;

  virtual std::vector<double> runif2(const int& n,
				     const double& min = 0,
				     const double& max = 1) = 0;

  virtual std::vector<double> pinkNoise(const int& n,
					const double& alpha) = 0;
  
  virtual void permuteSubblocks(std::vector<int>& v,
				const int& nSubblocks) = 0;

  
  std::vector<double> expDecay(const int& n,
			       const double& lambda,
			       // little extra noise (sqrt(.1) = .3162278)
			       const double& sdExtraNoise = 0.3162278);
  
  std::vector<double> randomWalk(const int& n,
				 // little extra noise (sqrt(.1) = .3162278)
				 const double& sdExtraNoise = 0.3162278);
};

#endif
