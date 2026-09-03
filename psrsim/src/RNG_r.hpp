#ifndef RNG_r_INCLUDED
#define RNG_r_INCLUDED

#include "RNG_base.hpp"

// this version is for using with R
class RNG_r : public RNG_base {
public:
  RNG_r();
  ~RNG_r();
  
  std::vector<double> rnorm2(const int& n,
			     const double& mean = 0,
			     const double& sd = 1) override;

  std::vector<double> runif2(const int& n,
			     const double& min = 0,
			     const double& max = 1) override;

  std::vector<double> pinkNoise(const int& n,
				const double& alpha) override;
  
  void permuteSubblocks(std::vector<int>& v,
			const int& nSubblocks) override;
};

#endif
