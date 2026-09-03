#include "RNG_r.hpp"
#include <Rcpp.h>

RNG_r::RNG_r() {
}

RNG_r::~RNG_r() {
}

void RNG_r::permuteSubblocks(std::vector<int>& v,
			     const int& nSubblocks) {

  const int sbSize = v.size() / nSubblocks;

  std::vector<int> vCopy(v);
  int ixFrom{0}; // index to span subblock: [ixFrom, ixFrom + sbSize)
  
  for (int s = 0; s < nSubblocks; s++) {
    ixFrom = s * sbSize;

    auto vsb = Rcpp::sample(sbSize, sbSize, false, R_NilValue, false);

    for (int i = 0; i < sbSize; i++) {
      v[ixFrom + i] = vCopy[ixFrom + vsb[i]];
    }
  }
}

std::vector<double> RNG_r::rnorm2(const int& n,
				  const double& mean,
				  const double& sd) {

  std::vector<double> vf =
    Rcpp::as<std::vector<double>>(Rcpp::rnorm(n, mean, sd));

  return vf;
}

std::vector<double> RNG_r::runif2(const int& n,
				  const double& min,
				  const double& max) {
  std::vector<double> vf =
    Rcpp::as<std::vector<double>>(Rcpp::runif(n, min, max));

  return vf;
}

std::vector<double> RNG_r::pinkNoise(const int& n,
				     const double& alpha) {

  Rcpp::Function f("pink");

  Rcpp::NumericVector nv = f(n, alpha);

  std::vector<double> vd(nv.begin(), nv.end());

  normalize(vd);

  return vd;
}
