#include <numbers>

#include <Rcpp.h>
using namespace Rcpp;

#include <memory>

#include "RNG_r.hpp"
#include "SimData.hpp"

double mean(const NumericVector& v) {
  double dSum{std::accumulate(v.begin(), v.end(), 0.0)};
  return dSum / v.size();
}

double stdDev(const NumericVector& v) {
  double dMean{mean(v)};
  double dSS{0.0};
  for (const auto& e : v) {
    dSS += (e - dMean) * (e - dMean);
  }
  return sqrt(dSS / (v.size() - 1));
}

void normalize(NumericVector& v) {
  double dMean{mean(v)};
  double dStdDev{stdDev(v)};
  for (auto& e : v) {
    e = (e - dMean) / dStdDev;
  }
}

// [[Rcpp::export]]
std::vector<double> test_pink(const int& n, const double& d) {
  std::unique_ptr<RNG_r> pRNG = std::make_unique<RNG_r>();

  return pRNG->pinkNoise(n, d);
}

// [[Rcpp::export]]
std::vector<int> test_psr(const int& n_obs, const int& n_subblocks) {
  std::unique_ptr<RNG_r> pRNG = std::make_unique<RNG_r>();
  std::vector<int> v(n_obs);
  std::iota(v.begin(), v.end(), 1);

  pRNG->permuteSubblocks(v, n_subblocks);

  return v;
}

// [[Rcpp::export]]
DataFrame simdata_1f_int(int nPart, int nLevels, int nReps, int nSubblocks,
			 double mu, double eta,
			 double prop_rint, double prop_rslp,
			 int estr_id) {

  std::unique_ptr<RNG_r> pRNG = std::make_unique<RNG_r>();
  SimDataOneFactor sd(nPart, nLevels, nReps, std::move(pRNG));
  sd.generate(mu, eta, prop_rint, prop_rslp,
	      static_cast<errstr>(estr_id),
	      nSubblocks);

  const auto tData = sd.getData();

  DataFrame df = DataFrame::create( Named("id") = std::get<0>(tData),
				    Named("cond") = std::get<1>(tData),
				    Named("order") = std::get<2>(tData),
				    Named("dv") = std::get<3>(tData) );

  return df;
}

// [[Rcpp::export]]
DataFrame simdata_2x2ww_int(int nPart, int nReps, int nSubblocks,
			    double mu,
			    double eta_row, double eta_col, double eta_rcx,
			    double prop_rint,
			    double prop_rslp_row,
			    double prop_rslp_col,
			    double prop_rslp_rcx,
			    int estr_id) {

  std::unique_ptr<RNG_r> pRNG = std::make_unique<RNG_r>();
  SimData2x2 sd(nPart, nReps, std::move(pRNG));
  sd.generate(mu,
	      eta_row, eta_col, eta_rcx,
	      prop_rint,
	      prop_rslp_row, prop_rslp_col, prop_rslp_rcx,
	      static_cast<errstr>(estr_id),
	      PSR,
	      nSubblocks);

  const auto tData = sd.getData();

  DataFrame df = DataFrame::create( Named("id") = std::get<0>(tData),
				    Named("A") = std::get<1>(tData),
				    Named("B") = std::get<2>(tData),
				    Named("order") = std::get<3>(tData),
				    Named("dv") = std::get<4>(tData) );

  return df;
}
