#include "RNG_base.hpp"
#include <complex>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>

double RNG_base::mean(const std::vector<double>& v) const {
  double dSum{std::accumulate(v.begin(), v.end(), 0.0)};
  return dSum / v.size();
}

double RNG_base::stdDev(const std::vector<double>& v) const {
  double dMean{mean(v)};
  double dSS{0.0};
  for (const auto& e : v) {
    dSS += (e - dMean) * (e - dMean);
  }
  return sqrt(dSS / (v.size() - 1));
}

void RNG_base::normalize(std::vector<double>& v) {
  double dMean{mean(v)};
  double dStdDev{stdDev(v)};
  for (auto& e : v) {
    e = (e - dMean) / dStdDev;
  }
}

std::vector<double> RNG_base::expDecay(const int& n,
				       const double& lambda,
				       const double& sdExtraNoise) {
  std::vector<double> vf(n, 0.0);

  for (int i = 0; i < n; i++) {
    double step = (double) i / n;
    vf[i] = exp(-exp(lambda) * step);
  }

  normalize(vf);

  std::vector<double> gn{rnorm2(n, 0, sdExtraNoise)};

  for (size_t i = 0; i < vf.size(); i++) {
    vf[i] += gn[i];
  }

  normalize(vf);
  
  return vf;
}

std::vector<double> RNG_base::randomWalk(const int& n,
					 const double& sdExtraNoise) {
  // implementation of R autocorr::stat_gp() with
  // sigma = 1, gamma = 2, and dt_GP = 0.8
  std::vector<double> vf{};

  double sigma = 1.0;
  double gamma = 2.0;
  
  // double Nt_cor = 100;
  // double dt = gamma / Nt_cor;
  // double eps_cov = .0001;
  // double Tmax = 430.0 * dt;
  // double Nt = Tmax / dt;
  // double int_0 = Nt + 1.0;

  const int Nw = 128;
  const int Mw = 550;
  const double dts = 1.72;
  const double dw2 = 0.01426959;

  double cov_t[2 * Nw];
  for (int i = 0; i < 2 * Nw; i++) {
    cov_t[i] = -128 + i;
  }
  for (int i = 0; i < 2 * Nw; i++) {
    double d = cov_t[i] * dts;
    cov_t[i] = sigma * sigma * exp(-(d * d) / (2 * gamma * gamma));
  }

  // run FFT on cov_t
  gsl_fft_real_wavetable * gsl_wavt =
    gsl_fft_real_wavetable_alloc(2 * Nw);
  gsl_fft_real_workspace * gsl_work =
    gsl_fft_real_workspace_alloc(2 * Nw);

  gsl_fft_real_transform(cov_t, 1, 2 * Nw, gsl_wavt, gsl_work);
  
  gsl_fft_real_wavetable_free(gsl_wavt);
  gsl_fft_real_workspace_free(gsl_work);

  const double d_pi = 2 * acos(0.0);
  
  std::vector<double> S2(Nw, 0.0);
  
  for (int i = 0; i < Nw; i++) {
    if (i > 1) {
      S2[i] = abs(cov_t[(i - 2) * 2 + 3]) / (2 * d_pi);
    } else {
      S2[i] = abs(cov_t[i]) / (2 * d_pi);
    }
  }

  std::vector<double> A(Mw, 0.0);
  for (int i = 0; i < Nw; i++) {
    if (i > 0) {
      A[i] = sqrt(2 * S2[i] * dw2);
    } else {
      A[i] = 0.0;
    }
  }

  std::vector<double> vRUnif{runif2(Nw)};
  std::vector<double> phi(Mw, 0.0);
  for (int i = 0; i < Nw; i++) {
    phi[i] = 2 * d_pi * vRUnif[i];
  }

  using namespace std::complex_literals;
  std::vector<std::complex<double>> B(Mw, 0.0+0.0i);
  for (size_t i = 0; i < B.size(); i++) {
    B[i] = sqrt(2) * A[i] * std::exp(0.0+1.0i * phi[i]);
  }
  // for (int i = 0; i < Nw; i++) {
  //   v[i] = cov_t[i];
  // }

  // create packed array for complex inverse FFT
  double * dpa = new double[B.size() * 2];

  // even indices will be the real part
  // odd indices are imaginary
  for (size_t i = 0; i < B.size(); i++) {
    dpa[i * 2] = B[i].real();
    dpa[i * 2 + 1] = B[i].imag();
  }

  gsl_fft_complex_wavetable * g_wavtbl = gsl_fft_complex_wavetable_alloc(B.size());
  gsl_fft_complex_workspace * g_work = gsl_fft_complex_workspace_alloc(B.size());

  gsl_fft_complex_inverse(dpa, 1, B.size(), g_wavtbl, g_work);
  
  gsl_fft_complex_wavetable_free(g_wavtbl);
  gsl_fft_complex_workspace_free(g_work);

  std::vector<double> vResult(n, 0.0);
  for (int i = 0; i < n; i++) {
    vResult[i] = dpa[i * 2];
  }

  normalize(vResult);
  
  // Need to add a small amount of noise
  std::vector<double> vExtra{rnorm2(n, 0, sdExtraNoise)};
  for (int i = 0; i < n; i++) {
    vResult[i] += vExtra[i];
  }

  normalize(vResult);
  
  delete [] dpa;
  
  return vResult;
}
