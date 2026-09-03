#include "SimData.hpp"
#include <gsl/gsl_cdf.h>
#include <iostream>

SimDataBase::SimDataBase(const int& nParticipants,
			 const int& nReps,
			 std::unique_ptr<RNG_base> pRNG) :
  m_nParticipants{nParticipants},
  m_nReps{nReps},
  m_pRNG{std::move(pRNG)} {
}

SimDataBase::~SimDataBase() {
}

double SimDataBase::sumOfSquares(const std::vector<double>& vd) const {
  double vdResult{0.0};

  for (auto& e : vd) {
    vdResult += e * e;
  }

  return vdResult;
}

std::vector<std::vector<double>> SimDataBase::errors(const enum errstr& estr,
						     const double& sdExtraNoise) {
    
  // ERRORS
  std::vector<std::vector<double>> vErr(m_nParticipants,
					std::vector<double>(NTrials(), 0));

  switch (estr) {
  case noAuto :
    for (int p = 0; p < m_nParticipants; p++) {
      vErr[p] = m_pRNG->rnorm2(vErr[p].size(), 0.0, 1.0);
    }
    break;
  case expDecay :
    {
      std::vector<double> vLambdaBase{m_pRNG->runif2(1, log(2), log(15))};
      // add a random effect for each participant
      std::vector<double> vOffset{m_pRNG->rnorm2(m_nParticipants, 0.0, 0.8)};

      for (int p = 0; p < m_nParticipants; p++) {
	double dLambda = exp(vLambdaBase[0] + vOffset[p]);
	vErr[p] = m_pRNG->expDecay(vErr[p].size(), dLambda, sdExtraNoise);
      }
    }
    break;
  case pinkNoise :
    {
      double dAlpha = m_pRNG->runif2(1, .8, 1.2)[0];
      for (int p = 0; p < m_nParticipants; p++) {
	vErr[p] = m_pRNG->pinkNoise(vErr[p].size(), dAlpha);
      }
    }
    break;
  case randomWalk :
    {
      for (int p = 0; p < m_nParticipants; p++) {
	vErr[p] = m_pRNG->randomWalk(vErr[p].size(), sdExtraNoise);
      }
    }
    break;
  case mixed :
    {
      std::vector<std::vector<double>> vDecay(errors(expDecay, 0.0));
      std::vector<std::vector<double>> vPink(errors(pinkNoise));
      std::vector<std::vector<double>> vWalk(errors(randomWalk, 0.0));
      for (int p = 0; p < m_nParticipants; p++) {
	// need positive weights that sum to one
	std::vector<double> vWeights{m_pRNG->runif2(3)};
        const double dSum{std::accumulate(vWeights.begin(),
					  vWeights.end(), 0.0)};
	for (int i = 0; i < 3; i++) {
	  vWeights[i] /= dSum;
	}
	for (size_t i = 0; i < vDecay[p].size(); i++) {
	  vErr[p][i] = vDecay[p][i] * vWeights[0] +
	    vPink[p][i] * vWeights[1] +
	    vWalk[p][i] * vWeights[2];
	}
	m_pRNG->normalize(vErr[p]);

	// now add some noise
	std::vector<double> dvExtra = m_pRNG->rnorm2(vErr[p].size(),
						     0, sdExtraNoise);
	for (size_t i = 0; i < vErr[p].size(); i++) {
	  vErr[p][i] += dvExtra[i];
	}
	m_pRNG->normalize(vErr[p]);
      }
    }
    break;
  }

  return vErr;
}

SimDataOneFactor::SimDataOneFactor(const int& nParticipants,
				   const int& nLevels,
				   const int& nReps,
				   std::unique_ptr<RNG_base> pRNG) :
  SimDataBase(nParticipants, nReps, std::move(pRNG)),
  m_nLevels{nLevels},
  m_pData{nullptr},
  m_pOrder{nullptr},
  m_nSignificant{0},
  m_vdPValues{0.0} {
  
  m_pData = new double **[nParticipants];
  m_pOrder = new int **[nParticipants];
  for (int i = 0; i < nParticipants; i++) {
    m_pData[i] = new double *[nLevels];
    m_pOrder[i] = new int *[nLevels];
    for (int j = 0; j < nLevels; j++) {
      m_pData[i][j] = new double[nReps];
      m_pOrder[i][j] = new int[nReps];
    }
  }

}

SimDataOneFactor::~SimDataOneFactor() {
  for (int i = 0; i < m_nParticipants; i++) {
    for (int j = 0; j < m_nLevels; j++) {
      delete[] m_pData[i][j];
      delete[] m_pOrder[i][j];
    }
    delete[] m_pData[i];
    delete[] m_pOrder[i];
  }
  delete[] m_pData;
  delete[] m_pOrder;
}

std::vector<double> SimDataOneFactor::MainEffOneFactor(const double& targSS /* = 1 */) {

  // make random numbers, center (subtract mean)
  std::vector<double> vf = m_pRNG->runif2(m_nLevels);

  double sum = std::accumulate(vf.begin(), vf.end(), 0.0);
  double mean = sum / vf.size();
  
  for (size_t i = 0; i < vf.size(); i++) {
    vf[i] = vf[i] - mean;
  }

  // calculate sum of squares, and figure out correction
  // factor so that sum of squares matches targSS
  double ss = sumOfSquares(vf);
  double fac = targSS / ss;
  
  for (size_t i = 0; i < vf.size(); i++) {
    vf[i] *= sqrt(fac);
  }
  
  return vf;
}


void SimDataOneFactor::generate(const double& mu,
		       const double& eta,
		       const double& prop_rint,
		       const double& prop_rslp,
		       const enum errstr& estr,
		       const int& nSubblocks) {
  
  //if ((eta > 1) || (eta < 0)) stop("eta must be between 0 and 1")

  randomize(nSubblocks);
  
  // we are using the property that
  // expected value of SS = N * sigma^2
  //
  // error total is:
  //   SS_errtot = SS_trialerr + SS_randint + SS_randslp
  
  // -- calculations
  double SS_trialerr = m_nLevels * m_nReps * m_nParticipants; // sigma = 1
  double SS_randint = prop_rint * SS_trialerr;
  double SS_randslp = prop_rslp * SS_trialerr;

  // MAIN EFFECT
  // calculate target SS_x needed to get eta
  double SS_x = (eta * (SS_trialerr + SS_randslp)) / (1 - eta);
  double targ_ss = SS_x / (m_nParticipants * m_nReps);

  std::vector<double> mainEff{MainEffOneFactor(targ_ss)};

  // RANDOM EFFECTS
  double targ_ss_rint = SS_randint / (m_nReps * m_nLevels);
  double rint_sd = sqrt(targ_ss_rint / m_nParticipants);
  std::vector<double> randomIntercepts{m_pRNG->rnorm2(m_nParticipants, 0.0, rint_sd)};
  
  double targ_ss_rslp = SS_randslp / (m_nParticipants * m_nReps);
  
  std::vector<std::vector<double>> randomSlopes(m_nParticipants,
						std::vector<double>(m_nLevels *
								    m_nReps, 0));

  for (int i = 0; i < m_nParticipants; i++) {
    randomSlopes[i] = MainEffOneFactor(targ_ss_rslp);
  }

  std::vector<std::vector<double>> vErr{errors(estr, sqrt(.1))};

  // apply linear model to calculate the DV
  for (int p = 0; p < m_nParticipants; p++) {
    for (int k = 0; k < m_nLevels; k++) {
      for (int r = 0; r < m_nReps; r++) {
	m_pData[p][k][r] = mu + mainEff[k] + // fixed effects
	  randomIntercepts[p] + randomSlopes[p][k] + // random effects
	  vErr[p][m_pOrder[p][k][r]];
      }
    }
  }
}

void SimDataOneFactor::generate2(const double& mu,
				 const double& eta,
				 const double& prop_rint,
				 const double& prop_rslp,
				 const enum errstr& estr_int,
				 const enum errstr& estr_slp,
				 const int& nSubblocks) {
  
  //if ((eta > 1) || (eta < 0)) stop("eta must be between 0 and 1")

  randomize(nSubblocks);
  
  // we are using the property that
  // expected value of SS = N * sigma^2
  //
  // error total is:
  //   SS_errtot = SS_trialerr + SS_randint + SS_randslp
  
  // -- calculations
  double SS_trialerr = m_nLevels * m_nReps * m_nParticipants; // sigma = 1
  double SS_randint = prop_rint * SS_trialerr;
  double SS_randslp = prop_rslp * SS_trialerr;

  // MAIN EFFECT
  // calculate target SS_x needed to get eta
  double SS_x = (eta * (SS_trialerr + SS_randslp)) / (1 - eta);
  double targ_ss = SS_x / (m_nParticipants * m_nReps);

  std::vector<double> mainEff{MainEffOneFactor(targ_ss)};

  // RANDOM EFFECTS
  double targ_ss_rint = SS_randint / (m_nReps * m_nLevels);
  double rint_sd = sqrt(targ_ss_rint / m_nParticipants);
  std::vector<double> randomIntercepts{m_pRNG->rnorm2(m_nParticipants, 0.0, rint_sd)};
  
  double targ_ss_rslp = SS_randslp / (m_nParticipants * m_nReps);
  
  std::vector<std::vector<double>> randomSlopes(m_nParticipants,
						std::vector<double>(m_nLevels *
								    m_nReps, 0));

  for (int i = 0; i < m_nParticipants; i++) {
    randomSlopes[i] = MainEffOneFactor(targ_ss_rslp);
  }

  std::vector<std::vector<double>> vErr{errors(estr_int, sqrt(.1))};

  // apply linear model to calculate the DV
  for (int p = 0; p < m_nParticipants; p++) {
    for (int k = 0; k < m_nLevels; k++) {
      for (int r = 0; r < m_nReps; r++) {
	m_pData[p][k][r] = mu + mainEff[k] + // fixed effects
	  randomIntercepts[p] + randomSlopes[p][k] + // random effects
	  vErr[p][m_pOrder[p][k][r]];
      }
    }
  }
}

void SimDataOneFactor::writeCSV(const std::string& fname) const {
  auto strNP = std::to_string(m_nParticipants);

  std::ofstream fout(fname);
  fout << "id,cond,tnum,dv" << std::endl;
  
  for (int p = 0; p < m_nParticipants; p++) {
    auto strID = std::to_string(p + 1);
    if (strID.size() < strNP.size()) {
      strID = std::string(strNP.size() - strID.size(), '0') + strID;
    }
    for (int k = 0; k < m_nLevels; k++) {
      auto strCond = std::string(1, 'A') + std::to_string(k + 1);
      for (int r = 0; r < m_nReps; r++) {
	fout << "\"P" << strID << "\",\"" << strCond << "\","
	     << (m_pOrder[p][k][r] + 1) << ","
	     << m_pData[p][k][r] << std::endl;
      }
    }
  }
}

anovaStats SimDataOneFactor::anova() {
  anovaStats result(m_nLevels);

  double grandSum{0.0};
  std::vector<double> mainSum(m_nLevels, 0.0);
  std::vector<double> partEff(m_nParticipants, 0.0);
  std::vector<std::vector<double>> partCondEff(m_nParticipants,
					       std::vector<double>(m_nLevels, 0.0));

  for (int p = 0; p < m_nParticipants; p++) {
    for (int k = 0; k < m_nLevels; k++) {
      for (int r = 0; r < m_nReps; r++) {
	partEff[p] += m_pData[p][k][r];
	partCondEff[p][k] += m_pData[p][k][r];
	mainSum[k] += m_pData[p][k][r];
	grandSum += m_pData[p][k][r];
	// result.ssTot += m_pData[p][k][r] * m_pData[p][k][r];
      }
      partCondEff[p][k] = partCondEff[p][k] / m_nReps;
    }
  }

  // calculate main effects
  result.grandMean = grandSum / (m_nParticipants * m_nLevels * m_nReps);
  
  for (int k = 0; k < m_nLevels; k++) {
    
    result.mainEff[k] = mainSum[k] / (m_nParticipants * m_nReps) -
      result.grandMean;
    
    result.ssEff += result.mainEff[k] * result.mainEff[k] *
      m_nParticipants * m_nReps;
  }

  // finish calculating participant and condition effects
  for (int p = 0; p < m_nParticipants; p++) {
    partEff[p] = partEff[p] / (m_nLevels * m_nReps) - result.grandMean;
    for (int k = 0; k < m_nLevels; k++) {
      partCondEff[p][k] = partCondEff[p][k] - partEff[p] - result.mainEff[k] -
	result.grandMean;
      result.ssError += (partCondEff[p][k] * partCondEff[p][k]) * m_nReps;
    }
  }
  
  result.dfNum = m_nLevels - 1;

  result.dfDen = (m_nParticipants - 1) * (m_nLevels - 1);

  double msEff = result.ssEff / result.dfNum;
  double msErr = result.ssError / result.dfDen;

  result.FRatio = msEff / msErr;
  result.pValue = 1.0 - gsl_cdf_fdist_P(result.FRatio,
					result.dfNum, result.dfDen);
  
  return result;
}

void SimDataOneFactor::run(const int& nmc, const double& eta2,
			   const enum errstr& estr, const int& nSubblocks) {
  
  m_nSignificant = 0;
  anovaStats result(m_nLevels);
  
  for (int i = 0; i < nmc; i++) {
    generate(0.0, eta2, .35, .11, estr, nSubblocks);
    result = anova();
    if (result.pValue < .05) {
      m_nSignificant++;
    }
  }
}

void SimDataOneFactor::run2(const int& nmc, // no. Monte Carlo runs
			    const double& eta2, // eta-squared
			    const enum errstr& estr_int, // err structure icept
			    const enum errstr& estr_slp, // err structure slope
			    const int& nSubblocks, 
			    const double& dMu, //= 0.0, // intercept (grand mean)
			    const double& prop_rint, // = .35,
			    const double& prop_rslp) { // = .11) 

  m_vdPValues.resize(nmc);

  anovaStats result(m_nLevels);

  for (int i = 0; i < nmc; i++) {
    generate2(dMu, eta2, prop_rint, prop_rslp,
	      estr_int, estr_slp, nSubblocks);
    result = anova();
    m_vdPValues[i] = result.pValue;
  }
  
}

void SimDataOneFactor::randomize(const int& nSubblocks) {

  const int sbSize = (m_nLevels * m_nReps) / nSubblocks;

  std::vector<int> v(m_nLevels * m_nReps, 0);
  std::iota(v.begin(), v.end(), 0);
  
  for (int p = 0; p < m_nParticipants; p++) {
    m_pRNG->permuteSubblocks(v, nSubblocks);

    for (int s = 0; s < nSubblocks; s++) {
      int ix = s * sbSize;
      for (int se = 0; se < sbSize; se++) {
	int k = se / (sbSize / m_nLevels);
	int r = se % (sbSize / m_nLevels) + (sbSize / m_nLevels) * s;
	m_pOrder[p][k][r] = v[ix + se];
      }
    }
  }
  
}

void SimDataOneFactor::writePValues(std::ofstream& out) {
  for (auto d : m_vdPValues) {
    out << d << "\n";
  }
}

SimData2x2::SimData2x2(const int& nParticipants,
		       const int& nReps,
		       std::unique_ptr<RNG_base> pRNG) :
  SimDataBase(nParticipants, nReps, std::move(pRNG)), 
  m_pData{nullptr},
  m_pOrder{nullptr} {

  m_vSignificant = std::vector<int>(3, 0);
  
  m_pData = new double ***[nParticipants];
  m_pOrder = new int ***[nParticipants];
  for (int p = 0; p < nParticipants; p++) {
    m_pData[p] = new double **[2];
    m_pOrder[p] = new int **[2];
    for (int r = 0; r < 2; r++) {
      m_pData[p][r] = new double *[2];
      m_pOrder[p][r] = new int *[2];
      for (int c = 0; c < 2; c++) {
	m_pData[p][r][c] = new double[nReps];
	m_pOrder[p][r][c] = new int[nReps];
      }
    }
  }
}

SimData2x2::~SimData2x2() {
  for (int p = 0; p < m_nParticipants; p++) {
    for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) { 
	delete[] m_pData[p][r][c];
	delete[] m_pOrder[p][r][c];
      }
      delete[] m_pData[p][r];
      delete[] m_pOrder[p][r];
    }
    delete[] m_pData[p];
    delete[] m_pOrder[p];
  }
  delete[] m_pData;
  delete[] m_pOrder;
}

void SimData2x2::generate(const double& mu,
			  const double& eta_row,
			  const double& eta_col,
			  const double& eta_rcx,
			  const double& prop_rint,
			  const double& prop_rslp_row,
			  const double& prop_rslp_col,
			  const double& prop_rslp_rcx,
			  const enum errstr& estr,
			  const enum facRandStrategy& randStrategy,
			  const int& nSubblocks) {

  if (randStrategy == PSR) {
    randomize(nSubblocks);
  } else if (randStrategy == circular) {
    circularWalk();
  } else if (randStrategy == fig8) {
    fig8Walk();
  }
  
  // we are using the property that
  // expected value of SS = N * sigma^2
  //
  // error total is:
  //   SS_errtot = SS_trialerr + SS_M + SS_A + SS_B + SS_AB
  
  // -- calculations
  double SS_trialerr = 2 * 2 * m_nReps * m_nParticipants; // sigma = 1
  double SS_randint = prop_rint * SS_trialerr;
  double SS_randA = prop_rslp_row * SS_trialerr;
  double SS_randB = prop_rslp_col * SS_trialerr;
  double SS_randAB = prop_rslp_rcx * SS_trialerr;

  // calculate target SS_x needed to get eta_x
  double SS_a = (eta_row * (SS_trialerr + SS_randA)) / (1 - eta_row);
  double SS_b = (eta_col * (SS_trialerr + SS_randB)) / (1 - eta_col);
  double SS_ab = (eta_rcx * (SS_trialerr + SS_randAB)) / (1 - eta_rcx);

  // (true) main effects; work backwards from SS
  double val_a = sqrt(SS_a / (2 * 2 * m_nReps * m_nParticipants));
  double val_b = sqrt(SS_b / (2 * 2 * m_nReps * m_nParticipants));
  double val_ab = sqrt(SS_ab / (2 * 2 * m_nReps * m_nParticipants));
  std::vector<double> va{-val_a, val_a};
  std::vector<double> vb{-val_b, val_b};
  std::vector<std::vector<double>> vab{ {val_ab, -val_ab}, {-val_ab, val_ab} };

  // random effects
  // random intercepts: work backwards
  double targ_ss_rint = SS_randint / (m_nReps * 2 * 2);
  double rint_sd = sqrt(targ_ss_rint / m_nParticipants);
  std::vector<double> vri{m_pRNG->rnorm2(m_nParticipants, 0.0, rint_sd)};

  // random slope for A
  double targ_ss_randA = SS_randA / (m_nReps * 2 * 2);
  double randA_sd = sqrt(targ_ss_randA / m_nParticipants);
  std::vector<std::vector<double>> vvra(m_nParticipants,
					std::vector<double>(2, 0.0));
  std::vector<double> randomSlopes{m_pRNG->rnorm2(m_nParticipants, 0.0, randA_sd)};
  for (int p = 0; p < m_nParticipants; p++) {
    vvra[p][0] = -randomSlopes[p];
    vvra[p][1] = randomSlopes[p];
  }

  // random slope for B
  double targ_ss_randB = SS_randB / (m_nReps * 2 * 2);
  double randB_sd = sqrt(targ_ss_randB / m_nParticipants);
  std::vector<std::vector<double>> vvrb(m_nParticipants,
					std::vector<double>(2, 0.0));
  randomSlopes = m_pRNG->rnorm2(m_nParticipants, 0.0, randB_sd);
  for (int p = 0; p < m_nParticipants; p++) {
    vvrb[p][0] = -randomSlopes[p];
    vvrb[p][1] = randomSlopes[p];
  }

  // random slope for AB
  double targ_ss_randAB = SS_randAB / (m_nReps * 2 * 2);
  double randAB_sd = sqrt(targ_ss_randAB / m_nParticipants);
  std::vector<std::vector<std::vector<double>>> vvrab;
  randomSlopes = m_pRNG->rnorm2(m_nParticipants, 0.0, randAB_sd);
    
  for (int p = 0; p < m_nParticipants; p++) {
    vvrab.push_back(std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0)));
    vvrab[p][0][0] = randomSlopes[p];
    vvrab[p][0][1] = -randomSlopes[p];
    vvrab[p][1][0] = -randomSlopes[p];
    vvrab[p][1][1] = randomSlopes[p];
  }
  
  // errors
  std::vector<std::vector<double>> vErr{errors(estr, sqrt(.1))};

  // apply linear model to calculate the DV
  for (int p = 0; p < m_nParticipants; p++) {
    for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) {
	for (int i = 0; i < m_nReps; i++) {
	  m_pData[p][r][c][i] = mu + va[r] + vb[c] + vab[r][c] + // fixed effects
	    vri[p] + vvra[p][r] + vvrb[p][c] + vvrab[p][r][c] + // random effects
	    vErr[p][m_pOrder[p][r][c][i]]; // error
	}
      }
    }
  }
}

anovaStats2x2 SimData2x2::anova() {
  anovaStats2x2 astats;

  double grandSum{0.0};
  std::vector<double> vRowEff(2, 0.0);
  std::vector<double> vColEff(2, 0.0);
  std::vector<std::vector<double>> vRCEff(2, std::vector<double>(2, 0.0));
  std::vector<double> vPEff(m_nParticipants, 0.0);
  std::vector<std::vector<double>>
    vPRowEff(m_nParticipants, std::vector<double>(2, 0.0));
  std::vector<std::vector<double>>
    vPColEff(m_nParticipants, std::vector<double>(2, 0.0));
  std::vector<std::vector<std::vector<double>>>
    vPRCEff(m_nParticipants,
	    std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0)));

  for (int p = 0; p < m_nParticipants; p++) {
    for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) {
	for (int i = 0; i < m_nReps; i++) {
	  grandSum += m_pData[p][r][c][i];
	  vRowEff[r] += m_pData[p][r][c][i];
	  vColEff[c] += m_pData[p][r][c][i];
	  vRCEff[r][c] += m_pData[p][r][c][i];
	  vPEff[p] += m_pData[p][r][c][i];
	  vPRowEff[p][r] += m_pData[p][r][c][i];
	  vPColEff[p][c] += m_pData[p][r][c][i];
	  vPRCEff[p][r][c] += m_pData[p][r][c][i];
	}
	// finished one cell for participant p
	vPRCEff[p][r][c] = vPRCEff[p][r][c] / m_nReps;
	if (r == 1) { // finished one column for participant p
	  vPColEff[p][c] = vPColEff[p][c] / (2 * m_nReps);
	}
      }
      // finished one row for one participant
      vPRowEff[p][r] = vPRowEff[p][r] / (2 * m_nReps);
    } // finished participant p
    vPEff[p] = vPEff[p] / (2 * 2 * m_nReps);
  }

  astats.grandMean = grandSum / (m_nParticipants * 2 * 2 * m_nReps);

  // calculate main effects  
  for (size_t r = 0; r < 2; r++) {
    vRowEff[r] = vRowEff[r] / (m_nParticipants * m_nReps * 2) - astats.grandMean;
    astats.SSRow += vRowEff[r] * vRowEff[r];
  }

  for (size_t c = 0; c < 2; c++) {
    vColEff[c] = vColEff[c] / (m_nParticipants * m_nReps * 2) - astats.grandMean;
    astats.SSCol += vColEff[c] * vColEff[c];
  }
  
  for (size_t r = 0; r < 2; r++) {
    for (size_t c = 0; c < 2; c++) {
      vRCEff[r][c] = vRCEff[r][c] / (m_nParticipants * m_nReps) -
	vRowEff[r] - vColEff[c] - astats.grandMean;
      astats.SSRC += vRCEff[r][c] * vRCEff[r][c];
    }
  }

  // finish off participant effects
  for (int p = 0; p < m_nParticipants; p++) {
    vPEff[p] = vPEff[p] - astats.grandMean;
    for (int r = 0; r < 2; r++) {
      vPRowEff[p][r] = vPRowEff[p][r] - vRowEff[r] - vPEff[p] - astats.grandMean;
      astats.SSPRow += vPRowEff[p][r] * vPRowEff[p][r];
    }
    for (int c = 0; c < 2; c++) {
      vPColEff[p][c] = vPColEff[p][c] - vColEff[c] - vPEff[p] - astats.grandMean;
      astats.SSPCol += vPColEff[p][c] * vPColEff[p][c];
    }
    for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) {
	vPRCEff[p][r][c] = vPRCEff[p][r][c] - vRowEff[r] - vColEff[c] - vRCEff[r][c] -
	  vPRowEff[p][r] - vPColEff[p][c] - vPEff[p] - astats.grandMean;
	astats.SSPRC += vPRCEff[p][r][c] * vPRCEff[p][r][c];
      }
    }
  }

  // calculate sums of squares, mean squares, and F-ratios
  astats.SSRow = astats.SSRow * m_nParticipants * NTrials() / 2;
  astats.SSCol = astats.SSCol * m_nParticipants * NTrials() / 2;
  astats.SSRC = astats.SSRC * m_nParticipants * NTrials() / 4;
  astats.SSPRow = astats.SSPRow * NTrials() / 2;
  astats.SSPCol = astats.SSPCol * NTrials() / 2;
  astats.SSPRC = astats.SSPRC * NTrials() / 4;
  astats.F_Row = astats.SSRow / (astats.SSPRow / (m_nParticipants - 1));
  astats.F_Col = astats.SSCol / (astats.SSPCol / (m_nParticipants - 1));
  astats.F_RC = astats.SSRC / (astats.SSPRC / (m_nParticipants - 1));
  astats.p_row = 1.0 - gsl_cdf_fdist_P(astats.F_Row, 1, m_nParticipants - 1);
  astats.p_col = 1.0 - gsl_cdf_fdist_P(astats.F_Col, 1, m_nParticipants - 1);
  astats.p_rc = 1.0 - gsl_cdf_fdist_P(astats.F_RC, 1, m_nParticipants - 1);

  // uncomment to check calculations
  // std::cout << astats << std::endl;
  
  return astats;
}

std::vector<int> possibleSubblocks(const int& nReps) {
  std::vector<int> vResult{};
  for (int i = 1; i <= nReps; i++) {
    if ((nReps % i) == 0) {
      vResult.push_back(i);
    }
  }
  return vResult;
}

void SimData2x2::randomize(const int& nSubblocks /* = 1 */) {
  const int sbSize = (2 * 2 * m_nReps) / nSubblocks;

  std::vector<int> v(2 * 2 * m_nReps, 0);
  std::iota(v.begin(), v.end(), 0);
  
  for (int p = 0; p < m_nParticipants; p++) {
    m_pRNG->permuteSubblocks(v, nSubblocks);

    for (int s = 0; s < nSubblocks; s++) {
      int ix = s * sbSize;
      for (int se = 0; se < sbSize; se++) {
	int k = se / (sbSize / (2 * 2)); // which level (there are 4)
	int r = k / 2;  // which row
	int c = k % 2;  // which column
	int i = se % (sbSize / (2 * 2)) + (sbSize / (2 * 2)) * s; // which rep
	m_pOrder[p][r][c][i] = v[ix + se];
      }
    }
  }

}

void SimData2x2::walkInternal(const std::array<int, 4>& an_cellSeq) {
  for (int p = 0; p < m_nParticipants; p++) {
    int n_startingPos = p % 4;
    for (int i = 0; i < NTrials(); i++) {
      int n_thisCell = an_cellSeq[(n_startingPos + i) % 4];
      int n_r = n_thisCell / 2;
      int n_c = n_thisCell % 2;
      m_pOrder[p][n_r][n_c][i / 4] = i;
    }
  }
}

void SimData2x2::circularWalk() {
  // simple 'circular' walk of conditions, with starting position counterbalanced
  // across participants
  // P1: R1C1 -> R1C2 -> R2C2 -> R2C1 -> ...
  // P2: R1C2 -> R2C2 -> R2C1 -> R1C1 -> ...
  // P3: R2C1 -> ... etc
  std::array<int, 4> an_cellSeq{0, 1, 3, 2};
  walkInternal(an_cellSeq);
}

void SimData2x2::fig8Walk() {
  // 'z-shaped' walk of conditions, with starting position counterbalanced
  // across participants
  // P1: R1C1 -> R1C2 -> R2C1 -> R2C2 -> ...
  // P2: R1C2 -> R2C1 -> R2C2 -> R1C1 -> ...
  // P3: R2C1 -> ... etc
  std::array<int, 4> an_cellSeq{0, 1, 2, 3};
  walkInternal(an_cellSeq);
}

void SimData2x2::writeCSV(const std::string& fname) const {
  auto strNP = std::to_string(m_nParticipants);
  
  std::ofstream fout(fname);
  fout << "id,A,B,tnum,dv" << std::endl;
  
  for (int p = 0; p < m_nParticipants; p++) {
    auto strID = std::to_string(p + 1);
    if (strID.size() < strNP.size()) {
      strID = std::string(strNP.size() - strID.size(), '0') + strID;
    }
    for (int r = 0; r < 2; r++) {
      auto strCondA = std::string(1, 'A') + std::to_string(r + 1);
      for (int c = 0; c < 2; c++) {
	auto strCondB = std::string(1, 'B') + std::to_string(c + 1);
	for (int i = 0; i < m_nReps; i++) {
	  fout << "\"P" << strID << "\",\""
	       << strCondA << "\",\""
	       << strCondB << "\","
	       << (m_pOrder[p][r][c][i] + 1) << ","
	       << m_pData[p][r][c][i] << std::endl;
	}
      }
    }
  }
}

void SimData2x2::run(const int& nmc,
		     const double& eta2_row,
		     const double& eta2_col,
		     const double& eta2_rc,
		     const enum errstr& estr,
		     const enum facRandStrategy& randStrategy,
		     const int& nSubblocks) {
  
  m_vSignificant = {0, 0, 0};
  anovaStats2x2 result;
  
  for (int i = 0; i < nmc; i++) {
    generate(0.0, eta2_row, eta2_col, eta2_rc,
	     .11, .11, .11, .11,
	     estr, randStrategy, nSubblocks);
    result = anova();
    if (result.p_row < .05) {
      m_vSignificant[0]++;
    }
    if (result.p_col < .05) {
      m_vSignificant[1]++;
    }
    if (result.p_rc < .05) {
      m_vSignificant[2]++;
    }
  }
}

std::tuple<std::vector<std::string>,
	   std::vector<std::string>,
	   std::vector<int>, std::vector<double>> SimDataOneFactor::getData() {

  std::vector<std::string> vID;
  std::vector<std::string> vCond;
  std::vector<int> vOrder;
  std::vector<double> vDV;

  auto strNP = std::to_string(m_nParticipants);
  
  for (int p = 0; p < m_nParticipants; p++) {
    auto strID = std::to_string(p + 1);
    strID = "P" + std::string(strNP.size() - strID.size(), '0') + strID;
    for (int k = 0; k < m_nLevels; k++) {
      auto strCond = std::string(1, 'A') + std::to_string(k + 1);
      for (int r = 0; r < m_nReps; r++) {
	vID.push_back(strID);
	vCond.push_back(strCond);
	vOrder.push_back(m_pOrder[p][k][r] + 1);
	vDV.push_back(m_pData[p][k][r]);
      }
    }
  }

  return {vID, vCond, vOrder, vDV};
};

std::tuple<std::vector<std::string>,
	   std::vector<std::string>,
	   std::vector<std::string>,
	   std::vector<int>, std::vector<double>> SimData2x2::getData() {

  std::vector<std::string> vID;
  std::vector<std::string> vA;
  std::vector<std::string> vB;
  std::vector<int> vOrder;
  std::vector<double> vDV;

  auto strNP = std::to_string(m_nParticipants);
  
  for (int p = 0; p < m_nParticipants; p++) {
    auto strID = std::to_string(p + 1);
    strID = "P" + std::string(strNP.size() - strID.size(), '0') + strID;
    for (int r = 0; r < 2; r++) {
      auto strCondA = std::string(1, 'A') + std::to_string(r + 1);
      for (int c = 0; c < 2; c++) {
	auto strCondB = std::string(1, 'B') + std::to_string(c + 1);
	for (int i = 0; i < m_nReps; i++) {
	  // fout << "\"P" << strID << "\",\""
	  //      << strCondA << "\",\""
	  //      << strCondB << "\","
	  //      << (m_pOrder[p][r][c][i] + 1) << ","
	  //      << m_pData[p][r][c][i] << std::endl;
	  vID.push_back(strID);
	  vA.push_back(strCondA);
	  vB.push_back(strCondB);
	  vOrder.push_back(m_pOrder[p][r][c][i] + 1);
	  vDV.push_back(m_pData[p][r][c][i]);	  
	}
      }
    }
  }

  return {vID, vA, vB, vOrder, vDV};
};
