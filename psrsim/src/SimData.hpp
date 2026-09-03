#ifndef SIMDATA_INCLUDED
#define SIMDATA_INCLUDED

#include "RNG_base.hpp"

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>

enum errstr {
  noAuto = 1,
  expDecay = 2,
  pinkNoise = 3,
  randomWalk = 4,
  mixed = 5
};

enum facRandStrategy {
  PSR,
  circular,
  fig8
};

class anovaStatsBase {
public:
  double grandMean;
  explicit anovaStatsBase() : grandMean{0.0} {};
  anovaStatsBase(const anovaStatsBase&) = default;
  anovaStatsBase& operator=(const anovaStatsBase& rhs) = default;
};

class anovaStats : public anovaStatsBase {
public:
  std::vector<double> mainEff;
  double ssEff;
  double ssError;
  double FRatio;
  double dfNum;
  double dfDen;
  double pValue;
  explicit anovaStats(int nLevels) :
    anovaStatsBase(),
    mainEff{std::vector<double>(nLevels, 0.0)},
    ssEff{0.0}, ssError{0.0},
    FRatio{0.0}, dfNum{0.0}, dfDen{0.0},
    pValue{0.0} {};
  anovaStats(const anovaStats&) = default;
  anovaStats& operator=(const anovaStats& rhs) = default;

  friend std::ostream& operator<<(std::ostream& os, const anovaStats& as) {
    os << "grand mean : " << as.grandMean << std::endl;
    os << "main eff   : (";
    for (size_t k = 0; k < as.mainEff.size(); k++) {
      os << as.mainEff[k];
      if (k < (as.mainEff.size() - 1)) {
	os << ", ";
      }
    }
    os << ")" << std::endl;
    os << "SS_eff     : " << as.ssEff << std::endl;
    os << "SS_err     : " << as.ssError << std::endl;
    os << "df_num     : " << as.dfNum << std::endl;
    os << "df_den     : " << as.dfDen << std::endl;
    os << "F-ratio    : " << as.FRatio << std::endl;
    os << "p-value    : " << as.pValue << std::endl;
    return os;
  }
};

class anovaStats2x2 : public anovaStatsBase {
public:
  double SSRow;
  double SSCol;
  double SSRC;
  double SSPRow;
  double SSPCol;
  double SSPRC;
  double F_Row;
  double F_Col;
  double F_RC;
  double p_row;
  double p_col;
  double p_rc;
  explicit anovaStats2x2() : anovaStatsBase(),
			     SSRow{0.0}, SSCol{0.0}, SSRC{0.0},
			     SSPRow{0.0}, SSPCol{0.0}, SSPRC{0.0},
			     F_Row{0.0}, F_Col{0.0}, F_RC{0.0},
			     p_row{0.0}, p_col{0.0}, p_rc{0.0} {};
  anovaStats2x2(const anovaStats2x2&) = default;
  anovaStats2x2& operator=(const anovaStats2x2& rhs) = default;
  
  friend std::ostream& operator<<(std::ostream& os, const anovaStats2x2& as) {
    os << "grand mean : " << as.grandMean << std::endl;
    os << "SSRow  : " << as.SSRow << std::endl;
    os << "SSCol  : " << as.SSCol << std::endl;
    os << "SSRC   : " << as.SSRC << std::endl;
    os << "SSPRow : " << as.SSPRow << std::endl;
    os << "SSPCol : " << as.SSPCol << std::endl;
    os << "SSPRC  : " << as.SSPRC << std::endl;
    os << "F_Row  : " << as.F_Row << ", p = " << as.p_row << std::endl;
    os << "F_Col  : " << as.F_Col << ", p = " << as.p_col << std::endl;
    os << "F_RC  : " << as.F_RC << ", p = " << as.p_rc << std::endl;
    return os;
  }
};

class SimDataBase {
protected:
  int m_nParticipants;
  int m_nReps;
  std::unique_ptr<RNG_base> m_pRNG;
public:
  SimDataBase(const int& nParticipants,
	      const int& nReps,
	      std::unique_ptr<RNG_base> pRNG);
  virtual ~SimDataBase();

  // void seed(unsigned long value);
  virtual int NTrials() const = 0;
  
  double sumOfSquares(const std::vector<double>& vd) const;
  std::vector<std::vector<double>> errors(const enum errstr& estr,
					  const double& sdExtraNoise = 0.0);
};

class SimDataOneFactor : public SimDataBase {
private:
  int m_nLevels;
  double *** m_pData;
  int *** m_pOrder;
  int m_nSignificant;
  std::vector<double> m_vdPValues;
  
public:
  SimDataOneFactor(const int& nParticipants,
		   const int& nLevels,
		   const int& nReps,
		   std::unique_ptr<RNG_base> pRNG);
  
  virtual ~SimDataOneFactor();

  inline int getNSig() { return m_nSignificant; };
  void writePValues(std::ofstream& out);

  inline int NTrials() const override {return m_nLevels * m_nReps; };

  std::vector<double> MainEffOneFactor(const double& targSS = 1);
  
  void generate(const double& mu,
		const double& eta,
		const double& prop_rint,
		const double& prop_rslp,
		const enum errstr& estr,
		const int& nSubblocks);

  void generate2(const double& mu,
		 const double& eta,
		 const double& prop_rint,
		 const double& prop_rslp,
		 const enum errstr& estr_int,
		 const enum errstr& estr_slp,
		 const int& nSubblocks);
  
  virtual void writeCSV(const std::string& fname) const;

  anovaStats anova();

  void run(const int& nmc, const double& eta2,
	   const enum errstr& estr, const int& nSubblocks);

  void run2(const int& nmc, // no. Monte Carlo runs
	    const double& eta2, // eta-squared
	    const enum errstr& estr_int, // err structure icept
	    const enum errstr& estr_slp, // err structure slope
	    const int& nSubblocks,
	    const double& dMu = 0.0, // intercept (grand mean)
	    const double& prop_rint = 0.35,
	    const double& prop_rslp = 0.11);
  
  void randomize(const int& nSubblocks = 1);

  std::tuple<std::vector<std::string>,
	     std::vector<std::string>,
	     std::vector<int>, std::vector<double>> getData();
};

class SimData2x2 : public SimDataBase {
private:
  double **** m_pData;
  int **** m_pOrder;
  std::vector<int> m_vSignificant;

  void walkInternal(const std::array<int, 4>& an_cellSeq);
  
public:

  SimData2x2(const int& nParticipants,
	     const int& nReps,
	     std::unique_ptr<RNG_base> pRNG);
  virtual ~SimData2x2() override;

  inline int NTrials() const override { return 2 * 2 * m_nReps; };
  inline int getNSig(int ix) const { return m_vSignificant[ix]; };
  
  void generate(const double& mu,
		const double& eta_row,
		const double& eta_col,
		const double& eta_rcx,
		const double& prop_rint,
		const double& prop_rslp_row,
		const double& prop_rslp_col,
		const double& prop_rslp_rcx,
		const enum errstr& estr,
		const enum facRandStrategy& randStrategy,
		const int& nSubblocks);

  void run(const int& nmc,
	   const double& eta2_row,
	   const double& eta2_col,
	   const double& eta2_rc,
	   const enum errstr& estr,
	   const enum facRandStrategy& randStrategy,
	   const int& nSubblocks);
  
  virtual void randomize(const int& nSubblocks = 1);
  void circularWalk();
  void fig8Walk();
  virtual void writeCSV(const std::string& fname) const;

  std::tuple<std::vector<std::string>,
	     std::vector<std::string>,
	     std::vector<std::string>,
	     std::vector<int>, std::vector<double>> getData();
  
  anovaStats2x2 anova();
};

// HELPERS
std::vector<int> possibleSubblocks(const int& nReps);

#endif
