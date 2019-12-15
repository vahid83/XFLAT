#pragma once
#include "global.h"

#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include <string>

#include <sys/time.h>

#include "NBGroup.h"

namespace Util
{
	const double SQRT2 = 1.41421356237;                 ///< SQRT2
	const double pi    = 3.14159265359;                 ///< PI
	const double Pi    = 3.14159265359;                 ///< PI
	const double Pix2  = 6.28318530718;                 ///< PIx2
	const double pi2   = 1.57079632679;                 ///< PI/2
	const double Pih   = 1.57079632679;                 ///< PI/2
	const double e     = 2.71828182846;                 ///< e
	const double GF    = 1.16637E-11 * MEV_2;/*MeV^-2*/	///< Fermi coupling constant

	bool& isBench();
	int& newFilestep();
	int& syncstep();

	double& rstep1();
	int& tstep1();
	int& itrstep1();

	double& rstep2();
	int& tstep2();
	int& itrstep2();

	int& startBeam();
	int& endBeam();

	int& multiNodeBench();
        int& minNodes();
        int& maxNodes();

        int& hasMatter();

	double& eps0();
	double& ksi();
	double& mu();

	int& Ts();
	int& Tn();
	double& R0();
	double& Rn();
	double& dr();
	double& maxDr();
	double& E0();
	double& E1();

	int& SPoints();
	int& Pbins();
	int& Abins();
	int& Ebins();
	int& Flvs();

	double& Ye();
	double& Rv();
	double& Mns();
	double& gs();
	double& S();
	double& hNS();
	double& nb0();

	double& theta();
	double& theta2();
	double& dm2();

	std::string& filePrefix();
	int& dumpMode();
	char* inVals();
	void inVals(char* s);

/// ************************ General helper methods ****************************

	/*!
	 * Calculates multiplication of two complex numbers
	 */
	inline void mul_cmplx(const double ar, const double ai, const double br, const double bi, double& res_r, double& res_i)
	{
		res_r = (ar * br) - (ai * bi);
		res_i = (ar * bi) + (ai * br);
	}
	
	/*!
	 * Calculates norm2 of a and b
	 */
	inline double norm2(const double a, const double b) { return a*a + b*b; }

	double benchmark();

	void beginBenchmark();
	void endBenchmark();
} /// End of namespace

