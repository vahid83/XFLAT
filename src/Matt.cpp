/*
 * Matt.cpp
 *
 *  Created on: Jul 20, 2013
 *      Author: vahid
 */

#include <cmath>
#include "Util.h"
#include "Matt.h"

namespace matt {

double nb0, Rv, hNS;

double mat_coeff, nbfar_coeff;

double& matt_nb0()	{ return nb0; } //!< To be set from Parser
double& matt_Rv()	{ return Rv;  } //!< To be set from Parser
double& matt_hNS()	{ return hNS; } //!< To be set from Parser

void init()
{
	mat_coeff	= Util::SQRT2 * Util::GF * Util::Ye();

	double M	= Util::Mns() / 1.4;
	double s	= 100. / Util::S();
	nbfar_coeff 	= 4.2E+45 * Util::gs() * (M*M*M) * (s*s*s*s);
}

inline double nb_far(double r)
{
	double R = 10. / r;

	return nbfar_coeff * (R*R*R);	///< 4.2E+30 x 1E+15 to convert to 1/Km3
}

inline double nb_near(double r)
{
	return nb0 * exp((Rv - r) / hNS);
}

inline double nb(double r)
{
	return nb_far(r) + nb_near(r);
}

void getHmatt(	const double radius,		///< input
		double ret_r[4] ) throw()	///< output
{
	double cf = mat_coeff * nb(radius);
	ret_r[0] = cf;
	ret_r[1] = 0.;
	ret_r[2] = 0.;
	ret_r[3] = -cf;
}

void getHm( const double radius, double& ret_r ) throw()
{
	if (!Util::hasMatter())	///< No matter profile is included
		ret_r = 0;
	else
		ret_r = mat_coeff * nb(radius) * .5; ///< there is a 1/2 Hamiltonian factor
}

}


