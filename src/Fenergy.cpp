/*
 * Fenergy.cpp
 *
 *  Created on: Sep 12, 2013
 *      Author: vahid
 */

#include "Fenergy.h"
#include "Util.h"
#include "NBeam.h"
namespace fed {

double dE;

/// This is the result of the integral for F2 and eta = 3
inline double fk(double eta = 3., int k = 2) throw()
{ return 18.9686; }

double fq(double q, double Tv, double eta_v) throw()	
{
	//return 1. / dE;
	return (1. / fk()) * (1. / (Tv*Tv*Tv)) * ( (q*q) / (exp(q/Tv - eta_v)+1.) );
}

double avgEnergy(int prtcl) throw()
{
        const int Eb = 800;//4096;///< high resolution integral!
        double de = (Util::E1() - Util::E0()) / Eb;
        double frst = Util::E0() + de*.5;
        double avg = 0., Eavg = 0.;
        for (int e = 0; e < Eb; ++e)
        {
                double q = frst + e * de;
                double fe = fq(q, nbm::Tv(prtcl), nbm::eta(prtcl));

                avg += fe;// * de;
                Eavg += fe * q;// * de;
        }
	//printf("The Average energy for neutrino# %d is %f MeV\n",prtcl,  Eavg / avg / MEV);
        return Eavg / avg;
}


void init(int ebins)
{
	dE = (Util::E1() - Util::E0()) / ebins;
}

}


