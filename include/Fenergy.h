/*
 * Fenergy.h
 *
 *  Created on: Sep 12, 2013
 *      Author: vahid
 */

#ifndef FENERGY_H_
#define FENERGY_H_

namespace fed {

extern double dE;

/**
 * init. energy bins and dE
 */
void init(int ebins);

/**
 * Energy spectra function
 * @q enery bin
 * @Tv temprature
 * @eta_v we take it as 3
 * @return initial energy spectra for neutrinos
 */
double fq(double q, double Tv, double eta_v) throw();

/**
 * Calculates the average energy for each particle
 * @prtcl the particle index
 * @return the average energy
 */
double avgEnergy(int prtcl) throw();
}


#endif /* FENERGY_H_ */
