/*
 * Nmr.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifndef NMR_H_
#define NMR_H_

namespace nbgrp {

namespace nmr {

/**
 * Initializes NBeam's
 * @len total size of bins
 * @return number of step size
 */
int init(int len);

/**
 * free ups memory
 */
void freemem();

/**
 * Evolves neutrinos from R0 to Rn
 * @return the number of computed steps
 */
int evolutionLoop() throw();

} /// Nmr namespace

} /// NBGroup namespace


#endif /* NMR_H_ */
