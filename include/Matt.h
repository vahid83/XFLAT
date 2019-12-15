/*
 * Matt.h
 *
 *  Created on: Jul 17, 2013
 *      Author: vahid
 */

#ifndef MATT_H_
#define MATT_H_


namespace matt {

double& matt_nb0();
double& matt_Rv();
double& matt_hNS();

/**
 * init. coefficients
 */
void init();

/**
 * returns matter profile
 * @radius current radius
 */
inline double nb(double r);

/**
 * returns matter density
 * @radius current radius
 * @ret_r return value
 */
void getHm( const double radius, double& ret_r ) throw();

}


#endif /* MATT_H_ */
