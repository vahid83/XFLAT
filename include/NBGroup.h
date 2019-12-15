/*
 * NBGroup.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifndef NBGROUP_H_
#define NBGROUP_H_

namespace nbgrp {

/**
 * Initializes all the lower layer modules
 */
void init();

/**
 * Calls into main evolution loop
 */
void particleLoop();

/**
 * Free up memory
 */
void freemem();


}


#endif /* NBGROUP_H_ */
