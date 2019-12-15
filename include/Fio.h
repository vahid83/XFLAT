/*
 * Fio.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifndef FIO_H_
#define FIO_H_

#include "NBGroup.h"
#include "NBeam.h"

namespace nbgrp {

namespace IO_MODULE {

/**
 * Initialize netcdf handles
 * @file_counter indicates the current # of generated files
 */
void init(int file_counter = 0);

/**
 * Free caches
 */
void freemem();

/**
 * if the input file provided, this function handle that!
 * @nubeam the beam array
 */
void fillInitData(nbm::NBeam *REST nubeam);

/**
 * Based on verbose mode decides which function to be called
 */
void dumpToFile(const nbm::NBeam *REST nubeam, const int itr, const double r, const double dr = 1, int theta = 0, int phi = 0) throw();

} /// Fio namespace

namespace fio=IO_MODULE;

} /// NBGroup namespace


#endif /* FIO_H_ */
