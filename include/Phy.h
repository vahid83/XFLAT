/*
 * Phy.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifndef PHY_H_
#define PHY_H_

#include <string>

#include "NBeam.h"
#include "Fio.h"

namespace nbgrp {

namespace PHY_MODULE {

/**
 * Returns the length of the NBeam object array for the current node
 * @retval beam array len. in the module
 */
int beamLen();

/**
 * Return the number of data dimentions
 * @retval count_dim[] array
 */
int getDim();

/**
 * Returns naming info for the geometry -- to be used by IO module
 * @param str[] the array of strings indicating the name of each dimention
 */
void getDimInfo(std::string str[]);

/**
 * Returns an array of starting points for data dimentions -- to be used by IO
 * @return the start[]
 */
size_t* startDim();

/**
 * Returns an array containing the length of each dimention -- to be used by IO
 * @retval the count[]
 */
size_t* countDim();

/**
 * The setter and getter for the beam index that the current node start with
 * @retval the reference to starting beam index
 */
int& startBeamIdx();

/**
 * The setter and getter for the beam index that the current node end with
 * @retval the reference to ending beam index
 */
int& endBeamIdx();

/**
 * Returns the dimention length on which the job will be distributed over nodes
 * To be used in NBGroup module
 * @retval the length to the first dimention
 */
int& firstDimLen();

/**
 * Initializes ithe module and allocates memory for arrays
 */
void init();

/**
 * freeups allocated arrays
 */
void freemem();

/**
 * Allocates memory for the hamiltonian array -- to be used by Numeric module
 * @hvv the pointer to the allocated array
 */
void newHvv(double*& hvv);

/**
 * Deletes the memory of hamiltonian array -- to be used by Num. module
 * @hvv array handle
 */
void deleteHvv(double*& hvv);

/**
 * Initializes an array of NBeam objects by calling the constructor
 * @beam array of neutrino beams to be init.
 */
void initBeam(nbm::NBeam* beam);

/**
 * Call destructors of NBeam objects
 * @beam the pointer to the NBeam array
 */
void freeBeam(nbm::NBeam* beam);

/**
 * Calculates and initializes angle bins
 * @r the current radius at which the bins are computed
 * @step_num indicates the current step number (first, mid, or last poiint)
 */
void calcAngleBins(const double r, const int step_num) throw();

/**
 * Calculates dl for each angle bin
 * @dr the current delta radius
 * @cur_pnt the current point on which dls are calculated
 * @s_pnt the start point -- source point
 * @e_pnt the end point -- destination point
 */
void calcDeltaLs(const double dr, const int cur_pnt, const int s_pnt, const int e_pnt) throw();

/**
 * Calculates the hamiltonians 
 * @beam the NBeam array on which the Hamiltonian is calculated
 * @step_num the point at which the Hamiltonian is calculated
 * @hvv the result array
 */
void calc_Hvv(const nbm::NBeam *REST beam, const int step_num, double *REST hvv) throw();

/**
 * Calculates the hamiltonian for the point for two points simultaneously
 * @beam1 the first NBeam array on which the first Hamiltonian is calculated
 * @step_num1 the first point at which the first Hamiltonian is calculated
 * @hvv1 the first Hamiltonian result array
 * @beam2 the second NBeam array on which the second Hamiltonian is calculated
 * @step_num2 the second point at which the second Hamiltonian is calculated
 * @hvv2 the second Hamiltonian result array
 *
 */
void calc_Hvv(const nbm::NBeam *REST beam1, const int step_num1, double *REST hvv1, const nbm::NBeam *REST beam2, const int step_num2, double *REST hvv2) throw();

/**
 * Calls the evolution function of NBeam class
 * @ibeam the first beam
 * @pnt the point at which the evolution is done
 * @hvv the hamiltonian
 * @hmatt the matter
 * @obeam the result beam
 */
void evolve(const nbm::NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST obeam) throw();
void evolveHvv(const nbm::NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST obeam) throw();
void evolveAvg(const nbm::NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST obeam, nbm::NBeam *REST obeamAvg) throw();
void evolveHvvAvg(const nbm::NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST obeam, nbm::NBeam *REST obeamAvg) throw();
void evolveAvgErr(const nbm::NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST obeam, nbm::NBeam *REST obeamAvg, nbm::NBeam *REST obeamErr) throw();
void evolveAvg(const int pnt, const double *REST hvv, const double hmatt, nbm::NBeam *REST iobeam, nbm::NBeam *REST obeamAvg) throw();

/**
 * Takes the average of two beams
 * @ibeam the first beam
 * @obeam the second beam -- the average stores in this beam
 */
void avgBeam(const nbm::NBeam *REST ibeam, nbm::NBeam *REST obeam) throw();

}

namespace phy=PHY_MODULE;

} /// NBGroup namespace

#endif /* PHY_H_ */
