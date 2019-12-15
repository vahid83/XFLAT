/*
 * Phy_SA.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#include <stdlib.h>
#include <new>
#include <cstring>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Fenergy.h"
#include "Fio.h"
#include "Nmr.h"

#include "Phy.h"

#include <stdio.h>

namespace nbgrp {

namespace phy_sa {

using namespace fed;
using namespace nbm;

///===================== Variables ========================

/// Starting and ending points for beams, and total #of theta bins (theta=1)
int start_beam, end_beam, abins;

/// Holds the lenght of the nubeam array
int beam_len;

/// Holds the current radius, this array's length is equal to number of middle points in Nmr mod.
double* rarr;
/// Holds the current delta radius, this array's length is equal to number of middle points in Nmr mod.
double* drarr;

/// cache for delta-r array -- one per theta angle
double* delta_l;

/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

double* L_E;	///< the result of L / <E> for each particle

const int Ndims = 4;    ///< [r,partcl,cmpn,ebin]
size_t start[Ndims], count[Ndims];
int last_idx;           ///< to prevent recalculating all of the indecis when only component index is changed!

///=================== Implementation =====================
/**
 * Returns the length of the NBeam object array for the current node
 * @retval beam array len. in the module
 */
int beamLen()	{ return beam_len; }

///----------------------------------------------
/**
 * Return the number of data dimentions
 * @retval count_dim[] array
 */
int getDim()    { return Ndims; }

///----------------------------------------------
/**
 * Returns naming info for the geometry -- to be used by IO module
 * @param str[] the array of strings indicating the name of each dimention
 */
void getDimInfo(std::string str[])
{
        str[0] = "r";
        str[1] = "prtcl";
        str[2] = "comp";
        str[3] = "ebin";
}

///----------------------------------------------
/**
 * Returns an array of starting points for data dimentions -- to be used by IO
 * @return the start[]
 */
size_t* startDim()      { return start; }

///----------------------------------------------
/**
 * Returns an array containing the length of each dimention -- to be used by IO
 * @retval the count[]
 */
size_t* countDim()      { return count; }

///----------------------------------------------
/**
 * The setter and getter for the beam index that the current node start with
 * @retval the reference to starting beam index
 */
int& startBeamIdx()     { return start_beam; }

///----------------------------------------------
/**
 * The setter and getter for the beam index that the current node end with
 * @retval the reference to ending beam index
 */
int& endBeamIdx()       { return end_beam; }

///----------------------------------------------
/**
 * Returns the dimention length on which the job will be distributed over nodes
 * To be used in NBGroup module
 * @retval the length to the first dimention
 */
int& firstDimLen()       { return Util::Abins(); }

///----------------------------------------------
void calcAngleBins(const double r, const int step_num) throw()
{
	rarr[step_num] = r;
}

///----------------------------------------------
/**
 * Returns the delta-l for a point
 * @pnt the step point at which the result is returned
 * @idx the current index of the angle bin
 * @return the dl
 */
inline double getDeltaL(const int pnt, const int tet) throw()
{
        return delta_l[pnt];
}

///----------------------------------------------
/**
 * Calculates dl for each angle bin
 * @dr the current delta radius
 * @cur_pnt the current point on which dls are calculated
 * @s_pnt the start point -- source point
 * @e_pnt the end point -- destination point
 */
void calcDeltaLs(const double dr, const int cur_pnt, const int s_pnt, const int e_pnt) throw()
{
	delta_l[cur_pnt] = dr;
}

///----------------------------------------------
void init()
{
	if (Util::isBench()) ///< if the code is running in benchmark mode
        {
                abins = 1;
                start_beam = 0;
                end_beam = abins-1;
                printf("Bnechmarking using %d trajectory beams...\n", abins);
        }
        else    /// start_beam and end_beam are set before!
        {
                abins = Util::Abins();
        }

	last_idx = -1; /// beam index is never neg.

	/// Hvv integral coeff init.
	hvv_coeff  = (Util::SQRT2 * Util::GF) / (2. * Util::pi * Util::Rv() * Util::Rv());
	integ_coef = hvv_coeff ;//* dE;

	/// Needed by IO module -- should be done before nmr::init->initBeams
	start[0] = start[1] = start[2] = start[3] = 0;
        count[0] = 1;
        count[1] = nbm_prtcls;
        count[2] = nbm_cmpn;
        count[3] = nbm_ebins;
	
	beam_len = nbm_prtcls;
	int pnts_num = nmr::init(beam_len);
	posix_memalign((void**)&rarr, ALIGN_LEN, pnts_num*sizeof(double));

	posix_memalign((void**)&L_E,            ALIGN_LEN, nbm_prtcls*sizeof(double));
	posix_memalign((void**)&delta_l,	ALIGN_LEN, pnts_num*sizeof(double));
	for (int i = 0; i < nbm_prtcls; ++i)
		L_E[i] = Lv(i)  / Ev(i);

}

///----------------------------------------------
void freemem()
{
	free(rarr);
	nmr::freemem();

	free(L_E);
}

///----------------------------------------------
/**
 * Allocates memory for the hamiltonian array -- to be used by Numeric module
 * @hvv the pointer to the allocated array
 */
void newHvv(double*& hvv)
{
        posix_memalign((void**)&hvv, ALIGN_LEN, nbm_size*sizeof(double));
}

///----------------------------------------------
/**
 * Deletes the memory of hamiltonian array -- to be used by Num. module
 * @hvv array handle
 */
void deleteHvv(double*& hvv)
{
        free(hvv);
}

///----------------------------------------------
/**
 * Initializes an array of NBeam objects by calling the constructor
 * @beam array of neutrino beams to be init.
 */
void initBeam(NBeam* beam)
{
	for (int j = 0; j < nbm_prtcls; ++j)
	{
		/// Call the constructor!
		new (&beam[j]) NBeam(j);
	}	

	if (Util::inVals() == NULL || Util::isBench()) return;
        fio::fillInitData(beam);	 ///< if a data file is provided for resumption

}

///----------------------------------------------
/**
 * Call destructors of NBeam objects
 * @beam the pointer to the NBeam array
 */
void freeBeam(NBeam* beam)
{
	for (int j = 0; j < nbm_prtcls; ++j)
	{
		/// Call destructor!
		beam[j].~NBeam();
	}
}

///----------------------------------------------
inline double D(double r)
{
	double R = Util::Rv() / r;

	double tm = 1. - sqrt( 1. - (R*R) );
	return 0.5 * (tm*tm);
}

///----------------------------------------------
void calc_Hvv(const NBeam* beam, const int step_num, double *REST hvv) throw()
{
	memset(hvv, 0., nbm_size*sizeof(double));

	for (int n = 0; n < nbm_flvs; ++n)
	{
		int neu  = n  << 1;	///< even indices are neu
		int aneu = neu + 1;	///< odds are anti-neu

		Res_t  res_neu = {};
		Res_t res_aneu = {};

		beam[neu].getESum(res_neu);
		beam[aneu].getESum(res_aneu);

		double neu_coeff  = L_E[neu];//Lv(neu)  / Ev(neu);
		double aneu_coeff = L_E[aneu];//Lv(aneu) / Ev(aneu);

		upd_nu_coef(res_neu, res_aneu, neu_coeff, aneu_coeff, hvv);
	}

	double coeff = integ_coef  * D(rarr[step_num]);
	for (int i = 0; i < nbm_size; ++i)	hvv[i] *= coeff;
}

void evolve(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
	double dl = getDeltaL(pnt,0);
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int prc = 0; prc < nbm_prtcls; ++prc)
	{
		obeam[prc].evolveBins(ibeam[prc], prc, dl, hvv, hmatt);
	}
}

void evolveHvv(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
        double dl = getDeltaL(pnt,0);
#ifdef OMP
#       pragma omp parallel for
#endif
        for (int prc = 0; prc < nbm_prtcls; ++prc)
        {
                obeam[prc].evolveBinsHvv(ibeam[prc], prc, dl, hvv, hmatt);
        }
}

void evolveAvg(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg) throw()
{
        double dl = getDeltaL(pnt,0);
#ifdef OMP
#       pragma omp parallel for
#endif
        for (int prc = 0; prc < nbm_prtcls; ++prc)
        {
		obeam[prc].evolveBinsAvg(ibeam[prc], prc, dl, hvv, hmatt, obeamAvg[prc]);
        }
}

void evolveAvg(const int pnt, const double *REST hvv, const double hmatt, NBeam *REST iobeam, NBeam *REST obeamAvg) throw()
{
        double dl = getDeltaL(pnt,0);
#ifdef OMP
#       pragma omp parallel for
#endif
        for (int prc = 0; prc < nbm_prtcls; ++prc)
        {
		iobeam[prc].evolveBinsAvg(prc, dl, hvv, hmatt, obeamAvg[prc]);
        }
}

///----------------------------------------------
/**
 * Takes the average of two beams
 * @ibeam the first beam
 * @obeam the second beam -- the average stores in this beam
 */
void avgBeam(const NBeam *REST ibeam, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int prc = 0; prc < nbm_prtcls; ++prc)
	{
		obeam[prc].addAvg(ibeam[prc]);
	}
}

} /// Phy namespace

} /// NBGroup namespace

