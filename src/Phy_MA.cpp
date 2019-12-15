/*
 * Phy_MA.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */
#ifdef MMPI
#include <mpi.h>
#endif

#include <cstdio>
#include <stdlib.h>
#include <new>
#include <cstring>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fenergy.h"
#include "Fio.h"

#include "Phy.h"

#include <iostream>

namespace nbgrp {


namespace phy_ma {

using namespace fed;
using namespace nbm;

///===================== Variables ========================

/// Starting and ending points for beams, and total #of theta bins
int start_beam, end_beam, abins, pbins;

int beam_len;///< length of the nubeam array

int tetbins; ///< number of beams on the node based on the computational capability! 

/// place holder for partial and Hvv result
Res_t part_hvv;
Res_t part_hvv_cos;
//double* Hvv;

/// cache for delta-r array -- one per theta angle
double* delta_l;

/// Arrays for manual reduction - omp doesn't support reduction over dynamic arrays


/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

/// Arrays for cos^2(x) and dcos^2(x) bins at surface -- has to be computed once
double* Cos2_R;
double* dCos2_R;

/// Arrays for cos^2(x) and dcos^2(x) bins at 'r' -- recalculated at each step
double* cos_tbin;
double* dcos_tbin;

double* L_E;	///< the result of L / <E> for each particle

const int Ndims = 5;    ///< [r,theta,partcl,cmpn,ebin]
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
        str[1] = "theta";
        str[2] = "prtcl";
        str[3] = "comp";
        str[4] = "ebin";
}

///----------------------------------------------
/**
 * Returns an array of starting points for data dimentions -- to be used by IO
 * @return the start[]
 */
size_t* startDim()	{ return start; }

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
/**
 * Init. cos2 bins of the surface
 */
inline void init_cos2R() throw()
{
/*
///warm up!
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int i = 0; i < abins; ++i)
		Cos2_R[i] = dCos2_R[i] = 0.;
*/
	double dcos_theta = (1. - 0.) / abins;

	/// delta theta -- since we're creating bins on (0-1) uniformly, the distance is fixed
	Cos2_R[0] = 1. - dcos_theta*.5;
	dCos2_R[0] = dcos_theta;//1. - Cos2_R[0];
	double cos_theta = Cos2_R[0] - dcos_theta;

	for (int i = 1; i < abins; ++i, cos_theta -= dcos_theta)
	{
		/// fill up the Cos2 and dCos2 bins
		Cos2_R[i]  = cos_theta;
		dCos2_R[i] = dcos_theta;//Cos2_R[i-1] - Cos2_R[i];
	}
}

///----------------------------------------------
/**
 * Initializes cosine angle bins
 * @r the current radius at which the bins are computed
 * @step_num indicates the current step number (first, mid, or last poiint)
 */
void calcAngleBins(const double r, const int step_num) throw()
{
	double r2 = Util::Rv() / r;
	r2 *= r2;
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		int idx = step_num*tetbins + ang;
		cos_tbin [idx] = sqrt(1. - r2 *  (1. - Cos2_R[start_beam+ang]));
		dcos_tbin[idx] = .5 * (1./cos_tbin[idx]) * r2 * dCos2_R[start_beam+ang];
	}
/*
//#ifdef OMP
//#       pragma omp parallel for
//#endif
        for (int ang = 0; ang < tetbins-1; ++ang)
	{
		int idx = step_num*tetbins + ang;
		dcos_tbin[idx] = cos_tbin [idx+1] - cos_tbin [idx];
	}
	dcos_tbin[step_num*tetbins + tetbins-1] = dcos_tbin[step_num*tetbins + tetbins-2];
*/
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
	return delta_l[pnt*tetbins+tet];
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
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		delta_l[cur_pnt*tetbins+ang] = dr / ( (cos_tbin[s_pnt*tetbins+ang] ));//+ cos_tbin[e_pnt*tetbins+ang])*.5 );
	}
}

///----------------------------------------------
/**
 * Initializes ithe module and allocates memory for arrays
 */
void init()
{
	if (Util::isBench()) ///< if the code is running in benchmark mode
        {
                abins = 1800;
                start_beam = 0;
                end_beam = abins-1;
		printf("Bnechmarking using %d trajectory beams...\n", abins);
        }
        else	/// start_beam and end_beam are set before!
        {
                abins = Util::Abins();
        }

	last_idx = -1; /// beam index is never neg.

        /// start_beam and end_beam are set from NBGroup module
	tetbins = end_beam - start_beam + 1;// abins;

	/// Needed by IO module -- should be done before nmr::init->initBeams
	start[0] = start[1] = start[2] = start[3] = start[4] = 0;
        count[0] = 1;
        count[1] = tetbins;
        count[2] = nbm_prtcls;
        count[3] = nbm_cmpn;
        count[4] = nbm_ebins;

	//posix_memalign((void**)&Hvv, 				ALIGN_LEN, tetbins*nbm_size*sizeof(double));



	posix_memalign((void**)&Cos2_R,    			ALIGN_LEN, abins*sizeof(double));
	posix_memalign((void**)&dCos2_R,   			ALIGN_LEN, abins*sizeof(double));

	/// Init. cosine^2(R0)
	init_cos2R();

	/// Hvv integral coeff init.
	hvv_coeff  = (Util::SQRT2 * Util::GF) / (2. * Util::pi * Util::Rv() * Util::Rv());
	integ_coef = hvv_coeff ;//* dE;

	beam_len = tetbins*nbm_prtcls;
	int pnts_num = nmr::init(beam_len);	///< calls initBeam
	posix_memalign((void**)&cos_tbin,     		ALIGN_LEN, tetbins*pnts_num*sizeof(double));
	posix_memalign((void**)&dcos_tbin,    		ALIGN_LEN, tetbins*pnts_num*sizeof(double));
	posix_memalign((void**)&delta_l,    		ALIGN_LEN, tetbins*pnts_num*sizeof(double));

        posix_memalign((void**)&L_E,            ALIGN_LEN, nbm_prtcls*sizeof(double));
	for (int i = 0; i < nbm_prtcls; ++i)
		L_E[i] = Lv(i) / Ev(i);

}

///----------------------------------------------
/**
 * freeups allocated arrays
 */
void freemem()
{
	nmr::freemem();

	free(cos_tbin);
	free(dcos_tbin);
	free(delta_l);

	free(Cos2_R);
	free(dCos2_R);

	free(L_E);
}

///----------------------------------------------
/**
 * Allocates memory for the hamiltonian array -- to be used by Numeric module
 * @hvv the pointer to the allocated array
 */
void newHvv(double*& hvv)
{
	posix_memalign((void**)&hvv, ALIGN_LEN, tetbins*nbm_size*sizeof(double));
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
#if defined(OMP) && !defined(__MIC__)
#       pragma omp parallel for
#endif
	for (int i = 0; i < tetbins; ++i)
	{
		for (int j = 0; j < nbm_prtcls; ++j)
		{
			/// Call constructor!
			new (&beam[i*nbm_prtcls+j]) NBeam(j);
		}
	}

	if (Util::inVals() == NULL || Util::isBench()) return;
	fio::fillInitData(beam);	///< if a data file is provided for resumption
}

///----------------------------------------------
/**
 * Call destructors of NBeam objects
 * @beam the pointer to the NBeam array
 */
void freeBeam(NBeam* beam)
{
	for (int i = 0; i < tetbins; ++i)
	{
		for (int j = 0; j < nbm_prtcls; ++j)
		{
			/// Call destructor!
			beam[i*nbm_prtcls+j].~NBeam();
		}
	}
}

///----------------------------------------------
/**
 * Calculates Hvv matrix over all bins
 * @Cos cos(x) coeff
 * @partial_nu partial integral result
 * @partial_nu_cos partial integral result with cos(x) term
 * @ret the whole integral result over angle/energy bins for the current theta
 */
inline void getHvv( const double Cos, const double *REST partial_nu, const double *REST partial_nu_cos,	///< inputs
				double *REST ret) throw()																///< output
{
	for (int z = 0; z < nbm_size; ++z)
		ret[z] = ( partial_nu[z] - (partial_nu_cos[z] * Cos) ) * integ_coef;
}

///----------------------------------------------
/**
 * Calculates partial Hvv matrix over all angle bins
 * @beam a pointer to the beam array
 * @Cos pre-computed cos(x)
 * @dCos pre-computed dcos(x)
 * @partial_nu returns partial sum over angle/energy bins
 * @partial_nu_cos returns partial sum with cos(a) coeff. over angle/energy bins
 */
void getHvv_partial( const NBeam *REST beam, const double *REST Cos, const double *REST dCos,	///< inputs
						double *REST partial_nu, double *REST partial_nu_cos) throw()			///< output
{
	double* neu_partial;
	double* neu_partial_cos;

	posix_memalign((void**)&neu_partial,  		ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_partial_cos,  	ALIGN_LEN, tetbins*nbm_size*sizeof(double));

	/// Reset the buffer!
	memset(neu_partial, 0., tetbins*nbm_size*sizeof(double));
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		/// loop over flavors
		for (int n = 0; n < nbm_flvs; ++n)
		{
			int neu  = n  << 1;	///< even indices are neu
			int aneu = neu + 1;	///< odds are anti-neu

			Res_t res_neu  AlignAs(ALIGN_LEN) = {};
			Res_t res_aneu AlignAs(ALIGN_LEN) = {};

			///----------------------Calculates Sum_(E)------------------------
			beam[(ang)*nbm_prtcls+neu ].getESum(res_neu );
			beam[(ang)*nbm_prtcls+aneu].getESum(res_aneu);

			///-----------Calculates partial summation over angles-------------
			upd_nu_coef(res_neu, res_aneu, L_E[neu], L_E[aneu], &neu_partial[ang*nbm_size]);
		}	/// End of flavor loop

		for (int z = 0; z < nbm_size; ++z)
		{
			neu_partial[ang*nbm_size+z] *= dCos[ang];

			neu_partial_cos[ang*nbm_size+z] = neu_partial[ang*nbm_size+z]*Cos[ang];
		}

	}	/// End of angle loop

	///Manual reduction!
	for (int a = 0; a < tetbins; ++a)
	{
		for (int z = 0; z < nbm_size; ++z)
		{
			int idx = a*nbm_size+z;

			partial_nu[z]     += neu_partial[idx];
			partial_nu_cos[z] += neu_partial_cos[idx];
		}
	}


	free(neu_partial);
	free(neu_partial_cos);

}

///----------------------------------------------
/**
 * Calculates the hamiltonians 
 * @beam the NBeam array on which the Hamiltonian is calculated
 * @step_num the point at which the Hamiltonian is calculated
 * @hvv the result array
 */
void calc_Hvv(const NBeam *REST beam, const int step_num, double *REST hvv) throw()
{
	Res_t part_Hvv AlignAs(ALIGN_LEN)={}, part_Hvv_cos AlignAs(ALIGN_LEN)={};
#ifdef MMPI
        Res_t mpi_part_Hvv AlignAs(ALIGN_LEN) = {}, mpi_part_Hvv_cos AlignAs(ALIGN_LEN) = {};
        getHvv_partial(beam, &cos_tbin[step_num*tetbins] /*cos(r)*/, &dcos_tbin[step_num*tetbins] /*dcos(r)*/, mpi_part_Hvv, mpi_part_Hvv_cos);
        MPI_Allreduce(mpi_part_Hvv, part_Hvv, nbm_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(mpi_part_Hvv_cos, part_Hvv_cos, nbm_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        getHvv_partial(beam, &cos_tbin[step_num*tetbins] /*cos(r)*/, &dcos_tbin[step_num*tetbins] /*dcos(r)*/, part_Hvv, part_Hvv_cos);
#endif

#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		getHvv(cos_tbin[step_num*tetbins+ang], part_Hvv, part_Hvv_cos, &hvv[ang*nbm_size]);
	}
}

void evolve(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt,ang);
		//double dl = dr / ( (cos_tbin[p1*tetbins+ang] + cos_tbin[p2*tetbins+ang])*.5 );
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].evolveBins(ibeam[idx], prc, dl, &hvv[ang*nbm_size], hmatt);
                        //ibeam[idx].evolveBins(prc, dl, &hvv[ang*nbm_size], hmatt, obeam[idx]);
		}
	}
}

void evolveHvv(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt,ang);
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].evolveBinsHvv(ibeam[idx], prc, dl, &hvv[ang*nbm_size], hmatt);
		}
	}
}

void evolveAvg(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt,ang);
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].evolveBinsAvg(ibeam[idx], prc, dl, &hvv[ang*nbm_size], hmatt, obeamAvg[idx]);
		}
	}
}

void evolveAvg(const int pnt, const double *REST hvv, const double hmatt, NBeam *REST iobeam, NBeam *REST obeamAvg) throw()
{
#ifdef OMP
#       pragma omp parallel for
#endif
        for (int ang = 0; ang < tetbins; ++ang)
        {
                double dl = getDeltaL(pnt,ang);
                for (int prc = 0; prc < nbm_prtcls; ++prc)
                {
                        int idx = ang * nbm_prtcls + prc;
                        iobeam[idx].evolveBinsAvg(prc, dl, &hvv[ang*nbm_size], hmatt, obeamAvg[idx]);
                }
        }
}

void evolveAvgErr(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg, NBeam *REST obeamErr) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt,ang);
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].evolveBinsAvgErr(ibeam[idx], prc, dl, &hvv[ang*nbm_size], hmatt, obeamAvg[idx], obeamErr[idx]);
		}
	}
}

void evolveHvvAvg(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt,ang);
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].evolveBinsHvvAvg(ibeam[idx], prc, dl, &hvv[ang*nbm_size], hmatt, obeamAvg[idx]);
		}
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
	for (int ang = 0; ang < tetbins; ++ang)
	{
		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			obeam[idx].addAvg(ibeam[idx]);
		}
	}
}

} /// Phy namespace

} /// NBGroup namespace
