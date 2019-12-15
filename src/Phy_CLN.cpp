/*
 * Phy_CLN.cpp
 *
 *  Created on: Aug 30, 2014
 *      Author: vahid
 */
#ifdef MMPI
#include <mpi.h>
#endif

#include <cstdio>
#include <stdlib.h>
#include <new>
#include <cstring>
#include <algorithm>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fenergy.h"
#include "Fio.h"

#include "Phy.h"

#include <iostream>

namespace nbgrp {

extern int proc_rank, proc_size;
extern double timeReduceSum;

namespace phy_cln {

using namespace fed;
using namespace nbm;

///===================== Variables ========================

/// Starting and ending points for beams, and total #of points on the surface and theta bins per point
int start_beam, end_beam, source_points, theta_bins;

int node_points;	///< number of source points on the current compute node, based on its computational capability!

int beam_len;		///< holds the lenght of the nubeam array

double* R_r;		///< Holds the precalculated results of R/r

/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

/// cache for delta-r array
double* delta_l;

/// Arrays for alpha, theta, sin(x), cos(x) bins at the surface -- need to be computed only once
double* Alpha_R;
double* Theta_R;
double* Sin_R;
double* Cos_R;
double dAlpha, dTet;	///< delta alpha and theta in radian
double* L_E;		///< the result of L / <E> for each particle

bool first_beamSet;     ///< indicate if we're dealing with beam0 set

/// IO related variables
const int Ndims = 6; 	///< [r,surface_point,theta,partcl,cmpn,ebin]
size_t start[Ndims], count[Ndims];
int last_idx;		///< to prevent recalculating all of the indecis when only component index is changed!

///=================== Implementation =====================

/**
 * Returns the length of the NBeam object array for the current node
 * @retval beam array len. in the module
 */
int beamLen()   { return beam_len; }

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
        str[1] = "source_points";
        str[2] = "theta";
        str[3] = "prtcl";
        str[4] = "comp";
        str[5] = "ebin";
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
size_t* countDim()	{ return count; }

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
int& firstDimLen()       { return Util::SPoints(); }

///----------------------------------------------
/**
 * Init. angle bins at the surface
 */
inline void init_bins_R() throw()
{
/*
/// warm up!
#ifdef OMP
#       pragma omp parallel for
#endif
	for (int i = 0; i < surface_points; ++i)
		Alpha_R[i] = 0.;

#ifdef OMP
#	pragma omp parallel for
#endif
	for (int i = 0; i < theta_bins; ++i)
		Theta_R[i] = Cos_R[i] = Sin_R[i] = 0.;
*/
	/// --- Init alpha bins ---
	dAlpha = Util::Pix2 / source_points;
	double a = dAlpha*.5;
	for (int i = 0; i < source_points; ++i, a+=dAlpha)
	{
		Alpha_R[i] = a;
	}
	/// --- Init theta, sin and cos bins section ---
	dTet = Util::Pi / theta_bins;
	double t = -Util::Pih + dTet*.5; 
	for (int i = 0; i < theta_bins; ++i, t+=dTet)
	{
		Theta_R[i] = t;
		Cos_R[i] = cos(t);
		Sin_R[i] = sin(t);
	}

}

///----------------------------------------------
/**
 * In this module it only calculates and caches the R/r
 * @r the current radius at which the bins are computed
 * @step_num indicates the current step number (first, mid, or last poiint)
 */
void calcAngleBins(const double r, const int step_num) throw()
{
	/// It's not supernova! we only need to have R/r
	R_r[step_num] = Util::Rv() / r;
/*
	double ra = Util::Rv() / r;
	double r2 = ra * ra;
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		int idx = step_num*tetbins + ang;
		int tet = start_beam+ang;

		cos_tetbins [idx] = sqrt(1. - r2 *  (1. - Cos2_R[tet]));
		sin_tetbins [idx] = ra * Sin_R[tet];
		dcos_tetbins[idx] = .5 * (1./cos_tetbins[idx]) * r2 * dCos2_R[tet];
	}
*/
}

///----------------------------------------------
/// Returns the dl array's pointer for the input point
inline double getDeltaL(const int pnt) throw()
{
	return delta_l[pnt];
//	return delta_l[pnt*tbins+idx];
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
	/// It's not supernova! so dl = dr
	delta_l[cur_pnt] = dr;
/*
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		delta_l[cur_pnt*tbins+ang] = dr / ( (cos_tetbins[s_pnt*tbins+ang] + cos_tetbins[e_pnt*tbins+ang])*.5 );
	}
*/
}

///----------------------------------------------
/**
 * Initializes ithe module and allocates memory for arrays
 */
void init()
{
	if (Util::isBench()) ///< if the code is running in benchmark mode
	{
		source_points = 360;
		start_beam = 0; 
		end_beam = source_points-1;
        	theta_bins = 300;
		printf("Bnechmarking using %dx%d trajectory beams...\n", source_points, theta_bins);
	}
	else	/// start_beam and end_beam are set before!
	{
		source_points = Util::SPoints();
	        theta_bins = Util::Abins();
	}

	last_idx = -1; /// beam index is never neg.

	/// start_beam and end_beam are set from NBGroup module
	node_points = end_beam - start_beam + 1;

	/// Needed by IO module -- should be done before nmr::init->initBeams
	start[0] = start[1] = start[2] = start[3] = start[4] = start[5] = 0;
        count[0] = 1;
        count[1] = node_points;
        count[2] = theta_bins;
        count[3] = nbm_prtcls;
        count[4] = nbm_cmpn;
        count[5] = nbm_ebins;

	posix_memalign((void**)&Alpha_R,	ALIGN_LEN, source_points*sizeof(double));
	posix_memalign((void**)&Theta_R,	ALIGN_LEN, theta_bins*sizeof(double));
	posix_memalign((void**)&Sin_R, 		ALIGN_LEN, theta_bins*sizeof(double));
	posix_memalign((void**)&Cos_R, 		ALIGN_LEN, theta_bins*sizeof(double));

	/// Init. alpha, theta, sin and cos bins and calc. dTet 
	init_bins_R();

	first_beamSet = true;   ///< should be done before nmr::init!
	beam_len = node_points*theta_bins*nbm_prtcls;
	int pnts_num = nmr::init(beam_len);		///< calls initBeam
	posix_memalign((void**)&R_r,		ALIGN_LEN, pnts_num*sizeof(double));
	
	///< FixMe! the len. should be proportional to number of angels, but for now we assume they are all the same!!
	posix_memalign((void**)&delta_l,	ALIGN_LEN, pnts_num*sizeof(double));

	/// Hvv integral coeff init.	
	hvv_coeff  = 1E-50; //(Util::SQRT2 * Util::GF) / (4. * Util::pi * Util::pi * Util::Rv() * Util::Rv());
        integ_coef = hvv_coeff  * dTet;//* dE;

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

	free(R_r);
	free(delta_l);

	free(Alpha_R);
	free(Theta_R);
	free(Sin_R);
	free(Cos_R);

	free(L_E);
}

///----------------------------------------------
/**
 * Allocates memory for the hamiltonian array -- to be used by Numeric module
 * @hvv the pointer to the allocated array
 */
void newHvv(double*& hvv)
{
	posix_memalign((void**)&hvv, ALIGN_LEN, node_points*theta_bins*nbm_size*sizeof(double));
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
///---------------------------------------------
#define COS_PERT	///< If defined the perturbation distribution is based on cosine
/**
 * Initializes an array of NBeam objects by calling the constructor
 * @beam array of neutrino beams to be init.
 */
void initBeam(NBeam* beam)
{
#ifdef OMP
#       pragma omp parallel for
#endif
	for (int a = 0; a < node_points; ++a)
		for (int t = 0; t < theta_bins; ++t)
			for (int n = 0; n < nbm_prtcls; ++n)
				new (&beam[(a*theta_bins+t)*nbm_prtcls+n]) NBeam(n);	///< Call constructor for every beam

/// It is better the rest of the func. only exec. in serial

	/// Add perturbation to the first set of beams as the other sets derive from it
	if (first_beamSet)
	{
	double al = start_beam*dAlpha + dAlpha*.5;
	for (int a = 0; a < node_points; ++a, al+=dAlpha)
	{
		/// To make sure that all surface points have the same perturbation
		srand(0);

		for (int t = 0; t < theta_bins; ++t)
		{
			for (int n = 0; n < nbm_prtcls; ++n)
			{
                                for (int e = 0; e <  nbm_ebins; ++e)
                                {
					int idx = (a*theta_bins+t)*nbm_prtcls+n;

                                        double eps = (double( rand() )/RAND_MAX) * (1.e-8) - (.5*1.e-8);
                                        double eta = (double( rand() )/RAND_MAX) * (1.e-8) - (.5*1.e-8);
#ifdef COS_PERT
					eps *= cos(al);
                                        eta *= cos(al);
#endif
                                        double one_sq = sqrt(1 - eps*eps - eta*eta);
                                        beam[idx].Ar(e) = (n < 2) ? one_sq : eps;
                                        beam[idx].Ai(e) = 0;
                                        beam[idx].Br(e) = (n < 2) ? eps : one_sq;
                                        beam[idx].Bi(e) = eta;
                                
                                } /// End energy bins

			} /// End number of particles

		} /// End theta_bins

	} /// End source_points

	} /// End if

	if (Util::inVals() == NULL || Util::isBench() || !first_beamSet) return;
        fio::fillInitData(beam);	///< if a data file is provided for resumption
        first_beamSet = false;
}

///----------------------------------------------
/**
 * Call destructors of NBeam objects
 * @beam the pointer to the NBeam array
 */
void freeBeam(NBeam* beam)
{
	for (int a = 0; a < node_points; ++a)
        {
                for (int t = 0; t < theta_bins; ++t)
                {
                        for (int n = 0; n < nbm_prtcls; ++n)
                        {
				/// Call destructor!
				beam[(a*theta_bins+t)*nbm_prtcls+n].~NBeam();
			}
		}
	}
}

///----------------------------------------------
/**
 * Calculates the hamiltonians 
 * @beam the NBeam array on which the Hamiltonian is calculated
 * @step_num the point at which the Hamiltonian is calculated
 * @hvv the result array
 */
void calc_Hvv(const NBeam *REST beam, const int step_num, double *REST hvv) throw()	///< FixMe! should the hvv be passed by restrict?
{
// **** Extract the summation over energy bins from all NBeam classes *********

	double* beam_ESum;
	//double* node_ESum;
	size_t tot_size = source_points*theta_bins*nbm_prtcls*sizeof(Res_t);
	posix_memalign((void**)&beam_ESum, ALIGN_LEN, tot_size);
	//size_t node_size = node_points*theta_bins*nbm_prtcls*sizeof(Res_t);
	//posix_memalign((void**)&node_ESum, ALIGN_LEN, node_size);
	
#ifdef OMP
//#       pragma omp parallel for
#endif
	for (int srcPnt = 0; srcPnt < node_points; ++srcPnt)
		for (int tet = 0; tet < theta_bins; ++tet)
			for (int p = 0; p < nbm_prtcls; ++p)
			{
				int ind = ( srcPnt*theta_bins + tet )*nbm_prtcls + p;		///< Index into the current compute node's beams
				int idx = ( ( (srcPnt+start_beam)*theta_bins + tet )*nbm_prtcls + p )*nbm_size;	///< each compute node, write their result on their own section
				beam[ind].getESum(&beam_ESum[idx]);	///< The class' method directly copies the return to the buffer
			}

#ifdef MMPI
if (!Util::isBench())
{
	/// Gather starting points for each node -- to be used to in the beam_ESum buffer
	int* proc_starts = new int[proc_size]();	
#ifdef DBG
double t1 = MPI_Wtime();
#endif
	MPI_Allgather(&start_beam, 1, MPI_INT, proc_starts, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef DBG
double t2 = MPI_Wtime();
timeReduceSum += t2-t1;
#endif
	for (int i = 0; i < proc_size; ++i)
		proc_starts[i] *= theta_bins*nbm_prtcls*nbm_size;

	int* proc_chunks = new int[proc_size]();
	/// Calculated each node's data chunk size
	for (int i = 0; i < proc_size; ++i)
		/// There last node is a special case!
		proc_chunks[i] = (i+1 < proc_size) ? proc_starts[i+1] - proc_starts[i] : source_points*theta_bins*nbm_prtcls*nbm_size - proc_starts[i];

#ifdef DBG
t1 = MPI_Wtime();
#endif
	/// Compute nodes will exchange their 'node_ESum' data buffers with eachother
	MPI_Allgatherv(MPI_IN_PLACE, /*ignored*/0, /*ignored*/MPI_DATATYPE_NULL, beam_ESum, proc_chunks, proc_starts, MPI_DOUBLE, MPI_COMM_WORLD);
#ifdef DBG
t2 = MPI_Wtime();
timeReduceSum += t2-t1;
#endif

	delete[] proc_chunks;
	delete[] proc_starts;
}		
#endif /// End of MMPI

// ****************************************************************************

	const double* alpha_c = Alpha_R; 
#ifdef OMP
#       pragma omp parallel for //collapse(2)
#endif	
	for (int srcPnt = 0; srcPnt < node_points; ++srcPnt)
	{
		for (int tet = 0; tet < theta_bins; ++tet)
		{
			Res_t partial_H AlignAs(ALIGN_LEN) = {};

			/// Loop over theta'
			for (int tp = 0; tp < theta_bins; ++tp)
			{

				double alpha_prime = alpha_c[srcPnt+start_beam] + Theta_R[tp] - Theta_R[tet] - asin(R_r[step_num]*Sin_R[tp]) + asin(R_r[step_num]*Sin_R[tet]);
				alpha_prime += alpha_prime < 0. ? Util::Pix2 : (alpha_prime > Util::Pix2 ? -Util::Pix2 : 0.); ///< add 2xPi to neg. rad. and sub 2xPi from those that are more than 2xPi

				int ap = std::upper_bound(alpha_c, alpha_c+source_points, alpha_prime) - alpha_c;
				int y = ap % source_points;
				int x = y == 0 ? source_points-1 : y-1;

				double dRad = fabs(alpha_c[x] - alpha_prime);
				dRad = dRad > dAlpha ? Util::Pix2-dRad : dRad; ///< when one angle is above 0-line and one is below!
				double intrpltn_coeff = dRad / dAlpha;

				Res_t neu_ESum0 AlignAs(ALIGN_LEN) = {};
				Res_t neu_ESum1 AlignAs(ALIGN_LEN) = {};

				/// Loop over flavors
				for (int n = 0; n < nbm_flvs; ++n)
			        {
        			        int neu  = n  << 1;     ///< even indices are neu
                			int aneu = neu + 1;     ///< odds are anti-neu

        			        double *res_neu0, *res_aneu0, *res_neu1, *res_aneu1;

			                ///----------------------Calculates Sum_(E)------------------------
					res_neu0  = &beam_ESum[((x * theta_bins + tp) * nbm_prtcls + neu ) * nbm_size];
					res_aneu0 = &beam_ESum[((x * theta_bins + tp) * nbm_prtcls + aneu) * nbm_size];

                                        res_neu1  = &beam_ESum[((y * theta_bins + tp) * nbm_prtcls + neu ) * nbm_size];
					res_aneu1 = &beam_ESum[((y * theta_bins + tp) * nbm_prtcls + aneu) * nbm_size];

					upd_nu_coef(res_neu0, res_aneu0, (1.-intrpltn_coeff) * L_E[neu], (1.-intrpltn_coeff) * L_E[aneu], neu_ESum0);
					upd_nu_coef(res_neu1, res_aneu1, intrpltn_coeff * L_E[neu], intrpltn_coeff * L_E[aneu], neu_ESum1);
				}
				/// Add other components!
				for (int z = 0; z < nbm_size; ++z)
				{
					double sin2 = Sin_R[tet] - Sin_R[tp];
					partial_H[z] += (neu_ESum0[z] + neu_ESum1[z]) * Cos_R[tp] * sin2*sin2;
				}

			}
			/// Add (R/r)^3 and save the Hamiltonian!
			for (int z = 0; z < nbm_size; ++z)
				hvv[(srcPnt*theta_bins+tet)*nbm_size+z] = partial_H[z] * R_r[step_num]*R_r[step_num]*R_r[step_num] * integ_coef;
			
		} /// End of theta beams

        } /// End of source points                     

	free(beam_ESum);
}

///----------------------------------------------
/**
 * Calls the evolution function of NBeam class
 * @ibeam the first beam
 * @pnt the point at which the evolution is done
 * @hvv the hamiltonian
 * @hmatt the matter
 * @obeam the result beam
 */
void evolve(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
	double dl = getDeltaL(pnt);
#ifdef OMP
#	pragma omp parallel for //collapse(2)
#endif
	for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
			int h_idx = (a*theta_bins + tet)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].evolveBins(ibeam[idx], prc, dl, &hvv[h_idx], hmatt);
			}
		}
	}
}

void evolveHvv(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam) throw()
{
	double dl = getDeltaL(pnt);
#ifdef OMP
#	pragma omp parallel for //collapse(2)
#endif
	for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        int h_idx = (a*theta_bins + tet)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].evolveBinsHvv(ibeam[idx], prc, dl, &hvv[h_idx], hmatt);
			}
		}
	}
}

void evolveAvg(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg) throw()
{
	double dl = getDeltaL(pnt);
#ifdef OMP
#       pragma omp parallel for //collapse(2)
#endif
        for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        int h_idx = (a*theta_bins + tet)*nbm_size;

                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].evolveBinsAvg(ibeam[idx], prc, dl, &hvv[h_idx], hmatt, obeamAvg[idx]);
			}
		}
	}
}

void evolveAvg(const int pnt, const double *REST hvv, const double hmatt, NBeam *REST iobeam, NBeam *REST obeamAvg) throw()
{
        double dl = getDeltaL(pnt);
#ifdef OMP
#       pragma omp parallel for //collapse(2)
#endif
	for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        int h_idx = (a*theta_bins + tet)*nbm_size;

                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
                                iobeam[idx].evolveBinsAvg(prc, dl, &hvv[h_idx], hmatt, obeamAvg[idx]);
                        }
                }
        }
}


void evolveAvgErr(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg, NBeam *REST obeamErr) throw()
{
	double dl = getDeltaL(pnt);
#ifdef OMP
#       pragma omp parallel for //collapse(2)
#endif
        for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        int h_idx = (a*theta_bins + tet)*nbm_size;

                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].evolveBinsAvgErr(ibeam[idx], prc, dl, &hvv[h_idx], hmatt, obeamAvg[idx], obeamErr[idx]);
			}
		}
	}
}

void evolveHvvAvg(const NBeam *REST ibeam, const int pnt, const double *REST hvv, const double hmatt, NBeam *REST obeam, NBeam *REST obeamAvg) throw()
{
	double dl = getDeltaL(pnt);
#ifdef OMP
#       pragma omp parallel for //collapse(2)
#endif
        for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        int h_idx = (a*theta_bins + tet)*nbm_size;

                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].evolveBinsHvvAvg(ibeam[idx], prc, dl, &hvv[h_idx], hmatt, obeamAvg[idx]);
			}
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
#       pragma omp parallel for //collapse(2)
#endif
        for (int a = 0; a < node_points; ++a)
        {
                for (int tet = 0; tet < theta_bins; ++tet)
                {
                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (a*theta_bins + tet)*nbm_prtcls + prc;
				obeam[idx].addAvg(ibeam[idx]);
			}
		}
	}
}

}
} /// NBGroup namespace

