/*
 * Phy_PNT.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */
#ifdef MMPI
#include <mpi.h>
#endif

#include <stdlib.h>
//#include <malloc.h>
#include <new>
#include <cstring>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fenergy.h"
#include "Fio.h"

#include "Phy.h"

namespace nbgrp {

extern double timeReduceSum;
extern int proc_rank;

namespace phy_pnt {

using namespace fed;
using namespace nbm;


///===================== Variables ========================

/// Starting and ending points for beams, and total #of theta bins
int start_beam, end_beam, abins, pbins;

int beam_len;		///< holds the lenght of the nubeam array

int tetbins; 		///< number of beams on the node based on the computational capability!

/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

double dPhi;		///< will be init in init_sincos() function

/// Arrays for sin(theta), cos(theta),  sin(phi), cos(phi),  and dcos(theta) bins at surface -- have to be computed once
double* sin_tetR;
double* cos_tetR;
double* dcos_tetR;
double* sin_phiR;
double* cos_phiR;

/// cache for delta-r array -- one per theta angle
double* delta_l;

double* L_E;		///< the result of L / <E> for each particle
double jq_aen;
double jq_en ;
double jq_atn;
double jq_tn ;

bool first_beamSet;	///< indicate if we're dealing with beam0 set

/// IO related variables
const int Ndims = 6; 	///< [r,theta,phi,partcl,cmpn,ebin]
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
int getDim()	{ return Ndims; }

///----------------------------------------------
/**
 * Returns naming info for the geometry -- to be used by IO module
 * @param str[] the array of strings indicating the name of each dimention
 */
void getDimInfo(std::string str[])
{
	str[0] = "r";
	str[1] = "theta";
	str[2] = "phi";
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
int& firstDimLen()	{ return Util::Abins(); }

/// Based on the 'i' return neutrino density for a flavor
double Jq(int i)
{
	switch (i) {
		case 0:
			return jq_en;
		case 1:
			return jq_aen;
		case 2:
			return jq_tn;
		case 3:
			return jq_atn;
		default:
			return -1;
	}
}

///----------------------------------------------
/// Init. sin and cos bins at the surface
inline void init_sincos_R() throw()
{
	/// *** Theta bins calculations ***
	double dtet = Util::Pih / tetbins; 		///< create theta bins over [0:pi/2]
	double dtet2 = dtet*.5;
	double tet = dtet2;

#ifdef OMP
//#	pragma omp parallel for
#endif
	for (int t = 0; t < tetbins; ++t, tet+=dtet)
	{
		sin_tetR[t]  = sin(tet);
		cos_tetR[t]  = cos(tet);
		dcos_tetR[t] = cos(tet-dtet2) - cos(tet+dtet2);
	}

	/// *** Phi bins calculations ***
	double dphi = (Util::Pix2) / pbins;	///< create phi bins over [0:2pi]
	dPhi = dphi; ///< init dPhi as well!
	double phi = dphi*.5;

#ifdef OMP
//#	pragma omp parallel for
#endif
	for (int p = 0; p < pbins; ++p, phi+=dphi)
	{
		sin_phiR[p] = sin(phi);
		cos_phiR[p] = cos(phi);
	}

}

///----------------------------------------------
void calcAngleBins(const double r, const int step_num) throw()
{
	return;	/// No usage for this function in this module!
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
	/// In this module there is no usage other arguments after 'dr'
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		delta_l[cur_pnt*tetbins+ang] = dr / cos_tetR[ang];
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
		abins = 720;
        pbins = 100;
		start_beam = 0; 
		end_beam = abins-1;
		printf("Bnechmarking using %dx%d trajectory beams...\n", abins, pbins);
	}
	else	/// start_beam and end_beam are set before!
	{
		abins = Util::Abins();
		pbins = Util::Pbins();
	}

	last_idx = -1; /// beam index is never neg.

	/// start_beam and end_beam are set from NBGroup module
	tetbins = end_beam - start_beam + 1;

	/// Needed by IO module -- should be done before nmr::init->initBeams
	start[0] = start[1] = start[2] = start[3] = start[4] = start[5] = 0;
	count[0] = 1;
	count[1] = tetbins;
	count[2] = pbins;
	count[3] = nbm_prtcls;
	count[4] = nbm_cmpn;
	count[5] = nbm_ebins;

	double mu = 0.249233 * Util::mu();
	jq_en = mu/(Util::SQRT2*Util::GF*2.*Util::Pi);
	jq_aen = jq_en * .8;
	jq_atn = 0.;
	jq_tn = 0.;

	posix_memalign((void**)&sin_tetR,    	ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&cos_tetR,    	ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&dcos_tetR,   	ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&sin_phiR,    	ALIGN_LEN, pbins*sizeof(double));
	posix_memalign((void**)&cos_phiR,    	ALIGN_LEN, pbins*sizeof(double));

	first_beamSet = true;	///< should be done before nmr::init!
	beam_len = tetbins*pbins*nbm_prtcls;
	int pnts_num = nmr::init(beam_len);	///< calls initBeam

	posix_memalign((void**)&delta_l,    	ALIGN_LEN, pnts_num*tetbins*sizeof(double));

	init_sincos_R();
	/// Hvv integral coeff init.	
	hvv_coeff = Util::SQRT2 * Util::GF;
	integ_coef = hvv_coeff * dPhi;//* dE;

	posix_memalign((void**)&L_E,		ALIGN_LEN, nbm_prtcls*sizeof(double));
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

	free(delta_l);

	free(sin_tetR);
	free(cos_tetR);
	free(dcos_tetR);
	free(sin_phiR);
	free(cos_phiR);

	free(L_E);
}

///----------------------------------------------
/**
 * Allocates memory for the hamiltonian array -- to be used by Numeric module
 * @hvv the pointer to the allocated array
 */
void newHvv(double*& hvv)
{
	posix_memalign((void**)&hvv, ALIGN_LEN, pbins*tetbins*nbm_size*sizeof(double));
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
	for (int a = 0; a < tetbins; ++a)
		for (int p = 0; p < pbins; ++p)
			for (int n = 0; n < nbm_prtcls; ++n)
				new (&beam[(a*pbins+p)*nbm_prtcls+n]) NBeam(n);	///< Call constructor for each Beam

/// It is better the rest of this func. exe. in serial!

	/// We just want to add noise to the first set of beam!
	/// As the other beams derived from the first one
	if (first_beamSet)
	{

	for (int a = 0; a < tetbins; ++a)
	{
		//unsigned int seed;
		for (int p = 0; p < pbins; ++p)
		{
			//seed += p;
			for (int n = 0; n < nbm_prtcls; ++n)
			{
				//seed += n;

				int idx = (a*pbins+p)*nbm_prtcls+n;
	
				/// Add some noise!
				for (int e = 0; e <  nbm_ebins; ++e)
				{
					double eps = 1.e-6;//(double( rand() )/RAND_MAX) * (1.e-8) - (.5*1.e-8);
					double eta = 1.e-6;//(double( rand() )/RAND_MAX) * (1.e-8) - (.5*1.e-8);
					
					double one_sq = sqrt(1 - eps*eps - eta*eta);
					
					beam[idx].Ar(e) = (n < 2) ? one_sq : 0./*eps*/;
					beam[idx].Ai(e) = (n < 2) ? 0. : 0.;
					beam[idx].Br(e) = (n < 2) ? eps : 0./*one_sq*/;
					beam[idx].Bi(e) = (n < 2) ? eta : 0.;

				} /// End of energy bins

			} /// End of particle number

		} /// End of pbins

	} /// End tetbins

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
	for (int a = 0; a < tetbins; ++a)
	{
		for (int p = 0; p < pbins; ++p)
		{
			for (int n = 0; n < nbm_prtcls; ++n)
			{
				/// Call destructor!
				beam[(a*pbins+p)*nbm_prtcls+n].~NBeam();
			}
		}
	}
}

///----------------------------------------------
/**
 * Calculates Hvv matrix over all bins
 * @sin_t sin(theta) coeff
 * @cos_t cos(theta) coeff
 * @sin_p sin(phi) coeff
 * @cos_p cos(phi) coeff
 * @part_Hvv_e_p_t partial integral result
 * @part_Hvv_e_cosp_sint partial integral result with cos(phi) * sin(theta)
 * @part_Hvv_e_sinp_sint partial integral result with sin(phi) * sin(theta)
 * @part_Hvv_e_p_cost partial integral result with cos(theta)
 * @ret the whole integral result over angle/energy bins for the current theta
 */
inline void getHvv( const double sin_t, 						///< inputs
					const double cos_t,							///< "
					const double sin_p,							///< "
					const double cos_p,							///< "
					const double *REST part_Hvv_e_p_t, 			///< "
					const double *REST part_Hvv_e_cosp_sint,	///< "
					const double *REST part_Hvv_e_sinp_sint,	///< "
					const double *REST part_Hvv_e_p_cost,		///< "
					double *REST ret ) throw()					///< output
{
	for (int z = 0; z < nbm_size; ++z)
		ret[z]= ( 	part_Hvv_e_p_t[z] -
					( sin_t * cos_p * part_Hvv_e_cosp_sint[z]
					+ sin_t * sin_p * part_Hvv_e_sinp_sint[z]
					+ cos_t * part_Hvv_e_p_cost[z] )
				) * integ_coef;
}

///----------------------------------------------
/**
 * Calculates partial Hvv matrix over all angle bins
 * @beam a pointer to the beam array
 * @step_num indicates the at which step we are doing the calculations
 * @part_Hvv_e_p_t partial integral result
 * @part_Hvv_e_cosp_sint partial integral result with cos(phi) * sin(theta)
 * @part_Hvv_e_sinp_sint partial integral result with sin(phi) * sin(theta)
 * @part_Hvv_e_p_cost partial integral result with cos(theta)
 */
void getHvv_partial(const NBeam *REST beam, const int step_num,	///< inputs
		double *REST part_e_p_t, double *REST part_e_Cosp_Sint, double *REST part_e_Sinp_Sint, double *REST part_e_p_Cost) throw()	///< output
{
	/// Arrays for manual reduction - omp doesn't support reduction over dynamic arrays
	double* neu_partialEngPhiTet;
	double* neu_partialEngCospSint;
	double* neu_partialEngSinpSint;
	double* neu_partialEngPhiCost;

	posix_memalign((void**)&neu_partialEngPhiTet,  		ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_partialEngCospSint, 	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_partialEngSinpSint, 	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_partialEngPhiCost, 		ALIGN_LEN, tetbins*nbm_size*sizeof(double));

	/// Reset the buffer!
	memset(neu_partialEngPhiTet, 0., tetbins*nbm_size*sizeof(double));
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		Res_t neu_partialE_cosP AlignAs(ALIGN_LEN) = {};
		Res_t neu_partialE_sinP AlignAs(ALIGN_LEN) = {};

		for (int phi = 0; phi < pbins; ++phi)
		{
			Res_t neu_partialE AlignAs(ALIGN_LEN) = {};

			/// loop over flavors
			for (int n = 0; n < nbm_flvs; ++n)
			{
				int neu  = n  << 1;	///< even indices are neu
				int aneu = neu + 1;	///< odds are anti-neu

				Res_t res_neu  AlignAs(ALIGN_LEN) = {};
				Res_t res_aneu AlignAs(ALIGN_LEN) = {};

				///----------------------Calculates Sum_(E)------------------------
				beam[(ang*pbins+phi)*nbm_prtcls+neu ].getESum(res_neu );
				beam[(ang*pbins+phi)*nbm_prtcls+aneu].getESum(res_aneu);

				upd_nu_coef(res_neu, res_aneu, Jq(neu), Jq(aneu), neu_partialE);

			}	/// --- End of the flavor loop

			/// Manual reduction!
			for (int z = 0; z < nbm_size; ++z)
			{
				neu_partialEngPhiTet[ang*nbm_size+z] += neu_partialE[z];

				neu_partialE_cosP[z] += neu_partialE[z] * cos_phiR[phi];
				neu_partialE_sinP[z] += neu_partialE[z] * sin_phiR[phi];
			}

		}	/// --- End of the Phi loop

		/// Manual reduction!
		for (int z = 0; z < nbm_size; ++z)
		{
			int idx = ang*nbm_size+z;
			int aid = step_num*tetbins+ang;

			neu_partialEngPhiTet[idx] *= dcos_tetR[aid];

			neu_partialEngCospSint[idx] = neu_partialE_cosP[z] * sin_tetR[aid] * dcos_tetR[aid];
			neu_partialEngSinpSint[idx] = neu_partialE_sinP[z] * sin_tetR[aid] * dcos_tetR[aid];

			neu_partialEngPhiCost[idx] = neu_partialEngPhiTet[idx] * cos_tetR[aid]; ///< there is an implicit dcos(t) here
		}

	}	/// --- End of the angle loop

	/// Manual reduction!
	for (int z = 0; z < nbm_size; ++z)
	{
		for (int ang = 1; ang < tetbins; ++ang)
		{
			int idx = ang*nbm_size+z;

			neu_partialEngPhiTet[z] 	+= neu_partialEngPhiTet[idx];
			neu_partialEngCospSint[z] 	+= neu_partialEngCospSint[idx];
			neu_partialEngSinpSint[z] 	+= neu_partialEngSinpSint[idx];
			neu_partialEngPhiCost[z] 	+= neu_partialEngPhiCost[idx];
		}

		part_e_p_t[z] 		= neu_partialEngPhiTet[z];
		part_e_Cosp_Sint[z]	= neu_partialEngCospSint[z];
		part_e_Sinp_Sint[z]	= neu_partialEngSinpSint[z];
		part_e_p_Cost[z] 	= neu_partialEngPhiCost[z];
	}

	free(neu_partialEngPhiTet);
	free(neu_partialEngCospSint);
	free(neu_partialEngSinpSint);
	free(neu_partialEngPhiCost);
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
	Res_t partHvv_e_p_t AlignAs(ALIGN_LEN)={}, partHvv_e_Cosp_Sint AlignAs(ALIGN_LEN)={}, partHvv_e_Sinp_Sint AlignAs(ALIGN_LEN)={}, partHvv_e_p_Cost AlignAs(ALIGN_LEN)={};
#ifdef MMPI

if (!Util::isBench())
{
        Res_t mpi_partHvv_e_p_t AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Cosp_Sint AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Sinp_Sint AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_p_Cost AlignAs(ALIGN_LEN) = {};
        getHvv_partial(beam, step_num, mpi_partHvv_e_p_t, mpi_partHvv_e_Cosp_Sint, mpi_partHvv_e_Sinp_Sint, mpi_partHvv_e_p_Cost);

        double mpi_sendBuffer[4*nbm_size]  AlignAs(ALIGN_LEN) = {}, mpi_recvBuffer[4*nbm_size]  AlignAs(ALIGN_LEN) = {}; ///< packing four mpi calls into one!
        memcpy(&mpi_sendBuffer[0*nbm_size], mpi_partHvv_e_p_t, 		nbm_size*sizeof(double));
        memcpy(&mpi_sendBuffer[1*nbm_size], mpi_partHvv_e_Cosp_Sint, 	nbm_size*sizeof(double));
        memcpy(&mpi_sendBuffer[2*nbm_size], mpi_partHvv_e_Sinp_Sint, 	nbm_size*sizeof(double));
        memcpy(&mpi_sendBuffer[3*nbm_size], mpi_partHvv_e_p_Cost, 	nbm_size*sizeof(double));

#ifdef DBG
double t1 = MPI_Wtime();
#endif
        MPI_Allreduce(mpi_sendBuffer, mpi_recvBuffer, 4*nbm_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef DBG
double t2 = MPI_Wtime();
timeReduceSum += t2-t1;
#endif

        memcpy(partHvv_e_p_t, 		&mpi_recvBuffer[0*nbm_size], nbm_size*sizeof(double));
        memcpy(partHvv_e_Cosp_Sint, 	&mpi_recvBuffer[1*nbm_size], nbm_size*sizeof(double));
        memcpy(partHvv_e_Sinp_Sint, 	&mpi_recvBuffer[2*nbm_size], nbm_size*sizeof(double));
        memcpy(partHvv_e_p_Cost, 	&mpi_recvBuffer[3*nbm_size], nbm_size*sizeof(double));
}
else
        getHvv_partial(beam, step_num, partHvv_e_p_t, partHvv_e_Cosp_Sint, partHvv_e_Sinp_Sint, partHvv_e_p_Cost);
#else
        getHvv_partial(beam, step_num, partHvv_e_p_t, partHvv_e_Cosp_Sint, partHvv_e_Sinp_Sint, partHvv_e_p_Cost);
#endif


#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = step_num*tetbins+ang;
			getHvv(sin_tetR[aid], cos_tetR[aid], sin_phiR[phi], cos_phiR[phi],
					partHvv_e_p_t, partHvv_e_Cosp_Sint, partHvv_e_Sinp_Sint, partHvv_e_p_Cost, &hvv[(ang*pbins+phi)*nbm_size]);
		}
	}

}


void calc_Hvv(const NBeam *REST beam1, const int step_num1, double *REST hvv1, const NBeam *REST beam2, const int step_num2, double *REST hvv2) throw()
{
	Res_t partHvv_e_p_t1 AlignAs(ALIGN_LEN)={}, partHvv_e_Cosp_Sint1 AlignAs(ALIGN_LEN)={}, partHvv_e_Sinp_Sint1 AlignAs(ALIGN_LEN)={}, partHvv_e_p_Cost1 AlignAs(ALIGN_LEN)={};
	Res_t partHvv_e_p_t2 AlignAs(ALIGN_LEN)={}, partHvv_e_Cosp_Sint2 AlignAs(ALIGN_LEN)={}, partHvv_e_Sinp_Sint2 AlignAs(ALIGN_LEN)={}, partHvv_e_p_Cost2 AlignAs(ALIGN_LEN)={};

#ifdef MMPI
if (!Util::isBench())
{
	/// Compute the partial result of the first Hamiltonian
	Res_t mpi_partHvv_e_p_t1 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Cosp_Sint1 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Sinp_Sint1 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_p_Cost1 AlignAs(ALIGN_LEN) = {};
	getHvv_partial(beam1, step_num1, mpi_partHvv_e_p_t1, mpi_partHvv_e_Cosp_Sint1, mpi_partHvv_e_Sinp_Sint1, mpi_partHvv_e_p_Cost1);

	/// Compute the partial result of the second Hamiltonian	
	Res_t mpi_partHvv_e_p_t2 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Cosp_Sint2 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_Sinp_Sint2 AlignAs(ALIGN_LEN) = {}, mpi_partHvv_e_p_Cost2 AlignAs(ALIGN_LEN) = {};
	getHvv_partial(beam2, step_num2, mpi_partHvv_e_p_t2, mpi_partHvv_e_Cosp_Sint2, mpi_partHvv_e_Sinp_Sint2, mpi_partHvv_e_p_Cost2);

	/// Packing MPI calls into one!
	double mpi_sendBuffer[8*nbm_size]  AlignAs(ALIGN_LEN) = {}, mpi_recvBuffer[8*nbm_size]  AlignAs(ALIGN_LEN) = {};

	/// The first half of the buffer relates to the first Hamiltonian
	memcpy(&mpi_sendBuffer[0*nbm_size],	mpi_partHvv_e_p_t1,		nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[1*nbm_size],	mpi_partHvv_e_Cosp_Sint1,	nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[2*nbm_size],	mpi_partHvv_e_Sinp_Sint1,	nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[3*nbm_size],	mpi_partHvv_e_p_Cost1,		nbm_size*sizeof(double));
	/// Packing the second Hamiltonian into the second half of the buffer
	memcpy(&mpi_sendBuffer[4*nbm_size],	mpi_partHvv_e_p_t2,         nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[5*nbm_size],	mpi_partHvv_e_Cosp_Sint2,   nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[6*nbm_size],	mpi_partHvv_e_Sinp_Sint2,   nbm_size*sizeof(double));
	memcpy(&mpi_sendBuffer[7*nbm_size],	mpi_partHvv_e_p_Cost2,      nbm_size*sizeof(double));

#ifdef DBG
double t1 = MPI_Wtime();
#endif
	MPI_Allreduce(mpi_sendBuffer, mpi_recvBuffer, 8*nbm_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef DBG
double t2 = MPI_Wtime();
timeReduceSum += t2-t1;
#endif

	/// Unpacking the buffer
	memcpy(partHvv_e_p_t1,		&mpi_recvBuffer[0*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_Cosp_Sint1,	&mpi_recvBuffer[1*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_Sinp_Sint1,	&mpi_recvBuffer[2*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_p_Cost1,	&mpi_recvBuffer[3*nbm_size],	nbm_size*sizeof(double));

	memcpy(partHvv_e_p_t2,          &mpi_recvBuffer[4*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_Cosp_Sint2,    &mpi_recvBuffer[5*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_Sinp_Sint2,    &mpi_recvBuffer[6*nbm_size],	nbm_size*sizeof(double));
	memcpy(partHvv_e_p_Cost2,       &mpi_recvBuffer[7*nbm_size],	nbm_size*sizeof(double));

}
else
{
	getHvv_partial(beam1, step_num1, partHvv_e_p_t1, partHvv_e_Cosp_Sint1, partHvv_e_Sinp_Sint1, partHvv_e_p_Cost1);
	getHvv_partial(beam2, step_num2, partHvv_e_p_t2, partHvv_e_Cosp_Sint2, partHvv_e_Sinp_Sint2, partHvv_e_p_Cost2);
}
#else
	getHvv_partial(beam1, step_num1, partHvv_e_p_t1, partHvv_e_Cosp_Sint1, partHvv_e_Sinp_Sint1, partHvv_e_p_Cost1);
	getHvv_partial(beam2, step_num2, partHvv_e_p_t2, partHvv_e_Cosp_Sint2, partHvv_e_Sinp_Sint2, partHvv_e_p_Cost2);
#endif


#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		for (int phi = 0; phi < pbins; ++phi)
		{
				/// First Hvv
				int id1 = step_num1*tetbins+ang;
				getHvv(sin_tetR[id1], cos_tetR[id1], sin_phiR[phi], cos_phiR[phi],
					partHvv_e_p_t1, partHvv_e_Cosp_Sint1, partHvv_e_Sinp_Sint1, partHvv_e_p_Cost1, &hvv1[(ang*pbins+phi)*nbm_size]);
				/// Second Hvv
				int id2 = step_num2*tetbins+ang;
				getHvv(sin_tetR[id2], cos_tetR[id2], sin_phiR[phi], cos_phiR[phi],
					partHvv_e_p_t2, partHvv_e_Cosp_Sint2, partHvv_e_Sinp_Sint2, partHvv_e_p_Cost2, &hvv2[(ang*pbins+phi)*nbm_size]);
		}
	}

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
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = getDeltaL(pnt, ang);

		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = (ang*pbins+phi)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].evolveBins(ibeam[idx], prc, dl, &hvv[aid], hmatt);
			}
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
		double dl = getDeltaL(pnt, ang);

		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = (ang*pbins+phi)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].evolveBinsHvv(ibeam[idx], prc, dl, &hvv[aid], hmatt);
			}
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
		double dl = getDeltaL(pnt, ang);

		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = (ang*pbins+phi)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].evolveBinsAvg(ibeam[idx], prc, dl, &hvv[aid], hmatt, obeamAvg[idx]);
			}
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
                double dl = getDeltaL(pnt, ang);

                for (int phi = 0; phi < pbins; ++phi)
                {
                        int aid = (ang*pbins+phi)*nbm_size;

                        for (int prc = 0; prc < nbm_prtcls; ++prc)
                        {
                                int idx = (ang*pbins + phi) * nbm_prtcls + prc;
                                iobeam[idx].evolveBinsAvg(prc, dl, &hvv[aid], hmatt, obeamAvg[idx]);
                        }
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
		double dl = getDeltaL(pnt, ang);

		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = (ang*pbins+phi)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].evolveBinsAvgErr(ibeam[idx], prc, dl, &hvv[aid], hmatt, obeamAvg[idx], obeamErr[idx]);
			}
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
		double dl = getDeltaL(pnt, ang);

		for (int phi = 0; phi < pbins; ++phi)
		{
			int aid = (ang*pbins+phi)*nbm_size;

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].evolveBinsHvvAvg(ibeam[idx], prc, dl, &hvv[aid], hmatt, obeamAvg[idx]);
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
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		for (int phi = 0; phi < pbins; ++phi)
		{
			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (ang*pbins + phi) * nbm_prtcls + prc;
				obeam[idx].addAvg(ibeam[idx]);
			}
		}
	}
}

}
} /// NBGroup namespace

