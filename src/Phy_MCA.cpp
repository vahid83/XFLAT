/*
 * Phy_MA.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

/// ***************** THIS IS AN OLD AND DEPREACATED MODULE *******************
/// ************ IT IS NOT COMPATIBLE WITH THE REST OF THE CODE ***************
/// ********* THE CODE IS MAINTAINED HERE FOR COMPARISON USAGE ONLY ***********

#include <stdlib.h>
//#include <malloc.h>
#include <new>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fenergy.h"
#include "Fio.h"

#include "Phy.h"

//#define ASYMM

namespace nbgrp {

//extern int abins;

namespace phy_mca {

using namespace fed;
using namespace nbm;

///===================== Variables ========================
/// Starting and ending points for beams, and total #of theta bins
int start_beam, end_beam, abins, pbins;


int tetbins;
int tbins;

/// place holder for partial and Hvv result
Res_t part_hvv;
Res_t part_hvv_cos;
double* Hvv;

double *neu_partial, *aneu_partial, *neu_partial_cos, *aneu_partial_cos;

/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

/// Arrays for cos^2(x) and dcos^2(x) bins at surface -- has to be computed once
double* Cos2_R;
double* dCos2_R;

/// Arrays for cos^2(x) and dcos^2(x) bins at 'r' -- recalculated at each step
double* cos_tbin;
double* dtbin;

double* L_E;	///< the result of L / <E> for each particle

///=================== Implementation =====================

/// Init. cos2 bins
inline void init_cos2R() throw()
{
///warm up!
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int i = 0; i < tetbins; ++i)
	{
		Cos2_R[i] = dCos2_R[i] = 0.;
	}

	/// delta theta -- since we're creating bins on (0-1) uniformly, the distance is fixed
	double dcos_theta = (1. - 0.) / tbins;
	Cos2_R[0] = 1. - dcos_theta*.5;
	dCos2_R[0] = dcos_theta;//1. - Cos2_R[0];
#ifdef ASYMM
	Cos2_R[tbins] = Cos2_R[0];
	dCos2_R[tbins] = dcos_theta;//1. - Cos2_R[tetbins];
#endif
	double cos_theta = Cos2_R[0] - dcos_theta;
	for (int i = 1; i < tbins; ++i, cos_theta -= dcos_theta)
	{
		/// fill up the Cos2 and dCos2 bins for the first half [0:pi/2]
		Cos2_R[i]  = cos_theta;
		dCos2_R[i] = dcos_theta;
#ifdef ASYMM
		Cos2_R[tbins+i] = cos_theta;
		dCos2_R[tbins+i] = dcos_theta;
#endif
	}
}

void calcAngleBins(const double r, const int step_num) throw()
{
	double r1 = Util::Rv() / r;
	double r2 = r1*r1;
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		int idx = step_num*tetbins + ang;
		double c2 = 1. - Cos2_R[ang];

		cos_tbin [idx] = sqrt(1. - r2 * (c2));
		dtbin[idx] = (-.5 * r1 * dCos2_R[ang]) / sqrt(c2 - r2*(c2*c2));
	}
}

void init()
{
	if (Util::isBench()) ///< if the code is running in benchmark mode
        {
                abins = 1800;
                start_beam = 0;
                end_beam = abins-1;
        }
        else
        {
                abins = Util::Abins();
        }

	tetbins = abins;
#ifdef ASYMM
	tbins = tetbins * .5; ///< number of bins in [-pi/2:0] and [0:pi/2]
#else
	tbins = tetbins;
#endif
	posix_memalign((void**)&Hvv, 				ALIGN_LEN, tetbins*nbm_size*sizeof(double));

	posix_memalign((void**)&neu_partial,  		ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_partial,  		ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_partial_cos,  	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_partial_cos,  	ALIGN_LEN, tetbins*nbm_size*sizeof(double));

	posix_memalign((void**)&Cos2_R,    			ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&dCos2_R,   			ALIGN_LEN, tetbins*sizeof(double));

	/// Init. cosine^2(R0)
	init_cos2R();

	/// Hvv integral coeff init.
	hvv_coeff	= (Util::SQRT2 * Util::GF) / (Util::pi * Util::Rv());
#ifndef ASYMM
	hvv_coeff *= 2.;
#endif
	integ_coef = hvv_coeff ;//* dE;

	int pnts_num = nmr::init(tetbins*nbm_prtcls);
	posix_memalign((void**)&cos_tbin,			ALIGN_LEN, tetbins*pnts_num*sizeof(double));
	posix_memalign((void**)&dtbin,    			ALIGN_LEN, tetbins*pnts_num*sizeof(double));

	L_E = new double[nbm_prtcls];
	for (int i = 0; i < nbm_prtcls; ++i)
		L_E[i] = Lv(i)  / Ev(i);

#ifdef FIRST_TOUCH
#	pragma omp parallel for
	for (int i = 0; i < tetbins; ++i)
		for (int z = 0; z < nbm_size; ++z)
			Hvv[i*nbm_size+z] = neu_partial[i*nbm_size+z] = aneu_partial[i*nbm_size+z] =
					neu_partial_cos[i*nbm_size+z] = aneu_partial_cos[i*nbm_size+z] = 0.;
#endif
}

void freemem()
{
	nmr::freemem();

	free(cos_tbin);
	free(dtbin);

	free(Hvv);

	free(neu_partial);
	free(aneu_partial);
	free(neu_partial_cos);
	free(aneu_partial_cos);

	free(Cos2_R);
	free(dCos2_R);

	delete[] L_E;
}

void initBeam(NBeam* beam)
{
	for (int i = 0; i < tbins; ++i)
	{
		for (int j = 0; j < nbm_prtcls; ++j)
		{
			/// Call constructor!
			new (&beam[i*nbm_prtcls+j]) NBeam(j);
		}
	}
#ifdef ASYMM
	for (int i = tbins; i < tetbins; ++i)
	{
		for (int j = 0; j < nbm_prtcls; ++j)
		{
			/// Call constructor!
			new (&beam[i*nbm_prtcls+j]) NBeam(j);
			for (int e = 0; e <  nbm_ebins; ++e)
			{
				double eps = (double(rand())/RAND_MAX) * (1.e-6) - (.5*1.e-6);
				double one_sq = sqrt(1 - eps*eps);
				beam[i*nbm_prtcls+j].Ar(e) = (j < 2) ? one_sq : eps;
				beam[i*nbm_prtcls+j].Ai(e) = 0.;
				beam[i*nbm_prtcls+j].Br(e) = (j < 2) ? eps : one_sq;
				beam[i*nbm_prtcls+j].Bi(e) = 0.;
			}
		}
	}
#endif
}

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

/**
 * Calculates Hvv matrix over all bins
 * @Cos cos(x) coeff
 * @partial_nu partial integral result
 * @partial_nu_cos partial integral result with cos(x) term
 * @ret the whole integral result over angle/energy bins for the current theta
 */
void getHvv( const double Cos, const double *REST partial_nu, const double *REST partial_nu_cos,	///< inputs
				double *REST ret) throw()															///< output
{
	for (int i = 0; i < nbm_size; ++i)
		ret[i] = ( partial_nu[i] - (partial_nu_cos[i] * Cos) ) * integ_coef;
}

/**
 * Calculates partial Hvv matrix over all angle bins
 * @beam a pointer to the beam array
 * @Cos pre-computed cos(x)
 * @dCos pre-computed dcos(x)
 * @partial_nu returns partial sum over angle/energy bins
 * @partial_nu_cos returns partial sum with cos(a) coeff. over angle/energy bins
 */
void getHvv_partial( const NBeam *REST beam, const double *REST Cos, const double *REST dtet,	///< inputs
						double *REST partial_nu, double *REST partial_nu_cos) throw()			///< output
{
	for (int n = 0; n < nbm_flvs; ++n)
	{
		int neu  = n  << 1;	///< even indices are neu
		int aneu = neu + 1;	///< odds are anti-neu

#ifdef OMP
#		pragma omp parallel for
#endif
		for (int ang = 0; ang < tetbins; ++ang)
		{
			Res_t  res_neu = {};
			Res_t res_aneu = {};

			///----------------------Calculates Sum_(E)------------------------
			beam[ang*nbm_prtcls+neu].getESum(res_neu);
			beam[ang*nbm_prtcls+aneu].getESum(res_aneu);
			///-----------Calculates partial summation over angles-------------
			for (int z = 0; z < nbm_size; ++z)
			{
				double cfn  = res_neu [z] * dtet[ang];
				double cfan = res_aneu[z] * dtet[ang];

				int idx = ang*nbm_size+z;

				neu_partial [idx] = cfn;
				aneu_partial[idx] = cfan;

				neu_partial_cos [idx] = cfn  * Cos[ang];
				aneu_partial_cos[idx] = cfan * Cos[ang];
			}

		}	/// End of angle loop

		/// Manual reduction!!!
		Res_t n_part = {}, an_part = {}, n_part_cos = {}, an_part_cos = {};
		for (int a = 0; a < tetbins; ++a)
			for (int z = 0; z < nbm_size; ++z)
			{
				int idx = a*nbm_size+z;

				n_part[z]      += neu_partial[idx];
				an_part[z]     += aneu_partial[idx];

				n_part_cos[z]  += neu_partial_cos[idx];
				an_part_cos[z] += aneu_partial_cos[idx];
			}

		double neu_coeff  = L_E[neu];//Lv(neu)  / Ev(neu);
		double aneu_coeff = L_E[aneu];//Lv(aneu) / Ev(aneu);

		upd_nu_coef(n_part, an_part, neu_coeff, aneu_coeff, partial_nu);
		upd_nu_coef(n_part_cos, an_part_cos, neu_coeff, aneu_coeff, partial_nu_cos);
	}	/// End of flavor loop
}

void calc_Hvv(const NBeam* beam, const int step_num) throw()
{
	/// reset the tmp
	memset(part_hvv, 0., nbm_size*sizeof(double));
	memset(part_hvv_cos, 0., nbm_size*sizeof(double));

	getHvv_partial(beam, &cos_tbin[step_num*tetbins] /*cos(r)*/, &dtbin[step_num*tetbins] /*dtheta(r)*/, part_hvv, part_hvv_cos);
}

void evolve(const NBeam *REST ibeam, const double dr, const double hmatt, const int p1, const int p2, bool calc, const int q, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int ang = 0; ang < tetbins; ++ang)
	{
		double dl = dr / ( (cos_tbin[p1*tetbins+ang] + cos_tbin[p2*tetbins+ang])*.5 );

		if (calc) /// if already calculated, no need to recalculate
			getHvv(cos_tbin[q*tetbins+ang], part_hvv, part_hvv_cos, &Hvv[ang*nbm_size]);

		for (int prc = 0; prc < nbm_prtcls; ++prc)
		{
			int idx = ang * nbm_prtcls + prc;
			ibeam[idx].evolveBins(prc, dl, &Hvv[ang*nbm_size], hmatt, obeam[idx]);
		}
	}
}

void avgBeam(const NBeam *REST ibeam, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for// collapse(2)
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

