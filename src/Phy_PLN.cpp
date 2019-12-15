/*
 * Phy_PLN.cpp
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

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fenergy.h"
#include "Fio.h"

#include "Phy.h"

namespace nbgrp {

extern int abins, pbins;

namespace phy_pln {

using namespace fed;
using namespace nbm;

///===================== Variables ========================
/// Starting and ending points for beams, and total #of theta bins
int start_beam, end_beam, abins, pbins;


int tetbins, phibins; ///< to be accessible from Fio

/// global holder for partial Hvv result
Res_t partial_nu;
Res_t partial_nu_costet_phi;
Res_t partial_nu_sintet_sinphi;
Res_t partial_nu_sintet_cosphi;

/// global holder for Hvv
double* Hvv_costet_phi_cos;
double* Hvv_sintet_sinphi_sin;
double* Hvv_sintet_cosphi_sin;
double* Hvv;

double* sin_tet;
double* cos_tet;
double* sin_phi;
double* cos_phi;
double* dcos_tet;

/// Arrays for manual reduction - omp doesn't support reduction over dynamic arrays
double* neu_partial_E;
double* aneu_partial_E;
double* neu_costet_phi;
double* aneu_costet_phi;
double* neu_sintet_sinphi;
double* aneu_sintet_sinphi;
double* neu_sintet_cosphi;
double* aneu_sintet_cosphi;

double dPhi;	///< will be init in sincos function

/// some coefficients which are calculated once
double hvv_coeff, integ_coef;

double jq_aen;
double jq_en ;
double jq_atn;
double jq_tn ;

double* L_E;	///< the result of L / <E> for each particle

///=================== Implementation =====================

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
			return jq_en;
	}
}

/// Init. sin and cos bins
void init_sincos() throw()
{
	double dtet = Util::pi2 / tetbins; 		///< create theta bins over [0:pi/2]
	double dtet2 = dtet*.5;
	double tet = dtet2;
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int t = 0; t < tetbins; ++t)
	{
		sin_tet[t] = sin(tet);
		cos_tet[t] = cos(tet);
		dcos_tet[t] = cos(tet-dtet2) - cos(tet+dtet2);
		tet+=dtet;
	}

	double dphi = (Util::pi*2.) / phibins;	///< create phi bins over [0:2pi]
	dPhi = dphi; ///< init dPhi as well!
	double phi = dphi*.5;
#ifdef OMP
#	pragma omp parallel for
#endif
	for (int p = 0; p < phibins; ++p)
	{
		sin_phi[p] = sin(phi);
		cos_phi[p] = cos(phi);
		phi+=dphi;
	}
}

void calcAngleBins(const double r, const int step_num) throw()
{
	return;
}

void init()
{
	if (Util::isBench()) ///< if the code is running in benchmark mode
        {
                abins = 720;
                pbins = 100;
                start_beam = 0;
                end_beam = abins-1;
        }
        else	/// start_beam and end_beam are set before!
        {
                abins = Util::Abins();
                pbins = Util::Pbins();
        }

	tetbins = abins;
	phibins = pbins;

	double mu = 0.249233 * Util::mu();
	jq_aen = mu/(Util::SQRT2*2.*Util::pi*Util::GF);
	jq_en = jq_aen * 1.5;
	jq_atn = jq_aen * .5;
	jq_tn = jq_aen * .5;

	posix_memalign((void**)&Hvv_costet_phi_cos, 	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&Hvv_sintet_sinphi_sin,	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&Hvv_sintet_cosphi_sin,	ALIGN_LEN, tetbins*nbm_size*sizeof(double));
	posix_memalign((void**)&Hvv,					ALIGN_LEN, phibins*tetbins*nbm_size*sizeof(double));

	posix_memalign((void**)&sin_phi,   				ALIGN_LEN, phibins*sizeof(double));
	posix_memalign((void**)&cos_phi,   				ALIGN_LEN, phibins*sizeof(double));
	posix_memalign((void**)&sin_tet,   				ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&cos_tet,   				ALIGN_LEN, tetbins*sizeof(double));
	posix_memalign((void**)&dcos_tet,  				ALIGN_LEN, tetbins*sizeof(double));

	posix_memalign((void**)&neu_partial_E,    		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_partial_E,    		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_costet_phi,    		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_costet_phi,   		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_sintet_sinphi,		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_sintet_sinphi,		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&neu_sintet_cosphi, 		ALIGN_LEN, phibins*nbm_size*sizeof(double));
	posix_memalign((void**)&aneu_sintet_cosphi,		ALIGN_LEN, phibins*nbm_size*sizeof(double));

	/// Init. sin and cos bins and also init dPhi
	init_sincos();

	/// Hvv integral coeff init.
	hvv_coeff = Util::SQRT2 * Util::GF;
	integ_coef = hvv_coeff * dPhi;//* dE;

	nmr::init(phibins*tetbins*nbm_prtcls);

	L_E = new double[nbm_prtcls];
	for (int i = 0; i < nbm_prtcls; ++i)
		L_E[i] = Lv(i)  / Ev(i);
}

void freemem()
{
	nmr::freemem();

	free(Hvv_costet_phi_cos);
	free(Hvv_sintet_sinphi_sin);
	free(Hvv_sintet_cosphi_sin);
	free(Hvv);

	free(sin_tet);
	free(cos_tet);
    free(sin_phi);
	free(cos_phi);
	free(dcos_tet);

	free(neu_partial_E);
	free(aneu_partial_E);
	free(neu_costet_phi);
	free(aneu_costet_phi);
	free(neu_sintet_sinphi);
	free(aneu_sintet_sinphi);
	free(neu_sintet_cosphi);
	free(aneu_sintet_cosphi);

	delete[] L_E;
}

void initBeam(NBeam* beam)
{
#if defined(OMP) && !defined(__MIC__)
#       pragma omp parallel for
#endif
	for (int i = 0; i < phibins; ++i)
		for (int j = 0; j < tetbins; ++j)
			for (int k = 0; k < nbm_prtcls; ++k)
				/// Call constructor!
				new (&beam[(i*tetbins+j)*nbm_prtcls+k]) NBeam(k);
}

void freeBeam(NBeam* beam)
{
	for (int i = 0; i < phibins; ++i)
		for (int j = 0; j < tetbins; ++j)
			for (int k = 0; k < nbm_prtcls; ++k)
				/// Call destructor!
				beam[(i*tetbins+j)*nbm_prtcls+k].~NBeam();
}

/**
 * Calculates Hvv matrix over all bins
 * @Cos cos(x) coeff
 * @partial_nu partial integral result
 * @partial_nu_cos partial integral result with cos(x) term
 * @ret the whole integral result over angle/energy bins for the current theta
 */
void getHvv( const int phi,
			const double *REST part_nu_E, const double *REST nu_costet_phi_cos, const double *REST nu_sintet_sinphi_sin, const double *REST nu_sintet_cosphi_sin,	///< inputs
			double *REST hvv) throw()	///< output
{
	for (int i = 0; i < nbm_size; ++i)
		hvv[i] = ( part_nu_E[i] - ( (nu_sintet_cosphi_sin[i] * cos_phi[phi]) + (nu_sintet_sinphi_sin[i] * sin_phi[phi]) + nu_costet_phi_cos[i] ) ) * integ_coef;
}

/**
 * Calculates partial cos(tet) and sin(tet) multiplication for Hvv
 * @theta current angle
 * @part_costet partial result from previous step
 * @part_sintet_sinphi partial result from previous step
 * @part_sintet_cosphi partial result from previous step
 * @ret* the partial result after coefficient multiplication
 */
void getHvv_partial_tet( const int theta, const double *REST part_costet_phi, const double *REST part_sintet_sinphi, const double *REST part_sintet_cosphi,	///< inputs
						double *REST ret_costet_phi_cos, double *REST ret_sintet_sinphi_sin, double *REST ret_sintet_cosphi_sin) throw()	///< output
{
	for (int i = 0; i < nbm_size; ++i)
	{
		ret_costet_phi_cos[i]    = part_costet_phi[i]    * cos_tet[theta];
		ret_sintet_sinphi_sin[i] = part_sintet_sinphi[i] * sin_tet[theta];
		ret_sintet_cosphi_sin[i] = part_sintet_cosphi[i] * sin_tet[theta];
	}
}

/**
 * Calculates partial Hvv matrix over all angle bins
 * @beam a pointer to the beam array
 * @Cos pre-computed cos(x)
 * @dCos pre-computed dcos(x)
 * @partial_nu returns partial sum over angle/energy bins
 * @partial_nu_cos returns partial sum with cos(a) coeff. over angle/energy bins
 */
void getHvv_partial( const NBeam *REST beam, ///< inputs
					double *REST part_nu, double *REST part_nu_costet_phi, double *REST part_nu_sintet_sinphi, double *REST part_nu_sintet_cosphi) throw()			///< output
{
	/// flavor loop
	for (int n = 0; n < nbm_flvs; ++n)
	{
		int neu  = n  << 1;	///< even indices are neu
		int aneu = neu + 1;	///< odds are anti-neu

#ifdef OMP
#		pragma omp parallel for
#endif
		for (int i = 0; i < phibins; ++i)
		{
			Res_t  nu_sintet = {};
			Res_t anu_sintet = {};

			for (int j = 0; j < tetbins; ++j)
			{

				Res_t nu_E  = {};
				Res_t anu_E = {};
				///----------------------Calculates Sum_(E)------------------------
				beam[(i*tetbins+j)*nbm_prtcls+neu].getESum(nu_E);
				beam[(i*tetbins+j)*nbm_prtcls+aneu].getESum(anu_E);
				///----------------------------------------------------------------

				///==========Calculates partial summation over angles==============
				for (int z = 0; z < nbm_size; ++z)
				{
					double cfn  = nu_E [z] * dcos_tet[j];
					double cfan = anu_E[z] * dcos_tet[j];

					int idx = i*nbm_size+z;
					/// Each omp thread writes to different index
					neu_partial_E [idx] = cfn;
					aneu_partial_E[idx] = cfan;

					neu_costet_phi [idx] = cfn  * cos_tet[j];
					aneu_costet_phi[idx] = cfan * cos_tet[j];

					/// Since the following arrays are tmp, we use +=
					/// then they assigned to their corresponding array in the next loop
					nu_sintet [z] += cfn  * sin_tet[j];
					anu_sintet[z] += cfan * sin_tet[j];

				}
				///================================================================

			}	/// End of theta angle loop

			for (int z = 0; z < nbm_size; ++z)
			{
				int idx = i*nbm_size+z;

				neu_sintet_sinphi [idx] = nu_sintet [z] * sin_phi[i];
				aneu_sintet_sinphi[idx] = anu_sintet[z] * sin_phi[i];

				neu_sintet_cosphi [idx] = nu_sintet [z] * cos_phi[i];
				aneu_sintet_cosphi[idx] = anu_sintet[z] * cos_phi[i];
			}
		}	/// End of phi angle loop

		/// Manual reduction!
		Res_t n_part_E = {}, an_part_E = {}, n_part_cos_phi = {}, an_part_cos_phi = {};
		Res_t n_part_sin_sin = {}, an_part_sin_sin = {}, n_part_sin_cos = {}, an_part_sin_cos = {};
		for (int p = 0; p < phibins; ++p)
			for (int z = 0; z < nbm_size; ++z)
			{
				int idx = p*nbm_size+z;

				n_part_E [z] += neu_partial_E [idx];
				an_part_E[z] += aneu_partial_E[idx];

				n_part_cos_phi [z] += neu_costet_phi [idx];
				an_part_cos_phi[z] += aneu_costet_phi[idx];

				n_part_sin_sin [z] += neu_sintet_sinphi [idx];
				an_part_sin_sin[z] += aneu_sintet_sinphi[idx];

				n_part_sin_cos [z] += neu_sintet_cosphi [idx];
				an_part_sin_cos[z] += aneu_sintet_cosphi[idx];
			}

		double neu_coeff  = Jq(neu);//L_E[neu];
		double aneu_coeff = Jq(aneu);//L_E[aneu];

		upd_nu_coef(n_part_E, an_part_E, neu_coeff, aneu_coeff, part_nu);
		upd_nu_coef(n_part_cos_phi, an_part_cos_phi, neu_coeff, aneu_coeff, part_nu_costet_phi);
		upd_nu_coef(n_part_sin_sin, an_part_sin_sin, neu_coeff, aneu_coeff, part_nu_sintet_sinphi);
		upd_nu_coef(n_part_sin_cos, an_part_sin_cos, neu_coeff, aneu_coeff, part_nu_sintet_cosphi);
	}	/// End of flavor loop

}

void calc_Hvv(const NBeam* beam, const int step_num) throw()
{
	/// reset the tmp
	memset(partial_nu, 0., nbm_size*sizeof(double));
	memset(partial_nu_costet_phi, 0., nbm_size*sizeof(double));
	memset(partial_nu_sintet_sinphi, 0., nbm_size*sizeof(double));
	memset(partial_nu_sintet_cosphi, 0., nbm_size*sizeof(double));

	getHvv_partial(beam, partial_nu, partial_nu_costet_phi, partial_nu_sintet_sinphi, partial_nu_sintet_cosphi);
}

void evolve(const NBeam *REST ibeam, const double dr, const double hmatt, const int p1, const int p2, bool calc, const int q, NBeam *REST obeam) throw()
{

#ifdef OMP
#	pragma omp parallel for
#endif
	for (int tet = 0; tet < tetbins; ++tet)
	{
		double dl = dr / cos_tet[tet];

		if(calc) /// if already calculated, no need to recalculate
			getHvv_partial_tet(tet, partial_nu_costet_phi, partial_nu_sintet_sinphi, partial_nu_sintet_cosphi,
					&Hvv_costet_phi_cos[tet*nbm_size], &Hvv_sintet_sinphi_sin[tet*nbm_size], &Hvv_sintet_cosphi_sin[tet*nbm_size]);

		for (int phi = 0; phi < phibins; ++phi)
		{
			if(calc) /// if already calculated, no need to recalculate
				getHvv(phi, partial_nu, &Hvv_costet_phi_cos[tet*nbm_size], &Hvv_sintet_sinphi_sin[tet*nbm_size], &Hvv_sintet_cosphi_sin[tet*nbm_size],
						&Hvv[(phi*tetbins+tet)*nbm_size]);

			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (phi*tetbins+tet)*nbm_prtcls+prc;
				ibeam[idx].evolveBins(prc, dl, &Hvv[(phi*tetbins+tet)*nbm_size], hmatt, obeam[idx]);
			}
		}
	}
}

void avgBeam(const NBeam *REST ibeam, NBeam *REST obeam) throw()
{
#ifdef OMP
#	pragma omp parallel for// collapse(3)
#endif
	for (int tet = 0; tet < tetbins; ++tet)
	{
		for (int phi = 0; phi < phibins; ++phi)
		{
			for (int prc = 0; prc < nbm_prtcls; ++prc)
			{
				int idx = (phi*tetbins+tet)*nbm_prtcls+prc;
				obeam[idx].addAvg(ibeam[idx]);
			}
		}
	}
}

} /// Phy namespace

} /// NBGroup namespace

