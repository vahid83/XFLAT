/*
 * NBeam.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */
#include <cstdlib>
#include <cmath>

#include "global.h"
#include "Util.h"
#include "NBGroup.h"
#include "Fenergy.h"

#include "NBeam.h"
/// if we want to add perturbation to wavefunctions!
//#define	PERTERB

namespace nbm {

///===================== Variables ========================

int nbm_flvs;	///< flavor number
int nbm_prtcls;	///< particle number = flv x 2
int nbm_cmpn;	///< Psi component number = flv x 2
int nbm_ebins;	///< place holder for number of energy bins
double* h_vacc; ///< an array of H0's, needs to calculate once -- size= ebin x 2

//double dE;//, dE2;

/// Energy spectra arrays
double* fq_ne;
double* fq_ane;
double* fq_nt;
double* fq_ant;

/// energy luminosity
double Lve, Lvae, Lvt, Lvat;

/// average energy
double Ev_e, Ev_ae, Ev_t, Ev_at;

/// Temperature
double Tve, Tvae, Tvt, Tvat;

/// eta
double eta_ve, eta_vae, eta_vt, eta_vat;

///=================== Implementation =====================

inline double* Fq(int prtcl_idx) throw()
{
	switch (prtcl_idx) {
		case 0:
			return fq_ne;
		case 1:
			return fq_ane;
		case 2:
			return fq_nt;
		case 3:
			return fq_ant;
		default:
			return NULL;
	}
}

double error = -4.;
double& Lv(int i)
{
	switch(i)
	{
		case 0:
			return Lve;
		case 1:
			return Lvae;
		case 2:
			return Lvt;
		case 3:
			return Lvat;
		default:
			return error;
	}
}

double& Ev(int i)
{
	switch(i)
	{
		case 0:
			return Ev_e;
		case 1:
			return Ev_ae;
		case 2:
			return Ev_t;
		case 3:
			return Ev_at;
		default:
			return error;
	}
}

double& Tv(int i)
{
	switch(i)
	{
		case 0:
			return Tve;
		case 1:
			return Tvae;
		case 2:
			return Tvt;
		case 3:
			return Tvat;
		default:
			return error;
	}
}

double& eta(int i)
{
	switch(i)
	{
		case 0:
			return eta_ve;
		case 1:
			return eta_vae;
		case 2:
			return eta_vt;
		case 3:
			return eta_vat;
		default:
			return error;
	}
}

///=================== implementations ====================

void init(int flavors, int ebins)
{
	/// init the energy module
	fed::init(ebins);

	/// Calculates the average evergy for each particle
	Ev_e  = fed::avgEnergy(0);
	Ev_ae = fed::avgEnergy(1);
	Ev_t  = fed::avgEnergy(2);
	Ev_at = fed::avgEnergy(3);	

	nbm_flvs = flavors;
	nbm_prtcls = flavors << 1;
	nbm_cmpn = nbm_prtcls;
	nbm_ebins = ebins;

	posix_memalign((void**)&fq_ne,  ALIGN_LEN, ebins*sizeof(double));
	posix_memalign((void**)&fq_ane, ALIGN_LEN, ebins*sizeof(double));
	posix_memalign((void**)&fq_nt,  ALIGN_LEN, ebins*sizeof(double));
	posix_memalign((void**)&fq_ant, ALIGN_LEN, ebins*sizeof(double));

	/// For each energy bin we need to keep two double for its hamiltonian
	posix_memalign((void**)&h_vacc, ALIGN_LEN, 2*ebins*sizeof(double));

	double sum_ne = 0., sum_ane = 0., sum_nt = 0., sum_ant = 0.;
	double frst = Util::E0() + fed::dE*.5;
#ifdef	SIMD
//#	pragma omp simd aligned(fq_ne, fq_ane, fq_nt, fq_ant, h_vacc: ALIGN_LEN) reduction(+: sum_ne, sum_ane, sum_nt, sum_ant)
#endif
	for (int e = 0; e < ebins; ++e)
	{
		double q = frst + e * fed::dE;
		///energy spectra pre-calculations
		fq_ne [e] = fed::fq(q, Tve, eta_ve);
		fq_ane[e] = fed::fq(q, Tvae, eta_vae);
		fq_nt [e] = fed::fq(q, Tvt, eta_vt);
		fq_ant[e] = fed::fq(q, Tvat, eta_vat);
		sum_ne  += fq_ne [e];
		sum_ane += fq_ane[e];
		sum_nt  += fq_nt [e];
		sum_ant += fq_ant[e];
		/// H_vaccume calculations
		double delta = .5 * Util::dm2() / q;
		int j = e << 1;
		double theta = Util::theta2();	///< theta x 2 of Vaccume mixing angle
		h_vacc[j  ] = -cos(theta) * delta * .5;	///< there is a 1/2 Hamiltonian factor
		h_vacc[j+1] =  sin(theta) * delta * .5;	///< " "
		//h_vacc[3] = -h_vacc[0];
	}
//	sum_ne  *= fed::dE; //< excluded here
//	sum_ane *= fed::dE;
//	sum_nt  *= fed::dE;
//	sum_ant *= fed::dE;

	/// normalizing fv functions
	sum_ne  = 1. / sum_ne;
	sum_ane = 1. / sum_ane;
	sum_nt  = 1. / sum_nt;
	sum_ant	= 1. / sum_ant;

#ifdef	SIMD
#	pragma omp simd aligned(fq_ne, fq_ane, fq_nt, fq_ant: ALIGN_LEN)
#endif
	for (int e = 0; e < ebins; ++e)
	{
		fq_ne [e] *= sum_ne ;//* fed::dE; //< and excluded here
		fq_ane[e] *= sum_ane;// * fed::dE;
		fq_nt [e] *= sum_nt ;// * fed::dE;
		fq_ant[e] *= sum_ant;// * fed::dE;

	}
}

void freemem()
{
	free(h_vacc);

	free(fq_ne);
	free(fq_ane);
	free(fq_nt);
	free(fq_ant);
}

/**
 * Returns H0(ebin)
 * @ebin the current energy bin
 * @ret return H0
 */
inline void getH0(	const int ebin, ///< input
					double ret[2])	///< output
{
	int i  = ebin << 1;		///< x2
	ret[0] = h_vacc[i  ];	///< h00
	ret[1] = h_vacc[i+1];	///< h01
}

//========== NBeam implementations =============

NBeam::NBeam(int prtc)
{
	fv = Fq(prtc);
	hSumMat[0] = hSumMat[1] = hSumMat[2] = hSumMat[3] = Err = 0.;

	posix_memalign((void**)&ar, ALIGN_LEN, nbm_ebins*sizeof(double));
	posix_memalign((void**)&ai, ALIGN_LEN, nbm_ebins*sizeof(double));
	posix_memalign((void**)&br, ALIGN_LEN, nbm_ebins*sizeof(double));
	posix_memalign((void**)&bi, ALIGN_LEN, nbm_ebins*sizeof(double));

	double a_r = (prtc < 2) ? 1. : 0.;
	double a_i = 0.;
	double b_r = (prtc < 2) ? 0. : 1.;
	double b_i = 0.;

	double* ar = this->ar;
	double* ai = this->ai;
	double* br = this->br;
	double* bi = this->bi;
#ifdef	SIMD
//#	pragma omp simd aligned(ar, ai, br, bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
#ifdef	PERTERB
		double e1 = rand() * (1.e-8) - (.5*1.e-8);
		double e2 = rand() * (1.e-8) - (.5*1.e-8);

		double eps = (double( e1 )/RAND_MAX);
		double eta = (double( e2 )/RAND_MAX);
#else
		double eps = 0.;
		double eta = 0.;
#endif
		double one_sq = sqrt(1 - eps*eps - eta*eta);

		ar[e] = (prtc < 2) ? one_sq : eps;
		ai[e] = 0.;
		br[e] = (prtc < 2) ? eps : one_sq;
		bi[e] = eta;
	}
}

NBeam::~NBeam()
{
	free(ar);
	free(ai);
	free(br);
	free(bi);
}

///	for particle the matrix is:
///  __                              __
/// | res0,res1(=0)    res2,res3       |
/// | res2,-res3       -res0,-res1(=0) |
/// |__                              __|
///
///	the conjugate of anti-particle is:
///  __                              __
/// | res0,res1(=0)    res2,-res3      |
/// | res2,res3        -res0,-res1(=0) |
/// |__                              __|
///
void NBeam::density(const int e, 				///< inputs
					Res_t ret) const throw()	///< output
{
	const double *REST ar = this->ar;
	const double *REST ai = this->ai;
	const double *REST br = this->br;
	const double *REST bi = this->bi;

	/// A(0,0): (r,i)
	ret[0] = .5 * ( Util::norm2(ar[e],ai[e]) - Util::norm2(br[e],bi[e]) );
	ret[1] = 0.;	///< if this line is commented out, make sure that res = {0.} before call to this function!
	/// A(0,1): (r,i)
	Util::mul_cmplx(ar[e], ai[e], br[e], -bi[e], ret[2], ret[3]);
}

void NBeam::density(const double ar, const double ai, const double br, const double bi,  ///< inputs
                                        Res_t ret) const throw()        ///< output
{
	/// A(0,0): (r,i)
        ret[0] = .5 * ( Util::norm2(ar,ai) - Util::norm2(br,bi) );
	ret[1] = 0.;    ///< if this line is commented out, make sure that res = {0.} before call to this function!
	/// A(0,1): (r,i)
	Util::mul_cmplx(ar, ai, br, -bi, ret[2], ret[3]);
}

/// Old style API
void NBeam::U(	const int e, const double dr, const double h_r0, const double h_i0, const double h_r1, const double h_i1,
				const double n_ar, const double n_ai, const double n_br, const double n_bi) throw()
{
	double mat0, mat1, mat2, mat3, mat4, mat5, mat6, mat7;

	double lambda = sqrt( h_r0*h_r0 + Util::norm2(h_r1, h_i1) );
	double ldr = lambda * dr;

	double cosCoef = cos(ldr);
	double sinCoef = sin(ldr) / lambda;

	mat0 =  cosCoef;
	mat1 = -h_r0 * sinCoef;

	mat2 =  h_i1 * sinCoef;
	mat3 = -h_r1 * sinCoef;

	mat4 = -mat2;
	mat5 =  mat3;

	mat6 =  mat0;
	mat7 = -mat1;

	double res00_r, res00_i, res01_r, res01_i, res10_r, res10_i, res11_r, res11_i;
	Util::mul_cmplx(mat0, mat1, n_ar, n_ai, res00_r, res00_i);
	Util::mul_cmplx(mat2, mat3, n_br, n_bi, res01_r, res01_i);
	Util::mul_cmplx(mat4, mat5, n_ar, n_ai, res10_r, res10_i);
	Util::mul_cmplx(mat6, mat7, n_br, n_bi, res11_r, res11_i);

	ar[e] = res00_r + res01_r;
	ai[e] = res00_i + res01_i;
	br[e] = res10_r + res11_r;
	bi[e] = res10_i + res11_i;
}

void NBeam::U(	const double dr, const double h_r0, const double h_i0, const double h_r1, const double h_i1,
				const double a_r, const double a_i, const double b_r, const double b_i,
				double& n_ar, double& n_ai, double& n_br, double& n_bi	) const	throw()
{
	double mat0, mat1, mat2, mat3, mat4, mat5, mat6, mat7;

	double lambda = sqrt( h_r0*h_r0 + Util::norm2(h_r1, h_i1) );
	double ldr = lambda * dr;
	double cosCoef = cos(ldr);
	double sinCoef = sin(ldr) / lambda;

	mat0 =  cosCoef;
	mat1 = -h_r0 * sinCoef;

	mat2 =  h_i1 * sinCoef;
	mat3 = -h_r1 * sinCoef;

	mat4 = -mat2;
	mat5 =  mat3;

	mat6 =  mat0;
	mat7 = -mat1;

	double res00_r, res00_i, res01_r, res01_i, res10_r, res10_i, res11_r, res11_i;
	Util::mul_cmplx(mat0, mat1, a_r, a_i, res00_r, res00_i);
	Util::mul_cmplx(mat2, mat3, b_r, b_i, res01_r, res01_i);
	Util::mul_cmplx(mat4, mat5, a_r, a_i, res10_r, res10_i);
	Util::mul_cmplx(mat6, mat7, b_r, b_i, res11_r, res11_i);

	n_ar = res00_r + res01_r;
	n_ai = res00_i + res01_i;
	n_br = res10_r + res11_r;
	n_bi = res10_i + res11_i;
}

void NBeam::evolveBinsAvgErr( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg, NBeam& beamErr ) throw()
{
	const double *REST fv = this->fv;

	double *REST ar = this->ar;
	double *REST ai = this->ai;
	double *REST br = this->br;
	double *REST bi = this->bi;

	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	double *REST bmav_ar = beamAvg.ar;
	double *REST bmav_ai = beamAvg.ai;
	double *REST bmav_br = beamAvg.br;
	double *REST bmav_bi = beamAvg.bi;

	const double *REST bmer_ar = beamErr.ar;
	const double *REST bmer_ai = beamErr.ai;
	const double *REST bmer_br = beamErr.br;
	const double *REST bmer_bi = beamErr.bi;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

	double eps = 0.;

#ifdef SIMD
#	pragma omp simd reduction(+: eps) aligned(ar,ai,br,bi, bm_ar,bm_ai,bm_br,bm_bi, bmav_ar,bmav_ai,bmav_br,bmav_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		//U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
		//	bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                        bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e],
                        ar[e], ai[e], br[e], bi[e] );

///---------------Calc avg.--------------------
		bmav_ar[e] = (ar[e] + bmav_ar[e]) * .5;
		bmav_ai[e] = (ai[e] + bmav_ai[e]) * .5;
		bmav_br[e] = (br[e] + bmav_br[e]) * .5;
		bmav_bi[e] = (bi[e] + bmav_bi[e]) * .5;

///--------------Calc err.---------------------
		eps += ( Util::norm2(bmer_ar[e] - bmav_ar[e], bmer_ai[e] - bmav_ai[e]) + Util::norm2(bmer_br[e] - bmav_br[e], bmer_bi[e] - bmav_bi[e]) ) * (fv[e]);
	}

	beamErr.getErr() = eps;

}

void NBeam::evolveBinsAvg( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw()
{
	double *REST ar = this->ar;
	double *REST ai = this->ai;
	double *REST br = this->br;
	double *REST bi = this->bi;

	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	double *REST bmav_ar = beamAvg.ar;
	double *REST bmav_ai = beamAvg.ai;
	double *REST bmav_br = beamAvg.br;
	double *REST bmav_bi = beamAvg.bi;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

#ifdef SIMD
#	pragma omp simd aligned(ar,ai,br,bi, bm_ar,bm_ai,bm_br,bm_bi, bmav_ar,bmav_ai,bmav_br,bmav_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		//U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
		//	bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                        bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e],
                        ar[e], ai[e], br[e], bi[e] );

///---------------Calc avg.--------------------
		bmav_ar[e] = (ar[e] + bmav_ar[e]) * .5;
		bmav_ai[e] = (ai[e] + bmav_ai[e]) * .5;
		bmav_br[e] = (br[e] + bmav_br[e]) * .5;
		bmav_bi[e] = (bi[e] + bmav_bi[e]) * .5;

	}

}

void NBeam::evolveBinsHvvAvg( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw()
{
	double *REST ar = this->ar;
	double *REST ai = this->ai;
	double *REST br = this->br;
	double *REST bi = this->bi;

	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	double *REST bmav_ar = beamAvg.ar;
	double *REST bmav_ai = beamAvg.ai;
	double *REST bmav_br = beamAvg.br;
	double *REST bmav_bi = beamAvg.bi;

	const double *REST fv = this->fv;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

	/// Reset the summation result
	double r0 = 0., r1 = 0., r2 = 0., r3 = 0.;


#ifdef SIMD
#	pragma omp simd reduction(+: r0,r1,r2,r3) aligned(ar,ai,br,bi, bm_ar,bm_ai,bm_br,bm_bi, bmav_ar,bmav_ai,bmav_br,bmav_bi, fv: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		//U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
		//	bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                        bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e],
                        ar[e], ai[e], br[e], bi[e] );

///---------------Calc. Hvv-------------------
		Res_t neu_dens;
		density(ar[e], ai[e], br[e], bi[e], neu_dens);

		r0 += neu_dens[0] * fv[e];
		r1 += neu_dens[1] * fv[e];	///< it's = 0.;
		r2 += neu_dens[2] * fv[e];
		r3 += neu_dens[3] * fv[e];

///---------------Calc avg.--------------------
		bmav_ar[e] = (ar[e] + bmav_ar[e]) * .5;
		bmav_ai[e] = (ai[e] + bmav_ai[e]) * .5;
		bmav_br[e] = (br[e] + bmav_br[e]) * .5;
		bmav_bi[e] = (bi[e] + bmav_bi[e]) * .5;
	}

	hSumMat[0] = r0;
	hSumMat[1] = r1;	///< zero!
	hSumMat[2] = r2;
	hSumMat[3] = r3;
}

void NBeam::evolveBinsHvv( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw()
{
	double *REST ar = this->ar;
        double *REST ai = this->ai;
        double *REST br = this->br;
        double *REST bi = this->bi;

	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	const double *REST fv = this->fv;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

	/// Reset the summation result
	double r0 = 0., r1 = 0., r2 = 0., r3 = 0.;


#ifdef SIMD
#	pragma omp simd reduction(+: r0,r1,r2,r3) aligned(ar,ai,br,bi,bm_ar,bm_ai,bm_br,bm_bi, fv: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		//U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
		//	bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                        beam.ar[e], beam.ai[e], beam.br[e], beam.bi[e],
			ar[e], ai[e], br[e], bi[e] );
///-------------------------------------
		Res_t neu_dens;
		density(ar[e], ai[e], br[e], bi[e], neu_dens);

		r0 += neu_dens[0] * fv[e];
		r1 += neu_dens[1] * fv[e];	///< it's = 0.;
		r2 += neu_dens[2] * fv[e];
		r3 += neu_dens[3] * fv[e];
	}

	hSumMat[0] = r0;
	hSumMat[1] = r1;	///< zero!
	hSumMat[2] = r2;
	hSumMat[3] = r3;
}

void NBeam::evolveBins( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw()														///< outputs
{
	double *REST ar = this->ar;
        double *REST ai = this->ai;
        double *REST br = this->br;
        double *REST bi = this->bi;

	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

#ifdef SIMD
#	pragma omp simd aligned(ar,ai,br,bi,bm_ar,bm_ai,bm_br,bm_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		//U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
		//	bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
			bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e],
			ar[e], ai[e], br[e], bi[e] );

	}
}

void NBeam::evolveBinsAvg( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw()
{
	double *REST ar = this->ar;
        double *REST ai = this->ai;
        double *REST br = this->br;
        double *REST bi = this->bi;

	double *REST bmav_ar = beamAvg.ar;
        double *REST bmav_ai = beamAvg.ai;
        double *REST bmav_br = beamAvg.br;
        double *REST bmav_bi = beamAvg.bi;

        double hmvv[4] AlignAs(ALIGN_LEN);
        /// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
        hmvv[1] = hvv[1];
        hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
        hmvv[3] = hvv[3];

#ifdef SIMD
#       pragma omp simd aligned(ar,ai,br,bi, bmav_ar,bmav_ai,bmav_br,bmav_bi: ALIGN_LEN)
#endif
        for (int e = 0; e < nbm_ebins; ++e)
        {
                double h0[2] AlignAs(ALIGN_LEN);
                getH0(e, h0);

                double hamilt[4] AlignAs(ALIGN_LEN);
                hamilt[0] = h0[0] + hmvv[0];
                hamilt[1] = hmvv[1];
                hamilt[2] = h0[1] + hmvv[2];
                hamilt[3] = hmvv[3];

                //U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                //        bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
                U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
                        ar[e], ai[e], br[e], bi[e],
                        ar[e], ai[e], br[e], bi[e] );

///---------------Calc avg.--------------------
		bmav_ar[e] = (ar[e] + bmav_ar[e]) * .5;
                bmav_ai[e] = (ai[e] + bmav_ai[e]) * .5;
                bmav_br[e] = (br[e] + bmav_br[e]) * .5;
                bmav_bi[e] = (bi[e] + bmav_bi[e]) * .5;
        }

}

void NBeam::evolveBins( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw()														///< outputs
{
	double *REST bm_ar = this->ar;
	double *REST bm_ai = this->ai;
	double *REST bm_br = this->br;
	double *REST bm_bi = this->bi;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

#ifdef SIMD
#	pragma omp simd aligned(bm_ar,bm_ai,bm_br,bm_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		U( e, dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
			bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
	}
}

/// Old style API
void NBeam::evolveBins( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt,	///< inputs
						NBeam& beam ) const	throw()														///< outputs
{
	const double *REST ar = this->ar;
	const double *REST ai = this->ai;
	const double *REST br = this->br;
	const double *REST bi = this->bi;
	double *REST bm_ar = beam.ar;
	double *REST bm_ai = beam.ai;
	double *REST bm_br = beam.br;
	double *REST bm_bi = beam.bi;

	double hmvv[4] AlignAs(ALIGN_LEN);
	/// the (ptc_idx%2) determines particle/anti-particle
	hmvv[0] = (ptc_idx%2) ? -(hvv[0] + hmatt) : (hvv[0] + hmatt);
	hmvv[1] = hvv[1];
	hmvv[2] = (ptc_idx%2) ? -hvv[2] : hvv[2];
	hmvv[3] = hvv[3];

#ifdef SIMD
#	pragma omp simd aligned(ar,ai,br,bi, bm_ar,bm_ai,bm_br,bm_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		double h0[2] AlignAs(ALIGN_LEN);
		getH0(e, h0);

		double hamilt[4] AlignAs(ALIGN_LEN);
		hamilt[0] = h0[0] + hmvv[0];
		hamilt[1] = hmvv[1];
		hamilt[2] = h0[1] + hmvv[2];
		hamilt[3] = hmvv[3];

		U( dr, hamilt[0], hamilt[1], hamilt[2], hamilt[3],
			ar[e], ai[e], br[e], bi[e],
			bm_ar[e], bm_ai[e], bm_br[e], bm_bi[e] );
	}

}

void NBeam::calcESum() throw()
{
	const double *REST ar = this->ar;
	const double *REST ai = this->ai;
	const double *REST br = this->br;
	const double *REST bi = this->bi;

	const double *REST fv = this->fv;

	double r0 = 0., r1 = 0., r2 = 0., r3 = 0.;
#ifdef SIMD
#	pragma omp simd reduction(+: r0,r1,r2,r3) aligned(ar,ai,br,bi,fv: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		Res_t neu_dens;
		density(ar[e],ai[e],br[e],bi[e],neu_dens);

		r0 += neu_dens[0] * fv[e];
		r1 += neu_dens[1] * fv[e];	///< it's = 0.;
		r2 += neu_dens[2] * fv[e];
		r3 += neu_dens[3] * fv[e];
	}

	hSumMat[0] = r0;
	hSumMat[1] = r1;	///< zero!
	hSumMat[2] = r2;
	hSumMat[3] = r3;
}

void NBeam::getESum( Res_t ret ) const throw()	///< output
{
/*
	const double *REST ar = this->ar;
        const double *REST ai = this->ai;
        const double *REST br = this->br;
        const double *REST bi = this->bi;

	const double* fv = this->fv;

	double r0 = 0., r1 = 0., r2 = 0., r3 = 0.;
#ifdef SIMD
#	pragma omp simd reduction(+: r0,r1,r2,r3) aligned(ar,ai,br,bi,fv: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		Res_t neu_dens;
		//density(e, neu_dens);
		density(ar[e],ai[e],br[e],bi[e],neu_dens);

		r0 += neu_dens[0] * fv[e];
		r1 += neu_dens[1] * fv[e];	///< it's = 0.;
		r2 += neu_dens[2] * fv[e];
		r3 += neu_dens[3] * fv[e];
	}

	ret[0] = r0;
	ret[1] = r1;	///< zero!
	ret[2] = r2;
	ret[3] = r3;
*/
	ret[0] = hSumMat[0];
	ret[1] = hSumMat[1];	///< zero!
	ret[2] = hSumMat[2];
	ret[3] = hSumMat[3];

}

void NBeam::addAvg( const NBeam& beam ) throw()
{
	double *REST ar = this->ar;
	double *REST ai = this->ai;
	double *REST br = this->br;
	double *REST bi = this->bi;
	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

#ifdef SIMD
#	pragma omp simd aligned(ar, ai, br, bi, bm_ar, bm_ai, bm_br, bm_bi: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		ar[e] = (ar[e] + bm_ar[e]) * .5;
		ai[e] = (ai[e] + bm_ai[e]) * .5;
		br[e] = (br[e] + bm_br[e]) * .5;
		bi[e] = (bi[e] + bm_bi[e]) * .5;
	}
}

double NBeam::calcErrHvv( const NBeam& beam ) throw()
{
        double *REST ar = this->ar;
        double *REST ai = this->ai;
        double *REST br = this->br;
        double *REST bi = this->bi;
        const double *REST bm_ar = beam.ar;
        const double *REST bm_ai = beam.ai;
        const double *REST bm_br = beam.br;
        const double *REST bm_bi = beam.bi;

        const double *REST fv = this->fv;
        double eps = 0., r0 = 0., r1 = 0., r2 = 0., r3 = 0.;
#ifdef SIMD
#       pragma omp simd reduction(+: eps,r0,r1,r2,r3) aligned(ar, ai, br, bi, bm_ar, bm_ai, bm_br, bm_bi, fv: ALIGN_LEN)
#endif
        for (int e = 0; e < nbm_ebins; ++e)
        {

		eps += ( Util::norm2(ar[e] - bm_ar[e], ai[e] - bm_ai[e]) + Util::norm2(br[e] - bm_br[e], bi[e] - bm_bi[e]) ) * (fv[e]);
		//--------------------------------
		Res_t neu_dens;
		density(ar[e],ai[e],br[e],bi[e],neu_dens);

		r0 += neu_dens[0] * fv[e];
		r1 += neu_dens[1] * fv[e];      ///< it's = 0.;
		r2 += neu_dens[2] * fv[e];
		r3 += neu_dens[3] * fv[e];
        }
	hSumMat[0] = r0;
	hSumMat[1] = r1;        ///< zero!
	hSumMat[2] = r2;
	hSumMat[3] = r3;

        return eps;
}

double NBeam::calcErr( const NBeam& beam ) const throw()
{
	const double *REST ar = this->ar;
	const double *REST ai = this->ai;
	const double *REST br = this->br;
	const double *REST bi = this->bi;
	const double *REST bm_ar = beam.ar;
	const double *REST bm_ai = beam.ai;
	const double *REST bm_br = beam.br;
	const double *REST bm_bi = beam.bi;

	const double* fv = this->fv;

	double eps = 0.;
#ifdef SIMD
#	pragma omp simd reduction(+: eps) aligned(ar, ai, br, bi, bm_ar, bm_ai, bm_br, bm_bi, fv: ALIGN_LEN)
#endif
	for (int e = 0; e < nbm_ebins; ++e)
	{
		eps += ( Util::norm2(ar[e] - bm_ar[e], ai[e] - bm_ai[e]) + Util::norm2(br[e] - bm_br[e], bi[e] - bm_bi[e]) ) * (fv[e]);
	}

	return eps;
}

double& NBeam::getErr() throw()
{
	return Err;
}


}



