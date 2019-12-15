/*
 * NBeam.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifndef NBEAM_H_
#define NBEAM_H_

#include "global.h"

namespace nbm {

/// Place holder for temp results
const int nbm_size = 4;
typedef double Res_t[nbm_size];

extern int nbm_flvs;	///< flavor number
extern int nbm_prtcls;	///< particle number = flv x 2
extern int nbm_cmpn;	///< psi component number = flv x 2
extern int nbm_ebins;	///< number of energy bins

double& Lv(int i);
double& Ev(int i);
double& Tv(int i);
double& eta(int i);

/**
 * Initializes h_vacc values - needs to be called once
 */
void init(int flavors, int ebins);

/**
 * Free h_vacc and hamiltonian arrays
 */
void freemem();

/**
 * Multiply nu and anu mat by their coefficients
 * @nu mat of neu
 * @anu mat of aneu
 * @n_cf neu coef.
 * @an_cf aneu_coef.
 * @ret result
 */
inline void upd_nu_coef( const double *REST nu, const double *REST anu, const double n_cf, const double an_cf, double *REST ret ) throw() ///< Must be inlined!
{
	ret[0] += nu[0] * n_cf - anu[0] * an_cf;
	ret[1] += nu[1] * n_cf - anu[1] * an_cf; ///< they are = 0.;
	ret[2] += nu[2] * n_cf - anu[2] * an_cf;
	ret[3] += nu[3] * n_cf + anu[3] * an_cf;
}

/// ============== NBeam class ================

class NBeam
{
protected:

	double *fv;					///< pointer to energy spectra
	double *ar, *ai, *br, *bi;	///< different component of the wave functions
	NBeam ()
	{ fv = ar = ai = br = bi = 0; hSumMat[0] = hSumMat[1] = hSumMat[2] = hSumMat[3] = Err = 0.;};

	Res_t hSumMat;
	double Err;

public:

	NBeam (int prtcl);
	~NBeam();

	/// Setters/Getters
	inline double& Ar(const int e) { return ar[e]; }
	inline double& Ai(const int e) { return ai[e]; }
	inline double& Br(const int e) { return br[e]; }
	inline double& Bi(const int e) { return bi[e]; }
	inline const double& Ar(const int e) const { return ar[e]; }
	inline const double& Ai(const int e) const { return ai[e]; }
	inline const double& Br(const int e) const { return br[e]; }
	inline const double& Bi(const int e) const { return bi[e]; }

	inline double* Ar()	{ return ar; }
	inline double* Ai()	{ return ai; }
	inline double* Br()	{ return br; }
	inline double* Bi()	{ return bi; }
	inline const double* Ar() const	{ return ar; }
	inline const double* Ai() const	{ return ai; }
	inline const double* Br() const	{ return br; }
	inline const double* Bi() const	{ return bi; }

	inline double Fv(int e) const { return fv[e]; }
	inline double* Fv() const { return fv; }

	/// Return components of the wave function
	inline const double* psi(int cmpn) const throw()
	{
		switch (cmpn) {
			case 0:
				return ar;
			case 1:
				return ai;
			case 2:
				return br;
			case 3:
				return bi;
			default:
				return 0;
		}
	}

	inline double* psi(int cmpn) throw()
	{
		switch (cmpn) {
			case 0:
				return ar;
			case 1:
				return ai;
			case 2:
				return br;
			case 3:
				return bi;
			default:
				return 0;
		}
	}

	/**
	 * Calculates density matrix for enerybin = ebin
	 * @ebin current energy bin
	 * @ret return matrix
	 */
	inline void density(const int ebin, Res_t ret) const throw();
	inline void density(const double ar, const double ai, const double br, const double bi, Res_t ret) const throw();

	/**
	 * Evolves an neutrino bins with Hamiltonian
	 * h_* are hamiltonian matrix components, dr is delta r
	 */
	inline void U(	const int e, const double dr, const double h_r0, const double h_i0, const double h_r1, const double h_i1,
					const double n_ar, const double n_ai, const double n_br, const double n_bi) throw();

	inline void U(	const double dr, const double h_r0, const double h_i0, const double h_r1, const double h_i1,
					const double a_r, const double a_i, const double b_r, const double b_i,
					double& n_ar, double& n_ai, double& n_br, double& n_bi) const throw();

	/**
	 * Calculates the summation over energy bins
	 */
	void calcESum() throw();

	/**
	 * Returns the summation over energy bins
	 * @ret the return value
	 */
	void getESum( Res_t ret ) const throw();

	/**
	 * Evolves nu arrays with hamiltonian matrix
	 * @ptc_idx indicates type of particle
	 * @dr delta radius
	 * @beam the result will be saved into this array
	 */
	void evolveBinsAvgErr( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg, NBeam& beamErr ) throw();
	void evolveBinsHvvAvg( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw();
	void evolveBinsAvg( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw();
	void evolveBinsAvg( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beamAvg ) throw();
	void evolveBinsHvv( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw();
	void evolveBins( const NBeam& beam, const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw();
	void evolveBins( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt ) throw();

	void evolveBins( const int ptc_idx, const double dr, const double *REST hvv, const double hmatt, NBeam& beam ) const throw();

	/**
	 * Add the value of passed NBeam's bins to the current bins and take the average
	 * @beam the second NBeam's bins to be added to the current bins
	 */
	void addAvg( const NBeam& beam ) throw();

	/**
	 * Calculates maximum err between energy bins of two Beam
	 * @beam the second NBeam's bins to be taken difference
	 * @return the maximum err
	 */
	double calcErr( const NBeam& beam ) const throw();
	double calcErrHvv( const NBeam& beam ) throw();
	double& getErr() throw();
};

}

#endif /* NBEAM_H_ */
