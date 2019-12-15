#ifndef GLOBAL_H_
#define GLOBAL_H_

/// Indicates the already plugged in physics module
#ifdef	MAA
#define PHY_MODULE	phy_maa
#elif	CLN
#define PHY_MODULE	phy_cln
#elif	LIN
#define PHY_MODULE	phy_lin
#elif	MA
#define PHY_MODULE	phy_ma
#elif	SA
#define PHY_MODULE	phy_sa
#elif   PLN
#define PHY_MODULE      phy_pln
#else
#define PHY_MODULE	phy_maa
#endif

/// Indicates the already plugged in IO module
#ifdef	IOFI
#define IO_MODULE	io_fi
#elif   IOF
#define IO_MODULE	io_f
#else
#define	IO_MODULE	io_f
#endif

#define FIRST_TOUCH

///Defined in Makefile
//#define SIMD											///< To enable vectorization e.g. AVX
//#define OMP											///< To enable OpenMP mode
//#define MMPI											///< To enable MPI mode

#define AVX_LEN			64								///< Length of alignment for AVX inst.
#define XPHI_LEN		64								///< Length of alignment for Xeon Phi
#define ALIGN_LEN		AVX_LEN

////=============================

const double ERG_S			= 1.053613095584080e+16;	///< For converting erg/s to 1/Km2
const double MEV			= 5.067731161976211e+15;	///< For converting MeV to 1/Km
const double MEV_2			= 3.893793036626876e-32;	///< For converting MeV^-2 to Km2
const double EV2			= 2.568189913006477e+19;	///< For converting eV^2 to 1/Km2

#ifdef WIN32
#define AlignAs(i)			__declspec(align(i))		///< To be replaced by C++11 alignas(i)
#else
#define AlignAs(i)			__attribute__((aligned(i)))	///< To be replaced by C++11 alignas(i)
#endif

#define REST				__restrict

#define DBG

#endif /* GLOBAL_H_ */
