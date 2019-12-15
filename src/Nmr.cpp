/*
 * Nmr.cpp
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
#include <cmath>
#include <sys/time.h>

#include "global.h"
#include "Util.h"
#include "Nmr.h"
#include "NBeam.h"
#include "NBGroup.h"
#include "Fio.h"
#include "Matt.h"
#include "Phy.h"

namespace nbgrp {

#ifdef DBG
double totaltime, timeReduceMax, timeReduceSum, timeBcast;
double ioTime;
#endif

extern int proc_rank;

namespace nmr {

using namespace nbm;

int length; 		///< length of Beam array = phibins x thetabins x particles
const int pnts_num = 3;	///< number of total steps (points) in order to advance

NBeam* nubm_0;
NBeam* nubm_1;
//NBeam* nubm_2;
NBeam* nubm_3;
NBeam* nubm_4;

int init(int len)
{
#ifdef DBG
	totaltime = timeReduceMax= timeReduceSum= timeBcast= 0.;
	ioTime = 0;
#endif
	printf("Allocating memory for beams...\n");
	posix_memalign((void**)&nubm_0, ALIGN_LEN, len*sizeof(NBeam));
	posix_memalign((void**)&nubm_1, ALIGN_LEN, len*sizeof(NBeam));
//	posix_memalign((void**)&nubm_2, ALIGN_LEN, len*sizeof(NBeam));
	posix_memalign((void**)&nubm_3, ALIGN_LEN, len*sizeof(NBeam));
	posix_memalign((void**)&nubm_4, ALIGN_LEN, len*sizeof(NBeam));

	phy::initBeam(nubm_0);
	phy::initBeam(nubm_1);
//	phy::initBeam(nubm_2);
	phy::initBeam(nubm_3);
	phy::initBeam(nubm_4);

	length = len;

	return pnts_num;
}

void freemem()
{
	nbm::freemem();

	phy::freeBeam(nubm_0);
	phy::freeBeam(nubm_1);
//	phy::freeBeam(nubm_2);
	phy::freeBeam(nubm_3);
	phy::freeBeam(nubm_4);

	free(nubm_0);
	free(nubm_1);
//	free(nubm_2);
	free(nubm_3);
	free(nubm_4);
}

int evolutionLoop() throw()
{
#ifdef MMPI
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Node %d entering the evolution loop...\n", proc_rank);
#endif
	double r  = Util::R0();
	double dr = Util::dr();

	NBeam* neuBeam0 = nubm_0;
	NBeam* neuBeam1 = nubm_1;
//	NBeam* neuBeam2 = nubm_2;
	NBeam* neuBeam3 = nubm_3;
	NBeam* neuBeam4 = nubm_4;

	/// Four hamiltonians are needed
	double *h0, *h1, *h3, *h4;
	phy::newHvv(h0);
	phy::newHvv(h1);
	phy::newHvv(h3);
	phy::newHvv(h4);
	
	/// indicates starting point of execution
	timeval start_time, end_time;
        start_time.tv_sec = end_time.tv_sec = 0;
        gettimeofday(&start_time, NULL);

        /// Calculating Hamiltonian for the first array
#ifdef OMP
#       pragma omp parallel for
#endif
        for (int i = 0; i < length; ++i)
                neuBeam0[i].calcESum();

	int counter_tot = 0, counter = 0;
	/// Main Loop!
	while (true)
	{
		/// Check for Termination!
		double break_flag = 0;
		if (!Util::isBench())
			break_flag = ( (r < Util::Rn()) && ((end_time.tv_sec - start_time.tv_sec) < Util::Tn()) && (counter_tot < Util::Ts()) );
		else	/// in benchmark mode the benchmark must always continue for the specific amount of time!
			break_flag = (end_time.tv_sec - start_time.tv_sec) < Util::Tn();
#ifdef MMPI
if (!Util::isBench())
{

		double mpi_buffer[3] AlignAs(ALIGN_LEN) = {break_flag, r, dr};
#ifdef DBG
double t1 = MPI_Wtime();
#endif
		MPI_Bcast(mpi_buffer, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef DBG
double t2 = MPI_Wtime();
timeBcast += t2-t1; 
#endif
		break_flag = mpi_buffer[0];
		r = mpi_buffer[1];
		dr = mpi_buffer[2];
}
#endif
		if (!break_flag) break;

		double dr2 = .5 * dr;

		double rdr = r + dr;
		double rdr2 = r + dr2;

		/// step lengths: [r...r+dr/2...r+dr]
		/// step numbers: [0.......1.......2]
			
		/// pre-compute H_mat
		double hm_r, hm_rdr, hm_rdr2;
		matt::getHm(r, hm_r);		/// at r
		matt::getHm(rdr2, hm_rdr2);	/// at r+dr/2
		matt::getHm(rdr, hm_rdr);	/// at r+dr

		phy::calcAngleBins(r,    0); 	/// at (r)
		phy::calcAngleBins(rdr2, 1); 	/// at (r+dr/2)
		phy::calcAngleBins(rdr,  2); 	/// at (r+dr)

		phy::calcDeltaLs(dr2, 0, 0, 1);	/// calc length=dr/2, start point=0 (r), end point=1 (r+dr/2)
		phy::calcDeltaLs(dr,  1, 0, 2);	/// calc length=dr, start point=0 (r), end point=2 (r+dr)
		phy::calcDeltaLs(dr2, 2, 1, 2);	/// calc length=dr/2, start point=1 (r+dr/2), end point=2 (r+dr)



                phy::calc_Hvv (neuBeam0, 0, h0);				/// S1
                phy::evolveHvv(neuBeam0, 1, h0, hm_r, neuBeam1);		/// S2
                phy::evolveHvv(neuBeam0, 0, h0, hm_r, neuBeam3);		/// S6

		phy::calc_Hvv (neuBeam1, 2, h1);				/// S3
//phy::calc_Hvv (neuBeam3, 1, h3);				/// S7

		//phy::calc_Hvv (neuBeam1, 2, h1, neuBeam3, 1, h3);		/// S3/S7
                phy::evolveAvg(neuBeam0, 1, h1, hm_rdr, neuBeam4, neuBeam1);	/// S4/S5
              
		phy::calc_Hvv (neuBeam3, 1, h3);				/// S7
                phy::evolveHvv(neuBeam0, 1, h3, hm_rdr2, neuBeam4);		/// S8

                phy::calc_Hvv (neuBeam4, 2, h4);				/// S9
                phy::evolveAvg(2, h4, hm_rdr, neuBeam3, neuBeam4);		/// S10/S11

		//evolveAvgErr(neuBeam3, 2, h4, hm_rdr, neuBeam2, neuBeam4, neuBeam1);
		/// =========== calc. max. err. =============
		double max_err = -1.;
#ifdef OMP
#		pragma omp parallel for reduction(max: max_err)
#endif
		for (int i = 0; i < length; ++i)
		{
			//double e = neuBeam4[i].calcErr(neuBeam1[i]);
			double e = neuBeam4[i].calcErrHvv(neuBeam1[i]);
			//double e = neuBeam1[i].getErr();
			max_err = (e > max_err) ? e : max_err;
		}
		max_err = sqrt(max_err);
		double Max_err;
#ifdef MMPI
if (!Util::isBench())
{
#ifdef DBG
double t1 = MPI_Wtime();
#endif
		MPI_Allreduce(&max_err, &Max_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#ifdef DBG
double t2 = MPI_Wtime();
timeReduceMax += t2-t1;
#endif

}
else
		Max_err = max_err;
#else
		Max_err = max_err;
#endif

		if (Max_err < Util::eps0())
		{
			r += dr;
//printf("%f\n",r);
			NBeam* tmp = neuBeam0;
			neuBeam0 = neuBeam4; ///< FixME! if they now point to the same location, can we elimiate the 4?
			neuBeam4 = tmp;
#ifdef DBG
			timeval t0;
			gettimeofday(&t0, NULL);
#endif
			fio::dumpToFile(neuBeam0, counter_tot, r, dr);
#ifdef DBG
                        timeval t1;
                        gettimeofday(&t1, NULL);
			ioTime += double(t1.tv_sec + t1.tv_usec*1.0e-6) - double(t0.tv_sec + t0.tv_usec*1.0e-6);
#endif
			counter++;
		}
		counter_tot++;
		dr *= Util::ksi() * sqrt(Util::eps0() / Max_err);	///< Adjust the delta-r
		dr = ( dr > Util::maxDr() ) ? Util::maxDr() : dr;	///< Trim it to make sure that it's not too large
		dr = ( r+dr > Util::Rn() ) ? Util::Rn() - r : dr;	///< Make sure that it won't go beyond the maximum radius
		gettimeofday(&end_time, NULL);

	} /// End of while loop
#ifdef DBG
	printf("$Rank: %d, brodcast time:%f, reduction_sum time:%f, reduction_max time:%f, total MPI time:%f, IO time:%f total time:%f\n", proc_rank, timeBcast, timeReduceSum, timeReduceMax, timeBcast+timeReduceSum+timeReduceMax, ioTime, double(end_time.tv_sec + end_time.tv_usec*1.0e-6) - double(start_time.tv_sec + start_time.tv_usec*1.0e-6));
#endif
	phy::deleteHvv(h0);
	phy::deleteHvv(h1);
	phy::deleteHvv(h3);
	phy::deleteHvv(h4);
	printf("Total iterations: %d\tCorrect iterations: %d\t at radius(Km): %f\n",counter_tot, counter, r);
	return counter_tot;
}

} /// Nmr namespace

} /// NBGroup namespace



