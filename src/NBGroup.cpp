/*
 * NBGroup.cpp
 *
 *  Created on: Jul 20, 2013
 *      Author: vahid
 */
#ifdef MMPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cstring>
#include "global.h"
#include "NBGroup.h"
#include "Nmr.h"
#include "Fio.h"
#include "Matt.h"
#include "Phy.h"
#include "Util.h"

namespace nbgrp {

/// Get the number of processes
int proc_size;
/// Get the rank of the process
int proc_rank;

void node_benchmark()
{
        Util::isBench() = false; ///< by default we're not benchmarking
        double totalflops = 0., partialflops = 0., benchResult = 0.;

#ifdef MMPI
        MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        /// if the start and end beam indices did not provided, they deteremine with benchmarking
        if (Util::startBeam() < 0 || Util::endBeam() < 0)
        {
                /// Quick benchmark!
                printf("Benchmarking node:#%d ...\n", proc_rank);
		Util::isBench() = true;
                Util::beginBenchmark();
                        nbm::init(Util::Flvs(), 100/*Ebins*/);
                        phy::init();
			matt::init();
                        benchResult = nmr::evolutionLoop();
			phy::freemem();
                Util::endBenchmark();
		Util::isBench() = false;
                printf("End of benchmark: %f ...\n", benchResult);

                MPI_Allreduce(&benchResult, &totalflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Exscan(&benchResult, &partialflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                partialflops = (proc_rank) ? partialflops : 0;
                printf("Rank %d out of %d capable of: %f\n",proc_rank,proc_size,benchResult);
        }
#else
        proc_size = 1;
        proc_rank = 0;
        partialflops = 0.;
        totalflops = benchResult = 1.;
#endif

        /// Distribute NBeams over nodes accordingly
	phy::startBeamIdx() = (Util::startBeam() < 0 || Util::endBeam() < 0) ? int((partialflops / totalflops) * phy::firstDimLen()) : Util::startBeam();
	phy::endBeamIdx()   = (Util::startBeam() < 0 || Util::endBeam() < 0) ? int(( (partialflops+benchResult)/totalflops ) * phy::firstDimLen() -1) : Util::endBeam()-1;

	if (!Util::multiNodeBench())
		return;	///< If we're not in the multi-node benchmaraing modek, no need to go beyond this line!

Util::beginBenchmark();

	int start = phy::startBeamIdx();
        int end = phy::endBeamIdx();
	int Len = phy::firstDimLen();

	const int Num_of_Points = 9;
	double (*bench_res)[Num_of_Points] = new double[Util::maxNodes()-Util::minNodes()+1][Num_of_Points];
	int i = 0;
	for (int n = Util::minNodes(); n <= Util::maxNodes(); ++n, ++i)
	{
	
		double d = -0.04;
		for (int j = 0; j < Num_of_Points; d+=0.01, ++j)
		{
			phy::firstDimLen()  = Len   / n;
			phy::startBeamIdx() = start / n + (proc_rank ? phy::firstDimLen()*d : 0);
			phy::endBeamIdx()   = end   / n + (proc_rank ? 0 : phy::firstDimLen()*d);
			

			nbm::init(Util::Flvs(), Util::Ebins());
			phy::init();
			matt::init();
			bench_res[i][j] = nmr::evolutionLoop();
                        phy::freemem();
		
/*
			printf("Node:%d : totalLen=%d $ start=%d , end=%d # bench/total=%f >> adjustedStart=%d , adjustedEnd=%d # adjustedBench/total=%f\n", 
				proc_rank, 
				Len,
				start, 
				end, 
				benchResult/totalflops, 
				proc_rank ? int(start+Len*d) : start, 
				proc_rank ? end : int(end+Len*d), 
				benchResult/totalflops+(proc_rank ? d : -d)
				);
*/
		}
		printf("----next set of nodes----\n");
	}
	for (int i = 0; i < Util::maxNodes()-Util::minNodes()+1; ++i)
	{
		printf("For %d nodes:\n", Util::minNodes()+i);
		double d = -.04;
		for (int j = 0; j < Num_of_Points; ++j, d+=.01)
			printf("For ratio:%f :: res=%f\n", benchResult/totalflops+(proc_rank ? d : -d), bench_res[i][j]);
	}
	delete[] bench_res;

	phy::startBeamIdx() = start;
	phy::endBeamIdx() = end;
	phy::firstDimLen() = Len;

Util::endBenchmark();
}

/// Initialize other modules
void init()
{
	node_benchmark();

	/// !!! The order of calls is important !!!
	nbm::init(Util::Flvs(), Util::Ebins());	///< it also calls into Fed init function
	phy::init();				///< it also calls into Nmr init function
	fio::init();
	matt::init();

}

void particleLoop()
{
	nmr::evolutionLoop();
}

void freemem()
{
	fio::freemem();
	phy::freemem(); ///< it calls into Nmr freemem function and Nmr calls into nbeam freemem function!
}

}



