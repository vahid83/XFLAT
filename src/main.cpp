#ifdef MMPI
#include <mpi.h>
#endif

#include "global.h"
#include "Parser.h"
#include "Util.h"
#include "NBGroup.h"
#include "Nmr.h"

int main(int argc, char* argv[])
{

#ifdef MMPI
#ifdef OMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, NULL);
#else
	MPI_Init(&argc, &argv);
#endif
#endif

	Parser parser(argv[1]);

	/// if a data file is presented, save the name which is used in Fio mod.
	Util::inVals( (argc == 3) ? argv[2] : NULL );

	///=======
	nbgrp::init();
	nbgrp::particleLoop();
	nbgrp::freemem();

#ifdef MMPI
	MPI_Finalize();
#endif

	return 0;
}
