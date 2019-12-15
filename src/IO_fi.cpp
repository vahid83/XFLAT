/*
 * Fio.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vahid
 */

#ifdef MMPI
#include <mpi.h>
#endif

#include <netcdf.h>
#include <stdlib.h>
//#include <malloc.h>
#include <cstring>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "global.h"
#include "NBeam.h"
#include "Fenergy.h"
#include "Util.h"
#include "Phy.h"
#include "Fio.h"

#define PCI

namespace nbgrp {

extern int proc_rank, proc_size;

namespace io_fi {

int pbins, start_beam, end_beam;

using namespace std;
using namespace nbm;

///===================== Variables ========================
int tbins; 				///< theta bins number
//int beam_length; 			///< total number of beams -- phixtheta

///-----variables for when there is no IO!-----
double last_epoch_r0;			///< when we don't have IO it keeps the last 'r' that the progress is shown
const double rstep0 = 0.5;		///< in no-IO case, we show the radius every 'rstep0' kilometers

/// TODO! To be removed!
const int Ndims = 6;			///< rad, phi, theta, particle, psi_component, energy bins - FixMe! has to be flexible
size_t start[Ndims], count[Ndims];

///-----verbose1 variables-----
int nc1id, rid, drid, psid, beamid;//, dimids[Ndims];
size_t itr;				///< after some specific number of iterations we're creating a new file
size_t ctr;				///< holds the number of generated files
/// To see if it's IO time!
double last_epoch_r1;
timeval start_epoch1;
/// Remote process section
int nc1id_mic, rid_mic, drid_mic, psid_mic;//, beamid_mic, dimids_mic[Ndims];
size_t* count_mic;			///< The MIC dimensions might be different due to difference in load
int buff_mic_len;			///< Holds the lenght of the mic's buffer
//int tet_mic, phi_mic;
int remote_rank;			///< Holds either the rank of remote CPU for a MIC, or remote MIC for a CPU
///----------------------------

///-----verbose2 variables-----
const int Ndims2 = 3;
int nc2id, r2id, dimids2[Ndims2];
size_t start2[Ndims2], count2[Ndims2];
int* psi_ids2;
/// nc id's
int flat_tot_fid, flat_r_fid, rad_dm_id, r_dm_id;
int psi_id;
int r_id;
/// To see if it's IO time!
double last_epoch_r2;
timeval start_epoch2;
///----------------------------

///-----verbose4 variables-----
/// To see if it's IO time!
double last_epoch_r4;
timeval start_epoch4;
///----------------------------

int cmpn;				///< number of components for each particle i.e. ar,ai,br,bi,etc.
int retval;				///< returned error code
//size_t start_psi;			///< start point for appending data in the array
size_t one;
size_t len_e;//start_e, start_r;
//size_t current;

timeval end_epoch;

///=================== Implementation =====================
#define NCRUN(exe)	{int err=exe;if(err){printf("NetCDF Status msg: %s\nError in %s line %d\n",nc_strerror(err),__FILE__,__LINE__);exit(-10);}}

#ifdef MMPI
bool doneRemoteRanks = false;	///< shows that we've already found the remote ranks, so there is no need to redo it again
/// If there are MIC installed, and the current MPI process is on CPU, it should handle MIC's IO data
void findRemoteMateRank()
{
	if (doneRemoteRanks) return;
	doneRemoteRanks = true;

	/// hostnames are usually in hostname.domain for CPUs and hostname-mic.domain for MICs
	char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

	string proc_name(processor_name);
        stringstream ss_name;
        ss_name << proc_name;
        string s_name;
        getline(ss_name, s_name, '.');		///< take the first part of the hostname (without domain)
        size_t pcn = s_name.find("-mic");	///< check the name to see if it's mic ('-mic' is included)
        bool isMic = (pcn == string::npos) ? false : true;

	char* processor_list = new char[MPI_MAX_PROCESSOR_NAME*proc_size];
        MPI_Allgather(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, processor_list, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);

	remote_rank = proc_rank;

	/// CPUs need not to find any other nodes!
	if (isMic)
	{
		int prevMic = 1;		///< indicates how many MICs (including itself) are installed before the current one!
		while (true)
		{
			--remote_rank;		///< decreasing the rank until the proper CPU rank is found
			if (remote_rank < 0) break;
		        stringstream ss;
		        ss << &processor_list[remote_rank*MPI_MAX_PROCESSOR_NAME];
			string s;
		        getline(ss, s, '.');
			size_t pos = s.find("-mic");
		        if (pos != string::npos)///< Another Mic with lower rank is detected!
        		{
				++prevMic;
                		continue;
		        }

			if (s.substr(0,pos) == s_name.substr(0,pcn))
				continue;
			break;
		}
		remote_rank += prevMic;		///< it maps the first MIC to the first CPU, the second MIC to the second CPU, etc.
	}

	/// Each node must know its mate's rank!
	/// Index shows the rank of the node, and the item shows the remote_rank
	int* remote_rank_list = new int[proc_size];
	MPI_Allgather(&remote_rank, 1, MPI_INT, remote_rank_list, 1, MPI_INT, MPI_COMM_WORLD);

	/// If we're on MIC, we've already found the CPU rank
	if (proc_rank == remote_rank)		///< we're on a CPU -- try to find the MIC rank
	{
	        for (int i = 0; i < proc_size; ++i)
                	if (remote_rank_list[i] == proc_rank && proc_rank != i)
			{
        	                remote_rank = i;
				break;
			}
	}
printf("%d - %d\n",proc_rank, remote_rank);
	
	delete[] remote_rank_list;
	delete[] processor_list;
}
#endif

void init(int file_counter)
{
	pbins = Util::Pbins();	///FixMe! temp solution!
	tbins = phy::beamLen()/(pbins*nbm_prtcls);//end_beam-start_beam+1;
//tbins = Util::SPoints();
//pbins = Util::Abins();


	last_epoch_r0 = Util::R0();	///< when it's no-IO, it shows the progress!

	ctr = file_counter;		///< Keeps the number of files that are generated up to now!
	itr = 0;			///< Keeps the number of file IO jobs
	//beam_length = pbins * tbins;
	//start_psi = start_e = start_r = current = 0;
	
	/// number of psi's components which is the same as number of particles i.e. ar, ai, br, bi..
	cmpn  = nbm_prtcls;
	len_e = nbm_ebins;
//	psi_bins = fv_bins = NULL;

	one = 1;			///< to write one 'r' at each step

	stringstream ss;
	ofstream ofs;

	///--------------------- netCDF file1 init ----------------------
	if (Util::dumpMode() & 1)
	{
		ss.str(std::string());
		ss << Util::filePrefix() << "Snapshot" << ctr << "_" << proc_rank << "_.nc";
		NCRUN( nc_create(ss.str().c_str(), NC_CLOBBER, &nc1id) );
		size_t* dim_length = phy::countDim();
		const int D = phy::getDim();
		string  dim_names[D];
		phy::getDimInfo(dim_names);
		int dimids[D];
		NCRUN( nc_def_dim(nc1id, dim_names[0].c_str(),	NC_UNLIMITED, 	&dimids[0]) );
		for (int i = 1; i < D; ++i)
			NCRUN( nc_def_dim(nc1id, dim_names[i].c_str(),	dim_length[i],	&dimids[i]) );				

		NCRUN( nc_def_var(nc1id, "r", 		NC_DOUBLE,	1, 		&dimids[0],	&rid)	);
		NCRUN( nc_def_var(nc1id, "dr", 		NC_DOUBLE,	1, 		&dimids[0],	&drid)	);
//		NCRUN( nc_def_var(nc1id, "theta_bm",	NC_INT, 	1, 		&dimids[1],	&beamid));
		NCRUN( nc_def_var(nc1id, "psi", 	NC_DOUBLE,	D, 		dimids,		&psid)	);

		NCRUN( nc_put_att_double(nc1id, psid, "E0", 	NC_DOUBLE, 1, &Util::E0())	);
		NCRUN( nc_put_att_double(nc1id, psid, "E1", 	NC_DOUBLE, 1, &Util::E1())	);
		for (int i = 0; i < nbm_flvs; ++i)
		{
			int n = i << 1;
			NCRUN( nc_put_att_double(nc1id, psid, string("Teprature for particle"+i).c_str(),    NC_DOUBLE, 1, &nbm::Tv(i))	 );
			NCRUN( nc_put_att_double(nc1id, psid, string("Teprature for anti-particle"+i).c_str(),    NC_DOUBLE, 1, &nbm::Tv(i+1)));
			NCRUN( nc_put_att_double(nc1id, psid, string("eta for particle"+i).c_str(), NC_DOUBLE, 1, &nbm::eta(i))	);
			NCRUN( nc_put_att_double(nc1id, psid, string("eta for anti-particle"+i).c_str(), NC_DOUBLE, 1, &nbm::eta(i+1)));
		}	
		NCRUN( nc_enddef(nc1id) );

//		/// Writes theta beams indices, from star_beam to end_beam -- info to be used by merge scripts
//		size_t s = 0, c = 1;
//		for (int bm = start_beam; bm <= end_beam; ++bm, s++)
//			NCRUN( nc_put_vara_int(nc1id, beamid, &s, &c, &bm) );
#ifdef MMPI
		findRemoteMateRank();	///< Check to see whether or not should the node recieve data from other nodes

#ifdef __MIC__	///< on MIC
		MPI_Send(dim_length, D, MPI_UNSIGNED_LONG, remote_rank, 0, MPI_COMM_WORLD);
#else		///< on CPU
		if (remote_rank != proc_rank)	///< if there is a mate!
		{

		buff_mic_len = 1;
		int msg_len;
		/// Probe for an incoming message from process zero
		MPI_Status status;
		MPI_Probe(remote_rank, 0, MPI_COMM_WORLD, &status);		
		MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_len);
                count_mic = new size_t[msg_len];
		int* dim_ids = new int[msg_len];
		MPI_Recv(count_mic, msg_len, MPI_UNSIGNED_LONG, remote_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		ss.str(std::string());
                ss << Util::filePrefix() << "Snapshot" << ctr << "_" << remote_rank << "_.nc";
		NCRUN( nc_create(ss.str().c_str(), NC_CLOBBER, &nc1id_mic) );
		NCRUN( nc_def_dim(nc1id_mic, dim_names[0].c_str(),  NC_UNLIMITED,   &dim_ids[0]) );
                for (int i = 1; i < msg_len; ++i)
		{
                        NCRUN( nc_def_dim(nc1id_mic, dim_names[i].c_str(),  count_mic[i],  &dim_ids[i]) );
			buff_mic_len *= count_mic[i];
		}
                NCRUN( nc_def_var(nc1id_mic, "r",           NC_DOUBLE,      1,          &dim_ids[0],     &rid_mic)   );
                NCRUN( nc_def_var(nc1id_mic, "dr",          NC_DOUBLE,      1,          &dim_ids[0],     &drid_mic)  );
          //    NCRUN( nc_def_var(nc1id_mic, "theta_bm",    NC_INT,         1,          &dim_ids[1],     &beamid_mic));
                NCRUN( nc_def_var(nc1id_mic, "psi",         NC_DOUBLE,      D,		dim_ids,         &psid_mic)  );

		NCRUN( nc_put_att_double(nc1id_mic, psid_mic, "E0",     NC_DOUBLE, 1, &Util::E0())      );
                NCRUN( nc_put_att_double(nc1id_mic, psid_mic, "E1",     NC_DOUBLE, 1, &Util::E1())      );
		for (int i = 0; i < nbm_flvs; ++i)
                {
                        int n = i << 1;
                        NCRUN( nc_put_att_double(nc1id_mic, psid_mic, string("Teprature for particle"+i).c_str(),    NC_DOUBLE, 1, &nbm::Tv(i))  );
                        NCRUN( nc_put_att_double(nc1id_mic, psid_mic, string("Teprature for anti-particle"+i).c_str(),    NC_DOUBLE, 1, &nbm::Tv(i+1)));
                        NCRUN( nc_put_att_double(nc1id_mic, psid_mic, string("eta for particle"+i).c_str(), NC_DOUBLE, 1, &nbm::eta(i)) );
                        NCRUN( nc_put_att_double(nc1id_mic, psid_mic, string("eta for anti-particle"+i).c_str(), NC_DOUBLE, 1, &nbm::eta(i+1)));
                }
                NCRUN( nc_enddef(nc1id_mic) );
		delete[] dim_ids;
//		size_t s = 0, c = 1;
//		for (int bm = start_beam; bm <= end_beam; ++bm, s++)
//			NCRUN( nc_put_vara_int(nc1id_mic, beamid_mic, &s, &c, &bm) );
		}
#endif

#endif
		last_epoch_r1 = Util::R0();
		gettimeofday(&start_epoch1, NULL); ///< start the clock for IO intervals
	}
	///-------------------------end of file1-------------------------

	///--------------------- netCDF file2 init ----------------------
	if (Util::dumpMode() & 2)
	{
		psi_ids2 = new int[nbm_prtcls];
	
		ss.str(std::string());
                ss << Util::filePrefix() << "Probability" << ctr << "_" << proc_rank << "_.nc";
                NCRUN( nc_create(ss.str().c_str(), NC_CLOBBER, &nc2id) );

                NCRUN( nc_def_dim(nc2id, "r",		NC_UNLIMITED,	&dimids2[0]) );
		NCRUN( nc_def_dim(nc2id, "theta",       tbins,          &dimids2[1]) );
                NCRUN( nc_def_dim(nc2id, "phi",         pbins,          &dimids2[2]) );

                NCRUN( nc_def_var(nc2id, "r",	NC_DOUBLE,	1,	&dimids2[0], &r2id) );
		for (int i = 0; i < nbm_prtcls; ++i)
		{
			ss.str(std::string());
			ss << "prtcl" << i;
			NCRUN( nc_def_var(nc2id, ss.str().c_str(), NC_DOUBLE, Ndims2, dimids2, &psi_ids2[i] ) );
		}
		NCRUN( nc_enddef(nc2id) );

		start2[0] = start2[1] = start2[2] = 0;
                count2[0] = count2[1] = count2[2] = 1;

		last_epoch_r2 = Util::R0();
		gettimeofday(&start_epoch2, NULL); ///< start the clock for IO intervals
	}
	///-------------------------end of file2-------------------------

	///--------------------- netCDF file4 init ----------------------
	if (Util::dumpMode() & 4)
	{
		last_epoch_r4 = Util::R0();	
		gettimeofday(&start_epoch4, NULL); ///< start the clock for IO intervals
	}
	///-------------------------end of file4------------------------

}

/// Get the init. values from a file and set psi's with them -- if no file, use 0's and 1's
void fillInitData(NBeam *REST nubeam)
{
	int nc_inputid, nc_psid, nc_rid, nc_drid, nc_rdimid;
        size_t rlen;                            ///< number of iterations (length of 'r')

	printf("Provided data file: %s\n", Util::inVals());
        NCRUN( nc_open(Util::inVals(), NC_NOWRITE, &nc_inputid) );
        NCRUN( nc_inq_varid(nc_inputid, "psi", &nc_psid) );
        NCRUN( nc_inq_varid(nc_inputid, "r", &nc_rid) );
        NCRUN( nc_inq_varid(nc_inputid, "dr", &nc_drid) );

        NCRUN( nc_inq_unlimdim(nc_inputid, &nc_rdimid) );

        NCRUN( nc_inq_dimlen(nc_inputid, nc_rdimid, &rlen) );
        size_t st = rlen - 1, cnt = 1;

        double rdr;             ///< the last 'r' and 'dr' values from the previous calculation restored from an input data file
        NCRUN( nc_get_vara_double(nc_inputid, nc_rid, &st, &cnt, &rdr) );
        Util::R0() = rdr;       ///< set the starting radius
        NCRUN( nc_get_vara_double(nc_inputid, nc_drid, &st, &cnt, &rdr) );
        Util::dr() = rdr;

	printf("Loading the provided data file...\n");
        double* buffer;
        posix_memalign((void**) &buffer, ALIGN_LEN, phy::beamLen()*nbm_cmpn*nbm_ebins*sizeof(double));

        size_t* start_phys = phy::startDim();
        size_t* count_phys = phy::countDim();

	size_t s0tmp = start_phys[0];
        start_phys[0] = st; /// Accesing the last set of data in the file

        /// Always only the first dimension is distributed over nodes, so this is the only dimensionality difference from node to node
	size_t s1tmp = start_phys[1];
	start_phys[1] = phy::startBeamIdx();

	NCRUN( nc_get_vara_double(nc_inputid, nc_psid, start_phys, count_phys, buffer) );

        for (size_t i = 0; i < phy::beamLen(); ++i)
        {
                for (size_t c = 0; c < nbm_cmpn; ++c)
                {
                        memcpy( nubeam[i].psi(c), &buffer[(i*nbm_cmpn+c)*nbm_ebins], nbm_ebins*sizeof(double) );
                }
        }
        free(buffer);

        ///Restore starts
        start_phys[0] = s0tmp;
        start_phys[1] = s1tmp;

        NCRUN( nc_close(nc_inputid) );
}

void freemem()
{
	if  (Util::dumpMode() & 1)
	{
		NCRUN( nc_close(nc1id) );
#if defined(MMPI) && !defined(__MIC__)
		if(proc_rank != remote_rank)
		{
			delete[] count_mic;
			NCRUN( nc_close(nc1id_mic) );
		}
#endif
	}
	if  (Util::dumpMode() & 2)
	{
		NCRUN( nc_close(nc2id) );
		delete[] psi_ids2;
	}
}

/// IO mode = 001
void writeToFile(const NBeam *REST nubeam, const double r, const double dr) throw()
{
	/// Create a new file if necessary, when the current file is enough large
	if (itr >= Util::newFilestep())
	{
		freemem();
		init(++ctr); ///< creates a new file with ctr+1 name
	}

	double* buffer;
	posix_memalign((void**) &buffer, ALIGN_LEN, phy::beamLen()*nbm_cmpn*nbm_ebins*sizeof(double));
//--------------
	size_t* count_phys = phy::countDim();
	size_t* start_phys = phy::startDim();
	for (size_t i = 0; i < phy::beamLen(); ++i)
	{

		for (size_t c = 0; c < nbm_cmpn; ++c)
                {
			memcpy( &buffer[(i*nbm_cmpn+c)*nbm_ebins], nubeam[i].psi(c), nbm_ebins*sizeof(double) );
                }

	}
	start_phys[0] = itr;


#ifdef MMPI
#ifdef __MIC__
	MPI_Send(buffer, phy::beamLen()*nbm_cmpn*nbm_ebins, MPI_DOUBLE, remote_rank, 0, MPI_COMM_WORLD);
	/// FixMe! do we need a barrier here?
	free(buffer);
	return;

#else	/// on CPU

#ifdef OMP
#pragma omp parallel
{
#pragma omp single
{
#endif

	if (proc_rank != remote_rank)
	{

#ifdef OMP
#	pragma omp task	/// task 1
{
#endif
		double* buffer_mic;
	        posix_memalign((void**) &buffer_mic, ALIGN_LEN, buff_mic_len*sizeof(double));
		MPI_Recv(buffer_mic, buff_mic_len, MPI_DOUBLE, remote_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		NCRUN( nc_put_vara_double(nc1id_mic, psid_mic, start_phys, count_mic, buffer_mic) );
		NCRUN( nc_put_vara_double(nc1id_mic, rid_mic, &itr, &one, &r));
	        NCRUN( nc_put_vara_double(nc1id_mic, drid_mic, &itr, &one, &dr));
		free(buffer_mic);
#ifdef OMP
}	/// task 1 end
#endif
	}/// if

#ifdef OMP
#       pragma omp task	/// task 2
{
#endif
	NCRUN( nc_put_vara_double(nc1id, psid, start_phys, count_phys, buffer) );
	NCRUN( nc_put_vara_double(nc1id, rid, &itr, &one, &r));
	NCRUN( nc_put_vara_double(nc1id, drid, &itr, &one, &dr));
	free(buffer);
#ifdef OMP
}	/// task 2 end
#pragma omp taskwait
#endif

	
	itr++;

	if (!(itr % Util::syncstep()))	///< sync file every 'n' itr
		NCRUN( nc_sync(nc1id) );

#ifdef OMP
}	/// single end
}	/// parallel end
#endif


#endif	/// __MIC__
#endif	/// MMPI
}

/// IO mode = 010
void writeToFile(const NBeam *REST nubeam, const double r) throw()
{
	double* part_sum = new double[nbm_prtcls]();


	for (int i = 0; i < tbins*pbins; ++i)
	{
		int i0 = i / (pbins);
                int j0 = i % (pbins);
                start2[1] = i0; start2[2] = j0; 
		for (int p = 0; p < nbm_prtcls; ++p)
                {
                                part_sum[p] = 0.;
                                const double* ar = nubeam[(i)*nbm_prtcls+p].psi(0);
                                const double* ai = nubeam[(i)*nbm_prtcls+p].psi(1);
                                const double* fv = nubeam[(i)*nbm_prtcls+p % nbm_flvs].Fv();
				for (int e = 0; e < len_e; ++e)
                                {
                                        part_sum[p] += Util::norm2(ar[e], ai[e]) * fv[e];
                                }
                                NCRUN( nc_put_vara_double(nc2id, psi_ids2[p], start2, count2, &part_sum[p] ) );
		}
	}
/*
	for (int tet = 0; tet < tbins; ++tet)
        {
		start2[1] = tet;
                for (int phi = 0; phi < pbins; ++phi)
                {
			start2[2] = phi;
                        for (int p = 0; p < prtcl; ++p)
			{
				part_sum[p] = 0.;
				const double* ar = nubeam[(tet*pbins+phi)*prtcl+p].psi(0);
				const double* ai = nubeam[(tet*pbins+phi)*prtcl+p].psi(1);
				const double* fv = nubeam[(tet*pbins+phi)*prtcl+p % nbm_flvs].Fv(); ///< FixME! how this line should be in three flv. case
				/// TODO! has to be vectorized!	
				for (int e = 0; e < len_e; ++e)
				{
					part_sum[p] += Util::norm2(ar[e], ai[e]) * fv[e];
				}
				NCRUN( nc_put_vara_double(nc2id, psi_ids2[p], start2, count2, &part_sum[p] ) );				
			}
		}
	}
*/
	NCRUN( nc_put_vara_double(nc2id, r2id, &start2[0], &one, &r));
	start2[0]++;

	delete part_sum;

	if (!(itr % Util::syncstep()))  ///< sync file every 'n' itr
                NCRUN( nc_sync(nc2id) );
}

/// General function which based on verbose mode decides which function to call
void dumpToFile(const NBeam *REST nubeam, const int it, const double r, const double dr, int theta, int phi) throw()
{
	gettimeofday(&end_epoch, NULL);

	if (!Util::dumpMode())	/// if NO IO
	{
		if ( r - last_epoch_r0 > rstep0 )
		{
			last_epoch_r0 = r;
	               	printf("%f\n",r);
		}
	}

	if  (Util::dumpMode() & 1) /// if IO: 001
	{
		if ( r == Util::Rn() || r - last_epoch_r1 > Util::rstep1() || (end_epoch.tv_sec - start_epoch1.tv_sec) > Util::tstep1() || !(it % Util::itrstep1()) )
		{
                	printf("%f\n",r);
                        last_epoch_r1 = r;
                        start_epoch1.tv_sec = end_epoch.tv_sec;
			writeToFile(nubeam, r, dr);
		}
	}
	if  (Util::dumpMode() & 2) /// if IO: 010
	{
		if ( r == Util::Rn() || r - last_epoch_r2 > Util::rstep2() || (end_epoch.tv_sec - start_epoch2.tv_sec) > Util::tstep2() || !(it % Util::itrstep2()) )
                {
                        printf("%f\n",r);
                        last_epoch_r2 = r;
                        start_epoch2.tv_sec = end_epoch.tv_sec;
			writeToFile(nubeam, r);
                }
	}
	if  (Util::dumpMode() & 4) /// if IO: 100
	{
		return;
	}
}


} /// io_fi namespace

} /// NBGroup namespace



