#include "global.h"
#include "Util.h"

#ifdef OMP
#include <omp.h>
#endif

namespace Util
{
	bool m_isBench;		///< indicates if we're just benchmarking

	int m_syncstep;		///< interval for sync. with disk
	int m_newFilestep;	///< interval for creating new file

	double m_rstep1;	///< radius interval for doing IO in IO=1 config
	int m_tstep1;		///< time interval for doing IO in IO=1 config
	int m_itrstep1;		///< iteration interval for doing IO in IO=1 config

	double m_rstep2;	///< radius interval for doing IO in IO=2 config
	int m_tstep2;		///< time interval for doing IO in IO=2 config
	int m_itrstep2;		///< iteration interval for doing IO in IO=2 config

	int m_start_beam;	///< starting beam index for the current process node
	int m_end_beam;		///< end beam index and its max value is Abins

	int m_multi_node;	///< true if we're doing a multi node benchmark, otherwise false
	int m_min_nodes;	///< the minimum number of nodes for the multinode benchmark
	int m_max_nodes;	///< the maximum number of nodes for the multinode benchmark

	int m_has_matter;	///< indicates whether or not the matter profile is included

	double m_ksi;		///< empirical constant for estimating the error in evolve
	double m_eps0;		///< the prescribed error tolerance
	double m_mu;		///< ae flux

	int m_Ts;		///< Max. number of computed steps
	int m_Tn;		///< Max. time of execution
	double m_R0;		///< Start Radius
	double m_Rn;		///< Final Radius
	double m_dr;		///< Initial step size
	double m_MaxDr;		///< Maximum step size
	double m_E0;		///< Start energy bin
	double m_E1;		///< Final energy bin

	int m_SPoints;		///< The number of points on the star surface
	int m_Pbins;		///< Phi angle bins number
	int m_Abins;		///< Angle bins number
	int m_Ebins;		///< Energy bins number
	int m_Flvs;		///< Flavors number

	double m_Ye;		///< Electron fraction
	double m_Rv;		///< Radius of NP star in Km
	double m_Mns;		///< NP star mass in unit of solar mass
	double m_gs;		///< Statistical weight in relativistic particles. 11/2
	double m_S;		///< Entropy per Baryon
	double m_hNS;		///< The scale height
	double m_nb0;		///< Initial baryion density

	double m_theta;		///< Vaccume mixing angle
	double m_theta2;	///< Vaccume mixing angle x 2
	double m_dm2;		///< Neutrino mass-squared differences
	double m_Lve;		///< Energy Luminosity for V-el
	double m_Lv_e;		///< Energy Luminosity for anti V-el
	double m_Lvt;		///< Energy Luminosity for V-tau
	double m_Lv_t;		///< Energy Luminosity for anti V-tau

	double m_Tve;		///< Temprature for V-el
	double m_Tv_e;		///< Temprature for anti V-el
	double m_Tvt;		///< Temprature for V-tau
	double m_Tv_t;		///< Temprature for anti V-tau
	double m_eta_ve;	///< Degeneracy parameter for V-el
	double m_eta_v_e;	///< Degeneracy parameter for anti V-el
	double m_eta_vt;	///< Degeneracy parameter for V-tau
	double m_eta_v_t;	///< Degeneracy parameter for anti V-tau

	std::string m_prefix;	///< Dump file prefix name
	int m_dumpMode;		///< verbosity mode
	char* m_inVals;		///< if presented, psi's are init based on file values

///============function impl. section=================
	bool& isBench()				{ return m_isBench;	}

	int& syncstep()				{ return m_syncstep;	}
	int& newFilestep()			{ return m_newFilestep;	}

	double& rstep1()			{ return m_rstep1;	}
	int& tstep1()				{ return m_tstep1;	}
	int& itrstep1()				{ return m_itrstep1;	}
	double& rstep2()			{ return m_rstep2;	}
	int& tstep2()				{ return m_tstep2;	}
	int& itrstep2()				{ return m_itrstep2;	}

	int& startBeam()			{ return m_start_beam;	}
	int& endBeam()				{ return m_end_beam;	}

	int& multiNodeBench()			{ return m_multi_node;	}
	int& minNodes()				{ return m_min_nodes;	}
	int& maxNodes()				{ return m_max_nodes;	}

	int& hasMatter()			{ return m_has_matter;	}

	double& eps0()				{ return m_eps0;	}
	double& ksi()				{ return m_ksi;		}
	double& mu()				{ return m_mu;		}

	int& Ts()				{ return m_Ts;		}
	int& Tn()				{ return m_Tn;		}
	double& R0()				{ return m_R0;		}
	double& Rn()				{ return m_Rn;		}
	double& dr()				{ return m_dr;		}
	double& maxDr()				{ return m_MaxDr;	}
	double& E0()				{ return m_E0;		}
	double& E1()				{ return m_E1;		}

	int& SPoints()				{ return m_SPoints;	}
	int& Abins()				{ return m_Abins;	}
	int& Pbins()				{ return m_Pbins;	}
	int& Ebins()				{ return m_Ebins;	}
	int& Flvs()				{ return m_Flvs;	}

	double& Ye()				{ return m_Ye;		}
	double& Rv()				{ return m_Rv;		}
	double& Mns()				{ return m_Mns;		}
	double& gs()				{ return m_gs;		}
	double& S()				{ return m_S;		}
	double& hNS()				{ return m_hNS;		}
	double& nb0()				{ return m_nb0;		}

	double& theta()				{ return m_theta;	}
	double& theta2()			{ return m_theta2;	}
	double& dm2()				{ return m_dm2;		}

	std::string& filePrefix()		{ return m_prefix;	}
	int& dumpMode()				{ return m_dumpMode;	}
	char* inVals()				{ return m_inVals;	}
	void inVals(char* s)			{ m_inVals = s;		}

///============benchmark section=================
int t_tmp, s_tmp, d_tmp;
double r_tmp;

void beginBenchmark()
{
	t_tmp = Util::Tn();
	Util::Tn() = 10;	///< total time of the benchmark
	s_tmp = Util::Ts();
	Util::Ts() = 100000;	///< large enough to prevent early termination during benchmark
	r_tmp = Util::Rn();
	Util::Rn() = 1000;	///< large enough to prevent early termination during benchmark
	d_tmp = Util::dumpMode();
	Util::dumpMode() = 0;
}

void endBenchmark()
{
	/// Restore to prev. values
	Util::Tn() = t_tmp;
	Util::Ts() = s_tmp;
	Util::Rn() = r_tmp;
	Util::dumpMode() = d_tmp;
}

inline double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

} /// End of namespace
