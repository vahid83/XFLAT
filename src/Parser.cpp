#include "Parser.h"
#include "global.h"
#include <sstream>
#include <iostream>

#include "Matt.h"
#include "NBeam.h"

using namespace std;

#include "Util.h"

Parser::Parser(const char* configFile) 
{
	ifstream fs;
	fs.open(configFile);
	if ( !fs.is_open() ) 
		cout << "Cannot open the input file...\n";

	readLines(fs);
	
	fs.close();
}

void Parser::readLines(ifstream& fs) 
{	
	string line; ///< Hold the current read line from the file

	while ( getline(fs, line) ) 
	{	
		stringstream ssLine;
		ssLine << line;

		string sFirst;
		ssLine >> sFirst;

		if ( sFirst == "" || sFirst[0] == '#' ) 
			continue;	

		else if ( sFirst == "filePrefix="	)	{ ssLine >> Util::filePrefix();		}
		else if ( sFirst == "dumpMode=" 	) 	{ ssLine >> Util::dumpMode();		}

		else if ( sFirst == "newFile_step="	)	{ ssLine >> Util::newFilestep();	}
		else if ( sFirst == "sync_step="	)	{ ssLine >> Util::syncstep();		}

		else if ( sFirst == "r_step1=" 		)	{ ssLine >> Util::rstep1();		}
		else if ( sFirst == "t_step1=" 		)	{ ssLine >> Util::tstep1();		}
		else if ( sFirst == "itr_step1="	)	{ ssLine >> Util::itrstep1();		}

		else if ( sFirst == "r_step2=" 		)	{ ssLine >> Util::rstep2();		}
		else if ( sFirst == "t_step2=" 		)	{ ssLine >> Util::tstep2();		}
		else if ( sFirst == "itr_step2="	)	{ ssLine >> Util::itrstep2();		}

		else if ( sFirst == "start_beam="	)	{ ssLine >> Util::startBeam();		}
		else if ( sFirst == "end_beam="		)	{ ssLine >> Util::endBeam();		}

		
		else if ( sFirst == "multiNodeBench="	)	{ ssLine >> Util::multiNodeBench();	}
		else if ( sFirst == "minNodes="		)	{ ssLine >> Util::minNodes();		}
		else if ( sFirst == "maxNodes="		)	{ ssLine >> Util::maxNodes();		}

		else if ( sFirst == "hasMatter="	)	{ ssLine >> Util::hasMatter();		}

		else if ( sFirst == "eps0=" 		)	{ ssLine >> Util::eps0();		}
		else if ( sFirst == "ksi=" 		)	{ ssLine >> Util::ksi();		}
		else if ( sFirst == "mu=" 		)	{ ssLine >> Util::mu();			}

		else if ( sFirst == "Tn=" 		)	{ ssLine >> Util::Tn();			}
		else if ( sFirst == "Ts=" 		)	{ ssLine >> Util::Ts();			}
		else if ( sFirst == "R0=" 		)	{ ssLine >> Util::R0();			}
		else if ( sFirst == "Rn=" 		)	{ ssLine >> Util::Rn();			}
		else if ( sFirst == "dr=" 		)	{ ssLine >> Util::dr();			}
		else if ( sFirst == "max_dr=" 		)	{ ssLine >> Util::maxDr();		}
		else if ( sFirst == "E0=" 		)	{ ssLine >> Util::E0();		
								Util::E0() *= MEV;			}
		else if ( sFirst == "E1=" 		)	{ ssLine >> Util::E1();		
								Util::E1() *= MEV;			}

		else if ( sFirst == "SPoints=" 		)	{ ssLine >> Util::SPoints();		}
		else if ( sFirst == "Pbins=" 		)	{ ssLine >> Util::Pbins();		}
		else if ( sFirst == "Abins=" 		)	{ ssLine >> Util::Abins();		}
		else if ( sFirst == "Ebins=" 		)	{ ssLine >> Util::Ebins();		}
		else if ( sFirst == "Flvs=" 		)	{ ssLine >> Util::Flvs();		}
		
		else if ( sFirst == "Ye=" 		)	{ ssLine >> Util::Ye();			}
		else if ( sFirst == "Rv=" 		)	{ ssLine >> matt::matt_Rv();
								Util::Rv() = matt::matt_Rv();		}
		else if ( sFirst == "Mns=" 		)	{ ssLine >> Util::Mns();		}
		else if ( sFirst == "gs=" 		)	{ ssLine >> Util::gs();			}
		else if ( sFirst == "S=" 		)	{ ssLine >> Util::S();			}
		else if ( sFirst == "hNS=" 		)	{ ssLine >> matt::matt_hNS();
								 Util::hNS() = matt::matt_hNS();	}
		else if ( sFirst == "nb0=" 		)	{ ssLine >> matt::matt_nb0();
								matt::matt_nb0() *= 1.E+15;
								Util::nb0() = matt::matt_nb0();		}

		else if ( sFirst == "theta=" 		)	{ ssLine >> Util::theta();	
								Util::theta2() = Util::theta()*2;	}
		else if ( sFirst == "dm2=" 		)	{ ssLine >> Util::dm2();	
								Util::dm2() *= EV2;			}

		else if ( sFirst == "Lve=" 		)	{ ssLine >> nbm::Lv(0);
								nbm::Lv(0) *= ERG_S;			}
		else if ( sFirst == "Lv_e=" 		)	{ ssLine >> nbm::Lv(1);
								nbm::Lv(1) *= ERG_S;			}
		else if ( sFirst == "Lvt=" 		)	{ ssLine >> nbm::Lv(2);
								nbm::Lv(2) *= ERG_S;			}
		else if ( sFirst == "Lv_t=" 		)	{ ssLine >> nbm::Lv(3);
								nbm::Lv(3) *= ERG_S;			}

		else if ( sFirst == "Tve=" 		)	{ ssLine >> nbm::Tv(0);
								nbm::Tv(0) *= MEV;			}
		else if ( sFirst == "Tv_e=" 		)	{ ssLine >> nbm::Tv(1);
								nbm::Tv(1) *= MEV;			}
		else if ( sFirst == "Tvt=" 		)	{ ssLine >> nbm::Tv(2);
								nbm::Tv(2) *= MEV;			}
		else if ( sFirst == "Tv_t=" 		)	{ ssLine >> nbm::Tv(3);
								nbm::Tv(3) *= MEV;			}

		else if ( sFirst == "eta_ve=" 		)	{ ssLine >> nbm::eta(0);		}
		else if ( sFirst == "eta_v_e=" 		)	{ ssLine >> nbm::eta(1);		}
		else if ( sFirst == "eta_vt=" 		)	{ ssLine >> nbm::eta(2);		}
		else if ( sFirst == "eta_v_t="		)	{ ssLine >> nbm::eta(3);		}

	}

}
