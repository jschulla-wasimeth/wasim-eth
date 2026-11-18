/*   $Revision: 21 $ 
	$Date: 16.05.10 0:09 $  */
/*

Rahmenprogramm fuer das Wasserhaushalts- und Flussgebietsmodell
WaSiM-ETH, Version 2 (Program version see file wasim.h) 

(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

*/

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include "modapuzr.h"  /* Definition des Anwendungsobjectes */
#include "wasim.h"
#include "wsstream.h"

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


ModellApplication *Anwendung;  /* deklaration des Anwendungsobjectes */

int controll_c_flag;

void fpe(int sig)
{
	signal(SIGFPE,(functionpointer)fpe);
	if (Anwendung->unsatzonmodell->fpeflag==0)
	{
		Anwendung->expokorrektur->fpeflag=1;
		Anwendung->verdunstungsmodell->fpeflag=1;
		Anwendung->schneemodell->fpeflag=1;
		Anwendung->interzeptionsmodell->fpeflag=1;
		Anwendung->unsatzonmodell->fpeflag=1;
		Anwendung->irrigationmodell->fpeflag=1;
		Anwendung->gwflow->fpeflag=1;
		Anwendung->abflussmodell->fpeflag=1;
		cout << "floating point error\n";
	}
}

void controll_c(int handlerflag)	/* Signalhandler/Abfangfunktion */
{
	/* ignoriere alle weiteren Signale des gleichen Typs: */

	signal(handlerflag, SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif
	if (Anwendung->controll_c_flag==0)
	{
		Anwendung->controll_c_flag=1;
		cout << "\nCntr-C pressed! please wait...\n\n";
	}
}


void SetSignalAbfangfunktion()
{
	/* catching Ctrl-C and Ctrl-Break */
	signal(SIGINT,(functionpointer)controll_c);
#if 	defined _BCC_
	signal(SIGBREAK,(functionpointer)controll_c);
#endif
	//  printf("Ctrl-C active again.\n");
}

void SetSignalDefaults()
{
	/* setting Ctrl-C interupts to the system-default-handler */
	signal(SIGINT,SIG_DFL);
#if 	defined _BCC_
	signal(SIGBREAK,SIG_DFL);
#endif
	//  printf("Ctrl-C set to system default.\n");
}

std::string toUpper(std::string in)
{
	for (std::string::iterator iter=in.begin();iter<in.end();iter++)
	{
		if (((iter[0])>='a')&&((iter[0])<='z'))
		{
			iter[0]-=32;
		}
	}
	return in;
}


// Usage text for -h/--help. Printed before MPI or OpenMP are started, so that asking
// for help (or calling without arguments) never launches parallel processes.
static void PrintUsage(const char *program)
{
	printf(
		"WaSiM %s\n"
		"\n"
		"usage: %s <control file> [<control file> ...] [options]\n"
		"\n"
		"Runs the model for each control file in turn.\n"
		"\n"
		"options:\n"
		"  -h, --help                 show this help and exit\n"
		"  --version                  show the version and exit\n"
		"  -diagfile:<file>           write diagnostics (XML) to <file> instead of WaSiM_Diag.xml\n"
		"  -diaglevel:<level>         FATAL, ERROR, WARNING, INFO or DEBUG (default)\n"
		"  -mode:routing              only do channel routing, from the given control file\n"
		"  -dumpmode[:<n>]            write a dump file when soil moisture, soil temperature or\n"
		"                             hydraulic head leave their expected range (at most n files)\n"
		"  -closeAfterDump            stop after the first dump\n"
		"  -repairAfterDump           try to repair the values after a dump and continue\n"
		"\n"
		"environment:\n"
		"  OMP_NUM_THREADS            number of OpenMP threads\n"
		"  OMP_SCHEDULE               schedule of the parallel loops (see README)\n"
		"\n"
		"The model description (manual) is docs/wasim_2025_en.pdf in the source repository.\n",
		_wasimuzr_version_english_, program);
}

int main(int argc,char *argv[])
{
	int steuerdateinummer;
	int anfang;
	char *revisionstring,*revisionstring2;

	int mpi_rank=0;
	int mpi_size=1;

	// handle help and version requests before starting MPI
	for (int i = 1; i < argc; i++)
	{
		if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-help") == 0))
		{
			PrintUsage(argv[0]);
			return 0;
		}
		if (strcmp(argv[i], "--version") == 0)
		{
			printf("WaSiM %s\n", _wasimuzr_version_english_);
			return 0;
		}
	}
	if (argc < 2)
	{
		PrintUsage(argv[0]);
		return 2;
	}

#ifdef _MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);



#ifdef _MPIDEBUG
if (mpi_rank == 0)
{
	cout << "option _MPIDEBUG was set during compiling the code. After attaching all ranks processes to the debugger, start the debug run by entering an integer and hitting Enter..." << endl;
	int x;
	scanf("%d",&x);
}
#endif
MPI_Barrier(MPI_COMM_WORLD);

#endif

	wasim::mystream *WSstream;
	int diaglevel;
	int diaglevelparam=0;
	int diagfileparam=0;
	std::string paramstrlvl("-diaglevel:");
	std::string paramstrname("-diagfile:");
	std::string diagfilename;

	int modeparam=0;
	std::string paramstrmode("-mode:");
	std::string modename;


	try
	{

		// default: DEBUG (FATAL=0, ERROR=1, WARNING=2, INFO=3, DEBUG=4, SKIP=5)
		diaglevel=4;
		for (int i=1;i<argc;i++)
		{
			std::string actparam(argv[i]);
			if (paramstrlvl.compare(actparam.substr(0,11))==0) 	{
				diaglevelparam=i;
				std::string sLevel = actparam.substr(11);
				sLevel=toUpper(sLevel);
				if (sLevel.compare("FATAL")==0) {
					diaglevel=0;
				} else if (sLevel.compare("ERROR")==0) {
					diaglevel=1;
				} else if (sLevel.compare("WARNING")==0) {
					diaglevel=2;
				} else if (sLevel.compare("INFO")==0) {
					diaglevel=3;
				} else if (sLevel.compare("DEBUG")==0) {
					diaglevel=4;
				} else {
					diaglevel=4;
				}
			} else if (paramstrname.compare(actparam.substr(0,10))==0) 	{
				diagfileparam=i;
				diagfilename = actparam.substr(10);
			} else if (paramstrmode.compare(actparam.substr(0, 6)) == 0) {
				modeparam = i;
				modename = actparam.substr(6);
			}
		}

		WSstream = new wasim::mystream(diagfilename,diaglevel);

		Anwendung=NULL;
		revisionstring=(char*)calloc(_stringsize_,sizeof(char));
		revisionstring2=(char*)calloc(_stringsize_,sizeof(char));

		sprintf(revisionstring2,"%s",_wasimuzr_version_english_);
		strcpy(revisionstring,revisionstring2);
		while (strstr(revisionstring2,"$"))
		{
			revisionstring2[strstr(revisionstring2,"$")-revisionstring2]='\0';
			sprintf(revisionstring,"%s%s",revisionstring2,strstr(revisionstring,"$")+1);
			strcpy(revisionstring2,revisionstring);
		}

		(*WSstream) << "----------------------------------------------------------------------------"<<"[end-entry]";

#ifdef _MPI
		(*WSstream) << "WaSiM-ETH (MPI) " << "[end-entry]";
#elif _OPENMP
		(*WSstream) << "WaSiM-ETH (OpenMP) " << "[end-entry]";
#else
		(*WSstream) << "WaSiM-ETH " << "[end-entry]";
#endif

#ifdef WIN64
		(*WSstream) << "64-bit Version on Windows "<< "[end-entry]";
#elif WIN32
		(*WSstream) << "32-bit Version on Windows "<< "[end-entry]";
#elif LINUX32
		(*WSstream) << "32-bit Version on Linux "<< "[end-entry]";
#elif LINUX64
		(*WSstream) << "64-bit Version on Linux "<< "[end-entry]";
#elif MACOS64
		(*WSstream) << "64-bit Version on macOS "<< "[end-entry]";
#else
		(*WSstream) << " unspecified version "<< "[end-entry]";
#endif

		(*WSstream) << "Version: " << revisionstring << "[end-entry]";
		(*WSstream) << "Runoff and Water Balance Simulation Model\n"
		"Copyright (C) 1995-2019 Joerg Schulla, 2000-2004 Karsten Jasper\n"
		"2016-2018 Kristian Foerster and Matthias Kopp\n"
		"Any use for commercial purposes needs a suitable licence agreement!\n"
		"All changes in the sources must be reported to the authors!\n"
		"Best wishes for your application!\n\n"
		"Contact: j.schulla@wasim.ch\n"
		"---------------------------------------------------------------------------"<<"[end-entry]";
#ifdef _MPI
		(*WSstream).ssc->flush();
		if (mpi_rank==0) {
			(*WSstream) << "\n this is the MPI-Version of WaSiM, here is instance " << mpi_rank+1 << " of " << mpi_size << "[end-entry]";
	        int th_id, nthreads;
#ifdef _OPENMP
			#pragma omp parallel private(th_id)
            { 
                th_id = omp_get_thread_num();
                printf("This is OpenMP thread %d\n",th_id);
                #pragma omp barrier
                if ( th_id == 0 ) {
                    nthreads = omp_get_num_threads();
                    printf("*** There are %d OpenMP threads total for MPI process %d\n",nthreads,mpi_rank);
                }
            }
#endif
			fflush(NULL);
		}
		else
		{
			printf("\n this is the MPI-Version of WaSiM, here is instance %d of %d\n",mpi_rank+1,mpi_size);
	        int th_id, nthreads;
#ifdef _OPENMP
            #pragma omp parallel private(th_id)
            { 
                th_id = omp_get_thread_num();
                printf("This is OpenMP thread %d\n",th_id);
                #pragma omp barrier
                if ( th_id == 0 ) {
                    nthreads = omp_get_num_threads();
                    printf("*** There are %d OpenMP threads total for MPI process %d\n",nthreads,mpi_rank);
                }
            }
#endif
			fflush(NULL);
		} 
		(*WSstream).ssc->flush();
#endif
		free(revisionstring); revisionstring=NULL;
		free(revisionstring2); revisionstring2=NULL;
		if (argc<2) anfang=0; else anfang=1;

		if (argc<2)
		{
			(*WSstream) << "FATAL: \nusage: wasimvz[o]<32|64>|wasimmpi64 <1-st cntrl-file> [<2-nd cntrl-file> ... [n-th cntrl-file]] [-diagfile:<diagnostics_file> [-diaglevel:<FATAL|ERROR|WARNING|INFO|DEBUG>]] [-dumpmode [-closeAfterDump|-repairAfterDump]]\n" << "[end-entry]";
			(*WSstream) << "INFO: \n-dumpmode: if soil moisture or soil temperatures or hydraulic head (all in Richards approach only) leave the expected range, a dumpfile is written containing the actual internal registers (can be used for debugging) "
				<< "\n-closeAfterDump: WaSiM stops the program execution immediately (usefull to prevent the hard disk from being filled up with dumpfiles)"
				<< "\n-repairAfterDump: let WaSiM try to repair the unexpected results, e.g. a soil temperature that is higher or lower than that of both neighboring soil layers will be set to the maxumum/minimum of the neighboring layers etc. Although the error is repaired it may well happen that the hard disk is flooded with other dump files." << "[end-entry]";
			throw 2;
		}

		for (steuerdateinummer=anfang;steuerdateinummer<argc;steuerdateinummer++)
		{
			std::string sActParam(argv[steuerdateinummer]);
			if ((diaglevelparam!=steuerdateinummer) &&
				(diagfileparam != steuerdateinummer) &&
				(modeparam != steuerdateinummer) &&
				(sActParam.find("-dumpmode") == std::string::npos) &&
				(sActParam.find("-closeAfterDump") == std::string::npos) &&
				(sActParam.find("-repairAfterDump") == std::string::npos))
			{
				Anwendung = new ModellApplication(WSstream,argc,argv,steuerdateinummer);
				
				/* catching Ctrl-C and Ctrl-Break */
				controll_c_flag=0;
				SetSignalAbfangfunktion();
				/* catching floating point errors */
				signal(SIGFPE,(functionpointer)fpe);

				controll_c_flag=Anwendung->run();
				Anwendung->ClearMemory();
				delete Anwendung;
			}
		}
		(*WSstream) << "exitcode " << controll_c_flag << "[end-entry]";
		(*WSstream) << "o.k.\n\nthe end."<<"[end-entry]";
#ifdef _MPI
		MPI_Finalize();
#endif
		delete WSstream;
		return(controll_c_flag);
	}
	catch (int ex)
	{
		(*WSstream) << "FATAL: Exception " << ex << " caught at top level!" << "[end-entry]";
		(*WSstream) << "FATAL: Program aborted with exit code " << ex << "[end-entry]";
#ifdef _MPI
		delete WSstream;
		MPI_Abort(MPI_COMM_WORLD, ex);
#else
		delete WSstream;
		exit(ex);
#endif
	}
}


