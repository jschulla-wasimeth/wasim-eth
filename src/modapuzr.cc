/*   $Revision: 61 $ 
	$Date: 15.09.10 1:07 $  */

/*

(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

Von wasimuzr aus aufgerufenes Hauptprogramm des Modells WaSiM-ETH, Version 3.5

Stand: Juli 1997 einige Routinen, wie Infiltrationsmodell, sind rausgefallen, unsatzon neu hinzugekommen

Stand: 19. September 1997
		(Verduntungs- und Interzeptionsberechnung sind in ihrer Aufrufreihenfolge mit neuem Bodenmodell abgestimmt)

Stand: 6.Januar 1998
		Bodenmodell als FDM nach der Richards-Gleichung mit Interflow, Infiltartionsmodell (ungeschichteter Boden),
		und eigenstaendiger Basisabflussbildung noch ohne Kopplung an ASM

State: 14. January 1998
		simple Irrigation-submodel included
State: February 1998
		Ctrl-C-Handler implemented
State: February/March 1998
		substance-transport (program control basics) implemented

date : 16.6.1998  error in end-date for irrigation module is fixed
date : 9.9.1998   groundwater model implemented

date       : march, 9 1999
version    : hour-entry can be given for start or for end of time step via min_to_int_end

date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.

*/

#include <errno.h>
#include <limits.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <time.h>
#ifdef WIN32
#include <windows.h>
#elif WIN64
#include <windows.h>
#else
#include <stdio.h>
#endif
#include "modapuzr.h"
#include "wasim.h"
#include "wsstream.h"

#ifdef _MPI
#include <mpi.h>
#endif


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#include <cmath>
#include <cstdlib>
#include <vector>

/* Debug aid for porting processes: dump named grids at fixed points of the time loop.
   Only active when WASIM_DUMP_GRIDS (comma-separated grid identifications, as in modconst.h,
   e.g. "precipitation,snow_storage_solid") and WASIM_DUMP_DIR are set; WASIM_DUMP_POINTS
   optionally restricts the points (default: all). Each (point, grid) pair is appended to
   <dir>/<point>__<grid>.f64 as raw float64, row-major over rows 1..nrows and columns 1..ncols
   (stacks: levels+1 values per cell); rows not held by this process are NaN. <dir>/manifest.txt
   gets one line "point grid nrows ncols levels" per file, <dir>/identifications.txt lists all
   grids known at the first point. Nothing is changed in the model. */
static void griddump(const char *point, gridinfotype *gridinfo, int nrows, int ncols)
{
	static int active = -1;
	static std::vector<std::string> names, points;
	static std::string dir, rank;
	static std::vector<std::string> seen;
	if (active == 0) return;
	if (active < 0)
	{
		const char *g = getenv("WASIM_DUMP_GRIDS"), *d = getenv("WASIM_DUMP_DIR"), *pt = getenv("WASIM_DUMP_POINTS");
		active = (g && d && *g && *d) ? 1 : 0;
		if (!active) return;
		dir = d;
		auto split = [](const char *s, std::vector<std::string> &out) {
			std::stringstream ss(s); std::string item;
			while (std::getline(ss, item, ',')) if (!item.empty()) out.push_back(item);
		};
		split(g, names);
		if (pt) split(pt, points);
#ifdef _MPI
		int r = 0, n = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		MPI_Comm_size(MPI_COMM_WORLD, &n);
		if (n > 1) rank = ".rank" + std::to_string(r);
#endif
		FILE *l = fopen((dir + "/identifications" + rank + ".txt").c_str(), "w");
		if (l)
		{
			for (int i = 1; i <= gridmaxanzahl; i++)
				if (gridinfo[i].identification && *gridinfo[i].identification)
					fprintf(l, "%s %d\n", gridinfo[i].identification, gridinfo[i].stackpointer ? gridinfo[i].levels : 0);
			fclose(l);
		}
	}
	if (!points.empty())
	{
		bool wanted = false;
		for (auto &p : points) if (p == point) wanted = true;
		if (!wanted) return;
	}
	for (auto &name : names)
	{
		int k = 0;
		for (int i = 1; i <= gridmaxanzahl; i++)
			if (gridinfo[i].identification && name == gridinfo[i].identification) { k = i; break; }
		if (!k || (!gridinfo[k].pointer && !gridinfo[k].stackpointer)) continue;
		gridinfotype &g = gridinfo[k];
		int nvals = g.stackpointer ? g.levels + 1 : 1;
		std::string safe = name;
		for (auto &c : safe) if (c == '/' || c == ' ') c = '_';
		std::string file = dir + "/" + point + "__" + safe + rank + ".f64";
		bool first = true;
		for (auto &f : seen) if (f == file) first = false;
		if (first)
		{
			seen.push_back(file);
			FILE *m = fopen((dir + "/manifest" + rank + ".txt").c_str(), "a");
			if (m) { fprintf(m, "%s %s %d %d %d\n", point, safe.c_str(), nrows, ncols, nvals - 1); fclose(m); }
		}
		FILE *f = fopen(file.c_str(), first ? "wb" : "ab");
		if (!f) continue;
		std::vector<double> row((size_t)ncols * nvals);
		for (int i = 1; i <= nrows; i++)
		{
			bool have = i >= g.nrl && i <= g.nrh && (g.stackpointer ? g.stackpointer[i] != NULL : g.pointer[i] != NULL);
			for (int j = 1; j <= ncols; j++)
				for (int l = 0; l < nvals; l++)
					row[(size_t)(j - 1) * nvals + l] = !have ? NAN : (g.stackpointer ? g.stackpointer[i][j][l] : g.pointer[i][j]);
			fwrite(row.data(), sizeof(double), row.size(), f);
		}
		fclose(f);
	}
}


ModellApplication::ModellApplication(wasim::mystream *Stream, int argc, char *argv[], int steuerdateinummer) :
	memory_handling(),
	controll_c_flag(0), i(0), taste(0), meteonamen(0), progaus(0), warning_level(0), runoff_measure(0),
	minutes_to_interval_end(0), MinutesToInvEndRouting(0), numsubb(0), listsubb(0), anzahlgroessen(0), breite(0), laenge(0), bezugsmeridian(0),
	programmende(0), steuerdateinamearg(0), steuerdateiname(0), steuerdateiarg(0), steuerdatei(0),
	satz(0), tsatz(0), tag(0), monat(0), stunde(0), jahr(0), minute(0), tga(0), sta(0), mina(0), mta(0), jra(0),
	tge(0), ste(0), mine(0), mte(0), jre(0), LUtable(0), SOILtable(0), irrigationtable(0), zone(0),
	zarray(0), zonmax(0), zonmaxwert(0), ncols(0), nrows(0), mpistartrow(0), mpiendrow(0), rarray(0), regiomax(0), regiomaxwert(0),
	center_of_mass(0), xll(0), yll(0), cellsize(0), nodata(0), FirstLastCol(0), gridinfo(0), gridanzahl(0),
	resultnr(0), timekorr(0), tracer_count(0), tracer_info(0), 
	initial(0), 
	meteo(0), 
	precicorr(0), 
	data_assimilation_model(0),
	superposition(0),
	schneemodell(0),
	interzeptionsmodell(0),
	siltingupmodell(0),
	surfaceroutingmodell(0),
	unsatzonmodell(0), 
	irrigationmodell(0), 
	abflussmodell(0), 
	verdunstungsmodell(0),
	expokorrektur(0), 
	gwflow(0), 
	tracer_ini(0),
    specialoutput(0),
	lastCouplingHour(-1),
	lastCouplingDay(-1),
	lastCouplingMonth(-1),
	lastCouplingYear(-1),
	lastCouplingInterval(-1),
	CouplingIntervalCounter(0),
	CouplingMode('I'),
	WriteASCIIGrids(false),
	WriteDoubleGrids(false),
	WriteNetcdfGrids(false),
	DefaultInputDirectory(0),
	DefaultOutputDirectory(0),
	routingmode(false),
	vHRU(0),
	doSpinUp(0), 
	doTemperaturePreSpinUp(0), 
	repeatCnt(0),
	routingTimestep(0),
	routingTimestepFactor(0)
{
	int j,retvalue;
	try      
	{
		WSstream=Stream;
#ifdef __TIME__
		(*WSstream) << "file " << __FILE__ << ", compiled at " << __DATE__ << ", " << __TIME__ << "[end-entry]";
#endif
		controll_c_flag = 0;
		programmende = false;
		satz               = (char *)calloc(_stringsize_+1,sizeof(char)); tsatz=satz;
		steuerdateiname    = (char *)calloc(_stringsize_+1,sizeof(char));
		steuerdateinamearg = (char *)calloc(_stringsize_+1,sizeof(char));
		int substepcnt = 1;

		/* Einlesen des Namens des Steuerfiles oder Abbruch bei Fehlen des Namens als Programmparameter */
		if (argc<2)
		{
			(*WSstream) << "FATAL: \nusage: wasimvzo<32|64> <1.cntrl-file> [<2.cntrl-file> ... [n.cntrl-file]] [-mode:routing] [-diagfile:<diagnostics_file> [-diaglevel:<FATAL|ERROR|WARNING|INFO|DEBUG>]] [-dumpmode[:<maxdumpfiles>]] [-closeAfterDump|-repairAfterDump]\n"<<"[end-entry]";
			(*WSstream) << "INFO: \n-mode:routing: only channel routing is done after reading the given control file\n-dumpmode: if soil moisture or soil temperatures or hydraulic head (all in Richards approach only) leave the expected range, a dumpfile is written containing the actual internal registers (can be used for debugging) "
				<< "\n-closeAfterDump: WaSiM stops the program execution immediately (usefull to prevent the hard disk from beeing filled up with dumpfiles)"
				<<"\n-repairAfterDump: let WaSiM try to repair the unexpected results, e.g. a soil temperature that is higher or lower than that of both neighboring soil layers will be set to the maxumum/minimum of the neighboring layers etc. Although the error is repaired it may well happen that the hard disk is flooded with other dump files." << "[end-entry]";
			throw 2;
		}
		else strcpy(steuerdateinamearg,argv[steuerdateinummer]);

		std::string paramstrmode("-mode:");
		std::string modename;
		for (int i = 1; i < argc; i++)
		{
			std::string actparam(argv[i]);
			if (paramstrmode.compare(actparam.substr(0, 6)) == 0) {
				modename = actparam.substr(6);
				routingmode = true;
			}
		}

		gridinfo = (gridinfotype *)calloc(gridmaxanzahl+1,sizeof(gridinfotype));
		if (!gridinfo) 
		{
			(*WSstream) << "FATAL: no free memory for creating grids "<<"[end-entry]";
			throw 1;
		}

		/* Initialisieren sicherheitshalber trotz calloc. Die Teilmodelle bekommen gridinfo
	als Adresse uebergeben und belegen ihrerseits Speicher, indem sie suchen,
	welcher index als erster NULL ist
	*/

		for (i=1;i<=gridmaxanzahl;i++)
		{
			gridinfo[i].pointer        = gridinfo[i].minpointer      = gridinfo[i].maxpointer   = NULL;
			gridinfo[i].summenpointer  = NULL;
			gridinfo[i].name           = gridinfo[i].summenname      = gridinfo[i].minname      = gridinfo[i].maxname      = NULL;
			gridinfo[i].writegrid      = gridinfo[i].summenwritegrid = gridinfo[i].minwritegrid = gridinfo[i].maxwritegrid = 0;
			gridinfo[i].stackpointer   = NULL;
			gridinfo[i].levels         = 1;
			gridinfo[i].identification = NULL;
			gridinfo[i].anzahl         = 1;
			gridinfo[i].minmaxflag     = 0;
			// variables for vegetation period dependent balances
			gridinfo[i].vegWritegrid   = 0;
			gridinfo[i].vegSumPointer  = 0;   // pointer to a grid with summed up values for active Vegetation 
			gridinfo[i].vegName        = 0;   // output name for grid with vegetation period balance
			gridinfo[i].vegWriteFlag   = 0;   // is set externally to 1, grids must be written at next call
			gridinfo[i].vegResetFlag   = 0;   // if set externally to 1, all Values in vegSumPointer are reset to 0, vegCnt=0, vegWriteFlag=0
			gridinfo[i].VEGSTART       = NULL;
			gridinfo[i].VEGSTOP        = NULL;
			gridinfo[i].VEGDURAT	     = NULL;
			gridinfo[i].iPer			 = -1; // periodicity of input grid (to be read each time the period is over
			gridinfo[i].iPerAct		 = -1; // actual counter of elapsed periods as given in cPerType (H, D, M or Y)
			gridinfo[i].iSubInt		 = 0; // number of subintervals of the periodtype which must be over until reading a new grid (e.g. 15 for 15 days when cPerType is M
			gridinfo[i].iSubIntAct	 = 0; // actual counter of elapsed sub intervals until the grid must be imported again
			gridinfo[i].iLastSubInt	 = 0; // used to determine if the actual sub interval is less than the old one (when a new month starts, the day will be smaller than the last day of the old month)
			gridinfo[i].iLastPer		 = 0; // same for period counter
			gridinfo[i].cPerType		 = NULL; // type of periodicity: H, D, M, Y
			gridinfo[i].cPerGridName   = NULL; // gridname to read from (will be constructed each time a grid is read from the actual dates)
			gridinfo[i].bPersistency   = true; // set to true by default. If set to false (by persistent=false in the control file) 
			// the original grid will be stored in a temporary grid and after one timestep the 
			// new read in grid will be replaced by the old one again
			gridinfo[i].bReplaceGridNow = false; // is set to true when bPersitency=false and a new grid is read in, so in the next call to the read function it will recognize the flag and replace the grid with the default grid again
			gridinfo[i].TempGrid        = NULL;   // here, the old grid will be stored when the bPersistency-Flag is set to false
			gridinfo[i].iReplaceAfter   = 0;
			gridinfo[i].iReplaceCnt     = 0;
			gridinfo[i].defaultValue    = 0;
			gridinfo[i].outname         = NULL;
			gridinfo[i].statcode        = 0;
			gridinfo[i].statname        = NULL;
			gridinfo[i].zone            = NULL;
			gridinfo[i].stat            = NULL;
			gridinfo[i].zn              = NULL;
			gridinfo[i].zarray          = NULL;
			gridinfo[i].readcode        = 1; // default: read grid content
			gridinfo[i].outnameinternally = 0;
			gridinfo[i].subIntervalCnt = 0;
			gridinfo[i].lastYear = 0;
			gridinfo[i].lastMonth = 0;
			gridinfo[i].lastDay = 0;
			gridinfo[i].lastHour = 0;
			gridinfo[i].lastMinute = 0;
			gridinfo[i].ncl = 0;
			gridinfo[i].nch = 0;
			gridinfo[i].nrl = 0;
			gridinfo[i].nrh = 0;
		}

		maxLastResults = _maxResults_;
		LastResults = (TLastResults*)calloc(_maxResults_+1,sizeof(TLastResults));
		for (i=0;i<=_maxResults_;i++)
		{
			LastResults[i].wert      = NULL;
			LastResults[i].refvalue  = NULL;
			LastResults[i].rel       = NULL;
			LastResults[i].ZoneCodes = NULL;
			LastResults[i].Name      = NULL;
			LastResults[i].NumZones  = 0;
			LastResults[i].j         = 0;
			LastResults[i].m         = 0;
			LastResults[i].s         = 0;
			LastResults[i].t         = 0;
		}

		/* Lesen der allgemeinen Parameter und initialisieren der allgemeinen Datenstrukturen */

		initial = new defaultini(retvalue,WSstream,steuerdateinamearg,LastResults,maxLastResults,2,&gridinfo);
		if (retvalue!=0) throw retvalue;
	
		int nextParamsFunctionState = initial->ParamsFunctionState;
		bool ready=false;
		do
		{
			if (!routingmode || (routingmode == true && (nextParamsFunctionState <= 22 || nextParamsFunctionState > 22)))
			{
				nextParamsFunctionState = initial->Parameters(&ready);
			}
			else
			{
				initial->ParamsFunctionState++;
				nextParamsFunctionState++; // set the next iteration manually
			}
		} 
		while (!ready);

		initial->CreateStandardGridsStatistics();
		initial->GetParameters(steuerdateiname, &iLogLevel,
			&progaus, &warning_level, &runoff_measure,
			&minutes_to_interval_end, &MinutesToInvEndRouting, &timestep,
			&anzahlgroessen, &breite, &laenge, &bezugsmeridian,
			&timekorr, &meteonamen, &numsubb, &listsubb,
			&gridanzahl, &zarray, &zone, &FirstLastCol, &zonmax, &zonmaxwert,
			&center_of_mass, &regiomax, &regiomaxwert, &rarray,
			&ncols, &nrows, &mpistartrow, &mpiendrow, &xll, &yll, &cellsize, &nodata,
			&mina, &sta, &tga, &mta, &jra, &mine, &ste, &tge, &mte, &jre,
			&LUtable, &SOILtable, &irrigationtable, &MultiLayerLUTable,
			&WriteASCIIGrids, &WriteDoubleGrids, &WriteNetcdfGrids, &DefaultInputDirectory, &DefaultOutputDirectory,
			&v_mpiregioncnt, &sendHalo, &recvHalo);
		initial->SetInitialized(1); // now the module is really initialized

		resultnr = 2 + gridanzahl; /* die ersten beiden Grids sind hoehenmodell und Zonengrid
				dann folgen <gridanzahl> Standardgrids (Aspect, Slope ..),
				dann erst die Meteogrids, danach irgendwann die Speichergrids */

				/* lesen der Parameter fuer jede Metegroesse und Initialisieren der Datenstrukturen */
				// new place for making MPI regions
		// now in defualtini before any grid is read in in order to define the offsets into the grid files for each MPI-rank (collective I/O from binary files)
		//		initial->MakeMPIRegions(&v_mpiregioncnt,&sendHalo,&recvHalo);

		if (!routingmode) initial->MakeHRUs(doSpinUp, doTemperaturePreSpinUp, repeatCnt, &vHRU);

		if (!routingmode)
		{
			meteo = new meteointerpol*[anzahlgroessen + 1];

			for (i = 1; i <= anzahlgroessen; i++)
			{
				meteo[i] = NULL;
			}

			for (i = 1; i <= anzahlgroessen; i++)
			{
				j = 1;
				while (gridinfo[j].identification != NULL)
				{
					j++;
					if (j == gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no more free grids in memeory, (max. " << gridmaxanzahl << ")!" << "[end-entry]";
						throw 12;
					}
				}
				resultnr = j;
				/* Erzeugen dyn. Instanzen mit Constructor */
				meteo[i] = new meteointerpol(retvalue, WSstream, LastResults, maxLastResults, meteonamen[i], steuerdateiname, iLogLevel, progaus, minutes_to_interval_end, timestep,
					numsubb, listsubb, gridinfo, resultnr, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
					center_of_mass, regiomax, regiomaxwert, rarray,
					ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
					mina, sta, tga, mta, jra, mine, ste, tge, mte, jre, LUtable, MultiLayerLUTable,
					WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, vHRU);
				if (retvalue != 0) throw retvalue;
				meteo[i]->set_warning_level(warning_level);
				meteo[i]->set_hour_offset(minutes_to_interval_end / 60.0);
				if (strstr(meteonamen[i], _irrigation_amounts_)) meteo[i]->set_warning_level(0);
				meteo[i]->SetInitialized(1); // now the module is really initialized
			}
		}
		else
		{
			anzahlgroessen = 0;
		}

		tracer_ini = new tracer_ini_class(retvalue, WSstream, steuerdateiname, iLogLevel, &tracer_count, &tracer_info,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, vHRU);
		if (retvalue != 0) throw retvalue;
		tracer_ini->SetInitialized(1); // now the module is really initialized
		//    tracer_count=0;

		// should be initialized before the precipitation correction because precicorr will use the output of the regional superposition
		superposition = new RegionalSuperposition(retvalue, WSstream, LastResults, maxLastResults, breite, laenge,
			bezugsmeridian, (char *)_RegionalSuperposition_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			regiomax, regiomaxwert, rarray,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);

		precicorr = new precicorrclass(retvalue, WSstream, meteo, meteonamen, steuerdateiname,
			iLogLevel, anzahlgroessen, regiomax, regiomaxwert, rarray, gridinfo, FirstLastCol,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			LUtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;

		/* Uebergeben der maximalen Hoehe aus defaultini an Niederschlagskorrektur. Das ist noetig,
	wenn P nur oder teilweise hoehenabhaengig berechnet wird, weil dann in den regressdaten nur igu und igo stehen,
	aber nicht, wie hoch die Regressionsgerade gilt. will not be needed on grid based corrections */
		precicorr->put_hmax(initial->get_hmax());
		precicorr->SetInitialized(1); // now the module is really initialized

		expokorrektur = new expokorrklasse(retvalue, WSstream, LastResults, maxLastResults, breite, laenge,
			bezugsmeridian, (char *)_Expositionskorrektur_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;
		expokorrektur->SetInitialized(1); // now the module is really initialized


		verdunstungsmodell = new verdunstungsmodellklasse(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_Verdunstungsmodell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;
		verdunstungsmodell->SetExpoKorrAdress(expokorrektur->GetExpoKorrAdress());


		schneemodell = new schneemodellklasse(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_Schneemodell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;
		schneemodell->SetExpoKorrAdress(expokorrektur->GetExpoKorrAdress());

		/* Damit das Verdunstungsmodell vom Schnee verdunsten kann, wird dessen Adresse geholt.
		ist das Schneemodell abgeschalten, dann wird NULL geliefert, das wird bei der
		Verdunstungsberechnung abgefragt */
		verdunstungsmodell->SetSnowStorAdress(schneemodell->GetSnowStorAdress());
		verdunstungsmodell->setSnowModelMethod(schneemodell->getMethod());
		verdunstungsmodell->SetInitialized(1); // now the module is really initialized


		interzeptionsmodell = new interzeptionsmodellklasse(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_Interzeptionsmodell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;

		interceptionMethod = interzeptionsmodell->getMethod();
		if (interceptionMethod > 10)
		{
			bInterceptionBeforeSnow = true;
			// $kf 2018-06-15
// in this case we want to make use of the rain/snow separation parameters from the snow model --> ALSO SURFTEMP???
			interzeptionsmodell->transzone = schneemodell->transzone;
			interzeptionsmodell->T0R = schneemodell->T0R;
			//			//SURFTEMP AUS SNOWMOD --> Wie geht grid Uebergabe???
			//			interzeptionsmodell->SURFTEMP  = schneemodell->SURFTEMP;
		}
		else
		{
			bInterceptionBeforeSnow = false;
		}
		schneemodell->setInterceptionBeforeSnow(bInterceptionBeforeSnow);

		siltingupmodell = new siltingupclass(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_SiltingUpModell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;
		siltingupmodell->setInterceptionBeforeSnow(bInterceptionBeforeSnow);

		// move this (the complete reading of sub time step count parameter) 
		// into the section where we read from the control file. This fix setting is for testing only
		substepcnt = initial->getSubStepCnt();

		RivLnkCellsToSend = &(initial->RivLnkCellsToSend);
		RivLnkCellsToReceive = &(initial->RivLnkCellsToReceive);
		RiverLinkVector = &(initial->RiverLinkVector);

		unsatzonmodell = new unsatzonmodellklasse(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, _Unsatzonmodell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable, irrigationtable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU, RivLnkCellsToSend, RivLnkCellsToReceive, RiverLinkVector, substepcnt);
		if (retvalue != 0) throw retvalue;
		unsatzonmodell->get_max_pumping(); /* initializing the MAXPUMP-array for the first time (maximum amount of irrigation water from groundwater for the next time-step) */
		unsatzonmodell->setInterceptionBeforeSnow(bInterceptionBeforeSnow);

		bool dumpmode = false;
		bool closeafterdump = false;
		bool repairafterdump = false;
		long maxDumpFiles = 0; // default number of dumpfiles is unlimited
		for (int c = 1; c < argc; c++)
		{
			std::string actArgument(argv[c]);

			if (actArgument.find("-dumpmode") != std::string::npos) {
				dumpmode = true;
				// looking for optional number of dumpfiles
				if (actArgument.find(":") != std::string::npos)
				{
#if __cplusplus > 201100L
					// for newer compilers 
					try
					{
						maxDumpFiles = std::stol(actArgument.substr(actArgument.find(":") + 1));
					}
					catch (...)
					{
						// do nothing here, because just the conversion went wrong, so the ddefault for maxDumpFiles is taken as the initial 0
						(*WSstream) << "FATAL: Error in argument list: invalid number of max dumpfiles given in argument " << actArgument << "[end-entry]";
						throw 1061;
					}
					(*WSstream) << "INFO: max dumpfiles given in argument " << actArgument << " read as " << maxDumpFiles << "(__cplusplus defined as " << __cplusplus << ")" << "[end-entry]";
#else
					// optional for compilers that have no C++11 standard yet
					std::string sMaxDumpfiles = actArgument.substr(actArgument.find(":") + 1);
					if (sMaxDumpfiles.length() > 0)
					{
						char sTemp[1000];
						char *endptr = NULL;
						sprintf(sTemp, "%s", sMaxDumpfiles.c_str());
						errno = 0;
						maxDumpFiles = strtol(sTemp, &endptr, 10);
						if ((errno == ERANGE && (maxDumpFiles == LONG_MAX || maxDumpFiles == LONG_MIN))
							|| (errno != 0 && maxDumpFiles == 0) || (endptr == sTemp))
						{
							(*WSstream) << "FATAL: Error (C-lib errno: " << errno << ") in argument list: invalid number of max dumpfiles given in argument " << actArgument << ", stopped at " << endptr << " LONG_MAX is defined as" << LONG_MAX
								<< "[end-entry]";
							throw 1061;
						}
					}
					(*WSstream) << "INFO: max dumpfiles given in argument " << actArgument << " read as " << maxDumpFiles << "(__cplusplus defined as " << __cplusplus << ")" << "[end-entry]";
#endif
				}
			}
			else if (actArgument.find("-closeAfterDump") != std::string::npos) {
				closeafterdump = true;
			}
			else if (actArgument.find("-repairAfterDump") != std::string::npos) {
				repairafterdump = true;
			}

			//			if (strcmp(argv[c],"-dumpmode")==0)	{
			//				dumpmode=true;
			//			}	else if (strcmp(argv[c], "-closeAfterDump") == 0)	{
			//				closeafterdump=true;
			//			}	else if (strcmp(argv[c], "-repairAfterDump") == 0)	{
			//				repairafterdump = true;
			//			}
		}
		unsatzonmodell->setDumpMode(dumpmode, closeafterdump, repairafterdump, maxDumpFiles);

		// because POND-grid will be known after unsatzonmodel only, the initialize_lakes-function is called here and not from the verdunstungsmodellklasse-constructor itself
		verdunstungsmodell->initialize_lakes(_lake_model_);
		schneemodell->initialize_lakes(_lake_model_);

		// also because the state of heat transfer model and layered snow model are known only after unsatzonmodel was initialized, we set these states now in snow model
		schneemodell->setHeatTransferActive(unsatzonmodell->getHeatTransferActive());
		schneemodell->setLayeredSnowActive(unsatzonmodell->getLayeredSnowActive(), unsatzonmodell->getLayeredSnowProperties(), unsatzonmodell->getTOffsetFreshSnow());
		schneemodell->setDumpMode(dumpmode, closeafterdump, repairafterdump, maxDumpFiles);

		surfaceroutingmodell = new surfaceroutingclass(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_SurfaceRoutingModell_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU, substepcnt);
		if (retvalue != 0) throw retvalue;
		surfaceroutingmodell->setDumpMode(dumpmode, closeafterdump, repairafterdump, maxDumpFiles);

		irrigationmodell = new irrigationclass(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, _irrigation_parameters_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo,
			unsatzonmodell->zarray,
			unsatzonmodell->zone, FirstLastCol,
			unsatzonmodell->zonmax,
			unsatzonmodell->zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable, irrigationtable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU);
		if (retvalue != 0) throw retvalue;

		gwflow = new gwflowclass(retvalue, WSstream, LastResults, maxLastResults, breite,
			laenge, bezugsmeridian, _gwflow_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo,
			unsatzonmodell->zarray,
			unsatzonmodell->zone, FirstLastCol,
			unsatzonmodell->zonmax,
			unsatzonmodell->zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mine, ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable, irrigationtable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU, substepcnt);
		if (retvalue != 0) throw retvalue;

		// the commented submodels got their substepcnt-parameter already when creating them.
		// but verdunstungsmodell is not part of substeps except the real_uzr-method, so we have to set substepcnt manually to be used by this method only
//		unsatzonmodell->setSubStepCnt(substepcnt);
//		surfaceroutingmodell->setSubStepCnt(substepcnt);
		verdunstungsmodell->setSubStepCnt(substepcnt);
		//		gwflow->setSubStepCnt(substepcnt);

				/* if groundwaterflow is calculated by the groundwater flow module, than the baseflow is calculated within
			the groundwater module and not in the unsaturated zone module */
		if (gwflow->do_it != 0) unsatzonmodell->gwflow_active = 1; else unsatzonmodell->gwflow_active = 0;

		/* some initializiations for substance transports (discharge_module: subst. transp. init. is done
	within this module, no subststance transport within evaporation and interpolations */

		irrigationmodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		schneemodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		interzeptionsmodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		interzeptionsmodell->SetInitialized(1); // now the module is really initialized

		siltingupmodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		siltingupmodell->SetInitialized(1); // now the module is really initialized

		gwflow->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		/* unsatzonmodell-tracer are initialized after groundwater module, because in unsatzonmodell the concentrations of
	the groundwater module are required, if the groundwater module is processed */
		unsatzonmodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		unsatzonmodell->SetInitialized(1); // now the module is really initialized

		surfaceroutingmodell->substance_transport_ini(retvalue, tracer_count, tracer_info, LastResults, maxLastResults);
		if (retvalue != 0) throw retvalue;
		surfaceroutingmodell->SetInitialized(1); // now the module is really initialized

		int externalTimeStep = unsatzonmodell->getTimeStep();
		timestep = externalTimeStep;
		routingTimestep = unsatzonmodell->getTimeStep("routing_model");

		initial->timestep = timestep; // important for storing time step dependend standard grids

		// if runoff model runs with another time step we need to get the first time step here explcitely based on the other modules timestep and minutes to interval end etc.
		std::vector<std::pair<int, int>>* timestepvectorbegin = getSubSteps(timestep, routingTimestep, minutes_to_interval_end, MinutesToInvEndRouting, jra, mta, tga, sta, mina);
		int locmina = timestepvectorbegin->begin()->first;
		int locsta = timestepvectorbegin->begin()->second;
		// also the last time step is required
		std::vector<std::pair<int, int>>* timestepvectorend = getSubSteps(timestep, routingTimestep, minutes_to_interval_end, MinutesToInvEndRouting, jre, mte, tge, ste, mine);
		int locmine = timestepvectorend->back().first;
		int locste = timestepvectorend->back().second;

		abflussmodell = new abflussmodellklasse(retvalue,WSstream,LastResults,maxLastResults,
							breite,laenge,bezugsmeridian, (char *)_Abflussmodell_,steuerdateiname,iLogLevel,progaus,
							minutes_to_interval_end, MinutesToInvEndRouting, externalTimeStep,
							numsubb,listsubb,
							gridinfo,
							(routingmode == false) ? unsatzonmodell->zarray:zarray,
							(routingmode == false) ? unsatzonmodell->zone: zone,FirstLastCol,
							(routingmode == false) ? unsatzonmodell->zonmax : zonmax,
							(routingmode == false) ? unsatzonmodell->zonmaxwert : zonmaxwert,
							ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
							locmina,locsta,tga,mta,jra, locmine, locste, tge, mte, jre, timekorr,LUtable,SOILtable,MultiLayerLUTable,
							tracer_count,tracer_info,
							WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
							v_mpiregioncnt,sendHalo,recvHalo,doSpinUp, routingmode, vHRU,
							(routingmode == false) ? unsatzonmodell->gridslesen:0);

		if (retvalue!=0) throw retvalue;
		//    abflussmodell->SetRunoffMeasure(runoff_measure);
		routingTimestepFactor = timestep / routingTimestep;
		// together with minutes_to_interval_end and MinutesToInvEndRouting we are now able to run routing in a loop for smaller timesteps

		if (surfaceroutingmodell->do_it > 0)
		{
			surfaceroutingmodell->SetQBAddress(unsatzonmodell->GetQBAdresse());
			surfaceroutingmodell->SetKelsQDAddress(unsatzonmodell->GetKelsQDAdresse());
			surfaceroutingmodell->SetKelsQIAddress(unsatzonmodell->GetKelsQIAdresse());
			surfaceroutingmodell->ReadSlrContents(unsatzonmodell->GetStorageFileName(), unsatzonmodell->gridslesen);
			unsatzonmodell->SetSurfaceRoutingState(true);
			if (abflussmodell->do_it > 0) 
			abflussmodell->SetQgesAdresse(surfaceroutingmodell->GetQgesAdresse());
		}
		else
		{
			if (abflussmodell->do_it > 0) 
			abflussmodell->SetQgesAdresse(unsatzonmodell->GetQgesAdresse());
		}
		abflussmodell->CountMPICells(bInterceptionBeforeSnow); // do_it is zero for all MPI ranks > 0, so this must be called regardless do_it

		if (abflussmodell->do_it > 0)
		{
			abflussmodell->SetRoutingFluxAdresses();
			abflussmodell->SetZoneSWAdress(irrigationmodell->GetZoneSWAdress());    /* now an information exchange between irrigatiion and routing can start */
			/* now the riverflow_table adress is known (allocated in the unsatzon-substance_transport_ini() */
			gwflow->SetRiverflow_tableAdress(unsatzonmodell->GetRiverflow_tableAdress()); // not really used yet
			abflussmodell->SetRiverflow_tableAdress(unsatzonmodell->GetRiverflow_tableAdress());
			abflussmodell->ReadStorageFile(unsatzonmodell->GetStorageFileName(),unsatzonmodell->gridslesen);
			abflussmodell->SetInitialized(1); // now the module is really initialized
		}
		/* some adresses for the irrigation model were not known until unsatzon module and routing module are initialized */
		irrigationmodell->SetRiverflow_tableAdress(unsatzonmodell->GetRiverflow_tableAdress());
		irrigationmodell->SetQgesAdresse(unsatzonmodell->GetQgesAdresse()); // not yet really used 
		irrigationmodell->GetPumpingConcentrations(); /* needs adresses, that are initialized from within unsatzonmodell */
		irrigationmodell->SetMaxIrrigAdress(abflussmodell->GetMaxIrrigAdress(),abflussmodell->GetMaxIrrigFromReservoirsAddress());
		irrigationmodell->SetRtargetAdress(abflussmodell->GetRtargetAdress());
		irrigationmodell->SetAbfrelAdress(abflussmodell->GetAbfrelAdress());
		irrigationmodell->SetInitialized(1); // now the module is really initialized

		gwflow->SetQgesAdresse(unsatzonmodell->GetQgesAdresse()); // not yet really used
		gwflow->SetInitialized(1); // now the module is really initialized
		
		// Qges is used here for glacerized areas only. For those areas, the unsatzonmodel is bypassed
		// when using silting up model, we have to check, what address must be given to the snowmodel to still work properly
		if (schneemodell->do_it > 0)
		{
			schneemodell->SetQgesAdresse(unsatzonmodell->GetQgesAdresse());
			schneemodell->ReadSlrContents(unsatzonmodell->GetStorageFileName(),unsatzonmodell->gridslesen);
			schneemodell->SetInitialized(1); // now the module is really initialized
		}

		specialoutput = new SpecialOutputClass(retvalue,WSstream,LastResults,maxLastResults,
							breite,laenge,bezugsmeridian, (char *)_special_output_,steuerdateiname,iLogLevel,progaus,
							minutes_to_interval_end,numsubb,listsubb,
							gridinfo,zarray,zone,FirstLastCol,zonmax,zonmaxwert,
							ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
							mine,ste,tge,mte,jre,timekorr,LUtable,SOILtable,MultiLayerLUTable,
							WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
							v_mpiregioncnt,sendHalo,recvHalo,doSpinUp, routingmode, vHRU,substepcnt);
		if (retvalue!=0) throw retvalue;

		if (abflussmodell->do_it != 0)
		{
			abflussmodell->get_reservoir_vector(reservoirs);
		}
		if (abflussmodell->do_it_mpi0 != 0)
		{
			lakeCellsPerMPIRank = abflussmodell->lakeCellsPerMPIRank;
		}
		else
		{
			lakeCellsPerMPIRank = NULL;
		}

		data_assimilation_model = new data_assimilation(retvalue, WSstream, LastResults, maxLastResults,
			breite, laenge, bezugsmeridian, (char *)_data_assimilation_, steuerdateiname, iLogLevel, progaus,
			minutes_to_interval_end, numsubb, listsubb,
			gridinfo, zarray, zone, FirstLastCol, zonmax, zonmaxwert,
			ncols, nrows, mpistartrow, mpiendrow, xll, yll, cellsize, nodata,
			mina,sta,tga,mta,jra,mine,ste, tge, mte, jre, timekorr, LUtable, SOILtable, MultiLayerLUTable,
			WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
			v_mpiregioncnt, sendHalo, recvHalo, doSpinUp, routingmode, vHRU, reservoirs, lakeCellsPerMPIRank, substepcnt);
		if (retvalue != 0) throw retvalue;

		if (unsatzonmodell->do_it > 0)
		{
			(*WSstream) << "initializing re-infiltration in-order-processing..." << "[end-entry]";
			unsigned long numdbl = unsatzonmodell->prepareMPIReinfiltration();
		}

		initial->InitializeExternalCoupling(LastResults,maxLastResults);
		(*WSstream) << "end of initialization reached..."<<"[end-entry]";

		// old place
		//	initial->MakeMPIRegions();

		/* This can be used for testing the output to LastResults   
	ofstream zwf;
	zwf.open("d:\\temp\\zwf.dat",ios::app);
	i=0;
	while ((i<maxLastResults)&&(LastResults[i].wert!=NULL))
	{
	int kk;
	zwf << "------------------\n" << 
			"Name      : " << LastResults[i].Name << "[end-entry]";
	zwf << "NumZones  : " << LastResults[i].NumZones << "[end-entry]";
	zwf << "Codes     :\t";
	for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].ZoneCodes[kk] << "\t" ;
	zwf << "\nAktWert :\t";
	for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].wert[kk] << "\t" ;
	zwf << "\nrefValue:\t";
	for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].refvalue[kk] << "\t" ;
	zwf << "\nRelArea :\t";
	for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].rel[kk] << "\t" ;
	zwf << "[end-entry]";
	i++;
	}
	zwf.close();
*/  
	}
	catch (int ex)
	{
		(*WSstream) << "FATAL: exceptions in ModellApplication::init() caught!" << "[end-entry]";
		throw ex;
	}
	initialized=1;
}

void ModellApplication::RunPreSpinUp()
{
	int temperatureindex=0;
	// looking for temperature interpolation
	char key;
	int mpi_rank=0;
#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#endif

	for (int i=1;i<=anzahlgroessen;i++)
	{
		if (strcmp(meteo[i]->GetInstanzname(),"temperature")==0)
		{
			temperatureindex=i;
			break;
		}
		(*WSstream) << "FATAL: temperature interpolation not found. This is required for pre-SpinUp! (no regions should be used for pre-SpinUp)" << "[end-entry]";
		throw 1035;
	}
	// run the spinup period for temperature interpolation and unsatzon::temperaturePreSpinup once
	while (programmende == false)
	{
		int global_controlc;
#ifdef _MPI // synchronizing external keyboard events with all instances
		MPI_Reduce(&controll_c_flag,&global_controlc,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&global_controlc,1,MPI_INT,0, MPI_COMM_WORLD);
		controll_c_flag=global_controlc;
#endif
//			(*WSstream) <<"global_controlc is "<<global_controlc<<"[end-entry]";
			
		if ((controll_c_flag==1)&&(mpi_rank==0)) // only one instance should communicate with the user
		{
			controll_c_flag=0;
			(*WSstream) <<"SKIP: \n\nCtrl-C: cancel the model? \n"
			"0 = no (continue)\n"
			"1 = cancel with next interval (writing all grids and buffers)\n"
			"2 = cancel now without saving buffers!\n"
			"3 = flush buffers and go on\n"
			"\nchoose one of the above options: "<<"[end-entry]";
			cin >> key;
			char skey[4];
			skey[0]=key;
			skey[1]=0;
			skey[2]=0;
			skey[3]=0;
			controll_c_flag=atoi(&(skey[0]));
			//        while (cin.peek()!='\n') cin.get(key);
		}
		else
		{
			controll_c_flag=0;
		}

#ifdef _MPI // synchronizing external keyboard events with all instances
		MPI_Reduce(&controll_c_flag,&global_controlc,1,MPI_INT,MPI_MAX,0, MPI_COMM_WORLD);
		MPI_Bcast(&global_controlc,1,MPI_INT,0, MPI_COMM_WORLD);
		controll_c_flag=global_controlc;
#endif
			
		if (controll_c_flag>0)
		{
			if ((controll_c_flag<0)||(controll_c_flag>3)) controll_c_flag=0;
			if (controll_c_flag==1)
			{
				for (int i=1;i<=anzahlgroessen;i++) meteo[i]->set_end_date();
				meteo[1]->getdatetime(&minute,&stunde,&tag,&monat,&jahr);
				expokorrektur->set_end(jahr,monat,tag,stunde, minute);
				superposition->set_end(jahr,monat,tag,stunde, minute);
				verdunstungsmodell->set_end(jahr,monat,tag,stunde, minute);
				schneemodell->set_end(jahr,monat,tag,stunde, minute);
				interzeptionsmodell->set_end(jahr,monat,tag,stunde, minute);
				siltingupmodell->set_end(jahr,monat,tag,stunde, minute);
				irrigationmodell->set_end(jahr,monat,tag,stunde, minute);
				unsatzonmodell->set_end(jahr,monat,tag,stunde, minute);
				surfaceroutingmodell->set_end(jahr,monat,tag,stunde, minute);
				gwflow->set_end(jahr,monat,tag,stunde, minute);
				abflussmodell->set_end(jahr,monat,tag,stunde, minute);
				initial->set_end(jahr,monat,tag,stunde, minute);
				specialoutput->set_end(jahr,monat,tag,stunde, minute);
				controll_c_flag=99;
			}
			else
			{
				if (controll_c_flag==2)
				{
					(*WSstream) << "WARNING: program manually killed by pressing Ctrl-C and selecting 2!"<<"[end-entry]";
					delete WSstream;
					exit(98);
				}
				else
				{
					if (controll_c_flag==3)
					{
						for (int i=1;i<=anzahlgroessen;i++) meteo[i]->flushbuffers();
						superposition->flushbuffers();
						verdunstungsmodell->flushbuffers();
						schneemodell->flushbuffers();
						interzeptionsmodell->flushbuffers();
						siltingupmodell->flushbuffers();
						irrigationmodell->flushbuffers();
						unsatzonmodell->flushbuffers();
						surfaceroutingmodell->flushbuffers();
						gwflow->flushbuffers();
						abflussmodell->flushbuffers();
				  		initial->FlushBuffersGridStatistics();
						specialoutput->flushbuffers();
						controll_c_flag=0;
					}
				}
			}
			/* Es soll also weitergerechnet werden.
	Erlaube nun wieder Signale, und fange Sie erneut ab: */
			SetSignalAbfangfunktion();
		}
		meteo[temperatureindex]->getdatetime(&minute,&stunde,&tag,&monat,&jahr);     /* Die Modell-Zeit wird vom ersten Meteoobject geholt */
		if ((meteo[temperatureindex]->getmethode()==10)||(meteo[temperatureindex]->getmethode()==11))
		{
			meteo[temperatureindex]->CalculateRegressionParameters(0); // only one set of regressionparameters --> one region only in direct regression
		}
		if (stunde+(minute+minutes_to_interval_end)/60.0==24) (*WSstream) << wasim::setw(20) << "temperature pre-SpinUp "<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde<<":"<<minute<<" h "<<"[end-entry]";
		meteo[temperatureindex]->run(minute,stunde,tag,monat,jahr);
		meteo[temperatureindex]->WriteResultGrid();
		unsatzonmodell->temperaturePreSpinUp(minute, stunde,tag,monat,jahr);
		meteo[temperatureindex]->nextrecord();
		if (meteo[temperatureindex]->isend()) programmende = true;

	}
	(*WSstream) << "INFO: pre-SpinUp reached end of defined time range. Soil-temperatures are going to be set now." << "[end-entry]";
	unsatzonmodell->initializeSoilTemperatures();
	programmende=false; // this was only pre-SpinUp
	// reset meteorologic input and look for the starting date again
	meteo[temperatureindex]->resetInputFiles();
	meteo[temperatureindex]->nextrecord();

}

int ModellApplication::run()
{
	char key;
	int actSpinUpRepeatCnt=1;
	//  int i;
	try
	{
		int mpi_rank=0;
#ifdef _MPI
		MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#endif
		if (!routingmode)
		{

			if (doSpinUp > 0 && doTemperaturePreSpinUp > 0)
			{
				RunPreSpinUp();
			}

			if (doSpinUp > 0)
			{
				if (actSpinUpRepeatCnt > repeatCnt)
				{
					programmende = true;
				}
			}
		}

		while (programmende == false)
		{
			int global_controlc;
#ifdef _MPI // synchronizing external keyboard events with all instances
			MPI_Reduce(&controll_c_flag,&global_controlc,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&global_controlc,1,MPI_INT,0, MPI_COMM_WORLD);
			controll_c_flag=global_controlc;
#endif
//			(*WSstream) <<"global_controlc is "<<global_controlc<<"[end-entry]";
			
			if ((controll_c_flag==1)&&(mpi_rank==0)) // only one instance should communicate with the user
			{
				controll_c_flag=0;
				(*WSstream) <<"SKIP: \n\nCtrl-C: cancel the model? \n"
				"0 = no (continue)\n"
				"1 = cancel with next interval (writing all grids and buffers)\n"
				"2 = cancel now without saving buffers!\n"
				"3 = flush buffers and go on\n"
				"\nchoose one of the above options: "<<"[end-entry]";
				cin >> key;
				char skey[4];
				skey[0]=key;
				skey[1]=0;
				skey[2]=0;
				skey[3]=0;
				controll_c_flag=atoi(&(skey[0]));
				//        while (cin.peek()!='\n') cin.get(key);
			}
			else
			{
				controll_c_flag=0;
			}

#ifdef _MPI // synchronizing external keyboard events with all instances
//			MPI_Barrier(MPI_COMM_WPRLD);
			MPI_Reduce(&controll_c_flag,&global_controlc,1,MPI_INT,MPI_MAX,0, MPI_COMM_WORLD);
			MPI_Bcast(&global_controlc,1,MPI_INT,0, MPI_COMM_WORLD);
			controll_c_flag=global_controlc;
#endif

			/* fuer jeden Zeitschritt werden erst einmal die noetigen Eingangsdaten
			auf die in der Steuerdatei mit Namen benannten Grids interpoliert.
			Anzahl der Eingangsdaten, deren Anfangs- und Endzeitpunkt sowie die
			speziellen Parameter wurden bei der Initialisierung der meteorologischen
			Interpolations-Instanzen aus der Steuerdatei gelesen */

			if (routingmode)
			{
				abflussmodell->getdatetime(&minute, &stunde, &tag, &monat, &jahr);
			}
			else if (anzahlgroessen >= 1) {
				meteo[1]->getdatetime(&minute, &stunde, &tag, &monat, &jahr);     /* Die Modell-Zeit wird vom ersten Meteoobject geholt */
			}
			else {
				minute = stunde = tag = monat = jahr = -9999;
			}

			// now wait for the external gw table grid if the model runs as coupled model
			if (initial->IsCoupled() && !routingmode)
			{
				CouplingMode = initial->GetCouplingMode();
				if ((CouplingMode == 'I') ||
					((CouplingMode == 'H') && (lastCouplingHour != stunde)) ||
					((CouplingMode == 'D') && (lastCouplingDay != tag)) ||
					((CouplingMode == 'M') && (lastCouplingMonth != monat)) ||
					((CouplingMode == 'Y') && (lastCouplingYear != jahr)))
				{
					if (mpi_rank == 0)
					{
						WaitForFile(initial->GetSemaphoreFileNameRead(), initial->GetSemaphoreWaitInterval());
					}
#ifdef _MPI
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Reduce(&controll_c_flag, &global_controlc, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
					MPI_Bcast(&global_controlc, 1, MPI_INT, 0, MPI_COMM_WORLD);
					controll_c_flag = global_controlc;
#endif
					lastCouplingHour = stunde;
					lastCouplingDay = tag;
					lastCouplingMonth = monat;
					lastCouplingYear = jahr;
					initial->InitializeFromExternalGrids(jahr, monat, tag, stunde);
					if (initial->ExistsExternalGWTable())
					{
						unsatzonmodell->GetExternalGWTable(initial->GetExternalTimeStep());
					}
				}
			}

			
			if (controll_c_flag>0)
			{
				if ((controll_c_flag<0)||(controll_c_flag>3)) controll_c_flag=0;
				if (controll_c_flag==1)
				{
					for (i=1;i<=anzahlgroessen;i++) meteo[i]->set_end_date();
					meteo[1]->getdatetime(&minute,&stunde,&tag,&monat,&jahr);
					expokorrektur->set_end(jahr,monat,tag,stunde,minute);
					superposition->set_end(jahr,monat,tag,stunde, minute);
					verdunstungsmodell->set_end(jahr,monat,tag,stunde, minute);
					schneemodell->set_end(jahr,monat,tag,stunde, minute);
					interzeptionsmodell->set_end(jahr,monat,tag,stunde, minute);
					siltingupmodell->set_end(jahr,monat,tag,stunde, minute);
					irrigationmodell->set_end(jahr,monat,tag,stunde, minute);
					unsatzonmodell->set_end(jahr,monat,tag,stunde, minute);
					surfaceroutingmodell->set_end(jahr,monat,tag,stunde, minute);
					gwflow->set_end(jahr,monat,tag,stunde, minute);
					abflussmodell->set_end(jahr,monat,tag,stunde, minute);
					initial->set_end(jahr,monat,tag,stunde, minute);
					specialoutput->set_end(jahr,monat,tag,stunde, minute);
					data_assimilation_model->set_end(jahr, monat, tag, stunde, minute);
					controll_c_flag=99;
				}
				else
				{
					if (controll_c_flag==2)
					{
						(*WSstream) << "WARNING: program manually killed by pressing Ctrl-C and selecting 2!"<<"[end-entry]";
						delete WSstream;
						exit(98);
					}
					else
					{
						if (controll_c_flag==3)
						{
							for (i=1;i<=anzahlgroessen;i++) meteo[i]->flushbuffers();
							superposition->flushbuffers();
							verdunstungsmodell->flushbuffers();
							schneemodell->flushbuffers();
							interzeptionsmodell->flushbuffers();
							siltingupmodell->flushbuffers();
							irrigationmodell->flushbuffers();
							unsatzonmodell->flushbuffers();
							surfaceroutingmodell->flushbuffers();
							gwflow->flushbuffers();
							abflussmodell->flushbuffers();
				  		    initial->FlushBuffersGridStatistics();
							specialoutput->flushbuffers();
							data_assimilation_model->flushbuffers();
							controll_c_flag=0;
						}
					}
				}
				/* Es soll also weitergerechnet werden.
		Erlaube nun wieder Signale, und fange Sie erneut ab: */
				SetSignalAbfangfunktion();
			}

			for (int i=1;i<=anzahlgroessen;i++)
			{
				if ((meteo[i]->getmethode()==10)||(meteo[i]->getmethode()==11))
				meteo[i]->CalculateRegressionParameters(0); // only one set of regressionparameters --> one region only in direct regression
			}

			griddump("start", gridinfo, nrows, ncols);
			if (superposition->do_it==0)  precicorr->run();

			// if this pragma is uncommented, then the parallel pragmas in meteoint should be commented if the compiler 
			// does not support nested omp-directives -> it's faster to user omp in the meteoint-for routine then here 
			// (estimated by experience)
			//#pragma omp parallel for schedule(dynamic,1) private(i)

			//#pragma omp parallel for schedule(runtime) 
			// attention: tests with this for-loop-parallelized version showed that something is still wrong with it. 
			// is is NOT save to interpolate in parallel here. Rather do it within meteointerpol
			for (int i=1;i<=anzahlgroessen;i++)
			{
				meteo[i]->run(minute, stunde,tag,monat,jahr);
			}

			superposition->run(minute, stunde,tag,monat,jahr);
			if (superposition->do_it || precicorr->doGriddedCorrection)
			{
				precicorr->runOnSuperposedRegions();
			}

			// because temperatures are changed in expokorrektur and in 
			// verdunstungsmodell (+changes for day and -changes for night!) 
			// but the phenology needs original (or real average daily) temperatures, we 
			// calculate the dynamic phenology now - so at the 24th hour of a day 
			// it may happen that the parameters will suddenly change. 
			// But for daily input values it's ok
			int Interval;
			if (!routingmode)
			{
				if (verdunstungsmodell->doit())
				{
					Interval = verdunstungsmodell->GetInV();
				}
				else
				{
					Interval = 1440; // default 1 day
				}
				// reading new standard grids, if required
				// note: this must be done before interpolating variable grids (here after superposing meteorologic values
				// in order to avoid situations where a read in grid is already in use but 
				// has missing or undefined values for cells which will be set by interpolation hereafter
				initial->readPeriodicStandardGrids(minute, stunde, tag, monat, jahr, Interval); // for all standard grids with a cPerType of H, D, M or Y
				initial->int_tvar_params(minute, stunde, tag, monat, jahr, Interval); // interpolate_zeitvariable_parameter
			}

			data_assimilation_model->run(minute, stunde, tag, monat, jahr);

			// the results written should be the changed ones (e.g. changed radiation or temperature by aspect dependent correction)
			for (int i=1;i<=anzahlgroessen;i++)
			{
				meteo[i]->WriteResultGrid();
			}

			superposition->WriteResultGrids();
			griddump("after_meteo", gridinfo, nrows, ncols);

//			// reading new standard grids, if required
//			initial->readPeriodicStandardGrids(minute,stunde,tag,monat,jahr,Interval); // for all standard grids with a cPerType of H, D, M or Y

			/* Errechnen der expositions-, neigungs- und
		abschattungsbedingten Korrekturfaktoren fuer die direkte Sonnenstrahlung */

			expokorrektur->run(minute, stunde,tag,monat,jahr);
			griddump("after_expokorr", gridinfo, nrows, ncols);

			/* Errechnen der potentiellen und realen Verdunstung (fuer Boden/Pflanzen)
		als Eingabegroessen fuer die weiteren Teilmodelle, die reale Verdunstung
		wird mit den Meteogroessen des aktuellen Intervalles und mit den Bodendaten
		des vorhergehenden Intervalles berechnet!
	*/
			verdunstungsmodell->potentiell(minute, stunde,tag,monat,jahr);
			griddump("after_etp", gridinfo, nrows, ncols);

			/* applying irrigation water according to the irrigation table or an optimizing instruction */
			irrigationmodell->run(minute, stunde,tag,monat,jahr);
			griddump("after_irrigation", gridinfo, nrows, ncols);

//			/* Schneeverdunstung von ETP abziehen und GMB aktualisieren*/
//			verdunstungsmodell->snow(minute, stunde,tag,monat,jahr); // this component could be switched off when usig energy balance approach for snow melt. Or at least the approach should be worked over since usually the available energy cannot be used for pre evaporation but the melting energy must be satisfied first

			if (bInterceptionBeforeSnow)
			{
				/* Bilanzieren des Interzeptionsspeichers einfangen, Abtropfen, Verdunsten */
				interzeptionsmodell->run(minute, stunde, tag, monat, jahr);
				griddump("after_interception", gridinfo, nrows, ncols);
				/* Akkumulation, Schmelze, Rueckgefrieren, Ausfluss, evtl. Verdunstung */
				schneemodell->run1(minute, stunde, tag, monat, jahr);
				griddump("after_snow1", gridinfo, nrows, ncols);
			}
			else
			{
				/* Akkumulation, Schmelze, Rueckgefrieren, Ausfluss, evtl. Verdunstung */
				schneemodell->run1(minute, stunde, tag, monat, jahr);
				griddump("after_snow1", gridinfo, nrows, ncols);
				/* Bilanzieren des Interzeptionsspeichers einfangen, Abtropfen, Verdunsten */
				interzeptionsmodell->run(minute, stunde, tag, monat, jahr);
				griddump("after_interception", gridinfo, nrows, ncols);
			}
			/* da die reale Verdunstung in verschiedenen Teilmodellen noch veraendert werden kann,
				wird sie erst als letzte Groesse endgueltig berechnet */

			// Schneeverdunstung von ETP abziehen und GMB aktualisieren
			// hierher von oben verschoben, da erst im Schneemodell EVRS berechnet wird, die hier dann statistisch erfasst wird
			verdunstungsmodell->snow(minute, stunde, tag, monat, jahr); // this component could be switched off when usig energy balance approach for snow melt. Or at least the approach should be worked over since usually the available energy cannot be used for pre evaporation but the melting energy must be satisfied first
			griddump("after_snowevap", gridinfo, nrows, ncols);

			siltingupmodell->run(minute, stunde,tag,monat,jahr);
			griddump("after_siltingup", gridinfo, nrows, ncols);
			
			// dummy frame for time step splitting 
			// since groundwater module, unsaturated zone and surface routing affect each other
			// these submodels (as well as accumulation of real evapotranspiration) should be part of 
			// a greater loop. All we have to do is
			// - defining a common time step ("InV") for each of the affetced submodels
			// - initializing statistics only for the first sub step
			// - accumulating statistics for each sub time step as usual
			// - calculate statistics in the last sub time step only
			// - manging other connections to other submodels e.g. riverflowtable from routing module
			// - make sure the tracer module still works
			unsatzonmodell->pre_run(minute, stunde,tag,monat,jahr,0);
			surfaceroutingmodell->pre_run(minute, stunde,tag,monat,jahr,0);
			verdunstungsmodell->pre_real_uzr(minute, stunde,tag,monat,jahr,0);
			gwflow->pre_run(minute, stunde,tag,monat,jahr,0);
			
			int substepcnt=unsatzonmodell->getSubStepCnt();
			for (int substep=1;substep<=substepcnt;substep++)
			{
				unsatzonmodell->run(minute, stunde,tag,monat,jahr,substep);
				griddump("after_unsatzon", gridinfo, nrows, ncols);
				// surface routing is not yet MPI or OpenMP aware, so all grids used by Surface routing must be synchronized now or at least when starting surface routing
				surfaceroutingmodell->run(minute, stunde,tag,monat,jahr,substep);
				griddump("after_surfrout", gridinfo, nrows, ncols);
				verdunstungsmodell->real_uzr(minute, stunde,tag,monat,jahr,substep);
				griddump("after_etreal", gridinfo, nrows, ncols);
				// gwflow is not yet MPI aware, so all grids used by GWflow must be syncronised now or at least when starting wflow module
				gwflow->run(minute, stunde,tag,monat,jahr,substep);
				griddump("after_gwflow", gridinfo, nrows, ncols);
				specialoutput->run(minute, stunde,tag,monat,jahr,true,substep);
			}
//			(*WSstream) << "back in main loop"<< "[end-entry]";

			schneemodell->run2(minute, stunde, tag, monat, jahr);
			griddump("after_snow2", gridinfo, nrows, ncols);

			unsatzonmodell->get_max_pumping();

			unsatzonmodell->post_run(minute, stunde,tag,monat,jahr,0);
			surfaceroutingmodell->post_run(minute, stunde,tag,monat,jahr,0);
			verdunstungsmodell->post_real_uzr(minute, stunde,tag,monat,jahr,0);
			gwflow->post_run(minute, stunde,tag,monat,jahr,0);
			griddump("after_post_run", gridinfo, nrows, ncols);



#ifdef _MPI
//		MPI_Barrier(MPI_COMM_WPRLD);
#endif

			// if the timestep for routing and unsatzon is different (i.e. routing has a shorter timestep), we need to cann the routing several times now
			// routingTimestep = time step in minutef routing
			// timestep = time step in unsatzonmodell
			// routingTimestepFactor = externalTimeStep / routingTimestep;
			// minutes_to_interval_end = 0 or time step length in minutes for normal modules
			// MinutesToInvEndRouting   = 0 or time step length in minutes for routing

			std::vector<std::pair<int, int>>* timestepvector = getSubSteps(timestep, routingTimestep, minutes_to_interval_end, MinutesToInvEndRouting, jahr, monat, tag, stunde, minute);
//			std::vector<std::pair<int, int>>* timestepvector = getSubSteps(1440, 60, 1440, 0, 1984, 12, 31, 0, 0);
			// sample input: 1984 01 01 12 (00), i.e. 12 is the hour value of the middle of the day 
			// sample output: 
			// 1984 01 01 00 30
			// 1984 01 01 01 00
			// ...
			// 1984 01 01 23 30
			// 1984 01 01 24 00
			int actsubtimestep = 0;
			int numsubtimesteps = timestepvector->size();
			for (std::vector<pair<int, int>>::iterator iter = timestepvector->begin(); iter != timestepvector->end(); iter++)
			{
				int minute = iter->first;
				int stunde = iter->second;
				abflussmodell->run(minute, stunde, tag, monat, jahr, actsubtimestep, numsubtimesteps);
				griddump("after_routing", gridinfo, nrows, ncols);
				actsubtimestep++;
			}
			delete timestepvector;

			if (!routingmode)
			{
				initial->write_Standard_Grids(minute, stunde, tag, monat, jahr);
				initial->write_Variable_Grids(minute, stunde, tag, monat, jahr);

				// now all standard grids which have a statistics output defined are processed, i.e. 
				// the statistics registers are updated
				initial->InitGridsStatistics();
				initial->AddValuesGridStatistics();
				initial->CalcAndDoneGridStatistics(minute, stunde, tag, monat, jahr);
			}
			specialoutput->run(minute, stunde,tag,monat,jahr,false,0);

			if (!routingmode)
			{
				for (int i = 1; i <= anzahlgroessen; i++)
				{
					meteo[i]->nextrecord();
				}

				// Das Programmende ist erreicht, wenn in allen Eingangsgroessen
				// das in der Steuerdatei angegebene Enddatum erreicht ist 
				programmende = true;
				//#pragma omp parallel for schedule(runtime) 
				for (int i = 1; i <= anzahlgroessen; i++)
				{
					if (!meteo[i]->isend())
					{
						programmende = false;
					}
				}
			}
			else
			{
				programmende = abflussmodell->routingmode_isend();
				abflussmodell->routingmode_nextrecord();
			}

			if (initial->IsCoupled() && !routingmode)
			{
				if (anzahlgroessen>=1) meteo[1]->getdatetime(&minute,&stunde,&tag,&monat,&jahr); /* Die Modell-Zeit wird vom ersten Meteoobject geholt */
				CouplingMode = initial->GetCouplingMode();
				// this is already the time of the next timestep, so if an hour, a day, month or year has changed and the coupling mode is set
				// to the corresponding value, the grids are printed out, otherwise the values are accumulated according to the settings
				if ( (CouplingMode=='I')|| (programmende==true) ||
						((CouplingMode=='H')&&(lastCouplingHour != stunde  )) ||
						((CouplingMode=='D')&&(lastCouplingDay  != tag     )) ||
						((CouplingMode=='M')&&(lastCouplingMonth!= monat   )) ||
						((CouplingMode=='Y')&&(lastCouplingYear != jahr    ))    )
				{
					initial->WriteExternalResults(stunde,tag,monat,jahr,
					(lastCouplingHour != stunde  ) || (programmende == true),
					(lastCouplingDay  != tag     ) || (programmende == true),
					(lastCouplingMonth!= monat   ) || (programmende == true),
					(lastCouplingYear != jahr    ) || (programmende == true));
					if (mpi_rank == 0)
					{
						if (!programmende)
						{
							WriteSemaphore(initial->GetSemaphoreFileNameWrite(), initial->GetSemaphoreContent());
						}
						else
						{
							WriteSemaphore(initial->GetSemaphoreFileNameWrite(), (char *)"Command=Finish;");
							std::string extendedname(initial->GetSemaphoreFileNameWrite());
							extendedname.append("_Done");
							WriteSemaphore(extendedname.c_str(), (char *)"Command=Finish;");
						}
					}
#ifdef _MPI
					MPI_Barrier(MPI_COMM_WORLD);
#endif
				}
				else
				{
					initial->AccumulateExternalEntities();
				}
			}

			if ((doSpinUp==1) && (programmende==true) && !routingmode)
			{
				// check if the program should really close or if the spin uip period must be repeated
				if (actSpinUpRepeatCnt<repeatCnt)
				{
					(*WSstream) << "INFO: SpinUp reached end of defined time range for run #"<<actSpinUpRepeatCnt<<", resetting inputs and starting again now (repetition "<<actSpinUpRepeatCnt+1<<" of "<<repeatCnt<<")" << "[end-entry]";
					programmende=false;
					// reset meteorologic input and look for the starting date again
					for (int i=1;i<=anzahlgroessen;i++) 
					{
						meteo[i]->resetInputFiles();
						meteo[i]->nextrecord();
					}
					actSpinUpRepeatCnt++;
				}
			}


			/*   this can be used for testing output in LastResults
	ofstream zwf;
	zwf.open("d:\\temp\\zwf.dat",ios::app);
	i=0;
	while ((i<maxLastResults)&&(LastResults[i].wert!=NULL))
	{
		int kk;
		zwf << "------------------\n" << 
				"Name      : " << LastResults[i].Name << "[end-entry]";
		zwf << "NumZones  : " << LastResults[i].NumZones << "[end-entry]";
		zwf << "Codes     :\t";
		for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].ZoneCodes[kk] << "\t" ;
		zwf << "\nAktWert :\t";
		for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].wert[kk] << "\t" ;
		zwf << "\nrefValue:\t";
		for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].refvalue[kk] << "\t" ;
		zwf << "\nRelArea :\t";
		for (kk=1;kk<=LastResults[i].NumZones;kk++)
		zwf << LastResults[i].rel[kk] << "\t" ;
		zwf << "[end-entry]";
		i++;
	}
	zwf.close();
*/ 
		}
		/* ------------------ Ende der Modellroutine -------------------- */
	}
	catch (int ex)
	{
		(*WSstream) << "FATAL: Exception in Modellapplication::run() caught!" << "[end-entry]";
		throw ex;
	}
	return(controll_c_flag);
}


/* Destructor des Modells */

void ModellApplication::ClearMemory()
{
	// int i;

	try
	{
		if (WSstream) (*WSstream) << "reached end of simulation sucessfully"<<"[end-entry]";

		/* flushing all buffers in order to avoid data loss, if in one of the following ClearMemory-functions
	a memory access violates security and the program is halted from the system. This may be due to errors
	of the optimiszing compiler, not due to programmers errors! */
		if (initialized)
		{
			for (int i=1;i<=anzahlgroessen;i++) 
			{
				if (meteo[i]) meteo[i]->flushbuffers();
			}
			if (superposition)	   superposition->flushbuffers();
			if (verdunstungsmodell)  verdunstungsmodell->flushbuffers();
			if (schneemodell)        schneemodell->flushbuffers();
			if (interzeptionsmodell) interzeptionsmodell->flushbuffers();
			if (irrigationmodell)    irrigationmodell->flushbuffers();
			if (siltingupmodell)     siltingupmodell->flushbuffers();
			if (unsatzonmodell)      unsatzonmodell->flushbuffers();
			if (surfaceroutingmodell) surfaceroutingmodell->flushbuffers();
			if (gwflow)              gwflow->flushbuffers();
			if (abflussmodell)       abflussmodell->flushbuffers();
			if (specialoutput)		 specialoutput->flushbuffers();
			if (data_assimilation_model) data_assimilation_model->flushbuffers();
			if (initial)			 initial->FlushBuffersGridStatistics();
		}
		else
		{
			for (int i=1;i<=anzahlgroessen;i++) 
			if (meteo[i]) meteo[i]->SetInitialized(0);
			if (superposition)	   superposition->SetInitialized(0);
			if (verdunstungsmodell)  verdunstungsmodell->SetInitialized(0);
			if (schneemodell)        schneemodell->SetInitialized(0);
			if (interzeptionsmodell) interzeptionsmodell->SetInitialized(0);
			if (irrigationmodell)    irrigationmodell->SetInitialized(0);
			if (surfaceroutingmodell) surfaceroutingmodell->SetInitialized(0);
			if (unsatzonmodell)      unsatzonmodell->SetInitialized(0);
			if (siltingupmodell)     siltingupmodell->SetInitialized(0);
			if (gwflow)              gwflow->SetInitialized(0);
			if (abflussmodell)       abflussmodell->SetInitialized(0);
			if (expokorrektur)       expokorrektur->SetInitialized(0);
			if (precicorr)           precicorr->SetInitialized(0);
			if (tracer_ini)          tracer_ini->SetInitialized(0);
			if (specialoutput)       specialoutput->SetInitialized(0);
			if (data_assimilation_model) data_assimilation_model->SetInitialized(0);
			if (initial)             initial->SetInitialized(0);
		}

		for (int i=anzahlgroessen;i>=1;i--)
		{
			if (meteo[i]) meteo[i]->ClearMemory();
			if(meteonamen[i]) free(meteonamen[i]); 
			meteonamen[i]=NULL;
		}
		if (meteonamen) free(meteonamen); 
		meteonamen=NULL;

		if (schneemodell)         schneemodell->substance_transport_delete();
		if (interzeptionsmodell)  interzeptionsmodell->substance_transport_delete();
		if (siltingupmodell)      siltingupmodell->substance_transport_delete();
		// next line: is done in clearmemory now because QDd must be available for writing	SLR contents 
		// if (surfaceroutingmodell) surfaceroutingmodell->substance_transport_delete();
		if (irrigationmodell)     irrigationmodell->substance_transport_delete();
		if (gwflow)               gwflow->substance_transport_delete();
		if (superposition)        superposition->substance_transport_delete();
		if (data_assimilation_model) data_assimilation_model->substance_transport_delete();

		bool appendmode=false;
		if (abflussmodell)        
		{
			if (abflussmodell->do_it) appendmode=true;
			abflussmodell->ClearMemory();
		}
		if (unsatzonmodell)
		{
			bool unsatzon_appendmode = false;
			if (unsatzonmodell->do_it) unsatzon_appendmode = true;
			unsatzonmodell->ClearMemory(appendmode);
			if (unsatzon_appendmode)
			{
				appendmode = true;
			}
		}

		if (siltingupmodell)      siltingupmodell->ClearMemory();
		if (surfaceroutingmodell) surfaceroutingmodell->ClearMemory(appendmode);
		if (interzeptionsmodell)  interzeptionsmodell->ClearMemory();
		if (irrigationmodell)     irrigationmodell->ClearMemory();
		if (schneemodell)         schneemodell->ClearMemory();
		if (verdunstungsmodell)   verdunstungsmodell->ClearMemory();
		if (expokorrektur)        expokorrektur->ClearMemory();
		if (superposition)		  superposition->ClearMemory();
		if (precicorr)            precicorr->ClearMemory();
		if (gwflow)               gwflow->ClearMemory();
		if (tracer_ini)           tracer_ini->ClearMemory();
		if (specialoutput)        specialoutput->ClearMemory();
		if (data_assimilation_model) data_assimilation_model->ClearMemory();
		if (initial)              initial->ClearMemory(); 

		if (gridinfo)           free(gridinfo);           gridinfo=NULL;
		if (steuerdateinamearg) free(steuerdateinamearg); steuerdateinamearg=NULL;
		if (steuerdateiname)    free(steuerdateiname);    steuerdateiname=NULL;
		satz=tsatz; 
		if (satz)               free(satz);               
		satz=NULL; 
		tsatz=NULL;

		if (numsubb>0) 
		{
			if (listsubb) 
			free(listsubb);
		}
		numsubb=0;
		listsubb=NULL; 

		if (WSstream) (*WSstream) << "o.k.\n" << "[end-entry]";

		if (WSstream) (*WSstream) << "deleting modules..."<<"[end-entry]";
		for (int i=anzahlgroessen;i>=1;i--)
		{
			if (meteo[i]) 
			{
				delete meteo[i];
				meteo[i]=NULL;
			}
		}
		if (meteo) delete[] meteo;
		meteo=NULL;
		anzahlgroessen=0;

		if (tracer_ini)          delete tracer_ini;          tracer_ini          = NULL;
		if (abflussmodell)       delete abflussmodell;       abflussmodell       = NULL;
		if (gwflow)              delete gwflow;              gwflow              = NULL;
		if (irrigationmodell)    delete irrigationmodell;    irrigationmodell    = NULL;
		if (siltingupmodell)     delete siltingupmodell;     siltingupmodell     = NULL;
		if (surfaceroutingmodell)delete surfaceroutingmodell; surfaceroutingmodell= NULL;
		if (unsatzonmodell)      delete unsatzonmodell;      unsatzonmodell      = NULL;
		if (interzeptionsmodell) delete interzeptionsmodell; interzeptionsmodell = NULL;
		if (schneemodell)        delete schneemodell;        schneemodell        = NULL;
		if (verdunstungsmodell)  delete verdunstungsmodell;  verdunstungsmodell  = NULL;
		if (expokorrektur)       delete expokorrektur;       expokorrektur       = NULL;
		if (superposition)       delete superposition;       superposition       = NULL;
		if (precicorr)           delete precicorr;           precicorr           = NULL; 
		if (specialoutput)       delete specialoutput;       specialoutput       = NULL;
		if (data_assimilation_model) delete data_assimilation_model; data_assimilation_model = NULL;
		if (initial)             delete initial;             initial             = NULL;
		if (LastResults) free(LastResults); LastResults=NULL;

		if (v_mpiregioncnt) delete[] v_mpiregioncnt;
		if (sendHalo) delete sendHalo;
		if (recvHalo) delete recvHalo; 

		if (vHRU)
		{
			for (int i = 0; i < vHRU->size(); i++)
			{
				delete(vHRU->at(i));
			}
			vHRU->clear();
			delete vHRU;
			vHRU = NULL;
		}

	}
	catch (int ex)
	{
		if (WSstream) (*WSstream) << "FATAL: Exception in Modellapplication::ClearMemory() caught!" << "[end-entry]";
		throw ex;
	}
}

void ModellApplication::m_sleep(clock_t wait)
{
	clock_t goal;
	goal = wait + clock();
	while( goal > clock() );
}

bool ModellApplication::WaitForFile(char *filename,long waitinterval)
{
	try
	{
		char key;
		FILE *file;
		if (WSstream) (*WSstream) << "waiting for external synchronisation file " << filename << "...";
		while ((file=fopen(filename,"r"))==NULL)
		{
			if (controll_c_flag!=0)
			{
				controll_c_flag=0;
				(*WSstream) <<"SKIP: \n\nCtrl-C pressed while waiting for " << filename << "\nCancel the model? \n"
				"0 = no (continue)\n"
				"2 = cancel now without saving buffers!\n"
				"\nchoose one of the above options: "<<"[end-entry]";
				cin >> key;
				controll_c_flag=atoi(&key);
				//        while (cin.peek()!='\n') cin.get(key);
				if ((controll_c_flag<0)||(controll_c_flag>3)) controll_c_flag=0;
				if (controll_c_flag==2)
				{
					throw(98);
				}
				SetSignalAbfangfunktion();
			}
#ifdef WIN32
			//		Sleep(waitinterval);	
			//_sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
			Sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
			// however, _sleep is a deprecated C++ feature, so it must be replaced somewhen in the future
#elif WIN64
			//_sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
			Sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
#else
			//m_sleep(waitinterval);
			timespec tim;
			if (waitinterval>1000)
			{
				tim.tv_sec = (unsigned long)waitinterval/(unsigned long)1000;
				tim.tv_nsec = 1000000 * (waitinterval%1000);
			}
			else
			{
				tim.tv_sec = 0;
				tim.tv_nsec = 1000000 * waitinterval;
			}
			nanosleep(&tim, NULL);
#endif
			SetSignalAbfangfunktion();
		};
		if (file)
		{
			if (WSstream) (*WSstream) << "found it! " << "[end-entry]";
			char *content = new char[_stringsize_ + 1];
			fgets(content, _stringsize_, file);
			std::string commandstring(content);
			if (commandstring.find("Command=") != std::string::npos)
			{
				int index1 = commandstring.find("Command=") + 8;
				int index2 = commandstring.substr(index1).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_");
				commandstring = commandstring.substr(index1, index2);
				if (commandstring.find("SaveAndContinue") != std::string::npos)	{
					controll_c_flag = 3;
				}
				else if (commandstring.find("Continue") != std::string::npos) {
					controll_c_flag = 0;
				}
				else if (commandstring.find("Finish") != std::string::npos) {
					controll_c_flag = 1;
				}
				else if (commandstring.find("StopWithoutSaving") != std::string::npos) {
					controll_c_flag = 2;
				}
			}
			delete[] content;
			fclose(file);
		}
		return true;
	}
	catch (int ex)
	{
		(*WSstream) << "FATAL: Exception in Modellapplication::WaitForFile() caught!" << "[end-entry]";
		throw ex;
	}
}

void ModellApplication::WriteSemaphore(const char *filename,char *content)
{
	ofstream file;
	file.open("WaSiM_ETH_temporary_semaphore.file");
	file << content;
	file.close();
	rename("WaSiM_ETH_temporary_semaphore.file",filename);
}

// routingTimestep = time step in minutef routing
// timestep = time step in unsatzonmodell
// routingTimestepFactor = externalTimeStep / routingTimestep;
// minutes_to_interval_end = 0 or time step length in minutes for normal modules
// MinutesToInvEndRouting   = 0 or time step length in minutes for routing
// the caller should clean up the memory by calling delete
std::vector<std::pair<int,int>>* ModellApplication::getSubSteps(int timestep, int routingTimestep, int minutes_to_interval_end,int MinutesToInvEndRouting,int year,int month,int day,int hour,int minute)
{
	std::vector<std::pair<int,int>> *timestepvector = new std::vector<pair<int,int>>;

	// todo: calculate all minute and hour values as pairs for the inbetween time steps based on the 4 parameters transferred to this method
	// 1. get start time value of time step end

	time_t routingtime;
	time(&routingtime); // initialising routingtime with the actual UTC time in order to prevent a call to gmtime with an invalid time argument
	tm* tm1 = gmtime(&routingtime); // initialize tm1 with the actual calendar date and time (this call is to get a valid pointer to tm1 without having to allocate it seperately
	tm1->tm_hour = (int)((60.0*hour + minutes_to_interval_end - timestep)/60.0); // 0...23 --> if hour = 24, then (60*24 + 0 - 1440)/60 = 0, should usually be 0 for daily time steps.
	tm1->tm_min = minute; // 0..59
	tm1->tm_sec = 0; // 0..59
	tm1->tm_mon = month - 1;  // 0..11
	tm1->tm_year = year - 1900; // e.g. 2010 = 110
	tm1->tm_mday = day;
	tm1->tm_wday = 0;
	tm1->tm_yday = 0;
	tm1->tm_isdst = 0;

#ifdef WIN64
	routingtime = _mkgmtime(tm1);
#else
	routingtime = timegm(tm1);
#endif

	tm1 = gmtime(&routingtime);

	// 2. get number of sub steps
	int steps = timestep / routingTimestep;

	// 3. iterate till the end of the time step and create new elements for minute and hour in the result vector push_back(pair<int,int>(x,y))
	// sample input: 1984 01 01 12 (00), i.e. 12 is the hour value of the middle of the day 
	// sample output: 
	// 1984 01 01 00 30
	// 1984 01 01 01 00
	// ...
	// 1984 01 01 23 30
	// 1984 01 01 24 00

	for (int i = 0; i < steps; i++)
	{
		int rminute;
		int rhour;
		if (routingTimestep < 60)
		{
			rminute = tm1->tm_min + routingTimestep - MinutesToInvEndRouting;
			rhour = tm1->tm_hour;
		}
		else
		{   // we assume no other time steps than whole hours (one or multiple, but not e.g. 90 min)
			rminute = 0;
			rhour = tm1->tm_hour + (routingTimestep - MinutesToInvEndRouting) / 60;
		}

		timestepvector->push_back(std::pair<int, int>(rminute, rhour));
#ifdef _WASIMDEBUG_
		std::stringstream ss;
		ss << "time in seconds since 1970-01-01 UTC: " << routingtime << " year: " << tm1->tm_year << " month: " << tm1->tm_mon << " day: " << tm1->tm_mday << " hour: " << tm1->tm_hour << " timestepvector[1]: " << timestepvector->back().first << " timestepvector[1]: " << timestepvector->back().second;
		(*WSstream) << ss.str() << "[end-entry]";
#endif

		routingtime += routingTimestep * 60.0;

		tm1 = gmtime(&routingtime);

	}
	return timestepvector;
}