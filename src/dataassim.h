
#ifndef _DATAASSIMH_
#define _DATAASSIMH_

#include <iostream>
#include <ctime>
#include "teilmod.h"


class data_assimilation;

class data_assimilation_entity
{
public:	
	wasim::mystream *WSstream;
	data_assimilation *parent;
	std::string assimtype;

	std::string entityString;
	gridinfotype *gridinfo;
	double xll, yll, cs, ncols, nrows, mpistartrow, mpiendrow;
	StringUtils su;
	int mina, sta, tga, mta, jra; // start time (usually only available in meteo interpolators, not in derived clases stemming from teilmodellklasse, but required here for scanning the input list file foer the starting date
	int minute, hour, day, month, year; // actual values read from input list file
	int oldmonth;
	std::ifstream infile;
	std::string infilename;
	std::string listFileEntry;
	time_t startdatetime,actdatetime, nextdatetime;
	double min_to_int_end;
	int InV;
	data_assimilation_entity(data_assimilation *p_parent, wasim::mystream *Stream,
		std::string entityString, gridinfotype *gridinfo,
		double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int mina,int sta, int tga, int mta, int jra, double min_to_int_end,int InV);
	virtual ~data_assimilation_entity();
	// pure virtual function; only derived classes that implement run() can be instanciated
	virtual void run(long me,long se, long tg, long mt, long jr) = 0;

};

class data_assimilation_grid_entity :
	public data_assimilation_entity
{
public:
	std::string entityName;
	double defaultLiquidPerc;
	bool doSplitSWE;
	bool doWeighByGlacierizedFraction;
	double xcoord, ycoord, elevation;
	std::string entityNameInFile;
	int gridindex,
		gridindexliquid,
		gridindexsolid,
		gridindexTotalSWEonglac,
		gridindexsolidSWEOnglac,
		gridindexliquidSWEOnglac,
		gridindexglacfrac;

		double **grid,
		**SWELIQUID,
		**SWESOLID,
		**SWELIQGL,
		**SWESOLGL,
		**SWETOTGL,
		**GLACFRAC;

	firstlast *FirstLastCol;
	double **zn;
	double **MPIREG;
	int mpi_rank;
	int *zarray;
	double nodata;

	data_assimilation_grid_entity(data_assimilation *p_parent, wasim::mystream *Stream,
		std::string entityString, gridinfotype *gridinfo,
		double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow,int mina, int sta, int tga, int mta, int jra,double min_to_int_end,int InV);
	~data_assimilation_grid_entity();
	virtual void run(long me,long se, long tg, long mt, long jr);
};

class data_assimilation_reservoir_entity :
	public data_assimilation_entity
{
public:
	int lakeID;
	SP_class *SP;
	double AverageWaterTable;

	double xcoord, ycoord, elevation;
	std::string entityNameInFile;

	firstlast *FirstLastCol;
	double **zn;
	double **MPIREG;
	int mpi_rank;
	int *zarray;
	double nodata;
	std::string sFormat;
	std::string sLine;
	char *tstr;
	int skiprows;

	double waterlevel;

	data_assimilation_reservoir_entity(data_assimilation *p_parent, wasim::mystream *Stream,
		std::string entityString, gridinfotype *gridinfo,
		double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int mina, int sta, int tga, int mta, int jra, double min_to_int_end,int InV,std::vector<SP_class*> &reservoirs);
	~data_assimilation_reservoir_entity();
	virtual void run(long me,long se, long tg, long mt, long jr);
};

class data_assimilation:
public teilmodellklasse
{
public:

	int numrules;
	std::vector<data_assimilation_entity*> vdaRules;
	long *lakeCellsPerMPIRank;
	int pondindex, lakeindex;

	// constructor of the class
	data_assimilation(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
						double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
						char *steuername,int LogLevel,int progaus,
						double p_min_to_int_end,int p_numsubb,int *p_listsubb,
						gridinfotype *p_gridinfo,
						int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
						int p_zonmax,  int p_zonmaxwert,
						int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
						double p_xll ,  double p_yll,
						double p_cellsize    , double      p_nodata,
						long mina, long sta, long tga, long mta, long jra,
						long mine, long ste, long tge, long mte, long jre,
						double p_timekorr,
						LUSEtabletype *p_LUtable,
						SOILtabletype *p_SOILtable,
						MultiLayerLUSEType *p_MultiLayerLUTable,
						bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
						long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
						int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU, std::vector<SP_class*> &reservoirs, long *lakeCellsPerMPIRank, int substepcnt=1);

	
	virtual void   run(long me, long se,long tg,long mt,long jr);

	virtual void   ClearMemory();
	virtual void   flushbuffers();

	/* overwriting the inherited functions for substance transport from teilmodellklasse */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */

};

#endif