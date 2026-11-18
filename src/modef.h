/*   $Revision: 71 $ 
	$Date: 7.03.11 0:36 $  */
//   (c) 1995-2007 by Joerg Schulla
//   more version history: see source safe comments


#ifndef MODEFH
#define MODEFH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <map>
#include <list>
#include <math.h>
#include <string>
#include "wsstream.h"
#include "modconst.h"
#ifdef _MPI
#include <mpi.h>
#endif
#ifdef _NETCDF
#include "netcdfhlp.h"
#endif

typedef void (*functionpointer) (int);

struct HRUMapElement
{
	// holds information on hydroligic response units used for spin up mode for soil temperatures (and other entities)
	unsigned long ID; // unique ID, independent on sID
	std::string sID; // string representation of the chained values of each (reclassified) grid used to create the HRUs, e.g. 0001'0002'0003'0001 for DHM class=01, Soil=02, landuse=03, slope class=01 (without ')
	// for fast reference to the representative HRU when filling in the other cells
	unsigned int MasterRow; 
	unsigned int MasterCol;
	bool done; // flag that indicates that this HRU was already processed
};

struct mpihalo 
{
	int row; int col; int cnt; int rank;
};

struct  statrecord    
{ 
	double n;
	double sx;
	double zonenwert;
	double rel;
	double RelCellSum;
}; /* for zonal statistics */

typedef int     Tzarray;

struct TLastResults
{
	long j,m,t,s,min;
	long NumZones;
	double *wert;
	double *refvalue;
	double *rel;
	char *Name;
	unsigned long *ZoneCodes;
};

class zonalstatistic;
/* ein Typ, der Variableninformationen ueber die Speichergrids aufnimmt */
struct gridinfotype
{
	double **pointer;              /* Zeiger auf das Grid */
	double **summenpointer;       /* Zeiger auf das Summengrid */
	double  **minpointer;          /* Zeiger auf ein Minimumgrid */
	double  **maxpointer;          /* Zeiger auf ein Maximumgrid */
	double ***stackpointer;        /* Zeiger auf den Gridstapel bzw. ein Grid
								aus double-Arrays, die in senkrechter Richtung
								die Eigenschaften des Bodens speichern */
	int   levels;                 /* Dimension der 1-D-Arrays im Gridstapel */
	char  *name;                  /* Name des Grids */
	char  *summenname;            /* Name des Summengrids */
	char  *minname;               /* Name des Minimumgrids */
	char  *maxname;               /* Name des Maximumgrids */
	char  *identification;        /* Bezeichnung der Groesse im Klartext */
	int    writegrid;             /* Schreibcode fuer normales Grid */
	int    summenwritegrid;       /* Schreibcode fuer Summengrid */
	int    minwritegrid;          /* schreibcode fuer minimumgrid */
	int    maxwritegrid;          /* schreibcode fuer maximumgrid */
	int    minmaxflag;            /* if thif flag is 0, initialization with grid values must be done */
	int    fillcode;              /* 0: kein Auffuellen von Luecken, 1=Auffuellen von Luecken moeglich beim Lesen */
	int    getfirstlevel;         /* usually zero, only 1, if the first level from a stack should not be copied into the pointer field (e.g., if there are already other values stored, which should be printed) */
	double  anzahl;               /* Anzahl an Zeitintervallen zur Mittelwertbildung
								beim Summengrid, wenn Summen-Schreibcode > 3 */

	int    vegWritegrid;          // writecode for vegetation period balance grid
	double **vegSumPointer;       // pointer to a grid with summed up values for active Vegetation 
	char   *vegName;              // output name for grid with vegetation period balance
	int    vegWriteFlag;          // is set externally to 1, grids must be written at next call
	int    vegResetFlag;          // if set externally to 1, all Values in vegSumPointer are reset to 0, vegCnt=0, vegWriteFlag=0
	double  **VEGSTART,**VEGSTOP,**VEGDURAT; // pointers to the grids which are handled in defaultuzr but used in gridinfotype (depending on vegetation layer, it may be the address of VEGSTART1..n, VEGSTOP1..n, VEGDURAT1..n, resp.)

	int    iPer; // periodicity of input grid (to be read each time the period is over
	int    iPerAct; // actual counter of elapsed periods as given in cPerType (H, D, M or Y)
	int    iSubInt; // number of subintervals of the periodtype which must be over until reading a new grid (e.g. 15 for 15 days when cPerType is M
	int    iSubIntAct; // actual counter of elapsed sub intervals until the grid must be imported again
	int    iLastSubInt; // used to determine if the actual sub interval is less than the old one (when a new month starts, the day will be smaller than the last day of the old month)
	int    iLastPer; // same for period counter
	char   *cPerType; // type of periodicity: H, D, M, Y
	char   *cPerGridName; // gridname to read from (will be constructed each time a grid is read from the actual dates)
	bool   bPersistency; // set to true by default. If set to false (by persistent=false in the control file) 
	// the original grid will be stored in a temporary grid and after one timestep the 
	// new read in grid will be replaced by the old one again
	bool   bReplaceGridNow; // is set to true when bPersitency=false and a new grid is read in, so in the next call to the read function it will recognize the flag and replace the grid with the default grid again
	double  **TempGrid;   // here, the old grid will be stored when the bPersistency-Flag is set to false
	int    iReplaceAfter; // number of intervals, the new default grid should be valid until it is replaced by the old default grid (to avoid oscillations)
	int    iReplaceCnt; // actual counter of intervals until iReplaceAfter is reached

	double defaultValue; // value for nodata celss when grid is created or no neighbor was found (depending on the fillcode: 1 = fill with nearest neighbor, 2=fill with default (also when reading in a grid which is not completely filled and has some holes)
	char *outname; // name or directory for writing the grid. This is not identical to the read location as was the case until now  with variable grids
	int statcode; // statistics write code (e.g. 2001 for writing the sum for each time step and each sub basin)
	char *statname; // statistics output file, when configured to write statistics, too
	statrecord *zone; // statistics regristers per zone (.sx, .n, .zonenwert, .rel, .RelCellSum)
	zonalstatistic *stat; // handling the output file, including the file header, buffering etc.
	double **zn; // zone grid reference
	Tzarray *zarray; // reference to the zone index array
	int readcode; // if grid should be read in or otherwise created internally
	int outnameinternally; // flag if the outname was reserved internally in order to free it accordingly
	long int subIntervalCnt; // if time steps of less than 1 hour are used, this counter is used to extend the filename by an additional index which is incremented until a new time stamp is encountered
	long int lastYear, lastMonth, lastDay, lastHour, lastMinute; // stores the date/time of the last written grid. If the next write command encounters the same values, subIntervalCnt is incremented, alse set to 0.

	// memory minimization for MPI related variables (Nov. 2015)
	long nrl, nrh, ncl, nch; // usually, these numbers will be given as parameters. But it may be a good idea to keep them here
	bool withHalo;

#ifdef _NETCDF
	// If given in the input file, read latitude and longitude grids to
	// be able to write them into the netcdf output grids. However, only
	// grid 1 will store this information.
	bool has_latgrid_, has_longrid_;
	int latgrid_index_, longrid_index_;

	// Constructor
	gridinfotype() {
		has_latgrid_ = false;
		has_longrid_ = false;
	}
#endif

};

/* for storing the centers of mass for all regions in order to be able to look for the nearest region for stations outside the model domain */
struct  center_of_mass_type 
{ 
	int row;
	int col; 
	double ycoord;
	double xcoord;
	double n;
};

struct  Tstatist      
{ 
	double ndiff;
	double sxdiff;
	double sx2diff;
	double ngem;
	double sxgem;
	double sx2gem;
};

class StringUtils
{
public:
	void stripInt(std::string &row, std::string sWhatParam, int *param);
	void stripDouble(std::string &row, std::string sWhatParam, double *param);
	void stripString(std::string &row,std::string sWhatParam,std::string &param);
	std::string stripSectionBetweenBrackets(std::string &row);
	int stripDoubleVector(std::string &row,std::string sWhatParam,std::vector<double> &X);
	int stripIntVector(std::string &row,std::string sWhatParam,std::vector<int> &X);
	size_t StringSplit (const std::string &text, const std::string &separators, std::vector< std::string > &words, bool bAddEmptyStrings );
};

class CommonTableType
{
protected:
	char *Name;
	int Code;
	CommonTableType() {Name=new char[200];};
	CommonTableType(int code) {Code = code; Name=new char[200];};
	~CommonTableType();
	StringUtils su;
};

class MultiLayerLUSEType:
public CommonTableType
{
public:
	int *LULayers;
	double *ArcVegLUFractions;
	double* ArcVegBiomass;
	double* ArcVegRetrans;
private:
	int iNumLayers;
	int iNumFractions;
	int iNumBiomass;
	int iNumRetrans;
	int *OriLULayers;
	int OriCode;
	double K_lambert; // extinction coefficient of d'Lambert-Beer's law for reducing radiation in its way through the canopies
	double LAI_scale; // Scaling factor for calculating the aerodynamic resistencies of layer 2..n dependent on the cumulated leaf area index (after COUPmodel, Jansson and Karlberg)

public:
	void setCode(int i) {Code = i;};
	int getCode() {return Code;};
	int getOriCode() {return OriCode;};
	int getNumLayers() { return iNumLayers; };
	int getNumFractions() { return iNumFractions; };
	int getNumBiomass() { return iNumBiomass; };
	int getNumRetrans() { return iNumRetrans; };

	int* getLULayers() { return LULayers; };
	int getLULayer(int i) { if ((i <= iNumLayers) && (i > 0)) return LULayers[i]; else return -9999; };
	void setLULayer(int j, int k) { if ((j > 0) && (j <= iNumLayers)) { LULayers[j] = k; } };

	double* getArcVegLUFractions() { return ArcVegLUFractions; };
	double getArcVegLUFraction(int i) { if ((i <= iNumLayers) && (i > 0)) return ArcVegLUFractions[i]; else return -9999; };
	void setArcVegLUFractions(int j, double k) { if ((j > 0) && (j <= iNumLayers)) { ArcVegLUFractions[j] = k; } };

	double* getArcVegBiomass() { return ArcVegBiomass; };
	double getArcVegBiomass(int i) { if ((i <= iNumBiomass) && (i > 0)) return ArcVegBiomass[i]; else return -9999; };
	void setArcVegBiomass(int j, double k) { if ((j > 0) && (j <= iNumBiomass)) { ArcVegBiomass[j] = k; } };

	double* getArcVegRetrans() { return ArcVegRetrans; };
	double getArcVegRetrans(int i) { if ((i <= iNumRetrans) && (i > 0)) return ArcVegRetrans[i]; else return -9999; };
	void setArcVegRetrans(int j, double k) { if ((j > 0) && (j <= iNumRetrans)) { ArcVegRetrans[j] = k; } };

	int *getOriLULayers() {return OriLULayers;};
	int getOriLULayer(int i) {if ((i<=iNumLayers)&&(i>0)) return OriLULayers[i]; else return -9999;};
	void setK_lambert(double kl){K_lambert = kl;};
	void setLAI_scale(double LAIsc){LAI_scale = LAIsc;};
	double getK_lambert() {return K_lambert;};
	double getLAI_scale() {return LAI_scale;};

	MultiLayerLUSEType();
	MultiLayerLUSEType(int code);
	~MultiLayerLUSEType();
	void initialize(wasim::mystream *WSstream,std::string row,int version=1);
	void stripParameters(std::string& row, std::string sWhatParam, int** param, int& NumLayers);
	void stripParameters(std::string& row, std::string sWhatParam, double** param, int& NumLayers);

};


class LUSEtabletype:
public CommonTableType
{
private:
	int method;
	int iNumFixPoints;
	int iNumMTStartLevels;
	int StartVegetationPeriod;
	int StopVegetationPeriod;
	int JDVegetationReset;
	int JDVegetationWrite;

	unsigned int dynflags; // flags: used for swithcing on/off the dynamic calculation of 
	std::string sDynFlags;
	// 0x01 Albedo	--> sDynFlags = 00000001 when set to 1, the entity is modelled internally, when set to 0, it should be read from a grid periodically
	// 0x02 RSE		--> sDynFlags = 00000010 --"-- dito...
	// 0x04 RSC		--> sDynFlags = 00000100 --"-- dito...
	// 0x08 RSI		--> sDynFlags = 00001000 --"-- dito...
	// 0x10 LAI		--> sDynFlags = 00010000 --"-- dito...
	// 0x20 VCF		--> sDynFlags = 00100000 --"-- dito...
	// 0x40	ROOT	--> sDynFlags = 01000000 --"-- dito...
	// 0x80 Z0		--> sDynFlags = 10000000 --"-- dito...
	// all togeter  --> sDynFlags = 11111111 --> eill be translated into an integer for easier matching
	// or in any combination. Note: if a flag is set to 0, the corresponding grid should be read in dynamically

	// MTSTart, TSTART and FSTAR are grids (and one stack) which are handled by defaultini (int_tvar_params)
	// but used in LUtable-objects. 
	double ****MTSTART;    // Start days for subsequent phenologic phases (each sampling point has its own starting day)
	double ***TSTART;      // Start day t1 for dynamic_phenology_1, each layer gets one grid
	double ***FSTAR;       // F* for dynamic phenology 2 (F3 will be calculated there depending on chilling units, so its dynamically for each cell as well as FSTART
	double ***FORC;        // actual forcing units for dynamic phenology 3 for interpolating actual phenologic values depending on forcing units and last and next sample days forcing thresholds

	int Tstart;			// julian day which was calculated using dynamic phenology
	int Tstart_col;		// which column contains TStart?
	int *days;	        // Stützstellen der Jahresgang-Kurve 
	int *daystemp;        // nach Höhe korrigierten aktuellen Stützstellen 
	double	*rsc,			// minimal canopy-resistance der Landnutzung als Monatswerte 
			*rsi,			// minimal interception-resistance der Landnutzung als Monatswerte 
			*rse,			// minimal soil-resistance der Landnutzung als Monatswerte 
			*albedo,		// albedo if no snow is covering the crops                     
			*LAI,			// Blattflaechenindizes an den 4 Knickpunkten des Jahresganges 
			*z0,			// effektive Bewuchshoehe der Landnutzung zu den Knickpunkten  
			*Vbg,			// Vegetationsbedeckung der Landnutzung an den Knickpunkten    
			*Wurz,			// Wurzeltiefe der Landnutzung                                 
			*AltDep,		// Altitude dependency in days per m                           
			*ForcingThreshold; // ForcingThreshold for each 
	double form_par,		// Wurzeldichte fuer Entnahmefunktion 0...                     
		redu_parDry,    // ab welcher Bodenfeuchte im Trockenbereich beginnt Verdunstungsreduktion einlesen als suction in m, Umrechnen auf absolten Wassergehalt Theta 
		redu_parWet,    // ab welcher Bodenfeuchte im Nassbereich beginnt Verdunstungsreduktion, einlesen als relativer Wassergehalt (0..1), Umrechnen auf absoluten Wassergehalt Theta 
		LimitReduWet,   // minimal relative reduction of real transpiration when water content reaches saturation. The reduction factor value will go down linearly starting at 1.0 when relative Theta equals TReduWet (e.g. 0.95) to LimitReduWet when the soil is saturated (Theta rel = 1.0) 
		ETPHamonFactor; // some evaporation methods underestimate potential evaporation, so it is possible here per landuse to set a correction factor
	double F_star;         // "Temperatursum" which must be exceeded for starting the phenological cycle (unfolding leaves)       
	int   JDReset_TStart; // Julian Day when TStart is reset to -1 and Forcing untis are reset to 0 for a new vegetation period  
	int   DP1_t1_dorm;    // starting day (julian day number), forcing units will be summed up after this day of year, thermal time model after equation 12b
	double DP1_T_Bf;	    // threshold temperatur for a positive forcing unit after Model 12b (thermal time model = DP1) 
	int   DP2_t0_dorm;    // starting day (julian day number), chilling units will be summed up after this day of year, sequential model 2 after equation 24b
	int   DP2_t1_dorm;    // starting day (julian day number), forcing units will be summed up after this day of year, sequential model 2 after equation 24b
	double DP2_T_Bc;       // threshold temperatur for a chilling unit after Model 24b (Sequential model 2 = DP2)
	double DP2_T_Bf;       // threshold temperatur for a positive forcing unit after Model 24b (Sequential model 2 = DP2)
	double DP2_Par_a;		// Parameter a in F*=a*exp(bC*) after Model 24b (sequential model 2)
	double DP2_Par_b;		// Parameter b in F*=a*exp(bC*) after Model 24b (sequential model 2)
	double DP2_Offset_1;	// used for DynamicPhenology_2 only: value for z1 in R_c(T_i)=(T_i-z1)/(T_Bc-z1) when z1 < T_i < T_Bc
	double DP2_Offset_2;	// used for DynamicPhenology_2 only: value for z2 in R_c(T_i)=(T_i-z2)/(T_Bc-z2) when T_Bc < T_i < z2

	int   DP4_t0_dorm;    //  starting day (julian day number), chilling units will be summed up after this day of year
	int   DP4_t1_dorm;	//  starting day (julian day number), forcing units will be summed up after this day of year
	double DP4_T0;			// threshold temperatur for chilling units (if T < T0)
	double DP4_T1;			// threshold temperatur for forcing units (T > T1)
	double DP4_Par_a;		// tree specific parameter
	double DP4_Par_b;		// tree specific parameter
	double DP4_T_xylstop;	// threshold temperature, the moving average of temperature has to stay below for 5 days in order to start leave fall
	int   DP4_t_xs1;		// time of year (julian day number), the moving average will be analysed from
	int   DP4_SPSLF;		// sample point ID which is replaced by the day of start leave fall

	int   maxStartJDforDP1; // latest start day for the model run to use DynamicPhenology_1. 
	// If start date is after this date, then TStart is set to maxStartJDforDP1 
	// minus the delta of the next column (e.g. maxStartJDforDP1(150) - Delta(18) = 132), so we assume 
	// that this start date meets a fully developed vegetation. If start day is 
	// even after DP2_t0_dorm, then the next year will use DP1 only
	bool  useMethod3;		// when starting a new model run, method 4 (sequential model) cannot be used for 
	// calculating TStart when the model start date is less than maxStartJDforDP1 because Chilling units calculation
	// will usually start in autumn of the old year. Then, this bool is set to true until DP2_t0_dorm is reached, which sets this bool to false.
	// by the way: if model start day is after maxStartJDforDP1, then the actual year will be run without TStart calculation at all, neither Method 3 nor method 4,
	// so the method 4 can be applied in the next year only (assuming the model run encompasses at least 2 years)
	double m_fStressFactor;   // used for dynamicPhenology_3: a value of 0 means no changes in dynamic phenology due to water stress (dry soil), 
	// values above 1 mean that the effect of dry soils accellerates the development and 
	// values between 0 and 1 do mean that the stress slows the development down (with 0 stopping the development at all)
	// stress starts from redu_parDry with increasing severity towards zero.
	
	double m_fIntercepStor; // specific interception storage capacity in mm -> will be multiplied by LAI to get the effektive storage capacity

	int *m_SoilTillage; // array of Julian Days containig the Julian Days when the soil is reworked 
	int m_numSoilTillages; // number of Julian Days with soil tillage

	// parameters used for surface routing only
	int SR_cultivation; // default: other; type of cultivating the surface; possible values: conserving (c)=1, nonconserving (nc)=2, other (o)=3
	double *SR_StemDiameter;	// (nc) vector for each sample of the JD-table, default: 0.01 surface routing: diameter of the stalks or effective diameter of crop/tree stems (perpendicular to the streamlines); used for nonconserving land use (using ploughs)
	double *SR_ResidCovrg;	// (c)  vector for each sample of the JD-table, default: 0.75 fraction of soil covered by vegetation residuals (mulch); used for conserving land use only
	double SR_CropDistX;		// (nc) default: 0.3 distance of crops in flow direction; used for nonconserving land use (using ploughs) only
	double SR_CropDistY;      // (nc) default: 0.3 distance of crops perpendicular to flow direction; used for nonconserving land use (using ploughs) only
	double SR_K_ST;           // (nc/o) default: 10 Manning-Strickler-roughness for the soil surface in (m^1/3)/s; used for nonconserving and other landuse only
	double SR_C_WR;           // (nc) default: 1.5 form drag coefficient for a circular cylinder in a group; used for nonconserving landuse only
	double *SR_rough_coeff1;  // (c) default = 0.0198 0.712 0.142; coefficients for k_st-calculation for conserving landuse, Reynolds number < 20000
	double *SR_rough_coeff2;  // (c) default = 0.0198 0.712 0.142; coefficients for k_st-calculation for conserving landuse, Reynolds number >= 20000

	// heat transfer parameters
	double *nfactor;	// a vector of coefficients which are applied to the air temperature to get the surface temperature as boundary condition for the heat transport model
	double nf_freezing,nf_thawing; // alternative n-factors per land use type: specific values for temperatures above and below 0°C which will be used if no monthly values are available
	
	// ArcVeg Parameters
	double ArcVegInFactFull; // Inundation scaling factor for full inundation
	double ArcVegInFactSome; // inundation scaling factor for some inundation
	double ArcVegInUDIceWdg; // inundation scaling factor for other cells that are not on an ice wedge
	double ArcVegBNratio; // biomass to N ration (each in g/m^2)
	double ArcVegSenesce; // fraction of senescing biomass
	double ArcVegRetran;  // plant type specific factor for retranslocated N
	double ArcVegFoliarProp; // foliar prop for each plant type (used for senesce)
	double ArcVegPropN; // PropN is another factor used in ArcVeg
	double ArcVegAGProp; // AGProp is yet another factor used in ArcVeg
	double ArcVegMaxBio; // maxium biomass for a m^2 (used in ArcVeg)
	double *ArcVegSeedMaxProb; // maximum probability in each year that a seedling will establish with low biomass (0.1 g/m2) for all seven climate zones (7 values are read in)
	double *ArcVegUptakeEfficiency; // uptake efficiency for every plant type and every PGP
	int mon[14];

public:
	LUSEtabletype(int index);
	LUSEtabletype(); // defaultconstructor
	~LUSEtabletype();

private:
	double interpolate(double x1,double y1,double x2,double y2,double x);
	void AltitudeCorrection(wasim::mystream *WSstream,double dhm,int tstart);
	int GetFixPoint(wasim::mystream *WSstream,int iJulDay,double height,int tstart);
	void stripParameters(std::string &row,std::string sWhatParam,double *param,int maxIndex);
	void stripParametersJulDays(wasim::mystream *WSstream,std::string &row);
	void stripParametersSoilTillageDays(wasim::mystream *WSstream,std::string &row);

public:

	int GetJulDay(int iMonthDay,int iMonth,int iYear);
	void initialize(wasim::mystream *WSstream,std::string row,int version);
	void initialize(wasim::mystream *WSstream,std::string row,std::string sVersion);

	unsigned int getDynFlags();
	std::string getSDynFlags();
	void	setDynFlags(unsigned int flags);
	unsigned int isFlagSet(unsigned int flag); // returns the "and"-ed result of dynflags and the input parameter. So if a flag is set in the dynflags-array and in the parameter, the matching bit in the result will be set, too 

	int   getStartVegetationPeriod();
	int   getStopVegetationPeriod();
	int   getJDVegetationReset();
	int   getJDVegetationWrite();
	void  setJDVegetationReset(int jd);
	void  setJDVegetationWrite(int jd);

	int   getNumFixPoints();
	double getFormPar();
	double getReduParDry();
	double getReduParWet();
	double getETPHamonFactor();
	double getLimitReduWet();
	double getIntercepStor();
	char* getName();
	int   getCode();
	void  setCode(int code);
	int   GetMethod();

	int   GetTstart();
	void  SetTstart(int day);
	void  SetFStarGridVector(double*** &_FSTAR);
	void  SetTStartGridVector(double*** &_TSTART);
	void  SetFORCGridVector(double*** &_FORC);
	void  SetMTStartStackVector(double**** &_MTSTART,int _iNumMTStartLevels);
	int   getTstartcol();
	double GetFstar();    
	int   GetDP1_t1dorm();   
	double GetDP1_TBf();      
	int   GetDP2_t0dorm();   
	int   GetDP2_t1dorm();   
	double GetDP2_TBf();      
	double GetDP2_TBc();      
	double GetDP2_Par_a();      
	double GetDP2_Par_b();      
	double GetDP2_Offset_1();
	double GetDP2_Offset_2();

	int   GetDP4_t0dorm();   // usually November, 1st; here we start the chilling units to add up
	int   GetDP4_t1dorm();   // usually February, 1st; here we start the forcing units to add up
	double GetDP4_T0();       // threshold temperature for chilling units (T < T0)
	double GetDP4_T1();       // threshold temperature for forcing untis (T > T1)
	double GetDP4_Par_a();    // empirical parameter a
	double GetDP4_Par_b();    // empirical parameter b
	double GetDP4_Txylstop(); // temperature of the moving average which must be under-run for at least 5 days to stop xylem growth
	int   GetDP4_txs1();  // time of year as julian day, the evaluation of the moving average of temperature starts to estimate stop of xylem growth
	int   GetDP4_SPSLF(); // SamplePontStartLeaveFall

	int   GetJDReset_TStart();
	int   GetMaxStartJDforDP1();
	void  SetUseMethod3(bool mode);
	bool  GetUseMethod3();
	void  InitDynamicPhenology(int taga,int monata,int jahra);
	double* GetForcingThresholds();
	double GetForcingThreshold(int i);
	double* GetAltDeps();
	double GetAltDep(int i);
	int*  GetJulianDays();
	int   GetJulianDay(int i);
	double GetStressFactor();

	void GetSoilTillage(int &SoilTillageCount,int* &SoilTillageDays);

	// SR=SurfaceRouting-parameters 
	int getSR_cultivation(); // default: other; type of cultivating the surface; possible values: conserving (c)=1, nonconserving (nc)=2, other (o)=3
	double GetSR_CropDistX();		// (nc) default: 0.3 distance of crops in flow direction; used for nonconserving land use (using ploughs) only
	double GetSR_CropDistY();      // (nc) default: 0.3 distance of crops perpendicular to flow direction; used for nonconserving land use (using ploughs) only
	double GetSR_K_ST();           // (nc/o) default: 10 Manning-Strickler-roughness for the soil surface in (m^1/3)/s; used for nonconserving and other landuse only
	double GetSR_C_WR();           // (nc) default: 1.5 form drag coefficient for a circular cylinder in a group; used for nonconserving landuse only

	void GetSR_rough_coeff1(double &c1,double &c2,double &c3);  // (c) default = 0.0198 0.712 0.142; coefficients for k_st-calculation for conserving landuse, Reynolds number < 20000
	void GetSR_rough_coeff2(double &c1,double &c2,double &c3);  // (c) default = 0.0198 0.712 0.142; coefficients for k_st-calculation for conserving landuse, Reynolds number >= 20000

	double GetSR_StemDiameter(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);	// (nc) vector for each sample of the JD-table, default: 0.01 surface routing: diameter of the stalks or effective diameter of crop/tree stems (perpendicular to the streamlines); used for nonconserving land use (using ploughs)
	double GetSR_StemDiameter(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);	// (nc) vector for each sample of the JD-table, default: 0.01 surface routing: diameter of the stalks or effective diameter of crop/tree stems (perpendicular to the streamlines); used for nonconserving land use (using ploughs)

	double GetSR_ResidCovrg(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);  // (c)  vector for each sample of the JD-table, default: 0.75 fraction of soil covered by vegetation residuals (mulch); used for conserving land use only
	double GetSR_ResidCovrg(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);	// (c)  vector for each sample of the JD-table, default: 0.75 fraction of soil covered by vegetation residuals (mulch); used for conserving land use only

	double GetRsc(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetRsc(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetRsi(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetRsi(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetRse(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetRse(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetZ0(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetZ0(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetVegStart(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetVegStart(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetVegStop(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetVegStop(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetLAI(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetLAI(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetVbg(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetVbg(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetWurz(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetWurz(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetAlbedo(wasim::mystream *WSstream,int iJulDay, double height=0.0,int layer=1,int i=0,int j=0);
	double GetAlbedo(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear, double height=0.0,int layer=1,int i=0,int j=0);

	double GetNfactorJD(wasim::mystream *WSstream,int iJulDay, int i=0,int j=0);
	double GetNfactor(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,int i=0,int j=0);
	double GetNFactorFreezing(wasim::mystream *WSstream,int i=0,int j=0);
	double GetNFactorThawing(wasim::mystream *WSstream,int i=0,int j=0);

	double GetArcVegInFactFull(); // Inundation scaling factor for full inundation
	double GetArcVegInFactSome(); // inundation scaling factor for some inundation
	double GetArcVegInUDIceWdg(); // inundation scaling factor for other cells that are not on an ice wedge
	double GetArcVegBNratio(); // biomass to N ration (each in g/m^2)
	double GetArcVegSenesce(); // fraction of senescing biomass
	double GetArcVegRetran();  // plant type specific factor for retranslocated N
	double GetArcVegFoliarProp(); // foliar prop for each plant type (used for senesce)
	double GetArcVegPropN(); // PropN is another factor used in ArcVeg
	double GetArcVegAGProp(); // AGProp is yet another factor used in ArcVeg
	double GetArcVegMaxBio(); // maxium biomass for a m^2 (used in ArcVeg)

	double *GetArcVegSeeedMaxProb(); // maximum probability in each year that a seedling will establish with low biomass (0.1 g/m2) for all seven climate zones (7 values are read in)
	double *GetArcVegUptakeEfficiency(); // returns vector of uptake efficiencies for all possible 5 plant grow periods
};

// struct  LUSEtabletype
// {
// char *Name;     /* Name der Landnutzungsarten z.B. "Nadelwald"                 */
// int   Code;     /* Code dieser Nutzungsart im Landuse-Grid                     */
// double Albedo;   /* Albedo der entsprechenden Landnutzungsart                   */
// double *rsc,     /* minimaler canopy-resistance der Landnutzung als Monatswerte */
// *d,       /* 4 Knickpunkte der Jahresgang-Kurve fuer LAI, z0 und Vbg     */
// *LAI,     /* Blattflaechenindizes an den 4 Knickpunkten des Jahresganges */
// *z0,      /* effektive Bewuchshoehe der Landnutzung zu den Knickpunkten  */
// *Vbg,     /* Vegetationsbedeckung der Landnutzung an den Knickpunkten    */
// *Wurz;    /* Wurzeltiefe der Landnutzung                                 */
// double form_par, /* Wurzeldichte fuer Entnahmefunktion 0...                     */
// redu_par; /* ab welcher Bodenfeuchte beginnt Verdunstungsreduktion 0..1  */
// };

class SOILtabletype:
public CommonTableType
{ 
private:
	int		parmethode; /* Methode der Parameteruebergabe 1=van-Genuchten, 2=look-up-tables        */

	double dIsArcVegType; // 0=no, other=yes
	bool bIsArcVegType;

	double	nFK,        /* nutzbare Feldkapazitaet in Vol.%                                        */
			mSB,        /* maximaler Bodenspeicher in Vol.% = Feuchte bei SAettigung - Welkepunkt  */ 
			drain,      /* der durch Gravitation entwaesserbare Porenraum mSB-nFK in rel. Anteil   */
			k,          /* vert. gesaettigte hydraulische Wasserleitfaehigkeit des Bodens in m/s   */
			psi;        /* Saugspannung in mm                                                      */
	int		anz_k,      /* Anzahl von folgenden Wertepaaren fuer die k von theta - Beziehung       */
			anz_h,      /* Anzahl von folgenden Wertepaaren fuer die h von theta - Beziehung       */
			levels;     /* Anzahl an Schichten fuer diesen Boden. Alle Schichten haben identische k-theta und h-theta-Beziehungen, aber untershciedliche k_sat-Werte */
	/* die erste Schicht ist im alten WaSiM bis version 6.4 0.05 m dick, die zweite 0.2 m  und ab der dritten alle Schichten je 0.25 m */
	int		iNumHorizons; /* Anzahl Horizonte - momentan nur für das Einlesen der Parameter genutzt, da diese horizontweise eingelesen werden */	
	int		*LevelPerHorizon; 
	double	**theta_k,   /* Werte fuer Wassergehalt in der k-Theta-Tabelle */
			**k_tab,     /* zu den entsprechenden Indizes in theta_k passende relativb Leitfaehigkeiten (0..1) */
			**theta_h,   /* Werte fuer den Wassergehalt in der Saugspannungs-Saettigungsbez. (h-theta) */
			**h_tab,     /* zu den entsprechenden Indizes in theta_h passende Saugspannungen in m */
			*k_sat,     /* gesaettigte hydraulische Leitfaehigkeit fuer die einzelnen Schichten des Bodens in m/s */
			*dzui;      /* Dicke der einzelnen Bodenschichten (in Index 0 wird die allgemeine Schichtdicke eingelesen, dann wird die erste Schicht auf 5 cm, die zweite auf 20, alle weiteren auf 25 cm gesetzt */
	double 
			*thetas,     /* Saettigungswassergehalt */
			thetafc,    /* Feldkapazitaet */ 
			*thetar,     /* residualer Wassergehalt (gilt fuer alle Schichten des Bodens) */
			*psi_FC; // suction at field capacity. Used for new statistics on useable soil water (nFK/uFC)
	bool useLocalDefaultPsiFC; // true, if there are no local psi_FC values defined and the default of 200 hPa will be used. That value can then be overidden by the global psi_FC_global 
	bool useGlobalDefaultPsiFC; // true, if the gloabl value was set
	double
		*n,
		*m,
		*alpha,
		*tau, // mualem-parameter, per default = 0.5 in van-Genuchten-approach (optional parameter).
		*k_rec, // k recession for new method (per horizon)
		k_abnahme,/* parameter fuer van-geneuchten und Abnahme des k-Wertes pro Meter als ko/ku */
		maxratio, /* maximum ratio of ko/ku or ku/ko for each soil type */
		dEvapMaxDepth, // maximum depth for evaporation. Default = first layer
		dMaxInterflowDepth, // maximum depth for interflow: down to this depth the hydr. conductivity is decreased in order to be able to create interflow at all (works only if k_recession < 1)
		dThicknessScaling, // scale factor for soil layewr thickness in *dzui. default=1.0
		dTotalThickness; // for easier access to set aquifer thickness

	// some parameters for macro pores --> only one set for a soil profile, NOT for each layer!
	double	dPMacroThreshold;   // threshold of precipitation intensity for activating the macro pore system
	double	dPMacroCapacity;    // capacity of the macropores in m/s
	double	dMacroDepth;        // depth of the macropores in m below surface
	double	dReductionPerMeter; // ration of dPMacroCapacity in 1m depth to the one at the surface 
	// --> 1.0 means, there is no decreasing capacity at all until the dMacroDepth is reached
	// 0.5 means, in a depth of 1m the macropres can transport only 50% of the water flowing 
	// into the macropore system at the surface. 
	// also numbers > 1 are allowed, thus allowing e.g. reinfiltration from rivers to 
	// flow rapidly into deeper soil layers it can be 
	// parameters for silting up model
	double	dFractionSand, // fraction of sand in the simplified grain size distribution (0..1)
			dFractionSilt, // fraction of silt in the simplified grain size distribution (0..1)
			dFractionClay, // fraction of clay in the simplified grain size distribution (0..1)
			dFractionStones1, // dFractionSand+dFractionSilt+dFractionClay+dFractionStones1..4 should be 1.0
			dFractionStones2, 
			dFractionStones3, 
			dFractionStones4; 

	// optional dynamic parameters for the expression parser
	double SU_PAR[11];

	// parameters for temperature model
	double *DryHeatCapacity;   // dry heat capacity in J/(Kg*K)
	double *DryDensity;		  // dry density in m^3/m^3
	double *DryThermalConduct; // dry thermal conductivity in W/(m*K) (or J/(m*s*K) --> 1J = 1Ws)
	double KMinFrozenSoil;    // reduced k_sat (minimum hydraulic conductivity for fully frozrn soils)
	double lnKMinFrozenSoil;    // log. value of KMinFrozenSoil as used in the calculation of the effective k-value
	double PhaseChngMinTemp;  // minimum temperature, we want to be accounted for for phase change calculation, usually -10.0°C
	double drySoilThermalCondEff; // effective (constant) thermal conductivity for the soil particles part for saturated soil
	double WaterThermalCondEff; // effective (constant) thermal conductivity for the water part for saturated soil of > 0°C
	
	int anz_T; // number of look-Up-table entries for SE(T) table and dEoverdT-Table
public:
	double **SE_T; // lookup table for relative saturation as a function of temperature (per soil layer!)
	double **T_SE; // for reverse lookup, i.e. getting the temperatur for a given SE (used for adapting the temperature to a partly drained snow layer in order to protect the solid snow from melting away too fast due to missing liquid fraction)
private:
	double **dEOverdT; // lookup table for internal energy change (per soil layer!)
	double **Ts_RelTheta; // lookup table for temperature shift depending on water saturation, i.e. (theta_rel-theta_res)/(theta_sat-theta_res)
//	double **ksat_SE; // ksat as function of SE

	void LookUpTableGeneration(int layer);
	void stripSingleParameter(std::string &row,std::string sWhatParam,double *param);
	void stripGrainSizeParameters(std::string &row);
	int stripParameters(std::string &row,std::string sWhatParam,double *param,int &horizons);

public:
	void LookUpTableGenerationTMod(int layer,double factor);

	SOILtabletype(int index);
	SOILtabletype(); // defaultconstructor
	~SOILtabletype();
	void initialize(wasim::mystream *WSstream,std::string row,int version);
	void initialize(wasim::mystream *WSstream,std::string row,std::string sVersion);

	char  *getName() {return Name;}; /* Name der Bodenart, z.B.: "anlehmiger_Sand"                              */
	
	int    getCode() {return Code;}; /* Code dieser Bodenart im Bodenarengrid                                   */
	void   setCode(int code) {Code=code;}; /* Code dieser Bodenart im Bodenarengrid                                   */
	bool getIsArcVegType() { return bIsArcVegType; }; // if this soiltype should be handelt as a soil type with variable uppermost soil layers (3 of them) is defined by this optional attribite.
	int    getParmethode() {return parmethode;}; /* Methode der Parameteruebergabe 1=van-Genuchten, 2=look-up-tables        */
	double  getNFK() {return nFK;};        /* nutzbare Feldkapazitaet in Vol.%                                        */
	double  getMSB() {return mSB;};        /* maximaler Bodenspeicher in Vol.% = Feuchte bei SAettigung - Welkepunkt  */
	double  getDrain() {return drain;};      /* der durch Gravitation entwaesserbare Porenraum mSB-nFK in rel. Anteil   */
	double  getK() {return k;};          /* vert. gesaettigte hydraulische Wasserleitfaehigkeit des Bodens in m/s   */
	double  getPsi() { return psi; };        /* Saugspannung in mm                                                      */
	void    setPsiFCDefault(double psiFCdefault, int usePsiFCglobalForUnassignedOnly);
	double  *getPsi_FC() { return psi_FC; };        /* Saugspannung bei FC in mm                                                      */
	int    getAnz_k() {return anz_k;};      /* Anzahl von folgenden Wertepaaren fuer die k von theta - Beziehung       */
	int    getAnz_h() {return anz_h;};      /* Anzahl von folgenden Wertepaaren fuer die h von theta - Beziehung       */

	int    getLevels() {return levels;};     /* Anzahl an Schichten fuer diesen Boden. Alle Schichten haben identische k-theta und h-theta-Beziehungen, aber untershciedliche k_sat-Werte */
	void   setLevels(int _levels) {levels=_levels;};     /* Anzahl an Schichten fuer diesen Boden. Alle Schichten haben identische k-theta und h-theta-Beziehungen, aber untershciedliche k_sat-Werte */

	//	  double *getTheta_k() {return theta_k;}; /* Werte fuer Wassergehalt in der k-Theta-Tabelle */
	//    double *getK_tab() {return k_tab;};     /* zu den entsprechenden Indizes in theta_k passende relativb Leitfaehigkeiten (0..1) */
	//    double *getTheta_h() {return theta_h;};   /* Werte fuer den Wassergehalt in der Saugspannungs-Saettigungsbez. (h-theta) */
	//    double *getH_tab() {return h_tab;};     /* zu den entsprechenden Indizes in theta_h passende Saugspannungen in m */
	void CalcSE(double *SE_sl,double *temper_sl,double factor,double *theta_sl); // calculating new values for relative saturation with water (only ice = 0, only water = 1)
	void CalcdEOverdT(double *dEOverdT_sl,double *temper_sl,double *theta_sl);
	void CalcSEanddEOverdT(int offset,double* SE_sl, double *dEOverdT_sl, double *temper_sl, double factor, double *theta_sl, double *thetares,double *thetasat);
	void CalcHeatCapacity(int offset,double *temper_sl, double *heatcap_sl, double *SE_sl, double *dEOverdT_sl, double *theta_sl, double dTModRhoWater, double dTModCWater, double dTModCIce0, double dTModCIce20); // calculating effective heat capacity for the given temperature and absolute and effective saturation, regarding the effect of phase change
	void CalcKOfSE(int offset,double *k_sat_of_SE,double *SE_sl,double *theta_sl,double *thetasat,double *thetares);
	double GetTemperatureShift(int level,double actTheta);

	double *getK_sat(); /* gesaettigte hydraulische Leitfaehigkeit fuer die einzelnen Schichten des Bodens in m/s */
	double *getDzui() {return dzui;};      /* Dicke der einzelnen Bodenschichten (in Index 0 wird die allgemeine Schichtdicke eingelesen, dann wird die erste Schicht auf 5 cm, die zweite auf 20, alle weiteren auf 25 cm gesetzt */
	
	double getThetas(int layer) {return thetas[layer];};     /* Saettigungswassergehalt */
	double *getThetas() {return thetas;};     /* Saettigungswassergehalt als Vektor */

	double getThetafc() {return thetafc;};    /* Feldkapazitaet */ 
	
	double getThetar(int layer) {return thetar[layer];};     /* residualer Wassergehalt (gilt fuer alle Schichten des Bodens) */
	double *getThetar() {return thetar;};     /* residualer Wassergehalt (gilt fuer alle Schichten des Bodens)  als Vektor */
	
	double getN(int layer) {return n[layer];};
	double getM(int layer) {return m[layer];};
	double getAlpha(int layer) {return alpha[layer];};
	double getK_abnahme() {return k_abnahme;}; /* parameter fuer van-geneuchten und Abnahme des k-Wertes pro Meter als ko/ku */
	double getMaxratio() {return maxratio;};  /* maximum ratio of ko/ku or ku/ko for each soil type */

	double k_calc(int layer,double theta,double SE, double thetarex=-1, double thetasex=-1);
	double h_calc(int layer,double theta,double thetarex=-1, double thetasex=-1);
	double get_theta(int layer,double h, double thetarex=-1, double thetasex=-1);

	double	get_PMacroThreshold() {return dPMacroThreshold;};		// threshold of precipitation intensity for activating the macro pore system
	double	get_PMacroCapacity() {return dPMacroCapacity;};		// capacity of the macropores in m/s
	double	get_MacroDepth() {return dMacroDepth;};					// depth of the macropores in m below surface
	double	get_ReductionPerMeter() {return dReductionPerMeter;};	// ration of dPMacroCapacity in 1m depth to the one at the surface 

	double get_SU_PAR(int i) {
		if ((i>=0)&&(i<=10))
		{
			return SU_PAR[i];
		}
		else
		{
			return 0;
		}
	};

	void getGrainSizeDist(double &sand,double &silt,double &clay,double &stones1,double &stones2,double &stones3,double &stones4)
	{
		sand = dFractionSand;
		silt = dFractionSilt;
		clay = dFractionClay;
		stones1 = dFractionStones1;
		stones2 = dFractionStones2;
		stones3 = dFractionStones3;
		stones4 = dFractionStones4;
	}

	// returning pointers to the arrays
	double *getDryHeatCapacity()	{return DryHeatCapacity;};  
	double *getDryDensity()			{return DryDensity;};		 
	double *getDryThermalConduct()	{return DryThermalConduct;};

	// returning single layer values
	double getDryHeatCapacity(int i)	{return DryHeatCapacity[i];};  
	double getDryDensity(int i)			{return DryDensity[i];};		 
	double getDryThermalConduct(int i)	{return DryThermalConduct[i];};

	double getKMinFrozenSoil()      {return KMinFrozenSoil;};
	double getLnKMinFrozenSoil()    {return lnKMinFrozenSoil;};
	double getAnz_T()               {return anz_T;};

	void setDryHeatCapacity(int i,double _p)	{DryHeatCapacity[i]=_p;};  
	void setDryDensity(int i,double _p)			{DryDensity[i]=_p;};		 
	void setDryThermalConduct(int i,double _p)	{DryThermalConduct[i]=_p;};
	void setKMinFrozenSoil(double _p)		{KMinFrozenSoil=_p;};
	void setLnKMinFrozenSoil(double _p)		{lnKMinFrozenSoil=_p;};

	double getDrySoilThermalCondEff(int i) // effective (constant) thermal conductivity for the soil particles part for saturated soil
	{
		if (drySoilThermalCondEff <= 0) 
		{
			drySoilThermalCondEff =  pow(this->getDryThermalConduct(i),1-this->getThetas()[i]);
		} 

		return drySoilThermalCondEff;
	};

	double getWaterThermalCondEff(int i,double LambdaWater) // effective (constant) thermal conductivity for the water part for saturated soil of > 0°C
	{
		if (WaterThermalCondEff <= 0) 
		{
			WaterThermalCondEff  =  pow(LambdaWater,this->getThetas()[i]);
		} 

		return WaterThermalCondEff;
	}; 

	double getEvapMaxDepth() {return dEvapMaxDepth;};
	double getThicknessScaling() { return dThicknessScaling; };
	double getTotalThickness() { return dTotalThickness; };


};

// struct SOILtabletype  
// { 
//   char  *Name;       /* Name der Bodenart, z.B.: "anlehmiger_Sand"                              */
//   int    Code,       /* Code dieser Bodenart im Bodenarengrid                                   */
//          parmethode; /* Methode der Parameteruebergabe 1=van-Genuchten, 2=look-up-tables        */
//   double  nFK,        /* nutzbare Feldkapazitaet in Vol.%                                        */
//          mSB,        /* maximaler Bodenspeicher in Vol.% = Feuchte bei SAettigung - Welkepunkt  */
//          drain,      /* der durch Gravitation entwaesserbare Porenraum mSB-nFK in rel. Anteil   */
//          k,          /* vert. gesaettigte hydraulische Wasserleitfaehigkeit des Bodens in m/s   */
//          psi;        /* Saugspannung in mm                                                      */
//   int    anz_k,      /* Anzahl von folgenden Wertepaaren fuer die k von theta - Beziehung       */
//          anz_h,      /* Anzahl von folgenden Wertepaaren fuer die h von theta - Beziehung       */
//          levels;     /* Anzahl an Schichten fuer diesen Boden. Alle Schichten haben identische k-theta und h-theta-Beziehungen, aber untershciedliche k_sat-Werte */
//                      /* die erste Schicht ist 0.05 m dick, die zweite 0.2 m  und ab der dritten alle Schichten je 0.25 m */
//   double *theta_k,   /* Werte fuer Wassergehalt in der k-Theta-Tabelle */
//          *k_tab,     /* zu den entsprechenden Indizes in theta_k passende relativb Leitfaehigkeiten (0..1) */
//          *theta_h,   /* Werte fuer den Wassergehalt in der Saugspannungs-Saettigungsbez. (h-theta) */
//          *h_tab,     /* zu den entsprechenden Indizes in theta_h passende Saugspannungen in m */
//          *k_sat,     /* gesaettigte hydraulische Leitfaehigkeit fuer die einzelnen Schichten des Bodens in m/s */
//          *dzui;      /* Dicke der einzelnen Bodenschichten (in Index 0 wird die allgemeine Schichtdicke eingelesen, dann wird die erste Schicht auf 5 cm, die zweite auf 20, alle weiteren auf 25 cm gesetzt */
//   double thetas,     /* Saettigungswassergehalt */
//          thetafc,    /* Feldkapazitaet */ 
//          thetar;     /* residualer Wassergehalt (gilt fuer alle Schichten des Bodens) */
//   double n,
//          m,
//          alpha,
//          k_abnahme,/* parameter fuer van-geneuchten und Abnahme des k-Wertes pro Meter als ko/ku */
//          maxratio; /* maximum ratio of ko/ku or ku/ko for each soil type */
// };


struct irrigationtabletype 
{ 
	int     Code;    /* Code for the corresponding grid containing irrigation codes */
	char   *Name;    /* name of the land use/crop type or other to be irrigated */
	int     methode; /* 0=no irrigation, 1=irrigation from table, 2=irrigation according to demand (optimal management) */
	int     origin;  /* 1=from groundwater, 2=from rivers/lakes/reservoirs (surface water) */ 
	double  psi_dry, /* maximum allowed suction */
	psi_wet, /* suction to which soil is irrigated if psi was below psi_dry */
	MinEtr2Etp; // minimum quotient of ETR/ETP to start irrigation 
	int     tabcnt;  /* if irrigation control according to a irrigation table, this is the number of entries in this table */
	long   *month,   
	*day;     /* date of irrigation, month and day are arrays with dimension of "tabcnt" */
	double *amount;  /* array for the amounts of irrigation water for each irrigation termin */
};

struct tracer_info_type  /* strcture for a table, containing infos about tracer */
{
	char *name;      /* name of the tracer for identification in error messages */
	int radioactive; /* 0 if nonradioactive, 1 if radioactive */
	int evaporating; /* 0 if non evaporating (like salt), 1 if evaporating (like tritium) */
	double min_conc, /* minimum concentration (z.B. 0 for NaCl, -9999 for 18O, 0 for 3H) */
			max_conc, /* maximum concentration (physically possible) -> e.g. 0.03 for NaCl (in kg/kg) or "30 g/l" */
			half_time, /* half_time of radioactive tracers in hours -> Tritium: "107820 h" or "4493 d" */
			gwini,     /* initial tracer concentration within groundwater (is used only, if no concentration grids are read in) */
			soilini,   /* initial tracer concentration within soil (equal distribution down to groundwater; only valid, if no grids are read in) */
			decfact,   /* factor to multiply the old value to get the new value after radioactive decay (decfact=exp(-(InV/60.0/24.0)/(half_time*(1/ln(2)))))*/
			routdecfact; /* decay factor for routing module (timestep for sub-routing intervals is much shorter than usual timestep) */
	char *path,*statext,*gridext; /* parts to built filenames from */
	int   statcode,writegrid; /* output code for the tracer statistic file and the grid-writecode */
};

/* a table used for the communication between the discharge model and the 
unsaturated zone model when estimating infiltration and exfiltration */
struct riverflow_table_type 
{
	double flow_mm[10];      /* <- from routing  : the flow, entering the actual routing channel from upstream in mm */
	double flowm3s;          /* <- from routing  : same variable but in m3/s */
	double speed;            /* <- from routing  : velocity of the flow, to get depth: depth =  */
	double upstream_area;    /* <- from routing  : total sum of the upstream tributaries */
	double infilt_amount;    /* -> from unsatzon : total infitration amount in m3 */
	double infilt_fraction;  /* <-> from both    : fraction of infiltrating water on river flow in a channel */
	int    ROUT_index;       /* <- from routing  : in wich routing description is this subbasin included as target subbasin (e.g. as x in "TG x (..) from ..." */
	int    OL_index;         /* <- from routing  : what OL-index has this subbasin within the routing description, important for SUMTRIBs, since they have in around 50% of all cases no own channel */
};

struct  firstlast 
{
	int First;
	int Last;
};

struct buffertype
{
	long j,m,t,s,min;
	double *wert;
};

struct doublebuffertype
{
	long j,m,t,s,min;
	double *wert;
};

// for the moment, the filebuffer will be a read buffer only - later, write capabilities may be added

#define FB_STATE_UNINITIALIZED  0
#define FB_STATE_INITIALIZED    1
#define FB_STATE_FILEENDREACHED 2
#define FB_STATE_FILENOTFOUND   3
#define FB_STATE_LASTRECORDREAD 4
#define FB_STATE_FILENOTGOOD    5


class FileBuffer
{
private:
	std::fstream m_file;
	std::string m_filename;
	buffertype *m_buffer;
	long m_memsize;
	long m_bufpos;
	long m_bufmax;
	long m_valcnt;
	
	bool m_useMinutes;

	unsigned long m_filesize;
	unsigned long m_filepos;
	int m_mode; // 0=read, 1=write (1 is not supported in the actual version)
	int m_status; // 0=uninitialized, no filename given (default constructor)
	// 1=initialized successfully
	// 2=file end reached
	// 3=file not found
	// 4=last record already read
	char *stringbuffer1;
	char *stringbuffer2;

public:
	FileBuffer(); // default constructor: empty filename, mode read=0,memory size 1MB
	FileBuffer(std::string filename,int mode=0,long memorysize=1000000);
	~FileBuffer();
	void ClearMemory();
	long GetMemsize() {return m_memsize;};
	long GetBufpos() {return m_bufpos;};
	long GetBufmax() {return m_bufmax;};
	long GetValcnt() {return m_valcnt;};
	std::string GetFilename() {return m_filename;};
	unsigned long GetFilepos() {return m_filepos;};
	int GetMode() {return m_mode;};
	int GetStatus() {return m_status;};

	void reset(); // clear all buffer elements, reset filepos, reset bufpos,clear filename
	int initialize(std::string filename,int mode,long memorysize);

	void SetBufPos(long pos); // for free access to the buffer, if required
	void SetBufPosRel(long relpos); // for free access to the buffer, if required; relpos is relative to the actual position, can be positive or negative

	void NextRecord();
	void ReadBufferFromFile();
	void WriteBufferToFile();
	void GetActdate(int &year,int &month,int &day,int &hour,int &minute);
	void GetActValues(double *value);
	double GetActValue(int i);
};

class zonalstatistic
/* Routinen fuer die zonale Statistik */
{
public:
	int initialized;
	bool useMinutes;
	wasim::mystream *WSstream; 
	statrecord *zone;
	buffertype *buffer;
	TLastResults lastresults;
	char *titel;
	/* depending on buffersize (the maximum number of bytes, used for each buffer, in modconst.h)
	and on the number of zones a maximum indef for buffer is calculated at the allocation time */
	unsigned long bufmax;
	unsigned long bufcnt;          /* actual index for buffer (1...bufmax) */
	unsigned long wertcnt;
	double **RelCellArea; /* relative real area of a cell of cellsize*cellsize (usually smaller values possible due to effects of projections) */
	int nrows,ncols;
	double nodata;
	int zonmax;
	int prec;
	int code;
	int intervalle;
	int numsubb,*listsubb;
	char *dateiname;
	char *DefaultOutputDirectory;
	std::FILE *zonestat;
	double mittelwert;
	int cnt; /* zaehlt die Intervalle und wird immer dann auf null gesetzt, wenn es == intervalle ist */
	double gesamtsumme;
	double *sumn,*sumsum,*sumRelCellSum; /* zur Speicherung von Zwischensummen bei Mittelwertbildung fuer mehrere Intervalle */
	/* Routinen fuer die zonale Statistik */
	Tzarray    *zarray;
	virtual void fileopen(const char *mode);
	virtual void fileclose();
	virtual void flushbuffer();
	virtual void wertzufuegen(int i,int j,int index,double x,double weight=1.0);
	virtual void init();
	virtual void calc_and_done(long jahr,long monat,long tag,long stunde,long minute,int gather=1);
	virtual int schreibzonenwerte(const char *titel,const char *instanzname,TLastResults *LastResults,int maxLastResults);
	zonalstatistic(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double **p_RelCellArea,int p_nrows,int p_ncols,statrecord **p_zone,int p_zonmax,int p_code,
	int p_prec,char *p_dateiname,char *p_DefaultOutputDirectory,int p_numsubb,int *p_listsubb,bool useMinutes);
	virtual void ClearMemory();
};

class routstatistic:
public zonalstatistic
{
public:
	double *abfrel;
	TLastResults observedRunoff;
	/* Routinen fuer die zonale Statistik */
	virtual void Rsquare(Tstatist **p_eff,Tstatist **p_leff);
	virtual void calc_and_done(long jahr,long monat,long tag,long stunde,long minute);
	virtual int schreibzonenwerte(const char *titel,const char *instanzname,TLastResults *LastResults,int maxLastResults);
	routstatistic(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double **p_RelCellArea,int p_nrows,int p_ncols,statrecord **p_zone,double *p_abfrel,int p_zonmax,int p_code,int p_prec,
	char *p_dateiname,char *p_DefaultOutputDirectory,int p_numsubb,int *p_listsubb, bool useMinutes);
	virtual void ClearMemory();
};

class memory_handling
{
public:
	int initialized;
	wasim::mystream   *WSstream; // for output in a specified stream which is given to the constructor of a module
	char *instanzname;
	char* GetInstanzname() {return instanzname;};
	/* Routinen zur Belegung und Freigabe von Speicher fuer Vektoren und Matrizen */
	/*   virtual double *vector(long nl,long nh);
virtual void free_vector(double *v,long nl,long nh); */
	memory_handling():initialized(0){};

	virtual double  ***drei_d_matrix(long nrl, long nrh, long ncl, long nch, long levels, long nrows=0);
	virtual void   free_drei_d_matrix(double ***m, long nrl, long nrh, long ncl, long nch, long levels,long nrows=0);
	
//	virtual double **d_matrix(long nrl, long nrh, long ncl, long nch, long nrows = 0);
//	virtual void   free_d_matrix(double **m, long nrl, long nrh, long ncl, long nch, long nrows = 0);
	
	virtual double  **matrix(long nrl, long nrh, long ncl, long nch, long nrows = 0);
	virtual void   free_matrix(double **m, long nrl, long nrh, long ncl, long nch, long nrows = 0);
	
	virtual char   *newstring(int count,int errorcode);
	virtual void   speichererror(int fehler);
	virtual double log_2(double x);
};

class modelobject:
public memory_handling
{

	/* diese Klasse ist virtuel, d.h., sie wird nie benutzt, nur Ableitungen verwenden sie
Die hier definierten Typen und functionen sind auch in allen abgeleiteten Typen verfuegbar
abgeleietet werden z.B. die Classen der Interpolatoren und der Teilmodelle,
welche alle mit den Grunddaten des Zonengrids und der Gridbeschreibenden Daten arbeiten
sowie alle eine zonale Statistik berechnen und in das Ergebnisgrid schreiben */

public:
	/* nbumsubb: number of subbasins, for which an output should be written (if output code for a statistics file is 5xxx) */
	/* listsubb: contains as much entries as specified by numsubb, each containing a number of a subbasin, output should be generated */
	int iLogLevel; // 0=no output, 1=information level 2 = debug level
	// will be read from control file if paragraph [logging] is present as first parameter
	int   numsubb,*listsubb;
	double min_to_int_end, /* if the hour-entry in the date from the meoet-files is given for interval start, the min_to_int_end should be the minutes of a time step (like InV), otherwise 0 */
	hour_offset; /* converted value min_to_int_end into hours */

	int mon[13];  /* wird in init-routine initialisiert
				die meissten Instanzen werden eine Datums-Zaehlung eingebaut haben,
				deshalb wird das Array mit den Tagesanzahlen jedes Monats gleich in
				der Mutterklasse gespeichert */

	/* satz dient immer als Puffer zum Einlesen von ASCII-Daten */
	/*tsatz und ttstr sind Speicher fuer die ursprueglichen Adessen von Satz und tstr
			um beim Freigeben des Speichers nicht Fehler durch geaenderte Zeiger zu bekommen, da
			z.B. fgets beim Dateiende eine Null liefert fuer die Adresse des Satzes */
	char       *satz,*tsatz,*tstr,*ttstr;

	/* Alle Instanzen werden auf eine Steuerdatei zugriefen */
	char       *steuerdateiname;
	std::FILE       *steuerdatei;

	//LAURA set up request array for JustSend and JustRecv
#ifdef _MPI
	MPI_Request requestVectorLAURA[1000]; //LAURA
	MPI_Status statusVector[1000]; //Joerg
	int request_cnt_LAURA; //LAURA set up counter
#endif

	/* Die Standard-Statistikspeicher enthalten alle Instanzen */
	double statanz,statmini,statmaxi,sx,sx2,stat_nodata;

	/* sowohl fuer ein aktuell zu bearbeitendes Grid als auch fuer
				das Zonengrid koennen die Parameter eingelesen werden */
	int        ncols, nrows, mpistartrow, mpiendrow;
	double      xll,yll,yol,cellsize,nodata,fehlwert;
	double      **RelCellArea,**CellSizeX,**CellSizeY,**CoordX,**CoordY;
	int        cellsizexindex,cellsizeyindex,relcellareaindex,coordxindex,coordyindex;

	double		xllcorner_ext,yllcorner_ext,cellsize_ext; // when GIFs are read from external files, these values are taken from a meteorologic inpup file and must be set before reading the GIF from file with readgridfromfile
	double      **Mask; // reference to a mask grid that can be used to mask GIF-cuts
	double		gmx1, gmf1, gmx2, gmf2; // mapping factors for GIF-Value Mapping
	int			gmcut; // threshold for using gmx2 and gmf2 for mapping GIF Values (if GIF value is larger than gmcut, use *2, else *1)

	/* zum Einlesen der ersten 12 Datensaetze (beschreibende Daten) eines Grids */
	double      zz[13];

	/* fuer die zonale Statistik wird in zonalstatistik_init das Array
				zone mit <maxzonenzahl> Elementen eroeffnet */
	statrecord *zone;

	/* im Array zarray wird fuer jeden der hoechstens <maxzonenwert>
				grossen Codierungen des Zonengrids der Index im Array zone festgehalten */
	Tzarray    *zarray,*rarray;

	firstlast  *FirstLastCol;

	/* speichert die im aktuellen Zonengrid maximal vorliegende Anzahl an Zonen */
	int        zonmax;
	int        zonmaxwert;
	/* wenn mehrere Regionen verwendet werden sollen, dann muessen die Regionen
				wie uebergeordnete Zonen auch gezaehlt werden */
	double     transdist; /* transition distance in m used to calculate weights of neighboring regions which affect the regression results in the actual cell */
	int        regiomax;  /* how many regions? */
	int        regiomaxwert;  /* what is the largest code of the regions */
	center_of_mass_type *center_of_mass; /* a structure containing for each region the center of mass in real and in row/column coordinates */
	int        maxregionenzahl; /* a "constant" like for the zones, set to 32767 (so no region is allowed to  have a code <1 or >32767) */

	/* Standard-Datumsregister */
	long	minute,stunde,tag,monat,jahr,
			minutea,stundea,taga,monata,jahra,
			minutee,stundee,tage,monate,jahre;

	long       jf,mf,df,hf,minf,endtime; /* to compare dates assembled from year, month day and hour, a factor for each part is generated,
									so that the sum of the products jf*jahre + mf*monate + df*tage + hf*stundee is always comparable to
									other dates */

	MultiLayerLUSEType *MultiLayerLUTable;
	LUSEtabletype *LUtable;
	SOILtabletype *SOILtable;
	irrigationtabletype *irrigationtable;

	int  writegrid;
	bool WriteASCIIGrids; // when set to 0, Grids are written in binary format, otherwise in ESRI ASCII format
	bool WriteDoubleGrids; // when set to 0, Grids are written in double format, otherwise in double format, if not WriteASCIIGrids is set to values other than 0
	bool WriteNetcdfGrids; // when set to 0, Grids are written in double format, otherwise in double format, if not WriteASCIIGrids is set to values other than 0
	char *DefaultOutputDirectory, *DefaultInputDirectory;
	int timestep; // timestep in minutes, optional parameter, default will be set depending on meteorologic input files 



	/* wird aus der Steuerdatei gelesen. Ist writegrid > 0, dann werden die
						Daten aus dem Speicher double **gridinfo[i].pointer in eine Datei mit dem Namen
						char *gridinfo[i].name geschrieben. Der Name wird aus den ersten 4 buchstaben
						des gridinfo[i].name gebildet plus den aktuellen Daten fuer monat, tag und stunde
						z.B. gridinfo[i].name = regenidw.grd am 14.1. um 7 Uhr, -> schreiben in rege0114.07 */

	gridinfotype *gridinfo; /* ein array von records, in welchem jeweils der Name (char *name)
							des records (sollte einem eindeutigen Dateinamen
							entsprechen), eine Identification (z.B. "[Niederschlag]")
							sowie die Adresse des Speicherbereiches mit den Daten (double **pointer)
							gespeichert sind. Das Array (mit den pointern!) wird vom Haupprogramm
							initialisiert und dessen Adresse in der init-routine einer Instanz uebergeben.
							Diese belegt dann erst Speicher fuer die Namen und Grids */

	double **result; /* Ergebnisgrid im Memory (Adresse wird auf das richtige Element in gridinfo gesetzt,
					um den Zugriff zu beschleunigen */

	// MPI variables
	double **MPIREG,**MPIHLP; // MPIREG: coding the region to be used for each instance (with its rank); MPIHLP: contains a compressed code for fast addressing of the MPI zones
	int mpi_rank,mpi_size;
	long *v_mpiregioncnt;
	std::vector<mpihalo> *sendHalo,*recvHalo;

	int doSpinUp;
	std::vector<HRUMapElement*> *vHRU;
	double **HRU,**HRUMaster;

	int   resultnr; /* Speichert die Nummer im gridinfo-array, unter welcher das Grid direkt angesprochen werden kann,
					wenn nicht die indirekte Adressierung ueber **result genutzt wird (z.B., um den Namen zu
					erfragen) */
	int nachricht;  /* ist standardmaessig 1, dann wird zu jedem Intervall eine Nachricht ausgegeben
					andrenfalls nur dann, wenn die Stunde 24 ist (also bei Stundentakt einmal am Tag) */

	virtual void LogEntry(int iLevel,std::string sMessage);

	virtual void set_end(long y,long m,long d,long h,long min);                      /* sets the end-date (jahre,monate,tage,stundee) to the actual value */
	virtual void putnachrichtcode(int code);
	virtual int GetNachricht() {return nachricht;};
	virtual double GetHour_offset() {return hour_offset;};
	virtual int ist_vorhanden(const char *name,bool ignorepointer=false); /* liefert eine 1, wenn der uebergebene Parameter im
												gridinfo-array.identification gefunden werden konnte */
	/* Routinen zum Suchen von Grids im Arbeitsspeicher, die erste Variante sucht nach einem einzigen Grid,
			die zweite Variante sucht in der angegebenen Reihenfolge nach dem ersten Auftreten der Grids im
			gridinfo-array und liefert den index ins gridinfo-aray zurueck, unter dem das Grid angesprochen werden
			kann (wenn in der zweiten Variante das erste Grid nicht gefunden werden konnte, dann wird das zweite Grid
			gesucht, wenn diese auch nicht da ist, wird das dritte gesucht, erst wenn dieses ebenfalls nicht gefunden
			wurde, gibt es einen Fehler. Die erste Variante bricht mit einem Fehler ab, wenn das erste Grid nicht
			gefunden wurde*/
	virtual int getindex(const char *name);
	virtual int getindex_first_of_3(const char *name1,const char *name2,const char *name3);
	virtual int getindex_first_of_4(const char *name1,const char *name2,const char *name3,const char *name4);

	/* Standardfehlerbehandlungen */
	virtual void steuererror(const char *text,const char *message);
	virtual void writeerror(char *dateiname,int fehler);
	virtual void readerror(char *dateiname,int fehler);

	/* Leseroutinen fuer das Lesen von Daten aus der Steuerdsatei */
	virtual int  fscanfln(std::FILE *,const char *format,void *var,const char *errorstring,char *string);
	virtual int  fscanfln(std::FILE *,const char *format,void *var1,void *var2,const char *errorstring,char *string);
	virtual bool findeabschnitt(const char *Abschnittsname,const bool bErrorOnNoFound=true);


	virtual void clearstat(double nodata);
	virtual void addstat(double x,double weight);
	virtual void statistics(double *anzahl,double *mini,double *maxi,double *summe,double *mittel,double *stddev);
	virtual int  readstackfromfile(char *dateiname, double ****p_stack, bool MPIFullGrid = true);
	// To enable reading netcdfs instead of grid files, replace the file name by a netcdf and
	// pass the netcdf layer to this routine
	virtual void readgridfromfile(char *dateiname, double ***p_grid, int fillmiss,
								bool bAllocateGridMemory = _AllocateGridMemoryOn_, bool showmessages = _ShowMessagesOn_,
								double defaultValue = 0, bool limitAutoFill = _LimitAutoFillOn_,
								bool MPIFullGrid = _MPIFullGridOn_, bool MPIRank0Only = _MPIRank0OnlyOff_,
								bool MPISyncHalo = _MPISyncHaloOn_, bool isTimeSeries=false);
	virtual void printstackfile(char *dateiname, double ***stack, int levels, bool WriteASCIIGrids, bool bDoubleMode, bool WriteNetcdfGrids);
	virtual void printgridfile(char *dateiname, double **grid, bool bTextMode, bool bDoubleMode, bool WriteNetcdfGrids,int year=-1,int month = -1,int day = -1,int hour = -1,int minute = -1, double hour_offset=0);
	virtual void printsummengridfile(char *dateiname, double **grid, bool bTextMode, bool bDoubleMode, bool WriteNetcdfGrids);
	virtual void writestacktofile(int idx);
	virtual void writegridtofile(int idx, int gather = 1, bool doNotPrintAverageYet=false);
	virtual void initialize_grid(int i,const char *dateiname,double **grid,double ***stack,int Levels,const char *Bezeichnung,int writegrid,int pnrl,int pnrh);
	virtual void delete_grid(int i);

	/* lineare Interpolation */
	virtual double interpolate(double x1,double y1,double x2,double y2,double x);
	virtual void   SetInitialized(int p_initialized);

	/* Constructor, Destructor */
	modelobject(); // default constructor
	modelobject(wasim::mystream *Stream, char *steuername, int LogLevel, bool p_WriteASCIIGrids, bool p_WriteDoubleGrids, bool p_WriteNetcdfGrids, char *p_DefaultInputDirectory, char *p_DefaultOutputDirectory,
				long *pv_mpiregioncnt,std::vector<mpihalo> *pv_sendHalo,std::vector<mpihalo> *pv_recvHalo,int p_doSpinUp,std::vector<HRUMapElement*> *p_vHRU);
	virtual void ClearMemory();

	void MPIGatherGrid(int idx, double **grid = NULL);
	double** MPIGatherIntoNewGrid(int idx, double **grid = NULL);
	void MPIGatherStack(int idx, double ***stack = NULL);
	double*** MPIGatherIntoNewStack(int idx, double ***stack = NULL);

	void MPIScatterGrid(int idx,double **grid=NULL);
	void MPIScatterStack(int idx,double ***stack=NULL,int levels=0);
	void MPIGatherZone(statrecord *localzone, int zonecnt, bool gatherrelcellsum = true);
	void MPIScatterZone(statrecord *localzone, int zonecnt, bool scatterrelcellsum = true);
	void MPIGatherGridWithMask(int idx, int maskidx, long *CellsPerRank, double **OptionalTargetGrid = NULL);
	void MPIScatterGridWithMask(int idx,int maskidx,long *CellsPerRank);
	void MPISynchronizeGridHalo(int idx,double **grid = NULL);
	void MPISynchronizeStackHalo(int idx,int levels,double ***grid=NULL);
	void MPISynchronizeGridHaloJustRecv(int idx, double **grid = NULL); //LAURA
	void MPISynchronizeGridHaloJustSend(int idx, double **grid = NULL); //LAURA
	void MPISynchronizeGridHaloJustWait(int idx, double **grid = NULL); //LAURA

	void ExtendFileName(char* directory,char *dateiname);
	void SetExternalGridProperties(double xll,double yll,double cs,double **maskgrid);

#ifdef _NETCDF
    
protected:
	netcdfWriter *netcdfOut;
	netcdfReader *netcdfIn;

#endif

};

/* die Klasse defaultini reserviert Speicher fuer das Hoehenmodell
	und fuer das Zonengrid und gibt die Adressen dieser Arrays an den Aufrufer
	des Constructors zurueck (**dhm,**zn). Ausserdem untersucht es die Zonenfolge
	des Zonengrids und speichert diese im array zarray. Weiterhin wird das array
	*zone initialisiert, in welchem die statistikregister der zonalen Statistik fuer
	jede Zone zu finden sind (index des Registers ist in zarray fuer jeden zonenwert
	(index in zarray) gespeichert. Die Adressen von zarray und von zone werden
	ebenfalls an den Aufrufer zurueckgegeben. Wenn der Speicherplatz nicht
	ausreicht, sind die Adressen NULL */


class RiverLinkCell
{
public:
	int i; // original position of the cells row in the model grid
	int j; // original position of the cells column in the model grid
	int link; // river link ID, i.e. which river the cell belongs to
	int mpirank; // which rank does the cell belong to
	int processingrank; // which mpi rank is scheduled for processing that cell later
	double z; // elevation of the cell
	int processingorder; // defines the location of that cells buffer data. 1 is the first location, then numbers increase
	
	// default constructor
	RiverLinkCell() { i = 0; j = 0;  z = 0; link = 0; mpirank = 0; processingrank = 0; processingorder = 0; };
	
	// constructor with i, j, z, linm as parameters and optional mpirank and processingrank, although these values will and must be set later
	RiverLinkCell(int I, int J, double Z, int Link, int Mpirank = 0, int Processingrank = 0, int Processingorder=0) 
	{ 
		i = I; 
		j = J; 
		z = Z; 
		link = Link; 
		mpirank = Mpirank; 
		processingrank = Processingrank; 
		processingorder = Processingorder;
	};

	// copy constructor
	RiverLinkCell(const RiverLinkCell &c)
	{
		i = c.i;
		j = c.j;
		z = c.z;
		link = c.link;
		mpirank = c.mpirank;
		processingrank = c.processingrank;
		processingorder = c.processingorder;
	};

	// for sort function of a container class (like list) of these objects:
	// sort the river link cells by river link ID first, then by elevation
	bool operator<(const RiverLinkCell &t2) const
	{
		if (this->link < t2.link)
		{
			return true; // if link ID is smaller, then return true in any case
		}
		else
		{
			// if river link ID is equal, evaluate elevation (we are within the same river link)
			if (this->link == t2.link)
			{
				if (this->z > t2.z)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{ // link ID is greater than that of t2, so must return false here regardless of elevation
				return false;
			}
		}
	};

	bool operator>(const RiverLinkCell &t2) const
	{
		if (this->link > t2.link)
		{
			return true; // if link ID is smaller, then return true in any case
		}
		else
		{
			// if river link ID is equal, evaluate elevation (we are within the same river link)
			if (this->link == t2.link)
			{
				if (this->z < t2.z)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{ // link ID is greater than that of t2, so must return false here regardless of elevation
				return false;
			}
		}
	};

};

class RiverLink {
public:
	int linkID; // river link ID, i.e. which river the link belongs to
	int linkProcessingrank; // which mpi rank is scheduled for processing that links cells later
	int cellcnt; // numebr of cells that belong to that link
	std::vector<RiverLinkCell> *linkCellVector; // sub list for this link only, ordered from top to bottom, i.e. by elevation
	// default constructor
	RiverLink() :linkID(0), linkProcessingrank(0), cellcnt(0), linkCellVector(0)
	{
		linkCellVector = new std::vector<RiverLinkCell>;
	};

	RiverLink(int LinkID, int LinkProcessingrank, int Cellcnt, std::vector<RiverLinkCell> *LinkCellVector) :
		linkID(LinkID), linkProcessingrank(LinkProcessingrank), cellcnt(Cellcnt), linkCellVector(LinkCellVector)
	{
		if (LinkCellVector == NULL)
		{
			linkCellVector = new std::vector<RiverLinkCell>;
		}
	};

	~RiverLink()
	{
		linkCellVector->clear();
		delete linkCellVector;
	}

	virtual void setCellVector(std::vector<RiverLinkCell> *LinkCellVector)
	{
		// delete old linkCellList
		if (linkCellVector)
		{
			linkCellVector->clear();
			delete linkCellVector;
		}
		// set new linkCellList
		linkCellVector = LinkCellVector;
	};

	virtual std::vector<RiverLinkCell>* getCellVector() { return linkCellVector; };

	virtual std::vector<RiverLinkCell>::iterator getCellVectorLastElement() { return linkCellVector->end() - 1; }; // return iterator to the last element of the vector

	virtual void addRiverLinkCell(RiverLinkCell cell)
	{
		linkCellVector->push_back(cell);
		cellcnt++;
	};

	bool operator<(const RiverLink* t2) const
	{
		if (this->cellcnt > t2->cellcnt)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

	bool operator<=(const RiverLink* t2) const
	{
		if (this->cellcnt > t2->cellcnt)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

	bool operator>(const RiverLink* t2) const
	{
		if (this->cellcnt < t2->cellcnt)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

	bool operator>=(const RiverLink* t2) const
	{
		if (this->cellcnt < t2->cellcnt)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

};

class defaultini:
public modelobject
{
public:
	std::FILE *grid;
	char *gridname;
	double **zn,**nFK,**Na,**SB,**DHM,**SOIL,**IRRIG,**SLP,**ASP,**kWert,**ICE,**ALBEDO,**RSE; // some grids which are created or read in this instance rather than in any other instance (they will be used in multiple other instances usually like the zone grid or slope etc.)
	double ***LUSE,***RSC,***RSI,***Z0,***WURZ,***VBG,***LAI,***VEGSTART,***VEGSTOP,***VEGDURAT; // vectors of grids for multilayer Landuser -> this are not stacks!
	double*** ARCVEGFRACTIONS, *** AVBIOM, ***AVRETRANS;

	double hmax;
	int use;

	TLastResults *LastResults; // taking a pointer from an external calling program only (e.g. using wasim for windows)
	int maxLastResults;

	int naindex,kindex,nfkindex,soilindex,
	slopeindex,aspectindex,dhmindex,
	cellsizexread,cellsizeyread,relcellarearead;

	int rseindex;		// only needed once, not per layer - soil is there onyl once
	int albedoindex;	// only needed once, not per layer, albedo will be taken from the uppermost layer
	int *wurzindex;		// per vegetation layer
	int *vbgindex;		// per vegetation layer
	int *laiindex;		// per vegetation layer
	int *rscindex;		// per vegetation layer
	int *rsiindex;		// per vegetation layer
	int *z0index;		// per vegetation layer
	int *vegstartindex; // per vegetation layer
	int *vegstopindex;  // per vegetation layer
	int *vegduratindex; // per vegetation layer
	int *forcindex;     // per vegetation layer
	int *chillindex;    // per vegetation layer
	int *tstartindex;   // per vegetation layer
	int *fstarindex;    // per vegetation layer
	int *mtstart_index,mtstart_index1; // actual gridindex for MTSTART-stack and the index for the stack in the first vegetation layer
	int *arcvegfractionsindex; // per vegetation layer
	int* arcvegbiomassindex;
	int* arcvegretransindex;

	int tmovavg7index; // only needed once, not per layer
	int tavgindex;// only needed once, not per layer
	int startfallindex;// only needed once, not per layer
	int movavgcnt;// only needed once, not per layer

	int gridanzahl,do_expokorr,do_irrigation;
	int varGridsCount,varGridsStart,*varGridsRead;
	int numvegetationlayers;

	int regio1index,regio2index,regio3index,weight1index,weight2index,weight3index;
	double **REGION1,**REGION2,**REGION3,**WEIGHT1,**WEIGHT2,**WEIGHT3;

	int progaus_m;                
	int warning_level_m;          
	int runoff_measure_m;         
	double minutes_to_interval_end_m, MinutesToInvEndRouting_m;
	int timestep_m;
	int anzahlgroessen_m;         
	double breite_m;                 
	double laenge_m;                 
	double bezugsmeridian_m;         
	double timekorr_m;               
	char **meteonamen_m;             
	int numsubb_m;                
	int *listsubb_m;               

	// probably never really used
//	int ncols_m;   
//	int nrows_m;   
//	double xll_m;     
//	double yll_m;     
//	double cellsize_m;
//	double nodata_m;  

	int   version;
	char *steuername_m;

	// variables used for coupling with external models (e.g. groundwater models or external interpolations)
	int ExternalSyncCounter;
	int iDoExternalCoupling;
	char *SemaphoreFileNameRead;
	long lSemaphoreWaitInterval; // wait time between directory scans for the gwtableexternFilename
	char cCouplingMode; // 'I' for Interval, 'H' for hour, 'D' for day, 'M' for Month and 'Y' for yearly coupling
	int iExternalTimeStep; // time step used by the external model; important for converting e.g. groundwater level changes into fluxes in m/s

	int iNumExternalGridsRead; // number of external grids to read in for each coupling step
	char **ExternalGridNamesRead; // a vector holding the name for each of the external grids
	char **ExternalGridEntitiesRead; // a vector holding the wasim-internal entity name (e.g. "directflow" for the direct flow)
	double ***ExternalGridsRead; // vector holding references to internal grids
	int *iRenameMode; // vector of rename modes for each grid 0=delete external grid after reading in, 1=rename with new extension YYYYMMDD
	int *iFillMiss; // vector with fill-parameter for external grids: if 0, no missing values are filled, 1 when missing values (-9999) should be filled with nearest neighbors value

	int   iNumExternalGridsWrite;
	char  **ExternalGridNamesWrite; // a vector holding the name for each of the external grids to be written by wasim during synchronization
	char  **ExternalGridEntitiesWrite; // a vector holding the wasim-internal entity name for grid to write during syncgronisation (e.g. "directflow" for the direct flow)
	double ***ExternalGridsWrite; // vector holding references to internal grids
	double ***ExternalGridsWriteSum; // vector holding references to internal grids
	char  **ExternalGridsWriteCodes; // vector holding abbreviations on how the grids should be accumulated over time (Mean, Sum or Last)

	int iNumExternalStatsWrite;
	char **ExternalStatNamesWrite; // a vector holding the name for each of the external grids to be written by wasim during synchronization
	char **ExternalStatEntitiesWrite; // a vector holding the wasim-internal entity name for grid to write during syncgronisation (e.g. "directflow" for the direct flow)
	TLastResults *ExternalStatLastResults; // a vector holding references (will be initialized at the end of the initialization)
	int ExternalStatsMaxLastResults;
	char  **ExternalStatsWriteCodes; // vector holding abbreviations on how the grids should be accumulated over time (Mean, Sum or Last)

	char *SemaphoreFileNameWrite; // name of the file wasim will write after all writes for coupling are written
	char *SemaphoreContent; // content of the semaphoreWriteFile


	int ParamsFunctionState;
	/* Constructors */

	int substepcnt;
	int getSubStepCnt();

	// for use in a Windows-Environment when themodel is used as a DLL
	// in order to allow shorter subroutine calls. If This Constructor is used, 
	// the function Parameters() has to be called until the stateparameter is false
	// the parameters itself must then be retrived by callling GetParameters(...)
	defaultini(int &retvalue,wasim::mystream *Stream,char *steuernamearg,
		TLastResults *LastResults,int maxLastResults,
		int p_version,gridinfotype **p_gridinfo);

	int Parameters(bool *state);
	
	int GetParameters(char *steuername,int *LogLevel,
		int *progaus,int *warning_level,
		int *runoff_measure,double *minutes_to_interval_end, double *MinutesToInvEndRouting,int *timestep,
		int *anzahlgroessen,double *breite,double *laenge,
		double *bezugsmeridian,double *timekorr,
		char ***meteonamen,int *numsubb,int **listsubb,
		int *p_gridanzahl,
		int **p_zarray,  statrecord **p_zone, firstlast **p_FirstLastCol,
		int *p_zonmax,  int *p_zonmaxwert,center_of_mass_type **p_center_of_mass,
		int *p_regiomax,  int *p_regiomaxwert,int **p_rarray,
		int *p_ncols, int *p_nrows, int *p_mpistartrow, int *p_mpiendrow,
		double *p_xll ,  double *p_yll,
		double *p_cellsize  ,  double *p_nodata,
		long *mina,long *sta , long *tga , long *mta , long *jra,
		long *mine,long *ste , long *tge , long *mte , long *jre,
		LUSEtabletype **p_LUtable,SOILtabletype **p_SOILtable,
		irrigationtabletype **p_irrigationtable,
		MultiLayerLUSEType **p_MultiLayerLUTable,
		bool *WriteASCIIGrids,
		bool *WriteDoubleGrids,
		bool *WriteNetcdfGrids,
		char **DefaultInputDirectory, char **DefaultOutputDirectory,
		long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo);

	// this contructor does the same like the three prev. functions, but it takes 
	// longer until the control is givenback to the main program (usefull for text based applications)
	defaultini(int &retvalue,wasim::mystream *Stream,char *steuerdateinamearg,char *steuerdateiname,
		TLastResults *LastResults,int maxLastResults,
		int *LogLevel,
		int version,int *progaus,int *warning_level,
		int *runoff_measure,double *minutes_to_interval_end,int *timestep,
		int *anzahlgroessen,double *breite,double *laenge,
		double *bezugsmeridian,double *timekorr,
		char ***meteonamen,int *numsubb,int **listsubb,
		gridinfotype **p_gridinfo,int *p_gridanzahl,
		int **p_zarray,  statrecord **p_zone, firstlast **p_FirstLastCol,
		int *p_zonmax,  int *p_zonmaxwert,center_of_mass_type **p_center_of_mass,
		int *p_regiomax,  int *p_regiomaxwert,int **p_rarray,
		int *p_ncols, int *p_nrows, int *p_mpistartrow, int *p_mpiendrow,
		double *p_xll ,  double *p_yll,
		double *p_cellsize  ,  double *p_nodata,
		long *mina, long *sta , long *tga , long *mta , long *jra,
		long *mine, long *ste , long *tge , long *mte , long *jre,
		LUSEtabletype **p_LUtable,SOILtabletype **p_SOILtable,
		irrigationtabletype **p_irrigationtable,
		MultiLayerLUSEType **p_MultiLayerLUTable,
		bool *WriteASCIIGrids,
		bool *WriteDoubleGrids,
		bool *WriteNetcdfGrids,
		char **DefaultInputDirectory, char **DefaultOutputDirectory,
		long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo);

	/* Destructor */
	virtual void free_pointers(char *p1,char *p2,char *p3,char *p4,char *p5,char **p6);
	virtual void ClearMemory();

	virtual void replacesymbols(char *dateinameein,char *dateinameaus);

	void CreateStandardGridsStatistics();
	void InitGridsStatistics();
	void CalcAndDoneGridStatistics(long minute,long stunde,long tag,long monat,long jahr);
	void AddValuesGridStatistics();
	void FlushBuffersGridStatistics();
	void ClearMemoryGridStatistics();

	virtual void GetBasicControls(char *steuerdateiname,int *progaus,int *warning_level,
		int *runoff_measure,double *minutes_to_interval_end, double *MinutesToInvEndRouting, int *p_timestep,
		int *anzahlgroessen,double *breite,double *laenge,
		double *bezugsmeridian,double *timekorr,
		char ***meteonamen,int *numsubb,int **listsubb);

	virtual char *compilestring(char *satz,int zeile,char *symbol[],
	char *eintrag[],int anzahleintraege);


	virtual void irrigationtable_error(int i);
	virtual void write_Standard_Grids(long me,long se,long tg,long mt,long jr);
	virtual void write_Variable_Grids(long me,long se,long tg,long mt,long jr);
	virtual void GetMultiLayerLUSETable(int version);
	virtual void GetLuseTable(int version);
	virtual void GetSoilTable(int version);
	virtual void GetExpokorrIgnore();
	virtual void GetDates();
	virtual void GetGridNames();
	virtual void GetIrrigTable(int version);
	virtual void GetGridDimensions(int *p_ncols, int *p_nrows, double *p_xll, double *p_yll, double *p_cellsize, double *p_nodata);
//	virtual void readzonengrid(int *p_ncols, int *p_nrows, double *p_xll, double *p_yll, double *p_cellsize, double *p_nodata);
	virtual void zonenanalyse(gridinfotype **p_gridinfo, double **p_RelCellArea,
								int *p_ncols ,  int *p_nrows, double *p_nodata);
	virtual void  check_LU_Codes(int l);
	virtual void  check_SOIL_Codes(int l);

	// 2026-10-31 new methods for computing soil type from grain size distributions
	virtual int getUSDASoilTexture(double clay, double silt, double sand);
	virtual void  makeSoilGridFromClaySiltSand(int l);

	virtual void  check_irrig_Codes(int l);
	virtual void  readPeriodicStandardGrids(long minute,long stunde,long tag,long monat,long jahr,long Interval);
	virtual void  read_variable_Grids();
	virtual void  make_nFK_kWert_and_Na();
	virtual void  getnachbar(int i,int j,double **Grid);
	virtual void  make_Slope_and_Aspect();
	virtual void  check_ice_firn_grid();
	virtual void  MakeRegionWeightGrids();
	virtual void  SetRegionWeights();
	virtual void  analyze_regions(int regioindex);
	virtual void  CheckRelCellArea();
	virtual void  MakeCoordGrids();
	virtual double get_hmax();
	virtual double get_h_maximum(int i);
	// two versions of int_tvar_params are provided, the first one for Richards, the second one for Topmodel
	virtual void  int_tvar_params(long me,long st, long tg, long mt, long jr,int InV); /* fuer wasimuzr */
	virtual void  int_tvar_params(long me,long st, long tg, long mt, long jr,int InV,  /* fuer wasim */
		double **SBalt,double **SSP,
		double **nFKia,double **SBia);
	virtual void generateMultiLayerLanduse();
	virtual void generateMultiLayerParameterGrids(int numlayers);
	//  virtual void make_Wurz();
	// routines for external coupling
	virtual bool IsCoupled();
	virtual void GetCouplingParameters();
	virtual void InitializeFromExternalGrids(long jahr,long monat,long tag,long stunde);
	virtual char *GetSemaphoreFileNameRead();
	virtual char *GetSemaphoreFileNameWrite();
	virtual char *GetSemaphoreContent();
	virtual long GetSemaphoreWaitInterval();
	virtual char GetCouplingMode();
	virtual void InitializeExternalCoupling(TLastResults *p_LastResult,int p_maxLastResults);
	virtual void AccumulateExternalEntities();
	virtual void WriteExternalResults(int actCouplingHour,int actCouplingDay, int actCouplingMonth, int actCouplingYear,
		bool newHour,bool newDay,bool newMonth,bool newYear);
	//  virtual void WriteExternalResults();
	virtual bool ExistsExternalGWTable();
	virtual int  GetExternalTimeStep();

	// dynamic phenology grids and stacks
	// these grids are managed by defaultini but must be accessible within LUtable-objects
	// therefor the addresses of TSTART, FSTAR and MTSTART will be given to all LUTable entries 
	// after initialization
	double **T,**TAVG,**STARTFALL,**TMOVAVG0,***TMOVAVG7,**THROOT,***FORC,***TSTART,***CHILL,***FSTAR,****MTSTART;
	int iNumMTStartLevels; // used as info for LUtable-objects

	bool bUseForcAndTStart;
	bool bUseMTStart;
	void CalculateDynamicPhenology(int tagnr,int numlayers);

	double **MPIREG, **MPIHLP, **MPIHALO;
	int mpiregindex, mpihlpindex, mpihaloindex;
	void MakeMPIRegions(long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo, int *p_mpistartrow, int *p_mpiendrow);

	double **HRU,**HRUMaster; // result grids
	int hruindex,hrumasterindex;

	void MakeHRUs(int &doSpinUp,int &doTemperaturePreSpinUp,int &repeatCnt,std::vector<HRUMapElement*> **vHRU);

	std::list<RiverLinkCell> RiverLinkCellList;
	std::list<RiverLink*> RiverLinkList; // while creating the river link strukture for MPI-reinfiltration, the links need to be sorted several times, so we use a list
	std::vector<RiverLink*> RiverLinkVector; // after sorting, direct access is needed, so the list will be converted to a vector
	std::vector< std::vector<RiverLinkCell> > RivLnkCellsToSend;    // Array for each instance that needs to send cells to the processing rank. first dimension: remote processing rank, second dimension: single cells
	std::vector< std::vector<RiverLinkCell> > RivLnkCellsToReceive; // Array for each processing instance. First dimension: rank to receive from; second dimension: single cells in processing order
};

struct tracer_module_type
{
	statrecord      *zone,        /* initializing a statistc register set for groundwater level [m below surface] */
	**zone_stack;  /* if a module has more than onelayer (like layered groundwater module) */
	zonalstatistic  *stat,        /* output stream for tracer concentrations */
	**stat_stack;  /* for layered modules */
	char  *statfilename,          /* name of the statistic output file for unlayered modules */
	**statfilename_stack;    /* and for layered modules like groundwater module */
	char  *gridname,              /* name of the output grid for unlayered modules */
	**gridname_stack;        /* and for layered modules (like layered groundwater module) */
	int   storindex,fluxindex;    /* the index of the storage grid in the gridinfo-array for unlayered modules */
	int  *storindex_stack,
	*fluxindex_stack;        /* and for layered modules */
	/* grids for tracer concentrations, one for input, one for conc in storage, two other for conc. in fluxes */
	double **c_temp,
	**concstor,
	**cinput,
	**concflux,
	***conc_stack;    /* a stack of concentrations, if a module (like groundwater) has a layered structure */
};

class tracer_ini_class:
public modelobject
{
public:
	int tracer_count;
	tracer_info_type *tracer_info;

	tracer_ini_class(int &retvalue,wasim::mystream *Stream,char *steuername,int LogLevel,int *p_tracer_count,
		tracer_info_type **p_tracer_info,bool WriteASCIIGrids,bool WriteDoubleGrids,bool WriteNetcdfGrids,char *DefaultInputDirectory,char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,int p_doSpinUp,std::vector<HRUMapElement*> *p_vHRU);
	virtual void ClearMemory();
};
#endif
