/*   $Revision: 107 $:  $ 
	$Date: 7.03.11 0:36 $  */

/*

(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

header file, containing all deklarations for abstracted classes used in wasimuzr/wasim
other deklarations may be found in
modconst.h            (program constants, strings aso.
modappli.h/modapuzr.h (deklarations for main program variables)
modef.h               (deklarations for the origin classes and for utilities used by the abstracted classes in this file)
stringp.h, statist.h  (some global variables and function deklarations for string-functions and for statistics)


*/

#ifndef _TEILMODH_
#define _TEILMODH_

#include <iostream>
#include <stdio.h>
#include <map>
#include <queue>
#include "cspline.h"
#include "snow_canopy.h"

// used for statistics with respect to subbasins --> that's why there are bare and glacerized cells
struct glacierlisttype {
	bool isGlacierized;
	double CellCntGlacerized;
	double CellCntBare;
};

class elevationbelt {
public:
	double Hmin,
		  Hmax,
		  Hrange;
	long  NumMarginCells;
	std::vector<double*> marginCells; // stores direct access pointers to the margin cells in ICE-grids
	std::vector<double*> marginCellsIDs; // stores direct access pointers to the GLACID-grid
	std::vector<int> inactiveCells;
	void clearMarginCells();
	double getChangeableMarginCells(double change);
	double addChangesToMarginCells(double changePerCell);
	void addMarginCell(double *cellvalue,double *cellid);
	elevationbelt(double hmin,double hmax,double hrange):Hmin(hmin),Hmax(hmax),Hrange(hrange),NumMarginCells(0) {};
	bool isInRange(double elev);
};

// used for statistics on the separate glaciers and for V-A-relation updates 
// (there may be multiple glaciers in each basin, but the V-A-relation relies upon seperate glacier bodies)
class glacierclass {
public:
	// members
	int ID; // which glacier is it?
	double cellsize;
	long cnt; // number of glacierized cells
	double Area,      // area in km^2
		   Volume,    // volume in 1000000 m^3
		   oldArea,   // backup of old area
		   oldVolume; // backup of old volume 
	double ParB; // 28.5; in Volume = ParB*Area*exp(ParExp)
	double ParExp; // 1.36; in Volume = ParB*Area*exp(ParExp)
	double Hmin, // minimum elevation of the glacier tonge
		   Hmax; // maximum elevation of the glaciers accumulation area
	int numIter; // count of elevation belts for iterative grow/erode
	double totalChange; // in cells
	double remainingChanges; // in cells (for iterative changes)
    long  totalMarginCells; // count of margin cells
	int extraWeightFactorBand0; // if lowest elevation band should become more weight for growing/shrinking (recommended: 1 or 2)
	int subbasin; // which subbasin does the glacier belongs to. This will be checked for every cell when creating the statistic belts. If a glacier is located in more than one subbasins, an error will be thrown
	double BalanceCarryOver; // when glaciers are resized, the balance of the parts that are removed needs to be redistributed to all other remaining cells. 
	// Also, when the glacier grows, the balance of all glacier cells needs to be a adopted to the new amount of glacierized cells. The value of BalanceCarryOver is used to add up all these balance leftovers for redistribution to the remainig cells
	double TotalCellsBeforeResize; // weighed (be glacierized fraction) number of cells for the entire glacier before resizing
	double TotalCellsAfterResize; // weighed (be glacierized fraction) number of cells for the entire glacier after resizing
	double BalanceRedistribution; // amount to redistribute to all cells


	std::vector<elevationbelt*> elevBands;
	// methods
	glacierclass(); // default contructor 
	~glacierclass(); // default destructor 
	glacierclass(int id,double b,double exponent,double V,double A,int iter,double cs,int f=0,int subbasin=0); // constructor which sets Volume and area
	void addCell(double fraction,double elevation,double cellarea,bool isIce); // when creating the elevation belts, here we add the values
	void addBalanceChangeToOldVolume(double dV); // for initialization of oldVolume and oldArea
	void createElevationBelts(); // elevation belts must be created for each V-A-recalculation
	void resetMarginCounter();
	void addMarginCell(double *cell,double *id,double elevation);
	void calcVolume(); // doing the re-balancing for the entire Glacier
	void calcArea();
	void calcOldArea();
	void serialize(); // 
	double getArea();
	double getVolume();
	double getOldArea();
	double getOldVolume();
	int getExtraWeightFactorBand0();
	void setArea(double a);
	void setVolume(double v);
	void setOldArea(double a);
	void setOldVolume(double v);
	void setExtraWeightFactorBand0(int f);
	// get the amount in cellsize units which must be added or 
	// subtracted from the glacier
	double getTotalChangeInCells(); 
	void calcTotalChangeInCells(); 
	void addBalanceChangeToVolume(double dV);
	double changePerCell(double remainingChange,bool flag);
	double addChangesToMarginCells(double fraction,int iter); // takes the change per margin cell and returns the remaining changes (if there were not enough margin cells - another run must follow) 
	long getChangeableCellCount(double Change);

};

class meteointerpol:
public modelobject
{
public:

	int statfilecnt;

	int coordxrindex,coordyrindex;
	double **CoordXR,**CoordYR;

	long *stunde_a,*tag_a,*monat_a,*jahr_a,*minute_a;
	int   *ascii; /* is set to 0 (zero), if regress-parameter file is binaery, otherwise "ascii" is set to 1 */
	long  *regfilepos,*regfilepos_start,*regfilepos_preread; /* offset to the start of the file for regression files */
	long  *anetzfilepos,*anetzfilepos_start,*anetzfilepos_preread; /* offset to the start of the station data file if using methods 1, 3, 4 */

	double xmin, xmax, ymin, ymax;

	int   warning_level;
	buffertype **bufferanetz,**bufferregress;
	/* depending on buffersize (the maximum number of bytes, used for each buffer, in modconst.h)
and on the number of zones a maximum indef for buffer is calculated at the allocation time */
	long *bufmaxanetz,*bufmaxregress;
	long *bufcntanetz,*bufcntregress;          /* actual index for buffer row (1...bufmax) */
	int  *valcntanetz,*valcntregress;         /* index forthe buffer column */
	char  **regresssatz;
	int  *numAdditionalColumns;
	char	**anetzstat,**AdditionalColumns,**dateinamekoef,
	*ergebnisgridname,*zonestatname;
	int    zonestatcode;
	FILE   **anetz,**dateikoef,*ergebnisgrid;
	zonalstatistic *zonestat;
	/* Ergebnisgrid wird nur dann genutzt, wenn die Ergebnisse auch auf die
		Plate geschrieben werden sollen. Dazu muss der Parameter writegrids > 0 sein
		(in der Steuerdatei) */
	int    methode; /* methode, nach welcher interpoliert wird 1=idw, 2=regress, 3=idhw*/
	int    *anzahl,fehlwerti;
	bool *useMinutes; // if time step is shorter than one hour and/or the input file contains a mm columns, then thsi flag should be set tot true
	/* parameters used by all methods */
	double  grenzwertu,grenzwerto,fixu,fixo;

	// parameters for regressions done within WaSiM instead of in regress or regr
	double lowerinversionfix,upperinversionfix,tolerance,clusterlimit;
	int overlap;

	/* parameters used by distance weighting methods */
	double  maxentf,maxentfinv,wichtung,aniso,anisoneig,summedh,danteil; /* fuer idhw: wert = summe wert*(danteil*(entfgew)+(1-danteil)*(entfgew)) */

	/* private variables for OMP -> will be dclared private in a omp parallel directive or must be moved into the function calls therein */
	//  double *pe,*pt,*px,*py; 
	//  double  dx,dy,wert;
	//  double summeentf;

	/* parameters for regression methods */
	double *af,*bf1,*bf2,*igu,*igo,**a,**b;

	// parameters for lapse rate
	double lapserate,upperlapserate;
	char datatype[1000]; // P-type for data which are discontinuous (like precip or sunshine duration) or T-type for continuous data (like temperature)
	int iDataType;


	/* coordinate transformation identifiers/anisotropy parameters */
	double  alpha,beta; /* Drehwinkel bei Anisotropie */
	int     wi; /* Wichtung */

	/* identifiers for bilinera interpolation */
	//  double u,v,x[5],y[5],z[5],d[5],dmin,dsum;
	//  unsigned int stat[5];
	unsigned int ***STAT; /* Grid-array containing for each cell the four nearest stations as numbers of the input file */
	double *mindist,*direction,*slopeclass;
	double *dhmin,*dvmin,AvgHorDist,AvgVertDist;
	int    *minstation;

	/* identifiers for bicubic spline interpolation */
	TSplines *Spline;

	int *station;

	//  vector<int> StatArray;

	double **StationGrid,**StationIndexGrid,**SecondDerivative,**result2,**CoordXStatR,**CoordYStatR;
	int StatMeshNumCols,StatMeshNumRows;

	double **dhm,**zn;
	//double **RelCellArea; -> is already defined in modelobject or memory_handling
	/* Matrizen fuer das Hoehenmodell und fuer das Zonengrid
					beide werden als symbolische Adressen genutzt, denn die
					eigentlichen Daten stehen im Array gridinfo. um eine schnellere
					Adressierung zu gewaehrleisten, werden im Constructor
					jeweils die Adressen der entsprechenden Elemente des gridinfo-
					arrays an dhm und zn uebergeben. So kann ueber dhm und zn
					direkt das Array gridinfo angesprochen werden */
	double **xrot,**yrot; /* Matrizen, welche die x- bzw. y-Koordinaten des gedrehten
						DHM aufnehmen, wenn Anisotropie beruecksichtigt werden muss */
	double **value,**hoehe,**xkoord,**ykoord,**entf;
	double **xkoord_rot,**ykoord_rot;     /* Auch die Koordinaten der Stationen muessen
											bei Anisotropie gedreht werden, damit die
											Entfernungen zu den Stationen geometrisch stimmen */
	unsigned short int    *gleichewerte;
	unsigned short int    ende;
	char       *js,*ms,*ts,*ss,*mins;
	char       *meteosatz;         /* zum Lesen der Daten aus den EIngabedateien */
	char       *t_meteosatz;

	double kor; /* zur Korrektur der in das Grid zu schreibenden Daten
			mit einem bestimmten Faktor, z.B., um aus 1/10 C
			ganze Grade zu machen */

	double **szenario; /* enthaelt fuer jeden Monat einen Korrekturwert, der entweder addiert (methode 1), multipliziert (methode 2) oder
					prozentual addiert (methode 3) wird. Die Indizies sind so vergeben, dass der erste Index fuer
					den Rasterpunkt (entspricht einer Station im herkoemmlichen Sinn), der zweite fuer den Monat genutzt wird.
					die Indizies [*][1] und [*][1] sind die Koordinaten. Es wird eine Matrix (1,anzahl,1,15) gebildet.
					Im letzen Index wird jeweils der zeitlich interpolierte Szenariowert fuer eine Zelle gespeichert. */
	double **SCEN;     /* if monthly scenario grids in high resolution should be used for scaling the results -> like a weighting grid */
	char  **scengridname,*scengrididentifier; /* the 12 names of the scenario grids */
	int   scenmonth;  /* month to which the actual SCEN-grid refers */
	int   scenindex;

	int    szenmeth;  /* 0,1,2,3 je nach oben erklaerter Methode der ANwendung der Korrektur. 0=keine Korrektur, 4=Korrektur aus monatlichen Grids */
	int    szennum;   /* Anzahl an Szenariozellen */
	int    szenuse;   /* 0=keine Szenarien nutzen, 1=Szenarien nutzen */

	int regio1index,regio2index,regio3index,weight1index,weight2index,weight3index;
	double ***REGION,***WEIGHT;


	double dxmin,dymin,dist,delta,sumdhmin,sumdvmin,anzdhmin,anzdvmin;
	double xcrot,yrrot,csroty,csrotx; /* wie xc und yr, nur dass die Koordinaten die Drehung beruecksichtigen */

	char *gridlist_actname;
	char *gridlist_lastvalidname;

	double **scalinggrid;

	/* Constructor */
	meteointerpol(int &retcode,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,char *meteonamennummer,
		char *steuername,int LogLevel,int progaus,double p_min_to_int_end,int p_timestep,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,int p_resultnr,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,center_of_mass_type *p_center_of_mass,
		int p_regiomax,int p_regiomaxwert,int *p_rarray,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize  ,  double p_nodata,
		long mina, long sta , long tga , long mta , long jra,
		long mine, long ste , long tge , long mte , long jre,
		LUSEtabletype *p_LUtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp,std::vector<HRUMapElement*> *p_vHRU);
	/* dummy-version */
	//  meteointerpol(char *steuername);

	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
	virtual void set_warning_level(int warnlevel);
	virtual void set_hour_offset(double p_hour_offset);
	virtual void getdatetime(long *me,long *se,long *tg,long *mt,long *jr); /* holt die aktuelle Zeit von der Instanz */
	/* verzweigt je nach Methode in regressrasterizing, idw- oder idhwinterpolation */
	virtual void run(long vgl_minute,long vgl_stunde,long vgl_tag,long vgl_monat,long vgl_jahr);
	virtual void nextrecord(); /* verzweigt je nach Methode in nextrecordanetz und/oder nextrecordregress */
	virtual int  isend();  /* gibt 1 urueck, wenn das Dateieinde der EIngabedateien erreicht ist */
	virtual void set_end_date();  /* setzt das End-Datum auf den aktuellen Wert, d.h., das Modell soll anschliessend kontrolliert beendet werden */

	double *gethoehen()  { return hoehe[0];  }
	double *getxkoords() { return xkoord[0]; }
	double *getykoords() { return ykoord[0]; }
	double *getwert()    { return value[0]; }

	/* mit der folgenden Methode kann die Interpolation
	sowie die Berechnung des Korrekturwertes der Klimaszenarien nur fuer eine einzige Stelle
	anstatt fuer das gesamte Grid durchgefuehrt werden. Das wird genutzt, um bei der
	Niederschlagskorrektur die Daten von Wind und Temperatur auf die Standorte der Nieder-
	schlagsstationen zu interpolieren. Wenn der Niederschlag ganz oder teilweise durch
	hoehenabhaengige Regression berechnet wird, dann erfolgt die Korrektur fuer den
	hoehenabhaengigen Teil so, dass die Anstiege und/oder die Konstanten der
	Regressionsgeraden korrigiert werden. (Es wird fuer jeden Geradenabschnitt fuer jeweils
	10 Punkte der korrigierte Niederschlag berechnet. Aus diesen korrigierten Werten
	werden unter Beibehaltung der alten Inversionsgrenzen neue Geraden berechnet, wobei
	die mittlere Gerade nur als Verbindung der unteren und der oberen Geraden gezogen und
	nicht neu berechnet wird */
	double getsinglepoint(double x,double y,double z,double *entf);

	int getmethode() { return methode; }
	int getanzahl()  { return anzahl[0]; }

	/* liefert die Koeffizienten einer hoehenabhaengigen Regression an den Aufrufer zurueck */
	void getregresscoeffs(int nr,double *p_igu,double *p_igo,double p_a[4],double p_b[4]);
	/* schreibt die veraenderten Regressionskoeffizienten in die Felder der Instanz zurueck (P-Korrektur!) */
	void putregresscoeffs(int nr,double p_a[4],double p_b[4],double p_igo);

	virtual void flushbuffers();
	virtual void resetInputFiles();

private:

	virtual void readbufferanetz(int nr);      /* reading data from the anetz input file in an internal buffer */
	virtual void readbufferregress(int nr);    /* reading data from the regressinput file in an internal buffer */
	virtual void nextrecordanetz(int nr);  /* liesst naechsten Satz aus der Anetz-Datei */
	virtual void nextrecordgridlist();
	virtual void nextrecordregress(int nr); /* liesst naechsten Satz aus der Regress-Datei */
	virtual void parametereingabe(void); /* liesst die Parameter aus der STeuerdatei */
	virtual void stationswerte(FILE *statfile,char *statfilename,
								int numAddCol,bool *useMinutes,int *anzahl,double **hoehe,
								double **xkoord,double **ykoord,double **entf,
								double **value);
	//  virtual void stationswerte(void);    /* liesst aus der ANetz-Datei die Stationswerte */
	virtual void finddateanetz(FILE *statfile,char *statfilename,int numAddCol,bool useMinutes,int anzahl,
								double *value,long *statfilepos,long *statfilepos_preread,unsigned short *gw,
								long *me,long *se,long *tg,long *mt,long *jr);
	virtual void finddategridlist(FILE *statfile,char *statfilename,
								int numAddCol, bool useMinutes, char *gridname,long *statfilepos,long *statfilepos_preread,
								long *me,long *se,long *tg,long *mt,long *jr);
	//  virtual void finddateanetz(void);    /* sucht in der ANetz-Datei nach dem Anfangsdatum */
	virtual void finddateregress(int nr);  /* sucht in der Regress-Datei nach dem Anfangsdatum */
	virtual void GetMeanStationDist(int nr);
	virtual void ConstructGridFromStationCoordinates(int nr);
	virtual void MakeBilinearStationGrid(int nr); /* generates 4 grids containing for each cell the stations to use for bilinear interpolation (all station files must have the same stations) */
	virtual void GetBilinearSkewness();
	virtual void RotateCoordGrids();
	virtual void RotateStationCoordiantes();
	virtual int  GetNextStation(double x,double y,double dxmin,double dxmax,double dymin,double dymax);
	virtual int  StationsExist(int *p_stations, int p_anzahl);
	/* Routinen fuer die verschiedenen Interpolationsmethoden */

	virtual double getszenkorr(int i,int j,double x,double y);
	virtual double getidw(int nr,double x,double y,double *entf); /* liefert idw-Ergebnis fuer einen Punkt */
	virtual double getidw_anisotropisch(int nr,double x,double y,double *entf); /* liefert idw-Ergebnis fuer einen Punkt mit Beruecksichtigung der Anisotropie */
	virtual double getthiessen(int nr,double x,double y,double *entf); /* liefert idw-Ergebnis fuer einen Punkt */
	virtual double getThiessenWithLapseRate(double x,double y,double elevation,double lapserate,double upperinversionfix=-9999.0,double upperlapserate=-9999.0);
	virtual double getregress(int i,int j,double x,double y,double z); /* liefert regress-Ergebnis fuer einen Punkt */
	virtual double getbilinear(int nr,int i,int j,double x,double y,double *entf); /* liefert bilineare Interpolation fuer Stationen innerhalb eines festen Rasters-> Kopplung mit Klimamodellen, Methode 5 */
	virtual double determinante_drei(double *a,double *b,double *c);
public:
	virtual void  CalculateRegressionParameters(int nr);
	virtual void  WriteResultGrid();

};

class precicorrclass: /* eine Klasse zur Korrektur der Niederschlaege */
public modelobject
{
public:
	int do_it;
	meteointerpol **meteo;
//	double Tgrenz,ar,as,br,bs; /* Grenztemperatur Regen/Schnee, Anstiege und Konstanten der Korrekturfunktion */
	double *Tgrenz, *ar, *as, *br, *bs; /* Grenztemperatur Regen/Schnee, Anstiege und Konstanten der Korrekturfunktion */
	double *Phoehe;  /* Hoehen der Messstationen */
	double *Px,*Py;  /* x- und y-Koordinaten der Messstationen */
	int    Pmeth;                     /* Methode der P-Interpolation */
	int    Panz;
	int    regiomax;
	double  *igu,*igo,**a,**b,hmax;
	double *Twert,*Uwert,*Pwert;
	int    Pindex,Tindex,Uindex;
	double  **TEMP,**WIND,**PRECIP; // will be used in grid based precipitation correction only
	int    nowind,notemp;            /* wird als logische Variable abgefragt, um die Art der Korrektur zu steuern, wenn
									kein Wind und/oder keine Temperatur gemessen wurden */
	bool   doGriddedCorrection;
	
	bool   routingmode;

	// for regional correction, a weighted superposition of the correctionresults is required
	int regio1index, regio2index, regio3index, weight1index, weight2index, weight3index;
	double*** REGION, *** WEIGHT;

	precicorrclass(int &retvalue,wasim::mystream *p_Stream,meteointerpol **p_meteo, /* an array of pointers to dynamic instances of meteointerpol */
		char **p_meteonamen,char *steuername,
		int LogLevel,int anzahlgroessen,
		int p_regiomax, int p_regiomaxwert, int* p_rarray,
		gridinfotype *p_gridinfo,
		firstlast *p_FirstLastCol ,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll, double p_yll,
		double p_cellsize, double p_nodata,
		LUSEtabletype *p_LUtable,
		MultiLayerLUSEType *p_MultiLayerLUTable, bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,int p_doSpinUp,bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU);

	virtual void run();
	virtual void runOnSuperposedRegions();
	virtual void put_hmax(double p_hmax);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
};

// a class dealing with different interpolationj results wich are used for superpositioning
class teilmodellklasse:
public modelobject
{
public:

	//  double **RelCellArea;
	int tracer_count;
	bool routingmode; // if this flag is set, wasim runs for runoff routing only and all modulues will set do_it to 0 (except those needed for routing)
	tracer_info_type *tracer_info; /* contains infos for all tracers  */
	tracer_module_type *tracer_module; /* and a structure containing the variables for statistics (zonalstatistics, filenames etc.) */

	int    do_it;    /* wird als erster Parameter aus der Steuerdatei fuer jedes
			Teilmodell eingelesen und bestimmt die Abarbeitung des Modells
			0 = ignorieren, 1 = abarbeiten */
	int    gridslesen; /* 0=read grids not from disk but initialize them with a standard value; 1=read from disk */
	int    methode;  /* wenn mehrere Methoden angeboten werden, wird hierin die
			Nummer der zu nutzenden Methode gespeichert */
	int    fpeflag;
	double  InV,InVtag;   /* Intervallaenge in Minuten sowie Anzahl der Intervalle pro Tag*/
	int    substepcnt; // for unsatzonmodel, etpetr(real_uzr only), surface routing, and gw-model, the time step may be splitted

	double  *stack;   /* ein einfacher Stack (array) */
	double  **stackr; /* ein array (row of stacks) */
	int    stacktiefe;
	double  sunazimuthangle, sunaltitudeangle, sunzenithangle, sunazimuthangleAtNoon, sunaltitudeangleAtNoon;
	double sunweighingfactor; // used for weighing the temperature addition/subtraction depending on the relative sun elevation with respect to the altitude at noon (highest)
	double  cosincidence,breite,laenge,bezugsmeridian;

	double Anteil_hell,Anteil_dunkel, /* relative Anteile von Nacht und Tag am Zeitintervall */
	dark_a,dark_e;             /* Zeitdauer in Stunden, die es am Anfang und Ende des Intervalls dunkel ist */
	int   nostrl;
	double sunr,suns; /* Sonnenauf- und Sonnenuntergangszeiten in hh.hh (dezimal) */

	double pi,rad;

	double timeequation[14];

	/* gibt die Verschiebung der realen Sonnenzeit gegenueber der
		mittleren Sonnenzeit in min an (Elliptische Bahn) das Array steht
		an jeder Seite um eis ueber, d.h., es kann auch mit monat-1
		bzw. monat+1 gerechnet werden*/
	double  **zn; /* nimmt das zonengrid auf (bzw. Zeiger auf Array, welches das Zonengrid enthaelt */
	double  timekorr; /* Zeitverschiebung der Zeitzone gegebueber den Zeitangaben der gemessenen Daten,
					wenn die Daten in UTC angegeben werden, die Zeitzone aber MEZ ist, dann ist
					timekorr=1.0 zu setzen usw.  */

	int        *zzarray;   // temporaeres Indizierungsarray fuer die Zonananlyse 
	statrecord *zzone;     // temporaeres Zonanarray (statistikregister der einzelnen Zonen) fuer Zonenanalyse 

	virtual int doit();
	virtual int GetInV() {return (int)InV;};

	virtual int getSubStepCnt() {return substepcnt;};
	virtual void setSubStepCnt(int step) {substepcnt=step;};

	virtual void   fpe_check(int i,int j,int k,int z,double a1,double a2,double a3,double a4,double a5);
	virtual double deklination(long tag,long monat);
	virtual double gettimeequation(long tag,long monat);
	virtual double sonnenhoehe(long tag,long monat,double breite,double laenge,double stunde);
	virtual void   getsuncoordinates(long tag,long monat,double breite,double laenge,double stunde,
	double *sunazimuth,double *sunaltitude);
	virtual double getcosincidence(double aspect,double slope,double azimuth,double zenith);
	virtual void   sunrise(long tag,long monat,double breite,double laenge,double *sr,double *ss);
	virtual void   get_Tageslicht_Daten();

	/* substance_transport algorithms */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults); /* initializing substance_transport algorithms (getting tracer_info table adress) */
	virtual double decay(double conc, double decfact);               /* for calculating Ci = Ci-1 * decfact; decfact = exp(-(InV/60/24)/half_time);
																	returns the new concentration after dte decay */
	virtual double mixing(double c1,double c2, double v1,double v2); /* calc. mixing concentrations for volume1 (c1) with volume 2 (c2);
																	returns mixing concentration */
	virtual double evaporating_changes(double c,double v1,double ve); /* calc. and returning the concentration c after evaporating the amount ve from v */
	virtual void   substance_transport_delete(); /* free memeory */

	void zonenidentification(double ***zn,statrecord **zone,Tzarray **zarray,int &zonmax,int &zonmaxwert);
	
	/* Desctructor */
	~teilmodellklasse();
	/* Constructor */
	teilmodellklasse(wasim::mystream *p_Stream,
			double breite,double laenge,double bezugsmeridianarg,const char *bezeichnung,
			char *steuername,int LogLevel,int progaus,
			double p_min_to_int_end,int numsubb,int *listsubb,
			gridinfotype *p_gridinfo,
			int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
			int p_zonmax, int p_zonmaxwert,
			int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
			double p_xll , double p_yll,
			double p_cellsize , double p_nodata,
			long mine, long ste , long tge , long mte , long jre, double p_timekorr,
			LUSEtabletype *p_LUtable,
			SOILtabletype *p_SOILtable,
			MultiLayerLUSEType *p_MultiLayerLUTable,
			bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
			long *pv_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
			int p_doSpinUp,bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU,
			int psubstepcnt=1,bool SectionNotFoundRethrow=true);

	virtual void run(long me,long se,long tg,long mt,long jr,char *name, char* additionalInfo);

	virtual int getTimeStep(std::string section = std::string(""));

	/* es wird anstatt virtual ~teilmodellklasse(); ein Standarddestruktor genutzt */

};

// forward decrlaration of RegionalSuperposition because the members cSuperposition will use a poiter to its parent
class RegionalSuperposition;

class cSuperposition
{
	/* example for the control file structure to parse
temperature {
	entityinputgrid = temperature_reg ;
		regions = 1 ,   2,   3;
		weights = 0 , 0.5, 1.0;
	entityinputgrid = temperature_idw ;
		regions = 1 ,   2,   3;
		weights = 1.0 , 0.5, 0;
	outputgrid  = $outpath//$tempegrid ; 
		writecode = $Writegrid ;        
	outputtable = $outpath//temp//$grid//.//$code//$year ;
		statcode  = $hour_mean ;
}
*/
public:

	struct inputstruct {
		std::string sInputEntity; // name of the input grid as used as name in the section for meteointerpolation
		std::vector<int> region;
		std::vector<double> weight;
		double **INGRID;
		int ingridindex;
	};

	int jahr,monat,tag,stunde,minute;

	int ncols, nrows, mpistartrow, mpiendrow;
	double xll,yll,cellsize,nodata;
	Tzarray    *zarray,*rarray;
	firstlast  *FirstLastCol;
	int regiomax;
	int regiomaxwert;
	int regio1index,regio2index,regio3index,weight1index,weight2index,weight3index;
	double ***REGION,***WEIGHT,**zn;

	StringUtils su;

	std::string sEntityName; // how will this entity refered to in the model, e.g. temperature 

	std::vector<inputstruct*> inputs;

	RegionalSuperposition *parent;
	std::string initials;

	double **RESULT;
	int gridwritecode,resultindex;
	std::string resultgridname;

	// variables for zonal statistics
	statrecord *stats;
	zonalstatistic *zonestat;
	int  statwritecode;
	std::string statfilename;
	char *cStatfilename;

	double **MPIREG,**MPIHLP;
	int mpi_rank,mpi_size;

	cSuperposition(RegionalSuperposition *p_parent,std::string p_initials);
	~cSuperposition();
	virtual void WriteResultGrids();
	virtual void flushbuffers();
	virtual void run(long me, long se,long tg,long mt,long jr);
	virtual void ClearMemory();

};

class RegionalSuperposition:
public teilmodellklasse
{
public:

	int iNumberEntities;
	TLastResults *LastResults;
	int maxLastResults;
	int regiomax;
	int regiomaxwert;
	int *rarray;

	std::vector<cSuperposition*> vSuperpos;

	RegionalSuperposition(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_regiomax,int p_regiomaxwert,int *p_rarray,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre,double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void WriteResultGrids();
	virtual void flushbuffers();
	virtual void run(long me, long se,long tg,long mt,long jr);
	virtual void ClearMemory();
};

class expokorrklasse:
public teilmodellklasse
{
	/* in der Instanz dieser Klasse werden je nach  */
public:

	double hmax;
	double gewicht,sp,spa,spe,strlsumme;
	double cossunazim;
	int   korrparameter; /* 1 nur Strahlung, 2 Strahlung und Temperatur */
	double **EXKORR, /* Ergebis der Expositionsbedingten Korrektur */
	**TKORR,  /* korrigierte Temperaturen */
	**ASP,    /* Aspectgrid */
	**SL,     /* Slopegrid */
	** SSDABS,    /* Sonnenscheindauer-Grid absolute values in minutes per time step */
	** SSD, /* Sonnenscheindauer-Grid relative values*/
	**dhm,    /* Hoehenmodell */
	**SHAPE,  /* Schattengrid 0=Sonne, 1=Schatten */
	**T,      /* Temperatur (Tag bei InV > 720, sonst alleinige T für ges. Intervall) */
	**Tnacht, /* Temperatur Nacht bei InV > 720 */
	**IEFFMAX;/* theoretic maximum radiation regarding also the local slope, i.e. the Angle of incidence between the normal to the slope and the sun angle */

	int ssdIsAbsolute;

	double Isenk,spreiz;
	int   *indexx,*indexy; /* zwei Zeiger auf ein Array aus Integern  */
	//  double cosTheta;
	int   dx,dy,vx,vy; /* fuer die Geradengleichungen */
	int optim;
	double dz;
	double  *cos_za;  /* Ein array, das fuer alle Teilintervalle eines Intervalles InV > 60 min den sonnen-zenithwinkel speichert */
	int   EXKORRwritegrid,SHAPEwritegrid,TKORRwritegrid,ieffmaxwritegrid;
	int   exkorrindex,shapeindex,tnachtindex,tindex,tkorrindex,ieffmaxindex,ssdindex,ssdabsindex;
	int   cn,count;    /* nach je count Intervallen werden neue Werte berechnet */
	int   split;    /* es werden jeweils split Teilintervalle gebildet. Fuer jedes dieser
					Teilintervalle wird die Mittelzeit genommen, fuer diese Mittelzeit wird dann
					die Expositionskorrektur gerechnet. Der Mittelwert aus allen Berechnungen
					wird als repraesentativ fuer das Intervall angenommen. Es werden nur Intervalle
					beruecksichtigt, an denen die Sonne ueber dem Horizont steht.
					Bei einem Zeitintervall von 24 Stunden kann split z.B. 24 sein, es werden dann
					nach jeweils count Rechenaufrufen (also count Tagen) fuer alle 24 Stunden
					die Korrekturwerte ermittelt un dann ain repraesentetiver Tageswert errechnet. */
	double split_temp; /* wenn das Intervall nicht vollstaendig besonnt ist, wird split_temp so berechnet,
					dass eine minimale Anzahl von hoechstens der Laenge eines normalen Teilintervalls
					gebildet werden (also wenn z.B. split = 4 ist und die Sonne nur in den
					letzten 40% des Intervalls scheint, werden 2 Teilintervalle a 20 min gebildet */

	char  *EXKORRgridname,
	*SHAPEgridname,
	*TKORRgridname;

	int sr, er; // startrow and endrow for MPI applications

	expokorrklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre,double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void run(long me, long se,long tg,long mt,long jr);
	virtual double **GetExpoKorrAdress();
	virtual void convertAbsoluteToRelativeSSD(int minute,int hour,int day, int month,int houroffset);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();

};

class verdunstungsmodellklasse:
public teilmodellklasse
{
public:
	/* some constants for calculating evaporation (initialized in constructor, see there for explanation */
	double sealev1013   ,   /* = 1013.0    air pressure at see level */
	c7991        ,   /* = 7991.0;   constant for calculating pressure from elevation */
	c29_33       ,   /* = 29.33;    constant for calculating pressure from elevation */
	c273_15      ,   /* = 273.15;   absolute temperature at 0 degree celsius */
	c0_48        ,   /* = 0.48;     1. constant for calulating longwave radiation */
	c0_065       ,   /* = 0.065;    2. constant for calulating longwave radiation */
	c694_67      ,   /* = 694.67;   1. factor for konverting Wh/m2 to mm/Interval */
	c0_6589      ,   /* = 0.6589;   2. factor for konverting Wh/m2 to mm/Interval */
	c0_9         ,   /* = 0.90;     reduciong net radiation by accounting for soil heat flux */
	mag237_3     ,   /* = 237.3;    factor from the Magnus-formula for saturation vapour pressure  */
	mag25029     ,   /* = 25029.0;  derivation factor of magnus formula  */
	mag17_26     ,   /* = 17.260;   factor from Magnus formula */
	latheat2500_8,   /* = 2500.8;   latenbt heat of water */
	psychro1_6158,   /* = 1.6158;   1. factor for psychrometric constant */
	latheat2_372 ,   /* = 2.372;    2. factor for psychrometric constant */
	c64_0        ,   /* = 64.0;     */
	c1_0         ,   /* = 1.0;      */
	c0_54        ,   /* = 0.54;     */
	c20_0        ,   /* = 20.0;     */
	c4_72        ,   /* = 4.72;     */
	c0_0004      ,   /* = 0.0004;   */
	c0_004       ,   /* = 0.004;    */
	c0_6         ,   /* = 0.6;      */
	c0_0067      ,   /* = 0.0067;   */
	c0_216876    ,   /* = 0.216876; */
	c60_0        ;   /* = 60.0;     */

	int iSnowModelMethod;
	int mcode,netto;
	int m_iInterceptionIsActive;
	int lake_active,lake_ini_method;
	double rsc_OpenWater,z0_OpenWater,LAI_OpenWater,Albedo_OpenWater,VCF_OpenWater; // default parameter values for evapotranspiration for open water (lakes)
	int p_T,p_E,p_RLF,p_U,p_SSD,p_RG,p_LUSE,p_LDR; /* zur Berechnung eines Bit-Codes fuer
						die einzelnen Methoden */
	double  **Albedo,**dhm,**RLF,**E,**LDR,**RG,**SSD,**U,**SSDTemp, //(ggf. modellierte SSD) in SSDTemp
	**T,      /* Temperatur (Tag bei InV > 720, sonst alleinige T für ges. Intervall) */
	**Tnacht, /* Temperatur Nacht bei InV > 720 */
	**NFK,**SB,**SSNOW,**EXKORR,**SSNOWGL,
	**CL; // cloudiness, used for new snow model extensions by Warscher
	double  **EI; // alias for the gridinfo-element holding EI[0] (sum of all interzeption evaporation) for calculating the total real evapotranspiration

	// grids for layer specific entities
	double  ***ETP,***ETR,***EIP,***LUSE,***WURZ,***LAI,***VBG,***RSC,***RSI,***Z0; 
	double **EVAP,**EVAR,**EVRS,**EVRSGL,**RGEX,**RSE; // potential and real evaporation from soil (only one grid is needed because there is one evaporating soil layer only); EVRS: Snow-evaporation
	double **LAKE,**POND; // used for lake model: if a cell is within a lake, and the POND-value is positive, the potential evaporation will be calculated after Penman (without monteith)
	double **ICE; // for statistics: ice covered fraction has to be regarded for ETR, EI and EVAR (but not for EVRS)
	int iNumVegetationLayers;
	bool firsttime;

	int do_ArcVeg;
	int arcVegNumLayers;
	double*** ARCVEGFRACTIONS; // used for keepting track of the fraction of each landuse type when ArcVeg is used: there may be multiple (usually 5) land use types per cell, these grids keep track of their fraction (created in defaultini)

//	// stacks for layered snow model, declared in unsatzonmodellklasse. Will be used for subtracting snow evaporation for simple snow melt methods
//	double ***refSWES;
//	double ***refSWEL;

	/* Memorygrids zur Aufnahme der Variablen, bis auf
	**ETP und **ETR sind alle schon belegt, **ETP und **ETR sind
	die Ergebnisse dieser Instanz */
	double nfkratio,tagnr;
	double *PEC,fk,*FH;    /* empirische Faktoren fuer HAMOIN, WENDLING und HAUDE-Verdunstung */
	//  double deltaT,    /* Temperaturamplitude zwischen an e gewichteten Tages- und
	//                      Nachtmittel der Temperatur, Ergebnis einer e-Funktion
	//                      f(dhm,tag im Jahr,relRG) */
	double dtmax_akt, /* aus Monatswerten auf den Tag interpolierter Wert der max. T-Differenz der Tag-Nacht-T-Mittel als f(tag im Jahr)  */
	dtmax[14], /* Monatswerte der max. T-Differenz (Tag-Nacht-T-Mittel) fuer 0 m Meereshoehe */
	dt_koeff,  /* Abklingkonstante der e-Funktion zur Ber. der Hoehenabh. der T-Diff. in m */
	Ampl_tag,
	Ampl_nacht,/* Aufteilugn der Amplitude der Temperatur zwischen Tag und Nacht (gewoehnlich Tags etwas mehr -> Ampl_tag = 0.63 */
	Ampl_korr, /* der Amplitudenaufteilungsfaktor hat einen Jahresgang,
					(im Winter Ampl_tag groesser, bis 0.75, im Sommer dagegen kleiner, bis 0.45, also eine Spanne von 0.3,
					Ampl_korr wird fuer jeden Tag aus Ampl_korr = Ampl_Jg *cos(tagnr/265.0*2.0*Pi) berechnet und dann zu Ampl_tag
					addiert und von Ampl_nacht subtrahiert)
					*/
	Ampl_Jg,   /* Jahresgang-maxima der Abweichung der Temp.Amplitudenaufteilung zwischen Tag/Nacht-Temperatur
					bezogen auf den Jahresmittelwert (also z.B. 0.15, wenn die mittlere Aufteilung der T-Amplitude
					zwischen 0.75 im Sommer und 0.45 im Winter schwanken kann */
	minrelRG;  /* Ein minimalwert für die relative Strahlung, unterhalb dessen SSD=0 angenommen wird.
					Dieser Parameter wird als Nullstelle aus dem Polynom 3. Ordnung für die Modellierung der
					Sonnenscheindauer aus der relativen Strahlung nach dem Sekantenverfahren berechnet
					und muss zwischen 0.1 und 0.5 liegen, ansonsten wird ein Fehler vermutet. */
	
	int		ssdCnt; // counter of intervals in light day in order to calculate mean values for the following night
	int     ssdFlag; // will be set to 0 or 1 according to the next action required (0=no action, 1=to calculate mean value)

	// will be stored in SSDTemp now (because avoiding synchronisation between tasks in MPI-Versions
	//  statrecord *ssdstat; // fuehrt fuer alle Zeitintervalle mit Anteil_hell > 0.05 eine teilgebietsweise
	// Statistik der SSD mit, um in anschliessenden dunklen Intervallen die SSD abschaetzen
	// zu koennen (fuer langwellige Ausstrahlung in der Nacht)

	statrecord *zoneRGEX; /* statistikregistersatz fuer die korrigierte Globalstrahlung */
	statrecord *zoneCL; // statistics for Cloudiness
	statrecord **zoneETP;
	statrecord **zoneEIP;
	statrecord **zoneETR;
	statrecord *zoneEVAR; // real and potential avaopration from soil is handled seperately 
	statrecord *zoneEVAP; // real and potential avaopration from soil is handled seperately 
	statrecord *zoneEVRS; // real avaopration from snow is handled seperately 

	double *cos_za,split;
	double CS[4],cr[4],extink,phi,dphi;

	double sigma; /* Boltzmann-Konstante in Wh/m^2) */

	double ldr,rlf,e; /* je nach Methode werden diese Variablen erst errechnet oder
				direkt aus den Eingabegrids entnommen */
	double d_dhm[5];  /* fuer die aus der LUtable interpolierten Werte */
	int Nadelwald; /* bekommt den Code, den Nadelwald hat (in LUtabelle in Steuerdatei codiert),
					da die Widerstaende von Nadelwald speziell behandelt werden muss */

	double globalETPHamonFactor;

	char           **ETPstatfilename, **EIPstatfilename, **ETRstatfilename, *RGEXstatfilename, *EVARstatfilename, *EVAPstatfilename, *EVRSstatfilename, *CLstatfilename;
	zonalstatistic **ETPstat,         **EIPstat,         **ETRstat,         *RGEXstat,         *EVARstat,         *EVAPstat,         *EVRSstat,         *CLstat;

	char           **ETPgridname,     **EIPgridname,     **ETRgridname,     *RGEXgridname,     *EVARgridname,     *EVAPgridname,     *EVRSgridname,     *CLgridname;

	int    *etpindex,  *eipindex, *etrindex,  sbindex, evarindex,evapindex,evrsindex,rgexindex,ssdtempindex,clindex;/* um beim Freigeben des Speichers den Index in Gridinfo zu kennen */
	int    ETPwritegrid, EIPwritegrid, ETRwritegrid, EVAPwritegrid, EVARwritegrid,EVRSwritegrid,RGEXwritegrid,CLwritegrid;
	int    ETPstatcode,  EIPstatcode, ETRstatcode,  EVAPstatcode,  EVARstatcode,EVRSstatcode, RGEXstatcode,CLstatcode;

	verdunstungsmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre,double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void setSubStepCnt(int step); // inherited from teilmodellklasse and changed here

	virtual void flushbuffers();
	virtual int nullstellen(double a,double b,double *n1,double *n2);
	virtual double Esatt(double t);
	virtual double getETP(int layer,int i,int j,double tm,double Anteil,int tags,double &evap,double &eip,double &ra_local,double &Rrest_night,double &Rrest_day,double &ssd,double &RK,double &LAIkum);
	//  virtual double getETP(int layer,int i,int j,double tm,double Anteil,int tags);
	virtual void find_SB_and_NFK();
	virtual void SetSnowStorAdress(double **adresse);
	virtual void SetExpoKorrAdress(double **adresse);
	virtual void potentiell(long me, long se,long tg,long mt,long jr);
	virtual void snow(long me, long se,long tg,long mt,long jr);
	virtual void real(long me, long se,long tg,long mt,long jr);

	virtual void pre_real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void post_real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep);
	
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();

	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	/* initializing substance_transport algorithms (getting tracer_info table adress) */
	virtual void substance_transport_delete(); /* free memeory */
	virtual void initialize_lakes(std::string instanzname);
	virtual void setSnowModelMethod(int sMethod) {this->iSnowModelMethod = sMethod;};
};

class schneemodellklasse:
public teilmodellklasse
{
public:

	zonalstatistic *qschnee,*qsnow_out,*snowdays,*snow_age,*p_snow,*p_rain,*snowstor,*snowstor2,*albedo;
	statrecord *zoneQSNOW,*zoneQSNOW_OUT,*zoneSSNOW,*zoneSSNOW2,*zoneSNOWDAYS,*zoneSNOW_AGE,*zoneP_SNOW,*zoneP_RAIN,*zoneALBEDO;
	int  qschneestatcode,qsnow_outstatcode,snowstorstatcode,snowstor2statcode,snowdaysstatcode,snow_agestatcode,p_snowstatcode,p_rainstatcode,albedostatcode;
	char *qschneegridname,*qschneestatname,*qsnow_outgridname,*qsnow_outstatname,*snowstorstatname,*snowstorstatname2,
	*SSNOWgridname,*SLIQgridname,*SSTOgridname,*SNOWDAYSgridname,*SNOW_AGEgridname,*snowdaysstatname,*snow_agestatname,*albedostatname,
	*P_SNOWgridname,*p_snowstatname,*P_RAINgridname,*p_rainstatname;


	// new snow model extensions ================================================
	zonalstatistic 	*snowtemp,		*surftemp,		*ebal_snow,		*radflux,		*sensflux,		*latflux,		*adflux,		*soilflux,		*input_mass,		*mobile_mass,		*deposition;
	statrecord 	*zoneSNOWTEMP,	*zoneSURFTEMP,	*zoneEBAL_SNOW,	*zoneRADFLUX,	*zoneSENSFLUX,	*zoneLATFLUX,	*zoneADFLUX,	*zoneSOILFLUX,	*zoneINPUT_MASS,	*zoneMOBILE_MASS,	*zoneDEPOSITION;
	int	snowtempstatcode,surftempstatcode,ebal_snowstatcode,radfluxstatcode,sensfluxstatcode,latfluxstatcode,adfluxstatcode,soilfluxstatcode,input_massstatcode,mobile_massstatcode,depositionstatcode;
    char *SNOWTEMPgridname,*SURFTEMPgridname,
	  *INPUT_MASSgridname,*MOBILE_MASSgridname,*DEPOSITIONgridname,
	  *snowtempstatname, *surftempstatname, *ebal_snowstatname, *radfluxstatname, *sensfluxstatname, *latfluxstatname, *adfluxstatname, *soilfluxstatname, *input_massstatname,
	  *mobile_massstatname, *depositionstatname;
	double **RLF, **CL, **SNOWTEMP, **SURFTEMP, **SURFTEMPINTERPOL,	**RGEX, **SL,
		**SLIDEFRACTION1,	**SLIDEFRACTION2,	**SLIDEFRACTION3,	**SLIDEFRACTION4,	**DEPOINDEX,
		**ELEV_ORDER_ROW,	**ELEV_ORDER_COL,	**INPUT_MASS,	**MOBILE_MASS,	**DEPOSITION, **DEPOSITION_MAX,
		**EVRS; // snow evaporation
	int snowtempindex,	surftempindex,	input_massindex,	mobile_massindex,	depositionindex, depositionmaxindex, elevorderrowindex,elevordercolindex;
	int SNOWTEMPwritegrid,SURFTEMPwritegrid, INPUT_MASSwritegrid,MOBILE_MASSwritegrid, DEPOSITIONwritegrid;
	double maxDepositionSlope; // maximum slope for a cell in order to allwo deposition of moving snow masses (avalanches)
	double maxDepositionFactor; // factor for scaling maximum deposition amount from avalanches
	double minSlideSlope; // minimum slope that may lead to slides
	double slideMassFrac; // part of the snow pack that will be moved away in a slide (if steep enough)
	double LWINcorr; // factor for adjusting energy balance by reducing/increasing incoming long wave radiation
	double LWOUTcorr; // factor for adjusting energy balance by reducing/increasing incoming long wave radiation
	// end of new snow model extensions ================================================

	double  T0R,T0,CWH,CRFR,actC0,actC1,actC2,z0,RMFMIN,RMFMAX,SCA,ALBsnow_min,ALBsnow_max;
	double *C0, *C1, *C2; // monthly values for melt coefficients
	/* Grids from the gridinfo-array */
	double  **DHM,
	**P,
	**T,      /* Temperatur (Tag bei InV > 720, sonst alleinige T für ges. Intervall) */
	**Tnacht, /* Temperatur Nacht bei InV > 720 */
	**U,**E,**QSNOW,**QSNOW_OUT,RMFakt,
	**SSNOW,**SNOWDAYS,**SNOW_AGE,**P_SNOW,**P_RAIN,**SLIQ,**SSTO,**Albedo,**EXKORR;
	double ***LUSE; // LUSE in snowmodel will be the uppermost landuse layer (is most important for albedo, anyway)
	double **LAKE,**POND;
	int lake_active;
	int iNumVegetationLayers;
	double **ETRS;
	int etrsindex;
	bool ETRSchecked;

	// extension december 2014 for correct mass balaning
	double **SSNOWGL		,**SLIQGL		 ,**SSTOGL		  ,**QSNOWGL		,**SNOWTEMPGL		  ,**SURFTEMPGL		    ,**EVRSGL; // the three snow storage grids for glacierized fractions (must be read in when using readgrids=1) and surface temperature and snow temperature as well as snow evaporation
	int  ssnowglindex		,sliqglindex	 , sstoglindex	  ,qsnowglindex		,snowtempglindex	  ,surftempglindex      ,evrsglindex ;
	int  SSNOWGLwritegrid	,SLIQGLwritegrid ,SSTOGLwritegrid ,QSNOWGLwritegrid ,SNOWTEMPGLwritegrid  ,SURFTEMPGLwritegrid  ,EVRSGLwritegrid ;
	char *SSNOWGLgridname	,*SLIQGLgridname ,*SSTOGLgridname ,*QSNOWGLgridname ,*SNOWTEMPGLgridname  ,*SURFTEMPGLgridname  ,*EVRSGLgridname ;


	/* Niederschlag, Temperatur,Windgeschw.,Dampfdruck,Ausfluss,
		aktueller Radiation Melt Factor, Wasseraequivalent,
		Schneespeicher (fest), Schneespeicher (fluessig), Albedo,
			Landnutzung, Expositionskorrektur fuer Strahlungs-Schmelzfaktor */

	double **Qges; /* array of flows and tracer concentrations for all subbasins and all tracers */
	double r_a,r_b,r_c,r_d; /* Parameter fuer die RMF-Interpolation mit Sinusgang */
	double transzone;
	int  ssnowindex , sliqindex , sstoindex, weindex , albindex, qsnowindex,qsnow_outindex,
	snowdaysindex,snow_ageindex,p_snowindex,p_rainindex,albedoindex;
	int  SSNOWwritegrid,Qwritegrid,qsnow_outwritegrid,SLIQwritegrid,SSTOwritegrid, /* ,WEwritegrid */ 
	Albedowritegrid,SNOWDAYSwritegrid,SNOW_AGEwritegrid,P_SNOWwritegrid,P_RAINwritegrid;

	/* for glaciermodel */
	int  glacmethode;
	char *qfirngridname,*qicegridname,*firnstackname,
	*qfirnstatname,*qicestatname,*qglacstatname,*qsnowonglacstatname,
	*firnstatname,*gmbgridname,*gmbogridname,*gmbstatname,*gmb2statname;
	
	int  qfirnwritegrid,qicewritegrid,firnwritestack,gmbwritegrid,gmbowritegrid;
	int  qfirnindex,qiceindex,firnindex,gmbindex,gmboindex;
	int  iceoldindex;
	int  glacidoldindex;
	int ice_gridslesen,firn_levels;
//	double firn_base, deltafirnper100m;
	double *firn_base, *deltafirnper100m;
	int GMBEvalDay,GMBEvalMonth; // Evaluation day and month for glacier mass balance )e.g. month = 9 and day = 30: each model year on September 30, 24h, the Volume-Area-Relation is applied
	double	 **ICE,**ICEOLD,**GLACID,**GLACIDOLD,**QFIRN,**QICE,**RAD,**GMB,**GMBO,**GLACDEB; // GMB=glacier mass balance, GMBO = glacier mass balance as backup copy when last calculating Volume-area relation)
	double ***FIRN;// stack with dynamically modelled firn areas
	double **FIRN1; // sum of all firn layers
	int  qfirnstatcode,qicestatcode,qglacstatcode,qsnowonglacstatcode,firnstatcode,gmbstatcode,gmb2statcode;
	zonalstatistic *qfirn,*qice,*qglac,*qsnowonglac,*firnstat,*gmbstat,*gmb2stat;
	statrecord *zoneQICE,*zoneQFIRN,*zoneQGLAC,*zoneQSnowOnGlac,*zoneFIRN,*zoneGMB,*zoneGMB2;
	double Kice, Kfirn, Ksnow;

	//	double Asnow_a, Aice_a, Asnow_b, Aice_b, Asnow_akt, Aice_akt;
//	double MF, DDFsnow, DDFfirn, DDFice, Kice, Kfirn, Ksnow;
//	double Asnow_a, Aice_a, Asnow_b, Aice_b, Asnow_akt, Aice_akt, Asnow_min, Aice_min, Asnow_max, Aice_max;
//	double kfirn_els, ksnow_els, kice_els; /* storage coefficients for single lienear storages (glacier melt routing model) */

	// new definition as arrays since these parameters should be available per sub basin
	int DetailedGlacierParams; // 0 if method is <below>= 100, 1 if method is > 100
	int* Codes;
	double *MF, *DDFsnow, *DDFfirn, *DDFice;
	double *Asnow_min, *Aice_min, *Asnow_max, *Aice_max;
	double *Asnow_a, *Aice_a, *Asnow_b, *Aice_b, *Asnow_akt, *Aice_akt;
	double *kfirn_els, *ksnow_els, *kice_els; /* storage coefficients for single lienear storages (glacier melt routing model) */

	double *Qfirn, *Qsnow,
	//*Qsnow_out, 
	*Qice; /* runoff components for each Glacier as areal average (not yet routed!)*/
	//  double *Snowdays,*Snow_age,*P_snow,*P_rain, *ALBEDO;
	double *Qfirn_els, *Qsnow_els, *Qice_els; /* routed runoff components for each Glacier */
	glacierlisttype *glacierlist;
	void GlacierRetreatAndAdvance();
	void CreateGlacierStatisticsZones();
	void CreateNewGlacierCells(int id); // when the glaciers grow, here we add cells to the margins
	void UpdateGlacierStatisticsZones(int id=-1,bool useFirnasWell=false);
	void CountMarginsGroupedByElevationBelts(int id=-1);
	void ComputeOldBalance(int id=-1);
	void RedistributeOldBalance(int id=-1);
	std::vector<glacierclass*> glaciers;
	std::map<long,long> glacierindexmap;
	double VAscaling,VAexponent; // parameters for the Volume-Area-relation in V = VAscaling*A^VAexponent
	int numElevBands; // number of iterations (=elevaltion bands) for grow and erode
	int extraWeightFactorBand0; // for overemphasizing grow and erode for lowest elevation band (default=0, recommended values: 1 ... 3)
	double rho_ice; // density of ice = 850 Kg/m^3 or 0.85Kg/(mm*m^2)
	double dGloabDebrisScaling;

	// when using other zones for unsatzonmodel and routing than for snow model, the addresses of the zoen-connected variables from the unsatzonmodel will be used instead of the normal values
	double      **zn_routing;
	statrecord *zone_routing;
	Tzarray    *zarray_routing;
	int        zonmax_routing;
	int        zonmaxwert_routing;

	// for permafrost thawing model
	double  **ALPHA,		 //	a grid which may be read in as standard grid (if not present: no permafrost model). containes alpha value for thawing depth
	**SNOWFREEDAYS, //  a newly generated grid containing the days without snow (or days with snow under a certain limit). will be used for calculating thaw depth (in conjunction with alpha)
	**SNOWCOVERDAYS, //  a newly generated grid containing the days with snow (or days with snow over a certain limit). will be used for calculating thaw depth (in conjunction with alpha)
	**THAWDEPTH,	 //  the actual thaw depth in the soil
	**ROOTDEPTH;    //  get a reference to the root depth grid in order to limit the root zone to the thaw depth
	int    alphaindex,snowfreedaysindex,snowcoverdaysindex,thawdepthindex,rootdepthindex; 
	int    permfrostmethod;
	double  MinSAgeRefreeze; // minimum snow age to let the soil freeze again (e.g. 10)
	double  MaxSWE_SnowFree; // maximum SWE to be counted as snow-free

	zonalstatistic *thawdepth;
	statrecord *zoneTHAWDEPTH;
	int  thawdepthstatcode;
	char *thawdepthstatname;
	double reduced_ksat; // when using frozen soils in permafrost model, this value defines the ksat for frozen soil layers
						// to which the real ksat is reduced.

	/* variables for the storage-richards or storage.ftz file */
	FILE *storagefile; /* a text file containing the water content of all glacier single linear storages */
	char *storagefilename;

	int sr, er; // startrow and endrow, if rank 0 should do slide calculation
	bool bFullGridRequired;

	bool interceptionBeforeSnow;

	// $kf 2016-07-20: add solid throughfall to snow model input (relevant if interception is computed prior to snow)
	double** SCN_snow_throuhfall;
	double** SCN_ir_radiation;

	double maxSnowPackToCreateGlacier;
	double relSizeNewGlacier;
	double EnBalScaling;

	bool bHeatTransferActive;
	bool bLayeredSnowActive;
	double dTOffsetFreshSnow;
	double ***refSWES; // reference to the layered solid snow stack
	double ***refSWEL; // reference to the layered water content stack
	double **refUpperIndexSnow; // reference to the helper grid holding the index into the swes and swel stacks in which layer the first snow/water can be found
	double ***refThickSnow; // thickness of the snow layers
	double ***refDensSnow; // density of the snow in each layer
	double ***refTemperSnow; // temperature of the snow in each layer

	double **EBAL; // new grid that will be forwarded to unsatzonmodell in order to be able there to estimate the upper boundary condition
	int		ebalindex;
	int  EBALwritegrid;
	char *EBALgridname;
	
	double **EBALGL; // new grid that will be forwarded to unsatzonmodell in order to be able there to estimate the upper boundary condition
	int		ebalglindex;
	int  EBALGLwritegrid;
	char *EBALGLgridname;

	SOILtabletype *LayeredSnowProperties;

	bool dumpmode;
	bool closeAfterDump;
	bool repairAfterDump;
	long maxDumpFiles;
	long actDumpFiles;

	void readSubBasinCodes(int *param);
	void readGlacierParam(double *param);

	schneemodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);
	
//	void zonenidentification(double ***zn_routing,statrecord **zone_routing,Tzarray **zarray_routing,int &zonmax_routing,int &zonmaxwert_routing);

	virtual double calc_Accu_and_Melt(int i, int j, double temperatur, double Anteil, double &schneeanteil, double &Ssnow, double &Sliq, double &Ssto, double &Snowtemp, double &Surftemp, double &Evrs, double &EnergyBalance, double solid_throughfall, bool onGlacier = false);
	virtual void run1(long me, long se, long tg, long mt, long jr);
	virtual void run2(long me, long se, long tg, long mt, long jr);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();

	virtual void  flushbuffers();
	virtual void  SetExpoKorrAdress(double **adresse);
	virtual double **GetSnowStorAdress();
	virtual void  SetQgesAdresse(double **p_Qges);
	virtual double RMF(long tg,long mt);
	virtual double Aice(int subbasin,long tg,long mt);
	virtual double Asnow(int subbasin,long tg,long mt);
	virtual double esatt(double t); /* Berechnet den Saettigungsdampfdruck */
	virtual void   initialize_lakes(std::string instanzname);

	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */

	virtual void ReadSlrContents(char *_storagefilename,int flag);
	virtual void WriteSlrContents(bool appendmode);
	virtual void StorageFileError(const char *message,double c1,double c2);
	virtual void ReadFloatFromStoragefile(double *FValue,const char *ErrorMessage);

	virtual int getMethod() { if (do_it > 0) return methode; else return 0; };
	virtual void setInterceptionBeforeSnow(bool b_interceptionBeforeSnow);

	virtual void setHeatTransferActive(bool flag) 
	{ 
		this->bHeatTransferActive = flag; 
	};
	virtual bool getHeatTransferActive() 
	{ 
		return this->bHeatTransferActive; 
	};

	virtual void setLayeredSnowActive(bool flag,SOILtabletype *LayeredSnowProperties,double dTOffsetFreshSnow);
	virtual bool getLayeredSnowActive() 
	{ 
		return this->bLayeredSnowActive; 
	};

	virtual void UpdateSWEStack(int i, int j, double Anteil, double ChangeSolid, double ChangeLiquid,bool onGlacier, bool fromRefreezing,double externalTemperature);
	virtual void setDumpMode(bool p_dumpmode, bool p_closeAfterDump, bool p_repairAfterDump, long p_maxDumpFiles) { this->dumpmode = p_dumpmode; this->closeAfterDump = p_closeAfterDump; this->repairAfterDump = p_repairAfterDump; this->maxDumpFiles = p_maxDumpFiles; };

};


/* A GridStack is a vector of Grid objects. By adding a grid,
* the lowermost grid is dropped if the GridStack already
* includes 24 grids. In this way, a GridStack object always
* comprises the last 24 grids. A smaller number is also
* possible but it never exceeds the capacity of 24 grids.
* This simple type of data allows one to perform computations
* with data processed in the past.
*/
class GridStack : modelobject
{
public:
	GridStack(int capacity, int mpistartrow, int mpiendrow, int rows, int cols);
	~GridStack();
	bool addToStack(double** grid_in);
	double getAvg(int x, int y);
private:
	int n;
	std::vector<double**>* stack;
	int mpistartrow;
	int mpiendrow;
	int nrows;
	int ncols;
	int capacity;
};



class interzeptionsmodellklasse:
public teilmodellklasse
{
public:

//	friend class GridStack;

	/* Felder der Interzeptionsklasse */
	zonalstatistic  **EIstat,**Qstat,**SIstat;
	char  **EIstatfilename,**SIstatfilename,**Qstatfilename;
	char  **EIgridname,**SIgridname,**Qgridname;
	double **QInput,**dhm,**ICE,**SLIQ; /* SLIQ added for determination of snow surface temperature*/
	double  ***LUSE,***WURZ,***LAI,***VBG,***ETP,***EIP,***ETR,***EI,***SI,***Q,***Z0; /* $kf 2017-06-09: z0 added for threshold detection*/
	int   *siindex,*eiindex,*qindex;
        int sliqindex;
	int iNumVegetationLayers;

	int do_ArcVeg;
	int arcVegNumLayers;
	double*** ARCVEGFRACTIONS; // used for keepting track of the fraction of each landuse type when ArcVeg is used: there may be multiple (usually 5) land use types per cell, these grids keep track of their fraction (created in defaultini)

	statrecord **zoneEI;
	statrecord **zoneSI;
	statrecord **zoneQ;

	int   EIstatcode,SIstatcode,Qstatcode;
	int   EIwritegrid,SIwritegrid,Qwritegrid,SLIQwritegrid;
	double SImax,SIdicke;
	double tagnr,d_dhm[5],lai,vbg; /* Tagnummer, Knicktage im Jahrgang,
								Blattflaechenindex und Vegetationsbedeckungsgrad */
	int use;

	// $kf 2016-06-15: Additional variables for Snow Canopy Modell (SCN)
        // $mk 2018-04-15: More variables for Canopy Energy Balance (CEB)
	double **TInput, **RInput, **WSInput, **RHInput; // grids meteo input
	double **CLInput; // $kf 2017-04-13
	double **SCN_Albedo, **SCN_LOAD, **SCN_AGE; // grids state variables
	double **SCN_SUBL, **SCN_Out_solid, **SCN_Out_liquid;
	double **SCN_Temp, **SCN_RH, **SCN_RAD, **SCN_Wind;
	double **SCN_IR;
        double **CEB_CANENBAL, **CEB_DUMMY, **CEB_CANTEMP, **CEB_CANTEMP_OLD;
	int grd_idx_scn_temp, grd_idx_scn_rh, grd_idx_scn_rad, grd_idx_scn_wind;
	int grd_idx_scn_ir;
	int grd_idx_scn_alb, grd_idx_scn_load, grd_idx_scn_age, grd_idx_scn_subl;
	int grd_idx_scn_out_s, grd_idx_scn_out_l;
        int grd_idx_ceb_canenbal, grd_idx_ceb_dummy, grd_idx_ceb_cantemp, grd_idx_ceb_cantemp_old;
	zonalstatistic *SCN_Temp_stats, *SCN_RH_stats, *SCN_RAD_stats, *SCN_Wind_stats;
	zonalstatistic *SCN_Alb_stats, *SCN_Age_stats, *SCN_Load_stats, *SCN_Sub_stats, *SCN_Out_solid_stats, *SCN_Out_liq_stats;
	statrecord *SCN_Temp_zone, *SCN_RH_zone, *SCN_RAD_zone, *SCN_Wind_zone;
	statrecord *SCN_Alb_zone, *SCN_Age_zone, *SCN_Load_zone, *SCN_Sub_zone, *SCN_Out_solid_zone, *SCN_Out_liq_zone;
	//ceb
	zonalstatistic *CEB_CANENBAL_stats, *CEB_DUMMY_stats, *CEB_CANTEMP_stats, *CEB_CANTEMP_OLD_stats;
	statrecord *CEB_CANENBAL_zone, *CEB_DUMMY_zone, *CEB_CANTEMP_zone, *CEB_CANTEMP_OLD_zone;
	// parameters read from the control file
	bool bSCN_model, CEB_model;
	char *SCN_Temp_filename, *SCN_RH_filename, *SCN_RAD_filename, *SCN_Wind_filename;
	char *SCN_ALB_filename, *SCN_LOAD_filename, *SCN_AGE_filename, *SCN_SUBL_filename, *SCN_Out_solid_filename;
	char *SCN_Out_liq_filename;
        char *CEB_CANENBAL_filename, *CEB_DUMMY_filename, *CEB_CANTEMP_filename, *CEB_CANTEMP_OLD_filename;
	char *SCN_Temp_gridname, *SCN_RH_gridname, *SCN_RAD_gridname, *SCN_Wind_gridname;
	char *SCN_ALB_gridname, *SCN_LOAD_gridname, *SCN_AGE_gridname, *SCN_SUBL_gridname, *SCN_Out_solid_gridname;
	char *SCN_Out_liq_gridname;
        char *CEB_CANENBAL_gridname, *CEB_DUMMY_gridname, *CEB_CANTEMP_gridname, *CEB_CANTEMP_OLD_gridname;
	int SCN_Temp_code, SCN_RH_code, SCN_RAD_code, SCN_Wind_code;
	int SCN_ALB_code, SCN_LOAD_code, SCN_AGE_code, SCN_SUBL_code, SCN_Out_solid_code, SCN_Out_liq_code;
        int CEB_CANENBAL_code, CEB_DUMMY_code, CEB_CANTEMP_code, CEB_CANTEMP_OLD_code;
	int SCN_Temp_writegrid, SCN_RH_writegrid, SCN_RAD_writegrid, SCN_Wind_writegrid;
	int SCN_ALB_writegrid, SCN_LOAD_writegrid, SCN_AGE_writegrid, SCN_SUBL_writegrid, SCN_Out_solid_writegrid;
	int SCN_Out_liq_writegrid;
        int CEB_CANENBAL_writegrid, CEB_DUMMY_writegrid, CEB_CANTEMP_writegrid, CEB_CANTEMP_OLD_writegrid;
	double SCN_factor_lai;
	double SCN_min_lai;
	double SCN_min_z0;
        // $kf 2018-06-15: same rain/snow separation approach in interception model
	double transzone, T0R;
        //CEB: PARAMETERS FOR CEB
        //SURFTEMP
        double **SURFTEMP;
        //ALBEDO_SURF
        double **Albedo;
        //SNOW STORAGE
        double **SSTO;
        //RSC + RSI for LE Transfer Resistance
        double ***RSC, ***RSI;
	
	// $kf 2016-06-20: Add the temperature stack for calculating averages
	GridStack* TempGrids;
	int minute_old,stunde_old, tag_old, monat_old,jahr_old;

	interzeptionsmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void flushbuffers();
	virtual void run(long me,long se,long tg,long mt,long jr);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void  substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void substance_transport_delete(); /* free memeory */
	virtual int getMethod() { if (do_it > 0) return methode; else return 0; };
};

class bodenmodellklasse:
public teilmodellklasse
{
public:

	char *storagefilename;
	FILE *storagefile;

	/* Grids */
	double  **K,**ATB,**SB,**SBalt,**nFK,**SD,**SUZ,**SIF,**P,**ICE,
	**SSNOW,**SSP,**Peak,**QDgrid,**QIgrid,**SOIL,
	**nFKia,**SBia,**FZT;
	double  ***WURZ;
	double **LAKE,**MAXPOND,**POND; // used for lake model: if a cell is within a lake, and the POND-value is positive, the potential evaporation will be calculated after Penman (without monteith)
	int iNumVegetationLayers;

	/* fuer Saettigungsflaechenmodell
				ATB    = Topographischer Faktor ln(A/(T*tan(beta)))
				SB     = Bodenspeicher (Bodenfeuchte als relativer Wert),
				SBalt  = Bodenspeicher aus dem vorhergehenden Intervall, interne Variable ohne explizit
						anzugebenden Namen, aber mit Anspruch auf einen Platz im gridinfo-Array, deshalb
						wird ein konstanter Name vergeben ("$$sb2_grid$$")
				nFK    = nutzbare Feldkapazitaet,
				SD     = Saettigungsdefizit
				K      = K-Wert in m/s,
				P      = Belastung (Niederschlag)
				Smax   = entwaesserbare Porositaet (durch Gravitation), vorerst wird der
						Zeiger auf nFK gesetzt, also nFK = Smax
				SUZ    = Speicher der ungesaettigten Bodenzone
				SIF    = Speicher fuer Interflow
				SSP    = Saturated Soil Partition 0..1 (Anteil an vom GW erreichbaren Bodenspeicher
				Peak   = Peakflow jeder Gridzelle waehrend des gesamten Modellzeitraumes
				QDgrid = Direktabflussgrid
				QIgrid = Interflowgrid
				LUSE   = Landnutzung
				SOIL   = Bodenarten
				WURZ   = Wurzeltiefen
				nFKia  = nutzbare Feldkapazitaet in mm der inaktiven Bodenzone
				SBia   = Bodenfeuchte der inaktiven Bodenzone (relativwert)
				FZT    = Fliesszeitengrid (in Sekunden die Fliesszeit fuer jedes Pixel bis zum naechsten Pegel (
						nur Landoberflaechenabfluss, evtl. mit Faktor auch fuer Interflow anwendbar */

	double    **INFEX;      /* der errechnete Infiltrationsexzess, also Oberflaechenabfluss, wird erst in Infiltrationsmodell genutzt */
	double Ex,Sif,UZ,SDgrenz,kapillarflow,SnowRunoff,SSPneu,Qbb;
	/* Ex           : Saettigungsdirektabfluss,
					Sif          : lokaler Wert des Interflowspeichers
					UZ           : Drainageabfluss aus der unges. Zone in die gesaettigte
					SDgrenz      : Grenztiefe fuer SD (in rueck*nFK), ab der kein Kapillaraufstieg mehr stattfindet
					kapillarflow : aufstiegsrate des GW in den Boden
					SnowRunoff   : Anteil am Schmelzwasser, das wegen Bodenfrost nicht in den Boden eindringen kann
					SSPneu       : Verhaeltnis SD/SDgrenz im aktuellen Intervall
					Qbb          : langsamer Basisabfluss in einem Element */

	/* Grid-Schreib-Codes fuer die wichtigsten Grids */
	int SDwritegrid, SUZwritegrid, SIFwritegrid,
	SSPwritegrid, QDwritegrid, QIwritegrid,
	SBiawritegrid, nFKiawritegrid,fztindex;

	double *ATM; /* ATM: mittlerer Wert von ln(a/tanB) fuer jedes Teilgebiet */
	double *m;    /* Rezessionsparameter */
	double *anffuell,*anfdefi; /* Anfangwerte fuer Fuellung des Drainspeichers und des Saett.defizits in n*nFK [mm] */
	double *code,*Pgrenz,*Tkorr,*Kkorr,*k_ELS_QD,*k_ELS_QI,
	*SD_QI_lim,*rueck,*SR_Anteil,*k_ELS_QBB,*alpha;
	/* code            :
					Pgrenz          : Grenzniederschlag fuer Makroporenansprung,
					Tkorr,Kkorr     : KorrekturFaktor fuer Transmissvitaeten bzw. vertikale Versickerung,
					k_ELS_QD        : Verzoegerung des Direktabflusses im ELS (k) (1..10h),
					SD_QI_lim       : Grenze des Saettigungsdefizites, unterhalb derer SIF gefuellt werden kann
					k_ELS_QI        : Verzoegerung des Interflows im ELS (k) (100..1000 h),
					rueck           : bei Unterschreiten von rueck*nFK durch SD erfolgt Auffuellung des Bodenspeichers aus dem GW (0.1..1.0)
					SR_Anteil       : Anteil an der effektiven Schneeschmelze, der auch bei ungesaettigtem Boden nicht aufgenommen wird,
										sondern direkt abfliesst
					k_ELS_QBB       : Rueckgangskonstante fuer langsamen Basisabfluss
					alpha           : maximale Versickerungsrate ins Tiefengrundwasser bei Saettigung, wird nach exp(-S/m) reduziert
				*/
	double *QUZ,
	*QD,*QDd, /* this arrays are 2-d-arrays (each tracer x subbasin) direct flow */
	*QI,*QIi, /* this arrays are 2-d-arrays (each tracer x subbasin) interflow */ 
	**Qg,       /* total flow and tracer concentration for each tracer and each subbasin */
	*QB,*QBB,*QBBb;  /* die (maximal) vier Abflusskomponenten QD (QDd ist QD mit ELS-Wirkung),
						QI, QB und QBB pro Teilgebiet, vertikale Versickerung QUZ als Teilgebietsmittel
							und Gesamtabfluss Qg (QBBb ist verz. Basisabfluss vor ELS, QBB nach ELS) */
	double *Sm,*Q0;         /* mittl. Speicherdefizit pro Gebiet, max. Basisabfluss pro Teilgebiet */

	char  *SDgridname,     /* saturation defizit mm        */
	*SUZgridname,    /* unsaturated zone storage mm  */
	*SIFgridname,    /* interflow storage mm         */
	*SSPgridname,    /* saturated soil-water-storage partition - (0..1) */
	*QDgridname,     /* direct runoff (surface)      */
	*QIgridname,     /* interflow                    */
	*Peakgridname,   /* max. surface runoff          */
	*SBiagridname,   /* inaktive Bodenfeuchte */
	*nFKiagridname;  /* nutzbare Feldkapazitaet inaktiver Teil */


	zonalstatistic
	*QDstat,*QIstat,
	*QBstat,*QBBstat, /* Statistikstreams fuer die drei (vier) ABflusskomponenten */
	*Qgstat,          /* dito fuer Gesamtabfluss */
	*SBstat,*SUZstat,
	*SIFstat,*SDstat; /* Statistikfiles fuer die Zustaende der Speicher der
							Bodenfeuchte (SB) der ungesaettigten Zone (SUZ), des Interflowspeichers (SIF)
							und des Saettigungsdefizites SD */
	double       *maxzeit;
	statrecord **zoneFZT;  /* fuer die Translation des Oberflaechenabflusses wird
							pro Teilgebiet pro Zeitflaechenzone ein Statistikregistersatz
							benoetigt, ausserdem ein Array mit der Information ueber die
							maximale Anzahl an Registern in jedem Teilgebiet (maxzeit) */


	int  QDstatcode,QIstatcode,QBstatcode,QBBstatcode,
	Qgstatcode,SBstatcode,SUZstatcode,SIFstatcode,SDstatcode; /* Ausgabecodes fuer Statistikfiles */

	char *QDstatfilename,  /* Teilgebietsmittelwert des Direktabflusses             */
	*QIstatfilename,  /* Teilgebietsmittelwert des Interflows                  */
	*QBstatfilename,  /* Teilgebietsmittelwert des Basisabflusses              */
	*QBBstatfilename, /* Teilgebietsmittelwert des verzoegerten Basisabflusses */
	*Qgstatfilename,  /* Teilgebietsmittelwert des Gesamtabflusses             */
	*SBstatfilename,  /* Teilgebietsmittelwert des Bodenspeichers              */
	*SUZstatfilename, /* Teilgebietsmittelwert des Drainspeichers              */
	*SIFstatfilename, /* Teilgebietsmittelwert des Interflowspeichers          */
	*SDstatfilename;  /* Teilgebietsmittelwert des Saettigungsdefizits         */

	int sdindex,suzindex,sifindex,sspindex,
	sbaltindex,qdindex,qiindex,peakindex,
	sbiaindex,nfkiaindex;

//	glacierlisttype *glacierlist;
	int *glacierlist;

	int methode,lake_active,lake_ini_method;

	statrecord *zoneSD,  /* Statistik fuer Saettigungsdefizit  */
	*zoneQI,  /* Statistik fuer Zwischenabfluss     */
	*zoneQBB; /* Statistik der Tiefenversickerung   */

	bool interceptionBeforeSnow;

	bodenmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,const char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

//	virtual void zonenidentification();
	virtual void   ReadFloatFromStoragefile(double *FValue,const char *ErrorMessage);
	virtual void   StorageFileError(const char *message,double c1,double c2);
	virtual void   WriteSlrContents();
	virtual void   ReadSlrContents();
	virtual char   *GetStorageFileName();

	virtual double **GetQgesAdresse();
	virtual void run(long me, long se,long tg,long mt,long jr);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
	virtual void flushbuffers();
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */
	virtual void	initialize_lakes(std::string instanzname);

	virtual void setInterceptionBeforeSnow(bool b_interceptionBeforeSnow);

	// some grids used for soil stability index
	double **CROOT,		// _root_cohesion_ [Pa]
		  **CSOIL,		// _soil_cohesion_ [Pa]
		  **RHOWET,		// _wet_soil_density_ [kg/m^3]
		  **IFA,		// _internal_friction_angle_ [-]
		  **SDEP,		// _soil_depth_ [m]
		  **SLP,        // _slope_angle_ [0..90]
		  **SSI;		// _soil_stability_index_ [-]
	bool doSSI;


};

class infiltrationsmodellklasse:
public bodenmodellklasse
{
public:

	char      *instanznameinfil;
	int        do_infil;

	int        INFEXstatcode;
	char      *INFEXstatfilename;
	zonalstatistic   *INFEXstat;
	char      *INFEXgridname,*SATTgridname;
	int        INFEXwritegrid,SATTwritegrid;
	int        infexindex,sattindex;
	double    **Na;         /* auffuellbares Porenvolumen, wird vor jeder Berechnung aus Bodenspeicher, Parameter */
	double    **SATT;       /* ist auf 1, wenn das Element am Ende des Zeitintervalles gesaettigt war,
							wenn im neuen Zeitintervall dann PI>Ks, dann wird weiterhin Dirketabfluss gebildet */
	double    ***LUSE;      
	double sicker;          /* Anteil an vom Infiltrationsueberschuss doch in den Boden versickernden Wassers */



	infiltrationsmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void flushbuffers();
	virtual void run(long me, long se,long tg,long mt,long jr);
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
};

class SP_class; 

class AL_class {
public:
	wasim::mystream *WSstream; // for stream based output. The stream is given to the constructor as a pointer-parameter
	bool useMinutes;
	double hour_offset;

	int mon[13];
	int    tracer_count;
	tracer_info_type *tracer_info;

	double **AL_buffer;
	double AL_anteil,AL_schwelle,AL_max;
	int AL_art;
	unsigned long AL_bufmax,AL_bufcnt,AL_valcnt;
	char *AL_filename;
	FILE *AL_file;
	char *InFileName;
	FILE *InFile;
	int year,month,day,hour,minute;
	int startyear,startmonth,startday,starthour,startminute;

	double *Abakt; // actual abstraction 
	int TG,AL_code; // counter of ALs for the actual routing description

	int   tabmax;  // Anzahl der Werte, die in Abgaberegeltabelle enthalten sind 
	int	  rulecnt; // Anzahl Stützstellen (Julian-Day-Spalten)
	long *WeekDays;
	int   *JD,*start_hour,*stop_hour,* start_hour2, * stop_hour2;
	double *ARTabIndex,  // Tabelle mit <tabmax> Punkten fuer Inhalt des Speichers oder Abfluss in einem Fluss
	**ARTabValue, // Tabelle mit <tabmax> Punkten fuer Abgabe pro Stützstelle erster Index: Stützstelle, zweiter Index: 1..Tabmax
	**ARChangeRate,  // Steigung der Verbindungsgerade zwischen 2 Inhalt-Abgabe-Punkten fuer jede Stützstelle
	*targetCapacity; // target capacity in m^3/3

	double TgAe,TgAe_rel;
	double InV;
	FileBuffer *FB;
	int column;

	int ZL_target;
	double ZL_TgAe;
	double ZL_TgAe_rel;
	int ZL_flussindex;
	double ZL_SourceAE;
	
	char *DefaultInputDirectory;
	char *DefaultOutputDirectory;

	bool bReverseLookUp;

	// vector for referenced AL-IDs when AL_art==9 (external file that is read in and that replaces other ALs
	std::vector<int> vrefAL_localIDs;
	std::vector<int> vrefAL_globalIDs;
	bool bUseReplaceValue;
	double dReplaceValue;

	AL_class();
	~AL_class();
	virtual void  ExtendFileName(char *directory,char *dateiname);
	virtual void  AL_flushbuffer(int code);
	virtual void  WriteAbstraction(long p_jr,long p_mt,long p_tg,long p_se,long p_me,double *flowm3s);
	virtual void  SetAbstractionMode(FILE *steuerdatei,int p_ALnr,std::stringstream &mystream,int subint,int jahra,int monata,int taga,int stundea,int minutea);
	//	virtual void  SetAbstractionMode(FILE *steuerdatei,int p_ALnr,char *satz1,int subint);
	//	virtual void  SetAbstractionMode(int p_ALnr,int p_AL_art,const char *p_AL_filename,
	//				  double p_AL_schwelle,double p_AL_anteil,double p_AL_max,const char *p_AL_einheit,int p_subint);

	virtual void  ReadAbstractionRule(FILE *steuerdatei);
	virtual void  findeabschnitt(FILE *steuerdatei,char *Abschnittsname);
	virtual void  initialize(wasim::mystream *p_WSstream,tracer_info_type *p_tracer_info,int p_tracer_count,int p_TG,double p_TgAe,double p_TgAe_rel,double p_InV,
	char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,double hour_offset);
	virtual double CalcWaterBalance(double **Qges,double &abstraction_sum,double &flowmax,int p_year,int p_month,int p_day,int p_hour,int p_minute);
	virtual double m3s_to_mmh(double p_m3s,double p_TgAe,double p_InV);
	virtual double mmh_to_m3s(double p_mmh,double p_TgAe,double p_InV);
	virtual void  SetTargetData(int p_ZL_target,double p_ZL_TgAe,double p_ZL_TgAe_rel,int flussindex,double SourceAE);
	virtual double GetActualAbstraction(double **Qges,double Qact,int p_year,int p_month,int p_day,int p_hour,int p_minute);
	virtual double GetActualAbstractionFromReservoir(double ***fluss,double Vnew,int p_year,int p_month,int p_day,int p_hour,int p_minute);
	virtual int get_weekday(int y, int m, int d);
};

class SP_class {
public:
	double hour_offset;
	bool useMinutes;
	int timestep;

	int mon[13];
	wasim::mystream *WSstream; // for stream based output. The stream is given to the constructor as a pointer-parameter
	// variables from modelobject, which will be required for analyzing the lakes statistics
	gridinfotype *gridinfo;
	int *zarray;
	statrecord *zone;
	firstlast *FirstLastCol;
	int zonmax;
	int zonmaxwert;
	int ncols, nrows, mpistartrow, mpiendrow;
	double xll,yll,cellsize,nodata;
	int lake_active;
	int year,month,day,hour,minute;

	int    tabmax;  /* Anzahl der Werte, die in Abgaberegeltabelle enthalten sind */
	int    tracer_count;
	tracer_info_type *tracer_info;
	double Vmax,    /* maximaler Inhalt -> darueber geht alles in Abgabe */
	*Vakt,   /* aktueller Inhalt und aktuelle Konzentrationen von bis zu 9 tracern */
	*Vold,   // volume and tracer concentrations from the last time step -> stored at the beginning of a routing step to handle balances
	*Vtemp,  // temporal values as function result of GetVolumeFromPonds(double *V)
	Abakt;   /* aktuelle Abgabe */
	std::list<double>	V_hist; // for running average of volume over 24 hours. Number of elemenst in this queue is defined by time step lenght an may be by an additional parameter of the time window required (usually 24 hours for now)

	double WaterlevelInit, // optional initialization water level when using routing mode and a hydraulic connection between two lakes should be regarded -> lake mnodel ist switched off then and no water levels are known
		AreaInit; // initialization area when using routing mode 

	//           *Inhalt, /* Tabelle mit <tabmax> Punkten fuer Inhalt des Speichers */
	//           *Abgabe, /* Tabelle mit <tabmax> Punkten fuer Abgabe */
	//           *ARChangeRate;  /* Steigung der Verbindungsgerade zwischen 2 Inhalt-Abgabe-Punkten fuer */
	int	  rulecnt; // Anzahl Stützstellen (Julian-Day-Spalten)
	long *WeekDays;
	int   *JD,*start_hour,*stop_hour, * start_hour2, * stop_hour2;
	double *ARTabIndex,  // Tabelle mit <tabmax> Punkten fuer Inhalt des Speichers oder Abfluss in einem Fluss
	**ARTabValue, // Tabelle mit <tabmax> Punkten fuer Abgabe pro Stützstelle erster Index: Stützstelle, zweiter Index: 1..Tabmax
	**ARChangeRate,  // Steigung der Verbindungsgerade zwischen 2 Inhalt-Abgabe-Punkten fuer jede Stützstelle
	*targetCapacity; // target capacity in m^3/3

	double AREAm2, // balance entities for later output to the baance file
	ABSTRmm,ABSTRm3s,ABSTRm3, 
	IRRIGmm,IRRIGm3s,IRRIGm3, 
	PRECmm,PRECm3s,PRECm3,
	QINmm,QINm3s,QINm3,
	QOUTmm,QOUTm3s,QOUTm3,
	QGWmm,QGWm3s,QGWm3,
	ETRmm,ETRm3s,ETRm3,
	DVmm,DVm3s,DVm3,
	ERRmm,ERRm3s,ERRm3;

	double deltaPond; // the change in PONDs caused by PREC, ETR and GWN during an external time step 

	double dTmin; // internal sub time step, can be as low as less than a second for e.g. hydraulic connections between lakes or for small lakes with large inflow/outflow
	int UseRunningAvg; // if <=1, the outflwo of the lake is calculated from the actual content using dTmin (i.e. several times per interval, the effective outflow is then summed up for all sub time steps)
	// if > 1, the outflow is determined by the running average of the reservoirs volume over the last n time steps, with n as the value of this parameter

	double **SP_buffer;
	unsigned long SP_bufmax,SP_bufcnt,SP_valcnt;
	char *SP_filename;
	FILE *SP_file;

	// only used, when LAKE grid is present
	double **LAKE,**MAXPOND,**POND,**zn,**dhm,**ETR,**EI,**qu,****conc,**DELTAGW,**PREC;
	double *WaterLevel; // depth for double referenced depth-volume-table
	double *Volume;     // volume for double referenced depth-volume-table
	int W_V_total_steps;
	int startcol,startrow,endcol,endrow; // grid rows and columns which must be scanned for this lake (it's a performance issue only)
	int SP_code; // code given in the control file in the routing description. Must match the code in the LAKE-grid
	int TG; // code for the subcatchment -> must match the zonegrid for unsatzon-model where the actual lake is located
	int readgrids;
	char *DefaultInputDirectory;
	char *DefaultOutputDirectory;

	bool bReverseLookUp;
	bool useHydraulicConnection;
	int 	TargetReservoir;
	double 	ChannelWidth;
	double 	ChannelBottom;
	double 	ChannelLength;
	double 	RoughnessCoeff;
	//double  ActFluxm3s;
	double  TargetWaterLevel;
	double  TargetWaterLevelChangeRate;
	double  actWaterLevel;
	std::map<int, double**> *mapReservoirLevels; // array of three elements of pointers to double: [0]=water level, [1]=VolumeChange, [2]=reserve
	std::map<int, bool> *mapTributariesWithHydCon;

	SP_class();
	~SP_class();
	virtual void  ExtendFileName(char *directory,char *dateiname);
	virtual void  SP_flushbuffer(int code);
	virtual void  WriteSpeicherinhalt(long p_jr,long p_mt,long p_tg,long p_se,long p_me);
	virtual void  SetSpeicherParams(FILE *steuerdatei,int p_SPnr,std::stringstream &mystream);
	virtual void  ReadReservoirRule(FILE *steuerdatei);
	virtual void  findeabschnitt(FILE *steuerdatei,char *Abschnittsname);
	virtual void  initialize(wasim::mystream *p_WSstream,tracer_info_type *p_tracer_info,int p_tracer_count,
							int p_TG,gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,
							firstlast *p_FirstLastCol ,int p_zonmax,int p_zonmaxwert,
							int p_ncols,int p_nrows,int p_mpistartrow,int p_mpiendrow,double p_xll,double p_yll,
							double p_cellsize,double p_nodata,int p_lake_active,
							char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,double hour_offset, int timestep);
	virtual void  initialize_lake();
	virtual double GetWaterLevelByVolume(double V);
	virtual double GetVolumeByWaterLevel(double W);
	virtual double GetVolumeFromPonds(double *V); // an array of 0..tracer_count doubles which will hold the total Volume as well as the tracer concentrations in its elements
	virtual	double GetActualOutflow(double ***fluss,AL_class *Abstraction,int ALs,double Qin,double InV,double hour_offfset,int p_year,int p_month,int p_day,int p_hour,int p_minute,double dOldOutflow=0,double iMinTimeStep=-1);
	virtual double SetPonds(double Volume);
	virtual void  SetVaktFromPonds();
	virtual void  CreateStatFile();
	virtual void  CalcWaterBalance(double InV,double irr,double abstraction_summ3,double ***fluss,double **Qges,double *flowm3s,AL_class* Abstraction,int ALs,int p_year,int p_month,int p_day,int p_hour,int p_minute,int p_timestepfactor,int actsubtimestep,int numsubtimesteps);
	virtual double mixing(double c1,double c2, double v1,double v2);
	virtual void  FillMapReservoirLevels(std::map<int, double**> *p_mapReservoirLevels);
	virtual void  FillMapTributariesWithHydCon(std::map<int, bool> *p_mapTributariesWithHydCon);
	virtual int get_weekday(int y, int m, int d);
};

/* Eine Struktur, die fuer jedes Teilgebiet, das sich aus mehreren Oberliegern und dem letzten Zwischengebiet zusammensetzt,
Parameter fuer das Routing speichert. Die Arraydimension von routing ist gleich der Anzahl an zu routenden Tilgebieten */
class routingclass:
//public memory_handling
public modelobject
{
public:
	int tracer_count;
	int readgrids;
	tracer_info_type *tracer_info;
	wasim::mystream *WSstream; // for stream based output. The stream is given to the constructor as a pointer-parameter
	int lake_active;
	riverflow_table_type *riverflow_table; /* for communication between unsaturated zone and rivers -> infiltration */
	double v1,v2,weightsum,weight,*frachtsum,*fracht,frachtweightsum,frachtweight,discharge,dischargesum;

	int timestepfactor; // if routing is done several times per model time step (=externalTimeStep), then this value contains the factor. E.g. 24, if unsatzonmodel is run in 1d time step but routing must run in 1h time step because of abstraction rules

	double TgAe; /* Einzugsgebiet des Gesamten Teilgebietes inclusive Oberlieger, wird zur
				Berechnung der Spenden aus Zuleitungen und zur Berechnung der m3/s aus ABleitungen in mm benutzt */
	double TgAe_rel; /* relativer Anteil der als Flaechen der als Oberlieger gereouteten Teilgebiete an der
					Gesamteinzugsgebietsflaeche des Teilgebietes. Wenn z.B. ein einzelner Oberlieger nicht geroutet wird,
					weil dessen Zulaeufe ungenau modelliert werden, dieses TG aber im Modell auch modelliert wird, so muss
					TgAe_rel das Beruecksichtigen (Bsp.: TgAe_rel ist 0.8279, wenn ein 84 km2 grosser Oberlieger
					an einem insgesamt 488.1 km2 grossen Teilgebiet nicht mit geroutet wird, sondern diese Daten
					von einer Datei eingelesen werden. Wenn die Daten nicht eingelesen werden, dann MUSS TgAe_rel == 1.0 sein,
					weil ja von vornherein mit Spenden gerechnet wurde. */
	int zonmax,TG,OLs,*OL,SUMTRIBs,**indexST,ets,*et;       /* OLs: Anzahl Oberlieger, TG: fuer welches Gebiet, OL eine Liste mit den Nummern der Oberlieger */
	int SPs,*SP,ZLs,*ZL,ALs,*AL;  /* Anzahl Speicher, Zulaeufe, Ablaeufe und Nummern als Liste */
	int *ZL_art;//,*AL_art;  /* fuer jeden Zu- oder Ablauf wird mit einem Code die ARt festgelegt:
	//                           1 constanter Wert in m3/s
	//                           2 relativer Wert am Gesamtabfluss (0..1) (nur Ableitung)
	//                           3 file Werte aus Datei (nur Zuleitung)
	//                           4 intern, Werte werden aus vorher berechneter Ableitung uebernommen (nur Zuleitung)
	//                           5 extern, Werte werden in Datei geschrieben) */

	//  double ***ZL_buffer;
	FileBuffer *ZL_buffer; // array of Zulauf-buffers for 1..n ZL-definitions
	//  double ***AL_buffer;
	//  double ***SP_buffer;
	//  unsigned long *ZL_bufmax,*AL_bufmax,*SP_bufmax;
	//  unsigned long *ZL_bufcnt,*AL_bufcnt,*SP_bufcnt;         /* actual index for buffer row (1...bufmax) */
	//  unsigned long *ZL_valcnt,*AL_valcnt,*SP_valcnt;         /* index forthe buffer column */

	//  char  
	//		 **ZL_filename,
	//       **SP_filename,
	//       **AL_filename;   /* fuer jede der theor. beliebig vielen Zu- oder
	//                                         Ableitungen kann ein DAteiname angegeben
	//                                         werden (bei ABleitungen nur, wenn "extern"
	//                                         angegeben ist */
	//  int   *ZL_filecol;
	//  long  *ZL_filerow;      /* fuer jedes Zulauffeil muss angegeben werden, in welcher Spalte die Werte stehen und ab welcher Reihe sie gelten */
	//  FILE **ZL_file,**AL_file,**SP_file; /* File, in welchem die Zulaeufe festgehalten werden. Muss nicht, kann aber immer dasselbe File sein */
	//  FILE **AL_file; /* File, in welchem die Zulaeufe festgehalten werden. Muss nicht, kann aber immer dasselbe File sein */
	double **ZL_wert;//,*AL_anteil,
	//        *AL_schwelle,*AL_max; /* ZL_wert: der Wert, der fuer jeden Zulauf aus seiner DAtei gelesen wurde,
	//                               AL_anteil: relativer ANteil des die Schwelle AL_schwelle uebersteigenden Wassers, dass abgeleitet wird
	//                               AL_max; maximales Fassungsvermoegen der Ableitung */
	//  char **ZL_format,*ZL_tempstr;       /* Leseformat fuer die entsprechende Zulaufdatei, wird durch ZL_spalte bestimmt, Rest der Zeile wird in ZL_tempstr gelesen */
	int   ZLjr,ZLmt,ZLtg,ZLse,ZLme;         /* aktuelles Lesen aus der Zulaufdatei->Speichern des Datums */
	double  *flowmmh,*flowm3s,flowmax;     /* fuer Speicherinhalts- und Ableitungsberechnungen */

	double InV; /* Zeitschrittweite in h */
	double *kv,*kh,*gew,**total; /* t Translationszeiten, kv/kh Speicherkonstante Vorland, Hauptbett, *gew relativer Flaechenanteil (pro OL) */
	double *Bh,*Bv,*Mh,*Mv,*Tm,*I,*L,*AE,*minspende,*maxspende,*intzahl,*dq;
	int   subint;
	/* fuer hydraulische Berechnung:
				*Bh, *Bv Liste mit den Breiten der Gerinne fuer jeden Oberlieger (Hauptbett, Vorlaender)
				*Mh, *Mv Liste mit den Manningbeiwerten der Gerinne fuer jeden Oberlieger (Hauptbett, Vorlaender)
				*Th Liste mit den Hauptgerinnetiefen fuer jede Gerinnestrecke
				*I Liste mit den Gefaellewerten (in m/m) fuer die Gerinne
				*L Liste mit den Laengen (in km) jedes Gerinnes
				*AE Liste mit den absoluten EZG-Groessen (eventuell als Mittelwert angebbar) zur Berechnung des
					Gesamtabflusses zur Abschaetzung der Fliessgeschwindigkeit (kann groesser sein, als das Einzugsgebiet
					eines Oberliegers, da es auch das EZG bis zur Muendung in den naechsten bedeutenden Fluss erfassen soll
				*minspende minimale Abflussspende, fuer die noch eine Translationszeit in wqtabelle geschrieben werden soll
				*maxspende maximale Spende ------------"-----------------"----------------"----------------"---------------
				*intzahl   Anzahl Intervalle, in die die Abflusskurve unterteilt werden soll, um zwischen minspende und
							maxspende zu interpolieren
				*dq        Aenderung des Abflusslogarithmus zwischen zwei in der wq-tabelle enthaltenen Spenden
				subint    Anzahl an Teilzeitintervallen, in die ein Zeitintervall unterteilt werden soll
			*/

	/* Tabelle fuer eine Spenden-Geschwindigkeits-Beziehung pro Gerinnestrecke */
	struct w_q_tabelle {
		double *sp_lskm2,
		*sp_mm,
		*Q;         /* Spende in mm/InV, l/s/km2 bzw. Abfluss in m3/s */
		double *v,*vh,*vv; /* Fliessgeschwindigkeit in m/s (gesamt sowie fuer Hauptbett und Vorland*/
		double *t;         /* Fliesszeit in h */
		double n;          /* Unterteilungsanzahl gemaess des Verhaeltnisses aus groesster zu kleinster Translationszeit */
		double *Rh,*Rv;    /* hydraulische Radien in m (Hauptbett, Vorland */
		double *Th,*Tv;    /* Wassertiefe in m fuer Hauptbett und Vorland */
		int n12t;         /* Anzahl von 1/12 Intervallen, die die laengste Translationszeit fuellen */
		int *delta;       /* Anzahl an 1/12-tel Intervallen, die ein bestimmter Abfluss vom Auslass
						entfernt abgelegt werden muss, damit er beim konstanten Zeitschritt in
						der richtigen Zeit dort ankommt */
		double *QvzuQ;     /* Verhaeltnis des Vorlandabflusses zum Gesamtabfluss (ist null, wenn keine Ausuferung eintritt */
		int *last_index;  /* merkt sich fuer jede (eventuell geteilte) Strecke den Index der 1/12-tel Abflussarrays,
						in dem der Abfluss der Vorstunde endet.
						So koennen Uebrrollen und leere Zwischenraeume beim Routen der Abfluesse der
						naechsten Stunde festgestellt werden.
						Muss mit dem delta-Wert des Initialabflusses initialisiert werden */
		double **ueberrollt; /* speichert sich fuer jede Gerinneteilstrecke die Menge von einer schnellen HW-Welle
						ueberholten Abflusses */
		double *qvzuq;         /* haelt Verhaeltnis von Vorland-zu Gesamtabfluss vom vorigen Zeitintervall fest */
		double **Qaus;       /* Array, das den Ausfluss im vorigen Zeitintervall speichert; hat fuer jede
						Teilstrecke einer Routingstrecke <subint> Elemente */
	} *wqtab;

	SP_class *Speicher;
	AL_class *Abstraction;

	/* temporaere Variablen fuer die W-Q-Tabellen-Berechnung,
	teilweise auch fuer das Routing (Qh, Qv, Qh_alt,Qv_alt) */
	double vh,vv,Rh,Rv,Qh,Qh_alt,Qv,Qv_alt,Th,Th_alt,Tv,Tv_alt,qh,qv,lnq;

	/* Variablen fuer das eigentliche Routing */
	int geb, /* = zarray[OL[i]] */
	mi,  /* maximale Anzahl Intervalle einer Teilstrecke */
	ber, /* Abflussbereich (Index in wqtab) */
	ver; /* aktuelle verschiebung */
	double diff,Qref; /* diff: fehlendes Wasser zwischen zwei gerouteten Intervallabfluessen,
					das dann aus ueberrollt und ggf. der vorhergehenden oder
					der nachfolgenden Welle ausgeglichen wird;
					Qref welcher Abfluss soll als referenz fuer Auffuellen/Speichern
					von ueberrolltem Wasser dienen, der letzte oder der aktuelle
					Abfluss? in Qref wird der entsprechende Wert je nach aktuellem
					Anstieg der v-Q-Beziehung festgehalten */
	double  *CQref; /* concentrations of max. 9 tracers in Qref */
	int     gem_Zul; /* 0, wenn keine gemessenen Zulaufganglinien genutzt werden solln, sonst 1 */
	double **Spende;  /* gemessene Abfluesse, um eventuell gemessene Zulaufganglinien zu routen */
	int    *Spalte;  /* fuer Spenden ist die Zuordnung der Spalten in der Spendendatei zu den Teilgebieten noetig */
	double  **abfrel;  /* relative Anteile des Teilgebietes am Gesamtgebiet */
	int *zarray;     /* Umindizierung der Teilgebietscodes (1...16384 in Codes 1...zonmax) zum Speichersparen */
	/* Uebernahme der Adressen aus dem abflussklassenobject,
	um auch in routingclass darauf zugreifen zu koennen */
	double **Qges,***Q12ges,***fluss,**Qakt;
	statrecord *zoneSW;
	double **maxirrig;
	double **maxirrigFromReservoirs;
	double abstraction_sum;

	std::map<int, double**> *mapReservoirLevels; // array of three elements of pointers to double: [0]=water level, [1]=VolumeChange, [2]=reserve
	std::map<int, bool> *mapTributariesWithHydCon;

	/* for storing actual contentof river channels and single linear storages for the next model run */
	virtual void WriteStorageFile(FILE *storagefile);
	virtual void ReadStorageFile(FILE *storagefile,char *storagefilename,double timestepstor, int syntaxVersion);
	virtual void ConvertStorageFile(FILE *storagefile, int z, int tribindex, int geb, int src_subint, int src_maxdelta, int src_numchannels, double src_timestep,int mode=1);

	virtual void SetZoneSWAdress(statrecord *p_zoneSW);
	//  virtual void findeabschnitt(FILE *steuerdatei,char *Abschnittsname);
	virtual void riverflow_table_ini(riverflow_table_type *adresse,int nr);
	virtual void riverflow_table_update();
	virtual int   getmaxdelta(int tmax);
	virtual int   getbereich(int ol,double q);
	virtual double maxflow(double *flow);
	virtual double GetAbstractionSum();

	virtual double decay(double conc, double decfact);               /* for calculating Ci = Ci-1 * decfact; decfact = exp(-(InV/60/24)/half_time);
																	returns the new concentration after dte decay */
	virtual double mixing(double c1,double c2, double v1,double v2); /* calc. mixing concentrations for volume1 (c1) with volume 2 (c2);
																	returns mixing concentration */
	virtual double evaporating_changes(double c,double v1,double ve); /* calc. and returning the concentration c after evaporating the amount ve from v */

	//  virtual void  WriteSpeicherinhalt(int i,double *pwert,long p_jr,long p_mt,long p_tg,long p_se);
	//  virtual void  WriteAbleitung(int i,double *pwert,long p_jr,long p_mt,long p_tg,long p_se);
	virtual void  ReadZuleitung(double **p_Qges,long p_jr,long p_mt,long p_tg,long p_se,long p_me);

	virtual void  flushbuffers();
	//  virtual void  AL_flushbuffer(int code,FILE *datei,char *dateiname,double **buffer,unsigned long bufcnt,unsigned long wertcnt);
	//  virtual void  SP_flushbuffer(int code,FILE *datei,char *dateiname,double **buffer,unsigned long bufcnt,unsigned long wertcnt);

	virtual void  SetSpeicherParams(FILE *steuerdatei,int p_index,std::stringstream &mystream,int p_subint);
	virtual void  LocalizeRefALs();
	//  virtual void  SetSpeicherParams(int p_index,int p_SPnr,double *p_V0,char *p_SP_filename,int p_subint);
	virtual void  ReadReservoirRule(FILE *steuerdatei);
	virtual void  FillMapReservoirLevels(std::map<int, double**> *p_mapReservoirLevels);
	virtual void  FillMapTributariesWithHydCon(std::map<int, bool> *p_mapTributariesWithHydCon);
	virtual void  SetFirstAbstractionValue();

	virtual void  ReadAbstractionRule(FILE *steuerdatei);
	virtual void  SetAblaufModus(FILE *steuerdatei,int p_index,std::stringstream &mystream,int p_subint,long jahra,long monata,long taga,long stundea,long minutea);
	//  virtual void     SetAblaufModus(int p_index,int p_ALnr,int p_AL_art,const char *p_AL_filename,double p_AL_schwelle,double p_AL_anteil,double AL_max,const char *p_AL_einheit,int p_subint);
	virtual void     SetZulaufModus(int p_index,int p_ZLnr,int p_ZL_art,const char *p_ZL_filename,
									long p_ZL_reihe,int p_ZL_spalte,double p_ZL_wert,
									long jahra,long monata,long taga,long stundea,long minutea);
	virtual void  SetRoutingParams(int p_index,int p_OL,int *p_numST,double p_kh,double p_kv,double p_Bh,double p_Bv,double p_Tm,double p_Mh,double p_Mv,
	double p_I,double p_L,double p_AE,double p_minspende,double p_maxspende,double p_intzahl,int p_subint,int p_gem_Zul,bool bUseTgAe_rel=true);
	virtual void  SetWeights();
	virtual void  Print_W_Q_Table();
	virtual void  IterateManning(double *Q,double *B,double *T,double *R,double *v,double M,double I,int maxiter,double maxdiff);
	virtual void  Routing(long jahr, long monat, long tag, long stunde, long minute, int actsubtimestep, int numsubtimesteps);
	virtual void  SetFluxAdresses(double **p_Qges,double ***p_Q12ges,double ***p_fluss,double **Spende,int *Spalte);
	virtual void  set_initial_last_index();
	virtual void  SetMaxIrrigAdress(double **p_maxirrig,double **p_maxirrigFromReservoirs);

	virtual double m3s_to_mmh(double p_m3s,double p_TgAe,double p_InV);
	virtual double mmh_to_m3s(double p_mmh,double p_TgAe,double p_InV);
	virtual void add_reservoir(std::vector<SP_class*> &reservoirs);
	/* Constructor und Destructor */

	routingclass(int &retvalue,wasim::mystream *p_Stream,char *steuername,int LogLevel,double p_InV,
		int p_zonmax,int anz_OL,int anz_SUMTRIB,int anz_et,int anz_ZL,int anz_AL,int anz_SP,int tg,
		double p_TgAe,double p_TgAe_rel,double **p_abfrel,int *p_zarray,
		int p_tracer_count, tracer_info_type *p_tracer_info_table,
		gridinfotype *p_gridinfo,statrecord *p_zone,
		firstlast *p_FirstLastCol ,int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata, int p_lake_active, int p_readgrids,
		char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,
		long *pv_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, std::vector<HRUMapElement*> *p_vHRU,
		double hour_offset, int timestep,int timestepfactor,
		int minutea,int stundea,int taga,int monata,int jahra
	);
	virtual void  ClearMemory();
	/*  virtual      ~routingclass(); */ /* es wird ein Standarddestructor verwendet */

};

class abflussmodellklasse:
public teilmodellklasse
{
public:
	char *storagefilename;
	FILE *storagefile;
	int lastresultsindex;
	TLastResults *LastResults;
	riverflow_table_type *riverflow_table; /* for communication between unsaturated zone and rivers -> infiltration */
	int anz_OL,anz_ZL,anz_AL,anz_SP,anz_SUMTRIB,anz_et;  /* Anzahl OL, ZL, AL uns SP pro Routingvorschrift */
	int p__OL,p__ZL,p__AL,p__SP,p__et,p__SUMTRIB; /* temporaere Speicher fuer anz_OL usw. */
	int  *num_ST; /* enthaelt die Nummern der Oberlieger, die fuer ein 
					SUMTRIB-routing addiert werden muessen -> array wird 
					als Parameter an SetSUMTRIBRoutingParams() uebergeben */
	char *Qstatfilename,*spendenfilename;
	char *satz1,*tsatz1,*satz2,*tsatz2,*satz3,*tsatz3;
//	char *kombisatz,*kombisatztemp,*tkombisatz;
	int numcoldescriptors; //,kombisatzlaenge;
	int Qstatcode;
	statrecord *zoneSW;
	double **maxirrig;
	double **maxirrigFromReservoirs;
	double *maxirrig_temp;
	int **rtarget; /* target catchment for tha actual routing description */
	int timeoffset,time;  /* timeoffset ist die Zeit, die benoetigt wird, um am Ausfluss echt geroutete Abfluesse zu landen,
						wenn die Summe aller Translationszeiten im Gebiet bis zum Gesamtauslass z.B.
						14 h ist, dann sollte timeoffset=14 gewaehlt werden, time zaehlt die Intervalle,
						wenn time > timeoffset, dann wird R-square gerechnet */
	FILE *spendenfile;
	FileBuffer *SpendenBuffer;

	char *qgesfileRMname; // name of the input file for total runoffs from subbasins when in routing mode
	FileBuffer *QgesBufferRM; // buffer to hold the read in values from external total runoffs when in routing mode (has a member that represents the external file)
	int qgesfileend;
	int lastday; // used for printing a message, which should be printed for each new day, so everytime, the actual day is different from the lastday, the message will be printed and the lastday will be be updated
	char *additionalinfo;

	routstatistic *Qstat; /* Instanz einer von zonalstatistic abgeleiteten Klasse:
						Es wird kein GEsamtgebietsmittelwert ausgegeben,
						Die Relativanteile sind Relativanteile der Pegelgebiete am Gesamtgebiet,
						deren Summe groesser 1 ist */

	int gem_Zul,anzrout;  /* gem_Zul: wenn gemessene Zulaufganglinien genutzt werden sollen, wird hier eine 1, sonst eine 0 gespeichert
						anzrout: Anzahl an nacheinander zu berechnenden Routings, jedes Routing kann aus beliebig vielen
						Oberliegern zusammengesetzt sein, die aber alle nur zu einem gemeinsamen Ziel geroutet werden */
	double **abfrel;        /* relative Anteile der Pegelgebiete (nicht der Teilgebiete) am Gesamtgebiet */
	double ***fluss,**Qges,**Qges_uebernommen,***Q12ges,**Spende;
	int Ys,Ms,Ds,Hs,ms;      /* Jahr,monat,Tag,STunde fuer Spenden (zur Kontrolle mit Modellzeit) */
	/* Zuordnung der Spalten dr Spendendatei zu den Teilgebieten: in spalte[1]
	steht die Spalte der Spendendatei, die fuer TG 1 gilt usw. */
	int *Spalte;
	double tmax;
	int spendenende,spendenanzahl,p_subint;
	double p_minspende,p_maxspende,p_intzahl,p_kv,p_kh,p_Bv,p_Bh,p_Mv,p_Mh,p_Th,p_I,p_L,p_AE,p_AErel;
	routingclass  **routing;
	int runoff_measure; /* 0=mm/interval, 1=m3/s */

	int MinToInvEnd;
	int extTimestep;
	double timestepfactor;

	/* Fuer jeden Pegel aus der Spendendatei, bzw. fuer jedes einer Spalte in der Spendendatei zugeordnete Teilgebiet
	wird ein Effizienzmass ausgerechnet (R-square). Dafuer werden einige Statistikregister benoetigt. Da
	dieses Effizienzmass auf +- 10 Intervalle bezogen wird, wird gleich ein ganzes Feld an Statistikregistern
	reserviert. Es wird ausser dem linearen (effizienz) noch das logarithmische R-square berechnet -> log_effi */
	Tstatist **effizienz,**log_effi;

	int lake_active;
	int readgrids;
	int lakeindex;
	int pondindex;
	int etrindex, eiindex,quindex, qbindex, pindex;
	int tracerindex[10];
	long *lakeCellsPerMPIRank;
	int do_it_mpi0;
	double **PRECONLAKES;
	int preconlakesindex;

	bool interceptionBeforeSnow;

	std::map<int, double**> *mapReservoirLevels; // array of three elements of pointers to double: [0]=water level, [1]=VolumeChange, [2]=reserve
	//std::map<int, double*> *mapReservoirLevels; // used to allow any reservoir in any routing class instance to access the water levels of any other reservoir
	std::map<int, bool> *mapTributariesWithHydCon; // indicates for each potential tributary if the reservoir there has defined a hydraulic connection to the actual 

	/* der Constructor initialisiert alle Datenstrukturen, erzeugt die noetigen Ausgabedateien, errechnet Gewichte,
		reserviert Speicher usw.*/
	abflussmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *p_LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end, double MinutesToInvEndRouting, int externalTimeStep,
		int numsubb,int *listsubb,
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
		int p_tracer_count, tracer_info_type *p_tracer_info_table,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int p_readgrids);

	virtual int count_substr(char *satz,const char *substring);
	virtual void StorageFileError(const char *message,double c1,double c2);
	virtual void ReadStorageFile(char *filename,int flag);
	virtual void WriteStorageFile();

	/* Um die Adresse des Arrays mit den Teilgebietsabfluessen zu erhalten wird eine eigene Methode benoetigt,
		nur so kann die Kapselung der Objekte erreicht werden */
	virtual void SetRoutingFluxAdresses();
	virtual void SetQgesAdresse(double **adresse);
	virtual void SetRiverflow_tableAdress(riverflow_table_type *adresse);
	virtual int **GetRtargetAdress();
	virtual double **GetAbfrelAdress();
	virtual void SetRunoffMeasure(int p_runoff_measure);

	/* fuer Bewaesserung: Adresse des Zonenarrays mit dem mittleren Bewaesserungsbedarf aus Oberflaechenwasser */
	virtual void SetZoneSWAdress(statrecord *p_zoneSW);
	virtual double **GetMaxIrrigAdress();
	virtual double **GetMaxIrrigFromReservoirsAddress();

	virtual void flushbuffers();
	/* Ausfuehrungsroutine des Abflussmodells */
	virtual void run(long me, long se,long tg,long mt,long jr,int actsubtimestep, int numsubtimesteps);
	/* Zerstoert die Instanz der Klasse, gibt Speicher frei */
	/* anstelle des Destructors wird Speicher freigegeben und Standarddestructor vom Compiler erzeugt  */
	virtual void ClearMemory();
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void substance_transport_delete(); /* free memeory */
	virtual void initialize_lakes(std::string instanzname);
	virtual void CountMPICells(bool b_interceptionBeforeSnow);
	virtual void get_reservoir_vector(std::vector<SP_class*> &reservoirs);

	virtual void routingmode_initialize();
	virtual void routingmode_nextrecord();
	virtual void routingmode_cleanup();
	virtual bool routingmode_isend();
	virtual void getdatetime(long *me, long *se, long *tg, long *mt, long *jr);
};

class unsatzonmodellklasse:
public teilmodellklasse
{
public:

	// variables for coupling the model with a groundwater model 
	double **GWTableExtern,**EXTBALANCE,**FILLPOR; // gw table as imported from the external groundwater model, a balance grid and the fillable porosity also imported
	int gwtableexternindex,balanceindex; 
	bool m_bCoupled;
	bool m_useSurfaceRouting;
	bool m_bAddQdsuToP;

	bool dumpmode;
	bool closeAfterDump;
	bool repairAfterDump;
	long maxDumpFiles;
	long actDumpFiles;

	/* variables for the storage.ftz file */
	FILE *storagefile; /* a text file containing the water content of all flow time zones and 
					the actual content of the single linear reservoirs for all subbasins 
					for both, surface runoff and interflow */
	char *storagefilename;
	
	//  int fpeflag; already declared in teilmodellklasse
	//  double *balance,*maxbalance,deltapond,deltasb;

	int  pond_use,drain_use,clay_use,surface_interaction,gwflow_active,lake_active,lake_ini_method; /* flags for internal use -> controling branching of program flow */
	glacierlisttype *glacierlist;
	riverflow_table_type *riverflow_table; /* for communication between unsaturated zone and rivers -> infiltration */

	buffertype *qubuffer,*thetabuffer,*hhydbuffer,*otherdatabuffer,*etrdetailbuffer,*interflowbuffer,**concbuffer; /* concentrations are written for up to 9 tracer species */
	char *additionals; /* for display output */
	unsigned long  qubufmax,    qubufcnt,
	thetabufmax, thetabufcnt,
	hhydbufmax,  hhydbufcnt,
	otherdatabufmax,  otherdatabufcnt,
	etrdetailbufmax,etrdetailbufcnt,
	interflowbufmax,interflowbufcnt,
	*concbufmax,  *concbufcnt;

	FILE *logdatei;
	FILE *thetadatei, *qudatei, *hhyddatei, *otherdatadatei, *etrdetaildatei, *interflowdatei, **concdatei, *tempdatei;
	char *thetadateiname,*qudateiname,*hhyddateiname,*otherdatadateiname,**concdateiname,*etrdetaildateiname,*interflowdateiname;
	int control_row,control_col;

	SOILtabletype *ST;
	/* Grids */
	double  **P, **SOIL,      **ICE,       **ICEOLD,   **THAWDEPTH, // the THAWDEPTH is calculated for perma frost soils. if not present, then this option was not set in the control file (by commenting the alpha grid there or not reading it in, alpha is the parameter for calculating the thaw depth)
	**FZT,       **SNOW,      **IRRIG,     **HDRAIN,
	**DRAINLEVEL,**DRAINDIST, **SLP,       **GWSTAND,   **GWN,       **TANBETA, **UPRISE, **PERCOL,
	**QDRAIN,    **QDgrid,    **QIgrid,    **DELTAGW,   **oldDELTAGW,**GWLEVEL,
	**GWtheta,   **INFEX,     **INFILT,    **RIVDEP,    **RIVWIT,    **RIVLNK,
	**MAXPOND,   **POND,      **CLAY,      **SB_1_,      **SB_2_,    **ROOT,  **SBLOW,     **SBHIGH, **ETR2ETP,
	**IRRIGCODE, **MAXPUMP,   **GWIN,      **GWEX,      **KOLM,      **S01,	**MACRO,    **THROOT,
	**SUBSTEPS,  **QDSU, // QDSU is the direct runof fraction from silting up model and is used only when silting up model runs
	**LAKE,      **S0; // LAKE: unique codes for each lake, used for setting some defaults (DRAIN=0, SLOPE=0, GWSTAND=0) 
	// as those cells will be handled as lakes with a saturated soil column and without lateral outflows 
	// like interflow, drainage or baseflow



	// this are vectors of grids for the seperate vegetation layers
	double  ***LUSE,***ETP,***ETR,***VBG,***WURZ,***LAI,***RSC,***RSI,***Z0;
	double **ALBEDO,**RSE;
	double **EVAP,**EVAR; // potential and real evaporation from soil
	int iNumVegetationLayers;

	double **HEAD1; // is taken from groudnwater model to synchronize the exact groundwater table between both sub models (this is used when there are numerical oscillations in gwmodel, 
		// so the use of DELTAGW for synchronizing the grouindwater table between gw model and richards model is balance-safe, but the use of the real HEAD is mathematically safe 
		// (and the balance may even be hold when applied for longer times when the groundwater table crosses layer borders several times in both directions, so plus and minus balance errors cancel each other out)
	bool headinitialized;

	bool mpiinitialized;


	/* alias fuer werte der oben deklarierten Grids */
	//  double  gw_depth,p,wurz,vbg,evap,evar,etp,etr,q_u,draindist,deltagw,maxrootdepth;
	double k_clay,ini_gw_level; // parameters unique for the model domain -> no parallelization required

	double reduced_ksat; // when using frozen soils in permafrost model, this value defines the ksat for frozen soil layers
	// to which the real ksat is reduced.

	double max_Error,max_ziter; // global values for maximum iteration error and maximum number of iterations

	/*
double *a_etr,
		*a_evar,
		*etrdetail,
		*a_dicke,
		*theta_temp,
		*h_temp1,
		*theta_average,
		*theta_sl_old, // fuer Speicherung des alten theta-Wertes und der Zwischenwerte in der Iteration im finite differences scheme 
		*q_temp,
		*extractionweight,
		*qlat,
		*gwin,
		*gwex,
		*interflow,
		*lateralinflow,
		*thetasat,
		*thetares,
		maxlateralflow, // fluxes as vertical 1-d-arrays (to be summerized to get total fluxes for e.g. interflow) 
*/         
	//  double old_gwstand,old_GWtheta;

	//  int    clay_layer,old_start_level,maxlevel; 
	double  ***theta,  // 3-D-Matrize mit Bodenwassergehalten pro Schicht und Gridzelle 
	**theta1,  // grid of the first layer (for statistics) 
	//*theta_sl, // Substitutionsvariable fuer theta[i][j][1..levels] 
	//*hgeo_sl,  // Substitutionsvariable fuer hgeo[i][j][1..levels] 
	//*qu_sl,    // Substitutionsvariable fuer qu[i][j][1..levels] 
	// *hh_sl,    // Substitutionsvariable fuer hh[i][j][1..levels] 
	//**conc_sl,		// substitution variable for a single column of ***conc -> conc_sl=conc[i][j] -> for each tracer a vector 
	***hh,     // 3-D-Matrize mit hydraulischer Hoehe (Potential) jeder Schichtmitte 
	**hh1,     // grid of first layer of ***hh 
	***hgeo,   // geodaetische Hoehe der Schichtmitten -> 3D-Matrix
	**hgeo1,   // grid of first layer of ***hgeo 
	***qu,     // 3-D-Matrize mit Fluessen aus den Schichten iu die folgende Schicht unterhalb (+) beziehungsweise aus dieser in die aktuelle Schicht (-) 
	**qu1,     // grid of first layer of qu (for statistics) 
	**qu1_ext, // complete qu1 grid (mpi rank0) to be used for lake model in runoff routing (needs to be allocated when lake model is on)
	// for tracers, there are max. 9 species arranged in arrays, thats why the pointers have one indirection in addition 
	****conc,   // concentration of a tracer (salt, O18...) in each soil layer -> a vector of stacks 
	***conc1,   // grid of first layer of conc (for statistics) a vector of grids 
	//         **conc_sl_old, // for storing old tracer concentrations 
	**gw_thick,    // thickenss of the first qauifer in m from soil surface to aquifer bottom 
	***gw_conc,    // mixed concentrations for the first layer from Groundwater module 
	***CONC_QD,	// vector of grids, containing tracer concentrations in direct discharge (only internal, will not be written)
	***CONC_QI;	// vector of grids, containing tracer concentrations in interflow (only internal, will not be written)

	double QDsum;
	double count;

	zonalstatistic	
	*iteration_stat,    /* statistic for iteration counts */
	*QDstat,
	*QIstat,
	*QBstat,          /* Statistikstreams fuer die drei ABflusskomponenten */
	*Qgstat,          /* dito fuer Gesamtabfluss */
	*SB_1_stat,       /* fuer Bodenfeuchte 0..0.5 m */
	*SB_2_stat,       /* fuer Bodenfeuchte 0..1 m */
	*SB_2_IMMO_stat,       /* fuer Bodenfeuchte 0..1 m */
	*ROOTstat,       /* Wurzeltiefe */
	*QDRAINstat,      /* Drainageabfluss */
	*GWSTANDstat,     /* Grundwasserstand */
	*GWNstat,         /* Grundwasserneubildung */
	*GWINstat,        /* infiltration from rivers into unsat. zone or groundwater */
	*GWEXstat,        /* exfiltration from groundwater into surface water channels */
	*MACROstat,       // statistics for infoltration into macro pores
	*INFILTstat,      // statistics for total infiltration into the first soil layer
	*PONDstat,        /* ponding water storage content */
	*SUBSTEPSstat;	 // ststistics for the number of substeps for each timestep in the richards-model

	int  QDstatcode,    QIstatcode,     QBstatcode, Qgstatcode,
	GWINstatcode,  GWEXstatcode,
	QDRAINstatcode,GWSTANDstatcode,GWNstatcode,PONDstatcode,
	SB_1_statcode, SB_2_statcode,SB_2_IMMO_statcode, ROOTstatcode, iteration_statcode, MACROstatcode,INFILTstatcode,
	SUBSTEPSstatcode; /* Ausgabecodes fuer Statistikfiles */

	char *iteration_statfilename,
	*QDstatfilename,      /* Teilgebietsmittelwert des Direktabflusses              */
	*QIstatfilename,      /* Teilgebietsmittelwert des Interflows                   */
	*QBstatfilename,      /* Teilgebietsmittelwert des Basisabflusses               */
	*Qgstatfilename,      /* Teilgebietsmittelwert des Gesamtabflusses              */
	*GWINstatfilename,    /* Teilgebietsmittelwert der re-infiltration aus Fluessen */
	*GWEXstatfilename,    /* Teilgebietsmittelwert der exfiltration (Basisabfluss)  */
	*QDRAINstatfilename,  /* Drainageabfluss */
	*GWSTANDstatfilename, /* Grundwasserstand*/
	*GWNstatfilename,     /* Grundwasserneubildung */
	*PONDstatfilename,     /* Grundwasserneubildung */
	*ROOTstatfilename,   /* Wurzeltiefe */
	*MACROstatfilename,  // macropore infiltration
	*INFILTstatfilename,  // macropore infiltration
	*SB_2_statfilename,   /* Bodenfeuchte 0..1 m */
	*SB_2_IMMO_statfilename,   /* Bodenfeuchte 0..1 m released or immobilised by glacier growth or shrinking */
	*SB_1_statfilename,   /* Bodenfeuchte 0..0.5 m */
	*SUBSTEPSstatfilename; // number of substeps for each timestep as areal average value for each timestep

	double *maxzeit; // no parallelization, will be used in fzt_grid_auswerten and in ReadWrite SLR-Content only

	/* einzulesende Namen fuer Grids und Stacks */
	char *theta_stackname,
	*hh_stackname,
	*hgeo_stackname,
	*qu_stackname,
	*GWSTANDgridname,
	*GWthetagridname,
	*GWNgridname,
	*UPRISEgridname,
	*PERCOLgridname,
	*PONDgridname,
	*GWLEVELgridname,
	*QDRAINgridname,
	*QDgridname,
	*QIgridname,
	*DELTAGWname,
	*GWINgridname,
	*GWEXgridname,
	*SB_1_gridname,
	*SB_2_gridname,
	*MACROgridname,
	*IRRIGgridname,
	*ROOTgridname,
	*SUBSTEPSgridname;

	/* Schreibcodes fuer Grids und Stacks */
	int  theta_writestack,    theta_index,
	hh_writestack,       hh_index,
	hgeo_writestack,     hgeo_index,
	qu_writestack,       qu_index,qu1extindex,
	GWSTANDwritegrid,    gwstandindex,
	GWNwritegrid,        gwnindex,
	UPRISEwritegrid,     upriseindex,
	PERCOLwritegrid,     percolindex,
	PONDwritegrid,       pondindex,
	GWLEVELwritegrid,    gwlevelindex,
	QDRAINwritegrid,     qdrainindex,
	QDwritegrid,         qdindex,
	QIwritegrid,         qiindex,
	DELTAGWwritegrid,    deltagwindex,  olddeltagwindex,
	GWINwritegrid,       gwinindex,
	GWEXwritegrid,       gwexindex,
	GWthetawritegrid,    gwthetaindex,
	DRAINLEVELwritegrid, drainlevelindex,
	tanbetaindex,
	SB_1_writegrid,		sb_1_index,
	SB_2_writegrid,		sb_2_index,
	THROOTwritegrid,		throotindex,
	ROOTwritegrid,		rootindex,
	MACROwritegrid,      macroindex,
	SUBSTEPSwritegrid,   substepsindex,

	/* some indices for grids for irrigation only */
	maxpumpindex,
	sblowindex,
	sbhighindex,
	threduindex,
	IRRIGwritegrid, irrigindex, etr2etpindex;
	int *conc_qdindex, // indices for each element in CONC_QD, a vector of grids, containing tracer concentrations in direct discharge (internal use only, will not be written)
	*conc_qiindex; // indices for each element in CONC_QD, a vector of grids, containing tracer concentrations in interflow (internal use only, will not be written)


	statrecord ***zoneFZT_QD,
	***zoneFZT_QI;  /* fuer die Translation des Oberflaechenabflusses wird
							pro Teilgebiet pro Zeitflaechenzone ein Statistikregistersatz
							benoetigt, ausserdem ein Array mit der Information ueber die
							maximale Anzahl an Registern in jedem Teilgebiet (maxzeit) */


	//  int   TG;
	//  double qlattemp,SnowRunoff;
	double dParameterMinStep;   // for dynamic time steps this value gives the minimum timestep in secounds
	/* Arrays fuer pro TG einzulesende Parameter */
	double *code,                /* Nummern der Teilgebiete (werden in der STeuerdatei fuer die Parameterlisten zu Kontrollzwecken eingelesen) */
	*flussdichte,         /* Flussdichte zur Berechnung des Interflows */
	*SR_Anteil,           /* Anteil der Schneeschmelze, der direkt abfliesst (bzw. bei Regen auf Schnee:
								Anteil des gesamten Niederschlages + Schmelze, der direkt abfleisst im Schnee)) */
	*kelsqi,              /* kelsqi: ELS-Konstante (h) fuer Interflow */
	//        *conc_lat,*conc_gwin,*conc_gwex,
	**QD,**QDd, /* this arrays are 2-d-arrays (each tracer x subbasin) direct flow */
	**QI,**QIi, /* this arrays are 2-d-arrays (each tracer x subbasin) interflow */
	**QB,**QBb, /* this arrays are 2-d-arrays (each tracer x subbasin) baseflow */
	**Qg,       /* total flow and tracer concentration for each tracer and each subbasin */
	*kelsqd,     /* dito fuer Direktabfluss              */
	qb,*kelsqlat,  /* dito fuer Basisabfluss               */
	*Q0;             /* Gesamtabfluss,maximale Basisabfluss  */
	

	statrecord *zoneQD,         /* Statistik fuer Direktabfluss         */
	*zoneQI,         /* Statistik fuer Zwischenabfluss       */
	*zoneQB,         /* Statistik fuer Basisabfluss          */
	*zoneQg,         /* Statistik fuer Gesamtabfluss         */
	*zoneGWIN,       /* Statistik fuer reinfiltration aus Fluessen */
	*zoneGWEX,       /* Statistik fuer exfiltration aus GW (Baseflow) */
	*zoneQDRAIN,     /* Statistik fuer Drainage              */
	*zoneGWSTAND,    /* Statistik fuer Grundwasserstand      */
	*zoneGWN,        /* Statistik fuer Grundwasserneubildung */
	*zonePOND,       /* Statistik fuer Inhalt des ponding water storages */
	*zoneSB_1_,      /* Statistik fuer Bodenfeuchte 0..0.5 m */
	*zoneSB_2_,      /* Statistik fuer Bodenfeuchte 0..1 m   */
	*zoneSB_2_IMMO,      /* Statistik fuer immobilised or mobilised (by glacier resizing) Bodenfeuchte 0..1 m   */
	*zoneROOT,      /* Statistik fuer Wurzeltiefe   */
	*zoneINFEX,      /* Statistik fuer Infiltrationsexzess    */
	*zoneINFILT,     /* Statistik fuer totale Infiltration in erste Bodenschicht */
	*zoneMACRO,      /* Statistik fuer Macroporeninfiltration */
	*zoneSUBSTEPS;   // Statistics for sub-steps per time steps for dynamic time steps in method 3
	

	int methode,levels; // nur Lesezugriff, kann global bleiben
	//  long weight;
	//  double DT,timestep; 
	//  double m_err,max_Error;
	//  double *dicken; // nur Lesezugriff, kann global bleiben
	//  int start_level;

	/* variables for concentration mixing, used in function void new_concentrations() */
	/*
double *maxfromlr,  // what amount of water could flow out of the layer itself 
		*qufromlr,   // what amount of water could flow out of the layer as downward flux (given by ratio*qu_sl[k]) 
		*iflfromlr,  // what amount of water could flow out of the layer as interflow (given by ration*interflow[k]) 
		*qlatfromlr, // what amount of water could flow out of the layer as lateral flow (given by ratio*qlat[k]) 
		*gwexfromlr, // what amount of water could flow out of the layer as exfiltration (given by ratio*gwex[k])
		*ratio;      // ratio of maxfromlr/(qu_ls[k]+interflow[k]+gwex[k]+qlat[k]) 
double *macrocapacity; // capacity of the macropores in any level
*/
	//  double content_upper,content_lower; /* masses and water volumes in all saturated layers (and the first groundwater layer) for concentration mixing with groundwater */


	/* nun folgen die Variablen fuer das Infiltrationsmodell, dass hier voll in das Bodenmodell integriert ist */

	int        INFEXstatcode;
	char      *INFEXstatfilename;
	zonalstatistic   *INFEXstat;
	char      *INFEXgridname,*SATTgridname;
	int        INFEXwritegrid,SATTwritegrid;
	int        infexindex,sattindex,infiltindex;
	double    **SATT;       /* ist auf 1, wenn das Element am Ende des Zeitintervalles gesaettigt war,
							wenn im neuen Zeitintervall dann PI>Ks, dann wird weiterhin Dirketabfluss gebildet */
	//  double    Fsatt;        /* speichert die kumulative Infiltration waehrend des Saettigungsphase. */
	//  double    F;            /* gesamte Infiltration nach Saettigung. Fsatt und F addiert ergeben die gesamte
	//                            Infiltrationsmene in einem Zeitintervall. Diese Summe wird ins **P-Grid geschrieben, vorher die
	//                            Differenz **P-(F-Fsatt) ins **INFEX-Grid. */
	//  double ts;              /* Saettigungszeit, wird pro Pixel berechnet */
	//  double saug;            /* Saugspannung an der Feuchtefront (gut: 100mm) */
	//  double sicker;          /* Anteil an vom Infiltrationsueberschuss doch in den Boden versickernden Wassers */

	/* Ende Infiltrationsmodell */

	// variables for temperature model
	int do_TMod; // 0 = do not run, 1=run
	int iTModMethod; // method for temperature model: 0 = off, 1 = vertical heat transfer (diffusion only), 2 = diffusion and advection (energy transport by water flow)
	double dTModLowBoundBase; // base temperature for lower boundary condition (in °C, with regard to sea level)
	double dTModLapseRate; // lapse rate for calculating the lower boundary condition according to the elevation of each cell
	bool do_initializeTemper;
	// default soil "constants" --> soiltable values will be set if they are assigned default values of <0
	double dTModCDrySoil; // default heat capacity of dry soil in J/(Kg*K), default 800
	double dTModRhoDrySoil; // default density of dry soil in Kg/m^3 , default 1500
	double dTModLamDrySoil; // default thermal conductivity for dry soil in J/(m*s*K), default 0.58
	double dTModKmin; // minimum hydraulic conductivity for frozen soils in m/s
	double dTModLamWat; // thermal conductivity of liquid water (0.5562)
	double dTModLamIce; // thermal conductivity of ice (2.33)
	double dTModCWater; // heat capacity of water in J/(Kg*K) (4187)
	double dTModCIce20; // heat capacity of ice in J/(Kg*K) at -20°C (1940)
	double dTModCIce0; // heat capacity of ice in J/(Kg*K) at 0°C (2090)
	double dTModLatHeat; // latent heat of freezing in J/Kg (334000)
	double dTModRhoWater; // density of water in Kg/m^3 (1000)
	double dTModFactor; // Clapeyron solution dH/Tm in m/K or J/(m^3*K) (1.22)
	double dTModSEForStat; // relative saturation where we asusme ice (for thaw depth statistics)
	double dTModShortestStep; // minimum sub time step allowed for heat transfer model (if the required time step would be shorter, numeric errors like extrem temperature fluctuations are possible)
	double dTModLongestStep;  // maximum sub time step allowed for heat transfer model (to avoid instabilites induced by the nonlinearity of the processes)
	double dTModNFactorFreezing; // factor to multiply the radiation with to get the real surface temperature during negative tempertures
	double dTModNFactorThawing; // factor to multiply the radiation with to get the real surface temperature during positive tempertures
//	double actNFactor; // value from Landuse table or from the above two defaults (if not specified in landuse table)
	int useHeatFluxAsLowerBC; // is set to 1 if a grid with ID HeatFlux_Lower_Boundary_Condition is found
	
	buffertype *temperbuffer; // buffer for writing level temperatures of the control location
	unsigned long  temperbufmax,    temperbufcnt; // buffer indices and maximum counter
	FILE *temperfile; // file for output of temperatures per layer
	char *temperfilename; // file name for temperature file

	double  ***temper;  // 3-D-Matrix mit temperatures per layer and grid cell
	double **temper1; // dummy grid for first layer (for initialization)
	char *temper_stackname;
	int  temper_writestack;
	int temperindex; // index into the gridinfo-array

	double  ***SEstack;  // 3-D-Matrix with relative saturation of ice per layer and grid cell
	double **SE1; // dummy grid for first layer (for initialization)
	char *se_stackname;
	int  se_writestack;
	int seindex; // index into the gridinfo-array

	double **KX,**KY,**FROZENDEPTH; // hydraulic conductivities from groundwater: will be set to an effective value of dTModKmin if the layer is frozen (effective transmissivity will be calculated using each layers temperature)
	int frozendepthingwindex;
	// this means, the k-values are controlled by the unsatzon model and thus by the soil table rather than by the groundwater input (later, I may want to store the original hydraulic conductivities in temporary
	// grids to set it back to normal values when the ground is thawed again
	double **GW1COLMFROZEN, **GW1COLMUNFROZEN; // references to two standard grids for colmation between first and second groundwater layer, one for unfrozen and one for frozen conditions (lowest soil layer of unsaturated zone is either frozen or not)
	double **GW1COLMATION; // the original colmation grid (will be overwritten when using the above  tow grids)

	double **T; // air temperature as upper boudary condition
	int tlowbndindex; // index into the gridinfo array
	double **TLOWBND; // lower boundary condition as standard grid in _T_Lower_Boundary_Condition_
	double **QLOWBND; // lower boundary condition as constant heat flux as standard grid in _HeatFlux_Lower_Boundary_Condition_
	double **THAWTMOD; // thaw depth for physical heat transfer model (in contrast to the old conceptual THAWDEPTH)
	int thawtmodindex;

	zonalstatistic	*thawdepth_stat;    // statistic for soil temperature (Thaw depth)
	char *thawdepth_gridname;
	int  thawdepth_writegrid;
	int  thawdepthstatcode;
	char *thawdepth_statfilename;
	statrecord *zonethawdepth;

	double **FDD,**TDD,**TDAY,**TYEAR; // Freezing degree days and thawing degree days for initialization of the active layer bottom temperature
	int fddindex,tddindex,tdayindex,tyearindex;
	int numIntPreSpinUpAll; // counts intervals for pre-spinup to calculate annual average temperature
	int numIntPreSpinUp; // counts intervals for pre-spinup to calculate daily average temperature

	bool interceptionBeforeSnow;

	////////////////////////////////////////////////////////////////////////////////////////
	// variables for LayeredSnowModel, in addition to the heat transfer parameters
	////////////////////////////////////////////////////////////////////////////////////////

	int  do_LayeredSnow;
	bool do_initializeTemperSnow;
	int iLayeredSnowMethod; // method of dealing with snow layers : 1 = constant thickness, 2 = constant density(still under construction)
	int levelsSnow;
	double dLayeredSnowMinDensity; // minimum density of freshly fallen dry snow (will later be modeled in more detail, i.e. depending on temperatur and wind or so)
	double dLayeredSnowMaxDensity; // maximum density of old snow (will later be modeled dynamically depending on snow age, temperature and water content etc.)
	double dTOffsetFreshSnow;
	double dDrainageDensitiySnow;
	double dWaterHoldingCapacity;
	double dSEeffectiveness;
	double dSEeffectivenessSoil;
	SOILtabletype *LayeredSnowProperties;

	// stacks: SWE solid
	double  ***swes;  // 3-D-Matrix snow water equivalent solid (layered snow model)
	double **swes1; // dummy grid for first layer (for initialization)
	char *swes_stackname;
	int  swes_writestack;
	int swesindex; // index into the gridinfo-array

	// stacks: SWE liquid (water content in snow layers)
	double  ***swel;  // 3-D-Matrix snow water equivalent liquid (layered snow model)
	double **swel1; // dummy grid for first layer (for initialization)
	char *swel_stackname;
	int  swel_writestack;
	int swelindex; // index into the gridinfo-array

	// stacks: temperatures in snow stack
	double  ***temperSnow;  // temperatures in each snow layer (3-D-matrix)
	double **temperSnow1; // dummy grid for first layer (for initialization)
	char *temperSnow_stackname;
	int  temperSnow_writestack;
	int temperSnowindex; // index into the gridinfo-array

	// stacks: theta in snow stack
	double  ***thetaSnow;  // theta (relative volumetric water content) in each snow layer (3-D-matrix)
	double **thetaSnow1; // dummy grid for first layer (for initialization)
	char *thetaSnow_stackname;
	int  thetaSnow_writestack;
	int thetaSnowindex; // index into the gridinfo-array

	// stacks: hydraulic head in snow stack
	double  ***hhydSnow;  // hydraulic head in m above sea level for each snow layer (center) (3-D-matrix)
	double **hhydSnow1; // dummy grid for first layer (for initialization)
	char *hhydSnow_stackname;
	int  hhydSnow_writestack;
	int hhydSnowindex; // index into the gridinfo-array

	// stacks: geodetic head in snow stack
	double  ***hgeoSnow;  // geodetic head in m above sea level for each snow layer bottom (3-D-matrix)
	double **hgeoSnow1; // dummy grid for first layer (for initialization)
	char *hgeoSnow_stackname;
	int  hgeoSnow_writestack;
	int hgeoSnowindex; // index into the gridinfo-array

	// stacks: bottom flow from each layer in the snow stack
	double  ***quSnow;  // bottom flow from each layer (3-D-matrix)
	double **quSnow1; // dummy grid for first layer (for initialization)
	char *quSnow_stackname;
	int  quSnow_writestack;
	int quSnowindex; // index into the gridinfo-array

	// stacks: bottom flow from each layer in the snow stack
	double  ***thickSnow;  // bottom flow from each layer (3-D-matrix)
	double **thickSnow1; // dummy grid for first layer (for initialization)
	char *thickSnow_stackname;
	int  thickSnow_writestack;
	int thickSnowindex; // index into the gridinfo-array

	// stacks: bottom flow from each layer in the snow stack
	double  ***densSnow;  // bottom flow from each layer (3-D-matrix)
	double **densSnow1; // dummy grid for first layer (for initialization)
	char *densSnow_stackname;
	int  densSnow_writestack;
	int densSnowindex; // index into the gridinfo-array

	// a helper grid that holds the index of the uppermnost layer that contains some snow/water
	double **upperIndexSnow; // dummy grid for first layer (for initialization)
	char *upperIndexSnow_name;
	int  upperIndexSnow_writegrid;
	int upperIndexSnowindex; // index into the gridinfo-array

	// grids that will be initialized in snow model and need to be referenced in unsatzonmodel only
	double **SURFTEMP; // snow surface temperature
	double **EBAL; // energy balance
	double **SURFTEMPGL; // snow surface temperature on glacierized fractions
	double **EBALGL; // energy balance snow on glacierized fractions
	double **P_RAIN; // liquid fraction of precipitation in mm/time-step
	double **SSNO; // pointer to the snow pack from snow model
	double **SLIQ; // pointer to the water in snow pack from snow model
	double **SSTO; // total snow storage (un-layered snow model)
	double **QSNOW_OUT; // returns the effective snow melt (change in SWES) from unsatzon model to the snow model

//////////////////////////////////////////////////////////////////////////////////
// variables for ordered processing of re-infiltration, especially in MPI mode ///
//////////////////////////////////////////////////////////////////////////////////

	std::vector<RiverLink*> *RiverLinkVector;
	std::vector< std::vector<RiverLinkCell> > *RivLnkCellsToSend;    // Array for each instance that needs to send cells to the processing rank. first dimension: remote processing rank, second dimension: single cells
	std::vector< std::vector<RiverLinkCell> > *RivLnkCellsToReceive; // Array for each processing instance. First dimension: rank to receive from; second dimension: single cells in processing order
	unsigned long numAllocatedDoubles;
	double **tempGridRows;
	double ***tempStackRows;
	unsigned int numTempGridRows;
	unsigned int numTempStackRows;
	firstlast tempFirstLastCol;

	double **RecvBuffer;
	unsigned long *RecvBufferSize;
	unsigned long *RecvBufferPos; // used in processing ranks. Holds for each potential rank the position into the buffer (i.e. RecvBuffer index)

	double **SendBuffer;
	unsigned long *SendBufferSize;

	//////////////////////////////////////////////////////////////////////////////////
	// variables for integration of ArcVeg ///
	//////////////////////////////////////////////////////////////////////////////////

	int do_ArcVeg;
	bool bArcVegActive; // set internally to true when start JD for arcVeg is reached and end Julian Day is not yet reached, i.e. JD between start and end JD -> or on southern hemisphere if outside the range between start end end JD
	int ArcVegStartJD, ArcVegEndJD, ArcVegWaitPeriod;
	bool bArcVegWaiting,bArcVegSnowFree;

	// grids that holds the depth of the upper ice wedge elevation (in m below soil surface)
	double **ICEWEDGEUPPERDEPTH; // 
	// char *iceWedgeUpperDepth_name;
	// int  iceWedgeUpperDepth_writegrid;
	// int iceWedgeUpperDepthindex; // index into the gridinfo-array

	// grids that holds the depth of the lower ice wedge elevation (in m below soil surface)
	double **ICEWEDGELOWERDEPTH; // 
	// char *iceWedgeLowerDepth_name;
	// int  iceWedgeLowerDepth_writegrid;
	// int iceWedgeLowerDepthindex; // index into the gridinfo-array

	// Growing degree days for ArcVeg
	double **AVGDD; 
	char *AVGDD_name;
	int  AVGDD_writegrid;
	int AVGDDindex; // index into the gridinfo-array
	// grids that holds maximum value of the 7-day-moving average of the inundation depth (in m)
	double **AVINMAX7DAVG; // maximum 7-day-moving average of inundation after start Julian day and waiting period until end julian day fpr ArcVeg
	char *AVINMAX7DAVG_name;
	int  AVINMAX7DAVG_writegrid;
	int AVINMAX7DAVGindex; // index into the gridinfo-array

	double **AVTDAVG; // ArcVeg daily average Temperature grid
	char *AVTDAVG_name;
	int  AVTDAVG_writegrid;
	int AVTDAVGindex; // index into the gridinfo-array

	double **AVINDAVG; // ArcVeg daily average Inundation grid
	char *AVINDAVG_name;
	int  AVINDAVG_writegrid;
	int AVINDAVGindex; // index into the gridinfo-array

	double** AVSNOWFREE; // ArcVeg Snow Free Days grid
	char* AVSNOWFREE_name;
	int  AVSNOWFREE_writegrid;
	int AVSNOWFREEindex; // index into the gridinfo-array

	double** AVPGP; // ArcVeg number of plant groing periods
	char* AVPGP_name;
	int  AVPGP_writegrid;
	int AVPGPindex; // index into the gridinfo-array

	// stack for internal inundation scaling when calculating actual changes once a year (value per land use fraction)
	double*** AVINSCAL; // ArcVeg Inundation Scaling
	double** AVINSCAL1; // ArcVeg Inundation Scaling
	char* AVINSCAL_stackname;
	int  AVINSCAL_writestack;
	int AVINSCALindex; // index into the gridinfo-array

	// stack for computing the moving average of the inundation depth (in m) over 7 days
	double  ***AVINHLP;  // a stack that holds 7 day values for inundation and a running average in level 0
	double **AVINHLP1; // dummy grid for first layer (for initialization)
	char *AVINHLP_stackname;
	int  AVINHLP_writestack;
	int AVINHLPindex; // index into the gridinfo-array
	
	// arcveg total soil organic mineralization 
	double** AVSON; // ArcVeg soil organic nitrogen
	char* AVSON_name;
	int  AVSON_writegrid;
	int AVSONindex; // index into the gridinfo-array

	// stack for nitrogen mineralization per plant growth period
	double*** AVNMIN; // ArcVeg nitrogen mineralization stack
	double** AVNMIN1; // ArcVeg nitrogen mineralization grid 
	char* AVNMIN_stackname;
	int  AVNMIN_writestack;
	int AVNMINindex; // index into the gridinfo-array

	// arcveg mineralization total
	double** AVNMINTOT; // ArcVeg soil organic nitrogen
	char* AVNMINTOT_name;
	int  AVNMINTOT_writegrid;
	int AVNMINTOTindex; // index into the gridinfo-array

	int AVNumClimZon; // number of climate zones for which N-mineralization rates for each PGP must be read in 
	int AVNumPGPs; // number of plant growing periods (PGPs) for which N-mineralization parameters must be read in
	int AVActClimZon; // actual climate zone used for the actual model setup
	double** AVNMinParam; // 2-D-Array for N-mineralization. first dimension is climate zones count, second is PGPs count. 1-based with an empty 0-element 

	// variables that are needed for translating the biomass change into soil thickness change
	int AVLUSEType; // land use type that will be used for changing the soil thickness (i.e. when the biomass of moss is changing, this soil type represents the moss)
	int AVLayerCnt; // COunt of soil layers that should be changed in thickenss when arcveg runs
	double AVSoilDryDensity; // density of dry MOSS (or whatever land use type is used for thickness change)

	// stack for biomass for all regarded plant types 
//	double*** AVBIOM;  // a stack that holds 7 day values for inundation and a running average in level 0
//	double** AVBIOM1; // dummy grid for first layer (for initialization)
//	char* AVBIOM_stackname;
//	int  AVBIOM_writestack;
//	int AVBIOMindex; // index into the gridinfo-array

	int arcVegNumLayers;
	int arcVegNumBiomass;
	int arcVegNumRetrans;
	double ***ARCVEGFRACTIONS; // used for keepting track of the fraction of each landuse type when ArcVeg is used: there may be multiple (usually 5) land use types per cell, these grids keep track of their fraction (created in defaultini)
	double ***ARCVEGBIOMASS; // used for keepting track of the biomass of each landuse type when ArcVeg is used: there may be multiple (usually 5) land use types per cell, these grids keep track of their biomass (created in defaultini)
	double ***ARCVEGRETRANS; // used for keepting track of the retrans amount (from senescense modul) of each landuse type when ArcVeg is used: there may be multiple (usually 5) land use types per cell, these grids keep track of their retrans amount (created in defaultini)

	// stacks: ExcessIceFactor
	double  ***ExcessIceFactor;  // a stack that holds excess ice factors for each soil layer. for initialisation, ICEWEDGEUPPERDEPTH and ICEWEDGELOWERDEPTH will be used 
	double **ExcessIceFactor1; // dummy grid for first layer (for initialization)
	char *ExcessIceFactor_stackname;
	int  ExcessIceFactor_writestack;
	int ExcessIceFactorindex; // index into the gridinfo-array

	int indexFCmin_mm; // gridindex for resiual water content for root zone 
	int indexFCmax_mm; // gridindex for usable water content at field capacity = ttoal pore water at field capacity minus residual water content
	int indexFCact_mm; // grid index for usable water content actual 
	int indexFCact_perc; // same as actual usable water content but in percent of the maximum available water content
	double **FCMIN, **FCMAX, **FCACT, **FCACTPERC; // grids for the mentioned various water contents, will be handled by standard grids in module initial, but need to be available here as references
	bool doFCcalc;

	int sb_1_rel_index; // gridindex for soil moisture as relative value on total pore volume instead as relative value on total soil volume as in SB_1
	double **SB_1_REL; // soil moisture as relative value on total pore volume instead as relative value on total soil volume as in SB_1
	bool doSB1RELcalc; // switch, if SB_1_REL was found as standard grid

///////////////////////////////////////////////////////////////
// methods ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

	unsatzonmodellklasse(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
						double breite,double laenge,double bezugsmeridianarg,const char *bezeichnung,
						char *steuername,int LogLevel,int progaus,
						double p_min_to_int_end,int numsubb,int *listsubb,
						gridinfotype *p_gridinfo,
						int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
						int p_zonmax,  int p_zonmaxwert,
						int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
						double p_xll ,  double p_yll,
						double p_cellsize    , double      p_nodata,
						long mine, long ste , long tge , long mte , long jre, double p_timekorr,
						LUSEtabletype *p_LUtable,
						SOILtabletype *p_SOILtable,
						MultiLayerLUSEType *p_MultiLayerLUTable,
						irrigationtabletype *p_irrigationtable,
						bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
						long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
						int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,
						std::vector< std::vector<RiverLinkCell> > *p_RivLnkCellsToSend,
						std::vector< std::vector<RiverLinkCell> > *p_RivLnkCellsToReceive,
						std::vector<RiverLink*> *p_RiverLinkVector, int substepcnt = 1);

	virtual void createFieldCapacityGrids(); // creating internal helper grids to keep track of field capacity and its usable fraction during the model run (mostly defined in standard_grids but handled in unsatzonmodel)
	virtual void updateFieldCapacityGrids(); // updating the usable field capacity

	virtual void createSB1RELGrid(); // checking for relative soil moisture grid and set doSB1RELcalc to true;
	virtual void updateSB1RELGrid(); // updating the filled pores relative tot total pore volume grid

	virtual void setSubStepCnt(int step); // inherited from teilmodellklasse and changed here

	virtual void   ReadFloatFromStoragefile(double *FValue,const char *ErrorMessage);
	virtual void   StorageFileError(const char *message,double c1,double c2);
	virtual void   WriteSlrContents(bool appendmode);
	virtual void   ReadSlrContents();
	virtual char   *GetStorageFileName();
	virtual double quflow(double kf,double h1,double h2,double dz1,double dz2);
	virtual double gradient_restriction(int j,int s,double thetao,double thetau,
					double qup,double qmid,double qlow,double *h,double* &interflow,double* &hh_sl,double DT,double *SE_sl,double *thetar, double *thetas);
	virtual double  GetMeanTheta(int &s,double* &dicken,double* &theta_sl,double upper=0.0,double lower=9999.0,int AsAbsoluteAmount=0);
	
	virtual int    evapotranspiration(int i,int j,double a,int s,double dReduParDry,double dMaxThetaWet,double dLimitReduWet,
					double *a_dicke,double *a_etr,double *a_evar,
					double &wurz,double &vbg,double &etp,double evap,
					double* &dicken,double* &thetasat,double* &thetares,double* &theta_sl,double* &SE_sl,
					double* &hgeo_sl,double* &etrdetail,double* &extractionweight,int gwlevel,double DT,double evapmaxdepth);
	
	virtual void   hydraulic_heads_iterative(int offset,int i,int j,int z,
					int &maxlevel,int &s,double &h_temp,
					double* &gwin,double* &gwex,double* &h_temp1,
					double* &theta_sl_old,double* &theta_average,
					double* &hh_sl,double* &hgeo_sl,double * SE_sl, double *thetar, double *thetas);
	
	virtual void specific_fluxes(int offset,int i,int j,int s,int maxlevel,double draindist,int clay_layer,double p,
					double* &interflow,double* &qu_sl,double* &gwin,double* &gwex,double* &conc_gwin,
					double* &thetares, double* &thetasat,double* &theta_sl_old,double* &theta_average,
					double* &hgeo_sl,double* &hh_sl,double* &dicken,double &qdrain,double &GWINtotal,double &GWEXtotal,double DT,int do_TMod=0,double *SE_sl=NULL);
	
	virtual void   flux_restrictions(int offset, int i, int j, int anz,
		int &maxlevel, double* &gwin, double* &gwex,
		double* &theta_sl_old, double* &theta_temp, double* &thetasat, double* &thetares,
		double* &qu_sl, double* &interflow, double* &dicken, double &GWINtotal, double &GWEXtotal, double DT, double *SE_sl);

	virtual void   flux_restrictions_opt(int offset, int i, int j, int anz,
		int &maxlevel, double* &gwin, double* &gwex,
		double* &theta_sl_old, double* &theta_temp, double* &thetasat, double* &thetares,
		double* &qu_sl, double* &interflow, double* &dicken, double &GWINtotal, double &GWEXtotal, double DT, double *SE_sl);

	virtual double new_thetas(int offset,int i,int j,int z,
					int &maxlevel,double* &gwin,double* &gwex,double* &theta_sl,
					double* &theta_average,double* &theta_sl_old,double* &qu_sl,
					double* &interflow,double* &dicken,long &weight,double DT);
	
	virtual void   calc_surface_interaction(int i,int j,int lr,
					double* &gwin,double* &gwex,double* &conc_gwin,
					double* &thetares, double* &thetasat,double* &theta_sl_old,double* &theta_average,
					int &s,double* &hgeo_sl,double* &hh_sl,double* &dicken,double &GWINtotal,double &GWEXtotal,double DT);

	virtual double  groundwater_recharge(int i,int j,double gw_exact,
					int &old_start_level,int &start_level,double &old_gwstand,
					double* &hgeo_sl,double* &theta_sl,double* &thetasat,
					double* &theta_sl_old,double* &lateralinflow,double* &dicken,double* &qu_sl,
					double &GWINtotal,double &GWEXtotal,double DT,double timestep);
	
	virtual void   infiltration(int i,int j,int &clay_layer,double &gw_depth,
					int &s,double* &dicken,double* &hgeo_sl,double* &theta_sl,double &P,double &INFEX,double timestep);
	
	virtual void   new_concentrations(int i,int j,
					int &start_level,double** &conc_sl,double** &conc_sl_old,double* &qu_sl,double* &theta_sl,double* &hgeo_sl,
					double* &lateralinflow,double* &gwin,double* &gwex,double* &theta_sl_old,double* &thetasat,double* &maxfromlr,
					double* &dicken,double* &interflow,double* &qlat,double* &ratio,double* &iflfromlr,double* &qlatfromlr,
					double* &gwexfromlr,double* &qufromlr,double* &conc_lat,double* &conc_gwex,double* &conc_gwin,
					double &QDgrid,double &QIgrid,double GWEXtotal,double DT,double timestep);

	virtual void   flushbuffer(FILE *datei,char *dateiname,unsigned long *bufcnt,buffertype *buffer,int start=1,int maxlevels=0);
	virtual void   flushbuffers();

	virtual void   get_max_pumping();
	virtual void   fzt_grid_auswerten();
//	virtual void   zonenidentification();
	virtual double  **GetQgesAdresse();
	virtual double  *GetKelsQIAdresse();
	virtual double  *GetKelsQDAdresse();
	virtual double  **GetQDAdresse();
	virtual double  **GetQIAdresse();
	virtual double  **GetQBAdresse();
	virtual riverflow_table_type *GetRiverflow_tableAdress();
	virtual void   pre_run(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void   run(long me, long se,long tg,long mt,long jr,long actsubstep);
	virtual void   post_run(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void   ClearMemory(bool appendmode);
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */

	virtual int  GetReadGrid();
	virtual void GetExternalGWTable(int iExternalTimeStep);	
	virtual double CalculateMacroporeInfiltration(double &p,int &s,double* &dicken,double* &macrocapacity,double* &thetasat,double* &theta_sl,double timestep,double* &macroinfiltration);
	virtual void initialize_lakes(std::string instanzname);
	virtual void SetSurfaceRoutingState(bool bUseSurfaceRouting);

	virtual void CalculateLambdaEff(int offset,double *lamda_sl,int soiltype,double *SE_sl,double *theta_sl, double *thetas, double LambdaWater,double LambdaIce);
	virtual void Expand_dEoverdT(int offset,double *theta_sl,double *theta_res,double *dEOverdT_sl);
	virtual void CalcHeatFluxes(int offset,double *qHeat_sl,double *lambda_sl,double *dicken,double *temper_sl,double TLOWBND,double Tair,double actNFactor,double Tsurf);
	virtual void CalcDeltaT(int offset,double *deltaT_sl, double *qHeat_sl,double *theta_sl,double *dicken,double *heatcap_sl,int s,int timestep);
	virtual void CalcNewTemperature(int offset,double *temper_sl_new,double *temper_sl, double *deltaT_sl,double *qHeat_sl,double Tbndlow,double Tbndup,bool withlimitation,double Tsurf);
	virtual void SolveTriDiag(int offset,int row,int col,double *lambda_sl,double *heatcap_sl,double *dicken,double *temper_sl,double timestep_heat,
		                      double *coeff_a,double *coeff_b,double *coeff_c,double *right_v,double *result_x,double actNFactor,
							  double *macroinfiltration,double* qu_sl,double Tsurf);

	virtual void temperaturePreSpinUp(long me, long se,long tg,long mt,long jr);
	virtual void initializeSoilTemperatures();
	virtual void setDumpMode(bool p_dumpmode, bool p_closeAfterDump, bool p_repairAfterDump, long p_maxDumpFiles) { this->dumpmode = p_dumpmode; this->closeAfterDump = p_closeAfterDump; this->repairAfterDump = p_repairAfterDump; this->maxDumpFiles = p_maxDumpFiles; };
	virtual void setInterceptionBeforeSnow(bool b_interceptionBeforeSnow);

	virtual int  readLayeredSnowModelParameters();
	virtual void readLayeredSnowProperties();
	virtual void ReadOrCreateLayeredSnowGridsAndStacks();

	virtual bool getHeatTransferActive() 	{ 		return (this->do_TMod>0); 	};
	virtual bool getLayeredSnowActive() 	{ 		return (this->do_LayeredSnow>0); 	};
	virtual SOILtabletype* getLayeredSnowProperties()	{		return this->LayeredSnowProperties;	};
	virtual double getTOffsetFreshSnow()	{		return this->dTOffsetFreshSnow;};
	virtual double updateSnowStorageLayers(double *swes_sl,double *swel_sl,double *theta_sl,double *thetasat,double *dicken,double *SE_sl,double *densSnow_sl,double *temper_sl);

	virtual unsigned long prepareMPIReinfiltration();
	virtual int allocateGridElementsInRow0(double **grid, unsigned long cnt);
	virtual int allocateStackElementsInRow0(double ***stack, unsigned long cnt,int level);
	virtual unsigned long putCellDataIntoBuffer(int i,int j,double *buffer, unsigned long startpos);
	virtual unsigned long getCellDataFromBuffer(int i, int j, double *buffer, unsigned long nextpos);
	virtual unsigned long fillSendBuffer(std::vector<RiverLinkCell> &RivLnkVector, double *buffer, unsigned long startpos);
	virtual void MPISyncRivLnkCells1(); // sending data from ranks 0 ... (mpi_size-1) to the processing ranks
	virtual void MPISyncRivLnkCells2(); // sending data from processing ranks back to the ranks 0 ... (mpi_size-1)
	virtual void SetTemporaryRowPtrs();
	virtual void RecoverOrigRowPtrs();
	virtual void MPIBCastRiverFlowTable();
	virtual void MPIGatherRiverFlowTable();
	virtual void MPIPrepareRiverFlowTable();
	virtual void run_arcveg(long me, long se, long tg, long mt, long jr);
};

class irrigationclass:
public teilmodellklasse
{
public:

	/* Grids */
	double  **P,**SOIL,**IRRIG,**SB,**SBLOW,**SBHIGH,**MAXPUMP,**IRRIGCODE,**ICE,**ETR2ETP;
	double **WURZ;
	int iNumVegetationLayers;
	double ***PUMPCONC; /* tracer concentrations within pumping water (for salts: masses!) */
	double  p,luse,soil,wurz,irrig,sb,maxpump; /* alias fuer die oben deklarierten Grids */
	double **Qges; /* array of flows and tracer concentrations for all subbasins and all tracers */
	double **maxirrig; /* array with maximum allowed irrigation water amount from surface water (calculated in routing_module at the last time step for the next (now the actual) timestep) */
	double **maxirrigFromReservoirs; /* array with maximum allowed irrigation water amount from reservoirs (calculated in routing_module at the last time step for the next (now the actual) timestep) */
	zonalstatistic  *IRRIGgwstat,*IRRIGswstat; /* statistics for groundwater and surface water irrigation amounts */

	riverflow_table_type *riverflow_table; /* for communication between unsaturated zone and rivers -> infiltration */
	int  IRRIGgwstatcode,IRRIGswstatcode;

	char *IRRIGgwstatfilename,    /* Teilgebietsmittelwert des IRRIGATIONwassers aus Grundwasser */
	*IRRIGswstatfilename;    /* Teilgebietsmittelwert des IRRIGATIONwassers aus Grundwasser */

	/* Schreibcodes fuer Grids und Stacks */
	int  TG,
	methode;
	int **rtarget; /* target of the routing subbasins for information on how much irrigation water is used and how much can be given to the target basin */
	double **abfrel; /* relative fractions of routed subbasins area on the total basin area */
	/* Arrays fuer pro TG einzulesende Parameter */
	double *code;                /* Nummern der Teilgebiete (werden in der STeuerdatei fuer die Parameterlisten zu Kontrollzwecken eingelesen) */

	statrecord *zoneGW,*zoneSW;      /* Statistik fuer Bewaesserungswasserauftrag  */

	irrigationclass(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
		double breite,double laenge,double bezugsmeridianarg,const char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int numsubb,int *listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
		int p_zonmax,  int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll ,  double p_yll,
		double p_cellsize    , double      p_nodata,
		long mine, long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		irrigationtabletype *p_irrigationtable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void   run(long me, long se,long tg,long mt,long jr);
	virtual void   SetQgesAdresse(double **p_Qges);
	virtual statrecord *GetZoneSWAdress();
	virtual void   SetRiverflow_tableAdress(riverflow_table_type *adresse);
	virtual void   SetMaxIrrigAdress(double **p_maxirrig,double **p_maxirrigFromReservoirs);
	virtual void   SetRtargetAdress(int **p_rtarget);
	virtual void   SetAbfrelAdress(double **p_abfrel);
	virtual void   GetPumpingConcentrations();
	virtual void   ClearMemory();
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void   flushbuffers();
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */

};

class gwflowclass:
public teilmodellklasse
{
public:

	/* Grids */
	double   **DHM,      /* digital elevation model */
	**GWN,      /* recharge rate in m/s perpendicular to the grid surface */
	**GWSTAND,  /* groundwater table below surface [m] */
	**SOIL   ,  /* soil types to get information about porosities */
	**DELTAGW,  /* rate of change in groundwater level, given as filter velocity in m/s perpendicular to the grid surface */
	**PONDold,  // for calculating the balance after a time step
	***BOUNDH,   /* boundary conditions (fix head) for layers 1..n                       */
	***BOUNDQ,   /* boundary conditions (fix flux q) for layers 1..n                     */
	***HEAD,     /* hydraulic heads for layers 1 .. n                                    */
	***HEADold,  /* old hydraulic head (used for interpolation)                          */
	***HEADolditer,  /* hydraulic head from the forelast iteration step (used for estimating optimal SOR-factor */
	***HEADoldold, /* for internal use: automatic calculation of best sor-factor: forelast head is needed */
	***BALANCE,  /* balance of in- and outflows with storage change (check of mass conservation and influence of boundary consitions) */
	***THICK,    /* aquifer thicknes for all layers                                      */
	***KX,       /* saturated hydraulic conductivity for all layers in x direction (m/s) */
	***KY,       /* saturated hydraulic conductivity for all layers in y direction (m/s) */
	***TRANSX,   /* Transmissivity for all 3 layers in x-direction (for layers 2 and lower, it can be used instaed of KX) */
	***TRANSY,   /* Transmissivity for all 3 layers in y-direction (for layers 2 and lower, it can be used instaed of KY) */
	***FLOWX,    /* output flow in x direction in m2/s for all layers                    */
	***FLOWY,    /* output flow in y direction in m2/s for all layers                    */
	***S0,       /* storage coefficients                                                 */
	***LEAKAGE,  /* leakage between groundwater layers */
	***UPLIM,    /* upper norder of aquifer in m a.s.l. */
	***KOLM;     /* leakage factors for all layers (with respect to the next lower layer) */

	// double precision versions for internal use
//	double ***dHEAD,     /* hydraulic heads for layers 1 .. n                                    */
//	***dHEADold,  /* old hydraulic head (used for interpolation)                          */
//	***dHEADolditer,  /* hydraulic head from the forelast iteration step (used for estimating optimal SOR-factor */
//	***dHEADoldold; /* for internal use: automatic calculation of best sor-factor: forelast head is needed */

	double  **ICE; /* codes for glaciers if such a grid is used */
	double  **POND; /* used for calculating the correct hydraulic head for cells with a pond in order to model lake growth by groundwater exfiltration */
	double  **LAKE; // indicator, that POND MUST be used, since lakes are models using PONDs
	double  **THAWDEPTH; // indicator, that permafrost soils are modelled and groundwater movement should be limited to the active layer only (for deep permafrost)
	double  **FROZENDEPTH; // thickness of the uppermost aquifer which is frozen and in water saturated conditions (in m), unsaturated zone does not count
	double  **QDSU; // surface storage from surface routing model

	buffertype *glogbuffer;
	unsigned long  glogbufmax,glogbufcnt;
	FILE *glogfile;
	char *glogfilename;
	int control_row,control_col;
	int lake_active,lake_ini_method;


	int       layers,auto_sor,s;   /* number of layers (max 3), auto_sor: boolean variable to indicate automatic estimation of an relaxation factor, s: soil type */
	int       iUsePond,			// a flag read in in the control file if pond grids should be used at all
	iPondAvailable;     // an internal flag if the pond grid is available at all. Both flags must be > 0
	long      num_oscill;        /* on how many points doesthe result oszillate */
	double    dx,dy,    /* dx,dy: cellsize in x and y direction; */
	alpha,    /* weighting of forward differences in Crank-Nicholson-approach */
	sor,      /* relaxation coefficient for successive over relaxation method standard=1.0 */
	sor_sum,sor_anz,sum_diff,maxsum_diff;

	double    maxerror,     /* actual maximum iteration error between 2 iteration steps (maximum of all cells of an entire grid)) */
	maxerrorf,    /* allowed maximum iteration error bwtween two iteration steps */
	errorsum,     /* sum of all errors in a single iterationstep of an entire grid */
	errorcount,   /* number of added errors in errorsum */
	meanerror,    /* actual mean error */
	meanerrorf;   /* allowed mean error */

	double    mass0,     /* tracer mass within the actual cell                */
	massflux1, /* massflux from/into the back cell (lokal index 1)  */
	massflux2, /* massflux from/into the right cell (lokal index 2) */
	massflux3, /* massflux from/into the front cell (lokal index 3) */
	massflux4, /* massflux from/into the left cell (lokal index 4)  */
	massflux5, /* massflux from/into the upper layer (by laekage)   */
	massflux6; /* massflux from/into the lower layer (by laekage)   */

	int       iter,         /* actual iteration step count */
	maxiter;      /* allowed number of iteration steps */

	char    **FLOWXgridname,
	**FLOWYgridname,
	**BALANCEgridname,
	**HEADgridname,
	**HEADoldgridname,
	**HEADoldoldgridname,
	**HEADolditergridname;

	int      *FLOWXwritegrid,
	*FLOWYwritegrid,
	*HEADwritegrid,
	*HEADoldwritegrid,
	*HEADoldoldwritegrid,
	*HEADolditerwritegrid,
	*BALANCEwritegrid,
	*transxindex,
	*transyindex,
	*flowxindex,
	*flowyindex,
	*leakageindex,
	*uplimindex,
	*headindex,
	*headoldindex,
	*headolditerindex,
	*headoldoldindex,
	*balanceindex;

	int s0index,
		gwstandindex,
		deltagwindex,
		pondindex,
		pondoldindex,
		thawdepthindex,
		frozendepthindex;


	double   **Qges; /* array of flows and tracer concentrations for all subbasins and all tracers */

	riverflow_table_type *riverflow_table; /* for communication between unsaturated zone and rivers -> infiltration */

	/* Schreibcodes fuer Grids und Stacks */
	int	TG;
	int	methode;

	double frozenDepthScaling;

	gwflowclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
				double breite,double laenge,double bezugsmeridianarg,const char *bezeichnung,
				char *steuername,int LogLevel,int progaus,
				double p_min_to_int_end,int numsubb,int *listsubb,
				gridinfotype *p_gridinfo,
				int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
				int p_zonmax,  int p_zonmaxwert,
				int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
				double p_xll ,  double p_yll,
				double p_cellsize    , double      p_nodata,
				long mine, long ste , long tge , long mte , long jre, double p_timekorr,
				LUSEtabletype *p_LUtable,
				SOILtabletype *p_SOILtable,
				MultiLayerLUSEType *p_MultiLayerLUTable,
				irrigationtabletype *p_irrigationtable,
				bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
				long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
				int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int substepcnt=1);

	virtual void setSubStepCnt(int step); // inherited from teilmodellklasse and changed here

	virtual void   initialize_lakes(std::string instanzname);
	virtual void   transmissivities();
	virtual void   flushbuffer(FILE *datei,char *dateiname,unsigned long *bufcnt,buffertype *buffer);
	virtual void   fluxfield();
	virtual double gauss_seidel_iteration(int lr);
	virtual double HorizontalEquationSystem(int lr);
	virtual double get_oscillations(double **H, double **Hdold, double **Holdold);
	virtual void   pre_run(long me, long se,long tg,long mt,long jr,long actsubstep);
	virtual void   run(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void   post_run(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void   SetQgesAdresse(double **p_Qges);
	virtual void   SetRiverflow_tableAdress(riverflow_table_type *adresse);
	virtual void   ClearMemory();
	/* overwriting the inherited funstions for substance transport from teilmodellklasse */
	virtual void   flushbuffers();
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */

};

enum type_SubStepsCumulationType {NONE = 0,ADD = 1,AVG = 2};
class SpecialOutputClass;
class SpecialOutputFileClass;

class SpecialOutputEntityClass
{
wasim::mystream *WSstream;
public:
	int numPoints;
	int gridinfoindex;
	double xll,yll,cs;
	int ncols, nrows,mpistartrow,mpiendrow;
	std::string ID;
	std::string Symbol;
	std::vector<double> Xcoord; 
	std::vector<double> Ycoord; 
	std::vector<int> Level; 
	std::vector<int> col; 
	std::vector<int> row; 
	StringUtils su;
	gridinfotype *gridinfo;
	SpecialOutputFileClass *parent;
	
	bool bUsesSubSteps;
	int iNumSubSteps;
	type_SubStepsCumulationType iSubStepsCumulationType;

	int ist_vorhanden(const char *name);

	SpecialOutputEntityClass(SpecialOutputFileClass *p_parent, wasim::mystream *Stream, std::string entityString, gridinfotype *gridinfo, double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int numsubsteps);
	~SpecialOutputEntityClass();
	int addValues(buffertype *buffer,unsigned long bufcnt,int startindex,bool bInSubStepsLoop,int actSubStep); // returns the last index written to in the buffer array
	int getNumPoints() {return numPoints;};
	int getGridinfoindex() {return gridinfoindex;};
	std::string getSymbol() {return Symbol;};
	std::vector<double> getXcoord() {return Xcoord;};
	std::vector<double> getYcoord() {return Ycoord;};
	std::vector<int> getCol() {return col;};
	std::vector<int> getRow() {return row;};
	std::vector<int> getLevel() {return Level;};
};

class SpecialOutputFileClass
{
wasim::mystream *WSstream;
public:
	int numEntities;
	buffertype *buffer;
	unsigned long bufmax,bufcnt;
	FILE *outputfile;
	std::string header;
	std::string filename;
	std::vector<SpecialOutputEntityClass*> entities;
	StringUtils su;
	gridinfotype *gridinfo;
	double xll,yll,cs;
	SpecialOutputClass *parent;
	int iNumSubSteps;
	bool incrementBufCntNextTime;
	
	SpecialOutputFileClass(SpecialOutputClass *p_parent, wasim::mystream *Stream, std::string record, gridinfotype *p_gridinfo, double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, char *DefaultOutputDirectory, int substepcnt);
	~SpecialOutputFileClass();
	void ExtendFileName(char *directory,char *dateiname);
	void flushbuffer();
	void addValues(long me, long se,long tg,long mt,long jr,bool bInSubStepsLoop,int actSubStep);
};

class SpecialOutputClass:
public teilmodellklasse
{
public:
	std::vector<SpecialOutputFileClass*> files;

	int numfiles;
	int iNumSubSteps;
	StringUtils su;

	void flushbuffers();
	void ClearMemory();
	void stripInt(std::string &row,std::string sWhatParam,int *param);

	SpecialOutputClass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
						double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
						char *steuername,int LogLevel,int progaus,
						double p_min_to_int_end,int p_numsubb,int *p_listsubb,
						gridinfotype *p_gridinfo,
						int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
						int p_zonmax,  int p_zonmaxwert,
						int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
						double p_xll ,  double p_yll,
						double p_cellsize    , double      p_nodata,
						long mine,long ste , long tge , long mte , long jre, double p_timekorr,
						LUSEtabletype *p_LUtable,
						SOILtabletype *p_SOILtable,
						MultiLayerLUSEType *p_MultiLayerLUTable,
						bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
						long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
						int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int substepcnt=1);

//	~SpecialOutputClass();
	void run(long me, long se,long tg,long mt,long jr,bool bInSubStepsLoop,int actSubStep);
};

typedef struct {
	int maxlayer;
	double h_temp;
	double p;
	double DELTAGW;
	double *ETR;
	double *gwin;
	double *gwex;
	double *h_temp1;
	double *theta_sl_old;
	double *theta_average;
	double *hh_sl;
	double *hgeo_sl;
	double *SE_sl;
	double *qu_sl;
	double *interflow;
	double *thetasat;
	double *dicken;
	double *theta_temp;
	double *thetares;
	double *temper_sl;
	double *theta_sl;
	double *lambda_sl;
	double *dEOverdT_sl;
	double *heatcap_sl;
	double *qHeat_sl;
	double *temper_sl_new;
	double *swes;
	double *swel;
	double *thickSnow;
	double *densSnow;
	double upperIndexSnow;
	double SSNO;
	double SLIQ;
	double SSTO;
	int s;
	double *LUSE;
	double timestep_heat;
	int stunde;
	int tag;
	int monat;
	int jahr;
	double splitfactor;
	double DT;
	int swesoffset;
	double splitfactor_heat;
	int heatmultiple;
	double T;
	double Tsurf;
} Dumptype;

typedef struct {
	double *refSWES;
	double *refSWEL;
	double *refThickSnow;
	double *refDensSnow;
	double refUpperIndexSnow;
	double EVRS;
	double SSNOW;
	double SLIQ;
	double SSTO;
	double T;
	double Albedo;
	double *LUSE;
	double DHM;
	double ICE;
	double P;
	double SSNOWGL;
	double SLIQGL;
	double SNOWTEMP;
	double SURFTEMP;
	double SNOWTEMPGL;
	double SURFTEMPGL;
	double EVRSGL;
	double EBAL;
	double EBALGL;
	double *FIRN;
	double LAKE;
	double POND;
	double RLF;
	double CL;
	double RGEX;
	double U;
	double EXKORR;
	int stunde;
	int tag;
	int monat;
	int jahr;
} DumptypeSnow;


#endif
