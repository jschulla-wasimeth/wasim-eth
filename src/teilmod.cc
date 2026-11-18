
/*   $Revision: 18 $ 
	$Date: 8.06.10 0:05 $  */

/*
(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <sstream>
#include <csignal>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "teilmod.h"   /* Definition der Teilmodelle */
#include "wasim.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

/* ----------------------- virtuelle Teilmodellklasse --------------------- */

/*
Diese Klasse enthaelt einige Routinen, welche allen Teilmodellen gemein ist
unter anderem die Zuweisung des Namens der Instanz, die Uebergabe des Namens der
Steuerdatei sowie diverse Nachrichten ueber den Fortschritt der Abarbeitung
des Programmes

wird waehrend der Modellentwicklung erkannt, dass bestimmte Routinen von vielen
oder allen Teilmodellen genutzt werden, so koennen diese Routinen in diese Klasse
aufgenommen werden

*/

teilmodellklasse::teilmodellklasse(wasim::mystream *Stream,double breitearg,double laengearg,double bezugsmeridianarg,const char *bezeichnung,
	char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo, int *p_zarray,  statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, 
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll ,  double p_yll, double p_cellsize  ,  double p_nodata,
	long mine, long ste , long tge , long mte , long jre,double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *pv_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp,bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU,int psubstepcnt,bool SectionNotFoundRethrow):
	modelobject(Stream, steuername, LogLevel, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, pv_mpiregioncnt, p_sendHalo, p_recvHalo, p_doSpinUp, p_vHRU),
tracer_count(0),tracer_info(0),tracer_module(0),
do_it(0),gridslesen(0),methode(0),fpeflag(0),InV(0),InVtag(0),stack(0),stackr(0),
stacktiefe(0),sunazimuthangle(0),sunaltitudeangle(0),sunzenithangle(0),
cosincidence(0),breite(0),laenge(0),bezugsmeridian(0),Anteil_hell(0),Anteil_dunkel(0),
dark_a(0),dark_e(0),nostrl(0),sunr(0),suns(0),pi(0),rad(0),zn(0),timekorr(0),substepcnt(psubstepcnt),routingmode(p_routingmode)
{
	/* Initialisieren der Adreseen und Daten mit den aus der
allgemeinen Initialisierungsroutine uebergebenen Daten */
	if (p_zarray==NULL)
	{
		// this is the case, if unsatzonmodel or soilmodel where not run and runoff-ignore>=1
		do_it=0;
		return;
	}
	nachricht  = progaus;
	min_to_int_end = p_min_to_int_end;
	numsubb    = p_numsubb;
	listsubb   = p_listsubb;
	if (listsubb) listsubb[0]=numsubb; /* number of subbasins to be considered for statistic files is stored in the zero-index */
	gridinfo   = p_gridinfo;
	zn         = gridinfo[1].pointer;
	zarray     = p_zarray;
	zone       = p_zone;
	FirstLastCol = p_FirstLastCol;
	zonmax     = p_zonmax;
	zonmaxwert = p_zonmaxwert;
	ncols      = p_ncols;
	nrows      = p_nrows;
	mpistartrow = p_mpistartrow;
	mpiendrow = p_mpiendrow;
	xll        = p_xll;
	yll        = p_yll;
	cellsize   = p_cellsize;
	nodata     = p_nodata;
	minutee    = mine;
	stundee    = ste;
	tage       = tge;
	monate     = mte;
	jahre      = jre;
	endtime = minf * minutee + jf * jahre + mf * monate + df * tage + hf * stundee;
	timekorr   = p_timekorr;
	LUtable    = p_LUtable;
	SOILtable  = p_SOILtable;
	MultiLayerLUTable = p_MultiLayerLUTable;

	if (ist_vorhanden(_RelCellArea_))
	RelCellArea=gridinfo[getindex(_RelCellArea_)].pointer;
	else
	RelCellArea=NULL;


	pi = 3.141592653;
	rad = 180/pi;

	/* Initialisieren des Arrays fuer die Abweichung der wahren Sonnenzeit von der Ortszeit
	in Abhaengigkeit vom Monat */
	timeequation[0]  =  11.3;
	timeequation[1]  =  -3.2;
	timeequation[2]  = -13.6;
	timeequation[3]  = -12.6;
	timeequation[4]  =  -4.2;
	timeequation[5]  =   2.8;
	timeequation[6]  =   2.5;
	timeequation[7]  =  -3.5;
	timeequation[8]  =  -6.3;
	timeequation[9]  =  -0.3;
	timeequation[10] =  10.0;
	timeequation[11] =  16.4;
	timeequation[12] =  11.3;
	timeequation[13] =  -3.2;

	/* geographische Koordinaten */
	breite         = breitearg;
	laenge         = laengearg;
	bezugsmeridian = bezugsmeridianarg;

	strcpy(instanzname,bezeichnung);
	(*WSstream) << "DEBUG: initializing module "<<instanzname<<"[end-entry]";

	if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}

	strcpy(tstr,"[");
	strcat(tstr,instanzname);
	strcat(tstr,"]");

	/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
	try
	{
		findeabschnitt(tstr);
	}
	catch (int ex)
	{
		// only if SectionNotFoundRethrow is false, the model shouild go on (which is only for optional sections and submodels
		if (SectionNotFoundRethrow)
		{
			throw ex;
		}
		else
		{
			do_it=0;
			return;
		}
	}

	fscanfln(steuerdatei,"%d",&do_it,"flow control parameter for ",tstr);
	fscanfln(steuerdatei,"%lf",&InV,"time interval for ",tstr);
	timestep = InV;

	hour_offset = min_to_int_end/60.0;

	if ((hour_offset>InV/60.0+0.00001)&&(strcmp(bezeichnung,"routing_model")!=0))
	{
		(*WSstream) << "FATAL: \n"
		"-------- error --------\n"
		"\n"
		"the time step "<<InV<<" minutes as given in the section for the module "
		<<bezeichnung<<" \n"
		"is shorter than the parameter min_to_int_end with "<<min_to_int_end<<" minutes\n"
		"from the [program_output]-section which gives the shift from \n"
		"the hour entry of meteo-data to the time step end\n"
		"this parameter (the 4-th in the section [program_output])\n"
		"should be set to 0 if the hour markes the end of the time step\n"
		"e.g. like in \"1984 01 01 24\" for daily time steps -> parameter set to 0\n"
		"and it should set to the variable $time, if the hour marks the time step start\n"
		"e.g. like in \"1984 01 01 00\" for daily time steps -> parameter set to 1440 or $time\n"
		"\n-------- error --------"<<"[end-entry]";
		throw(63);
	}

	/* Umrechnung der Intervallaenge in die Anzahl der Intervalle pro Tag */
	/* Diese Angabe wird oft benoetigt, z.B. bei der Schneeschmelze       */
	InVtag = 1440.0 / InV;

	if (do_it == 0)
	{
		/* Modell soll ignoriert werden, als gleich wieder verlassen, noch bevor Daten gelesen werden */
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
	}
}



void teilmodellklasse::run(long me, long se,long tg,long mt,long jr,char *name,char* additionalInfo)
{
	if (do_it == 0) return;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	char message[200];
	char additionalMessage[200];
	if (name)
	{
		sprintf(message, "%s", name);
	}
	else
	{
		sprintf(message, "%s", instanzname);
	}
	if (additionalInfo)
	{
		sprintf(additionalMessage, "%s", additionalInfo);
	}
	else
	{
		sprintf(additionalMessage, "%s", "");
	}
#ifdef _OPENMP
	//#pragma omp critical
	if (!omp_in_parallel())
	{
#endif
		if (nachricht == 1)
		{
			if ((InV > 0) && (InV < 60))
			{
				(*WSstream) << wasim::setw(20) << message << " " << tag << "." << monat << "." << jahr << " " << stunde << ":"<<minute<<" h " << _Meteodatenzeit_ << additionalMessage << "[end-entry]";
			}
			else // traditional message without minutes
			{
				(*WSstream) << wasim::setw(20) << message << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << _Meteodatenzeit_ << additionalMessage << "[end-entry]";
			}
		}
		else
		{
			if ((InV > 0) && (InV < 60))
			{
				if (((60*stunde + (int)((minute+60.0*hour_offset))) % (60*nachricht) == 0) && (stunde != 0))
				{
					(*WSstream) << wasim::setw(20) << message << " " << tag << "." << monat << "." << jahr << " " << stunde << ":"<<minute<<" h " << _Meteodatenzeit_ << additionalMessage << "[end-entry]";
				}
			}
			else
			{	// traditional message without minutes
				if ((stunde + (int)hour_offset) % nachricht == 0)
				{
					(*WSstream) << wasim::setw(20) << message << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << _Meteodatenzeit_ << additionalMessage << "[end-entry]";
				}
			}
		}
#ifdef _OPENMP
	}
#endif
}

int teilmodellklasse::doit()
{
	return do_it; /* gibt den wert des Ignoranzparameters do_it zurueck,
				do_it = 0 bedeutet, Teilmodell nicht ausfuehren */
}

double teilmodellklasse::deklination(long tag,long monat)
{
	long tagnr;
	int i;
	tagnr = 0;
	for (i=1;i<(int)monat;i++)
	{
		tagnr += mon[i];
	}
	tagnr += tag;
	return -23.4 * cos( (tagnr+10)/58.09155423); /* 58.09 = 2*pi/365  */
}

double teilmodellklasse::gettimeequation(long tag,long monat)
/* gibt die Anzahl Stunden an, die zur aktuellen mittleren Sonnenzeit
hinzuadduert werden muessen, um auf die reale SOnnenzeit zu kommen */
{
	return (double)(mon[(int)monat]-tag+1)/(double)mon[(int)monat]*(timeequation[(int)monat]/60.0) +
	(double)(tag-1)/(double)mon[(int)monat]*(timeequation[(int)monat+1]/60.0);
	/* einfache lineare Interpolation */
}



double teilmodellklasse::sonnenhoehe(long tag,long monat,double breite,double laenge,double stunde)
{ /* Angabe der exacten Stunde mit Bruchteil noetig, also z.B. 11.25 fuer 11:15 gesetzliche Zeit */
	double dekli,timeeq,
	stuwi,cosz,sonnenstand;
	/* int    i;
long tagnr;

tagnr = 0;
for (i=1;i<=monat;i++)
{
	tagnr += mon[i];
}
tagnr += tag; */
	timeeq = gettimeequation(tag,monat);
	/* Deklination der Sonne gegen die Aequatorialebene haengt nur vom Tag ab (1.Naeherung) */
	dekli  = deklination(tag,monat);

	stunde = stunde-(bezugsmeridian-laenge)/15.0+timeeq;
	stuwi  = 15.0*(12.0-stunde)/180.0*pi;
	cosz   = sin(breite*pi/180.0)*sin(dekli*pi/180.0)+cos(breite*pi/180.0)*cos(dekli*pi/180.0)*cos(stuwi);
	sonnenstand = acos(cosz)*180.0/pi;
	if (cosz<0.0) sonnenstand = sonnenstand-90.0; else sonnenstand = 90.0-sonnenstand;
	return sonnenstand;
}



void teilmodellklasse::getsuncoordinates(long tag,long monat,double breite,double laenge,double stunde,
double *sunazimuth,double *sunaltitude)
{   /* liefert bei Eingabe der geographischen Koordinaten in Grad und des Tages, Monats und der Stunde (in exacten Bruchteilen!)
	die Koordinaten der Sonne azimuth und altitude zurueck (in RADIANT!) */
	double dekli,timeeq,stuwi,cosz,cosOM;
	/* int    i;
long tagnr;

tagnr = 0;
for (i=1;i<monat;i++)
{
	tagnr += mon[i];
}
tagnr += tag; */
	breite = fmod(breite,360.0);
	laenge = fmod(laenge,360.0);
	timeeq = gettimeequation(tag,monat);
	dekli  = deklination(tag,monat);
	stunde = fmod(stunde-(bezugsmeridian-laenge)/15.0+timeeq+24.0,24.0);
	/*  stunde = stunde-(bezugsmeridian-laenge)/15.0+timeeq; */
	stuwi  = (12.0-stunde)/12.0 * pi;
	breite *= pi/180.0;
	dekli  *= pi/180.0;
	cosz   = sin(breite)*sin(dekli)+cos(breite)*cos(dekli)*cos(stuwi);
	cosOM  = (sin(dekli) * cos(breite) - cos(dekli)*sin(breite)*cos(stuwi))/sqrt(1-cosz*cosz);
	*sunaltitude = acos(cosz); /* sunaltitude ist eigentlich der zenithwinkel! */
	if (*sunaltitude <= pi/2.0) *sunaltitude = pi/2.0 - *sunaltitude; /* Sonnenstand ueber Horizont */
	else *sunaltitude = -(*sunaltitude - pi/2.0); /* Sonnenstand unter Horizont */
	if ((cosOM>-1.0)&&(cosOM<1.0)) *sunazimuth = acos(cosOM); else *sunazimuth=Pi;
	if (stunde > 12.0) *sunazimuth = 2*Pi - *sunazimuth;
}

double teilmodellklasse::getcosincidence(double aspect,double slope,double azimuth,double zenith)
{
	/* berechnet den cosinus des WInkels zwischen Sonnenvektor und Ebenennormale
	alle Angaben in Radiant, auch die Eingabedaten! */
	return cos(slope)*cos(zenith) + sin(slope)*sin(zenith)*cos(azimuth - aspect);
}

void teilmodellklasse::sunrise(long tag,long monat,double breite,double laenge,double *sr,double *ss)
{
	double z,dekli,timeeq;
	/* long   tagnr;
int    i;
*/
	/*  tagnr = 0;
for (i=1;i<monat;i++)
{
	tagnr += mon[i];
}
tagnr += tag; */
	timeeq = gettimeequation(tag,monat);
	dekli  = deklination(tag,monat);
	breite = fmod(breite,360.0);
	laenge = fmod(laenge,360.0);

	z = sin(breite*pi/180.0)/cos(breite*pi/180.0)*sin(dekli*pi/180.0)/cos(dekli*pi/180.0);
	if (fabs(z)<=1.0)
	{ /* da es ansonsten keinen Sonnenuntergang gibt! */

		if (z!=0.0) *sr = 12.0/pi * acos(z); else *sr = 6;
		*ss = 24.0 - *sr;
		*sr = *sr + (bezugsmeridian-laenge)/15.0-timeeq;
		*ss = *ss + (bezugsmeridian-laenge)/15.0-timeeq;
	}
	else
	{
		*ss = 0;
		*sr = 24;
	}
}

int teilmodellklasse::getTimeStep(std::string section)
{ 
	int locInV;
	if ((section.length() == 0)||(section.compare(this->instanzname)==0))
	{
		locInV=this->InV;
	}
	else
	{
		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}
		strcpy(tstr, (std::string("[") + section + std::string("]")).c_str());
		

		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		try
		{
			findeabschnitt(tstr);
		}
		catch (int ex)
		{
			throw ex;
		}
		int do_it = 0;
		fscanfln(steuerdatei, "%d", &do_it, "flow control parameter for ", tstr);
		fscanfln(steuerdatei, "%d", &locInV, "time interval for ", tstr);
		
		fclose(steuerdatei);
	}
	return locInV;
}


void teilmodellklasse::get_Tageslicht_Daten()
{
	double inta,inte;
	/* Berechnen des Sonnenauf- und Unterganges, um herauszufinden, ob es in diesem Zeitintervall
	ueberhaupt einen Sinn macht, die Strahlungsbilanz zu rechnen (da sie bei Nacht immer negativ
	ist und somit die Verdunstung nur vom Wind abhaengt ausser vom Saettigungsdefizit) */
	sunrise(tag,monat,breite,laenge,&sunr,&suns);
	inta=InV<=1440?(stunde+hour_offset)-InV/60.0:(stunde+hour_offset)-24; /* if timestep is greater than 24hours, the timestep is assumed to be a multiple of 1440 min (24h, 1day), so the end hour is 24 in each interval */
	inte=(stunde+hour_offset);
	sunr = fmod(sunr+24.0,24.0); /* Sonnenaufgang zwischen 0 ... 24 Uhr fixieren */
	suns = fmod(suns+24.0,24.0); /* Sonnenuntergang ebenfalls fixieren */
	/* wie lange ist es am Anfang eines Zeitintervalls noch dunkel ? */

	/* dark_a */
	if (sunr<inta)
	dark_a = 0.0;
	else
	if (sunr<inte)
	dark_a = sunr - inta;
	else /* Stunde komplett vor Sonnenaufgang */
	dark_a = 0.0;

	/* dark_e */
	if (suns>sunr)
	{ /* der Normalfall */
		if ((suns>inte)&&(sunr<inte))
		dark_e = 0.0;
		else
		if ((suns>inta)&&(suns<inte))
		dark_e = inte - suns;
		else
		dark_e = InV/60.0;
	}
	else
	{ /* der Fall, dass die Sonne erst nach Mitternacht gesetzl. Zeit untergeht */
		if (suns>inte)
		dark_e = 0.0;
		else
		if (suns>inta)
		dark_e = inte - suns;
		else
		dark_e = sunr > inte ? InV/60.0 : 0.0;
	}

	Anteil_dunkel = (dark_a + dark_e)/((InV>1440?1440:InV)/60.0);
	/* Wenn die Sonne im Intervall erst unter, dann aufgeht, muss zum Anfang und Ende
	des Intervalls das Komplement zu 1 von dark_a und dark_e gebildet und addiert werden
	(1-dark_a)+(1-dark_e) = 2-(dark_a+dark_e) = 2.0 - Anteil_dunkel */
	if (Anteil_dunkel>1.0) Anteil_dunkel--;
	/* wenn die Sonne nur ganz kurz auftaucht, kann sie auch gleich
	untergetaucht bleiben */
	if (Anteil_dunkel > 0.95) Anteil_dunkel=1.0;
	Anteil_hell   = 1.0 - Anteil_dunkel;

	/* wenn das aktuelle Berechnungsintervall komplett zu mehr als 95% vor Sonnenaufgang
	oder nach Sonnenuntergang liegt, wird keine Strahlung beruecksichtigt */
	if (Anteil_hell<0.05)
	nostrl = 1;
	else
	nostrl = 0;
}

void teilmodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count,tracer_info_type *p_tracer_info,TLastResults *LastResults,int maxLastResults) /* initializing substance_transport algorithms (getting tracer_info table adress) */
{
	try
	{
		retvalue=0;
		int i,j,l,k;
		char *identify;
		identify=newstring(300,1);
		tracer_count = p_tracer_count;
		if (tracer_count>0)
		{
			tracer_info = p_tracer_info;
			tracer_module = (tracer_module_type*)calloc(tracer_count+1,sizeof(tracer_module_type));
			for (i=1;i<=tracer_count;i++)
			{
				if (tracer_info[i].radioactive==1)
				/* this is the factor, with which the actual concentration of a radioactive
			tracer has to be multiplied to get the radioactive decay for one time step */
				tracer_info[i].decfact = exp(-(InV/60/24)/(tracer_info[i].half_time*_1_ln2_));
				/* initializing a statistc register set for this tracer */
				tracer_module[i].zone = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
				if (!tracer_module[i].zone) speichererror(1);
				for (j=0;j<=zonmax;j++)
				{ tracer_module[i].zone[j].n         = tracer_module[i].zone[j].sx        = 0;
					tracer_module[i].zone[j].zonenwert = zone[j].zonenwert;
					tracer_module[i].zone[j].rel       = zone[j].rel;
					tracer_module[i].zone[j].RelCellSum = 0;
				}
				/* initializing the output stream for the tracer statistics */
				tracer_module[i].statfilename=newstring(_stringsize_+1,1);
				sprintf(tracer_module[i].statfilename,"%s%.3s%.5s.%.3s",tracer_info[i].path,tracer_info[i].name,instanzname,tracer_info[i].statext);
				tracer_module[i].stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&tracer_module[i].zone,(int)zonmax,(int)tracer_info[i].statcode,(int)99,(char*)tracer_module[i].statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
				if (retvalue!=0) throw retvalue;
				sprintf(identify,"tracer concentration for tracer %d (%s) ",i,tracer_info[i].name);
				tracer_module[i].stat->schreibzonenwerte(identify,instanzname,LastResults,maxLastResults);
				/* allocating memory for the concentration in a storage */
				j=1;
				while (gridinfo[j].identification!=NULL)
				{
					j++;
					if (j==gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no free grid left. Please contact the author to set the number of available grids to a larger value "<<gridmaxanzahl<<"[end-entry]";
						throw(12);
					}
				}
				sprintf(identify,"concstor_tracer_%d_%s",i,instanzname);
				tracer_module[i].gridname=newstring(_stringsize_+1,1);
				sprintf(tracer_module[i].gridname,"%s%.3s%.5s.%.3s",tracer_info[i].path,tracer_info[i].name,instanzname,tracer_info[i].gridext);
				if (gridslesen!=1)
				/* Grid wird nicht von Platte eingelesen */
				{
					tracer_module[i].concstor=matrix(1,nrows,1,ncols);
					initialize_grid(j, tracer_module[i].gridname, tracer_module[i].concstor, NULL, 0, identify, tracer_info[i].writegrid, 1, nrows);
					for (k = mpistartrow; k <= mpiendrow; k++)
					{
						for (l = 1; l <= ncols; l++)
						{
							tracer_module[i].concstor[k][l] = zn[k][l] == nodata ? nodata : 0;
						}
					}
					tracer_module[i].storindex = j;
				}
				else
				{
					/* alte Dateien einlesen */
					readgridfromfile(tracer_module[i].gridname, &(tracer_module[i].concstor), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_);
					initialize_grid(j, tracer_module[i].gridname, tracer_module[i].concstor, NULL, 0, identify, tracer_info[i].writegrid, 1, nrows);
					tracer_module[i].storindex = j;
				}
				/* allocating a grid for the concentration in the outflow */
				j++;
				sprintf(identify,"concflux_tracer_%d_%s",i,instanzname);
				tracer_module[i].concflux=matrix(1,nrows,1,ncols);
				initialize_grid(j, tracer_module[i].gridname, tracer_module[i].concflux, NULL, 0, identify, 0, 1, nrows);
				tracer_module[i].fluxindex = j;
			}
		}
		if (identify) free(identify); identify=NULL;
	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	{
		retvalue=103;
	}
}

void teilmodellklasse::substance_transport_delete() /* deleting all the grids for substance concentrations */
{
	int i;
	if (tracer_count>0)
	{
		for (i=1;i<=tracer_count;i++)
		{
			free_matrix(tracer_module[i].concstor,1,nrows,1,ncols);
			delete_grid(tracer_module[i].storindex);
			free_matrix(tracer_module[i].concflux,1,nrows,1,ncols);
			delete_grid(tracer_module[i].fluxindex);
			if (tracer_module[i].gridname) free(tracer_module[i].gridname);
			tracer_module[i].gridname=NULL;
			tracer_module[i].stat->ClearMemory();
			if (tracer_module[i].stat) delete(tracer_module[i].stat);
			if (tracer_module[i].statfilename) free(tracer_module[i].statfilename);
			tracer_module[i].statfilename=NULL;
			if (tracer_module[i].zone) 
			free(tracer_module[i].zone);
			tracer_module[i].zone=NULL;
		}
	}
}

double teilmodellklasse::decay(double conc, double decfact)
/* for calculating Ci = Ci-1 * decfact; decfact = exp(-(InV/60/24)/half_time); returns the new concentration after dte decay */
{
	return conc*decfact;
}

double teilmodellklasse::mixing(double c1,double c2, double v1,double v2)
/* calc. mixing concentrations for volume1 (c1) with volume 2 (c2); returns mixing concentration */
{
	if (v1+v2!=(double)0.0)
	return (v1*c1+v2*c2)/(v1+v2); /* return the mixing concentration */
	else
	return (double)0.0;
}


double teilmodellklasse::evaporating_changes(double c,double v1,double ve) /* calc. and returning the concentration c after evaporating the amount ve from v */
{
	if (v1-ve!=(double)0.0)
	return (v1*c)/(v1-ve);  /* agglomerating e.g. salts due to evaporation / transpiration */
	else
	return (double)0.0;
}

void teilmodellklasse::fpe_check(int i,int j,int k,int z,double a1,double a2,double a3,double a4,double a5)
{
	if (fpeflag==1)
	{
		(*WSstream) << "INFO: fpe in i1: "<<i<<", i2: "<<j<<", i3: "<<k<<", i4: "<<z
		<<"\na1: "<<a1<<" a2: "<<a2<<" a3: "<<a3<<" a4: "<<a4<<" a5: "<<a5<<" "<<"[end-entry]";
		fpeflag=2;
	}
	else fpeflag=0;
}

teilmodellklasse::~teilmodellklasse()
{
	//  (*WSstream) << "INFO: deleting module "<<instanzname<<"[end-entry]";
}

void teilmodellklasse::zonenidentification(double ***zn,statrecord **zone,Tzarray **zarray,int &zonmax,int &zonmaxwert)
{
	int  i,j;
	long zonenid;
	statrecord tempzone;
	unsigned short int  weiter;

	// identifies the zone codes and rankes them from lowest to highest, generating an index without nodata-value in the index-table 

	// the foll. variables are allocated and initializied outside of this module 
	*zn         = NULL;
	*zarray     = NULL;
	*zone       = NULL;
	zonmaxwert = zonmax=0;

	// searching the zonegrid (zn)
	*zn     = gridinfo[getindex(_Zonengrid_Bodenmodell_)].pointer;

	if (mpi_rank == 0)
	{
		(*WSstream) << "ranking zones from " << gridinfo[getindex(_Zonengrid_Bodenmodell_)].name << " for " << instanzname << "[end-entry]";
		zzone = (statrecord *)calloc(maxzonenzahl + 1, sizeof(statrecord)); /* allocating maximum allowed memeory for a temporary statistics-array */
		if (!zzone) speichererror(1);
		zzarray = (Tzarray *)calloc((size_t)maxzonenwert + 1, sizeof(Tzarray)); /* allocating maximum allowed memory for temporary zone-array */
		if (!zzarray) speichererror(1);
		zzone[0].n = zzone[0].RelCellSum = 0.0;
		zzone[0].sx = 0.0;
		zzone[0].rel = 1.0;
		zzone[0].zonenwert = 1;

		zonmax = 0;
		zonmaxwert = 0;

		for (i = 1; i < maxzonenwert; i++)   zzarray[i] = 0;
		for (i = 1; i <= nrows; i++)
		{
			for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((*zn)[i][j] != nodata)
				{
					zonenid = (unsigned long)floor((*zn)[i][j] + 0.5);
					if ((zonenid<0) || (zonenid>maxzonenwert))
					{
						(*WSstream) << "FATAL: wrong zone index in " << gridinfo[getindex(_Zonengrid_Bodenmodell_)].name
							<< "!\n Zones have to range from 0 to " << maxzonenwert << "!" << "[end-entry]";
						throw(3);
					}
					if (zzarray[zonenid] == 0)
					{ /*gibt es die Zone mit der ID zonenid noch nicht? */
						zonmax = zonmax + 1;           /* nein, also in die Liste einfuegen */
						if (zonmax > maxzonenzahl)
						{
							(*WSstream) << "FATAL: too much zones in " << gridinfo[getindex(_Zonengrid_Bodenmodell_)].name
								<< ", (max. " << maxzonenzahl << "!)" << "[end-entry]";
							throw(4);
						}
						if (zonenid > zonmaxwert) zonmaxwert = zonenid;
						zzarray[zonenid] = zonmax;
						zzone[zzarray[zonenid]].zonenwert = (*zn)[i][j];
						zzone[zzarray[zonenid]].n = 0;
						zzone[zzarray[zonenid]].RelCellSum = 0;
						zzone[zzarray[zonenid]].sx = 0;
					}
				}
			}
		}
		/* now ranking the zones for zone-value */
		for (i = 1; i <= zonmax - 1; i++)
		{
			weiter = true;
			for (j = 2; j <= zonmax; j++)
			{
				if ((zzone[j].zonenwert) < (zzone[j - 1].zonenwert))
				{
					tempzone = zzone[j - 1];
					zzone[j - 1] = zzone[j];
					zzone[j] = tempzone;
					zzarray[(int)floor(zzone[j].zonenwert + 0.5)] = j;
					zzarray[(int)floor(zzone[j - 1].zonenwert + 0.5)] = j - 1;

					weiter = false;
				}
			}
			if (weiter == true) break;
		}

		/* for using not to much memory, now a re-allocation of zone and zarray is done
	(with only as many records as needed (number of different zones)) */
		*zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
		if (!*zone) speichererror(1);
		for (i = 0; i <= zonmax; i++)
		{
			(*zone)[i].n = zzone[i].n;         /* number of values for zonal statistics */
			(*zone)[i].sx = zzone[i].sx;        /* sum of values for zonal statistcs     */
			(*zone)[i].zonenwert = zzone[i].zonenwert; /* what is the code of the zone?  */
			(*zone)[i].rel = zzone[i].rel;       /* relativ fraction of the zone from the entire grid */
			(*zone)[i].RelCellSum = zzone[i].RelCellSum;
		}
		if (zzone) free(zzone);

		*zarray = (Tzarray *)calloc(zonmaxwert + 1, sizeof(Tzarray));
		if (!*zarray) speichererror(1);
		for (i = 0; i <= zonmaxwert; i++) (*zarray)[i] = zzarray[i];
		if (zzarray) free(zzarray);

		/* calculating relative fraction of all zones from the entire grid */
		/* scanning the zone grid and adding to all registers a "1" per valid cell */
		for (i = 1; i <= nrows; i++)
		{
			for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((*zn)[i][j] != nodata)
				{
					(*zone)[(*zarray)[(int)((*zn)[i][j])]].sx++;
					if (RelCellArea)
					{
						(*zone)[(*zarray)[(int)((*zn)[i][j])]].n += RelCellArea[i][j];
						(*zone)[(*zarray)[(int)((*zn)[i][j])]].RelCellSum += RelCellArea[i][j];
					}
					else
					{
						(*zone)[(*zarray)[(int)((*zn)[i][j])]].n++;
						(*zone)[(*zarray)[(int)((*zn)[i][j])]].RelCellSum++;
					}
				}
			}
		}

		/* total number of valid values */
		unsigned int cellcount = 0;
		for (i = 1; i <= zonmax; i++)
		{
			cellcount += (*zone)[i].n;
		}

		if (cellcount == 0)
		{
			(*WSstream) << "FATAL: all cells in " << gridinfo[getindex(_Zonengrid_Bodenmodell_)].name
				<< " are NODATA!" << "[end-entry]";
			throw(7);
		}

		/* relativ fraction per zone */
		for (i = 1; i <= zonmax; i++)
		{
			(*zone)[i].rel = (*zone)[i].n / cellcount;
		}
		(*zone)[0].rel = cellcount*cellsize*cellsize / 1000000.0; /* in zone[0].rel: area of the entire basin */

		/* not more zones than "maxteilgebiete" are supported */

		if (zonmax > maxteilgebiete)
		{
			(*WSstream) << "FATAL: too much subbasins (" << zonmax << ") in " << instanzname
				<< " (max. " << maxteilgebiete << "). Ende!" << "[end-entry]";
			throw(4);
		}
	}

	// only if there are more than one instances, we need to send the zone structure to other ranks
#ifdef _MPI
	if (mpi_size > 1)
	{
		MPI_Bcast(&zonmax, 1, MPI_INT, 0,MPI_COMM_WORLD);
		MPI_Bcast(&zonmaxwert, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (mpi_rank > 0)
		{
			(*zone) = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			if (!zone) speichererror(1);
			(*zarray) = (Tzarray *)calloc(zonmaxwert + 1, sizeof(Tzarray));
			if (!zarray) speichererror(1);
		}
		MPI_Bcast(&((*zarray)[0]), zonmaxwert + 1, MPI_INT, 0, MPI_COMM_WORLD);
		int cnt = int(sizeof(statrecord) / sizeof(double))*(zonmax + 1);

		MPI_Datatype zonevectortype;
		MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
		MPI_Type_commit(&zonevectortype);

		for (int i = 0; i <= zonmax; i++)
		{
			MPI_Bcast(&((*zone)[i]), 1, zonevectortype, 0,MPI_COMM_WORLD);
		}
		MPI_Type_free(&zonevectortype);

	}
#endif

	/* ###################### ende of zone analysis #################### */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////  SpecialOutputClass with its helper classes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StringUtils::stripInt(std::string &row, std::string sWhatParam, int *param)
{
	int index1, index2;
	std::string buffer;

	buffer = row.substr(row.find(sWhatParam));
	index1 = buffer.find("=") + 1;
	index2 = buffer.find(";");
	buffer = buffer.substr(index1, index2 - index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> *param;
}

void StringUtils::stripDouble(std::string &row, std::string sWhatParam, double *param)
{
	int index1, index2;
	std::string buffer;

	buffer = row.substr(row.find(sWhatParam));
	index1 = buffer.find("=") + 1;
	index2 = buffer.find(";");
	buffer = buffer.substr(index1, index2 - index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> *param;
}

void StringUtils::stripString(std::string &row,std::string sWhatParam,std::string &param)
{
	int index1,index2;
	std::string buffer;

	buffer=row.substr(row.find(sWhatParam));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> param;
}

std::string StringUtils::stripSectionBetweenBrackets(std::string &row)
{
	int index1,index2;
	std::string buffer;

	buffer=row.substr(row.find("{"));
	index2 = buffer.find("}");
	buffer=buffer.substr(1,index2);
	row = row.substr(row.find("}")+1);
	return buffer;
}

int StringUtils::stripDoubleVector(std::string &row,std::string sWhatParam,std::vector<double> &X)
{
	int index1,index2;
	std::string buffer;
	std::vector<std::string> words;

	if (row.find(sWhatParam)==std::string::npos) return 0;
	buffer=row.substr(row.rfind(sWhatParam));

	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	
	buffer=buffer.substr(index1,index2-index1);
	words.clear();
	StringSplit(buffer,"\t ;,",words,false);
	for (int i=0;i<words.size();i++)
	{
		if (words[i].size()>0)
		{
			std::stringstream ss;
			ss << words[i];
			double value;
			ss >> value;
			X.push_back(value);
		}
	}
	return X.size();
}

int StringUtils::stripIntVector(std::string &row,std::string sWhatParam,std::vector<int> &X)
{
	int index1,index2;
	std::string buffer;
	std::vector<std::string> words;

	if (row.find(sWhatParam)==std::string::npos) return 0;

	buffer=row.substr(row.rfind(sWhatParam));

	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	
	buffer=buffer.substr(index1,index2-index1);
	words.clear();
	StringSplit(buffer,"\t ;,",words,false);
	for (int i=0;i<words.size();i++)
	{
		if (words[i].size()>0)
		{
			std::stringstream ss;
			ss << words[i];
			int value;
			ss >> value;
			X.push_back(value);
		}
	}
	return X.size();
}


//////////////////////////////////////////// SpecialOutputEntityClass

SpecialOutputEntityClass::SpecialOutputEntityClass(SpecialOutputFileClass *p_parent,wasim::mystream *Stream,std::string entityString,gridinfotype *gridinfo,
	double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow,
	int numsubsteps):
parent(p_parent),WSstream(Stream),gridinfo(gridinfo),xll(xll),yll(yll),cs(cs),ncols(ncols),nrows(nrows),
iNumSubSteps(numsubsteps),
bUsesSubSteps(false),
iSubStepsCumulationType(NONE)
{
	// read the entities ID, Smbol, X, Y and Level values
	su.stripString(entityString,"ID",ID);
	su.stripString(entityString,"Symbol",Symbol);
	su.stripDoubleVector(entityString,"Xcoords",Xcoord);	
	su.stripDoubleVector(entityString,"Ycoords",Ycoord);	
	su.stripIntVector(entityString,"Level",Level);	

	if (Xcoord.size()!=Ycoord.size())
	{
		(*WSstream) << "FATAL: Error in section [special_output]: count of X coordinates dos not match the number of Y coordinates for entity "<<ID<<"[end-entry]";
		throw(1006);
	}

	if (Level.size()!=Xcoord.size() && Level.size()!=0)
	{
		(*WSstream) << "FATAL: Error in section [special_output]: number of levels dos not match the number of X and Y coordinates for entity "<<ID<<". EIther each point has a level or none has."<<"[end-entry]";
		throw(1007);
	}

	// calculating the indices for the column and row of the gridinfo.pointer
	for (int i=0;i<Xcoord.size();i++)
	{
		col.push_back(ceil((Xcoord[i]-xll)/cs));
		if ((col.back()<1)||(col.back()>ncols))
		{
			(*WSstream) << "FATAL: Error in section [special_output]: x-coordinate "<<Xcoord[i]<<", defined as point "<<i<<" for entity "<<ID<<", is located outside the actual grid."<<"[end-entry]";
			throw(1008);
		}
	}
	for (int i=0;i<Ycoord.size();i++)
	{
		row.push_back(nrows-floor((Ycoord[i]-yll)/cs));
		if ((row.back()<1)||(row.back()>nrows))
		{
			(*WSstream) << "FATAL: Error in section [special_output]: y-coordinate "<<Ycoord[i]<<", defined as point "<<i<<" for entity "<<ID<<", is located outside the actual grid."<<"[end-entry]";
			throw(1008);
		}
	}
	this->numPoints=col.size();
	// suchen des Grids mit der identiy ID
	gridinfoindex = ist_vorhanden(ID.c_str());

	// looking for special entities which may require summing up or averaging the content of the grid when sub steps are used for unsatzin model, surface routing and griundwater model
	if (this->iNumSubSteps>1)
	{
		//unsatzon:
		if (
			(ID.find(_Basisabfluss_)!=std::string::npos) ||
			(ID.find(_irrigation_amounts_)!=std::string::npos) ||
			(ID.find(_macropore_infiltration_)!=std::string::npos) ||
			(ID.find(_Grundwasserneubildung_)!=std::string::npos) ||
			(ID.find(_Capillar_Uprise_)!=std::string::npos) ||
			(ID.find(_Percolation_)!=std::string::npos) ||
			(ID.find(_DIRECTRUNOFFSILTINGUP_)!=std::string::npos) ||
			(ID.find(_Direktabfluss_)!=std::string::npos) ||
			(ID.find(_Interflow_)!=std::string::npos) ||
			(ID.find(_Infiltrationsexzess_)!=std::string::npos) ||
			(ID.find(_Total_Infiltration_)!=std::string::npos) ||
			(ID.find(_Drainageabfluss_)!=std::string::npos) ||
			(ID.find(_river_infiltration_)!=std::string::npos) ||
			(ID.find(_river_exfiltration_)!=std::string::npos)
		   )
		{
			bUsesSubSteps=true;
			iSubStepsCumulationType=ADD;
		} else if (
		//unsatzon:
			(ID.find(_soil_bottom_fluxes_)!=std::string::npos) ||
			//surface routing:
			(ID.find(_SURFACEROUTING_FLUXES_)!=std::string::npos) ||
			(ID.find(_SURFACEROUTING_HEAD_)!=std::string::npos) ||
			(ID.find(_SURFACEROUTING_QIN_)!=std::string::npos) ||
			(ID.find(_SURFACEROUTING_QOUT_)!=std::string::npos) ||
			// gw-model:
			(ID.find(_gw_flow_x_)!=std::string::npos) ||
			(ID.find(_gw_flow_y_)!=std::string::npos) ||
			(ID.find(_gw_leakage_)!=std::string::npos)
			)
		{
			bUsesSubSteps=true;
			iSubStepsCumulationType=ADD;
		}
	}
}

SpecialOutputEntityClass::~SpecialOutputEntityClass()
{
	Xcoord.clear();
	Ycoord.clear();
	Level.clear();
	ID.clear();
	Symbol.clear();
	row.clear();
	col.clear();
}

int SpecialOutputEntityClass::addValues(buffertype *buffer,unsigned long bufcnt,int startindex,bool bInSubStepsLoop,int actSubStep)
{
	int actindex=startindex;
	if ( !bUsesSubSteps && bInSubStepsLoop)
	{
		// for all entities, which are either not affected by the submodels that are in the sub steps loop or that doesn't 
		// have to be handled specially when in this loop, the actual call must be ignored, since we are in the sub steps loop
		return actindex+col.size();
	}
	else
	{
		// if either in sub steps loop and the entity needs to be handled specially or if not in sub steps loop 
		// (all entities will be handled there, either by the normal procedure of adding a value or by calculating 
		// the average and setting the new index for specially handled values from the sub steps loop)
	
		for (int index=0;index<col.size();index++)
		{
			int ri=row[index];
			int ci=col[index];
			int mpir = (parent->parent->MPIREG[ri]) ? parent->parent->MPIREG[ri][ci]:-1;
			int mpir2=parent->parent->mpi_rank;
			if ((gridinfoindex>0)&&(mpir==mpir2))
	//		if ((gridinfoindex>0)&&(parent->parent->MPIREG[row[index]][col[index]]==parent->parent->mpi_rank))
			{
				if ((gridinfo[gridinfoindex].stackpointer!=NULL) && (Level.size()>0))
				{ // a stack
					if ((row[index] >= parent->parent->mpistartrow) && (row[index] <= parent->parent->mpiendrow) && (col[index]>0) && (col[index] <= ncols) && (Level[index] <= gridinfo[gridinfoindex].levels))
					{
						if (bUsesSubSteps && bInSubStepsLoop)
						{
							// add value to the value from the last sub step
							if (actSubStep>1)
							{	// add to value from last sub step
								buffer[bufcnt].wert[actindex] += gridinfo[gridinfoindex].stackpointer[row[index]][col[index]][Level[index]];
							}
							else
							{  // initialize for a new time step with the actual value
								buffer[bufcnt].wert[actindex] = gridinfo[gridinfoindex].stackpointer[row[index]][col[index]][Level[index]];
							}
							if ((iNumSubSteps==actSubStep) && (iSubStepsCumulationType==AVG))
							{
								buffer[bufcnt].wert[actindex] /= iNumSubSteps;
							}
						}
						else if (!bUsesSubSteps)
						{
							// add new value only once for entities that are not touched in sub steps loop or that doesn't need to be handled specially
							buffer[bufcnt].wert[actindex] = gridinfo[gridinfoindex].stackpointer[row[index]][col[index]][Level[index]];
						}
					}
					else
					{
						buffer[bufcnt].wert[actindex] = -9999;
					}
				}
				else // only a 2D-grid
				{
					if ((row[index] >= parent->parent->mpistartrow) && (row[index] <= parent->parent->mpiendrow) && (col[index]>0) && (col[index] <= ncols))
					{
						if (bUsesSubSteps && bInSubStepsLoop)
						{
							// add value to the value from the last sub step
							if (actSubStep>1)
							{	// add to value from last sub step
								buffer[bufcnt].wert[actindex] += gridinfo[gridinfoindex].pointer[row[index]][col[index]];
							}
							else
							{  // initialize for a new time step with the actual value
								buffer[bufcnt].wert[actindex] = gridinfo[gridinfoindex].pointer[row[index]][col[index]];
							}
							if ((iNumSubSteps==actSubStep) && (iSubStepsCumulationType==AVG))
							{
								buffer[bufcnt].wert[actindex] /= iNumSubSteps;
							}
						}
						else if (!bUsesSubSteps)
						{
							// add new value only once for entities that are not touched in sub steps loop or that doesn't need to be handled specially
							buffer[bufcnt].wert[actindex] = gridinfo[gridinfoindex].pointer[row[index]][col[index]];
						}
					}
					else
					{
						buffer[bufcnt].wert[actindex] = -9999;
					}
				}
			}
			else
			{
				buffer[bufcnt].wert[actindex] = -9999;
			}
			actindex++;
		}
	}
	return actindex;
}

int SpecialOutputEntityClass::ist_vorhanden(const char *name)
{
	int i;
	i=1;
	while (strcmp(gridinfo[i].identification,name)!=0)
	{
		i++;
		if ((i==gridmaxanzahl)||(gridinfo[i].pointer==NULL)) return 0;
	}
	return i;
}


//////////////////////////////////////////// SpecialOutputFileClass

SpecialOutputFileClass::SpecialOutputFileClass(SpecialOutputClass *p_parent,wasim::mystream *Stream,std::string record,
	gridinfotype *p_gridinfo, double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow,
	char *DefaultOutputDirectory,int substepcnt):
parent(p_parent),WSstream(Stream),gridinfo(p_gridinfo),
	numEntities(0),buffer(NULL),bufmax(0),bufcnt(0),
	outputfile(NULL),iNumSubSteps(substepcnt),incrementBufCntNextTime(true)
{
//	gridinfo = p_gridinfo;
//	WSstream = Stream;
//	parent = p_parent;

	// for each entity in this file: creating an instance of the class SpecialOutputEntityClass
	su.stripString(record,"header",header);
	su.stripString(record,"filename",filename);
	char *tempname;
	tempname= new char[_stringsize_+1];
	sprintf(tempname,"%s",filename.c_str());
	ExtendFileName(DefaultOutputDirectory,tempname);
	filename=std::string(tempname);
	delete[] tempname;
	std::string entityString=record;
	// extract entityString
	while (entityString.find("entity")!=std::string::npos)
	{
		record = su.stripSectionBetweenBrackets(entityString);
		entities.push_back(new SpecialOutputEntityClass(this,WSstream,record,gridinfo,xll,yll,cs,ncols,nrows,mpistartrow,mpiendrow,this->iNumSubSteps));
	}
	// how many data colmuns are there for all entities of this file?

	numEntities=0;
	for (std::vector<SpecialOutputEntityClass*>::iterator iter=entities.begin();iter<entities.end();iter++)
	{
		numEntities+=(*iter)->getNumPoints();
	}

	// create the buffer structure for the SpecialOutputFileClass
	bufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(numEntities+1)*sizeof(double)));
	if (bufmax<1) bufmax=1;
	bufcnt      = 0; /* initialising the buffer index counter with 0 */
	buffer      = (buffertype*)calloc(bufmax+1,sizeof(buffertype));
	for (unsigned long bufrow=1;bufrow<=bufmax;bufrow++) buffer[bufrow].wert=(double*)calloc(numEntities+1,sizeof(double));
	
	if (parent->mpi_rank==0)
	{
		// create the outputfile and write the header
		if ((outputfile=fopen(filename.c_str(),"wt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening file for special output of file\n"<<filename<<"\n-> check file name, disk space or user rights"<<"[end-entry]";
			throw(1007);
		}
		// 1st row
		fprintf(outputfile,"YY\tMM\tDD\tHH");
		for (std::vector<SpecialOutputEntityClass*>::iterator iter=entities.begin();iter<entities.end();iter++)
		{
			int num=(*iter)->getNumPoints();
			for (int c=0;c<num;c++)
			{
				fprintf(outputfile,"\t%s/%d/%d",(*iter)->getSymbol().c_str(),(*iter)->getCol()[c],(*iter)->getRow()[c]);
				if ((*iter)->getLevel().size()>0)
				{
					fprintf(outputfile,"/%d",(*iter)->getLevel()[c]);
				}
			}
		}
		// 2nd row
		fprintf(outputfile,"\n%s\nYY\tMM\tDD\tHH",header.c_str());
		// 3rd row: contains coordinates instead of columns/rows
		for (std::vector<SpecialOutputEntityClass*>::iterator iter=entities.begin();iter<entities.end();iter++)
		{
			int num=(*iter)->getNumPoints();
			for (int c=0;c<num;c++)
			{
				fprintf(outputfile,"\t%s/%lg/%lg",(*iter)->getSymbol().c_str(),(*iter)->getXcoord()[c],(*iter)->getYcoord()[c]);
				if ((*iter)->getLevel().size()>0)
				{
					fprintf(outputfile,"/%d",(*iter)->getLevel()[c]);
				}
			}
		}

		fprintf(outputfile,"\n");
		fclose(outputfile);
	}
}
	
void SpecialOutputFileClass::ExtendFileName(char *directory,char *dateiname)
{
	bool bIsComplete=false;
	std::string temp(dateiname);
	// check if the dateiname contains a drive letter (Windows) or starts with a / or ~ (Linux)
#if defined(WIN32)
	if (temp.find(":\\")!=std::string::npos) bIsComplete=true;
#elif defined(WIN64)
	if (temp.find(":\\")!=std::string::npos) bIsComplete=true;
#else
	if (temp.find_first_of("/~")==0) bIsComplete=true; // starts with a / or a ~ meaning there are absolute paths defined
#endif
	
	if ((!bIsComplete)&&(directory!=NULL))
	{
		char *completename=new char[_stringsize_+1];
		sprintf(completename,"%s%s",directory,dateiname);
		sprintf(dateiname,"%s",completename);
		delete[] completename;
	}
}


SpecialOutputFileClass::~SpecialOutputFileClass()
{
	this->flushbuffer();
	// deleting all entities data
	for (std::vector<SpecialOutputEntityClass*>::iterator iter = entities.begin(); iter<entities.end(); iter++)
	{
		delete *iter;
	}
	entities.clear();
	// deleting buffer data
	if (buffer) 
	{
		for (int i=1;i<=bufmax;i++) 
		{
			if (buffer[i].wert) 
			{
				free(buffer[i].wert);
				buffer[i].wert=NULL;
			}
		}
		free(buffer);
		buffer=NULL;
	}

}


void SpecialOutputFileClass::flushbuffer()
{
	// implement here the flushing of the buffers structure into a file and resetting the counters
	int i;
	unsigned long j;

	if (bufcnt>0)
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		if (parent->mpi_rank==0)
		{
			(*WSstream) << "write "<<filename.c_str()<<"[end-entry]";
			if ((outputfile=fopen(filename.c_str(),"a+"))==NULL)
			{
				(*WSstream) << "FATAL: \nerror opening file "<<filename.c_str()<<"[end-entry]";
				throw(65);
			}
		}

		// reducing the buffers from all instances into a large buffer-array in rank 0
#ifdef _MPI
		// Variable Array-size 
		int rows = bufcnt;
		int cols = numEntities;
		int total_elements = rows * cols;

		// Lokale 1D-Array 
		std::vector<double> local_array(total_elements,-9999); // all elements are initialized with -9999
		// now copy the values into the 
		long position = 0;
		for (j = 1; j <= bufcnt; j++)
		{
			std::memmove(local_array.data()+(j-1)*numEntities,
						&(buffer[j].wert[1]),
						sizeof(double)*numEntities);
		}

		// allocating the result vector on root only
		std::vector<double> global_result;
		if (parent->mpi_rank == 0) {
			global_result.resize(total_elements);
		}

		// calling the MPI-reduction
		MPI_Reduce(local_array.data(), global_result.data(), total_elements, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#endif


		for (j=1;j<=bufcnt;j++)
		{
			if (parent->mpi_rank==0)
			{
				fprintf(outputfile,"%ld\t%ld\t%ld\t%ld",buffer[j].j,buffer[j].m,buffer[j].t,buffer[j].s);
			}

			// Gather the maximum of all instances (if a cell is not in a instance, the value is -9999, so we need to get the maximum simply
			double *dmax = new double[numEntities+1];
#ifdef _MPI
//			MPI_Reduce(&(buffer[j].wert[1]), &(dmax[1]), numEntities, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);
			if (parent->mpi_rank == 0)
			{
				for (int i = 1; i <= numEntities; i++)
				{
					dmax[i] = global_result[(j - 1) * numEntities + (i - 1)];
				}
			}
#else
			for (int i=1;i<=numEntities;i++)
			{
				dmax[i] = buffer[j].wert[i];
			}
#endif

			if (parent->mpi_rank == 0)
			{
				for (int i = 1; i <= numEntities; i++)
				{
					fprintf(outputfile,"\t%g",dmax[i]);
				}
			}
			delete[] dmax;

			if (parent->mpi_rank==0) 
			{
				fprintf(outputfile,"\n");
			}
		}

		if (parent->mpi_rank==0)
		{
			fclose(outputfile);
			(*WSstream) << "(" << (long)((bufcnt)*((long)4*sizeof(long)+(numEntities)*sizeof(buffer[1].wert[1])))
			<<" of "<<(long)buffersize<<" buffer-bytes)"<<"[end-entry]";
		}

		bufcnt=0;
		SetSignalAbfangfunktion();
	}

}

void SpecialOutputFileClass::addValues(long me, long se,long tg,long mt,long jr,bool bInSubStepsLoop,int actSubStep)
{
	if (incrementBufCntNextTime) 
	{
		bufcnt++;
		buffer[bufcnt].j=jr;
		buffer[bufcnt].m=mt;
		buffer[bufcnt].t=tg;
		buffer[bufcnt].s=se;
		buffer[bufcnt].min = me;
		incrementBufCntNextTime=false;
	}

	int startindex=1;
	for (std::vector<SpecialOutputEntityClass*>::iterator iter=entities.begin();iter<entities.end();iter++)
	{
		startindex=(*iter)->addValues(buffer,bufcnt,startindex,bInSubStepsLoop,actSubStep);
	}

	if (actSubStep==0)
	{
		if (bufcnt>= bufmax) 
		{
			flushbuffer();
		}
		incrementBufCntNextTime=true;
	}

}


//////////////////////////////////////////// SpecialOutputClass itself

SpecialOutputClass::SpecialOutputClass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
						double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
						char *steuername,int LogLevel,int progaus,
						double p_min_to_int_end,int p_numsubb,int *p_listsubb,
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
						int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int substepcnt):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine, ste, tge, mte, jre, p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
	v_mpiregioncnt,p_sendHalo,p_recvHalo,p_doSpinUp, p_routingmode, p_vHRU,1,false),
numfiles(0),iNumSubSteps(substepcnt)
{
	try
	{
        retvalue=0;

		if (do_it == 0) return;
		/* reservieren von Speicher fuer die Dateinamen */

		if (MPIREG==NULL)
		{
			if (ist_vorhanden(_MPIRegions_))
				MPIREG=gridinfo[getindex(_MPIRegions_)].pointer;
			else
				MPIREG=NULL;
		}

		if (MPIHLP==NULL)
		{
			if (ist_vorhanden(_MPIHelper_))
				MPIHLP=gridinfo[getindex(_MPIHelper_)].pointer;
			else
				MPIHLP=NULL;
		}


#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif
		// ToDo's: 1. looking for the section [special_output] in the control file -> is done in teilmodellklasse's constructor. 
		//	control file is open and positioned at the correct line already if we are here
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		
		std::string row(tstr);

		su.stripInt(row,"numfiles",&numfiles);

		// 2. for each output file: creating an instance of the class SpecialOutputFileClass
		for (int i=1;i<=numfiles;i++)
		{
			int depth=1;
			tstr=fgets(tstr,_stringsize_,steuerdatei);
			std::string record(tstr);
			while (record.find("}")==std::string::npos || depth>0)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				if (std::string(tstr).find("{")!=std::string::npos) depth++;
				if (std::string(tstr).find("}")!=std::string::npos) depth--;
				record.append(std::string(tstr));
			}
			files.push_back(new SpecialOutputFileClass(this,WSstream,record,gridinfo,xll,yll,cellsize,ncols,nrows,mpistartrow,mpiendrow,DefaultOutputDirectory,this->iNumSubSteps));
		}


	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	{
		retvalue=104;
	}
}

void SpecialOutputClass::run(long me,long se,long tg,long mt,long jr,bool bInSubStepsLoop,int actSubStep)
{
	// add values to buffer by calling SpecialOutputFileClass's addValues routine
	for (std::vector<SpecialOutputFileClass*>::iterator iter=files.begin();iter<files.end();iter++)
	{
		(*iter)->addValues(me,se,tg,mt,jr,bInSubStepsLoop,actSubStep);
	}
}

void SpecialOutputClass::flushbuffers()
{
	// call the flushbuffer method of all instances of the class SpecialOutputFileClass
	for (std::vector<SpecialOutputFileClass*>::iterator iter=files.begin();iter<files.end();iter++)
	{
		(*iter)->flushbuffer();
	}
}

void SpecialOutputClass::ClearMemory()
{
	// destroy all instances of the SpecialOutputFileClass -> the destructors of the vectors elements are called automatically
	for (std::vector<SpecialOutputFileClass*>::iterator iter = files.begin(); iter<files.end(); iter++)
	{
		delete *iter;
	}
	files.clear();
	modelobject::ClearMemory();
}

/* ------------- Ende der Definition der Teilmodell-Klasse ---------------------------- */

