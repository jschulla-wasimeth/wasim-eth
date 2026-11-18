
/*
(c) 1995-2016 by Joerg Schulla

Autor      : Joerg Schulla
more version history: see comments in SVN

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "teilmod.h"   /* Definition der Teilmodelle */
#include "wasim.h"

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


/* '''''''''''''''''''' Beginn der Expositionskorrekturroutinen ''''''''''''''''''''' */

expokorrklasse::expokorrklasse(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double breitearg,double laengearg,double bezugsmeridianarg,
	char *bezeichnung,char *steuername,int LogLevel,int progaus,
	double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax,int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode, p_vHRU),
hmax(0),gewicht(0),sp(0),spa(0),spe(0),strlsumme(0),cossunazim(0),korrparameter(0),
EXKORR(0),IEFFMAX(0),TKORR(0),ASP(0),SL(0),SSD(0),SSDABS(0),dhm(0),SHAPE(0),T(0),Tnacht(0),Isenk(0),spreiz(0),
indexx(0),indexy(0),dx(0),dy(0),vx(0),vy(0),
optim(0),dz(0),cos_za(0),EXKORRwritegrid(0),ieffmaxwritegrid(0),SHAPEwritegrid(0),TKORRwritegrid(0),
exkorrindex(0),ieffmaxindex(0),shapeindex(0),tnachtindex(0),tindex(0),tkorrindex(0),
cn(0),count(0),split(0),split_temp(0),EXKORRgridname(0),SHAPEgridname(0),TKORRgridname(0), ssdIsAbsolute(0)
{
	try
	{
		retvalue=0;
		int i,j,k;

		if (MPIREG == NULL)
		{
			if (ist_vorhanden(_MPIRegions_))
				MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
			else
				MPIREG = NULL;
		}
		if (MPIHLP == NULL)
		{
			if (ist_vorhanden(_MPIHelper_))
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			else
				MPIHLP = NULL;
		}


		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it != 0)
		{
#ifdef __TIME__
			(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
			//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif
			EXKORRgridname     = newstring(_stringsize_+1,1);
			SHAPEgridname      = newstring(_stringsize_+1,1);
			TKORRgridname      = newstring(_stringsize_+1,1);

			/* reading names and parameters */
			fscanfln(steuerdatei,"%d",&korrparameter    ,"methode of correction for "  ,tstr);
			fscanfln(steuerdatei,"%s",TKORRgridname     ,"name of the TKORR-Grid for"  ,tstr);
			fscanfln(steuerdatei,"%d",&TKORRwritegrid   ,"Writegrid EXKORR for "       ,tstr);
			fscanfln(steuerdatei,"%lf",&spreiz           ,"scaling factor for Tkorr for " ,tstr);
			fscanfln(steuerdatei,"%s",EXKORRgridname    ,"name of the EXKORR-Grids for",tstr);
			fscanfln(steuerdatei,"%d",&EXKORRwritegrid  ,"Writegrid EXKORR for "       ,tstr);
			fscanfln(steuerdatei,"%s",SHAPEgridname     ,"name of the SHAPE-Grid for"  ,tstr);
			fscanfln(steuerdatei,"%d",&SHAPEwritegrid   ,"Writegrid SHAPE for "        ,tstr);
			fscanfln(steuerdatei,"%d",&count            ,"number if timesteps for "    ,tstr);
			fscanfln(steuerdatei,"%d",&split            ,"number of splittings for "   ,tstr);

			/* als Korrekturparameter gelten nur die 1 (nur Strahlungskorrektur)
		und die 2 (auch Temperaturkorrektur) bzw. die 3 und 4, jeweils mit Abschattung */
			if ((korrparameter<1)||(korrparameter>4)) korrparameter=2;

			sr = mpistartrow;
			er = mpiendrow;
			if ((mpi_rank == 0) && (korrparameter <= 2))
			{
				sr = 1;
				er = nrows;
			}

			if (fclose(steuerdatei)==EOF)
			{
				(*WSstream) << "FATAL: error closing control file "<<steuerdateiname<<"[end-entry]";
				throw(11);
			}

			dhm    = gridinfo[getindex(_Hoehenmodell_)].pointer;
			ASP    = gridinfo[getindex(_Exposition_)].pointer;
			SL     = gridinfo[getindex(_Gefaelle_)].pointer;

			if (doSpinUp!=0 && HRU==NULL)
			{
				HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
				HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
			}

			/* SSD wird extra gesucht, da das Programm nicht abbrechen soll, wenn SSD nicht gefunden wird, getindex aber
		bei erfolgloser Suche das Programm mit einer Fehlermeldung abbricht */
			

//			i=1;
//			while (strcmp(gridinfo[i].identification,_Sonnenscheindauer_)!=0)
//			{
//				i++;
//				if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL)) break;
//			}
//
//			if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
//			{
//				SSD = NULL;
//				ssdindex = -1;
//			}
//			else
//			{
//				SSD = gridinfo[i].pointer;
//				ssdindex = i;
//			}

			/* auch die Temperatur soll nur verwendet werden, wenn sie als Eingabegroesse interpoliert wird */
		} /* if (doit!=0) */

		if (ist_vorhanden(_Sonnenscheindauer_))
		{
			SSD = gridinfo[getindex(_Sonnenscheindauer_)].pointer;
		}
		else
		{
			ssdabsindex = ist_vorhanden(_SunshineDurationAbsolute_);
			if (ssdabsindex > 0)
			{
				SSDABS = gridinfo[getindex(_SunshineDurationAbsolute_)].pointer;
				// create new SSD grid for ETPETR and other sub models
				ssdIsAbsolute = 1;
			}
			else
			{
				SSD = SSDABS = NULL;
			}
		}

		/* die Suche nach T wird auch ausgefuehrt, wenn doit==0, weil
	nachher Tnacht gebildet werden muss, wofuer T genutzt wird */
		i=1;
		while (strcmp(gridinfo[i].identification,_Temperatur_)!=0)
		{
			i++;
			if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL)) break;
		}

		if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
		{
			i=1;
			while (strcmp(gridinfo[i].identification,_Temperatur_14_)!=0)
			{
				i++;
				if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL)) break;
			}
			if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
			{
				T = NULL;
				tindex = -1;
			}
			else
			{
				T = gridinfo[i].pointer;
				tindex = i;
			}
		}
		else
		{
			T = gridinfo[i].pointer;
			tindex = i;
		}

		if (do_it!=0)
		{
			/* wo ist Platz fuer EXKORR  im array gridinfo? suchen des ersten freien Zeigers */
			i=1;
			while (gridinfo[i].identification!=NULL)
			{
				i++;
				if (i==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no free grid left"<<gridmaxanzahl<<"[end-entry]";
					throw(12);
				}
			}

			EXKORR = matrix(mpistartrow, mpiendrow, 1, ncols,nrows); // will only be written for actual mpi rank
			initialize_grid(i, EXKORRgridname, EXKORR, NULL, 0, _Expositionskorrektur_, EXKORRwritegrid, mpistartrow, mpiendrow);
			for (k = mpistartrow; k <= mpiendrow; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					EXKORR[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			exkorrindex = i;

			i++;
			SHAPE  = matrix(sr,er,1,ncols,nrows); // will be written in rank 0 for all rows and then scattered to other ranks, sowe need to reserve memory for sr to er and not mpistartrow to mpiendrow
			initialize_grid(i, SHAPEgridname, SHAPE, NULL, 0, _Schatten_, SHAPEwritegrid, sr, er);
			for (k = sr; k <= er; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					SHAPE[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			shapeindex = i;

			i++;
			IEFFMAX = matrix(mpistartrow, mpiendrow, 1, ncols, nrows); // will also be written for actual rank only
			initialize_grid(i, _IEFF_MAX_, IEFFMAX, NULL, 0, _IEFF_MAX_, 0, mpistartrow, mpiendrow);
			for (k = mpistartrow; k <= mpiendrow; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					IEFFMAX[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			ieffmaxindex = i;

			cn=count; /* das erste Mal soll auf jeden Fall die Abschattung berechnet werden */

			/* Suchen der groessten Hoehe, um den Schattensuchalgorithmus abbrechen zu koennen,
		wenn die Gerade ueber hmax steigt, und immer noch kein Schatten berechnet wurde */

			double globalhmax = nodata;
			hmax=-1e10;
			for (i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j = 1; j <= ncols; j++)
				{
					if ((zn[i][j] != nodata) && (dhm[i][j] > hmax) && (MPIREG[i][j]==mpi_rank))
					{
						hmax = dhm[i][j];
					}
				}
			}
#ifdef _MPI
			MPI_Reduce(&hmax, &globalhmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Bcast(&globalhmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			hmax = globalhmax;
#endif

			/* wenn T korrigiert werden soll,dann wird nun ein Array eingerichtet, welches die
		korrigierte Temperatur aufnimmt. Die Module zur Schmelze, Verdunstung usw. muessen
		dann diese Temperatur suchen! */
			if (((korrparameter==2)||(korrparameter==4))&&(T))
			{
				/* wo ist Platz fuer Tkorrt  im array gridinfo? suchen des ersten freien Zeigers */
				i=1;
				while (gridinfo[i].identification!=NULL)
				{
					i++;
					if (i==gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no free grid left"<<gridmaxanzahl<<"[end-entry]";
						throw(12);
					}
				}
				TKORR  = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
				initialize_grid(i, TKORRgridname, TKORR, NULL, 0, _Korrigierte_Temperatur_, TKORRwritegrid, mpistartrow, mpiendrow);
				for (k = mpistartrow; k <= mpiendrow; k++)
				{
					for (j = 1; j <= ncols; j++)
					{
						TKORR[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				tkorrindex = i;
			}
			else
			{
				TKORR=NULL;
				tkorrindex=-1;
			}
			if (EXKORRgridname ) free(EXKORRgridname ); EXKORRgridname =NULL;
			if (SHAPEgridname  ) free(SHAPEgridname  ); SHAPEgridname  =NULL;
			if (TKORRgridname  ) free(TKORRgridname  ); TKORRgridname  =NULL;

		}   /* if (do_it!=0), der Rest muss gemacht werden, weil ind er etpetr-Routine bei Tagesschritten Tnacht gesucht wird */

		/*
	Es muss ein zweites Temperaturgrid initialisiert werden
	wenn InV > 720 min ist, da dann im Verdunstungsmodell sowie im Schneemodell
	die Temperatur eines Intervalles aufgrund empirischer Beziehungen zwischen
	der Differenz der am Saettigungsdampfdruck gewichteten Tages- und Nachttemperaturmittelwerte
	einerseits sowie der Hoehe, der Jahreszeit und dem Verhaeltnis aus interpolierter/modellierter Strahlung
	zur maximal moeglichen Strahlung andererseits fuer Tag und Nacht modifiziert wird.
	Das macht die Verwendung zweier Temperaturgrids noetig, da im Verdunstungsmodell die Tagtemperatur
	aufgrund der bereits im Expositionskorrekturmodul korrigierten Temperatur erhoeht
	und die Nachttemperatur von den hier nicht korrigierten Ausgangstemperaturen abgesenkt werden.
	Da ausserdem beide Temperaturen auch im Schneemodell genutzt werden, muessen beide Temeraturen
	gespeichert werden. Wenn das Zeitintervall <720 min. ist, dann wird kein zweites T-Grid
	erzeugt. Die Adresse von Tnacht wird dann auf T gesetzt.
	*/

		if (InV > 720.0)
		{
			/* nur wenn T auch existiert, kann und soll es dupliziert werden, ansonsten werden
		Verdunstungs- und Schneemodell ja nicht ausgefuehrt */
			if (T)
			{
				/* wo ist Platz fuer Tnacht  im array gridinfo? suchen des ersten freien Zeigers */
				i=1;
				while (gridinfo[i].identification!=NULL)
				{
					i++;
					if (i>=gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no free grid left"<<gridmaxanzahl<<"[end-entry]";
						throw(12);
					}
				}
				Tnacht  = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
				initialize_grid(i, __Nacht_Temperatur__, Tnacht, NULL, 0, _Nacht_Temperatur_, 0, mpistartrow, mpiendrow);
				for (k = mpistartrow; k <= mpiendrow; k++)
				{
					for (j = 1; j <= ncols; j++)
					{
						Tnacht[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				tnachtindex = i;
			}
		}
		else
		{
			Tnacht = T;
			tnachtindex = tindex;
		}

		// SSD ahouls be created regardles is expocorr runs or not since ot is required in evaporation and unsatzon as well
		if (ssdIsAbsolute == 1)
		{
			// create SSD grid
			i++;
			SSD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows); // will also be written for actual rank only
			initialize_grid(i, _Sonnenscheindauer_, SSD, NULL, 0, _Sonnenscheindauer_, 0, mpistartrow, mpiendrow);
			for (k = mpistartrow; k <= mpiendrow; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					SSD[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			ssdabsindex = i;
		}



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


double **expokorrklasse::GetExpoKorrAdress()
{
	if (do_it==1)
	{
		return EXKORR;
	}
	else
	{
		return NULL;
	}
}


// methid convertAbsoluteToRelativeSSD
// this method is called only once per time step. 
// it converts the sunshine duration from absolute values in minutes per time step into relative values in 1/1.
// It should work for all comibinations of sunset and sunrise with respect to the time step start and end as well as in polar day and polar night (northern and southern hemisphere) 
// and also for time zones where the sun rises after the sun set during a 24 hours day.
void expokorrklasse::convertAbsoluteToRelativeSSD(int minute,int hour, int day, int month, int houroffset)
{
	
	double sunr =0, suns = 0;

// The relative sunshine duration is valid for all following time steps, even if there are 24 hours as sub steps within a day (daily time step)
// However, this is handled in the algorithms later. They also consider sunrise and sunset, so we can concentrate here on the pure conversion
// call the sunrise function here 
	sunrise(day, month, breite, laenge, &sunr, &suns);
	sunr = fmod(sunr + 24.0, 24.0); // fixes the sunrise between 0 ... 24
	suns = fmod(suns + 24.0, 24.0); // fixes the sunset between 0 ... 24 

	bool polarday = false;
	bool polarnight = false;

	// for areas with no sunset or no sunrise within a day, suns will be returned with 0.0 and sunr with 24.0 by the sunrise function. This may be polar regions during polar day and polar night
	// so lets check
	if ((suns == 0.0) && (sunr == 24.0))
	{
		if ((breite > 66.0) && ((month >= 9) || (month <= 3))) 		 { 
			// polarnight northern hemisphere
			polarnight = true;
		}
		else if ((breite > 66.0) && ((month < 9) || (month > 3))) 	 { 
			// polarday northern hemisphere
			polarday = true;
		}
		else if ((breite < -66.0) && ((month >= 9) || (month <= 3))) { 
			// polarday southern hemisphere
			polarday = true;
		}
		else if ((breite < -66.0) && ((month < 9) || (month > 3)))   { 
			// polarnight southern hemisphere
			polarnight = true;
		}
		else {
			// something is wrong, since nowhere on earth the sun shines either 0 or 24 hours outside the polar regions, which begin above 66.33 degree north or below 66.33 degree south
		}
	}

	int maxdaylightminutes;
	if (polarday) {
		maxdaylightminutes = InV;
	}
	else if (polarnight) {
		maxdaylightminutes = 0;
	} 

	sunr *= 60.0; // conversion into minutes
	suns *= 60.0; // conversion into minutes
	double tstart = (hour * 60) + houroffset - InV + minute; // Start off time subste in Minutes since midnight
	double tend = tstart + InV; // Time step end in minutes since midnight

	// now calculate the interval lenght with the sun over the horizon which gives then the maximum sunshine duration in minutes for the interval
	if (sunr < suns)
	{ // normal case
		if ((tstart < sunr) && (tend < sunr)) {
			maxdaylightminutes = 0; // interval starts and ends before sunrise
		} 
		else if ((tstart<sunr) && (tend>sunr) && (tend<suns)) {
			maxdaylightminutes = tend - sunr; // sun rises within the time step but doesn't set within the same time step
		}
		else if ((tstart>sunr) && (tend>sunr) && (tend<suns)) {
			maxdaylightminutes = InV; // sunset before the time step starts and sunrise after the time step ends, i.e. gright day always
		}
		else if ((tstart>sunr) && (tend>sunr) && (tend>suns)) {
			maxdaylightminutes = suns - tstart; // sunriose before the time step starts and sunset during the time step
		}
		else if ((tstart>suns) && (tend>suns)) {
			maxdaylightminutes = 0; // sunrise and sunset befor the time step starts, everything is dark
		}
		else if ((tstart<sunr) && (tend>suns)) {
			maxdaylightminutes = suns - sunr; // e.g. for daily time steps, the time between sunrise and sunset if taken as maximum possible sunshine duration in minutes
		}
	}
	else
	{ // sunset before sunrise, can be due to wide time zones
		if ((tstart < suns) && (tend < suns)) {
			maxdaylightminutes = InV; // interval starts and end before sunset
		}
		else if ((tstart<suns) && (tend>suns) && (tend<sunr)) {
			maxdaylightminutes = suns - tstart; // sunset is within the time step
		}
		else if ((tstart>suns) && (tend>suns) && (tend<sunr)) {
			maxdaylightminutes = 0; // dark durgin the entire time step
		}
		else if ((tstart<sunr) && (tend>sunr) && (tstart>suns)) {
			maxdaylightminutes = tend - sunr; // sunrise within the timestep
		}
		else if ((tstart > sunr) && (tend > sunr)) {
			maxdaylightminutes = InV; // complete time step in day light
		}
		else if ((tstart<suns) && (tend>sunr)) {
			maxdaylightminutes = (suns - tstart) + (tend - sunr); // light at the beginning of the time step and at the  end of the time step is added up
		}

	}



	// the relative duration is simply the fraction of the observed value on the theoretical maximum (with limits of 0 and 1).

	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
			{
				SSD[i][j] = (maxdaylightminutes > 0) ? SSDABS[i][j] / maxdaylightminutes : 0.0; // first, we only copy the values for development testing the correct calls
				if (SSD[i][j] > 1.0) SSD[i][j] = 1.0;
				
			}
		}
	}

}


/* Ausfuehren der expositions-, gefaelle- und
abschattungsabhaengigen Korrektur fuer die direkte Sonnenstrahlung
sowie für die Temperatur */



void expokorrklasse::run(long me, long se,long tg,long mt,long jr)
{
	int i,j,k,l,m;
	int nrl,nrh,ncl,nch; /* werden jeweil für die untere und obere Grenze der Schleifen zum
						Durchlaufen der Zeilen (nrl, nrh) und der Spalten (ncl, nch) berechnet */
	double cosZ;
	double inta,inte;

	/* auch wenn do_it == 0 ist, muss die Temperatur, wenn T!=NULL und wenn InV>720.0 ist,
	in Tnacht dupliziert werden, da das Verdunstungs- sowie das Schneemodell beide bei
	zeitlicher Aufloesung > 720 min eine Unterteilung in Tag und Nachttemperatur vornehmen.
	Erst danach kann mit do_it==0 abgebrochen werden */

	if (T) 
	{
		if (InV>720.0)
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
					{
						Tnacht[i][j] = T[i][j];
					}
				}
			}
		}
	}

	if (ssdIsAbsolute)
	{
		this->convertAbsoluteToRelativeSSD(me,se,tg,mt,hour_offset);
	}

	if (do_it == 0) return;

	/* wenn z.B. die Faktorbestimmumg nur alle 3 Stunden getan werden soll,
	dann ist count==3, es werden jeweils 2 Durchgaenge an dieser Stelle abgebrochen,
	nur jeder dritte Durchgang wird durchgefuehrt */
	if (cn<count)
	{
		cn++;
		return;
	}
	cn=1;

	/* Ausgabe einer Nachricht */
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);
	gewicht=1.0/(double)split;
	strlsumme=0.0;

	/* Berechnen von sunazimuthangle, sunaltitudeangle, sunzenithangle fuer diese Stunde
	Angaben der Ergebnisse in rad, der Eingangsdaten aber in Grad! */
	sunrise(tag,monat,breite,laenge,&sunr,&suns);
	/*  if (suns<0) suns = 24.0; */
	/* wenn es keinen Sonnenaufgang gibt, wird der Beginn des Tages angenommen als Sonnenaufgang */
	/*  if (sunr<0) sunr = 0.0; */
	sunr = fmod(sunr + 24.0, 24.0); /* Sonnenaufgang zwischen 0 ... 24 Uhr fixieren */
	suns = fmod(suns + 24.0, 24.0); /* Sonnenuntergang ebenfalls fixieren */

	getsuncoordinates(tag, monat, breite, laenge,(sunr+suns)/2.0, &sunazimuthangleAtNoon, &sunaltitudeangleAtNoon);

	/* wenn die Temperatur korrigiert werden soll, dann mussen die T-Werte in TKORR geschrieben werden,
	weil etpetr und snowmod in diesem Fall nach diesem Grid suchen anstatt nach T. TKORR wird aber nicht
	berechnet (unten), wenn SSD fehlt, so dass es sonst immer 0 Grad haette. */
	if (TKORR)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
				{
					//        if (zn[i][j]!=nodata) 
					TKORR[i][j] = T[i][j];
				}
				else
				{
					TKORR[i][j] = nodata;
				}
			}
		}
	}
	
	/* Nullsetzen der Korrekturfaktoren erfolgt explizit, da die Zuweisung in der
	Schleife zur Korrekturfaktorermittlung incrementell erfolt (es wird dann fuer
	jedes Teilintervall eine um den Korrekturfaktor fuer dieses Zeitintervall korrigierte
	Strahlung addiert, so dass zum Intervallende genau <split> Teilintervall-Ergebnisse
	addiert sind. Der Wert in EXKORR entspricht dann dem Korrigierten Stundenwert der Strahlung
	(Korrekturfaktoren mit der astronomisch und topographisch moeglichen STrahlung
	gewichtet). Ganz am Ende wird dann durch eine separat ermittelte unkrrigierte
	Intervallsumme der STrahlung geteilt, um so den effektiven Korrekturfaktor zu erhalten */

	for (int i=mpistartrow;i<=mpiendrow;i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
			//      if (zn[i][j]!=nodata)
			{
				EXKORR[i][j]=0.0;
				IEFFMAX[i][j]=0.0;
			}
		}
	}

	/* Da fuer mehrere Intervalle die Strahlung unterschiedlich gross sein kann, muessen die
	Korrekturfaktoren mit den jeweils astronomisch moeglichen Strahlungen
	gewichtet werden. Es wird also das Integral der korrigierten maximal
	moeglichen Strahlungen errechnet, anschliessend durch das Integral der
	unkorrigierten maximal moeglichen Strahlung geteilt (strlsumme) und so
	der mittlere Korrekturfaktor erhalten. */

	/* Teilintervall am Anfang der Integration (Mitte des 1. der <split> Teilintervalle) */
	spa = (se+hour_offset + (me / 60.0))-InV/60.0*(1-0.5/(double)split)+timekorr;

	/* Teilintervall am Ende der Integration (Mitte des letzten der <split> Teilintervalle*/
	spe = (se+hour_offset + (me / 60.0))+0.01+timekorr;

	split_temp = InV/60.0/(double)split; /* Dauer eines Teilintervalles in h */

	l=0;

	for (int sp=(double)spa;sp<(double)spe;sp+=gewicht*(InV/60.0))
	/* sp ist die Zeit (Bruchteile voller Stunden erlaubt und noetig) des aktuellen Teilintervalles,
	korrigiert mit der Zeitverschiebung der Messzeit gegen die gesetzliche Zeit (darf nur Winterzeit sein);
	so ist spa bei einem Zeitintervall von 2 Stunden und split=4 (jede halbe Stunde neu berechnen)
	nach obiger Berechnung fuer se=10 (die an die Routine uebergebene Stundennummer [fuer das Ende des Intervalls])
	10 - 120/60 * (1 - 0.5/4) + korr = 10 - 2 * (1-0.25) = 10 - 2 + 0.5 = 8.5.
	Das erste Teiintervall ist also 0.5 Stunden lang und endet um 8:30.
	Als Increment wird jeweils die mit 1/<split> gewichtete Laenge eines
	Intervalles in Stunden (hier 0.25*120/60 = 0.5 addiert */
	{
		inta = sp - split_temp / 2.0;
		inte = sp + split_temp / 2.0;

		/* Wie lange ist es am Anfang eines Teilintervalls noch dunkel und am Ende schon dunkel? */

		/* dark_a */
		if (sunr < inta)
		{
			dark_a = 0.0;
		}
		else
		{
			if (sunr < inte)
			{
				dark_a = sunr - inta;
			}
			else /* Stunde komplett vor Sonnenaufgang */
			{
				dark_a = 0.0;
			}
		}
		/* dark_e */
		if (suns > sunr)
		{ /* der Normalfall */
			if ((suns>inte) && (sunr < inte))
			{
				dark_e = 0.0;
			}
			else
			{
				if ((suns > inta) && (suns < inte))
				{
					dark_e = inte - suns;
				}
				else
				{
					dark_e = split_temp;
				}
			}
		}
		else
		{ /* der Fall, dass die Sonne erst nach Mitternacht gesetzl. Zeit untergeht */
			if (suns > inte)
			{
				dark_e = 0.0;
			}
			else
			{
				if (suns > inta)
				{
					dark_e = inte - suns;
				}
				else
				{
					dark_e = sunr > inte ? split_temp : 0.0;
				}
			}
		}

		Anteil_dunkel = (dark_a + dark_e) / split_temp;
		/* Wenn die Sonne im Intervall erst unter, dann aufgeht, muss zum Anfang und Ende
	des Intervalls das Komplement zu 1 von dark_a und dark_e gebildet und addiert werden
	Anteil_hell=(1-dark_a)+(1-dark_e) = 2-(dark_a+dark_e) = 2.0 - Anteil_dunkel
	was für Anteil_dunkel zu Anteil_dunkel = 1-Anteil_hell = 1 - 2 + Anteil_hell = Anteil_dunkel-1 fuehrt */
		if (Anteil_dunkel > 1.0)
		{
			Anteil_dunkel--;
		}
		/* wenn die Sonne nur ganz kurz auftaucht, nuetzt es gar nichts!
	Da kann sie auch gleich untergetaucht bleiben */
		if (Anteil_dunkel > 0.95)
		{
			Anteil_dunkel = 1.0;
		}
		Anteil_hell = 1.0 - Anteil_dunkel;

		/* wenn das aktuelle Berechnungsintervall komplett (zu 95%) vor Sonnenaufgang
	oder nach Sonnenuntergang liegt, wird keine Strahlung beruecksichtigt */
		if (Anteil_hell < 0.05)
		{
			nostrl = 1;
		}
		else
		{
			nostrl = 0;
		}

		if (!nostrl)
		{
			/* wo steht die Sonne zur Mitte des hellen Abschnittes des aktuellen
		Berechnungsintervalles? Die Stunde wurde schon oben bei spa und spe mit dem hour_offset verrechnet*/
			getsuncoordinates(tag, monat, breite, laenge,
				/* Mitte des helle Bereichs des Teilintervalls */
				/* = sp - split_temp/2.0 + dark_a + split_temp/2.0*Anteil_hell */
				sp + split_temp / 2.0*(Anteil_hell - 1.0) + dark_a,
				&sunazimuthangle, &sunaltitudeangle);
			sunzenithangle = pi / 2 - sunaltitudeangle;
			cosZ = cos(sunzenithangle);

			// new for testing smoother temperature offsets 2025-11-09
			sunweighingfactor = sunaltitudeangle / sunaltitudeangleAtNoon;
			if (sunweighingfactor >= 0)
			{
				sunweighingfactor = sqrt(sunweighingfactor);
			}
			else
			{
				sunweighingfactor = 0.0;
			}

			if (sunaltitudeangle > 0.017) /* 1 Grad */
			{
				/* Berechnung der bei wolkenfreiem Himmel und klarer Atmosphaere auf einen
		senkrecht zur Einstrahlungsrichtung liegenden Quadratmeter
		auftreffender Strahlungsenergie [Isenk] = Wh/m^2 */
				/* cos(sunazimuthangle) wirkt hier auf den Extinktionskoeffizienten!  */
				/* auf 1000m bezogen, 0.9 = (1-1000.0/10000.0) */

				Isenk = split_temp*Anteil_hell*1367.0 * pow(0.86, 0.9 / cosZ);

				/*        Isenk *= cos(sunzenithangle); */ /* beziehen auf waagerechte Flaeche */

			}
			else
			{
				Isenk = 0.0; /* Nachts gibt es keine kurzwellige Einstrahlung */
				sunazimuthangle = 0.0;
				cosZ = 1.0;
				sunweighingfactor = 0.0;
			}

		}
		else
		{
			Isenk = 0.0;
			sunazimuthangle = 0.0;
			cosZ = 1.0;
			sunweighingfactor = 0.0;
		}

		strlsumme += (Isenk*cosZ);

		/* SHAPE auf 0 setzen */
		for (int i = sr; i <= er; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && ((MPIREG[i][j] == mpi_rank) || (mpi_rank == 0 && korrparameter<=2)) && !(doSpinUp && HRUMaster[i][j] == nodata))
				// if (zn[i][j]!=nodata)
				{
					SHAPE[i][j] = 0.0;
				}
			}
		}


		/* Berechnung der Abschattung will be skipped if in spinUp-mode and will only be done in mpi_rank 0 and then scattered to all other ranks */
		if ((!nostrl) && ((korrparameter == 1) || (korrparameter == 2)) && !doSpinUp)
		{
			if (mpi_rank == 0)
			{
				/* Vorzeichen der Geradenanstiege in x- und y-Richtung */
				if (sunazimuthangle <= Pi / 2.0)
				{
					vx = 1;
					vy = 1;
				}
				else
				{
					if (sunazimuthangle <= Pi)
					{
						vx = 1;
						vy = -1;
					}
					else
					{
						if (sunazimuthangle <= Pi*1.5)
						{
							vx = -1;
							vy = -1;
						}
						else
						{
							vx = -1;
							vy = 1;
						}
					}
				}


				if (vx > 0) 
				{ 
					ncl = 1;     
					nch = ncols + 1; 
				}
				else 
				{ 
					ncl = ncols; 
					nch = 0; 
				}

				if (vy > 0) 
				{ 
					nrl = er; 
					nrh = 0; 
				}
				else 
				{ 
					nrl = 1;     
					nrh = er + 1; 
				}

				// in x Richtung kleiner 45 Grad Neigung 
				if (fabs(cos(sunazimuthangle)) <= 0.707106781)
					// x-schritt ist immer 1, y-schritt wird sukzessive nachgefuehrt 
				{
					dy = (int)fabs((ncols*cos(sunazimuthangle) / sin(sunazimuthangle)));
					dx = ncols;
					dz = fabs((double)cellsize / sin(sunazimuthangle))*tan(sunaltitudeangle);

					//#pragma omp parallel for schedule(runtime) private(i,j)
					for (int i = nrl; i != nrh; i -= vy)
					{
						for (int j = ncl; j != nch; j += vx)
						{
							if (i != nrh + vy)
							{
								// only the last (or first) column should be processed
								if (vx > 0) { j = ncols; }
								else { j = 1; }
							}
							// must be done if it's a nodata cell or not --> follow the ray and you will eventually stumble across a valiud cell
							// this algorithm is NOT suitable for MPI parallelization, but probably for OpenMP
							//				(*WSstream) << "row " << i <<", col " << j << "[end-entry]";
							int fx = -dx / 2;
							int ii = i;
							int jj = j;
							double z = dhm[i][j];
							while ((ii != nrl + vy) && (jj != ncl - vx))
							{
								if (z > nodata) z -= dz;
								jj -= vx;
								fx += dy;
								if (fx >= 0)
								{
									fx -= dx;
									ii += vy;
								}
								if ((ii != nrl + vy) && (jj != ncl - vx))
								{
									if (zn[ii][jj] != nodata)
									{
										if (dhm[ii][jj] < z)
										{
											SHAPE[ii][jj] = 1;
										}
										else
										{
											SHAPE[ii][jj] = 0;
											z = dhm[ii][jj]; // equal or higher altitudes must be updated now in the comparison value
										}
									}
								}
							}
						}
					}
				}
				else
					// y-schritt ist immer 1, x-schritt wird nachgefuehrt 
				{
					dx = (int)fabs(nrows*sin(sunazimuthangle) / cos(sunazimuthangle));
					dy = nrows;
					dz = fabs((double)cellsize / cos(sunazimuthangle))*tan(sunaltitudeangle);

					//#pragma omp parallel for schedule(runtime) private(i,j)
					for (int i = nrl; i != nrh; i -= vy)
					{
						for (int j = ncl; j != nch; j += vx)
						{
							if (i != nrh + vy)
							{
								// only the last (or first) column should be processed
								if (vx > 0) { j = ncols; }
								else { j = 1; }
							}
							//				(*WSstream) << "row " << i <<", col " << j << "[end-entry]";
							int fy = -dy / 2;
							int ii = i;
							int jj = j;
							double z = dhm[i][j];
							//              while ((ii!=nrh+vy)&&(jj!=nch-vx)&&(z<hmax)&&(SHAPE[i][j]<0.5))
							while ((ii != nrl + vy) && (jj != ncl - vx))
							{
								z -= dz;
								ii += vy;
								fy += dx;
								if (fy >= 0)
								{
									fy -= dy;
									jj -= vx;
								}
								if ((ii != nrl + vy) && (jj != ncl - vx))
								{
									if (zn[ii][jj] != nodata)
									{
										if (dhm[ii][jj] < z)
										{
											SHAPE[ii][jj] = 1;
										}
										else
										{
											SHAPE[ii][jj] = 0;
											z = dhm[ii][jj]; // equal or higher altitudes must be updated now in the comparison value
										}
									}
								}
							}
						}
					}
				}
			}

#ifdef _MPI
			MPIScatterGrid(shapeindex);
#endif
		} /* ((!nostrl)&&((korrparameter==1)||(korrparameter==2))) */

#ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		/* expositions- und neigungsbedingte Strahlunsgkorrektur */
		double cosTheta;
#pragma omp parallel for schedule(runtime) shared(cosZ) private(i,j,cosTheta)
		for (i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) /* Spalten des Grids */
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
		        //if (zn[i][j]!=nodata)
				{
					if ((!nostrl) && (SHAPE[i][j] < 0.5))
					{
						cosTheta = getcosincidence(ASP[i][j] * pi / 180.0, SL[i][j] * pi / 180.0, sunazimuthangle, sunzenithangle);
						if (cosTheta <= 0)
						{
							cosTheta = 0.0;
						}
						else
						{
							if (SSD)
							{
								if (SSD[i][j] < 0.1)
								{
									cosTheta = cosZ; /* wenn keine Sonne scheint, soll keine Korrektur erfolgen */
								}
							}
						}
					}
					else
					{
						cosTheta = 0;
					}

					/* Korrektur fuer direkte Sonnenstrahlung:
						wenn Sonne unter Horizont      : 1 (keine Korrektur
						(alle gemessene STrahlung ist
						diffuse Strahlung)
						wenn Sonne Flaeche nicht trifft: 0 (keine Korrektur)
						wenn Sonne Flaeche trifft, aber
						SSD zu klein war          : 1 (keine Korrektur -> alles diffuse STrahlung)
						wenn Flaeche im Schatten liegt : 0 (keine Korrektur)
						wenn Sonne < 12 Grad           : cosTheta = 5
						sonst                          : echte Korrektur     */

					EXKORR[i][j] += (cosTheta > 5.0*cosZ ? 5.0*cosZ*Isenk : cosTheta*Isenk);

				}
			}
		}
//		stunde=se*10+l;
		writegridtofile(shapeindex, 1, sp + gewicht*(InV / 60.0)<(double)spe);
		l++;

	} /* Schleife fuer die <split> Teil-(zeit-)intervalle */

	/* Strahlungssumme aller Teilintervalle wieder wegdividieren,
	also den effektiven Faktor, das eigentliche Ergebniss dieses Teilmodells,
	errechnen */

#pragma omp parallel for schedule(runtime) private(i,j)
	for (i=mpistartrow;i<=mpiendrow;i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
			//      if (zn[i][j]!=nodata)
			{
				IEFFMAX[i][j] = EXKORR[i][j];
				EXKORR[i][j] = strlsumme > 0.0 ? EXKORR[i][j]/strlsumme : 0 ;
				/* nun die Temperaturkorrektur:
			Die interpolierte Temperatur wird dann korrigiert,
			wenn es Sonne gab. Als Korrektur wird der natuerliche Logarithmus
			des Strahlungskorrekturkoeffizienten verwendet, da EXKORR im Wertebereich
			zwischen 0 und (theoretisch) +Inf liegt sowie bei 1.0 seinen Normalwert hat.
			Die Logarithmen haben dann den Wertebereich von -INF ... +INF, wobei durch die
			Begrenzung des EXKORR-Faktors auf 5.0 die Korrekturwerte für die Temperatur
			im Bereich -INF ... +1.6K liegen. Wird für die Temperaturkorrektur
			noch ein unterer Grenzwert für EXKORR von 0.2 festgelegt, dann liegt
			der additive Temperaturkorrekturwert im Bereich +- 1.6 Kelvin.
			Der sich so ergebende Korrekturwert wird noch an der relativen Sonnenscheindauer
			gewichtet, um zu erreichen, dass bei wolkigem Himmel die Temperaturen nicht
			so stark schwanken wie bei klarem Himmel.
		*/

				double security_factor = (strlsumme>0.0)?1.0:0.0;
				if (((korrparameter==2)||(korrparameter==4))&&(TKORR)&&(T))
				{
					double ssdfact=0.5;
					if (SSD) ssdfact=SSD[i][j];
					/* am Anfang von expositionskorrektur::run wurde TKORR auf T gesetzt
						(damit etpetr auch dann gueltige Werte bekommt, wenn zwar korrekturparameter==2 oder 4,
						aber SSD==NULL ist, was ja bedeutet, dass die T-Korrektur nicht ausgefuehrt wird */
					if (EXKORR[i][j]>0.2)
					{
						TKORR[i][j] = (split==1) ? 
							T[i][j] + (spreiz * log(EXKORR[i][j]) * ssdfact * sunweighingfactor):
							T[i][j] + (spreiz * log(EXKORR[i][j]) * ssdfact);
					}
					else
					{
						/* diese Zeile erfasst auch den Fall, dass es Nacht oder bewölkt ist,
							weil dann SSD 0 ist. Bei Tageswerten wird hier nur die Tagesmitteltewmperatur
							korrigiert -1.609 ist ln(0.2)*/
						/* (1.0-nostrl) will be implemented as security factor: if SSD is >0 for any reason 
							(even if it is impossible by physics, one could have given a ssd value of > 0 in the input),
							then the temperature would be depleted in every dark interval */
						TKORR[i][j] = (split == 1) ?
							T[i][j] - (spreiz * 1.6094 * ssdfact * security_factor * sunweighingfactor):
							T[i][j] - (spreiz * 1.6094 * ssdfact * security_factor);
					}
				}
			}
		}
	}
	stunde = se;
	if (EXKORR) writegridtofile(exkorrindex,1);
	if (TKORR)  writegridtofile(tkorrindex,1);
//	if (SHAPE)  writegridtofile(shapeindex);
}

void expokorrklasse::ClearMemory()
{
	/* Tnacht muss freigegeben werden, wenn es erzeugt wurde -> (T!=NULL)&&(InV>720.0) */
	if (T)
	{
		if (InV>720.0)
		{
			free_matrix(gridinfo[tnachtindex].pointer, gridinfo[tnachtindex].nrl, gridinfo[tnachtindex].nrh, 1, ncols, nrows);
			delete_grid(tnachtindex);
		}
	}

	if (do_it != 0)
	{
		if (((korrparameter==2)||(korrparameter==4))&&(T))
		{
			if (TKORR)
			/* TKORR freigeben */
			{
				free_matrix(gridinfo[tkorrindex].pointer, gridinfo[tkorrindex].nrl, gridinfo[tkorrindex].nrh, 1, ncols, nrows);
				delete_grid(tkorrindex);
			}
		}

		free_matrix(gridinfo[shapeindex].pointer, gridinfo[shapeindex].nrl, gridinfo[shapeindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[exkorrindex].pointer, gridinfo[exkorrindex].nrl, gridinfo[exkorrindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[ieffmaxindex].pointer, gridinfo[ieffmaxindex].nrl, gridinfo[ieffmaxindex].nrh, 1, ncols, nrows);

		if (ssdIsAbsolute == 1)
		{
			free_matrix(gridinfo[ssdabsindex].pointer, gridinfo[ssdabsindex].nrl, gridinfo[ssdabsindex].nrh, 1, ncols, nrows);;
			delete_grid(ssdabsindex);
		}

		delete_grid(shapeindex);
		delete_grid(exkorrindex);
		delete_grid(ieffmaxindex);

	}
	modelobject::ClearMemory();
}

/* '''''''''''''''''''' Ende der Exkpositionskorrekturroutinen ''''''''''''''''''''' */

