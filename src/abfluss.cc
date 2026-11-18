/*   $Revision: 61 $ 
	$Date: 24.08.10 23:11 $  */

/*

(c) 1995-2007 by Joerg Schulla

Version: 1.0 (Juni 94) einfaches Routing mit Zeitkonstante und Speicherkonstante nur lin. r^2
			1.1 1996 mit r^2 und EV, jeweils linear und logarithmisch
			2.0 Februar 1997: variable Translationszeiten durch adaptiertes Translations-Retentions-Verfahren
							in Anlehnung an Ludwig, 1979 (FGMOD) aber bilanzrein, etwas einfachere Retentionsberechnung
							mit einem nichtlinearen Stufenmodell getrennt fuer Hauptbett und Vorlandabfluss bzw. Abfluss
							mit oder ohne Ueberflutung der Vorlaender.
							Ausserdem besseres Handling des Spendenfiles
							(es muss nun keine Spendendatei mehr vorhanden sein, wenn keine
							vorhanden ist, dann wird kein Bestimmtheitsmass berechnet, in der
							Spendendatei koennen nun die Spalten durch beliebig viele
							Trennzeichen (ausser \n) getrennt sein).
		3.0 19. Juni 1997: Beruecksichtigung von Ueberleitungen im Gebiet, Zu- und Ableitungen sowie von Speichern.
							externe Dateien enthalten den ABfluss immer in m^3/s (fuer Speicher: Inhalt in m^3)
						- Zuleitungen koennen exern von Datei (Stunden- oder Tageswerte) oder intern von
							einer ebenfalls deklarierten Ableitung aus einem anderen Teilgebiet uebernommen werden.
						- Ableitungen koennen oberhalb einer Schwelle mit einem bestimmten Anteil des die Schwelle uebersteigenden
							Abflusses bis zum maximalen Fassungsvermoegen der Ableitungseinrichtung vom modellierten Afluss
							vom Abfluss abgezweigt werden.
						- Speicher sollten als eigenes Teilgebiet fuer die Modellierung behandelt werden. Sie koennen als
							Start- bzw. Endpunkte von beliebig vielen Zu- und Ableitungen genutzt werden.
							Abgaberegeln: * Tabelle Inhalt <-> Abgabe
						- Routingregeln muessen nun nicht mehr alle OLs, ZLs, ALS und SPs in einer einzigen Zeile pro TG enthalten,
							sondern je OL, AL, ZL bzw. SP kann eine eigene Zeile genutzt werden, sofern das ganze TG ohne Leerzeile
							nacheinander beschrieben wird. Nach wie vor ist die erste Leerzeile das Erkennungszeichen dafuer, dass
							keine neue Routingvorschrift kommt.

		4.0 Maerz 1998: - Stofftransportalgorithmen parallel zum Wasserfluss integriert
						- Initialisierung der Speicher: neue Routingvorschrift derart:
							und SP <n> (file = <filename> , V0 = <ini-volume> , C0 = <ini-conc-1> ... <ini-conc-9> )

		4.1 16.6.98   : - irrigation amouts can now also got from reservoirs in a head basin. Reservoirs can now also placed in a head basin 
							also if there is no other routing (inflow, abstraction, tributary)
		
		5.0 04.08.98  : - until now routing was carried out for each tributary in a single routing channel, 
							now its possible to superpose tributaries prior to routing (they are flowing into another at the beginning of a routing channel)
							the Symbol used for this is SUMTRIB <n1>+<n2>+<n3>+... with n1 ... nx the codes of tributaries
							tanalys was already changed in order to generate an appropriate routing structure file

date       : march, 9 1999
version    : 5.0a hour-entry can be given for start or for end of time step via min_to_int_end

5.1 June 5, 1999
- It's now possible to save and restrore the contents of the runoff which is already in the runoff concentration
	process, as is true for the content of the flow time zones and their single linear reservoirs as well as for
	the content of all routing channels and their single linear reservoirs. Also the tracer concentrations are stored
	in the same way. The storing/restoring mechanism works on a file base. The filename is given as the last
	parameter in the [unsatzon_model] section of the control file and it's valid also for the [routing_model] section

date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.

more version history: see source safe comments

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>
#include <csignal>
#include <ctime>
#include <string>
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

/*~~~~~~~~~~~~~~~~~~~~ Beginn der Routing - routinen ~~~~~~~~~~~~~~~~~ */
void routingclass::riverflow_table_ini(riverflow_table_type *adresse,int nr)
{
	int i,j,z,ber,ind;
	riverflow_table=adresse;
	if (riverflow_table)
	{
		for (i=1;i<=OLs+SUMTRIBs;i++)
		{
			if (indexST[i])
			{
				for (j=1;j<=indexST[i][0];j++)
				{
					ind=zarray[indexST[i][j]];
					riverflow_table[ind].OL_index = i; 
					riverflow_table[ind].ROUT_index    = nr; 
					riverflow_table[ind].upstream_area = AE[i];  /* area of the first of the SUMTRIBs is already the entire basin area of all tribs. */
					ber = getbereich(i,Qges[0][ind]);
					riverflow_table[ind].infilt_amount=0;
					riverflow_table[ind].infilt_fraction=0;
					for (z=0;z<=tracer_count;z++)
					riverflow_table[ind].flow_mm[z]=Qges[z][ind];                      
					riverflow_table[ind].flowm3s=mmh_to_m3s(Qges[0][ind],AE[i],InV);
					riverflow_table[ind].speed  =wqtab[i].v[ber];
				}
			}
			ind=zarray[OL[i]];
			riverflow_table[ind].OL_index = i; 
			riverflow_table[ind].ROUT_index    = nr; 
			riverflow_table[ind].upstream_area = AE[i];  /* area of the first of the SUMTRIBs is already the entire basin area of all tribs. */
			ber = getbereich(i,Qges[0][ind]);
			riverflow_table[ind].infilt_amount=0;
			riverflow_table[ind].infilt_fraction=0;
			for (z=0;z<=tracer_count;z++)
			riverflow_table[ind].flow_mm[z]=Qges[z][ind];                      
			riverflow_table[ind].flowm3s=mmh_to_m3s(Qges[0][ind],AE[i],InV);
			riverflow_table[ind].speed  =wqtab[i].v[ber];
		}
	}
}

void routingclass::riverflow_table_update()
{
	int i,j,z,ber,ind,ind_ol;
	for (i=1;i<=OLs+SUMTRIBs;i++)
	{
		ind_ol=zarray[OL[i]];
		if (indexST[i])
		{
			for (j=1;j<=indexST[i][0];j++)
			{
				ind=zarray[indexST[i][j]];
				ber = getbereich(i,Qges[0][ind]);
				for (z = 0; z <= tracer_count; z++)
				{
					riverflow_table[ind].flow_mm[z] = Qges[z][ind] * timestepfactor;
//					riverflow_table[ind].flow_mm[z] = z == 0 ? Qges[z][ind] * timestepfactor : Qges[z][ind];
//					riverflow_table[ind].flow_mm[z] = Qges[z][ind];
				}
				riverflow_table[ind].infilt_amount   = 0.0;
				riverflow_table[ind].infilt_fraction = 0.0;
				riverflow_table[ind].flowm3s         = mmh_to_m3s(Qges[0][ind],AE[i],InV);
				riverflow_table[ind].speed           = wqtab[i].v[ber];
			}
		}
		ber = getbereich(i,Qges[0][ind_ol]);
		for (z=0;z<=tracer_count;z++)
		{
			riverflow_table[ind_ol].flow_mm[z] = Qges[z][ind_ol] * timestepfactor;
//			riverflow_table[ind_ol].flow_mm[z] = (z == 0) ? Qges[z][ind_ol] * timestepfactor : Qges[z][ind_ol];
//			riverflow_table[ind_ol].flow_mm[z] = Qges[z][ind_ol];
		}
		riverflow_table[ind_ol].infilt_amount   = 0.0;
		riverflow_table[ind_ol].infilt_fraction = 0.0;
		riverflow_table[ind_ol].flowm3s         = mmh_to_m3s(Qges[0][ind_ol],AE[i],InV);
		riverflow_table[ind_ol].speed           = wqtab[i].v[ber];
	}
}

double routingclass::m3s_to_mmh(double p_m3s,double p_TgAe,double p_InV)
{
	/* wandelt den absoluten Abfluss in m^3/s in mm/Zeitintervall um  */
	return p_m3s*(double)0.06*p_InV/p_TgAe;
	/* 0.06 = 3.6 / 60 */
}

double routingclass::mmh_to_m3s(double p_mmh,double p_TgAe,double p_InV)
{
	/* wandelt mm/Zeitintervall in m^3/s um */
	return p_mmh*p_TgAe*(double)16.6666667/p_InV;
}

void routingclass::add_reservoir(std::vector<SP_class*> &reservoirs)
{
	for (int i = 1; i <= SPs; i++)
	{
		reservoirs.push_back(&Speicher[i]);
	}
}

double routingclass::decay(double conc, double decfact)
/* for calculating Ci = Ci-1 * decfact; decfact = exp(-(InV/60/24)/half_time); returns the new concentration after dte decay */
{
	return conc*decfact;
}

double routingclass::mixing(double c1,double c2, double v1,double v2)
/* calc. mixing concentrations for volume1 (c1) with volume 2 (c2); returns mixing concentration */
{
	/* absolute values are necessary in mixing of storages for overtaken routing water, because this storage can get negative values */
	v1=fabs(v1);
	v2=fabs(v2);
	if ((v1+v2)!=(double)0.0)
	return (v1*c1+v2*c2)/(v1+v2); /* return the mixing concentration */
	else
	return (double)0.0;
}

double routingclass::evaporating_changes(double c,double v1,double ve) /* calc. and returning the concentration c after evaporating the amount ve from v */
{
	if ((v1-ve)!=(double)0.0)
	return (v1*c)/(v1-ve);  /* agglomerating e.g. salts due to evaporation / transpiration */
	else
	return (double)0.0;
}

void  routingclass::flushbuffers()
{
	for (int i=1;i<=ALs;i++)
	{
		Abstraction[i].AL_flushbuffer(1);
	}
	for (int i=1;i<=SPs;i++)
	{
		Speicher[i].SP_flushbuffer(1);
	}
}

//void routingclass::AL_flushbuffer(int code,FILE *datei,char *dateiname,double **buffer,unsigned long bufcnt,unsigned long wertcnt)
//{
//unsigned long i,j;
//  if (code>0)
//  {
//    signal(SIGINT,SIG_IGN);
//#if defined _BCC_
//    signal(SIGBREAK,SIG_IGN);
//#endif
//	(*WSstream) << "flushing to disk: "<<dateiname<<" ("<<buffersize<<"d buffersize)"<<"[end-entry]";
//    if ((datei=fopen(dateiname,"a+"))==NULL)
//    {
//      (*WSstream) << "FATAL: error opening abstraction-file "<<dateiname<<"[end-entry]";
//      throw(61);
//    }
//    for (i=1;i<=bufcnt;i++)
//    {
////      fprintf(datei,"%2ld\t%2ld\t%2ld\t%2ld",buffer[i].j,buffer[i].m,buffer[i].t,buffer[i].s);
//      fprintf(datei,"%2.0f\t%2.0f\t%2.0f\t%2.0f",buffer[i][0],buffer[i][1],buffer[i][2],buffer[i][3]);
//      for (j=4;j<=3+wertcnt;j++)
//      {
////        fprintf(datei,"\t%lg",buffer[i].wert[j]);
//        fprintf(datei,"\t%lg",buffer[i][j]);
//      }
//      fprintf(datei,"\n");
//    }
//    fclose(datei);
//    SetSignalAbfangfunktion();
//  }
//}

void routingclass::ReadZuleitung(double **p_Qges,long p_jr,long p_mt,long p_tg,long p_se,long p_me)
{
	/* Wenn Zufluss-modus = extern ist, dann werden aus den vorher geoeffneten Dateien Werte eingelesen */
	int i,idx;
	double wert;
	for (i=1;i<=ZLs;i++)
	{
		if (ZL_art[i]==4)
		{
			idx=ZL[i]+zonmax;
			if (ZL_buffer[i].GetStatus()==FB_STATE_LASTRECORDREAD)
			{
				(*WSstream) << "FATAL: error reading from "<<ZL_buffer[i].GetFilename()<<" (or from its buffer) for inflow "<<ZL[i]
				<< "\nreached end of file!"<<"[end-entry]";
				throw(39);
			}
			ZL_buffer[i].GetActdate(ZLjr,ZLmt,ZLtg,ZLse,ZLme);
			if (ZLjr<1000) ZLjr+=1900;
			if ((ZLjr!=p_jr)||(ZLmt!=p_mt)||(ZLtg!=p_tg)||(ZLse!=p_se) || (ZLme != p_me))
			{
				(*WSstream) << "FATAL: error comparing dates in "<<ZL_buffer[i].GetFilename()<<":\n"
				"     read: "<<ZLjr<<"."<<ZLmt<<"."<<ZLtg<<", "<<ZLse<<":"<<ZLme<<" h,\n"
				"should be: "<<p_jr<<"."<<p_mt<<"."<<p_tg<<", "<<p_se<<":"<<p_me<<" h.\n"<<"[end-entry]";
				throw(36);
			}
			wert=ZL_buffer[i].GetActValue(0);
			if (wert<=-9998) wert=0; // nodata values are handled as zero values for runoff.
			for (int k=1;k<=tracer_count;k++)
			{
				p_Qges[1][idx]=ZL_buffer[i].GetActValue(k);
				if (p_Qges[1][idx]<=-9998)
				{
					p_Qges[1][idx]=0; // nodata values are handled as 0 values
				}
			}
			// Umwandeln der in m^3/s gelesenen Werte in mm/Zeitintervall, es wird die Fläche angenommen,
			// die im Parametesatz als AE=... steht. Die Konzentrationen bleiben, wie sie eingelesen wurden
			p_Qges[0][idx]=m3s_to_mmh(wert,AE[OLs+SUMTRIBs+i],InV);
			ZL_buffer[i].NextRecord();


		}
	}
}

void routingclass::SetSpeicherParams(FILE *steuerdatei,int p_index,std::stringstream &mystream,int p_subint)
{
	int z;

	int p_SPnr=0;
	mystream.seekg(mystream.str().find(_SP_,mystream.tellg())+strlen(_SP_));
	mystream >> p_SPnr;

	if (((OLs+SUMTRIBs+ZLs)==0)&&(total==NULL))
	{
		// reservieren von Speicher fuer Gesamtabfluss im total-array, wenn es in einem Quellgebiet nur einen Speicher, aber kein Routing gibt (OLs==0, ZLs==0, SUMTRIBs==0)
		subint=p_subint;
		if ((total=(double**)calloc(tracer_count+1,sizeof(double*)))==NULL)
		{
			(*WSstream) << "FATAL: memory allocation error in SetSpeicherParams()\n"<<"[end-entry]";
			throw(1);
		}
		for (z=0;z<=tracer_count;z++) 
		total[z]=(double*)calloc(subint+1,sizeof(double));
	}

	p_index=p_index-OLs-SUMTRIBs-ZLs-ALs;
	SP[p_index] = p_SPnr;

	Speicher[p_index].SetSpeicherParams(steuerdatei,p_SPnr,mystream);
	Speicher[p_index].initialize_lake();
	Speicher[p_index].CreateStatFile();
	if (readgrids==0)
	{   // if intialized by grids, no POND-manipulation is required here
		Speicher[p_index].SetPonds(Speicher[p_index].Vakt[0]);
	}
	else
	{    // but the Volume must be taken from the POND-fillage
		// set Vakt[0..z] from the read POND-value
		Speicher[p_index].SetVaktFromPonds();
		// but for initialization of the actWaterLevel, we also re-set the ponds from the actually retrieved reservoir volume
		Speicher[p_index].SetPonds(Speicher[p_index].Vakt[0]);
	}
}

void  routingclass::LocalizeRefALs()
{

	// std::vector<int> vrefAL_localIDs;
	// std::vector<int> vrefAL_globalIDs;
	for (int k = 1; k <= ALs; k++)
	{
		if (Abstraction[k].AL_art == 9)
		{
			for (int i = 0; i < Abstraction[k].vrefAL_globalIDs.size(); i++)
			{
				int globalID = Abstraction[k].vrefAL_globalIDs[i];
				
				if (globalID == 0)
				{
					Abstraction[k].vrefAL_localIDs.push_back(0);
				}
				else
				{
					for (int j = 1; j <= ALs; j++)
					{
						if (Abstraction[j].AL_code == globalID)
						{
							Abstraction[k].vrefAL_localIDs.push_back(j);
							break;
						}
					}
				}
			}
		}
	}

}

void routingclass::ReadReservoirRule(FILE *steuerdatei)
{
	// Lesen der Abgaberegeln fuer die Speicher
	for (int i = 1; i <= SPs; i++)
	{
		Speicher[i].ReadReservoirRule(steuerdatei);
	}
}

void routingclass::FillMapReservoirLevels(std::map<int, double**> *p_mapReservoirLevels)
{
	this->mapReservoirLevels = p_mapReservoirLevels;

	for (int i = 1; i <= SPs; i++)
	{
		Speicher[i].FillMapReservoirLevels(mapReservoirLevels);
	}

}

void  routingclass::FillMapTributariesWithHydCon(std::map<int, bool> *p_mapTributariesWithHydCon)
{
	this->mapTributariesWithHydCon = p_mapTributariesWithHydCon;
	if (SPs > 0)
	{
		Speicher[SPs].FillMapTributariesWithHydCon(mapTributariesWithHydCon); // only the last reservoir may have an hydraulic connection to the next subbasin (the target subbasin of this routing description)
	}
}


void routingclass::SetFirstAbstractionValue()
{
	// Lesen der Abgaberegeln fuer die Speicher
	for (int i = 1; i <= SPs; i++)
	{
		// year must be set to 0 in order to prevent the method to reduce storage volume after calculating the abstraction
		Speicher[i].Abakt = Speicher[i].GetActualOutflow(NULL, NULL, 0, 0, InV, hour_offset,this->jahr, this->monata, this->taga, this->stundea,this->minutea);
		if ((Speicher[i].Abakt*InV*60.0)>Speicher[i].Vakt[0])
		{
			Speicher[i].Abakt = Speicher[i].Vakt[0] / 60.0 / InV;
		}
	}
}

void routingclass::ReadAbstractionRule(FILE *steuerdatei)
{
	// Lesen der Abgaberegeln fuer die Speicher
	for (int i=1;i<=ALs;i++)
	{
		Abstraction[i].ReadAbstractionRule(steuerdatei);
	}
}

routingclass::routingclass(int &retvalue,wasim::mystream *p_Stream,char *steuername,int LogLevel,double p_InV,int p_zonmax,int anz_OL,int anz_SUMTRIB,int anz_et,int anz_ZL,int anz_AL,
	int anz_SP,int tg,double p_TgAe,double p_TgAe_rel,double **p_abfrel,int *p_zarray,
	int p_tracer_count, tracer_info_type *p_tracer_info_table,
	gridinfotype *p_gridinfo,statrecord *p_zone,
	firstlast *p_FirstLastCol ,int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata, int p_lake_active, int p_readgrids,
	char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,long *pv_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, std::vector<HRUMapElement*> *p_vHRU,double hour_offset,int timestep, int timestepfactor,int minutea, int stundea, int taga, int monata, int jahra):
	//memory_handling(),
modelobject(p_Stream,steuername,LogLevel,0,0,0,p_DefaultInputDirectory,p_DefaultOutputDirectory,pv_mpiregioncnt,p_sendHalo,p_recvHalo,p_doSpinUp,p_vHRU),
tracer_count(0),tracer_info(0),WSstream(0),riverflow_table(0),
v1(0),v2(0),weightsum(0),weight(0),frachtsum(0),fracht(0),frachtweightsum(0),frachtweight(0),
discharge(0),dischargesum(0),TgAe(0),TgAe_rel(0),zonmax(0),TG(0),OLs(0),OL(0),SUMTRIBs(0),
indexST(0),ets(0),et(0),SPs(0),SP(0),ZLs(0),ZL(0),ALs(0),AL(0),ZL_art(0),
ZL_buffer(0),
ZL_wert(0),
ZLjr(0),ZLmt(0),ZLtg(0),ZLse(0),flowm3s(0),flowmax(0),
InV(0),kv(0),kh(0),gew(0),total(0),Bh(0),Bv(0),Mh(0),Mv(0),Tm(0),I(0),L(0),AE(0),minspende(0),
maxspende(0),intzahl(0),dq(0),subint(0),vh(0),vv(0),Rh(0),Rv(0),Qh(0),Qh_alt(0),Qv(0),
Qv_alt(0),Th(0),Th_alt(0),Tv(0),Tv_alt(0),qh(0),qv(0),lnq(0),
geb(0),mi(0),ber(0),ver(0),diff(0),Qref(0),CQref(0),gem_Zul(0),Spende(0),Spalte(0),
abfrel(0),zarray(0),Qges(0),Q12ges(0),fluss(0),Qakt(0),zoneSW(0),maxirrig(0),
maxirrigFromReservoirs(0),abstraction_sum(0)
{
	try
	{
		retvalue=0;
		int i,z;
		zoneSW   = NULL; /* initializing pointers with NULL (IMPORTANT for "flowrout.exe") */
		maxirrig = NULL; /* dito ------""------ */
		maxirrigFromReservoirs = NULL; /* dito ------""------ */

		WSstream =  p_Stream;

		InV      =  p_InV;
		this->timestep = timestep;
		this->timestepfactor = timestepfactor;
		zonmax   =  p_zonmax;
		OLs      =  anz_OL;
		SUMTRIBs =  anz_SUMTRIB;
		ets      =  anz_et;
		ZLs      =  anz_ZL;
		ALs      =  anz_AL;
		SPs      =  anz_SP;
		TG       =  tg;

		this->minutea = minutea;
		this->stundea = stundea;
		this->taga = taga;
		this->monata = monata;
		this->jahra = jahra;

		/* EZG des Teilgebietes wird exakt benoetigt, um die Zulaeufe und ABgaben aus dem
	Gebiet von m3/s in mm/Intervall umrechnen zu koennen */
		TgAe   = p_TgAe;

		TgAe_rel = p_TgAe_rel;
		/* fuer einen Sonderfall muss TgAe_rel, der relative Anteil der Einzugsgebietsflaechen
	aller Oberlieger und des Zwischengebietes, die am Routing in einem TG beteiligt sind,
	kleiner 1 sein: wenn naemlich ein Oberlieger, der zwar im Bodenmodell modelliert wurde,
	im Routing nicht geroutet wird, sondern durch eine DAtei ersetzt wird. Da im Bodenmodell
	(bzw. im WaSiM ueberhaupt) die relativen Flaechenanteile der Teilgebiete an der Gesamt-
	flaeche berechnet werden, die sich ja auch aus der nun nicht gerouteten Flaeche zusam-
	mensetzt, gaebe es sonst zu grosse Werte durch die zusaetzliche Einleitung mit Hilfe
	einer Datei. Wenn der OL einfach unberuecksichtigt bleiben soll, dann muss TgAe_rel im
	betrachteten Gebiet aber 1.0 sein, da alle Gewichte der OL und des Zwischengebietes
	eines TG auf die Summe ihrer Gewichte bezogen werden, das TG also automatisch das Gesamt-
	gewicht 1.0 hat. In den folgenden Unterliegenden TG hat diese Gewichtsbildung keinen
	Einfluss mehr (es kann das betrachtete TG also mit seinem normalen Gewicht laut WaSiM
	behandelt werden, da davon ausgegangen wird, dass im Routing des Teilgebietes wirklich
	aller Abfluss erfasst wurde */

		strcpy(this->instanzname,"routingmodel");

		gridinfo     = p_gridinfo;
		zarray       = p_zarray;
		zone         = p_zone;
		FirstLastCol = p_FirstLastCol;
		zonmax       = p_zonmax;
		zonmaxwert   = p_zonmaxwert;
		ncols        = p_ncols;
		nrows        = p_nrows;
		mpistartrow  = p_mpistartrow;
		mpiendrow    = p_mpiendrow;
		xll          = p_xll;
		yll          = p_yll;
		cellsize     = p_cellsize;
		nodata       = p_nodata;
		abfrel		 = p_abfrel;
		lake_active  = p_lake_active;
		readgrids    = p_readgrids;
		tracer_count = p_tracer_count;
		tracer_info  = p_tracer_info_table;



		fracht    = (double*)calloc(tracer_count+1,sizeof(double)); /* for mixing freights in SUMTRIB-channels */
		frachtsum = (double*)calloc(tracer_count+1,sizeof(double)); /* for mixing freights in SUMTRIB-channels */
		//    flowmmh   =  (double*)calloc(tracer_count+1,sizeof(double));
		flowm3s   =  (double*)calloc(tracer_count+1,sizeof(double));
		CQref     =  (double*)calloc(tracer_count+1,sizeof(double));

		total=NULL;
		Qakt=NULL;

		if ((OL     = (int*)  calloc(OLs+SUMTRIBs+ZLs+1,sizeof(int)))==NULL)      {(*WSstream) << "FATAL: routing: memory allocation error 2 "<<"[end-entry]";throw(1);}
		if ((indexST= (int**) calloc(OLs+SUMTRIBs+ZLs+1,sizeof(int*)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 2 "<<"[end-entry]";throw(1);};
		/* die Gewichte muessen fuer alle Oberlieger, auch die am Anfang der Routingstrecke addierten und gemeinsam gerouteten bekannt sein,
	deshalb muss fuer das gew-array die ANzahl Oberlieger + SUMTRIBs + ets + Anzahl Zulaeufe ZLs reserviert werden*/
		if ((gew   = (double*)calloc(OLs+SUMTRIBs+ets+ZLs+1,sizeof(double)))==NULL) {(*WSstream) << "FATAL: routing: memory allocation error 3 "<<"[end-entry]";throw(1);}

		if ((kh    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 0 "<<"[end-entry]";throw(1);}
		if ((kv    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 1 "<<"[end-entry]";throw(1);}
		if ((Bh    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 4 "<<"[end-entry]";throw(1);}
		if ((Bv    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 5 "<<"[end-entry]";throw(1);}
		if ((Mh    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 6 "<<"[end-entry]";throw(1);}
		if ((Mv    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 7 "<<"[end-entry]";throw(1);}
		if ((Tm    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 8 "<<"[end-entry]";throw(1);}
		if ((I     = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 9 "<<"[end-entry]";throw(1);}
		if ((L     = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 10 "<<"[end-entry]";throw(1);}
		if ((AE    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 11 "<<"[end-entry]";throw(1);}
		if ((dq    = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL)     {(*WSstream) << "FATAL: routing: memory allocation error 12 "<<"[end-entry]";throw(1);}
		if ((minspende = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL) {(*WSstream) << "FATAL: routing: memory allocation error 13 "<<"[end-entry]";throw(1);}
		if ((maxspende = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL) {(*WSstream) << "FATAL: routing: memory allocation error 14 "<<"[end-entry]";throw(1);}
		if ((intzahl   = (double*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(double)))==NULL) {(*WSstream) << "FATAL: routing: memory allocation error 15 "<<"[end-entry]";throw(1);}
		/* das reservieren der Spalten fuer diese Tabell erfolgt in
	den einzelnen AUfrufen von SetRoutingParams(), weil erst
	dann die Zeilenzahl bekannt sein wird */
		if ((wqtab = (w_q_tabelle*)calloc(OLs+SUMTRIBs+ZLs+1,sizeof(w_q_tabelle)))==NULL) {(*WSstream) << "FATAL: routing: memory allocation error 16 "<<"[end-entry]";throw(1);}
		/* das erste Gebiet ist immer das Zwischengebiet des Pegels, fuer das das Routing durchgefuehrt wird.
	Dieser Pegel hat keine Translationzeit und auch keine Speicherverzoegerung (die Ganglinie fuer ein
	Teilgebiet wird ja schon im Bodenmodell errechnet). */
		OL[0]  = TG;
		kh[0]  = 0;
		kv[0]  = 0;
		gew[0] = abfrel[0][zarray[TG]];

		/* einrichten der Speicher fuer die Indizierung der Zuleitungen, Ableitungen und Speicher sowie der  */
		if (ZLs>0)
		{
			ZL_buffer   = new FileBuffer[ZLs+1];
			ZL          = (int*)  calloc(ZLs+1,sizeof(int));
			ZL_art      = (int*)  calloc(ZLs+1,sizeof(int));
			ZL_wert     = (double**)calloc(tracer_count+1,sizeof(double*));
			for (z=0;z<=tracer_count;z++) ZL_wert[z]=(double*)calloc(ZLs+1,sizeof(double));
		}
		if (ALs>0)
		{
			AL          = new int[ALs+1];
			Abstraction = new AL_class[ALs+1];
			for (i=1;i<=ALs;i++) 
			{
				Abstraction[i].initialize(WSstream,tracer_info,tracer_count,TG,TgAe,TgAe_rel,InV,DefaultInputDirectory,DefaultOutputDirectory,hour_offset);
			}
			
		}
		if (SPs>0)
		{
			SP=new int[SPs+1];
			//      Speicher = (SP_struct*)calloc(SPs+1,sizeof(SP_struct));
			Speicher = new SP_class[SPs+1];
			for (i=1;i<=SPs;i++) 
			{
				Speicher[i].initialize(WSstream,tracer_info,tracer_count,
				TG,gridinfo,zarray,zone,FirstLastCol,zonmax,zonmaxwert,
				ncols,nrows,mpistartrow,mpiendrow,xll,yll,cellsize,nodata,lake_active,p_DefaultInputDirectory,p_DefaultOutputDirectory,hour_offset,timestep);
				
				// check for the switch for lake model first!
				if (lake_active)
				{
					Speicher[i].LAKE    = gridinfo[getindex(_LAKE_CODES_)].pointer;
					if (ist_vorhanden(_ponding_storage_))
					{
						Speicher[i].POND    = gridinfo[getindex(_ponding_storage_)].pointer;
					}
					else
					{
						Speicher[i].POND    = gridinfo[getindex(_ponding_storage_top_)].pointer;
					}
					Speicher[i].MAXPOND = gridinfo[getindex(_max_ponding_storage_)].pointer;
					Speicher[i].zn      = gridinfo[getindex(_Zonengrid_Bodenmodell_)].pointer;
					Speicher[i].dhm     = gridinfo[getindex(_Hoehenmodell_)].pointer;
					char name[200];
					sprintf(name,"%s%d",_Reale_Verdunstung_,0);
					if (ist_vorhanden(name))
					{
						Speicher[i].ETR = gridinfo[getindex(name)].pointer;
					}
					else
					{
						Speicher[i].ETR = NULL;
					}

					sprintf(name,"%s%d", _Interzeptionsverdunstung_, 0);
					if (ist_vorhanden(name))
					{
						Speicher[i].EI = gridinfo[getindex(name)].pointer;
					}
					else
					{
						Speicher[i].EI = NULL;
					}

					sprintf(name,"%s%s",_soil_bottom_fluxes_,"ext");
					if (ist_vorhanden(name))
					{
						// we don't need the stack since only the exchange between the lake and the uppermost layer is important. This value was copied into the 2D-grid in unsatzonmodel
						Speicher[i].qu      = gridinfo[getindex(name)].pointer; 
					}
					else
					{
						Speicher[i].qu=NULL;
					}
					if (ist_vorhanden(_Basisabfluss_))
					{
						Speicher[i].DELTAGW = gridinfo[getindex(_Basisabfluss_)].pointer;
					}
					else
					{
						Speicher[i].DELTAGW=NULL;
					}
					Speicher[i].PREC    = gridinfo[getindex(_PrecOnLakes_)].pointer;
					Speicher[i].readgrids = readgrids;
					for (int z=1;z<=tracer_count;z++)
					{
						sprintf(name,"concstor_tracer_%d_%s",z,_Unsatzonmodell_);
						Speicher[i].conc[z] = gridinfo[getindex(name)].stackpointer;
					}
					//			next command: must be done later when other values are initialized
					//			Speicher[i].initialize_lake();
				}
			}
		}
		initialized=0; // darf erst nach der letzten Adressuebergabe gesetzt werden
	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	// catching all exceptions in order to allow the object to be constructed correct. This is required 
	// in order to allow a dextructor call in case of an error during construction
	{
		retvalue=103;
	}
}

//void routingclass::SetAblaufModus(int p_index,int p_ALnr,int p_AL_art,const char *p_AL_filename,
//                                  double p_AL_schwelle,double p_AL_anteil,double p_AL_max,const char *p_AL_einheit,int p_subint)

void routingclass::SetAblaufModus(FILE *steuerdatei,int p_index,std::stringstream &mystream,int p_subint,
long jahra,long monata,long taga,long stundea,long minutea)								  
{
	p_index=p_index-OLs-SUMTRIBs-ZLs;

	int p_ALnr=0;
	// Nummer des Ablaufes 

	mystream.seekg(mystream.str().find(_AL_,mystream.tellg())+strlen(_AL_));
	mystream >> p_ALnr;

	AL[p_index] = p_ALnr;

	subint=p_subint;
	if (((OLs+SUMTRIBs+ZLs)==0)&&(total==NULL))
	{
		// reservieren von Speicher fuer Gesamtabfluss im total-array, wenn es in einem Quellgebiet nur eine Ableitung, aber kein Routing gibt
		// (OLs==0, ZLs==0, SUMTRIBs==0, evtl. SPs==0)
		if ((total=(double**)calloc(tracer_count+1,sizeof(double*)))==NULL)
		{
			(*WSstream) << "FATAL: memory allocation error in SetAblaufModus()"<<"[end-entry]";
			throw(1);
		}
		for (int z=0;z<=tracer_count;z++) total[z]=(double*)calloc(subint+1,sizeof(double));
	}
	Abstraction[p_index].SetAbstractionMode(steuerdatei,p_ALnr,mystream,subint,jahra,monata,taga,stundea,minutea);

}

void routingclass::SetZulaufModus(int p_index,int p_ZLnr,int p_ZL_art,const char *p_ZL_filename,
long p_ZL_reihe,int p_ZL_spalte,double p_ZL_wert,
long jahra,long monata,long taga,long stundea,long minutea)
{
	long i;

	//  jahra = jahra%100; /* uebergebenes Anfangsjahr modulo 100 -> Werte 0 .. 99 */
	p_index-=(OLs+SUMTRIBs);
	ZL[p_index] = p_ZLnr;
	ZL_art[p_index] = p_ZL_art;
	if (ZL_art[p_index]==4) /* 4: Werte von File lesen */
	{
		int Header_Rows=ZL_buffer[p_index].initialize(p_ZL_filename,0,100000);

		if (ZL_buffer[p_index].GetStatus()==FB_STATE_FILENOTFOUND)
		{
			(*WSstream) << "FATAL: error opening "<<ZL_buffer[p_index].GetFilename()<<" for inflow "<<p_index<<"[end-entry]";
			throw(37);
		}

		for (i=Header_Rows;i<p_ZL_reihe-1;i++)
		{
			ZL_buffer[p_index].NextRecord();
			if (ZL_buffer[p_index].GetStatus()==FB_STATE_LASTRECORDREAD)
			{
				(*WSstream) << "FATAL: error searching in "<<ZL_buffer[p_index].GetFilename()<<" for row "
				<< p_ZL_reihe<<"\nactual reached row: "<<i<<"[end-entry]";
				throw(38);
			}
		}

		// wind froward to the actual starting date
		ZLjr=ZLmt=ZLtg=ZLse=ZLme=0;
		ZL_buffer[p_index].SetBufPosRel(-1); // auf eine Position vor Beginn setzen
		while (!((ZLjr==jahra)&&(ZLmt==monata)&&(ZLtg==taga) && (ZLse == stundea) && (ZLme == minutea)))
		{
			ZL_buffer[p_index].NextRecord();
			ZL_buffer[p_index].GetActdate(ZLjr,ZLmt,ZLtg,ZLse,ZLme);
			if (ZLjr<1000) ZLjr+=1900;
			if (ZL_buffer[p_index].GetStatus()==FB_STATE_LASTRECORDREAD)
			{
				(*WSstream) << "FATAL: starting date "<<taga<<"."<<monata<<"."<<jahra<<", "<<stundea<<":"<<minutea<<" h in "
				<< ZL_buffer[p_index].GetFilename()<<" not found (" << "search started at row " << p_ZL_reihe << ")" << "[end-entry]";
				throw(38);
			}
		}
		if (ZL_buffer[p_index].GetValcnt()<p_ZL_spalte+tracer_count-4)
		{
			(*WSstream) << "FATAL: there is not the required number of columns in " 
			<< ZL_buffer[p_index].GetFilename() << 
			" for reading ZL["<<p_index<<"] and the follwoing "
			<<tracer_count<<" tracer values starting from column " 
			<<p_ZL_spalte << "[end-entry]";
			throw(38);
		}
		else
		{
			for (int i=0;i<=tracer_count;i++)
			{
				ZL_wert[i][p_index]=ZL_buffer[p_index].GetActValue(i);
			}
		}

	}
	else
	ZL_wert[0][p_index]=p_ZL_wert;
}

/* In dieser Routine werden fuer jede Gerinnestrecke (fuer jede Strecke muss die Routine aufgerufen werden)
die Parameter uebergeben und gespeichert. Aus den Gerinnedaten wird eine w-q-Beziehung berechnet, die dann
in einer Tabell gespeichert wird. Spaeter werden dann beim Routen die in dieser Tabelle festgehaltenen
Translationszeiten und ggf. Speicherkoeffizienten einfach abgerufen bzw. linear interpoliert. */

void routingclass::IterateManning(double *Q,double *B,double *T,double *R,double *v,double M,double I,int maxiter,double maxdiff)
{
	double valt,A,U;
	int k;
	k=0;
	do
	{
		valt =*v;
		A   = *v > 0.0 ? (*Q) / (*v) : 0.0;         /* Querschnitssflaeche ergibt sich aus Abfluss und Geschwindigkeit                          */
		*T   = A/(*B);                               /* Tiefe ergibt sich aus Querschnitssflaeche durch Breite                                   */
		U   = 2.0*(*T)+*B;                          /* benetzter Umfang = 2*Tiefe + Breite (hier liesse sich auch Trapezprofil beruecksichtigen */
		*R   = U> 0.0 ? A/U : 0.0;                   /* hydraulischer Radius aus Flaeche durch benetzten Umfang                                  */
		*v   = *R > 0.0 ? M*pow((double)*R,(double)0.66667)*I : 0.0; /* Berechnung der sich daraus ergebenden Fliessgeschwindigkeit                              */
		k++;
	}
	while ((k<maxiter)&&(fabs(valt-*v) > maxdiff)); /* iterieren, solange sich v um mehr als 0.01 m/s aendert, aber max. 10 mal */
}

void routingclass::SetRoutingParams(int p_index,int p_OL,int *p_numST,double p_kh,double p_kv,double p_Bh,double p_Bv,double p_Tm,double p_Mh,double p_Mv,
	double p_I, double p_L, double p_AE, double p_minspende, double p_maxspende, double p_intzahl, int p_subint, int p_gem_Zul, bool bUseTgAe_rel)
{
	int i,k,z;
	double Qdiffh,Qdiffv;

	OL[p_index]  = p_OL;        /* Nummer des Oberliegers, p_index ist die laufende Nummer des Oberliegers in der aktuellen Routingvorschrift */

	if ((p_index>OLs)&&(p_index<=OLs+SUMTRIBs)) /* wenn die Routingstrecken aus mehreren Oberliegern gespeist wird (SUMTRIB) */
	{
		indexST[p_index]=(int*)calloc(p_numST[0]+1,sizeof(int)); /* reservieren von Speicher fuer die Indizes der gemeinsam gerouteten Oberlieger */
		for (i=0;i<=p_numST[0];i++) /* initialisieren der Indizes; im Element 0 ist die Anzahl enthalten */
		indexST[p_index][i]=p_numST[i];
	}
	else
	{
		indexST[p_index]=NULL; /* wenn der Routingkanal nur aus einem Oberlieger oder aus einem Zulauf gespeist wird, dann soll indexST NULL sein */
	}

	kh[p_index]  = p_kh  > 0 ? p_kh      : 0.01;     /* Speicherkonstante Hauptbett in h, minimal 0.01 h */
	kv[p_index]  = p_kv  > 0 ? p_kv      : 0.01;     /* Speicherkonstante Vorland in h, minimal 0.01 h */
	Bh[p_index]  = p_Bh  > 0 ? p_Bh      : 0.1;      /* Gerinnebreite in m minimal 0.1 m */
	Bv[p_index]  = p_Bv  > 0 ? p_Bv      : 0.1;      /* Vorlandbreite in m minimal 0.1 m */
	Tm[p_index]  = p_Tm  > 0 ? p_Tm      : 0.1;      /* Hauptbettiefe in m minimal 0.1 m */
	Mh[p_index]  = p_Mh  > 0 ? p_Mh      : 1;        /* Gerinnerauhigkeit in (m^(1/3))/s minimal 1 m^(1/3)/s*/
	Mv[p_index]  = p_Mv  > 0 ? p_Mv      : 1;        /* Vorlandrauhigkeit in (m^(1/3))/s minimal 1 m^(1/3)/s*/
	I[p_index]   = p_I  > 0 ? sqrt(p_I): 0.0316;   /* Gerinnegefaelle in m/m minimal 0.001 = 1 m/km*/
	L[p_index]   = p_L  > 0 ? p_L      : 1;        /* Gerinnelaenge in m minimal 1 m*/
	AE[p_index]  = p_AE > 0 ? p_AE     : 0.000001; /* effektives EZG fuer Gerinne in km^2 minimal 1 m^2 */
	minspende[p_index]   = p_minspende;           /* minimale Abflussspende fuer t-q-Beziehung */
	maxspende[p_index]   = p_maxspende;           /* maximale Abflussspende fuer t-q-Beziehung */
	intzahl[p_index]     = p_intzahl;             /* intervallanzahl, in der zwischen minspende und maxspende Tabellenwerte berechnet werden sollen */
	subint      = p_subint;                      /* Anzahl Teilintervalle eines Zeitintervalls (z.B. 12 Teilintervalle a 5 Minuten) */
	gem_Zul     = p_gem_Zul;
	if (p_index==1)
	{ /* nur beim ersten Durchlauf soll der Platz feur diese beiden Arrays belegt werden,
	dass sie ueberhaupt erst hier und nicht schon im Construktor belegt wurden liegt daran, dass
	sie vom Umfang erst durch subint bestimmbar werden, welches an diese function uebergeben wurde */
		total = (double**)calloc(tracer_count+1,sizeof(double*));
		for (z=0;z<=tracer_count;z++)
		total[z]=(double*)calloc(subint+1,sizeof(double));

		Qakt  = (double**)calloc(tracer_count+1,sizeof(double*));
		for (z=0;z<=tracer_count;z++)
		Qakt[z]=(double*)calloc(subint+1,sizeof(double));
	}
	/* nun ist die Intervallanzahl (Abflussintervalle WQ-Tabelle) bekannt, es
	kann fuer die einzelnen wq-Tabellenspalten Speicher reservier werden */
	if ((wqtab[p_index].sp_lskm2 = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1017 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].sp_mm    = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1018 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].Q        = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1019 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].v        = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1020 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].vh       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1021 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].vv       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1022 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].t        = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1023 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].Rh       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1024 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].Rv       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1025 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].Th       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1026 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].Tv       = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1027 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].delta    =   (int*)calloc((int)intzahl[p_index]+2,sizeof(int  )))==NULL) {(*WSstream) << "FATAL: memory allocation error 1028 "<<"[end-entry]";throw(1);}
	if ((wqtab[p_index].QvzuQ    = (double*)calloc((int)intzahl[p_index]+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1029 "<<"[end-entry]";throw(1);}
	/* nun berechnen der Tabelle mit den Translationszeiten */
	/* jeweiliges logarithmisches Increment der Spenden fuer die dann die Fliesszustaende berechnet werden */
	dq[p_index] = log(maxspende[p_index]/minspende[p_index])/intzahl[p_index];
	lnq=log(maxspende[p_index])+dq[p_index];
	vh=1.0; /* erste Geschwindigkeit geschaetzt */
	for (i=(int)intzahl[p_index]+1;i>=1;i--) /* von oben nach unten Fuellen der Tabelle durch Berechnen der Fliessgeschwindigkeit und zugehoeriger DAten*/
	{
		lnq-=dq[p_index];        /* dekrementieren der Spende des logarithmus l/s/km2 */
		wqtab[p_index].sp_lskm2[i] = qh   = exp(lnq);          /* Berechnen der Spende in l/s/km2 */
		wqtab[p_index].sp_mm[i]    = 0.0036 * qh * InV/60;     /* Berechnen der Spende in mm/Zeitintervall */
		if (bUseTgAe_rel)
		{
			wqtab[p_index].Q[i] = Qh = qh / 1000.0*AE[p_index] / TgAe_rel; /* Berechnen des Abflusses in m^3/s */
		}
		else
		{
			wqtab[p_index].Q[i] = Qh = qh / 1000.0*AE[p_index]; /* Berechnen des Abflusses in m^3/s */
		}
		Tv = 0.0;
		vv = 0.0;
		Qv = 0.0; /* Anfangs wurde nicht ueberflutetes Vorland angenommen */
		Rv = 0.0; /* Vorland wird zunaechst als nicht ueberflutet angesehen */

		/* zunaechst Ermitteln von v, T und R fuer Hauptbett */

		IterateManning(&Qh,&(Bh[p_index]),&Th,&Rh,&vh,Mh[p_index],I[p_index],50,(double)0.005);

		/* wenn Th > Tm ist (Ausuferung tritt auf), muss Vorlandabfluss beruecksichtigt werden */
		if (Th>Tm[p_index])
		{
			vv = vh;  /* damit keine Division durch null auftritt beim Iterieren */
			k=0;
			do
			{
				/* Abfluss nach Hauptbett und Vorland aufteilen */
				Th_alt = Th;
				Qh_alt = Qh;
				/* an den Breiten gewichteter Mittelwert aus den Tiefen in Vorland
		und Hauptbett (bezogen auf Hauptbett!) ist neue Hauptbettiefe */
				Th = (Th + Tv + Tm[p_index])/2.0;
				Tv = Th - Tm[p_index]; /* neue Vorlandtiefe ist um dieselbe Differenz gestiegen, wie Hauptbettiefe abnimmt */

				if (Th_alt>Th) /* Vorland bekommt nun mehr Abfluss */
				{
					Qdiffh = (Th_alt-Th)*Bh[p_index]*vh; /* was fuer eine Q-Aenderung bringt es im Hauptbett (gleiche Geschwindigkeit angenommen) */
					Qdiffv = (Th_alt-Th)*Bv[p_index]*vv; /* dito erhoffte Vergroesserung des ABflusses im Vorland */
					/* der Mittelwert aus beiden Q-Aenderungen wird nun abgegeben */
					Qh -= Qdiffh < Qdiffv ? Qdiffh : Qdiffv;
					Qv += (Qh_alt-Qh);
				}
				else /* Hauptbett bekommt mehr */
				{
					Qdiffh = (Th-Th_alt)*Bh[p_index]*vh; /* was fuer eine Q-Aenderung bringt es im Hauptbett (gleiche Geschwindigkeit angenommen) */
					Qdiffv = (Th-Th_alt)*Bv[p_index]*vv; /* dito erhoffte Vergroesserung des ABflusses im Vorland */
					/* der kleinere der beiden Q-Aenderungen wird nun abgegeben */
					Qh += Qdiffh < Qdiffv ? Qdiffh : Qdiffv;
					Qv -= (Qh-Qh_alt);
				}
				/*      Qh = Th/Th_alt * Qh;
		Qv = Qv + (Qh_alt - Qh); */
				/* nun beide Geschwindigkeiten, R's und Tiefen iterieren */
				IterateManning(&Qh,&(Bh[p_index]),&Th,&Rh,&vh,Mh[p_index],I[p_index],50,(double)0.005);
				IterateManning(&Qv,&(Bv[p_index]),&Tv,&Rv,&vv,Mv[p_index],I[p_index],50,(double)0.005);
				k++;
			}
			while ((k<50)&&(fabs(Th-Tv-Tm[p_index]) > 0.01));
		}
		/* Nun fuellen der Tabelle */
		wqtab[p_index].vh[i]    = vh;
		wqtab[p_index].vv[i]    = vv;
		/* wqtab[_index].v ist die wichtige MITTLERE FLIESSGESCHWINDIGKEIT,
	hier aber nicht ueber v=Q/A ermittelt (was dann v = (Qh+Qv)/(Ah+Av) lauten muesste),
	sondern als Abflussanteilsgewichtetes Mittel aus vh und vv. */
		wqtab[p_index].v[i]     = vh*(Qh/(Qh+Qv)) + vv*(Qv/(Qh+Qv));
		//wqtab[p_index].v[i] = (Qh + Qv) / (Bh[p_index] *Th + Bv[p_index]*Tv);
		wqtab[p_index].Rh[i]    = Rh;
		wqtab[p_index].Rv[i]    = Rv;
		wqtab[p_index].Th[i]    = Th;
		wqtab[p_index].Tv[i]    = Tv;
		wqtab[p_index].QvzuQ[i] = Qv/(Qv+Qh);
		/* Fliesszeit in h */
		wqtab[p_index].t[i]     = L[p_index]/wqtab[p_index].v[i]/3600.0;
	}
	/* die richtige Anzahl Unterteilungen der Gesamtstrecke eines Gerinnes, damit sich Wellen nicht ueberholen koennen */
	wqtab[p_index].n=1;
	while ((wqtab[p_index].t[(int)intzahl[p_index]+1]/wqtab[p_index].n+InV/60.0)<(wqtab[p_index].t[1]/wqtab[p_index].n))
	{
		/* das Flowrouting muss dann fuer jede Gerinnestrecke in n Teilen ausgefuehrt werden,
	wobei jedesmal dieselben Parameter gelten. So kann es nicht zum Ueberholen einer
	langsameren durch eine schnelle Welle kommen */
		wqtab[p_index].n++;
	}
	/* je nach Anzahl der Unterteilungen der Strecke muss nun Speicher fuer die Zeiger auf die pro Teilstrecke
	aktuell belegten Abflussarrayplaetze reserviert werden */
	/* speichert nochmal den AUsfluss des letzten Zeitintervalles fuer ELS Berechnung */
	if ((wqtab[p_index].Qaus       = matrix(0,(int)wqtab[p_index].n+1,0,subint+1))==NULL) {(*WSstream) << "FATAL: memory allocation error 1030 "<<"[end-entry]"; throw(1);};
	if ((wqtab[p_index].qvzuq      = (double*)calloc((int)wqtab[p_index].n+2,sizeof(double)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1031 "<<"[end-entry]"; throw(1);};
	if ((wqtab[p_index].last_index =   (int*)calloc((int)wqtab[p_index].n+2,  sizeof(int)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1032 "<<"[end-entry]"; throw(1);};
	if ((wqtab[p_index].ueberrollt = (double**)calloc(tracer_count+1,sizeof(double*)))==NULL) {(*WSstream) << "FATAL: memory allocation error 1033 "<<"[end-entry]"; throw(1);};
	for (z=0;z<=tracer_count;z++)
	wqtab[p_index].ueberrollt[z]=(double*)calloc((int)wqtab[p_index].n+2,sizeof(double));


	kh[p_index] = exp(-1.0/kh[p_index]);
	kv[p_index] = exp(-1.0/kv[p_index]);

	/* fuer jede Gerinnestrecke wird nun fuer jede der intzahl+1 vorgegebenen Spenden die Anzahl von 1/subint-tel Intervallen
	berechnet, die der entsprechende Abfluss beim Routing vom Auslass entfernt abgelegt werden muss */
	int maxdelta=0;
	for (i=1;i<=((int)intzahl[p_index]+1);i++)
	{
		wqtab[p_index].delta[i]=(int)(wqtab[p_index].t[i]/wqtab[p_index].n/InV*60.0*(double)subint);
		if (wqtab[p_index].delta[i]>maxdelta)
		{
			maxdelta=wqtab[p_index].delta[i]; // es kann sein, dass beim Überfluten des Vorlandes die Gesamtfliessgeschwindigkeit kleiner wird
		}
	}

	/* nun wird errechnet, wieviele 1/subint-tel Zeitintervalle fuer jede der maximal n Routingstrecken
	maximal gespeichert werden muessen (also bei einer maximalen Verschiebungszeit von 1.75 Stunden
	und InV = 1h sind es 21 Intervalle a 5 min, bei t = 5h und InV = 24h waeren es 3 Intervalle a 2 h,
	dabei ist die in wqtab[_index].delta[1] gespeicherte Anzahl Intervalle die maximale Zahl pro Teilstrecke,
	fuer die gesamte Strecke muss also der Inhalt von wqtab[].delta[1] mit der Anzahl der Teilstrecken
	multipliziert werden. */
	//  wqtab[p_index].n12t = (int)((wqtab[p_index].delta[1]+(double)subint)*wqtab[p_index].n);
	wqtab[p_index].n12t = (int)((maxdelta+(double)subint)*wqtab[p_index].n);
}

void routingclass::SetWeights()
{
	int i,j;
	double ftotal;
	ftotal=0.0;
	/* Gewichte aus abfrel entnehmen  A C H T U N G: Teilgebiete, die schon aus mehreren OLs bestehen,
	und die in der Routing-Vorschrift vor dem aktuellen TG stehen, sind bereits mit dem Gewicht
	versehen, das der Summe aller seiner Untergebiete am Gesamtgebiet entspricht!  */
	for (j=0;j<=OLs;j++)
	{
		/* die erste ABfrage ist noetig, wenn ein Teilgebiet sich aus sich selbst 
		ergibt (z.B. TG 10 aus OL 10 (...). Dann soll durchGewicht 0 erreicht werden, 
		dass es trotzdem keine Fehler gibt */
		gew[j] = (j>0)&&(zarray[OL[j]]==zarray[TG]) ? 0.0 : abfrel[0][zarray[OL[j]]]; 
		ftotal += gew[j];
	}
	/* gewichte der SUMTRIB-Oberlieger -> Summe der Gewichte ist im Gewicht fuer den ersten in indexST auftretenden Oberlieger gespeichert */
	for (i=1;i<=SUMTRIBs;i++)
	{
		gew[OLs+i]=abfrel[0][zarray[OL[OLs+i]]];
		for (j=1;j<=indexST[i+OLs][0];j++)
		{
			gew[OLs+i]+=abfrel[0][zarray[indexST[i+OLs][j]]];
		}
		ftotal += gew[OLs+i];
	}
	/* Berechnen der Gewichte als Anteile am Teilgebiet TG */
	for (j=0;j<=OLs+SUMTRIBs;j++)
	gew[j] = gew[j]/ftotal*TgAe_rel;
	/* Ersetzen des Gewichtes des Zwischengebietes durch das Gewicht des gesamten
	Teilgebietes am Gesamtgebiet. Somit kann fuer das Routing das Teilgebiet behandelt werden,
	als ob es keine weiteren Untergebiete haette. Die Abflusswerte muessen vorher natuerlich durch
	die fuer dieses Teilgebiet bereits gerouteten Abfluesse ersetz werden. */
	abfrel[0][zarray[TG]]=ftotal;

	/* in abfrel[1][...] ist die jeweilige relative Einzugsgebietsflaeche gespeichert, die an der 
	tatsaechlichen nur Flaeche Modelliert wird. Wird ein Gebiet als Oberlieger genutzt, dessen relative 
	Einzugsgebietsflaeche auch schon kleiner als 1 war, so muss dieser Zufluss wieder mit dem Reziproken 
	des im Oberliegergebiet geltenden relativen Gebietsanteiles multipliziert werden, um fuer das aktuelle
	Gebiet zu einem richtig gewichteten Wert zu gelangen. Dies geschieht jeweils im Routing am Ende, bei der 
	Ueberlagerung der Teilgebietszufluesse. */
	abfrel[1][zarray[TG]]=TgAe_rel;

	if (abfrel[0][0]>0)
	if ((abfrel[0][zarray[TG]]*abfrel[0][0]-0.001>TgAe)||(abfrel[0][zarray[TG]]*abfrel[0][0]+0.001<TgAe))
	{
		(*WSstream) << "INFO: Basin "<<TG<<": area in contr.file: "<<TgAe
		<<" <-> zone-grid "<<abfrel[0][zarray[TG]]*abfrel[0][0]
		<<" (adjusted to the latter)"<<"[end-entry]";

		for (int al=1;al<ALs;al++)
		{
			// convert AL_schwelle from mm with regard to the old are to mm/h with regard to the new area
			// this is exactly the quotient of old and new area
			Abstraction[al].AL_schwelle *= TgAe/(abfrel[0][zarray[TG]]*abfrel[0][0]);
			Abstraction[al].AL_max *= TgAe/(abfrel[0][zarray[TG]]*abfrel[0][0]);
			Abstraction[al].TgAe = abfrel[0][zarray[TG]]*abfrel[0][0];
		}
		TgAe=abfrel[0][zarray[TG]]*abfrel[0][0];
	}
}

double routingclass::maxflow(double *flow)
{
	int i;
	double max;
	max=flow[1];
	for (i=2;i<=subint;i++)
	if (flow[i]>max) max=flow[i];
	return max;
}              

void routingclass::Routing(long jahr,long monat,long tag,long stunde,long minute,int actsubtimestep,int numsubtimesteps)
{
	int i,j,k,geb,ind;

	/* entnehmen der Bewaesserungsmenge aus den Zufluessen, welche Quellfluesse sind 
	(in maxirrig durch 0 gekennzeichnet). Die Bewaesserungsmenge ist schon im Modul irrigation
	entnommen worden, bei Oberliegern auch im Teilgebiet des Oberliegers selbst */
	if ((zoneSW)&&(monat>0))
	{
		for (i=1;i<=OLs+SUMTRIBs;i++)
		{
			if (indexST[i])
			{
				for (k=1;k<=indexST[i][0];k++)
				{
					geb=zarray[indexST[i][k]];
					/* if this is a 1-st level subbasin, irrigation 
			is abstracted now from the inflow (is applied already in 
			the source subbasin of the inflow) */
					if (maxirrig[0][geb]==0.0) 
					{
						for (int z=0;z<=tracer_count;z++)
						maxirrig[z][geb]=Q12ges[z][1][geb];
						/* abziehen der in zoneSW gespeicherten Bewaesserungsmenge -> beide Angaben sind in mm; 
			keine Aenderung der Concentrations noetig */
						for (j=1;j<=subint;j++)
						Q12ges[0][j][geb]-=(zoneSW[geb].sx); /* zoneSW[].sx is already the mean value (in abflussmodellklasse::run() calculated as one of the first commands)) */
					}
				}
			}
			geb=zarray[OL[i]];
			if (maxirrig[0][geb]==0.0) /* if this is a 1-st level subbasin, irrigation is abstracted now from the inflow (is applied already in the source subbasin of the inflow) */
			{
				for (int z=0;z<=tracer_count;z++)
				maxirrig[z][geb]=Q12ges[z][1][geb];
				/* abziehen der in zoneSW gespeicherten Bewaesserungsmenge -> beide Angaben sind in mm; 
		keine Aenderung der Concentrations noetig */
				for (j=1;j<=subint;j++)
				Q12ges[0][j][geb]-=(zoneSW[geb].sx); /* zoneSW[].sx is already the mean value (in abflussmodellklasse::run() calculated as one of the first commands)) */
			}
		}
	}

	/* der Gesamtausfluss aus dem Gebiet einschliesslich der gerouteten Oberlieger ergibt sich dann als Summe
	aus diesen Oberliegerabfluessen und dem am relativen Gebietsanteil des Zwischengebietes
	gewichteten Abfluss des Zwischengebietes */
	for (i=1;i<=subint;i++)
	{
		total[0][i]=gew[0]*Q12ges[0][i][zarray[TG]];
		for (int z=1;z<=tracer_count;z++)
		total[z][i]=Q12ges[z][i][zarray[TG]]; /* concentrations are not weighed */
	}   

	/* nun wird jede Routingstrecke einzeln abgearbeitet.
	Zulaeufe werden ebenfalls geroutet, haben aber ein Gewicht von 0.0 (gew[OLs+1]...gew[OLs+ZLs])
	und muessen nachher extra zum Gesamtabfluss addiert werden */
	for (i=1;i<=OLs+SUMTRIBs+ZLs;i++)
	{
		geb = i<=OLs+SUMTRIBs ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */

		if ((monat > 0) && (i <= OLs + SUMTRIBs) && (riverflow_table))
		{
			if (indexST[i])
			{
				for (k=1;k<=indexST[i][0];k++)
				{
					ind= zarray[indexST[i][k]];
					/* converting infiltration amount from m3/s into mm/intervall */
					if (actsubtimestep == 0)
					{
						riverflow_table[ind].infilt_amount = m3s_to_mmh(riverflow_table[ind].infilt_amount, riverflow_table[ind].upstream_area, InV);
					}
/* now reducing inflow by the infiltration amount -> negative numbers may occure due to transition of one time intervall, if inflow is less than infiltration amount */
					for (j = 1; j <= subint; j++)
					{
						Q12ges[0][j][ind] -= riverflow_table[ind].infilt_amount;
					}
				}
			}
			/* converting infiltration amount from m3/s into mm/intervall */
			if (actsubtimestep == 0)
			{
				riverflow_table[geb].infilt_amount = m3s_to_mmh(riverflow_table[geb].infilt_amount, riverflow_table[geb].upstream_area, InV);
			}
			/* now reducing inflow by the infiltration amount -> negative numbers may occure due to transition of one time intervall, if inflow is less than infiltration amount */
			for (j = 1; j <= subint; j++)
			{
				Q12ges[0][j][geb] -= riverflow_table[geb].infilt_amount;
			}
		}

		/* if the actual routing channel is a SUMTRIB-channel, the discharge of the tributaries has to be added as well as the concentrations have to be mixed */
		if ((indexST[i])&&(i<=OLs+SUMTRIBs)) /* indexST!=NULL, if this is a SUMTRIB-channel */
		//    if ((indexST[i])&&(i<OLs+SUMTRIBs)) /* indexST!=NULL, if this is a SUMTRIB-channel */
		/* !!!! this block is run instead of the entire routing for all SUMTRIB basins except the first one!!!!  */
		{
			for (j=1;j<=subint;j++)
			{
				weightsum=frachtweightsum=dischargesum=0;
				for (int z=1;z<=tracer_count;z++) frachtsum[z]=0;
				for (k=1;k<=indexST[i][0];k++)
				{
					ind              = zarray[indexST[i][k]];
					weight = abfrel[0][ind] / abfrel[1][ind]; 
					weightsum       += weight;

					
					if ((gem_Zul >= 1) && (ind <= zonmax) && (Spalte[ind] != 0) && (Spende != NULL) && (Spende[10][Spalte[ind]] > nodata))
					{
						if (gem_Zul == 2)
						{
							discharge = m3s_to_mmh(Spende[10][Spalte[ind]], abfrel[0][0]*weight, InV);
						}
						else
						{
							discharge = Spende[10][Spalte[ind]];
						}
					}
					else
					{
						discharge = Q12ges[0][j][ind];
					}
					dischargesum    += discharge*weight;
					frachtweight     = weight*discharge;
					frachtweightsum += frachtweight;
					for (int z=1;z<=tracer_count;z++)
					{
						fracht[z]        = frachtweight*Q12ges[z][j][ind];
						frachtsum[z]    += fracht[z];
					}
				}
				/* the first one of the tribs. is a normal subbasin */
				weight           = abfrel[0][geb]/abfrel[1][geb];
				//        weight           = abfrel[0][geb]-weightsum; weiss auch nicht, wieso das so aussah
				weightsum       += weight;

				if ((gem_Zul >= 1) && (geb <= zonmax) && (Spalte[geb] != 0) && (Spende != NULL) && (Spende[10][Spalte[geb]] > nodata))
				{
					
					if (gem_Zul == 2)
					{
						discharge = m3s_to_mmh(Spende[10][Spalte[geb]], abfrel[0][0] * weight, InV);
					}
					else
					{
						discharge = Spende[10][Spalte[geb]];
					}
				}
				else
				{
					discharge = Q12ges[0][j][geb];
				}
//				discharge = Q12ges[0][j][geb];
				dischargesum    += discharge*weight;
				frachtweight     = weight*discharge;
				frachtweightsum += frachtweight;
				for (int z=1;z<=tracer_count;z++)
				{
					fracht[z]        = frachtweight*Q12ges[z][j][geb];
					frachtsum[z]    += fracht[z];
				}

				dischargesum  = dischargesum/weightsum; /* weighted average of inflows -> specific discharge */
				Q12ges[0][j][geb] = dischargesum;
				for (int z=1;z<=tracer_count;z++)
				Q12ges[z][j][geb] = frachtsum[z]/frachtweightsum;
			}
		}

		if (i>OLs+SUMTRIBs)
		if (ZL_art[i-OLs-SUMTRIBs]==3)   /* fuer Zulaeufe aus internen Ableitungen muessen nun die m3/s in mm/Zeitiuntervall umgewandelt werden */
		{
			for (j=1;j<=subint;j++)
			{
				/* konverting discharge from m3/s to mm/timestep; but: tracer concentrations remain constant! */
				Q12ges[0][j][geb]=m3s_to_mmh(Q12ges[0][j][geb],AE[i],InV);
			}
			for (int z=0;z<=tracer_count;z++) 
			Qges[z][geb]=Q12ges[z][1][geb]; /* getting the inflow for the next timestep (if possibly this inflow is get from a routing element laying behind the actual element -> reversala of flow!) */
		}


		/* Anzahl der maximalen Intervalle einer Teilstrecke errechnen, wird oft benoetigt und spart so Zeit */
		mi = wqtab[i].delta[1]+subint;

		/* fuer die oberste der Teilstrecken sind die Abfluesse aus der Eingabedatei bzw dem Bodenmodell in Qges gegeben,
	oder, wenn gemessene Zulaufganglinien verwendet werden, sie werden aus der Spendendatei geliefert
	was vor Aufruf der Routing-Routinen in Q12ges kopiert wurde (jeweils <subint> gleiche Werte fuer ein Intervall.
	Diese Werte werden nun in Qakt uebernommen, welches fuer jede unterliegende Teilstrecke dann jeweils
	aus den gerouteten Werten der oberen Teilstrecken als Eingabeabfluss dient (Zufluss) */

		if ((gem_Zul >= 1) && (geb <= zonmax) && !((indexST[i]) && (i <= OLs + SUMTRIBs)) && (Spalte[geb] != 0) && (Spende != NULL) && (Spende[10][Spalte[geb]] > nodata))
		{
			if (gem_Zul == 2)
			{
				for (k = 1; k <= subint; k++)
				{
					Qakt[0][k] = m3s_to_mmh( Spende[10][Spalte[geb]],AE[i],InV); /* nun gehen gemessene Ganglinien in die Translation */
				}
			}
			else
			{
				for (k = 1; k <= subint; k++)
				{
					Qakt[0][k] = Spende[10][Spalte[geb]]; /* nun gehen gemessene Ganglinien in die Translation */
				}
			}
		}
		else
		{
			/* wenn geb > zonmax ist, dann ist also ein Zulauf mit der Berechnung dran,
		welcher NICHT aus der Spendendatei gelesen werden darf; es kann ein externer oder ein interner Zulauf sein,
		in beiden Faellen wird davon ausgegangen dass in Q12ges auch die Tracerkonzentrationen enthalten sind */
			for (k = 1; k <= subint; k++)
			{
				for (int z = 0; z <= tracer_count; z++)
				{
					Qakt[z][k] = Q12ges[z][k][geb];
				}
			}
		}


		std::map<int, bool>::iterator p;
		p = (mapTributariesWithHydCon->find(geb));
		bool useHydCon = false;
		if (p != mapTributariesWithHydCon->end()) useHydCon = p->second;

		if (useHydCon)
		{
//			int sourceSubBasin = p->first;
			for (int z = 0; z <= tracer_count; z++)
			{
				for (int k = 1; k <= subint; k++)
				{
					fluss[z][geb][k] = Q12ges[z][k][geb];
				}
			}
		}
		else
		{
			/* innerhalb der einzelnen Strecken gibt es je j Teilstrecken,
			dabei wird von oben nach unten geroutet */
			for (j = (int)wqtab[i].n; j >= 1; j--)
			{
				/* welcher Abflussbereich ist fuer  den maximalen Abfluss in diesem Intervall gueltig?
			(Index in w-q-tabelle wqtab, daraus berechnet sich die Verschiebung fuer den gesamten
			Zufluss der <subint> Teilintervalle) */
				ber = getbereich(i, maxflow(Qakt[0]));

				/* wie gross ist die aktuelle Verschiebung? */
				ver = wqtab[i].delta[ber];
				/* calculating the actual decay factor for the time, the water is running down the channel */
				for (int z = 1; z <= tracer_count; z++)
				{
					if (tracer_info[z].half_time <= 0)
						tracer_info[z].routdecfact = 0;
					else
						tracer_info[z].routdecfact = exp(-(ver*InV / subint / 60 / 24) / (tracer_info[z].half_time*_1_ln2_));
				}

				/* Verschiebung auf Index im fluss-array umrechnen */
				ver = subint + 1 + (j - 1)*mi + ver;

				/* hier muss nun ab dem Index "ver" die Werte von Qakt rein,
			wobei je nach letztem last_index entweder ueberolltes Wasser
			zwischengespeichert werden muss (Lastindex >= ver) oder,
			bei Luecken zwischen Lastindex und ver, diese Luecken gefuellt
			werden muessen */
				if (ver > (wqtab[i].last_index[j] + 1))
				{
					/* Auffuellen von Luecken aus dem Speicher */

					if (fluss[0][geb][wqtab[i].last_index[j]] >= Qakt[0][1])
					{
						/* der Normalfall: Luecken entstehen beim Ablaufen eines Hochwassers,
				also ist der aktuelle Abfluss kleiner als der vorhergehende */
						Qref = Qakt[0][1];
						for (int z = 1; z <= tracer_count; z++) /* for (z=1 ist richtig, weil CQref[0] nicht gebraucht wird! */
							CQref[z] = Qakt[z][1];
					}
					else
					{
						/* beim Austreten des Wassers ins Vorland kann aber die Geschwindigkeit nach dem
				in der Translationszeitenberechnung genutzten, am jeweiligen Abfluss
				in Hauptbett und Vorlaendern gewichteten Fliessgeschwindigkeitsmittel
				ploetzlich groesser werden, also muss dann der Abfluss des
				vorherigen Intervalles, der ja noch im letzten (mit last_index bezeichneten
				Arrayelement liegt, aufgefuellt werden.
				*/
						k = wqtab[i].last_index[j];
						Qref = fluss[0][geb][k];
						for (int z = 1; z <= tracer_count; z++) /* for (z=1 ist richtig, weil CQref[0] nicht gebraucht wird! */
							CQref[z] = fluss[z][geb][k];
					}

					diff = (ver - wqtab[i].last_index[j] - 1)*Qref;

					/* in dieser Programmversion wird zunaechst erstmal davon ausgegangen,
			dass ueberollt auch negativ werden kann. Bei spaeterem Ueberollen neuer
			Wellen wird es dann wieder aufgefuellt. So bleibt die Bilanzreinheit auf lange
			Sicht gewaehrt */

					wqtab[i].ueberrollt[0][j] -= diff; /* beim Entnehmen von Wasser aus dem Puffer muessen die Konzentrationen nicht neu berechnet werden */

					for (k = wqtab[i].last_index[j] + 1; k < ver; k++)
					{
						fluss[0][geb][k] = Qref;
						for (int z = 1; z <= tracer_count; z++)
							fluss[z][geb][k] = CQref[z];
					}
				}
				for (k = ver; k <= ver + subint - 1; k++)
				{
					if (k <= wqtab[i].last_index[j])
					{
						/* speichern ueberrollten Wassers, auch dabei ist es normal, dass ein ansteiegendes
				Hochwasser kleinere ABfluesse ueberrollt, also werden die kleineren ABfluesse des vorherigen
				Intervalles in ueberrollt gespeichert. Beim Rueckzug des Wassers ins Hauptbett
				wird v nach dem gewaehlten Verfahren der am relativen Anteil am Gesamtabfluss
				gewichteten Mittelung der Geschwindigkeiten im Hauptbett und auf dem Vorland dagegen groesser,
				was dazu fuehrt, das der aktuelle Abfluss gespeichert werden muss, und der alte Abfluss
				ungestoert seinen Weg fortsetzen kann.
				*/
						if (fluss[0][geb][k] <= Qakt[0][k - ver + 1])
						{
							/* mixing of tracers in buffer <ueberrollt> */
							v1 = wqtab[i].ueberrollt[0][j];
							v2 = fluss[0][geb][k];
							for (int z = 1; z <= tracer_count; z++)
								wqtab[i].ueberrollt[z][j] = (double)mixing(wqtab[i].ueberrollt[z][j], fluss[z][geb][k], v1, v2);
							/* now adding the exceeding amount of water to the buffer amount */
							wqtab[i].ueberrollt[0][j] += v2;
						}
						else
						{
							/* mixing of tracers in buffer <ueberrollt> */
							v1 = wqtab[i].ueberrollt[0][j];
							v2 = Qakt[0][k - ver + 1];
							for (int z = 1; z <= tracer_count; z++)
								wqtab[i].ueberrollt[z][j] = (double)mixing(wqtab[i].ueberrollt[z][j], Qakt[z][k - ver + 1], v1, v2);
							/* now adding the exceeding amount of water to the buffer amount */
							wqtab[i].ueberrollt[0][j] += Qakt[0][k - ver + 1];
						}
					}
					/* nun noch einfaches kopieren der Qakt-Abfluesse in die Fluss-array-Felder */
					for (int z = 0; z <= tracer_count; z++)
						fluss[z][geb][k] = Qakt[z][k - ver + 1];
				}

				/* Werte verschieben, das eigentliche Routing, aber nun um konstant <subint> Werte */
				for (k = 1 + (j - 1)*mi; k <= j*mi; k++)
				{
					for (int z = 0; z <= tracer_count; z++)
						fluss[z][geb][k] = fluss[z][geb][k + subint];
					/* radioactive decay with a decay factor correspondent to the real flow time (not <subint>-Intervals but <ver>-Intervals */
					for (int z = 1; z <= tracer_count; z++)
						fluss[z][geb][k] *= tracer_info[z].routdecfact;
				}
				/* neuen last_index berechnen */
				wqtab[i].last_index[j] = ver - 1;

				/* nach dem Routen werden die ersten <subint> Abfluesse aus der ersten Strecke
			die letzten <subint> in der zweiten Teilstrecke sein (muessen aber eventuell
			nach vorn geschoben werden, weshalb sie nun gleichzeitig Qakt[1..subint] in
			dieser Teilstrecke werden (auuserdem wird nun die Retention berechnet) */

				/* zunaechst Feststellen des Abflussbereichs des Ausflusses im aktuellen
			Intervall in der wq-Tabelle der aktuellen Strecke i */
//				ber = getbereich(i, maxflow(&(fluss[0][geb][1 + (j - 1)*subint])));
				ber = getbereich(i, maxflow(&(fluss[0][geb][1 + (j - 1)*mi]))); // new 20210505: corrected index for ber calculation

				for (k = 1 + (j - 1)*mi; k <= subint + (j - 1)*mi; k++)
				{
					/* nun uebernehmen der verschobenen Abfluesse in Qakt
			(als Eingabe in naechstes Teilelement)
			sowie in wqtab[i].Qaus, weil diese Abfluesse noch fuer den
			Auslaufterm der Retentionsberechnung im folgenden Intervall
			benoetigt werden */
					Qh = fluss[0][geb][k] * (1.0 - wqtab[i].QvzuQ[ber]);                /* aktueller Hauptbettabfluss */
					Qv = fluss[0][geb][k] * wqtab[i].QvzuQ[ber];                 /* aktueller Vorlandabfluss   */
					Qh_alt = wqtab[i].Qaus[j][k - (j - 1)*mi] * (1.0 - wqtab[i].qvzuq[j]); /* letzter Hauptbettabfluss   */
					Qv_alt = wqtab[i].Qaus[j][k - (j - 1)*mi] * wqtab[i].qvzuq[j];  /* letzter Vorlandabfluss     */

					Qakt[0][k - (j - 1)*mi] = Qh * (1.0 - kh[i]) + Qh_alt*kh[i] + Qv * (1.0 - kv[i]) + Qv_alt*kv[i];

					/* Qaus = Qakt, Qakt wird pro Teilstrecke veraendert,
			Qaus wird im naechsten Intervall aber noch benoetigt */
					wqtab[i].Qaus[j][k - (j - 1)*mi] = Qakt[0][k - (j - 1)*mi];

				}

				/* Verhaeltnis des Vorland- zum Gesamtabfluss fuer aktuellen Ausfluss
			aus der Teilstrecke j im i-ten Routingelement wird im sonst nicht
			benutzten 0-ten Element des aktuellen Ausflusses gespeichert, um
			im naechsten Intervall die Retentionberechnung fuer die Vorland- und
			Hauptbettkomponenten jeweils richtig berechnen zu koennen (was den
			Leerlaufterm angeht). "ber" wurde bereits fuer den (unretentionierten)
			aber gerouteten Ausfluss bestimmt. Wenn nur noch weniger als 0.1 %
			des Gesamtabflusses vom Vorland kommen, wird QvzuQ auf 0 gesetzt, dann
			geht es anschliessend schneller */

				wqtab[i].qvzuq[j] = wqtab[i].QvzuQ[ber];
			}
		}

		if (i<=OLs+SUMTRIBs)
		{ /* als letzte Aktion in jeder Routingstrecke wird der geroutete Abfluss zum Gesamtabfluss addiert */
			/* bei der letzten Teilstrecke ist Qakt[k]==fluss[zarray[OL[i]]][k] */
			for (k=1;k<=subint;k++)
			{
				/* mixing tracer concentrations in total-array */
				v1 = total[0][k];
				if (indexST[i])
				{
					// This branch is reached only if SUMTRIBS ar routed. In such a case, the variable weightsum contains the
					// sum of the single weights of all participants of the sumtrib-cannel on the modelled basin normalized by
					// their respective weight on total real catchment area of the single participants.
					// example: if a routing channel is a sumtrib-channel for subbasins 4, 14 and 15, weightsum is calculated
					// (at the beginning of routing) by
					// weightsum=abfrel[0][4]/abfrel[1][4]+abfrel[0][14]/abfrel[1][14]+abfrel[0][15]/abfrel[1][15];
					// with abfrel[0][..] the weight of the subbasin on the total modelled area
					// abfrel[1][..] the weight of the modelled part of a subbasin on the real catchment area of this basin
					// (it is 1.0 for most cases, except if a river basin is not modelled as a whole but only in its lower parts
					// so that a influx has to be considered and the relative fraction of all affected subbasins has to be set to
					// values less than 1.0.
					// The expression is ont simple explainable, but it is comparable to the analogous expression
					// for conventional routing channels: the weighting for Qakt[0][k] is the sum of the single weightings of all
					// participants of the SUMTRIB.....
					v2 = Qakt[0][k] * TgAe_rel / abfrel[0][zarray[TG]] * weightsum;
				}
				else
				{
					v2 = Qakt[0][k] * (gew[i] / abfrel[1][zarray[OL[i]]]);
				}
				for (int z = 1; z <= tracer_count; z++)
				{
					total[z][k] = (double)mixing(total[z][k], Qakt[z][k], v1, v2);
				}
				/* adding the water amount v2 to v1 */
				total[0][k] += v2;
			}
		}
		else
		{ 
			// AE[i]/TgAe ist die relative Groesse des Zuflussgebietes am Gebiet des
			// betrachteten TG's, wobei darin die Groesse des Zuflussgebietes nicht enthalten sein darf.
			// Wenn AE[i] z.B. 0.5*TgAe ist, dann wird der aus dem Zulauf stammende
			// Abfluss in Bezug auf das gesamte TG nur noch zu 50% gewichtet. Beim Einlesen derZufluesse wurde
			// der Zufluss von m3/s in mm/Zeitintervall umgewandelt, dabei wird Q durch AE geteilt. Hier wird daswieder
			// rueckgangig gemachtund aus der Spende bezogen auf das fiktive Einzugsgebiet des Zuflusses wird eine
			// Spende bezogen auf das modell-Teilgebiet (in m3s_to_mmh:  mm/h = _m3s*3600000/60*_InV/AE/1000000;
			// durch mm/h/AE*TgAe wird daraus: mm/h = m3/s*3600000/60*InV/TgAe/1000000, also eine Spende, bezogen auf TgAe
			for (k=1;k<=subint;k++)
			{
				v1 = total[0][k];
				v2 = Qakt[0][k]*AE[i]/TgAe*TgAe_rel;   /* TgAe_rel is the relative fraction of the modelled catchment area to the total upstream catchment area. So specific discharges are referring to the total area, like measurements at strem gauges do, too */
				for (int z=1;z<=tracer_count;z++)
				total[z][k] = (double)mixing(total[z][k],Qakt[z][k],v1,v2);
				total[0][k] += v2;
			}
		}
	}


	/* Das Teilgebietsergebnis wird in des Q12ges-array geschrieben, von dort wird spaeter das
	naechste Teilgebiet geroutet und es werden die Ergebnisse geschrieben */
	for (int z=0;z<=tracer_count;z++) 
	{
		total[z][0]=0.0; /* fuer mittleren Abfluss und mittlere Konzentrationen */
	}

	for (i=1;i<=subint;i++)
	{
		for (int z=0;z<=tracer_count;z++)  
		{
			Q12ges[z][i][zarray[TG]]=total[z][i];
		}
		total[0][0]+=total[0][i]; /* nur Abflussmenge wird zunaechst addiert, Konzentrationen werden anschliessend gewichtet gemittelt */
	}
	/* nun steht in total[0][0] die Summe der Abfluesse und in total[0][i] der jeweilige Intervallabfluss */
	for (i=1;i<=subint;i++)
	{
		/* nun gewichtete Mittelung der Konzentrationen der tracer */
		for (int z=1;z<=tracer_count;z++)
		total[z][0]+=total[0][0]>0?(total[z][i]*(total[0][i]/total[0][0])):0.0; /* gewichtet am Anteil des Abflusses im Teil-Intervall gemessen am Gesamtintervall */
	}

	/* nun noch den Mittelwert fuer die Ausgabe in die Statistikdatei errechnen */
	Qges[0][zarray[TG]] = total[0][0]/(double)subint;
	for (int z=1;z<=tracer_count;z++) Qges[z][zarray[TG]] = total[z][0];

	/* nun festlegen der maximalen Entnahmemenge im naechsten Zeitintervall (in mm/zone)
	sowie Entnahme der im Bewaesserungsmodul berechneten Menge (in zoneSW gespeichert in mm/Zone)
*/
	if ((zoneSW)&&(monat>0))
	{
		if (maxirrig[0][zarray[TG]]==0.0) /* maxirrig is zero, if it wasnt calculated for the actual basin yet */
		{
			maxirrig[0][zarray[TG]]=Qges[0][zarray[TG]]*(abfrel[0][zarray[TG]]/zoneSW[zarray[TG]].rel)/TgAe_rel; /* new max. value for irrigation */
			for (int z=1;z<=tracer_count;z++)
			maxirrig[z][zarray[TG]]=Qges[z][zarray[TG]]; /* concentrations for irrigation water is the same like the mixed concentration of all inflows */
			for (j=1;j<=subint;j++)
			{
				Q12ges[0][j][zarray[TG]]-=(zoneSW[zarray[TG]].sx)*(zoneSW[zarray[TG]].rel/abfrel[0][zarray[TG]])*TgAe_rel; /* zoneSW[].sx is already the mean value (in abflussmodellklasse::run() calculated as one of the first commands)) */
				if (Q12ges[0][j][zarray[TG]]<0)
				{
					Q12ges[0][j][zarray[TG]]=0;
				}
			}
			/* reducing also the total amount that counts for balances with other subasins and for output (Q12ges is only an internal variable) */
			Qges[0][zarray[TG]]-=(zoneSW[zarray[TG]].sx)*(zoneSW[zarray[TG]].rel/abfrel[0][zarray[TG]])*TgAe_rel;
			if (Qges[0][zarray[TG]]<0)
			{
				Qges[0][zarray[TG]]=0;
			}
		}
	}

	abstraction_sum=0;
	/* nun koennen aus dem Speicher oder aus dem Abfluss in Qges noch die Abgaben in Fremdgebiete
	abgezogen werden, wenn ALs > 0 ist */
	if (ALs > 0)
	{ /* in die entsprechenden Qges-Elemente, die sich hinter den <zonmax> Elementen fuer echte Teilgebiete
	aufsteigend entsprechend der Nummer der ALs einordnen */
		double irr=0;
		// Bewässerungsentnahme aus dem letzten Zeitschritt im irrigation-model muss noch entnommen werden können, dieser Rest muss also im Fluss bzw. im Speicher bleiben
		if ((zoneSW)&&(monat>0))
		{
			irr = InV*60.0*mmh_to_m3s(maxirrigFromReservoirs[0][zarray[TG]],TgAe/TgAe_rel,InV); // Umrechnen mm/Interval in m^3
		}
		flowmax = SPs>0 ? m3s_to_mmh((Speicher[SPs].Vakt[0]-irr)/InV/60.0,TgAe/TgAe_rel,InV) : Qges[0][zarray[TG]]; /* maximaler Abfluss in einem Intervall in m3/s -> dann laeuft Speicher aber leer! */
		for (i=1;i<=ALs;i++) /* first: outflows (abstractios), second: satisfy irrigation water demand */
		{
			if (Abstraction[i].AL_art!=8) // abstractions from reservoir with AL_art==8 will be handled in SP_class 
			{
				double flowmmh = Abstraction[i].CalcWaterBalance(Qges,abstraction_sum,flowmax,jahr,monat,tag,stunde,minute);
				flowm3s[0] = mmh_to_m3s(flowmmh,TgAe/TgAe_rel,InV);
				for (int z=1;z<=tracer_count;z++) 
				{
					if (SPs==0)
					{
						flowm3s[z]=Qges[z][zarray[TG]]; /* uebernehmen der Konzentrationen fuer das Schreiben der Abgaben */
					}
					else
					{
						flowm3s[z]=Speicher[SPs].Vakt[z];
					}
				}
				if (monat>0)  /* soll in der Initialisierung nicht geschrieben werden */
				{
					Abstraction[i].WriteAbstraction(jahr,monat,tag,stunde,minute,flowm3s);
				}
				if (SPs==0) // otherwise, AL amounts will be taken from SP-elements later (using abstraction_sum)
				{
					/* wenn kein Speicher da ist, dann muss die Ableitung nur vom aktuellen Abfluss abgezogen werden
			(ist nie groesser als Qges, da flowmax maximal Qges ist) */
					/* Konzentrationen BLEIBEN GLEICH!! */
					Qges[0][zarray[TG]]-=flowmmh;
					for (j=1;j<=subint;j++)
					{
						Q12ges[0][j][zarray[TG]]-=flowmmh;
					}
				}
				/* bei internen Ableitungen muss die Ableitung nun noch im Qges-Array gespeichert werden als Zulauf derselben Nummer */
				if (     (Abstraction[i].AL_art==3) // _intern_
					  || (Abstraction[i].AL_art==6) // _intern_with_rule
				   ) 
				{
					for (j=1;j<=subint;j++)
					{
						Q12ges[0][j][zonmax+AL[i]]=flowm3s[0];
						for (int z=1;z<=tracer_count;z++) /* Schleife ab z=1 nur fuer die Konzentrationen!! */
						Q12ges[z][j][zonmax+AL[i]]=Qges[z][zarray[TG]];
					}
					Qges[0][zonmax+AL[i]]=flowm3s[0];
				}
			}
		}
		
	}
	// der in Qges gespeicherte Abfluss kann nun in den Speicher geleitet werden,
	// um die Abgaberegel anzuwenden. 
	if (SPs > 0)
	{
		for (i=1;i<=SPs;i++)
		{
			// zunächst entnehmen der Bewässerungsmenge aus den Speichern, aber nur aus dem letzten Speicher einer Serie
			double irr=0;
			if ((zoneSW)&&(monat>0)&&(i==SPs))
			{
				irr = InV*60.0*mmh_to_m3s(maxirrigFromReservoirs[0][zarray[TG]],TgAe/TgAe_rel,InV); // Umrechnen mm/Interval in m^3
			}
			flowm3s[0] = mmh_to_m3s(Qges[0][zarray[TG]],TgAe/TgAe_rel,InV); // Umrechnen mm/h in m^3/s
			double abstraction_summ3;
			if ((monat>0)&&(i==SPs)) // take abstractions from the last reservoir only and not in initialization routing (when monat==0)
			{
				abstraction_summ3 = InV*60.0*mmh_to_m3s(abstraction_sum,TgAe/TgAe_rel,InV); // Umrechnen mm/h in m^3/s
			}
			else
			{
				abstraction_summ3=0.0; 
			}
			Speicher[i].CalcWaterBalance(InV,irr,abstraction_summ3,fluss,Qges,flowm3s,Abstraction,ALs,jahr,monat,tag,stunde,minute,timestepfactor,actsubtimestep,numsubtimesteps);

			for (int al=1;al<=ALs;al++)
			{
				// here we do not need to add AL_art == 9 abstractions (because, if read in abstraction was already used to copy it's value to "AL_art == 8"-abstractions, then the values will be taken into account already here
				if (Abstraction[al].AL_art==8)
				{
					flowm3s[0] = Abstraction[al].Abakt[0];
					for (int z=1;z<=tracer_count;z++) 
					{
						flowm3s[z]=Speicher[i].Vakt[z];
					}
					// Werte aus Abstraction from Reservoir in Qges-Array kopieren 
					for (j=1;j<=subint;j++)
					{
						Q12ges[0][j][zonmax+AL[al]]=flowm3s[0];
						for (int z=1;z<=tracer_count;z++) /* Schleife ab z=1 nur fuer die Konzentrationen!! */
						Q12ges[z][j][zonmax+AL[al]]=Qges[z][zarray[TG]];
					}
					Qges[0][zonmax+AL[i]]=flowm3s[0];

					// if (monat>0) Abstraction[i].WriteAbstraction(jahr,monat,tag,stunde,flowm3s);
					if (jahr>0)
					{
						Abstraction[al].WriteAbstraction(jahr,monat,tag,stunde,minute,flowm3s);
					}
				}
			}

			Qges[0][zarray[TG]]=m3s_to_mmh(Speicher[i].Abakt,TgAe/TgAe_rel,InV);

			// Abgabekonzentrationen muessen nun noch mit dem Wert der Konzentrationen im Speicher belegt werden 
			for (int z=1;z<=tracer_count;z++) Qges[z][zarray[TG]]=Speicher[i].Vakt[z];
			// fuer das weitere Routing muessen die Werte noch in Q12ges uebernommen werden
			for (j=1;j<=subint;j++)
			{
				for (int z=0;z<=tracer_count;z++)
				Q12ges[z][j][zarray[TG]]=Qges[z][zarray[TG]];
			}
			if (monat>0)  /* soll in der Initialisierung nicht geschrieben werden */
			{
				// verfügbare Bewässerungsmenge aus Speicher berechnen in mm bezogen auf aktuelles Teilgebiet
				if (i==SPs)
				{
					maxirrigFromReservoirs[0][zarray[TG]] = m3s_to_mmh(Speicher[i].Vakt[0]/60.0/InV,TgAe/TgAe_rel,InV);
					for (int z=1;z<=tracer_count;z++)
					maxirrigFromReservoirs[z][zarray[TG]]=Speicher[i].Vakt[z]; /* concentrations for irrigation water is the same like the mixed concentration of all inflows */
				}
				Speicher[i].WriteSpeicherinhalt(jahr,monat,tag,stunde,minute);
			}
		}
	}
}

double routingclass::GetAbstractionSum()
{
	return abstraction_sum;
}


int routingclass::getmaxdelta(int tmax)
{
	int i,t;
	t=0;
	for (i=1;i<=OLs+SUMTRIBs+ZLs;i++)
	if (wqtab[i].n12t>t) t=wqtab[i].n12t;
	if (t>tmax) return t; else return tmax;
}

void routingclass::SetFluxAdresses(double **p_Qges,double ***p_Q12ges,double ***p_fluss,double **p_Spende,int *p_Spalte)
{
	int i,j,z,k;
	Qges   = p_Qges;
	Q12ges = p_Q12ges;
	fluss  = p_fluss;
	Spende = p_Spende;
	Spalte = p_Spalte;
	for (i=1;i<=ZLs;i++)
	Qges[0][(int)(ZL[i]+zonmax)]=m3s_to_mmh(ZL_wert[0][i],AE[OLs+SUMTRIBs+i],InV);
	for (i=1;i<=OLs+SUMTRIBs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		if (indexST[i])
		{
			for (k=1;k<=indexST[i][0];k++) /* now initializing the rest of the SUMTRIB-tributaries */
			{
				geb = zarray[indexST[i][k]];
				for (z=0;z<=tracer_count;z++)
				{
					if ((fluss[z][geb]=(double*)calloc(1+subint+wqtab[i].n12t,sizeof(double)))==NULL)
					{
						(*WSstream) << "FATAL: error allocating memory for flow-arrays"<<"[end-entry]";
						throw(1);
					}
					if ((gem_Zul>=1)&&(Spalte[geb]!=0)&&(Spende!=NULL))
					{
						if (gem_Zul == 2)
						{
							for (j = 1; j <= subint + wqtab[i].n12t; j++)
							{
								fluss[z][geb][j] = m3s_to_mmh(Spende[10][Spalte[geb]], (abfrel[0][geb] / abfrel[1][geb])* abfrel[0][0] , InV);
							}
						}
						else
						{
							for (j = 1; j <= subint + wqtab[i].n12t; j++)
							{
								fluss[z][geb][j] = Spende[10][Spalte[geb]];
							}
						}
					}
					else
					{
						for (j = 1; j <= subint + wqtab[i].n12t; j++)
						{
							fluss[z][geb][j] = Qges[z][geb];
						}
					}
				}
			}
		}
		geb = i<=OLs+SUMTRIBs ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
			if ((fluss[z][geb]=(double*)calloc(1+subint+wqtab[i].n12t,sizeof(double)))==NULL)
			{
				(*WSstream) << "FATAL: error allocating memory for flow-arrays"<<"[end-entry]";
				throw(1);
			}
			if ((gem_Zul>=1)&&(Spalte[geb]!=0)&&(Spende!=NULL))
			{
				if (gem_Zul == 2)
				{
					for (j = 1; j <= subint + wqtab[i].n12t; j++)
					{
						fluss[z][geb][j] = m3s_to_mmh(Spende[10][Spalte[geb]], (abfrel[0][geb] / abfrel[1][geb])* abfrel[0][0], InV);
					}
				}
				else
				{
					for (j = 1; j <= subint + wqtab[i].n12t; j++)
					{
						fluss[z][geb][j] = Spende[10][Spalte[geb]];
					}
				}
			}
			else
			{
				for (j = 1; j <= subint + wqtab[i].n12t; j++)
				{
					fluss[z][geb][j] = Qges[z][geb];
				}
			}
		}
	}
	for (i=1;i<=ZLs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		geb = ZL[i]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
			fluss[z][geb]=(double*)calloc(1+subint+wqtab[OLs+SUMTRIBs+i].n12t,sizeof(double));
			for (j=1;j<=subint+wqtab[OLs+SUMTRIBs+i].n12t;j++)
			fluss[z][geb][j] =  Qges[z][geb];
		}
	}
}

int routingclass::getbereich(int ol,double q)
/* ol ist nicht die Nummer von OL (also Teilgebietsnummer aus dem Zonengrid),
sondern der Index im OL-Array von 1...OLs */
{
	int b;
	/* umrechnen von mm/Zeitintervall in l/s/km^2,
	dabei ist 1 mm/Stunde = 1000 m^3/km^2/h = 1000000/3600 l/s/km^2 = 277.777778 l/s/km^2 */
	q*=(277.7778*60.0/InV);
	if (q<(double)0.0000001) q=(double)0.0000001; /* damit es beim Logarithmieren keinen "floating point underflow error" gibt */
	/* fuer jede Strecke kann ein anderes Tripel aus
	minspende, maxspende und intervallanzahl intzahl angegeben worden sein,
	weshalb alle drei Werte als Array in Routingclass definiert worden
	sind (jeweils fuer jede Strecke ein Element) */
	b = dq[ol] > 0.0 ? (int)ceil(log(q/minspende[ol])/dq[ol]) : 1;
	if (b<1) b=1; /* weniger als der in minspende angegebene Abfluss soll mit minimaler Geschw. von minspende fliessen */
	if (b>(int)(intzahl[ol]+1)) b=(int)(intzahl[ol]+1); /* ebenso sollen groessere Hochwasser mit max. Geschw. von maxspende fliessen */
	return b;
}

void routingclass::set_initial_last_index()
{
	int i,j,k,ber;
	/* hier werden nun die zum Initialabfluss passenden Verschiebungsintervalle ausgerechnet
	und anschliessend dieser Wert (jeweils die Anzahl von 1/<subint>-tel Zeitintervallen) in
	dem fuer die entsprechende Teilstrecke der gerade bearbeiteten Routingstrecke
	zustaendigen Element von wqtab[strecke].last_index[teilstrecke] abgelegt, so dass
	beim ersten Routen festgestellt werden kann, ob der nach dem Initialabfluss
	kommende erste Modellabfluss die vorher berechneten Abfluesse ueberrollt oder
	hinter ihnen zurueckbleibt. */

	for (i=1;i<=OLs+SUMTRIBs+ZLs;i++) /* einzelne Strecken */
	{
		geb = i<=OLs+SUMTRIBs ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		ber=getbereich(i,Qges[0][geb]);
		for (j=1;j<=wqtab[i].n;j++) /* Teilstrecken innerhalb der einzelnen Strecken */
		{
			for (k=1;k<=subint;k++) wqtab[i].Qaus[j][k]=Qges[0][geb];
			/* beim Unterscheiden in Vorland- und Hauptbettretention muss
		in jedem Zeitintervall fuer den AUslaufterm bekannt sein, wieviel
		Abfluss im vorigen Zeitintervall im Vorland bzw. im Hauptbett abfloss.
		In Qaus ist nur der Gesamtabfluss gespeichert. In
		wqtab[Teilstrecke].QvzuQ[ber] wird fuer jeden Abfluss das Verhaeltnis von Vorland- zu
		gesamtabfluss festgehalten, dieser Wert muss von einem Teilintervall
		ins naechste gerettet werden, um das Auslaufen von Hauptbett und Vorland getrennt
		berechnen zu koennen. Als initialwert wird erstmal null gesetzt (wenn man
		in einem Hochwasser anfaengt zu modelleiren, koennte es anfaenglich zusaetzliche
		Retention durch das AUffuellen der Vorlaender geben, die hier also leer angenommen
		werden*/
			wqtab[i].qvzuq[j]=0.0;

			/* letztes belegtes Element pro Teilstrecke: die ersten subint Elemente sind die
		Breite eines Abflusses, ((j-1)*wqtab[i].delta[1]+subint) gibt den offset zur naechsten
		Teilstrecke an, wqtab[i].delta[ber] ist schliesslich die Verschiebung, die
		fuer den Initialabfluss berechnet wurde (in getbereich()). Der last_index gilt
		fuer eine erfolgte Translation um ein Zeitintervall, d.h., es wird kein
		Zeigerdrehen mehr nach diesen Festlegungen ausgefuehrt, die Speicher
		sind nun im abarbeitungsfaehigen Zustand. */
			wqtab[i].last_index[j] = subint + (j-1)*(wqtab[i].delta[1]+subint) + wqtab[i].delta[ber];
		}
	}
}

void routingclass::Print_W_Q_Table()
{
	char *dateiname;
	int i,j;
	FILE *datei;
	/* In dieser function werden die in SetRoutingParams() berechneten
W_Q_Tabellen fuer jede Strecke in eine Datei geschrieben. Der Dateiname
setzt sich aus wq_table.<streckennummer> zusammen */
	dateiname = (char*)calloc(101,sizeof(char));
	for (i=1;i<=OLs+SUMTRIBs+ZLs;i++)
	{
		geb = i<=(OLs+SUMTRIBs) ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs];
		if (i<=(OLs+SUMTRIBs))
		{
			sprintf(dateiname,"wq_table.%.3d",geb);
		}
		else
		{
			sprintf(dateiname,"wq_tabZL.%.3d",geb);
		}
		if ((datei=fopen(dateiname,"wt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening file "<<dateiname<<" (check user rights or disk space)"<<"[end-entry]";
		}
		else
		{
			if (i<=OLs)
			fprintf(datei,"waterlevel-discharge-table for the channel of subbasin %d to the outlet of subbasin %d\n",OL[i],TG);
			else
			{
				if (i>OLs+SUMTRIBs) /* inflow */
				fprintf(datei,"waterlevel-discharge-table for the channel of inflow %d to the outlet of subbasin %d\n",OL[i],TG);
				else /* SUMTRIB */
				{
					fprintf(datei,"waterlevel-discharge-table for the channel of subbasins %d",OL[i]);
					for (j=1;j<indexST[i][0];j++) fprintf(datei,", %d",indexST[i][j]);
					fprintf(datei," and %d to the outlet of subbasin %d\n",indexST[i][indexST[i][0]],TG);
				}
			}
			fprintf(datei,"input parameters from the control file:\n");
			fprintf(datei,"------------------------------------------\n");
			fprintf(datei,"storage constant mean channel [h] : %7.2f\n",kh[i] != 1.0 ? ( kh[i] > 0.0 ? -1.0/(log(kh[i])) : -9999.0 ): 0.0);
			fprintf(datei,"storage constant flood plain [h]  : %7.2f\n",kv[i] != 1.0 ? ( kv[i] > 0.0 ? -1.0/(log(kv[i])) : -9999.0 ): 0.0);
			fprintf(datei,"width of the mean channel [m]     : %7.2f\n",Bh[i]);
			fprintf(datei,"width of flood plains  [m]        : %7.2f\n",Bv[i]);
			fprintf(datei,"depth of mean channel  [m]        : %7.2f\n",Tm[i]);
			fprintf(datei,"Manning-coef. mean chan[m^(1/3)/s]: %7.2f\n",Mh[i]);
			fprintf(datei,"Manning-ceof. flood pl.[m^(1/3)/s]: %7.2f\n",Mv[i]);
			fprintf(datei,"slope angle            [1/1]      : %7.4f\n",I[i]*I[i]);
			fprintf(datei,"channel length         [m]        : %7.1f\n",L[i]);
			fprintf(datei,"catchment area         [km^2]     : %7.2f\n",AE[i]);
			fprintf(datei,"min. spec. discharge   [l/s/km^2] : %7.1f\n",minspende[i]);
			fprintf(datei,"max. spec. discharge   [l/s/km^2] : %7.1f\n",maxspende[i]);
			fprintf(datei,"intervals for w-q-table[-]        : %7.0f\n",intzahl[i]+1);
			fprintf(datei,"sub-intervals (time)   [-]        : %7d  \n",subint);
			fprintf(datei,"number of sub-channels            : %7.0f\n",wqtab[i].n);
			fprintf(datei,"total number of storage elements  : %7d  \n",wqtab[i].n12t+subint);
			/* Tabellenkopf */
			fprintf(datei,"-------------------------------------------------------------------------------\n"
			"  #    q      q       Q     dm    df     v      vm     vf     t   delta   Qm/Q \n"
			"     [mm]  l/s/km2 [m^3/s]  [m]   [m]   [m/s]  [m/s]  [m/s]  [h]  [Int]   [1/1]\n"
			"-------------------------------------------------------------------------------\n");
			
			for (j=1;j<=(int)(intzahl[i])+1;j++)
			{
				/* Ausfuellen der Tabelle */
				fprintf(datei,"%3d  %5.3f %6.1f %7.2f  %5.2f %5.2f %5.2f  %5.2f  %5.2f  %5.2f  %3d   %5.3f\n",
				j,wqtab[i].sp_mm[j],wqtab[i].sp_lskm2[j],wqtab[i].Q[j],
				wqtab[i].Th[j],wqtab[i].Tv[j],wqtab[i].v[j],wqtab[i].vh[j],wqtab[i].vv[j],
				wqtab[i].t[j],wqtab[i].delta[j],wqtab[i].QvzuQ[j]);
			}
			fclose(datei);
		}
	}
	if (dateiname)
	{
		free(dateiname);
		dateiname = NULL;
	}
}

void routingclass::ClearMemory()
{
	int i,j,z;
	(*WSstream) << " free memory TG "<<TG;

	for (i=1;i<=OLs+SUMTRIBs+ZLs;i++)
	{
		if (indexST[i])
		{
                    for (int k=1;k<=indexST[i][0];k++) /* now initializing the rest of the SUMTRIB-tributaries */
                    {
                        geb = zarray[indexST[i][k]];
                        for (z=0;z<=tracer_count;z++)
                        {
							if (fluss[z][geb])
							{
								free(fluss[z][geb]);
								fluss[z][geb] = NULL;
							}
                        }
                    }
		}

		geb = i<=(OLs+SUMTRIBs) ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		if (i<=OLs+SUMTRIBs)
		{
			(*WSstream) << ", OL "<<i<<" (Subb. "<<OL[i];
			if (indexST[i])
			for (j=1;j<=indexST[i][0];j++)
			(*WSstream) << ", "<<indexST[i][j];
			(*WSstream) << ")";
		}
		else
		(*WSstream) << ", (only ZL's)";

		if (i==1)
		{
			if (Qakt)
			{
				for (z=0;z<=tracer_count;z++)
				{
					if (Qakt[z])
					{
						free(Qakt[z]);
						Qakt[z] = NULL;
					}
				}
				free(Qakt);
				Qakt  = NULL;
			}

			if (total)
			{
				for (z=0;z<=tracer_count;z++)
				{
					if (total[z])
					{
						free(total[z]);
						total[z] = NULL;
					}
				}
				free(total);
				total = NULL;
			}
		}

                if (fluss)
                {
                        for (z=0;z<=tracer_count;z++)
                        {
							if (fluss[z][geb])
							{
								free(fluss[z][geb]);
								fluss[z][geb] = NULL;
							}
                        }
                }
                // fluss darf hier nicht freigegebenw werden, weil es nur ein übergebener Parameter war!
                // die einlenen Array-Elemente wurden jedoch von den einzelnen Instanzen reserviert und muessen also freigegeben werden

		if (wqtab)
		{
			if (wqtab[i].ueberrollt)
			{
				for (z=0;z<=tracer_count;z++)
				{
					if (wqtab[i].ueberrollt[z])
					{
						free(wqtab[i].ueberrollt[z]);
						wqtab[i].ueberrollt[z] = NULL;
					}
				}
				free(wqtab[i].ueberrollt);
				wqtab[i].ueberrollt  = NULL;
			}
			if (wqtab[i].last_index) { free(wqtab[i].last_index); wqtab[i].last_index  = NULL; }
			if (wqtab[i].qvzuq)      { free(wqtab[i].qvzuq);      wqtab[i].qvzuq       = NULL; }
			free_matrix(wqtab[i].Qaus,0,(int)wqtab[i].n+1,0,subint+1); wqtab[i].Qaus = NULL;
			if (wqtab[i].QvzuQ)    { free(wqtab[i].QvzuQ);      wqtab[i].QvzuQ      = NULL;}
			if (wqtab[i].delta)    { free(wqtab[i].delta);      wqtab[i].delta      = NULL;}
			if (wqtab[i].Tv)       { free(wqtab[i].Tv);         wqtab[i].Tv         = NULL;}
			if (wqtab[i].Th)       { free(wqtab[i].Th);         wqtab[i].Th         = NULL;}
			if (wqtab[i].Rv)       { free(wqtab[i].Rv);         wqtab[i].Rv         = NULL;}
			if (wqtab[i].Rh)       { free(wqtab[i].Rh);         wqtab[i].Rh         = NULL;}
			if (wqtab[i].t)        { free(wqtab[i].t);          wqtab[i].t          = NULL;}
			if (wqtab[i].vv)       { free(wqtab[i].vv);         wqtab[i].vv         = NULL;}
			if (wqtab[i].vh)       { free(wqtab[i].vh);         wqtab[i].vh         = NULL;}
			if (wqtab[i].v)        { free(wqtab[i].v);          wqtab[i].v          = NULL;}
			if (wqtab[i].Q)        { free(wqtab[i].Q);          wqtab[i].Q          = NULL;}
			if (wqtab[i].sp_mm)    { free(wqtab[i].sp_mm);      wqtab[i].sp_mm      = NULL;}
			if (wqtab[i].sp_lskm2) { free(wqtab[i].sp_lskm2);   wqtab[i].sp_lskm2   = NULL;}

			if (indexST)
			{
				if (indexST[i])
				{
					free(indexST[i]);
					indexST[i] = NULL;
				}
			}
		}
	}

	if (ZLs>0)
	{
		(*WSstream) << ", --> ZLs...";
		for (i=1;i<=ZLs;i++)
		{
			if (ZL)
			{
				(*WSstream) << " "<<ZL[i]<<",";
			}
		}
		if (ZL_buffer)
		{
			ZL_buffer->ClearMemory();
			delete[] ZL_buffer;
		}
		if (ZL) { free(ZL); ZL = NULL; }
		if (ZL_art) { free(ZL_art); ZL_art = NULL; }
		if (ZL_wert)
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (ZL_wert[z])
				{
					free(ZL_wert[z]);
					ZL_wert[z] = NULL;
				}
			}
			free(ZL_wert);
			ZL_wert=NULL;
		}
	}

	if (ALs>0)
	{
//		if ((OLs+SUMTRIBs+ZLs)==0)
//		{
			if (total)
			{
				for (z=0;z<=tracer_count;z++)
				{
					if (total[z])
					{
						free(total[z]);
						total[z] = NULL;
					}
				}
				free(total);
				total=NULL;
			}
//		}
		(*WSstream) << ", --> ALs...";
		for (i=1;i<=ALs;i++)
		{
			Abstraction[i].AL_flushbuffer(1);
			(*WSstream) << " "<<AL[i]<<",";
		}
		if (AL)          {delete[] AL;          AL=NULL;		 }
		if (Abstraction) {delete[] Abstraction; Abstraction=NULL;}
	}
	if (SPs>0)
	{
		(*WSstream) << ", --> SPs...";
		if (total)
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (total[z])
				{
					free(total[z]);
					total[z] = NULL;
				}
			}
			if (total) { free(total); total = NULL; }
		}
		for (i=1;i<=SPs;i++)
		{

			Speicher[i].SP_flushbuffer(1);
			(*WSstream) << " " << i << " ";

		}
		if (SP)          {delete[] SP;          SP=NULL;	  }
		if (Speicher)    {delete[] Speicher;    Speicher=NULL;}
	}

	(*WSstream) << " Rest...";

	if (wqtab    ) {free(wqtab);     wqtab     = NULL;}
	if (intzahl  ) {free(intzahl);   intzahl   = NULL;}
	if (maxspende) {free(maxspende); maxspende = NULL;}
	if (minspende) {free(minspende); minspende = NULL;}
	if (dq       ) {free(dq);        dq        = NULL;}
	if (AE       ) {free(AE);        AE        = NULL;}
	if (L        ) {free(L);         L         = NULL;}
	if (I        ) {free(I);         I         = NULL;}
	if (Tm       ) {free(Tm);        Tm        = NULL;}
	if (Mv       ) {free(Mv);        Mv        = NULL;}
	if (Mh       ) {free(Mh);        Mh        = NULL;}
	if (Bv       ) {free(Bv);        Bv        = NULL;}
	if (Bh       ) {free(Bh);        Bh        = NULL;}
	if (gew      ) {free(gew);       gew       = NULL;}
	if (OL       ) {free(OL);        OL        = NULL;}
	if (indexST  ) {free(indexST);   indexST   = NULL;}
	if (kv       ) {free(kv);        kv        = NULL;}
	if (kh       ) {free(kh);        kh        = NULL;}

	//  if (flowmmh  ) free(flowmmh);   flowmmh   = NULL;
	if (flowm3s  ) {free(flowm3s);   flowm3s   = NULL;}
	if (CQref    ) {free(CQref);     CQref     = NULL;}
	if (fracht   ) {free(fracht);    fracht    = NULL;}
	if (frachtsum) {free(frachtsum); frachtsum = NULL;}
        modelobject::ClearMemory();
	(*WSstream) << " ok"<<"[end-entry]";
	return;

}

void routingclass::SetZoneSWAdress(statrecord *p_zoneSW)
{
	zoneSW=p_zoneSW; /* contains for each subbasin the amount of irrigation water demand in mm/zone */
}

void routingclass::SetMaxIrrigAdress(double **p_maxirrig,double **p_maxirrigFromReservoirs)
{
	maxirrig=p_maxirrig;
	maxirrigFromReservoirs=p_maxirrigFromReservoirs;
}

void routingclass::WriteStorageFile(FILE *storagefile)
{
	int i,j,k,z,n;
	fprintf(storagefile,"%d (routing subbasin with tributaries",TG);
	for (i=1;i<=OLs;i++) fprintf(storagefile," %d",OL[i]);
	fprintf(storagefile,")\n");

	// new from 10.06.00: some additional parameters to enable the model to converse inputs from different time step configurations


	/* now the content of the river channels and subchannels are stored for each tributary and for each tracer */
	for (i=1;i<=OLs+SUMTRIBs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		if (indexST[i])
		{
			for (k=1;k<=indexST[i][0];k++) /* now initializing the rest of the SUMTRIB-tributaries */
			{
				geb = zarray[indexST[i][k]]; 
				for (z=0;z<=tracer_count;z++)
				{
//					fprintf(storagefile, "%d %d (tracer/flow index and tributary subbasin)\n", z, geb);
					fprintf(storagefile, "\n%d %d %d %d %.0lf (tracer/flow index and tributary subbasin, subint, maxdelta, num-subchannels)\n", z,geb,subint, wqtab[i].delta[1], wqtab[i].n);
					for (j=1;j<=subint+wqtab[i].n12t;j++)
					{
						fprintf(storagefile,"%g ",fluss[z][geb][j]);
//						if (!(j%10)) fprintf(storagefile,"\n");
						if ((j == subint) || ((j - subint) % (wqtab[i].delta[1] + subint) == 0)) fprintf(storagefile, "\n");
					}
					fprintf(storagefile,"\n");
				}
				for (n=1;n<=wqtab[i].n;n++)
				{
					for (j=1;j<=subint;j++)
					fprintf(storagefile,"%g ",wqtab[i].Qaus[n][j]);
					fprintf(storagefile," (wqtab[local tributary %d].Qaus[subchannel %d][1..%d])\n",i,n,subint);
					fprintf(storagefile,"%g (wqtab[local tributary %d].qvzuq[subchannel %d])\n",wqtab[i].qvzuq[n],i,n);
				}
				/* store content of riverflow table */
				if (riverflow_table)
				{
					for (z=0;z<=tracer_count;z++)
					fprintf(storagefile,"%lg ",riverflow_table[geb].flow_mm[z]);
					fprintf(storagefile,"(riverflow_table[%d].flow_mm[%d])\n",geb,z);
					fprintf(storagefile,"%lg %lg (riverflow_table[%d].flowm3s/speed\n",
					riverflow_table[geb].flowm3s,riverflow_table[geb].speed,geb);
				}
			}
		}
		geb = i<=OLs+SUMTRIBs ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
//			fprintf(storagefile,"%d %d (tracer/flow index and tributary subbasin)\n",z,geb);
			fprintf(storagefile, "\n%d %d %d %d %.0lf (tracer/flow index and tributary subbasin, subint, maxdelta, num-subchannels)\n", z, geb, subint, wqtab[i].delta[1], wqtab[i].n);
			for (j=1;j<=subint+wqtab[i].n12t;j++)
			{
				fprintf(storagefile,"%g\t",fluss[z][geb][j]);
//				if (!(j % 10)) fprintf(storagefile, "\n");
				if ((j == subint)||((j-subint) % (wqtab[i].delta[1]+subint)==0)) fprintf(storagefile, "\n");
			}
			fprintf(storagefile,"\n");
		}
		for (n=1;n<=wqtab[i].n;n++)
		{
			for (j=1;j<=subint;j++)
			fprintf(storagefile,"%g ",wqtab[i].Qaus[n][j]);
			fprintf(storagefile," (wqtab[local tributary %d].Qaus[subchannel %d][1..%d])\n",i,n,subint);
			fprintf(storagefile,"%g (wqtab[local tributary %d].qvzuq[subchannel %d])\n",wqtab[i].qvzuq[n],i,n);
		}
		/* store content of riverflow table */
		if (riverflow_table)
		{
			for (z=0;z<=tracer_count;z++)
			fprintf(storagefile,"%lg ",riverflow_table[geb].flow_mm[z]);
			fprintf(storagefile,"(riverflow_table[%d].flow_mm[%d])\n",geb,z);
			fprintf(storagefile,"%lg %lg (riverflow_table[%d].flowm3s/speed\n",
			riverflow_table[geb].flowm3s,riverflow_table[geb].speed,geb);
		}
	}
	for (i=1;i<=ZLs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		geb = ZL[i]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
//			fprintf(storagefile,"%d %d (tracer/flow index and external inflow)\n",z,geb);
			fprintf(storagefile, "\n%d %d %d %d %.0lf (tracer/flow index and external inflow, subint, maxdelta, num-subchannels)\n", z, geb, subint, wqtab[OLs + SUMTRIBs + i].delta[1], wqtab[OLs + SUMTRIBs + i].n);
			for (j=1;j<=subint+wqtab[OLs+SUMTRIBs+i].n12t;j++)
			{
				fprintf(storagefile,"%g ",fluss[z][geb][j]);
//				if (!(j%10)) fprintf(storagefile,"\n");
				if ((j == subint) || ((j - subint) % (wqtab[OLs + SUMTRIBs + i].delta[1] + subint) == 0)) fprintf(storagefile, "\n");
			}
		}
		fprintf(storagefile,"\n");
		for (n=1;n<=wqtab[OLs + SUMTRIBs + i].n;n++)
		{
			for (j=1;j<=subint;j++)
			fprintf(storagefile,"%g ",wqtab[OLs + SUMTRIBs + i].Qaus[n][j]);
			fprintf(storagefile," (wqtab[inflow %d].Qaus[subchannel %d][1..%d])\n",i,n,subint);
			fprintf(storagefile,"%g (wqtab[inflow %d].qvzuq[subchannel %d])\n",wqtab[OLs + SUMTRIBs + i].qvzuq[n],i,n);
		}
	}
	/* store the actual content of reservoirs */
	fprintf(storagefile,"%d (number of reservoirs)\n",SPs);
	for (i=1;i<=SPs;i++)
	{
		for (z = 0; z <= tracer_count; z++)
		{
			fprintf(storagefile, "%lg ", Speicher[i].Vakt[z]);
		}
		fprintf(storagefile," (reservoir %d)\n",i);
		fprintf(storagefile,"%d (useRunningAvg reservoir %d)\n",Speicher[i].UseRunningAvg, i);
		if (Speicher[i].UseRunningAvg > 0)
		{
			for (std::list<double>::iterator iter = Speicher[i].V_hist.begin(); iter != Speicher[i].V_hist.end(); iter++)
			{
				fprintf(storagefile, "%lg ", double(*iter));
			}
		}
	}
	fprintf(storagefile,"\n\n");
}

void routingclass::ReadStorageFile(FILE *storagefile, char *storagefilename, double timestepstor ,int syntaxVersion = 2)
{
	int i,j,k,z,n,st_int1,st_int2;
	char *tstr;
	long filepos;

	int src_subint, src_maxdelta, src_numchannels;
	bool convertsource = false; // if set to true during analysing the input storage file: the read storage file values are not read in 1:1 into the fluss array etc but the actual structure is filles based on a mapping of the source input to the actual structure

	filepos=0;
	tstr=(char*)calloc(_stringsize_+1,sizeof(char));
	fscanf(storagefile,"%d",&st_int1); //routing subbasin with tributaries
	tstr=fgets(tstr,_stringsize_,storagefile);
	filepos++;
	if (st_int1!=TG) 
	{
		(*WSstream) << "FATAL: error reading subbasin code in file "<<storagefilename
		<<" (position "<<filepos<<")"<<"[end-entry]";
		throw(71);
	}

	/* now the content of the river channels and subchannels are stored for each tributary and for each tracer */
	for (i=1;i<=OLs+SUMTRIBs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		if (indexST[i])
		{
			for (k=1;k<=indexST[i][0];k++) /* now initializing the rest of the SUMTRIB-tributaries */
			{
				geb = zarray[indexST[i][k]]; 
				for (z=0;z<=tracer_count;z++)
				{
					fscanf(storagefile,"%d%d",&st_int1,&st_int2); //(tracer/flow index and tributary subbasin)

					if (st_int1!=z) 
					{
						(*WSstream) << "FATAL: error reading tracer code (expected " << z << ", read " << st_int1 << ") in file "<<storagefilename	<<" (position "<<filepos<<")"<<"[end-entry]";
						throw(71);
					} 
					if (st_int2!=geb) 
					{
						(*WSstream) << "FATAL: error reading tributary index (expected " << geb << ", read " << st_int2 << ") in file "<<storagefilename	<<" (position "<<filepos<<")"<<"[end-entry]";
						throw(71);
					} 

					if (syntaxVersion == 2)
					{
						if (fscanf(storagefile, "%d%d%d", &src_subint, &src_maxdelta, &src_numchannels) != 3)
						{
							(*WSstream) << "FATAL: error reading routing_model_v2 extensions src_subint, src_maxdelta and src_numchannels for routing channel "<<geb<<" in file " << storagefilename << " (position " << filepos << ")" << "[end-entry]";
							throw(71);
						}
						if ((subint!=src_subint) || (wqtab[i].delta[1]!=src_maxdelta) || (src_numchannels != (int)wqtab[i].n))
						{
							convertsource = true;
						}
					}
					tstr = fgets(tstr, _stringsize_, storagefile);

					filepos++;

					if (convertsource)
					{
						ConvertStorageFile(storagefile, z, i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor);
					}
					else
					{
						// reading the values into the fluss array: first subint values that represent Qaus, then for each sub-channel maxdelta+subint values (wqtab[i].n12t has as many elements as maxdelt+subint times number of subchannels)
						for (j = 1; j <= subint + wqtab[i].n12t; j++)
						{
							fscanf(storagefile, "%lf", &(fluss[z][geb][j]));
						}
					}
				}
				if (convertsource)
				{
					ConvertStorageFile(storagefile, z, i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor,2);
				}
				else
				{
					for (n = 1; n <= wqtab[i].n; n++)
					{
						for (j = 1; j <= subint; j++)
						{
							fscanf(storagefile, "%lf", &(wqtab[i].Qaus[n][j]));
						}
						tstr = fgets(tstr, _stringsize_, storagefile);
						filepos++;
						fscanf(storagefile, "%lf", &(wqtab[i].qvzuq[n])); //(wqtab[local tributary].qvzuq[subchannel])
						tstr = fgets(tstr, _stringsize_, storagefile);
						filepos++;
					}
				}
				/* store content of riverflow table */
				if (riverflow_table)
				{
					for (z = 0; z <= tracer_count; z++)                   //(riverflow_table[%d].flow_mm[%d])\n
					{
						fscanf(storagefile, "%lf", &(riverflow_table[geb].flow_mm[z]));
					}
					tstr=fgets(tstr,_stringsize_,storagefile);
					filepos++;
					fscanf(storagefile,"%lf%lf",&(riverflow_table[geb].flowm3s),&(riverflow_table[geb].speed)); //(riverflow_table[%d].flowm3s/speed\n
					tstr=fgets(tstr,_stringsize_,storagefile);
					filepos++;
				}
			}
		}
		geb = i<=OLs+SUMTRIBs ? zarray[OL[i]] : ZL[i-OLs-SUMTRIBs]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
			fscanf(storagefile,"%d%d",&st_int1,&st_int2); //(tracer/flow index and tributary subbasin)

			if (st_int1!=z) 
			{
				(*WSstream) << "FATAL: error reading tracer code (expected "<< z <<", read " << st_int1 << ") in file "<<storagefilename<<" (position "<<filepos<<")"<<"[end-entry]"; 
				throw(71);
			} 
			if (st_int2!=geb) 
			{
				(*WSstream) << "FATAL: error reading tributary index (expected " << geb << ", read " << st_int2 << ") in file "<<storagefilename<<" (position "<<filepos<<")"<<"[end-entry]";
				throw(71);
			} 

			if (syntaxVersion == 2)
			{
				if (fscanf(storagefile, "%d%d%d", &src_subint, &src_maxdelta, &src_numchannels) != 3)
				{
					(*WSstream) << "FATAL: error reading routing_model_v2 extensions src_subint, src_maxdelta and src_numchannels for routing channel " << geb << " in file " << storagefilename << " (position " << filepos << ")" << "[end-entry]";
					throw(71);
				}
				if ((subint != src_subint) || (wqtab[i].delta[1] != src_maxdelta) || (src_numchannels != (int)wqtab[i].n))
				{
					convertsource = true;
				}
			}
			tstr = fgets(tstr, _stringsize_, storagefile);

			filepos++;

			if (convertsource)
			{
				ConvertStorageFile(storagefile, z, i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor);
			}
			else
			{
				// reading the values into the fluss array: first subint values that represent Qaus, then for each sub-channel maxdelta+subint values (wqtab[i].n12t has as many elements as maxdelt+subint times number of subchannels)
				for (j = 1; j <= subint + wqtab[i].n12t; j++)
				{
					fscanf(storagefile, "%lf", &(fluss[z][geb][j]));
				}
			}
		}
		if (convertsource)
		{
			ConvertStorageFile(storagefile, z, i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor, 2);
		}
		else
		{
			for (n = 1; n <= wqtab[i].n; n++)
			{
				for (j = 1; j <= subint; j++)
				{
					fscanf(storagefile, "%lf", &(wqtab[i].Qaus[n][j]));
				}
				tstr = fgets(tstr, _stringsize_, storagefile);
				filepos++;
				fscanf(storagefile, "%lf", &(wqtab[i].qvzuq[n])); //(wqtab[local tributary].qvzuq[subchannel])
				tstr = fgets(tstr, _stringsize_, storagefile);
				filepos++;
			}
		}
		/* store content of riverflow table */
		if (riverflow_table)
		{
			for (z = 0; z <= tracer_count; z++)                   //(riverflow_table[%d].flow_mm[%d])\n
			{
				fscanf(storagefile, "%lf", &(riverflow_table[geb].flow_mm[z]));
			}
			tstr=fgets(tstr,_stringsize_,storagefile);
			filepos++;
			fscanf(storagefile,"%lf%lf",&(riverflow_table[geb].flowm3s),&(riverflow_table[geb].speed)); //(riverflow_table[%d].flowm3s/speed\n
			tstr=fgets(tstr,_stringsize_,storagefile);
			filepos++;
		}
	}
	for (i=1;i<=ZLs;i++) /* Initialabfluss in fluss-array setzen fuer jeden Oberlieger */
	{
		geb = ZL[i]+zonmax; /* wenn ein Zulauf geroutet wird, soll keine echte TG-Nummer, sondern die um die TG-Anzahl erhoehte ZL-Nummer genommen werden */
		for (z=0;z<=tracer_count;z++)
		{
			fscanf(storagefile,"%d%d",&st_int1,&st_int2); //(tracer/flow index and tributary subbasin)

			if (st_int1!=z) 
			{
				(*WSstream) << "FATAL: error reading tracer code (expected " << z << ", read " << st_int1 << ") in file "<<storagefilename<<" (position "<<filepos<<")"<<"[end-entry]";
				throw(71);
			} 
			if (st_int2!=geb) 
			{
				(*WSstream) << "FATAL: error reading tributary index (expected " << geb << ", read " << st_int2 << ") in file "<<storagefilename<<" (position "<<filepos<<")"<<"[end-entry]";
				throw(71);
			} 

			if (syntaxVersion == 2)
			{
				if (fscanf(storagefile, "%d%d%d", &src_subint, &src_maxdelta, &src_numchannels) != 3)
				{
					(*WSstream) << "FATAL: error reading routing_model_v2 extensions src_subint, src_maxdelta and src_numchannels for routing channel " << geb << " in file " << storagefilename << " (position " << filepos << ")" << "[end-entry]";
					throw(71);
				}
				if ((subint != src_subint) || (wqtab[OLs + SUMTRIBs + i].delta[1] != src_maxdelta) || (src_numchannels != (int)wqtab[OLs + SUMTRIBs + i].n))
				{
					convertsource = true;
				}
			}
			tstr = fgets(tstr, _stringsize_, storagefile);

			filepos++;

			if (convertsource)
			{
				ConvertStorageFile(storagefile, z, OLs + SUMTRIBs + i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor);
			}
			else
			{
				// reading the values into the fluss array: first subint values that represent Qaus, then for each sub-channel maxdelta+subint values (wqtab[i].n12t has as many elements as maxdelt+subint times number of subchannels)
				for (j = 1; j <= subint + wqtab[OLs + SUMTRIBs + i].n12t; j++)
				{
					fscanf(storagefile, "%lf", &(fluss[z][geb][j]));
				}
			}
		}
//		fprintf(storagefile,"\n");
		if (convertsource)
		{
			ConvertStorageFile(storagefile, z, OLs + SUMTRIBs + i, geb, src_subint, src_maxdelta, src_numchannels, timestepstor, 2);
		}
		else
		{
			if (syntaxVersion == 2)
			{
				for (n = 1; n <= wqtab[OLs + SUMTRIBs + i].n; n++)
				{
					for (j = 1; j <= subint; j++)
					{
						fscanf(storagefile, "%lf", &(wqtab[OLs + SUMTRIBs + i].Qaus[n][j]));
					}
					tstr = fgets(tstr, _stringsize_, storagefile);
					filepos++;
					fscanf(storagefile, "%lf", &(wqtab[OLs + SUMTRIBs + i].qvzuq[n])); //(wqtab[local tributary].qvzuq[subchannel])
					tstr = fgets(tstr, _stringsize_, storagefile);
					filepos++;
				}
			}
			else // to be comatible with old storage files, we will read the wrong content here, despite it's erronous, but it's at least compatible with old wasim versions
			{
				for (n = 1; n <= wqtab[i].n; n++)
				{
					for (j = 1; j <= subint; j++)
					{
						fscanf(storagefile, "%lf", &(wqtab[i].Qaus[n][j]));
					}
					tstr = fgets(tstr, _stringsize_, storagefile);
					filepos++;
					fscanf(storagefile, "%lf", &(wqtab[i].qvzuq[n])); //(wqtab[local tributary].qvzuq[subchannel])
					tstr = fgets(tstr, _stringsize_, storagefile);
					filepos++;
				}
			}

		}
	}
	/* store the actual content of reservoirs */
	fscanf(storagefile,"%d",&st_int1); // (number of reservoirs)\n
	tstr=fgets(tstr,_stringsize_,storagefile);
	filepos++;
	if (st_int1!=SPs) 
	{
		(*WSstream) << "FATAL: error reading reservoir code in file "<<storagefilename<<" (position "<<filepos<<")"<<"[end-entry]"; 
		throw(71);
	} 
	for (i=1;i<=SPs;i++)
	{
		for (z=0;z<=tracer_count;z++)
		fscanf(storagefile,"%lf ",&(Speicher[i].Vakt[z]));
		tstr=fgets(tstr,_stringsize_,storagefile);
		filepos++;
		unsigned long actfilepos = ftell(storagefile);
		// sample for running average volume list for a big lake
		// 48 (useRunningAvg reservoir 1)
		// 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.698e+10 3.69805e+10 3.69812e+10 3.69826e+10 3.69851e+10 3.69892e+10 3.69955e+10 3.70049e+10 3.7016e+10 3.70287e+10 3.70424e+10 3.7057e+10 3.70721e+10 3.70876e+10 3.71033e+10 3.71191e+10 3.71349e+10 3.71507e+10 3.71662e+10 3.71817e+10 3.71969e+10 3.72121e+10 3.7227e+10 3.72418e+10
		fscanf(storagefile, "%d", &st_int1); // (number of reservoirs)\n
		tstr = fgets(tstr, _stringsize_, storagefile);
		if (std::string(tstr).find("useRunningAvg") != std::string::npos)
		{
			sscanf(tstr, "%d", &(Speicher[i].UseRunningAvg));
			if (Speicher[i].UseRunningAvg > 0)
			{
				tstr = fgets(tstr, _stringsize_, storagefile);
				std::stringstream ss;
				ss << tstr;
				int i = 0;
				double tempval = 0;
				while (ss.good())
				{
					ss >> tempval;
					if (ss.good())
					{
						Speicher[i].V_hist.push_back(tempval);
					}
				}
//				for (int j = 1; j < Speicher[i].UseRunningAvg; j++)
//				{
//					if (sscanf(tstr, "%lg", &tempval) == 1)
//					{
//						Speicher[i].V_hist.push_back(tempval);
//						std::string mystring(tstr);
//						
//					}
//				}
			}
			
		}
		else
		{
			fseek(storagefile, actfilepos, 0);
		}


	}
	if (tstr) { free(tstr); tstr = NULL; }
}

void routingclass::ConvertStorageFile(FILE * storagefile, int z, int tribindex, int geb,int src_subint, int src_maxdelta, int src_numchannels, double src_timestep,int mode)
{
	double value;
	
	int trg_numchannels = (int)wqtab[tribindex].n;
	int trg_subint = subint;
	int trg_maxdelta = wqtab[tribindex].delta[1];

	if (mode == 1)
	{

		// first line in v2 syntax contains irrelevant values only (routed results for the last sub channel), hence we can skip that line
		fgets(satz, _stringsize_, storagefile);

		// now each line contains src_maxdelta" + "src_subint" values. We need to loop through the source and target subchannels and distribute the values
		// easy approach as first estimate: take the first value of each input row and populate as many target sub channels as appropriate for an even redistribution
		int numValues = 0;
		double sumValues = 0;
		double subchannelFactor = double(trg_numchannels) / double(src_numchannels);
		for (int i = 1; i <= src_numchannels; i++)
		{
			fscanf(storagefile, "%lf", &value);
			tstr = fgets(tstr, _stringsize_, storagefile);
			sumValues += value;
			numValues++;

			int startchannel = (int)floor((double)(i - 1) * subchannelFactor);
			int stopchannel = (int)floor((double)(i)* subchannelFactor);

			int startindex = subint + 1 + startchannel * (trg_subint + trg_maxdelta);
			int stopindex = subint + stopchannel * (trg_subint + trg_maxdelta);

			if (startindex <= stopindex)
			{
				// only if the stopindex is larger than the startindex, the initialization is done now
				for (int j = startindex; j <= stopindex; j++)
				{
					fluss[z][geb][j] = (sumValues / numValues) * (timestep / src_timestep);
				}
				sumValues = 0;
				numValues = 0;
			}
		}
	}
	else if (mode == 2)
	{
		// reading Qaus-Values and accumulating oder spreading, whatever is needed by the different counters for subchannels and subint
		// every other line contains "src_subint" values for Qaus. We need to loop through the source and target subchannels and distribute the values
		// easy approach as first estimate: take the first value of each input row and populate as many target Qaus-Values as appropriate for an even redistribution
		double src_qvzuq = 0;
		int numValues = 0;
		double sumValues = 0;
		double subchannelFactor = double(trg_numchannels) / double(src_numchannels);
		for (int i = 1; i <= src_numchannels; i++)
		{
			fscanf(storagefile, "%lf", &value);
			tstr = fgets(tstr, _stringsize_, storagefile);
			sumValues += value;
			numValues++;
			fscanf(storagefile, "%lf", &src_qvzuq);
			tstr = fgets(tstr, _stringsize_, storagefile);

			int startchannel = (int)floor((double)(i - 1) * subchannelFactor)+1; // first index is 1
			int stopchannel = (int)floor((double)(i)* subchannelFactor);

//			int startindex = subint + 1 + startchannel * (trg_subint + trg_maxdelta);
//			int stopindex = subint + stopchannel * (trg_subint + trg_maxdelta);

			if (startchannel <= stopchannel)
			{
				// only if the stopindex is larger than the startindex, the initialization is done now
				for (int j = startchannel; j <= stopchannel; j++)
				{
					for (int k = 1; k <= trg_subint; k++)
					{
						wqtab[tribindex].Qaus[j][k] = (sumValues / numValues) * (timestep / src_timestep);
					}
				}
				sumValues = 0;
				numValues = 0;
			}
		}

	}
}

/*
routingclass::~routingclass()
{
ClearMemory();
}
*/
/*~~~~~~~~~~~~~~~~~~~~ Ende der Routing - routinen ~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~ Beginn der Abflussmodellroutinen ~~~~~~~~~~~~~~~~~ */

void abflussmodellklasse::SetRiverflow_tableAdress(riverflow_table_type *adresse)
{
	int i;
	if (do_it==0) 
	{
		return;
	}
	else
	{
		riverflow_table=adresse;
		for (i=1;i<=anzrout;i++)
		routing[i]->riverflow_table_ini(adresse,i);
	}
}

int abflussmodellklasse::count_substr(char *satz,const char *substring)
{
	int i;
	char *temp,*temp2;
	temp=newstring(strlen(satz)+2,1);
	temp2=newstring(strlen(satz)+2,1);
	strcpy(temp,satz);
	/* Anzahl substring in satz zaehlen: */
	i=0;
	while (strstr(temp,substring))
	{
		strcpy(temp2,temp);
		strcpy(temp,&(temp2[(int)(strstr(temp2,substring)-temp2)+1]));
		i++;
	}
	if (temp)  { free(temp); temp=NULL;	  }
	if (temp2) { free(temp2); temp2=NULL; }
	return i;
}

void abflussmodellklasse::CountMPICells(bool bLinterceptionBeforeSnow)
{
	if (do_it_mpi0>0)
	{
		if (mpi_rank>0)
		{
			// this call was not done for any instance with mpi_rank!=0, so we need to do it here in order to get the status of the lake model (lake_active), even if the instance has do_it==0
			initialize_lakes(_lake_model_);
		}
		if (lake_active)
		{
			lakeCellsPerMPIRank = new long[mpi_size + 1];
			for (int i=0;i<=mpi_size;i++)
			{
				lakeCellsPerMPIRank[i]=0;
			}
			// getting the indices of grids which must be gathered in each time step (when running in MPI mode)
			if (ist_vorhanden(_LAKE_CODES_)) 
			{
				lakeindex=getindex(_LAKE_CODES_);
				double **LAKE = gridinfo[lakeindex].pointer;

				pondindex=getindex(_ponding_storage_);
				char name[200];
				sprintf(name,"%s%d",_Reale_Verdunstung_,0);
				etrindex=getindex(name);

				sprintf(name, "%s%d", _Interzeptionsverdunstung_, 0);
				eiindex = getindex(name);

				sprintf(name,"%s%s",_soil_bottom_fluxes_,"ext");
				if (ist_vorhanden(name))
				{
					quindex=getindex(name);
				}

				if (ist_vorhanden(_Basisabfluss_))
				{
					qbindex=getindex(_Basisabfluss_);
				}

				
				if (interceptionBeforeSnow)
				{
					sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
					pindex = getindex_first_of_3(_Schneeausfluss_, name, _Niederschlag_);
				}
				else
				{
					sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
					pindex = getindex_first_of_3(name, _Schneeausfluss_, _Niederschlag_);
				}

				// todo: hier noch alle conc_stor_unsatzon-n gathern
				for (int z=1;z<=tracer_count;z++)
				{
					sprintf(name,"concstor_tracer_%d_%s",z,_Unsatzonmodell_);
					tracerindex[z]=getindex(name);
				}

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

//				if (mpi_rank == 0) can be used when the entire grid is available. Here, we reduce the active portion to the actual MPI rank's cells
//				{
					for (int i = mpistartrow; i <= mpiendrow; i++)
					{
						for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
						{
//							if ((zn[i][j] > nodata) && (MPIREG[i][j] > nodata) && (LAKE[i][j] > nodata)) would be used if the entire grid is processed
							if ((zn[i][j] > nodata) && (MPIREG[i][j] == mpi_rank) && (LAKE[i][j] > nodata))
							{
								lakeCellsPerMPIRank[(long)MPIREG[i][j]]++;
							}
						}
					}
//				}


#ifdef _MPI
				long *rank0_lakeCellsPerMPIRank = new long[mpi_size + 1];
				MPI_Reduce(&(lakeCellsPerMPIRank[0]), &(rank0_lakeCellsPerMPIRank[0]), mpi_size, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				for (int i = 0; i < mpi_size; i++)
				{
					lakeCellsPerMPIRank[i] = rank0_lakeCellsPerMPIRank[i];
				}
				MPI_Bcast(&(lakeCellsPerMPIRank[0]), mpi_size, MPI_LONG, 0,MPI_COMM_WORLD);
				delete[] rank0_lakeCellsPerMPIRank;
#endif
//cout << "---------------------- rank " << mpi_rank << " nach Broadcast lakeCellsPerMPIRank" << endl;
				for (int i = 0; i<mpi_size; i++)
				{
					(*WSstream) << "DEBUG: lakeCellsPerMPIRank["<<i<<"] "<<lakeCellsPerMPIRank[i]<<"[end-entry]";
//cout << "rank " << mpi_rank << " lakeCellsPerMPIRank[" << i << "]=" << lakeCellsPerMPIRank[i] << endl;
				}
				MPIScatterGridWithMask(pondindex,lakeindex,lakeCellsPerMPIRank);
			}
		}
	}
}

/* Initialisierung (Constructor) der Abflusskonzentrationsroutine */
abflussmodellklasse::abflussmodellklasse(int &retvalue,wasim::mystream *Stream,TLastResults *p_LastResults,int maxLastResults,double breitearg,double laengearg,double bezugsmeridianarg,
	char *bezeichnung,char *steuername,int LogLevel,int progaus,double p_min_to_int_end, double MinutesToInvEndRouting, int externalTimeStep,
	int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,
	firstlast *p_FirstLastCol ,int p_zonmax,int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata,
	long mina, long sta, long tga, long mta, long jra,
	long mine, long ste, long tge, long mte, long jre,
	double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	int p_tracer_count, tracer_info_type *p_tracer_info_table,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int p_readgrids):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols,p_nrows,p_mpistartrow,p_mpiendrow,p_xll,p_yll,p_cellsize,p_nodata,mine,ste,tge,mte,jre,p_timekorr,p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids,WriteDoubleGrids,WriteNetcdfGrids, DefaultInputDirectory,DefaultOutputDirectory,v_mpiregioncnt,p_sendHalo,p_recvHalo,
	p_doSpinUp, p_routingmode, p_vHRU),
storagefilename(0),storagefile(0),riverflow_table(0),anz_OL(0),anz_ZL(0),anz_AL(0),anz_SP(0),
anz_SUMTRIB(0),anz_et(0),p__OL(0),p__ZL(0),p__AL(0),p__SP(0),p__et(0),p__SUMTRIB(0),num_ST(0),
Qstatfilename(0),spendenfilename(0),satz1(0),tsatz1(0),satz2(0),tsatz2(0),satz3(0),tsatz3(0),
//kombisatz(0),kombisatztemp(0),tkombisatz(0),numcoldescriptors(0),kombisatzlaenge(0),Qstatcode(0),
numcoldescriptors(0),Qstatcode(0),
zoneSW(0),maxirrig(0),maxirrigFromReservoirs(0),maxirrig_temp(0),rtarget(0),timeoffset(0),time(0),spendenfile(0),Qstat(0),
gem_Zul(0),anzrout(0),abfrel(0),fluss(0),Qges(0),Qges_uebernommen(0),Q12ges(0),Spende(0),
Ys(0),Ms(0),Ds(0),Hs(0),Spalte(0),tmax(0),spendenende(0),spendenanzahl(0),p_subint(0),
p_minspende(0),p_maxspende(0),p_intzahl(0),p_kv(0),p_kh(0),p_Bv(0),p_Bh(0),p_Mv(0),p_Mh(0),
p_Th(0),p_I(0),p_L(0),p_AE(0),p_AErel(0),routing(0),runoff_measure(0),effizienz(0),log_effi(0),
SpendenBuffer(0),lake_active(0),lakeCellsPerMPIRank(0),lakeindex(0), pondindex(0), etrindex(0), quindex(0), qbindex(0), pindex(0),do_it_mpi0(0),
preconlakesindex(0), PRECONLAKES(NULL), interceptionBeforeSnow(false), lastday(-1),additionalinfo(NULL),
MinToInvEnd(MinutesToInvEndRouting),extTimestep(externalTimeStep)
{
	try
	{
		retvalue=0;
		int k,tg,z,n;

		char *ZL_modus,*ZL_filename,*SP_filename; //*AL_modus,*AL_filename,*AL_einheit,
		long  ZL_reihe;
		int  ZL_spalte,ZL_art;//,AL_art;
		double ZL_wert;// ,AL_max,AL_schwelle,AL_anteil;

		additionalinfo = new char[255];

		if (mpi_rank!=0)
		{
			do_it_mpi0=do_it;
			do_it=0;
		}
		else
		{
			do_it_mpi0=1;
		}

		if (do_it == 0) return;
		/* Einlesen der Parameter und Namen aus der Steuerdatei */
#ifdef __TIME__
		// (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
#endif
		gem_Zul = do_it-1;

		timestepfactor = extTimestep / InV;

		//		gem_Zul = 0;
		//		if (do_it == 2) gem_Zul = 1;
		//		if (do_it == 3) gem_Zul = 2;

		/* die relativen Anteile der Teilgebiete am Gesamtgebiet werden extra gefuehrt, da sie fuer ganze Pegelgebiete
	groesser werden. Die Gesamtsumme der Relativanteile wird etwa 2..3 sein, da viele Pegel doppelt und
	dreifach zaehlen.
	*/

		num_ST=(int*)calloc(10,sizeof(int));

		abfrel = (double**)calloc(2,sizeof(double*));
		abfrel[0] = (double*)calloc(zonmax+1,sizeof(double));
		abfrel[1] = (double*)calloc(zonmax+1,sizeof(double));
		for (int i=0;i<=zonmax;i++) 
		{
			abfrel[0][i]=zone[i].rel;    /* in element 0: total area of the entire basin in km^2 */
			abfrel[1][i]=1.0;            /* all relative subbasin fractions are set to initial 1 */
		}

		rtarget=(int**)calloc(3,sizeof(int*));
		rtarget[0]=(int*)calloc(2*zonmax+1,sizeof(int)); /* order of routing */
		rtarget[1]=(int*)calloc(zonmax+1,sizeof(int)); /* targets of subbasins */
		rtarget[2]=(int*)calloc(zonmax+1,sizeof(int)); /* index of the routing object instance for getting relative area for m3/s output */

		/* Spalte haelt fuer jedes Teilgebiet fest, welche Spalte der Spendendatei mit ddem Teilgebiet korrespondiert,
	diese Werte nuessen in der STeuerdatei durch blanks getrennt fuer jedes Teilgebiet angegeben werden.
	Hat ein Teilgebiet keine korrespondierende Spalte in der Spendendatei, dann wird 0 eingetragen */
		Spalte = (int*)calloc(zonmax+1,sizeof(int));
		for (int i=1;i<=zonmax;i++) Spalte[i]=0;

		LastResults = p_LastResults;
		
		/* Dieses Teilmodell nimmt Anfangszeiten entgegen, die es zum Suchen in der Spendendatei braucht, alle anderen Teilmodelle
	dagegen nehmen Endzeiten entgegen, um ihre Berechnung beenden zu koennen */
		jahra   = jra;
		monata  = mta;
		stundea = sta;
		taga    = tga;
		minutea = mina;

		readgrids = p_readgrids;

		storagefilename=NULL; /* this filename is initialized and allocated in the routine "ReadStorageFile". If the routine is not 
	called it must be initialized with zero in order to avoid the writing of this file */
		Qstatfilename   = newstring(_stringsize_+1,1);
		spendenfilename = newstring(_stringsize_+1,1);
		satz1           = newstring(10*_stringsize_+1,1);  tsatz1 = satz1;
		satz2           = newstring(10*_stringsize_+1,1);  tsatz2 = satz2;
		satz3			= newstring(10*_stringsize_+1,1);  tsatz3 = satz3;

		if (fscanf(steuerdatei,"%lf%lf%lf%d",&p_minspende,&p_maxspende,&p_intzahl,&p_subint)!=4)
		{
			(*WSstream) << "FATAL: values for min/max specific discharge and/or discharge and time interval splitting   \n"
			"not found in control file"<<"[end-entry]";
			throw(17);
		}
		tstr=fgets(tstr,_stringsize_,steuerdatei);

		fscanfln(steuerdatei,"%s%d",Qstatfilename,&Qstatcode  ,"name of the Qstatfiles for ",tstr);
		fscanfln(steuerdatei,"%s",spendenfilename,"name of spec.discharge files for "       ,tstr);
		fscanfln(steuerdatei,"%d",&numcoldescriptors,"number column descriptors for"        ,tstr);

		for (int i=1;i<=numcoldescriptors;i++)
		{
			if (i > zonmax)
			{
				(*WSstream) << "FATAL: too much column descriptors ("<<numcoldescriptors<<") in routing model (more than sub basins: "<< zonmax <<")" << "[end-entry]";
				throw(17);
			}

			if (fscanf(steuerdatei,"%d",&k)!=1)
			{
				(*WSstream) << "FATAL: error reading a subbasin index for column indexing in routing model"<<"[end-entry]";
				throw(17);
			}
			if (fscanf(steuerdatei,"%d",&(Spalte[zarray[k]]))!=1)
			{
				(*WSstream) << "FATAL: error reading a column number for column indexing in routing model"<<"[end-entry]";
				throw(17);
			}
		}
		if (numcoldescriptors>0) satz1=fgets(satz1,_stringsize_,steuerdatei);
		fscanfln(steuerdatei,"%d",&timeoffset,"Zeitoffset fuer Beginn der R-square-Berechnung fuer ",tstr);

		/* routing ist eine Struktur, die fuer jede Gruppe von Teilgebieten, die zusammengefasstwerden sollen,
	einen Parametersatz speichert (Anzahl Teilgebiete (OLs), Nummern der Teilgebiete (OL), Nummer des Gesamt(Teil)gebietes (TG),
	die Routingparameter fuer jedes OL (kh, kv, Th, Mh, Mv, Bh, Bv, I, L, AE usw.) sowie die Gewichte der OL und des TG relativ zum Gesamt-TG.
	Beim Routing muss in einer Reihenfolge vorgegangen werden, dass immer zuerst die obersten Teilgebiete geroutet werden.
	Die Ergebnisse dieses Routings gehen dann in die unteren Teilgebiets-routings ein.
	Zunaechst muss die Anzahl von Routingvorschriften gezaehlt werden (alle Zeilen, die nach den ersten 7 Eintraegen im
	Abflussmodell stehen, sind Routingvorschriften, bis eine Leerzeile kommt. Ein neues Teilgebiet wird an "TG" erkannt
	*/
		anzrout=0;
		anz_OL=0; anz_ZL=0; anz_AL=0; anz_SP=0; anz_SUMTRIB=0; anz_et=0;
		if (!feof(steuerdatei)) tstr=fgets(tstr,_stringsize_,steuerdatei);
		else
		{
			(*WSstream) << "FATAL: error reading a routing description from control file! (eof reached)."<<"[end-entry]";
			throw(38);
		}
		while ((strstr(tstr,_TG_))||(strstr(tstr,_OL_))||(strstr(tstr,_SUMTRIB_))||(strstr(tstr,_ZL_))||(strstr(tstr,_AL_))||(strstr(tstr,_SP_)))
		{
			if (strstr(tstr,_TG_)) anzrout++; /* Routingvorschriften fuer ein neues TG werden nun definiert */
			anz_OL+=count_substr(tstr,_OL_);
			anz_et+=count_substr(tstr,"&");
			anz_SUMTRIB+=count_substr(tstr,_SUMTRIB_);
			anz_ZL+=count_substr(tstr,_ZL_);
			anz_AL+=count_substr(tstr,_AL_);
			anz_SP+=count_substr(tstr,_SP_);
			if (!feof(steuerdatei)) tstr=fgets(tstr,_stringsize_,steuerdatei);
			else
			{
				(*WSstream) << "FATAL: error reading a routing description from control file! (eof reached)."<<"[end-entry]";
				throw(38);
			}
		}
		(*WSstream) << "number inbetween areas (routings)    : "<<anzrout<<"[end-entry]";;
		(*WSstream) << "number tributaries  total      (OL)  : "<<(anz_OL+anz_et+anz_SUMTRIB)<<"[end-entry]";
		(*WSstream) << " - flowing into another at the end   : "<<anz_OL<<"[end-entry]";
		(*WSstream) << " - flowing into another at the start : "<<(anz_et+anz_SUMTRIB)<<"[end-entry]";
		(*WSstream) << "so the total of routingchannels is   : "<<(anz_OL+anz_SUMTRIB)<<"[end-entry]";
		(*WSstream) << "number of inflows              (ZL)  : "<<anz_ZL<<"[end-entry]";
		(*WSstream) << "number of abstractions         (AL)  : "<<anz_AL<<"[end-entry]";
		(*WSstream) << "number reservoirs              (SP)  : "<<anz_SP<<"[end-entry]";

		substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (retvalue!=0) throw retvalue;

		/* maxirrig: maximum amount of irrigation water from surface water sources (rivers or reservoirs) in m3/s per zone */
		maxirrig    = (double**)calloc(tracer_count+1,sizeof(double*));
		for (z=0;z<=tracer_count;z++)
		maxirrig[z] = (double*)calloc(zonmax+1,sizeof(double));

		/* maxirrig: maximum amount of irrigation water from surface water sources (rivers or reservoirs) in m3/s per zone */
		maxirrigFromReservoirs    = (double**)calloc(tracer_count+1,sizeof(double*));
		for (z=0;z<=tracer_count;z++)
		maxirrigFromReservoirs[z] = (double*)calloc(zonmax+1,sizeof(double));

		zoneSW=NULL;
		maxirrig_temp = (double*)calloc(zonmax+1,sizeof(double));

		routing   = (routingclass**)calloc(anzrout+1,sizeof(routingclass*));

//cout << "rank: " << mpi_rank << " in Abflussmodellklasse vor initialize lakes" << endl;

		initialize_lakes(_lake_model_);

//cout << "rank: " << mpi_rank << " in Abflussmodellklasse nach initialize lakes" <<endl;

		// this grid will be required in rank 0 when using the lake model. Remember: we are here in rank 0 only, all other ranks skip the initialization of routing classes
		preconlakesindex = 1;
		while (gridinfo[preconlakesindex].identification != NULL)
		{
			preconlakesindex++;
			if (preconlakesindex == gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no more free grids in memory, (max. " << gridmaxanzahl << ")!" << "[end-entry]";
				throw(12);
			}
		}
		PRECONLAKES = matrix(1, nrows, 1, ncols, nrows);
		initialize_grid(preconlakesindex, _PrecOnLakes_, PRECONLAKES, NULL, 0, _PrecOnLakes_, 0, 1, nrows); /* because theta is not a grid but a stack, the adress-parameter for the grid is NULL */

		// search conrol file
		strcpy(tstr,"[");
		strcat(tstr,bezeichnung);
		strcat(tstr,"]");
		findeabschnitt(tstr);
		/* steuerdatei bis zur ersten Routingvorschrift vorspulen */
		for (int i=1;i<=7+numcoldescriptors;i++)
		{
			if (!feof(steuerdatei))
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
			else
			{
				(*WSstream) << "FATAL: error searching the second time for routing descriptions\n"
				"the first run was ok., but not the second"<<"[end-entry]";
				throw(38);
			}
		}

		/* Es wird aus der einzulesenden Datei der Abflussspenden und aus den simulierten Abfluessen ein Effizienz-
	Kriterium gebildet (R-square). Dazu muss angegeben werden, welche Spalte in der Spenden-Datei mit welchem Teilgebiet
	kommuniziert. */
		effizienz = (Tstatist**)calloc(21,sizeof(Tstatist*));
		log_effi  = (Tstatist**)calloc(21,sizeof(Tstatist*));
		for (int i=0;i<=20;i++)
		{
			effizienz[i] = (Tstatist*)calloc(zonmax+1,sizeof(Tstatist));
			log_effi[i]  = (Tstatist*)calloc(zonmax+1,sizeof(Tstatist));
		}

		for (int i=0;i<=20;i++)
		{
			for (int j=1;j<=zonmax;j++)
			{
				effizienz[i][j].sxdiff = effizienz[i][j].sx2diff = effizienz[i][j].ngem  =
				effizienz[i][j].sxgem  = effizienz[i][j].sx2gem  = effizienz[i][j].ndiff = 0.0;
				log_effi[i][j].sxdiff = log_effi[i][j].sx2diff = log_effi[i][j].ngem  =
				log_effi[i][j].sxgem  = log_effi[i][j].sx2gem  = log_effi[i][j].ndiff = 0.0;
			}
		}

		/* Absuchen der Steuerdatei nach routing-Vorschriften, "TG" wird als Schluesselwort benutzt */
		if (!feof(steuerdatei))
		{
			satz1=fgets(satz1,_stringsize_,steuerdatei);
		}
		else
		{
			(*WSstream) << "FATAL: rerror reading routing description \n"
			"2nd run, (1st run: o.k.)"<<"[end-entry]";
			throw(38);
		}
		/* Jeder Satz, in dem "TG" vorkommt, wird auf Oberlieger "OL" untersucht. */
		anzrout=0;
		tmax=0;
		strcpy(satz3,satz1);

		while (strstr(satz3,_TG_))
		{
			(*WSstream) << "reading routing description "<<(anzrout+1)<<"[end-entry]";
			/* Ein Satz mit einer Routing - Vorschrift ist gefunden */
			anzrout++;

			std::stringstream mystream;
			std::string ss;
			strcpy(satz1,satz3); /* nun wieder uebernehmen des Inhaltes von satz3 in Satz1, um diesen auszuwerten */
			do
			{
				if (satz1[strlen(satz1)-1]=='\n') satz1[strlen(satz1)-1]='\0';
				else
				{
					(*WSstream) << "FATAL: string is too long (max. "<<_stringsize_<<" chars):\n"<<satz1<<"[end-entry]";
					throw(40);
				}
			ss.append(satz1);

				if (!feof(steuerdatei))
				{
					satz1=fgets(satz1,_stringsize_,steuerdatei);
				}
				else
				{
					(*WSstream) << "FATAL: unexpected end of control file reached reading routing descriptions"<<"[end-entry]";
					throw(38);
				}
			}
			/* wenn satz1 "TG" enthaelt, so ist das der Beginn einer neuen
		Routingvorschrift, die spaeter ausgewertet wird. satz1 bleibt erhalten in satz3 */
			while (  (!strstr(satz1,_TG_)) && ( (strstr(satz1,_SUMTRIB_)) || (strstr(satz1,_OL_)) || (strstr(satz1,_AL_)) || (strstr(satz1,_ZL_)) || (strstr(satz1,_SP_)) ));
			mystream << ss;

			strcpy(satz3,satz1); /* Inhalt von satz1 wird in satz3 zwischengespeichert */

			strcpy(satz1,ss.c_str());
			strcpy(satz2,ss.c_str());

			/* Feststellen, wieviele Oberlieger _OL_, Zulaeufe _ZL_,
		Abgaben _AL_ und Speicher _SP_ im Gebiet sind */

			p__OL      = count_substr((char*)mystream.str().c_str(),_OL_);
			p__SUMTRIB = count_substr((char*)mystream.str().c_str(),_SUMTRIB_);
			p__et      = count_substr((char*)mystream.str().c_str(),"&");
			p__ZL      = count_substr((char*)mystream.str().c_str(),_ZL_);
			p__AL      = count_substr((char*)mystream.str().c_str(),_AL_);
			p__SP      = count_substr((char*)mystream.str().c_str(),_SP_);
			
			/* Nummer des Teilgebietes lesen */
			mystream.seekg(mystream.str().find(_TG_)+2);
			mystream >> tg;
			if ((tg > zonmaxwert)||(zarray[tg]<=0)||(zarray[tg]>zonmax))
			{
				(*WSstream) << "FATAL: there is no subbasin with code "<<tg<<"[end-entry]";
				throw(41);
			}

			/* TgAe Einzugsgebietsgroesse des Teilgebietes lesen */
			mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
			mystream >> p_AE;

			/* TgAe_rel relative Einzugsgebietsgroesse des Teilgebietes lesen
		(wenn ein OL nicht normal mitgeroutet, sondern als Zulauf aus
		einer externen Datei geroutet wird) */
			mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
			mystream >> p_AErel;

			/* eine neue Routing-Instanz einrichten */
			routing[anzrout] = new routingclass(
				retvalue,WSstream,steuername,LogLevel,InV,
				zonmax,p__OL,p__SUMTRIB,p__et,p__ZL,p__AL,
				p__SP,tg,p_AE,p_AErel,abfrel,zarray,
				tracer_count,tracer_info,
				gridinfo,zone,FirstLastCol ,zonmaxwert,
				ncols,nrows,mpistartrow,mpiendrow,xll,yll,cellsize,nodata,lake_active,readgrids,
				DefaultInputDirectory,DefaultOutputDirectory,
				v_mpiregioncnt,p_sendHalo,p_recvHalo,doSpinUp,vHRU, hour_offset, timestep,timestepfactor,minutea,taga,stundea,monata,jahra);
			if (retvalue!=0) throw retvalue;
			rtarget[2][zarray[tg]]=anzrout;
			k=0;
			/* Lesen der Routingparameter fuer die Oberlieger nacheinander */
			while (mystream.str().find(_OL_,mystream.tellg())!=std::string::npos)
			{
				for (int i=0;i<=9;i++) num_ST[i]=0;
				k++;

				/* Nummer (Zone) des Oberliegers */
				mystream.seekg(mystream.str().find(_OL_,mystream.tellg())+2);
				mystream >> p__OL;

				if ((p__OL > zonmaxwert)||(zarray[p__OL]<=0)||(zarray[p__OL]>zonmax))
				{
					(*WSstream) << "FATAL: there cant be a tributary with code "
					<<p__OL<<" for no subbasin (according to zones)\n"<<"[end-entry]";
					throw(41);
				}
				if (p__OL==tg)
				{
					(*WSstream) << "FATAL: subbasin "<<tg<<": code of tributary "<<k<<" ("<<p__OL<<") is equal to subbasin code!\n"
					"It is not allowed to rout the outflow of a subbasin back into \n"
					"the same subbasin. Use a reservoir instead.\n"<<"[end-entry]";
					throw(41);
				}
				/* Retentionskonstante Hauptbett dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kh;

				/* Retentionskonstante Vorland dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kv;

				/* Breite des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bh;

				/* Breite des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bv;

				/* Tiefe des Hauptbettes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Th;

				/* Manningbeiwert des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mh;

				/* Manningbeiwert des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mv;

				/* Gefaelle des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_I;

				/* Laenge des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_L;

				/* effektives Einzugsgebiet dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_AE;

				routing[anzrout]->SetRoutingParams(k,p__OL,num_ST,p_kh,p_kv,p_Bh,p_Bv,p_Th,p_Mh,p_Mv,p_I,p_L,p_AE,p_minspende,p_maxspende,p_intzahl,p_subint,gem_Zul);
				n=1;
				while (rtarget[0][n]>0) n++;
				rtarget[0][n]=zarray[p__OL];
				rtarget[1][zarray[p__OL]]=zarray[tg];

				//(*WSstream) << "DEBUG: read OL parameters for OL "<<p__OL<<"[end-entry]";
			}
			//(*WSstream) << "DEBUG: done: reading OL parameters "<<"[end-entry]";

			/* Lesen der Routingparameter fuer die SUMTRIBS nacheinander */
			mystream.clear();
			mystream.seekg(0,std::ios::beg); // reset the read position in order to start from the very beginning when looking for SUMTRIB
			while (mystream.str().find(_SUMTRIB_,mystream.tellg())!=std::string::npos)
			{
				for (int i=0;i<=9;i++) num_ST[i]=0; /* zunaechst nullsetzen der Teilgebietsnummern fuer alle potentiellen SUMTRIB-Oberlieger */
				k++; /* k zaehlt weiter von OLs! */

				/* Nummer (Zone) des 1. Oberliegers aller gemeinsam gerouteten Oberlieger */
				mystream.seekg(mystream.str().find(_SUMTRIB_,mystream.tellg())+strlen(_SUMTRIB_));
				mystream >> p__OL;
				if ((p__OL > zonmaxwert)||(zarray[p__OL]<=0)||(zarray[p__OL]>zonmax))
				{
					(*WSstream) << "FATAL: there cant be a tributary with code "<<p__OL<<" for no subbasin (according to zones)"<<"[end-entry]";
					throw(41);
				}
				if (p__OL==tg)
				{
					(*WSstream) << "FATAL: subbasin "<<tg<<": code of tributary "<<k<<" ("<<p__OL<<") is equal to subbasin code!\n"
					"It is not allowed to rout the outflow of a subbasin back into \n"
					"the same subbasin. Use a reservoir instead."<<"[end-entry]";
					throw(41);
				}
				/* nun lesen der restlichen Oberlieger-Nummern und Belegen des num_ST-arraya */
				
				n=1;
				while (rtarget[0][n]>0) n++;
				rtarget[0][n]=zarray[p__OL];
				rtarget[1][zarray[p__OL]]=zarray[tg];
				unsigned long pos1 = mystream.tellg();
				unsigned long pos2 = mystream.str().find_first_of(" (k",pos1);
				std::string SUMTRIBstring=mystream.str().substr(pos1,pos2-pos1);
				{
					int i=0;
					while (SUMTRIBstring.find("&")!=std::string::npos)
					{
						i++;
						SUMTRIBstring=SUMTRIBstring.substr(SUMTRIBstring.find("&")+1);
						sscanf(SUMTRIBstring.c_str(),"%d",&(num_ST[i]));
						n++; 
						rtarget[0][n]=zarray[num_ST[i]];
						rtarget[1][zarray[num_ST[i]]]=zarray[tg];
						num_ST[0]++;
					}
				}

				/* Retentionskonstante Hauptbett dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kh;

				/* Retentionskonstante Vorland dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kv;

				/* Breite des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bh;

				/* Breite des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bv;

				/* Tiefe des Hauptbettes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Th;

				/* Manningbeiwert des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mh;

				/* Manningbeiwert des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mv;

				/* Gefaelle des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_I;

				/* Laenge des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_L;

				/* effektives Einzugsgebiet dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_AE;

				routing[anzrout]->SetRoutingParams(k,p__OL,num_ST,p_kh,p_kv,p_Bh,p_Bv,p_Th,p_Mh,p_Mv,p_I,p_L,p_AE,p_minspende,p_maxspende,p_intzahl,p_subint,gem_Zul);
				//(*WSstream) << "DEBUG: read SUMTRIB parameters for SUMTRIB "<<p__OL<<"[end-entry]";
			}
			//(*WSstream) << "DEBUG: done: reading SUMTRIB parameters "<<"[end-entry]";

			/* nocheinmal Wiederherstellen des Satzes 1, um nun die Zulaufparameter lesen zu koennen */
			mystream.clear();
			mystream.seekg(0,std::ios::beg); // reset the read position in order to start from the very beginning when looking for SUMTRIB

			ZL_modus = (char*)calloc(20,sizeof(char));
			ZL_filename = (char*)calloc(300,sizeof(char));
			/* Lesen der Routingparameter fuer die Zuleitungen nacheinander */
			while (mystream.str().find(_ZL_,mystream.tellg())!=std::string::npos)
			{
				for (int i=0;i<=9;i++) num_ST[i]=0; /* nullsetzen der Teilgebietsnummern fuer alle potentiellen SUMTRIB-Oberlieger */
				k++;
				ZL_reihe=0;
				/* Nummer des Zulaufes */
				mystream.seekg(mystream.str().find(_ZL_,mystream.tellg())+strlen(_ZL_));
				mystream >> p__ZL;


				if (p__ZL > anz_ZL)
				{
					(*WSstream) << "FATAL: there is/are only "<<anz_ZL
					<<" inflow(s) for the entire basin,\n"
					" but the codes ranges from 1 to at least "<<p__ZL<<"! in "<<instanzname<<"[end-entry]";
					throw(42);
				}
				/* bei Zulaeufen erst den Modus und gegebenenfalls den Dateinamen und die Spalte der Datei lesen */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> ZL_modus;

				if (strcmp(ZL_modus,_constant_)==0)
				{
					ZL_art = 1;
					strcpy(ZL_filename,"");
					ZL_spalte=0;
					mystream >> ZL_wert;
				}
				else
				{
					if (strcmp(ZL_modus,_extern_)==0)
					{
						ZL_art=4;
						mystream >> ZL_filename >> ZL_reihe >> ZL_spalte;
						ZL_wert=0;
					}
					else
					{
						if (strcmp(ZL_modus,_intern_)==0)
						{
							ZL_art=3;
							strcpy(ZL_filename,"");
							ZL_spalte=0;
							ZL_reihe=0;
							ZL_wert=0;
						}
						else
						{
							(*WSstream) << "FATAL: the modus \""<<ZL_modus<<"\" is not supported for inflow "<<p__ZL<<"[end-entry]";
							throw(42);
						}
					}
				}
				routing[anzrout]->SetZulaufModus(k,p__ZL,ZL_art,ZL_filename,ZL_reihe,ZL_spalte,ZL_wert,jahra,monata,taga,stundea,minutea);

				/* Retentionskonstante Hauptbett dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kh;

				/* Retentionskonstante Vorland dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_kv;

				/* Breite des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bh;

				/* Breite des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Bv;

				/* Tiefe des Hauptbettes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Th;

				/* Manningbeiwert des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mh;

				/* Manningbeiwert des Vorlandes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_Mv;

				/* Gefaelle des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_I;

				/* Laenge des Gerinnes dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_L;

				/* effektives Einzugsgebiet dieses Oberliegers */
				mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
				mystream >> p_AE;

				routing[anzrout]->SetRoutingParams(k,p__ZL,num_ST,p_kh,p_kv,p_Bh,p_Bv,p_Th,p_Mh,p_Mv,p_I,p_L,p_AE,p_minspende,p_maxspende,p_intzahl,p_subint,gem_Zul,false);
				//(*WSstream) << "DEBUG: read ZL parameters for ZL "<<p__ZL<<"[end-entry]";
			}
			//(*WSstream) << "DEBUG: done: reading ZL parameters "<<"[end-entry]";

			/* nun lesen der Parameter fuer die Ableitungen */
			mystream.clear();
			mystream.seekg(0,std::ios::beg);
			while (mystream.str().find(_AL_,mystream.tellg())!=std::string::npos)
			{
				/* die Nummern der Ableitungen koennen mit denen der Zuleitungen identisch sein, dann
		gilt: wenn die Ableitung im Abarbeitungsschema vor der gleichlautenden Zuleitung
		auftaucht, dann wird noch im gleichen Zeitintervall die Ableitung beruecksichtigt,
		anderenfalls erst im folgenden Intervall (wenn also die Zuleitung vor der Ableitung
		definiert wurde. Wenn bei einer Zuleitung das Schluesselwort "intern"
		definiert wurde, die Werte also nicht constant sind oder aus einer Datei gelesen
		werden, muss es auch eine Ableitung mit diesem Schluesselwort geben, wobei dann
		natuerlich die Ableitung die gleiche Nummer wie die angezielte Zuleitung haben muss.
		Ableitung und Zuleitung muessen sich jeweils auf ein anderes Gebiet beziehen, duerfen also,
		wenn sie in einer einzigen Routingvorschrift gleichzeitig auftreten, nicht diegleiche
		Nummer haben, weil Ableitungen am AUslass, Zuleitungen aber beliebig im Gebiet erfolgen
		koennen! */
				k++;
				routing[anzrout]->SetAblaufModus(steuerdatei,k,mystream,p_subint,jahra,monata,taga,stundea,minutea);
				//(*WSstream) << "DEBUG: read AL parameters for AL "<<p__AL<<"[end-entry]";
			}
			//(*WSstream) << "DEBUG: done: reading AL parameters "<<"[end-entry]";

			/* nun lesen der Parameter fuer die Speicher */
			mystream.clear();
			mystream.seekg(0,std::ios::beg); // reset the read position in order to start from the very beginning when looking for SUMTRIB
			while (mystream.str().find(_SP_,mystream.tellg())!=std::string::npos)
			{
				k++;
				routing[anzrout]->SetSpeicherParams(steuerdatei,k,mystream,p_subint);
				//(*WSstream) << "DEBUG: read SP parameters for SP "<<k<<"[end-entry]";
			}
			//(*WSstream) << "DEBUG: done: reading SP parameters "<<"[end-entry]";

			if (routing[anzrout]->ALs > 0)
			{
				routing[anzrout]->LocalizeRefALs();
			}

			if (ZL_modus)    {free(ZL_modus); ZL_modus=NULL;	  }
			if (ZL_filename) {free(ZL_filename); ZL_filename=NULL;}
			//(*WSstream) << "DEBUG:done (TG: "<<routing[anzrout]->TG<<", OLs: "<<routing[anzrout]->OLs<<", ZLs: "<<routing[anzrout]->ZLs<<", ALs: "<<routing[anzrout]->ALs<<", SPs: "<<routing[anzrout]->SPs<<")"<<"[end-entry]";
		}
		satz3 = tsatz3; if (satz3) {free(satz3); satz3 = NULL;}
		satz1 = tsatz1; if (satz1) {free(satz1); satz1 = NULL;}
		satz2 = tsatz2; if (satz2) {free(satz2); satz2 = NULL;}

		satz1 = (char*)calloc(10*_stringsize_,sizeof(char)); tsatz1=satz1;


		mapReservoirLevels = new map<int, double**>();
		mapTributariesWithHydCon = new map<int, bool>();

		// this is called after initializing all routing elements because the control files file pointer will be 
		// changed. Same is true for reading the abstractionrules for AL-elements
		for (int i = 1; i <= anzrout; i++)
		{
			routing[i]->ReadReservoirRule(steuerdatei);
			routing[i]->FillMapReservoirLevels(mapReservoirLevels);
			routing[i]->FillMapTributariesWithHydCon(mapTributariesWithHydCon);
		}
		// cannot be included in the above loop because all elements of the mapReservoirLevels must be 
		// known before the first actual outflow of a reservoir can be calculated, that is 
		// hydraulically connected to another reservoir
		for (int i = 1; i <= anzrout; i++)
		{
			routing[i]->SetFirstAbstractionValue();
			routing[i]->ReadAbstractionRule(steuerdatei);
		}



		/* suchen der maximalen Anzahl an zu verschiebenden 1/subint-tel Zeitintervallen */
		tmax = 0;
		for (int i=1;i<=anzrout;i++)
		tmax=routing[i]->getmaxdelta((int)tmax);

		/* Gewichte neu berechnen (fuer Summe Teilgebiete=1, daraus wird jeweils das
	Gewicht eines Pegeleinzugsgebietes am Gesamtgebiet berechnet) */
		for (int i=1;i<=anzrout;i++)
		{
			routing[i]->SetWeights();
			routing[i]->Print_W_Q_Table();
			routing[i]->SetMaxIrrigAdress(maxirrig,maxirrigFromReservoirs);
		}
		Qstat = new routstatistic(retvalue,WSstream,LastResults,maxLastResults,NULL,1,1,&zone,abfrel[0],(int)zonmax,(int)Qstatcode,
		(int)5,(char*)Qstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));

		strcpy(tstr,"[output_interval]");
		findeabschnitt(tstr);
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		fscanf(steuerdatei,"%d",&runoff_measure); 
		tstr=fgets(tstr,_stringsize_,steuerdatei);

		if (retvalue!=0) throw retvalue;
#if defined _ENGLISH_
		//    lastresultsindex=Qstat->schreibzonenwerte("routed discharge [mm per subbasin] ",instanzname,LastResults,maxLastResults);
		if ( runoff_measure == 0)  {
			lastresultsindex=Qstat->schreibzonenwerte("routed discharge [mm per subbasin] ",instanzname,LastResults,maxLastResults);
		} else if ( runoff_measure == 1) {
			lastresultsindex=Qstat->schreibzonenwerte("routed discharge [ m3 per sec ] ",instanzname,LastResults,maxLastResults);
		} else {
			lastresultsindex=Qstat->schreibzonenwerte("routed discharge (unknown measure)",instanzname,LastResults,maxLastResults);
		}
#else
		lastresultsindex=Qstat->schreibzonenwerte("gerouteter Abfluss je (Gesamt)Teilgebiet ",instanzname,LastResults,maxLastResults);
#endif

		/* Einrichten eines Speicherarrays, welches jeden Fluss fuer maxinal tmax Intervalle zwischenspeichert
	(wegen der Translation, wobei jedes Zeitintervall in <subint> Unterintervalle aufgeteilt ist!).
	im fluss[1..subint,x] steht fuer den x. Fluss der geroutete Abfluss im aktuellen Intervall
	In fluss[1+subint..2*subint,x] steht dann der Abfluss aus dem naechsten Intervall usw. Die Werte werden am Beginn des Routings
	aus Qges uebernommen und soweit, wie es ihrer Translationszeit (aus wqtab.delta[Q]) entspricht hinter dem aktuellen
	Wert eingeordnet.
	Die (Zwischengebiets)-Abfluesse in Routing-Elemente werden auch erst uebernommen,
	jedoch dann durch das Routing-Ergebnisse ueberschrieben (das sich aus der Ueberlagerung der gerouteten Oberlieger
	abfluesse und dem Zwischengebietsabfluss ergibt, weshalb die Uebernahme zunaechst auch noetig ist. */

		/* wenn bei einer Zeitintervallunterteilung in 12 Teilintervalle die laengste moegliche Translationszeit etwa 10 h ist,
	wobei diese Strecke in 5 Teilstrecken a max. 2 Stunden aufgeteilt ist,so werden fuer alle Gerinnestrecken demnach
	5*(2*12+12)+12 = 192 Elemente reserviert. Der Platzbedarf ist relativ gering:
	bei 100 Teilgebieten bzw. 100 Routingstrecken mit einer Teilung in jeweils 10 Strecken und
	12 zeitintervalle bei max. 2h Fliesszeit pro Teilstrecke sind das:
	100*10*(2*12+12) = 36000 Werte, entspricht
	etwa 140 Kbyte, bei 200 Gebieten und max. 40 h Laufzeit werden es 144000 Elemente oder
	562 KByte, was im Vergleich zur Gridgroesse (Thur, 7000 Punkte: etwa 58 Kbyte,
	Schweiz, 250-m-Raster: 4 Mbyte!) doch noch relativ wenig ist.
	Dass noch 12 Werte pro Strecke mehr reserviert werden liegt daran, dass sie das aktuelle Intervall
	repraesentieren. Fuer jede Teilstrecke werden soviele Elemente benoetigt, wie die maximale Verschiebung
	zuzueglich der Laenge eines Zeitintervalls in Stunden multipliziert mit der Anzahl der Teilzeitintervalle ergibt,
	also bei max 2h Fliesszeit in einer Teilstrecke und 12 Intervallen a 5 min. in einem Stundenintervall des Modells
	sind das 3*12 = 36 Elemente, eine gesamte Routingstrecke benoetigt dann die Anzahl der Elemente einer Teilstrecke
	mal die Anzahl der Unterteilungen zuzueglich der Anzahl der Teilzeitintervalle, also bei einer Einteilung in
	sieben Teilstrecken 7*36 + 12 = 264 Elemente (=1056 Byte) */

		fluss    = (double***)calloc(tracer_count+1,sizeof(double**));
		for (z=0;z<=tracer_count;z++)
		fluss[z] = (double**)calloc(1+zonmax+anz_ZL,sizeof(double*));

		/* da die Abfluesse in 12-tel Intervalle eingeteilt werden, wird nun noch eine Matrix aus 12 Feldern pro Gebiet
	eingerichtet, in die die Original- und die gerouteten Abfluesse  (wie frueher Qges) kopiert werden */
		Q12ges=(double***)calloc(tracer_count+1,sizeof(double**));
		for (z=0;z<=tracer_count;z++)
		{
			Q12ges[z] = (double**)calloc(p_subint+1,sizeof(double*));
			for (k=0;k<=p_subint;k++)
			Q12ges[z][k]=(double*)calloc(zonmax+anz_ZL+1,sizeof(double));
		}
		//    Q12ges[z] = matrix(0,p_subint,0,zonmax+anz_ZL);

		// nun noch Oeffnen des Spendenfiles und Vorspulen bis zum benoetigten Datum 
		// -> Spenden werden in einen Spendenbuffer eingelesen, um Lesezugriffe sicherer zu machen
		SpendenBuffer=new FileBuffer(spendenfilename,0,1000000);
		if (SpendenBuffer->GetStatus()==FB_STATE_FILENOTFOUND)
		{
			(*WSstream) << "INFO: error opening "<<spendenfilename
			<<"\nRouting will be done without comparison to observed discharge data!"<<"[end-entry]";
			spendenende=10;
			Spende=NULL;
			gem_Zul = 0; /* wenn keine Spendendatei angegeben werden wurde, koennen auch keine gemessenen Zulaufganglinien genutzt werden */
		} else if (SpendenBuffer->GetStatus()==FB_STATE_FILENOTGOOD) {
			(*WSstream) <<"FATAL: error analysing or reading "<<spendenfilename
			<<" at row "<<SpendenBuffer->GetBufpos() 
			<<" or before \n(e.g. due to wrong formatting or wrong number of columns or start date not found)" 
			<<"[end-entry]";
			throw 203;
		}
		else /* Nur dann, wenn auch eine Spendendatei gefunden wurde, werden die folgenden Sequenzen noetig sein */
		{

			(*WSstream) << "INFO: looking for starting date in " << spendenfilename	<< "[end-entry]";
			/* Spende wird als Feld mit dem aktuellen und den letzten und naechsten 10 Datensaetzen im Speicher gehalten
		Beim Einlesen des naechsten Wertes werden jeweils die Zeiger verdreht */
			spendenanzahl=SpendenBuffer->GetValcnt();
			// if there are less columns in the spendenfile than there are zones, or subcatchments to rout, that is, 
			// the elements in Spalte having a too large value have to be set to 0 in order to avoid invalid array access later
			for (int i=1;i<=zonmax;i++)
			{
				if (Spalte[i]>spendenanzahl)
				{
					Spalte[i]=0;
				}
			}
			Spende = (double**)calloc(21,sizeof(double*));
			for (int i=0;i<=20;i++)
			Spende[i] = (double*)calloc(spendenanzahl+1,sizeof(double));

			spendenende=0;
			/* Vorspulen der Spendendatei */
			Ys=0;
			int oldDs = 0;
			SpendenBuffer->SetBufPos(-1); // auf eine Position vor Beginn setzen
			while ((Ys!=jahra)||(Ms!=monata)||(Ds!=taga)||(Hs!=stundea))
			{
				SpendenBuffer->NextRecord();
				SpendenBuffer->GetActdate(Ys,Ms,Ds,Hs,ms);
				if (Ds < oldDs)
				{
					(*WSstream) << ".";
					cout.flush();
				}
				oldDs = Ds;
				if (Ys<1000) Ys+=1900;
				if (SpendenBuffer->GetStatus()==FB_STATE_LASTRECORDREAD) 
				{
					spendenende++;
					break;
				}

			}
			(*WSstream) << "[end-entry]";

			/*    satz=(char*)calloc(_stringsize_+1,sizeof(char)); satz wird schon in modelobjekt::modelobjekt belegt */
			/* Da die Spenden als Feld je 10 Werte im Voraus und im Nachhinein gehalten werden, werden die naechsten
		10 Saetze noch gelesen */
			for (int i=10;i<=20;i++)
			{
				for (int j=1;j<=spendenanzahl;j++)
				{
					Spende[i][j]=SpendenBuffer->GetActValue(j-1);
				}
				if (SpendenBuffer->GetStatus()==FB_STATE_LASTRECORDREAD) 
				spendenende++;
				SpendenBuffer->NextRecord();

				if (Ys<1000) Ys+=1900;
			}
			/* fuer den Initialisierungsabfluss werden die Werte aus Spende[10] in die Elemente [0..9] kopiert */
			for (int i=0;i<10;i++)
			{
				for (int j=1;j<=spendenanzahl;j++)
				{
					Spende[i][j]=Spende[10][j];
				}
			}
		}
		satz1 = tsatz1; if (satz1) { free(satz1); satz1 = NULL; }
		time=0;
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	// catching all exceptions in order to allow the object to be constructed correct. This is required 
	// in order to allow a dextructor call in case of an error during construction
	{
		retvalue=103;
	}
}

/* das setzen der Adressen des Abflussarrays erfolgt erst nach der
Initialisierung der Abflussmodellklassen-Instanz, weil in verschiedenen
Programmen auf verschiendene Weise die ABfluesse bereitgestellt werden */
void abflussmodellklasse::SetQgesAdresse(double **adresse)
{
	int z;
	if (do_it==0) return;
	
	if (routingmode && adresse == NULL)
	{ // 
		// ToDo: when in routing only mode, the Qges values must be read from a file, so we will open the file here and set the data structure of the values read from file as address here
		// this should be set of separate methods rather than an implemenation in this method
		routingmode_initialize();
	}
	else
	{
		Qges_uebernommen = adresse;
	}

	Qges = (double**)calloc(tracer_count+1,sizeof(double*));
	for (z = 0; z <= tracer_count; z++)
	{
		Qges[z] = (double*)calloc(zonmax + 1 + anz_ZL + 1 + anz_AL + 1, sizeof(double));
	}
	/* Die initial-Abfluesse werden in das Array aller gerouteten Abfluesse geschrieben,
	es wird anfangs ein stationaeres Fliessen mit dem Initialabfluss vorausgesetzt */
	if (Qges_uebernommen==NULL)
	{
		(*WSstream) << "FATAL: No valid discharges for discharge routing!\n"
		"You have to run the model including at least precipitation interpolation and soil model\n"
		"if discharges are wished."<<"[end-entry]";
		throw(43);
	}
	return ;
}

void abflussmodellklasse::StorageFileError(const char *message,double c1,double c2)
{
	(*WSstream) << "FATAL: error in storage file content reading "<<message<<":\n"
	"value in actual control is:   "<<c1<<"\n"
	"but value in storage file is: "<<c2<<"\n"
	"WaSiM throws..."<<"[end-entry]";
	throw(72);
}

void abflussmodellklasse::ReadStorageFile(char *filename,int flag)
{
	char *satz;
	int syntaxVersion = 2; // default syntax expected is for version 2 of the storage file

	storagefilename=(char*)calloc(_stringsize_,sizeof(char));
	if (filename) strcpy(storagefilename,filename); else storagefilename=NULL;

	char *completedname = new char[_stringsize_+1];
	sprintf(completedname,"%s",storagefilename);
	ExtendFileName(DefaultInputDirectory,completedname);
	//  delete[] completedname;

	if ((filename)&&(flag!=2)&&(flag!=4))
	{
		// if readgrid in Unsatzonmodel is 2 or 4, then only the groundwater table is read in, but no storage file.
		if (flag) // read the storage file for unsatzonmodel.readgrid==1 or 3
		{
			satz=(char*)calloc(_stringsize_,sizeof(char));
			if ((storagefile = fopen(completedname,"rt"))==NULL)
			{
				(*WSstream) << "FATAL: error opening storage content file "<<completedname<<"[end-entry]";
				throw(64);
			}

			char cmpString[100];
			char* tmp = satz;
			double timestepstor = this->timestep;
			sprintf(cmpString, "[routing_model_v2]");

			do
			{ 
				satz=fgets(satz,_stringsize_,storagefile); 
				if (feof(storagefile))
				{
					satz = tmp; // because fgets sets satz to NULL when feof has been reached, we reset the pointer to it's original value
					(*WSstream) << "INFORMATION: couldn't find section [unsatzon_model_v2], looking now for [unsatzon_model] in " << completedname << "[end-entry]";
					syntaxVersion = 1;
					fseek(storagefile, 0, 0); // reset storagefile
					sprintf(cmpString, "[routing_model]");
					do
					{
						satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
						if (feof(storagefile))
						{
							(*WSstream) << "FATAL: error locking for the section [routing_model_v2] or [routing_model] in " << completedname << "[end-entry]";
							throw(71);
						}
					} while (!strstr(satz, cmpString) && !feof(storagefile));
				}
			} 
			while (!strstr(satz, cmpString) && !feof(storagefile));

			if (syntaxVersion > 1)
			{
				fscanf(storagefile, "%lf", &timestepstor);
				satz = fgets(satz, _stringsize_, storagefile);
			}

			for (int i=1;i<=anzrout;i++)
			{
				routing[i]->ReadStorageFile(storagefile,completedname,timestepstor,syntaxVersion);
			}

			fclose(storagefile);
			if (satz) { free(satz); satz = NULL; }
		}

	}
	delete[] completedname;
}

void abflussmodellklasse::WriteStorageFile()
{  /* writing contents of overland flow and interflow single linear reservoirs to a file named storage.ftz */
	int i;
	FILE *tempfile;
	if (storagefilename)
	{
		char *completedname = new char[_stringsize_+1];
		sprintf(completedname,"%s",storagefilename);
		ExtendFileName(DefaultOutputDirectory,completedname);
		//  delete[] completedname;

		if ((storagefile = fopen(completedname,"wt"))!=NULL)
		{
//			fprintf(storagefile, "[" _Abflussmodell_ "]\n");
			fprintf(storagefile, "[routing_model_v2]\n");
			fprintf(storagefile, "%.3lf (timestep in minutes)\n", InV);
			/* hier muss nun der Inhalt aller Flussstrecken rein */
			for (i=1;i<=anzrout;i++)
			{
				routing[i]->WriteStorageFile(storagefile);
			}
			fprintf(storagefile,"\n");
			fclose(storagefile);
		}
		else
		{
			(*WSstream) << "WARNING: couldn't open file \""<<completedname<<"\" \n"
			"(model version doesn't support this feature or write access error)"<<"[end-entry]";
		}
		delete[] completedname;
	}
	else /* trying to write a default file */
	{
		char *completedname = new char[_stringsize_+1];
		sprintf(completedname,"storage.ftz");
		ExtendFileName(DefaultOutputDirectory,completedname);
		//  delete[] completedname;
		if ((tempfile = fopen(completedname,"wt"))!=NULL)
		{
//			fprintf(tempfile, "[" _Abflussmodell_ "]\n");
			fprintf(tempfile, "[routing_model_v2]\n");
			/* hier muss nun der Inhalt aller Flussstrecken rein */
			for (i=1;i<=anzrout;i++)
			{
				routing[i]->WriteStorageFile(tempfile);
			}
			fprintf(tempfile,"\n");
			fclose(tempfile);
		}
		else
		{
			(*WSstream) << "WARNING: couldn't open file "<<completedname<<" \n"
			"(no filename was given or model version \n"
			"doesn't support this feature or the current user has no write access)"<<"[end-entry]";
		}
	}
	/* -----------ENDE STORAGEFILE -------------- */
}

void abflussmodellklasse::SetRoutingFluxAdresses()
{
	int i,j,n,max,z,k;
	double **Qtemp;
	if (do_it==0) return;

	/* formales Uebernehmen der Abfluesse in das erweiterte Array Qges aus Qges_uebernommen (Qges_uebernommen darf nocht direkt
	geroutet werden, weil im Bodenmodell Qges im naechsten Zeitschritt noch gebraucht werden kann) */
	for (z = 0; z <= tracer_count; z++)
	{
		for (i = 1; i <= zonmax; i++)
		{
			Qges[z][i] = Qges_uebernommen[z][i];
		}
	}

	for (z = 0; z <= tracer_count; z++)
	{
		for (i = zonmax + 1; i <= zonmax + anz_ZL; i++)
		{
			Qges[z][i] = 0;
		}
	}

	for (i = 1; i <= anzrout; i++)
	{
		routing[i]->SetFluxAdresses(Qges, Q12ges, fluss, Spende, Spalte);
	}
	/* darin werden auch die Qges-Werte ins fluss-array uebernommen */

	for (z = 0; z <= tracer_count; z++)
	{
		for (i = 1; i <= p_subint; i++)
		{
			for (j = 1; j <= zonmax + anz_ZL; j++)
			{
				Q12ges[z][i][j] = Qges[z][j];
			}
		}
	}

	// check, if all internal ZL/AL-combinations are correctly connected
	(*WSstream) << "checking connection of internal inflows and abstractions" << "[end-entry]";
	bool connected;
	for (int i=1;i<=anzrout;i++)
	{
		if (routing[i]->ZLs>0)
		{ // do this check for routing elements with inflows only
			for (int zl=1;zl<=routing[i]->ZLs;zl++)
			{ // step through all inflows
				if (routing[i]->ZL_art[zl]==3)
				{  // do connection check for internal inflows only
					connected = false;
					int   ZL_code     = routing[i]->ZL[zl];
					double ZL_TgAe     = routing[i]->TgAe;
					double ZL_TgAe_rel = routing[i]->TgAe_rel;
					for (int j=1;j<=anzrout;j++)
					{
						if (routing[j]->ALs>0)
						{
							for (int al=1;al<=routing[j]->ALs;al++)
							{
								if ((
											(routing[j]->Abstraction[al].AL_art==3) || // _intern_
											(routing[j]->Abstraction[al].AL_art==6) || // _intern_with_rule
											(routing[j]->Abstraction[al].AL_art==8)    // _intern_with_rule_from_reservoir
											)&&
										(routing[j]->Abstraction[al].AL_code==ZL_code))
								{
									if (connected)
									{
										(*WSstream) << "FATAL: \n ---> Error in connecting internal abstractions with internal inflows:\n" <<
										"in routing description "<<j<<" for TG "<<routing[j]->TG<<":\n"<<
										"Inflow "<<routing[i]->ZL[zl]<<" is already connected to an internal abstraction.\n"<<
										"Each abstaction (AL-element) must be conencted to exactly one internal inflow(ZL-element)\n"<<"[end-entry]";
										throw(300);
									}
									if (routing[i]->TG==routing[j]->TG)
									{
										(*WSstream) << "FATAL: \n ---> Error in connecting internal abstractions with internal inflows:\n" <<
										"Abstraction AL "<<routing[j]->Abstraction[al].AL_code<<" from routing description "<<j<<
										" and inflow ZL "<<routing[i]->ZL[zl]<<" are both defined for the same routing-subbasin TG "<<routing[j]->TG<<".\n"<<
										"this is not allowed\n"<<"[end-entry]";
										throw(301);
									}
									connected=true;
//									routing[j]->Abstraction[al].SetTargetData(zarray[routing[i]->TG], ZL_TgAe, ZL_TgAe_rel, zonmax + routing[j]->AL[al], routing[j]->TgAe);
									routing[j]->Abstraction[al].SetTargetData(zarray[routing[i]->TG], ZL_TgAe, ZL_TgAe_rel, routing[i]->OLs>0?zarray[routing[i]->OL[1]]:0, routing[i]->TgAe);
								}
							}
						}
					}
					if (!connected)
					{
						(*WSstream) << "FATAL: \n ---> Error in connecting internal abstractions with internal inflows:\n" <<
						"in routing description "<<i<<" for TG "<<routing[i]->TG<<":\n"<<
						"Internal Inflow ZL "<<routing[i]->ZL[zl]<<" has no conenction to an internal abstraction.\n"<<
						"this is not allowed\n"<<"[end-entry]";
						throw(302);
					}
				}
			}
		}
	}
	// now check, if there are abstractions which still are not connected 
	// (if there are more internal abstractions than internal inflows, this may be the case)
	for (int i=1;i<=anzrout;i++)
	{
		if (routing[i]->ALs>0)
		{ // do this check for routing elements with inflows only
			for (int al=1;al<=routing[i]->ALs;al++)
			{ // step through all inflows
				if ((
							(routing[i]->Abstraction[al].AL_art==3)|| // _intern_
							(routing[i]->Abstraction[al].AL_art==6)|| // _intern_with_rule_
							(routing[i]->Abstraction[al].AL_art==8)   // _intern_with_rule_from_reservoir_
							)&&
						(routing[i]->Abstraction[al].ZL_target<=0))
				{  // do connection check for internal inflows only
					(*WSstream) << "FATAL: \n ---> Error in connecting internal abstractions with internal inflows:\n" <<
					"in routing description "<<i<<" for TG "<<routing[i]->TG<<":\n"<<
					"Internal abstraction AL "<<routing[i]->Abstraction[al].AL_code<<" is not connected to any internal inflow.\n"<<"[end-entry]";
					throw(303);
				}
			}
		}
	}

	/* nun muessen fuer alle Routinginstanzen die Initialwerte
	fuer wqtab[*].last_index[*] gesetzt werden (pro Teilstrecke, wo noetig) */
	for (i = 1; i <= anzrout; i++)
	{
		routing[i]->set_initial_last_index();
	}
	/* Ermitteln der laengsten Translationszeit
	jeder Routingvorschrift und Addieren dieser Zeiten.
	Das Ergebnis gibt an, wie oft das Routing vor Beginn
	des Modells mit konstanten Abfluessen schon ausgefuehrt
	werden muss, um einen stabilen Anfangsabfluss zu erhalten */
	n=0;
	for (i=1;i<=anzrout;i++)
	{
		max=0;
		for (j=1;j<=routing[i]->OLs+routing[i]->SUMTRIBs+routing[i]->ZLs;j++)
		{
			if ((int)(routing[i]->wqtab[j].t[1]) > max) max = (int)(routing[i]->wqtab[j].t[1]);
		}
		n+=max;
	}

	/* Anzahl der Initialisierungsroutingsequenzen von Stunden auf Zeitintervalle umrechnen */
	n = (int)(ceil(n/InV*60.0));
	if (n<1) n=1;

	if (readgrids==0)
	{
		(*WSstream) << "routing the initialization discharges ("<<n<<" times)"<<"[end-entry]";
		Qtemp = matrix(0,tracer_count+1,0,zonmax+anz_ZL+1);
		for (int j=1;j<=n;j++)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				for (int i=1;i<=zonmax+anz_ZL;i++) 
				{
					Qtemp[z][i]=Qges[z][i];
				}
			}

			for (int i=1;i<=p_subint;i++)
			{
				for (int k=1;k<=zonmax+anz_ZL;k++)
				{
					for (int z=0;z<=tracer_count;z++)
					{
						Q12ges[z][i][k]=Qges[z][k]; /* die Zulaeufe werden erst in den Routingelementen in Qges und Q12ges uebernommen */
					}
				}
			}


			for (int i=1;i<=anzrout;i++) 
			{
				routing[i]->Routing(0,monata,taga,stundea,0,0,1);
			}
			
			(*WSstream) << ".";
			
			for (int z=0;z<=tracer_count;z++) 
				for (int i=1;i<=zonmax;i++) 
			//     for (i=1;i<=zonmax+anz_ZL;i++) 
					Qges[z][i]=Qtemp[z][i];
		}
		(*WSstream) << "[end-entry]";
		free_matrix(Qtemp,0,tracer_count+1,0,zonmax+anz_ZL+1); Qtemp = NULL;
	}
	return;
}

void abflussmodellklasse::SetZoneSWAdress(statrecord *p_zoneSW)
{
	int i;
	if (do_it==0) return;
	else
	{
		zoneSW=p_zoneSW; /* if no irrigation is applied (1st parameter in control file for irrigation is 0), p_zoneSW is set to NULL */
		for (i=1;i<=anzrout;i++)
		routing[i]->SetZoneSWAdress(zoneSW); /* NULL, if no irrigation is applied, otherwise != 0  */
	}
}

double **abflussmodellklasse::GetMaxIrrigAdress()
{
	if (do_it==0) return NULL;
	else return maxirrig;
}

double **abflussmodellklasse::GetMaxIrrigFromReservoirsAddress()
{
	if (do_it==0) return NULL;
	else return maxirrigFromReservoirs;
}

int **abflussmodellklasse::GetRtargetAdress()
{
	if (do_it==0) 
	return NULL;
	else
	return rtarget;
}

double **abflussmodellklasse::GetAbfrelAdress()
{
	if (do_it==0)
	return NULL;
	else
	return abfrel;
}


void abflussmodellklasse::SetRunoffMeasure(int p_runoff_measure)
{
	if (do_it==0) return;
	else
	runoff_measure=p_runoff_measure;
}


void abflussmodellklasse::run(long me,long se, long tg, long mt, long jr, int actsubtimestep,int numsubtimesteps)
{
	int i,j,z,routindex,actindex;
	double temp,logwert;
	double *zeiger;
	double irr_diff;

	jahr=jr;
	monat=mt;
	tag=tg;
	stunde=se;

	if (routingmode)
	{
		if (lastday != tg)
		{
			sprintf(additionalinfo, "routing only for day %.04ld-%.02ld-%.02ld",jr,mt,tg);
			lastday = tg;
		}
	}
	else {
		sprintf(additionalinfo, "");
	}

	// some grids must be gathered if we are in MPI mode
	// this will be called regardless do_it, since all mpi_ranks > 0 have do_it set to 0.
	// but since lake_active is only > 0 if the instance should run in mpi_rank 0 it is save to do the grid gathering here
	if ((this->lake_active)&&(do_it_mpi0>0))
	{
		MPIGatherGridWithMask(pondindex,lakeindex,lakeCellsPerMPIRank);
		
		if (etrindex > 0)
		{
			MPIGatherGridWithMask(etrindex, lakeindex, lakeCellsPerMPIRank);
		}
		
		if (eiindex > 0)
		{
			MPIGatherGridWithMask(eiindex, lakeindex, lakeCellsPerMPIRank);
		}
		if (quindex>0)
		{
			MPIGatherGridWithMask(quindex,lakeindex,lakeCellsPerMPIRank);
		}

// QB was already gathered after unsatzonmodel
		if (qbindex>0)
		{
			MPIGatherGridWithMask(qbindex,lakeindex,lakeCellsPerMPIRank);
		}

		MPIGatherGridWithMask(pindex,lakeindex,lakeCellsPerMPIRank,PRECONLAKES); // PRECONLAKES is a local grid known to runoff routing only which is the target for gathering the prec grids from all other ranks

		// todo: hier noch alle conc_stor_unsatzon-n gathern
		for (int z=1;z<=tracer_count;z++)
		{
			MPIGatherGridWithMask(tracerindex[z],lakeindex,lakeCellsPerMPIRank);
		}

//		if (mpi_rank==0)
//		{
//			char name[300];
//			sprintf(name,"MPI_Pond_%04d-%02d-%02d_%02d.grd",jahr,monat,tag,stunde);
//			printgridfile(name,gridinfo[pondindex].pointer,true,true);
//			sprintf(name,"MPI_ETR__%04d-%02d-%02d_%02d.grd",jahr,monat,tag,stunde);
//			printgridfile(name,gridinfo[etrindex].pointer,true,true);
//			sprintf(name,"MPI_QU___%04d-%02d-%02d_%02d.grd",jahr,monat,tag,stunde);
//			printgridfile(name,gridinfo[quindex].pointer,true,true);
//			sprintf(name,"MPI_QB___%04d-%02d-%02d_%02d.grd",jahr,monat,tag,stunde);
//			printgridfile(name,gridinfo[qbindex].pointer,true,true);
//			sprintf(name,"MPI_PREC_%04d-%02d-%02d_%02d.grd",jahr,monat,tag,stunde);
//			printgridfile(name,gridinfo[pindex].pointer,true,true);
//		}
	}

	if (mpi_rank>0) 
	{
		// MPIScatterGridWithMask will be called for rank 0 at the end of the ::run method in order to synchronize again with all ranks
		MPIScatterGridWithMask(pondindex,lakeindex,lakeCellsPerMPIRank); 
//		cout << "rank " << mpi_rank << " skipping routing" << endl;
		return; // this module is currently not able to run on a cluster. Only the instance running mpi rank 0 will do channel routing in order to keep the model running at all
	}
	if (do_it == 0) return;

	/* Ausgabe einer Nachricht */
#ifdef _OPENMP
	if (!omp_in_parallel())
	{
		teilmodellklasse::run(me,se,tg,mt,jr,NULL,additionalinfo);
	}
#else
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,additionalinfo);
#endif

	/* set maxirrig to 0. This is an indicator within the routing-method, 
	that the tributary is a 1st-level-river (irrigation water has to 
	be abstracted from the inflow, not from routed fluxes) */
	if (zoneSW)
	{
		for (i=1;i<=zonmax;i++)
		{
			for (z=0;z<=tracer_count;z++)
			{
				maxirrig[z][i]=0;
				// do not set maxirrigFromReservoirs to 0 because it contains the 
				// amount of water, which must be subtracted 
				// from the reservoir content in the routing model right now (in routingclass::Routing()
				// maxirrigFromReservoirs[z][i]=0;
			}
			if (zoneSW[i].RelCellSum>0)
			zoneSW[i].sx=(zoneSW[i].sx/zoneSW[i].RelCellSum); /* mean irrigation demand for each zone from surface water in mm */
			else
			zoneSW[i].sx=0.0;
		}
		for (i=1;i<=zonmax;i++)
		{
			zoneSW[i].sx-=maxirrigFromReservoirs[0][i];
			if (zoneSW[i].sx<0)
			{
				zoneSW[i].sx=0;
			}
		}
	}

//	Drehen der Zeiger muss erst dann erfolgen, wenn Effizienzkriterien upgedated wurden, also nach dem Routing und nicht davor
//	if (SpendenBuffer->GetStatus()!=FB_STATE_FILENOTFOUND)
//	{
//		if (spendenende>=10) 
//		(*WSstream) << "INFO: observed runoff: actual date "<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde
//		<<"h not found (end of discharge file at "<<Ds<<"."<<Ms<<"."<<Ys<<", "<<Hs<<"h)"<<"[end-entry]";
//		zeiger=Spende[0];
//		for (i=0;i<20;i++) Spende[i]=Spende[i+1];
//		Spende[20]=zeiger;
//
//		/* Einlesen der nächsten Spenden */
//		if (spendenende<10)
//		{
//			for (i=1;i<=spendenanzahl;i++)
//			{
//				Spende[20][i]=SpendenBuffer->GetActValue(i-1);
//			}
//			SpendenBuffer->NextRecord();
//			if (SpendenBuffer->GetStatus()==FB_STATE_LASTRECORDREAD)
//			{
//				(*WSstream) << "INFO: end of discharge file reached!"<<"[end-entry]";
//				spendenende++;
//			}
//			else
//			{
//				SpendenBuffer->GetActdate(Ys,Ms,Ds,Hs);
//				if (Ys<1000) Ys+=1900;
//			}
//		}
//		else
//		{
//			for (i=1;i<=spendenanzahl;i++) Spende[20][i]=0.0;
//		}
//	}

	/* Uebernehmen der Abflussdaten aus dem Qges-Array aus dem Bodenmodell
	bzw. aus der eigens (im Programm flowrout) eingelesenen Datei mit
	Gebietsgesamtabfluessen) und das eigentliche Routing */

	/* eigentlich ist Qges eine Duplizierung von Qges_uebernommen,
	aber es hat noch <ZLs> zusaetzliche Felder, die spaeter gefuellt werden.
	Um das Programm uebersichtlich zu halten, sollen Zulaeufe genauso gereoutet werden wie
	normale Oberlieger, dazu muss aber Qges so gross sein, dass es alle Oberlieger (Teilgebiete)
	und Externen Zufluesse (ZLs) aufnehmen kann */
	for (j = 1; j <= zonmax; j++)
	{
		for (z = 0; z <= tracer_count; z++)
		{
			Qges[z][j] = (z==0)? Qges_uebernommen[z][j] / timestepfactor : Qges_uebernommen[z][j];
		}
	}
	for (i=1;i<=anzrout;i++)
	{
		routing[i]->ReadZuleitung(Qges,jahr,monat,tag,stunde,minute);
	}

	for (i = 1; i <= p_subint; i++)
	{
		for (j = 1; j <= zonmax + anz_ZL; j++)
		{
			for (z = 0; z <= tracer_count; z++)
			{
				Q12ges[z][i][j] = Qges[z][j]; /* die Zulaeufe werden erst in den Routingelementen in Qges und Q12ges uebernommen */
			}
		}
	}

//#ifdef _DEBUG
//	ofstream rivflow;
//	rivflow.open("runoff_routing_riv_flow.dat",std::ios::app);
//	rivflow << "step " << jr << "-" << mt << "-" << tg << "-" << se << endl;
//	for (int i = 1; i <= zonmax; i++)
//	{
//		rivflow << "zone:\t" << i;
//		rivflow << "\tflowm3s:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].flowm3s << "\tinfilt_amount:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].infilt_amount << endl;
//	}
//	rivflow.close();
//#endif
// #pragma omp parallel for schedule(runtime) private(ij)
	for (i=1;i<=anzrout;i++) 
	{
		routing[i]->Routing(jahr,monat,tag,stunde,minute,actsubtimestep,numsubtimesteps);
		if (riverflow_table && (actsubtimestep == numsubtimesteps-1))
		{
			routing[i]->riverflow_table_update();
		}
	}

	/* updating maxirrig: in a first step, the maximum allowed irrigation abstractions should be added to the downstream subbasin in order to allow there a 
	additionally abstraction if the actual subbasin does not need all possible irrigation water 
rtarget[0][1..n] order of routings 1..n 
rtarget[1][1..zonmax] targets of subbasins 
rtarget[2][1..zonmax] index of the routing object instance for the actual subbasin 
in Qges[0][1..zonmax] are the routed specific fluxes (with respect to their real catchment area) in mm
in maxirrig[0][1..zonmax] are the potential values for irrigation as specific fluxes with respect to the subbasin area only, but
there are no considerations of exceeding amounts from upstream areas that could be used in downstream areas This is done now! */

	//  for (i=1;i<=zonmax;i++)
	//    maxirrig_temp[i]=maxirrig[0][i]; /* maxirirg will be changed recursively, but in each recursion, the original value is also required. Thats why it is stored in a temporal array */

	if (zoneSW)
	{
		for (i=1;i<=zonmax;i++)
		{
			routindex=rtarget[0][i];
			actindex=routindex;
			if (routindex>0)
			{
				/* irr_diff is the difference between irrmax for the actual subbasin and the actual runoff in Qges, 
		both as specific discharges with respect to the area of the actual subbasin */
				if (rtarget[2][actindex]>0)
				{
					/* if this subbasin is the target of other subbasins (i.e., this basin has a relevant relative subbasin area) */
					//        irr_diff=maxirrig_temp[actindex]-(routing[rtarget[2][actindex]]->GetAbstractionSum()+Qges[0][actindex])*
					//                                        abfrel[0][actindex]/zoneSW[actindex].rel/routing[rtarget[2][actindex]]->TgAe_rel;
					irr_diff=zoneSW[actindex].sx; /* the irrigation-amount was subtracted from the runoff. Since it is assumed in 
										a first estimate that the next interval will bring the same runoff, the maxirrig-array can be 
										corrected by adding the irrigation-amount from the last timestep to the expected runoff */
				}
				else
				{
					irr_diff=maxirrig_temp[actindex]-Qges[0][actindex]*abfrel[0][actindex]/zoneSW[actindex].rel; /* a 1. order subbasin without inflows and tributaries */
				}
				if (irr_diff<0) 
				{
					irr_diff=0; 
				}
				while (rtarget[1][routindex]>0)
				{
					routindex=rtarget[1][routindex];
					if (routindex>0)
					{
						maxirrig[0][routindex]+=irr_diff*zoneSW[actindex].rel/zoneSW[routindex].rel;
					}
				}
			}
		}
	}

	/* Effizienzkriterium updaten */
	if ((spendenende<10)&&(SpendenBuffer->GetStatus()!=FB_STATE_FILENOTFOUND))
	{
		time++;
		for (i=1;i<=zonmax;i++)
		{

			int wertcnt=0;
			if ((Qstatcode/1000==2)||(Qstatcode/1000==4))
			{
				LastResults[lastresultsindex].wert[i]=Spende[11][Spalte[i]];
			}
			else
			{
				if ((Qstatcode/1000==5)||(Qstatcode/1000==6))
				{
					if (wertcnt<numsubb)
					{
						if (i==listsubb[wertcnt+1])
						{
							LastResults[lastresultsindex].wert[i]=Spende[11][Spalte[i]];
							wertcnt++;
						}
					}
				}
			}
			

			if ((Spalte[i]>0)&&(time>timeoffset))
			for (j=0;j<=20;j++)
			{
				if (Spende[j][Spalte[i]]>0.0)
				{
					double ModVal = Qges[0][i];
					double ObsVal = Spende[j][Spalte[i]];
					if ((do_it == 1) || (do_it == 2)) // do_it = 1 is routing with modelled runoff in mm/time step, do_it = 2 is routing with observed runoff in mm/time step
					{ // so we need to convert both, mpodelled and observed runoff into m3/s, if runoff_measure==1
						if (rtarget[2][i] > 0)
						{
							if (runoff_measure == 1)
							{
								ModVal = (Qges[0][i]		   * (abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV)) / routing[rtarget[2][i]]->TgAe_rel;
								ObsVal = (Spende[j][Spalte[i]] * (abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV)) / routing[rtarget[2][i]]->TgAe_rel;
							}
						}
						else
						{
							if (runoff_measure == 1)
							{
								ModVal = (Qges[0][i]           * abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV);
								ObsVal = (Spende[j][Spalte[i]] * abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV);
							}
						}
					}
					else // do_it = 3 means routing with observed runoff in m3/s, so if runoff_measure is 1, only ModVal has to be converted into m3/s
					{    // nut if runoff_measure is not 1, ObsVal has to be coinverted into mm/time step
						if (rtarget[2][i] > 0)
						{
							if (runoff_measure == 1)
							{
								ModVal = (Qges[0][i] * (abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV)) / routing[rtarget[2][i]]->TgAe_rel;
							}
							else
							{
								ObsVal = (Spende[j][Spalte[i]] * InV * routing[rtarget[2][i]]->TgAe_rel) / (abfrel[0][0] * abfrel[0][i] * (double)16.6666667);
							}
						}
						else
						{
							if (runoff_measure == 1)
							{
								ModVal = (Qges[0][i] * abfrel[0][0] * abfrel[0][i] * (double)16.6666667 / InV);
							}
							else
							{
								ObsVal = (Spende[j][Spalte[i]] * InV) / (abfrel[0][0] * abfrel[0][i] * (double)16.6666667);
							}
						}
					}
					/* fuer lineares rsquare */
					effizienz[j][i].ngem    ++;
					effizienz[j][i].sxgem   += ObsVal;
					effizienz[j][i].sx2gem  += ObsVal * ObsVal;

					effizienz[j][i].ndiff   ++;
//					temp                     = Qges[0][i] - ObsVal;
					temp                     = ModVal - ObsVal;
					effizienz[j][i].sxdiff  += temp;
					effizienz[j][i].sx2diff += temp*temp;
					/* fuer logarithmisches rsquare */
//					if (Qges[0][i]>0.0)
					if ((ModVal>0.0)&&(ObsVal>0))
					{
						logwert = log(ObsVal);
						log_effi[j][i].ngem    ++;
						log_effi[j][i].sxgem   += logwert;
						log_effi[j][i].sx2gem  += logwert*logwert;

						log_effi[j][i].ndiff   ++;
//						temp                    = log(Qges[0][i]) - logwert;
						temp                    = log(ModVal) - logwert;
						log_effi[j][i].sxdiff  += temp;
						log_effi[j][i].sx2diff += temp*temp;
					}
				}
			}
		}
	}

	/* Drehen der Zeiger */
	if (SpendenBuffer->GetStatus()!=FB_STATE_FILENOTFOUND)
	{
		if (spendenende>=10) 
		(*WSstream) << "INFO: observed runoff: actual date "<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde
		<<"h not found (end of discharge file at "<<Ds<<"."<<Ms<<"."<<Ys<<", "<<Hs<<"h)"<<"[end-entry]";
		zeiger=Spende[0];
		for (i=0;i<20;i++) Spende[i]=Spende[i+1];
		Spende[20]=zeiger;

		/* Einlesen der nächsten Spenden */
		if (spendenende<10)
		{
			for (i=1;i<=spendenanzahl;i++)
			{
				Spende[20][i]=SpendenBuffer->GetActValue(i-1);
			}
			SpendenBuffer->NextRecord();
			if (SpendenBuffer->GetStatus()==FB_STATE_LASTRECORDREAD)
			{
				(*WSstream) << "INFO: end of discharge file reached!"<<"[end-entry]";
				spendenende++;
			}
			else
			{
				SpendenBuffer->GetActdate(Ys,Ms,Ds,Hs,ms);
				if (Ys<1000) Ys+=1900;
			}
		}
		else
		{
			for (i=1;i<=spendenanzahl;i++) Spende[20][i]=0.0;
		}
	}

	/* Schreiben der Ausgabewerte in die Statistikdatei */
	Qstat->init();
	for (j=1;j<=zonmax;j++) 
	{ 
		zone[j].n++;
		zone[j].RelCellSum+=1.0;
		if (rtarget[2][j]>0)
		{
			zone[j].sx += runoff_measure==1 ?
			(Qges[0][j]*(abfrel[0][0]*abfrel[0][j]*(double)16.6666667/InV))/routing[rtarget[2][j]]->TgAe_rel
			:  Qges[0][j];
		}
		else
		{
			zone[j].sx += runoff_measure==1 ? (Qges[0][j]*abfrel[0][0]*abfrel[0][j]*(double)16.6666667/InV):Qges[0][j];
		}
	}
	Qstat->calc_and_done(jahr,monat,tag,stunde,minute);
	/* Schreiben der Konzentrationen in die Statistikdatei(en) */
	for (z=1;z<=tracer_count;z++)
	{
		tracer_module[z].stat->init();
		for (j=1;j<=zonmax;j++)
		{
			tracer_module[z].zone[j].n++;
			tracer_module[z].zone[j].RelCellSum+=1.0;
			tracer_module[z].zone[j].sx += Qges[z][j];
		}
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde,minute);
	}
	if (mpi_rank==0) 
	{
		MPIScatterGridWithMask(pondindex,lakeindex,lakeCellsPerMPIRank); // will be called for rank >0 at the beginning of the ::run method
	}
}

/* neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void abflussmodellklasse::ClearMemory()
{
	int i,z,k;
	if (do_it > 0)
	{

		if (initialized) WriteStorageFile();
		if (storagefilename) { free(storagefilename); storagefilename = NULL; }
		substance_transport_delete();
		/* Schreiben der R-square in Ergebnisdatei an die letzte Stelle */
		if (initialized)
		{
			if (SpendenBuffer->GetStatus()!=FB_STATE_FILENOTFOUND)
			if (mpi_rank==0) 
			{
				Qstat->Rsquare(effizienz,log_effi);
			}
			SpendenBuffer->ClearMemory();
			delete SpendenBuffer;
		}
		Qstat->ClearMemory();
		if (Qstat)
		{
			delete(Qstat);
			Qstat=NULL;
		}

		if (routing) 
		{
			for (i=1;i<=anzrout;i++)
			routing[i]->ClearMemory();
			for (i=1;i<=anzrout;i++)
			{
				if (routing[i])
				{
					delete (routing[i]);
					routing[i] = NULL;
				}
			}
			free(routing); 
			routing = NULL;
		}

		if (Q12ges) 
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (Q12ges[z])
				{
					for (k=0;k<=p_subint;k++)
					{
						if (Q12ges[z][k])
						{
							free(Q12ges[z][k]);
							Q12ges[z][k]=NULL;
						}
					}
					free(Q12ges[z]);
					Q12ges[z]=NULL;
				}
			}
			free(Q12ges);
			Q12ges=NULL;
		}

// comes from unsatzonmodel?
		if (Qges)
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (Qges[z])
				{
					free(Qges[z]);
					Qges[z]=NULL;
				}
			}
			free(Qges);
			Qges=NULL;
		}

		if (fluss) 
		{
			for (z=0;z<=tracer_count;z++) 
			{ 
				if (fluss[z]) 
				{
					free(fluss[z]); 
					fluss[z]=NULL; 
				}
			}
			free(fluss); 
			fluss=NULL;
		}
		
		if (effizienz)
		{
			for (i=0;i<=20;i++)
			{
				if (effizienz[i]) 
				{
					free(effizienz[i]); 
					effizienz[i] = NULL;
				}
			}
			free(effizienz); 
			effizienz = NULL;
		}

		if (log_effi)
		{
			for (i=0;i<=20;i++)
			{
				if (log_effi[i]) 
				{
					free(log_effi[i]); 
					log_effi[i] = NULL;
				}
			}
			free(log_effi); 
			log_effi = NULL;
		}

		if (spendenfilename) 
		{
			free(spendenfilename); 
			spendenfilename = NULL;
		}
		
		if (Qstatfilename) 
		{
			free(Qstatfilename);   
			Qstatfilename   = NULL;
		}
		
		if (Spende) 
		{
			for (i=0;i<=20;i++) 
			{
				if (Spende[i]) 
				{
					free(Spende[i]); 
					Spende[i] = NULL;
				}
			}
			free(Spende); 
			Spende = NULL;
		}
		if (Spalte) { free(Spalte); Spalte = NULL; }
		
		if (rtarget)
		{
			for (i=0;i<3;i++)
			{
				if (rtarget[i]) 
				{
					free(rtarget[i]); 
					rtarget[i]=NULL;
				}
			}
			free(rtarget); 
			rtarget=NULL;
		}

		if (abfrel) 
		{
			if (abfrel[0]) {			free(abfrel[0]); 			abfrel[0]=NULL;}
			if (abfrel[1]) {			free(abfrel[1]); 			abfrel[1]=NULL;}
			free(abfrel); 
			abfrel = NULL;
		}

		if (num_ST) { free(num_ST); num_ST = NULL; }
		
		if (maxirrig) 
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (maxirrig[z])
				{
					free(maxirrig[z]);
					maxirrig[z] = NULL;
				}
			}
			free(maxirrig); 
			maxirrig=NULL;
		}
		if (maxirrigFromReservoirs) 
		{
			for (z=0;z<=tracer_count;z++)
			{
				if (maxirrigFromReservoirs[z])
				{
					free(maxirrigFromReservoirs[z]);
					maxirrigFromReservoirs[z] = NULL;
				}
			}
			free(maxirrigFromReservoirs); 
			maxirrigFromReservoirs=NULL;
		}
		if (maxirrig_temp) { free(maxirrig_temp); maxirrig_temp = NULL; }
		if (lakeCellsPerMPIRank) delete[] lakeCellsPerMPIRank;

		mapReservoirLevels->clear();
		mapTributariesWithHydCon->clear();
		delete mapReservoirLevels;
		delete mapTributariesWithHydCon;
		mapReservoirLevels = NULL;
		mapTributariesWithHydCon = NULL;

		if (additionalinfo)
		{
			delete[] additionalinfo;
			additionalinfo = NULL;
		}

	}
	if ((PRECONLAKES)&&(mpi_rank==0))
	{
		delete_grid(preconlakesindex);
		free_matrix(PRECONLAKES,gridinfo[preconlakesindex].nrl,gridinfo[preconlakesindex].nrh,1,ncols,nrows);
	}
	modelobject::ClearMemory();
}

void abflussmodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	try
	{
		int i,j;
		char *identify;
		if (do_it==0) return;
		identify=newstring(300,1);
		tracer_count = p_tracer_count;
		if (tracer_count>0)
		{
			tracer_info = p_tracer_info_table;
			tracer_module = (tracer_module_type*)calloc(tracer_count+1,sizeof(tracer_module_type));
			for (i=1;i<=tracer_count;i++)
			{
				if (tracer_info[i].radioactive==1)
				{
					/* this is the factor, with which the actual concentration of a radioactive
			tracer has to be multiplied to get the radioactive decay for one time step */
					if (tracer_info[i].half_time<=0) { tracer_info[i].decfact = 0.0; tracer_info[i].routdecfact = 0.0; }
					else
					{ tracer_info[i].decfact = exp(-(InV/60/24)/(tracer_info[i].half_time*_1_ln2_));
						tracer_info[i].routdecfact = exp(-(InV/p_subint/60/24)/(tracer_info[i].half_time*_1_ln2_));
					}
				}
				/* initializing a statistc register set for this tracer */
				tracer_module[i].zone = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
				if (!tracer_module[i].zone) speichererror(1);
				for (j=0;j<=zonmax;j++)
				{ tracer_module[i].zone[j].n          = tracer_module[i].zone[j].sx        = 0;
					tracer_module[i].zone[j].zonenwert  = zone[j].zonenwert;
					tracer_module[i].zone[j].rel        = zone[j].rel;
					tracer_module[i].zone[j].RelCellSum = zone[j].RelCellSum;
				}

				/* initializing the output stream for the tracer statistics */
				tracer_module[i].statfilename=newstring(_stringsize_+1,1);
				sprintf(tracer_module[i].statfilename,"%s%.3s%.5s.%.3s",tracer_info[i].path,tracer_info[i].name,instanzname,tracer_info[i].statext);
				tracer_module[i].stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,NULL,1,1,&tracer_module[i].zone,(int)zonmax,(int)tracer_info[i].statcode,(int)99,(char*)tracer_module[i].statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
				if (retvalue!=0) throw retvalue;
				sprintf(identify,"tracer concentration for tracer %d (%s) ",i,tracer_info[i].name);
				tracer_module[i].stat->schreibzonenwerte(identify,instanzname,LastResults,maxLastResults);

			}
		}
		if (identify) { free(identify); identify = NULL; }
	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	// catching all exceptions in order to allow the object to be constructed correct. This is required 
	// in order to allow a dextructor call in case of an error during construction
	{
		retvalue=103;
	}
}

void abflussmodellklasse::substance_transport_delete()
{
	int i;
	if (do_it==0) return;
	if (tracer_count>0)
	{
		for (i=1;i<=tracer_count;i++)
		{
			tracer_module[i].stat->ClearMemory();
			if (tracer_module[i].stat) delete(tracer_module[i].stat);
			if (tracer_module[i].statfilename) 
			{
				free(tracer_module[i].statfilename);
				tracer_module[i].statfilename=NULL;
			}
			if (tracer_module[i].zone) 
			{
				free(tracer_module[i].zone);
				tracer_module[i].zone=NULL;
			}
		}
	}
}

void abflussmodellklasse::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
		if (Qstat) Qstat->flushbuffer();
		for (z=1;z<=tracer_count;z++)
		{
			tracer_module[z].stat->flushbuffer();
		}
		for (int i=1;i<=anzrout;i++)
		{
			routing[i]->flushbuffers();
		}
	}
}

void abflussmodellklasse::initialize_lakes(std::string instanzname)
{
	char *tstr;

	if (!ist_vorhanden(_LAKE_CODES_))
	{
		(*WSstream) << "INFO: no LAKE-grid was found in "<<this->instanzname<<"\n"<<
		"Lake modelling is switched off for "<<this->instanzname<<".\n"<<"[end-entry]";
		//	throw(223);
	}
	else
	{
		tstr = new char[_stringsize_];
		strcpy(tstr,instanzname.c_str());
		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		findeabschnitt(tstr);
		fscanfln(steuerdatei,"%d",&lake_active,"flow control parameter for ",tstr);
		// switch off lake model when in Spinup-mode
		if (doSpinUp || routingmode) {
			lake_active=0;
		}
		if (lake_active)
		{
			for (int i=1;i<=3;i++) 
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
			fscanfln(steuerdatei,"%d",&readgrids,"readgrids for ",tstr);
		}
		else
		{
			readgrids=0;
		}
		delete[] tstr;
	}
}

void abflussmodellklasse::get_reservoir_vector(std::vector<SP_class*> &reservoirs)
{
	for (int i = 1; i <= anzrout; i++)
	{
		routing[i]->add_reservoir(reservoirs);
	}
}

void abflussmodellklasse::routingmode_initialize()
{
	Qges_uebernommen = new double*[tracer_count + 1];
	for (int z = 0; z <= tracer_count; z++)
	{
		Qges_uebernommen[z] = new double[zonmax + 1];
	}
	// read filename for qges from control file into qgesfileRMname
	if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
	{
		(*WSstream) << "FATAL: error opening " << steuerdateiname << " (in initialize routingmode)"<<"[end-entry]";
		throw(28);
	}
	strcpy(tstr, "[");
	strcat(tstr, _Unsatzonmodell_);
	strcat(tstr, "]");
	findeabschnitt(tstr);

	// over-reading lines until qges filename must be read in
	qgesfileRMname = newstring(_stringsize_ + 1, 1);

	for (int i = 1; i <= 20; i++)
	{
		tstr = fgets(tstr, _stringsize_, steuerdatei);
	}
	int tcode;
	fscanfln(steuerdatei, "%s%d", qgesfileRMname, &tcode, "name of the Qgstatfiles in ", tstr); /* total flow mm/interval */
	ExtendFileName(DefaultOutputDirectory, qgesfileRMname);

	if (fclose(steuerdatei) == EOF)
	{
		(*WSstream) << "FATAL: control file \"" << steuerdateiname << "\" cannot be closed (in initialize routingmode)" << "[end-entry]";
		throw(11);
	}

	// create new filebuffer instance
	QgesBufferRM = new FileBuffer(qgesfileRMname, 0, 1000000);
	if (QgesBufferRM->GetStatus() == FB_STATE_FILENOTFOUND || QgesBufferRM->GetStatus() == FB_STATE_FILENOTGOOD)
	{
		(*WSstream) << "ERROR: error opening " << qgesfileRMname
			<< "\nseparate routing with -mode:routing can only be done with valid input data as result of a former model run (total runoff for each subbasin)!" << "[end-entry]";
		throw 1400;
	}
	else /* Nur dann, wenn auch eine Spendendatei gefunden wurde, werden die folgenden Sequenzen noetig sein */
	{

		(*WSstream) << "INFO: looking for starting date in " << qgesfileRMname << "[end-entry]";
		int QgesCount = QgesBufferRM->GetValcnt();
		// if there are less columns in the external file than there are zones, or subcatchments to rout, that is, 
		// an erorr must be thrown (note that the total value for all subbasins is expected here, too, so there must be one column more than sub basins)
		if (QgesCount != zonmax+1)
		{
			(*WSstream) << "ERROR: error: number of expected columns in " << qgesfileRMname << " doesn't match the number of subbasins defined in the control file!" << "[end-entry]";
			throw 1410;
		}

		qgesfileend = 0;
		// wind fast forward to the starting date and read the first values
		Ys = 0;
		int oldDs = 0;
		QgesBufferRM->SetBufPos(-1); // auf eine Position vor Beginn setzen
		while ((Ys != jahra) || (Ms != monata) || (Ds != taga) || (Hs != stundea) || (ms != minutea))
		{
			QgesBufferRM->NextRecord();
			QgesBufferRM->GetActdate(Ys, Ms, Ds, Hs, ms);
			if (Ds < oldDs)
			{
				(*WSstream) << ".";
				cout.flush();
			}
			oldDs = Ds;
			if (Ys < 1000) Ys += 1900;
			if (QgesBufferRM->GetStatus() == FB_STATE_LASTRECORDREAD)
			{
				qgesfileend++;
				break;
			}

		}
		(*WSstream) << "[end-entry]";
		// if starting date is not found, throw an error
		if ((Ys != jahra) || (Ms != monata) || (Ds != taga) || (Hs != stundea) || (ms != minutea))
		{
			(*WSstream) << "ERROR: error looking for start date " 
				<< wasim::setw(4) << wasim::setfill('0') << jahra   <<"-"
				<< wasim::setw(2) << wasim::setfill('0') << monata  <<"-" 
				<< wasim::setw(2) << wasim::setfill('0') << taga    <<" "
				<< wasim::setw(2) << wasim::setfill('0') << stundea <<":"
				<< wasim::setw(2) << wasim::setfill('0') << minutea <<" (run into error at actual position: "
				<< wasim::setw(2) << wasim::setfill('0') << Ys <<"-"
				<< wasim::setw(2) << wasim::setfill('0') << Ms <<"-"
				<< wasim::setw(2) << wasim::setfill('0') << Ds <<" "
				<< wasim::setw(2) << wasim::setfill('0') << Hs <<":"
				<< wasim::setw(2) << wasim::setfill('0') << ms <<") in " << qgesfileRMname << "[end-entry]";
			throw 1420;
		}

	}
	for (int i = 1; i <= zonmax; i++)
	{
		Qges_uebernommen[0][i] = QgesBufferRM->GetActValue(i - 1);
	}
}

void abflussmodellklasse::routingmode_nextrecord()
{
	// read next record from filebuffer and, if required, refill the buffer from the input file
	QgesBufferRM->NextRecord();
	for (int i = 1; i <= zonmax; i++)
	{
		Qges_uebernommen[0][i] = QgesBufferRM->GetActValue(i-1);
	}
}

void abflussmodellklasse::routingmode_cleanup()
{
	// free all memory in the filebuffer
	if (qgesfileRMname) 
	{
		free(qgesfileRMname);
		qgesfileRMname = NULL;
	}

}

bool abflussmodellklasse::routingmode_isend()
{
	int min, h, d, m, y;
	QgesBufferRM->GetActdate(y, m, d, h, min);
	if ((QgesBufferRM->GetStatus() >= 2 && QgesBufferRM->GetBufpos()>= QgesBufferRM->GetBufmax()) || 
		((y == jahre) && (m == monate) && (d == tage) && (h == stundee) && (min == minutee)))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void abflussmodellklasse::getdatetime(long *me, long *se, long *tg, long *mt, long *jr)
{
	/* Die uebergebenen Parameter werden mit dem aktuell in der
	Instanz gueltigen Zeit (jahr,Monat,tag und stunde) geladen.
	So koennen andere Teilmodelle mit den fuer die Interpolation gueltigen Zeiten
	versorgt werden */
	int ime, ise, itg, imt, ijr;
	QgesBufferRM->GetActdate(ijr, imt, itg, ise, ime);
	*me = ime;
	*se = ise;
	*tg = itg;
	*mt = imt;
	*jr = ijr;
}



/* ~~~~~~~~~~~~~~~~~~~~ Ende der Abflussmodellroutinen ~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~ Speicher/Reservoirs/Lakes ~~~~~~~~~~~~~~~~~ */

SP_class::SP_class() :
WSstream(), tabmax(0), Vmax(0), Vakt(0), Vold(0), Vtemp(0), conc(0), qu(0), rulecnt(0), JD(0),
Abakt(0), ARTabIndex(0), ARTabValue(0), ARChangeRate(0),
SP_buffer(0), SP_bufmax(0), SP_bufcnt(0), SP_valcnt(0),
SP_filename(0), SP_file(0),
LAKE(0), MAXPOND(0), POND(0),
WaterLevel(0), Volume(0),
startcol(0), startrow(0), endcol(0), endrow(0),
SP_code(0), TG(0), tracer_count(0), readgrids(0),
year(0), month(0), day(0), hour(0), dTmin(60), UseRunningAvg(0), targetCapacity(0),
WeekDays(0), start_hour(0), stop_hour(0), start_hour2(0), stop_hour2(0), bReverseLookUp(false),
useHydraulicConnection(false), TargetReservoir(0), ChannelWidth(0),
ChannelBottom(0), ChannelLength(0), RoughnessCoeff(0), 
TargetWaterLevel(0), TargetWaterLevelChangeRate(0), actWaterLevel(0),
ABSTRmm(0), ABSTRm3s(0), ABSTRm3(0),
IRRIGmm(0), IRRIGm3s(0), IRRIGm3(0),
PRECmm(0), PRECm3s(0), PRECm3(0),
QINmm(0), QINm3s(0), QINm3(0),
QOUTmm(0), QOUTm3s(0), QOUTm3(0),
QGWmm(0), QGWm3s(0), QGWm3(0),
ETRmm(0), ETRm3s(0), ETRm3(0),
DVmm(0), DVm3s(0), DVm3(0),
ERRmm(0), ERRm3s(0), ERRm3(0),
useMinutes(false),hour_offset(0),timestep(0), deltaPond(0)
{
	mon[1]=31;
	mon[2]=28;
	mon[3]=31;
	mon[4]=30;
	mon[5]=31;
	mon[6]=30;
	mon[7]=31;
	mon[8]=31;
	mon[9]=30;
	mon[10]=31;
	mon[11]=30;
	mon[12]=31;

}

SP_class::~SP_class()
{
	if (SP_code) // delete initialized objects only
	{
	// allocate memory for output buffers
		for (int k=0;k<=SP_bufmax;k++)
		{
			if (SP_buffer[k])
			delete[] SP_buffer[k];
			SP_buffer[k]=NULL;
		}
		if (SP_buffer) delete[] SP_buffer;
		if (Vakt) free(Vakt);
		if (Vold) free(Vold);
		if (Vtemp) free(Vtemp);
		if (conc) free(conc);
		if (ARTabIndex) delete[] ARTabIndex;
		if (ARTabValue)
		{
			for (int i = 0; i <= rulecnt; i++)
			{
				delete[] ARTabValue[i];
			}
			delete[] ARTabValue;
		}
		if (ARChangeRate)
		{
			for (int i = 0; i <= rulecnt; i++)
			{
				delete[] ARChangeRate[i];
			}
			delete[] ARChangeRate;
		}
		if (WaterLevel) delete[] WaterLevel;
		if (Volume) delete[] Volume;
		if (SP_filename) delete[] SP_filename;
		if (JD) delete[] JD;
		if (WeekDays) delete[] WeekDays;
		if (start_hour) delete[] start_hour;
		if (stop_hour) delete[] stop_hour;
		if (start_hour2) delete[] start_hour2;
		if (stop_hour2) delete[] stop_hour2;
		if (targetCapacity) delete[] targetCapacity;
	}
}

void SP_class::SP_flushbuffer(int code)
{
	unsigned long i,j;
	if ((code>0)&&(SP_bufcnt>0))
	{

		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		(*WSstream) << "flushing to disk: "<<SP_filename<<" ("<<buffersize<<" buffersize)"<<"[end-entry]";
		if ((SP_file=fopen(SP_filename,"a+"))==NULL)
		{
			(*WSstream) << "FATAL: error opening reservoir file "<<SP_filename<<"[end-entry]";
			throw(65);
		}
		for (i=1;i<=SP_bufcnt;i++)
		{
			if (useMinutes)
			{
				fprintf(SP_file, "%2.0lf\t%2.0lf\t%2.0lf\t%2.0lf\t%2.0lf", SP_buffer[i][0], SP_buffer[i][1], SP_buffer[i][2], SP_buffer[i][3], SP_buffer[i][4]);
			}
			else
			{
				fprintf(SP_file, "%2.0lf\t%2.0lf\t%2.0lf\t%2.0lf", SP_buffer[i][0], SP_buffer[i][1], SP_buffer[i][2], SP_buffer[i][3]);
			}
			for (j = 5; j <= 4 + SP_valcnt; j++) // also when useMinutes is false, the values are stored in element 5 and higher
			{
				fprintf(SP_file, "\t%lg", SP_buffer[i][j]);
			}
			fprintf(SP_file, "\n");
		}
		fclose(SP_file);
		SP_bufcnt=SP_valcnt=0;
		SetSignalAbfangfunktion();
	}
}

void SP_class::FillMapReservoirLevels(std::map<int, double**> *p_mapReservoirLevels)
{
	mapReservoirLevels = p_mapReservoirLevels;
	double** temppointers = new double*[3];
	temppointers[0] = &(this->actWaterLevel);
	temppointers[1] = &(this->DVm3s);
	temppointers[2] = &(this->AREAm2);

	mapReservoirLevels->insert(pair<int, double**>(this->SP_code, temppointers));
}

void SP_class::FillMapTributariesWithHydCon(std::map<int, bool> *p_mapTributariesWithHydCon)
{
	mapTributariesWithHydCon = p_mapTributariesWithHydCon;
	mapTributariesWithHydCon->insert(pair<int, bool>(this->TG, this->useHydraulicConnection==true));
}



void SP_class::ReadReservoirRule(FILE *steuerdatei)
{
	char *tstr;
	tstr = (char*)calloc(_stringsize_, sizeof(char));
	sprintf(tstr, "[abstraction_rule_reservoir_%d]", SP_code);
	findeabschnitt(steuerdatei, tstr);
	if (!feof(steuerdatei)) tstr = fgets(tstr, _stringsize_, steuerdatei);
	// Anzahl der Eintraege der Abgaberegel
	sscanf(tstr, "%d", &(tabmax));
	if (tabmax <= 0)
	{
		// looking for String "HydraulicConnection"
		if (!feof(steuerdatei)) tstr = fgets(tstr, _stringsize_, steuerdatei);
		std::string s(tstr);
		if (s.find("HydraulicConnection") != std::string::npos)
		{
			this->useHydraulicConnection = true;
			// read in the entire definition with all parameters and analyse the concatanated string like in a landuse table entry
			while (s.find("}") == std::string::npos)
			{
				tstr = fgets(tstr, _stringsize_, steuerdatei);
				s.append(std::string(tstr));
			}
			// now extracting parameters
			std::string buffer;
			int index1;
			int index2;
			try
			{
				buffer = s.substr(s.find("TargetReservoir"));
				index1 = buffer.find("=") + 1;
				index2 = buffer.find_first_of(";\n");
				buffer = buffer.substr(index1, index2 - index1);
				TargetReservoir = atof(buffer.c_str());

				buffer = s.substr(s.find("ChannelWidth"));
				index1 = buffer.find("=") + 1;
				index2 = buffer.find_first_of(";\n");
				buffer = buffer.substr(index1, index2 - index1);
				ChannelWidth = atof(buffer.c_str());

				buffer = s.substr(s.find("ChannelBottom"));
				index1 = buffer.find("=") + 1;
				index2 = buffer.find_first_of(";\n");
				buffer = buffer.substr(index1, index2 - index1);
				ChannelBottom = atof(buffer.c_str());

				buffer = s.substr(s.find("ChannelLength"));
				index1 = buffer.find("=") + 1;
				index2 = buffer.find_first_of(";\n");
				buffer = buffer.substr(index1, index2 - index1);
				ChannelLength = atof(buffer.c_str());

				buffer = s.substr(s.find("RoughnessCoeff"));
				index1 = buffer.find("=") + 1;
				index2 = buffer.find_first_of(";\n");
				buffer = buffer.substr(index1, index2 - index1);
				RoughnessCoeff = atof(buffer.c_str());
			}
			catch (...)
			{
				(*WSstream) << "FATAL: reading HydraulicConnection parameters for SP " << SP_code << ": some values could not be read\n" << s.c_str() << "[end-entry]";
				throw(29);
			}

		}
	}
	else
	{
		// analysing the number of julian days
		std::string s(tstr);
		rulecnt = 0;
		int pos = s.find_first_not_of(" \t\n\r\0");
		while ((s.length() > 0) && (pos != std::string::npos))
		{
			// remove leading whitespaces
			if (s.find_first_not_of(" \t") > 0)
			{
				s = s.substr(s.find_first_not_of(" \t"));
			}
			s = s.substr(s.find_first_not_of("0123456789.eE+-"));
			rulecnt++;
			pos = s.find_first_not_of(" \t\n\r\0");
		}
		// since the first entry was not a Julian Day, the rulecnt is 1 to high (except for old style rules with no JD at all)
		if (rulecnt > 1) rulecnt--;
		// now the arrays can be allocated
		ARTabIndex = new double[tabmax + 1];
		ARTabValue = new double*[rulecnt + 1];
		ARChangeRate = new double*[rulecnt + 1];
		for (int i = 0; i <= rulecnt; i++)
		{
			ARTabValue[i] = new double[tabmax + 1];
			ARChangeRate[i] = new double[tabmax + 1];
		}
		targetCapacity = new double[rulecnt + 1];
		JD = new int[rulecnt + 1];
		WeekDays = new long[rulecnt + 1];
		start_hour = new int[rulecnt + 1];
		stop_hour = new int[rulecnt + 1];
		start_hour2 = new int[rulecnt + 1];
		stop_hour2 = new int[rulecnt + 1];
		// read Julian Days from the first row (from std::string s)
		s = std::string(tstr);
		// overide first value (is counter of rows)
		if (rulecnt == 1)
		{
			JD[1] = 1;
		}
		else
		{
			if (s.find_first_not_of(" \t") > 0)
			{
				s = s.substr(s.find_first_not_of(" \t"));
			}
			s = s.substr(s.find_first_not_of("0123456789.eE+-"));
			for (int i = 1; i <= rulecnt; i++)
			{
				std::stringstream ss(s.c_str());
				ss >> JD[i];
				if (s.find_first_not_of(" \t") > 0)
				{
					s = s.substr(s.find_first_not_of(" \t"));
				}
				s = s.substr(s.find_first_not_of("0123456789.eE+-"));
			}
		}
		// now reading the table
		for (int t = 1; t <= tabmax; t++)
		{
			if (!feof(steuerdatei)) tstr = fgets(tstr, _stringsize_, steuerdatei);
			s = std::string(tstr);
			for (int i = 0; i <= rulecnt; i++)
			{
				std::stringstream ss(s.c_str());
				if (i == 0)
				{
					ss >> ARTabIndex[t];
				}
				else
				{
					ss >> ARTabValue[i][t];
				}
				if (s.find_first_not_of(" \t") > 0)
				{
					s = s.substr(s.find_first_not_of(" \t"));
				}
				s = s.substr(s.find_first_not_of("0123456789.eE+-"));
			}
		}
		Vmax = ARTabIndex[tabmax];

		// initializing additional parameters in case they are not present in the control file
		for (int i = 1; i <= rulecnt; i++)
		{
			targetCapacity[i] = -1.0;
			WeekDays[i] = -1;
			start_hour[i] = -1;
			stop_hour[i] = -1;
			start_hour2[i] = -1;
			stop_hour2[i] = -1;
		}

		// reading target capacities
		if (!feof(steuerdatei))
		{
			tstr = fgets(tstr, _stringsize_, steuerdatei);
			s = std::string(tstr);
		}
		while ((s.find("[") == std::string::npos) && (s != "\n") && (s.length() != 0) && !feof(steuerdatei))
		{
			if (s.find("TargetCap") != std::string::npos)
			{
				// reading target capacities
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> targetCapacity[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("0123456789.eE+-"));
				}
			}
			else if (s.find("WeekDays") != std::string::npos) {
				// Reading WeekDays
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> WeekDays[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567"));
				}
			}
			else if (s.find("start_hour") != std::string::npos) {
				// Reading start_hour
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> start_hour[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("stop_hour") != std::string::npos) {
				// Reading stop_hour
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> stop_hour[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("start2_hour") != std::string::npos) {
				// Reading start_hour2
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> start_hour2[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("stop2_hour") != std::string::npos) {
				// Reading stop_hour2
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> stop_hour2[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("TableIndex") != std::string::npos) {
				s = s.substr(s.find("=") + 1);
				if (s.find("Flow") != std::string::npos)
				{
					bReverseLookUp = true;
				}
				else
				{
					bReverseLookUp = false;
				}
			}
			if (!feof(steuerdatei))
			{
				tstr = fgets(tstr, _stringsize_, steuerdatei);
				s = std::string(tstr);
			}
		}

		// interpolating between the table entries or table indices, depending on the state of breverseLookUp
		for (int i = 1; i <= rulecnt; i++)
		{
			if (bReverseLookUp)
			{
				for (int t = 1; t < tabmax; t++)
				{

					ARChangeRate[i][t] = (ARTabValue[i][t + 1] - ARTabValue[i][t]) != 0 ?
						(ARTabIndex[t + 1] - ARTabIndex[t]) /
						(ARTabValue[i][t + 1] - ARTabValue[i][t]) : 0.0;
				}
				ARChangeRate[i][tabmax] = 0.0;
			}
			else // default
			{
				for (int t = 1; t < tabmax; t++)
				{

					ARChangeRate[i][t] = (ARTabIndex[t + 1] - ARTabIndex[t]) != 0 ?
						(ARTabValue[i][t + 1] - ARTabValue[i][t]) /
						(ARTabIndex[t + 1] - ARTabIndex[t]) : 0.0;
				}
				ARChangeRate[i][tabmax] = 0.0;
			}
		}
	}

	if (tstr) { free(tstr); tstr = NULL; }
}

// 1. Wochentag direkt mathematisch berechnen (Sakamoto's Algorithmus)
// Liefert: 0 = Sonntag, 1 = Montag, ..., 6 = Samstag
int SP_class::get_weekday(int y, int m, int d) {
	static int t[] = { 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 };
	if (m < 3) y -= 1;
	return (y + y / 4 - y / 100 + y / 400 + t[m - 1] + d) % 7;
}

double SP_class::GetActualOutflow(double ***fluss,AL_class *Abstraction,int ALs,double Qin,double InV,double hour_offset,int p_year,int p_month,int p_day,int p_hour,int p_minute,double dOldOutflow,double iMinTimeStep)
{

	// ToDo
	// here we can cycle through ALs. If an AL_art==9 is found, then we read the actual value. If this is valid (i.e. no nodata), and the replace-option includes abstraction 0 (i.e. the downstream outflow out of the reservoir itself)
	// then the table lookup can be skipped for each sub-step and instead the read in value can be used as Abakt
	// for all other replace-options, only the first listed abstraction to be replaced will be replaced by the respective AL_art==9 abstratcion
	// example: 
	// 
	// TG 21 (AE = 155.0000000, AErel = 1.0)
	//  from OL 20 (kh = 0.1, kv = 0.4, Bh = 9.7, Bv = 38.7, Th = 0.97, Mh = 40.0, Mv = 15.0, I = 0.0010, L = 11414.2, AE = 32.0000000)
	//   and OL 6 (kh = 0.1, kv = 0.4, Bh = 12.1, Bv = 48.6, Th = 1.21, Mh = 40.0, Mv = 15.0, I = 0.0010, L = 11414.2, AE = 58.7500000)
	//   and SP 3 (file = $outpath//spv__03_sub21.//$code//$year , V0 = 92.2E6, C0 = 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0, dTmin = $dtmin1 )  # Sihlsee C0: initial concentration for all of max. 9 tracers in the reservoir water; optional parameter dTmin is in seconds (60s is default even without this parameter given)
	//   and AL $28(modus = intern_from_file_from_reservoir, externalAbstractionFile = $inpath_hydro//spv__03_sub21.s91 31, replaced_Abstractions = { 0 $1 $2 } ) # that means, that this abstraction is read in from an external file, reading in column 5 and the values will replace the abstraction $1 and set abstraction $2 to 0
	//   and AL $1(modus = intern_with_rule_from_reservoir) # used for residual water --> will be replaced with 0, since $1 is the second replace_Abstractions list entry; the abstraction for the downstream area is used for the read in value here (0 in replace_Abstractions list)
	//   and AL $2(modus = intern_with_rule_from_reservoir) # used for residual water --> will be replaced with 0, since $1 is the 3rd replace_Abstractions list entry
	//   and AL $29(modus = intern_from_file_from_reservoir, externalAbstractionFile = $inpath_hydro//givenAbstractionsSihlsee2.dat 5, replaced_Abstractions = { $21 $22 } ) # that means, that this abstraction is read in from an external file, reading in 5 and the values will replace the abstraction $21 and set abstraction $22 to 0
	//   and AL $21(modus = intern_with_rule_from_reservoir) # used for hydro power generation weekdays Sihlsee --> will be replaced with the read in abstraction since it is the first entry in the replace_Abstractions list
	//   and AL $22(modus = intern_with_rule_from_reservoir) # used for hydro power generation weekends Sihlsee --> will be replaced with 0, since $1 is the second replace_Abstractions list entry

	bool bUseExternalAbstraction = false;
	double dExternalAbstraction = 0;
	for (int i = 1; i <= ALs; i++)
	{
		if (Abstraction[i].AL_art == 9)
		{
			for (int j = 0; j < Abstraction[i].vrefAL_localIDs.size(); j++)
			{
				if (Abstraction[i].vrefAL_localIDs[j] == 0)
				{
					bUseExternalAbstraction = true;
					dExternalAbstraction = Abstraction[i].FB->GetActValue(Abstraction[i].column-1); // FB class ignores date and time columns and starts with index 0, so the real column is the one from the control file minus 1 
					break;
				}
			}
		}
		if (bUseExternalAbstraction)
		{
			break;
		}
	}

	double dAbTotal = 0; // this is the return value of the method

	if (!useHydraulicConnection)
	{
		year = p_year;
		month = p_month;
		day = p_day;
		hour = p_hour;
		minute = p_minute;

		int actJD = 0;
		for (int i = 1; i < month; i++) actJD += mon[i];
		actJD += day;
		int rule = 1;
		while ((actJD > JD[rule]) && (rule < rulecnt)) rule++;
		double WeightNextRule = 1;
		double WeightLastRule = 0;
		double dLastAbakt;
		int nextJD = JD[rule];
		int LastJD;
		if (rule > 1)
		{
			LastJD = JD[rule - 1];
			WeightNextRule = double(actJD - LastJD) / double(nextJD - LastJD);
			WeightLastRule = 1.0 - WeightNextRule;
		}
		else
		{
			WeightNextRule = 1.0;
			WeightLastRule = 0.0;
		}

		// optional parameter UseRunningAvg defines the number of hours the running average fo reserboir content should be calculated from. So we need to convert that hours into numbe rof time steps. 
		// which is equal to UseRunningAvg when time step is one hour, but smalle when time step is larger and larger when time step is smaller
		unsigned long numVhistSteps = (60.0*UseRunningAvg) / InV;
		if (numVhistSteps <= 0)    {
			numVhistSteps = 1; // at least one step should be stored in the queue for the running average of volumes, even if the time step if longer than one day
		}
//		else if (numVhistSteps >144) { // 10 minutes is shortest time step for using 24 hour running average, all shorter itme steps will use shorter running average windows
//			numVhistSteps = 144;
//		}
		unsigned long steps = 60.0*InV / dTmin;
		//double dAbTotal = 0;
		double *dAbFromResTotal;
		double *dAbFromResAct;
		dAbFromResTotal = new double[ALs + 1];
		dAbFromResAct = new double[ALs + 1];
		for (int i = 0; i <= ALs; i++)
		{
			dAbFromResTotal[i] = 0.0;
			dAbFromResAct[i] = 0.0;
		}
		// internally used volume for subtracting sub time step outflows in order to account for changing outflow rates if the volume decreases considerably
		// the real balancing will be done outside of this function when returning the outflow as result
		double Vnew = Vakt[0];

		V_hist.push_back(Vakt[0]);
		if (V_hist.size() > numVhistSteps)
		{
			V_hist.pop_front();
		}
		// computing running average rom all values in the list, could also be done by using a stored value V, where the first list entry is subtarcted from and the last list entry is added to.
		double VrunningAvg = 0;
		double cnt = 0;
		if (V_hist.size() > 1)
		{
			for (std::list<double>::iterator iter = V_hist.begin(); iter != V_hist.end(); iter++)
			{
				VrunningAvg += (*iter);
				cnt += 1.0;
			}
			VrunningAvg /= cnt;
		}
		else
		{
			VrunningAvg = Vnew;
		}
		
		bool bApplyRule = false;

//		if (month > 0) // valid month, day, hour and minute must be provided, otherwise the weekday cannot be calculated
//		{
//			time_t acttime;
//			tm* tm1 = new(tm);
//			tm1->tm_hour = hour - (int)hour_offset; // 0...23
////			tm1->tm_hour = hour; // 0...23
//			tm1->tm_min = useMinutes?minute:0; // 0..60 
//			tm1->tm_sec = 0; // 0..60
//			tm1->tm_mon = month - 1;  // 0..11
//			tm1->tm_year = year - 1900; // e.g. 2010 = 110
//			tm1->tm_mday = day;
//			tm1->tm_wday = 0;
//			tm1->tm_yday = 0;
//			tm1->tm_isdst = 0;
//
//#ifdef WIN64
//			acttime = _mkgmtime(tm1);
//#else
//			acttime = timegm(tm1);
//#endif
//			//acttime = mktime(&tm1);
//
//			std::stringstream s1, s2;
//			if (tm1->tm_wday > 0)
//			{
//				s1 << tm1->tm_wday;
//			}
//			else
//			{
//				s1 << 7; // since Sunday is not stored as 7 but as 0 (days since Sunday=0)
//			}
//			s2 << WeekDays[rule];
//			std::string wday;
//			std::string weekdays;
//			s1 >> wday;
//			s2 >> weekdays;
//			if (WeekDays[rule] == -1)
//			{
//				bApplyRule = true; // if no weekdays ared efined
//			}
//			else
//			{
//				if (weekdays.find(wday.c_str()) != std::string::npos)
//				{
//					if (
//							((start_hour[rule] <= stop_hour[rule]) && ((hour >= start_hour[rule]) && (hour <= stop_hour[rule]))) // conventional, e.g. start=7, stop=15: all hours between 7 and 15 will apply the rule
//							||
//							((start_hour[rule] > stop_hour[rule]) && ((hour >= start_hour[rule]) || (hour <= stop_hour[rule]))) // over bight, e.g. start = 16, stop = 6: all hours between 16 and 24 and all hours between 1 and 6 will apply
//					   )
//					{
//						bApplyRule = true;
//					}
//				}
//			}
//			delete tm1;
//		}

		// regard WeekDays, start_hour, stop_hour, start_hour2 and stop_hour2
		// 	   new method which is much faster (found by asking AI to optimize the above commented section)
		if (month > 0) // valid month, day, hour and minute must be provided, otherwise the weekday cannot be calculated
		{
			if (WeekDays[rule] == -1)
			{
				bApplyRule = true; // if no weekdays ared defined
			}
			else
			{
				// get week day without timoe consuming mkdate function calls
				int wday = get_weekday(year, month, day);

				// ISO-week day (1 = Monday, ..., 7 = Sunday)
				int iso_wday = (wday == 0) ? 7 : wday;

				// much more efficient week day check (instead of using stringstream)
				// precondition here: WeekDays[rule] contains e.g. "12345" for Mo-Fr.
				// If WeekDays is a std::string, we use char-searching instead of ttring-conversion
				char wday_char = '0' + iso_wday;
				const std::string& rule_days = std::to_string(WeekDays[rule]);

				if (rule_days.find(wday_char) != std::string::npos) {
					// 4. Stunden-Regel (Bleibt performant, da reine CPU-Register-Vergleiche)
					int start = start_hour[rule];
					int stop = stop_hour[rule];
					int start2 = start_hour2[rule];
					int stop2 = stop_hour2[rule];

					if (
						(start <= stop && hour >= start && hour <= stop) ||
						(start > stop && (hour >= start || hour <= stop)) ||
						(start2 <= stop2 && hour >= start2 && hour <= stop2) ||
						(start2 > stop2 && (hour >= start2 || hour <= stop2))
						)
					{
						bApplyRule = true;
					}
				}
			}
		}


		//	if (year>0) (*WSstream) << "\n" << this->SP_code <<": ";
		// do the balancing for each subtimestep
		for (unsigned long st = 1; st <= steps; st++)
		{
			double dAbakt = 0;
			if (year > 0)
			{
				Vnew += Qin * dTmin;
			}
			// Ermitteln der aktuellen Speicherabgabe 
			if (!bUseExternalAbstraction)
			{
				if (bReverseLookUp == false) // old, default version with a table that has volumes as index values and fluxes as table values
				{
					// check if the reservoir is more than full. If so, do not use any restricted rule
					if ((Vnew) > ARTabIndex[tabmax]) // Vnew is still computed with the actual volume and inflow. If, for any reason, the volume goes beyond the reservoirs maximum capacity, the outflow should not be calculated by any rule but by the maximum flow possible, i.e. like inflow to protect the reservoir from becoming numerically unstable or irrationale
					{
						//		dAbakt=Qin;
						double dWlevel = GetWaterLevelByVolume(Vnew); // get water level for the volume plus inflow
						// extrapolate from the W-Q-Table
						dAbakt = ARTabValue[rule][tabmax] +
							(Vnew - ARTabIndex[tabmax]) * ARChangeRate[rule][tabmax - 1];
						if (rule > 1) // interpolate between extrapolated abstraction from last ruel and actual rule
						{
							dLastAbakt = ARTabValue[rule - 1][tabmax] +
								(Vnew - ARTabIndex[tabmax]) * ARChangeRate[rule - 1][tabmax - 1];
							// interpolating between abstractions in time
							dAbakt = dAbakt * WeightNextRule + dLastAbakt * WeightLastRule;
						}
						// waterlevel way beyound the controllable level (i.e. more then 0.5m above any table entry for the volume-outflow-relation)
						if (lake_active)
						{
							if ((dWlevel > (WaterLevel[tabmax] + 0.5)) && (dAbakt < Qin))
							{
								dAbakt = Qin;
							}
						}
						else // conventional lake model: always set outflow to inflow if maximum volume is reached
						{
							dAbakt = Qin;
						}
					}
					else
					{
						for (int j = tabmax - 1; j >= 1; j--)
						{
							double Vtemp = (UseRunningAvg > 1) ? VrunningAvg : Vnew;

							if (Vtemp > ARTabIndex[j])
							{
								dAbakt = ARTabValue[rule][j] +
									(Vtemp - ARTabIndex[j]) * ARChangeRate[rule][j];
								if (rule > 1)
								{
									dLastAbakt = ARTabValue[rule - 1][j] +
										(Vtemp - ARTabIndex[j]) * ARChangeRate[rule - 1][j];
									// interpolating between abstractions in time
									dAbakt = dAbakt * WeightNextRule + dLastAbakt * WeightLastRule;
								}
								break;
							}
						}
						// prevent abstraction from beeing taken from reservoir when the rule is not aktiv regarding the schedule
						if (bApplyRule == false)
						{
							dAbakt = 0;
						}
					}
				}
				else
				{ // reverse LookUp, contains volumes as table values and fluxes as table indices

					for (int j = tabmax - 1; j >= 1; j--)
					{
						double ARTabValueInterpol = (rule > 1) ?
							ARTabValue[rule][j] * WeightNextRule + ARTabValue[rule - 1][j] * WeightLastRule :
							ARTabValue[rule][j];

						double Vtemp = (UseRunningAvg > 1) ? VrunningAvg : Vnew;

						if (Vtemp > ARTabValueInterpol)
						{
							double ARChangeRateInterpol = (rule > 1) ?
								ARChangeRate[rule][j] * WeightNextRule + ARChangeRate[rule - 1][j] * WeightLastRule :
								ARChangeRate[rule][j];

							dAbakt = ARTabIndex[j] +
								(Vtemp - ARTabValueInterpol) * ARChangeRateInterpol;

							if ((Vnew > ARTabValueInterpol) || (Vtemp > ARTabValueInterpol))
							{
								if (j == tabmax)
								{
									double maxV = max(Vnew, Vtemp);
									double dWlevel = GetWaterLevelByVolume(maxV); // get water level for the volume plus inflow
									// waterlevel way beyound the controllable level
									if (lake_active)
									{
										if ((dWlevel > (WaterLevel[tabmax] + 0.5)) && (dAbakt < Qin))
										{
											dAbakt = Qin;
										}
									}
									else // conventional lake model: always set outflow to inflow if maximum volume is reached
									{
										dAbakt = Qin;
									}
								}
							}
							break;
						}
					}
					// prevent abstraction from beeing taken from reservoir when the rule is not aktiv regarding the schedule
					if (bApplyRule == false)
					{
						dAbakt = 0;
					}
				}
			}
			else
			{
				dAbakt = dExternalAbstraction;
			}

			// if target capacity restricts the outflow
			if (targetCapacity[rule] > 0)
			{
				if (dAbakt > targetCapacity[rule])
				{
					dAbakt = targetCapacity[rule];
				}
			}
			//		  if (year>0) (*WSstream) << dAbakt << " ";
			if (year > 0)
			{
				Vnew -= (dAbakt*dTmin); // reducing the value by the amount of outflow during this interval
			}
			if (Vnew < 0) // if the actual abstraction leads to an empty reservoir, numerically it may even have a negative content, so the abstraction must be reduced according to give exactly a zero content
			{
				dAbakt -= (-Vnew) / dTmin; // Vnew is negativ, so (-Vnew) will be positive in order to really reduce Abakt 
			}
			if (dAbakt < 0) dAbakt = 0;
			dAbTotal = dAbTotal + dAbakt / steps; // result is in m^3/s

			// now cycling through the abstractions from reservoir
			for (int i = 1; i <= ALs; i++)
			{
				// Todo 
				// if one or more abstractions of type AL_art == 9 are defined and the actual value from the input file is valid, then the first Abstrations with AL_art == 8 that was listed in the 
				// replace-option of the respective AL_art == 9 abstractions will be replaced by that read-in value
				if (Abstraction[i].AL_art == 8)
				{
					// abstraction are in m^3/s, but dAbFromResTotal is in m^3 (multiplying abstraction withdTmin here!)
					double Vtemp = (UseRunningAvg > 1) ? VrunningAvg : Vnew;
					dAbFromResAct[i] = Abstraction[i].GetActualAbstractionFromReservoir(fluss, Vtemp, year, month, day, hour,minute)*dTmin;

					// now check if this abstraction needs to be replaced by an external value (from an abstractionof type AL_art == 9 (modus = intern_from_file_from_reservoir)
					for (int j = 1; j <= ALs; j++)
					{
						if (Abstraction[j].AL_art == 9)
						{
							for (int k = 0; k < Abstraction[j].vrefAL_localIDs.size(); k++)
							{
								if (Abstraction[j].vrefAL_localIDs[k] == i) // there is an entry in the other abstractions reference vector holding the index of the actual abstraction, so the actual sbatractions value needs to be replaced by either 0 or the other abstractions value
								{
									double actExtAbstractionValue = Abstraction[j].FB->GetActValue(Abstraction[j].column - 1);
									if (actExtAbstractionValue != nodata)
									{
										if (k == 0) // it is the first entry of the reference vector, so the others abstraction value is taken to replace the actual abstractions value with
										{
											dAbFromResAct[i] = Abstraction[j].FB->GetActValue(Abstraction[j].column - 1) * dTmin;
											// cloumn-1-4 means: since the control file specifies the absolute column including date columns and also starting with index 1 and 
											// but the FB class starts to count at 0 and ommits the date columns, we need to adjust for this difference here
										}
										else
										{
											dAbFromResAct[i] = 0; // all other abstractions that should be replaced by the external abstraction will be set to 0 because the first abstractions in the replace list gets all
										}
									}
								}
							}
						}
					}


					Vnew -= dAbFromResAct[i]; // volume change in m^3
					if (Vnew < 0.9*Vtemp)
					{
						dAbFromResAct[i] -= (0.9*Vtemp-Vnew); // the negative Vnew is added to actual abstraction, reducing it effectively
					}
					if (dAbFromResTotal[i] < 0) dAbFromResTotal[i] = 0;
					dAbFromResTotal[i] += dAbFromResAct[i]; // in m^3
				}
			}

		}

		// setting the Abakt values of each abstraction (in m3/s)
		for (int i = 1; i <= ALs; i++)
		{
			if (Abstraction[i].AL_art == 8)
			{
				Abstraction[i].Abakt[0] = dAbFromResTotal[i] / (InV*60.0); // absolute abstraction value is now converted into m^3/s
			}
		}

		//	if (year>0) (*WSstream) << ">> mean: "<<dAbTotal<< "[end-entry]";
		delete[] dAbFromResTotal;
		delete[] dAbFromResAct;
//		return dAbTotal; // return actual discharge into downstream water in m^3/s; this is not the total of downstram discharge and abstraction! Abstractions have to be balanced later!
	}
	else
	{

		unsigned long steps = 60.0*InV / dTmin;
		//double dAbTotal = 0;
		double *dAbFromResTotal;
		double *dAbFromResAct;
		dAbFromResTotal = new double[ALs + 1];
		dAbFromResAct = new double[ALs + 1];
		for (int i = 0; i <= ALs; i++)
		{
			dAbFromResTotal[i] = 0.0;
			dAbFromResAct[i] = 0.0;
		}
		// internally used volume for subtracting sub time step outflows in order to account for changing outflow rates if the volume decreases considerably
		// the real balancing will be done outside of this function when returning the outflow as result
		double Vnew = Vakt[0];
		double localActWaterLevel = actWaterLevel;
		TargetWaterLevel = *(mapReservoirLevels->find(TargetReservoir)->second[0]);


		if (*(mapReservoirLevels->find(TargetReservoir)->second[2]) != 0)
		{
			double tempDVm3s = (*(mapReservoirLevels->find(TargetReservoir)->second[1]));
			double tempArea = (*(mapReservoirLevels->find(TargetReservoir)->second[2]));
				TargetWaterLevelChangeRate = (tempDVm3s * dTmin) / tempArea;
		}
		else
		{
			TargetWaterLevelChangeRate = 0;
		}

		//	if (year>0) (*WSstream) << "\n" << this->SP_code <<": ";
		// do the balancing for each subtimestep
		for (unsigned long st = 1; st <= steps; st++)
		{
			double dAbakt = 0;

			Vnew += Qin*dTmin;
			if (this->AREAm2 > 0)
			{
				localActWaterLevel += (Qin*dTmin) / this->AREAm2;
			}
			// Ermitteln der aktuellen Speicherabgabe 

			// calculate flux between the lakes based on the hydraulic properties given in the AbstractionRule definition
// the consideration of volume change in the targes reservoir is sometimes causing oscillations
//			TargetWaterLevel += TargetWaterLevelChangeRate/2.0;
			double actSlope = 0;
			double actChannelDepth = 0;
			double actHydraulicRadius = 0;
			double actChannelArea = 0;
			if ((localActWaterLevel > ChannelBottom) && (TargetWaterLevel > ChannelBottom))
			{
				actSlope = (localActWaterLevel - TargetWaterLevel) / ChannelLength;
				if (actSlope > 0)
				{
					actChannelDepth = localActWaterLevel - ChannelBottom;
				}
				else
				{
					actChannelDepth = TargetWaterLevel - ChannelBottom;
				}
			}
			else if ((localActWaterLevel > ChannelBottom) && (TargetWaterLevel <= ChannelBottom))
			{
				actSlope = (localActWaterLevel - ChannelBottom) / ChannelLength;
				actChannelDepth = localActWaterLevel - ChannelBottom;
			}
			else if ((localActWaterLevel <= ChannelBottom) && (TargetWaterLevel > ChannelBottom))
			{
				actSlope = (ChannelBottom - TargetWaterLevel) / ChannelLength;
				actChannelDepth = TargetWaterLevel - ChannelBottom;
			}

			if ((actChannelDepth > 0) && (actSlope != 0))
			{
				actHydraulicRadius = (actChannelDepth*ChannelWidth) / (ChannelWidth + 2.0*actChannelDepth);
				actChannelArea = ChannelWidth*actChannelDepth;
				// flux equation after Manning/Gauckler/Strickler Q = M*A*R^2/3*S^1/2
				if (actSlope > 0)
				{
					dAbakt = RoughnessCoeff*actChannelArea*pow(actHydraulicRadius, double(2.0 / 3.0))*sqrt(actSlope);
				}
				else
				{
					dAbakt = RoughnessCoeff*actChannelArea*pow(actHydraulicRadius, double(2.0 / 3.0))*sqrt(-actSlope);
					dAbakt = -dAbakt;
				}
			}
			else
			{
				dAbakt = 0;
			}

			// replace dAbakt with externally read value if AL_art==9 and 0 is part of the replaceALs-list (was evaluated at the beginning of this routine
			if (bUseExternalAbstraction)
			{
				dAbakt = dExternalAbstraction;
			}

			if (this->AREAm2 > 0)
			{
				localActWaterLevel -= (dAbakt*dTmin) / this->AREAm2;
			}
			if ((localActWaterLevel < ChannelBottom)&&(actWaterLevel > ChannelBottom))
			{
				dAbakt -= this->AREAm2*(ChannelBottom - localActWaterLevel) / dTmin;
				localActWaterLevel = ChannelBottom;
			}
			
			//		  if (year>0) (*WSstream) << dAbakt << " ";

			Vnew -= (dAbakt*dTmin); // reducing the value by the amount of outflow during this interval
			if (Vnew < 0) // if the actual abstraction leads to an empty reservoir, numerically it may even have a negative content, so the abstraction must be reduced according to give exactly a zero content
			{
				dAbakt -= (-Vnew) / dTmin; // Vnew is negativ, so (-Vnew) will be positive in order to really reduce Abakt 
			}
			//if (dAbakt < 0) dAbakt = 0;
			dAbTotal = dAbTotal + dAbakt / steps; // result is in m^3/s

			// now cycling through the abstractions from reservoir
			for (int i = 1; i <= ALs; i++)
			{
				// Todo 
				// if one or more abstractions of type AL_art == 9 are defined and the actual value from the input file is valid, then the first Abstrations with AL_art == 8 that was listed in the 
				// replace-option of the respective AL_art == 9 abstractions will be replaced by that read-in value
				if (Abstraction[i].AL_art == 8)
				{
					// abstraction are in m^3/s, but dAbFromResTotal is in m^3 (multiplying abstraction withdTmin here!)
					dAbFromResAct[i] = Abstraction[i].GetActualAbstractionFromReservoir(fluss, Vnew, year, month, day, hour,minute)*dTmin;

					// now check if this abstraction needs to be replaced by an external value (from an abstractionof type AL_art == 9 (modus = intern_from_file_from_reservoir)
					for (int j = 1; j <= ALs; j++)
					{
						if (Abstraction[j].AL_art == 9)
						{
							for (int k = 0; k < Abstraction[j].vrefAL_localIDs.size(); k++)
							{
								if (Abstraction[j].vrefAL_localIDs[k] == i) // there is an entry in the other abstractions reference vector holding the index of the actual abstraction, so the actual sbatractions value needs to be replaced by either 0 or the other abstractions value
								{
									double actExtAbstractionValue = Abstraction[j].FB->GetActValue(Abstraction[j].column - 1);
									if (actExtAbstractionValue != nodata)
									{
										if (k == 0) // it is the first entry of the reference vector, so the others abstraction value is taken to replace the actual abstractions value with
										{
											dAbFromResAct[i] = Abstraction[j].FB->GetActValue(Abstraction[j].column - 1) * dTmin;
											// cloumn-1-4 means: since the control file specifies the absolute column including date columns and also starting with index 1 and 
											// but the FB class starts to count at 0 and ommits the date columns, we need to adjust for this difference here
										}
										else
										{
											dAbFromResAct[i] = 0; // all other abstractions that should be replaced by the external abstraction will be set to 0 because the first abstractions in the replace list gets all
										}
									}
								}
							}
						}
					}

					Vnew -= dAbFromResAct[i]; // volume change in m^3
					if (Vnew < 0)
					{
						dAbFromResAct[i] -= (-Vnew); // the negative Vnew is added to actual abstraction, reducing it effectively
					}
					if (dAbFromResTotal[i] < 0) dAbFromResTotal[i] = 0;
					dAbFromResTotal[i] += dAbFromResAct[i]; // in m^3
				}
			}

		}

		// setting the Abakt values of each abstraction (in m3/s)
		for (int i = 1; i <= ALs; i++)
		{
			if (Abstraction[i].AL_art == 8)
			{
				Abstraction[i].Abakt[0] = dAbFromResTotal[i] / (InV*60.0); // absolute abstraction value is now converted into m^3/s
			}
		}

		//	if (year>0) (*WSstream) << ">> mean: "<<dAbTotal<< "[end-entry]";
		delete[] dAbFromResTotal;
		delete[] dAbFromResAct;
//		return dAbTotal; // return actual discharge into downstream water in m^3/s; this is not the total of downstram discharge and abstraction! Abstractions have to be balanced later!
	}

	// read next value from AL_art == 9 abstractions for next global time step
	if (year > 0)
	{
		for (int i = 1; i <= ALs; i++)
		{
			if (Abstraction[i].AL_art == 9)
			{
				Abstraction[i].FB->NextRecord();
			}
		}
	}

	return dAbTotal; // return actual discharge into downstream water in m^3/s; this is not the total of downstram discharge and abstraction! Abstractions have to be balanced later!

}


void SP_class::findeabschnitt(FILE *steuerdatei,char *Abschnittsname)
{
	char *tempstr,*satz,*ttempstr,*tsatz;

	tempstr=(char*)calloc(_stringsize_,sizeof(char)); ttempstr=tempstr;
	satz=(char*)calloc(_stringsize_,sizeof(char)); tsatz=satz;

	fseek(steuerdatei,0,0);
	strcpy(satz,"");
	while ((strcmp(satz,Abschnittsname)!=0)&&(!feof(steuerdatei))&&tempstr)
	{
		fscanf(steuerdatei,"%s",satz);
		fgets(tempstr,_stringsize_,steuerdatei);
	}
	if (feof(steuerdatei))
	{
		(*WSstream) << "FATAL: section "<<Abschnittsname<<" not found in control file"<<"[end-entry]";
		throw(34);
	}
	tempstr=ttempstr; if (tempstr) {free(tempstr); tempstr=NULL;}
	satz=tsatz;       if (satz)    {free(satz);    satz=NULL;	}
}

void SP_class::initialize(wasim::mystream *p_WSstream,tracer_info_type *p_tracer_info,int p_tracer_count,
int p_TG,gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,
firstlast *p_FirstLastCol ,int p_zonmax,int p_zonmaxwert,
int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata, int p_lake_active,
char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,double hour_offset,int timestep)
{
	WSstream = p_WSstream;
	tracer_info = p_tracer_info;
	tracer_count = p_tracer_count;
	TG = p_TG;
	gridinfo = p_gridinfo;
	zarray = p_zarray;
	zone = p_zone;
	FirstLastCol = p_FirstLastCol;
	zonmax = p_zonmax;
	zonmaxwert = p_zonmaxwert;
	ncols = p_ncols;
	nrows = p_nrows;
	mpistartrow = p_mpistartrow;
	mpiendrow = p_mpiendrow;
	xll = p_xll;
	yll = p_yll;
	cellsize = p_cellsize;
	nodata = p_nodata;
	lake_active = p_lake_active;
	Vakt=(double*)calloc(tracer_count+1,sizeof(double));
	Vold=(double*)calloc(tracer_count+1,sizeof(double));
	Vtemp=(double*)calloc(tracer_count+1,sizeof(double));
	conc=(double****)calloc(tracer_count+1,sizeof(double***));
	DefaultInputDirectory=p_DefaultInputDirectory;
	DefaultOutputDirectory=p_DefaultOutputDirectory;
	this->hour_offset = hour_offset;
	this->timestep = timestep;
	if (this->timestep < 60 && this->timestep>0)
	{
		useMinutes = true;
	}
}

// now, we can develop the depth-volume-relation into a lookup-table
// and do some other performance optimizing things 
void SP_class::initialize_lake()
{
	// do this step when LAKE grid is there and [lake_model] is activated, only
	if 	(!lake_active) return;

	// looking for the start end end row and column in the LAKE grid for the actual lake
	startcol=ncols+1;
	startrow=nrows+1;
	endcol=0;
	endrow=0;
	double DHM_ini;
	bool dhm_initialized=false;

	for (int i=1;i<=nrows;i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if (zn[i][j]!=nodata)
			{
				if (LAKE[i][j]==SP_code)
				{
					if (!dhm_initialized)
					{
						DHM_ini=dhm[i][j]+MAXPOND[i][j]/1000.0;
						dhm_initialized=true;
					}
					if (fabs(DHM_ini-(dhm[i][j]+MAXPOND[i][j]/1000.0))>0.01)
					{
						(*WSstream) << "FATAL: \n============== ERROR ====================\n"
						<< "Error in initialization for lake "<<SP_code<<"."
						<< "All Cells for a lake must have the same Altitude value for DEM + MAXPOND\n"
						<< "Value at row "<<i<<", col "<<j<<" is "<<dhm[i][j]<<"\n"
						<< "but first value of dhm for this lake was "<<DHM_ini<<"\n"
						<< "\n\n=======================================\n"<< "[end-entry]";
						throw(213);
					}

					if (zarray[TG]!=zarray[(int)zn[i][j]])
					{
						(*WSstream) << "FATAL:  error in routing description or LAKE grid:\n"
						<< " Lake "<<SP_code<<" must be located completely within subbasin "<<TG<<"." <<"[end-entry]";
						throw(224);
					}
					if (i<startrow) startrow=i;
					if (i>endrow) endrow=i;
					if (j<startcol) startcol=j;
					if (j>endcol) endcol=j;
				}
			}
		}
	}
	if (endcol==0) endcol=startcol;
	if (endrow==0) endrow=startrow;
	// now creating the depth-volume-relationship
	double min_dhm=9999;
	double max_depth=0;
	double max_WaterLevel=0;
	double min_WaterLevel=0;
	double act_depth=0;
	double depth_range;
	
	int steps1=0;
	int steps2=0;
	int steps3=0;
	
	double stepThick1 = 0.01;
	double stepThick2 = 0.02;
	double stepThick3 = 0.05;

	double LayerMaxThick1 = 4.0;
	double LayerMaxThick2 = 10.0;

	double upper_layer,mid_layer,lower_layer; // actual thickness of the layers
	// 1. looking for the lowest value in the dhm within the lake
	// 2. looking for the maximum value of dhm+maxpond -> the difference of max(dhm+maxpond)-min(dhm) will be the range of the depth
	for (int i=startrow;(i<=endrow) && (i<=nrows);i++)
	{
		for (int j=startcol;(j<=endcol)&&(j<ncols);j++)
		{
			if (LAKE[i][j]==SP_code)
			{
				if (dhm[i][j]<min_dhm)
				{
					min_dhm = dhm[i][j];
				}
				if ((dhm[i][j] + MAXPOND[i][j] / 1000.0) > max_WaterLevel)
				{
					max_WaterLevel = dhm[i][j] + MAXPOND[i][j] / 1000.0;
				}
			}
		}
	}
	min_WaterLevel=min_dhm;
	// 3. dividing the depth-range into steps of 1cm starting from max-depth for the first 2m, then 10cm for the next 10m and then 0.5m for the rest of the range
	depth_range=max_WaterLevel-min_WaterLevel;
	if (depth_range >= LayerMaxThick1)
	{
		upper_layer = LayerMaxThick1;
	}
	else
	{
		upper_layer = depth_range;
	}

	if (depth_range - upper_layer >= LayerMaxThick2)
	{
		mid_layer = LayerMaxThick2;
	}
	else
	{
		mid_layer = depth_range - upper_layer;
	}

	lower_layer = depth_range-upper_layer-mid_layer;
	
	steps1 = (int)(upper_layer/stepThick1); // 1cm levels in the uppermost 2m of possible water levels
	steps2 = (int)(mid_layer/stepThick2); // 2cm levels in the next 10m of possible water levels
	steps3 = (int)ceil(lower_layer/stepThick3); // 5cm levels in the rest down to the bottom of the LAKE
	if (steps1==0) steps1=1; // at minimum, 1 step is required in order to get a table with at least 2 entries
	W_V_total_steps = steps1+steps2+steps3+1;
	// allocating memory for the table
	WaterLevel = new double[W_V_total_steps+2];
	Volume = new double[W_V_total_steps+2];
	// 4. for each depth
	//  4.1 get the corresponding altitude of the water level for the actual depth as wlevel
	//  4.2 running through all cells of the lake, 
	//   4.2.1 calculating the difference between the wlevel and the dhm-value
	//   4.2.2 calculating the corresponding water volume, summing up the value in a total sum
	//  4.3 make the entry into the depth-volume-table
	// 0.01m-steps for the uppermost layer of the water body
	for (int z=0;z<steps1;z++)
	{
		WaterLevel[z]=max_WaterLevel-z*stepThick1;
		double V=0;
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if (LAKE[i][j]==SP_code)
				{
					if (WaterLevel[z]>dhm[i][j])
					{
						V+=cellsize*cellsize*(WaterLevel[z]-dhm[i][j]);
					}
				}
			}
		}
		Volume[z]=V;
	}
	// 0.1m-steps for the middle layer of the water body
	for (int z=0;z<steps2;z++)
	{
		WaterLevel[steps1+z]=max_WaterLevel-upper_layer-z*stepThick2;
		double V=0;
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if (LAKE[i][j]==SP_code)
				{
					if (WaterLevel[steps1+z]>dhm[i][j])
					{
						V+=cellsize*cellsize*(WaterLevel[steps1+z]-dhm[i][j]);
					}
				}
			}
		}
		Volume[steps1+z]=V;
	}
	// 0.5m-steps for the lower layer of the water body
	for (int z=0;z<steps3;z++)
	{
		WaterLevel[steps1+steps2+z]=max_WaterLevel-upper_layer-mid_layer-z*stepThick3;
		double V=0;
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if (LAKE[i][j]==SP_code)
				{
					if (WaterLevel[steps1+steps2+z]>dhm[i][j])
					{
						V+=cellsize*cellsize*(WaterLevel[steps1+steps2+z]-dhm[i][j]);
					}
				}
			}
		}
		Volume[steps1+steps2+z]=V;
	}
	Volume[W_V_total_steps-1]=0;
	WaterLevel[W_V_total_steps-1]=min_WaterLevel;

	// 5. when ready, write the table to the logfile -> at least for the first time.
	std::ofstream W_V_file;
	char *filename;
	filename=new char[255];
	sprintf(filename,"W_V_Table_%02d.dat",SP_code);
	W_V_file.open(filename,std::ios::out);
	if (!W_V_file.good() || W_V_file.eof())
	{
		(*WSstream) << "FATAL: error opening " << filename << " for writing the waterlevel-volume-relation to disk."<<"[end-entry]";
		throw (225);
	}
	else
	{
		W_V_file << "waterlevel-volume-table lake "<<SP_code<<"\n";
		W_V_file << "Level[m a.s.l.]\tVolume[m^3]\n";
		for (int i=0;i<W_V_total_steps;i++)
		{
			W_V_file << WaterLevel[i] << "\t" << Volume[i] << "\n";
		}
		W_V_file.close();
	}
	delete[] filename;

}

double SP_class::GetWaterLevelByVolume(double V)
{
	if (lake_active)
	{
		// 6. now, by looking for an entry in the WaterDepth-array, the corresponding Volume may be get from the same index in the Volume-array
		//    vice versa, any WaterDepth may be get by looking into the Volume-array for the proper index and get the WaterDepth from the same index in the WaterDepth-array
		int k=0;
		while ((Volume[k+1]>V)&&(k<W_V_total_steps-2)) k++;
		return ((Volume[k]-V)*WaterLevel[k+1]+(V-Volume[k+1])*WaterLevel[k])/(Volume[k]-Volume[k+1]);

	}
	else 
	return 0;
}

double SP_class::GetVolumeByWaterLevel(double W)
{
	if (lake_active)
	{
		int k=0;
		while ((WaterLevel[k+1]>W)&&(k<W_V_total_steps-2)) k++;
		return ((WaterLevel[k]-W)*Volume[k+1]+(W-WaterLevel[k+1])*Volume[k])/(WaterLevel[k]-WaterLevel[k+1]);

	}
	else
	return 0;
}

// V has the type double[tracer_count+1] to hold all tracer masses
// result: a vector with the total Volume and each tracer concentrations as medium value for all LAKE-cells whith more than 0mm water depth (in POND)
// this function must be called at the beginning of the lake/reservoir handling to redistribute the water of the POND cells onto all valid POND cells
// for the actual water table and to mix all tracers in the lake with eventually changed values in the single POND-cells.
// vector V may then be copied into Vakt, and when balancing the Volumes, the difference between V and Vakt(old-Values before routing) may be compared
// but the old Vakt must be stored somewehere alse in order to be able to make a correct balance (e.g. Vtemp)
// after routing, the next function call is: SetPonds(), which redistributes the volumes and tracer concentrations across all POND cells
// Differences between Vold and Vakt are an indicator for the loss of water by evaporation and groundwater fluxes, resp.
double SP_class::GetVolumeFromPonds(double *V)
{
	double cnt = 0;
	if (lake_active)
	{
		double Vt=0;
		for (int i=0;i<=tracer_count;i++) V[i]=0;
		// Volumes and concentrations from all LAKE-cells which are actually filled with water
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				// a pond belongs to the lake only if there are more than 10mm water in it
				// otherwise, the cell is still a normal cell outside the lake
				if ((LAKE[i][j]==SP_code)&&(POND[i][j]>0.0)) 
				{
					Vt=((POND[i][j]/1000.0)*cellsize*cellsize);
					V[0] += Vt;
					for (int z=1;z<=tracer_count;z++)
					{
						// tracer-mass must be added up for later redistribution
						V[z]+=conc[z][i][j][0]*Vt;
					}
					cnt += 1;
				}
			}
		}
		for (int z=1;z<=tracer_count;z++)
		{
			// calc tracer-concentrations by dividing with the volume
			if (V[0]>0)
			{
				V[z]=V[z]/V[0];
			}
			else
			{
				V[z]=0;
			}
		}

		return cnt;
	}
	else
	return cnt;
}

// setting POND-values and conc-values from Vakt[0..z]-values
double SP_class::SetPonds(double Volume)
{
	double cnt = 0;
	if (lake_active)
	{
		// 7. Then, if the initialization was not read from grid, a function may change all actual POND-Values with the WaterDepth-value for the corresponding cell
		//    but this is yet another function, which will be called externally
		actWaterLevel=GetWaterLevelByVolume(Volume);

		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if (LAKE[i][j]==SP_code)
				{
					if (actWaterLevel>dhm[i][j])
					{
						POND[i][j] = (actWaterLevel - dhm[i][j]) * 1000;
						for (int z=1;z<=tracer_count;z++)
						{
							conc[z][i][j][0]=Vakt[z];
						}
						cnt += 1.0;
					}
					else
					{
						POND[i][j]=0;
						for (int z=1;z<=tracer_count;z++)
						{
							conc[z][i][j][0]=0;
						}

					}
				}
			}
		}
	}
	return cnt;
}

// function is needed for initialization only
void SP_class::SetVaktFromPonds()
{
	if (lake_active)
	{
		double W=0;
		int i_f,j_f;
		bool cellfound=false;
		Vakt[0] = 0;
		double cellarea = cellsize*cellsize/1000.0; // when multiplied with a pond value in mm, the result is in m3 for that cell
		for (int i = startrow; i <= endrow; i++)
		{
			for (int j = startcol; j <= endcol; j++)
			{
				// a pond belongs to the lake only if there are more than 10mm water in it
				// otherwise, the cell is still a normal cell outside the lake
				if ((LAKE[i][j] == SP_code) && (POND[i][j]>0.0))
				{
					Vakt[0] += cellarea*POND[i][j];
					for (int z = 1; z <= tracer_count; z++)
					{
						Vakt[z] += conc[z][i][j][0];
					}
				}
			}
		}
		for (int z = 1; z <= tracer_count; z++)
		{
			if (Vakt[0] > 0)
			{
				Vakt[z] /= Vakt[0];
			}
			else
			{
				Vakt[z] = 0;
			}
		}

// old version which resulted in wrong volumina when a lake had shallower rim cells which were not yet connected to the lake but contained some water in the pond anyway because there was some precip in the last time step
//		// get i,j and WaterLevel from a LAKE-cell which is actually filled with water
//		for (int i = startrow; i <= endrow; i++)
//		{
//			for (int j=startcol;j<=endcol;j++)
//			{
//				// a pond belongs to the lake only if there are more than 10mm water in it
//				// otherwise, the cell is still a normal cell outside the lake
//				if ((LAKE[i][j]==SP_code)&&(POND[i][j]>0.0)) 
//				{
//					W = dhm[i][j]+POND[i][j]/1000.0;
//					i_f=i;
//					j_f=j;
//					cellfound=true;
//					break;
//				}
//			}
//			if (cellfound) break;
//		}
//		if (cellfound)
//		{
//			Vakt[0]=GetVolumeByWaterLevel(W);
//			for (int z=1;z<=tracer_count;z++)
//			{
//				// better for for non evaporating tracers: summing up mass of tracer and divide it by the volume in Vakt[0] to get each Vakt[z]
//				Vakt[z]=conc[z][i_f][j_f][0];
//			}
//		}
//		else
//		{
//			Vakt[0]=0;
//			for (int z=1;z<=tracer_count;z++)
//			{
//				// better for for non evaporating tracers: summing up mass of tracer and divide it by the volume in Vakt[0] to get each Vakt[z]
//				Vakt[z]=0;
//			}
//		}
	}
	else
	{
		return;
	}
}


void SP_class::CalcWaterBalance(double InV,double irr,double abstraction_summ3,double ***fluss,double **Qges,double *flowm3s,AL_class* Abstraction,int ALs,int p_year,int p_month,int p_day,int p_hour,int p_minute,int p_timestepfactor, int actsubtimestep, int numsubtimesteps)
{
	double v1=0,v2=0;
	//  Variables used in SP_class for balance will be initialized with 0
	if (lake_active)
	{
		AREAm2 = 0; // if lake_active=false, AREAm2 was set from routingdefinition after SP X (... A0 = Y...). This Value must not be reset, since it cannot be recalculated from any water covered cells exept when lake_active si true
	}
	ABSTRmm = ABSTRm3s = ABSTRm3 =
	PRECmm  = PRECm3s  = PRECm3  =
	IRRIGmm = IRRIGm3s = IRRIGm3 =
	QINmm   = QINm3s   = QINm3   =
	QOUTmm  = QOUTm3s  = QOUTm3  =
	QGWmm   = QGWm3s   = QGWm3   =
	ETRmm   = ETRm3s   = ETRm3   =
	DVmm    = DVm3s    = DVm3    =
	ERRmm   = ERRm3s   = ERRm3   = 0.0;

	double abstraction_summ3s = 0; // for abstractions from reservoir (except the generic abstraction into downwater), given in m^3/s
	double abstraction_mm = 0;
	
	IRRIGm3  = irr;
	ABSTRm3  = abstraction_summ3;
	QINm3    = flowm3s[0]*(InV*60.0);
	QINm3s   = QINm3/(InV*60.0);

	if (lake_active)
	{
		double cellcnt1 = 0;
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if ((LAKE[i][j]==SP_code)&&(POND[i][j]>0.0)) 
				{
					if (ETR)
					{
						// if interception is modelled, it must be subtracted from ETR, since POND storage in unsatzonmodel is only reduced by ETR (EI was taken from interception storage instead)
						// --> note: EI can only be > 0 when there is a multi layer landuse defined on a lake cell, e.g. for lakes with varying surface area
						ETRm3 += (cellsize * cellsize) / 1000.0 / p_timestepfactor * (EI ? (ETR[i][j] - EI[i][j]) : ETR[i][j]);
					}
					if ((qu)&&(DELTAGW))
					{
						if (fabs(qu[i][j])>1e-10)
							QGWm3  += (cellsize*cellsize)*qu[i][j]*(InV*60.0); // qu is already gioven in m/s, so we have no reason to divide when running in other time steps than unsatzonmodel
						else
							QGWm3  -= (cellsize*cellsize)*DELTAGW[i][j]*(InV*60.0); // DELTAWG is already gioven in m/s, so we have no reason to divide when running in other time steps than unsatzonmodel
					}
					else
					{
						QGWm3=0;
					}
					PRECm3 += (cellsize*cellsize)*PREC[i][j]/1000.0 / p_timestepfactor;
					cellcnt1 += 1.0;
				}
			}
		}
		// if the POND-cells where filled from rain or groundwater saturation, then the water table may not even reach the cells,
		// so we will regard the water from that cells as direct inflow into the lake, resulting in a bit larger volume
		// and in empty pond-cells

		// usual process: ponds were filed in unsatzon with rain and GW-Exchange, emptied by ETR and GWexchng. The total pond for all cells is now distributed obver ll cells of the lake
		// evenly.
		// process when the time step of unsatzon and routing are not identical:
		// POND values will have been updated externally already by ETR, PREC and GWN.
		// in order to get the most likely value for the first sub step, a new temp variable (e.g. deltaPond) must be created that contains the total amout of PREC-ETR-GWN
		// This must be subtracted from the Vtemp value which is calculated as first value hereunder.
		// then, only the fraction that represents the actual substep should be added for each sub step to either Vtemp or POND (or both when using GetVolumeFromPonds and/or SetPonds
		// At the last substep, the full amount of deltaPond was added to Vtemp and all POND cells, thus making sure that the volume of the lake 
		// is smoothly changing also at the next main time step

		double cellcnt2 = this->GetVolumeFromPonds(Vtemp);
		if (numsubtimesteps > 1)
		{
			if (actsubtimestep == 0)
			{
				// if this is the first of several sub-timesteps, the entire amount of PREC - ETR + GWN - IRRIG - ABSTR must be subtractzed again from the PONDS
				// only the fraction for the furst sub time step will be added again
				deltaPond = QGWm3 - ETRm3 - IRRIGm3 - ABSTRm3 + PRECm3;
				Vtemp[0] = Vtemp[0] - deltaPond*(numsubtimesteps-1);
			}
			else
			{
				// in every remaining sub time step, the same deltaPond amount will be added
				Vtemp[0] = Vtemp[0] + deltaPond;
			}
		}
		double cellcnt3 = this->SetPonds(Vtemp[0]); // this call also sets actWaterLevel when lake_active. For -mode:routing, actWaterLevel must be set elsewhere!
		double cellcnt4 = 0;
		for (int i=startrow;i<=endrow;i++)
		{
			for (int j=startcol;j<=endcol;j++)
			{
				if ((LAKE[i][j]==SP_code)&&(POND[i][j]>0.0)) 
				{
					AREAm2 += (cellsize*cellsize);
					cellcnt4 += 1.0;
				}
			}
		}
		for (int z=0;z<=tracer_count;z++)
		{
			Vold[z]=Vakt[z];  // saving the old Vakt-values for later calculation of balances
			Vakt[z]=Vtemp[z]; // Vtemp contains the reservoir-volume and tracer concentraions after redistributing the POND-cells
			// Vtemp will stay in memory for later balancing as well as Vold
		}
	}


	// now that Vakt was saved in Vold, Vakt may be changed
	Vakt[0] -= IRRIGm3;
	Vakt[0] -= ABSTRm3;
	if (Vakt[0] < 0)
	{
		Vakt[0]=0;
	} 

	if (lake_active)
	{
		// calculation of storage change by irrigation water extraction, evaporation, precipitation and groundwater fluxes and abstractions
		DVm3 = QGWm3-ETRm3-IRRIGm3-ABSTRm3+PRECm3;
		ERRm3 = Vakt[0]-Vold[0]-DVm3; 
	}
	else
	{
		DVm3 -= ABSTRm3;
	}

	// aktuelle Abgabe aus dem Speicher: ist Speicher voll? dann maximale Abgabe laut Regel (HW-Entlastungsleistung)
	// Wenn Zufluss < maximale Abgabe: maximale Abgabe, sonst Zufluss (wird ja nicht mehr gespeichert)
	year=p_year;
	month=p_month;
	day=p_day;
	hour=p_hour;
	minute = p_minute;
	
	if (year>0) // no tracer mixing during initializing routing 
	{ 
		// radioactive decay of radioactive tracers
		for (int z=1;z<=tracer_count;z++)
		{
			if (tracer_info[z].radioactive==1) 
			{
				Vakt[z]*=tracer_info[z].decfact;
			}
		}
		// mixing tracer concentrations
		v1=Vakt[0];
		v2=InV*60.0*flowm3s[0];
		for (int z=1;z<=tracer_count;z++)
		{
			Vakt[z]=mixing(Vakt[z],Qges[z][zarray[TG]],v1,v2);
		}
		// aktueller Speicherinhalt: addieren des Zuflusses
		//		Vakt[0] += v2; -> do not add here, because Abakt will be calculated with a sub time step and still relies upon the old volume
	}
	
	// Abakt=0; we use the old value now!
	Abakt = this->GetActualOutflow(fluss,Abstraction,ALs,QINm3s,InV,hour_offset,year,month,day,hour,minute,Abakt,dTmin);
	// GetActualOutflow also calculates abstractions from reservoirs. The abstracted water is not returned in Abakt but must be subtracted now (it is stored in Abstraction[i].Abakt)
	if (year>0) // do not balance the content for initialisation routing
	{
		for (int i=1;i<=ALs;i++)
		{
			if (Abstraction[i].AL_art==8) // only abstraction with type intern_with_rule_from_reservoir should now be processed
			{
				abstraction_summ3 += Abstraction[i].Abakt[0]*InV*60.0; // total abstractions from reservoir in m^3/s
				Vakt[0] -= Abstraction[i].Abakt[0]*InV*60.0;
			}
		}
	}

	// bei falsch aufgestellten Abgaberegeln soll nun verhindert werden, dass mehr aus dem Speicher
	// abgegeben werden kann, als drin ist, indem die ABgabemenge auf den Speicherinhalt begrenzt wird
	// now adding inflow to original volume, since this is now possible and required (GetActualOutflow only changed the volume internally and thus temporarily)
	Vakt[0] += v2;
	if ((Abakt*InV*60.0)>Vakt[0])
	{
		Abakt=Vakt[0]/60.0/InV;
	}

	if (year>0) // no volume change during initialzing routing 
	{ 
		Vakt[0]-=(InV*60.0*Abakt);
		this->SetPonds(Vakt[0]);
	}

	if (lake_active)
	{
		// calculation of storage change by irrigation water extraction, evaporation, precipitation and groundwater fluxes and abstractions
		DVm3 = QINm3 + QGWm3 + PRECm3 - ETRm3 - IRRIGm3 - abstraction_summ3 - Abakt*(InV*60);
		ERRm3 = Vakt[0]-Vold[0]-DVm3; 
	}
	else
	{
		DVm3 = QINm3 - abstraction_summ3 - Abakt * (InV * 60);
		if ((AREAm2 > 0) && (year>0))
		{
			actWaterLevel += DVm3 / AREAm2;
		}
	}
	ABSTRm3  = abstraction_summ3;


	
	// QINm3
	QINmm    = AREAm2>0?QINm3/AREAm2*1000.0:0;

	QOUTm3   = InV*60.0*Abakt;
	QOUTm3s  = Abakt;
	QOUTmm   = AREAm2>0?QOUTm3/AREAm2*1000.0:0;
	// IRRIGm3
	IRRIGmm  = AREAm2>0?IRRIGm3/AREAm2*1000.0:0;
	IRRIGm3s = IRRIGm3/(InV*60.0);
	// ABSTRm3
	ABSTRmm  = AREAm2>0?ABSTRm3/AREAm2*1000.0:0;
	ABSTRm3s = ABSTRm3/(InV*60.0);
	// QGWm3
	QGWmm    = AREAm2>0?QGWm3/AREAm2*1000.0:0;
	QGWm3s   = QGWm3/(InV*60.0);
	// ETRm3
	ETRmm    = AREAm2>0?ETRm3/AREAm2*1000.0:0;
	ETRm3s   = ETRm3/(InV*60.0);
	// PRECm3
	PRECmm   = AREAm2>0?PRECm3/AREAm2*1000.0:0;
	PRECm3s  = PRECm3/(InV*60.0);
	// DVm3
	DVmm     = AREAm2>0?DVm3/AREAm2*1000.0:0;
	DVm3s    = DVm3/(InV*60.0);
	// ERRm3
	ERRmm    = AREAm2>0?ERRm3/AREAm2*1000.0:0;
	ERRm3s   = ERRm3/(InV*60.0);
}

double SP_class::mixing(double c1,double c2, double v1,double v2)
/* calc. mixing concentrations for volume1 (c1) with volume 2 (c2); returns mixing concentration */
{
	/* absolute values are necessary in mixing of storages for overtaken routing water, because this storage can get negative values */
	v1=fabs(v1);
	v2=fabs(v2);
	if ((v1+v2)!=(double)0.0)
	return (v1*c1+v2*c2)/(v1+v2); /* return the mixing concentration */
	else
	return (double)0.0;
}

void SP_class::CreateStatFile()
{
	if ((SP_file=fopen(SP_filename,"w"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<SP_filename<<" for reservoir "<<SP_code<<"[end-entry]";
		throw(37);
	}
	fprintf(SP_file, "YY\tMM\tDD\tHH");
	if (useMinutes) fprintf(SP_file, "\tmm");
	fprintf(SP_file,"\tW[m]\tA[m2]\tV[m3]\tQin[m3]\tQin[mm]\tQin[m3s]\tQout[m3]\tQout[mm]\tQout[m3s]"
	"\tIrrig[m3]\tIrrig[mm]\tIrrig[m3s]\tQGW[m3]\tQGW[mm]\tQGW[m3s]\tETR[m3]\tETR[mm]\tETR[m3s]"
	"\tPREC[m3]\tPREC[mm]\tPREC[m3s]\tdV[m3]\tdV[mm]\tdV[m3s]\tAbstr[m3]\tAbstr[mm]\tAbstr[m3s]\tErr[m3]\tErr[mm]\tErr[m3s]");

	for (int z=1;z<=tracer_count;z++) fprintf(SP_file,"\t%s",tracer_info[z].name);
	
	fprintf(SP_file,"\nreservoir_%d\n",SP_code);
	
	fprintf(SP_file, "YY\tMM\tDD\tHH");
	if (useMinutes) fprintf(SP_file, "\tmm");
	fprintf(SP_file, "\tW[m]\tA[m2]\tV[m3]\tQin[m3]\tQin[mm]\tQin[m3s]\tQout[m3]\tQout[mm]\tQout[m3s]"
	"\tIrrig[m3]\tIrrig[mm]\tIrrig[m3s]\tQGW[m3]\tQGW[mm]\tQGW[m3s]\tETR[m3]\tETR[mm]\tETR[m3s]"
	"\tPREC[m3]\tPREC[mm]\tPREC[m3s]\tdV[m3]\tdV[mm]\tdV[m3s]\tAbstr[m3]\tAbstr[mm]\tAbstr[m3s]\tErr[m3]\tErr[mm]\tErr[m3s]");
	
	for (int z=1;z<=tracer_count;z++) fprintf(SP_file,"\t%s",tracer_info[z].name);
	
	fprintf(SP_file,"\n");
	fclose(SP_file);
}


void SP_class::WriteSpeicherinhalt(long p_jr,long p_mt,long p_tg,long p_se,long p_me)
{
	int z,k;
	if (SP_bufcnt>=SP_bufmax)
	{
		SP_flushbuffer(1);
		SP_bufcnt=SP_valcnt=0;
	}
	SP_bufcnt++;
	k=SP_bufcnt;
	SP_buffer[k][0]=(double)p_jr;
	SP_buffer[k][1]=(double)p_mt;
	SP_buffer[k][2]=(double)p_tg;
	SP_buffer[k][3]=(double)p_se;
	SP_buffer[k][4]=(double)p_me;
	SP_buffer[k][5]=GetWaterLevelByVolume(Vakt[0]);
	SP_buffer[k][6]=AREAm2;
	SP_buffer[k][7]=Vakt[0];
	SP_buffer[k][8]=QINm3;
	SP_buffer[k][9]=QINmm;
	SP_buffer[k][10]=QINm3s;
	SP_buffer[k][11]=QOUTm3;
	SP_buffer[k][12]=QOUTmm;
	SP_buffer[k][13]=QOUTm3s;
	SP_buffer[k][14]=IRRIGm3;
	SP_buffer[k][15]=IRRIGmm;
	SP_buffer[k][16]=IRRIGm3s;
	SP_buffer[k][17]=QGWm3;
	SP_buffer[k][18]=QGWmm;
	SP_buffer[k][19]=QGWm3s;
	SP_buffer[k][20]=ETRm3;
	SP_buffer[k][21]=ETRmm;
	SP_buffer[k][22]=ETRm3s;
	SP_buffer[k][23]=PRECm3;
	SP_buffer[k][24]=PRECmm;
	SP_buffer[k][25]=PRECm3s;
	SP_buffer[k][26]=DVm3;
	SP_buffer[k][27]=DVmm;
	SP_buffer[k][28]=DVm3s;
	SP_buffer[k][29]=ABSTRm3;
	SP_buffer[k][30]=ABSTRmm;
	SP_buffer[k][31]=ABSTRm3s;
	SP_buffer[k][32]=ERRm3;
	SP_buffer[k][33]=ERRmm;
	SP_buffer[k][34]=ERRm3s;
	SP_valcnt=30;

	for (z=1;z<=tracer_count;z++)
	{
		SP_valcnt++;
		SP_buffer[k][z+33]=Vakt[z];
	}
}

void SP_class::ExtendFileName(char *directory,char *dateiname)
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

void SP_class::SetSpeicherParams(FILE *steuerdatei,int p_SPnr,std::stringstream &mystream)
{
	int z;

	SP_filename=new char[_stringsize_+1];
	SP_code=p_SPnr;
	// reading the output file name for reservoir storage
	if (mystream.str().find("file")!=std::string::npos) {
		mystream.seekg(mystream.str().find("file",mystream.tellg())+strlen("file"));
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> SP_filename;
	} else  {
		(*WSstream) << "FATAL: for reservoir "<<p_SPnr<<": no \"file =\" found!"<<"[end-entry]";
		throw(35);
	}
	
	ExtendFileName(DefaultOutputDirectory,SP_filename);

	double *V0;
	V0 = new double[10];
	// reading initialization storage 
	if (mystream.str().find("V0")!=std::string::npos) {
		mystream.seekg(mystream.str().find("V0",mystream.tellg())+strlen("V0"));
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> V0[0];
	} else {
		(*WSstream) << "FATAL: for reservoir "<<p_SPnr<<": no \"V0 =\" found!"<<"[end-entry]";
		throw(35);
	}

	AreaInit = 0.0;
	if (mystream.str().find("A0") != std::string::npos) {
		mystream.seekg(mystream.str().find("A0", mystream.tellg()) + strlen("A0"));
		mystream.seekg(mystream.str().find("=", mystream.tellg()) + 1);
		mystream >> AreaInit;
		this->AREAm2 = this->AreaInit;
	}
	WaterlevelInit = 0.0; // when in -mode:routing, water levels are unknown. In case of a hydraulic connection between two lakes we need to track the water level change therefore explizitely, so we need the initialiation values for water level and areas
	if (mystream.str().find("W0") != std::string::npos) {
		mystream.seekg(mystream.str().find("W0", mystream.tellg()) + strlen("W0"));
		mystream.seekg(mystream.str().find("=", mystream.tellg()) + 1);
		mystream >> WaterlevelInit;
		this->actWaterLevel = this->WaterlevelInit;
	}

	// raeding initialization concentrations for ALL 9 tracers
	if (mystream.str().find("C0")!=std::string::npos) {
		mystream.seekg(mystream.str().find("C0",mystream.tellg())+strlen("C0"));
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> V0[1] >> V0[2] >> V0[3] >> V0[4] >> V0[5] >> V0[6] >> V0[7] >> V0[8] >> V0[9]; 
	} else {
		(*WSstream) << "FATAL: for reservoir "<<p_SPnr<<": no \"C0 =\" found!"<<"[end-entry]";
		throw(35);
	}
	for (z=0;z<=tracer_count;z++) Vakt[z]=V0[z];
	delete[] V0; V0=NULL;

	// reading minimum timestep for dynamic time step splitting (usefull for smaller reservoirs which may run empty within a single interval ptherwise) 
	if (mystream.str().find("dTmin")!=std::string::npos) {
		mystream.seekg(mystream.str().find("dTmin",mystream.tellg())+strlen("dTmin"));
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> dTmin; 
	}
	else
	{
		dTmin=60.0; // minimum time step in secounds (if not specified in the control file
	}

	// reading optional parameter UseRunningAvg for smoother outflow calculation, usefull for larger lakes
	if (mystream.str().find("UseRunningAvg") != std::string::npos) {
		mystream.seekg(mystream.str().find("UseRunningAvg", mystream.tellg()) + strlen("UseRunningAvg"));
		mystream.seekg(mystream.str().find("=", mystream.tellg()) + 1);
		mystream >> UseRunningAvg;
	}
	else
	{
		UseRunningAvg = 0; // minimum time step in secounds (if not specified in the control file
	}


	// allocate memory for output buffers
	SP_bufmax = (unsigned long)(buffersize/((tracer_count+35)*sizeof(double)));
	if (SP_bufmax<1) SP_bufmax=1;
	SP_valcnt   = SP_bufcnt   = 0;
	SP_buffer   = new double*[SP_bufmax+1];
	for (int i=0;i<=SP_bufmax;i++)
	{
		SP_buffer[i]=new double[tracer_count+36]; /* tracer_count + yy + mm + dd + hh + discharge = +5 */
	}

}


/* ~~~~~~~~~~~~~~~~~~~~ Ende Speicher/Reservoirs/Lakes ~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~ Ableitungen/Abstractions ~~~~~~~~~~~~~~~~~ */

AL_class::AL_class():
tabmax(0),AL_buffer(0),AL_bufmax(0),AL_bufcnt(0),AL_valcnt(0),
AL_filename(0),AL_file(0),
AL_code(0),TG(0),tracer_count(0),tracer_info(0),ARTabIndex(0),ARTabValue(0),Abakt(0),ARChangeRate(0),targetCapacity(0),
AL_anteil(0),AL_schwelle(0),AL_max(0),AL_art(0),TgAe(0),TgAe_rel(0),InV(0),InFileName(0),InFile(0),
FB(0),year(0),month(0),day(0),hour(0),startyear(0),startmonth(0),startday(0),starthour(0),
ZL_target(-1), ZL_TgAe(-1), ZL_TgAe_rel(1.0), rulecnt(0), JD(0), WeekDays(0), start_hour(0), stop_hour(0), start_hour2(0), stop_hour2(0), bReverseLookUp(false), useMinutes(false)
{
	mon[1]=31;
	mon[2]=28;
	mon[3]=31;
	mon[4]=30;
	mon[5]=31;
	mon[6]=30;
	mon[7]=31;
	mon[8]=31;
	mon[9]=30;
	mon[10]=31;
	mon[11]=30;
	mon[12]=31;
}

AL_class::~AL_class()
{
	if (AL_code) // delete initialized objects only 
	{
		if (AL_buffer)
		{
			for (int k=0;k<=AL_bufmax;k++)
			{
				if (AL_buffer[k])
				{
					free(AL_buffer[k]);
					AL_buffer[k] = NULL;
				}
			}
		}
		if (ARTabIndex) delete[] ARTabIndex;
		if (Abakt) delete[] Abakt;
		if (AL_buffer) free(AL_buffer);
		if (AL_filename) delete[] AL_filename;
		if (InFileName) delete[] InFileName;
		if (FB) 
		{
			FB->ClearMemory();
			delete FB;
		}
		if (rulecnt>0)
		{
			for (int i=0;i<=rulecnt;i++)
			{
				if (ARTabValue[i]) delete[] ARTabValue[i];
				if (ARChangeRate[i])  delete[] ARChangeRate[i];
			}
		}
		if (ARTabValue) delete[] ARTabValue;
		if (ARChangeRate)  delete[] ARChangeRate;
		if (targetCapacity) delete[] targetCapacity;
		if (JD) delete[] JD;
		if (WeekDays) delete[] WeekDays;
		if (start_hour) delete[] start_hour;
		if (stop_hour) delete[] stop_hour;
		if (start_hour2) delete[] start_hour2;
		if (stop_hour2) delete[] stop_hour2;
	}
}

void AL_class::ExtendFileName(char *directory,char *dateiname)
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

void AL_class::AL_flushbuffer(int code)
{
	unsigned long i,j;
	if ((AL_art==4)||(AL_art==5)||(AL_art==7))
	{
		if ((code>0)&&(AL_bufcnt>0))
		{
			signal(SIGINT,SIG_IGN);
#if defined _BCC_
			signal(SIGBREAK,SIG_IGN);
#endif
			(*WSstream) << "flushing to disk: "<<AL_filename<<" ("<<buffersize<<" buffersize)"<<"[end-entry]";
			if ((AL_file=fopen(AL_filename,"a+"))==NULL)
			{
				(*WSstream) << "FATAL: error opening reservoir file "<<AL_filename<<"[end-entry]";
				throw(65);
			}
			for (i=1;i<=AL_bufcnt;i++)
			{
				fprintf(AL_file,"%2.0lf\t%2.0lf\t%2.0lf\t%2.0lf",AL_buffer[i][0],AL_buffer[i][1],AL_buffer[i][2],AL_buffer[i][3]);
				if (useMinutes)
				{
					fprintf(AL_file, "\t%2.0lf", AL_buffer[i][4]);
				}
				for (j = 5; j <= 4 + AL_valcnt; j++) // even is useMinutes is false, values are stored from element j=5 upwards (element 4 is then only not written
				{
					fprintf(AL_file, "\t%lg", AL_buffer[i][j]);
				}
				fprintf(AL_file,"\n");
			}
			fclose(AL_file);
			AL_bufcnt=AL_valcnt=0;
			SetSignalAbfangfunktion();
		}
	}
}
void AL_class::findeabschnitt(FILE *steuerdatei,char *Abschnittsname)
{
	char *tempstr,*satz,*ttempstr,*tsatz;

	tempstr=(char*)calloc(_stringsize_,sizeof(char)); ttempstr=tempstr;
	satz=(char*)calloc(_stringsize_,sizeof(char)); tsatz=satz;

	fseek(steuerdatei,0,0);
	strcpy(satz,"");
	while ((strcmp(satz,Abschnittsname)!=0)&&(!feof(steuerdatei))&&tempstr)
	{
		fscanf(steuerdatei,"%s",satz);
		fgets(tempstr,_stringsize_,steuerdatei);
	}
	if (feof(steuerdatei))
	{
		(*WSstream) << "FATAL: section "<<Abschnittsname<<" not found in control file"<<"[end-entry]";
		throw(34);
	}
	tempstr=ttempstr; if (tempstr) {free(tempstr); tempstr=NULL;}
	satz=tsatz;       if (satz)    {free(satz);    satz=NULL;	}
}

void AL_class::initialize(wasim::mystream *p_WSstream,tracer_info_type *p_tracer_info,int p_tracer_count,int p_TG,double p_TgAe,double p_TgAe_rel,double p_InV,
char *p_DefaultInputDirectory,char *p_DefaultOutputDirectory,double hour_offset)
{
	WSstream = p_WSstream;
	tracer_info = p_tracer_info;
	tracer_count = p_tracer_count;
	TG = p_TG;
	TgAe = p_TgAe;
	TgAe_rel = p_TgAe_rel;
	InV = p_InV;
	Abakt  = new double[tracer_count+1];
	DefaultOutputDirectory=p_DefaultOutputDirectory;
	DefaultInputDirectory=p_DefaultInputDirectory;
	this->hour_offset = hour_offset;
	if (InV < 60) useMinutes = true;
}


double AL_class::CalcWaterBalance(double **Qges,double &abstraction_sum,double &flowmax,int p_year,int p_month,int p_day,int p_hour,int p_minute)
{
	year=p_year;
	month=p_month;
	day=p_day;
	hour=p_hour;

	double flowmmh=0;
	if ((AL_art==3)||(AL_art==4)) {
		flowmmh = (flowmax-abstraction_sum) >= AL_schwelle ? (flowmax - AL_schwelle)*AL_anteil : 0.0;
		if (flowmax-abstraction_sum-flowmmh<AL_schwelle) 
		{
			flowmmh=flowmax-abstraction_sum-AL_schwelle; /* limiting each abstraction to the value, that the threshold exceeds after subtracting all former abstractions */
		}
		if (flowmmh<0.0) 
		{
			flowmmh=0.0;
		}
		if (flowmmh>AL_max) 
		{
			flowmmh=AL_max;
		}
	} else if (AL_art==5) {
		flowmmh = this->m3s_to_mmh(FB->GetActValue(column-1),TgAe/TgAe_rel,InV);
		if ((flowmax-abstraction_sum) < flowmmh)
		{
			// if the required abstraction cannot be satisfied because of to less 
			// weter in the river or reservoir, even the abstraction gieven from 
			// file will be reduced
			flowmmh=flowmax-abstraction_sum; 
		}
		if (year>0) // do not wind forward in initialization mode
		{
			FB->NextRecord(); 
		}
	} else if ((AL_art==6)||(AL_art==7)) {
		// evaluation the time variable abstraction table here
		flowmmh = GetActualAbstraction(Qges,mmh_to_m3s(flowmax,TgAe/TgAe_rel,InV),year,month,day,hour,minute);
		if ((flowmax-abstraction_sum) < flowmmh)
		{
			// if the required abstraction cannot be satisfied because of to less 
			// weter in the river or reservoir, even the abstraction gieven from 
			// file will be reduced
			flowmmh=flowmax-abstraction_sum; 
		}
	} else if (AL_art==8) {
		// evaluation the time variable abstraction table here like for 6 and 7, but from rerservoir volume instead from discharge
		// to be changed...
		flowmmh = GetActualAbstraction(Qges,mmh_to_m3s(flowmax,TgAe/TgAe_rel,InV),year,month,day,hour,minute);
		if ((flowmax-abstraction_sum) < flowmmh)
		{
			// if the required abstraction cannot be satisfied because of to less 
			// weter in the river or reservoir, even the abstraction gieven from 
			// file will be reduced
			flowmmh=flowmax-abstraction_sum; 
		}
	}
	// ToDo:
	// here, we need to read the actual value like for AL_art == 5 when AL_art==9
	abstraction_sum+=flowmmh;
	return flowmmh;
}

void AL_class::WriteAbstraction(long p_jr,long p_mt,long p_tg,long p_se,long p_me,double *flowm3s)
{
	if ((AL_art==4)||(AL_art==5)||(AL_art==7))
	{
		if (AL_bufcnt>=AL_bufmax)
		{
			AL_flushbuffer(1);
			AL_bufcnt=AL_valcnt=0;
		}
		AL_bufcnt++;
		AL_valcnt=0;
		int k=AL_bufcnt;
		AL_buffer[k][0]=(double)p_jr;
		AL_buffer[k][1]=(double)p_mt;
		AL_buffer[k][2]=(double)p_tg;
		AL_buffer[k][3]=(double)p_se;
		if (useMinutes)
		{
			AL_buffer[k][4] = (double)p_me;
		}
		else
		{
			AL_buffer[k][4] = 0;
		}
		for (int z = 0; z <= tracer_count; z++)
		{
			AL_valcnt++;
			AL_buffer[k][z + 5] = flowm3s[z]; // even if useMinutes is false, values are stored from element 5 upwards
		}
	}
}


double AL_class::m3s_to_mmh(double p_m3s,double p_TgAe,double p_InV)
{
	/* wandelt den absoluten Abfluss in m^3/s in mm/Zeitintervall um  */
	return p_m3s*(double)0.06*p_InV/p_TgAe;
	/* 0.06 = 3.6 / 60 */
}

double AL_class::mmh_to_m3s(double p_mmh,double p_TgAe,double p_InV)
{
	/* wandelt mm/Zeitintervall in m^3/s um */
	return p_mmh*p_TgAe*(double)16.6666667/p_InV;
}


//void AL_class::SetAbstractionMode(int p_ALnr,int p_AL_art,const char *p_AL_filename,double p_AL_schwelle,
//								  double p_AL_anteil,double p_AL_max,const char *p_AL_einheit,int p_subint)
void AL_class::SetAbstractionMode(FILE *steuerdatei,int p_ALnr,std::stringstream &mystream,int subint,int jahra,int monata,int taga,int stundea,int minutea)
{

	char  AL_modus[100];
	char  AL_einheit[100];
	//	char  p_AL_filename[256];
	//  double AL_schwelle,AL_anteil,AL_max;

	AL_code = p_ALnr;

//	strcpy(satz1,&(satz1[(int)(strstr(satz1,"=")-satz1)+1]));
	mystream.seekg(mystream.str().find("=",mystream.tellg())+1);

	// Ableitungsmodus: intern: filename wird nicht gelesen, extern: filename lesen (jede Ableitung eine eigene Datei)
	mystream >> AL_modus;

	if (strstr(AL_modus,")"))
	{
		AL_modus[index(AL_modus,")")]='\0';
	}
	// override the modus word (which length is not known)
	AL_filename=new char[_stringsize_+1];       
        sprintf(AL_filename,"");
	InFileName=new char[_stringsize_+1];
        sprintf(InFileName,"");
	if (strcmp(AL_modus,_intern_)==0)  {
		AL_art=3;
		strcpy(AL_filename,"");
	} else if (strcmp(AL_modus,_extern_)==0) {
		AL_art=4;
		mystream >> AL_filename;
	} else if (strcmp(AL_modus,_from_file_)==0) {
		AL_art=5;

		mystream.seekg(mystream.str().find("input_file",mystream.tellg())+strlen("input_file"));
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> InFileName;
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> column;
		mystream.seekg(mystream.str().find("=",mystream.tellg())+1);
		mystream >> AL_filename;

	} else if (strcmp(AL_modus,_intern_with_rule_)==0) {
		AL_art=6;
		strcpy(AL_filename,"");
	} else if (strcmp(AL_modus,_extern_with_rule_)==0) {
		AL_art=7;
		mystream >> AL_filename;

	} else if (strcmp(AL_modus, _intern_with_rule_from_reservoir_) == 0) {
		AL_art = 8;
		strcpy(AL_filename, "");
	} else if (strcmp(AL_modus, _intern_from_file_from_reservoir_) == 0) {
		AL_art = 9;
		
		// here we need to read the filename and the replace-list (into an array of integers, which must be made available to callers ba specific get functions in the AL-class
		// so, the AL-class needs to be extended by a std::vector<int> for ALs to be replaced, a getter function should return the n-th element of that array when called 

		
		mystream.seekg(mystream.str().find("externalAbstractionFile", mystream.tellg()) + strlen("externalAbstractionFile"));
		mystream.seekg(mystream.str().find("=", mystream.tellg()) + 1);
		mystream >> InFileName >> column;

		// the control file should contain a replacement list of AL-IDs that will be replaced by this external data like this:
		// replaced_Abstractions = {0, 1, 2}
		// this list will now be evaluated and a std::vector is generated that holds the values
		mystream.seekg(mystream.str().find("replaced_Abstractions", mystream.tellg()) + strlen("replaced_Abstractions"));
		mystream.seekg(mystream.str().find("=", mystream.tellg()) + 1);
		mystream.seekg(mystream.str().find("{", mystream.tellg()) + 1);
		int actAL_ID;
		do
		{
			mystream >> actAL_ID;
			if (mystream.good())
			{
				// add value to the refAL-vector
				vrefAL_globalIDs.push_back(actAL_ID);
			}
		}
		while (mystream.good());
		if (mystream.bad() || mystream.fail())
		{
			mystream.clear(mystream.goodbit);
		}

	} else {
		(*WSstream) << "FATAL: the modus \""<<AL_modus
		<<"\" is not supported for abstraction "<<p_ALnr<<"[end-entry]";
		throw(44);
	}
	

	if ((strlen(AL_filename)>0) && strstr(AL_filename,")"))
	{
		AL_filename[index(AL_filename,")")]='\0';
	}
	if (strlen(InFileName)>0)
	{
                if (strstr(InFileName,")"))
                {
		        InFileName[index(InFileName,")")]='\0';
	        }
	}

	// the old style abstraction only allowed threshold, fraction and max capacity parameters, the new style allows a table only or a filename
	if ((AL_art==4)||(AL_art==3))
	{
		// 4: Werte in File schreiben 
		// Art der Ableitung als sub-modus: constant oder Schwelle mit Anteil
		// reading ebstraction threshold, percentage of abstraction relative to the flow exceeding the threshold, 
		//   the capacity of the abstaction channel and the units in which the numbers are given 
		mystream >> AL_schwelle >> AL_anteil >> AL_max >> AL_einheit;
		if (strstr(AL_einheit,")")) // removing all trailing chars after the closing bracket
		AL_einheit[strstr(AL_einheit,")")-AL_einheit]='\0';
		AL_schwelle = ((strstr(AL_einheit,"m^3/s"))||(strstr(AL_einheit,"m3/s"))||(strstr(AL_einheit,"m3s"))) ?
		m3s_to_mmh(AL_schwelle,TgAe/TgAe_rel,InV) : AL_schwelle;
		AL_max      = ((strstr(AL_einheit,"m^3/s"))||(strstr(AL_einheit,"m3/s"))||(strstr(AL_einheit,"m3s"))) ?
		m3s_to_mmh(AL_max,TgAe/TgAe_rel,InV) : AL_max;
	}
	else
	{
		AL_schwelle=0;
		AL_anteil=0;
		AL_max=0;
	}

	// if abstraction will be read from a file, this file must be opened as a FileBuffer
	// and the starting date must be looked for
	// also for AL_art==9 lookup the first date and time now
	if ((AL_art==5)||(AL_art==9))
	{
		FB = new FileBuffer(InFileName,0,1000000);

		if (FB->GetStatus()==FB_STATE_FILENOTFOUND)
		{
			(*WSstream) << "FATAL: error opening "<<FB->GetFilename()<<" for abstraction "<<AL_code<<"[end-entry]";
			throw(37);
		}

		// wind froward to the actual starting date
		startyear=startmonth=startday=starthour=startminute=0;
		FB->SetBufPosRel(-1); // auf eine Position vor Beginn setzen
		while (!((startyear==jahra)&&(startmonth==monata)&&(startday==taga) && (starthour == stundea) && (startminute == minutea)))
		{
			FB->NextRecord();
			FB->GetActdate(startyear,startmonth,startday,starthour,startminute);
			if (startyear<1000) startyear+=1900;
			if (FB->GetStatus()==FB_STATE_LASTRECORDREAD)
			{
				(*WSstream) << "FATAL: starting date "<<taga<<"."<<monata<<"."<<jahra<<", "<<stundea<<". h in "
				<< FB->GetFilename()<<" not found."  << "[end-entry]";
				throw(38);
			}
		}
	}

	if ((AL_art==4)||(AL_art==5)||(AL_art==7)) // modus 4,5 and 7: Werte in File schreiben 
	{
		//		AL_filename=new char[_stringsize_+1];
		//		strcpy(AL_filename,p_AL_filename);
		ExtendFileName(DefaultOutputDirectory,AL_filename);
		if ((AL_file=fopen(AL_filename,"w"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<AL_filename<<" for abstraction "<<AL_code<<"[end-entry]";
			throw(37);
		}
		// first header row
		if ((InV < 60) && (InV > 0))
		{
			fprintf(AL_file, "YY\tMM\tDD\tHH\tmm\tm^3/s");
		}
		else
		{
			fprintf(AL_file, "YY\tMM\tDD\tHH\tm^3/s");
		}

		for (int z=1;z<=tracer_count;z++) fprintf(AL_file,"\t%s",tracer_info[z].name);
		fprintf(AL_file,"\n");

		// second header row
		fprintf(AL_file,"abstraction_%d\n",AL_code);
		
		// third header row (= 1. row)
		if ((InV < 60) && (InV > 0))
		{
			fprintf(AL_file, "YY\tMM\tDD\tHH\tmm\tm^3/s");
		}
		else
		{
			fprintf(AL_file, "YY\tMM\tDD\tHH\tm^3/s");
		}
		for (int z=1;z<=tracer_count;z++) fprintf(AL_file,"\t%s",tracer_info[z].name);
		fprintf(AL_file,"\n");
		
		fclose(AL_file);
		// allocate memory for output buffers
		AL_bufmax = (unsigned)(buffersize/(5*sizeof(long)+(tracer_count+1)*sizeof(double))); // 5 date/time columns plus one data column plus one column for each active tracer
		if (AL_bufmax<1) AL_bufmax=1;
		AL_valcnt = AL_bufcnt = 0;

		AL_buffer = (double**)calloc(AL_bufmax+1,sizeof(double*));
		for (int i=0;i<=AL_bufmax;i++)
		AL_buffer[i]=(double*)calloc(tracer_count+6,sizeof(double)); // 5 date/time columns plus one data column plus one column for each active tracer

	}
}

void AL_class::ReadAbstractionRule(FILE *steuerdatei)
{
	if ((AL_art==6)||(AL_art==7)||(AL_art==8))
	{
		char *tstr;
		tstr=(char*)calloc(_stringsize_,sizeof(char));
		sprintf(tstr,"[abstraction_rule_abstraction_%d]",AL_code);
		findeabschnitt(steuerdatei,tstr);
		if (!feof(steuerdatei)) tstr=fgets(tstr,_stringsize_,steuerdatei);
		// Anzahl der Eintraege der Abgaberegel 
		sscanf(tstr,"%d",&(tabmax));
		// analysing the number of julian days
		std::string s(tstr);
		rulecnt=0;
		int pos=s.find_first_not_of(" \t\n\r\0");
		while ((s.length()>0)&&(pos!=std::string::npos))
		{
			// remove leading whitespaces
			if (s.find_first_not_of(" \t")>0) 
			{
				s=s.substr(s.find_first_not_of(" \t"));
			}
			s=s.substr(s.find_first_not_of("0123456789.eE+-"));
			rulecnt++;
			pos=s.find_first_not_of(" \t\n\r\0");
		}
		// since the first entry was not a Julian Day, the rulecnt is 1 to high (except for old style rules with no JD at all)
		if (rulecnt>1) rulecnt--;
		// now the arrays can be allocated
		ARTabIndex = new double[tabmax+1];
		ARTabValue = new double*[rulecnt+1];
		ARChangeRate  = new double*[rulecnt+1];
		for (int i=0;i<=rulecnt;i++)
		{
			ARTabValue[i]= new double[tabmax+1];
			ARChangeRate[i] = new double[tabmax+1];
		}
		targetCapacity = new double[rulecnt+1];
		JD = new int[rulecnt+1];
		WeekDays= new long[rulecnt+1];
		start_hour = new int[rulecnt + 1];
		stop_hour = new int[rulecnt + 1];
		start_hour2 = new int[rulecnt + 1];
		stop_hour2 = new int[rulecnt + 1];

		// read Julian Days from the first row (from std::string s)
		s=std::string(tstr);
		// overide first value (is counter of rows)
		if (rulecnt==1)
		{
			JD[1]=1;
		}
		else
		{
			if (s.find_first_not_of(" \t")>0) 
			{
				s=s.substr(s.find_first_not_of(" \t"));
			}
			s=s.substr(s.find_first_not_of("0123456789.eE+-"));
			for (int i=1;i<=rulecnt;i++)
			{
				std::stringstream ss(s.c_str());
				ss >> JD[i];
				if (s.find_first_not_of(" \t")>0) 
				{
					s=s.substr(s.find_first_not_of(" \t"));
				}
				s=s.substr(s.find_first_not_of("0123456789.eE+-"));
			}
		}
		// now reading the table
		for (int t=1;t<=tabmax;t++)
		{
			if (!feof(steuerdatei)) tstr=fgets(tstr,_stringsize_,steuerdatei);
			s=std::string(tstr);
			for (int i=0;i<=rulecnt;i++)
			{
				std::stringstream ss(s.c_str());
				if (i==0)
				{
					ss >> ARTabIndex[t];
					//					ARTabIndex[t]=m3s_to_mmh(ARTabIndex[t],TgAe/TgAe_rel,InV);
				}
				else
				{
					ss >> ARTabValue[i][t];
					//					Abgabe[i][t]=m3s_to_mmh(Abgabe[i][t],TgAe/TgAe_rel,InV);
				}
				if (s.find_first_not_of(" \t")>0) 
				{
					s=s.substr(s.find_first_not_of(" \t"));
				}
				s=s.substr(s.find_first_not_of("0123456789.eE+-"));
			}
		}

		// initializing additional parameters in case they are not present in the control file
		for (int i = 1; i <= rulecnt; i++)
		{
			targetCapacity[i] = -1.0;
			WeekDays[i] = -1;
			start_hour[i] = -1;
			stop_hour[i] = -1;
			start_hour2[i] = -1;
			stop_hour2[i] = -1;
		}

		// reading target capacities
		if (!feof(steuerdatei))
		{
			tstr = fgets(tstr, _stringsize_, steuerdatei);
			s = std::string(tstr);
		}
		while ((s.find("[") == std::string::npos) && (s != "\n") && (s.length() != 0) && !feof(steuerdatei))
		{
			if (s.find("TargetCap") != std::string::npos)
			{
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> targetCapacity[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("0123456789.eE+-"));
				}
			}
			else if (s.find("WeekDays") != std::string::npos) {
				// Reading WeekDays
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> WeekDays[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567"));
				}
			}
			else if (s.find("start_hour") != std::string::npos) {
				// Reading start_hour
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> start_hour[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("stop_hour") != std::string::npos) {
				// Reading stop_hour
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> stop_hour[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("start2_hour") != std::string::npos) {
				// Reading start_hour2
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> start_hour2[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("stop2_hour") != std::string::npos) {
				// Reading stop_hour2
				s = s.substr(s.find("=") + 1);
				for (int i = 1; i <= rulecnt; i++)
				{
					std::stringstream ss(s.c_str());
					ss >> stop_hour2[i];
					if (s.find_first_not_of(" \t") > 0)
					{
						s = s.substr(s.find_first_not_of(" \t"));
					}
					s = s.substr(s.find_first_not_of("1234567890"));
				}
			}
			else if (s.find("TableIndex") != std::string::npos) {
				s = s.substr(s.find("=") + 1);
				if (s.find("Flow") != std::string::npos)
				{
					bReverseLookUp = true;
				}
				else
				{
					bReverseLookUp = false;
				}
			}

			if (!feof(steuerdatei))
			{
				{
					tstr = fgets(tstr, _stringsize_, steuerdatei);
					s = std::string(tstr);
				}
			}
		}

		// interpolating between the table entries or table indices, depending on the state of breverseLookUp
		for (int i = 1; i <= rulecnt; i++)
		{
			if (bReverseLookUp)
			{
				for (int t = 1; t < tabmax; t++)
				{

					ARChangeRate[i][t] = (ARTabValue[i][t + 1] - ARTabValue[i][t]) != 0 ?
						(ARTabIndex[t + 1] - ARTabIndex[t]) /
						(ARTabValue[i][t + 1] - ARTabValue[i][t]) : 0.0;
				}
				ARChangeRate[i][tabmax] = 0.0;
			}
			else // default
			{
				for (int t = 1; t < tabmax; t++)
				{

					ARChangeRate[i][t] = (ARTabIndex[t + 1] - ARTabIndex[t]) != 0 ?
						(ARTabValue[i][t + 1] - ARTabValue[i][t]) /
						(ARTabIndex[t + 1] - ARTabIndex[t]) : 0.0;
				}
				ARChangeRate[i][tabmax] = 0.0;
			}
		}
		if (tstr) { free(tstr); tstr = NULL; }
	}
}

void  AL_class::SetTargetData(int p_ZL_target,double p_ZL_TgAe,double p_ZL_TgAe_rel,int flussindex,double SourceAE)
{
	ZL_target = p_ZL_target;
	ZL_TgAe = p_ZL_TgAe;
	ZL_TgAe_rel = p_ZL_TgAe_rel;
	ZL_flussindex=flussindex;
	ZL_SourceAE=SourceAE;
}

double AL_class::GetActualAbstraction(double **Qges,double Qact,int p_year,int p_month,int p_day,int p_hour,int p_minute)
{
	year=p_year;
	month=p_month;
	day=p_day;
	hour=p_hour;
	minute = p_minute;
	int actJD=0;
	for (int i=1;i<month;i++) actJD+=mon[i];
	actJD+=day;
	int rule=1;
	while ((actJD>JD[rule])&&(rule<rulecnt)) rule++;

	double WeightNextRule = 1;
	double WeightLastRule = 0;
	double dLastAbakt;
	int nextJD = JD[rule];
	int LastJD;
	if (rule>1)
	{
		LastJD = JD[rule - 1];
		WeightNextRule = double(actJD - LastJD) / double(nextJD - LastJD);
		WeightLastRule = 1.0 - WeightNextRule;
	}
	else
	{
		WeightNextRule = 1.0;
		WeightLastRule = 0.0;
	}

	double flowmmh=0;
	double flowm3s=0;
	// Ermitteln der aktuellen SpeicherARTabValue 
	if (bReverseLookUp == false) // old, default version with a table that has volumes as index values and fluxes as table values
	{
		if (Qact > ARTabIndex[tabmax])
		{
			flowm3s = ARTabValue[rule][tabmax];
		}
		else
		{
			for (int j = tabmax - 1; j >= 1; j--)
			{
				if (Qact > ARTabIndex[j])
				{
					flowm3s = ARTabValue[rule][j] +
						(Qact - ARTabIndex[j])*ARChangeRate[rule][j];
					if (rule > 1)
					{
						dLastAbakt = ARTabValue[rule - 1][j] +
							(Qact - ARTabIndex[j])*ARChangeRate[rule - 1][j];
						// interpolating between abstractions in time
						flowm3s = flowm3s * WeightNextRule + dLastAbakt * WeightLastRule;
					}
					break;
				}
			}
		}
	}
	else
	{
		for (int j = tabmax - 1; j >= 1; j--)
		{
			double ARTabValueInterpol = (rule > 1) ?
				ARTabValue[rule][j] * WeightNextRule + ARTabValue[rule - 1][j] * WeightLastRule :
				ARTabValue[rule][j];

			if (Qact > ARTabValueInterpol)
			{
				double ARChangeRateInterpol = (rule > 1) ?
					ARChangeRate[rule][j] * WeightNextRule + ARChangeRate[rule - 1][j] * WeightLastRule :
					ARChangeRate[rule][j];

				flowm3s = ARTabIndex[j] +
					(Qact - ARTabValueInterpol)*ARChangeRateInterpol;
				break;
			}
		}
	}

	// and take into account the target capacity, if possible
	if ((ZL_target>0) && (ZL_flussindex>0))
	{
//		double Qt=mmh_to_m3s(Qges[0][ZL_target],ZL_TgAe/ZL_TgAe_rel,InV);
//		if (Qt+flowm3s>targetCapacity[rule])
//		{
//			flowm3s=targetCapacity[rule]-Qt;
//			if (flowm3s<0) flowm3s=0.0;
//		}
		double Qin_MainOL = mmh_to_m3s(Qges[0][ZL_flussindex], ZL_SourceAE, InV);
		if (Qin_MainOL + flowm3s>targetCapacity[rule])
		{
			flowm3s = targetCapacity[rule] - Qin_MainOL;
			if (flowm3s<0) flowm3s = 0.0;
		}
	}
	flowmmh = m3s_to_mmh(flowm3s,TgAe/TgAe_rel,InV);
	return flowmmh;
}

// 1. Wochentag direkt mathematisch berechnen (Sakamoto's Algorithmus)
// Liefert: 0 = Sonntag, 1 = Montag, ..., 6 = Samstag
int AL_class::get_weekday(int y, int m, int d) {
	static int t[] = { 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 };
	if (m < 3) y -= 1;
	return (y + y / 4 - y / 100 + y / 400 + t[m - 1] + d) % 7;
}


double AL_class::GetActualAbstractionFromReservoir(double ***fluss,double Vnew,int p_year,int p_month,int p_day,int p_hour,int p_minute)
{
	year=p_year;
	month=p_month;
	day=p_day;
	hour=p_hour;
	minute = p_minute;
	int actJD=0;
	for (int i=1;i<month;i++) actJD+=mon[i];
	actJD+=day;
	int rule=1;
	while ((actJD>JD[rule])&&(rule<rulecnt)) rule++;

	double WeightNextRule = 1;
	double WeightLastRule = 0;
	double dLastAbakt;
	int nextJD = JD[rule];
	int LastJD;
	if (rule>1)
	{
		LastJD = JD[rule - 1];
		WeightNextRule = double(actJD - LastJD) / double(nextJD - LastJD);
		WeightLastRule = 1.0 - WeightNextRule;
	}
	else
	{
		WeightNextRule = 1.0;
		WeightLastRule = 0.0;
	}

	bool bApplyRule = false;

// regard WeekDays, start_hour and stop_hour
// 	   old method which was unbeleavable slow due to string stream conversion
//	time_t acttime;
//	tm *tm1 = new(tm);
//	if (month>0)
//	{
//		tm1->tm_hour = hour - (int)hour_offset; // 0...23
////		tm1->tm_hour = hour; // 0...23
//		tm1->tm_min = useMinutes ? minute : 0; // 0..60 
//		tm1->tm_sec = 0; // 0..60
//		tm1->tm_mon = month-1;  // 0..11
//		tm1->tm_year = year - 1900; // e.g. 2010 = 110
//		tm1->tm_mday = day; 
//		tm1->tm_wday = 0; 
//		tm1->tm_yday = 0; 
//		tm1->tm_isdst = 0; 
//
//#ifdef WIN64
//		acttime = _mkgmtime(tm1);
//#else
//		acttime = timegm(tm1);
//#endif
//		//acttime = mktime(&tm1);
//	
//		std::stringstream s1,s2;
//		if (tm1->tm_wday>0)
//		{
//			s1 << tm1->tm_wday;
//		}
//		else
//		{
//			s1 << 7; // since Sunday is not stored as 7 but as 0 (days since Sunday=0)
//		}
//		s2 << WeekDays[rule];
//		std::string wday;
//		std::string weekdays;
//		s1 >> wday;
//		s2 >> weekdays;
//		if (weekdays.find(wday.c_str())!=std::string::npos) 
//		{
//			if (
//				((start_hour[rule] <= stop_hour[rule]) && ((hour >= start_hour[rule]) && (hour <= stop_hour[rule]))) // conventional, e.g. start=7, stop=15: all hours between 7 and 15 will apply the rule
//				||
//				((start_hour[rule] > stop_hour[rule]) && ((hour >= start_hour[rule]) || (hour <= stop_hour[rule]))) // over bight, e.g. start = 16, stop = 6: all hours between 16 and 24 and all hours between 1 and 6 will apply
//				)
//			{
//				bApplyRule = true;
//			}
//		}
//	}

// regard WeekDays, start_hour and stop_hour
// 	   new method which is much faster (found by asking AI to optimize the above commented section)
	if (month > 0) // valid month, day, hour and minute must be provided, otherwise the weekday cannot be calculated
	{
		// get week day without timoe consuming mkdate function calls
		int wday = get_weekday(year, month, day);

		// ISO-week day (1 = Monday, ..., 7 = Sunday)
		int iso_wday = (wday == 0) ? 7 : wday;

		// much more efficient week day check (instead of using stringstream)
		// precondition here: WeekDays[rule] contains e.g. "12345" for Mo-Fr.
		// If WeekDays is a std::string, we use char-searching instead of ttring-conversion
		char wday_char = '0' + iso_wday;
		const std::string& rule_days = std::to_string(WeekDays[rule]);

		if (rule_days.find(wday_char) != std::string::npos) {
			// 4. Stunden-Regel (Bleibt performant, da reine CPU-Register-Vergleiche)
			int start = start_hour[rule];
			int stop = stop_hour[rule];
			int start2 = start_hour2[rule];
			int stop2 = stop_hour2[rule];

			if (
				(start <= stop && hour >= start && hour <= stop) ||
				(start > stop && (hour >= start || hour <= stop)) ||
				(start2 <= stop2 && hour >= start2 && hour <= stop2) ||
				(start2 > stop2 && (hour >= start2 || hour <= stop2))
				)
			{
				bApplyRule = true;
			}
		}
	}


	double flowmmh=0;
	double flowm3s=0;
	if (bApplyRule==true)
	{
		// Ermitteln der aktuellen Speicherabgabe 
		if (bReverseLookUp == false) // old, default version with a table that has volumes as index values and fluxes as table values
		{
			if (Vnew > ARTabIndex[tabmax])
			{
				flowm3s = ARTabValue[rule][tabmax];
			}
			else
			{
				for (int j = tabmax - 1; j >= 1; j--)
				{
					if (Vnew > ARTabIndex[j])
					{
						flowm3s = ARTabValue[rule][j] +
							(Vnew - ARTabIndex[j])*ARChangeRate[rule][j];
						if (rule > 1)
						{
							dLastAbakt = ARTabValue[rule - 1][j] +
								(Vnew - ARTabIndex[j])*ARChangeRate[rule - 1][j];
							// interpolating between abstractions in time
							flowm3s = flowm3s * WeightNextRule + dLastAbakt * WeightLastRule;
						}
						break;
					}
				}
			}
		}
		else
		{
			for (int j = tabmax - 1; j >= 1; j--)
			{
				double ARTabValueInterpol = (rule > 1) ?
					ARTabValue[rule][j] * WeightNextRule + ARTabValue[rule - 1][j] * WeightLastRule :
					ARTabValue[rule][j];

				if (Vnew > ARTabValueInterpol)
				{
					double ARChangeRateInterpol = (rule > 1) ?
						ARChangeRate[rule][j] * WeightNextRule + ARChangeRate[rule - 1][j] * WeightLastRule :
						ARChangeRate[rule][j];

					flowm3s = ARTabIndex[j] +
						(Vnew - ARTabValueInterpol)*ARChangeRateInterpol;
					break;
				}
			}

		}


		// and take into account the target capacity, if possible; Attention: here, we interpret the target capacity as the capacity of the abstraction channel itself 
		// (so it is not really the target channel, but the abstraction itself, limiting the abstracted amount of water per time step)
		// new in 9.00.13: target capacity should be the flow in the target subbasin, since the abstraction itself is restricted by the rule already.
		// but attention: if the target subbasin has no routing into other basins (e.g. because it is the final target or there are multiple subbasins acting as sink or only as source for abstractions etc.
		// then fluss[0][ZL-target] will be NULL. So this must be checked before trying to balance th abstraction against any target restrictions
		if ((ZL_target>0) && (ZL_flussindex>0) && (fluss[0][ZL_flussindex]!=NULL))
		{
//			double Qt = mmh_to_m3s(fluss[0][ZL_target][1], ZL_TgAe / ZL_TgAe_rel, InV);
//			double Qin_old = mmh_to_m3s(fluss[0][ZL_flussindex][1], ZL_SourceAE, InV);
//			if (Qt - Qin_old + flowm3s>targetCapacity[rule])
//			{
//				flowm3s = targetCapacity[rule] - (Qt - Qin_old);
//				if (flowm3s<0) flowm3s = 0.0;
//			}

			double Qin_MainOL = mmh_to_m3s(fluss[0][ZL_flussindex][1], ZL_SourceAE, InV);
			if (Qin_MainOL + flowm3s>targetCapacity[rule])
			{
				flowm3s=targetCapacity[rule]-Qin_MainOL;
				if (flowm3s<0) flowm3s=0.0;
			}
		}

// old version which is limiting the flux by a constant values (rather the divice3 itself, e.g. a tube or so)
//		if (flowm3s>targetCapacity[rule])
//		{
//			flowm3s=targetCapacity[rule];
//			if (flowm3s<0) flowm3s=0.0;
//		}
	}
//	flowmmh = m3s_to_mmh(flowm3s,TgAe/TgAe_rel,InV);
//	this->Abakt[0] = flowm3s;
//	return flowmmh;

//	delete tm1;

	return flowm3s;
}
// ---------------------------------- end of Abstractions -----------------------------
