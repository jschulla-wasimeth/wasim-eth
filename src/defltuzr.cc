/*   $Revision: 114 $ 
	$Date: 22.02.11 21:03 $  */

/*
(c) 1995, 96, 97, 2006, 2007 by Joerg Schulla

	Version history: see source safe comments
	
*/

// comment added from Linux for test of write access
#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <signal.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stddef.h>
#include <map>
#include <algorithm>
#include "stringp.h"
#include "modconst.h"
#include "modef.h"
#include "wasim.h"

using namespace std;

char* defaultini::compilestring(char *satz,int zeile,char *symbol[],char *eintrag[],int anzahleintraege)
{
	char *satzaus,*tsatzaus;
	char *muster,*tmuster;
	int  i,start,stop;

	satzaus = newstring(_stringsize_+1,1); tsatzaus=satzaus;
	muster  = newstring(_stringsize_+1,1); tmuster=muster;

	/* wenn ein '$' in diesem Satz steht, dann wird das Symbol ausgetauscht
	mit dem zugehoerigen Eintrag */
	i = index(satz,"#");
	if (i>0)
	{
		satz[i-1]='\n';
		satz[i]='\0';
	}
	strcpy(satzaus,satz);
	if (!strstr(satz,"$set")) /* variablendeklarationen selbst sollen in der uebersetzten Datei erhalten bleiben */
	{
		while (strstr(satz,"$"))
		{
			/* sollten mehere Symbole in einem Satz stehen, muss der Satz
		so wiederhergestellt werden, dass das i. Symbol ersetzt ist.
		Normalerweise wird bei der substr-function der Eingabesatz
		voellig deformiert und nur als das Ergebnis zurueckgeliefert */
			sscanf(strstr(satz,"$"),"%s",muster);
			if (strstr(muster,"//")) muster[index(muster,"//")]='\0';
			if (strstr(muster,";")) muster[index(muster,";")]='\0';
			for (i=1;i<=anzahleintraege;i++)
			if (strcmp(muster,symbol[i])==0)
			{
				if (index(satz,symbol[i])==0)
				strcpy(satzaus,eintrag[i]);
				else
				strcat(substr(satzaus,0,index(satz,symbol[i])-1),eintrag[i]);
				start = index(satz,symbol[i])+strlen(symbol[i]);
				stop  = strlen(satz);
				strcat(satzaus,substr(satz,start,stop));
				strcpy(satz,satzaus);
				break;
			}
			else
			{
				if ((i>=anzahleintraege)&&(index(satz,"$")>0))
				{
					(*WSstream) << "FATAL: \n-----------------------------------------------------------"
					<< "\npossible error in control file!\n"
					<< "undefined symbol at line "<<zeile<<":\n"<<satz
					<< "\n\n-----------------------------------------------------------\n"<<"[end-entry]";
					throw(17);
				}
			}
		}
		/* String-Concatenation */
		while (strstr(satz,"//")/*!=NULL*/)
		{
			i=index(satz,"//");
			/* Bedeutung der folgenden Zeile: wenn in satz "//" auftaucht, dann gibt i den
		Index im Satz an, an der der erste "/" steht. Da die Zaehlung mit 0 beginnt, liegen also
		i Zeichen vor diesem ersten "/". die Formatanweisung %.*s schreibt maximal soviele
		Zeichen in den Ausgabesatz, wie durch das naechste Argument der Argumentliste (das * steht
		also fuer i in der Argumentenliste) angegeben werden. das naechste Argument beginnt erst
		ab dem i+2. Index in satz, also nach dem zweiten "/". so werden beide "/" aus dem Satz
		eleminiert */
			sprintf(satzaus,"%.*s%s",i,satz,&(satz[i+2]));
			strcpy(satz,satzaus);
		}
	}
	satzaus=tsatzaus; if (satzaus) free(satzaus); satzaus = NULL;
	muster=tmuster;   if (muster) free(muster);  muster = NULL;
	return satz;
}

void defaultini::free_pointers(char *p1,char *p2,char *p3,char *p4,char *p5,char **p6)
{
	if (p1) free(p1);  p1=NULL;
	if (p2) free(p2);  p2=NULL;
	if (p3) free(p3);  p3=NULL;
	if (p4) free(p4);  p4=NULL;
	if (p5) free(p5);  p5=NULL;
	if (p6) free(p6);  p6=NULL;
}

void defaultini::replacesymbols(char *dateinameein,char *dateinameaus)
/* Diese Routine sucht in der eingabedatei nach Eintraegen, welche die Syntax

set $<symbol> = <eintrag>

aufweisen (Bsp.: set $Anfangstag = 1) und ersetzt alle weiteren Eintraege
von $<symbol> mit <eintrag> (sozusagen eine Entsprechung fuer die
#define - Compiler Anweisung von C.
Die Ausgabedatei wird dann wie gewohnt durch die Einleseprogramme weiterverarbeitet.

Eingabeparameter: dateinameein: Name der Datei, welche die Symbole noch enthaelt, also i.d.R. der
								Name der Steuerdatei, wie sie vom Anwender als Parameter fuer das
								modell angegeben wird
					dateinameaus: Name der Datei, welche vom modell dann als gueltige Steuerdatei
						genutzt wird (Rueckgabeparameter, Dateinameein durch ".$$$" erweitert,
*/
{
	FILE       *steuerdateiarg;
	char       *satz,*tempsatz2;
	char *symbol[maxsymbolanzahl],*eintrag[maxsymbolanzahl]; /* maximal 200 Symbole duerfen definiert werden, dafuer muessen
								naemlich je 130 Byte Speicherplatz zur Verfuegung gestellt werden,
					was maximal 26000 Byte + 1600 Byte fuer die Zeiger verlangt */
	char *set,*tempset,*gleich,*tempgleich;
	int  k,anzahleintraege,i,zeile;
	char *satzaus,*tempsatz;
	char **errpntr,**temperrpntr;
	char *muster,*tempmuster;
	char *zeiger,dirz;
	int  wert,minpos;

	satz    = newstring(2*_stringsize_+1,1);
	satzaus = newstring(2*_stringsize_+1,1);
	set     = newstring(100,1);
	gleich  = newstring(100,1);
	errpntr = (char **)calloc(2,sizeof(char*));
	muster  = newstring(2*_stringsize_+1,1);

	tempsatz2=satz;
	tempsatz=satzaus;
	tempset=set;
	tempgleich=gleich;
	temperrpntr=errpntr;
	tempmuster=muster;

	anzahleintraege = 0;

	strcpy(dateinameaus,dateinameein);
	/* die Extension wird mit '$$$' ersetzt. Hat der Dateiname keinen Punkt, dann wird '.$$$' an den
	Dateinamen angefuegt */

	zeiger=strrchr(dateinameaus,dirzeichen);
	if (zeiger)  minpos=(int)(zeiger-dateinameaus); else minpos=0;
	/* nach dem letzten directory-Zeichen (entweder '\\' fuer Dos oder '/' fuer UNIX wird das
letzte Auftreten eines Punktes gesucht, wenn es sich um DOS handelt, da bei UNIX auch mehrere Punkte
im Daeinamen vorkommen koennen, bei DOS aber nicht */
//	dirz=dirzeichen;
//	if (dirz == '\\')
//	{
//		zeiger=strrchr(&(dateinameaus[minpos+1]),'.');
//		if (zeiger) dateinameaus[(int)(zeiger-dateinameaus)]='\0';
//	}

	/*  if (index(dateinameaus,".")>0)
	dateinameaus = substr(dateinameaus,0,index(dateinameaus,".")-1); */

//	sprintf(dateinameaus, "%s_%d.$$$", dateinameaus, mpi_rank);
	sprintf(dateinameaus, "%s.intern", dateinameein);

	(*WSstream) << "compiling symbols in control file " << dateinameein << " to " << dateinameaus << "[end-entry]";

//	strcat(dateinameaus,".$$$");

	/* Suchen von 'set'- eintraegen in der Eingabedatei */

	if ((steuerdateiarg = fopen(dateinameein,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<dateinameein<<"[end-entry]";
		free_pointers(satz,satzaus,set,muster,gleich,errpntr);
		throw(28);
	}

	zeile=0;
	while (! feof(steuerdateiarg))
	{
		/* Satz fuer Satz lesen */
		if (! feof(steuerdateiarg)) satzaus = fgets(satzaus,_stringsize_,steuerdateiarg);
		zeile++;
		if (satzaus==NULL) satzaus = tempsatz;
		strcpy(satz," ");
		strcat(satz,satzaus);
		if ((satz[1]!='#')&&(index(satz,"$set")>0))
		{
			anzahleintraege += 1;
			eintrag[anzahleintraege] = newstring(300,1);
			symbol[anzahleintraege]  = newstring(300,1);
			sscanf(satz,"%s%s%s%s",set,symbol[anzahleintraege],gleich,eintrag[anzahleintraege]);
			eintrag[anzahleintraege]=compilestring(eintrag[anzahleintraege],zeile,symbol,eintrag,anzahleintraege);
			k=anzahleintraege;
			/* nun folgt eine Ueberpruefung, ob in zwei Variabeln zufaellig dergleiche Wert steht. Wenn ja, dann wird eine Warnung ausgegeben
		Wenn aber die Variable mit INF, Inf, inf, NAN, Nan oder nan anfaengt, dann gibt es auf verschiedenen MAschinen floating-errors, weshalb diese
		von der function strtod als gueltige Zahlen interpretierte Werte hier abgefangen werden muessen, die Warnung dann also auch nicht ausgegeben werden soll */
			if (!(strstr(eintrag[k],"inf")||strstr(eintrag[k],"nan")||strstr(eintrag[k],"INF")||strstr(eintrag[k],"NAN")||strstr(eintrag[k],"Inf")||strstr(eintrag[k],"Nan")))
			{
				wert=(unsigned long)strtod(eintrag[anzahleintraege],errpntr);
				if (**errpntr!=0)
				{
					for (i=1;i<anzahleintraege;i++)
					{
						if (strcmp(eintrag[i],eintrag[anzahleintraege])==0)
						{
							(*WSstream) << "INFO: \n------------------------ WARNING --------------------------"
							<< "\npossible error in control file!\n"
							<< "content of symbol "<<symbol[i]<<" is the same as symbol "
							<< symbol[anzahleintraege]<<", \n (which is: "<<eintrag[i]<<").\n"
							<< "this may lead to collusions!"
							<< "\n-----------------------------------------------------------\n"<<"[end-entry]";
						}
					}
				}
			}
		}
	}
	if (satzaus==NULL) satzaus = tempsatz;

	/* Kopieren von Steuerdateiarg nach steuerdate und dabei ersetzen der Symbole */

	if (fseek(steuerdateiarg,0,0)!=0)
	{
		(*WSstream) << "FATAL: error rewinding "<<dateinameein<<"[end-entry]";
		free_pointers(satz,satzaus,set,muster,gleich,errpntr);
		throw(28);
	}

	if ((steuerdatei = fopen(dateinameaus,"wt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<dateinameaus<<"[end-entry]";
		free_pointers(satz,satzaus,set,muster,gleich,errpntr);
		throw(28);
	}

	zeile=0;
	while (! feof(steuerdateiarg))
	{
		if ( ! feof(steuerdateiarg)) satzaus = fgets(satzaus,_stringsize_,steuerdateiarg);
		zeile++;
		if (satzaus==NULL) satzaus = tempsatz;
		strcpy(satz," ");
		strcat(satz,satzaus);
		strcpy(satzaus,satz);
		/*    if (satz[1]!='#') */ /* Kommentarzeilen mit einem # am Anfang werden sofort ausgelassen */
		satz=compilestring(satz,zeile,symbol,eintrag,anzahleintraege);
		if (satz[1]!='$') /* $set-Anweisungen sollen nicht mehr in der einzulesenden Steuerdatei erscheinen */
		fputs(&(satz[1]),steuerdatei);
	}
	if (satzaus==NULL) satzaus = tempsatz;

	if ((fclose(steuerdateiarg))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<dateinameein<<"[end-entry]";
		free_pointers(satz,satzaus,set,muster,gleich,errpntr);
		throw(11);
	}

	if ((fclose(steuerdatei))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<dateinameaus<<"[end-entry]";
		free_pointers(satz,satzaus,set,muster,gleich,errpntr);
		throw(11);
	}


	for (i=1;i<=anzahleintraege;i++)
	{
		if (symbol[i]) free(symbol[i]); symbol[i]=NULL;
		if (eintrag[i]) free(eintrag[i]); eintrag[i]=NULL;
	}
	/* sicherstellen, dass alle Zeiger auf die urspruenglich reservierte ADresse zeigen */
	free_pointers(satz,satzaus,set,muster,gleich,errpntr);
}

void defaultini::GetBasicControls(char *steuerdateiname,int *progaus,int *warning_level,int *runoff_measure,double *minutes_to_interval_end,double *MinutesToInvEndRouting,int *p_timestep,
int *anzahlgroessen,double *breite,double *laenge,double *bezugsmeridian,double *timekorr,
char ***meteonamen,int *numsubb,int **listsubb)
/* liest aus der Steuerdatei die Anzahl der Meteogroessen sowie deren
Namen und die geographischen Koordinaten ein */
{
	int i;
	char *satz,*tempsatz;

	satz    = (char *)calloc(_stringsize_+1,sizeof(char));
	tempsatz=satz;
	if ((steuerdatei = fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: Error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}

	*warning_level=0;
	findeabschnitt(_Programmausgabe_);
	if (fscanf(steuerdatei,"%d",progaus)!=1) steuererror("output control","");
	fgets(satz,_stringsize_,steuerdatei);
	fscanf(steuerdatei,"%d",warning_level);
	if ((*warning_level <0)||(*warning_level>1)) *warning_level=1;
	fgets(satz,_stringsize_,steuerdatei);
	fscanf(steuerdatei,"%d",runoff_measure);
	if ((*runoff_measure<0)||(*runoff_measure>1)) *runoff_measure=0;
	fgets(satz,_stringsize_,steuerdatei);
	fscanf(steuerdatei,"%lf",minutes_to_interval_end);
	if (*minutes_to_interval_end<0) *minutes_to_interval_end=0;
	fgets(satz,_stringsize_,steuerdatei);
	*MinutesToInvEndRouting = *minutes_to_interval_end;
	bool DoContinue=true;
	do
	{
		fgets(satz,_stringsize_,steuerdatei);
		std::string s(satz);
		s=s.substr(s.find_first_not_of(" \t"));
		sprintf(satz,"%s",s.c_str());
		if ((!feof(steuerdatei))&&(index(satz,"\n")>0)) // skip all empty lines
		{
			DoContinue=true;
			std::string parameter(satz);
			if (parameter.find("WriteAsciiGrids")!=std::string::npos) {
				std::stringstream ss;
				parameter=parameter.substr(parameter.find("=")+1);
				ss << parameter.c_str();
				ss >> WriteASCIIGrids;
			}
			else if (parameter.find("WriteDoubleGrids") != std::string::npos) {
				std::stringstream ss;
				parameter = parameter.substr(parameter.find("=") + 1);
				ss << parameter.c_str();
				ss >> WriteDoubleGrids;
			}
			else if (parameter.find("WriteNetcdfGrids") != std::string::npos) {
				std::stringstream ss;
				parameter = parameter.substr(parameter.find("=") + 1);
				ss << parameter.c_str();
				ss >> WriteNetcdfGrids;
			}
			else if (parameter.find("InitialStateDirectory") != std::string::npos) {
				std::stringstream ss;
				parameter=parameter.substr(parameter.find("=")+1);
				ss << parameter.c_str();
				DefaultInputDirectory = newstring(_stringsize_+1,1001);
				ss >> DefaultInputDirectory;
			}
			else if (parameter.find("DefaultOutputDirectory") != std::string::npos) {
				std::stringstream ss;
				parameter = parameter.substr(parameter.find("=") + 1);
				ss << parameter.c_str();
				DefaultOutputDirectory = newstring(_stringsize_ + 1, 1001);
				ss >> DefaultOutputDirectory;
			}
			else if (parameter.find("MinutesToInvEndRouting") != std::string::npos) {
				std::stringstream ss;
				parameter = parameter.substr(parameter.find("=") + 1);
				ss << parameter.c_str();
				ss >> MinutesToInvEndRouting_m;
				*MinutesToInvEndRouting = MinutesToInvEndRouting_m;
			}
			else if (parameter.find("timestep")!=std::string::npos) {
				std::stringstream ss;
				parameter=parameter.substr(parameter.find("=")+1);
				ss << parameter.c_str();
				ss >> timestep_m;
				timestep = timestep_m;
				this->timestep = timestep_m;
			}
		}
		else 
		{
			DoContinue=false;
		}
	} while (DoContinue);

	/* looking for the output-list */
	findeabschnitt(_output_list_);
	if (fscanf(steuerdatei,"%d",numsubb)!=1) steuererror("number of subbasins in output list","");
	if (*numsubb>0)
	{
		*listsubb=(int*)calloc(*numsubb+1,sizeof(int));
		for (i=1;i<=*numsubb;i++)
		{
			if (fscanf(steuerdatei,"%d",&((*listsubb)[i]))!=1) steuererror("codes in output-list","");
			fgets(satz,_stringsize_,steuerdatei);
		}
		if (i>1)
		for (i=1;i<*numsubb;i++)
		if ((*listsubb)[i]>=(*listsubb)[i+1])
		{
			(*WSstream) << "FATAL: \nError in control file: order of subbasin codes in output_list!\n"
			"if output should be written for a specified set of subbasins,\n"
			"(specified by a statcode of 5 for averaging in time or 6 for summing in time)\n"
			"than the codes has to be arranged in ascendend order\n."
			"look in [output_list]\n"<<"[end-entry]";
			throw(97);
		}
	}
	else
	listsubb=NULL;


	/* time relevant data like coordinates, meridians */
	findeabschnitt(_Koordinaten_);
	if (fscanf(steuerdatei,"%lf",breite)!=1) steuererror("geographic latitude","");
	fgets(satz,_stringsize_,steuerdatei);
	if (fscanf(steuerdatei,"%lf",laenge)!=1) steuererror("geographic longitude","");
	fgets(satz,_stringsize_,steuerdatei);
	if (fscanf(steuerdatei,"%lf",bezugsmeridian)!=1) steuererror("reference meridian","");
	fgets(satz,_stringsize_,steuerdatei);
	if (fscanf(steuerdatei,"%lf",timekorr)!=1) steuererror("time shift of data time to lokal time","");
	fgets(satz,_stringsize_,steuerdatei);

	/* reading the number of meteo-data */
	findeabschnitt(_anzahl_);
	if (fscanf(steuerdatei,"%d",anzahlgroessen)!=1) steuererror("number of meteorological variables","");

	/* reading names of paragraphs for meteorological input */
	findeabschnitt(_meteonamen_);
	*meteonamen         = (char**)calloc(*anzahlgroessen+1,sizeof(char*));
	for (i=1;i<=*anzahlgroessen;i++)
	{
		(*meteonamen)[i] = (char *)calloc(201,sizeof(char));
		if ((*meteonamen)[i]==NULL) speichererror(1);

		fscanf(steuerdatei,"%s",(*meteonamen)[i]);
		fgets(satz,_stringsize_,steuerdatei);
	}

	if ((fclose(steuerdatei))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
		throw(11);
	}

	satz=tempsatz;
	if (satz) free(satz);satz=NULL;
}


void defaultini::GetGridNames()
{
	int i;

	/* suchen des Abschnittes _Zonengrid_ in der Steuerdatei */
	findeabschnitt(_Zonengrid_);
	gridinfo[1].name = newstring(_stringsize_+1,1);
	fscanfln(steuerdatei,"%s",gridinfo[1].name," Name for ",tstr);
	gridinfo[1].identification = newstring(101,1);
	strcpy(gridinfo[1].identification,_Zonengrid_);

	/* suchen des Abschnittes _Hoehenmodell_ in der Steuerdatei */
	findeabschnitt(_Hoehenmodell_);
	gridinfo[2].name = newstring(_stringsize_+1,1);
	fscanfln(steuerdatei,"%s",gridinfo[2].name," Name for ",tstr);
	gridinfo[2].identification = newstring(101,1);
	strcpy(gridinfo[2].identification,_Hoehenmodell_);

	/* suchen des Abschnittes "Standardgrids" in der Steuerdatei */
	findeabschnitt(_Standardgrids_);

	fscanfln(steuerdatei,"%d",&gridanzahl,"counter/number of grids ",tstr);
	for (i=1;i<=gridanzahl;i++)
	{
		gridinfo[i+2].name = newstring(_stringsize_+1,1);
		gridinfo[i+2].identification = newstring(101,1);

		// the new format of a standard grid entry is (all in one row!):
		// <[path]inputname> <identifier> fillcode = <fillcode> 
		// [readcode = <readcode>] [defaultValue = <defaultValue>] [writecode = <writecode>]
		// [periodicity = <periodicity> <type> <sub-intervals> persistent = <persistency>]
		// [outgrid = <[path]outputname>] [statname = <statisticsname> statcode <statcode>]
		// 
		// if using the new version, the fillcode will not be read correctly, so the error 
		// handling should consider this fact and should try again looking for the new syntax

		bool mustFindFillCode = false;
		unsigned long actpos = ftell(steuerdatei);
		if (fscanf(steuerdatei,"%s%s%d",gridinfo[i+2].name,gridinfo[i+2].identification,&(gridinfo[i+2].fillcode))!=3)
		{
			fseek(steuerdatei,actpos,0);
			if (fscanf(steuerdatei,"%s%s",gridinfo[i+2].name,gridinfo[i+2].identification)!=2)
			{
				sprintf(tstr,"\nerror reading a grid name or identifier %d in\n",i);
				steuererror(tstr,instanzname);
			}
			mustFindFillCode = true;
		}

		// new in version 8.3.x: checking, if there is an additional parameter controlling 
		// periodically reading of standardgrids
		fgets(tstr,_stringsize_,steuerdatei);
		std::string s(tstr);

		// looking for fillcode = <fillcode> default: 0
		if (s.find("fillcode")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("fillcode"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].fillcode;
		}
		else
		{
			if (mustFindFillCode)
			{
				// it's a new style definition but no fillcode = <code> was found
				sprintf(tstr,"\nerror reading the fillcode for standard grid %d (found neither in old nor in new style standard grid definition) in\n",i);
				steuererror(tstr,instanzname);
			}
			// gridinfo[i+2].fillcode = 0; was either set by reading in an old style standard grid or read by searching for fillcode = <code>
		}

		// looking for readcode = <readcode> default: 1
		if (s.find("readcode")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("readcode"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].readcode;
		}
		else
		{
			gridinfo[i+2].readcode = 1;
		}

		// looking for defaultValue = <defaultValue> default: 0
		if (s.find("defaultValue")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("defaultValue"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].defaultValue;
		}
		else
		{
			gridinfo[i+2].defaultValue = 0;
		}

		// looking for writecode = <writecode> default: 0
		if (s.find("writecode")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("writecode"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].writegrid;
		}
		else
		{
			gridinfo[i+2].writegrid = 0;
		}

		// looking for outname = <outname> default: 0
		if (s.find("outname")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("outname"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			gridinfo[i+2].outname = newstring(_stringsize_+1,1);
			ss >> gridinfo[i+2].outname;
		}
		else
		{
			// use a copy of the input name -- thus, the new output defaults to the same as the old one but by using a special output we can implement a new behaviour
			gridinfo[i+2].outname = newstring(_stringsize_+1,1);
			strcpy(gridinfo[i+2].outname,gridinfo[i+2].name);
		}

		// looking for statfile = <statfile> default: 0
		if (s.find("statfile")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("statfile"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			gridinfo[i+2].statname = newstring(_stringsize_+1,1);
			ss >> gridinfo[i+2].statname;
		}
		else
		{
			gridinfo[i+2].statname = NULL;
		}

		// looking for statcode = <statcode> default: 0
		if (s.find("statcode")!=std::string::npos)
		{
			int index1;
			std::string buffer;
			buffer=s.substr(s.find("statcode"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].statcode;
		}
		else
		{
			gridinfo[i+2].statcode = 0;
		}

		if (s.find("periodicity")!=std::string::npos)
		{
			// sample: "periodicity = 1 M 15" --> every (1) month (M) on the 15th day
			int index1;
			gridinfo[i+2].cPerType = new char[100];
			gridinfo[i+2].cPerGridName = newstring(_stringsize_+1,1);
			std::string buffer;
			buffer=s.substr(s.find("periodicity"));
			index1 = buffer.find("=")+1;
			buffer=buffer.substr(index1);
			std::stringstream ss;
			ss << buffer.c_str();
			ss >> gridinfo[i+2].iPer >> gridinfo[i+2].cPerType >> gridinfo[i+2].iSubInt;
			if (s.find("persistent")!=std::string::npos)
			{
				buffer=s.substr(s.find("persistent"));
				index1 = buffer.find("=")+1;
				buffer=buffer.substr(index1);
				std::stringstream ss1;
				ss1 << buffer.c_str();
				ss1 >> gridinfo[i+2].bPersistency >> gridinfo[i+2].iReplaceAfter;
			}
			else
			{
				gridinfo[i+2].bPersistency = true;
			}
			gridinfo[i+2].bReplaceGridNow = false;

		}
		else
		{
			gridinfo[i+2].iPer = -1;
			gridinfo[i+2].iSubInt = -1;
			gridinfo[i+2].cPerType = NULL;
			gridinfo[i+2].cPerGridName = NULL;
			gridinfo[i+2].bPersistency = true;
			gridinfo[i+2].bReplaceGridNow = false;
		}

		// for the moment, we will not allow the feature "periodicity" and the feature "output as grid and/or statistics file" together
//		if ((gridinfo[i+2].iPer >=0)&&((gridinfo[i+2].outname!=NULL) || (gridinfo[i+2].statname !=NULL)))
//		{
//			(*WSstream) << "FATAL: \nError in control file: \n"<<
//							"in section [standard_grids]:\n"<<
//							"feature \"periodicity\" and either one of the features \"writing to statistics\" (using \"statfile = <statfile>\")\n"<<
//							"or \"writing to grids\" (using \"outname = <outname>\") are invalid combined in this version\n"<<
//							"here: entry "<<i<<" ("<<gridinfo[i+2].name<<"/"<<gridinfo[i+2].identification<<")"<<"[end-entry]";
//			throw(198);
//		}

		// initialization of counters for periods and subintervals
		gridinfo[i+2].iPerAct=gridinfo[i+2].iSubIntAct=0;
		gridinfo[i+2].iLastSubInt=gridinfo[i+2].iLastPer=9999; // big numbers, because we will check the new value against the old ones if they are smaller

#ifdef _NETCDF
		// Try to find latitude and longitude grids in the standard grids so
		// we can write them into the netcdf files
		if (strcmp(gridinfo[i + 2].identification, "latitude") == 0) {
			(*WSstream) << "found latitude grid on position: " << i << "[end-entry]";
			//for (int j=1;j<=gridanzahl;j++) {
			gridinfo[1].has_latgrid_ = true;
			gridinfo[1].latgrid_index_ = i + 2;
			//}
		}
		if (strcmp(gridinfo[i + 2].identification, "longitude") == 0) {
			(*WSstream) << "found longitude grid on position: " << i << "[end-entry]";
			//for (int j=1;j<=gridanzahl;j++) {
			gridinfo[1].has_longrid_ = true;
			gridinfo[1].longrid_index_ = i + 2;
			//}
		}
#endif

	}
}

void defaultini::GetDates()
{
	/* Suchen des Abschnitts mit der Modellzeit */
	findeabschnitt(_Modellzeit_);

	if (timestep == 0) // if the optional entry "timestep = n" was not found in section [output_interval], timestep is not yet set and hence the old starttime and endtime is used. Otherwise, the new version must be present in the control file
	{
		int retcode = 0;
		retcode += fscanfln(steuerdatei, "%ld", &stundea, "starting hour in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &taga, "starting day in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &monata, "starting month in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &jahra, "starting year in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &stundee, "end hour in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &tage, "and day in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &monate, "end month in ", tstr);
		retcode += fscanfln(steuerdatei, "%ld", &jahre, "end year in ", tstr);
		if (retcode != 8)
		{
			(*WSstream) << "ERROR: one or more values of start and end time couldn't be read in section " << _Modellzeit_ << "\n" 
				<< "when NOT using the optional parameter timestep in section [output_interval], the start and end date/time must be provided in the legacy format one single parameter per line like this:\n"
				<< _Modellzeit_ << "\n"
				<< "$starthour\n"
				<< "$startday\n" 
				<< "$startmonth\n" 
				<< "$startyear\n" 
				<< "$endhour\n" 
				<< "$endday\n" 
				<< "$endmonth\n" 
				<< "$endyear\n" << "\n(or, of course, by using the actual values without using variables)"<< "[end-entry]";
				
			throw(1300);
		}
	}
	else
	{
		bool DoContinue = true;
		bool bminutea = false, bstundea = false, btaga = false, bmonata = false, bjahra = false, bminutee = false, bstundee = false, btage = false, bmonate = false, bjahre = false;
		do
		{
			fgets(satz, _stringsize_, steuerdatei);
			std::string s(satz);
			s = s.substr(s.find_first_not_of(" \t"));
//			sprintf(satz, "%s", s.c_str());
			if ((!feof(steuerdatei)) && (s.find_first_not_of(" \t\r\n")!=std::string::npos)) // stop at the first empty line
			{
				if (s.find("startminute") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> minutea;
					bminutea = true;
				}
				if (s.find("starthour") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> stundea;
					bstundea = true;
				}
				if (s.find("startday") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> taga;
					btaga = true;
				}
				if (s.find("startmonth") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> monata;
					bmonata = true;
				}
				if (s.find("startyear") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> jahra;
					bjahra = true;
				}
				if (s.find("endminute") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> minutee;
					bminutee = true;
				}
				if (s.find("endhour") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> stundee;
					bstundee = true;
				}
				if (s.find("endday") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> tage;
					btage = true;
				}
				if (s.find("endmonth") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> monate;
					bmonate = true;
				}
				if (s.find("endyear") != std::string::npos) {
					std::stringstream ss;
					s = s.substr(s.find("=") + 1);
					ss << s.c_str();
					ss >> jahre;
					bjahre = true;
				}
			}
			else
			{
				DoContinue = false;
			}
		} while (DoContinue);
		if (!(bminutea && bstundea && btaga && bmonata && bjahra && bminutee && bstundee && btage && bmonate && bjahre))
		{
			(*WSstream) << "ERROR: one or more values of start and end time couldn't be read in section [" << _Modellzeit_ << "]\n"
				<< "when using the optional parameter timestep in section [output_interval], the start and end date/time must be provided in the new format of key = value pairs like this:\n"
				<< _Modellzeit_<<"\n"
				<< "startminute = $startminute\n"
				<< "starthour = $starthour\n"
				<< "startday = $startday\n"
				<< "startmonth = $startmonth\n"
				<< "startyear = $startyear\n"
				<< "endminute = $endminute\n"
				<< "endhour = $endhour\n"
				<< "endday = $endday\n"
				<< "endmonth = $endmonth\n"
				<< "endyear = $endyear\n" << "\n(or, of course, by using the actual values without using variables)" << "[end-entry]";

			throw(1301);
		}

	}

	if ( jahra<1000 ) jahra=jahra+1900;
	if ( jahre<1000 ) jahre=jahre+1900;
	if ( monata<1  )  monata=1;
	if ( monata>12 )  monata=12;
	if ( monate<1  )  monate=1;
	if ( monate>12 )  monate=12;

	if ( ((jahre % 4) == 0) && ((jahre % 100) != 0) )  mon[2]=29; else mon[2]=28;
	if ( tage < 1  ) tage = 1;
	if ( tage > mon[(int)monate])  tage = mon[(int)monate];
	if ( taga < 1  ) taga = 1;
	if ( taga > mon[(int)monata])  taga = mon[(int)monata];
}

void defaultini::GetMultiLayerLUSETable(int version)
{
	int i,usenr;
	/* Einlesen der Tabelle der Landnutzungsspezifischen Daten aus der Steuerdatei */
	if (findeabschnitt(_MultiLayerLUSETable_,false)==true)
	{
		fscanf(steuerdatei,"%d",&usenr);
		tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* Erzeugen und initialisieren der Landnutzungstabelle. Die Werte werden in
		defaultini aus der Steuerdatei eingelesen */
		MultiLayerLUTable = new MultiLayerLUSEType[usenr+1];


		for (i=1;i<=usenr;i++)
		{
			//	  LUtable[i]=new LUSEtabletype(i);
			MultiLayerLUTable[i].setCode(i);
		}
		
		/* im nicht genutzten ersten Element des LUtable-Arrays wird die gesamtzahl an Landnutzungen
		festgehalten um beim spaeteren Zugriff diesen Wert noch zu kennen und um beim Freigeben
		keinen Fehler zu provozieren */
		//  LUtable[0] = new LUSEtabletype(usenr);
		MultiLayerLUTable[0].setCode(usenr);

		for (i=1;i<=usenr;i++)
		{
			tstr=fgets(tstr,_stringsize_,steuerdatei);
			std::string record(tstr);

			if (record.find("Landuse_Layers")!=std::string::npos)
			{
				while (record.find("}")==std::string::npos)
				{
					tstr=fgets(tstr,_stringsize_,steuerdatei);
					record.append(std::string(tstr));
				}
				MultiLayerLUTable[i].initialize(WSstream,record);

				if (MultiLayerLUTable[i].ArcVegLUFractions != NULL)
				{
					if (MultiLayerLUTable[i].getNumFractions() != MultiLayerLUTable[i].getNumLayers())
					{
						(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of landuse layers and ArcVeg_Initial_Fractions has to be identical if ArcVeg_Initial_Fractions are used. \n"
							"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "), number of layers is " << MultiLayerLUTable[i].getNumLayers()
							<< " , but number of fractions for ArcVeg is " << MultiLayerLUTable[i].getNumFractions() << ")" << "[end-entry]";
						throw(111);
					}
				}

				if (MultiLayerLUTable[i].ArcVegBiomass != NULL)
				{
					if (MultiLayerLUTable[i].getNumBiomass() != MultiLayerLUTable[i].getNumLayers())
					{
						(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of landuse layers and ArcVeg_Initial_Biomass has to be identical if ArcVeg_Initial_Biomass is used. \n"
							"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "), number of layers is " << MultiLayerLUTable[i].getNumLayers()
							<< " , but number of biomasses for ArcVeg is " << MultiLayerLUTable[i].getNumBiomass() << ")" << "[end-entry]";
						throw(111);
					}
				}

				if (MultiLayerLUTable[i].ArcVegRetrans != NULL)
				{
					if (MultiLayerLUTable[i].getNumRetrans() != MultiLayerLUTable[i].getNumLayers())
					{
						(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of landuse layers and ArcVeg_Initial_Retrans has to be identical if ArcVeg_Initial_Retrans is used. \n"
							"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "), number of layers is " << MultiLayerLUTable[i].getNumLayers()
							<< " , but number of Retrans for ArcVeg is " << MultiLayerLUTable[i].getNumRetrans() << ")" << "[end-entry]";
						throw(111);
					}
				}

			}
			if (i>1)
			{
				if (MultiLayerLUTable[i].getNumLayers() != MultiLayerLUTable[i - 1].getNumLayers())
				{
					(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of landuse layers has to be equal for all multilayer-landuse types\n"
						"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "): " << MultiLayerLUTable[i].getNumLayers()
						<< " layers, but multilayer-landuse type " << MultiLayerLUTable[i - 1].getCode() << ": "
						<< MultiLayerLUTable[i - 1].getNumLayers() << " layers)\n\n" << "[end-entry]";
					throw(111);
				}

				if (MultiLayerLUTable[i].getNumFractions() != MultiLayerLUTable[i - 1].getNumFractions())
				{
					(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of ArcVeg_Initial_Fractions has to be equal for all multilayer-landuse types if ArcVeg_Initial_Fractions are used.\n"
						"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "): " << MultiLayerLUTable[i].getNumFractions()
						<< " fractions, but multilayer-landuse type " << MultiLayerLUTable[i - 1].getCode() << ": "
						<< MultiLayerLUTable[i - 1].getNumFractions() << " fractions)\n\n" << "[end-entry]";
					throw(111);
				}

				if (MultiLayerLUTable[i].getNumBiomass() != MultiLayerLUTable[i - 1].getNumBiomass())
				{
					(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of ArcVeg_Biomass has to be equal for all multilayer-landuse types if ArcVeg_Initial_Biomass are used.\n"
						"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "): " << MultiLayerLUTable[i].getNumBiomass()
						<< " fractions, but multilayer-landuse type " << MultiLayerLUTable[i - 1].getCode() << ": "
						<< MultiLayerLUTable[i - 1].getNumBiomass() << " biomases)\n\n" << "[end-entry]";
					throw(111);
				}

				if (MultiLayerLUTable[i].getNumRetrans() != MultiLayerLUTable[i - 1].getNumRetrans())
				{
					(*WSstream) << "FATAL: \nError in [multilayer_landuse]: The numbers of ArcVeg_Retrans has to be equal for all multilayer-landuse types if ArcVeg_Initial_Retrans are used.\n"
						"(here: multilayer-landuse type " << MultiLayerLUTable[i].getCode() << " (entry number " << i << "): " << MultiLayerLUTable[i].getNumRetrans()
						<< " fractions, but multilayer-landuse type " << MultiLayerLUTable[i - 1].getCode() << ": "
						<< MultiLayerLUTable[i - 1].getNumRetrans() << " biomases)\n\n" << "[end-entry]";
					throw(111);
				}
			}
		}
	}
	else
	{
		MultiLayerLUTable=NULL;
	}
}

void defaultini::GetLuseTable(int version)
{
	int i,usenr;

	/* Einlesen der Tabelle der Landnutzungsspezifischen Daten aus der Steuerdatei */
	findeabschnitt(_Landnutzungstabelle_);
	fscanf(steuerdatei,"%d",&usenr);
	tstr=fgets(tstr,_stringsize_,steuerdatei);

	/* Erzeugen und initialisieren der Landnutzungstabelle. Die Werte werden in
	defaultini aus der Steuerdatei eingelesen */
	LUtable = new LUSEtabletype[usenr+1];


	for (i=1;i<=usenr;i++)
	{
		//	  LUtable[i]=new LUSEtabletype(i);
		LUtable[i].setCode(i);
	}
	
	/* im nicht genutzten ersten Element des LUtable-Arrays wird die gesamtzahl an Landnutzungen
	festgehalten um beim spaeteren Zugriff diesen Wert noch zu kennen und um beim Freigeben
	keinen Fehler zu provozieren */
	//  LUtable[0] = new LUSEtabletype(usenr);
	LUtable[0].setCode(usenr);

	for (i=1;i<=usenr;i++)
	{
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		std::string record(tstr);

		if ((record.find("VariableDayCount")!=std::string::npos)&&(record.find("method")!=std::string::npos)&&(record.find("{")!=std::string::npos))
		{
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			LUtable[i].initialize(WSstream,record,"VariableDayCount");
		} else if ((record.find("DynamicPhenology_1")!=std::string::npos)&&(record.find("method")!=std::string::npos)&&(record.find("{")!=std::string::npos)) {
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			LUtable[i].initialize(WSstream,record,"DynamicPhenology_1");
		} else if ((record.find("DynamicPhenology_2")!=std::string::npos)&&(record.find("method")!=std::string::npos)&&(record.find("{")!=std::string::npos)) {
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			LUtable[i].initialize(WSstream,record,"DynamicPhenology_2");
		} else if ((record.find("DynamicPhenology_3")!=std::string::npos)&&(record.find("method")!=std::string::npos)&&(record.find("{")!=std::string::npos)) {
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			LUtable[i].initialize(WSstream,record,"DynamicPhenology_3");
		} else if ((record.find("DynamicPhenology_4")!=std::string::npos)&&(record.find("method")!=std::string::npos)&&(record.find("{")!=std::string::npos)) {
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			LUtable[i].initialize(WSstream,record,"DynamicPhenology_4");
		} else	{
			LUtable[i].initialize(WSstream,std::string(tstr),version);
		}
	}

	int JDVegWrite=-1;
	int JDVegReset=-1;
	// looking for JDVegWrite in any land use code. The first one wins and will be set to all other land use codes
	for (i=1;i<=usenr;i++)
	{
		int temp=LUtable[i].getJDVegetationWrite();
		if (temp>0)
		{
			JDVegWrite=temp;
			break;
		}
	}
	if (JDVegWrite<0) JDVegWrite=365; // default: at day 365 the vegetation period dependent grids will be written

	// looking for JDVegReset in any land use code. The first one wins and will be set to all other land use codes
	for (i=1;i<=usenr;i++)
	{
		int temp=LUtable[i].getJDVegetationReset();
		if (temp>0)
		{
			JDVegReset=temp;
			break;
		}
	}
	if (JDVegReset<0) JDVegReset=1; // default: first day of the year is the da, the vegetation period dependent grids will be reset

	// LUtable[0] holds common parameters like the number of landuse codes (in code) and now the write and reset JDs, too
	LUtable[0].setJDVegetationReset(JDVegReset);
	LUtable[0].setJDVegetationWrite(JDVegWrite);

}


void defaultini::GetSoilTable(int version)
{
	int i,usenr;

	/*  --------------------------------------------------------  */
	/* |                                                        | */
	/* | Erzeugen und Einlesen der Werte der Bodenarten-Tabelle | */
	/* |                                                        | */
	/*  --------------------------------------------------------  */

	findeabschnitt(_Bodenartentabelle_);
	fscanf(steuerdatei,"%d",&usenr);
	tstr=fgets(tstr,_stringsize_,steuerdatei);

	/* Erzeugen und initialisieren der Bodenartentabelle. Die Werte werden in
	defaultini aus der Steuerdatei eingelesen */
	SOILtable = new SOILtabletype[usenr+1];
	for (i=1;i<=usenr;i++)
	{
		SOILtable[i].setCode(i);
	}

	/* im nicht genutzten ersten Element des SOILtable-Arrays wird die gesamtzahl an Bodenarten
	festgehalten um beim spaeteren Zugriff diesen WErt noch zu kennen und um beim Freigeben
	keinen Fehler zu provozieren */
	SOILtable[0].setCode(usenr);

	for (i=1;i<=usenr;i++)
	{
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		std::string record(tstr);
		if ((record.find("MultipleHorizons")!=std::string::npos)||(record.find("method")!=std::string::npos)||(record.find("{")!=std::string::npos))
		{
			while (record.find("}")==std::string::npos)
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				record.append(std::string(tstr));
			}
			SOILtable[i].initialize(WSstream,record,"MultipleHorizons");
		}
		else
		{
			SOILtable[i].initialize(WSstream,std::string(tstr),version); // legacy method; not supported any more
		}
		if ((i>1)&&(SOILtable[i].getLevels()!=SOILtable[i-1].getLevels()))
		{
			(*WSstream) << "FATAL: the numbers of soil layers has to be equal for all soil types\n"
			"(here: soil type "<<SOILtable[i].getCode()<<": "<<SOILtable[i].getLevels()
			<<" layers, but soil type "<<SOILtable[i-1].getCode()<<": "
			<<SOILtable[i-1].getLevels()<<" layers)\n\n"<<"[end-entry]";
			throw(30);
		}
	}

	double psi_FC_global = 0;
	int usePsiFCglobalForUnassignedOnly = 0;

	if (findeabschnitt("[GlobalSoilParameters]", false)==true)
	{
		tstr = fgets(tstr, _stringsize_, steuerdatei);
		std::string row(tstr);

		std::string buffer = row.substr(row.find("psi_FC_global"));
		int index1 = buffer.find("=") + 1;
		int index2 = buffer.find(";");
		buffer = buffer.substr(index1, index2 - index1);
		psi_FC_global = atof(buffer.c_str());

		tstr = fgets(tstr, _stringsize_, steuerdatei);
		row = std::string(tstr);

		buffer = row.substr(row.find("usePsiFCglobalForUnassignedOnly"));
		index1 = buffer.find("=") + 1;
		index2 = buffer.find(";");
		buffer = buffer.substr(index1, index2 - index1);
		usePsiFCglobalForUnassignedOnly = atoi(buffer.c_str());

		for (i = 1; i <= usenr; i++) // now set the global default value for all soil table entries
		{
			SOILtable[i].setPsiFCDefault(psi_FC_global, usePsiFCglobalForUnassignedOnly);
		}

	}
	else
	{
		double psi_FC_global = 222;
		int usePsiFCglobalForUnassignedOnly = 0;
	}

}


void defaultini::GetExpokorrIgnore()
{

	/* nun noch Suchen des Abschnittes "Expositionskorrektur",
	um den Ignoranzparameter zu lesen, was spaeter bei der
	Erzeugung von Aspect und Slope genutzt wird. */
	findeabschnitt("[" _Expositionskorrektur_ "]");
	fscanf(steuerdatei,"%d",&do_expokorr);

}

void defaultini::irrigationtable_error(int i)
{
	(*WSstream) << "FATAL: irrigation description for irrigation type "<<i<<" is wrong "
	<< "there are 4 possible ways to define irrigation rules:\n"
	"two schedule variants:\n"
	"1.: \"starting from MM DD with XX mm to MM DD with YY mm every ZZ day\"\n"
	"2.: \"starting from MM DD with XX mm AA times every ZZ days\"\n"
	" where MM and DD are Month and day, XX is the first irrigation amount \n"
	" per irrigation event in mm/day, ZZ is the time step (e.g. every 4 days..),\n"
	" YY is the last irrigation amount in mm/day (it is linearly interpolated)\n"
	" and AA is a repeat counter (only in the secound format).\n"
	" The model recognises both formats by the different formats by analysing the\n"
	" key-words \"starting from\", \"to\" and \"times\"\n"
	" note: there may be defined multiple irrigation schemes in one row,\n"
	" but they have to be of the same type (1. or 2.) \n"
	" and they have to be seperated by a semicolon (;)\n"
	" 3.: irrigation control by demand (1): suction (psi) for start and end of irrigation must be provided\n"
	" 4.: irrigation control by demand (2): ratio between ETR and ETP for irrigation start must be provided\n"<<"[end-entry]";
	throw(29);
}

void defaultini::GetIrrigTable(int version)
{
	int i,j,usenr;
	int startmonth,startday,endmonth,endday,timestep,times,actindex,startindex,endindex;
	double startamount,endamount;
	char *scheme;

	if (version>1)
	{
		//    ***************************************************
		//    *       reading irrigation table                  *
		//    *************************************************** 
		
		findeabschnitt("[" _irrigation_parameters_ "]");
		fscanf(steuerdatei,"%d",&do_irrigation);

		if (do_irrigation)
		{
			/* Einlesen der Tabelle der Landnutzungsspezifischen Daten aus der Steuerdatei */
			findeabschnitt(_irrigation_table_);
			fscanf(steuerdatei,"%d",&usenr);
			
			/* Erzeugen und initialisieren der Bewässerungstabelle. Die Werte werden in
	defaultini aus der Steuerdatei eingelesen */
			irrigationtable = (irrigationtabletype*)calloc(usenr+1,sizeof(irrigationtabletype));
			for (i=1;i<=usenr;i++)
			{
				irrigationtable[i].Name = (char*)calloc(50,sizeof(char));
				irrigationtable[i].Code = i;
			}
			
			/* im nicht genutzten ersten Element des irrigationtable-Arrays wird die gesamtzahl an codes
	festgehalten um beim spaeteren Zugriff diesen Wert noch zu kennen und um beim Freigeben
	keinen Fehler zu provozieren */
			irrigationtable[0].Code = usenr;
			
			for (i=1;i<=usenr;i++)
			{
				/* Einlesen des Codes, des Namens und der Methode */
				if (fscanf(steuerdatei,"%d%50s%d%d",&(irrigationtable[i].Code),irrigationtable[i].Name,&(irrigationtable[i].methode),
							&(irrigationtable[i].origin))!=4)
				{
					(*WSstream) << "FATAL: for irrigation code "<<i<<": code, name, method, or origin could not be read correctly\n"<<"[end-entry]";
					throw(29);
				}
				if (irrigationtable[i].methode==4)
				{
					if (fscanf(steuerdatei,"%lf",&(irrigationtable[i].MinEtr2Etp))!=1)
					{
						(*WSstream) << "FATAL: for irrigation code "<<i<<": min(ETR/ETP) could not be read correctly\n"<<"[end-entry]";
						throw(29);
					}
					// default values, although never really used, except during initialization of unsatzonmodel, so it is save to set a value in the expected range here
					// to avoid accidental crashes by uninitialized floating point values.
					irrigationtable[i].psi_dry=20;
					irrigationtable[i].psi_wet=3.45;
				}
				else
				{
					if (fscanf(steuerdatei,"%lf%lf",&(irrigationtable[i].psi_dry),&(irrigationtable[i].psi_wet))!=2)
					{
						(*WSstream) << "FATAL: for irrigation code "<<i<<": psi_dry or psi_wet could not be read correctly\n"<<"[end-entry]";
						throw(29);
					}
				}
				if (irrigationtable[i].methode==1)
				{
					if (fscanf(steuerdatei,"%d",&(irrigationtable[i].tabcnt))!=1)
					{
						(*WSstream) << "FATAL: for irrigation code "<<i<<": tabcnt could not be read correctly\n"<<"[end-entry]";
						throw(29);
					}
					irrigationtable[i].month  = (long*)calloc(irrigationtable[i].tabcnt+1,sizeof(long));
					irrigationtable[i].day    = (long*)calloc(irrigationtable[i].tabcnt+1,sizeof(long));
					irrigationtable[i].amount = (double*)calloc(irrigationtable[i].tabcnt+1,sizeof(double));
					/* Einlesen der amounts zu den Terminen */
					for (j=1;j<=irrigationtable[i].tabcnt;j++)
					{ 
						if (fscanf(steuerdatei,"%ld%ld%lf",&(irrigationtable[i].month[j]),&(irrigationtable[i].day[j]),&(irrigationtable[i].amount[j]))!=3)
						{
							(*WSstream) << "FATAL: for irrigation code "<<i<<" month["<<j<<"] or day["<<j<<"] or amount["<<j<<"] could not be read\n"<<"[end-entry]";
							throw(29);
						}
					} 
				}
				else
				{
					if (irrigationtable[i].methode==2) 
					/* verbal description of irrigation schedule in 2 different formats: 
			1: "starting from MM DD with XX mm to MM DD with YY mm every ZZ days" here, the start end end date are explicitly given 
			2: "starting from MM DD with XX mm YY times every ZZ days" Here, the number of irrigation events is given explicitly */
					{
						scheme=(char*)calloc(_stringsize_+1,sizeof(char));
						scheme=fgets(scheme,_stringsize_,steuerdatei); /* now the irrigation scheme description is read in in scheme */
						actindex=0; /* first index is 1 for the irrigation table */
						irrigationtable[i].month  = (long*)calloc(2,sizeof(long));
						irrigationtable[i].day    = (long*)calloc(2,sizeof(long));
						irrigationtable[i].amount = (double*)calloc(2,sizeof(double));
						if (!strstr(scheme,"starting from")) irrigationtable_error(i);
						while (strstr(scheme,"starting from")) /* there may be more then one irrigation schemes for a single landuse (e.g. if the amount are increasing and then decreasing) */
						{
							actindex++;
							if (!strstr(scheme,"to") && !strstr(scheme,"times")) /* checking, if there are different types of irrigation schemes defines in one single line */
							irrigationtable_error(i);
							/* raeding startdate */
							if (sscanf(strstr(scheme,"starting from")+strlen("starting from")+1,"%d%d",&startmonth,&startday)!=2) 
							irrigationtable_error(i);
							if (mon[startmonth]<startday) startday=mon[startmonth];
							if (startday<1) startday=1;
							if (startmonth<1) startmonth=1;
							if (startmonth>12) startmonth=12;
							sprintf(scheme,"%s",strstr(scheme,"with"));
							if (strstr(scheme,"with"))
							{ /* rading startamount */
								if (sscanf(strstr(scheme,"with")+strlen("with")+1,"%lf",&startamount)!=1) 
								irrigationtable_error(i);
								/* now two possibilities are open: 1. a given enddate and a timestep or
				2.: a given repeat count and a tiemstep (without altering irrigation amount) */
								if (strstr(scheme,"to")) /* a given enddate -> to */
								{
									sprintf(scheme,"%s",strstr(scheme,"to"));
									if (sscanf(strstr(scheme,"to")+strlen("to")+1,"%d%d",&endmonth,&endday)!=2) 
									irrigationtable_error(i);
									if (mon[endmonth]<endday) endday=mon[endmonth];
									if (endday<1) endday=1;
									if (endmonth<1) endmonth=1;
									if (endmonth>12) endmonth=12;
									if (!strstr(scheme,"with")) irrigationtable_error(i);
									sprintf(scheme,"%s",strstr(scheme,"with")+strlen("with")+1);
									sscanf(scheme,"%lf",&endamount);
									if (!strstr(scheme,"every")) irrigationtable_error(i);
									sprintf(scheme,"%s",strstr(scheme,"every")+strlen("every")+1);
									sscanf(scheme,"%d",&timestep);
									/* now filling the rrigation table using an internal calendar (without 02/29) */
									irrigationtable[i].month[actindex] =startmonth;
									irrigationtable[i].day[actindex]   =startday;
									irrigationtable[i].amount[actindex]=startamount;
									startindex=actindex;
									while ((irrigationtable[i].month[actindex]!=endmonth)||((irrigationtable[i].day[actindex]<=endday)&&(irrigationtable[i].month[actindex]==endmonth)))
									{
										irrigationtable[i].tabcnt = actindex+1; /* first: there are only two known entries (start and end date) */
										irrigationtable[i].month  = (long*)realloc(irrigationtable[i].month,(irrigationtable[i].tabcnt+1)*sizeof(long));
										irrigationtable[i].day    = (long*)realloc(irrigationtable[i].day,(irrigationtable[i].tabcnt+1)*sizeof(long));
										irrigationtable[i].amount = (double*)realloc(irrigationtable[i].amount,(irrigationtable[i].tabcnt+1)*sizeof(double));
										actindex++;
										irrigationtable[i].day[actindex]=irrigationtable[i].day[actindex-1]+timestep;
										if ((irrigationtable[i].day[actindex]>endday)&&(irrigationtable[i].month[actindex-1]==endmonth)) 
										{
											irrigationtable[i].tabcnt--;
											actindex--;
											break;
										}
										if (irrigationtable[i].day[actindex]>mon[irrigationtable[i].month[actindex-1]])
										{
											irrigationtable[i].day[actindex]-=mon[irrigationtable[i].month[actindex-1]];
											irrigationtable[i].month[actindex]=irrigationtable[i].month[actindex-1]+1;
											if (irrigationtable[i].month[actindex]>12) irrigationtable[i].month[actindex]=1;
										}
										else
										irrigationtable[i].month[actindex]=irrigationtable[i].month[actindex-1];
										if ((irrigationtable[i].day[actindex]>endday)&&(irrigationtable[i].month[actindex]==endmonth)) 
										{
											actindex--;
											break;
										}
									}
									endindex=actindex;
									irrigationtable[i].amount[endindex]=endamount;
									/* now interpolating the irrigation amounts */
									for (j=startindex+1;j<endindex;j++)
									irrigationtable[i].amount[j]=((double)j-(double)startindex)/((double)endindex-(double)startindex)*irrigationtable[i].amount[endindex]+
									((double)endindex-(double)j)/((double)endindex-(double)startindex)*irrigationtable[i].amount[startindex];
								}
								else /* a given number of irrigation events -> times*/
								{ 
									if (strstr(scheme,"times"))
									{
										endamount=startamount;
										if (!strstr(scheme,"mm")) irrigationtable_error(i);
										sprintf(scheme,"%s",strstr(scheme,"mm")+strlen("mm")+1);
										sscanf(scheme,"%d",&times);
										if (!strstr(scheme,"every")) irrigationtable_error(i);
										sprintf(scheme,"%s",strstr(scheme,"every")+strlen("every")+1);
										sscanf(scheme,"%d",&timestep);
										/* now filling the rrigation table using an internal calendar (without 02/29) */
										irrigationtable[i].tabcnt = actindex+times-1; /* first: there are only two known entries (start and end date) */
										irrigationtable[i].month  = (long*)realloc(irrigationtable[i].month,(irrigationtable[i].tabcnt+1)*sizeof(long));
										irrigationtable[i].day    = (long*)realloc(irrigationtable[i].day,(irrigationtable[i].tabcnt+1)*sizeof(long));
										irrigationtable[i].amount = (double*)realloc(irrigationtable[i].amount,(irrigationtable[i].tabcnt+1)*sizeof(double));
										irrigationtable[i].month[actindex] =startmonth;
										irrigationtable[i].day[actindex]   =startday;
										irrigationtable[i].amount[actindex]=startamount;
										startindex=actindex;
										for (j=2;j<=times;j++)
										{
											actindex++;
											irrigationtable[i].day[actindex]=irrigationtable[i].day[actindex-1]+timestep;
											if (irrigationtable[i].day[actindex]>mon[irrigationtable[i].month[actindex-1]])
											{
												irrigationtable[i].day[actindex]-=mon[irrigationtable[i].month[actindex-1]];
												irrigationtable[i].month[actindex]=irrigationtable[i].month[actindex-1]+1;
												if (irrigationtable[i].month[actindex]>12) irrigationtable[i].month[actindex]=1;
											}
											else
											irrigationtable[i].month[actindex]=irrigationtable[i].month[actindex-1];
										}
										endindex=actindex;
										irrigationtable[i].amount[endindex]=endamount;
										/* now interpolating the irrigation amounts */
										for (j=startindex+1;j<endindex;j++)
										irrigationtable[i].amount[j]=((double)j-(double)startindex)/((double)endindex-(double)startindex)*irrigationtable[i].amount[endindex]+
										((double)endindex-(double)j)/((double)endindex-(double)startindex)*irrigationtable[i].amount[startindex];
									}
									else
									{
										(*WSstream) << "FATAL: none of the keywords \"to\" or \"times\" could be found!\n"<<"[end-entry]";
										irrigationtable_error(i);
									}
								}
							}
							else irrigationtable_error(i);
						}
						if (scheme) free(scheme);
					}
					else
					{
						irrigationtable[i].month  = NULL;
						irrigationtable[i].day    = NULL;
						irrigationtable[i].amount = NULL;
					}
				}
			}
		}
	}
}

void defaultini::GetGridDimensions(int *p_ncols, int *p_nrows, double *p_xll, double *p_yll, double *p_cellsize, double *p_nodata)
{
	int i, j, format;
	FILE *datei;
	char *FormatTest;

	signal(SIGINT, SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif
	/* Oeffnen des grids und lesen der Parameter */
	if ((datei = fopen(gridinfo[1].name, "rb")) == NULL)
	{
		(*WSstream) << "FATAL: error opening " << gridinfo[1].name << "[end-entry]";
		throw(21);
	}
	FormatTest = (char*)calloc(100, sizeof(char));
	fread(FormatTest, sizeof(char), 4, datei);
	FormatTest[4] = '\0';

	if (strcmp(FormatTest, "ncol") == 0)
	{
		format = 1;
	}
	else
	{
		// analyse binary format: is it a double or a double grid?
		fseek(datei, 0, 0);
		fread(&zz[1], sizeof(double), 12, datei);
		if (zz[6] == -9999) {
			format = 2; // double
		}
		else {
			float fzz[7];
			fseek(datei, 0, 0);
			fread(&fzz[1], sizeof(float), 6, datei);
			if (fzz[6] == -9999)
			{
				format = 0; // float
			}
			else
			{
				(*WSstream) << "FATAL: could not recognise file format (ascii, float or double) for " << gridinfo[1].name << "[end-entry]";
				throw(22);
			}
		}
	}

	fseek(datei, 0, 0);
	// rewind(datei); not sure if this is really ANSII?

	if (format == 1)	{ // ascii
		for (i = 1; i <= 6; i++)
			if (fscanf(datei, "%*s%lf", &(zz[i])) != 1)
			{
				(*WSstream) << "FATAL: error reading " << gridinfo[1].name << "[end-entry]";
				throw(22);
			}
	}
	else if (format == 2)	{ // double
		if ((fread(&zz[1], sizeof(double), 12, datei)) != 12)
		{
			(*WSstream) << "FATAL: error reading double value from the header of " << gridinfo[1].name << "[end-entry]";
			throw(22);
		}
	}
	else { // float (compatibility mode for old grids)
		float fzz[13];
		for (i = 1; i <= 12; i++)
		{
			if ((fread(&(fzz[i]), sizeof(float), 1, datei)) != 1)
			{
				(*WSstream) << "FATAL: error reading double value #" << i << " from the header of " << gridinfo[1].name << "[end-entry]";
				throw(22);
			}
			zz[i] = (double)fzz[i];
		}
	}


	*p_ncols = (int)floor(zz[1]);
	*p_nrows = (int)floor(zz[2]);
	*p_xll = zz[3];
	*p_yll = zz[4];
	*p_cellsize = zz[5];
	*p_nodata = zz[6];


	if (fclose(datei) == EOF)
	{
		(*WSstream) << "FATAL: error closing grid " << gridinfo[1].name << "[end-entry]";
		throw(23);
	}
	SetSignalDefaults();
	if (FormatTest) free(FormatTest);
	return;
}


void defaultini::zonenanalyse(gridinfotype **p_gridinfo, double **RelCellArea,
int *p_ncols, int *p_nrows, double *p_nodata)
{
	int i,j;
	int        zonenid;
	statrecord tempzone;
	unsigned short int  weiter;
	int        *zzarray;
	statrecord *zzone;
	double cellcount;

	/* zone, zarray und zonmax sind lokale Variablen. Ihr Inhalt wird
an die Aufrufparameter mit denselben Namen mit fuerendem Unterstrich
am Ende des defaultini::constructors per Adressuebergabe uebergeben.
Deshalb wird diese Routine ohne Parameter aufgerufen. */

	if (mpi_rank == 0)
	{
		(*WSstream) << "ranking zones from " << (*p_gridinfo)[1].name << " for " << instanzname << "[end-entry]";
		zzone = (statrecord *)calloc(maxzonenzahl + 1, sizeof(statrecord));  /* maxzonenzahl zones (16000), index 0 unused */
		if (!zzone) speichererror(1);
		zzarray = (Tzarray *)calloc((size_t)maxzonenwert + 1, sizeof(Tzarray)); /* zones codes up to max. 32767 (maxzonenwert) */
		if (!zzarray) speichererror(1);
		FirstLastCol = (firstlast *)calloc(nrows + 1, sizeof(firstlast));
		zzone[0].n = 0.0;
		zzone[0].sx = 0.0;
		zzone[0].RelCellSum = 0.0;
		zonmax = 0;
		zonmaxwert = 0;
		for (i = 1; i < maxzonenwert; i++)   zzarray[i] = 0;
		for (i = 1; i <= *p_nrows; i++)
		{
			FirstLastCol[i].First = 0;
			FirstLastCol[i].Last = -1;
			for (j = 1; j <= *p_ncols; j++)
			{
				if (zn[i][j] != *p_nodata)
				{
					if (FirstLastCol[i].First == 0) FirstLastCol[i].First = j;
					FirstLastCol[i].Last = j;
					zonenid = (int)(zn[i][j] + 0.5);
					if ((zonenid<0) || (zonenid>maxzonenwert))
					{
						(*WSstream) << "FATAL: invalid zone code in " << gridinfo[2].name
							<< "!\n (codes has to be between 0 and " << maxzonenwert << ")!" << "[end-entry]";
						throw(3);
					}
					if (zzarray[zonenid] == 0)
					{ /*gibt es die Zone mit der ID zonenid noch nicht? */
						zonmax = zonmax + 1;           /* nein, also in die Liste einfuegen */
						if (zonmax > maxzonenzahl)
						{
							(*WSstream) << "FATAL: too much zones in " << gridinfo[2].name << ", (max. " << maxzonenzahl << ")!" << "[end-entry]";
							throw(4);
						}
						if (zonenid > zonmaxwert) zonmaxwert = zonenid;
						zzarray[zonenid] = zonmax;
						zzone[zonmax].zonenwert = zn[i][j];
						zzone[zonmax].n = 0;
						zzone[zonmax].sx = 0;
						zzone[zonmax].RelCellSum = 0;
					}
				}
			}
		}
		/* Ordnen der Zonen nach dem Zonenwert */
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

		/* zone und zarray brauchen nur so gross zu sein, wie Zonen da sind,
		also findet ein Neubelegen von Speicher statt  */

		zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
		if (!zone) speichererror(1);
		for (i = 0; i <= zonmax; i++)
		{
			zone[i].n = zzone[i].n;
			zone[i].sx = zzone[i].sx;
			zone[i].zonenwert = zzone[i].zonenwert;
			zone[i].rel = zzone[i].rel;
			zone[i].RelCellSum = zzone[i].RelCellSum;
		}
		if (zzone) free(zzone); zzone = NULL;
		/* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */

		zarray = (Tzarray *)calloc(zonmaxwert + 1, sizeof(Tzarray));
		if (!zarray) speichererror(1);
		for (i = 0; i <= zonmaxwert; i++) zarray[i] = zzarray[i];
		if (zzarray) free(zzarray); zzarray = NULL;

		/* Ermitteln des relativen Flaechenanteils der Zonen an der gesamten Modellflaeche */
		/* zunaechst Zonengrid durchsuchen und die register jeder Zone auffuellen (mit 1) */
		for (i = 1; i <= nrows; i++)
			for (j = 1; j <= ncols; j++)
			{
				if (zn[i][j] != nodata)
				{
					zone[zarray[(int)(zn[i][j])]].sx++;
					if (RelCellArea)
					{
						zone[zarray[(int)(zn[i][j])]].n += RelCellArea[i][j];
						zone[zarray[(int)(zn[i][j])]].RelCellSum += RelCellArea[i][j];
					}
					else
					{
						zone[zarray[(int)(zn[i][j])]].n++;
						zone[zarray[(int)(zn[i][j])]].RelCellSum++;
					}
				}
			}

		/* Feststellen der Gesamtzahl an gueltigen Werten */
		cellcount = 0;
		for (i = 1; i <= zonmax; i++)
		{
			cellcount += zone[i].n;
		}

		if (cellcount == 0)
		{
			(*WSstream) << "FATAL: all cells in the zone grid " << gridinfo[2].name << " are NODATA!" << "[end-entry]";
			throw(7);
		}

		/* berechnen der relativen Anteile pro Zone */
		for (i = 1; i <= zonmax; i++)
		{
			zone[i].rel = zone[i].n / cellcount;
		}
		zone[0].rel = cellcount*cellsize*cellsize / 1000000.0; /* in zone[0].rel: area of the entire basin */
	}

// only if there are more than one instances, we need to send the zone structure to other ranks
#ifdef _MPI
	if (mpi_size > 1)
	{
		MPI_Bcast(&zonmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&zonmaxwert, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (mpi_rank > 0)
		{
			FirstLastCol = (firstlast *)calloc(nrows + 1, sizeof(firstlast));
			zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			if (!zone) speichererror(1);
			zarray = (Tzarray *)calloc(zonmaxwert + 1, sizeof(Tzarray));
			if (!zarray) speichererror(1);
		}
		MPI_Bcast(&(zarray[0]), zonmaxwert+1, MPI_INT, 0, MPI_COMM_WORLD);
		int cnt = int(sizeof(statrecord) / sizeof(double))*(zonmax + 1);

		MPI_Datatype zonevectortype;
		MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
		MPI_Type_commit(&zonevectortype);
		for (int i = 0; i <= zonmax; i++)
		{
			MPI_Bcast(&(zone[i]), 1, zonevectortype, 0,MPI_COMM_WORLD);
		}
		MPI_Type_free(&zonevectortype);

		MPI_Datatype firstlastcoltype;
		MPI_Type_contiguous(2,MPI_INT,&firstlastcoltype);
		MPI_Type_commit(&firstlastcoltype);
		for (int i = 0; i <= nrows; i++)
		{
			MPI_Bcast(&(FirstLastCol[i]), 1, firstlastcoltype, 0,MPI_COMM_WORLD);
		}
		MPI_Type_free(&firstlastcoltype);

		if (mpi_rank > 0)
		{
			for (int i = 0; i < this->mpistartrow; i++)
			{
				FirstLastCol[i].First = 0;
				FirstLastCol[i].Last = -1;
			}
			for (int i = this->mpiendrow+1; i <= nrows; i++)
			{
				FirstLastCol[i].First = 0;
				FirstLastCol[i].Last = -1;
			}
		}

	}
#endif

}


void defaultini::check_LU_Codes(int l)
{
	int i,j,k,flag;
	if (l==0) return;

	if (mpi_rank == 0)
	{
		if (MultiLayerLUTable == NULL)
		{

			/* Umwandeln der Codes in LUSE in die Nummern 1..LUtable[0].Code, wie sie in der
			Steuerdatei als Reihenfolge der Landnutzungen in der Landnutzungstabelle
			auftreten. Die effektive Reihenfolge ist nicht wichtig, da die Parameter
			jeder Nutzungsform in einem eigenen, durch den Code 1..LUtable[0].Code zu
			adressierenden Record abgelegt sind (in LUtable), jedoch duerfen nur die
			Codes 1..LUtable[0].Code im Landnutzungsgrid auftreten, da nur soviele
			Records reserviert sind (erster Parameter im Abschnitt Landnutzungstabelle).
			Wird in der Modellierung auf eine bestimmte Landnutzung ein bestimmter
			Algorithmus angewandt, so muss anhand des Feldes "Name" des LUtable-arrays
			entschieden werden, um welche Landnutzung es sich handelt. Die vorgegebenen
			Landnutzungsformen in der Steuerdatei duerfen deshalb nicht veraendert
			werden, soweit die Namen im Modell gebraucht werden. Das betrifft zur Zeit
			nur die Landnutzungsform "Nadelwald". */
			LUSE[0] = gridinfo[l].pointer;
			double **TEMP = matrix(1, nrows, 1, ncols);

			(*WSstream) << "1.: checking codes in land use grid" << "[end-entry]";
			for (i = 1; i <= nrows; i++)
			{
				//ARNOLD2006_EXTEND TO FULL AREA!
				// Reason: In other function (constructor unsatuzr - Constructor, after line:
				//    after initialize_grid(i,_SBLOW_grid_,SBLOW,NULL,0,_SBLOW_grid_,0);
				//  the translation is called for COMPLETE area, not only for the cut-out part!
				for (j = 1; j <= ncols; j++)
					//		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					TEMP[i][j] = LUSE[0][i][j]; /* Erklaerung analog zur Erklaerung in der naechsten Routine check_SOIL_Codes() */
					if (zn[i][j] > nodata)
					{
						flag = 0;
						for (k = 1; k <= LUtable[0].getCode(); k++)
						{
							if ((int)LUSE[0][i][j] == LUtable[k].getCode()) flag = 1;
						}
						if (!flag)
						{
							(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
								"land use " << LUSE[0][i][j] << " in grid " << gridinfo[l].name
								<< " (" << gridinfo[l].identification << ")\n"
								"is not present in the land use table in the control file \n"
								"-----------------------------------------------------------------------\n\n" << "[end-entry]";
							throw(31);
						}
					}
				}
			}

			(*WSstream) << "2.:replacing codes in land use grid " << gridinfo[l].name << "[end-entry]";
			for (k = 1; k <= LUtable[0].getCode(); k++)
			{
				/* Jede Landnutzungsform wird abgefragt */
				if (LUtable[k].getCode() != k)
					/* Wenn die Zeilennummer (Index im LUtable-array) nicht mit
				dem Code fuer diese Nutzung in LUSE uebereinstimmt
				wird im LUSE-Grid Code durch k ersetzt */
				{
					(*WSstream) << LUtable[k].getCode() << " replaced by " << k << " ..." << "[end-entry]";
					for (i = 1; i <= nrows; i++)
						for (j = 1; j <= ncols; j++)
							//			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
							if (TEMP[i][j] == LUtable[k].getCode())
								LUSE[0][i][j] = k;
				}
			}
			free_matrix(TEMP, 1, nrows, 1, ncols);
			TEMP = NULL;
		}
		else // multiple layers for LUSE, so a set of grids holding parameters for all layers must be 
			// created. Also, not only the codes in the landuse table but also those in the multilayer-landuse-table
			// has to be checked for matching codes in the LUSE grid.
		{
			LUSE[0] = gridinfo[l].pointer;
			double **TEMP = matrix(1, nrows, 1, ncols);

			(*WSstream) << "1.: checking codes in land use grid" << "[end-entry]";
			for (i = 1; i <= nrows; i++)
			{
				for (j = 1; j <= ncols; j++)
					//		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					TEMP[i][j] = LUSE[0][i][j]; /* Erklaerung analog zur Erklaerung in der naechsten Routine check_SOIL_Codes() */
					if (zn[i][j] > nodata)
					{
						flag = 0;
						for (k = 1; k <= MultiLayerLUTable[0].getCode(); k++)
						{
							if ((int)LUSE[0][i][j] == MultiLayerLUTable[k].getCode()) flag = 1;
						}
						if (!flag)
						{
							(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
								"land use " << LUSE[0][i][j] << " in grid " << gridinfo[l].name
								<< " (" << gridinfo[l].identification << ")\n"
								"is not present in the multilayer-landuse-table in the control file \n"
								"-----------------------------------------------------------------------\n\n" << "[end-entry]";
							throw(311);
						}
					}
				}
			}
			(*WSstream) << "2.:replacing codes in land use grid " << gridinfo[l].name << "[end-entry]";
			for (k = 1; k <= MultiLayerLUTable[0].getCode(); k++)
			{
				/* Jede Landnutzungsform wird abgefragt */
				if (MultiLayerLUTable[k].getCode() != k)
					/* Wenn die Zeilennummer (Index im LUtable-array) nicht mit
				dem Code fuer diese Nutzung in LUSE uebereinstimmt
				wird im LUSE-Grid Code durch k ersetzt */
				{
					(*WSstream) << MultiLayerLUTable[k].getCode() << " replaced by " << k << " ..." << "[end-entry]";
					for (i = 1; i <= nrows; i++)
					{
						for (j = 1; j <= ncols; j++)
							//			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						{
							if (TEMP[i][j] == MultiLayerLUTable[k].getCode())
								LUSE[0][i][j] = k;
						}
					}
					// ? MultiLayerLUTable[k].setCode(k); ?
				}
			}
			free_matrix(TEMP, 1, nrows, 1, ncols);
			TEMP = NULL;

			// no check the MultiLayerLUTable codes for each layer against all codes in LUtable -> 
			// there may be more codes in LUTable, but at least all codes appearing in MultiLayerLUTable 
			// must be present in LUtable, too
			int code;
			for (i = 1; i <= MultiLayerLUTable[0].getCode(); i++)
			{
				for (j = 1; j <= MultiLayerLUTable[i].getNumLayers(); j++)
				{
					code = MultiLayerLUTable[i].getLULayer(j);
					if (code > nodata)
					{
						flag = 0;
						for (k = 1; k <= LUtable[0].getCode(); k++)
						{
							if (code == LUtable[k].getCode())
							{
								flag = 1;
								break;
							}
						}
						if (!flag)
						{
							(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
								"Error: not all layer codes from landuse " << MultiLayerLUTable[i].getCode() <<
								" from MultiLayer_Landuse-Table\n"
								"are present in the landuse-table in the control file: Here, code " <<
								MultiLayerLUTable[i].getLULayer(j) << " is missing.\n"
								"-----------------------------------------------------------------------\n\n" << "[end-entry]";
							throw(312);
						}
					}
					else
					{
						if (j == 1)
						{
							(*WSstream) << "FATAL: \n\n-----------------------------------------------------------------------\n"
								"Error: the first Layer code of each MultiLayer-Landuse Entry \n"
								"must be a valid entry (not 9999 or larger):\nHere, line " << i << " contains "
								<< MultiLayerLUTable[i].getLULayer(j) << " as first layer code.\n"
								"-----------------------------------------------------------------------\n\n" << "[end-entry]";
							throw(312);
						}
					}
				}
			}

			// make a reguar sequence from the LUTable-Codes (i.e. 1...n step 1) in order to be able to adress the LUTable-entries 
			// by their index
			for (int k = 1; k <= LUtable[0].getCode(); k++)
			{
				// for each LUTable-entry the entire MultiLayer-Landuse-table will be scanned for codes which are not in correct sequence
				int code = LUtable[k].getCode();
				if (code != k)
				{
					// if a LUTable code is not in it's correct sequence, then we replace the Layer code by the corresponding index into the LUTable
					for (int i = 1; i <= MultiLayerLUTable[0].getCode(); i++)
					{
						for (int j = 1; j <= MultiLayerLUTable[i].getNumLayers(); j++)
						{
							if (MultiLayerLUTable[i].getOriLULayer(j) == code) // use ORI-Code in order to avoid multiple changes
							{
								MultiLayerLUTable[i].setLULayer(j, k);
								(*WSstream) << "Code replacement in MultiLayer-Landuse " << i << ": changed Layer-Code #" << j << " from " << code <<
									" to " << k << "[end-entry]";
							}
						}
					}
				}
			}
		}
	}
	// MultiLayerLUTable needs to be broadcasted to all ranks
#ifdef _MPI
	for (int i = 1; i <= MultiLayerLUTable[0].getCode(); i++)
	{
		MPI_Bcast(&(MultiLayerLUTable[i].LULayers[0]), MultiLayerLUTable[i].getNumLayers()+1, MPI_INT, 0, MPI_COMM_WORLD);
	}
#endif
	// as well as the changed land use codes
	MPIScatterGrid(l);
	MPISynchronizeGridHalo(l);
}

// gereates LUSE stack with codes for LUSE table for every landuse layer from multi layer landuese table
void defaultini::generateMultiLayerLanduse()
{
	delete[] LUSE;
	int cnt=MultiLayerLUTable[1].getNumLayers();
	LUSE = new double**[cnt+1]; // index 0 holds the original LUSE grid as it is used for the old compatibility mode
	LUSE[0] = gridinfo[ist_vorhanden(_Landnutzung_)].pointer;
	if (LUSE[0]!=NULL)
	{
		// suchen des ersten freien Grids im gridinfo-array
		int k=1;
		while (gridinfo[k].name!=NULL)
		{ 
			k++;
			if (k==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (at position "<<k<<")\n"<<"[end-entry]";
				throw(12);
			}
		}
		char *layername;
		layername=new char[200];
		for (int z=1;z<=cnt;z++)
		{
			int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
			int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
			if (mpi_rank == 0)
			{
				sr = 1;
				er = nrows;
			}

			sprintf(layername,"%s%d",_Landnutzung_,z);
			initialize_grid(k, layername, NULL, NULL, 0, layername, 0, sr, er);
			gridinfo[k].pointer = matrix(sr,er,1,ncols,nrows);
			LUSE[z] = gridinfo[k].pointer; 
			for (int i=sr;i<=er;i++)
			{
				for (int j=1;j<=ncols;j++)
				//			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata) && (mpi_rank==MPIREG[i][j]))
					{
						LUSE[z][i][j]=MultiLayerLUTable[(int)LUSE[0][i][j]].getLULayer(z);
					}
					else LUSE[z][i][j] = nodata;
				}
			}
#ifdef _MPI
			MPISynchronizeGridHalo(k, LUSE[z]);
#endif
			#ifdef _WASIMDEBUG_
				char *aktname = new char[_stringsize_ + 1];
				sprintf(aktname, "%s_%d", layername, mpi_rank);
				printgridfile(aktname, LUSE[z], WriteASCIIGrids, WriteDoubleGrids,WriteNetcdfGrids);
				delete[] aktname;
			#endif
			gridanzahl++;
			k++;
		}
		delete[] layername;
	}
}

void defaultini::generateMultiLayerParameterGrids(int numlayers)
/* erzeugt aus der Landnutzung die folgenden Grids:
Wurzeltiefe			WURZ, 
Vegetationsbedeckungsgrad VBG,
Blattflaechenindex	LAI,
SurfaceResistance	RSC,
SurfaceResistance	RSI,
Rahigkeitslänge		Z0
*/
{
	int i,j;
	(*WSstream) << "1"<<"[end-entry]";

	if (LUSE)
	{
		if (MultiLayerLUTable[1].getNumFractions() > 0)
		{
			ARCVEGFRACTIONS = new double**[numlayers + 1];
			arcvegfractionsindex = new int[numlayers + 1];
		}
		else
		{
			ARCVEGFRACTIONS = NULL;
			arcvegfractionsindex = NULL;
		}

		if (MultiLayerLUTable[1].getNumBiomass() > 0)
		{
			AVBIOM = new double** [numlayers + 1];
			arcvegbiomassindex = new int[numlayers + 1];
		}
		else
		{
			AVBIOM = NULL;
			arcvegbiomassindex = NULL;
		}

		if (MultiLayerLUTable[1].getNumRetrans() > 0)
		{
			AVRETRANS = new double** [numlayers + 1];
			arcvegretransindex = new int[numlayers + 1];
		}
		else
		{
			AVRETRANS = NULL;
			arcvegretransindex = NULL;
		}

		WURZ      = new double**[numlayers+1];
		VBG       = new double**[numlayers+1]; 
		LAI       = new double**[numlayers+1];
		RSC       = new double**[numlayers+1];
		RSI       = new double**[numlayers+1];
		Z0        = new double**[numlayers+1];
		VEGSTART  = new double**[numlayers+1];
		VEGSTOP   = new double**[numlayers+1];
		VEGDURAT  = new double**[numlayers+1];

		wurzindex	  = new int[numlayers+1];
		vbgindex	  = new int[numlayers+1];
		laiindex	  = new int[numlayers+1];
		rscindex	  = new int[numlayers+1];
		rsiindex	  = new int[numlayers+1];
		z0index		  = new int[numlayers+1];
		vegstartindex = new int[numlayers+1];
		vegstopindex  = new int[numlayers+1];
		vegduratindex = new int[numlayers+1];

		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;

		iNumMTStartLevels=0;

		for (int luz=1;luz<=LUtable[0].getCode();luz++)
		{
			if (iNumMTStartLevels<LUtable[luz].getNumFixPoints())
			{
				iNumMTStartLevels=LUtable[luz].getNumFixPoints();
			}
			if ((LUtable[luz].GetMethod()==3)||(LUtable[luz].GetMethod()==4)||(LUtable[luz].GetMethod()==6))
			{
				// bUseForcAndTStart is a flag used later (also in run-mode) for switching on some special 
				// initializations and calculations: calculating the average temperature for each day as well as 
				// reserve memeory for all of the extra grids and calling the appropriate methods is controlled by this flag
				bUseForcAndTStart=true;
			}
			if ((LUtable[luz].GetMethod()==5)||(LUtable[luz].GetMethod()==6))
			{
				// bUseMTStart indicates the usage of dynamic phenology 3 in at least one land use table entry
				// that means extra grids (resp. stacks) has to be reserved and managed which is flagged by this boolean
				bUseMTStart=true;
				// since most of the other grids will also be used, we set the bUseForcAndTStart-flag to true here, too.
				bUseForcAndTStart=true;
			}

		}

		for (int luz=1;luz<=LUtable[0].getCode();luz++)
		{
			LUtable[luz].InitDynamicPhenology(taga,monata,jahra);
		}

		if (bUseForcAndTStart)
		{
			FORC   = new double**[numlayers+1];
			TSTART = new double**[numlayers+1]; 
			CHILL  = new double**[numlayers+1];
			FSTAR  = new double**[numlayers+1];
			forcindex	  = new int[numlayers+1];
			chillindex	  = new int[numlayers+1];
			tstartindex	  = new int[numlayers+1];
			fstarindex	  = new int[numlayers+1];
			mtstart_index = new int[numlayers+1];

			if (bUseMTStart)
			{
				// MTStart means Multiple Start Times, so its a stack of grids (for each layer one stack --> means a lot of memory)
				// iNumMTStartLevels holds the levels in each stack (was calculates earlier when setting bUseMTStart)
				MTSTART = new double***[numlayers+1];
			}
			else
			{
				MTSTART = NULL;
			}
		}
		else
		{
			FORC   = NULL;
			TSTART = NULL;
			CHILL  = NULL;
			FSTAR  = NULL;
			MTSTART = NULL; // MTSTART is a stack of grids 
			forcindex	  = NULL;
			chillindex	  = NULL;
			tstartindex	  = NULL;
			fstarindex	  = NULL;
			mtstart_index = NULL;
		}

		int albedo_index   = ist_vorhanden(_Albedo_);
		if (albedo_index<=0)
		{
			(*WSstream) << "FATAL: Albedo Grid not found for use within Landuse-Table.\n"
			" This grid must be read in as standard or variable grid (there, it can be created automatically)"<<"[end-entry]";
			throw(12);
		}
		else
		{	
			ALBEDO=gridinfo[albedo_index].pointer;
		}


		char *layername;
		layername=new char[200];
		mtstart_index1=0;
		// suchen des ersten freien Grids
		int gridindex=1;
		while (gridinfo[gridindex].name!=NULL)
		{ 
			gridindex++;
			if (gridindex==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (at position "<<gridindex<<")\n"<<"[end-entry]";
				throw(12);
			}
		}

		// reading or creating the grid for rse
		sprintf(layername,_SurfaceEvaporationResistance_);
		if (!ist_vorhanden(layername)) // could also be read in as variable grid
		{
			initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
			gridinfo[gridindex].pointer = matrix(sr,er,1,ncols,nrows);
			RSE = gridinfo[gridindex].pointer;
			rseindex=gridindex;
			gridanzahl++;
			gridindex++;
		}
		else
		{
			rseindex=getindex(layername);
			RSE = gridinfo[rseindex].pointer;
		}

		// for dynamic phenology, only one grid will be needed for averaging temperatures
		// all other grids will be required per layer
		if (bUseForcAndTStart)
		{
			sprintf(layername,_TemperatureAVG_);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				TAVG = gridinfo[gridindex].pointer;/* nun ist TAVG ansprechbar in Klasse defaultini */
				tavgindex=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				tavgindex=getindex(layername);
				TAVG = gridinfo[tavgindex].pointer;
			}

			sprintf(layername,_StartDayOfLeaveFall_);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				STARTFALL = gridinfo[gridindex].pointer;/* nun ist TAVG ansprechbar in Klasse defaultini */
				gridanzahl++;
				gridindex++;
			}
			else
			{
				STARTFALL = gridinfo[getindex(layername)].pointer;
			}

			sprintf(layername,_TemperatureMOVAVG7_);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				TMOVAVG7 = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, 8, nrows);
				TMOVAVG0 = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(gridindex, layername, TMOVAVG0, TMOVAVG7, 8, layername, gridinfo[tavgindex].writegrid, mpistartrow, mpiendrow);
				tmovavg7index=gridindex;
				for (i=mpistartrow;i<=mpiendrow;i++) 
				{
					for (j=1;j<=ncols;j++)
					{
						if (zn[i][j]!=nodata)
						{
							TMOVAVG0[i][j]=0.0;
							for (int lay=0;lay<=7;lay++)
							{
								TMOVAVG7[i][j][lay]=0.0;
							}
						}
						else
						{
							TMOVAVG0[i][j]=nodata;
							for (int lay=0;lay<=7;lay++)
							{
								TMOVAVG7[i][j][lay]=nodata;
							}
						}
					}
				}
				gridanzahl++;
				gridindex++;
			}
			else
			{
				tmovavg7index=getindex(layername);
				TMOVAVG7 = gridinfo[tmovavg7index].stackpointer;
				TMOVAVG0 = gridinfo[tmovavg7index].pointer;
			}
		}
		else
		{
			TAVG=NULL;
			TMOVAVG7=NULL;
			STARTFALL=NULL;
		}

		for (int layer=1;layer<=numlayers;layer++)
		{

			if (ARCVEGFRACTIONS)
			{
				sprintf(layername, "%s%d", "ArcVeg_Fraction_", layer);
				ARCVEGFRACTIONS[layer]= gridinfo[getindex(layername)].pointer;
				arcvegfractionsindex[layer] = getindex(layername);
			}

			if (AVBIOM)
			{
				sprintf(layername, "%s%d", "ArcVeg_Biomass_", layer);
				AVBIOM[layer] = gridinfo[getindex(layername)].pointer;
				arcvegbiomassindex[layer] = getindex(layername);
			}

			if (AVRETRANS)
			{
				sprintf(layername, "%s%d", "ArcVeg_Retrans_", layer);
				AVRETRANS[layer] = gridinfo[getindex(layername)].pointer;
				arcvegretransindex[layer] = getindex(layername);
			}

			sprintf(layername,"%s%d",_Wurzeltiefe_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				WURZ[layer] = gridinfo[gridindex].pointer;/* nun ist WURZ ansprechbar in Klasse defaultini */
				wurzindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				WURZ[layer]=gridinfo[getindex(layername)].pointer;
				wurzindex[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_Vegetationsbedeckungsgrad_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				VBG[layer] = gridinfo[gridindex].pointer;                            /* nun ist VBG ansprechbar in Klasse defaultini */
				vbgindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				VBG[layer]=gridinfo[getindex(layername)].pointer;
				vbgindex[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_Blattflaechenindex_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				LAI[layer] = gridinfo[gridindex].pointer;                            /* nun ist LAI ansprechbar in Klasse defaultini */
				laiindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				LAI[layer]=gridinfo[getindex(layername)].pointer;
				laiindex[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_SurfaceCanopyResistance_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				RSC[layer] = gridinfo[gridindex].pointer;                            /* nun ist RSC ansprechbar in Klasse defaultini */
				rscindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				RSC[layer]=gridinfo[getindex(layername)].pointer;
				rscindex[layer]=getindex(layername);
			}


			sprintf(layername,"%s%d",_SurfaceInterzepResistance_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				RSI[layer] = gridinfo[gridindex].pointer;                            /* nun ist RSC ansprechbar in Klasse defaultini */
				rsiindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				RSI[layer]=gridinfo[getindex(layername)].pointer;
				rsiindex[layer]=getindex(layername);
			}


			sprintf(layername,"%s%d",_RoughnessLength_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				Z0[layer] = gridinfo[gridindex].pointer;                            /* nun ist Z0 ansprechbar in Klasse defaultini */
				z0index[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				Z0[layer]=gridinfo[getindex(layername)].pointer;
				z0index[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_VegetationStart_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				VEGSTART[layer] = gridinfo[gridindex].pointer;                     /* nun ist VEGSTART ansprechbar in Klasse defaultini */
				vegstartindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				VEGSTART[layer]=gridinfo[getindex(layername)].pointer;
				vegstartindex[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_VegetationStop_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				VEGSTOP[layer] = gridinfo[gridindex].pointer;                     /* nun ist VEGSTOP ansprechbar in Klasse defaultini */
				vegstopindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				VEGSTOP[layer]=gridinfo[getindex(layername)].pointer;
				vegstopindex[layer]=getindex(layername);
			}

			sprintf(layername,"%s%d",_VegetationDuration_,layer);
			if (!ist_vorhanden(layername)) // could also be read in as variable grid
			{
				initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
				gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
				VEGDURAT[layer] = gridinfo[gridindex].pointer;                     /* nun ist VEGDURAT ansprechbar in Klasse defaultini */
				vegduratindex[layer]=gridindex;
				gridanzahl++;
				gridindex++;
			}
			else
			{
				VEGDURAT[layer]=gridinfo[getindex(layername)].pointer;
				vegduratindex[layer]=getindex(layername);
			}

			if (bUseForcAndTStart)
			{
				// grids used for dynamic phenology
				std::stringstream ss2;
				ss2 << _SumOfForcingUnits_ << layer;
				ss2 >> layername;
				if (!ist_vorhanden(layername))
				{
					// creation of the new grid
					// registering the pointer as FORC[i]
					initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
					gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
					FORC[layer] = gridinfo[gridindex].pointer;                            /* nun ist FORC ansprechbar in Klasse defaultini */
					forcindex[layer]=gridindex;
					gridanzahl++;
					gridindex++;
				}
				else
				{	
					FORC[layer]=gridinfo[getindex(layername)].pointer;
					forcindex[layer]=getindex(layername);
				}

				// TStart grids used for dynamic phenology
				std::stringstream ss3;
				ss3 << _Pheno_start_ << layer;
				ss3 >> layername;
				if (!ist_vorhanden(layername))
				{
					// creation of the new grid
					// registering the pointer as TSTART[i]
					initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
					gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
					TSTART[layer] = gridinfo[gridindex].pointer;                            /* nun ist TSTART ansprechbar in Klasse defaultini */
					tstartindex[layer]=gridindex;
					gridanzahl++;
					gridindex++;
				}
				else
				{	
					TSTART[layer]=gridinfo[getindex(layername)].pointer;
					tstartindex[layer]=getindex(layername);
				}

				// CHILL grids used for dynamic phenology in sequential model (used in method 4 only, but created anyway)
				std::stringstream ss4;
				ss4 << _SumOfChillingUnits_ << layer;
				ss4 >> layername;
				if (!ist_vorhanden(layername))
				{
					// creation of the new grid
					// registering the pointer as CHILL[i]
					initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
					gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
					CHILL[layer] = gridinfo[gridindex].pointer;                            /* nun ist CHILL ansprechbar in Klasse defaultini */
					chillindex[layer]=gridindex;
					gridanzahl++;
					gridindex++;
				}
				else
				{	
					CHILL[layer]=gridinfo[getindex(layername)].pointer;
					chillindex[layer]=getindex(layername);
				}

				// FSTAR grids used for dynamic phenology in sequential model (used in method 4 only, but created anyway)
				std::stringstream ss5;
				ss5 << _FStar_ForcingThreshold_ << layer;
				ss5 >> layername;
				if (!ist_vorhanden(layername))
				{
					// creation of the new grid
					// registering the pointer as FSTAR[i]
					initialize_grid(gridindex, layername, NULL, NULL, 0, layername, 0, sr, er);
					gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
					FSTAR[layer] = gridinfo[gridindex].pointer;                            /* nun ist FSTAR ansprechbar in Klasse defaultini */
					fstarindex[layer]=gridindex;
					gridanzahl++;
					gridindex++;
				}
				else
				{	
					FSTAR[layer]=gridinfo[getindex(layername)].pointer;
					fstarindex[layer]=getindex(layername);
				}

				if (bUseMTStart)
				{	
					int starttagnr=0; // get julian day for the start date
					for (int mc=1;mc<monata;mc++)
					{
						starttagnr+=mon[mc];
					}
					starttagnr+=taga;

					// reserving memory for a stack of TStart-Grids, the MTSTART-Stack
					// iNumMTStartLevels holds the levels in each stack
					std::stringstream ss6;
					ss6 << _MultipleStartTimes_ << layer;
					ss6 >> layername;
					mtstart_index[layer]   = ist_vorhanden(layername);
					if (!ist_vorhanden(layername))
					{
						MTSTART[layer]= drei_d_matrix(1,nrows,1,ncols,iNumMTStartLevels+1);
						initialize_grid(gridindex,layername,NULL,MTSTART[layer],iNumMTStartLevels+1,layername,3,1,nrows);
						gridinfo[gridindex].pointer = matrix(sr, er, 1, ncols, nrows);
						mtstart_index[layer]=gridindex;
						if (layer==1) {
							mtstart_index1=mtstart_index[layer];
						}
						gridanzahl++;
						gridindex++;
						// the stack must be initialized using the sample days from the land use table
						// when calculating the dynamic phenology, values for LAI, VBG, Wurz a.s.o. must be interpolated between
						// the sample days using the actual forcing units compared to the forcing thresholds for the last and the next sample day

						for (int i=sr;i<=er;i++) 
						{
							for (int j = 1; j <= ncols; j++)
							{
								use = (int)LUSE[layer][i][j];
								if ((zn[i][j]!=nodata)&&(use!=nodata)&&(mpi_rank==MPIREG[i][j]))
								{
									MTSTART[layer][i][j][0]=0; // later we may want to store special tags here, so level 0 is reserved
									for (int lv=1;lv<=iNumMTStartLevels;lv++)
									{
										if (lv<=LUtable[use].getNumFixPoints())
										{
											// hier wird MTSTART mit -1 für alle dynamischen Punkte belegt, ansonsten mit dem gültigen Julian Day
											if ((LUtable[use].GetForcingThreshold(lv)<0)||(starttagnr >= LUtable[use].GetJulianDay(lv)))
											{
												MTSTART[layer][i][j][lv] = LUtable[use].GetJulianDay(lv)+(int)(LUtable[use].GetAltDep(lv)*(DHM[i][j]-400.0));
											}
											else
											{
												MTSTART[layer][i][j][lv] = -1;
											}
											// if starting day is someweher within the year we have to store the last valid FSTAR-Value in the 0-level of MTSTART
											// in order to be able to interpolate the phenologic parameters between a dynamic or static day and a dynamic day
											if ((starttagnr >= LUtable[use].GetJulianDay(lv))&&(LUtable[use].GetForcingThreshold(lv)>0))
											{
												MTSTART[layer][i][j][0]=LUtable[use].GetForcingThreshold(lv);
											}
										}
										else
										{
											MTSTART[layer][i][j][lv] = nodata;
										}
									}
								}
								else
								{
									for (int lv=1;lv<=iNumMTStartLevels;lv++) 
									{
										MTSTART[layer][i][j][lv]=nodata;
									}
								}
							}
						}
					}
				}
			}

			if (layer==1)
			{
				for (int luz=1;luz<=LUtable[0].getCode();luz++)
				{
					if ((LUtable[luz].GetMethod()==3)||
							(LUtable[luz].GetMethod()==4)||
							(LUtable[luz].GetMethod()==5)||
							(LUtable[luz].GetMethod()==6))
					{
						LUtable[luz].SetTStartGridVector(TSTART);
						LUtable[luz].SetFStarGridVector(FSTAR);
						LUtable[luz].SetFORCGridVector(FORC);
						LUtable[luz].SetMTStartStackVector(MTSTART,iNumMTStartLevels);
					}
				}
			}

			/* fuer jede Landnutzungsform werden nun die aktuellen Parameter berechnet.
		Dabei werden die Knickpunkttage auf die Hoehe umgerechnet.
		Anschliessend wird daraus die nutzbare Feldkapazitaet und der aktuelle
		Inhalt des Bodenspeichers neu berechnet */

			for (int i=sr;i<=er;i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if ((zn[i][j]>nodata) && (mpi_rank==MPIREG[i][j]))
					{
						use = (int)LUSE[layer][i][j];
						if (use>nodata)
						{
							if ((LUtable[use].GetMethod()==3)||(LUtable[use].GetMethod()==4)||(LUtable[use].GetMethod()==6))
							{
								// if the grids were read as variable grids, they may be initialized already. -1 is the sign for not initialized grids, so do it here
								if (FORC[layer][i][j]<0) {
									FORC[layer][i][j]    = 0;
								}
								if (TSTART[layer][i][j]<0) {
									TSTART[layer][i][j]  = LUtable[use].GetTstart();
								}
								if (CHILL[layer][i][j]<0)  {
									CHILL[layer][i][j]   = 0;
								}
								else
								{
									LUtable[use].SetUseMethod3(false); // if a valid chilling units grid was read, set useMethod3 to false, since the sequential model can well be applied in this model year
								}
								if (FSTAR[layer][i][j]<0)  {
									FSTAR[layer][i][j]   = LUtable[use].GetFstar(); // using FStar for Method 3 when initialising the model
								}
							}
							// actual parameters will be calculated for land uses codes with the matching dnyflags only, all other parameters will be read in from the standardgrids/variable grid section
							if (LUtable[use].isFlagSet(ROOTFLAG))
								WURZ[layer][i][j]     = LUtable[use].GetWurz(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (LUtable[use].isFlagSet(LAIFLAG))
								LAI[layer][i][j]      = LUtable[use].GetLAI(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (LUtable[use].isFlagSet(VCFFLAG))
								VBG[layer][i][j]      = LUtable[use].GetVbg(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (LUtable[use].isFlagSet(RSCFLAG))
								RSC[layer][i][j]      = LUtable[use].GetRsc(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (LUtable[use].isFlagSet(RSIFLAG))
								RSI[layer][i][j]      = LUtable[use].GetRsi(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (LUtable[use].isFlagSet(Z0FLAG))
								Z0[layer][i][j]       = LUtable[use].GetZ0(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							if (VEGSTART[layer][i][j]<=0)
							{
								VEGSTART[layer][i][j] = LUtable[use].GetVegStart(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							}
							VEGDURAT[layer][i][j] = 0;
							VEGSTOP[layer][i][j] = -1;

							if (layer==1)
							{
								if (LUtable[use].isFlagSet(RSEFLAG))
									RSE[i][j]    = LUtable[use].GetRse(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
								if (LUtable[use].isFlagSet(ALBEDOFLAG))
									ALBEDO[i][j] = LUtable[use].GetAlbedo(WSstream,taga,monata,jahra,DHM[i][j],layer,i,j);
							}

							/* 1 cm Wurzeltiefe wird als Minimum erachtet (Einfluss der Evaporation)
					in nFK ist immer soviel Wasser, wie 1 cm Wurzeltiefe entspricht */
							if (WURZ[layer][i][j] < double(0.01)) WURZ[layer][i][j] = double(0.01);
							if (VBG[layer][i][j]  < double(0))    VBG[layer][i][j]  = double(0.0); 
							else if (VBG[layer][i][j] > double(1)) VBG[layer][i][j] = double(1.0);
							if (LAI[layer][i][j]  < double(0))    LAI[layer][i][j]  = double(0.0);
							if (RSC[layer][i][j]  < double(0))    RSC[layer][i][j]  = double(0.0);
							if (RSI[layer][i][j]  < double(0))    RSI[layer][i][j]  = double(0.0);
							if (layer==1)
							{
								if (RSE[i][j]  < double(0))    RSE[i][j]  = double(0.0);
							}
							if (Z0[layer][i][j]   < double(0.01))    Z0[layer][i][j]   = double(0.01);

							if (ARCVEGFRACTIONS)
							{
								if (ARCVEGFRACTIONS[layer][i][j] < 0)
								{
									ARCVEGFRACTIONS[layer][i][j] = MultiLayerLUTable[(int)LUSE[0][i][j]].getArcVegLUFractions()[layer];
									gridinfo[arcvegfractionsindex[layer]].writegrid = 8;
								}
							}

							if (AVBIOM)
							{
								if (AVBIOM[layer][i][j] < 0)
								{
									AVBIOM[layer][i][j] = MultiLayerLUTable[(int)LUSE[0][i][j]].getArcVegBiomass()[layer];
									gridinfo[arcvegbiomassindex[layer]].writegrid = 8;
								}
							}

							if (AVRETRANS)
							{
								if (AVRETRANS[layer][i][j] < 0)
								{
									AVRETRANS[layer][i][j] = MultiLayerLUTable[(int)LUSE[0][i][j]].getArcVegRetrans()[layer];
									gridinfo[arcvegretransindex[layer]].writegrid = 8;
								}
							}


						}
						else
						{
							WURZ[layer][i][j] = nodata;
							LAI[layer][i][j]  = nodata;
							VBG[layer][i][j]  = nodata;
							RSC[layer][i][j]  = nodata;
							RSI[layer][i][j]  = nodata;
							if (layer==1)
							{
								RSE[i][j]    = nodata;
								ALBEDO[i][j] = nodata;
							}
							Z0[layer][i][j]       = nodata;
							VEGSTART[layer][i][j] = nodata;
							VEGSTOP[layer][i][j] = nodata;
							VEGDURAT[layer][i][j] = nodata;
							if (bUseForcAndTStart)
							{
								// even for read variable grids we will ensure that cells outside the model domain will set to nodata
								FORC[layer][i][j]   = nodata;
								TSTART[layer][i][j] = nodata;
								CHILL[layer][i][j]  = nodata;
								FSTAR[layer][i][j]  = nodata;
							}
							if (ARCVEGFRACTIONS)
							{
								if (ARCVEGFRACTIONS[layer][i][j] < 0)
								{
									ARCVEGFRACTIONS[layer][i][j] = nodata;
								}
							}
							if (AVBIOM)
							{
								if (AVBIOM[layer][i][j] < 0)
								{
									AVBIOM[layer][i][j] = nodata;
								}
							}
						}
					}
					else 
					{
						WURZ[layer][i][j] = nodata;
						LAI[layer][i][j]  = nodata;
						VBG[layer][i][j]  = nodata;
						RSC[layer][i][j]  = nodata;
						RSI[layer][i][j]  = nodata;
						if (layer==1)
						{
							RSE[i][j]  = nodata;
						}
						Z0[layer][i][j]   = nodata;
						VEGSTART[layer][i][j] = nodata;
						VEGSTOP[layer][i][j] = nodata;
						VEGDURAT[layer][i][j] = nodata;

						if (bUseForcAndTStart)
						{ // even for read variable grids we will ensure that cells outside the model domain will set to nodata
							FORC[layer][i][j]   = nodata;
							TSTART[layer][i][j] = nodata;
							CHILL[layer][i][j]  = nodata;
							FSTAR[layer][i][j]  = nodata;
						}
						if (ARCVEGFRACTIONS)
						{
							if (ARCVEGFRACTIONS[layer][i][j] < 0)
							{
								ARCVEGFRACTIONS[layer][i][j] = nodata;
							}
						}
						if (AVBIOM)
						{
							if (AVBIOM[layer][i][j] < 0)
							{
								AVBIOM[layer][i][j] = nodata;
							}
						}
					}
				}
			}

			MPISynchronizeGridHalo(0, WURZ[layer]);
			MPISynchronizeGridHalo(0, VBG[layer]);
			MPISynchronizeGridHalo(0, LAI[layer]);
			MPISynchronizeGridHalo(0, RSC[layer]);
			MPISynchronizeGridHalo(0, RSI[layer]);
			if (layer == 1)
			{
				MPISynchronizeGridHalo(0, RSE);
			}
			MPISynchronizeGridHalo(0, Z0[layer]);
			 
#ifdef _WASIMDEBUG_
			char *aktname = new char[_stringsize_ + 1];

			sprintf(aktname, "%s_%d", gridinfo[wurzindex[layer]].name, mpi_rank);
			printgridfile(aktname,WURZ[layer],WriteASCIIGrids,WriteDoubleGrids,WriteNetcdfGrids);

			sprintf(aktname, "%s_%d", gridinfo[vbgindex[layer]].name, mpi_rank);
			printgridfile(aktname, VBG[layer], WriteASCIIGrids, WriteDoubleGrids,WriteNetcdfGrids);

			sprintf(aktname, "%s_%d", gridinfo[laiindex[layer]].name, mpi_rank);
			printgridfile(aktname, LAI[layer], WriteASCIIGrids, WriteDoubleGrids,WriteNetcdfGrids);

			sprintf(aktname, "%s_%d", gridinfo[rscindex[layer]].name, mpi_rank);
			printgridfile(aktname, RSC[layer], WriteASCIIGrids, WriteDoubleGrids,WriteNetcdfGrids);

			sprintf(aktname, "%s_%d", gridinfo[rsiindex[layer]].name, mpi_rank);
			printgridfile(aktname, RSI[layer], WriteASCIIGrids, WriteDoubleGrids,WriteNetcdfGrids);
			if (mtstart_index!=NULL) 
			{
				if (mtstart_index[layer]>0)
				{
					sprintf(aktname, "%s_%d_Start", gridinfo[mtstart_index[layer]].name, mpi_rank);
					printstackfile(aktname,MTSTART[layer],iNumMTStartLevels+1,WriteASCIIGrids,WriteDoubleGrids,WriteNetcdfGrids);
				}
			}
			if (layer==1)
			{
				sprintf(aktname, "%s_%d", gridinfo[rseindex].name, mpi_rank);
				printgridfile(aktname, RSE, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}

			sprintf(aktname, "%s_%d", gridinfo[z0index[layer]].name, mpi_rank);
			printgridfile(aktname, Z0[layer], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			delete[] aktname;

			if (ARCVEGFRACTIONS)
			{
				if (gridinfo[arcvegfractionsindex[layer]].writegrid == 8)
				{
					writegridtofile(arcvegfractionsindex[layer]);
				}
				
			}
				
			if (AVBIOM)
			{
				if (gridinfo[arcvegbiomassindex[layer]].writegrid == 8)
				{
					writegridtofile(arcvegbiomassindex[layer]);
				}

			}

			if (AVRETRANS)
			{
				if (gridinfo[arcvegretransindex[layer]].writegrid == 8)
				{
					writegridtofile(arcvegretransindex[layer]);
				}

			}

#endif
		}
		delete[] layername;
	}
	return;
}

void defaultini::check_SOIL_Codes(int l)
{
	int i,j,k,flag;
	if (l==0) return;
	/* Umwandeln der Codes in SOIL in die Nummern 1..SOILtable[0].Code, wie sie in der
	Steuerdatei als Reihenfolge der Bodenarten in der Bodenartentabelle
	auftreten. Die effektive Reihenfolge ist nicht wichtig, da die Parameter
	jeder Bodenart in einem eigenen, durch den Code 1..SOILtable[0].Code zu
	adressierenden Record abgelegt sind (in SOILtable), jedoch duerfen nur die
	Codes 1..SOILtable[0].Code im Bodenartengrid auftreten, da nur soviele
	Records reserviert sind (erster Parameter im Abschnitt "[Bodentabelle]").
	Wird in der Modellierung auf eine bestimmte Bodenart ein bestimmter
	Algorithmus angewandt, so muss anhand des Feldes "Name" des SOILtable-arrays
	entschieden werden, um welche Bodenart es sich handelt. Die vorgegebenen
	Bodenarten in der Steuerdatei duerfen deshalb nicht veraendert
	werden, soweit die Namen im Modell gebraucht werden.
*/
	if (mpi_rank == 0)
	{
		SOIL = gridinfo[l].pointer;
		double **TEMP = matrix(1, nrows, 1, ncols);

		(*WSstream) << "1.: checking soil type codes" << "[end-entry]";
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
				//    for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				/* beim Ersetzen braucht man nachher ein temporaeres Grid, aus dem die Originaldaten gelesen werden koennen.
				Wenn naemlich im Originalgrid z.B. die Codes 8 und 9 vorkommen und in der Bodenartentabelle aber der Code
				8 in der 9. Zeile und der Code 9 in der 8. Zeile vorkommen, dann wird zunaechst die 9 durch die 8 ersetzt
				und anschliessend die 8 durch die 9 - es gibt also keinen Unterschied mehr zwischen 8 und 9, alle achter
				sind nun neuner. Im Tempraergrid bleiben aber die achter achter und die nuener neuner. Deshalb werden auch
				nur die jeweils wirklich benoetigten Codes umgewandelt. */
				TEMP[i][j] = SOIL[i][j];
				if (zn[i][j] > nodata)
				{
					flag = 0;
					for (k = 1; k <= SOILtable[0].getCode(); k++)
					{
						if ((int)SOIL[i][j] == SOILtable[k].getCode()) {
							flag = 1;
							break;
						}
					}
					if (!flag)
					{
						(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
							<< "soil type " << SOIL[i][j] << " in grid " << gridinfo[l].name
							<< " (" << gridinfo[l].identification << ")\n"
							<< "is not present in the control file \n"
							<< "-----------------------------------------------------------------------\n\n" << "[end-entry]";
						throw(32);
					}
				}
			}
		}
		(*WSstream) << "2.:replacing soil type codes " << gridinfo[l].name << "[end-entry]";
		for (k = 1; k <= SOILtable[0].getCode(); k++)
		{
			/* Jede Bodenart wird abgefragt */
			if (SOILtable[k].getCode() != k)

				/* Wenn die Zeilennummer (Index im SOILtable-array) nicht mit
				dem Code fuer diese Bodenart in SOIL uebereinstimmt
				wird im SOIL-Grid Code durch k ersetzt */
			{
				(*WSstream) << SOILtable[k].getCode() << " replaced by " << k << " ..." << "[end-entry]";
				for (i = 1; i <= nrows; i++)
				{
					for (j = 1; j <= ncols; j++)
					{
						//        for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						if (TEMP[i][j] == SOILtable[k].getCode())
						{
							SOIL[i][j] = k;
						}
					}
				}
			}
		}
		free_matrix(TEMP, 1, nrows, 1, ncols);
		TEMP = NULL;
	}
	MPIScatterGrid(l);
	MPISynchronizeGridHalo(l);
}

int defaultini::getUSDASoilTexture(double clay, double silt, double sand)
{
	const int CsiltyClay = 1;
	const int Cclay = 2;
	const int CsiltyClayLoam = 3;
	const int CclayLoam = 4;
	const int Csilt = 5;
	const int CsiltLoam = 6;
	const int CsandyClay= 7;
	const int CLoam = 8;
	const int CsandyClayLoam = 9;
	const int CSandyLoam = 10;
	const int CLoamySand=11;
	const int CSand=12;

	// implementing the USDA soil tetxure triangle

	// convert fractions into percentages
	clay *= 100.0; 
	sand *= 100.0;
	silt *= 100.0;

	// uppermost part with clay above 40
	if (clay >= 40.0) {
		if (sand >= 45.0) {
			return CsandyClay; // 7
		}
		if (silt >= 40.0) {
			return CsiltyClay; // 1
		}
		return Cclay; // 2
	}
	// special case sandyClay above 35% clay
	else if (clay >= 35.0 && sand >= 45.0) {
		return CsandyClay; // 7
	}
	// medium clay range (27% tp 40% clay)
	else if (clay >= 27.0 && clay < 40.0) {
		if (sand >= 45.0)
		{
			return CsandyClayLoam; // 9
		}
		if (sand < 20.0)
		{
			return CsiltyClayLoam; // 3
		}
		return CclayLoam; // 4
	}
	// low clay range (20% to 27% clay)
	else if (clay >= 20.0 && clay < 27.0) {
		if (silt >= 50.0)
		{
			return CsiltLoam; // 6
		}
		if (silt < 28.0) {
			return CsandyClayLoam; // 9
		}
		return CLoam; // 8
	}
	// minimum clay range (clay < 20.0)
	else {
		if (clay <= 12.0 && silt > 80.0 ) {
			return Csilt; // 5
		}
		if (silt >= 50.0)
		{
			return CsiltLoam; // 6
		}
		if (clay >= 7.0 && sand <= 52.0) {
			return CLoam; // 8
		}

		if (sand >= 85.0) {
			if (silt + 1.5 * clay < 15.0) {
				return CSand; // 12
			}
			return CLoamySand; // 11
		}
		else if (sand >= 70.0) {
			if (silt + 2.0 * clay <= 30.0) {
				return CLoamySand; // 11
			}
			return CSandyLoam; // 10
		}
		return CSandyLoam; // 10
	}



	return CsiltyClay;
}

void defaultini::makeSoilGridFromClaySiltSand(int l)
{

	// for testing the usda texture triangle method
	// int res;
	// res = getUSDASoilTexture(0.52, 0.06, 0.42); // expected 2
	// res = getUSDASoilTexture(0.52, 0.01, 0.47); // expected 7
	// res = getUSDASoilTexture(0.45, 0.50, 0.05); // expected 1
	// res = getUSDASoilTexture(0.40, 0.10, 0.50); // expected 7
	// res = getUSDASoilTexture(0.37, 0.30, 0.33); // expected 4
	// res = getUSDASoilTexture(0.37, 0.55, 0.08); // expected 3
	// res = getUSDASoilTexture(0.30, 0.10, 0.60); // expected 9
	// res = getUSDASoilTexture(0.30, 0.30, 0.40); // expected 4
	// res = getUSDASoilTexture(0.30, 0.60, 0.10); // expected 3
	// res = getUSDASoilTexture(0.25, 0.20, 0.55); // expected 9
	// res = getUSDASoilTexture(0.25, 0.40, 0.35); // expected 8
	// res = getUSDASoilTexture(0.25, 0.60, 0.15); // expected 6
	// res = getUSDASoilTexture(0.12, 0.02, 0.86); // expected 11
	// res = getUSDASoilTexture(0.10, 0.08, 0.82); // expected 11
	// res = getUSDASoilTexture(0.12, 0.10, 0.78); // expected 10
	// res = getUSDASoilTexture(0.12, 0.40, 0.48); // expected 8
	// res = getUSDASoilTexture(0.12, 0.60, 0.28); // expected 6
	// res = getUSDASoilTexture(0.13, 0.84, 0.03); // expected 6
	// res = getUSDASoilTexture(0.11, 0.84, 0.05); // expected 5
	// res = getUSDASoilTexture(0.05, 0.05, 0.90); // expected 12
	// res = getUSDASoilTexture(0.00, 0.12, 0.88); // expected 12
	// res = getUSDASoilTexture(0.05, 0.15, 0.80); // expected 11
	// res = getUSDASoilTexture(0.05, 0.25, 0.70); // expected 10
	// res = getUSDASoilTexture(0.05, 0.45, 0.50); // expected 10
	// res = getUSDASoilTexture(0.05, 0.55, 0.40); // expected 6
	// res = getUSDASoilTexture(0.05, 0.85, 0.10); // expected 5

	if (l == 0) return;

	if (!(ist_vorhanden(_clayfraction_) && ist_vorhanden(_siltfraction_) && ist_vorhanden(_sandfraction_)))
	{
		(*WSstream) << "WARNING: no clay, sand and silt grids found, no soil types will be created fr the grain size fraction grids..." << "[end-entry]";
		return;
	}

	if (mpi_rank == 0)
	{

		SOIL = gridinfo[l].pointer;

		double** CLAY = gridinfo[getindex(_clayfraction_)].pointer;
		double** SILT = gridinfo[getindex(_siltfraction_)].pointer;
		double** SAND = gridinfo[getindex(_sandfraction_)].pointer;

		double** TEMP = matrix(1, nrows, 1, ncols);

		(*WSstream) << "analysing clay, sand and silt grids..." << "[end-entry]";

		// using default values for shallow, medium and deept soils
		double mindepth    = 0;
		double mediumdepth = 0;
		double maxdepth = 0;
		int levels = SOILtable[1].getLevels();

		// defaul shallow soils
		double* depths = SOILtable[1].getDzui();
		for (int k = 1; k <= levels; k++)
		{
			mindepth += depths[k];
		}

		// default medium soils
		depths = SOILtable[2].getDzui();
		for (int k = 1; k <= levels; k++)
		{
			mediumdepth += depths[k];
		}

		// default thick soils
		depths = SOILtable[3].getDzui();
		for (int k = 1; k <= levels; k++)
		{
			maxdepth += depths[k];
		}

		if ((mindepth * 1.2 >= mediumdepth) || (mediumdepth * 1.2 >= maxdepth))
		{
			(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
				<< "Error: Depths of soil table entries must be at least differ by factor 1.2 between shallow, medium and thick soils. \n"
				<< "the first three entries of the soil table must be a shallow, a medium and a thick soil, in this order. \n"
				<< "The actual soil table doesn't fullfill this requirement:\n"
				<< "Shallow soil (entry 1, "<< SOILtable[1].getName() <<") is " << mindepth << "m deep\n"
				<< "Medium soil (entry 2, "<< SOILtable[2].getName() <<") is " << mediumdepth << "m deep\n"
				<< "Thick soil (entry 3, "<< SOILtable[3].getName() <<") is " << maxdepth << "m deep\n"
				<< "-----------------------------------------------------------------------\n\n" << "[end-entry]";
			throw(1601);
		}

		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if (zn[i][j] > nodata)
				{
					bool isvalidcell = false;
					if ((CLAY[i][j] > nodata) && (SILT[i][j] > nodata) && (SAND[i][j] > nodata))
					{
						// we have a valid cell 
						isvalidcell = true;
					}
					else
					{
						// todo:
						// get nearest neighbors value here, if the  cell is not rock or water

					}
					if (isvalidcell)
					{
						double clay = CLAY[i][j];
						double silt = SILT[i][j];
						double sand = SAND[i][j];
						double sum = clay + sand + silt;
						int actsoiltype = (int)SOIL[i][j];

						// estimating from the first 3 entries of the soil table if the actual soil is shallo, medium, or thick
						double* actdepths = SOILtable[actsoiltype].getDzui();

						double actThickness = 0;

						for (int k = 1; k <= levels; k++)
						{
							actThickness += actdepths[k];
						}

						bool isShallow = false;
						bool isMedium = false;
						bool isThick = false;
						int thicknessOffset = 1;
						if (actThickness < 1.2 * mindepth) {
							isShallow = true;
						}
						else if (actThickness < 1.2 * mediumdepth) {
							isMedium = true;
							thicknessOffset = 1;
						}
						else {
							isThick = true;
							thicknessOffset = 2;
						}


						if (sum <= 0)
						{
							(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
								<< "at row " << i << ", column " << j << " the sum of clay, sand and silt is 0"
								<< "-----------------------------------------------------------------------\n\n" << "[end-entry]";
							throw(1600);
						}
						else
						{
							// normalizing weight of the fractions to 1.0
							clay = clay / sum;
							silt = silt / sum;
							sand = sand / sum;
								
							// getting the Soiltype from the soil triangle (after US
							int usdaSoil = getUSDASoilTexture(clay,silt,sand);
							TEMP[i][j] = usdaSoil + thicknessOffset;
							
						}
					}
					else
					{
						TEMP[i][j] = SOIL[i][j];
					}
				}
				else
				{
					TEMP[i][j] = nodata;
				}
			}
		}

		

		// replacing soil types by the ones, computed before, but only for the non-nodata cells
		(*WSstream) << "replacing default soil codes with the computed ones ..." << "[end-entry]";
		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if (TEMP[i][j] > nodata)
				{
					SOIL[i][j] = TEMP[i][j];
				}
				else
				{
					SOIL[i][j] = nodata;
				}
			}
		}
		free_matrix(TEMP, 1, nrows, 1, ncols);
		TEMP = NULL;
	}
	
	MPIScatterGrid(l);
	MPISynchronizeGridHalo(l);

	// writing the result grid to default output
	if (mpi_rank == 0)
	{
		char *aktname = new char[500];
		sprintf(aktname, "ComputedSoils_%d", mpi_rank);
		printgridfile(aktname, SOIL, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
		delete[] aktname;
	}
}

void defaultini::check_irrig_Codes(int l)
{
	int i,j,k,flag;
	if (l==0) return;
	/* Umwandeln der Codes in IRRIG in die Nummern 1..irrigationtable[0].Code, wie sie in der
	Steuerdatei als Reihenfolge der Landnutzungen in der Landnutzungstabelle
	auftreten. Die effektive Reihenfolge ist nicht wichtig, da die Parameter
	jeder Nutzungsform in einem eigenen, durch den Code 1..irrigationtable[0].Code zu
	adressierenden Record abgelegt sind, jedoch duerfen nur die
	Codes 1..irrigationtable[0].Code im irrigationgrid auftreten, da nur soviele
	Records reserviert sind (erster Parameter im Abschnitt irrigation_table).
	Wird in der Modellierung auf eine bestimmte Bewaesserungsart ein bestimmter
	Algorithmus angewandt, so muss anhand des Feldes "Name" des irrigationtable-arrays
	entschieden werden, um welche Landnutzung es sich handelt. Die vorgegebenen
	Namen in der Steuerdatei duerfen deshalb nicht veraendert
	werden, soweit die Namen im Modell gebraucht werden. */
	if (mpi_rank == 0)
	{
		IRRIG = gridinfo[l].pointer;
		double **TEMP = matrix(1, nrows, 1, ncols);

		(*WSstream) << "1.: checking codes in irrigation grid" << "[end-entry]";
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
				//    for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				TEMP[i][j] = IRRIG[i][j]; /* Erklaerung analog zur Erklaerung in der naechsten Routine check_SOIL_Codes() */
				if ((zn[i][j] > nodata) && (IRRIG[i][j] > 0))
				{
					flag = 0;
					for (k = 1; k <= irrigationtable[0].Code; k++)
					{
						if ((int)IRRIG[i][j] == irrigationtable[k].Code) flag = 1;
					}
					if (!flag)
					{
						(*WSstream) << "FATAL: \n-----------------------------------------------------------------------\n"
							"irrigation type " << IRRIG[i][j] << " in grid " << gridinfo[l].name
							<< " (" << gridinfo[l].identification << ")\n"
							"is not present in the irrigation table in the control file \n"
							"-----------------------------------------------------------------------\n\n" << "[end-entry]";
						throw(31);
					}
				}
			}
		}

		(*WSstream) << "2.:replacing codes in irrigation grid " << gridinfo[l].name << "[end-entry]";
		for (k = 1; k <= irrigationtable[0].Code; k++)
			/* Jede Landnutzungsform wird abgefragt */
			if (irrigationtable[k].Code != k)
				/* Wenn die Zeilennummer (Index im irrigationtable-array) nicht mit
				dem Code fuer diese Nutzung in IRRIG uebereinstimmt
				wird im IRRIG-Grid Code durch k ersetzt */
			{
				(*WSstream) << irrigationtable[k].Code << " replaced by " << k << " ..." << "[end-entry]";
				for (i = 1; i <= nrows; i++)
				{
					for (j = 1; j <= ncols; j++)
					{
						//        for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						if (TEMP[i][j] == irrigationtable[k].Code)
						{
							IRRIG[i][j] = k;
						}
					}
				}
			}
		free_matrix(TEMP, 1, nrows, 1, ncols);
		TEMP = NULL;
	}
	// distribute corrected irrigation codes to all other instances
	MPIScatterGrid(l);
	MPISynchronizeGridHalo(l);
}


void defaultini::readPeriodicStandardGrids(long minute,long stunde,long tag,long monat,long jahr,long Interval)
{

	//  int    iPer; // periodicity of input grid (to be read each time the period is over
	//  int   iPerAct; // actual counter of elapsed periods as given in cPerType (H, D, M or Y)
	//  int    iSubInt; // number of subintervals of the periodtype which must be over until reading a new grid (e.g. 15 for 15 days when cPerType is M
	//  int   iSubIntAct; // actual counter of elapsed sub intervals until the grid must be imported again
	//  int   iLastSubInt; // used to determine if the actual sub interval is less than the old one (when a new month starts, the day will be smaller than the last day of the old month)
	//  int   iLastPer; // same for period counter, if usefull 
	//  char   *cPerType; // type of periodicity: H, D, M, Y
	//  char   *cPerGridName; // gridname to read from (will be constructed each time a grid is read from the actual dates)
	

	for (int i=3;i<=gridanzahl+2;i++)
	{
		bool readnow=false;
		if (gridinfo[i].iPer>0)
		{
			std::string pertype(gridinfo[i].cPerType);
			if (pertype.find("H")!=std::string::npos) {
				if (stunde>gridinfo[i].iLastSubInt) // if multiple hours are used to identify e.g. 15-minutes intervals, only real hour-steps are regarded here (and only the first of the subintervals)
				{
					gridinfo[i].iSubIntAct++;
				}
				else if (stunde<gridinfo[i].iLastSubInt) // a new day starts
				{
					gridinfo[i].iPerAct++;
					gridinfo[i].iSubIntAct=1;
				}
				gridinfo[i].iLastSubInt=stunde;
			} else if (pertype.find("D")!=std::string::npos) {
				if (stunde>gridinfo[i].iLastSubInt) // if multiple hours are used to identify e.g. 15-minutes intervals, only real hour-steps are regarded here (and only the first of the subintervals)
				{
					gridinfo[i].iSubIntAct++;
				}
				else if (stunde<gridinfo[i].iLastSubInt) // a new day starts
				{
					gridinfo[i].iPerAct++;
					gridinfo[i].iSubIntAct=1;
				}
				gridinfo[i].iLastSubInt=stunde;
			} else if (pertype.find("M")!=std::string::npos) {
				if (tag>gridinfo[i].iLastSubInt) // new day has started
				{
					gridinfo[i].iSubIntAct++;
				}
				else if (tag<gridinfo[i].iLastSubInt) // a new month starts
				{
					gridinfo[i].iPerAct++; // month is the item of periodicity
					gridinfo[i].iSubIntAct=1;
				}
				gridinfo[i].iLastSubInt=tag;
			} else if (pertype.find("Y")!=std::string::npos) {
				if (monat>gridinfo[i].iLastSubInt) // new day has started
				{
					gridinfo[i].iSubIntAct++;
				}
				else if (monat<gridinfo[i].iLastSubInt) // a new month starts
				{
					gridinfo[i].iPerAct++; // month is the item of periodicity
					gridinfo[i].iSubIntAct=1;
				}
				gridinfo[i].iLastSubInt=monat;
			}
			if ((gridinfo[i].iPerAct==gridinfo[i].iPer)&&(gridinfo[i].iSubInt==gridinfo[i].iSubIntAct))
			{
				readnow=true;
				// now creating the filename and reading the grid
				sprintf(gridinfo[i].cPerGridName,"%s_%04ld%02ld%02ld%02ld",gridinfo[i].name,jahr,monat,tag,stunde);
				readgridfromfile(gridinfo[i].cPerGridName, &(gridinfo[i]).pointer, gridinfo[i].fillcode, _AllocateGridMemoryOff_, _ShowMessagesOn_, gridinfo[i].defaultValue, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
				// special treatment for landuse grid since some other grids depend on that grid (the layers!)
				if (std::string(gridinfo[i].name).compare("landuse")==0)
				{
					// set all teh layer grids with the same algorithms as in generateMultiLayerLanduse()
					int cnt=MultiLayerLUTable[1].getNumLayers();
//					LUSE[0] = gridinfo[ist_vorhanden(_Landnutzung_)].pointer;
					if (LUSE[0]!=NULL)
					{
						char *layername;
						layername=new char[200];
						for (int z=1;z<=cnt;z++)
						{
//							sprintf(layername,"%s%d",_Landnutzung_,z);
//							LUSE[z] = gridinfo[ist_vorhanden(layername)].pointer;; 
							for (int i=mpistartrow;i<=mpiendrow;i++)
							{
								for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
								//			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
								{
									if ((zn[i][j]>nodata) && (zn[i][j]==MPIREG[i][j]))
									{
										LUSE[z][i][j]=MultiLayerLUTable[(int)LUSE[0][i][j]].getLULayer(z);
									}
									else LUSE[z][i][j] = nodata;
								}
							}
						}
						delete[] layername;
					}
				}
				MPISynchronizeGridHalo(i);
				// resetting iPerAct and iSubIntAct 
				gridinfo[i].iSubIntAct=0;
				gridinfo[i].iPerAct=0; 
				// when the newly read in grid should be valid for the next interval only
				if (gridinfo[i].bPersistency == false)
				{
					gridinfo[i].iReplaceCnt=0;
					gridinfo[i].bReplaceGridNow=true;
				}
			} else if (gridinfo[i].bReplaceGridNow == true) {
				// replace default grid by its original read in at model start
				gridinfo[i].iReplaceCnt++;
				if (gridinfo[i].iReplaceCnt==gridinfo[i].iReplaceAfter)
				{
					if (gridinfo[i].TempGrid!=NULL)
					{
						for (int r = mpistartrow; r <= mpiendrow; r++)
						{
							for (int c=FirstLastCol[r].First;c<=FirstLastCol[r].Last;c++)
							{
								gridinfo[i].pointer[r][c]=gridinfo[i].TempGrid[r][c];
							}
						}
					}
					MPISynchronizeGridHalo(i);
					gridinfo[i].bReplaceGridNow=false;
				}
			}
		}
	}

}

void defaultini::read_variable_Grids()
{
	int aktgrid;
	int i,j,k,wrgr,aktfillcode;
	double defaultvalue;
	char *aktname,*aktident;
	if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}
	findeabschnitt(_variable_Grids_);
	fscanfln(steuerdatei,"%d",&varGridsCount," number of variable grids ",tstr);

	if (varGridsCount > 0)
	{
		/* zunaechst suchen des ersten freien EIntrages */
		i=1;
		while (gridinfo[i].name!=NULL)
		{
			i++;
			if (i==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (an Position "<<i<<")\n"<<"[end-entry]";
                                fclose(steuerdatei);
				throw(12);
			}
		}
		varGridsStart=i;

		varGridsRead  = (int *)calloc(varGridsCount+1,sizeof(int));
		aktname=newstring(_stringsize_,1);
		aktident=newstring(_stringsize_,1);

		for (j = 1;j <= varGridsCount; j++)
		{
			aktgrid=j+2+gridanzahl;
			/* Gridnamen und Bezeichnung lesen */
			if (fscanf(steuerdatei,"%s%s%d%lf",aktname,aktident,&aktfillcode,&defaultvalue)!=4)
			{
				sprintf(tstr,"\nfile for vargrid %d or the identifier or the fillcode or the defaultvalue (which is new in version 7.7) in ",j);
				steuererror(tstr,instanzname);
			}
			fgets(tstr,_stringsize_,steuerdatei);
			fscanfln(steuerdatei,"%d",&wrgr,"Writegrid for a variable grid ",tstr);
			fscanfln(steuerdatei,"%d",&(varGridsRead[j]),"Readgrid for a variable grid ",tstr);

			int sr = (mpistartrow>1) ? mpistartrow - 1 : mpistartrow;
			int er = (mpiendrow<nrows) ? mpiendrow + 1 : mpiendrow;
			if (mpi_rank == 0)
			{
				sr = 1;
				er = nrows;
			}

			initialize_grid(aktgrid, aktname, gridinfo[aktgrid].pointer, NULL, 0, aktident, wrgr,sr,er);
			if (varGridsRead[j] == 1)
			{
				readgridfromfile(gridinfo[aktgrid].name, &(gridinfo[aktgrid]).pointer, aktfillcode, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
			}
			else
			{
				gridinfo[aktgrid].pointer=matrix(sr,er,1,ncols,nrows);
				for (i = sr; i <= er; i++)
				{
					for (k = 1; k <= ncols; k++)
					{
						gridinfo[aktgrid].pointer[i][k] = zn[i][k] > nodata ? defaultvalue : nodata;
					}
				}
			}

			/* Bodenspeicher als relative Bodenfeuchte entweder einlesen oder auf 1 festlegen  */
			if (strcmp(gridinfo[aktgrid].identification,_Bodenspeicher_)==0)
			{
				SB = gridinfo[aktgrid].pointer; /* nun kann auf SB zugegriffen werden in der Klasse defaultini  */
			}
		}

		if (aktname) free(aktname); aktname=NULL;
		if (aktident) free(aktident); aktident=NULL;
	}
	/* da gridanzahl an das Hauptprogramm uebergeben wird und die maximal bereits belegte Anzahl von Grids
	anzeigt, wird nun varGridsCount zu gridanzahl addiert. */
        fclose(steuerdatei);
	gridanzahl += varGridsCount;
}

void defaultini::make_nFK_kWert_and_Na()
{
	int i,j,k;
	char *aktname,*aktident;

	nfkindex=naindex=kindex=0;

	/* Wenn ein Bodenartengrid eingelesen wurde, wird nun geprueft,
	ob nFK und SB schon als Standardgrids eingelesen wurden.
	Wenn nicht, so werden beide Grids nun aus den Daten der Bodenartentabelle
	erzeugt. Wird weder "Bodenarten" noch "nutzbare_Feldkapazitaet" noch "Bodenspeicher"
	im gridinfo-array gefunden, so gibt es keinen Fehler, da ja auch eine
	Modellierung ohne das Bodenmodell moeglich sein kann, es wird erst dann ein
	Fehler ausgegeben, wenn ein Teilmodell, welches eines der Grids benoetigt, dieses nicht
	finden kann */
	soilindex =ist_vorhanden(_Bodenarten_);

	if (soilindex)
	{
		aktname=newstring(_stringsize_,1);
		aktident=newstring(_stringsize_,1);

		/*    SOIL=gridinfo[soilindex].pointer;  */              /* von nun an ist SOIL ansprechbar, wurde schon in check_SOIL_codes belegt */
		k=ist_vorhanden(_Nutzbare_Feldkapazitaet_);
		if (k==0)
		{ /* ist nFK noch nicht eingelesen, wird sie nun aus Landnutzung und Bodenartenparametern erzeugt */
			nfkindex=1;
			while (gridinfo[nfkindex].name!=NULL)
			{ nfkindex++;
				if (nfkindex==gridmaxanzahl)
				{ 
					(*WSstream) << "FATAL: no free grid found in memory (at position "<<nfkindex<<"[end-entry]";
					throw(12);
				}
			}

			if (strrchr(gridinfo[soilindex].name,'.'))
			{ i=(int)(strrchr(gridinfo[soilindex].name,'.')-gridinfo[soilindex].name);
				sprintf(aktname,"%.*s.intern.fk",i,gridinfo[soilindex].name);
			} else  sprintf(aktname,"%s.intern.fk",gridinfo[soilindex].name);
			strcpy(aktident,_Nutzbare_Feldkapazitaet_);
			initialize_grid(nfkindex,aktname,NULL,NULL,0,aktident,0,1,nrows);
			gridinfo[nfkindex].pointer = matrix(1,nrows,1,ncols);
			nFK=gridinfo[nfkindex].pointer;
			for (i=1;i<=nrows;i++)
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if (zn[i][j]>nodata)
				{
					/* nFK in SOILtable ist in Vol.% angegeben, Umrechnung in mm: *10.0
			Beruecksichtigung der Wurzeltiefe durch LUtable.Wurz oder 0.4m,
			wenn keine Landnutzung eingelesen wurde*/
					nFK[i][j] = WURZ[1][i][j] * SOILtable[(int)SOIL[i][j]].getNFK()*10.0;
				}
				else nFK[i][j] = nodata;
			}

			/* die genauere Berechnung der nFK erfolgt jeweils anhand der aktuellen Wurzeltiefe in defaultini::int_tvar_params() */
			gridanzahl++;
			gridinfo[nfkindex].writegrid=1;
			writegridtofile(nfkindex);
			gridinfo[nfkindex].writegrid=0;
		}
		else
		{
			nfkindex=k;
			nFK=gridinfo[nfkindex].pointer;                            /* nun ist nFK ansprechbar in Klasse defaultini */
		}
		k=ist_vorhanden(_Hydraulische_Leitfaehigkeit_);
		if (k==0)
		{ /* max. Bodenspeicher (auffuellbare Porositaet nun erzeugen,
		da es noch nicht existiert */
			kindex=1;
			while (gridinfo[kindex].name!=NULL)
			{ kindex++;
				if (kindex==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no free grid found in memeory (at position "<<kindex<<")"<<"[end-entry]";
					throw(12);
				}
			}

			if (strrchr(gridinfo[soilindex].name,'.'))
			{ i=(int)(strrchr(gridinfo[soilindex].name,'.')-gridinfo[soilindex].name);
				sprintf(aktname,"%.*s.intern.k",i,gridinfo[soilindex].name);
			} else  sprintf(aktname,"%s.intern.k",gridinfo[soilindex].name);
			strcpy(aktident,_Hydraulische_Leitfaehigkeit_);
			initialize_grid(kindex, aktname, NULL, NULL, 0, aktident, 0, 1, nrows);
			gridinfo[kindex].pointer = matrix(1,nrows,1,ncols);
			kWert=gridinfo[kindex].pointer;                        /* von hier an ist kWert in der Klasse defaultini ansprechbar */
			for (i=1;i<=nrows;i++)
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			kWert[i][j] = zn[i][j] != nodata ? SOILtable[(int)SOIL[i][j]].getK() : nodata;
			gridanzahl++;
			gridinfo[kindex].writegrid=1;
			writegridtofile(kindex);
			gridinfo[kindex].writegrid=0;
		}
		else
		{
			kindex=k;
			kWert=gridinfo[kindex].pointer;                        /* von hier an ist kWert in der Klasse defaultini ansprechbar */
		}
		k=ist_vorhanden(_Auffuellbare_Porositaet_);
		if (k==0)
		{ /* max. Bodenspeicher (auffuellbare Porositaet nun erzeugen,
		da es noch nicht existiert */
			naindex=1;
			while (gridinfo[naindex].name!=NULL)
			{ naindex++;
				if (naindex==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no free grid found in memeory (at position "<<naindex<<")"<<"[end-entry]";
					throw(12);
				}
			}
			if (strrchr(gridinfo[soilindex].name,'.'))
			{ i=(int)(strrchr(gridinfo[soilindex].name,'.')-gridinfo[soilindex].name);
				sprintf(aktname,"%.*s.intern.na",i,gridinfo[soilindex].name);
			} else  sprintf(aktname,"%s.intern.na",gridinfo[soilindex].name);
			strcpy(aktident,_Auffuellbare_Porositaet_);
			initialize_grid(naindex, aktname, NULL, NULL, 0, aktident, 0, 1, nrows);
			gridinfo[naindex].pointer = matrix(1,nrows,1,ncols);
			Na=gridinfo[naindex].pointer;                          /* ab jetzt wird Na in der Klasse defaultini ansprechbar */
			for (i=1;i<=nrows;i++)
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			/* auffuellbare Porositaet ist gleich maximaler Bodenspeicher-Welkepunkt
				Umrechnung von Vol.% in relativen Anteil: *0.01 bzw. /100.0*/
			Na[i][j] = zn[i][j] != nodata ? SOILtable[(int)SOIL[i][j]].getMSB()*0.01 : nodata;
			gridanzahl++;
			gridinfo[naindex].writegrid=1;
			writegridtofile(naindex);
			gridinfo[naindex].writegrid=0;
		}
		else
		{
			naindex=k;
			Na=gridinfo[naindex].pointer;                          /* ab jetzt wird Na in der Klasse defaultini ansprechbar */
		}
		if (aktname) free(aktname); aktname=NULL;
		if (aktident) free(aktident); aktident=NULL;
	}
}

/* fuer Zellen, die mindestens einen nodata-Wert als Rand haben kann keine
3x3-Matrix zur Berechnugn von slope, aspect, curvature aufgestellt werden.
Diese Kennwerte werden mit der folgenden Funktion aus den Werten der
direkten Nachbarn mit gueltigen Werten berechnet */
void defaultini::getnachbar(int i,int j,double **Grid)
{
	double  summe=0.0;
	double anzahl=0.0;
	if ((i>1)    &&(j<ncols)) if (Grid[i-1][j+1]!=nodata) {summe+=Grid[i-1][j+1];anzahl++;}
	if             (j<ncols)  if (Grid[i  ][j+1]!=nodata) {summe+=Grid[i  ][j+1];anzahl++;}
	if ((i<nrows)&&(j<ncols)) if (Grid[i+1][j+1]!=nodata) {summe+=Grid[i+1][j+1];anzahl++;}
	if  (i<nrows)             if (Grid[i+1][j  ]!=nodata) {summe+=Grid[i+1][j  ];anzahl++;}
	if ((i<nrows)&&(j>1))     if (Grid[i+1][j-1]!=nodata) {summe+=Grid[i+1][j-1];anzahl++;}
	if                 (j>1)  if (Grid[i  ][j-1]!=nodata) {summe+=Grid[i  ][j-1];anzahl++;}
	if ((i>1)    &&(j>1))     if (Grid[i-1][j-1]!=nodata) {summe+=Grid[i-1][j-1];anzahl++;}
	if  (i>1)                 if (Grid[i-1][j  ]!=nodata) {summe+=Grid[i-1][j  ];anzahl++;}
	if (anzahl>0)
	{
		Grid[i][j]=summe/anzahl;
	}
}

void defaultini::make_Slope_and_Aspect()
{
	int i,j,asp,slp;
	double H,G;
	char *aktname,*aktident;

	aspectindex=slopeindex=0;

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	if (ist_vorhanden(_Hoehenmodell_))
	{
		aktname=newstring(_stringsize_,1);
		aktident=newstring(_stringsize_,1);

		dhmindex=ist_vorhanden(_Hoehenmodell_); /* eigentlich redundant, DHM ist immer in Position 2 des gridinfo-arrays */
		/*    DHM=gridinfo[dhmindex].pointer; */
		asp=ist_vorhanden(_Exposition_);
		slp=ist_vorhanden(_Gefaelle_);
		if ((asp==0)||(slp==0))
		{ /* ist Exposition oder Gefaelle noch nicht eingelesen, wird sie nun erzeugt */
			slopeindex=1;
			if (slp==0)
			{
				while (gridinfo[slopeindex].name!=NULL)
				{ slopeindex++;
					if (slopeindex==gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no free grid found in memeory (at position "<<slopeindex<<")"<<"[end-entry]";
						throw(12);
					}
				}
				if (strrchr(gridinfo[dhmindex].name,'.'))
				{ 
					i=(int)(strrchr(gridinfo[dhmindex].name,'.')-gridinfo[dhmindex].name);
					sprintf(aktname,"%.*s.intern.sl",i,gridinfo[dhmindex].name);
				} 
				else
				{
					sprintf(aktname, "%s.intern.sl", gridinfo[dhmindex].name);
				}
				strcpy(aktident,_Gefaelle_);
				initialize_grid(slopeindex, aktname, NULL, NULL, 0, aktident, 0, sr, er);
				gridinfo[slopeindex].pointer = matrix(sr,er,1,ncols,nrows);
				SLP=gridinfo[slopeindex].pointer;                             /* SLP nun ansprechbar in der Klasse defaultini */
				gridanzahl++;
			}
			else
			{
				slopeindex=slp;
				SLP = gridinfo[slp].pointer;
			}
			aspectindex=1;
			if (asp==0)
			{
				while (gridinfo[aspectindex].name!=NULL)
				{ aspectindex++;
					if (aspectindex==gridmaxanzahl)
					{
						(*WSstream) << "FATAL: no free grid left (at position "<<aspectindex<<")\n"<<"[end-entry]";
						throw(12);
					}
				}
				if (strrchr(gridinfo[dhmindex].name,'.'))
				{
					i=(int)(strrchr(gridinfo[dhmindex].name,'.')-gridinfo[dhmindex].name);
					sprintf(aktname,"%.*s.intern.ex",i,gridinfo[dhmindex].name);
				} 
				else
				{
					sprintf(aktname, "%s.intern.ex", gridinfo[dhmindex].name);
				}
				strcpy(aktident,_Exposition_);
				initialize_grid(aspectindex,aktname,NULL,NULL,0,aktident,0,sr,er);
				gridinfo[aspectindex].pointer = matrix(sr,er,1,ncols,nrows);
				ASP=gridinfo[aspectindex].pointer;                            /* ASP nun ansprechbar in der Klasse defaultini */
				gridanzahl++;
			}
			else
			{
				aspectindex=asp;
				ASP = gridinfo[asp].pointer;
			}

			(*WSstream) << "calculating the lokal topography (aspect, slope ...).. "<<"[end-entry]";
			for (j=1;j<=ncols;j++)
			{
				if (!slp) 
				{
					if (sr<=1) SLP[1][j] = nodata; 
					if (er>=nrows) SLP[nrows][j] = nodata;
				}
				if (!asp) 
				{ 
					if (sr <= 1) ASP[1][j] = nodata;
					if (er >= nrows) ASP[nrows][j] = nodata;
				}
			}

			for (i=sr+1;i<er;i++)
			{
				if (!slp) { SLP[i][1] = nodata; SLP[i][ncols] = nodata;}
				if (!asp) { ASP[i][1] = nodata; ASP[i][ncols] = nodata;}
				for (j=2;j<ncols;j++)
				{
					/* es wird eine Flaeche 2. Ordnung durch die 9 Punkte gelegt. Fuer die Berechnung des Gefaelles
			bedeutet das, dass es sich (als 1. Ableitung) als geometrisches Mittel der Ableitungen
			nach x und nach y ergibt. Folglich wird die Flaeche 2. Ordung in 2 lineare Gleichungen 1.Ordnung
			zerlegt, was dazu fuehrt, dass nur die 4 direkt an die Zelle grenzenden 4 Punkte
			zur Berechnung benoetigt werden (fuer Aspect gilt dasselbe)  */
					if ( (DHM[i-1][j  ]!=nodata)&&
							(DHM[i  ][j-1]!=nodata)&&
							(DHM[i  ][j+1]!=nodata)&&
							(DHM[i+1][j  ]!=nodata))
					{
						/* Ermitteln der Koefizienten des Gleichungssystems */
						H = (DHM[i][j+1]-DHM[i][j-1])/(cellsize+cellsize);
						G = (DHM[i-1][j]-DHM[i+1][j])/(cellsize+cellsize);
						if (!slp) { if ((G!=0)||(H!=0)) SLP[i][j] = 180.0/Pi*atan(sqrt(G*G + H*H));            else SLP[i][j]=0;}
						if (!asp) { if ( H!=0)          ASP[i][j] = 180.0-180.0/Pi*atan(G/H)+90.0*(H/fabs(H)); else ASP[i][j] = G>0 ? 180.0 : 0;}
					}
					else
					{
						if (!asp) ASP[i][j]=nodata;
						if (!slp) SLP[i][j]=nodata;
					}
					if (!slp) { SLP[i][1] = nodata; SLP[i][ncols] = nodata;}
					if (!asp) { ASP[i][1] = nodata; ASP[i][ncols] = nodata;}
				}
			}
			/* Fehlwerte dort auffuellen, wo im Hoehenmodell Daten stehen -> naechster Nachbar oder Mittelwert der Nachbarn */
			for (i=sr;i<=er;i++)
			{
				for (j=1;j<=ncols;j++)
				{
					if ((!asp)&&(DHM[i][j]!=nodata)&&(ASP[i][j]==nodata))  
					getnachbar(i,j,ASP);
					if ((!slp)&&(DHM[i][j]!=nodata)&&(SLP[i][j]==nodata))  
					getnachbar(i,j,SLP);
				}
			}
			// second iteration 
			for (i=sr;i<=er;i++)
			{
				for (j=1;j<=ncols;j++)
				{
					if ((!asp)&&(DHM[i][j]!=nodata)&&(ASP[i][j]==nodata))  
					getnachbar(i,j,ASP);
					if ((!slp)&&(DHM[i][j]!=nodata)&&(SLP[i][j]==nodata))  
					getnachbar(i,j,SLP);
				}
			}
			if (!asp)
			{
				MPISynchronizeGridHalo(aspectindex, ASP);
				gridinfo[aspectindex].writegrid=1;
				writegridtofile(aspectindex);
				gridinfo[aspectindex].writegrid=0;
			}
			if (!slp)
			{
				MPISynchronizeGridHalo(slopeindex, SLP);
				gridinfo[slopeindex].writegrid = 1;
				writegridtofile(slopeindex);
				gridinfo[slopeindex].writegrid=0;
			}
		}
		else
		{
			aspectindex = asp;
			ASP = gridinfo[asp].pointer;
			slopeindex = slp;
			SLP = gridinfo[slp].pointer;
		}
		if (aktname) free(aktname); aktname=NULL;
		if (aktident) free(aktident); aktident=NULL;
	}
}

void defaultini::check_ice_firn_grid()
{
	int *list;
	int i,j,iceindex;
	iceindex=1;
	while (strcmp(gridinfo[iceindex].identification,_ice_firn_)!=0)
	{
		iceindex++;
		if ((iceindex>=gridmaxanzahl)||(gridinfo[iceindex].name==NULL))
		{
			(*WSstream) << "INFO: ice-grid not found in memory (at position "<<iceindex<<")"<<"[end-entry]";
			iceindex=-1;
			break;
		}
	}

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	if (iceindex>0)
	{
		ICE=gridinfo[iceindex].pointer;
		list=(int*)calloc(zonmax+1,sizeof(int));
		/* part 1: building a list with subbasin codes which are glaciers */
		for (i=sr;i<=er;i++)
		{
			for (j=1;j<=ncols;j++)
			{
				if (zn[i][j]>nodata)
				{
					if (ICE[i][j]>0) /* allowed values: 1 for ice, 2 for firn */
					{
						list[zarray[(int)zn[i][j]]]=1;
					}
				}
			}
		}
		/* part 2: checking subasins, if somewhere are cells which are not ice or firn, if at least one cell of the subbasin is ice or firn */
		/* the new glacier model will no longer require 1 to 1 matching for a glacier and it's catchment
	for (i=1;i<=nrows;i++)
	{
	for (j=1;j<=ncols;j++)
	{
		if (zn[i][j]>nodata)
		{
		if (list[zarray[(int)zn[i][j]]]&&((ICE[i][j]<1)||(ICE[i][j]>2)))
		{
			(*WSstream) << "FATAL: error in subbasin codes with respect to the codes in the ice_firn_grid:\n"
				"glaciers (which have the codes 1 for ice and 2 for firn) have to be defined \n"
				"as an own subbasin (each glacier may ba a single subbasin or several \n"
				"glaciers may form a subbasin, this is free). It is not allowed that a \n"
				"glacier-subbasin (defined within the zone-grid!) has codes in the \n"
				"ice_firn_grid which are different from 1 (for ice) or 2 (for firn).\n"
				"It is also not allowed, that ice or firn in the glacier grid are somewhere\n"
				"outside of a glacier-subbasin (e.g., the matching cells in the zone grid and\n"
				"in the ice_firn_grid have to coincide).\n"<<"[end-entry]";
			throw(69);
		}
		}
	}
	}
	*/
		if (list) free(list); list=NULL;
	}
}

/* sucht im angegebenen Grid (index i im Gridinfo-array) nach dem
Maximalwert und uebergibt ihn als Rueckgabeparameter */
double defaultini::get_h_maximum(int i)
{
	double locHmax = nodata;
	double globalhmax=nodata;
	for (int j = mpistartrow; j <= mpiendrow; j++)
	{
		for (int k = FirstLastCol[j].First; k <= FirstLastCol[j].Last; k++)
		{
			if (locHmax < gridinfo[i].pointer[j][k])
			{
				locHmax = gridinfo[i].pointer[j][k];
			}
		}
	}
#ifdef _MPI // synchronizing external keyboard events with all instances
//cout << "hmax = " << locHmax << " in rank " << mpi_rank << endl;
	MPI_Reduce(&locHmax, &globalhmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&globalhmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	locHmax = globalhmax;
#endif
	return locHmax;
}


void defaultini::analyze_regions(int regioindex)
{
	int i,j;
	int        zonenid;
	statrecord tempzone;
	unsigned short int   weiter;
	int        *zzarray;
	statrecord *zzone;
	int maxregionenwert;

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	if (mpi_rank == 0)
	{
		(*WSstream) << "ranking regions from " << gridinfo[regioindex].name << " for " << instanzname << "[end-entry]";
		zzone = (statrecord *)calloc(maxzonenzahl + 1, sizeof(statrecord));
		if (!zzone) speichererror(1);
		zzarray = (Tzarray *)calloc((size_t)maxzonenwert + 1, sizeof(Tzarray));
		if (!zzarray) speichererror(1);
		zzone[0].n = 0.0;
		zzone[0].sx = 0.0;
		regiomax = 0;
		regiomaxwert = 0;
		maxregionenwert = maxzonenwert; /* Regionen-codes koennen bis 32768 reichen */
		maxregionenzahl = _max_regionen_zahl_; /* aber insgesamt sollen nichtmehr als z.Zt.20 Regionen verwendet werden */
		for (i = 1; i < maxregionenwert; i++)   zzarray[i] = 0;
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
			{
				if (zn[i][j] != nodata)
				{
					zonenid = (int)(REGION1[i][j] + 0.5);
					if ((zonenid<1) || (zonenid>maxregionenwert))
					{
						(*WSstream) << "FATAL: invalid region code in " << gridinfo[regioindex].name
							<< "!\n (codes has to be between 1 and " << maxregionenwert << ")!\n" << "[end-entry]";
						throw(3);
					}
					if (zzarray[zonenid] == 0)
					{ /*gibt es die Zone mit der ID zonenid noch nicht? */
						regiomax = regiomax + 1;           /* nein, also in die Liste einfuegen */
						if (regiomax > maxregionenzahl)
						{
							(*WSstream) << "FATAL: too much zones in " << gridinfo[regioindex].name
								<< ", (max. " << maxregionenzahl << ")!\n" << "[end-entry]";
							throw(4);
						}
						if (zonenid > regiomaxwert) regiomaxwert = zonenid;
						zzarray[zonenid] = regiomax;
						zzone[zzarray[zonenid]].zonenwert = REGION1[i][j];
						zzone[zzarray[zonenid]].n = 0;
						zzone[zzarray[zonenid]].sx = 0;
					}
				}
			}
		}
		/* Ordnen der Zonen nach dem Zonenwert */
		for (i = 1; i <= regiomax - 1; i++)
		{
			weiter = true;
			for (j = 2; j <= regiomax; j++)
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
		rarray = (Tzarray *)calloc(regiomaxwert + 1, sizeof(Tzarray));
		if (!rarray) speichererror(1);
		for (i = 0; i <= regiomaxwert; i++) rarray[i] = zzarray[i];

		if (zzarray) free(zzarray); zzarray = NULL;
		if (zzone) free(zzone); zzone = NULL;
		/* getting the centers of mass for all regions in grid coordinates
		(rows/colums) and real coordinates */
		center_of_mass = (center_of_mass_type*)calloc(regiomax + 1, sizeof(center_of_mass_type));
		for (i = 1; i <= regiomax; i++)
		{
			/* initializing with 0 */
			center_of_mass[i].row = 0;
			center_of_mass[i].col = 0;
			center_of_mass[i].ycoord = 0;
			center_of_mass[i].xcoord = 0;
			center_of_mass[i].n = 0;
		}
		/* adding up the values of row and column in the yccord and xcoord fields */
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
			{
				//	  if ((REGION1[i][j]!=nodata))
				if ((zn[i][j] != nodata) && (REGION1[i][j] != nodata))
				{
					center_of_mass[rarray[(int)REGION1[i][j]]].row += i;
					center_of_mass[rarray[(int)REGION1[i][j]]].col += j;
					center_of_mass[rarray[(int)REGION1[i][j]]].ycoord += i;
					center_of_mass[rarray[(int)REGION1[i][j]]].xcoord += j;
					center_of_mass[rarray[(int)REGION1[i][j]]].n++;
				}
			}
		}
		/* getting the average row and column number of each region */
		for (i = 1; i <= regiomax; i++)
		{
			center_of_mass[i].row = center_of_mass[i].row / center_of_mass[i].n;   /* in n the number of cells for each region is stored */
			center_of_mass[i].col = center_of_mass[i].col / center_of_mass[i].n;
			center_of_mass[i].ycoord = center_of_mass[i].ycoord / center_of_mass[i].n;
			center_of_mass[i].xcoord = center_of_mass[i].xcoord / center_of_mass[i].n;
		}
	}
	// now bcasting regiomaxwert and regiomax, 
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&regiomax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&regiomaxwert, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// creating rarray in all other ranks and bcast rarray from rank 0
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank > 0)
	{
		rarray = (int *)calloc(regiomaxwert + 1, sizeof(int));
		if (!rarray) speichererror(1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Datatype vectortype;
	MPI_Type_contiguous(regiomaxwert + 1,MPI_INT,&vectortype);
	MPI_Type_commit(&vectortype);
	MPI_Bcast(&(rarray[0]), 1, vectortype, 0,MPI_COMM_WORLD);
	MPI_Type_free(&vectortype);
	
	// bcastig center_of_mass to all other ranks
	MPI_Barrier(MPI_COMM_WORLD);
	const int nitems = 5;
	int blocklengths[5] = { 1, 1, 1, 1, 1 };
	MPI_Datatype types[5] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	MPI_Aint offsets[5] = { offsetof(center_of_mass_type, row), 
							 offsetof(center_of_mass_type, col), 
							 offsetof(center_of_mass_type, ycoord), 
							 offsetof(center_of_mass_type, xcoord), 
							 offsetof(center_of_mass_type, n) };
	MPI_Datatype MPICenterOfMassType;
	MPI_Type_create_struct(nitems, blocklengths, offsets, types,&MPICenterOfMassType);
	MPI_Type_commit(&MPICenterOfMassType);

	if (mpi_rank > 0)
	{
		center_of_mass = (center_of_mass_type*)calloc(regiomax + 1, sizeof(center_of_mass_type));
		if (!center_of_mass) speichererror(1);
	}
	MPI_Bcast(&(center_of_mass[0]), regiomax + 1, MPICenterOfMassType, 0,MPI_COMM_WORLD);

	MPI_Type_free(&MPICenterOfMassType);

#endif

}

#pragma optimize ("p",on)
void defaultini::SetRegionWeights()
{
	int i,j,ii,jj,m,n,o,p,r;
	double radius,r2,dist,inv2r,wmin,
	x2,x3,y2,y3,d1,d2,d3,dsum,dmin;

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	if (mpi_rank == 0)
	{

		(*WSstream) << "generating region weight grids" << "[end-entry]";
		if (transdist <= 0) transdist = 0.01;
		radius = transdist / cellsize; /* search radius in cells */
		r2 = radius*radius;
		inv2r = 1.0 / (2.0*radius);
		// first, make sure, no distant regions will affect the actual model area
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
			{
				if (zn[i][j] <= nodata)
				{
					REGION1[i][j] = nodata;
				}
			}
		}
		for (i = 1; i <= nrows; i++)
		{
			for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if (zn[i][j] > nodata)
				{
					/* loking around to find other regions near the actual cell */
					WEIGHT1[i][j] = radius;
					if (regiomax > 1)
					{
						WEIGHT2[i][j] = radius;
						REGION2[i][j] = 0;
						if (regiomax > 2)
						{
							WEIGHT3[i][j] = radius;
							REGION3[i][j] = 0;
						}
						/* suche in cellweiten von der aktuellen Zeile nach oben (bis n) und unten (bis m) */
						n = (i - (int)radius) >= 1 ? (i - (int)radius) : 1;
						m = (i + (int)radius) <= nrows ? (i + (int)radius) : nrows;
						for (ii = n; ii <= m; ii++)
						{
							/* die Suche in x richtung kann nun nach der Kreisformel begrenzt werden */
							r = (int)sqrt(r2 - (i - ii)*(i - ii));
							o = (j - r) >= 1 ? (j - r) : 1;
							p = (j + r) < ncols ? (j + r) : ncols;
							/* nur wenn im Regionengrid an der Stelle ii,jj ein gueltiger Wert ist,
				soll eine Neugewichtung erfolgen */
							for (jj = o; jj <= p; jj++)
							{
								/* wenn an der Suchstelle eine andere Region gefunden wurde */
								if ((REGION1[ii][jj] != REGION1[i][j]) && (REGION1[ii][jj] > nodata))
								{
									if (REGION2[i][j] == 0) REGION2[i][j] = REGION1[ii][jj]; /* a new region is found */
									dist = (int)sqrt(double((i - ii)*(i - ii) + (j - jj)*(j - jj)));
									//                  weight=(dist/(2.0*radius)+0.5); /* minimum weight of 0.5 for the own region is assumed? */
									if (REGION1[ii][jj] == REGION2[i][j])
									{
										if (dist <= WEIGHT1[i][j])
										{
											WEIGHT1[i][j] = dist; /* the same region is found but nearer to the actual cell */
											WEIGHT2[i][j] = dist;
											x2 = jj;
											y2 = ii;
										}
									}
									else
									{    /* another region was found than given in REGION2[i][j] */
										if (REGION3[i][j] == 0) REGION3[i][j] = REGION1[ii][jj];
										if (REGION1[ii][jj] == REGION3[i][j])
										{
											if (dist <= WEIGHT3[i][j])
											{
												WEIGHT3[i][j] = dist; /* the same region is found but nearer to the actual cell */
												x3 = jj;
												y3 = ii;
											}
										}
										else
										{ /* although at least 2 other regions are already known, there hasbeen found another one,
						so if it is nearer than the others it should replace the until now farest region */
											if (WEIGHT2[i][j] < WEIGHT3[i][j]) /* weight3: region 3 is farer away,s check against region3 now */
											{
												if (WEIGHT3[i][j] > dist)
												{
													REGION3[i][j] = REGION1[ii][jj];
													WEIGHT3[i][j] = dist;
													x3 = jj; y3 = ii;
												}
											}
											else
											{ /* if region 2 is farer away than region 3, check against weight 2 */
												if (WEIGHT2[i][j] > dist)
												{
													REGION2[i][j] = REGION1[ii][jj];
													WEIGHT2[i][j] = dist;
													x2 = jj; y2 = ii;
												}
											}
										}
									}
								}
							}
						}
					}
					/* now converting distances into relative weights */
					if (regiomax > 2)
					{
						WEIGHT1[i][j] = (WEIGHT1[i][j] / (2.0*radius) + 0.5);
						if (WEIGHT1[i][j] > 1.0) WEIGHT1[i][j] = 1.0; else if (WEIGHT1[i][j] < 0.5) WEIGHT1[i][j] = 0.5;
						if (REGION3[i][j] > 0)
						{ /* the actual cell is affected by 3 regions */
							/* distance to each region */
							d2 = sqrt((i - y2)*(i - y2) + (j - x2)*(j - x2));
							d3 = sqrt((i - y3)*(i - y3) + (j - x3)*(j - x3));
							/* for the weights, the distance to the inner area (where weight is 100%) is needed */
							/* for the first region, this distance is the maximum radius minus the distance to the border */
							if (d2 < d3)
								d1 = radius - d2;
							else
								d1 = radius - d3;
							/* for the other ergions, this distance is given by the distance to the border of this
				region plus the maximum search radius */
							d2 += radius;
							d3 += radius;
							/* no distances of 0 are allowed (because of inverting the distances) */
							if (d1 <= 0) d1 = 1;
							if (d2 <= 0) d2 = 1;
							if (d3 <= 0) d3 = 1;
							/* the minimum distance to any region is used for setting the correct weights */
							if (d1 < d2) dmin = d1; else dmin = d2;
							if (dmin > d3) dmin = d3;

							wmin = cos(1.5707963*(dmin / radius));
							wmin = 0.9 - 0.5*(wmin);

							//            wmin=0.4+0.5*(dmin/radius)*(dmin/radius); 
							/* dmin may vary between 1 and radius. If dmin is small, the weight of the
				inverted distance should be reduced by a larger offset than for farer minimum distances.
				the limits are scaling factor 1 for small distances and 0.5 for far distances */
							d1 = 1.0 / d1 - wmin*inv2r;
							d2 = 1.0 / d2 - wmin*inv2r;
							d3 = 1.0 / d3 - wmin*inv2r;
							if (d1 < 0) d1 = 0;
							if (d2 < 0) d2 = 0;
							if (d3<0) d3 = 0;
							dsum = d1 + d2 + d3;
							WEIGHT1[i][j] = d1 / dsum;
							WEIGHT2[i][j] = d2 / dsum;
							if ((WEIGHT1[i][j] + WEIGHT2[i][j])>1.0)
								WEIGHT2[i][j] = 1.0 - WEIGHT1[i][j];
							WEIGHT3[i][j] = 1.0 - (WEIGHT1[i][j] + WEIGHT2[i][j]);
						}
						else
						{
							WEIGHT3[i][j] = 0.0;
							if (REGION2[i][j]>0)
							{ /* although there are at least 3 regions, the actual cell is only affected by 2 regions,
					so a linear interpolation is enough -> weight 2 is 1-weight 1 because of weight 1 is
					already calculated */
								WEIGHT2[i][j] = 1.0 - WEIGHT1[i][j];
							}
							else
							{ /* if the actual cell is only affected by one region, the weights of regions 2 and 3 are 0 */
								WEIGHT2[i][j] = 0.0;
							}
						}
					}
					else
					{
						if (regiomax > 1)
						{
							WEIGHT1[i][j] = (WEIGHT1[i][j] / (2.0*radius) + 0.5);
							if (WEIGHT1[i][j] > 1.0) WEIGHT1[i][j] = 1.0; else if (WEIGHT1[i][j] < 0.5) WEIGHT1[i][j] = 0.5;
							WEIGHT2[i][j] = 1.0 - WEIGHT1[i][j];
						}
						else
							WEIGHT1[i][j] = 1.0;
					}
				}
				else
				{
					REGION1[i][j] = nodata;
					WEIGHT1[i][j] = nodata;
					if (regiomax > 1)
					{
						REGION2[i][j] = nodata;
						WEIGHT2[i][j] = nodata;
						if (regiomax > 2)
						{
							REGION3[i][j] = nodata;
							WEIGHT3[i][j] = nodata;
						}
					}
				}
			}
			(*WSstream) << ".";
		}
		(*WSstream) << "[end-entry]";
	}
}
#pragma optimize ("p",off)


void  defaultini::MakeRegionWeightGrids()
{
	int i;
	char *aktname;

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	aktname=newstring(_stringsize_,1);
	REGION1=REGION2=REGION3=WEIGHT1=WEIGHT2=WEIGHT3=NULL;
	regiomax=0;regiomaxwert=0;rarray=NULL;center_of_mass=NULL;
	regio1index=ist_vorhanden(_regression_regions_1_);

	// new since 8.05 (10090912): if no region grid was read in, it is created now (with exactly one region with code 1)
	// Thus, the superposition code is much easier to implement
	if (regio1index<=0)
	{
		regio1index=1;
		while (gridinfo[regio1index].name!=NULL)
		{ 
			regio1index++;
			if (weight1index==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (at position "<<regio1index<<")\n"<<"[end-entry]";
				throw(12);
			}
		}
		if (strrchr(gridinfo[1].name,'.'))
		{ 
			i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
			sprintf(aktname,"%.*s.intern.r1",i,gridinfo[1].name);
		} 
		else
		{
			sprintf(aktname, "%s.intern.r1", gridinfo[1].name);
		}
		initialize_grid(regio1index, aktname, NULL, NULL, 0, _regression_regions_1_, 0, sr, er);
		gridinfo[regio1index].pointer = matrix(sr,er,1,ncols,nrows);
		REGION1=gridinfo[regio1index].pointer;
		for (int i=sr;i<=er;i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				if (zn[i][j] != nodata)
				{
					REGION1[i][j] = 1;
				}
				else
				{
					REGION1[i][j] = nodata;
				}
			}
		}
		gridanzahl++;
	}

	if (regio1index>0) /* ansonsten gibt es nur eine Region und kein regio-grid wurde eingelesen */
	{
		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		findeabschnitt(_region_transition_distance_);
		fscanfln(steuerdatei,"%lf",&transdist,"transition distance (if or if not using a region grid) ",tstr);
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		REGION1=gridinfo[regio1index].pointer;

		analyze_regions(regio1index);
		
		weight1index=1;
		while (gridinfo[weight1index].name!=NULL)
		{ 
			weight1index++;
			if (weight1index==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (at position "<<weight1index<<")"<<"[end-entry]";
				throw(12);
			}
		}
		if (strrchr(gridinfo[regio1index].name,'.'))
		{ 
			i=(int)(strrchr(gridinfo[regio1index].name,'.')-gridinfo[regio1index].name);
			sprintf(aktname,"%.*s.intern.w1",i,gridinfo[regio1index].name);
		} 
		else
		{
			sprintf(aktname, "%s.intern.rw1", gridinfo[regio1index].name);
		}
		initialize_grid(weight1index, aktname, NULL, NULL, 0, _regression_weights_1_, 0, sr, er);
		gridinfo[weight1index].pointer = matrix(sr, er, 1, ncols, nrows);
		WEIGHT1 = gridinfo[weight1index].pointer;
		gridanzahl++;

		/* 2 regions should be thereif a region grid is read in. However, if not, then the following
	section is skipped */
		if (regiomax>1)
		{
			/* if there are 2 regions, then we need 2 new grids: a weight grid and a region grid,
		the latter stores the second region influencing the actual cell and the weight grid
		stores the weight of that second region (depends on the distance of the actual cell
		to that region) */
			regio2index=weight1index+1;
			if (strrchr(gridinfo[regio1index].name,'.'))
			{ 
				i=(int)(strrchr(gridinfo[regio1index].name,'.')-gridinfo[regio1index].name);
				sprintf(aktname,"%.*s.intern.r2",i,gridinfo[regio1index].name);
			} 
			else
			{
				sprintf(aktname, "%s.intern.r2", gridinfo[regio1index].name);
			}
			initialize_grid(regio2index, aktname, NULL, NULL, 0, _regression_regions_2_, 0, sr, er);
			gridinfo[regio2index].pointer = matrix(sr, er, 1, ncols, nrows);
			REGION2 = gridinfo[regio2index].pointer;
			gridanzahl++;

			weight2index=regio2index+1;
			if (strrchr(gridinfo[regio1index].name,'.'))
			{ 
				i=(int)(strrchr(gridinfo[regio1index].name,'.')-gridinfo[regio1index].name);
				sprintf(aktname,"%.*s.intern.w2",i,gridinfo[regio1index].name);
			} 
			else
			{
				sprintf(aktname, "%s.intern.w2", gridinfo[regio1index].name);
			}
			initialize_grid(weight2index, aktname, NULL, NULL, 0, _regression_weights_2_, 0, sr, er);
			gridinfo[weight2index].pointer = matrix(sr, er, 1, ncols, nrows);
			WEIGHT2 = gridinfo[weight2index].pointer;
			gridanzahl++;

			if (regiomax>2)
			{
				/* if there are three regions, than it may happen, that three or more regions are neighboring each other,
		so to consider at least 3 regions, there are 2 more grids needed, the region code of the 3rd considered region
		and its weight for the actual cell */
				regio3index=weight2index+1;
				if (strrchr(gridinfo[regio1index].name,'.'))
				{ 
					i=(int)(strrchr(gridinfo[regio1index].name,'.')-gridinfo[regio1index].name);
					sprintf(aktname,"%.*s.intern.r3",i,gridinfo[regio1index].name);
				} 
				else
				{
					sprintf(aktname, "%s.intern.r3", gridinfo[regio1index].name);
				}
				initialize_grid(regio3index, aktname, NULL, NULL, 0, _regression_regions_3_, 0, sr, er);
				gridinfo[regio3index].pointer = matrix(sr, er, 1, ncols, nrows);
				REGION3 = gridinfo[regio3index].pointer;
				gridanzahl++;

				weight3index=regio3index+1;
				if (strrchr(gridinfo[regio1index].name,'.'))
				{ 
					i=(int)(strrchr(gridinfo[regio1index].name,'.')-gridinfo[regio1index].name);
					sprintf(aktname,"%.*s.intern.w3",i,gridinfo[regio1index].name);
				} 
				else
				{
					sprintf(aktname, "%s.intern.w3", gridinfo[regio1index].name);
				}
				initialize_grid(weight3index, aktname, NULL, NULL, 0, _regression_weights_3_, 0, sr, er);
				gridinfo[weight3index].pointer = matrix(sr, er, 1, ncols, nrows);
				WEIGHT3 = gridinfo[weight3index].pointer;
				gridanzahl++;

			}
		}
		if (aktname) free(aktname); aktname=NULL;

		/* calculating the weights and region codes for the other regions */

		SetRegionWeights();

		// scatter grids from rank 0 to all instances

		/* rausschreiben der Gewichte, damit man sich dieDinger mal ansehen kann */
		MPIScatterGrid(weight1index, WEIGHT1);
		MPISynchronizeGridHalo(weight1index, WEIGHT1);

		gridinfo[weight1index].writegrid=1;
		writegridtofile(weight1index);
		gridinfo[weight1index].writegrid=0;
		if (regiomax>1)
		{
			gridinfo[regio2index].writegrid=1;
			gridinfo[weight2index].writegrid=1;

			MPIScatterGrid(regio2index, REGION2);
			MPISynchronizeGridHalo(regio2index, REGION2);

			MPIScatterGrid(weight2index, WEIGHT2);
			MPISynchronizeGridHalo(weight2index, WEIGHT2);

			writegridtofile(regio2index);
			writegridtofile(weight2index);
			gridinfo[regio2index].writegrid=0;
			gridinfo[weight2index].writegrid=0;
			if(regiomax>2)
			{
				gridinfo[regio3index].writegrid=1;
				gridinfo[weight3index].writegrid=1;

				MPIScatterGrid(regio2index, REGION3);
				MPISynchronizeGridHalo(regio2index, REGION3);

				MPIScatterGrid(weight2index, WEIGHT3);
				MPISynchronizeGridHalo(weight2index, WEIGHT3);

				writegridtofile(regio3index);
				writegridtofile(weight3index);
				gridinfo[regio3index].writegrid=0;
				gridinfo[weight3index].writegrid=0;
			}
		}
	}
}

void defaultini::CheckRelCellArea()
{
	int i,j,k;
	char *aktname;

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	aktname=newstring(_stringsize_,1);
	/* relcellarea und cellsizex, cellsizey sollen immer berechnet werden, damit spaeter damit die 
	idw-interpolation mit echten Koordinaten berechnet werden kann */
	/*  if ((RelCellArea!=NULL) ||
	((RelCellArea==NULL) && ist_vorhanden(_CellSizeX_) && ist_vorhanden(_CellSizeY_))
	)
*/
	{

		k=1;
		while (gridinfo[k].name!=NULL)
		{ 
			k++;
			if (k==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid found in memory (at position "<<k<<")"<<"[end-entry]";
				throw(12);
			}
		}

		/* is there a cell size grid for cell sizes in x direction? */
		if (ist_vorhanden(_CellSizeX_))
		{
			CellSizeX=gridinfo[getindex(_CellSizeX_)].pointer;
			cellsizexread=true;
			cellsizexindex=k;
		}
		else
		{
			/* if not, then generate the grid (values are added later from RelCellArea) */
			if (strrchr(gridinfo[1].name,'.'))
			{ 
				i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
				sprintf(aktname,"%.*s.intern.csx",i,gridinfo[1].name);
			} 
			else
			{
				sprintf(aktname, "%s.intern.csx", gridinfo[1].name);
			}
			initialize_grid(k, aktname, NULL, NULL, 0, _CellSizeX_, 0, sr, er);
			//gridinfo[k].pointer = matrix(1,nrows,1,ncols);
			gridinfo[k].pointer = matrix(sr,er,1,ncols,nrows);
			CellSizeX = gridinfo[k].pointer;
			gridanzahl++;
			cellsizexread=false;
			cellsizexindex=k;
			k++;
		}
		/* is there a CellSize grid for cell sizes in y direction? */
		if (ist_vorhanden(_CellSizeY_))
		{
			CellSizeY=gridinfo[getindex(_CellSizeY_)].pointer;
			cellsizeyread=true;
			cellsizeyindex=k;
		}
		else
		{
			/* if not, then generate the grid (values are added later from RelCellArea) */
			if (strrchr(gridinfo[1].name,'.'))
			{
				i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
				sprintf(aktname,"%.*s.intern.csy",i,gridinfo[1].name);
			} 
			else
			{
				sprintf(aktname, "%s.intern.csy", gridinfo[1].name);
			}
			initialize_grid(k, aktname, NULL, NULL, 0, _CellSizeY_, 0, sr, er);
			//gridinfo[k].pointer = matrix(1,nrows,1,ncols);
			gridinfo[k].pointer = matrix(sr, er, 1, ncols, nrows); 
			CellSizeY = gridinfo[k].pointer;
			gridanzahl++;
			cellsizeyread=false;
			cellsizeyindex=k;
			k++;
		}

		if (RelCellArea==NULL)
		{
			if (strrchr(gridinfo[1].name,'.'))
			{ 
				i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
				sprintf(aktname,"%.*s.intern.rca",i,gridinfo[1].name);
			} 
			else
			{
				sprintf(aktname, "%s.intern.rca", gridinfo[1].name);
			}
			initialize_grid(k, aktname, NULL, NULL, 0, _RelCellArea_, 0, sr, er);
			//gridinfo[k].pointer = matrix(1,nrows,1,ncols);
			gridinfo[k].pointer = matrix(sr, er, 1, ncols, nrows);
			RelCellArea = gridinfo[k].pointer;
			gridanzahl++;
			relcellarearead=false;
			relcellareaindex=k;
			k++;
			if ((!cellsizeyread)&&(!cellsizexread))
			{ /* if both rel cell size grids and the relcellareagrid 
		are not read in, then relcellarea is set to a default value of 1.0 */
				for (i=sr;i<=er;i++)
				{
					for (j=1;j<=ncols;j++)
					{
						if (zn[i][j]!=nodata)
						{
							RelCellArea[i][j]=1.0;
						}
						else
						{
							RelCellArea[i][j]=nodata;
						}
					}
				}
			}
		}
		else
		{
			relcellarearead=true;
			relcellareaindex=getindex(_RelCellArea_);
		}

		(*WSstream) << "Checking cell sizes for variable Cells..."<<"[end-entry]";
		for (i=sr;i<=er;i++)
		{
			for (j=1;j<=ncols;j++)
			{
				if (zn[i][j]!=nodata)
				{
					if (cellsizexread)
					{
						if (CellSizeX[i][j]==0)
						{
							(*WSstream) << "FATAL: cell size in x direction is zero at row "<<i
							<<", column "<<j<<",\nread from file "<<gridinfo[cellsizexindex].name
							<<"\n"<<"[end-entry]";
							throw(66);
						}
						else
						{
							if (cellsizeyread)
							{
								if (relcellarearead)
								{ /* all three grids are read, so they will be checked on plausibility */
									if (RelCellArea[i][j]-CellSizeX[i][j]*CellSizeY[i][j]/(cellsize*cellsize)>0.00001)
									{
										(*WSstream) << "FATAL: Cell size Y ("<<CellSizeX[i][j]<<", from "
										<< gridinfo[cellsizexindex].name<<") times\n"
										"cell size X ("<<CellSizeY[i][j]<<", from "<<gridinfo[cellsizeyindex].name<<") \n"
										"divided by the mean cell area is "<<CellSizeX[i][j]*CellSizeY[i][j]/(cellsize*cellsize)
										<<" and differs from \n"
										"the relative grid cell size "<<RelCellArea[i][j]
										<<" (from "<<gridinfo[getindex(_RelCellArea_)].name<<") by more than 0.001%%\n"<<"[end-entry]";
									}
									throw(66);
								}
								else
								RelCellArea[i][j]=CellSizeX[i][j]*CellSizeY[i][j]/(cellsize*cellsize);
							}
							else
							/* CellSizeY not read, so RelCellArea and CellSizeX must be there,
				otherwise this statement would never be reached */
							{
								CellSizeY[i][j]=RelCellArea[i][j]*cellsize*cellsize/CellSizeX[i][j];
							}
						}
					}
					else
					{
						if (cellsizeyread)
						{ /* CellSizeX was not read, so this branch was reached */
							if (CellSizeY[i][j]==0)
							{
								(*WSstream) << "FATAL: cell size in y direction is zero at row "<<i<<", column "<<j
								<< ",\nread from file "<<gridinfo[cellsizeyindex].name<<"[end-entry]";
								throw(66);
							}
							else
							{
								CellSizeX[i][j]=RelCellArea[i][j]*cellsize*cellsize/CellSizeY[i][j];
							}
						}
						else
						{  /* wether CellSizeX nor CellSizeY are read, so both values have to be derived from RelCellArea */
							CellSizeX[i][j]=sqrt(cellsize*cellsize*RelCellArea[i][j]);
							CellSizeY[i][j]=CellSizeX[i][j];
						}
					}
				}
				else
				{
					RelCellArea[i][j]=CellSizeX[i][j]=CellSizeY[i][j]=nodata;
				}
			}
		}
	}
	if (aktname) free(aktname); aktname=NULL;
	if (!relcellarearead)
	{
		gridinfo[relcellareaindex].writegrid=1;
		writegridtofile(relcellareaindex);
		gridinfo[relcellareaindex].writegrid=0;
	}
	if (!cellsizexread)
	{
		gridinfo[cellsizexindex].writegrid=1;
		writegridtofile(cellsizexindex);
		gridinfo[cellsizexindex].writegrid=0;
	}
	if (!cellsizeyread)
	{
		gridinfo[cellsizeyindex].writegrid=1;
		writegridtofile(cellsizeyindex);
		gridinfo[cellsizeyindex].writegrid=0;
	}
}

void defaultini::MakeCoordGrids()
{
	int i,j,k;
	char *aktname;
	aktname=newstring(_stringsize_,1);
	k=1;
	while (gridinfo[k].name!=NULL)
	{ k++;
		if (k==gridmaxanzahl)
		{
			(*WSstream) << "FATAL: no free grid found in memory (at position "<<k<<")"<<"[end-entry]";
			throw(12);
		}
	}
	/* generate coordx and coordy-grids */
	if (strrchr(gridinfo[1].name,'.'))
	{ 
		i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
		sprintf(aktname,"%.*s.intern.cx",i,gridinfo[1].name);
	} 
	else
	{
		sprintf(aktname, "%s.intern.cx", gridinfo[1].name);
	}
	initialize_grid(k,aktname,NULL,NULL,0,_CoordX_,0,mpistartrow,mpiendrow);
	gridinfo[k].pointer = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
	CoordX = gridinfo[k].pointer;
	gridanzahl++;
	coordxindex=k;
	k++;
	/* now coordy */
	if (strrchr(gridinfo[1].name,'.'))
	{ 
		i=(int)(strrchr(gridinfo[1].name,'.')-gridinfo[1].name);
		sprintf(aktname,"%.*s.intern.cy",i,gridinfo[1].name);
	} 
	else
	{
		sprintf(aktname, "%s.intern.cy", gridinfo[1].name);
	}
	initialize_grid(k,aktname,NULL,NULL,0,_CoordY_,0,mpistartrow,mpiendrow);
	gridinfo[k].pointer = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
	CoordY = gridinfo[k].pointer;
	gridanzahl++;
	coordyindex=k;
	k++;

	/* storing coordinates of the cells into the coord-grids */
	for (i=mpiendrow;i>=mpistartrow;i--)
	{
		for (j=1;j<=ncols;j++)
		{
			CoordX[i][j]=(j-1)*cellsize+xll;
			CoordY[i][j]=((nrows-i)*cellsize)+yll;
		}
	}
	if (RelCellArea)
	{
		for (j=FirstLastCol[mpiendrow].First+1;j<=ncols;j++)
		{
			CoordX[mpiendrow][j] = CoordX[mpiendrow][j - 1] + (CellSizeX[mpiendrow][j - 1]>nodata ? CellSizeX[mpiendrow][j - 1] : cellsize); /* lower side is oriented from west to east */
		}
		for (i=mpiendrow-1;i>=mpistartrow;i--)
		{
			CoordX[i][1]=xll; /* the left side is oriented from south to north */
			for (j=2;j<=ncols;j++)
			{
				CoordX[i][j]=CoordX[i][j-1]+(CellSizeX[i][j-1]>nodata?CellSizeX[i][j-1]:cellsize);
				CoordY[i][j]=CoordY[i+1][j]+(CellSizeY[i+1][j]>nodata?CellSizeY[i+1][j]:cellsize);
			}
		}
	}
	for (i=mpistartrow;i<=mpiendrow;i++)
	{
		for (j=1;j<=ncols;j++)
		{
			if (zn[i][j]==nodata)
			{
				CoordX[i][j]=nodata;
				CoordY[i][j]=nodata;
			}
		}
	}
	if (aktname) free(aktname); aktname=NULL;
	gridinfo[coordxindex].writegrid=1;
	writegridtofile(coordxindex);
	gridinfo[coordxindex].writegrid=0;
	gridinfo[coordyindex].writegrid=1;
	writegridtofile(coordyindex);
	gridinfo[coordyindex].writegrid=0;
}

/* defaultini uebernimmt die Aufgabe, das Hoehenmodell in einen
zu reservierenden Speicherbereich einzulesen, ebenfalls das
Zonengrid. Das Zonengrid wird auf Anzahl der Zonen untersucht, die Arrays
zone und zarray werden erzeugt und belegt. Die Adressen der
Grids und der beiden Arrays werden an den Aufrufer uebergeben, ausserdem
die Anzahl der Zonen und der hoechste auftretende Zonenwert.
das record-Array zone wird mit den Zonen-IDs sowie mit den relativen Anteilen
jeder Zone am Gesamtgrid (ohne NODATA) belegt.
Gleichzeitig werden die grid-Daten (erste 6 saetze)
an den Aufrufer zurueckgeliefert, sowie die Daten von Beginn und
Ende der Modellierung. All diese Daten werden an alle Meteoobjecte und
Teilobjecte zur Initialisierung uebergeben */

// input parameters: Stream, steuernamearg, steuername, version (1 for Topmodel, or 2 for Richards)
// all other variables are return values which are initialized (and partly allocated) within the constructor
defaultini::defaultini(int &retvalue,wasim::mystream *Stream,char *steuernamearg,char *steuername,
	TLastResults *LastResults,int maxLastResults,
	int *LogLevel,int p_version,
	int *progaus,int *warning_level,
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
	irrigationtabletype **p_irrigationtable,MultiLayerLUSEType **p_MultiLayerLUTable,
	bool *WriteASCIIGrids, bool *WriteDoubleGrids, bool *WriteNetcdfGrids, char **DefaultInputDirectory, char **DefaultOutputDirectory,
	long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo) :
	modelobject(Stream, steuernamearg, *LogLevel, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, NULL, NULL, NULL, NULL, NULL, 0, NULL),
LastResults(LastResults),maxLastResults(maxLastResults),
grid(0),gridname(0),zn(0),nFK(0),Na(0),SB(0),DHM(0),LUSE(0),SOIL(0),IRRIG(0),SLP(0),
ASP(0),kWert(0),ARCVEGFRACTIONS(0),AVBIOM(0),WURZ(0),VBG(0),LAI(0),RSC(0),RSE(0),RSI(0),Z0(0),ICE(0),hmax(0),
use(0),naindex(0),kindex(0),nfkindex(0),soilindex(0),slopeindex(0),
aspectindex(0),dhmindex(0),cellsizexread(0),cellsizeyread(0),relcellarearead(0),
gridanzahl(0),do_expokorr(0),do_irrigation(0),
varGridsStart(0),varGridsCount(0),varGridsRead(0),regio1index(0),regio2index(0),regio3index(0),
weight1index(0),weight2index(0),weight3index(0),
REGION1(0),REGION2(0),REGION3(0),WEIGHT1(0),WEIGHT2(0),WEIGHT3(0),
progaus_m(0),warning_level_m(0),runoff_measure_m(0),minutes_to_interval_end_m(0),
anzahlgroessen_m(0),breite_m(0),laenge_m(0),bezugsmeridian_m(0),timekorr_m(0),               
meteonamen_m(0),numsubb_m(0),listsubb_m(0),
//ncols_m(0),nrows_m(0),xll_m(0),yll_m(0),cellsize_m(0),nodata_m(0),
version(0),steuername_m(0),ParamsFunctionState(0),
SemaphoreFileNameRead(NULL),
ExternalSyncCounter(0),
iDoExternalCoupling(0),
lSemaphoreWaitInterval(100),
cCouplingMode('I'),
iExternalTimeStep(0),
iNumExternalGridsRead(0),
ExternalGridNamesRead(NULL),
ExternalGridEntitiesRead(NULL),
ExternalGridsRead(NULL),
iRenameMode(NULL),
iFillMiss(NULL),
iNumExternalGridsWrite(0),
ExternalGridNamesWrite(NULL),
ExternalGridEntitiesWrite(NULL),
ExternalGridsWrite(NULL),
ExternalGridsWriteSum(NULL),
ExternalGridsWriteCodes(NULL),
iNumExternalStatsWrite(0),
ExternalStatNamesWrite(NULL),
ExternalStatEntitiesWrite(NULL),
SemaphoreFileNameWrite(NULL),
SemaphoreContent(NULL),
ExternalStatLastResults(NULL),
ExternalStatsMaxLastResults(0),
ExternalStatsWriteCodes(NULL),
T(NULL),
TAVG(NULL),
TMOVAVG0(NULL),
TMOVAVG7(NULL),
movavgcnt(0),
FORC(NULL),
TSTART(NULL),
CHILL(NULL),
MTSTART(NULL),
FSTAR(NULL),
THROOT(NULL),
VEGSTART(0),
VEGSTOP(0),
VEGDURAT(0),
bUseForcAndTStart(false),
bUseMTStart(false),
iNumMTStartLevels(0),
substepcnt(1),
wurzindex(0),vbgindex(0),laiindex(0),rscindex(0),rsiindex(0),z0index(0),
vegstartindex(0), vegstopindex(0), vegduratindex(0), forcindex(0), chillindex(0), tstartindex(0), fstarindex(0), mtstart_index(0),
MPIREG(0), MPIHLP(0), MPIHALO(0), mpiregindex(0), mpihlpindex(0), mpihaloindex(0),
HRU(0),HRUMaster(0),hruindex(0),hrumasterindex(0), arcvegfractionsindex(0), arcvegbiomassindex(0), arcvegretransindex(0)
{
	try
	{
		retvalue=0;
		int k;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif
		strcpy(instanzname,_Defaults_);
		version=p_version;
		ParamsFunctionState=1000;
		steuername_m = newstring(_stringsize_+1,1);

		strcpy(steuername_m,steuerdateiname);
		replacesymbols(steuername_m,steuerdateiname); // hier wird steuerdateiname um dieEndung $$$ erweitert
		strcpy(steuername_m,steuerdateiname);

		GetBasicControls(steuerdateiname,&progaus_m,&warning_level_m,
		&runoff_measure_m,&minutes_to_interval_end_m,&MinutesToInvEndRouting_m, &timestep_m,
		&anzahlgroessen_m,&breite_m,&laenge_m,
		&bezugsmeridian_m,&timekorr_m,
		&meteonamen_m,&numsubb_m,&listsubb_m);

		min_to_int_end=minutes_to_interval_end_m;
		hour_offset=minutes_to_interval_end_m/60.0;
		gridinfo = *p_gridinfo;

		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		if (index(instanzname," ")>0)  
		instanzname = substr(instanzname,1,index(instanzname," ")-1);

		GetGridNames(); // liesst aus der Steuerdatei die Namen von Zonen- und DHM-Grid sowie dieNamen aller Standardgrids, setzt gridanzahl
		GetDates(); /* die Anfangs- und Endzeit werden hier gelesen und zurueckgegeben  */
		endtime    = jf*jahre+mf*monate+df*tage+hf*stundee+minf*minutee;
		GetMultiLayerLUSETable(version); // liesst die Multilayer-Landuse Table ein, welche auf die normale LU-Tabelle referenziert
		GetLuseTable(version); // liesst die Landnutzungstabelle ein
		GetSoilTable(version); // liesst die Bodenartentabelle ein
		do_irrigation=0;
		GetIrrigTable(version);
		GetExpokorrIgnore(); // liesst den ignore-Parameter des expokorr-Modules
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		DHM = zn = SOIL = nFK = SLP = ASP 
		= SB = Na = kWert = RSE = NULL;

		LUSE = WURZ = VBG = LAI = RSC = RSI = Z0 = NULL;

		LUSE = new double**[1]; // until we know, how many layers there are, space for the first pointer must be reserved

		/* das erste Grid, also immer das Zonengrid, wird diskret eingelesen, da die Daten
		des Grids (nrows, ncols, xll, yll, cellsize, nodata) erst belegt werden muessen, bevor
		readgrdifromfile aufgerufen wird */
		GetGridDimensions(&ncols,&nrows,&xll,&yll,&cellsize,&nodata);
		readgridfromfile(gridinfo[1].name, &(gridinfo[1]).pointer, gridinfo[1].fillcode, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[1].defaultValue,_LimitAutoFillOn_,_MPIFullGridOn_,_MPIRank0OnlyOff_,_MPISyncHaloOn_);
		initialize_grid(1, gridinfo[1].identification, gridinfo[1].pointer, NULL, 0, gridinfo[1].identification, gridinfo[1].writegrid,1,nrows);
		zn = gridinfo[1].pointer; /* von nun an ist das Zonengrid immer ueber zn ansprechbar (in dieser Instanz) */

		MakeMPIRegions(&v_mpiregioncnt, sendHalo, recvHalo, &mpistartrow,&mpiendrow);

		for (k=2;k <= 2+gridanzahl;k++)
		{
			if (gridinfo[k].readcode > 0)
			{
//				--> add additional parameter for defaultValue: if fillcode = 1, use nearestNeighbor, if fillcode = 2, use default value
//				readgridfromfile(gridinfo[k].name,&(gridinfo[k]).pointer,gridinfo[k].fillcode);
				readgridfromfile(gridinfo[k].name,&(gridinfo[k]).pointer,gridinfo[k].fillcode,true,true,gridinfo[k].defaultValue);
				initialize_grid(k,gridinfo[k].identification ,gridinfo[k].pointer,NULL,0,gridinfo[k].identification ,gridinfo[k].writegrid,1,nrows);
			}
			else
			{
				// --> create the grid and fill it with the defaultValue
				gridinfo[k].pointer = matrix(1,nrows,1,ncols);
				initialize_grid(k,gridinfo[k].identification ,gridinfo[k].pointer,NULL,0,gridinfo[k].identification ,gridinfo[k].writegrid,1,nrows);
				double **p = gridinfo[k].pointer; // used as short cut reference only
				for (int i=1;i<=nrows;i++)
				{
					for (int j=1;j<=ncols;j++)
					{
						/* es interessieren die zwei Faelle: 1) ist grid kein Fehlwert, aber zone
							2) ist grid Fehlwert, aber zone nicht */
						if ((p[i][j]!=zz[6])&&(zn[i][j]==nodata)) p[i][j]=nodata;
						if ((p[i][j]==zz[6])&&(zn[i][j]!=nodata)) /* Fehlwert gefunden */
						{
							p[i][j]=gridinfo[k].defaultValue;
						}
					}
				}
			}

			// when this standard grid is periodically read in and is not persistent 
			// then we have to make a copy into TempGrid
			if (gridinfo[k].bPersistency==false)
			{
				gridinfo[k].TempGrid=matrix(1,nrows,1,ncols);
				for (int i=1;i<=nrows;i++)
				{
					for (int j=1;j<=ncols;j++)
					{
						gridinfo[k].TempGrid[i][j]=gridinfo[k].pointer[i][j];
					}
				}
			}
		}
		// testhalber hier schon lesen  
		read_variable_Grids();  

		for (k = 1; k <= 2 + gridanzahl; k++)
		{
			if (gridinfo[k].writegrid == 8)
			{
				writegridtofile(k);
			}
		}


		if (ist_vorhanden(_RelCellArea_))
		RelCellArea=gridinfo[getindex(_RelCellArea_)].pointer;
		else
		RelCellArea=NULL;

		/* Untersuchen der Anzahl und der Reihenfolge der Zonen, gleichzeitig Festhalten der Indizes
		der jeweils ersten und letzten gueltigen Spalte einer Zeile im Array *_FirstLastCol */
		zonenanalyse(&gridinfo,RelCellArea,&ncols,&nrows,&nodata);

		CheckRelCellArea();
		MakeCoordGrids();

		DHM = gridinfo[2].pointer; /* Hoehenmodell von jetzt an unter DHM ansprechbar */
		/* ermitteln der maximalen Hoehe im Hoehenmodell */
		hmax = get_h_maximum(2);

		check_LU_Codes(ist_vorhanden(_Landnutzung_)); /* anschliessend ist LUSE ansprechbar, wenn es existiert */

		if (ist_vorhanden(_Landnutzung_))
		{
			// generate the VBG, Z0, RSC, RSI, RSE, LAI, WURZ and LUSE-grids for all Layers here and initialize them
			int numvegetationlayers;
			if (MultiLayerLUTable) 
			{
				generateMultiLayerLanduse();
				numvegetationlayers=MultiLayerLUTable[1].getNumLayers();
			}
			else
			{ 
				delete[] LUSE;
				LUSE = new double**[2]; // index 0 holds the original LUSE grid as it is used for the old compatibility mode, index 1 ppoints to the same grid for easy use in the submodels
				int luindex=ist_vorhanden(_Landnutzung_);
				LUSE[0] = gridinfo[luindex].pointer;
				LUSE[1] = LUSE[0];
				sprintf(gridinfo[luindex].identification,"%s1",_Landnutzung_);
				numvegetationlayers=1;
			}

			generateMultiLayerParameterGrids(numvegetationlayers); /* wenn ein Landnutzungsgrid eingelesen wurde,
						so wird nun das Wurzeltiefengrid, die LAI, vbg usw. mit Hilfe der LUSE-Daten
						initialisiert */
		}
		else
		{
			delete[] LUSE;
			LUSE=NULL;
		}
		if (do_irrigation) 
		check_irrig_Codes(ist_vorhanden(_irrigation_codes_)); 
		else 
		irrigationtable=NULL;

		check_SOIL_Codes(ist_vorhanden(_Bodenarten_)); /* nun ist auch SOIL ansprechbar, wenn es existiert (der Name _Boenarten_ eingelesen wurde) */

		// either before or after check_SOIL_Codes, the new method for converting the clay-, sand- and silt grids into soil types must be placed.
		makeSoilGridFromClaySiltSand(ist_vorhanden(_Bodenarten_));

		/* Wenn nFK und/oder Na (auffuellbare Porosititaet) NICHT(!!!) als Standardgrids eingelesen wurden,
		so werden nun diese Daten aus dem Bodenartengrid erzeugt, wenn es eingelesen wurde.
		Wurde das Bodenartengrid nicht eingelesen und auch nFK und Na nicht, dann
		gibt es erst dann eine Fehlermeldung, wenn ein Teilmodell nFK oder Na benoetigt
		und nicht findet. */

		if (version==1)
		make_nFK_kWert_and_Na(); /* nFK, kWert, Na sind jetzt bei Bedarf verfuegbar */

		/* Ebenso koennen Exposition und Gefaelle aus dem Hoehenmodell hergeleitet werden, wenn sie
		fuer die Strahlungskorrektur benoetigt werden sollten und nicht als Grids auf der Platte
		zur Verfuegung stehen. Es kann u.U. sinnvoll sein, diese Datensaetze nicht hier erzeugen zu lassen,
		z.B., wenn die Zellgroesse zu gross ist, um eine einigermassen vertrauenswuerdige
		Berechnung zu gewaehrleisten (speziell bei Gefaelle). Das Gefaelle kann dann z.B. aus einem
		Hoehenmodell mit wesentlichn hoeherer Aufloesung berechnet und auf die im Modell benutzte
		Zellgroesse umgewandelt werden, wobei ein Mittelwert zu bilden ist.
		Eine kleine Einschraenkung ergibt sich aus der Begrenzung des Hoehenmodells nach dem Einlesen
		auf die im Zonengrid gueltigen Werte: die Werte der Exposition und des Gefaelles fuer die
		Raender des Einzugsgebietes werden nicht aus dem 3x3-Fenster heraus berechnet (das geht nicht, wenn
		nicht alle 9 Werte da sind), sondern sie werden spaeter aus den berechneten Nachbarwerten
		gemittelt. Auch das kann ein Grund sein, beide Grids extra berechnen zu lassen und als fertige
		Grids einzulesen.

		Ob Exposition und/oder Gefaelle hier berechnet werden, haengt davon ab, ob in der Steuerdatei
		im ABschnitt Standardgrids die Grids mit den Identifikationen "Exposition" bzw. "Gefaelle"
		enthalten sind. Ist das nicht der Fall, werden die Grids erzeugt, ansonsten nicht.

	*/

		make_Slope_and_Aspect(); /* SLP und ASP nun ansprechbar in der Klasse defualtini */

		check_ice_firn_grid();

		MakeRegionWeightGrids();



		/* Lesen der variablen Grids, wenn Leseparameter==0, dann erzeugen von Albedo
		aus Landnutzung und Bodenspeicher aus Feldkapazitaet */
		//    read_variable_Grids();  /* SB ansprechbar, wenn es eingelesen wurde oder erzeugt wurde, sonst natuerlich nicht! */
		
		GetCouplingParameters();

		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
		
		/* Adressen werden uebergeben */
		strcpy(steuername,steuerdateiname);
		*progaus                 = progaus_m;
		*warning_level           = warning_level_m;
		*runoff_measure          = runoff_measure_m;
		*minutes_to_interval_end = minutes_to_interval_end_m;
		*timestep                = timestep_m;
		*anzahlgroessen          = anzahlgroessen_m;
		*breite                  = breite_m;
		*laenge                  = laenge_m;
		*bezugsmeridian          = bezugsmeridian_m;
		*timekorr                = timekorr_m;
		*meteonamen              = meteonamen_m;
		*numsubb                 = numsubb_m;
		*listsubb                = listsubb_m;

		*p_ncols    = ncols;     *p_nrows    = nrows;
		*p_mpistartrow = mpistartrow; *p_mpiendrow = mpiendrow;
		*p_xll      = xll;       *p_yll      = yll;
		*p_cellsize = cellsize;  *p_nodata   = nodata;

		*mina = minutea; *sta = stundea;  *tga = taga;  *mta = monata;  *jra = jahra;
		*mine = minutee; *ste = stundee;  *tge = tage;  *mte = monate;  *jre = jahre;

		*p_FirstLastCol = FirstLastCol;  *p_zone       = zone;
		*p_zarray       = zarray;        *p_zonmax     = zonmax;  
		*p_zonmaxwert   = zonmaxwert;

		*p_rarray       = rarray;        *p_regiomax       = regiomax;
		*p_regiomaxwert = regiomaxwert;  *p_center_of_mass = center_of_mass;

		*p_gridanzahl      = gridanzahl;
		*p_LUtable         = LUtable;
		*p_SOILtable       = SOILtable;
		*p_irrigationtable = version==2?irrigationtable:NULL; /* no irrigation is assumed for the old program version */
		*p_MultiLayerLUTable = MultiLayerLUTable;
		if (steuername_m) free(steuername_m); steuername_m=NULL;
		*WriteASCIIGrids = this->WriteASCIIGrids;
		*WriteDoubleGrids = this->WriteDoubleGrids;
		*WriteNetcdfGrids = this->WriteNetcdfGrids;
		*DefaultInputDirectory = this->DefaultInputDirectory;
		*DefaultOutputDirectory = this->DefaultOutputDirectory;
		nachricht = 1;
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


defaultini::defaultini(int &retvalue,wasim::mystream *Stream,char *steuernamearg,
			TLastResults *LastResults,int maxLastResults,
			int p_version,gridinfotype **p_gridinfo):
modelobject(Stream,steuernamearg,0,false,false,false,NULL,NULL,NULL,NULL,NULL,0,NULL),
LastResults(LastResults),maxLastResults(maxLastResults),
grid(0),gridname(0),zn(0),nFK(0),Na(0),SB(0),DHM(0),LUSE(0),SOIL(0),IRRIG(0),SLP(0),
ASP(0),kWert(0),ARCVEGFRACTIONS(0),AVBIOM(),WURZ(0),VBG(0),LAI(0),RSC(0),RSI(0),RSE(0),Z0(0),ICE(0),hmax(0),
use(0),naindex(0),kindex(0),nfkindex(0),soilindex(0),slopeindex(0),
aspectindex(0),dhmindex(0),cellsizexread(0),cellsizeyread(0),relcellarearead(0),
gridanzahl(0),do_expokorr(0),do_irrigation(0),
varGridsStart(0),varGridsCount(0),varGridsRead(0),regio1index(0),regio2index(0),regio3index(0),
weight1index(0),weight2index(0),weight3index(0),
REGION1(0),REGION2(0),REGION3(0),WEIGHT1(0),WEIGHT2(0),WEIGHT3(0),
progaus_m(0),warning_level_m(0),runoff_measure_m(0),minutes_to_interval_end_m(0),
anzahlgroessen_m(0),breite_m(0),laenge_m(0),bezugsmeridian_m(0),timekorr_m(0),               
meteonamen_m(0),numsubb_m(0),listsubb_m(0),
//ncols_m(0),nrows_m(0),xll_m(0),yll_m(0),cellsize_m(0),nodata_m(0),
version(0),steuername_m(0),ParamsFunctionState(0),
SemaphoreFileNameRead(NULL),
ExternalSyncCounter(0),
lSemaphoreWaitInterval(100),
iDoExternalCoupling(0),
cCouplingMode('I'),
iExternalTimeStep(0),
iNumExternalGridsRead(0),
ExternalGridNamesRead(NULL),
ExternalGridEntitiesRead(NULL),
ExternalGridsRead(NULL),
iRenameMode(NULL),
iFillMiss(NULL),
iNumExternalGridsWrite(0),
ExternalGridNamesWrite(NULL),
ExternalGridEntitiesWrite(NULL),
ExternalGridsWrite(NULL),
ExternalGridsWriteSum(NULL),
ExternalGridsWriteCodes(NULL),
iNumExternalStatsWrite(0),
ExternalStatNamesWrite(NULL),
ExternalStatEntitiesWrite(NULL),
SemaphoreFileNameWrite(NULL),
SemaphoreContent(NULL),
ExternalStatLastResults(NULL),
ExternalStatsMaxLastResults(0),
ExternalStatsWriteCodes(NULL),
T(NULL),
TAVG(NULL),
TMOVAVG0(NULL),
TMOVAVG7(NULL),
movavgcnt(0),
FORC(NULL),
TSTART(NULL),
CHILL(NULL),
MTSTART(NULL),
FSTAR(NULL),
THROOT(NULL),
VEGSTART(0),
VEGSTOP(0),
VEGDURAT(0),
bUseForcAndTStart(false),
bUseMTStart(false),
iNumMTStartLevels(0),
substepcnt(1),
wurzindex(0),vbgindex(0),laiindex(0),rscindex(0),rsiindex(0),z0index(0),
vegstartindex(0), vegstopindex(0), vegduratindex(0), forcindex(0), chillindex(0), tstartindex(0), fstarindex(0), mtstart_index(0),
MPIREG(0), MPIHLP(0), MPIHALO(0), mpiregindex(0), mpihlpindex(0), mpihaloindex(0),
HRU(0),HRUMaster(0),hruindex(0),hrumasterindex(0), arcvegfractionsindex(0), arcvegbiomassindex(0), arcvegretransindex(0)
{
	try
	{
		retvalue=0;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif
		strcpy(instanzname,_Defaults_);
		version=p_version;
		gridinfo = *p_gridinfo;
		steuername_m = newstring(_stringsize_+1,1);
		ParamsFunctionState=0;
		LUSE = new double**[1]; // until we know, how many layers there are, space for the first pointer must be reserved
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

int defaultini::getSubStepCnt()
{
	int i;
	char *satz,*tempsatz;

	satz    = (char *)calloc(_stringsize_+1,sizeof(char));
	tempsatz=satz;
	if ((steuerdatei = fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: Error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}

	if (findeabschnitt(_soil_surface_groundwater_substeps_,false))
	{
		if (fscanf(steuerdatei,"%d",&substepcnt)!=1) 
			steuererror("sub step counter from [soil_surface_groundwater_substeps] not found","");
	}
	else
	{
		substepcnt=1;
	}

	if ((fclose(steuerdatei))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
		throw(11);
	}
	satz=tempsatz;
	if (satz) free(satz);satz=NULL;
	return substepcnt;
}

int defaultini::Parameters(bool *state)
{
	int k;
	char *aktname;
	int gridcnt = 1;
	int sr, er;
	int rivlnkindex=0;
	int dhmindex = 0;

	switch (ParamsFunctionState)
	{
	case 0:
		strcpy(steuername_m, steuerdateiname);
		if (mpi_rank == 0)
		{
			replacesymbols(steuername_m, steuerdateiname); // hier wird steuerdateiname um dieEndung $$$ erweitert
			strcpy(steuername_m, steuerdateiname);
		}
		else
		{
			sprintf(steuername_m, "%s.intern", steuerdateiname);
			strcpy(steuerdateiname, steuername_m);
		}

#ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		GetBasicControls(steuerdateiname, &progaus_m, &warning_level_m,
			&runoff_measure_m, &minutes_to_interval_end_m, &MinutesToInvEndRouting_m, &timestep_m,
			&anzahlgroessen_m, &breite_m, &laenge_m,
			&bezugsmeridian_m, &timekorr_m,
			&meteonamen_m, &numsubb_m, &listsubb_m);
		min_to_int_end = minutes_to_interval_end_m;
		hour_offset = minutes_to_interval_end_m / 60.0;
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 1:

		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}
		if (index(instanzname, " ") > 0)
			instanzname = substr(instanzname, 1, index(instanzname, " ") - 1);

		GetGridNames(); // liesst aus der Steuerdatei die Namen von Zonen- und DHM-Grid sowie dieNamen aller Standardgrids, setzt gridanzahl
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 2:
		GetDates(); /* die Anfangs- und Endzeit werden hier gelesen und zurueckgegeben  */
		endtime = jf*jahre + mf*monate + df*tage + hf*stundee + minf*minutee;
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 3:
		GetMultiLayerLUSETable(version);
		GetLuseTable(version); // liesst die Landnutzungstabelle ein
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 4:
		GetSoilTable(version); // liesst die Bodenartentabelle ein
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 5:
		do_irrigation = 0;
		GetIrrigTable(version);
		GetExpokorrIgnore(); // liesst den ignore-Parameter des expokorr-Modules
		if ((fclose(steuerdatei)) == EOF)
		{
			(*WSstream) << "FATAL: error closing " << steuerdateiname << "[end-entry]";
			throw(11);
		}

		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;
	case 6:
		DHM = zn = SOIL = nFK = SLP = ASP = SB = Na = kWert = RSE = NULL;

		RSC = RSI = Z0 = LAI = WURZ = VBG = 0;

		/* das erste Grid, also immer das Zonengrid, wird diskret eingelesen, da die Daten
		des Grids (nrows, ncols, xll, yll, cellsize, nodata) erst belegt werden muessen, bevor readgridfromfile aufgerufen wird */
		if (mpi_rank == 0)
		{
			GetGridDimensions(&ncols, &nrows, &xll, &yll, &cellsize, &nodata);
		}
#ifdef _MPI
		MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&xll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&yll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cellsize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nodata, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

		// read zone grid
		if (mpi_rank == 0)
		{
			readgridfromfile(gridinfo[1].name, &(gridinfo[1]).pointer, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[1].defaultValue, _LimitAutoFillOff_, _MPIFullGridOn_, _MPIRank0OnlyOn_, _MPISyncHaloOn_);
			initialize_grid(1, gridinfo[1].name, gridinfo[1].pointer, NULL, 0, gridinfo[1].identification, gridinfo[1].writegrid, 1, nrows);
			zn = gridinfo[1].pointer; /* von nun an ist das Zonengrid immer ueber zn ansprechbar (in dieser Instanz) */

			dhmindex = 2;
			readgridfromfile(gridinfo[dhmindex].name, &(gridinfo[dhmindex]).pointer, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[dhmindex].defaultValue, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOn_, _MPISyncHaloOn_);
			initialize_grid(dhmindex, gridinfo[dhmindex].name, gridinfo[dhmindex].pointer, NULL, 0, gridinfo[dhmindex].identification, gridinfo[dhmindex].writegrid, 1, nrows);
			DHM = gridinfo[dhmindex].pointer;

			if (ist_vorhanden(_river_links_, true))
			{ // rivr link grid is required by function MakeMPIRegions (also when not using MPI regions with one instance only). So it has to be read in here already like the zone grid for rank 0 only
				rivlnkindex = getindex(_river_links_);
				readgridfromfile(gridinfo[rivlnkindex].name, &(gridinfo[rivlnkindex]).pointer, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[rivlnkindex].defaultValue, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOn_, _MPISyncHaloOn_);
				initialize_grid(rivlnkindex, gridinfo[rivlnkindex].name, gridinfo[rivlnkindex].pointer, NULL, 0, gridinfo[rivlnkindex].identification, gridinfo[rivlnkindex].writegrid, 1, nrows);
			}
		}
		MakeMPIRegions(&v_mpiregioncnt, &sendHalo, &recvHalo, &mpistartrow, &mpiendrow);
		modelobject::MPIHLP = defaultini::MPIHLP;
		modelobject::MPIREG = defaultini::MPIREG;

		//cout << "mpistartrow " << mpistartrow << " in rank " << mpi_rank << endl;
		//cout << "mpiendrow " << mpiendrow << " in rank " << mpi_rank << endl;
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;

	case 7:
		/* Untersuchen der Anzahl und der Reihenfolge der Zonen, gleichzeitig Festhalten der Indizes
		der jeweils ersten und letzten gueltigen Spalte einer Zeile im Array *_FirstLastCol */

		zonenanalyse(&gridinfo, RelCellArea, &ncols, &nrows, &nodata);

		// print the MPIREG-Grids
#ifdef _MPI

		aktname = newstring(_stringsize_, 1);
		if (mpi_rank == 0)
		{
			sprintf(aktname, "MPIRegions_%d", mpi_rank);
			printgridfile(aktname, MPIREG, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
		}
#ifdef _MPIDEBUG
		sprintf(aktname, "MPIHelper_%d", mpi_rank);
		printgridfile(aktname, MPIHLP, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);

		sprintf(aktname, "MPIHalo_%d", mpi_rank);
		printgridfile(aktname, MPIHALO, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
#endif
		free(aktname);

#endif

		// release grid in instance 0 and re-read the grid, now for all instances
		//		if (mpi_rank == 0)
		//		{
		//			free_matrix(gridinfo[1].pointer, 1, nrows, 1, ncols,nrows);
		//			gridinfo[1].pointer = NULL;
		//		}

		if (ist_vorhanden(_river_links_, true))
		{
			rivlnkindex = getindex(_river_links_);
		}
		else
		{
			rivlnkindex = 0;
		}

		if (mpi_rank > 0)
		{
			sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
			er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;

			gridinfo[1].pointer = matrix(sr, er, 1, ncols, nrows);
			double **p = gridinfo[1].pointer; // used as short cut reference only
			for (int i = sr; i <= er; i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					p[i][j] = nodata;
				}
			}
			initialize_grid(1, gridinfo[1].identification, gridinfo[1].pointer, NULL, 0, gridinfo[1].identification, gridinfo[1].writegrid, sr, er);
			zn = gridinfo[1].pointer;

			gridinfo[2].pointer = matrix(sr, er, 1, ncols, nrows);
			p = gridinfo[2].pointer; // used as short cut reference only
			for (int i = sr; i <= er; i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					p[i][j] = nodata;
				}
			}
			initialize_grid(2, gridinfo[2].identification, gridinfo[2].pointer, NULL, 0, gridinfo[2].identification, gridinfo[2].writegrid, sr, er);
			DHM = gridinfo[2].pointer;

			if (rivlnkindex>0)
			{ // rivr link grid is required by function MakeMPIRegions (also when not using MPI regions with one instance only). So it has to be read in here already like the zone grid for rank 0 only
				gridinfo[rivlnkindex].pointer = matrix(sr, er, 1, ncols, nrows);
				double **p = gridinfo[rivlnkindex].pointer; // used as short cut reference only
				for (int i = sr; i <= er; i++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						p[i][j] = nodata;
					}
				}
				initialize_grid(rivlnkindex, gridinfo[rivlnkindex].identification, gridinfo[rivlnkindex].pointer, NULL, 0, gridinfo[rivlnkindex].identification, gridinfo[rivlnkindex].writegrid, sr, er);
			}

		}
		MPIScatterGrid(1, gridinfo[1].pointer);
		MPISynchronizeGridHalo(1, gridinfo[1].pointer);

		MPIScatterGrid(2, gridinfo[2].pointer);
		MPISynchronizeGridHalo(2, gridinfo[2].pointer);

		if (rivlnkindex > 0)
		{
			MPIScatterGrid(rivlnkindex, gridinfo[rivlnkindex].pointer);
			MPISynchronizeGridHalo(rivlnkindex, gridinfo[rivlnkindex].pointer);
		}

//cout << "before reading zonegrid a 2nd time in rank " << mpi_rank << endl;
//
//cout << "MPIHLP before calling readgridfromfile " << MPIHLP << " in rank " << mpi_rank << endl;
//cout << "MPIREG before calling readgridfromfile " << MPIREG << " in rank " << mpi_rank << endl;

//aktname = newstring(_stringsize_, 1);
//sprintf(aktname, "MPIHLP_%d", mpi_rank);
//printgridfile(aktname, MPIHLP, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//delete[] aktname;
//
//aktname = newstring(_stringsize_, 1);
//sprintf(aktname, "MPIREG_%d", mpi_rank);
//printgridfile(aktname, MPIREG, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//delete[] aktname;

		// here we can read the optimized stripe only to save memory
//cout << "before reading zone grid a 2nd time in rank " << mpi_rank << endl;
//		readgridfromfile(gridinfo[1].name, &(gridinfo[1]).pointer, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[1].defaultValue, _LimitAutoFillOff_, _MPIFullGridOn_);
//		zn = gridinfo[1].pointer; /* von nun an ist das Zonengrid immer ueber zn ansprechbar (in dieser Instanz) */
//cout << "stop after reading zonegrid a 2nd time in rank " << mpi_rank << endl;
//int xxx;
//cin >> xxx;
//		initialize_grid(1, gridinfo[1].name, gridinfo[1].pointer, NULL, 0, gridinfo[1].identification, gridinfo[1].writegrid);

		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;

	case 8:
//cout << " gridanzahl " << gridanzahl << "in rank " << mpi_rank << endl;
#ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		if (ist_vorhanden(_river_links_, true))
		{
			rivlnkindex = getindex(_river_links_);
		}
		else
		{
			rivlnkindex = 0;
		}
		for (k = 3; k <= 2 + gridanzahl; k++)
		{	
			if ((gridinfo[k].pointer == NULL) && (k!=rivlnkindex))
			{
				int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
				int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
				if (mpi_rank == 0)
				{
					sr = 1;
					er = nrows;
				}
				if (gridinfo[k].readcode > 0)
				{
					//				--> add additional parameter for defaultValue: if fillcode = 1, use nearestNeighbor, if fillcode = 2, use default value
					//				readgridfromfile(gridinfo[k].name,&(gridinfo[k]).pointer,gridinfo[k].fillcode);
//cout << " before reading in grid " << gridinfo[k].name << "in rank " << mpi_rank << endl;
					readgridfromfile(gridinfo[k].name, &(gridinfo[k]).pointer, gridinfo[k].fillcode, _AllocateGridMemoryOn_, _ShowMessagesOn_, gridinfo[k].defaultValue, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
					initialize_grid(k, gridinfo[k].name, gridinfo[k].pointer, NULL, 0, gridinfo[k].identification, gridinfo[k].writegrid,sr,er);
//cout << " after reading in grid " << gridinfo[k].name << "in rank " << mpi_rank << endl;

				}
				else
				{
//cout << " creating grid " << gridinfo[k].name << "in rank " << mpi_rank << endl;
					// --> create the grid and fill it with the defaultValue
					gridinfo[k].pointer = matrix(sr, er, 1, ncols,nrows);
					double **p = gridinfo[k].pointer; // used as short cut reference only
					for (int i = sr; i <= er; i++)
					{
						for (int j = 1; j <= ncols; j++)
						{
							/* es interessieren die zwei Faelle: 1) ist grid kein Fehlwert, aber zone 2) ist grid Fehlwert, aber zone nicht */
							if (zn[i][j] == nodata)
							{
								p[i][j] = nodata;
							}
							else
							{
								p[i][j] = gridinfo[k].defaultValue;
							}
						}
					}
					initialize_grid(k, gridinfo[k].identification, gridinfo[k].pointer, NULL, 0, gridinfo[k].identification, gridinfo[k].writegrid,sr,er);
					MPISynchronizeGridHalo(k, gridinfo[k].pointer);
				}

				// when this standard grid is periodically read in and is not persistent 
				// then we have to make a copy into TempGrid
				if (gridinfo[k].bPersistency == false)
				{
					gridinfo[k].TempGrid = matrix(sr, er, 1, ncols,nrows);
					for (int i = sr; i <= er; i++)
					{
						for (int j = 1; j <= ncols; j++)
						{
							gridinfo[k].TempGrid[i][j] = gridinfo[k].pointer[i][j];
						}
					}
				}
//aktname = newstring(_stringsize_, 1);
//sprintf(aktname, "%s_%d", gridinfo[k].identification, mpi_rank);
//printgridfile(aktname, gridinfo[k].pointer, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//free(aktname);
			}
		}
		// testhalber hier schon lesen
		(*WSstream) << "before reading variable grids" << mpi_rank << "[end-entry]";
#ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		read_variable_Grids();

		for (k = 1; k <= 2 + gridanzahl; k++)
		{
			if (gridinfo[k].writegrid == 8)
			{
				writegridtofile(k);
			}
		}

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 9:
		if (ist_vorhanden(_RelCellArea_))
			RelCellArea=gridinfo[getindex(_RelCellArea_)].pointer;
		else
			RelCellArea=NULL;

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 10:
		CheckRelCellArea();

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 11:
		MakeCoordGrids();

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 12:
		DHM = gridinfo[2].pointer; /* Hoehenmodell von jetzt an unter DHM ansprechbar */
		/* ermitteln der maximalen Hoehe im Hoehenmodell */
		hmax = get_h_maximum(2);
		
		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 13:
		check_LU_Codes(ist_vorhanden(_Landnutzung_)); /* anschliessend ist LUSE ansprechbar, wenn es existiert */
		// generate the VBG, Z0, RSC, RSE, RSI, LAI, WURZ and LUSE-grids for all Layers here and initialize them
		if (ist_vorhanden(_Landnutzung_))
		{
			if (MultiLayerLUTable) 
			{
				generateMultiLayerLanduse();
				numvegetationlayers=MultiLayerLUTable[1].getNumLayers();
			}
			else
			{ 
				delete[] LUSE;
				LUSE = new double**[2]; // index 0 holds the original LUSE grid as it is used for the old compatibility mode, index 1 ppoints to the same grid for easy use in the submodels
				int luindex=ist_vorhanden(_Landnutzung_);
				LUSE[0] = gridinfo[luindex].pointer;
				LUSE[1] = LUSE[0];
				sprintf(gridinfo[luindex].identification,"%s1",_Landnutzung_);
				numvegetationlayers=1;
			}
		}
		else
		{
			delete[] LUSE;
			LUSE=NULL;
		}

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;


	case 14:
		if (ist_vorhanden(_Landnutzung_))
		{
			(*WSstream) << "generateMultiLayerParameterGrids..." << "[end-entry]";
			generateMultiLayerParameterGrids(numvegetationlayers); 
			(*WSstream) << "...ok"<<"[end-entry]";
		}
		
		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 15:
		ParamsFunctionState++;
		if (do_irrigation) 
		{
			(*WSstream) << "check_irrig_Codes"<<"[end-entry]";
			check_irrig_Codes(ist_vorhanden(_irrigation_codes_));
			(*WSstream) << "...ok" << "[end-entry]";
		}
		else 
		irrigationtable=NULL;
		state=false;
		return ParamsFunctionState;
		break;

	case 16:
		(*WSstream) << "check_SOIL_Codes..." << "[end-entry]";
		check_SOIL_Codes(ist_vorhanden(_Bodenarten_)); /* nun ist auch SOIL ansprechbar, wenn es existiert (der Name _Boenarten_ eingelesen wurde) */
		(*WSstream) << "...ok" << "[end-entry]";

		/* Wenn nFK und/oder Na (auffuellbare Porosititaet) NICHT(!!!) als Standardgrids eingelesen wurden,
		so werden nun diese Daten aus dem Bodenartengrid erzeugt, wenn es eingelesen wurde.
		Wurde das Bodenartengrid nicht eingelesen und auch nFK und Na nicht, dann
		gibt es erst dann eine Fehlermeldung, wenn ein Teilmodell nFK oder Na benoetigt
		und nicht findet. */
		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;

	case 17:
		(*WSstream) << "make soil grid from clay, sand and silt..." << "[end-entry]";
		makeSoilGridFromClaySiltSand(ist_vorhanden(_Bodenarten_));
		(*WSstream) << "...ok" << "[end-entry]";

		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;

	case 18:
		if (version==1)
		{
			(*WSstream) << "make_nFK_kWert_and_Na..."<<"[end-entry]";
			make_nFK_kWert_and_Na(); /* nFK, kWert, Na sind jetzt bei Bedarf verfuegbar */
			(*WSstream) << "...ok"<<"[end-entry]";
		}
		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 19:
		/* Ebenso koennen Exposition und Gefaelle aus dem Hoehenmodell hergeleitet werden, wenn sie
		fuer die Strahlungskorrektur benoetigt werden sollten und nicht als Grids auf der Platte
		zur Verfuegung stehen. Es kann u.U. sinnvoll sein, diese Datensaetze nicht hier erzeugen zu lassen,
		z.B., wenn die Zellgroesse zu gross ist, um eine einigermassen vertrauenswuerdige
		Berechnung zu gewaehrleisten (speziell bei Gefaelle). Das Gefaelle kann dann z.B. aus einem
		Hoehenmodell mit wesentlichn hoeherer Aufloesung berechnet und auf die im Modell benutzte
		Zellgroesse umgewandelt werden, wobei ein Mittelwert zu bilden ist.
		Eine kleine Einschraenkung ergibt sich aus der Begrenzung des Hoehenmodells nach dem Einlesen
		auf die im Zonengrid gueltigen Werte: die Werte der Exposition und des Gefaelles fuer die
		Raender des Einzugsgebietes werden nicht aus dem 3x3-Fenster heraus berechnet (das geht nicht, wenn
		nicht alle 9 Werte da sind), sondern sie werden spaeter aus den berechneten Nachbarwerten
		gemittelt. Auch das kann ein Grund sein, beide Grids extra berechnen zu lassen und als fertige
		Grids einzulesen.
		Ob Exposition und/oder Gefaelle hier berechnet werden, haengt davon ab, ob in der Steuerdatei
		im ABschnitt Standardgrids die Grids mit den Identifikationen "Exposition" bzw. "Gefaelle"
		enthalten sind. Ist das nicht der Fall, werden die Grids erzeugt, ansonsten nicht.
	*/

		(*WSstream) << "make_Slope_and_Aspect..."<<"[end-entry]";
		make_Slope_and_Aspect(); /* SLP und ASP nun ansprechbar in der Klasse defualtini */
		(*WSstream) << "...ok"<<"[end-entry]";

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 20:

		(*WSstream) << "check_ice_firn_grid..." << "[end-entry]";
		check_ice_firn_grid();
		(*WSstream) << "...ok"<<"[end-entry]";

		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;
		break;

	case 21:

		(*WSstream) << "MakeRegionWeightGrids..."<<"[end-entry]";
		MakeRegionWeightGrids();
		(*WSstream) << "...ok"<<"[end-entry]";

		ParamsFunctionState++;
		state = false;
		return ParamsFunctionState;
		break;

	case 22:
		/* Lesen der variablen Grids, wenn Leseparameter==0, dann erzeugen von Albedo
		aus Landnutzung und Bodenspeicher aus Feldkapazitaet */
		//    read_variable_Grids();  /* SB ansprechbar, wenn es eingelesen wurde oder erzeugt wurde, sonst natuerlich nicht! */
		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;

	case 23:
		// reading parameters for external coupling
		(*WSstream) << "GetCouplingParameters..."<<"[end-entry]";
		GetCouplingParameters();
		(*WSstream) << "...ok"<<"[end-entry]";
		ParamsFunctionState++;
		state=false;
		return ParamsFunctionState;

	default:
		break;
	}
	if ((fclose(steuerdatei))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
		throw(11);
	}
	nachricht = 1;
	*state=true;
	return ParamsFunctionState;
}

void defaultini::CreateStandardGridsStatistics()
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			gridinfo[k].zn=this->zn;
			gridinfo[k].zarray = this->zarray;
			// create statistics output if requested
			// using: 
			// int statcode; // statistics write code (e.g. 2001 for writing the sum for each time step and each sub basin)
			// char *statname; // statistics output file, when configured to write statistics, too
			// statrecord *zone; // statistics regristers per zone (.sx, .n, .zonenwert, .rel, .RelCellSum)
			// zonalstatistic *stat; // handling the output file, including the file header, buffering etc.

			gridinfo[k].zone = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zone) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{
				gridinfo[k].zone[i].n          = gridinfo[k].zone[i].sx         = 0;
				gridinfo[k].zone[i].zonenwert  = zone[i].zonenwert;
				gridinfo[k].zone[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				gridinfo[k].zone[i].RelCellSum = zone[i].RelCellSum;
			}

			int retvalue=0;
			gridinfo[k].stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
				&gridinfo[k].zone,(int)zonmax,(int)gridinfo[k].statcode,(int)5,(char*)gridinfo[k].statname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, timestep<60 && timestep > 0);
			if (retvalue!=0) throw retvalue;
			
			gridinfo[k].stat->schreibzonenwerte(gridinfo[k].identification,"(default)",LastResults,maxLastResults);
		}
	}
}

void defaultini::InitGridsStatistics()
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			gridinfo[k].stat->init();
		}
	}
}

void defaultini::AddValuesGridStatistics()
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					if (zn[i][j]!=nodata)
					{
						gridinfo[k].stat->wertzufuegen(i,j,gridinfo[k].zarray[(int)(gridinfo[k].zn[i][j])],gridinfo[k].pointer[i][j]);
					}
				}
			}
		}
	}
}

void defaultini::CalcAndDoneGridStatistics(long minute,long stunde,long tag,long monat,long jahr)
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			gridinfo[k].stat->calc_and_done(jahr,monat,tag,stunde,minute);
		}
	}
}

void defaultini::FlushBuffersGridStatistics()
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			gridinfo[k].stat->flushbuffer();
		}
	}
}

void defaultini::ClearMemoryGridStatistics()
{
	for (int k=2;k <= 2+gridanzahl;k++)
	{
		if (gridinfo[k].statname != NULL)
		{
			gridinfo[k].stat->ClearMemory();
		}
	}
}

int defaultini::GetParameters(char *steuername,int *LogLevel,
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
	bool *WriteASCIIGrids, bool *WriteDoubleGrids, bool *WriteNetcdfGrids,
	char **DefaultInputDirectory, char **DefaultOutputDirectory, 
	long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo)
{
	/* Adressen werden uebergeben */
	strcpy(steuername,steuername_m);
	*progaus                 = progaus_m;
	*warning_level           = warning_level_m;
	*runoff_measure          = runoff_measure_m;
	*minutes_to_interval_end = minutes_to_interval_end_m;
	*MinutesToInvEndRouting  = MinutesToInvEndRouting_m;
	*timestep                = timestep_m;
	*anzahlgroessen          = anzahlgroessen_m;
	*breite                  = breite_m;
	*laenge                  = laenge_m;
	*bezugsmeridian          = bezugsmeridian_m;
	*timekorr                = timekorr_m;
	*meteonamen              = meteonamen_m;
	*numsubb                 = numsubb_m;
	*listsubb                = listsubb_m;

	*p_ncols    = ncols;     *p_nrows    = nrows;
	*p_mpistartrow = mpistartrow; *p_mpiendrow = mpiendrow;
	*p_xll      = xll;       *p_yll      = yll;
	*p_cellsize = cellsize;  *p_nodata   = nodata;

	*mina = minutea;  *sta = stundea;  *tga = taga;  *mta = monata;  *jra = jahra;
	*mine = minutee; *ste = stundee;  *tge = tage;  *mte = monate;  *jre = jahre;

	*p_FirstLastCol = FirstLastCol;  *p_zone       = zone;
	*p_zarray       = zarray;        *p_zonmax     = zonmax;  
	*p_zonmaxwert   = zonmaxwert;

	*p_rarray       = rarray;        *p_regiomax       = regiomax;
	*p_regiomaxwert = regiomaxwert;  *p_center_of_mass = center_of_mass;

	*p_gridanzahl      = gridanzahl;
	*p_LUtable         = LUtable;
	*p_SOILtable       = SOILtable;
	*p_irrigationtable = version==2?irrigationtable:NULL; /* no irrigation is assumed for the old program version */
	*p_MultiLayerLUTable = MultiLayerLUTable;
	if (steuername_m) free(steuername_m); steuername_m=NULL;
	*WriteASCIIGrids = this->WriteASCIIGrids;
	*WriteDoubleGrids = this->WriteDoubleGrids;
	*WriteNetcdfGrids = this->WriteNetcdfGrids;
	*DefaultInputDirectory = this->DefaultInputDirectory;
	*DefaultOutputDirectory = this->DefaultOutputDirectory;

	*v_regioncnt = this->v_mpiregioncnt;
	*sendHalo = this->sendHalo;
	*recvHalo = this->recvHalo;

	return 0;
}


void defaultini::write_Standard_Grids(long me,long se,long tg,long mt,long jr)
{
	int i;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	// the for loop takes all standard grids and a few MPI grids. 
	// gridanzahl is not the number of standard grids, since all internally generated grids like MPIhalo and -helper and -regions, all phaeno-helper grids etc.
	// are added to gridanzahl. so the nearest value to the number of standardgrids is varGridsStart, sionce the variable gruds were read in directly after reading in the standard grids and creating the MPI helper grids,
	// so varGridsStart is probably only a few grids away from the real number of the standardgrids.
//	for (i = 2; i <= 2+gridanzahl; i++)
	int endindex = gridanzahl;
	if (varGridsStart > 0) endindex = varGridsStart;
	// so if there are vargrids, that index is used, otherwise all grids that are not part of a submodel are looped through
	for (i = 2; i <= endindex; i++)
	{
		writegridtofile(i);
	}
}

void defaultini::write_Variable_Grids(long me,long se,long tg,long mt,long jr)
{
	int i;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	if (varGridsCount>0)
	{
		//    for (i=1;i<=varGridsCount;i++)
		//      writegridtofile(2+gridanzahl-i+1);
		for (i=0;i<varGridsCount;i++)
		writegridtofile(varGridsStart+i);
	}
}

double defaultini::get_hmax()
{
	return hmax;
}


/* Berechnen der aktuellen Pflanzenparameter je nach Vegetationszustand in Abhängigkeit 
vom Julianischen Tag und der Landnutzung u ggf. der Höhe oder Temperatur.
Das ist wichtig zur Berechnung der Verdunstung ueber die nutzbare
Feldkapazitaet nFK und ueber den Inhalt des Bodenspeichers SB.
nach der Berechnung der Wurzeltiefe werden nFK und SB neu berechnet,
ebenso SSP uns SBalt. Auch der unterhalb der ehemals durchwurzelten Zone
(bei Ernte) noch verbliebene Bodenspeicher mit SB<1 wird separat
neu in nFKia und SBia verwaltet */

void defaultini::int_tvar_params(long me,long st, long tg, long mt, long jr, int InV,
double **SBalt,double **SSP,
double **nFKia,double **SBia)
{
	int i,j,k;
	double nFKalt;

	if (bUseForcAndTStart)
	{ // if time step is not 1 day, we must calculate the mean temperatur
		if ((T==NULL)||(TSTART==NULL)||(FORC==NULL))
		{
			int t_index=0;
			char name[200];
			std::stringstream ss,ss2,ss3;
			ss << "temperature";
			ss >> name;

			t_index = ist_vorhanden(name);
			if ((t_index<1)||(t_index>=gridmaxanzahl))
			{
				(*WSstream) << "FATAL: could not find grid with label \""<<name<<"\". \nThis grid must be present (either by interpolating "<<
				"or by reading or as standard grid) for using dynamic phenology"<<"[end-entry]";
				throw(33);
			}
			else
			{
				// now set the pointer to the correct value. T can now be used inside defaultini
				T=gridinfo[t_index].pointer;
			}
		}
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if (zn[i][j]>nodata)
				{
					TAVG[i][j]+=T[i][j]*InV/1440.0;
				}
			}
		}
	}    

	int numlayers;
	if (MultiLayerLUTable)
	{
		numlayers = MultiLayerLUTable[1].getNumLayers();
	}
	else
	{
		numlayers = 1;
	}

	if ((st+(int)hour_offset)==24)
	{

#ifdef _OPENMP
		if (!omp_in_parallel())
		{
			if (nachricht==1)
			(*WSstream) << "INFO: "<< wasim::setw(20) << _actual_parameters_<<" "<<tg<<"."<<mt<<"."<<jr<<" "<<st<<". h"<<"[end-entry]";
			else
			{
				if ((stunde+(int)hour_offset)%nachricht==0) 
				(*WSstream) << "INFO: "<< wasim::setw(20) << _actual_parameters_<<" "<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde<<". h "<<_Meteodatenzeit_<<"[end-entry]";
			}
		}
#else
		if (nachricht==1)
		(*WSstream) << "INFO: "<< wasim::setw(20) << _actual_parameters_<<" "<<tg<<"."<<mt<<"."<<jr<<" "<<st<<". h"<<"[end-entry]";
		else
		{
			if ((stunde+(int)hour_offset)%nachricht==0) 
			(*WSstream) << "INFO: "<< wasim::setw(20) << _actual_parameters_<<" "<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde<<". h "<<_Meteodatenzeit_<<"[end-entry]";
		}
#endif

		// updating moving average for temperature
		if (this->bUseForcAndTStart)
		{ // only if at least one cell has DP2 or DP4 specific land use, the moving average of temperatures is needed
			movavgcnt++;
			if (movavgcnt>7) movavgcnt=7;
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if (zn[i][j]>nodata)
					{ 
						double sum=0.0;
						for (int level=7;level>1;level--)
						{
							TMOVAVG7[i][j][level]=TMOVAVG7[i][j][level-1];
							sum+=TMOVAVG7[i][j][level];
						}
						TMOVAVG7[i][j][1]=TAVG[i][j];
						sum+=TMOVAVG7[i][j][1];
						// calculating moving average now
						TMOVAVG0[i][j] = (double)(sum/movavgcnt); // 7-day-moving average
					}
				}
			}
			writegridtofile(tmovavg7index);
		}

		if ((LUSE)&&(SBalt)&&(SSP)&&(nFKia)&&(SBia)) /* nur dann, wenn auch ein Landnutzungsgrid eingelesen wurde, sollen die nFK und SB-Werte neu berechnet werden */
		{

			/* fuer jede Landnutzungsform wird nun die aktuelle Wurzeltiefe berechnet.
		Dabei werden die Knickpunkttage auf die Hoehe umgerechnet.
		Anschliessend wird daraus die nutzbare Feldkapazitaet und der aktuelle
		Inhalt des Bodenspeichers neu berechnet */

			int tagnr=0;
			for (k=1;k<monat;k++) tagnr += mon[k]; tagnr += (int)tag;

			monat=mt;
			tag=tg;
			jahr=jr;
			stunde=st;

			// check, if TStart was already calculated resp. if the forcingunits were summed up (TStartCalcDone==true)
			// this must be done for land uses with code 3 (method=DynamicPhenology_1)
			if (bUseForcAndTStart)
			{
				// we have found a land use code using methiod 3
				CalculateDynamicPhenology(tagnr,numlayers);
			}

			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if (zn[i][j]>nodata)
					{
						use = (int)LUSE[1][i][j];
						if (use!=nodata)
						{
							ALBEDO[i][j] = LUtable[use].GetAlbedo(WSstream,tag,monat,jahr,DHM[i][j],1,i,j);
							RSE[i][j]    = LUtable[use].GetRse(WSstream,tag,monat,jahr,DHM[i][j],1,i,j);
						}
						else
						{
							ALBEDO[i][j] = nodata;
							RSE[i][j]    = nodata;
						}
						int wurzindex=1;
						double Wurz=WURZ[1][i][j]; // save old value
						for (int layer=1;layer<=numlayers;layer++)
						{
							use = (int)LUSE[layer][i][j];
							if (use!=nodata)
							{
								if (Wurz<WURZ[layer][i][j])
								{
									Wurz = WURZ[layer][i][j]; // if a layer has a deeper root depth, save that value instead
									wurzindex=layer;
								}
								WURZ[layer][i][j]=LUtable[use].GetWurz(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								LAI[layer][i][j] =LUtable[use].GetLAI(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								VBG[layer][i][j] =LUtable[use].GetVbg(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								RSC[layer][i][j] =LUtable[use].GetRsc(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								RSI[layer][i][j] =LUtable[use].GetRsi(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								Z0[layer][i][j]  =LUtable[use].GetZ0(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								
								if (LUtable[use].getJDVegetationReset()==LUtable[use].GetJulDay(tag,monat,jahr))
								{
									VEGSTART[layer][i][j]=-1;
									VEGSTOP[layer][i][j]=-1;
									VEGDURAT[layer][i][j]=0;
								}
								
								if (VEGSTART[layer][i][j]<0)
								{
									VEGSTART[layer][i][j] = LUtable[use].GetVegStart(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								}

								if (VEGSTOP[layer][i][j]<0)
								{
									VEGSTOP[layer][i][j] = LUtable[use].GetVegStop(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								}
								
								/* 1 cm Wurzeltiefe wird als Minimum erachtet (Einfluss der Evaporation)
						in nFK ist immer soviel Wasser, wie 1 cm Wurzeltiefe entspricht */
								if (WURZ[layer][i][j] < double(0.01)) WURZ[layer][i][j]  = double(0.01);
								if (VBG[layer][i][j]  < double(0))    VBG[layer][i][j]   = double(0); 
								else if (VBG[layer][i][j] > double(1)) VBG[layer][i][j] = double(1.0);
								if (LAI[layer][i][j]  < double(0))    LAI[layer][i][j]   = double(0);
								if (RSC[layer][i][j]  < double(0))    RSC[layer][i][j]   = double(0);
								if (RSI[layer][i][j]  < double(0))    RSI[layer][i][j]   = double(0);
								if (Z0[layer][i][j]   < double(0.01))    Z0[layer][i][j]    = double(0.01);

							}
							else
							{
								WURZ[layer][i][j] = nodata;
								LAI[layer][i][j]  = nodata;
								VBG[layer][i][j]  = nodata;
								RSI[layer][i][j]  = nodata;
								Z0[layer][i][j]   = nodata;
								VEGSTART[layer][i][j] = nodata;
								VEGSTOP[layer][i][j] = nodata;
							}

						}

						/* nun Berechnen der Werte für nFK und SB, eventuell noch SSP ? */
						/* nFK entspricht dem maximalen Speicherinhalt des Bodenspeichers */
						nFKalt  = nFK[i][j];

						nFK[i][j] = (double)nFK[i][j]*((double)Wurz/(double)WURZ[wurzindex][i][j]);

						/* wenn nFK groesser wird (Wurzeln wachsen), dann muss untersucht werden,
				ob unterhalb der alten Wurzelgrenze Wassergehalt bei Feldkapazitaet herrschte,
				oder ob, aus alten, bereits geernteten  Pflanzenbestaenden etwa,
				auch dieser Bodenbereich noch nicht bis nFK gesaettigt ist.
				Aus der Bodenfeuchte im alten Wurzelbereich und im neuen Bereich wird dann
				eine neue Bodenfeuchte berechnet. */
						if (nFK[i][j] > nFKalt) /* Wurzeln wachsen */
						{
							if (nFK[i][j] <= (nFKalt+nFKia[i][j]))
							{
								/* SB (rel. Bodenfeuchte bezogen auf nFK) ist eine Mischkonzentration des
						durchwurzelten und des nicht durchwurzelten Bodens */
								SB[i][j]    = (nFKalt*SB[i][j]+(nFK[i][j]-nFKalt)*SBia[i][j])/nFK[i][j];
								nFKia[i][j]-= nFK[i][j]-nFKalt;
								/* SBia[i][j] bleibt so, wie es ist, es ist nur noch auf die kleiner gewordene nFKia bezogen */
							}
							else
							{
								/* Wurzeln wachsen bis unterhalb des nicht bis FC gesaettigten nFKia-Teils */
								SB[i][j]    = (nFKalt*(SB[i][j]-1.0)+nFKia[i][j]*(SBia[i][j]-1.0)+nFK[i][j])/nFK[i][j];
								nFKia[i][j] = 0.0; /* der Speicher beachtenswerter inaktiver Bodenfeuchte verschwindet */
								SBia[i][j]  = 1.0; /* sicherheitshalber wird inaktive relative Bodenfeuchte SBia (bez. auf nFKia) auf 1.0 gesetzt */
							}
						}
						else /* Wurzeln "schrumpfen" (Ernte, Ende der Vegetationsperiode) */
						{
							/* neue nFK ist bereits berechnet worden,
				neue SB bleibt gleich der alten, hier werden nun aber
				SBia und nFKia veraendert */
							if (nFK[i][j]<nFKalt)
							/* im Falle der Gleichheit von nFKalt und nFK
				wird bei nFKia[i][j]==0 der Divisor 0, das darf nicht passieren.
				da nFKia >= 0 und nFKalt - nFK immer >0 ist, wird der Divisor nie <= 0 */
							{
								SBia[i][j] = (nFKia[i][j]*SBia[i][j]+(nFKalt-nFK[i][j])*SB[i][j])/
								(nFKia[i][j]+nFKalt-nFK[i][j]);
								nFKia[i][j] += (nFKalt-nFK[i][j]);
							}
						}
						SSP[i][j]   = 1.0 + nFKalt*(SSP[i][j]-1.0)/nFK[i][j];
						SBalt[i][j] = SB[i][j];
					}
				}
			}
		}
		if (bUseForcAndTStart)
		{ // the Average must be initialized to 0 again for the next day
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if (zn[i][j]>nodata)
					{
						TAVG[i][j]=0;
					}
				}
			}
		}
		if (LUtable[0].getJDVegetationReset()==LUtable[0].GetJulDay(tag,monat,jahr))  {
			int i=0;
			while (gridinfo[i].identification!=NULL)
			{
				gridinfo[i].vegResetFlag=1;
			}
			// set in all grids the flag vegResetFlag (in writegridtofile all grids will be reset and then this flag will be reset to 0)
		} else if (LUtable[0].getJDVegetationWrite()==LUtable[0].GetJulDay(tag,monat,jahr))  {
			// set in all grids the flag vegWriteFlag (will be reset by the writegridtofile-routine itself when grid was written)
			int i=0;
			while (gridinfo[i].identification!=NULL)
			{
				gridinfo[i].vegWriteFlag=1;
			}
		}
	}
	// counting of intervals of the active vegetation period should be done in each interval
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if (zn[i][j]>nodata)
			{
				for (int layer=1;layer<=numlayers;layer++)
				{
					if ((VEGSTART[layer][i][j]>=0)&&(VEGSTOP[layer][i][j]<=0))
					{
						VEGDURAT[layer][i][j]+=1.0; // vegetation period duration is counts one day more
					}
				}
			}
			else
			{
				for (int layer=1;layer<=numlayers;layer++)
				{
					VEGDURAT[layer][i][j] = nodata;
				}
			}
		}
	}
}


// here, we have to do several things:
// 1) check if the required grids exist --> this was NOT done already during the initialisation because 
//    the temperature was not known then (defaultini was initialized before meteorological interpolators
// 3) check for each cell in the TSTART-grid if the value is already set to values of > 0. If yes, skip 
//    summing up forcing units because the phenophase has started already
//    otherwise: check for all other cells in the LUSE grid, if the actual day given in the parameter 
//    tagnr is equal to or greater than the DP1_t1_dorm parameter in the respective LUtable-record
//    if yes, start calculating Forcing units
// 5) check, if the forcing units exceed the maximum F* given as parameter in the LUtable. 
//    If yes, set TSTART-Parameter
// ACHTUNG: funktioniert nur mit Tagesschrittweite, denn diese Berechnung wird nur im letzten Intervall 
// des Tages mit dem dann geltenden Temperaturwert durchgeführt. Für Stundenwerte oder längere Zeitintervalle
void defaultini::CalculateDynamicPhenology(int tagnr,int numlayers)
{
	int i,j;

	// check all landuses if DP2_t0_dorm was reached for method 4 and switch useMethod3 to false
	for (int luz=1;luz<=LUtable[0].getCode();luz++)
	{
		if (LUtable[luz].GetMethod() ==4)
		{
			if (tagnr==LUtable[luz].GetDP2_t0dorm())
			{
				LUtable[luz].SetUseMethod3(false);
			}
		} else if (LUtable[luz].GetMethod() ==6) {// method 6 is DP4, a variation of DP2
			if (tagnr==LUtable[luz].GetDP4_t0dorm())
			{
				LUtable[luz].SetUseMethod3(false);
			}
		}
	}


	double Rf=0; // forcing unit
	double Rc=0; // chilling unit
	// doing the normal calculation of forcing units if necessary
	// using algorithm 12b
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if (zn[i][j]>nodata)
			{
				for (int layer=1;layer<=numlayers;layer++)
				{
					use = (int)LUSE[layer][i][j];
					if (use!=nodata)
					{
						double thetascaling=0.0;
						if (THROOT)
						{
							if (THROOT[i][j]>=LUtable[use].getReduParDry()) // it is more dry than Field capacity
							{
								thetascaling=(THROOT[i][j]-LUtable[use].getReduParDry())/(150.0 - LUtable[use].getReduParDry());
							}
							else
							{
								thetascaling=0; // no stress effects if soil is wetter than fieldcapacity
							}
						}
						if ((LUtable[use].GetMethod()==3)||(LUtable[use].GetUseMethod3()==true))
						{ // only for cells using dynamic phenology we will do the calculation of forcing units
							if (   (TSTART[layer][i][j]<0)			 // TSTART is not reached yet, so we have to go on
									&&(tagnr>=LUtable[use].GetDP1_t1dorm())) // start calculation of forcing units after DP1_t1_dorm was reached only
							{
								// now we know, that we are after DP1_t1_dorm and TSTART was not set, so this means, that the forcing untis must still be calculated
								if (TAVG[i][j]>LUtable[use].GetDP1_TBf())
								{
									Rf=28.4/(1.0+exp(-0.185*(TAVG[i][j]-LUtable[use].GetDP1_TBf()-18.4)));
								}
								else
								{
									Rf=0;
								}
								double sf = LUtable[use].GetStressFactor();
								if (sf>0)
								{
									FORC[layer][i][j]=FORC[layer][i][j]+(Rf*(1.0+thetascaling*sf));
								}
								else
								{
									FORC[layer][i][j]+=Rf;
								}
								if (FORC[layer][i][j]>LUtable[use].GetFstar())
								{ // now we can set TSTART
									TSTART[layer][i][j]=tagnr;
									if ((LUtable[use].GetMethod()==5)||(LUtable[use].GetMethod()==6))
									{
										int nextJDindex=1;
										// ermitteln der Stützstelle, für die als nächstes ein STarttag berechnet werden muss
										// for method 6 we use both, TSTART and MTSTART
										while ((MTSTART[layer][i][j][nextJDindex]<tagnr)
										&&(MTSTART[layer][i][j][nextJDindex]>0)
										&&(MTSTART[layer][i][j][nextJDindex+1]>nodata+1))
										{
											nextJDindex++;
										}
										MTSTART[layer][i][j][nextJDindex]=tagnr;
										MTSTART[layer][i][j][nextJDindex+1]=tagnr+LUtable[use].GetJulianDays()[nextJDindex+1];
									}
								}
							}
							if (tagnr==LUtable[use].GetJDReset_TStart())	{
								// z.B. m letzten Tag des Jahres oder dem ersten Tag werden die Forcingunits auf 0 gesetzt und TStart auf -1, damit wasim erkennt, dass dieser Zeitpunkt wieder dynamisch bestimmt werden muss
								FORC[layer][i][j]=0.0;
								TSTART[layer][i][j]=-1.0;

								if (STARTFALL) STARTFALL[i][j]=0.0;
								if (MTSTART)
								{
									for (int lv=1;lv<=LUtable[use].getNumFixPoints();lv++)
									{
										if ((LUtable[use].GetForcingThreshold(lv)<0)||(tagnr >= LUtable[use].GetJulianDay(lv)))
										{
											MTSTART[layer][i][j][lv] = LUtable[use].GetJulianDay(lv);
										}
										else
										{
											MTSTART[layer][i][j][lv] = -1;
										}
										// initialise the FSTAR-Value for interpolating the parameters to zero at the beginning of a new year
										// MTSTART[layer][i][j][0]=0;
									}
								}
							}
						}
						else if (LUtable[use].GetMethod()==4)
						{ // only for cells using dynamic phenology 2 we will do the calculation of chilling and forcing units here
							if ((TSTART[layer][i][j]<0)			 // TSTART is not reached yet, so we have to go on
									&&(tagnr>=LUtable[use].GetDP2_t1dorm()) // start calculation of forcing units after DP2_t1_dorm was reached only -> DP2_t1_dorm is usually different from DP1_t1_dorm!
									) 
							{
								// now we know, that we are after DP1_t1_dorm and TSTART was not set, so this means, that the forcing untis must still be calculated
								if (TAVG[i][j]>LUtable[use].GetDP2_TBf()) // DP2_T_Bf is different from DP1_T_Bf (possibly)
								{
									Rf=28.4/(1.0+exp(-0.185*(TAVG[i][j]-LUtable[use].GetDP2_TBf()-18.4)));
								}
								else
								{
									Rf=0;
								}
								double sf = LUtable[use].GetStressFactor();
								if (sf>0)
								{
									FORC[layer][i][j]=FORC[layer][i][j]+(Rf*(1.0+thetascaling*sf));
									//									FORC[layer][i][j]+=(Rf*thetascaling*sf);
								}
								else
								{
									FORC[layer][i][j]+=Rf;
								}
								if (FORC[layer][i][j]>FSTAR[layer][i][j]) // here we use FSTAR as an array because it is set dynamically and not fix for a landuse
								{ // now we can set TSTART
									TSTART[layer][i][j]=tagnr;
								}
							}
							
							// calculate chilling untis for northern Hemisphere 
							// (southern Hemisphere: if ((tagnr>=LUtable[use].GetDP2_t0dorm())&&(tagnr<LUtable[use].GetDP2_t1dorm()))
							if ((tagnr>=LUtable[use].GetDP2_t0dorm())||(tagnr<LUtable[use].GetDP2_t1dorm()))
							{
								if (TAVG[i][j]<=LUtable[use].GetDP2_Offset_1()) {
									Rc=0.0;
								} else if ((TAVG[i][j]<=LUtable[use].GetDP2_TBc())&&
										(TAVG[i][j]<LUtable[use].GetDP2_Offset_2())) {
									Rc=(TAVG[i][j]-LUtable[use].GetDP2_Offset_1())/
									(LUtable[use].GetDP2_TBc()-LUtable[use].GetDP2_Offset_1());
								} else if (TAVG[i][j]<LUtable[use].GetDP2_Offset_2()) {
									Rc=(TAVG[i][j]-LUtable[use].GetDP2_Offset_2())/
									(LUtable[use].GetDP2_TBc()-LUtable[use].GetDP2_Offset_2());
								} else {
									Rc=0.0;
								}
								CHILL[layer][i][j]+=Rc;
							}

							// at the specified day of year, the forcing units and theTStart-value must be reset
							if (tagnr==LUtable[use].GetJDReset_TStart())	
							{
								FORC[layer][i][j]=0.0;
								TSTART[layer][i][j]=-1.0;
							} 
							
							// at the last day before DP2_t0_dorm the Chilling units must be reset to 0
							if (tagnr==LUtable[use].GetDP2_t0dorm()-1)
							{
								CHILL[layer][i][j]=0.0;
							} 
							
							// at the last day before DP2_t1_dorm, the FSTAR-values must be estimated
							// no chilling units are calculated afterwards until DP2_t0_dorm
							// nach F_star=a*exp(b*C_star)
							if (tagnr==LUtable[use].GetDP2_t1dorm()-1)
							{
								FSTAR[layer][i][j]=LUtable[use].GetDP2_Par_a()*exp(LUtable[use].GetDP2_Par_b()*CHILL[layer][i][j]);
							}

						}
						else if (LUtable[use].GetMethod()==5)
						{	
							// Adding up Forcingunits
							// looking for the actual day in the MTSTART-Stack: is it somewhere 
							// in a sequence between at least one dynamic day and another dynamic or fix day?
							// if yes, check if Forcing units are sufficent to set the day, if not check if maximum
							// day is reached
							int nextJDindex=1;
							// ermitteln der Stützstelle, für die als nächstes ein STarttag berechnet werden muss
							// --> entfällt, wenn dieser Wert nicht -1 ist oder tagnr > dem JD für die Stützstelle ist
							while ((MTSTART[layer][i][j][nextJDindex]<tagnr)
							&&(MTSTART[layer][i][j][nextJDindex]>0)
							&&(MTSTART[layer][i][j][nextJDindex+1]>nodata+1))
							{
								nextJDindex++;
							}
							// was, wenn aktueller Tag >= Stützstelle, aber MTSTART noch -1 --> dann MTSTART auf JD setzen
							if ((LUtable[use].GetJulianDay(nextJDindex)<=tagnr)
									&&(LUtable[use].GetForcingThreshold(nextJDindex)>0))
							{
								MTSTART[layer][i][j][nextJDindex]=LUtable[use].GetJulianDay(nextJDindex);
								MTSTART[layer][i][j][0]=LUtable[use].GetForcingThreshold(nextJDindex);
							}

							// now we know, that we are after DP1_t1_dorm and TSTART was not set, so this means, that the forcing untis must still be calculated
							if (TAVG[i][j]>LUtable[use].GetDP1_TBf())
							{
								Rf=TAVG[i][j]-LUtable[use].GetDP1_TBf();
								// Rf=28.4/(1.0+exp(-0.185*(TAVG[i][j]-LUtable[use].GetDP1_TBf()-18.4)));
							}
							else
							{
								Rf=0;
							}
							double sf = LUtable[use].GetStressFactor();
							if (sf>0)
							{
								FORC[layer][i][j]=FORC[layer][i][j]+(Rf*(1.0+thetascaling*sf));
								//								FORC[layer][i][j]+=(Rf*thetascaling*sf);
							}
							else
							{
								FORC[layer][i][j]+=Rf;
							}

							if (MTSTART[layer][i][j][nextJDindex]<0)
							{
								if (FORC[layer][i][j]>=LUtable[use].GetForcingThreshold(nextJDindex))
								{ // now we can set TSTART for the nextJDindex-day
									MTSTART[layer][i][j][nextJDindex]=tagnr;
									// and store the actual threshold for interpolating between a dynamic or static day and a dynamic day
									MTSTART[layer][i][j][0]=LUtable[use].GetForcingThreshold(nextJDindex);
								}
							}
							if (tagnr==1)	{
								// am ersten Tag des Jahres werden die Forcingunits auf 0 gesetzt und alle MTStart auf -1 
								// oder den entsprechenden JulDay, damit wasim erkennt, dass diese Zeitpunkte wieder 
								// dynamisch bestimmt werden muessen
								FORC[layer][i][j]=0.0;
								MTSTART[layer][i][j][0]=FORC[layer][i][j]; // because the interpolation between the first dynamic day and the origin must start from 0 as FORC starts from 0, too
								for (int lv=1;lv<=LUtable[use].getNumFixPoints();lv++)
								{
									if ((LUtable[use].GetForcingThreshold(lv)<0)||(tagnr >= LUtable[use].GetJulianDay(lv)))
									{
										MTSTART[layer][i][j][lv] = LUtable[use].GetJulianDay(lv);
									}
									else
									{
										MTSTART[layer][i][j][lv] = -1;
									}
									// initialise the FSTAR-Value for interpolating the parameters to zero at the beginning of a new year
									// MTSTART[layer][i][j][0]=0;
								}
							} // tagnr==1
						} // method = 5
						else if (LUtable[use].GetMethod()==6)
						{ // only for cells using dynamic phenology 2 we will do the calculation of chilling and forcing units here

							if ((TSTART[layer][i][j]<0)			 // TSTART is not reached yet, so we have to go on
									&&(tagnr>=LUtable[use].GetDP4_t1dorm()) // start calculation of forcing units after DP2_t1_dorm was reached only -> DP2_t1_dorm is usually different from DP1_t1_dorm!
									) 
							{
								int nextJDindex=1;
								// ermitteln der Stützstelle, für die als nächstes ein STarttag berechnet werden muss
								// for method 6 we use both, TSTART and MTSTART
								while ((MTSTART[layer][i][j][nextJDindex]<tagnr)
								&&(MTSTART[layer][i][j][nextJDindex]>0)
								&&(MTSTART[layer][i][j][nextJDindex+1]>nodata+1))
								{
									nextJDindex++;
								}
								// now we know, that we are after DP4_t1_dorm and TSTART was not set, so this means, that the forcing untis must still be calculated
								if (TAVG[i][j]>=LUtable[use].GetDP4_T1()) // DP2_T_Bf is different from DP1_T_Bf (possibly)
								{
									Rf=TAVG[i][j]-LUtable[use].GetDP4_T1();
								}
								else
								{
									Rf=0;
								}
								double sf = LUtable[use].GetStressFactor();
								if (sf>0)
								{
									FORC[layer][i][j]=FORC[layer][i][j]+(Rf*(1.0+thetascaling*sf));
									//									FORC[layer][i][j]+=(Rf*thetascaling*sf);
								}
								else
								{
									FORC[layer][i][j]+=Rf;
								}
								if (FORC[layer][i][j]>FSTAR[layer][i][j]) // here we use FSTAR as an array because it is set dynamically and not fix for a landuse
								{ // now we can set TSTART
									TSTART[layer][i][j]=tagnr;
									MTSTART[layer][i][j][nextJDindex]=tagnr;
									MTSTART[layer][i][j][nextJDindex+1]=tagnr+LUtable[use].GetJulianDays()[nextJDindex+1];
								}
							}
							
							// calculate chilling untis for northern Hemisphere 
							// (southern Hemisphere: if ((tagnr>=LUtable[use].GetDP2_t0dorm())&&(tagnr<LUtable[use].GetDP2_t1dorm()))

							if ((tagnr>=LUtable[use].GetDP4_t0dorm())||(TSTART[layer][i][j]<0)) // correctly implemented, CD must be summed up until TSTART (no stop at DP4_t_dorm)
							//							if ((tagnr>=LUtable[use].GetDP4_t0dorm())||(tagnr<LUtable[use].GetDP4_t1dorm()))
							{
								if (TAVG[i][j]<=LUtable[use].GetDP4_T0()) {
									Rc=1.0;
								} else {
									Rc=0.0;
								}
								CHILL[layer][i][j]+=Rc;
							}

							// at the specified day of year, the forcing units and theTStart-value must be reset
							if (tagnr==LUtable[use].GetJDReset_TStart())	
							{
								FORC[layer][i][j]=0.0;
								TSTART[layer][i][j]=-1.0;
								STARTFALL[i][j]=0.0;
								for (int lv=1;lv<=LUtable[use].getNumFixPoints();lv++)
								{
									if ((LUtable[use].GetForcingThreshold(lv)<0)||(tagnr >= LUtable[use].GetJulianDay(lv)))
									{
										MTSTART[layer][i][j][lv] = LUtable[use].GetJulianDay(lv);
									}
									else
									{
										MTSTART[layer][i][j][lv] = -1;
									}
									// initialise the FSTAR-Value for interpolating the parameters to zero at the beginning of a new year
									// MTSTART[layer][i][j][0]=0;
								}
							} 
							
							// at the last day before DP2_t0_dorm the Chilling units must be reset to 0
							if (tagnr==LUtable[use].GetDP4_t0dorm()-1)
							{
								CHILL[layer][i][j]=0.0;
							} 
							
							// at the last day before DP4_t1_dorm, the FSTAR-values must be estimated
							// no chilling units are calculated afterwards until DP4_t0_dorm
							// nach F_star=a+b*ln(TScrit) with TSCrit calculated here as FSTAR
							if ((tagnr>=LUtable[use].GetDP4_t1dorm()-1)&&(TSTART[layer][i][j]<0))
							{
								//								FSTAR[layer][i][j]=LUtable[use].GetDP2_Par_a()*exp(LUtable[use].GetDP2_Par_b()*CHILL[layer][i][j]);
								if (CHILL[layer][i][j]>0)
								{
									FSTAR[layer][i][j]=LUtable[use].GetDP4_Par_a()+LUtable[use].GetDP4_Par_b()*log(CHILL[layer][i][j]);
								}
								else
								{
									FSTAR[layer][i][j]=LUtable[use].GetDP4_Par_a(); // CHILL will be assumed to be 1, log(CHILL) will be 0, so the Par_a will remain
								}
							}

							// DP4 has an extended functionallity: leave fall is estimated by a moving average falling below a threshold for at least 5 days
							if ((tagnr>=LUtable[use].GetDP4_txs1()) && (STARTFALL[i][j]<=0)&&(TSTART[layer][i][j]>0)) // will be executed after vegetation start and before leave fall start only
							{ // we are at or after the date, the temperature moving average for start of leave fall must be evaluated
								if (TMOVAVG0[i][j]<LUtable[use].GetDP4_Txylstop())
								{
									TMOVAVG7[i][j][0]+=1.0;
								}
								else
								{
									TMOVAVG7[i][j][0]=0.0; // reset if moving average is above threshold temperature
								}

								// threshold value for start leave fall is reached
								// or the latest possible day (fix) was reached
								int actsp=LUtable[use].GetDP4_SPSLF();
								if ((TMOVAVG7[i][j][0]>=5.0)||(tagnr>=LUtable[use].GetJulianDay(actsp))) 
								{ // now set sample points to the actual day
									if (actsp<=LUtable[use].getNumFixPoints())
									{
										MTSTART[layer][i][j][actsp+1] = tagnr+MTSTART[layer][i][j][actsp+1]-MTSTART[layer][i][j][actsp];
										MTSTART[layer][i][j][actsp] = tagnr;
										for (int spcount=1;spcount<actsp;spcount++)
										{
											if ((MTSTART[layer][i][j][spcount]>tagnr)||
													(MTSTART[layer][i][j][spcount]<0)||
													(MTSTART[layer][i][j][spcount]<MTSTART[layer][i][j][spcount-1]))
											{
												MTSTART[layer][i][j][spcount]=tagnr-1;
											}
										}
									}
									STARTFALL[i][j]=tagnr;
								}
							}

						} // method=6
					} // use != nodata
				} // layer
			} // zone != nodata
		} // j
	} // i
}


//   Berechnen der aktuellen Pflanzenparameter je nach Vegetationszustand in Abhängigkeit 
//   vom Julianischen Tag und der Landnutzung u ggf. der Höhe oder Temperatur.
//   hier in der Version für Richards.
void defaultini::int_tvar_params(long me,long st, long tg, long mt, long jr,int InV)
{
	int i,j,k,tstart;
	//bool TStartCalcDone=false;

	if (bUseForcAndTStart)
	{ // if time step is not 1 day, we must calculate the mean temperatur
		if ((T==NULL)||(TSTART==NULL)||(FORC==NULL))
		{
			int t_index=0;
			char name[200];
			std::stringstream ss,ss2,ss3;
			ss << "temperature";
			ss >> name;

			t_index = ist_vorhanden(name);
			if ((t_index<1)||(t_index>=gridmaxanzahl))
			{
				(*WSstream) << "FATAL: could not find grid with label \""<<name<<"\". \nThis grid must be present (either by interpolating "<<
				"or by reading or as standard grid) for using dynamic phenology"<<"[end-entry]";
				throw(33);
			}
			else
			{
				// now set the pointer to the correct value. T can now be used inside defaultini
				T=gridinfo[t_index].pointer;
			}

			if (ist_vorhanden(_THROOT_grid_))
			THROOT=gridinfo[getindex(_THROOT_grid_)].pointer;
			else
			THROOT=NULL;

		}
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank))
				{
					TAVG[i][j]+=T[i][j]*InV/1440.0;
				}
			}
		}
	}    

	int numlayers;
	if (MultiLayerLUTable)
	{
		numlayers = MultiLayerLUTable[1].getNumLayers();
	}
	else
	{
		numlayers = 1;
	}

	if ((60 * st + (int)((me + 60.0*hour_offset)))==1440) // days end
	//if ((st+(int)hour_offset)==24)
	{

		if (((60 * st + (int)((me + 60.0*hour_offset))) % (60 * nachricht) == 0) && (st != 0))
		{
			//if ((stunde+(int)hour_offset)%nachricht==0)
			if ((InV > 0) && (InV < 60))
			{
				(*WSstream) << "INFO: " << wasim::setw(20) << _actual_parameters_ << " " << tg << "." << mt << "." << jr << " " << st << ":" << me << " h " << _Meteodatenzeit_ << "[end-entry]";
			}
			else
			{
				(*WSstream) << "INFO: " << wasim::setw(20) << _actual_parameters_ << " " << tg << "." << mt << "." << jr << " " << st << ". h " << _Meteodatenzeit_ << "[end-entry]";
			}
		}


		// updating moving average for temperature
		if (bUseForcAndTStart)
		{
			movavgcnt++;
			if (movavgcnt>7) movavgcnt=7;
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank))
					//				  if (zn[i][j]>nodata)
					{ 
						double sum=0.0;
						for (int level=7;level>1;level--)
						{
							TMOVAVG7[i][j][level]=TMOVAVG7[i][j][level-1];
							sum+=TMOVAVG7[i][j][level];
						}
						TMOVAVG7[i][j][1]=TAVG[i][j];
						sum+=TMOVAVG7[i][j][1];
						// calculating moving average now
						TMOVAVG0[i][j] = (double)(sum/movavgcnt); // 7-day-moving average
					}
				}
			}
			writegridtofile(tmovavg7index);
		}

		if (LUSE) /* nur dann, wenn auch ein Landnutzungsgrid eingelesen wurde, sollen die Parameter neu berechnet werden */
		{
			/* fuer jede Landnutzungsform wird nun die aktuelle Wurzeltiefe berechnet.
		Dabei werden die Knickpunkttage auf die Hoehe umgerechnet. Diese Information wird
		im Bodenmodell benoetigt, um die Evapotranspiration aus den Bodenschichten richtig
		entnehmen zu koennen */
			monat=mt;
			tag=tg;
			jahr=jr;
			stunde=st;

			int tagnr=0;
			for (k=1;k<monat;k++) tagnr += mon[k]; tagnr += (int)tag;

			// check, if TStart was already calculated resp. if the forcingunits were summed up (TStartCalcDone==true)
			// this must be done for land uses with code 3 (method=DynamicPhenology_1)
			if (bUseForcAndTStart)
			{
				// we have found a land use code using methiod 3
				// no parameters like LAI or so are touched in this function
				CalculateDynamicPhenology(tagnr,numlayers);
			}

			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank))
					//			if (zn[i][j]>nodata)
					{
						for (int layer=1;layer<=numlayers;layer++)
						{
							use = (int)LUSE[layer][i][j];
							if (use!=nodata)
							{
								if (layer==1) 
								{ 
									// only if the flag for this particular entity and the particular land use is set in the flag-array 
									// (sDynFlags=11111111 or dynflag=0xFF = default: all entities should be calculated here)
									if (LUtable[use].isFlagSet(ALBEDOFLAG))
										ALBEDO[i][j] = LUtable[use].GetAlbedo(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
									if (LUtable[use].isFlagSet(RSEFLAG))
										RSE[i][j]    = LUtable[use].GetRse(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								}
								if (LUtable[use].isFlagSet(ROOTFLAG))
									WURZ[layer][i][j]=LUtable[use].GetWurz(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
								if (LUtable[use].isFlagSet(LAIFLAG))
									LAI[layer][i][j] =LUtable[use].GetLAI(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
								if (LUtable[use].isFlagSet(VCFFLAG))
									VBG[layer][i][j] =LUtable[use].GetVbg(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
								if (LUtable[use].isFlagSet(RSCFLAG))
									RSC[layer][i][j] =LUtable[use].GetRsc(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
								if (LUtable[use].isFlagSet(RSIFLAG))
									RSI[layer][i][j] =LUtable[use].GetRsi(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 
								if (LUtable[use].isFlagSet(Z0FLAG))
									Z0[layer][i][j]  =LUtable[use].GetZ0(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j); 

								if (LUtable[use].getJDVegetationReset()==LUtable[use].GetJulDay(tag,monat,jahr))
								{
									VEGSTART[layer][i][j]=-1;
									VEGSTOP[layer][i][j]=-1;
									VEGDURAT[layer][i][j]=0;
								}

								if (VEGSTART[layer][i][j]<=0)
								{
									VEGSTART[layer][i][j] = LUtable[use].GetVegStart(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								}

								if (VEGSTOP[layer][i][j]<0)
								{
									VEGSTOP[layer][i][j] = LUtable[use].GetVegStop(WSstream,tag,monat,jahr,DHM[i][j],layer,i,j);
								}
								
								/* 1 cm Wurzeltiefe wird als Minimum erachtet (Einfluss der Evaporation)
					in nFK ist immer soviel Wasser, wie 1 cm Wurzeltiefe entspricht */
								if (WURZ[layer][i][j] < double(0.01)) WURZ[layer][i][j]  = double(0.01);
								if (VBG[layer][i][j]  < double(0))    VBG[layer][i][j]   = double(0); 
								else if (VBG[layer][i][j] > double(1)) VBG[layer][i][j] = double(1.0);
								if (LAI[layer][i][j]  < double(0))    LAI[layer][i][j]   = double(0);
								if (RSI[layer][i][j]  < double(0))    RSI[layer][i][j]   = double(0);
								if (Z0[layer][i][j]   <= double(0))    Z0[layer][i][j]    = double(0.01);
							}
							else
							{
								if (layer==1) 
								{
									ALBEDO[i][j]=nodata;
									RSE[i][j]=nodata;
								}
								WURZ[layer][i][j] = nodata;
								LAI[layer][i][j]  = nodata;
								VBG[layer][i][j]  = nodata;
								RSI[layer][i][j]  = nodata;
								RSC[layer][i][j]  = nodata;
								Z0[layer][i][j]   = nodata;
								VEGSTART[layer][i][j] = nodata;
								VEGSTOP[layer][i][j] = nodata;
							}
						}
					}
				}
			}
		}
		if (bUseForcAndTStart)
		{ // the Average must be initialized to 0 again for the next day
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank))
					//				  if (zn[i][j]>nodata)
					{
						TAVG[i][j]=0;
					}
				}
			}
		}

		if (LUtable[0].getJDVegetationReset()==LUtable[0].GetJulDay(tag,monat,jahr))  
		{
			int i=1;
			while (gridinfo[i].identification!=NULL)
			{
				gridinfo[i].vegResetFlag=1;
				i++;
			}
			// set in all grids the flag vegResetFlag (in writegridtofile all grids will be reset and then this flag will be reset to 0)
		} 
		
		if (LUtable[0].getJDVegetationWrite()==LUtable[0].GetJulDay(tag,monat,jahr))  
		{
			// set in all grids the flag vegWriteFlag (will be reset by the writegridtofile-routine itself when grid was written)
			int i=1;
			while (gridinfo[i].identification!=NULL)
			{
				gridinfo[i].vegWriteFlag=1;
				i++;
			}
		}
	}

	// counting of intervals of the active vegetation period should be done in each interval
	if (LUSE!=NULL)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank))
				//		  if (zn[i][j]>nodata)
				{
					for (int layer=1;layer<=numlayers;layer++)
					{
						if ((VEGSTART[layer][i][j]>=0)&&(VEGSTOP[layer][i][j]<=0))
						{
							VEGDURAT[layer][i][j]+=1.0; // vegetation period duration is counts one day more
						}
					}
				}
				else
				{
					if (MPIREG[i][j]==mpi_rank)
					{
						for (int layer=1;layer<=numlayers;layer++)
						{
							VEGDURAT[layer][i][j] = nodata;
						}
					}
				}
			}
		}
	}
}

/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void defaultini::ClearMemory()
{
	int i;

	if  (iDoExternalCoupling)
	{
		if (SemaphoreFileNameRead) free(SemaphoreFileNameRead);   SemaphoreFileNameRead = NULL;
		if (SemaphoreFileNameWrite) free(SemaphoreFileNameWrite); SemaphoreFileNameWrite = NULL;
		if (SemaphoreContent) free(SemaphoreContent);             SemaphoreContent = NULL;

		for (i=0;i<iNumExternalGridsRead;i++)
		{
			if (ExternalGridNamesRead[i])    delete[] ExternalGridNamesRead[i]; ExternalGridNamesRead[i]=NULL;
			if (ExternalGridEntitiesRead[i]) delete[] ExternalGridEntitiesRead[i]; ExternalGridEntitiesRead[i]=NULL;
		}
		if (ExternalGridNamesRead)    delete[] ExternalGridNamesRead; ExternalGridNamesRead=NULL;
		if (ExternalGridEntitiesRead) delete[] ExternalGridEntitiesRead; ExternalGridEntitiesRead=NULL;
		if (ExternalGridsRead)        delete[] ExternalGridsRead; ExternalGridsRead=NULL;
		if (iRenameMode)              delete[] iRenameMode; iRenameMode=NULL;
		if (iFillMiss)                delete[] iFillMiss; iFillMiss=NULL;

		for (i=0;i<iNumExternalGridsWrite;i++)
		{
			if (ExternalGridNamesWrite[i]) delete[] ExternalGridNamesWrite[i]; ExternalGridNamesWrite[i]=NULL;
			if (ExternalGridEntitiesWrite[i]) delete[] ExternalGridEntitiesWrite[i]; ExternalGridEntitiesWrite[i]=NULL;
			if (ExternalGridsWriteCodes[i]) delete[] ExternalGridsWriteCodes[i]; ExternalGridsWriteCodes[i]=NULL;
			if (ExternalGridsWriteSum[i]) 
			{
				free_matrix(ExternalGridsWriteSum[i],1,nrows,1,ncols);
			}
		}
		if (ExternalGridNamesWrite) delete[] ExternalGridNamesWrite; ExternalGridNamesWrite=NULL;
		if (ExternalGridEntitiesWrite) delete[] ExternalGridEntitiesWrite; ExternalGridEntitiesWrite=NULL;
		if (ExternalGridsWrite) delete[] ExternalGridsWrite; ExternalGridsWrite=NULL;
		if (ExternalGridsWriteSum) delete[] ExternalGridsWriteSum; ExternalGridsWriteSum=NULL;
		if (ExternalGridsWriteCodes) delete[] ExternalGridsWriteCodes; ExternalGridsWriteCodes=NULL;
		

		for (i=0;i<iNumExternalStatsWrite;i++)
		{
			if (ExternalStatNamesWrite[i]) delete[] ExternalStatNamesWrite[i]; ExternalStatNamesWrite[i]=NULL;
			if (ExternalStatEntitiesWrite[i]) delete[] ExternalStatEntitiesWrite[i]; ExternalStatEntitiesWrite[i]=NULL;
			if (ExternalStatsWriteCodes[i]) delete[] ExternalStatsWriteCodes[i]; ExternalStatsWriteCodes[i]=NULL;
		}
		if (ExternalStatNamesWrite) delete[] ExternalStatNamesWrite; ExternalStatNamesWrite=NULL;
		if (ExternalStatEntitiesWrite) delete[] ExternalStatEntitiesWrite; ExternalStatEntitiesWrite=NULL;
		if (ExternalStatsWriteCodes) delete[] ExternalStatsWriteCodes; ExternalStatsWriteCodes=NULL;
		

		if (SemaphoreFileNameWrite) delete[] SemaphoreFileNameWrite; SemaphoreFileNameWrite=NULL;
		if (SemaphoreContent) delete[] SemaphoreContent; SemaphoreContent=NULL;
	} 

	if (RiverLinkCellList.size() > 0) { RiverLinkCellList.clear(); }
	if (RiverLinkList.size() > 0) { RiverLinkList.clear(); }
	if (RiverLinkVector.size() > 0) 
	{ 
		for (std::vector<RiverLink*>::iterator iter = RiverLinkVector.begin(); iter != RiverLinkVector.end(); iter++)
		{
			(*iter)->linkCellVector->clear();
			delete (*iter)->linkCellVector;
		}
		RiverLinkVector.clear();
	}
	
	for (std::vector<std::vector<RiverLinkCell> >::iterator iter1 = RivLnkCellsToSend.begin(); iter1 != RivLnkCellsToSend.end(); iter1++)
	{
		iter1->clear();
	}
	RivLnkCellsToSend.clear();

	for (std::vector<std::vector<RiverLinkCell> >::iterator iter1 = RivLnkCellsToReceive.begin(); iter1 != RivLnkCellsToReceive.end(); iter1++)
	{
		iter1->clear();
	}
	RivLnkCellsToReceive.clear();

	if (MultiLayerLUTable)
	{
		delete[] MultiLayerLUTable;
		MultiLayerLUTable = NULL;
	}

	if (LUtable) 
	{
		//	  int usecnt=LUtable[0]->getCode();
		//	  for (int z=1;z<=usecnt;z++)
		//	  {
		//		  delete LUtable[z];
		//		  LUtable[z]=NULL;
		//	  }
		delete[] LUtable;
		LUtable=NULL;
	}

	//  if (SOILtable)
	//  {
	//    for (i=1;i<=SOILtable[0].Code;i++)
	//    {
	//      if (SOILtable[i].Name) free(SOILtable[i].Name);
	//      SOILtable[i].Name=NULL;
	//
	//      if (version==2)
	//      {
	//        if (SOILtable[i].dzui    ) free(SOILtable[i].dzui    ); SOILtable[i].dzui    =NULL;
	//        if (SOILtable[i].k_sat   ) free(SOILtable[i].k_sat   ); SOILtable[i].k_sat   =NULL;
	//        if (SOILtable[i].theta_k ) free(SOILtable[i].theta_k ); SOILtable[i].theta_k =NULL;
	//        if (SOILtable[i].k_tab   ) free(SOILtable[i].k_tab   ); SOILtable[i].k_tab   =NULL;
	//        if (SOILtable[i].theta_h ) free(SOILtable[i].theta_h ); SOILtable[i].theta_h =NULL;
	//        if (SOILtable[i].h_tab   ) free(SOILtable[i].h_tab   ); SOILtable[i].h_tab   =NULL;
	//      }
	//    }
	//  }
	//  if (SOILtable) free(SOILtable); SOILtable=NULL;

	if (SOILtable) 
	{
		delete[] SOILtable; 
		SOILtable=NULL;
	}

	if (irrigationtable)
	{
		for (i=1;i<=irrigationtable[0].Code;i++)
		{
			if (irrigationtable[i].Name)   free(irrigationtable[i].Name);   
			irrigationtable[i].Name   = NULL;
			if (irrigationtable[i].month)  free(irrigationtable[i].month);  
			irrigationtable[i].month  = NULL;
			if (irrigationtable[i].day)    free(irrigationtable[i].day);    
			irrigationtable[i].day    = NULL;
			if (irrigationtable[i].amount) free(irrigationtable[i].amount); 
			irrigationtable[i].amount = NULL;
		}
		if (irrigationtable) free(irrigationtable); irrigationtable=NULL;
	}

	if (ARCVEGFRACTIONS) { delete[] ARCVEGFRACTIONS;   ARCVEGFRACTIONS = NULL; }
	if (AVBIOM) { delete[] AVBIOM;   AVBIOM = NULL; }
	if (LUSE)   {delete[] LUSE;   LUSE   = NULL;}
	if (WURZ)   {delete[] WURZ;   WURZ   = NULL;}
	if (VBG)	  {delete[] VBG;    VBG    = NULL;}
	if (LAI)	  {delete[] LAI;    LAI    = NULL;}
	if (RSC)	  {delete[] RSC;    RSC    = NULL;}
	if (RSI)	  {delete[] RSI;    RSI    = NULL;}
	if (Z0)	  {delete[] Z0;     Z0     = NULL;}
	if (VEGSTART){delete[] VEGSTART; VEGSTART = NULL;}
	if (VEGSTOP){delete[] VEGSTOP; VEGSTOP = NULL;}
	if (VEGDURAT){delete[] VEGDURAT; VEGDURAT = NULL;}
	if (FORC)	  {delete[] FORC;   FORC   = NULL;}
	if (TSTART) {delete[] TSTART; TSTART = NULL;}
	if (FSTAR)  {delete[] FSTAR;  FSTAR  = NULL;}
	if (CHILL)  {delete[] CHILL;  CHILL  = NULL;}

	if (arcvegfractionsindex) { delete[] arcvegfractionsindex;   arcvegfractionsindex = NULL; }
	if (arcvegbiomassindex) { delete[] arcvegbiomassindex;   arcvegbiomassindex = NULL; }
	if (arcvegretransindex) { delete[] arcvegretransindex;   arcvegretransindex = NULL; }
	if (wurzindex) {delete[] wurzindex;wurzindex = NULL;}
	if (vbgindex)  {delete[] vbgindex ;vbgindex  = NULL;}
	if (laiindex)  {delete[] laiindex ;laiindex  = NULL;}
	if (rscindex)  {delete[] rscindex ;rscindex  = NULL;}
	if (rsiindex)  {delete[] rsiindex ;rsiindex  = NULL;}
	if (z0index)   {delete[] z0index  ;z0index   = NULL;}
	if (vegstartindex) {delete[] vegstartindex ;vegstartindex = NULL;}
	if (vegstopindex)  {delete[] vegstopindex  ;vegstopindex  = NULL;}
	if (vegduratindex) {delete[] vegduratindex ;vegduratindex = NULL;}
	if (forcindex)     {delete[] forcindex	   ;forcindex     = NULL;}
	if (chillindex)    {delete[] chillindex    ;chillindex    = NULL;}
	if (tstartindex)   {delete[] tstartindex   ;tstartindex   = NULL;}
	if (fstarindex)    {delete[] fstarindex    ;fstarindex    = NULL;}
	if (mtstart_index) {delete[] mtstart_index ;mtstart_index = NULL;}

	if (mtstart_index1>0)
	{
		char n1[255];
		sprintf(n1,"%s_End",gridinfo[mtstart_index1].name);
		printstackfile(n1, MTSTART[1], iNumMTStartLevels + 1, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	}

	if (MTSTART) {delete[] MTSTART; MTSTART = NULL;}

	if (varGridsCount > 0)
	{
		if (varGridsRead) free(varGridsRead); 
		varGridsRead=NULL;
	}
	this->ClearMemoryGridStatistics();
	for (int i=1;i<=2+gridanzahl;i++)
	{
		// check, ob Stack belegt ist, wenn ja, dann auch den Stack löschen
		free_matrix(gridinfo[i].pointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
		if (gridinfo[i].stackpointer != NULL)
		{
			this->free_drei_d_matrix(gridinfo[i].stackpointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, gridinfo[i].levels, nrows);
		}
		delete_grid(i);
		if (gridinfo[i].name)     free (gridinfo[i].name);
		if (gridinfo[i].outname)  free (gridinfo[i].outname);
		if (gridinfo[i].statname) free (gridinfo[i].statname);
        if (gridinfo[i].identification) free (gridinfo[i].identification);
	}

	if (FirstLastCol) free(FirstLastCol); 
	FirstLastCol=NULL;
	if (zarray) free(zarray);
	zarray=NULL;
	if (zone) free(zone);
	zone=NULL;
	if (rarray)
	{
		if (rarray)         free(rarray);
		rarray=NULL;
		if (center_of_mass) free(center_of_mass);
		center_of_mass=NULL;
	}
	if (DefaultInputDirectory)  free(DefaultInputDirectory);  DefaultInputDirectory=NULL;
	if (DefaultOutputDirectory) free(DefaultOutputDirectory); DefaultOutputDirectory=NULL;
//	if (v_mpiregioncnt) delete[] v_mpiregioncnt;
	modelobject::ClearMemory();
}

/* defining tracer_ini */

tracer_ini_class::tracer_ini_class(int &retvalue,wasim::mystream *Stream,char *steuername, int LogLevel,
	int *p_tracer_count,tracer_info_type **p_tracer_info,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,int p_doSinUp,std::vector<HRUMapElement*> *vHRU):
	modelobject(Stream, steuername, LogLevel, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
	v_mpiregioncnt,p_sendHalo,p_recvHalo,p_doSinUp,vHRU),
tracer_count(0),tracer_info(0)
{
	try
	{
		retvalue=0;
		signed long i,k;
		char *equal_sign,*testsatz;
		char *measure;
		strcpy(instanzname,_tracer_ini_);
		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		tracer_count=0;
		tracer_info=NULL;
		findeabschnitt("[" _substance_transport_ "]");
		fscanfln(steuerdatei,"%i",&tracer_count," number of tracers in ",tstr);
		if (tracer_count>9)
		{
			(*WSstream) << "FATAL: not more than 9 tracer species are allowed at one time (they need a very lot of memory!)"<<"[end-entry]";
			throw(58);
		}
		if (tracer_count>0)
		{
			tracer_info = (tracer_info_type*)calloc(tracer_count+1,sizeof(tracer_info_type)); /* allocating memory for the tracer-info-table */
			equal_sign=(char*)calloc(_stringsize_,sizeof(char));
			measure=(char*)calloc(_stringsize_,sizeof(char));
			for (i=1;i<=tracer_count;i++)
			{ /* reading parameters for each tracer type */
				(*WSstream) << "reading parameters for tracer "<<i<<" ..."<<"[end-entry]";
				tracer_info[i].name=(char*)calloc(100,sizeof(char));
				satz=fgets(satz,_stringsize_,steuerdatei);
				sscanf(satz,"%3s ",tracer_info[i].name); /* getting the name of the tracer into the tracer_info-array */
				if (strlen(tracer_info[i].name) > 3)
				{
					tracer_info[i].name[3] = '\0'; /* limiting the name to the first 3 chars */
				}
				else
				{
					for (k = strlen(tracer_info[i].name); k < 3; k++)
					{
						tracer_info[i].name[k] = '_';
					}
				}
				/* making sure, no tracer has the same name like another (in the first 3 digits) */
				for (k = 0; k < (signed long)strlen(tracer_info[i].name); k++)
				{
					if ((tracer_info[i].name[k] >= 'A') && (tracer_info[i].name[k] <= 'Z'))
					{
						tracer_info[i].name[k] = (char)((unsigned int)'a' + (unsigned int)tracer_info[i].name[k] - (unsigned int)'A');
					}
				}
				for (k = 1; k < i; k++)
				{
					if (!strcmp(tracer_info[i].name, tracer_info[k].name))
					{
						(*WSstream) << "FATAL: \n----------- ERROR in substance_transport ---------\n"
							"\nsubstance #" << i << " and substance #" << k << " have the same\n"
							"first 3 chars (lowercase) in its name (" << tracer_info[i].name
							<< " <-> " << tracer_info[k].name << ")\n"
							"\n--------------------------------------------------\n" << "[end-entry]";
						throw(59);
					}
				}

				/* NOW SEARCHING FOR KEY-WORDS */
				/* radioactive? 1 = yes, 2 = no, 3 = no entry */
				if ((strstr(satz,_radioactive_))&&(!strstr(satz,_non_radioactive_))) tracer_info[i].radioactive=1;
				else
				{
					if (strstr(satz,_non_radioactive_)) tracer_info[i].radioactive=0;
					else
					{
						(*WSstream) << "FATAL: error reading parameters for tracer "<<i<<":\n need a key word \""<<_radioactive_
						<<"\" or \""<<_non_radioactive_<<"\"\n"
						"(may be, the number of tracers ("<<tracer_count<<") is wrong?\n"<<"[end-entry]";
						throw(56); }
				}
				/* half time */
				if (tracer_info[i].radioactive)
				{
					if (strstr(satz,_half_time_))
					{
						testsatz=strstr(strstr(satz,_half_time_),"=");
						if (!testsatz)
						{
							(*WSstream) << "FATAL: error reading half_time for tracer "<<i<<": no \""<<_half_time_<<" =\" found!"<<"[end-entry]";
							throw(56);
						}
						else testsatz=testsatz+1;
						k=sscanf(testsatz,"%lf%s",&(tracer_info[i].half_time),measure);
						if (k!=2)
						{
							(*WSstream) << "FATAL: error reading "<<_half_time_<<" for tracer "<<i<<"[end-entry]";
							throw(56);
						}
						if (strcmp(measure,"d")!=0)
						{
							(*WSstream) << "FATAL: error reading "<<_half_time_<<" for tracer "<<i<<": no \"d\" for measure found!"<<"[end-entry]";
							throw(56);
						}
						if (tracer_info[i].half_time<0)
						{
							(*WSstream) << "FATAL: error in "<<_half_time_<<" for tracer "<<i<<": must be greater than 0!"<<"[end-entry]";
							throw(56);
						}
					}
					else
					{
						(*WSstream) << "FATAL: error in tracer_info_table tracer "<<i<<": need a \""<<_half_time_<<" = xxxx d\" sequence!"<<"[end-entry]";
						throw(56);
					}
				}
				else tracer_info[i].half_time=-9999.0;

				/* evaporating ? */
				if ((strstr(satz,_evaporating_))&&(!strstr(satz,_non_evaporating_)))
				{
					tracer_info[i].evaporating=1;
				}
				else
				{
					if (strstr(satz,_non_evaporating_)) tracer_info[i].evaporating=0;
					else
					{
						(*WSstream) << "FATAL: error reading parameters for tracer "<<i<<":\n need a key word \""<<_evaporating_
						<<"\" or \""<<_non_evaporating_<<"\"\n"
						"(may be, the number of tracers ("<<tracer_count<<") is wrong?"<<"[end-entry]";
						throw(56);
					}
				}
				sprintf(measure,"");

				/* minimum concentrations */
				if (strstr(satz,_min_conc_))
				{
					testsatz=strstr(strstr(satz,_min_conc_),"=");
					if (!testsatz)
					{
						(*WSstream) << "FATAL: error reading "<<_min_conc_<<" for tracer "<<i<<": no \""<<_min_conc_<<" =\" found!"<<"[end-entry]";
						throw(56);
					}
					else testsatz=testsatz+1;
					if (sscanf(testsatz,"%lf",&(tracer_info[i].min_conc))!=1)
					{
						(*WSstream) << "FATAL: error reading "<<_min_conc_<<" for tracer "<<i<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading parameters for tracer "<<i<<":\nneed a key word \""<<_min_conc_<<" = xxxx\""<<"[end-entry]";
					throw(56);
				}

				/* maximum concentrations */
				if (strstr(satz,_max_conc_))
				{
					testsatz=strstr(strstr(satz,_max_conc_),"=");
					if (!testsatz)
					{
						(*WSstream) << "FATAL: error reading max_conc for tracer "<<i<<": no \""<<_max_conc_<<" =\" found!"<<"[end-entry]";
						throw(56);
					}
					else testsatz=testsatz+1;
					if (sscanf(testsatz,"%lf%s",&(tracer_info[i].max_conc),measure)!=2)
					{ sprintf(measure,"kg/kg"); }
				}
				else
				{
					(*WSstream) << "FATAL: error reading parameters for tracer "<<i<<":\nneed a key word \""<<_max_conc_<<" = xxxx\""<<"[end-entry]";
					throw(56);
				}
				if (strcmp(measure,"mg/l")==0)
				{
					//        tracer_info[i].min_conc*=1e-6;
					//        tracer_info[i].max_conc*=1e-6;
				} /* because standard concentration measure is 1/1 and 1 mg/l is only around 1e-6 */
				else
				if (strcmp(measure,"g/l")==0)
				{
					//          tracer_info[i].min_conc*=1e-3;
					//          tracer_info[i].max_conc*=1e-3;
				} /* because standard concentration measure is 1/1 and 1 g/l is only around 1e-3 */
				else
				if ((strcmp(measure,"kg/kg")!=0)&&(strcmp(measure,"Kg/Kg")!=0)) /* kg/kg is the standard measuer (relative mass fraction) */
				{
					(*WSstream) << "FATAL: error reading max_conc for tracer "<<i
					<<":\nwrong measure (only mg/l, g/l or kg/kg resp. no measure are allowed!)\n"
					"read here: "<<measure<<"[end-entry]";
					throw(56);
				}

				/* initial concentrations within the soil (down to groundwater-1 layer distribution) */
				if (strstr(satz,_soilini_))
				{
					testsatz=strstr(strstr(satz,_soilini_),"=");
					if (!testsatz)
					{
						(*WSstream) << "FATAL: error reading "<<_soilini_<<" for tracer "<<i<<": no \""<<_soilini_<<" =\" found!"<<"[end-entry]";
						throw(56);
					}
					else testsatz=testsatz+1;
					if (sscanf(testsatz,"%lf",&(tracer_info[i].soilini))!=1)
					{
						(*WSstream) << "FATAL: error reading "<<_soilini_<<" for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading parameters for tracer "<<i
					<<":\nneed a key word \""<<_soilini_<<" = xxxx\""<<"[end-entry]";
					throw(56);
				}

				/* initial concentrations within the soil (down to bedrock equal distribution) */
				if (strstr(satz,_gwini_))
				{
					testsatz=strstr(strstr(satz,_gwini_),"=");
					if (!testsatz)
					{
						(*WSstream) << "FATAL: error reading "<<_gwini_<<" for tracer "<<i
						<<": no \""<<_gwini_<<" =\" found!"<<"[end-entry]";
						throw(56);
					}
					else testsatz=testsatz+1;
					if (sscanf(testsatz,"%lf",&(tracer_info[i].gwini))!=1)
					{
						(*WSstream) << "FATAL: error reading "<<_gwini_<<" for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading parameters for tracer "<<i
					<<":\nneed a key word \""<<_gwini_<<" = xxxx\""<<"[end-entry]";
					throw(56);
				}

				/* statistic codes and grid write codes */
				if (strstr(satz,"statcode ="))
				{
					testsatz=strstr(satz,"statcode =");
					testsatz=testsatz+10;
					if (sscanf(testsatz,"%d",&(tracer_info[i].statcode))!=1)
					{
						(*WSstream) << "FATAL: error reading statcode for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading statcode for tracer "<<i
					<<":\nneed a key word \"statcode = xxxx\""<<"[end-entry]";
					throw(56);
				}
				if (strstr(satz,"gridcode ="))
				{
					testsatz=strstr(satz,"gridcode =");
					testsatz=testsatz+10;
					if (sscanf(testsatz,"%d",&(tracer_info[i].writegrid))!=1)
					{
						(*WSstream) << "FATAL: error reading gridcode for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading gridcode for tracer "<<i
					<<":\nneed a key word \"gridcode = xxxx\""<<"[end-entry]";
					throw(56);
				}
				if (strstr(satz,"path ="))
				{
					testsatz=strstr(satz,"path =");
					testsatz=testsatz+6;
					tracer_info[i].path=newstring(300,1);
					if (sscanf(testsatz,"%299s",tracer_info[i].path)!=1)
					{
						(*WSstream) << "FATAL: error reading path for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading path for tracer "<<i
					<<":\nneed a key word \"path = xxxx\""<<"[end-entry]";
					throw(56);
				}
				if (strstr(satz,"statext ="))
				{
					testsatz=strstr(satz,"statext =");
					testsatz=testsatz+9;
					tracer_info[i].statext=newstring(300,1);
					if (sscanf(testsatz,"%299s",tracer_info[i].statext)!=1)
					{
						(*WSstream) << "FATAL: error reading statext for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading statext for tracer "<<i
					<<":\nneed a key word \"statext = xxxx\""<<"[end-entry]";
					throw(56);
				}
				if (strstr(satz,"gridext ="))
				{
					testsatz=strstr(satz,"gridext =");
					testsatz=testsatz+9;
					tracer_info[i].gridext=newstring(300,1);
					if (sscanf(testsatz,"%299s",tracer_info[i].gridext)!=1)
					{
						(*WSstream) << "FATAL: error reading gridext for tracer "<<i<<"!"<<"[end-entry]";
						throw(56);
					}
				}
				else
				{
					(*WSstream) << "FATAL: error reading gridext for tracer "<<i
					<<":\nneed a key word \"gridext = xxxx\""<<"[end-entry]";
					throw(56);
				}
				(*WSstream) << "done."<<"[end-entry]";
			}
			if (measure) if (measure) free(measure);
			measure=NULL;
			if (equal_sign) if (equal_sign) free(equal_sign);
			equal_sign=NULL;
		}
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
		/* passing values to return parameters */
		*p_tracer_count = tracer_count;
		*p_tracer_info  = tracer_info;
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


void tracer_ini_class::ClearMemory()
{
	int i;
	if (tracer_count>0)
	{
		for (i=1;i<=tracer_count;i++)
		{
			if (tracer_info[i].name) free(tracer_info[i].name);
			tracer_info[i].name=NULL;
			if (tracer_info[i].path) free(tracer_info[i].path);
			tracer_info[i].path=NULL;
			if (tracer_info[i].statext) free(tracer_info[i].statext);
			tracer_info[i].statext=NULL;
			if (tracer_info[i].gridext) free(tracer_info[i].gridext);
			tracer_info[i].gridext=NULL;
		}
		if (tracer_info) free(tracer_info);
	}
	modelobject::ClearMemory();
}


void defaultini::GetCouplingParameters()
{    
	char *tstr;

	tstr=new char[_stringsize_];
	
	// some parameters for coupling the model with an external groundwatermodel
	iDoExternalCoupling=0; // already initialized, but we do it a second time here
	if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
		throw(14);
	}
	if (findeabschnitt(_ExternalCoupling_,false))
	{
		if (fscanfln(steuerdatei,"%d",&iDoExternalCoupling,"Coupling mode ",tstr)!=1)
		{
			(*WSstream) << "FATAL: -----------------error------------\n"
			"error reading coupling mode for "<<instanzname<<"\n"
			"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block\n"<<"[end-entry]";
			throw(10);
		}
		if (iDoExternalCoupling)
		{
//			if (mpi_size > 1)
//			{
//				(*WSstream) << "FATAL: -----------------error------------\n"
//					"external coupling does not work when running in MPI mode with more than 1 instance (i.e. you should use the non-MPI-version)" << "[end-entry]";
//				throw(1009);
//			}
			SemaphoreFileNameRead  = newstring(_stringsize_+1,1);
			SemaphoreFileNameWrite = newstring(_stringsize_+1,1);
			SemaphoreContent       = newstring(_stringsize_+1,1);
			
			if (fscanfln(steuerdatei,"%s",SemaphoreFileNameRead,"name of the semaphore file ",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading name of the external groundwater table grid for "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			if (fscanfln(steuerdatei,"%ld",&lSemaphoreWaitInterval,"wait interval for scanning for external grids ",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading of wait interval for scanning for external gw table grid "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			fscanfln(steuerdatei,"%c",&cCouplingMode,"Coupling mode for external grids (I, H, D, M or Y) ",tstr);
			if ((cCouplingMode!='I')&&
					(cCouplingMode!='H')&&
					(cCouplingMode!='D')&&
					(cCouplingMode!='M')&&
					(cCouplingMode!='Y'))
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading Coupling Mode (I, H, D, M or Y) for external gw table grid "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			if (fscanfln(steuerdatei,"%d",&iExternalTimeStep,"the external model's timestep in minutes ",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading time step of the external model (should be given in minutes after coupling mode) "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			if (fscanfln(steuerdatei,"%d",&iNumExternalGridsRead,"number external grids (to read)",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading number of external grids to read in "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			// allocating memory for the external gridnames, entities and grid pointers, rename modes and fillmiss codes
			if (iNumExternalGridsRead>0)
			{
				ExternalGridNamesRead    = new char*[iNumExternalGridsRead];
				ExternalGridEntitiesRead = new char*[iNumExternalGridsRead];
				ExternalGridsRead        = new double**[iNumExternalGridsRead];
				iRenameMode              = new int[iNumExternalGridsRead];
				iFillMiss                = new int[iNumExternalGridsRead];
				for (int i=0;i<iNumExternalGridsRead;i++)
				{
					ExternalGridNamesRead[i] = new char[_stringsize_+1];
					ExternalGridEntitiesRead[i] = new char[_stringsize_+1];
					if (fscanf(steuerdatei,"%s%s%d%d",&(ExternalGridNamesRead[i][0]),&(ExternalGridEntitiesRead[i][0]),&(iFillMiss[i]),&(iRenameMode[i]))!=4)
					{
						(*WSstream) << "FATAL: -----------------error------------\n"
						"error reading external (read) grid name/entity/fill code/rename code "<<instanzname<<"\n"
						"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
						throw(10);
					}
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
			if (fscanfln(steuerdatei,"%d",&iNumExternalGridsWrite,"number external grids (to write)",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading number of external grids to write in "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			// allocating memory for the external gridnames, entities and grid pointers, rename modes and fillmiss codes
			if (iNumExternalGridsWrite>0)
			{
				ExternalGridNamesWrite    = new char*[iNumExternalGridsWrite];
				ExternalGridEntitiesWrite = new char*[iNumExternalGridsWrite];
				ExternalGridsWrite        = new double**[iNumExternalGridsWrite];
				ExternalGridsWriteSum     = new double**[iNumExternalGridsWrite];
				ExternalGridsWriteCodes    = new char*[iNumExternalGridsWrite];
				for (int i=0;i<iNumExternalGridsWrite;i++)
				{
					ExternalGridNamesWrite[i] = new char[_stringsize_+1];
					ExternalGridEntitiesWrite[i] = new char[_stringsize_+1];
					ExternalGridsWriteCodes[i] = new char[50];
					if (fscanf(steuerdatei,"%s%s%s",
								&(ExternalGridNamesWrite[i][0]),
								&(ExternalGridEntitiesWrite[i][0]),
								&(ExternalGridsWriteCodes[i][0]))!=3)
					{
						(*WSstream) << "FATAL: -----------------error------------\n"
						"error reading external (write) grid name/entity/Writecode"<<instanzname<<"\n"
						"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
						throw(10);
					}
					if ((strcmp(ExternalGridsWriteCodes[i],"Mean")==0)||
							(strcmp(ExternalGridsWriteCodes[i],"Sum")==0)||
							(strncmp(ExternalGridsWriteCodes[i],"SumTotal",8)==0))
					{
						ExternalGridsWriteSum[i]=matrix(1,nrows,1,ncols);
						for (int m = mpistartrow; m <= mpiendrow; m++)
						{
							for (int n = 1; n <= ncols; n++)
							{
								ExternalGridsWriteSum[i][m][n] = (zn[m][n] == nodata) ? nodata : 0;
							}
						}
					}
					else
					{
						ExternalGridsWriteSum[i]=NULL;
					}
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
			if (fscanfln(steuerdatei,"%d",&iNumExternalStatsWrite,"number external Stats (to write)",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading number of external Stats to write in "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			// allocating memory for the external gridnames, entities and grid pointers, rename modes and fillmiss codes
			if (iNumExternalStatsWrite>0)
			{
				ExternalStatNamesWrite    = new char*[iNumExternalStatsWrite];
				ExternalStatEntitiesWrite = new char*[iNumExternalStatsWrite];
				ExternalStatsWriteCodes     = new char*[iNumExternalStatsWrite];
				for (int i=0;i<iNumExternalStatsWrite;i++)
				{
					ExternalStatNamesWrite[i] = new char[_stringsize_+1];
					ExternalStatEntitiesWrite[i] = new char[_stringsize_+1];
					ExternalStatsWriteCodes[i]  = new char[50];
					if (fscanf(steuerdatei,"%s%s%s",
								&(ExternalStatNamesWrite[i][0]),
								&(ExternalStatEntitiesWrite[i][0]),
								&(ExternalStatsWriteCodes[i][0]))!=3)
					{
						(*WSstream) << "FATAL: -----------------error------------\n"
						"error reading external (write) Statfile name/entity/writecode"<<instanzname<<"\n"
						"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
						throw(10);
					}
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
			
			if (fscanfln(steuerdatei,"%s",SemaphoreFileNameWrite,"semaphore file name for informing the external model ",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading semaphore file name for file to be written by wasim "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			if (fscanfln(steuerdatei,"%s",SemaphoreContent,"semaphore Content for informing the external model ",tstr)!=1)
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
				"error reading semaphore Content for file to be written by wasim "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
		}
	}
	delete[] tstr;
}

bool defaultini::IsCoupled()
{
	return (iDoExternalCoupling==1)?true:false;
}

char* defaultini::GetSemaphoreFileNameRead()
{
	return SemaphoreFileNameRead;
}

char* defaultini::GetSemaphoreFileNameWrite()
{
	return SemaphoreFileNameWrite;
}

char* defaultini::GetSemaphoreContent()
{
	return SemaphoreContent;
}

long defaultini::GetSemaphoreWaitInterval()
{
	return lSemaphoreWaitInterval;
}


char defaultini::GetCouplingMode()
{
	return cCouplingMode;
}

void defaultini::InitializeExternalCoupling(TLastResults *p_LastResults,int p_maxLastResults)
{
	// find adresses and of grids and statistics 
	ExternalStatLastResults = p_LastResults;
	ExternalStatsMaxLastResults = p_maxLastResults;
	int i=0;
	int index;

	// find addresses of the grids (for reading) which must already be created in gridinfo
	for (i=0;i<iNumExternalGridsRead;i++)
	{
		if (ist_vorhanden(ExternalGridEntitiesRead[i]))
		{
			index=getindex(ExternalGridEntitiesRead[i]);
			ExternalGridsRead[i]=gridinfo[index].pointer;
		}
		else
		{
			(*WSstream) << "FATAL: -----------------error------------\n"
			"error looking for internal grid to read from "<<ExternalGridNamesRead[i]<<", with internal name \""<<ExternalGridEntitiesRead[i]<<"\"\n"
			"possible reason: internal name of the grid is wrong (case sensitive!)"<<"[end-entry]";
			throw(1003);
		}
	}

	// find addresses of the grids (for writing) which, too, must already be created in gridinfo
	for (i=0;i<iNumExternalGridsWrite;i++)
	{
		if (ist_vorhanden(ExternalGridEntitiesWrite[i]))
		{
			index=getindex(ExternalGridEntitiesWrite[i]);
			ExternalGridsWrite[i]=gridinfo[index].pointer;
		}
		else
		{
			(*WSstream) << "FATAL: -----------------error------------\n"
			"error looking for internal grid to write to "<<ExternalGridNamesWrite[i]<<" with internal name \""<<ExternalGridEntitiesWrite[i]<<"\"\n"
			"possible reason: internal name of the grid is wrong (case sensitive!)"<<"[end-entry]";
			throw(1004);
		}
	}

	// find addresses of the zonal statistics in TLastResults holding the statistic mean value for each zone
	if (mpi_rank == 0)
	{
		for (i = 0; i < iNumExternalStatsWrite; i++)
		{
			std::string entity;
			std::string entityFromLastResults;
			entity = std::string(ExternalStatEntitiesWrite[i]);
			int k;
			int z = 0;
			bool found = false;
			while ((p_LastResults[z].wert != NULL) && (z < ExternalStatsMaxLastResults))
			{
				entityFromLastResults = std::string(p_LastResults[z].Name);
				// the names in the control file must contain "_" when more than one word is specified there. So here we replace all Spaces by "_"
				for (k = 0; k < entityFromLastResults.size(); k++)
				{
					if (entityFromLastResults[k] == ' ')
						entityFromLastResults[k] = '_';
				}
				if (entityFromLastResults.find(entity) != std::string::npos)
				{
					found = true;
					break; // z will contain the index to LastResults holding the required information

				}
				z++;
			}
			if (found)
			{
				ExternalStatLastResults[i] = p_LastResults[z];
			}
			else
			{
				(*WSstream) << "FATAL: -----------------error------------\n"
					"error looking for zonal statistic with internal name \"" << ExternalStatEntitiesWrite[i] << "\"\n"
					"possible reason: name of the statistic entity is wrong (case sensitive!)" << "[end-entry]";
				throw(1005);
			}
		}
	}

}

void defaultini::AccumulateExternalEntities()
{
	int i;
	// accumulate grids for external coupling
	for (i=0;i<iNumExternalGridsWrite;i++)
	{
		if ((strcmp(ExternalGridsWriteCodes[i],"Sum")==0) || (strcmp(ExternalGridsWriteCodes[i],"Mean")==0))
		{
			int m,n;
			int idx = getindex(ExternalGridEntitiesWrite[i]);
			for (m = gridinfo[idx].nrl; m <= gridinfo[idx].nrh; m++)
			{
				for (n = FirstLastCol[m].First; n <= FirstLastCol[m].Last; n++)
				{
					if ((zn[m][n] != nodata) && (mpi_rank == MPIREG[m][n]))
					{
						ExternalGridsWriteSum[i][m][n]+=ExternalGridsWrite[i][m][n];
					}
					else
					{
						ExternalGridsWriteSum[i][m][n]=nodata;
					}
				}
			}
		}
	}

	// accumulating statistic values for each zone
	for (i=0;i<iNumExternalStatsWrite;i++)
	{
		if ((strcmp(ExternalStatsWriteCodes[i],"Sum")==0) 
				|| (strncmp(ExternalGridsWriteCodes[i],"SumTotal",8)==0)
				|| (strcmp(ExternalStatsWriteCodes[i],"Mean")==0))
		{
			for (int k=1; k<ExternalStatLastResults[i].NumZones;k++)
			{
				ExternalStatLastResults[i].refvalue[k]+=ExternalStatLastResults[i].wert[k];
			}
		}
	}
	ExternalSyncCounter++;
}

void defaultini::WriteExternalResults(int actCouplingHour,int actCouplingDay, int actCouplingMonth, int actCouplingYear,
bool newHour,bool newDay,bool newMonth,bool newYear)
{

	int n,m;
	int i;
	// write grids for external coupling --> no sums or averages are calculated by now!
	ExternalSyncCounter++;
	for (i=0;i<iNumExternalGridsWrite;i++)
	{
		int idx = getindex(ExternalGridEntitiesWrite[i]);
		if ((strcmp(ExternalGridsWriteCodes[i], "Sum") == 0) || (strncmp(ExternalGridsWriteCodes[i], "SumTotal", 8) == 0))
		{
			for (m = gridinfo[idx].nrl; m <= gridinfo[idx].nrh; m++)
			{
				for (n = FirstLastCol[m].First; n <= FirstLastCol[m].Last; n++)
				{
					if ((zn[m][n] != nodata) && (mpi_rank == MPIREG[m][n]))
					{
						ExternalGridsWriteSum[i][m][n] += ExternalGridsWrite[i][m][n];
					}
					else
					{
						ExternalGridsWriteSum[i][m][n]=nodata;
					}
				}
			}
			// new:
			double** tempgrid;
			if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWriteSum[i]); else tempgrid = ExternalGridsWriteSum[i];
			if (mpi_rank == 0)
			{
				printgridfile(ExternalGridNamesWrite[i], tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
			// old: 
			// printgridfile(ExternalGridNamesWrite[i], ExternalGridsWriteSum[i], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);

			// writing an extra grid with the total sum for a given Period, here the last day
			if ((newDay==true)&&(strncmp(ExternalGridsWriteCodes[i],"SumTotal_D",10)==0))
			{
				char name[1024];
				std::string Sname=std::string(ExternalGridNamesWrite[i]);
				sprintf(name,"%sSumTotal_EndOf_Day_%.04d_%.02d_%.02d_%s",
				Sname.substr(0,Sname.find_last_of(dirzeichen)+1).c_str(),
				actCouplingYear,actCouplingMonth,actCouplingDay,
				Sname.substr(Sname.find_last_of(dirzeichen) + 1).c_str());

				// new:
				double** tempgrid;
				if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWriteSum[i]); else tempgrid = ExternalGridsWriteSum[i];
				if (mpi_rank == 0)
				{
					printgridfile(name, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
				if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
				// old: 
				// printgridfile(name, ExternalGridsWriteSum[i], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			// writing an extra grid with the total sum for the last month
			if ((newMonth==true)&&
					((strncmp(ExternalGridsWriteCodes[i],"SumTotal_M",10)==0) ||
						(strncmp(ExternalGridsWriteCodes[i],"SumTotal_DM",11)==0))
					)
			{
				char name[1024];
				std::string Sname=std::string(ExternalGridNamesWrite[i]);
				sprintf(name,"%sSumTotal_EndOf_Month_%.04d_%.02d_%s",
					Sname.substr(0, Sname.find_last_of(dirzeichen) + 1).c_str(),
				actCouplingYear,actCouplingMonth,
				Sname.substr(Sname.find_last_of(dirzeichen) + 1).c_str());
				
				// new:
				double** tempgrid;
				if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWriteSum[i]); else tempgrid = ExternalGridsWriteSum[i];
				if (mpi_rank == 0)
				{
					printgridfile(name, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
				if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
				// old: 
				// printgridfile(name, ExternalGridsWriteSum[i], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			// writing an extra grid with the total sum for the last year
			if ((newYear==true)&&
					((strncmp(ExternalGridsWriteCodes[i],"SumTotal_Y",12)==0) ||
						(strncmp(ExternalGridsWriteCodes[i],"SumTotal_MY",12)==0) ||
						(strncmp(ExternalGridsWriteCodes[i],"SumTotal_DMY",11)==0))
					)
			{
				char name[1024];
				std::string Sname=std::string(ExternalGridNamesWrite[i]);
				sprintf(name,"%sSumTotal_EndOf_Year_%.04d_%s",
				Sname.substr(0,Sname.find_last_of(dirzeichen)+1).c_str(),
				actCouplingYear,
				Sname.substr(Sname.find_last_of(dirzeichen) + 1).c_str());

				// new:
				double** tempgrid;
				if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWriteSum[i]); else tempgrid = ExternalGridsWriteSum[i];
				if (mpi_rank == 0)
				{
					printgridfile(name, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
				if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
				// old: printgridfile(name, ExternalGridsWriteSum[i], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
		}
		else if (strcmp(ExternalGridsWriteCodes[i],"Mean")==0)
		{
			for (m = gridinfo[idx].nrl; m <= gridinfo[idx].nrh; m++)
			{
				for (n = FirstLastCol[m].First; n <= FirstLastCol[m].Last; n++)
				{
					if ((zn[m][n] != nodata) && (mpi_rank == MPIREG[m][n]))
					{
						ExternalGridsWriteSum[i][m][n] += ExternalGridsWrite[i][m][n];
						ExternalGridsWriteSum[i][m][n]=ExternalGridsWriteSum[i][m][n]/ExternalSyncCounter;
					}
					else
					{
						ExternalGridsWriteSum[i][m][n]=nodata;
					}
				}
			}
			// new:
			double** tempgrid;
			if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWriteSum[i]); else tempgrid = ExternalGridsWriteSum[i];
			if (mpi_rank == 0)
			{
				printgridfile(ExternalGridNamesWrite[i], tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
			// old:
			// printgridfile(ExternalGridNamesWrite[i],ExternalGridsWriteSum[i],WriteASCIIGrids,WriteDoubleGrids, WriteNetcdfGrids);
		}
		else
		{
			// new:
			double** tempgrid;
			if (mpi_size>1) tempgrid = MPIGatherIntoNewGrid(idx, ExternalGridsWrite[i]); else tempgrid = ExternalGridsWrite[i];
			if (mpi_rank == 0)
			{
				printgridfile(ExternalGridNamesWrite[i], tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if (mpi_size>1) free_matrix(tempgrid, 1, nrows, 1, ncols);
			// old:
			// printgridfile(ExternalGridNamesWrite[i], ExternalGridsWrite[i], WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
		}
		// set all values back to 0, except writecode was SumTotal
		if ((strcmp(ExternalGridsWriteCodes[i],"Mean")==0)||(strcmp(ExternalGridsWriteCodes[i],"Sum")==0))
		{
			for (m = gridinfo[idx].nrl; m <= gridinfo[idx].nrh; m++)
			{
				for (n = FirstLastCol[m].First; n <= FirstLastCol[m].Last; n++)
				{
					if ((zn[m][n] != nodata) && (mpi_rank == MPIREG[m][n]))
					{
						ExternalGridsWriteSum[i][m][n] = 0;
					}
					else
					{
						ExternalGridsWriteSum[i][m][n]=nodata;
					}
				}
			}
		}
	}

	// gather all statistics now before 
	// writing statistic values for each zone --> no sums or averages are calculated by now!
	if (mpi_rank == 0)
	{
		for (i = 0; i < iNumExternalStatsWrite; i++)
		{
			std::ofstream StatFile;
			StatFile.open(ExternalStatNamesWrite[i], ios::trunc);

			if (std::string(ExternalStatLastResults[i].Name).size()>150)
			{
				StatFile 
					<< "date: " << setfill('0') << setw(4) << actCouplingYear << "-" << setfill('0') << setw(2) << actCouplingMonth << "-" << setfill('0') << setw(2) << actCouplingDay << "\n"
					<< "hour: " << setfill('0') << setw(2) << actCouplingHour << "h\n"
					<< std::string(ExternalStatLastResults[i].Name).substr(0, 150).c_str() << "..." << "; " << zonmax << " zones" << endl;
			}
			else
			{
				StatFile 
					<< "date: " << setfill('0') << setw(4) << actCouplingYear << "-" << setfill('0') << setw(2) << actCouplingMonth << "-" << setfill('0') << setw(2) << actCouplingDay << "\n"
					<< "hour: " << setfill('0') << setw(2) << actCouplingHour << "h\n" << ExternalStatLastResults[i].Name << "; " << zonmax << " zones" << endl;
			}

			int k = 0;
			// todo: write mean or sum when such a writecde was defined
			if (strcmp(ExternalStatsWriteCodes[i], "Sum") == 0)
			{
				double sum = 0;
				for (k = 1; k < ExternalStatLastResults[i].NumZones; k++)
				{
					// one last accumulation must be done
					ExternalStatLastResults[i].refvalue[k] += ExternalStatLastResults[i].wert[k];
					StatFile << ExternalStatLastResults[i].ZoneCodes[k] << "\t" << ExternalStatLastResults[i].refvalue[k] << endl;
					sum += ExternalStatLastResults[i].rel[k] * ExternalStatLastResults[i].refvalue[k];
					ExternalStatLastResults[i].refvalue[k] = 0;
				}
				StatFile << "average" << "\t" << sum << endl;
			}
			else if (strncmp(ExternalStatsWriteCodes[i], "SumTotal", 8) == 0)
			{
				double sum = 0;
				double sumrel = 0;
				for (k = 1; k < ExternalStatLastResults[i].NumZones; k++)
				{
					// one last accumulation must be done
					ExternalStatLastResults[i].refvalue[k] += ExternalStatLastResults[i].wert[k];
					StatFile << ExternalStatLastResults[i].ZoneCodes[k] << "\t" << ExternalStatLastResults[i].refvalue[k] << endl;
					sum += ExternalStatLastResults[i].rel[k] * ExternalStatLastResults[i].refvalue[k];
					sumrel += ExternalStatLastResults[i].rel[k];
					// ExternalStatLastResults[i].refvalue[k]=0; because SumTotal, no reset is done
				}
				StatFile << "average" << "\t" << sum << endl;
			}
			else if (strcmp(ExternalStatsWriteCodes[i], "Mean") == 0)
			{
				double sum = 0;
				for (k = 1; k < ExternalStatLastResults[i].NumZones; k++)
				{
					// one last accumulation must be done
					ExternalStatLastResults[i].refvalue[k] += ExternalStatLastResults[i].wert[k];
					StatFile << ExternalStatLastResults[i].ZoneCodes[k] << "\t" << ExternalStatLastResults[i].refvalue[k] / ExternalSyncCounter << endl;
					sum += ExternalStatLastResults[i].rel[k] * ExternalStatLastResults[i].refvalue[k];
					ExternalStatLastResults[i].refvalue[k] = 0;
				}
				StatFile << "average" << "\t" << sum / ExternalSyncCounter << endl;
			}
			else
			{
				double sum = 0;
				for (k = 1; k < ExternalStatLastResults[i].NumZones; k++)
				{
					StatFile << ExternalStatLastResults[i].ZoneCodes[k] << "\t" << ExternalStatLastResults[i].wert[k] << endl;
					sum += ExternalStatLastResults[i].rel[k] * ExternalStatLastResults[i].wert[k];
				}
				StatFile << "average" << "\t" << sum << endl;
			}
			StatFile.close();
		}
	}
	ExternalSyncCounter=0;
}

void defaultini::InitializeFromExternalGrids(long jahr,long monat,long tag,long stunde)
{

	// reading the grids from external models
	// the specific initialization for each submodel has to be done there (but may be called from outside, 
	// because the model itself does not need to know if it was coupled or not)
	int i;
	char *newname;
	newname = new char[_stringsize_];
	for (i=0;i<iNumExternalGridsRead;i++)
	{
		try
		{
			readgridfromfile(ExternalGridNamesRead[i], &(ExternalGridsRead[i]), iFillMiss[i], false, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
			SetSignalAbfangfunktion();
		}
		catch (exception &e)
		{
			(*WSstream) << "FATAL: error reading external grid "<<ExternalGridNamesRead[i]<<"[end-entry]";
			throw(1001);
		}
		if (mpi_rank == 0)
		{
			if (iRenameMode[i] == _RenameMode_Rename_)
			{
				sprintf(newname, "%s_%.4ld_%.2ld_%.2ld_%.2ld_%ld", ExternalGridNamesRead[i], jahr, monat, tag, stunde, clock());
				if (rename(ExternalGridNamesRead[i], newname) != 0)
				{
					(*WSstream) << "FATAL: error renaming external grid " << ExternalGridNamesRead[i] << " nach " << newname << "[end-entry]";
					throw(1002);
				}
			}
			else if (iRenameMode[i] == _RenameMode_Delete_)
			{
				if (remove(ExternalGridNamesRead[i]) != 0)
				{
					(*WSstream) << "FATAL: error deleting external grid " << ExternalGridNamesRead[i] << "[end-entry]";
					throw(1002);
				}
			}
		}
	}
	delete[] newname;
	if (mpi_rank == 0)
	{
		if (remove(SemaphoreFileNameRead) != 0) // first try
		{
			clock_t wait=100;
#ifdef WIN32
			_sleep(wait);
#elif WIN64
			_sleep(wait);
#else
			clock_t goal;
			goal = wait + clock();
			while (goal > clock());
#endif
			if (remove(SemaphoreFileNameRead) != 0) // second try
			{
				(*WSstream) << "FATAL: error deleting semaphore file " << SemaphoreFileNameRead << "[end-entry]";
				throw(1001);
			}
		}
	}
}

bool defaultini::ExistsExternalGWTable()
{
	bool retval=false;
	for (int z=0;z<iNumExternalGridsRead;z++)
	{
		if (strcmp(ExternalGridEntitiesRead[z],_GWTableExtern_)==0)
		{
			retval = true;
			break;
		}
	}
	return retval;
}

int defaultini::GetExternalTimeStep()
{
	return iExternalTimeStep;
}

void defaultini::MakeMPIRegions(long **v_regioncnt, std::vector<mpihalo> **sendHalo, std::vector<mpihalo> **recvHalo, int *p_mpistartrow, int *p_mpiendrow)
{
	int i;
	char *aktname;

	// default: mpistartrow is 1 and mpiendrow is nrows (for non-MPI applications)
	(*p_mpistartrow) = 1;
	(*p_mpiendrow) = nrows;
	int *locMpiStartRow;
	int *locMpiEndRow;

	locMpiStartRow = new int[mpi_size+1];
	locMpiEndRow = new int[mpi_size+1];

	if (mpi_rank == 0)
	{ 
		// MPI grids will now be created for rank 0 only since all the other instances will get their optimized grids later when we know how many rows are required for each rank 
		// (note: each MPI region grid will include the halo on each side, so the grids to reserve will be  abit larger then the MPI region itself (i.e. 2 additional rows) in order to allow the application 
		// to get the neighboring MPI ranks values from the halos
		MPIREG = NULL;
		MPIHLP = NULL;
		MPIHALO = NULL;
		aktname = newstring(_stringsize_, 1);
		mpiregindex = 1;
		while (gridinfo[mpiregindex].name != NULL)
		{
			mpiregindex++;
			if (mpiregindex == gridmaxanzahl)
			{
				(*WSstream) << "FATAL: when trying to create MPI-Regions: no free grid found in memory (at position " << mpiregindex << ")" << "[end-entry]";
				throw(12);
			}
		}

		initialize_grid(mpiregindex, aktname, NULL, NULL, 0, _MPIRegions_, 0,1,nrows);
		gridinfo[mpiregindex].pointer = matrix(1, nrows, 1, ncols,nrows);
		MPIREG = gridinfo[mpiregindex].pointer;
		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIREG[i][j] = nodata;
			}
		}
		gridanzahl++;

		mpihlpindex = mpiregindex + 1;
		initialize_grid(mpihlpindex, aktname, NULL, NULL, 0, _MPIHelper_, 0, 1, nrows);
		gridinfo[mpihlpindex].pointer = matrix(1, nrows, 1, ncols,nrows);
		MPIHLP = gridinfo[mpihlpindex].pointer;
		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIHLP[i][j] = 0;
			}
		}
		gridanzahl++;

		mpihaloindex = mpihlpindex + 1;
		initialize_grid(mpihaloindex, aktname, NULL, NULL, 0, _MPIHalo_, 0, 1, nrows);
		gridinfo[mpihaloindex].pointer = matrix(1, nrows, 1, ncols, nrows);
		MPIHALO = gridinfo[mpihaloindex].pointer;
		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIHALO[i][j] = nodata;
			}
		}
		gridanzahl++;

		free(aktname);
	}
	// MPI performance analysis for every node
	double starttime,endtime,wtime;
	double *v_wtime = new double[mpi_size+1];
	double *v_mpiweight = new double[mpi_size+1];
	*v_regioncnt = new long[mpi_size+1];
	double **m1,**m2,**m3;
#ifdef _MPI
	// this is an autmatic load balancer, doing a mathematic computation on every machine (node) beeing part of the MPI call
	// the number of cells allocated to each MPI regions is in direct proportion to the speed, this loop will be done
	// However, if a file with name MpiSpeedFactor.dat is found in the working directory, the speed test is bypassed and the value read from the file is taken instead, thus providing a way
	// to define the work loads for each node (assuming there are independent file systems, 
	// otherwise the workload will be evenly distributed across all nodes, as is usualy the best for supercomputers)
	std::fstream mpispeedfile;
	mpispeedfile.open("MpiSpeedFactor.dat",std::ios_base::in);
	if (!mpispeedfile.good())
	{
		int mpitest_nrows=1000;
		int mpitest_ncols=1000;

		m1 = matrix(1,mpitest_nrows,1,mpitest_ncols);
		m2 = matrix(1,mpitest_nrows,1,mpitest_ncols);
		m3 = matrix(1,mpitest_nrows,1,mpitest_ncols);

		wtime=0;

		for (int k=1;k<=30;k++)
		{
			starttime = MPI_Wtime();
			int i,j;
#pragma omp parallel for schedule(runtime) private(i,j)
			for (i=1;i<=mpitest_nrows;i++)
			{
				for (j=1;j<=mpitest_ncols;j++)
				{
					m1[i][j]=cos(double(i+j));
					m2[i][j]=sin(double(i+j));
					m3[i][j]=atan(m1[i][j]*m2[i][j]);
				}
			}
			endtime = MPI_Wtime();
			wtime += endtime-starttime;
		}

		free_matrix(m1,1,mpitest_nrows,1,mpitest_ncols);
		free_matrix(m2,1,mpitest_nrows,1,mpitest_ncols);
		free_matrix(m3,1,mpitest_nrows,1,mpitest_ncols);

		(*WSstream) << "MPI wtime: " << wtime << "for instance: " << mpi_rank << "[end-entry]";
	}
	else
	{
		mpispeedfile >> wtime;
		mpispeedfile.close();
	}
	// now sending all times to the master process. There, the fraction to be used for each node will be calculated
	if (mpi_rank==0)
	{
		v_wtime[0]=wtime;
		for (int i=1;i<mpi_size;i++)
		{
			MPI_Status s;
			MPI_Recv(&(v_wtime[i]),1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&s);
		}
	}
	else
	{ // all other ranks should send their results to instance 0
		MPI_Send(&(wtime),1,MPI_DOUBLE,0,mpi_rank, MPI_COMM_WORLD);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank==0)
	{
		for (int i=0;i<mpi_size;i++)
		{
			(*WSstream) << "received MPI wtime: " << v_wtime[i] << " for instance: " << i << "[end-entry]";
		}
	}
	// now weighing the factors
	double wtimesum=0;
	for (int i=0;i<mpi_size;i++)
	{
		if (v_wtime[i]>0)
		{
			v_wtime[i]=1.0/v_wtime[i];
		}
		else
		{
			v_wtime[i]=1.0;
		}
		wtimesum+=v_wtime[i];
	}
	// actual weight
	for (int i=0;i<mpi_size;i++)
	{
		v_mpiweight[i]=v_wtime[i]/wtimesum;
	}
	if (mpi_rank==0)
	{
		for (int i=0;i<mpi_size;i++)
		{
			(*WSstream) << "MPI performance weight: " << v_mpiweight[i] << " for instance: " << i << "[end-entry]";
		}
	}


#endif
	// creating regions for MPI, will be done in the master instance only and then send to the other instances
	if (mpi_rank==0)
	{
		// all cells below a river MUST be calculated by the same instance, so we will assign those cells to instance with rank 0
		// which will be available even with no MPI
		double **RIVLNK;
		int rivlnkindex=ist_vorhanden(_river_links_,true);

#ifdef _MPI
		// for release 10.03: switch on river link cell counting again
		// but only if interaction code is set to 2 (i.e. exfiltration and infiltration possible)
		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: Error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}

		if (findeabschnitt("[" _Unsatzonmodell_ "]", false))
		{
			int param=0;
			if (fscanf(steuerdatei, "%d", &param) != 1) steuererror("looking for re-infiltration parameter in [" _Unsatzonmodell_ "]", "");
			fgets(satz, _stringsize_, steuerdatei);
			if (fscanf(steuerdatei, "%d", &param) != 1) steuererror("looking for re-infiltration parameter in [" _Unsatzonmodell_ "]", "");
			fgets(satz, _stringsize_, steuerdatei);
			if (fscanf(steuerdatei, "%d", &param) != 1) steuererror("looking for re-infiltration parameter in [" _Unsatzonmodell_ "]", "");
			fgets(satz, _stringsize_, steuerdatei);
			if (fscanf(steuerdatei, "%d", &param) != 1) steuererror("looking for re-infiltration parameter in [" _Unsatzonmodell_ "]", "");
			if (param != 2)
			{
				rivlnkindex = 0;
			}
		}
		else
		{
			rivlnkindex = 0;
		}

		if ((fclose(steuerdatei)) == EOF)
		{
			(*WSstream) << "FATAL: error closing " << steuerdateiname << "[end-entry]";
			throw(11);
		}
		// when ignoring the riverlinks, it will NOT be possible to reinfiltrate water in another MPI-Region than 0 since because the river flow table will be avalable in rank 0 only (
		// and routing will be done in rank 0 only, anyway. The advantage is, that the application is much faster because only some messages (in stripes) have to be exchanged.
#endif
		if (rivlnkindex>0) RIVLNK=gridinfo[rivlnkindex].pointer;
		// now counting active cells and calculating the fraction for each instance
		long cellcnt=0;
		long rivcnt=0;
		
		// counting
		for (int i=1;i<=nrows;i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				if (zn[i][j]>nodata)
				{
					if ((rivlnkindex>0)&&(RIVLNK[i][j]>0))
					{
						rivcnt++;
						// adding each river linkk cell to the list of river link cells to be sorted and rearranged later
						RiverLinkCellList.push_back(RiverLinkCell(i, j, DHM[i][j], (int)RIVLNK[i][j]));
						
					}
					else
					{
						cellcnt++;
					}
				}
			}
		}
		RiverLinkCellList.sort();

		// constructing a list of river links, each holding a vector of river cells belonging to this link, ordered from highest to lowest elevation, taken from the sorted list above (the reaso we used a list for the cells was the easy sorting)
		// in case the list has been filled somewhere else... just to make sure we start with an empty list
		RiverLinkList.clear();
		std::list<RiverLinkCell>::iterator RLCL_iter = RiverLinkCellList.begin();
		
		while (RLCL_iter != RiverLinkCellList.end())
		{
			// a new river link starts here, initialize with the links ID and the processingrank, although this rank is still unknown
			RiverLinkList.push_back(new RiverLink(RLCL_iter->link, RLCL_iter->processingrank,0,NULL));
			int actlink = RLCL_iter->link;
			do
			{	// add all list elements from the RiverLinkCellList as new elements of a vector (later, the elements needs to be directly accessed for performance reasons)
				RiverLinkList.back()->addRiverLinkCell(*RLCL_iter); // this will also increment the cell counter for each RiverLinkList-element, so at the end of this process we know the number of cells for each link not only from the vectors size but also from the element cellcnt of a RiverLinkList element
				RLCL_iter++;
			} while (RLCL_iter != RiverLinkCellList.end() && (RLCL_iter->link == actlink));
		}
		// using a so called lambda to sort the list
		RiverLinkList.sort([](const RiverLink* f, const RiverLink* s) {return f->cellcnt > s->cellcnt; });
		// this call will not work since we have pointers in the list and an overloaded operator< will not sort the objects after any internal members but by address		
		//		RiverLinkList.sort();

		// now converting the list into a vector for random access
		std::list<RiverLink*>::iterator RLL_iter = RiverLinkList.begin();
		while (RLL_iter != RiverLinkList.end())
		{
			RiverLinkVector.push_back(*RLL_iter);
			RLL_iter++;
		}
		// deleting the riverLinkList and the RiverLinkCellList now (it's all part of the RiverLinkVector now)
		RiverLinkCellList.clear(); // we do not need the RiverLinkCellList anymore
		RiverLinkList.clear();

		// ToDo: distribute links to ranks 0...n until all links are covered somehow. If the number of free rows doesn't suffice, we have still a problem...

		long CellsPerRank = (long)floor((cellcnt + rivcnt + 0.4999)/mpi_size);
		int actRiverLinkVectorIndex = 0;

		(*WSstream) << "MPI regular cells: " << cellcnt << "\nMPI river cells: " << rivcnt << "[end-entry]";
		// setting rank codes
		long totalcellcnt = 0;
		long totalrivercells = 0;
		for (int i=0;i<mpi_size;i++)
		{
			if (i < mpi_size - 1)
			{
				// first guess for cell count for the actual region depending on mpiweight
				(*v_regioncnt)[i] = (long)floor(v_mpiweight[i] * (cellcnt + rivcnt) + 0.4999);
			}
			else
			{
				// last rank gets all rounding error cells
				(*v_regioncnt)[i] = cellcnt + rivcnt - totalcellcnt - totalrivercells;
			}
			// now subtracting the number of link cells that will be processed by this rank (starting with 0 usually means that most link cells will be processes in that rank)
			// by the way: the number of link cells will be subtracted from the active cells of all ranks, and when assigning an MPI rank to all cells, 
			// the rive rlink cells will also be assigned the actual rank but the counter will not be incremented, so in the end the counter was not incremented as often as there are river link cells, 
			// and on the other hand, that is the same number of cells that will be additionally processed by ranks 0 (and sometimes 1, 2 etc.) in addition to their reduced number of cells.
			// the MPI region grid will look strange (with no cells at all or very few cells for rank 0 sometimes). To avoid numnerical problems, the minimum fraction of cells should still be > 0
			// so there is a problem when the longest link has more cells that the average MPI rank will have cells assigned to. That will cause a warning message (user should add additional 
			// pour points to the long link in order to get shorter links, otherwise the processing of the link cells will slow down all ranks.
			bool stop = (RiverLinkVector.size() <= actRiverLinkVectorIndex); // if no more river links are available, we should right here set stop to true
			do
			{
				if (RiverLinkVector.size() > actRiverLinkVectorIndex)
				{
					(*v_regioncnt)[i] -= RiverLinkVector[actRiverLinkVectorIndex]->cellcnt;
					totalrivercells += RiverLinkVector[actRiverLinkVectorIndex]->cellcnt;
					RiverLinkVector[actRiverLinkVectorIndex]->linkProcessingrank = i;
					
					for (std::vector<RiverLinkCell>::iterator iter = RiverLinkVector[actRiverLinkVectorIndex]->linkCellVector->begin();	iter != RiverLinkVector[actRiverLinkVectorIndex]->linkCellVector->end(); iter++)
					{
						iter->processingrank = i;
					}
					if ((*v_regioncnt)[i] <= 0)
					{
						(*v_regioncnt)[i] = 1;
					}
					actRiverLinkVectorIndex++;
					if (actRiverLinkVectorIndex < RiverLinkVector.size())
					{
						if (((*v_regioncnt)[i] - RiverLinkVector[actRiverLinkVectorIndex]->cellcnt) <= 0) // one original cell should at least be processed by  every rank
						{
							stop = true;
						}
					}
					else
					{
						stop = true;
					}
				}
			} while (stop == false);  

			totalcellcnt += (*v_regioncnt)[i];
			int cnt = (*v_regioncnt)[i];
			(*WSstream) << "MPI num cells: " << (*v_regioncnt)[i] << " for instance: " << i << "[end-entry]";
			for (int k = 0; k < actRiverLinkVectorIndex; k++)
			{
				if (RiverLinkVector[k]->linkProcessingrank == i)
				{
					(*WSstream) << " + link cells: " << RiverLinkVector[k]->cellcnt << " for river link " << RiverLinkVector[k]->linkID << "[end-entry]";
					cnt += RiverLinkVector[k]->cellcnt;
				}
			}
			if (cnt > (*v_regioncnt)[i])
			{
				(*WSstream) << "        total: " << cnt << " for instance " << i << "[end-entry]";
			}
		}

		int actinstance=mpi_size-1;
		for (int i = 0; i < mpi_size; i++)
		{
			locMpiStartRow[i] = -1;
			locMpiEndRow[i] = -1;
		}

		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if (zn[i][j]>nodata)
				{
					if (locMpiStartRow[actinstance] < 0)
					{
						// this must be the first valid cell after initializing the start/end row arrays, so we can set the actual row for each
						locMpiStartRow[actinstance] = i;
						locMpiEndRow[actinstance] = i;
					}
// new for re-ifiltration in MPI: river link cells will be skipped at normal computation and then collected, processed in order in a single rank and 
// then scattered back to their residing rank, so the MPI region assignment is much more simple here
//					if ((rivlnkindex>0)&& (RIVLNK[i][j]>0))
//					{
//						MPIREG[i][j]=0;
//						rivcnt--;
//						if (locMpiStartRow[0] <= 0)
//						{
//							locMpiStartRow[0] = i;
//							locMpiEndRow[0] = i;
//						}
//					}
//					else
//					{
//						MPIREG[i][j] = actinstance;
//					}
					// new: river link cells will get the same MPI region but they don't count
					MPIREG[i][j] = actinstance;
					if ((rivlnkindex > 0) && (RIVLNK[i][j] != nodata))
					{
						// do nothing for now, later we could try to populate some fields in the RiverLinkVector, like the actual Region for each link cell
					}
					else
					{
						(*v_regioncnt)[actinstance]--;
					}
					if ((*v_regioncnt)[actinstance]==0)
					{
						locMpiEndRow[actinstance] = i;
						actinstance--;
						if (actinstance < 0)
						{
							actinstance = 0;
						}
					}
				}
			}
		}

//		std::vector<RiverLinkCell> *linkCellVector; // sub list for this link only, ordered from top to bottom, i.e. by elevation
//		std::vector<RiverLink*> RiverLinkVector;

		// now assigning the MPI rank to the river link vectors member linkCellVector
		for (int lnk = 0; lnk < RiverLinkVector.size(); lnk++)
		{
			for (int cell = 0; cell < RiverLinkVector[lnk]->linkCellVector->size(); cell++)
			{
				int i = RiverLinkVector[lnk]->linkCellVector->at(cell).i;
				int j = RiverLinkVector[lnk]->linkCellVector->at(cell).j;
				RiverLinkVector[lnk]->linkCellVector->at(cell).mpirank = MPIREG[i][j];
				RiverLinkVector[lnk]->linkCellVector->at(cell).processingorder = cell;
			}
		}


#ifdef _MPIDEBUG
//		printgridfile("MPIRegions_0",MPIREG,true,false);	
#endif
	}

	// Now creating a structure that holds river link cells per rank. One structure per remote rank and one per processing rank, 
	// so each rank does know where to send the cells to or from where to receive the cells 
	// buffers need to be created as well, but not within these structures (only untyped buffers will be sent by MPI, it's faster
	
#ifdef _MPI
	int numLinks;
	if (mpi_rank == 0)	
	{		
		numLinks = RiverLinkVector.size();	
	}
	MPI_Bcast(&numLinks, 1, MPI_INT, 0, MPI_COMM_WORLD);
	(*WSstream) << "After broadcasting numLinks: " << numLinks << " for instance: " << mpi_rank << "[end-entry]";
	for (int i = 0; i < numLinks; i++)
	{
		int vectorsize;
		if (mpi_rank == 0)
		{
			vectorsize = RiverLinkVector[i]->linkCellVector->size();
		}
		MPI_Bcast(&vectorsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		int elementsize = 6 * sizeof(int) + sizeof(double);

		// now all ranks now how many entries the actual vector will have
		char *buffer = new char[elementsize*vectorsize];
		if (mpi_rank == 0)
		{	// coding the vectors elements into the buffer
			// i.e. serializing RiverLinkVector
			for (int k = 0; k < vectorsize; k++)
			{
				memmove(&(buffer[k*elementsize + 0 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).i), sizeof(int));
				memmove(&(buffer[k*elementsize + 1 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).j), sizeof(int));
				memmove(&(buffer[k*elementsize + 2 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).link), sizeof(int));
				memmove(&(buffer[k*elementsize + 3 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).mpirank), sizeof(int));
				memmove(&(buffer[k*elementsize + 4 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).processingrank), sizeof(int));
				memmove(&(buffer[k*elementsize + 5 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).processingorder), sizeof(int));
				memmove(&(buffer[k*elementsize + 6 * sizeof(int)]), &(RiverLinkVector[i]->linkCellVector->at(k).z), sizeof(double));
			}
		}
		// broadcasting the actual serialized river link to all instances
		MPI_Bcast(buffer, elementsize*vectorsize, MPI_CHAR, 0, MPI_COMM_WORLD);
		if (mpi_rank != 0)
		{
			// create a new riverlink in the RiverLinkVector in all istances except 0
			RiverLinkVector.push_back(new RiverLink());
			for (int k = 0; k < vectorsize; k++)
			{
				RiverLinkVector.back()->linkCellVector->push_back(RiverLinkCell()); // create a new empty element, we know the number of elemenst that should follow from the broadcasted variable vectorsize
				// deserializing the buffer als replacing the defaultvalues of the new link cell element 
				memmove(&(RiverLinkVector.back()->linkCellVector->back().i					),&(buffer[k*elementsize + 0 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().j					),&(buffer[k*elementsize + 1 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().link				),&(buffer[k*elementsize + 2 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().mpirank			),&(buffer[k*elementsize + 3 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().processingrank		),&(buffer[k*elementsize + 4 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().processingorder	),&(buffer[k*elementsize + 5 * sizeof(int)]), sizeof(int));
				memmove(&(RiverLinkVector.back()->linkCellVector->back().z					),&(buffer[k*elementsize + 6 * sizeof(int)]), sizeof(double));
				if (k == 0)
				{
					// the riverlink itself needs to be initialized as well with a few common values that are also part of each river link cell
					RiverLinkVector.back()->cellcnt = vectorsize;
					RiverLinkVector.back()->linkID = RiverLinkVector.back()->linkCellVector->back().link;
					RiverLinkVector.back()->linkProcessingrank = RiverLinkVector.back()->linkCellVector->back().processingrank;
				}
			}
		}
		
	}

	MPI_Barrier(MPI_COMM_WORLD);

#ifdef _MPIDEBUG
	// only for debugging: is the received river link vector the same as the sent one?  
	if (mpi_rank == 1)
	{
		(*WSstream) << "received and decoded RiverLinkVector in rank "<<mpi_rank<<" has " << (int)RiverLinkVector.size() << " links." << "[end-entry]";
		for (int i = 0; i < RiverLinkVector.size(); i++)
		{
			(*WSstream) << " RiverLinkVector[" << i << "] in rank 1: cellcount is " << RiverLinkVector[i]->cellcnt << " (vector size " << (int)RiverLinkVector[i]->linkCellVector->size() << ")" << "[end-entry]";
			(*WSstream) << " RiverLinkVector[" << i << "] in rank 1: linkID is " << RiverLinkVector[i]->linkID << "[end-entry]";
			(*WSstream) << " RiverLinkVector[" << i << "] in rank 1: linkProcesingRank is " << RiverLinkVector[i]->linkProcessingrank << "[end-entry]";
			for (int k = 0; k < RiverLinkVector[i]->linkCellVector->size(); k++)
			{
				(*WSstream) << " RiverLinkVector[" << i << "] in rank 1: linkCellVector["<<k<<"].processingorder  " << (int)RiverLinkVector[i]->linkCellVector->at(k).processingorder << "[end-entry]";
			}
			(*WSstream) << "------------------ end of actual link "<<i<<"-------------------------------" << "[end-entry]";
		}
	}
#endif // _MPIDEBUG

#endif // _MPI

	// Step 3: creating localRiverLinkVector with sorted cells of the actual rank only for each processing rank
	// can have another name, like, e.g., 
	for (std::vector<RiverLink*>::iterator iter = RiverLinkVector.begin(); iter != RiverLinkVector.end(); iter++)
	{
		for (std::vector<RiverLinkCell>::iterator celliter = (*iter)->linkCellVector->begin(); celliter != (*iter)->linkCellVector->end(); celliter++)
		{
			if (celliter->mpirank == mpi_rank)
			{
				while ((signed int)(*iter)->linkProcessingrank > (signed int)RivLnkCellsToSend.size() - 1) // add a new remote rank if it doesn't exist yet
				{
					RivLnkCellsToSend.push_back(std::vector<RiverLinkCell>());
				}
				// add the actrual celliter cell to the RivLnkCellsToSend vector
				RivLnkCellsToSend.back().push_back(RiverLinkCell(*celliter));
			}
		}
	}

	// Step 4: creating another lcal RiverLinkVector for each processing rank as well
	// first creating empty vector elements
	for (int i = 0; i < mpi_size; i++)
	{
		RivLnkCellsToReceive.push_back(std::vector<RiverLinkCell>()); // one vector for each rank
	}
	// fill the receive vector. Later, this vector will be processed in every time step to receive the remote link cells data in a bufer and to send results back
	// each remote rank will have the same cells in its RivLnkCellsToSend for the processing rank as the processing rank expects cells from this remote rank. THis is made sure 
	// by using the same RiverLinkTable to create the vectors RivLnkCellsToSend and RivLnkCellsToReceive. So even the same instance can send and receive cells to itself (but it may be faster to copy it simply within the same instance)
	for (std::vector<RiverLink*>::iterator iter = RiverLinkVector.begin(); iter != RiverLinkVector.end(); iter++)
	{
		for (std::vector<RiverLinkCell>::iterator celliter = (*iter)->linkCellVector->begin(); celliter != (*iter)->linkCellVector->end(); celliter++)
		{
			if (celliter->processingrank == mpi_rank) // the cell in celliter should be processed in the actual instance, i.e. processingrank equals the mpi_rank of the actual instance
			{
				RivLnkCellsToReceive[celliter->mpirank].push_back(RiverLinkCell(*celliter));
			}
		}

	}


#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	// Step 5: exchanging data as a test several times: each cell does need e.g. 4000 byte of data, so there could be several MBytes to be sent and received in each time step

	int sr = 1;
	int er = nrows;

	// send the regions to the other instances
#ifdef _MPI
	MPI_Bcast(&(locMpiStartRow[0]), mpi_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(locMpiEndRow[0]), mpi_size, MPI_INT, 0, MPI_COMM_WORLD);

	(*p_mpistartrow) = locMpiStartRow[mpi_rank];
	(*p_mpiendrow) = locMpiEndRow[mpi_rank];

	sr = ((*p_mpistartrow) > 1) ? (*p_mpistartrow) - 1 : (*p_mpistartrow);
	er = ((*p_mpiendrow) < nrows) ? (*p_mpiendrow) + 1 : (*p_mpiendrow);

	if (mpi_rank == 0)
	{

		for (int remoterank = 1; remoterank < mpi_size; remoterank++)
		{
			int locsr = (locMpiStartRow[remoterank] > 1) ? locMpiStartRow[remoterank] - 1 : locMpiStartRow[remoterank];
			int locer = (locMpiEndRow[remoterank] < nrows) ? locMpiEndRow[remoterank] + 1 : locMpiEndRow[remoterank];

			for (int remoterow = locsr; remoterow <= locer; remoterow++)
			{
				// if (remoterow == locsr)
				// cout << "sending regions from startrow " << locsr << " to rank " << remoterank << " (global sr = " << sr << ")" << endl;
				// if (remoterow == locer)
				// cout << "sending regions till endrow " << locer << " to rank " << remoterank << " (global er = " << er << ")" << endl;
				MPI_Send(&(MPIREG[remoterow][1]), ncols, MPI_DOUBLE, remoterank, remoterow, MPI_COMM_WORLD);
			}
		}
	}
	else
	{
		// reserving optimized space for each rank
		MPIREG = NULL;
		MPIHLP = NULL;
		MPIHALO = NULL;
		aktname = newstring(_stringsize_, 1);
		mpiregindex = 1;
		while (gridinfo[mpiregindex].name != NULL)
		{
			mpiregindex++;
			if (mpiregindex == gridmaxanzahl)
			{
				(*WSstream) << "FATAL: when trying to create MPI-Regions: no free grid found in memory (at position " << mpiregindex << ")" << "[end-entry]";
				throw(12);
			}
		}

//		cout << "creating grids, startrow " << sr << " in rank " << mpi_rank << endl;
//		cout << "creating grids, endrow " << er << " in rank " << mpi_rank << endl;

		initialize_grid(mpiregindex, aktname, NULL, NULL, 0, _MPIRegions_, 0,sr,er);
		gridinfo[mpiregindex].pointer = matrix(sr, er, 1, ncols, nrows);
		MPIREG = gridinfo[mpiregindex].pointer;
		for (int i = sr; i<= er; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIREG[i][j] = nodata;
			}
		}
		gridanzahl++;

		mpihlpindex = mpiregindex + 1;
		initialize_grid(mpihlpindex, aktname, NULL, NULL, 0, _MPIHelper_, 0, sr, er);
		gridinfo[mpihlpindex].pointer = matrix(sr, er, 1, ncols, nrows);
		MPIHLP = gridinfo[mpihlpindex].pointer;
		for (int i = sr; i <= er; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIHLP[i][j] = 0;
			}
		}
		gridanzahl++;

		mpihaloindex = mpihlpindex + 1;
		initialize_grid(mpihaloindex, aktname, NULL, NULL, 0, _MPIHalo_, 0, sr, er);
		gridinfo[mpihaloindex].pointer = matrix(sr, er, 1, ncols, nrows);
		MPIHALO = gridinfo[mpihaloindex].pointer;
		for (int i = sr; i <= er; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				MPIHALO[i][j] = nodata;
			}
		}
		gridanzahl++;

		// receiving the relevant rows
		for (int localrow = sr; localrow <= er; localrow++)
		{
			//if (localrow == sr)
			//	cout << "receiving regions from startrow " << sr << " in rank " << mpi_rank << endl;
			//if (localrow == er)
			//	cout << "receiving regions till endrow " << er << " in rank " << mpi_rank << endl;
			MPI_Status s;
			MPI_Recv(&(MPIREG[localrow][1]), ncols, MPI_DOUBLE, 0, localrow, MPI_COMM_WORLD,&s);
		}

	}

#endif

	//filling the helper grid with informations on how long each strip in each row is (for send and recv-commands)
	// attention: it's not possible to simply take (FirstLastCol[i].Last-FirstLastCol[i].First+1) as runlength
	// because some cells will belong to rank 0 (river cells of a routing link when using re-infiltration will be handled by rank 0 only)
	if (mpi_rank>0)
	{
		for (int i = (*p_mpistartrow); i <= (*p_mpiendrow); i++)
		{
			int cnt=0;
			double lastj = nodata;
			for (int j = 1; j <= ncols; j++)
			{
				if (MPIREG[i][j] == mpi_rank)
				{
					if (MPIREG[i][j] == lastj)
					{
						MPIHLP[i][2 * (cnt - 1) + 2]++;
					}
					else
					{
						cnt++;
						MPIHLP[i][2 * (cnt - 1) + 1] = j;
						MPIHLP[i][2 * (cnt - 1) + 2]++;
					}
				}
				lastj = MPIREG[i][j];
			}
		}
	}
	else
	{ // for the master instance, ALL strips must be stored because this instance will receive all messages
		for (int i=1;i<=nrows;i++)
		{
			int cnt=0;
			double lastj = nodata;
			for (int j = 1; j <= ncols; j++)
			{
				if (MPIREG[i][j] == lastj)
				{
					MPIHLP[i][2 * (cnt - 1) + 2]++;
				}
				else
				{
					cnt++;
					MPIHLP[i][2 * (cnt - 1) + 1] = MPIREG[i][j];
					MPIHLP[i][2 * (cnt - 1) + 2]++;
				}
				lastj = MPIREG[i][j];
			}
		}
	}

	// creating MPIHALO 
	int startr, endr;
	if (mpi_rank == 0)
	{
		startr = 1;
		endr = nrows;
	}
	else
	{
		startr = (*p_mpistartrow);
		endr = (*p_mpiendrow);
	}

	for (int i = startr; i <= endr; i++)
	{
		for (int j = 1; j <= ncols; j++)
		{
			if (MPIREG[i][j] == mpi_rank)
			{
				// we set the actual and the real halo cell to the value of the other zone. This is because each rank will receive the other ranks values
				// at the positions where there is the rank-value of the other instance on a cell with the actual rank.
				// on the other hand, each cell on a halo where only the other instances rank is there will be cells which vaue must be send to that other ranks instance
				// 2016-06-19: only the Halo cells to be received from each surrounding region can be assigned to a singe grid (because the actual cell may be in the halo of multiple other regions, there is no single send-Halo grid possible)
				// and the send-Halo vector will be filled by sending the halo vector elements from instance to instance later
				if ((i >     1) && (j > 1)     && (MPIREG[i - 1][j - 1] != MPIREG[i][j]) && (MPIREG[i - 1][j - 1] != nodata))	{ MPIHALO[i - 1][j - 1] = MPIREG[i - 1][j - 1]; }
				if ((i >     1) &&                (MPIREG[i - 1][j    ] != MPIREG[i][j]) && (MPIREG[i - 1][j    ] != nodata))	{ MPIHALO[i - 1][j    ] = MPIREG[i - 1][j    ]; }
				if ((i >     1) && (j < ncols) && (MPIREG[i - 1][j + 1] != MPIREG[i][j]) && (MPIREG[i - 1][j + 1] != nodata))	{ MPIHALO[i - 1][j + 1] = MPIREG[i - 1][j + 1]; }
				if (               (j < ncols) && (MPIREG[i    ][j + 1] != MPIREG[i][j]) && (MPIREG[i    ][j + 1] != nodata))	{ MPIHALO[i    ][j + 1] = MPIREG[i    ][j + 1]; }
				if ((i < nrows) && (j < ncols) && (MPIREG[i + 1][j + 1] != MPIREG[i][j]) && (MPIREG[i + 1][j + 1] != nodata))	{ MPIHALO[i + 1][j + 1] = MPIREG[i + 1][j + 1]; }
				if ((i < nrows) &&                (MPIREG[i + 1][j    ] != MPIREG[i][j]) && (MPIREG[i + 1][j    ] != nodata))	{ MPIHALO[i + 1][j    ] = MPIREG[i + 1][j    ]; }
				if ((i < nrows) && (j > 1)     && (MPIREG[i + 1][j - 1] != MPIREG[i][j]) && (MPIREG[i + 1][j - 1] != nodata))	{ MPIHALO[i + 1][j - 1] = MPIREG[i + 1][j - 1]; }
				if (               (j > 1)     && (MPIREG[i    ][j - 1] != MPIREG[i][j]) && (MPIREG[i    ][j - 1] != nodata))	{ MPIHALO[i    ][j - 1] = MPIREG[i    ][j - 1]; }
			}
		}
	}

	// updating the array containing the number of cells per MPI region (will be used later for dimensioning the MPI send and receive buffers
	for (int i=0;i<mpi_size;i++)
	{
		(*v_regioncnt)[i]=0;
	}

	if (mpi_rank == 0)
	{
		for (int i = 1; i <= nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[i][j]>nodata) && (MPIREG[i][j] > nodata))
				{
					// cell must belong to a MPI region, so increment the respective counter
					(*v_regioncnt)[(int)MPIREG[i][j]]++;
				}
			}
		}
	}


#ifdef _MPI

	MPI_Bcast(&((*v_regioncnt)[0]), mpi_size, MPI_LONG, 0, MPI_COMM_WORLD);

#endif

	// mpihalo = struct (int row,col,cnt;};
	*sendHalo = new vector<mpihalo>();
	*recvHalo = new vector<mpihalo>();
	
	map<int, int> haloMapR; // map that will be used for storing the number of halo sections we will need to exchange for each unique rank ID with the remote ranks. 
    // Key: rank-ID of the remote rank; Value: number of vector elements in recvHalo
	map<int, int> haloMapS;

	double actRegion=-1;
	double actHalo=-1;
#ifdef _MPI
	// creating elements in the sendHalo and recvHalo vectors
	for (int i = sr; i <= er; i++)
	{ // an element in sendHalo or recvHalo encompasses at maximum one row
		// if the halo lies within the actual region, the cells value must be used for sending to another region. 
		// If the halo is outside the actual region, the cells value must be used for receiving from another region
		actRegion=-1;
		actHalo=-1;
		for (int j = 1; j <= ncols; j++)
		{
			if (MPIHALO[i][j]>=0) // we have found the beginning of a halo
			{
				if ((MPIREG[i][j]==actRegion)&&(MPIHALO[i][j]==actHalo))
				{
					(*recvHalo)->back().cnt++;
				}
				else
				{
					// first check if there is a halo with the same region code for this grid row already. 
					// If so, connect the two halo sections
					if (((*recvHalo)->size() > 0) && ((*recvHalo)->back().rank == MPIHALO[i][j]) && ((*recvHalo)->back().row == i))
					{
						(*recvHalo)->back().cnt = j - (*recvHalo)->back().col + 1;
					}
					else
					{
						mpihalo temp = { i, j, 1, int(MPIHALO[i][j]) };
						(*recvHalo)->push_back(temp);
						if (haloMapR.count(temp.rank) == 0)
						{
							haloMapR.insert(pair<int, int>(temp.rank, 1));
						}
						else
						{
							haloMapR.find(temp.rank)->second++;
						}
					}
				}
				actRegion=MPIREG[i][j];
				actHalo=MPIHALO[i][j];
			}
			else
			{
				actRegion=-1;
				actHalo=-1;
			}
		}
	}

	// send to each neighbor the count of sendHalo elements that this neighbor should send (because this actual rank expects as many elements)
	// so that the next operation is to exchange the dimensions of these halo elements
	MPI_Request *requestVector;
	requestVector = (MPI_Request*)calloc(2 * haloMapR.size(), sizeof(MPI_Request));

	typedef int mpihalomaptype[2];

	mpihalomaptype *haloMapFlatR = new mpihalomaptype[haloMapR.size()];
	mpihalomaptype *haloMapFlatS = new mpihalomaptype[haloMapR.size()];
	int halocntR = 0;
	int mapcnt = 0;
	// content of haloMap will be put into a buffer haloMapFlatS that will then used for sending messages to the receiving rank and vice versa with haloMapFlatR (initialized to 0)
	for (std::map<int,int>::iterator iter = haloMapR.begin(); iter != haloMapR.end(); iter++)
	{
		haloMapFlatR[mapcnt][0] = iter->first;
		haloMapFlatR[mapcnt][1] = iter->second;
		halocntR += iter->second;

		haloMapFlatS[mapcnt][0] = 0;
		haloMapFlatS[mapcnt][1] = 0;

		mapcnt++;
	}

	int request_cnt = 0;
	mapcnt = 0;
	for (std::map<int,int>::iterator iter = haloMapR.begin(); iter != haloMapR.end(); iter++)
	{
		MPI_Isend(&(haloMapFlatR[mapcnt][0]), 2, MPI_INT, haloMapFlatR[mapcnt][0], 0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
		request_cnt++;
		MPI_Irecv(&(haloMapFlatS[mapcnt][0]), 2, MPI_INT, haloMapFlatR[mapcnt][0], 0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
		request_cnt++;
		mapcnt++;
	}
	MPI_Waitall(request_cnt, requestVector, statusVector);
	free(requestVector);
	
	int halocntS = 0;
	mapcnt = 0;
	for (std::map<int, int>::iterator iter = haloMapR.begin(); iter != haloMapR.end(); iter++)
	{
		// when crating the haloMapS, we set the rank not from the received elements (since that value will always be the actual rank) but tot the rank of the corresponding haloMapR
		haloMapS.insert(pair<int, int>(haloMapFlatR[mapcnt][0], haloMapFlatS[mapcnt][1]));
		halocntS += haloMapFlatS[mapcnt][1];
		mapcnt++;
	}

	// now that the recv-vector is complete we know, that each instance we are receiving a halo from should also receive a halo from this instance, 
	// so we can construct the send-vector be receiving messgaes from the other instances 
	// (and send messages to the other instances that expect the halo from this instance)

	requestVector = (MPI_Request*)calloc(halocntR+halocntS, sizeof(MPI_Request));
	request_cnt = 0;

	int lastrank = -1;

	typedef int mpihalovectortype[4];

	mpihalovectortype *vHaloS = new mpihalovectortype[halocntS];
	mpihalovectortype *vHaloR = new mpihalovectortype[halocntR];
	
	// copying Halo metadata for halos to be received into a buffer that will be used for sending the requirements to the remote rank
	int totalcount = 0;
	for (std::vector<mpihalo>::iterator iter = (*recvHalo)->begin(); iter < (*recvHalo)->end(); iter++)
	{
		vHaloR[totalcount][0] = iter->row;
		vHaloR[totalcount][1] = iter->col;
		vHaloR[totalcount][2] = iter->cnt;
		vHaloR[totalcount][3] = mpi_rank;
		totalcount++;
	}

	// now sending the halo metadata information (on halos that this rank requires from the remote rank) to that rank.
	// the remote rank knows about the number of requests because we exchanged the haloMaps before where for all remote ranks the expected count of halos was listed
	totalcount = 0;
	for (std::vector<mpihalo>::iterator iter = (*recvHalo)->begin(); iter < (*recvHalo)->end(); iter++)
	{
		MPI_Isend(&(vHaloR[totalcount][0]), 4, MPI_INT, iter->rank, 0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
		request_cnt++;
		totalcount++;
	}

	// now receiving the halo metadata for the halos the remote rank expects from this rank
	totalcount = 0;
	for (std::map<int,int>::iterator iter = haloMapS.begin(); iter != haloMapS.end(); iter++)
	{
		for (int i = 1; i <= iter->second; i++)
		{
			MPI_Irecv(&(vHaloS[totalcount][0]), 4, MPI_INT, iter->first,0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
			request_cnt++;
			totalcount++;
		}
	}
	MPI_Waitall(request_cnt, requestVector,statusVector);
	free(requestVector);

	// putting the sendHalo elements into the vector that is used for halo exchange later
	for (int i = 0; i < totalcount;i++)
	{
		mpihalo temp;
		temp.row = vHaloS[i][0];
		temp.col = vHaloS[i][1];
		temp.cnt = vHaloS[i][2];
		temp.rank = vHaloS[i][3];
	
		(*sendHalo)->push_back(temp);
	}

	delete[] haloMapFlatR;
	delete[] haloMapFlatS;
	delete[] vHaloR;
	delete[] vHaloS;

//	for debugging HALO creation	
#ifdef _MPIDEBUG
	std::fstream ofs;
	std::stringstream ss;
	ss << "halo"<<mpi_rank<<".txt";
	std::string name;
	ss >> name;
	ofs.open(name.c_str(),std::ios::out);
	ofs << "DEBUG: " << "sendHalo-vector: "<<int((*sendHalo)->size())<<" elements.\n";
	for (std::vector<mpihalo>::iterator iter=(*sendHalo)->begin();iter<(*sendHalo)->end();iter++)
	{
		ofs << "row: " << int(iter->row) << " col: " << int(iter->col) << " cnt: " << int(iter->cnt) << " rank: " << int(iter->rank) << "\n";
	}

	ofs << "\nDEBUG: " << "recvHalo-vector: "<<int((*recvHalo)->size())<<" elements.\n";
	for (std::vector<mpihalo>::iterator iter=(*recvHalo)->begin();iter<(*recvHalo)->end();iter++)
	{
		ofs << "row: " << int(iter->row) << " col: " << int(iter->col) << " cnt: " << int(iter->cnt) << " rank: " << int(iter->rank) << "\n";
	}
	ofs.close();

	MPI_Barrier(MPI_COMM_WORLD);
#endif

#endif
	
// print the MPIREG-Grids --> must be done later when zones are analysed and FirstLastCol exists
//#ifdef _MPI
//
//	aktname = newstring(_stringsize_, 1);
//	if (mpi_rank == 0)
//	{
//		sprintf(aktname, "MPIRegions_%d", mpi_rank);
//		printgridfile(aktname, MPIREG, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//	}
//	sprintf(aktname, "MPIHelper_%d", mpi_rank);
//	printgridfile(aktname, MPIHLP, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//
//	sprintf(aktname, "MPIHalo_%d", mpi_rank);
//	printgridfile(aktname, MPIHALO, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//	free(aktname);
//	
//#endif
//  // also deleting these grids cannot be done here since the grids should still be printed later
//	if (mpi_rank == 0)
//	{
//		// the MPIHALO grid is not used any longer. SInce it was the last grid that was reserved, we can free it here without problems (there will be no gap in the gridinfo araay)
//		free_matrix(MPIHALO, 1, nrows, 1, ncols);
//		delete_grid(mpihaloindex);
//		gridanzahl--;
//	}
//	else
//	{
//		free_matrix(MPIHALO, sr, er, 1, ncols, nrows);
//		delete_grid(mpihaloindex);
//		gridanzahl--;
//	}

	delete[] v_wtime;
	delete[] v_mpiweight;
//	delete[] v_regioncnt; // will be used in other modules, so no deletion here but at model end
	delete[] locMpiStartRow;
	delete[] locMpiEndRow;

}

void defaultini::MakeHRUs(int &doSpinUp,int &doTemperaturePreSpinUp,int &repeatCnt,std::vector<HRUMapElement*> **vHRU)
{
	char *satz=NULL,*tempsatz=NULL;
	int numGrids=0;
	char *sSpacing=NULL;
	char *sName=NULL;
	double Classwidth=1;
	std::vector<double**> grids;
	std::vector<std::string> name; // grid identifier like 
	std::vector<double> classwidth; // equidistance between class boundaries for remapping values prior to creating HRUs
	std::vector<std::string> spacing; // lin or sqrt
	std::map<std::string,HRUMapElement*> mapHRU;
	char *cID;
	cID = new char[200];

	(*WSstream) << "DEBUG: Creating HRUs"<<"[end-entry]";

	mapHRU.clear();
	(*vHRU) = new vector<HRUMapElement*>();
	(*vHRU)->clear();
	// add a zero index element to the vector in order to have indices started from 1 to n rather from 0 to n-1 later.
	HRUMapElement *DummyMapElem;
	DummyMapElem=new HRUMapElement();
	DummyMapElem->ID = mapHRU.size();
	DummyMapElem->done = false;
	DummyMapElem->sID = "";
	DummyMapElem->MasterCol = 0;
	DummyMapElem->MasterRow = 0;
	(*vHRU)->push_back(DummyMapElem);

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
	mpi_rank=0;
	mpi_size=1;
#endif

	int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
	int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
	if (mpi_rank == 0)
	{
		sr = 1;
		er = nrows;
	}

	// 1. Looking for [SpinUp] in the control file
	satz    = (char *)calloc(_stringsize_+1,sizeof(char));
	tempsatz=satz;
	if ((steuerdatei = fopen(steuerdateiname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: Error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}

	if (!findeabschnitt(_SpinUp_,false))
	{   // if the new section [SpinUp] is not found in the control file, no error should be thrown and the program should go on without doing the SpinUp routines
		doSpinUp=0;
		doTemperaturePreSpinUp=0;
		repeatCnt=0;

		satz = tempsatz;
		if (satz) free(satz); satz = NULL;
		delete[] cID;

		return;
	}

	// 2. reading the parameters
	if (fscanf(steuerdatei,"%d",&doSpinUp)!=1) 
	{
		steuererror("SpinUp-switch not found","");	
	}

	if (doSpinUp == 0) 
	{
		doTemperaturePreSpinUp = 0;
		repeatCnt = 0;

		satz = tempsatz;
		if (satz) free(satz); satz = NULL;
		delete[] cID;

		return;
	}
	
	fgets(satz,_stringsize_,steuerdatei);
	
	if (fscanf(steuerdatei,"%d",&doTemperaturePreSpinUp)!=1) 
	{
		steuererror("doTemperaturePreSpinUp not found","");	fgets(satz,_stringsize_,steuerdatei);
	}

	if (fscanf(steuerdatei,"%d",&repeatCnt)!=1) 
	{
		steuererror("repeatCnt not found","");	
	}
	fgets(satz,_stringsize_,steuerdatei);

	if (fscanf(steuerdatei,"%d",&numGrids)!=1)
	{
		steuererror("numGrids not found","");
	}
	fgets(satz,_stringsize_,steuerdatei);

	// 4a: Creating or reading grid HRU
	if (numGrids==0)
	{
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
		if (!ist_vorhanden(_HydrologicResponseUnits_))
		{
			(*WSstream) << "FATAL: standard grid with identifier "<<_HydrologicResponseUnits_<<" not found (Hydrologic Respond Units, HRU; used for SpinUp mode)"<<"[end-entry]";
			throw(1030);
		}
		else
		{
			hruindex=getindex(_HydrologicResponseUnits_);
			HRU=gridinfo[hruindex].pointer;

			hrumasterindex=1;
			while (gridinfo[hrumasterindex].name!=NULL)
			{ 
				hrumasterindex++;
				if (hrumasterindex==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: when trying to create HRUMaster grid: no free grid found in memory (max position "<<hruindex<<")"<<"[end-entry]";
					throw(12);
				}
			}

			initialize_grid(hrumasterindex,_MasterHydrologicResponseUnits_,NULL,NULL,0,_MasterHydrologicResponseUnits_,0,sr,er);
			gridinfo[hrumasterindex].pointer = matrix(sr,er,1,ncols,nrows);
			HRUMaster = gridinfo[hrumasterindex].pointer;
			gridanzahl++;

			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					HRUMaster[i][j]=nodata;
					sprintf(cID,"%0.lf",HRU[i][j]);
					if (mapHRU.find(std::string(cID))!=mapHRU.end())
					{
						HRU[i][j]=mapHRU.find(cID)->second->ID; // already in the map, so simply set the HRU to the ID from the map
					}
					else
					{
						// new HRU
						HRUMapElement *mapElem;
						mapElem=new HRUMapElement();
						mapElem->ID = mapHRU.size()+1;
						mapElem->done = false;
						mapElem->sID = cID;
						mapElem->MasterCol = j;
						mapElem->MasterRow = i;
						(*vHRU)->push_back(mapElem);
						mapHRU.insert(std::pair<std::string,HRUMapElement*>(std::string(cID),mapElem));
						HRU[i][j]=HRUMaster[i][j]=mapHRU.size();

					}
				}
			}

		}
	}
	else
	{
		// read grid identifiers and other parameters, locate the grid in the gridinfo array, put all this values to the respective vectors
		sSpacing = new char[200];
		sName = new char[200];
		for (int i=1;i<=numGrids;i++)
		{
			if (fscanf(steuerdatei,"%s%lf%s",sName,&Classwidth,sSpacing)!=3)
			{
				steuererror(" error reading HRU defintions ","");
			}
			fgets(satz,_stringsize_,steuerdatei);

			name.push_back(std::string(sName));
			classwidth.push_back(Classwidth);
			spacing.push_back(std::string(sSpacing));
			grids.push_back(gridinfo[getindex_first_of_3(sName,(std::string("[")+sName+"]").c_str(),sName)].pointer);
		}
		delete[] sSpacing; sSpacing=NULL;
		delete[] sName; sName=NULL;

		hruindex=1;
		while (gridinfo[hruindex].name!=NULL)
		{ 
			hruindex++;
			if (hruindex==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: when trying to create HRU grid: no free grid found in memory (max position "<<hruindex<<")"<<"[end-entry]";
				throw(12);
			}
		}

		initialize_grid(hruindex,"All_HydrologicResponseUnits",NULL,NULL,0,_HydrologicResponseUnits_,0,sr,er);
		gridinfo[hruindex].pointer = matrix(sr, er, 1, ncols, nrows);
		HRU = gridinfo[hruindex].pointer;
		gridanzahl++;

		hrumasterindex=hruindex+1;
		initialize_grid(hrumasterindex,"Master_HydrologicResponseUnits",NULL,NULL,0,_MasterHydrologicResponseUnits_,0,sr,er);
		gridinfo[hrumasterindex].pointer = matrix(sr, er, 1, ncols, nrows);
		HRUMaster = gridinfo[hrumasterindex].pointer;
		gridanzahl++;

		for (int i=sr;i<=er;i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				HRU[i][j]=nodata;
				HRUMaster[i][j]=nodata;
				cID[0]='\0';

				for (int z=0;z<grids.size();z++)
				{
					if ((grids[z][i][j]==nodata)||(zn[i][j]==nodata)||(mpi_rank!=MPIREG[i][j])) // each MPI instance has its own set of HRUs
					{
						break; // leave Value at -9999
					}
					else
					{
						// 3. Creating the map and/or vector that will hold the HRU Infos
						// 4. Creating HRUs
						HRU[i][j]=0;
						double transformedValue=grids[z][i][j];
						if (spacing[z].find("sqrt")!=std::string::npos)
						{
							if (transformedValue<0) {
								transformedValue = sqrt(-transformedValue);
							} else if (transformedValue ==0) {
								transformedValue = 0;
							} else {
								transformedValue = sqrt(transformedValue);
							}

						}
						transformedValue = floor(transformedValue/classwidth[z]);
						if (transformedValue>9999)
						{
							(*WSstream) << "FATAL: Creating HRUs resulted in a re-classify result of "<<transformedValue
								<<" for grid "<<name[i]<<" at row="<<i<<" and col="<<j
								<<". Only 10'000 Classes are allowed for each grid. Choose a larger class width or use lg2 spacing to reduce the amount of HRUs."
								<<"[end-entry]";
							throw(1031);
						}
						sprintf(cID + strlen(cID), "%04.0lf'", transformedValue);
					}
				}
				// try to add the sID to the Map
				if (HRU[i][j]==0)
				{
					if (mapHRU.find(std::string(cID))!=mapHRU.end())
					{
						HRU[i][j]=mapHRU.find(cID)->second->ID; // already in the map, so simply set the HRU to the ID from the map
					}
					else
					{
						// new HRU
						HRUMapElement *mapElem;
						mapElem=new HRUMapElement();
						mapElem->ID = mapHRU.size()+1;
						mapElem->done = false;
						mapElem->sID = cID;
						mapElem->MasterCol = j;
						mapElem->MasterRow = i;
						(*vHRU)->push_back(mapElem);
						mapHRU.insert(std::pair<std::string,HRUMapElement*>(std::string(cID),mapElem));
						HRU[i][j]=HRUMaster[i][j]=mapHRU.size();

					}
				}
			}
		}

		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}
	}

	satz=tempsatz;
	if (satz) free(satz);satz=NULL;
	delete[] cID;

	mapHRU.clear();

	// writing the newly created grids
	char *gname;
	gname = new char[500];
	
#ifdef _MPI
	sprintf(gname,"%s_%03d",gridinfo[hruindex].name,mpi_rank);
	printgridfile(gname,gridinfo[hruindex].pointer,WriteASCIIGrids,WriteDoubleGrids, WriteNetcdfGrids);

	sprintf(gname,"%s_%03d",gridinfo[hrumasterindex].name,mpi_rank);
	printgridfile(gname,gridinfo[hrumasterindex].pointer,WriteASCIIGrids,WriteDoubleGrids, WriteNetcdfGrids);

	modelobject::MPIREG=gridinfo[getindex(_MPIRegions_)].pointer;
	modelobject::MPIHLP=gridinfo[getindex(_MPIHelper_)].pointer;

	MPIGatherGrid(hruindex);
	MPIGatherGrid(hrumasterindex);
#endif	

	if (mpi_rank==0)
	{
		sprintf(gname,"%s",gridinfo[hruindex].name);
		printgridfile(gname, gridinfo[hruindex].pointer, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);

		sprintf(gname,"%s",gridinfo[hrumasterindex].name);
		printgridfile(gname, gridinfo[hrumasterindex].pointer, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	}

	delete[] gname;
	// returning HRU-map or vector
}

LUSEtabletype::LUSEtabletype(int index):
CommonTableType(index),
method(0),iNumFixPoints(0),iNumMTStartLevels(0),StartVegetationPeriod(0),StopVegetationPeriod(1000),
JDVegetationReset(1),JDVegetationWrite(0),MTSTART(0),TSTART(0),FSTAR(0),FORC(0),Tstart(0),Tstart_col(0),days(0),
daystemp(0),rsc(0),rsi(0),rse(0),albedo(0),LAI(0),z0(0),Vbg(0),Wurz(0),AltDep(0),ForcingThreshold(0),form_par(0),
redu_parDry(0),redu_parWet(0),LimitReduWet(0),F_star(0),JDReset_TStart(0),DP1_t1_dorm(0),DP1_T_Bf(0),DP2_t0_dorm(0),
DP2_t1_dorm(0),DP2_T_Bc(0),DP2_T_Bf(0),DP2_Par_a(0),DP2_Par_b(0),DP2_Offset_1(0),DP2_Offset_2(0),DP4_t0_dorm(0),
DP4_t1_dorm(0),DP4_T0(0),DP4_T1(0),DP4_Par_a(0),DP4_Par_b(0),DP4_T_xylstop(0),DP4_t_xs1(0),DP4_SPSLF(0),
maxStartJDforDP1(0),useMethod3(0),m_fStressFactor(0),m_fIntercepStor(-1.0),m_SoilTillage(0),m_numSoilTillages(0),
SR_cultivation(0),SR_StemDiameter(0),SR_ResidCovrg(0),SR_CropDistX(0),SR_CropDistY(0),SR_K_ST(0),SR_C_WR(0),
SR_rough_coeff1(0),SR_rough_coeff2(0),sDynFlags("11111111"),dynflags(255),nfactor(0)
{
	mon[0]=31;
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

// defautlconstructor, setting code to 0
LUSEtabletype::LUSEtabletype():
CommonTableType(),
method(0),iNumFixPoints(0),iNumMTStartLevels(0),StartVegetationPeriod(0),StopVegetationPeriod(1000),
JDVegetationReset(1),JDVegetationWrite(0),MTSTART(0),TSTART(0),FSTAR(0),FORC(0),Tstart(0),Tstart_col(0),days(0),
daystemp(0),rsc(0),rsi(0),rse(0),albedo(0),LAI(0),z0(0),Vbg(0),Wurz(0),AltDep(0),ForcingThreshold(0),form_par(0),
redu_parDry(0),redu_parWet(0),LimitReduWet(0),F_star(0),JDReset_TStart(0),DP1_t1_dorm(0),DP1_T_Bf(0),DP2_t0_dorm(0),
DP2_t1_dorm(0),DP2_T_Bc(0),DP2_T_Bf(0),DP2_Par_a(0),DP2_Par_b(0),DP2_Offset_1(0),DP2_Offset_2(0),DP4_t0_dorm(0),
DP4_t1_dorm(0),DP4_T0(0),DP4_T1(0),DP4_Par_a(0),DP4_Par_b(0),DP4_T_xylstop(0),DP4_t_xs1(0),DP4_SPSLF(0),
maxStartJDforDP1(0),useMethod3(0),m_fStressFactor(0),m_fIntercepStor(-1.0),m_SoilTillage(0),m_numSoilTillages(0),
SR_cultivation(0),SR_StemDiameter(0),SR_ResidCovrg(0),SR_CropDistX(0),SR_CropDistY(0),SR_K_ST(0),SR_C_WR(0),
SR_rough_coeff1(0),SR_rough_coeff2(0),sDynFlags("11111111"),dynflags(255),nfactor(0)
{
	mon[0]=31;
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

LUSEtabletype::~LUSEtabletype()
{
	if (days)     delete[] days;
	if (daystemp) delete[] daystemp;
	if (rsc)      delete[] rsc;
	if (rse)      delete[] rse;
	if (rsi)      delete[] rsi;
	if (albedo)   delete[] albedo;
	if (LAI)      delete[] LAI;
	if (z0)       delete[] z0;
	if (Vbg)      delete[] Vbg;
	if (Wurz)     delete[] Wurz;
	if (AltDep)   delete[] AltDep;
	if (ForcingThreshold) delete[] ForcingThreshold;
	if (m_SoilTillage) delete[] m_SoilTillage;
	if (SR_rough_coeff1) delete[] SR_rough_coeff1;
	if (SR_rough_coeff2) delete[] SR_rough_coeff2;
	if (SR_StemDiameter) delete[] SR_StemDiameter;
	if (SR_ResidCovrg) delete[] SR_ResidCovrg;
	if (nfactor) delete[] nfactor;

}

void  LUSEtabletype::SetFStarGridVector(double*** &_FSTAR)
{
	FSTAR=_FSTAR;
}

void  LUSEtabletype::SetFORCGridVector(double*** &_FORC)
{
	FORC=_FORC;
}

void  LUSEtabletype::SetTStartGridVector(double*** &_TSTART)
{
	TSTART=_TSTART;
}

void  LUSEtabletype::SetMTStartStackVector(double**** &_MTSTART,int _iNumMTStartLevels)
{
	MTSTART=_MTSTART;
	iNumMTStartLevels=_iNumMTStartLevels;
}

double LUSEtabletype::getFormPar()
{
	return form_par;
}

double LUSEtabletype::getReduParDry()
{
	return redu_parDry;
}

double LUSEtabletype::getReduParWet()
{
	return redu_parWet;
}

double LUSEtabletype::getETPHamonFactor()
{
	return ETPHamonFactor;
}

double LUSEtabletype::getLimitReduWet()
{
	return LimitReduWet;
}

double LUSEtabletype::getIntercepStor()
{
	return m_fIntercepStor;
}

char* LUSEtabletype::getName()
{
	return Name;
}

int   LUSEtabletype::GetMethod()
{
	return method;
}

int   LUSEtabletype::GetTstart()
{
	return Tstart;
}

double* LUSEtabletype::GetForcingThresholds()
{
	return ForcingThreshold;
}

double LUSEtabletype::GetForcingThreshold(int i)
{
	if ((i>=0) && (i<=iNumFixPoints) && (method==5))
	return ForcingThreshold[i];
	else
	return -1;
}

double* LUSEtabletype::GetAltDeps()
{
	return AltDep;
}

double LUSEtabletype::GetAltDep(int i)
{
	if ((i>=0) && (i<=iNumFixPoints))
	return AltDep[i];
	else
	return 0;
}

int* LUSEtabletype::GetJulianDays()
{
	return days;
}

int LUSEtabletype::GetJulianDay(int i)
{
	if ((i>=0) && (i<=iNumFixPoints))
	return days[i];
	else
	return -1;
}

double LUSEtabletype::GetStressFactor()
{
	return m_fStressFactor;
}

void LUSEtabletype::GetSoilTillage(int &SoilTillageCount,int* &SoilTillageDays)
{
	SoilTillageCount = m_numSoilTillages;
	SoilTillageDays  = m_SoilTillage;
}

// will be used in int_tvar_params to set the starting day for method 3 when 
// summing up the temperature or for setting the Tstart to 0 at the beginning of a new year
void  LUSEtabletype::SetTstart(int day)
{
	Tstart=day; 
}

void LUSEtabletype::SetUseMethod3(bool mode)
{
	useMethod3=mode;
}

bool LUSEtabletype::GetUseMethod3()
{
	return useMethod3;
}


int   LUSEtabletype::getTstartcol()
{
	return Tstart_col;
}

double LUSEtabletype::GetFstar()
{
	return F_star;
}

int   LUSEtabletype::GetDP1_t1dorm()
{
	return DP1_t1_dorm;
}

double LUSEtabletype::GetDP1_TBf()
{
	return DP1_T_Bf;
}

int   LUSEtabletype::GetDP2_t1dorm()
{
	return DP2_t1_dorm;
}

int   LUSEtabletype::GetDP2_t0dorm()
{
	return DP2_t0_dorm;
}

double LUSEtabletype::GetDP2_TBf()
{
	return DP2_T_Bf;
}

double LUSEtabletype::GetDP2_TBc()     
{
	return DP2_T_Bc;
}

double LUSEtabletype::GetDP2_Par_a()      
{
	return DP2_Par_a;
}

double LUSEtabletype::GetDP2_Par_b()
{
	return DP2_Par_b;
}

double LUSEtabletype::GetDP2_Offset_1()
{
	return DP2_Offset_1;
}

double LUSEtabletype::GetDP2_Offset_2()
{
	return DP2_Offset_2;
}

int   LUSEtabletype::GetDP4_t0dorm()   // usually November, 1st; here we start the chilling units to add up
{
	return DP4_t0_dorm;
}

int   LUSEtabletype::GetDP4_t1dorm()   // usually February, 1st; here we start the forcing units to add up
{
	return DP4_t1_dorm;
}

double LUSEtabletype::GetDP4_T0()       // threshold temperature for chilling units (T < T0)
{
	return DP4_T0;
}

double LUSEtabletype::GetDP4_T1()       // threshold temperature for forcing untis (T > T1)
{
	return DP4_T1;
}

double LUSEtabletype::GetDP4_Par_a()    // empirical parameter a
{
	return DP4_Par_a;
}

double LUSEtabletype::GetDP4_Par_b()    // empirical parameter b
{
	return DP4_Par_b;
}

double LUSEtabletype::GetDP4_Txylstop() // temperature of the moving average which must be under-run for at least 5 days to stop xylem growth
{
	return DP4_T_xylstop;
}

int   LUSEtabletype::GetDP4_txs1()  // time of year as julian day, the evaluation of the moving average of temperature starts to estimate stop of xylem growth
{
	return DP4_t_xs1;
}

int   LUSEtabletype::GetDP4_SPSLF() // SamplePontStartLeaveFall
{
	return DP4_SPSLF;
}

int LUSEtabletype::GetJDReset_TStart()
{
	return JDReset_TStart;
}

int LUSEtabletype::GetMaxStartJDforDP1()
{
	return maxStartJDforDP1;
}

void LUSEtabletype::InitDynamicPhenology(int taga,int monata,int jahra)
{
	if ((GetMethod()==3)||(GetMethod()==4)||(GetMethod()==6))
	{
		SetUseMethod3(true); // the first year has no chilling units to calculate FStar from (after squential model 2), so we start with thermal time model
		// this will be set to false when DP2_t0_dorm or DP4_t0_dorm are reached for the first time
		if (GetJulDay(taga,monata,jahra)>GetMaxStartJDforDP1())
		{ // start day is later than the latest date allowed for dynamic phenology, so set 
			SetTstart(GetMaxStartJDforDP1()-this->days[getTstartcol()+1]);
		}
	}
}

double LUSEtabletype::interpolate(double x1,double y1,double x2,double y2,double x)
{
	if ((x2-x1)==0.0) return y1; else return (x-x1)/(x2-x1)*y2 + (x2-x)/(x2-x1)*y1;
}

int LUSEtabletype::GetJulDay(int iMonthDay,int iMonth,int iYear)
{
	int i=0;
	if ((iYear%4==0)&&(iYear%400!=0)) mon[2]=29; else mon[2]=28;
	for (int j=1;j<iMonth && j<12;j++)
	{
		i+=mon[j];
	}
	i+=iMonthDay;
	if ((i<0)||(i>366))
	{
		throw 150;
	}
	return i;
}

int LUSEtabletype::GetFixPoint(wasim::mystream *WSstream,int iJulDay,double height,int tstart)
{
	AltitudeCorrection(WSstream,height,tstart); // for Method 5 no days are calculated because they are already in MTSTART
	int k=1;
	for (int i=1;i<=iNumFixPoints;i++)
	{
		if ((iJulDay>daystemp[i])&&(daystemp[i]>=0))
		{
			k++;
		}
		else
		{
			break;
		}
	}
	if ((k<1)||(k>iNumFixPoints+1))
	{
		(*WSstream) << "FATAL: Julian-day-fixpoint at "<<k<<" is outside of the allowed range of "<<1<<" to "<<iNumFixPoints+1<<"[end-entry]";
		throw 151;
	}
	return k;
}

int  LUSEtabletype::getNumFixPoints()
{
	return this->iNumFixPoints;
}

unsigned int LUSEtabletype::getDynFlags()
{
	return dynflags;
}

std::string LUSEtabletype::getSDynFlags()
{
	return sDynFlags;
}

void	LUSEtabletype::setDynFlags(unsigned int flags)
{
	dynflags=flags;
}

unsigned int LUSEtabletype::isFlagSet(unsigned int flag)
{
	return dynflags&flag;
}

int LUSEtabletype::getStartVegetationPeriod()
{
	return StartVegetationPeriod;
}

int LUSEtabletype::getStopVegetationPeriod()
{
	return StopVegetationPeriod;
}


// optimized version (by Laura Carrington, staff of xseed)
void LUSEtabletype::AltitudeCorrection(wasim::mystream *WSstream, double dhm, int tstart)
{

	//      if (this->GetMethod()!=5) // war bis 21.04.2012 <5, hatte also auf DP4 die Wirkung, dass AltDep nie wirkte
	//      {
	//LAURA make FOR-loop easy to optimize by compiler...move if-then outside
	for (int i = 1; i <= iNumFixPoints; i++) //LAURA MOVED
	{ //LAURA MOVED
		daystemp[i] = days[i] + (int)(AltDep[i] * (dhm - 400.0));
	}
	if ((this->GetMethod() == 3) || (this->GetMethod() == 4) || (this->GetMethod() == 6))
	{
		int i = this->getTstartcol();
		// no altitude correction for TSTART and the next position
		// because the forcing units already depend on altitude via the temperature)
		daystemp[i] = tstart;

		i = this->getTstartcol() + 1;
		// next position after tstart is a relative value and must be added to tstart
		daystemp[i] = tstart + days[i];
	}

	for (int i = 1; i <= iNumFixPoints; i++) //LAURA MOVED
	{
		if ((daystemp[i]<-1) || (daystemp[i]>366))
		{
			(*WSstream) << "FATAL: \n---ERROR--\naltitude correction resulted in a Julian day of "
				<< daystemp[i] << " for fix point " << i << " (JD should be between 1 and 366)" << "[end-entry]";
			throw 152;
		}
	}
	daystemp[0] = days[iNumFixPoints] - 366;
}

// original version not optimized
//void LUSEtabletype::AltitudeCorrection(wasim::mystream *WSstream,double dhm,int tstart)
//{
////	if (this->GetMethod()!=5) // war bis 21.04.2012 <5, hatte also auf DP4 die Wirkung, dass AltDep nie wirkte
////	{
//		for (int i=1;i<=iNumFixPoints;i++)
//		{
//			if ((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))
//			{
//				if (i==this->getTstartcol())
//				{
//					// no altitude correction for TSTART and the next position 
//					// (because the forcing units already depend on altitude via the temperature)
//					daystemp[i] = tstart; 
//					if ((daystemp[i]<-1)||(daystemp[i]>366))
//					{
//						(*WSstream) << "FATAL: \n---ERROR--\ndynamic phenology (method 1 or 2 at TStart) resulted in a Julian day of "<<daystemp[i]<<" for fix point "<<i<<" (JD should be between 1 and 366)"<< "[end-entry]";
//						throw 152;
//					}
//				} else if (i==this->getTstartcol()+1) {
//					daystemp[i]=tstart+days[i]; // next position after tstart is a relative value and must be added to tstart
//					if ((daystemp[i]<-1)||(daystemp[i]>366))
//					{
//						(*WSstream) << "FATAL: \n---ERROR--\ndynamic phenology (method 1 or 2 after TStart) resulted in a Julian day of "<<daystemp[i]<<" for fix point "<<i<<" (JD should be between 1 and 366)"<< "[end-entry]";
//						throw 152;
//					}
//				} else {
//					daystemp[i] = days[i] + (int)(AltDep[i]*(dhm-400.0));
//					if ((daystemp[i]<-1)||(daystemp[i]>366))
//					{
//						(*WSstream) << "FATAL: \n---ERROR--\naltitude correction resulted in a Julian day of "<<daystemp[i]<<" for fix point "<<i<<" (JD should be between 1 and 366)"<< "[end-entry]";
//						throw 152;
//					}
//				}
//			}
//			else
//			{
//				daystemp[i] = days[i] + (int)(AltDep[i]*(dhm-400.0));
//				if ((daystemp[i]<-1)||(daystemp[i]>366))
//				{
//					(*WSstream) << "FATAL: \n---ERROR--\naltitude correction resulted in a Julian day of "<<daystemp[i]<<" for fix point "<<i<<" (JD should be between 1 and 366)"<< "[end-entry]";
//					throw 152;
//				}
//			}
//		}
//		daystemp[0] = days[iNumFixPoints] - 366;
////	}
//}

double LUSEtabletype::GetRsc(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],rsc[k-1],this->GetForcingThreshold(k),rsc[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],rsc[k-1],MTSTART[layer][i][j][k],rsc[k],iJulDay);
		}
	}
	else if (method>1) // all new methods using variable number of sample days as well as dynamic phenology 1 and 2
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==4))&&(k==this->getTstartcol()))
		{
			return rsc[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],rsc[k-1],daystemp[k],rsc[k],iJulDay);
		}
	}
	else // old method 1 using a fix number of sample days
	{
		if (iJulDay>365)
		{
			iJulDay=iJulDay%365;
		}
		int k=1;
		int z=mon[1];
		int remainder=iJulDay;

		while (iJulDay>z)
		{
			remainder-=mon[k];
			k++;
			z+=mon[k];
		}
		if (remainder==0) remainder=1;
		if ((remainder<1)||(remainder>mon[k])||(k<1)||(k>12))
		{
			throw 153;
		}
		return GetRsc(WSstream,remainder,k,2000,0,layer,i,j); // day=remainder, k=month, year=2000, altitude correction not required, so altitude=0
	}
}

double LUSEtabletype::GetRsc(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height,int layer,int i,int j)
{
	if (method>1) 
	{
		return GetRsc(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
	}
	else
	{
		int k;
		if (iMonthDay<15)
		{
			k=iMonth-1;
			iMonthDay+=(mon[iMonth-1]-1);
		}
		else
		{
			k=iMonth;
		}
		if (iMonth<0) iMonth=0;
		if ((iMonth<1)||(iMonth>12))
		{
			throw 154;
		}

		if ((k<0)||(k>12))
		{
			throw 155;
		}
		return interpolate(15,rsc[k],15+mon[iMonth],rsc[k+1],iMonthDay);
	}
}

double LUSEtabletype::GetRsi(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],rsi[k-1],this->GetForcingThreshold(k),rsi[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],rsi[k-1],MTSTART[layer][i][j][k],rsi[k],iJulDay);
		}
	}
	else if (method>1) // all new methods using variable number of sample days as well as dynamic phenology 1 and 2
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return rsi[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],rsi[k-1],daystemp[k],rsi[k],iJulDay);
		}
	}
	else
	{
		if (iJulDay>365)
		{
			iJulDay=iJulDay%365;
		}
		int k=1;
		int z=mon[1];
		int remainder=iJulDay;

		while (iJulDay>z)
		{
			remainder-=mon[k];
			k++;
			z+=mon[k];
		}
		if (remainder==0) remainder=1;
		if ((remainder<1)||(remainder>mon[k])||(k<1)||(k>12))
		{
			throw 153;
		}
		return GetRsi(WSstream,remainder,k,2000,0,layer,i,j); // day=remainder, k=month, year=2000, altitude correction not required, so altitude=0
	}
}

double LUSEtabletype::GetRsi(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height,int layer,int i,int j)
{
	if (method>1) 
	{
		return GetRsi(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
	}
	else
	{
		int k;
		if (iMonthDay<15)
		{
			k=iMonth-1;
			iMonthDay+=(mon[iMonth-1]-1);
		}
		else
		{
			k=iMonth;
		}
		if (iMonth<0) iMonth=0;
		if ((iMonth<1)||(iMonth>12))
		{
			throw 154;
		}

		if ((k<0)||(k>12))
		{
			throw 155;
		}
		return interpolate(15,rsi[k],15+mon[iMonth],rsi[k+1],iMonthDay);
	}
}

double LUSEtabletype::GetRse(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],rse[k-1],this->GetForcingThreshold(k),rse[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],rse[k-1],MTSTART[layer][i][j][k],rse[k],iJulDay);
		}
	}
	else if (method>1) // all new methods using variable number of sample days as well as dynamic phenology 1 and 2
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return rse[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],rse[k-1],daystemp[k],rse[k],iJulDay);
		}
	}
	else
	{
		if (iJulDay>365)
		{
			iJulDay=iJulDay%365;
		}
		int k=1;
		int z=mon[1];
		int remainder=iJulDay;

		while (iJulDay>z)
		{
			remainder-=mon[k];
			k++;
			z+=mon[k];
		}
		if (remainder==0) remainder=1;
		if ((remainder<1)||(remainder>mon[k])||(k<1)||(k>12))
		{
			throw 153;
		}
		return GetRse(WSstream,remainder,k,2000,0,layer,i,j); // day=remainder, k=month, year=2000, altitude correction not required, so altitude=0
	}
}

double LUSEtabletype::GetRse(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height,int layer,int i,int j)
{
	if (method>1) 
	{
		return GetRse(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
	}
	else
	{
		int k;
		if (iMonthDay<15)
		{
			k=iMonth-1;
			iMonthDay+=(mon[iMonth-1]-1);
		}
		else
		{
			k=iMonth;
		}
		if (iMonth<0) iMonth=0;
		if ((iMonth<1)||(iMonth>12))
		{
			throw 154;
		}

		if ((k<0)||(k>12))
		{
			throw 155;
		}
		return interpolate(15,rse[k],15+mon[iMonth],rse[k+1],iMonthDay);
	}
}

double LUSEtabletype::GetZ0(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],z0[k-1],this->GetForcingThreshold(k),z0[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],z0[k-1],MTSTART[layer][i][j][k],z0[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return z0[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],z0[k-1],daystemp[k],z0[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetZ0(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,double height,int layer,int i,int j)
{
	return GetZ0(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

// prerequisites: The value of "StartVegetationPeriodForBalance" and "StopVegetationPeriodForBalance" must be set 
// to non-zero values in the control file for every landuse, the vegetation period should be explicitely used for balances
// get the Julian day for the start of the vegetation period. 
// If the actual iJulDay is at the starting point of the phenological cycle (start of vegetation period), 
// then that day will be returned, otherwise -1
// works for all types of landuse table entries: fix sample nodes, dynamic phenology 1, 2 and 3
double LUSEtabletype::GetVegStart(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (this->getStartVegetationPeriod()==0) 
	{
		return 0; // if no "StartVegetationPeriodForBalance"-parameter was given, the intire year will be valid as vegetation period
	}
	else
	{
		if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
		{
			// get the last and the next fixpoint
			int k=1;
			while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
			// since k will be the next fixpoint after the last sample point, MTSTART was already calculated, we have to subtract 1 now
			k--;
			if (
					((k>=this->getStartVegetationPeriod())&&(k<this->getStopVegetationPeriod())) ||
					((k>=this->getStartVegetationPeriod())&&(this->getStartVegetationPeriod()>this->getStopVegetationPeriod())) // for southern hemisphere, start is later than stop within a year expressed as JD; however: dynamicPhenology_3 will probably not work for southern Hemisphere anyway
					)
			{
				if (iJulDay>=MTSTART[layer][i][j][k]) 
				{
					return iJulDay; // will return the actual JD even if VEGSTART was set to an earlier date! So this function should not be called once VEGSTART is > 0
				}
				else 
				{
					return -1;
				}
			}
			else
			{
				return -1;
			}
		}
		else 
		{
			int tstart;
			if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART)) 
			{ // DynamicPhenology_1 resp. DynamicPhenology_2
				tstart=TSTART[layer][i][j];
			}
			else
			{ // fixed sample points are only dependent altitude, but not on temperature, so tstart is always -1
				// meaning that it will not be used in GetFixPoint/AltitudeCorrection
				tstart=-1;
			}

			int k=GetFixPoint(WSstream,iJulDay,height,tstart);
			if (k>iNumFixPoints)
			{
				k=1;
				iJulDay-=366;
			}

			if (
					((k>=this->getStartVegetationPeriod())&&(k<this->getStopVegetationPeriod())) ||
					((k>=this->getStartVegetationPeriod())&&(this->getStartVegetationPeriod()>this->getStopVegetationPeriod())) // for southern hemisphere, start is later than stop within a year expressed as JD; will work for variableDaycount only, since dynamicPhenology may not work in sourhern hemisphere as implemented here
					)
			{
				if ((iJulDay>=daystemp[k])&&(daystemp[k]>0)) // will be true at the end of the interpolation interval only, thus we do not need to use k-1 here
				{
					return iJulDay; // will return the actual JD even if VEGSTART was set to an earlier date! So this function should not be called once VEGSTART is > 0
				}
				else 
				{
					return -1;
				}
			}
			else
			{
				return -1;
			}
		}
	}
}

double LUSEtabletype::GetVegStart(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,double height,int layer,int i,int j)
{
	return GetVegStart(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

// prerequisites: The value of "StartVegetationPeriodForBalance" and "StopVegetationPeriodForBalance" must be set 
// to non-zero values in the control file for every landuse, the vegetation period should be explicitely used for balances
// get the Julian day for the end of the vegetation period. 
// If the actual iJulDay is at the end point of the phenological cycle (end of vegetation period), 
// then that day will be returned, otherwise -1
// works for all types of landuse table entries: fix sample nodes, dynamic phenology 1, 2 and 3

int LUSEtabletype::getSR_cultivation()
{
	return SR_cultivation;
}

double LUSEtabletype::GetSR_CropDistX()
{
	return SR_CropDistX;
}

double LUSEtabletype::GetSR_CropDistY()
{
	return SR_CropDistY;
}

double LUSEtabletype::GetSR_K_ST()
{
	return SR_K_ST;
}

double LUSEtabletype::GetSR_C_WR()
{
	return SR_C_WR;
}

void LUSEtabletype::GetSR_rough_coeff1(double &c1,double &c2,double &c3)
{
	c1=SR_rough_coeff1[1];
	c2=SR_rough_coeff1[2];
	c3=SR_rough_coeff1[3];
}

void LUSEtabletype::GetSR_rough_coeff2(double &c1,double &c2,double &c3)
{
	c1=SR_rough_coeff2[1];
	c2=SR_rough_coeff2[2];
	c3=SR_rough_coeff2[3];
}

double LUSEtabletype::GetSR_StemDiameter(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],SR_StemDiameter[k-1],this->GetForcingThreshold(k),SR_StemDiameter[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],SR_StemDiameter[k-1],MTSTART[layer][i][j][k],SR_StemDiameter[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return SR_StemDiameter[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],SR_StemDiameter[k-1],daystemp[k],SR_StemDiameter[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetSR_StemDiameter(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height,int layer,int i,int j)
{
	return GetSR_StemDiameter(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetSR_ResidCovrg(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],SR_ResidCovrg[k-1],this->GetForcingThreshold(k),SR_ResidCovrg[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],SR_ResidCovrg[k-1],MTSTART[layer][i][j][k],SR_ResidCovrg[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return SR_ResidCovrg[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],SR_ResidCovrg[k-1],daystemp[k],SR_ResidCovrg[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetSR_ResidCovrg(wasim::mystream *WSstream,int iMonthDay,int iMonth,int iYear, double height,int layer,int i,int j)
{
	return GetSR_ResidCovrg(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetVegStop(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (this->getStopVegetationPeriod()==0) 
	{
		return 0; // if no "StartVegetationPeriodForBalance"-parameter was given, the entire year will be valid as vegetation period
	}
	else
	{
		if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
		{
			// get the last and the next fixpoint
			int k=1;
			while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
			// since k will be the next fixpoint after the last sample point, MTSTART was already calculated, we have to subtract 1 now
			k--;
			if ((k)==this->getStopVegetationPeriod())
			{
				if (iJulDay==MTSTART[layer][i][j][k]) 
				{
					return iJulDay;
				}
				else 
				{
					return -1;
				}
			}
			else
			{
				return -1;
			}
		}
		else 
		{
			int tstart;
			if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART)) 
			{ // DynamicPhenology_1 resp. DynamicPhenology_2
				tstart=TSTART[layer][i][j];
			}
			else
			{ // fixed sample points are only dependent altitude, but not on temperature, so tstart is always -1
				// meaning that it will not be used in GetFixPoint/AltitudeCorrection
				tstart=-1;
			}

			int k=GetFixPoint(WSstream,iJulDay,height,tstart);
			if (k>iNumFixPoints)
			{
				k=1;
				iJulDay-=366;
			}

			if (k==this->getStopVegetationPeriod())
			{
				if (iJulDay==daystemp[k]) // will be true at the end of the interpolation interval only, thus we do not need to use k-1 here
				{
					return iJulDay;
				}
				else 
				{
					return -1;
				}
			}
			else
			{
				return -1;
			}
		}
	}
}

double LUSEtabletype::GetVegStop(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,double height,int layer,int i,int j)
{
	return GetVegStop(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

void LUSEtabletype::setJDVegetationReset(int jd)
{
	JDVegetationReset=jd;
}

void LUSEtabletype::setJDVegetationWrite(int jd)
{
	JDVegetationWrite=jd;
}

int   LUSEtabletype::getJDVegetationReset()
{
	return JDVegetationReset;
}

int   LUSEtabletype::getJDVegetationWrite()
{
	return JDVegetationWrite;
}

double LUSEtabletype::GetLAI(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],LAI[k-1],this->GetForcingThreshold(k),LAI[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],LAI[k-1],MTSTART[layer][i][j][k],LAI[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return LAI[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],LAI[k-1],daystemp[k],LAI[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetLAI(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,double height,int layer,int i,int j)
{
	return GetLAI(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetVbg(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],Vbg[k-1],this->GetForcingThreshold(k),Vbg[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],Vbg[k-1],MTSTART[layer][i][j][k],Vbg[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return Vbg[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],Vbg[k-1],daystemp[k],Vbg[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetVbg(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,double height,int layer,int i,int j)
{
	return GetVbg(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetWurz(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	try
	{
		if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
		{
			// get the last and the next fixpoint
			int k=1;
			while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
			//		k++;
			// interpolate
			// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
			if (this->GetForcingThreshold(k)>0)
			{
				return interpolate(MTSTART[layer][i][j][0],Wurz[k-1],this->GetForcingThreshold(k),Wurz[k],FORC[layer][i][j]);
			}
			// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
			else
			{
				return interpolate(MTSTART[layer][i][j][k-1],Wurz[k-1],MTSTART[layer][i][j][k],Wurz[k],iJulDay);
			}
		}
		else 
		{
			int tstart;
			if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
			{
				tstart=TSTART[layer][i][j];
			}
			else
			{
				tstart=-1;
			}
			int k=GetFixPoint(WSstream,iJulDay,height,tstart);
			if (k>iNumFixPoints)
			{
				k=1;
				iJulDay-=366;
			}
			if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
			{
				return Wurz[k-1];
			}
			else
			{
				return interpolate(daystemp[k-1],Wurz[k-1],daystemp[k],Wurz[k],iJulDay);
			}
		}
	}
	catch (int ex)
	{
		(*WSstream) << "FATAL: Coordinates: i="<< i << ", j=" << j << "\nLayer=" << layer <<"\nAltitude="<<height<<"m\nJulday="<<iJulDay<<"\nLanduse Code="<<Code<<"\nTStart="<<this->GetTstart()<<"[end-entry]";
		throw;
	}
}

double LUSEtabletype::GetWurz(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear, double height,int layer,int i,int j)
{
	return GetWurz(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetAlbedo(wasim::mystream *WSstream,int iJulDay, double height,int layer,int i,int j)
{
	if (((method==5)||(method==6))&&(MTSTART!=NULL)) // dynamic phenology 3 using multiple FStar-Thresholds for multiple dynamic phenologic phases
	{
		// get the last and the next fixpoint
		int k=1;
		while ((iJulDay>=MTSTART[layer][i][j][k])&&(k<this->getNumFixPoints())&&(MTSTART[layer][i][j][k]>0)) k++;
		//		k++;
		// interpolate
		// a) using the MTSTART[layer][i][j][0] and the this->GetForcingThreshold(nextFixPoint) and the actual FORC-value
		if (this->GetForcingThreshold(k)>0)
		{
			return interpolate(MTSTART[layer][i][j][0],albedo[k-1],this->GetForcingThreshold(k),albedo[k],FORC[layer][i][j]);
		}
		// b) using the usual days from MTSTART[layer][i][j][k] and MTSTART[layer][i][j][k-1] and the actual julian day
		else
		{
			return interpolate(MTSTART[layer][i][j][k-1],albedo[k-1],MTSTART[layer][i][j][k],albedo[k],iJulDay);
		}
	}
	else 
	{
		int tstart;
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(TSTART))
		{
			tstart=TSTART[layer][i][j];
		}
		else
		{
			tstart=-1;
		}
		int k=GetFixPoint(WSstream,iJulDay,height,tstart);
		if (k>iNumFixPoints)
		{
			k=1;
			iJulDay-=366;
		}
		if (((this->GetMethod()==3)||(this->GetMethod()==4)||(this->GetMethod()==6))&&(k==this->getTstartcol()))
		{
			return albedo[k-1];
		}
		else
		{
			return interpolate(daystemp[k-1],albedo[k-1],daystemp[k],albedo[k],iJulDay);
		}
	}
}

double LUSEtabletype::GetAlbedo(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear, double height,int layer,int i,int j)
{
	return GetAlbedo(WSstream,GetJulDay(iMonthDay,iMonth,iYear),height,layer,i,j);
}

double LUSEtabletype::GetNFactorFreezing(wasim::mystream *WSstream,int i,int j)
{
	return this->nf_freezing;
}

double LUSEtabletype::GetNFactorThawing(wasim::mystream *WSstream,int i,int j)
{
	return this->nf_thawing;
}

double LUSEtabletype::GetArcVegInFactFull()
{
	return ArcVegInFactFull;
}

double LUSEtabletype::GetArcVegInFactSome()
{
	return ArcVegInFactSome;
}

double LUSEtabletype::GetArcVegInUDIceWdg()
{
	return ArcVegInUDIceWdg;
}

double LUSEtabletype::GetArcVegBNratio()
{
	return ArcVegBNratio;
}

double LUSEtabletype::GetArcVegSenesce()
{
	return ArcVegSenesce;
}

double LUSEtabletype::GetArcVegRetran()
{
	return ArcVegRetran;
}

double LUSEtabletype::GetArcVegFoliarProp()
{
	return ArcVegFoliarProp;
}

double LUSEtabletype::GetArcVegPropN()
{
	return ArcVegPropN;
}

double LUSEtabletype::GetArcVegAGProp()
{
	return ArcVegAGProp;
}

double LUSEtabletype::GetArcVegMaxBio()
{
	return ArcVegMaxBio;
}

double* LUSEtabletype::GetArcVegSeeedMaxProb()
{
	return ArcVegSeedMaxProb;
}

double* LUSEtabletype::GetArcVegUptakeEfficiency()
{
	return ArcVegUptakeEfficiency;
}

double LUSEtabletype::GetNfactorJD(wasim::mystream *WSstream,int iJulDay,int i,int j)
{
	// calculate the index of the sample point
	int k=1;
	for (int ii=1;ii<=iNumFixPoints;ii++)
	{
		if ((iJulDay>days[ii])&&(days[ii]>=0))
		{
			k++;
		}
		else
		{
			break;
		}
	}
	if ((k<1)||(k>iNumFixPoints+1))
	{
		(*WSstream) << "FATAL: Julian-day-fixpoint for n-factor at "<<k<<" is outside of the allowed range of "<<1<<" to "<<iNumFixPoints+1<<"[end-entry]";
		throw 151;
	}
	if (k>iNumFixPoints)
	{
		k=1;
		iJulDay-=366;
	}

	return interpolate(days[k-1],nfactor[k-1],days[k],nfactor[k],iJulDay);
}

double LUSEtabletype::GetNfactor(wasim::mystream *WSstream,int iMonthDay,int iMonth, int iYear,int i,int j)
{
	return GetNfactorJD(WSstream,GetJulDay(iMonthDay,iMonth,iYear),i,j);
}

int   LUSEtabletype::getCode()
{
	return Code;
}

void  LUSEtabletype::setCode(int code)
{
	Code = code;
}

// this is the intialization method for the variable landuse table using more than 4 fix points and for all dynamic phenology models
void LUSEtabletype::initialize(wasim::mystream *WSstream,std::string row,std::string sVersion)
{

	if (sVersion.compare("VariableDayCount")==0) {
		method=2; // the new method: using an open number of days for defining the phenology
	}	else if (sVersion.compare("DynamicPhenology_1")==0) {
		method=3; // for dynamic phenology 1 (thermal time model) we use a new method index
		//	}	else if (row.find("DynamicPhenology_2")!=std::string::npos)	{
	}	else if (sVersion.compare("DynamicPhenology_2")==0) {
		method=4; // for dynamic phenology 2 (sequential model 2) we use a new method index
	}	else if (sVersion.compare("DynamicPhenology_3")==0) {
		method=5; // for dynamic phenology 3 (thermal time model with multiple thresholds) we use a new method index
	}	else if (sVersion.compare("DynamicPhenology_4")==0) {
		method=6; // for dynamic phenology 4 (thermal time model similar to DP 2 with additional leave fall estimation) we use a new method index
	} else {
		(*WSstream) << "FATAL: landuse table method unknown "<< row.c_str() << "[end-entry]";
		throw(29);
	}

	useMethod3=false;
	Tstart=-1;

	std::stringstream ss; 
	ss << row.c_str();
	std::string buffer;
	int index1;
	int index2;

	//	std::vector<std::string> words;

	try
	{
		ss >> Code >> Name;

		row = row.substr(row.find("{"));

		buffer=row.substr(row.find("RootDist"));
		index1 = buffer.find("=")+1;
		index2 = buffer.find(";");
		buffer=buffer.substr(index1,index2-index1);
		form_par = atof( buffer.c_str());

		buffer=row.substr(row.find("TReduWet"));
		index1 = buffer.find("=")+1;
		index2 = buffer.find(";");
		buffer=buffer.substr(index1,index2-index1);
		redu_parWet = atof( buffer.c_str());

		if (row.find("useInternalParameters")==std::string::npos)
		{
			sDynFlags = "11111111";
			dynflags = 255;
		}
		else
		{
			buffer=row.substr(row.find("useInternalParameters"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			buffer = buffer.substr(buffer.find_first_of("01"));
			sDynFlags = buffer;
			dynflags=0;
			// no adding the single flags
			for (int cnt=0;cnt<8 && cnt<sDynFlags.length();cnt++)
			{
				if (sDynFlags.substr(cnt,1)=="1")
				{
					dynflags+=0x80 >> cnt;
				}
			}
		}

		if (row.find("LimitReduWet") == std::string::npos)
		{
			LimitReduWet = 0;
		}
		else
		{
			buffer = row.substr(row.find("LimitReduWet"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			LimitReduWet = atof(buffer.c_str());
		}

		if (row.find("ETPHamonFactor") == std::string::npos)
		{
			ETPHamonFactor = 1.0;
		}
		else
		{
			buffer = row.substr(row.find("ETPHamonFactor"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			ETPHamonFactor = atof(buffer.c_str());
		}


		buffer=row.substr(row.find("HReduDry"));
		index1 = buffer.find("=")+1;
		index2 = buffer.find(";");
		buffer=buffer.substr(index1,index2-index1);
		redu_parDry = atof( buffer.c_str());

		if (row.find("StartVegetationPeriodForBalance")==std::string::npos)
		{
			StartVegetationPeriod = 0;
		}
		else
		{
			buffer=row.substr(row.find("StartVegetationPeriodForBalance"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			StartVegetationPeriod = atoi( buffer.c_str());
		}

		if (row.find("StopVegetationPeriodForBalance")==std::string::npos)
		{
			StopVegetationPeriod = 0;
		}
		else
		{
			buffer=row.substr(row.find("StopVegetationPeriodForBalance"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			StopVegetationPeriod = atoi( buffer.c_str());
		}

		if (row.find("JDVegetationResetForBalance")==std::string::npos)
		{
			JDVegetationReset = -1;
		}
		else
		{
			buffer=row.substr(row.find("JDVegetationResetForBalance"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			JDVegetationReset = atoi( buffer.c_str());
		}

		if (row.find("JDVegetationWriteForBalance")==std::string::npos)
		{
			JDVegetationWrite = -1;
		}
		else
		{
			buffer=row.substr(row.find("JDVegetationWriteForBalance"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			JDVegetationWrite = atoi( buffer.c_str());
		}

		if (row.find("IntercepCap")==std::string::npos)
		{
			m_fIntercepStor = -1;
		}
		else
		{
			buffer=row.substr(row.find("IntercepCap"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			m_fIntercepStor = atof( buffer.c_str());
		}

		if (row.find("StressFactorDynPhen")==std::string::npos)
		{
			m_fStressFactor = 0;
		}
		else
		{
			buffer=row.substr(row.find("StressFactorDynPhen"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			m_fStressFactor = atof( buffer.c_str());
		}

		stripParametersJulDays(WSstream,row);
		
		// if dynamic phenology is used, a few new parameter have to be read
		if ((method==3)||(method==4)||(method==6))	{ // for methods 4 and 6 the parameters are used in fall back situations only (when model run starts after DP2_t0_dorm or DP4_t0_dorm, resp.)
			Tstart=-1;
			Tstart_col=-1;
			if (row.find("F*")!=std::string::npos)	{
				buffer=row.substr(row.find("F*"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();	ss >> F_star;
			}
			else {
				(*WSstream) << "FATAL: parameter \"F*\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP1_t1_dorm")!=std::string::npos)	{
				buffer=row.substr(row.find("DP1_t1_dorm"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();		ss >> DP1_t1_dorm;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"DP1_t1_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP1_T_Bf")!=std::string::npos)	{
				buffer=row.substr(row.find("DP1_T_Bf"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();	ss >> DP1_T_Bf;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"DP1_T_Bf\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("JDReset_TStart")!=std::string::npos)	{
				buffer=row.substr(row.find("JDReset_TStart"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();	ss >> JDReset_TStart;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"JDReset_TStart\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("maxStartJDforDP1")!=std::string::npos) {
				buffer=row.substr(row.find("maxStartJDforDP1"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> maxStartJDforDP1;
			}
			else {
				(*WSstream) << "FATAL: parameter \"maxStartJDforDP1\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			for (int i=1;i<iNumFixPoints;i++)
			{
				if (days[i]<0)
				{
					Tstart_col=i;
					break;
				}
			}

		}

		// some additional parameters required in sequential model for dynamic phenology
		// attantion: the parameters for method 3 are required as well, because in the first year 
		// the thermal time method (3) must be used (when starting on January 1st, the DP2_t0_dorm is 
		// over long ago, so FSTAR could not estiamted correctly thus leaving the only alternative to use 
		// the thermal time model for the first year.
		if (method==4) {
			if (row.find("DP2_t0_dorm")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_t0_dorm"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_t0_dorm;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_t0_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_t1_dorm")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_t1_dorm"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_t1_dorm;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_t1_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_T_Bf")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_T_Bf"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_T_Bf;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_T_Bf\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_T_Bc")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_T_Bc"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_T_Bc;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_T_Bc\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_Par_a")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_Par_a"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_Par_a;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_Par_a\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_Par_b")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_Par_b"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_Par_b;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP2_Par_b\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP2_Offset_1")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_Offset_1"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_Offset_1;
			}
			else {
				DP2_Offset_1=-3.4;
				(*WSstream) << "INFO: parameter \"DP2_Offset_1\" not found for land use "<<getCode() << "; Using default of " << DP2_Offset_1 << "[end-entry]";
			}

			if (row.find("DP2_Offset_2")!=std::string::npos) {
				buffer=row.substr(row.find("DP2_Offset_2"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP2_Offset_2;
			}
			else {
				DP2_Offset_2=12.0; // changed from 10.4 to 12.0 on March, 8, 2008 because of an errata from Prof. Chmielewski
				(*WSstream) << "INFO: parameter \"DP2_Offset_2\" not found for land use "<<getCode() << "; Using default of " << DP2_Offset_2 << "[end-entry]";
			}

		}
		// reading Parameters for method 5 (a subset of those needed for method 3 only)
		if (method==5)
		{
			if (row.find("DP1_t1_dorm")!=std::string::npos)	{
				buffer=row.substr(row.find("DP1_t1_dorm"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();		ss >> DP1_t1_dorm;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"DP1_t1_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP1_T_Bf")!=std::string::npos)	{
				buffer=row.substr(row.find("DP1_T_Bf"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();	ss >> DP1_T_Bf;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"DP1_T_Bf\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("JDReset_TStart")!=std::string::npos)	{
				buffer=row.substr(row.find("JDReset_TStart"));	index1 = buffer.find("=")+1;	index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str();	ss >> JDReset_TStart;
			}
			else	{
				(*WSstream) << "FATAL: parameter \"JDReset_TStart\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("maxStartJDforDP1")!=std::string::npos) {
				buffer=row.substr(row.find("maxStartJDforDP1"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> maxStartJDforDP1;
			}
			else {
				(*WSstream) << "FATAL: parameter \"maxStartJDforDP1\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}
		}
		if (method==6) // dynamic phenology 4 similar to DP2 but with other equations and parameters and with leave fall estimation for trees
		{
			if (row.find("DP4_t0_dorm")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_t0_dorm"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_t0_dorm;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_t0_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_t1_dorm")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_t1_dorm"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_t1_dorm;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_t1_dorm\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_T0")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_T0"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_T0;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_T0\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_T1")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_T1"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_T1;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_T1\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_Par_a")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_Par_a"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_Par_a;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_Par_a\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_Par_b")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_Par_b"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_Par_b;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_Par_b\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_T_xylstop")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_T_xylstop"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_T_xylstop;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_T_xylstop\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_t_xs1")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_t_xs1"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_t_xs1;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_t_xs1\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}

			if (row.find("DP4_SPSLF")!=std::string::npos) {
				buffer=row.substr(row.find("DP4_SPSLF"));	index1 = buffer.find("=")+1; index2 = buffer.find(";");
				buffer=buffer.substr(index1,index2-index1);
				std::stringstream ss;	ss << buffer.c_str(); ss >> DP4_SPSLF;
			}
			else {
				(*WSstream) << "FATAL: parameter \"DP4_SPSLF\" not found for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}
		}

		if (method==5) {
			ForcingThreshold=new double[iNumFixPoints+1];
		}
		rsc  = new double[iNumFixPoints+1];                    
		rsi  = new double[iNumFixPoints+1];                    
		rse  = new double[iNumFixPoints+1];                    
		albedo = new double[iNumFixPoints+1];
		LAI  = new double[iNumFixPoints+1];
		z0   = new double[iNumFixPoints+1];
		Vbg  = new double[iNumFixPoints+1];
		Wurz = new double[iNumFixPoints+1];
		AltDep = new double[iNumFixPoints+1];
		SR_rough_coeff1 = new double[3+1];
		SR_rough_coeff2 = new double[3+1];
		SR_StemDiameter = new double[iNumFixPoints+1];
		SR_ResidCovrg   = new double[iNumFixPoints+1];
		nfactor = new double[iNumFixPoints+1];

		if (method==5) {
			stripParameters(row,"ForcingThreshold",ForcingThreshold,iNumFixPoints);
		}
		stripParameters(row,"Albedo",albedo,iNumFixPoints);
		stripParameters(row,"rsc",rsc,iNumFixPoints);
		stripParameters(row,"LAI",LAI,iNumFixPoints);
		stripParameters(row,"Z0",z0,iNumFixPoints);
		stripParameters(row,"VCF",Vbg,iNumFixPoints);
		stripParameters(row,"RootDepth",Wurz,iNumFixPoints);
		stripParameters(row,"AltDep",AltDep,iNumFixPoints);

		int i;
		if (row.find("rs_interception")!=std::string::npos)
		{
			stripParameters(row,"rs_interception",rsi,iNumFixPoints);
		}
		else
		{
			for (i=0;i<=iNumFixPoints;i++)
			{
				rsi[i] = 20;
			}
		}

		if (row.find("rs_evaporation")!=std::string::npos)
		{
			stripParameters(row,"rs_evaporation",rse,iNumFixPoints);
		}
		else
		{
			for (i=0;i<=iNumFixPoints;i++)
			{
				rse[i] = 200;
			}
		}

		if (row.find("SoilTillage")!=std::string::npos)
		{
			stripParametersSoilTillageDays(WSstream,row);
		}
		else
		{
			m_numSoilTillages=0;
			m_SoilTillage = new int[m_numSoilTillages+1];
			m_SoilTillage[0]=-9999;
		}

		// parameters used in surface routing only, with default values
		if (row.find("SR_cultivation")!=std::string::npos)
		{
			buffer=row.substr(row.find("SR_cultivation"));
			index1 = buffer.find("=")+1;	
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			if (buffer.find("other")!=std::string::npos) {
				SR_cultivation=3;
			} else if (buffer.find("nonconserving")!=std::string::npos) {
				SR_cultivation=2;
			} else if (buffer.find("conserving")!=std::string::npos) {
				SR_cultivation=1;
			} else {
				(*WSstream) << "FATAL: no valid value for parameter \"SR_cultivation\" for land use "<<getCode()<<"[end-entry]";
				throw(17);
			}
		}
		else
		{
			this->SR_cultivation=3; // 1=conserving, 2=non_conserving, 3=other
		}

		if (row.find("SR_CropDistX")!=std::string::npos) {
			buffer=row.substr(row.find("SR_CropDistX"));	
			index1 = buffer.find("=")+1; 
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> SR_CropDistX;
		}
		else {
			SR_CropDistX=0.3;

		}

		if (row.find("SR_CropDistY")!=std::string::npos) {
			buffer=row.substr(row.find("SR_CropDistY"));	
			index1 = buffer.find("=")+1; 
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> SR_CropDistY;
		}
		else {
			SR_CropDistY=0.3;
		}

		if (row.find("SR_K_ST")!=std::string::npos) {
			buffer=row.substr(row.find("SR_K_ST"));	
			index1 = buffer.find("=")+1; 
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> SR_K_ST;
		}
		else {
			SR_K_ST=10.0;
		}

		if (row.find("SR_C_WR")!=std::string::npos) {
			buffer=row.substr(row.find("SR_C_WR"));	
			index1 = buffer.find("=")+1; 
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> SR_C_WR;
		}
		else {
			SR_C_WR=1.5;
		}

		if (row.find("SR_rough_coeff1")!=std::string::npos)
		{
			stripParameters(row,"SR_rough_coeff1",SR_rough_coeff1,3);
		}
		else
		{
			SR_rough_coeff1[0]= 0.0;
			SR_rough_coeff1[1]= 0.0198;
			SR_rough_coeff1[2]= 0.712;
			SR_rough_coeff1[3]= 0.142;
		}

		if (row.find("SR_rough_coeff2")!=std::string::npos)
		{
			stripParameters(row,"SR_rough_coeff2",SR_rough_coeff2,3);
		}
		else
		{
			SR_rough_coeff2[0]= 0.0;
			SR_rough_coeff2[1]= 0.0198;
			SR_rough_coeff2[2]= 0.712;
			SR_rough_coeff2[3]= 0.142;
		}

		if (row.find("SR_StemDiameter")!=std::string::npos)
		{
			stripParameters(row,"SR_StemDiameter",SR_StemDiameter,iNumFixPoints);
		}
		else
		{
			for (i=0;i<=iNumFixPoints;i++)
			{
				SR_StemDiameter[i] = 0.01;
			}
		}

		if (row.find("SR_ResidCovrg")!=std::string::npos)
		{
			stripParameters(row,"SR_ResidCovrg",SR_ResidCovrg,iNumFixPoints);
		}
		else
		{
			for (i=0;i<=iNumFixPoints;i++)
			{
				SR_ResidCovrg[i] = 0.01;
			}
		}

		if (row.find("n-factor")!=std::string::npos)
		{
			stripParameters(row,"n-factor",nfactor,iNumFixPoints);
		}
		else
		{
			for (i=0;i<=iNumFixPoints;i++)
			{
				nfactor[i] = -9999.9;
			}
		}

		if (row.find("nf_freezing")!=std::string::npos) {
			buffer=row.substr(row.find("nf_freezing"));	
			index1 = buffer.find("=")+1; 
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> nf_freezing;
		}
		else {
			nf_freezing=-9999.9;
		}

		if (row.find("nf_thawing") != std::string::npos) {
			buffer = row.substr(row.find("nf_thawing"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> nf_thawing;
		}
		else {
			nf_thawing = -9999.9;
		}

		// ArcVeg parameters
		if (row.find("ArcVegInFactFull") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegInFactFull"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegInFactFull;
		}
		else {
			ArcVegInFactFull = 1.0;
		}

		if (row.find("ArcVegInFactSome") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegInFactSome"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegInFactSome;
		}
		else {
			ArcVegInFactSome = 1.0;
		}

		if (row.find("ArcVegInUDIceWdg") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegInUDIceWdg"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegInUDIceWdg;
		}
		else {
			ArcVegInUDIceWdg = 1.0;
		}

		// biomass to N ration (each in g/m^2)
		if (row.find("ArcVegBNratio") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegBNratio"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegBNratio;
		}
		else {
			ArcVegBNratio = 1.0;
		}

		// fraction of senescing biomass
		if (row.find("ArcVegSenesce") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegSenesce"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegSenesce;
		}
		else {
			ArcVegSenesce = 1.0;
		}

		// plant type specific factor for retranslocated N
		if (row.find("ArcVegRetran") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegRetran"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegRetran;
		}
		else {
			ArcVegRetran = 1.0;
		}

		// foliar prop for each plant type (used for senesce)
		if (row.find("ArcVegFoliarProp") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegFoliarProp"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegFoliarProp;
		}
		else {
			ArcVegFoliarProp = 1.0;
		}

		// PropN is another factor used in ArcVeg
		if (row.find("ArcVegPropN") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegPropN"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegPropN;
		}
		else {
			ArcVegPropN = 1.0;
		}

		// AGProp is yet another factor used in ArcVeg
		if (row.find("ArcVegAGProp") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegAGProp"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegAGProp;
		}
		else {
			ArcVegAGProp = 1.0;
		}

		// maxium biomass for a m^2 (used in ArcVeg)
		if (row.find("ArcVegMaxBio") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegMaxBio"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss >> ArcVegMaxBio;
		}
		else {
			ArcVegMaxBio = 1.0;
		}

		ArcVegSeedMaxProb = new double[8];
		ArcVegSeedMaxProb[0] = 0;
		if (row.find("ArcVegSeedMaxProb") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegSeedMaxProb"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss
				>> ArcVegSeedMaxProb[1]
				>> ArcVegSeedMaxProb[2]
				>> ArcVegSeedMaxProb[3]
				>> ArcVegSeedMaxProb[4]
				>> ArcVegSeedMaxProb[5]
				>> ArcVegSeedMaxProb[6]
				>> ArcVegSeedMaxProb[7];
		}
		else {
			for (int k = 1; k <= 7; k++)
			{
				ArcVegSeedMaxProb[k] = 0;
			}
		}

		ArcVegUptakeEfficiency = new double[7];
		ArcVegUptakeEfficiency[0] = 0;
		if (row.find("ArcVegUptakeEfficiency") != std::string::npos) {
			buffer = row.substr(row.find("ArcVegUptakeEfficiency"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			std::stringstream ss;	ss << buffer.c_str(); ss
				>> ArcVegUptakeEfficiency[1]
				>> ArcVegUptakeEfficiency[2]
				>> ArcVegUptakeEfficiency[3]
				>> ArcVegUptakeEfficiency[4]
				>> ArcVegUptakeEfficiency[5];
		}
		else {
			for (int k = 1; k <= 5; k++)
			{
				ArcVegUptakeEfficiency[k] =	0;
			}
		}



	}
	catch (...)
	{
		(*WSstream) << "FATAL: landuse table code "<<Code<<": some values could not be read\n" << row.c_str() << "[end-entry]";
		throw(29);
	}
}

void LUSEtabletype::stripParameters(std::string &row,std::string sWhatParam,double *param,int maxIndex)
{
	int index1,index2;
	std::string buffer;
	buffer=row.substr(row.find(sWhatParam));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();
	for (int i=1;i<=maxIndex;i++)
	{
		ss >> param[i];
	}
	param[0]=param[maxIndex];
}

void LUSEtabletype::stripParametersJulDays(wasim::mystream *WSstream,std::string &row)
{
	int index1,index2;
	std::string buffer;
	std::vector<std::string> words;

	buffer=row.substr(row.find("JulDays"));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	words.clear();
	su.StringSplit(buffer,"\t ;,",words,false);
	iNumFixPoints=0;
	int i;
	for (i=0;i<words.size();i++)
	{
		if (words[i].size()>0)
		{
			iNumFixPoints++;
		}
	}
	days     = new int[iNumFixPoints+1];
	daystemp = new int[iNumFixPoints+1];
	std::stringstream ss;
	ss << buffer.c_str();
	for (i=1;i<=iNumFixPoints;i++)
	{
		ss >> days[i];
		daystemp[i]=days[i];
	}
	days[0]=days[iNumFixPoints]-366;
	daystemp[0]=daystemp[iNumFixPoints]-366;
	if (this->GetMethod()==5)
	{
		if (days[1]!=1)
		{
			(*WSstream) << "FATAL: Error in Landuse Table (Code: " << this->Code << "): for Method DynamicPhenology_3 the first sample day (JulDay) must be 1 and not " << days[1] <<"[end-entry]";
			throw(29);
		}
		if (days[iNumFixPoints]!=366)
		{
			(*WSstream) << "FATAL: Error in Landuse Table (Code: " << this->Code << "): for Method DynamicPhenology_3 the last sample day (JulDay) must be 366 and not " << days[iNumFixPoints] <<"[end-entry]";
			throw(29);
		}

	}
}

void LUSEtabletype::stripParametersSoilTillageDays(wasim::mystream *WSstream,std::string &row)
{
	int index1,index2;
	std::string buffer;
	std::vector<std::string> words;

	buffer=row.substr(row.find("SoilTillage"));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	words.clear();
	su.StringSplit(buffer,"\t ;,",words,false);
	m_numSoilTillages=0;
	int i;
	for (i=0;i<words.size();i++)
	{
		if (words[i].size()>0)
		{
			m_numSoilTillages++;
		}
	}
	m_SoilTillage = new int[m_numSoilTillages+1];
	std::stringstream ss;
	ss << buffer.c_str();
	for (i=1;i<=m_numSoilTillages;i++)
	{
		ss >> m_SoilTillage[i];
	}
	m_SoilTillage[0]=m_SoilTillage[m_numSoilTillages]-366;
}

void LUSEtabletype::initialize(wasim::mystream *WSstream,std::string row,int version)
{

	method=1; // the old static table

	std::stringstream ss; 
	ss << row.c_str();

	try
	{
		iNumFixPoints=4;
		days     = new   int[iNumFixPoints+1];
		daystemp = new   int[iNumFixPoints+1];
		rsc      = new double[14]; // it's special for the old landuse table to use 12 monthly surface resistances instead of 4 as for the other parameters
		rse      = new double[14]; // it's special for the old landuse table to use 12 monthly surface resistances instead of 4 as for the other parameters
		rsi      = new double[14]; // it's special for the old landuse table to use 12 monthly surface resistances instead of 4 as for the other parameters
		albedo   = new double[iNumFixPoints+1];
		LAI      = new double[iNumFixPoints+1];
		z0       = new double[iNumFixPoints+1];
		Vbg      = new double[iNumFixPoints+1];
		Wurz     = new double[iNumFixPoints+1];
		AltDep   = new double[iNumFixPoints+1];


		ss >> Code >> Name >> albedo[0];
		if (!ss.good())
		{
			(*WSstream) << "FATAL: landuse table code "<<Code<<": some values could not be read\n" << row.c_str() << "[end-entry]";
			throw(29);
		}
		albedo[1]=albedo[2]=albedo[3]=albedo[4]=albedo[0];
		int i;
		for (i=1;i<=12;i++)
		{
			ss >> rsc[i];
			rse[i] = 200;    // the old LU-table did not contain this parameter. So let's use a default value.
			rse[i] = rsc[i]; // the old LU-table did not contain this parameter. So let's use a default value.
		}
		rsc[0]=rsc[12];
		rsc[13]=rsc[1];
		rse[0]=rse[12];
		rse[13]=rse[1];
		rsi[0]=rsi[12];
		rsi[13]=rsi[1];

		int day;
		for (i=1;i<=iNumFixPoints;i++)
		{
			ss >> day;
			days[i]=day;
			daystemp[i]=day;
		}
		days[0]=days[iNumFixPoints]-366;
		daystemp[0]=daystemp[iNumFixPoints]-366;

		for (i=1;i<=iNumFixPoints;i++)
		{
			ss >> LAI[i];
		}
		LAI[0]=LAI[iNumFixPoints];

		for (i=1;i<=iNumFixPoints;i++)
		{
			ss >> z0[i];
		}
		z0[0]=z0[iNumFixPoints];

		for (i=1;i<=iNumFixPoints;i++)
		{
			ss >> Vbg[i];
		}
		Vbg[0]=Vbg[iNumFixPoints];

		for (i=1;i<=iNumFixPoints;i++)
		{
			ss >> Wurz[i];
		}
		Wurz[0]=Wurz[iNumFixPoints];

		if (version==2)
		ss >> form_par >> redu_parDry;

		/* default-Value for reducing transpiration under anaerob conditions in the old model version */
		redu_parWet = 0.95; /* ab welcher Bodenfeuchte im Nassbereich beginnt Verdunstungsreduktion, einlesen als relativer Wassergehalt (0..1), Umrechnen auf absoluten Wassergehalt Theta */
		LimitReduWet = 0; /* minimal relative reduction of real transpiration when water content reaches saturation. The reduction factor value will go down linearly starting at 1.0 when relative Theta equals TReduWet (e.g. 0.95) to LimitReduWet when the soil is saturated (Theta rel = 1.0) */

		AltDep[1] = 0.025;
		AltDep[2] = 0.025;
		AltDep[3] = -0.025;
		AltDep[4] = -0.025;
		AltDep[0] = AltDep[iNumFixPoints];
	}
	catch (...)
	{
		(*WSstream) << "FATAL: landuse table code "<<Code<<": some values could not be read\n" << row.c_str() << "[end-entry]";
		throw(29);
	}

}

// now the new methods for the soil table class

SOILtabletype::SOILtabletype(int index):
CommonTableType(index),
parmethode(0),nFK(0),
mSB(0),drain(0),k(0),psi(0),anz_k(0),anz_h(0),
levels(0),theta_k(NULL),k_tab(NULL),theta_h(NULL),
h_tab(NULL),k_sat(NULL),dzui(NULL),thetas(0),
thetafc(0),thetar(0),n(0),m(0),alpha(0),tau(0),k_rec(NULL),k_abnahme(0),
maxratio(0),LevelPerHorizon(0),
dPMacroThreshold(1.0),dPMacroCapacity(0.0),dMacroDepth(0.0),dReductionPerMeter(1.0),
dFractionSand(0.3333),dFractionSilt(0.3333),dFractionClay(0.3334),
dFractionStones1(0),dFractionStones2(0),dFractionStones3(0),dFractionStones4(0),
DryHeatCapacity(0),DryDensity(0),DryThermalConduct(0),KMinFrozenSoil(-1),lnKMinFrozenSoil(-20),
dEOverdT(NULL),Ts_RelTheta(NULL),SE_T(NULL),anz_T(1000),PhaseChngMinTemp(-10.0),dEvapMaxDepth(0.1),
drySoilThermalCondEff(-1),WaterThermalCondEff(-1)
//,ksat_SE(NULL)
{
}

SOILtabletype::SOILtabletype():
CommonTableType(),parmethode(0),nFK(0),
mSB(0),drain(0),k(0),psi(0),anz_k(0),anz_h(0),
levels(0),theta_k(NULL),k_tab(NULL),theta_h(NULL),
h_tab(NULL),k_sat(NULL),dzui(NULL),thetas(0),
thetafc(0),thetar(0),n(0),m(0),alpha(0),tau(0),k_rec(NULL),k_abnahme(0),
maxratio(0),LevelPerHorizon(0),
dPMacroThreshold(1.0),dPMacroCapacity(0.0),dMacroDepth(0.0),dReductionPerMeter(1.0),
dFractionSand(0.3333),dFractionSilt(0.3333),dFractionClay(0.3334),
dFractionStones1(0),dFractionStones2(0),dFractionStones3(0),dFractionStones4(0),
DryHeatCapacity(0),DryDensity(0),DryThermalConduct(0),KMinFrozenSoil(-1),lnKMinFrozenSoil(-20),
dEOverdT(NULL),Ts_RelTheta(NULL),SE_T(NULL),anz_T(1000),PhaseChngMinTemp(-10.0),dEvapMaxDepth(0.1),
drySoilThermalCondEff(-1),WaterThermalCondEff(-1)
//,ksat_SE(NULL)
{
}

SOILtabletype::~SOILtabletype()
{
	if (dzui    ) delete[] dzui    ; 
	if (k_sat   ) delete[] k_sat   ; 
	int i;
	if (theta_k ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] theta_k[i];
		}
		delete[] theta_k ; 
	}
	if (k_tab   ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] k_tab[i];
		}
		delete[] k_tab   ; 
	}
	if (theta_h ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] theta_h[i];
		}
		delete[] theta_h ; 
	}
	if (h_tab   ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] h_tab[i];
		}
		delete[] h_tab   ; 
	}
//	if (ksat_SE   ) 
//	{
//		for (int i=0;i<=levels;i++)
//		{
//			delete[] ksat_SE[i];
//		}
//		delete[] ksat_SE   ; 
//	}
	if (SE_T   ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] SE_T[i];
		}
		delete[] SE_T   ; 
	}
	if (dEOverdT   ) 
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] dEOverdT[i];
		}
		delete[] dEOverdT; 
	}
	if (Ts_RelTheta)
	{
		for (int i=0;i<=levels;i++)
		{
			delete[] Ts_RelTheta[i];
		}
		delete[] Ts_RelTheta; 
	}

	
	if (thetas  ) delete[] thetas  ;
	if (thetar  ) delete[] thetar  ;
	if (m       ) delete[] m       ;
	if (n       ) delete[] n       ;
	if (alpha   ) delete[] alpha   ;
	if (tau     ) delete[] tau     ;
	if (k_rec   ) delete[] k_rec   ;
	if (LevelPerHorizon)	delete[] LevelPerHorizon;
	if (DryDensity)			delete[] DryDensity;
	if (DryHeatCapacity)	delete[] DryHeatCapacity;
	if (DryThermalConduct)	delete[] DryThermalConduct;

}

void SOILtabletype::LookUpTableGenerationTMod(int layer,double factor)
{
	// temperatures are linearly distributed between -10°C and 0°C with 0°C beeing the 0-index already which is valid for temperatures at and above 0°+C
	SE_T[layer][0]=1.0; // water saturation is 1 at and above 0°C
	dEOverdT[layer][0]=0; // no energy change above 0°C
	Ts_RelTheta[layer][0]=abs(PhaseChngMinTemp); // temperature shift for rel theta at or below below 0 limited to the actual regarded temperature
	double dT=PhaseChngMinTemp/anz_T;

	double MaxdEOverdT=0; // temporary value for the maximum dE/dT to scale temperature shift
	double T_maxdEdT=0;

	for (int i=1;i<=anz_T;i++)
	{
		// SE_T-table
		double T = i*dT;
		SE_T[layer][i] = pow(double(1.0/(1.0+pow(double(this->getAlpha(layer)*fabs(factor*T)),double(this->getN(layer))))),double(this->getM(layer)));
		// dEOverdT table --> attention: the value calculated here must be multiplied by dTModRhoDrySoil*dTModLatHeat*(Theta_act-Theta_res) in order to accounting for the volumetric water or ice content of the soil
		// so this is only the scaling factor (integrated by T it should result in 1.0)
		double AlphaFactorT = fabs(T*factor*this->getAlpha(layer)); // alpha*f*T
		double temp1=pow(AlphaFactorT,double(this->getN(layer)-1));
		double temp2=pow(double(1.0+pow(AlphaFactorT,double(this->getN(layer)))),double(-1.0*this->getM(layer)-1));
		double temp3=this->getM(layer)*this->getN(layer)*this->getAlpha(layer)*factor;
		dEOverdT[layer][i]=double(temp1*temp2*temp3);
		if (dEOverdT[layer][i]>MaxdEOverdT)
		{
			MaxdEOverdT = dEOverdT[layer][i];
			T_maxdEdT   = T;
		}

	}

	// now the temperature shift for unsaturated soils
	double dTheta=1.0/anz_T;
	for (int i=0;i<=anz_T;i++)
	{
		// first getting the maximum dE/dT value for the actual reduced theta
		double relTheta=i*dTheta;
		double reducedMaxdE = relTheta * MaxdEOverdT;
		double T_redMaxdE = 0; // temperature for that reduced maximum in the unsaturated dE/dT curve
		// now searching the corresponding temperature for that energy change value
		int j=anz_T; // j is the temperature index
		if (dEOverdT[layer][j]>reducedMaxdE)
		{
			T_redMaxdE = j*dT;
		}
		else
		{
			while ((j>0)&&(dEOverdT[layer][j]<reducedMaxdE))
			{
				j--;
			}
			double weight=(reducedMaxdE-dEOverdT[layer][j+1])/(dEOverdT[layer][j]-dEOverdT[layer][j+1]);
			T_redMaxdE = weight*j*dT + (1.0-weight)*(j+1)*dT;
		}
		double fc=1.15;
		if (T_redMaxdE>=T_maxdEdT)
		{
			T_redMaxdE=T_maxdEdT;
			Ts_RelTheta[layer][i]=0;
		}
		// left of the double offset from the maximum dE/dT (e.g. offset is 0.15 when factor fc is 1.15, so the double offset is 2*0.15=0.3*T(of maxdE/dT) the factor is applied as a constant
		// but between that offset and the maximum dE/dT itself, the factor is linearly degraded to 0
		// now the branch when the temperature of reduced dE/dT-maximum is still lower than that offset 
		else if (T_redMaxdE<T_maxdEdT*(1+2*(fc-1))) 
		{
			Ts_RelTheta[layer][i] = T_redMaxdE - fc*T_maxdEdT;
		}
		else
		// and now the linear degrading of the factor to 1.0
		{
			Ts_RelTheta[layer][i] = T_redMaxdE - T_maxdEdT*(1.0+(fc-1.0)*(T_redMaxdE - T_maxdEdT)/(2*(fc-1)*T_maxdEdT));
		}
	}

	// ksat_SE will not be calculated anymore by this lookup-table but analytically in k_calc since we changed the algorithms for unsaturated soils
	// in k_calc, theta is replaced by SE if SE becomes smaller than theta
//	ksat_SE[layer][0]=this->getKMinFrozenSoil();
//	double dSE = 1.0/anz_T; // change rate per index increment in the look up table
//	for (int i=1;i<=anz_T;i++)
//	{
//		// ksat_SE-table
//		double SE_act=i*dSE;
//		double ln_ksat=k_sat[layer]<=0?-20.0:log(k_sat[layer]);
//		// contains now the logarithmic average fo the saturated hydraulic conductivity when all water is liquid and the saturated hydraulic conductivity when all water is ice
//		// weighed by the relative water saturation SE
//		// in fact, ksat is here a function of T, since SE is a function of T only (and soil parameters, of course, which are constant for each soil horizon)
//		ksat_SE[layer][i] = exp(SE_act*ln_ksat+(1.0-SE_act)*this->getLnKMinFrozenSoil());
//	}

}

void SOILtabletype::LookUpTableGeneration(int layer)
{

	theta_k[layer] = new double[anz_k+1];
	k_tab[layer]   = new double[anz_k+1];
	theta_h[layer] = new double[anz_h+1];
	h_tab[layer]   = new double[anz_h+1];

	theta_k[layer][0] = theta_k[layer][1] = theta_h[layer][0] = thetas[layer];
	theta_k[layer][anz_k] = theta_h[layer][anz_h] = thetar[layer];

	/* Berechnung der theta_h, theta_k und k_tab (h_tab wird vorgegeben) */
	h_tab[layer][1]  = 0;     /* Saugspannung in m */
	h_tab[layer][2]  = 0.0005;
	h_tab[layer][3]  = 0.001;
	h_tab[layer][4]  = 0.002;
	h_tab[layer][5]  = 0.005;
	h_tab[layer][6]  = 0.007;
	h_tab[layer][7]  = 0.01;
	h_tab[layer][8]  = 0.02;
	h_tab[layer][9]  = 0.03;
	h_tab[layer][10] = 0.05;
	h_tab[layer][11] = 0.07;
	h_tab[layer][12] = 0.1;
	h_tab[layer][13] = 0.15;
	h_tab[layer][14] = 0.2;
	h_tab[layer][15] = 0.3;
	h_tab[layer][16] = 0.4;
	h_tab[layer][17] = 0.6;
	h_tab[layer][18] = 0.8;
	h_tab[layer][19] = 1.2;
	h_tab[layer][20] = 1.6;
	h_tab[layer][21] = 2.6;
	h_tab[layer][22] = 3.45;
	h_tab[layer][23] = 6.9;
	h_tab[layer][24] = 20;
	h_tab[layer][25] = 35;
	h_tab[layer][26] = 50;
	h_tab[layer][27] = 75;
	h_tab[layer][28] = 100;
	h_tab[layer][29] = 125;
	h_tab[layer][30] = 150;

	int ak;
	double temp1,temp2,temp3;
	for (ak=1;ak<anz_h;ak++)
	{
		/* Berechnung des relativen Theta-Wertes */
		theta_h[layer][ak]=pow(1.0/(1.0+pow(alpha[layer]*h_tab[layer][ak],n[layer])),m[layer]);
		/* umrechnen in realen theta-Wert */
		theta_h[layer][ak]*=(thetas[layer]-thetar[layer]);
		theta_h[layer][ak]+=thetar[layer];
		if (ak<anz_k)
		{
			theta_k[layer][ak]=theta_h[layer][ak];
			/* Berechnen des k-Wertes */

			// old
			//		temp1 = 1.0-pow(alpha[layer]*h_tab[layer][ak],n[layer]-1.0)*pow(1.0+pow(alpha[layer]*h_tab[layer][ak],n[layer]),-m[layer]);
			//		temp1*=temp1;
			//		temp2 = pow(1.0+pow(alpha[layer]*h_tab[layer][ak],n[layer]),m[layer]*tau[layer]);
			//		k_tab[layer][ak] = temp1/temp2;

			// new
			temp1 = pow( 1.0/(1.0+pow(alpha[layer]*h_tab[layer][ak],n[layer])),m[layer]); // relative theta (fraction of Theta_sat-theta_res)
			temp2 = pow(temp1,tau[layer]);
			temp3 = pow(1.0-pow(1.0-pow(temp1,1.0/m[layer]),m[layer]),2.0);
			k_tab[layer][ak] = temp2*temp3;
			
		}
	}
	k_tab[layer][anz_k]=0.0; /* bei Erreichen der Restsaettigung ist Leitfaehigkeit 0 */
}

// initializing for the old soiltable layout without multiple horizons
void SOILtabletype::initialize(wasim::mystream *WSstream,std::string row,int version)
{
	int ak,i,layer;

	std::stringstream ss; 
	ss << row.c_str();
	double thetas_temp,thetar_temp,alpha_temp,n_temp,m_temp;

	try
	{
		/* Einlesen des Codes, des Namens und der Parameter fuer Topmodel version*/
		ss >> Code >> Name >> nFK >> mSB >> k >> psi;
		/* Berechnen des durch gravitation entwaesserbaren Porenraumes n_e*/
		drain=(mSB-nFK)*0.01;

		if (version==2)
		{
			//-------------------------------------------------------------------

			/* ########### Nun erfolgt das Einlesen der fuer wasimuzr noetigen Parameter der einzelnen Schichten des Bodens ########## */
			ss >> parmethode;
			if (parmethode==1)
			{  /* erzeugen der k-theta- und h-theta-Tabellen durch van-Genuchten Parameter */
				ss >> thetas_temp >> thetar_temp >> alpha_temp >> n_temp >> levels;
				m_temp    = 1.0 - 1.0/n_temp; /* parameter m fuer van-genuechten-Pedotransfer-function */
				
				dzui   = new double[levels+1];
				k_sat  = new double[levels+1];
				thetas = new double[levels+1];
				thetar = new double[levels+1];
				m      = new double[levels+1];
				n      = new double[levels+1];
				alpha  = new double[levels+1];
				tau    = new double[levels+1];

				for (i = 0;i<=levels;i++)
				{
					thetas[i] = thetas_temp; 
					thetar[i] = thetar_temp; 
					m     [i] = m_temp; 
					n     [i] = n_temp; 
					alpha [i] = alpha_temp; 
					tau   [i] = 0.5;
				}
				

				try { ss >> dzui[0];  } catch (...) {(*WSstream) << "FATAL: error reading layer thickness for soiltype "<<i<<"[end-entry]"; throw(30);}
				try { ss >> maxratio; } catch (...) {(*WSstream) << "FATAL: error reading maxratio for soil type "      <<i<<"[end-entry]"; throw(30);}
				try { ss >> k_abnahme;} catch (...) {(*WSstream) << "FATAL: error reading k reduction for soil type "<<i<<"[end-entry]"; throw(30);}

				anz_k=30;
				anz_h=30;

				typedef double* dp ;
				theta_k = new dp[levels+1];
				k_tab   = new dp[levels+1];
				theta_h = new dp[levels+1];
				h_tab   = new dp[levels+1];

// set in constructor already				anz_T=10000; // 0.01°C between 0 and -10°C
				SE_T      = new dp[levels+1]; // lookup table for relative saturation as a function of temperature (per soil layer!)
				dEOverdT  = new dp[levels+1]; // lookup table for internal energy change (per soil layer!)
				Ts_RelTheta = new dp[levels+1]; // lookup table for temperature shifts
//				ksat_SE   = new dp[levels+1]; // lookup table for internal energy change (per soil layer!)

				for (layer=0;layer<=levels;layer++)
				{
					// calling tab generation here
					LookUpTableGeneration(layer);
					// allocation is done here in order to make sure the values exist (because they will be deleted from this class, as well, and then no query is made if all vectors were allocated correctly
					// even for index 0!
					SE_T[layer] = new double[anz_T+1];
					dEOverdT[layer] = new double[anz_T+1];
					Ts_RelTheta[layer] = new double[anz_T+1];
//					ksat_SE[layer] = new double[anz_T+1];
					// will be called from unsatzon-constructor since we know the factor only there
					// LookUpTableGenerationTMod(layer);
				}

				/* Schichtdicken berechnen (first Layer: constant value from modconst.h) */
				if (dzui[0]<2*_Dicke_Schicht_1_old_)
				{
					(*WSstream) << "FATAL: \nerror soil type "<<Code<<":\n"
					"thickness of single layers (here "<<dzui[0]<<" m)\n"
					"has to be at minimum twice as thick as the\n"
					"thickness of the top layer (here "<<_Dicke_Schicht_1_old_<<" m)"<<"[end-entry]";
					throw(30);
				}
				dzui[1]=_Dicke_Schicht_1_old_; /* erste Schicht ist 5 cm dick */
				/* zweite Schicht ist 5 cm duenner als die Standardschichtdicke, um insgesamt
		in der Summe der ersten und zweiten Schicht auf dzui zu kommen */
				dzui[2]=dzui[0]-dzui[1];
				for (ak=3;ak<=levels;ak++)
				dzui[ak]=dzui[0];

				/* calculation of k_sat - values depending on the relative reduction per m (k_abnahme) */
				k_sat[1] = k*pow(k_abnahme,_Dicke_Schicht_1_old_/2.0); /* erste Schicht ist immer _Dicke_Schicht_1_ m dick (in modconst.h) */
				k_sat[2] = k*pow(k_abnahme,_Dicke_Schicht_1_old_+(dzui[2]-_Dicke_Schicht_1_old_)/2.0);
				for (ak=3;ak<=levels;ak++)
				{
					k_sat[ak] = ((double)ak-1.5)*dzui[0] <= 2.0 ? /* below 2.0 meter depth k_sat is assumed to be constant */
					k*pow(k_abnahme,((double)ak-1.5)*dzui[0]) :
					k_sat[ak-1];
				}
			}
			else // if (parmethode==2)
			{
				(*WSstream) << "FATAL: error reading soil parameters for soil type "<<Code<<": only old method 1 (van-Genuchten Parameter) is supported in WaSiM newer than 6.4.\n"
				<< "Use the new soil type table format instead."<<"[end-entry]";
				throw(131);
			}
		}
	}
	catch (...)
	{
		(*WSstream) << "FATAL: soil table code "<<Code<<": some values could not be read\n" << row.c_str() << "[end-entry]";
		throw(30);
	}
}

// if the temperature model is active, the relative saturation must be taken into account for calculating the hydraulic conductivity
// this function is obsolete and will never be called
void SOILtabletype::CalcKOfSE(int offset,double *k_sat_of_SE,double *SE_sl,double *theta_sl,double *thetasat,double *thetares)
{
//	double x = 1/0.0; // make sure an error is thrown when this function is called accidentally
	double kfrozen=this->getLnKMinFrozenSoil();
	double thetarel;
	double k_of_theta;
	for (int i=1;i<=this->levels;i++)
	{
		if (SE_sl)
		{
//			double ln_ksat=k_sat[i]<=0?-20.0:log(k_sat[i]);
			
			// contains now the logarithmic average fo the saturated hydraulic conductivity when all water is liquid and the saturated hydraulic conductivity when all water is ice
			// weighed by the relative water saturation SE
			
			// is done now with using a lookup table 
			// the commented section uses interpolations between sample points. It needs less sample points (100...1000) but is slightly slower than the direct adressing 
			// (total heat transfer model speed will be reduced by 5% when doing the 3 lookups by interpolating between sparse sample points
			// same algorithm like below but slower
			// version 1.0 explicitely solved (without pre-computed lookup-tables)
//			k_sat_of_SE[i] = exp(SE_sl[i]*ln_ksat+(1.0-SE_sl[i])*this->getLnKMinFrozenSoil());

			// version 1.0


			// this type of direct index addressing is dangerous but a little bit faster than the interpolation above, it works best for anz_T = 10000 (0.01°C steps in the array)
			// k_sat_of_SE[i] = ksat_SE[i][int(anz_T*SE_sl[i])];

			// version 1.1
//			double indexf=(SE_sl[i]*anz_T);
//			int indexlow=int(indexf);
//			double fraction = (indexf-indexlow);
//			if (indexlow < 0) indexlow=0; else if (indexlow >= anz_T) indexlow=anz_T-1;
//			int indexup=indexlow+1;
//			// is done now with a lookup table
//			k_sat_of_SE[i] = ksat_SE[i][indexup]*fraction + ksat_SE[i][indexlow]*(1.0-fraction);
			// 17 sek.

			// version 2.0: only SE and the normal Richards-aproach will be taken into account
//			k_sat_of_SE[i] = k_calc(i,theta_sl[i]*SE_sl[i]); // this will result in less dynamic k-values with temperature and in extremely small k values when theta is very low
			// 22 sek.

			// version 2.1: using both SE and Theta (the lookuptable ksat_SE is not used here since we have not only SE but also Theta as degree of freedom)
//////			thetarel=(theta_sl[i]-thetar[i])/(thetas[i]-thetar[i]);
			thetarel=((double)theta_sl[i])/(thetasat[i]);
			if (thetarel<=0) thetarel = 0.001; // to avoid undefined values
			else if (thetarel>1) thetarel = 1; 
			k_of_theta=(double)k_calc(i,theta_sl[i],1)*k_sat[i];
			if (k_of_theta<=0) k_of_theta=1e-23;
			double exponent=(1.0-(1.0-SE_sl[i])*thetarel)*log(k_of_theta)+(1.0-(double)SE_sl[i])*thetarel*kfrozen;
			k_sat_of_SE[i] = exp(exponent);
			// 25 sek.
		}
		else
		{
			k_sat_of_SE[i]=k_sat[i];
		}

	}
}

double SOILtabletype::GetTemperatureShift(int level,double actTheta)
{
	double Tshift=0;
	double indexf=anz_T*actTheta;
	int indexlow=int(indexf);
	if (indexlow<0) {
		indexlow=0; 
		Tshift = Ts_RelTheta[level][indexlow];
		return abs(Tshift);
	} else if (indexlow>=anz_T) {
		indexlow=anz_T; 
		Tshift = Ts_RelTheta[level][indexlow];
		return abs(Tshift);
	} else {
		double fraction = (indexf-indexlow);
		int indexup=indexlow+1;
		Tshift = Ts_RelTheta[level][indexup]*fraction + Ts_RelTheta[level][indexlow]*(1.0-fraction);
		return abs(Tshift);
	}
}

void SOILtabletype::CalcdEOverdT(double *dEOverdT_sl,double *temper_sl,double *theta_sl)
{
	if (temper_sl)
	{
		for (int i=1;i<=this->levels;i++)
		{
			if (temper_sl[i]>0)
			{
				dEOverdT_sl[i]=0.0; // alle ice thawed when soil temperature si > 0, so no phase change energy is needed
			}
			else if (temper_sl[i]<PhaseChngMinTemp)
			{
				dEOverdT_sl[i]=0.0; // below PhaseChngMinTemp(e.g. -10.0°C) we assume total ice saturation, so no phase change energy si needed here
			}
			else
			{
				// to account for unsaturated soils, an effective temperature is now calculated, see equ. 2.24.6 in documentation
				double t_eff = 0;
				double relTheta=(theta_sl[i]-this->thetar[i])/(thetas[i]-thetar[i]);
				t_eff=temper_sl[i]+GetTemperatureShift(i,relTheta);
				
				if (t_eff>=0)
				{
					dEOverdT_sl[i]=0.0;
				}
				else
				{
					double indexf = t_eff*this->anz_T/PhaseChngMinTemp;
					int indexlow=int(indexf);
					double fraction = (indexf-indexlow);
					if (indexlow < 0) indexlow=0; else if (indexlow >= anz_T) indexlow=anz_T-1;
					int indexup=indexlow+1;
					dEOverdT_sl[i] = dEOverdT[i][indexup]*fraction + dEOverdT[i][indexlow]*(1.0-fraction);
				
					// this type of direct index addressing is dangerous but a little bit faster than the interpolation above, it works best for anz_T = 10000 (0.01°C steps in the array)
	//				dEOverdT_sl[i] = dEOverdT[i][int(-(anz_T/10.0)*temper_sl[i])];
				}

			}
		}
	}
}

void SOILtabletype::CalcSEanddEOverdT(int offset,double* SE_sl, double *dEOverdT_sl, double *temper_sl, double factor, double *theta_sl, double *thetares, double *thetasat)
{
	if (temper_sl)
	{
		for (int i = 1; i <= this->levels; i++)
		{
			if (temper_sl[i+offset]>0)
			{
				SE_sl[i+offset] = 1.0; // alle ice thawed when soil temperature si > 0, so water saturation is 1.0
				dEOverdT_sl[i+offset] = 0.0; // alle ice thawed when soil temperature si > 0, so no phase change energy is needed
			}
			else if (temper_sl[i+offset] <= this->PhaseChngMinTemp)
			{
				SE_sl[i+offset] = 0.0; // below -10.0 °C we assume total ice saturation
				dEOverdT_sl[i+offset] = 0.0; // alle ice thawed when soil temperature si > 0, so no phase change energy is needed
			}
			else
			{
				// to account for unsaturated soils, an effective temperature is now calculated, see equ. 2.24.6 in documentation
				double t_eff = temper_sl[i+offset];
				double relTheta = (theta_sl[i+offset] - thetares[i]) / (thetasat[i] - thetares[i]);
				t_eff = temper_sl[i+offset] + GetTemperatureShift(i, relTheta);

//				if (offset>=0) // should only be corrected for soil layers, not for snow
//				{
//					t_eff = temper_sl[i+offset] + GetTemperatureShift(i, relTheta);
//				}

				if (t_eff >= 0)
				{
					SE_sl[i+offset] = 1.0; // alle ice thawed when soil temperature si > 0, so water saturation is 1.0
					dEOverdT_sl[i+offset] = 0.0;
				}
				else
				{
					double indexf = t_eff*this->anz_T / PhaseChngMinTemp;
					int indexlow = int(indexf);
					double fraction = (indexf - indexlow);
					if (indexlow < 0) indexlow = 0; else if (indexlow >= anz_T) indexlow = anz_T - 1;
					int indexup = indexlow + 1;
					SE_sl[i+offset] = SE_T[i][indexup] * fraction + SE_T[i][indexlow] * (1.0 - fraction);
					dEOverdT_sl[i+offset] = dEOverdT[i][indexup] * fraction + dEOverdT[i][indexlow] * (1.0 - fraction);
					// this type of direct index addressing is dangerous but a little bit faster than the interpolation above, it works best for anz_T = 10000 (0.01°C steps in the array
					// SE_sl[i] = SE_T[i][int(-(anz_T/10.0)*temper_sl[i])]; // range from 0 to 10000, since temper may vary from 0 to -10
					// this type of direct index addressing is dangerous but a little bit faster than the interpolation above, it works best for anz_T = 10000 (0.01°C steps in the array)
					// dEOverdT_sl[i] = dEOverdT[i][int(-(anz_T/10.0)*temper_sl[i])];

					// old without lookup table: 
					// SE_sl[i] = pow(double(1.0/(1.0+pow(double(this->getAlpha(i)*fabs(factor*temper_sl[i])),double(this->getN(i))))),double(this->getM(i)));
				}
			}
		}
	}
}

void SOILtabletype::CalcSE(double* SE_sl,double *temper_sl,double factor,double *theta_sl)
{
	if (temper_sl)
	{
		for (int i=1;i<=this->levels;i++)
		{
			if (temper_sl[i]>0)
			{
				SE_sl[i]=1.0; // alle ice thawed when soil temperature si > 0, so water saturation is 1.0
			}
			else if (temper_sl[i]<=this->PhaseChngMinTemp)
			{
				SE_sl[i]=0.0; // below -10.0 °C we assume total ice saturation
			}
			else
			{
				// to account for unsaturated soils, an effective temperature is now calculated, see equ. 2.24.6 in documentation
				double t_eff = 0;
				double relTheta=(theta_sl[i]-this->thetar[i])/(thetas[i]-thetar[i]);
				t_eff=temper_sl[i]+GetTemperatureShift(i,relTheta);
				if (t_eff>=0)
				{
					SE_sl[i]=1.0; // alle ice thawed when soil temperature si > 0, so water saturation is 1.0
				}
				else
				{
				
					double indexf = t_eff*this->anz_T/PhaseChngMinTemp;
					int indexlow=int(indexf);
					double fraction = (indexf-indexlow);
					if (indexlow < 0) indexlow=0; else if (indexlow >= anz_T) indexlow=anz_T-1;
					int indexup=indexlow+1;
					// is done now with a lookup table
					SE_sl[i] = SE_T[i][indexup]*fraction + SE_T[i][indexlow]*(1.0-fraction);
				
					// this type of direct index addressing is dangerous but a little bit faster than the interpolation above, it works best for anz_T = 10000 (0.01°C steps in the array
	//				SE_sl[i] = SE_T[i][int(-(anz_T/10.0)*temper_sl[i])]; // range from 0 to 10000, since temper may vary from 0 to -10
				
					// old without lookup table: 
					// SE_sl[i] = pow(double(1.0/(1.0+pow(double(this->getAlpha(i)*fabs(factor*temper_sl[i])),double(this->getN(i))))),double(this->getM(i)));
				}
			}
		}
	}
}

// calculating effective heat capacity for the given temperature and absolute and effective saturation, regarding the effect of phase change
void SOILtabletype::CalcHeatCapacity(int offset,double *temper_sl,double *heatcap_sl,double *SE_sl,double *dEOverdT_sl,double *theta_sl,double dTModRhoWater,double dTModCWater,double dTModCIce0,double dTModCIce20)
{
	double *drydensity=this->getDryDensity();
	double *dryheatcap=this->getDryHeatCapacity();
	for (int i=1;i<=levels; i++)
	{
		double CIceAct = dTModCIce0;
		if (temper_sl[i+offset]<-20)		{
			CIceAct = dTModCIce20;
		} else if (temper_sl[i+offset]<0)	{
			CIceAct = temper_sl[i+offset]/(-20.0) * dTModCIce20 + (1.0-temper_sl[i]/(-20.0))*dTModCIce0;
		}
		double temp1 = drydensity[i]*dryheatcap[i];
		double temp2 = SE_sl[i+offset]*dTModCWater + (1.0-SE_sl[i+offset])*CIceAct;
		heatcap_sl[i+offset] = temp1 + theta_sl[i+offset]*dTModRhoWater*temp2 + dEOverdT_sl[i+offset];
	}
}

/* gesaettigte hydraulische Leitfaehigkeit fuer die einzelnen Schichten des Bodens in m/s */
double* SOILtabletype::getK_sat()
{
	return k_sat;
}     

//#pragma optimize ("p",on)
double SOILtabletype::h_calc(int layer,double theta, double thetarex, double thetasex)
/* searching in the suction look up table for the corresponding value for the given parameter theta */
{
	const double *theta_h_layer = theta_h[layer];
	const double *h_tab_layer = h_tab[layer];
//	double theta_h_layer[31];
//	double h_tab_layer[31];
//	int sz = sizeof(double) * 31;
//	memmove(theta_h_layer, theta_h[layer], sz);
//	memmove(h_tab_layer, h_tab[layer], sz);

	// convert theta in case fo differently thetas and thetar than when creating the lookup table
	if ((thetasex > 0) && (thetarex > 0))
	{
		double rel_theta = (theta - thetarex) / (thetasex - thetarex);
		theta = this->thetar[layer] + (this->thetas[layer] - this->thetar[layer])*rel_theta;
	}

	if (theta<=theta_h_layer[anz_h])
	return h_tab_layer[anz_h];  /* if theta was below wilting point, maximal suction is returned */

	if (theta>=theta_h_layer[1])
	return 0.0;

	// find index by binary search
	unsigned short int min=1;
	unsigned short int max=anz_h;
//	unsigned short int offset = min + (max - min) / 2;
	unsigned short int offset = min + ((max - min) >> 1);

	while(min <= max)   {
		if(theta_h_layer[offset] <= theta) 	{
			max = offset-1;
		} else if(theta_h_layer[offset] > theta) {
			min = offset+1;
		} else {
			break;
		}
//		offset = min + (max - min) / 2;
		offset = min + ((max - min) >> 1);
	}
	offset=max;

	const double a=(theta-theta_h_layer[offset+1])/(theta_h_layer[offset]-theta_h_layer[offset+1]); /* a is the weight for the left (first) value, 1-a the weight for the second (right) value */
	return h_tab_layer[offset]*a + h_tab_layer[offset+1]*(1.0-a);   /* weighted average from the first and the second entry of the look up table */
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double SOILtabletype::k_calc(int layer,double theta,double SE, double thetarex, double thetasex)
// searching in the kf look up table for the corresponding value for the given parameter theta
{
	const double *theta_k_layer=theta_k[layer];
	const double *k_tab_layer=k_tab[layer];

//	double theta_k_layer[31];
//	double k_tab_layer[31];
//	int sz = sizeof(double) * 31;
//	memmove(theta_k_layer, theta_k[layer], sz);
//	memmove(k_tab_layer, k_tab[layer], sz);


	
	// convert theta in case fo differently thetas and thetar than when creating the lookup table
	if ((thetasex > 0) && (thetarex > 0))
	{
		double rel_theta = (theta - thetarex) / (thetasex - thetarex);
		// multiplied by SE 2023-10-06
		theta = this->thetar[layer] + (this->thetas[layer] - this->thetar[layer])*rel_theta*SE;
	}
	else
	{
		// added 2023-10-06
		theta = this->thetar[layer] + (theta - this->thetar[layer])*SE;
	}
	
	// added 2023-10-06
	if (theta < this->thetar[layer])
	{
		theta = this->thetar[layer]; // limit theta to the residual value, which is usually 0 in the lookup-table
	}
		
	// possible variation: 
	// if (SE < theta) theta = SE; // here, we set the minimum of SE or theta as theta in order to compute k-values for partly frozen soil
	// theta *= SE; // changed 2023-04-14: reduce the effective theta by factor SE (i.e. if no ice and SE=1, theta remains as is, if everything is frozen, theta converges to 0)
	// commented out 2023-10-06

	if (theta<theta_k_layer[anz_k])
		return k_tab_layer[anz_k];  // returning minimum conductivity for theta below wilting point
	if (theta>=theta_k_layer[1])
		return 1.0;

	unsigned short int min=1;
	unsigned short int max=anz_h;
//	unsigned short int offset = min + (max - min) / 2;
	unsigned short int offset = min + ((max - min) >> 1);

	while(min <= max)   {
		if(theta_k_layer[offset] <= theta) 	{
			max = offset-1;
		} else if(theta_k_layer[offset] > theta) {
			min = offset+1;
		} else {
			break;
		}
//		offset = min + (max - min) / 2;
		offset = min + ((max - min) >> 1);
	}
	offset=max;

	const double a=(theta-theta_k_layer[offset+1])/(theta_k_layer[offset]-theta_k_layer[offset+1]);/* a is the weight for the left (first) value, 1-a the weight for the second (right) value */
	return k_tab_layer[offset]*a + k_tab_layer[offset+1]*(1.0-a); /* weighted average from the first and the second entry of the look up table */

}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double SOILtabletype::get_theta(int layer,double h, double thetarex, double thetasex)
{ /* calculates from the h-theta-trable a theta-value for a given h-value (h=suction in m) */

	if (h>=h_tab[layer][anz_h])
	return theta_h[layer][anz_h];  /* if h is to large, theta is at or below wilting point */
	const double *theta_h_layer=theta_h[layer];
	const double *h_tab_layer=h_tab[layer];

	int min=1;
	int max=anz_h;
	int offset= min + (max - min) / 2;

	while(min <= max)   {
		if(h_tab_layer[offset] >= h) 	{
			max = offset-1;
		} else if(h_tab_layer[offset] < h) {
			min = offset+1;
		} else {
			break;
		}
		offset = min + (max - min) / 2;
	}
	offset=max;

	double a;
	if ((h_tab_layer[offset+1]-h_tab_layer[offset]) < 0.001)
	{
		a = 0;
	}
	else
	a=(h-h_tab_layer[offset])/(h_tab_layer[offset+1]-h_tab_layer[offset]);
	double res = theta_h_layer[offset+1]*a + theta_h_layer[offset]*(1.0-a);

	// convert theta in case fo differently thetas and thetar than when creating the lookup table
	if ((thetasex > 0) && (thetarex > 0))
	{
		double rel_theta = (res - this->thetar[layer]) / (this->thetas[layer] - this->thetar[layer]);
		res = thetarex + (thetasex - thetarex)*rel_theta;
		return res;
	}
	else
	{
		return res;
	}
}
//#pragma optimize ("p",off)

// initializing for the new soiltable layout with multiple horizons
void SOILtabletype::initialize(wasim::mystream *WSstream,std::string row,std::string sVersion)
{

	if (row.find("MultipleHorizons")!=std::string::npos)
	{
		//		method=2; // the new method: using an open number of horizons
	}
	else
	{
		(*WSstream) << "FATAL: soil table definition method unknown "<< row.c_str() << "[end-entry]";
		throw(29);
	}

	std::stringstream ss; 
	ss << row.c_str();
	std::string buffer;
	int index1;
	int index2;

	//	std::vector<std::string> words;

	try
	{
		ss >> Code >> Name;

		row = row.substr(row.find("{"));

		iNumHorizons=-1;
		levels=stripParameters(row,"layers",NULL,iNumHorizons);

		dzui   = new double[levels+1];
		k_sat  = new double[levels+1];
		thetas = new double[levels+1];
		thetar = new double[levels+1];
		psi_FC = new double[levels+1];
		m      = new double[levels+1];
		n      = new double[levels+1];
		alpha  = new double[levels+1];
		tau    = new double[levels+1];
		k_rec  = new double[levels+1];
		DryHeatCapacity	  = new double[levels+1];
		DryDensity		  = new double[levels+1];
		DryThermalConduct = new double[levels+1];


		stripParameters(row,"ksat",k_sat,iNumHorizons);
		stripParameters(row,"theta_sat",thetas,iNumHorizons);
		stripParameters(row, "theta_res", thetar, iNumHorizons);
		stripParameters(row,"alpha",alpha,iNumHorizons);
		stripParameters(row,"Par_n",n,iNumHorizons);
		stripParameters(row,"thickness",dzui,iNumHorizons);

		// defaultvalue for maximum evaporation depth: first numeric layers depth
		dEvapMaxDepth=dzui[1];

		if (row.find("isArcVegType") != std::string::npos)
		{
			stripSingleParameter(row, "isArcVegType", &dIsArcVegType);
			bIsArcVegType = (dIsArcVegType != 0);
		}
		else
		{
			bIsArcVegType = false;
		}


		if (row.find("EvapMaxDepth") != std::string::npos)
		{
			stripSingleParameter(row, "EvapMaxDepth", &dEvapMaxDepth);
		}

		if (row.find("ThicknessScaling") != std::string::npos)
		{
			stripSingleParameter(row, "ThicknessScaling", &dThicknessScaling);
			if (dThicknessScaling != 1.0)
			{
				(*WSstream) << "INFO: Soil layers for soil type " << Code << " will be scaled by a factor of " << dThicknessScaling << "[end-entry]";
			}
		}
		else
		{
			dThicknessScaling = 1.0;
		}
		dTotalThickness = 0;
		for (int i = 0; i <= levels; i++)
		{
			dzui[i] *= dThicknessScaling;
			if (i > 0) dTotalThickness += dzui[i];
		}
		
		// optional parameters for silting up model
		if (row.find("GrainSizeDist")!=std::string::npos)
		{
			stripGrainSizeParameters(row);
			if (dFractionSand<=0)
			{
				(*WSstream) << "FATAL: Grainsize distribution: fraction of sand must be > 0 \nsoil type: "<<Code<<"\n"<< row.c_str() << "[end-entry]";
				throw(110);
			}

			if ((dFractionSilt<0)||(dFractionClay<0)||(dFractionStones1 <0)||(dFractionStones2 <0)||(dFractionStones3 <0)||(dFractionStones4<0))
			{
				(*WSstream) << "FATAL: Grainsize distribution: all grain size fractions must be >= 0 (sand must be >0)\nsoil type: "<<Code<<"\n"<< row.c_str() << "[end-entry]";
				throw(110);
			}

			if (fabs(dFractionSand + dFractionSilt + dFractionClay + dFractionStones1 + dFractionStones2 + dFractionStones3 + dFractionStones4 -1.0)>1e-4)
			{
				(*WSstream) << "FATAL: Grainsize distribution: sum of all fractions (sand, silt, clay, stones1..4) must be 1.0\n soil type: "<<Code<<"\n"<< row.c_str() << "[end-entry]";
				throw(110);
			}
		}
		else
		{
			dFractionSand=0.3333;
			dFractionSilt=0.3333;
			dFractionClay=0.3334;
			dFractionStones1=0.0;
			dFractionStones2=0.0;
			dFractionStones3=0.0;
			dFractionStones4=0.0;
		}

		// optional parameters for silting up model with method 3 (dynamic expressions)
		for (int su_par=0;su_par<=10;su_par++)
		{
			char parname[20];
			sprintf(parname,"SU_PAR%.02d",su_par);
			SU_PAR[su_par]=0.0; // initial value
			// if the parameter si defined, the value will be overwritten now
			if (row.find(parname)!=std::string::npos)
			{
				stripSingleParameter(row,std::string(parname),&(SU_PAR[su_par]));
			}
		}

		// optional parameter psi_FC
		useGlobalDefaultPsiFC = false; // set false by default. Can be changed externally after initializing all soil table types, if useLocalDefaultPsiFC is true
		if (row.find("psi_FC") != std::string::npos)
		{
			stripParameters(row, "psi_FC", psi_FC, iNumHorizons);
			useLocalDefaultPsiFC = false;
		}
		else
		{
			for (int i = 0; i <= levels; i++)
			{
				psi_FC[i] = 200; // default is 200 cm (negative pressure, i.e. suction)
			}
			useLocalDefaultPsiFC = true;
		}

		// optional parameter tau
		if (row.find("Par_tau")!=std::string::npos)
		{
			stripParameters(row,"Par_tau",tau,iNumHorizons);
		}
		else
		{
			int i;
			for (i=0;i<=levels;i++) 
			{
				tau[i]=0.5;
			}
		}

		// optional parameter k_recession
		if (row.find("k_recession")!=std::string::npos)
		{
			stripParameters(row,"k_recession",k_rec,iNumHorizons);
		}
		else
		{
			int i;
			for (i=1;i<=levels;i++) 
			{
				k_rec[i]=1.0;
			}
		}

		if (row.find("PMacroThresh")!=std::string::npos)		{
			stripSingleParameter(row,"PMacroThresh",&dPMacroThreshold);
		}		else		{
			dPMacroThreshold=1000; // 1000 mm/h means no macro pore runoff at all probably
		}

		if (row.find("MacroCapacity")!=std::string::npos)		{
			stripSingleParameter(row,"MacroCapacity",&dPMacroCapacity);		
		}		else		{
			dPMacroCapacity=0.0; // 1000 mm/h means no macro pore runoff at all probably
		}

		if (row.find("CapacityRedu")!=std::string::npos)		{
			stripSingleParameter(row,"CapacityRedu",&dReductionPerMeter);		
		}		else		{
			dReductionPerMeter=1.0; // 1000 mm/h means no macro pore runoff at all probably
		}

		if (row.find("MacroDepth")!=std::string::npos)		{
			stripSingleParameter(row,"MacroDepth",&dMacroDepth);
		}		else		{
			dMacroDepth=1.0; // 1000 mm/h means no macro pore runoff at all probably
		}

		if (row.find("DryHeatCapacity")!=std::string::npos)		{
			stripParameters(row,"DryHeatCapacity",DryHeatCapacity,iNumHorizons);
//			stripSingleParameter(row,"DryHeatCapacity",&DryHeatCapacity);
		} // default was set when calling the constructor; global defaults will be set for each non-assigned soil type after initialization of soil model
		else
		{
			for (int l=0;l<=levels;l++) DryHeatCapacity[l]=-1;
		}

		if (row.find("DryDensity")!=std::string::npos)		{
			stripParameters(row,"DryDensity",DryDensity,iNumHorizons);
//			stripSingleParameter(row,"DryDensity",&DryDensity);
		} // default was set when calling the constructor; global defaults will be set for each non-assigned soil type after initialization of soil model
		else
		{
			for (int l=0;l<=levels;l++) DryDensity[l]=-1;
		}

		if (row.find("DryThermalConduct")!=std::string::npos)		{
			stripParameters(row,"DryThermalConduct",DryThermalConduct,iNumHorizons);
//			stripSingleParameter(row,"DryThermalConduct",&DryThermalConduct);
		} // default was set when calling the constructor; global defaults will be set for each non-assigned soil type after initialization of soil model
		else
		{
			for (int l=0;l<=levels;l++) DryThermalConduct[l]=-1;
		}

		if (row.find("KMinFrozenSoil")!=std::string::npos)		{
			stripSingleParameter(row,"KMinFrozenSoil",&KMinFrozenSoil);
			if (KMinFrozenSoil>0) lnKMinFrozenSoil=log(KMinFrozenSoil); else lnKMinFrozenSoil=-27.631; // natural logarithm -27.631 equals 1e-12 
		} // default was set when calling the constructor; global defaults will be set for each non-assigned soil type after initialization of soil model

		// converting mm/h into m/s
		dPMacroThreshold=dPMacroThreshold/3.6e6; 
		dPMacroCapacity=dPMacroCapacity/3.6e6;
		
		// should be read later, too
		if (row.find("maxratio")!=std::string::npos)		{
			stripSingleParameter(row,"maxratio",&maxratio);
		}		else		{
			maxratio=1000; // maximum proportion of effektive k-values between two layers is limited by this factor, using the higher value as reference
		}

		double dnFK=0;
                double dmSB=0;
                double dk=0;
                double dpsi=0;
		// for compatibility with topmodel version...
		if (row.find("FCap")!=std::string::npos)		{
			stripSingleParameter(row,"FCap",&dnFK);
		}		else		{
			dnFK=25; // nutzbare Feldkapazitaet/useable field capacity in Vol%
		}

		if (row.find("mSB")!=std::string::npos)		{
			stripSingleParameter(row,"mSB",&dmSB);
		}		else		{
			dmSB=35; // maximum soil storage in Vol%
		}
		
		if (row.find("ksat_topmodel")!=std::string::npos)		{
			stripSingleParameter(row,"ksat_topmodel",&dk);
		}		else		{
			k=1e-5; // hydraulic conductivity in m/h
		}

		if (row.find("suction")!=std::string::npos)		{
			stripSingleParameter(row,"suction",&dpsi);
		}		else		{
			dpsi=300; // suction at infiltration front
		}

		nFK=(double)dnFK;
		mSB=(double)dmSB;
		k=(double)dk;
		psi=(double)dpsi;

		drain=(mSB-nFK)*0.01;

		if (row.find("Par_m")!=std::string::npos)		{
			stripParameters(row,"Par_m",m,iNumHorizons);
		}
		else
		{
			// generating m from alpha
			for (int k=0;k<=levels;k++)
			{
				if (n[k]!=0)
				{
					m[k] = 1.0 - 1.0/n[k];
				}
				else
				m[k]=0;
			}
		}

		if (row.find("MaxInterflowDepth")!=std::string::npos)		{
			stripSingleParameter(row,"MaxInterflowDepth",&dMaxInterflowDepth);		
		}		else		{
			dMaxInterflowDepth=2.0; // default value for backwards compatibility
		}

		anz_k=30;
		anz_h=30;

		typedef double* dp ;
		theta_k = new dp[levels+1];
		k_tab   = new dp[levels+1];
		theta_h = new dp[levels+1];
		h_tab   = new dp[levels+1];

// set in constructor already		anz_T=10000; // 0.01°C between 0 and -10°C
		SE_T      = new dp[levels+1]; // lookup table for relative saturation as a function of temperature (per soil layer!)
		dEOverdT  = new dp[levels+1]; // lookup table for internal energy change (per soil layer!)
		Ts_RelTheta = new dp[levels+1]; // lookup table for temperature shifts
//		ksat_SE   = new dp[levels+1]; // lookup table for internal energy change (per soil layer!)

		for (int layer=0;layer<=levels;layer++)
		{
					// calling tab generation here
					LookUpTableGeneration(layer);
					// allocation is done here in order to make sure the values exist (because they will be deleted from this class, as well, and then no query is made if all vectors were allocated correctly
					// even for index 0!
					SE_T[layer] = new double[anz_T+1];
					dEOverdT[layer] = new double[anz_T+1];
					Ts_RelTheta[layer]=new double[anz_T+1];
//					ksat_SE[layer] = new double[anz_T+1];
					// will be called from unsatzon-constructor since we know the factor only there
					// LookUpTableGenerationTMod(layer);
		}

		//		double depth=0;
		//		k_abnahme=0.5;
		//        for (int ak=2;ak<=levels;ak++)
		//        {
		//			depth+=dzui[ak]/2.0+dzui[ak-1]/2.0;
		//			if (depth>2) break;
		//			k_sat[ak] = k_sat[1]*pow(k_abnahme,depth);
		//        }

		// calculating recession in saturated hydraulic conductivity if this parameter was given
		int lev=0;
		double depthSumTotal=0;
		for (int i=1;i<=iNumHorizons;i++)
		{
			int start=lev+1;
			int stop=lev+LevelPerHorizon[i];
			double depthSum=0;
			for (int k=start;k<=stop;k++)
			{
				if (k>start)
				{
					depthSum+=dzui[k]/2.0+dzui[k-1]/2.0;
					depthSumTotal+=dzui[k]/2.0+dzui[k-1]/2.0;
				}
				else
				{
					depthSum=dzui[start]/2.0;
					if (k==1)
					{
						depthSumTotal=dzui[start]/2.0;
					}
					else
					{
						depthSumTotal+=dzui[k]/2.0+dzui[k-1]/2.0;
					}
				}
				if (depthSumTotal>dMaxInterflowDepth) break;
				k_sat[k] = k_sat[k]*pow(k_rec[i],depthSum);
				//				(*WSstream) << "ksat["<<k<<"]: " << k_sat[k] <<"[end-entry]";
				lev++;
			}

		}


	}
	catch (...)
	{
		(*WSstream) << "FATAL: Soil table code "<<Code<<": some values could not be read\n" << row.c_str() << "[end-entry]";
		throw(29);
	}
}

void SOILtabletype::stripSingleParameter(std::string &row,std::string sWhatParam,double *param)
{
	int index1,index2;
	std::string buffer;

	buffer=row.substr(row.find(sWhatParam));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> *param;
}

void SOILtabletype::stripGrainSizeParameters(std::string &row)
{
	int index1,index2;
	std::string buffer;

	buffer=row.substr(row.find("GrainSizeDist"));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> dFractionSand >> dFractionSilt >> dFractionClay >> dFractionStones1 >> dFractionStones2 >> dFractionStones3 >> dFractionStones4;

}

int SOILtabletype::stripParameters(std::string &row,std::string sWhatParam,double *param,int &horizons)
{
	if (horizons<=0) 
	{	// read parameters from an unknown number of arguments and return horizons 
		int index1,index2;
		std::string buffer;
		std::vector<std::string> words;

		buffer=row.substr(row.rfind(sWhatParam));
		index1 = buffer.find("=")+1;
		index2 = buffer.find(";");
		buffer=buffer.substr(index1,index2-index1);
		words.clear();
		su.StringSplit(buffer,"\t ;,",words,false);
		iNumHorizons=0;
		int i;
		for (i=0;i<words.size();i++)
		{
			if (words[i].size()>0)
			{
				iNumHorizons++;
			}
		}
		LevelPerHorizon = new int[iNumHorizons+1];
		LevelPerHorizon[0]=0;
		std::stringstream ss;
		ss << buffer.c_str();
		for (i=1;i<=horizons;i++)
		{
			ss >> LevelPerHorizon[i];
			LevelPerHorizon[0]+=LevelPerHorizon[i];
		}
		horizons = iNumHorizons;
		return LevelPerHorizon[0];
	}
	else
	{ // read the other parameters with a known number of read-arguments and copy each "LevelPerHorizon" times
		//		param[0]=param[maxIndex]; // the zero-element is always holding the number of following parameter values
		double tempParam;
		int index1,index2,lev;
		std::string buffer;
		buffer=row.substr(row.rfind(sWhatParam));
		index1 = buffer.find("=")+1;
		index2 = buffer.find(";");
		buffer=buffer.substr(index1,index2-index1);
		std::stringstream ss;
		ss << buffer.c_str();
		lev=0;
		for (int i=1;i<=horizons;i++)
		{
			// if there are not sufficient parameters provided, the list will be filled with the last known values
			if (!ss.eof())	ss >> tempParam;
			int start=lev+1;
			int stop=lev+LevelPerHorizon[i];
			for (int k=start;k<=stop;k++)
			{
				param[k]=tempParam;
				lev++;
			}

		}
		param[0]=param[1];
		return lev;
	}
}

void SOILtabletype::setPsiFCDefault(double psiFCdefault, int usePsiFCglobalForUnassignedOnly)
{
	// only if the actual soil table entry has not assigned any psi_FC values yet OR if all values should be overwritten because usePsiFCglobalForUnassignedOnly is not set (i.e. equals 0, otherwise the assigned entries will be jumped over)
	if (((usePsiFCglobalForUnassignedOnly == 1) && (useLocalDefaultPsiFC == true)) || (usePsiFCglobalForUnassignedOnly == 0))
	{
		for (int i = 0; i <= levels; i++)
		{
			psi_FC[i] = psiFCdefault;
		}
	}
	
}

MultiLayerLUSEType::MultiLayerLUSEType():
CommonTableType(0),
iNumLayers(1),
iNumFractions(0),
iNumBiomass(0),
iNumRetrans(0),
LULayers(NULL),
ArcVegLUFractions(NULL), 
ArcVegBiomass(NULL),
ArcVegRetrans(NULL),
OriLULayers(NULL),
K_lambert(0.3),
LAI_scale(10),
OriCode(0)
{
}

MultiLayerLUSEType::MultiLayerLUSEType(int code):
CommonTableType(code),
iNumLayers(1),
iNumFractions(0),
iNumBiomass(0),
iNumRetrans(0),
LULayers(NULL),
ArcVegLUFractions(NULL),
ArcVegBiomass(NULL),
ArcVegRetrans(NULL),
OriLULayers(NULL),
K_lambert(0.3),
LAI_scale(10),
OriCode(code)
{
}

MultiLayerLUSEType::~MultiLayerLUSEType()
{
	if (LULayers)
	{
		delete[] LULayers;
		LULayers=NULL;
	}
	if (ArcVegLUFractions)
	{
		delete[] ArcVegLUFractions;
		ArcVegLUFractions = NULL;
	}
	if (ArcVegBiomass)
	{
		delete[] ArcVegBiomass;
		ArcVegBiomass = NULL;
	}
	if (ArcVegRetrans)
	{
		delete[] ArcVegRetrans;
		ArcVegRetrans = NULL;
	}
	if (OriLULayers)
	{
		delete[] OriLULayers;
		OriLULayers=NULL;
	}
}

void MultiLayerLUSEType::initialize(wasim::mystream *WSstream,std::string row,int version)
{
	if (row.find("Landuse_Layers")==std::string::npos)
	{
		(*WSstream) << "FATAL: multiple landuse table: parameter Landuse_Layers not found "<< row.c_str() << "[end-entry]";
		throw(112);
	}

	std::stringstream ss; 
	ss << row.c_str();
	std::string buffer;
	int index1;
	int index2;

	try
	{
		ss >> Code >> Name;
		OriCode=Code;
		stripParameters(row,"Landuse_Layers",&LULayers,iNumLayers);
		OriLULayers = new int[iNumLayers+1];
		for (int i=0;i<=iNumLayers;i++)
		{
			OriLULayers[i]=LULayers[i];
		}

		// optionally read ArcVegLUFractions
		if (row.find("ArcVeg_Initial_Fractions") == std::string::npos)
		{
			ArcVegLUFractions = NULL; // if later this variable is found to be 0 and ArcVeg is active, an error message will be thrown and the program will be aborted
		}
		else
		{
			stripParameters(row, "ArcVeg_Initial_Fractions", &ArcVegLUFractions, iNumFractions);
		}

		// optionally read ArcVegLUBiomass
		if (row.find("ArcVeg_Initial_Biomass") == std::string::npos)
		{
			ArcVegBiomass = NULL; // if later this variable is found to be 0 and ArcVeg is active, an error message will be thrown and the program will be aborted
		}
		else
		{
			stripParameters(row, "ArcVeg_Initial_Biomass", &ArcVegBiomass, iNumBiomass);
		}

		// optionally read ArcVeg_Initial_Retrans
		if (row.find("ArcVeg_Initial_Retrans") == std::string::npos)
		{
			ArcVegRetrans = NULL; // if later this variable is found to be 0 and ArcVeg is active, an error message will be thrown and the program will be aborted
		}
		else
		{
			stripParameters(row, "ArcVeg_Initial_Retrans", &ArcVegRetrans, iNumRetrans);
		}

		// optionally read extinction coefficient
		if (row.find("k_extinct") == std::string::npos)
		{
			K_lambert = 0.3;
		}
		else
		{
			buffer = row.substr(row.find("k_extinct"));
			index1 = buffer.find("=") + 1;
			index2 = buffer.find(";");
			buffer = buffer.substr(index1, index2 - index1);
			K_lambert = atof(buffer.c_str());
		}

		// optionally read LAI-scaling parameter for aerodynamic resistencies in lower layers
		if (row.find("LAI_scale")==std::string::npos)
		{
			LAI_scale  = 10;
		}
		else
		{
			buffer=row.substr(row.find("LAI_scale"));
			index1 = buffer.find("=")+1;
			index2 = buffer.find(";");
			buffer=buffer.substr(index1,index2-index1);
			LAI_scale = atof( buffer.c_str());
		}


	}
	catch (...)
	{
		(*WSstream) << "FATAL: multilayer_landuse_table code "<<Code<<": some parameters could not be read\n" << row.c_str() << "[end-entry]";
		throw(110);
	}
}

void MultiLayerLUSEType::stripParameters(std::string &row,std::string sWhatParam,int **param,int &NumLayers)
{
	int i,index1,index2;
	std::string buffer;
	buffer=row.substr(row.find(sWhatParam));
	index1 = buffer.find("=")+1;
	index2 = buffer.find(";");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();

	std::vector<std::string> words;
	std::vector<int> vi;
	int temp;
	words.clear();
	su.StringSplit(buffer,",;",words,false);
	
	for (i=0;i<words.size();i++) 
	{
		std::stringstream ss2;
		ss2 << words[i].c_str() << " ";
		ss2 >> temp;
		if (ss2.good())
		{
			vi.push_back(temp);
		}
		ss2.clear();
	}
	*param = new int[vi.size()+1];
	for (i=0;i<vi.size();i++)
	{
		(*param)[i+1]=vi[i];
	}
	(*param)[0]=vi.size();
	NumLayers=(*param)[0];
	vi.clear();
	words.clear();
}

void MultiLayerLUSEType::stripParameters(std::string& row, std::string sWhatParam, double** param, int& NumLayers)
{
	int i, index1, index2;
	std::string buffer;
	buffer = row.substr(row.find(sWhatParam));
	index1 = buffer.find("=") + 1;
	index2 = buffer.find(";");
	buffer = buffer.substr(index1, index2 - index1);
	std::stringstream ss;
	ss << buffer.c_str();

	std::vector<std::string> words;
	std::vector<double> vi;
	double temp;
	words.clear();
	su.StringSplit(buffer, ",;", words, false);

	for (i = 0; i < words.size(); i++)
	{
		std::stringstream ss2;
		ss2 << words[i].c_str() << " ";
		ss2 >> temp;
		if (ss2.good())
		{
			vi.push_back(temp);
		}
		ss2.clear();
	}
	*param = new double[vi.size() + 1];
	for (i = 0; i < vi.size(); i++)
	{
		(*param)[i + 1] = vi[i];
	}
	(*param)[0] = vi.size();
	NumLayers = vi.size();
	vi.clear();
	words.clear();
}

CommonTableType::~CommonTableType() 
{
	if (Name) 
	{
		delete[] Name;
		Name=NULL;
	}
}


