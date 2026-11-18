
/*
(c) 1995-2015 by Joerg Schulla
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include "teilmod.h"
#include "stringp.h"
#include "modconst.h"
#include "modef.h"
#include "wasim.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

void meteointerpol::set_warning_level(int warnlevel)
{
	warning_level=warnlevel;
}

void meteointerpol::set_hour_offset(double p_hour_offset)
{
	hour_offset=p_hour_offset;
}

void meteointerpol::getdatetime(long *me,long *se,long *tg,long *mt,long *jr)
{
	/* Die uebergebenen Parameter werden mit dem aktuell in der
	Instanz gueltigen Zeit (jahr,Monat,tag und stunde) geladen.
	So koennen andere Teilmodelle mit den fuer die Interpolation gueltigen Zeiten
	versorgt werden */
	*me = minute;
	*se = stunde;
	*tg = tag;
	*mt = monat;
	*jr = jahr;

}


void meteointerpol::parametereingabe(void)
{
	int i,j,s;

	if ((steuerdatei=fopen(steuerdateiname,"rb"))==NULL)
	{
		(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
		throw(28);
	}

	if (index(instanzname," ")>0)  instanzname = substr(instanzname,1,index(instanzname," ")-1);
	sprintf(tstr,"[%s]",instanzname);

	/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
	findeabschnitt(tstr);
	fscanfln(steuerdatei,"%d",&methode,"1 [idw], 2 [regress], [3] [idw+regr], 4 [thiessen] 5 [bilinear] 6 [BIGRES] 7 [bicubic spline] 8 [SPLINERES] 9 [ReadGrid] 10 [EDR from station data] 11 [IDW + EDR from station data] 12 [Thiessen combines with elevation gradient] for ",tstr);
	if ((methode==6)||(methode==8)) statfilecnt=2; else statfilecnt=1;
	for (s=0;s<statfilecnt;s++) 
	{
		numAdditionalColumns[s]=0;
		anetzstat[s]           = newstring(_stringsize_+1,1);
		AdditionalColumns[s]   = newstring(_stringsize_+1,1);
	}

	//    fscanfln(steuerdatei,"%s%s",(anetzstat[0]),AdditionalColumns,"name of the station file in ",tstr);
	for (s=0;s<statfilecnt;s++) 
	{
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		if (sscanf(tstr,"%s",&(anetzstat[s][0]))!=1)
		{
			(*WSstream) << "FATAL: error reading station data file name nr. "<<s<<" for interpolation of "<<instanzname<<"[end-entry]";
			throw(49);
		}
		strcpy(AdditionalColumns[s],tstr);
		if (strstr(AdditionalColumns[s],"AdditionalColumns=")) /* if a number for additional columns after the date columns is given, this number is claused in a AdditionalColumns=<n> statement */
		{
			sscanf((char*)(strstr(AdditionalColumns[s],"AdditionalColumns=")+18),"%d",&(numAdditionalColumns[s]));
		}
		std::string parameter(tstr);
		if (parameter.find("xmin") != std::string::npos) {
			std::stringstream ss;
			parameter = parameter.substr(parameter.find("xmin") + 1);
			parameter = parameter.substr(parameter.find("=") + 1);
			ss << parameter.c_str();
			ss >> xmin;
		}
		parameter = std::string(tstr);
		if (parameter.find("xmax") != std::string::npos) {
			std::stringstream ss;
			parameter = parameter.substr(parameter.find("xmax") + 1);
			parameter = parameter.substr(parameter.find("=") + 1);
			ss << parameter.c_str();
			ss >> xmax;
		}
		parameter = std::string(tstr);
		if (parameter.find("ymin") != std::string::npos) {
			std::stringstream ss;
			parameter = parameter.substr(parameter.find("ymin") + 1);
			parameter = parameter.substr(parameter.find("=") + 1);
			ss << parameter.c_str();
			ss >> ymin;
		}
		parameter = std::string(tstr);
		if (parameter.find("ymax") != std::string::npos) {
			std::stringstream ss;
			parameter = parameter.substr(parameter.find("ymax") + 1);
			parameter = parameter.substr(parameter.find("=") + 1);
			ss << parameter.c_str();
			ss >> ymax;
		}

	}
	if ((methode == 9) || (methode == 10) || (methode == 11) || (methode == 12))
	{
		if (methode == 9) {
			// read mapping factors
			tstr = fgets(tstr, _stringsize_, steuerdatei);
			int numread = sscanf(tstr, "%lf%lf%d%lf%lf", &gmf1, &gmx1, &gmcut, &gmf2, &gmx2);
			if (numread<5)
			{
				(*WSstream) << "FATAL: error reading the required 5 GIF-Mapping parameters f1 (default=8.916787e-5), x1 (0.0767528404), threshold (127), f2 (default=6.00e-006), x2 (0.098) for interpolation of " << instanzname << "[end-entry]";
				throw(49);
			}
		}
		else if ((methode == 10) || (methode == 11)) {
			// read regression parameters
			tstr=fgets(tstr,_stringsize_,steuerdatei);
			int numread=sscanf(tstr,"%lf%lf%lf%d%lf%lf",&lowerinversionfix,&upperinversionfix,&tolerance,&overlap,&clusterlimit,&lapserate);
			if (numread<5)
			{
				(*WSstream) << "FATAL: error reading the required 5 or 6 regression parameters <lower inversion> <upper inversion> <tolerance> <overlap> <clusterrange> [<lapse rate above upper inversion>] for interpolation of "<<instanzname<<"[end-entry]";
				throw(49);
			}
			else
			{
				if (numread==5) 
				{ // set lapserate to nodata, otherwise the readin value will be used
					lapserate=-9999; // so we will check the value of lapserate and decide dynamically, if the other regression below the upper layer should be still valid or if we apply the new lapserate based on the value at the layer boundary (above 
				}
			}
		} else	{
			// read elevation gradient as a single double value
			// sprintf(datatype,"P-type");
			tstr=fgets(tstr,_stringsize_,steuerdatei);
			int numread=sscanf(tstr,"%lf%lf%lf%s",&lapserate,&upperinversionfix,&upperlapserate,datatype);
			if ((numread<1)||(numread==2))
			{
				(*WSstream) << "FATAL: error reading the required parameter \"lapse rate\" and/or optional parameters inversion_hight and fix upper lapse rate for interpolation of "<<instanzname<<"[end-entry]";
				throw(49);
			}
			else
			{
				if (numread==1) 
				{ // set lapserate to nodata, otherwise the readin value will be used
					upperinversionfix=-9999; // so we will check the value of lapserate and decide dynamically, if the other regression below the upper layer should be still valid or if we apply the new lapserate based on the value at the layer boundary (above 
					upperlapserate=-9999; // so we will check the value of lapserate and decide dynamically, if the other regression below the upper layer should be still valid or if we apply the new lapserate based on the value at the layer boundary (above 
				}
			}
			if (strstr(datatype,"P-type"))
			{
				iDataType=PTYPE;
			}
			else
			{
				iDataType=TTYPE;
			}
		}
	}
	else
	{
		for (i=1;i<=regiomax;i++)
		fscanfln(steuerdatei,"%s",dateinamekoef[i],"name of the regression file in ",tstr);
	}
	fscanfln(steuerdatei,"%s",ergebnisgridname,"name of the result grid in ",tstr);
	fscanfln(steuerdatei,"%d",&writegrid,"grid write code for result grid in ",tstr);
	fscanfln(steuerdatei,"%lf",&kor,"scaling factor for ",tstr);
	fscanfln(steuerdatei,"%s%d",zonestatname,&zonestatcode,"name of the statistic result file in ",tstr);
	fscanfln(steuerdatei,"%lf",&fehlwert,"nodata value in ",tstr);
	fehlwerti = (int)fehlwert;

	fscanfln(steuerdatei,"%lf",&wichtung  ,"idw-weight in ",tstr);
	fscanfln(steuerdatei,"%lf",&danteil   ,"fraction of idw on weighted interpolation (method 3)",tstr);
	fscanfln(steuerdatei,"%lf",&maxentf   ,"maximum distance in ",tstr);
	fscanfln(steuerdatei,"%lf",&anisoneig ,"slope of the main axis of the anisotropy ellipse in ",tstr);
	fscanfln(steuerdatei,"%lf",&aniso     ,"anisotropy b/a in ",tstr);
	fscanfln(steuerdatei,"%lf",&grenzwertu,"lower limit of results in ",tstr);
	fscanfln(steuerdatei,"%lf",&fixu      ,"upper limit of results in ",tstr);
	fscanfln(steuerdatei,"%lf",&grenzwerto,"lower replacement for limit-exceeding results in ",tstr);
	fscanfln(steuerdatei,"%lf",&fixo      ,"upper replacement for limit-exceeding results in ",tstr);
	fscanfln(steuerdatei,"%d",&szenuse   ,"SzenUse in ",tstr);

	if (szenuse==1)
	{
		fscanfln(steuerdatei,"%d",&szenmeth  ,"SzenMeth in ",tstr);
		if ((szenmeth==5)||(szenmeth==6)) /* monthly grids are used for scaling the results */
		{
			scengridname=(char**)calloc(13,sizeof(char*));
			scengrididentifier=(char*)calloc(200,sizeof(char));
			for (i=1;i<=12;i++)
			{
				scengridname[i]=(char*)calloc(_stringsize_,sizeof(char));
				if ((j=fscanf(steuerdatei,"%s",scengridname[i])) != 1) 
				{
					(*WSstream) << "FATAL: error in control file: scenario grid for month "<<i<<" in "
					<<instanzname<<" not found!\n"
					"There have to be exactly 12 filenames after the sceario application mode\n"
					"referring to each of the 12 months if using scenario method 5 or 6.\n"
					"However, if methods 0..4 are used, NO FILENAMES are allowed there\n"
					" (for backwards compatibility reasons)" << "[end-entry]";
					throw(17);
				}
				fgets(tstr,_stringsize_,steuerdatei);
			}
		}
		else
		{
			fscanfln(steuerdatei,"%d",&szennum   ,"SzenNum in ",tstr);
			szenario = matrix(1,szennum,1,15);
			for (i=1;i<=szennum;i++)
			{
				for (j=1;j<=14;j++)
				{
					if (fscanf(steuerdatei,"%lf",&(szenario[i][j]))!=1)
					{
						(*WSstream) << "FATAL: error! scenario-data row "<<i
						<<", column "<<j<<" in "<<instanzname<<": read error!" << "[end-entry]";
						throw(48);
					}
					if ((szenmeth==3)&&(j>2)) /* prozentuale Aenderung */
					szenario[i][j] = 1.0 + szenario[i][j]/100.0;
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
			}
		}
	}

	if (wichtung==2) maxentf = pow(maxentf,wichtung);

	if ((fclose(steuerdatei))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
		throw(11);
	}
}

void meteointerpol::stationswerte(FILE *statfile,char *statfilename,
int numAddCol,bool *useMinutes,int *anz,
double **h,double **xk,double **yk,
double **dist,double **v)
{
	int i;
	double wert;
	(*WSstream) << "analyzing station data for interpolator "<<instanzname<<" from "<<statfilename<<" ..."<<"[end-entry]";

	std::ifstream tfile(statfilename); 
	std::string line;
	getline(tfile, line); // overread first line
	getline(tfile, line); // read the line that will be analyzed later by old C commands
	int linelength = line.size()+100;
	tfile.close();

	if (linelength > _meteosatzlaenge_)
	{
		if (meteosatz) free(meteosatz);
		meteosatz = (char*)calloc(linelength, sizeof(char));
		t_meteosatz = meteosatz;
	}
	else
	{
		linelength = _meteosatzlaenge_;
	}

	if ((statfile = fopen(statfilename,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: \nerror opening "<<statfilename<<"[end-entry]";
		throw(49);
	}

	do
	{ /* ueberlesen des ersten Satzes, der die Datenart kennzeichnet */
		meteosatz=fgets(meteosatz, linelength,statfile);
		i=strlen(meteosatz); if (i==0) i=1;
	} 
	while (meteosatz[i-1]!='\n'); /* wenn der reservierte Speicher fuer den Satz zu kurz war, weiter lesen */

	/* hier stehen die Hoehen drin */
	meteosatz=fgets(meteosatz, linelength,statfile);
	(*anz)=0;
	/* spaltenanalyse anhand der Hoehenangaben */
	/* Weglesen der ersten 4 Spalten oder 5 fuer das Datum */
	std::string row(meteosatz);
	
	int idx = 0;
	int len = 0;
	// row = row.find_first_not_of("YMDH \t"); will not work if additional columns start with other values than YMDH
	for (int i = 1; i <= 4; i++)
	{
		idx = row.find_first_not_of("YMDH"); 
		len = row.length() - idx;
		// optional: row.find_first_of("+-0123456789.")
		row = row.substr(idx, len);
	
		idx = row.find_first_not_of(" \t");
		len = row.length() - idx;
		// optional: row.find_first_of("+-0123456789.")
		row = row.substr(idx, len);
	}

	if ((row.find("mm") == 0) && (row.find_first_of(" \t") == 2))
	{
		(*useMinutes) = true;
		idx = row.find_first_not_of("m \t"); // looking for the first numeric value
		len = row.length() - idx;
		row = row.substr(idx, len);
	}

	for (int i = 1; i <= numAddCol; i++)
	{
		idx = row.find_first_of(" \t");
		len = row.length() - idx;
		// optional: row.find_first_of("+-0123456789.")
		row = row.substr(idx, len);

		idx = row.find_first_not_of(" \t");
		len = row.length() - idx;
		// optional: row.find_first_of("+-0123456789.")
		row = row.substr(idx, len);
	}

	row = row + std::string("\t");
	while (row.find_first_not_of(" \t\n\r") == 0)
	{
		(*anz)++;
		row = row.substr(row.find_first_of(" \t"), row.length() - row.find_first_of(" \t"));
		while (row.find_first_of(" \t") == 0) row = row.substr(1, row.length() - 1);
	}


	/* fuer meteosatz reicht auch ein kuerzerer Satz aus */
	//  meteosatz=t_meteosatz;
	if (meteosatz) free(meteosatz); 
	meteosatz=(char*)calloc((20*((*anz)+5+numAddCol))> linelength ?20*((*anz)+5+numAddCol): linelength,sizeof(char));
	t_meteosatz=meteosatz;

	/* nun belegen der arrays mit der richtigen Anzahlvon Feldern */
	*h  = (double *)calloc((*anz)+1,sizeof(double)); if (!*h) speichererror(1);
	*xk = (double *)calloc((*anz)+1,sizeof(double)); if (!*xk) speichererror(1);
	*yk = (double *)calloc((*anz)+1,sizeof(double)); if (!*yk) speichererror(1);
	*dist   = (double *)calloc((*anz)+1,sizeof(double)); if (!*dist) speichererror(1);
	*v  = (double *)calloc((*anz)+1,sizeof(double)); if (!*v) speichererror(1);
	for (i=0;i<=(*anz);i++)  (*v)[i] = -fehlwert -1;

	(*WSstream) << " ("<<(*anz)<<") Stationen "<<"[end-entry]";

	fclose(statfile);
	if ((statfile = fopen(statfilename,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening file "<<statfilename<<"[end-entry]";
		throw(49);
	}

	do
	{ /* ueberlesen des ersten Satzes, der die Datenart kennzeichnet */
		meteosatz=fgets(meteosatz, linelength,statfile);
		i=strlen(meteosatz); if (i==0) i=1;
	} 
	while (meteosatz[i-1]!='\n'); /* wenn der reservierte Speicher fuer den Satz zu kurz war, weiter lesen */

	/* Einlesen der Stationshoehen aus der 2. anetzzeile */
	fscanf(statfile,"%s%s%s%s",js,ms,ts,ss);
	if (*useMinutes) fscanf(statfile, "%s", mins);
	for (i=1;i<=numAddCol;i++) fscanf(statfile,"%s",js); /* reading the headers of add. columns away */
	for (i=1;i<=(*anz);i++)
	{
		fscanf(statfile,"%s",ss);
		(*h)[i] = atof(ss);
	}

	/* Einlesen der x-Koordinate aus der 3. anetzzeile */
	fscanf(statfile,"%s%s%s%s",js,ms,ts,ss);
	if (*useMinutes) fscanf(statfile, "%s", mins);
	for (i=1;i<=numAddCol;i++) fscanf(statfile,"%s",js); /* reading the headers of add. columns away */
	for (i=1;i<=(*anz);i++)
	{
		fscanf(statfile,"%s",ss);
		(*xk)[i] = atof(ss);
	}

	/* Einlesen der y-Koordinate aus der 4. anetzzeile */
	fscanf(statfile,"%s%s%s%s",js,ms,ts,ss);
	if (*useMinutes) fscanf(statfile, "%s", mins);
	for (i=1;i<=numAddCol;i++) fscanf(statfile,"%s",js); /* reading the headers of add. columns away */
	for (i=1;i<=(*anz);i++)
	{
		fscanf(statfile,"%s",ss);
		(*yk)[i] = atof(ss);
	}

	//  meteosatz=fgets(meteosatz,(20*(*anz))>_stringsize_?20*(*anz):_stringsize_,statfile);
	if (fclose(statfile)==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<statfilename<<"[end-entry]";
		throw(50);
	}
}

void meteointerpol::finddateanetz(FILE *statfile,char *statfilename,
					int numAddCol,bool useMinutes,int anz,double *v,long *statfilepos,long *statfilepos_preread,
					unsigned short *gw,
					long *me,long *se,long *tg,long *mt,long *jr)
{
	int i,j;
	char dummy[100];

	/* Suchen des Anfangsdatums in der statfile
datei */
	(*WSstream) << "looking for starting date "<< wasim::setw(2) << taga<<"."<< wasim::setw(2) << monata<<"."<<jahra<<" "<< wasim::setw(2) << stundea<<":"<< wasim::setw(2) << minutea<<" for "<<instanzname<<" in "<<statfilename<<"[end-entry]";

	if ((statfile = fopen(statfilename,"rb"))==NULL)
	{
		(*WSstream) << "FATAL: \nerror opening file "<<statfilename<<"[end-entry]";
		throw(49);
	}

	/* es wird ein meteosatz weniger vorgespult, als Kopfzeilen da sind, weil dann beim Suchen
	immer erst ein meteosatz weggelesen wird (der naemlich schon fuer das Datum angelesen wurde
	aber beim erstenmal noch gar nicht gelesen wurde */
	for (i=1;i<=4;i++)  
	{
		do
		{ /* ueberlesen der ersten 4 Saetze mit */
			meteosatz=fgets(meteosatz,(20*anz)>_stringsize_?20*anz:_stringsize_,statfile);
			j=strlen(meteosatz); if (j==0) j=1;
		} 
		while (meteosatz[j-1]!='\n'); /* wenn der reservierte Speicher fuer den Satz zu kurz war, weiter lesen */
	}
	// meteosatz=fgets(meteosatz,(20*anz)>_stringsize_?20*anz:_stringsize_,statfile);
	for (j=1;j<=anz;j++) v[j] = -fehlwert - 1.0;

	while ( (! feof(statfile)) && ! ( (*tg == taga) && (*mt == monata)
	&& (*jr == jahra) && (*se == stundea) && (*me == minutea)))
	{
		if (! feof(statfile))
		{
			meteosatz=fgets(meteosatz,(20*anz)>_stringsize_?20*anz:_stringsize_,statfile);
			*statfilepos_preread=ftell(statfile);
			fscanf(statfile,"%ld%ld%ld%ld",jr,mt,tg,se);
			if (useMinutes) fscanf(statfile, "%ld", me); else *me = 0;
			if (*jr < 100)  *jr = *jr+1900;
			/* reading the contents of add. columns into dummy variables */
			for (i=1;i<=numAddCol;i++) fscanf(statfile,"%s",dummy); 
		}
		if (((*se >= 24) || ((*se == 0) && (*me == 0))) && (*tg == 1))
		{
			(*WSstream) << ".";
		}
	}
	if (((*tg!=taga) || (*se != stundea) || (*me != minutea) ||(*mt!=monata)||(*jr!=jahra)) && (feof(statfile)))
	{
		(*WSstream) << "\nfile "<<statfilename<<" dosnt contain date "<< wasim::setw(2) << taga	<<"."<< wasim::setw(2) << monata<<"."<<jahra<<", " << wasim::setw(2) << stundea<<":"<< wasim::setw(2) << minutea<<"(feof)!"<<"[end-entry]";
		throw(51);
	}

	for (i=1;i<=anz;i++)
	{
		fscanf(statfile,"%lf",&(v[i]));
		if ((v[i] < -fehlwert)||(v[i] > fehlwert)) v[i] = -fehlwert -1;
	}
	meteosatz=fgets(meteosatz,(20*anz)>_stringsize_?20*anz:_stringsize_,statfile);

	*gw = false;
	(*WSstream) << "[end-entry]";
	*statfilepos=ftell(statfile);
	if (fclose(statfile)==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<statfilename<<"[end-entry]";
		throw(50);
	}
}


// for method 9 (read grids according to a grid list in Wasim-format)

void meteointerpol::finddategridlist(FILE *statfile,char *statfilename,
					int numAddCol,bool useMinutes,char *gridname,long *statfilepos,long *statfilepos_preread,
					long *me,long *se,long *tg,long *mt,long *jr)
{
	int i,j;
	char dummy[100];

	// Suchen des Anfangsdatums in der statfile datei 
	(*WSstream) << "looking for starting date "<<taga<<"."<<monata<<"."<<jahra<<" for "<<instanzname<<" in "<<statfilename<<"[end-entry]";

	if ((statfile = fopen(statfilename,"rb"))==NULL)
	{
		(*WSstream) << "FATAL: \nerror opening file "<<statfilename<<"[end-entry]";
		throw(49);
	}

	// es wird ein meteosatz weniger vorgespult, als Kopfzeilen da sind, weil dann beim Suchen
	// immer erst ein meteosatz weggelesen wird (der naemlich schon fuer das Datum angelesen wurde
	// aber beim erstenmal noch gar nicht gelesen wurde 

	for (i=1;i<=4;i++)  
	{
		do
		{ // ueberlesen der ersten 4 Saetze mit 
			meteosatz=fgets(meteosatz,_stringsize_,statfile);
			j=strlen(meteosatz); if (j==0) j=1;
		} 
		while (meteosatz[j-1]!='\n'); // wenn der reservierte Speicher fuer den Satz zu kurz war, weiter lesen 
	}

	// von ganz früher meteosatz=fgets(meteosatz,(20*anz)>_stringsize_?20*anz:_stringsize_,statfile);
	// 20.05.2007 for (j=1;j<=anz;j++) v[j] = -fehlwert - 1.0;

	while ( (! feof(statfile)) && ! ( (*tg == taga) && (*mt == monata)	&& (*jr == jahra) && (*se == stundea) && (*me == minutea)))
	{
		if (! feof(statfile))
		{
			meteosatz=fgets(meteosatz,_stringsize_,statfile);
			*statfilepos_preread=ftell(statfile);
			fscanf(statfile,"%ld%ld%ld%ld",jr,mt,tg,se);
			if (useMinutes) fscanf(statfile, "%ld", me); else *me = 0;
			if (*jr < 100)  *jr = *jr+1900;
			// reading the contents of add. columns into dummy variables 
			for (i=1;i<=numAddCol;i++) fscanf(statfile,"%s",dummy); 
		}
		if (((*se >= 24)||((*se==0)&&(*me==0)))&&(*tg==1))
		{
			(*WSstream) << ".";
		}
	}
	if (((*tg!=taga)||(*se!=stundea)||(*mt!=monata)||(*jr!=jahra)) && (feof(statfile)))
	{
		(*WSstream) << "FATAL: \nfile "<<statfilename<<" dosnt contain date "<<taga
		<<"."<<monata<<"."<<jahra<<", hour "<<stundea<<"(feof)!"<<"[end-entry]";
		throw(51);
	}

	// 20.05.2007  for (i=1;i<=anz;i++)
	// 20.05.2007  {
	// 20.05.2007    fscanf(statfile,"%lf",&(v[i]));
	// 20.05.2007    if ((v[i] < -fehlwert)||(v[i] > fehlwert)) v[i] = -fehlwert -1;
	// 20.05.2007  }

	// neu 20.05.2007
	fscanf(statfile,"%s",gridname);
	meteosatz=fgets(meteosatz,_stringsize_,statfile);

	(*WSstream) << "[end-entry]";
	*statfilepos=ftell(statfile);
	if (fclose(statfile)==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<statfilename<<"[end-entry]";
		throw(50);
	}
}


void meteointerpol::finddateregress(int nr)
{
	char c[5];
	/* Suchen des Anfangsdatums in der Koeffizientendatei bei Regression */
	(*WSstream) << "looking for starting date "<<taga<<"."<<monata<<"."<<jahra<<" for "
	<<instanzname<<" in "<<&(dateinamekoef[nr][0])<<"[end-entry]";
	/* check, if the file is binary or ascii */
	ascii[nr]=0; /* initializing the binary/ascii-flag with binary (0) */
	regresssatz[nr]=newstring(_stringsize_+1,1);
	if ((dateikoef[nr] = fopen(dateinamekoef[nr],"rb")) == NULL)
	{
		(*WSstream) << "FATAL: \nerror opening file " << &(dateinamekoef[nr][0])<<"[end-entry]";
		throw(49);
	}
	fread(c,sizeof(char),5,dateikoef[nr]);
	if ((c[0]=='j')&&(c[1]=='r')&&(c[3]=='m')&&(c[4]=='o'))
	{
		ascii[nr]=1;
		regresssatz[nr]=fgets(regresssatz[nr],_stringsize_,dateikoef[nr]);
	}
	else
	{
		if (regresssatz[nr]) free(regresssatz[nr]);
		regresssatz[nr]=NULL;
		fseek(dateikoef[nr],0,0);
	}
	monat = tag = stunde = jahr = 0;
	while (((tag!=taga)||(stunde!=stundea)||(monat!=monata)||(jahr!=jahra)) && (!feof(dateikoef[nr])))
	{
		if (ascii[nr])
		{
			regfilepos_preread[nr]=ftell(dateikoef[nr]);
			regresssatz[nr]=fgets(regresssatz[nr],_stringsize_,dateikoef[nr]);
			if (!regresssatz[nr]) break;
			sscanf(regresssatz[nr],"%ld%ld%ld%ld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&(jahr),&(monat),&(tag),&(stunde),
			&(a[nr][1]),&(b[nr][1]),&(igu[nr]),&(a[nr][2]),&(b[nr][2]),&(igo[nr]),&(a[nr][3]),&(b[nr][3]),&(af[nr]),&(bf1[nr]),&(bf2[nr]));
		}
		else
		{
			regfilepos_preread[nr]=ftell(dateikoef[nr]);
			fread(&(jahr)  ,sizeof(jahr),1,dateikoef[nr]);
			fread(&(monat) ,sizeof(monat),1,dateikoef[nr]);
			fread(&(tag)   ,sizeof(tag),1,dateikoef[nr]);
			fread(&(stunde),sizeof(stunde),1,dateikoef[nr]);

			fread(&(a[nr][1])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(b[nr][1])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(igu[nr])   ,sizeof(double),1,dateikoef[nr]);
			fread(&(a[nr][2])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(b[nr][2])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(igo[nr])   ,sizeof(double),1,dateikoef[nr]);
			fread(&(a[nr][3])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(b[nr][3])  ,sizeof(double),1,dateikoef[nr]);
			fread(&(af[nr])    ,sizeof(double),1,dateikoef[nr]);
			fread(&(bf1[nr])   ,sizeof(double),1,dateikoef[nr]);
			fread(&(bf2[nr])   ,sizeof(double),1,dateikoef[nr]);
		}
		if (jahr<1000)  jahr=jahr+1900;
		if (((stunde >= 24)||(stunde==0))&&(tag==1))
		{
			(*WSstream) << ".";
		}
	}
	(*WSstream) << "[end-entry]";
	if (((tag!=taga)||(stunde!=stundea)||(monat!=monata)||(jahr!=jahra)) && (feof(dateikoef[nr])))
	{
		(*WSstream) << "FATAL: \nEOF reached: "<<dateinamekoef[nr]<<", date "<<taga
		<<"."<<monata<<"."<<jahra<<", "<<stundea<<". h not found!"<<"[end-entry]";
		throw(51);
	}
	gleichewerte[nr] = false;
	regfilepos[nr]=ftell(dateikoef[nr]);
	if (fclose(dateikoef[nr])==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<dateinamekoef[nr]<<"[end-entry]";
		throw(50);
	}
}

void meteointerpol::nextrecord()
{
	int nr;
	long int stundet,tagt,monatt,jahrt;
	/* die idhw-Methode bildet aus idw und regress ein gewichtetes Mittel
	Deshalb werden beide Eingabedateien benoetigt */
	if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12)) 
	{
		for (nr=0;nr<statfilecnt;nr++) 
		{
			nextrecordanetz(nr);
			if (nr==0)
			{
				minute = minute_a[0];
				stunde=stunde_a[0];
				tag=tag_a[0];
				monat=monat_a[0];
				jahr=jahr_a[0];
			}
		}
		for (nr=0;nr<statfilecnt;nr++) 
		{
			if ((minute != minute_a[nr]) || 
				(stunde!=stunde_a[nr])||
				(tag!=tag_a[nr])||
				(monat!=monat_a[nr])||
				(jahr!=jahr_a[nr]))
			{
				(*WSstream) << "FATAL: error reading date from stationfiles:"<<"[end-entry]";
				for (nr=0;nr<statfilecnt;nr++)
				(*WSstream) << "file "<<anetzstat[nr]<<": YYYY="<<jahr_a[nr]
				<<", MM="<<monat_a[nr]<<", DD="<<tag_a[nr]
				<<", HH="<<stunde_a[nr]<<" mm="<<minute_a[nr]<<"[end-entry]";
				(*WSstream) << "[end-entry]";
				throw(73);
			}
		}
	}
	if ((methode==2)||(methode==3))
	{
		stundet=stunde;
		tagt=tag;
		monatt=monat;
		jahrt=jahr;
		for (nr=1;nr<=regiomax;nr++) 
		{
			nextrecordregress(nr);
			if (((methode==3)||((methode==2)&&(nr>1)))&&((jahr!=jahrt)||(monat!=monatt)||(tag!=tagt)||(stunde!=stundet)))
			{
				(*WSstream) << "FATAL: -----------------------------------------------------------------------------\n"
				"error in date column of "<<instanzname<<" for IDW*Regress \n\n"
				"Date from station/regress file "
				<<(methode==3?(anetzstat[0]):dateinamekoef[nr])<<":\n"
				<<tagt<<"."<<monatt<<"."<<jahrt<<" "<<stundet<<" h\n\n"
				"but date from regres file "<<dateinamekoef[nr]<<":\n"
				<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde<<" h\n\n"
				"-----------------------------------------------------------------------------\n" << "[end-entry]";
				throw(52);
			}
			stundet=stunde;
			tagt=tag;
			monatt=monat;
			jahrt=jahr;
		}
	}
	if (methode==9)
	{
		nextrecordgridlist();
		minute = minute_a[0];
		stunde=stunde_a[0];
		tag=tag_a[0];
		monat=monat_a[0];
		jahr=jahr_a[0];
	}
}

void meteointerpol::readbufferanetz(int nr)
{
	char *tempsatz;
	int i;

	if (((anetz[nr]) = fopen((anetzstat[nr]),"rb"))==NULL)
	{
		(*WSstream) << "FATAL: \nerror opening file "<<(anetzstat[nr])<<"[end-entry]";
		throw(49);
	}
	tempsatz=meteosatz;
	fseek(anetz[nr],anetzfilepos[nr],0);
	if (!feof((anetz[nr])))
	{
		bufcntanetz[nr]=0;
		while ((!feof((anetz[nr])))&&(bufcntanetz[nr]<bufmaxanetz[nr]))
		{
			bufcntanetz[nr]++;
			fscanf((anetz[nr]),"%ld%ld%ld%ld",
			&(bufferanetz[nr][(int)bufcntanetz[nr]].j),&(bufferanetz[nr][(int)bufcntanetz[nr]].m),&(bufferanetz[nr][(int)bufcntanetz[nr]].t),&(bufferanetz[nr][(int)bufcntanetz[nr]].s));
			if (useMinutes[nr]) fscanf((anetz[nr]), "%ld", &(bufferanetz[nr][(int)bufcntanetz[nr]].min));
			/* reading the contents of add. columns into dummy variables */
			for (i=1;i<=numAdditionalColumns[nr];i++) fscanf((anetz[nr]),"%s",js); 
			for (i=1;i<=(anzahl[nr]);i++)
			{
				fscanf((anetz[nr]),"%lf",&(bufferanetz[nr][(int)bufcntanetz[nr]].wert[i]));
				if ((bufferanetz[nr][(int)bufcntanetz[nr]].wert[i] < -fehlwert)||(bufferanetz[nr][(int)bufcntanetz[nr]].wert[i] > fehlwert))
				bufferanetz[nr][(int)bufcntanetz[nr]].wert[i] = -fehlwert - 1;
			}
			meteosatz=fgets(meteosatz,(20*(anzahl[nr]))>_stringsize_?20*(anzahl[nr]):_stringsize_,(anetz[nr]));
		}
		if (feof((anetz[nr]))) 
		{
			for (i=bufcntanetz[nr]+1;i<=bufmaxanetz[nr];i++)
			{  
				free(bufferanetz[nr][i].wert);
				bufferanetz[nr][i].wert=NULL;
			}
			bufmaxanetz[nr]=bufcntanetz[nr];  /* if end of input file is reached, the max-indicator is set to the number of valid records in the buffer */
		}
		bufcntanetz[nr] = bufmaxanetz[nr]>=1? 1 : 0;  /* pointing to the actual valid record */
	}
	else
	ende=true;
	meteosatz=tempsatz;
	anetzfilepos[nr]=ftell((anetz[nr]));
	if (fclose((anetz[nr]))==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<(anetzstat[nr])<<"[end-entry]";
		throw(50);
	}
}

void meteointerpol::nextrecordgridlist()
{
	int i;
	double min,max;
	if ((jf*jahr_a[0]+mf*monat_a[0]+df*tag_a[0]+hf*stunde_a[0]+minf*minute_a[0])>=endtime)
	{
		ende=true;
		if (anetz[0])
		{
			fclose(anetz[0]);
			anetz[0] = NULL;
		}
	}	
	else
	{
		if (! ende)
		{
			if (!anetz[0])
			{
				if (((anetz[0]) = fopen((anetzstat[0]),"rb"))==NULL)
				{
					(*WSstream) << "FATAL: \nerror opening file "<<(anetzstat[0])<<"[end-entry]";
					throw(49);
				}
				fseek(anetz[0],anetzfilepos[0],0); // positioning the file pointer
			}
			if (! feof(anetz[0]))
			{
				fscanf(anetz[0], "%ld%ld%ld%ld", &(jahr_a[0]), &(monat_a[0]), &(tag_a[0]), &(stunde_a[0]));
				if (useMinutes[0]) 
//				if ((timestep < 60) && (timestep > 0))
				{
					fscanf(anetz[0], "%ld", &(minute_a[0]));
				}
				else
				{
					minute_a[0] = 0;
				}
				if (jahr_a[0] < 100)  jahr_a[0] = jahr_a[0]+1900;
				// reading the contents of add. columns into dummy variables 
				for (i=1;i<=numAdditionalColumns[0];i++) fscanf((anetz[0]),"%s",js); 
				fscanf(anetz[0],"%s",this->gridlist_actname);
				meteosatz=fgets(meteosatz,_stringsize_,anetz[0]);
			}
		}
	}
}

void meteointerpol::nextrecordanetz(int nr)
{
	int i;
	double min,max;
	if ((jf*jahr_a[nr] + mf*monat_a[nr] + df*tag_a[nr] + hf*stunde_a[nr] + minf*minute_a[nr]) >= endtime)
	{
		ende = true;
	}
	else
	{
		if (bufcntanetz[nr] > bufmaxanetz[nr])
		{
			readbufferanetz(nr);
		}

		if (! ende)
		{
			jahr_a[nr]   = bufferanetz[nr][(int)bufcntanetz[nr]].j;
			monat_a[nr]  = bufferanetz[nr][(int)bufcntanetz[nr]].m;
			tag_a[nr]    = bufferanetz[nr][(int)bufcntanetz[nr]].t;
			stunde_a[nr] = bufferanetz[nr][(int)bufcntanetz[nr]].s;
			minute_a[nr] = bufferanetz[nr][(int)bufcntanetz[nr]].min;
			if ((jahr_a[nr]==0)&&(monat_a[nr]==0)&&(tag_a[nr]==0))
			{
				ende=true;
			}
			else
			{
				if (jahr_a[nr] < 100)  jahr_a[nr] = jahr_a[nr]+1900;
				for (i=1;i<=(anzahl[nr]);i++)
				value[nr][i]=bufferanetz[nr][(int)bufcntanetz[nr]].wert[i];
				gleichewerte[nr] = true;
				min=max=-fehlwert-1;
				for (i=1;i<=(anzahl[nr]);i++) if (value[nr][i]>= -fehlwert) {min=max=value[nr][i]; break;}
				for (i=1;i<=(anzahl[nr]);i++)
				{
					if ((value[nr][i]<min)&&(value[nr][i] >= -fehlwert)) min=value[nr][i];
					else
					if ((value[nr][i]>max)&&(value[nr][i] >= -fehlwert)) max=value[nr][i];
				}
				if (min!=max) gleichewerte[nr]=false;
			}
		}
		if (!ende) bufcntanetz[nr]++;
	}
}

void meteointerpol::readbufferregress(int nr)
{
	int i;
	if (ascii[nr])
	{
		if ((dateikoef[nr] = fopen(dateinamekoef[nr],"rb"))==NULL)
		{
			(*WSstream) << "FATAL: \nerror opening file "<<dateinamekoef[nr]<<"[end-entry]";
			throw(49);
		}
	}
	else
	{ /* opening in binary mode */
		if ((dateikoef[nr]=fopen(dateinamekoef[nr],"rb"))==NULL)
		{
			(*WSstream) << "FATAL: error opening regression file "<<dateinamekoef[nr]<<" in binary mode"<<"[end-entry]";
			throw(49);
		}
	}
	fseek(dateikoef[nr],regfilepos[nr],0);

	if (!feof(dateikoef[nr]))
	{
		bufcntregress[nr]=0;
		while ((!feof(dateikoef[nr]))&&(bufcntregress[nr]<bufmaxregress[nr]))
		{
			bufcntregress[nr]++;
			if (ascii[nr])
			{
				fscanf(dateikoef[nr],"%ld%ld%ld%ld",
				&(bufferregress[nr][(int)(bufcntregress[nr])].j),
				&(bufferregress[nr][(int)(bufcntregress[nr])].m),
				&(bufferregress[nr][(int)(bufcntregress[nr])].t),
				&(bufferregress[nr][(int)(bufcntregress[nr])].s));
				// setting the minute value to 0 until a minute column is supported by regr.exe
				bufferregress[nr][(int)(bufcntregress[nr])].min = 0;
				for (i=1;i<=11;i++)
				{
					fscanf(dateikoef[nr],"%lf",&(bufferregress[nr][(int)(bufcntregress[nr])].wert[i]));
				}
				regresssatz[nr]=fgets(regresssatz[nr],_stringsize_,dateikoef[nr]);
			}
			else
			{
				/* reading 4 long values (date) */
				fread(&(bufferregress[nr][(int)(bufcntregress[nr])].j),sizeof(long),4,dateikoef[nr]);
				// setting the minute value to 0 until a minute column is supported by regr.exe
				bufferregress[nr][(int)(bufcntregress[nr])].min = 0;
				/* reading the 11 regression-parameters */
				fread(&(bufferregress[nr][(int)(bufcntregress[nr])].wert[1]),sizeof(double),11,dateikoef[nr]);
			}
		}
		if (feof(dateikoef[nr])) 
		{
			for (i=bufcntregress[nr]+1;i<=bufmaxregress[nr];i++)
			{ 
				free(bufferregress[nr][i].wert);
				bufferregress[nr][i].wert=NULL;
			}
			bufmaxregress[nr]=bufcntregress[nr];  /* if end of input fileis reached, the max-indicator is set to the number of valid records in the buffer */
		}
		bufcntregress[nr]=bufmaxregress[nr]>=1?1:0;  /* pointing to the actual valid record */
	}
	else
	ende=true;

	regfilepos[nr]=ftell(dateikoef[nr]);
	if (fclose(dateikoef[nr])==EOF)
	{
		(*WSstream) << "FATAL: error closing "<<dateinamekoef[nr]<<"[end-entry]";
		throw(50);
	}
}

void meteointerpol::nextrecordregress(int nr)
{
	/* check for the end date and time: jf, mf, df and hf are factors for multiplying the year, month, 
	day and hour (defined in modelobject) in order to be able to compare different dates with 
	"<" and ">" operators */
	if (((jf*jahr+mf*monat+df*tag+hf*stunde+minf*minutee)>=endtime)&&(nr==1))
	{
		if (methode!=3)
		/* wenn methode == 3 ist, dann wird vorher schon nextrecordanetz durchgefuehrt,
	wobei die stunde schon weitergesetzt wurde, was hier zum vorzeitigen Erkennen
	des Endes fuehren wurde. Deshalb braucht und darf diese Abfrage im modus 3 nicht
	durchgefuehrt (zu) werden */
		ende = true;
	}
	if (!ende)
	{
		if (bufcntregress[nr]>bufmaxregress[nr])
		readbufferregress(nr);

		if (! ende)
		{
			jahr   = bufferregress[nr][(int)(bufcntregress[nr])].j;
			monat  = bufferregress[nr][(int)(bufcntregress[nr])].m;
			tag    = bufferregress[nr][(int)(bufcntregress[nr])].t;
			stunde = bufferregress[nr][(int)(bufcntregress[nr])].s;
			minute = 0;
			if ((jahr==0)&&(monat==0)&&(tag==0))
			{
				ende=true;
			}
			else
			{
				a[nr][1]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[1];
				b[nr][1]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[2];
				igu[nr]    = bufferregress[nr][(int)(bufcntregress[nr])].wert[3];
				a[nr][2]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[4];
				b[nr][2]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[5];
				igo[nr]    = bufferregress[nr][(int)(bufcntregress[nr])].wert[6];
				a[nr][3]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[7];
				b[nr][3]   = bufferregress[nr][(int)(bufcntregress[nr])].wert[8];
				af[nr]     = bufferregress[nr][(int)(bufcntregress[nr])].wert[9];
				bf1[nr]    = bufferregress[nr][(int)(bufcntregress[nr])].wert[10];
				bf2[nr]    = bufferregress[nr][(int)(bufcntregress[nr])].wert[11];
				if (jahr<1000)  jahr=jahr+1900;
			}
		}
		if (!ende) bufcntregress[nr]++;
	}
}

int meteointerpol::isend()
{
	return (ende>0);
}

void meteointerpol::set_end_date()
{
	/* if Ctrl-C was pressed, this interval should be the last one, so the end date is set to actual date */
	minutee = minute;
	stundee = stunde;
	tage    = tag;
	monate  = monat;
	jahre   = jahr;
	endtime = jf*jahre+mf*monate+df*tage+hf*stundee+minf*minutee;
}




double meteointerpol::getregress(int row,int column,double x,double y,double z)
{
	int nr,region,nearestregion_index;
	double dx,dy,dist,mindist;
	double wert[4];

	/* --------------- without region transition -> if no region grid was read in 
or if the actual location doesnt belong to any region -------------- */
	if (WEIGHT)
	{
		region=nearestregion_index=rarray[(int)REGION[1][row][column]]; 
		if (REGION[1][row][column]<=0)
		{
			mindist=1e20;
			for (nr=1;nr<=regiomax;nr++)
			{
				dx=(x!=nodata)?x-center_of_mass[nr].xcoord:center_of_mass[nr].xcoord;
				dy=(y!=nodata)?y-center_of_mass[nr].ycoord:center_of_mass[nr].ycoord;
				dist=dx*dx+dy*dy;
				if (dist<mindist)
				{
					mindist=dist;
					nearestregion_index=nr;
				}
			}
			row=center_of_mass[nearestregion_index].row;
			column=center_of_mass[nearestregion_index].col;
		}
	}
	else
	{
		region=1;
		nearestregion_index=0;
	}

	if ((!WEIGHT)||((WEIGHT)&&(nearestregion_index<=0)))
	{
		if (a[region][1] < -fehlwert) return -fehlwert - 1.0;
		/* Das Regress-Programm schreibt einen Fehlwert fuer den Anstieg der untersten Geradein die
	Regress-Datei, wenn alle Daten Fehlwerte waren, deshalb kann hier dann sofort ein Fehlwert
	zurueckgeliefert werden, wenn a[1] ein Fehlwert ist */
		if ((af[region] >= -fehlwert)&&(x!=nodata)&&(y!=nodata))
		{
			if (z > igo[region])
			/* oberste Gerade */
			wert[0] = (double)( af[region] + bf1[region]*x + bf2[region]*y - a[region][1] + a[region][3] + b[region][3]*z);
			else
			{
				if (z < igu[region])
				/* unterste Gerade */
				{
					wert[0] = (double)( af[region] + bf1[region]*x + bf2[region]*y /* - a[1] + a[1]*/ + b[region][1]*z);
				}
				else
				{
					/* mittlere Gerade */
					wert[0] = (double)( af[region] + bf1[region]*x + bf2[region]*y - a[region][1] + a[region][2] + b[region][2]*z);
				}
			}
		}
		else /* wenn in der Regression alle Werte gleich waren (z.B. alle == 0),
			dann stehen in den Parametern fuer die areale Regression nur
			-999 - Werte drin -> dann sollen keine Interpolationen aus
			arealen Regressionen berechnet werden  */
		{
			if (z > igo[region])
			/* oberste Gerade */
			wert[0] = (double)( a[region][3] + b[region][3]*z);
			else
			{
				if (z < igu[region])
				/* unterste Gerade */
				{
					wert[0] = (double)( a[region][1] + b[region][1]*z);
				}
				else
				{
					/* mittlere Gerade */
					wert[0] = (double)( a[region][2] + b[region][2]*z);
				}
			}
		}
	}
	/* ---------------------------- with different regions -> transition of results -------------- */
	else
	{
		wert[0]=0;
		for (nr=1;nr<=regiomax&&nr<=3;nr++)  /* regiomax may be 2 or 3, thats why the loop condition reads nr<=regiomax&&nr<=3 */
		{
			/* getting the region code by using the row/colum data which where get as parameters or, 
		if there was no row/column given, are calculated at first in this function (e.g., if getsinglepoint 
		is called, which does not know its coordinates ) */
			region=rarray[(int)REGION[nr][row][column]]; 
			if (WEIGHT[nr][row][column]>0.0)
			{
				if (a[region][1] < -fehlwert) return -fehlwert - 1.0;
				if ((af[region] >= -fehlwert)&&(x!=nodata)&&(y!=nodata))
				{ /* with areal regression of residuals */
					if (z > igo[region])
					wert[nr] = (double)( af[region] + bf1[region]*x + bf2[region]*y - a[region][1] + a[region][3] + b[region][3]*z);
					else
					{
						if (z < igu[region])
						{
							wert[nr] = (double)( af[region] + bf1[region]*x + bf2[region]*y /* - a[1] + a[1]*/ + b[region][1]*z);
						}
						else
						{
							wert[nr] = (double)( af[region] + bf1[region]*x + bf2[region]*y - a[region][1] + a[region][2] + b[region][2]*z);
						}
					}
				}
				else
				{ /* without areal regression */
					if (z > igo[region])
					wert[nr] = (double)( a[region][3] + b[region][3]*z);
					else
					{
						if (z < igu[region])
						{
							wert[nr] = (double)( a[region][1] + b[region][1]*z);
						}
						else
						{
							wert[nr] = (double)( a[region][2] + b[region][2]*z);
						}
					}
				}
				wert[0]+=WEIGHT[nr][row][column]*wert[nr];
			}
		}
	}
	return wert[0];
}

double meteointerpol::determinante_drei(double *a,double *b,double *c)
{
	return ( a[1]*b[2]*c[3]+b[1]*c[2]*a[3]+c[1]*a[2]*b[3]-
	a[3]*b[2]*c[1]-b[3]*c[2]*a[1]-c[3]*a[2]*b[1]);
}

void meteointerpol::CalculateRegressionParameters(int nr)
{
	// these values have to be set by this funktion: double *af,*bf1,*bf2,*igu,*igo,**a,**b;
	// and this are the parameters used to calculate the resulting regression parameters
	//  double lowerinversionfix,upperinversionfix,tolerance,clusterlimit;
	//  int overlap;

	double     sx[3],sy[3],sx2[3],sy2[3],sxy[3],n[3];
	double     syi,sx1i,sx2i,sx1i2,sx1ix2i,sx2i2,sx1iyi,sx2iyi;
	double     minh[3],maxh[3],hdiff[3];

	/* statistische berechnung und Einzeinen der Punkte und regressionsgeraden */
	for (int ber=1;ber<=2;ber++)
	{
		sx[ber]   = 0.0 ;
		sy[ber]   = 0.0 ;
		sx2[ber]  = 0.0 ;
		sy2[ber]  = 0.0 ;
		sxy[ber]  = 0.0 ;
		n[ber]    = 0.0 ;
		minh[ber] = 0.0 ;
		maxh[ber] = 0.0 ;
	}

	for (int i=1;i<=anzahl[nr];i++)
	{
		/* Fehlwerte werden nicht mit gerechnet */
		if (value[nr][i] > -fehlwert) 
		{
			int ber=1;
			if (hoehe[nr][i] > lowerinversionfix)  ber = 2; else ber = 1;
			sx[ber]  += hoehe[nr][i];
			sy[ber]  += value[nr][i];
			n[ber]   += 1;
			sx2[ber] += hoehe[nr][i]*hoehe[nr][i];
			sy2[ber] += value[nr][i]*value[nr][i];
			sxy[ber] += hoehe[nr][i]*value[nr][i];
			if (n[ber]==1) /* der erste Wert initialisiert die Register zum feststellen der Hoehenausdehnung */
			{
				minh[ber] = hoehe[nr][i];
				maxh[ber] = hoehe[nr][i];
			}
			else
			{
				if (hoehe[nr][i] < minh[ber])  
				{
					minh[ber] = hoehe[nr][i]; 
				}
				else 
				{
					if (hoehe[nr][i]>maxh[ber]) 
					maxh[ber] = hoehe[nr][i];
				}
			}
		}
	}

	double xminoben=1000000000.0;
	double xmaxunten=-100000000.0;

	for (int i=1;i<=anzahl[nr];i++)
	if ((hoehe[nr][i] < xminoben) && (hoehe[nr][i] > lowerinversionfix))  xminoben= hoehe[nr][i];

	for (int i=1;i<=anzahl[nr];i++)
	if ((hoehe[nr][i] > xmaxunten) && (hoehe[nr][i] < lowerinversionfix))  xmaxunten=hoehe[nr][i];

	// looking for overlapping stations
	for (int i=1;i<=anzahl[nr];i++)
	{
		/* Fehlwerte werden nicht mit gerechnet */
		if  (value[nr][i] > -fehlwert) 
		{
			if ((hoehe[nr][i] == xmaxunten)&&((overlap>0)||(n[2]==1)) )
			{
				int ber=2;
				sx[ber]  += hoehe[nr][i];
				sy[ber]  += value[nr][i];
				n[ber]   += 1;
				sx2[ber] += hoehe[nr][i]*hoehe[nr][i];
				sy2[ber] += value[nr][i]*value[nr][i];
				sxy[ber] += hoehe[nr][i]*value[nr][i];
				if (hoehe[nr][i] < minh[ber])  
				{
					minh[ber] = hoehe[nr][i]; 
				}
				else 
				{
					if (hoehe[nr][i]>maxh[ber]) 
					maxh[ber] = hoehe[nr][i];
				}
			}

			if ((hoehe[nr][i] == xminoben)&&((overlap>0)||(n[1]==1)) )
			{
				int ber=1;
				sx[ber]  += hoehe[nr][i];
				sy[ber]  += value[nr][i];
				n[ber]   += 1;
				sx2[ber] += hoehe[nr][i]*hoehe[nr][i];
				sy2[ber] += value[nr][i]*value[nr][i];
				sxy[ber] += hoehe[nr][i]*value[nr][i];
				if (hoehe[nr][i] < minh[ber])  
				{
					minh[ber] = hoehe[nr][i]; 
				}
				else 
				{
					if (hoehe[nr][i]>maxh[ber]) 
					maxh[ber] = hoehe[nr][i];
				}
			}
		}
	}

	for (int i=1;i<=2;i++)
	{
		/* wenn die Stationen zu eng aneinander liegen, werden Cluster gebildet */
		hdiff[i] = maxh[i]-minh[i];
		if ((hdiff[i] < clusterlimit)&&(n[i]>1)) 
		{
			sx[i]   = sx[i] / n[i];
			sy[i]   = sy[i] / n[i];
			n[i]    = 1;
			sx2[i]  = sx[i] * sx[i];
			sy2[i]  = sy[i] * sy[i];
			sxy[i]  = sx[i] * sy[i];
			minh[i] = sx[i];
			maxh[i] = sx[i];
		}
	}

	if ((n[1]==1)&&(n[2]==1)) 
	{
		if ((maxh[2]-minh[1])>=clusterlimit)
		{
			sx[1]  = sx[1]+sx[2];    sx[2]  = sx[1];
			sy[1]  = sy[1]+sy[2];    sy[2]  = sy[1];
			n[1]   = 2;              n[2]   = n[1];
			sx2[1] = sx2[1]+sx2[2];  sx2[2] = sx2[1];
			sy2[1] = sy2[1]+sy2[2];  sy2[2] = sy2[1];
			sxy[1] = sxy[1]+sxy[2];  sxy[2] = sxy[1];
		} 
		else
		{
			sx[2]=sy[2]=sx2[2]=sy2[2]=sxy[2]=n[2]=0;
		}
	}

	// now the partial sums for calculating the regression parameters are known...
	// estimating parameters now

	/* 0 = nur Anetz, 1=Anetz und Klimastationen */
	if ((n[1]==0)&&(n[2]==0))
	{ /* zu wenige Werte */
		a[nr+1][1] = -fehlwert-1;
		b[nr+1][1] = -fehlwert-1;
		a[nr+1][2] = -fehlwert-1;
		b[nr+1][2] = -fehlwert-1;
		igu[nr+1]  = -fehlwert-1;
		igo[nr+1]  = -fehlwert-1;
	}
	else
	{
		for (int m=1;m<=2;m++)
		{
			double quotient=0;
			if (n[m]>1)  quotient=(sx2[m]-((sx[m])*(sx[m]))/n[m]);
			if (quotient!=0)
			{ /* quotient wird null, wenn entweder
		nur ein Wert da ist, durch dan man natuerlich keine Gerade legen kann,
		oder wenn beide Werte exakt die Gleiche Abszisse haben, was zu einem
		unendlich grossen Anstieg fuehren wuerde */
				b[nr+1][m] = (sxy[m]-(sx[m]*sy[m])/n[m])/quotient;
				a[nr+1][m] = (sy[m]-b[nr+1][m]*sx[m])/n[m];
			}
			else
			{
				b[nr+1][m] = -fehlwert-1;
				a[nr+1][m] = -fehlwert-1;
			}
		}
		/* Wenn ein Bereich nur einen Wert hat, der andere aber mehr, wird eine Verbindungsgerade */
		if (n[1]==0) 
		{ /* unten kein Wert */
			if (n[2]==1) 
			{
				/* nur ein Punkt ueberhaupt, Senkrechte definieren */
				a[nr+1][2]=sy[2]; 
				b[nr+1][2]=0; 
				a[nr+1][1]=sy[2]; 
				b[nr+1][1]=0; 
			} 
			else
			{
				if (n[2]>1) 
				{
					/* oben eine echte Gerade, Parameter unten uebernehmen */
					a[nr+1][1]=a[nr+1][2];
					b[nr+1][1]=b[nr+1][2];
				}
			}
		}
		else 
		{
			if (n[1]==1) 
			{
				if (n[2]==0) 
				{
					/* nur ein Punkt unten, oben nichts, eine Senkrechte definieren*/
					a[nr+1][2]=sy[1]; 
					b[nr+1][2]=0; 
					a[nr+1][1]=sy[1]; 
					b[nr+1][1]=0; 
				}
				else
				{
					if (n[2]>1) 
					{
						/* obere Gerade bis igufix legen, dann verbinden zu unterem Punkt */
						b[nr+1][1]=(a[nr+1][2]+lowerinversionfix*b[nr+1][2]-sy[1])/(lowerinversionfix-sx[1]);
						a[nr+1][1]=a[nr+1][2]+lowerinversionfix*(b[nr+1][2]-b[nr+1][1]);
					}
				}
			}
			else
			{
				if (n[1]>1) 
				{
					if (n[2]==0) 
					{
						/* unten eine echte Gerade, Parameter oben uebernehmen */
						a[nr+1][2]=a[nr+1][1];
						b[nr+1][2]=b[nr+1][1];
					}
					else
					{
						if (n[2]==1) 
						{
							/* obere Gerade bis igufix legen, dann verbinden zu unterem Punkt */
							b[nr+1][2]=(a[nr+1][1]+lowerinversionfix*b[nr+1][1]-sy[2])/(lowerinversionfix-sx[2]);
							a[nr+1][2]=a[nr+1][1]+lowerinversionfix*(b[nr+1][1]-b[nr+1][2]);
						}
					}
				}
			}
		}

		/* Errechnen der inversionsgrenze */
		if ((a[nr+1][2]>-fehlwert) && ((b[nr+1][2]-b[nr+1][1])!=0))  igo[nr+1] = (a[nr+1][2]-a[nr+1][1]) / (b[nr+1][1]-b[nr+1][2]);
		else igo[nr+1]=lowerinversionfix;
		igu[nr+1] = 0;
		if ((igo[nr+1] < (lowerinversionfix - tolerance) ) || ( igo[nr+1] > (upperinversionfix + tolerance)))
		{
			igo[nr+1] = upperinversionfix ;
			igu[nr+1] = lowerinversionfix ;
		}
	}

	// moving resulting parameters to the correct location, depending on the existence of two or three valid regressions
	if (igu[nr+1]>0)
	{ /* es gibt 3 Geraden */
		a[nr+1][3] = a[nr+1][2]; // moving regression parameters 2 to regression 3
		b[nr+1][3] = b[nr+1][2];
		b[nr+1][2] = (( a[nr+1][3]+igo[nr+1]*b[nr+1][3] - a[nr+1][1]-igu[nr+1]*b[nr+1][1] ) / (igo[nr+1] - igu[nr+1]));
		a[nr+1][2] = (a[nr+1][1]  + igu[nr+1]*(b[nr+1][1]-b[nr+1][2]));
	}
	else
	{ /* es gibt nur 2 Geraden, die erste ist 0,0 */
		a[nr+1][3] = a[nr+1][2];
		b[nr+1][3] = b[nr+1][2];
		a[nr+1][2] = a[nr+1][1];
		b[nr+1][2] = b[nr+1][1];
		b[nr+1][1] = 0;
	}
	if (lapserate>-9998)
	{
		// set the uppermost regression tangens to the read value
		a[nr+1][3] = a[nr+1][3]+b[nr+1][3]*igo[nr+1] - lapserate*igo[nr+1];
		b[nr+1][3] = lapserate;
	}

	// now areal regression parameters
	double *residual;
	residual = new double[anzahl[nr]+1];
	sx2[nr] = 0;
	n[nr]   = 0;
	double delta=0;
	for (int i=1;i<=anzahl[nr];i++)
	{
		if ( value[nr][i]> -fehlwert)
		{
			if (hoehe[nr][i] < igu[nr+1])  
			{
				delta = value[nr][i]-(a[nr+1][1]+b[nr+1][1]*(hoehe[nr][i]));
			}
			else
			{
				if (hoehe[nr][i] > igo[nr+1])  
				{
					delta=value[nr][i]-(a[nr+1][3]+b[nr+1][3]*(hoehe[nr][i]));
				}
				else 
				{
					delta = value[nr][i]-(a[nr+1][2]+b[nr+1][2]*(hoehe[nr][i]));
				}
			}
			residual[i]=a[nr+1][1]+delta;

			/* zum Ermitteln der Reststreuung in bezug auf die Hoehenregression */

			n[nr]   += 1;
			sx2[nr] += delta*delta;

		}
	}

	if ((n[nr]>3)) // only for real regression without fix lapse rate we will calculate the areal regression
//	if ((n[nr]>3)&&(lapserate<-9998)) // only for real regression without fix lapse rate we will calculate the areal regression
	{ /* Ermitteln der Ausgleichsflaeche */
		n[nr]=0;  sx1i=0; sx2i=0; sx1i2=0; sx1ix2i=0; sx2i2=0; syi=0; sx1iyi=0; sx2iyi=0;
		for (int i=1;i<=anzahl[nr];i++)
		{ /* Summenbildung */
			/* nur die gueltigen Werte werden gezaehlt */
			if (value[nr][i] > -fehlwert)
			{
				n[nr]   += 1;
				sx1i    += xkoord[nr][i];
				sx2i    += ykoord[nr][i];
				syi     += residual[i];
				sx1i2   += xkoord[nr][i]*xkoord[nr][i];
				sx2i2   += ykoord[nr][i]*ykoord[nr][i];
				sx1ix2i += xkoord[nr][i]*ykoord[nr][i];
				sx1iyi  += xkoord[nr][i]*residual[i];
				sx2iyi  += ykoord[nr][i]*residual[i];
			}
		}
		/* Parameterschaetzung von af, bf[1] und bf[2] */
		/* belegung der Vektoren (spaltenweise)*/
		double va[4],vb[4],vc[4];
		va[1] = n[nr];    vb[1] = sx1i;      vc[1] = sx2i;
		va[2] = sx1i;     vb[2] = sx1i2;     vc[2] = sx1ix2i;
		va[3] = sx2i;     vb[3] = sx1ix2i;   vc[3] = sx2i2;
		double det=determinante_drei(va,vb,vc);
		/* Parameter b1 */
		vb[1] = syi;
		vb[2] = sx1iyi;
		vb[3] = sx2iyi;
		bf1[nr+1]=(double)(determinante_drei(va,vb,vc)/det);
		/* Parameter b2 */
		vb[1] = sx1i;
		vb[2] = sx1i2;
		vb[3] = sx1ix2i;
		vc[1] = syi;
		vc[2] = sx1iyi;
		vc[3] = sx2iyi;

		bf2[nr+1]=(double)(determinante_drei(&(va[0]),&(vb[0]),&(vc[0]))/det);
		// parameter a
		af[nr+1] = (double)(( syi - bf1[nr+1]*sx1i - bf2[nr+1]*sx2i ) / n[nr]);
	}
	else
	{ // if less than 3 stations or if a fix lapse rate for the upper regression is defined
		af[nr+1]       = -fehlwert-1;
		bf1[nr+1]    = af[nr+1];
		bf2[nr+1]    = af[nr+1];
	}

	delete[] residual;

}

double meteointerpol::getidw_anisotropisch(int nr,double x,double y,double *entf)
{
	int k;
	double *pe,*pt,*px,*py; 
	double  dx,dy,wert;
	double summeentf;

	summeentf = wert = 0;
	/* schleife je nach Wichtung */
	if (wi==2000)
	{
		pt = &(value[nr][0]);
		px = &xkoord_rot[nr][0];
		py = &ykoord_rot[nr][0];
		//  for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
		for (pe=&(entf[1]);pe<=&(entf[anzahl[nr]]);pe++)
		{
			pt++; px++; py++;
			if (*pt >= -fehlwert)
			{
				dx = (*px-x)*aniso;
				dy = (*py-y);
				*pe = dx*dx + dy*dy;
				if (*pe==0.0) *pe = 1.0;
				if (*pe<=maxentf)
				{
					summeentf += *pe = 1.0/(*pe);
				}
				else *pe = 0;
			}
			else *pe = 0;                                 
		}
	}
	else
	{
		if (wi==1000)
		{
			pt = &(value[nr][0]);
			px = &xkoord_rot[nr][0];
			py = &ykoord_rot[nr][0];
			for (pe=&(entf[1]);pe<=&(entf[anzahl[nr]]);pe++)
			{
				pt++; px++; py++;
				if (*pt >= -fehlwert)
				{
					dx = (*px-x)*aniso;
					dy = (*py-y);
					*pe = dx*dx + dy*dy;
					if (*pe==0.0) *pe = 1.0;
					*pe = sqrt(*pe);
					if (*pe<=maxentf)
					{
						summeentf += *pe = 1.0 / *pe;
					}
					else *pe = 0;
				}
				else *pe = 0;
			}
		}
		else
		{
			pt = &(value[nr][0]);
			px = &xkoord_rot[nr][0];
			py = &ykoord_rot[nr][0];
			for (pe=&(entf[1]);pe<=&(entf[anzahl[nr]]);pe++)
			//      for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
			{
				pt++; px++; py++;
				if (*pt >= -fehlwert)
				{
					dx = (*px-x)*aniso;
					dy = (*py-y);
					*pe = dx*dx + dy*dy;
					if (*pe<=0.0) *pe = 1.0;
					*pe = sqrt(*pe);
					if (*pe<=maxentf)
					{
						summeentf += *pe = pow((double)*pe,(double)-wichtung);
					}
					else *pe = 0;
				}
				else *pe = 0;
			}
		}
	}
	/* da alle Entfernungen reziprokiert wurden, sind sie nun nur dann zugelassen,
	wenn sie groesser als das entsprechend der selben Formel reziprokierte
	maxentf sind. war value[k] ein Fehlwert oder entf[k]>maxentf, dann ist nun
	entf[k]=0, was ein Ausscheiden nach sich zieht*/
	if (summeentf==0) wert = -fehlwert - 1.0; /* waren alle Werte Fehlwerte, so wird ein Fehlwert zurueckgeliefert */
	else
	{
		for (k=1;k<=(anzahl[nr]);k++)
		wert += entf[k]>0.0 ? value[nr][k]*entf[k] : 0.0;
		wert /= summeentf;
	}
	return(wert);
}


double meteointerpol::getidw(int nr,double x,double y,double *entf)
{
	int k;
	double *pe,*pt,*px,*py; 
	double  dx,dy,wert;
	double summeentf;
	summeentf = wert = 0;
	/* schleife je nach Wichtung */
	if (wi==2000)
	{
		pt = &(value[nr][0]);
		px = &(xkoord[nr][0]);
		py = &(ykoord[nr][0]);
		for (pe=&(entf[1]);pe<=&(entf[(anzahl[nr])]);pe++)
		{
			pt++; px++; py++;
			if (*pt >= -fehlwert)
			{
				dx = (*px-x);
				dy = (*py-y);
				*pe = dx*dx + dy*dy;
				if (*pe==0.0) *pe = 1.0;
				if (*pe<=maxentf)
				{
					summeentf += *pe = 1.0/(*pe);
				}
				else *pe = 0;
			}
			else *pe = 0;
		}
	}
	else
	{
		if (wi==1000)
		{
			pt = &(value[nr][0]);
			px = &(xkoord[nr][0]);
			py = &(ykoord[nr][0]);
			for (pe=&(entf[1]);pe<=&(entf[(anzahl[nr])]);pe++)
			//      for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
			{
				pt++; px++; py++;
				if (*pt >= -fehlwert)
				{
					dx = (*px-x);
					dy = (*py-y);
					*pe = dx*dx + dy*dy;
					if (*pe==0.0) *pe = 1.0;
					*pe = sqrt(*pe);
					if (*pe<=maxentf)
					{
						summeentf += *pe = 1.0 / *pe;
					}
					else *pe = 0;
				}
				else *pe = 0;
			}
		}
		else
		{
			pt = &(value[nr][0]);
			px = &(xkoord[nr][0]);
			py = &(ykoord[nr][0]);
			for (pe=&(entf[1]);pe<=&(entf[(anzahl[nr])]);pe++)
			//      for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
			{
				pt++; px++; py++;
				if (*pt >= -fehlwert)
				{
					dx = (*px-x);
					dy = (*py-y);
					*pe = dx*dx + dy*dy;
					if (*pe<=0.0) *pe = 1.0;
					*pe = sqrt(*pe);
					if (*pe<=maxentf)
					{
						summeentf += *pe = pow((double)*pe,(double)-wichtung);
					}
					else *pe = 0;
				}
				else *pe = 0;
			}
		}
	}
	/* da alle Entfernungen reziprokiert wurden, sind sie nun nur dann zugelassen,
	wenn sie groesser als das entsprechend der selben Formel reziprokierte
	maxentf sind. war value[k] ein Fehlwert oder entf[k]>maxentf, dann ist nun
	entf[k]=0, was ein Ausscheiden nach sich zieht*/
	if (summeentf==0) wert = -fehlwert - 1.0; /* waren alle Werte Fehlwerte, so wird ein Fehlwert zurueckgeliefert */
	else
	{
		for (k=1;k<=(anzahl[nr]);k++)
		wert += entf[k]>0.0 ? value[nr][k]*entf[k] : 0.0;
		wert /= summeentf;
	}
	return(wert);
}

double meteointerpol::getthiessen(int nr,double x,double y,double *entf)
{
	double *pe,*pt,*px,*py; 
	double  dx,dy,wert;
	double summeentf;
	/* Wenn in der Interpolation kein richtiger Wert gefunden wird,
	wird der Fehlwert zurueckgeliefert */
	wert = -fehlwert - 1.0;
	summeentf=1e20;
	/* schleife je nach Wichtung */
	pt = &(value[nr][0]);
	px = &(xkoord[nr][0]);
	py = &(ykoord[nr][0]);

	//  for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
	for (pe=&(entf[1]);pe<=&(entf[(anzahl[nr])]);pe++)
	{
		pt++; px++; py++;
		if (*pt >= -fehlwert)
		{
			dx = (*px-x);
			dy = (*py-y);
			*pe = dx*dx + dy*dy;
			if (*pe<=summeentf)
			{
				wert=*pt;
				summeentf=*pe;
			}
		}
	}
	return(wert);
}

double meteointerpol::getThiessenWithLapseRate(double x,double y,double elevation,double lapserate,double upperinversionfix,double upperlapserate)
{
	double  dx,dy,dist,result;
	double mindist=1e20;
	/* Wenn in der Interpolation kein richtiger Wert gefunden wird,
	wird der Fehlwert zurueckgeliefert */
	result = -fehlwert - 1.0;
	for (int i=1;i<=anzahl[0];i++)
	{
		if (value[0][i] >= -fehlwert)
		{
			dx = (xkoord[0][i]-x);
			dy = (ykoord[0][i]-y);
			dist = dx*dx + dy*dy;
			if (dist<mindist)
			{
				if (upperinversionfix>-9998.0)
				{
					if (hoehe[0][i]<=upperinversionfix) // Stationshöhe unter igo
					{
						if (elevation<=upperinversionfix)  // aktuelle Zelle auch unter igo => entspricht Methode 12 ohne upperlapserate
						{
							result=value[0][i]+(elevation-hoehe[0][i])*lapserate;
						}
						else // aktuelle Zelle oberhalb igo
						{
							result=value[0][i]+(upperinversionfix-hoehe[0][i])*lapserate + (elevation-upperinversionfix)*upperlapserate;
						}
					}
					else // stationshöhe über igo
					{
						if (elevation>upperinversionfix)  // aktuelle Zelle auch über igo 
						{
							result=value[0][i]+(elevation-hoehe[0][i])*upperlapserate;
						}
						else // aktuelle Zelle unterhalb igo
						{
							result=value[0][i] + (upperinversionfix-hoehe[0][i])*upperlapserate + (elevation-upperinversionfix)*lapserate;
						}
					}
				}
				else
				{
					result=value[0][i]+(elevation-hoehe[0][i])*lapserate;
				}
				if ((iDataType==PTYPE)&&(value[0][i]==0))
				{
					result=0; // do not change 0 values. Only non-zero precipitation should be changed by lapse rates
				}
				mindist=dist;
			}
		}
	}
	return (result);
}


double meteointerpol::getszenkorr(int i,int j,double x,double y)
/* getszenkorr dient dazu, den fuer den aktuellen Punkt und den aktuellen Zeitpunkt
gueltigen Szenarium-Korrekturwert aus den Klimaszenarien zu extrahieren. Es findet
in dieser Routine keine zeitliche Interpolation statt, das wird in run() erledigt, da die
zeitliche Interpolation fuer alle Rasterpunkte des Interpolationsfeldes gilt. Die zeitlich
interpolierten Szenarienwerte werden jeweils im 16. Element (index 15) der Szenarienarrays gespeichert.
Es wird jeweils der der Interpolationsstelle am naechsten liegende Szenariumwert zurueckgeliefert. */
{
	int k;
	double entf;
	//double *pe,*pt,*px,*py; 
	double  dx,dy,wert;
	double summeentf;

	wert = 0;
	summeentf=1e15;
	/* schleife je nach Wichtung */

	if ((szenmeth==5)||(szenmeth==6))
	{
		wert=SCEN[i][j];
	}
	else
	{
		for (k=1;k<=szennum;k++)
		{
			dx = (szenario[k][0]-x);
			dy = (szenario[k][1]-y);
			entf = dx*dx + dy*dy;
			if (entf<=summeentf)
			{
				wert=szenario[k][15];
				summeentf=entf;
			}
		}
	}
	return(wert);
}

/* liefert bilineare Interpolation fuer Stationen innerhalb eines festen Rasters-> Kopplung mit Klimamodellen, Methode 5 */
double meteointerpol::getbilinear(int nr,int i,int j,double xi,double yi,double *entf)
{
	int k;
	double u,v,x[5],y[5],z[5],d[5],dmin,dsum;
	unsigned int stat[5];
	/* the station indizes are locally ordered in mathematically positive orientation: 
	index 1 is used for the upper right station,
	index 2 is used for the upper left station,
	index 3 is used for the lower left station,
	index 4 is used for the lower reight station,
	*/
	/* getting the indizes of the four stations which are nearest to the actual cell in each of the four quadrants */
	for (k=1;k<=4;k++)
	{
		stat[k]=STAT[k][i][j];
	}
	/* if there is no station in one or more quadrants, or if one or more stations have nodata values,
	the actual station is indexed with a zero to indicate that this station is not available for 
	the actual interpolation */
	for (k=1;k<=4;k++)
	{
		x[k]=xkoord_rot[nr][stat[k]];
		y[k]=ykoord_rot[nr][stat[k]];
		z[k]=value[nr][stat[k]];
		if (z[k]<-fehlwert) stat[k]=0;
	}
	/* if the actual point is on a statio location */
	if ((stat[1]==stat[2])&&(stat[2]==stat[3])&&(stat[3]==stat[4])&&(stat[1]!=0)) return z[1];
	/* if one or more stations are not available for a bilinear interpolation, 
	the inverse distance interpolation without any anisotropy is used instead
	in order to get reasonable interpolation results for all grid cells despite of 
	a temprally lack of data */
	if ((stat[1]==0)||(stat[2]==0)||(stat[3]==0)||(stat[4]==0))
	return getidw_anisotropisch(nr,xi,yi,entf);
	/* the distances in x and y direction in a regular mesh */
	d[1]=x[1]-x[2];
	d[2]=y[3]-y[2];
	/* if, for any reason, one of these distances is zero, the 
	idw-interpolation is called. Since the stat-index of at least one station
	should be zero, this statement should never be reached, but is is surely safer this way 
	(e.g. if future program changes changes the codes of stations) */
	if ((d[1]==0)||(d[2]==0)) getidw_anisotropisch(nr,xi,yi,entf);
	u = (xi-x[2])/d[1]; /* weight of data of stations 1 and 4 in the "horizontal" x direction (stations 2 and 3 have the weight 1-u) */
	v = (yi-y[2])/d[2]; /* weight of data of stations 3 and 4 in the "vertical" y direction (stations 1 and 2 have the weight 1-v) */

	/* the actual interpolation */
	z[0] = (1.0-u)*(1.0-v)*z[2]+
	(    u)*(1.0-v)*z[1]+
	(1.0-u)*(    v)*z[3]+
	(    u)*(    v)*z[4];

	return (double)z[0];
}

int meteointerpol::StationsExist(int *p_stations,int p_anzahl)
{
	int i;
	for (i=1;i<=p_anzahl;i++)
	{
		if (p_stations[i]!=0)
		{
			return 1;
		}
	}
	return 0;
}

typedef std::vector<int> vect;

void meteointerpol::ConstructGridFromStationCoordinates(int nr)
{
	int i,j,actstat,actrow,startstat;
	vector<int>::iterator stat_iter;
	vector<double> mean_x(0);
	vector<double>::iterator mean_x_iter;
	double xmin,xmean,ymean,xsum,ysum;
	vector<vect> StatArray(0);

	StatMeshNumCols=0;
	station=new int[anzahl[0]+1];
	// a station which is linked to the list will get a zero-flag, 
	// i.e., the station may be linked while it's flag is 1
	for (i=1;i<=anzahl[0];i++) station[i]=1; 
	do
	{
		/* Station mit kleinster x-Koordinate */
		xmin=1e40;
		for (i=1;i<=anzahl[0];i++)
		{
			if ((xkoord_rot[nr][i]<xmin)&&(station[i]>0))
			{
				xmin=xkoord_rot[nr][i];
				startstat=i;
			}
		}
		xmean=ymean=xsum=ysum=0;
		actstat=startstat;
		StatArray.push_back(vector<int>(0)); // new column with actually one row
		j=StatArray.size()-1;
		actrow=-1;
		// searching downwards first
		do                                 
		{
			actrow++;
			station[actstat]=0;
			StatArray[j].push_back(actstat);
			xsum+=xkoord_rot[nr][actstat];
			actstat=GetNextStation(xkoord_rot[nr][actstat],ykoord_rot[nr][actstat],
			-AvgHorDist*(0.0+delta) , AvgHorDist*(0.0+delta),
			-AvgVertDist*(1.0+delta),-AvgVertDist*(1.0-delta));
		} 
		while ((actstat!=0) && (station[actstat]!=0));

		// now searching upwards
		actrow=0;
		
		//starting from the same station as the downwards search loop before did
		actstat=GetNextStation(xkoord_rot[nr][startstat],ykoord_rot[nr][startstat],
		-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
		AvgVertDist*(1.0-delta),AvgVertDist*(1.0+delta));
		while ((actstat!=0) && (station[actstat]!=0))
		{
			stat_iter = StatArray[j].begin();
			StatArray[j].insert(stat_iter,actstat);
			xsum+=xkoord_rot[nr][actstat];
			station[actstat]=0;
			actstat=GetNextStation(xkoord_rot[nr][actstat],ykoord_rot[nr][actstat],
			-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
			AvgVertDist*(1.0-delta),AvgVertDist*(1.0+delta));
		} 

		mean_x.push_back(xsum/(double)(StatArray[StatArray.size()-1].size()));
		mean_x_iter=mean_x.begin();
	} 
	while (StationsExist(station,anzahl[0])==1);

	// generating a temporary Grid according to the number of rows and columns in the statoin mesh
	StatMeshNumRows=0;
	for (vector<vect>::iterator p=StatArray.begin();p<StatArray.end();p++)
	{
		if ((*p).size()!=StatArray[0].size())
		{
			(*WSstream) << "FATAL: Error converting the station coordinates into a grid:\n"
			"The actual version of WaSiM-ETH accepts only station meshs\n"
			"with an equal number of more or less equally spaced stations\n"
			"in each column, whereas the total number of both, rows \n"
			"and columns, are each virtually unlimited\n" << "[end-entry]";
			throw(70);
		}
	}
	StatMeshNumRows=StatArray[0].size();
	StatMeshNumCols=StatArray.size();

	StationIndexGrid = matrix(1,StatMeshNumRows,1,StatMeshNumCols);
	StationGrid      = matrix(1,StatMeshNumRows,1,StatMeshNumCols);
	SecondDerivative = matrix(1,StatMeshNumRows,1,StatMeshNumCols);
	CoordXStatR      = matrix(1,StatMeshNumRows,1,StatMeshNumCols);
	CoordYStatR      = matrix(1,StatMeshNumRows,1,StatMeshNumCols);
	/* hier müssen dann die Koordinaten entsprechend der Spalten/Reihenzuordnung der Stationen aus der Liste rein! */
	for (i=1;i<=StatMeshNumRows;i++)
	for (j=1;j<=StatMeshNumCols;j++)
	{
		CoordXStatR[i][j]=xkoord_rot[0][StatArray[j-1][i-1]];
		CoordYStatR[i][j]=ykoord_rot[0][StatArray[j-1][i-1]];
		StationIndexGrid[i][j]=StatArray[j-1][i-1];
	}
	if (methode==8) result2=matrix(mpistartrow,mpiendrow,1,ncols,nrows);

	/* constructor for the bicubic spline class */
	Spline = new TSplines(CoordXStatR,CoordYStatR,CoordXR,CoordYR,
	StatMeshNumRows,StatMeshNumCols,CoordXStatR[StatMeshNumRows][1],CoordYStatR[StatMeshNumRows][1],
	AvgHorDist,AvgVertDist,-9999.0,nrows,ncols,xll,yll,cellsize,nodata,anisoneig,1,0,zn);
}

void meteointerpol::GetMeanStationDist(int nr)
{
	int i,j,s;

	s=0; /* only for station file 0 required, because all other files have to have identical stations */
	delta=0.25;

	(*WSstream) << "calculating average horizontal and vertical distance to the next station..."<<"[end-entry]";
	for (i=1;i<=anzahl[s];i++)
	{ // the search for the nearest station is performed for each station 
		dist=dhmin[i]=dvmin[i]=-9999;
		for (j=1;j<=anzahl[s];j++)
		{
			if (i!=j)
			{
				/* distance to the actual station as square */
				//        if (((xkoord_rot[s][j]-xkoord_rot[s][i])>0)&& 
				//            ((ykoord_rot[s][j]-ykoord_rot[s][i])>0))   /* only stations to the rightand to the upper direction are considered */
				{
					dist=sqrt((xkoord_rot[s][j]-xkoord_rot[s][i])*(xkoord_rot[s][j]-xkoord_rot[s][i])+
					(ykoord_rot[s][j]-ykoord_rot[s][i])*(ykoord_rot[s][j]-ykoord_rot[s][i]));
					
					if (fabs(xkoord_rot[s][j]-xkoord_rot[s][i])< /* station in south nord direction */
							fabs(ykoord_rot[s][j]-ykoord_rot[s][i])) 
					{
						if (((ykoord_rot[s][j]-ykoord_rot[s][i])>0)&&((dist<dvmin[i])||(dvmin[i]<-1000)))
						dvmin[i]=dist;
					}
					else /* station in west east direction */
					{
						if (((xkoord_rot[s][j]-xkoord_rot[s][i])>0)&&((dist<dhmin[i])||(dhmin[i]<-1000)))
						dhmin[i]=dist;
					}
				}   
			}
		}
	}

	/* getting the average horizontal and vertical distances to the nearest stations */
	anzdhmin=anzdvmin=sumdhmin=sumdvmin=0.0;
	for (i=1;i<=anzahl[s];i++)
	{
		if (dhmin[i]>-1000) /* if there is a neighbouring station in the east */
		{
			anzdhmin+=1.0;
			sumdhmin+=dhmin[i];
		}
		if (dvmin[i]>-1000) /* if there is a neighbouring station in the north */
		{
			anzdvmin+=1.0;
			sumdvmin+=dvmin[i];
		}
	}

	/* now getting the mean distance between 2 stations in East West direction */
	if (anzdhmin>0)
	AvgHorDist=sumdhmin/anzdhmin;
	else 
	AvgHorDist=-9999;

	/* now getting the mean distance between 2 stations in North-South direction */
	if (anzdvmin>0)
	AvgVertDist=sumdvmin/anzdvmin;
	else 
	AvgVertDist=-9999;

	(*WSstream) << "avg. horizontal distance: "<<AvgHorDist<<" m"<<"[end-entry]";
	(*WSstream) << "avg. vertical   distance: "<<AvgVertDist<<" m"<<"[end-entry]";
}

void meteointerpol::MakeBilinearStationGrid(int nr)
{
	int i,j,k,minstat,s;
	double u,v,x[5],y[5],z[5],d[5],dmin,dsum;
	unsigned int stat[5];
	double  dx,dy;

	if (MPIREG == NULL)
	{
		if (ist_vorhanden(_MPIRegions_))
		{
			MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
		}
		else
		{
			MPIREG = NULL;
		}
	}

	/* generating the grids for station number look-up */
	if ((methode==5)||(methode==6))
	{
		STAT=(unsigned int ***)calloc(1+4*(nr+1),sizeof(unsigned int **));
		for (k=1;k<=4*(nr+1);k++)
		{
			STAT[k]=(unsigned int **)calloc(nrows+1,sizeof(unsigned int *));
			for (int i = 1; i < mpistartrow; i++)
			{
				STAT[k][i] = NULL;
			}
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				STAT[k][i]=(unsigned int*)calloc(ncols+1,sizeof(unsigned int));
				for (j=1;j<=ncols;j++)
				STAT[k][i][j]=0;
			}
			for (int i = mpiendrow + 1; i <= nrows; i++)
			{
				STAT[k][i] = NULL;
			}
		}
		(*WSstream) << "choosing look-up-stations...";
		/* now making the station look-up grids */

		for (s=0;s<nr;s++)
		{
			for (i=mpistartrow;i<=mpiendrow;i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j]==mpi_rank))
					{
						dmin=dxmin=dymin=1e30;
						/* looking for the nearest station */
						for (k=1;k<=(anzahl[s]);k++)
						{                                
							dx=xkoord_rot[s][k]-CoordXR[i][j];
							dy=ykoord_rot[s][k]-CoordYR[i][j];
							dist = sqrt(dx*dx+dy*dy);
							if (dist<dmin) 
							{ 
								dmin=dist;
								dxmin=dx;
								dymin=dy;
								minstat=k;
							}
						}
						dx=dxmin;
						dy=dymin;
						k=minstat;
						/* now the nearest station was found, so the next step is 
			to look for the other 3 stations in the station grid, 
			depending on the quadrant */
						if (dx>0)
						{ /* quadrants 1 and 4 */
							if (dy>0)
							{  /* quadrant 1 */
								STAT[4*s+1][i][j]=k; 
								/* looking for a station to the west of station 1 */
								STAT[4*s+2][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(1.0+delta),-AvgHorDist*(1.0-delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
								/* looking for a station to the south of station 2 */
								STAT[4*s+3][i][j]=GetNextStation(xkoord_rot[s][STAT[4*s+2][i][j]],ykoord_rot[s][STAT[4*s+2][i][j]],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								-AvgVertDist*(1.0+delta),-AvgVertDist*(1.0-delta));
								/* looking for a station to the south of station 1 */
								STAT[4*s+4][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								-AvgVertDist*(1.0+delta),-AvgVertDist*(1.0-delta));
							}
							else
							{  /* dy<0: quadrant 4 */
								STAT[4*s+4][i][j]=k; 
								/* looking for a station to the west of station 1 */
								STAT[4*s+3][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(1.0+delta),-AvgHorDist*(1.0-delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
								/* looking for a station to the north of station 3 */
								STAT[4*s+2][i][j]=GetNextStation(xkoord_rot[s][STAT[4*s+3][i][j]],ykoord_rot[s][STAT[4*s+3][i][j]],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								AvgVertDist*(1.0-delta),AvgVertDist*(1.0+delta));
								/* looking for a station to the north of station 4 */
								STAT[4*s+1][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								AvgVertDist*(1.0-delta),AvgVertDist*(1.0+delta));
							}
						}
						else
						{ /* dx<0 quadrants 2 and 3 */
							if (dy>0)
							{  /* quadrant 2 */
								STAT[4*s+2][i][j]=k; 
								/* looking for a station to the south of station 2 */
								STAT[4*s+3][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								-AvgVertDist*(1.0+delta),-AvgVertDist*(1.0-delta));
								/* looking for a station to the east of station 3 */
								STAT[4*s+4][i][j]=GetNextStation(xkoord_rot[s][STAT[4*s+3][i][j]],ykoord_rot[s][STAT[4*s+3][i][j]],
								AvgHorDist*(1.0-delta),AvgHorDist*(1.0+delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
								/* looking for a station to the east of station 2 */
								STAT[4*s+1][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								AvgHorDist*(1.0-delta),AvgHorDist*(1.0+delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
							}
							else
							{  /* quadrant 3 */
								STAT[4*s+3][i][j]=k; 
								/* looking for a station to the north of station 3 */
								STAT[4*s+2][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								-AvgHorDist*(0.0+delta),AvgHorDist*(0.0+delta),
								AvgVertDist*(1.0-delta),AvgVertDist*(1.0+delta));
								/* looking for a station to the east of station 2 */
								STAT[4*s+1][i][j]=GetNextStation(xkoord_rot[s][STAT[4*s+2][i][j]],ykoord_rot[s][STAT[4*s+2][i][j]],
								AvgHorDist*(1.0-delta),AvgHorDist*(1.0+delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
								/* looking for a station to the east of station 3 */
								STAT[4*s+4][i][j]=GetNextStation(xkoord_rot[s][k],ykoord_rot[s][k],
								AvgHorDist*(1.0-delta),AvgHorDist*(1.0+delta),
								-AvgVertDist*(0.0+delta),AvgVertDist*(0.0+delta));
							}
						}
						/* if there is no station at all (should never be reached because 
			that would cause an error earlier in the meteoint-constructor) */
						if ((STAT[4*s+1][i][j]&&0)||(STAT[4*s+2][i][j]&&0)||
								(STAT[4*s+3][i][j]&&0)||(STAT[4*s+4][i][j]&&0))
						{
							(*WSstream) << "FATAL: \n---------------------- error ------------------\n"
							"in data file "<<(anetzstat[s])<<":\n"
							"For bilinear interpolation the entire model domain\n"
							"must be covered by a mesh of stations/UCM-outputs!\n"
							"The cell at row "<<i<<", column "<<j<<" is outside of such a mesh!\n"
							"-------------------------------------------------\n" << "[end-entry]";
							throw(67);
						}
					}
				}
				(*WSstream) << ".";
			}
			(*WSstream) << "[end-entry]";
		}

#ifdef _WASIMDEBUG_
		/* writing the station grids */
		FILE *statfile;
		statfile=fopen("stat1.asc","wt");
		fprintf(statfile,"ncols       %d\n",ncols);
		fprintf(statfile,"nrows       %d\n",nrows);
		fprintf(statfile,"xll         %g\n",xll);
		fprintf(statfile,"yll         %g\n",yll);
		fprintf(statfile,"cellsize    %g\n",cellsize);
		fprintf(statfile,"nodatavalue %g\n",nodata);
		for (int i = 1; i < mpistartrow; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		for (int i=mpistartrow;i<=mpiendrow;i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				fprintf(statfile,"%g\t",(double)STAT[1][i][j]);
			}
			fprintf(statfile,"\n");
		}
		for (int i = mpiendrow+1; i < nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		fclose(statfile);

		statfile=fopen("stat2.asc","wt");
		fprintf(statfile,"ncols       %d\n",ncols);
		fprintf(statfile,"nrows       %d\n",nrows);
		fprintf(statfile,"xll         %g\n",xll);
		fprintf(statfile,"yll         %g\n",yll);
		fprintf(statfile,"cellsize    %g\n",cellsize);
		fprintf(statfile,"nodatavalue %g\n",nodata);
		for (int i = 1; i < mpistartrow; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				fprintf(statfile,"%g\t",(double)STAT[2][i][j]);
			}
			fprintf(statfile,"\n");
		}
		for (int i = mpiendrow + 1; i < nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		fclose(statfile);

		statfile=fopen("stat3.asc","wt");
		fprintf(statfile,"ncols       %d\n",ncols);
		fprintf(statfile,"nrows       %d\n",nrows);
		fprintf(statfile,"xll         %g\n",xll);
		fprintf(statfile,"yll         %g\n",yll);
		fprintf(statfile,"cellsize    %g\n",cellsize);
		fprintf(statfile,"nodatavalue %g\n",nodata);
		for (int i = 1; i < mpistartrow; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		for (int i = 1; i <= nrows; i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				fprintf(statfile,"%g\t",(double)STAT[3][i][j]);
			}
			fprintf(statfile,"\n");
		}
		for (int i = mpiendrow + 1; i < nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		fclose(statfile);

		statfile=fopen("stat4.asc","wt");
		fprintf(statfile,"ncols       %d\n",ncols);
		fprintf(statfile,"nrows       %d\n",nrows);
		fprintf(statfile,"xll         %g\n",xll);
		fprintf(statfile,"yll         %g\n",yll);
		fprintf(statfile,"cellsize    %g\n",cellsize);
		fprintf(statfile,"nodatavalue %g\n",nodata);
		for (int i = 1; i < mpistartrow; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		for (int i = 1; i <= nrows; i++)
		{
			for (int j=1;j<=ncols;j++)
			{
				fprintf(statfile,"%g\t",(double)STAT[4][i][j]);
			}
			fprintf(statfile,"\n");
		}
		for (int i = mpiendrow + 1; i < nrows; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				fprintf(statfile, "%g\t", nodata);
			}
			fprintf(statfile, "\n");
		}
		fclose(statfile);
#endif
	}
}

int meteointerpol::GetNextStation(double x,double y,double dxmin,double dxmax,double dymin,double dymax)
{
	int i,s;
	s=0;

	for (i=1;i<=anzahl[s];i++)
	{
		if (((x+dxmin)<xkoord_rot[s][i])&&((x+dxmax)>xkoord_rot[s][i]) && /* station within x range */
				((y+dymin)<ykoord_rot[s][i])&&((y+dymax)>ykoord_rot[s][i]))
		{
			return i;
		}
	}

	return 0;
}

/* this function calculates the skewness (anisotropie slope) of the 
station grid for bilinear interpolation on a statistical basis.
the result is a single scalar value which is the average anisotropie 
slope of the station mesh when tried to constructed automatically */
void meteointerpol::GetBilinearSkewness()
{
	int i,j,s,MostFreqSlp,MostFreqSlpClass,absminstation;
	double dist,absmin;
	double fk,anz,sum,val;

	if ((methode==5)||(methode==6)||(methode==7)||(methode==8))
	{
		s=0; 
		/* anisotropy slope is calculated only for the stations in grid 0, since all 
	other stations in the other files are identical to those in file 0 */
		for (i=1;i<=anzahl[s];i++)
		{ // the search for the nearest station is performed for each station 
			dist=mindist[i]=minstation[i]=-9999;
			for (j=1;j<=anzahl[s];j++)
			{
				if (i!=j)
				{
					/* distance to the actual station as square */
					dist=(xkoord[s][j]-xkoord[s][i])*(xkoord[s][j]-xkoord[s][i])+
					(ykoord[s][j]-ykoord[s][i])*(ykoord[s][j]-ykoord[s][i]);
					if (mindist[i]<0) 
					{
						absmin=dist;
						absminstation=j;
						if (minstation[j]!=i)
						{
							mindist[i]=dist;
							minstation[i]=j;
						}
					}
					else
					{
						if (dist<mindist[i])
						{
							absmin=dist;
							absminstation=j;
						}
						if ((dist<mindist[i])&&(minstation[j]!=i))
						{ /* if for station j the actual station is already 
				found as nearest station, another station should be taken */
							mindist[i]=dist;
							minstation[i]=j;
						}
					}
				}
			}
			if (mindist[i]>0) 
			{
				/* if the found station was already cross related to the actual station 
		(i.e., the actual station was already the nearest station for the found station),
		then another station was looked for. But if this other station is too far away, 
		then the cross related station should be the nearest station again. This mens,
		if the nearest station was found to the left of the actual station but this station was
		already cross related to the actual station, then usuallay the station to the right of
		the actual station would be in the same distance and hence another good candidate
		as the nearest station. But if there is no other station to the right, and if also
		the grid size in the other direction is much larger, then the compromise must be taken. */
				if (mindist[i]>1.01*absmin)
				{
					mindist[i]=absmin;
					minstation[i]=absminstation;
				}
				mindist[i]=sqrt(mindist[i]);
			}
		}
		/* now getting the direction */
		for (i=1;i<=anzahl[s];i++)
		{
			direction[i]=atan2((double)(ykoord[s][minstation[i]]-ykoord[s][i]),
			(double)(xkoord[s][minstation[i]]-xkoord[s][i]));
			/* only directions between -90 and +90 are allowed -> direction is not a vector, so e.g. 120=-60! */
			if (direction[i]>(Pi/2.0)) direction[i]-=Pi;
			if (direction[i]<(-Pi/2.0)) direction[i]+=Pi;

			if ((direction[i]*180.0/Pi)>355.0)
			slopeclass[0]++;
			else
			{
				fk=direction[i]*180.0/Pi+365.0;
				fk=fmod(fk,360.0);
				fk=fmod((double)((int)(fk/10)),18.0);
				slopeclass[(int)(fk)]++;
			}
		}
		/* getting the most frequent slope class */
		MostFreqSlp=MostFreqSlpClass=0;
		for (i=0;i<=35;i++)
		{
			if (slopeclass[i]>MostFreqSlp)
			{
				MostFreqSlp=slopeclass[i];
				MostFreqSlpClass=i;
			}
		}
		anz=sum=0.0;
		/* getting the mean slope in that slope class */
		for (i=1;i<=anzahl[s];i++)
		{
			fk=direction[i]*180.0/Pi+365.0;
			fk=fmod(fk,360.0);
			fk=fmod((double)((int)(fk/10)),18.0);
			//      fk=(int)(fk/10);
			if (MostFreqSlpClass==fk)
			{
				val=fmod(direction[i]*180.0/Pi+360.0,180.0);
				if ((val > 5.0)&&((int)fk==0))
				{
					val-=180.0;
				}
				else
				{ 
					if ((val < -5.0)&&((int)fk==0))
					val+=180.0;
				}
				if (((fk*10-5.0)<=val)&&((fk*10.0+5.0)>=val))
				{
					anz++;
					sum+=val;
				}
			}
		}

		anisoneig=anz>0?sum/anz:0;
		/* now the anisotropy slope should be converted into a range between -45 
	degree to +45 degree */
		if (anisoneig>45)
		{
			if (anisoneig>45)
			{ /* >135 ... 180 */
				anisoneig-=90.0;
			}
			else
			{ /* >45 ... 135 */
				anisoneig-=180.0;
			}
		}
		/* aniso is set to a value less than 1 in order to assure that the anisotropic versions of getidw
	are called (even if its not needed yet, but may be later) */
		if (fabs(anisoneig)>1e-3) 
		{
			aniso=0.995; 
		}
		else 
		{
			aniso=1.0; 
			anisoneig=0.0;
		}
	}
}

void meteointerpol::RotateCoordGrids()
{
	int j,i;
	char *CoordXRName,*CoordYRName;
	double entf;
	double  dx,dy;

	if ((aniso<1)&&((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12)))
	{ /* generating rotated grids of CoordX and CoordY called then CoordXR and CoordYR */
		j=1;
		while (gridinfo[j].identification!=NULL)
		{ 
			j++;
			if (j==gridmaxanzahl)
			{ 
				(*WSstream) << "FATAL: no more free grids in memeory, (max. "<<gridmaxanzahl<<")!"<<"[end-entry]";  
				throw(12); 
			}
		}
		coordxrindex=j;
		coordyrindex=j+1;
		CoordXRName=newstring(_stringsize_+1,1);
		CoordYRName=newstring(_stringsize_+1,1);

		sprintf(CoordXRName,"CoordXRot_%s",instanzname);
		CoordXR = matrix(mpistartrow,mpiendrow,1,ncols,nrows);   
		initialize_grid(coordxrindex,CoordXRName,CoordXR,NULL,0,CoordXRName,0,mpistartrow,mpiendrow);

		sprintf(CoordYRName,"CoordYRot_%s",instanzname);
		CoordYR = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(coordyrindex, CoordYRName, CoordYR, NULL, 0, CoordYRName, 0, mpistartrow, mpiendrow);

		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j=1;j<=ncols;j++) 
			{
				if (CoordX[i][j]!=nodata)
				{
					/* berechnen der Entfernung */
					dx = (CoordX[i][j]-xll);
					dy = (CoordY[i][j]-yll);
					entf = sqrt(dx*dx+dy*dy);
					/* actual angle of the points vector wiht respect to xll and yll */
					if (dx==0.0)
					alpha = dy > 0.0 ? Pi/2.0 : -Pi/2.0;
					else alpha = atan2(dy,dx);
					/* Winkel nach rechts drehen um Aniotropiewinkel */
					beta = alpha - anisoneig;
					/* neue Koordinaten */
					CoordXR[i][j] = xll+cos(beta)*entf;
					CoordYR[i][j] = yll+sin(beta)*entf;
				}
				else
				{
					CoordXR[i][j]=nodata;
					CoordYR[i][j]=nodata;
				}
			}
		}
		if (CoordXRName) free(CoordXRName); CoordXRName=NULL;
		if (CoordYRName) free(CoordYRName); CoordYRName=NULL;
	}
}

// rotating the coordinates of meteo stations (is done once for all possible files 
//   (if method == 6), since the stations have to be identical in all files
void meteointerpol::RotateStationCoordiantes()
{
	int i,s;
	double  dx,dy;
	wi = (int)floor(1000.0*wichtung);

	/* Wenn eine Anisotropie im IDW-Verfahren beruecksichtigt werden soll, dann muessen die Stationen
	und die Koordinaten aller Gridpunkte nun um den Winkel der Neigung der Anisotropieellipsenhauotachse
	gedreht werden. Dieser Winkel anisoneig muss zwischen -90 und +90 Grad liegen. Ob eine Anistotropie vorliegt,
	wird anhand aniso, dem Achsenverhaeltnis der Anisotropieellipse bestimt, welches 1.0 ist, wenn keine
	Anisotropie vorliegt */
	xkoord_rot = (double**)calloc(statfilecnt+1,sizeof(double*)); if (!xkoord_rot) speichererror(1);
	ykoord_rot = (double**)calloc(statfilecnt+1,sizeof(double*)); if (!ykoord_rot) speichererror(1);
	for (s=0;s<statfilecnt;s++)
	{
		xkoord_rot[s]=NULL;
		ykoord_rot[s]=NULL;
		if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12))
		{
			if ((anisoneig < -(Pi/2.0))||(anisoneig > (Pi/2.0)))
			{
				(*WSstream) << "FATAL: anisotropy slope has to range from -90 to +90 degree \n" << "[end-entry]";
				throw(54);
			}
			if ((aniso<=0.0)||(aniso>1.0))
			{
				(*WSstream) << "FATAL: ratio of anisotropy has to range from +1 to >0 liegen"<<"[end-entry]";
				throw(54);
			}
			if (aniso<1.0) /* nun erzeugen der beiden Koordinatenvektoren */
			{
				xkoord_rot[s] = (double *)calloc(anzahl[s]+1,sizeof(double)); if (!xkoord_rot[s]) speichererror(1);
				ykoord_rot[s] = (double *)calloc(anzahl[s]+1,sizeof(double)); if (!ykoord_rot[s]) speichererror(1);
				/* Drehen der Koordinaten der Stationen */
				/*      double yol=yll+cellsize*(double)(nrows-1); */
				for (i=1;i<=anzahl[s];i++)
				{
					/* berechnen der Entfernung */
					dx = (xkoord[s][i]-xll);
					dy = (ykoord[s][i]-yll);
					entf[s][i] = sqrt(dx*dx+dy*dy);
					/* aktueller Winkel der Station bezogen auf xll,yll */
					if (dx==0.0)
					alpha = dy > .0 ? Pi/2.0 : -Pi/2.0;
					else alpha = atan2(dy,dx);
					/* Winkel nach rechts drehen um Aniotropiewinkel */
					beta = alpha - anisoneig;
					/* neue Koordinaten */
					xkoord_rot[s][i] = xll+cos(beta)*entf[s][i];
					ykoord_rot[s][i] = yll+sin(beta)*entf[s][i];
				}
			}
			else
			{
				xkoord_rot[s]=xkoord[s];
				ykoord_rot[s]=ykoord[s];
			}
		}
	}
}

void meteointerpol::run(long vgl_minute,long vgl_stunde,long vgl_tag,long vgl_monat,long vgl_jahr)
{
	double  wert,wert2;
	int    i,j,k,mi;
	double gradient; /* xcolumne, yrow */
	double  szenkorr;  /* der fuer einen Punkt zeitlich und raeumlich aus den Szenariowerten interpolierte Wert */
	bool NewGrid = true;

	if ((jahr!=vgl_jahr)||(monat!=vgl_monat)||(tag!=vgl_tag)||(stunde!=vgl_stunde) || (minute != vgl_minute))
	{
		(*WSstream) << "FATAL: -----------------------------------------------------------------------------\n"
		"error in date columns of "<<instanzname<<"\n\n"
		"first Meteo data set (this is the set from which data is get:\n"
		<<vgl_tag<<"."<<vgl_monat<<"."<<vgl_jahr<<" "<<vgl_stunde<<":"<<vgl_minute<<" h\n\n"
		"but in "<<instanzname<<":\n"
		<<tag<<"."<<monat<<"."<<jahr<<" "<<stunde<<":"<<minute<<" h\n\n"
		"-----------------------------------------------------------------------------\n" << "[end-entry]";
		throw(52);
	}

	if (MPIREG == NULL)
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

	//  LogEntry(1,"Test");

#ifdef _OPENMP
	//#pragma omp critical
	if (!omp_in_parallel())
	{
#endif
		if (nachricht == 1)
		{
			if (szenuse == 1)
			{
				if ((timestep > 0) && (timestep < 60))
				{
					(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ":" << minute << " h " << _Scenarios_ << "[end-entry]";
				}
				else
				{
					(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << _Scenarios_ << "[end-entry]";
				}
			}
			else
			{
				if ((timestep > 0) && (timestep < 60))
				{
					(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ":" << minute << " h " << "[end-entry]";
				}
				else
				{
					(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << "[end-entry]";
				}
			}
		}
		else
		{
			if (((60 * stunde + (int)((minute + 60.0*hour_offset))) % (60 * nachricht) == 0) && (stunde != 0))
			// if ((((stunde + (int)hour_offset) % nachricht == 0) && (stunde > 0)) || ((minute + timestep == 60) && (stunde == 23)))
			{
				if (szenuse == 1)
				{
					if ((timestep > 0) && (timestep < 60))
					{
						(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ":" << minute << " h " << _Scenarios_ << "[end-entry]";
					}
					else
					{
						(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << _Scenarios_ << "[end-entry]";
					}
				}
				else
				{
					if ((timestep > 0) && (timestep < 60))
					{
						(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ":" << minute << " h " << "[end-entry]";;
					}
					else
					{
						(*WSstream) << wasim::setw(20) << instanzname << " " << tag << "." << monat << "." << jahr << " " << stunde << ". h " << "[end-entry]";;
					}
				}
			}
		}
#ifdef _OPENMP
		//#pragma omp critical
	}
#endif


	/* Zeitliche Interpolation der Szenariodaten */
	if (szenuse==1)
	{
		if ((szenmeth==5)||(szenmeth==6))
		{
			if (scenmonth!=monat)
			{
				scenmonth=monat;
				if (scenmonth>12) scenmonth=(scenmonth-1)%12+1;
				free_matrix(gridinfo[scenindex].pointer, mpistartrow, mpiendrow, 1, ncols, nrows);
				delete_grid(scenindex);
				// grid is deleted every time and re-initialized
				readgridfromfile(scengridname[scenmonth], &SCEN, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);  /* includes filling of gaps -> paramter is 1 */
				sprintf(scengrididentifier,"scenariogrid_%s",instanzname);
				initialize_grid(scenindex, scengridname[scenmonth], SCEN, NULL, 0, scengrididentifier, 0, mpistartrow, mpiendrow);
			}
		}
		else
		{
			if (tag<=15)
			{
				j=((int)monat+10)%12+3;
				k=((int)monat+11)%12+3;
				for (i=1;i<=szennum;i++)
				{
					szenario[i][15]=(double)(mon[j-2]-16+tag)/mon[j-2]*szenario[i][k] + (double)(16-tag)/mon[j-2]*szenario[i][j];
				}
			}
			else
			{
				j=((int)monat+11)%12+3;
				k=((int)monat+12)%12+3;
				for (i=1;i<=szennum;i++)
				{
					szenario[i][15]=(double)(tag-15)/mon[j-2]*szenario[i][k] + (double)(mon[j-2]-tag+15)/mon[j-2]*szenario[i][j];
				}
			}
		}
	}

	/* doing the spline interpolation (must be done for the entire Grid in a single call) */
	if (methode==7)
	{
		/* wird spaeter anhand der Mapping-Liste Stationen->Stationsgrid  geregelt */
		for (i=1;i<=StatMeshNumRows;i++)
		{
			for (j=1;j<=StatMeshNumCols;j++)
			{
				StationGrid[i][j]=value[0][(int)StationIndexGrid[i][j]];
			}
		}
		Spline->BicubicColumnSpline(StationGrid,SecondDerivative);
		Spline->BicubicRowSpline(StationGrid,SecondDerivative,result);
	}

	/* doing the spline interpolation for linear combinations of two grids (must be done for the entire Grid in a single call) */
	if (methode==8)
	{

		/* wird spaeter anhand der Mapping-Liste Stationen->Stationsgrid  geregelt */
		for (i=1;i<=StatMeshNumRows;i++)
		{
			for (j=1;j<=StatMeshNumCols;j++)
			{
				StationGrid[i][j]=value[0][(int)StationIndexGrid[i][j]];
			}
		}
		//    for (i=1;i<=10;i++)
		//      for (j=1;j<=10;j++)
		//        StationGrid[11-i][j]=value[0][10*(i-1)+j];
		
		Spline->BicubicColumnSpline(StationGrid,SecondDerivative);
		Spline->BicubicRowSpline(StationGrid,SecondDerivative,result);

		/* wird spaeter anhand der Mapping-Liste Stationen->Stationsgrid  geregelt */
		for (i=1;i<=StatMeshNumRows;i++)
		{
			for (j=1;j<=StatMeshNumCols;j++)
			{
				StationGrid[i][j]=value[1][(int)StationIndexGrid[i][j]];
			}
		}
		//    for (i=1;i<=10;i++)
		//      for (j=1;j<=10;j++)
		//        StationGrid[11-i][j]=value[1][10*(i-1)+j];
		
		Spline->BicubicColumnSpline(StationGrid,SecondDerivative);
		Spline->BicubicRowSpline(StationGrid,SecondDerivative,result2);
	}

	// read grids according to a gridlist
	if (methode==9)
	{
		SetSignalDefaults();
		if (strcmp(gridlist_actname,"-9999")!=0)
		{ // if there is no grid, the old grid is still valid and nothing is read in
			readgridfromfile(this->gridlist_actname, &result, 1, _AllocateGridMemoryOff_, _ShowMessagesOff_, 0, _LimitAutoFillOff_, _MPIFullGridOn_, _MPIRank0OnlyOff_,_MPISyncHaloOff_,true);  // includes filling of gaps without limits -> fillmiss paramter is 1 and last parameter is set to false (limitAutoFill) 
			strcpy(gridlist_lastvalidname, gridlist_actname);
		}
		else
		{
			if (strcmp(gridlist_lastvalidname, "") != 0) // if there was a valid gridname already, use that one again
			{
				readgridfromfile(this->gridlist_lastvalidname, &result, 1, _AllocateGridMemoryOff_, _ShowMessagesOff_, 0, _LimitAutoFillOff_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOff_, true);  // includes filling of gaps without limits -> fillmiss paramter is 1 and last parameter is set to false (limitAutoFill) 
			}
			else // otherwise throw an error
			{
				(*WSstream) << "FATAL: interpolation method 9 (read from grid) no valid grid name was found (-9999 as input) and no last valid grid name present (i.e., simulation starts, when input list file contains -9999 only)" << "[end-entry]";
				throw(1050);
			}
		}
		SetSignalAbfangfunktion();
	}

	/* --> this is done in modapuzr already because for precipitation correction the parameters must be known, so we can skip this step here
// estimate regression parameters as those parameters are not read in in methods 10 and 11
if ((methode==10)||(methode==11))
{
	CalculateRegressionParameters(0);
}
*/

	zonestat->init();
#pragma omp parallel for schedule(runtime) private(i,j,wert,wert2,k,mi,gradient,szenkorr)
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		double *entf;
		entf=new double[anzahl[0]+1];
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			//      if (zn[i][j] != nodata)
			if (((CoordX[i][j]>=xmin) && (CoordX[i][j] <= xmax) && (CoordY[i][j] >= ymin) && (CoordY[i][j] <= ymax))&&(zn[i][j] != nodata)&&(mpi_rank==MPIREG[i][j]) && !(doSpinUp && HRUMaster[i][j]==nodata))
			{
				if ((methode==2)||(methode==10)) { /* hoehenabhaengige Regression, hat als einzige Methode nicht die Abfrage auf gleichewerte und kommt deshalbals erste */
					if (WEIGHT)	{
						gradient = a[rarray[(int)(REGION[1][i][j])]][1];
					} else {
						gradient = a[1][1];
					}

					if (gradient >= -fehlwert) {/* interpolate only, if the gradient is valid in the actual region (else: take old values) */
						result[i][j] = getregress(i, j,CoordX[i][j], CoordY[i][j], dhm[i][j]);
					} else {
						NewGrid = false;
					}
				}
				else
				{
					if (!gleichewerte[0])
					{
						if (methode==1) { /* Inverse distance Interpolation */
							result[i][j] = aniso < 1.0 ? getidw_anisotropisch(0,CoordXR[i][j],CoordYR[i][j],entf) : getidw(0,CoordXR[i][j],CoordYR[i][j],entf);
							if ((warning_level) && (result[i][j] < -fehlwert))
							{
#pragma omp critical
								(*WSstream) << "WARNING:  at x=" << CoordXR[i][j] << ", y=" << CoordYR[i][j] << " " << _suchradius_ << "[end-entry]";
							}
						} else if ((methode==3)||(methode==11)) { /* Kombination aus IDW und hoehenabhaengiger Regression */
							/* wenn in Methode 3 gleichewerte true ist, dann heisst das, dass auch die Gradienten 0 sind, also nur die Ergebnisse von IDW zaehlen */
							if (WEIGHT) gradient = a[rarray[(int)(REGION[1][i][j])]][1];
							else        gradient = a[1][1];
							if (gradient >= -fehlwert) /* interpolate only, if the gradient is valid in the actual region (else: take old values) */
							{
								wert   =  getregress(i,j,CoordX[i][j],CoordY[i][j],dhm[i][j]);
								wert2  =  aniso < 1.0 ? getidw_anisotropisch(0,CoordXR[i][j],CoordYR[i][j],entf) : getidw(0,CoordXR[i][j],CoordYR[i][j],entf);
								if ((warning_level) && (wert2 < -fehlwert))
								{
#pragma omp critical
									(*WSstream) << "WARNING: at x=" << CoordXR[i][j] << ", y=" << CoordYR[i][j] << " " << _suchradius_ << "[end-entry]";
								}
								wert  *=  1.0-danteil;
								wert2 *=  danteil;
								result[i][j]    = wert + wert2;
							}
							else
							{
								result[i][j] = aniso < 1.0 ? getidw_anisotropisch(0,CoordXR[i][j],CoordYR[i][j],entf) : getidw(0,CoordXR[i][j],CoordYR[i][j],entf);
								if ((warning_level) && (result[i][j] < -fehlwert))
								{
#pragma omp critical
									(*WSstream) << "WARNING: at x=" << CoordXR[i][j] << ", y=" << CoordYR[i][j] << " " << _suchradius_ << "[end-entry]";
								}
							}
						} else if (methode==4) {/* Thiessen-Polygone */
							result[i][j] = getthiessen(0,CoordX[i][j],CoordY[i][j],entf);
						} else if (methode==5)	{ /* methode == 5 bilineare Interpolation */
							result[i][j] = getbilinear(0,i,j,CoordXR[i][j],CoordYR[i][j],entf);
						} else if (methode==6)	{ /* methode == 6 kombinierte bilineare Interpolation fuer Residuen und Gradienten */
							result[i][j] = getbilinear(0,i,j,CoordXR[i][j],CoordYR[i][j],entf)+
							dhm[i][j]*getbilinear(1,i,j,CoordXR[i][j],CoordYR[i][j],entf);
						} else if (methode==7)	{
							if (result[i][j]==nodata)
							{
								result[i][j] = 
								aniso < 1.0 ? 
								getidw_anisotropisch(0,CoordXR[i][j],CoordYR[i][j],entf) : 
								getidw(0,CoordXR[i][j],CoordYR[i][j],entf);
							}
						} else if (methode==8)	{
							if ((result[i][j]==nodata)||(result2[i][j]==nodata))
							{
								result[i][j] = 
								aniso < 1.0 ? 
								getidw_anisotropisch(0,CoordXR[i][j],CoordYR[i][j],entf) : 
								getidw(0,CoordXR[i][j],CoordYR[i][j],entf);

							}
							else
								result[i][j] = result[i][j]+dhm[i][j]*result2[i][j];
						} else if (methode==12) {/* Thiessen-Polygone with applied elevation lapse rate */
							if (upperlapserate>-9998)
							{
								result[i][j] = getThiessenWithLapseRate(CoordX[i][j],CoordY[i][j],dhm[i][j],lapserate,upperinversionfix,upperlapserate);
							}
							else
							{
								result[i][j] = getThiessenWithLapseRate(CoordX[i][j],CoordY[i][j],dhm[i][j],lapserate);
							}
						}
					}
					else  /* wenn alle Werte Fehlwerte sind, ist gleichewerte==1, dann soll keine Interpolation stattfinden,
				es gelten die alten Daten, ansonsten, wenn alle Werte der Eingabedatei (bis auf die Stationen mit Fehlwerten)
				gleich sind (für alle Stationen), dann wird jetzt gesucht, welche Station keinen Fehlwert hat (und
				angenommen, dass dann alle anderen Stationen mit gueltigen Werten auch diesen Wert haben, denn gleichewerte
				wurde beim Einlesen gesetzt als Ergebnis der Pruefung auf Gleichheit aller gueltigen Werte */
					{
						for (mi=1;mi<=(anzahl[0]);mi++) /* alle Stationen durchsuchen */
						{
							if (value[0][mi]>=-fehlwert)  /* erste gueltige Station gefunden */
							{
								value[0][1]=value[0][mi]; /* Wert uebernehmen */
								break;                    /* Suche kann abgebrochen werden, da davon ausgegangen wird, dass ale anderen STationen den gleichen Wert haben, so sie keinen Fehlwert haben */
							}
						}
						if (value[0][1] >= -fehlwert) /* alle anderen Werte koennen gleich dem ersten Wert gesetzt werden */
						{
							if (methode==12) /* Thiessen-Polygone with applied elevation lapse rate */
							{
								if (upperlapserate>-9998)
								{
									result[i][j] = getThiessenWithLapseRate(CoordX[i][j],CoordY[i][j],dhm[i][j],lapserate,upperinversionfix,upperlapserate);
								}
								else
								{
									result[i][j] = getThiessenWithLapseRate(CoordX[i][j],CoordY[i][j],dhm[i][j],lapserate);
								}
							}
							else
							{
								result[i][j] = value[0][1];
							}
						}
					}
				}
				// correction and scenario application can only be done when a new grid was created. I.e. if there where only nodata values in the input or no grid was read in, the 
				// old values should still be valid, including the already applied correction and scenario adjustion
				if (NewGrid)
				{
					result[i][j] *= kor;
					if (result[i][j] < grenzwertu) {
						result[i][j] = fixu;
					} 
					else if (result[i][j] > grenzwerto)  {
							result[i][j] = fixo; 
					}
					if (szenuse==1)
					{
						szenkorr= getszenkorr(i,j,CoordXR[i][j],CoordYR[i][j]);
						if ((szenmeth==1)||(szenmeth==5))
						result[i][j] += szenkorr; 
						else 
						result[i][j] *= szenkorr;
					}
					if (scalinggrid)
					{
						result[i][j] *= scalinggrid[i][j];
					}
				}
				//#pragma omp critical
				//		zonestat->wertzufuegen(i,j,zarray[(int)zn[i][j]],result[i][j]);
			}
			else
			{
				result[i][j] = nodata;
			}
		}
		delete[] entf;
	} // ende parallel-section
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if ((zn[i][j] != nodata)&&(mpi_rank==MPIREG[i][j]))
			{
				if (doSpinUp!=0 && HRUMaster[i][j]==nodata)
				{
					int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
					int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
					result[i][j] = result[imaster][jmaster];
				}
				zonestat->wertzufuegen(i,j,zarray[(int)zn[i][j]],result[i][j]);
			}
		}
	}
	zonestat->calc_and_done(jahr,monat,tag,stunde, minute);
	// will be written in an extra call now because 
	// the start and end of the vegetation period will be calculated after interpolation
	// (temperature data may be beeded for that calculation)
	// writegridtofile(resultnr); 
}

void meteointerpol::WriteResultGrid()
{
	writegridtofile(resultnr);
}

double meteointerpol::getsinglepoint(double xc,double yr,double zh,double *entf)
/* Diese Methode liefert das Interpolationsergebnis nur fuer eine einzige
Stelle mit den Koordinaten x=xc, y=yr, z=zh zurueck. Das wird bei der Korrektur der
Niederschlagsstationswerte genutzt */
{
	double xw,temp,entfernung;
	double xcrot,yrrot;
	int i,j,k;
	double  dx,dy;
	double  szenkorr;  /* der fuer einen Punkt zeitlich und raeumlich aus den Szenariowerten interpolierte Wert */

	/* Zeitliche Interpolation der Szenariodaten */
	if (szenuse==1)
	{
		if (tag<=15)
		{
			j=((int)monat+10)%12+3;
			k=((int)monat+11)%12+3;
			for (i=1;i<=szennum;i++)
			{
				szenario[i][15]=(double)(mon[j-2]-16+tag)/mon[j-2]*szenario[i][k] + (double)(16-tag)/mon[j-2]*szenario[i][j];
			}
		}
		else
		{
			j=((int)monat+11)%12+3;
			k=((int)monat+12)%12+3;
			for (i=1;i<=szennum;i++)
			{
				szenario[i][15]=(double)(tag-15)/mon[j-2]*szenario[i][k] + (double)(mon[j-2]-tag+15)/mon[j-2]*szenario[i][j];
			}
		}
	}

	if ((aniso<1.0)&&(xc!=nodata)&&(yr!=nodata))
	{
		dx = xc-xll;
		dy = yr-yll;
		entfernung = sqrt(dx*dx+dy*dy);
		/* aktueller Winkel des Punktes bezogen auf xll,yll */
		if (dx==0.0)
		alpha = dy > .0 ? Pi/2.0 : -Pi/2.0;
		else alpha = atan2(dy,dx);
		/* Winkel nach rechts drehen um Aniotropiewinkel */
		beta = alpha - anisoneig;
		/* neue Koordinaten */
		xcrot = xll+cos(beta)*entfernung;
		yrrot = yll+sin(beta)*entfernung;
	}
	else
	{
		xcrot=xc;
		yrrot=yr;
	}

	/* grid column */
	i=nrows-(int)((yr-yll)/cellsize);
	if  (i<1)  i=1; else if  (i>nrows)  i=nrows;
	/* grid row */
	j=(int)((xc-xll)/cellsize)+1;

	if  (j<1)  j=1; else if  (j>ncols)  j=ncols;

	if (methode==1) {
		xw = aniso < 1.0 ? getidw_anisotropisch(0,xcrot,yrrot,entf) : getidw(0,xc,yr,entf);
	} else if ((methode==2)||(methode==10)) {
		xw = getregress(i,j,xc,yr,zh);
	} else	if ((methode==3)||(methode==11)) {
		xw = aniso < 1.0 ? getidw_anisotropisch(0,xcrot,yrrot,entf) : getidw(0,xc,yr,entf); /* Ergebnis aus IDW */
		temp = getregress(i,j,xc,yr,zh); /* Ergebnis aus hoehenabhaengiger Regression */
		if ((xw >= -fehlwert)&&(temp >= -fehlwert)) 
		xw = danteil*xw + (1.0-danteil)*temp;
		else
		{
			/* das bedeutet, dass fabs(xw)>fehlwert ist, da ansonsten der else-Block nicht aufgerufen worden waere */
			if (temp >= -fehlwert) /* es gilt nur der Regressionswert */
			xw = temp;
		}
	} else  if (methode==4) {/* methode==4 */
		xw = getthiessen(0,xc,yr,entf);
	} else  if (methode==12) {/* methode==12 */
		if (upperlapserate>-9998)
		{
			xw = getThiessenWithLapseRate(xc,yr,zh,lapserate,upperinversionfix,upperlapserate);
		}
		else
		{
			xw = getThiessenWithLapseRate(xc,yr,zh,lapserate);
		}
	} else	if (methode==5) { /* bilinear */
		xw = getbilinear(0,i,j,xcrot,yrrot,entf);
	} else	if (methode==6)	{ /* methode==6 -> mehrfach bilninear */
		xw = getbilinear(0,i,j,xcrot,yrrot,entf)+dhm[i][j]*getbilinear(1,i,j,xcrot,yrrot,entf);
	} else if ((methode == 7)||(methode==8))	{ /* methode==7 oder 8 -> bicubic spline interpolation */
		xw = getthiessen(0,CoordXR[i][j],CoordYR[i][j],entf);
	} else { // methode == 9 read grids
		xw = result[i][j];
	}

	if (xw >= -fehlwert)
	{
		xw *= kor;
		if (xw < grenzwertu)  xw = fixu;
		else if (xw > grenzwerto)  xw = fixo;
		if (szenuse==1)
		{
			szenkorr= getszenkorr(i,j,xc,yr);
			if ((szenmeth==1)||(szenmeth==5))
			xw += szenkorr; 
			else 
			xw *= szenkorr;
		}
	}
	return xw;
}


meteointerpol::meteointerpol(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,
	int maxLastResults,char *meteonamennummer,char *steuername,int LogLevel,int progaus,
	double p_min_to_int_end,int p_timestep,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int p_resultnr,
	int *p_zarray,  statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax,  int p_zonmaxwert,center_of_mass_type *p_center_of_mass,
	int p_regiomax,int p_regiomaxwert,int *p_rarray,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll ,  double p_yll,
	double p_cellsize  ,  double p_nodata,
	long mina,long sta , long tga , long mta , long jra,
	long mine,long ste , long tge , long mte , long jre,
	LUSEtabletype *p_LUtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp,std::vector<HRUMapElement*> *p_vHRU):
	modelobject(Stream, steuername, LogLevel, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo, p_doSpinUp, p_vHRU),
statfilecnt(0),coordxrindex(0),coordyrindex(0),CoordXR(0),CoordYR(0),
stunde_a(0),tag_a(0),monat_a(0),minute_a(0),jahr_a(0),ascii(0),
regfilepos(0),anetzfilepos(0),regfilepos_start(0),anetzfilepos_start(0),regfilepos_preread(0),anetzfilepos_preread(0),
warning_level(0),bufferanetz(0),bufferregress(0),
bufmaxanetz(0),bufmaxregress(0),bufcntanetz(0),bufcntregress(0),
valcntanetz(0),valcntregress(0),regresssatz(0),numAdditionalColumns(0),
anetzstat(0),AdditionalColumns(0),dateinamekoef(0),ergebnisgridname(0),zonestatname(0),
zonestatcode(0),anetz(0),dateikoef(0),ergebnisgrid(0),zonestat(0),
methode(0),anzahl(0),useMinutes(NULL),fehlwerti(0),aniso(0),anisoneig(0),summedh(0),danteil(0),
grenzwertu(0),grenzwerto(0),fixu(0),fixo(0),maxentf(0),maxentfinv(0),wichtung(0),
af(0),bf1(0),bf2(0),igu(0),igo(0),a(0),b(0),alpha(0),beta(0),wi(0),STAT(0),
mindist(0),direction(0),slopeclass(0),dhmin(0),dvmin(0),AvgHorDist(0),AvgVertDist(0),
minstation(0),Spline(0),station(0),StationGrid(0),StationIndexGrid(0),
SecondDerivative(0),result2(0),CoordXStatR(0),CoordYStatR(0),StatMeshNumCols(0),StatMeshNumRows(0),
dhm(0),zn(0),xrot(0),yrot(0),value(0),hoehe(0),xkoord(0),ykoord(0),xkoord_rot(0),ykoord_rot(0),
gleichewerte(0),ende(0),js(0),ms(0),ts(0),ss(0),mins(0),meteosatz(0),t_meteosatz(0),kor(0),szenario(0),
SCEN(0),scengridname(0),scengrididentifier(0),scenmonth(0),scenindex(0),szenmeth(0),
szennum(0),szenuse(0),regio1index(0),regio2index(0),regio3index(0),weight1index(0),weight2index(0),
weight3index(0),REGION(0),WEIGHT(0),dxmin(0),dymin(0),dist(0),delta(0),sumdhmin(0),sumdvmin(0),
anzdhmin(0), anzdvmin(0), xcrot(0), yrrot(0), csroty(0), csrotx(0), gridlist_actname(0), gridlist_lastvalidname(0),
lowerinversionfix(0), upperinversionfix(0), tolerance(0), clusterlimit(0), overlap(0),scalinggrid(0),
xmin(p_xll),ymin(p_yll),xmax(p_xll+ p_cellsize *(p_ncols+1)),ymax(p_yll+ p_cellsize *(p_nrows+1))
{
	try
	{
		int i,j,s;

		/* Initialisieren der Adreseen und Daten mit den aus der
	allgemeinen Initialisierungsroutine uebergebenen Daten */
		retvalue=0;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif
		sprintf(instanzname,"%s",meteonamennummer);
		//    strcpy(instanzname,bezeichnung);
		(*WSstream) << "DEBUG: initializing interpolator for "<<instanzname<<"[end-entry]";

		gridinfo   = p_gridinfo;
		resultnr   = p_resultnr;
		nachricht  = progaus;
		min_to_int_end = p_min_to_int_end;
		this->timestep = p_timestep;
		numsubb    = p_numsubb;
		listsubb   = p_listsubb;
		if (listsubb) listsubb[0]=numsubb; /* number of subbasins to be considered for statistic files is stored in the zero-index */
		dhm        = gridinfo[2].pointer;
		zn         = gridinfo[1].pointer;
		zarray     = p_zarray;
		zone       = p_zone;
		FirstLastCol = p_FirstLastCol;
		zonmax     = p_zonmax;
		zonmaxwert = p_zonmaxwert;
		center_of_mass = p_center_of_mass;
		regiomax   = p_regiomax;
		regiomaxwert = p_regiomaxwert;
		rarray     = p_rarray;
		ncols      = p_ncols;
		nrows      = p_nrows;
		mpistartrow = p_mpistartrow;
		mpiendrow = p_mpiendrow;
		xll        = p_xll;
		yll        = p_yll;
		cellsize   = p_cellsize;
		nodata     = p_nodata;
		minutea    = mina;
		stundea    = sta;
		taga       = tga;
		monata     = mta;
		jahra      = jra;
		minutee    = mine;
		stundee    = ste;
		tage       = tge;
		monate     = mte;
		jahre      = jre;
		LUtable		= p_LUtable,
		MultiLayerLUTable = p_MultiLayerLUTable;

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

		if (doSpinUp!=0 && HRU==NULL)
		{
			HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
			HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
		}

		coordxrindex=coordyrindex=-1;
		endtime = jf*jahre+mf*monate+df*tage+hf*stundee+minf*minutee; /* for comparing the time from input file with the end-date */

		if (ist_vorhanden(_RelCellArea_))
		{
			RelCellArea=gridinfo[getindex(_RelCellArea_)].pointer;
			CoordX=gridinfo[getindex(_CoordX_)].pointer;
			CoordY=gridinfo[getindex(_CoordY_)].pointer;
			/* if the coordinate system will not be rotated, then the original 
		coordinates can be used (here done by using CoordXR as alias for CoordX */
			CoordXR=CoordX;
			CoordYR=CoordY;
		}
		else
		{
			(*WSstream) << "FATAL: \n-----------------------------------------------------------\n"
			"WaSiM-ETH no longer works without a RelCellArea-Grid.\n"
			"If this message appears, then something is wrong with the code,\n"
			"so please contact the author.\n"
			"-----------------------------------------------------------\n\n" <<"[end-entry]";
			throw(74);
		}

		js = newstring(50,1);
		ms = newstring(50,1);
		ts = newstring(50,1);
		ss = newstring(50, 1);
		mins = newstring(50, 1);

		/* reservieren von Speicher fuer Parameter-Arrays (11.6.1999) */
		hoehe=(double**)calloc(2,sizeof(double*));
		value=(double**)calloc(2,sizeof(double*));
		xkoord=(double**)calloc(2,sizeof(double*));
		ykoord=(double**)calloc(2,sizeof(double*));
		entf=(double**)calloc(2,sizeof(double*));
		tag_a=(long*)calloc(2,sizeof(long));
		minute_a = (long*)calloc(2, sizeof(long));
		stunde_a = (long*)calloc(2, sizeof(long));
		monat_a=(long*)calloc(2,sizeof(long));
		jahr_a=(long*)calloc(2,sizeof(long));
		anetzfilepos=(long*)calloc(2,sizeof(long));
		anetzfilepos_start=(long*)calloc(2,sizeof(long));
		anetzfilepos_preread=(long*)calloc(2,sizeof(long));
		anetz=(FILE**)calloc(2,sizeof(FILE*));
		anetzstat=(char**)calloc(2,sizeof(char*));
		AdditionalColumns=(char**)calloc(2,sizeof(char*));
		numAdditionalColumns=(int*)calloc(2,sizeof(int));
		anzahl=(int*)calloc(2,sizeof(int));
		useMinutes = (bool*)calloc(2, sizeof(bool));
		bufferanetz=(buffertype**)calloc(2,sizeof(buffertype*));
		bufcntanetz=(long*)calloc(2,sizeof(long));
		bufmaxanetz=(long*)calloc(2,sizeof(long));
		valcntanetz=(int*)calloc(2,sizeof(int));
		
		regio1index=regio2index=regio3index=weight1index=weight2index=weight3index=0;

		// supressing the use of a region grid by simply not allowing the if branch to be processed 
		// (as if no rarray was created and no regions were used)
		if ((1==2)&&(rarray)) /* es wurde ein Regionengrid eingelesen und ausgewertet, es gibt also mindestens eine Region */
		{
			WEIGHT=(double***)calloc(4,sizeof(double**));
			REGION=(double***)calloc(4,sizeof(double**));
			REGION[1]=REGION[2]=REGION[3]=WEIGHT[1]=WEIGHT[2]=WEIGHT[3]=NULL;
			regio1index= ist_vorhanden(_regression_regions_1_);
			weight1index=ist_vorhanden(_regression_weights_1_);
			REGION[1]=gridinfo[regio1index].pointer;
			WEIGHT[1]=gridinfo[weight1index].pointer;
			if (regiomax>1)
			{
				regio2index= ist_vorhanden(_regression_regions_2_);
				weight2index=ist_vorhanden(_regression_weights_2_);
				REGION[2]=gridinfo[regio2index].pointer;
				WEIGHT[2]=gridinfo[weight2index].pointer;
			}
			if (regiomax>2)
			{
				regio3index= ist_vorhanden(_regression_regions_3_);
				weight3index=ist_vorhanden(_regression_weights_3_);
				REGION[3]=gridinfo[regio3index].pointer;
				WEIGHT[3]=gridinfo[weight3index].pointer;
			}
		}
		else
		{
			regiomax=1; /* ansonsten wird regiomax, dass bisher 0 war, zu eins gesetzt (als maximaler Zaehler fuer die regressdateien) */
			REGION=WEIGHT=NULL;
		}

		dateinamekoef = (char**)calloc(regiomax+1,sizeof(char*));
		for (i=1;i<=regiomax;i++)
		dateinamekoef[i] = newstring(_stringsize_+1,1);

		gleichewerte=(unsigned short*)calloc(regiomax>2?regiomax+1:2,sizeof(unsigned short));
		bufferregress=(buffertype**)calloc(regiomax+1,sizeof(buffertype*));
		ascii        = (int*)calloc(regiomax+1,sizeof(int));
		regfilepos   = (long*)calloc(regiomax+1,sizeof(long));
		regfilepos_start   = (long*)calloc(regiomax+1,sizeof(long));
		regfilepos_preread = (long*)calloc(regiomax+1,sizeof(long));
		regresssatz  = (char**)calloc(regiomax+1,sizeof(char*));
		bufmaxregress=(long*)calloc(regiomax+1,sizeof(long));
		bufcntregress=(long*)calloc(regiomax+1,sizeof(long));
		valcntregress= (int*)calloc(regiomax+1,sizeof(int));
		dateikoef  = (FILE**)calloc(regiomax+1,sizeof(FILE*));
		af  = (double*)calloc(regiomax+1,sizeof(double));
		bf1 = (double*)calloc(regiomax+1,sizeof(double));
		bf2 = (double*)calloc(regiomax+1,sizeof(double));
		igu = (double*)calloc(regiomax+1,sizeof(double));
		igo = (double*)calloc(regiomax+1,sizeof(double));
		a   = (double**)calloc(regiomax+1,sizeof(double*));
		
		for (i=0;i<=regiomax;i++)
        {
			a[i]=(double*)calloc(5,sizeof(double));
        }
		
		b   = (double**)calloc(regiomax+1,sizeof(double*));
		
		for (i=0;i<=regiomax;i++)
        {
			b[i]=(double*)calloc(5,sizeof(double));
        }

		ergebnisgridname = newstring(_stringsize_+1,1);
		zonestatname     = newstring(_stringsize_+1,1);

		meteosatz=newstring(8192,1);
		t_meteosatz=meteosatz;

		parametereingabe();

		// when reading in grids, rank 0 should allocate the complete grid in order to read ascii grids and to fill up missing values before scattering the files to the other ranks
		if ((mpi_rank == 0) && (methode==9))
		{
			/* belegen von Speicher fuer name und inhalt des Ergebnisgrids auf dem Speicher */
			result = matrix(1, nrows, 1, ncols, nrows);
			initialize_grid(resultnr, ergebnisgridname, result, NULL, 0, instanzname, writegrid, 1, nrows);
		}
		else
		{
			/* belegen von Speicher fuer name und inhalt des Ergebnisgrids auf dem Speicher */
			result = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(resultnr, ergebnisgridname, result, NULL, 0, instanzname, writegrid, mpistartrow, mpiendrow);
		}

		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (j=1;j<=ncols;j++) 
			{
				result[i][j]=nodata;
			}
		}

		/* if scenario grids should be used the first grid have to be read in now */
		if ((szenmeth==5)||(szenmeth==6))
		{
			scenindex=resultnr+1;
			scenmonth=monata;
			// allocating memeory when reading scenario grids for the first time
			readgridfromfile(scengridname[scenmonth], &SCEN, 1, _AllocateGridMemoryOn_, _ShowMessagesOff_, 0, _LimitAutoFillOff_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			sprintf(scengrididentifier,"scenariogrid_%s",instanzname);
			initialize_grid(scenindex,scengridname[scenmonth],SCEN,NULL,0,scengrididentifier,0,mpistartrow,mpiendrow);
#ifdef _MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
		else
		{
			SCEN=NULL;
		}

		/* Datenspeicher mit dem korrekten Fehlwert belegen */
		for (s=0;s<statfilecnt;s++) 
		{
			minute_a[s] = -fehlwerti - 1;
			stunde_a[s] = -fehlwerti - 1;
			tag_a[s]    = -fehlwerti-1;
			monat_a[s]  = -fehlwerti-1;
			jahr_a[s]   = -fehlwerti-1;
		}
		stunde=-fehlwerti-1;
		tag=-fehlwerti-1;
		monat=-fehlwerti-1;
		jahr=-fehlwerti-1;

		if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12))
		{
			for (s=0;s<statfilecnt;s++) 
			{
				stationswerte(anetz[s],anetzstat[s],numAdditionalColumns[s],&(useMinutes[s]),&(anzahl[s]),&(hoehe[s]),&(xkoord[s]),&(ykoord[s]),&(entf[s]),&(value[s]));
				if (s>0)
				{
					if (anzahl[s]!=anzahl[0])
					{
						(*WSstream) << "FATAL: \n------- Error 73 in Interpolation input data ---------\n"
						<< "meteorological station files "<<anetzstat[s]
						<<" \nand "<<anetzstat[0]
						<< " contain a different number of stations \n("<<anzahl[s]
						<<" and "<<anzahl[0]<<", resp.)\n" 
						<< "these data has to be identical in both files for interpolation method 6\n"
						<< "------- Program terminated due to input data error ---------\n\n" << "[end-entry]";
						throw(73);
					}
					else
					{
						for (int b=1;b<=anzahl[0];b++)
						{
							if ((hoehe[s][b]!=hoehe[0][b])||
									(xkoord[s][b]!=xkoord[0][b])||
									(ykoord[s][b]!=ykoord[0][b]))
							{
								(*WSstream) << "FATAL: \n------- Error 73 in Interpolation input data ---------\n"
								<< "Altitude and/or coordinates of station "<<b<<" are different in\n"
								"meteorological station files "<<anetzstat[s]
								<<"\nH: "<<hoehe[s][b]<<",\tX: "<<xkoord[s][b]<<"\tY: "<<ykoord[s][b]
								<<"\nand\n"<<anetzstat[0]<<"\nH: "<<hoehe[0][b]
								<<",\tX: "<<xkoord[0][b]<<"\tY: "<<ykoord[0][b]<<"\n"
								<< "These data has to be identical in both files for interpolation method 6\n"
								<< "------- Program terminated due to input data error ---------\n\n" << "[end-entry]";
								throw(73);
							}
						}
					}
				}
				//    stationswerte(); /* Einlesen der Stationswerte Hoehe, x- und y-Koordinaten */
				finddateanetz(anetz[s],anetzstat[s],numAdditionalColumns[s],useMinutes[s],anzahl[s],value[s],&(anetzfilepos[s]),&(anetzfilepos_preread[s]),
								&(gleichewerte[s]), &(minute_a[s]), &(stunde_a[s]), &(tag_a[s]),&(monat_a[s]),&(jahr_a[s]));
				anetzfilepos_start[s]=anetzfilepos_preread[s]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
				//    finddateanetz(); /* Suchen des Anfangsdatums in der EIngabedatei */
				/* allocate buffer memory for input files */
				bufmaxanetz[s]      = (unsigned long)(buffersize/(5*sizeof(long)+(anzahl[s]+1)*sizeof(double)));
				if (bufmaxanetz[s]<1) bufmaxanetz[s]=1;
				bufcntanetz[s]      = bufmaxanetz[s]+1; /* initialising the anetz-buffer index counter with max+1 (next check must result in call to readbufferanetz()) */
				bufferanetz[s]      = (buffertype*)calloc((size_t)bufmaxanetz[s]+1,sizeof(buffertype));
				for (i=1;i<=bufmaxanetz[s];i++) 
				bufferanetz[s][i].wert=(double*)calloc(anzahl[s]+1,sizeof(double)); /* allocating memory for each station */
			}
			s=0;
			
			mindist=(double*)calloc(anzahl[s]+1,sizeof(double));
			dhmin=(double*)calloc(anzahl[s]+1,sizeof(double));
			dvmin=(double*)calloc(anzahl[s]+1,sizeof(double));
			direction=(double*)calloc(anzahl[s]+1,sizeof(double));
			slopeclass=(double*)calloc(37,sizeof(double)); /* 10 degree intervals */
			minstation=(int*)calloc(anzahl[s]+1,sizeof(int)); /* 10 degree intervals */
			
			/* general slope of the "station" mesh, that is, the minimal slope to 
		rotate the grid into a cartesian coordinate system with horizontal 
		and vertical axis' */

			GetBilinearSkewness(); 
			anisoneig *= Pi/180.0; /* converting the annisotropy slope into radiant */
			/* transform the grid cell coordinates into the rotated system */
			RotateCoordGrids();
			/* rotate the station coordinates into the rotated system */
			RotateStationCoordiantes();
			GetMeanStationDist(0);
			if ((methode==5)||(methode==6))
			{
				/* output for method 5 and 6: 4 Grids serving as look-up-grids for the 
		four stations used for bilinear interpolation for each grid cell,
		output for method 7: an array of "station"-columns which contain each 
		the indizes of all stations belonging to a column */
				MakeBilinearStationGrid(statfilecnt);
			}
			else
			if ((methode==7)||(methode==8))
			{ /* building the station-columns for spline interpolation */
				ConstructGridFromStationCoordinates(s);
			}


			if (mindist) free(mindist); mindist=NULL;
			if (dhmin) free(dhmin); dhmin=NULL;
			if (dvmin) free(dvmin); dvmin=NULL;
			if (direction) free(direction); direction=NULL;
			if (slopeclass) free(slopeclass); slopeclass=NULL;
			if (minstation) free(minstation); minstation=NULL;

			/* stunde, tag monat und Jahr sind Variablen in modelobject. Das Array fuer eine 
		ganze Reihe Stationsdateien heisst daher etwas anders. Geprüft wird jeweils beim Einlesen der 
		Daten, ob in allen Dateien die Werte gleich sind. */
			minute = minute_a[0];
			stunde = stunde_a[0];
			tag=tag_a[0];
			monat=monat_a[0];
			jahr=jahr_a[0];
		}
		if ((methode==2)||(methode==3))
		{
			for (j=1;j<=regiomax;j++)
			{
				if (timestep < 60)
				{
					(*WSstream) << "FATAL: interpolation method 2 or 3 with externally prepared regression files does not work for time steps less than one hour - please use interpolation methods 10 or 11 instead" << "[end-entry]";
					throw(1302);
				}
				finddateregress(j);
				regfilepos_start[j]=regfilepos_preread[j]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
				/* allocate buffer memory for input files */
				bufmaxregress[j]      = (unsigned long)(buffersize/(5*sizeof(long)+12*sizeof(double)));  /* regress-file format is fix: 4 long and 11 double */
				if (bufmaxregress[j]<1) bufmaxregress[j]=1;
				bufcntregress[j]      = bufmaxregress[j]+1; /* initialising the regress-buffer index counter with max+1 (next check must result in a call to readbufferregress()) */
				bufferregress[j]      = (buffertype*)calloc((size_t)bufmaxregress[j]+1,sizeof(buffertype));
				for (i=1;i<=bufmaxregress[j];i++) 
				bufferregress[j][i].wert=(double*)calloc(12,sizeof(double)); /* allocating memory for the regression-parameters */
			}
		}
		if (methode==9)
		{
			// read grids according to the names contained in a list file
			// this method works withput a buffer since disk access is required anyway
			gridlist_actname = new char[_stringsize_ + 1];
			sprintf(gridlist_actname, "");
			gridlist_lastvalidname = new char[_stringsize_ + 1];
			sprintf(gridlist_lastvalidname, "");
			stationswerte(anetz[0],anetzstat[0],numAdditionalColumns[0],&(useMinutes[0]),&(anzahl[0]),&(hoehe[0]),&(xkoord[0]),&(ykoord[0]),&(entf[0]),&(value[0]));
			finddategridlist(anetz[0],anetzstat[0],numAdditionalColumns[0],useMinutes[0],gridlist_actname,
							&(anetzfilepos[0]),&(anetzfilepos_preread[0]), &(minute_a[0]), &(stunde_a[0]), &(tag_a[0]),&(monat_a[0]),&(jahr_a[0]));
			anetzfilepos_start[0]=anetzfilepos_preread[0]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
			modelobject::SetExternalGridProperties(xkoord[0][1],ykoord[0][1],hoehe[0][1],zn);
			gleichewerte[0]=false;
			minute = minute_a[0];
			stunde = stunde_a[0];
			tag=tag_a[0];
			monat=monat_a[0];
			jahr=jahr_a[0];

		}
		ende = false; /* wird in nextrecord auf false gesetz, wenn die Datei nicht zu Ende ist */
		if (methode == 0)  ende = true; else ende = false;
		if (methode == 1) {
			sprintf(meteosatz, "%s interpolated with idw from %s, zones from %s, weighting = %4.1f,  maxentf = %8.0f, control file %s, additional columns: %d",
				instanzname, (anetzstat[0]), gridinfo[1].name, wichtung, wichtung == 2 ? sqrt(maxentf) : maxentf, steuername, numAdditionalColumns[0]);
		}
		else if (methode == 2) {
			sprintf(meteosatz, "%s interpolated with regression from ", instanzname);
			for (i = 1; i <= regiomax; i++)
				sprintf(meteosatz + strlen(meteosatz), " %s,", dateinamekoef[i]);
			sprintf(meteosatz + strlen(meteosatz), " zones from %s, control file: %s", gridinfo[1].name, steuername);
		}
		else if (methode == 3) {
			sprintf(meteosatz, "%s interpolated with idw+regression from %s, zones from %s,weighting = %4.1f, fraction distance: = %4.1f max distance = %8.0f, control file %s, additional columns: %d",
				instanzname, (anetzstat[0]), gridinfo[1].name, wichtung, danteil, wichtung == 2 ? sqrt(maxentf) : maxentf, steuername, numAdditionalColumns[0]);
		}
		else if (methode == 4) {
			sprintf(meteosatz, "%s interpolated with Thiessen from %s, zones from %s, control file %s, additional columns: %d",
				instanzname, (anetzstat[0]), gridinfo[1].name, steuername, numAdditionalColumns[0]);
		}
		else if (methode == 5) {
			sprintf(meteosatz, "%s interpolated with bilinear interpolation from %s, zones from %s, control file %s, additional columns: %d",
				instanzname, (anetzstat[0]), gridinfo[1].name, steuername, numAdditionalColumns[0]);
		}
		else if (methode == 6) {
			sprintf(meteosatz, "%s interpolated with bilinear interpolation of gradients (%s, add. col.: %d) and residuals (%s, add. col.: %d), linearly combined, zones from %s, control file %s",
				instanzname, anetzstat[0], numAdditionalColumns[0], anetzstat[1], numAdditionalColumns[1], gridinfo[1].name, steuername);
		}
		else if (methode == 7) {
			sprintf(meteosatz, "%s interpolated with bicubic spline interpolation (%s, add. col.: %d), zones from %s, control file %s",
				instanzname, anetzstat[0], numAdditionalColumns[0], gridinfo[1].name, steuername);
		}
		else if (methode == 8) {
			sprintf(meteosatz, "%s interpolated with bicubic spline interpolation of gradients (%s, add. col.: %d) and residuals (%s, add. col.: %d), linearly combined, zones from %s, control file %s",
				instanzname, anetzstat[0], numAdditionalColumns[0], anetzstat[1], numAdditionalColumns[1], gridinfo[1].name, steuername);
		}
		else if (methode == 9) {
			sprintf(meteosatz, "%s read in from grids containing already interpolated data (%s, add. col.: %d) zones from %s, control file %s",
				instanzname, anetzstat[0], numAdditionalColumns[0], gridinfo[1].name, steuername);
		}
		else if (methode == 10) {
			sprintf(meteosatz, "%s interpolated with elevation dependent regression for values taken from %s, control file %s", instanzname,
				anetzstat[0], steuername);
		}
		else if (methode == 11) {
			sprintf(meteosatz, "%s interpolated with elevation dependent regression + idw, all values taken from %s, control file %s", instanzname,
				anetzstat[0], steuername);
		}
		else if (methode == 12) {
			sprintf(meteosatz, "%s interpolated with Thiessen plus elevation lapserate, all values taken from %s, control file %s", instanzname,
				anetzstat[0], steuername);
		}
		else {
			(*WSstream) << "FATAL: error! interpolation method has to range from 1..12 (here: " << methode << ")\n" << "[end-entry]";
			throw(53);
		}

		zonestat = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
			&zone, (int)zonmax, (int)zonestatcode, (int)3, (char*)zonestatname, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, useMinutes[0]);
		if (retvalue != 0) throw retvalue;
		zonestat->schreibzonenwerte(meteosatz, "", LastResults, maxLastResults);

		meteosatz = t_meteosatz;
		if (meteosatz) free(meteosatz);
		meteosatz = newstring(_meteosatzlaenge_ + 1, 1);
		t_meteosatz = meteosatz;

		// looking for a scalinggrid
		char *scalinggridname = new char[1024];
		sprintf(scalinggridname, "%s_%s", _ScalingGrid_, this->instanzname);
		if (ist_vorhanden(scalinggridname))
		{
			scalinggrid = gridinfo[getindex(scalinggridname)].pointer;
		}
		delete[] scalinggridname;
	}
	catch (int ex) 
	// all exceptions must be caught within the constructor, because throwing an exception 
	// outside off this constructor would cause the object to remain uninitialized
	// resulting in the inpossibility of destructing (i.e., the memory clearing) of the 
	// instance and hence in memory leaks
	{
		retvalue=ex;
		throw ex;
	}
	catch (...)
	// also unknown exceptions have to be caught...
	{
		retvalue=UNKNOWN_EXCEPTION;
		throw UNKNOWN_EXCEPTION;
	}
}

// used for resetting the input file date when model is in SpinUp mode and multiple cycles through the input data must be done
void meteointerpol::resetInputFiles()
{
	this->ende=0;
	if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12))
	{
		for (int s=0;s<statfilecnt;s++) 
		{
			anetzfilepos[s]=anetzfilepos_start[s]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
			bufcntanetz[s] = bufmaxanetz[s]+1; /* initialising the anetz-buffer index counter with max+1 (next check must result in call to readbufferanetz()) */
			jahr_a[s]   = jahra;
			monat_a[s]  = monata;
			tag_a[s]    = taga;
			stunde_a[s] = stundea;
			minute_a[s] = minutea;
		}
	}
	if ((methode==2)||(methode==3))
	{
		for (int j=1;j<=regiomax;j++)
		{
			regfilepos[j]=regfilepos_start[j]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
			bufcntregress[j]      = bufmaxregress[j]+1; /* initialising the regress-buffer index counter with max+1 (next check must result in a call to readbufferregress()) */
			jahr_a[j]   = jahra;
			monat_a[j]  = monata;
			tag_a[j]    = taga;
			stunde_a[j] = stundea;
			minute_a[j] = minutea;
		}
	}
	if (methode==9)
	{
		anetzfilepos[0]=anetzfilepos_start[0]; // will be used in SpinUp mode to reset the file pointer when doing multiple cycles
		jahr_a[0]   = jahra;
		monat_a[0]  = monata;
		tag_a[0]    = taga;
		stunde_a[0] = stundea;
		minute_a[0] = minutea;
	}
}

//meteointerpol::meteointerpol(char *steuername):
//modelobject(steuername)
//{
//}


/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void meteointerpol::ClearMemory()
{
	int i,j,k,s;

	(*WSstream) << "DEBUG: deleting interpolator for "<<instanzname<<"[end-entry]";

	/* befreien des Speichers fuer das Ergebnisgrid muss innerhalb des Objektes erfolgen,
	da auch innerhalb des Objektes der Speicher reserviert wurde. Es kann aber
	auch ausserhalb erfolgen, da die Adresse von gridinfo als Parameter erst an das
	Object uebergeben wurde */
	if (szenuse==1) 
	{
		if ((szenmeth==5)||(szenmeth==6))
		{
			for (i=1;i<=12;i++)
			{
				if (scengridname[i]) free(scengridname[i]);
				scengridname[i]=NULL;
			}
			if (scengrididentifier) free(scengrididentifier);
			scengrididentifier=NULL;
			free_matrix(gridinfo[scenindex].pointer, mpistartrow, mpiendrow, 1, ncols, nrows);
			delete_grid(scenindex);
		}
		else
		free_matrix(szenario,1,szennum,1,15); szenario=NULL;
	}

	if (station) delete station; station=NULL;

	free_matrix(gridinfo[resultnr].pointer, gridinfo[resultnr].nrl, gridinfo[resultnr].nrh, 1, ncols, nrows);
	delete_grid(resultnr);

	if (coordxrindex>0)
	{ /* then, the Coordinate Grids are extra reserved grids and have to be freed */
		free_matrix(gridinfo[coordxrindex].pointer, gridinfo[coordxrindex].nrl, gridinfo[coordxrindex].nrh, 1, ncols, nrows);
		delete_grid(coordxrindex);
	}

	if (coordyrindex>0)
	{
		free_matrix(gridinfo[coordyrindex].pointer, gridinfo[coordyrindex].nrl, gridinfo[coordyrindex].nrh, 1, ncols, nrows);
		delete_grid(coordyrindex);
	}

	if (zonestat) 
	{
		zonestat->ClearMemory();
		delete zonestat;
		zonestat=NULL;
	}

	if ((methode==5)||(methode==6)) /* bilinear interpolation */
	{
		if (STAT) 
		{
			for (k=1;k<=4*(statfilecnt+1);k++)
			{
				if (STAT[k])
				{
					for (i=mpistartrow;i<=mpiendrow;i++)
					{
						if (STAT[k][i]) free(STAT[k][i]);
						STAT[k][i]=NULL;
					}
					free(STAT[k]);
					STAT[k]=NULL;
				}
			}
			free(STAT);
			STAT=NULL;
		}
	}

	if ((methode==7)||(methode==8))
	{
		free_matrix(StationGrid,1,StatMeshNumRows,1,StatMeshNumCols);
		free_matrix(StationIndexGrid,1,StatMeshNumRows,1,StatMeshNumCols);
		free_matrix(SecondDerivative,1,StatMeshNumRows,1,StatMeshNumCols);
		free_matrix(CoordXStatR,1,StatMeshNumRows,1,StatMeshNumCols);
		free_matrix(CoordYStatR,1,StatMeshNumRows,1,StatMeshNumCols);
		if (Spline) delete Spline;
		if (methode == 8)
		{
			free_matrix(result2, mpistartrow, mpiendrow, 1, ncols, nrows);
		}
	}



	if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12))
	{
		if (bufferanetz)
		{
			for (s=0;s<statfilecnt;s++)
			{
				if (bufferanetz[s]) 
				{
					for (i=1;i<=bufmaxanetz[s];i++)
					{
						if (bufferanetz[s][i].wert) 
						free(bufferanetz[s][i].wert);
						bufferanetz[s][i].wert=NULL;
					}
					free(bufferanetz[s]);
					bufferanetz[s]=NULL;
				}
			}
		}
	}

	if ((methode==2)||(methode==3))
	{
		if (bufferregress) 
		{
			for (j=1;j<=regiomax;j++)
			{
				if (ascii[j]) 
				{
					if (regresssatz[j]) 
					free(regresssatz[j]);
					regresssatz[j]=NULL;
				}
				if (bufferregress[j]) 
				{
					for (i=1;i<=bufmaxregress[j];i++)
					{
						if (bufferregress[j][i].wert) 
						free(bufferregress[j][i].wert);
						bufferregress[j][i].wert=NULL;
					}
					free(bufferregress[j]);
					bufferregress[j]=NULL;
				}
			}
			free(bufferregress);
			bufferregress=NULL;
		}
	}

	if ((methode==1)||(methode==3)||(methode==4)||(methode==5)||(methode==6)||(methode==7)||(methode==8)||(methode==10)||(methode==11)||(methode==12))
	{
		for (s=0;s<statfilecnt;s++)
		{
			if (xkoord_rot) 
			{
				if (xkoord_rot[s]&&(aniso<1.0)) 
				{
					free(xkoord_rot[s]);
					xkoord_rot[s]=NULL;
				}
			}
			if (ykoord_rot)
			{
				if (ykoord_rot[s]&&(aniso<1.0)) 
				{
					free(ykoord_rot[s]);
					ykoord_rot[s]=NULL;
				}
			}
			if (value)  if (value[s])     { free(value[s]); value[s]=NULL;   }
			if (entf)   if (entf[s])      { free(entf[s]);  entf[s]=NULL;    } 
			if (ykoord) if (ykoord[s])    { free(ykoord[s]); ykoord[s]=NULL; }
			if (xkoord) if (xkoord[s])    { free(xkoord[s]); xkoord[s]=NULL; }
			if (hoehe)  if (hoehe[s])     { free(hoehe[s]);  hoehe[s]=NULL;  }
		}
		if (xkoord_rot) free(xkoord_rot); xkoord_rot=NULL;
		if (ykoord_rot) free(ykoord_rot); ykoord_rot=NULL;
	}

	if (ascii) free(ascii); ascii=NULL;
	if (regresssatz) free(regresssatz); regresssatz=NULL;
	if (regfilepos) free(regfilepos); regfilepos=NULL;
	if (regfilepos_start) free(regfilepos_start); regfilepos_start = NULL; if (regfilepos_preread) free(regfilepos_preread); regfilepos_preread = NULL;
	if (bufmaxregress) free(bufmaxregress); bufmaxregress=NULL;
	if (bufcntregress) free(bufcntregress); bufcntregress=NULL;
	if (valcntregress) free(valcntregress); valcntregress=NULL;
	if (dateikoef) free(dateikoef); dateikoef=NULL;
	if (af) free(af); af=NULL;
	if (bf1) free(bf1); bf1=NULL;
	if (bf2) free(bf2); bf2=NULL;
	if (igu) free(igu); igu=NULL;
	if (igo) free(igo); igo=NULL;

	if (a) 
	{
		for (i=0;i<=regiomax;i++)
		{
			if (a[i]) 
			free(a[i]);
			a[i]=NULL;
		}
		free(a);
		a=NULL;
	}

	if (b) 
	{
		for (i=0;i<=regiomax;i++)
		{
			if (b[i]) 
			free(b[i]);
			b[i]=NULL;
		}
		free(b);
		b=NULL;
	}

	if (dateinamekoef    ) 
	{
		for (i=1;i<=regiomax;i++) 
		{
			if (dateinamekoef[i]) 
			free(dateinamekoef[i]);
			dateinamekoef[i]=NULL;
		}
		free(dateinamekoef    ); 
		dateinamekoef    =NULL;
	}

	if (zonestatname     ) free(zonestatname     ); zonestatname     =NULL;
	if (ergebnisgridname ) free(ergebnisgridname ); ergebnisgridname =NULL;

	if (AdditionalColumns) 
	{
		for (s=0;s<statfilecnt;s++)
		{
			if (AdditionalColumns[s]) 
			free(AdditionalColumns[s]);
			AdditionalColumns[s]=NULL;
		}
		free(AdditionalColumns); 
		AdditionalColumns=NULL;
	}

	if (anetzstat) 
	{
		for (s=0;s<statfilecnt;s++)
		{
			if (anetzstat[s]) 
			free(anetzstat[s]);
			anetzstat[s]=NULL;
		}
		free(anetzstat); 
		anetzstat=NULL;
	}

	if (hoehe)        free(hoehe);               hoehe = NULL;
	if (value)        free(value);               value = NULL;
	if (xkoord)       free(xkoord);             xkoord = NULL;
	if (ykoord)       free(ykoord);             ykoord = NULL;
	if (entf)         free(entf);                 entf = NULL;
	if (tag_a)        free(tag_a);               tag_a = NULL;
	if (minute_a)     free(minute_a);         minute_a = NULL;
	if (stunde_a)     free(stunde_a);         stunde_a = NULL;
	if (monat_a)      free(monat_a);           monat_a = NULL;
	if (jahr_a)       free(jahr_a);             jahr_a = NULL;
	if (anetzfilepos) free(anetzfilepos); anetzfilepos = NULL;
	if (anetzfilepos_start) free(anetzfilepos_start); anetzfilepos_start = NULL; 	if (anetzfilepos_preread) free(anetzfilepos_preread); anetzfilepos_preread = NULL;
	if (anetz)        free(anetz);               anetz = NULL;
	if (gleichewerte) free(gleichewerte); gleichewerte = NULL;
	if (bufferregress) free(bufferregress); bufferregress = NULL;
	if (numAdditionalColumns) free(numAdditionalColumns); numAdditionalColumns=NULL;
	if (anzahl)       free(anzahl);        anzahl      = NULL;
	if (useMinutes)	  free(useMinutes);	   useMinutes  = NULL;
	if (bufferanetz)  free(bufferanetz);   bufferanetz = NULL;
	if (bufcntanetz)  free(bufcntanetz);   bufcntanetz = NULL;
	if (bufmaxanetz)  free(bufmaxanetz);   bufmaxanetz = NULL;
	if (valcntanetz)  free(valcntanetz);   valcntanetz = NULL;

	if (mins) free(mins); mins = NULL;
	if (ss) free(ss); ss = NULL;
	if (ts) free(ts); ts=NULL;
	if (ms) free(ms); ms=NULL;
	if (js) free(js); js=NULL;
	//  meteosatz=t_meteosatz; /* poiting to the correct memory array for meteosatz */
	if (meteosatz) free(meteosatz); meteosatz=NULL;
	
	if (gridlist_actname) delete[] gridlist_actname; gridlist_actname = NULL;
	if (gridlist_lastvalidname) delete[] gridlist_lastvalidname; gridlist_lastvalidname = NULL;

	modelobject::ClearMemory();
}


void meteointerpol::getregresscoeffs(int nr,double *p_igu,double *p_igo,double *p_a,double *p_b)
{
	if ((methode!=2)&&(methode!=3)&&(methode!=10)&&(methode!=11))
	{
		*p_igu = *p_igo = p_a[0] = p_a[1] = p_a[2] = p_a[3] = p_b[0] = p_b[1] = p_b[2] = p_b[3] = -fehlwert -1.0;
	}
	else
	{
		*p_igu = igu[nr];
		*p_igo = igo[nr];
		p_a[0] = a[nr][0];
		p_a[1] = a[nr][1];
		p_a[2] = a[nr][2];
		p_a[3] = a[nr][3];
		p_b[0] = b[nr][0];
		p_b[1] = b[nr][1];
		p_b[2] = b[nr][2];
		p_b[3] = b[nr][3];
	}
	return;
}

void meteointerpol::putregresscoeffs(int nr,double *p_a,double *p_b,double p_igo)
{
	for (int j=0; j<=3;j++) {a[nr][j]=p_a[j]; b[nr][j]=p_b[j];}
	igo[nr] = p_igo;
}


void meteointerpol::flushbuffers()
{
	if (zonestat) zonestat->flushbuffer();
}

precicorrclass::precicorrclass(int &retvalue,wasim::mystream *Stream,meteointerpol **p_meteo,
	char **meteonamen,char *steuername,int LogLevel,int anzahlgroessen,
	int p_regiomax, int p_regiomaxwert, int* p_rarray,
	gridinfotype *p_gridinfo,
	firstlast *p_FirstLastCol ,
	int   p_ncols, int        p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll   , double      p_yll,
	double p_cellsize    , double      p_nodata,
	LUSEtabletype *p_LUtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory, long *v_mpiregioncnt, std::vector<mpihalo> *p_sendHalo, std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
	modelobject(Stream, steuername, LogLevel, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo, p_doSpinUp, p_vHRU),
do_it(0),meteo(0),Tgrenz(0),ar(0),as(0),br(0),bs(0),Phoehe(0),Px(0),Py(0),Pmeth(0),
Panz(0),regiomax(0),igu(0),igo(0),a(0),b(0),hmax(0),Twert(0),Uwert(0),Pwert(0),
Pindex(0),Tindex(0),Uindex(0),nowind(0),notemp(0),doGriddedCorrection(0),WIND(0),TEMP(0),routingmode(p_routingmode)
{
	try
	{
		retvalue=0;
		int nr,i;
		FILE *datei;

		do_it=1;
		/* Lesen der Korrektur-Parameter aus der Steuerdatei */

		gridinfo = p_gridinfo;
		FirstLastCol = p_FirstLastCol;
		ncols = p_ncols; 
		nrows = p_nrows;
		mpistartrow = p_mpistartrow;
		mpiendrow = p_mpiendrow;
		xll = p_xll;
		yll = p_yll;
		cellsize = p_cellsize; 
		nodata = p_nodata;
		
		meteo = p_meteo;
		regiomax=p_regiomax;
		if (regiomax<=0) regiomax=1;
		LUtable=p_LUtable;
		MultiLayerLUTable=p_MultiLayerLUTable;

		sprintf(instanzname,"%s",_Niederschlagskorrektur_);

		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}

		if (index(instanzname," ")>0)  instanzname = substr(instanzname,1,index(instanzname," ")-1);
		sprintf(tstr,"[%s]",instanzname);

		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		findeabschnitt(tstr);

		/* lesen der EInstellungsparameter */

		fscanfln(steuerdatei,"%d",&do_it   ,"do_it fuer "           ,tstr);
		if (routingmode)
		{
			do_it = 0; // in routingmode, no other modules than routing are permitted to run
		}
		if (do_it == 0) return;

		// new in 11.02.00: parameters per region

//		fscanfln(steuerdatei,"%lf",&Tgrenz  ,"Grenztemperatur fuer " ,tstr);
//		fscanfln(steuerdatei,"%lf",&br      ,"Konstante Regen fuer " ,tstr);
//		fscanfln(steuerdatei,"%lf",&ar      ,"Anstieg Regen fuer "   ,tstr);
//		fscanfln(steuerdatei,"%lf",&bs      ,"Konstante Schnee fuer ",tstr);
//		fscanfln(steuerdatei,"%lf",&as      ,"Anstieg Schnee fuer "  ,tstr);

		Tgrenz = new double[regiomax + 1];
		br = new double[regiomax + 1];
		ar = new double[regiomax + 1];
		bs = new double[regiomax + 1];
		as = new double[regiomax + 1];

		if (do_it == 1)
		{
			// traditional one-parameter correction
			fscanfln(steuerdatei,"%lf",&Tgrenz[0], "temperatur limit snow/rain for ", tstr);
			fscanfln(steuerdatei,"%lf",&br[0]      ,"correction constant rain for " ,tstr);
			fscanfln(steuerdatei,"%lf",&ar[0]      ,"correction factor rain for "   ,tstr);
			fscanfln(steuerdatei,"%lf",&bs[0]      ,"correction constant snow for ",tstr);
			fscanfln(steuerdatei,"%lf",&as[0]      ,"correction factor snow for "  ,tstr);

			for (int i = 1; i <= regiomax; i++)
			{
				Tgrenz[i] = Tgrenz[0];
				br[i] = br[0];
				ar[i] = ar[0];
				bs[i] = bs[0];
				as[i] = as[0];
			}

		}
		else
		{ // do_it=2, per region one parameter
			// Tgrenz
			for (nr = 1; nr < regiomax; nr++)
			{
				fscanf(steuerdatei, "%lf", &(Tgrenz[nr]));
			}
			// last number with an optional comment or line feed after the value
			fscanfln(steuerdatei, "%lf", &Tgrenz[regiomax], "temperatur limit snow/rain for all regions for", tstr);
			Tgrenz[0] = 0.0;

			// br - correction constant rain
			for (nr = 1; nr < regiomax; nr++)
			{
				fscanf(steuerdatei, "%lf", &(br[nr]));
			}
			// last number with an optional comment or line feed after the value
			fscanfln(steuerdatei, "%lf", &br[regiomax], "correction constant rain  for all regions for", tstr);
			br[0] = 0.0;

			// ar - correction factor rain
			for (nr = 1; nr < regiomax; nr++)
			{
				fscanf(steuerdatei, "%lf", &(ar[nr]));
			}
			// last number with an optional comment or line feed after the value
			fscanfln(steuerdatei, "%lf", &ar[regiomax], "correction factor rain  for all regions for", tstr);
			ar[0] = 1.0;

			// bs - correction constant snow
			for (nr = 1; nr < regiomax; nr++)
			{
				fscanf(steuerdatei, "%lf", &(bs[nr]));
			}
			// last number with an optional comment or line feed after the value
			fscanfln(steuerdatei, "%lf", &bs[regiomax], "correction constant snow for all regions for", tstr);
			bs[0] = 0.0;

			// as correction factor snow
			for (nr = 1; nr < regiomax; nr++)
			{
				fscanf(steuerdatei, "%lf", &(as[nr]));
			}
			// last number with an optional comment or line feed after the value
			fscanfln(steuerdatei, "%lf", &as[regiomax], "correction factor snow  for all regions for", tstr);
			as[0] = 1.0;
		}



		// all parameters read, closing control file for this class
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing"<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		igo = (double*)calloc(regiomax+1,sizeof(double));
		igu = (double*)calloc(regiomax+1,sizeof(double));
		a   = (double**)calloc(regiomax+1,sizeof(double*));
		for (nr = 1; nr <= regiomax; nr++)
		{
			a[nr] = (double*)calloc(4, sizeof(double));
		}
		b   = (double**)calloc(regiomax+1,sizeof(double*));
		for (nr = 1; nr <= regiomax; nr++)
		{
			b[nr] = (double*)calloc(4, sizeof(double));
		}

		Pindex=1;
		/* suchen des Indexes fuer Niederschlag im meteo[1..meteoanzahl]-Array */
		while ((strcmp(_Niederschlag_,meteonamen[Pindex])!=0)&&(Pindex<anzahlgroessen)) Pindex++;
		if  (strcmp(_Niederschlag_,meteonamen[Pindex])!=0)
		{
			// looking for a grid with this name, because it may be that superpositioning is used, so there is no interpolation instance with precipitation but an outputgrid of superposition
			doGriddedCorrection=true;
			if ((Pindex<=0)&&(mpi_rank==0))
			{
				do_it=0;
				char logfilename[1024];
				sprintf(logfilename,"modell_%d.log",mpi_rank);
				if ((datei = fopen(logfilename,"a"))==NULL)
				{
					(*WSstream) << "FATAL: error opening log-file "<<logfilename << "[end-entry]";
					throw(65);
				}
				fprintf(datei,"no precipitation correction will be done since %s was not found in %s!\n\n",_Niederschlag_,steuername);
				fclose(datei);
				return;
			}
		}

		// for grid based correction, these data should be available
		PRECIP = gridinfo[getindex(_Niederschlag_)].pointer;

		Tindex=1;
		/* suchen des Indexes fuer Temperatur im meteo[1..meteoanzahl]-Array */
		notemp=0;
		while ((strcmp(_Temperatur_,meteonamen[Tindex])!=0)&&(Tindex<anzahlgroessen)) Tindex++;
		if  (strcmp(_Temperatur_,meteonamen[Tindex])!=0)
		{
			// looking for a grid with this name, because it may be that superpositioning is used, so there is no interpolation instance with precipitation but an outputgrid of superposition
			if (ist_vorhanden(_Temperatur_)) 
			{
				doGriddedCorrection=true;
			}
			else	
			{
				Tindex=-1;
				notemp=1;
			}
		}
		if (!notemp)
		{
			TEMP = gridinfo[getindex(_Temperatur_)].pointer;
		}
		else
		{
			TEMP = NULL;
		}

		Uindex=1;
		nowind=0;
		/* suchen des Indexes fuer Windgeschwindigkeit im meteo[1..meteoanzahl]-Array */
		while ((strcmp(_Windgeschwindigkeit_,meteonamen[Uindex])!=0)&&(Uindex<anzahlgroessen)) Uindex++;
		if  (strcmp(_Windgeschwindigkeit_,meteonamen[Uindex])!=0)
		{
			// looking for a grid with this name, because it may be that superpositioning is used, so there is no interpolation instance with precipitation but an outputgrid of superposition
			if (ist_vorhanden(_Windgeschwindigkeit_)) 
			{
				doGriddedCorrection=true;
			}
			else 
			{
				Uindex=-1;
				nowind=1;
			}
		}
		if (!nowind)
		{
			WIND = gridinfo[getindex(_Windgeschwindigkeit_)].pointer;
		}
		else
		{
			WIND = NULL;
		}

		if (Uindex > 0)
		{
			if (meteo[Uindex]->getmethode() == 9)
			{
				doGriddedCorrection = true;
			}
		}

		if (Tindex > 0)
		{
			if (meteo[Tindex]->getmethode() == 9)
			{
				doGriddedCorrection = true;
			}
		}


		// parameter for location corrections may be taken from the interpolators now 
		// but only if no superposition of any variable used here is done (the distributed correction will then be called anyway)
		if (!doGriddedCorrection)
		{
			/* holen der Methoden und der (anzahl[0]) der Werte */
			Pmeth = meteo[Pindex]->getmethode();
			if (Pmeth == 9)
			{
				doGriddedCorrection = true;
			}
			else
			{
				if ((Pmeth != 2) && (Pmeth != 10))
				{
					Panz = meteo[Pindex]->getanzahl();
					/* holen der Adressen der Koordinaten und der Wertearrays, wenn diese gebraucht werden koennen */
					Phoehe = meteo[Pindex]->gethoehen();
					Px = meteo[Pindex]->getxkoords();
					Py = meteo[Pindex]->getykoords();
					Pwert = meteo[Pindex]->getwert();
				}
				else
				{
					Panz = 1;
				}
				Twert = (double *)calloc(Panz + 1, sizeof(double)); if (!Twert) speichererror(1);
				Uwert = (double *)calloc(Panz + 1, sizeof(double)); if (!Uwert) speichererror(1);
				for (i = 1; i <= Panz; i++)
				{
					Twert[i] = 0.0;
					Uwert[i] = 0.0;
				}
			}
		}

		// getting regions and weights
		regiomax = p_regiomax;
		regiomaxwert = p_regiomaxwert;
		rarray = p_rarray;

		regio1index = regio2index = regio3index = weight1index = weight2index = weight3index = 0;
		WEIGHT = NULL;
		REGION = NULL;
		if (rarray) /* es wurde ein Regionengrid eingelesen und ausgewertet, es gibt also mindestens eine Region */
		{
			WEIGHT = (double***)calloc(4, sizeof(double**));
			REGION = (double***)calloc(4, sizeof(double**));
			REGION[1] = REGION[2] = REGION[3] = WEIGHT[1] = WEIGHT[2] = WEIGHT[3] = NULL;
			regio1index = ist_vorhanden(_regression_regions_1_);
			weight1index = ist_vorhanden(_regression_weights_1_);
			REGION[1] = gridinfo[regio1index].pointer;
			WEIGHT[1] = gridinfo[weight1index].pointer;
			if (regiomax > 1)
			{
				regio2index = ist_vorhanden(_regression_regions_2_);
				weight2index = ist_vorhanden(_regression_weights_2_);
				REGION[2] = gridinfo[regio2index].pointer;
				WEIGHT[2] = gridinfo[weight2index].pointer;
			}
			if (regiomax > 2)
			{
				regio3index = ist_vorhanden(_regression_regions_3_);
				weight3index = ist_vorhanden(_regression_weights_3_);
				REGION[3] = gridinfo[regio3index].pointer;
				WEIGHT[3] = gridinfo[weight3index].pointer;
			}
		}
		else
		{
			regiomax = 0;
			REGION = WEIGHT = NULL;
		}

	}
	catch (int ex) 
	// all exceptions must be caught within the constructor, because throwing an exception 
	// outside off this constructor would cause the object to remain uninitialized
	// resulting in the inpossibility of destructing (i.e., the memory clearing) of the 
	// instance and hence in memory leaks
	{
		retvalue=ex;
	}
	catch (...)
	// also unknown exceptions have to be caught...
	{
		retvalue=UNKNOWN_EXCEPTION;
	}
}

/* Korrektur des Niederschlages anhand der Daten von Wind und Temperatur */
void precicorrclass::run()
{
	int i,nr;
	double dh;
	double hr;
	double Sx,Sy,Sxy,Sx2,Sy2,N;
	double P3,P2,P1;
	double xcen,ycen;
	double *entf;

	if (do_it==0) return;
	if (doGriddedCorrection) return;
	/* ##################################################################################### */
	/*                                                                                       */
	/* Korrektur der Stationswerte des Niederschlages, wenn die Niederschlagsinterpolation   */
	/* nicht ausschliesslich durch hoehenabhaengige Regression durchgefuehrt werden          */
	/*                                                                                       */
	/* ##################################################################################### */

	int Uanz;
	if (Uindex>0)
		Uanz=meteo[Uindex]->getanzahl();
	else
		Uanz=0;
	int Tanz=meteo[Tindex]->getanzahl();

	int maxstat=Panz;
	if (Tanz>maxstat) maxstat=Tanz;
	if (Uanz>maxstat) maxstat=Uanz;
	entf=new double[maxstat+1];

	if ((Pmeth!=2)&&(Pmeth!=10))
	{
		/* Berechnen der Stationswerte von Wind und Temperatur an den Niederschlags-Stationen */
		if (!notemp)
		{
			for (i=1;i<=Panz;i++)
			{
				/* mit meteo[index]->getsinglepoint(x,y,h) kann ein einzelner Wert einer durch meteo[<index>]
		spezifizierten Instanz der Interpolationsklasse an der STelle mit den Koordinaten x,y und
		der Hoehe h berechnet werden. <index> ist dabei der Index im array meteo, der fuer die
		zu interpolierende Groesse gilt (je nach Anordnung der Interpolationsabschnitte in der
		Steuerdatei). Die Interpolation wird unter Beruecksichtigung aller in der Steuerdatei
		festgelegten Parameter nur fuer diesen einen Punkt ausgefuehrt. */

				Twert[i]=meteo[Tindex]->getsinglepoint(Px[i],Py[i],Phoehe[i],entf);
				/* liefert Temperatur-Interpolation einen Fehlwert, wird Regen angenommen, also value > Tgrenz gesetzt */
				if (Twert[i] < -meteo[Tindex]->fehlwert) Twert[i]=Tgrenz[0] + 1.0;
			}
		}
		if (!nowind)
		{
			for (i=1;i<=Panz;i++)
			{
				Uwert[i]=meteo[Uindex]->getsinglepoint(Px[i],Py[i],Phoehe[i],entf);
				/* liefert die Windinterpolation einen Fehlwert, so wird Wind=0 gesetzt, somit erfolgt die Korrektur,
		als ob gar keine Windinterpolation stattgefunden haette */
				if (Uwert[i] < -meteo[Uindex]->fehlwert) Uwert[i]=0.0;
			}
		}

		/* die eigentliche Korrektur */
		if (!notemp) /* Temperatur steht zur Verfuegung */
		{
			/* Temperatur und Wind stehen zur Verfuegung, wenn kein Wind gemessen wurde,
		ist Uwert immer 0, die Korrektur erfolgt dann nur anhand der konstanten Faktoren br/bs */
			for (i=1;i<=Panz;i++)
			if (Twert[i] > Tgrenz[0])
			{ 
				if (Pwert[i] >= -meteo[Pindex]->fehlwert)
				{
					Pwert[i] = Pwert[i] * (br[0] + ar[0] * Uwert[i]);
				}
			}
			else
			{ 
				if (Pwert[i] >= -meteo[Pindex]->fehlwert)
				{
					Pwert[i] = Pwert[i] * (bs[0] + as[0] * Uwert[i]);
				}
			}
		}
		else
		{ /* Wenn keine Temperatur zur Verfuegung steht, dann wird immer angenommen, dass es sich um Regen handelt
		die Korrektur erfolgt dann nur durch die Windwerte. Steht auch der Wind nicht zur Verfuegung, sind alle
		Uwerte == 0, so dass nur der Korrekturfaktor br zur Anwendung kommt */
			for (i=1;i<=Panz;i++)
			if (Pwert[i] >= -meteo[Pindex]->fehlwert) Pwert[i] = Pwert[i]*(br[0] + ar[0] * Uwert[i]);
		}
	}

	/* ##################################################################################### */
	/*                                                                                       */
	/* Korrektur der Niederschlagsregressionsparameter, wenn die Niedserschlagsinterpolation */
	/* teilweise oder ganz durch hoehenabhaengige Regression durchgefuehrt wird              */
	/*                                                                                       */
	/* ##################################################################################### */

	// since 8.05.00 the regions are used differently. The precipitation correction
	// may use one single region only (since there are no regions in an interpolator
	// but only multiple interpolators for precipitation. This can only be handled by another
	// correction method (see below runOnSuperposedRegions()). Here, we have to make sure that only 
	// one region is used (for future reference on how the correction worked with multiple regions, we 
	// do not rewrite the entire procedure
	regiomax=1; 

	for (nr=1;nr<=regiomax;nr++)
	{
		if ((Pmeth==2)||(Pmeth==3)||(Pmeth==10)||(Pmeth==11))
		{
			/* holen der Parameter  */
			meteo[Pindex]->getregresscoeffs(nr,&(igu[nr]),&(igo[nr]),a[nr],b[nr]);
			/* zentrale Koordinaten für die hoehenabhaengige Regression */
			xcen = meteo[Pindex]->xll+meteo[Pindex]->cellsize*meteo[Pindex]->ncols/2.0;
			ycen = meteo[Pindex]->yll+meteo[Pindex]->cellsize*meteo[Pindex]->nrows/2.0;


			/* Wenn die Regression ungueltig ist, sollen die Parameter nicht veraendert werden
		Das sollte aber eigentlich nie vorkommen */
			if (a[nr][1] >= -meteo[Pindex]->fehlwert)
			{
				/* ------------------------------------------------------------------ */
				/* Abschnittsweise Korrektur des Hoehenprofiles, unterster Abschnitt, */
				/* es wird jeweils die Temperatur und die Windgeschwindigkeit fuer    */
				/* ein Hoehenprofil in der Mitte des Grids errechnet (je Geradenab-   */
				/* schnitt 10 Unterteilungen) und anhand dieser Werte korrigiert.     */
				/* ------------------------------------------------------------------ */
				if ((igu[nr]>0.0)&&(b[nr][1]!=0.0)) /* Wenn der Anstieg 0 ist, faellt kein Niederschlag, es muss und darf nicht korrigiert werden */
				{
					dh = igu[nr]/10.0;   /* Unterteilung des Bereiches in 10 einzelbereiche, um den Einfluss der Temperatur etwas besser zu erfassen */
					/* Statistikregister fuer die aus den Regressionsgeraden berechneten und korrigierten Niederschlaege 0 setzen */
					Sx = Sy = Sxy = Sx2 = Sy2 = N = 0;
					/* Statistikregister fuer jeden der 11 Teilbereiche fuellen */
					for (hr = 0.0; hr <= igu[nr]; hr = hr + dh)
					{
						/* Berechnen des Niederschlages anhand der uebergebenen Parameter a und b */
						P1   = a[nr][1] + b[nr][1]*hr;
						/* Berechnen der Temperatur un des Windes an dieser Stelle wie bei P-Interpolation durch IDW, Erlaeuterung siehe oben */
						if (!notemp) {
//							Twert[1]=meteo[Tindex]->getsinglepoint(xcen,ycen,hr,entf);
							Twert[1]=meteo[Tindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert Temperatur-Interpolation einen Fehlwert, wird Regen angenommen, also value > Tgrenz gesetzt */
							if (Twert[1] < -meteo[Tindex]->fehlwert) Twert[1]=Tgrenz[0] + 1.0;
						}
						else Twert[1]=Tgrenz[0] + 1.0;
						if (!nowind) {
//							Uwert[1]=meteo[Uindex]->getsinglepoint(xcen,ycen,hr,entf);
							Uwert[1]=meteo[Uindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert die Windinterpolation einen Fehlwert, so wird Wind=0 gesetzt, somit erfolgt die Korrektur,
				als ob gar keine Windinterpolation stattgefunden haette */
							if (Uwert[1] < -meteo[Uindex]->fehlwert) Uwert[1]=0.0;
						} else Uwert[1]=0.0;
						/* Korrektur des Niederschlagswertes */
						if (Twert[1] > Tgrenz[0])  P1 = P1 * (br[0] + ar[0] * Uwert[1]); else P1 = P1 * (bs[0] + as[0] * Uwert[1]);
						/* Registerupdate */
						Sx  += hr;
						Sy  += P1;
						Sxy += P1*hr;
						Sx2 += hr*hr;
						Sy2 += P1*P1;
						N++;
					}
					/* Parameter a und b fuer untersten Abschnitt nach linearer Regression Modell 1 neu berechnen */
					b[nr][1] = (double)((N*Sxy-Sx*Sy)/(N*Sx2-Sx*Sx));
					if (N>0) a[nr][1] = (double)((Sy-b[nr][1]*Sx)/N);
				}
				else /* wenn es nur eine oder keine Inversion gibt, dann gilt die mittlere Gerade bis zu igo, welches dann anstelle von igu tritt */
				{    /* igu ist 0 und b[1] ist auch 0 */
					dh = igo[nr]/10.0;   /* Unterteilung des Bereiches in 10 einzelbereiche, um den Einfluss der Temperatur etwas besser zu erfassen */
					/* Statistikregister fuer die aus den Regressionsgeraden berechneten und korrigierten Niederschlaege 0 setzen */
					Sx = Sy = Sxy = Sx2 = Sy2 = N = 0;
					/* Statistikregister fuer jeden der 11 Teilbereiche fuellen */
					for (hr = 0.0; hr <= igo[nr]; hr = hr + dh)
					{
						/* Berechnen des Niederschlages anhand der uebergebenen Parameter a und b */
						P1   = a[nr][2] + b[nr][2]*hr;
						/* Berechnen der Temperatur un des Windes an dieser Stelle wie bei P-Interpolation durch IDW, Erlaeuterung siehe oben */
						if (!notemp) {
//							Twert[1]=meteo[Tindex]->getsinglepoint(xcen,ycen,hr,entf);
							Twert[1]=meteo[Tindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert Temperatur-Interpolation einen Fehlwert, wird Regen angenommen, also value > Tgrenz gesetzt */
							if (Twert[1] < -meteo[Tindex]->fehlwert) Twert[1]=Tgrenz[0] + 1.0;
						}
						else Twert[1]=Tgrenz[0] + 1.0;
						if (!nowind) {
//							Uwert[1]=meteo[Uindex]->getsinglepoint(xcen,ycen,hr,entf);
							Uwert[1]=meteo[Uindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert die Windinterpolation einen Fehlwert, so wird Wind=0 gesetzt, somit erfolgt die Korrektur,
				als ob gar keine Windinterpolation stattgefunden haette */
							if (Uwert[1] < -meteo[Uindex]->fehlwert) Uwert[1]=0.0;
						} else Uwert[1]=0.0;
						/* Korrektur des Niederschlagswertes */
						if (Twert[1] > Tgrenz[0])  P1 = P1 * (br[0] + ar[0] * Uwert[1]); else P1 = P1 * (bs[0] + as[0] * Uwert[1]);
						/* Registerupdate */
						Sx  += hr;
						Sy  += P1;
						Sxy += P1*hr;
						Sx2 += hr*hr;
						Sy2 += P1*P1;
						N++;
					}
					/* Parameter a und b fuer mittleren Abschnitt nach linearer Regression Modell 1 neu berechnen */
					b[nr][2] = (double)((N*Sxy-Sx*Sy)/(N*Sx2-Sx*Sx));
					if (N>0) a[nr][2] = (double)((Sy-b[nr][2]*Sx)/N);
				}

				/* ---------------------------------------------------------------- */
				/*                                                                  */
				/* Abschnittsweise Korrektur des Hoehenprofiles, oberster Abschnitt */
				/*                                                                  */
				/* ---------------------------------------------------------------- */

				if ((b[nr][3]!=0.0)&&(hmax>igo[nr]))
				{
					dh = (hmax-igo[nr])/10.0;
					Sx = Sy = Sxy = Sx2 = Sy2 = N = 0;
					/* Statistikregister fuellen */
					for (hr = igo[nr]; (hr <= 2.0*igo[nr]) || (hr <= hmax); hr = hr + dh)
					{
						P1   = a[nr][3] + b[nr][3]*hr;
						/* Berechnen der Temperatur un des Windes an dieser Stelle */
						if (!notemp) {
//							Twert[1]=meteo[Tindex]->getsinglepoint(xcen,ycen,hr,entf);
							Twert[1]=meteo[Tindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert Temperatur-Interpolation einen Fehlwert, wird Regen angenommen, also value > Tgrenz gesetzt */
							if (Twert[1] < -meteo[Tindex]->fehlwert) Twert[1]=Tgrenz[0] + 1.0;
						}
						else Twert[1]=Tgrenz[0] +1.0;
						if (!nowind) {
//							Uwert[1]=meteo[Uindex]->getsinglepoint(xcen,ycen,hr,entf);
							Uwert[1]=meteo[Uindex]->getsinglepoint(nodata,nodata,hr,entf);
							/* liefert die Windinterpolation einen Fehlwert, so wird Wind=0 gesetzt, somit erfolgt die Korrektur,
				als ob gar keine Windinterpolation stattgefunden haette */
							if (Uwert[1] < -meteo[Uindex]->fehlwert) Uwert[1]=0.0;
						} else Uwert[1]=0.0;
						/* Korrektur des Niederschlagswertes */
						if (Twert[1] > Tgrenz[0])  P1 = P1*(br[0] + ar[0] * Uwert[1]); else P1 = P1*(bs[0] + as[0] * Uwert[1]);
						Sx  += hr;
						Sy  += P1;
						Sxy += P1*hr;
						Sx2 += hr*hr;
						Sy2 += P1*P1;
						N++;
					}
					/* Parameter a und b fuer obersten Abschnitt neu berechnen */
					b[nr][3] = (double)((N*Sxy-Sx*Sy)/(N*Sx2-Sx*Sx));
					if (N>0) a[nr][3] = (double)((Sy-b[nr][3]*Sx)/N);
				}

				/* --------------------------------------------------- */
				/*                                                     */
				/* Verbinden des obersten und des untersen Abschnittes */
				/* (zwei Inversionen, alle drei Geraden gueltig)      */
				/* oder suchen des neuen Schnittpunktes, wenn nur zwei */
				/* Geraden galten (max. eine Inversion)                */
				/*                                                     */
				/* --------------------------------------------------- */

				if (igu[nr]>0.0)
				{ /* verbinden obere und untere Gerade */
					P1 = a[nr][1] + b[nr][1]*igu[nr];
					P3 = a[nr][3] + b[nr][3]*igo[nr];
					P2 = P3 - P1;
					b[nr][2] = P2/(igo[nr]-igu[nr]);
					a[nr][2] = a[nr][1] + (b[nr][1]-b[nr][2])*igu[nr];
				}
				else
				{ /* finden Schnittpunkt igo */
					if (fabs(b[nr][3]-b[nr][2])>0.0001)
					igo[nr] = (a[nr][2]-a[nr][3])/(b[nr][3]-b[nr][2]);
				}
			}
			/* zurueckschreiben der Parameter in die Felder des Niederschlagsinterpolators */
			meteo[Pindex]->putregresscoeffs(nr,a[nr],b[nr],igo[nr]);
		}
		// attention: methods 10 and 11 are for one region only, so these methods should not loop through all regions but the first one (default region if no region grid is defined)
		if ((Pmeth==10)||(Pmeth==11)) nr=regiomax;
	}
	if (entf) delete[] entf;
}

void precicorrclass::runOnSuperposedRegions()
{
	if (do_it==0) return;

	if (MPIREG == NULL)
	{
		if (ist_vorhanden(_MPIRegions_))
			MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
		else
			MPIREG = NULL;
	}

	int i,j;
#pragma omp parallel for schedule(runtime) private(i,j)
	for (int i=mpistartrow;i<=mpiendrow;i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if (MPIREG[i][j] == mpi_rank)
			{
				/* do the correction here */
				// get the regions, at max a point is assigned to three regions
				int* reg = new int[4] {-9999, rarray[(int)REGION[1][i][j]], 0, 0};
				double* weight = new double[4] {-9999.0, WEIGHT[1][i][j], 0.0, 0.0};

				if (regiomax>1)
				{ 
					if ((int)REGION[2][i][j] > 0) reg[2] = rarray[(int)REGION[2][i][j]];
					if (reg[2] > 0)
					{
						weight[2] = WEIGHT[2][i][j];
					}

					if (regiomax > 2)
					{
						if ((int)REGION[3][i][j] > 0) reg[3] = rarray[(int)REGION[3][i][j]];
						if (reg[3] > 0)
						{
							weight[3] = WEIGHT[3][i][j];
						}
					}
				}
				// no we have reg1, reg2 and reg3 with their respective weights weight1, weight2 and weight3. If a cell has only one (or two) regions assigned, then reg2 and reg3 (or only reg3) are 0 and weight2 and weight3 (or only weight3) are 0
				// thus, we can use all reg-indices and weights within a single formula
				// -----------------------------------------------------------------------------
				// next step: getting temperature weight (Tgrenz minus 1 is only snow, Tgrenz+1 is only rain, 
				// everything between is a certain amount rain and a certain amount snow, 
				// assuming a linear transition we need to interpolate linear between both points
				
				// default values for rain and snow fraction weights (per region, since Tgrenz may vary from region to region)
				double* ws = new double[4] {-9999.0, 0.0, 0.0, 0.0}; // weight for snow fraction
				double* wr = new double[4] {-9999.0, 1.0, 1.0, 1.0}; // weight for rain fraction
				
				// constants for better readibility of the following formulas
				const double maxweightrain = 1.0;
				const double minweightrain = 0.0;

				double T = notemp ? 0.0: TEMP[i][j];
				double rainoffset = 1.0; // one degree celsius more than Tgrenz is the limit for any snow
				double snowoffset = -1.0; // one degree celsius less than Tgrenz is the limit for any rain
				
				double transitionrange = rainoffset - snowoffset;

				if (!notemp)
				{
					for (int r = 1; r <= 3; r++)
					{

						if (T < Tgrenz[reg[r]] - 1.0)
						{ // temperature below threshold for pure snow
							ws[r] = 1.0;
							wr[r] = 0.0;
						}
						else if (T < Tgrenz[reg[r]] + 1.0)
						{ // temperature between snow and rain, do a linear interpolation 
							wr[r] = (T - (Tgrenz[reg[r]] + snowoffset)) / (rainoffset-snowoffset);
							ws[r] = 1.0 - wr[r];
						}
					}
				}
				else
				{ // no temperature information, we assume always rain
					for (int r = 1; r <= 3; r++)
					{
						ws[r] = 0.0;
						wr[r] = 1.0;
					}
				}
				double nodata = -meteo[Pindex]->fehlwert; // e.g. 9998 is defined in the control file, so the actual value is -9998, i.e. -9999 counts as nodata
				double u = (WIND ? WIND[i][j] : 0.0);

				double f1 = (br[reg[1]] + ar[reg[1]] * u) * weight[1] * wr[1]; // weight for rain for region 1 with parameters for region 1
				double f2 = (br[reg[2]] + ar[reg[2]] * u) * weight[2] * wr[2]; // weight for rain for region 1 with parameters for region 2
				double f3 = (br[reg[3]] + ar[reg[3]] * u) * weight[3] * wr[3]; // weight for rain for region 1 with parameters for region 3
				double f4 = (bs[reg[1]] + as[reg[1]] * u) * weight[1] * ws[1]; // weight for rain for region 1 with parameters for region 1
				double f5 = (bs[reg[2]] + as[reg[2]] * u) * weight[2] * ws[2]; // weight for rain for region 1 with parameters for region 2
				double f6 = (bs[reg[3]] + as[reg[3]] * u) * weight[3] * ws[3]; // weight for rain for region 1 with parameters for region 3

				PRECIP[i][j] = (PRECIP[i][j] >= nodata ? PRECIP[i][j] * (f1+f2+f3+f4+f5+f6) : PRECIP[i][j]);

//				if (!notemp) /* Temperatur steht zur Verfuegung */
//				{
//					/* Temperatur und Wind stehen zur Verfuegung, wenn kein Wind gemessen wurde,
//				ist Uwert immer 0, die Korrektur erfolgt dann nur anhand der konstanten Faktoren br/bs */
//					if (TEMP[i][j] > Tgrenz[0])
//					{
//						if (WIND)
//						{
//							if (PRECIP[i][j] >= -meteo[Pindex]->fehlwert) PRECIP[i][j] = PRECIP[i][j] * (br[0] + ar[0] * WIND[i][j]);
//						}
//						else
//						{
//							if (PRECIP[i][j] >= -meteo[Pindex]->fehlwert) PRECIP[i][j] = PRECIP[i][j] * (br[0] + ar[0] * 0);
//						}
//					}
//					else
//					{
//						if (WIND)
//						{
//							if (PRECIP[i][j] >= -meteo[Pindex]->fehlwert) PRECIP[i][j] = PRECIP[i][j] * (bs[0] + as[0] * WIND[i][j]);
//						}
//						else
//						{
//							if (PRECIP[i][j] >= -meteo[Pindex]->fehlwert) PRECIP[i][j] = PRECIP[i][j] * (bs[0] + as[0] * 0);
//						}
//					}
//				}
//				else
//				{ /* Wenn keine Temperatur zur Verfuegung steht, dann wird immer angenommen, dass es sich um Regen handelt
//				die Korrektur erfolgt dann nur durch die Windwerte. Steht auch der Wind nicht zur Verfuegung, sind alle
//				Uwerte == 0, so dass nur der Korrekturfaktor br zur Anwendung kommt */
//					if (PRECIP[i][j] >= -meteo[Pindex]->fehlwert) PRECIP[i][j] = PRECIP[i][j] * (br[0] + ar[0] * WIND[i][j]);
//				}
			}
		}
	}
}


void precicorrclass::put_hmax(double p_hmax)
{
	hmax = p_hmax;
}

/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void precicorrclass::ClearMemory()
{
	int nr;
	if (do_it > 0)
	{
		if (Twert) free(Twert); Twert=NULL;
		if (Uwert) free(Uwert); Uwert=NULL;
		if (igo) free(igo);     igo=NULL;
		if (igu) free(igu);     igu=NULL;
		for (nr=1;nr<=regiomax;nr++)
		{
			if (a[nr]) 
			free(a[nr]); 
			a[nr]=NULL;
			if (b[nr]) 
			free(b[nr]); 
			b[nr]=NULL;
		}
		if (a) free(a); 
		a=NULL;
		if (b) free(b); 
		b=NULL;
		if (WEIGHT) free(WEIGHT); WEIGHT = NULL;
		if (REGION) free(REGION); REGION = NULL;

		if (Tgrenz) { delete[] Tgrenz;	Tgrenz = NULL; }
		if (ar) delete[] ar;	ar = NULL;
		if (br) delete[] br;	br = NULL;
		if (as) delete[] as;	as = NULL;
		if (bs) delete[] bs;	bs = NULL;

	}
	modelobject::ClearMemory();
}

// implementations for regional superposition of interpolation results
// this has the status of a sub model but must be run before any other submodel after interpolation 
// (and before new precipitation correction using interpolated temperatures, wind and precip data)

RegionalSuperposition::RegionalSuperposition(int &retvalue,wasim::mystream *p_Stream,TLastResults *LastResults,int maxLastResults,
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
	long mine,long ste , long tge , long mte , long jre,double p_timekorr,
	LUSEtabletype *p_LUtable,
	SOILtabletype *p_SOILtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory, long *v_mpiregioncnt, std::vector<mpihalo> *p_sendHalo, std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(p_Stream,breite,laenge,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,numsubb,listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols,p_nrows,p_mpistartrow,p_mpiendrow,p_xll,p_yll,p_cellsize,p_nodata,mine,ste,tge,mte,jre,p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo, p_doSpinUp, p_routingmode,p_vHRU)
{
	int index1,index2;

	if (this->do_it==0) return;

	this->LastResults    = LastResults;
	this->maxLastResults = maxLastResults;
	regiomax    = p_regiomax;
	regiomaxwert = p_regiomaxwert;
	rarray       = p_rarray;

	// read number of entities
	tstr=fgets(tstr,_stringsize_,steuerdatei);
	std::string record(tstr);
	if (record.find("NumberOfEntities")!=std::string::npos)
	{
		index1 = record.find("=")+1;
		index2 = record.find(";");
		record=record.substr(index1,index2-index1);
		iNumberEntities = atoi( record.c_str());
	}
	else
	{
		(*WSstream) << "FATAL: could not read \"NumberOfEntities\" keyword in section ["<<this->instanzname<<"]\n\n"<<"[end-entry]";
		throw 133;
	}


	for (int i=1;i<=iNumberEntities;i++)
	{
		tstr=fgets(tstr,_stringsize_,steuerdatei);
		std::string record(tstr);
		while (record.find("}")==std::string::npos)
		{
			tstr=fgets(tstr,_stringsize_,steuerdatei);
			record.append(std::string(tstr));
		}
		this->vSuperpos.push_back(new cSuperposition(this,record));
	}
	
}

void RegionalSuperposition::run(long me,long se,long tg,long mt,long jr)
{
	jahr=jr;
	monat=mt;
	tag=tg;
	stunde=se;
	minute = me;

	if (this->do_it==0) return;
	for (int i=0;i<this->vSuperpos.size();i++)
	{
		vSuperpos[i]->run(me,se,tg,mt,jr);
	}
}

void RegionalSuperposition::ClearMemory()
{
	if (this->do_it==0)

        {
            teilmodellklasse::ClearMemory();
            return;
        }
	for (int i=0;i<this->vSuperpos.size();i++)
	{
		vSuperpos[i]->ClearMemory();
	}
	vSuperpos.clear();
	teilmodellklasse::ClearMemory();

}

void RegionalSuperposition::flushbuffers()
{
	if (this->do_it==0) return;
	for (int i=0;i<this->vSuperpos.size();i++)
	{
		vSuperpos[i]->flushbuffers();
	}
}

void RegionalSuperposition::WriteResultGrids()
{
	if (this->do_it==0) return;
	for (int i=0;i<this->vSuperpos.size();i++)
	{
		this->vSuperpos[i]->WriteResultGrids();
	}
}




cSuperposition::cSuperposition(RegionalSuperposition *p_parent,std::string p_initials):
parent(p_parent),
initials(p_initials),
ncols(p_parent->ncols),
nrows(p_parent->nrows),
mpistartrow(p_parent->mpistartrow),
mpiendrow(p_parent->mpiendrow),
xll(p_parent->xll),
yll(p_parent->yll),
cellsize(p_parent->cellsize),
nodata(p_parent->nodata),
regiomax(p_parent->regiomax),
regiomaxwert(p_parent->regiomaxwert),
rarray(p_parent->rarray),
zarray(p_parent->zarray),
FirstLastCol(p_parent->FirstLastCol),
zn(p_parent->zn),
jahr(0),
monat(0),
tag(0),
stunde(0),
MPIREG(0),MPIHLP(0)
{
	// now the 
	// parameters have to be parsed
	std::string buffer(p_initials);
	int index1,index2;
	index1 = buffer.find_first_not_of(" \t");
	index2 = buffer.find("{");
	buffer=buffer.substr(index1,index2-index1);
	std::stringstream ss;
	ss << buffer.c_str();
	ss >> sEntityName;

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#else
	mpi_rank=0;
#endif

	if (MPIREG==NULL)
	{
		if (parent->ist_vorhanden(_MPIRegions_))
			MPIREG=parent->gridinfo[parent->getindex(_MPIRegions_)].pointer;
		else
			MPIREG=NULL;
	}

	if (MPIHLP==NULL)
	{
		if (parent->ist_vorhanden(_MPIHelper_))
			MPIHLP=parent->gridinfo[parent->getindex(_MPIHelper_)].pointer;
		else
			MPIHLP=NULL;
	}

	buffer = p_initials.substr(index2+1,p_initials.find("}")-index2+2);
	std::string subbuffer;
	while (buffer.find("entityinputgrid")!=std::string::npos)
	{
		inputstruct *is = new inputstruct;
		buffer = buffer.substr(buffer.find("entityinputgrid")+15);
		index1=buffer.find("=")+1;
		index2=buffer.find(";");
		subbuffer=buffer.substr(index1,index2-index1);
		std::stringstream ss; 
		ss.clear();
		ss << subbuffer;
		ss >> is->sInputEntity;
		is->ingridindex=parent->getindex(is->sInputEntity.c_str());
		is->INGRID=parent->gridinfo[is->ingridindex].pointer;

		index1=buffer.find("regions")+8; // cut off the entity name and the regions keyword
		buffer = buffer.substr(index1);
		index1=buffer.find("=")+1;
		index2=buffer.find(";");
		subbuffer = buffer.substr(index1,index2-index1);
		std::vector<std::string> words;
		int nReg=su.StringSplit(subbuffer,"\n\t ;,",words,false);
		if (nReg!=regiomax)
		{
			(*parent->WSstream) << "\nerror reading regions in "<< parent->instanzname << 
			" entity: " << sEntityName.c_str() << ", inputgrid " << is->sInputEntity.c_str() 
			<< "\nThe number of regions here ("<<nReg<<") does not equal the number of"
			<< "\nregions in the regiongrid ("<<regiomax<<")\n"<< "[end-entry]";
			throw 134;
		}
		for (int i=0;i<nReg;i++)
		{
			int reg = atoi(words[i].c_str());
			if (rarray[reg]!=(i+1))
			{
				if (rarray[i+1]==0)
				{
					(*parent->WSstream) << "\nerror reading regions in "<< parent->instanzname << 
					" entity: " << sEntityName.c_str() << ", inputgrid " << is->sInputEntity.c_str() 
					<< "\nregion on position "<<(i+1)<<" was read as "<<reg
					<<"\nThere is now region with this code in the region grid\n"
					<<"\nOnly the region codes from the region grid must be listet in ascending order here\n"
					<<"regardless of their actual weight for this superposition instance\n"<<"[end-entry]";
				}
				else
				{
					(*parent->WSstream) << "\nerror reading regions in "<< parent->instanzname << 
					" entity: " << sEntityName.c_str() << ", inputgrid " << is->sInputEntity.c_str() 
					<< "\nregion on position "<<(i+1)<<" in the region list should have value "<<rarray[i+1]
					<<"\nbut was read as "<<reg<<"\nAll regions must be listet in ascending order here\n"
					<<"regardless of their actual weight for this superposition instance\n"<<"[end-entry]";
				}
				throw 134;
			}
			is->region.push_back(reg);
		}

		index1=buffer.find("weights")+8; // cut off the weights keyword
		buffer = buffer.substr(index1);
		index1=buffer.find("=")+1;
		index2=buffer.find(";");
		subbuffer = buffer.substr(index1,index2-index1);
		words.clear();
		int nW=su.StringSplit(subbuffer,"\n\t ;,",words,false);
		if (nReg!=nW)
		{
			(*parent->WSstream) << "\nerror reading weights for regions in "<< parent->instanzname << 
			" entity: " << sEntityName.c_str() << ", inputgrid " << is->sInputEntity.c_str() 
			<< "\nThe number of regions ("<<nReg<<") does not equal the number of weights ("<<nW<<")\n"<< "[end-entry]";
			throw 135;
		}
		for (int i=0;i<nW;i++)
		{
			is->weight.push_back(atof(words[i].c_str()));
		}

		inputs.push_back(is);
	}

	// check if the single weights will add up to 1.0 for each region
	for (int rc=1;rc<=regiomax;rc++)
	{
		double regweight=0.0;
		for (int ic=0;ic<inputs.size();ic++)
		{
			regweight+=inputs[ic]->weight[rc-1];
		}
		//if (fabs(regweight-1.0)>0.0001)
		//{
		//	(*parent->WSstream) << "\nerror reading weights for regions in "<< parent->instanzname << 
		//	"\nentity: " << sEntityName.c_str() << 
		//	"\nregion "<<rc<<": The sum of all weights is not 1.0 but "<<regweight<<" \n"<< "[end-entry]";
		//	throw 135;
		//}
	}

	// now reading the output grid name
	index1=buffer.find("outputgrid")+11;
	subbuffer = buffer.substr(index1);
	index1=subbuffer.find("=")+1;
	index2=subbuffer.find(";");
	subbuffer = subbuffer.substr(index1,index2-index1);
	ss.clear();
	ss << subbuffer.c_str();
	ss >> resultgridname;

	// now reading the output grid write code
	index1=buffer.find("writecode")+10;
	subbuffer = buffer.substr(index1);
	index1=subbuffer.find("=")+1;
	index2=subbuffer.find(";");
	subbuffer = subbuffer.substr(index1,index2-index1);
	ss.clear();
	ss << subbuffer.c_str();
	ss >> gridwritecode;

	// now reading the output table name
	index1=buffer.find("outputtable")+12;
	subbuffer = buffer.substr(index1);
	index1=subbuffer.find("=")+1;
	index2=subbuffer.find(";");
	subbuffer = subbuffer.substr(index1,index2-index1);
	ss.clear();
	ss << subbuffer.c_str();
	ss >> statfilename;

	// now reading the output table statcode
	index1=buffer.find("statcode")+9;
	subbuffer = buffer.substr(index1);
	index1=subbuffer.find("=")+1;
	index2=subbuffer.find(";");
	subbuffer = subbuffer.substr(index1,index2-index1);
	ss.clear();
	ss << subbuffer.c_str();
	ss >> statwritecode;

	// the result grid must be created
	resultindex=1;
	while (parent->gridinfo[resultindex].identification!=NULL)
	{
		resultindex ++;
		if (resultindex==gridmaxanzahl)
		{
			(*parent->WSstream) << "no free grid anymore, allowed number of grids is "<<gridmaxanzahl<<"!"<<"[end-entry]";
			throw(12);
		}
	}

	//	double **RESULT;
	//	int gridwritecode,resultindex;
	//	std::string resultgridname;

	RESULT = parent->matrix(mpistartrow,mpiendrow,1,ncols,nrows);
	parent->initialize_grid(resultindex,resultgridname.c_str(),RESULT,NULL,0,sEntityName.c_str(),gridwritecode,mpistartrow,mpiendrow);
	for (int k = mpistartrow; k <= mpiendrow; k++)
	{
		for (int j = 1; j <= ncols; j++)
		{
//			if (MPIREG[k][j] == mpi_rank)
//			{
//				RESULT[k][j] = parent->zn[k][j] == nodata ? nodata : 0;
//			}
			RESULT[k][j] = parent->zn[k][j] == nodata ? nodata : 0;
		}
	}

	// the statistics file must be created and initialized
	stats = (statrecord *)calloc(parent->zonmax+1,sizeof(statrecord));
	if (!stats) parent->speichererror(1);
	for (int i=0;i<=parent->zonmax;i++)
	{ 
		stats[i].n           = stats[i].sx        = 0;
		stats[i].zonenwert   = parent->zone[i].zonenwert;
		stats[i].rel         = parent->zone[i].rel;
		stats[i].RelCellSum  = parent->zone[i].RelCellSum;
	}

	// initializing the statistics class for each register set we use in this module
	int retvalue;
	cStatfilename=new char[_stringsize_];
	sprintf(cStatfilename,"%s",statfilename.c_str());
	zonestat = new zonalstatistic(retvalue,parent->WSstream,parent->LastResults,parent->maxLastResults,parent->RelCellArea,
		parent->nrows,parent->ncols,&stats,(int)parent->zonmax,(int)statwritecode,(int)4,cStatfilename,
		parent->DefaultOutputDirectory,(int)parent->numsubb,(int*)parent->listsubb,(parent->InV < 60));
	if (retvalue!=0) throw retvalue;
	zonestat->schreibzonenwerte(sEntityName.c_str(),parent->instanzname,parent->LastResults,parent->maxLastResults);

	// getting regions and weights
	regio1index=regio2index=regio3index=weight1index=weight2index=weight3index=0;
	WEIGHT=NULL;
	REGION=NULL;
	if (rarray) /* es wurde ein Regionengrid eingelesen und ausgewertet, es gibt also mindestens eine Region */
	{
		WEIGHT=(double***)calloc(4,sizeof(double**));
		REGION=(double***)calloc(4,sizeof(double**));
		REGION[1]=REGION[2]=REGION[3]=WEIGHT[1]=WEIGHT[2]=WEIGHT[3]=NULL;
		regio1index= parent->ist_vorhanden(_regression_regions_1_);
		weight1index=parent->ist_vorhanden(_regression_weights_1_);
		REGION[1]=parent->gridinfo[regio1index].pointer;
		WEIGHT[1]=parent->gridinfo[weight1index].pointer;
		if (regiomax>1)
		{
			regio2index= parent->ist_vorhanden(_regression_regions_2_);
			weight2index=parent->ist_vorhanden(_regression_weights_2_);
			REGION[2]=parent->gridinfo[regio2index].pointer;
			WEIGHT[2]=parent->gridinfo[weight2index].pointer;
		}
		if (regiomax>2)
		{
			regio3index= parent->ist_vorhanden(_regression_regions_3_);
			weight3index=parent->ist_vorhanden(_regression_weights_3_);
			REGION[3]=parent->gridinfo[regio3index].pointer;
			WEIGHT[3]=parent->gridinfo[weight3index].pointer;
		}
	}
	else
	{
		regiomax=0; 
		REGION=WEIGHT=NULL;
	}
}

/* example for the control file structure to parse
temperature {
	entityinputgrid = temperature_reg ;
		regions = 1    2   3;
		weights = 0  0.5 1.0;
	entityinputgrid = temperature_idw ;
		regions = 1    2   3;
		weights = 1.0  0.5 0;
	outputgrid  = $outpath//$tempegrid ; 
		writecode = $Writegrid ;        
	outputtable = $outpath//temp//$grid//.//$code//$year ;
		statcode  = $hour_mean ;
}
*/

/* members of cSuperposition

	struct inputstruct {
		std::string sInputEntity; // name of the input grid as used as name in the section for meteointerpolation
		std::vector<int> region;
		std::vector<double> weight;
		double **INGRID;
	};

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
*/

cSuperposition::~cSuperposition()
{
}

// writing the results will be called after the precipitation correction ran in order to write out the corrected values
void cSuperposition::WriteResultGrids()
{

	parent->writegridtofile(resultindex);
	// statistics in a non-parallel section only
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				zonestat->wertzufuegen(i,j,zarray[(int)zn[i][j]],RESULT[i][j]);
			}
		}
	}
	zonestat->calc_and_done(jahr,monat,tag,stunde, minute);
}

void cSuperposition::flushbuffers()
{
	if (zonestat) zonestat->flushbuffer();
}

void cSuperposition::run(long me,long se,long tg,long mt,long jr)
{
	jahr = jr;
	monat=mt;
	tag=tg;
	stunde=se;
	minute = me;
	zonestat->init(); 
	int incnt=inputs.size();
	int i,j;

#pragma omp parallel for schedule(runtime) shared(incnt) private(i,j)
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				//
				// finally, do the superposition here
				//
				// loop through all inputgrids, get the actual points region and its weight and superpose the results
				RESULT[i][j]=0;
				int r1=REGION[1][i][j];
				for (int is=0;is<incnt;is++)
				{
					double regionweight1=WEIGHT[1][i][j];
					double value=inputs[is]->INGRID[i][j];
					double methodweight=inputs[is]->weight[rarray[r1]-1];
					RESULT[i][j]+=value*methodweight*regionweight1;
				}
				if (regiomax>1)
				{
					int r2=REGION[2][i][j];
					if (r2>0)
					{
						for (int is=0;is<incnt;is++)
						{
							double regionweight2=WEIGHT[2][i][j];
							double value=inputs[is]->INGRID[i][j];
							double methodweight=inputs[is]->weight[rarray[r2]-1];
							RESULT[i][j]+=value*methodweight*regionweight2;
						}
					}
				}
				if (regiomax>2)
				{
					int r3=REGION[3][i][j];
					if (r3>0)
					{
						for (int is=0;is<incnt;is++)
						{
							double regionweight3=WEIGHT[3][i][j];
							double value=inputs[is]->INGRID[i][j];
							double methodweight=inputs[is]->weight[rarray[r3]-1];
							RESULT[i][j]+=value*methodweight*regionweight3;
						}
					}
				}
			}
		}
	}
}

void cSuperposition::ClearMemory()
{
	parent->free_matrix(parent->gridinfo[resultindex].pointer, parent->gridinfo[resultindex].nrl, parent->gridinfo[resultindex].nrh, 1, ncols, nrows);
	parent->delete_grid(resultindex);
	zonestat->ClearMemory();
	if (zonestat) 
	delete(zonestat);
	zonestat=NULL;
	if (stats)    
	free(stats);     
	stats=NULL;
}
