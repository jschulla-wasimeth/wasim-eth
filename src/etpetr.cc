/*   $Revision: 47 $ 
	$Date: 8.06.10 0:04 $  */

/*

(c) 1995-2007 by Joerg Schulla
erste Version: Juni 1995

date       : 02.07.97
version    : 2   Einfuegen der Ausgabe einer Datei mit den korrigierten Strahlungen fuer jede Zone (Statistikdatei)

date       : 19.09.97
version    : 2.5 Neue function "real_uzr" eingefuegt -> soll nur noch die im Bodenmodell ins ETR-Array geschriebene
					reale Verdunstung in eine STatistikdatei schreiben.
Gleichzeitig:    Aenderung der Berechnung der realen Verdunstung: Schneeverdunstung wird getrennt von ::real in ::snow
					berechnet und nach dem Schneemodell aufgerufen,
					Fuer wasim:    reale Verdunstung wird erst nach dem Bodenmodell aufgrund der
								neuen Bodenfeuchte aus dem Bodenmodell berechnet und vom Bodenspeicher abgezogen (in ::real).
					Fuer wasimuzr: reale Verdunstung wird in unsatzon berechnet und in ::real_uzr() nur noch in Statistikdatei
								ausgegeben

date       : 11.03.98
version    : 3   Anpassung an MS-Visual C++: In get_ETP(...) explizites Typecasting beim Aufruf der exp() function erzwungen,
					indem alle Constanten vorher an double-Variablen zugewiesen werden. So kann in MS-VC++ das Flag /Op ausgeschaltet werden
					(improve double conssitency), was den Code extrem schnell gemacht hat (+50%!!!) 

date       : 16.6.1998     
version    : 4   surface resitances (rsc) for Landuse-table are now interpolated daily.
					aerodynamical resistance for vegetation height of > 2m is lowered (25/() instead of 64/()), so it is a continuous function
					using both approaches, the one for vegetation heights below and the one for vegetation heights above 2m

date       : march, 9 1999
version    : 4.1 hour-entry can be given for start or for end of time step via min_to_int_end

date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.

date       : January, 2006
version    : implementing multi layer evaporation model

more version history: see source safe comments

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "teilmod.h"   /* Definition der Teilmodelle */
#include "wasim.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/* oooooooooooooooooooo Beginn der evapotrans_pot - routinen 	oooooooooooo */

int verdunstungsmodellklasse::nullstellen(double a,double b,double *n1,double *n2)
{
	/* liefert in n1 und n2 die Nullstellen der quadratischen Gleichung y=x^2 + ax + b.
Die Werte gelten folgendermassen:
wenn Rueckgabewert eine 0: keine Nullstelle (n1 = n2 =0, beide ungueltig)
						1: doppelte relle Nullstelle (n1 =n2, beide gueltig)
						2: zwei relle Nullstellen (n1 != n2, beide gueltig) */
	double diskriminante;

	diskriminante = a*a/4.0-b;
	if (diskriminante>0)
	{ *n1 = -a/2.0 - sqrt(diskriminante);
		*n2 = -a/2.0 + sqrt(diskriminante);
		return 2;
	}
	else
	if (diskriminante==0.0)
	{ *n1 = *n2 = -a/2.0;
		return 1;
	}
	else
	{ *n1 = *n2 = 0.0;
		return 0;
	}
}

void verdunstungsmodellklasse::initialize_lakes(std::string instanzname)
{
	if (do_it == 0) return;
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

		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}

		strcpy(tstr,instanzname.c_str());

		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		findeabschnitt(tstr);

		fscanfln(steuerdatei,"%d",&lake_active,"flow control parameter for ",tstr);
		
		// switch off lake model when in Spinup-mode
		if (doSpinUp) {
			lake_active=0;
		} 

		if (lake_active != 0)
		{
			fscanfln(steuerdatei,"%d",&lake_ini_method,"initialization method for ",tstr);
			fscanfln(steuerdatei,"%lf",&Albedo_OpenWater,"Albedo_OpenWater for ",tstr);
			//		  fscanfln(steuerdatei,"%lf",&rsc_OpenWater,"rsc_OpenWater for ",tstr);
			rsc_OpenWater=0;
			fscanfln(steuerdatei,"%lf",&z0_OpenWater,"z0_OpenWater for ",tstr);
			//		  fscanfln(steuerdatei,"%lf",&LAI_OpenWater,"LAI_OpenWater for ",tstr);
			LAI_OpenWater=0;
			//		  fscanfln(steuerdatei,"%lf",&VCF_OpenWater,"VCF_OpenWater for ",tstr);
			VCF_OpenWater=0;

			LAKE=gridinfo[getindex(_LAKE_CODES_)].pointer;
			if (ist_vorhanden(_ponding_storage_))
			{
				POND=gridinfo[getindex(_ponding_storage_)].pointer;
			}
			else
			{
				if (ist_vorhanden(_ponding_storage_top_))
				{
					POND=gridinfo[getindex(_ponding_storage_top_)].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n=======================================\n"<<
					"!!!! ERROR !!!! \n\n"<<
					"\nLAKE-grid was found, but no POND-Grid. For using lakes,"<<
					"\nboth grids as well as the MAXPOND grid must be present."<<
					"\n\n=======================================\n"<< "[end-entry]";
					throw(28);
				}
			}

			if ((lake_ini_method<1)||(lake_ini_method>2))
			{
				(*WSstream) << "FATAL: \n=======================================\n"<<
				"!!!! ERROR !!!! \n\n "
				<< "Error in initialization method for "<<instanzname<<".\n"
				<< "Allowed values are 1 (no change of grids) or 2 (adjust DEM, SLOPE, DRAIN).\n"
				<< "In either case, theta_soil will be set to saturated, if the cell is marked as lake and POND is > 0." <<
				"\n\n=======================================\n"<< "[end-entry]";
				throw(211);
			}

		}
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"after reading initialization method for "<<instanzname<<"[end-entry]";
			throw(11);
		}
		delete[] tstr;
	}
}

verdunstungsmodellklasse::verdunstungsmodellklasse(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double breitearg,double laengearg,double bezugsmeridianarg,
	char *bezeichnung,char *steuername,int LogLevel,int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids,bool WriteDoubleGrids,bool WriteNetcdfGrids,char *DefaultInputDirectory,char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow,p_xll, p_yll, p_cellsize, p_nodata, mine, ste, tge, mte, jre, p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode, p_vHRU),
sealev1013(0),c7991(0),c29_33(0),c273_15(0),c0_48(0),c0_065(0),c694_67(0),c0_6589(0),c0_9(0),
mag237_3(0),mag25029(0),mag17_26(0),latheat2500_8(0),psychro1_6158(0),latheat2_372(0),c64_0(0),
c1_0(0),c0_54(0),c20_0(0),c4_72(0),c0_0004(0),c0_004(0),c0_6(0),c0_0067(0),c0_216876(0),c60_0(0),
mcode(0),netto(0),p_T(0),p_E(0),p_RLF(0),p_U(0),p_SSD(0),p_RG(0),p_LUSE(0),p_LDR(0),
Albedo(0),dhm(0),RLF(0),E(0),LDR(0),RG(0),SSD(0),U(0),T(0),Tnacht(0),CL(0),
ETP(0),EIP(0),ETR(0),EVAR(0),EVRS(0),RGEX(0),EVAP(0),RSE(0),RSI(0),LUSE(0),NFK(0),SB(0),SSNOW(0),EXKORR(0),WURZ(0),LAI(0),VBG(0),RSC(0),Z0(0),
nfkratio(0),tagnr(0),PEC(0),fk(0),FH(0),dtmax_akt(0),dt_koeff(0),
Ampl_tag(0),Ampl_nacht(0),Ampl_korr(0),Ampl_Jg(0),minrelRG(0),
zoneRGEX(0),zoneETP(0),zoneETR(0),zoneCL(0),cos_za(0),split(0),extink(0),phi(0),dphi(0),
sigma(0),ldr(0),Nadelwald(0),ssdCnt(0),ssdFlag(0),
EIPstatfilename(0), ETPstatfilename(0), ETRstatfilename(0), EVARstatfilename(0), EVAPstatfilename(0),EVRSstatfilename(0),RGEXstatfilename(0),CLstatfilename(0),
EIPstat(0),ETPstat(0),ETRstat(0),RGEXstat(0), EVARstat(0), EVAPstat(0),EVRSstat(0),CLstat(0),
EIPgridname(0),ETPgridname(0),ETRgridname(0),EVARgridname(0),EVAPgridname(0),EVRSgridname(0),
eipindex(0),etpindex(0),etrindex(0),sbindex(0),evarindex(0),evapindex(0),evrsindex(0),
EIPwritegrid(0),ETPwritegrid(0),ETRwritegrid(0),EVARwritegrid(0),EVAPwritegrid(0),EVRSwritegrid(0),
EIPstatcode(0),ETPstatcode(0),ETRstatcode(0),EVARstatcode(0),EVAPstatcode(0),EVRSstatcode(0),RGEXstatcode(0),CLstatcode(0),
m_iInterceptionIsActive(0),EI(0),
lake_active(0),LAKE(0),POND(0),lake_ini_method(0),rsc_OpenWater(0),z0_OpenWater(0),LAI_OpenWater(0),Albedo_OpenWater(0),VCF_OpenWater(0),ICE(0),
firsttime(true),EVRSGL(0), globalETPHamonFactor(1.0),
do_ArcVeg(0),ARCVEGFRACTIONS(NULL), arcVegNumLayers(0)
{
	try
	{
		retvalue=0;
		int i,j,k;
		double n1,n2;


		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;
		/* reservieren von Speicher fuer die Dateinamen */

#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif

		if (MPIREG==NULL)
		{
			if (ist_vorhanden(_MPIRegions_))
			{
				MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
			}
			else
			{
				MPIREG = NULL;
			}
			if (ist_vorhanden(_MPIHelper_))
			{
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			}
			else
			{
				MPIHLP = NULL;
			}

		}

		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}

		LUSE   = new double**[iNumVegetationLayers+1];
		LAI    = new double**[iNumVegetationLayers+1];
		WURZ   = new double**[iNumVegetationLayers+1];
		VBG    = new double**[iNumVegetationLayers+1];
		RSC    = new double**[iNumVegetationLayers+1];
		RSI    = new double**[iNumVegetationLayers+1];
		Z0     = new double**[iNumVegetationLayers+1];
		ETP    = new double**[iNumVegetationLayers+1];
		EIP    = new double**[iNumVegetationLayers+1];
		ETR    = new double**[iNumVegetationLayers+1];
		

		EIPstatfilename = new char*[iNumVegetationLayers+1];
		ETPstatfilename = new char*[iNumVegetationLayers+1];
		ETRstatfilename = new char*[iNumVegetationLayers+1];
		EIPstat         = new zonalstatistic*[iNumVegetationLayers+1];
		ETPstat         = new zonalstatistic*[iNumVegetationLayers+1];
		ETRstat         = new zonalstatistic*[iNumVegetationLayers+1];
		EIPgridname     = new char*[iNumVegetationLayers+1];
		ETPgridname     = new char*[iNumVegetationLayers+1];
		ETRgridname     = new char*[iNumVegetationLayers+1];
		eipindex        = new int[iNumVegetationLayers+1];
		etpindex        = new int[iNumVegetationLayers+1];
		etrindex        = new int[iNumVegetationLayers+1];

		zoneEIP			= new statrecord*[iNumVegetationLayers+1];
		zoneETP			= new statrecord*[iNumVegetationLayers+1];
		zoneETR			= new statrecord*[iNumVegetationLayers+1];
		
		int lay;
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			ETRgridname[lay]        = newstring(_stringsize_+1,1);
			ETRstatfilename[lay]    = newstring(_stringsize_+1,1);
			ETPgridname[lay]        = newstring(_stringsize_+1,1);
			ETPstatfilename[lay]    = newstring(_stringsize_+1,1);
			EIPgridname[lay]        = newstring(_stringsize_+1,1);
			EIPstatfilename[lay]    = newstring(_stringsize_+1,1);
		}
		EVAPgridname       = newstring(_stringsize_+1,1);
		EVARgridname       = newstring(_stringsize_+1,1);
		EVRSgridname       = newstring(_stringsize_+1,1);
		RGEXgridname       = newstring(_stringsize_+1,1);
		CLgridname         = newstring(_stringsize_+1,1);

		EVARstatfilename   = newstring(_stringsize_+1,1);
		EVAPstatfilename   = newstring(_stringsize_+1,1);
		EVRSstatfilename   = newstring(_stringsize_+1,1);
		RGEXstatfilename   = newstring(_stringsize_+1,1);
		CLstatfilename     = newstring(_stringsize_+1,1);

		/* Parameter einlesen */
		fscanfln(steuerdatei,"%d",&methode          ,"Methode fuer ",tstr);

		if ((methode>1)&&(iNumVegetationLayers>1))
		{
			(*WSstream) << "FATAL: \n---------------------- error in evapotranspiration model ------------------- \n"
			"\nModelling evapotranspiration after HAMON, WENDLING or HAUDE is only possible\n"
			"for unlayered Vegetation. Actual parameters: method: "<<methode<<", number of layers: "<<iNumVegetationLayers<<"\n\n"
			"-----------------------------------------------------------------------\n"<<"[end-entry]";
			throw (44);
		}


		if ((methode>1)&&(InV<1440))
		{
			(*WSstream) << "WARNING: \n---------------------- error in evapotranspiration model ------------------- \n"
			"\nModelling evapotranspiration after HAMON, WENDLING or HAUDE should only be used\n"
			"for daily time steps or larger. The actual time step is less than 1440 min.\n"
			"actual parameters: method: "<<methode<<", timestep: "<<InV<<" minutes\n\n"
			"-----------------------------------------------------------------------\n"<<"[end-entry]";
			// throw(44);
		}
		PEC=(double*)calloc(14,sizeof(double));
		if (fscanf(steuerdatei,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&PEC[1],&PEC[2],&PEC[3],&PEC[4],
					&PEC[5],&PEC[6],&PEC[7],&PEC[8],&PEC[9],&PEC[10],&PEC[11],&PEC[12])!=12)
		{
			(*WSstream) << "FATAL: \n---------------------- error in evapotranspiration model ------------------- \n"
			"there have to be 12 (monthly) values for the PEC-parameter for \n"
			"HAMON-evaporation in the control file.\n"
			"\n-----------------------------------------------------------------------\n"<<"[end-entry]";
			throw(45);
		}
		PEC[0]=PEC[12];
		PEC[13]=PEC[1];

		FH=(double*)calloc(14,sizeof(double));
		if (fscanf(steuerdatei,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&FH[1],&FH[2],&FH[3],&FH[4],
					&FH[5],&FH[6],&FH[7],&FH[8],&FH[9],&FH[10],&FH[11],&FH[12])!=12)
		{
			(*WSstream) << "FATAL: \n---------------------- error in evapotranspiration model ------------------- \n"
			"there have to be 12 (monthly) values for the fh-parameter for \n"
			"HAUDE-evaporation in the control file.\n"
			"\n-----------------------------------------------------------------------\n"<<"[end-entry]";
			throw(45);
		}
		FH[0]=FH[12];
		FH[13]=FH[1];

		fscanfln(steuerdatei,"%lf",&fk               ,"fk - Coastfactor (WENDLING) for",tstr);
		
		fscanfln(steuerdatei,"%s",&(ETPgridname[0][0])    ,"Name of the ETP-Grids for "       ,tstr);
		fscanfln(steuerdatei,"%d",&ETPwritegrid    ,"Writegrid ETP for "            ,tstr);
		fscanfln(steuerdatei,"%s%d",&(ETPstatfilename[0][0]),&ETPstatcode ,"Name of ETP-stat.file for "  ,tstr);
		
		fscanfln(steuerdatei,"%s",&(ETRgridname[0][0])       ,"Name of the ETR-Grids for "       ,tstr);
		fscanfln(steuerdatei,"%d",&ETRwritegrid    ,"Writegrid ETR for    "         ,tstr);
		fscanfln(steuerdatei,"%s%d",&(ETRstatfilename[0][0]),&ETRstatcode ,"Name of the ETR-stat.file for "  ,tstr);
		
		fscanfln(steuerdatei,"%s",EVAPgridname       ,"Name of the evaporation_potential-grid for "       ,tstr);
		fscanfln(steuerdatei,"%d",&EVAPwritegrid    ,"Writegrid EVAP for    "         ,tstr);
		fscanfln(steuerdatei,"%s%d",EVAPstatfilename,&EVAPstatcode ,"Name of the evaporation_potential-stat.files for "  ,tstr);

		fscanfln(steuerdatei,"%s",EVARgridname       ,"Name of the evaporation_real grid for "       ,tstr);
		fscanfln(steuerdatei,"%d",&EVARwritegrid    ,"Writegrid evaporation real for    "         ,tstr);
		fscanfln(steuerdatei,"%s%d",EVARstatfilename,&EVARstatcode ,"Name of the evaporation_real-stat.files for "  ,tstr);
		
		fscanfln(steuerdatei,"%s",EVRSgridname       ,"Name of the snow evaporation grid for "       ,tstr);
		fscanfln(steuerdatei,"%d",&EVRSwritegrid    ,"Writegrid snow evaporation for    "         ,tstr);
		fscanfln(steuerdatei,"%s%d",EVRSstatfilename,&EVRSstatcode ,"Name of the evaporation_real-stat.files for "  ,tstr);
		
		fscanfln(steuerdatei,"%s",&(EIPgridname[0][0])    ,"Name of the EIP-Grids for "       ,tstr);
		fscanfln(steuerdatei,"%d",&EIPwritegrid    ,"Writegrid ETP for "            ,tstr);
		fscanfln(steuerdatei,"%s%d",&(EIPstatfilename[0][0]),&EIPstatcode ,"Name of EIP-stat.file for "  ,tstr);
		

		//	fscanfln(steuerdatei,"%s",RGEXgridname       ,"Name of the corrected radiation grid for "       ,tstr);
		//     fscanfln(steuerdatei,"%d",&RGEXwritegrid    ,"Writegrid corrected radiation for    "         ,tstr);
		sprintf(RGEXgridname,"%s","correctedradiation");
		RGEXwritegrid=0;
		fscanfln(steuerdatei,"%s%d",RGEXstatfilename,&RGEXstatcode ,"Name of the RGEX-stat.file for "  ,tstr);
		
		/* einlesen der Koeffizienten 0. bis 3. Ordung fuer die Berechnung der Strahlung aus der SSD */
		fscanf(steuerdatei,"%lf%lf%lf%lf",&(cr[0]),&(cr[1]),&(cr[2]),&(cr[3]));   satz=fgets(satz,_stringsize_,steuerdatei);
		/* einlesen der Koeffizienten 0. bis 3. Ordung fuer die Berechnung der SSD aus der Globalstrahlung */
		fscanf(steuerdatei,"%lf%lf%lf%lf",&(CS[0]),&(CS[1]),&(CS[2]),&(CS[3]));   satz=fgets(satz,_stringsize_,steuerdatei);
		/* einlesen des Extinktionskoeffizienten und seines Aenderungswertes im Jahresverlauf */
		fscanf(steuerdatei,"%lf%lf",&phi,&dphi);                                satz=fgets(satz,_stringsize_,steuerdatei);
		/* einlesen des Koeffizienten zur Beschr. der Abnahme der Hoehenabhaengigkeit der max. taegl.
	Differenz der Mittelwerte der Tages und der NAchttemperaturen */
		fscanf(steuerdatei,"%lf",&dt_koeff);                                   satz=fgets(satz,_stringsize_,steuerdatei);
		/* einlesen der monatlichen Werte der max. taegl. Differenzen der Tages- und Nachttemperaturmittel (an e gewichtet)  */
		for (i=1;i<=12;i++) fscanf(steuerdatei,"%lf",&(dtmax[i]));             satz=fgets(satz,_stringsize_,steuerdatei);
		dtmax[0]=dtmax[12];
		dtmax[13]=dtmax[1];
		if (fscanf(steuerdatei,"%lf%lf",&Ampl_tag,&Ampl_Jg)!=2)
		{
			(*WSstream) << "FATAL: error reading fraction for temperature amplitude for day/night\n" << 
			"or the annual variability\n"<<"[end-entry]";
			throw(10);
		}
		satz=fgets(satz,_stringsize_,steuerdatei);

		satz=fgets(satz,_stringsize_,steuerdatei);
		sprintf(CLgridname,"%s",_Cloudiness_);
		CLwritegrid=0;
		if (strlen(satz)>2)
		{
			sscanf(satz,"%s%d",CLstatfilename,&CLstatcode);
		}
		else
		{
			sprintf(CLstatfilename,"cloudiness");
			CLstatcode=0;
		}

		satz=fgets(satz, _stringsize_, steuerdatei);
		std::string ts(satz);
		if ((ts.find("HamonCorrectionOpenWater") != std::string::npos)&&(ts.find("=") != std::string::npos))
		{
			std::stringstream ss;
			ts = ts.substr(ts.find("=")+1);
			ss << ts.c_str();
			ss >> globalETPHamonFactor;
		}

		if ((Ampl_tag>=0.0)&&(Ampl_tag<=1.0))
		{    Ampl_nacht = 1.0 - Ampl_tag;  }
		else
		{
			(*WSstream) << "FATAL: \n---------------------------- error in control file --------------------\n\n"
			"The splitting of the temperature amplitude between day and night\n"
			"has to be done in a way that the fraction for the bright day is between 0..1\n"
			"read: amplitude_fraction_day: "<<Ampl_tag<<"\n"
			"\n-----------------------------------------------------------------------\n"<< "[end-entry]";
			throw(10);
		}

		if (cr[0]<0.05)
		{
			(*WSstream) << "FATAL: fraction of diffuse radiation is too small: "<<cr[0]<<" (should be 0.05 ... 0.5)\n"<<"[end-entry]";
			throw(10);
		}
		else
		if (cr[0]>0.5)
		{
			(*WSstream) << "FATAL: fraction of diffuse radiation is too large: "<<cr[0]<<" (should be 0.05 ... 0.5)\n"<<"[end-entry]";
			throw(10);
		}

		if (fabs((cr[0]+cr[1]+cr[2]+cr[3])-1.0)>0.01)
		{
			(*WSstream) << "FATAL: die Koeffizienten fuer die Modellierung der "
			"Globalstrahlung muessen in ihrer Summe 1 ergeben\n"
			<< " derzeitige Summe: "<<cr[0]<<" + "<<cr[1]<<" + "<<cr[2]
			<<" + "<<cr[3]<<" = "<< (cr[0]+cr[1]+cr[2]+cr[3]) <<"[end-entry]";
			throw(10);
		}

		/* Ermitteln von minrelRG, dem Wert der relativen Strahlung, unterhalb dessen SSD=0 ist. Da die SSD durch
	ein Polygon 3. Ordnung aus relRG ermittelt werden kann, wobei aber zwischen 0 und 1 mindestens eine Nullstelle
	(2.Ordnung oder CS[3]==0) oder ein lokales Minimum des Polygons (3.Ordnung, CS[3]!=0) liegen muss, wird
	nun diese Nullstelle bzw. das Minimum (das eine Nullstelle der ersten Ableitung ist) gesucht. Nur relRG-Werte
	ueber dieser Nullstelle werden zur Modelleirung genutzt.
	*/
		if (CS[3]!=0.0)
		{ /* echtes Polynom 3.Ordnung */
			/* 1. Ableitung bilden */
			i=nullstellen(2.0*CS[2]/(3.0*CS[3]),CS[1]/(3.0*CS[3]),&n1,&n2);
			if (i==2)
			{ 
				if ((n2<=1.0)&&(n2>=0.0)) 
				{
					minrelRG=n2; 
					j=1;
				}
				else
				{
					if ((n1 >= 0.0) && (n1 <= 1.0))
					{
						minrelRG = n1;
						j = 1;
					}
					else
					{
						j = 0;
					}
				}
			}
			else
			{
				if (i == 1)
				{
					if ((n1 >= 0.0) && (n1 <= 1.0))
					{
						minrelRG = n1;
						j = 1;
					}
					else
					{
						j = 0;
					}
				}
				else
				{
					j = 0; /* keine Nullstellen zwischen 0 und 1 gefunden */
				}
			}
		}
		else
		{ /* nur Polynom geringerer Ordnung */
			if (CS[2]!=0.0)
			{ /* 2. Ordnung */
				i=nullstellen(CS[1]/CS[2],CS[0]/CS[2],&n1,&n2);
				if (i==2)
				{ 
					if ((n2<=1.0)&&(n2>=0.0)) 
					{
						minrelRG=n2; 
						j=1;
					}
					else 
					{ 
						if ((n1>=0.0)&&(n1<=1.0)) 
						{
							minrelRG=n1; 
							j=1;
						} 
						else
						{
							j = 0;
						}
					}
				}
				else
				{ 
					if (i==1) 
					{ 
						if ((n1>=0.0)&&(n1<=1.0)) 
						{ 
							minrelRG=n1; 
							j=1;
						} 
						else
						{
							j = 0;
						}
					}
					else {  j=0; /* keine Nullstellen zwischen 0 und 1 gefunden */  }
				}
			}
			else
			{
				if (CS[1]!=0)
				{
					minrelRG = -CS[0]/CS[1];
					if ((minrelRG<0.0) || (minrelRG>1.0))
					{
						j = 0;
					}
					else
					{
						j = 1;
					}
				}
				else
				{
					(*WSstream) << "FATAL: The coefficients for sunshine duration modelling are wrong!\n"
					"There must be at minimum a linear relationship\n"
					"between radiation and sunshine duration\n"
					"CS[0]="<<CS[0]<<", CS[1]="<<CS[1]<<", CS[2]="<<CS[2]<<", CS[3]="<<CS[3]<<"\n\n"<<"[end-entry]";
					throw(10);
				}
			}
		}
		if (j==0)
		{
			(*WSstream) << "FATAL: \n----------------------------!! Error !! -----------------------\n"
			"The coefficients for sunshine duration modelling are wrong!\n"
			"No null point (root) was found between 0 and 1!\n"
			"CS[0]="<<CS[0]<<", CS[1]="<<CS[1]<<", CS[2]="<<CS[2]<<", CS[3]="<<CS[3]<<"\n\n"
			"use e.g.: CS[0]=0.072, CS[1]=-0.808, CS[2]=1.749, CS[3]=-0.239\n"
			"----------------------------------------------------------------"<<"[end-entry]";
			throw(10);
		}


		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		// now reading the ArcVeg state (on or off)
		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}

		if (findeabschnitt(_ArcVegIntegration_, false)) // looks for a section  [ArcVegIntegration]
		{
			fscanfln(steuerdatei, "%d", &do_ArcVeg, "", tstr);
			
			if (do_ArcVeg > 0)
			{
				arcVegNumLayers = MultiLayerLUTable[1].getNumFractions();
				if (arcVegNumLayers != iNumVegetationLayers)
				{ // second try: assume the grids may be there but not the initial fractions in the multi layer landuse table
					// check for grids based on num layers
					arcVegNumLayers = iNumVegetationLayers;
				}
				
				bool bAllArcVegFractionsFound = true;
				if (arcVegNumLayers > 0)
				{
					char* layername;
					layername = new char[200];
					arcVegNumLayers = iNumVegetationLayers;
					ARCVEGFRACTIONS = new double** [arcVegNumLayers + 1];

					for (int layer = 1; layer <= arcVegNumLayers; layer++)
					{
						sprintf(layername, "%s%d", "ArcVeg_Fraction_", layer);
						if (ist_vorhanden(layername))
						{
							ARCVEGFRACTIONS[layer] = gridinfo[getindex(layername)].pointer;
						}
						else
						{
							bAllArcVegFractionsFound = false;
						}
					}
					
					delete[] layername;
				}
				else
				{
					bAllArcVegFractionsFound = false;
				}
				if (!bAllArcVegFractionsFound)
				{
					(*WSstream) << "FATAL: ERROR in initializing ETPETR: [ArcVegIntegration] is active but there are no or not all fraction grids found as standard grids (ArcVeg_Fraction_1..n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}
			}
		}

		if (fclose(steuerdatei) == EOF)
		{
			(*WSstream) << "FATAL: error closing " << steuerdateiname << "[end-entry]";
			throw(11);
		}



		/* Methoden zur Verdunstungsberechnung je nach Datenlage (Reihenfolge:
	1,2,3,4,11,12,13,14,21,22,23,24

	Methode  dez. T   e  RLF  U  SSD RG  LUSE LDR
	----------------------------------------------
	11011111  223  x   x   -   x   x   x   x   x
	10111111  191  x   -   x   x   x   x   x   x
	11011110  222  x   x   -   x   x   x   x   -
	10111110  190  x   -   x   x   x   x   x   -
	11010111  215  x   x   -   x   -   x   x   x
	10110111  183  x   -   x   x   -   x   x   x
	11010110  214  x   x   -   x   -   x   x   -
	10110110  182  x   -   x   x   -   x   x   -
	11011011  219  x   x   -   x   x   -   x   x
	10111011  187  x   -   x   x   x   -   x   x
	11011010  218  x   x   -   x   x   -   x   -
	10111010  186  x   -   x   x   x   -   x   -

		gute Methoden 13/14 (nur 4 meteorologische Groessen)
	*/

		p_LDR = 1; p_LUSE = 2; p_RG = 4; p_SSD = 8; p_U = 16; p_RLF = 32; p_E = 64; p_T = 128;

		mcode=0;
		if (methode==1)
		{ /* nur Penman benoetigt diese Groessen  */
			mcode  = p_LDR  * (0!=ist_vorhanden(_Luftdruck_));
			mcode += p_LUSE * (0!=ist_vorhanden(_Landnutzung_"0"));
			mcode += p_U    * (0!=ist_vorhanden(_Windgeschwindigkeit_));
			mcode += p_RLF  * (0!=ist_vorhanden(_Luftfeuchte_));
			mcode += p_E    * (0!=ist_vorhanden(_Dampfdruck_));
		}

		if (methode == 2)
		{ /*  Hamon needs access to land use to define ETPHamonFactor and to calculate EVAP like all other methods --> important for e.g. bare soil, when LAI=0 orVBG=0 */
			mcode += p_LUSE * (0 != ist_vorhanden(_Landnutzung_"0"));
		}

		if (methode == 4)
		{ /*  Haude braucht Feuchte und/oder Dampfdruck von 14 Uhr */
			mcode += p_LUSE * (0 != ist_vorhanden(_Landnutzung_"0"));
			mcode += p_RLF * (0 != ist_vorhanden(_Luftfeuchte_14_));
			mcode += p_E * (0 != ist_vorhanden(_Dampfdruck_14_));
		}

		if ((methode==1)||(methode==3))
		{ /* RG oder SSD werden von PENMAN und WENDLING benoetigt */
			if (ist_vorhanden(_Nettostrahlung_))
			{
				mcode += p_RG   * (0!=ist_vorhanden(_Nettostrahlung_));
				netto=1;
			}
			else
			{
				netto=0;
				mcode += p_RG   * (0!=ist_vorhanden(_Globalstrahlung_));
			}
			mcode += p_SSD  * (0!=ist_vorhanden(_Sonnenscheindauer_));
			mcode += p_LUSE * (0 != ist_vorhanden(_Landnutzung_"0"));
		}


		if ((!(mcode&p_SSD)) && (!(mcode&p_RG)) && (methode!=2) && (methode!=4) )
		/* bei Methode 2 (HAMON) und Methode 4 (HAUDE) wird keine Strahlung benoetigt */
		{
			(*WSstream) << "FATAL: --------------------------- ERROR ----------------------------------------\n"
			<< "Evapotranspiration requires at least one of the following entities\n"
			<< "erlative sunshine duration (SSD, in 1/1)  or global radiation  (RG, in Wh/m^2).\n"
			<< "or net radiation (RN, in Wh/m^2) (either RG or RN!)\n"
			<< "-------------------------------------------------------------------\n"<<"[end-entry]";
			throw(46);
		}

		if ((!(mcode&p_E)) && (!(mcode&p_RLF)) && ((methode==1)||(methode==4)))
		/* bei Methode 1 (PENMAN/MONTEITH) oder 4 (HAUDE) werden e oder rLF benoetigt */
		/* nur bei Methode 1 (PENMAN/MONTEITH) werden e oder rLF benoetigt */
		{
			(*WSstream) << "FATAL: ---------------------------- ERROR ---------------------------------------\n"
			<< "Evapotranspiration requires at least one of the following entities\n"
			<< "vapor pressure (e, in mbar) or relative air humidity (RLF, in 1/1).\n"
			<< "-------------------------------------------------------------------\n"<<"[end-entry]";
			throw(46);
		}

		/* die Temperatur wird von allen Verdunstungsformeln benoetigt */
		if (methode!=4)
		mcode += p_T    * (0!=getindex_first_of_3(_Korrigierte_Temperatur_,_Temperatur_,_Temperatur_));
		else
		mcode += p_T    * (0!=ist_vorhanden(_Temperatur_14_));

		/* zuordnen der Zeiger zu den stellvertretenden Variablen */
		/* T, SB und NFK werden immer benoetigt (von PENMAN/MONTEITH, HAMON und WENDLING) */
		T      = methode!=4? gridinfo[getindex_first_of_3(_Korrigierte_Temperatur_,_Temperatur_,_Temperatur_)].pointer:
		gridinfo[getindex(_Temperatur_14_)].pointer;

		dhm    = gridinfo[getindex(_Hoehenmodell_)].pointer; // dhm is used for radiation modelling, used by Penman (1) and Wendling(3)
		
		/* aber SB und NFK nur, wenn nicht WaSiM-uzr/ETH sondern WaSiM-ETH das aufrufende programm ist:
	deshalb wird die Suche nach diesen beiden Pointern in einer eigenen function durchgefuehrt,
	die nur in WaSiM-ETH aufgerufen wird */
		/*  SB     = gridinfo[getindex(_Bodenspeicher_)].pointer;
	NFK    = gridinfo[getindex(_Nutzbare_Feldkapazitaet_)].pointer; */

		char name[500];

		// LUSE will be used also in snow evaporation in other evaporation methods than Penman
		if (iNumVegetationLayers>1)
		{
			sprintf(name,"%s",_Landnutzung_);
			LUSE[0]   = gridinfo[getindex(name)].pointer;
		}
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;
		}

		
		for (lay = 1; lay <= iNumVegetationLayers; lay++)
		{
			sprintf(name, "%s%d", _Blattflaechenindex_, lay);
			LAI[lay] = gridinfo[getindex(name)].pointer;

			sprintf(name, "%s%d", _Wurzeltiefe_, lay);
			WURZ[lay] = gridinfo[getindex(name)].pointer;

			sprintf(name, "%s%d", _Vegetationsbedeckungsgrad_, lay);
			VBG[lay] = gridinfo[getindex(name)].pointer;

			sprintf(name, "%s%d", _SurfaceCanopyResistance_, lay);
			RSC[lay] = gridinfo[getindex(name)].pointer;

			sprintf(name, "%s%d", _SurfaceInterzepResistance_, lay);
			RSI[lay] = gridinfo[getindex(name)].pointer;

			sprintf(name, "%s%d", _RoughnessLength_, lay);
			Z0[lay] = gridinfo[getindex(name)].pointer;
		}

		if (methode==1) /* dhm, LUSE und U werden nur von PENMAN/MONTEITH benoetigt */
		{
// used by wendling as well --> retrieved earlier			dhm    = gridinfo[getindex(_Hoehenmodell_)].pointer;
			U   = gridinfo[getindex(_Windgeschwindigkeit_)].pointer;
			RSE = gridinfo[getindex(_SurfaceEvaporationResistance_)].pointer;
		}

		/* Albedo wird sowohl von Wendling als auch von Penman/Monteith benoetigt */
		if ((methode==1)||(methode==3)) Albedo = gridinfo[getindex(_Albedo_)].pointer;

		/* bei InV > 720 min. wird die Temperatur aufgrund emprisch gewonnener Beziehungen
	in Tag- und Nachtwerte aufgeteilt. das Temperaturgrid mit den Intervallmittelwerten
	wurde vor der expositionsabhaengigen Temperaturkorrektur in Tnacht dupliziert, so
	dass in T die korrigierten Temperaturen und in Tnacht die unkorrigierten Temperaturen
	zu finden sind. Von Tnacht wird spaeter pro Intervall die Haelfte der fuer die
	entsprechende Hoehe, Jahreszeit und SSD berechneten Temperaturamplitude abgezogen,
	zu T wird dieser Wert addiert. */

		if ((InV>720.0)&&(methode==1))
		/* nur bei PENMAN/MONTEITH wird die 24-h-Verdunstung in Tag und Nacht aufgeteilt */
		Tnacht = gridinfo[getindex(_Nacht_Temperatur_)].pointer;
		else
		Tnacht = T;

		/* von der genutzten Methode abhaengige Zeiger (Datenlage) */
		if (netto)
		{
			if (mcode&p_RG)    RG = gridinfo[getindex(_Nettostrahlung_)].pointer;
		}
		else
		{
			if (mcode&p_RG)    RG = gridinfo[getindex(_Globalstrahlung_)].pointer;
		}
		if (mcode&p_SSD)  SSD = gridinfo[getindex(_Sonnenscheindauer_)].pointer;
		if (mcode&p_E)      E = methode!=4? gridinfo[getindex(_Dampfdruck_)].pointer:gridinfo[getindex(_Dampfdruck_14_)].pointer;
		if (mcode&p_RLF)  RLF = methode!=4? gridinfo[getindex(_Luftfeuchte_)].pointer:gridinfo[getindex(_Luftfeuchte_14_)].pointer;
		if (mcode&p_LDR)  LDR = gridinfo[getindex(_Luftdruck_)].pointer;

		sigma = InV/60*0.0000000567; /* Strahlung eines Koerpers in Wh/m^2 (also schon auf Intervallaenge in h bezogen),
									Stephan Boltzmann Konstante */
		ldr = 1013.0; /* Luftdruckwert in etwa 400..500m Hoehe, wird als Standardwert benutzt, wenn keine PENMAN-Methode
					zur Verdunstungsberechnung genutzt wird, um die theor. Strahlung zu berechnen (Extinktion braucht
					eigentlich den Luftdruck zur Bestimmumg der Dichte der Luft) */

		/* wo ist Platz fuer ETP und ETR im array gridinfo? suchen des ersten freien Zeigers */
		i=1;
		while (gridinfo[i].identification!=NULL)
		{
			i++;
		}

		int ii;
		char layername[2000];
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			if (lay>0) 
			{
				int p;
				//			p=index(ETPgridname[0],"."); 
				std::string s;

				s=std::string(ETPgridname[0]);
				p=s.find_last_of("."); 
				sprintf(ETPgridname[lay],"%.*s_Layer_%d%s",p,ETPgridname[0],lay,&(ETPgridname[0][p]));

				//			p=index(EIPgridname[0],"."); 
				s=std::string(EIPgridname[0]);
				p=s.find_last_of("."); 
				sprintf(EIPgridname[lay],"%.*s_Layer_%d%s",p,EIPgridname[0],lay,&(EIPgridname[0][p]));

				//			p=index(ETRgridname[0],"."); 
				s=std::string(ETRgridname[0]);
				p=s.find_last_of("."); 
				sprintf(ETRgridname[lay],"%.*s_Layer_%d%s",p,ETRgridname[0],lay,&(ETRgridname[0][p]));
				
				//			p=index(ETPstatfilename[0],"."); 
				s=std::string(ETPstatfilename[0]);
				p=s.find_last_of("."); 
				sprintf(ETPstatfilename[lay],"%.*s_Layer_%d%s",p,ETPstatfilename[0],lay,&(ETPstatfilename[0][p]));

				//			p=index(EIPstatfilename[0],"."); 
				s=std::string(EIPstatfilename[0]);
				p=s.find_last_of("."); 
				sprintf(EIPstatfilename[lay],"%.*s_Layer_%d%s",p,EIPstatfilename[0],lay,&(EIPstatfilename[0][p]));
				
				//			p=index(ETRstatfilename[0],"."); 
				s=std::string(ETRstatfilename[0]);
				p=s.find_last_of("."); 
				sprintf(ETRstatfilename[lay],"%.*s_Layer_%d%s",p,ETRstatfilename[0],lay,&(ETRstatfilename[0][p]));

			}
			if (i==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid left"<<gridmaxanzahl<<"[end-entry]";
				throw(12);
			}
			
			sprintf(layername,"%s%d",_Potentielle_Verdunstung_,lay);
			ETP[lay] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(i, ETPgridname[lay], ETP[lay], NULL, 0, layername, ETPwritegrid, mpistartrow, mpiendrow);
			for (k = mpistartrow; k <= mpiendrow; k++)
			{
				for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						ETP[lay][k][j] = 0;
					}
					else
					{
						ETP[lay][k][j] = nodata;
					}
				}
			}
			etpindex[lay] = i;


			zoneETP[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneETP[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneETP[lay][ii].n          = zoneETP[lay][ii].sx         = 0;
				zoneETP[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneETP[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				zoneETP[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}
			
			ETPstat[lay] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&(zoneETP[lay]),(int)zonmax,(int)ETPstatcode,(int)5,(char*)ETPstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb,(InV<60));
			if (retvalue!=0) throw retvalue;

			if (methode==1)
			{ 
				sprintf(layername + strlen(layername), "%s", "(Penman/Monteith)_");
				ETPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults); 
			}
			else
			{
				if (methode==2)
				{ 
					sprintf(layername + strlen(layername), "%s", "(HAMON)_");
					ETPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults); 
				}
				else
				if (methode==3)
				{ 
					sprintf(layername + strlen(layername), "%s", "(WENDLING)_");
					ETPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults); 
				}
				else
				if (methode==4)
				{ 
					sprintf(layername + strlen(layername), "%s", "(HAUDE)_");
					ETPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults); 
				}
				else
				{ 
					sprintf(layername + strlen(layername), "%s", "(unknown_method)_");
					ETPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults); 
				}
			}

			i++;
			if (i==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid left of "<<gridmaxanzahl<<"!\n"<< "[end-entry]";
				throw(12);
			}

			// when running lake model, ETR of layer 0 is required as full grid in rank 0
			int sretr0=(mpi_rank==0)?1:mpistartrow;
			int eretr0=(mpi_rank==0)?nrows:mpiendrow;

			ETR[lay] = matrix(sretr0, eretr0, 1, ncols,nrows);
			sprintf(layername,"%s%d",_Reale_Verdunstung_,lay);
			initialize_grid(i, ETRgridname[lay], ETR[lay], NULL, 0, layername, ETRwritegrid, sretr0, eretr0);
			for (k = sretr0; k <= eretr0; k++)
			{
				for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						ETR[lay][k][j] = 0;
					}
					else
					{
						ETR[lay][k][j] = nodata;
					}
				}
			}
			etrindex[lay] = i;

			zoneETR[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneETR[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneETR[lay][ii].n          = zoneETR[lay][ii].sx         = 0;
				zoneETR[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneETR[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				zoneETR[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}

			ETRstat[lay] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&(zoneETR[lay]),(int)zonmax,(int)ETRstatcode,(int)5,(char*)ETRstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			if (lay==0)
			{
				ETRstat[lay]->schreibzonenwerte("total_Evapotranspiration_as_sum_of_ETR1..n_EI1..n_EVAR_and_EVRS_",instanzname,LastResults,maxLastResults);
			}
			else
			{
				ETRstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults);
			}


			i++;
			if (i==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free ghrid left of "<<gridmaxanzahl<<"!\n"<< "[end-entry]";
				throw(12);
			}
			EIP[lay] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			
			sprintf(layername,"%s%d",_Potential_InterceptionEvap_,lay);
			initialize_grid(i, EIPgridname[lay], EIP[lay], NULL, 0, layername, EIPwritegrid, mpistartrow, mpiendrow);
			for (k = mpistartrow; k <= mpiendrow; k++)
			{
				for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						EIP[lay][k][j] = 0;
					}
					else
					{
						EIP[lay][k][j] = nodata;
					}
				}
			}
			eipindex[lay] = i;

			zoneEIP[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneEIP[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneEIP[lay][ii].n          = zoneEIP[lay][ii].sx         = 0;
				zoneEIP[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneEIP[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				zoneEIP[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}

			EIPstat[lay] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&(zoneEIP[lay]),
			(int)zonmax,(int)EIPstatcode,(int)5,(char*)EIPstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			EIPstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults);
			i++;
		}

		////////////////////////////////////////// now potential evaporation

		EVAP = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, EVAPgridname, EVAP, NULL, 0, _Potential_Evaporation_, EVAPwritegrid, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					EVAP[k][j] = 0;
				}
				else
				{
					EVAP[k][j] = nodata;
				}
			}
		}
		evapindex = i;
		i++;


		////////////////////////////////////////// now real evaporation

		EVAR = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, EVARgridname, EVAR, NULL, 0, _Reale_Evaporation_, EVARwritegrid, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					EVAR[k][j] = 0;
				}
				else
				{
					EVAR[k][j] = nodata;
				}
			}
		}
		evarindex = i;
		i++;

		////////////////////////////////////////// now snow evaporation

		EVRS = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, EVRSgridname, EVRS, NULL, 0, _Snow_Evaporation_, EVRSwritegrid, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					EVRS[k][j] = 0;
				}
				else
				{
					EVRS[k][j] = nodata;
				}
			}
		}
		evrsindex = i;
		i++;

		////////////////////////////////////////// grid with corrected global radiation (using correction factors of radiation correction module)

		RGEX = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, RGEXgridname, RGEX, NULL, 0, _AspectCorrectedGlobalRadiation_, RGEXwritegrid, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					RGEX[k][j] = 0;
				}
				else
				{
					RGEX[k][j] = nodata;
				}
			}
		}
		rgexindex = i;
		i++;

		CL = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, CLgridname, CL, NULL, 0, _Cloudiness_, CLwritegrid, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					CL[k][j] = 0;
				}
				else
				{
					CL[k][j] = nodata;
				}
			}
		}
		clindex = i;
		i++;

		SSDTemp = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(i, "SSDTemp", SSDTemp, NULL, 0, "SSDTemp", 0, mpistartrow, mpiendrow);
		for (k = mpistartrow; k <= mpiendrow; k++)
		{
			for (j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					SSDTemp [k][j] = 0.5;
				}
				else
				{
					SSDTemp [k][j] = nodata;
				}
			}
		}
		ssdtempindex = i;
		i++;


		/* Einrichten eines weiteren Statistikspeichersatzes fuer potentielle evaporation */
		zoneEVAP = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneEVAP) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneEVAP[i].n          = 0;
			zoneEVAP[i].sx         = 0;
			zoneEVAP[i].zonenwert  = zone[i].zonenwert;
			zoneEVAP[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
			zoneEVAP[i].RelCellSum = zone[i].RelCellSum;
		}

		/* Einrichten eines weiteren Statistikspeichersatzes fuer reale evaporation */
		zoneEVAR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneEVAR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneEVAR[i].n          = 0;
			zoneEVAR[i].sx         = 0;
			zoneEVAR[i].zonenwert  = zone[i].zonenwert;
			zoneEVAR[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
			zoneEVAR[i].RelCellSum = zone[i].RelCellSum;
		}

		/* Einrichten eines weiteren Statistikspeichersatzes fuer reale evaporation */
		zoneEVRS = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneEVRS) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneEVRS[i].n          = 0;
			zoneEVRS[i].sx         = 0;
			zoneEVRS[i].zonenwert  = zone[i].zonenwert;
			zoneEVRS[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
			zoneEVRS[i].RelCellSum = zone[i].RelCellSum;
		}

		/* Einrichten eines weiteren Statistikspeichersatzes fuer die korrigierte Globalstrahlung */
		zoneRGEX = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneRGEX) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneRGEX[i].n          = 0;
			zoneRGEX[i].sx         = 0;
			zoneRGEX[i].zonenwert  = zone[i].zonenwert;
			zoneRGEX[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
			zoneRGEX[i].RelCellSum = zone[i].RelCellSum;
		}

		zoneCL = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneCL) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
		    zoneCL[i].n          = 0;
		    zoneCL[i].sx         = 0;
		    zoneCL[i].zonenwert  = zone[i].zonenwert;
		    zoneCL[i].rel        = zone[i].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
		    zoneCL[i].RelCellSum = zone[i].RelCellSum;
		}

		EVAPstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneEVAP,(int)zonmax,(int)EVAPstatcode,(int)5,(char*)EVAPstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		EVAPstat->schreibzonenwerte(_Potential_Evaporation_,"",LastResults,maxLastResults);

		EVARstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneEVAR,(int)zonmax,(int)EVARstatcode,(int)5,(char*)EVARstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		EVARstat->schreibzonenwerte(_Reale_Evaporation_,"",LastResults,maxLastResults);

		EVRSstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneEVRS,(int)zonmax,(int)EVRSstatcode,(int)5,(char*)EVRSstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		EVRSstat->schreibzonenwerte(_Snow_Evaporation_,"",LastResults,maxLastResults);

		RGEXstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneRGEX,(int)zonmax,(int)RGEXstatcode,(int)2,(char*)RGEXstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		RGEXstat->schreibzonenwerte("Global_Radiation_corrected_by_Exposition_",instanzname,LastResults,maxLastResults);

		CLstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&zoneCL,(int)zonmax,(int)CLstatcode,(int)2,(char*)CLstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		CLstat->schreibzonenwerte(_Cloudiness_,instanzname,LastResults,maxLastResults);
    

		/* max. 24 Intevalle innnerhalb eines Tages zur Strahlungsmodellierung */
		cos_za = (double*)calloc(25,sizeof(double));

		for (i=1;i<=LUtable[0].getCode();i++)
		{
			if (strcmp(LUtable[i].getName(),_Nadelwald_)==0)
			{
				Nadelwald=i;
				break;
			}
		}

		/* looking for ICE grid */
		if (ist_vorhanden(_ice_firn_)||ist_vorhanden(_GlacierizedCells_))
		{
			if (ist_vorhanden(_GlacierizedCells_))
				ICE=gridinfo[getindex(_GlacierizedCells_)].pointer;
			else
				ICE=gridinfo[getindex(_ice_firn_)].pointer;
		}
		else
		{
			ICE=NULL;
		}


		/* um fuer dunkle Intervalle die Bewoelkung zur
	Berechnung der langwelligen Strahlungsbilanz abschaetzen zu koennen,
	muss die Sonnenscheindauer ueber alle vorangegangenen hellen Intervalle
	eines Tages gemittelt werden. Dies geschieht pro Teilgebiet, um eine
	gewisse raumliche Variabilitaet zu erlauben. Dazu wird ein Struktur eingerichtet,
	wie sie zur Berechnung der zonalen Statistiken genutzt wird (statrecord).
	Die Zonenwerte werden durch den Mittelwert ersetzt, der fuer eine gamnze Nacht gilt.
	Anstelle der relativen Anteile wird im Feld rel ein Code abgelegt, der darueber
	Auskunft gibt, ob das vorige Intervall hell oder dunkel war (1=hell, 0=dunkel).
	Davon kann dann eine Neuinitialisierung der Register abhaengig gemacht werden,
	denn wenn ein neuer Tag (heller Tag) beginnt, soll die Mittelwertberechnung der
	Sonnenscheindauern ebenfalls neu erfolgen. */
		//    ssdstat=(statrecord*)calloc(zonmax+1,sizeof(statrecord));
		//    for (i=1;i<=zonmax;i++)
		//    {
		//      /* wenn ssd = ssdstat[(int)zn[i][j]].sx/ssdstat[(int)zn[i][j]].n gerechnet
		//         wird, ergibt sich als STandard ssd = 0.5 mit n = 2.0 und sx = 1.0 */
		//      ssdstat[i].n         = 2.0;
		//      ssdstat[i].sx        = 1.0;
		//      ssdstat[i].zonenwert = 0.5;
		//      ssdstat[i].rel       = 0.0;
		//      ssdstat[i].RelCellSum = 2.0;
		//    }

		sealev1013    = 1013.0; /* air pressure at see level */
		c7991         = 7991.0; /* constant for calculating pressure from elevation */
		c29_33        = 29.33;  /* constant for calculating pressure from elevation */
		c273_15       = 273.15; /* absolute temperature at 0 degree celsius */
		c0_48         = 0.48;   /* 1. constant for calulating longwave radiation */
		c0_065        = 0.065;  /* 2. constant for calulating longwave radiation */
		c694_67       = 694.67; /* 1. factor for konverting Wh/m2 to mm/Interval */
		c0_6589       = 0.6589; /* 2. factor for konverting Wh/m2 to mm/Interval */
		c0_9          = 0.90;     /* reduciong net radiation by accounting for soil heat flux */

		mag237_3      = 237.3;   /* factor from the Magnus-formula for saturation vapour pressure  */
		mag25029      = 25029.0; /* derivation factor of magnus formula  */
		mag17_26      = 17.260;  /* factor from Magnus formula */
		latheat2500_8 = 2500.8;  /* latenbt heat of water */
		psychro1_6158 = 1.6158;  /* 1. factor for psychrometric constant */
		latheat2_372  = 2.372;   /* 2. factor for psychrometric constant */

		c64_0         = 25.0;    /* ratio of measurement height and roughness length is set to a constant value for forest */
		c1_0          = 1.0; 
		c0_54         = 0.54;
		c20_0         = 20.0; 
		c4_72         = 4.72; 
		c0_0004       = 0.0004; 
		c0_004        = 0.004; 
		c0_6          = 0.6; 
		c0_0067       = 0.0067; 
		c0_216876     = 0.216876; 
		c60_0         = 60.0; 

		if (doSpinUp!=0 && HRU==NULL)
		{
			HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
			HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
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

void verdunstungsmodellklasse::setSubStepCnt(int step)
{
	substepcnt=step;
}

void verdunstungsmodellklasse::find_SB_and_NFK()
{
	if (do_it == 0) return;
	SB     = gridinfo[getindex(_Bodenspeicher_)].pointer;
	NFK    = gridinfo[getindex(_Nutzbare_Feldkapazitaet_)].pointer; 
}

double verdunstungsmodellklasse::Esatt(double t)
{
	return 6.1078*exp(17.269*t/(237.3+t));
}

/* for Microsoft Visual C++ it is necessary to use the "improve floating point consistency" flag /Op
because of undefined results if not using it (heuristic results). 
Other compilers are not affected by this pragma. (the /Op flag makes the code much faster, 
but it is necessary to carefully test the results compared to versions not using this optimization flag) */
//#pragma optimize ("p",on)
double verdunstungsmodellklasse::getETP(int layer,int i,int j,double tm,double Anteil,int tags,double &evap,double &eip,double &ra_local,double &Rrest_night,double &Rrest_day,double &ssd,double &RK,double &LAIkum)
/* diese Routine berechnet die Verdunstung anhand der Temperatur, es wird anhand dees Parameters tag
entschieden, ob die Widerstaende fuer Tag oder Nacht berechnet werden sollen. Wenn ein Zeitintervall
einen Wert Anteil_hell < 0.95 hat, wird diese Routine zweimal aufgerufen. Besonders bei Tageswerten
ist das erforderlich. */
{
	double a;
	double K_lambert=0.3;
	double LAIscale=10;
	double defizit,ldr,RN,RL,Vv,delta,gamma,dg,rs_evap,rs,ra,rsca,rsi,e;
	double  Ruse_day,Ruse_night;
	double lai_local,vbg_local,rsc_local,z0_local;
	int usePondParameters;


	bool lastlayer=false;
	if (layer < iNumVegetationLayers)
	{
		if (LAI[layer+1][i][j]==nodata) lastlayer=true;
	}

	rs_evap=0.0; // initialize to 0

	vbg_local = VBG[layer][i][j];
	lai_local = LAI[layer][i][j];
	rsc_local = RSC[layer][i][j];
	z0_local  = Z0[layer][i][j];
	usePondParameters = 0;

	// temporary replacement of phenologic parameters by the default parameters for open water in ponds and lakes
	// when the last layer is calculated
	if (POND && LAKE && do_ArcVeg==0) // todo: check if ArcVeg requires normal parameters on open water or the ones for open water as in normal configurations
	{
		if (((iNumVegetationLayers==1)||(layer==iNumVegetationLayers)||lastlayer)&&(LAKE[i][j]>0)&&(POND[i][j]>0))
		{
			if ((LAKE[i][j]>0)&&(POND[i][j]>0))
			{
				vbg_local = VCF_OpenWater;
				lai_local = LAI_OpenWater;
				rsc_local = rsc_OpenWater;
				z0_local  = z0_OpenWater;
				usePondParameters = 1;
			}
		}
	}

	/* Luftdruck errechnen, wenn er nicht interpoliert wird, was meisst unnoetig ist */
	if (!(mcode&p_LDR)) ldr = sealev1013 * exp(-dhm[i][j]/(c7991+c29_33*tm)); else ldr = LDR[i][j];

	//(*WSstream) << "layer\ti\tj\ttm\tAnteil\ttags" << "[end-entry]";
	//(*WSstream) << layer << "\t" << i << "\t" << j << "\t" << tm << "\t" << Anteil << "\t" << tags << "[end-entry]";

	defizit=Esatt(tm); /* um Esatt() nur einmal aufzurufen, wird Esatt in defizit gespeichert */
	/* Dampfdruck errechnen aus RLF, wenn er nicht interpoliert wurde */
	if (!(mcode&p_E))   e   = defizit * RLF[i][j]; else e = E[i][j];

	if ((iNumVegetationLayers>1)&&(do_ArcVeg==0)) // if ArcVeg is active, it should work as if there is only one layer, i.e. the default value is loaded
	{
		K_lambert= MultiLayerLUTable[(int)LUSE[0][i][j]].getK_lambert();
		LAIscale = MultiLayerLUTable[(int)LUSE[0][i][j]].getLAI_scale();
	}

	//(*WSstream) << "K_lambert\tLAIscale" << "[end-entry]";
	//(*WSstream) << K_lambert << "\t" << LAIscale << "[end-entry]";

	if ((layer==1)||(do_ArcVeg>0)) // in arcVeg, the fractions of land use are not stacked but beside each other, so input is the  same for all
	{
		if (!netto)
		{
			a = tm+c273_15;
			RL = a*a*a*a;
			RL *= (sigma*Anteil);
			RL *= ( c0_48 - c0_065*sqrt(e))*(cr[0]+(c1_0-cr[0])*ssd); /* in Wh/m^2 */
			/* Berechnung der Netto-Strahlungsbilanz aus Kurz- und Langwellenstrahlung
		die 0.90 beruecksichtigt den Bodenwaermestrom, Ergebnis in mm!
		wenn das Intervall in einen hellen und in einen dunklen Anteil
		unterteilt ist, wird die Kurzwellige Strahlung nicht nach Anteil_hell/Anteil_dunkel aufgeteilt,
		sondern kommt voll zur Wirkung bei tags==1, bei tags==0 wird keine Kurzwellenstrahlung
		beruecksichtigt. So kann die Verdunstung tags ueber groesser, nachts ueber kleiner gehalten werden. */
			if (tags==0)
			RN = -RL / (c694_67-c0_6589*tm);
			else
			{
				RN = (RK - RL);
				RN = c0_9 * RN;
				RN = RN / (c694_67 - c0_6589*tm);
			}

			/* die 694.67 - 0.6589 ergeben sich aus (2500.8 - 2.372*T) KJ/Kg, der
		Verdunstungswaerme von Wasser nach Umrechnen  in Wh/m^2.
		Da 1 KJ/Kg = 1000 Ws / Kg = 1/3600 * 1000 Wh/Kg = 1/3.6 Wh/Kg
		bedeutet das: aus 2500.8 KJ/Kg werden 2500.8/3.6 = 694.67 Wh/Kg.
		analog werden aus 2.372 KJ/Kg/K eben 0.6589 Wh/Kg/K
		RN ist dann in (Wh/m^2)/(Wh/Kg) = Kg/m^2 = mm angegeben, es bezieht sich also
		auf ein gesamtes Zeitintervall als absolute Verdunstungsgroesse in mm ( noch
		durch den Nenner (delta/gamma +1 +rs/ra) skaliert)
		*/
		}
		else RN=RG[i][j]/(c694_67-c0_6589*tm);
		
		// for multilayer vegetation the effektive net radiation for each layer has to be calculated
		if (tags==0)
		{
			Ruse_night = (do_ArcVeg >= 0) ? RN * vbg_local * (1.0 - exp(-K_lambert * lai_local)) : RN * vbg_local ; // no reduction for ArcVeg != 0, since all land use fractions are beside each other rather than stacked
			Rrest_night = RN-Ruse_night; // this amount of radiation will be available in the lower layers or for evaporation from soil (if the actual layer is the last one)
			if (Rrest_night<0) Rrest_night=0;
		}
		else
		{
			Ruse_day  = (do_ArcVeg >= 0) ? RN*vbg_local*(1.0-exp(-K_lambert*lai_local)): RN * vbg_local; // no reduction for ArcVeg != 0, since all land use fractions are beside each other rather than stacked
			Rrest_day = RN-Ruse_day; // this amount of radiation will be available in the lower layers or for evaporation from soil (if the actual layer is the last one)
			if (Rrest_day<0) Rrest_day=0;
		}
	}
	else
	{
		if (tags==0)
		{
			Ruse_night = Rrest_night * vbg_local*(1.0-exp(-K_lambert*lai_local));
			Rrest_night = Rrest_night - Ruse_night;
			if (Rrest_night<0) Rrest_night=0;
		}
		else
		{
			Ruse_day = Rrest_day * vbg_local*(1.0-exp(-K_lambert*lai_local));
			Rrest_day = Rrest_day - Ruse_day;
			if (Rrest_day<0) Rrest_day=0;
		}
		//	  if (firsttime==true) // if getETP is called multiple times for an interval, then do not change net radiation, except for the first layer (which calculates net radiation)
		//	  {
		//		  Ruse = Rrest * vbg_local*(1.0-exp(-K_lambert*lai_local));
		//		  Rrest = Rrest - Ruse;
		//		  if (Rrest<0) Rrest=0;
		//	  }
	}

	//(*WSstream) << "RK\tRL\tRN\tRuse_day\tRuse_night\tRrest_day\tRrest_night" << "[end-entry]";
	//(*WSstream) << RK << "\t" << RL << "\t" << RN << "\t" << Ruse_day << "\t" << Ruse_night << "\t" << Rrest_day << "\t" << Rrest_night << "[end-entry]";


	/* Berechnung der Steigung der Dampfdruckkurve bei Saettigung */
	a = c1_0/(mag237_3+ tm);
	delta = exp(mag17_26*tm*a)*mag25029*a*a;

	/* Berechnung der Psychrometerkonstante mit Hilfe von Druck und Temperatur als hPa/K */
	gamma = psychro1_6158 * ldr / (latheat2500_8 - latheat2_372*tm);
	defizit -= e;      /* in defizit war vorher der Wert fuer Esatt() gespeichert (Aufruf vor der Berechnung von e */
	if (defizit < 0) defizit=0;

	dg=delta/gamma;

	/* Berechnung des Aerodynamischen Widerstandes */
	if ((layer==1)||(do_ArcVeg>0))
	{
		if (z0_local>2.0)
		ra_local = c64_0/(c1_0+c0_54*U[i][j]);
		else
		{
			a  = log(c20_0/z0_local); /* anstatt z0=0.1*z0 (Rauhigkeitslaenge) -> a= log(2.0/(0.1*z)) */
			ra_local = c4_72*a*a/(c1_0+c0_54*U[i][j]);
		}
		ra=ra_local;
	}
	else
	{
		ra = ra_local + LAIscale*LAIkum; // --> factor 2 should be made a variable given in the control file
		//(*WSstream) << "ra\tra_local\tLAIscale\tLAIKum" << "[end-entry]";
		//(*WSstream) << ra << "\t" << ra_local << "\t" << LAIscale << "\t" << LAIkum << "\t" << "[end-entry]";
	}

	/* Stomatawiderstand */
	//	rsca = LUtable[use].GetRsc(tag,monat,jahr);
	rsca = rsc_local;
	rsi  = RSI[layer][i][j];
	//    rsca = (tag < 16) ? interpolate(-15.5,LUtable[use].rsc[(int)monat-1],15                ,LUtable[use].rsc[(int)monat]  ,tag) :
	//                        interpolate(   16,LUtable[use].rsc[(int)monat]  ,mon[(int)monat]+15,LUtable[use].rsc[(int)monat+1],tag);
	//  rsca = LUtable[use].rsc[(int)monat];


	/* bei Nadelwald rsc korrigieren */
	if (LUtable[(int)LUSE[layer][i][j]].getCode()==Nadelwald)
	{
		/* 1. Temperaturkorrektur */
		/*
	if (tm<= -4.99) rsca=500.0;
	else
	if (tm<20.0) rsca = 25.0*rsca/(tm+5.0);
*/
		/* 2. defizit-Korrektur. Der Stomata-Widerstand steigt von geringen Temperaturen
	gegen 20 mbar steil an (rsc=70: 1.0-0.043*defizit): 5/89 10/123 12/145 15/197 17/260 18/310 19/382 20/500
							(rsc=50: 1.0-0.045*defizit): 5/65 10/91  12/109 15/154 17/213 18/263 19/345 20/500 */
		/*
	if (defizit<=20.0) rsca=rsca/(1.0-0.043*defizit); else rsca=500.0;
*/
	}

	//  bool lastlayer=false;
	if (layer < iNumVegetationLayers)
	{
		if (LAI[layer+1][i][j]==nodata) lastlayer=true;
	}

	rs_evap=0.0; // initialize to 0

	if ((iNumVegetationLayers==1)||(layer==iNumVegetationLayers)||lastlayer||do_ArcVeg>0)
	{ // unstructured vegetation or the bottom layer in structured vegetation or when ArcVeg is active
		/* Berechnen von rs fuer die Nacht (Stomataspalten geschlossen) oder fuer den Tag */
		if ((tags==0)&&(usePondParameters==0))
		{  // Nacht
			//		  rs = c1_0/(lai_local*vbg_local*c0_0004 + c0_004);    /* 1.0/(LAI*Vbg/2500.0 + 1.0/250.0) */
			if ((lai_local>0)&&(vbg_local>0)&&(rsca>0))
			{ // all parameters within valid range
				rs = c1_0/(lai_local*vbg_local*c0_0004);    /* 1.0/(LAI*Vbg/2500.0) -> without evaporation */
			} else if (((lai_local<=0)||(vbg_local<=0))&&(rsca>0))  {
				rs = 2500.0; // no leafs at all, set resitance to maximum value
			} else {
				rs=0.0; // if rsca < = 0 for lake evaporation
			}
		}
		else
		{ // Tag or usePondParameters==1 (when open water is in the POND grid and it is a LAKE)
			if ((lai_local>0)&&(vbg_local>0)&&(rsca>0))
			{
				a=pow(c0_6,(double)(lai_local*vbg_local));
				rs  = c1_0/((c1_0-a)/rsca);   /* rs=1.0/((1.0-a)/rsca);  -> without evaporation */
			} else if (((lai_local<=0)||(vbg_local<=0))&&(rsca>0))  {
				rs=2500;
			} else {
				rs=0.0; // for lakes, rs should be 0 -> rsca is set to 0 for lake cells
			}
		}
		rs_evap = RSE[i][j];
	}
	else
	{	// vertically structured vegetation, all layers above the bottom layer
		if (tags==0)
		{  // Nacht
			if ((lai_local>0)&&(vbg_local>0)&&(rsca>0)) {
				rs = c1_0/(lai_local*vbg_local*c0_0004);    /* 1.0/(LAI*Vbg/2500.0)  --> ohne + 1.0/250.0, weil Evaporation vom Boden nur in unterster Schicht berücksichtigt wird */
			} else if (((lai_local<=0)||(vbg_local<=0))&&(rsca>0))  {
				rs = 2500.0; // no leafs at all, set resitance to maximum value
			} else {
				rs=0.0; // if rsca < = 0 for lake evaporation
			}
		}
		else
		{ // Tag
			//		rs=rsca/(0.5*lai_local); /* r = rsca/(0.5*LAI) */
			if ((lai_local>0)&&(vbg_local>0)&&(rsca>0))
			{
				a=pow(c0_6,(double)(lai_local*vbg_local));
				rs=c1_0/((c1_0-a)/rsca);   /* rs=1.0/((1.0-a)/rsca + a/150.0); */
			} else if (((lai_local<=0)||(vbg_local<=0))&&(rsca>0))  {
				rs=2500;
			} else {
				rs=0.0; // for lakes, rs should be 0 -> rsca is set to 0 for lake cells
			}
		}
		//	  rsi=c1_0/((c1_0-a)/rsi);   /* rs=1.0/((1.0-a)/rsca + a/150.0); */
	}

	//(*WSstream) << "a\tdelta\tgamma\tdefizit\tdg\tra\trsca\trsi\trs\trs_evap" << "[end-entry]";
	//(*WSstream) << a << "\t" << delta << "\t" << gamma << "\t" << defizit << "\t" << dg << "\t" << 
	//ra << "\t" << rsca << "\t" << rsi << "\t" << rs << "\t" << rs_evap << "[end-entry]";

	/* Penman 1948 potentiell
	Vv = InV/60*0.010833*(1+0.54*U[i][j])*defizit;
	ETP[i][j] = (dg * RN + Vv)/(dg + 1.0);
*/

	// Penman monteith real rs >= 0 
	Vv = Anteil*InV*c60_0 * c0_216876 * defizit / (c273_15+tm) / ra;

	/*
	Herleitung dieser Gleichung:

	in E = [(delta/gamma * RN/lambda) + Vv/lambda ] / (delta/gamma + 1 + rs/ra)
	ist Vv/lambda das Ventilationsfeuchteglied in mm pro Zeitintervall.

	Originalgleichung: Vv/lambda = rho*cp/gamma * (es-e)/ra               (1)

	fuer rho (Dichte tr. Luft) wird nach der Gasgleichung eingesetzt:

	rho = 100*p/(Rl*T)                                                  (2)
		mit p     Luftdruck in hPa
			Rl    Gaskonstante fuer Luft = 286.8 KJ/Kg/K
			T     Kelvintemperatur

	gamma = f(p,T) = (cp*p)/(0.622*lamda)                               (3)
		mit cp    spez. Waermekap. der tr. Luft bei konst. Druck = 1.005 KJ/Kg/K
			p     Luftdruck in hPa
			0.622 Verhaeltnis der molaren Massen von Wasserdampf und trockener Luft
			lamda latente Verdunstungswaerme = (2500.8 - 2.372*T) KJ/Kg [T in °C]

	Durch Einsetzen von (2) und (3) in (1) erhaelt man:

	Vv/lambda = (100*p*cp*0.622*lambda)/(Rl*Tkelvin*cp*p*lambda) * (es-e)/ra
				= 62.2/286.8 / T[kelvin] * (es-e)/ra
				= 0.216876 / (273.15 + T) * (es-e)/ra    [Kg/m^2/s]       (4)

	Um die Einheit von Kg/m^2/s in Kg/m^2 zu wandeln, muss das Ergebnis (4) noch mit
	der Anzahl der zum Zeitintervall (mit der relativen Laenge "Anteil")
	gehoerenden Sekunden "InV*60.0" (InV in Minuten) multipliziert werden.
*/
	double Ruse;
	double Rrest;
	if (usePondParameters==0)
	{
		Ruse = (tags==0)?Ruse_night:Ruse_day;
		Rrest = (tags==0)?Rrest_night:Rrest_day;
	}
	else
	{
		Ruse = (tags==0)?Ruse_night+Rrest_night:Ruse_day+Rrest_day;
		Rrest = 0;
	}


	if (((iNumVegetationLayers==1)||(layer==iNumVegetationLayers)||lastlayer||do_ArcVeg>0)&&usePondParameters==0)
	{ // for PONDs in a LAKE (usePondParameters==1), no evaporation must be extracted from the soil
		evap = (dg * Rrest + Vv)/(dg + c1_0 + rs_evap/ra);
	}
	else
	{
		evap=0;
	}

	if (usePondParameters==0)
	{
		eip = (dg * Ruse + Vv)/(dg + c1_0 + rsi/ra); // for the moment, EIP-formula = ETP-formula except differont resistancies
	}
	else
	{ // open water will not have interception evaporation
		eip = 0;
	}

	//(*WSstream) << "Ruse\tRrest\tevap\teip" << "[end-entry]";
	//(*WSstream) << Ruse << "\t" << Rrest << "\t" << evap << "\t" << eip << "[end-entry]";

	return   (dg * Ruse + Vv)/(dg + c1_0 + rs/ra);
}
//#pragma optimize ("p",off)

void verdunstungsmodellklasse::potentiell(long me,long se,long tg,long mt,long jr)
{
	int i,j,k,m;
	double a;
	if (do_it == 0) return;

	minute = me;
	stunde = (int)(se+timekorr);
	tag    = tg;
	monat  = mt;
	jahr   = jr;

	teilmodellklasse::run(me, se, tg, mt, jr, (char *)_ETP_, NULL);
	tagnr=0;

	for (k=1;k<monat;k++) tagnr += mon[k]; tagnr += tag;


	int lay;
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETPstat[lay]->init();
		EIPstat[lay]->init();
	}
	EVAPstat->init();

	RGEXstat->init();
    CLstat->init();

	/* Aufruf der Routine zur Berechnung von Sonnenaufgang (sunr), Sonnenuntergang (suns),
	Dauer der Dunkelheit am Beginn (dark_a) und am Ende des Intervalls (dark_e)
	sowie des relativen Anteils von Helligkeit (Anteil_hell)  und von
	Dunkelheit (Anteil_dunkel) am Intervall */
	get_Tageslicht_Daten();

	/* bei HAMON u. HAUDE wird und darf keine Strahlung genutzt werden, deshalb werden die 
	Strahlungsberechnungen durch ausdrueckliches Angeben von nostrl=1 ausgeschalten */
	if (methode==2 || methode==4) nostrl=1;
	/* wo steht die Sonne zur Mitte des hellen Abschnittes des aktuellen
	Berechnungsintervalles? */
	if ((!nostrl)&&(!((mcode&(p_SSD+p_RG)) == (p_SSD+p_RG))))
	{
		/* die an die Routine uebergebene Zeit muss mit der Korrektur der Daten bezueglich
	der Zeit des Bezugsmeridianes (also der aus der Steuerdatei eingelesenen
	Zeitkorrektur) addiert werden. Da in der Schweiz die ANETZ-Daten in UTC angegeben
	werden, ist die echte Zeit der Meteodatenzeit also um 1h voraus, die Zeitkorrektur
	ist also in diesem Fall 1h (in timekorr) */

		m = (int)ceil(InV/60.0 * Anteil_hell-0.001);
		split = (InV/60.0*Anteil_hell)/ceil(InV/60.0 * Anteil_hell-0.001);
		/* in split ist die Laenge eines Teilintervalles zur Integration der theor. Strahlung
	gemerkt (in h). split soll maximal eine Stunde lang sein. Es wird so gewaehlt, dass
	der helle Anteil eines Intervalles in eine minimale Anzahl gleichlanger,
	jeweils maximal eine Stunde lange Teilintervalle eingeteilt wird. */
		for (k=1;k<=m;k++)
		/* fuer jedes Teilintervall werden nun die Sonnenstaende berechnet, um es nicht
	beim zellweisen Durchlauf fuer jede Zelle tun zu muessen. Die Strahlungsberechnung selbst
	ist von der Hoehe abhaengig und muss deshalb an jedem Punkt durchgefuehrt werden. */
		{
			getsuncoordinates(tag,monat,breite,laenge,
			stunde+hour_offset - InV/60.0 + dark_a + (k-0.5)*split,
			&sunazimuthangle,&sunaltitudeangle);
			sunzenithangle = pi/2 - sunaltitudeangle;
			cos_za[k] = cos(sunzenithangle);
		}
		extink = ((monat-1)*30.5+10.0+tag)*2.0*Pi/365.0;
		extink = phi + dphi*cos(extink);

	}

	if (InV<720.0)
	{
		if (Anteil_hell>0.05)
		{
			if (ssdFlag==0)
			{
				for (i = mpistartrow; i <= mpiendrow; i++)
				{
					for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
					{
						if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
						{
							SSDTemp[i][j]=0;
						}
					}
				}
				// a new day started, so all values for statistics must be set to 0 and the flag must be set to 1 in order to 
				// go over this step in the next interval when it's still bright day
				ssdFlag=1;/* im naechsten Intervall wird daran erkannt, dass es sich
						nicht um das erste helle Intervall des aktuellen Tages handelt. Bei Tageswerten wird
						diese Berechnung nicht benoetigt und auch nicht verwendet, da nostrl immer 0 ist und
						somit die Berechnung der langwelligen Strahlungsbilanz nicht mit der hier berechneten
						SSD erfolgt. */
				ssdCnt=0;
			}
			ssdCnt++;
		}
		else
		{
			// will be reached in the first intervall with Anteil_hell < 0.05, so we will have ssdCnt>0 already and some values in SSDTemp
			if (ssdFlag==1)
			{
				for (i = mpistartrow; i <= mpiendrow; i++)
				{
					for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
					{
						if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
						{
							if (ssdCnt>0)
							SSDTemp[i][j]=SSDTemp[i][j]/ssdCnt;
							else
							SSDTemp[i][j]=0.5;
						}
					}
				}
				/* ssdFlag wird 0 gesetzt, damit im naechsten Intervall erkannt wird, dass die in
			zonenwert abgelegten Werte bereits fuer diese Nacht gueltige, aus den SSD-Werten
			des vorangegangenen Tages berechnete Mittelwerte der SSD sind */
				ssdFlag=0;
			}
		}
	}

	if ((InV>720.0)&&(methode==1))
	{
		dtmax_akt     = (tag < 16) ? interpolate(-15.5,dtmax[(int)monat-1],15  ,dtmax[(int)monat]  ,tag) :
		interpolate(   16,dtmax[(int)monat]  ,mon[(int)monat]+15,dtmax[(int)monat+1],tag);
		Ampl_korr = Ampl_Jg * cos(tagnr/365.0*2*Pi);
	}

	/* ###################### Zellweise Berechnung der Verdunstung ############################ */

	double Isenk,
	alb,
	relRG,     /* aus ssd berechnete zu erwartende relative Strahlung */
	LAIkum,
	Rrest_day,
	Rrest_night; // LAIkum: akkumulated Leaf Area Index, Rrest: remaining radiation for the actual layer; Ruse useable radiation for the actual layer
	double ssd,RK,e;
	double deltaT;    /* Temperaturamplitude zwischen an e gewichteten Tages- und
				Nachtmittel der Temperatur, Ergebnis einer e-Funktion
				f(dhm,tag im Jahr,relRG) */
	double etpTag, etpNacht, evapTag, evapNacht, eipTag, eipNacht;

#pragma omp parallel for schedule(runtime)\
		shared(m) private(i,j,k,Isenk,relRG,ssd,RK,\
		deltaT,LAIkum,lay,Rrest_day,Rrest_night,\
		etpTag,etpNacht,evapTag,evapNacht,eipTag,eipNacht,alb,a,e)
	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank) && !(doSpinUp && HRUMaster[i][j]==nodata))
			//      if (zn[i][j]!=nodata)
			{
				/* Berechnung der max. moeglichen Strahlung auf einer Flaeche senkrecht zur Sonne */
				if (!nostrl)  /* also bei mind. 5% des Intervalls mit Helligkeit */
				{
					/* --------------------- Isenk ------------------ */
					if ((!((mcode&(p_SSD+p_RG))==(p_SSD+p_RG)))&&(cos_za[m/2+1]>0.017)) /* Sonne in der Mitte des Intervalls hoeher als 1 Grad ueber Horizont */
					{
						/* Berechnung der bei Wolkenfreiem Himmel und klarer Atmosphaere auf einen Quadratmeter
			ebene Flaeche auftreffenden Strahlungsenergie [Isenk] = Wh/m^2 */
						Isenk=0.0;
						double ldr=0;
						if (!(mcode&p_LDR)) ldr = sealev1013 * exp(-dhm[i][j]/(c7991+c29_33*T[i][j])); else ldr = LDR[i][j];
						for (k=1;k<=(int)ceil(InV/60.0 * Anteil_hell-0.001);k++)
						Isenk += split*1367.0 * pow((double)extink,(double)((ldr/1013.0)/cos_za[k])) * cos_za[k];
					}
					else Isenk=0; /* Nachts gibt es keine kurzwellige Einstrahlung */

					/* --------------------- SSD ------------------- */
					if (!(mcode&p_SSD))
					{
						if (Isenk>0)
						{ /* wenn SSD modelliert wird, muss RG bekannt sein. Also kann (und muss)
				die relative Strahlung aus der theor. moegl. Strahlung berechnet werden.
				relRG wird spaeter noch fuer die Temperaturkorrektur genutzt */
							relRG = RG[i][j]/Isenk;
							if (relRG<0) 
							{
								relRG=0.0; // no less than 0 allowed since the theoretical incoming shortwave radiation may not be below zero
							} else if (relRG>1) {
								relRG=1.0; // also for the algorithms no more than a relative radiation of 1.0 is allowed in order to limit ssd and temperature change (between day and night) to reasonable limits
							}
							/* minrelRG ist die Nullstelle oder das lok. Minimum des Polynoms, welches
				durch CS beschrieben wird, minrelRG muss zwischen 0 und 1 liegen und wurde im preprocessing bestimmt.
				Die Sonnenscheindauer wird nun aus diesen Koeffizienten und relRG hergeleitet. */
							if (relRG>minrelRG) ssd = CS[0]+relRG*(CS[1]+relRG*(CS[2]+relRG*CS[3])); else ssd=0.0;
						}
						else
						{ ssd   = 0.0; relRG = 0.0; }
						if (ssd<0.0) ssd=0.0; else if (ssd>1.0) ssd=1.0;
					}
					else
					{
						ssd = SSD[i][j];
						/* relRG ist der Anteil an zu erwartender RG anhand der SSD berechnet,
			wird noch fuer T-Korrektur benoetigt und steht deshalb for der if-Abfrage */
						relRG = cr[0]+ssd*(cr[1]+ssd*(cr[2]+ssd*cr[3]));
					}

					/* --------------------  RG  ------------------- */
					if (!(mcode&p_RG))
					{ /* RG wird aus astronomisch moeglicher Strahlung und SSD modelliert und mit expositionskorrektur modelliert */
						if (EXKORR)
						/* es wird nur der direkte Anteil korrigiert, der durch cr[0], dem diffusen Anteil, von der gesamten Strahlung abgetrennt wird */
						RGEX[i][j]  = Isenk*(cr[0]+ssd*(cr[1]+ssd*(cr[2]+ssd*cr[3]))*EXKORR[i][j]);
						else
						RGEX[i][j]  = Isenk*relRG;
					} else
					{ /* gemessenen Globalstrahlung wird korrigiert in Abhaengigkeit von der Sonnenscheindauer */
						if (EXKORR)
						RG[i][j] = RGEX[i][j] = RG[i][j]*(1.0+(1.0-cr[0])*ssd*(EXKORR[i][j]-1));
						else
						RGEX[i][j] = RG[i][j];
					}
					/* uebernehmen der korrigierten Globalstrahlung in ihre Statistikregister */

					/* die modellierte oder interpolierte Sonnenscheindauer wird nun in die Statistikregister uebernommen */
					SSDTemp[i][j]+=ssd; // later, we will calculate the mean ssd for the dark part of the day (i.e. the night)

					/* ist RG groesser, als es theoretisch moeglich ist? Mit 1300 W/m2
			ist eine absolute Obergrenze erreicht die bei einem Extinktions-
			koeffizienten von 0.9 und senkrecht stehender Sonne in 5000m
			erreicht werden wuerde. In gemaessigten Breiten hier:48 Grad, kann
			am 21.6. im Hochgebirge auf 3000m bei sehr klarer Luft die
			maximale Strahlung noch ca. 1260 W/m2 betragen */

					if (RGEX[i][j] > 1260.0*InV/60.0*Anteil_hell) 
					RGEX[i][j] = 1260.0*InV/60.0*Anteil_hell;

					/* Kurzwellige Strahlungsbilanz */
					if (POND && LAKE)
						alb = ((LAKE[i][j]>0)&&(POND[i][j]>0))?Albedo_OpenWater:Albedo[i][j];
					else
						alb = Albedo[i][j];

					RK = (1.0 - alb)*RGEX[i][j]; /* in Wh/m^2 */
				    if (relRG > 1.0) 
					{
						CL[i][j]=0.0;
					}
					else 
					{
						CL[i][j]= 0.997+0.4473*relRG-1.4059*relRG*relRG;
					}

					if (CL[i][j]>1.0) 
					{
						CL[i][j]= 1.0;
					}

				}
				else
				{
					/* Nachts gibt es keine Kurzwellige Strahlung (es ist ja dunkel!)
			die langwellige Strahlung dagegen kann durchaus bis zu -100 W/m^2
			erreichen. Dadurch kann die Verdunstung reduziert werden. Erst
			die Summe (delta/gamma*RN+Vv) muss positiv sein. Negative Werte
			sollten 0 gesetzt werden */
					RK = RGEX[i][j] = relRG = 0.0;  /* relRG wird zur T-Korrektur verwendet. Da diese nur fuer InV>720 (also 1440) stattfindet und
									so nostrl immer 0 ist (es sei denn, es ist Polarnacht!), wird das hier berechnete relRG
									nie beruecksichtigt. Es steht nur sicherheitshalber hier. */
					ssd = SSDTemp[i][j];
					//		  ssd = ssdstat[zarray[(int)zn[i][j]]].zonenwert;
					CL[i][j]=0.5; // default cloudiness by night is 0.7
				}

				//		RGEXstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],RGEX[i][j]);

				EVAP[i][j]=0.0;   // initializing evaporation to 0
				EVAR[i][j]=0.0;   // and do this for real evaporation, too --> only Penman/Monteith will calculate evaporation seperately, so it's 0 for all other methods
				// but it will be set when calculating the real evaporation for the other methods (as remainig amount of water to be extracted from soil if 
				// not all energy reserved for transpiration could be consumed by the plants
				EVRS[i][j]=0.0;
				/* ---------------------- Temperaturkorrektur bei InV > 720.0 auf Tag- und Nachtwert ---------------- */
				if (methode==1)
				{
					if (InV>720.0)
					{
						deltaT        = dtmax_akt*relRG*exp(-dhm[i][j]/dt_koeff);
						T[i][j]      += (Ampl_tag  +Ampl_korr)*deltaT;
						Tnacht[i][j] -= (Ampl_nacht-Ampl_korr)*deltaT;
					}

					ETP[0][i][j]=0.0; // initializing transpiration to 0
					EIP[0][i][j]=0.0; // initializing transpiration to 0
					ETR[0][i][j]=0.0;


					double ra_local=0;
					LAIkum = 0;
					Rrest_day=Rrest_night=0; // will be changed in getETP!

					for (lay=1;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]!=nodata)
						{
							if (lay>1) LAIkum+=LAI[lay-1][i][j]; // only the accumulated LAI ABOVE the current layer should be used here
							//				  firsttime=true;
							if (Anteil_dunkel >= 0.01) 
							{
								etpNacht   = getETP(lay,i,j,(double)Tnacht[i][j],(double)Anteil_dunkel,0,evapNacht,eipNacht,ra_local,Rrest_night,Rrest_day,ssd,RK,LAIkum); 
								//					(*WSstream) << wasim::setw(20) << "Nacht, Layer " << lay << ": " << etpNacht <<"[end-entry]";

								//	                firsttime=false;
							}
							else 
							{
								etpNacht  = 0.0;
								evapNacht = 0.0;
								eipNacht  = 0.0;
							}
							
							if (Anteil_hell   >= 0.01)
							{
								etpTag     = getETP(lay,i,j,(double)T[i][j]     ,(double)Anteil_hell  ,1,evapTag,eipTag,ra_local,Rrest_night,Rrest_day,ssd,RK,LAIkum); 
								//					(*WSstream) << wasim::setw(20) << "Tag, Layer " << lay << ": " << etpTag <<"[end-entry]";
								//	                firsttime=false;
							}
							else 
							{
								etpTag   = 0.0;
								eipTag   = 0.0;
								evapTag  = 0.0;
							}

							if (do_ArcVeg > 0)
							{
								ETP[lay][i][j] = ARCVEGFRACTIONS[lay][i][j] * (etpTag + etpNacht); // this is the potential transpiration only
								EIP[lay][i][j] = ARCVEGFRACTIONS[lay][i][j] * (eipTag + eipNacht); // this is the potential transpiration only
								EVAP[i][j]    += ARCVEGFRACTIONS[lay][i][j] * (evapTag + evapNacht);  // and this is the potential evaporation (from first soil layer)
							}
							else
							{
								ETP[lay][i][j] = etpTag + etpNacht; // this is the potential transpiration only
								EIP[lay][i][j] = eipTag + eipNacht; // this is the potential interception evaporation only
								EVAP[i][j] += evapTag + evapNacht;  // this is the potential evaporation from first soil layer
							}

							if (ETP[lay][i][j]<0.0) ETP[lay][i][j]=0.0;
							if (EIP[lay][i][j]<0.0) EIP[lay][i][j]=0.0;
							ETR[lay][i][j]=0.0; /* Sicherheitshalber, damit nicht die alte Verdunstung noch irgendwo genutzt wird */
							/* Entnehmen der verdunsteten Wassermenge erfolgt in verdunstungsmodell::real, da die
					Interzeption noch ETR reduzieren Kann */

							// for balance, the total potential evaporation should be written out
							// this must be done as a simple addition for layered vegetation and for ArcVeg type cacluation, since the raw ETP[lay][i][j] values were already fraction-weighed in getETP
							ETP[0][i][j] += ETP[lay][i][j];
							EIP[0][i][j] += EIP[lay][i][j];
							
						}
						else
						{
							ETP[lay][i][j]=0;
							EIP[lay][i][j]=0;
							//				  ETPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0); // to weigh all values regarding the total number of valid cells, also a nodata
							//				  EIPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0); // to weigh all values regarding the total number of valid cells, also a nodata
							// cell for a specific layer has to be added to the statistics
						}
					}
					ETP[0][i][j] += EVAP[i][j]; // since EVAP is incremented for each layer, it must be added only once to ETP, which is after all layers are processed (was before in the layer loop, leading to too high increments and possibly to a larger ETP than the sum of all layers ETP plus EVAP)
				}
				else // other methods will only give results for the uppermost layer -> 
				// should be made sure during startup that if another method for 
				// evaporation than penmen is used, no multiple layers are allowed
				{
					
					if (methode==2)
					{
						/* Verdunstung nach HAMON */
						if (tag<16)
						{
							a=interpolate(-15,PEC[(int)monat-1],15,PEC[(int)monat],tag);
						}
						else
						{
							a=interpolate(16,PEC[(int)monat],16+mon[(int)monat],PEC[(int)monat+1],tag);
						}
						ETP[0][i][j]  = 0.1651 * a * 2.0 * Anteil_hell * 216.7 * Esatt(T[i][j]);
						ETP[0][i][j] /= (T[i][j] + 273.3);
						ETP[0][i][j] *= InV / 1440.0;
						double f = LUtable[(int)LUSE[1][i][j]].getETPHamonFactor();
						if (f != 1.0)
						{
							ETP[0][i][j] *= f; // optional correction factor per land use. Important e.g. for open watre, i.e. lake evaporation
						}
						else
						{
							// check if the actual POND is > 5mm or so and set a separate correction factor valid for all land use types (read from evapotranspiration section)
							if (lake_active && (POND[i][j] > 5.0 && LAKE[i][j] != nodata))
							{
								ETP[0][i][j] *= globalETPHamonFactor;
							}
						}
						EVAP[i][j] = (1 - VBG[1][i][j]) * ETP[0][i][j];
					}
					else
					{
						if (methode==3)
						{
							/* Verdunstung nach WENDLNG */
							ETP[0][i][j]  = (RGEX[i][j]*0.36*(1.1-alb)) + 93.0*fk;
							ETP[0][i][j] *= (T[i][j]+22.0);
							ETP[0][i][j] /= (150.0*(T[i][j]+123.0));
							ETP[0][i][j] *= InV / 1440.0;
							EVAP[i][j] = (1 - VBG[1][i][j]) * ETP[0][i][j];
						}
						else
						{
							if (methode==4)
							{
								/* Verdunstung nach HAUDE */
								/* Berechnung es=f(T) */
								if (T[i][j]<0) 
								{
									e=6.11*exp(22.46*T[i][j]/(272.46+T[i][j]));
								}
								else
								e=6.11*exp(17.62*T[i][j]/(243.12+T[i][j])); 

								if (!(mcode&p_RLF))
								{
									/* if no rel. humidity but vapour pressure was read in, the rel. hum. 
					can be calculated from e (temperature dependend saturation pressure)
					and from the observed vapour pressure
				*/
									a=E[i][j]/e;
									if (a<0) a=0;
									if (a>1) a=1;
								}
								else
								{
									a=RLF[i][j];
								}

								ETP[0][i][j] = FH[int(monat)]*e*(1-a);
								if (ETP[0][i][j]>7) { ETP[0][i][j]=7; } // ETP muss < 7 mm/d sein.
								ETP[0][i][j] *= InV / 1440.0;
								EVAP[i][j] = (1 - VBG[1][i][j]) * ETP[0][i][j];
								/*    (*WSstream) << "Monat=%i, T=%3.2f, es=%10.4f, fh=%3.2f, RLF=%10.4f, ETP=%10.4f\n",
	monat,T[i][j],a,FH[int(monat)-1],RLF[i][j],ETP[i][j]); */  //Ausgabekontrolle
							}
						}
					}
					if (ETP[0][i][j]<0.0) ETP[0][i][j]=0.0;
					EIP[0][i][j]=ETP[0][i][j];
					ETR[0][i][j]=0.0; /* Sicherheitshalber, damit nicht die alte Verdunstung noch irgendwo genutzt wird */
					ETR[1][i][j]=0.0;
					// Entnehmen der verdunsteten Wassermenge erfolgt in verdunstungsmodell::real, da die
					// Interzeption noch ETR reduzieren Kann 
					// weil dort aber ETP[1..n] genutzt werden, muss hier ETP[1]=ETP[0] gesetzt werden.
					ETP[1][i][j]=ETP[0][i][j];
					EIP[1][i][j]=EIP[0][i][j];
				}
				//        ETPstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETP[0][i][j]);
				//        EIPstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EIP[0][i][j]);
				//        EVAPstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVAP[i][j]);
			}
			else 
			{
				for (lay=0;lay<=iNumVegetationLayers;lay++)
				{
					ETP[lay][i][j] = nodata;
					EIP[lay][i][j] = nodata;
					EVAP[i][j] = nodata;
				}
			}
		}
	}
	// filling statistics is stripped out of the parallel section because of commonly used variables in the called statistics functions
	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//      if (zn[i][j]!=nodata)
			{  
				// when in SpinUp mode, all cells that were not computed must be filled with the HRUMaster data now
				if (doSpinUp!=0 && HRUMaster[i][j]==nodata)
				{ 
					// getting the coordinates of the master cell
					int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
					int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
					// filling the grid cells
					for (lay=0;lay<=iNumVegetationLayers;lay++)
					{
						ETP[lay][i][j] = ETP[lay][imaster][jmaster];
						EIP[lay][i][j] = EIP[lay][imaster][jmaster];
					}
					EVAP[i][j] = EVAP[imaster][jmaster];
					RGEX[i][j] = RGEX[imaster][jmaster];
					CL[i][j] = CL[imaster][jmaster];
				}

				if (methode==1)
				{
					for (lay=1;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]!=nodata)
						{
							ETPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETP[lay][i][j]);
							EIPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EIP[lay][i][j]);
						}
						else
						{
							ETPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0); // to weigh all values regarding the total number of valid cells, also a nodata
							EIPstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0); // to weigh all values regarding the total number of valid cells, also a nodata
						}
					}
				}
				ETPstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETP[0][i][j]);
				EIPstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EIP[0][i][j]);
				EVAPstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVAP[i][j]);
				RGEXstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],RGEX[i][j]);
				CLstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],CL[i][j]);
			}
		}
	}


	stunde=se;
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETPstat[lay]->calc_and_done(jahr,monat,tag,stunde,minute);
		EIPstat[lay]->calc_and_done(jahr,monat,tag,stunde,minute);
		writegridtofile(etpindex[lay]);
		writegridtofile(eipindex[lay]);
	}
	EVAPstat->calc_and_done(jahr,monat,tag,stunde, minute);
	writegridtofile(evapindex);
	RGEXstat->calc_and_done(jahr,monat,tag,stunde, minute);
    CLstat->calc_and_done(jahr,monat,tag,stunde, minute);

}


void verdunstungsmodellklasse::SetExpoKorrAdress(double **adresse)
{
	if (do_it == 0) return;
	EXKORR = adresse;
	return;
}


void verdunstungsmodellklasse::SetSnowStorAdress(double **adresse)
{
	if (do_it == 0) return;
	/* Da zum Zeitpunkt des Erzeugens und Initialisierens des Verdunstungsmodells
das Schneemodell noch nicht existiert, wird diese Abfrage der
Adresse des Schneespeichers separat ausgefuehrt (erst, wenn schneemodel existiert) */
	SSNOW = adresse;
	return;
}

void verdunstungsmodellklasse::snow(long me,long se,long tg,long mt,long jr)
{
	int i,j;
	double etrs;
	if (do_it == 0) return;
	if (!SSNOW) return;

	/* no change in tracer concentrations occurs by evaporation (non-evaporating tracers are handled as masses, 
	they remain constant, too!) */
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	teilmodellklasse::run(me, se, tg, mt, jr, (char *)_Schneeverdunstung_ , NULL);

	if (firsttime)
	{
		if (ist_vorhanden(_SnowEvaporation_onGlac_)&&ist_vorhanden(_TotalSnowStorage_onGlacier_))
		{
			EVRSGL = gridinfo[getindex(_SnowEvaporation_onGlac_)].pointer;
			SSNOWGL = gridinfo[getindex(_SolidSnowStorage_onGlacier_)].pointer;
		}
		else
		{
			EVRSGL=NULL;
			SSNOWGL=NULL;
		}
		firsttime=false;
	}

	EVRSstat->init();

	double *ETPtemp,EVAPtemp;
	ETPtemp = new double[iNumVegetationLayers+1];

//#pragma omp parallel for schedule(runtime) private(i,j,etrs)
	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank) && !(doSpinUp && HRUMaster[i][j]==nodata))
			//      if (zn[i][j]!=nodata)
			{
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

				/* Wenn Das Schneemodell auch gerechnet wird, dann soll,
		wenn mehr als 0 mm Schnee liegt, ETR zunaechst vom Schnee und erst dannach
		vom Boden entnommen werden */
				// bei Nutzung der neuen energybalance Methoden aber wird EVRS aus Schneemodell übernommen
				if (this->iSnowModelMethod==5 || this->iSnowModelMethod==6 || this->iSnowModelMethod==7 || this->iSnowModelMethod==11)
				{
					// EVRS was already set in snowmodel. If there is still snow, there is no energy left for evaporation (EVAR) or transpiration (ETR)
					// so we need to set EVAP and ETP to 0 for all layers. Only when there is no snow cover, ETP and EVAP may be used for extraction of water from the soil. 
					// But this is not handled here, since SNOW[i][j] must be > 0 anyway
					if (SSNOW[i][j] > 0.0) /* liegen mehr als 0 mmWE Schnee ? */
					{
						for (int lay=1;lay<=iNumVegetationLayers;lay++)
						{
							if (LUSE[lay][i][j]>nodata)
							{
								ETP[lay][i][j]=0;
							}
							EVAP[i][j]=0;
						}
					} // otherwise let EVAR and ETR calculation use the potential values of EVAP and ETP[lay]
				}
				else
				{
					EVRS[i][j]=0.0; // if there is no snow, we have no snow evaporation. Otherwise, EVRS will be filled now
					if (EVRSGL) EVRSGL[i][j]=0.0;

					if (icevalue<1.0) // the unglacierized fraction of a cell is handled as a normal landuse
					{
						for (int lay=1;lay<=iNumVegetationLayers;lay++)
						{
							ETPtemp[lay]=ETP[lay][i][j];
							EVAPtemp=EVAP[i][j];
						}
						if (SSNOW[i][j] > 0.0) /* liegen mehr als 0 mmWE Schnee ? */
						{
							for (int lay=1;lay<=iNumVegetationLayers;lay++)
							{
								if (LUSE[lay][i][j]>nodata)
								{
									if (SSNOW[i][j]<ETP[lay][i][j])
									{
										etrs=SSNOW[i][j]; /* aller Schnee verdunstet und wird zu ETR */
										EVRS[i][j]+=etrs; 
										SSNOW[i][j]=0;         /* Schneespeicher ist leer */
										ETP[lay][i][j]-=etrs;  /* fuer ETP bleibt noch der Rest ETP-ETR zur Entnahme aus dem Boden */
									}
									else /* es bleibt noch Schnee uebrig */
									{
										etrs=ETP[lay][i][j];
										EVRS[i][j]+=etrs;
										SSNOW[i][j]-=etrs;
										ETP[lay][i][j]=0;
									}
								}
							}
							// and now the same for EVAP, because in case of snow there will be now evapporation from soil
							if (SSNOW[i][j]<EVAP[i][j])
							{
								etrs = SSNOW[i][j];
								EVRS[i][j]+=etrs; /* aller Schnee verdunstet und wird zu EVRS */
								SSNOW[i][j]=0;         /* Schneespeicher ist leer */
								EVAP[i][j]-=etrs;  /* fuer EVAP bleibt noch der Rest EVAP-EVAR zur Entnahme aus dem Boden */
							}
							else /* es bleibt noch Schnee uebrig */
							{
								etrs = EVAP[i][j];
								EVRS[i][j]+=etrs;
								SSNOW[i][j]-=etrs;
								EVAP[i][j]=0;
							}
						}
					}

					if (icevalue > 0.0) // the glacierized fraction of a cell can always evaporate with potential rate, the result must be stored in EVRSGL
					{
						for (int lay = 1; lay <= iNumVegetationLayers; lay++)
						{
							ETPtemp[lay] = ETP[lay][i][j];
							EVAPtemp = EVAP[i][j];
						}
						for (int lay = 1; lay <= iNumVegetationLayers; lay++)
						{
							if (LUSE[lay][i][j]>nodata)
							{
								if (SSNOWGL[i][j]<ETPtemp[lay])
								{
									etrs=SSNOWGL[i][j]; /* aller Schnee verdunstet und wird zu ETR */
									EVRSGL[i][j]+=etrs; 
									SSNOWGL[i][j]=0;         /* Schneespeicher ist leer */
									ETPtemp[lay]-=etrs;  
								}
								else /* es bleibt noch Schnee uebrig */
								{
									etrs=ETPtemp[lay];
									EVRSGL[i][j]+=etrs;
									SSNOWGL[i][j]-=etrs;
									ETPtemp[lay]=0;
								}
							}
						}
						// and now the same for EVAP, because in case of snow there will be now evapporation from soil
						if (SSNOWGL[i][j]<EVAPtemp)
						{
							etrs = SSNOWGL[i][j];
							EVRSGL[i][j]+=etrs; /* aller Schnee verdunstet und wird zu EVRS */
							SSNOWGL[i][j]=0;         /* Schneespeicher ist leer */
							EVAPtemp-=etrs;  /* fuer EVAP bleibt noch der Rest EVAP-EVAR zur Entnahme aus dem Boden */
						}
						else /* es bleibt noch Schnee uebrig */
						{
							etrs = EVAPtemp;
							EVRSGL[i][j]+=etrs;
							SSNOWGL[i][j]-=etrs;
							EVAPtemp=0;
						}
					}
				}
				EVRSstat->wertzufuegen(i, j, zarray[(int)zn[i][j]], EVRS[i][j]);
			}
		}
	}

	delete[] ETPtemp;

	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//      if (zn[i][j]!=nodata)
			{
				// when in SpinUp mode, all cells that were not computed must be filled with the HRUMaster data now
				if (doSpinUp!=0 && HRUMaster[i][j]==nodata)
				{ 
					// getting the coordinates of the master cell
					int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
					int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
					// filling the grid cells
					for (int lay=0;lay<=iNumVegetationLayers;lay++)
					{
						ETP[lay][i][j] = ETP[lay][imaster][jmaster];
						EIP[lay][i][j] = EIP[lay][imaster][jmaster];
					}
					EVAP[i][j]    = EVAP[imaster][jmaster];
					EVRS[i][j]    = EVRS[imaster][jmaster];
					if (EVRSGL) 
					{
						EVRSGL[i][j]  = EVRSGL[imaster][jmaster];
					}
					SSNOW[i][j]   = SSNOW[imaster][jmaster];
					if (SSNOWGL) 
					{
						SSNOWGL[i][j] = SSNOWGL[imaster][jmaster];
					}
				}

				// no weighing for glacierized and unglazierized parts of a cell is done here because 
				// snow evaporation will be taken from every cell (even from Firn or Ice, if no snow is there, 
				// which may bring a small numerical error into the snow balance but not into the total balance 
				// (glacier balance will then correct the snow balance)
				EVRSstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVRS[i][j]);
			}
		}
	}
	EVRSstat->calc_and_done(jahr,monat,tag,stunde, minute);
	writegridtofile(evrsindex);
}

void verdunstungsmodellklasse::real(long me, long se,long tg,long mt,long jr)
/* wird in der alten Version wasim aufgerufen (von maoappli aus). In der neuen Version erfolgt die
Berechnung der realen Verdunstung im Bodenmodell und im Verdunstungsmodell wird nur noch die
Statistikdatei geschrieben -> function verdunstungsmodell::real_uzr() */
{
	int i,j;
	double ETR_temp,EVAR_temp;
	if (do_it == 0) return;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;

	teilmodellklasse::run(me, se, tg, mt, jr, (char *)_ETR_, NULL); 

	// if this is the first time, the evaporation and transpiration components are summarized, look for interception evaporation
	if (m_iInterceptionIsActive==0)
	{
		char layername[2000];
		sprintf(layername,"%s%d",_Interzeptionsverdunstung_,0);
		int pos=ist_vorhanden(layername);
		if (pos>0)
		{
			EI=gridinfo[getindex(layername)].pointer;
			m_iInterceptionIsActive=1;
		}
	}

	int lay;
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETRstat[lay]->init();
	}
	EVARstat->init();

	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//      if (zn[i][j]!=nodata)
			{
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

				if (icevalue<1.0)
				{
					/* in ETP steht nun alles, was noch verdunstet werden kann
						(und zwar ist das die in ::potentiell berechnete Menge ETP
						abzueglich Interzeptions- und Schneeverdunstung */

					/* Beruecksichtigung der Fuellung des Bodenspeichers. Die Landnutzung "Wasser" (mit nFK==0) fuehrt zu ETP=ETR,
						da nfkratio nun bei nfk==0 zu 1.0 gesetzt wird. */
					if (NFK[i][j]>0) nfkratio=SB[i][j]; else nfkratio = 1.0;

					/*      Reduzieren der Widerstaende
						if (nfkratio < 0.6)
						if (nfkratio>=0.2)
						rs = (nfkratio-0.2)/0.4*rs + (0.6-nfkratio)/0.4*1000.0;
						else rs=1000.0;
						ETR[i][j] = (dg * RN + Vv)/(dg + 1.0 + rs/ra);
					*/
					/* Reduzieren der potentiellen Verdunstung */
					ETR[0][i][j]=0;
					for (lay=1;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]>nodata)
						{
							if ((lake_active==1)&&(LAKE[i][j]!=nodata)&&(POND[i][j]>1.0))
							{   // if this is a lake cell, then ETR should be equal to potential ETP
								// and ETR should be taken from POND instead from soil

								if (ETP[lay][i][j]>=POND[i][j])  /* more ETP than is available in pond? (both storages are given in mm) */
								{
									/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
									ETR_temp = POND[i][j]; /* adding the POND-content to the actual evaporation (there may be already some actual evap., because of snowevaporation) */
									ETP[lay][i][j] -= POND[i][j]; /* remaining evaporation is reduced by the POND-evaporation */
									POND[i][j] = 0.0;        /* POND evaporates complete */
								}
								else                       /* in POND is more water than necessary for evaporation */
								{
									/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
									ETR_temp  = ETP[lay][i][j]; /* adding all of the potential to the actual evaporation */
									POND[i][j] -= ETP[lay][i][j]; /* evaporation is taken from the pond */
									ETP[lay][i][j] = 0.0;       /* there is nothing left for further evaporation */
								}

							}
							else
							{
								if (nfkratio < 0.6)
								ETR_temp = nfkratio/0.6 * ETP[lay][i][j];
								else /* nfkratio > 0.6, alles kann verdunsten, ETP==ETR */
								ETR_temp = ETP[lay][i][j];

								SB[i][j] -= (ETR_temp/NFK[i][j]);
								if (SB[i][j]<0.0)
								{
									/* von ETR_temp wird, wenn soviel gar nicht mehr im Boden war, der zuviel entnommene Anteil abgezogen,
							indem die nun negative Bodenfeuchte nach Umrechnen in mm zu ETR_temp addiert wird */
									ETR_temp += (SB[i][j]*NFK[i][j]);
									SB[i][j]=0; /* dann ist der Boden wirklich leer */
								}
							}

							ETR[lay][i][j]+=ETR_temp; /* zu der in ::snow berechneten ETR muss nun die echte Entnahme aus dem Boden addiert werden,
					um auf eine gesamte reale Verdunstung zu kommen */

							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETR[lay][i][j],1.0-icevalue);
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);

							if (!((lake_active==1)&&(LAKE[i][j]!=nodata)&&(POND[i][j]>1.0))) // NO Lake
							{
								EVAR_temp = EVAP[i][j];
								SB[i][j] -= (EVAR_temp/NFK[i][j]);
								if (SB[i][j]<0.0)
								{
									/* von ETR_temp wird, wenn soviel gar nicht mehr im Boden war, der zuviel entnommene Anteil abgezogen,
							indem die nun negative Bodenfeuchte nach Umrechnen in mm zu ETR_temp addiert wird */
									EVAR_temp += (SB[i][j]*NFK[i][j]);
									SB[i][j]=0; /* dann ist der Boden wirklich leer */
								}
								EVAR[i][j]+=EVAR_temp; /* zu der in ::snow berechneten ETR muss nun die echte Entnahme aus dem Boden addiert werden,
						um auf eine gesamte reale Verdunstung zu kommen */
							}
							else
							{
								EVAR[i][j]=0.0;
							}
							EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVAR[i][j],1.0-icevalue);
							EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
							ETR[0][i][j]+=ETR[lay][i][j];

						}
						else
						{
							ETR[lay][i][j]=0;
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0); // 0 is valid for the entore cell, so no break up for glarierized and unglazierized fraction of a cell is required here
						}
					}
					if (m_iInterceptionIsActive>0)
					{
						ETR[0][i][j]=ETR[0][i][j]+EI[i][j]+EVRS[i][j]+EVAR[i][j];
					}
					else
					{
						ETR[0][i][j]=ETR[0][i][j]+EVRS[i][j]+EVAR[i][j];
					}
					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETR[0][i][j],1.0-icevalue);
					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue); // add 0 for the unglazierîzed part of the actual cell in order to get statistics for the entire subbasin
				}
//				if (icevalue>0) // add 0 for glacierized fractions of a cell in order to get weighed statistics for the entire basin
//				{
//					for (lay=0;lay<=iNumVegetationLayers;lay++)
//					{
//						ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
//					}
//					EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
//				}
			}
			else
			{
				for (lay=0;lay<=iNumVegetationLayers;lay++)
				{
					ETR[lay][i][j] = nodata;
				}
			}
		}
	}
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETRstat[lay]->calc_and_done(jahr,monat,tag,stunde, minute);
		writegridtofile(etrindex[lay]);
	}
	EVARstat->calc_and_done(jahr,monat,tag,stunde, minute);
	writegridtofile(evarindex);
}



void verdunstungsmodellklasse::pre_real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep)
{
	int i,j;
	if (do_it == 0) return;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	
	teilmodellklasse::run(me, se, tg, mt, jr, (char *)_ETR_, NULL);

	// if this is the first time, the evaporation and transpiration components are summarized, look for interception evaporation
	if (m_iInterceptionIsActive==0)
	{
		char layername[2000];
		sprintf(layername,"%s%d",_Interzeptionsverdunstung_,0);
		int pos=ist_vorhanden(layername);
		if (pos>0)
		{
			EI=gridinfo[getindex(layername)].pointer;
			m_iInterceptionIsActive=1;
		}
	}

	int lay;
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETRstat[lay]->init();
	}
	EVARstat->init();

}

/* hier wird nur noch die bereits in unsatzonmodell berechnete reale Verdunstung in die
Statistikdatei geschrieben */
void verdunstungsmodellklasse::real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep)
{
	int i,j;
	if (do_it == 0) return;
	// accumulation is done for all substeps only once since the etp/etr/avar values are accumulated already and set to 0 in etpetr::potentiell only
}

void verdunstungsmodellklasse::post_real_uzr(long me, long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;

	int i,j;

	for (i = mpistartrow; i <= mpiendrow; i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			// no special treatment for HRU and non-HRU Master cells neccessary since ETR will be distributed to all cells already in Unsatzon and this routine is only for updating statistics registers
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//      if (zn[i][j]!=nodata)
			{
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

//				if (icevalue<1.0)
//				{
					for (int lay=1;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]>nodata)
						{
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETR[lay][i][j],1.0-icevalue);
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
						}
						else
						{
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,1.0-icevalue);
							ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
						}

					}
					// ETRstat[0] will contain the total evapotranspiration of interception, transpiration, soil evaporation and snow evaporation
					// but ETR[0] alrady contains EVAR (was added in unsatzonmodel)
					if (m_iInterceptionIsActive>0)
					{
						ETR[0][i][j] = EI[i][j]+EVRS[i][j]+EVAR[i][j];
						for (int lay=1;lay<=iNumVegetationLayers;lay++)
						{
							ETR[0][i][j]+=ETR[lay][i][j];
						}
						
					}
					else
					{
						ETR[0][i][j]=ETR[0][i][j]+EVRS[i][j];
					}
//					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETR[0][i][j],1.0-icevalue);
//					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVRS[i][j],icevalue); 
					double statval=(EVRSGL!=NULL)?EVRSGL[i][j]*icevalue:EVRS[i][j]*icevalue;
					statval+=ETR[0][i][j]*(1.0-icevalue);
					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],statval);

					// old:	ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],ETR[0][i][j]); // attention: EVAR was added in unsatzonmodel arleady
					EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVAR[i][j],1.0-icevalue);
					EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
//				}
//				if (icevalue>0) // add 0 for glacierized fractions of a cell in order to get weighed statistics for the entire basin
//				{
//					for (int lay=1;lay<=iNumVegetationLayers;lay++)
//					{
//						ETRstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
//					}
//					// add only EVRS for ETR since this entity is also taken from glacierized fractions
//					ETRstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EVRS[i][j],icevalue);
//					EVARstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
//				}
			}
		}
	}

	for (int lay=0;lay<=iNumVegetationLayers;lay++)
	{
		ETRstat[lay]->calc_and_done(jahr,monat,tag,stunde, minute);
		writegridtofile(etrindex[lay]);
	}
	EVARstat->calc_and_done(jahr,monat,tag,stunde, minute);
	writegridtofile(evarindex);
}

/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void verdunstungsmodellklasse::ClearMemory()
{
	if (do_it > 0)
	{
		//    if (ssdstat) free(ssdstat); ssdstat=NULL;
		if (cos_za) free(cos_za); cos_za=NULL;

		if (LUSE ) delete[] LUSE; LUSE=NULL;
		if (LAI  ) delete[] LAI;  LAI =NULL;
		if (WURZ ) delete[] WURZ; WURZ=NULL;
		if (VBG  ) delete[] VBG;  VBG =NULL;
		if (RSI  ) delete[] RSI;  RSI =NULL;
		if (RSC  ) delete[] RSC;  RSC =NULL;
		if (Z0   ) delete[] Z0  ; Z0  =NULL;
		if (ARCVEGFRACTIONS) delete[] ARCVEGFRACTIONS; ARCVEGFRACTIONS = NULL;

		RGEXstat->ClearMemory();
		if (RGEXstat) delete(RGEXstat); RGEXstat=NULL;
		if (zoneRGEX) free(zoneRGEX); zoneRGEX=NULL;
		if (RGEXstatfilename )  free(RGEXstatfilename);  RGEXstatfilename=NULL;

		CLstat->ClearMemory();
		if (CLstat) delete(CLstat); CLstat=NULL;
		if (zoneCL) free(zoneCL); zoneCL=NULL;
		if (CLstatfilename )  free(CLstatfilename);  CLstatfilename=NULL;

		EVAPstat->ClearMemory();
		if (EVAPstat) delete(EVAPstat); EVAPstat=NULL;
		if (zoneEVAP) free(zoneEVAP); zoneEVAP=NULL;
		if (EVAPstatfilename )  free(EVAPstatfilename);  EVAPstatfilename=NULL;

		EVARstat->ClearMemory();
		if (EVARstat) delete(EVARstat); EVARstat=NULL;
		if (zoneEVAR) free(zoneEVAR); zoneEVAR=NULL;
		if (EVARstatfilename )  free(EVARstatfilename);  EVARstatfilename=NULL;

		EVRSstat->ClearMemory();
		if (EVRSstat) delete(EVRSstat); EVRSstat=NULL;
		if (zoneEVRS) free(zoneEVRS); zoneEVRS=NULL;
		if (EVRSstatfilename )  free(EVRSstatfilename);  EVRSstatfilename=NULL;

		int i;
		for (i=0;i<=iNumVegetationLayers;i++)
		{
			ETPstat[i]->ClearMemory();
			EIPstat[i]->ClearMemory();
			ETRstat[i]->ClearMemory();

			if (ETPstat[i]) delete(ETPstat[i]); ETPstat[i]=NULL;
			if (EIPstat[i]) delete(EIPstat[i]); EIPstat[i]=NULL;
			if (ETRstat[i]) delete(ETRstat[i]); ETRstat[i]=NULL;
			
			if (zoneETP[i]) free(zoneETP[i]); zoneETP[i]=NULL;
			if (zoneEIP[i]) free(zoneEIP[i]); zoneEIP[i]=NULL;
			if (zoneETR[i]) free(zoneETR[i]); zoneETR[i]=NULL;

			if (ETPstatfilename[i]) free(ETPstatfilename[i]); ETPstatfilename[i]=NULL;
			if (EIPstatfilename[i]) free(EIPstatfilename[i]); EIPstatfilename[i]=NULL;
			if (ETRstatfilename[i]) free(ETRstatfilename[i]); ETRstatfilename[i]=NULL;
		}

		for (i=0;i<=iNumVegetationLayers;i++)
		{
			free_matrix(ETP[i], gridinfo[etpindex[i]].nrl, gridinfo[etpindex[i]].nrh, 1, ncols,nrows); ETP[i] = NULL;
			free_matrix(EIP[i], gridinfo[eipindex[i]].nrl, gridinfo[eipindex[i]].nrh, 1, ncols, nrows); EIP[i] = NULL;
			free_matrix(ETR[i], gridinfo[etrindex[i]].nrl, gridinfo[etrindex[i]].nrh, 1, ncols, nrows); ETR[i] = NULL;

			if (ETPgridname[i]) free(ETPgridname[i]); ETPgridname[i]=NULL;
			if (EIPgridname[i]) free(EIPgridname[i]); EIPgridname[i]=NULL;
			if (ETRgridname[i]) free(ETRgridname[i]); ETRgridname[i]=NULL;

			delete_grid(etpindex[i]);
			delete_grid(eipindex[i]);
			delete_grid(etrindex[i]);

		}
		if (RGEXstatfilename) free(RGEXstatfilename); RGEXstatfilename=NULL;
		if (CLstatfilename) free(CLstatfilename); CLstatfilename=NULL;

		free_matrix(EVAP, gridinfo[evapindex].nrl, gridinfo[evapindex].nrh, 1, ncols, nrows); EVAP = NULL;
		free_matrix(EVAR, gridinfo[evarindex].nrl, gridinfo[evarindex].nrh, 1, ncols, nrows); EVAR = NULL;
		free_matrix(EVRS, gridinfo[evrsindex].nrl, gridinfo[evrsindex].nrh, 1, ncols, nrows); EVRS = NULL;
		free_matrix(RGEX, gridinfo[rgexindex].nrl, gridinfo[rgexindex].nrh, 1, ncols, nrows); RGEX = NULL;
		free_matrix(CL, gridinfo[clindex].nrl, gridinfo[clindex].nrh, 1, ncols, nrows); RGEX = NULL;
		free_matrix(SSDTemp, gridinfo[ssdtempindex].nrl, gridinfo[ssdtempindex].nrh, 1, ncols, nrows); SSDTemp = NULL;

		if (EVAPgridname) free(EVAPgridname); EVAPgridname=NULL;
		if (EVARgridname) free(EVARgridname); EVARgridname=NULL;
		if (EVRSgridname) free(EVRSgridname); EVRSgridname=NULL;
		if (RGEXgridname) free(RGEXgridname); RGEXgridname=NULL;
		if (CLgridname) free(CLgridname); CLgridname=NULL;
		
		delete_grid(evapindex);
		delete_grid(evarindex);
		delete_grid(evrsindex);
		delete_grid(rgexindex);
		delete_grid(ssdtempindex);
		delete_grid(clindex);

		if (ETPstatfilename )  delete[] ETPstatfilename;  ETPstatfilename=NULL;
		if (EIPstatfilename )  delete[] EIPstatfilename;  EIPstatfilename=NULL;
		if (ETRstatfilename )  delete[] ETRstatfilename;  ETRstatfilename=NULL;

		if (ETPstat         ) delete[] ETPstat;     ETPstat=NULL;
		if (EIPstat         ) delete[] EIPstat;     EIPstat=NULL;
		if (ETRstat         ) delete[] ETRstat;     ETRstat=NULL;

		if (ETPgridname     ) delete[] ETPgridname;  ETPgridname=NULL;     
		if (EIPgridname     ) delete[] EIPgridname;  EIPgridname=NULL;     
		if (ETRgridname     ) delete[] ETRgridname;  ETRgridname=NULL;     

		if (etpindex        ) delete[] etpindex;    etpindex=NULL;
		if (eipindex        ) delete[] eipindex;    eipindex=NULL;
		if (etrindex        ) delete[] etrindex;    etrindex=NULL;

		if (zoneETP         ) delete[] zoneETP;     zoneETP=NULL;
		if (zoneEIP         ) delete[] zoneEIP;     zoneEIP=NULL;
		if (zoneETR         ) delete[] zoneETR;     zoneETR=NULL;

                if (EIP             ) delete[] EIP; EIP=NULL;
                if (ETP             ) delete[] ETP; ETP=NULL;
                if (ETR             ) delete[] ETR; ETR=NULL;

		if (PEC) free(PEC); PEC=NULL;
		if (FH) free(FH); FH=NULL;
	}
	modelobject::ClearMemory();

}

void verdunstungsmodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
	if (tracer_count<=0) return;
}

void verdunstungsmodellklasse::substance_transport_delete()
{
	teilmodellklasse::substance_transport_delete();
	if (tracer_count<=0) return;
}

void verdunstungsmodellklasse::flushbuffers()
{
	if (do_it==0) return;
	else
	{
		int lay;
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			if (ETPstat[lay])  ETPstat[lay]->flushbuffer();
			if (EIPstat[lay])  EIPstat[lay]->flushbuffer();
			if (ETRstat[lay])  ETRstat[lay]->flushbuffer();
		}
		if (RGEXstat) RGEXstat->flushbuffer();
		if (CLstat)   CLstat->flushbuffer();
		if (EVARstat) EVARstat->flushbuffer();
		if (EVAPstat) EVAPstat->flushbuffer();
		if (EVRSstat) EVRSstat->flushbuffer();
	}
}

/* oooooooooooooooooooo Ende der Verdunstungsmodell-routinen oooooooooooooooooo */

