/*   $Revision: 44 $
	$Date: 13.02.11 16:32 $  */

/*
(c) 1995-2013 by Joerg Schulla
version history: see source safe / SVN comments

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <sstream>
#include <stddef.h>
#include "stringp.h"   
#include "modconst.h"  
#include "modef.h"     
#include "teilmod.h"   
#include "wasim.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

// ++++++++++++++++ snow model methods +++++++++++++++++++++++++++++++++++ 

class tupel
{
public:
	int i;
	int j;
	double z;
	tupel() {i=0;j=0;z=0;};
	tupel(int I,int J,double Z) {i=I;j=J;z=Z;};
	bool operator<(const tupel &t2) const;
};

bool tupel::operator<(const tupel &t2) const
{
  if (this->z < t2.z) 
	  return true;
  else 
	  return false;
}


void schneemodellklasse::initialize_lakes(std::string instanzname)
{	
	if (do_it == 0) return;
	char *tstr;

	if (!ist_vorhanden(_LAKE_CODES_))
	{
		(*WSstream) << "INFO: no LAKE-grid was found in "<<this->instanzname<<"\n"<<
		"Lake modelling is switched off for "<<this->instanzname<<"."<<"[end-entry]";
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

		// suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse 
		findeabschnitt(tstr);

		fscanfln(steuerdatei,"%d",&lake_active,"flow control parameter for ",tstr);
		// switch off lake model when in Spinup-mode
		if (doSpinUp) {
			lake_active=0;
		} 

		if (lake_active != 0)
		{
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
					(*WSstream) << "FATAL: \n\n=======================================\n"<<
					"!!!! ERROR !!!! \n\n"<<
					"\nLAKE-grid was found, but no POND-Grid. For using lakes,"<<
					"\nboth grids as well as the MAXPOND grid must be present."<<
					"\n\n=======================================\n"<< "[end-entry]";
					throw(28);
				}
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

double **schneemodellklasse::GetSnowStorAdress()
{	if (do_it==0)
	{
		return NULL;
	}
	else
	{
		return SSNOW;
	}
}

void schneemodellklasse::SetExpoKorrAdress(double **adresse)
{	if (do_it == 0)
	{
		return;
	}
	else
	{
		EXKORR = adresse;
		return;
	}
}

void schneemodellklasse::SetQgesAdresse(double **p_Qges)
{	if (do_it==0)
	return;
	else
	{
		Qges=p_Qges;
		initialized=1;
	}
}

void schneemodellklasse::StorageFileError(const char *message,double c1,double c2)
{
	(*WSstream) << "FATAL: (snow model) error in storage file content reading "<<message<<":\n"
	"value in actual control is:   "<<c1<<"\n"
	"but value in storage file is: "<<c2<<"\n"
	"WaSiM throws..."<<"[end-entry]";
	throw(71);
}

void schneemodellklasse::ReadFloatFromStoragefile(double *FValue,const char *ErrorMessage)
{
	double st_float;
	fscanf(storagefile,"%lf",&st_float);
	if (st_float!=0)
	{
		if (fabs((st_float- *FValue)/st_float)>1e-5) StorageFileError(ErrorMessage,(double)(*FValue),(double)st_float);
	}
	else
	{
		if (*FValue!=0) StorageFileError(ErrorMessage,(double)(*FValue),(double)st_float);
	}
}


void schneemodellklasse::ReadSlrContents(char *_storagefilename,int flag)
{
	int i,j,z,k;
	int st_int,st_int2;
	char *satz;
	int syntaxVersion = 2; // default syntax expected is for version 2 of the storage file
	double origTimestep = InV;
	double dMultipleOfOrigTimestep = 1.0;

	if ((_storagefilename!=NULL)&&(flag!=2)&&(flag!=4))
	{
		storagefilename=(char*)calloc(_stringsize_,sizeof(char));
	    sprintf(this->storagefilename,"%s",_storagefilename);
        // if readgrid in Unsatzonmodel is 2 or 4, then only the groundwater table is read in, but no storage file.
        if (flag!=0) // read the storage file for unsatzonmodel.readgrid==1 or 3
        {
            if ((Qice_els!=NULL)&&(Qfirn_els!=NULL)&&(Qsnow_els!=NULL))
            {
                satz=(char*)calloc(_stringsize_,sizeof(char));

                char *completedname = new char[_stringsize_+1];
                sprintf(completedname,"%s",storagefilename);
                ExtendFileName(DefaultInputDirectory,completedname);
                if ((storagefile = fopen(completedname,"rt"))==NULL)
                {
					(*WSstream) << "WARNING: in snow/glacier model: error opening storage content file "<<completedname<<"[end-entry]";
					(*WSstream) << "WARNING: glacier model will be run with cold start initialization values"<<"[end-entry]";
					if (satz!=NULL) free(satz);
					delete[] completedname;
					return;
                }

//                do
//                {
//					satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
//					if (feof(storagefile))
//					{
//						(*WSstream) << "WARNING: in snow/glacier model: [snow_model] not found in storage content file "<<completedname<<"[end-entry]";
//						(*WSstream) << "WARNING: glacier model will be run with cold start initialization values"<<"[end-entry]";
//						delete[] completedname;
//						fclose(storagefile);
//						if (satz!=NULL) free(satz);
//						return;
//					}
//                }
//                while (!strstr(satz,_Schneemodell_));
//
//                delete[] completedname;

				char cmpString[100];
				char* tmp = satz;
				sprintf(cmpString, "[snow_model_v2]");
				do
				{
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					if (feof(storagefile))
					{
						satz = tmp; // because fgets sets satz to NULL when feof has been reached, we reset the pointer to it's original value
						(*WSstream) << "INFORMATION: couldn't find section [snow_model_v2], looking now for [snow_model] in " << completedname << "[end-entry]";
						syntaxVersion = 1;
						fseek(storagefile, 0, 0); // reset storagefile
						sprintf(cmpString, "[snow_model]");
						do
						{
							satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
							if (feof(storagefile))
							{
								(*WSstream) << "FATAL: error locking for the section [unsatzon_model_v2] or [unsatzon_model] in " << completedname << "[end-entry]";
								throw(71);
							}
						} while (!strstr(satz, cmpString) && !feof(storagefile));
					}
				} while (!strstr(satz, cmpString) && !feof(storagefile));
				delete[] completedname;


				/* TG_number */
				fscanf(storagefile, "%d", &st_int); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (st_int != zonmax_routing) StorageFileError("number of subbasins", (double)zonmax_routing, (double)st_int);

				if (syntaxVersion >= 2)
				{
					/* timestep used to create that storagefile */
					fscanf(storagefile, "%lf", &origTimestep); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					dMultipleOfOrigTimestep = origTimestep / InV;
				}

				/* reading subbasin codes */
                for (j=1;j<=zonmax_routing;j++)
                {
					fscanf(storagefile,"%d",&st_int);
					if (st_int!=j)
					{
							StorageFileError("subbasin code",(double)j,(double)st_int);
					}

					std::string errmsg("Qice_slr for subbasin ");
					char basincode[100];
					sprintf(basincode,"%d",st_int);
					errmsg.append(basincode);
					fscanf(storagefile,"%lf",&(Qice_els[j]));
					Qice_els[j] /= dMultipleOfOrigTimestep;
    //		ReadFloatFromStoragefile(&(Qice_els[j]),errmsg.c_str());
                    satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

                    fscanf(storagefile,"%d",&st_int);
                    if (st_int!=j)
                    {
                            StorageFileError("subbasin code",(double)j,(double)st_int);
                    }
                    errmsg = "Qfirn_slr for subbasin ";
                    errmsg.append(basincode);
                    fscanf(storagefile,"%lf",&(Qfirn_els[j]));
					Qfirn_els[j] /= dMultipleOfOrigTimestep;
    //		ReadFloatFromStoragefile(&(Qfirn_els[j]),errmsg.c_str());
                    satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

                    fscanf(storagefile,"%d",&st_int);
                    if (st_int!=j)
                    {
						StorageFileError("subbasin code",(double)j,(double)st_int);
                    }
                    errmsg = "Qsnow_slr for subbasin ";
                    errmsg.append(basincode);
                    fscanf(storagefile,"%lf",&(Qsnow_els[j]));
					Qsnow_els[j] /= dMultipleOfOrigTimestep;
            //		ReadFloatFromStoragefile(&(Qsnow_els[j]),errmsg.c_str());

					satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
                }

                fclose(storagefile);
				if (satz!=NULL) free(satz);
	        }
		}
    }
	/* -----------ENDE STORAGEFILE -------------- */
}

void schneemodellklasse::WriteSlrContents(bool appendmode)
{  /* wrinting content of the three single linear reservoirs for the glacier runoff components to the storage file (either Richards version or the topmodel version) */
	int i,j,z;
//	if (mpi_rank==0)
//	{
		if (storagefilename!=NULL) // so it was initialized somewhen during start of the model, otherwise default is NULL (initialized in constructor)
		{
			char *completedname = new char[_stringsize_+1];
			sprintf(completedname,"%s",storagefilename);
			ExtendFileName(DefaultOutputDirectory,completedname);

			if (appendmode==true) // if routing model or unsatzon was active
				storagefile = fopen(completedname,"a");
			else
				storagefile = fopen(completedname,"w");

			delete[] completedname;

			fprintf(storagefile,"[" _Schneemodell_ "_v2]\n");
			fprintf(storagefile,"%d (number of subasins)\n",zonmax_routing);
			fprintf(storagefile, "%.3lf (timestep in minutes)\n", InV);

			for (j=1;j<=zonmax_routing;j++)
			{ 
				fprintf(storagefile,"%d\t%lf",j,Qice_els!=NULL?Qice_els[j]:0); fprintf(storagefile,"(TG code and Qice_slr value)\n");		
				fprintf(storagefile,"%d\t%lf",j,Qfirn_els!=NULL?Qfirn_els[j]:0); fprintf(storagefile,"(TG code and Qfirn_slr value)\n");		
				fprintf(storagefile,"%d\t%lf",j,Qsnow_els!=NULL?Qsnow_els[j]:0); fprintf(storagefile,"(TG code and Qsnow_slr value)\n");		
			} 

			fprintf(storagefile,"\n");
			fclose(storagefile);
		}
		/* -----------ENDE STORAGEFILE -------------- */
//	}
}

void schneemodellklasse::readSubBasinCodes(int *param)
{
//todo: extract codes from control file, fill missing ranges
	int myindex = 0;
	// int kzone = (int)(zone[1].zonenwert); 
	int kzone = (int)(zone_routing[1].zonenwert); /* first code is minmal for the lowest basin number */
	do /* reading ranges for subbasin codes and parameters */
	{
		if (fscanf(steuerdatei, "%d", &myindex) != 1)
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error reading code[" << (int)myindex << "] for " << instanzname << "\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n5.0\t6.0\n4.0\t4.5\n3.0\t3.4\n1.4\t1.45\n-0.0001\t-0.0000970\n0.0003\t0.00032\n0.00006\t0.00007\n0.00025\t0.00025\n12.0\t13.0\n240.0\t220.0\n24.0\t23.0\n0.01\t0.012\n0.01\t0.012\n0.01\t0.012\n3200\t2800\n1.80\t1.9\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!" << "[end-entry]";
			throw(10);
		}

		if ((myindex > zonmaxwert_routing) || (myindex < 1) || (zarray_routing[myindex] == 0))
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error! there is no code " << myindex << " (from " << steuerdateiname
				<< ") for a subbasin in " << instanzname << "\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n5.0\t6.0\n4.0\t4.5\n3.0\t3.4\n1.4\t1.45\n-0.0001\t-0.0000970\n0.0003\t0.00032\n0.00006\t0.00007\n0.00025\t0.00025\n12.0\t13.0\n240.0\t220.0\n24.0\t23.0\n0.01\t0.012\n0.01\t0.012\n0.01\t0.012\n3200\t2800\n1.80\t1.9\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!" << "[end-entry]";
			throw(10);
		}
		else
			param[zarray_routing[myindex]] = myindex;

		if (zarray_routing[myindex] > zarray_routing[kzone] + 1) /* filling codes inbetween with -1 */
			if ((kzone > zone_routing[1].zonenwert) || ((kzone == zone_routing[1].zonenwert) && (param[1] > 0)))
				for (int j = zarray_routing[kzone] + 1; j < zarray_routing[myindex]; j++) param[j] = -1;
			else
			{
				param[zarray_routing[kzone]] = zone_routing[1].zonenwert;
				for (int j = zarray_routing[kzone] + 1; j <= zarray_routing[myindex]; j++) param[j] = -1; /* moving the first valid subbasin code to the first subbasin, if first code was higer */
			}
		kzone = myindex;
	} while (myindex < zonmaxwert_routing);

	tstr = fgets(tstr, _stringsize_, steuerdatei);



}

void schneemodellklasse::readGlacierParam(double *param)
{
//todo: read parameters per subbasin from control file
	for (int i = 1; i <= zonmax_routing; i++)
	{
		if (Codes[i] > 0)
		{
			if (fscanf(steuerdatei, "%lf", &(param[i])) != 1)
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading kelsqd[" << (int)Codes[i] << "] for " << instanzname << "\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n5.0\t6.0\n4.0\t4.5\n3.0\t3.4\n1.4\t1.45\n-0.0001\t-0.0000970\n0.0003\t0.00032\n0.00006\t0.00007\n0.00025\t0.00025\n12.0\t13.0\n240.0\t220.0\n24.0\t23.0\n0.01\t0.012\n0.01\t0.012\n0.01\t0.012\n3200\t2800\n1.80\t1.9\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!" << "[end-entry]";
				throw(10);
			}
		}
		else param[i] = param[i - 1];
	}  
	tstr = fgets(tstr, _stringsize_, steuerdatei);


}


schneemodellklasse::schneemodellklasse(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,
	int maxLastResults,double breitearg,double laengearg,double bezugsmeridianarg,
	char *bezeichnung,char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, 
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll,double p_yll,double p_cellsize,double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols,p_nrows,p_mpistartrow,p_mpiendrow,p_xll,p_yll,p_cellsize,p_nodata,mine,ste,tge,mte,jre,p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode, p_vHRU),
qschnee(0),snowstor(0),snowstor2(0),zoneQSNOW(0),zoneSSNOW(0),zoneSSNOW2(0),qschneestatcode(0),snowstorstatcode(0),
qschneegridname(0),qschneestatname(0),snowstorstatname(0),snowstorstatname2(0),SSNOWgridname(0),SLIQgridname(0),SSTOgridname(0),
qsnow_out(0),qsnow_outgridname(0),qsnow_outstatname(0),qsnow_outstatcode(0),qsnow_outwritegrid(0),
qsnow_outindex(0),QSNOW_OUT(0),zoneQSNOW_OUT(0),
SNOW_AGE(0),SNOW_AGEgridname(0),SNOW_AGEwritegrid(0),snow_ageindex(0),
snow_age(0),snow_agestatname(0),snow_agestatcode(0),zoneSNOW_AGE(0),
SNOWDAYS(0),SNOWDAYSgridname(0),SNOWDAYSwritegrid(0),snowdaysindex(0),
snowdays(0),snowdaysstatname(0),snowdaysstatcode(0),zoneSNOWDAYS(0),
SNOWFREEDAYS(0),snowfreedaysindex(0),ALPHA(0),alphaindex(0),THAWDEPTH(0),thawdepthindex(0),ROOTDEPTH(0),
SNOWCOVERDAYS(0),snowcoverdaysindex(0),
P_SNOW(0),P_SNOWgridname(0),P_SNOWwritegrid(0),p_snowindex(0),
p_snow(0),p_snowstatname(0),p_snowstatcode(0),zoneP_SNOW(0),
P_RAIN(0),P_RAINgridname(0),P_RAINwritegrid(0),p_rainindex(0),
p_rain(0),p_rainstatname(0),p_rainstatcode(0),zoneP_RAIN(0),
Albedo(0),albedostatname(0),albedostatcode(0),zoneALBEDO(0), Albedowritegrid(0),albedoindex(0),
T0R(0),T0(0),CWH(0),CRFR(0),actC0(0),actC1(0),actC2(0),C0(0),C1(0),C2(0),
z0(0),RMFMIN(0),RMFMAX(0),SCA(0),
ALBsnow_min(0),ALBsnow_max(0),
P(0),T(0),Tnacht(0),U(0),E(0),QSNOW(0),RMFakt(0),
SSNOW(0),SLIQ(0),SSTO(0),LUSE(0),EXKORR(0),Qges(0),r_a(0),r_b(0),r_c(0),r_d(0),
Asnow_a(0),Asnow_b(0),Asnow_min(0),Asnow_max(0), Asnow_akt(0),
Aice_a(0),Aice_b(0),Aice_akt(0),Aice_min(0),Aice_max(0),
transzone(0),ssnowindex(0),sliqindex(0),sstoindex(0),weindex(0),albindex(0),qsnowindex(0),
SSNOWwritegrid(0),Qwritegrid(0),SLIQwritegrid(0),SSTOwritegrid(0),glacmethode(0),
qfirngridname(0),qicegridname(0),
qfirnstatname(0),qicestatname(0),qglacstatname(0),qsnowonglacstatname(0),
firnstackname(0),firnstatname(0),gmbgridname(0),gmbogridname(0),
gmbwritegrid(0),gmbowritegrid(0),gmbindex(0),gmboindex(0),
gmbstatcode(0),gmbstatname(0),
gmb2statcode(0),gmb2statname(0),
qfirnwritegrid(0),qicewritegrid(0),
qfirnindex(0),qiceindex(0),
ICE(0),QFIRN(0),QICE(0),RAD(0),GLACID(0),GLACDEB(0),
qfirnstatcode(0),qicestatcode(0),qglacstatcode(0),qsnowonglacstatcode(0),
qfirn(0),qice(0),qglac(0),qsnowonglac(0),
zoneQICE(0),zoneQFIRN(0),zoneQGLAC(0),zoneQSnowOnGlac(0),
zoneFIRN(0),zoneGMB(0),zoneGMB2(0),firnstat(0),gmbstat(0),gmb2stat(0),
MF(0),DDFsnow(0),DDFfirn(0),DDFice(0),
Kice(0),Kfirn(0),Ksnow(0),kfirn_els(0), ksnow_els(0), kice_els(0),
Qfirn(0), Qsnow(0), Qice(0),Qfirn_els(0), Qsnow_els(0), Qice_els(0),glacierlist(0),
rho_ice(0.85),
LAKE(0),POND(0),lake_active(0),
thawdepth(0),zoneTHAWDEPTH(0),thawdepthstatcode(0),thawdepthstatname(0),
storagefilename(0),storagefile(0),
zn_routing(0),zone_routing(0),zarray_routing(0),zonmax_routing(0),zonmaxwert_routing(0),dGloabDebrisScaling(1.0),
//ETRS(0),etrsindex(0),ETRSchecked(false),
RLF(0), RGEX(0),SL(0),CL(0), 
SNOWTEMP(0),SNOWTEMPgridname(0),SNOWTEMPwritegrid(0),snowtempindex(0), 
snowtemp(0),snowtempstatname(0),snowtempstatcode(0),zoneSNOWTEMP(0),
SURFTEMP(0),SURFTEMPgridname(0),SURFTEMPwritegrid(0),surftempindex(0), 
surftemp(0),surftempstatname(0),surftempstatcode(0),zoneSURFTEMP(0),
SLIDEFRACTION1(0),SLIDEFRACTION2(0),SLIDEFRACTION3(0),SLIDEFRACTION4(0),
ELEV_ORDER_ROW(0),ELEV_ORDER_COL(0),
INPUT_MASS(0),INPUT_MASSgridname(0),INPUT_MASSwritegrid(0),input_massindex(0),
input_mass(0),input_massstatname(0),input_massstatcode(0),zoneINPUT_MASS(0),
MOBILE_MASS(0),MOBILE_MASSgridname(0),MOBILE_MASSwritegrid(0),mobile_massindex(0),
mobile_mass(0),mobile_massstatname(0),mobile_massstatcode(0),zoneMOBILE_MASS(0),
DEPOSITION(0),DEPOSITIONgridname(0),DEPOSITIONwritegrid(0),depositionindex(0),
DEPOSITION_MAX(0),depositionmaxindex(0),maxDepositionSlope(55),maxDepositionFactor(1),
deposition(0),depositionstatname(0),depositionstatcode(0),zoneDEPOSITION(0),
minSlideSlope(30),slideMassFrac(0.01),LWINcorr(1.0),LWOUTcorr(1.0),EVRS(NULL),
ICEOLD(0),iceoldindex(0),GLACIDOLD(0),glacidoldindex(0),
SSNOWGL(0),SLIQGL(0),SSTOGL(0),QSNOWGL(0),SNOWTEMPGL(0),SURFTEMPGL(0),EVRSGL(0),
ssnowglindex(0),sliqglindex(0),sstoglindex(0),qsnowglindex(0),snowtempglindex(0),surftempglindex(0),evrsglindex(0),
SSNOWGLwritegrid(0),SLIQGLwritegrid(0),SSTOGLwritegrid(0),QSNOWGLwritegrid(0),SNOWTEMPGLwritegrid(0),SURFTEMPGLwritegrid(0),EVRSGLwritegrid(0),
SSNOWGLgridname(0),SLIQGLgridname(0),SSTOGLgridname(0),QSNOWGLgridname(0),SNOWTEMPGLgridname(0),SURFTEMPGLgridname(0),EVRSGLgridname(0),
SURFTEMPINTERPOL(0), bFullGridRequired(false), sr(mpistartrow), er(mpiendrow), interceptionBeforeSnow(false), 
maxSnowPackToCreateGlacier(1000), relSizeNewGlacier(0.25), EnBalScaling(1.0), bHeatTransferActive(false), bLayeredSnowActive(false),
EBAL(NULL), ebalindex(0), EBALwritegrid(0), EBALgridname(NULL),
EBALGL(NULL), ebalglindex(0), EBALGLwritegrid(0), EBALGLgridname(NULL), 
refSWES(NULL), refSWEL(NULL), LayeredSnowProperties(NULL),
refUpperIndexSnow(NULL), refThickSnow(NULL), refDensSnow(NULL), refTemperSnow(NULL),
SCN_ir_radiation(NULL), SCN_snow_throuhfall(NULL), DetailedGlacierParams(0)
{																	   
	try
	{
		retvalue=0;
		int i,j,k,sti;


		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;

		// belegen von Speicher 

#ifdef __TIME__
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
#endif

		DHM = gridinfo[2].pointer; // Hoehenmodell von jetzt an unter DHM ansprechbar 

		zonenidentification(&zn_routing,&zone_routing,&zarray_routing,zonmax_routing,zonmaxwert_routing);

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

			if (ist_vorhanden(_MPIHelper_))
			{
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			}
			else
			{
				MPIHLP = NULL;
			}
		}


		// old values (for compatibility testing)
//		zn_routing=zn;
//		zone_routing=zone;
//		zarray_routing=zarray;
//		zonmax_routing=zonmax;
//		zonmaxwert_routing=zonmaxwert;

		if (MultiLayerLUTable!=NULL)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}
		LUSE   = new double**[iNumVegetationLayers+1];
		
		qschneegridname     = newstring(_stringsize_+1,1);
		qschneestatname     = newstring(_stringsize_+1,1);
		snowstorstatname    = newstring(_stringsize_+1,1);
		SSNOWgridname       = newstring(_stringsize_+1,1);
		SNOW_AGEgridname    = newstring(_stringsize_+1,1);
		snow_agestatname    = newstring(_stringsize_+1,1);
		SNOWDAYSgridname    = newstring(_stringsize_+1,1);
		snowdaysstatname    = newstring(_stringsize_+1,1);
		P_SNOWgridname      = newstring(_stringsize_+1,1);
		p_snowstatname      = newstring(_stringsize_+1,1);
		P_RAINgridname      = newstring(_stringsize_+1,1);
		p_rainstatname      = newstring(_stringsize_+1,1);
		albedostatname      = newstring(_stringsize_+1,1);
		SLIQgridname        = newstring(_stringsize_+1,1);
		SSTOgridname        = newstring(_stringsize_+1,1);
		qsnow_outgridname   = newstring(_stringsize_+1,1);
		qsnow_outstatname   = newstring(_stringsize_+1,1);

		// lesen der Werte aus der Steuerdatei 

		// ACHTUNG !!!!!: die Steuerdatei wird im Constructor der Mutterklasse (teilmodell)
		// geoeffnet und bis zum Erreichen des
		// ignoranzparameters
		// +  Intervallaenge
		// gelesen.
		// ist der Ignoranzparameter = 0, dann wird die Datei nach diesem Parameter (InV) nicht
		// weitergelesen und gleich geschlossen.
		// So kann der gesamte Ablauf der ini-Prozedur abgebrochen werden, da der
		// Ignoranzparameter schon bekannt ist!
		

		fscanfln(steuerdatei,"%d",&methode   ,"method for ",tstr);
		fscanfln(steuerdatei,"%lf",&transzone ,"transzone for ",tstr);
		fscanfln(steuerdatei,"%lf",&T0R   ,"T0R for "     ,tstr);
		fscanfln(steuerdatei,"%lf",&T0    ,"T0 for "     ,tstr);
		fscanfln(steuerdatei,"%lf",&CWH   ,"CWH for "    ,tstr);
		fscanfln(steuerdatei,"%lf",&CRFR  ,"CRFR for "   ,tstr);

		int cnt = 0;
		C0 = new double[14];
		C1 = new double[14];
		C2 = new double[14];

		satz = fgets(satz, _stringsize_, steuerdatei);
		cnt = sscanf(satz, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &C0[1], &C0[2], &C0[3], &C0[4], &C0[5], &C0[6], &C0[7], &C0[8], &C0[9], &C0[10], &C0[11], &C0[12]);
		if (cnt==12)
		{
			// all 12 values where read
			C0[0]  = C0[12];
			C0[13] = C0[1];
		}
		else if (cnt>=1) {
			// at least one value was read; all values get the value from the first element
			C0[0] = C0[2] = C0[3] = C0[4] = C0[5] = C0[6] = C0[7] = C0[8] = C0[9] = C0[10] = C0[11] = C0[12] = C0[13] = C0[1];
		}

		satz = fgets(satz, _stringsize_, steuerdatei);
		cnt = sscanf(satz, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &C1[1], &C1[2], &C1[3], &C1[4], &C1[5], &C1[6], &C1[7], &C1[8], &C1[9], &C1[10], &C1[11], &C1[12]);
		if (cnt == 12)
		{
			// all 12 values where read
			C1[0]  = C1[12];
			C1[13] = C1[1];
		}
		else if (cnt >= 1) {
			// at least one value was read; all values get the value from the first element
			C1[0] = C1[2] = C1[3] = C1[4] = C1[5] = C1[6] = C1[7] = C1[8] = C1[9] = C1[10] = C1[11] = C1[12] = C1[13] = C1[1];
		}

		satz = fgets(satz, _stringsize_, steuerdatei);
		cnt = sscanf(satz, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &C2[1], &C2[2], &C2[3], &C2[4], &C2[5], &C2[6], &C2[7], &C2[8], &C2[9], &C2[10], &C2[11], &C2[12]);
		if (cnt == 12)
		{
			// all 12 values where read
			C2[0] = C2[12];
			C2[13] = C2[1];
		}
		else if (cnt >= 1) {
			// at least one value was read; all values get the value from the first element
			C2[0] = C2[2] = C2[3] = C2[4] = C2[5] = C2[6] = C2[7] = C2[8] = C2[9] = C2[10] = C2[11] = C2[12] = C2[13] = C2[1];
		}


		// old
		//fscanfln(steuerdatei, "%lf", &C0, "C0 for ", tstr);
		//fscanfln(steuerdatei, "%lf", &C1, "C1 for ", tstr);
		//fscanfln(steuerdatei, "%lf", &C2, "C2 for ", tstr);
		
		fscanfln(steuerdatei,"%lf",&z0    ,"z0 for "     ,tstr);
		fscanfln(steuerdatei,"%lf",&RMFMIN,"RMFMIN for " ,tstr);
		fscanfln(steuerdatei,"%lf",&RMFMAX,"RMFMAX for " ,tstr);
		fscanfln(steuerdatei,"%lf",&ALBsnow_min,"ALBsnow_min for ",tstr);
		fscanfln(steuerdatei,"%lf",&ALBsnow_max,"ALBsnow_max for ",tstr);

		fscanfln(steuerdatei,"%s",P_RAINgridname     ,"Name of the rain-grid for "  ,tstr);
		fscanfln(steuerdatei,"%d",&P_RAINwritegrid   ,"Writegrid of rain-grid for " ,tstr);
		fscanfln(steuerdatei,"%s%d",p_rainstatname,&p_rainstatcode ,"Name of the rain statistics for ",tstr);
		fscanfln(steuerdatei,"%s",P_SNOWgridname     ,"Name of the snow grid for "      ,tstr);
		fscanfln(steuerdatei,"%d",&P_SNOWwritegrid   ,"Writegrid of the snow grid for " ,tstr);
		fscanfln(steuerdatei,"%s%d",p_snowstatname,&p_snowstatcode ,"Name of the snow statistics for ",tstr);
		fscanfln(steuerdatei,"%s",SNOWDAYSgridname     ,"Name of the snowdays grid for "      ,tstr);
		fscanfln(steuerdatei,"%d",&SNOWDAYSwritegrid   ,"Writegrid of the snowdays grid for " ,tstr);
		fscanfln(steuerdatei,"%s%d",snowdaysstatname,&snowdaysstatcode ,"Name of the snowdays statistics for " ,tstr);
		fscanfln(steuerdatei,"%s",SNOW_AGEgridname     ,"Name of the snow age grid for "      ,tstr);
		fscanfln(steuerdatei,"%d",&SNOW_AGEwritegrid   ,"Writegrid of the snow age grid for " ,tstr);
		fscanfln(steuerdatei,"%s%d",snow_agestatname,&snow_agestatcode ,"Name of the snow age statistics for " ,tstr);
		fscanfln(steuerdatei,"%s%d",albedostatname,&albedostatcode ,"name of the albedo statistics for ",tstr);
		fscanfln(steuerdatei,"%s",qschneegridname      ,"Name of the snow outflow grid for "            ,tstr);
		fscanfln(steuerdatei,"%d",&Qwritegrid          ,"Writegrid of the snow outflow grid for "       ,tstr);
		fscanfln(steuerdatei,"%s%d",qschneestatname,&qschneestatcode ,"Name of the snow outflow statistics for ",tstr);
		fscanfln(steuerdatei,"%s",qsnow_outgridname    ,"Name of the snow melt only outflow grid for "  ,tstr);
		fscanfln(steuerdatei,"%d",&qsnow_outwritegrid  ,"Writegrid of the snow melt only outflow grid for "     ,tstr);
		fscanfln(steuerdatei,"%s%d",qsnow_outstatname,&qsnow_outstatcode ,"Name of the snow melt only outflow statistics for ",tstr);
		fscanfln(steuerdatei,"%s",SSNOWgridname        ,"Name of the snow storage grids for "           ,tstr);
		fscanfln(steuerdatei,"%d",&SSNOWwritegrid      ,"Writegrid of the snow storage grids for "      ,tstr);
		fscanfln(steuerdatei,"%s",SLIQgridname         ,"Name of the liquid storage grids for "         ,tstr);
		fscanfln(steuerdatei,"%d",&SLIQwritegrid       ,"Writegrid of the liquid storage grids for "    ,tstr);
		fscanfln(steuerdatei,"%s%d",snowstorstatname,&snowstorstatcode  ,"Name of the total snow storage statistics for " ,tstr);
		fscanfln(steuerdatei,"%s",SSTOgridname         ,"name of the total snow storage grid for "      ,tstr);
		fscanfln(steuerdatei,"%d",&SSTOwritegrid       ,"Writegrid of the total snow storage grid for " ,tstr);
		fscanfln(steuerdatei,"%d",&gridslesen   ,"readgrid forr storage grids for ",tstr);

		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14))
		{
			// attention: sr and er and bFullGridRequired may be set also depending on glacmethod later
			if (mpi_rank == 0)
			{
				bFullGridRequired = true; // in mpi_rank 0, the full grid is required for a series of entities in order to do the slide model
				sr = 1;
				er = nrows;
			}
		}

		// new snow model extensions (Warscher 2013)
		if (methode>4) // all parameters should be there, but the grids will only be created or looked for or read in if it is required by hte methode's value 
		{
			SNOWTEMPgridname          = newstring(_stringsize_+1,1);
			snowtempstatname          = newstring(_stringsize_+1,1);
			SURFTEMPgridname          = newstring(_stringsize_+1,1);
			surftempstatname          = newstring(_stringsize_+1,1);
			
			EBALgridname			  = newstring(_stringsize_ + 1, 1);
			
			INPUT_MASSgridname        = newstring(_stringsize_+1,1);
			input_massstatname        = newstring(_stringsize_+1,1);
			MOBILE_MASSgridname       = newstring(_stringsize_+1,1);
			mobile_massstatname       = newstring(_stringsize_+1,1);
			DEPOSITIONgridname        = newstring(_stringsize_+1,1);
			depositionstatname        = newstring(_stringsize_+1,1);

			fgets(satz,_stringsize_,steuerdatei);
			if (strlen(satz)<=2)
			{
				(*WSstream) << "FATAL: cannot read additional parameters for snow model methods 5 to 11 from "<<steuerdateiname<<"[end-entry]";
				throw(404);
			}
			sscanf(satz,"%s",SNOWTEMPgridname);

			fscanfln(steuerdatei,"%d",&SNOWTEMPwritegrid   ,"Writegrid of the snow pack temperature grid for "    ,tstr);
			fscanfln(steuerdatei,"%s%d",snowtempstatname,&snowtempstatcode ,"Name of the snow pack temperature statistics for ",tstr);
			fscanfln(steuerdatei,"%s",SURFTEMPgridname     ,"Name of the snow surface temperature grid for "      ,tstr);
			fscanfln(steuerdatei,"%d",&SURFTEMPwritegrid   ,"Writegrid of the snow surface temperature grid for " ,tstr);
			fscanfln(steuerdatei,"%s%d",surftempstatname,&surftempstatcode ,"Name of the snow surface temperature statistics for ",tstr);
			fscanfln(steuerdatei,"%s",INPUT_MASSgridname     ,"Name of the input mass grid for "        ,tstr);
			fscanfln(steuerdatei,"%d",&INPUT_MASSwritegrid   ,"Writegrid of the input mass grid for "   ,tstr);
			fscanfln(steuerdatei,"%s%d",input_massstatname,&input_massstatcode ,"Name of the input mass statistics for ",tstr);
			fscanfln(steuerdatei,"%s",MOBILE_MASSgridname     ,"Name of the mobile mass grid for "      ,tstr);
			fscanfln(steuerdatei,"%d",&MOBILE_MASSwritegrid   ,"Writegrid of the mobile mass grid for " ,tstr);
			fscanfln(steuerdatei,"%s%d",mobile_massstatname,&mobile_massstatcode ,"Name of the mobile mass statistics for ",tstr);
			fscanfln(steuerdatei,"%s",DEPOSITIONgridname     ,"Name of the deposition grid for "        ,tstr);
			fscanfln(steuerdatei,"%d",&DEPOSITIONwritegrid   ,"Writegrid of the deposition grid for "   ,tstr);
			fscanfln(steuerdatei,"%s%d",depositionstatname,&depositionstatcode ,"Name of the deposition statistics for ",tstr);
			fscanfln(steuerdatei,"%lf",&maxDepositionSlope,"maximum slope for deposition for ",tstr);
			if ((maxDepositionSlope<0)||(maxDepositionSlope>90))
			{
				(*WSstream) << "FATAL: maxDepositionSlope not within parameter bounds of [0...90] Here: "<<maxDepositionSlope<<"[end-entry]";
				throw(405);
			}
			fscanfln(steuerdatei,"%lf",&maxDepositionFactor,"maximum deposition factor for ",tstr);
			if (maxDepositionFactor<=0)
			{
				(*WSstream) << "FATAL: maxDepositionFactor not within parameter bounds of (0...inf). Here: "<<maxDepositionFactor<<"[end-entry]";
				throw(406);
			}
			fscanfln(steuerdatei,"%lf",&minSlideSlope,"minimum slope for slides for ",tstr);
			if ((minSlideSlope<=0)||(minSlideSlope>90))
			{
				(*WSstream) << "FATAL: minSlideSlope not within parameter bounds of [0...90]. Here: "<<minSlideSlope<<"[end-entry]";
				throw(407);
			}
			fscanfln(steuerdatei,"%lf",&slideMassFrac,"slide mass fraction for ",tstr);
			if ((slideMassFrac<=0)||(slideMassFrac>1.0))
			{
				(*WSstream) << "FATAL: slideMassFrac not within parameter bounds of (0...1.0]. Here: "<<slideMassFrac<<"[end-entry]";
				throw(408);
			}
			fscanfln(steuerdatei,"%lf",&LWINcorr,"incoming long wave radiation correction factor for ",tstr);
			if ((LWINcorr<=0)||(LWINcorr>3.0))
			{
				(*WSstream) << "FATAL: incoming long wave radiation correction factor LWINcorr not within parameter bounds of (0...3.0]. Here: "<<LWINcorr<<"[end-entry]";
				throw(408);
			}
			fscanfln(steuerdatei,"%lf",&LWOUTcorr,"outgoing long wave radiation correction factor for ",tstr);
			if ((LWOUTcorr <= 0) || (LWOUTcorr>3.0))
			{
				(*WSstream) << "FATAL: outgoing long wave radiation correction factor LWOUTcorr not within parameter bounds of (0...3.0]. Here: "<<LWOUTcorr<<"[end-entry]";
				throw(408);
			}
			if (fscanf(steuerdatei, "%lf", &EnBalScaling) != 1)
			{
				(*WSstream) << "WARNING: \n------------- warning: -------------\n"
					<< "Could not read the optional parameter EnBalScaling (double),\n"
					<< "for snow melt methods with energy balance that will be used for scaling the energy balance (i.e. affect melt rates)\n"
					<< "Setting default values of 1.0" << "[end-entry]";
			}
			tstr = fgets(tstr, _stringsize_, steuerdatei);

		}

		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: cannot close "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		zoneALBEDO = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneALBEDO) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneALBEDO[i].n           = zoneALBEDO[i].sx        = 0;
			zoneALBEDO[i].zonenwert   = zone[i].zonenwert;
			zoneALBEDO[i].rel         = zone[i].rel;
			zoneALBEDO[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneQSNOW = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQSNOW) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneQSNOW[i].n           = zoneQSNOW[i].sx        = 0;
			zoneQSNOW[i].zonenwert   = zone[i].zonenwert;
			zoneQSNOW[i].rel         = zone[i].rel;
			zoneQSNOW[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneQSNOW_OUT = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQSNOW_OUT) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneQSNOW_OUT[i].n           = zoneQSNOW_OUT[i].sx        = 0;
			zoneQSNOW_OUT[i].zonenwert   = zone[i].zonenwert;
			zoneQSNOW_OUT[i].rel         = zone[i].rel;
			zoneQSNOW_OUT[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneSSNOW = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSSNOW) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneSSNOW[i].n           = zoneSSNOW[i].sx        = 0;
			zoneSSNOW[i].zonenwert   = zone[i].zonenwert;
			zoneSSNOW[i].rel         = zone[i].rel;
			zoneSSNOW[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneSNOWDAYS = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSNOWDAYS) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneSNOWDAYS[i].n           = zoneSNOWDAYS[i].sx        = 0;
			zoneSNOWDAYS[i].zonenwert   = zone[i].zonenwert;
			zoneSNOWDAYS[i].rel         = zone[i].rel;
			zoneSNOWDAYS[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneSNOW_AGE = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSNOW_AGE) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneSNOW_AGE[i].n           = zoneSNOW_AGE[i].sx        = 0;
			zoneSNOW_AGE[i].zonenwert   = zone[i].zonenwert;
			zoneSNOW_AGE[i].rel         = zone[i].rel;
			zoneSNOW_AGE[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneP_SNOW = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneP_SNOW) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneP_SNOW[i].n           = zoneP_SNOW[i].sx        = 0;
			zoneP_SNOW[i].zonenwert   = zone[i].zonenwert;
			zoneP_SNOW[i].rel         = zone[i].rel;
			zoneP_SNOW[i].RelCellSum  = zone[i].RelCellSum;
		}

		zoneP_RAIN = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneP_RAIN) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{
			zoneP_RAIN[i].n           = zoneP_RAIN[i].sx        = 0;
			zoneP_RAIN[i].zonenwert   = zone[i].zonenwert;
			zoneP_RAIN[i].rel         = zone[i].rel;
			zoneP_RAIN[i].RelCellSum  = zone[i].RelCellSum;
		}

		// new snow model extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			zoneSNOWTEMP = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			if (!zoneSNOWTEMP) speichererror(1);
			for (i = 0; i <= zonmax; i++)
			{
				zoneSNOWTEMP[i].n = zoneSNOWTEMP[i].sx = 0;
				zoneSNOWTEMP[i].zonenwert = zone[i].zonenwert;
				zoneSNOWTEMP[i].rel = zone[i].rel;
				zoneSNOWTEMP[i].RelCellSum = zone[i].RelCellSum;
			}

			zoneSURFTEMP = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			if (!zoneSURFTEMP) speichererror(1);
			for (i = 0; i <= zonmax; i++)
			{
				zoneSURFTEMP[i].n = zoneSURFTEMP[i].sx = 0;
				zoneSURFTEMP[i].zonenwert = zone[i].zonenwert;
				zoneSURFTEMP[i].rel = zone[i].rel;
				zoneSURFTEMP[i].RelCellSum = zone[i].RelCellSum;
			}
		}

		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			// for methodes with slides
			zoneINPUT_MASS = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneINPUT_MASS) speichererror(1);
			for (i=0;i<=zonmax;i++)
			{
			  zoneINPUT_MASS[i].n           = zoneINPUT_MASS[i].sx        = 0;
			  zoneINPUT_MASS[i].zonenwert   = zone[i].zonenwert;
			  zoneINPUT_MASS[i].rel         = zone[i].rel;
			  zoneINPUT_MASS[i].RelCellSum  = zone[i].RelCellSum;
			}

			zoneMOBILE_MASS = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneMOBILE_MASS) speichererror(1);
			for (i=0;i<=zonmax;i++)
			{
			  zoneMOBILE_MASS[i].n           = zoneMOBILE_MASS[i].sx        = 0;
			  zoneMOBILE_MASS[i].zonenwert   = zone[i].zonenwert;
			  zoneMOBILE_MASS[i].rel         = zone[i].rel;
			  zoneMOBILE_MASS[i].RelCellSum  = zone[i].RelCellSum;
			}
   

			zoneDEPOSITION = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneDEPOSITION) speichererror(1);
			for (i=0;i<=zonmax;i++)
			{
			  zoneDEPOSITION[i].n           = zoneDEPOSITION[i].sx        = 0;
			  zoneDEPOSITION[i].zonenwert   = zone[i].zonenwert;
			  zoneDEPOSITION[i].rel         = zone[i].rel;
			  zoneDEPOSITION[i].RelCellSum  = zone[i].RelCellSum;
			}
		}

		qschnee = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQSNOW,(int)zonmax,(int)qschneestatcode,(int)4,(char*)qschneestatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		qschnee->schreibzonenwerte("snow_melt_(effective_outflow_from_snow_cover)___",instanzname,LastResults,maxLastResults);

		qsnow_out = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQSNOW_OUT,(int)zonmax,(int)qsnow_outstatcode,(int)4,(char*)qsnow_outstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		qsnow_out->schreibzonenwerte("snow_melt_only(effective_outflow_from_snow_cover)___",instanzname,LastResults,maxLastResults);

		snowstor = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSSNOW,(int)zonmax,(int)snowstorstatcode,(int)3,(char*)snowstorstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		snowstor->schreibzonenwerte("snow_storage_in_mm___",instanzname,LastResults,maxLastResults);

		snowdays = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSNOWDAYS,(int)zonmax,(int)snowdaysstatcode,(int)3,(char*)snowdaysstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		snowdays->schreibzonenwerte("snow_days____",instanzname,LastResults,maxLastResults);

		snow_age = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSNOW_AGE,(int)zonmax,(int)snow_agestatcode,(int)3,(char*)snow_agestatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		snow_age->schreibzonenwerte("snow_age____",instanzname,LastResults,maxLastResults);

		p_snow = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneP_SNOW,(int)zonmax,(int)p_snowstatcode,(int)3,(char*)p_snowstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		p_snow->schreibzonenwerte("p_snow____",instanzname,LastResults,maxLastResults);

		p_rain = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneP_RAIN,(int)zonmax,(int)p_rainstatcode,(int)3,(char*)p_rainstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		p_rain->schreibzonenwerte("p_rain____",instanzname,LastResults,maxLastResults);

		albedo = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneALBEDO,(int)zonmax,(int)albedostatcode,(int)3,(char*)albedostatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		albedo->schreibzonenwerte("albedo____",instanzname,LastResults,maxLastResults);

		// new snow model extensions (Warscher, 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			snowtemp = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols, &zoneSNOWTEMP, (int)zonmax, (int)snowtempstatcode, (int)3, (char*)snowtempstatname, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			snowtemp->schreibzonenwerte("snowtemp____", instanzname, LastResults, maxLastResults);

			surftemp = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols, &zoneSURFTEMP, (int)zonmax, (int)surftempstatcode, (int)3, (char*)surftempstatname, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			surftemp->schreibzonenwerte("surftemp____", instanzname, LastResults, maxLastResults);
		}

		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			input_mass = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols, &zoneINPUT_MASS, (int)zonmax, (int)input_massstatcode, (int)3, (char*)input_massstatname, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			input_mass->schreibzonenwerte("input_mass____",instanzname,LastResults,maxLastResults);

			mobile_mass = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneMOBILE_MASS,(int)zonmax,(int)mobile_massstatcode,(int)3,(char*)mobile_massstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			mobile_mass->schreibzonenwerte("mobile_mass____",instanzname,LastResults,maxLastResults);

			deposition = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneDEPOSITION,(int)zonmax,(int)depositionstatcode,(int)3,(char*)depositionstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			deposition->schreibzonenwerte("deposition____",instanzname,LastResults,maxLastResults);
		}

		// Parameter zur Interpolation des Radiation-Melt-Factors RMF(t) ueber eine Sinusfunktion 
		r_b = (RMFMAX - RMFMIN)/2.0; // Amplitude der Sinusfunktion 
		r_a = RMFMIN + r_b;  // Verschiebung der Sinusfunktion entlang der Ordinate 
		r_c = 0.017214206; // 2*pi/365 Anpassung der Periode der Sinusschwingung 
		r_d = -80.0; // Differenz zwischen dem 1.1. und dem 21.3, wenn sin(c*(d+tag)) == 0 

		Albedo = gridinfo[getindex(_Albedo_)].pointer;

		// is evaporation model is run, the snow evaporation (sublimation) of the new energy balance approaches must be taken into account in evaporation model
		if (ist_vorhanden(_Snow_Evaporation_))
		{
			EVRS=gridinfo[getindex(_Snow_Evaporation_)].pointer;
		}

		char name[500];
		int lay;
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;
		}

		// looking for the ice_firn_grid: if this grid is found, the model has to consider glaciers
		if (ist_vorhanden(_ice_firn_) || ist_vorhanden(_GlacierizedCells_))
		{
			steuerdatei=fopen(steuerdateiname,"r");
			findeabschnitt("[" _ice_firn_ "]");

			fscanfln(steuerdatei,"%d",&glacmethode,"Glacier-Method for ",tstr);

			if (glacmethode > 100)
			{
				glacmethode -= 100;
				DetailedGlacierParams = 1;
			}

			if ((glacmethode!=1)&&(glacmethode!=2)&&(glacmethode!=11)&&(glacmethode!=12)&&(glacmethode!=13))
			{
				(*WSstream) << "FATAL: \n-----------Fatal Error:----------\nMethod for glacier melt must be 1(resp 11 for dynamic model; t-index)\n"
				" or 2 (resp. 12 for dynamic model; t-index with correction by radiation) or 13 (dynamic glaciers with T-Index for ice and firn and the melt method used for snow as in the snow model)"<<"[end-entry]";
				throw(80);
			}

			if (ist_vorhanden(_debris_on_glaciers_))
			{
				GLACDEB=gridinfo[getindex(_debris_on_glaciers_)].pointer;
			}
			else
			{
				GLACDEB=NULL;
			}

			if (glacmethode > 10)
			{
				bFullGridRequired = true; // in mpi_rank 0, the full grid is required for a series of entities in order to do the slide model
				if (mpi_rank == 0)
				{   // that also means, that the rows to read must be extended to the entire grid
					sr = 1;
					er = nrows;
				}

				ICE = gridinfo[getindex(_GlacierizedCells_)].pointer;
				GLACID = gridinfo[getindex(_GlacierCodes_)].pointer;
				// check if all cells of the ICE-Grid have correspondend entries in the GLACID-grid (otherwise these cells would not be part of the grow/erode mechanism)
				// but do this check in rank 0 only since all the other instances will not have the entire grids loaded
				if (mpi_rank == 0)
				{
					for (int i = sr; i <= er; i++)
					{
						for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
						{
							if ((ICE[i][j] > 0) && (GLACID[i][j] == nodata))
							{
								(*WSstream) << "FATAL: \n-----------Fatal Error:----------\nThere is a cell in the GlacierizedCells-grid having no valid value in \n"
									<< _GlacierCodes_ << " (must be an integer ID for each separate glacier) at row " << i << " and col " << j << "[end-entry]";
								throw(80);

							}
						}
					}
					// check if the ICE grid was read in correctly or if it was created and initialized with -9999
					bool icegridisvalid = false;
					for (int i = sr; i <= er; i++)
					{
						for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
						{
							if (ICE[i][j] > 0)
							{
								icegridisvalid = true;
								break; // at least one cell was found with a valid
							}
							if (icegridisvalid) break;
						}
					}
					// try to create the GlacierizedCells grid from ice_firn_grid
					if (!icegridisvalid)
					{
						(*WSstream) << "FATAL: \n-----------Fatal Error:----------\nThere are no glacierized cells in GlacierizedCells-grid (all nodata or equal to 0)!" << "[end-entry]";
						throw(80);


						//					double **ICETEMP=gridinfo[getindex(_ice_firn_)].pointer;
						//					for (int i=1;i<=nrows;i++)
						//					{
						//						for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						//						{
						//							ICE[i][j]=ICETEMP[i][j];
						//							if (ICE[i][j]>1.0) ICE[i][j]=1.0;
						//						}
						//					}
					}
					// check if the GLACID grid was read in correctly or if it was created and initialized with -9999
					bool glacidgridisvalid = false;
					for (int i = sr; i <= er; i++)
					{
						for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
						{
							if ((zn_routing[i][j] != nodata) && (GLACID[i][j] > 0))
							{
								glacidgridisvalid = true;
								break; // at least one cell was found with a valid value
							}
							if (glacidgridisvalid) break;
						}
					}
					if (!glacidgridisvalid)
					{
						(*WSstream) << "FATAL: \n----------- Fatal Error: ----------\nThe " << _GlacierCodes_ << " grid (" << gridinfo[getindex(_GlacierCodes_)].name << ")contains nodata values only. \n"
							<< "Dynamic glacier modelling will not be possible without \n"
							<< "distinct glacier codes for each glacier" << "[end-entry]";
						throw(80);
					}
					// check if each glacier is located in one sub basin only
				}
			}
			else
			{
				ICE=gridinfo[getindex(_ice_firn_)].pointer;
			}

			// lesen der Parameter für Gletschermodell 
			qfirngridname       = newstring(_stringsize_+1,1);
			qfirnstatname       = newstring(_stringsize_+1,1);
			qicegridname        = newstring(_stringsize_+1,1);
			qicestatname        = newstring(_stringsize_+1,1);
			qglacstatname       = newstring(_stringsize_+1,1);
			qsnowonglacstatname = newstring(_stringsize_+1,1);
			snowstorstatname2   = newstring(_stringsize_+1,1);
			SSNOWGLgridname		= newstring(_stringsize_+1,1);
			SLIQGLgridname		= newstring(_stringsize_+1,1);
			SSTOGLgridname		= newstring(_stringsize_+1,1);
			QSNOWGLgridname		= newstring(_stringsize_+1,1);
			SNOWTEMPGLgridname	= newstring(_stringsize_+1,1);
			SURFTEMPGLgridname	= newstring(_stringsize_+1,1);
			EBALGLgridname		= newstring(_stringsize_+1,1);
			EVRSGLgridname		= newstring(_stringsize_+1,1);
			

			if (glacmethode>10) // we need those grids and statistics for dynamic glacier model only
			{
				firnstatname      = newstring(_stringsize_+1,1);
				gmbgridname		  = newstring(_stringsize_+1,1);
				gmbogridname	  = newstring(_stringsize_+1,1);
				firnstackname	  = newstring(_stringsize_+1,1);
				gmbstatname	      = newstring(_stringsize_+1,1);
				gmb2statname	  = newstring(_stringsize_+1,1);

			}

			Qice_els  = (double*)calloc(zonmax_routing+1,sizeof(double));
			Qfirn_els = (double*)calloc(zonmax_routing+1,sizeof(double));
			Qsnow_els = (double*)calloc(zonmax_routing+1,sizeof(double));
			
			Codes = (int*)calloc(zonmax_routing + 1, sizeof(int));
			DDFice = (double*)calloc(zonmax_routing + 1, sizeof(double));
			DDFfirn = (double*)calloc(zonmax_routing + 1, sizeof(double));
			DDFsnow = (double*)calloc(zonmax_routing + 1, sizeof(double));
			MF = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Asnow_a = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Aice_a = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Asnow_b = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Aice_b = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Asnow_akt = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Aice_akt = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Aice_min = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Aice_max = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Asnow_min = (double*)calloc(zonmax_routing + 1, sizeof(double));
			Asnow_max = (double*)calloc(zonmax_routing + 1, sizeof(double));
			kice_els = (double*)calloc(zonmax_routing + 1, sizeof(double));
			kfirn_els = (double*)calloc(zonmax_routing + 1, sizeof(double));
			ksnow_els = (double*)calloc(zonmax_routing + 1, sizeof(double));
			deltafirnper100m = (double*)calloc(zonmax_routing + 1, sizeof(double));
			firn_base = (double*)calloc(zonmax_routing + 1, sizeof(double));

			Qice  = (double*)calloc(zonmax_routing+1,sizeof(double));
			Qfirn = (double*)calloc(zonmax_routing+1,sizeof(double));
			Qsnow = (double*)calloc(zonmax_routing+1,sizeof(double));
			
			//      Qsnow_out = (double*)calloc(zonmax+1,sizeof(double));
			//      Snowdays = (double*)calloc(zonmax+1,sizeof(double));
			//      Snow_age = (double*)calloc(zonmax+1,sizeof(double));
			//      P_snow = (double*)calloc(zonmax+1,sizeof(double));
			//      P_rain = (double*)calloc(zonmax+1,sizeof(double));
			//      ALBEDO = (double*)calloc(zonmax+1,sizeof(double));

			// hier muessen die Parameter eingelesen werden!! 
			if (DetailedGlacierParams == 0)
			{	// legacy method: only one global parameter isread in in first array element
				fscanfln(steuerdatei, "%lf", &(DDFice[0]), "DDFice  for ", tstr);
				fscanfln(steuerdatei, "%lf", &(DDFfirn[0]), "DDFfirn for ", tstr);
				fscanfln(steuerdatei, "%lf", &(DDFsnow[0]), "DDFsnow for ", tstr);
				fscanfln(steuerdatei, "%lf", &(MF[0]), "MF  for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Aice_min[0]), "Aice_min for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Aice_max[0]), "Aice_max for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Asnow_min[0]), "Asnow_min for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Asnow_max[0]), "Asnow_max for ", tstr);
				fscanfln(steuerdatei, "%lf", &(kice_els[0]), "kice_els for ", tstr);
				fscanfln(steuerdatei, "%lf", &(kfirn_els[0]), "kfirn_els ", tstr);
				fscanfln(steuerdatei, "%lf", &(ksnow_els[0]), "ksnow_els for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Qice_els[0]), "initial Qice_els for ", tstr);
				fscanfln(steuerdatei, "%lf", &(Qfirn_els[0]), "initial Qfirn_els ", tstr);
				fscanfln(steuerdatei, "%lf", &(Qsnow_els[0]), "initial Qsnow_els for ", tstr);
				for (i = 1; i <= zonmax_routing; i++)
				{
					DDFice[i] = DDFice[i - 1];
					DDFfirn[i] = DDFfirn[i - 1];
					DDFsnow[i] = DDFsnow[i - 1];
					MF[i] = MF[i - 1];
					Aice_min[i] = Aice_min[i - 1];
					Aice_max[i] = Aice_max[i - 1];
					Asnow_min[i] = Asnow_min[i - 1];
					Asnow_max[i] = Asnow_max[i - 1];
					kice_els[i] = kice_els[i - 1];
					kfirn_els[i] = kfirn_els[i - 1];
					ksnow_els[i] = ksnow_els[i - 1];
					Qice_els[i] = Qice_els[i - 1];
					Qfirn_els[i] = Qfirn_els[i - 1];
					Qsnow_els[i] = Qsnow_els[i - 1];
				}
			}
			else
			{
				// parameters will be read per subbasin (or range iof subbasins like in unsatzonmodel)
				readSubBasinCodes(Codes);
				
				readGlacierParam(DDFice);
				readGlacierParam(DDFfirn);
				readGlacierParam(DDFsnow);
				readGlacierParam(MF);
				readGlacierParam(Aice_min);
				readGlacierParam(Aice_max);
				readGlacierParam(Asnow_min);
				readGlacierParam(Asnow_max);
				readGlacierParam(kice_els);
				readGlacierParam(kfirn_els);
				readGlacierParam(ksnow_els);
				readGlacierParam(Qice_els);
				readGlacierParam(Qfirn_els);
				readGlacierParam(Qsnow_els);
				readGlacierParam(firn_base);
				readGlacierParam(deltafirnper100m);

			}

			fscanfln(steuerdatei,"%s",qfirngridname    ,"name of the firn grid for "     ,tstr);
			fscanfln(steuerdatei,"%d",&qfirnwritegrid  ,"write code for firn grid for "  ,tstr);
			fscanfln(steuerdatei,"%s%d",qfirnstatname,&qfirnstatcode  ,"name of the qfirn-statistic-file for "    ,tstr);

			fscanfln(steuerdatei,"%s",qicegridname    ,"name of the ice grid for "     ,tstr);
			fscanfln(steuerdatei,"%d",&qicewritegrid  ,"write code for ice grid for "  ,tstr);
			fscanfln(steuerdatei,"%s%d",qicestatname,&qicestatcode  ,"name of the qice-statistic-file for "    ,tstr);

			fscanfln(steuerdatei,"%s%d",qglacstatname,&qglacstatcode  ,"name of the qglac-statistic-file for "    ,tstr);
			if (glacmethode>2) // for the new glaciermodel, a specific file for snow runoff from the glacier is implemented
			{
				fscanfln(steuerdatei,"%s%d",qsnowonglacstatname,&qsnowonglacstatcode  ,"name of the qSnowOnGlac-statistic-file for "    ,tstr);
				fscanfln(steuerdatei,"%d",&ice_gridslesen,"read grid parameter for glacier model "    ,tstr);
				if (fscanf(steuerdatei,"%d%lf%lf",&firn_levels,&(firn_base[0]), &(deltafirnper100m[0])) != 3)
				{
					(*WSstream) << "FATAL: \n------------- Fatal Error: -------------\n"
							    <<"Could not read the three parameters firn_levels (int),\n"
								<<"firn_base (double), and deltafirnper100m (double)  \n"
								<<"for dynamic glacier modelling "<<"[end-entry]";
					throw(80);
				}
				if (DetailedGlacierParams == 0)
				{ // otherwise, these parameters were already read in for detailed parameters
					for (i = 1; i <= zonmax_routing; i++)
					{
						firn_base[i] = firn_base[i - 1];
						deltafirnper100m[i] = deltafirnper100m[i - 1];
					}
				}

				tstr=fgets(tstr,_stringsize_,steuerdatei);
				if (fscanf(steuerdatei,"%d%d",&GMBEvalMonth,&GMBEvalDay)!=2)
				{
					(*WSstream) << "FATAL: \n------------- Fatal Error: -------------\n"
							    << "Could not read the two parameters GMBEvalDay (int 1...31)\n"
							    << "and GMBEvalMonth (int 1..12) for dynamic glacier modelling "<<"[end-entry]";
					throw(80);
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				if (fscanf(steuerdatei,"%lf%lf%d%d",&VAscaling,&VAexponent,&numElevBands,&extraWeightFactorBand0)!=4)
				{
					(*WSstream) << "FATAL: \n------------- Fatal Error: -------------\n"
							    << "Could not read the four parameters VAscaling, VAexponent, numElevBands and extraWeightFactorBand0\n"
							    << "for Volume-Area-Relation in dynamic glacier modelling "<<"[end-entry]";
					throw(80);
				}
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				fscanfln(steuerdatei,"%s",firnstackname,"name of the firn stack for ",tstr);
				fscanfln(steuerdatei,"%d",&firnwritestack,"write code for the firn stack in ",tstr);
				fscanfln(steuerdatei,"%s%d",firnstatname,&firnstatcode  ,"name of the firn-statistic-file for "    ,tstr);
				fscanfln(steuerdatei,"%s",gmbgridname,"name of the glacier mass balance grid for "    ,tstr);
				fscanfln(steuerdatei,"%d",&gmbwritegrid,"write code for glacier mass balance grid in ",tstr);
				fscanfln(steuerdatei,"%s",gmbogridname,"name of the glacier mass balance grid (old reference) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&gmbowritegrid,"write code for glacier mass balance grid (old reference) in ",tstr);
				fscanfln(steuerdatei,"%s%d",gmbstatname,&gmbstatcode  ,"name of the glacier mass balance statistics file (1) for "    ,tstr);
				fscanfln(steuerdatei,"%s%d",gmb2statname,&gmb2statcode  ,"name of the glacier mass balance statistics file (2) for "    ,tstr);
//				if (GLACDEB!=NULL)
//				{
					fscanfln(steuerdatei,"%lf",&dGloabDebrisScaling,"global debris scaling factor (requirede when using a debris_on_glaciers grid) ",tstr);
//				}
//				else dGloabDebrisScaling=1.0;
				fscanfln(steuerdatei,"%s%d",snowstorstatname2,&snowstor2statcode,"name of the snow statistics for snow on glacier for ", tstr);
				fscanfln(steuerdatei,"%s",SSNOWGLgridname,"name of the snow storage grid (solid fraction of WE, glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&SSNOWGLwritegrid,"write code for snow storage grid (solid fraction of WE, glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",SLIQGLgridname,"name of the snow storage grid (liquid fraction of WE, glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&SLIQGLwritegrid,"write code for snow storage grid (liquid fraction of WE, glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",SSTOGLgridname,"name of the snow storage grid (total on glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&SSTOGLwritegrid,"write code for snow storage grid (total on glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",QSNOWGLgridname,"name of the snow snow melt grid (incl. rain, total on glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&QSNOWGLwritegrid,"write code for snow melt grid (incl. rain, total on glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",SURFTEMPGLgridname,"name of the snow surface temperature grid (incl. rain, total on glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&SURFTEMPGLwritegrid,"write code for snow surface temperature  (incl. rain, total on glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",SNOWTEMPGLgridname,"name of the snow temperature grid (incl. rain, total on glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&SNOWTEMPGLwritegrid,"write code for snow temperature grid (incl. rain, total on glacierized fraction of cell) in ",tstr);
				fscanfln(steuerdatei,"%s",EVRSGLgridname,"name of the snow evaporation grid (incl. rain, total on glacierized fraction of cell) for "    ,tstr);
				fscanfln(steuerdatei,"%d",&EVRSGLwritegrid,"write code for snow evaporation grid (incl. rain, total on glacierized fraction of cell) in ",tstr);
				// optional parameters for creating glacier cells
				if (fscanf(steuerdatei, "%lf%lf", &maxSnowPackToCreateGlacier, &relSizeNewGlacier) != 2)
				{
					(*WSstream) << "WARNING: \n------------- warning: -------------\n"
						<< "Could not read the two optional parameters maxSnowPackToCreateGlacier (double),\n"
						<< "and relSizeNewGlacier (double) for dynamic glacier modelling that will be used for creating new glaciers\n"
						<< "Setting default values of 1000 and 0.25, i.e. 1000 mm SWE will rmain on a cell that is partly glacierized already, the rest is put into the glacier.\n"
						<< "For cells with no existing glacier, the amount > 1000 will be used to create a glacier that encompasses 25% of the cell" << "[end-entry]";
				}
				tstr = fgets(tstr, _stringsize_, steuerdatei);
			}
			else
			{	
				qsnowonglacstatcode = qglacstatcode;
				
				std::string part1(qglacstatname);
				part1=part1.substr(0,part1.find_last_of("."));
				
				std::string part2(qglacstatname);
				part2=part2.substr(part2.find_last_of("."));

				sprintf(qsnowonglacstatname,"%s_qsgl_%s",part1.c_str(),part2.c_str());

			}

			if ((fclose(steuerdatei))==EOF)
			{
				(*WSstream) << "FATAL: cannot close "<<steuerdateiname<<"[end-entry]";
				throw(11);
			}
			
			
			zoneQFIRN = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
			if (!zoneQFIRN) speichererror(1);
			for (i=0;i<=zonmax_routing;i++)
			{
				zoneQFIRN[i].n           = zoneQFIRN[i].sx        = 0;
				zoneQFIRN[i].zonenwert   = zone_routing[i].zonenwert;
				zoneQFIRN[i].rel         = zone_routing[i].rel;
				zoneQFIRN[i].RelCellSum  = zone_routing[i].RelCellSum;
			}
			
			qfirn = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&zoneQFIRN,(int)zonmax_routing,(int)qfirnstatcode,(int)5,(char*)qfirnstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			qfirn->schreibzonenwerte("glacier_melt_from_firn_areas__",instanzname,LastResults,maxLastResults);

			zoneQICE = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
			if (!zoneQICE) speichererror(1);
			for (i=0;i<=zonmax_routing;i++)
			{
				zoneQICE[i].n           = zoneQICE[i].sx        = 0;
				zoneQICE[i].zonenwert   = zone_routing[i].zonenwert;
				zoneQICE[i].rel         = zone_routing[i].rel;
				zoneQICE[i].RelCellSum  = zone_routing[i].RelCellSum;
			}

			qice = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&zoneQICE,(int)zonmax_routing,(int)qicestatcode,(int)5,(char*)qicestatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			qice->schreibzonenwerte("glacier_melt_from_ice_areas__",instanzname,LastResults,maxLastResults);

			zoneQSnowOnGlac = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
			if (!zoneQSnowOnGlac) speichererror(1);
			for (i=0;i<=zonmax_routing;i++)
			{
				zoneQSnowOnGlac[i].n           = zoneQSnowOnGlac[i].sx        = 0;
				zoneQSnowOnGlac[i].zonenwert   = zone_routing[i].zonenwert;
				zoneQSnowOnGlac[i].rel         = zone_routing[i].rel;
				zoneQSnowOnGlac[i].RelCellSum  = zone_routing[i].RelCellSum;
			}
			
			qsnowonglac = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&zoneQSnowOnGlac,(int)zonmax_routing,(int)qsnowonglacstatcode,(int)5,(char*)qsnowonglacstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			qsnowonglac->schreibzonenwerte("snow_melt_from_ice_areas__",instanzname,LastResults,maxLastResults);

			zoneQGLAC = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
			if (!zoneQGLAC) speichererror(1);
			for (i=0;i<=zonmax_routing;i++)
			{
				zoneQGLAC[i].n           = zoneQGLAC[i].sx        = 0;
				zoneQGLAC[i].zonenwert   = zone_routing[i].zonenwert;
				zoneQGLAC[i].rel         = zone_routing[i].rel;
				zoneQGLAC[i].RelCellSum  = zone_routing[i].RelCellSum;
			}
			
			qglac = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&zoneQGLAC,(int)zonmax_routing,(int)qglacstatcode,(int)5,(char*)qglacstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			qglac->schreibzonenwerte("glacier_melt_total(firn+ice)__",instanzname,LastResults,maxLastResults);

			zoneSSNOW2 = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneSSNOW2) speichererror(1);
			for (i=0;i<=zonmax;i++)
			{
				zoneSSNOW2[i].n           = zoneSSNOW2[i].sx        = 0;
				zoneSSNOW2[i].zonenwert   = zone[i].zonenwert;
				zoneSSNOW2[i].rel         = zone[i].rel;
				zoneSSNOW2[i].RelCellSum  = zone[i].RelCellSum;
			}

//			std::string tmp;
//			tmp=std::string(snowstorstatname);
//			tmp=tmp.substr(0,tmp.find_last_of("."))+"_2"+tmp.substr(tmp.find_last_of("."));
//			sprintf(snowstorstatname2,"%s",tmp.c_str());
			snowstor2 = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
				&zoneSSNOW2,(int)zonmax,(int)snowstor2statcode,(int)3,(char*)snowstorstatname2,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			
			if (retvalue!=0) throw retvalue;
			snowstor2->schreibzonenwerte("snow_storage_in_mm_WITHOUT_GLACIERS_",instanzname,LastResults,maxLastResults);

			// Parameter zur Interpolation des Snow-Melt-Factors Asnow(t) und Ice-Melt-Factors Aice(t) ueber eine Sinusfunktion 
			for (i = 1; i <= zonmax_routing; i++)
			{
				Asnow_b[i] = (Asnow_max[i] - Asnow_min[i]) / 2.0; // Amplitude der Sinusfunktion 
				Asnow_a[i] =  Asnow_min[i] + Asnow_b[i];  // Verschiebung der Sinusfunktion entlang der Ordinate 

				Aice_b[i] = (Aice_max[i] - Aice_min[i]) / 2.0; // Amplitude der Sinusfunktion 
				Aice_a[i] = Aice_min[i] + Aice_b[i];  // Verschiebung der Sinusfunktion entlang der Ordinate 
			}
			
			if ((glacmethode==2)||(glacmethode==12)||(glacmethode==13)) RAD = gridinfo[this->getindex_first_of_3(_Globalstrahlung_,_IEFF_MAX_,NULL)].pointer;

			if (mpi_rank == 0)
			{
				glacierlist = (glacierlisttype*)calloc(zonmax_routing + 1, sizeof(glacierlisttype));
				for (int i = 0; i <= zonmax_routing; i++)
				{
					glacierlist[i].isGlacierized = 0;
					glacierlist[i].CellCntBare = 0;
					glacierlist[i].CellCntGlacerized = 0;
				}

				for (int i = 1; i <= nrows; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if (zn_routing[i][j] > nodata)
						{
							double icevalue = ICE[i][j]; // valid range of icevalue now between 0...2
							if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
							else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
							if (icevalue > 0)
							{
								glacierlist[zarray_routing[(int)zn_routing[i][j]]].isGlacierized = 1;
							}
							glacierlist[zarray_routing[(int)zn_routing[i][j]]].CellCntBare += (1.0 - icevalue);
							glacierlist[zarray_routing[(int)zn_routing[i][j]]].CellCntGlacerized += icevalue;
						}
					}
				}
			}

#ifdef _MPI
			// bcastig center_of_mass to all other ranks
			MPI_Barrier(MPI_COMM_WORLD);
			const int nitems = 3;
			int blocklengths[3] = { 1, 1, 1};
			MPI_Datatype types[3] = { MPI_C_BOOL, MPI_DOUBLE, MPI_DOUBLE };
			MPI_Aint offsets[3] = { offsetof(glacierlisttype, isGlacierized),
				offsetof(glacierlisttype, CellCntGlacerized),
				offsetof(glacierlisttype, CellCntBare) };

			MPI_Datatype MPIglacierlisttype;
			MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPIglacierlisttype);
			MPI_Type_commit(&MPIglacierlisttype);

			if (mpi_rank > 0)
			{
				glacierlist = (glacierlisttype*)calloc(zonmax_routing + 1, sizeof(glacierlisttype));
				if (!glacierlist) speichererror(1);
			}
			MPI_Bcast(&(glacierlist[0]), zonmax_routing + 1, MPIglacierlisttype, 0,MPI_COMM_WORLD);

			MPI_Type_free(&MPIglacierlisttype);
#endif

//for (int i=0;i<=zonmax_routing;i++)
//{
//	cout << "DEBUG: glacierlist["<<i<<"].isGlacierized: "<<glacierlist[i].isGlacierized<<", glacierlist["<<i<<"].CellCntBare: "<<glacierlist[i].CellCntBare<<", glacierlist[i].CellCntGlacerized: "<<glacierlist[i].CellCntGlacerized<<" in rank "<<mpi_rank<<endl;
//}

			// statistics for dynamic glacier model only
			if (glacmethode>10)
			{
				// todo: create statistics for Firn and GWBalance
				zoneFIRN = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
				if (!zoneFIRN) speichererror(1);
				for (int i=0;i<=zonmax_routing;i++)
				{
					zoneFIRN[i].n           = zoneFIRN[i].sx        = 0;
					zoneFIRN[i].zonenwert   = zone_routing[i].zonenwert;
					zoneFIRN[i].rel         = zone_routing[i].rel;
					zoneFIRN[i].RelCellSum  = zone_routing[i].RelCellSum;
				}
				
				firnstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
				&zoneFIRN,(int)zonmax_routing,(int)firnstatcode,(int)4,(char*)firnstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
				if (retvalue!=0) throw retvalue;
				firnstat->schreibzonenwerte("WE_in_firn_[mm_per_subbasin]",instanzname,LastResults,maxLastResults);

				zoneGMB = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
				if (!zoneGMB) speichererror(1);
				for (int i=0;i<=zonmax_routing;i++)
				{
					zoneGMB[i].n           = zoneGMB[i].sx        = 0;
					zoneGMB[i].zonenwert   = zone_routing[i].zonenwert;
					zoneGMB[i].rel         = zone_routing[i].rel;
					zoneGMB[i].RelCellSum  = zone_routing[i].RelCellSum;
				}
				
				gmbstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
				&zoneGMB,(int)zonmax_routing,(int)gmbstatcode,(int)5,(char*)gmbstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
				if (retvalue!=0) throw retvalue;
				gmbstat->schreibzonenwerte("glacier_mass_balance_[mm_per_subbasin]",instanzname,LastResults,maxLastResults);

				zoneGMB2 = (statrecord *)calloc(zonmax_routing+1,sizeof(statrecord));
				if (!zoneGMB2) speichererror(1);
				for (int i=0;i<=zonmax_routing;i++)
				{
					zoneGMB2[i].n           = zoneGMB2[i].sx        = 0;
					zoneGMB2[i].zonenwert   = zone_routing[i].zonenwert;
					zoneGMB2[i].rel         = zone_routing[i].rel;
					zoneGMB2[i].RelCellSum  = zone_routing[i].RelCellSum;
				}
				
				gmb2stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
				&zoneGMB2,(int)zonmax_routing,(int)gmb2statcode,(int)5,(char*)gmb2statname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
				if (retvalue!=0) throw retvalue;
				gmb2stat->schreibzonenwerte("glacier_mass_balance_[mm_GlacierPartOfSubbasin]",instanzname,LastResults,maxLastResults);
			}
		}
		else
		{
			ICE=NULL;
			FIRN=NULL;
			glacierlist=NULL;
		}

		// $kf 2016-08-22: If the snow canopy model is applied, read modified meteorology for the sub-canopy.
		// If there is no canopy, the standard values will be applied instead. 
		// see InterceptionBeforeSnow function
		if (methode == 4) E = gridinfo[getindex(_Dampfdruck_)].pointer; //$kf: ?? Dow we need to re-calculate vapor pressure from relative humidity?

// since some grids may be not known when snow model is initialized, the following grids will be assigned only after interception model was created
//		if ((methode != 1) && (methode != 8) && (methode != 9) && (methode != 10))
//		{
//			U = gridinfo[getindex_first_of_3(_Snow_Canopy_Wind_, _Windgeschwindigkeit_, _ungueltiger_Name_)].pointer; // needed for all but T-Index methods
//		}
//
//		T = gridinfo[getindex_first_of_4(_Snow_Canopy_Temperature_,_Korrigierte_Temperatur_, _Temperatur_, _Temperatur_14_)].pointer;

//		P    = gridinfo[getindex(_Niederschlag_)].pointer; must be initialized during the first time step since the initialization here may need the interception output and this sub model is initialized after the snow model


		// Die vom Expositionskorrekturmodul duplizierte und vom Verdunstungsmodell modifizierten Tag/Nacht-
		// Temperaturmittelwerte werden auch von der Schneeschmelzroutine bei Intervallaenge Tageswert genutzt,
		// um eine Unterscheidung in Tag und Nacht durchzufuehren 

		if (InV>720.0)
		Tnacht = gridinfo[getindex(_Nacht_Temperatur_)].pointer;
		else Tnacht = T;

		// suchen des ersten freien Eintrages 
		sti=1;
		while (gridinfo[sti].identification!=NULL)
		{
			sti ++;
			if (sti==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid left. Please contact the author to set the number of available grids to a larger number "<<gridmaxanzahl <<"[end-entry]";
				throw(12);
			}
		}

		// new snow model extensions (Warscher 2013)
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			SLIDEFRACTION1 = gridinfo[getindex(_slidefraction1_)].pointer;
			SLIDEFRACTION2 = gridinfo[getindex(_slidefraction2_)].pointer;
			SLIDEFRACTION3 = gridinfo[getindex(_slidefraction3_)].pointer;
			SLIDEFRACTION4 = gridinfo[getindex(_slidefraction4_)].pointer;
			// elevation order grids will be created internally to make sure that the correct basin delineation is used (to avoid wrong indirections to cells outside the basin)
			//			ELEV_ORDER_ROW = gridinfo[getindex(_elev_order_row_)].pointer;
			//			ELEV_ORDER_COL = gridinfo[getindex(_elev_order_col_)].pointer;
			SL = gridinfo[getindex(_Gefaelle_)].pointer;
		}

		if ((methode == 7) || (methode == 9) || (methode == 10) || (methode == 11) || (methode == 13) || (methode == 14)) // wind resistribution
		{
			DEPOINDEX = gridinfo[getindex(_depoindex_)].pointer; // wind redistribution grid, nothing to do with slides deposition!
		}

		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			CL = gridinfo[getindex(_Cloudiness_)].pointer;
//$kf 2016-08-22: moved to setInterceptionBeforeSnow
// since some grids may be not known when snow model is initialized, the following grids will be assigned only after interception model was created
//			RLF = gridinfo[getindex_first_of_3(_Snow_Canopy_RH_, _Luftfeuchte_, _ungueltiger_Name_)].pointer;
//			RGEX = gridinfo[getindex_first_of_3(_Snow_Canopy_RAD_, _AspectCorrectedGlobalRadiation_, _ungueltiger_Name_)].pointer;
			if (ist_vorhanden(_SurfaceTemperature_))
			{
				SURFTEMPINTERPOL = gridinfo[getindex(_SurfaceTemperature_)].pointer;
			}
		}

		if (gridslesen != 1)
		{
			// neu erstellen, also nur Speicher reservieren 
			if (!ist_vorhanden(_snowstorage_solid_))
			{
				SSNOW = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(sti, SSNOWgridname, SSNOW, NULL, 0, _snowstorage_solid_, SSNOWwritegrid, sr, er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SSNOW[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				ssnowindex = sti;
				sti++;
			}
			else
			{
				ssnowindex = getindex(_snowstorage_solid_);
				SSNOW = gridinfo[ssnowindex].pointer;
			}


			if (!ist_vorhanden(_snowstorage_liquid_))
			{
				SLIQ = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(sti, SLIQgridname, SLIQ, NULL, 0, _snowstorage_liquid_, SLIQwritegrid, sr, er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SLIQ[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				sliqindex = sti;
				sti++;
			}
			else
			{
				sliqindex = getindex(_snowstorage_liquid_);
				SLIQ = gridinfo[sliqindex].pointer;
			}

			if (!ist_vorhanden(_Schneespeicher_Total_))
			{
				SSTO = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(sti, SSTOgridname, SSTO, NULL, 0, _Schneespeicher_Total_, SSTOwritegrid, sr,er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SSTO[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				sstoindex = sti;
				sti++;
			}
			else
			{
				sstoindex = getindex(_Schneespeicher_Total_);
				SSTO = gridinfo[sstoindex].pointer;
			}

			if (!ist_vorhanden(_Schneetage_))
			{
				SNOWDAYS = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(sti, SNOWDAYSgridname, SNOWDAYS, NULL, 0, _Schneetage_, SNOWDAYSwritegrid, sr,er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SNOWDAYS[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				snowdaysindex = sti;
				sti++;
			}
			else
			{
				snowdaysindex = getindex(_Schneetage_);
				SNOWDAYS = gridinfo[snowdaysindex].pointer;
			}

			if (!ist_vorhanden(_Schneealter_))
			{
				SNOW_AGE = matrix(sr,er,1,ncols,nrows);
				initialize_grid(sti, SNOW_AGEgridname, SNOW_AGE, NULL, 0, _Schneealter_, SNOW_AGEwritegrid, sr,er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SNOW_AGE[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				snow_ageindex = sti;
				sti++;
			}
			else
			{
				snow_ageindex = getindex(_Schneealter_);
				SNOW_AGE = gridinfo[snow_ageindex].pointer;
			}

			// new snow model extensions (Warscher 2013)
			if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // enbal
			{
				if (!ist_vorhanden(_snowtemp_))
				{
					SNOWTEMP = matrix(sr,er,1,ncols,nrows);
					initialize_grid(sti, SNOWTEMPgridname, SNOWTEMP, NULL, 0, _snowtemp_, SNOWTEMPwritegrid, sr,er);
					for (int k = sr; k <= er; k++)
					{
						for (int j = 1; j <= ncols; j++)
						{
							SNOWTEMP[k][j] = nodata;
						}
					}
					snowtempindex = sti;
					sti++;
				}
				else
				{
					snowtempindex = getindex(_snowtemp_);
					SNOWTEMP = gridinfo[snowtempindex].pointer;
				}

				if (!ist_vorhanden(_surftemp_))
				{
					SURFTEMP = matrix(sr, er, 1, ncols, nrows);
					initialize_grid(sti, SURFTEMPgridname, SURFTEMP, NULL, 0, _surftemp_, SURFTEMPwritegrid, sr, er);
					for (int k = sr; k <= er; k++)
					{
						for (int j = 1; j <= ncols; j++)
						{
							SURFTEMP[k][j] = nodata;
						}
					}
					surftempindex = sti;
					sti++;
				}
				else
				{
					surftempindex = getindex(_surftemp_);
					SURFTEMP = gridinfo[surftempindex].pointer;
				}

// not yet read from a file, so kept uotside of this section yet (see after nect else branch)
//				if (!ist_vorhanden(_EnergyBalance_))
//				{
//					EBAL = matrix(sr, er, 1, ncols, nrows);
//					sprintf(EBALgridname, "%s", _EnergyBalance_);
//					initialize_grid(sti, EBALgridname, EBAL, NULL, 0, _EnergyBalance_, EBALwritegrid, sr, er);
//					for (int k = sr; k <= er; k++)
//					{
//						for (int j = 1; j <= ncols; j++)
//						{
//							EBAL[k][j] = 0;
//						}
//					}
//					ebalindex = sti;
//					sti++;
//				}
//				else
//				{
//					ebalindex = getindex(_EnergyBalance_);
//					EBAL = gridinfo[ebalindex].pointer;
//				}

			}
		}
		else
		{
			// alte Dateien einlesen 
			// for snow storages: if a method with slides is choosen, those grids need to be read completely in in rank 0 (and fractionally in other ranks)
			if (!ist_vorhanden(_snowstorage_solid_))
			{
				readgridfromfile(SSNOWgridname, &SSNOW, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(sti, SSNOWgridname, SSNOW, NULL, 0, _snowstorage_solid_, SSNOWwritegrid, sr,er);
				ssnowindex = sti;
				sti++;
			}
			else
			{
				ssnowindex = getindex(_snowstorage_solid_);
				SSNOW = gridinfo[ssnowindex].pointer;
			}

			if (!ist_vorhanden(_snowstorage_liquid_))
			{
				readgridfromfile(SLIQgridname, &SLIQ, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(sti, SLIQgridname, SLIQ, NULL, 0, _snowstorage_liquid_, SLIQwritegrid, sr,er);
				sliqindex = sti;
				sti++;
			}
			else
			{
				sliqindex = getindex(_snowstorage_liquid_);
				SLIQ = gridinfo[sliqindex].pointer;
			}

			if (!ist_vorhanden(_Schneespeicher_Total_))
			{
				readgridfromfile(SSTOgridname, &SSTO, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(sti, SSTOgridname, SSTO, NULL, 0, _Schneespeicher_Total_, SSTOwritegrid, sr, er);
				sstoindex = sti;
				sti++;
			}
			else
			{
				sstoindex = getindex(_Schneespeicher_Total_);
				SSTO = gridinfo[sstoindex].pointer;
			}

			if (!ist_vorhanden(_Schneetage_))
			{
				readgridfromfile(SNOWDAYSgridname, &SNOWDAYS, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(sti, SNOWDAYSgridname, SNOWDAYS, NULL, 0, _Schneetage_, SNOWDAYSwritegrid, mpistartrow, mpiendrow);
				snowdaysindex = sti;
				sti++;
			}
			else
			{
				snowdaysindex = getindex(_Schneetage_);
				SNOWDAYS = gridinfo[snowdaysindex].pointer;
			}

			if (!ist_vorhanden(_Schneealter_))
			{
				readgridfromfile(SNOW_AGEgridname, &SNOW_AGE, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(sti, SNOW_AGEgridname, SNOW_AGE, NULL, 0, _Schneealter_, SNOW_AGEwritegrid, mpistartrow, mpiendrow);
				snow_ageindex = sti;
				sti++;
			}
			else
			{
				snow_ageindex = getindex(_Schneealter_);
				SNOW_AGE = gridinfo[snow_ageindex].pointer;
			}

			// new snow model extensions (Warscher 2013)
			if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // enbal
			{
				if (!ist_vorhanden(_snowtemp_))
				{
					readgridfromfile(SNOWTEMPgridname, &SNOWTEMP, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 273.16, _LimitAutoFillOff_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, SNOWTEMPgridname, SNOWTEMP, NULL, 0, _snowtemp_, SNOWTEMPwritegrid, mpistartrow, mpiendrow);
					snowtempindex = sti;
					sti++;
				}
				else
				{
					snowtempindex = getindex(_snowtemp_);
					SNOWTEMP = gridinfo[snowtempindex].pointer;
				}

				if (!ist_vorhanden(_surftemp_))
				{
					readgridfromfile(SURFTEMPgridname, &SURFTEMP, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 273.16, _LimitAutoFillOff_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, SURFTEMPgridname, SURFTEMP, NULL, 0, _surftemp_, SURFTEMPwritegrid, mpistartrow, mpiendrow);
					surftempindex = sti;
					sti++;
				}
				else
				{
					surftempindex = getindex(_surftemp_);
					SURFTEMP = gridinfo[surftempindex].pointer;
				}

// not yet read from a file
//				if (!ist_vorhanden(_EnergyBalance_))
//				{
//					sprintf(EBALgridname, "%s", _EnergyBalance_);
//					readgridfromfile(EBALgridname, &EBAL, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOff_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
//					initialize_grid(sti, EBALgridname, EBAL, NULL, 0, _EnergyBalance_, EBALwritegrid, mpistartrow, mpiendrow);
//					ebalindex = sti;
//					sti++;
//				}
//				else
//				{
//					ebalindex = getindex(_EnergyBalance_);
//					EBAL = gridinfo[ebalindex].pointer;
//				}
			}
		}

		// new snow model extensions (Warscher, 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // energy balance methods
		{
			if (!ist_vorhanden(_EnergyBalance_))
			{
				EBAL = matrix(sr, er, 1, ncols, nrows);
				sprintf(EBALgridname, "%s", _EnergyBalance_);
				initialize_grid(sti, EBALgridname, EBAL, NULL, 0, _EnergyBalance_, EBALwritegrid, sr, er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						EBAL[k][j] = 0;
					}
				}
				ebalindex = sti;
				sti++;
			}
			else
			{
				ebalindex = getindex(_EnergyBalance_);
				EBAL = gridinfo[ebalindex].pointer;
			}
		}

		// new snow model extensions (Warscher, 2013)
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			if (mpi_rank == 0) // slides will be modeled in rank 0 only
			{
				if (!ist_vorhanden(_max_deposition_))
				{
					DEPOSITION_MAX = matrix(1, nrows, 1, ncols);
					initialize_grid(sti, "", DEPOSITION_MAX, NULL, 0, _max_deposition_, 0, 1, nrows);
					for (int k = 1; k <= nrows; k++)
					{
						for (int j = 1; j <= ncols; j++)
						{
							if (zn[k][j] == nodata)
							{
								DEPOSITION_MAX[k][j] = nodata;
							}
							else
							{
								if (SL[k][j] >= maxDepositionSlope)
								{
									DEPOSITION_MAX[k][j]=0;
								}
								else
								{
									DEPOSITION_MAX[k][j] = zn[k][j] == nodata ? nodata : (1 - SL[k][j] / maxDepositionSlope)*maxDepositionFactor;
								}
								if (DEPOSITION_MAX[k][j] < 0) DEPOSITION_MAX[k][j] = 0;
							}
						}
					}
					depositionmaxindex = sti;
					sti++;
				}
				else
				{
					depositionmaxindex = getindex(_max_deposition_);
					DEPOSITION_MAX = gridinfo[depositionmaxindex].pointer;
				}

				ELEV_ORDER_ROW = matrix(1, nrows, 1, ncols);
				initialize_grid(sti, "", ELEV_ORDER_ROW, NULL, 0, _elev_order_row_, 0, 1, nrows);
				elevorderrowindex = sti;
				sti++;

				ELEV_ORDER_COL = matrix(1, nrows, 1, ncols);
				initialize_grid(sti, "", ELEV_ORDER_COL, NULL, 0, _elev_order_col_, 0, 1, nrows);
				elevordercolindex = sti;
				sti++;

				// now calculating the ordered vector of elevations and write it back to the two newly created grids
				std::list<tupel> Liste;
				for (int k = 1; k <= nrows; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						if (zn[k][j] != nodata)
						{
							Liste.push_back(tupel(k, j, DHM[k][j]));
						}
					}
				}
				Liste.sort();
				long cnt = 0;
				std::list<tupel>::iterator iter;
				iter = Liste.end();
				// now writing the sorted coordinates back to the EOC and EOR grids
				// EOC = Elevation Order Columns
				// EOR = Elevation Order Rows
				// this is used in WaSiM to start the slides at the highest possible point to get the amount of moving masses correctly
				for (int k = 1; k <= nrows; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						if (zn[k][j] != nodata)
						{
							iter--;
							ELEV_ORDER_COL[k][j] = iter->j;
							ELEV_ORDER_ROW[k][j] = iter->i;
						}
						else
						{
							ELEV_ORDER_COL[k][j] = nodata;
							ELEV_ORDER_ROW[k][j] = nodata;
						}
					}
				}
			}
			INPUT_MASS = matrix(sr, er, 1, ncols,nrows);
			initialize_grid(sti, INPUT_MASSgridname, INPUT_MASS, NULL, 0, _input_mass_, INPUT_MASSwritegrid, sr,er);
			for (int k = sr; k <= er; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					INPUT_MASS[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			input_massindex = sti;
			sti++;

			MOBILE_MASS = matrix(sr, er, 1, ncols, nrows);
			initialize_grid(sti, MOBILE_MASSgridname, MOBILE_MASS, NULL, 0, _mobile_mass_, MOBILE_MASSwritegrid, sr, er);
			for (int k = sr; k <= er; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					MOBILE_MASS[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			mobile_massindex = sti;
			sti++;

			DEPOSITION = matrix(sr, er, 1, ncols, nrows);
			initialize_grid(sti, DEPOSITIONgridname, DEPOSITION, NULL, 0, _deposition_, DEPOSITIONwritegrid, sr, er);
			for (int k = sr; k <= er; k++)
			{
				for (j = 1; j <= ncols; j++)
				{
					DEPOSITION[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			depositionindex = sti;
			sti++;

		}

		P_SNOW = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
		initialize_grid(sti, P_SNOWgridname, P_SNOW, NULL, 0, _Schneeanteil_, P_SNOWwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				P_SNOW[k][j] = zn[k][j] == nodata ? nodata : 0;
			}
		}
		p_snowindex = sti;
		sti ++;

		P_RAIN = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(sti, P_RAINgridname, P_RAIN, NULL, 0, _Regenanteil_, P_RAINwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				P_RAIN[k][j] = zn[k][j] == nodata ? nodata : 0;
			}
		}
		p_rainindex = sti;
		sti++;

		// Reservieren von Platz fuer Schneedeckenausfluss erfolgt immer ohne Dateieinlesen 
		QSNOW = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(sti, qschneegridname, QSNOW, NULL, 0, _Schneeausfluss_, Qwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				QSNOW[k][j] = zn[k][j] == nodata ? nodata : 0;
			}
		}
		qsnowindex = sti;
		sti++;

		// Reservieren von Platz fuer Schneedeckenausfluss (ohne Regen) 
		QSNOW_OUT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(sti, qsnow_outgridname, QSNOW_OUT, NULL, 0, _Snow_out_, qsnow_outwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				QSNOW_OUT[k][j] = zn[k][j] == nodata ? nodata : 0;
			}
		}
		qsnow_outindex = sti;
		sti++;

		if (ICE != NULL)
		{
			// for all glacier grids: sr and er should be set depending on the methode and glac_method before
			// runoff grids doesn't need to be allocated for the entire grid in rank 0, only thise grids that are part of the resizing algorithm
			QFIRN = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
			initialize_grid(sti, qfirngridname, QFIRN, NULL, 0, _melt_from_firn_, qfirnwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					QFIRN[k][j] = zn_routing[k][j] == nodata ? nodata : 0;
				}
			}
			qfirnindex = sti;
			sti++;

			QICE = matrix(mpistartrow,mpiendrow, 1, ncols,nrows);
			initialize_grid(sti, qicegridname, QICE, NULL, 0, _melt_from_ice_, qicewritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					QICE[k][j] = zn_routing[k][j] == nodata ? nodata : 0;
				}
			}
			qiceindex = sti;
			sti++;

			if (glacmethode>10)
			{
				if (ice_gridslesen !=1)
				{
					FIRN  = drei_d_matrix(sr,er,1,ncols,firn_levels,nrows);
					FIRN1 = matrix(sr, er, 1, ncols,nrows);
					initialize_grid(sti, firnstackname, FIRN1, FIRN, firn_levels, _firn_WE_stack_, firnwritestack, sr, er);
					for (int k=sr;k<=er;k++)
					{
						for (int j=1;j<=ncols;j++)
						{
							if ((zn_routing[k][j]!=nodata)&&(ICE[k][j]>0.0))
							{
								FIRN1[k][j]=FIRN[k][j][0]=0;
								for (int lv=1;lv<=firn_levels;lv++)
								{
									FIRN[k][j][lv] = (DHM[k][j]-firn_base[zarray_routing[(int)zn_routing[k][j]]])*deltafirnper100m[zarray_routing[(int)zn_routing[k][j]]];
									if (FIRN[k][j][lv]<0.0)
									FIRN[k][j][lv]=0; // only above firn_base there should be any firn...
									FIRN[k][j][0] += FIRN[k][j][lv];
									FIRN1[k][j] = FIRN[k][j][0];
								}
							}
							else
							{
								for (int lv=0;lv<=firn_levels;lv++)
								{
									FIRN[k][j][lv]=nodata;
								}
								FIRN1[k][j]=nodata;
							}
						}
					}
					firnindex = sti;
					sti++;

					// Glacier Mass Balance in mm
					GMB = matrix(sr,er,1,ncols,nrows);
					initialize_grid(sti, gmbgridname, GMB, NULL, 0, _GlacierMassBalance_, gmbwritegrid, sr, er);
					for (int k=sr;k<=er;k++)
					{
						for (int j=1;j<=ncols;j++)
						{
							GMB[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
						}
					}
					gmbindex = sti;
					sti++;

					// Glacier Mass Balance for the year before (old reference) in mm
					GMBO = matrix(sr, er, 1, ncols,nrows);
					initialize_grid(sti, gmbogridname, GMBO, NULL, 0, _OldGlacierMassBalance_, gmbowritegrid, sr, er);
					for (int k=sr;k<=er;k++)
					{
						for (int j=1;j<=ncols;j++)
						{
							GMBO[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
						}
					}
					gmboindex = sti;
					sti++;


					// new 2014-12-21: correction of snow balance bug:
					// three snow storages will be created that contain solid, liquid and total snow storage for the glacierized fraction of cells
					SSNOWGL = matrix(sr, er, 1, ncols,nrows);
					initialize_grid(sti, SSNOWGLgridname, SSNOWGL, NULL, 0, _SolidSnowStorage_onGlacier_, SSNOWGLwritegrid, sr, er);
					for (int k=sr;k<=er;k++)	
					{
						for (int j=1;j<=ncols;j++)	
						{
							SSNOWGL[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
						}
					}
					ssnowglindex = sti;
					sti++;

					SLIQGL = matrix(sr, er, 1, ncols,nrows);
					initialize_grid(sti, SLIQGLgridname, SLIQGL, NULL, 0, _LiquidSnowStorage_onGlacier_, SLIQGLwritegrid, sr, er);
					for (int k=sr;k<=er;k++)	
					{
						for (int j=1;j<=ncols;j++)	
						{
							SLIQGL[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
						}
					}
					sliqglindex = sti;
					sti++;

					SSTOGL = matrix(sr, er, 1, ncols,nrows);
					initialize_grid(sti, SSTOGLgridname, SSTOGL, NULL, 0, _TotalSnowStorage_onGlacier_, SSTOGLwritegrid, sr, er);
					for (int k=sr;k<=er;k++)	
					{
						for (int j=1;j<=ncols;j++)	
						{
							SSTOGL[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
						}
					}
					sstoglindex = sti;
					sti++;

					if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
					{

						if (!ist_vorhanden(_surftemp_onGlac_))
						{
							SURFTEMPGL = matrix(sr, er, 1, ncols, nrows);
							initialize_grid(sti, SURFTEMPGLgridname, SURFTEMPGL, NULL, 0, _surftemp_onGlac_, SURFTEMPGLwritegrid, sr, er);
							for (int k = sr; k <= er; k++)
							{
								for (int j = 1; j <= ncols; j++)
								{
									SURFTEMPGL[k][j] = (zn_routing[k][j] == nodata) ? nodata : 0;
								}
							}
							surftempglindex = sti;
							sti++;
						}
						else
						{
							surftempglindex = getindex(_surftemp_onGlac_);
							SURFTEMPGL = gridinfo[surftempglindex].pointer;
						}

						if (!ist_vorhanden(_snowtemp_onGlac_))
						{
							SNOWTEMPGL = matrix(sr, er, 1, ncols, nrows);
							initialize_grid(sti, SNOWTEMPGLgridname, SNOWTEMPGL, NULL, 0, _snowtemp_onGlac_, SNOWTEMPGLwritegrid, sr, er);
							for (int k = sr; k <= er; k++)
							{
								for (int j = 1; j <= ncols; j++)
								{
									SNOWTEMPGL[k][j] = (zn_routing[k][j] == nodata) ? nodata : 0;
								}
							}
							snowtempglindex = sti;
							sti++;
						}
						else
						{
							snowtempglindex = getindex(_snowtemp_onGlac_);
							SNOWTEMPGL = gridinfo[snowtempglindex].pointer;
						}

// not yet a grid that can be read in from file
//						if (!ist_vorhanden(_EnergyBalance_onGlac_))
//						{
//							EBALGL = matrix(sr, er, 1, ncols, nrows);
//							sprintf(EBALGLgridname, "%s", _EnergyBalance_onGlac_)
//							initialize_grid(sti, EBALGLgridname, EBALGL, NULL, 0, _EnergyBalance_onGlac_, EBALGLwritegrid, sr, er);
//							for (int k = sr; k <= er; k++)
//							{
//								for (int j = 1; j <= ncols; j++)
//								{
//									EBALGL[k][j] = (zn_routing[k][j] == nodata) ? nodata : 0;
//								}
//							}
//							ebalglindex = sti;
//							sti++;
//						}
//						else
//						{
//							ebalglindex = getindex(_EnergyBalance_onGlac_);
//							EBALGL = gridinfo[ebalglindex].pointer;
//						}
					}
				}
				else
				{
					firn_levels = readstackfromfile(firnstackname, &FIRN, bFullGridRequired); // firn_levels is the number of levels from the stack to be read in
					FIRN1=matrix(sr,er,1,ncols,nrows);
					initialize_grid(sti,firnstackname,FIRN1,FIRN,firn_levels,_firn_WE_stack_,firnwritestack,sr,er); // because theta is not a grid but a stack, the adress-parameter for the grid is NULL 
					firnindex=sti;
					for (int k=mpistartrow;k<=mpiendrow;k++)
					{
						for (int j=1;j<=ncols;j++)
						{
							FIRN1[k][j]=FIRN[k][j][0];
						}
					}
					sti++;

					// read GMB Glacier Mass Balance grid
					readgridfromfile(gmbgridname, &GMB, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, gmbgridname, GMB, NULL, 0, _GlacierMassBalance_, gmbwritegrid, sr, er);
					gmbindex = sti;
					sti++;

					// read GMBO grid (Glacier Mass Balance for the last main balanceing day on October 1st)
					readgridfromfile(gmbogridname, &GMBO, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, gmbogridname, GMBO, NULL, 0, _OldGlacierMassBalance_, gmbowritegrid, sr,er);
					gmboindex = sti;
					// make sure no cell neither in GMB nor in GMBO has nodata values where the other is valid -> replace with other value automatically
					for (int k=sr;k<=er;k++)
					{
						for (int j=1;j<=ncols;j++)
						{
							if ((GMB[k][j]==nodata)&&(GMBO[k][j]!=nodata))	{
								GMB[k][j]=GMBO[k][j];
							} else if ((GMBO[k][j]==nodata)&&(GMB[k][j]!=nodata))	{
								GMBO[k][j]=GMB[k][j];
							}
						}
					}
					sti++;

					// new 2014-12-21: correction of snow balance bug:
					// three snow storages will be read on from disk that contaoin sloid, liquid and total snow storage for the glacierized fraction of cells
					readgridfromfile(SSNOWGLgridname, &SSNOWGL, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, SSNOWGLgridname, SSNOWGL, NULL, 0, _SolidSnowStorage_onGlacier_, SSNOWGLwritegrid, sr,er);
					ssnowglindex = sti;
					sti ++;

					readgridfromfile(SLIQGLgridname, &SLIQGL, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, SLIQGLgridname, SLIQGL, NULL, 0, _LiquidSnowStorage_onGlacier_, SLIQGLwritegrid, sr,er);
					sliqglindex = sti;
					sti ++;

					readgridfromfile(SSTOGLgridname, &SSTOGL, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, bFullGridRequired, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
					initialize_grid(sti, SSTOGLgridname, SSTOGL, NULL, 0, _TotalSnowStorage_onGlacier_, SSTOGLwritegrid, sr, er);
					sstoglindex = sti;
					sti++;

					if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
					{
						readgridfromfile(SURFTEMPGLgridname, &SURFTEMPGL, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, 273.16, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
						initialize_grid(sti, SURFTEMPGLgridname, SURFTEMPGL, NULL, 0, _surftemp_onGlac_, SURFTEMPGLwritegrid, mpistartrow, mpiendrow);
						surftempglindex = sti;
						sti++;

						readgridfromfile(SNOWTEMPGLgridname, &SNOWTEMPGL, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, 273.16, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
						initialize_grid(sti, SNOWTEMPGLgridname, SNOWTEMPGL, NULL, 0, _snowtemp_onGlac_, SNOWTEMPGLwritegrid, mpistartrow, mpiendrow);
						snowtempglindex = sti;
						sti++;

// not a grid that can be read in from a file yet.
//						readgridfromfile(EBALGLgridname, &EBALGL, 0, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
//						initialize_grid(sti, EBALGLgridname, EBALGL, NULL, 0, _EnergyBalance_onGlac_, EBALGLwritegrid, mpistartrow, mpiendrow);
//						snowtempglindex = sti;
//						sti++;


					}

				}


				if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
				{
					if (!ist_vorhanden(_EnergyBalance_onGlac_))
					{
						EBALGL = matrix(sr, er, 1, ncols, nrows);
						sprintf(EBALGLgridname, "%s", _EnergyBalance_onGlac_);
						initialize_grid(sti, EBALGLgridname, EBALGL, NULL, 0, _EnergyBalance_onGlac_, EBALGLwritegrid, sr, er);
						for (int k = sr; k <= er; k++)
						{
							for (int j = 1; j <= ncols; j++)
							{
								EBALGL[k][j] = (zn_routing[k][j] == nodata) ? nodata : 0;
							}
						}
						ebalglindex = sti;
						sti++;
					}
					else
					{
						ebalglindex = getindex(_EnergyBalance_onGlac_);
						EBALGL = gridinfo[ebalglindex].pointer;
					}
				}


				// for keeping the glacier mass balance even when resizing the glacier, 
				// the old glacierized fraction and old glacier-ID must be stored for calculating the balance to be shifted to the new glacier cells after resizing.
				ICEOLD = matrix(sr,er,1,ncols,nrows);
				initialize_grid(sti, _GlacierizedCellsOld_, ICEOLD, NULL, 0, _GlacierizedCellsOld_, 0, sr,er);
				for (int k=sr;k<=er;k++)
				{
					for (int j=1;j<=ncols;j++)
					{
						if (ICE)
							ICEOLD[k][j]=(ICE[k][j]==nodata) ? nodata : ICE[k][j];
						else
							ICEOLD[k][j]=0;
					}
				}
				iceoldindex = sti;
				sti ++;

				GLACIDOLD = matrix(sr, er, 1, ncols, nrows);
				initialize_grid(sti, _GlacierCodesOld_, GLACIDOLD, NULL, 0, _GlacierCodesOld_, 0, sr,er);
				for (int k=sr;k<=er;k++)
				{
					for (int j=1;j<=ncols;j++)
					{
						if (ICE)
							GLACIDOLD[k][j]=(GLACID[k][j]==nodata) ? nodata : GLACID[k][j];
						else
							GLACIDOLD[k][j]=0;
					}
				}
				glacidoldindex = sti;
				sti ++;

				QSNOWGL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(sti, QSNOWGLgridname, QSNOWGL, NULL, 0, _SnowMelt_onGlacier_, QSNOWGLwritegrid, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j=1;j<=ncols;j++)	
					{
						QSNOWGL[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
					}
				}
				qsnowglindex = sti;
				sti++;

				EVRSGL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(sti, EVRSGLgridname, EVRSGL, NULL, 0, _SnowEvaporation_onGlac_, EVRSGLwritegrid, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j=1;j<=ncols;j++)	
					{
						EVRSGL[k][j]=(zn_routing[k][j]==nodata) ? nodata : 0;
					}
				}
				evrsglindex = sti;
				sti++;

				if (mpi_rank == 0)
				{
					CreateGlacierStatisticsZones();
					// calculate the Volume and create the elevation belts for each glacier
					for (int glc = 0; glc < glaciers.size(); glc++)
					{
						glaciers[glc]->calcVolume(); // using the V-A-relation to determine the actual volume
						//					glaciers[glc]->setOldVolume(glaciers[glc]->getOldVolume() + glaciers[glc]->getVolume()); // using the balance difference to determine the old volume
						glaciers[glc]->setOldVolume(glaciers[glc]->getVolume()); // since the old volume is always based on the actual area as long as no resizing was done, we can and must ignore the balance change since tha last resizing event
						glaciers[glc]->calcOldArea(); // determine the old area from the old volume
						glaciers[glc]->createElevationBelts(); // create elevation bands, they will be filled later with margin cells
					}
					CountMarginsGroupedByElevationBelts();
				}
			} // for dynamic glacier model only
		}

		// for Permafrost, only this grid has to be created as internal grid
		if ((ist_vorhanden(_ThawCoeffPermaFrost_))&&(ist_vorhanden(_PermafrostThawDepth_)))
		{
			// here, we set the number of snow free days to 0. This would mean that the soil thawing starts right now.
			// if for any reason it may be more suitable to start with an already thawed soil, this grid must be changed into a grid which
			// must be either read in as variable grid or as a grid specific to the snowmodel, thus changign the control file
			if (!ist_vorhanden(_SnowFreeDaysGrid_))
			{
				SNOWFREEDAYS = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				initialize_grid(sti, "", SNOWFREEDAYS, NULL, 0, _SnowFreeDaysGrid_, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SNOWFREEDAYS[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				snowfreedaysindex = sti;
				sti++;
			}
			else
			{
				snowfreedaysindex=getindex(_SnowFreeDaysGrid_);
				SNOWFREEDAYS = gridinfo[snowfreedaysindex].pointer;
			}

			if (!ist_vorhanden(_SnowFreeDaysGrid_))
			{
				SNOWCOVERDAYS = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				initialize_grid(sti, "", SNOWCOVERDAYS, NULL, 0, _SnowCoverDaysGrid_, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SNOWCOVERDAYS[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				snowcoverdaysindex = sti;
				sti++;
			}
			else
			{
				snowcoverdaysindex = getindex(_SnowCoverDaysGrid_);
				SNOWCOVERDAYS = gridinfo[snowcoverdaysindex].pointer;
			}

			ALPHA     = gridinfo[getindex(_ThawCoeffPermaFrost_)].pointer;
			THAWDEPTH = gridinfo[getindex(_PermafrostThawDepth_)].pointer;

			sprintf(name,"%s%d",_Wurzeltiefe_,1);
			ROOTDEPTH = gridinfo[getindex(name)].pointer;

			steuerdatei=fopen(steuerdateiname,"r");
			findeabschnitt("[" _permafrost_ "]");

			fscanfln(steuerdatei,"%d",&permfrostmethod,"Permafrost-Method for ",tstr);
			if ((permfrostmethod<1)||(permfrostmethod>1))
			{
				(*WSstream) << "FATAL: method for permafrost melt must be 1 (alpha*sqrt(time))"<<"[end-entry]";
				throw(80);
			}
			// hier muessen die Parameter eingelesen werden!! 

			thawdepthstatname   = newstring(_stringsize_+1,1);

			fscanfln(steuerdatei,"%lf",&MinSAgeRefreeze  ,"Min_Snow_Age_Refreeze for "    ,tstr);
			fscanfln(steuerdatei,"%lf",&MaxSWE_SnowFree ,"Max_SWE_SnowFree for "    ,tstr);
			fscanfln(steuerdatei,"%s%d",thawdepthstatname,&thawdepthstatcode  ,"name of thaw depth statistics file "    ,tstr);
			fscanfln(steuerdatei,"%lf",&reduced_ksat,"reduced k_sat (in m/s)",tstr);

			zoneTHAWDEPTH = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneTHAWDEPTH) speichererror(1);
			for (i=0;i<=zonmax;i++)
			{
				zoneTHAWDEPTH[i].n           = zoneTHAWDEPTH[i].sx        = 0;
				zoneTHAWDEPTH[i].zonenwert   = zone[i].zonenwert;
				zoneTHAWDEPTH[i].rel         = zone[i].rel;
				zoneTHAWDEPTH[i].RelCellSum  = zone[i].RelCellSum;
			}
			thawdepth = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneTHAWDEPTH,
			(int)zonmax,(int)thawdepthstatcode,(int)4,(char*)thawdepthstatname,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			thawdepth->schreibzonenwerte("thaw_depth_in_m_",instanzname,LastResults,maxLastResults);

			if ((fclose(steuerdatei))==EOF)
			{
				(*WSstream) << "FATAL: cannot close "<<steuerdateiname<<"[end-entry]";
				throw(11);
			}

		}

		if (T0R<T0)
		{
			(*WSstream) << "ERROR: \n-----------------------------------------------------------------------\n"
			<< "   T0R (threshold temperature rain/snow) < T0 (threshold temperature snow melt)  \n"
			<< "   because these thresholds are contradicting each other, T0R was set to T0!                    \n"
			<< "-----------------------------------------------------------------------"<<"[end-entry]";
			T0R=T0;
		}
		initialized=0; // darf erst nach SetQgesAdresse auf 1 gesetzt werden

		if (doSpinUp!=0 && HRU==NULL)
		{
			HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
			HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
		}
	}
	catch (int ex)
	{
		retvalue = ex;
	}
	catch (...)
	{
		retvalue = 103;
	}
}

double schneemodellklasse::RMF(long tg,long mt)
{	long t;
	int i;
	// errechnet aus der uebergebenen Tagnummer tag den aktuellen
	// Radiation-Melt-Factor, wobei a, b, c und d folgende Parameter sind:

	// r_b = (RMFMAX - RMFMIN)/2.0  Amplitude
	// r_a = RMFMIN + r_b   Additative Konstante der Sinusfunction
	// r_c = 2*pi/365 = 0.017214206 zur Anpassung der Periode der Sinusschwingung
	// r_d = -80 (Differenz zwischen dem 1.1. und dem 21.3, wenn sin(c*(d+tag)) == 0)

	// Errechnen der Tagesnummer aus dem an die function uebergebenen Datums-Angaben 
	t = tg;
	for (i=1;i<mt;i++) t += mon[(int)mt];
	return (r_a + r_b*sin(r_c*(r_d+t)));
}

double schneemodellklasse::Asnow(int subbasin,long tg,long mt)
{	long t;
	int i;
	t = tg;
	for (i=1;i<mt;i++) t += mon[(int)mt];
	return (Asnow_a[subbasin] + Asnow_b[subbasin]*sin(r_c*(r_d+t)));
}

double schneemodellklasse::Aice(int subbasin,long tg,long mt)
{	long t;
	int i;
	t = tg;
	for (i=1;i<mt;i++) t += mon[(int)mt];
	return (Aice_a[subbasin] + Aice_b[subbasin] *sin(r_c*(r_d+t)));
}

double schneemodellklasse::esatt(double t)
{	return (6.1078 * exp((21.874*t)/(265.5+t)));
	// gilt fuer E(eis), E(wasser) = 6.1078 * exp((17.269*t)/(237.3+t))
	// Angaben nach Muench, 1994; dagegen in Weischet, Einfuehrung in die Klimatologie:
	// E(wasser) = 6.1078 * exp((17.08085*t)/(234.175+t))  
}

// snow stack must be update layer by layer for added snow and also when snow melts (+/- sign of newSnow)
// ChangeLiquid is actually not used. May be later.
void schneemodellklasse::UpdateSWEStack(int i, int j, double Anteil, double ChangeSolid, double ChangeLiquid ,bool onGlacier,bool fromRefreezing,double externalTemperature)
{
	if (onGlacier) return; // on glacier we will not do the layered snow model at all for now.

	int levelsSnow = LayeredSnowProperties[1].getLevels();
	int actLevelSnow = (int)refUpperIndexSnow[i][j];

	if (ChangeSolid >= 0) // new fallen snow or refreezing water in the snow pack that needs to be added to swes and subtracted from swel
	{

		// if there is no snow yet, start a new level (which is levelsSnow, the lowermost layer
		if (actLevelSnow > levelsSnow)
		{
			actLevelSnow = levelsSnow;
		}

		double SnowInExcess = ChangeSolid * Anteil;
		while ((SnowInExcess > 0) && (actLevelSnow > 0)) // actLevelSnow is always > 0, so this last condition is not really required here, just makes the code more robust for future changes
		{
			// if the amount of snow is too much for this layer, reserve the exceeding snow for the next layer
			// actual content: SWES in mm = kg/m^2
			// new Snow in mm = kg/m^2
			// maximum content of a layer: thickness[m] * density [kg/m^3] (each from LayeredSnowProperties or of the actual stack) -> result in kg/m^2 or mm
			double defaultDensity = LayeredSnowProperties[1].getDryDensity()[actLevelSnow]/1000.0; 
			double defaultThickness = LayeredSnowProperties[1].getDzui()[actLevelSnow];

			// old
			// double maxSWE = refDensSnow[i][j][actLevelSnow] * 1000.0 * refThickSnow[i][j][actLevelSnow];
			// double MaxSnowAmountToAdd = maxSWE - refSWES[i][j][actLevelSnow];

			// new
			double fillableThickness = 0;
			double MaxSnowAmountToAdd = 0;
			if (actLevelSnow>1)
			{
				fillableThickness = (refThickSnow[i][j][actLevelSnow] <= -9999) ? defaultThickness : (defaultThickness - refThickSnow[i][j][actLevelSnow]);
				MaxSnowAmountToAdd = fillableThickness * defaultDensity*1000.0;
			}
			else
			{
				MaxSnowAmountToAdd = SnowInExcess; 
			}

			if (MaxSnowAmountToAdd > 0) // there is less free SWE left in the actual layer than we have new snow fall
			{
				// first mixing the temperatures of new snow and old snow
				// note: the temperature of rain will be accounted for in the heat transfer model itself (when convection is on)
				double actSnowAmountToAdd = (MaxSnowAmountToAdd>SnowInExcess)?SnowInExcess:MaxSnowAmountToAdd;
				if (bHeatTransferActive && (!onGlacier))
				{
					double tsnow = (externalTemperature+dTOffsetFreshSnow <= 0) ? externalTemperature+dTOffsetFreshSnow : 0.0;
					double Tnew = (refSWES[i][j][actLevelSnow] * refTemperSnow[i][j][actLevelSnow] + actSnowAmountToAdd*tsnow) / (refSWES[i][j][actLevelSnow] + actSnowAmountToAdd);
					refTemperSnow[i][j][actLevelSnow] = Tnew;
				}
				refSWES[i][j][actLevelSnow] += actSnowAmountToAdd;

				double oldThickness = (refThickSnow[i][j][actLevelSnow] <= -9999) ? 0 : (refThickSnow[i][j][actLevelSnow]);
				refThickSnow[i][j][actLevelSnow] = oldThickness+(actSnowAmountToAdd/1000.0/defaultDensity); // result in m, so actSnowAmountToAdd/1000.0 is conversion form mm to m, since density is used in Kg/dm^3 and not in kg/m^3, i.e. there is no factor 1000 alreday included
				refDensSnow[i][j][actLevelSnow] = (refThickSnow[i][j][actLevelSnow]!=0) ? (refSWES[i][j][actLevelSnow]+refSWEL[i][j][actLevelSnow])/(1000.0*refThickSnow[i][j][actLevelSnow]) : 0; // in kg/dm^3

				SnowInExcess -= actSnowAmountToAdd;

				if ((actLevelSnow>1)&&(SnowInExcess>0)) 
				{
					actLevelSnow--; // if actLevelSnow reaches zero, this loop is broken and the remainig nsow must be added to the uppermost layer (or distributed to all layers or whatever..., not sure yet)
				}
			}
			else
			{
				actLevelSnow--;
			}
		}
	}
	else // snow melt
	{
		double SnowToMelt = -ChangeSolid * Anteil; // newsnow is negative!
		while ((SnowToMelt > 0) && (actLevelSnow <= levelsSnow))
		{
			double SWESAvailableInLayer = refSWES[i][j][actLevelSnow];
			if (SWESAvailableInLayer <= SnowToMelt) // there is more snow to be molten than what is available in the actual layer, so make the actual layer empty and move to the next layer
			{   // no temperature mixing here
				SnowToMelt -= SWESAvailableInLayer;
				refSWES[i][j][actLevelSnow] = 0;
				refSWEL[i][j][0] += refSWEL[i][j][actLevelSnow]; // content of swel will be transferred to the melt and thus be available as liquid input in unsatzonmodel again, even if there was no rain fall
				refSWEL[i][j][actLevelSnow] = 0;
				// since this level is empty now, reset all other entities to nodata
				refDensSnow[i][j][actLevelSnow] = -9999; 
				refThickSnow[i][j][actLevelSnow] = -9999;
				if (bHeatTransferActive)
				{
					refTemperSnow[i][j][actLevelSnow] = -9999;
				}
				actLevelSnow++; // if actLevelSnow reaches levelsSnow, this loop is broken and the remainig melt must be returned 
			}
			else
			{   // no temperature mixing here
				refSWES[i][j][actLevelSnow] -= SnowToMelt;
				// snow density remains constant when melt occurs but thickness decreases
				// in order to keep thickness and density synchronous with swes, thickness is calculated from those variables now
				// but only for the uppermost snow layer, all layers below will not be touched by this routine except in unsatzon model (if heat transfer suggests melting in a layer)
				refThickSnow[i][j][actLevelSnow] = ((refSWES[i][j][actLevelSnow]+refSWEL[i][j][actLevelSnow])/refDensSnow[i][j][actLevelSnow])/1000.0;
				SnowToMelt = 0;
			}
		}
		if (SnowToMelt > 0)
		{
			SnowToMelt = 0;
		}
	}

	refSWES[i][j][0] = 0;
	for (int z = 1; z <= levelsSnow; z++)
	{
		refSWES[i][j][0] += refSWES[i][j][z];
	}

	refUpperIndexSnow[i][j] = actLevelSnow;

#ifdef _WASIMDEBUG_
	if (fabs(refSWES[i][j][0] - SSNOW[i][j]) > 0.001)
	{
//		(*WSstream) << "DEBUG: hour "<<this->stunde<<" difference in SSNOW "<<SSNOW[i][j]<<"and SWES "<<refSWES[i][j][0]<<" in row "<<i<<" and col "<<j<< "[end-entry]";
	}
#endif
}

// calling the snow melt routine: i and j are coordinates, temperature is the air temperature,
// Anteil is the fraction of the timestep, for which the melt will be calculeted
// (if in a daily timestep two calculations for day and night should be done) 
// if onGlacier is true, then this is always the second call and all changes of input (rain/snow etc) must not be made in the second call again
double schneemodellklasse::calc_Accu_and_Melt(int i, int j, double temperatur, double Anteil, double &schneeanteil, double &Ssnow, double &Sliq, double &Ssto, double &Snowtemp, double &Surftemp, double &Evrs, double &EnergyBalance, double solid_throughfall, bool onGlacier)
{	double faktor; // zur Expositionsabhaengigen Korrektur des Strahlungsschmelzfaktors 
	double RMFkorr; // zur Expositionsabhaengigen Korrektur des Strahlungsschmelzfaktors 
	double M,Mneg,Mh,Me,Mr,Ms;
	// variables for new snoe model extensions (Warscher, 2103)
	double Ebal,Radflux,Sensflux,Latflux,Adflux,Soilflux,Mass_lat,bolz_clouds;
	double LWin, LWout, SurfaceTemperature, SnowTemperature, AirTemperature;
	double VapPresAir, VapPresSnow, AirPres;
	double WindCorr;
	double temp1,temp2;
	int nr_iterations;
	//							 M    = Gesamtschmelze als Summe der Teilschmelzen oder im T-Index-Verfahren
	//                           Mneg = negative Schmelze,
	//                           Mh   = Schmelze durch sensible Waerme
	//                           Me   = Schmelze durch latente Waerme
	//                           Mr   = Schmelze durch Regen
	//                           Ms   = Schmelze durch Strahlung 
	// ----------- Aufteilen des Niederschlages auf die --------------- 
	// ----------- beiden Schneespeicher fest und fluessig ------------ 

	// for LAKEs, all snow is assumed to melt immediately (later, we could model the accumulation of ice, e.g. by a thermal-time-approach (degree-days)
	if ((lake_active)&&(LAKE[i][j]>0)&&(POND[i][j]>0))
	{
		schneeanteil=0.0;
	}
	else
	{
		if (temperatur <= (T0R+transzone))
		{
			// Lufttemperatur unter der Grenztemperatur T0R 
			// Ermitteln des Schneeanteils am Niederschlag 
			if (temperatur>(T0R-transzone))
			schneeanteil=((T0R+transzone)-temperatur)/(transzone+transzone);
			else
			schneeanteil = 1.0;
		}
		else
		schneeanteil = 0.0;
	}
	if (!onGlacier)
	{
		// Niederschlagsstatistik, doesn't need to be updated in second call
		P_SNOW[i][j] = schneeanteil*P[i][j];
		P_RAIN[i][j] = (1-schneeanteil)*P[i][j];

		// Correction wind redistribution -> do it only for the first call with onGlacier=false!
		if (methode == 7 || methode == 9 || methode == 10 || methode == 11 || methode == 13 || methode == 14)
		{
			P[i][j] += (DEPOINDEX[i][j]*P_SNOW[i][j]);
			P_SNOW[i][j] = P_SNOW[i][j] + (DEPOINDEX[i][j]*P_SNOW[i][j]);
		}
		// $kf 2016-07-20: add solid throughfall to snow accumulation
		P_SNOW[i][j] += solid_throughfall;
	}
		// Akkumulation: wenn es zu warm fuer Schnee ist, dann wird der Regen in SLIQ
		//gespeichert. Bei der Bilanzierung ergibt dieser Anteil den Niederschlag fuer das
		//folgende Teilmodell. 
	Ssnow += P_SNOW[i][j]*Anteil;
	Sliq  += P_RAIN[i][j]*Anteil;
	Ssto = Ssnow+ Sliq;
	Evrs = 0;
	// now updating the SWES stack, if present 
	if (bLayeredSnowActive && P_SNOW[i][j] > 0)
	{
		// during development, we assume that the rain fraction is addeed to the swel as well. This will be replaced later by the following (commented) line with 0 rain input (which is then the input to the Richardsmodel)
		// UpdateSWEStack(i, j, Anteil, P_SNOW[i][j], P_RAIN[i][j], onGlacier, _fromRefreezingFalse_);
		UpdateSWEStack(i, j, Anteil, P_SNOW[i][j], 0, onGlacier, _fromRefreezingFalse_,T[i][j]);
	}

	// -------- Rueckgefrieren und Schmelze werden nun berechnet ------ 

	// Initialize variables
	M=0.0;
	Mneg=0.0;
	if (methode>4)
	{

		Mass_lat=0.0;
		SnowTemperature = Snowtemp;
		AirTemperature=temperatur+273.16;
		if ((bHeatTransferActive && bLayeredSnowActive)&& (!onGlacier))
		{
			Soilflux = 1.0; // will later be calculated in heat transfer depending on the temperature of the lowest snow layer and the uppermost soil layer and returned in grid to the snow model
		}
		else
		{
			Soilflux = 1.0;
		}
		

		if ((Surftemp == nodata) || (Surftemp == 0.0))
		{
			Surftemp=AirTemperature-3;
		} 
		if (SURFTEMPINTERPOL)
		{
			SurfaceTemperature = SURFTEMPINTERPOL[i][j];
			if (SurfaceTemperature < 100) // input is probably not in Kelvin but in °C
			{
				SurfaceTemperature += 273.16; // convert into Kelvin
			}
		}
		else
		{
			SurfaceTemperature = Surftemp;
		}
		
		if (SurfaceTemperature>273.16)
		{
			SurfaceTemperature=273.16;
		}

		if ((SnowTemperature == nodata) || (SnowTemperature == 0.0)) 
		{
			if (AirTemperature < 273.16) 
			{
				SnowTemperature = AirTemperature;
			}  
			else 
			{
				SnowTemperature = 273.16;
			}
		}
	}

	// Rueckgefrieren und Schmelze haengt von T0 und nicht von T0R ab 
	if (temperatur <= T0)
	// Rueckgefrieren 
	{
		if ((methode == 3) || (methode == 4)) // Kombinierte T-U-Index - Energiebilanzmethoden 
		{
			// Es wird das Rueckgefrieren berechnet, Vorzeichen negativ! 
			Mneg = fabs(CRFR * RMFakt * (temperatur - T0) / InVtag);
			Mneg *= Anteil; // mit Laenge des der Routine uebergebenen Zeitintervalls wichten 

			if (Mneg <= Sliq) // ist weniger Wasser rueckgefroren, als im Schnee steckte? 
			{ // ja 
				Sliq -= Mneg;  // Schneewasserspeicher um M reduzieren 
				Ssnow += Mneg; // Scheespeicher um M erhoehen 
			}
			else
			{ // nein, wasser reichte nicht aus 
				Ssnow += Sliq; // Schneespeicher nur um den bisher fluessigen Anteil erhoehen 
				Sliq = 0;            // Fluessigspeicher 0 setzen 
			}
			if (!bHeatTransferActive && !onGlacier)
			{ // if heat transfer is active, Mneg is calculated in soil model and the balance is handled there, too
				if (bLayeredSnowActive && (Mneg != 0))
				{
					// will be called when heat transfer is off only, so the negative melt must be handled here
					UpdateSWEStack(i, j, Anteil, Mneg, 0, onGlacier, _fromRefreezingTrue_,T[i][j]);
				}
			}
		}
		//new energy balance calculation (Warscher, 2013)
		else if (methode == 5 || methode == 6 || methode == 7 || methode == 11)
		{   
			if (Ssnow > 0.0)
			{
				// ToDo: (2015-02-14) sub time step for more realistic energy balance calculation (since surface temperature can change during long time steps significantly and thus the energy balance changes as well

				if (AirTemperature > 273.16) 
				{
					temp1=(17.08085*(AirTemperature-273.16))/(234.175+(AirTemperature-273.16));
					VapPresAir = 6.1078 * exp(temp1)*RLF[i][j];
				}
				else 
				{
					temp1=(22.4429*(AirTemperature-273.16))/(172.44+(AirTemperature-273.16));
					VapPresAir = 6.1071 * exp(temp1)*RLF[i][j];
				}	
				if (P[i][j] == 0.0) {
					bolz_clouds = 0.17;
					if (CL[i][j] <= 0.2) {bolz_clouds = 0.08;}
					if (CL[i][j] >= 0.6) {bolz_clouds = 0.20;}
				}
				else {bolz_clouds = 0.24;}	     
				temp1=AirTemperature*AirTemperature*AirTemperature*AirTemperature;
				temp1=temp1*(0.610+0.050*sqrt(VapPresAir));
				temp2=(1.0+bolz_clouds*CL[i][j]*CL[i][j]);
				LWin = 0.0000000567*LWINcorr*temp1*temp2;
				// $kf 2016-08-22: use canopy conditions for incoming longwave radiation if available
				if (SCN_ir_radiation)
				{
					double LWin_scn = this->SCN_ir_radiation[i][j] * LWINcorr;
					if (LWin_scn > 0.)
					{
						// use these values only in case valid values are available
						LWin = LWin_scn;
					}
				}
				//mk: Adflux in CEB??
				Adflux = ((4180.0 * (AirTemperature - 273.16)) * P_RAIN[i][j] / (InV*60.0)) + ((2100.0 * (AirTemperature - 273.16)) * P_SNOW[i][j] / (InV*60.0));

				if ((Sliq > 0) || (AirTemperature > 275.16)) // if the snow pack contains liquid water, we assume the surface temperature to be at 0°C. Also, if Air temperature is above 2°C
				{
					SurfaceTemperature = 273.16; // first guess for iterative calculation of surface temperature or limiting when surface temperature is known 
				}
				else // otherwise, first guess surface temperature is 2K below air temperature (which is usually below zero because we are in the refreezing branch)
				{
					if (!SURFTEMPINTERPOL)
					{
						SurfaceTemperature = AirTemperature - 2;
					}
				}
					
				// iteration to find a surface temperature that brings the energy balance into equilibrium
				double EBALold=0;
				Ebal=0;
				double deltaT=-10;
				double reducefact=0.5;

				nr_iterations=0;
				bool signchanged=false;
				do
				{
					EBALold=Ebal;
					nr_iterations++;
					temp1=SurfaceTemperature*SurfaceTemperature;
					temp1=temp1*temp1;
					LWout= -0.0000000567*LWOUTcorr*temp1;	               	
					temp1=(22.4429*(SurfaceTemperature-273.16))/(172.44+(SurfaceTemperature-273.16));
					VapPresSnow = 6.1071 * exp(temp1);
					Radflux=((1.0-Albedo[i][j])*RGEX[i][j]*60.0/InV) + LWin + LWout ;
					Sensflux=18.85*(0.18+0.098*U[i][j])*(AirTemperature-SurfaceTemperature);
					Latflux=32.82*(0.18+0.098*U[i][j])*(VapPresAir-VapPresSnow);
						
					Ebal=Radflux+Sensflux+Latflux+Adflux+Soilflux;

					if ((fabs(Ebal) > 0.25) && (nr_iterations <= 150) && ((Sliq <= 0.0) || (bHeatTransferActive && bLayeredSnowActive)) && !SURFTEMPINTERPOL)
					{
						if (Ebal < 0)
						{ // meaning: temperaturwas too high, we are loosing too much energy, so we need to set SurfaceTemperature to lower values
							// now we will compare against the previous value
							if (EBALold>0) // is only possible after the second iteration, since initial EBAL
							{ // old value was with other sign, so we need to go back a bit, i.e. we set deltaT to -0.5*deltaT
								deltaT = (fabs(Ebal/(Ebal-EBALold))<1?-fabs(Ebal/(Ebal-EBALold)):-reducefact)*deltaT; // here we assume that the step to be taken backward is in a first approximation linearly connected to the relation of new and old Ebal
								signchanged=true;
							}
							else
							{ // EBALold was already < 0 and we are still < 0, so we need to push the temperature further in the same direction
								deltaT = (signchanged?0.5:1.0)*deltaT;
							}
						}
						else // Ebal>0
						{ // meaning: surface temperature too low, we need to lose more energy, so we need to rise the surface temperature 
							if (EBALold>0)
							{
								deltaT = (signchanged?0.5:1.0)*deltaT;
							}
							else
							{
								if (nr_iterations>1)
								{
									deltaT = (fabs(Ebal/(Ebal-EBALold))<1?-fabs(Ebal/(Ebal-EBALold)):-reducefact)*deltaT; // here we assume that the step to be taken backward is in a first approximation linearly connected to the relation of new and old Ebal
									signchanged=true;
								}
								else
								{
									deltaT = (signchanged?-0.5:-1.0)*deltaT;
								}
							}
						}
						SurfaceTemperature += deltaT;
					}
				}
				while ((fabs(Ebal) > 0.25) && (nr_iterations <= 150) && ((Sliq <= 0.0) || (bHeatTransferActive && bLayeredSnowActive)) && !SURFTEMPINTERPOL); // loop is done only once when SURFTEMPINTERPOL was found 
				
				EnergyBalance = 0; // for external use, assume Ebal to be 0 now (e.g. in heat transfer model)

				if (!SURFTEMPINTERPOL)
				{
					if ((AirTemperature - SurfaceTemperature) > 13.0)
					{
						SurfaceTemperature = AirTemperature - 13.0;
						temp1 = (22.4429*(SurfaceTemperature - 273.16)) / (172.44 + (SurfaceTemperature - 273.16));
						VapPresSnow = 6.1071 * exp(temp1);
						Latflux = 32.82*(0.18 + 0.098*U[i][j])*(VapPresAir - VapPresSnow);
					}

					if (SurfaceTemperature > 273.16)
					{
						SurfaceTemperature = 273.16;
						Sensflux = 18.85*(0.18 + 0.098*U[i][j])*(AirTemperature - SurfaceTemperature);
						temp1 = (22.4429*(SurfaceTemperature - 273.16)) / (172.44 + (SurfaceTemperature - 273.16));
						VapPresSnow = 6.1071 * exp(temp1);
						Latflux = 32.82*(0.18 + 0.098*U[i][j])*(VapPresAir - VapPresSnow);
					}
				}

				// computing temperature change of the snow pack
				// we are in the path T<T0 and there is no liquid fraction present
				// The change in snow temperature can be estimated from the ols snow temperature, the surface temperature and the thermal properties of snow
				// we take the surface temperature as the upper boundary condition for the heat transport equation for a one-layer snow pack (lower boundary condition: soil heat flux)
				// the heat flux is defined by 
				// q_snow = lambda_snow * dT/dx, i.e. heat conduction of snow * temperature gradient
				// this energy flux in J/(s*m^2) is taken to heat up/cool down the complete snow pack:
				// rho*V*C_snow*dT/dt = q_snow
				// with SWE=1000*V und x=0.003*SWE (measure of WE in mm, x in m) we get
				// Tsnow,new = Tsnow,old + lambda_snow * (Tsurf-Tsnow,old)*dt/(SWE*SWE*C_snow)*1000*2*0.33
				// with 0.33 as density of snow (which varies between 0.05 und 0.6, but for old snow it is something around 0.33)
				// lamda is assumed to be 0.3 (likewise for older snow)

				double lambda = 0.3;
				double SnowDensity = 0.3;
				double Csnow= 2100; // heat capacity of snow, rough estimation only; [Csnow]=J/(Kg*K)
				temp1=(Anteil*InV*60.0)*1000*2*SnowDensity/(Ssto*Ssto*Csnow);
				double oldSnowTemp = SnowTemperature;
				SnowTemperature = SnowTemperature + lambda * (SurfaceTemperature - oldSnowTemp)*temp1;

				// limiting thermal flow: surface temperature cannot be over/-underrun in either direction, depending on the original gradient.
				if (oldSnowTemp >= SurfaceTemperature)
				{ // flux from snow to surface, i.e. snow is cooling down
					// snow temperature may not be lower than the surface temperature
					if (SnowTemperature < SurfaceTemperature)
					{
						SnowTemperature = SurfaceTemperature;
					}
				}
				else // snow was cooler than surface
				{ // flux from surface into snow, i.e. snow is warming up
					// snow temperature may not rise abve surface temperature
					if (SnowTemperature > SurfaceTemperature)
					{
						SnowTemperature = SurfaceTemperature;
					}
				}

//				if (SnowTemperature < SurfaceTemperature-3.0) 
//				{
//					SnowTemperature=SurfaceTemperature-3.0;
//				}
//				else if (SnowTemperature > 273.16)
//				{
//					SnowTemperature = 273.16;
//				}

				if (SnowTemperature > 273.16)
				{
					SnowTemperature = 273.16;
				}

				Mass_lat=(Latflux*(InV*60.0))/2835500.0;
				Mass_lat *= Anteil; // only for the valid fraction of the time step

				Ssnow = Ssnow + Mass_lat;
				if (Ssnow<0) // when Mass_lat is negative
				{
// for compatibility with layered snow model we remove snow evaporation from solid snow only now (20170326)
//					Sliq += Ssnow; // trying to take evaporation from liquid fraction first before reducing evaporation
//					Ssnow = 0;
//					if (Sliq<0)
//					{
//						Mass_lat -= Sliq; // reducing evaporation loss from snow to the maximum amount possible
//						Sliq = 0;
//					}
					Mass_lat -= Ssnow;
					Ssnow = 0;
					Ssto = Ssnow + Sliq;
				}

				// now updating the SWES stack, if present: add to or remove from the SWES-stack the amount of latent heat flux (evaporation or sublimation) 
				if (bLayeredSnowActive && (Mass_lat!=0))
				{
					UpdateSWEStack(i, j, Anteil, Mass_lat, 0, onGlacier, _fromRefreezingFalse_,T[i][j]);
				}

				if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier)
				{ // if heat transfer is active, Ebal is always 0 and refreezing is computed in the heat transfer model. This old way of estimating the refreezing is 
				  // required only because there is no simple way to estimate the refreezing amount. So the usage of the energy balance is only a work around 
				  // (and by the way: setting the surface temperature to 0°C for conditions when there is any liquid water in the snow pack is a bit awkward in any case)
					Mneg = -(Ebal*(Anteil*InV*60.0)) / 337500.0;
					if ((Mneg > 0.0) && (Sliq>0.0))
					{
						if (Mneg <= Sliq) /* ist weniger Wasser rueckgefroren, als im Schnee steckte? */
						{ /* ja */
							Sliq = Sliq - Mneg;  /* Schneewasserspeicher um M reduzieren */
							Ssnow = Ssnow + Mneg; /* Scheespeicher um M erhoehen */
						}
						else
						{ /* nein, wasser reichte nicht aus */
							Mneg = Sliq;
							Ssnow += Mneg; /* Schneespeicher nur um den bisher fluessigen Anteil erhoehen */
							// until 2017 03 10		SnowTemperature = SnowTemperature - ((((Mneg - Sliq)*337500.0)) / (Ssnow*2100.0));
							// when water freezes, energy is set freee, so the temperature should rise a bit instead of fall
							SnowTemperature = SnowTemperature + (((Mneg*337500.0)) / (Ssnow*2100.0));
							Sliq = 0;            /* Fluessigspeicher 0 setzen */
							if (SnowTemperature < SurfaceTemperature - 3.0)
							{
								SnowTemperature = SurfaceTemperature - 3.0;
							}
							else if (SnowTemperature > 273.16)
							{
								SnowTemperature = 273.16;
							}
						}

						if (bLayeredSnowActive && (Mneg != 0))
						{
							// will be called when heat transfer is off only, so the negative melt must be handled here
							UpdateSWEStack(i, j, Anteil, Mneg, 0, onGlacier, _fromRefreezingTrue_,T[i][j]);
						}

					}
				}
				// Ssnow and Sliq may have been changed by sublimation and (if heat transfer is not active) by refreezing, so Ssto has to be updated as well
				Ssto = Ssnow+ Sliq;

				// sublimation from air to snow: negative amount must be subtracted from EVRS as well to get the balance correctly
				// the amount of Mass_lat which was added to SSNOW must be subtracted from EVRS, the snow evaporation/sublimation
				Evrs = -Mass_lat;

				Surftemp = SurfaceTemperature;
				if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier)
				{ // change snow temperature when heat transfer is not done in unsatzonmodell
					Snowtemp = SnowTemperature;
				}
				if (Ssto <= 0.0) {Snowtemp=nodata;}
				if (Ssto <= 0.0) {Surftemp=nodata;}
			}	
			else
			{
				Surftemp = nodata;
				Snowtemp = nodata;
			}
		}
	}
	else // Schmelze berechnen
	{
		// Lufttemperatur T ueber der Grenztemperatur 
		if (Ssnow>0) // ist noch schmelzbarer Schnee da?
		{
			if (methode==1 || methode == 8 || methode==9 || methode==10)
			{
				M = actC0*(temperatur-T0)/InVtag;
			}
			else
			{
				if ((methode==2) || (methode==12) || (methode==13) || (methode==14))	{
					M = (actC1 + actC2*U[i][j])*(temperatur-T0)/InVtag;
				} else if (methode==3)	{
					if (P[i][j] > (2.4/InVtag)) // konvektives Ereignis mit Niederschlagswaerme 
					{
						//	Ms = 1.2*temperatur/InVtag;              // Radiation Melt
						Ms = 0.6*temperatur/InVtag;              // Radiation Melt
						Mh = (actC1 + actC2*U[i][j])*temperatur/InVtag;   // infolge fuehlbarer Waerme 
						// nun infolge latenter Waerme; 0.662: Psychrometerkonstante in hPa/K 
						Me = (actC1 + actC2*U[i][j])*(esatt(temperatur) - 6.11)/(0.662*InVtag);
						//	Mr = 0.0125*P[i][j]*temperatur;
						Mr = 0.006*P[i][j]*temperatur;
						M  = Ms + Mh + Me + Mr;
					}
					else
					{
						// direkte Sonnenstrahlung wird korrigiert:
						//0, wenn es Nacht ist, oder die Flaeche nicht getroffen wird,
						//echte Korrektur, wenn Tag ist
						if (EXKORR)
						faktor  = 0.5*(1.0+EXKORR[i][j]);
						else faktor = 1.0;
						RMFkorr = RMFakt * faktor;
						M = RMFkorr*(temperatur-T0)/InVtag;
					}
				} else if (methode==4) {
					if (P[i][j] > (2.4/InVtag)) // konvektives Ereignis mit Niederschlagswaerme 
					{
						//	Ms = 1.2*temperatur/InVtag;              // Radiation Melt 
						Ms = 0.6*temperatur/InVtag;              // Radiation Melt 
						Mh = (actC1 + actC2*U[i][j])*temperatur/InVtag;   // infolge fuehlbarer Waerme 
						// nun infole latenter Waerme; 0.662: Psychrometerkonstante in hPa/K 
						Me = (actC1 + actC2*U[i][j])*(E[i][j] - 6.11)/(0.662*InVtag);
						//	Mr = 0.0125*P[i][j]*temperatur;
						Mr = 0.006*P[i][j]*temperatur;
						M  = Ms + Mh + Me + Mr;
					}
					else
					{
						// direkte Sonnenstrahlung wird korrigiert:
						// 0, wenn es Nacht ist, oder die Flaeche nicht getroffen wird,
						// echte Korrektur, wenn Tag ist
						if (EXKORR)
						faktor = 0.5*(1.0+EXKORR[i][j]);
						else
						faktor=1.0;
						RMFkorr = RMFakt * faktor;
						M = RMFkorr*(temperatur-T0)/InVtag;
					}
				} else if (methode == 5 || methode == 6 || methode == 7 || methode == 11) { // new snow model extensions (Warscher, 2013)
					// no sub time steps are required here since surface temperature is kept constant at 0°C and also (positive) energy balances will lead to snow melt and no change in any temperature (wether snow surface nor snow pack temperature)
					if (!SURFTEMPINTERPOL)
					{ 
						SurfaceTemperature = 273.16; // melting conditions, so the surface temperature is constant at 0°C or 273.16 K
					} // otherwise the observed surface temperature will be used
					double VapPresAirSat;
					if (AirTemperature > 273.16) 
					{
						VapPresAirSat = 6.1078 * exp((17.08085*temperatur)/(234.175+temperatur));
					}
					else 
					{
						VapPresAirSat = 6.1071 * exp((22.4429*temperatur)/(172.44+temperatur));
					}
					VapPresAir = VapPresAirSat*RLF[i][j];
					// VapPresSnow = 6.1071 * exp((22.4429*(SurfaceTemperature-273.16))/(172.44+(SurfaceTemperature-273.16)))*1.0;
					VapPresSnow = 6.1071; // since surface temperature is defined as 273.16, vapor pressure can be set to a constant value as well
					// since SurfaceTemperature is constant for melting conditions at 0°C, this expression can be pre-calculated instead of using the power function or doing four multiplications each time
					LWout= -0.0000000567*LWOUTcorr*SurfaceTemperature*SurfaceTemperature*SurfaceTemperature*SurfaceTemperature;
					//LWout = -312.2564; // if surface temperature is const. at 273.16K
					if (P[i][j] == 0.0) {
						bolz_clouds = 0.17;
						if (CL[i][j] <= 0.2) {bolz_clouds = 0.08;}
						if (CL[i][j] >= 0.6) {bolz_clouds = 0.20;}
					}
					else {bolz_clouds = 0.24;}	               	
					LWin = 0.0000000567*LWINcorr*AirTemperature*AirTemperature*AirTemperature*AirTemperature*(0.610+0.050*sqrt(VapPresAir))*(1.0+bolz_clouds*CL[i][j]*CL[i][j]);
					// $kf 2016-08-22: use canopy conditions for incoming longwave radiation if available
					if (SCN_ir_radiation)
					{
						double LWin_scn = this->SCN_ir_radiation[i][j] * LWINcorr;
						if (LWin_scn > 0.) // use these values only in case valid values are available
							LWin = LWin_scn;
					}
					Radflux = ((1.0 - Albedo[i][j])*RGEX[i][j] * 60.0 / InV) + LWin + LWout;
					Sensflux=18.85*(0.18+0.098*U[i][j])*(AirTemperature-SurfaceTemperature);
					Latflux=32.82*(0.18+0.098*U[i][j])*(VapPresAir-VapPresSnow);
					Adflux=((4180.0 * (AirTemperature-273.16)) * P_RAIN[i][j]/(InV*60.0))+((2100.0 * (AirTemperature-273.16)) * P_SNOW[i][j]/(InV*60.0));
					Mass_lat=Anteil*(Latflux*(InV*60.0))/2835500.0; 
					// Ssnow = Ssnow + Mass_lat;
					// Ssto = Ssnow+ Sliq;	
					Ebal=Radflux+Sensflux+Latflux+Adflux+Soilflux;
					// new 2016-10-17 for scaling energy balance
					Ebal *= EnBalScaling;

					if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier) // otherwise the snow temperatures will be calculated in the unsaturated zones layered snow model
					{
						if (Ebal > 0) // for negative energy balances the snow temperatur would still fall and nothing would melt. So this is handled as an exception and nothing is done here
						{
							if ((Sliq > 0.0)||(Ssnow > 0.0))
							{
								// SnowTemperature=SnowTemperature+(((Ssnow/Ssto)*Ebal*(Anteil*InV*60.0))/(Ssnow*2100.0))+(((Sliq/Ssto)*Ebal*(Anteil*InV*60.0))/(Sliq*4180.0));
								// this last line is probably wrong. The single components are added without real weighing
								// this approach here is probably better and fits for both cases: sliq=0 and sliq>0 as well as Ssnow==0 or Ssnow>0
								SnowTemperature = SnowTemperature + ((Ebal*(Anteil*InV*60.0)) / (Ssnow*2100.0 + (Sliq*4180.0)));
							}
						}
					}
					

					// calculating wet bulb temperature in order to differentiate between snow melt and sublimation (vaporization), 
					// which would require much more energy than melt and takes place when wet bulb temperature is below 0°C
					// rough estimation for actual air pressure, using the international elevation equation
					double airpressure = 1013.25 * pow((1.0 - (0.0065*DHM[i][j])/288.15),5.255);
					double mixingratio = 0.622*VapPresAir/(airpressure-VapPresAir);
					double temp = log((mixingratio*airpressure)/((0.622+mixingratio)*6.112));
					double dewpoint = 272.62*temp/(22.46-temp);
					double gamma = 0.00066; // psychrometric constant in 1/K
					// iterative estimation of Tf
					// first: initializing Tf
					double Tf=dewpoint;
					Tf = (dewpoint+temperatur+temperatur)/3.0;
					int cnt=0;
					double ESat_Tf=0;
					double VapPresAir2=0;
					double dTf= (temperatur+Tf)/2.0; // average between T and Tf is taken as first guess for delta Tf in the following iteration
					do
					{
						cnt++;
						ESat_Tf = 6.1115*exp(22.46*Tf/(272.62+Tf));
						VapPresAir2 = ESat_Tf - gamma * (temperatur - Tf) * airpressure * (1.0 + 0.00115*Tf);
						if (VapPresAir2 > VapPresAir)
						{ // too high temperature
							dTf = -dTf;
						}
						else
						{
							dTf = fabs(dTf)/2.0;
						}
						Tf+=dTf;
					}
					while ((cnt<15) && (fabs(dTf)>0.025));

					if (Tf>0.0)
					{
						if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier)
						{ // can only be calculated after updating the snow temperatrue in heat transfer model
							M = (SnowTemperature - 273.16)*2100.0*Ssnow / 337500.0;
						}
					}
					else
					{
//						tempml = (SnowTemperature - 273.16)*2100.0*Ssnow / 2835500.0; 
						// new: sublilmation calculated directly from energy balance, so we do not need the new snow temperature yet
						// slight differences to the old version are possible since the old version used the changed temperature but the new version does not change the temperature of the remaining snow
						double tempml = (Ebal*(Anteil*InV*60.0)) / 2835500.0;
						Mass_lat -= tempml;
						Mass_lat *= Anteil;
						M=0; // no melt, i.e. all energy was needed for sublimation
						Ebal = 0; // when used again in heat transfer in unsatzonmodel as upper boundary condition
					}

					if (M < 0) M = 0; // when Ebal is negative (because too much outgoing radiation even for temperatures above 0°C), the negative melt is set to 0

					if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier)
					{ // if heta transfer is active, we can much better estimate the snow pack temperature
						SnowTemperature = 273.16;
					}

					// subtracting Mass_lat (snow evaporation) from snow storage
					if (Mass_lat < 0)
					{
						Ssnow = Ssnow + Mass_lat;

						if (Ssnow<0) // when Mass_lat is negative
						{
							// for compatibility with layered snow model we remove snow evaporation from solid snow only now (20170326)
							//					Sliq += Ssnow; // trying to take evaporation from liquid fraction first before reducing evaporation
							//					Ssnow = 0;
							//					if (Sliq<0)
							//					{
							//						Mass_lat -= Sliq; // reducing evaporation loss from snow to the maximum amount possible
							//						Sliq = 0;
							//					}
							Mass_lat -= Ssnow;
							Ssnow = 0;
						}

						Ssto = Ssnow + Sliq;
						// sublimation from air to snow: negative amount must be subtracted from EVRS as well to get the balance correctly
						// the amount of Mass_lat which was added to SSNOW must be subtracted from EVRS, the snow evaporation/sublimation
						Evrs = -Mass_lat;

						// now updating the SWES stack, if present: add to or remove from the SWES-stack the amount of latent heat flux (evaporation or sublimation) 
						if (bLayeredSnowActive)
						{
							UpdateSWEStack(i, j, Anteil, Mass_lat, 0, onGlacier, _fromRefreezingFalse_,T[i][j]);
						}
					}


					Surftemp = SurfaceTemperature;
					if (!(bHeatTransferActive && bLayeredSnowActive) || onGlacier)
					{
						Snowtemp = SnowTemperature;
					}
					if (Ssto <= 0.0) {Snowtemp=nodata;}
					if (Ssto <= 0.0) {Surftemp=nodata;}
					EnergyBalance = Ebal;

				}
			}
		}
		else
		{
			if (SURFTEMP) Surftemp = nodata;
			if (SNOWTEMP) Snowtemp = nodata;
		}
	}
	return M * Anteil; // berechnete Schmelze gilt nur fuer den als Parameter
					   // uebergeben Anteil am Intervall (tags oder nachts)
}

void schneemodellklasse::setInterceptionBeforeSnow(bool b_interceptionBeforeSnow)
{
	if (do_it == 0) return;
	interceptionBeforeSnow = b_interceptionBeforeSnow;
	if (!P)
	{
		if (interceptionBeforeSnow)
		{
			char name[500];
			sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
			P = gridinfo[getindex_first_of_3(name, _Niederschlag_, _ungueltiger_Name_)].pointer;
		}
		else
		{
			P = gridinfo[getindex(_Niederschlag_)].pointer;
		}

		// $kf 2016-08-22: apply modified meteorological values for the sub-canopy
		// additional grids that are not known when snow model was initialized primarily.
		if (ist_vorhanden(_Snow_Canopy_Unload_solid_))
		{
			this->SCN_snow_throuhfall = gridinfo[getindex(_Snow_Canopy_Unload_solid_)].pointer;
		}
		else
		{
			this->SCN_snow_throuhfall = NULL;
		}

		if (ist_vorhanden(_Snow_Canopy_IR_))
		{
			this->SCN_ir_radiation = gridinfo[getindex(_Snow_Canopy_IR_)].pointer;
		}
		else
		{
			this->SCN_ir_radiation = NULL;
		}

		if ((methode != 1) && (methode != 8) && (methode != 9) && (methode != 10))
		{
			U = gridinfo[getindex_first_of_3(_Snow_Canopy_Wind_, _Windgeschwindigkeit_, _ungueltiger_Name_)].pointer; // needed for all but T-Index methods
		}

		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			RLF = gridinfo[getindex_first_of_3(_Snow_Canopy_RH_, _Luftfeuchte_, _ungueltiger_Name_)].pointer;
			RGEX = gridinfo[getindex_first_of_3(_Snow_Canopy_RAD_, _AspectCorrectedGlobalRadiation_, _ungueltiger_Name_)].pointer;
		}

		T = gridinfo[getindex_first_of_4(_Snow_Canopy_Temperature_, _Korrigierte_Temperatur_, _Temperatur_, _Temperatur_14_)].pointer;
	}
}


// first part of the run method. It runs until snow accumulation and melt is calculated for all methods, except when bHeatTransferActive is set to true.
// In that case, for energy balance methods the melt is not returned but rather the energy balance and surface tempertatures are returned. Unsatzon model can then use the heat transfer 
// algorithms to use that values as boudary conditions as well as the liquid precipitation as input to the richards approach and heat trasfer equation.
// for other methods, melt is output as usual. The second part of the run method will gather the results as usual and do all the other stuff like glacier resizing, statistics update, spin-up etc.
void schneemodellklasse::run1(long me,long se, long tg, long mt, long jr)
{

	// gibt lediglich eine Nachricht aus
	teilmodellklasse::run(me,se, tg, mt, jr,NULL,NULL);

	if (do_it == 0) return;

	//	if (ETRSchecked==false)
	//	{
	//		ETRSchecked=true;
	//		if (ist_vorhanden(_Snow_Evaporation_))
	//		{
	//			etrsindex=getindex(_Snow_Evaporation_);
	//			ETRS=gridinfo[etrsindex].pointer;
	//		}
	//	}

	// interpolate the actual C0, C1 and C2 parameters
	if (tag < 16)
	{
		actC0 = interpolate(-15, C0[(int)monat - 1], 15, C0[(int)monat], tag);
		actC1 = interpolate(-15, C1[(int)monat - 1], 15, C1[(int)monat], tag);
		actC2 = interpolate(-15, C2[(int)monat - 1], 15, C2[(int)monat], tag);
	}
	else
	{
		actC0 = interpolate(16, C0[(int)monat], 16 + mon[(int)monat], C0[(int)monat + 1], tag);
		actC1 = interpolate(16, C1[(int)monat], 16 + mon[(int)monat], C1[(int)monat + 1], tag);
		actC2 = interpolate(16, C2[(int)monat], 16 + mon[(int)monat], C2[(int)monat + 1], tag);
	}


	RMFakt = RMF(tg, mt);

	if (glacmethode > 0)
	{
		for (int i = 1; i <= zonmax_routing; i++)
		{
			Asnow_akt[i] = Asnow(i, tg, mt);
			Aice_akt[i] = Aice(i, tg, mt);
		}
	}

	// Berechnen von sunrise, sunset, dark_a, dark_e, Anteil_hell, Anteil_dunkel -> in teilmod.cc
	get_Tageslicht_Daten();

	snowstor->init(); // oeffnen des Statistikfiles fuer Schneespeicher -> modef.cc
	qschnee->init();
	qsnow_out->init();
	snowdays->init();
	snow_age->init();
	p_snow->init();
	p_rain->init();
	albedo->init();

	// new snow model extensions (Warscher 2013)
	if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
	{
		snowtemp->init();
		surftemp->init();
	}

	if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
	{
		input_mass->init();
		mobile_mass->init();
		deposition->init();
	}

	if (ICE)
	{
		snowstor2->init();
		qfirn->init();
		qice->init();
		qglac->init();
		qsnowonglac->init();
		if (glacmethode>10)
		{
			firnstat->init();
			gmbstat->init();
			gmb2stat->init();
		}
	}

	// if thawing depth for the active layer should be modelled, too)
	if (SNOWFREEDAYS) // only if SNOWFREEDAYS exists
	{
		thawdepth->init();
	}

	for (int z = 1; z <= tracer_count; z++) tracer_module[z].stat->init();

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

#pragma omp parallel for schedule(runtime)\
		shared(tg,mt,jr)  
	for (int i = mpistartrow; i <= mpiendrow; i++) // rows of the grid
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) // columns of the grid
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && !(doSpinUp && HRUMaster[i][j] == nodata))
				//      if (zn[i][j]!=nodata)
			{
				double icevalue, trms, trmsmax, schneeanteil, M, Mgl, Sliq_max_new/*,Sliq_max_old*/;
				int startlayer;
				startlayer = 1; // index in LUSE-array -> take the uppermost layer for albed calculation
				icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

				int levelsSnow=0;
				if (bLayeredSnowActive)
				{
					levelsSnow=LayeredSnowProperties[1].getLevels();
				}

				int subbasin = (int)zarray_routing[(int)zn_routing[i][j]];

				DumptypeSnow dump;
				if (this->dumpmode==true)
				{
					if (bLayeredSnowActive)
					{
						dump.refSWES		= new double[levelsSnow+1];
						dump.refSWEL		= new double[levelsSnow+1];
						dump.refThickSnow   = new double[levelsSnow+1];
						dump.refDensSnow	= new double[levelsSnow+1];
					}
					if ((ICE)&&(glacmethode>10))
					{
						dump.FIRN			= new double[firn_levels+1];
					}
					dump.LUSE           = new double[iNumVegetationLayers+1];

#ifdef _READDUMPFILESNOW
				bool readdump=true;
#else
				bool readdump=false;
#endif
					if (readdump)
					{
						std::ifstream dumpstream;
						dumpstream.open("dumpsnow.dmp",std::ios::binary);

						if (bLayeredSnowActive)
						{
							dumpstream.read((char *)(&refSWES[i][j][0]),		sizeof(double)*(1+levelsSnow));
							dumpstream.read((char *)(&refSWEL[i][j][0]),		sizeof(double)*(1+levelsSnow));
							dumpstream.read((char *)(&refThickSnow[i][j][0]),	sizeof(double)*(1+levelsSnow));
							dumpstream.read((char *)(&refDensSnow[i][j][0]),	sizeof(double)*(1+levelsSnow));
							dumpstream.read((char *)(&refUpperIndexSnow[i][j]),	sizeof(double));
						}
						dumpstream.read((char *)(&EVRS[i][j]),			sizeof(double));
						dumpstream.read((char *)(&SSNOW[i][j]),			sizeof(double));
						dumpstream.read((char *)(&SLIQ[i][j]),			sizeof(double));
						dumpstream.read((char *)(&SSTO[i][j]),			sizeof(double));
						dumpstream.read((char *)(&T[i][j]),				sizeof(double));
						dumpstream.read((char *)(&Albedo[i][j]),		sizeof(double));
						for (int y=1;y<=iNumVegetationLayers;y++) 
						{
							dumpstream.read((char *)(&(LUSE[y][i][j])),sizeof(double));
						}
						dumpstream.read((char *)(&DHM[i][j]),			sizeof(double));
						dumpstream.read((char *)(&P[i][j]),				sizeof(double));
						if (ICE)
						{
							dumpstream.read((char *)(&ICE[i][j]),			sizeof(double));
						}
						if (SSNOWGL)
						{
							dumpstream.read((char *)(&SSNOWGL[i][j]),		sizeof(double));
							dumpstream.read((char *)(&SLIQGL[i][j]),		sizeof(double));
							dumpstream.read((char *)(&EVRSGL[i][j]),		sizeof(double));
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								dumpstream.read((char *)(&SNOWTEMPGL[i][j]),	sizeof(double));
								dumpstream.read((char *)(&SURFTEMPGL[i][j]),	sizeof(double));
								dumpstream.read((char *)(&EBALGL[i][j]),		sizeof(double));
							}
						}
						if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
						{
							dumpstream.read((char *)(&SNOWTEMP[i][j]),		sizeof(double));
							dumpstream.read((char *)(&SURFTEMP[i][j]),		sizeof(double));
							dumpstream.read((char *)(&EBAL[i][j]),			sizeof(double));
							dumpstream.read((char *)(&RLF[i][j]),			sizeof(double));
							dumpstream.read((char *)(&CL[i][j]),			sizeof(double));
							dumpstream.read((char *)(&RGEX[i][j]),			sizeof(double));
						}
						if ((ICE)&&(glacmethode>10))
						{
							dumpstream.read((char *)(&FIRN[i][j][0]),		sizeof(double)*(firn_levels+1));
						}
						if (lake_active)
						{
							dumpstream.read((char *)(&LAKE[i][j]),			sizeof(double));
							dumpstream.read((char *)(&POND[i][j]),			sizeof(double));
						}
						if ((methode != 1) && (methode != 8) && (methode != 9) && (methode != 10))
						{
							dumpstream.read((char *)(&U[i][j]),				sizeof(double));
						}
						if (EXKORR)
						{
							dumpstream.read((char *)(&EXKORR[i][j]),		sizeof(double));
						}
						dumpstream.read((char *)(&me), sizeof(int));
						dumpstream.read((char *)(&se), sizeof(int));
						dumpstream.read((char *)(&tg),sizeof(int));
						dumpstream.read((char *)(&mt),sizeof(int));
						dumpstream.read((char *)(&jr),sizeof(int));

						dumpstream.close();
					}
					else // read values should be moved to te dump-array in order to allow a proper comparison after core routines that simulates the decision if a dump file will be written
					{
						if (bLayeredSnowActive)
						{
							memmove(&(dump.refSWES[0]		),&(refSWES[i][j][0]		),sizeof(double)*(levelsSnow+1));
							memmove(&(dump.refSWEL[0]		),&(refSWEL[i][j][0]		),sizeof(double)*(levelsSnow+1));
							memmove(&(dump.refThickSnow[0]	),&(refThickSnow[i][j][0]	),sizeof(double)*(levelsSnow+1));
							memmove(&(dump.refDensSnow[0]	),&(refDensSnow[i][j][0]	),sizeof(double)*(levelsSnow+1));
							dump.refUpperIndexSnow=refUpperIndexSnow[i][j];
						}
						dump.EVRS=EVRS[i][j];
						dump.SSNOW=SSNOW[i][j];
						dump.SLIQ=SLIQ[i][j];
						dump.SSTO=SSTO[i][j];
						dump.T=T[i][j];
						dump.Albedo=Albedo[i][j];
						for (int y=1;y<=iNumVegetationLayers;y++)
						{
							dump.LUSE[y]=LUSE[y][i][j];
						}
						dump.DHM=DHM[i][j];
						dump.P=P[i][j];
						if (ICE)
						{
							dump.ICE=ICE[i][j];
						}
						if (SSNOWGL)
						{
							dump.SSNOWGL=SSNOWGL[i][j];
							dump.SLIQGL=SLIQGL[i][j];
							dump.EVRSGL=EVRSGL[i][j];
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								dump.SNOWTEMPGL=SNOWTEMPGL[i][j];
								dump.SURFTEMPGL=SURFTEMPGL[i][j];
								dump.EBALGL=EBALGL[i][j];
							}
						}
						if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
						{
							dump.SNOWTEMP=SNOWTEMP[i][j];
							dump.SURFTEMP=SURFTEMP[i][j];
							dump.EBAL=EBAL[i][j];
							dump.RLF=RLF[i][j];
							dump.CL=CL[i][j];
							dump.RGEX=RGEX[i][j];
						}
						if ((ICE)&&(glacmethode>10))
						{
							memmove(&(dump.FIRN[0]		),&(FIRN[i][j][0]					),sizeof(double)*(firn_levels+1));
						}
						if (lake_active)
						{
							dump.LAKE=LAKE[i][j];
							dump.POND=POND[i][j];
						}
						if ((methode != 1) && (methode != 8) && (methode != 9) && (methode != 10))
						{
							dump.U=U[i][j];
						}
						if (EXKORR)
						{
							dump.EXKORR=EXKORR[i][j];
						}
						dump.stunde = stunde;
						dump.tag=tag;
						dump.monat=monat;
						dump.jahr=jahr;
					}
				}


				if (bLayeredSnowActive)
				{
					refSWEL[i][j][0]=0; // reset the sum value of liquid snow storage since it is used now as transportation vehicle to communicate the total melt to unsatzon in case a melt method was used without energy balance
#ifdef _DEBUG
					if (bHeatTransferActive)
					{
						for (int k=1;k<=LayeredSnowProperties[1].getLevels();k++)
						{
							if (refTemperSnow[i][j][k]!=refTemperSnow[i][j][k])
							{
								(*WSstream) << "ERROR: NaN value found for snow temperature at row "<<i<<", column "<<j<<" with other values:\n";
								for (int z=1;z<=LayeredSnowProperties[1].getLevels();z++)
								{
									(*WSstream) <<"refTemperSnow["<<z<<"] = "<<refTemperSnow[i][j][z];
								}
								(*WSstream) <<"[end-entry]";
								throw(2001);
							}
						}
					}
#endif
				}
				
				// now updating the SWES stack, if present: remove from the SWES-stack the amount of snow evaporation if the melt method was a simple method (non-energy balance)
				// This was computed for non-energy melt methods in evaporation model and subtracted from normal bulk snow storage in verdunstungsmodell->snow() already
				// i.e. for energy balance methods EVRS will always be 0 since the EVRS computation in snow model is skipped
				if (bLayeredSnowActive && (EVRS[i][j] > 0))
				{
					UpdateSWEStack(i, j, 1.0, -EVRS[i][j], 0, _onGlacierFalse_, _fromRefreezingFalse_,T[i][j]); 
				}

				QSNOW[i][j] = 0.0;
				if (QSNOWGL) QSNOWGL[i][j] = 0;
				QSNOW_OUT[i][j] = 0.0;
				M = 0.0;        // wichtig fuer spaetere Bilanzierung
				Mgl = 0.0;

				if ((SSNOW[i][j] + SLIQ[i][j]) > 2.5)
				{
					// Ansatz aus Corps of Engineers (1956) 
					if (T[i][j]>0)
						Albedo[i][j] = ALBsnow_min + (ALBsnow_max - ALBsnow_min)*exp(-0.05 * SNOW_AGE[i][j]);
					else
						Albedo[i][j] = ALBsnow_min + (ALBsnow_max - ALBsnow_min)*exp(-0.12 * SNOW_AGE[i][j]);
				}
				else
				{
					// ueber die Landnutzungstabelle wird hier bei abgetauter Schneedecke die ALbedo fuer Firn oder Eis eingesetzt 
					Albedo[i][j] = LUtable[(int)(LUSE[startlayer][i][j])].GetAlbedo(WSstream, tg, mt, jr, DHM[i][j], 1, i, j);
				}

				if (ICE) // is there a glacier grid at all? if so, then we have to take glaciers into account
				{
					// Gletscherschmelzmodell: these values need to be set to 0 because otherwise cells that will be not glacierized after glacier resizing 
					// keep draining with that amount for ever - leading to constant or rising runoff values over time
					QFIRN[i][j] = 0.0;
					QICE[i][j] = 0.0;
					ICEOLD[i][j] = ICE[i][j];
					// old: do it for ice or firn cells
					//		  if ((ICE[i][j]>=1)&&(ICE[i][j]<=2)) /* 1 for ice, 2 for firn 
					// new: basically the same procedure, but glacier cells are now identified by an ice value of larger than zero
					if (icevalue>0.0) // performance is better when we can skip this section for pure ice cells
					{
						// to put this cells out of the statistics
						//SNOWTEMP[i][j]=nodata;
						//SURFTEMP[i][j]=nodata;

						double moltenFirn = 0;
						double moltenSnow = 0;

						if ((glacmethode == 1) || (glacmethode == 2) || (glacmethode == 11) || (glacmethode == 12)) // traditionally, the melt routines for glaciers were implemented separately
						{
							// aus Funktion: calc_Accu_and_Melt entnommen (s.o.)
							// Schneeakkumulation Aufteilen des Niederschlages
							if (T[i][j] <= (T0R + transzone))
							{
								if (T[i][j]>(T0R - transzone))
									schneeanteil = ((T0R + transzone) - T[i][j]) / (transzone + transzone);
								else
									schneeanteil = 1.0;
							}
							else
								schneeanteil = 0.0;
							//Sliq_max_old = SSNOW[i][j] * CWH;
							// Niederschlagsstatistik
							P_SNOW[i][j] = schneeanteil*P[i][j];
							P_RAIN[i][j] = (1 - schneeanteil)*P[i][j];

							// Akkumulation 
							SSNOW[i][j] += schneeanteil*P[i][j];
							SLIQ[i][j] += (1 - schneeanteil)*P[i][j];
							if (SSNOWGL)
							{
								// for glaciers, separate storages are used
								SSNOWGL[i][j] += schneeanteil*P[i][j];
								SLIQGL[i][j] += (1 - schneeanteil)*P[i][j];
							}

							// now updating the SWES stack, if present: add solid precip to SWES layer(s)
							if (bLayeredSnowActive && (schneeanteil*P[i][j] > 0))
							{
								UpdateSWEStack(i, j, 1.0, schneeanteil*P[i][j], 0, _onGlacierFalse_, _fromRefreezingFalse_,T[i][j]); // M needs to be converted to a negative number since we need to subtract the snow fomr the storage
							}

							// snow melt off glacier
							if ((T[i][j] > T0) && (SSTO[i][j] > 0))
							{
								if ((glacmethode == 1) || (glacmethode == 11))		{
									M = DDFsnow[subbasin] *(T[i][j] - T0) / InVtag;
								}
								else { // if ((glacmethode ==2) ||(glacmethode ==12))
									M = (MF[subbasin] / InVtag + Asnow_akt[subbasin] *RAD[i][j])*(T[i][j] - T0);
								}
							}
							if (M > SSNOW[i][j]) M = SSNOW[i][j]; // es kann nur das schmelzen, was auch da ist!
							SLIQ[i][j] += M;
							SSNOW[i][j] -= M;
							Sliq_max_new = SSNOW[i][j] * CWH;  // Groesse des neuen Wasserspeichers im Schnee 
							if (Sliq_max_new < SLIQ[i][j]) // Wasserhaltekapazitaet des Schnees erschoepft? There will be real outflow out of the snow
							{
								QSNOW[i][j] = SLIQ[i][j] - Sliq_max_new;  // Abfluss ist Ueberschuss des Wasserspeichers
								if (M > 0) QSNOW_OUT[i][j] = M;
								SLIQ[i][j] = Sliq_max_new;             // neuer Wasserspeicher ist gefuellt
							}

							// now updating the SWES stack, if present: remove from the SWES-stack the amount of melt water (off glacier, when using methods 1, 2, 11 or 12, but not 13)
							if (bLayeredSnowActive && (M > 0))
							{
								UpdateSWEStack(i, j, 1.0, -M, 0, _onGlacierFalse_, _fromRefreezingFalse_,T[i][j]); // M needs to be converted to a negative number since we need to subtract the snow fomr the storage
								refSWEL[i][j][0] += M; // addig the pure melt to layer 0 of the SWEL stack. This value is taken in unsatzonmodel and added to the P_RAIN input as total input into the uppermost snow layer.
							}

							// snow melt on glacier
							if ((T[i][j] > T0) && (SSTOGL[i][j] > 0))
							{
								if ((glacmethode == 1) || (glacmethode == 11))		{
									Mgl = DDFsnow[subbasin] *(T[i][j] - T0) / InVtag;
								}
								else { // if ((glacmethode ==2) ||(glacmethode ==12))
									Mgl = (MF[subbasin] / InVtag + Asnow_akt[subbasin] *RAD[i][j])*(T[i][j] - T0);
								}
							}
							if (Mgl > SSNOWGL[i][j]) Mgl = SSNOWGL[i][j]; // es kann nur das schmelzen, was auch da ist!
							SLIQGL[i][j] += Mgl;
							SSNOWGL[i][j] -= Mgl;
							Sliq_max_new = SSNOWGL[i][j] * CWH;  // Groesse des neuen Wasserspeichers im Schnee 
							if (Sliq_max_new < SLIQGL[i][j]) // Wasserhaltekapazitaet des Schnees erschoepft? There will be real outflow out of the snow
							{
								QSNOWGL[i][j] = SLIQGL[i][j] - Sliq_max_new;  // Abfluss ist Ueberschuss des Wasserspeichers
								//								if (Mgl > 0) QSNOW_OUT[i][j] = SLIQGL[i][j]-Sliq_max_new; --> not needed for pure complete outflow on glaciers. Is part of glacier mass balance already
								SLIQGL[i][j] = Sliq_max_new;             // neuer Wasserspeicher ist gefuellt
							}


						}
						else // glacmethode==13: use snow model method
						{
							// first for the unglacierized part
							double snowtemp = 0;
							double surftemp = 0;
							double etrs = 0;
							double energybalance = 0;

							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								snowtemp = SNOWTEMP[i][j];
								surftemp = SURFTEMP[i][j];
								etrs = EVRS[i][j];
								if (bHeatTransferActive)
								{
									energybalance = EBAL[i][j];
								}
							}

							if (InV > 720.0)
							{
								// Aufruf der Routine zur Berechnung der Scheeakkumulation bzw. -schmelze
								// fuer den hellen und den dunklen Teil des Zeitintervalls. Bei Intervalllaenge
								// 1440 min wird dabei die Temperaturamplitude zwischen Tag und Nacht (aus
								// dem Verdunstungsmodell) angewandt. 
								
								double local_etrs = 0;
								double local_energybalance = 0;
								M = calc_Accu_and_Melt(i, j, Tnacht[i][j], Anteil_dunkel, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, 0, false);
								local_etrs = etrs;
								local_energybalance = energybalance;
								M += calc_Accu_and_Melt(i, j, T[i][j], Anteil_hell, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, 0, false);
								etrs += local_etrs;
								energybalance += local_energybalance;
							}
							else
							{
								// bei kuerzeren Zeitintervallen wird eh keine Trennung in Tag- und Nachttemperatur
								// vorgenommen, weshalb der doppelte Aufruf erspart werden kann 
								M = calc_Accu_and_Melt(i, j, T[i][j], 1.0, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, 0, false);
							}

							// return values will be used to update the entities depending on used methods
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								SNOWTEMP[i][j] = snowtemp;
								SURFTEMP[i][j] = surftemp;
								EVRS[i][j] = etrs;
								if (bHeatTransferActive)
								{
									EBAL[i][j] = energybalance;
								}
							}

							if (M > SSNOW[i][j]) M = SSNOW[i][j]; // es kann nur das schmelzen, was auch da ist!
							SLIQ[i][j] += M;
							SSNOW[i][j] -= M;
							Sliq_max_new = SSNOW[i][j] * CWH;  // Groesse des neuen Wasserspeichers im Schnee 
							if (Sliq_max_new < SLIQ[i][j]) // Wasserhaltekapazitaet des Schnees erschoepft? There will be real outflow out of the snow
							{
								QSNOW[i][j] = SLIQ[i][j] - Sliq_max_new;  // Abfluss ist Ueberschuss des Wasserspeichers
								if (M > 0) QSNOW_OUT[i][j] = M;
								SLIQ[i][j] = Sliq_max_new;             // neuer Wasserspeicher ist gefuellt
							}

							// now updating the SWES stack, if present: remove from the SWES-stack the amount melt water for glacier melt method method 13 (which is like glacier method 12 combined with the normal off-glacier snow melt method)
							if (bLayeredSnowActive && (M > 0))
							{
								UpdateSWEStack(i, j, 1.0, -M, 0, _onGlacierFalse_, _fromRefreezingFalse_,T[i][j]); // M needs to be converted to a negative number since we need to subtract the snow fomr the storage
								refSWEL[i][j][0] += M; // addig the pure melt to layer 0 of the SWEL stack. This value is taken in unsatzonmodel and added to the P_RAIN input as total input into the uppermost snow layer.
							}

							// now again for the glacierized part
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								snowtemp = SNOWTEMPGL[i][j];
								surftemp = SURFTEMPGL[i][j];
								etrs = EVRSGL[i][j];
								if (bHeatTransferActive)
								{
									energybalance = EBALGL[i][j];
								}
							}

							if (InV > 720.0)
							{
								// Aufruf der Routine zur Berechnung der Scheeakkumulation bzw. -schmelze
								// fuer den hellen und den dunklen Teil des Zeitintervalls. Bei Intervalllaenge
								// 1440 min wird dabei die Temperaturamplitude zwischen Tag und Nacht (aus
								// dem Verdunstungsmodell) angewandt. 
								double local_etrs = 0;
								double local_energybalance = 0;
								Mgl = calc_Accu_and_Melt(i, j, Tnacht[i][j], Anteil_dunkel, schneeanteil, SSNOWGL[i][j], SLIQGL[i][j], SSTOGL[i][j], snowtemp, surftemp, etrs, energybalance, 0, true);
								local_etrs = etrs;
								local_energybalance = energybalance;
								Mgl += calc_Accu_and_Melt(i, j, T[i][j], Anteil_hell, schneeanteil, SSNOWGL[i][j], SLIQGL[i][j], SSTOGL[i][j], snowtemp, surftemp, etrs, energybalance, 0, true);
								etrs += local_etrs;
								energybalance += local_energybalance;
							}
							else
							{
								// bei kuerzeren Zeitintervallen wird eh keine Trennung in Tag- und Nachttemperatur
								// vorgenommen, weshalb der doppelte Aufruf erspart werden kann 
								Mgl = calc_Accu_and_Melt(i, j, T[i][j], 1.0, schneeanteil, SSNOWGL[i][j], SLIQGL[i][j], SSTOGL[i][j], snowtemp, surftemp, etrs, energybalance, 0, true);
							}
							// return values will be used to update the entities depending on used methods
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								SNOWTEMPGL[i][j] = snowtemp;
								SURFTEMPGL[i][j] = surftemp;
								EVRSGL[i][j] = etrs;
								if (bHeatTransferActive)
								{
									EBALGL[i][j] = energybalance;
								}
							}

							if (Mgl > SSNOWGL[i][j]) Mgl = SSNOWGL[i][j]; // es kann nur das schmelzen, was auch da ist!
							SLIQGL[i][j] += Mgl;
							SSNOWGL[i][j] -= Mgl;
							Sliq_max_new = SSNOWGL[i][j] * CWH;  // Groesse des neuen Wasserspeichers im Schnee 
							if (Sliq_max_new < SLIQGL[i][j]) // Wasserhaltekapazitaet des Schnees erschoepft? There will be real outflow out of the snow
							{
								QSNOWGL[i][j] = SLIQGL[i][j] - Sliq_max_new;  // Abfluss ist Ueberschuss des Wasserspeichers
								// if (Mgl > 0) QSNOW_OUT[i][j] = SLIQGL[i][j]-Sliq_max_new; -> not needed for pure glacier mass balance
								SLIQGL[i][j] = Sliq_max_new;             // neuer Wasserspeicher ist gefuellt
							}

							// now updating the SWES stack, if present: add to or remove from the SWES-stack the amount of latent heat flux (evaporation or sublimation) 
							// not possible on a glacier since we do not use the layered snow model on a glacier (no unsatzon model active there)
							// if (bLayeredSnowActive && (M > 0))
							// {
							// 	UpdateSWEStack(i, j, 1.0, M, 0, _onGlacierTrue_,_fromRefreezingFalse_,T[i][j]);
							//	refSWEL[i][j][0] += M; // addig the pure melt to layer 0 of the SWEL stack. This value is taken in unsatzonmodel and added to the P_RAIN input as total input into the uppermost snow layer.
							// }

						}
						// Achtung; Zufluss zu firn, schnee oder eis speicher richtet sich danach, welche Unterlage unter dem Schnee ist:
						// Schmelze aus Schnee über Firn geht in den Firnspeicher (wird also eigentlich firnschmelze)
						// Schmelze aus Schnee ueber Eis bleibt aber im Schneespeicher.
						// removed on 12.12.2014; Qsnow is added to the firn and qsnowonglac statistics only since it must be left constant for the unglacierized fraction for the soil model as precipitation input
						//						if (((glacmethode>10)&&(FIRN1[i][j]>0.0))||(ICE[i][j]==2))
						//						{
						//							// auf Firnflächen wird die Schneeschmelze als Firnschmelze betrachtet
						//							QFIRN[i][j] = QSNOW[i][j]; // auf Firnflächen wird die Schneeschmelze als Firnschmelze betrachtet
						//							QSNOW[i][j]=0; // und dafuer gibt es keine Schneeschmelze an diesen Flaechen -> wichtig fuer
						//							// Parallelspeicherkonzept, Firn ist am langsamsten!! und dominiert daher bei einer Schnee-auf-Firn-Situation
						//							QSNOW_OUT[i][j]=0;
						//						}

						moltenSnow = Mgl; // in order to limit the total melt to the maximum of snowmelt, icemelt and firnmelt

						// Es liegt kein Schnee mehr und es ist noch nichts geschmolzen worden -> Schmelze von Firn oder Eis
						if ((T[i][j] > T0) && (SSNOWGL[i][j] <= 0))
						{
							double debrisfactor = 1.0;
							if (GLACDEB != NULL)
							{
								if (GLACDEB[i][j]>nodata + 1)
								{
									debrisfactor = GLACDEB[i][j];
								}
							}
							if (glacmethode == 1)  // einfaches T-Index-Verfahren
							{
								if (ICE[i][j] == 1)  // Eisflaechen 
								{
									Mgl = dGloabDebrisScaling*debrisfactor*DDFice[subbasin] *(T[i][j] - T0) / InVtag - moltenSnow;
									if (Mgl>0) QICE[i][j] = Mgl;
								}
								if (ICE[i][j] == 2)
								{
									Mgl = DDFfirn[subbasin] *(T[i][j] - T0) / InVtag - moltenSnow;    // Firnflaechen
									if (Mgl>0) QFIRN[i][j] += Mgl; // since snow melt was already converted to firn melt
								}
							}
							else if (glacmethode == 2) {  // Verfahren Hock 3 (T-Index mit Korrektur durch korrigierte Globalstrahlung)
								if (ICE[i][j] == 1)
								{
									Mgl = dGloabDebrisScaling*debrisfactor*(MF[subbasin] / InVtag + Aice_akt[subbasin] *RAD[i][j])*(T[i][j] - T0) - moltenSnow;
									if (Mgl>0) QICE[i][j] = Mgl;
								}
								if (ICE[i][j] == 2)
								{
									Mgl = (MF[subbasin] / InVtag + Asnow_akt[subbasin] *RAD[i][j])*(T[i][j] - T0) - moltenSnow;
									if (Mgl>0) QFIRN[i][j] += Mgl;
								}
							}
							else if (glacmethode == 11) { // simple temperature method but for dynamic glacier model
								if (FIRN1[i][j] <= 0.0) // ice melt only
								{
									Mgl = dGloabDebrisScaling*debrisfactor*DDFice[subbasin] *(T[i][j] - T0) / InVtag - moltenSnow;
									if (Mgl>0) QICE[i][j] = Mgl;
								}
								else // firn melt --> FIRN1[i][j]>0
								{
									Mgl = DDFfirn[subbasin] *(T[i][j] - T0) / InVtag - moltenSnow;    // Firnflaechen
									if (Mgl>0) QFIRN[i][j] += Mgl; else Mgl = 0;
									moltenFirn = Mgl;
								}
							}
							else if ((glacmethode == 12) || (glacmethode == 13)) { // Approach Hock 3 (t-index with correction by global radiation) but for dynamic glacier model
								if (FIRN1[i][j] <= 0.0) // ice melt only
								{
									Mgl = dGloabDebrisScaling*debrisfactor*(MF[subbasin] / InVtag + Aice_akt[subbasin] *RAD[i][j])*(T[i][j] - T0) - moltenSnow;
									if (Mgl>0) QICE[i][j] = Mgl;
								}
								else // firn melt
								{
									Mgl = (MF[subbasin] / InVtag + Asnow_akt[subbasin] *RAD[i][j])*(T[i][j] - T0) - moltenSnow;
									if (Mgl>0) QFIRN[i][j] += Mgl; else Mgl = 0;
									moltenFirn = Mgl;
								}
							}
							if (glacmethode>10) // for dynamic glacier model only: balancing firn layers
							{
								for (int lv = 1; lv <= firn_levels; lv++)
								{
									// when the glacier grows it may happen, that not all firn cells are initialized correctly in some test cases...
									if (FIRN[i][j][lv] == nodata) FIRN[i][j][lv] = 0;

									FIRN[i][j][lv] -= moltenFirn;
									if (FIRN[i][j][lv]<0.0) // the actual firn layer does not contain as much water as required by the melted amount
									{
										moltenFirn = -FIRN[i][j][lv]; // remaining melt for lower layers
										FIRN[i][j][lv] = 0;
									}
									else // balance ok
									{
										moltenFirn = 0;
										break;
									}
								}
								if (moltenFirn>0) // still some melt left --> let it be ice melt
								{
									QICE[i][j] += moltenFirn;
								}
								// balancing FIRN1 grid to have a newly calculated statistics for each time step
								FIRN1[i][j] = 0;
								for (int lv = 1; lv <= firn_levels; lv++)
								{
									FIRN1[i][j] += FIRN[i][j][lv];
								}
								FIRN[i][j][0] = FIRN1[i][j];
							}
						}
					}
				}

				// old:		icevalue=ICE?ICE[i][j]:nodata;
				//		icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				//		if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				//  	    else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				// normal snow areas without ice and firn below
				if (icevalue == 0.0) // performance is better when we can skip this section for non-iced cells
					//		if ((!ICE)||((ICE)&&(icevalue==nodata)))
				{
					for (int z = 1; z <= tracer_count; z++) // radioactive decay for radioactive tracers in the storage 
						if (tracer_info[z].radioactive)
							tracer_module[z].concstor[i][j] *= tracer_info[z].decfact;

					// $kf 2016-07-20: get snow (solid) throughfall from interception model if applicable
					double snow_throughfall = 0.;
					if (SCN_snow_throuhfall)
					{
						snow_throughfall = this->SCN_snow_throuhfall[i][j];
					}

					if (P[i][j]>0)
					{
						for (int z = 1; z <= tracer_count; z++)
						{ // updating tracer concentrations in the storages 
							tracer_module[z].concstor[i][j] = tracer_info[z].evaporating ?
								mixing(tracer_module[z].concstor[i][j], tracer_module[z].cinput[i][j],
								SSNOW[i][j] + SLIQ[i][j], P[i][j]) :
								tracer_module[z].concstor[i][j] + tracer_module[z].cinput[i][j] * P[i][j];
						}
					}

					// first for the unglacierized part
					double snowtemp = 0;
					double surftemp = 0;
					double etrs = 0;
					double energybalance = 0;

					if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
					{
						snowtemp = SNOWTEMP[i][j];
						surftemp = SURFTEMP[i][j];
						etrs = EVRS[i][j];
						if (bHeatTransferActive)
						{
							energybalance = EBAL[i][j];
						}
					}

					if (InV > 720.0)
					{
						// Aufruf der Routine zur Berechnung der Scheeakkumulation bzw. -schmelze
						// fuer den hellen und den dunklen Teil des Zeitintervalls. Bei Intervalllaenge
						// 1440 min wird dabei die Temperaturamplitude zwischen Tag und Nacht (aus
						// dem Verdunstungsmodell) angewandt. 
						double local_etrs=0;
						double local_energybalance = 0;
						double local_snow_throughfall = 0;
						M = calc_Accu_and_Melt(i, j, Tnacht[i][j], Anteil_dunkel, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, snow_throughfall, false);
						local_etrs = etrs;
						local_energybalance = energybalance;
						local_snow_throughfall = snow_throughfall;
						M += calc_Accu_and_Melt(i, j, T[i][j], Anteil_hell, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, snow_throughfall, false);
						etrs+= local_etrs;
						energybalance+= local_energybalance;
						snow_throughfall+= local_snow_throughfall;
					}
					else
					{
						// bei kuerzeren Zeitintervallen wird eh keine Trennung in Tag- und Nachttemperatur
						// vorgenommen, weshalb der doppelte Aufruf erspart werden kann 
						M = calc_Accu_and_Melt(i, j, T[i][j], 1.0, schneeanteil, SSNOW[i][j], SLIQ[i][j], SSTO[i][j], snowtemp, surftemp, etrs, energybalance, snow_throughfall, false);
					}

					// return values will be used to update the entities depending on used methods
					if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
					{
						SNOWTEMP[i][j] = snowtemp;
						SURFTEMP[i][j] = surftemp;
						EVRS[i][j] = etrs;
						if (bHeatTransferActive)
						{
							EBAL[i][j] = energybalance;
						}
					}

					// nun wird bilanziert: es gilt, den Schneespeicher SSNOW, den Fluessigspeicher
					// SLIQ und die berechnete Schmelze miteinander in Einklang zu bringen
					// (SLIQ ist max. CWH*SSNOW). P spielt keine Rolle mehr, da er am Anfang der
					// Routine auf SSNOW und SLIQ aufgeteilt wurde! 

					if (M > SSNOW[i][j]) M = SSNOW[i][j]; // es kann nur das schmelzen, was auch da ist! 
					SLIQ[i][j] += M;
					SSNOW[i][j] -= M;
					Sliq_max_new = SSNOW[i][j] * CWH;  // Groesse des neuen Wasserspeichers im Schnee 
					if (Sliq_max_new < SLIQ[i][j]) // Wasserhaltekapazitaet des Schnees erschoepft? 
					{
						QSNOW[i][j] = SLIQ[i][j] - Sliq_max_new;  // Abfluss ist Ueberschuss des Wasserspeichers 
						if (M > 0) QSNOW_OUT[i][j] = M;
						SLIQ[i][j] = Sliq_max_new;             // neuer Wasserspeicher ist gefuellt 
					}

					if (SSNOWGL && SLIQGL && SSTOGL)
					{
						// for absolutely unglacierized cells, all three snow storage components for the glacierized fraction must be 0 (when a cell lost all glaciers, thius is especially important)
						SSNOWGL[i][j] = 0;
						SLIQGL[i][j] = 0;
						SSTOGL[i][j] = 0;
					}

					// now updating the SWES stack, if present: remove from the SWES-stack the amount of melt water like we did for the bulk snow pack above
					if (bLayeredSnowActive && (M > 0))
					{
						UpdateSWEStack(i, j, 1.0, -M, 0, _onGlacierFalse_, _fromRefreezingFalse_,T[i][j]); // M needs to be converted to a negative number since we need to subtract the snow from the storage, UpdateSWEStack will handle positive and negative inputs differently
						refSWEL[i][j][0] += M; // addig the pure melt to layer 0 of the SWEL stack. This value is taken in unsatzonmodel and added to the P_RAIN input as total input into the uppermost snow layer.
					}

					// Niederschlagsstatistik 
					// do not update here, since in calc_Accu_and_Melt this was already done and P[i][j] was changed for the factor of the snow redistribution grid (as well as P_SNOW is based on thatgrid as well already)
					// P_SNOW[i][j] = schneeanteil*P[i][j];
					// P_RAIN[i][j] = (1 - schneeanteil)*P[i][j];

					// balancing tracer concentrations for output-flux 
					double temp = QSNOW[i][j] + SLIQ[i][j] + SSNOW[i][j];
					for (int z = 1; z <= tracer_count; z++)
					{
						if (tracer_info[z].evaporating != 0)
							tracer_module[z].concflux[i][j] = tracer_module[z].concstor[i][j];
						else
						{ // subtracting massflux from storage and add it to the outflow 
							trms = (temp)>0.0 ? tracer_module[z].concstor[i][j] * QSNOW[i][j] / temp : 0.0;
							trmsmax = tracer_info[z].max_conc*QSNOW[i][j];
							if (trmsmax<trms) trms = trmsmax;
							tracer_module[z].concstor[i][j] -= trms;
							tracer_module[z].concflux[i][j] = QSNOW[i][j]>0 ? trms / QSNOW[i][j] : 0;
						}
					}
				}

				if ((SSNOW[i][j] + SLIQ[i][j]) > 5.0)
					SNOWDAYS[i][j] = (1 / InVtag);
				else
					SNOWDAYS[i][j] = 0;

				if ((T[i][j] <= (T0R + transzone)) && (P[i][j]>1.0) && ((SSNOW[i][j] + SLIQ[i][j]) > 5.0))
					SNOW_AGE[i][j] = 0;
				else
					SNOW_AGE[i][j] += (1 / InVtag);

				if ((SSNOW[i][j] + SLIQ[i][j]) > 5.0)
				{
					// Ansatz aus Corps of Engineers (1956) 
					if (T[i][j]>0)
						Albedo[i][j] = ALBsnow_min + (ALBsnow_max - ALBsnow_min)*exp(-0.05 * SNOW_AGE[i][j]);
					else
						Albedo[i][j] = ALBsnow_min + (ALBsnow_max - ALBsnow_min)*exp(-0.12 * SNOW_AGE[i][j]);
				}
				else
				{
					// ueber die Landnutzungstabelle wird hier bei abgetauter Schneedecke die ALbedo fuer Firn oder Eis eingesetzt 
					Albedo[i][j] = LUtable[(int)(LUSE[startlayer][i][j])].GetAlbedo(WSstream, tg, mt, jr, DHM[i][j], 1, i, j);
				}

				SSTO[i][j] = SSNOW[i][j] + SLIQ[i][j];
				if (SSTOGL)
				{
					SSTOGL[i][j] = SSNOWGL[i][j] + SLIQGL[i][j];
				}

				// modelling thaw depth
				if (SNOWFREEDAYS) // only if SNOWFREEDAYS exists
				{
					if ((SSNOW[i][j] + SLIQ[i][j]) > MaxSWE_SnowFree)
					{
						// if timestep is not 1 day, InvTag will be the number of intervals per day,
						// so after one day, the sum of all intervals of this day will be 1.0
						SNOWCOVERDAYS[i][j] += (1.0 / InVtag);
					}
					else
					{ // no snow cover
						SNOWCOVERDAYS[i][j] = 0;
					}

					// snow free days: will increase if there is no snow cover.
					// will be reset to 0 if snow cover is older than 30 days (this is the refreezing model, very simple, but sufficient
					if (SNOWCOVERDAYS[i][j] == 0)
					{
						SNOWFREEDAYS[i][j] += (1.0 / InVtag);
						if (SNOWFREEDAYS[i][j]>0.0)
							THAWDEPTH[i][j] = ALPHA[i][j] * sqrt(SNOWFREEDAYS[i][j]);
					}
					else if (SNOWCOVERDAYS[i][j] >= MinSAgeRefreeze)
					{
						SNOWFREEDAYS[i][j] = 0;
						// finally, the soil freezes again...
						THAWDEPTH[i][j] = 0.0;
					}

					if (ROOTDEPTH[i][j]>THAWDEPTH[i][j])
						ROOTDEPTH[i][j] = THAWDEPTH[i][j];
				}

				if (this->dumpmode==true)
				{
					bool writedumpswes=false;
					bool writedumpswel=false;

					int dumplayerswes = 0;
					int dumplayerswel = 0;

					if (bLayeredSnowActive)
					{
						for (int k = 1; k < levelsSnow; k++)
						{
							if ((refSWES[i][j][k]>0)&&(refSWES[i][j][k+1]<=0))
							{
								writedumpswes = true;
								dumplayerswes = k;
								break;
							}
						}

						for (int k = 1; k < levelsSnow; k++)
						{
							if ((refSWEL[i][j][k]>0)&&(refSWES[i][j][k]<=0))
							{
								writedumpswel = true;
								dumplayerswel = k;
								break;
							}
						}
					}

					if (writedumpswes || writedumpswel)
					{
						if ((actDumpFiles < maxDumpFiles) || maxDumpFiles==0)
						{
							std::ofstream dumpstream;
							std::stringstream fname;
							fname << "MPI_" << mpi_rank <<"_";
							if (writedumpswes)
							{
								fname << "swes(" << dumplayerswes << ")_";
							}
							if (writedumpswel)
							{
								fname << "swel(" << dumplayerswel << ")_";
							}

							fname << "dump_row_" << i << "_col_" << j << "_" << jr << setw(2) << setfill('0') << mt << setw(2) << setfill('0') << tg << setw(2) << setfill('0') << se << ".dmp";
							std::string sname(fname.str());
							dumpstream.open(sname.c_str(), std::ios::binary);

							if (bLayeredSnowActive)
							{
								dumpstream.write((char *)(&refSWES[i][j][0]),			sizeof(double)*(1+levelsSnow));
								dumpstream.write((char *)(&refSWEL[i][j][0]),			sizeof(double)*(1+levelsSnow));
								dumpstream.write((char *)(&refThickSnow[i][j][0]),		sizeof(double)*(1+levelsSnow));
								dumpstream.write((char *)(&refDensSnow[i][j][0]),		sizeof(double)*(1+levelsSnow));
								dumpstream.write((char *)(&refUpperIndexSnow[i][j]),	sizeof(double));
							}
							dumpstream.write((char *)(&EVRS[i][j]),			sizeof(double));
							dumpstream.write((char *)(&SSNOW[i][j]),			sizeof(double));
							dumpstream.write((char *)(&SLIQ[i][j]),			sizeof(double));
							dumpstream.write((char *)(&SSTO[i][j]),			sizeof(double));
							dumpstream.write((char *)(&T[i][j]),				sizeof(double));
							dumpstream.write((char *)(&Albedo[i][j]),		sizeof(double));
							for (int y=1;y<=iNumVegetationLayers;y++) 
							{
								dumpstream.write((char *)(&(LUSE[y][i][j])),sizeof(double));
							}
							dumpstream.write((char *)(&DHM[i][j]),			sizeof(double));
							dumpstream.write((char *)(&P[i][j]),				sizeof(double));
							if (ICE)
							{
								dumpstream.write((char *)(&ICE[i][j]),			sizeof(double));
							}
							if (SSNOWGL)
							{
								dumpstream.write((char *)(&SSNOWGL[i][j]),		sizeof(double));
								dumpstream.write((char *)(&SLIQGL[i][j]),		sizeof(double));
								dumpstream.write((char *)(&EVRSGL[i][j]),		sizeof(double));
								if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
								{
									dumpstream.write((char *)(&SNOWTEMPGL[i][j]),	sizeof(double));
									dumpstream.write((char *)(&SURFTEMPGL[i][j]),	sizeof(double));
									dumpstream.write((char *)(&EBALGL[i][j]),		sizeof(double));
								}
							}
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								dumpstream.write((char *)(&SNOWTEMP[i][j]),		sizeof(double));
								dumpstream.write((char *)(&SURFTEMP[i][j]),		sizeof(double));
								dumpstream.write((char *)(&EBAL[i][j]),			sizeof(double));
								dumpstream.write((char *)(&RLF[i][j]),			sizeof(double));
								dumpstream.write((char *)(&CL[i][j]),			sizeof(double));
								dumpstream.write((char *)(&RGEX[i][j]),			sizeof(double));
							}
							if ((ICE)&&(glacmethode>10))
							{
								dumpstream.write((char *)(&FIRN[i][j][0]),		sizeof(double)*(firn_levels+1));
							}
							if (lake_active)
							{
								dumpstream.write((char *)(&LAKE[i][j]),			sizeof(double));
								dumpstream.write((char *)(&POND[i][j]),			sizeof(double));
							}
							if ((methode != 1) && (methode != 8) && (methode != 9) && (methode != 10))
							{
								dumpstream.write((char *)(&U[i][j]),				sizeof(double));
							}
							if (EXKORR)
							{
								dumpstream.write((char *)(&EXKORR[i][j]),		sizeof(double));
							}
							dumpstream.write((char *)(&me), sizeof(int));
							dumpstream.write((char *)(&se), sizeof(int));
							dumpstream.write((char *)(&tg),sizeof(int));
							dumpstream.write((char *)(&mt),sizeof(int));
							dumpstream.write((char *)(&jr),sizeof(int));

							dumpstream.close();
						}
						actDumpFiles++;


						if (closeAfterDump) 	{
								(*WSstream) << "FATAL: \n\n -------- Snowm odel: command line option -closeAfterDump active! ---------------------------------------------\n\n"
								<< " a dumpfile has been written and WaSiM aborts now.\n"
								<<      "-------------------------------------------------------------\n\n" <<"[end-entry]";

								throw(199);
						}

						if ((actDumpFiles < maxDumpFiles)&&(repairAfterDump))
						{
							(*WSstream) << "DEBUG: \n ---- Snow model: command line option -repairAfterDump active! --------\n"
								<< "nothing to here" << "[end-entry]";
						}
						else if (actDumpFiles == maxDumpFiles)
						{
							(*WSstream) << "DEBUG: \n ---- Snow model: command line option -dumpmode:<maxDumpFiles> with maxDumpFiles="<<maxDumpFiles<<" active! --------\n"
								<< "\n "<<maxDumpFiles<<" has now be written. Continuing from now on without writing dump files for this instance and MPI-rank..." 
								<< "[end-entry]";
						}
					}
					if (bLayeredSnowActive)
					{
						delete[] dump.refSWES;
						delete[] dump.refSWEL;
						delete[] dump.refThickSnow;
						delete[] dump.refDensSnow;
					}
					if ((ICE)&&(glacmethode>10))
					{
						delete[] dump.FIRN;
					}
					delete[] dump.LUSE;
				}
			}
		}
	}
}

// here, the old run method was cut into two partsin order to allow another sub model (in this case unsatzon model) to manipulate snow model entities like snow pack itself and also the melt
// process. 
void schneemodellklasse::run2(long me,long se, long tg, long mt, long jr)
{

	// gibt lediglich eine Nachricht aus
	teilmodellklasse::run(me,se, tg, mt, jr,NULL,NULL);

	if (do_it == 0) return;

	// when in SpinUp mode, all cells that were not computed must be filled with the HRUMaster data now
#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) /* Spalten des Grids */
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
			{
				if (doSpinUp != 0 && HRUMaster[i][j] == nodata)
				{
					// getting the coordinates of the master cell
					int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
					int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
					// filling the grid cells
					SSTO[i][j] = SSTO[imaster][jmaster];
					SSNOW[i][j] = SSNOW[imaster][jmaster];
					SLIQ[i][j] = SLIQ[imaster][jmaster];
					SNOWDAYS[i][j] = SNOWDAYS[imaster][jmaster];
					SNOW_AGE[i][j] = SNOW_AGE[imaster][jmaster];
					P_SNOW[i][j] = P_SNOW[imaster][jmaster];
					P_RAIN[i][j] = P_RAIN[imaster][jmaster];
					Albedo[i][j] = Albedo[imaster][jmaster];
					QSNOW[i][j] = QSNOW[imaster][jmaster];
					QSNOW_OUT[i][j] = QSNOW_OUT[imaster][jmaster];
					if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // enbal
					{
						SNOWTEMP[i][j] = SNOWTEMP[imaster][jmaster];
						SURFTEMP[i][j] = SURFTEMP[imaster][jmaster];
					}
					if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
					{

						INPUT_MASS[i][j] = INPUT_MASS[imaster][jmaster];
						MOBILE_MASS[i][j] = MOBILE_MASS[imaster][jmaster];
						DEPOSITION[i][j] = DEPOSITION[imaster][jmaster];
					}
					if (ICE)
					{
						QFIRN[i][j] = QFIRN[imaster][jmaster];
						QICE[i][j] = QICE[imaster][jmaster];
						if (glacmethode>10)
						{
							FIRN[i][j] = FIRN[imaster][jmaster];
							GMB[i][j] = GMB[imaster][jmaster];
							GMBO[i][j] = GMBO[imaster][jmaster];
							SSNOWGL[i][j] = SSNOWGL[imaster][jmaster];
							SLIQGL[i][j] = SLIQGL[imaster][jmaster];
							SSTOGL[i][j] = SSTOGL[imaster][jmaster];
							QSNOWGL[i][j] = QSNOWGL[imaster][jmaster];
							if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
							{
								SURFTEMPGL[i][j] = SURFTEMPGL[imaster][jmaster];
								SNOWTEMPGL[i][j] = SNOWTEMPGL[imaster][jmaster];
							}
							EVRSGL[i][j] = EVRSGL[imaster][jmaster];
						}
					}
				}
			}
		}
	}

	// reset counters for weights
	//	if (glacierlist)
	//	{
	//		for (i=1;i<=zonmax_routing;i++)
	//		{
	//			if (glacierlist[i].isGlacierized)
	//			{
	//				glacierlist[i].CellCntBare=0;
	//				glacierlist[i].CellCntGlacerized =0;
	//			}
	//		}
	//	}

	// new snow model extensiond (Warscher, 2013); here: Slides
	double icevalue;
	if (doSpinUp == 0)
	{
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14))
		{
			double Flow1, Flow2, Flow3, Flow4, liq, liq1, liq2, liq3, liq4;
			int row, col;

			//		SL = gridinfo[getindex(_Gefaelle_)].pointer; // has already be done in the constructor when initializing the snow model
			int mpi_rank = 0;
			int mpi_size = 1;

#ifdef _MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
			this->MPIGatherGrid(sstoindex);
			this->MPIGatherGrid(sliqindex);
			this->MPIGatherGrid(ssnowindex);
			if (glacmethode > 10)
			{
				this->MPIGatherGrid(sstoglindex);
				this->MPIGatherGrid(sliqglindex);
				this->MPIGatherGrid(ssnowglindex);
			}

			if (mpi_rank == 0) // slides will be calculated in rank 0 only since the order of cells is completely independent on any MPI regions
			{

				for (int i = 1; i <= nrows; i++)
				{

					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if (zn[i][j] != nodata)
						{
							DEPOSITION[i][j] = 0.0;
							INPUT_MASS[i][j] = 0.0;
							MOBILE_MASS[i][j] = 0.0;
						}
						else
						{
							DEPOSITION[i][j] = nodata;
							INPUT_MASS[i][j] = nodata;
							MOBILE_MASS[i][j] = nodata;
						}
					}
				}
				// no OpenMP can be used either since the order of processing is importatnt (from highest to lowest without parallelization!)
				for (int i = 1; i <= nrows; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if (zn[i][j] != nodata)
						{
							row = (int)ELEV_ORDER_ROW[i][j];
							col = (int)ELEV_ORDER_COL[i][j];

							if (ICE&&SSTOGL) icevalue = ICE[row][col]; else icevalue = 0; // valid range of icevalue here between 0...2 (for old glacier model ICE can be 1 for ice or 2 for firn)
							if (icevalue>1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
							else if (icevalue<0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

							double slfr1, slfr2, slfr3, slfr4;
							double slfrSum = 0;

							// checking all slide fractions that they do not allow a slide to end in another sub-basin
							// north
							if ((row >     1) && (zn[row - 1][col] != nodata) && (zn[row - 1][col] == zn[row][col])) { slfr1 = SLIDEFRACTION1[row][col]; }
							else { slfr1 = 0; }
							// west
							if ((col >     1) && (zn[row][col - 1] != nodata) && (zn[row][col - 1] == zn[row][col])) { slfr2 = SLIDEFRACTION2[row][col]; }
							else { slfr2 = 0; }
							// east
							if ((col < ncols) && (zn[row][col + 1] != nodata) && (zn[row][col + 1] == zn[row][col])) { slfr3 = SLIDEFRACTION3[row][col]; }
							else { slfr3 = 0; }
							// south
							if ((row < nrows) && (zn[row + 1][col] != nodata) && (zn[row + 1][col] == zn[row][col])) { slfr4 = SLIDEFRACTION4[row][col]; }
							else { slfr4 = 0; }

							slfrSum = slfr4 + slfr3 + slfr2 + slfr1;

							// if some direction(s) is/are blocked because it leads to another sub-basin, the weights must be recalculated
							if (slfrSum>0)
							{
								slfr1 /= slfrSum;
								slfr2 /= slfrSum;
								slfr3 /= slfrSum;
								slfr4 /= slfrSum;
							}

							// compute fot the actual cell the amount of snow that will moved downwards, which is either the amount coming from another cell (in MOBILE_MASS) and/or a fraction of the actual SWE
							if ((SSTO[row][col] > 0.0) || (SSTOGL && SSTOGL[row][col] > 0.0))
							{
								if (SL[row][col] >= minSlideSlope)
								{
									if ((icevalue < 1.0) && (SSTO[row][col] > 0.0))
									{
										INPUT_MASS[row][col] = (slideMassFrac*SSTO[row][col] * (1.0 - icevalue)); // unglacierized faction of the cell can add its part to the slide
										liq = SLIQ[row][col] / SSTO[row][col]; // fraction of liquid water compared to total SWE
										SSTO[row][col] -= (slideMassFrac*SSTO[row][col]); // subtract the mobilized fraction from the SWE total 
										SLIQ[row][col] = liq * SSTO[row][col]; // new amount of liquid water in the snow pack, according to the old fraction of liquid water on SWE before subtracting the mobilized mass
										SSNOW[row][col] = SSTO[row][col] - SLIQ[row][col]; // new amount of solid snow as difference between total SWE and liquid water
									}
									if (ICE && SSTOGL && SSTOGL[row][col]>0 && icevalue > 0.0) // only if there is some ice, the slide input will also be taken from there, according to different snow pack masses on the glacier
									{
										INPUT_MASS[row][col] += (slideMassFrac*SSTOGL[row][col] * icevalue); // glacierized fraction of area may add its part to the slide
										liq = SLIQGL[row][col] / SSTOGL[row][col]; // fraction of liquid water compared to total SWE
										SSTOGL[row][col] -= (slideMassFrac*SSTOGL[row][col]); // subtract the mobilized fraction from the SWE total 
										SLIQGL[row][col] = liq * SSTOGL[row][col]; // new amount of liquid water in the snow pack, according to the old fraction of liquid water on SWE before subtracting the mobilized mass
										SSNOWGL[row][col] = SSTOGL[row][col] - SLIQGL[row][col]; // new amount of solid snow as difference between total SWE and liquid water
									}
								}
								else
								{
									//13.12.2014								INPUT_MASS[row][col]=MOBILE_MASS[row][col];
									INPUT_MASS[row][col] = 0;
								}
							}
							else
							{
								//13.12.2014							INPUT_MASS[row][col]=MOBILE_MASS[row][col];
								INPUT_MASS[row][col] = 0;
							}

							if (slfrSum>0)
							{
								Flow1 = slfr1 * (INPUT_MASS[row][col] + MOBILE_MASS[row][col]); // north
								Flow2 = slfr2 * (INPUT_MASS[row][col] + MOBILE_MASS[row][col]); // west
								Flow3 = slfr3 * (INPUT_MASS[row][col] + MOBILE_MASS[row][col]); // east
								Flow4 = slfr4 * (INPUT_MASS[row][col] + MOBILE_MASS[row][col]); // south
								MOBILE_MASS[row][col] = 0; // this mass is leaving the actual cell (but only, if slfrSum is > 0, otherwise the mobile mass cannot leave the cell and will stay there). Later, when handling the cell where this mass is flowing to, this mass can create a new mobile mass in that other cell
								//13.12.2014							Flow1 = slfr1 * INPUT_MASS[row][col]; // north
								//13.12.2014							Flow2 = slfr2 * INPUT_MASS[row][col]; // west
								//13.12.2014							Flow3 = slfr3 * INPUT_MASS[row][col]; // east
								//13.12.2014							Flow4 = slfr4 * INPUT_MASS[row][col]; // south
							}
							else
							{	// lsfrSum == 0.0, i.e. no way out of here...
								// the amount that was computed to leave the cell cannot leave the cell and must be added back to the snow pack
								// since we set INPUT_MASS to MOBILE_MASS above for cells that were either not stee enough or without snow, 
								// the following line is not required for that cells but for cells where we subtracted a certain amount of 
								// snow from the actual snow pack to INPUT_MASS in addition to the incoming MOBILE_MASS. Later in the next grid-loop,
								// all depositions and mobile masses will be added to the snow pack, so it is just fine to define the Mobile mas to 
								// be the above calculated Input Mass
								//13.12.2014							MOBILE_MASS[row][col]=INPUT_MASS[row][col];
								MOBILE_MASS[row][col] += INPUT_MASS[row][col];
								Flow1 = Flow2 = Flow3 = Flow4 = 0;
							}

							// north facing slopes (or fraction that is north facing)
							if ((Flow1>0) && (row > 1) && (zn[row - 1][col] != nodata) && (zn[row - 1][col] == zn[row][col]))
							{
								if (DEPOSITION[row - 1][col] + Flow1 <= DEPOSITION_MAX[row - 1][col])
								{
									DEPOSITION[row - 1][col] += Flow1;
								}
								else
								{
									MOBILE_MASS[row - 1][col] += Flow1 + DEPOSITION[row - 1][col] - DEPOSITION_MAX[row - 1][col];
									DEPOSITION[row - 1][col] = DEPOSITION_MAX[row - 1][col];
								}
							}

							// west facing slopes (or fraction that is west facing)
							if ((Flow2>0) && (col>1) && (zn[row][col - 1] != nodata) && (zn[row][col - 1] == zn[row][col]))
							{
								if (DEPOSITION[row][col - 1] + Flow2 <= DEPOSITION_MAX[row][col - 1])
								{
									DEPOSITION[row][col - 1] += Flow2;
								}
								else
								{
									MOBILE_MASS[row][col - 1] += Flow2 + DEPOSITION[row][col - 1] - DEPOSITION_MAX[row][col - 1];
									DEPOSITION[row][col - 1] = DEPOSITION_MAX[row][col - 1];
								}
							}

							// east
							if ((Flow3>0) && (col<ncols) && (zn[row][col + 1] != nodata) && (zn[row][col + 1] == zn[row][col]))
							{
								if (DEPOSITION[row][col + 1] + Flow3 <= DEPOSITION_MAX[row][col + 1])
								{
									DEPOSITION[row][col + 1] += Flow3;
								}
								else
								{
									MOBILE_MASS[row][col + 1] += Flow3 + DEPOSITION[row][col + 1] - DEPOSITION_MAX[row][col + 1];
									DEPOSITION[row][col + 1] = DEPOSITION_MAX[row][col + 1];
								}
							}

							// south
							if ((Flow4>0) && (row<nrows) && (zn[row + 1][col] != nodata) && (zn[row + 1][col] == zn[row][col]))
							{
								if (DEPOSITION[row + 1][col] + Flow4 <= DEPOSITION_MAX[row + 1][col])
								{
									DEPOSITION[row + 1][col] += Flow4;
								}
								else
								{
									MOBILE_MASS[row + 1][col] += Flow4 + DEPOSITION[row + 1][col] - DEPOSITION_MAX[row + 1][col];
									DEPOSITION[row + 1][col] = DEPOSITION_MAX[row + 1][col];
								}
							}

						}
					}
				}

				// finally balancing the snow storage grids
				for (int i = 1; i <= nrows; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if (zn[i][j] != nodata)
						{
							if (ICE&&SSTOGL) icevalue = ICE[i][j]; else icevalue = 0; // valid range of icevalue here between 0...2 (for old glacier model ICE can be 1 for ice or 2 for firn)
							if (icevalue>1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
							else if (icevalue<0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

							// first for the unglacierized fraction
							if (icevalue<1.0)
							{
								if (SSTO[i][j] > 0.0)
								{
									liq = SLIQ[i][j] / SSTO[i][j];
									SSTO[i][j] = SSTO[i][j] + DEPOSITION[i][j] + MOBILE_MASS[i][j];
									SLIQ[i][j] = liq * SSTO[i][j];
									SSNOW[i][j] = SSTO[i][j] - SLIQ[i][j];
								}
								else
								{   // no snow yet on that cell
									SLIQ[i][j] = 0.0;
									SSTO[i][j] = DEPOSITION[i][j] + MOBILE_MASS[i][j];
									SSNOW[i][j] = SSTO[i][j]; // no liquid fraction, so we can set ssnow to ssto
								}
							}
							if (ICE && SSTOGL && icevalue>0.0) // same for glacierized fraction (will count with weight of icevalue in the statistics!)
							{
								if (SSTOGL[i][j] > 0.0)
								{
									liq = SLIQGL[i][j] / SSTOGL[i][j];
									SSTOGL[i][j] = SSTOGL[i][j] + DEPOSITION[i][j] + MOBILE_MASS[i][j];
									SLIQGL[i][j] = liq * SSTOGL[i][j];
									SSNOWGL[i][j] = SSTOGL[i][j] - SLIQGL[i][j];
								}
								else
								{   // no snow yet on that cell
									SLIQGL[i][j] = 0.0;
									SSTOGL[i][j] = DEPOSITION[i][j] + MOBILE_MASS[i][j];
									SSNOWGL[i][j] = SSTOGL[i][j]; // no liquid fraction, so we can set ssnow to ssto
								}
							}
						}
					}
				}
			}
			this->MPIScatterGrid(sstoindex);
			this->MPIScatterGrid(sliqindex);
			this->MPIScatterGrid(ssnowindex);
			if (glacmethode > 10)
			{
				this->MPIScatterGrid(sstoglindex);
				this->MPIScatterGrid(sliqglindex);
				this->MPIScatterGrid(ssnowglindex);
			}
			this->MPIScatterGrid(input_massindex);
			this->MPIScatterGrid(mobile_massindex);
			this->MPIScatterGrid(depositionindex);
			
			// since the mass and hence the heat transport is quite turbulent in a snow slide it doesn't make sense to exactly 
			// calculate the mixing temperature of the remaining snow pack after a slide has been stride the actual cell.
			// However, the accumulation or source area of a slide does lose only the upper layer(s), so that can be done quite deterministic. 
			// The deposition will be done in a way that the deposited snow is simply added on top of the old snow.
			// later, we can consider a mixing temperature of the snow in motion to set the deposition temperature as well. 
			// For now, the temperature can be set to the value of the uppermost snow layer before the slide stroke the cell.
			if (bLayeredSnowActive)
			{ // we are now loacl again. Foreach cell a comparison between the bulk snow storage and the refSWES sum (in layer 0) indicates a 
			  // change during the slide module. If that happened, the refSWES must be changed now, too
				for (int i = mpistartrow; i <= mpiendrow; i++) // rows of the grid
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) // columns of the grid
					{
						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						{
							double diff = SSNOW[i][j] - refSWES[i][j][0];
							if (fabs(diff) > 0.0005)
							{
								int levelsSnow = LayeredSnowProperties[1].getLevels();
								int actLevelSnow = (int)refUpperIndexSnow[i][j];
								double snowTemp = (T[i][j]>0)?0:T[i][j]; // if there was no snow, we take the air temperature as first approximation
								if (bHeatTransferActive)
								{
									if (actLevelSnow<=levelsSnow) // if there was snow, we take the average from the existing snow cover
									{
										snowTemp = refTemperSnow[i][j][int((levelsSnow + actLevelSnow)/2.0)];
									}
								}
								UpdateSWEStack(i, j, 1.0, diff, 0, _onGlacierFalse_, _fromRefreezingFalse_,snowTemp); 
							}
						}
					}
				}
			}
		}
	}

	// statistics must be done outside the parallel section because common variables will be accessed by wertzufuegen
	// as well as by tracer_module[z].stat->zone[].sx and other operators
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				//      if (zn[i][j]!=nodata)
			{
				icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

				// all of the following statistics do not depend on glaciers, so the weight will be 1.0 by default (no icevalue-weighing is done for fractionally glacerized cells)
				snowstor->wertzufuegen(i, j, zarray[(int)zn[i][j]], SSNOW[i][j] + SLIQ[i][j]);
				qschnee->wertzufuegen(i, j, zarray[(int)zn[i][j]], QSNOW[i][j]);
				qsnow_out->wertzufuegen(i, j, zarray[(int)zn[i][j]], QSNOW_OUT[i][j]);
				snowdays->wertzufuegen(i, j, zarray[(int)zn[i][j]], SNOWDAYS[i][j]);
				snow_age->wertzufuegen(i, j, zarray[(int)zn[i][j]], SNOW_AGE[i][j]);
				p_snow->wertzufuegen(i, j, zarray[(int)zn[i][j]], P_SNOW[i][j]);
				p_rain->wertzufuegen(i, j, zarray[(int)zn[i][j]], P_RAIN[i][j]);
				albedo->wertzufuegen(i, j, zarray[(int)zn[i][j]], Albedo[i][j]);

				// new snowmodel extensions (Warscher 2013)
				if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
				{
					if (SNOWTEMP[i][j] > nodata + 1) snowtemp->wertzufuegen(i, j, zarray[(int)zn[i][j]], SNOWTEMP[i][j]);
					if (SURFTEMP[i][j] > nodata + 1) surftemp->wertzufuegen(i, j, zarray[(int)zn[i][j]], SURFTEMP[i][j]);
				}

				if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
				{
					input_mass->wertzufuegen(i, j, zarray[(int)zn[i][j]], INPUT_MASS[i][j]);
					mobile_mass->wertzufuegen(i, j, zarray[(int)zn[i][j]], MOBILE_MASS[i][j]);
					deposition->wertzufuegen(i, j, zarray[(int)zn[i][j]], DEPOSITION[i][j]);
				}

				if (ICE) // it is possible thatthere is no glaciermodel active!
				{
					snowstor2->wertzufuegen(i, j, zarray[(int)zn[i][j]], SSNOW[i][j] + SLIQ[i][j], 1.0 - icevalue);
					snowstor2->wertzufuegen(i, j, zarray[(int)zn[i][j]], 0.0, icevalue);
					if (icevalue>0.0)
					{
						double snowOnFirnWeight = (FIRN1[i][j]>0) ? 1.0 : 0.0;
						// for the glacierized part -- only the glacerized fraction is important
						//						qfirn->wertzufuegen(i,j,zarray_routing[(int)zn_routing[i][j]],QFIRN[i][j]+QSNOW[i][j]*icevalue); // icevalue fraction for QFIRN is already accounted for since melt rates were calculated using a icevalue factor

						qfirn->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], QFIRN[i][j] + snowOnFirnWeight*QSNOWGL[i][j], icevalue);
						qice->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], QICE[i][j], icevalue);
						qsnowonglac->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], (1.0 - snowOnFirnWeight)*QSNOWGL[i][j], icevalue);

						qfirn->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0, 1 - icevalue); // since QFIRN was already scaled, it is not required to add the unglacierized fraction woth 0 to the statistics
						qice->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0, 1 - icevalue); // dito
						qsnowonglac->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0, 1 - icevalue); // here, we need the QSNOW from the unglacierized fraction of the cell as 0-value

						if (glacmethode>10)
						{
							// mass balance for the glacerized part:
							// QSNOW was already reduced, so it has tobe re-reduced here
							// QFIRN and QICE where scaled to icevalue when calculating melt rates
							// So also P has to be scaled and GMB is consequently a scaled value too

							//13.12.2014							GMB[i][j]+=(P[i][j]-QSNOW[i][j]-QFIRN[i][j]-QICE[i][j]); 
							// BALANCE FOR GLACIER: 
							// P_SNOW accounts for wind driven accumulations/reductions already, 
							// P_RAIN must be added because it will be subtracted as runoff in QSNOW and/or QFIRN and/or QICE
							// DEOPSIION and MOBILE_MASS are masses that could be brought onto the glacier by avalanches, 
							// INPUT_MASS is the mass that could be brought away from the glacier by avalanches
							GMB[i][j] += (P_SNOW[i][j] + P_RAIN[i][j] - QSNOWGL[i][j] - QFIRN[i][j] - QICE[i][j]);
							if (DEPOSITION) GMB[i][j] += DEPOSITION[i][j] + MOBILE_MASS[i][j] - INPUT_MASS[i][j];
							// if evaporation model was on, then snow evaporation must be subtracted, too		
							if (EVRSGL) GMB[i][j] -= EVRSGL[i][j];

							firnstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], FIRN1[i][j], icevalue);
							gmbstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], GMB[i][j], icevalue);
							gmb2stat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], GMB[i][j], icevalue); // this is done here only, since unglacierized cells will not be part of this statistics

							firnstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0, 1.0 - icevalue);
							gmbstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0, 1.0 - icevalue);
						}

						// 12.12.2014 removed this section
						//						// Achtung; Zufluss zu firn, schnee oder eis speicher richtet sich danach, welche Unterlage unter dem Schnee ist:
						//						// Schmelze aus Schnee über Firn geht in den Firnspeicher (wird also eigentlich firnschmelze)
						//						// Schmelze aus Schnee ueber Eis bleibt aber im Schneespeicher.
						//						if (((glacmethode>10)&&(FIRN1[i][j]>0.0))||(ICE[i][j]==2))
						//						{
						//							// auf Firnflächen wird die Schneeschmelze als Firnschmelze betrachtet
						//							//	QFIRN[i][j] += QSNOW[i][j]*icevalue; // take into account the glacierized fraction only
						//							QFIRN[i][j] += QSNOW[i][j]; // take into account the glacierized fraction only
						//							// for remaining snow storage balance, do take into account the unglacierized fraction only
						//							// --> no, it's not important, since it counts for statistics only. subsequent modules will handle the glacierized fraction by their own
						//							// QSNOW[i][j] *= (1.0-icevalue); // und dafuer gibt es keine Schneeschmelze an diesen Flaechen -> wichtig fuer
						//							// Parallelspeicherkonzept, Firn ist am langsamsten!! und dominiert daher bei einer Schnee-auf-Firn-Situation
						//							// QSNOW_OUT[i][j] *= (1.0-icevalue);
						//						}

					}
					else
					{
						// totally unglacerized cells will have a 0-input for the entire cell, so this statistics will be valid for the entire basin (e.g. mm flow from firn for the entire basin)
						qfirn->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0);
						qice->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0);
						qsnowonglac->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0);
						if (glacmethode>10)
						{
							firnstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0);
							gmbstat->wertzufuegen(i, j, zarray_routing[(int)zn_routing[i][j]], 0.0);
						}
					}
				}

				//		icevalue=ICE?ICE[i][j]:nodata;
				//      if (ICE)
				//      {
				//        if ((ICE[i][j]>=1)&&(ICE[i][j]<=2)) /* 1 for ice, 2 for firn 
				//        {
				//          qfirn->wertzufuegen(i,j,zarray_routing[(int)zn_routing[i][j]],QFIRN[i][j]);
				//          qice->wertzufuegen(i,j,zarray_routing[(int)zn_routing[i][j]],QICE[i][j]);
				//	      }
				// 	    }

				if (SNOWFREEDAYS) // only if SNOWFREEDAYS exists
				{
					thawdepth->wertzufuegen(i, j, zarray[(int)zn[i][j]], THAWDEPTH[i][j]);
				}

				//        if ((!ICE)||((ICE)&&(icevalue==nodata)))
				if (icevalue == 0)
				{
					for (int z = 1; z <= tracer_count; z++)
					{
						if (tracer_info[z].evaporating != 0)
						{
							if (RelCellArea != NULL)
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx += tracer_module[z].concstor[i][j] * RelCellArea[i][j] * (SSNOW[i][j] + SLIQ[i][j]);
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n += (SSNOW[i][j] + SLIQ[i][j]);
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum += RelCellArea[i][j] * (SSNOW[i][j] + SLIQ[i][j]);
							}
							else
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx += tracer_module[z].concstor[i][j] * (SSNOW[i][j] + SLIQ[i][j]);
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n += (SSNOW[i][j] + SLIQ[i][j]);
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum += (SSNOW[i][j] + SLIQ[i][j]);
							}
						}
						else
						{
							if (RelCellArea != NULL)
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx += RelCellArea[i][j] * tracer_module[z].concstor[i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum += RelCellArea[i][j];
							}
							else
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx += tracer_module[z].concstor[i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum += 1.0;
							}
						}
					}
				}
			}
		}
	}

	snowstor->calc_and_done(jahr, monat, tag, stunde, minute);
	qschnee->calc_and_done(jahr, monat, tag, stunde, minute);
	qsnow_out->calc_and_done(jahr, monat, tag, stunde, minute);
	snowdays->calc_and_done(jahr, monat, tag, stunde, minute);
	snow_age->calc_and_done(jahr, monat, tag, stunde, minute);
	p_snow->calc_and_done(jahr, monat, tag, stunde,minute);
	p_rain->calc_and_done(jahr, monat, tag, stunde,minute);
	albedo->calc_and_done(jahr, monat, tag, stunde,minute);

	// new snowmodel extensions (Warscher 2013)
	if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
	{
		snowtemp->calc_and_done(jahr, monat, tag, stunde, minute);
		surftemp->calc_and_done(jahr, monat, tag, stunde, minute);
	}

	if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
	{
		input_mass->calc_and_done(jahr, monat, tag, stunde, minute);
		mobile_mass->calc_and_done(jahr, monat, tag, stunde, minute);
		deposition->calc_and_done(jahr, monat, tag, stunde, minute);
	}

	if (ICE)
	{
		snowstor2->calc_and_done(jahr, monat, tag, stunde, minute);
		qfirn->calc_and_done(jahr, monat, tag, stunde, minute);
		qice->calc_and_done(jahr, monat, tag, stunde, minute);
		qsnowonglac->calc_and_done(jahr, monat, tag, stunde, minute);
		if (glacmethode>10)
		{
			firnstat->calc_and_done(jahr, monat, tag, stunde, minute);
			gmbstat->calc_and_done(jahr, monat, tag, stunde, minute);
			gmb2stat->calc_and_done(jahr, monat, tag, stunde, minute);
		}
	}

	if (SNOWFREEDAYS) // only if SNOWFREEDAYS exists
	{
		thawdepth->calc_and_done(jahr, monat, tag, stunde, minute);
	}

	if (mpi_rank == 0)
	{
		if (ICE)
		{
			// ToDo 2013-05-21:
			// Gather all records of zoneQICE, zoneQFIRN and zoneQSnowOnGlac
			// -> should be done already by the calc_and_done methods of qice, qfirn and qsnowonglac earlier, so it is not necessary to do it here twice
			for (int i = 1; i <= zonmax_routing; i++)
			{
				//(*WSstream) << "DEBUG: glacierlist["<<i<<"].isGlacierized: "<<glacierlist[i].isGlacierized<<", glacierlist["<<i<<"].CellCntBare: "<<glacierlist[i].CellCntBare<<", glacierlist[i].CellCntGlacerized: "<<glacierlist[i].CellCntGlacerized<<"[end-entry]";
				if (glacierlist)
				{
					if (glacierlist[i].isGlacierized)
					{
						double glacierizedPart = glacierlist[i].CellCntGlacerized / (glacierlist[i].CellCntGlacerized + glacierlist[i].CellCntBare);

						Qice[i] = zoneQICE[i].sx / zoneQICE[i].RelCellSum; // -> Gesamtschmelze vom EIS fuer das Teilgebiet i (bezogen auf Gesamtfläche des Teilgebietes, nicht nur auf die Gletscher!)
						Qfirn[i] = zoneQFIRN[i].sx / zoneQFIRN[i].RelCellSum; // -> Gesamtschmelze vom FIRN fuer das Teilgebiet i (bezogen auf Gesamtfläche des Teilgebietes, nicht nur auf die Gletscher!)
						//			  Qsnow[i] = zoneQSNOW[i].sx/zoneQSNOW[i].RelCellSum; // -> Gesamtschmelze vom SNOW fuer das Teilgebiet i (bezogen auf Gesamtfläche des Teilgebietes, nicht nur auf die Gletscher!)/
						// now with a separate statistics for snow melt on glaciers
						Qsnow[i] = zoneQSnowOnGlac[i].sx / zoneQSnowOnGlac[i].RelCellSum; // -> Gesamtschmelze vom SNOW fuer das Teilgebiet i (bezogen auf Gesamtfläche des Teilgebietes, nicht nur auf die Gletscher!)/

						double step = InV / 60.0;
						// since a disappeared glacier would result in division by zero, we set the components to 0 in that cases
						Qice_els[i] = (glacierizedPart>0.0) ? (Qice_els[i] * exp(-step / kice_els[i]) + (1.0 - exp(-step / kice_els[i])) * Qice[i] / glacierizedPart) : 0.0;
						Qfirn_els[i] = (glacierizedPart>0.0) ? (Qfirn_els[i] * exp(-step / kfirn_els[i]) + (1.0 - exp(-step / kfirn_els[i]))*Qfirn[i] / glacierizedPart) : 0.0;
						Qsnow_els[i] = (glacierizedPart>0.0) ? (Qsnow_els[i] * exp(-step / ksnow_els[i]) + (1.0 - exp(-step / ksnow_els[i]))*Qsnow[i] / glacierizedPart) : 0.0;

						zoneQGLAC[i].sx = (Qice_els[i] + Qfirn_els[i] + Qsnow_els[i])*glacierizedPart*zoneQICE[i].n;
						zoneQGLAC[i].n = zoneQICE[i].n;
						zoneQGLAC[i].RelCellSum = zoneQICE[i].RelCellSum;

						// if no unsatzonmodel runs, then Qges is unavailable
						if (Qges) Qges[0][i] = zoneQGLAC[i].sx / zoneQGLAC[i].RelCellSum;
					}
					else
					{
						zoneQGLAC[i].sx = 0;
						zoneQGLAC[i].n = zoneQICE[i].n;
						zoneQGLAC[i].RelCellSum = zoneQICE[i].RelCellSum;
					}
				}
			}
		}
	}

	if (ICE)
	{
		qglac->calc_and_done(jahr, monat, tag, stunde, minute);
	}

	// balancing the glaciers area for the dynamic glacier model and switching status of the unmelted snow to firn for glacierized areas
	if (ICE && glacmethode>10)
	{
		if ((monat == GMBEvalMonth) && (tag == GMBEvalDay) && (((60 * se + int(me + 60 * hour_offset)) == 1440)))
		{
			// now we will at first shift down the firn layers, change snow to firn and set snow to 0 and save the old glacierized fraction
			for (int i = mpistartrow; i <= mpiendrow; i++) // rows of the grid
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) // columns of the grid
				{
					if ((zn[i][j] != nodata) && (ICE[i][j]>0.0) && (MPIREG[i][j] == mpi_rank)) // for glacierized cells only
					{
						// -------------------------------- start new code fragment 2016-10-13 ----------------------------------------
						// 2016-10-13 using maxSnowPackToCreateGlacier to decide wether part of the snow pack (i.e. the exceeding part above that threshold)
						// should be moved onto the glacier to make sure that no mounts of snow can grow beside tha glaciers
						if (ICE[i][j] < 1.0)
						{
							if (SSTO[i][j] > maxSnowPackToCreateGlacier)
							{
								if (SSTOGL)
								{
									if (ICE[i][j] < relSizeNewGlacier)
									{
										ICE[i][j] = relSizeNewGlacier;
										ICEOLD[i][j] = relSizeNewGlacier;
									}

									// total SWE to be transferred from the non-glacierized part to the glacierized part of the cell
									double transferamount = (SSTO[i][j] - maxSnowPackToCreateGlacier);
									// liquid fraction is weighed on the liquid fraction of the total SWE 
									double transferamountsliq = transferamount * SLIQ[i][j] / SSTO[i][j];
									// solid fraction as well 
									double transferamountsnow = transferamount * SSNOW[i][j] / SSTO[i][j];
									// adding the changes to the glacierized part, weighed by the ratio of the respective relative fractions to get the correct scaling factor for the glacierized fraction, 
									// e.g. if the glacierized fraction is 0.2 and the unglacierized is 0.8, the amount from the unglacierized fraction will be 0.8/0.2 = 4.0 times higher on the glacierized fraction
									SSTOGL[i][j] += transferamount*((1.0 - ICE[i][j]) / ICE[i][j]);
									SSNOWGL[i][j] += transferamountsnow*((1.0 - ICE[i][j]) / ICE[i][j]);
									SLIQGL[i][j] += transferamountsliq*((1.0 - ICE[i][j]) / ICE[i][j]);
									// subtracting the transferamount without weighing, since the concentration takes place on the glacierized fraction
									// only and the three storages here do not "know" about the glacier a all, 
									// i.e. it is as taking away the exceeding SWE in mm from the entire cell as an even layer
									SSTO[i][j] -= transferamount;
									SSNOW[i][j] -= transferamountsnow;
									SLIQ[i][j] -= transferamountsliq;
									GMB[i][j] += transferamount;
									if (bLayeredSnowActive)
									{
										int levelsSnow = LayeredSnowProperties[1].getLevels();
										int actLevelSnow = (int)refUpperIndexSnow[i][j];
										double snowTemp = 0;
										if (bHeatTransferActive)
										{
											refTemperSnow[i][j][int((levelsSnow + actLevelSnow)/2.0)];
										}
										UpdateSWEStack(i, j, 1.0, -transferamount, 0, _onGlacierFalse_, _fromRefreezingFalse_,snowTemp);
									}
								}
							}
						}
						else
						{ // ICE[i][j]==1.0, i.e. no unglacierized fraction in this cell, so all snow storages can be set to 0 in order to have grids without extreme snow accumulations
							SSTO[i][j] = 0;
							SLIQ[i][j] = 0;
							SSNOW[i][j] = 0;
						}
						// -------------------------------- end new code fragment 2016-10-13 ----------------------------------------

						// next line is done already right at the beginning of ::run in all timesteps to ensure unsatzonmodel finds out about the difference        
						//						ICEOLD[i][j]=ICE[i][j]; // ICEOLD will be used after resizing the glacier to redistribute the old balance to the remaining glacierized cells evenly
						GLACIDOLD[i][j] = GLACID[i][j]; // GLACIDOLD will also be used after resizing the glacier to redistribute the old balance to the remaining glacierized cells evenly
						FIRN1[i][j] = 0; // rest helper grid values
						for (int lv = firn_levels; lv>1; lv--)
						{
							FIRN[i][j][lv] = FIRN[i][j][lv - 1];
							FIRN1[i][j] += FIRN[i][j][lv];
						}
						FIRN[i][j][1] = SSTOGL ? SSTOGL[i][j] : SSTO[i][j];
						FIRN1[i][j] += SSTOGL ? SSTOGL[i][j] : SSTO[i][j];
						FIRN[i][j][0] = FIRN1[i][j];
						// all water from the snow storage will be firn now, even the liquid fraction, which is minimal anyway but important for the total balance
						if (SSTOGL) SSTOGL[i][j] = 0; else SSTO[i][j] = 0;
						if (SSNOWGL) SSNOWGL[i][j] = 0; else SSNOW[i][j] = 0;
						if (SLIQGL) SLIQGL[i][j] = 0; else SLIQ[i][j] = 0;
					}
					else if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank) && (ICE[i][j] <= 0.0) && (SSTO[i][j] > maxSnowPackToCreateGlacier)) // for unglacierized cells only where the snow pack doesn't melt for a long time
					{
						if (GLACID[i][j] < 0)
						{
							GLACID[i][j] = -1; // flag this cell to be added as a new glacier cell later in updateGlacierStatistics
						}
						else
						{
							GLACIDOLD[i][j] = GLACID[i][j];
						}
						// defining ICE[i][j]
						ICE[i][j] = relSizeNewGlacier;
						ICEOLD[i][j] = relSizeNewGlacier;
						if (SSTOGL) SSTOGL[i][j] = SSTO[i][j];
						if (SLIQGL) SLIQGL[i][j] = SLIQ[i][j];
						if (SSNOWGL) SSNOWGL[i][j] = SSNOW[i][j];

						// todo: moving snow from the unglacierized part of the cell to the glacierized part as above for the other cells that already were markes as glaciers
						if (SSTO[i][j] > maxSnowPackToCreateGlacier)
						{
							if (SSTOGL)
							{
								// total SWE to be transferred from the non-glacierized part to the glacierized part of the cell
								double transferamount = (SSTO[i][j] - maxSnowPackToCreateGlacier);
								// liquid fraction is weighed on the liquid fraction of the total SWE 
								double transferamountsliq = transferamount * SLIQ[i][j] / SSTO[i][j];
								// solid fraction as well 
								double transferamountsnow = transferamount * SSNOW[i][j] / SSTO[i][j];
								// adding the changes to the glacierized part, weighed by the ratio of the respective relative fractions to get the correct scaling factor for the glacierized fraction, 
								// e.g. if the glacierized fraction is 0.2 and the unglacierized is 0.8, the amount from the unglacierized fraction will be 0.8/0.2 = 4.0 times higher on the glacierized fraction
								SSTOGL[i][j] += transferamount*((1.0 - ICE[i][j]) / ICE[i][j]);
								SSNOWGL[i][j] += transferamountsnow*((1.0 - ICE[i][j]) / ICE[i][j]);
								SLIQGL[i][j] += transferamountsliq*((1.0 - ICE[i][j]) / ICE[i][j]);
								// subtracting the transferamount without weighing, since the concentration takes place on the glacierized fraction
								// only and the three storages here do not "know" about the glacier a all, 
								// i.e. it is as taking away the exceeding SWE in mm from the entire cell as an even layer
								SSTO[i][j] -= transferamount;
								SSNOW[i][j] -= transferamountsnow;
								SLIQ[i][j] -= transferamountsliq;
								GMB[i][j] += transferamount;
								if (bLayeredSnowActive)
								{
									int levelsSnow = LayeredSnowProperties[1].getLevels();
									int actLevelSnow = (int)refUpperIndexSnow[i][j];
									double snowTemp = 0;
									if (bHeatTransferActive) 
									{
										snowTemp = refTemperSnow[i][j][int((levelsSnow + actLevelSnow)/2.0)];
									}
									UpdateSWEStack(i, j, 1.0, -transferamount, 0, _onGlacierFalse_, _fromRefreezingFalse_,snowTemp);
								}
							}
						}

						FIRN1[i][j] = 0; // rest helper grid values
						for (int lv = firn_levels; lv>1; lv--)
						{
							FIRN[i][j][lv] = 0;
						}
						FIRN[i][j][1] = SSTOGL ? SSTOGL[i][j] : SSTO[i][j];
						FIRN1[i][j] += SSTOGL ? SSTOGL[i][j] : SSTO[i][j];
						FIRN[i][j][0] = FIRN1[i][j];
						// all water from the snow storage will be firn now, even the liquid fraction, which is minimal anyway but important for the total balance
						if (SSTOGL) SSTOGL[i][j] = 0; else SSTO[i][j] = 0;
						if (SSNOWGL) SSNOWGL[i][j] = 0; else SSNOW[i][j] = 0;
						if (SLIQGL) SLIQGL[i][j] = 0; else SLIQ[i][j] = 0;
						//GMBO[i][j] = 0;
						GMBO[i][j] = GMB[i][j] = FIRN1[i][j]; // new glacier cell: balance starts with the freshly added firn
					}
				}
			}

			// ToDo 2013-05-21:
			// Gather all Grids and other structures here which will be used in the following method
			// FIRN, GMB,
			MPIGatherGrid(iceoldindex);
			MPIGatherGrid(glacidoldindex);
			MPIGatherGrid(gmbindex);
			MPIGatherGrid(gmboindex);
			MPIGatherGrid(firnindex);
			MPIGatherStack(firnindex);
			MPIGatherGrid(getindex(_GlacierizedCells_));
			MPIGatherGrid(getindex(_GlacierCodes_));
			// now grow and erode the ice using the GMB and GMBO grids
			// this must be done in a single thread for the entire grid since lateral relations are important
			if (mpi_rank == 0)
			{
				GlacierRetreatAndAdvance();
				RedistributeOldBalance(-1);
			}
			MPIScatterGrid(getindex(_GlacierizedCells_));
			MPIScatterGrid(getindex(_GlacierCodes_));
			MPIScatterGrid(gmbindex);
			MPIScatterGrid(gmboindex);
			MPIScatterStack(firnindex);
			// iceoldgrid will not be changed and thus doen't need to be scattered

			// reset counters for weights
			(*WSstream) << "DEBUG: glaciers changed, updating glaicerlist now..." << "[end-entry]";
			if (glacierlist)
			{
				for (int i = 1; i <= zonmax_routing; i++)
				{
					if (glacierlist[i].isGlacierized)
					{
						glacierlist[i].CellCntBare = 0;
						glacierlist[i].CellCntGlacerized = 0;
					}
				}
				// in the end, the fields of glacierlist will contain the number of glacierized and unglacierized cells as floating point numbers
				for (int i = 1; i <= nrows; i++) // rows of the grid
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) // columns of the grid
					{
						if (zn[i][j] != nodata)
						{
							double icevalue = ICE[i][j]; // valid range of icevalue now between 0...2
							if (icevalue>1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
							else if (icevalue<0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
							glacierlist[zarray_routing[(int)zn_routing[i][j]]].CellCntBare += (1.0 - icevalue);
							glacierlist[zarray_routing[(int)zn_routing[i][j]]].CellCntGlacerized += icevalue;
						}
					}
				}
			}

#ifdef _MPI
			// bcastig glacierlist to all other ranks
			MPI_Barrier(MPI_COMM_WORLD);
			const int nitems = 3;
			int blocklengths[3] = { 1, 1, 1 };
			MPI_Datatype types[3] = { MPI_C_BOOL, MPI_DOUBLE, MPI_DOUBLE };
			MPI_Aint offsets[3] = { offsetof(glacierlisttype, isGlacierized),
				offsetof(glacierlisttype, CellCntGlacerized),
				offsetof(glacierlisttype, CellCntBare) };

			MPI_Datatype MPIglacierlisttype;
			MPI_Type_create_struct(nitems, blocklengths, offsets, types,&MPIglacierlisttype);
			MPI_Type_commit(&MPIglacierlisttype);

			if (mpi_rank > 0)
			{
				glacierlist = (glacierlisttype*)calloc(zonmax_routing + 1, sizeof(glacierlisttype));
				if (!glacierlist) speichererror(1);
			}
			MPI_Bcast(&(glacierlist[0]), zonmax_routing + 1, MPIglacierlisttype, 0,MPI_COMM_WORLD);

			MPI_Type_free(&MPIglacierlisttype);
#endif

			// ToDo 2013-05-21:
			// Scatter all Grids and other structures here which will be used in other MPI-Instances
			// Scatter also glacierlist

			// at the end, set the old balance grid to the actual balance grid value and set cells to 0 which are no longer part of the glacier
			for (int i = mpistartrow; i <= mpiendrow; i++) // rows of the grid
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) // columns of the grid
				{
					if ((zn_routing[i][j] != nodata) && (ICE[i][j]>0))
					{
						GMBO[i][j] = GMB[i][j]; // GMB is not changed, so we have a cumulative series
					}
					else
					{
						GMB[i][j] = 0;
						GMBO[i][j] = 0;
					}
				}
			}
		} // end of special balance day
	}

	// Schreiben der Ergebnisgrids mit den Daten des Ausflusses aus der Schneedecke
	// und den aktuellen Speicherinhalten 
	writegridtofile(sstoindex);
	writegridtofile(ssnowindex);
	writegridtofile(sliqindex);
	writegridtofile(snowdaysindex);
	writegridtofile(snow_ageindex);
	writegridtofile(p_snowindex);
	writegridtofile(p_rainindex);
	writegridtofile(albedoindex);
	writegridtofile(qsnowindex);
	writegridtofile(qsnow_outindex);

	// new snowmodel extensions (Warscher 2013)
	if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
	{
		writegridtofile(input_massindex);
		writegridtofile(mobile_massindex);
		writegridtofile(depositionindex);
	}

	if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
	{
		writegridtofile(snowtempindex);
		writegridtofile(surftempindex);
	}


	if (ICE)
	{
		writegridtofile(qfirnindex);
		writegridtofile(qiceindex);
		if (glacmethode>10)
		{
			writestacktofile(firnindex); // this is a stack! (WE for each layer of the firn grid)
			writegridtofile(gmbindex);  // glacier mass balance
			writegridtofile(gmboindex); // glacier mass balance for the last Voluma-Area-Modelling time (usually October 1st)
			writegridtofile(ssnowglindex); // solid snow storage for glacierized cell fractions
			writegridtofile(sliqglindex); // liquid snow storage for glacierized cell fractions
			writegridtofile(sstoglindex); // total snow storage for glacierized cell fractions
			writegridtofile(qsnowglindex); // total snow storage for glacierized cell fractions
			if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
			{
				writegridtofile(surftempglindex); // total snow storage for glacierized cell fractions
				writegridtofile(snowtempglindex); // total snow storage for glacierized cell fractions
			}
			writegridtofile(evrsglindex); // total snow storage for glacierized cell fractions
		}
	}

	for (int z = 1; z <= tracer_count; z++)
	{
		tracer_module[z].stat->calc_and_done(jahr, monat, tag, stunde, minute);
		writegridtofile(tracer_module[z].storindex);
	}
}


// 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
// C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
// Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. 
void schneemodellklasse::ClearMemory()
{	
	if (do_it > 0)
	{
		// freigeben der belegten Speicher 
		if (initialized!=0) WriteSlrContents(true);

		qschnee->ClearMemory();
		qsnow_out->ClearMemory();
		snowstor->ClearMemory();
		snowdays->ClearMemory();
		snow_age->ClearMemory();
		p_snow->ClearMemory();
		p_rain->ClearMemory();
		albedo->ClearMemory();

		if (C0) delete[] C0; C0 = NULL;
		if (C1) delete[] C1; C1 = NULL;
		if (C2) delete[] C2; C2 = NULL;

		// new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			snowtemp->ClearMemory();
			surftemp->ClearMemory();
		}

		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			input_mass->ClearMemory();
			mobile_mass->ClearMemory();
			deposition->ClearMemory();
		}

		if (qschnee) delete(qschnee);
		if (qsnow_out) delete(qsnow_out);
		if (snowstor) delete(snowstor);
		if (snowdays) delete(snowdays);
		if (snow_age) delete(snow_age);
		if (p_snow) delete(p_snow);
		if (p_rain) delete(p_rain);
		if (albedo) delete(albedo);
		// new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			if (snowtemp) delete(snowtemp);
			if (surftemp) delete(surftemp);
		}
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			if (input_mass) delete(input_mass);
			if (mobile_mass) delete(mobile_mass);
			if (deposition) delete(deposition);
		}

		if (zoneQSNOW) free(zoneQSNOW); zoneQSNOW=NULL;
		if (zoneSSNOW) free(zoneSSNOW); zoneSSNOW=NULL;
		if (zoneSSNOW2) free(zoneSSNOW2); zoneSSNOW2=NULL;
		if (zoneQSNOW_OUT) free(zoneQSNOW_OUT);zoneQSNOW_OUT=NULL;
		if (zoneSNOWDAYS) free(zoneSNOWDAYS);zoneSNOWDAYS=NULL;
		if (zoneSNOW_AGE) free(zoneSNOW_AGE);zoneSNOW_AGE=NULL;
		if (zoneP_SNOW) free(zoneP_SNOW);zoneP_SNOW=NULL;
		if (zoneP_RAIN) free(zoneP_RAIN);zoneP_RAIN=NULL;
		if (zoneALBEDO) free(zoneALBEDO);zoneALBEDO=NULL;
		// new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			if (zoneSNOWTEMP) free(zoneSNOWTEMP); zoneSNOWTEMP = NULL;
			if (zoneSURFTEMP) free(zoneSURFTEMP); zoneSURFTEMP = NULL;
		}
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			if (zoneINPUT_MASS) free(zoneINPUT_MASS); zoneINPUT_MASS = NULL;
			if (zoneMOBILE_MASS) free(zoneMOBILE_MASS);zoneMOBILE_MASS=NULL;
			if (zoneDEPOSITION) free(zoneDEPOSITION);zoneDEPOSITION=NULL;
		}

		free_matrix(gridinfo[ssnowindex].pointer, gridinfo[ssnowindex].nrl, gridinfo[ssnowindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[snowdaysindex].pointer, gridinfo[snowdaysindex].nrl, gridinfo[snowdaysindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[snow_ageindex].pointer, gridinfo[snow_ageindex].nrl, gridinfo[snow_ageindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[p_snowindex].pointer, gridinfo[p_snowindex].nrl, gridinfo[p_snowindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[p_rainindex].pointer, gridinfo[p_rainindex].nrl, gridinfo[p_rainindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[albedoindex].pointer, gridinfo[albedoindex].nrl, gridinfo[albedoindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[sliqindex].pointer, gridinfo[sliqindex].nrl, gridinfo[sliqindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[sstoindex].pointer, gridinfo[sstoindex].nrl, gridinfo[sstoindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[qsnowindex].pointer, gridinfo[qsnowindex].nrl, gridinfo[qsnowindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[qsnow_outindex].pointer, gridinfo[qsnow_outindex].nrl, gridinfo[qsnow_outindex].nrh, 1, ncols, nrows);

		// new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			free_matrix(gridinfo[snowtempindex].pointer, gridinfo[snowtempindex].nrl, gridinfo[snowtempindex].nrh, 1, ncols, nrows);
			free_matrix(gridinfo[surftempindex].pointer, gridinfo[surftempindex].nrl, gridinfo[surftempindex].nrh, 1, ncols, nrows);
		}
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			free_matrix(gridinfo[input_massindex].pointer, gridinfo[input_massindex].nrl, gridinfo[input_massindex].nrh, 1, ncols, nrows);
			free_matrix(gridinfo[mobile_massindex].pointer, gridinfo[mobile_massindex].nrl, gridinfo[mobile_massindex].nrh, 1, ncols, nrows);
			free_matrix(gridinfo[depositionindex].pointer, gridinfo[depositionindex].nrl, gridinfo[depositionindex].nrh, 1, ncols, nrows);
			if (mpi_rank == 0)
			{
				free_matrix(gridinfo[depositionmaxindex].pointer, gridinfo[depositionmaxindex].nrl, gridinfo[depositionmaxindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[elevorderrowindex].pointer, gridinfo[elevorderrowindex].nrl, gridinfo[elevorderrowindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[elevordercolindex].pointer, gridinfo[elevordercolindex].nrl, gridinfo[elevordercolindex].nrh, 1, ncols, nrows);
			}
		}

		delete_grid(ssnowindex);
		delete_grid(sliqindex);
		delete_grid(sstoindex);
		delete_grid(snowdaysindex);
		delete_grid(snow_ageindex);
		delete_grid(p_snowindex);
		delete_grid(p_rainindex);
		delete_grid(albedoindex);
		delete_grid(qsnowindex);
		delete_grid(qsnow_outindex);
		// new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			delete_grid(snowtempindex);
			delete_grid(surftempindex);
		}

		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			delete_grid(input_massindex);
			delete_grid(mobile_massindex);
			delete_grid(depositionindex);
			if (mpi_rank == 0)
			{
				delete_grid(depositionmaxindex);
				delete_grid(elevorderrowindex);
				delete_grid(elevordercolindex);
			}
		}

		if (SNOWDAYSgridname   ) free(SNOWDAYSgridname   ); SNOWDAYSgridname=NULL;
		if (SNOW_AGEgridname   ) free(SNOW_AGEgridname   ); SNOW_AGEgridname=NULL;
		if (P_SNOWgridname     ) free(P_SNOWgridname     ); P_SNOWgridname=NULL;
		if (P_RAINgridname     ) free(P_RAINgridname     ); P_RAINgridname=NULL;
		if (qschneegridname    ) free(qschneegridname    ); qschneegridname=NULL;
		if (qsnow_outgridname  ) free(qsnow_outgridname  ); qsnow_outgridname=NULL;
		if (SSNOWgridname      ) free(SSNOWgridname      ); SSNOWgridname=NULL;
		if (SLIQgridname       ) free(SLIQgridname       ); SLIQgridname=NULL;
		if (SSTOgridname       ) free(SSTOgridname       ); SSTOgridname=NULL;
		if (ICE)
		{
			if (qicegridname) free(qicegridname); qicegridname=NULL;
			if (qfirngridname) free(qfirngridname); qfirngridname=NULL;
		}

		// new snow model extensions (Warscher 2013)
//		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		if (methode > 4)
		{
			if (SNOWTEMPgridname) free(SNOWTEMPgridname); SNOWTEMPgridname = NULL;
			if (SURFTEMPgridname) free(SURFTEMPgridname); SURFTEMPgridname = NULL;
			if (EBALgridname) free(EBALgridname); EBALgridname = NULL;
//		}
//		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
//		{
			if (INPUT_MASSgridname) free(INPUT_MASSgridname); INPUT_MASSgridname = NULL;
			if (MOBILE_MASSgridname   ) free(MOBILE_MASSgridname   ); MOBILE_MASSgridname=NULL;
			if (DEPOSITIONgridname   ) free(DEPOSITIONgridname   ); DEPOSITIONgridname=NULL;
		}

		if (SNOWFREEDAYS)
		{
			free_matrix(gridinfo[snowfreedaysindex].pointer, gridinfo[snowfreedaysindex].nrl, gridinfo[snowfreedaysindex].nrh, 1, ncols, nrows);
			delete_grid(snowfreedaysindex);
			free_matrix(gridinfo[snowcoverdaysindex].pointer, gridinfo[snowcoverdaysindex].nrl, gridinfo[snowcoverdaysindex].nrh, 1, ncols, nrows);
			delete_grid(snowcoverdaysindex);

			thawdepth->ClearMemory();
			if (thawdepth) delete(thawdepth);
			
			if (zoneTHAWDEPTH)
			{
				free(zoneTHAWDEPTH);
				zoneTHAWDEPTH=NULL;
			}
			if (thawdepthstatname)
			{
				free(thawdepthstatname);
				thawdepthstatname=NULL;
			}
		}

		// freigeben von Speicher 
		if (ICE)
		{
			if (QICE)
			{
				free_matrix(gridinfo[qiceindex].pointer, gridinfo[qiceindex].nrl, gridinfo[qiceindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[qfirnindex].pointer, gridinfo[qfirnindex].nrl, gridinfo[qfirnindex].nrh, 1, ncols, nrows);
				delete_grid(qiceindex);
				delete_grid(qfirnindex);
			}

			if (snowstor2) {
				snowstor2->ClearMemory();
				delete(snowstor2);
				snowstor2=NULL;
			}
			if (qglac)  {
                qglac->ClearMemory();
                delete(qglac);
                qglac=NULL;
            }
			if (qice)   {
                qice->ClearMemory();
                delete(qice);
                qice=NULL;
            }
			if (qfirn)  {
                qfirn->ClearMemory();
                delete(qfirn);
                qfirn=NULL;
            }
			if (qsnowonglac) {
                qsnowonglac->ClearMemory();
                delete(qsnowonglac);
                qsnowonglac=NULL;
            }

			if (zoneQFIRN) free(zoneQFIRN);     zoneQFIRN=NULL;
			if (zoneQICE) free(zoneQICE);       zoneQICE=NULL;
			if (zoneQGLAC) free(zoneQGLAC);     zoneQGLAC=NULL;
			if (zoneQSnowOnGlac) free(zoneQSnowOnGlac); zoneQSnowOnGlac=NULL;

			if (Qice_els ) free(Qice_els);  Qice_els=NULL;
			if (Qfirn_els) free(Qfirn_els); Qfirn_els=NULL;
			if (Qsnow_els) free(Qsnow_els); Qsnow_els=NULL;

			if (Codes)	    free(Codes);      Codes     = NULL;
			if (DDFice)	    free(DDFice);     DDFice    = NULL;
			if (DDFfirn)	free(DDFfirn);	  DDFfirn   = NULL;
			if (DDFsnow)	free(DDFsnow);	  DDFsnow   = NULL;
			if (MF)	        free(MF);		  MF	    = NULL;
			if (Asnow_a)	free(Asnow_a);    Asnow_a   = NULL;
			if (Aice_a) 	free(Aice_a);	  Aice_a    = NULL;
			if (Asnow_b) 	free(Asnow_b);	  Asnow_b   = NULL;
			if (Aice_b)		free(Aice_b);	  Aice_b    = NULL;
			if (Asnow_akt)	free(Asnow_akt);  Asnow_akt = NULL;
			if (Aice_akt)	free(Aice_akt);	  Aice_akt  = NULL;
			if (Aice_min)	free(Aice_min);	  Aice_min  = NULL;
			if (Aice_max)	free(Aice_max);	  Aice_max  = NULL;
			if (Asnow_min)	free(Asnow_min);  Asnow_min = NULL;
			if (Asnow_max)	free(Asnow_max);  Asnow_max = NULL;
			if (kice_els)	free(kice_els);	  kice_els  = NULL;
			if (kfirn_els)	free(kfirn_els);  kfirn_els = NULL;
			if (ksnow_els)	free(ksnow_els);  ksnow_els = NULL;
			if (firn_base)	free(firn_base);  firn_base = NULL;
			if (deltafirnper100m)	free(deltafirnper100m);  deltafirnper100m = NULL;

			if (Qice) free(Qice); Qice=NULL;
			if (Qfirn) free(Qfirn); Qfirn=NULL;
			if (Qsnow) free(Qsnow); Qsnow=NULL;

			if (qicestatname) free(qicestatname); qicestatname=NULL;
			if (qfirnstatname) free(qfirnstatname); qfirnstatname=NULL;
			if (qglacstatname) free(qglacstatname); qglacstatname=NULL;
			if (qsnowonglacstatname) free(qsnowonglacstatname); qsnowonglacstatname=NULL;
			if (snowstorstatname2  ) free(snowstorstatname2  ); snowstorstatname2  =NULL;
			if (SSNOWGLgridname)	 free(SSNOWGLgridname	); SSNOWGLgridname	   = NULL;
			if (SLIQGLgridname)		 free(SLIQGLgridname	); SLIQGLgridname	   = NULL;
			if (SSTOGLgridname)		 free(SSTOGLgridname	); SSTOGLgridname	   = NULL;
			if (QSNOWGLgridname)	 free(QSNOWGLgridname	); QSNOWGLgridname	   = NULL;
			if (EBALGLgridname)		 free(EBALGLgridname    ); EBALGLgridname	   = NULL;
			if (SURFTEMPGLgridname)	 free(SURFTEMPGLgridname); SURFTEMPGLgridname  = NULL;
			if (SNOWTEMPGLgridname)	 free(SNOWTEMPGLgridname); SNOWTEMPGLgridname  = NULL;
			if (EVRSGLgridname)		 free(EVRSGLgridname	); EVRSGLgridname	   = NULL;

			if (glacmethode>10)
			{
				if (zoneFIRN) free(zoneFIRN);    zoneFIRN=NULL;
				if (zoneGMB)  free(zoneGMB);     zoneGMB =NULL;
				if (zoneGMB2) free(zoneGMB2);    zoneGMB2=NULL;

				if (firnstat)
                                {
                                    firnstat->ClearMemory();
                                    delete firnstat;
                                    firnstat=NULL;
                                }
				if (gmbstat)
                                {
                                    gmbstat->ClearMemory();
                                    delete gmbstat;
                                    gmbstat=NULL;
                                }
				if (gmb2stat)
                                {
                                    gmb2stat->ClearMemory();
                                    delete gmb2stat;
                                    gmb2stat=NULL;
                                }

				if (firnstatname) free(firnstatname); firnstatname=NULL;
				if (gmbstatname) free(gmbstatname); gmbstatname=NULL;
				if (gmb2statname) free(gmb2statname); gmb2statname=NULL;

				if (gmbgridname) free(gmbgridname); gmbgridname=NULL;
				if (gmbogridname) free(gmbogridname); gmbogridname=NULL;
				if (firnstackname) free(firnstackname); firnstackname=NULL;

				free_matrix(gridinfo[gmbindex].pointer, gridinfo[gmbindex].nrl, gridinfo[gmbindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[gmboindex].pointer, gridinfo[gmboindex].nrl, gridinfo[gmboindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[iceoldindex].pointer, gridinfo[iceoldindex].nrl, gridinfo[iceoldindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[glacidoldindex].pointer, gridinfo[glacidoldindex].nrl, gridinfo[glacidoldindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[sstoglindex].pointer, gridinfo[sstoglindex].nrl, gridinfo[sstoglindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[ssnowglindex].pointer, gridinfo[ssnowglindex].nrl, gridinfo[ssnowglindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[sliqglindex].pointer, gridinfo[sliqglindex].nrl, gridinfo[sliqglindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[qsnowglindex].pointer, gridinfo[qsnowglindex].nrl, gridinfo[qsnowglindex].nrh, 1, ncols, nrows);
				free_matrix(gridinfo[evrsglindex].pointer, gridinfo[evrsglindex].nrl, gridinfo[evrsglindex].nrh, 1, ncols, nrows);

				if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
				{
					free_matrix(gridinfo[surftempglindex].pointer, gridinfo[surftempglindex].nrl, gridinfo[surftempglindex].nrh, 1, ncols, nrows);
					free_matrix(gridinfo[snowtempglindex].pointer, gridinfo[snowtempglindex].nrl, gridinfo[snowtempglindex].nrh, 1, ncols, nrows);
				}

				
				delete_grid(gmbindex);
				delete_grid(gmboindex);
				delete_grid(iceoldindex);
				delete_grid(glacidoldindex);
				delete_grid(sstoglindex);
				delete_grid(ssnowglindex);
				delete_grid(sliqglindex);
				delete_grid(qsnowglindex);
				delete_grid(evrsglindex);
				if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
				{
					delete_grid(surftempglindex);
					delete_grid(snowtempglindex);
				}

				free_drei_d_matrix(FIRN, 1, nrows, 1, ncols, firn_levels, nrows);
				FIRN = NULL;
	
				free_matrix(gridinfo[firnindex].pointer, gridinfo[firnindex].nrl, gridinfo[firnindex].nrh, 1, ncols, nrows);
				FIRN1 = NULL;

                delete_grid(firnindex);

                for (int i=0;i<glaciers.size();i++)
                {
	                delete glaciers[i]; // deletign dynamically allocated instance of the glacier and calling its destrcutor here
                }
                glaciers.clear();

			}

			if (glacierlist) free(glacierlist); glacierlist=NULL;
		}

		if (qschneestatname    ) free(qschneestatname    ); qschneestatname    =NULL;
		if (snowstorstatname   ) free(snowstorstatname   ); snowstorstatname   =NULL;
		if (qsnow_outstatname  ) free(qsnow_outstatname  ); qsnow_outstatname  =NULL;
		if (snowdaysstatname   ) free(snowdaysstatname   ); snowdaysstatname   =NULL;
		if (snow_agestatname   ) free(snow_agestatname   ); snow_agestatname   =NULL;
		if (p_snowstatname     ) free(p_snowstatname     ); p_snowstatname     =NULL;
		if (p_rainstatname     ) free(p_rainstatname     ); p_rainstatname     =NULL;
		if (albedostatname     ) free(albedostatname     ); albedostatname     =NULL;
		// new snowmodel extensions (Warscher 2013)
//		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		if (methode>4)
		{
			if (snowtempstatname) free(snowtempstatname); snowtempstatname = NULL;
			if (surftempstatname) free(surftempstatname); surftempstatname = NULL;
//		}
//		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
//		{
			if (input_massstatname) free(input_massstatname); input_massstatname = NULL;
			if (mobile_massstatname   ) free(mobile_massstatname   ); mobile_massstatname   =NULL;
			if (depositionstatname   ) free(depositionstatname   ); depositionstatname   =NULL;
		}

		if (LUSE) delete[] LUSE; LUSE=NULL;

		if (zone_routing)
		{
			free(zone_routing);
			zone_routing = NULL;
		}
		if (zarray_routing)
		{
			free(zarray_routing);
			zarray_routing = NULL;
		}

		if (storagefilename) free (storagefilename); storagefilename = NULL;

	}
	modelobject::ClearMemory();

}

void schneemodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{	int z;
	char *tempstr1,*tempstr2,*tempstr3;
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (retvalue!=0) throw retvalue;
		if (tracer_count<=0) return;
		else
		{
			tempstr1=newstring(100,sizeof(char));
			tempstr2=newstring(100,sizeof(char));
			tempstr3=newstring(100,sizeof(char));

			for (z=1;z<=tracer_count;z++)
			{
				sprintf(tempstr3,"!!_no_entry_!!");
				// if the tracer is evaporating, the input concentration is always the same as in the last storage,
				// but if it is not evaporating (salts), the input-concentration may not rise above the maximum
				// concentration, therefore another grid containing the concentration in the flux is allocated and
				// consequently searched for here in such cases 
				sprintf(tempstr2,"concflux_tracer_%1d_%s",z,"input");
				sprintf(tempstr1,"concflux_tracer_%1d_%s",z,_irrigation_parameters_); // if irrigation is considered, tracer concentrations are get from irrigation submodel 
				tracer_module[z].cinput=gridinfo[getindex_first_of_3(tempstr1,tempstr2,tempstr3)].pointer;
			}
			if (tempstr1) free(tempstr1); tempstr1=NULL;
			if (tempstr2) free(tempstr2); tempstr2=NULL;
			if (tempstr3) free(tempstr3); tempstr3=NULL;
			return;
		}
	}
}

void schneemodellklasse::substance_transport_delete()
{	if (do_it==0) return;
	else
	teilmodellklasse::substance_transport_delete();
}

void schneemodellklasse ::flushbuffers()
{	int z;
	if (do_it==0) return;
	else
	{
		if (qschnee) qschnee->flushbuffer();
		if (qsnow_out) qsnow_out->flushbuffer();
		if (snowstor) snowstor->flushbuffer();
		if (snowdays) snowdays->flushbuffer();
		if (snow_age) snow_age->flushbuffer();
		if (p_snow) p_snow->flushbuffer();
		if (p_rain) p_rain->flushbuffer();
		if (albedo) albedo->flushbuffer();
	    // new snowmodel extensions (Warscher 2013)
		if ((methode == 5) || (methode == 6) || (methode == 7) || (methode == 11)) // Enbal
		{
			if (snowtemp) snowtemp->flushbuffer();
			if (surftemp) surftemp->flushbuffer();
		}
		if ((methode == 6) || (methode == 7) || (methode == 8) || (methode == 9) || (methode == 12) || (methode == 14)) // slides
		{
			if (input_mass) input_mass->flushbuffer();
			if (mobile_mass) mobile_mass->flushbuffer();
			if (deposition) deposition->flushbuffer();
		}

		if (ICE)
		{
			if (snowstor2) snowstor2->flushbuffer();
			if (qfirn) qfirn->flushbuffer();
			if (qice)  qice->flushbuffer();
			if (qglac) qglac->flushbuffer();
			if (qsnowonglac) qsnowonglac->flushbuffer();
			if (glacmethode>10)
			{
				firnstat->flushbuffer();
				gmbstat->flushbuffer();
				gmb2stat->flushbuffer();
			}
		}
		
		if (SNOWFREEDAYS)
		{
			if (thawdepth) thawdepth->flushbuffer();
		}

		for (z=1;z<=tracer_count;z++)
		tracer_module[z].stat->flushbuffer();
	}
}

// --------------------------------------
// --------------------------------------
// glacier related methods for snowmodel
// --------------------------------------
// --------------------------------------

void schneemodellklasse::GlacierRetreatAndAdvance()
{
	// after a certain time (usually one year) the balance changes are distributed across the area of the glacier and certain margin cells will loose or gain material
	// the statistics zones are already created, but they have to be estimated again after changing the glaciers extend. Only the old volume and Area will remain

	if (ICE && glacmethode>10)
	{

		// MPI: gather all grids needed here

		for (int glc=0;glc<glaciers.size();glc++)
		{
			glaciers[glc]->cnt=0;
			glaciers[glc]->Area=0;
			glaciers[glc]->Hmax=-10000;
			glaciers[glc]->Hmin=10000;
			// make sure the Volume is the one that was the result of the last balance 
			glaciers[glc]->setVolume(glaciers[glc]->getOldVolume());
		}

		UpdateGlacierStatisticsZones();

		for (int glc=0;glc<glaciers.size();glc++)
		{
			if (glaciers[glc]->Hmin==10000)
			{
				// no firn free ice cell found for this glacier
				glaciers[glc]->cnt=0;
				glaciers[glc]->Area=0;
				glaciers[glc]->Hmax=-10000;
				glaciers[glc]->Hmin=10000;
				UpdateGlacierStatisticsZones(glaciers[glc]->ID,true);
			}
			glaciers[glc]->createElevationBelts();
		}

		CountMarginsGroupedByElevationBelts();

		// make a statistic for the Balance per glacier (separated by GLACID-codes)
		// the result of this procedure will be the new Volume
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				if ((zn_routing[i][j]!=nodata)&&(ICE[i][j]>=0.0)&&(GLACID[i][j]>=0)) // for glacierized cells only
				{
					double icevalue = ICE[i][j];
					if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
					{
						int index=glacierindexmap.find((int)GLACID[i][j])->second;
						// dV must be given in 1'000'000 m^3, e.g. if cellsize is 100m, 
						// 1mm balance change is 10000m^3*0.001m = 10 m^3 = 0.00001 * 10^6 m^3
						// rho_ice is the density of ice (0.85 kg/l), resulting in a larger volume for ice than for water
						glaciers[index]->addBalanceChangeToVolume(((GMB[i][j]-GMBO[i][j])/rho_ice)*icevalue/1000.0*cellsize*cellsize/1000000.0);
					}
				}
			}
		}
		(*WSstream) << "DEBUG: Glacierbalance\n" <<"[end-entry]";
		for (int glc=0;glc<glaciers.size();glc++)
		{
			glaciers[glc]->calcArea(); // now the new Volume is calculated and we can get the new Area
			glaciers[glc]->calcTotalChangeInCells(); // also resets the internal value of remainingChange (member of glaciers) to 0
			double change=glaciers[glc]->getTotalChangeInCells();
			double changePerCell=glaciers[glc]->changePerCell(0,false);
			(*WSstream) << "DEBUG: " 
				<< "Glacier:\t"<<glc<< "\n"
				<< "Old Volume:\t"<<glaciers[glc]->getOldVolume() << "\n"
				<< "New Volume:\t"<<glaciers[glc]->getVolume() << "\n"
				<< "Old Area:\t"<<glaciers[glc]->getOldArea() << "\n"
				<< "New Area:\t"<<glaciers[glc]->getArea() << "\n"
		        <<"[end-entry]";

			// core of the algorithm: add or remove the changes in Area from the marginal cells$
			double remainingChange=change;
			int cnt=0;
			while (fabs(remainingChange)>0.00001)
			{
				long changeableCellCount = glaciers[glc]->getChangeableCellCount(remainingChange);
				if (changeableCellCount > 0)
				{
					if (cnt==0)
					{
						changePerCell=glaciers[glc]->changePerCell(0,false);
					}
					else
					{
						changePerCell=glaciers[glc]->changePerCell(remainingChange,true);
					}
					(*WSstream) << "DEBUG: " 
						<< "number of cells to add: \t" << remainingChange << "\n"
						<< "number of changeable cells: \t" << changeableCellCount  << "\n"
						<< "fraction of cells to add per iteration: \t" << changePerCell
						<<"[end-entry]";
	
					double oldRemainingChange=remainingChange; // store old value to chack after adding values to margin cells if the remainingChange value has changes at all
					remainingChange=0;
					for (int iter=0;iter<glaciers[glc]->numIter;iter++)
					{
						// internally, the addChangesToMarginCells function checks, if the total change is reached and will then return.
						// Thus, remainingChange will always be synchronous with the internal 
						remainingChange+=glaciers[glc]->addChangesToMarginCells(changePerCell,iter);
					}
				}
				if (fabs(remainingChange)>0.00001)
				{
					glaciers[glc]->cnt=0;
					glaciers[glc]->Area=0;
					glaciers[glc]->Hmax=-10000;
					glaciers[glc]->Hmin=10000;
					if (remainingChange>0) // positive change left: at least on elevation band must be extended by new margin cells
					{
						(*WSstream) << "DEBUG: " 
							<< "growing the glacier..."
							<<"[end-entry]";
						CreateNewGlacierCells(glaciers[glc]->ID);
					}
					else
					{
						(*WSstream) << "DEBUG: " 
							<< "shrinking the glacier..."
							<<"[end-entry]";
					}
// when all margin cells must be filled before a new margin is created
//					if (fabs(remainingChange-oldRemainingChange)<0.00001) // only if all changes were distributed and no margin cell is left...
//					{
//						CreateNewGlacierCells(glaciers[glc]->ID);
//					}
					UpdateGlacierStatisticsZones(glaciers[glc]->ID);
					if (glaciers[glc]->Hmin==10000) // no margin cells with no firn found (special case!)
					{
						// no firn free ice cell found for this glacier
						glaciers[glc]->cnt=0;
						glaciers[glc]->Area=0;
						glaciers[glc]->Hmax=-10000;
						glaciers[glc]->Hmin=10000;
						UpdateGlacierStatisticsZones(glaciers[glc]->ID,true); // shrink even the firn covered cells (it may be old firn from some years ago)
					}
					glaciers[glc]->createElevationBelts();
					CountMarginsGroupedByElevationBelts(glaciers[glc]->ID);
				}
				// emergency exit
				cnt++;
				if (cnt>30) 
				{
					(*WSstream) << "WARNING: \n----------- WARNING: ----------\nThere were too many iterations when trying to remove or add glacier areas for glacier" 
						        << glc << "\n-----------------------------\n" << "[end-entry]";
					break;
				}
			}
			(*WSstream) << "DEBUG: " 
				<< "---------------------------------\n\n"
				<<"[end-entry]";
			glaciers[glc]->setOldArea(glaciers[glc]->getArea());
			glaciers[glc]->setOldVolume(glaciers[glc]->getVolume());
			
		}
		// ToDo: MPI Scatter Grids or other records which were changed here
	}
}

// calculating elevation ranges for all or only a specific glacier
void schneemodellklasse::UpdateGlacierStatisticsZones(int id,bool useFirnasWell)
{
	if (ICE && glacmethode>10)
	{
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				// for glacierized cells only, the -1 is set when a new glacier cell should be created earlier in the FIRN moving sequence, 
				// so it is a flag for regarding this cell here as glacierized although there is no valid glacid yet
				if ((zn_routing[i][j]!=nodata)&&(GLACID[i][j]>=-1)) // changed at 2016-10-13 from GLACID[i][j]>=0 to GLACID[i][j]>=-1 
				{
					if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
					{
						int index=glacierindexmap.find((int)GLACID[i][j])->second;
						if ((id<0)||((int)GLACID[i][j]==id))
						{
								// if the actual ICE-value is -9999 but the GLACID is a known one, that can mean that the glacier is expanded into an origininally already occupied cell (that was part of this glacier). 
								// this can happen if the GLID grid from the initial start has a valid entry but the glc grid has not (nodata)
								// so we will set ICE[i][j] to the initzial value here, too, and change also the OldVolume etc.
								if (ICE[i][j] == nodata)
								{
									GLACIDOLD[i][j] = GLACID[i][j];
									// defining ICE[i][j] here
									ICE[i][j] = relSizeNewGlacier;

									if (useFirnasWell)
									{
										glaciers[index]->addCell(ICE[i][j], DHM[i][j], cellsize*cellsize, true); // always assume that this cell is firn free for creating elevation bands	
									}
									else
									{
										// here, we account for old firn only, all firn which was created from old snow of the actual year does not count for 
										// firn free or covered cells, since otherwise some fresh snow would also count for firn (all leftover snow from the actual year 
										// has already been converted to firn, so only the older layers of FIRN do count here)
										// glaciers[index]->addCell(ICE[i][j],DHM[i][j],cellsize*cellsize,(FIRN1[i][j]-FIRN[i][j][1])<=0);	

										// new for MPI optimized version: we avoid to Gather the firn stack by simply assuming that firn will be old usually
										glaciers[index]->addCell(ICE[i][j], DHM[i][j], cellsize*cellsize, false);
									}

									if (glaciers[index]->getVolume() > 0) // i.e. the glacier already exists
									{ // addig cells to an existing glacier means that the volume should also be extended by an amount equivalent to the added cells.
										// for a new glacier, this is done during the calculation of the Volume/Area relation, but for existing glaciers, the update must be daone "manually"
										double Vold = glaciers[index]->getOldVolume();
										double Aold = glaciers[index]->getOldArea();
										double ParB = glaciers[index]->ParB;
										double ParExp = glaciers[index]->ParExp;
										double V = glaciers[index]->getVolume();
										double expo = Aold + ICE[i][j] * (cellsize*cellsize / (double)1.0e6);
										double Vtemp = ParB*pow(expo>0?expo:0, ParExp);
										double Vextra = Vtemp - Vold;
										glaciers[index]->setOldArea(Aold + ICE[i][j] * (cellsize*cellsize / (double)1.0e6));
										glaciers[index]->setOldVolume(Vold + Vextra);
										glaciers[index]->setVolume(V + Vextra);
									}
								}
								else
								{
									if (useFirnasWell)
									{
										glaciers[index]->addCell(ICE[i][j], DHM[i][j], cellsize*cellsize, true); // always assume that this cell is firn free for creating elevation bands	
									}
									else
									{
										// here, we account for old firn only, all firn which was created from old snow of the actual year does not count for 
										// firn free or covered cells, since otherwise some fresh snow would also count for firn (all leftover snow from the actual year 
										// has already been converted to firn, so only the older layers of FIRN do count here)
										// glaciers[index]->addCell(ICE[i][j],DHM[i][j],cellsize*cellsize,(FIRN1[i][j]-FIRN[i][j][1])<=0);	

										// new for MPI optimized version: we avoid to Gather the firn stack by simply assuming that firn will be old usually
										glaciers[index]->addCell(ICE[i][j], DHM[i][j], cellsize * cellsize, false);
									}
// like for new glacier cells, the volume and area should now be made larger in order to welcome the old snow pack as a new part of the glacier
									if (glaciers[index]->getVolume() > 0) // i.e. the glacier already exists
									{ // addig cells to an existing glacier means that the volume should also be extended by an amount equivalent to the added cells.
										// for a new glacier, this is done during the calculation of the Volume/Area relation, but for existing glaciers, the update must be daone "manually"
										double Vold = glaciers[index]->getOldVolume();
										double Aold = glaciers[index]->getOldArea();
										double ParB = glaciers[index]->ParB;
										double ParExp = glaciers[index]->ParExp;
										double V = glaciers[index]->getVolume();
										double expo = Aold + ICE[i][j] * (cellsize * cellsize / (double)1.0e6);
										double Vtemp = ParB * pow(expo > 0 ? expo : 0, ParExp);
										double Vextra = Vtemp - Vold;
										glaciers[index]->setOldArea(Aold + ICE[i][j] * (cellsize * cellsize / (double)1.0e6));
										glaciers[index]->setOldVolume(Vold + Vextra);
										glaciers[index]->setVolume(V + Vextra);
									}

								}
						}
					}
					// -------------------------------- start new code fragment 2016-10-13 ----------------------------------------
					else
					{	// glacier not yet in glacierlist, so a new glacier is created --> an ICE[i][j] value must be set before calling this function, 
						// which is done when changing snow to firn and switching the firn layers and moving firn from non glacierized parts of glacier cells to the glacier 
						// defining GLACID[i][j] -> may be the next free ID of glacierlist or a glid of the neighbourhood
//					todo: if GLACID is - 1, then get the highest glacier id from the glaicerlist and add the new glacier with the next ID if no other glid call is within a cell or two radius around the actual cell
						if (GLACID[i][j] == -1)
						{
							if ((i > 1) && (GLACID[i - 1][j] > 0) && (zn_routing[i][j] == zn_routing[i-1][j])) { // north
								GLACID[i][j] = GLACID[i - 1][j];
							}
							else if ((j < ncols) && (GLACID[i][j + 1] > 0) && (zn_routing[i][j] == zn_routing[i][j+1])) { // east
								GLACID[i][j] = GLACID[i][j+1];
							}
							else if ((i < nrows) && (GLACID[i + 1][j] > 0) && (zn_routing[i][j] == zn_routing[i + 1][j])) { //south
								GLACID[i][j] = GLACID[i + 1][j];
							}
							else if ((j > 1) && (GLACID[i][j - 1] > 0) && (zn_routing[i][j] == zn_routing[i][j-1])) { // west
								GLACID[i][j] = GLACID[i][j - 1];
							}
							else if ((i>1) && (j < ncols) && (GLACID[i - 1][j + 1] > 0) && (zn_routing[i][j] == zn_routing[i - 1][j+1])) { // north east
								GLACID[i][j] = GLACID[i-1][j + 1];
							}
							else if ((j < nrows) && (j<ncols) && (GLACID[i + 1][j + 1] > 0) && (zn_routing[i][j] == zn_routing[i+1][j+1])) { // south east
								GLACID[i][j] = GLACID[i+1][j + 1];
							}
							else if ((j < nrows) && (j>1) && (GLACID[i + 1][j - 1] > 0) && (zn_routing[i][j] == zn_routing[i + 1][j-1])) { // south west
								GLACID[i][j] = GLACID[i+1][j - 1];
							}
							else if ((j > 1) && (i>1) && (GLACID[i - 1][j - 1] > 0) && (zn_routing[i][j] == zn_routing[i - 1][j-1])) { // north west
								GLACID[i][j] = GLACID[i-1][j - 1];
							}
							else
							{ // no glacier in the immediate neighbourhood, so a new glacier is created
								GLACID[i][j] = glacierindexmap.rbegin()->first + 1;
							}
						}
						GLACIDOLD[i][j] = GLACID[i][j];
						// defining ICE[i][j]
						ICE[i][j] = relSizeNewGlacier;
						// add glacier to list
						int totalids = glacierindexmap.size();
						// inserts a new element in the indexmap only when the GLACID is not yet found in the map
						if (glacierindexmap.find((int)GLACID[i][j]) == glacierindexmap.end())
						{
							glacierindexmap.insert(std::pair<int, int>((int)GLACID[i][j], (int)totalids));
							glaciers.push_back(new glacierclass((int)GLACID[i][j], VAscaling, VAexponent,
								0,   // actual Value
								0,   // actual Area
								numElevBands,   // iterations
								cellsize,
								extraWeightFactorBand0, // if growing and shrinking of the glacier should be overemphasized in the lowest elevation band. Default: 0
								(int)zn_routing[i][j]));
						}
						int index = glacierindexmap.find((int)GLACID[i][j])->second;

						if (glaciers[index]->getVolume() > 0) // i.e. the glacier already exists
						{ // addig cells to an existing glacier means that the volume should also be extended by an amount equivalent to the added cells.
							// for a new glacier, this is done during the calculation of the Volume/Area relation, but for existing glaciers, the update must be daone "manually"
							double Vold = glaciers[index]->getOldVolume();
							double Aold = glaciers[index]->getOldArea();
							double ParB = glaciers[index]->ParB;
							double ParExp = glaciers[index]->ParExp;
							double V = glaciers[index]->getVolume();
							double expo = Aold + ICE[i][j] * (cellsize*cellsize / (double)1.0e6);
							double Vtemp = ParB*pow(expo>0?expo:0, ParExp);
							double Vextra = Vtemp-Vold;
							glaciers[index]->setOldArea(Aold + ICE[i][j] * (cellsize*cellsize / (double)1.0e6));
							glaciers[index]->setOldVolume(Vold + Vextra);
							glaciers[index]->setVolume(V + Vextra);

						}

//							if (glaciers[index]->subbasin != zn_routing[i][j])
//							{
//								(*WSstream) << "FATAL: \n\n=======================================\n" <<
//									"!!!! ERROR !!!! \n\n" <<
//									"glacier " << (int)GLACID[i][j] << " is located in more than one subbasin:\n" <<
//									"First cell was found in subbasin " << glaciers[index]->subbasin << ",\n"
//									<< "whereas the cell in row " << i << " and col " << j << " belongs to subbasin " << zn_routing[i][j] << ".\n"
//									"Please make sure that each single glacier belongs to one specific subbasin only." << "[end-entry]";
//							}
//							throw(84);

						glaciers[index]->addCell(ICE[i][j], DHM[i][j], cellsize*cellsize, FIRN1[i][j] == 0);
						// glaciers[index]->addBalanceChangeToOldVolume((GMBO[i][j] - GMB[i][j])*ICE[i][j] / 1000.0*cellsize*cellsize / 1000000.0);
					}
					// -------------------------------- end new code fragment 2016-10-13 ----------------------------------------
				}
			}
		}
	}
}

// this method adds up all glacier mass balances from before resizing the glacier and redistributes it to all remaining cells after resizing the glacier in order to keep the balance correctly for the entire basin
void schneemodellklasse::RedistributeOldBalance(int id)
{
	if (ICE && glacmethode>10)
	{
		// initializing sums
		for (int glc=0;glc<glaciers.size();glc++)
		{
			if ((glaciers[glc]->ID==id)||(id<0)) // do it for the specific glacier only or for all glaciers (if no id was given or id was < 0)
			{
				glaciers[glc]->BalanceCarryOver=0;
				glaciers[glc]->TotalCellsBeforeResize=0; 
				glaciers[glc]->TotalCellsAfterResize=0; 
			}
		}

		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				if ((zn_routing[i][j]!=nodata)&&(GLACIDOLD[i][j]>=0)) // for glacierized cells only
				{
					if (glacierindexmap.find((int)GLACIDOLD[i][j])!=glacierindexmap.end())
					{
						int index=glacierindexmap.find((int)GLACIDOLD[i][j])->second;
						if ((id<0)||((int)GLACIDOLD[i][j]==id))
						{
							if (GLACIDOLD[i][j]==GLACID[i][j])
							{
								glaciers[index]->BalanceCarryOver		+= GMB[i][j]*(ICEOLD[i][j]-ICE[i][j]);
							}
							else // GLACID was set to 0 or even another glacier advanced into the cell formerly populated by the actual glacier, so we assume ICE[i][j] to be 0 for the actual cell and the actual glacier, i.e. all the old balance must be redistributed
							{
								glaciers[index]->BalanceCarryOver		+= GMB[i][j]*ICEOLD[i][j];
							}
							glaciers[index]->TotalCellsBeforeResize += ICEOLD[i][j];
						}
					}
				}
			}
		}

		// counting cells after resizing, using GLACID instead of GLACIDOLD
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				if ((zn_routing[i][j]!=nodata)&&(GLACID[i][j]>=0)) // for glacierized cells only
				{
					if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
					{
						int index=glacierindexmap.find((int)GLACID[i][j])->second;
						if ((id<0)||((int)GLACID[i][j]==id))
						{
							glaciers[index]->TotalCellsAfterResize  += ICE[i][j];
						}
					}
				}
			}
		}

		// now we now for every glacier the total amount of balance to redistribute to all cells and we can now define the redistribution amount per cell
		for (int glc=0;glc<glaciers.size();glc++)
		{
			if ((glaciers[glc]->ID==id)||(id<0)) // do it for the specific glacier only or for all glaciers (if no id was given or id was < 0)
			{
				if (glaciers[glc]->TotalCellsAfterResize > 0)
					glaciers[glc]->BalanceRedistribution = glaciers[glc]->BalanceCarryOver / glaciers[glc]->TotalCellsAfterResize;
				else
					glaciers[glc]->BalanceRedistribution = 0;
			}
		}

		// now redistribute the amount to all cells
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				if ((zn_routing[i][j]!=nodata)&&(GLACID[i][j]>=0)) // for glacierized cells only
				{
					if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
					{
						int index=glacierindexmap.find((int)GLACID[i][j])->second;
						if ((id<0)||((int)GLACID[i][j]==id))
						{
							GMB[i][j]+=glaciers[index]->BalanceRedistribution;
						}
					}
				}
			}
		}
	}
}


void schneemodellklasse::CreateGlacierStatisticsZones()
{
	glaciers.clear();
	if (ICE && glacmethode>10)
	{
		// creating map of glacierindices
		glacierindexmap.clear();
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
//				if ((zn_routing[i][j] != nodata) && (ICE[i][j]>0.0) && (GLACID[i][j]>0.0)) // for glacierized cells only
				if ((zn_routing[i][j] != nodata) && (GLACID[i][j]>0.0)) // for glacierized cells only
				{
					if (glacierindexmap.find((int)GLACID[i][j])==glacierindexmap.end())
					{ // a new glacier ID is found
						int totalids = glacierindexmap.size();
						glacierindexmap.insert(std::pair<int,int>((int)GLACID[i][j],(int)totalids));
						glaciers.push_back(new glacierclass((int)GLACID[i][j],VAscaling,VAexponent,
							0,   // actual Value
							0,   // actual Area
							numElevBands,   // iterations
							cellsize,
							extraWeightFactorBand0, // if growing and shrinking of the glacier should be overemphasized in the lowest elevation band. Default: 0
							(int)zn_routing[i][j]));
					}
					int index=glacierindexmap.find((int)GLACID[i][j])->second;
					if (glaciers[index]->subbasin != zn_routing[i][j])
					{
						(*WSstream) << "FATAL: \n\n=======================================\n"<<
						"!!!! ERROR !!!! \n\n"<<
						"glacier "<<(int)GLACID[i][j]<<" is located in more than one subbasin:\n"<<
						"First cell was found in subbasin "<<glaciers[index]->subbasin<<",\n"
						<<"whereas the cell in row "<<i<<" and col "<<j<<" belongs to subbasin "<<zn_routing[i][j]<<".\n"
						"Please make sure that each single glacier belongs to one specific subbasin only."<<"[end-entry]";
						throw(84);
					}
					glaciers[index]->addCell(ICE[i][j],DHM[i][j],cellsize*cellsize,FIRN1[i][j]==0);
					glaciers[index]->addBalanceChangeToOldVolume((GMBO[i][j]-GMB[i][j])*ICE[i][j]/1000.0*cellsize*cellsize/1000000.0);
				}
			}
		}
	}
}

void schneemodellklasse::CreateNewGlacierCells(int id)
{
	if (glacierindexmap.find(id)!=glacierindexmap.end())
	{
		return;
	}
	else
	{
		int glid=glacierindexmap.find(id)->second;
		glaciers[glid]->resetMarginCounter();
		bool oldglaciercell=false;

		// first run: grow into old glacier cells (which have ICE[i][j]=0 instead of -9999) and only for those cells which are firn free
		for (int i=1;i<=nrows;i++) // rows of the grid
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
			{
				bool margincell=false;
				if ((zn_routing[i][j]!=nodata)&& // valid cell
					(GLACID[i][j]<0)&& // and the actual cell is NOT part of any glacier yet --> IMPORTANT: clear this ID when removing a cell from the glacier
					(ICE[i][j]==0.0))   // which is not glacierized actually but was formerly (since ICE is not nodata)
				{
					// check if there is the glacier with code id to either one side of this cell. If yes, mark this cell for that glacier
					// and add the cell to the marginsCell-vector
					// conditions (here for northern cell):
					// cell to the north belongs to the same subbasin
					// cell to the north is glacierized
					// and belongs to the glacier of interest
					// and has no firn on it
					if ((i>1) && (zn_routing[i][j]==zn_routing[i-1][j]) && (ICE[i-1][j]==1.0) && (GLACID[i-1][j]==id) && (FIRN1[i-1][j]<0.0001)) { 
						margincell=true; // flag it for beeing registered as a margin cell later
					} else if ((j<ncols) && (zn_routing[i][j]==zn_routing[i][j+1]) && (ICE[i][j+1]==1.0) && (GLACID[i][j+1]==id) && (FIRN1[i][j+1]<0.0001))  { 
						margincell=true; 
					} else if ((i<nrows) && (zn_routing[i][j]==zn_routing[i+1][j]) && (ICE[i+1][j]==1.0) && (GLACID[i+1][j]==id) && (FIRN1[i+1][j]<0.0001))  {
						margincell=true; 
					} else if ((j>1) && (zn_routing[i][j]==zn_routing[i][j-1]) && (ICE[i][j-1]==1.0) && (GLACID[i][j-1]==id) && (FIRN1[i][j-1]<0.0001))  { 
						margincell=true;
					}
					if (margincell)
					{
						// even if no ICE=1 will be found as neighbor of any cell, we know that at least one cell is 
						// still to be filled, so NO other cell should be filled first					
						// this will not be reset for this glacier until the next grow command is executed
						if (ICE[i][j]==0.0) 
							oldglaciercell=true; 
						ICE[i][j]=0; // initialize the new glacierized cll with an icevalue of 0
						GLACID[i][j]=id; // make it part of the glacier of interest
						FIRN1[i][j]=0;
						for (int fl=0;fl<=firn_levels;fl++)
						{
							FIRN[i][j][fl]=0.0;
						}
						glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
					}
				}
			}
		}

		// second run (if first run retruned 0 cells) run: grow into old glacier cells (which have ICE[i][j]=0 instead of -9999), regardless if there is firn or not on the actual cell
		if (glaciers[glid]->totalMarginCells==0)
		{
			oldglaciercell=false;
			for (int i=1;i<=nrows;i++) // rows of the grid
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
				{
					bool margincell=false;
					double oldfirn=0;
					if ((zn_routing[i][j]!=nodata)&& // valid cell
						(ICE[i][j]==0.0))   // which is not glacierized actually
					{
						// check if there is the glacier with code id to either one side of this cell. If yes, mark this cell for that glacier
						// and add the cell to the marginsCell-vector
						if ((i>1) && (zn_routing[i][j]==zn_routing[i-1][j]) && (ICE[i-1][j]==1.0) && (GLACID[i-1][j]==id)) {  
							margincell=true; 
							oldfirn=FIRN1[i-1][j];
						} else if ((j<ncols) && (zn_routing[i][j]==zn_routing[i][j+1]) && (ICE[i][j+1]==1.0) && (GLACID[i][j+1]==id)) {  
							margincell=true; 
							oldfirn=FIRN1[i][j+1];
						} else if ((i<nrows) && (zn_routing[i][j]==zn_routing[i+1][j]) && (ICE[i+1][j]==1.0) && (GLACID[i+1][j]==id)) { 
							margincell=true; 
							oldfirn=FIRN1[i+1][j];
						} else if ((j>1) && (zn_routing[i][j]==zn_routing[i][j-1]) && (ICE[i][j-1]==1.0) && (GLACID[i][j-1]==id)) { 
							margincell=true; 
							oldfirn=FIRN1[i][j-1];
						}
						if (margincell)
						{
							if (ICE[i][j]==0.0) 
								oldglaciercell=true; 
							ICE[i][j]=0; // initialize the new glacierized cll with an icevalue of 0
							GLACID[i][j]=id; // make it part of the glacier of interest
							FIRN1[i][j]=oldfirn;
							for (int fl=0;fl<=firn_levels;fl++)
							{
								FIRN[i][j][fl]=oldfirn/firn_levels;
							}
							glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
						}
					}
				}
			}
		}

		// third run: grow anywhere into the same subbasin, regardless if the new cells were once glacier cells or not, do this for actual firn free cells only
		if ((glaciers[glid]->totalMarginCells==0)&&(oldglaciercell==false))
		{
			// only if no cells were found wich once were part of the glacier
			for (int i=1;i<=nrows;i++) // rows of the grid
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
				{
					bool margincell=false;
					if ((zn_routing[i][j]!=nodata)&& // valid cell
						(GLACID[i][j]<0)&& // and the actual cell is NOT part of any glacier yet --> IMPORTANT: clear this ID when removing a cell from the glacier
						(ICE[i][j]<=0.0))   // which is not glacierized actually
					{
						// check if there is the glacier with code id to either one side of this cell. If yes, mark this cell for that glacier
						// and add the cell to the marginsCell-vector
						// conditions (here for northern cell):
						// cell to the north belongs to the same subbasin
						// cell to the north is glacierized
						// and belongs to the glacier of interest
						// and has no firn on it
						if ((i>1) && (zn_routing[i][j]==zn_routing[i-1][j]) && (ICE[i-1][j]==1.0) && (GLACID[i-1][j]==id) && (FIRN1[i-1][j]<0.0001)) { 
							margincell=true; // flag it for beeing registered as a margin cell later
						} else if ((j<ncols) && (zn_routing[i][j]==zn_routing[i][j+1]) && (ICE[i][j+1]==1.0) && (GLACID[i][j+1]==id) && (FIRN1[i][j+1]<0.0001))  { 
							margincell=true; 
						} else if ((j<nrows) && (zn_routing[i][j]==zn_routing[i+1][j]) && (ICE[i+1][j]==1.0) && (GLACID[i+1][j]==id) && (FIRN1[i+1][j]<0.0001))  {
							margincell=true; 
						} else if ((j>1) && (zn_routing[i][j]==zn_routing[i][j-1]) && (ICE[i][j-1]==1.0) && (GLACID[i][j-1]==id) && (FIRN1[i][j-1]<0.0001))  { 
							margincell=true;
						}
						if (margincell)
						{
							ICE[i][j]=0; // initialize the new glacierized cll with an icevalue of 0
							GLACID[i][j]=id; // make it part of the glacier of interest
							FIRN1[i][j]=0;
							for (int fl=0;fl<=firn_levels;fl++)
							{
								FIRN[i][j][fl]=0.0;
							}
							glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
						}
					}
				}
			}
		}
		// this glacier has only firn covered cells, so we may create elevation bands in the firn covered area as a back door to keep the algorithm working
		// also theer are no old glacier cells, so let the glacier grow whereever possible
		if ((glaciers[glid]->totalMarginCells==0)&&(oldglaciercell==false))
		{
			for (int i=1;i<=nrows;i++) // rows of the grid
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
				{
					bool margincell=false;
					double oldfirn=0;
					if ((zn_routing[i][j]!=nodata)&& // valid cell
						(ICE[i][j]<=0.0))   // which is not glacierized actually
					{
						// check if there is the glacier with code id to either one side of this cell. If yes, mark this cell for that glacier
						// and add the cell to the marginsCell-vector
						if ((i>1) && (zn_routing[i][j]==zn_routing[i-1][j]) && (ICE[i-1][j]==1.0) && (GLACID[i-1][j]==id)) {  
							margincell=true; 
							oldfirn=FIRN1[i-1][j];
						} else if ((j<ncols) && (zn_routing[i][j]==zn_routing[i][j+1]) && (ICE[i][j+1]==1.0) && (GLACID[i][j+1]==id)) {  
							margincell=true; 
							oldfirn=FIRN1[i][j+1];
						} else if ((i<nrows) && (zn_routing[i][j]==zn_routing[i+1][j]) && (ICE[i+1][j]==1.0) && (GLACID[i+1][j]==id)) { 
							margincell=true; 
							oldfirn=FIRN1[i+1][j];
						} else if ((j>1) && (zn_routing[i][j]==zn_routing[i][j-1]) && (ICE[i][j-1]==1.0) && (GLACID[i][j-1]==id)) { 
							margincell=true; 
							oldfirn=FIRN1[i][j-1];
						}
						if (margincell)
						{
							ICE[i][j]=0; // initialize the new glacierized cll with an icevalue of 0
							GLACID[i][j]=id; // make it part of the glacier of interest
							FIRN1[i][j]=oldfirn;
							for (int fl=0;fl<=firn_levels;fl++)
							{
								FIRN[i][j][fl]=oldfirn/firn_levels;
							}
							glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
						}
					}
				}
			}
		}
	}
}

void schneemodellklasse::CountMarginsGroupedByElevationBelts(int id)
{
	for (int glc=0;glc<glaciers.size();glc++)
	{
		if ((glaciers[glc]->ID==id)||(id<0)) // do it for the specific glacier only or for all glaciers (if no id was given or id was 0)
			glaciers[glc]->resetMarginCounter();
	}

	for (int i=1;i<=nrows;i++) // rows of the grid
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
		{
			if ((zn_routing[i][j]!=nodata)&&(ICE[i][j]>0.0)&&(GLACID[i][j]>0.0)) // for glacierized cells only
			{
				// choose the correct glacier from the glaciers-array and check if this cell is a margin cell and if it is covered by Firn or not
				bool margincell=false;

				if (FIRN1[i][j]<=0) // only pure ice cells should be srink and grow
				// note: here, we possibly could use also the GMB grid: a negative mass balance would mean that there could be a fraction removed from this cell
				{
					if ((i>1)     && (GLACID[i-1][j  ]!=GLACID[i][j])) margincell=true;
					if ((j>1)     && (GLACID[i  ][j-1]!=GLACID[i][j])) margincell=true;
					if ((i<nrows) && (GLACID[i+1][j  ]!=GLACID[i][j])) margincell=true;
					if ((j<ncols) && (GLACID[i  ][j+1]!=GLACID[i][j])) margincell=true;
					if ((ICE[i][j]<1.0)&&(ICE[i][j]>=0.0)) margincell=true;
				
					if (margincell)
					{
						if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
						{
							int glid=glacierindexmap.find((int)GLACID[i][j])->second;
							if (((int)GLACID[i][j]==id)||(id<0)) // do it for the specific glacier only or for all glaciers (if no id was given or id was 0)
							{
								glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
							}
						}
						else
						{
//							(*WSstream) << "WARNING: Glacier with ID "<<(int)GLACID[i][j]<<" was not found in glacier list. This may happen if the GLID grid has a valid entry where the Glacierized fraction is 0."<<"[end-entry]";
//							throw(83);
						}
					}
				}
			}
		}
	}

	for (int glc=0;glc<glaciers.size();glc++)
	{
		if ((glaciers[glc]->totalMarginCells==0))
		{
			// this glacier has only firn covered cells, so we may create elevation bands in the firn covered area as a back door to keep the algorithm working
			for (int i=1;i<=nrows;i++) // rows of the grid
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) // columns of the grid
				{
					if ((zn_routing[i][j]!=nodata)&&(GLACID[i][j]>0.0)) // for glacierized cells only
					{
						// choose the correct glacier from the glaciers-array and check if this cell is a margin cell and if it is covered by Firn or not
						bool margincell=false;

						if ((i>1)     && (GLACID[i-1][j  ]!=GLACID[i][j])) margincell=true;
						if ((j>1)     && (GLACID[i  ][j-1]!=GLACID[i][j])) margincell=true;
						if ((i<nrows) && (GLACID[i+1][j  ]!=GLACID[i][j])) margincell=true;
						if ((j<ncols) && (GLACID[i  ][j+1]!=GLACID[i][j])) margincell=true;
						if ((ICE[i][j]<1.0)&&(ICE[i][j]>=0.0)) margincell=true; // all cells with < 1.0 icevalue will be handled as cells which may shrink or grow (important for growing since some new margin cells may be encircled from other valid cells after after growing)
					
						if (margincell)
						{
							if (glacierindexmap.find((int)GLACID[i][j])!=glacierindexmap.end())
							{
								int glid=glacierindexmap.find((int)GLACID[i][j])->second;
								if (((int)GLACID[i][j]==id)||(id<0)) // do it for the specific glacier only or for all glaciers (if no id was given or id was 0)
								{
									glaciers[glid]->addMarginCell(&(ICE[i][j]),&(GLACID[i][j]),DHM[i][j]); // we store the address of the cell for later access when changing it directly
								}
							}
							else
							{
//								(*WSstream) << "WARNING: Glacier with ID "<<(int)GLACID[i][j]<<" was not found in glacier list. This may happen if the GLID grid has a valid entry where the Glacierized fraction is 0."<<"[end-entry]";
								//throw(83);
							}
						}
					}
				}
			}
		}
	}
}

void schneemodellklasse::setLayeredSnowActive(bool flag, SOILtabletype *LayeredSnowProperties,double dTOffsetFreshSnow)
{ 
	if (do_it == 0) return;

	this->bLayeredSnowActive = flag; 
	if (bLayeredSnowActive)
	{
		this->LayeredSnowProperties = LayeredSnowProperties;
		
		// getting references to some stacks and grids from unsatzonmodel
		refSWES = gridinfo[getindex(_SnowWaterEquivalentSolid_)].stackpointer;
		refSWEL = gridinfo[getindex(_SnowWaterEquivalentLiquid_)].stackpointer;
		refUpperIndexSnow = gridinfo[getindex(_upperIndexSnow_)].pointer;
		refThickSnow = gridinfo[getindex(_SnowLayerThickness_)].stackpointer;
		refDensSnow = gridinfo[getindex(_SnowDensity_)].stackpointer;
		if (this->bHeatTransferActive)
		{
			refTemperSnow = gridinfo[getindex(_Soil_Temperature_Stack_Snow_)].stackpointer;
		}
		this->dTOffsetFreshSnow = dTOffsetFreshSnow;
	}
}

// ++++++++++++++++++++ End Snow model ++++++++++++++++++++++++++++++ 

// ++++++++++++++++++++ Start Glacierclass ++++++++++++++++++++++++++++++ 

// int    ID; // which glacier is it?
// double Area,      // area in km^2
//	      Volume,    // volume in 1000000 m^3
//	      oldArea,   // backup of old area
//	      oldVolume; // backup of old volume 
// double ParB; // 28.5; in Volume = ParB*Area*exp(ParExp)
// double ParExp; // 1.36; in Volume = ParB*Area*exp(ParExp)
// double Hmin, // minimum elevation of the glacier tonge
//	      Hmax, // maximum elevation of the glaciers ablation area = maximum elevation of equilibrium line
// int    numIter; // count of elevation belts for iterative grow/erode
// std::vector<elevationbelttype> elevBands;

glacierclass::~glacierclass()
{
	for (int i=0;i<elevBands.size();i++)
	{
        	delete elevBands[i];
	}
        elevBands.clear();
}


glacierclass::glacierclass():
ID(0),Area(0),Volume(0),oldArea(0),oldVolume(0),
ParB(28.5),ParExp(1.36),Hmin(10000),Hmax(-10000),numIter(1),cnt(0),cellsize(100),extraWeightFactorBand0(0),subbasin(0)
{
}

glacierclass::glacierclass(int id,double b,double exponent,double V,double A,int iter,double cs,int f,int sb):
ID(id),Area(A),Volume(V),oldArea(A),oldVolume(V),
ParB(b),ParExp(exponent),Hmin(10000),Hmax(-10000),numIter(iter),cnt(0),cellsize(cs),extraWeightFactorBand0(f),subbasin(sb)
{
	if (iter<1) iter=1; // at least one iteration must be done for grow/erode
}

void glacierclass::addCell(double fraction,double elevation,double cellarea,bool isIce) // cellarea must be given in m^2, isIce is used for calculating Hmax of firn free ice
{
	cnt++;
	Area+=fraction*cellarea/(double)1000000.0; // result in km^2
	if (cnt==1)
	{
		Hmax=Hmin=elevation;
	}
	else
	{
		if ((elevation>Hmax)&&(isIce==true)) {
			Hmax=elevation;
		}	else if ((elevation<Hmin)&&(isIce==true)) {
			Hmin=elevation;
		}
	}
}

// for initialization of the glaciers we must know the old area and volume, if possible
// Thus, we create the difference between the actual mass balance and the old mass balance (valid for the last area calculation)
// so we can determine the old value and thus the old area indirectly
void glacierclass::addBalanceChangeToOldVolume(double dV) 
{
	// remember: gmb and gmbo have measures of mm whereas oldVolume and Volume have measures of 1000000 m^3
	this->oldVolume+=dV; 
}

void glacierclass::createElevationBelts()
{
	// we know now the minimum and maximum elevation for the firn-free glacier and thus we can create elevation belts for grow/erode algorithm
	elevBands.clear();
	double Hrange=(Hmax-Hmin)/numIter;
	for (int i=0;i<numIter;i++)
	{
		elevBands.push_back(new elevationbelt(Hmin+Hrange*i,Hmin+Hrange*(i+1),Hrange));
	}
}

void glacierclass::resetMarginCounter()
{
	for (int cnt=0;cnt<numIter;cnt++)
	{
		elevBands[cnt]->clearMarginCells();
	}
	totalMarginCells=0;
}

void glacierclass::addMarginCell(double *cell,double *id,double elevation)
{
	// look for the correct elevation band
	for (int band=0;band<elevBands.size();band++)
	{
		if ((elevBands[band]->isInRange(elevation))
			|| 
			((elevation>=Hmax)&&(band==elevBands.size()-1))
			||
			((elevation<=Hmin)&&(band==0))
			)
		{
			elevBands[band]->addMarginCell(cell,id);
			totalMarginCells++;
			break;
		}
	}
}

double glacierclass::getTotalChangeInCells()
{
	return totalChange;
}

// get the amount in cellsize units which must be added or 
// subtracted from the actual cell, depending on elevation and the iteration
void glacierclass::calcTotalChangeInCells()
{
	this->remainingChanges=0;
	totalChange=(this->Area-this->oldArea)*(double)1000000.0/((double)cellsize*(double)cellsize);
}


long glacierclass::getChangeableCellCount(double Change)
{
	long sumMarginCells=0; 
	for (int band=0;band<elevBands.size();band++)
	{
		sumMarginCells += elevBands[band]->getChangeableMarginCells(Change);
	}
	return sumMarginCells;
}


double glacierclass::changePerCell(double remainingChange,bool flag)
{
	double changeCnt;
	if (!flag) 
		changeCnt=getTotalChangeInCells(); 
	else 
		changeCnt=remainingChange;

	double sumMarginCells=0;
	
	for (int band=0;band<elevBands.size();band++)
	{
		sumMarginCells += elevBands[band]->getChangeableMarginCells(changeCnt)*(elevBands.size()-band);
	}
	sumMarginCells += elevBands[0]->getChangeableMarginCells(changeCnt)*extraWeightFactorBand0*elevBands.size(); // add much more extra weight for band 0 as read from control file

	return changeCnt/sumMarginCells;
}

// takes the change per margin cell and returns the remaining changes (if there were not enough margin cells - another run must follow) 
// iter defines, which elevation belts should be part of the actual action
double glacierclass::addChangesToMarginCells(double fraction,int iter)
{
	remainingChanges=0;
	for (int band=0;band<=iter;band++)
	{
		remainingChanges += elevBands[band]->addChangesToMarginCells((double)fraction);
// extra weight for band 0 as read from control file
		if (band==0)
		{
			for (int extra=1;extra<=extraWeightFactorBand0;extra++)
			{
				remainingChanges += elevBands[band]->addChangesToMarginCells((double)fraction);
			}
		}
	}
	return remainingChanges;
}


void glacierclass::addBalanceChangeToVolume(double dV)
{
	this->Volume += dV;
}


void glacierclass::calcVolume()
{
//	const double Amin = 0.001;
	const double Amin = 0.01;
	double Vmin = ParB*pow(Amin, ParExp);

	if (Area <= Amin)
	{
		Volume = Vmin; // equivalent to an area of Amin
		Area = Amin;
	}
	Volume = ParB*pow(((Area>0)?Area:Amin), ParExp); // result in 1000000 m^3

}

void glacierclass::calcArea()
{
	// if volume is 0, this is a newly created glacier and the vlume can now be computed from the initial area to prevent the resizing in this first year
	// It may also be a good idea to not limit the volume and area at all or at least to much smaller values and also to do it depending on the ParA and ParB-Parameters
//	const double Amin = 0.001;
	const double Amin = 0.0543;
	double Vmin = ParB*pow(Amin, ParExp);

	if (this->Volume == 0) // freshly created glacier, so set Volume depending on actual area which is a known entity
	{
		Volume = ParB*pow(Area>0?Area:Amin, ParExp); // result in 1000000 m^3
		oldVolume = Volume;
		oldArea = Area;
	}

	if (Volume <= Vmin)
	{
		Volume = Vmin; // equivalent to an area of Amin
		Area = Amin;
	}
	else
	{
		Area = pow((((Volume>0)?Volume:0) / ParB), (1.0 / ParExp)); // result in km^2
	}
}

void glacierclass::calcOldArea()
{
	// if volume is 0, this is a newly created glacier and the vlume can now be computed from the initial area to prevent the resizing in this first year
	// It may also be a good idea to not limit the volume and area at all or at least to much smaller values and also to do it depending on the ParA and ParB-Parameters
	//	const double Amin = 0.001;
	const double Amin = 0.0543;
	double Vmin = ParB*pow(Amin, ParExp);

	if (this->Volume == 0) // freshly created glacier, so set Volume depending on actual area which is a known entity
	{
		oldVolume = ParB*pow(oldArea>0?oldArea:0, ParExp); // result in 1000000 m^3
	}

	if (oldVolume <= Vmin)
	{
		oldVolume = Vmin; // equivalent to an area of Amin
		oldArea = Amin;
	}
	else
	{
		oldArea = pow((((oldVolume>0)?oldVolume:0) / ParB), (1.0 / ParExp)); // result in km^2
	}
}

void glacierclass::serialize() 
{
	// writing the content of variables to a file
	// still a dummy right now
}

double glacierclass::getArea()
{
	return Area;
}

double glacierclass::getVolume()
{
	return Volume;
}

double glacierclass::getOldArea()
{
	return oldArea;
}

double glacierclass::getOldVolume()
{
	return oldVolume;
}

int glacierclass::getExtraWeightFactorBand0()
{
	return extraWeightFactorBand0;
}

void glacierclass::setExtraWeightFactorBand0(int f)
{
	extraWeightFactorBand0=f;
}

void glacierclass::setArea(double a)
{
	Area = a;
}

void glacierclass::setVolume(double v)
{
	Volume=v;
}

void glacierclass::setOldArea(double a)
{
	oldArea=a;
}

void glacierclass::setOldVolume(double v)
{
	oldVolume=v;
}

void elevationbelt::clearMarginCells() 
{
	NumMarginCells=0;
	marginCells.clear();
	marginCellsIDs.clear();
	inactiveCells.clear();
}

double elevationbelt::getChangeableMarginCells(double change)
{
	double counter=0;
	if (change!=0)
	{
		for (int nc=0;nc<this->marginCells.size();nc++)
		{
			if (((*marginCells[nc])>0.0)&&(change<(double)0.0)) { // negative change: margin cells must have a value greater than 0
				counter+=1.0;
			} else if (((*marginCells[nc])<(double)1.0)&&(change>(double)0.0)) { // positive change: margin cells must have a value less than 1
				counter+=1.0;
			} else { 
				// flag this cell as not changeable until the next margin is created
				inactiveCells[nc]=1;
			}
		}
	}
	return counter;
}

double elevationbelt::addChangesToMarginCells(double changePerCell) 
{
	double remainingChanges=0;
	long activeCells=0;
	for (int mc=0;mc<inactiveCells.size();mc++)
	{
		if (!inactiveCells[mc]) activeCells++;
	}
	remainingChanges = activeCells*changePerCell;
	double oldRemainingChanges=remainingChanges;
	int iterations=0;
	// recursive: if some of the margin cells were filled completely (change > 0) or meptied and removed completely
	// then the changePerCel will not result in a remainingChange of 0
	// Thus, we try to distribute this remaining changes to all other cells of this elevation band
	// if ALL cells of that band are filled or removed, then the remainig changes will not change during one loop and 
	// the loop can be broken (additional cells will be added to the glacier or the new margin cells will be registerd 
	// and the this procedure is repeated
	int cnt=activeCells;
	do
	{ 
		oldRemainingChanges=remainingChanges;
		double internalChangePerCell;
		if (cnt>0)
			internalChangePerCell=remainingChanges/(double)cnt;
		else
			internalChangePerCell=0;
		cnt=0;
		for (int i=0;i<marginCells.size();i++)
		{
			if (inactiveCells[i]==0)
			{
				double temp=(*marginCells[i])+internalChangePerCell;
				if ((temp<0.0)&&((*marginCellsIDs[i])>=0)) { // do not remove cells twice
					(*marginCells[i])=0;        // direct access to ICE[i][j]
					(*marginCellsIDs[i])=-9999; // direct access to GLACID[i][j] -> remove this cell from the glacier so it will not be listed in the matgincells array again
					remainingChanges -= temp;
				} else if ((temp>1.0)&&((*marginCellsIDs[i])>=0))  {
					(*marginCells[i])=1.0;      // direct access to ICE[i][j]
					remainingChanges -= (internalChangePerCell-(temp-1.0));
				} else if ((*marginCellsIDs[i])>=0) 	{
					(*marginCells[i])=temp;    // direct access to ICE[i][j]
					remainingChanges -= internalChangePerCell;
					cnt++; // is the number of cells which have still the capacity of removing or adding a fraction of ice -> used for iterative runs when a recalculation of the change per cell is required
				}
			}
		}
		iterations++;
	} while ((fabs(remainingChanges)>0.00001)&&(fabs(oldRemainingChanges-remainingChanges)>0.00001)&&(iterations<5));
	if (iterations>4)
	{
		std::cerr << "error: could not remove or add the required number of cells for the glacier" << std::endl;
	}

	return remainingChanges;
}

void elevationbelt::addMarginCell(double *cellvalue,double *cellid) 
{
	NumMarginCells++; 
	marginCells.push_back(cellvalue); 
	marginCellsIDs.push_back(cellid); 
	inactiveCells.push_back(0);
}

bool elevationbelt::isInRange(double elev) 
{
	if ((elev<Hmax)&&(elev>=Hmin)) 
		return true; 
	else 
		return false; 
}


