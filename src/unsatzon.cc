/*

(c) 1995-2017 by Joerg Schulla
more version history: see source safe comments

unsaturated zone water and substance flows, discharge generation (direct flow, interflow, baseflow)
heat transfer, snow melt (when included in layered snow model and an energy balance approach for snow melt is used in snowmod.cc)

author     : Joerg Schulla
method     : Richardsequation, FDM, Pedotransferfunctions:  van-Genuchten


*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
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


/* ******************** some functions for internal organisation ********************* */
//#pragma optimize ("p",off)

/* a function wich analyses the flow travel time grid (the time, which
surface ruoff from a grid cell takes to reach the next subbasin outlet)
as a result there will allocated a set of statistic registers (for each
subbasin and each flow travel time zone one statistic register).
The distribution of travel zones is printed to the log file */
//#pragma optimize ("p",on)

#ifdef _MPI
void zonesumunsat(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype)
{
	for (int i = 0; i<*len; i++)
	{
		inout[i].n += in[i].n;
		inout[i].sx = 0;
		//		inout[i].rel+=in[i].rel; --> do not add the already calculated relative values which must remain constant
		//		inout[i].zonenwert+=in[i].zonenwert; --> also the zonenwert is a constant
		inout[i].RelCellSum += in[i].RelCellSum;
	}
}
#endif

void unsatzonmodellklasse::fzt_grid_auswerten()
{
	int logdateiopen;
	int i,j,id1,id2,z;
	FILE *logdatei;
	double cellcount;

	/* searching in the flow time grid for the longest flow time in order to reserve enough memory (in steps of time intervals) */
	(*WSstream) << "looking for longest flow travel time"<<"[end-entry]";
	maxzeit=(double*)calloc(zonmax+1,sizeof(double));
	for (int i=1;i<=zonmax;i++) maxzeit[i]=0.0;
	if (FZT)
	{
		for (int i=mpistartrow;i<=mpiendrow;i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					FZT[i][j] = ceil(FZT[i][j]/60.0/InV);
					if (FZT[i][j] > maxzeit[zarray[(int)zn[i][j]]]) 
					{
						maxzeit[zarray[(int)zn[i][j]]]=FZT[i][j];
					}
				}
			}
		}
	} 
	else
	{
		for (int i=1;i<=zonmax;i++) 
		{
			maxzeit[i]=1.0;
		}
	}

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	double *maxzeit2 = new double[zonmax + 1];
	MPI_Allreduce(&(maxzeit[1]), &(maxzeit2[1]), zonmax, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
	for (int i=1;i<=zonmax;i++) 
	{
		maxzeit[i]=maxzeit2[i];
	}
	delete[] maxzeit2;
	MPI_Barrier(MPI_COMM_WORLD);
	(*WSstream) << "done and reduced to rank 0" << "[end-entry]";
#else
	(*WSstream) << "done" << "[end-entry]";
#endif

	/* per subbasin: allocating an array of statistk-registers as many as specified in "maxzeit[i]" plus one as a buffer,
because the actual surface flow is routed through these registers by pointer rotating */

	zoneFZT_QD = (statrecord***)calloc(tracer_count+1,sizeof(statrecord**));
	zoneFZT_QI = (statrecord***)calloc(tracer_count+1,sizeof(statrecord**));
	for (int z=0;z<=tracer_count;z++)
	{
		zoneFZT_QD[z] = (statrecord**)calloc(zonmax+1,sizeof(statrecord*));
		zoneFZT_QI[z] = (statrecord**)calloc(zonmax+1,sizeof(statrecord*));
		for (int i=1;i<=zonmax;i++)
		{
			zoneFZT_QD[z][i]=(statrecord*)calloc((int)maxzeit[i]+1,sizeof(statrecord));
			zoneFZT_QI[z][i]=(statrecord*)calloc((int)maxzeit[i]+1,sizeof(statrecord));
		}
	}

	/* calculating relative fractions of each travel zone for all subbasins */
	if (FZT)
	{
		for (int i = 1; i <= zonmax; i++)
		{
			for (int j = (int)maxzeit[i]; j >= 1; j--)
			{
				zoneFZT_QD[0][i][j].n = 0;
				zoneFZT_QD[0][i][j].RelCellSum = 0;
				zoneFZT_QI[0][i][j].n = 0;
				zoneFZT_QI[0][i][j].RelCellSum = 0;
			}
		}
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			(*WSstream) << ".";
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					id1=zarray[(int)zn[i][j]];
					id2=(int)FZT[i][j];
					if (RelCellArea!=NULL)
					{
						zoneFZT_QD[0][id1][id2].n += RelCellArea[i][j];
						zoneFZT_QD[0][id1][id2].RelCellSum += RelCellArea[i][j];
						zoneFZT_QI[0][id1][id2].n += RelCellArea[i][j];
						zoneFZT_QI[0][id1][id2].RelCellSum += RelCellArea[i][j];
					}
					else
					{
						zoneFZT_QD[0][id1][id2].n++;
						zoneFZT_QD[0][id1][id2].RelCellSum++;
						zoneFZT_QI[0][id1][id2].n++;
						zoneFZT_QI[0][id1][id2].RelCellSum++;
					}
				}
			}
		}
		(*WSstream) << "[end-entry]";

#ifdef _MPI

		// reducing all the arrays to a common array and bcasting the result
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Datatype zonevectortype;
		MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
		MPI_Type_commit(&zonevectortype);

		void(*zonesumfunction)(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype) = &zonesumunsat;
		MPI_Op mySumOp;
		MPI_Op_create((MPI_User_function*)zonesumfunction, true, &mySumOp);

		for (int i = 1; i <= zonmax; i++)
		{

			statrecord *tempzone;
			tempzone = new statrecord[(int)maxzeit[i] + 1];

			for (int j = 0; j <= (int)maxzeit[i]; j++)
			{
				tempzone[j].n = tempzone[j].RelCellSum = 0;
			}

			MPI_Allreduce(&(zoneFZT_QD[0][i][0]), tempzone,(int)(maxzeit[i] + 1), zonevectortype, mySumOp, MPI_COMM_WORLD);

			for (int j = 0; j <= maxzeit[i]; j++)
			{
				zoneFZT_QD[0][i][j].n = tempzone[j].n;
				zoneFZT_QD[0][i][j].RelCellSum = tempzone[j].RelCellSum;
				zoneFZT_QI[0][i][j].n = tempzone[j].n;
				zoneFZT_QI[0][i][j].RelCellSum = tempzone[j].RelCellSum;
			}

			delete[] tempzone;
		}

		MPI_Op_free(&mySumOp);
		MPI_Type_free(&zonevectortype);

		MPI_Barrier(MPI_COMM_WORLD);
#endif
		
		logdateiopen=false;
		if (mpi_rank==0)
		{
			char logfilename[1024];
			sprintf(logfilename,"modell_%d.log",mpi_rank);
			if ((logdatei = fopen(logfilename,"a"))==NULL)
			{
				(*WSstream) << "WARNING: error opening log-file "<<logfilename<<"[end-entry]";
				logdateiopen=false;
			}
			else logdateiopen=true;
			fprintf(logdatei,"flow time distribution:\n");
		}

		for (int i=1;i<=zonmax;i++)
		{
			if (mpi_rank==0)
			{
				if (logdateiopen) fprintf(logdatei,"TG\t%4d\t(%3d\tZonen):",(int)zone[i].zonenwert,(int)maxzeit[i]);
			}
			/* how many valid values are in the flow time grid */
			cellcount=0;
			for (int j = 1; j <= (int)maxzeit[i]; j++)
			{
				cellcount += zoneFZT_QD[0][i][j].n;
			}

			if (cellcount==0)
			{
				(*WSstream) << "FATAL: all cells in "<<gridinfo[getindex(_Fliesszeiten_)].name
				<< ", subbasin "<<zone[i].zonenwert<<" are NODATA!"<<"[end-entry]";
				throw(6);
			}
			/* relativ fractions of each travel zone */
			for (int z=0;z<=tracer_count;z++)
			{
				for (int j=1;j<=(int)maxzeit[i];j++)
				{
					zoneFZT_QI[z][i][j].rel = zoneFZT_QD[z][i][j].rel = zoneFZT_QD[0][i][j].n/cellcount;
					zoneFZT_QI[z][i][j].zonenwert = zoneFZT_QD[z][i][j].zonenwert = j;
					if (z==0)
					{
						if (mpi_rank==0)
						{
							if (logdateiopen) fprintf(logdatei,"\t%5.3f",zoneFZT_QD[0][i][j].rel);
						}
						if (zoneFZT_QD[0][i][1].rel<=0)
						{
							(*WSstream) << "FATAL: in flow-travel-time grid: the first time zone of subbasin "<< zone[i].zonenwert
							<<" has no cells!"<<"[end-entry]";
							throw(61);
						}
						zoneFZT_QD[0][i][j].n = 0;
					}
				}
			}
			if (mpi_rank==0)
			{
				if (logdateiopen) fprintf(logdatei,"\n");
			}
		}
		if (mpi_rank==0)
		{
			if (logdateiopen) fclose(logdatei);
		}
	}
	else
	{
		for (int z = 0; z <= tracer_count; z++)
		{
			for (int i = 1; i <= zonmax; i++)
			{
				zoneFZT_QD[z][i][1].rel = 1.0;
				zoneFZT_QI[z][i][1].rel = 1.0;
			}
		}
		if (mpi_rank==0)
		{
			char logfilename[1024];
			sprintf(logfilename,"modell_%d.log",mpi_rank);
			if ((logdatei = fopen(logfilename,"a"))==NULL)
			{
				(*WSstream) << "WARNING: log file "<<logfilename<<" could not be opened"<<"[end-entry]";
			}
			else
			{
				fprintf(logdatei,"flow time distribution will not be calculated since no flow time grid was provided\n");
				fclose(logdatei);
			}
		}
	}

	/* kumulieren der relativen antiele der zonen, so dass die erste zone immer alle zellen enthält */
	for (int i = 1; i <= zonmax; i++)
	{
		for (int j = (int)maxzeit[i]; j>1; j--)
		{
			zoneFZT_QD[0][i][j - 1].n += zoneFZT_QD[0][i][j].n;
			zoneFZT_QD[0][i][j-1].RelCellSum += zoneFZT_QD[0][i][j].RelCellSum;
			zoneFZT_QI[0][i][j-1].n += zoneFZT_QI[0][i][j].n;
			zoneFZT_QI[0][i][j-1].RelCellSum += zoneFZT_QI[0][i][j].RelCellSum;
		}
	}

}

// if lakes are used, some grid values are changed (no interflow, no drainage, no slope a.s.o)
// during the model run, the lakes are not handled in any other way than the other cells: the POND-grid must be active and this 
// was implemented long before...
void unsatzonmodellklasse::initialize_lakes(std::string instanzname)
{
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
                        fclose(steuerdatei);
			throw(28);
		}

		strcpy(tstr,instanzname.c_str());

		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		findeabschnitt(tstr);

		fscanfln(steuerdatei,"%d",&lake_active,"flow control parameter for ",tstr);
		// switch off lake model when in Spinup-mode
		if (doSpinUp)
		{
			lake_active=0;
		} 

		if (lake_active != 0)
		{
			if (!pond_use)
			{
				(*WSstream) << "ERROR:\n\n=======================================\n"
				<< "!!!! ERROR !!!! \n\n"
				<< "Error in initialization method for "<<instanzname<<".\n"
				<< "pond_use must be set in "<<this->instanzname<<"\n"
				<< "for using the "<<instanzname
				<< "\n\n=======================================\n"<< "[end-entry]";
                                fclose(steuerdatei);
				throw(215);
			}

			fscanfln(steuerdatei,"%d",&lake_ini_method,"initialization method for ",tstr);
			if ((fclose(steuerdatei))==EOF)
			{
				(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"after reading initialization method for "<<instanzname<<"[end-entry]";
                                fclose(steuerdatei);
				throw(11);
			}

			LAKE=gridinfo[getindex(_LAKE_CODES_)].pointer;
			if ((lake_ini_method<1)||(lake_ini_method>2))
			{
				(*WSstream) << "FATAL: \n\n=======================================\n"<<
				"!!!! ERROR !!!! \n\n "
				<< "Error in initialization method for "<<instanzname<<".\n"
				<< "Allowed values for lake-ini-method are 1 (no change of grids) or 2 (adjust DEM, SLOPE, DRAIN).\n"
				<< "In either case, theta_soil will be set to saturation,\n"<<
				"if the cell is marked as lake and POND is > 0. Also, the groudwater level\n"<<
				"will be set to 0." <<
				"\n\n=======================================\n"<< "[end-entry]";
				fclose(steuerdatei);
				throw(215);
			}

			// change some values as initialization for lake use
			if (lake_ini_method==1) // only change saturation values
			{
				for (int i = mpistartrow; i <= mpiendrow; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						{
							if (LAKE[i][j] >= 0)
							{
								if ((gridslesen==0)||((gridslesen>0)&&(POND[i][j]>1.0)))
								{
									// if the lake will be initialized not by old grids but from scratch, we can set GWSTAND and GWLEVEL as well as thetas values as saturated
									GWSTAND[i][j]=0; // saturated soil
									GWLEVEL[i][j]=1; // GW is in first soil level
									for (int z=1;z<levels;z++)
									{
										theta[i][j][z]=ST[(int)SOIL[i][j]].getThetas(z);
									}
								}
								else
								{
									// when a lake cell is fallen dry (Pond contains no water any longer), 
									// we will take the values from old grids as they are without changes to 
									// GWSTAND, GWLEVEL and theta values
									// so we do nothing here...
								}
							}
						}
					}
				}
			}
			else if (lake_ini_method==2)
			{
				for (int i = mpistartrow; i <= mpiendrow; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						{
							if (LAKE[i][j] >= 0)
							{
								if ((gridslesen==0)||((gridslesen>0)&&(POND[i][j]>0.3)))
								{
									// if the lake will be initialized not by old grids but from scratch, we can set GWSTAND and GWLEVEL as well as thetas values as saturated
									if (HDRAIN)
									{
										DRAINLEVEL[i][j]=-1; // no drainage, since drainllevel -1 will never be reached
									}
									SLP[i][j]=0; // no interflow
									if (IRRIG)
									{
										IRRIG[i][j]=0; // no irrigation will take place in a lake
									}
									GWSTAND[i][j]=0; // saturated soil
									GWLEVEL[i][j]=1; // GW is in first soil level
									for (int z=1;z<levels;z++)
									{
										theta[i][j][z]=ST[(int)SOIL[i][j]].getThetas(z);
									}
									// correction for elevation model DEM/DHM, always stored in gridinfo[2]
									gridinfo[2].pointer[i][j]-=MAXPOND[i][j]/1000.0;
								}
								else
								{
									// when a lake cell is fallen dry (Pond contains no water any longer), 
									// we will take the values from old grids as they are without changes to 
									// GWSTAND, GWLEVEL and theta values
									// so we do nothing here...
									gridinfo[2].pointer[i][j]-=MAXPOND[i][j]/1000.0;
								}
							}
						}
					}
				}
				MPIGatherGrid(2,gridinfo[2].pointer); // changed DHM cells must be available in rank 0 for initialing lakes
			}
		}
		delete[] tstr;
	}
}

//#pragma optimize ("p",off)

/* this is the initialization of the unsaturated zone model. The function is called only one time at model initialization.
All Names and parameters from the control file are read in. All initialization is done here. */
//#pragma optimize ("p",on)
unsatzonmodellklasse::unsatzonmodellklasse(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
	double breitearg,double laengearg,double bezugsmeridianarg,const char *bezeichnung,
	char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, 
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll,double p_yll,
	double p_cellsize,double p_nodata,long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,irrigationtabletype *p_irrigationtable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,
	std::vector< std::vector<RiverLinkCell> > *p_RivLnkCellsToSend,
	std::vector< std::vector<RiverLinkCell> > *p_RivLnkCellsToReceive,
	std::vector<RiverLink*> *p_RiverLinkVector, int substepcnt) :
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr, p_LUtable, p_SOILtable, p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode,p_vHRU,substepcnt),
storagefile(0),storagefilename(0),
pond_use(0),drain_use(0),clay_use(0),surface_interaction(0),gwflow_active(0),glacierlist(0),
riverflow_table(0),qubuffer(0),thetabuffer(0),hhydbuffer(0),otherdatabuffer(0),concbuffer(0),etrdetailbuffer(0),interflowbuffer(0),
additionals(0),qubufmax(0),qubufcnt(0),thetabufmax(0),thetabufcnt(0),hhydbufmax(0),otherdatabufmax(0),hhydbufcnt(0),otherdatabufcnt(0),etrdetailbufcnt(0),interflowbufcnt(0),
concbufmax(0),concbufcnt(0),logdatei(0),thetadatei(0),qudatei(0),hhyddatei(0),otherdatadatei(0),etrdetaildatei(0),interflowdatei(0),concdatei(0),tempdatei(0),
thetadateiname(0),qudateiname(0),hhyddateiname(0),otherdatadateiname(0),etrdetaildateiname(0),interflowdateiname(0),concdateiname(0),control_row(0),control_col(0),
ST(0),P(0),LUSE(0),SOIL(0),WURZ(0),VBG(0),ICE(0),ICEOLD(0),FZT(0),ETP(0),ETR(0),EVAP(0),EVAR(0),SNOW(0),IRRIG(0),HDRAIN(0),
DRAINLEVEL(0),DRAINDIST(0),SLP(0),GWSTAND(0),GWN(0),PERCOL(0),UPRISE(0),TANBETA(0),QDRAIN(0),QDgrid(0),QIgrid(0),DELTAGW(0),
oldDELTAGW(0),GWLEVEL(0),GWtheta(0),INFEX(0),RIVDEP(0),RIVWIT(0),RIVLNK(0),MAXPOND(0),POND(0),QDSU(0),
CLAY(0),SB_1_(0),SB_2_(0),THROOT(0),SBLOW(0),SBHIGH(0),IRRIGCODE(0),ETR2ETP(0),MAXPUMP(0),GWIN(0),GWEX(0),KOLM(0),S01(0),
theta(0),theta1(0),hh(0),hh1(0),hgeo(0),hgeo1(0),qu(0),qu1(0),
conc(0),conc1(0),gw_thick(0),gw_conc(0),iteration_stat(0),QDstat(0),
QIstat(0),QBstat(0),Qgstat(0),SB_1_stat(0),SB_2_stat(0),SB_2_IMMO_stat(0),QDRAINstat(0),GWSTANDstat(0),GWNstat(0),GWINstat(0),
GWEXstat(0),PONDstat(0),QDstatcode(0),QIstatcode(0),QBstatcode(0),Qgstatcode(0),GWINstatcode(0),GWEXstatcode(0),
QDRAINstatcode(0),GWSTANDstatcode(0),GWNstatcode(0),PONDstatcode(0),SB_1_statcode(0),SB_2_statcode(0),SB_2_IMMO_statcode(0),
iteration_statcode(0),iteration_statfilename(0),QDstatfilename(0),QIstatfilename(0),QBstatfilename(0),
Qgstatfilename(0),GWINstatfilename(0),GWEXstatfilename(0),QDRAINstatfilename(0),GWSTANDstatfilename(0),
GWNstatfilename(0),PONDstatfilename(0),SB_2_statfilename(0),SB_2_IMMO_statfilename(0),SB_1_statfilename(0),maxzeit(0),theta_stackname(0),
hh_stackname(0),hgeo_stackname(0),qu_stackname(0),GWSTANDgridname(0),GWthetagridname(0),GWNgridname(0),UPRISEgridname(0),PERCOLgridname(0),
PONDgridname(0),GWLEVELgridname(0),QDRAINgridname(0),QDgridname(0),QIgridname(0),DELTAGWname(0),GWINgridname(0),
GWEXgridname(0),theta_writestack(0),theta_index(0),hh_writestack(0),hh_index(0),hgeo_writestack(0),hgeo_index(0),
qu_writestack(0),qu_index(0),GWSTANDwritegrid(0),gwstandindex(0),GWNwritegrid(0),gwnindex(0),PONDwritegrid(0),pondindex(0),
GWLEVELwritegrid(0),gwlevelindex(0),QDRAINwritegrid(0),qdrainindex(0),QDwritegrid(0),qdindex(0),QIwritegrid(0),qiindex(0),
DELTAGWwritegrid(0),deltagwindex(0),olddeltagwindex(0),GWINwritegrid(0),gwinindex(0),GWEXwritegrid(0),gwexindex(0),
GWthetawritegrid(0),gwthetaindex(0),DRAINLEVELwritegrid(0),drainlevelindex(0),tanbetaindex(0),maxpumpindex(0),
SB_1_gridname(0),SB_1_writegrid(0),sb_1_index(0),sblowindex(0),sbhighindex(0),
SB_2_gridname(0),SB_2_writegrid(0),sb_2_index(0),
THROOTwritegrid(0),throotindex(0),
ROOT(0),ROOTstat(0),ROOTstatcode(0),ROOTstatfilename(0),ROOTgridname(0),ROOTwritegrid(0),rootindex(0),zoneROOT(0),
threduindex(0),irrigindex(0),zoneFZT_QD(0),zoneFZT_QI(0),
code(0),flussdichte(0),SR_Anteil(0),kelsqi(0),
QD(0),QDd(0),QI(0),QIi(0),Qg(0),kelsqd(0),QB(0),QBb(0),kelsqlat(0),Q0(0),
zoneQD(0),zoneQI(0),zoneQB(0),zoneQg(0),zoneGWIN(0),zoneGWEX(0),zoneQDRAIN(0),zoneGWSTAND(0),zoneGWN(0),zonePOND(0),
zoneSB_1_(0),zoneSB_2_(0),zoneSB_2_IMMO(0),zoneINFEX(0),methode(0),levels(0),
INFEXstatcode(0),INFEXstatfilename(0),INFEXstat(0),INFEXgridname(0),SATTgridname(0),INFEXwritegrid(0),SATTwritegrid(0),
infexindex(0),sattindex(0),SATT(0),ini_gw_level(0),
GWTableExtern(NULL),m_bCoupled(false),FILLPOR(NULL),EXTBALANCE(NULL),
MACRO(0),MACROstat(0),MACROstatfilename(0),MACROstatcode(0),
IRRIGgridname(0),IRRIGwritegrid(0),MACROgridname(0),MACROwritegrid(0),macroindex(0),zoneMACRO(0),
INFILTstat(0),INFILTstatfilename(0),INFILTstatcode(0),zoneINFILT(0),
SUBSTEPS(0),SUBSTEPSstat(0),SUBSTEPSstatfilename(0),SUBSTEPSstatcode(0),
SUBSTEPSgridname(0),SUBSTEPSwritegrid(0),substepsindex(0),zoneSUBSTEPS(0),
CONC_QD(0),CONC_QI(0),conc_qdindex(0),conc_qiindex(0),
dParameterMinStep(0),lake_active(0),LAKE(0),lake_ini_method(0),m_useSurfaceRouting(false),S0(0),
reduced_ksat(1e-09),
do_TMod(0),iTModMethod(0),dTModLowBoundBase(12.0),dTModLapseRate(-0.007),dTModCDrySoil(800),dTModRhoDrySoil(1500),dTModLamDrySoil(0.58),dTModKmin(1e-10),
dTModLamWat(0.5562),dTModLamIce(2.33),dTModCWater(4187),dTModCIce20(1940),dTModCIce0(2090),dTModLatHeat(334000),dTModRhoWater(1000),dTModFactor(1.22),dTModSEForStat(0.5),
dTModShortestStep(60.0),dTModLongestStep(1800.0),dTModNFactorThawing(1),dTModNFactorFreezing(1),
temperbuffer(0),temperbufmax(0),temperbufcnt(0),temperfile(0),temperfilename(0),temper(0),SEstack(0),T(0),TLOWBND(0),tlowbndindex(0),temperindex(0),
thawdepth_stat(0),thawdepthstatcode(0),thawdepth_statfilename(0),temper_stackname(0),se_stackname(0),temper_writestack(0),zonethawdepth(0),do_initializeTemper(0),KX(0),KY(0),FROZENDEPTH(0),
THAWTMOD(0),thawdepth_gridname(0),thawdepth_writegrid(0),useHeatFluxAsLowerBC(0),max_Error(1.0e-5),max_ziter(100),
FDD(0),TDD(0),TDAY(0),TYEAR(0),numIntPreSpinUp(0),numIntPreSpinUpAll(0),fddindex(0),tddindex(0),tdayindex(0),tyearindex(0), 
interceptionBeforeSnow(false), dumpmode(false), maxDumpFiles(0), actDumpFiles(0), closeAfterDump(false), repairAfterDump(false),
headinitialized(false), HEAD1(NULL), mpiinitialized(false),
do_LayeredSnow(0),iLayeredSnowMethod(0),dLayeredSnowMinDensity(0),dLayeredSnowMaxDensity(0),
swes(NULL),swes1(NULL),swes_stackname(NULL),swes_writestack(0),swesindex(0),
swel(NULL),swel1(NULL),swel_stackname(NULL),swel_writestack(0),swelindex(0),
temperSnow(NULL), temperSnow1(NULL), temperSnow_stackname(NULL), temperSnow_writestack(0),temperSnowindex(0),
thetaSnow(NULL), thetaSnow1(NULL), thetaSnow_stackname(NULL), thetaSnow_writestack(0),thetaSnowindex(0),
hhydSnow(NULL), hhydSnow1(NULL), hhydSnow_stackname(NULL), hhydSnow_writestack(0),hhydSnowindex(0),
hgeoSnow(NULL), hgeoSnow1(NULL), hgeoSnow_stackname(NULL),hgeoSnow_writestack(0),hgeoSnowindex(0),
quSnow(NULL), quSnow1(NULL), quSnow_stackname(NULL), quSnow_writestack(0), quSnowindex(0),
thickSnow(NULL), thickSnow1(NULL), thickSnow_stackname(NULL), thickSnow_writestack(0), thickSnowindex(0),
densSnow(NULL), densSnow1(NULL), densSnow_stackname(NULL), densSnow_writestack(0), densSnowindex(0),
upperIndexSnow(NULL), upperIndexSnow_name(NULL), upperIndexSnow_writegrid(0),upperIndexSnowindex(0),
LayeredSnowProperties(NULL), SURFTEMP(NULL), EBAL(NULL), SURFTEMPGL(NULL), EBALGL(NULL), P_RAIN(NULL), SSNO(NULL), SLIQ(NULL), SSTO(NULL), QSNOW_OUT(NULL), QLOWBND(NULL),
levelsSnow(0), do_initializeTemperSnow(false),dTOffsetFreshSnow(0),dDrainageDensitiySnow(0),dWaterHoldingCapacity(0.1),dSEeffectiveness(0),dSEeffectivenessSoil(1),
GW1COLMATION(0), GW1COLMFROZEN(0), GW1COLMUNFROZEN(0),
RivLnkCellsToSend(p_RivLnkCellsToSend), RivLnkCellsToReceive(p_RivLnkCellsToReceive), RiverLinkVector(p_RiverLinkVector), 
numAllocatedDoubles(0), SendBuffer(NULL), RecvBuffer(0), SendBufferSize(NULL), RecvBufferSize(NULL), RecvBufferPos(NULL),
tempGridRows(NULL), tempStackRows(NULL), numTempGridRows(0), numTempStackRows(0),
do_ArcVeg(0), bArcVegActive(false), bArcVegWaiting(false), bArcVegSnowFree(false),
ExcessIceFactor(NULL), ExcessIceFactor1(NULL), ExcessIceFactor_stackname(NULL), ExcessIceFactor_writestack(0),ExcessIceFactorindex(0),
AVGDD(0),AVGDD_name(0),AVGDD_writegrid(0),AVGDDindex(0),
AVINMAX7DAVG(NULL),AVINMAX7DAVG_name(NULL),AVINMAX7DAVG_writegrid(0),AVINMAX7DAVGindex(0),
AVINHLP(NULL),AVINHLP1(NULL),AVINHLP_stackname(NULL),AVINHLP_writestack(0),AVINHLPindex(0),
AVTDAVG(NULL),AVTDAVG_name(NULL),AVTDAVG_writegrid(0),AVTDAVGindex(0),
AVINDAVG(NULL),AVINDAVG_name(NULL),AVINDAVG_writegrid(0),AVINDAVGindex(0),
AVSNOWFREE(NULL), AVSNOWFREE_name(NULL), AVSNOWFREE_writegrid(0), AVSNOWFREEindex(0),
AVPGP(NULL), AVPGP_name(NULL), AVPGP_writegrid(0), AVPGPindex(0),
AVINSCAL(NULL), AVINSCAL1(NULL), AVINSCAL_stackname(NULL), AVINSCAL_writestack(0), AVINSCALindex(0),
ArcVegStartJD(0), ArcVegEndJD(0), ArcVegWaitPeriod(0),
AVSON(NULL),AVSON_name(NULL),AVSON_writegrid(0),AVSONindex(0),
AVNMIN(NULL),AVNMIN1(NULL),AVNMIN_stackname(NULL),AVNMIN_writestack(0),AVNMINindex(0),
AVNMINTOT(NULL),AVNMINTOT_name(NULL),AVNMINTOT_writegrid(0),AVNMINTOTindex(0),
ARCVEGFRACTIONS(0),ARCVEGBIOMASS(0),ARCVEGRETRANS(0),
AVNumClimZon(0),AVNumPGPs(0),AVActClimZon(0),
AVLUSEType(0),AVLayerCnt(0),AVSoilDryDensity(0),
ICEWEDGELOWERDEPTH(NULL),
ICEWEDGEUPPERDEPTH(NULL),
indexFCmin_mm(0),indexFCmax_mm(0),indexFCact_mm(0),indexFCact_perc(0),
FCMIN(NULL), FCMAX(NULL), FCACT(NULL), FCACTPERC(NULL),doFCcalc(false),
sb_1_rel_index(0),SB_1_REL(NULL),doSB1RELcalc(false)
{
	try
	{
		retvalue=0;
		int     myindex,m,lv_temp,gw_depth_ini,s;
		unsigned long bufrow;
		double cellcount,pi_180;
		double th_start_sat;
		double *dicken; // alias for ST[s].getDzui()

		fpeflag=0;

		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif

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

		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow; // sr = start row for mpi with halo
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow; // er = end row for mpi with halo

		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}
		LUSE = new double**[iNumVegetationLayers+1];
		ETP = new double**[iNumVegetationLayers+1];
		ETR = new double**[iNumVegetationLayers+1];
		VBG = new double**[iNumVegetationLayers+1];
		WURZ = new double**[iNumVegetationLayers+1];
		// get the addresses of some parameter grids for output in otherdata only
		LAI = new double**[iNumVegetationLayers+1];
		RSC = new double**[iNumVegetationLayers+1];
		RSI = new double**[iNumVegetationLayers+1];
		Z0  = new double**[iNumVegetationLayers+1];
		// ALBEDO and RSE are single grids only but needed as well
		

		irrigationtable=p_irrigationtable;
//		double DT = InV*60.0;           /* interval in secounds (InV: interval in minutes) */
		double timestep = 6.0e4*InV/substepcnt;    /* factor for changing the measure of fluxes from m/s to mm/timestep */
		ST = SOILtable;

		riverflow_table=NULL;

		zonenidentification(&zn,&zone,&zarray,zonmax,zonmaxwert);

		/* allocating memory for statistik files */
		QDstatfilename      = newstring(_stringsize_+1,1);
		QIstatfilename      = newstring(_stringsize_+1,1);
		QBstatfilename      = newstring(_stringsize_+1,1);
		GWINstatfilename    = newstring(_stringsize_+1,1);
		GWEXstatfilename    = newstring(_stringsize_+1,1);
		Qgstatfilename      = newstring(_stringsize_+1,1);
		SB_1_statfilename   = newstring(_stringsize_+1,1);
		SB_2_statfilename   = newstring(_stringsize_+1,1);
		ROOTstatfilename    = newstring(_stringsize_+1,1);
		QDRAINstatfilename  = newstring(_stringsize_+1,1);
		GWSTANDstatfilename = newstring(_stringsize_+1,1);
		GWNstatfilename     = newstring(_stringsize_+1,1);
		PONDstatfilename    = newstring(_stringsize_+1,1);
		INFEXstatfilename   = newstring(_stringsize_+1,1);
		INFILTstatfilename  = newstring(_stringsize_+1,1);
		MACROstatfilename   = newstring(_stringsize_+1,1);
#ifdef _WASIMDEBUG_
		iteration_statfilename = newstring(_stringsize_+1,1);
#endif
		/* allocating memory for grids and stacks (stacks are stacks of grids) */
		theta_stackname   = newstring(_stringsize_+1,1);
		hh_stackname      = newstring(_stringsize_+1,1);
		hgeo_stackname    = newstring(_stringsize_+1,1);
		qu_stackname      = newstring(_stringsize_+1,1);
		GWSTANDgridname   = newstring(_stringsize_+1,1);
		GWthetagridname   = newstring(_stringsize_+1,1);
		GWNgridname       = newstring(_stringsize_+1,1);
		UPRISEgridname    = newstring(_stringsize_+1,1);
		PERCOLgridname    = newstring(_stringsize_+1,1);
		PONDgridname      = newstring(_stringsize_+1,1);
		GWLEVELgridname   = newstring(_stringsize_+1,1);
		QDRAINgridname    = newstring(_stringsize_+1,1);
		QDgridname        = newstring(_stringsize_+1,1);
		QIgridname        = newstring(_stringsize_+1,1);
		DELTAGWname       = newstring(_stringsize_+1,1);
		GWINgridname      = newstring(_stringsize_+1,1);
		GWEXgridname      = newstring(_stringsize_+1,1);
		INFEXgridname     = newstring(_stringsize_+1,1);
		SATTgridname      = newstring(_stringsize_+1,1);
		SB_1_gridname     = newstring(_stringsize_+1,1);
		SB_2_gridname     = newstring(_stringsize_+1,1);
		ROOTgridname      = newstring(_stringsize_+1,1);
		MACROgridname     = newstring(_stringsize_+1,1);
		IRRIGgridname     = newstring(_stringsize_+1,1);

		/* allocating memory for the names of poit-output-files (outputs for a single grid cell, but for all layers) */
		thetadateiname    = (char*)calloc(_stringsize_,sizeof(char));
		hhyddateiname     = (char*)calloc(_stringsize_,sizeof(char));
		otherdatadateiname = (char*)calloc(_stringsize_,sizeof(char));
		qudateiname       = (char*)calloc(_stringsize_,sizeof(char));
		etrdetaildateiname       = (char*)calloc(_stringsize_,sizeof(char));
		interflowdateiname       = (char*)calloc(_stringsize_,sizeof(char));
		/* filename for the textfile conataining the actual contents of all single linear reservoirs and of all flow travel zones at model end */
		storagefilename =(char*)calloc(_stringsize_,sizeof(char));

		/* reading parameters */
		fscanfln(steuerdatei,"%d",&methode            ," methode of "                            ,tstr);

		fscanfln(steuerdatei,"%d",&surface_interaction," surface water interaction parameter of ",tstr);
#ifdef _MPI
		int mpi_size_local;
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_local);
		if ((surface_interaction>1) && (mpi_size_local > 1))
		{ // is possible from version 10.03.00 on to use re-ifiltration in the MPI version with more than 1 rank. River link cells will be processed in flow order 
//			(*WSstream) << "FATAL: in modul " << bezeichnung << ":\n"
//				<< "parameter for surface_interaction (actual value: " << surface_interaction << ") can only be 0 or 1 for MPI version when running more than one instance!" << "[end-entry]";
//			throw(9);
		}
#endif

		fscanfln(steuerdatei,"%d",&pond_use           ," pond parameter of "                     ,tstr);
		fscanfln(steuerdatei,"%d",&drain_use          ," drainage parameter of "                 ,tstr);
		fscanfln(steuerdatei,"%d",&clay_use           ," clay parameter "                        ,tstr);
		fscanfln(steuerdatei,"%lf",&k_clay            ," k_sat. for clay  "                      ,tstr);
		fscanfln(steuerdatei,"%lf",&ini_gw_level      ," ini_gw_level  "                         ,tstr);

		if (methode>2)
		{
			SUBSTEPSstatfilename = newstring(_stringsize_+1,1);
			SUBSTEPSgridname  = newstring(_stringsize_+1,1);
		}


		/* name of statistic files */
		fscanfln(steuerdatei,"%s%d",QDRAINstatfilename, &QDRAINstatcode  ,"name of the QDRAINstatfiles in "  ,tstr); /* drainage discharge mm/interval */
		fscanfln(steuerdatei,"%s%d",GWSTANDstatfilename,&GWSTANDstatcode ,"name of the GWSTANDstatfiles in " ,tstr); /* groundwater level */
		fscanfln(steuerdatei,"%s%d",GWNstatfilename,    &GWNstatcode     ,"name of the GWNstatfiles in "     ,tstr); /* groundwater recharge mm/interval */
		fscanfln(steuerdatei,"%s%d",SB_1_statfilename,  &SB_1_statcode   ,"name of the SB_1_statfiles in "   ,tstr); /* rel. soilmoisture within root zone */
		fscanfln(steuerdatei,"%s%d",SB_2_statfilename,  &SB_2_statcode   ,"name of the SB_2_statfiles in "   ,tstr); /* abs. soilmoisture within unsaturated and saturated zone */
		fscanfln(steuerdatei,"%s%d",ROOTstatfilename,   &ROOTstatcode    ,"name of the ROOTstatfiles in "    ,tstr);  /* Root depth */
		fscanfln(steuerdatei,"%s%d",INFEXstatfilename,  &INFEXstatcode   ,"name of the INFEXstatfiles in "   ,tstr); /* infitration excess mm/interval */
		fscanfln(steuerdatei,"%s%d",PONDstatfilename,   &PONDstatcode    ,"name of the PONDstatfiles in "    ,tstr); /* contents of ponds mm */
		fscanfln(steuerdatei,"%s%d",QDstatfilename,     &QDstatcode      ,"name of the QDstatfiles in "      ,tstr); /* surface runoff mm/interval */
		fscanfln(steuerdatei,"%s%d",QIstatfilename,     &QIstatcode      ,"name of the QIstatfiles in "      ,tstr); /* interflow mm/intervall */
		fscanfln(steuerdatei,"%s%d",QBstatfilename,     &QBstatcode      ,"name of the QBstatfiles in "      ,tstr); /* baseflow mm/interval */
		fscanfln(steuerdatei,"%s%d",Qgstatfilename,     &Qgstatcode      ,"name of the Qgstatfiles in "      ,tstr); /* total flow mm/interval */
		fscanfln(steuerdatei,"%s%d",GWINstatfilename,   &GWINstatcode    ,"name of the GWINstatfiles in "    ,tstr); /* infiltration from rivers mm/interval */
		fscanfln(steuerdatei,"%s%d",GWEXstatfilename,   &GWEXstatcode    ,"name of the GWEXstatfiles in "    ,tstr); /* exfiltration into rivers (baseflow) mm/interval */
		fscanfln(steuerdatei,"%s%d",MACROstatfilename,  &MACROstatcode   ,"name of the MACROstatfiles in "   ,tstr); /* exfiltration into rivers (baseflow) mm/interval */
		fscanfln(steuerdatei,"%s%d",INFILTstatfilename, &INFILTstatcode  ,"name of the INFILTstatfiles in "  ,tstr); /* exfiltration into rivers (baseflow) mm/interval */

#ifdef _WASIMDEBUG_
		sprintf(iteration_statfilename,"iterate.log");iteration_statcode=2001;
#endif

		/* reading names for grids and/or stacks */
		fscanfln(steuerdatei,"%s",SB_1_gridname    ,"name of the SB05-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&SB_1_writegrid ,"Writegrid SB05-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",SB_2_gridname    ,"name of the SB1_-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&SB_2_writegrid ,"Writegrid SB1_-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",ROOTgridname    ,"name of the WURZ-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&ROOTwritegrid  ,"Writegrid WURZ-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",theta_stackname    ,"name of the Theta_stacks in "   ,tstr);
		fscanfln(steuerdatei,"%d",&theta_writestack ,"Writecode Theta in "           ,tstr);
		fscanfln(steuerdatei,"%s",hh_stackname       ,"name of the hh-stacks in "      ,tstr);
		fscanfln(steuerdatei,"%d",&hh_writestack    ,"Writecode hh in "              ,tstr);
		fscanfln(steuerdatei,"%s",hgeo_stackname     ,"name of the hgeo-stacks in "    ,tstr);
		fscanfln(steuerdatei,"%d",&hgeo_writestack  ,"Writecode hgeo  in "           ,tstr);
		fscanfln(steuerdatei,"%s",qu_stackname       ,"name of the q_unten-stacks in " ,tstr);
		fscanfln(steuerdatei,"%d",&qu_writestack    ,"Writecode qunten in "          ,tstr);

		fscanfln(steuerdatei,"%s",GWSTANDgridname    ,"name of the GWSTAND-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&GWSTANDwritegrid ,"Writegrid GWSTAND-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",GWthetagridname    ,"name of the GWtheta-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&GWthetawritegrid ,"Writegrid GWtheta-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",GWNgridname        ,"name of the GWN-Grids in "      ,tstr);
		fscanfln(steuerdatei,"%d",&GWNwritegrid     ,"Writegrid GWN-Grid in   "      ,tstr);
		fscanfln(steuerdatei,"%s",GWLEVELgridname    ,"name of the GWLEVEL-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&GWLEVELwritegrid ,"Writegrid GWLEVEL-Grid in   "  ,tstr);
		fscanfln(steuerdatei,"%s",QDRAINgridname     ,"name of the QDRAIN-Grids in "   ,tstr);
		fscanfln(steuerdatei,"%d",&QDRAINwritegrid  ,"Writegrid QDRAIN-Grid in   "   ,tstr);
		fscanfln(steuerdatei,"%s",SATTgridname       ,"name of the saturation grid in "          ,tstr);
		fscanfln(steuerdatei,"%d",&SATTwritegrid    ,"Writegrid Satt-Merker in "                        ,tstr);
		fscanfln(steuerdatei,"%s",INFEXgridname      ,"name of the infiltration excess grid in "   ,tstr);
		fscanfln(steuerdatei,"%d",&INFEXwritegrid   ,"Writegrid INFEX in "                              ,tstr);

		fscanfln(steuerdatei,"%s",QDgridname         ,"name of the QD-Grids in "           ,tstr);
		fscanfln(steuerdatei,"%d",&QDwritegrid      ,"Writegrid QD-Grid in   "           ,tstr);
		fscanfln(steuerdatei,"%s",QIgridname         ,"name of the QI-Grids in "           ,tstr);
		fscanfln(steuerdatei,"%d",&QIwritegrid      ,"Writegrid QI-Grid in   "           ,tstr);
		fscanfln(steuerdatei,"%s",DELTAGWname        ,"name of the DELTAGW-Grids in "      ,tstr);
		fscanfln(steuerdatei,"%d",&DELTAGWwritegrid ,"Writegrid DELTAGW-Grid in   "      ,tstr);
		fscanfln(steuerdatei,"%s",GWINgridname       ,"name of the GWIN-Grids in "         ,tstr);
		fscanfln(steuerdatei,"%d",&GWINwritegrid    ,"Writegrid GWIN-Grid in   "         ,tstr);
		fscanfln(steuerdatei,"%s",GWEXgridname       ,"name of the GWEX-Grids in "         ,tstr);
		fscanfln(steuerdatei,"%d",&GWEXwritegrid    ,"Writegrid GWEX-Grid in   "         ,tstr);
		fscanfln(steuerdatei,"%s",PONDgridname       ,"name of the POND-Grids in "         ,tstr);
		fscanfln(steuerdatei,"%d",&PONDwritegrid    ,"Writegrid POND-Grid in   "         ,tstr);

		fscanfln(steuerdatei,"%s",UPRISEgridname     ,"name of the GWN-Uprise-Grids in "   ,tstr);
		fscanfln(steuerdatei,"%d",&UPRISEwritegrid  ,"Writegrid GWN-Uprise-Grid in   "   ,tstr);
		fscanfln(steuerdatei,"%s",PERCOLgridname     ,"name of the GWN-PERCOL-Grids in "   ,tstr);
		fscanfln(steuerdatei,"%d",&PERCOLwritegrid  ,"Writegrid GWN-PERCOL-Grid in   "   ,tstr);
		fscanfln(steuerdatei,"%s",MACROgridname      ,"name of MACROPORE-Infiltration-Grids for "   ,tstr);
		fscanfln(steuerdatei,"%d",&MACROwritegrid  ,"Writegrid MACROPORE-Infiltration-Grid in   "   ,tstr);
		fscanfln(steuerdatei,"%s",IRRIGgridname    ,"name of the IRRIG-Grids in "  ,tstr);
		fscanfln(steuerdatei,"%d",&IRRIGwritegrid ,"Writegrid IRRIG-Grid in   "  ,tstr);

		/* here the row/column indices of a control location is given.
For this location, the three following files are generated */
		fscanfln(steuerdatei,"%d%d",&control_row,&control_col   ,"control_row/column for control point in ",tstr);
		fscanfln(steuerdatei,"%s",qudateiname        ,"Name of the Qbottom-file in "         ,tstr);
		fscanfln(steuerdatei,"%s",thetadateiname     ,"Name of the Theta-file in "          ,tstr);
		fscanfln(steuerdatei,"%s",hhyddateiname      ,"Name of the hhyd-file in "           ,tstr);
		fscanfln(steuerdatei,"%s",otherdatadateiname ,"Name of the Otherdata file in "      ,tstr);
		fscanfln(steuerdatei,"%s",etrdetaildateiname ,"Name of the etr-detail file in "     ,tstr);
		fscanfln(steuerdatei,"%s",interflowdateiname ,"Name of the interflow file in "      ,tstr);

		ExtendFileName(DefaultOutputDirectory,qudateiname);
		ExtendFileName(DefaultOutputDirectory,thetadateiname);
		ExtendFileName(DefaultOutputDirectory,hhyddateiname);
		ExtendFileName(DefaultOutputDirectory,otherdatadateiname);
		ExtendFileName(DefaultOutputDirectory,etrdetaildateiname);
		ExtendFileName(DefaultOutputDirectory,interflowdateiname);

		/* row/column-indizes are wrong */
		if ((control_row<0) || (control_row>nrows))
		{
			(*WSstream) << "FATAL: in modul "<<bezeichnung<<":\n"
			<< "row index ("<<control_row<<") for output-cell must be between "<<1<<" and "<<nrows<<" (Input as row col) or 0 0 "<<"[end-entry]";
			throw(8);
		}
		if ((control_col<0)||(control_col>ncols))
		{
			(*WSstream) << "FATAL: in modul "<<bezeichnung<<":\n"
			<< "column index ("<<control_col<<") for output-cell must be between "<<1<<" and "<<ncols<<" (Input as row col) "<<"[end-entry]";
			throw(9);
		}

		if ((methode<2)||(methode>3))
		{
			(*WSstream) << "FATAL: in modul "<<bezeichnung<<":\n"
			<< "invalid method ("<<methode<<")\n"
			<< "possible values are:\n 2 for Finite Differences with fixed timestep\n 3 for Finite Differences with dynamically calculated timestep\n"<<"[end-entry]";
			throw(91);
		}


		/* allocating memory for parameters (arrays with as many elements as zones) */
		code      = (double*)calloc(zonmax+1,sizeof(double)); /* codes of the zones corresponding to the zone grid */
		kelsqd    = (double*)calloc(zonmax+1,sizeof(double)); /* storage constant for QD */
		kelsqi    = (double*)calloc(zonmax+1,sizeof(double)); /* storage constant for interflow */
		flussdichte = (double*)calloc(zonmax+1,sizeof(double)); /* channel/stream density for interflow */
		//    QB        = (double*)calloc(zonmax+1,sizeof(double)); /* baseflow before SLS*/
		//    QBb       = (double*)calloc(zonmax+1,sizeof(double)); /* baseflow after SLS -> zonalmean */
		kelsqlat    = (double*)calloc(zonmax+1,sizeof(double)); /* storage constant for baseflow */
		/* !!!!! Qg, QD,QDd,QI,QIi are alloctaed in substance_transport_ini() because the numbers of tracers is unknown yet !!!!! */
		
		Qg        = NULL;
		QB		  = NULL;
		QBb       = NULL;
		QD        = NULL;
		QDd       = NULL;
		QI        = NULL;
		QIi       = NULL;

		Q0        = (double*)calloc(zonmax+1,sizeof(double)); /* maximum baseflow */
		SR_Anteil = (double*)calloc(zonmax+1,sizeof(double)); /* fraction of snowmelt, running off as surface runoff without infiltrating into the soil*/

		int kzone=(int)(zone[1].zonenwert); /* first code is minmal for the lowest basin number */
		do /* reading ranges for subbasin codes and parameters */
		{
			if (fscanf(steuerdatei,"%d",&myindex)!=1)
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading code[" << (int)myindex << "] for " << instanzname << "\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
				throw(10);
			}

			if ((myindex>zonmaxwert) || (myindex<1) || (zarray[myindex] == 0))
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error! there is no code " << myindex << " (from " << steuerdateiname
				<<") for a subbasin in "<<instanzname<<"\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
				throw(10);
			}
			else
				code[zarray[myindex]] = myindex;

			if (zarray[myindex]>zarray[kzone] + 1) /* filling codes inbetween with -1 */
				if ((kzone>zone[1].zonenwert) || ((kzone == zone[1].zonenwert) && (code[1]>0)))
					for (int j = zarray[kzone] + 1; j<zarray[myindex]; j++) code[j] = -1;
			else
			{
				code[zarray[kzone]] = zone[1].zonenwert;
				for (int j = zarray[kzone] + 1; j <= zarray[myindex]; j++) code[j] = -1; /* moving the first valid subbasin code to the first subbasin, if first code was higer */
			}
			kzone = myindex;
		}
		while (myindex<zonmaxwert);

		tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* reading storage constants for surface runoff */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(kelsqd[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading kelsqd["<<(int)code[i]<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
				/* Umrechnen kelsqd auf min */
				kelsqd[i]= kelsqd[i] >= 0.001 ? 60.0*kelsqd[i] : .06;
			} else kelsqd[i]=kelsqd[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* reading storage constants for interflow */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(kelsqi[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading kelsqi["<<(int)code[i]<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
				/* Umrechnen kelsqi auf min */
				kelsqi[i]= kelsqi[i] >= 0.001 ? 60.0*kelsqi[i] : .06;
			} else kelsqi[i]=kelsqi[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* reading stream density for interflow */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(flussdichte[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading channel network density["<<(int)code[i]
					<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
				//      flussdichte[i]*=cellsize/(double)1000.0; /* Umrechnen der Flussdichte von km/km in Flussdichte pro Zellen-m */
			} else flussdichte[i]=flussdichte[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* reading storage constants baseflow */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(kelsqlat[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading kelsqlat["<<(int)code[i]
					<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
				kelsqlat[i]= kelsqlat[i] >= 0 ? 1.0/kelsqlat[i] : 0.0;
			} else kelsqlat[i]=kelsqlat[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		/* reading maximum base flow */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(Q0[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading QB-Red.Faktor["<<(int)code[i]
					<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!\n"<<"[end-entry]";
					throw(10);
				}
			} else Q0[i]=Q0[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		/*  reading fraction of snow melt running off on the surface */
		for (int i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(SR_Anteil[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading SnowRunoff_fraction["<<(int)code[i]<<"] for "<<instanzname<<"\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
			} else SR_Anteil[i]=SR_Anteil[i-1];
		}  tstr=fgets(tstr,_stringsize_,steuerdatei);

		for (int i=1;i<=zonmax;i++) code[i]=zone[i].zonenwert;

		// reading the value of readgrids from control file
		if (fscanfln(steuerdatei,"%d",&gridslesen   ,"read control parameter for storage grids for ",tstr)!=1)
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
			"error reading readgrid-Parameter for "<<instanzname<<"\n"
			"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
			throw(10);
		}

		// file name of the storage file containing all content of the single linear storages for all flow time zones 
		if (fscanfln(steuerdatei,"%s",storagefilename,"name of the storage content file ",tstr)!=1)
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
			"error reading storage file name for "<<instanzname<<"\n"
			"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
			throw(10);
		}
        char tsn[1000];
		sprintf(tsn, "%s_%d", storagefilename, mpi_rank);
//		sprintf(tsn, "%s_%d", storagefilename, 0);
		sprintf(storagefilename, "%s", tsn);

		// file name of the storage file containing all content of the single linear storages for all flow time zones 
		if (methode>2)
		{
			if (fscanfln(steuerdatei,"%lf",&dParameterMinStep,"minimum dynamic time step in secounds ",tstr)!=1)
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error reading dParameterMinStep in "<<instanzname<<"\n"
				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
				throw(10);
			}
			fscanfln(steuerdatei,"%s%d",SUBSTEPSstatfilename,  &SUBSTEPSstatcode   ,"name of the SUBSTEPSstatfiles in "   ,tstr); /* exfiltration into rivers (baseflow) mm/interval */
			fscanfln(steuerdatei,"%s",SUBSTEPSgridname      ,"name of SUBSTEPS-Grids for "   ,tstr);
			fscanfln(steuerdatei,"%d",&SUBSTEPSwritegrid  ,"Writegrid SUBSTEPS-Grid in   "   ,tstr);
		}

		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: control file \""<<steuerdateiname<<"\" cannot be closed"<<"[end-entry]";
			throw(11);
		}

		// now reading the heat transfer parameters if this section exists and do_TMod = 1 gesetzt ist
		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		if (findeabschnitt(_heat_transfer_,false))
		{
			fscanfln(steuerdatei,"%d",&do_TMod,"",tstr);
			if (do_TMod>0)
			{
				int dummyi;
				fscanfln(steuerdatei,"%d",&iTModMethod," method of heat transfer modelling ",tstr);
				fscanfln(steuerdatei,"%d",&dummyi,"",tstr);
				fscanfln(steuerdatei,"%d",&dummyi,"",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLowBoundBase," dTModLowBoundBase ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLapseRate," dTModLapseRate ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModCDrySoil," dTModCDrySoil ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModRhoDrySoil," dTModRhoDrySoil ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLamDrySoil," dTModLamDrySoil ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModKmin," dTModKmin ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLamWat," dTModLamWat ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLamIce," dTModLamIce ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModCWater," dTModCWater ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModCIce20," dTModCIce20 ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModCIce0," dTModCIce0 ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLatHeat," dTModLatHeat ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModRhoWater," dTModRhoWater ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModFactor," dTModFactor ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModSEForStat," dTModSEForStat ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModShortestStep," dTModShortestStep ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModLongestStep," dTModLongestStep ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModNFactorFreezing," dTModNFactorFreezing ",tstr);
				fscanfln(steuerdatei,"%lf",&dTModNFactorThawing," dTModNFactorThawing ",tstr);
				
				thawdepth_statfilename   = newstring(_stringsize_+1,1);
				temper_stackname      = newstring(_stringsize_+1,1);
				se_stackname = newstring(_stringsize_ + 1, 1);
				temperfilename     = (char*)calloc(_stringsize_,sizeof(char));
				thawdepth_gridname = newstring(_stringsize_ + 1, 1);

				fscanfln(steuerdatei,"%s",temperfilename      ,"Name of the temperature-file for the control point location in "           ,tstr);
				fscanfln(steuerdatei,"%s%d",thawdepth_statfilename,  &thawdepthstatcode   ,"name of the temperstatfiles in "   ,tstr); 
				
				fscanfln(steuerdatei,"%s",temper_stackname       ,"name of the temper-stacks in "      ,tstr);
				fscanfln(steuerdatei,"%d",&temper_writestack    ,"Writecode temper in "              ,tstr);
				
				fscanfln(steuerdatei, "%s", thawdepth_gridname, "name of the thawdepth-grid (heat-tr.mod.) in ", tstr);
				fscanfln(steuerdatei, "%d", &thawdepth_writegrid, "Writecode thawdepth (heat-tr.mod.) in ", tstr);
				// --> new in 10.04.05, control file needs to be upgraded
				fscanfln(steuerdatei, "%s", se_stackname, "name of the SE-stack (heat-tr.mod.) in ", tstr);
				fscanfln(steuerdatei, "%d", &se_writestack, "Writecode SE-stack (heat-tr.mod.) in ", tstr);
				
				ExtendFileName(DefaultOutputDirectory,temperfilename);
				
			}
		}

		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: control file \""<<steuerdateiname<<"\" cannot be closed"<<"[end-entry]";
			throw(11);
		}



		// now reading the ARCVEG parameters if this section exists and do_TMod = 1 gesetzt ist
		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}
		if (findeabschnitt(_ArcVegIntegration_, false)) // looks for a section  [ArcVegIntegration]
		{
			fscanfln(steuerdatei, "%d", &do_ArcVeg, "", tstr);
			if (do_TMod == 0)
			{
				(*WSstream) << "FATAL: ArcVeg integration is only available if heat_transfer is switched on. " << "[end-entry]";
				throw(1501);
			}

			fscanfln(steuerdatei, "%d", &ArcVegStartJD, "ArcVeg start Julian Day in ", tstr);
			fscanfln(steuerdatei, "%d", &ArcVegEndJD, "ArcVeg end Julian Day in ", tstr);
			fscanfln(steuerdatei, "%d", &ArcVegWaitPeriod, "ArcVeg waiting period for inundation avg. after start day in ", tstr);
			fscanfln(steuerdatei, "%d%d", &AVNumClimZon, &AVNumPGPs,"ArcVeg count of climate zones and count of plant growing periods (PGPs) in ", tstr);
			
			if (AVNumClimZon <= 0 || AVNumPGPs <= 0)
			{
				(*WSstream) << "FATAL: ArcVeg climate zone count (AVNumClimZon) and PGP count (AVNumPGPs) must be both provided and be >= 0" << "[end-entry]";
				throw(1502);
			}

			
			AVNMinParam = new double*[AVNumClimZon+1];
			for (int k = 0; k <= AVNumClimZon; k++)
			{
				AVNMinParam[k] = new double[AVNumPGPs+1];
			}
			for (int i = 1; i <= AVNumClimZon; i++)
			{
				fgets(tstr, _stringsize_, steuerdatei);
				std::stringstream ss(tstr);
				for (int j = 1; j <= AVNumPGPs; j++)
				{
					ss >> AVNMinParam[i][j];
					if (!(ss.good()))
					{
						(*WSstream) << "FATAL: ArcVeg AVNMinParam for climate zone " << i << " and PGP " << j << " could not be read in correctly " << "[end-entry]";
						(*WSstream) << "This is the last line read in: " << tstr << "[end-entry]";
						throw(1503);
					}
				}
			}

			fscanfln(steuerdatei, "%d", &AVActClimZon, "ArcVeg actual climate zone (AVActClim) in ", tstr);
			if (AVActClimZon <= 0 || AVActClimZon > AVNumClimZon) 
			{
				(*WSstream) << "FATAL: ArvVeg actual climate zone \"" << AVActClimZon << "\" not in range of NMin params dimensions (" << AVNumClimZon << " climate zones and " << AVNumPGPs << " Plant Growing Periods PGPs)"<< "[end-entry]";
				throw(1504);
			}

			// land use type that will be used for changing the soil thickness (i.e. when the biomass of moss is changing, this soil type represents the moss)
			fscanfln(steuerdatei, "%d", &AVLUSEType, "ArcVeg land use type for changing soil thickness (AVLUSEType) in ", tstr);
			if (AVLUSEType <= 0 || AVLUSEType > LUtable[0].getCode())
			{
				(*WSstream) << "FATAL: ArvVeg land use type to be used for thickness change (AVLUSEType) \"" << AVLUSEType << "\" not in range of Landuse table entries (" << 1 << " to " << LUtable[0].getCode() << ")" << "[end-entry]";
				throw(1506);
			}

			// COunt of soil layers that should be changed in thickenss when arcveg runs
			fscanfln(steuerdatei, "%d", &AVLayerCnt, "ArcVeg Layer count of soil layers to be changed in thickness  (AVLayerCnt) in ", tstr);
			if (AVLayerCnt <= 0 || AVLUSEType > SOILtable[1].getLevels())
			{
				(*WSstream) << "FATAL: ArcVeg Layer count of soil layers to be changed in thickness  (AVLayerCnt) \"" << AVLayerCnt << "\" not in range of soil table layers (" << 1 << " to " << SOILtable[1].getLevels() << ")" << "[end-entry]";
				throw(1507);
			}


			// density of dry MOSS (or whatever land use type is used for thickness change)
			fscanfln(steuerdatei, "%lg", &AVSoilDryDensity, "ArcVeg dry density of soil layers to be changed in thickness  (AVSoilDryDensity) in ", tstr);
			if (AVSoilDryDensity <= 1 || AVLUSEType > 1000)
			{
				(*WSstream) << "FATAL: ArcVeg dry density of soil layers to be changed in thickness  (AVSoilDryDensity) \"" << AVSoilDryDensity << "\" not in a feasible range (" << 1.0 << " to " << 1000.0 << ")" << "[end-entry]";
				throw(1508);
			}

			ExcessIceFactor_stackname = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", ExcessIceFactor_stackname, "name of the ExcessIceFactor-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &ExcessIceFactor_writestack, "Writecode ExcessIceFactor-stack in ", tstr);
			
			AVGDD_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVGDD_name, "name of the ArcVeg-growing-degree-days grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVGDD_writegrid, "Writecode ArcVeg-growing-degree-days grid in ", tstr);

			AVINMAX7DAVG_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVINMAX7DAVG_name, "name of the ArcVeg-Inundation-Max-7-Day-Avg in ", tstr);
			fscanfln(steuerdatei, "%d", &AVINMAX7DAVG_writegrid, "Writecode ArcVeg-Inundation-Max-7-Day-Avg in ", tstr);

			AVINHLP_stackname = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVINHLP_stackname, "name of the ArcVeg-Inundation-7-Day-Avg-Helper-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &AVINHLP_writestack, "Writecode ArcVeg-Inundation-7-Day-Avg-Helper-stack in ", tstr);

			AVTDAVG_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVTDAVG_name, "name of the ArcVeg-Daily-Average-Temperature-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVTDAVG_writegrid, "Writecode ArcVeg-Daily-Average-Temperature-Grid in ", tstr);

			AVINDAVG_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVINDAVG_name, "name of the ArcVeg-Daily-Average-Inundation-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVINDAVG_writegrid, "Writecode ArcVeg-Daily-Average-Inundation-Grid in ", tstr);

			AVSNOWFREE_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVSNOWFREE_name, "name of the ArcVeg-SnowFreeDays-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVSNOWFREE_writegrid, "Writecode ArcVeg-SnowFreeDays-Grid in ", tstr);

			AVPGP_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVPGP_name, "name of the ArcVeg-AVPGP-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVPGP_writegrid, "Writecode ArcVeg-AVPGP-Grid in ", tstr);

			AVINSCAL_stackname = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVINSCAL_stackname, "name of the ArcVeg-AVINSCAL-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &AVINSCAL_writestack, "Writecode ArcVeg-AVINSCAL-stack in ", tstr);

			AVNMIN_stackname = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVNMIN_stackname, "name of the ArcVeg-AVNMIN-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &AVNMIN_writestack, "Writecode ArcVeg-AVNMIN-stack in ", tstr);

			AVNMINTOT_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVNMINTOT_name, "name of the ArcVeg-NMINTOT-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVNMINTOT_writegrid, "Writecode ArcVeg-NMINTOT-Grid in ", tstr);

			AVSON_name = newstring(_stringsize_ + 1, 1);
			fscanfln(steuerdatei, "%s", AVSON_name, "name of the ArcVeg-AVSON-Grid in ", tstr);
			fscanfln(steuerdatei, "%d", &AVSON_writegrid, "Writecode ArcVeg-AVSON-Grid in ", tstr);

		}

		if (fclose(steuerdatei) == EOF)
		{
			(*WSstream) << "FATAL: control file \"" << steuerdateiname << "\" cannot be closed" << "[end-entry]";
			throw(11);
		}



		// [unsatzon_tuning]
		// max_error = 1e-7 # maximum allowed error between the hydraulic heads of two susequent iterations (maximum pairwise difference of all layers)
		// max_iterations = 1000 # maximum number of iterations for Richards equation (will later be replaced by an implicit solution schema)
		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		strcpy(tstr,"[");
		strcat(tstr,_unsatzon_tuning_);
		strcat(tstr,"]");
	
		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		if (findeabschnitt(tstr,false))
		{
			bool DoContinue=true;
			do
			{
				fgets(tstr,_stringsize_,steuerdatei);
				std::string s(tstr);
				s=s.substr(s.find_first_not_of(" \t"));
				sprintf(tstr,"%s",s.c_str());
				if ((!feof(steuerdatei)) && (index(tstr, "\n")>0) && (index(tstr, "\r")>0)) // skip all empty lines
				{
					DoContinue=true;
					std::string parameter(tstr);
					if (parameter.find("max_error")!=std::string::npos) {
						std::stringstream ss;
						parameter=parameter.substr(parameter.find("=")+1);
						ss << parameter.c_str();
						ss >> max_Error;
					} else if (parameter.find("max_iterations")!=std::string::npos) {
						std::stringstream ss;
						parameter=parameter.substr(parameter.find("=")+1);
						ss << parameter.c_str();
						ss >> max_ziter;
					}
				}
				else 
				{
					DoContinue=false;
				}
			} while (DoContinue);

		}

		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: control file \""<<steuerdateiname<<"\" cannot be closed"<<"[end-entry]";
			throw(11);
		}

		// searching for the grid containing the thawing depth for permafrost soils.
		// if found, permafrost must be regarded at points where THAWDEPTH hat valid values
		// to be sure that the snow model correctly initialized also snowfreedays and thawcoefficient alpha, the other gridas must be found as well
		if (do_TMod==0 && ist_vorhanden(_PermafrostThawDepth_) && ist_vorhanden(_SnowFreeDaysGrid_) && ist_vorhanden(_ThawCoeffPermaFrost_))
		{
			THAWDEPTH=gridinfo[getindex(_PermafrostThawDepth_)].pointer;
			// getting reduced_ksat from control file
			if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
			{
				(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
				throw(28);
			}
			strcpy(tstr,"[");
			strcat(tstr,_permafrost_);
			strcat(tstr,"]");
	
			/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
			findeabschnitt(tstr);

			// overread first four lines
			int dummyi;
			double dummyf;
			char dummys[1000];
			fscanfln(steuerdatei,"%d",&dummyi,"",tstr);
			fscanfln(steuerdatei,"%lf",&dummyf,"",tstr);
			fscanfln(steuerdatei,"%lf",&dummyf,"",tstr);
			fscanfln(steuerdatei,"%s%d",dummys,&dummyi,"",tstr);
			fscanfln(steuerdatei,"%lf",&reduced_ksat,"reduced k_sat (in m/s)",tstr);
			
			if (fclose(steuerdatei)==EOF)
			{
				(*WSstream) << "FATAL: control file \""<<steuerdateiname<<"\" cannot be closed"<<"[end-entry]";
				throw(11);
			}
		}
		else
		{
			THAWDEPTH=NULL;
		}


		/* initializing a statistc register set for groundwater level [m below surface] */
		zoneGWSTAND = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneGWSTAND) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{
			zoneGWSTAND[i].n         = zoneGWSTAND[i].RelCellSum = zoneGWSTAND[i].sx        = 0;
			zoneGWSTAND[i].zonenwert = zone[i].zonenwert;
			zoneGWSTAND[i].rel       = zone[i].rel;
		}

		/* initializing a statistc register set for soil moisture [0...root depth] */
		zoneSB_1_ = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSB_1_) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneSB_1_[i].n         = zoneSB_1_[i].RelCellSum = zoneSB_1_[i].sx        = 0;
			zoneSB_1_[i].zonenwert = zone[i].zonenwert;
			zoneSB_1_[i].rel       = zone[i].rel;
		}

		/* initializing a statistc register set for soil moisturte [entire unsaturated and saturated zone] */
		zoneSB_2_ = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSB_2_) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneSB_2_[i].n         = zoneSB_2_[i].RelCellSum = zoneSB_2_[i].sx        = 0;
			zoneSB_2_[i].zonenwert = zone[i].zonenwert;
			zoneSB_2_[i].rel       = zone[i].rel;
		}

		/* initializing a statistc register set for root depth */
		zoneROOT = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneROOT) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneROOT[i].n         = zoneROOT[i].RelCellSum = zoneROOT[i].sx        = 0;
			zoneROOT[i].zonenwert = zone[i].zonenwert;
			zoneROOT[i].rel       = zone[i].rel;
		}

		/* initializing a statistc register set for infitration excess */
		zoneINFEX = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneINFEX) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneINFEX[i].n         = zoneINFEX[i].n         = zoneINFEX[i].sx        = 0;
			zoneINFEX[i].zonenwert = zone[i].zonenwert;
			zoneINFEX[i].rel       = zone[i].rel;
		}

		/*initializing a statistc register set for direct runoff */
		zoneQD = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQD) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneQD[i].n         = zoneQD[i].RelCellSum = zoneQD[i].sx        = 0;
			zoneQD[i].zonenwert = zone[i].zonenwert;
			zoneQD[i].rel       = zone[i].rel;
		}

		/*initializing a statistc register set for interflow */
		zoneQI = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQI) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneQI[i].n         = zoneQI[i].RelCellSum = zoneQI[i].sx        = 0;
			zoneQI[i].zonenwert = zone[i].zonenwert;
			zoneQI[i].rel       = zone[i].rel;
		}

		/*initializing a statistc register set for baseflow */
		zoneQB = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQB) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneQB[i].n         = zoneQB[i].RelCellSum = zoneQB[i].sx        = 0;
			zoneQB[i].zonenwert = zone[i].zonenwert;
			zoneQB[i].rel       = zone[i].rel;
		}

		/*initializing a statistc register set for total flow */
		zoneQg = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQg) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneQg[i].n         = zoneQg[i].RelCellSum = zoneQg[i].sx        = 0;
			zoneQg[i].zonenwert = zone[i].zonenwert;
			zoneQg[i].rel       = zone[i].rel;
		}


		/*initializing a statistc register set for groundwater recharge */
		zoneGWN = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneGWN) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneGWN[i].n         = zoneGWN[i].RelCellSum = zoneGWN[i].sx        = 0;
			zoneGWN[i].zonenwert = zone[i].zonenwert;
			zoneGWN[i].rel       = zone[i].rel;
		}

		/*initializing a statistc register set for macropore infiltration */
		zoneMACRO = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneMACRO) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneMACRO[i].n         = zoneMACRO[i].RelCellSum = zoneMACRO[i].sx        = 0;
			zoneMACRO[i].zonenwert = zone[i].zonenwert;
			zoneMACRO[i].rel       = zone[i].rel;
		}

		if (methode>2)
		{
			zoneSUBSTEPS = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneSUBSTEPS) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{ 
				zoneSUBSTEPS[i].n         = zoneSUBSTEPS[i].RelCellSum = zoneSUBSTEPS[i].sx        = 0;
				zoneSUBSTEPS[i].zonenwert = zone[i].zonenwert;
				zoneSUBSTEPS[i].rel       = zone[i].rel;
			}
		}

		/*initializing a statistc register set for infiltration */
		zoneINFILT = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneINFILT) speichererror(1);
		for (int i=0;i<=zonmax;i++)
		{ 
			zoneINFILT[i].n         = zoneINFILT[i].RelCellSum = zoneINFILT[i].sx        = 0;
			zoneINFILT[i].zonenwert = zone[i].zonenwert;
			zoneINFILT[i].rel       = zone[i].rel;
		}

		if (do_TMod>0)
		{
			/* initializing a statistc register set for infitration excess */
			zonethawdepth = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zonethawdepth) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{ 
				zonethawdepth[i].n         = zonethawdepth[i].n         = zonethawdepth[i].sx        = 0;
				zonethawdepth[i].zonenwert = zone[i].zonenwert;
				zonethawdepth[i].rel       = zone[i].rel;
			}

		}


		/* scanning the gridinfo-array for the adresses of some special grids */

		char name[500];
		int lay;
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;
			
			sprintf(name,"%s%d",_Wurzeltiefe_,lay);
			WURZ[lay]  =  gridinfo[getindex(name)].pointer;                       /* root depth */

			sprintf(name,"%s%d",_Vegetationsbedeckungsgrad_,lay);
			VBG[lay]   =  gridinfo[getindex(name)].pointer;         /* vegetation coverage degree */

			// the following entities will be used for output in otherdata only
			sprintf(name,"%s%d",_Blattflaechenindex_,lay);
			LAI[lay]   =  gridinfo[getindex(name)].pointer;         /* leaf are index */

			sprintf(name,"%s%d",_SurfaceCanopyResistance_,lay);
			RSC[lay]   =  gridinfo[getindex(name)].pointer;         /* canopy surface resistance */

			sprintf(name,"%s%d",_SurfaceInterzepResistance_,lay);
			RSI[lay]   =  gridinfo[getindex(name)].pointer;         /* interception surface resistance */

			sprintf(name,"%s%d",_RoughnessLength_,lay);
			Z0[lay]   =  gridinfo[getindex(name)].pointer;         /* roughess length */

		}
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Potentielle_Verdunstung_,lay);
			ETP[lay]    = gridinfo[getindex(name)].pointer;

			sprintf(name,"%s%d",_Reale_Verdunstung_,lay);
			ETR[lay]    = gridinfo[getindex(name)].pointer;

		}
		// the following entities will be used for output in otherdata only
		RSE = gridinfo[getindex(_SurfaceEvaporationResistance_)].pointer;
		ALBEDO = gridinfo[getindex(_Albedo_)].pointer;

		EVAP = gridinfo[getindex(_Potential_Evaporation_)].pointer;
		EVAR = gridinfo[getindex(_Reale_Evaporation_)].pointer;

		SLP   =  gridinfo[getindex(_Gefaelle_)].pointer;                          /* slope angle */
		SOIL  =  gridinfo[getindex(_Bodenarten_)].pointer;                        /* soil types (texture) */

		/* searching precipitation grid (input into soil model from other model components
-> if interception module was run: throughfall)
-> if not: snowmelt
-> if no snow model was run: precipitation*/

// is done in method setInterceptionBeforeSnow now because sequence of snow model and interception mode may be reversed
//		sprintf(name,"%s%d",_Bestandesniederschlag_,0);
//		P     =  gridinfo[getindex_first_of_3(name,_Schneeausfluss_,_Niederschlag_)].pointer;

		sprintf(name,"%s%d",_gw_storage_coeff_,1);
		if (ist_vorhanden(name)) // if no groundwtare model is used, S0 may not be there, then simply do not use it later
		{ 
			S0 = gridinfo[getindex(name)].pointer;  /* storage coefficients for layer 1 in groundwater module will be set in each interval dynamically */
		}

		sprintf(name,"%s%d",_gw_k_x_,1);
		if (ist_vorhanden(name)) // if no groundwtare model is used, S0 may not be there, then simply do not use it later
		{ 
			KX = gridinfo[getindex(name)].pointer;  /* storage coefficients for layer 1 in groundwater module will be set in each interval dynamically */
		}

		sprintf(name,"%s%d",_gw_k_y_,1);
		if (ist_vorhanden(name)) // if no groundwtare model is used, S0 may not be there, then simply do not use it later
		{ 
			KY = gridinfo[getindex(name)].pointer;  /* storage coefficients for layer 1 in groundwater module will be set in each interval dynamically */
		}

		if (pond_use==1)
		{
			MAXPOND = gridinfo[getindex(_max_ponding_storage_)].pointer;
			int index=getindex(_max_ponding_storage_);
			// make sure that in rank 0 all MAXPOND cells will be changes in a way that runoff routing gets mm instead of m
			if (mpi_rank == 0)
			{
				// entire grid is in memory for rank 0 except empty rows
				for (int j = 1; j <= nrows; j++)
				{
					for (int k = FirstLastCol[j].First; k <= FirstLastCol[j].Last; k++)
					{
						if (zn[j][k] != nodata)
						{
							MAXPOND[j][k] *= 1000.0; /* converting the pond depth unit from m into mm */
						}
					}
				}
			}
#ifdef _MPI
			MPI_Barrier(MPI_COMM_WORLD); // important since MAXPOND will be used later when calling initialize_lakes, so we need to make sure that 
#endif
			MPIScatterGrid(index);
#ifdef _MPI
			MPI_Barrier(MPI_COMM_WORLD); // important since MAXPOND will be used later when calling initialize_lakes, so we need to make sure that 
#endif
		}
		else
		{
			MAXPOND=NULL;
		}
		if (drain_use==1)
		{ /* drainage depth and densiti are searched for only if drainage is used */
			HDRAIN    = gridinfo[getindex(_Drainagehoehen_)].pointer;
			DRAINDIST = gridinfo[getindex(_Drainagedistance_)].pointer;
		}
		else
		{
			HDRAIN = NULL;
			DRAINDIST = NULL;
		}
		if (clay_use==1)
		CLAY = gridinfo[getindex(_clay_depth_)].pointer;
		else
		CLAY = NULL;

		/* SNOW is searched for seperately without ending the programm (may be, snowmodel was turned off) */
		int gridindex = 1;
		while (strcmp(gridinfo[gridindex].identification, _snowstorage_solid_) != 0)
		{ 
			gridindex++;
			if ((gridindex >= gridmaxanzahl) || (gridinfo[gridindex].identification == NULL)) break;
		}
		if ((gridindex >= gridmaxanzahl) || (gridinfo[gridindex].identification == NULL))
		{
			SNOW = NULL;
		}
		else
		{
			SNOW = gridinfo[gridindex].pointer;
		}

		/* flow travel time is searched for seperately without ending the programm (may be, there was no fzt-grid) */
		gridindex=1;
		while (strcmp(gridinfo[gridindex].identification, _Fliesszeiten_) != 0)
		{
			gridindex++;
			if ((gridindex >= gridmaxanzahl) || (gridinfo[gridindex].identification == NULL)) break;
		}
		if ((gridindex >= gridmaxanzahl) || (gridinfo[gridindex].identification == NULL))
		{
			FZT = NULL;
			(*WSstream) << "INFO: \n------------------------------------------------------------------------\n"
			"\nWarning! There was no flow travel time grid for surface runoff read in\n"
			"All surface runoof will be routed immediately to the outlet of the sub-basin\n"
			"(except the effects of a linear reservoir).\n"
			"A flow travel time grid may be created by TANALYS, writing the\n"
			"flow travel times for each grid cell to the sub basins outlet in secounds\n"
			"which are reclassified into model time steps in WaSiM.\n"
			"\n------------------------------------------------------------------------\n"<<"[end-entry]";
		}
		else  FZT = gridinfo[gridindex].pointer;

		if (do_TMod>0)
		{
			T = gridinfo[getindex(_Temperatur_)].pointer;
			if (ist_vorhanden(_T_Lower_Boundary_Condition_))
			{
				TLOWBND = gridinfo[getindex(_T_Lower_Boundary_Condition_)].pointer;
			}
			else
			{
				// create the lower boundary condition later, but initialize the index here with a specific value
				tlowbndindex=-99;
			}
			if (ist_vorhanden(_HeatFlux_Lower_Boundary_Condition_))
			{
				useHeatFluxAsLowerBC=1;
				QLOWBND = gridinfo[getindex(_HeatFlux_Lower_Boundary_Condition_)].pointer;
			}
			else
			{
				useHeatFluxAsLowerBC=0; // was already set by default during member initialization for unsatzonmodellklasse
				QLOWBND = NULL;
			}
		}

		/* constructing statistic files */
		GWSTANDstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneGWSTAND,(int)zonmax,(int)GWSTANDstatcode,(int)5,(char*)GWSTANDstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		GWSTANDstat->schreibzonenwerte("groundwater level [m below surface] (GWSTAND)",instanzname,LastResults,maxLastResults);
		GWNstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneGWN,(int)zonmax,(int)GWNstatcode,(int)5,(char*)GWNstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		GWNstat->schreibzonenwerte("groundwater recharge per zone [mm] (GWN)",instanzname,LastResults,maxLastResults);
		SB_1_stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSB_1_,(int)zonmax,(int)SB_1_statcode,(int)5,(char*)SB_1_statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		SB_1_stat->schreibzonenwerte("rel. soil moisture 0..root depth per zone ",instanzname,LastResults,maxLastResults);
		SB_2_stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSB_2_,(int)zonmax,(int)SB_2_statcode,(int)5,(char*)SB_2_statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		SB_2_stat->schreibzonenwerte("abs. soil moisture 0..profile depth per zone ",instanzname,LastResults,maxLastResults);
		ROOTstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneROOT,(int)zonmax,(int)ROOTstatcode,(int)5,(char*)ROOTstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		ROOTstat->schreibzonenwerte("root depth [mm per zone] ",instanzname,LastResults,maxLastResults);
		INFEXstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneINFEX,(int)zonmax,(int)INFEXstatcode,(int)5,(char*)INFEXstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		INFEXstat->schreibzonenwerte("Infiltration excess [mm per zone]  ",instanzname,LastResults,maxLastResults);
		QDstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQD,(int)zonmax,(int)QDstatcode,(int)6,(char*)QDstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		QDstat->schreibzonenwerte("direct discharge [mm per Zone] (QD)",instanzname,LastResults,maxLastResults);
		QIstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQI,(int)zonmax,(int)QIstatcode,(int)6,(char*)QIstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		QIstat->schreibzonenwerte("Interflow [mm per zone]",instanzname,LastResults,maxLastResults);
		QBstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQB,(int)zonmax,(int)QBstatcode,(int)6,(char*)QBstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		QBstat->schreibzonenwerte("base flow [mm per Zone]",instanzname,LastResults,maxLastResults);
		Qgstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQg,(int)zonmax,(int)Qgstatcode,(int)6,(char*)Qgstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		Qgstat->schreibzonenwerte("total runoff [mm per zone] (sum of QB, QI und QD)",instanzname,LastResults,maxLastResults);
		MACROstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneMACRO,(int)zonmax,(int)MACROstatcode,(int)6,(char*)MACROstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		MACROstat->schreibzonenwerte("macropore_infiltration_per_zone_[mm] ",instanzname,LastResults,maxLastResults);
		INFILTstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneINFILT,(int)zonmax,(int)INFILTstatcode,(int)6,(char*)INFILTstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		INFILTstat->schreibzonenwerte("total_infiltration_per_zone_[mm] ",instanzname,LastResults,maxLastResults);
		if (methode>2)
		{
			SUBSTEPSstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSUBSTEPS,(int)zonmax,(int)SUBSTEPSstatcode,(int)6,(char*)SUBSTEPSstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
			SUBSTEPSstat->schreibzonenwerte("Number_of_Substeps_per_Time_Step[-] ",instanzname,LastResults,maxLastResults);
		}

		if (do_TMod)
		{
			thawdepth_stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zonethawdepth,(int)zonmax,
							(int)thawdepthstatcode,(int)5,(char*)thawdepth_statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			thawdepth_stat->schreibzonenwerte("soil thaw depth [m] ",instanzname,LastResults,maxLastResults);
		}

#ifdef _WASIMDEBUG_
		iteration_stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zone,(int)zonmax,(int)iteration_statcode,(int)6,(char*)iteration_statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
		iteration_stat->schreibzonenwerte("mean number of iterations per subbasin ",instanzname,LastResults,maxLastResults);
#endif

		/* usually, the soil types from the soil table should have all the same number of levels.
now the soiltype with the most levels is searched for */
		levels=0;
		for (int i=1;i<=ST[0].getCode();i++) /* ST stands for soil table. In ST[0].code the number of soil types is stored */
		{
			if (ST[i].getLevels()>levels) levels=ST[i].getLevels(); /* finding the maximum number of levels for following checks of soil table validity */
			ST[0].setLevels(levels);
		}

		/* now some new grids are generated in memory */

		/* first: looking for the first empty entry in the gridinfo-array */
		gridindex =1;
		while (gridinfo[gridindex].identification != NULL)
		{ 
			gridindex++;
			if (gridindex == gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no more free grids in memeory, (max. "<<gridmaxanzahl<<")!"<<"[end-entry]";
				throw(12);
			}
		}

		// gridslesen: a variable which is 1 (read all) or 3 (read all and initialize gwtable for each time step from an external grid), 
		// if grids should be read from disk, otherwise it is 0 (from control file) or 2 (read only gwdepth) or 4 (read gwtable for 
		// each time step like 3 but not all the other grids)
		if ((gridslesen==1)||(gridslesen==3))
		{ /* reading stacks from disk */
			/* the theta-stack */
			lv_temp = readstackfromfile(theta_stackname, &theta, _MPIFullGridOff_); /* lv_temp is the number of levels from the stack, read in */
			if (lv_temp!=levels) /* comparing the read number of levels with the maximum number of levels from the soil table */
			{
				(*WSstream) << "FATAL: wrong number of levels in stack "<<gridinfo[gridindex].name
				<<" (there are "<<lv_temp<<" level, it must be "<<levels<<"[end-entry]";
				throw(13);
			}
			theta1 = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(gridindex, theta_stackname, theta1, theta, lv_temp, _Theta_, theta_writestack, mpistartrow, mpiendrow); /* because theta is not a grid but a stack, the adress-parameter for the grid is NULL */
			theta_index = gridindex;
			gridindex++;

			/* hgeo: the stack for topographic altitudes */
			lv_temp = readstackfromfile(hgeo_stackname, &hgeo, _MPIFullGridOff_);
			if (lv_temp!=levels)
			{
				(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<<" (there are "<<lv_temp<<" level, it must be "<<levels<<"[end-entry]";
				throw(13);
			}
			hgeo1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, hgeo_stackname, hgeo1, hgeo, lv_temp, _geodetic_height_, hgeo_writestack, mpistartrow, mpiendrow);
			hgeo_index = gridindex;
// initializing should not be done when the stack is read in anyway. This is important when hgeo has been changed by ice wedge melt			
//			for (int k = mpistartrow; k <= mpiendrow; k++)
//			{
//				for (int j = 1; j <= ncols; j++)
//				{
//					if (zn[k][j] != nodata)
//					{
//						hgeo[k][j][0] = gridinfo[2].pointer[k][j]; /* the zero-level contains then the elevation on the surface from the digital elevation model*/
//						/* all other levels contain the elevations of the lower layer boundary (needed for hydraulic heads calculation) */
//						for (int lv = 1; lv <= levels; lv++)
//						{
//							hgeo[k][j][lv] = hgeo[k][j][lv - 1] - ST[(int)SOIL[k][j]].getDzui()[lv];
//						}
//					}
//					else
//					{
//						for (int lv = 1; lv <= levels; lv++)
//						{
//							hgeo[k][j][lv] = nodata;
//						}
//					}
//				}
//			}
			
			gridindex++;

			/* hh: the stack for hydraulic heads */
			lv_temp = readstackfromfile(hh_stackname, &hh, _MPIFullGridOff_);
			if (lv_temp!=levels)
			{
				(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<<" (there are "<<lv_temp<<" level, it must be "<<levels<<"[end-entry]";
				throw(13);
			}

			hh1 = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(gridindex, hh_stackname, hh1, hh, lv_temp, _hydraulic_head_, hh_writestack, mpistartrow, mpiendrow);
			hh_index = gridindex;
			gridindex++;

			// creating the lower boundary condition for temperature modelling if it was not read in as standard grid
			// this is possible only after reading in hgeo or creating hgeo
			if (do_TMod>0 && tlowbndindex==-99)
			{
				TLOWBND = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, _T_Lower_Boundary_Condition_, TLOWBND, NULL, 0, _T_Lower_Boundary_Condition_, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
//					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							TLOWBND[k][j] = dTModLowBoundBase + dTModLapseRate*hgeo[k][j][0];
						}
						else
						{
							TLOWBND[k][j] = nodata;
						}
					}
				}
				tlowbndindex = gridindex;
				// writing the generated lower boundary condition
				char filename[500];
				sprintf(filename,"%s_%.03d",_T_Lower_Boundary_Condition_,mpi_rank);
				printgridfile(filename, TLOWBND, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				gridindex++;
			}

			if (do_TMod>0)
			{
				/* temper: the stack for soil temperatures */
				lv_temp = readstackfromfile(temper_stackname, &temper,_MPIFullGridOff_);
				if (lv_temp != levels)
				{
					(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
					<<" (there are "<<lv_temp<<" level, it must be "<<levels<<"[end-entry]";
					throw(13);
				}

				temper1 = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				initialize_grid(gridindex, temper_stackname, temper1, temper, lv_temp, _Soil_Temperature_Stack_, temper_writestack, mpistartrow, mpiendrow);
				temperindex = gridindex;
				gridindex++;
			}

			if (do_ArcVeg)
			{
				/* : the stack for hydraulic heads */
				lv_temp = readstackfromfile(ExcessIceFactor_stackname, &ExcessIceFactor, _MPIFullGridOff_);
				if (lv_temp != levels)
				{
					(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
						<< " (there are " << lv_temp << " level, it must be " << levels << "[end-entry]";
					throw(13);
				}
				ExcessIceFactor1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, ExcessIceFactor_stackname, ExcessIceFactor1, ExcessIceFactor, lv_temp, _ExcessIceFactor_, ExcessIceFactor_writestack, mpistartrow, mpiendrow);
				ExcessIceFactorindex = gridindex;
				gridindex++;

				// arcVeg inundation helper grid (for computing 7-day-moving average of inundation between start day (plus waiting period) and end day)
				lv_temp = readstackfromfile(AVINHLP_stackname, &AVINHLP, _MPIFullGridOff_);
				if (lv_temp != 8)
				{
					(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
						<< " (there are " << lv_temp << " levels, but expected are 8 levels" << "[end-entry]";
					throw(13);
				}
				AVINHLP1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVINHLP_stackname, AVINHLP1, AVINHLP, lv_temp, _AVINHLP_, AVINHLP_writestack, mpistartrow, mpiendrow);
				AVINHLPindex = gridindex;
				gridindex++;
//				(*WSstream) << "INFO: AVINHLP read in" << "[end-entry]";

				readgridfromfile(AVGDD_name, &AVGDD, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVGDD_name, AVGDD, NULL, 0, _AVGDD_, AVGDD_writegrid, mpistartrow, mpiendrow);
				AVGDDindex = gridindex;
				gridindex++;

				readgridfromfile(AVINMAX7DAVG_name, &AVINMAX7DAVG, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVINMAX7DAVG_name, AVINMAX7DAVG, NULL, 0, _AVINMAX7DAVG_, AVINMAX7DAVG_writegrid, mpistartrow, mpiendrow);
				AVINMAX7DAVGindex = gridindex;
				gridindex++;

				readgridfromfile(AVTDAVG_name, &AVTDAVG, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVTDAVG_name, AVTDAVG, NULL, 0, _AVTDAVG_, AVTDAVG_writegrid, mpistartrow, mpiendrow);
				AVTDAVGindex = gridindex;
				gridindex++;

				readgridfromfile(AVINDAVG_name, &AVINDAVG, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVINDAVG_name, AVINDAVG, NULL, 0, _AVINDAVG_, AVINDAVG_writegrid, mpistartrow, mpiendrow);
				AVINDAVGindex = gridindex;
				gridindex++;

				readgridfromfile(AVSNOWFREE_name, &AVSNOWFREE, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVSNOWFREE_name, AVSNOWFREE, NULL, 0, _AVSNOWFREE_, AVSNOWFREE_writegrid, mpistartrow, mpiendrow);
				AVSNOWFREEindex = gridindex;
				gridindex++;

				readgridfromfile(AVSON_name, &AVSON, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(gridindex, AVSON_name, AVSON, NULL, 0, _AVSON_, AVSON_writegrid, mpistartrow, mpiendrow);
				AVSONindex = gridindex;
				gridindex++;
//				(*WSstream) << "INFO: AVSON read in" << "[end-entry]";

			}

			/* now the grids */

			/* groundwater table [m a.s.l.) */
			
			readgridfromfile(GWSTANDgridname, &GWSTAND, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
			initialize_grid(gridindex, GWSTANDgridname, GWSTAND, NULL, 0, _Grundwasserstand_, GWSTANDwritegrid, (mpi_rank == 0) ? 1 : sr, (mpi_rank == 0) ? nrows : er);
			gwstandindex = gridindex;
			gridindex++;

			/* level index of the layer, in which the groundwater table is laying */
			readgridfromfile(GWLEVELgridname, &GWLEVEL, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(gridindex, GWLEVELgridname, GWLEVEL, NULL, 0, _Grundwasserlevel_, GWLEVELwritegrid, mpistartrow, mpiendrow);
			gwlevelindex = gridindex;
			gridindex++;

			/* the theta-value of the layer, in which the groundwater table is
	(needed for estimation of groundwater depth even between layer boundaries) */
			readgridfromfile(GWthetagridname, &GWtheta, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(gridindex, GWthetagridname, GWtheta, NULL, 0, _Theta_in_GWLEVEL_, GWthetawritegrid, mpistartrow, mpiendrow);
			gwthetaindex = gridindex;
			gridindex++;

			/* baseflow: it's the balance of lateral inflows and outflows of the groundwater model when using gwflow. 
	otherwise its the empiric baseflow from unsatzonmodel*/
			readgridfromfile(DELTAGWname, &DELTAGW, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOff_); // we need the full grid in instance 0, so set _MPIFullGridOn_
			initialize_grid(gridindex, DELTAGWname, DELTAGW, NULL, 0, _Basisabfluss_, DELTAGWwritegrid, (mpi_rank==0)?1:mpistartrow, (mpi_rank==0)?nrows:mpiendrow);
			deltagwindex = gridindex;
			// do not change the value that was stored in gw model 2014-07-18
//			for (int k=1;k<=nrows;k++) 
//			{
//				for (int j=1;j<=ncols;j++)
//				{
//					if (zn[k][j]!=nodata)
//					{
//						DELTAGW[k][j]/=timestep;
//					}
//				}
//			}
			gridindex++;

			if (MAXPOND) /* is a grid defined, containing maximum pond storage capacities? */
			{
				/* actual content of pond storage */
				readgridfromfile(PONDgridname, &POND, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
				initialize_grid(gridindex, PONDgridname, POND, NULL, 0, _ponding_storage_, PONDwritegrid, (mpi_rank==0)?1:sr, (mpi_rank==0)?nrows:er);
				pondindex = gridindex;
				gridindex++;
			}
		}
		else /* stacks and grids are not read from disk but generated internally */
		{
			if (gridslesen==2) /* only gwstand (gw-table) is read, from this data all the others will be derived if possible */
			{
				/* groundwater table [m a.s.l.) */
				readgridfromfile(GWSTANDgridname, &GWSTAND, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
				initialize_grid(gridindex, GWSTANDgridname, GWSTAND, NULL, 0, _Grundwasserstand_, GWSTANDwritegrid, (mpi_rank == 0) ? 1 : sr, (mpi_rank == 0) ? nrows : er);
				gwstandindex = gridindex;
				gridindex++;
				if (MAXPOND) /* is a grid defined, containing maximum pond storage capacities? */
				{
					/* actual content of pond storage */
					readgridfromfile(PONDgridname, &POND, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOn_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
					initialize_grid(gridindex, PONDgridname, POND, NULL, 0, _ponding_storage_, PONDwritegrid, (mpi_rank==0)?1:sr, (mpi_rank==0)?nrows:er);
					pondindex = gridindex;
					gridindex++;
				}
			}
			else
			{
				/* groundwater level [m a.s.l.] */
				GWSTAND = matrix((mpi_rank == 0) ? 1 : sr, (mpi_rank == 0) ? nrows : er, 1, ncols, nrows);
				initialize_grid(gridindex, GWSTANDgridname, GWSTAND, NULL, 0, _Grundwasserstand_, GWSTANDwritegrid, (mpi_rank == 0) ? 1 : sr, (mpi_rank == 0) ? nrows : er);
				gwstandindex = gridindex;
				gridindex++;
				
				/* pond storage actual value */
				if (MAXPOND)
				{
					// in channel routing we need the entire grid when using the lake model, so this grid should be read in completely in rank 0
					int loc_sr = (mpi_rank == 0) ? 1 : sr;
					int loc_er = (mpi_rank == 0) ? nrows : er;

					POND = matrix(loc_sr,loc_er, 1, ncols,nrows);
					initialize_grid(gridindex, PONDgridname, POND, NULL, 0, _ponding_storage_, PONDwritegrid, loc_sr,loc_er);
					for (int k = loc_sr; k <= loc_er; k++)
					{
						for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
						{
							if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
							{
								POND[k][j] = 0;
							}
							else
							{
								POND[k][j] = nodata;
							}
						}
					}
					pondindex = gridindex;
					gridindex++;
				}
			}

			/* 1. stacks */
			/* geographical altitude [m a.s.l.] */
			hgeo = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels,nrows);
			hgeo1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, hgeo_stackname, hgeo1, hgeo, levels, _geodetic_height_, hgeo_writestack, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						hgeo[k][j][0] = gridinfo[2].pointer[k][j]; /* the zero-level contains then the elevation on the surface from the digital elevation model*/
						/* all other levels contain the elevations of the lower layer boundary (needed for hydraulic heads calculation) */
						for (int lv = 1; lv <= levels; lv++)
						{
							hgeo[k][j][lv] = hgeo[k][j][lv - 1] - ST[(int)SOIL[k][j]].getDzui()[lv];
						}
					}
					else
					{
						for (int lv=1;lv<=levels;lv++) 
						{
							hgeo[k][j][lv]=nodata;
						}
					}
				}
			}
			hgeo_index = gridindex;
			gridindex++;

			/* a stack for water content for all levels (theta-values) */
			theta = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels, nrows);
			theta1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, theta_stackname, theta1, theta, levels, _Theta_, theta_writestack, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						if (gridslesen == 2)
						{
							gw_depth_ini=0;
							while ((hgeo[k][j][0]-hgeo[k][j][gw_depth_ini+1] < -GWSTAND[k][j]) && (gw_depth_ini<levels))
							{
								gw_depth_ini++;
							}

					
							//              if (-GWSTAND[k][j]<ST[(int)SOIL[k][j]].getDzui()[1])
							//              {
							//                if (-GWSTAND[k][j]<ST[(int)SOIL[k][j]].getDzui()[1]/2.0) gw_depth_ini=0;
							//  //                if (-GWSTAND[k][j]<_Dicke_Schicht_1_/2.0) gw_depth_ini=0;
							//                else
							//                  if (-GWSTAND[k][j]<((ST[(int)SOIL[k][j]].getDzui()[0]-_Dicke_Schicht_1_)/2.0+_Dicke_Schicht_1_))
							//                    gw_depth_ini=1;
							//                  else
							//                    gw_depth_ini=2;
							//              }                    
							//              else
							//                gw_depth_ini= +1-floor(GWSTAND[k][j]/ST[(int)SOIL[k][j]].getDzui()[0]+0.4999); /* groundwater level is rounded to the nearest layer boundary */
						}
						else
						//            gw_depth_ini= levels/4; /* groundwater level is assumed to stay in 1/5 of the depth of the soil column */
						gw_depth_ini= int(1.4999+floor(double(levels-1)/double(ini_gw_level))); /* initialization groundwater level: 1 .. levels (most 4) */
				
						/* th_start_sat: saturation water content */
						th_start_sat = ST[(int)SOIL[k][j]].getThetas(1);
						/* SOIL[k][j] is the code for the soiltable at the location k,j;
				the recored theta_h contains the theta-values for the theta-h-relation */
						theta[k][j][0]=th_start_sat; /* initializing the zero-level of theta-stack with saturation water content */
						for (int lv=1;lv<=gw_depth_ini;lv++)
						{
							//          theta[k][j][lv]=get_theta(ST[(int)SOIL[k][j]].anz_h,
							//                                    0.3333,   /* all levels have a water content equivalent to a suction of 0.3333 m */
							//                                    ST[(int)SOIL[k][j]].h_tab,
							//                                    ST[(int)SOIL[k][j]].theta_h);

							/* all levels are initialized with the same potential -> groundwater level potential: there are no fluxes at all at the beginning */
							theta[k][j][lv]=ST[(int)SOIL[k][j]].get_theta(lv,
							(hgeo[k][j][lv]+ST[(int)SOIL[k][j]].getDzui()[lv]*0.5)-
							(hgeo[k][j][gw_depth_ini+1]+ST[(int)SOIL[k][j]].getDzui()[gw_depth_ini+1]*0.5)
							);

						}
						for (int lv=gw_depth_ini+1;lv<=levels;lv++) 
						{
							//				theta[k][j][lv]=th_start_sat;
							theta[k][j][lv]=ST[(int)SOIL[k][j]].getThetas(lv);
						}
					}
					else
					{
						for (int lv=1;lv<=levels;lv++) 
						{
							theta[k][j][lv]=nodata;
						}
					}
				}
			}
			theta_index = gridindex;
			gridindex++;

			/* hydraulic head */
			hh = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels,nrows);
			hh1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, hh_stackname, hh1, hh, levels, _hydraulic_head_, hh_writestack, mpistartrow, mpiendrow);
			double h_temp;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if (zn[k][j] != nodata)
					{
						hh[k][j][0] = hgeo[k][j][0];
						for (int lv = 1; lv <= levels; lv++)
						{
							h_temp = ST[(int)SOIL[k][j]].h_calc(lv, /* calculating suction for the given theta-value */
								theta[k][j][lv]);
							hh[k][j][lv] = (hgeo[k][j][lv - 1] + hgeo[k][j][lv])*0.5 - h_temp;
							/* hydraulic head is calculated from mean geograohical elevation minus suction */
						}
					}
					else
					{
						for (int lv = 1; lv <= levels; lv++)
						{
							hh[k][j][lv] = nodata;
						}
					}
				}
			}
			hh_index = gridindex;
			gridindex++;

			// creating the lower boundary condition for temperature modelling if it was not read in as standard grid
			// this is possible only after reading in hgeo or creating hgeo
			if (do_TMod>0 && tlowbndindex==-99)
			{
				TLOWBND = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				initialize_grid(gridindex, _T_Lower_Boundary_Condition_, TLOWBND, NULL, 0, _T_Lower_Boundary_Condition_, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
//					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							TLOWBND[k][j] = dTModLowBoundBase + dTModLapseRate*hgeo[k][j][0];
						}
						else
						{
							TLOWBND[k][j] = nodata ;
						}
					}
				}
				tlowbndindex = gridindex;
				// writing the generated lower boundary condition
                // from netcdf code >>>
                // printgridfile(_T_Lower_Boundary_Condition_,TLOWBND,WriteASCIIGrids,WriteDoubleGrids,WriteNetcdfGrids);
                // <<< from netcdf code
				char filename[500];
				sprintf(filename,"%s_%.03d",_T_Lower_Boundary_Condition_,mpi_rank);
				printgridfile(filename, TLOWBND, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				gridindex++;
			}

			if (do_TMod>0)
			{
				/* soil temperatures */
				temper = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels,nrows);
				temper1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, temper_stackname, temper1, temper, levels, _Soil_Temperature_Stack_, temper_writestack, mpistartrow, mpiendrow);
				double temper_temp,actdiff,maxdiff;
				do_initializeTemper=true;
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							temper1[k][j] = 0;
						}
						else
						{
							temper1[k][j]=nodata;
						}
					}
				}
/* will be done in first time step, since temperature is known only then
				for (k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							temper[k][j][0]=TLOWBND[k][j];
							maxdiff = hgeo[k][j][0]-hgeo[k][j][levels];
							for (lv=1;lv<=levels;lv++)
							{
								// todo: linear interpolieren zwischen TLOWBND und aktueller Obergrenze (T Luft oder alles konstant)
								actdiff = hgeo[k][j][0]-(hgeo[k][j][lv-1]+hgeo[k][j][lv])/2; // elevation is taken from center of layer
								temper[k][j][lv]=TLOWBND[k][j]*(actdiff/maxdiff)+T[k][j]*(1.0-actdiff/maxdiff);
							}
						}
						else
						{
							for (lv=1;lv<=levels;lv++) 
							{
								temper[k][j][lv]=nodata;
							}
						}
					}
				}
*/
				temperindex = gridindex;
				gridindex++;
			}

			if (do_ArcVeg)
			{
				/* hydraulic head */
				ExcessIceFactor = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels, nrows);
				ExcessIceFactor1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, ExcessIceFactor_stackname, ExcessIceFactor1, ExcessIceFactor, levels, _ExcessIceFactor_, ExcessIceFactor_writestack, mpistartrow, mpiendrow);

				// reading helper grids with info in upper and lower ice wedge depth for creating excess ice stack
				ICEWEDGEUPPERDEPTH = gridinfo[getindex(_IceWedgeUpperDepth_)].pointer;
				ICEWEDGELOWERDEPTH = gridinfo[getindex(_IceWedgeLowerDepth_)].pointer;
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						for (int lv = 1; lv <= levels; lv++)
						{
							// initialilzing everything in the excessIceStack with nodata first
							ExcessIceFactor[k][j][lv] = nodata;
						}
					}
				}
				ExcessIceFactorindex = gridindex;
				gridindex++;

				// arcVeg inundation helper grid (for computing 7-day-moving average of inundation between start day (plus waiting period) and end day)
				AVINHLP = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, 8, nrows);
				AVINHLP1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVINHLP_stackname, AVINHLP1, AVINHLP, 8, _AVINHLP_, AVINHLP_writestack, mpistartrow, mpiendrow);
				AVINHLPindex = gridindex;
				gridindex++;

				AVGDD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVGDD_name, AVGDD, NULL, 0, _AVGDD_, AVGDD_writegrid, mpistartrow, mpiendrow);
				AVGDDindex = gridindex;
				gridindex++;

				AVINMAX7DAVG = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVINMAX7DAVG_name, AVINMAX7DAVG, NULL, 0, _AVINMAX7DAVG_, AVINMAX7DAVG_writegrid, mpistartrow, mpiendrow);
				AVINMAX7DAVGindex = gridindex;
				gridindex++;

				AVTDAVG = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVTDAVG_name, AVTDAVG, NULL, 0, _AVTDAVG_, AVTDAVG_writegrid, mpistartrow, mpiendrow);
				AVTDAVGindex = gridindex;
				gridindex++;

				AVINDAVG = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVINDAVG_name, AVINDAVG, NULL, 0, _AVINDAVG_, AVINDAVG_writegrid, mpistartrow, mpiendrow);
				AVINDAVGindex = gridindex;
				gridindex++;

				AVSNOWFREE = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVSNOWFREE_name, AVSNOWFREE, NULL, 0, _AVSNOWFREE_, AVSNOWFREE_writegrid, mpistartrow, mpiendrow);
				AVSNOWFREEindex = gridindex;
				gridindex++;

				AVSON = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, AVSON_name, AVSON, NULL, 0, _AVSON_, AVSON_writegrid, mpistartrow, mpiendrow);
				AVSONindex = gridindex;
				gridindex++;

				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					// set nodata to all cells outside the actual first-last range
					for (int j = 1; j <= ncols; j++)
					{
						if ((j< FirstLastCol[k].First)||(j> FirstLastCol[k].Last))
						{
							AVINHLP1[k][j] = nodata;
							AVGDD[k][j] = nodata;
							AVINMAX7DAVG[k][j] = nodata;
							AVTDAVG[k][j] = nodata;
							AVINDAVG[k][j] = nodata;
							AVSNOWFREE[k][j] = nodata;
							AVSON[k][j] = nodata;
							for (int lv = 0; lv <= 8; lv++)
							{
								// initialilzing everything in the excessIceStack with nodata first
								AVINHLP[k][j][lv] = nodata;
							}
						}
					}
					// initialize all cells from first to last in this row with 0. If a cell has no zone value or another mpi rank, set it to nodata
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							AVINHLP1[k][j] = 0;
							AVGDD[k][j] = 0;
							AVINMAX7DAVG[k][j] = 0;
							AVTDAVG[k][j] = 0;
							AVINDAVG[k][j] = 0;
							AVSNOWFREE[k][j] = 0;
							AVSON[k][j] = 0;
							for (int lv = 0; lv <= 8; lv++)
							{
								// initialilzing everything in the excessIceStack with nodata first
								AVINHLP[k][j][lv] = 0;
							}
						}
						else
						{
							AVINHLP1[k][j] = nodata;
							AVGDD[k][j] = nodata;
							AVINMAX7DAVG[k][j] = nodata;
							AVTDAVG[k][j] = nodata;
							AVINDAVG[k][j] = nodata;
							AVSNOWFREE[k][j] = nodata;
							AVSON[k][j] = nodata;
							for (int lv = 0; lv <= 8; lv++)
							{
								// initialilzing everything in the excessIceStack with nodata first
								AVINHLP[k][j][lv] = nodata;
							}
						}
					}
				}
			}

			/* now generating grids */

			/* baseflow: it's the balance of lateral inflows and outflows of the groundwater model when using gwflow. 
	otherwise its the empiric baseflow from unsatzonmodel*/

			// when running lake model, DELTAGW is required as full grid in rank 0
			int srdgw0=(mpi_rank==0)?1:mpistartrow;
			int erdgw0=(mpi_rank==0)?nrows:mpiendrow;

			DELTAGW = matrix(srdgw0, erdgw0, 1,ncols,nrows);
			initialize_grid(gridindex, DELTAGWname, DELTAGW, NULL, 0, _Basisabfluss_, DELTAGWwritegrid, srdgw0, erdgw0);
			for (int k = srdgw0; k <= erdgw0; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						DELTAGW[k][j] = 0;
					}
					else
					{
						DELTAGW[k][j]=nodata;
					}
				}
			}
			deltagwindex = gridindex;
			gridindex++;

			/* index of soil layer in which the groundwater table is */
			GWLEVEL = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(gridindex, GWLEVELgridname, GWLEVEL, NULL, 0, _Grundwasserlevel_, GWLEVELwritegrid, mpistartrow, mpiendrow);
			gwlevelindex = gridindex;
			gridindex++;

			/* theta limit of the layer containig groundwater level */
			GWtheta = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(gridindex, GWthetagridname, GWtheta, NULL, 0, _Theta_in_GWLEVEL_, GWthetawritegrid, mpistartrow, mpiendrow);
			gwthetaindex = gridindex;
			gridindex++;

			double temp;
			/* now initializing groundwater table (GWSTAND),  groundwater level index (GWLEVEL) and groundwater level theta limit (GWtheta) */
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						temp = 0;
						//            th_start_sat = ST[(int)SOIL[k][j]].getThetas(1);
						GWtheta[k][j] = theta[k][j][1];
						if (gridslesen!=2) GWSTAND[k][j] = hgeo[k][j][levels]-hgeo[k][j][0]; /* start condition: groundwater level is on the lower end of the soil column */
						for (int lv=levels-1;lv>0;lv--) /* from there the search is started: when is occuring the first level containing less water than saturation water content */
						{
							if (theta[k][j][lv]> (ST[(int)SOIL[k][j]].getThetas(lv)*0.9995) && temp!=2)
							{
								if (gridslesen!=2) GWSTAND[k][j] = hgeo[k][j][lv-1]-hgeo[k][j][0];
								GWLEVEL[k][j] = lv;  /* if g.w. is within the first level: GWLEVEL=1; within the 2nd level: GWLEVEL=2; within the 3rd: 3 a.s.o. */
								GWtheta[k][j] = theta[k][j][lv-1];
								temp=1;
							}
							else if (temp==1 || theta[k][j][levels]!=ST[(int)SOIL[k][j]].getThetas(lv))
							temp=2;
						}
					}
					else
					{
						GWSTAND[k][j] = GWLEVEL[k][j] = GWtheta[k][j] = nodata;
					}
				}
			}
		}

		// Todo 20241122
		// find grids ARCVEGBIOMASS analog zu interzeption module
			if (do_ArcVeg > 0)
			{
				arcVegNumLayers = MultiLayerLUTable[1].getNumFractions();
				arcVegNumBiomass = MultiLayerLUTable[1].getNumBiomass();
				arcVegNumRetrans = MultiLayerLUTable[1].getNumRetrans();

				if ((arcVegNumLayers != iNumVegetationLayers)||(arcVegNumLayers != arcVegNumBiomass)||(arcVegNumLayers != arcVegNumRetrans))
				{
					(*WSstream) << "FATAL: ERROR in initializing atcVeg: [ArcVegIntegration] is active but there are no or not all fraction and/or biomass and/or retrans grids found as standard grids (ArcVeg_Fraction_1..n and ArcVeg_Biomass_1...n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}

				bool bAllArcVegFractionsFound = true;
				bool bAllArcVegBiomassFound = true;
				bool bAllArcVegRetransFound = true;

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
					(*WSstream) << "FATAL: ERROR in initializing interception: [ArcVegIntegration] is active but there are no or not all fraction grids found as standard grids (ArcVeg_Fraction_1..n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}

				if (arcVegNumLayers > 0)
				{
					char* layername;
					layername = new char[200];
					arcVegNumLayers = iNumVegetationLayers;
					ARCVEGBIOMASS = new double** [arcVegNumLayers + 1];

					for (int layer = 1; layer <= arcVegNumLayers; layer++)
					{
						sprintf(layername, "%s%d", "ArcVeg_Biomass_", layer);
						if (ist_vorhanden(layername))
						{
							ARCVEGBIOMASS[layer] = gridinfo[getindex(layername)].pointer;
						}
						else
						{
							bAllArcVegBiomassFound = false;
						}
					}

					delete[] layername;
				}
				else
				{
					bAllArcVegBiomassFound = false;
				}
				if (!bAllArcVegBiomassFound)
				{
					(*WSstream) << "FATAL: ERROR in initializing ArcVeg: [ArcVegIntegration] is active but there are no or not all biomass grids found as standard grids (ArcVeg_Biomass_1..n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}

				if (arcVegNumLayers > 0)
				{
					char* layername;
					layername = new char[200];
					arcVegNumLayers = iNumVegetationLayers;
					ARCVEGRETRANS = new double** [arcVegNumLayers + 1];

					for (int layer = 1; layer <= arcVegNumLayers; layer++)
					{
						sprintf(layername, "%s%d", "ArcVeg_Retrans_", layer);
						if (ist_vorhanden(layername))
						{
							ARCVEGRETRANS[layer] = gridinfo[getindex(layername)].pointer;
						}
						else
						{
							bAllArcVegRetransFound = false;
						}
					}

					delete[] layername;
				}
				else
				{
					bAllArcVegRetransFound = false;
				}
				if (!bAllArcVegRetransFound)
				{
					(*WSstream) << "FATAL: ERROR in initializing ArcVeg: [ArcVegIntegration] is active but there are no or not all Retrans grids found as standard grids (ArcVeg_Retrans_1..n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}

			}

		



			/* now  some grids and a stack are newly generated even if grids and stacks are read in
(these are no storage grids, so they dont have to be stored!) */

		/* outflux from the layers (Q-bottom) */
		qu = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels,nrows);
		qu1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, qu_stackname, qu1, qu, levels, _soil_bottom_fluxes_, qu_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv = 1; lv <= levels; lv++)
					{
						qu[k][j][lv] = 0;
					}
					qu1[k][j] = 0;
				}
				else 
				{
					for (int lv = 1; lv <= levels; lv++)
					{
						qu[k][j][lv] = nodata;
					}
					qu1[k][j] = nodata;
				}
			}
		}
		qu_index = gridindex;
		gridindex++;

		// SE relative saturation
		if (do_TMod > 0)
		{
			/* soil temperatures */
			SEstack = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels, nrows);
			SE1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, se_stackname, SE1, SEstack, levels, _Soil_SE_Stack_, se_writestack, mpistartrow, mpiendrow);
			double temper_temp, actdiff, maxdiff;

			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						SE1[k][j] = 0;
					}
					else
					{
						SE1[k][j] = nodata;
					}
				}
			}
			seindex = gridindex;
			gridindex++;
		}

		// reserve a full grid for qu1_ext in rank 0
		sprintf(name,"%s%s",_soil_bottom_fluxes_,"ext");
		int srext = (mpi_rank==0)?1:mpistartrow;
		int erext = (mpi_rank==0)?nrows:mpiendrow;
		qu1_ext = matrix(srext,erext,1, ncols, nrows);
		initialize_grid(gridindex, name, qu1_ext , NULL, levels, name, 0, srext,erext);
		for (int k = srext; k <= erext; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					qu1_ext[k][j] = 0;
				}
				else 
				{
					qu1_ext[k][j] = nodata;
				}
			}
		}
		qu1extindex = gridindex;
		gridindex++;

		/* soil moisture in root zone is used also in irrigation initializing, thats why it's initialized before the other grids */
		SB_1_ = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(gridindex, SB_1_gridname, SB_1_, NULL, 0, _SB_1_grid_, SB_1_writegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					SB_1_[k][j] = 0;
				}
				else
				{
					SB_1_[k][j] = nodata;
				}
			}
		}
		sb_1_index = gridindex;
		gridindex++;

		// mean suction within root zone, calculated from mean soil moisture (SB_1_) for use in LUtable (dynamic phenology 3)
		THROOT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, _THROOT_grid_, THROOT, NULL, 0, _THROOT_grid_, THROOTwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					THROOT[k][j] = 0;
				}
				else
				{
					THROOT[k][j] = nodata;
				}
			}
		}
		throotindex = gridindex;
		gridindex++;


		/* grids only necessary for irrigation..
if irrigationtable is NULL, the "do_it" - parameter in [irrigation] was set to 0, 
this was checked within defaults::defaults (reading standard grids etc.). 
if irrigation should be considered (do=it in [irrigation]==1), a grid with 
identification string "irrigation_codes" must be available -> has to be read in as a standard grid
*/
		if (irrigationtable) 
		{
			IRRIGCODE = gridinfo[getindex(_irrigation_codes_)].pointer;

			IRRIG = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, IRRIGgridname, IRRIG, NULL, 0, _irrigation_amounts_, IRRIGwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						IRRIG[k][j] = 0;
					}
					else
					{
						IRRIG[k][j] = nodata;
					}
				}
			}
			irrigindex = gridindex;
			gridindex++;

			MAXPUMP = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, _Max_Pumping_, MAXPUMP, NULL, 0, _Max_Pumping_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						MAXPUMP[k][j] = 0;
					}
					else
					{
						MAXPUMP[k][j] = nodata;
					}
				}
			}
			maxpumpindex = gridindex;
			gridindex++;

			
			SBLOW = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, _SBLOW_grid_, SBLOW, NULL, 0, _SBLOW_grid_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						s = (int)SOIL[k][j];
						SBLOW[k][j] = IRRIGCODE[k][j] == nodata ? nodata :
							ST[s].get_theta(2, irrigationtable[(int)IRRIGCODE[k][j]].psi_dry);
					}
					else
					{
						SBLOW[k][j] = nodata;
					}
				}
			}
			sblowindex = gridindex;
			gridindex++;

			SBHIGH = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, _SBHIGH_grid_, SBHIGH, NULL, 0, _SBHIGH_grid_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						s = (int)SOIL[k][j];
						if (IRRIGCODE[k][j] != nodata)
						{
							if (irrigationtable[(int)IRRIGCODE[k][j]].methode == 4)
							{
								SBHIGH[k][j] = ST[s].get_theta(2, LUtable[(int)LUSE[1][k][j]].getReduParDry());
							}
							else
							{
								SBHIGH[k][j] = ST[s].get_theta(2, irrigationtable[(int)IRRIGCODE[k][j]].psi_wet);
							}
						}
						else
						{
							SBHIGH[k][j] = nodata;
						}
					}
					else
					{
						SBHIGH[k][j] = nodata;
					}
				}
			}
			sbhighindex = gridindex;
			gridindex++;

			// check if grid ETR2ETP exists (if read in as variable grid, otherwise create it here (but without output options!)
			if (ist_vorhanden(_ETR2ETP_))
			{
				etr2etpindex=getindex(_ETR2ETP_);
				ETR2ETP=gridinfo[etr2etpindex].pointer;
				etr2etpindex=-etr2etpindex; // when deleting the grid, a negative index indicates that the grid must be deleted within this module
			}
			else
			{
				ETR2ETP = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				initialize_grid(gridindex, "", ETR2ETP, NULL, 0, _ETR2ETP_, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							ETR2ETP[k][j] = 1.0;
						} else	
						{
							ETR2ETP[k][j] = nodata;
						}
					}
				}
				etr2etpindex = gridindex;
				gridindex++;
			}

			
			/* initializing first mean soil moisture for irrigation control */
			double  *theta_sl;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						s = (int)SOIL[k][j];
						dicken=ST[s].getDzui();
						theta_sl = &(theta[k][j][0]);
						if (iNumVegetationLayers==0)
						{
							if (WURZ[0][k][j] > 0)
							{
								SB_1_[k][j] = GetMeanTheta(s, dicken, theta_sl, 0.0, WURZ[0][k][j], 1) / WURZ[0][k][j]; //GetMeanTheta(0.0,WURZ[0][k][j]);
							}
							else
							{
								SB_1_[k][j] = theta_sl[1];
							}
							THROOT[k][j]=ST[(int)SOIL[k][j]].h_calc(1,SB_1_[k][j]);
						}
						else
						{
							if (WURZ[1][k][j] > 0)
							{
								SB_1_[k][j] = GetMeanTheta(s, dicken, theta_sl, 0.0, WURZ[1][k][j], 1) / WURZ[1][k][j]; //GetMeanTheta(0.0,WURZ[1][k][j]);
							}
							else
							{
								SB_1_[k][j] = theta_sl[1];
							}
							THROOT[k][j]=ST[(int)SOIL[k][j]].h_calc(1,SB_1_[k][j]);
						}
					}
					else
					{
						SB_1_[k][j] = nodata;
					}
				}
			}

		}
		else
		{
			// search for standard grid ETR2ETP anyway, since it should be possible to write out this value even without the irrigation module beeing active
			if (ist_vorhanden(_ETR2ETP_))
			{
				etr2etpindex = getindex(_ETR2ETP_);
				ETR2ETP = gridinfo[etr2etpindex].pointer;
				etr2etpindex = -etr2etpindex; // when deleting the grid, a negative index indicates that the grid must be deleted within this module
			}
			else
			{ // only if irrigation is switched off AND ETR2ETP is not listed in standard- or variable grids section, ETR2ETP should be ignored
				ETR2ETP = NULL;
			}

			IRRIG=NULL;
			MAXPUMP=NULL;
			//      SB_1_=NULL; is created earlier now
			SBLOW=NULL;
			SBHIGH=NULL;
		}

		/* baseflow: it's the balance of lateral inflows and outflows of the groundwater model when using gwflow. 
otherwise its the empiric baseflow from unsatzonmodel*/
		//    DELTAGW = matrix(1,nrows,1,ncols);
		//    initialize_grid(gridindex,DELTAGWname,DELTAGW,NULL,0,_Basisabfluss_,DELTAGWwritegrid,1,nrows);
		//    for (int k=1;k<=nrows;k++) for (int j=1;j<=ncols;j++)  DELTAGW[k][j] = zn[k][j] == nodata ? nodata : 0;
		//    deltagwindex = gridindex;
		//    gridindex++;
		

		GWTableExtern = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, _GWTableExtern_, GWTableExtern, NULL, 0, _GWTableExtern_, 0, mpistartrow, mpiendrow);
		gwtableexternindex = gridindex;
		gridindex++;

		EXTBALANCE = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, _BALANCE_, EXTBALANCE, NULL, 0, _BALANCE_, 0, mpistartrow, mpiendrow);
		balanceindex = gridindex;
		gridindex++;

		/* konverted slope angle from degree in tangent */
		TANBETA = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, _Tangens_Beta_, TANBETA, NULL, 0, _Tangens_Beta_, 0, mpistartrow, mpiendrow);
		pi_180=pi/180.0;
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					TANBETA[k][j] = SLP[k][j] * pi_180;
				}
				else
				{
					TANBETA[k][j] = nodata;
				}
			}
		}
		tanbetaindex = gridindex;
		gridindex++;

		/* soil moisture in entire unsaturated zone */
		SB_2_ = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, SB_2_gridname, SB_2_, NULL, 0, _SB_2_grid_, SB_2_writegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					SB_2_[k][j] = 0;
				}
				else
				{
					SB_2_[k][j] = nodata;
				}
			}
		}
		sb_2_index = gridindex;
		gridindex++;

		/* soil moisture in entire unsaturated zone */
		MACRO = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, MACROgridname, MACRO, NULL, 0, _macropore_infiltration_, MACROwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					MACRO[k][j] = 0;
				}
				else
				{
					MACRO[k][j] = nodata;
				}
			}
		}
		macroindex = gridindex;
		gridindex++;

		/* number of sub time steps */
		if (methode>2)
		{
			SUBSTEPS = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, SUBSTEPSgridname, SUBSTEPS, NULL, 0, _number_of_substeps_, SUBSTEPSwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
							SUBSTEPS[k][j] = 0;
					}
					else
					{
						SUBSTEPS[k][j] = nodata;
					}
				}
			}
			substepsindex = gridindex;
			gridindex++;
		}

		/* root depth */
		ROOT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, ROOTgridname, ROOT, NULL, 0, _ROOTgrid_, ROOTwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					ROOT[k][j] = 0;
				}
				else
				{
					ROOT[k][j] = nodata;
				}
			}
		}
		rootindex = gridindex;
		gridindex++;

		/* groundwater recharge */
		GWN = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, GWNgridname, GWN, NULL, 0, _Grundwasserneubildung_, GWNwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					GWN[k][j] = 0;
				}
				else
				{
					GWN[k][j] = nodata;
				}
			}
		}
		gwnindex = gridindex;
		gridindex++;

		/* groundwater recharge uprise fraction */
		UPRISE = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, UPRISEgridname, UPRISE, NULL, 0, _Capillar_Uprise_, UPRISEwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					UPRISE[k][j] = 0;
				}
				else
				{
					UPRISE[k][j] = nodata;
				}
			}
		}
		upriseindex = gridindex;
		gridindex++;

		/* groundwater recharge percolation fraction */
		PERCOL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, PERCOLgridname, PERCOL, NULL, 0, _Percolation_, PERCOLwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					PERCOL[k][j] = 0;
				}
				else
				{
					PERCOL[k][j] = nodata;
				}
			}
		}
		percolindex = gridindex;
		gridindex++;

		/* surface runoff */
		// the first component from silting up model will be available if the siltingup-model is used, only
		// to be able to mix the direct runoff components from unsatzonmodel and from silting up model, 
		// we get the address of the runoff grid here
		
		if (ist_vorhanden(_DIRECTRUNOFFSILTINGUP_))
		{
			QDSU=gridinfo[getindex(_DIRECTRUNOFFSILTINGUP_)].pointer;
			m_bAddQdsuToP=false;
		}
		else
		{
			// if silting up model is switched off, then the input into unsatzon model should be taken from 
			// surface routing, if that is switched on. Otherwise, set QDSU to NULL to indicate that no additional
			// input from silting up nor surface routing has to be regarded in unsatzon model
			if (ist_vorhanden(_SURFSTORSILTINGUP_))
			{
				QDSU = gridinfo[getindex(_SURFSTORSILTINGUP_)].pointer;
				m_bAddQdsuToP=true;
			}
			else
			{
				QDSU=NULL;
				m_bAddQdsuToP=false;
			}
		}

		QDgrid = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, QDgridname, QDgrid, NULL, 0, _Direktabfluss_, QDwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					QDgrid[k][j] = 0;
				}
				else
				{
					QDgrid[k][j] = nodata;
				}
			}
		}
		qdindex = gridindex;
		gridindex++;

		/* interflow */
		QIgrid = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, QIgridname, QIgrid, NULL, 0, _Interflow_, QIwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					QIgrid[k][j] = 0;
				}
				else
				{
					QIgrid[k][j] = nodata;
				}
			}
		}
		qiindex = gridindex;
		gridindex++;

//		oldDELTAGW = matrix(mpistartrow, mpiendrow,1,ncols,nrows);
//		initialize_grid(gridindex,"_oldDELTAGW_",oldDELTAGW,NULL,0,"_oldDELTAGW_",0,mpistartrow, mpiendrow);
//		for (int k=1;k<=nrows;k++) for (int j=1;j<=ncols;j++)  oldDELTAGW[k][j] = zn[k][j] == nodata ? nodata : 0;
//		olddeltagwindex = gridindex;
//		gridindex++;

		/* infitration excess */
		INFEX = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, INFEXgridname, INFEX, NULL, 0, _Infiltrationsexzess_, INFEXwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					INFEX[k][j] = 0.0;
				}
				else
				{
					INFEX[k][j] = nodata;
				}
			}
		}
		infexindex = gridindex;
		gridindex++;

		/* infitration excess */
		INFILT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, _Total_Infiltration_, INFILT, NULL, 0, _Total_Infiltration_, 0, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					INFILT[k][j] = 0.0;
				}
				else
				{
					INFILT[k][j] = nodata ;
				}
			}
		}
		infiltindex = gridindex;
		gridindex++;

		/* a grid used as a flag in the infitration routine, if at the interval end saturation occures */
		SATT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, SATTgridname, SATT, NULL, 0, _Infiltrationssaettigung_, SATTwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					SATT[k][j] = 0.0;
				}
				else
				{
					SATT[k][j] = nodata ;
				}
			}
		}
		sattindex = gridindex;
		gridindex++;

		/* a grid containing the minimum water content without reduction of etp */
		//    THREDU = matrix(1,nrows,1,ncols);
		//    initialize_grid(gridindex,_THETA_LIMIT_FOR_ETP_REDUCTION_,THREDU,NULL,0,_THETA_LIMIT_FOR_ETP_REDUCTION_,0);
		//    (*WSstream) << "calculating theta for evaporation reduction"<<"[end-entry]";
		//    for (int k=1;k<=nrows;k++) 
		//    {
		//      for (int j=1;j<=ncols;j++) 
		//      {
		//        if (zn[k][j]>nodata)
		//        {
		//          s=(int)SOIL[k][j];
		//          THREDU[k][j] = zn[k][j]==nodata ? nodata : 
		//            ST[s].get_theta(ST[s].getAnz_h(),LUtable[(int)LUSE[k][j]].getReduParDry(),ST[s].getH_tab(),ST[s].getTheta_h());
		//          if (THREDU[k][j]<=ST[s].getThetar())
		//          {
		//            THREDU[k][j]=(2*ST[s].getThetar()+ST[s].getThetas())/3.0; /* standard-value: weighted average to residual water content */
		//            (*WSstream) << "*";
		//          }
		//        }
		//        else THREDU[k][j]=nodata;
		//      }
		//      (*WSstream) << ".";
		//    }
		//    (*WSstream) << "[end-entry]";
		//    threduindex = gridindex;
		//    gridindex++;

		if (this->do_TMod>0)
		{
			int sr = (mpistartrow>1)?mpistartrow-1:mpistartrow;
			int er = (mpiendrow<nrows)?mpiendrow+1:mpiendrow;

			if (KX&&KY)
			{
				/* original KX value from groundwater module: if found, we calculate the frozen thickness of the soil in order to transfer this value to the  groundwater model to be accounted for transmissivity estimation */
				FROZENDEPTH = matrix(sr,er, 1, ncols, nrows);
				initialize_grid(gridindex, _frozenSoilInGroundwater_, FROZENDEPTH, NULL, 0, _frozenSoilInGroundwater_, 0, sr,er);
				for (int k = sr; k <= er; k++)
				{
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							FROZENDEPTH[k][j] = hgeo[k][j][0] - hgeo[k][j][levels];
						}
						else
						{
							FROZENDEPTH[k][j] = nodata;
						}
					}
				}
				frozendepthingwindex = gridindex;
				gridindex++;

				if (ist_vorhanden(_colmation_frozen_) && ist_vorhanden(_colmation_unfrozen_))
				{
					GW1COLMFROZEN = gridinfo[getindex(_colmation_frozen_)].pointer;
					GW1COLMUNFROZEN = gridinfo[getindex(_colmation_unfrozen_)].pointer;
					sprintf(tstr, "%s%d", _gw_kolmation_, 1);
					GW1COLMATION = gridinfo[getindex(tstr)].pointer;
				}

			}

			THAWTMOD = matrix(sr,er, 1, ncols, nrows);
			initialize_grid(gridindex, thawdepth_gridname, THAWTMOD, NULL, 0, _ThawDepthTMod_, thawdepth_writegrid, sr,er);
			for (int k = sr; k <= er; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						THAWTMOD[k][j] = 0;
					}
					else
					{
						THAWTMOD[k][j] = nodata ;
					}
				}
			}
			thawtmodindex = gridindex;
			gridindex++;

		}


		if (HDRAIN) /* if a drainage grid was read in (drainage should be considerd) now the layer indizes
			for the drainage networg are calculated and a drainage discharge grid is generated */
		{
			/* drainage discharge */
			QDRAIN = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, QDRAINgridname, QDRAIN, NULL, 0, _Drainageabfluss_, QDRAINwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						QDRAIN[k][j] = 0;
					}
					else
					{
						QDRAIN[k][j] = nodata;
					}
				}
			}
			qdrainindex = gridindex;
			gridindex++;

			/* soil level containing drainage pipes */
			DRAINLEVEL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			DRAINLEVELwritegrid = 0;
			initialize_grid(gridindex, _Drainagelevel_, DRAINLEVEL, NULL, 0, _Drainagelevel_, DRAINLEVELwritegrid, mpistartrow, mpiendrow);
			drainlevelindex = gridindex;

			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						m = 0;
						while ((hgeo[k][j][0] - hgeo[k][j][m]) < HDRAIN[k][j] - 0.01) m++; /* depth: differece between surface (hgeo[j][k][0]) and layer bottom hgeo[j][k][m] */
						/*  */
						DRAINLEVEL[k][j] = m;
					}
				}
			}

			/* allocating memeory for statistic registers for drainage flow */
			zoneQDRAIN = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneQDRAIN) speichererror(1);
			for (int i = 0; i <= zonmax; i++)
			{ zoneQDRAIN[i].n         = zoneQDRAIN[i].RelCellSum = zoneQDRAIN[i].sx        = 0;
				zoneQDRAIN[i].zonenwert = zone[i].zonenwert;
				zoneQDRAIN[i].rel       = zone[i].rel;
			}

			QDRAINstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneQDRAIN,(int)zonmax,(int)QDRAINstatcode,(int)5,(char*)QDRAINstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));       if (retvalue!=0) throw retvalue;
			QDRAINstat->schreibzonenwerte("drainage discharge [mm per zone] (QDRAIN)",instanzname,LastResults,maxLastResults);
		}
		else
		{
			QDRAIN=NULL;
			DRAINLEVEL = NULL;
			zoneQDRAIN = NULL;
			QDRAINstat = NULL;
		}

		if (MAXPOND)
		{
			/* allocating memory for the statistic registers for the content of the pond-storage (as subbasin average) */
			zonePOND = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zonePOND) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{ zonePOND[i].n         = zonePOND[i].RelCellSum = zonePOND[i].sx        = 0;
				zonePOND[i].zonenwert = zone[i].zonenwert;
				zonePOND[i].rel       = zone[i].rel;
			}
			PONDstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zonePOND,(int)zonmax,(int)PONDstatcode,(int)5,(char*)PONDstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));       if (retvalue!=0) throw retvalue;
			PONDstat->schreibzonenwerte("ponding water storage content [mm per zone]  ",instanzname,LastResults,maxLastResults);
		}
		else
		{
			zonePOND=NULL;
			POND=NULL;
			PONDstat=NULL;
		}
		
		initialize_lakes(_lake_model_);

		if (do_ArcVeg)
		{
			// generate internal grids and stacks that may not be read in 
			// AVPGP
			AVPGP = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, AVPGP_name, AVPGP, NULL, 0, _AVPGP_, AVPGP_writegrid, mpistartrow, mpiendrow);
			AVPGPindex = gridindex;
			gridindex++;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						AVPGP[k][j] = 0;
					}
				}
			}

			// AVINSCAL 
			int nuse = MultiLayerLUTable[1].getNumFractions();
			AVINSCAL = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, nuse, nrows);
			AVINSCAL1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, AVINSCAL_stackname, AVINSCAL1, AVINSCAL, nuse, _AVINSCAL_, AVINSCAL_writestack, mpistartrow, mpiendrow);
			AVINSCALindex = gridindex;
			gridindex++;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				//for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						for (int lv = 0; lv <= nuse; lv++)
						{
							AVINSCAL[k][j][lv] = 0;
						}
						AVINSCAL1[k][j] = 0;
					}
					else
					{
						for (int lv = 0; lv <= nuse; lv++)
						{
							AVINSCAL[k][j][lv] = nodata;
						}
						AVINSCAL1[k][j] = nodata;
					}
				}
			}

			// AVNMINTOT
			AVNMINTOT = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, AVNMINTOT_name, AVNMINTOT, NULL, 0, _AVNMINTOT_, AVNMINTOT_writegrid, mpistartrow, mpiendrow);
			AVNMINTOTindex = gridindex;
			gridindex++;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						AVNMINTOT[k][j] = 0;
					}
				}
			}

			// AVNMIN stack 
			AVNMIN = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, nuse, nrows);
			AVNMIN1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, AVNMIN_stackname, AVNMIN1, AVNMIN, nuse, _AVNMIN_, AVNMIN_writestack, mpistartrow, mpiendrow);
			AVNMINindex = gridindex;
			gridindex++;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
					//for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						for (int lv = 0; lv <= 5; lv++)
						{
							AVNMIN[k][j][lv] = 0;
						}
						AVNMIN1[k][j] = 0;
					}
					else
					{
						for (int lv = 0; lv <= 5; lv++)
						{
							AVNMIN[k][j][lv] = nodata;
						}
						AVNMIN1[k][j] = nodata;
					}
				}
			}

		}

		if ((do_ArcVeg)&&(gridslesen!=1)&&(gridslesen != 3))
		{
			double excessicedefaultfactor = 1.0;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if ((lake_active) && LAKE[k][j] != nodata)
					{
						// no ice wedges below a lake allowed
						ICEWEDGEUPPERDEPTH[k][j] = nodata;
						ICEWEDGELOWERDEPTH[k][j] = nodata;
					} 
					// now generating excess ice where required
					else if ((zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k][j] != nodata) && (ICEWEDGELOWERDEPTH[k][j] != nodata))
					{   // at the center lines of an ice wedge with full depth
						double depthsum = 0;
						for (int lv = 1; lv <= levels; lv++)
						{
							double thickness = SOILtable[(int)SOIL[k][j]].getDzui()[lv]; // in m
							depthsum += thickness;
							if (depthsum <= ICEWEDGEUPPERDEPTH[k][j])
							{
								ExcessIceFactor[k][j][lv] = nodata;
							}
							else if (depthsum < ICEWEDGELOWERDEPTH[k][j])
							{
								double excessicethickness = excessicedefaultfactor * thickness; // in m
								ExcessIceFactor[k][j][lv] = excessicethickness * 0.918*1000.0; // in mm WE

								hgeo[k][j][lv] = hgeo[k][j][lv - 1] - excessicethickness - thickness; // new hgeo is by excessicethickness deeper than before (accumulating over all alyers that are within the wedge)
								depthsum += excessicethickness;

								double extrawaterequivalent = ExcessIceFactor[k][j][lv] / 1000.0; //in m						
								double watercontentsat = thickness * SOILtable[(int)SOIL[k][j]].getThetas()[lv] + extrawaterequivalent; // satration water content in mm
								theta[k][j][lv] = watercontentsat / (thickness + excessicethickness);

							}
							else
							{   // below the ice wedge, all additional added thicknesses must also be added (or rather subtracted) from the heodetic height of the layers
								ExcessIceFactor[k][j][lv] = nodata;
								hgeo[k][j][lv] = hgeo[k][j][lv - 1] - thickness; // new hgeo is by excessicethickness deeper than before (accumulating over all alyers that are within the wedge)									
							}
						}
					}
					else
					{   // one cell to the left, right, top and bottom as well as diagonal we add half of the ice wedge center thickness as ice wedge
						int kk = k;
						int jj = j;
						bool onwedge = false;
						if ((k > 1) && (j > 1) && (zn[k - 1][j - 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k - 1][j - 1] != nodata) && (ICEWEDGELOWERDEPTH[k - 1][j - 1] != nodata))
						{
							kk = k - 1;
							jj = j - 1;
							onwedge = true;
						}
						else if ((k > 1) && (zn[k - 1][j] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k - 1][j] != nodata) && (ICEWEDGELOWERDEPTH[k - 1][j] != nodata))
						{
							kk = k - 1;
							jj = j;
							onwedge = true;
						}
						else if ((k > 1) && (j < ncols) && (zn[k - 1][j + 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k - 1][j + 1] != nodata) && (ICEWEDGELOWERDEPTH[k - 1][j + 1] != nodata))
						{
							kk = k - 1;
							jj = j + 1;
							onwedge = true;
						}
						else if ((j > 1) && (zn[k][j - 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k][j - 1] != nodata) && (ICEWEDGELOWERDEPTH[k][j - 1] != nodata))
						{
							kk = k;
							jj = j - 1;
							onwedge = true;
						}
						else if ((j < ncols) && (zn[k][j + 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k][j + 1] != nodata) && (ICEWEDGELOWERDEPTH[k][j + 1] != nodata))
						{
							kk = k;
							jj = j + 1;
							onwedge = true;
						}
						else if ((k < nrows) && (j > 1) && (zn[k + 1][j - 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k + 1][j - 1] != nodata) && (ICEWEDGELOWERDEPTH[k + 1][j - 1] != nodata))
						{
							kk = k + 1;
							jj = j - 1;
							onwedge = true;
						}
						else if ((k < nrows) && (zn[k + 1][j] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k + 1][j] != nodata) && (ICEWEDGELOWERDEPTH[k + 1][j] != nodata))
						{
							kk = k + 1;
							jj = j;
							onwedge = true;
						}
						else if ((k < nrows) && (j < ncols) && (zn[k + 1][j + 1] != nodata) && (zn[k][j] != nodata) && (ICEWEDGEUPPERDEPTH[k + 1][j + 1] != nodata) && (ICEWEDGELOWERDEPTH[k + 1][j + 1] != nodata))
						{
							kk = k + 1;
							jj = j + 1;
							onwedge = true;
						}

						if (onwedge)
						{
							double depthsum = 0;
							for (int lv = 1; lv <= levels; lv++)
							{
								double thickness = SOILtable[(int)SOIL[k][j]].getDzui()[lv]; // in m
								depthsum += thickness;
								if (depthsum <= ICEWEDGEUPPERDEPTH[kk][jj])
								{
									ExcessIceFactor[k][j][lv] = nodata;
								}
								else if (depthsum < ICEWEDGELOWERDEPTH[kk][jj] / 2.0)
								{
									double excessicethickness = excessicedefaultfactor * thickness; // in m
									ExcessIceFactor[k][j][lv] = excessicethickness * 0.918*1000.0; // in mm WE

									hgeo[k][j][lv] = hgeo[k][j][lv - 1] - excessicethickness - thickness; // new hgeo is by excessicethickness deeper than before (accumulating over all alyers that are within the wedge)
									depthsum += excessicethickness;

									double extrawaterequivalent = ExcessIceFactor[k][j][lv] / 1000.0; //in m						
									double watercontentsat = thickness * SOILtable[(int)SOIL[k][j]].getThetas()[lv] + extrawaterequivalent; // satration water content in mm
									theta[k][j][lv] = watercontentsat / (thickness + excessicethickness);

								}
								else
								{   // below the ice wedge, all additional added thicknesses must also be added (or rather subtracted) from the heodetic height of the layers
									ExcessIceFactor[k][j][lv] = nodata;
									hgeo[k][j][lv] = hgeo[k][j][lv - 1] - thickness; // new hgeo is by excessicethickness deeper than before (accumulating over all alyers that are within the wedge)									
								}
							}
						}
					}
				}
			}
			// printstackfile(ExcessIceFactor_stackname, ExcessIceFactor, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
		}


		/* allocating memory for array variables (for 1-D-columns of the 3-D-stacks) */
		/*
a_etr          = (double*)calloc(levels+1,sizeof(double));  // array for real transpiration 
a_evar         = (double*)calloc(levels+1,sizeof(double));  // array for real evaporation 
etrdetail      = (double*)calloc(levels+1,sizeof(double));  // array for real evaptranspiration summed up for all vegetation layers 
a_dicke        = (double*)calloc(levels+1,sizeof(double));  // array for the thickness of the layers 
theta_temp     = (double*)calloc(levels+1,sizeof(double));  // array for theta-values of the layers (for temporal storing) 
h_temp1        = (double*)calloc(levels+1,sizeof(double));  // array for hydraulic head values of the last timestep (for repeated use in iteration) 
theta_average  = (double*)calloc(levels+1,sizeof(double));  // array for averaged theta-values of the layers (for temporal storing during iteration, step 1..i-1) 
theta_sl_old   = (double*)calloc(levels+1,sizeof(double));  // array for theta-values at the begin of the iteration 
q_temp         = (double*)calloc(levels+1,sizeof(double));  // array for fluxes 
extractionweight = (double*)calloc(levels+1,sizeof(double));  // extraction density for the layers within root zone for evapotranspiration 
qlat           = (double*)calloc(levels+1,sizeof(double));  // array for baseflow per level (traditional, i.e. if gw-module is not active) 
gwin           = (double*)calloc(levels+1,sizeof(double));  // array for re-infiltration per level 
gwex           = (double*)calloc(levels+1,sizeof(double));  // array for exfiltration into rivers (baseflow if gwflow_active==1) per level 
interflow      = (double*)calloc(levels+1,sizeof(double));  // interflow per level 
lateralinflow  = (double*)calloc(levels+1,sizeof(double));  // lateral inflow (calculated as balance from groundwater module) per level 
*/
		//    sicker=0;    /* fraction of infiltration excess, re-infiltrating into the soil  */

		/* now some memory buffer for the control point output files is allocated */


		if (control_row>0) /* writing the files only, if control_row is greater than zero */
		{
			if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
			{
				/* file buffer for q-bottom fluxes */
				qubufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
				if (qubufmax<1) qubufmax=1;
				qubufcnt      = 0; /* initialising the buffer index counter with 0 */
				qubuffer      = (buffertype*)calloc(qubufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=qubufmax;bufrow++) qubuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
				/* allocating memory for each level the zero-index-element */
				/* file buffer for the theta values of all layers at */
				thetabufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
				if (thetabufmax<1) thetabufmax=1;
				thetabufcnt      = 0; /* initialising the buffer index counter with 0 */
				thetabuffer      = (buffertype*)calloc(thetabufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=thetabufmax;bufrow++) thetabuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
				/* allocating memory for each level and the mean + the zero-index-element */

				/* file buffer for the theta values of all layers at */
				hhydbufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
				if (hhydbufmax<1) hhydbufmax=1;
				hhydbufcnt      = 0; /* initialising the buffer index counter with 0 */
				hhydbuffer      = (buffertype*)calloc(hhydbufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=hhydbufmax;bufrow++) hhydbuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
				/* allocating memory for each level and the mean + the zero-index-element */

				/* file buffer for some other data like the crop parameters */
				otherdatabufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(16+2+this->iNumVegetationLayers*6+ 1)*sizeof(double))); // 16 variables will be stored at the moment - may be more in later releases --> new in 8.09.00: for each vegetation layer, LAI, VCF, RSC, RSI, Z0 and ROOT will be written as well as RSE and ALBEDO for all layers once (in common)
				if (otherdatabufmax<1) otherdatabufmax=1;
				otherdatabufcnt      = 0; /* initialising the buffer index counter with 0 */
				otherdatabuffer      = (buffertype*)calloc(otherdatabufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=otherdatabufmax;bufrow++) otherdatabuffer[bufrow].wert=(double*)calloc(16+2+this->iNumVegetationLayers*6+1,sizeof(double));
				/* allocating memory for each level and the mean + the zero-index-element */

				/* file buffer for the etr values of all layers at */
				etrdetailbufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
				if (etrdetailbufmax<1) etrdetailbufmax=1;
				etrdetailbufcnt      = 0; /* initialising the buffer index counter with 0 */
				etrdetailbuffer      = (buffertype*)calloc(etrdetailbufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=etrdetailbufmax;bufrow++) etrdetailbuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
				/* allocating memory for each level and the mean + the zero-index-element */

				/* file buffer for the etr values of all layers at */
				interflowbufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
				if (interflowbufmax<1) interflowbufmax=1;
				interflowbufcnt      = 0; /* initialising the buffer index counter with 0 */
				interflowbuffer      = (buffertype*)calloc(interflowbufmax+1,sizeof(buffertype));
				for (bufrow=1;bufrow<=interflowbufmax;bufrow++) interflowbuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
				/* allocating memory for each level and the mean + the zero-index-element */

				/* file buffer for the soil temperature values of all layers at */
				if (do_TMod>0)
				{
					temperbufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
					if (temperbufmax<1) temperbufmax=1;
					temperbufcnt      = 0; /* initialising the buffer index counter with 0 */
					temperbuffer      = (buffertype*)calloc(temperbufmax+1,sizeof(buffertype));
					for (bufrow=1;bufrow<=temperbufmax;bufrow++) temperbuffer[bufrow].wert=(double*)calloc(levels+1,sizeof(double));
					/* allocating memory for each level and the mean + the zero-index-element */
				}

				if ((hhyddatei=fopen(hhyddateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for hydraulic head (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(hhyddatei,"YY\tMM\tDD\tHH\thhyd[%d][%d]",control_row,control_col);
				for (int i=1;i<=levels;i++) fprintf(hhyddatei,"[%d]\t",i);
				fprintf(hhyddatei,"\n");
				fclose(hhyddatei);

				// initializing the file header for otherdata (some control point water balance data)
				if ((otherdatadatei=fopen(otherdatadateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for other detail data (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(otherdatadatei,"YY\tMM\tDD\tHH\totherdata[%d][%d]",control_row,control_col);
				fprintf(otherdatadatei,"QD\tQI\tDELTAGW/QB\tQINF\tGWIN\tGWEX\tGWN\tUPRISE\tPERCOL\tMACR\tQDRAIN\tACT_POND\tINFEX\tGWST\tGWLEVEL\tROOT");
				fprintf(otherdatadatei,"\tAlbedo\tRSE");
				for (int i=1;i<=this->iNumVegetationLayers;i++)
				{
					fprintf(otherdatadatei,"\tLAI_%d\tVCF_%d\tRSC_%d\tRSI_%d\tROOT_%d\tZ0_%d",i,i,i,i,i,i);
				}
				fprintf(otherdatadatei,"\n");
				fclose(otherdatadatei);

				if ((thetadatei=fopen(thetadateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for soil theta (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(thetadatei,"YY\tMM\tDD\tHH\ttheta[%d][%d]",control_row,control_col);
				for (int i=1;i<=levels;i++) fprintf(thetadatei,"[%d]\t",i);
				fprintf(thetadatei,"\n");
				fclose(thetadatei);

				if ((qudatei=fopen(qudateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for hydraulic head (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(qudatei,"YY\tMM\tDD\tHH\tq_bottom[%d][%d]",control_row,control_col);
				for (int i=1;i<=levels;i++) fprintf(qudatei,"[%d]\t",i);
				fprintf(qudatei,"\n");
				fclose(qudatei);

				if ((etrdetaildatei=fopen(etrdetaildateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for transpiration/evaporation per layer (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(etrdetaildatei,"YY\tMM\tDD\tHH\tetr[%d][%d](etr[0] is evaporation from first layer)",control_row,control_col);
				// for etr: etr[0] will be evaporation from first soil layer
				for (int i=0;i<=levels;i++) fprintf(etrdetaildatei,"[%d]\t",i);
				fprintf(etrdetaildatei,"\n");
				fclose(etrdetaildatei);

				if ((interflowdatei=fopen(interflowdateiname,"wt"))==NULL)
				{
					(*WSstream) << "FATAL: error opening file for interflow (single point) \n"
					"-> check disk space or user rights"<<"[end-entry]";
					throw(62);
				}
				fprintf(interflowdatei,"YY\tMM\tDD\tHH\tinterflow_per_layer[%d][%d]",control_row,control_col);
				for (int i=1;i<=levels;i++) fprintf(interflowdatei,"[%d]\t",i);
				fprintf(interflowdatei,"\n");
				fclose(interflowdatei);

				if (do_TMod>0)
				{
					if ((temperfile=fopen(temperfilename,"wt"))==NULL)
					{
						(*WSstream) << "FATAL: error opening file for soil temperatures (single point details) \n"
						"-> check disk space or user rights"<<"[end-entry]";
						throw(62);
					}
					fprintf(temperfile,"YY\tMM\tDD\tHH\ttemper[%d][%d]",control_row,control_col);
					for (int i=1;i<=levels;i++) fprintf(temperfile,"[%d]\t",i);
					fprintf(temperfile,"\n");
					fclose(temperfile);
				}
			}
		}

#ifdef _WASIMDEBUG_
		char iterfilename[1024];
		sprintf(iterfilename,"iter_%d.log",mpi_rank);
		if ((tempdatei=fopen(iterfilename,"a"))==NULL)
		{
			(*WSstream) << "FATAL: error opening log-file for iteration counts (debug-version only)"<<"[end-entry]";
			throw(63);
		}
		fprintf(tempdatei,"month\tday\thour\ttheta1\ttheta2\ttheta3\tqu1\tinter1\tqu2\tinter2\tqu3\tinter3\thh1\thh2\thh3\tdeltah1\tdeltah2\tthet_av1\tthet_av2\tthet_av3\n");
#endif

		/* some variables required for tracer mixing */
		/*
maxfromlr     = (double*)calloc(lv+1,sizeof(double));
qufromlr      = (double*)calloc(lv+1,sizeof(double));
iflfromlr     = (double*)calloc(lv+1,sizeof(double));
qlatfromlr    = (double*)calloc(lv+1,sizeof(double));
gwexfromlr    = (double*)calloc(lv+1,sizeof(double));
ratio         = (double*)calloc(lv+1,sizeof(double));

macrocapacity = (double*)calloc(lv+1,sizeof(double));
*/
		additionals=newstring(30,1); /* contains display output information for presence of IRRIG, CLAY, HDRAIN, DRAINDIST, MAXPOND */
		sprintf(additionals,"%s%s%s%s%s%s",FZT?"+ft,":"",IRRIG?"+ir,":"",surface_interaction>0?(surface_interaction>1?"+srf2,":"+srf1"):"",HDRAIN?"+dr,":"",CLAY?"+cl,":"",MAXPOND?"+po,":"");

		/* grid and stack names are not needed further bacause the names are stored (and memory was allocated, too) in gridinfo-array */
		if (SB_1_gridname   ) free(SB_1_gridname   );   SB_1_gridname = NULL;
		if (SB_2_gridname   ) free(SB_2_gridname   );   SB_2_gridname = NULL;
		if (ROOTgridname    ) free(ROOTgridname    );    ROOTgridname = NULL;
		if (theta_stackname ) free(theta_stackname ); theta_stackname = NULL;
		if (hh_stackname    ) free(hh_stackname    );    hh_stackname = NULL;
		if (hgeo_stackname  ) free(hgeo_stackname  );  hgeo_stackname = NULL;
		if (qu_stackname    ) free(qu_stackname    );    qu_stackname = NULL;
		if (GWSTANDgridname ) free(GWSTANDgridname ); GWSTANDgridname = NULL;
		if (GWthetagridname ) free(GWthetagridname ); GWthetagridname = NULL;
		if (GWNgridname     ) free(GWNgridname     );     GWNgridname = NULL;
		if (UPRISEgridname  ) free(UPRISEgridname  );  UPRISEgridname = NULL;
		if (PERCOLgridname  ) free(PERCOLgridname  );  PERCOLgridname = NULL;
		if (PONDgridname    ) free(PONDgridname    );    PONDgridname = NULL;
		if (GWLEVELgridname ) free(GWLEVELgridname ); GWLEVELgridname = NULL;
		if (QDRAINgridname  ) free(QDRAINgridname  );  QDRAINgridname = NULL;
		if (QDgridname      ) free(QDgridname      );      QDgridname = NULL;
		if (QIgridname      ) free(QIgridname      );      QIgridname = NULL;
		if (DELTAGWname     ) free(DELTAGWname      );      DELTAGWname = NULL;
		if (INFEXgridname   ) free(INFEXgridname   );   INFEXgridname = NULL;
		if (SATTgridname    ) free(SATTgridname    );    SATTgridname = NULL;
		if (MACROgridname   ) free(MACROgridname   );   MACROgridname = NULL;
		if (IRRIGgridname   ) free(IRRIGgridname   );   IRRIGgridname = NULL;
		if (methode>2)
		{
			if (SUBSTEPSgridname   ) free(SUBSTEPSgridname   );   SUBSTEPSgridname = NULL;
		}

		if (do_TMod>0)
		{
			if (temper_stackname) free(temper_stackname); temper_stackname = NULL;
			if (se_stackname) free(se_stackname); se_stackname = NULL;
			if (thawdepth_gridname) free(thawdepth_gridname); thawdepth_gridname = NULL;

			// set default soil properties if they are not set in the soil table yet
			for (int i = 1; i <= ST[0].getCode(); i++) /* ST stands for soil table. In ST[0].code the number of soil types is stored */
			{
				for (int k = 1; k <= ST[i].getLevels(); k++)
				{
					if (ST[i].getDryDensity(k)<0)
					{
						ST[i].setDryDensity(k, dTModRhoDrySoil);
					}

					if (ST[i].getDryHeatCapacity(k)<0)
					{
						ST[i].setDryHeatCapacity(k, dTModCDrySoil);
					}

					if (ST[i].getDryThermalConduct(k)<0)
					{
						ST[i].setDryThermalConduct(k, dTModLamDrySoil);
					}
				}

				if (ST[i].getKMinFrozenSoil()<0)
				{
					ST[i].setKMinFrozenSoil(dTModKmin);
				}

				for (int k = 1; k <= ST[i].getLevels(); k++)
				{
					ST[i].LookUpTableGenerationTMod(k, this->dTModFactor);
				}
			}

		}



		/*
balance=(double*)calloc(11,sizeof(double));
maxbalance=(double*)calloc(11,sizeof(double));
for (int i=0;i<=9;i++) 
{
balance[i]=0.0;
maxbalance[i]=0.0;
}
*/
		/* looking for ICE grid */
		if (ist_vorhanden(_ice_firn_)||ist_vorhanden(_GlacierMassBalance_))
		{
			if (ist_vorhanden(_GlacierizedCells_) && ist_vorhanden(_GlacierMassBalance_))
			{
				ICE=gridinfo[getindex(_GlacierizedCells_)].pointer; // dynamic glacier model: ICE points to another grid possibly
				ICEOLD=gridinfo[getindex(_GlacierizedCellsOld_)].pointer; // for adding immobilized/mobilized water to an extra statfile for balance purposes
				/* initializing a statistc register set for soil moisturte [entire unsaturated and saturated zone] */
				zoneSB_2_IMMO = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
				if (!zoneSB_2_IMMO) speichererror(1);
				for (int i=0;i<=zonmax;i++)
				{   zoneSB_2_IMMO[i].n         = zoneSB_2_IMMO[i].RelCellSum = zoneSB_2_[i].sx        = 0;
					zoneSB_2_IMMO[i].zonenwert = zone[i].zonenwert;
					zoneSB_2_IMMO[i].rel       = zone[i].rel;
				}
				SB_2_IMMO_statcode = SB_2_statcode;
				SB_2_IMMO_statfilename   = newstring(_stringsize_+1,1);
				std::string tmp;
				tmp=std::string(SB_2_statfilename);
				tmp=tmp.substr(0,tmp.find_last_of("."))+"_immo"+tmp.substr(tmp.find_last_of("."));
				sprintf(SB_2_IMMO_statfilename,"%s",tmp.c_str());
				SB_2_IMMO_stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneSB_2_IMMO,(int)zonmax,(int)SB_2_IMMO_statcode,(int)5,(char*)SB_2_IMMO_statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));     if (retvalue!=0) throw retvalue;
				SB_2_IMMO_stat->schreibzonenwerte("immobilised or mobilised soil water by glacier resizing 0..profile depth per zone ",instanzname,LastResults,maxLastResults);
			}
			else
			{
				ICE=gridinfo[getindex(_ice_firn_)].pointer;
			}
			glacierlist=(glacierlisttype*)calloc(zonmax+1,sizeof(glacierlisttype));
			for (int i=0;i<=zonmax;i++)
			{
				glacierlist[i].isGlacierized = 0;
				glacierlist[i].CellCntBare = 0;
				glacierlist[i].CellCntGlacerized = 0;
			}

			if (mpi_rank == 0)
			{
				for (int i = 1; i <= nrows; i++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						if (zn[i][j] != nodata)
						{
							if (ICE[i][j] > 0)
							{
								glacierlist[zarray[(int)zn[i][j]]].isGlacierized = 1;
								glacierlist[zarray[(int)zn[i][j]]].CellCntGlacerized += ICE[i][j];
								glacierlist[zarray[(int)zn[i][j]]].CellCntBare += (1.0 - ICE[i][j]);
							}
							else
							{
								glacierlist[zarray[(int)zn[i][j]]].CellCntBare++;
							}
						}
					}
				}
			}
#ifdef _MPI
			// bcastig center_of_mass to all other ranks
			MPI_Barrier(MPI_COMM_WORLD);
			const int nitems = 3;
			int blocklengths[3] = { 1, 1, 1 };
			MPI_Datatype types[3] = { MPI_C_BOOL, MPI_DOUBLE, MPI_DOUBLE };
			MPI_Aint offsets[3] = { offsetof(glacierlisttype, isGlacierized),
				offsetof(glacierlisttype, CellCntGlacerized),
				offsetof(glacierlisttype, CellCntBare) };

			MPI_Datatype MPIglacierlisttype;
			MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPIglacierlisttype);
			MPI_Type_commit(&MPIglacierlisttype);

			if (mpi_rank > 0)
			{
				glacierlist = (glacierlisttype*)calloc(zonmax + 1, sizeof(glacierlisttype));
				if (!glacierlist) speichererror(1);
			}
			MPI_Bcast(&(glacierlist[0]), zonmax + 1, MPIglacierlisttype, 0,MPI_COMM_WORLD);

			MPI_Type_free(&MPIglacierlisttype);
#endif
		}
		else
		{
			ICE=NULL;
			glacierlist=NULL;
		}


		double deltasb=0;
		double deltapond=0;

		if (control_row>0)
		{
			if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
			{
				dicken=ST[(int)SOIL[control_row][control_col]].getDzui();
				for (int k=1;k<=levels;k++)
				deltasb+=(theta[control_row][control_col][k]*dicken[k]*1000.0);
				char logfilename[1024];
				sprintf(logfilename,"modell_%d.log",mpi_rank);
				logdatei=fopen(logfilename,"a+");
				fprintf(logdatei,"soil water content (for control point r=%d c=%d) at model-start: %lf mm\n",control_row,control_col,deltasb);
				fclose(logdatei);
			}
		}
		double sum_sb1=0.0;
		double anz_sb1=0.0;
		double sum_sb2=0.0;
		double anz_sb2=0.0;
		double *thetasat;
		for (int i=mpistartrow;i<=mpiendrow;i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if ((zn[i][j]>0)&&(MPIREG[i][j]==mpi_rank))
				{
					int s=(int)SOIL[i][j];
					dicken=ST[s].getDzui();
					thetasat = ST[s].getThetas();
					for (int k=1;k<=levels;k++)
					{
						sum_sb1+=(theta[i][j][k]*dicken[k])*(1.0-icevalue);
						sum_sb2+=(thetasat[k]*dicken[k])*icevalue;
					}
					anz_sb1+=RelCellArea!=NULL?RelCellArea[i][j]*(1.0-icevalue):(1.0-icevalue);
					anz_sb1+=icevalue;
					anz_sb2+=RelCellArea!=NULL?RelCellArea[i][j]*icevalue:icevalue;
					anz_sb2+=(1.0-icevalue);
				}
			}
		}
		double sumtemp1=sum_sb1;
		double anztemp1=anz_sb1;
		double sumtemp2=sum_sb2;
		double anztemp2=anz_sb2;
#ifdef _MPI
		// reduce anz_sb and sum_sb from all instances
		MPI_Reduce(&anz_sb1, &anztemp1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&sum_sb1, &sumtemp1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&anz_sb2, &anztemp2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&sum_sb2, &sumtemp2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
		if (mpi_rank==0)
		{
			char logfilename[1024];
			sprintf(logfilename,"modell_%d.log",mpi_rank);
			logdatei=fopen(logfilename,"a+");
			fprintf(logdatei,"mean soil water content for all cells at model-start: %lf m^3/m^2\n",sumtemp1/anztemp1+sumtemp2/anztemp2);
			fprintf(logdatei,"mean soil water content for unglacierized or only partly unglacierized cells at model-start: %lf m^3/m^2\n",sumtemp1/anztemp1);
			fprintf(logdatei,"mean soil water content for partly or completely glacierized cells at model-start: %lf m^3/m^2\n",sumtemp2/anztemp2);
			fclose(logdatei);
		}

		if (doSpinUp!=0 && HRU==NULL)
		{
			HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
			HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
		}

		// if the snow model should be part of the soil model (i.e. heat transfer and/or water fluxes through snow should 
		// be calculated with the algorithms already implemented in unsatzon model, we need to read those parameters. Default is "off"
		if (readLayeredSnowModelParameters() != 0)
		{
			readLayeredSnowProperties();
			ReadOrCreateLayeredSnowGridsAndStacks();
		}

		createFieldCapacityGrids();
		createSB1RELGrid();

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
//#pragma optimize ("p",off)

void unsatzonmodellklasse::createFieldCapacityGrids()
{

	if ((indexFCmin_mm = ist_vorhanden(_minFieldCapacity_mm_) > 0) &&
		(indexFCmax_mm = ist_vorhanden(_maxFieldCapacity_mm_) > 0) &&
		(indexFCact_mm = ist_vorhanden(_usableFieldCapacity_mm_) > 0) &&
		(indexFCact_perc = ist_vorhanden(_usableFieldCapacity_perc_) > 0))
	{
		indexFCmin_mm   = getindex(_minFieldCapacity_mm_);
		indexFCmax_mm   = getindex(_maxFieldCapacity_mm_);
		indexFCact_mm   = getindex(_usableFieldCapacity_mm_);
		indexFCact_perc = getindex(_usableFieldCapacity_perc_);
		FCMIN = gridinfo[indexFCmin_mm].pointer;
		FCMAX = gridinfo[indexFCmax_mm].pointer;
		FCACT = gridinfo[indexFCact_mm].pointer; 
		FCACTPERC = gridinfo[indexFCact_perc].pointer;

		doFCcalc = true;
		// call updateFieldCapacityGrids() for the first time, only to be sure the initialization is done 
		updateFieldCapacityGrids();
	}

}

void unsatzonmodellklasse::updateFieldCapacityGrids()
{
	if (doFCcalc)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					int s = (int)SOIL[i][j]; // get soil type
					double *dicken = ST[s].getDzui();
					double *thetar = ST[s].getThetar();
					double *thetas = ST[s].getThetas();
					double rootdepth = WURZ[1][i][j];
					for (int k = 2; k <= iNumVegetationLayers; k++)
					{
						if (rootdepth < WURZ[k][i][j])
						{
							rootdepth = WURZ[k][i][j];
						}
					}

					// setting FCMIN as result of GetMeanTheta with theta_res
					FCMIN[i][j] = GetMeanTheta(s, dicken, thetar, 0, rootdepth, 1)*1000.0; // result in mm, converted via 1000.0 from m

					// setting FCACT as result of GetMeanTheta with theta_sl minus FCMIN
					FCACT[i][j] = GetMeanTheta(s, dicken, theta[i][j], 0, rootdepth, 1)*1000.0; // result in mm, converted via 1000.0 from m
					FCACT[i][j] -= FCMIN[i][j]; 
					if (FCACT[i][j] < 0)
					{
						FCACT[i][j] = 0; // if suction is below wilting point, i.e. soil moisture below residual water content -> happens by evaporation
					}

					// getting thetas for the suction as given in soil table for each layer -> store in theta_fc
					double *theta_fc = new double[levels+1];
					double *psi_FC = ST[s].getPsi_FC();
					for (int k = 1; k <= levels; k++)
					{
						theta_fc[k] = ST[s].get_theta(k, psi_FC[k]/100.0, thetar[k],thetas[k]); // converting psi_fc from cm or hPa as given in control file into m as needed in GetMeanTheta
					}
					// getting FCMAX as result of GetMeanTheta with theta_fc minus FCMIN
					FCMAX[i][j] = GetMeanTheta(s, dicken, theta_fc, 0, rootdepth, 1) * 1000; // converting in m
					FCMAX[i][j] -= FCMIN[i][j];
					if (FCMAX[i][j] < 0)
					{
						FCMAX[i][j] = 0;
					}
					// calculating the percentage of FCACT with regard to FCMAX. If > 100%, then set to 100%
					FCACTPERC[i][j] = FCMAX[i][j] > 0 ? FCACT[i][j] / FCMAX[i][j] * 100.0 : 0;
					if (FCACTPERC[i][j] > 100.0)
					{
						FCACTPERC[i][j] = 100.0;
					}

					delete[] theta_fc;
				}
				else
				{
					if (MPIREG[i][j] == mpi_rank)
					{
						// setting all the grids to nodata
						FCMIN[i][j] = nodata;
						FCACT[i][j] = nodata;
						FCMAX[i][j] = nodata;
						FCACTPERC[i][j] = nodata;
					}
				}
			}
		}
	}
}


void unsatzonmodellklasse::createSB1RELGrid()
{
	if (sb_1_rel_index = ist_vorhanden(_SB_1_REL_grid_) > 0)
	{
		sb_1_rel_index = getindex(_SB_1_REL_grid_);
		SB_1_REL = gridinfo[sb_1_rel_index].pointer;

		doSB1RELcalc = true;
		// call updateSB1RELGrid() for the first time, only to be sure the initialization is done 
		updateSB1RELGrid();
	}
}

void unsatzonmodellklasse::updateSB1RELGrid()
{
	if (doSB1RELcalc && SB_1_)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					int s = (int)SOIL[i][j]; // get soil type
					double *dicken = ST[s].getDzui();
					double *thetar = ST[s].getThetar();
					double *thetas = ST[s].getThetas();
					double rootdepth = WURZ[1][i][j];

					for (int k = 2; k <= iNumVegetationLayers; k++)
					{
						if (rootdepth < WURZ[k][i][j])
						{
							rootdepth = WURZ[k][i][j];
						}
					}

					// setting FCMIN as result of GetMeanTheta with theta_res
					if (rootdepth > 0)
					{
						double temp = GetMeanTheta(s, dicken, thetas, 0, rootdepth, 1) / rootdepth; // result in mm, converted via 1000.0 from m
						if (temp >= 0)
						{
							SB_1_REL[i][j] = SB_1_[i][j] / temp;
						}
						else
						{
							SB_1_REL[i][j] = 1.0;
						}
					}
					else
					{
						SB_1_REL[i][j] = 1.0;
					}

				}
				else
				{
					if (MPIREG[i][j] == mpi_rank)
					{
						SB_1_REL[i][j] = nodata;
					}
				}
			}
		}
	}
}


void unsatzonmodellklasse::setSubStepCnt(int step)
{
	substepcnt=step;
}

char *unsatzonmodellklasse::GetStorageFileName()
{
	if (do_it) return storagefilename;
	else return NULL;
}


void unsatzonmodellklasse::StorageFileError(const char *message,double c1,double c2)
{
	(*WSstream) << "FATAL: error in storage file content reading "<<message<<":\n"
	"value in actual control is:   "<<c1<<"\n"
	"but value in storage file is: "<<c2<<"\n"
	"WaSiM throws..."<<"[end-entry]";
	throw(71);
}

void unsatzonmodellklasse::ReadFloatFromStoragefile(double *FValue,const char *ErrorMessage)
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

void unsatzonmodellklasse::ReadSlrContents()
{
	int i,j,z,k;
	int st_int,st_int2;
	char *satz, *loctsatz;
	int syntaxVersion = 2; // default syntax expected is for version 2 of the storage file
	double origTimestep = InV;
	double dMultipleOfOrigTimestep = 1.0;
	double dCumFracOfOrigTimestep = 0.0;
	int sourceMaxzones;
	double ndr, sxdr, reldr, RelCellSumdr, zonenwertdr; // used for temporary reading (d=direct flow r=read from file)
	double nir, sxir, relir, RelCellSumir, zonenwertir; // used for temporary reading (i=direct flow r=read from file)

	satz=(char*)calloc(_stringsize_,sizeof(char));
	loctsatz = satz; // store actual address to be able to release the memory later

	char *completedname = new char[_stringsize_+1];
	sprintf(completedname,"%s",storagefilename);
	ExtendFileName(DefaultInputDirectory,completedname);
	if ((storagefile = fopen(completedname,"rt"))==NULL)
	{
		(*WSstream) << "FATAL: error opening storage content file "<<completedname<<"[end-entry]";
		throw(64);
	}

	// now looking for the section name.
	// new version >= 10.06.00 that could deal with different time steps uses [unsatzon_model_v2], old version uses [unsatzon_model]
	// if new version is found, it is used. If not, old version is searched for and no flexible time step can be used (backwards compatible)
	char cmpString[100];
	char* tmp = satz;
	sprintf(cmpString, "[unsatzon_model_v2]"); 
	do
	{  
		satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
		if (feof(storagefile))
		{
			satz = tmp; // because fgets sets satz to NULL when feof has been reached, we reset the pointer to it's original value
			(*WSstream) << "INFORMATION: couldn't find section [unsatzon_model_v2], looking now for [unsatzon_model] in " << completedname << "[end-entry]";
			syntaxVersion = 1;
			fseek(storagefile, 0, 0); // reset storagefile
			sprintf(cmpString, "[unsatzon_model]");
			do
			{
				satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (feof(storagefile))
				{
					(*WSstream) << "FATAL: error locking for the section [unsatzon_model_v2] or [unsatzon_model] in " << completedname << "[end-entry]";
					throw(71);
				}
			} 
			while (!strstr(satz, cmpString) && !feof(storagefile));
		}
	}
	while (!strstr(satz, cmpString) && !feof(storagefile));
	delete[] completedname;


	/* TG_number */
	fscanf(storagefile,"%d",&st_int); satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
	if (st_int!=zonmax) StorageFileError("number of subbasins",(double)zonmax,(double)st_int);

	if (syntaxVersion >= 2)
	{
		/* timestep used to create that storagefile */
		fscanf(storagefile, "%lf", &origTimestep); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
		dMultipleOfOrigTimestep = origTimestep / InV;
	}

	/* reading subbasin codes */
	for (int j=1;j<=zonmax;j++)
	{ fscanf(storagefile,"%d",&st_int);
		if (st_int!=code[j]) StorageFileError("subbasin code",(double)code[j],(double)st_int);
	} satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading kelsqd */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(kelsqd[j]),"m");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading kelsqi */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(kelsqi[j]),"kelsqi");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading flussdichte */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(flussdichte[j]),"flussdichte");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading kelsqlat */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(kelsqlat[j]),"kelsqlat");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading Q0 */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(Q0[j]),"Q0");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* reading SR_Anteil */
	for (int j=1;j<=zonmax;j++) ReadFloatFromStoragefile(&(SR_Anteil[j]),"SR_Anteil");
	satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string

	/* number of tracers */
	fscanf(storagefile,"%d",&st_int); satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
	if (st_int!=tracer_count) StorageFileError("number of tracers",(double)tracer_count,(double)st_int);

	for (z=0;z<=tracer_count;z++)
	{
		/* actual tracer (0=runoff) */
		fscanf(storagefile,"%d",&st_int); satz=fgets(satz,_stringsize_,storagefile); // read the comment into a dummy string
		if (st_int!=z) StorageFileError("actual tracers(0=runoff)",(double)z,(double)st_int);
		
		if (origTimestep == InV)
		{ // old but still used version, when no time step change occured
			for (int i = 1; i <= zonmax; i++)
			{
				/* actual tracer (0=runoff) */
				fscanf(storagefile, "%d%d", &st_int, &st_int2); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (st_int != i) StorageFileError("actual zone", (double)i, (double)st_int);
				if (st_int2 != (int)maxzeit[i])
				{
					StorageFileError("maxzeit of actual subbasin", (double)((int)maxzeit[i]), (double)st_int);
				}

				for (int j = 1; j <= (int)maxzeit[i]; j++)
				{
					k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &(zoneFZT_QD[z][i][j].n), &(zoneFZT_QD[z][i][j].sx),
						&(zoneFZT_QD[z][i][j].zonenwert), &(zoneFZT_QD[z][i][j].rel), &(zoneFZT_QD[z][i][j].RelCellSum));
					if (k != 5) StorageFileError("number of fields for zonal statistics QD", 5.0, (double)k);
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					if ((k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &(zoneFZT_QI[z][i][j].n), &(zoneFZT_QI[z][i][j].sx),
						&(zoneFZT_QI[z][i][j].zonenwert), &(zoneFZT_QI[z][i][j].rel), &(zoneFZT_QI[z][i][j].RelCellSum)))
						!= 5) StorageFileError("number of fields for zonal statistics QI", 5.0, (double)k);
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				}
				fscanf(storagefile, "%lf", &(QDd[z][i])); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				fscanf(storagefile, "%lf", &(QIi[z][i])); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
			}
		}
		else
		{ // new version that recomputes storage contents based on different number of zones
			for (int i = 1; i <= zonmax; i++) // zonmax is the number of sub basins of the actual configuration. This must still be iodentical even when the time step changed
			{
				/* actual tracer (0=runoff) */
				fscanf(storagefile, "%d%d", &st_int, &st_int2); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (st_int != i) StorageFileError("actual zone", (double)i, (double)st_int);

				// sourceMaxzones will now be compared to actual maxzeit[i]
				if (st_int2 != (int)maxzeit[i])
				{
					//StorageFileError("maxzeit of actual subbasin", (double)((int)maxzeit[i]), (double)st_int);
					sourceMaxzones = st_int2;
				}

				int actMaxzones = (int)maxzeit[i]; // for each subbasin (index for subbasin is i)

				if (actMaxzones < sourceMaxzones)
				{ // initializing a timestep that is larger than the timestep from the sourcefile, i.e. there are multiple source zones that must be added up into one actual zone
					// since all values are read anew from the file and will be added to the existing ones (in case of cumulating multiple zones from the source file into one actual zone), we set all original values to 0 now
					for (int j = 1; j <= actMaxzones; j++)
					{
						zoneFZT_QD[z][i][j].n = 0;
						zoneFZT_QD[z][i][j].sx = 0;
						//zoneFZT_QD[z][i][j].rel = 0;
						zoneFZT_QD[z][i][j].RelCellSum = 0;
						//zoneFZT_QD[z][i][j].zonenwert = 0;
						zoneFZT_QI[z][i][j].n = 0;
						zoneFZT_QI[z][i][j].sx = 0;
						//zoneFZT_QI[z][i][j].rel = 0;
						zoneFZT_QI[z][i][j].RelCellSum = 0;
						//zoneFZT_QI[z][i][j].zonenwert = 0;
					}
					// sourceMaxzones -> number of lines to be read from source files, i.e. as many zones there were in the original time step version, e.g. 47 in 10min, but 8 in 1h resolution only
					// actMaxzones -> number of zones to be populated from the sourceMaxzones zones from the storagefile 
					dCumFracOfOrigTimestep = 0; //index for output is +1 or the ceil function is used to compute the actual output index
					int iSourceZoneCntr = 0;
					while ((ceil(dCumFracOfOrigTimestep - 0.0001) <= actMaxzones)&&(iSourceZoneCntr < sourceMaxzones)) // -0.0001 should make sure that numerical rouding errors so not accumulate too much, so even the last zone from the source file can match the last actual zone 
					// when we have 47 input values but onla 8 output values, we will step 47 times throug the loop but add the index for the ouput only 7 times
					//for (int j = 1; j <= (int)maxzeit[i]; j++)
					{

						dCumFracOfOrigTimestep += dMultipleOfOrigTimestep;
						// read the first content from the storage file
						k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &ndr, &sxdr, &zonenwertdr, &reldr, &RelCellSumdr);
						//					k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &(zoneFZT_QD[z][i][j].n), &(zoneFZT_QD[z][i][j].sx),
						//						&(zoneFZT_QD[z][i][j].zonenwert), &(zoneFZT_QD[z][i][j].rel), &(zoneFZT_QD[z][i][j].RelCellSum));
						if (k != 5) {
							StorageFileError("number of fields for zonal statistics QD", 5.0, (double)k);
						}
						satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string

						k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &nir, &sxir, &zonenwertir, &relir, &RelCellSumir);
						//					k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &(zoneFZT_QI[z][i][j].n), &(zoneFZT_QI[z][i][j].sx),
						//						&(zoneFZT_QI[z][i][j].zonenwert), &(zoneFZT_QI[z][i][j].rel), &(zoneFZT_QI[z][i][j].RelCellSum));
						if (k != 5) {
							StorageFileError("number of fields for zonal statistics QI", 5.0, (double)k);
						}
						satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string

						// now populating the actual records depending on their position in original and actual records arrays
						// sample data: when original time step is 10min, actual time step is 60 min
						// dMultipleOfOrigTimestep 0.16666
						// dCumFracOfOrigTimestep = 0.0; after 6 sub steps we have reached a full new time step
						// --> will reach at maximum actMaxzones (maxzeit[i])
						// iterative up-counting will give the actual index orf the array to add the read values to

						double corrfact;
						corrfact = (dMultipleOfOrigTimestep < 1.0) ? 1.0 : dMultipleOfOrigTimestep;

						int j = (int)ceil(dCumFracOfOrigTimestep-0.0001); // to which local index we should add the values

						if (zoneFZT_QD[z][i][j].n == 0) zoneFZT_QD[z][i][j].n = ndr; // taken from first element matching a larger zone
						zoneFZT_QD[z][i][j].sx += sxdr;
						//zoneFZT_QD[z][i][j].zonenwert = i;
						//zoneFZT_QD[z][i][j].rel += reldr;
						if (zoneFZT_QD[z][i][j].RelCellSum == 0) zoneFZT_QD[z][i][j].RelCellSum = RelCellSumdr; // taken from first element matching a larger zone

						if (zoneFZT_QI[z][i][j].n == 0) zoneFZT_QI[z][i][j].n = nir; // taken from first element matching a larger zone
						zoneFZT_QI[z][i][j].sx += sxir;
						//zoneFZT_QI[z][i][j].zonenwert = i;
						//zoneFZT_QI[z][i][j].rel += relir;
						if (zoneFZT_QI[z][i][j].RelCellSum == 0) zoneFZT_QI[z][i][j].RelCellSum = RelCellSumir; // taken from first element matching a larger zone

						iSourceZoneCntr++;
					}
					// read actual values for output of the reservoir chain (is this really needed...?)
					fscanf(storagefile, "%lf", &(QDd[z][i])); 
					QDd[z][i] /= dMultipleOfOrigTimestep;
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					fscanf(storagefile, "%lf", &(QIi[z][i])); 
					QIi[z][i] /= dMultipleOfOrigTimestep;
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				}
				else
				{ // actual time step is shorter than the source time step, so each source zone must be split into multiple actual zones
					// since all values are read anew from the file and will be spreaded to the existing ones (in case of more actual zones than source zones, which is typical for a shorter time stap than the source file had), we set all original values to 0 now
					for (int j = 1; j <= actMaxzones; j++)
					{
						zoneFZT_QD[z][i][j].sx = 0;
						zoneFZT_QD[z][i][j].n = zoneFZT_QI[z][i][j].n; // zoneFZT_QD[z][i][j].n was set to 0 when creating the zoneFZT_QD array, so we need to set this here to the identical value like in QI-array
						zoneFZT_QI[z][i][j].sx = 0;
					}
					int actZoneCntr = 0;
					for (int jr = 1; jr <= sourceMaxzones; jr++)
					{

						// read the first content from the storage file
						k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &ndr, &sxdr, &zonenwertdr, &reldr, &RelCellSumdr);

						if (k != 5) {
							StorageFileError("number of fields for zonal statistics QD", 5.0, (double)k);
						}
						satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string

						k = fscanf(storagefile, "%lf %lf %lf %lf %lf", &nir, &sxir, &zonenwertir, &relir, &RelCellSumir);

						if (k != 5) {
							StorageFileError("number of fields for zonal statistics QI", 5.0, (double)k);
						}
						satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string

						int corrfact;
						corrfact = (int)dMultipleOfOrigTimestep;
						
						// computing weights for each new zone depending on the cumulated number of cells in each zone of the target record-array
						double n_sum = 0;
						int localCounter = 0;
						for (int j = (jr - 1)*corrfact + 1; j <= jr * corrfact && j <= actMaxzones; j++)
						{
							n_sum += zoneFZT_QI[z][i][j].n;
							localCounter++;
						}

						for (int j = (jr - 1)*corrfact + 1; j <= jr * corrfact && j <= actMaxzones; j++)
						{
							zoneFZT_QD[z][i][j].sx = (n_sum > 0) ? (sxdr) * (zoneFZT_QD[z][i][j].n / n_sum) : 0.0;
							zoneFZT_QI[z][i][j].sx = (n_sum > 0) ? (sxir) * (zoneFZT_QI[z][i][j].n / n_sum) : 0.0;
						}

					}

					// read actual values for output of the reservoir chain (is this really needed...?)
					fscanf(storagefile, "%lf", &(QDd[z][i])); 
					QDd[z][i] /= dMultipleOfOrigTimestep;
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					fscanf(storagefile, "%lf", &(QIi[z][i])); 
					QIi[z][i] /= dMultipleOfOrigTimestep;
					satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				}
			}
		}
	}
	fclose(storagefile);
	satz = loctsatz;
	if (satz) free(satz);
	/* -----------ENDE STORAGEFILE -------------- */
}

void unsatzonmodellklasse::WriteSlrContents(bool appendmode)
{  /* writing contents of overland flow and interflow single linear reservoirs to a file named storage.ftz */
	int i,j,z;
//	if (mpi_rank==0)
//	{

		char *completedname = new char[_stringsize_+1];
		sprintf(completedname,"%s",storagefilename);
		ExtendFileName(DefaultOutputDirectory,completedname);
		if (appendmode==true) // if routing model was active
		storagefile = fopen(completedname,"a");
		else
		storagefile = fopen(completedname,"w");

		delete[] completedname;

		fprintf(storagefile,"[" _Unsatzonmodell_ "_v2]\n"); // writing only new style information to the storage file
		fprintf(storagefile, "%d (number of subasins)\n", zonmax);
		fprintf(storagefile, "%.3lf (timestep in minutes)\n", InV);
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%.0f\t",code[j]); fprintf(storagefile,"(code)\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",kelsqd[j]); fprintf(storagefile,"(kelsqd in min)\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",kelsqi[j]); fprintf(storagefile,"(kelsqi in min)\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",flussdichte[j]); fprintf(storagefile,"(drainage density)\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",kelsqlat[j]); fprintf(storagefile,"(inverse recession constant for baseflow [needed when using no groundwater model])\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",Q0[j]); fprintf(storagefile,"(max. baseflow [needed when using no groundwater model])\n");
		for (int j=1;j<=zonmax;j++) fprintf(storagefile,"%g ",SR_Anteil[j]); fprintf(storagefile,"(surface runoff fraction of snowmelt)\n");
		fprintf(storagefile,"%d (number of tracers)\n",tracer_count);
		for (z=0;z<=tracer_count;z++)
		{
			fprintf(storagefile,"%d (0=runoff, 1...9 = tracers)\n",z);
			for (int i=1;i<=zonmax;i++)
			{
				fprintf(storagefile,"%d %d (zone %d and maxzeit)\n",i,(int)maxzeit[i],i);
				for (int j=1;j<=(int)maxzeit[i];j++)
				{
					fprintf(storagefile,"%lg %lg %lg %lg %lg (zoneFZT_QD.n/.sx/.zonenwert/.rel/.RelCellSum)\n",
					zoneFZT_QD[z][i][j].n,
					zoneFZT_QD[z][i][j].sx,
					zoneFZT_QD[z][i][j].zonenwert,
					zoneFZT_QD[z][i][j].rel,
					zoneFZT_QD[z][i][j].RelCellSum);
					fprintf(storagefile,"%lg %lg %lg %lg %lg (zoneFZT_QI.n/.sx/.zonenwert/.rel/.RelCellSum)\n",
					zoneFZT_QI[z][i][j].n,
					zoneFZT_QI[z][i][j].sx,
					zoneFZT_QI[z][i][j].zonenwert,
					zoneFZT_QI[z][i][j].rel,
					zoneFZT_QI[z][i][j].RelCellSum);
				}
				fprintf(storagefile,"%g (QDd[%d][%d])\n",QDd[z][i],z,i);
				fprintf(storagefile,"%g (QIi[%d][%d])\n",QIi[z][i],z,i);
			}
		}
		fprintf(storagefile,"\n");
		fclose(storagefile);
//	}
	/* -----------ENDE STORAGEFILE -------------- */
}

//#pragma optimize ("p",on)
double **unsatzonmodellklasse::GetQgesAdresse()
{
	if (do_it!=0)
	{
		if (Qg!=NULL)
		{
			if (Qg[0]!=NULL)
			return (double**)Qg;
			else
			{
				(*WSstream) << "FATAL: Qg[0] not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
				return NULL;
			}
		}
		else
		{
			(*WSstream) << "FATAL: Qg not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
			throw(57);
			return NULL;
		}
	}
	else
	return NULL;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double  **unsatzonmodellklasse::GetQDAdresse()
{
	if (do_it!=0)
	{
		if (QD!=NULL)
		{
			if (QD[0]!=NULL)
			return (double**)QD;
			else
			{
				(*WSstream) << "FATAL: QD[0] not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
				return NULL;
			}
		}
		else
		{
			(*WSstream) << "FATAL: QD not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
			throw(57);
			return NULL;
		}
	}
	else
	return NULL;
}

double  *unsatzonmodellklasse::GetKelsQIAdresse()
{
	if (do_it!=0)
	{
		return kelsqi;
	}
	else
	return NULL;
}

double  *unsatzonmodellklasse::GetKelsQDAdresse()
{
	if (do_it!=0)
	{
		return kelsqd;
	}
	else
	return NULL;
}

void unsatzonmodellklasse::SetSurfaceRoutingState(bool bUseSurfaceRouting)
{
	m_useSurfaceRouting = bUseSurfaceRouting;
}


//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double  **unsatzonmodellklasse::GetQIAdresse()
{
	if (do_it!=0)
	{
		if (QI!=NULL)
		{
			if (QI[0]!=NULL)
			return (double**)QI;
			else
			{
				(*WSstream) << "FATAL: QI[0] not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
				return NULL;
			}
		}
		else
		{
			(*WSstream) << "FATAL: QI not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
			throw(57);
			return NULL;
		}
	}
	else
	return NULL;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double  **unsatzonmodellklasse::GetQBAdresse()
{
	if (do_it!=0)
	{
		if (QBb!=NULL)
		{
			if (QBb[0]!=NULL)
			return (double**)QBb;
			else
			{
				(*WSstream) << "FATAL: QBb[0] not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
				return NULL;
			}
		}
		else
		{
			(*WSstream) << "FATAL: QBb not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
			throw(57);
			return NULL;
		}
	}
	else
	return NULL;
}
//#pragma optimize ("p",off)


//#pragma optimize ("p",on)
riverflow_table_type *unsatzonmodellklasse::GetRiverflow_tableAdress()
{
	if (do_it==0)
	return NULL;
	else
	{
		return riverflow_table;
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double unsatzonmodellklasse::quflow(double kf,double h1,double h2,double dz1,double dz2)
{
	return (kf*(h1-h2)/(-0.5*(dz1+dz2))); /* darcy-flow k(Theta)*delta(hydraulic heads)/delta(altitudes) */
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
/* gradient restriction: maximum allowed flux can be restricted in order to avoid a reversal of gradients (prevent from oscillate) 
or to avoid oversaturation or drying */
double unsatzonmodellklasse::gradient_restriction(int j,int s,double thetao,double thetau,
double qup,double qmid,double qlow,double *h,double* &interflow,double* &hh_sl,double DT,double *SE_sl, double *thetar, double *thetas)
{
	double help1,help2, /* help variables */
	h_up_old,h_up_new,h_low_old,h_low_new; /* hydraulic heads for upper and lower layer */
	double theta[3],delta,qmid_orig,sign,redfact,qmid_temp,ifl_orig;
	int k,flag,rev;

	h_up_old=hh_sl[j];                  /* hydraulic head upper layer, old  */
	h_low_old=hh_sl[j+1];                /* hydraulic head lower layer, old */

	theta[1]=thetao+(qmid-interflow[j]-qup)*DT/(h[j-1]-h[j]);  /* what would the water content be in the upper layer after the flow? */
	theta[2]=thetau+(qlow-interflow[j+1]-qmid)*DT/(h[j]-h[j+1]);  /* and what in the lower layer? */

	/* now reducing flux, if oversaturation or drying out would occure */
	double compval=(j>0) ? 
//		ST[s].getThetas(j) :
		thetas[j] :
		LayeredSnowProperties[1].getThetas(levelsSnow+j);

	if ((theta[1]>compval)&&(qmid>0)) /* inflow from lower layer (not too important but yet possible */
	{
		qmid_temp=fabs((compval-thetao)*(h[j-1]-h[j])/DT);
		redfact=qmid_temp/(qmid);
		if (redfact<1.0) {
			qmid *= redfact; /* reducing only with redfact, if maximum allowed flow would be less than the once calculated flow */
		} else {
			qmid=0.0; /* if inflow from below can not be reduced enough in order to avoid oversaturation, it should be set to 0.0 */
		}
		theta[1]=thetao+(qmid-interflow[j]-qup)*DT/(h[j-1]-h[j]);  /* what would the water content be in the upper layer after the flow? */
		theta[2]=thetau+(qlow-interflow[j+1]-qmid)*DT/(h[j]-h[j+1]);  /* and what in the lower layer? */
	}

	/* calculation of new hydraulic heads using the water contents after the flow */
	
	double SE1=(SE_sl!=NULL)?SE_sl[j]:1.0; // for (partly) frozen soils
	double SE2=(SE_sl!=NULL)?SE_sl[j+1]:1.0; // for (partly) frozen soils

	SE1=(j<1)?(1-dSEeffectiveness*(1-SE1)):(1-dSEeffectivenessSoil*(1-SE1));
	SE2=(j<0)?(1-dSEeffectiveness*(1-SE2)):(1-dSEeffectivenessSoil*(1-SE2));

	help1=(j>0)  ? 
			ST[s].h_calc(j,theta[1]*SE1,thetar[j],thetas[j])   : 
			LayeredSnowProperties[1].h_calc(levelsSnow+j,theta[1]*SE1);
	help2=(j>=0) ? 
			ST[s].h_calc(j+1,theta[2]*SE2,thetar[j],thetas[j]) : 
			LayeredSnowProperties[1].h_calc(levelsSnow+j+1,theta[2]*SE2);
	
	h_up_new=(h[j]+h[j-1])*0.5-help1;       /* hydraulic head upper layer, new  */
	h_low_new=(h[j+1]+h[j])*0.5-help2;       /* hydraulic head lower layer, new */

	if ((h_up_old>h_low_old && (h_up_new+0.001)<h_low_new) || (h_up_old<h_low_old && (h_up_new-0.001)>h_low_new))
	{ /* a gradient reversal must be avoided */
		k=0;
		delta=0.5;
		qmid_orig=qmid;
		ifl_orig=interflow[j];
		sign=-1;
		rev=1;
		do /* because of numeric errors, because of the nonlinearity of h-theta-relation and because of edges 
within the h-theta-look-up-table it is necessary to solve the
gradient restrictions in an iterative way  */
		{
			k++;                                       
			flag=0;
			qmid+=delta*sign*qmid_orig; /* first step: 1/2 of qmid is subtracted */
			interflow[j]+=delta*sign*ifl_orig;
			delta*=0.5; /* next step: 1/4 of qmid is added or subtracted (depends on the relation of h_up_new and h_low_new) */

			theta[1]=thetao+(qmid-interflow[j]-qup)*DT/(h[j-1]-h[j]);  /* what would the water content be in the upper layer after the flow? */
			theta[2]=thetau+(qlow-interflow[j+1]-qmid)*DT/(h[j]-h[j+1]);  /* and what in the lower layer? */
			
			/* calculation of new hydraulic heads using the water contents after the flow */
			help1=(j>0)  ? 
				ST[s].h_calc(j,theta[1]*SE1, thetar[j], thetas[j])   :
				LayeredSnowProperties[1].h_calc(levelsSnow+j,theta[1]*SE1);
			help2=(j>=0) ? 
				ST[s].h_calc(j+1,theta[2]*SE2, thetar[j], thetas[j]) :
				LayeredSnowProperties[1].h_calc(levelsSnow+j+1,theta[2]*SE2);

			h_up_new=(h[j]+h[j-1])*0.5-help1;       /* hydraulic head upper layer, new  */
			h_low_new=(h[j+1]+h[j])*0.5-help2;       /* hydraulic head lower layer, new */
			
			if ((h_up_old>h_low_old && (h_up_new+0.001)<h_low_new) || (h_up_old<h_low_old && (h_up_new-0.001)>h_low_new))  
			{
				flag=1; /* if still a gradient reversal occurs, while loop has to be continued, flux is still too large -> no sign change */
				if (!rev) 
				{ 
					rev=1; 
					sign*=(double)-1.0;
				}
			}
			else /* no gradient reversal: are new heads identical (+- 1 mm) or are they to different? */
			if (fabs(h_low_new-h_up_new)>0.001) 
			{
				flag=1; 
				if (rev)
				{
					rev=0;
					sign*=(double)-1.0;
				}
			} /* flux change was too large, 1/2 of the last change is now compensated */
		}
		while (flag && k<51);
	}
	return qmid;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double unsatzonmodellklasse::GetMeanTheta(int &s,double* &dicken,double* &theta_sl,double upper,double lower,int AsAbsoluteAmount)
//double unsatzonmodellklasse::GetMeanTheta(double upper,double lower,int AsAbsoluteAmount,int &s,double* &dicken,double* &theta_sl)
/* returns the mean relative water content for a given soil depthrange
upper: depth of the upper limit of soil depth range (e.g. 0 (surface) for 0..x m , 1 for 1..x m)
lower: depth of the lower limit ôf the soil depth range (e.g. 1 for range x..1 m, 2 for x..2 m)
*/
{
	int k;
	double temp,total_thickness,dicksum,mean_rel_theta;
	/* chscking for realistic values: upper must be greater than lower, both have to be greater than zero */
	if ((lower<=(double)0)||(upper<(double)-0.0001)) return nodata;
	double maxdepthforsoil=0;
	for (int zcnt=1;zcnt<=levels;zcnt++)
	{
		maxdepthforsoil+=dicken[zcnt];
	}
	if ((lower>maxdepthforsoil-1e-8)||(lower==9999.0)) lower=maxdepthforsoil-1e-8;
	total_thickness = lower-upper;

	if (total_thickness<=1e-8) return 1.0; /* soil is completely saturated */
	mean_rel_theta = 0.0;
	dicksum = upper;
	if (upper<1e-8) k=1; /* in which soil layers starts the calculation? */
	else
	{ /* searching the layer in which upper is within */
		temp = 0.0;
		k=1;
		while ((temp+dicken[k])<=upper)
		{
			temp+=dicken[k];
			k++;
		}
		/* calculating the relativ theta for the first layer
(extra calculation, because layer thickness is not necessarily the full thickness of the layer) */
		if ((temp+1e-8)<upper)
		{
			if ((temp+dicken[k])>lower) return (double)(theta_sl[k]/ST[s].getThetas(k)); /* if this is the only layer (just a portion of the layer is of interest) */
			dicksum = temp + dicken[k];
			/* if there are more layers to consider, the actual layer is weighted by the relative thickness
	of the part between "upper" and the layer lower limit */
			mean_rel_theta = (AsAbsoluteAmount==1)?
			((theta_sl[k])*(dicksum-upper)): // absolute
			((theta_sl[k]/ST[s].getThetas(k))*(dicksum-upper)/total_thickness); // relative
			k++; /* now the layer index is increased for using the next layer */
		}
	}
	/* adding of complete layers */
	while ((dicksum+dicken[k])<lower)
	{
		/* weighting the relative thetas with its fractions of upper ... lower */
		mean_rel_theta += (AsAbsoluteAmount==1)?
		(theta_sl[k]*dicken[k]):
		((theta_sl[k]/ST[s].getThetas(k))*(dicken[k]/total_thickness));
		dicksum+=dicken[k];
		k++;
	}
	/* if the last layer was to thick for "lower", again only a fraction of this layer is considerd */
	if ((dicksum+1e-8)<lower) /* 1e-10 is for avoiding numerical problems */
	{
		mean_rel_theta += (AsAbsoluteAmount==1)?
		(theta_sl[k]*(lower-dicksum)):
		((theta_sl[k]/ST[s].getThetas(k))*(1.0-(dicksum-upper)/total_thickness));
		/* 1.0 - (dicksum-upper)/total_thickness is the relative fraction of the last layer on the soil depth range given by (upper..lower) */
	}
	return (double)mean_rel_theta;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
int unsatzonmodellklasse::evapotranspiration(int i,int j,double a,int s,double dReduParDry,double dMaxThetaWet,double dLimitReduWet,
double *a_dicke,double *a_etr,double *a_evar,
double &wurz,double &vbg,double &etp,double evap,
double* &dicken,double* &thetasat,double* &thetares,double* &theta_sl,double* &SE_sl,
double* &hgeo_sl,double* &etrdetail,double* &extractionweight,int gwlevel,double DT,double evapmaxdepth)
//void unsatzonmodellklasse::evapotranspiration(double a,double theta_grenz,double dMaxThetaWet,double dLimitReduWet,double *a_dicke,double *a_etr)
/* extracts for the cell i,j from all layers within root depth a fraction of avaporation water
--> extraction profile of transpiration over root depth is calculated considering
root density and soil moisture
parameters:
a           parameter for root density with depth: -1 concav, around 1: linear decreasing with depth, >1 convex, >100...1000 very convex
theta_grenz lower limit of relative soils moisture allowing the transpiration to be equal the potential evapotranspiration
a_dicke     array, containig relative fractions of the (considered) layer thicknesses at the root depth -> return values
a_etr       array with return values of the extraction per layer

all the other parameters are used here for protecting the function in a multi threading environment of changes by code running in parallel 
*/
{
	/* etp ist already converted into m/s */
	int k,lrs,maxlayer;
	double temp,temp2,maxdepth,depth,sum_extractiondensity,redu,diff,fact,theta_grenz,tet;
	diff = (double)0.00001;  /* for avoiding numerical problems */
	//   theta_grenz *= theta_sl[0]; /* relative auf absolute Bodenfeuchte umrechnen
	//                                  (anstatt bezogen auf Porenraum nun bezogen auf absolutes Volumen) */
	if (wurz<dicken[1]) /* wurz: root depth */
	{
		/* if root depth < thickness of the first layer, extraction is done only from this layer*/
		wurz = dicken[1];
		//vbg=0.0;
	}
	maxdepth = hgeo_sl[0]-wurz-diff; /* geodaetic altitude of the lower root zone limit;
								does not have to coincide with a layer boundary! */

	/* 1.: extraction for the vegetation covered fraction of the call (vbg: vegetation coverage degree) */
	sum_extractiondensity=redu=0; /* initialization of some variables */
	k=0;
	do /* from the first layer a extraction should be possible in any way */
	{
		k++;
		a_dicke[k] = hgeo_sl[k] > maxdepth ? dicken[k] : hgeo_sl[k-1]-maxdepth-diff ; /* considered layer thickness within root depth */
		depth = hgeo_sl[k]>maxdepth ? (hgeo_sl[0]-hgeo_sl[k])-(a_dicke[k]*0.5) : wurz-a_dicke[k]*0.5; /* average depth of the actual layer */
		a_dicke[k]/=wurz; /* absolute -> relative fraction of layer thickness at root depth */

		extractionweight[k] = (1.0- depth /wurz)*(1.0+a* depth /wurz)*a_dicke[k]; /* calculating extraction density */
		sum_extractiondensity += extractionweight[k]; /* sum this density to get later a weighting with its sum to be 1.0 */
	}
	while ((hgeo_sl[k]>maxdepth)&&(k<levels));
	maxlayer=k;

	/* now weighting the extraction densities to get a sum of 1.0 */
	for (lrs=1;lrs<=maxlayer;lrs++) extractionweight[lrs]/=sum_extractiondensity;

	/* considering reductions due to too high or too less soil moistures */
	double eta=0;	
	for (int k=1;k<=maxlayer;k++)
	{
		theta_grenz = ST[s].get_theta(k,dReduParDry,thetares[k],thetasat[k]);
		tet = (SE_sl)?(SE_sl[k] < 1.0 ? SE_sl[k] * theta_sl[k] + thetares[k] : theta_sl[k]):theta_sl[k]; // reducing theta for frozen layers
		if (theta_grenz <=ST[s].getThetar(k))
		{
			theta_grenz=(2*ST[s].getThetar(k)+ST[s].getThetas(k))/3.0; /* standard-value: weighted average to residual water content */
		}
		
		fact=2.0; /* factor for extraction density in the next lower layer, which gives the maximum amount of reduced extraction density to be added to the next layer */
		temp=1.0; /* standard is no reduction! -> reduction factor is 1.0 */
		if (tet>dMaxThetaWet*thetasat[k]) /* because of anaerobe consitions reduction of transpiration between soil moistures of 0.9 ... 095 */
		{
			/* linear reduction of etr from etp to dLimitReduWet */
			temp = interpolate(dMaxThetaWet*thetasat[k],1.0,thetasat[k],dLimitReduWet,(tet>thetasat[k])?thetasat[k]:tet);
		}
		else if (tet<theta_grenz) /* for the lower end of soil moisture */
		{
			/* below soil moistures of "theta_grenz" ()relative soil moisture of e.g. 0.6) a linear reduction of etr is done */
			temp=tet>thetares[k]?(tet-thetares[k])/(theta_grenz-thetares[k]):0;
		}

		if (temp<1.0)
		{
			temp2=(1-temp)*extractionweight[k]; /* reduced amount in terms of extraction density */
			if ((k<maxlayer)&&(k<gwlevel))
			{ /* reduced amount of the actual extraction density is added to the next lower layer in order to get the water from there then */
				if (temp2>fact*extractionweight[k+1])
				{ /* effective reduction is now calculated in terms of extraction density */
					redu+=temp2-fact*extractionweight[k+1]; /* what cannot be added to the next layer is still possible to get from the surface */
					extractionweight[k+1]*=(1+fact);
				}
				else
				extractionweight[k+1]+=temp2;
			}
			else
			{
				redu+=temp2; /* reduced amount can still get from the first layer by evaporation */
			}
			extractionweight[k]-=temp2;
		}
	}

	/* extracting the correct fraction of etp from each layer within root depth */
	double dicksum=0;
	for (int k=1;k<=maxlayer;k++)
	{  /* only the vegetated fraction is considerd here!
	!!! layer thickness is also considered  !!!
	!!! weighting to 1.0 is done later !!!,
	because the sum of the extraction densities  (extractionweight) is calculated later */
		dicksum+=dicken[k];
		a_etr[k] = vbg*etp*extractionweight[k];
		//no need for further reduction
		//if (THAWDEPTH)
		//{	
		//	if (THAWDEPTH[i][j]>nodata)
		//	{
		//		if (dicksum>THAWDEPTH[i][j])
		//		{
		//			a_etr[k]; // no extraction by vegetation from frozen soil
		//		}
		//	}
		//}
		etrdetail[k]+=a_etr[k];
	}

	redu *= etp; /* redu was given in terms of extraction density, now it is converted to m/s considering potential evapotranspiration */

	/* 2.: now extraction of soil evaporation from the uppermost layer(s)
		if this level(s) is/are empty, no soil evaporation takes place, only transpiration from deeper layers by root extraction! */
	dicksum=0;
	k=1;
	double actMaxAvailable = evap+vbg*redu; // initial vlue: this is the total available amount for soil evaporation.
	while (dicksum<evapmaxdepth || k<2) // one loop must be guaranteed --> k is 1 in the first loop
	{
		if (!SE_sl) // classical quadratic reduction funktion for applications without heat transfer - to keep the model compatible with old setups
		{
			temp2 = (thetasat[k] - theta_sl[k]) / thetasat[k]; /* calculating a quadratic reduction function depending on soil mooisture */
			temp2 *= temp2;  /* quadratic weighting: only low moistures should reduce evaporation strongly */
			temp2 = 1 - temp2; /* to ensure that at theta=1 the weight is 1.0 and at theta=0 the weigth is 0.0 */
			//   temp = ((1-vbg)*etp+vbg*redu)*temp2; /* max. possib. extraction from layer one in %vol by evaporation */
			temp = actMaxAvailable * temp2; /* max. possib. extraction from layer one in %vol by evaporation */
			/* etr from bare soil is stored in a_evar */

			a_evar[k] = (temp + vbg * a_etr[k])*DT / dicken[k] < theta_sl[k] ? temp /* if there is enough moisture in layer 1, all potential extraction is extracted */
				: (1.0 - vbg)*theta_sl[k] * dicken[k] / DT; /* else: the soil is dryed to theta=0 */
		}
		else
		{

			temp2 = (theta_sl[k] - thetares[k]) / (thetasat[k] - thetares[k]); /* calculating a linear reduction function depending on soil mooisture */
			//temp2 *= temp2;  /* quadratic weighting: only low moistures should reduce evaporation strongly */
			if (SE_sl[k] < 1.0) temp2 = SE_sl[k] * temp2; // reduction of evap if frozen
			//   temp = ((1-vbg)*etp+vbg*redu)*temp2; /* max. possib. extraction from layer one in %vol by evaporation */
			temp = actMaxAvailable * temp2*(fmax(evapmaxdepth - dicksum, 0) / evapmaxdepth); /* max. possib. extraction from layer one in %vol by evaporation */
			/* etr from bare soil is stored in a_evar */

			a_evar[k] = theta_sl[k] - (temp + vbg * a_etr[k])*DT / dicken[k] > thetares[k] ? temp /* if there is enough moisture in layer 1, all potential extraction is extracted */
				: (1.0 - vbg)*(theta_sl[k] - thetares[k])*dicken[k] / DT; /* else: the soil is dryed to theta=thetares */
		}

		if (a_evar[k]>actMaxAvailable)
		{
			a_evar[k]=actMaxAvailable;
		}
		// reducing the available evaporation for the next layer according to the difference of the total avauilable for this layer and the amount actaually used for the actual layer
		actMaxAvailable -=a_evar[k]; 
		etrdetail[0]+=a_evar[k];

		dicksum+=dicken[k];
		k++;
	}
	//   a_etr[1]+=a_etr[0]; /* the extraction for the layer 1 is the sum of evaporation and transpiration */
	return maxlayer;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::hydraulic_heads_iterative(int offset,int i,int j,int z,
int &maxlevel,int &s,double &h_temp,
double* &gwin,double* &gwex,double* &h_temp1,
double* &theta_sl_old,double* &theta_average,
double* &hh_sl,double* &hgeo_sl,double * SE_sl, double *thetar, double *thetas)
/* this version is for implementation within the FDM-scheme. it considers the old and the new (actual iteration value) theta */
{
	int k,l;
	double h2;
	/* maxlevel: the level to which evapotranspiration took place -> this may be below groundwater table! */
	l=maxlevel > (int)GWLEVEL[i][j] ? maxlevel : (int)GWLEVEL[i][j];
	/* also groundwater exfiltration can occure only below or at groundwater table, therefore these layers have to be considered in the scheme */
	if (l<levels)
	{
		while (gwin[l+1]>0||gwex[l+1]>0) 
		{
			l++;
			if (l>=levels) break;
		}
	}
//	if (l>levels) l=levels; 
	l+=2; /* to avoid fluxes within the saturated zone, at least the first both saturated layers should also be included 
			(one layer for considering a change in groundwater level and one for considering the flux between this changed layer and the next layer */
	if (l>levels) l=levels; 
	for (int k=1+offset;k<=l && k<levels;k++)
	{
		double SE=(SE_sl!=NULL)?SE_sl[k]:1.0; // for (partly) frozen soils
		SE=(k<1)?(1-dSEeffectiveness*(1-SE)):(1-dSEeffectivenessSoil*(1-SE));
		/* if z!=0, the same h_temp1 value is used, so a new calculation would only waste time */
		if (z==0) h_temp1[k] = (k>0) ? 
				ST[s].h_calc(k,theta_sl_old[k]*SE,thetar[k],thetas[k]) : 
				LayeredSnowProperties[1].h_calc(levelsSnow+k,theta_sl_old[k]*SE); /* suction for the theta from the last timestep */
		
		h2 = (k>0) ? 
				ST[s].h_calc(k,theta_average[k]*SE, thetar[k], thetas[k]) :
				LayeredSnowProperties[1].h_calc(levelsSnow+k,theta_average[k]*SE); /* suction for the actual theta */

		h_temp=(h_temp1[k]+h2)/2.0;
		hh_sl[k]=(hgeo_sl[k-1]+hgeo_sl[k])*0.5; /* calculating hydraulic head */
		if (k>=(int)GWLEVEL[i][j])
		{
			if ((k==(int)GWLEVEL[i][j])&&(hgeo_sl[0]+GWSTAND[i][j]<=hgeo_sl[k-1])&&(hgeo_sl[0]+GWSTAND[i][j]>=hgeo_sl[k]))
			{
				hh_sl[k]=hgeo_sl[0]+GWSTAND[i][j];
			}
			else if (hh_sl[k]<hh_sl[(int)GWLEVEL[i][j]])
			{
				hh_sl[k]=hh_sl[(int)GWLEVEL[i][j]];
			}
			else 
			hh_sl[k]-=h_temp;
		}
		else hh_sl[k]-=h_temp;
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
/* calculation of infiltration into soil from rivers */
void   unsatzonmodellklasse::calc_surface_interaction(int i,int j,int lr,
double* &gwin,double* &gwex,double* &conc_gwin,
double* &thetares, double* &thetasat, double* &theta_sl_old,double* &theta_average,
int &s,double* &hgeo_sl,double* &hh_sl,double* &dicken,double &GWINtotal,double &GWEXtotal,
double DT)
{ /* lr is the Groundwater level */
	int r, /* r: index of the layer, which is at the same level as the river water table */
//	k, /* loop index */
	m, /* layer, which is in contact with the river bed */
	z; /* index for tracer species */
	int lnknr; /* alias for the value at cell i,j in the RIVLNK-grid */
	double A,       /* flow profile area */
	T,       /* flow depth in the river */ 
	Hriver,  /* hydraulic head in the river in m.a.s.l. */ 
	Hground, /* hydraulic head in the communicating layer in the soil (in m.a.s.l.) */ 
	dH,      /* potential difference */ 
	maxgwex, /* what amount of water can flow into the river until the gradient is zero */ 
	drain_theta, /* theta that is the result of a gravitative drainage of the actual layer */
	temp;       /* infiltration amount in m3/s */
	double RivDep=0,
	RivWit=0,
	Kolm=0;
	double act_theta;
	int zcnt;
	if (surface_interaction<1) 
	return;
	else
	{
		for (int k=1;k<=levels;k++) gwex[k]=gwin[k]=0.0; /* initializing */

		if ((lake_active==0)||
				((lake_active==1)&&(LAKE[i][j]==nodata))||
				((lake_active==1)&&(LAKE[i][j]!=nodata)&&(POND[i][j]<1.0)))
		{
			RivDep=RIVDEP[i][j];
			RivWit=RIVWIT[i][j];
			Kolm=KOLM[i][j];
		}
		else
		{
			RivDep=0;
			RivWit=0;
			Kolm=0;
		}
		
		if ((RivDep>0)&&(RivWit>0)) /* only for valid river network cells */
		{
			// this was the old method...:
			//  	m = RivDep>dicken[1]?2+(int)(RivDep/dicken[0]):1; /* index of the layer, which is in contact with the river bottom */
			double firstlayerincontact=dicken[1];
			m=1;
			for (zcnt=1;zcnt<=levels;zcnt++)
			{
				if (RivDep>firstlayerincontact)
				{
					firstlayerincontact+=dicken[zcnt];
					m++;
				}
				else
				{
					break;
				}
			}
			if (m>levels)
			{
				m = levels;
			}

			Hriver     = hgeo_sl[0]-RivDep; /* water level in the river is neglected for exfiltration */
			Hground    = hgeo_sl[0]+GWSTAND[i][j];
			/* 2 possibilities: 1. exfiltration (hground>Hriver), 2: infiltration (Hriver>Hground) */
			if (Hground>Hriver)
			{ /* exfiltration -> generating baseflow */
				dH=Hground-Hriver;
				GWINtotal=0;
				/* potential amount of exfiltrating water as given by the hydraulic gradient and the width of the river bed and its resitance to flows */
				GWEXtotal = Kolm>0? Kolm*dH*RivWit/cellsize : 0;
				int k=lr;
				act_theta=(theta_sl_old[k]+theta_average[k])/2.0; /* theta for the actual layer k as valid in this iteration step */
				/* assuming equal properties in all affacted layers, the drainable porosity act_theta-get_theta() times dH/DT is the rate of drainage, 
	which would empty the pore volume between Hground and Hriver in one single time step */
				drain_theta=ST[s].get_theta(k,3.45, thetares[k], thetasat[k]);
				temp=GWEXtotal;
				if (k<m) 
				{
					/* first layer */
					maxgwex=(act_theta-drain_theta)*(hgeo_sl[0]+GWSTAND[i][j]-hgeo_sl[k])/DT; /* what amount of water can drain out of the layer */
					if (maxgwex<0.0) maxgwex=0.0;
					if (maxgwex>temp) maxgwex=temp; /* if this is more than the allowed rate, then limit it */
					gwex[k]=maxgwex; /* gwex is set to this (ev. limited) rate  */
					temp-=gwex[k];   /* what amount is still available for other layers? */
				}
				else /* k==m */
				{
					maxgwex=(act_theta-drain_theta)*dH/DT;
					if (maxgwex<0.0) maxgwex=0.0;
					if (maxgwex>temp) maxgwex=temp; /* if this is more than the allowed rate, then limit it */
					gwex[m]=maxgwex; /* gwex is set to this (ev. limited) rate  */
					temp-=gwex[m];   /* what amount is still available for other layers? */
				}
				for (int k=lr+1;k<m;k++) 
				{
					maxgwex=(act_theta-drain_theta)*dicken[k]/DT;
					if (maxgwex<0.0) maxgwex=0.0;
					if (maxgwex>temp) maxgwex=temp; /* if this is more than the allowed rate, then limit it */
					gwex[k]=maxgwex; /* gwex is set to this (ev. limited) rate  */
					temp-=gwex[k];   /* what amount is still available for other layers? */
				}
				if (k<m)
				{
					/* last affected layer */
					maxgwex=(act_theta-drain_theta)*(hgeo_sl[m-1]-Hriver)/DT;
					if (maxgwex<0.0) maxgwex=0.0;
					if (maxgwex>temp) maxgwex=temp; /* if this is more than the allowed rate, then limit it */
					gwex[m]=maxgwex; /* gwex is set to this (ev. limited) rate  */
					temp-=gwex[m];   /* what amount is still available for other layers? */
				}
				if (temp>0) GWEXtotal-=temp;
			}
			else
			{ /* infiltration */
				GWEXtotal=0;
				if ((RIVLNK[i][j]>0)&&(surface_interaction>1)) /* only in such channels, which get water from upstream areas via river links */
				{
					lnknr=zarray[(int)RIVLNK[i][j]];
					A  = riverflow_table[lnknr].speed>0?riverflow_table[lnknr].flowm3s/riverflow_table[lnknr].speed:0.0; /* discharge profile area */
					T  = A/RivWit;
					Hriver     = hgeo_sl[0]-RivDep+T; /* with consideration of water level in the river */
					/* if no direct contact with groundwater: 
	hydraulic head of the layer which is in contact with the river bed is used, otherwise the head of the groundwater table is used */
					Hground    = m<=lr?hh_sl[m]:hgeo_sl[0]+GWSTAND[i][j];
					dH = Hriver-Hground;
					/* potential amount of infiltrating water as given by the hydraulic gradient and the width of the river bed and its resitance to flows */
					GWINtotal = Kolm>0? Kolm*dH*RivWit*cellsize : 0; /* in m3/s for this grid cell */
					if (GWINtotal<0) /* may be, due to the fact, that hydraulic head for the layers is valid for  */
					GWINtotal=0;
					// since other threads may asking for or even updating the riverflow_table registers, we must make this bloc critical
#pragma omp critical
					{
						if (GWINtotal>((riverflow_table[lnknr].flowm3s/substepcnt) - riverflow_table[lnknr].infilt_amount)) /* balancing with possible infiltration as given by the discharge in the channel */
						GWINtotal= (riverflow_table[lnknr].flowm3s / substepcnt) -riverflow_table[lnknr].infilt_amount;
					}
					GWINtotal/=(cellsize*cellsize); /* converting back from m3/s into m/s */
					if (GWINtotal>0)		
					{
						/* filling up all affected layers */
						//            r=RivDep+T>dicken[1]?2+(RivDep+T)/dicken[0]:1; /* index of the uppermost affected layer */
						double firstlayerincontact=0;
						r=1;
						for (zcnt=1;zcnt<=levels;zcnt++)
						{
							if (RivDep-T>firstlayerincontact)
							{
								firstlayerincontact+=dicken[zcnt];
								r++;
							}
							else
							{
								break;
							}
						}

						temp=GWINtotal;
						for (int k=m;k<=lr;k++)
						{ /* filling downwards */
							act_theta=(theta_sl_old[k]+theta_average[k])/2.0; /* theta for the actual layer k as valid in this iteration step */
							gwin[k]=(thetasat[k]-act_theta)*dicken[k]/DT;
							if (gwin[k]>temp) gwin[k]=temp;
							for (z=1;z<=tracer_count;z++)
							conc_gwin[z]=riverflow_table[lnknr].flow_mm[z];
							temp-=gwin[k];
						}
						if (temp>0)
						for (int k=m-1;k>r;k--)
						{ /* and now filling upwards */
							act_theta=(theta_sl_old[k]+theta_average[k])/2.0; /* theta for the actual layer k as valid in this iteration step */
							gwin[k]=(thetasat[k]-act_theta)*dicken[k]/DT;
							if (gwin[k]>temp) gwin[k]=temp;
							for (z=1;z<=tracer_count;z++)
							conc_gwin[z]=riverflow_table[lnknr].flow_mm[z];
							temp-=gwin[k];
						}
						if (temp>0) GWINtotal-=temp; /* if there was not enough free pore volume, GWIN is reduced */
					}
				}
			}
		}
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::specific_fluxes(int offset,int i,int j,int s,int maxlevel,double draindist,int clay_layer,double p,
double* &interflow,double* &qu_sl,double* &gwin,double* &gwex,double* &conc_gwin,
double* &thetares,double* &thetasat,double* &theta_sl_old,double* &theta_average,
double* &hgeo_sl,double* &hh_sl,double* &dicken,double &qdrain,double &GWINtotal,double &GWEXtotal,
double DT,int do_TMod,double *SE_sl)
{
	int k,l,m,start,filled,TG;
	double ifl_min_theta,max_ifl,qtmp,act_theta,act_theta_low; /* help variables for specific fluxes */
	double *ksat;     // Substitution variable for double *SOILtable[SOIL[i][j]].k_sat and LayeredSnowProperties[1].k_sat 
	double keff,ko,ku; /* keff=effektive,ko und ku = relative hydr. Leitfaehigkeiten in der betrachteten und der darunter (u) liegenden Schicht */

	// here, we will manipulate the ksat-array when a THAWDEPTH-grid is available
	// goal: soil layers below thaw depth should not allow water movement at all (or only at a minimum rate)
	double *m_ksat = new double[levels+1+levelsSnow];
	ksat    = &(m_ksat  [levelsSnow]); 

	memmove(&(m_ksat [levelsSnow]), ST[s].getK_sat(),sizeof(double)*(levels+1)); 
	if (do_LayeredSnow) 
	{
//		memmove(&(m_ksat[levelsSnow+offset+1]), this->LayeredSnowProperties[1].getK_sat(),sizeof(double)*(-offset));
		memmove(&(m_ksat[0]), this->LayeredSnowProperties[1].getK_sat(),sizeof(double)*(levelsSnow+1));
	}


	// old perma frost model, not required for new model since k_sat depends on SE_sl then which is a parameter to this method
	if (THAWDEPTH)
	{
		if (THAWDEPTH[i][j]>nodata)
		{
			int z=1;
			double dicksum=0;
			for (int z=1;z<=levels;z++)
			{
				dicksum+=dicken[z];
				if (dicksum>THAWDEPTH[i][j])
				{
					// for layers below thaw depth, set minimum hydraulic conductivity
					ksat[z]=reduced_ksat;
					if ((z==1)&&(ksat[z]<fabs(p)))
					{ // reduce "precipitation" p to the amount defined by ksat*timestep
						p=-ksat[z];
					}
				}
			}
		}
	}

	TG      = zarray[(int)zn[i][j]]; /* which subbasin */

	qu_sl[0+offset] = -p; /* inflow into 1st layer is outflow from "surface"
				-> precipitation + irrigation - infitration excess (in case of ponding: infitration excess is all water that remains inthe pond) */

	/* Interaction with surface water */
	if (gwflow_active)
	{
		// if ((i==10)&&(j==30)&&(stunde==12))
		// {
		// 	(*WSstream) << "\nStunde " << stunde << " GWEXtotal vor calc_surface_interaction" << GWEXtotal*DT*1000 << " timestep " << DT*1000 << " DT " << DT << "[end-entry]";
		// 	for (int xx=1;xx<=9;xx++)
		// 	{
		// 		(*WSstream) << "gwex["<<xx<<"]: " << gwex[xx] << " gwin["<<xx<<"]: " <<gwin[xx] << "[end-entry]";
		// 	}
		// }
		calc_surface_interaction(i,j,(int)GWLEVEL[i][j],	gwin,gwex,conc_gwin,thetares,thetasat,theta_sl_old,theta_average,s,hgeo_sl,hh_sl,dicken,GWINtotal,GWEXtotal,DT);
		// if ((i==10)&&(j==30)&&(stunde==12))
		// {
		// 	(*WSstream) << "\nStunde " << stunde << " GWEXtotal nach calc_surface_interaction" << GWEXtotal*DT*1000 << " timestep " << DT*1000 << " DT " << DT <<  "[end-entry]";
		// 	for (int xx=1;xx<=9;xx++)
		// 	{
		// 		(*WSstream) << "gwex["<<xx<<"]: " << gwex[xx] << " gwin["<<xx<<"]: " <<gwin[xx] << "[end-entry]";
		// 	}
		// }
	}

	int firstExcessIceLevel = -9999;
	if (do_ArcVeg > 0)
	{
		if (ExcessIceFactor)
		{
			k = 1;
			while (k<levels)
			{
				if (ExcessIceFactor[i][j][k] > 0)
				{
					firstExcessIceLevel = k;
					break;
				}
				k++;
			}
		}
	}

	/* maxlevel is the deepest level from which extraction by vegetation took place */
	l=maxlevel > (int)GWLEVEL[i][j]+1 ? maxlevel : (int)GWLEVEL[i][j]+1;
	if (l<levels) //20230801 l must reach to the first excessIce layer at least (if excessIce is present), otherwise no excessIce melt would be calculated in groundwater saturated layers
	{
		while ((gwin[l+1]>0) || (gwex[l+1]>0) || (l< firstExcessIceLevel))
		{
			l++;
			if (l>=levels) break;
		}
	}
	if (l>=levels) l=levels-1; /* no interflow below the forelast layer -> prevent numerical problems */

	//  m=(int)(gw_depth/dicken[0])+2;
	/* the following sequence is used only by the infiltration (GREEN-AMPT) Model, which sets gw_depth eventually to a new value
	since infiltration is not used any more, we comment this
	m=2;
	double gwcontact=0;
	for (int zcnt=1;zcnt<=levels;zcnt++)
	{
	if (gwcontact<gw_depth)
	{
		gwcontact+=dicken[zcnt];
		m++;
	}
	else
	{
		break;
	}
	}
	*/
	/* if the infiltration amount exceeds the water holding capacity of the first soil level,
a fraction of infitration has to infiltrate into the second layer and so on*/

//	k=start=1+offset;
//	filled=0;

//	k=m=l; /* the comment should be removed, if infiltration from infiltration submodel should NOT forced into soil (infiltration is then only the amount of flow following richards equation) */
//	qtmp=0.0;
//	while ((k<l)&&(k<m)&&(qu_sl[k-1]<0)&&((theta_sl_old[k]-(qu_sl[0]+interflow[k])*DT/dicken[k])>ST[s].getThetas(k)))
//	{
//		act_theta = ST[s].getThetas(k);
//		double SE=(SE_sl!=NULL)?SE_sl[k]:1.0; // for k-reduction in frozen soils
//		qu_sl[k]=qu_sl[0]+(act_theta-theta_sl_old[k])*dicken[k]/DT; /* difference between the infiltrating amount and the capacity
//																	of the actual layer is the input into the underlaying layer */
//		if (hgeo_sl[k]-hh_sl[k]<3.45)
//		/* calculation of interflow, if saturation is high enough (suction is below 3.45 m) */
//		{
//			ifl_min_theta =ST[s].get_theta(k,3.45); /* getting the theta, interflow may decrease actual theta maximal */
//			max_ifl = (act_theta-ifl_min_theta)>0?(act_theta-ifl_min_theta)*dicken[k]/DT:0; /* maximal interflow is the draiageable theta */
//			if (ksat[k]>ksat[k+1])
//			//      if (ST[s].getK_sat()[k]>ST[s].getK_sat()[k+1])
//			{
//				interflow[k]= ST[s].k_calc(k,act_theta,SE)*ksat[k]*dicken[k]*flussdichte[TG];
//				interflow[k]*= (TANBETA[i][j]>1?1:TANBETA[i][j]);
//			}
//			else interflow[k]=0.0;
//			/* if k is the drainage level, then a drainage flow may be added to the interflow */
//			if (HDRAIN)
//			if (k==DRAINLEVEL[i][j])
//			{ 
//				qdrain=ST[s].k_calc(k,act_theta,SE)*2*ksat[k]*dicken[k]/draindist*cellsize/draindist;
//				// flow clipping for drainage: no more drainage is allowed than is drainable from the actual layer
//				// and this is defined by the suctions
//				if (qdrain>max_ifl)
//				{
//					qdrain=max_ifl;
//				}
//				/* 2*dicken[k]/draindist is a kind of hydraulic gradient, cellsize/draindist is the number of drainage pipes within a cell */
//				interflow[k]+=qdrain;
//				/* by adding the drainage flow to the interflow, the FDM-scheme should still work properly */
//			}
//			if (max_ifl<interflow[k]) interflow[k]=max_ifl;
//		}
//		else interflow[k]=0.0;
//
//		/* substracting the interflow from the qu-flow */
//		if ((qu_sl[k]+interflow[k])<0)  /* but only if interflow is less than qu */
//		qu_sl[k]+=interflow[k];
//		else  /* otherwise interflow[k] should be zero */
//		interflow[k]=0;
//
//		qtmp=qu_sl[k]; /* value of qu_sl[k] is stored into a temporary variable in order to compare 
//				the later calculated value of the flow into the first not-satureted layer 
//				with this value (if calculated flow is larger than qtmp, qu_sl[k] is replaced by the larger value) */
//		qu_sl[0]-=(qu_sl[k]-interflow[k]); /* for the actual layer the remaining amount is the differenz between input and output
//						(matches the actual empty porosity) */
//		k++;
//		start=k-1;
//	}
	
	qu_sl[0+offset] = -p;

	for (int k=1+offset;k<=l && k< levels;k++)
	{
		/* estimation of k-value */
		/* in the following rows for each layer at first the effective conductivity, then the flux is calculated */

		act_theta=(theta_sl_old[k]+theta_average[k])/2.0; /* averaging thetas from last timestep, last iteration step and actual iteration step for smoother results */
		act_theta_low=(theta_sl_old[k+1]+theta_average[k+1])/2.0; /* averaging thetas of the next lower layer from last timestep, last iteration step and actual iteration step */
		double SE=(SE_sl!=NULL)?SE_sl[k]:1.0; // for k-reduction in frozen soils
		double SE_low=(SE_sl!=NULL)?SE_sl[k+1]:1.0; // for k-reduction in frozen soils

		SE    =(k<1)?(1-dSEeffectiveness*(1-SE)):(1-dSEeffectivenessSoil*(1-SE));
		SE_low=(k<0)?(1-dSEeffectiveness*(1-SE_low)):(1-dSEeffectivenessSoil*(1-SE_low));

		if ((k==clay_layer)||((k+1)==clay_layer))
		{
			keff=k_clay*ST[s].k_calc(k,(theta_sl_old[clay_layer]+theta_average[clay_layer])/2.0,SE, thetares[k], thetasat[k]);
		}
		else
		{
			ko = (k>0)? ST[s].k_calc(k,act_theta,SE, thetares[k], thetasat[k]) : LayeredSnowProperties[1].k_calc(levelsSnow+k,act_theta,SE);
			//k_relativ, oben
			//      if (k<GWLEVEL[i][j]-1)
			//      {
			//        ku=k_calc(ST[s].anz_k,act_theta_low,ST[s].k_tab,ST[s].theta_k); // k_relativ, unten
			//        /* damit bei scharfen Feuchtefronten trotzdem genug fliessen kann,
			//           werden die Unterschiede in den Leitfaehigkeiten begrenzt auf maxratio */
			//        if (ko>ST[s].maxratio*ku) 
			//          ku=ko/ST[s].maxratio;
			//        else if (ku>ST[s].maxratio*ko) 
			//          ko=ku/ST[s].maxratio;
			//        if (ko==0) ko=ku; else if (ku==0) ku=ko; /* wenn beide 0 sind, dann wird unten keff auch 0 */
			//        ko*=ksat[k];
			//        ku*=ksat[k+1];
			//        /*  */
			//        keff = ko==0 ? 0 : ko*ku*(dicken[k]+dicken[k+1])/(dicken[k]*ku+dicken[k+1]*ko);
			//      }
			//      else keff=ko*ksat[k]; /* wenn untere Schicht im GW steht, soll als k-Wert der der darueberliegenden Schicht genommen werden */
			//
			//  neue Version, bei der die Leitfaehigkeit in der letzten vollstaendig ungesaettigten Schicht wie immer 
			//  berechnet wird (also als gewichtetes Mittel aus dieser Scicht und der ersten teilgesaettigten Schicht)
			//  bisher wurde fuer diese vorletzte Schicht immer der k-Wert dieser Schicht genommen, um zu schnelles Absickern zu verhindern. 
			//  Das fuehrt aber bei Aufstiegsbedingungen zu einer hydraulischen Sperre, weil diese unges. Schicht zu trocken ist
			//
			if ((k<=GWLEVEL[i][j]) || (k<= firstExcessIceLevel)) //20230801 und mindestens bis erster Level mit Exessice > 0
			{
				ku= (k>=0) ? ST[s].k_calc(k+1,act_theta_low,SE_low, thetares[k], thetasat[k]) : LayeredSnowProperties[1].k_calc(levelsSnow+k+1,act_theta_low,SE_low); // k_relativ, unten
			
				// to make sure that there is some water movement even in sharp wetting -up fronts, the ratio between higher and lower conductivity is limited to the (optional) maxratio parameter.
				double locMaxRatio = ST[s].getMaxratio();
				if (k<0) {
					locMaxRatio = LayeredSnowProperties[1].getMaxratio();
				} else if (k==0) { // k==0, i.e. the flux is between lowermost snow layer and uppermost soil layer
					if (LayeredSnowProperties[1].getMaxratio()>ST[s].getMaxratio()) locMaxRatio = LayeredSnowProperties[1].getMaxratio();
				}

				if (ko>locMaxRatio*ku) 	{
					ku=ko/locMaxRatio;
				} else	if (ku>locMaxRatio*ko) 	{
					ko=ku/locMaxRatio;
				}

				if (ko==0) {
					ko=ku; 
				} else if (ku==0) {
					ku=ko; /* wenn beide 0 sind, dann wird unten keff auch 0 */
				}
				
				ko*=ksat[k];
				ku*=ksat[k+1];
				
				// geometrisches Mittel 
				// abschalten für die folgende Testsequenz  keff = ko==0 ? 0 : ko*ku*(dicken[k]+dicken[k+1])/(dicken[k]*ku+dicken[k+1]*ko);

				// Test, ob bei k-Wert immer von der besser gesättigten Schicht ein gleichmässigeres Versickern stattfinden kann
				// Das funktioniert in dem Sinn, dass die Abflüsse aus den Schichten nicht mehr so extrem stark schwanken, sondern schön glatt sind.
				// Allerdings reagiert der Boden insgesamt schneller mit Tiefenversickerung, da ja jede Schicht, in die versickert wird, sofort mit vollem
				// K-Wert der darüberliegenden Schicht beliefert wird. Will man die Reaktion verlangsamen, dann muss man den k-Wert des Bodens entweder 
				// über die Bodenartentabelle verringern, oder die Abnahme des k-Wertes mit der Tiefe muss vergrössert werden.
				if (k<GWLEVEL[i][j]-1)
				{
					if (ko>ku)	{
						keff=ko;
					} else	{
						keff=ku;
					}
				}
				else
				{
					if ((ku==0)&&(ko==0))	{
						keff=0;
					}	else	{
						keff = ko*ku*(dicken[k]+dicken[k+1])/(dicken[k]*ku+dicken[k+1]*ko);
					}
				}
			}
			else 
			{
				keff=ksat[k]; /* wenn untere Schicht im GW steht, soll als k-Wert der der darueberliegenden Schicht genommen werden */
			}
		}

		qu_sl[k]=quflow(keff,hh_sl[k],hh_sl[k+1],dicken[k],dicken[k+1]);
		// suppress upwards flux between soil and snow in order to avoid wetting up of the snow which is highly unprobable usually
		if ((k==0)&&(qu_sl[k]>0))
		{
			qu_sl[k]=0; // only uprise is suppressed
		}

//		if (k==0)
//		{
//			qu_sl[0] -= (extraSnowRunoffmm/(DT*1000));
//		}

		/* and now interflow and drainage */
		if ((k<(int)GWLEVEL[i][j]) || (k <= firstExcessIceLevel)) /* only in layers without contact to groundwater */ //20230801 und mindestens bis erstes Level mit ExcessIce > 0
		{
			/* and only, if suction is below 3.45 m (otherwise the relative conductivity is much to small to produce interflow) */
			if (hgeo_sl[k]-hh_sl[k]<3.45)
			{
				/* interflow */
				ifl_min_theta = (k>0)? ST[s].get_theta(k,3.45, thetares[k], thetasat[k]) : LayeredSnowProperties[1].get_theta(levelsSnow+k,3.45); /* getting the theta, interflow may decrease actual theta maximal */
				//        max_ifl = (theta_sl[k]-ifl_min_theta)>0?(theta_sl[k]-ifl_min_theta)*dicken[k]/DT:0; /* maximal interflow is the draiageable theta */
				max_ifl = (theta_average[k]-ifl_min_theta)>0?(theta_average[k]-ifl_min_theta)*dicken[k]/DT:0; /* maximal interflow is the draiageable theta */
				/* only if hydraulic conductivity of the lower layer is smaller than the conductivity
	of the actual layer (otherwise there can be no interflow!) */
				if ((ksat[k]>ksat[k+1]) || ((do_ArcVeg == 1)&&(k >= firstExcessIceLevel))) //20230801 oder ExcessIce > 0
				{
					if (k>0)
					{
						if (do_ArcVeg==0 || ((do_ArcVeg == 1) && ExcessIceFactor[i][j][k]<=0)) //20230801 bei excessIce > 0 TanBeta nicht anwenden, sondern Faktor 1.0
						{
							// interflow within soil depends on soil prperties and drainage density (flussdichte)
							interflow[k] = ST[s].k_calc(k, theta_average[k], SE, thetares[k], thetasat[k])*ksat[k] * dicken[k] * flussdichte[TG];
							/*	dicken[k] is used as a multipliant because interflow is a scpecific flux in the horizontal direction,
								so a profile thickness is needed */
								// if on an Ice wedge, all interflow should be dealed with as melt water that is taken away from the cell
							interflow[k] *= (TANBETA[i][j] > 1 ? 1 : TANBETA[i][j]); /* if the local slope angle is larger than 45 degree,
																the hydraulic gradient is restricted to 1 */
						}
						else //i.e: (do_ArcVeg==1 && ExcessIceFactor[i][j][k]>0)
						{
							// ice wedge melt is not dependent on Flussdichte parameter nor on local slope (gradient)
							double rel_k = ST[s].k_calc(k, theta_average[k], SE, thetares[k], thetasat[k]);
							interflow[k] = rel_k * ksat[k] * dicken[k];
//							interflow[k] = ksat[k] * dicken[k];
						}
					} 
					else
					{	// "interflow" within soil depends on snow properties and not on drainage density (flussdichte)
						interflow[k]= LayeredSnowProperties[1].k_calc(levelsSnow+k,theta_average[k],SE)*ksat[k]*dicken[k];
						interflow[k]*= dDrainageDensitiySnow;
					}
				}
				else interflow[k]=0.0;

				/* drainage */
				if (HDRAIN)
				if (k==DRAINLEVEL[i][j]) // applies to levels within real soil only, not to snow pack layers
				{
					qdrain=ST[s].k_calc(k,theta_average[k],SE)*2*ksat[k]*dicken[k]/draindist*cellsize/draindist;
					// flow clipping for drainage: no more drainage is allowed than is drainable from the actual layer
					// and this is defined by the suctions
					if (qdrain>max_ifl)
					{
						qdrain=max_ifl;
					}
					/* 2*dicken[k]/draindist is a kind of hydraulic gradient, cellsize/draindist is the number of drainage pipes within a cell */
					// flow clipping for total interflow+drainage -> drainage will suffer as last one
					interflow[k]+=qdrain;
				}
				if (max_ifl<interflow[k]) interflow[k]=max_ifl;
			}
			else interflow[k]=0.0;
		}
		else interflow[k]=0.0;

	}  /* Berechnen des Flusses in der naechsten Schicht */

	for (int k=1+offset;k<=l && k<levels;k++) /* Berechnen der Gradientenrestriktionen */
	{
		if (k!=0) // if layeredSnow is active, k starts from values at or below zero. But the inflow into the first layer should not be restricetd by this routine since all precipitation or melt should be "pressed" into the ground
		{
			qu_sl[k]=gradient_restriction(k,s,theta_sl_old[k],theta_sl_old[k+1],qu_sl[k-1],qu_sl[k],qu_sl[k+1],hgeo_sl,interflow,hh_sl,DT,SE_sl,thetares,thetasat);
		}
	}

	delete[] m_ksat;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::flux_restrictions(int offset,int i,int j,int anz,
int &maxlevel,double* &gwin,double* &gwex,
double* &theta_sl_old,double* &theta_temp,double* &thetasat,double* &thetares,
double* &qu_sl,double* &interflow,double* &dicken,
double &GWINtotal,double &GWEXtotal,double DT,double *SE_sl)
/* balancing the inflows and outflows into and out of the layers and checking,
if saturation would be exceeded (over saturation) or if moisture would fall below wilting point (drying out) */
{
	int k,l,z,start,flag;
	double q_temp1;
	/* again, maxlevel is the deepest layer, from which transpiration tookes places */
	l=maxlevel > (int)GWLEVEL[i][j]+3 ? maxlevel : (int)GWLEVEL[i][j]+3; /* if maxlevel is below groundwater level,
																	the whole schmenis done to maxlevel, elseo to groundwater depth*/
	if (l<levels)
	{
		while (gwin[l+1]>0||gwex[l+1]>0) 
		{
			l++;
			if (l>=levels) break;
		}
	}
	l+=2; /* to avoid fluxes within the saturated zone, at least the first both saturated layers should also be included 
(one layer for considering a change in groundwater level and one for considering the flux between this changed layer and the next layer */
	if (l>levels) start=levels; else start=l;

	z = 0;
	do
	{
		z++;
		flag=0;
		for (int k=start;k>=1+offset;k--) /* upwards is faster */
		{
			/* new theta as balance of the boundary fluxes and the sinks (interflow contains also drainage);
	extraction by evaporation is already done */
			theta_temp[k] = (qu_sl[k]-qu_sl[k-1]-interflow[k]-gwex[k]+gwin[k])*DT;
			theta_temp[k] /= dicken[k];
			theta_temp[k] += theta_sl_old[k];

			if (theta_temp[k]>thetasat[k]) /* avoiding oversaturation */
			{ /* reducing inflows -> qu_sl[k-1] is an inflow, if it is <0, 
	qu_sl[k] is an inflow, if it is >0, gwin[k] is an inflow, gwex[k] and interflow[k] are always outflows! */
				if ((qu_sl[k-1]<=0)&&(qu_sl[k]<=0)) /* flowing through downwards */
				{  /* reducing inflow from above */
					q_temp1 = (theta_temp[k]-thetasat[k])*dicken[k]/DT; /* amount of to much inflow */
					qu_sl[k-1]+=q_temp1;
					if (qu_sl[k-1]>0) 			                
					{
						// now, both fluxes flow out of the layer, so there must be some other inflow into the layer,	
						// which can come from GWIN only (or from numerical uncertainties...)
						if (gwflow_active && (surface_interaction>1) && gwin[k]>0)
						{
							if (gwin[k]>=qu_sl[k-1])
							{
								gwin[k]-=qu_sl[k-1]; 
								GWINtotal-=qu_sl[k-1]; 
								qu_sl[k-1]=0;
							}
							else
							{
								qu_sl[k-1]-=gwin[k]; 
								GWINtotal-=gwin[k]; 
								gwin[k]=0;
							}
						}
					}
					flag=1;
				}
				else
				if ((qu_sl[k-1]<=0)&&(qu_sl[k]>0)) /* inflow from upper and lower layer  */
				{ /* reducing both fluxes */
					q_temp1 = (theta_temp[k]-thetasat[k])*dicken[k]/DT; /* amount of to much inflow */
					qu_sl[k-1]+=q_temp1;  /* first reducing inflow from above */
					if (qu_sl[k-1]>0)         /* if it waszoomuch reduction */
					{
						qu_sl[k]-=qu_sl[k-1];   /* reducing the inflow from below by the overhead */
						qu_sl[k-1]=0;           /* and set the inflow from above to 0 */
						if (qu_sl[k]<0)         /* if this was still too much */
						{
							if (gwflow_active && (surface_interaction>1) && gwin[k]>0)
							{
								if (gwin[k]+=qu_sl[k]>=0)
								{
									gwin[k]+=qu_sl[k];    /* reducing the infiltration --> qu_sl[k] is negative! */
									GWINtotal+=qu_sl[k];
									qu_sl[k]=0;           /* and set the inflow from below to 0 */
								}
								else
								{
									qu_sl[k]+=gwin[k];
									GWINtotal+=gwin[k];
									gwin[k]=0;
								}
							}
						}
					}
					flag=1;
				}
				else
				if ((qu_sl[k-1]>0)&&(qu_sl[k]>=0)) /* flowing through upwards */
				{ /* inflow from lower layer is reduced, outflow remains constant */
					q_temp1 = (theta_temp[k]-thetasat[k])*dicken[k]/DT; /* amount of to much inflow */
					qu_sl[k]-=q_temp1;
					if (qu_sl[k]<0) 
					{
						if (gwflow_active && (surface_interaction>1) && gwin[k]>0)
						{
							if (gwin[k]+=qu_sl[k]>=0)
							{
								gwin[k]+=qu_sl[k];    /* reducing the infiltration --> qu_sl[k] is negative! */
								GWINtotal+=qu_sl[k];
								qu_sl[k]=0;           /* and set the inflow from below to 0 */
							}
							else
							{
								qu_sl[k]+=gwin[k];
								GWINtotal+=gwin[k];
								gwin[k]=0;
							}
						}
					}
					flag=1;
				}
				else
				if ((qu_sl[k-1]>=0)&&(qu_sl[k]<=0)) /* outflow up- and down, but also infiltration */
				{ /* only infiltration is reduced, outflow remains constant */
					q_temp1     = (theta_temp[k]-thetasat[k])*dicken[k]/DT; /* amount of to much inflow */
					if (gwflow_active && (surface_interaction>1) && gwin[k]>0)
					{
						if (gwin[k]>=q_temp1)
						{
							gwin[k]-=q_temp1; 
							GWINtotal-=q_temp1; 
						}
						else
						{
							GWINtotal-=gwin[k]; 
							gwin[k]=0;
						}
					}
					flag=1;
				}
			}
			else /* avoiding drying out */
			{
				if (theta_temp[k]<thetares[k])
				{
					double SE=(SE_sl!=NULL)?SE_sl[k]:1.0; // for (partly) frozen soils
					SE=(k<1)?(1-dSEeffectiveness*(1-SE)):(1-dSEeffectivenessSoil*(1-SE));

					double local_thetares= thetares[k]+(theta_sl_old[k]-thetares[k])*(1.0-SE); // if the soil is frozen or we are in a snow layer, the frozen fraction of water reduces the available water for fluxes. SO we put this frozen water temporarily on top of the residual water content
					
					if ((qu_sl[k-1]<=0)&&(qu_sl[k]<=0))  /* flowing through downwards */
					{ /* inflow remains constant, outflows are reduced */
						q_temp1   = local_thetares<theta_sl_old[k]?
						(local_thetares - theta_temp[k])*dicken[k]/DT: /* fraction of to much outflow */
						(((theta_sl_old[k] - theta_temp[k])>0) ?
						(theta_sl_old[k] - theta_temp[k])*dicken[k] / DT : -qu_sl[k]); /* fraction of to much outflow if theta was already below wilting point*/
						qu_sl[k] += q_temp1; /* reducing outflow downwards by adding a positive number to a negative outflow */
						if (qu_sl[k]>0) /* if it was too much compensation */
						{
							if (interflow[k]>0)
							{
								if (interflow[k]>=qu_sl[k])
								{
									interflow[k] -= qu_sl[k];  /* reducing interflow */
									qu_sl[k]=0;                /* and setting outflow downwards to 0 */
								}
								else
								{
									qu_sl[k]-=interflow[k];
									interflow[k]=0;
								}
							}
							else                         /* or, if gwex was calculated in this layer */
							{
								if (gwflow_active && (surface_interaction>0) && gwex[k]>0)
								{
									if (gwex[k]>=qu_sl[k])
									{   // there was more or equal as much exfiltration than the reduction results on additional inflow from below
										gwex[k] -= qu_sl[k];       /* doing the same with groundwater exfiltration */
										GWEXtotal -= qu_sl[k];
										qu_sl[k]=0;
									}
									else
									{   // inflow from below is larger than the exfiltration amount
										qu_sl[k]-=gwex[k];
										GWEXtotal -= gwex[k];
										gwex[k]=0;
									}
								}
//								else qu_sl[k]-=q_temp1;
								else qu_sl[k]=0; // 
							}
						}
						flag=1;
					}
					else
					{
						if ((qu_sl[k-1]>0)&&(qu_sl[k]<=0)) /* outflow upwards and downwards */
						{ /* reducing fluxes up- and downward and interflow */
							q_temp1   = local_thetares<theta_sl_old[k]?
							(local_thetares - theta_temp[k])*dicken[k]/DT: /* fraction of to much outflow if theta (before flux) was not below wilting point */
							(((theta_sl_old[k] - theta_temp[k])>0) ?
							(theta_sl_old[k] - theta_temp[k])*dicken[k] / DT : -qu_sl[k]); /* fraction of to much outflow if theta was already below wilting point*/
							qu_sl[k-1] -= q_temp1; /* reducing outflow upwards by subtracting a positive number from a positive outflow (upwards is positive) */
							if (qu_sl[k-1]<0)       /* if it was too much reduction -> flow is reversed now */
							{
								qu_sl[k]-=qu_sl[k-1]; /* reducing outflow downwards by the exceeding amount */
								qu_sl[k-1]=0;         /* and set the outflow upwards to 0 */
								if (qu_sl[k]>0) /* if it was too much */
								{
									if (interflow[k]>0)
									{
										if (interflow[k]>=qu_sl[k])
										{
											interflow[k] -= qu_sl[k];  /* reducing interflow */
											qu_sl[k]=0;                /* and setting outflow downwards to 0 */
										}
										else
										{
											qu_sl[k]-=interflow[k];
											interflow[k]=0;
										}
									}
									else                         /* or, if gwex was calculated in this layer */
									if (gwflow_active && (surface_interaction>0) && gwex[k]>0)
									{
										if (gwex[k]>=qu_sl[k])
										{   // there was more or equal as much exfiltration than the reduction results on additional inflow from below
											gwex[k] -= qu_sl[k];       /* doing the same with groundwater exfiltration */
											GWEXtotal -= qu_sl[k];
											qu_sl[k]=0;
										}
										else
										{   // inflow from below is larger than the exfiltration amount
											qu_sl[k]-=gwex[k];
											GWEXtotal -= gwex[k];
											gwex[k]=0;
										}
									}
								}
							}
							flag=1;
						}
						else
						{
							if ((qu_sl[k-1]<=0)&&(qu_sl[k]>0)) /* inflow from both, lower and upper layer */
							{ /* inflows from upper and lower layer, but interflow or exfiltration are too large -> reducing interflow or exfiltration to avoid drying out! */
								q_temp1   = local_thetares<theta_sl_old[k]?
								(local_thetares - theta_temp[k])*dicken[k]/DT: /* fraction of to much outflow if theta (before flux) was not below wilting point */
								(((theta_sl_old[k] - theta_temp[k])>0) ?
								(theta_sl_old[k] - theta_temp[k])*dicken[k] / DT : -qu_sl[k]); /* fraction of to much outflow if theta was already below wilting point*/
								//                q_temp1      = (thetares - theta_temp[k])*dicken[k]/DT; /* fraction of to much inflow */
								if (interflow[k]>0)
								{
									if (interflow[k]-q_temp1>0)
									{
										interflow[k] -= q_temp1;
									}
									else
									{
										interflow[k]=0;
									}
								}
								else 
								{
									if (gwflow_active && (surface_interaction>0) && gwex[k]>0)
									{
										if (gwex[k]>=q_temp1)
										{
											gwex[k] -= q_temp1;
											GWEXtotal -= q_temp1;
										}
										else
										{
											GWEXtotal -= gwex[k];
											gwex[k] = 0;
										}
									}
								}
								flag=1;
							}
							else
							{
								if ((qu_sl[k-1]>0)&&(qu_sl[k]>0)) /* flowing through upwards */
								{ /* upwards outflow is reduced, inflow remains const */
									q_temp1   = local_thetares<theta_sl_old[k]?
									(local_thetares - theta_temp[k])*dicken[k]/DT: /* fraction of to much outflow if theta (before flux) was not below wilting point */
									(((theta_sl_old[k] - theta_temp[k])>0) ?
									(theta_sl_old[k] - theta_temp[k])*dicken[k] / DT : -qu_sl[k]); /* fraction of to much outflow if theta was already below wilting point*/
									//                  q_temp1   = (thetares - theta_temp[k])*dicken[k]/DT; /* fraction of to much outflow */
									qu_sl[k-1] -= q_temp1; /* reducing outflow upwards by adding a negative number to a positive outflow */
									if (qu_sl[k-1]<0) /* if it was too much */
									{
										if (interflow[k]>0)
										{
											if (interflow[k] + qu_sl[k-1] >= 0)
											{
												interflow[k] += qu_sl[k-1];  /* reducing interflow */
												qu_sl[k-1]=0;                /* and setting outflow downwards to 0 */
											}
											else
											{
												qu_sl[k-1]+=interflow[k];
												interflow[k]=0;
											}
										}
										else                         /* or, if gwex was calculated in this layer */
										if (gwflow_active && (surface_interaction>0) && gwex[k]>0)
										{
											if (gwex[k]+qu_sl[k-1]>0)
											{
												gwex[k]   += qu_sl[k-1];       /* doing the same with groundwater exfiltration */
												GWEXtotal += qu_sl[k-1];
												qu_sl[k-1] = 0;
											}
											else // gwex[k]+qu_sl[k-1]<0
											{
												qu_sl[k-1]+=gwex[k];
												GWEXtotal -=gwex[k];
												gwex[k] = 0;
											}
										}
									}
									flag=1;
								}
							}
						}
					}
				}
				if (interflow[k]<0) interflow[k]=0;
			}
		}
	} while ((flag==1)&&(z<anz));
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::flux_restrictions_opt(int offset, int i, int j, int anz,
	int &maxlevel, double* &gwin, double* &gwex,
	double* &theta_sl_old, double* &theta_temp, double* &thetasat, double* &thetares,
	double* &qu_sl, double* &interflow, double* &dicken,
	double &GWINtotal, double &GWEXtotal, double DT, double *SE_sl)
	/* balancing the inflows and outflows into and out of the layers and checking,
	if saturation would be exceeded (over saturation) or if moisture would fall below wilting point (drying out) */
{
	int k, l, z, start, flag;
	double q_temp1;
	/* again, maxlevel is the deepest layer, from which transpiration tookes places */
	l = maxlevel > (int)GWLEVEL[i][j] + 3 ? maxlevel : (int)GWLEVEL[i][j] + 3; /* if maxlevel is below groundwater level,
																	the whole schmenis done to maxlevel, elseo to groundwater depth*/
	if (l < levels)
	{
		while (gwin[l + 1] > 0 || gwex[l + 1] > 0)
		{
			l++;
			if (l >= levels) break;
		}
	}
	l += 2; /* to avoid fluxes within the saturated zone, at least the first both saturated layers should also be included
(one layer for considering a change in groundwater level and one for considering the flux between this changed layer and the next layer */
	if (l > levels) start = levels; else start = l;

	double theta_temp_loc[200];
	double qu_sl_loc[200];
	double interflow_loc[200];
	double gwex_loc[200];
	double gwin_loc[200];
	double thetasat_loc[200];
	double dicken_loc[200];
	double theta_sl_old_loc[200];
	double thetares_loc[200];
	double SE_sl_loc[200];

	int sz = sizeof(double)*(start - offset + 1);

	// for performance optimization, heap arrays are moved to local stack arrays (works with up to 200 layers for soil (including snow layers) only
	memmove(&(theta_temp_loc[0]), &(theta_temp[offset]), sz);
	memmove(&(qu_sl_loc[0]), &(qu_sl[offset]), sz);
	memmove(&(interflow_loc[0]), &(interflow[offset]), sz);
	memmove(&(gwex_loc[0]), &(gwex[offset]), sz);
	memmove(&(gwin_loc[0]), &(gwin[offset]), sz);
	memmove(&(thetasat_loc[0]), &(thetasat[offset]), sz);
	memmove(&(dicken_loc[0]), &(dicken[offset]), sz);
	memmove(&(theta_sl_old_loc[0]), &(theta_sl_old[offset]), sz);
	memmove(&(thetares_loc[0]), &(thetares[offset]), sz);
	if (SE_sl) memmove(&(SE_sl_loc[0]), &(SE_sl[offset]), sz);

	int startorig = start;
	int offsetorig = offset;
	if (offset < 0)
	{
		start -= offset;
		offset = 0;
	}

	z = 0;
	do
	{
		z++;
		flag = 0;
		for (int k = start; k >= 1 + offset; k--) /* upwards is faster */
		{
			/* new theta as balance of the boundary fluxes and the sinks (interflow contains also drainage);
	extraction by evaporation is already done */
			theta_temp_loc[k] = (qu_sl_loc[k] - qu_sl_loc[k - 1] - interflow_loc[k] - gwex_loc[k] + gwin_loc[k])*DT;
			theta_temp_loc[k] /= dicken_loc[k];
			theta_temp_loc[k] += theta_sl_old_loc[k];

			if (theta_temp_loc[k] > thetasat_loc[k]) /* avoiding oversaturation */
			{ /* reducing inflows -> qu_sl_loc[k-1] is an inflow, if it is <0,
	qu_sl_loc[k] is an inflow, if it is >0, gwin_loc[k] is an inflow, gwex_loc[k] and interflow_loc[k] are always outflows! */
				if ((qu_sl_loc[k - 1] <= 0) && (qu_sl_loc[k] <= 0)) /* flowing through downwards */
				{  /* reducing inflow from above */
					q_temp1 = (theta_temp_loc[k] - thetasat_loc[k])*dicken_loc[k] / DT; /* amount of to much inflow */
					qu_sl_loc[k - 1] += q_temp1;
					if (qu_sl_loc[k - 1] > 0)
					{
						// now, both fluxes flow out of the layer, so there must be some other inflow into the layer,	
						// which can come from GWIN only (or from numerical uncertainties...)
						if (gwflow_active && (surface_interaction > 1) && gwin_loc[k] > 0)
						{
							if (gwin_loc[k] >= qu_sl_loc[k - 1])
							{
								gwin_loc[k] -= qu_sl_loc[k - 1];
								GWINtotal -= qu_sl_loc[k - 1];
								qu_sl_loc[k - 1] = 0;
							}
							else
							{
								qu_sl_loc[k - 1] -= gwin_loc[k];
								GWINtotal -= gwin_loc[k];
								gwin_loc[k] = 0;
							}
						}
					}
					flag = 1;
				}
				else
					if ((qu_sl_loc[k - 1] <= 0) && (qu_sl_loc[k] > 0)) /* inflow from upper and lower layer  */
					{ /* reducing both fluxes */
						q_temp1 = (theta_temp_loc[k] - thetasat_loc[k])*dicken_loc[k] / DT; /* amount of to much inflow */
						qu_sl_loc[k - 1] += q_temp1;  /* first reducing inflow from above */
						if (qu_sl_loc[k - 1] > 0)         /* if it waszoomuch reduction */
						{
							qu_sl_loc[k] -= qu_sl_loc[k - 1];   /* reducing the inflow from below by the overhead */
							qu_sl_loc[k - 1] = 0;           /* and set the inflow from above to 0 */
							if (qu_sl_loc[k] < 0)         /* if this was still too much */
							{
								if (gwflow_active && (surface_interaction > 1) && gwin_loc[k] > 0)
								{
									if (gwin_loc[k] += qu_sl_loc[k] >= 0)
									{
										gwin_loc[k] += qu_sl_loc[k];    /* reducing the infiltration --> qu_sl_loc[k] is negative! */
										GWINtotal += qu_sl_loc[k];
										qu_sl_loc[k] = 0;           /* and set the inflow from below to 0 */
									}
									else
									{
										qu_sl_loc[k] += gwin_loc[k];
										GWINtotal += gwin_loc[k];
										gwin_loc[k] = 0;
									}
								}
							}
						}
						flag = 1;
					}
					else
						if ((qu_sl_loc[k - 1] > 0) && (qu_sl_loc[k] >= 0)) /* flowing through upwards */
						{ /* inflow from lower layer is reduced, outflow remains constant */
							q_temp1 = (theta_temp_loc[k] - thetasat_loc[k])*dicken_loc[k] / DT; /* amount of to much inflow */
							qu_sl_loc[k] -= q_temp1;
							if (qu_sl_loc[k] < 0)
							{
								if (gwflow_active && (surface_interaction > 1) && gwin_loc[k] > 0)
								{
									if (gwin_loc[k] += qu_sl_loc[k] >= 0)
									{
										gwin_loc[k] += qu_sl_loc[k];    /* reducing the infiltration --> qu_sl_loc[k] is negative! */
										GWINtotal += qu_sl_loc[k];
										qu_sl_loc[k] = 0;           /* and set the inflow from below to 0 */
									}
									else
									{
										qu_sl_loc[k] += gwin_loc[k];
										GWINtotal += gwin_loc[k];
										gwin_loc[k] = 0;
									}
								}
							}
							flag = 1;
						}
						else
							if ((qu_sl_loc[k - 1] >= 0) && (qu_sl_loc[k] <= 0)) /* outflow up- and down, but also infiltration */
							{ /* only infiltration is reduced, outflow remains constant */
								q_temp1 = (theta_temp_loc[k] - thetasat_loc[k])*dicken_loc[k] / DT; /* amount of to much inflow */
								if (gwflow_active && (surface_interaction > 1) && gwin_loc[k] > 0)
								{
									if (gwin_loc[k] >= q_temp1)
									{
										gwin_loc[k] -= q_temp1;
										GWINtotal -= q_temp1;
									}
									else
									{
										GWINtotal -= gwin_loc[k];
										gwin_loc[k] = 0;
									}
								}
								flag = 1;
							}
			}
			else /* avoiding drying out */
			{
				if (theta_temp_loc[k] < thetares_loc[k])
				{
					double SE = (SE_sl != NULL) ? SE_sl_loc[k] : 1.0; // for (partly) frozen soils
					SE = (k < 1) ? (1 - dSEeffectiveness * (1 - SE)) : (1 - dSEeffectivenessSoil * (1 - SE));

					double local_thetares = thetares_loc[k] + (theta_sl_old_loc[k] - thetares_loc[k])*(1.0 - SE); // if the soil is frozen or we are in a snow layer, the frozen fraction of water reduces the available water for fluxes. SO we put this frozen water temporarily on top of the residual water content

					if ((qu_sl_loc[k - 1] <= 0) && (qu_sl_loc[k] <= 0))  /* flowing through downwards */
					{ /* inflow remains constant, outflows are reduced */
						q_temp1 = local_thetares < theta_sl_old_loc[k] ?
							(local_thetares - theta_temp_loc[k])*dicken_loc[k] / DT : /* fraction of to much outflow */
							(((theta_sl_old_loc[k] - theta_temp_loc[k]) > 0) ?
							(theta_sl_old_loc[k] - theta_temp_loc[k])*dicken_loc[k] / DT : -qu_sl_loc[k]); /* fraction of to much outflow if theta was already below wilting point*/
						qu_sl_loc[k] += q_temp1; /* reducing outflow downwards by adding a positive number to a negative outflow */
						if (qu_sl_loc[k] > 0) /* if it was too much compensation */
						{
							if (interflow_loc[k] > 0)
							{
								if (interflow_loc[k] >= qu_sl_loc[k])
								{
									interflow_loc[k] -= qu_sl_loc[k];  /* reducing interflow */
									qu_sl_loc[k] = 0;                /* and setting outflow downwards to 0 */
								}
								else
								{
									qu_sl_loc[k] -= interflow_loc[k];
									interflow_loc[k] = 0;
								}
							}
							else                         /* or, if gwex was calculated in this layer */
							{
								if (gwflow_active && (surface_interaction > 0) && gwex_loc[k] > 0)
								{
									if (gwex_loc[k] >= qu_sl_loc[k])
									{   // there was more or equal as much exfiltration than the reduction results on additional inflow from below
										gwex_loc[k] -= qu_sl_loc[k];       /* doing the same with groundwater exfiltration */
										GWEXtotal -= qu_sl_loc[k];
										qu_sl_loc[k] = 0;
									}
									else
									{   // inflow from below is larger than the exfiltration amount
										qu_sl_loc[k] -= gwex_loc[k];
										GWEXtotal -= gwex_loc[k];
										gwex_loc[k] = 0;
									}
								}
								//								else qu_sl_loc[k]-=q_temp1;
								else qu_sl_loc[k] = 0; // 
							}
						}
						flag = 1;
					}
					else
					{
						if ((qu_sl_loc[k - 1] > 0) && (qu_sl_loc[k] <= 0)) /* outflow upwards and downwards */
						{ /* reducing fluxes up- and downward and interflow */
							q_temp1 = local_thetares < theta_sl_old_loc[k] ?
								(local_thetares - theta_temp_loc[k])*dicken_loc[k] / DT : /* fraction of to much outflow if theta (before flux) was not below wilting point */
								(((theta_sl_old_loc[k] - theta_temp_loc[k]) > 0) ?
								(theta_sl_old_loc[k] - theta_temp_loc[k])*dicken_loc[k] / DT : -qu_sl_loc[k]); /* fraction of to much outflow if theta was already below wilting point*/
							qu_sl_loc[k - 1] -= q_temp1; /* reducing outflow upwards by subtracting a positive number from a positive outflow (upwards is positive) */
							if (qu_sl_loc[k - 1] < 0)       /* if it was too much reduction -> flow is reversed now */
							{
								qu_sl_loc[k] -= qu_sl_loc[k - 1]; /* reducing outflow downwards by the exceeding amount */
								qu_sl_loc[k - 1] = 0;         /* and set the outflow upwards to 0 */
								if (qu_sl_loc[k] > 0) /* if it was too much */
								{
									if (interflow_loc[k] > 0)
									{
										if (interflow_loc[k] >= qu_sl_loc[k])
										{
											interflow_loc[k] -= qu_sl_loc[k];  /* reducing interflow */
											qu_sl_loc[k] = 0;                /* and setting outflow downwards to 0 */
										}
										else
										{
											qu_sl_loc[k] -= interflow_loc[k];
											interflow_loc[k] = 0;
										}
									}
									else                         /* or, if gwex was calculated in this layer */
										if (gwflow_active && (surface_interaction > 0) && gwex_loc[k] > 0)
										{
											if (gwex_loc[k] >= qu_sl_loc[k])
											{   // there was more or equal as much exfiltration than the reduction results on additional inflow from below
												gwex_loc[k] -= qu_sl_loc[k];       /* doing the same with groundwater exfiltration */
												GWEXtotal -= qu_sl_loc[k];
												qu_sl_loc[k] = 0;
											}
											else
											{   // inflow from below is larger than the exfiltration amount
												qu_sl_loc[k] -= gwex_loc[k];
												GWEXtotal -= gwex_loc[k];
												gwex_loc[k] = 0;
											}
										}
								}
							}
							flag = 1;
						}
						else
						{
							if ((qu_sl_loc[k - 1] <= 0) && (qu_sl_loc[k] > 0)) /* inflow from both, lower and upper layer */
							{ /* inflows from upper and lower layer, but interflow or exfiltration are too large -> reducing interflow or exfiltration to avoid drying out! */
								q_temp1 = local_thetares < theta_sl_old_loc[k] ?
									(local_thetares - theta_temp_loc[k])*dicken_loc[k] / DT : /* fraction of to much outflow if theta (before flux) was not below wilting point */
									(((theta_sl_old_loc[k] - theta_temp_loc[k]) > 0) ?
									(theta_sl_old_loc[k] - theta_temp_loc[k])*dicken_loc[k] / DT : -qu_sl_loc[k]); /* fraction of to much outflow if theta was already below wilting point*/
									//                q_temp1      = (thetares - theta_temp_loc[k])*dicken_loc[k]/DT; /* fraction of to much inflow */
								if (interflow_loc[k] > 0)
								{
									if (interflow_loc[k] - q_temp1 > 0)
									{
										interflow_loc[k] -= q_temp1;
									}
									else
									{
										interflow_loc[k] = 0;
									}
								}
								else
								{
									if (gwflow_active && (surface_interaction > 0) && gwex_loc[k] > 0)
									{
										if (gwex_loc[k] >= q_temp1)
										{
											gwex_loc[k] -= q_temp1;
											GWEXtotal -= q_temp1;
										}
										else
										{
											GWEXtotal -= gwex_loc[k];
											gwex_loc[k] = 0;
										}
									}
								}
								flag = 1;
							}
							else
							{
								if ((qu_sl_loc[k - 1] > 0) && (qu_sl_loc[k] > 0)) /* flowing through upwards */
								{ /* upwards outflow is reduced, inflow remains const */
									q_temp1 = local_thetares < theta_sl_old_loc[k] ?
										(local_thetares - theta_temp_loc[k])*dicken_loc[k] / DT : /* fraction of to much outflow if theta (before flux) was not below wilting point */
										(((theta_sl_old_loc[k] - theta_temp_loc[k]) > 0) ?
										(theta_sl_old_loc[k] - theta_temp_loc[k])*dicken_loc[k] / DT : -qu_sl_loc[k]); /* fraction of to much outflow if theta was already below wilting point*/
										//                  q_temp1   = (thetares - theta_temp_loc[k])*dicken_loc[k]/DT; /* fraction of to much outflow */
									qu_sl_loc[k - 1] -= q_temp1; /* reducing outflow upwards by adding a negative number to a positive outflow */
									if (qu_sl_loc[k - 1] < 0) /* if it was too much */
									{
										if (interflow_loc[k] > 0)
										{
											if (interflow_loc[k] + qu_sl_loc[k - 1] >= 0)
											{
												interflow_loc[k] += qu_sl_loc[k - 1];  /* reducing interflow */
												qu_sl_loc[k - 1] = 0;                /* and setting outflow downwards to 0 */
											}
											else
											{
												qu_sl_loc[k - 1] += interflow_loc[k];
												interflow_loc[k] = 0;
											}
										}
										else                         /* or, if gwex was calculated in this layer */
											if (gwflow_active && (surface_interaction > 0) && gwex_loc[k] > 0)
											{
												if (gwex_loc[k] + qu_sl_loc[k - 1] > 0)
												{
													gwex_loc[k] += qu_sl_loc[k - 1];       /* doing the same with groundwater exfiltration */
													GWEXtotal += qu_sl_loc[k - 1];
													qu_sl_loc[k - 1] = 0;
												}
												else // gwex_loc[k]+qu_sl_loc[k-1]<0
												{
													qu_sl_loc[k - 1] += gwex_loc[k];
													GWEXtotal -= gwex_loc[k];
													gwex_loc[k] = 0;
												}
											}
									}
									flag = 1;
								}
							}
						}
					}
				}
				if (interflow_loc[k] < 0) interflow_loc[k] = 0;
			}
		}
	} while ((flag == 1) && (z < anz));

	// restore indices
	start = startorig;
	offset = offsetorig;

	// move array from stack to heap again
	memmove(&(theta_temp[offset])	, &(theta_temp_loc[0])	, sz);
	memmove(&(qu_sl[offset])		, &(qu_sl_loc[0])		, sz);
	memmove(&(interflow[offset])	, &(interflow_loc[0])	, sz);
	memmove(&(gwex[offset])			, &(gwex_loc[0])		, sz);
	memmove(&(gwin[offset])			, &(gwin_loc[0])		, sz);
	memmove(&(thetasat[offset])		, &(thetasat_loc[0])	, sz);
//	memmove(&(theta_sl_old[offset])	, &(theta_sl_old_loc[0]), sz); // values have not changed, so re-move is not neccessary
//	memmove(&(dicken[offset])		, &(dicken_loc[0])		, sz); // values have not changed, so re-move is not unneccessary
//	memmove(&(thetares[offset])		, &(thetares_loc[0])	, sz); // values have not changed, so re-move is not neccessary
	if (SE_sl) memmove(&(SE_sl[offset])		, &(SE_sl_loc[0])		, sz);

}
//#pragma optimize ("p",off)


//#pragma optimize ("p",on)
double unsatzonmodellklasse::new_thetas(int offset,int i,int j,int z,
int &maxlevel,double* &gwin,double* &gwex,double* &theta_sl,
double* &theta_average,double* &theta_sl_old,double* &qu_sl,
double* &interflow,double* &dicken,long &weight,double DT)
{
	int k,l,v;
	double th_error,temp,theta_last_iter;
	th_error=0;  v=3;

	if (z>v) weight=(z-v);

	/* maxlevel: the level to which evapotranspiration took place -> this may be below groundwater table! */
	l=maxlevel > (int)GWLEVEL[i][j] ? maxlevel : (int)GWLEVEL[i][j];
	/* also groundwater exfiltration can occure only below or at groundwater table, therefore these layers have to be considered in the scheme */
	if (l<levels)
	{
		while (gwin[l+1]>0||gwex[l+1]>0) 
		{
			l++;
			if (l>levels) break;
		}
	}
	l+=2; /* to avoid fluxes within the saturated zone, at least the first both saturated layers should also be included 
(one layer for considering a change in groundwater level and one for considering the flux between this changed layer and the next layer */
	if (l>levels) l=levels; 

	for (int k=1+offset;k<=l;k++)
	{
		theta_last_iter=theta_sl[k];
		if (z<=v)
		theta_average[k]=theta_sl[k]; /* in the first v iteration steps, the moving average is replaced by only the last theta value */

		theta_sl[k]=theta_sl_old[k]+(qu_sl[k]-qu_sl[k-1]-interflow[k]-gwex[k]+gwin[k])*DT/dicken[k]; /* balancing inflows and outflows */
		if (theta_sl[k]<0.0)
		{
			theta_sl[k]=0.0;
		}
		else if (theta_sl[k]>1.0)
		{
			theta_sl[k]=1.0;
		}
		temp=fabs((double)theta_last_iter-(double)theta_sl[k]); /* change from the old theta */
		if (temp>th_error) th_error=temp;
		if (z>v)  /* moving average */
		{
			theta_average[k]*=weight;     /* calculating new moving average for the iteration thetas */
			theta_average[k]+=theta_sl[k];
			theta_average[k]/=(weight+1);
		}
		else
		{
			theta_average[k]+=theta_sl[k];
			theta_average[k]/=2.0;
		}
	}
	return th_error;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double unsatzonmodellklasse::groundwater_recharge(int i,int j,double gw_exact,
int &old_start_level,int &start_level,double &old_gwstand,
double* &hgeo_sl,double* &theta_sl,double* &thetasat,
double* &theta_sl_old,double* &lateralinflow,double* &dicken,double* &qu_sl,
double &GWINtotal,double &GWEXtotal,double DT,double timestep)
{
	int s,k,n;
	double gwn;
	double thetaup_old,thetaup_new,hdiff;
	
	if (!m_bCoupled)
	{
		s=(int)SOIL[i][j];
		k=old_start_level;
		if (start_level==old_start_level) /* groundwater table remains within thesame level */
		{
			if (old_gwstand>=gw_exact) /* falling groundwater table within the same layer */
			{
				hdiff=(hgeo_sl[k-1]-(double)gw_exact);
				if (hdiff!=0.0)
				{
					//thetaup_new = (ST[s].getDzui()[k] * theta_sl[k] - (gw_exact - hgeo_sl[k])*thetasat[k]) / hdiff;
					thetaup_new = (dicken[k] * theta_sl[k] - (gw_exact - hgeo_sl[k])*thetasat[k]) / hdiff;
					gwn = (gw_exact-old_gwstand)*(thetasat[k]-thetaup_new)*1000.0; /* recharge in mm for statistik files */
				} else gwn=0.0;
			}
			else  /* rising groundwater table within the same layer */
			{
				hdiff=(hgeo_sl[k-1]-old_gwstand);
				if (hdiff!=0.0)
				{
//					thetaup_old = (ST[s].getDzui()[k] * theta_sl_old[k] - (old_gwstand - hgeo_sl[k])*thetasat[k]) / (hgeo_sl[k - 1] - old_gwstand);
					thetaup_old = (dicken[k] * theta_sl_old[k] - (old_gwstand - hgeo_sl[k])*thetasat[k]) / (hgeo_sl[k - 1] - old_gwstand);
					gwn = (gw_exact-old_gwstand)*(thetasat[k]-thetaup_old)*1000.0; /* recharge in mm for statistik files */
				} else gwn=0.0;
			}
		}
		else 
		if (start_level>old_start_level) /* falling groundwater table (passing a layer boundary) */
		{
			gwn  = (hgeo_sl[old_start_level]-old_gwstand)*(thetasat[old_start_level]-theta_sl[old_start_level])*1000.0;
			for (n = old_start_level + 1; n <= start_level; n++)
			{
//				gwn += -ST[s].getDzui()[n] * (theta_sl_old[n] - theta_sl[n])*1000.0;
				gwn += -dicken[n] * (theta_sl_old[n] - theta_sl[n])*1000.0;
			}
		}
		else /* groundwater table rises into a higher level */
		{
			gwn  = 0;
			for (n = old_start_level; n > start_level; n--)
			{
//				gwn += ST[s].getDzui()[n] * (theta_sl[n] - theta_sl_old[n])*1000.0;
				gwn += dicken[n] * (theta_sl[n] - theta_sl_old[n])*1000.0;
			}
			gwn += (gw_exact-hgeo_sl[start_level])*(thetasat[start_level]-(theta_sl[start_level]-DT*lateralinflow[start_level]/dicken[start_level]))*1000.0;
		}
		/* return value: net groundwater recharge (without lateral inflows and outflows, which therefore have to be eleminated from the balance
	lateralinflow has to be subtracted from DELTAGW, 
	DELTAGW is negative and lateralinflow positive, so the sum means the subtraction */
		if (gwflow_active)
			return (double)(gwn+((double)DELTAGW[i][j]+(double)lateralinflow[0]+GWEXtotal-GWINtotal)*timestep); 
		else
			return (double)(gwn+((double)DELTAGW[i][j]+(double)lateralinflow[0])*timestep); 
	}
	else
	{
		// if the model is coupled with an external gw-model, then the calculation of gwn from the change in gw level will be 
		// unstable and not in all cases correct (because theta and theta_old and the other old_*-Values may be invalidated by a new external gw level
		// so we take the flow from the last completely unsaturated level as groundwater recharge (may be some early in case 
		// of short time steps and/or thick layers, but for a long modelling period it should be better than the internal methoid
		if (start_level>1) {
			gwn = qu_sl[start_level-1]*timestep;
		}
		else {
			gwn=0;
		}
		return -gwn; // negative because flows downwards are negativ and we want the recharge as positive value
	}

}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::infiltration(int i,int j,int &clay_layer,double &gw_depth,
int &s,double* &dicken,double* &hgeo_sl,double* &theta_sl,double &P,double &INFEX,
double timestep)
/* for the calculation of infiltration the precipitation is used in mm, not in m/s.
This is for compatibility reasons with the WaSiM-ETH model using TOPMODEL-algorithms  */
{
	double Na,temp3,temp4,infmax,A,B;
	double    Fsatt;        /* speichert die kumulative Infiltration waehrend des Saettigungsphase. */
	double    F;            /* gesamte Infiltration nach Saettigung. Fsatt und F addiert ergeben die gesamte
						Infiltrationsmene in einem Zeitintervall. Diese Summe wird ins **P-Grid geschrieben, vorher die
						Differenz **P-(F-Fsatt) ins **INFEX-Grid. */
	double ts;              /* Saettigungszeit, wird pro Pixel berechnet */
	double saug;            /* Saugspannung an der Feuchtefront (gut: 100mm) */
	double sicker=0;    /* fraction of infiltration excess, re-infiltrating into the soil  */

	gw_depth = -GWSTAND[i][j]; /* depth of ground water */
	if (CLAY)
	if (clay_layer>0) 
	{
		gw_depth=hgeo_sl[0]-hgeo_sl[clay_layer-1]; /* depth of clay layer is limiting infitration depth */
//		if (gw_depth <= ST[s].getDzui()[1])
//		{
//			gw_depth = ST[s].getDzui()[1] + 0.01; /* clay layer must be at least the second layer, the first layer is not allowed to be a clay layer  */
//		}
		if (gw_depth <= dicken[1])
		{
			gw_depth = dicken[1] + 0.01; /* clay layer must be at least the second layer, the first layer is not allowed to be a clay layer  */
		}
	}
	if (gw_depth>1.0) gw_depth=1.0; /* but not deeper than 1.0 m should infitration reach */
	//  Na  = (double)((double)1.0-(double)GetMeanTheta(0,gw_depth)); /* fillable fraction of pores */
	Na  = (double)((double)1.0-(double)GetMeanTheta(s,dicken,theta_sl,0.0,gw_depth)); /* fillable fraction of pores */

	Na  *= ST[s].getThetas(2); /* fillable fraction of entire soil */
	if (Na<0.0) Na=0.0;
	infmax    = 1000.0*Na; /* fillable porosity in [mm] */

	/* saturation can only occure, if P[i][j]> ST[s].k_sat[1] */
	/* calculation can only carried out, if ksat > 0 and
if Na > 0, otherwise division by zero is possible */
	if ((P>ST[s].getK_sat()[1]*timestep)&&(ST[s].getK_sat()[1]*timestep>0.0)&&(Na>0.0))
	{
		/* calculation of saturation time */
		saug=SOILtable[(int)SOIL[i][j]].getPsi();
		//    saug = 1000.0*((hgeo_sl[0]+hgeo_sl[1])/2-hh_sl[1]); /* suction in mm as calculated for the first layer */
		/* if there was already saturation at the time step begin, -> ts=0, no calculation of new saturation time */
		if (SATT[i][j]<0.5)
		ts=saug*Na/P/(P/(ST[s].getK_sat()[1]*timestep)-1);
		else
		ts=0.0;
		/* is saturation time within or outside interval end? */
		if (ts<=1)
		{
			/* saturation within interval. bacause rain intensity is constant, saturation remains still intervall end.
	SATT[i][j] has to set to 1. cumulative infiltration is calculated. */
			/* ts is a relative fracztion, Fsatt is estimated directly from ts*P */
			Fsatt=ts*P;

			temp3 = 1-ts;
			temp4 = Na*saug;

			A=ST[s].getK_sat()[1]*timestep*temp3;
			B=Fsatt+2*temp4;

			F=A/2+sqrt(A*A/4.0+A*B+Fsatt*Fsatt); /* nach Grundlagen der Hydrologie */

			/* calculation of F iterative */
			/*      iter  = 0;
err   = 1.0;
maxerr=P[i][j]>1?0.002:0.002*P[i][j];
Falt  = P[i][j]/2.0;
temp3 = 1-ts;
temp4 = Na*saug;
while ((iter<250)&&(err>maxerr))
{
	F=ST[s].k_sat[1]*timestep*temp3 + temp4*log((Falt+temp4)/(Fsatt+temp4))+Fsatt;
	err = fabs(F-Falt);
	Falt=F;
	iter++;
}
*/
			/* check if F > P */
			if (F>P) F=P;  else if (F<0) F=0;
			/* is there enough room in the soil for infiltration water? */
			if (F>infmax)
			{
				F=infmax;
			}
			/* balancing the input an output */
			INFEX = (P-F)*(1-sicker);
			P = F+(P-F)*sicker;;
			SATT[i][j] = 1.0;
		}
		else
		{
			/* if precipitation < ksat or if no saturation occured */
			INFEX = 0.0;
			SATT[i][j]  = 0.0;
		}
	}
	else
	{ /* rain intensity is less than minimum infitration rate (k_sat) */
		if (infmax<P)
		{ /* if soil will become saturated because of filling it up instead of to high rain intensities */
			INFEX = P-infmax;
			P = infmax;
			SATT[i][j] = 1.0;
		}
		/* no saturation -> all water is infitrating, no direct discharge */
		else
		{
			INFEX = 0.0;
			SATT[i][j]  = 0.0;
		}
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::flushbuffer(FILE *datei,char *dateiname,unsigned long *bufcnt,buffertype *buffer,int start,int maxlevels)
/* flushing the file buffer */
{
	int i;
	unsigned long j;

	if (*bufcnt>0)
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		(*WSstream) << "write "<<dateiname<<"[end-entry]";
		if ((datei=fopen(dateiname,"a+"))==NULL)
		{
			(*WSstream) << "FATAL: \nerror opening file "<<dateiname<<"[end-entry]";
			throw(65);
		}
		for (int j=1;j<=*bufcnt;j++)
		{
			fprintf(datei,"%ld\t%ld\t%ld\t%ld",buffer[j].j,buffer[j].m,buffer[j].t,buffer[j].s);
			if (maxlevels==0)
			{
				for (int i=start;i<=levels;i++)
					fprintf(datei,"\t%g",buffer[j].wert[i]);
			}
			else
			{
				for (int i=start;i<=maxlevels;i++)
					fprintf(datei,"\t%g",buffer[j].wert[i]);
			}
			fprintf(datei,"\n");
		}
		fclose(datei);
		(*WSstream) << "(" << (long)((*bufcnt)*((long)4*sizeof(long)+levels*sizeof(buffer[1].wert[1])))
		<<" of "<<(long)buffersize<<" buffer-bytes)"<<"[end-entry]";
		*bufcnt=0;
		SetSignalAbfangfunktion();
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::new_concentrations(int i,int j,
int &start_level,double** &conc_sl,double** &conc_sl_old,double* &qu_sl,double* &theta_sl,double* &hgeo_sl,
double* &lateralinflow,double* &gwin,double* &gwex,double* &theta_sl_old,double* &thetasat,double* &maxfromlr,
double* &dicken,double* &interflow,double* &qlat,double* &ratio,double* &iflfromlr,double* &qlatfromlr,
double* &gwexfromlr,double* &qufromlr,double* &conc_lat,double* &conc_gwex,double* &conc_gwin,
double &QDgrid,double &QIgrid,double GWEXtotal,double DT,double timestep)
{
	int z,lv;
	double masssum,volsum;
	double  massflux1, /* massflux in the inflow from the upper layer */
	massflux2, /* massflux in the outflow to the lower layer */
	massflux3, /* massflux in the interflow */
	massflux4, /* massflux in the lateral in- and outflows from groundwater */
	massflux5, /* massflux in the exfiltration to the river */
	massflux6, /* massflux in the infiltration from the river */
	trmsmax,         /* maximum allowed mass of a tracer to be transported with the calculated flus (e.g. due to limited salt concentration) */
	oldmass,newmass; /* masses of a tracer species in a layer before and after the flows */

	/* calculation of concentration changes within soil layers due to
mixing, radioactive decay and mass transfer */
	/* saving old concentrations and radioactive decay */
	for (int k = 0; k <= levels; k++)
	{
		for (z = 1; z <= tracer_count; z++)
		{
			/* radioactive decay */
			//	decay should be done after balancing because of grondwater models own decay algorithm
			//      if (tracer_info[z].radioactive)
			//        conc_sl[z][k]*=tracer_info[z].decfact;
			conc_sl_old[z][k] = conc_sl[z][k];
		}
	}

	if (levels>start_level) lv=start_level+1; else lv=levels;
	/* the last layer of unsaturated flow may be not the one, that was valid before a negative baseflow
could fill one ore more layers by lateralinflow. Therefor these layers have to be added to the number
of layers for which concentrations are updated now. */
	if (lv<levels)
	{
		while ((lateralinflow[lv+1]>0)&&(lv<levels))  lv++;
	}
	if (lv<levels)
	{
		while ((gwin[lv+1]>0||gwex[lv+1]>0)&&(lv<levels))	lv++;
	}
	if (lv>levels) lv=levels;

	for (int k=1;k<=lv;k++)
	{
		//    if (theta_sl_old[k]>thetares)
		if (theta_sl_old[k]>0)
		//      maxfromlr[k]=(theta_sl_old[k]-thetares)*dicken[k]/DT; /* what part of flow can at maximum come from the layer itself */
		maxfromlr[k]=theta_sl_old[k]*dicken[k]/DT; /* what part of flow can at maximum come from the layer itself */
		else
		maxfromlr[k]=0;
		/* if more flow is calculated, then it will usually come from the layer(s) above (it means that there is throughflow from the layer(s) above)*/
		if (((interflow[k]+qlat[k]-qu_sl[k]+gwex[k])>maxfromlr[k])&&(maxfromlr[k]>0.0))     
		ratio[k] = maxfromlr[k]/(interflow[k]+qlat[k]-qu_sl[k]+gwex[k]);
		else
		ratio[k] = 1.0;

		iflfromlr[k]  = interflow[k]*ratio[k];
		qlatfromlr[k] = qlat[k]*ratio[k];
		gwexfromlr[k] = gwex[k]*ratio[k];
		qufromlr[k]   = qu_sl[k]*ratio[k];
	}

	for (z=1;z<=tracer_count;z++)
	{
		//    conc_ifl[z]  = 0.0;
		CONC_QI[z][i][j]  = 0.0;
		conc_lat[z]  = 0.0;
		conc_gwex[z] = 0.0;
		for (int k=1;k<=lv;k++)
		{
			if (tracer_info[z].evaporating)
			{ /* mixing concentrations */
				/* plugflow is assumed, if more water flows out than there was in the layer before */
				oldmass = conc_sl_old[z][k]*(1000.0*theta_sl_old[k]*dicken[k]); /* old mass of tracer (resp. relative parts) */
				/* 1. flow between upper and actual layer */
				if (qu_sl[k-1]<0) /* infolow from above */
				{
					if (k>1) massflux1=massflux2;
					else
					{ 
						if (qu_sl[0]<0)
						massflux1 = conc_sl_old[z][0] *(-qu_sl[0])*timestep; 
						else massflux1 = 0.0;
					}                                                    
				}
				else /* outflow upwards */
				{
					if (k>1) massflux1=massflux2;
					else  massflux1 = conc_sl_old[z][k]*(-qu_sl[k-1])*timestep;
				}
				
				/* 2. flow between actual and lower layer */
				if (qu_sl[k]>0) /* inflow from below */
				{
					if (k<levels) massflux2 = -conc_sl_old[z][k+1]*qu_sl[k]*timestep; /* attention! indizes: lower boundary flow comes from [k+1]but is flux [k] in the actual compartiment! */
					else massflux2=0.0;
				}
				else /* mass loss downwards */
				{
					if (k<levels)
					{
						if (ratio[k]>0.99999)
						massflux2 = -conc_sl_old[z][k]*qu_sl[k]*timestep;
						else
						{
							/* first part: old tracer content of actual layer is transported downwards (assuming, water flushes all pores) */
							massflux2 = -conc_sl_old[z][k]*qufromlr[k]*timestep; 
							/* massflux2 is not exactly the old tracer content of the layer, bacause there may have been also interflow or baseflow! */
							/* second part: other flushing water from above */
							massflux2 += -conc_sl[z][k-1]*(qu_sl[k]-qufromlr[k])*timestep;
						}
					}
				}
				
				/* considering interflow */
				if (ratio[k]>0.99999)
				massflux3 = conc_sl[z][k]*interflow[k]*timestep; 
				else
				{ 
					massflux3  = conc_sl[z][k]*iflfromlr[k]*timestep;
					massflux3 += conc_sl[z][k-1]*(interflow[k]-iflfromlr[k])*timestep;
				}
				CONC_QI[z][i][j]+=massflux3/timestep;
				//conc_ifl[z]+=massflux3/timestep;
				
				/* considering lateral inflows */
				if (ratio[k]>0.99999) /* qlat < 0 is not possible, because if DELTAGW is < 0, qlat is not calculated (qlat==0), but lateralinflows are calculated*/
				massflux4= conc_sl[z][k]*qlat[k]*timestep; 
				else
				{ 
					massflux4  = conc_sl[z][k]*qlatfromlr[k]*timestep;
					massflux4 += conc_sl[z][k-1]*(qlat[k]-qlatfromlr[k])*timestep;
				}
				/* considering lateral inflows calculated according to the fluxes from the groundwater module */
				if (gwflow_active)
				massflux4 -= gw_conc[z][i][j]*lateralinflow[k]*timestep; /* because sign of massflux 4 is negative in the balance inflowing mass has to be subtracted from maxxflux 4 */
				conc_lat[z]+= massflux4/timestep;
				
				/* groundwater exfiltration */
				if (ratio[k]>0.99999)
				massflux5 = conc_sl[z][k]*gwex[k]*timestep; 
				else
				{ 
					massflux5  = conc_sl[z][k]*gwexfromlr[k]*timestep;
					massflux5 += conc_sl[z][k-1]*(gwex[k]-gwexfromlr[k])*timestep;
				}
				conc_gwex[z]+=massflux5/timestep;

				/* infiltration from rivers */
				if (gwin[k]>0)
				/* concentration in the river was already calculated in the flowrouting 
		and stored in the riverflow_table[].flowmmh[1..9] for each tracer. During 
		calculation of infiltration from the river, this value was stored in conc_gwin[k] 
		(for each affected layer the same concentration) */
				massflux6 = conc_gwin[z]*gwin[k]*timestep; 
				else
				massflux6=0.0;
				
				newmass = oldmass + massflux1 - massflux2 - massflux3 - massflux4 - massflux5 + massflux6;
				conc_sl[z][k] = theta_sl[k]>0? newmass/(1000.0*dicken[k]*theta_sl[k]) : 0.0;
			}
			else
			{ /* adding masses : massflux 1 to 6 are masses! */
				if (qu_sl[k-1]<0) /* inflow from above */
				{
					if (k>1)
					{
						massflux1=massflux2; /* starting with 2nd layer: the flow between actual and upper layer is the 
								same as in the last layer the flow between the actual and the lower layer */
						// or the same: massflux1 = theta_sl_old[k-1]>0?conc_sl_old[z][k-1]*(-qu_sl[k-1]*DT/dicken[k-1])/theta_sl_old[k-1]:0;
						//            trmsmax   = -tracer_info[z].max_conc*qu_sl[k-1]*timestep; /* interflow*DT: flux in m; *1000: flux in liter or kg/m^2 */
						//            if (trmsmax<massflux1) massflux1=trmsmax;
					} /* a fraction of the upper mass is won */
					else
					{
						/* it is asuemd, that QDgrid-qu_sl[0] will never be zero, because QDgrid is always >= 0 
		and if it's 0, then this statement is only reached, if qu_sl[0] is larger than 0 */
						if ((qu_sl[0]<0)&&((QDgrid-qu_sl[0])>0.0))
						massflux1 = conc_sl_old[z][0] *(-qu_sl[0]/(QDgrid-qu_sl[0])); /* QDgrid-qu_sl[0] is the sum of surface flow and infiltrating water */
						else massflux1 = 0.0;
						trmsmax  = -tracer_info[z].max_conc*qu_sl[0]*timestep; /* interflow*DT: flux in m; *1000: flux in liter or kg/m^2 */
						if (trmsmax<massflux1) massflux1=trmsmax;
					} /* if this is the first layer, the upper layer is rain or irrigation water */
				}
				else /* outflow -> mass loss upwards */
				{
					if (k>1)
					massflux1 = massflux2;
					else
					{
						massflux1 = theta_sl_old[k]>0.001 ? conc_sl_old[z][k]*(qu_sl[k-1]*DT/dicken[k])/theta_sl_old[k]:0.0; /* fraction of the actual mass is lossed -> negative sign */
						trmsmax  = tracer_info[z].max_conc*qu_sl[k-1]*timestep;
						if (trmsmax<massflux1) massflux1=trmsmax;
						massflux1*=(double)-1.0;
					}
				}

				if (qu_sl[k]>0)  /* inflow from below */
				{
					if (k<levels)
					{
						massflux2 = theta_sl_old[k+1]>0.001 ? conc_sl_old[z][k+1]*(-qu_sl[k]*DT/dicken[k+1])/theta_sl_old[k+1] : 0; /* attention! indizes: lower boundary flow comes from [k+1]but is flux [k] in the actual compartiment! */
						trmsmax  = tracer_info[z].max_conc * -qu_sl[k]*timestep;
						if (-trmsmax < -massflux2) massflux2=trmsmax; /* negative signs, because inflow from below (but: inflow from above is positive, outflow downwards is also positive) */
					}
					else massflux2=0.0;
				}
				else /* outflow downwards -> loss */
				{
					if (k<levels)
					{
						if (ratio[k]>0.99999)
						/* if no troughflow trough a layer takes place */
						massflux2 = theta_sl_old[k]>0.001? -conc_sl_old[z][k]*(qu_sl[k]*DT/dicken[k])/theta_sl_old[k]:0.0; /* qu_sl[k] is < 0 -> loss */
						else
						{ /* more water flows out, than was in the layer, so other water with other concentrations enters the layer */
							/* first part: old tracer content of actual layer is transported downwards (assuming, water flushes ALL pores) */
							massflux2 = theta_sl_old[k]>0.001?-conc_sl_old[z][k]*(qufromlr[k]*DT/dicken[k])/theta_sl_old[k]:0.0; 
							/* it is not the old tracer content of the layer, bacause there may have been also interflow or baseflow! */
							/* second part: other flushing water from above */

							if (k>1)
							{
								if (theta_sl[k-1]>0)
								massflux2 += -conc_sl[z][k-1] * ((qu_sl[k]-qufromlr[k])*DT)/(1000.0*theta_sl[k-1]*dicken[k-1]);
							}
							else
							{
								if ((QDgrid-qu_sl[0])>0.0)
								massflux2 += -conc_sl[z][k-1] * ((qu_sl[k]-qufromlr[k])*DT)/((QDgrid-qu_sl[0])*DT); /* QDgrid-qu_sl[0] is precipitation+irrigation */
							}

							//              massflux2 += conc_sl[z][k-1]/(k>1?1000.0*theta_sl[k-1]*dicken[k-1]:(QDgrid-qu_sl[0])*DT)*(qufromlr[k]-qu_sl[k])*DT;
							/* the above line means: the tracer mass, that was within the actual layer, is flushed out, and now to this mass
			the mass of tracer contained in the additional flushing water is added. Because for non-evaporating tracers the masses instead of
			the centrations are stored, but for calculating inflowing mass the flux and the cocentrations are needed, 
			first the concentration for the water from the upper layer is calculated and then this concentration is 
			multiplied by the flow from this layer going trough the actual layer. If the actual layer is the first layer, the incoming
			flux has the concentration of rain resp. irigation water, which may be devided into direct flow (QDgrid) and infiltration qu_sl[0].
			therefor the concentration of throughflowing water is calculated from tracer mass at the surface divided by (QD-qu_sl[0]), wherby qu_sl[0] is negative 
		*/
						}
						trmsmax  = tracer_info[z].max_conc * -qu_sl[k] * timestep;
						if (trmsmax<massflux2) massflux2=trmsmax;
					}
				}

				/* considering interflow */
				if (ratio[k]>(double)0.99999)
				massflux3 = theta_sl_old[k]>(double)0.001 ? conc_sl[z][k]*(interflow[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
				else 
				{ massflux3  = theta_sl_old[k]>(double)0.001 ? conc_sl[z][k]*(iflfromlr[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
					/* comment for following row see above at comments for massflux2 */
					if (theta_sl_old[k-1]>(double)0.001)
					if (k>1)
					{
						if (theta_sl[k-1]>0)
						massflux3 += (conc_sl[z][k-1] * ((interflow[k]-iflfromlr[k])*DT)/(1000.0*theta_sl[k-1]*dicken[k-1]));
					}
					else
					{
						if ((QDgrid-qu_sl[0])>0.0)
						massflux3 += (conc_sl[z][k-1] * ((interflow[k]-iflfromlr[k])*DT)/((QDgrid-qu_sl[0])*DT));
					}
					//        massflux3 += theta_sl_old[k-1]>(double)0.0001 ? conc_sl[z][k-1]/(k>1?1000.0*theta_sl[k-1]*dicken[k-1]:(QDgrid-qu_sl[0])*DT)*(interflow[k]-iflfromlr[k])*DT : 0.0;
				}
				trmsmax   = tracer_info[z].max_conc*interflow[k]*timestep; /* interflow*DT: flux in m; *1000: flux in liter or kg/m^2 */
				if (trmsmax<massflux3) massflux3=trmsmax;
				//conc_ifl[z]+=massflux3;
				CONC_QI[z][i][j]+=massflux3;

				/* considering baseflow */
				if (ratio[k]>0.99999)
				massflux4 = theta_sl_old[k]>0.001 ? conc_sl[z][k]*(qlat[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
				else
				{ massflux4  = theta_sl_old[k]>0.001 ? conc_sl[z][k]*(qlatfromlr[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
					/* comment for following row see above at comments for massflux2 */
					if (theta_sl_old[k-1]>(double)0.001) 
					{
						if (k>1)
						{
							if (theta_sl[k-1]>0)
							massflux4 += (conc_sl[z][k-1] * ((qlat[k]-qlatfromlr[k])*DT)/ (1000.0*theta_sl[k-1]*dicken[k-1]));
						}
						else
						{
							if ((QDgrid-qu_sl[0])>0)
							massflux4 += (conc_sl[z][k-1] * ((qlat[k]-qlatfromlr[k])*DT)/  ((QDgrid-qu_sl[0])*DT));
						}
					}
					//        massflux4 += theta_sl_old[k-1]>0.0001 ? conc_sl[z][k-1]/(k>1?1000.0*theta_sl[k-1]*dicken[k-1]:(QDgrid-qu_sl[0])*DT)*(qlat[k]-qlatfromlr[k])*DT : 0.0;
				}
				trmsmax   = tracer_info[z].max_conc*qlat[k]*timestep; /* interflow*DT: flux in m; *1000: flux in liter or kg/m^2 */
				if (trmsmax<massflux4) massflux4=trmsmax;
				/* considering lateral inflows calculated according to the fluxes from the groundwater module -> also salts are handled by their concentrations */
				if (gwflow_active)
				massflux4 -= gw_conc[z][i][j]*lateralinflow[k]*timestep; /* because sign of massflux 4 is negative in the balance inflowing mass has to be subtracted from maxxflux 4 */
				conc_lat[z]+= massflux4;

				/* groundwater exfiltration into rivers */
				if (ratio[k]>(double)0.99999)
				massflux5 = theta_sl_old[k]>(double)0.001 ? conc_sl[z][k]*(gwex[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
				else 
				{ massflux5  = theta_sl_old[k]>(double)0.001 ? conc_sl[z][k]*(gwexfromlr[k]*DT/dicken[k])/theta_sl_old[k] : 0.0;
					/* comment for following row see above at comments for massflux2 */
					if (theta_sl_old[k-1]>(double)0.001)
					{
						if (k>1)
						{
							if (theta_sl[k-1]>0)
							massflux5 += (conc_sl[z][k-1] * ((gwex[k]-gwexfromlr[k])*DT)/(1000.0*theta_sl[k-1]*dicken[k-1]));
						}
						else
						{
							if ((QDgrid-qu_sl[0]))
							massflux5 += (conc_sl[z][k-1] * ((gwex[k]-gwexfromlr[k])*DT)/((QDgrid-qu_sl[0])*DT));
						}
					}
					//        massflux5 += theta_sl_old[k-1]>(double)0.0001 ? conc_sl[z][k-1]/(k>1?1000.0*theta_sl[k-1]*dicken[k-1]:(QDgrid-qu_sl[0])*DT)*(interflow[k]-iflfromlr[k])*DT : 0.0;
				}
				trmsmax   = tracer_info[z].max_conc*gwex[k]*timestep; /* gwex*timestep: flux in mm == kg/m^2; *conc: massflux in kg/m^2 */
				if (trmsmax<massflux5) massflux5=trmsmax;
				conc_gwex[z]+=massflux5;

				/* infiltration from the river -> concentration in the river is never above saturation, so trmsmax is not calculated here */
				if (gwin[k]>0.0)
				massflux6 = conc_gwin[z]*gwin[k]*timestep; /* conc_gwin[k] is given as concentration, not as mass, because in the river only concentratiuos are modeled! */
				else
				massflux6=0;

				conc_sl[z][k] = conc_sl_old[z][k] + massflux1 - massflux2 - massflux3 - massflux4 - massflux5 + massflux6;
				if (k==1) /* balance with surface storage */
				conc_sl[z][0]-=massflux1;
			}
			fpe_check(i,j,k,z,massflux1,massflux2,massflux3,massflux4,massflux5);
		}

		/* evaporating tracers: averaging of concentrations. non-evaporating tracers: only adding masses - this has already been done */
		if (tracer_info[z].evaporating)
		{
			//      if (QIgrid>0.0) conc_ifl[z]/=QIgrid; else conc_ifl[z]=0.0;
			if (QIgrid>0.0) CONC_QI[z][i][j]/=QIgrid; else CONC_QI[z][i][j]=0.0;

			if (DELTAGW[i][j]>0.0) conc_lat[z]/=DELTAGW[i][j]; else conc_lat[z]=0.0;
			if (gwflow_active)
			if (GWEXtotal>0.0) conc_gwex[z]/=GWEXtotal; else conc_gwex[z]=0.0;
			else
			conc_gwex[z]=0.0;
		}

		/* mixing ocncentrations in all groundwater levels at once with the inflow into groundwater (recharge) */
		masssum=0.0;volsum=0.0;
		if (tracer_info[z].evaporating)
		{
			for (int k=start_level;k<=levels;k++)
			{                                         
				masssum += conc_sl[z][k]*theta_sl[k]*dicken[k];
				volsum  += theta_sl[k]*dicken[k];
			}
			if (gwflow_active) 
			{
				masssum += gw_conc[z][i][j]*thetasat[levels]*(gw_thick[i][j]-(hgeo_sl[0]-hgeo_sl[levels]));
				volsum  += thetasat[levels]*(gw_thick[i][j]-(hgeo_sl[0]-hgeo_sl[levels]));
				//        masssum += gw_conc[z][i][j]*thetasat[start_level]*(gw_thick[i][j]-(levels-1)*dicken[levels]);
				//        volsum  += thetasat[start_level]*(gw_thick[i][j]-(levels-1)*dicken[levels]);
			}
			conc_sl[z][start_level]=volsum>0?masssum/volsum:0; /* mixing over all groundwater layers */
			for (int k=start_level+1;k<=levels;k++) conc_sl[z][k]=conc_sl[z][k-1];
			if (gwflow_active)
			gw_conc[z][i][j]=conc_sl[z][levels];
		}
		else
		{
			for (int k=start_level;k<=levels;k++)
			{
				masssum+=conc_sl[z][k];
				volsum+=theta_sl[k]*dicken[k];
			}
			if (gwflow_active)
			{
				masssum += gw_conc[z][i][j]*thetasat[start_level]*1000.0*(gw_thick[i][j]-(hgeo_sl[0]-hgeo_sl[levels])); /* in groundwater non-evaporating tracers are not longer handled as masses */
				volsum  += thetasat[start_level]*(gw_thick[i][j]-(hgeo_sl[0]-hgeo_sl[levels]));
				//        masssum += gw_conc[z][i][j]*thetasat[start_level]*1000.0*(gw_thick[i][j]-(levels-1)*dicken[levels]); /* in groundwater non-evaporating tracers are not longer handled as masses */
				//        volsum  += thetasat[start_level]*(gw_thick[i][j]-(levels-1)*dicken[levels]);
			}
			for (int k=start_level;k<=levels;k++)
			conc_sl[z][k]= volsum > 0 ? masssum*(dicken[k]*theta_sl[k])/volsum : 0.0;
			if (gwflow_active)
			gw_conc[z][i][j]=volsum>0?masssum/volsum/1000.0:0;
		}
		fpe_check(i,j,levels,z,massflux1,massflux2,massflux3,massflux4,massflux5);
	} /* next tracer */

	for (int k=0;k<=levels;k++)
	for (z=1;z<=tracer_count;z++)
	{
		/* radioactive decay */
		if (tracer_info[z].radioactive)
		conc_sl[z][k]*=tracer_info[z].decfact;
	}

}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::get_max_pumping()
{
	int i,j,s,k;
	double th;
	if (do_it==0) return;
	else
	{
		if (irrigationtable)
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						//          if (zn[i][j]!=nodata)
					{
						s = (int)SOIL[i][j];
						int z;
						MAXPUMP[i][j] = 0;
						for (z = (int)GWLEVEL[i][j]; z < ST[s].getLevels(); z++)
						{
							th = ST[s].get_theta(z, 3.45); /* a suction of 3.45 m is assumed to coincide with field capacity */
							MAXPUMP[i][j] += 1000.0*ST[s].getDzui()[z] * (ST[s].getThetas(z) - th);
//							MAXPUMP[i][j] += 1000.0*dicken[z] * (ST[s].getThetas(z) - th);
						}
					}
					else
					{
						MAXPUMP[i][j] = nodata;
					}
				}
			}
		}
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
// calculate macropore infiltration depending on soil properties and rain intensity
// the resulting amount is stored in a grid (MACRO) as well as in a statistics file
// the input for Richards is actualised (theta-values, rain a.s.o.)
double unsatzonmodellklasse::CalculateMacroporeInfiltration(double &p,int &s,double* &dicken,
double* &macrocapacity,double* &thetasat,double* &theta_sl,
double timestep,double* &macroinfiltration)
{
	int k;
	double maxCapacity;
	double maxAvailable;
	double dicksum=0;
	double infiltrationSum=0;
	double infiltAmount=0;
	double MacroDepth;   // maximum depth of the macropores                                                                                                                                   
	double PMacroThresh; // precipitation capacity thresholding macropore runoff in mm per hour (not in m/s, because it's more convenient than to write it down in m/s, e.g. 5mm/h = 1.38e-6) 
	double CapacityRedu; // reduction of the macropore capacity with depth -> pores become less dense. This Factor describes the reduction ratio per meter                                    

	PMacroThresh = ST[s].get_PMacroThreshold();
	
	if (p>PMacroThresh)
	{
		MacroDepth   = ST[s].get_MacroDepth();
		CapacityRedu = ST[s].get_ReductionPerMeter();
		
		// macrocapacity auf Schichten berechnen
		maxCapacity =	ST[s].get_PMacroCapacity();
		k=0;
		while ((dicksum<(MacroDepth-dicken[k+1]/2.0)) && (k<levels))
		{
			k++;
			dicksum+=dicken[k];
			macrocapacity[k]=maxCapacity*pow(CapacityRedu,dicksum-dicken[k]/2.0);
		}

		// nun von unten nach oben Auffüllen des Bodens mit Makroporeninfiltration...
		// maximum ist entweder die max. Kapazität oder, wenn weniger Regen fiel, die Differenz aus Regen und Threshold
		// in m/s
		maxAvailable = ((p-PMacroThresh)>maxCapacity)?maxCapacity:(p-PMacroThresh);
		// füllen von unten her -> k wird wieder runtergezählt.
		while ((k>0) && (maxAvailable>0))
		{
			// in eine Schicht passt entweder soviel Wasser, wie in macrocapacity[k] berechnet wurde, oder weniger (wenn sie schon zu feucht ist)
			double maxFill;
			infiltAmount=0; // initializing
			maxFill = (thetasat[k]-theta_sl[k])*dicken[k]*1000.0; // in mm Wassersäule umgerechneter auffüllbarer Porenraum
			if (maxFill>0)
			{
				if ((maxFill < maxAvailable*timestep)&&(maxFill < macrocapacity[k]*timestep))
				{ // layer will be filled to saturation, but both, the absolute amount in maxAvailable (which is updated after every layer) and the capacity for the actual layer itself must be meet
					infiltAmount=maxFill;
					theta_sl[k]=thetasat[k];
					maxAvailable-=(maxFill/timestep);
					macroinfiltration[k]=infiltAmount/timestep;
				}
				else
				{ 
					if (maxAvailable > macrocapacity[k])
					{	// macrocapacity[k] is limiting the macro pore infiltration, 
						// there will remain some water in maxAvailable 
						// and the layer is not filled completely, because maxFill is larger than macrocapacity[k] 
						// (otherwise it would have been handled above in the if-block of the enclosing if-else-statement)
						infiltAmount=macrocapacity[k]*timestep;
						theta_sl[k]+=infiltAmount/(1000.0*dicken[k]);
						maxAvailable-=macrocapacity[k];
						macroinfiltration[k]=infiltAmount/timestep;
					}
					else
					{	// not enough water left to even satisfy macrocapacity[k], so thats the end of the infiltration for this cell and this timestep
						// maxAvailable will be set to 0 and only this layer is filled up partly, all layers above get no macropore water at all
						infiltAmount=maxAvailable*timestep;
						theta_sl[k]+=infiltAmount/(1000.0*dicken[k]);
						maxAvailable=0.0;
						macroinfiltration[k]=infiltAmount/timestep;
					}
				}
			}
			infiltrationSum+=infiltAmount;
			k--;
		}
		p-=(double)(infiltrationSum/timestep);
		
	}
	return infiltrationSum;
}
//#pragma optimize ("p",off)


// expanding the change in inner energy per temperature change with latent heat, densitiy, and the effective water content Theta-Theta_res
void unsatzonmodellklasse::Expand_dEoverdT(int offset,double *theta_sl,double *theta_res,double *dEOverdT_sl)
{
	for (int i=1+offset;i<=levels;i++)
	{
		if (theta_sl[i]-theta_res[i] > 0)
		{
			dEOverdT_sl[i] *= (this->dTModLatHeat*this->dTModRhoWater*(theta_sl[i]-theta_res[i]));
		}
		else
		{
			dEOverdT_sl[i] = 0;
		}
	}
}

// double *deltaT_sl  result: temperature change
// double *qHeat_sl	 input: heat fluxes between levels
// double *theta_sl   input: theta-values (water content per level)
// double *heatcap_sl input: effective actual heat capacity per level depending on temperatur, absolute water content and relative saturation (water/ice)
// double *dicken	 input: level thickness
// int s			 input: soil type (for getting dry soil density)
// int timestep      input: (sub)-time step in seconds
void unsatzonmodellklasse::CalcDeltaT(int offset,double *deltaT_sl, double *qHeat_sl,double *theta_sl,double *dicken,double *heatcap_sl,int s,int timestep)
{
	for (int i=1+offset;i<=levels;i++)
	{
//		deltaT_sl[i] = (timestep / (this->dTModRhoWater*dicken[i]*heatcap_sl[i])) * (qHeat_sl[i]-qHeat_sl[i-1]);
		deltaT_sl[i] = (timestep / (dicken[i]*heatcap_sl[i])) * (qHeat_sl[i]-qHeat_sl[i-1]);
	}
}

// double *temper_sl_new		Input/Output: the temperature deltas in deltaT_sl will be added to the input values and returned as new temperatures
// double *temper_sl			Input: temperatures when the sub time step startet (not changed by iterations): used for comparison of old and new temperatures to limit the changes
// double *deltaT_sl			Input: the changes calculated from the heat fluxes
// double *qHeat_sl			Input: the heat fluxes itself, used to get the sign of change
void unsatzonmodellklasse::CalcNewTemperature(int offset,double *temper_sl_new,double *temper_sl,double *deltaT_sl,double *qHeat_sl,double Tbndlow,double Tbndup,bool withlimitation,double Tsurf)
{
	// starting from the buttom since the flux restrictions are better to handle this way (smaller amplitudes in deeper layers)
	for (int i=levels;i>offset;i--)
	{
		temper_sl_new[i]+=deltaT_sl[i];
		// limiting heat fluxes in order to prevent numerical instabilities
		if ((withlimitation)&&(i>offset)) // but only for real layers, not the 0-layer
		{
			double tlow;
			double tup;  
			if (i==levels) 
			{ // lowest layer: replace T(i+1) with lower boundary condition
				tlow = Tbndlow;
				tup = temper_sl[i-1]; // --> this is the old temperature from the last (sub) time step
			} 
			else 
			{
				if (i==1)
				{ // uppermost layer
					tlow = temper_sl[i+1]; // --> this is the old temperature from the last (sub) time step
					tup = (Tsurf>-9999)?Tsurf:Tbndup;
				}
				else
				{ // all normal layers
					tlow = temper_sl[i+1]; // --> this is the old temperature from the last (sub) time step
					tup = temper_sl[i-1]; // --> this is the old temperature from the last (sub) time step
				}
			}
			// case 1: both fluxes positive (upwards)
			if ((qHeat_sl[i-1]>=0)&&(qHeat_sl[i]>=0)) {
				if (temper_sl_new[i] > tlow)  
				{
					 temper_sl_new[i] = tlow;
				}
			} else if ((qHeat_sl[i-1]>=0)&&(qHeat_sl[i]<0)) {
				double tmin=tup;
				if (tmin>tlow) tmin=tlow;
				if (temper_sl_new[i] < tmin)
				{
					temper_sl_new[i] = tmin;
				}
			} else if ((qHeat_sl[i-1]<0)&&(qHeat_sl[i]<0)) {
				if (temper_sl_new[i] > tup)
				{
					temper_sl_new[i] = tup;
				}
			} else if ((qHeat_sl[i-1]<0)&&(qHeat_sl[i]>=0)) {
				double tmax=tup;
				if (tmax<tlow) tmax=tlow;
				if (temper_sl_new[i] > tmax)
				{
					temper_sl_new[i] = tmax;
				}
			}
		}
	}
}

void unsatzonmodellklasse::CalcHeatFluxes(int offset,double *qHeat_sl,double *lambda_sl,double *dicken,double *temper_sl,double TLOWBND,double Tair,double actNFactor,double Tsurf)
{
	// flux from upper boundary into 1st layer
	if (Tsurf>-9999)
	{
		qHeat_sl[offset] = (lambda_sl[offset+1]*-2.0)/(dicken[offset+1]) * (Tsurf - temper_sl[offset+1]);
	}
	else
	{
		qHeat_sl[offset] = (lambda_sl[offset+1]*-2.0)/(dicken[offset+1]) * (Tair*actNFactor - temper_sl[offset+1]);
	}
	// fluxes between regular layers (only between actual layer and next lower layer; flux is positive upwards, negative downwards
	for (int i=offset+1;i<levels;i++)
	{
		double deff = -0.5*(dicken[i]+dicken[i+1]);
		qHeat_sl[i] = double((lambda_sl[i]/deff) * (temper_sl[i]-temper_sl[i+1]));
	}
	// flux between lowermost level and lower boundary flux is positive upwards, negative downwards
	qHeat_sl[levels] = (lambda_sl[levels]*-2.0)/(dicken[levels]) * (temper_sl[levels] - TLOWBND);
}

void unsatzonmodellklasse::SolveTriDiag(int offset,int row,int col,double *lambda_sl,double *heatcap_sl,double *dicken,double *temper_sl,double timestep_heat,
		double *coeff_a,double *coeff_b,double *coeff_c,double *right_v,double *result_x,double actNFactor,
		double *macroinfiltration,double *qu_sl,double Tsurf)
{
	// defining the coefficients for the tridiagonal system of linear equations which can be used to solve the heat transfer equation for finite differences when using an implicit approach 
	// (where the resulting temperatures depend on the temperatures at time n+1 instead on temperatures at time n only as is the case for the explicit time stepping used before)

	double temp_a=0;
	double temp_c=0;
	double cw=4.187e6; // heat capacity of water in J/(m^3*K)
	double x=0;
	int i = 1+offset;
	// uppermost node: here we have a known boundary condition T[i][j] (air temperature) which is used to move the first term to the known reight side of the equation.
	if (iTModMethod>11)
	{ // including 1D vertical advection and inflow from macropore infiltration
		temp_a = (2.0*timestep_heat*lambda_sl[i])/(dicken[i]*dicken[i]*heatcap_sl[i]);
		if (qu_sl[0+offset]<0)
		{
			temp_a += (qu_sl[0+offset]*cw*timestep_heat)/(dicken[1]*heatcap_sl[1]);
		}
		coeff_a[i]  = 0;
		coeff_c[i]  = -(2.0*timestep_heat*lambda_sl[i])/(dicken[i]*(dicken[i]+dicken[i+1])*heatcap_sl[i]);
		if (qu_sl[i+1]>0)
		{
			coeff_c[i] += (qu_sl[i+1]*cw*timestep_heat)/(dicken[i]*heatcap_sl[i]);
		}
//				coeff_b[i]  = 1 - coeff_c[i] + temp_a;
		coeff_b[i]  = 1 - coeff_c[i] + temp_a - (macroinfiltration[i]*cw*timestep_heat)/(dicken[i]*heatcap_sl[i]);
//				double Tact = (SNOW && SNOW[row][col]>0) ? 0 : T[row][col];
		double Tact = T[row][col];
//				double Tact = (T[row][col]>0)?T[row][col]:0;
		if (SNOW)
		{
			if (SNOW[row][col] > 0) Tact = 0;
		}
//				right_v[i] = temper_sl[i] + T[row][col]*actNFactor*temp_a;
		right_v[i] = temper_sl[i] - (timestep_heat/(dicken[i]*heatcap_sl[i]))*macroinfiltration[i]*cw*Tact + Tsurf*actNFactor*temp_a;
//				right_v[i] = temper_sl[i] - (timestep_heat/(dicken[i]*heatcap_sl[i]))*macroinfiltration[i]*cw*Tact + T[row][col]*temp_a; 
//				right_v[i] = temper_sl[i] - (timestep_heat/(dicken[i]*heatcap_sl[i]))*macroinfiltration[i]*cw*Tact + Tact*temp_a; 
		result_x[i] = 0;
	}
	else
	{   // only heat diffusionion
		temp_a = (2.0*timestep_heat*lambda_sl[i])/(dicken[i]*dicken[i]*heatcap_sl[i]);
		coeff_a[i]  = 0;
		coeff_c[i]  = -(2.0*timestep_heat*lambda_sl[i])/(dicken[i]*(dicken[i]+dicken[i+1])*heatcap_sl[i]);
		coeff_b[i]  = 1 - coeff_c[i] + temp_a;
//				right_v[i]  = temper_sl[i] + T[row][col]*actNFactor*temp_a; // with air temperature, melt may be too low during the melting season
		right_v[i]  = temper_sl[i] + Tsurf*actNFactor*temp_a; // surface temperature is much higgher than 0°C for melting conditions since it mus treflect the energy balance as a heat flux only (instead of radiation that will reach deepr into the snow)
		result_x[i] = 0;
	}

	// now all normal layers
	if (iTModMethod>11)
	{ // including 1D vertical advection and inflow from macropore infiltration
		for (int i = 2 + offset; i < levels; i++)
		{
			coeff_a[i] = -(2.0*timestep_heat*lambda_sl[i - 1]) / (dicken[i] * (dicken[i] + dicken[i - 1])*heatcap_sl[i]);
			if (qu_sl[i - 1] < 0)
			{
				coeff_a[i] -= (qu_sl[i - 1] * cw*timestep_heat) / (dicken[i] * heatcap_sl[i]);
			}
			coeff_c[i] = -(2.0*timestep_heat*lambda_sl[i]) / (dicken[i] * (dicken[i] + dicken[i + 1])*heatcap_sl[i]);
			if (qu_sl[i + 1] > 0)
			{
				coeff_c[i] += (qu_sl[i + 1] * cw*timestep_heat) / (dicken[i] * heatcap_sl[i]);
			}
			coeff_b[i] = 1 - coeff_a[i] - coeff_c[i] - (macroinfiltration[i] * cw*timestep_heat) / (dicken[i] * heatcap_sl[i]);
			//				double Tact = (SNOW && SNOW[row][col]>0) ? 0 : T[row][col];
			double Tact = T[row][col];
			//				double Tact = (T[row][col]>0)?T[row][col]:0;
			if (SNOW)
			{
				if (SNOW[row][col] > 0) Tact = 0;
			}

			right_v[i] = temper_sl[i] - (timestep_heat / (dicken[i] * heatcap_sl[i]))*macroinfiltration[i] * cw*Tact;
			result_x[i] = 0;
		}
	}
	else
	{  // only heat diffusionion
		for (int i = 2 + offset; i < levels; i++)
		{
			coeff_a[i] = -(2.0*timestep_heat*lambda_sl[i - 1]) / (dicken[i] * (dicken[i] + dicken[i - 1])*heatcap_sl[i]);
			coeff_c[i] = -(2.0*timestep_heat*lambda_sl[i]) / (dicken[i] * (dicken[i] + dicken[i + 1])*heatcap_sl[i]);
			coeff_b[i] = 1 - coeff_a[i] - coeff_c[i];
			right_v[i] = temper_sl[i];
			result_x[i] = 0;
		}
	}
	i = levels;

	// and now the lowermost layer, having a lower boundary condition
	double Tlow=TLOWBND[row][col];
	if (useHeatFluxAsLowerBC==1)
	{
		Tlow = temper_sl[i]+dicken[i]*QLOWBND[row][col]/(2*lambda_sl[i]);
	}
	if (iTModMethod>11)
	{ // including 1D vertical advection and inflow from macropore infiltration
		temp_c = (2.0*timestep_heat*lambda_sl[i])/(dicken[i]*dicken[i]*heatcap_sl[i]);
		if (qu_sl[i+1]>0)
		{
			temp_c -= (qu_sl[i+1]*cw*timestep_heat)/(dicken[i]*heatcap_sl[i]);
		}
		coeff_a[i]  = -(2.0*timestep_heat*lambda_sl[i-1])/(dicken[i]*(dicken[i]+dicken[i-1])*heatcap_sl[i]);
		if (qu_sl[i-1]<0)
		{
			coeff_a[i] -= (qu_sl[i-1]*cw*timestep_heat)/(dicken[i]*heatcap_sl[i]);
		}
		coeff_c[i]  = 0;
		coeff_b[i]  = 1 - coeff_a[i] + temp_c  - (macroinfiltration[i]*cw*timestep_heat)/(dicken[i]*heatcap_sl[i]);
		right_v[i]  = temper_sl[i]   + Tlow*(temp_c  - (timestep_heat/(dicken[i]*heatcap_sl[i]))*macroinfiltration[i]*cw);
		result_x[i] = 0;
	}
	else
	{   // only heat diffusionion
		temp_c = (2.0*timestep_heat*lambda_sl[i])/(dicken[i]*dicken[i]*heatcap_sl[i]);
		coeff_a[i]  = -(2.0*timestep_heat*lambda_sl[i-1])/(dicken[i]*(dicken[i]+dicken[i-1])*heatcap_sl[i]);
		coeff_c[i]  = 0;
		coeff_b[i]  = 1 - coeff_a[i] + temp_c;
		right_v[i]  = temper_sl[i]   + Tlow*temp_c;
		result_x[i] = 0;
	}
	// now solving the equation system using a tri diagonal matrix solver using the Thomas algorithm
	for (int i=2+offset;i<=levels;i++)
	{
		double m=coeff_a[i]/coeff_b[i-1];
		coeff_b[i] -= (m*coeff_c[i-1]);
		right_v[i] -= (m*right_v[i-1]);
	}

	result_x[levels] = right_v[levels]/coeff_b[levels];

	for (int i=levels-1;i>0+offset;i--)
	{
		result_x[i]=(right_v[i]-coeff_c[i]*result_x[i+1])/coeff_b[i];
	}

}

// calculating the effective thermal conductivity for all soil layers
void unsatzonmodellklasse::CalculateLambdaEff(int offset,double *lambda_sl,int soiltype,double *SE_sl,double *theta_sl, double *thetas, double LambdaWater,double LambdaIce)
{
	// calculate effective thermal conductivity for each soil layer (also in saturated zone)
//	double *thetas = ST[soiltype].getThetas();
//	double *thetasSnow = (do_LayeredSnow)?LayeredSnowProperties[1].getDryDensity():NULL;
	double *thetasSnow = (do_LayeredSnow)?LayeredSnowProperties[1].getThetas():NULL;
//	double *drythermalconduct = ST[soiltype].getDryThermalConduct();
	double lambdathawed;
	double lambdafrozen;

	for (int i=offset;i>-levelsSnow;i--)
	{
		lambda_sl[i]=-9999;
	}
	for (int i=1+offset;i<=levels;i++)
	{
//		double porosity = (i>0)?thetas[i]:(1.0-thetasSnow[levelsSnow+i]/1000.0);
		double porosity = (i>0)?thetas[i]:thetasSnow[levelsSnow+i];
		double theta = theta_sl[i];
		double lambda_air = 0.0262;
		double lambda_soil = (i>0)?ST[soiltype].getDryThermalConduct(i):LayeredSnowProperties[1].getDryThermalConduct(levelsSnow+i);
//		double soilfrac = (i>0)?(1-thetas[i]):(1-thetasSnow[levelsSnow+i]/1000.0); // volume fraction for solid soil material
		double soilfrac = (i>0)?(1-thetas[i]):(1-thetasSnow[levelsSnow+i]); // volume fraction for solid soil material
		double SoilAndWater = theta + soilfrac;
		double lambda_soil_sat  = (i>0)? ST[soiltype].getDrySoilThermalCondEff(i) : LayeredSnowProperties[1].getDrySoilThermalCondEff(levelsSnow+i);
		double lambda_water_sat = (i>0)? ST[soiltype].getWaterThermalCondEff(i,LambdaWater) : LayeredSnowProperties[1].getWaterThermalCondEff(levelsSnow+i,LambdaWater);

		// for dedug: comparison with optimized version
//		lambda_sl[i] = (porosity-theta)*lambda_air + SoilAndWater * pow(lambda_soil,soilfrac/SoilAndWater)
//																  * pow(LambdaWater,(SE_sl[i]*theta)/SoilAndWater)
//																  * pow(LambdaIce,(1-SE_sl[i])*theta/SoilAndWater);
//
		if (SoilAndWater>1.0) SoilAndWater=1.0;
		if (SoilAndWater<0.00001) 
		{
			lambda_sl[i]=lambda_soil;
		}
		else
		{
			if ((porosity-theta)<0.01)	
			{ // almost or completely saturated
				if (SE_sl[i]<0.99)
				{ // saturated, partly ice, partly water -> SoilAndWater will be 1 and can thus be ommitted, air-fraction is zero, thus lambda_air part can be ommitted, soil part is a constzant value which is predefined in soil table
					lambda_sl[i] = lambda_soil_sat * pow(LambdaWater,SE_sl[i]*theta) * pow(LambdaIce,(1-SE_sl[i])*theta);
				}
				else
				{ // saturated, only water
					lambda_sl[i] = lambda_soil_sat * lambda_water_sat;
				}
			} 
			else 
			{ // unsaturated
				if (SE_sl[i]<0.99)
				{ // unsaturated, partly ice, partly water, full equation
					lambda_sl[i] = (porosity-theta)*lambda_air + SoilAndWater * pow(lambda_soil,soilfrac/SoilAndWater)
																			  * pow(LambdaWater,(SE_sl[i]*theta)/SoilAndWater)
																			  * pow(LambdaIce,(1-SE_sl[i])*theta/SoilAndWater);
				}
				else
				{ // unsaturated, only water, ice part (as pow(LambdaIce,(1-SE_sl[i])*theta/SoilAndWater)) is always 1.0 and can be ommitted
					lambda_sl[i] = (porosity-theta)*lambda_air + SoilAndWater * pow(lambda_soil,soilfrac/SoilAndWater)
																			  * pow(LambdaWater,(SE_sl[i]*theta)/SoilAndWater);
				}
			}
		}
	}
}

//#pragma optimize ("p",on)
void unsatzonmodellklasse::pre_run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;

	if (!mpiinitialized)
	{
		if (ist_vorhanden(_MPIRegions_))
			MPIREG=gridinfo[getindex(_MPIRegions_)].pointer;
		else
			MPIREG=NULL;	
		mpiinitialized = true;
	}
	if (!headinitialized)
	{
		char tempstr[200];
		sprintf(tempstr, "%s%d", _gw_head_, 1);

		if (ist_vorhanden(tempstr))
			HEAD1 = gridinfo[getindex(tempstr)].pointer;
		else
			HEAD1 = NULL;
		headinitialized = true;
	}

	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;

	teilmodellklasse::run(me, se, tg, mt, jr, NULL, additionals);

	// if temperature modelling is switched on and the initial conditions were not read in from grids, the actual air temperature is taken as upper boundary condition 
	// and a linear interpolation down to the lower boundary condition is applied now in the first time step
	// SE is not initialized here since it is a result of the heat transfer model and only available for writing the stack after the model ran 
	if (do_TMod)
	{
		if (do_initializeTemper)
		{
#pragma omp parallel for schedule(runtime) 
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j=1;j<=ncols;j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j]==mpi_rank))
					{
						temper[k][j][0]=TLOWBND[k][j];
						double maxdiff = hgeo[k][j][0]-hgeo[k][j][levels];
						for (int lv=1;lv<=levels;lv++)
						{
							// todo: linear interpolieren zwischen TLOWBND und aktueller Obergrenze (T Luft oder alles konstant)
							double actdiff = hgeo[k][j][0]-(hgeo[k][j][lv-1]+hgeo[k][j][lv])/2; // elevation is taken from center of layer
							temper[k][j][lv]=TLOWBND[k][j]*(actdiff/maxdiff)+T[k][j]*(1.0-actdiff/maxdiff);
						}
					}
					else
					{
						for (int lv=1;lv<=levels;lv++) 
						{
							temper[k][j][lv]=nodata;
						}
					}
				}
			}
			do_initializeTemper=false;
		}
		if (do_LayeredSnow)
		{ // for layered snow model, the temperatures of the single snow layers must be initialized independently of the initialization of the soil temperatures (because no pre-Spinup is done for snow)
			if (do_initializeTemperSnow)
			{
#pragma omp parallel for schedule(runtime) 
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							int lr = levelsSnow;
							upperIndexSnow[k][j] = levelsSnow + 1; // if the value of upperIndexSnow is larger than the number of snow levels, there is no snow at all
							while ((swes[k][j][lr] > 0) && (lr>0))
							{
								upperIndexSnow[k][j] = lr;
								lr--;
							}
// ToDo: initialize snow layer temperatures
							double maxdiff = (upperIndexSnow[k][j] <= levelsSnow) ? 
								hgeoSnow[k][j][(int)(upperIndexSnow[k][j])] + thickSnow[k][j][(int)(upperIndexSnow[k][j])] - hgeoSnow[k][j][levelsSnow]:0;
							for (int lv = (int)(upperIndexSnow[k][j]); lv <= levelsSnow; lv++)
							{
								// todo: linear interpolieren zwischen TLOWBND und aktueller Obergrenze (T Luft oder alles konstant)
//								double actdiff = hgeoSnow[k][j][(int)(upperIndexSnow[k][j])-1] - (hgeoSnow[k][j][lv - 1] + hgeoSnow[k][j][lv]) / 2; // elevation is taken from center of layer
								double actdiff = (lv == upperIndexSnow[k][j]) ?
												 thickSnow[k][j][(int)(upperIndexSnow[k][j])]/2.0:
												 hgeoSnow[k][j][(int)(upperIndexSnow[k][j])] + thickSnow[k][j][(int)(upperIndexSnow[k][j])] - hgeoSnow[k][j][lv] + (hgeoSnow[k][j][lv]+hgeoSnow[k][j][lv-1])/2.0;

								temperSnow[k][j][lv] = temper[k][j][1] * (actdiff / maxdiff) + T[k][j] * (1.0 - actdiff / maxdiff);
								if (temperSnow[k][j][lv] > 0) temperSnow[k][j][lv] = 0.0;
							}
							for (int lv = 0; lv < (int)(upperIndexSnow[k][j]); lv++)
							{
								temperSnow[k][j][lv] = nodata;
							}
						}
						else
						{ 
							upperIndexSnow[k][j] = nodata;
							for (int lv = 0; lv <= levelsSnow; lv++)
							{
								temperSnow[k][j][lv] = nodata;
							}
						}
					}
				}
				do_initializeTemperSnow = false; // initialize in first time step only
			}
		}
	}
	/* shifting registers for flow travel time zones */
	for (int z=0;z<=tracer_count;z++)
	{
		for (int i=1;i<=zonmax;i++)
		{
			for (int j=2;j<=(int)maxzeit[i];j++)
			{
				zoneFZT_QD[z][i][j-1].n				= zoneFZT_QD[z][i][j].n;
				zoneFZT_QD[z][i][j-1].RelCellSum	= zoneFZT_QD[z][i][j].RelCellSum;
				zoneFZT_QD[z][i][j-1].sx			= zoneFZT_QD[z][i][j].sx;
				zoneFZT_QI[z][i][j-1].n				= zoneFZT_QI[z][i][j].n;
				zoneFZT_QI[z][i][j-1].RelCellSum	= zoneFZT_QI[z][i][j].RelCellSum;
				zoneFZT_QI[z][i][j-1].sx			= zoneFZT_QI[z][i][j].sx;
			}
			/* only the most distant zone is set to 0 (the other zones contain the surface runoff from the last intervals */
			zoneFZT_QD[z][i][(int)maxzeit[i]].sx = 0;
			zoneFZT_QD[z][i][(int)maxzeit[i]].n  = 0;
			zoneFZT_QD[z][i][(int)maxzeit[i]].RelCellSum  = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].sx = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].n  = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].RelCellSum  = 0;
		}
	}
	/* initializing statistic registers */
	GWSTANDstat->init();
	GWNstat->init();
	MACROstat->init();
	if (methode>2)
	{
		SUBSTEPSstat->init();
	}
	INFILTstat->init();
	SB_1_stat->init();
	SB_2_stat->init();
	if (ICEOLD)	SB_2_IMMO_stat->init();
	ROOTstat->init();
	INFEXstat->init();
	QDstat->init();
	QIstat->init();
	QBstat->init();
	if (gwflow_active)
	{
		GWINstat->init();
		GWEXstat->init();
	}
	Qgstat->init();
	if (HDRAIN) QDRAINstat->init();
	if (MAXPOND) PONDstat->init();
	for (int z=1;z<=tracer_count;z++) tracer_module[z].stat->init();

	if (do_TMod>0)
	{
		thawdepth_stat->init();
	}
	
	// when re-infiltration is on, the riverflow_table must be broadcasted now to all ranks in order to provide all ranks with the proper routed runoff amounts, no matter what river they are working on
	MPIBCastRiverFlowTable();

#ifdef _WASIMDEBUG_
	iteration_stat->init();
#endif

//	QDsum=0;
//	count=0;
}
//#pragma optimize ("p",off)

void unsatzonmodellklasse::setInterceptionBeforeSnow(bool b_interceptionBeforeSnow)
{
	if (do_it == 0) return;
	interceptionBeforeSnow = b_interceptionBeforeSnow;
	char name[500];
	if (!P)
	{
		if (interceptionBeforeSnow) // new approach: first look for snow melt, then for interception model output and then for interpolated precipitation
		{
			sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
			P = gridinfo[getindex_first_of_3(_Schneeausfluss_, name, _Niederschlag_)].pointer;
		}
		else // classic approach: first look for interception model output, then for snow melt and then for interpolated precipitation
		{
			sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
			P = gridinfo[getindex_first_of_3(name, _Schneeausfluss_, _Niederschlag_)].pointer;
		}
	}
}



//#pragma optimize ("p",on)
void unsatzonmodellklasse::run(long me,long se,long tg,long mt,long jr,long actsubstep)
{
	if (do_it == 0) return;

//	{
//		char name[255];
//		sprintf(name,"temperSnow_0.asc");
//		printstackfile(name, gridinfo[int(getindex(_Soil_Temperature_Stack_Snow_))].stackpointer,LayeredSnowProperties[1].getLevels(),true,false);
//	}

	double	bilanz,
			gw_exact,
			deltapond,
			deltasb,
			h_temp,
			maxlateralflow,
			content_upper,
			content_lower, /* masses and water volumes in all saturated layers (and the first groundwater layer) for concentration mixing with groundwater */
			m_err,
		//	max_Error,
			old_gwstand,
			accumThick;
	int		ti,
			maxlayer,
			ziter,
			lay,
			i,
			j,
			k;
	double	icevalue;
	int		TG,
			fzt_class,
			s, // soil type
			start_level,
			old_start_level,
			clay_layer,
			maxlevel,
			startlayer;
	long    weight;
	double  ft,				/* fraction of tracer, which goes into direct runoff */
			trms,
			trmsmax,	/* tracermass - temporary variable for max.conc-calculation in fluxes */
			gw_depth,
			p,
			wurz,
			vbg,
			evap,
			evar,
			etp,
			etr,
			draindist,
			maxrootdepth,
			qlattemp,
			SnowRunoff,
			extraSnowRunoff,
			PreciFromSnowMod;
	double  *thetasat, // --> alias
			*thetares, // --> alias
			*dicken;   // --> alias
	double  *maxfromlr,  // what amount of water could flow out of the layer itself 
			*qufromlr,   // what amount of water could flow out of the layer as downward flux (given by ratio*qu_sl[k]) 
			*iflfromlr,  // what amount of water could flow out of the layer as interflow (given by ration*interflow[k]) 
			*qlatfromlr, // what amount of water could flow out of the layer as lateral flow (given by ratio*qlat[k]) 
			*gwexfromlr, // what amount of water could flow out of the layer as exfiltration (given by ratio*gwex[k])
			*ratio,      // ratio of maxfromlr/(qu_ls[k]+interflow[k]+gwex[k]+qlat[k]) 
			*macrocapacity, // capacity of the macropores in any level
			*macroinfiltration, // amount of infiltrating water in m/s to be used in heat transfer model
			*a_etr,
			*a_evar,
			*etrdetail,
			*a_dicke,
			*theta_temp,
			*h_temp1,
			*theta_average,
			*theta_sl_old, // fuer Speicherung des alten theta-Wertes und der Zwischenwerte in der Iteration im finite differences scheme 
			*q_temp,
			*extractionweight,
			*qlat,  // fluxes as vertical 1-d-arrays (to be summerized to get total fluxes for e.g. interflow) 
			*gwin,
			*gwex,
			*interflow,
			*lateralinflow, 
			*theta_sl,       // --> alias Substitutionsvariable fuer theta[i][j][1..levels] 
			*hgeo_sl,        // --> alias Substitutionsvariable fuer hgeo[i][j][1..levels] 
			*hh_sl,			 // --> alias Substitutionsvariable fuer hh[i][j][1..levels] 
			*qu_sl;          // --> alias Substitutionsvariable fuer qu[i][j][1..levels] 
	double  *balance,		 // not per soil layer but 11 fix values for tracer and water
			*maxbalance;	 // not per soil layer but 11 fix values for tracer and water
	double  *conc_lat,
			*conc_gwin,
			*conc_gwex;
	double	**conc_sl_old,
			**conc_sl;		 // --> alias substitution variable for a single column of ***conc -> conc_sl=conc[i][j] -> for each tracer a vector 
	bool	lastlayer;		 
	double  DT,	timestep, timestep_heat; 
	// variables for heat transport model
	double  *SE_sl,          //	effective saturation (water/ice)
			*T_SE08,	     //	temperature for SE = 0.8, used for Thaw depth calculation
			*lambda_sl,      //	thermal conductivity lambda in J/(m*K)
			*heatcap_sl,     //	heat capacity in J/(Kg*K)
			*qHeat_sl,       //	actual heat flow between soil layers in J/(s*m^2)
			*deltaT_sl,      //	temperature change
			*temper_sl,      // --> alias for temperatures
			*temper_sl_new,  // temperatures
			*dEOverdT_sl,    // change rate of inenr energy
			*coeff_a,		 
			*coeff_b,		 
			*coeff_c,		 
			*right_v,		 
			*result_x;		 // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer

	double evapmaxdepth;
	double actNFactor;
	int startrow = mpistartrow;
	
//	if (surface_interaction > 1)
//	{
//#ifdef _DEBUG
//		if (numTempGridRows > 0)
//		{
//			ofstream rivflow;
//
//			stringstream ss;
//			ss << "unsatzon_rank_" << mpi_rank << "_riv_flow_before.dat";
//			std::string fn("");
//			ss >> fn;
//			rivflow.open(fn.c_str(), std::ios::app);
//			rivflow << "step " << jr << "-" << mt << "-" << tg << "-" << se << endl;
//			for (int i = 1; i <= zonmax; i++)
//			{
//				rivflow << "zone:\t" << i;
//				rivflow << "\tflowm3s:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].flowm3s << "\tinfilt_amount:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].infilt_amount << endl;
//			}
//			rivflow.close();
//		}
//#endif
//	}

	if (surface_interaction > 1)
	{
		MPISyncRivLnkCells1(); // even if not in MPI mode, all RiverLink-Cells will be processed seperately in a new re-ordered grid row. This routine fills that row and sets it's address to mpistartrow-1
//#if defined _MPI && defined _MPIDEBUG
#if defined _MPI 
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		if (numTempGridRows > 0)
		{
			startrow = mpistartrow - 1;
		}
//#if defined _MPI && defined _MPIDEBUG
#if defined _MPI 
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

#pragma omp parallel for schedule(runtime) \
	private(ti,lastlayer,maxlayer,ziter,bilanz, gw_exact,\
	deltapond,deltasb,h_temp,maxlateralflow,content_upper,\
	content_lower,m_err,old_gwstand,lay,i,j,\
	k,icevalue,TG,s,start_level,old_start_level,clay_layer,\
	maxlevel,startlayer, weight,ft,trms,trmsmax,gw_depth, p,\
	wurz,vbg,evap,evar,etp,etr,draindist,maxrootdepth,\
	qlattemp, SnowRunoff,PreciFromSnowMod, thetasat, thetares, dicken,maxfromlr,qufromlr, iflfromlr,\
	qlatfromlr,gwexfromlr,ratio,macrocapacity,macroinfiltration,a_etr,a_evar,etrdetail,a_dicke,				\
	theta_temp,h_temp1,theta_average,theta_sl_old,q_temp,extractionweight,qlat,gwin,			\
	gwex,interflow,lateralinflow,balance,maxbalance,theta_sl,hgeo_sl,hh_sl,					\
	qu_sl,conc_lat,conc_gwin,conc_gwex,conc_sl_old,conc_sl,accumThick,DT,timestep,          \
	SE_sl,T_SE08,lambda_sl,qHeat_sl,heatcap_sl,deltaT_sl,temper_sl_new,temper_sl,dEOverdT_sl,\
	timestep_heat,coeff_a,coeff_b,coeff_c,right_v,result_x,evapmaxdepth,actNFactor)
	for (int i = startrow; i <= mpiendrow; i++)
//	for (int i = 0; i <= mpiendrow && FirstLastCol[i].Last>0; i++)
	{
		a_etr			= 
		a_evar			= 
		etrdetail		= 
		a_dicke			= 
		theta_temp		= 
		h_temp1			= 
		theta_average	= 
		theta_sl_old	= 
		q_temp			= 
		extractionweight = 
		qlat			= 
		gwin			= 
		gwex			= 
		interflow		= 
		lateralinflow	= NULL;

		SE_sl			= 
		T_SE08			= 
		qHeat_sl		= 
		lambda_sl		= 
		heatcap_sl		= 
		deltaT_sl		= 
		dEOverdT_sl		= 
		temper_sl_new	= 
		coeff_a			= 
		coeff_b			= 
		coeff_c			= 
		right_v			= 
		result_x		= 	NULL;

		int arraysize=2*(levels+1+levelsSnow);
		extraSnowRunoff=0;                                                                           
		double *m_thetasat  = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_thetasat[kc]=0;  thetasat  = &(m_thetasat[levelsSnow]);
		double *m_thetares	= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_thetares[kc]=0;  thetares  = &(m_thetares[levelsSnow]);
		double *m_dicken	= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_dicken[kc]=0;    dicken    = &(m_dicken  [levelsSnow]);
		double *m_theta_sl	= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_theta_sl[kc]=0;  theta_sl  = &(m_theta_sl[levelsSnow]);
		double *m_hgeo_sl	= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_hgeo_sl[kc]=0;   hgeo_sl   = &(m_hgeo_sl [levelsSnow]);
		double *m_hh_sl		= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_hh_sl[kc]=0;     hh_sl     = &(m_hh_sl   [levelsSnow]);
		double *m_qu_sl		= new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_qu_sl[kc]=0;     qu_sl     = &(m_qu_sl   [levelsSnow]);
		
		double *m_temper_sl = NULL;
		if (this->do_TMod > 0)
		{
			m_temper_sl = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_temper_sl[kc]=0;  temper_sl = &(m_temper_sl[levelsSnow]);
		}

		double *m_a_etr             = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_a_etr             [kc]=0; a_etr             = &(m_a_etr            [levelsSnow]); // array for real transpiration 
		double *m_a_evar            = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_a_evar            [kc]=0; a_evar            = &(m_a_evar           [levelsSnow]); // array for real evaporation 
		double *m_etrdetail         = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_etrdetail         [kc]=0; etrdetail         = &(m_etrdetail        [levelsSnow]); // array for real evaptranspiration summed up for all vegetation layers 
		double *m_a_dicke           = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_a_dicke           [kc]=0; a_dicke           = &(m_a_dicke          [levelsSnow]); // array for the thickness of the layers 
		double *m_theta_temp        = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_theta_temp        [kc]=0; theta_temp        = &(m_theta_temp       [levelsSnow]); // array for theta-values of the layers (for temporal storing) 
		double *m_h_temp1           = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_h_temp1           [kc]=0; h_temp1           = &(m_h_temp1          [levelsSnow]); // array for hydraulic head values of the last timestep (for repeated use in iteration) 
		double *m_theta_average     = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_theta_average     [kc]=0; theta_average     = &(m_theta_average    [levelsSnow]); // array for averaged theta-values of the layers (for temporal storing during iteration, step 1..i-1) 
		double *m_theta_sl_old      = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_theta_sl_old      [kc]=0; theta_sl_old      = &(m_theta_sl_old     [levelsSnow]); // array for theta-values at the begin of the iteration 
		double *m_q_temp            = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_q_temp            [kc]=0; q_temp            = &(m_q_temp           [levelsSnow]); // array for fluxes 
		double *m_extractionweight  = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_extractionweight  [kc]=0; extractionweight  = &(m_extractionweight [levelsSnow]); // extraction density for the layers within root zone for evapotranspiration 
		double *m_qlat              = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_qlat              [kc]=0; qlat              = &(m_qlat             [levelsSnow]); // array for baseflow per level (traditional, i.e. if gw-module is not active) 
		double *m_gwin              = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_gwin              [kc]=0; gwin              = &(m_gwin             [levelsSnow]); // array for re-infiltration per level 
		double *m_gwex              = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_gwex              [kc]=0; gwex              = &(m_gwex             [levelsSnow]); // array for exfiltration into rivers (baseflow if gwflow_active==1) per level 
		double *m_interflow         = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_interflow         [kc]=0; interflow         = &(m_interflow        [levelsSnow]); // interflow per level 
		double *m_lateralinflow     = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_lateralinflow     [kc]=0; lateralinflow     = &(m_lateralinflow    [levelsSnow]); // lateral inflow (calculated as balance from groundwater module) per level 

		double *m_maxfromlr         = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_maxfromlr         [kc]=0; maxfromlr  = &(m_maxfromlr [levelsSnow]); 
		double *m_qufromlr          = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_qufromlr          [kc]=0; qufromlr   = &(m_qufromlr  [levelsSnow]); 
		double *m_iflfromlr         = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_iflfromlr         [kc]=0; iflfromlr  = &(m_iflfromlr [levelsSnow]); 
		double *m_qlatfromlr        = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_qlatfromlr        [kc]=0; qlatfromlr = &(m_qlatfromlr[levelsSnow]); 
		double *m_gwexfromlr        = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_gwexfromlr        [kc]=0; gwexfromlr = &(m_gwexfromlr[levelsSnow]); 

		/* some variables required for tracer mixing */
		double *m_ratio             = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_ratio             [kc]=0; ratio				= &(m_ratio			   [levelsSnow]);     
		double *m_macrocapacity     = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_macrocapacity     [kc]=0; macrocapacity		= &(m_macrocapacity    [levelsSnow]);
		double *m_macroinfiltration = new double[arraysize]; for (int kc=0;kc<arraysize;kc++) m_macroinfiltration [kc]=0; macroinfiltration	= &(m_macroinfiltration[levelsSnow]);

		conc_lat          = new double[tracer_count+1];
		conc_gwin         = new double[tracer_count+1];
		conc_gwex         = new double[tracer_count+1];
		conc_sl_old       = new double*[tracer_count+1];
		for (int z = 0; z <= tracer_count; z++)
		{
			conc_sl_old[z] = new double[arraysize];  for (int kc=0;kc<arraysize;kc++) conc_sl_old[z][kc]=0;
		}

		conc_sl       = new double*[tracer_count+1];

		balance		  = new double[11];
		maxbalance    = new double[11];

		actNFactor = 1;

		// variables for heat transfer module will be allocated only when heat transfer is used
		double *m_SE_sl = NULL;
		double *m_T_SE08 = NULL;
		double *m_qHeat_sl = NULL;
		double *m_lambda_sl = NULL;
		double *m_heatcap_sl = NULL;
		double *m_deltaT_sl = NULL;
		double *m_dEOverdT_sl = NULL;
		double *m_temper_sl_new = NULL;
		double *m_coeff_a = NULL;
		double *m_coeff_b = NULL;
		double *m_coeff_c = NULL;
		double *m_right_v = NULL;
		double *m_result_x = NULL;

		if (this->do_TMod > 0)
		{
			
			m_SE_sl           = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_SE_sl           [kc]=0; SE_sl           = &(m_SE_sl        [levelsSnow]); // ice saturation
			m_T_SE08          = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_T_SE08          [kc]=0; T_SE08          = &(m_T_SE08       [levelsSnow]); // temperature for SE = 0.8, used for Thaw depth calculation
			m_qHeat_sl        = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_qHeat_sl        [kc]=0; qHeat_sl        = &(m_qHeat_sl     [levelsSnow]); // actual heat flow between soil layers in J/(s*m^2)
			m_lambda_sl       = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_lambda_sl       [kc]=0; lambda_sl       = &(m_lambda_sl    [levelsSnow]); // thermal conductivity lambda in J/(m*K)
			m_heatcap_sl      = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_heatcap_sl      [kc]=0; heatcap_sl      = &(m_heatcap_sl   [levelsSnow]); // heat capacity in J/(Kg*K)
			m_deltaT_sl       = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_deltaT_sl       [kc]=0; deltaT_sl       = &(m_deltaT_sl    [levelsSnow]); // temperature change
			m_dEOverdT_sl     = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_dEOverdT_sl     [kc]=0; dEOverdT_sl     = &(m_dEOverdT_sl  [levelsSnow]); // change rate of inenr energy
			m_temper_sl_new   = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_temper_sl_new   [kc]=0; temper_sl_new   = &(m_temper_sl_new[levelsSnow]); // temperatures
			m_coeff_a         = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_coeff_a         [kc]=0; coeff_a         = &(m_coeff_a      [levelsSnow]); // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer
			m_coeff_b         = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_coeff_b         [kc]=0; coeff_b         = &(m_coeff_b      [levelsSnow]); // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer
			m_coeff_c         = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_coeff_c         [kc]=0; coeff_c         = &(m_coeff_c      [levelsSnow]); // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer
			m_right_v         = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_right_v         [kc]=0; right_v         = &(m_right_v      [levelsSnow]); // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer
			m_result_x        = new double[arraysize];for (int kc=0;kc<arraysize;kc++) m_result_x        [kc]=0; result_x        = &(m_result_x     [levelsSnow]); // vectors for the tridiagonal system of equations to be solved when using the implicit solution for the heat transfer
		}
		else
		{
			SE_sl=NULL; // will be checked in specific_fluxes before calling k_calc
		}

		startlayer=1;
#ifdef _OPENMP
		/*	printf("i=%d --> Thread %i\tTracer_count=%d,levels=%d\n",i,omp_get_thread_num(),tracer_count,levels); */
#endif

		double *ETPsub=0;
		double *ETRsub=0;
		ETPsub=new double[iNumVegetationLayers+1];
		ETRsub=new double[iNumVegetationLayers+1];

		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			//      icevalue=ICE?ICE[i][j]:nodata;
			//      if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue==nodata))

			icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
			if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

			if ((zn[i][j] != nodata) && ((MPIREG[i][j] == mpi_rank) || (i == mpistartrow - 1)) && ((RIVLNK && (RIVLNK[i][j] == nodata)) || !RIVLNK || (i == mpistartrow - 1) || (surface_interaction < 2)) && (icevalue < 1.0) && !(doSpinUp && HRUMaster[i][j] == nodata)) // for unglacierized cells only (at least partly unglacierized)
			//if ((zn[i][j] != nodata) && ((MPIREG[i][j] == mpi_rank) || (i == mpistartrow - 1)) && ((RIVLNK && (RIVLNK[i][j] == nodata)) || !RIVLNK || (i == mpistartrow - 1) ) && (icevalue < 1.0) && !(doSpinUp && HRUMaster[i][j] == nodata)) // for unglacierized cells only (at least partly unglacierized)
			{

				// defining helper variables which will hold absoulte values (like GWEX, GWIN, runoff-components a.s.o
				// for each single subinterval. After a subinterval, these Values are added to the respective grid values
				double  QDRAINsub	= 0,
						GWEXsub		= 0,
						GWINsub		= 0,
						QDgridsub	= 0,
						QDSUsub		= 0,
						QIgridsub	= 0,
						MACROsub	= 0,
						Psub		= 0,
						GWNsub		= 0,
						PERCOLsub	= 0,
						UPRISEsub	= 0,
						EVARsub		= 0,
						IRRIGsub	= 0,
						INFEXsub	= 0,
						INFILTsub	= 0,
						SnowRunoffsub = 0,
						EVAPsub		= 0,
						extraSnowRunoffsub = 0;

				extraSnowRunoff=0;
				if (QDRAIN) QDRAIN[i][j]=0;
				if (GWEX) GWEX[i][j]=0;
				if (GWIN)GWIN[i][j]=0;
				QDgrid[i][j]=0;
				QIgrid[i][j]=0;
				if (MACRO) MACRO[i][j]=0;
				//		for (int v=0;v<=iNumVegetationLayers;v++)
				//		{
				//		    ETR[v][i][j]=0;
				//		}
				GWN[i][j]=0;
				PERCOL[i][j]=0;
				UPRISE[i][j]=0;
				//		EVAR[i][j]=0;
				if (IRRIG) IRRIG[i][j]=0;
				if (INFEX) INFEX[i][j]=0;
				if (INFILT) INFILT[i][j]=0;
				SnowRunoff=0;
				PreciFromSnowMod=0;

				// initializing some variables with 0
				for (ti=0;ti<=ST[0].getLevels();ti++)
				{ 
					a_etr[ti]=0; // set real transpiration to 0 
					a_evar[ti]=0; // set real evaporation to 0 
					a_dicke[ti]=0; // set effective layer thickness to 0 
					etrdetail[ti]=0;
				}
				TG=zarray[(int)zn[i][j]]; /* which subbasin */
				s        = (int)SOIL[i][j];
				if (HDRAIN)
				{
					draindist = DRAINDIST[i][j];
					QDRAINsub=0.0;
				}
				else draindist=0;

				int swesoffset=0;
				double splitfactor_heat=1.0; // subsequently, we will calculate the timestep depending on the actual temperature change
				double splitfactor=1.0; // subsequently, we will calculate the timestep depending on the Courant-Condition or on a parameter read in
				int heatmultiple=1; // this is a reset-counter: every n subtimesteps, the water transport must be computed only, so we need to set an integer max_counter for this maximum here
				double Tsurf;

				// getting adresses for pointer aliases -> for version after 2017-04-14: new arrays m_* are defined which will contain a copy (that needs to be copied back after the loop)
				// of the 3D-matrix values. SInce the snow model stacks will be put on top of the old stacks, the simple reference setting will not work any longer.

				// old versions:
				// dicken = ST[s].getDzui();
				// theta_sl = &(theta[i][j][0]);  // water content               
				// thetasat = ST[s].getThetas();  // saturation theta            
				// thetares = ST[s].getThetar();  // residual theta              
				// hh_sl = &(hh[i][j][0]);     // hydraulic heads             
				// hgeo_sl = &(hgeo[i][j][0]);   // geodaetic heads             
				// qu_sl = &(qu[i][j][0]);     // Fluxes between the layers   
				
				// first the soil properties
				// new versions: values that come from the soil table and doesn't need to be synchronized back to the soil table after that cell is done
				arraysize = sizeof(double)*(levels + 1);
				memmove(&(m_dicken  [levelsSnow]), ST[s].getDzui()  , arraysize); 
				memmove(&(m_thetasat[levelsSnow]), ST[s].getThetas(), arraysize); // saturation theta            
				memmove(&(m_thetares[levelsSnow]), ST[s].getThetar(), arraysize); // residual theta              
				// values that come from the 3d-matrices and needed to be synchronized back to the original stacks after each cell is processed
				memmove(&(m_theta_sl[levelsSnow]), &(theta[i][j][0]), arraysize); // water content               
				memmove(&(m_hh_sl   [levelsSnow]), &(hh[i][j][0])   , arraysize); // hydraulic heads             
				memmove(&(m_hgeo_sl [levelsSnow]), &(hgeo[i][j][0]) , arraysize); // geodaetic heads             
				memmove(&(m_qu_sl   [levelsSnow]), &(qu[i][j][0])   , arraysize); // Fluxes between the layers   
				if (this->do_TMod > 0)
				{
					memmove(&(m_temper_sl[levelsSnow]), &(temper[i][j][0]), arraysize); // Fluxes between the layers   
					memmove(&(m_SE_sl[levelsSnow]), &(SEstack[i][j][0]), arraysize); // SE values of each layer (relative saturation)
				}
				
				// adjusting dicken (thickness) for each layer depending on ExcessIceStack for this cell
				// also adjusting theta_res and theta_sat, which are given relative to the basic soil type but have to be changed according to the extra ice
				// note: hgeo values were already adjusted during initialization
				if (do_ArcVeg)
				{
					// if arvceg is active, not only cells on an ice wedge should have dynamic thicknesses, but also all other cells, because arcveg will change the vegetation mass of the moss layer later
					// this moss layer is represented by the uppermost soil layers (still needs to be defined if a special soil type will be used for this or a new boolean atribute is set to true for such soil types which whould be changed in thickness)
					int s = (int)SOIL[i][j];
					bool isArcVegType = ST[s].getIsArcVegType();
					
					for (int k = 1; k <= levels; k++)
					{
						if ((isArcVegType) || (ExcessIceFactor[i][j][k] != nodata))
						{
							double thickness = hgeo_sl[k - 1] - hgeo_sl[k]; // in m
							dicken[k] = thickness;

							if (ExcessIceFactor[i][j][k] != nodata)
							{
								// ExcessIceFactor contains the extra mount of frozen water in mm per layer
								// here, this extra (excess) value is now integrated into the normal layer by changing the theta-values according to the excess ice amount
								// basic avlues are taken from the soiltable. When all excess ice is molten, the layer uses the default theta-values again but has become thinner by the amount, 
								// the frozen water in excessicefactor encompassed
								// theta_res = 0.05 and theta_sat = 0.4 and dicken = 0.2m and excessice = 100 mm water equivalent means:
								// skeleton = 60%, water = 40% with respect to the remaining part after excess ice was integrated into the layer. 
								// 100mm excess ice WE = 110 mm thickness, i.e. 90mm remain for the normal layer. 5% of 90mm are 4.5 mm, which is only 0.0225% or 200mm
								// also, 40% of 90mm is 36mm, which is 18% of 200mm but must be added to 100mm from excess ise = 136mm -> thetasat = 136/200 = 0.68
								// So the resulting theta_res and theta_sat are 0.0225 and 0.68, resp.
								// the actual theta itself should also be converted this way
								// attention: thickness is taken from hgeo, since dicken-array is not valid for cells with excess-ice

								// limit dicken[k] to the minimum as given by the soil table. hgeo[k-1]-hgeo[k] should as well never become smaller than dicken[k]. If so, some code is wrong (and here, we try to catch that issue)
								if (dicken[k] < ST[s].getDzui()[k])
								{
									dicken[k] = ST[s].getDzui()[k];
								}

								double excessicethickness = ExcessIceFactor[i][j][k] / 0.918 / 1000.0; //in m						
								double extrawaterequivalent = ExcessIceFactor[i][j][k] / 1000.0; //in m						
								double watercontentres = (thickness - excessicethickness)*thetares[k];
								double watercontentsat = (thickness - excessicethickness)*thetasat[k] + extrawaterequivalent; // water content in mm
								double watercontentact = thickness * theta_sl[k]; // water content in mm
								thetares[k] = watercontentres / (thickness);
								thetasat[k] = watercontentsat / (thickness);
								theta_sl[k] = watercontentact / (thickness); // should always be the same or almost the same as theta_sat for ice wedges
							}
						}
					}
				}

				// second the snow properties
				if (do_LayeredSnow)
				{
					arraysize = sizeof(double)*(levelsSnow+1);
					// values that come from the LayeredSnowProperties table and doesn't need to be synchronized back to the soil table after that cell is done
					memmove(&(m_dicken  [0]), this->LayeredSnowProperties[1].getDzui()  , arraysize);
					memmove(&(m_thetasat[0]), this->LayeredSnowProperties[1].getThetas(), arraysize);
					memmove(&(m_thetares[0]), this->LayeredSnowProperties[1].getThetar(), arraysize);
					// values that come from the 3d-matrices and needed to be synchronized back to the original stacks after each cell is processed
					memmove(&(m_theta_sl[0]), &(thetaSnow[i][j][0]), arraysize);   
					memmove(&(m_hh_sl   [0]), &(hhydSnow[i][j][0]) , arraysize);   
					memmove(&(m_hgeo_sl [0]), &(hgeoSnow[i][j][0]) , arraysize);   
					memmove(&(m_qu_sl   [0]), &(quSnow[i][j][0])   , arraysize);   
					if (this->do_TMod > 0)
					{
						memmove(&(m_temper_sl[0]), &(temperSnow[i][j][0]), arraysize);
//						memmove(&(m_SE_sl[0]), &(SE[i][j][0]), arraysize);
					}

				}

				// second the snow properties
				if (do_LayeredSnow)
				{
					// set dicken for snow layers, in particular for the uppermost snow layer, since the thickness can be variable
					// existing stacks: 
					// thickSnow, densSnow, upperIndexSnow -> should be updated now
					upperIndexSnow[i][j]=levelsSnow+1;
					
					for (int k=levelsSnow;k>=1;k--)
					{
						if (swes[i][j][k]>0)
						{
							upperIndexSnow[i][j] = k;
							if (densSnow[i][j][k] <= -9999)
							{
								densSnow[i][j][k]    = LayeredSnowProperties[1].getDryDensity()[k]/1000.0; // result in kg/dm^3 or t/m^3
							}
							if (thickSnow[i][j][k] <= -9999)
							{
								thickSnow[i][j][k]   = ((swes[i][j][k]+swel[i][j][k])/densSnow[i][j][k])/1000.0; // result in m
							}
							dicken[k-levelsSnow] = thickSnow[i][j][k];
//							theta_sl[k-levelsSnow] = thetares[k-levelsSnow]+(swes[i][j][k]+swel[i][j][k])/(1000.0*thickSnow[i][j][k]); // liquid and solid put fraction together here
//							theta_sl[k-levelsSnow] = (swel[i][j][k])/(1000.0*thickSnow[i][j][k]); // liquid and solid put fraction together here
							theta_sl[k-levelsSnow] = (swel[i][j][k]+swes[i][j][k])/(1000.0*thickSnow[i][j][k]); // liquid and solid put fraction together here
							if (this->do_TMod > 0)
							{
								if (temper_sl[k-levelsSnow]<=-9999)
								{
									temper_sl[k-levelsSnow]=temper_sl[k-levelsSnow+1]; // if sublimation took place, temperature may not have been changed in snowmodel
									if (temper_sl[k-levelsSnow]>0) temper_sl[k-levelsSnow]=0;
								}
							}
						}
						else
						{
							densSnow[i][j][k]    = -9999;
							thickSnow[i][j][k]   = -9999;
							if (this->do_TMod > 0)
							{
								temperSnow[i][j][k] = -9999;
								temper_sl[k-levelsSnow] = -9999;
							}
							dicken[k-levelsSnow] = -9999;
							theta_sl[k-levelsSnow] = -9999;
						}
					}
				}

				evapmaxdepth = ST[s].getEvapMaxDepth();

				for (int z=1;z<=tracer_count;z++) 
				{
					conc_sl[z]=&(conc[z][i][j][0]); /* tracer concentrations for each tracer species */
				}

				// first usage of P after this line. So we need to get the correct input for both cases, either usual sub model sequence with input unsatzon = output snow model or interception model
				// or input unsatzon = sum of snow melt (stored in SWEL[i][j][0] if no heat transfer is active) and liquid precipitation fraction (P_RAIN)
				if ((do_LayeredSnow) && (actsubstep==1))
				{
					// swel[i][j][0] is the snow melt computed in snowmodel. It is 0 if the heat transfer model is used (then, here in unsatzon model the snow melt will be calculated and put in swel[i][j][0])
					// P_RAIN[i][j] is the liquid fraction of precipitation as input to the snow model
					P[i][j]=swel[i][j][0]+P_RAIN[i][j];
				}

				if ((this->do_TMod > 0)&&(do_LayeredSnow))
				{
					while ((swes[i][j][levelsSnow+swesoffset]>0)&&(-swesoffset<levelsSnow))
					{
						swesoffset--;
					}
				}

#ifdef _WASIMDEBUG_
				if (QDSU!=NULL)
					balance[0]=(P[i][j]+ETR[0][i][j]+EVAR[i][j])/substepcnt+QDSU[i][j]; /* if here ETRsub[0] and/or EVARsub are positive, than it is evaporation from Snow, which is subtracted again later */
				else
					balance[0]=(P[i][j]+ETR[0][i][j]+EVAR[i][j])/substepcnt; /* if here ETRsub[0] and/or EVARsub are positive, than it is evaporation from Snow, which is subtracted again later */
				
				for (int z = 1; z <= tracer_count; z++)
				{
					balance[z] = P[i][j] / substepcnt*tracer_module[z].cinput[i][j];
				}
				deltapond=0;
				deltasb=0;
				if (MAXPOND) 
				{
					deltapond+=POND[i][j];
					for (int z = 1; z <= tracer_count; z++)
					{
						balance[z] += tracer_info[z].evaporating ? POND[i][j] * conc_sl[z][0] : conc_sl[z][0];
					}
				}
				for (int kz=1+swesoffset;kz<=levels;kz++)
				{
					deltasb+=(theta_sl[kz]*dicken[kz]*1000.0);
					for (int z = 1; z <= tracer_count; z++)
					{
						balance[z] += tracer_info[z].evaporating ? theta_sl[kz] * dicken[kz] * 1000 * conc_sl[z][kz] : conc_sl[z][kz];
					}
				}
				for (int z = 1; z <= tracer_count; z++)
				{
					if (tracer_info[z].radioactive)
					{
						balance[z] *= tracer_info[z].decfact;
					}
				}

				/* because radioactive decay in groundwater is considered in the groundwater-module, the balance for groudnwater 
	has to be added after decay-calculation for the unsaturated zone, not before! */
				if (gwflow_active)
				{
					for (int z = 1; z <= tracer_count; z++)
					{
						balance[z] += thetasat[levels] * (gw_thick[i][j] - (hgeo_sl[0] - hgeo_sl[levels]))*gw_conc[z][i][j] * 1000;
					}
				}
				//            balance[z]+=thetasat[levels]*(gw_thick[i][j]-(levels-1)*dicken[levels])*gw_conc[z][i][j]*1000;
#endif

				// make a copy of the actual soil temperature profile in order to use the original values in each iteration 
				// the values in temper_sl_new will be returned from the heat transport model and copied into temper_sl after the last iteration
				if (this->do_TMod > 0)
				{
					// temper_sl = &(temper[i][j][0]);     was done above already when copying the actual stack values to the local variables for all the other entities like theta, qu etc.
					memmove(&(m_temper_sl_new[0]), &(m_temper_sl[0]), sizeof(double)*(levels + 1 + levelsSnow));
				}

				// the temperature model is now processed one single time to estimate the split factor with regard to the heat transport.
				// in case the time step does not have to be splitted, the results are already valid for the entire time step, 
				// otherwise they are valid for the first sub time step, so the core iteration will do the heat transfer model only beginning with the second sub step

				if (this->do_TMod > 0)
				{
					DT = InV*60.0; // without splitfactor
					ST[s].CalcSEanddEOverdT(0,SE_sl, dEOverdT_sl, temper_sl, this->dTModFactor, theta_sl, thetares, thetasat); // calculating new relative saturation values (water=1...ice=0)
					if (do_LayeredSnow) LayeredSnowProperties[1].CalcSEanddEOverdT(-levelsSnow,SE_sl, dEOverdT_sl, temper_sl, this->dTModFactor, theta_sl, thetares, thetasat); // calculating new relative saturation values (water=1...ice=0)
					CalculateLambdaEff(swesoffset,lambda_sl, s, SE_sl, theta_sl, thetasat,this->dTModLamWat, this->dTModLamIce); // effective thermal conductivity
					// dEOverdT is yet only a scaling factor. Must be multiplied by Rho*LatHeat*(theta-theta_res)
					Expand_dEoverdT(swesoffset,theta_sl,thetares,dEOverdT_sl);
					// calculating actual heat capacity
					ST[s].CalcHeatCapacity(0,temper_sl_new,heatcap_sl,SE_sl,dEOverdT_sl,theta_sl,this->dTModRhoWater,this->dTModCWater,this->dTModCIce0,this->dTModCIce20);
					if (do_LayeredSnow) LayeredSnowProperties[1].CalcHeatCapacity(-levelsSnow,temper_sl_new,heatcap_sl,SE_sl,dEOverdT_sl,theta_sl,this->dTModRhoWater,this->dTModCWater,this->dTModCIce0,this->dTModCIce20);
	
					// calculating heat fluxes
					int use = (int)LUSE[1][i][j];
					actNFactor = LUtable[use].GetNfactor(WSstream,tag,monat,jahr);
					if (actNFactor<-9999)
					{
						if (T[i][j] < 0) 
						{ 
							actNFactor = LUtable[use].GetNFactorFreezing(WSstream);
							if (actNFactor<-9999)
							{
								actNFactor = dTModNFactorFreezing;
							}
						} 
						else 
						{
							actNFactor = LUtable[use].GetNFactorThawing(WSstream);
							if (actNFactor<-9999)
							{
								actNFactor = dTModNFactorThawing;
							}
						}
					}
					Tsurf = T[i][j];
					if (EBAL && do_LayeredSnow)
					{
						if (upperIndexSnow[i][j]<=levelsSnow) // change Tsurf only if there is some snow on the ground
						{
							if ((EBAL[i][j]==0)&&(SURFTEMP[i][j]<273.16))
							{	// cooling conditions
								Tsurf = SURFTEMP[i][j]-273.16;
								if (Tsurf<T[i][j]-13.0)
								{
									Tsurf=T[i][j]-13.0;
								}
								actNFactor = dTModNFactorFreezing;
							}	else if (EBAL[i][j]>0)	{
								// melting conditions: EBAL > 0
								// mixing temperature of the uppermost two layers
								double Tmix;
								double thicksum;
								if (swesoffset<-1)
								{
									// at least two snow layers: in case the upper layer gets too thin, we look at the two uppermost layers to avoid unrealistic surface temperatures
									double swel1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*SE_sl[1+swesoffset]; //swes[i][j][1+swesoffset+levelsSnow];
									double swes1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*(1-SE_sl[1+swesoffset]); //swel[i][j][1+swesoffset+levelsSnow];
									double swel2 = 1000*dicken[2+swesoffset]*theta_sl[2+swesoffset]*SE_sl[2+swesoffset]; //swes[i][j][2+swesoffset+levelsSnow];
									double swes2 = 1000*dicken[2+swesoffset]*theta_sl[2+swesoffset]*(1-SE_sl[2+swesoffset]); //swel[i][j][2+swesoffset+levelsSnow];

									Tmix = ((swes1*dTModCIce0+swel1*dTModCWater) * temper_sl[1+swesoffset] 
											+ (swes2*dTModCIce0+swel2*dTModCWater) * temper_sl[2+swesoffset]) / 
											((swes1+swes2)*dTModCIce0 + (swel1+swel2)*dTModCWater);
									thicksum=dicken[1+swesoffset]+dicken[2+swesoffset];
								}
								else
								{
									// only one snow layer, so we may use that layers properties or mix it with the first soil layer or we reduce the energy balance to an amount that is sufficient to melt the actual content of the last snow layer only
									double swes1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*SE_sl[1+swesoffset]; //swes[i][j][1+swesoffset+levelsSnow];
									double swel1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*(1-SE_sl[1+swesoffset]); //swel[i][j][1+swesoffset+levelsSnow];
									double requiredMeltEnergy = (337500.0*swes1 - (temper_sl[0] * (2100.0 * swes1 + 4187.0 * swel1)))/(InV*60.0);
									if (EBAL[i][j]>requiredMeltEnergy) EBAL[i][j]=requiredMeltEnergy;
									Tmix = temper_sl[0];
									thicksum = dicken[0];
								}
								Tsurf = Tmix + EBAL[i][j]*thicksum/(2.0*lambda_sl[1+swesoffset]);
								if (Tsurf<0)
								{
									Tsurf=0;
								} else if (Tsurf>T[i][j]+13) {
									Tsurf = T[i][j]+13;
								}
//										Tsurf = temper_sl[1+swesoffset]+EBAL[i][j]*dicken[1+swesoffset]/(2.0*lambda_sl[1+swesoffset]);
								actNFactor = dTModNFactorThawing;
							} else {
								Tsurf = temper_sl[1+swesoffset]; // sublimation conditions, i.e. all energy is required for snow evaporation, no heat flux into snow possible
								actNFactor = 1.0;
							}

							// actNFactor=1;
						}
					}
					// old version of setting Tsurf
//					if (EBAL && do_LayeredSnow)
//					{
//						if (upperIndexSnow[i][j]<=levelsSnow) // change Tsurf only if there is some snow on the ground
//						{
//							if ((EBAL[i][j]==0)&&(SURFTEMP[i][j]<273.16))
//							{	// cooling conditions
//								Tsurf = SURFTEMP[i][j]-273.16;
//								actNFactor = dTModNFactorFreezing;
//							}
//							else
//							{
//								// melting conditions: EBAL > 0
//								Tsurf = temper_sl[1+swesoffset]+EBAL[i][j]*dicken[1+swesoffset]/(2.0*lambda_sl[1+swesoffset]);
//								actNFactor = dTModNFactorThawing;
//							}
//						}
//					}

//					SolveTriDiag(swesoffset,i,j,lambda_sl,heatcap_sl,dicken,temper_sl,timestep_heat,coeff_a,coeff_b,coeff_c,
//									right_v,result_x,actNFactor,macroinfiltration,qu_sl,Tsurf);
//					if (this->do_TMod > 0)
//					{
//						for (int z=1+swesoffset;z<=levels;z++)
//						{
//							temper_sl_new[z]=result_x[z];
//						}
//					}

					// calculating heat fluxes
					CalcHeatFluxes(swesoffset,qHeat_sl,lambda_sl,dicken,temper_sl_new,TLOWBND[i][j],T[i][j],actNFactor,Tsurf); // upper boundary condition: air temperature, lower boundary condition: given in boundary condition grid
					// calculating temperature change
					CalcDeltaT(swesoffset,deltaT_sl, qHeat_sl,theta_sl,dicken,heatcap_sl,s,(int)DT);
					// calculating new temperatures
					CalcNewTemperature(swesoffset,temper_sl_new,temper_sl,deltaT_sl,qHeat_sl,TLOWBND[i][j],T[i][j],false,Tsurf);

					// restricting temperature change to a value which does not break the second law of thermodynamics (no gradient sign switch!)
					// here, a split factor for the time step will be calculated which must then be compared with the split factor for water transport (take the maximum value)
					// we check the sub time step for layers which are no deeper than 2m only (dynamic is much slower below 0.5m, 
					// so the time step splitting is probably never neccessary for depth below 0.5m
					// 1st check: uppermost layer (with air contact)
					double depth=0;
					
					double maxDeltaT = (T[i][j]-temper_sl[1+swesoffset])/4.0;
					while ((fabs(maxDeltaT * splitfactor_heat) < fabs(deltaT_sl[1+swesoffset])) && (DT/splitfactor_heat>dTModShortestStep) || (DT/splitfactor_heat>dTModLongestStep))//&& fabs(deltaT_sl[1])>0.01)
					{
							splitfactor_heat ++; 
					}

					for (int i=1+swesoffset;i<levels;i++) 
					{
						depth+=dicken[i];
						maxDeltaT = (temper_sl[i]-temper_sl[i+1])/4.0;
						while ((fabs(maxDeltaT * splitfactor_heat) < fabs(deltaT_sl[i+1])) && (DT/splitfactor_heat>dTModLongestStep) || (DT/splitfactor_heat>dTModLongestStep))//&& fabs(deltaT_sl[i+1])>0.01)
						{
								splitfactor_heat ++; 
						}
					}

//					(*WSstream) << "Splitfactor_Fm: "<<splitfactor_heatFM<<", Splifatcor: "<<splitfactor_heat<<"[end-entry]";
//
//
					for (int i=swesoffset;i<=levels;i++)
					{
						temper_sl_new[i]=temper_sl[i]; 
					}
				}

				if (methode>2)
				{
					// calculating the actual splitfactor
					// for each layer, the optimal timestep is calculated. The shortest one is then taken (only if the 
					// setting here DT and timestep for subintervals for dynamic timestepping
					int cour_lev=1+swesoffset;
					double dt_min=InV/substepcnt*60.0; // normal time step, longer values are not possible
					do
					{
						double k_rel=1.0;
						// only if it's not raining, the actual k-value will be consistent for a Courant-condition-checking
						// otherwise take the saturated k-value (with a k_rel of 1.0) for estimating the Courant-Condition
						if (P[i][j]==0)
						{
							double SE=(SE_sl!=NULL)?SE_sl[cour_lev]:1.0;
							// limitation of the effect of SE. Since a low SE may lead to very long time steps but some soils with extremely high conductivity may react even stronger 
							// to small changes of theta, we need to get a time step short enough to regard that effects
							// So, an extreme measure could be to even not to regard SE at al here and always set the sub time step as if the soil were completely unfrozen, leading to longer run times but better results.
							if (SE<0.5) SE=0.5; 
							k_rel = (cour_lev>0)? ST[s].k_calc(cour_lev,theta_sl[cour_lev],SE, thetares[cour_lev], thetasat[cour_lev]) : LayeredSnowProperties[1].k_calc(levelsSnow+cour_lev,theta_sl[cour_lev],SE);   //k_relativ, oben
						}
						if (k_rel>0)
						{
							double dt = (cour_lev>0) ? dicken[cour_lev]/(k_rel*ST[s].getK_sat()[cour_lev]) : dicken[cour_lev]/(k_rel*LayeredSnowProperties[1].getK_sat()[levelsSnow+cour_lev]); // recommended time step in seconds
							if (dt < dt_min)
							{
								dt_min=dt;
							}
						}
						cour_lev++;
					}
					while (cour_lev<GWLEVEL[i][j]);

					if (dt_min<dParameterMinStep)
					{
						dt_min=dParameterMinStep;
					}

					splitfactor = ceil(InV/substepcnt*60.0/dt_min);
				}
// splitfactor = 10; Test for heat transfer model
//				int heatmultiple_counter=0;
				// now we try two different split_factors
				if (splitfactor_heat>2*splitfactor)
				{
					// old: splitfactor=splitfactor_heat;
					splitfactor_heat = ((int)((int)splitfactor_heat/(int)splitfactor)+1)*splitfactor;
					heatmultiple = (int)(splitfactor_heat/splitfactor);
				}
				else
				{ // no special splitting between water and heat transfer when the difference in splitfactors is small
					if (splitfactor_heat>splitfactor)
					{ // splitfactor for heat transport is larger than for water transport, and the larger one wins
						splitfactor=splitfactor_heat;
					}
					else
					{ // splitfactor for heat transport is smaller than for water transport, so the larger one wins here, too
						splitfactor_heat=splitfactor;
					}
					heatmultiple = 1;
				}

// actually, this is a shortcut: setting the temperatur to the value calculated in the sub-step-calculation neglects the impact of all water fluxes (Evaporation, infiltration etc.)
// and represents the temperature profile rather before than after the time step. The differences are small but visible, so it is recommended to do the calculation again in the main iteration cycle.
//				if ((do_TMod>0)&&(splitfactor == 1)) // if the interval is not splitted, then the heat transfer calculated above is already the final result, otherwise it will be recalculated in the core iteration
//				{
//					for (int i=0;i<=levels;i++)
//					{
//						// make a copy of the actual soil temperature profile in order to use the original values in each iteration 
//						// the values in temper_sl_new will be returned from the heat transport model and copied into temper_sl after the last iteration
//						temper_sl[i]=temper_sl_new[i]; 
//					}
//				}

				DT = InV/substepcnt*60.0/splitfactor;           /* interval in secounds (InV: interval in minutes) */
				timestep = 6.0e4*InV/substepcnt/splitfactor;    /* factor for changing the measure of fluxes from m/s to mm/timestep */
				timestep_heat = 60.0*InV/substepcnt/splitfactor_heat;    /* timestep for heat fluxes in seconds */

				// input-flow-Values must splitted into the sub intervals
				Psub=(P[i][j]/substepcnt)/splitfactor;
				EVAPsub=(EVAP[i][j]/substepcnt)/splitfactor;
				//must be done within the splitfactor-loop because ETPsub[layer] will be set to 0 if ETP is taken from POND
				//		for (int v=0;v<=iNumVegetationLayers;v++)
				//		{
				//			ETPsub[v]=ETP[v][i][j]/splitfactor;
				//			ETRsub[v]=ETR[v][i][j]; //e.g. due to snow evaporation there could already be a positive value in ETR
				//		}

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// start of a sub interval
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//				int printflag=0;

				for (int splitcnt=1;splitcnt<=splitfactor;splitcnt++)
				{
					Dumptype dump;
					if (this->dumpmode==true)
					{
						dump.ETR			= new double[iNumVegetationLayers+1];
						dump.gwin			= new double[levels+1+levelsSnow];
						dump.gwex			= new double[levels+1+levelsSnow];
						dump.h_temp1		= new double[levels+1+levelsSnow];
						dump.theta_sl_old	= new double[levels+1+levelsSnow];
						dump.theta_average	= new double[levels+1+levelsSnow];
						dump.hh_sl			= new double[levels+1+levelsSnow];
						dump.hgeo_sl		= new double[levels+1+levelsSnow];
						dump.SE_sl			= new double[levels+1+levelsSnow];
						dump.qu_sl			= new double[levels+1+levelsSnow];
						dump.interflow		= new double[levels+1+levelsSnow];
						dump.thetasat		= new double[levels+1+levelsSnow];
						dump.dicken			= new double[levels+1+levelsSnow];
						dump.theta_temp		= new double[levels+1+levelsSnow];
						dump.thetares		= new double[levels+1+levelsSnow];
						dump.temper_sl		= new double[levels+1+levelsSnow];
						dump.theta_sl		= new double[levels+1+levelsSnow];
						dump.lambda_sl		= new double[levels+1+levelsSnow];
						dump.dEOverdT_sl	= new double[levels+1+levelsSnow];
						dump.heatcap_sl		= new double[levels+1+levelsSnow];
						dump.qHeat_sl		= new double[levels+1+levelsSnow];
						dump.temper_sl_new	= new double[levels+1+levelsSnow];
						dump.swes			= new double[levelsSnow+1];
						dump.swel			= new double[levelsSnow+1];
						dump.densSnow		= new double[levelsSnow+1];
						dump.thickSnow		= new double[levelsSnow+1];
						dump.LUSE           = new double[iNumVegetationLayers+1];
					
	#ifdef _READDUMPFILE
						bool readdump=true;
	#else
						bool readdump=false;
	#endif

						if (readdump)
						{
							std::ifstream dumpstream;
							dumpstream.open("dump.dmp",std::ios::binary);
							dumpstream.read((char *)(&maxlayer),sizeof(int));
							dumpstream.read((char *)(&h_temp),sizeof(double));
							dumpstream.read((char *)(&p),sizeof(double));
							dumpstream.read((char *)(&(DELTAGW[i][j])),sizeof(double));
							for (int y=0;y<=iNumVegetationLayers;y++)	
								dumpstream.read((char *)(&(ETR[y][i][j])),sizeof(double));
							dumpstream.read((char *)(&gwin[0-levelsSnow]),				sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&gwex[0-levelsSnow]),				sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&h_temp1[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&theta_sl_old[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&theta_average[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&hh_sl[0-levelsSnow]),				sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&hgeo_sl[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&qu_sl[0-levelsSnow]),				sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&interflow[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&thetasat[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&dicken[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&theta_temp[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&thetares[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							dumpstream.read((char *)(&theta_sl[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
							if (do_TMod) {
								dumpstream.read((char *)(&SE_sl[0-levelsSnow]),			sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&temper_sl[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&lambda_sl[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&dEOverdT_sl[0-levelsSnow]),	sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&heatcap_sl[0-levelsSnow]),	sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&qHeat_sl[0-levelsSnow]),		sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&temper_sl_new[0-levelsSnow]),	sizeof(double)*(levels+1+levelsSnow));
								dumpstream.read((char *)(&timestep_heat),				sizeof(double));
							}
							if (do_LayeredSnow)
							{
								dumpstream.read((char *)(&swes[i][j][0]),		sizeof(double)*(levelsSnow+1));
								dumpstream.read((char *)(&swel[i][j][0]),		sizeof(double)*(levelsSnow+1));
								dumpstream.read((char *)(&densSnow[i][j][0]),	sizeof(double)*(levelsSnow+1));
								dumpstream.read((char *)(&thickSnow[i][j][0]),	sizeof(double)*(levelsSnow+1));
								dumpstream.read((char *)(&upperIndexSnow[i][j]),sizeof(double));
								dumpstream.read((char *)(&SSNO[i][j]),			sizeof(double));
								dumpstream.read((char *)(&SLIQ[i][j]),			sizeof(double));
								dumpstream.read((char *)(&SSTO[i][j]),			sizeof(double));
							}
							dumpstream.read((char *)(&s),sizeof(int));
							for (int y=1;y<=iNumVegetationLayers;y++) 
								dumpstream.read((char *)(&(LUSE[y][i][j])),sizeof(double));
							dumpstream.read((char *)(&me), sizeof(int));
							dumpstream.read((char *)(&se), sizeof(int));
							dumpstream.read((char *)(&tg),sizeof(int));
							dumpstream.read((char *)(&mt),sizeof(int));
							dumpstream.read((char *)(&jr),sizeof(int));
							dumpstream.read((char *)(&splitfactor),sizeof(double));
							dumpstream.read((char *)(&DT),sizeof(double));
							dumpstream.read((char *)(&swesoffset),sizeof(int));
							if (do_TMod) {
								if (!dumpstream.eof())
								{
									dumpstream.read((char *)(&splitfactor_heat), sizeof(double));
									dumpstream.read((char *)(&heatmultiple), sizeof(int));
									dumpstream.read((char *)(&(T[i][j])), sizeof(int));
									dumpstream.read((char *)(&Tsurf),sizeof(double));
								}
								else
								{
									splitfactor_heat = 60.0*InV / (substepcnt*timestep_heat);
									heatmultiple = (int)(splitfactor_heat / splitfactor);
								}
							}
							dumpstream.close();
						}
						// else read values should be moved to te dump-array in order to allow a proper comparison after core routines that simulates the decision if a dump file will be written
						{
							memmove(&(dump.maxlayer				),&(maxlayer				),sizeof(int));
							memmove(&(dump.h_temp				),&(h_temp					),sizeof(double));
							memmove(&(dump.p					),&(p						),sizeof(double));
							memmove(&(dump.DELTAGW				),&(DELTAGW[i][j]			),sizeof(double));
							for (int y=0;y<=iNumVegetationLayers;y++)
								memmove(&(dump.ETR[y]			),&(ETR[y][i][j]			),sizeof(double));
							memmove(&(dump.gwin[0]				),&(gwin[0-levelsSnow]					),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.gwex[0]				),&(gwex[0-levelsSnow]					),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.h_temp1[0]			),&(h_temp1[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.theta_sl_old[0]		),&(theta_sl_old[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.theta_average[0]		),&(theta_average[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.hh_sl[0]				),&(hh_sl[0-levelsSnow]					),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.hgeo_sl[0]			),&(hgeo_sl[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.qu_sl[0]				),&(qu_sl[0-levelsSnow]					),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.interflow[0]			),&(interflow[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.thetasat[0]			),&(thetasat[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.dicken[0]			),&(dicken[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.theta_temp[0]		),&(theta_temp[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.thetares[0]			),&(thetares[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							memmove(&(dump.theta_sl[0]			),&(theta_sl[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
							if (do_TMod)
							{
								memmove(&(dump.SE_sl[0]				),&(SE_sl[0-levelsSnow]				),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.temper_sl[0]			),&(temper_sl[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.lambda_sl[0]			),&(lambda_sl[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.dEOverdT_sl[0]		),&(dEOverdT_sl[0-levelsSnow]		),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.heatcap_sl[0]		),&(heatcap_sl[0-levelsSnow]		),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.qHeat_sl[0]			),&(qHeat_sl[0-levelsSnow]			),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.temper_sl_new[0]		),&(temper_sl_new[0-levelsSnow]		),sizeof(double)*(levels+1+levelsSnow));
								memmove(&(dump.timestep_heat		),&(timestep_heat			),sizeof(double));
								memmove(&(dump.splitfactor_heat		),&(splitfactor_heat		),sizeof(double));
								memmove(&(dump.heatmultiple			),&(heatmultiple			),sizeof(int));
								memmove(&(dump.T					),&(T[i][j]					),sizeof(double));
								memmove(&(dump.Tsurf				),&(Tsurf					),sizeof(double));

							}
							if (do_LayeredSnow)
							{
								memmove(&(dump.swes[0]),&(swes[i][j][0]),sizeof(double)*(levelsSnow+1));
								memmove(&(dump.swel[0]),&(swel[i][j][0]),sizeof(double)*(levelsSnow+1));
								memmove(&(dump.densSnow[0]),&(densSnow[i][j][0]),sizeof(double)*(levelsSnow+1));
								memmove(&(dump.thickSnow[0]),&(thickSnow[i][j][0]),sizeof(double)*(levelsSnow+1));
								memmove(&(dump.upperIndexSnow),&(upperIndexSnow[i][j])   ,sizeof(double));
								memmove(&(dump.SSNO),&(SSNO[i][j])   ,sizeof(double));
								memmove(&(dump.SLIQ),&(SLIQ[i][j])   ,sizeof(double));
								memmove(&(dump.SSTO),&(SSTO[i][j])   ,sizeof(double));
							}
							memmove(&(dump.s					),&(s						),sizeof(int));
							for (int y=1;y<=iNumVegetationLayers;y++)
								memmove(&(dump.LUSE[y]				),&(LUSE[y][i][j]			),sizeof(double));

							memmove(&(dump.stunde				),&(stunde					), sizeof(int));
							memmove(&(dump.tag					),&(tag						),sizeof(int));
							memmove(&(dump.monat				),&(monat					),sizeof(int));
							memmove(&(dump.jahr					),&(jahr					),sizeof(int));
							memmove(&(dump.splitfactor			),&(splitfactor				),sizeof(double));
							memmove(&(dump.DT					),&(DT						),sizeof(double));
							memmove(&(dump.swesoffset			),&(swesoffset				),sizeof(int));
						}
					}

					for (int v=0;v<=iNumVegetationLayers;v++)
					{
						ETPsub[v]=(ETP[v][i][j]/substepcnt)/splitfactor;
						ETRsub[v]=0; //e.g. due to snow evaporation there could already be a positive value in ETR, but this is added later to the already existing ETR, which is not changed anymore here compared to earlier versions
					}
					Psub=(P[i][j]/substepcnt)/splitfactor;
					EVAPsub=(EVAP[i][j]/substepcnt)/splitfactor;
					if (QDSU)
					{
//						QDSUsub=(QDSU[i][j]/substepcnt)/splitfactor;
						QDSUsub=(QDSU[i][j])/splitfactor; // since QDSU is assumed to be a storage and no flux, we use the entire amount here
					}
					else
					{
						QDSUsub=0;
					}
					EVARsub=0;
					QDRAINsub=0;
					GWEXsub=0;
					GWINsub=0;
					QDgridsub=0;
					QIgridsub=0;
					MACROsub=0;
					GWNsub=0;
					PERCOLsub=0;
					UPRISEsub=0;
					IRRIGsub=0;
					INFEXsub=0;
					INFILTsub=0;
					extraSnowRunoffsub = 0;

					//			if ((i==10)&&(j==30))
					//			(*WSstream) << "Stunde " << stunde << " splitcnt " <<splitcnt <<  " GWSTAND " << GWSTAND[i][j] << "[end-entry]";

					//oldplace of pointer substitution for theta_sl a.s.o.

					//		printf("Thread:%d, i:%d ,j:%d ,tsat:%lf ,tres:%lf ,d:%lf ,thet:%lf ,hgeo:%lf ,hh:%lf ,qu:%lf\n",omp_get_thread_num(),i,j,thetasat[1],thetares[1],dicken[1],theta_sl[1],hgeo_sl[1],hh_sl[1],qu_sl[1]);

					if (gwflow_active)
					{
						GWEXsub=0;
						GWINsub=0;
					}


					//#if 0
					// old use	for (z=1;z<=tracer_count;z++) conc_sl[z]=&(conc[z][i][j][0]); /* tracer concentrations for each tracer species */
					for (int kz=swesoffset;kz<=levels;kz++)        /* initializing some variables */
					{
						qu_sl[kz] = qlat[kz] = etrdetail[kz] = a_etr[kz] = a_evar[kz] = interflow[kz] = lateralinflow[kz] = 0.0;
						gwex[kz]=gwin[kz]=0.0;
						theta_temp[kz]=theta_average[kz]=theta_sl_old[kz]=theta_sl[kz];
						macroinfiltration[kz]=0;
					}
					QDgridsub=QIgridsub=MACROsub=0.0;

					/* setting all saturated groundwater layers to the concentration of the first aquifer from the groundwater module */
					if (gwflow_active)
					{ 
						/* first groundwater layer is stored in GWLEVEL[i][j],
							the last water content of this layer is stored in GWtheta[i][j],
							the last groundwater table is stored in GWSTAND[i][j] (in meter below soil surface) */
						/* at first: the partly saturated uppermost groundwater layer */
						for (int z=1;z<=tracer_count;z++)
						{
							if (tracer_info[z].evaporating)
							{
								int k=(int)GWLEVEL[i][j];
								content_upper = conc_sl[z][k]*theta_sl[k]*dicken[k]; /* old tracer mass in the entire layer */ 
								content_lower = gw_conc[z][i][j]*thetasat[k]*(GWSTAND[i][j]-(hgeo_sl[k]-hgeo_sl[0])); /* new tracer mass in the lower part */
								/* subtracting the old mass in the saturated part */
								//				  if (k<levels)
								//					content_upper -= conc_sl[z][k+1]*thetasat[k+1]*(GWSTAND[i][j]-(hgeo_sl[k]-hgeo_sl[0]));
								//				  else
								content_upper -= content_lower;
								/* now adding the new mass in the lower part and calculating the new concentrations */
								conc_sl[z][k]=(content_upper+content_lower)/(theta_sl[k]*dicken[k]);
								/* now all the other layers */
								for (int k=(int)GWLEVEL[i][j]+1;k<=levels;k++) /* tracer concentrations of all completely saturated layers are set to the value from the first groundwater aquifer */
								{
									conc_sl[z][k]=gw_conc[z][i][j];
								}
							}
							else /* salts are handeld by their masses in the unsaturated zone, but by concentrations within groundwater, soa convertion has to be done here */
							{
								int k=(int)GWLEVEL[i][j];
								content_upper = (double)(conc_sl[z][k]); /* old salt mass in the entire layer times thickness */ 
								content_lower = (double)(gw_conc[z][i][j])*(double)1000.0*(double)thetasat[k]*((double)(GWSTAND[i][j])-((double)(hgeo_sl[k])-(double)(hgeo_sl[0]))); /* new salt mass in the lower part times saturated thickness in this layer */
								/* subtracting the old mass in the saturated part to get the real mass in the upper (unsaturated) part. 
				It is assumed, that the concentration in the next layer is the same like the concentration in 
				the lower (saturated) part of the actual layer*/
								//				  if (k<levels) 
								//					content_upper -= (double)(conc_sl[z][k+1])/(double)(dicken[k+1])*((double)(GWSTAND[i][j])-((double)(hgeo_sl[k])-(double)(hgeo_sl[0])));
								//				  else
								content_upper -= content_lower;
								/* now adding the new mass in the lower part and calculating the new concentrations */
								conc_sl[z][k]=(content_upper+content_lower);

								/* now all the other layers */
								for (int k=(int)GWLEVEL[i][j]+1;k<=levels;k++) /* tracer concentrations of all completely saturated layers are set to the value from the first groundwater aquifer */
								{
									conc_sl[z][k]=(double)(gw_conc[z][i][j])*(double)(theta_sl[k])*(double)(dicken[k])*1000.0;
								}
							}
						}
					}

					clay_layer=-9999;
					if (CLAY)
					if (CLAY[i][j]>0)
					{
						//            clay_layer = CLAY[i][j] > dicken[0]?
						//                         (int)floor(CLAY[i][j]/dicken[0]+2.4999): /* if the clay is below the 2nd layer */
						//                         (CLAY[i][j]<dicken[1]?                   /* else: is the clay within the first layer? */
						//                            2:                                    /* then set it to 2nd layer */
						//                            (int)floor((CLAY[i][j]-dicken[1])/dicken[2]+2.49999)); /* else: the  */
						clay_layer=0;
						accumThick=0;
						while ((accumThick<CLAY[i][j])&&(clay_layer<levels))
						{
							clay_layer++;
							accumThick+=dicken[clay_layer];
						}
					}

					// snow internal discharge is still valid in mm/interval
					// if layered snow model is active, the runoff components will be known later only - after Richards approach ran over the snow layers
					if ((SNOW)&&(!do_LayeredSnow))
					{
						if ((SNOW[i][j] > 10.0)&&(Psub>0.0))
						{
							SnowRunoffsub = Psub*SR_Anteil[TG]; // only the fraction "SR_Anteil" of snow melt is infiltrating into the soil
							Psub   -= SnowRunoffsub;            // rest is discharge
						}
						else SnowRunoffsub=0.0;
					}
					else SnowRunoffsub=0.0;

					if (MAXPOND && MAXPOND[i][j]>0) /* if there are dams around the field (ponds): mixing of tracers with water in ponds, evaporating water, let precipitation for infiltration be pond-storage */
					{
						for (int z=1;z<=tracer_count;z++)
						{ /* first level (pond) is in conc_sl[z][0] */
							if (tracer_info[z].radioactive) /* radioactive decay for radioactive tracers for the pond-storage */
							conc_sl[z][0]*=tracer_info[z].decfact;
							// since the surface runoff from silting up model has the same tracer concentrations like the 
							// precipitation and the snow outflow, we do not have to mix the tracers here, explicitely
							// but can use the precipitation concentration
							if (tracer_info[z].evaporating) /* tritium, 18O, ... */
								conc_sl[z][0]=mixing(conc_sl[z][0],tracer_module[z].cinput[i][j],POND[i][j],Psub+SnowRunoffsub+QDSUsub);
							else  /* salts like NaCl ... since INPUT is given in concentration ALWAYS, it must be converted to absolute masses now */
								conc_sl[z][0]+=tracer_module[z].cinput[i][j]*(Psub+SnowRunoffsub+QDSUsub);
						}
						POND[i][j] = POND[i][j] + SnowRunoffsub + QDSUsub; // let the SnowRunoffsub and the fraction of surface runoff from the silting up model flowing into the pond
						QDSUsub=0; // QDSU is stored in POND, again
						for (lay=startlayer;lay<=iNumVegetationLayers;lay++)
						{
							if (LUSE[lay][i][j]>nodata)
							{
								/* evaporating from pond, if there is enough content (otherwise the water would only infiltrate and extracted by roots) */
								if (ETPsub[lay]>=POND[i][j])  /* more ETP than is available in pond? (both storages are given in mm) */
								{
									/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
									ETRsub[lay] += POND[i][j]; /* adding the POND-content to the actual evaporation (there may be already some actual evap., because of snowevaporation) */
									ETPsub[lay] -= POND[i][j]; /* remaining evaporation is reduced by the POND-evaporation */
									POND[i][j] = 0.0;        /* POND evaporates complete */
								}
								else                       /* in POND is more water than necessary for evaporation */
								{
									/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
									ETRsub[lay]  += ETPsub[lay]; /* adding all of the potential to the actual evaporation */
									POND[i][j] -= ETPsub[lay]; /* evaporation is taken from the pond */
									ETPsub[lay]   = 0.0;       /* there is nothing left for further evaporation */
								}
							}
						}
						Psub    += POND[i][j]; /* maximum possible infitration amount equals rain plus content of ponding water storage */
						POND[i][j] = Psub;     /* available water for infiltration is the sum of precip. and POND-content */
						SnowRunoffsub=0.0;
					}
					else
					{ // no PONDs, but still QDSU from silting up and surface routing may have to be taken into account as surface storages
						for (int z=1;z<=tracer_count;z++)
						{
							if (tracer_info[z].radioactive) /* radioactive decay for radioactive tracers for the pond-storage */
								conc_sl[z][0]*=tracer_info[z].decfact;
							if (tracer_info[z].evaporating) /* tritium, 18O, ... */
								conc_sl[z][0]=tracer_module[z].cinput[i][j];
							else  /* salts like NaCl ... INPUT is given in concentration ALWAYS !! and must be converted to absolute masses now */
							{
								if (this->m_bAddQdsuToP) // if surface routing is active but silting up is not and no POND is active, 
								// then QDSU (which is Surface storage of surface routing module) is added to Psub
								{
									conc_sl[z][0]+=tracer_module[z].cinput[i][j]*(SnowRunoffsub+Psub+QDSUsub);
								}
								else
								{
									conc_sl[z][0]+=tracer_module[z].cinput[i][j]*(SnowRunoffsub+Psub);
								}
							}
						}
						if (this->m_bAddQdsuToP) // if surface routing is active but silting up is not and no POND is active, 
						// then evaporation should be taken from QDSU and remaing QDSU (which is Surface storage of surface routing module) is added to Psub
						{
							for (lay=startlayer;lay<=iNumVegetationLayers;lay++)
							{
								if (LUSE[lay][i][j]>nodata)
								{
									/* evaporating from pond, if there is enough content (otherwise the water would only infiltrate and extracted by roots) */
									if (ETPsub[lay]>=QDSUsub)  /* more ETP than is available in pond? (both storages are given in mm) */
									{
										/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
										ETRsub[lay] += QDSUsub; /* adding the QDSUsub-content to the actual evaporation (there may be already some actual evap., because of snowevaporation) */
										ETPsub[lay] -= QDSUsub; /* remaining evaporation is reduced by the QDSUsub-evaporation */
										QDSUsub = 0.0;        /* QDSUsub evaporates complete */
									}
									else                       /* in QDSUsub is more water than necessary for evaporation */
									{
										/* NO tracer-concentration change for evaporating tracers and ALSO NO change for MASSES of non evaporating tracers */
										ETRsub[lay]  += ETPsub[lay]; /* adding all of the potential to the actual evaporation */
										QDSUsub		 -= ETPsub[lay]; /* evaporation is taken from QDSUsub */
										ETPsub[lay]   = 0.0;       /* there is nothing left for further evaporation */
									}
								}
							}
							if (this->do_LayeredSnow==false || this->do_TMod==false)
							{   // if layered snow is true but some conventional melt approach is used, the snow melt is given as PREC to unsatzonmodel, so we need to add QDSU to P
								// only if heat transfer is used in unsatzonmodel (and energy balance approach is used in snowmodel together with layered snow model),
								// the melt is NOT provided in PREC, so we know that when either do_LayeredSnow is false or do_TMod 9is false, that we need to add QDSU to Psub now.
								// at the end of the split step, we will check again for this condition and add QDSU there 
								Psub += QDSUsub;
								QDSUsub = 0;
							}
						}
					}
					/* tracers/substance are considered in masses for non-evaporating tracers, not in concentrations, because evaporation may be lead to fallouts */
					SnowRunoffsub/=timestep; /* converting into m/s */
					//      infiltration(i,j); /* calculating infiltration and infiltration excess -> excess will become surface runoff */
					//      infiltration(i,j,clay_layer,gw_depth,s,dicken,hgeo_sl,theta_sl,Psub,INFEXsub);
					/* in infiltration() Psub is reduced by infiltration excess INFEXsub. Change in masses/concentrations is calculated later */
					/* from now on all (fluxes) are given in SI-measures what means:  "mm/interval" is converted into "m/s" a.s.o; */

					p     = Psub/timestep; /* converting P into m/s. p is now the total amount of rain infiltrating into the soil */

					// now calculating macropore infiltration. p is the infiltration amount from rain, snow melt or pond storage, so we can now
					// force as much water into the macropores as will be available from all possible sources
					// Psub was used for the last time here, so it's save to change p now.
					// Exception: If layered snow model is active, macropore infiltration is not possible as long as there si a snow cover
					if	(((lake_active==0)||
						  ((lake_active==1)&&(LAKE[i][j]==nodata))||
						  ((lake_active==1)&&(LAKE[i][j]!=nodata)&&(POND[i][j]<1.0))
						 ) && ((do_LayeredSnow==false)||((do_LayeredSnow==true) && (upperIndexSnow[i][j]>levelsSnow))) 
						)
					{
						MACROsub=CalculateMacroporeInfiltration(p,s,dicken,macrocapacity,thetasat,theta_sl,timestep,macroinfiltration);
					}
					else
					{
						MACROsub=0;
					}

					GWNsub = PERCOLsub = UPRISEsub = 0.0;
					maxlevel=1;
					maxrootdepth=0;
					for (lay=startlayer;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]>nodata)
						{
							etp   = ETPsub[lay]/timestep;
							vbg   = VBG[lay][i][j];
							wurz  = WURZ[lay][i][j];
							if (wurz>maxrootdepth) maxrootdepth=wurz;

							lastlayer=false;
							if (lay < iNumVegetationLayers)
							{
								if (LUSE[lay+1][i][j]==nodata) lastlayer=true;
							}
							else
							{
								lastlayer=true;
							}
							if (lastlayer==true)
							{
								evap = EVAPsub/timestep; // evaporation darf nur einmal entnommen werden!
							}
							else
							{
								evap=0;
							}
							for (int lc=0;lc<=ST[0].getLevels();lc++)
							{ 
								a_etr[lc]=0; /* set real transpiration to 0 -> must be done for each layer, to avoid using old content later in "while (a_etr[k]>0)" */
								a_evar[lc]=0; /* set real transpiration to 0 -> must be done for each layer, to avoid using old content later in "while (a_etr[k]>0)" */
							}

							//				/* to be removed
							//				maxlayer=1;
							maxlayer = evapotranspiration(i,j,LUtable[(int)LUSE[lay][i][j]].getFormPar(),// parameter for the densitiy-depth-function (1-depth/wurz)*(1+depth*formpar/wurz) 
															s, // soiltype for calculating theta for ETR reduction for each layer 
															LUtable[(int)LUSE[lay][i][j]].getReduParDry(),	
															LUtable[(int)LUSE[lay][i][j]].getReduParWet(), // maximum soilmoisture for potential transiration 
															LUtable[(int)LUSE[lay][i][j]].getLimitReduWet(), // minimum reduction the potential transpiration is reduced to when the soil is saturated 
															a_dicke, // return value: relativ fraction of all layers within root depth on root depth 
															a_etr,
															a_evar,  // return value: a_etr contains for each level the extraction by evapotranspiration m/s 
															wurz,vbg,etp,evap,dicken,thetasat,thetares,
															theta_sl,SE_sl,hgeo_sl,etrdetail,extractionweight,(int)GWLEVEL[i][j],DT,evapmaxdepth);
							//				*/
							/* now adding the fractions of evapotranspiration for calculating the total amount of evapotranspiration
			and calculating new theta-values */


							etr=0;
							evar=0;
							int klvl=1;
							int gwlvl=(int)GWLEVEL[i][j];
							while (klvl<=maxlayer)
							{
								if (klvl>GWLEVEL[i][j]) // extraction of remaining evaporation from soil
								{	// try to extract the water from the first saturated layer only (if theta will be too low, then take the next lower layer)
// until 20100605					double theta_grenz = ST[s].get_theta(gwlvl,gwlvl==1 ? (dicken[gwlvl]+dicken[gwlvl+1])*0.5+h_temp : dicken[gwlvl]+h_temp);
									// now (20100605) we try the same value as for all other theta-estimations
									double theta_grenz = ST[s].get_theta(gwlvl,3.45, thetares[gwlvl], thetasat[gwlvl]);
									// this is another older version
									// double theta_grenz = ST[s].get_theta(gwlvl,LUtable[(int)LUSE[lay][i][j]].getReduParDry());
									double theta_diff  = a_etr[klvl]*DT/dicken[gwlvl]; /* theta difference with regard to the level the gw is within (<_etr*DT/dicken[k]*(dicken[k]/dicken[gwlvl]) */
									double extrarate=0;
									while (theta_diff>0)
									{
										if (theta_sl[gwlvl]-theta_diff>=theta_grenz)
										{
											// when extracting the evaporation from the gwlv-layer, there is still enough water in it
											extrarate=theta_diff*dicken[gwlvl]/DT;
											a_etr[gwlvl]+=extrarate;
											a_etr[klvl]-=extrarate;
											theta_sl[gwlvl]-=theta_diff;
											theta_diff=0;
										}
										else if (theta_sl[gwlvl]>theta_grenz)
										{
											// there will not be enough water in the actual gwlv-layer after extracting the actual etr
											// but part of it could be drawn out of the actual level, the rest must be taken from below
											theta_diff=theta_sl[gwlvl]-theta_grenz;
											theta_sl[gwlvl]=theta_grenz;
											extrarate=theta_diff*dicken[gwlvl]/DT;
											a_etr[gwlvl]+=extrarate;
											a_etr[klvl]-=extrarate;
											gwlvl++;
											if (gwlvl > maxlayer)
											{
												gwlvl = maxlayer;
												break;
											}
											if ((dicken[gwlvl]>0)&&(dicken[gwlvl]!=dicken[gwlvl-1]))
											{	// if the layer thickness changes, the relative theta must be adjusted
												theta_diff=theta_diff*dicken[gwlvl-1]/dicken[gwlvl];
											}
										}
										else
										{ // the actual level is already too dry, take the water from the next level below the actual level, if possible
											gwlvl++;
											theta_diff=0;
											if (gwlvl > maxlayer)
											{
												gwlvl = maxlayer;
												break;
											}
										}
										etr+=extrarate;
#ifdef _WASIMDEBUG_
										for (int z=1;z<=tracer_count;z++)
										{
											if (tracer_info[z].evaporating)
											{
												balance[z] -= extrarate*timestep*conc_sl[z][gwlvl];
											}
										}
#endif
									}
								}
								else // normal extraction from all soil layers within the root depth
								{
									theta_sl[klvl] -= a_etr[klvl]*DT/dicken[klvl]; /* new theta values */
									if (theta_sl[klvl]<0)
									{ 
#ifdef _WASIMDEBUG_
#pragma omp critical

										(*WSstream) << "WARNING: bad theta is calculated after taking transpiration from row "<<i
										<<", col "<<j<<",layer "<<klvl<<", Theta-value: "<<theta_sl[klvl]<<"[end-entry]";
#endif
										if (klvl>1)
										{
											a_etr[klvl]+=(theta_sl[klvl]-thetares[klvl])*dicken[klvl]/DT;
											theta_sl[klvl]=thetares[klvl];
										}
										else
										{
											a_etr[klvl]+=theta_sl[klvl]*dicken[klvl]/DT;
											theta_sl[klvl]=0;
										}
									}
									etr+=a_etr[klvl];
#ifdef _WASIMDEBUG_
									for (int z=1;z<=tracer_count;z++)
									{
										if (tracer_info[z].evaporating)
										{
											balance[z] -= a_etr[klvl] * timestep*conc_sl[z][klvl];
										}
									}
#endif
								}

								klvl++;
							}

							klvl=1;
							while ((klvl<=maxlayer)&&(a_evar[klvl]>0))
							{
								// taking evar from the first one or more layers
								theta_sl[klvl] -= a_evar[klvl]*DT/dicken[klvl]; /* new theta values */
								evar += a_evar[klvl];

								if (theta_sl[klvl]<0)
								{ 
#ifdef _WASIMDEBUG_
#pragma omp critical
									(*WSstream) << "WARNING: attention: bad theta is calculated after taking evaporation from row " << i
									<<", col "<<j<<",layer "<<i<<", Theta-value "<<theta_sl[klvl]<<"[end-entry]";
									evar += theta_sl[klvl]*dicken[klvl]/DT; // reduce evar by the negative amount 
									theta_sl[klvl]=0; // set theta to 0 as absolute minimum for theta.
#endif
								}
								klvl++;
							}

							if (maxlevel<klvl-1)
							{
								maxlevel=klvl-1; /* to the level maxlevel the flux-calculation have to be carried out later,
								even if groundwater table is higher, auch wenn GW-Stand viel hoeher war,
								in order to get percolating water into the emty spaces within the layers the roots extracted water from */
							}
							ETRsub[lay] += (etr*timestep); /*  mm/h -> m/s, adding to ETR, bacause there may be already evaporation from snow */
							EVARsub+=(evar*timestep);
						}
						else
						{
							ETRsub[lay]=0;
						}

					}
					ETRsub[0]=0.0; // ETR and EVAR may contain some water already if snow evaporation was calculated (first, ETP if used, then, EVAP is used: consequently, ETR and EVAR are filled)
					
					for (lay=1;lay<=iNumVegetationLayers;lay++)
					{
						if (LUSE[lay][i][j]>nodata)
						{
							ETRsub[0]+=ETRsub[lay];
						}
					}
					ETRsub[0]+=EVARsub; // adding Evaporation to the total of the transpiration of all layers

					// ########################################### DELTAGW-handling now before calculating inflitration  - starts here ##################

					/* new groundwater discharge (baseflow) GWSTAND is negative! 
		this is done only, if no groundwater flow module is active, 
		otherwise baseflow is calculated within groundwater flow module */
					if ((!gwflow_active) && (!m_bCoupled)) 
					{
						DELTAGW[i][j]=Q0[TG]*exp((GWSTAND[i][j]*kelsqlat[TG]))/(timestep*substepcnt*splitfactor);
//						DELTAGW[i][j]=Q0[TG]*exp((GWSTAND[i][j]*kelsqlat[TG]))/(timestep);
						/* DELTAGW is always positive -> no rising groundwater in this branch */
					}
					else
					{
//						// if DELTAGW is the change in groundwater table (measured in flow perpendicular to the cells surface, i.e. in m/s
//						// we may have to either limit DELTAGW or to enlarge it depending on the change in porosity when a soil layer border is crossed by the rising or falling groundwater table
//						// Step 1: calculate potential DELTAGW based on HEAD1 (also imported from groundwater model) --> this is no more than a completely new calculation of DELTAGW, so it may even be removed from the groundwater model once it is working here
//						int k = (int)GWLEVEL[i][j];
//						double scalingfactor = 1.0; // how much the DELTAGW will be scaled depending on the effective differences in porosity between groundwater and soil paraneters for the affected soil layers
//						double temp=0; // the amount of water that fits into the filled up space by rising the water table or that will be released by dropping the water table to the level indicated by the value of HEAD1
//						if (HEAD1[i][j] > hgeo_sl[k])
//						{   // rising gw level
//							if (HEAD1[i][j]<hgeo_sl[0]) // if HEAD1 were above groudn, the full DELTAGW is used here in order to have exfiltration modeled
//							{
//								while ((HEAD1[i][j]>hgeo_sl[k - 1]) && (k > 1))
//								{
//									temp -= (thetasat[k] - theta_sl[k])*dicken[k] / DT; // negative DELTAGW for rising groundwater!
//									k--;
//								};
//								// now we are in the layer where HEAD1[i][j] is located in
//								if (k == GWLEVEL[i][j])
//								{ // no change of GW-level compared to the last time step
//									temp -= (thetasat[k] - theta_sl[k])*(HEAD1[i][j] - (hgeo_sl[0] + GWSTAND[i][j])) / DT;
//								}
//								else
//								{ // a new level was entered, so the complete filled volume counts
//									temp -= (thetasat[k] - theta_sl[k])*(HEAD1[i][j] - hgeo_sl[k]) / DT;
//								}
//
//								DELTAGW[i][j] = temp;
//							}
//						}
//						else
//						{   // falling gw level
//							while ((HEAD1[i][j]<hgeo_sl[k]) && (k<levels))
//							{
//								double templocal = ST[s].get_theta(k, 3.45);
//								if (k == GWLEVEL[i][j])
//								{
//									temp += (theta_sl[k] - templocal)*(hgeo_sl[0]+GWSTAND[i][j]-hgeo_sl[k]) / DT; // positive DELTAGW for falling groundwater! If we are in the layer of the groundwater table, only the saturated part drains
//								}
//								else
//								{
//									temp += (theta_sl[k] - templocal)*dicken[k] / DT; // positive DELTAGW for falling groundwater!
//								}
//								k++;
//							};
//							// last, incomplete layer
//							if (k < levels)
//							{
//								double templocal = ST[s].get_theta(k, 3.45);
//								if (k == GWLEVEL[i][j]) // gw level stays in the same layer as before
//								{
//									temp += (theta_sl[k] - templocal)*(hgeo_sl[0] + GWSTAND[i][j] - HEAD1[i][j]) / DT; // positive DELTAGW for falling groundwater! If we are in the layer of the groundwater table, only the saturated part drains
//								}
//								else
//								{
//									temp += (theta_sl[k] - templocal)*(hgeo_sl[k-1]-HEAD1[i][j]) / DT; // positive DELTAGW for falling groundwater!
//								}
//							}
//							DELTAGW[i][j] = temp;
//						}
////						if (DELTAGW[i][j] != 0) scalingfactor = temp / DELTAGW[i][j];
//						// that is the easiest way to synchronize gw model and unsaturated zone model
					}
					if (IRRIG) 
					if (IRRIGCODE[i][j]!=nodata)
					if (irrigationtable[(int)IRRIGCODE[i][j]].origin==1) /* if irrigationwater should be get from groundwater */
//					DELTAGW[i][j]+=(IRRIGsub/timestep); /* getting the irrigation water from groundwater */
					DELTAGW[i][j]+=(IRRIGsub/timestep); /* getting the irrigation water from groundwater */
					// if (!gwflow_active) DELTAGW[i][j] = 0; --> commented in version 10.06.04. This line was introduced my Mat when changing the temperature distribution. But it leads to zero base flow when using no groundwater
					// because IRRIG may change DELTAGW from negative to positve values and vice versa, the else branch must be changed into a regular piece of code
					//		else // if (geflow_active==true || m_bCoupled==true)
					//        {

					
					qlattemp=DELTAGW[i][j];  /* DELTAGW was calculated in the groundwater flow model in the previous time step */
					/* DELTAGW is now the balance of gw-flows. If this is negative (rising groundwater!) the layers are 
		now filled up and the baseflow calculation is not done later until irrigation loss makes the balance positive again */
					if (qlattemp<0.0) // rising groundwater table: DELTAGW is negative (upwards=negative, downwards=positive)
					{
						for (k = levels; k >= 1; --k) /* starting from the lowest layer -> no hanging water lenses are considered */
						{
							if (theta_sl[k] <= thetasat[k] * 0.9995) /* if saturation is above 99.95%, it is assumed, that the layer is within groundwater */
							{
								break;
							}
						}
						// updating the actual theta for the level, the groundwater surface is loacted actually
						GWtheta[i][j] = theta_sl[k];

						/* starting with layer k the upper layers are filled with lateral inflows
			(tracer concentrations are updated later in new_concentrations) */
						do
						{
							if (k>0)
							{
								maxlateralflow = (thetasat[k]- theta_sl[k])*dicken[k]/DT; /* what amount of water fits into the layer? positive sign! */
								if (maxlateralflow <= -qlattemp) /* DELTAGW is negative! */
								{ /* layer can absorb less than DELTAGW gives */
									lateralinflow[k]=maxlateralflow; 
									qlattemp+=maxlateralflow; /* DELTAGW neg.sign, maxlatflow pos. sign -> DELTAGW is decreased */
									theta_sl[k]=thetasat[k];
									k--; /* next layer is affected */
									// updating the actual theta for the level, the groundwater surface is loacted actually
									if (k==0)
										GWtheta[i][j] = thetasat[1];
									else
										GWtheta[i][j] = theta_sl[k];
								}
								else /* maxflow greater or equal DELTAGW */
								{ /* inflow only in one layer */
									lateralinflow[k]=-qlattemp; /* all inflow is lateral flow for the actual layer */
									qlattemp=0.0;
									theta_sl[k]+=lateralinflow[k]*DT/dicken[k];
									GWtheta[i][j] = theta_sl[k];
								}
							}
							else
							{ /**/
								lateralinflow[0]=-qlattemp; /* all inflow is lateral flow for the actual layer */
								qlattemp=0.0;
							}
						} 
						while (qlattemp<0.0);
					}
					//        } 

					for (k = levels; k >= 1; k--) /* starting from the lowest layer -> no hanging water lenses are considered */
					{
						if (theta_sl[k] <= thetasat[k] * double(0.9995)) /* if saturation is above 99.95%, it is assumed, that the layer is within groundwater */
						{
							break;
						}
					}
					if (!gwflow_active) k = (int)GWLEVEL[i][j];
					/* now substracting the groundwater flow from thetas */
					if (k==0) k=1;
					GWLEVEL[i][j]=start_level=k;
					// updating the actual theta for the level, the groundwater surface is loacted actually
					GWtheta[i][j] = theta_sl[k];
					if (DELTAGW[i][j] >= 0) /* if groundwater flow balance is positive (or zero) for baseflow (falling or constant groundwater table -> counts positive downwards) */
					{
						qlattemp=DELTAGW[i][j];
						do
						{
							theta_temp[k]=ST[s].get_theta(k, 3.45, thetares[k], thetasat[k]);
							//				theta_temp[k]=ST[s].get_theta(k, // 3.45,  a constant suction of 3.45 m is assumed to refer to field capacity 
							//				  k==1 ? (dicken[k]+dicken[k+1])*0.5 : dicken[k]); /* getting the theta after gravitative drainage */
							/* if the actual soil moisture is less than the minimum moisture required for downward flow,
										the baseflow is taken from the layer below the actual layer
			(this is mostly an initialisation problem) */
							if (theta_sl[k]>theta_temp[k])
							{
								if (theta_sl[k]-qlattemp*DT/dicken[k]<theta_temp[k]) /* is more flux calculated, than the actual layer contains */
								{ /* baseflow is extracted also from the layer below the actual layer a.s.o. */
									qlat[k]=(theta_sl[k]-theta_temp[k])*dicken[k]/DT;
									qlattemp-=qlat[k];
									theta_sl[k]=theta_temp[k];
									k++;
									// updating the actual theta for the level, the groundwater surface is loacted actually
									GWtheta[i][j] = theta_sl[k];
									start_level++;
									GWLEVEL[i][j]+=1.0;
								}
								else /* groundwater table remains in the actual level */
								{
									theta_sl[k]-=qlattemp*DT/dicken[k];
									// updating the actual theta for the level, the groundwater surface is loacted actually
									GWtheta[i][j] = theta_sl[k];
									qlat[k] = qlattemp;
									qlattemp=0;
									k=levels+1; /* break criteria */
								}
							}
							else
							{
								k++;
								// updating the actual theta for the level, the groundwater surface is loacted actually
								GWtheta[i][j] = theta_sl[k];
								start_level++;
								GWLEVEL[i][j]+=1.0;
							}
						}
						while (k<=levels);
						DELTAGW[i][j]-=qlattemp; /* if not all baseflow could be extracted (because soil is empty)
									the remaining part has to be substracted from DELTAGW */
						if (start_level>levels) start_level=levels;
						if (GWLEVEL[i][j]>levels) GWLEVEL[i][j]=levels;
					}

					// ########################################### DELTAGW-handling now before calculating inflitration  - ends here ##################


					/* ##################################################################################################################### */
					/* ###################################  This is the core of the unsaturated zone flow modeling ######################### */
					/* ##################################################################################################################### */

					// test, if moving these cells to the rank 0 works. This section produces compete nonsense and must be deleted for production mode!!
//					if (do_ArcVeg && splitcnt==1)
//					{
//						for (int k = 1; k <= levels; k++)
//						{
//							if (ExcessIceFactor[i][j][k] != nodata)
//							{
//					//			ExcessIceFactor[i][j][k] += (double(rand()) / double(RAND_MAX)) / 10.0;
//							}
//						}
//					}

					// do the temperature model every time step.
					// heat transfer modelling is done after water flux calculation because of heat transfer by advection does need the water fluxed
					if (this->do_TMod > 0)
					{
						for (int heatmultiple_counter=0;heatmultiple_counter<heatmultiple;heatmultiple_counter++)
						{
							ST[s].CalcSEanddEOverdT(0,SE_sl, dEOverdT_sl, temper_sl, this->dTModFactor, theta_sl, thetares,thetasat); // calculating new relative saturation values (water=1...ice=0)
							if (do_LayeredSnow) LayeredSnowProperties[1].CalcSEanddEOverdT(-levelsSnow,SE_sl, dEOverdT_sl, temper_sl, this->dTModFactor, theta_sl, thetares, thetasat); // calculating new relative saturation values (water=1...ice=0)
							CalculateLambdaEff(swesoffset,lambda_sl, s, SE_sl, theta_sl, thetasat,this->dTModLamWat, this->dTModLamIce); // effective thermal conductivity
							// dEOverdT is yet only a scaling factor. Must be multiplied by Rho*LatHeat*(theta-theta_res)
							Expand_dEoverdT(swesoffset,theta_sl,thetares,dEOverdT_sl);
							// calculating actual heat capacity
							ST[s].CalcHeatCapacity(0,temper_sl,heatcap_sl,SE_sl,dEOverdT_sl,theta_sl,this->dTModRhoWater,this->dTModCWater,this->dTModCIce0,this->dTModCIce20);
							if (do_LayeredSnow) LayeredSnowProperties[1].CalcHeatCapacity(-levelsSnow,temper_sl,heatcap_sl,SE_sl,dEOverdT_sl,theta_sl,this->dTModRhoWater,this->dTModCWater,this->dTModCIce0,this->dTModCIce20);
							int use = (int)LUSE[1][i][j];
							actNFactor = LUtable[use].GetNfactor(WSstream,tag,monat,jahr);
							if (actNFactor<-9999)
							{
								if (T[i][j] < 0) 
								{ 
									actNFactor = LUtable[use].GetNFactorFreezing(WSstream);
									if (actNFactor<-9999)
									{
										actNFactor = dTModNFactorFreezing;
									}
								} 
								else 
								{
									actNFactor = LUtable[use].GetNFactorThawing(WSstream);
									if (actNFactor<-9999)
									{
										actNFactor = dTModNFactorThawing;
									}
								}
							}
							Tsurf = T[i][j];
							if (EBAL && do_LayeredSnow)
							{
								if (upperIndexSnow[i][j]<=levelsSnow) // change Tsurf only if there is some snow on the ground
								{
									if ((EBAL[i][j]==0)&&(SURFTEMP[i][j]<273.16))
									{	// cooling conditions
										Tsurf = SURFTEMP[i][j]-273.16;
										if (Tsurf<T[i][j]-13.0)
										{
											Tsurf=T[i][j]-13.0;
										}
										actNFactor = dTModNFactorFreezing;
									}	else if (EBAL[i][j]>0)	{ 
										// melting conditions: EBAL > 0
										// mixing temperature of the uppermost two layers
										double Tmix;
										double thicksum;
										if (swesoffset<-1)
										{
											// at least two snow layers: in case the upper layer gets too thin, we look at the two uppermost layers to avoid unrealistic surface temperatures
											double swel1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*SE_sl[1+swesoffset]; //swes[i][j][1+swesoffset+levelsSnow];
											double swes1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*(1-SE_sl[1+swesoffset]); //swel[i][j][1+swesoffset+levelsSnow];
											double swel2 = 1000*dicken[2+swesoffset]*theta_sl[2+swesoffset]*SE_sl[2+swesoffset]; //swes[i][j][2+swesoffset+levelsSnow];
											double swes2 = 1000*dicken[2+swesoffset]*theta_sl[2+swesoffset]*(1-SE_sl[2+swesoffset]); //swel[i][j][2+swesoffset+levelsSnow];

											Tmix = ((swes1*dTModCIce0+swel1*dTModCWater) * temper_sl[1+swesoffset] 
												  + (swes2*dTModCIce0+swel2*dTModCWater) * temper_sl[2+swesoffset]) / 
												   ((swes1+swes2)*dTModCIce0 + (swel1+swel2)*dTModCWater);
											thicksum=dicken[1+swesoffset]+dicken[2+swesoffset];
										}
										else
										{
											// only one snow layer, so we may use that layers properties or mix it with the first soil layer or we reduce the energy balance to an amount that is sufficient to melt the actual content of the last snow layer only
											double swes1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*SE_sl[1+swesoffset]; //swes[i][j][1+swesoffset+levelsSnow];
											double swel1 = 1000*dicken[1+swesoffset]*theta_sl[1+swesoffset]*(1-SE_sl[1+swesoffset]); //swel[i][j][1+swesoffset+levelsSnow];
											double requiredMeltEnergy = (337500.0*swes1 - (temper_sl[0] * (2100.0 * swes1 + 4187.0 * swel1)))/(InV*60.0);
											if (EBAL[i][j]>requiredMeltEnergy) EBAL[i][j]=requiredMeltEnergy;
											Tmix = temper_sl[0];
											thicksum = dicken[0];
										}
										Tsurf = Tmix + EBAL[i][j]*thicksum/(2.0*lambda_sl[1+swesoffset]);
										if (Tsurf<0)
										{
											Tsurf=0;
										} else if (Tsurf>T[i][j]+13) {
											Tsurf = T[i][j]+13;
										}
//										Tsurf = temper_sl[1+swesoffset]+EBAL[i][j]*dicken[1+swesoffset]/(2.0*lambda_sl[1+swesoffset]);
										actNFactor = dTModNFactorThawing;
									} else {
										Tsurf = temper_sl[1+swesoffset]; // sublimation conditions, i.e. all energy is required for snow evaporation, no heat flux into snow possible
										actNFactor = 1.0;
									}

									// actNFactor=1;
									if (dumpmode)
									{
										dump.Tsurf=Tsurf;
									}
								}
							}
							if (iTModMethod>10)
							{
								SolveTriDiag(swesoffset,i,j,lambda_sl,heatcap_sl,dicken,temper_sl,timestep_heat,coeff_a,coeff_b,coeff_c,
									         right_v,result_x,actNFactor,macroinfiltration,qu_sl,Tsurf);
								if (this->do_TMod > 0)
								{
									for (int z=1+swesoffset;z<=levels;z++)
									{
									//	if (z==0)
									//	{
									//		if ((result_x[z]<-30)||(result_x[z]>30)||printflag>0)
									//		{
									//			cout << "temperatue["<<i<<"]["<<j<<"]["<<z<<"], step "<<splitcnt<<": " << result_x[z] << " | old T[-1]: "<< temper_sl[-1] <<" | T[0]: "<< temper_sl[0] << " | T[1]" << temper_sl[1] << " Tsurf: "<< Tsurf << " swesoffset: " << swesoffset << " lambda_sl[0]: "<< lambda_sl[0] <<" heatcap_sl: "<< heatcap_sl[0] << " dicken: "<< dicken[0] << endl;
									//			printflag=1;
									//		}
									//	}
										temper_sl[z]=temper_sl_new[z]=result_x[z];
									}
								}
							}
							else
							{
								// calculating heat fluxes
								CalcHeatFluxes(swesoffset,qHeat_sl,lambda_sl,dicken,temper_sl,TLOWBND[i][j],T[i][j],actNFactor,Tsurf); // upper boundary condition: air temperature, lower boundary condition: given in boundary condition grid
								// calculating temperature change
								CalcDeltaT(swesoffset,deltaT_sl, qHeat_sl,theta_sl,dicken,heatcap_sl,s,(int)timestep_heat);
								// calculating new temperatures
								CalcNewTemperature(swesoffset,temper_sl_new,temper_sl,deltaT_sl,qHeat_sl,TLOWBND[i][j],T[i][j],true,Tsurf);
								if (this->do_TMod > 0)
								{
									for (int z=1+swesoffset;z<=levels;z++)
									{
										temper_sl[z]=temper_sl_new[z];
									}
								}
							}
//							(*WSstream) << tag <<"."<< monat << " "<< stunde << " T: " << temper_sl[1] << " C: " << heatcap_sl[1] << "[end-entry]";
						}
					}

					if ((do_LayeredSnow))
					{
						double *swes_sl = swes[i][j];
						double *swel_sl = swel[i][j];
						double *densSnow_sl = densSnow[i][j];
						if (do_TMod)
						{
							// update of the swes and swel stacks 
							extraSnowRunoffsub = updateSnowStorageLayers(swes_sl,swel_sl,theta_sl,thetasat,dicken,SE_sl,densSnow_sl,temper_sl);

							// press extraSnowRunoffsub into first soil layer, only the remaining amount will be real extraSnowRunoff
							double freeVolumeLayer1 =  (thetasat[1]-theta_sl[1])*dicken[1]*1000.0;
							if (extraSnowRunoffsub<freeVolumeLayer1)
							{
								theta_sl[1]+=(extraSnowRunoffsub/(dicken[1]*1000.0));
								extraSnowRunoffsub=0;
							}
							else
							{
								theta_sl[1]=thetasat[1];
								extraSnowRunoffsub-=freeVolumeLayer1;
								if (extraSnowRunoffsub<0) extraSnowRunoffsub=0;
								double freeVolumeLayer2 =  (thetasat[2]-theta_sl[2])*dicken[2]*1000.0;
								if (extraSnowRunoffsub<freeVolumeLayer2)
								{
									theta_sl[2]+=(extraSnowRunoffsub/(dicken[2]*1000.0));
									extraSnowRunoffsub=0;
								}
								else
								{
									theta_sl[2]=thetasat[2];
									extraSnowRunoffsub-=freeVolumeLayer2;
									if (extraSnowRunoffsub<0) extraSnowRunoffsub=0;
							
								}
							}
						}
						swes_sl[0]=0;
						swel_sl[0]=0;
						// new swesoffset should be computed here after updateing the snow layers
						swesoffset=0;
						while ((swes[i][j][levelsSnow+swesoffset]>0)&&(-swesoffset<levelsSnow))
						{ 
							swesoffset--;
						}
						upperIndexSnow[i][j]=levelsSnow+1+swesoffset;
					}

					/* now calculating fluxes iterative (method ==2) or single step, but also with iterative gradient restrictions (mehtod==1) */
					for (int kz = swesoffset; kz <= levels; kz++)
					{
						theta_average[kz] = theta_temp[kz] = theta_sl_old[kz] = theta_sl[kz];
					}

					ziter=0;
					weight=0; /* weight: contains the sum of 1..ziter for weighting moving averages of theta */
//					double qu0 = 0; // as temporary memory of what value qu_sl[0] had after normal flux calculation
					double max_Err=max_Error;
					do /* iterate the finite differences scheme */
					{
						// if splitfactor is 1, then the result was already calculated in the precondition calculation for the split interval 
						// (with the exception of not capturing altered conditions by infiltrating water etc, which usually should also led to a increase in sub time steps)
						// if ((this->do_TMod > 0)&&(splitfactor>1)) 

						hydraulic_heads_iterative(swesoffset,i,j,ziter,maxlevel,s,h_temp,gwin,gwex,
						                          h_temp1,theta_sl_old,theta_average,hh_sl,hgeo_sl,SE_sl,thetares,thetasat);
						if (QDRAIN)
						{
							specific_fluxes(swesoffset,i,j,s,maxlevel,draindist,clay_layer,p,
							                interflow,qu_sl,gwin,gwex,conc_gwin,
							                thetares,thetasat,theta_sl_old,theta_average,
							                hgeo_sl,hh_sl,dicken,QDRAINsub,GWINsub,GWEXsub,DT,do_TMod,SE_sl);  // flux-calculations including gradient restrictions 
						}
						else
						{
							double dummy;
								specific_fluxes(swesoffset,i,j,s,maxlevel,draindist,clay_layer,p,
							                interflow,qu_sl,gwin,gwex,conc_gwin,
											thetares,thetasat,theta_sl_old,theta_average,
							                hgeo_sl,hh_sl,dicken,dummy,GWINsub,GWEXsub,DT,do_TMod,SE_sl);  // flux-calculations including gradient restrictions 
						}
		
						// extraSnowRunoffsub should now be added to qu_sl[0] or provided as an extra parameter to flux_restrictions
//						qu0 = qu_sl[0];
//						if (extraSnowRunoffsub>0)
//						{
//							qu_sl[0] -= (extraSnowRunoffsub/(DT*1000)); // try to infiltrate all snow melt into the first soil layer
//						}
						if ((levels >= 200) || do_LayeredSnow)
						{
							flux_restrictions(swesoffset, i, j, 10, maxlevel, gwin, gwex, theta_sl_old, theta_temp,
								thetasat, thetares, qu_sl, interflow, dicken, GWINsub, GWEXsub, DT, SE_sl);
						}
						else
						{ // this version is faster, because internally, fixed reserved buffer is used. However, there has to be a limited number of soil layers, which is defined to be 200. 
							// for soils with more levels, the traditional method will be used
							flux_restrictions_opt(swesoffset, i, j, 10, maxlevel, gwin, gwex, theta_sl_old, theta_temp,
								thetasat, thetares, qu_sl, interflow, dicken, GWINsub, GWEXsub, DT, SE_sl);
						}

						ziter++;
						if (ziter==21) max_Err=max_Error*5.0; /* avoiding to long iteration cycles by setting the maximum allowed error in theta to a higher value */
						m_err=new_thetas(swesoffset,i,j,ziter,maxlevel,gwin,gwex,theta_sl,theta_average,theta_sl_old,qu_sl,interflow,dicken,weight,DT);


#ifdef _WASIMDEBUG_
						/* output of the results of each iteration step is only done if _DEBUG is defined (in MS Visual C++ it is the standard setting for debugging) */
						/*          if ((i==control_row)&&(j==control_col))
			fprintf(tempdatei,"%2d\t%2d\t%2d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",monat,tag,stunde,theta_sl[1],
				theta_sl[2],theta_sl[3],qu_sl[1],interflow[1],qu_sl[2],interflow[2],qu_sl[3],interflow[3],hh_sl[1],hh_sl[2],hh_sl[3],hh_sl[1]-hh_sl[2],hh_sl[2]-hh_sl[3],theta_average[1],theta_average[2],theta_average[3]);
*/
#endif

					} while ((m_err>max_Err)&&(ziter<max_ziter));

					if ((do_LayeredSnow)&&(do_TMod))
					{
						double *swes_sl = swes[i][j];
						double *swel_sl = swel[i][j];
						double *densSnow_sl = densSnow[i][j];
						// update of the swes and swel stacks 
						extraSnowRunoffsub += updateSnowStorageLayers(swes_sl,swel_sl,theta_sl,thetasat,dicken,SE_sl,densSnow_sl,temper_sl);

						swes_sl[0]=0;
						swel_sl[0]=0;
						// new swesoffset should be computed here after updateing the snow layers
						swesoffset=0;
						while ((swes[i][j][levelsSnow+swesoffset]>0)&&(-swesoffset<levelsSnow))
						{ 
							swesoffset--;
						}
						upperIndexSnow[i][j]=levelsSnow+1+swesoffset;
					}


					if (extraSnowRunoffsub>0)
					{
//						extraSnowRunoffsub = extraSnowRunoffsub + (qu_sl[0]-qu0)*(DT*1000.0);
//						if (extraSnowRunoffsub<0) extraSnowRunoffsub=0;
						extraSnowRunoff += extraSnowRunoffsub; // extraSnowRunoffsub may have beenreduced in specific_fluxes() or flux_restrictions(). So the absolute (or corrected) amount can be added to extraSnowRunoff only after the soil water fluxes were computed
					}

					/* ##################################################################################################################### */
					/* ##################################  This is the end of the core of unsat. zone flow modeling ######################## */
					/* ##################################################################################################################### */

					if (this->dumpmode==true)
					{
						bool writedumptheta=false;
						bool writedumpheat=false;
						bool writedumpdeltagw=false;
						bool writedumpthetasnow=false;

						int dumpthetalayer = 0;
						int dumpheatlayer = 0;
						int dumpthetasnowlayer=-9999;

						if (DELTAGW[i][j]<-10 || DELTAGW[i][j]>10) 
							writedumpdeltagw=true;
						for (int kk = 1-levelsSnow; kk <= levels; kk++)
						{
							if (((theta_sl[kk]<0)&&!(theta_sl[kk]<=-9998.9)) || (theta_sl[kk]>thetasat[kk] + 0.001))
							{
								writedumptheta = true;
								dumpthetalayer = kk;
								break;
							}
						}

						for (int kk = 1-levelsSnow+1; kk <= 0; kk++)
						{
							if ((theta_sl[kk-1]>0)&&(theta_sl[kk]<=-9998.9))
							{
								writedumpthetasnow = true;
								dumpthetasnowlayer = kk;
								break;
							}
						}

						for (int kk = 1-levelsSnow; kk <= levels; kk++)
						{

							if (do_TMod && temper_sl[kk]>-9998.9) // only for valid layers
							{
								if (kk==1+swesoffset)
								{
									if ((temper_sl[kk] < -70.0) || (temper_sl[kk] > +70.0) ||
										
										(temper_sl[kk] > dump.temper_sl[kk + levelsSnow] + 0.001) && // temperature increases
										(temper_sl[kk] > Tsurf*actNFactor + 0.001) && // but not more than the upper boundary condition
										(temper_sl[kk] > dump.temper_sl[kk + 1+levelsSnow] + 0.001) && (temper_sl[kk] > temper_sl[kk + 1] + 0.001) // and also not more than the higher value of either the start or the end temperature of the layer below
										||
										(temper_sl[kk] < dump.temper_sl[kk + levelsSnow] - 0.001) && // temperature decreases
										(temper_sl[kk] < Tsurf*actNFactor - 0.001) && // but not more than the upper boundary condition
										(temper_sl[kk] < dump.temper_sl[kk + 1+levelsSnow] - 0.001) && (temper_sl[kk] < temper_sl[kk + 1] - 0.001)) // and also not omre than the lower value of either the start orthe end value of the layer below 
									{
										writedumpheat = true;
										dumpheatlayer = kk;
										break;
									}
								}
								else
								{
									if ((temper_sl[kk] < -70.0) || (temper_sl[kk] > +70.0) ||

										(temper_sl[kk] > dump.temper_sl[kk + levelsSnow] + 0.001) && // temperature rises
										(temper_sl[kk] > dump.temper_sl[kk - 1+levelsSnow] + 0.001) && (temper_sl[kk] > temper_sl[kk - 1] + 0.001) && // but should neither rise above the old nor above the new value of the layer above
										(temper_sl[kk] > dump.temper_sl[kk + 1+levelsSnow] + 0.001) && (temper_sl[kk] > temper_sl[kk + 1] + 0.001)  // and should also neither rise above the old nor above the new value of the layer below
										||
										(temper_sl[kk] < dump.temper_sl[kk + levelsSnow] - 0.001) && // temperature decreses
										(temper_sl[kk] < dump.temper_sl[kk - 1+levelsSnow] - 0.001) && (temper_sl[kk] < temper_sl[kk - 1] - 0.001) && // but should neither sink below the old nor below the new value of the layer above
										(temper_sl[kk] < dump.temper_sl[kk + 1+levelsSnow] - 0.001) && (temper_sl[kk] < temper_sl[kk + 1] - 0.001)) // and should also neither sink below the old nor below the new value of the layer below
									{
										writedumpheat = true;
										dumpheatlayer = kk;
										break;
									}
								}
							}
						}
	//					writedumpdeltagw=writedumpheat=writedumptheta =true;
						if (writedumpheat || writedumptheta || writedumpthetasnow || writedumpdeltagw)
						{
							if ((actDumpFiles < maxDumpFiles) || maxDumpFiles==0)
							{
								std::ofstream dumpstream;
								std::stringstream fname;
								fname << "MPI_" << mpi_rank <<"_";
								if (writedumpheat)
								{
									fname << "heat(" << dumpheatlayer << ")_";
								}
								if (writedumptheta)
								{
									fname << "theta(" << dumpthetalayer << ")_";
								}
								if (writedumpdeltagw)
								{
									fname << "dgw_";
								}
								if (writedumpthetasnow)
								{
									fname << "thetasnow(" << dumpthetasnowlayer << ")_";
								}

								fname << "dump_row_" << i << "_col_" << j << "_substep_" << splitcnt << "_" << jr << setw(2) << setfill('0') << mt << setw(2) << setfill('0') << tg << setw(2) << setfill('0') << se << ".dmp";
								std::string sname(fname.str());
								dumpstream.open(sname.c_str(), std::ios::binary);
								dumpstream.write((char *)(&dump.maxlayer), sizeof(int));
								dumpstream.write((char *)(&dump.h_temp), sizeof(double));
								dumpstream.write((char *)(&dump.p), sizeof(double));
								dumpstream.write((char *)(&dump.DELTAGW), sizeof(double));
								dumpstream.write((char *)(&dump.ETR[0]), sizeof(double)*(iNumVegetationLayers + 1));
								dumpstream.write((char *)(&dump.gwin[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.gwex[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.h_temp1[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.theta_sl_old[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.theta_average[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.hh_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.hgeo_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.qu_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.interflow[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.thetasat[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.dicken[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.theta_temp[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.thetares[0]), sizeof(double)*(levels + 1+levelsSnow));
								dumpstream.write((char *)(&dump.theta_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
								if (do_TMod) {
									dumpstream.write((char *)(&dump.SE_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.temper_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.lambda_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.dEOverdT_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.heatcap_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.qHeat_sl[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.temper_sl_new[0]), sizeof(double)*(levels + 1+levelsSnow));
									dumpstream.write((char *)(&dump.timestep_heat), sizeof(double));
								}
								if (do_LayeredSnow)
								{
									dumpstream.write((char *)(&dump.swes[0]),		sizeof(double)*(levelsSnow+1));
									dumpstream.write((char *)(&dump.swel[0]),		sizeof(double)*(levelsSnow+1));
									dumpstream.write((char *)(&dump.densSnow[0]),	sizeof(double)*(levelsSnow+1));
									dumpstream.write((char *)(&dump.thickSnow[0]),	sizeof(double)*(levelsSnow+1));
									dumpstream.write((char *)(&dump.upperIndexSnow),sizeof(double));
									dumpstream.write((char *)(&dump.SSNO),			sizeof(double));
									dumpstream.write((char *)(&dump.SLIQ),			sizeof(double));
									dumpstream.write((char *)(&dump.SSTO),			sizeof(double));
								}
								dumpstream.write((char *)(&dump.s), sizeof(int));
								for (int y=1;y<=iNumVegetationLayers;y++) 
									dumpstream.write((char *)(&(LUSE[y][i][j])),sizeof(double));
								dumpstream.write((char *)(&dump.stunde), sizeof(int));
								dumpstream.write((char *)(&dump.tag), sizeof(int));
								dumpstream.write((char *)(&dump.monat), sizeof(int));
								dumpstream.write((char *)(&dump.jahr), sizeof(int));
								dumpstream.write((char *)(&dump.splitfactor), sizeof(double));
								dumpstream.write((char *)(&dump.DT), sizeof(double));
								dumpstream.write((char *)(&dump.swesoffset), sizeof(int));
								if (do_TMod) {
									dumpstream.write((char *)(&dump.splitfactor_heat), sizeof(double));
									dumpstream.write((char *)(&dump.heatmultiple), sizeof(int));
									dumpstream.write((char *)(&dump.T), sizeof(double));
									dumpstream.write((char *)(&Tsurf),sizeof(double));
								}
								dumpstream.close();

								if ((do_LayeredSnow)||(do_TMod))
								{
									std::ofstream dumpstream2;
									std::stringstream fname2;
									fname2 << "MPI_" << mpi_rank <<"_";
									if (writedumpheat)
									{
										fname2 << "heat(" << dumpheatlayer << ")_";
									}
									if (writedumptheta)
									{
										fname2 << "theta(" << dumpthetalayer << ")_";
									}
									if (writedumpdeltagw)
									{
										fname2 << "dgw_";
									}
									if (writedumpthetasnow)
									{
										fname2 << "thetasnow(" << dumpthetasnowlayer << ")_";
									}
	
									fname2 << "dump_row_" << i << "_col_" << j << "_substep_" << splitcnt << "_" << jr << setw(2) << setfill('0') << mt << setw(2) << setfill('0') << tg << setw(2) << setfill('0') << se << ".dmp_asc";
									std::string sname2(fname2.str());
									dumpstream2.open(sname2.c_str());
	
									dumpstream2 << "this file contains some values in human readabel format for layered snow model and heat transfer model" << endl;
	
									if (do_LayeredSnow)
									{
										dumpstream2 << "     swes ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << swes[i][j][z] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.swes ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << dump.swes[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     swel ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << swel[i][j][z] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.swel ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << dump.swel[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     thickSnow ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << thickSnow[i][j][z] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.thickSnow ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << dump.thickSnow[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     densSnow ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << densSnow[i][j][z] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.densSnow ";
										for (int z=1;z<=levelsSnow;z++) dumpstream2 << dump.densSnow[z] << "\t";
										dumpstream2 << endl;
									}
									if (do_TMod)
									{
										dumpstream2 << "     theta_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << theta_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.theta_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.theta_sl[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     SE_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << SE_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.SE_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.SE_sl[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     temper_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << temper_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.temper_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.temper_sl[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     lambda_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << lambda_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.lambda_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.lambda_sl[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     dEOverdT_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dEOverdT_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.dEOverdT_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.dEOverdT_sl[z] << "\t";
										dumpstream2 << endl;
	
										dumpstream2 << "     heatcap_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << heatcap_sl[z-levelsSnow] << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.heatcap_sl ";
										for (int z=1;z<=levels+levelsSnow;z++) dumpstream2 << dump.heatcap_sl[z] << "\t";
										dumpstream2 << endl;

										dumpstream2 << "    Tsurf ";
										dumpstream2 << Tsurf << "\t";
										dumpstream2 << endl;
										dumpstream2 << "dump.Tsurf ";
										dumpstream2 << dump.Tsurf << "\t";
										dumpstream2 << endl;

									}
									dumpstream2.close();
								}
							}
							actDumpFiles++;
							if (closeAfterDump) 	{
									(*WSstream) << "FATAL: \n\n -------- command line option -closeAfterDump active! ---------------------------------------------\n\n"
									<< " a dumpfile has been written and WaSiM aborts now.\n"
									<<      "-------------------------------------------------------------\n\n" <<"[end-entry]";

									throw(199);
							} else if (repairAfterDump) 	{

								if (writedumpdeltagw)
								{
									if (DELTAGW[i][j] < -10)		{
										DELTAGW[i][j] = -10.0;
									} else if (DELTAGW[i][j] > -10) {
										DELTAGW[i][j] = 10.0;
									}
								}

								for (int kk = 1; kk <= levels; kk++)
								{

									if (writedumptheta)
									{
										if (theta_sl[kk] < 0)		{
											theta_sl[kk] = 0.0;
										}
										else if (theta_sl[kk]>thetasat[kk] + 0.001) {
											theta_sl[kk] = thetasat[kk];
										}
									}

									if (writedumpheat)
									{
										if (kk == 1)
										{
											if (temper_sl[kk] < -70.0) {
												temper_sl[kk] = -70.0;
											}
											else if (temper_sl[kk] > +70.0) {
												temper_sl[kk] = 70.0;
											}
											else if ((temper_sl[kk] > dump.temper_sl[kk] + 0.001) && // temperature increases
												(temper_sl[kk] > T[i][j] * actNFactor + 0.001) && // but not more than the upper boundary condition
												(temper_sl[kk] > dump.temper_sl[kk + 1] + 0.001) && (temper_sl[kk] > temper_sl[kk + 1] + 0.001)) // and also not more than the higher value of either the start or the end temperature of the layer below
											{
												double maxt = dump.temper_sl[kk]; // maximum of all temperatures is set to old temperature of that layer
												if ((T[i][j] * actNFactor) > maxt) maxt = (T[i][j] * actNFactor); // if upper boundary condition is higher, set maxt to that value
												if (dump.temper_sl[kk + 1] > maxt) maxt = dump.temper_sl[kk + 1]; // if lower layer at the beginning is higher, set maxt to that value
												if (temper_sl[kk + 1] > maxt) maxt = temper_sl[kk + 1]; // if lower layer at the end of the (sub) time step is higher, set maxt to that value
												// now set actual result to maxt
												temper_sl[kk] = maxt;
											}
											else if ((temper_sl[kk] < dump.temper_sl[kk] - 0.001) && // temperature decreases
												(temper_sl[kk] < T[i][j] * actNFactor - 0.001) && // but not more than the upper boundary condition
												(temper_sl[kk] < dump.temper_sl[kk + 1] - 0.001) && (temper_sl[kk] < temper_sl[kk + 1] - 0.001)) // and also not omre than the lower value of either the start orthe end value of the layer below 
											{
												double mint = dump.temper_sl[kk]; // minimum of all temperatures is set to old temperature of that layer
												if ((T[i][j] * actNFactor) < mint) mint = (T[i][j] * actNFactor); // if upper boundary condition is lower, set mint to that value
												if (dump.temper_sl[kk + 1] < mint) mint = dump.temper_sl[kk + 1];// if T of the lower layer at the beginning is lower, set mint to that value
												if (temper_sl[kk + 1] < mint) mint = temper_sl[kk + 1]; // if T of the lower layer at the end of the sub time step is lower, set mint to that value
												// no setting the actual value to the minimum of all checked temperatures
												temper_sl[kk] = mint;
											}
										}
										else
										{
											if (temper_sl[kk] < -70.0) {
												temper_sl[kk] = -70.0;
											}
											else if (temper_sl[kk] > +70.0) {
												temper_sl[kk] = 70.0;
											}
											else if ((temper_sl[kk] > dump.temper_sl[kk] + 0.001) && // temperature rises
												(temper_sl[kk] > dump.temper_sl[kk - 1] + 0.001) && (temper_sl[kk] > temper_sl[kk - 1] + 0.001) && // but should neither rise above the old nor above the new value of the layer above
												(temper_sl[kk] > dump.temper_sl[kk + 1] + 0.001) && (temper_sl[kk] > temper_sl[kk + 1] + 0.001)) // and should also neither rise above the old nor above the new value of the layer below
											{
												double maxt = dump.temper_sl[kk]; // maximum of all temperatures is set to old temperature of that layer
												if (dump.temper_sl[kk - 1] > maxt) maxt = dump.temper_sl[kk - 1]; // if T of upper layer at the beginning of the sub time step is higher, set maxt to that value
												if (temper_sl[kk - 1] > maxt) maxt = temper_sl[kk - 1]; // if T of upper layer at the end of the sub time step is higher, set maxt to that value
												if (dump.temper_sl[kk + 1] > maxt) maxt = dump.temper_sl[kk + 1]; // if lower layer at the beginning is higher, set maxt to that value
												if (temper_sl[kk + 1] > maxt) maxt = temper_sl[kk + 1]; // if lower layer at the end of the (sub) time step is higher, set maxt to that value
												// now set actual result to maxt
												temper_sl[kk] = maxt;

											}
											else if ((temper_sl[kk] < dump.temper_sl[kk] - 0.001) && // temperature decreses
												(temper_sl[kk] < dump.temper_sl[kk - 1] - 0.001) && (temper_sl[kk] < temper_sl[kk - 1] - 0.001) && // but should neither sink below the old nor below the new value of the layer above
												(temper_sl[kk] < dump.temper_sl[kk + 1] - 0.001) && (temper_sl[kk] < temper_sl[kk + 1] - 0.001)) // and should also neither sink below the old nor below the new value of the layer below
											{
												double mint = dump.temper_sl[kk]; // minimum of all temperatures is set to old temperature of that layer
												if (dump.temper_sl[kk - 1] < mint) mint = dump.temper_sl[kk - 1]; // if T of upper layer at the beginning of the sub time step is lower, set maxt to that value
												if (temper_sl[kk - 1] < mint) mint = temper_sl[kk - 1]; // if T of upper layer at the end of the sub time step is lower, set maxt to that value
												if (dump.temper_sl[kk + 1] < mint) mint = dump.temper_sl[kk + 1]; // if T of lower layer at the beginning of the sub time step is lower, set maxt to that value
												if (temper_sl[kk + 1] < mint) mint = temper_sl[kk + 1]; // if T of lower layer at the end of the sub time step is lower, set maxt to that value
												// now set actual result to mint
												temper_sl[kk] = mint;
											}
										}
									}
								}

								if (actDumpFiles < maxDumpFiles)
								{
									(*WSstream) << "DEBUG: \n ---- command line option -repairAfterDump active! --------\n"
										<< "setting the misbehaving parameter back to what is expected of it :-)" << "[end-entry]";
								}
								else if (actDumpFiles == maxDumpFiles)
								{
									(*WSstream) << "DEBUG: \n ---- command line option -repairAfterDump active! --------\n"
										<< "setting the misbehaving parameter back to what is expected of it :-)" 
										<< "\n continuing from now on without that message..." 
										<< "[end-entry]";
								}
							}
						}


						delete[] dump.ETR;
						delete[] dump.gwin;
						delete[] dump.gwex;
						delete[] dump.h_temp1;
						delete[] dump.theta_sl_old;
						delete[] dump.theta_average;
						delete[] dump.hh_sl;
						delete[] dump.hgeo_sl;
						delete[] dump.SE_sl;
						delete[] dump.qu_sl;
						delete[] dump.interflow;
						delete[] dump.thetasat;
						delete[] dump.dicken;
						delete[] dump.theta_temp;
						delete[] dump.thetares;
						delete[] dump.temper_sl;
						delete[] dump.theta_sl;
						delete[] dump.lambda_sl;
						delete[] dump.dEOverdT_sl;
						delete[] dump.heatcap_sl;
						delete[] dump.qHeat_sl;
						delete[] dump.temper_sl_new;
						delete[] dump.LUSE;
						delete[] dump.swes;
						delete[] dump.swel;
						delete[] dump.densSnow;
						delete[] dump.thickSnow;

					}
					
					/* update the infiltration amount in the riverflow_table */
					if ((gwflow_active)&&(surface_interaction>1))
					if (RIVLNK[i][j]>0)
					{
						// because RIVLNK[i][j] will not be unique for a given grid cell, the update of riverflow_table.infilt_amount must not be 
						// processed by more than one thread at a time
#pragma omp critical
						{
							if (((int)RIVLNK[i][j]<=0) || ((int)RIVLNK[i][j]>zonmaxwert))
							{
								(*WSstream) << "FATAL: \n\n -------- error! ---------------------------------------------\n\n"
								"There is a river link \""<<(int)RIVLNK[i][j]<<"\" from the river link grid\n"
								"which is not valid for the current routing structure -> there is no routing \n"
								"description for a river with code "<<(int)RIVLNK[i][j]<<"\n"
								<< "Please check your River link grid if there are cells defined as routing link "<<(int)RIVLNK[i][j]<<" and replace those values by nodata (-9999).\n"
								<<      "-------------------------------------------------------------\n\n" <<"[end-entry]";

								throw(199);
							}
							else
							{
								int lnknr=zarray[(int)RIVLNK[i][j]];
								if ((lnknr<=0) || (lnknr>zonmaxwert))
								{
									(*WSstream) << "FATAL: \n\n -------- error! ---------------------------------------------\n\n"
									"There is a river link \""<<(int)RIVLNK[i][j]<<"\" from the river link grid\n"
									"which is not valid for the current routing structure -> there is no routing \n"
									"description for a river with code "<<(int)RIVLNK[i][j]<<"\n"
									<< "Please check your River link grid if there are cells defined as routing link "<<(int)RIVLNK[i][j]<<" and replace those values by nodata (-9999).\n"
									<<      "-------------------------------------------------------------\n\n" <<"[end-entry]";
									throw(199);
								}
								riverflow_table[lnknr].infilt_amount+=GWINsub*cellsize*cellsize; /* subtracting infiltration in m3/s */
							}
						}
					}

#ifdef _WASIMDEBUG_
					/* iteration counts statistics */
					// this sequence will produce correct output only if this file is compiled without #pragma omp parallel
					// or if it is running with environment OMP_NUM_THREADS to 0 (or on a single processor machine without OMP_NUM_THREAD set)
					// otherwise multiple threads may lead to errors in the statistics by overwriting the values, the other thread eventually wrote
					// iteration_stat->wertzufuegen(0,0,TG,ziter);
#endif
					if ((ziter>100)&&(mpi_rank==0))
					{ 
						char logfilename[1024];
						sprintf(logfilename,"modell_%d.log",mpi_rank);
						if ((logdatei=fopen(logfilename,"a"))==NULL)
						{
							(*WSstream) << "FATAL: error opening log-file "<< logfilename <<"[end-entry]";
							throw(5);
						}
						fprintf(logdatei,"%2ld.%2ld., %2ld h r=%3d, c=%3d : ziter=%3d, max_err = %lg \n",tag,monat,stunde,i,j,ziter,m_err);
						fclose(logdatei);
					}


					/* old values of groundwater table and baseflow to estimate new recharge (difference to new values) */
					old_start_level = (int)GWLEVEL[i][j];
					old_gwstand     = (double)GWSTAND[i][j]+(double)hgeo_sl[0]; /* because of numerical precision only the groundwater level relative to the surface is stored! */
// not used anywhere	oldDELTAGW[i][j] *= timestep; /* converting the groundwater storage change from m/s in mm/timestep */
// not used anywhere	old_GWtheta     = GWtheta[i][j];

					/* calculating total interflow */
					for (int kk=1;kk<old_start_level;kk++) /* Interflow is calculated only on such layer bottoms,	which do'nt have contact to the groundwater */
					{
						QIgridsub+=interflow[kk];
					}

					if (do_LayeredSnow) // otherwise SnowRunoffsub was already considered above around line 7116 ff.
					{
						if (MAXPOND && MAXPOND[i][j]>0)
						{
							for (int kk=1+swesoffset;kk<=0;kk++)
							{
								POND[i][j] = POND[i][j] + interflow[kk]*timestep;
							}
							SnowRunoffsub=0;
						}
						else
						{
							for (int kk=1+swesoffset;kk<=0;kk++)
							{
								SnowRunoffsub+=interflow[kk];
							}
						}
					}

// ####################################################  starting from here the DELTAGW calculation was located before it was moved before the core ################
// just move it here if something is not working anymore to restore the old state.
// ####################################################  down to here the DELTAGW calculation was located before it was moved before the core ################

//					if (S0)
//					{
//						S0[i][j]=ST[s].getThetas((int)GWLEVEL[i][j]) - ST[s].get_theta((int)GWLEVEL[i][j], 3.45);
//						if (S0[i][j]<0.01) S0[i][j]=0.01;
//					}

					/* now calculating average Theta for the layer, the groundwater is within */
					//        s=(int)SOIL[i][j];
					gw_exact = hgeo_sl[0];
					k=start_level;


					if (GWLEVEL[i][j]>old_start_level) /* actual groundwater level drops below a layer boundary */
					{
// until 20100605		GWtheta[i][j]=ST[s].get_theta(k, // 3.45,
//						k==1 ? (dicken[k]+dicken[k+1])*0.5+h_temp : dicken[k]+h_temp);
						// now (20100605) we try the same value as for all other theta-estimations
						GWtheta[i][j]=ST[s].get_theta(k,3.45, thetares[k],thetasat[k]);
					}
					else
					{
						if (GWLEVEL[i][j]<old_start_level) /* rising groundwater level into a new layer: theta within the actual layer counts */
						GWtheta[i][j] =  gwflow_active?theta_sl_old[k]+(qu_sl[k]-qu_sl[k-1]-interflow[k]-qlat[k]-gwex[k]+gwin[k])*DT/dicken[k]:
													   theta_sl_old[k]+(qu_sl[k]-qu_sl[k-1]-interflow[k]-qlat[k]+gwin[k])*DT/dicken[k]; 
						else /* if groundwaterv remains in the same level, it has to be checked, if it raises or sinks */
						{
//							bilanz = gwflow_active ? qu_sl[k - 1] + a_etr[k] + a_evar[k] + DELTAGW[i][j] + interflow[k] + gwex[k] - gwin[k] : /* calculating balance */
//								qu_sl[k - 1] + a_etr[k] + a_evar[k] + DELTAGW[i][j] + interflow[k] - gwin[k];
							bilanz = gwflow_active ? qu_sl[k - 1] + a_etr[k] + a_evar[k] + interflow[k] + gwex[k] - gwin[k] : /* calculating balance */
												     qu_sl[k - 1] + a_etr[k] + a_evar[k] + interflow[k] - gwin[k];
							if (bilanz>0) /* if theta is decreasing */
							{
								/* GWTHETA has to decrease down to the theta, which is calculated by get_theta(...,dicken[k],...)
			(which is the above mentioned theta, allowing no gravity flow to the saturated layer) */
								GWtheta[i][j]-=bilanz*DT/dicken[k];
								//              theta_temp[k]=get_theta((int)(ST[s].anz_h),3.45,
// until 20100605				theta_temp[k]=ST[s].get_theta(k,/*0.3333*/dicken[k]);
								// now (20100605) we try the same value as for all other theta-estimations
								theta_temp[k]=ST[s].get_theta(k,3.45,thetares[k],thetasat[k]);

								if (GWtheta[i][j]<theta_temp[k])
								GWtheta[i][j]=theta_temp[k];
							}
						}
					}

					/* calculating the groundwater using the relative theta between saturation and drainage theta */
					if ((k==1)&&(theta_sl[1]>0.9995*thetasat[1])) 
						gw_exact=hgeo_sl[0];
					else
					{
						if (((double)thetasat[k] - (double)GWtheta[i][j]) > 0.0001) /* in order to avoid divisions by zero */
						{
//							gw_exact =
//								(double)hgeo_sl[k] +
//								(double)dicken[k] * ((double)theta_sl[k] - (double)GWtheta[i][j]) /
//								((double)thetasat[k] - (double)GWtheta[i][j]);
							double drain_theta = ST[s].get_theta(k, 3.45,thetares[k],thetasat[k]);
							if (GWtheta[i][j] > drain_theta)
							{
								gw_exact =
									hgeo_sl[k] +
									dicken[k] *
									(GWtheta[i][j] - drain_theta) /
									(thetasat[k] - drain_theta);
							}
							else
							{
								gw_exact = (double)hgeo_sl[k];
							}
						}
						else
						{
							gw_exact = (double)hgeo_sl[k - 1]; /* if GWtheta is almost like thetasat, then the groundwater is at the layers upper boundary */
						}
					}

					if (gw_exact<(double)hgeo_sl[k])
					{
						gw_exact=(double)hgeo_sl[k];
					}
					else
					{
						if (k>1)
						{
							if (gw_exact>(double)hgeo_sl[k-1])
							{
								gw_exact=(double)hgeo_sl[k-1];
							}
						}
					}
					GWSTAND[i][j]=(double)((double)gw_exact-(double)hgeo_sl[0]); /* little errors due to single precision double in GWSTAND may occure */

					/* now getting new values of groundwater table and recharge */
					///* to be removed
					GWNsub = groundwater_recharge(i,j,gw_exact,
					old_start_level,start_level,old_gwstand,
					hgeo_sl,theta_sl,thetasat,
					theta_sl_old,lateralinflow,dicken,qu_sl,GWINsub,GWEXsub,DT,timestep);
					// 			if ((i==10)&&(j==30))
					// 			(*WSstream) << "\nStunde " << stunde << " splitcnt " <<splitcnt << " GWEXsub nach GWRecharge " << GWEXsub*timestep << "[end-entry]";
					//*/
					if (GWNsub<0)
					{
						UPRISEsub=GWNsub;
					}
					else
					{
						PERCOLsub=GWNsub;
					}

					if (IRRIG) 
					if (IRRIGCODE[i][j]!=nodata)
					if ((irrigationtable[(int)IRRIGCODE[i][j]].origin==1)&&(IRRIGsub>=qlattemp*timestep))
					IRRIGsub-=qlattemp*timestep;

					// since infiltration after Green&Ampt was switched off, INFEX was always 0.
					// Infiltration Excess will now be set as p - qu_sl[0] for all qu_sl[0] > 0 (exfiltrating groudnwater will not be counted as infiltration excess)
					double infiltrationamount = (do_LayeredSnow && (swesoffset<0))?-qu_sl[0]:p; // if there is 
					if (qu_sl[0]<0)
					{
						if (MAXPOND && MAXPOND[i][j]>0)
						{
							INFEXsub  = (infiltrationamount + qu_sl[0]+SnowRunoffsub)*timestep-(POND[i][j]-MACROsub); // result in mm
						}
						else
						{
							INFEXsub  = (infiltrationamount + qu_sl[0]+SnowRunoffsub)*timestep; // result in mm
						}
						if (INFEXsub<0) INFEXsub=0;
						INFILTsub = -qu_sl[0]*timestep;
					}
					else
					{
						// rising water -> no real infintltration and no account for qu for infex
						if (MAXPOND && MAXPOND[i][j]>0)
						{
							INFEXsub=(infiltrationamount+SnowRunoffsub)*timestep-(POND[i][j]-MACROsub);
						}
						else
						{
							INFEXsub=(infiltrationamount+SnowRunoffsub)*timestep;
						}
						if (INFEXsub<0) INFEXsub=0;
						INFILTsub=0.0;
					}

					/* ############################################################### */
					/* do not change the order of the following 4 commands resp. blocks */
					/* ############################################################### */

					/* surface runoff is p + SnowRunoffsub - fabs(qu_sl[0]), p is the sum of p and infex ! */
					// Zwischenlösung, auch nicht korrekt:
					//			QDgridsub = lateralinflow[0]+p+INFEXsub/timestep+qu_sl[0]+SnowRunoffsub; // weil beim Berechnen der Infiltration INFEX von P abgezogen wurde,
					//			QDgridsub = lateralinflow[0]+INFEXsub/timestep+SnowRunoffsub; // weil beim Berechnen der Infiltration INFEX von P abgezogen wurde,
					QDgridsub = lateralinflow[0]+INFEXsub/timestep; // Snowrunoff was already added to INFEX

					/* using old thetas (in theta_sl_old) and new thetas in theta_sl! USING ALSO QDgrid AND conc_sl[z][0] !!! */
					///* to be removed
					new_concentrations(i,j,start_level,conc_sl,conc_sl_old,qu_sl,theta_sl,hgeo_sl,
								lateralinflow,gwin,gwex,theta_sl_old,thetasat,maxfromlr,
								dicken,interflow,qlat,ratio,iflfromlr,qlatfromlr,
								gwexfromlr,qufromlr,conc_lat,conc_gwex,conc_gwin,
								QDgridsub,QIgridsub,GWEXsub,DT,timestep);
					//*/
					// new_concentrations(i,j); /* using old thetas (in theta_sl_old) and new thetas in theta_sl! USING ALSO QDgrid AND conc_sl[z][0] !!! */

					for (int z=1;z<=tracer_count;z++) /* balancing tracer concentration/masses for direct runoff */
					{ /* CONCENTRATIONS for evaporating tracers remain constant (no mixing!),
			only MASSES of non evaporating tracers may have changed by infiltration*/
						if (!tracer_info[z].evaporating)
						{
							if (QDgridsub>0) /* infiltrating water: flux has negativ sign! */
							{
								if ((QDgridsub-qu_sl[0])!=0.0)
									trms    = conc_sl[z][0] *(-qu_sl[0]/(QDgridsub-qu_sl[0])); /* positive sign */
								else 
									trms=0.0;
								trmsmax = -tracer_info[z].max_conc*qu_sl[0]*timestep;        /* positive sign */
									if (trmsmax<trms) trms=trmsmax;
								conc_sl[z][0] -= trms;
							}
							//            else
							//              conc_sl[z][0]=0; /* no direct runoff: no masses -> set conc to 0 */
							/* fraction of infiltrating tracer mass: conc_sl[z][0] * -qu_sl[0]/QDgridsub; */
							/* transportable mass fraction for direct runoff */
							trmsmax = tracer_info[z].max_conc*QDgridsub*timestep;
							//            qd[z]=trmsmax<conc_sl[z][0]?trmsmax:conc_sl[z][0];
							//            conc_sl[z][0]-=qd[z];
							CONC_QD[z][i][j]=trmsmax<conc_sl[z][0]?trmsmax:conc_sl[z][0];
							conc_sl[z][0]-=CONC_QD[z][i][j];

						}
						else 
						{
							//            qd[z]=conc_sl[z][0];
							CONC_QD[z][i][j]=conc_sl[z][0];
						}
					}

					// when silting up is active, QDSUsub is the direct runoff fraction of this module
					// in that case: add the QDSU-value to the QDgrid to get the total surface runoff as result of the unsatzon model
					// since it was not already added to the input (m_bAddQdsuToP==false)
					// also, when layered Snow model is active and snow melt method is 5 (Energy balance), P from snow model is not used, so even if QDSUsub was 
					if ((this->m_bAddQdsuToP==false) || (this->do_LayeredSnow == true && this->do_TMod == true))
					{
						for (int z=1;z<=tracer_count;z++) /* balancing tracer concentration/masses for direct runoff */
						{
							if (!tracer_info[z].evaporating)
							{
								if (tracer_module[z].concflux!=NULL)
									CONC_QD[z][i][j]+=tracer_module[z].concflux[i][j]*QDSUsub*timestep;
							}
							else
							{
								if (tracer_module[z].concflux!=NULL)
									CONC_QD[z][i][j]=mixing(tracer_module[z].concflux[i][j],CONC_QD[z][i][j],QDSUsub,QDgridsub);
							}
						}
						QDgridsub+=QDSUsub/timestep; // convert QDSU into m/s!
					}
					
					if (MAXPOND && MAXPOND[i][j]>0) /* balancing the ponding storage only, if a max-pond-grid exists */
					{
						/* balancing the ponding water storage: if its contents is more than the max. value
			(here the ponding water storage is represented by QDgrid -> the POND-Storage was added to the
			rain before infiltration, so direct runoff is all of the exceeding infitration
			amount including the old content of the storage) */
						// note: when a maxpond grid exists, the QDSU-amount was added to the maxpond and thus to the rain 
						// So it's not neccessary to regard it here, again
						//			  if ((QDgridsub*timestep+(POND[i][j]-MACROsub))>MAXPOND[i][j])
						// commented 20090910			  if ((QDgridsub*timestep)>MAXPOND[i][j])
						if (POND[i][j]>MAXPOND[i][j])
						{
							QDgridsub=POND[i][j]/timestep;
							//				ft = (QDgridsub+(((POND[i][j]-MACROsub)-MAXPOND[i][j])/timestep))/(QDgridsub+POND[i][j]/timestep); /* fraction of tracer, which goes into runoff */
							ft = (QDgridsub-(MAXPOND[i][j]/timestep))/QDgridsub; /* fraction of tracer, which goes into runoff */
							//				QDgridsub=QDgridsub+(((POND[i][j]-MACROsub)-MAXPOND[i][j])/timestep); /* only if there is more water in the POND-storage than allowed */
							QDgridsub-=(MAXPOND[i][j]/timestep); /* only if there is more water in the POND-storage than allowed */
							POND[i][j]=MAXPOND[i][j];
							INFEXsub=0.0; /* if ponding occures, no infitration excess should be written to output, because the values are unrealistically large
								(before infitration the POND-storage was added to the rain in each time step, so very large infitration excess-values are possible!) */
							//            qd[0]=QDgridsub;
							for (int z = 1; z <= tracer_count; z++)
							{
								if (!tracer_info[z].evaporating)
								{
									trms = conc_sl[z][0] * ft; /* fraction of direct runoff on sum of direct runoff and pond storage */
									trmsmax = tracer_info[z].max_conc*QDgridsub*timestep; /* limit: what amount of salt can QD transport */
									if (trmsmax < trms) trms = trmsmax; /* is this limit less than actual fraction on salt referring to the fraction of flux? if yes, then limit mass flux to this value*/
									//                qd[z] = trms; /* mass of salt in direct runoff */
									CONC_QD[z][i][j] = trms; /* mass of salt in direct runoff */
									conc_sl[z][0] -= trms; /* this amount of mass has to be subtracted from mass storage */
								}
								else
								{
									//				qd[z]=conc_sl[z][0];
									CONC_QD[z][i][j] = conc_sl[z][0];
								}
							}
						}
						else
						{
							POND[i][j]=(lateralinflow[0]+p+qu_sl[0])*timestep; /* the amount of surface water (not infiltrated but not running off) remains in the POND-storage */
							// 				POND[i][j]=POND[i][j]-MACROsub+QDgridsub*timestep; /* the amount of surface water (not infiltrated but not running off) remains in the POND-storage */
							//				POND[i][j]=QDgridsub*timestep; /* the amount of surface water (not infiltrated but not running off) remains in the POND-storage */
							//            qd[0]=QDgridsub=0.0; /* and therefore there will be no surface runoff */
							QDgridsub=0.0; /* and therefore there will be no surface runoff */
							for (int z=1;z<=tracer_count;z++)
							{
								//              conc_sl[z][0]=qd[z];
								//				qd[z]=0.0;
								CONC_QD[z][i][j]=0.0;
							}
						}
					}
					
					// mean soil moisture as fraction of water filled pores over root depth, used for irrigation control
					//			SB_1_[i][j]=GetMeanTheta(s,dicken,theta_sl,0.0,maxrootdepth,0); // --> this is the relative soil moisture (0...1) for the requested range (here: root deoth)
					if (maxrootdepth>0)
						SB_1_[i][j]=GetMeanTheta(s,dicken,theta_sl,0.0,maxrootdepth,1)/maxrootdepth; // --> this is the absolute soil moisture as relative water content (0...Theta_sat) for the requested range (here: root deoth)
					else
						SB_1_[i][j]=theta_sl[1];
					// mean suction in root zone, used in LUtable for estimating 
					// stress conditions for dynamic phenology 3
					THROOT[i][j]=ST[(int)SOIL[i][j]].h_calc(1,GetMeanTheta(s,dicken,theta_sl,0.0,maxrootdepth,0));
					SB_2_[i][j]=1000*GetMeanTheta(s,dicken,theta_sl,0.0,9999.0,1); // soilmoisture over the entire soil profile - converted to absolute value as mm 
					ROOT[i][j]=maxrootdepth*1000; /* root depth in mm */

					/* storing the actual change in groundwater content of a
		cell for calculating the balance in the next time step */
//					oldDELTAGW[i][j]=DELTAGW[i][j];

					//			DELTAGW[i][j] *= timestep;

					if (HDRAIN)
					{
						QDRAINsub *= timestep;
					}

					QDgridsub  *= timestep;
					QIgridsub  *= timestep;
					if (gwflow_active)
					{
						GWINsub    *= timestep;
						GWEXsub    *= timestep;
					}

					// now adding up sub interval-data into grid cells
					if (QDRAIN) QDRAIN[i][j]+=QDRAINsub;
					if (GWEX) GWEX[i][j]+=GWEXsub;
					//		if ((i==10)&&(j==30))
					//		(*WSstream) << "\nStunde " << stunde << " lfd. Summe GWEX[i][j] " << GWEX[i][j] << "[end-entry]";
					if (GWIN) GWIN[i][j]+=GWINsub;
					QDgrid[i][j]+=QDgridsub;
					QIgrid[i][j]+=QIgridsub;
					if (MACRO) MACRO[i][j]+=MACROsub;
					if (methode>2)
					{
						if (SUBSTEPS) SUBSTEPS[i][j]=splitfactor;
					}
					for (int v=0;v<=iNumVegetationLayers;v++)
					{
						ETR[v][i][j]+=ETRsub[v];
					}
					GWN[i][j]+=GWNsub;
					PERCOL[i][j]+=PERCOLsub;
					UPRISE[i][j]+=UPRISEsub;
					EVAR[i][j]+=EVARsub;
					if (IRRIG) IRRIG[i][j]+=IRRIGsub;
					if (INFEX) INFEX[i][j]+=INFEXsub;
					if (INFILT) INFILT[i][j]+=INFILTsub;

					if (do_LayeredSnow)
					{
						SnowRunoff+=SnowRunoffsub*timestep;

						if (qu_sl[0]<=0)
						{
							PreciFromSnowMod+=(extraSnowRunoffsub);
						}
					}

					if (do_ArcVeg)
					{
						bool onIceWedge = false;
						for (k = levels; k >= 1; k--)
						{
							if ((ExcessIceFactor[i][j][k] > 0) && (interflow[k] > 0) || onIceWedge)
							{
								onIceWedge = true; // from here on upward, all hgeo etc. need to be changed
								// now change 
								// hgeo, 
								// dhm, 
								//  theta (should be the same value as at time step beginning)
								//  and ExcessIceFactor (runoff as interflow or other mains is taken from ExcessIce)

								double thickness_change = 0;
								if (ExcessIceFactor[i][j][k] > 0)
								{
									ExcessIceFactor[i][j][k] -= interflow[k] * timestep; // timestep is a factor for converting m/s in mm/timestep, so the result is in mm (WE)
									thickness_change = interflow[k] * timestep / 1000.0 * 1.1; // thickness change in mm (taken from ice, 1.1 is roughly the density of ice)
								}
								if ((ExcessIceFactor[i][j][k] <= 0) && (ExcessIceFactor[i][j][k] > -9999))
								{
									ExcessIceFactor[i][j][k] = -9999; // all excess ice molten, so this layer becomes a normal layer now
									thickness_change = 0;
								}

								if (thickness_change > 0)
								{
									// dicken[k] -= (interflow[k] * timestep) / 0.918 / 1000.0; // 0.918 is average density of ice, 1000.0 is a factor converting mm to m, so the result is change of thickness in m
									dicken[k] -= thickness_change;
									int s = (int)SOIL[i][j];
									if (dicken[k] < ST[s].getDzui()[k])
									{
										dicken[k] = ST[s].getDzui()[k];
									}
								}

								theta_sl[k] = theta[i][j][k]; // reset theta to the old value because we reduce the layer thickness as well
								hgeo_sl[k - 1] = hgeo_sl[k] + dicken[k];
								hgeo[i][j][k - 1] = hgeo_sl[k - 1];
								if ((k == 1) && (do_LayeredSnow))
								{ // change hgeoSnow when hgeo was changed, too
									double diff = hgeoSnow[i][j][levelsSnow] - hgeo[i][j][0];
									if (diff > 0)
									{
										for (int lv = levelsSnow ; lv >= 0; lv--)
										{
											hgeoSnow[i][j][lv] -= diff;
										}
									}
									// version that uses original thickness definition from layered snow properties, is slower
									// hgeoSnow[i][j][levelsSnow] = hgeo[i][j][0]; /* the zero-level contains then the elevation on the surface from the digital elevation model*/
									// for (int lv = levelsSnow - 1; lv >= 0; lv--)
									// {
									// 	hgeoSnow[i][j][lv] = hgeoSnow[i][j][lv + 1] + LayeredSnowProperties[1].getDzui()[lv];
									// }
								}
							}
						}
						if (onIceWedge)
						{
							gridinfo[2].pointer[i][j] = hgeo[i][j][0];
							if (gw_exact > (double)hgeo_sl[0])
							{
								gw_exact = (double)hgeo_sl[0];
								GWSTAND[i][j] = 0;
							}
						}
					}

					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// end of a sub interval
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}

				if (S0)
				{
					if (old_gwstand > gw_exact)
					{ // groundwater table tends to lower heads, so the storage capacity is that of drainable porosity
						S0[i][j] = ST[s].getThetas((int)GWLEVEL[i][j]) - ST[s].get_theta((int)GWLEVEL[i][j], 3.45,thetares[k], thetasat[k]);
						if (S0[i][j]<0.01) S0[i][j] = 0.01;
					}
					else
					{ // use the actual layers fillable theta when groundwater table if rising
						S0[i][j] = ST[s].getThetas((int)GWLEVEL[i][j]) - theta_sl[(int)GWLEVEL[i][j]];
						if (S0[i][j]<0.01) S0[i][j] = 0.01;
					}
				}

				if (this->do_TMod>0)
				{
					// calculating Thawdepth and the effective soil thickness which is within unfrozen conditions for later use in groundwater model
					for (int lv=1;lv<=ST[s].getLevels();lv++)
					{
						long ind=1;
						// todo: implement binary search here or simply create a reverse lookup table during start up
						while ((ST[s].SE_T[lv][ind]>(double)dTModSEForStat) && (ind<ST[s].getAnz_T())) ind++;
						T_SE08[lv]=-ind/(ST[s].getAnz_T()/10.0);
					}
					T_SE08[0] = T_SE08[1]; // index 0 was not calculated above but may be needed below, so take the values from the first layer for the 0-index
					int lvl=-1;
					double oldDepth = THAWTMOD[i][j];
					THAWTMOD[i][j]=0;
					do
					{
						lvl++;
						// if it is warmer in a lower level (except first level because of daily fluctuations), examine the other layers as well
						if ((lvl>0) &&
							(lvl<levels) &&
							(temper_sl[lvl+1]<=T_SE08[lvl]) &&
							(temper_sl[lvl]>T_SE08[lvl])
						   )
						{
							// interpolating between centers of the layers
							double d1=(hgeo_sl[lvl]+hgeo_sl[lvl-1])/2.0; // depth of center of upper layer
							double d2=(hgeo_sl[lvl+1]+hgeo_sl[lvl])/2.0; // depth of center of lower layer
							double a = (temper_sl[lvl+1]-temper_sl[lvl])/(d2-d1); // slope of the temperatures line between the two layers
							double c = temper_sl[lvl]-a*d1; // offset at depth = 0
							double dSE08 = (T_SE08[lvl]-c)/a; // depth at T_SE08
							THAWTMOD[i][j] = double(dSE08)-hgeo_sl[0]; // thaw depth (SE<=0.8 depth)
						} else if ((lvl==levels) &&
								   (temper_sl[lvl]>T_SE08[lvl]) &&
								   (TLOWBND[i][j]<=T_SE08[lvl])) { // only the lowest layer is (partly) frozen or the lower boundary condition is below T_SE08
							THAWTMOD[i][j] = hgeo_sl[lvl]-hgeo_sl[0];
						}
					}
					while (lvl<levels);
					if (gwflow_active)
					{
// old code before Mat changed the FROZENDEPTH calculation						
//						FROZENDEPTH[i][j] = 0; // FROZENDEPTH will contain the absolute thickness of the frozen part of the upper groundwater layer which is actually in the gorundwater, i.e. unsaturated zone doesn't count for this thickness
//						for (int lvl = 1; lvl <= levels; lvl++)
//						{
//							if ((hgeo_sl[0] - hgeo_sl[lvl - 1]) > -THAWTMOD[i][j]) // we are below the lowest numeric layer which is part of the active layer, which was already accounted for when initializing FROZENDEPTH[i][j]
//							{
//								if ((temper_sl[lvl] <= T_SE08[lvl]) && (lvl >= GWLEVEL[i][j]))
//								{
//									FROZENDEPTH[i][j] += dicken[lvl];
//								}
//							}
//						}
//						// now setting colmation between first and second aquifer depending on the frozen/unfrozen state of the lowest soil layer
//						if (temper_sl[levels] <= T_SE08[levels])
//						{
//							if (GW1COLMFROZEN) GW1COLMATION[i][j] = GW1COLMFROZEN[i][j];
//						}
//						else
//						{
//							if (GW1COLMUNFROZEN) GW1COLMATION[i][j] = GW1COLMUNFROZEN[i][j];
//						}
//
// MATs code starting here
						FROZENDEPTH[i][j] = 0; // FROZENDEPTH will contain the absolute thickness of the frozen part of the upper groundwater layer which is actually in the gorundwater, i.e. unsaturated zone doesn't count for this thickness
						for (int lvl = levels-1; lvl >0; lvl--)
						{
							if ((hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0 <= hgeo_sl[0] + GWSTAND[i][j]) // center of the current layer lower than gwtable
							{
								if ((temper_sl[lvl] <= T_SE08[lvl]) && (temper_sl[lvl + 1] > T_SE08[lvl + 1])) //if lower level is warmer
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // height of center of upper layer
									double d2 = (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0; // height of center of lower layer
									double a = (temper_sl[lvl + 1] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double a1 = (T_SE08[lvl + 1] - T_SE08[lvl]) / (d2 - d1); //slope of the threshold
									double c = temper_sl[lvl] - a * d1; // offset of the temp
									double c1 = T_SE08[lvl] - a1 * d1; // offset of the threshold
									double dSE08 = (c1 - c) / (a - a1); // height of intersect 
									FROZENDEPTH[i][j] += d1 - dSE08; // distance from the center of the upper layer
								}
								else if ((temper_sl[lvl] <= T_SE08[lvl]) && (temper_sl[lvl + 1] <= T_SE08[lvl + 1]))
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // height of center of upper layer
									double d2 = (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0; // height of center of lower layer
									FROZENDEPTH[i][j] += d1 - d2; // distance from the center of the upper layer
								}
								else if ((temper_sl[lvl] > T_SE08[lvl]) && (temper_sl[lvl + 1] <= T_SE08[lvl + 1])) //lower level is colder
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // height of center of upper layer
									double d2 = (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0; // height of center of lower layer
									double a = (temper_sl[lvl + 1] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double a1 = (T_SE08[lvl + 1] - T_SE08[lvl]) / (d2 - d1); //slope of the threshold
									double c = temper_sl[lvl] - a * d1; // offset of the temp
									double c1 = T_SE08[lvl] - a1 * d1; // offset of the threshold
									double dSE08 = (c1 - c) / (a - a1); // height of intersect 
									FROZENDEPTH[i][j] += dSE08 - d2; // distance from the center of the lower layer
								}
							}
							else if ((hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0 > hgeo_sl[0] + GWSTAND[i][j] && (hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0 <= hgeo_sl[0] + GWSTAND[i][j]) //gw table betwen two layers
							{
								if (temper_sl[lvl + 1] < T_SE08[lvl + 1]) //if lower level is cold
								{
									if (temper_sl[lvl] <= T_SE08[lvl])
									{
										FROZENDEPTH[i][j] += hgeo_sl[0] + GWSTAND[i][j] - (hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0;
									}
									else
									{
										double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // height of center of upper layer
										double d2 = (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0; // height of center of lower layer
										double a = (temper_sl[lvl + 1] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
										double a1 = (T_SE08[lvl + 1] - T_SE08[lvl]) / (d2 - d1); //slope of the threshold
										double c = temper_sl[lvl] - a * d1; // offset of the temp
										double c1 = T_SE08[lvl] - a1 * d1; // offset of the threshold
										double dSE08 = (c1 - c) / (a - a1); // height of intersect 

										FROZENDEPTH[i][j] += fmin(hgeo_sl[0] + GWSTAND[i][j], dSE08) - d2;
									}


								}
								else
								{
									if (temper_sl[lvl] < T_SE08[lvl])
									{
										double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // height of center of upper layer
										double d2 = (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0; // height of center of lower layer
										double a = (temper_sl[lvl + 1] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
										double a1 = (T_SE08[lvl + 1] - T_SE08[lvl]) / (d2 - d1); //slope of the threshold
										double c = temper_sl[lvl] - a * d1; // offset of the temp
										double c1 = T_SE08[lvl] - a1 * d1; // offset of the threshold
										double dSE08 = (c1 - c) / (a - a1); // height of intersect 
										if (dSE08 <= hgeo_sl[0] + GWSTAND[i][j])
										{
											FROZENDEPTH[i][j] += hgeo_sl[0] + GWSTAND[i][j] - dSE08;
										}
									}
									else if (temper_sl[lvl] == T_SE08[lvl] && temper_sl[lvl + 1] == T_SE08[lvl + 1]) // just in case both layers are at SE08
									{
										FROZENDEPTH[i][j] += hgeo_sl[0] + GWSTAND[i][j] - (hgeo_sl[lvl + 1] + hgeo_sl[lvl]) / 2.0;
									}
								}

							}
						}
						// checking for lowermost layer in the soil column
						int lvl = levels;
						if ((hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0 <= hgeo_sl[0] + GWSTAND[i][j]) // if gw table is above the lowest layer
						{
							if (temper_sl[lvl] <= T_SE08[lvl])
							{
								if (TLOWBND[i][j] <= T_SE08[lvl])
								{
									FROZENDEPTH[i][j] += (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0 - hgeo_sl[lvl]; //whole half of the layer
								}
								else
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // depth of center of  layer
									double d2 = hgeo_sl[lvl]; // depth of boundary condition
									double a = (TLOWBND[i][j] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double c = temper_sl[lvl] - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += d1 - dSE08; //only the frozen part
								}

							}
							else
							{
								if (TLOWBND[i][j] <= T_SE08[lvl])
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // depth of center of  layer
									double d2 = hgeo_sl[lvl]; // depth of boundary condition
									double a = (TLOWBND[i][j] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double c = temper_sl[lvl] - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += dSE08-d2; //only the frozen part
								}
							}

						}
						else // if gw table is underneath the center of the lowest layer
						{
							if (temper_sl[lvl] <= T_SE08[lvl])
							{
								if (TLOWBND[i][j] <= T_SE08[lvl])
								{
									FROZENDEPTH[i][j] += hgeo_sl[0]+GWSTAND[i][j] - hgeo_sl[lvl]; // onlythe part below the gw table
								}
								else
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // depth of center of  layer
									double d2 = hgeo_sl[lvl]; // depth of boundary condition
									double a = (TLOWBND[i][j] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double c = temper_sl[lvl] - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += (hgeo_sl[0] + GWSTAND[i][j])- fmin(hgeo_sl[0] + GWSTAND[i][j], dSE08); //only the frozen part
								}

							}
							else
							{
								if (TLOWBND[i][j] <= T_SE08[lvl])
								{
									double d1 = (hgeo_sl[lvl] + hgeo_sl[lvl - 1]) / 2.0; // depth of center of  layer
									double d2 = hgeo_sl[lvl]; // depth of boundary condition
									double a = (TLOWBND[i][j] - temper_sl[lvl]) / (d2 - d1); // slope of the temperatures line between the two layers
									double c = temper_sl[lvl] - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += fmin(hgeo_sl[0] + GWSTAND[i][j], dSE08) - d2; // only the frozen part within gw
								}
							}
						}
						// checing for the uppermost layer
						lvl = 0;
						if ((hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0 <= hgeo_sl[0] + GWSTAND[i][j]) // if gw table higher than the center of the highest layer
						{
							if (temper_sl[lvl] <= T_SE08[lvl])
							{
								if (Tsurf*actNFactor <= T_SE08[lvl])
								{
									FROZENDEPTH[i][j] += hgeo_sl[0] + GWSTAND[i][j] - (hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0; //whole half of the layer
								}
								else
								{
									double d1 = hgeo_sl[lvl];  // depth of boundary condition
									double d2 = (hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0; // depth of center of  layer
									double a = (Tsurf * actNFactor - temper_sl[lvl]) / (d1 - d2); // slope of the temperatures line between the two layers
									double c = Tsurf * actNFactor - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += fmin(hgeo_sl[0] + GWSTAND[i][j], dSE08) - d2; // only the frozen part within gw
								}
							}
							else
							{
								if (Tsurf*actNFactor <= T_SE08[lvl])
								{
									double d1 = hgeo_sl[lvl];  // depth of boundary condition
									double d2 = (hgeo_sl[lvl] + hgeo_sl[lvl + 1]) / 2.0; // depth of center of  layer
									double a = (Tsurf*actNFactor - temper_sl[lvl]) / (d1 - d2); // slope of the temperatures line between the two layers
									double c = Tsurf * actNFactor - a * d1; // offset at depth = 0
									double dSE08 = (T_SE08[lvl] - c) / a; // depth at T_SE08
									FROZENDEPTH[i][j] += (hgeo_sl[0] + GWSTAND[i][j]) - fmin(hgeo_sl[0] + GWSTAND[i][j], dSE08);
								}
							}
						}
						// now setting colmation between first and second aquifer depending on the frozen/unfrozen state of the lowest soil layer
						if (temper_sl[levels] <= T_SE08[levels])
						{
							if (GW1COLMFROZEN) GW1COLMATION[i][j] = GW1COLMFROZEN[i][j];
						}
						else
						{
							if (GW1COLMUNFROZEN) GW1COLMATION[i][j] = GW1COLMUNFROZEN[i][j];
						}
					}
				}

				// conversion must be done after all subintervals only
				DELTAGW[i][j] *= timestep*splitfactor;

				// the ETR2ETP quotient is now calculated. This is used in the irrigation model later to 
				// calculate required irrigation (if eg. ETR2ETP < 0.8, then the difference of SBHIGH - SB will be irrigated.
				// SBHIGH should be defined by crop specific suctions around field capacity (parameter HReduDry of each land use class)
				// Only dry stress is important here, so stress induced by too wet conditions must be replaced 
				// by values of 1 for the ETR2ETP grid. We do this by a simple check of theta in layer 2 against the SBHIGH value (assuming that 
				// in dry conditions there will be no extreme steep theta gradient in the soil, so the value of the 2nd layer is represantative 
				// enough to determine the difference between wet and dry stress
				if (irrigationtable)
				{
					if (SB_1_[i][j]<SBHIGH[i][j])
					{
					//			if (theta_sl[2]<SBHIGH[i][j])
// old code						
//						ETR2ETP[i][j] = (ETR[0][i][j]/(ETP[0][i][j]+EVAP[i][j]));

// new code
						double etploc=0,etrloc=0;	
						for (int lrs=1;lrs<=iNumVegetationLayers;lrs++)
						{
							etploc+=ETP[lrs][i][j];
							etrloc+=ETR[lrs][i][j];
						}
// AVAP and EVAR should not count for irrigation: only real vegetation induced extraction should be regarded compared to the potential amount
//						etploc+=EVAP[i][j];
//						etrloc+=EVAR[i][j];
						if (etploc > 0)
						{
							ETR2ETP[i][j] = etrloc / etploc;
						}
						else
						{
							ETR2ETP[i][j] = 1.0;
						}
					}
					else
						ETR2ETP[i][j] = 1.0;
				}
				else
				{
					if (ETR2ETP)
					{
						// irrigation table does not exist, but ETR2ETP is anyway initialized because a standard grid was found, so it must be calculated anyway
						double etploc = 0, etrloc = 0;
						for (int lrs = 1; lrs <= iNumVegetationLayers; lrs++)
						{
							etploc += ETP[lrs][i][j];
							etrloc += ETR[lrs][i][j];
						}
						if (etploc > 0)
						{
							ETR2ETP[i][j] = etrloc / etploc;
						}
						else
						{
							ETR2ETP[i][j] = 1.0;
						}
						if (ETR2ETP[i][j] > 1.0)
						{
							ETR2ETP[i][j] = 1.0; // limit ratio of real evapotranspiration to potential evapotranspiration to 1
						}
					}

				}
#ifdef _WASIMDEBUG_
				for (int z = 0; z <= tracer_count; z++)
				{
					maxbalance[z] = balance[z];
				}

				if (gwflow_active)
				{
					balance[0]-=GWEX[i][j];
					balance[0]+=GWIN[i][j];
					for (int z=1;z<=tracer_count;z++)
					{
						balance[z]-=tracer_info[z].evaporating?GWEX[i][j]*conc_gwex[z]:conc_gwex[z];
						balance[z]+=GWIN[i][j]*conc_gwin[z];
					}
				}
				balance[0]-=DELTAGW[i][j];
				balance[0]-=QIgrid[i][j];
				balance[0]-=QDgrid[i][j];
				balance[0]-=ETR[0][i][j]/substepcnt;
				//		balance[0]-=EVAR[i][j]; EVAR[i][j] was already added to ETR[0][i][j], so it must not be subtracted here a second time
				for (int z=1;z<=tracer_count;z++)
				{
					if (gwflow_active) 
					balance[z]-=DELTAGW[i][j]*gw_conc[z][i][j];
					else
					balance[z]-=tracer_info[z].evaporating?DELTAGW[i][j]*conc_lat[z]:conc_lat[z];
					balance[z]-=tracer_info[z].evaporating?QIgrid[i][j]*CONC_QI[z][i][j]:CONC_QI[z][i][j];
					balance[z]-=tracer_info[z].evaporating?QDgrid[i][j]*CONC_QD[z][i][j]:CONC_QD[z][i][j];
				}
				if (MAXPOND && MAXPOND[i][j]>0) 
				{
					deltapond-=POND[i][j];
					for (int z=1;z<=tracer_count;z++)
					{
						balance[z]-=tracer_info[z].evaporating?POND[i][j]*conc_sl[z][0]:conc_sl[z][0];
					}
				}
				for (int k=1+swesoffset;k<=levels;k++)
				{
					deltasb-=(theta_sl[k]*dicken[k]*1000.0);
					for (int z=1;z<=tracer_count;z++)
					balance[z]-=tracer_info[z].evaporating? theta_sl[k]*dicken[k]*1000*conc_sl[z][k] : conc_sl[z][k];
				}
				balance[0]+=deltapond;
				balance[0]+=deltasb;

				if (gwflow_active)
				{
					for (int z = 1; z <= tracer_count; z++)
					{
						balance[z] -= thetasat[levels] * (gw_thick[i][j] - (hgeo_sl[0] - hgeo_sl[levels]))*gw_conc[z][i][j] * 1000.0;
					}
				}
				//            balance[z]-=thetasat[levels]*(gw_thick[i][j]-(levels-1)*dicken[levels])*gw_conc[z][i][j]*1000.0;

//				if (fabs(balance[0])>0.01)
//				{
//#pragma omp critical
//					{
//						char logfilename[1024];
//						sprintf(logfilename,"modell_%d.log",mpi_rank);
//						logdatei=fopen(logfilename,"a+");
//						if (QDSU!=NULL)
//						{
//							fprintf(logdatei,"%ld %ld %ld %ld, Zelle %d %d: balance %lf, P %lf, QB %lf, QI %lf, QD %lf, QDSU %lf, ETR %lf, dPOND %lf, dSB %lf\n",
//							jahr,monat,tag,stunde,i,j,balance[0],P[i][j],DELTAGW[i][j],QIgrid[i][j],QDgrid[i][j],QDSU[i][j],etr*timestep,deltapond,deltasb);
//						}
//						else
//						{
//							fprintf(logdatei,"%ld %ld %ld %ld, Zelle %d %d: balance %lf, P %lf, QB %lf, QI %lf, QD %lf, QDSU %lf, ETR %lf, dPOND %lf, dSB %lf\n",
//							jahr,monat,tag,stunde,i,j,balance[0],P[i][j],DELTAGW[i][j],QIgrid[i][j],QDgrid[i][j],0.0,etr*timestep,deltapond,deltasb);
//						}
//						(*WSstream) << "b";
//						fclose(logdatei);
//					}
//				}
				for (int z = 1; z <= tracer_count; z++)
				{

					//          if ((fabs(balance[z])>1e-4)&&(!tracer_info[z].radioactive))
					if (fabs(balance[z]) / fabs(maxbalance[z]) > 1e-1)
					{
						char logfilename[1024];
						sprintf(logfilename, "modell_%d.log", mpi_rank);
						logdatei = fopen(logfilename, "a+");
						fprintf(logdatei, "%ld %ld %ld %ld, Zelle %d %d: tracer%d %lf (relation %lf)\n",
							jahr, monat, tag, stunde, i, j, z, balance[z], balance[z] / maxbalance[z]);
						//            (*WSstream) << z;
						fclose(logdatei);
					}
				}
#endif

				// write values back to the real stacks
				arraysize = sizeof(double)*(levels);
				memmove(&(theta[i][j][1]), &(m_theta_sl[levelsSnow+1]), arraysize);  // water content               
				memmove(&(hh[i][j][1])   , &(m_hh_sl[levelsSnow+1])   , arraysize);  // hydraulic heads             
				memmove(&(hgeo[i][j][1]) , &(m_hgeo_sl[levelsSnow+1]) , arraysize);  // geodaetic heads             
				memmove(&(qu[i][j][1])   , &(m_qu_sl[levelsSnow+1])   , arraysize);  // Fluxes between the layers   
				if (this->do_TMod > 0)
				{
					memmove(&(temper[i][j][1]), &(m_temper_sl[levelsSnow + 1]), arraysize);  // Fluxes between the layers   
					memmove(&(SEstack[i][j][1]), &(m_SE_sl[levelsSnow + 1]), arraysize);  // SE-Values
				}

				if (do_LayeredSnow)
				{

//					double *swes_sl = swes[i][j];
//					double *swel_sl = swel[i][j];
//					double *densSnow_sl = densSnow[i][j];
//					double extraSnowRunoff = 0;
//
//					upperIndexSnow[i][j]=levelsSnow+1;
//					swes_sl[0]=0;
//					swel_sl[0]=0;
//					
//					// update of the swes and swel stacks 
//					extraSnowRunoff = updateSnowStorageLayers(swes_sl,swel_sl,theta_sl,thetasat,dicken,SE_sl,densSnow_sl,temper_sl);
					
					arraysize = sizeof(double)*(levelsSnow);
					// values that come from the 3d-matrices and needed to be synchronized back to the original stacks after each cell is processed
					memmove(&(thetaSnow[i][j][1]), &(m_theta_sl[1]), arraysize);   
					memmove(&(hhydSnow[i][j][1]) , &(m_hh_sl   [1]), arraysize);   
					memmove(&(hgeoSnow[i][j][1]) , &(m_hgeo_sl [1]), arraysize);   
					memmove(&(quSnow[i][j][1])   , &(m_qu_sl   [1]), arraysize);   
					memmove(&(thickSnow[i][j][1]), &(m_dicken  [1]), arraysize);   
					if (this->do_TMod > 0)
					{
						memmove(&(temperSnow[i][j][1]), &(m_temper_sl[1]), arraysize);
						memmove(&(SEstack[i][j][1]), &(m_SE_sl[1]), arraysize);
					}

					for (int k=levelsSnow;k>=1;k--) 
					{
						swes[i][j][0]+=swes[i][j][k];
						swel[i][j][0]+=swel[i][j][k];
						if (swes[i][j][k]>0)
						{
							upperIndexSnow[i][j] = k;
						}
					}
//					QSNOW_OUT[i][j] = SSNO[i][j]-swes[i][j][0]; // change in SWES si the real snow melt without storage effects or so
					QSNOW_OUT[i][j] = SSTO[i][j]-swes[i][j][0]-swel[i][j][0]; // bor balance reasons we take the liquid fraction into account as well
					if (QSNOW_OUT[i][j]<0) QSNOW_OUT[i][j]=0; // when refreezing occuered from infiltrating rain water, swes will be greater than the old SSNOW, so we would get negative effective melt
					SSNO[i][j]=swes[i][j][0];
					SLIQ[i][j]=swel[i][j][0];
					SSTO[i][j]=SSNO[i][j]+SLIQ[i][j];

//					P[i][j]=PreciFromSnowMod+SnowRunoff; // P was taken from snow model as melt output. Since this value is overwritten for layered snow model by P_RAIN and swel[0] and melt is computed in unsatzonmodel only, the real melt will no be written back to the grid that is later statistically accumulated in snowmodel again
					P[i][j]=P_RAIN[i][j]+PreciFromSnowMod+SnowRunoff; // P was taken from snow model as melt output. Since this value is overwritten for layered snow model by P_RAIN and swel[0] and melt is computed in unsatzonmodel only, the real melt will no be written back to the grid that is later statistically accumulated in snowmodel again
					QDgrid[i][j]+=extraSnowRunoff; // the sudden snow melt will become direct runoff
					INFEX[i][j]+=extraSnowRunoff; // and must be counted here as well
				}
			} /* code was for valid zone grid cells only except for */

			    if ((zn[i][j] != nodata) && ((MPIREG[i][j] == mpi_rank) || (i == mpistartrow - 1)) && ((RIVLNK && (RIVLNK[i][j] == nodata)) || !RIVLNK || (i == mpistartrow - 1) || (surface_interaction < 2)) && (icevalue == 1.0) && !(doSpinUp && HRUMaster[i][j] == nodata)) // for completely glacierized cells only (at least partly glacierized)
			    //if ((zn[i][j] != nodata) && ((MPIREG[i][j] == mpi_rank) || (i == mpistartrow - 1)) && ((RIVLNK && (RIVLNK[i][j] == nodata)) || !RIVLNK || (i == mpistartrow - 1) ) && (icevalue == 1.0) && !(doSpinUp && HRUMaster[i][j] == nodata)) // for completely glacierized cells only (at least partly glacierized)
				{ /* this branch is reached, if the cell is valid in the zone grid and at the same time valid in the glacier grid (e.g., it is a glacier cell) */
				
				//	all these variables should not be set to 0 for partly glacierized cells, because for partly glacierized cells they must be kept until the statistics are built and grids are written

				TG=zarray[(int)zn[i][j]]; /* which subbasin */
				if (HDRAIN) QDRAIN[i][j]=0;
				if (gwflow_active)
				{
					GWIN[i][j]=0;
					GWEX[i][j]=0;
				}
				for (lay=1;lay<=iNumVegetationLayers;lay++)
				{
					ETR[lay][i][j]=0;
				}
				ETR[0][i][j]=0;
				EVAR[i][j]=0;
				INFEX[i][j]=0;
				INFILT[i][j]=0;
				QDgrid[i][j]=0;
				QIgrid[i][j]=0;
				if (MAXPOND) POND[i][j]=0;
				PERCOL[i][j]=UPRISE[i][j]=GWN[i][j]=MACRO[i][j]=0;
				if (methode>2) SUBSTEPS[i][j]=0;
				DELTAGW[i][j]=0;

				if (do_TMod)
				{
					// now setting soil temperatures which will be used when the cell may become partly unglacierized.
					// since we do not model the soil temperatures below glaciers, tha values will be taken from the nearest (at least partly) unglacierized cell
					// and corrected by the elevation difference (0.7°C per 100m)
					// 1. looking for the nearest neighbor with partly glacierized cells
					int radius = 0;
					double icevalue_iijj = 1;
					int ii = i;
					int jj = j;
					bool cellfound = false;
					int processedcells = 0;
					int min_ii = i;
					int max_ii = i;
					int min_jj = j;
					int max_jj = j;
					int nearest_i = i;
					int nearest_j = j;
					do
					{
						radius++; // starts with radius=1, will be incremented each round
						ii = i - radius;
						jj = j;
						if (ii < mpistartrow) ii = mpistartrow;
						// first run does the second half of the upper row
						while ((!cellfound) && (min_ii > ii) && (jj <= (j + radius)) && (jj <= ncols))
						{
							if ((MPIREG[ii][jj]==mpi_rank) && (zn[ii][jj] != nodata) && !(doSpinUp && HRUMaster[ii][jj] == nodata))
							{
								processedcells++;
								icevalue_iijj = ICE ? ICE[ii][jj] : 0; // valid range of icevalue now between 0...2
								if (icevalue_iijj < 1.0)
								{
									nearest_i = ii;
									nearest_j = jj;
									cellfound = true;
									break;
								}
							}
							jj++;
						}

						ii = i - radius + 1; // corner was already handled by the first while loop
						jj = j + radius;
						if (ii < mpistartrow) ii = mpistartrow;
						// second run does the column to the right, going down the rows
						while ((max_jj < jj) && (!cellfound) && (ii <= (i + radius)) && (ii <= nrows) && (ii <= mpiendrow))
						{
							if ((MPIREG[ii][jj] == mpi_rank) && (zn[ii][jj] != nodata) && !(doSpinUp && HRUMaster[ii][jj] == nodata))
							{
								processedcells++;
								icevalue_iijj = ICE ? ICE[ii][jj] : 0; // valid range of icevalue now between 0...2
								if (icevalue_iijj < 1.0)
								{
									nearest_i = ii;
									nearest_j = jj;
									cellfound = true;
									break;
								}
							}
							ii++;
						}

						ii = i + radius;
						jj = j + radius-1; // corner was already handled by the first while loop
						if (ii > mpiendrow) ii = mpiendrow;

						// third run does the row below the actual cell from right to left
						while ((max_ii < ii) && (!cellfound) && (jj >= (j - radius)) && (jj > 0))
						{
							if ((MPIREG[ii][jj] == mpi_rank) && (zn[ii][jj] != nodata) && !(doSpinUp && HRUMaster[ii][jj] == nodata))
							{
								processedcells++;
								icevalue_iijj = ICE ? ICE[ii][jj] : 0; // valid range of icevalue now between 0...2
								if (icevalue_iijj < 1.0)
								{
									nearest_i = ii;
									nearest_j = jj;
									cellfound = true;
									break;
								}
							}
							jj--;
						}

						ii = i + radius - 1; // corner was already handled by the first while loop
						jj = j - radius; 
						if (ii > mpiendrow) ii = mpiendrow;
						// fourth run does the column to the left, going upwards
						while ((min_jj > jj) && (!cellfound) && (ii >= (i - radius)) && (ii > 0) && (ii >= mpistartrow))
						{
							if ((MPIREG[ii][jj] == mpi_rank) && (zn[ii][jj] != nodata) && !(doSpinUp && HRUMaster[ii][jj] == nodata))
							{
								processedcells++;
								icevalue_iijj = ICE ? ICE[ii][jj] : 0; // valid range of icevalue now between 0...2
								if (icevalue_iijj < 1.0)
								{
									nearest_i = ii;
									nearest_j = jj;
									cellfound = true;
									break;
								}
							}
							ii--;
						}

						ii = i - radius; 
						jj = j - radius+1; // corner was already handled by the first while loop
						if (ii < mpistartrow) ii = mpistartrow;
						// last part is the left half of the row "radius" rows above the actual cell
						while ((min_ii > ii) && (!cellfound) && (jj < j))
						{
							if ((MPIREG[ii][jj] == mpi_rank) && (zn[ii][jj] != nodata) && !(doSpinUp && HRUMaster[ii][jj] == nodata))
							{
								processedcells++;
								icevalue_iijj = ICE ? ICE[ii][jj] : 0; // valid range of icevalue now between 0...2
								if (icevalue_iijj < 1.0)
								{
									nearest_i = ii;
									nearest_j = jj;
									cellfound = true;
									break;
								}
							}
							jj++; // this last value musn't be decrmented because we are now at j axactly and can use this as starting point for the next radius increment in the outer do loop
						}

						min_ii--;
						min_jj--;
						max_ii++;
						max_jj++;

					} while ((processedcells>0) && (!cellfound));
					// now taking the soil temperatures from that cell ii/jj
					if (cellfound)
					{
						double dT = (hgeo[nearest_i][nearest_j][1] - hgeo[i][j][1])*0.7 / 100.0; // 0.7°C per 100 m temperature decrement
						for (int k = 0; k <= levels; k++)
						{
							temper[i][j][k] = temper[nearest_i][nearest_j][k] + dT;
							// no SE value required since SE is always re-calculated and only stored if a regular run should write this stack
						}
					}
				}
			}

			if (doSpinUp==0)
			{
				if ((i==control_row) && (j==control_col))
				{
					if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
					{
						etrdetailbufcnt++;
						etrdetailbuffer[etrdetailbufcnt].j=jahr;
						//    etrdetailbuffer[etrdetailbufcnt].j=jahr%100;
						etrdetailbuffer[etrdetailbufcnt].m=monat;
						etrdetailbuffer[etrdetailbufcnt].t=tag;
						etrdetailbuffer[etrdetailbufcnt].s=stunde;
						etrdetailbuffer[etrdetailbufcnt].min = minute;

						interflowbufcnt++;
						interflowbuffer[interflowbufcnt].j=jahr;
						//    interflowbuffer[interflowbufcnt].j=jahr%100;
						interflowbuffer[interflowbufcnt].m=monat;
						interflowbuffer[interflowbufcnt].t=tag;
						interflowbuffer[interflowbufcnt].s = stunde;
						interflowbuffer[interflowbufcnt].min = minute;

						etrdetailbuffer[ etrdetailbufcnt].wert[0] = etrdetail[0]*timestep;
						for (int lc=1;lc<=levels;lc++)
						{
							etrdetailbuffer[ etrdetailbufcnt].wert[lc] = etrdetail[lc]*timestep;
							interflowbuffer[ interflowbufcnt].wert[lc] = interflow[lc]*timestep;
						}
						otherdatabufcnt++;
						otherdatabuffer[otherdatabufcnt].j=jahr;
						//      otherdatabuffer[otherdatabufcnt].j=jahr%100;
						otherdatabuffer[otherdatabufcnt].m=monat;
						otherdatabuffer[otherdatabufcnt].t=tag;
						otherdatabuffer[otherdatabufcnt].s = stunde;
						otherdatabuffer[otherdatabufcnt].min = minute;

						otherdatabuffer[ otherdatabufcnt].wert[1] = QDgrid?(QDgrid[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[2] = QIgrid?(QIgrid[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[3] = DELTAGW?(DELTAGW[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[4] = -qu_sl[0]*timestep;
						otherdatabuffer[ otherdatabufcnt].wert[5] = GWIN?(GWIN[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[6] = GWEX?(GWEX[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[7] = GWN?(GWN[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[8] = UPRISE?(UPRISE[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[9] = PERCOL?(PERCOL[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[10] = MACRO?(MACRO[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[11] = QDRAIN?(QDRAIN[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[12] = POND?(POND[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[13] = INFEX?(INFEX[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[14] = GWSTAND?(GWSTAND[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[15] = GWLEVEL?(GWLEVEL[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[16] = ROOT?(ROOT[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[17] = ALBEDO?(ALBEDO[i][j]):nodata;
						otherdatabuffer[ otherdatabufcnt].wert[18] = RSE?(RSE[i][j]):nodata;
						for (int veglay=1;veglay<=this->iNumVegetationLayers;veglay++)
						{
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+1] = LAI[veglay]?(LAI[veglay][i][j]):nodata;
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+2] = VBG[veglay]?(VBG[veglay][i][j]):nodata;
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+3] = RSC[veglay]?(RSC[veglay][i][j]):nodata;
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+4] = RSI[veglay]?(RSI[veglay][i][j]):nodata;
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+5] = WURZ[veglay]?(WURZ[veglay][i][j]):nodata;
							otherdatabuffer[ otherdatabufcnt].wert[18+(veglay-1)*6+6] = Z0[veglay]?(Z0[veglay][i][j]):nodata;
						}
					}
				}
			}

			qu1[i][j]=qu1_ext[i][j]=qu[i][j][0];

			//#endif    
		} // next grid column 

		// delete some temporary variables used for local threads (each row is processed by a single thread

		if (this->do_TMod>0)
		{
			delete[] m_SE_sl;			m_SE_sl			= NULL;
			delete[] m_T_SE08;			m_T_SE08		= NULL;
			delete[] m_lambda_sl;		m_lambda_sl		= NULL;
			delete[] m_qHeat_sl;		m_qHeat_sl		= NULL;
			delete[] m_heatcap_sl;		m_heatcap_sl	= NULL;
			delete[] m_deltaT_sl;		m_deltaT_sl		= NULL;
			delete[] m_temper_sl_new;	m_temper_sl_new = NULL;
			delete[] m_dEOverdT_sl;		m_dEOverdT_sl	= NULL;
			delete[] m_coeff_a;			m_coeff_a		= NULL;
			delete[] m_coeff_b;			m_coeff_b		= NULL;
			delete[] m_coeff_c;			m_coeff_c		= NULL;
			delete[] m_right_v;			m_right_v		= NULL;
			delete[] m_result_x;		m_result_x		= NULL;
		}
		delete[] ETPsub;
		delete[] ETRsub;

		delete[] conc_sl; conc_sl=NULL;
		for (int z=0;z<=tracer_count;z++) 
		{
			delete[] conc_sl_old[z]; conc_sl_old[z]=NULL;
		}
		delete[] conc_sl_old; conc_sl_old=NULL;

		delete[] conc_lat;			conc_lat		  = NULL;
		delete[] conc_gwin;			conc_gwin		  = NULL;
		delete[] conc_gwex;			conc_gwex		  = NULL;
		delete[] m_macrocapacity;		m_macrocapacity	    = NULL;
		delete[] m_macroinfiltration;	m_macroinfiltration = NULL;
		delete[] m_maxfromlr;			m_maxfromlr		    = NULL;
		delete[] m_qufromlr;			m_qufromlr		    = NULL;
		delete[] m_qlatfromlr;			m_qlatfromlr		= NULL;
		delete[] m_gwexfromlr;			m_gwexfromlr		= NULL;
		delete[] m_iflfromlr;			m_iflfromlr			= NULL;
		delete[] m_ratio;				m_ratio				= NULL;
		delete[] m_a_etr;				m_a_etr				= NULL;
		delete[] m_a_evar;				m_a_evar			= NULL;
		delete[] m_etrdetail;			m_etrdetail			= NULL;
		delete[] m_interflow;			m_interflow			= NULL;
		delete[] m_lateralinflow;		m_lateralinflow		= NULL;
		delete[] m_a_dicke;				m_a_dicke			= NULL;
		delete[] m_theta_temp;			m_theta_temp		= NULL;
		delete[] m_h_temp1;				m_h_temp1			= NULL;
		delete[] m_theta_average;		m_theta_average		= NULL;
		delete[] m_theta_sl_old;		m_theta_sl_old		= NULL;
		delete[] m_q_temp;				m_q_temp			= NULL;
		delete[] m_extractionweight;	m_extractionweight  = NULL;
		delete[] m_qlat;				m_qlat				= NULL;
		delete[] m_gwin;				m_gwin				= NULL;
		delete[] m_gwex;				m_gwex				= NULL;

		delete[] m_thetasat; m_thetasat	= NULL;
		delete[] m_thetares; m_thetares	= NULL;
		delete[] m_dicken;   m_dicken	= NULL;
		delete[] m_theta_sl; m_theta_sl	= NULL;
		delete[] m_hgeo_sl;  m_hgeo_sl	= NULL;
		delete[] m_hh_sl;    m_hh_sl	= NULL;
		delete[] m_qu_sl;    m_qu_sl    = NULL;
		if (this->do_TMod > 0)
		{
			delete[] m_temper_sl; m_temper_sl = NULL;
		}

		delete[] balance;				balance				= NULL;
		delete[] maxbalance;			maxbalance			= NULL;
	} /* next grid row */

	if (surface_interaction > 1)
	{
//#ifdef _DEBUG
//		if (numTempGridRows > 0)
//		{
//			ofstream rivflow;
//
//			stringstream ss;
//			ss << "unsatzon_rank_" << mpi_rank << "_riv_flow.dat";
//			std::string fn("");
//			ss >> fn;
//			rivflow.open(fn.c_str(), std::ios::app);
//			rivflow << "step " << jr << "-" << mt << "-" << tg << "-" << se << "-" << me << endl;
//			for (int i = 1; i <= zonmax; i++)
//			{
//				rivflow << "zone:\t" << i;
//				rivflow << "\tflowm3s:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].flowm3s << "\tinfilt_amount:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].infilt_amount << endl;
//			}
//			rivflow.close();
//		}
//#endif
		MPISyncRivLnkCells2(); // get the values from the extra grid row for in-order-processing of the river link cells and put them into the original cells.
//#ifdef _DEBUG
//		{
//			ofstream rivflow;
//
//			stringstream ss;
//			ss << "unsatzon_rank_" << mpi_rank << "_riv_flow_after_gather.dat";
//			std::string fn("");
//			ss >> fn;
//			rivflow.open(fn.c_str(), std::ios::app);
//			rivflow << "step " << jr << "-" << mt << "-" << tg << "-" << se << "-" << me << endl;
//			for (int i = 1; i <= zonmax; i++)
//			{
//				rivflow << "zone:\t" << i;
//				rivflow << "\tflowm3s:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].flowm3s << "\tinfilt_amount:\t" << std::setw(11) << std::setprecision(9) << riverflow_table[i].infilt_amount << endl;
//			}
//			rivflow.close();
//		}
//#endif
	}


	if (doSpinUp!=0)
	{
	// now setting the values of all non-master-HRU cells
#pragma omp parallel for schedule(runtime) private(i,j)
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
			{
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
				{  
					if (HRUMaster[i][j]==nodata)
					{ 
						// getting the coordinates of the master cell
						int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
						int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
						// filling the grid cells
						if (POND)		 POND[i][j]			= POND[imaster][jmaster];
						if (DRAINLEVEL)  DRAINLEVEL[i][j]	= DRAINLEVEL[imaster][jmaster];
						if (IRRIG)		 IRRIG[i][j]		= IRRIG[imaster][jmaster];
						if (THAWDEPTH)	 THAWDEPTH[i][j]	= THAWDEPTH[imaster][jmaster];
						if (MAXPUMP)	 MAXPUMP[i][j]		= MAXPUMP[imaster][jmaster];
						if (THAWTMOD)	 THAWTMOD[i][j]		= THAWTMOD[imaster][jmaster];
						if (FROZENDEPTH) FROZENDEPTH[i][j]	= FROZENDEPTH[imaster][jmaster];
						if (ETR2ETP)	 ETR2ETP[i][j]		= ETR2ETP[imaster][jmaster];
						if (P)			 P[i][j]			= P[imaster][jmaster];
						if (EVAR)		 EVAR[i][j]			= EVAR[imaster][jmaster];
						if (QDSU)		 QDSU[i][j]			= QDSU[imaster][jmaster];

						for (int lay=0;lay<=iNumVegetationLayers;lay++)
						{
							if (ETR[lay]) ETR[lay][i][j]	= ETR[lay][imaster][jmaster];
						}
						
						GWSTAND [i][j]	=  GWSTAND [imaster][jmaster];
						GWLEVEL [i][j]	=  GWLEVEL [imaster][jmaster];
						DELTAGW [i][j]	=  DELTAGW [imaster][jmaster];
						if (QDRAIN) QDRAIN  [i][j]	=  QDRAIN  [imaster][jmaster];
						if (gwflow_active)
						{
							GWEX    [i][j]	=  GWEX    [imaster][jmaster];
							GWIN    [i][j]	=  GWIN    [imaster][jmaster];
						}
						QDgrid  [i][j]	=  QDgrid  [imaster][jmaster];
						QIgrid  [i][j]	=  QIgrid  [imaster][jmaster];
						MACRO   [i][j]	=  MACRO   [imaster][jmaster];
						GWN     [i][j]	=  GWN     [imaster][jmaster];
						PERCOL  [i][j]	=  PERCOL  [imaster][jmaster];
						UPRISE  [i][j]	=  UPRISE  [imaster][jmaster];
						INFEX   [i][j]	=  INFEX   [imaster][jmaster];
						INFILT  [i][j]	=  INFILT  [imaster][jmaster];
						GWtheta [i][j]	=  GWtheta [imaster][jmaster];
						SB_1_   [i][j]	=  SB_1_   [imaster][jmaster];
						SB_2_   [i][j]	=  SB_2_   [imaster][jmaster];
						ROOT    [i][j]	=  ROOT    [imaster][jmaster];
						SUBSTEPS[i][j]	=  SUBSTEPS[imaster][jmaster];

						if (temper)
						{
							for (int z=0;z<=levels;z++)
							{
								temper[i][j][z] = temper[imaster][jmaster][z];
								SEstack[i][j][z] = SEstack[imaster][jmaster][z];
							}
						}
						for (int z=0;z<=levels;z++)
						{
							theta[i][j][z] = theta[imaster][jmaster][z];
							hgeo [i][j][z] = hgeo [imaster][jmaster][z];
							hh   [i][j][z] = hh   [imaster][jmaster][z];
							qu   [i][j][z] = qu   [imaster][jmaster][z];
						}
						
						if (do_ArcVeg)
						{
							for (int z=0;z<=levels;z++)
							{
								ExcessIceFactor[i][j][z] = ExcessIceFactor[imaster][jmaster][z];
							}
							if (ICEWEDGEUPPERDEPTH) ICEWEDGEUPPERDEPTH[i][j] =	ICEWEDGEUPPERDEPTH[imaster][jmaster];
							if (ICEWEDGELOWERDEPTH) ICEWEDGELOWERDEPTH[i][j] =	ICEWEDGELOWERDEPTH[imaster][jmaster];
							if (AVGDD)				AVGDD[i][j]				 =	AVGDD[imaster][jmaster];
							if (AVINMAX7DAVG)		AVINMAX7DAVG[i][j]		 =	AVINMAX7DAVG[imaster][jmaster];
							if (AVTDAVG)			AVTDAVG[i][j]			 =	AVTDAVG[imaster][jmaster]; // ** ArcVeg daily average Temperature grid
							if (AVINDAVG)			AVINDAVG[i][j]			 =	AVINDAVG[imaster][jmaster]; // ** ArcVeg daily average Inundation grid
							if (AVSNOWFREE)			AVSNOWFREE[i][j]		 =	AVSNOWFREE[imaster][jmaster]; // ** ArcVeg Snow Free Days grid
							if (AVPGP)				AVPGP[i][j]				 =  AVPGP[imaster][jmaster]; // ** ArcVeg number of plant groing periods
							if (AVSON)				AVSON[i][j]				 =  AVSON[imaster][jmaster]; // ** ArcVeg soil organic nitrogen
							if (AVNMINTOT)			AVNMINTOT[i][j]			 =  AVNMINTOT[imaster][jmaster]; // ** ArcVeg soil organic nitrogen total

							int kmax = MultiLayerLUTable[1].getNumFractions();
							if (AVINSCAL) // *** ArcVeg Inundation Scaling per land use fraction
							{
								for (int k = 1; k <= kmax; k++)
								{
									AVINSCAL[i][j][k]		 = AVINSCAL[imaster][jmaster][k]; // real stack i,j,k
									ARCVEGFRACTIONS[k][i][j] = ARCVEGFRACTIONS[k][imaster][jmaster]; // *** per plant fraction (land use fraction) stack of grids k,i,j
									ARCVEGBIOMASS[k][i][j]	 = ARCVEGBIOMASS[k][imaster][jmaster];   // *** per plant fraction (land use fraction) stack of grids k,i,j
									ARCVEGRETRANS[k][i][j]	 = ARCVEGRETRANS[k][imaster][jmaster];   // *** per plant fraction (land use fraction) stack of grids k,i,j
								}
							}

							if (AVINHLP) // *** a stack that holds 7 day values for inundation and a running average in level 0
							{
								for (int k = 0; k <= 7; k++)
								{
									AVINHLP[i][j][k] = AVINHLP[imaster][jmaster][k];
								}
							}

							if (AVNMIN) // *** ArcVeg nitrogen mineralization stack per plant growth period
							{
								for (int k = 1; k <= 5; k++)
								{
									AVNMIN[i][j][k] = AVNMIN[imaster][jmaster][k]; // real stack i,j,k
								}
							}


						}


						if (do_LayeredSnow)
						{
							for (int z = 0; z <= levelsSnow; z++)
							{
								swes[i][j][z] = swes[imaster][jmaster][z];
								swel[i][j][z] = swel[imaster][jmaster][z];
								thetaSnow[i][j][z] = thetaSnow[imaster][jmaster][z];
								hgeoSnow[i][j][z] = hgeoSnow[imaster][jmaster][z];
								hhydSnow[i][j][z] = hhydSnow[imaster][jmaster][z];
								quSnow[i][j][z] = quSnow[imaster][jmaster][z];
								thickSnow[i][j][z] = thickSnow[imaster][jmaster][z];
								densSnow[i][j][z] = densSnow[imaster][jmaster][z];
							}
							if (temperSnow)
							{
								for (int z = 0; z <= levelsSnow; z++)
								{
									temperSnow[i][j][z] = temperSnow[imaster][jmaster][z];
								}
							}
							upperIndexSnow[i][j] = upperIndexSnow[imaster][jmaster];
						}
					}
				}
			}
		}
	}

//	{
//		char name[255];
//		sprintf(name,"temperSnow_1.asc");
//		printstackfile(name, gridinfo[int(getindex(_Soil_Temperature_Stack_Snow_))].stackpointer,LayeredSnowProperties[1].getLevels(),true,false);
//	}
	// now make cumulative statistics
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			//    icevalue=ICE?ICE[i][j]:nodata;
			//	  if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue==nodata))
			icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
			if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))// performance is better when we can skip this section for pure ice cells 
			{
				TG = zarray[(int)zn[i][j]]; /* which subbasin */
				s = (int)SOIL[i][j];

//				QDsum+=QDgrid[i][j];
//				count++;

				/* updating the statistics */
				/* getting quick surface flow into the statistic register */
				if ((FZT) && (m_useSurfaceRouting==false))
				{
					fzt_class=(int)FZT[i][j];
				}
				else
				{
					fzt_class=1;
				}
				// index 0 is for water flux, other indizes are for tracerss, which are averaged weighted with water amount per cell 
				if (RelCellArea!=NULL)
				{
					zoneFZT_QD[0][TG][fzt_class].sx += RelCellArea[i][j]*QDgrid[i][j]*(1.0-icevalue)*substepcnt;
					zoneFZT_QD[0][TG][fzt_class].n  += (1.0-icevalue);
					zoneFZT_QD[0][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*(1.0-icevalue);

					zoneFZT_QI[0][TG][fzt_class].sx += RelCellArea[i][j]*QIgrid[i][j]*(1.0-icevalue)*substepcnt;
					zoneFZT_QI[0][TG][fzt_class].n  +=(1.0-icevalue);
					zoneFZT_QI[0][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*(1.0-icevalue);
				}
				else
				{
					zoneFZT_QD[0][TG][fzt_class].sx += QDgrid[i][j]*(1.0-icevalue)*substepcnt;
					zoneFZT_QD[0][TG][fzt_class].n  +=(1.0-icevalue);
					zoneFZT_QD[0][TG][fzt_class].RelCellSum  += (1.0-icevalue);

					zoneFZT_QI[0][TG][fzt_class].sx += QIgrid[i][j]*(1.0-icevalue)*substepcnt;
					zoneFZT_QI[0][TG][fzt_class].n  +=(1.0-icevalue);
					zoneFZT_QI[0][TG][fzt_class].RelCellSum  += (1.0-icevalue);
				}
				for (int z=1;z<=tracer_count;z++)
				{
					// for tracers the averaging has to be carried out weigted by absolute values 
					if (RelCellArea!=NULL)
					{
						zoneFZT_QD[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
							RelCellArea[i][j]*QDgrid[i][j]*CONC_QD[z][i][j]*(1.0-icevalue)*substepcnt : RelCellArea[i][j]*CONC_QD[z][i][j]*(1.0-icevalue)*substepcnt; // if evaporating: handling masses, now konverted to concentrations 
						zoneFZT_QD[z][TG][fzt_class].n  += QDgrid[i][j]*(1.0-icevalue);
						zoneFZT_QD[z][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*QDgrid[i][j]*(1.0-icevalue);

						zoneFZT_QI[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
							RelCellArea[i][j]*CONC_QI[z][i][j]*QIgrid[i][j]*(1.0-icevalue)*substepcnt : RelCellArea[i][j]*CONC_QI[z][i][j]*(1.0-icevalue)*substepcnt;
						zoneFZT_QI[z][TG][fzt_class].n  += QIgrid[i][j]*(1.0-icevalue);
						zoneFZT_QI[z][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*QIgrid[i][j]*(1.0-icevalue);
					}
					else
					{
						zoneFZT_QD[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
							QDgrid[i][j]*CONC_QD[z][i][j]*(1.0-icevalue)*substepcnt : CONC_QD[z][i][j]*(1.0-icevalue)*substepcnt; // if evaporating: handling masses, now konverted to concentrations 
						zoneFZT_QD[z][TG][fzt_class].n  += QDgrid[i][j]*(1.0-icevalue);
						zoneFZT_QD[z][TG][fzt_class].RelCellSum  += QDgrid[i][j]*(1.0-icevalue);
						
						zoneFZT_QI[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
							CONC_QI[z][i][j]*QIgrid[i][j]*(1.0-icevalue)*substepcnt : CONC_QI[z][i][j]*(1.0-icevalue)*substepcnt;
						zoneFZT_QI[z][TG][fzt_class].n  += QIgrid[i][j]*(1.0-icevalue);
						zoneFZT_QI[z][TG][fzt_class].RelCellSum  += QIgrid[i][j]*(1.0-icevalue);
					}
				}

				SB_1_stat->wertzufuegen(i,j,TG,SB_1_[i][j],(1.0-icevalue)); /* soilmoisture within root zone, relative value */
				SB_2_stat->wertzufuegen(i,j,TG,SB_2_[i][j],(1.0-icevalue)); /* soil moisture within the entire soil profile, abs. value */
//14.12.2014
				if (ICEOLD)
				{
					if ((ICEOLD[i][j]!=ICE[i][j])&&(ICEOLD[i][j]>=0)&&(ICE[i][j]>=0)) // to avoid nodata cells, we check for both values to be > 0
					{
						SB_2_IMMO_stat->wertzufuegen(i,j,TG,SB_2_[i][j]*(ICEOLD[i][j]-ICE[i][j]),(1.0-icevalue)); /* soil moisture within the entire soil profile, abs. value */
					}
					else
					{
						SB_2_IMMO_stat->wertzufuegen(i,j,TG,0,(1.0-icevalue)); /* soil moisture within the entire soil profile, abs. value */
					}
				}

				ROOTstat->wertzufuegen(i,j,TG,ROOT[i][j],(1.0-icevalue)); /* root depth */

				INFEXstat->wertzufuegen(i,j,TG,INFEX[i][j]*substepcnt,(1.0-icevalue)); // with Surface routing active, the INFEX value may come from ponding water, so it must be counted only once (and not multiplied with substepcnt)
				GWNstat->wertzufuegen(i,j,TG,GWN[i][j]*substepcnt,(1.0-icevalue));
				MACROstat->wertzufuegen(i,j,TG,MACRO[i][j]*substepcnt,(1.0-icevalue));

				if (methode>2) SUBSTEPSstat->wertzufuegen(i,j,TG,SUBSTEPS[i][j],(1.0-icevalue));
				timestep = 6.0e4*InV/substepcnt;    /* factor for changing the measure of fluxes from m/s to mm/timestep */
				INFILTstat->wertzufuegen(i,j,TG,INFILT[i][j]*substepcnt,(1.0-icevalue));
				GWSTANDstat->wertzufuegen(i,j,TG,GWSTAND[i][j],(1.0-icevalue)); /* Statistics for groundwater level beneath the floor */
				if (!gwflow_active)
				{
					if (IRRIG)
					{
						if (IRRIGCODE[i][j]!=nodata)
						{
							if (irrigationtable[(int)IRRIGCODE[i][j]].origin==1)
								QBstat->wertzufuegen(i,j,TG,(DELTAGW[i][j]-IRRIG[i][j]/timestep)*substepcnt,(1.0-icevalue)); /* no GWEX is used in thios version, so the mean DELTAGW is equal to the baseflow */
							else
								QBstat->wertzufuegen(i,j,TG,DELTAGW[i][j]*substepcnt,(1.0-icevalue)); /* no GWEX is used in thios version, so the mean DELTAGW is equal to the baseflow */
						}
						else
							QBstat->wertzufuegen(i,j,TG,DELTAGW[i][j]*substepcnt,(1.0-icevalue)); /* no GWEX is used in thios version, so the mean DELTAGW is equal to the baseflow */
					}
					else
						QBstat->wertzufuegen(i,j,TG,DELTAGW[i][j]*substepcnt,(1.0-icevalue)); /* no GWEX is used in thios version, so the mean DELTAGW is equal to the baseflow */
				}
				else
					QBstat->wertzufuegen(i,j,TG,GWEX[i][j]*substepcnt,(1.0-icevalue)); /* GWEX is baseflow */

				if (gwflow_active)
				{
					GWINstat->wertzufuegen(i,j,TG,GWIN[i][j]*substepcnt,(1.0-icevalue));
					GWEXstat->wertzufuegen(i,j,TG,GWEX[i][j]*substepcnt,(1.0-icevalue));
				}

				if (MAXPOND) PONDstat->wertzufuegen(i,j,TG,POND[i][j],(1.0-icevalue));
				if (HDRAIN)
				{
					QDRAINstat->wertzufuegen(i,j,TG,QDRAIN[i][j]*substepcnt,(1.0-icevalue));
				}

				// baseflow concentrations instead of frst soil layer 
				if (gwflow_active)
				{
					if (GWEX[i][j]>0)
					{
						for (int z=1;z<=tracer_count;z++)
						{
							if (RelCellArea!=NULL)
							{
								tracer_module[z].stat->zone[TG].sx += 
								tracer_info[z].evaporating ? 
								RelCellArea[i][j]*conc[z][i][j][(int)GWLEVEL[i][j]]*GWEX[i][j]*(1.0-icevalue)*substepcnt : 
								RelCellArea[i][j]*conc[z][i][j][(int)GWLEVEL[i][j]]*(1.0-icevalue)*substepcnt;

								tracer_module[z].stat->zone[TG].n  += 
								tracer_info[z].evaporating ? 
								GWEX[i][j]*(1.0-icevalue):
								1.0*(1.0-icevalue);

								tracer_module[z].stat->zone[TG].RelCellSum  += 
								tracer_info[z].evaporating ? 
								RelCellArea[i][j]*GWEX[i][j]*(1.0-icevalue): 
								RelCellArea[i][j]*(1.0-icevalue) ;
							}
							else
							{
								tracer_module[z].stat->zone[TG].sx += 
								tracer_info[z].evaporating ? 
								conc[z][i][j][(int)GWLEVEL[i][j]]*GWEX[i][j]*(1.0-icevalue)*substepcnt : 
								conc[z][i][j][(int)GWLEVEL[i][j]]*(1.0-icevalue)*substepcnt;

								tracer_module[z].stat->zone[TG].n  += 
								tracer_info[z].evaporating ? 
								GWEX[i][j]*(1.0-icevalue): 
								(1.0-icevalue) ;

								tracer_module[z].stat->zone[TG].RelCellSum  += 
								tracer_info[z].evaporating ? 
								GWEX[i][j]*(1.0-icevalue): 
								(1.0-icevalue) ;
							}
						}
					}
				}
				else
				{
					for (int z=1;z<=tracer_count;z++)
					{
						if (RelCellArea!=NULL)
						{
							tracer_module[z].stat->zone[TG].sx += 
							tracer_info[z].evaporating ? 
							RelCellArea[i][j]*conc[z][i][j][(int)GWLEVEL[i][j]]*DELTAGW[i][j]*(1.0-icevalue)*substepcnt : 
							RelCellArea[i][j]*conc[z][i][j][(int)GWLEVEL[i][j]]*(1.0-icevalue)*substepcnt;
							
							tracer_module[z].stat->zone[TG].n  += 
							tracer_info[z].evaporating ? 
							DELTAGW[i][j]*(1.0-icevalue):
							(1.0-icevalue) ;

							tracer_module[z].stat->zone[TG].RelCellSum  += 
							tracer_info[z].evaporating ? 
							RelCellArea[i][j]*DELTAGW[i][j]*(1.0-icevalue): 
							RelCellArea[i][j]*(1.0-icevalue) ;
						}
						else
						{
							tracer_module[z].stat->zone[TG].sx += 
							tracer_info[z].evaporating ? 
							conc[z][i][j][(int)GWLEVEL[i][j]]*DELTAGW[i][j]*(1.0-icevalue)*substepcnt : 
							conc[z][i][j][(int)GWLEVEL[i][j]]*(1.0-icevalue)*substepcnt;
							
							tracer_module[z].stat->zone[TG].n  += 
							tracer_info[z].evaporating ? 
							DELTAGW[i][j]*(1.0-icevalue): 
							(1.0-icevalue) ;

							tracer_module[z].stat->zone[TG].RelCellSum  += 
							tracer_info[z].evaporating ? 
							DELTAGW[i][j]*(1.0-icevalue): 
							(1.0-icevalue) ;
						}
					}
				}

				if (do_TMod>0)
				{
					// double thawdepth=0;
					// todo: calculate the thaw depth: this may be the depth where the temperature is 0°C 
					// or the depth where the SE_sl value is 0.1 or 0.5 or 0.9 or so
					// value can be taken from THAWTMOD, the dummy grid for the temperature stack, whic his now no dummy any longer...
					thawdepth_stat->wertzufuegen(i,j,TG,THAWTMOD[i][j]);
				}


			} // code was for valid data cells only with icevalue  1.0 (also partly glacerized cells, but weighed with weight of icevalue)
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue>0))// performance is better when we can skip this section for pure ice cells 
			{ // this branch is reached, if the cell is at least partly glacerized. For partly glacerized cells, the weight is 
				// complementary to the weight used above, so the statistics will remain to be based on the total number of cells in a subbasin
				TG=zarray[(int)zn[i][j]]; /* which subbasin */
				if ((FZT) && (m_useSurfaceRouting==false))
				{
					fzt_class=(int)FZT[i][j];
				}
				else
				{
					fzt_class=1;
				}

				if (RelCellArea!=NULL)
				{
					zoneFZT_QD[0][TG][(int)FZT[i][j]].sx += 0;
					zoneFZT_QD[0][TG][(int)FZT[i][j]].n  += icevalue;
					zoneFZT_QD[0][TG][(int)FZT[i][j]].RelCellSum  += RelCellArea[i][j]*icevalue;

					zoneFZT_QI[0][TG][(int)FZT[i][j]].sx += 0;
					zoneFZT_QI[0][TG][(int)FZT[i][j]].n  += icevalue;
					zoneFZT_QI[0][TG][(int)FZT[i][j]].RelCellSum  += RelCellArea[i][j]*icevalue;
				}
				else
				{
					zoneFZT_QD[0][TG][(int)FZT[i][j]].sx += 0;
					zoneFZT_QD[0][TG][(int)FZT[i][j]].n  += icevalue;
					zoneFZT_QD[0][TG][(int)FZT[i][j]].RelCellSum  += icevalue;

					zoneFZT_QI[0][TG][(int)FZT[i][j]].sx += 0;
					zoneFZT_QI[0][TG][(int)FZT[i][j]].n  += icevalue;
					zoneFZT_QI[0][TG][(int)FZT[i][j]].RelCellSum  += icevalue;
				}

//				SB_1_stat->wertzufuegen(i,j,TG,SB_1_[i][j],icevalue); /* soilmoisture within root zone, relative value */
//				SB_2_stat->wertzufuegen(i,j,TG,SB_2_[i][j],icevalue); /* soil moisture within unsaturated and saturated zone, absolute value */
				SB_1_stat->wertzufuegen(i,j,TG,0,icevalue); /* soilmoisture within root zone, relative value */
				SB_2_stat->wertzufuegen(i,j,TG,0,icevalue); /* soil moisture within unsaturated and saturated zone, absolute value */
//14.12.2014
				if (ICEOLD)
				{
					if ((ICEOLD[i][j]!=ICE[i][j])&&(ICEOLD[i][j]>=0)&&(ICE[i][j]>=0)) // to avoid nodata cells, we check for both values to be > 0
					{
						SB_2_IMMO_stat->wertzufuegen(i,j,TG,SB_2_[i][j]*(ICEOLD[i][j]-ICE[i][j]),icevalue); /* soil moisture within the entire soil profile, abs. value */
					}
					else
					{
						SB_2_IMMO_stat->wertzufuegen(i,j,TG,0,icevalue); /* soil moisture within the entire soil profile, abs. value */
					}
				}
				ROOTstat->wertzufuegen(i,j,TG,0,icevalue);

				INFEXstat->wertzufuegen(i,j,TG,0,icevalue);
				GWNstat->wertzufuegen(i,j,TG,0,icevalue);
				INFILTstat->wertzufuegen(i,j,TG,0,icevalue);
				MACROstat->wertzufuegen(i,j,TG,0,icevalue);

				if (methode>2) SUBSTEPSstat->wertzufuegen(i,j,TG,0,icevalue);
				GWSTANDstat->wertzufuegen(i,j,TG,0,icevalue); /* Statistics for groundwater level beneath the floor */
				QBstat->wertzufuegen(i,j,TG,0,icevalue); /* GWEX is baseflow */
					
				if (gwflow_active)
				{
					GWINstat->wertzufuegen(i,j,TG,0,icevalue);
					GWEXstat->wertzufuegen(i,j,TG,0,icevalue);
				}
					
				if (MAXPOND) PONDstat->wertzufuegen(i,j,TG,0,icevalue);
				if (HDRAIN)
				{
					QDRAINstat->wertzufuegen(i,j,TG,0,icevalue);
				}
			}
		}
	}

//(*WSstream) << "DEBUG: MPI["<<mpi_rank<<"] zoneFZT_QD[0][1][1].sx: "<<zoneFZT_QD[0][1][1].sx<< " zoneFZT_QD[0][1][1].RelCellSum: "<<zoneFZT_QD[0][1][1].RelCellSum<<"[end-entry]";

	if (control_row>0)
	{
		if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
		{
			// 	    (*WSstream) << "writing control row points in instance "<<mpi_rank<<" at row="<<control_row<<", col="<<control_col<<"[end-entry]";
			qubufcnt++;
			qubuffer[qubufcnt].j=jahr;
			//    qubuffer[qubufcnt].j=jahr%100;
			qubuffer[qubufcnt].m=monat;
			qubuffer[qubufcnt].t=tag;
			qubuffer[qubufcnt].s = stunde;
			qubuffer[qubufcnt].min = minute;

			thetabufcnt++;
			//    thetabuffer[thetabufcnt].j=jahr%100;
			thetabuffer[thetabufcnt].j=jahr;
			thetabuffer[thetabufcnt].m=monat;
			thetabuffer[thetabufcnt].t=tag;
			thetabuffer[thetabufcnt].s = stunde;
			thetabuffer[thetabufcnt].min = minute;

			hhydbufcnt++;
			hhydbuffer[hhydbufcnt].j=jahr;
			//    hhydbuffer[hhydbufcnt].j=jahr%100;
			hhydbuffer[hhydbufcnt].m=monat;
			hhydbuffer[hhydbufcnt].t=tag;
			hhydbuffer[hhydbufcnt].s = stunde;
			hhydbuffer[hhydbufcnt].min = minute;

			if (do_TMod>0 && this->iTModMethod > 0)
			{
				temperbufcnt++;
				temperbuffer[temperbufcnt].j=jahr;
				temperbuffer[temperbufcnt].m=monat;
				temperbuffer[temperbufcnt].t=tag;
				temperbuffer[temperbufcnt].s = stunde;
				temperbuffer[temperbufcnt].min = minute;
			}

			for (int i=1;i<=levels;i++)
			{
				thetabuffer[     thetabufcnt].wert[i] =     theta[control_row][control_col][i];
				qubuffer[        qubufcnt].wert[i] =        qu[control_row][control_col][i];
				hhydbuffer[      hhydbufcnt].wert[i] =        hh[control_row][control_col][i];
				if (do_TMod>0 && this->iTModMethod > 0)
				{
					temperbuffer[      temperbufcnt].wert[i] =        temper[control_row][control_col][i];
				}
			}
			if (   qubufcnt>=   qubufmax) flushbuffer(   qudatei,   qudateiname,   &qubufcnt,   qubuffer);
			if (thetabufcnt>=thetabufmax) flushbuffer(thetadatei,thetadateiname,&thetabufcnt,thetabuffer);
			if ( hhydbufcnt>= hhydbufmax) flushbuffer(hhyddatei,hhyddateiname,&hhydbufcnt,hhydbuffer);
			if (temperbufcnt>= temperbufmax) flushbuffer(temperfile,temperfilename,&temperbufcnt,temperbuffer);
			if ( otherdatabufcnt>= otherdatabufmax) 
				flushbuffer(otherdatadatei,otherdatadateiname,&otherdatabufcnt,otherdatabuffer,1,16+2+this->iNumVegetationLayers*6);
			if ( etrdetailbufcnt>= etrdetailbufmax) 
				flushbuffer(etrdetaildatei,etrdetaildateiname,&etrdetailbufcnt,etrdetailbuffer,0);
			if ( interflowbufcnt>= interflowbufmax) 
				flushbuffer(interflowdatei,interflowdateiname,&interflowbufcnt,interflowbuffer);

			for (int z=1;z<=tracer_count;z++)
			{
				concbufcnt[z]++;
				concbuffer[z][concbufcnt[z]].j=jahr;
				//      concbuffer[z][concbufcnt[z]].j=jahr%100;
				concbuffer[z][concbufcnt[z]].m=monat;
				concbuffer[z][concbufcnt[z]].t=tag;
				concbuffer[z][concbufcnt[z]].s = stunde;
				concbuffer[z][concbufcnt[z]].min = minute;
				if (tracer_info[z].evaporating)
				{
					for (int i=1;i<=levels;i++)
					concbuffer[z][ concbufcnt[z]].wert[i] =  conc[z][control_row][control_col][i];
				}
				else
				{
					s=(int)SOIL[control_row][control_col];
					for (int i=1;i<=levels;i++)
					concbuffer[z][ concbufcnt[z]].wert[i] = theta[control_row][control_col][i]>0.00001?
//						conc[z][control_row][control_col][i] / (1000.0*ST[s].getDzui()[i] * theta[control_row][control_col][i]) : /* converting tracer (salt) mass to concentration */
// convert getDzui() to dicken by initializing dicken again for the single cell
						conc[z][control_row][control_col][i] / (1000.0*ST[s].getDzui()[i] * theta[control_row][control_col][i]) : /* converting tracer (salt) mass to concentration */
						1e9;
				}

				if ( concbufcnt[z]>= concbufmax[z]) flushbuffer( concdatei[z], concdateiname[z], &(concbufcnt[z]), concbuffer[z]);
			}
		}
	}

}
//#pragma optimize ("p",off)

double unsatzonmodellklasse::updateSnowStorageLayers(double *swes_sl,double *swel_sl,double *theta_sl,double *thetasat,double *dicken,double *SE_sl,double *densSnow_sl,double *temper_sl)
{
	// return 0;
	double extraSnowRunoff=0;
	double maxFractionLiquid=this->dWaterHoldingCapacity;
	for (int k=1;k<=levelsSnow;k++) // going from uppermost layer to the bottom. If a layer is empty below a filled layer (e.g. because snow was warmed up above 0°C and thus there is no snow anymore (SE_sl=1), that filled layer and all above will be moved to the lower, empty layer in order to close gaps
	{
		double thetaUpperLayer=(k>1)?theta_sl[k-levelsSnow-1]:-9999; // comparing the actual layer with the upper layer makes sure that the next else-statement is done only for layers where the overlaying layer is empty, too
		if ((theta_sl[k-levelsSnow]>=0) || (thetaUpperLayer>=0))
		{
			double swesold=swes_sl[k];
			double swelold=swel_sl[k];
			double temp=(theta_sl[k-levelsSnow])*1000.0*dicken[k-levelsSnow];
			swel_sl[k] = temp*SE_sl[k-levelsSnow];
			swes_sl[k] = temp*(1.0-SE_sl[k-levelsSnow]);
			if (swes_sl[k]>0.05) // everything below 0.05 mm SWE should become melt anyway
			{
				// make sure, the liquid fraction is not higher than thetasat. If so, generate extraSnowRunoff
				if ((swes_sl[k] < swesold)&&(theta_sl[k-levelsSnow]>=thetasat[k-levelsSnow])) // some melt took place, so we should make the layer thinner in order to keep the density constant (compacting is not modeled here) 
				{
					dicken[k-levelsSnow] *= swes_sl[k]/swesold;
				}

				// removing excess water (more than snow should be able to hold) 
				if (swel_sl[k]>maxFractionLiquid*swes_sl[k])
				{
					double oldSWE = swel_sl[k]+swes_sl[k];
					extraSnowRunoff += (swel_sl[k]-maxFractionLiquid*swes_sl[k]);
					swel_sl[k]=maxFractionLiquid*swes_sl[k];
					// and adjusting the temperature to keep solid/liquid fractions after removal in equlibrium
					if((swel_sl[k]+swes_sl[k])>0)
					{
						SE_sl[k-levelsSnow]=swel_sl[k]/(swel_sl[k]+swes_sl[k]);
					}
					else
					{
						SE_sl[k-levelsSnow]=1;
					}

					double a = LayeredSnowProperties[1].getAlpha(k);
					double m = LayeredSnowProperties[1].getM(k);
					double n = LayeredSnowProperties[1].getN(k);
					double f = dTModFactor;

					// reverse function to calculate a temperature based on a given SE and the soils van-Genuchten-parameters
					double newT=temper_sl[k-levelsSnow];
					try
					{
						if ((SE_sl[k-levelsSnow]>0.0)&&(SE_sl[k-levelsSnow]<1.0))
						{
							double temp1 = pow(1.0-SE_sl[k-levelsSnow],1.0/m);
							double temp2 = pow(SE_sl[k-levelsSnow],1.0/(m*n));
							newT = -pow(temp1,1.0/n)/(a*f*temp2);
						}
						else
						{
							newT -= 2.0;
						}
					}
					catch (...)
					{
						newT -= 2.0;
					}
					temper_sl[k-levelsSnow]=newT;
					double newSWE = swel_sl[k]+swes_sl[k];
					dicken[k-levelsSnow] *= newSWE/oldSWE;
				}

				double swelMax = dicken[k-levelsSnow]*thetasat[k-levelsSnow]*1000.0-swes_sl[k];
				if (swelMax<0) swelMax=0;
				if (swelMax<swel_sl[k])
				{
					extraSnowRunoff += (swel_sl[k]-swelMax);
					swel_sl[k]=swelMax;
				}
				theta_sl[k-levelsSnow] = (swel_sl[k]+swes_sl[k])/(1000.0*dicken[k-levelsSnow]); // liquid and solid put fraction together here
				if (theta_sl[k-levelsSnow] > 1.0) theta_sl[k-levelsSnow]=1.0;
				densSnow_sl[k]  = ((swes_sl[k]+swel_sl[k])/dicken[k-levelsSnow])/1000.0; // result in kg/dm^3
				if (densSnow_sl[k] > 1.0) densSnow_sl[k]=1.0;
			}
			else
			{ // no solid snow anymore, so that layer needs to be eleminated and the water put to Snowrunoff
				extraSnowRunoff+=swel_sl[k];
				// reset the actual layer to "empty" defaults
				densSnow_sl[k]          = -9999;
				dicken[k-levelsSnow]    = -9999;
				swes_sl[k]              = 0;
				swel_sl[k]              = 0;
				temper_sl[k-levelsSnow] = -9999;
				theta_sl[k-levelsSnow]  = -9999; 
				// if there are filled layers above the actual one, move them down one layer each up to level 2, level 1 will get a special handling since it can be much thicker than the default (collects all snow above the layered part of the snow pack)
				for (int z=k;z>1;z--)
				{
					densSnow_sl[z]          = densSnow_sl[z-1];
					dicken[z-levelsSnow]    = dicken[z-levelsSnow-1];
					temper_sl[z-levelsSnow] = temper_sl[z-levelsSnow-1];
					swes_sl[z]              = swes_sl[z-1];
					swel_sl[z]              = swel_sl[z-1];
					theta_sl[z-levelsSnow]  = theta_sl[z-levelsSnow-1]; 
				}

				// handling of the uppermost layer
				if (k>1) // there must be at least two layers, otherwise all the snow remains in layer 1
				{   
					// shifting layer 1 downwards if if is thinner than the default thickness
					if (dicken[2-levelsSnow]<=LayeredSnowProperties[1].getDzui()[1]) // do not chenge getDzui to dicken, since this is the snow pack here
					{  // uppermost layer becomes empty if it was thinner than the default thickness, otherwise the uppermost layer is divided into two layers
						densSnow_sl[1]          = -9999;
						dicken[1-levelsSnow]    = -9999;
						swes_sl[1]              = 0;
						swel_sl[1]              = 0;
						temper_sl[1-levelsSnow] = -9999;
						theta_sl[1-levelsSnow]  = -9999; 

					} else { // part layer one into two layers if the thickness was thicker than the default layer thickness
						// no change to density 
						dicken[1-levelsSnow]  -= LayeredSnowProperties[1].getDzui()[1]; // do not chenge getDzui to dicken, since this is the snow pack here
						// no change to temperature 
						swes_sl[1]             = swes_sl[2] * (dicken[1-levelsSnow]/(dicken[1-levelsSnow]+dicken[2-levelsSnow]));
						swes_sl[2]            -= swes_sl[1];
						swel_sl[1]             = swel_sl[2] * (dicken[1-levelsSnow]/(dicken[1-levelsSnow]+dicken[2-levelsSnow]));
						swel_sl[2]            -= swel_sl[1];
						theta_sl[1-levelsSnow] = (swel_sl[1])/(1000.0*dicken[1-levelsSnow]); // liquid and solid put fraction together here
						theta_sl[2-levelsSnow] = (swel_sl[2])/(1000.0*dicken[1-levelsSnow]); // liquid and solid put fraction together here
					}
				}
			}
		}
		else
		{  // shortcut for all other empty layers above the first non-emptylayer: set all snow related values back to their defaults
			densSnow_sl[k]          = -9999;
			dicken[k-levelsSnow]    = -9999;
			swes_sl[k]              = 0;
			swel_sl[k]              = 0;
			temper_sl[k-levelsSnow] = -9999;
			theta_sl[k-levelsSnow]  = -9999; 
		}
	}

	// For (numerical) stability reasons, all snow layers are now examined for not becoming too thin or containing too less snow (i.e. becoming too less dense).
	// if that happens, such layers will be merged with the upper/lower layer to a new common layer
	for (int k=levelsSnow;k>1;k--) // moving upwards from the bottom of the snow pack to the top layer
	{
		double defaultdensity=LayeredSnowProperties[1].getDryDensity()[k]/1000.0; // in kg/dm^3
		double defaultthickness = LayeredSnowProperties[1].getDzui()[k];
		if (theta_sl[k-levelsSnow]>=0)
		{
			if ((densSnow_sl[k]<(defaultdensity*0.8)) &&(theta_sl[k-levelsSnow-1]>0)) // occurs due to melt e.g. from soil heat flux
			{
				densSnow_sl[k]=defaultdensity;
				dicken[k-levelsSnow] = ((swes_sl[k]+swel_sl[k])/densSnow_sl[k])/1000.0; // result in m
				if (dicken[k-levelsSnow]<0.25*defaultthickness) // if the layer becomes too thin, merge it with the layer above
				{
					if (k>1)
					{ // mixing/adding the actual layer with the layer above
						// mixing temperatures --> todo 20170429: attention: we mix snow and water here together, although they have different hear capacities. This may lead to problems later -> check
						double Tnew = ((swes_sl[k]*dTModCIce0+swel_sl[k]*dTModCWater) * temper_sl[k-levelsSnow] + (swes_sl[k-1]*dTModCIce0+swel_sl[k-1]*dTModCWater) * temper_sl[k-levelsSnow-1]) 
									/ ((swes_sl[k]+swes_sl[k-1])*dTModCIce0 + (swel_sl[k]+swel_sl[k-1])*dTModCWater);
						temper_sl[k-levelsSnow]  = Tnew;
						dicken[k-levelsSnow]    += dicken[k-levelsSnow-1];
						swes_sl[k]              += swes_sl[k-1];
						swel_sl[k]              += swel_sl[k-1];
						// new density as a mixture of ols and new swes+swel and the accumulated thickness
						densSnow_sl[k]           = ((swes_sl[k]+swel_sl[k])/dicken[k-levelsSnow])/1000.0; // result in kg/dm^3
						theta_sl[k-levelsSnow]   = (swes_sl[k]+swel_sl[k])/(1000.0*dicken[k-levelsSnow]); // liquid and solid put fraction together here
					}
					// this is similar to the algorith for moving down layers, except that the masses are now added and the temperatures mixed
					for (int z=k-1;z>1;z--)
					{
						densSnow_sl[z]         = densSnow_sl[z-1];
						dicken[z-levelsSnow]   = dicken[z-levelsSnow-1];
						temper_sl[z-levelsSnow]= temper_sl[z-levelsSnow-1];
						swes_sl[z]             = swes_sl[z-1];
						swel_sl[z]             = swel_sl[z-1];
						theta_sl[z-levelsSnow] = theta_sl[z-levelsSnow-1]; 
					}
								
					if (k>1) // there must be at least two layers, otherwise all the snow remains in layer 1
					{   
						// shifting layers downwards if one layer was emtied completely
						if (dicken[2-levelsSnow]<=LayeredSnowProperties[1].getDzui()[1])
						{  // uppermost layer becomes empty if it was thinner than the default thickness, otherwise the uppermost layer is divided into two layers
							densSnow_sl[1]          = -9999;
							dicken[1-levelsSnow]    = -9999;
							temper_sl[1-levelsSnow] = -9999;
							swes_sl[1]              = 0;
							swel_sl[1]              = 0;
							theta_sl[1-levelsSnow]  = -9999;
						} else {
							// no change to density 
							dicken[1-levelsSnow]  -= LayeredSnowProperties[1].getDzui()[1];
							// no change to temperature 
							swes_sl[1]          = swes_sl[2] * (dicken[1-levelsSnow]/(dicken[1-levelsSnow]+dicken[2-levelsSnow]));
							swes_sl[2]         -= swes_sl[1];
							swel_sl[1]          = swel_sl[2] * (dicken[1-levelsSnow]/(dicken[1-levelsSnow]+dicken[2-levelsSnow]));
							swel_sl[2]         -= swes_sl[1];
							theta_sl[1-levelsSnow] = (swes_sl[1]+swel_sl[1])/(1000.0*dicken[1-levelsSnow]); // liquid and solid put fraction together here
							theta_sl[2-levelsSnow] = (swes_sl[2]+swel_sl[2])/(1000.0*dicken[2-levelsSnow]); // liquid and solid put fraction together here
						}
					}
				}
			}
		}
		else
		{
			break; // first layer without snow: break the for loop
		}
	}
	return extraSnowRunoff;
}

// will be called from post_run after all sub steps are processed - meaning, that also surface routing ran, so the actual state of inundation is always after the surface routing was done, even if we are in the soil model here
void unsatzonmodellklasse::run_arcveg(long me, long se, long tg, long mt, long jr)
{
	// defined in teilmod.h (class unsatzonmodellklasse)
	// int ArcVegStartJD, ArcVegEndJD, ArcVegWaitPeriod;
	// bool bArcVegWaiting, bArcVegSnowFree;
	// 1. calculate actual Julian Day
	bool runArcVegNow = false;
	int actJD = 0;
	for (int i = 1; i < mt; i++)
	{
		actJD += mon[i];
	}
	actJD += tg;

	// compute actual minute counter for the start of the interval
	// compute actual minute counter for the end of the interval

	if (!bArcVegActive)
	{
		if (ArcVegStartJD < ArcVegEndJD) // normal case in northern hemisphere
		{
			if ((actJD >= ArcVegStartJD) && (actJD < ArcVegEndJD))
			{
				bArcVegActive = true;
			}
		}
		else
		{ // southern hemisphere: user should set ArcVegStartJD to a day that is greater than the end day
			if ((actJD >= ArcVegStartJD) || (actJD < ArcVegEndJD))
			{
				bArcVegActive = true;
			}
		}

		if (bArcVegActive && (actJD == ArcVegStartJD))
		{
			// initialize all grids, because we start no into the active ArcVeg period
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				// initialize all cells from first to last in this row with 0. 
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
					{
						AVINHLP1[i][j] = 0;
						AVGDD[i][j] = 0;
						AVINMAX7DAVG[i][j] = 0;
						AVTDAVG[i][j] = 0;
						AVINDAVG[i][j] = 0;
						AVSNOWFREE[i][j] = 0;
						AVPGP[i][j] = 0;
						AVNMINTOT[i][j] = 0;
						for (int lv = 0; lv <= 8; lv++)
						{
							// initialilzing everything in the excessIceStack with nodata first
							AVINHLP[i][j][lv] = 0;
						}
						for (int lv = 0; lv <= 5; lv++)
						{
							// initialilzing everything in the mineralization stack with 0 first
							AVNMIN[i][j][lv] = 0;
						}
						for (int lv = 0; lv <= MultiLayerLUTable[1].getNumFractions(); lv++)
						{
							AVINSCAL[i][j][lv] = 0;
						}
					}
				}
			}
		}
		// otherwise do nothing
	}

	if (bArcVegActive)
	{
		if (ArcVegStartJD < ArcVegEndJD) // normal case in northern hemisphere
		{
			if ((actJD >= ArcVegEndJD) || (actJD < ArcVegStartJD))
			{
				bArcVegActive = false; // end of active period reached
				runArcVegNow = true;
			}
		}
		else
		{ // southern hemisphere: user should set ArcVegStartJD to a day that is greater than the end day
			if ((actJD >= ArcVegEndJD) && (actJD < ArcVegStartJD))
			{
				bArcVegActive = false; // end of active period reached
				runArcVegNow = true;
			}
		}

		// if active:
		if (bArcVegActive)
		{
			int actMinuteMin = 60 * (se - 1 + int(hour_offset)) + me;
			int actMinuteMax =(int)( actMinuteMin + InV);
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					// check if cell is snow free and AVSNOWFREE is still 0. If so, set AVSNOWFREE to the actual JD
					if ((SSTO[i][j] < 0.01) && (AVSNOWFREE[i][j] == 0))
					{
						AVSNOWFREE[i][j] = actJD;
					}

					if (AVSNOWFREE[i][j] > 0)
					{
						// start calculating the 
						if (actJD >= AVSNOWFREE[i][j] + ArcVegWaitPeriod - 7) // start 7 days before the wait period is over in order to have valid inundation averages over 7 days already
						{
							// calculate 7-day-moving average inundation -> per day, 
							// i.e. only if the new day value (tg) is different from the old day value of tg, 
							// using the value of AVINDAVG and setting levels of AVINHLP (1...7, then 0 as running average); reset, when a new day started
							if (actMinuteMax == 1440) // end of a day is reached, so daily statistics can be computed
							{
								AVINHLP[i][j][0] -= AVINHLP[i][j][7] / 7.0; // the 7th day before the actual day will become the 8th day and will thus be removed from the moving average stack
								AVINHLP[i][j][7] = AVINHLP[i][j][6];
								AVINHLP[i][j][6] = AVINHLP[i][j][5];
								AVINHLP[i][j][5] = AVINHLP[i][j][4];
								AVINHLP[i][j][4] = AVINHLP[i][j][3];
								AVINHLP[i][j][3] = AVINHLP[i][j][2];
								AVINHLP[i][j][2] = AVINHLP[i][j][1];
								AVINHLP[i][j][1] = AVINDAVG[i][j];
								AVINHLP[i][j][0] += AVINHLP[i][j][1] / 7.0; // adding the newest day to the moving average layer
							}
							// updates for Average Temperature and Inundation will be doine per time step
							// calculate daily average temperature -> per interval, 
							// update AVTDAVG (AVTDAVG = ArcVegTemperatureDailyAVeraGe), 
							// use actual value AVTDAVG and minutecounters 
							AVTDAVG[i][j] = (AVTDAVG[i][j] * actMinuteMin + T[i][j] * InV) / double(actMinuteMax);

							// calculate daily average inundation -> per interval, 
							// update AVINDAVG (= ArcVegINundationDailyAVeraGe), 
							// use actual value AVTDAVG and minutecounters 
							AVINDAVG[i][j] = (AVINDAVG[i][j] * actMinuteMin + QDSU[i][j] * InV) / double(actMinuteMax);

						}

						// if AVSNOWFREE+ArcVegWaitPeriod >= actJD and we are at the end of a day, calculate actual values for maximum inundation 7-day-average and GDD
						if ((actJD >= (AVSNOWFREE[i][j] + ArcVegWaitPeriod)) && (actMinuteMax == 1440))
						{
							// calculate max 7d inundation -> also only once per day, i.e. if actMinuteMax = 1440, 
							// using the value of AVINHLP[i][j][0], i.e. the previously computed moving average over the last 7 days
							if (AVINMAX7DAVG[i][j] < AVINHLP[i][j][0])
							{
								AVINMAX7DAVG[i][j] = AVINHLP[i][j][0];
							}
							// calulate AVGDD from the last AVTDAVG (add the value, if positive; then, reset AVTDAVG for the beginning new day)
							double GDDThreshold = 0.0;
							if (AVTDAVG[i][j] > GDDThreshold) // positive degree days only. 
							{
								AVGDD[i][j] += (AVTDAVG[i][j] - GDDThreshold);
							}
						}
					}
				}
			}
		}
		else
		{
			if (runArcVegNow) // should run only once directly after setting bArcVegActive = false
			{
				// core of arcveg could be implemented here or called here
				// parameters are among others the AVGDD grid and the AVINMAX7DAVG grid
				(*WSstream) << "INFO: ArcVeg Core runs now" << "[end-entry]";

				int numfractions = MultiLayerLUTable[1].getNumFractions();
				srand(1); // to make sure we will always have the same start value when comparing different runs
				// implementing ArcVeg
				for (int i = mpistartrow; i <= mpiendrow; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						{
							int s = (int)SOIL[i][j];
							bool isArcVegType = ST[s].getIsArcVegType();
							if (isArcVegType)
							{
								// saving the old BioMass of all fractions. Later, this is be used for change of soil thickness for e.g. mosses
								double* OldBioMass = new double[numfractions + 1];
								for (int k = 1; k <= numfractions; k++)
								{
									OldBioMass[k] = ARCVEGBIOMASS[k][i][j];
								}

								int multiyear = 0;
#ifdef _DEBUG
								
								int maxmultiyear = 100; 
								AVGDD[i][j] += 600.0; 
#else
								// to be save, we set the number of years in production mode to 1 and also to not shift the degree days sum to any higher value to get more growing periods.
								int maxmultiyear = 1; 
								AVGDD[i][j] += 0.0; 
#endif // DEBUG

								for (multiyear = 1; multiyear <= maxmultiyear; multiyear++)
								{

									// saving the old BioMass of all fractions. Later, this is be used for change of soil thickness for e.g. mosses
									double* OldBioMassActYear = new double[numfractions + 1];
									for (int k = 1; k <= numfractions; k++)
									{
										OldBioMassActYear[k] = ARCVEGBIOMASS[k][i][j];
									}
									AVNMINTOT[i][j] = 0;
									for (int lv = 0; lv <= 5; lv++)
									{
										// initialilzing everything in the mineralization stack with 0 first
										AVNMIN[i][j][lv] = 0;
									}
									for (int lv = 0; lv <= MultiLayerLUTable[1].getNumFractions(); lv++)
									{
										AVINSCAL[i][j][lv] = 0;
									}

									// hier beginnt der Code für ein normales Jahr.
									// 1. calculating number of growth periods between 2 and 5
									
									if (AVGDD[i][j] < 250.0) {
										AVPGP[i][j] = 2;
									}
									else if ((AVGDD[i][j] >= 250.0) && (AVGDD[i][j] < 450.0)) {
										AVPGP[i][j] = 3;
									}
									else if ((AVGDD[i][j] >= 450.0) && (AVGDD[i][j] < 750.0)) {
										double normgdd = (AVGDD[i][j] - 450.0) / (750.0 - 450.0); // scaling gdd to a value between 0 and 1
										// 450 scales to 0, 750 scales to 1
										// the actual AVGDD determines the probabilities of being one pgpmax vs.the other
										double temp = (double(rand()) / double(RAND_MAX)); // random value between 0 and 1
										if (temp < normgdd) // if random value is less than scaled gdd
										{
											AVPGP[i][j] = 3;
										}
										else // if random value is greater than scaled gdd
										{
											AVPGP[i][j] = 4;
										}
									}
									else if ((AVGDD[i][j] >= 750.0) && (AVGDD[i][j] < 1000.0)) {
										// values in this range were found in both subzones D and E
										double normgdd = (AVGDD[i][j] - 750.0) / (1000.0 - 750.0); // scaling gdd to a value between 0 and 1
										// 750 scales to 0, 1000 scales to 1
										// the actual gdd determines the probabilities of being one pgpmax vs.the other
										double temp = (double(rand()) / double(RAND_MAX)); // random value between 0 and 1
										if (temp < normgdd) // if random value is less than scaled gdd
										{
											AVPGP[i][j] = 4;
										}
										else // if random value is greater than scaled gdd
										{
											AVPGP[i][j] = 5;
										}

									}
									else if (AVGDD[i][j] > 1000.0) {
										AVPGP[i][j] = 5;
									}

									// 2. Inundation marker depending on max average 7-day inundation, is used in ArcVeg as scaling factor between 0 and 1
									// full inundation : wasim_ind > 30 regardless of wedgePresent
									if (AVINMAX7DAVG[i][j] > 300.0)
									{
										for (int k = 1; k <= numfractions; k++)
										{
											int use = (int)LUSE[k][i][j];
											AVINSCAL[i][j][k] = LUtable[use].GetArcVegInFactFull();
										}
									}
									else if ((AVINMAX7DAVG[i][j] >= 0.0) && (AVINMAX7DAVG[i][j] <= 300.0))
									{ // inundation in mm 

										// some inundation on ice wedge
										bool onwedge = false;
										for (int lv = 1; lv <= levels; lv++)
										{
											if (ExcessIceFactor[i][j][lv] > 0)
											{
												onwedge = true;
												break;
											}
										}
										if (onwedge)
										{
											for (int k = 1; k <= numfractions; k++)
											{
												int use = (int)LUSE[k][i][j];
												// (1 - wasim_ind / 30)* UDIceWdg + (wasim_ind / 30) * SomeInund (Scale the inundation growth scalar between UDIceWeg and SomeInund)
												double actUDIceWdg = LUtable[use].GetArcVegInUDIceWdg();
												double actFactSome = LUtable[use].GetArcVegInFactSome();
												AVINSCAL[i][j][k] = ((1.0 - (AVINMAX7DAVG[i][j] / 300.0)) * actUDIceWdg) + ((AVINMAX7DAVG[i][j] / 300.0) * actFactSome);
												//inundation(i, j, :) = ((1 - (wasim_ind(i, j) . / 30)).*Plants(:, 20)) + ((wasim_ind(i, j) . / 30).*Plants(:, 21));
											}
										}
										else
											// some inundation not on ice wedge
											// if (wedgePresent(i, j) == 0)
										{
											for (int k = 1; k <= numfractions; k++)
											{
												int use = (int)LUSE[k][i][j];
												double actFactSome = LUtable[use].GetArcVegInFactSome();
												AVINSCAL[i][j][k] = ((1.0 - (AVINMAX7DAVG[i][j]) / 300.0) * 1.0) + ((AVINMAX7DAVG[i][j] / 300.0) * actFactSome);
												// (1 - wasim_ind / 30) * 1 + (wasim_ind / 30) * SomeInund // (Scale the inundation growth scalar between 1 and SomeInund)
												// inundation(i, j, :) = (1 - (wasim_ind(i, j) . / 30)) + ((wasim_ind(i, j) . / 30).*Plants(:, 21));
											}
											// if No surface water, growth scalar will either be 1 (not on an ice wedge) or UDIceWdg(on an ice wedge)
										}
									}

									// 3. calculate mineralization
									for (int k = 1; k <= (int)AVPGP[i][j]; k++)
									{
										if (k <= AVNumPGPs)
										{
											double actNMin = AVNMinParam[AVActClimZon][k];
											AVNMIN[i][j][k] = AVSON[i][j] * actNMin; // ToDo: Parameter for Mineralization per climate zone and growth period from control file
										}
										else
										{
											(*WSstream) << "FATAL: ArcVeg actual PGP " << k << "count is larger than allowed in the control file (AVNumPGPs = " << AVNumPGPs << ")" << "[end-entry]";
											throw(1505);
										}
										if (AVNMIN[i][j][k] < 0)
										{
											AVNMIN[i][j][k] = 0;
										}
										AVNMINTOT[i][j] += AVNMIN[i][j][k];
									}
									AVSON[i][j] -= AVNMINTOT[i][j];

									// 4. Seeding
	//								for k = 1:planttypes
	//									if rand < Plants(k, clim + 7) % seedling establishment f(climate)
	//										biom(i, j, k) = biom(i, j, k) + 0.1;% add seed biomass to biomass for each grid cell
	//									end
	//								end
									int numFractions = MultiLayerLUTable[1].getNumFractions(); // get the number of fractions for that code from the multilayer land use table

									for (int k = 1; k <= numFractions; k++) // for each of the fractions...
									{
										int luse = (int)LUSE[k][i][j]; // what is the code of this fraction in the landuse table
										double temp = (double(rand()) / double(RAND_MAX)); // random value between 0 and 1
										if (temp < LUtable[luse].GetArcVegSeeedMaxProb()[AVActClimZon]) // if the random number is smaller than the maximum probablity read for that actual landuse table entry (i.e. the actual fraction) for the given climate
										{
											ARCVEGBIOMASS[k][i][j] += 0.1; // add seed biomass to actual biomass for each grid cell, so if any biomass fraction were zero, we have the chance (depending on climate and the random number) to start that fraction again
										}
									}

									// 5. growth
									for (int p = 1; p <= (int)AVPGP[i][j]; p++)
									{
										double nut = 0.0; // NUT is the total available N for each grid cell and PGP
										if (p > 1)
										{
											AVNMIN[i][j][p] += AVNMIN[i][j][p - 1];
										}
										nut = AVNMIN[i][j][p];

										double* need = new double[6];
										double needsum = 0.0;

										for (int k = 1; k <= numFractions; k++)
										{
											int luse = (int)LUSE[k][i][j];
											double uptakeefficiency = LUtable[luse].GetArcVegUptakeEfficiency()[p];
											need[k] = ARCVEGBIOMASS[k][i][j] * AVINSCAL[i][j][k] * uptakeefficiency;
											needsum += need[k];
										}

										double* Fraction2BiomassFactor = new double[numFractions + 1];
										double sumFraction2BiomassFactor = 0;
										for (int k = 1; k <= numFractions; k++) // Go through each planttype
										{
											Fraction2BiomassFactor[k] = ARCVEGFRACTIONS[k][i][j] / ARCVEGBIOMASS[k][i][j];
											sumFraction2BiomassFactor += Fraction2BiomassFactor[k];

											double uptake = 0;
											if (needsum > 0) // A check in case there is no nutrient demand (i.e. grid cell is barren)
											{
												if (needsum > nut) // IF the total NEED (demand) > available N (NUT)
												{
													double needprop = 0;
													needprop = need[k] / needsum;
													uptake = needprop * nut; // THEN plants obtain N in proportion to their demand
												}
												else
												{
													int luse = (int)LUSE[k][i][j];
													double uptakeefficiency = LUtable[luse].GetArcVegUptakeEfficiency()[p];
													uptake = ARCVEGBIOMASS[k][i][j] * uptakeefficiency * AVINSCAL[i][j][k]; // Otherwise, N uptake is equal to the NEED(demand)->like line 29
												}
											}
											else
											{
												uptake = 0;
											}
											// Plant types retranslocate nitrogen(from prior year)
											// Retrans is calculated in the Senesce module(from prior year)
											double retrans = 0;

											// while arcveg sets the retrans amount to 0 for the first year, wasim needs initialized grids (either read from file or initialized by using the default values and writecode 8 in standard grids section)
											retrans = ARCVEGRETRANS[k][i][j] / AVPGP[i][j];	  			// Retranslocated N is split evenly across PGPs

											int luse = (int)LUSE[k][i][j];

											// wrong formula						double uptakeefficiency = LUtable[luse].GetArcVegUptakeEfficiency()[p];
											//										ARCVEGBIOMASS[k][i][j] += uptakeefficiency * (uptake + retrans);  // Plant type biomass grows using N from Nmin and retrans,
											double bnratio = LUtable[luse].GetArcVegBNratio();
											ARCVEGBIOMASS[k][i][j] += bnratio * (uptake + retrans);  // Plant type biomass grows using N from Nmin and retrans,
											// debug without retrans										ARCVEGBIOMASS[k][i][j] += bnratio * (uptake);  // Plant type biomass grows using N from Nmin and retrans,
																					// we make a check here if actual biomass exceeds the maximum biomass for this plant type and limit it to this value
											if (ARCVEGBIOMASS[k][i][j] > LUtable[luse].GetArcVegMaxBio())
											{
												ARCVEGBIOMASS[k][i][j] = LUtable[luse].GetArcVegMaxBio();
											}

											// Multiplied by the plant type biomass : N ratio
											AVNMIN[i][j][p] -= uptake;
										}

										// change fractions
//										double sumArcVegFractions = 0;
//										for (int k = 1; k <= numFractions; k++) // Go through each planttype
//										{
//											ARCVEGFRACTIONS[k][i][j] = ARCVEGBIOMASS[k][i][j] * Fraction2BiomassFactor[k];
//											sumArcVegFractions += ARCVEGFRACTIONS[k][i][j];
//										}
//										// normalize fractions
//										double temp = 0;
//										for (int k = 1; k <= numFractions; k++) // Go through each planttype
//										{
//											if (sumArcVegFractions != 0)
//											{
//												ARCVEGFRACTIONS[k][i][j] /= sumArcVegFractions;
//											}
//											else
//											{
//												ARCVEGFRACTIONS[k][i][j] = 0;
//											}
//											temp += ARCVEGFRACTIONS[k][i][j];
//										}

										delete[] need;
										delete[] Fraction2BiomassFactor;
									}

									// 6. Senesce module
									double* retransNs = new double[numFractions + 1];
									double sumRetransN = 0;
									for (int k = 1; k <= numfractions; k++)
									{
										int use = (int)LUSE[k][i][j];
										double f1 = LUtable[use].GetArcVegRetran();
										double f2 = LUtable[use].GetArcVegFoliarProp();
										double f3 = LUtable[use].GetArcVegPropN();
										double scalingfactor = (LUtable[use].GetArcVegMaxBio() - ARCVEGBIOMASS[k][i][j]) / LUtable[use].GetArcVegMaxBio();
										retransNs[k] = ARCVEGBIOMASS[k][i][j] * f1 * f2 * f3 * scalingfactor;
										// for each plant type in each plot: biomass x foliar proportion x N prop.x prop.retranslocated
										// scale retranslocation by the maximum biomass(col T); as the plant type gets closer to its maximum biomass it retranslocates less
										ARCVEGRETRANS[k][i][j] = retransNs[k];
										sumRetransN += retransNs[k]; //sum for all planttypes retrans 
									}
									delete[] retransNs;
									// subtract the retranslocated nitrogen total for each grid cell from the total son pool
									// in other words, this nitrogen does NOT go back into the soil when the rest of the plant material senesces
									AVSON[i][j] -= sumRetransN;

									// Senesce biomass - the N not retranslocated goes back to the SON pool
									double sensumt = 0;

									double* growthamount = new double[numFractions + 1];
									double* balanceamount = new double[numFractions + 1];
									double* senescamount = new double[numFractions + 1];
									for (int k = 1; k <= numfractions; k++)
									{
										int use = (int)LUSE[k][i][j];
										// determine the nitrogen senesced from each plant type
										// biomass x proportion senesced divided by biomass : N ratio
										double sensum = LUtable[use].GetArcVegSenesce() * ARCVEGBIOMASS[k][i][j] / LUtable[use].GetArcVegBNratio();
										sensumt += sensum;
										//									ARCVEGBIOMASS[k][i][j] *= (1 - LUtable[use].GetArcVegSenesce()); //reduce the biomass by the senesced amount
										senescamount[k] = ARCVEGBIOMASS[k][i][j] * LUtable[use].GetArcVegSenesce();
										growthamount[k] = ARCVEGBIOMASS[k][i][j] - OldBioMassActYear[k];
										balanceamount[k] = growthamount[k] - senescamount[k];
										ARCVEGBIOMASS[k][i][j] -= senescamount[k];
									}
									delete[] senescamount;
									delete[] growthamount;
									delete[] balanceamount;
									AVSON[i][j] += sensumt; //add the total senesced nitrogen back to the soil nitrogen pool
									


									// 7. now change sol thickness according to the biomass of moss
									// a) get the land use table entry which is used as soil (i.e. moss) from the Multilayer landuse table or from another global value (not yet defined) or from a parameter set in the ArcVeg section of the control file
									double dens = this->AVSoilDryDensity;
									int numLayers = this->AVLayerCnt;
									int luseID = this->AVLUSEType;
									// b) make the first few layers (number should be defined in the multilayer landuse table as well) thicker or thinner according to the old and ne biomass
									double actBiomassChange = ARCVEGBIOMASS[luseID][i][j] - OldBioMassActYear[luseID];
									double totalThicknesToChange = actBiomassChange / dens; // result is change in mm. Example: when biomass changed by 100g and density is 10 Kg/m^3, the result is (0.1kg/m^2)/(10 kg/m^3) = 0.01 m. Since Biomass is given in g/m^2 the calculation is 100/10 = 10, which is in mm instead of m
									double totalThicknesToChangeM = totalThicknesToChange / 1000.0; // result in m
									if ((totalThicknesToChangeM > 0.0)||(1==1)) // actually since 1==1 is true, all changes take effect
									{
										double thicksum = 0;
										// actual sum of all layers that must be changed in thickness
										for (int k = numLayers; k >= 1; k--)
										{
											thicksum += hgeo[i][j][k - 1] - hgeo[i][j][k];
										}

										double* newThickness = new double[numLayers + 1];
										// first compute all layers new thickness in a loop from bottom to top
										for (int k = numLayers; k >= 1; k--)
										{
											double actThickness = hgeo[i][j][k - 1] - hgeo[i][j][k];
											double frac = actThickness / thicksum;
											// apply thiy to the thicknesschange
											double actThickChange = frac * totalThicknesToChangeM;
											newThickness[k] = actThickness + actThickChange;
										}
										// now apply the actual thickness change to the layers upper elevation
										for (int k = numLayers; k >= 1; k--)
										{
											hgeo[i][j][k - 1] = hgeo[i][j][k] + newThickness[k];
										}
										gridinfo[2].pointer[i][j] = hgeo[i][j][0]; // actualize DHM grid 

										if (do_LayeredSnow)
										{
											// increase hgeoSnow for all layers
											for (int k = 0; k <= levelsSnow; k++)
											{
												hgeoSnow[i][j][k] += totalThicknesToChangeM;
											}
										}

										delete[] newThickness;
									}
									delete[] OldBioMassActYear;
								}
								// c) delete OldBioMass vector (must be created at the beginning ogf the runArcVeg-routine for each cell with the biomasses of each plant layer or at least the one for moss according to the entry in the multilayer landuse table
								delete[] OldBioMass;
							}
						}
					}

				}
				// TEST ONLY -------------------------------------------------------------------------------------- TEST ONLY
				// test implementation: make uppermost soil layers thinner by 50%
//				for (int i = mpistartrow; i <= mpiendrow; i++)
//				{
//					// initialize all cells from first to last in this row with 0. 
//					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
//					{
//						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
//						{
//							int s = (int)SOIL[i][j];
//							bool isArcVegType = ST[s].getIsArcVegType();
//							if (isArcVegType)
//							{
//								double thick3 = (hgeo[i][j][2] - hgeo[i][j][3]) * 0.5;
//								if (thick3 < ST[s].getDzui()[3] * 0.7)
//								{
//									thick3 = ST[s].getDzui()[3] * 0.7;
//								}
//								double thick2 = (hgeo[i][j][1] - hgeo[i][j][2]) * 0.5;
//								if (thick2 < ST[s].getDzui()[2] * 0.7)
//								{
//									thick2 = ST[s].getDzui()[2] * 0.7;
//								}
//								double thick1 = (hgeo[i][j][0] - hgeo[i][j][1]) * 0.5;
//								if (thick1 < ST[s].getDzui()[1] * 0.7)
//								{
//									thick1 = ST[s].getDzui()[1] * 0.7;
//								}
//								hgeo[i][j][2] = hgeo[i][j][3] + thick3;
//								hgeo[i][j][1] = hgeo[i][j][2] + thick2;
//								hgeo[i][j][0] = hgeo[i][j][1] + thick1;
//								gridinfo[2].pointer[i][j] = hgeo[i][j][0]; // actualize DHM grid 
//							}
//						}
//					}
//				}
			}
		}
	}
}

//#pragma optimize ("p",on)
void unsatzonmodellklasse::post_run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
	double els_fact;

	if (do_ArcVeg)
	{
		run_arcveg(me, se, tg, mt, jr);
	}

	// gathering all re-infiltrated amounts from all ranks for further processing in routing model
	MPIGatherRiverFlowTable();

//	QDsum/=count;

	for (int z=0;z<=tracer_count;z++)
	{
		//	(*WSstream) << "gather flow time zones tracer "<<z<<"[end-entry]";
		for (int i=1;i<=zonmax;i++)
		{
			//		(*WSstream) << "gather flow time zones from zone "<<i<<" with maxzeit = "<<(int)maxzeit[i]<<"[end-entry]";
			MPIGatherZone(zoneFZT_QD[z][i],1,true);
			MPIGatherZone(zoneFZT_QI[z][i],1,true);
		}
	}
	// (*WSstream) << "gather zoneQB from "<<i<<" zones"<<"[end-entry]";
	MPIGatherZone(zoneQB,zonmax,true);
//	MPIGatherZone(zoneQI,zonmax,true); not neccessary, since zoneQD and zoneQI will be initialized hereafter from zoneFZT_QD and zoneFZT_QI, resp.
//	MPIGatherZone(zoneQD,zonmax,true);

	/* discharge components for subbasins */
	if (mpi_rank==0)
	{
		for (int i=1;i<=zonmax;i++)
		{
//			(*WSstream) << "DEBUG: zoneFZT_QD[0]["<<i<<"][1].sx: "<<zoneFZT_QD[0][i][1].sx<<" zoneFZT_QD[0]["<<i<<"][1].RelCellSum: "<<zoneFZT_QD[0][i][1].RelCellSum<<"[end-entry]";

			// (*WSstream) << "calculating runoff components for zone "<<i<<"[end-entry]";
			double t1,t2,t3;
			//    QD[0][i]  = zoneFZT_QD[0][i][1].sx/zoneFZT_QD[0][i][1].RelCellSum*timestep; // FZT[i][1] is defined, always 
			QD[0][i]  = zoneFZT_QD[0][i][1].sx/zoneFZT_QD[0][i][1].RelCellSum; // FZT[i][1] is defined, always 

//(*WSstream) << "DEBUG: MPI["<<mpi_rank<<"] QD[0]["<<i<<"] "<<QD[0][i]<<" sx:"<<zoneFZT_QD[0][i][1].sx<<" RelCellSum: "<<zoneFZT_QD[0][i][1].RelCellSum<<"[end-entry]";

			for (int z=1;z<=tracer_count;z++)
			{
				QD[z][i]  = zoneFZT_QD[z][i][1].RelCellSum>0 ? zoneFZT_QD[z][i][1].sx/zoneFZT_QD[z][i][1].RelCellSum : 0;
			}
			// when surface routing is done, all cells were already counted in zoneFZT_QD[0][i][1]; 
			// now, the SLR has to be disabled by setting the SLR-coefficient to 0
			if (m_useSurfaceRouting==false)	els_fact=exp(-InV/kelsqd[i]); else els_fact=0; // 
			
			t3=QD[0][i]*(1-els_fact) + QDd[0][i]*els_fact;
			for (int z=1;z<=tracer_count;z++)
			{
				QDd[z][i] = QD[z][i]*QD[0][i]*(1-els_fact) + QDd[z][i]*QDd[0][i]*els_fact; /* adding tracer masses for influx and storage outflow */
				QDd[z][i] = t3>0 ? QDd[z][i]/t3:0.0; /* konverting masses to concentrations */
			}
			QDd[0][i] = t3; /* single linear storage for surface runoff */
//			QDstat->wertzufuegen(0,0,i,QDd[0][i]);
			zoneQD[i].n = zoneQD[i].rel;  // initializing statistic register with this mean value, weighed at the relative part of this subbasin on the total basin
			zoneQD[i].RelCellSum = zoneQD[i].rel;
			zoneQD[i].sx = QDd[0][i]*zoneQD[i].rel;

			QBb[0][i]  = zoneQB[i].sx/zoneQB[i].RelCellSum;   // mean baseflow for a subbasin 
			if (QBb[0][i]<0) QBb[0][i]=0;
			zoneQB[i].n = zoneQB[i].rel;                               /* initializing statistic register with this mean value */
			zoneQB[i].RelCellSum = zoneQB[i].rel;
			zoneQB[i].sx = QBb[0][i]*zoneQB[i].rel;

			QI[0][i]  = zoneFZT_QI[0][i][1].sx/zoneFZT_QI[0][i][1].RelCellSum;    // average subbasin value of interflow is calculated to add it to the baseflow for getting later total flow
			for (int z=1;z<=tracer_count;z++)
			{
				QI[z][i]  = zoneFZT_QI[z][i][1].RelCellSum>0 ? zoneFZT_QI[z][i][1].sx/zoneFZT_QI[z][i][1].RelCellSum : 0.0;
			}
			// when surface routing is done, all cells were already counted in zoneFZT_QD[0][i][1]; 
			// now, the SLR has to be disabled by setting the SLR-coefficient to 0
			if (m_useSurfaceRouting==false)	els_fact=exp(-InV/kelsqi[i]); else els_fact=0; // 
			// old    els_fact=exp(-InV/kelsqi[i]);
			t2 = QI[0][i]*(1-els_fact) + QIi[0][i]*els_fact; /* single linear storage effect for interflow */
			for (int z=1;z<=tracer_count;z++)
			{
				QIi[z][i] = QI[z][i]*QI[0][i]*(1-els_fact) + QIi[z][i]*QIi[0][i]*els_fact; /* adding tracer masses for influx and storage outflow */
				QIi[z][i] = t2>0 ? QIi[z][i]/t2:0.0; /* konverting masses to concentrations */
			}
			QIi[0][i] = t2; /* single linear storage effect for interflow */
			zoneQI[i].sx = QIi[0][i]*zoneQI[i].rel;
			zoneQI[i].n  = zoneQI[i].rel;
			zoneQI[i].RelCellSum  = zoneQI[i].rel;

			if (Qg[0])
			{
				if (glacierlist)
				{
					/* if there are some glaciers */
					if (glacierlist[i].isGlacierized == 0) /* else: the total discharge is already in the Qg-array (done within snow model) */
					{
						Qg[0][i] = QDd[0][i] + QIi[0][i] + QBb[0][i];
					}
					else
					{
						// here we have to add the glacier runoff to the unsaturated zone runoff with proper weights (depending on the fraction of the area for the actual subbasin)$
						double glacierrunoff = Qg[0][i]; // this is the value from the glacier model (formerly, each glacier had have its own subbasin without any unglacierized part)
						// it's measure is in mm with regard to the entire subbasin the glacier is located in (in snowmodel we already built the statistics based on the entire subbasin)
						double nonGlacierrunoff = QDd[0][i] + QIi[0][i] + QBb[0][i]; // that's the same thing as for subcatchments without glaciers --> here, tto, the runoff is weighted on the total area of the subbasin
						// now we have to calculate simply the sum of both components
//if ((tg == 1)&&(i==31))
//{
//	cout << jr << mt << tg << se << me << " glacierrunoff " << glacierrunoff << " nonGlacierrunoff "<<nonGlacierrunoff <<endl;
//
//}
						Qg[0][i] = glacierrunoff + nonGlacierrunoff;

					}
				}
				else
				{ /* if there are no glaciers */
					Qg[0][i] = QDd[0][i] + QIi[0][i] + QBb[0][i];
				}
			}
			else
			{
				(*WSstream) << "FATAL: Qg["<<i<<"] not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
			}
			// Qgstat->wertzufuegen(0,0,i,Qg[0][i]);
			zoneQg[i].sx		 = Qg[0][i]*zoneQg[i].rel;
			zoneQg[i].n			 = zoneQg[i].rel;
			zoneQg[i].RelCellSum = zoneQg[i].rel;
		}
		/* writing theta, conc and qbottom-values for the control call into the file buffers */
	}
	//  (*WSstream) << "calling statistics in unsatzonmodel"<<"[end-entry]";

	/* calculating zonal statistics for the output files */
	GWSTANDstat->calc_and_done(jahr,monat,tag,stunde, minute);
	GWNstat->calc_and_done(jahr,monat,tag,stunde, minute);
	INFILTstat->calc_and_done(jahr,monat,tag,stunde, minute);
	MACROstat->calc_and_done(jahr,monat,tag,stunde, minute);
	if (methode>2) SUBSTEPSstat->calc_and_done(jahr,monat,tag,stunde, minute);
	SB_1_stat->calc_and_done(jahr,monat,tag,stunde,minute);
	SB_2_stat->calc_and_done(jahr,monat,tag,stunde,minute);
	if (ICEOLD)	SB_2_IMMO_stat->calc_and_done(jahr,monat,tag,stunde, minute);
	ROOTstat->calc_and_done(jahr,monat,tag,stunde, minute);
	INFEXstat->calc_and_done(jahr,monat,tag,stunde, minute);
	QDstat->calc_and_done(jahr,monat,tag,stunde,minute,0);
	QIstat->calc_and_done(jahr,monat,tag,stunde,minute,0);
	QBstat->calc_and_done(jahr,monat,tag,stunde,minute,0);
	if (gwflow_active) GWINstat->calc_and_done(jahr,monat,tag,stunde,minute);
	if (gwflow_active) GWEXstat->calc_and_done(jahr,monat,tag,stunde,minute);
	Qgstat->calc_and_done(jahr,monat,tag,stunde, minute,0);
	if (HDRAIN) QDRAINstat->calc_and_done(jahr,monat,tag,stunde, minute);
	if (MAXPOND) PONDstat->calc_and_done(jahr,monat,tag,stunde, minute);
	if (do_TMod>0) 
	{
		if (this->iTModMethod>0) 
		{
			thawdepth_stat->calc_and_done(jahr,monat,tag,stunde, minute);
		}
	}

#ifdef _WASIMDEBUG_
	iteration_stat->calc_and_done(jahr,monat,tag,stunde, minute);
#endif
	/* writing grids and stacks (depending on write code) to the disc */
	writegridtofile(sb_1_index);
	writegridtofile(sb_2_index);
	writegridtofile(rootindex);
	writestacktofile(theta_index);
	writestacktofile(hh_index);
	writestacktofile(hgeo_index);
	writestacktofile(qu_index);
	writegridtofile(gwstandindex);
	writegridtofile(gwthetaindex);
	writegridtofile(gwnindex);
	writegridtofile(macroindex);
	//  writegridtofile(irrigindex); will be done in irrigationmodel now
	if (methode>2) writegridtofile(substepsindex);
	writegridtofile(upriseindex);
	writegridtofile(percolindex);
	writegridtofile(gwlevelindex);
	writegridtofile(infexindex);
	writegridtofile(sattindex);
	writegridtofile(qdindex);
	writegridtofile(qiindex);
	writegridtofile(deltagwindex);
	if (gwflow_active) writegridtofile(gwinindex);
	if (gwflow_active) writegridtofile(gwexindex);
	if (MAXPOND) writegridtofile(pondindex);
	if (HDRAIN) writegridtofile(qdrainindex);
	if (do_TMod>0)
	{
		if (this->iTModMethod>0) 
		{
			writestacktofile(temperindex);
			writestacktofile(seindex);
			writegridtofile(thawtmodindex);
		}
	}

	if (do_LayeredSnow)
	{
		writegridtofile(upperIndexSnowindex);
		writestacktofile(swesindex);
		writestacktofile(swelindex);
		writestacktofile(thetaSnowindex);
		writestacktofile(hhydSnowindex);
		writestacktofile(hgeoSnowindex);
		writestacktofile(quSnowindex);
		writestacktofile(thickSnowindex);
		writestacktofile(densSnowindex);
		if (this->iTModMethod>0)
		{
			writestacktofile(temperSnowindex);
		}
	}

	if (do_ArcVeg)
	{
		writestacktofile(ExcessIceFactorindex);
		writestacktofile(AVINHLPindex);
		writegridtofile(AVGDDindex);
		writegridtofile(AVINMAX7DAVGindex);
		writegridtofile(AVTDAVGindex);
		writegridtofile(AVINDAVGindex);
		writegridtofile(AVSNOWFREEindex);
		writegridtofile(AVPGPindex);
		writestacktofile(AVINSCALindex);
		writestacktofile(AVNMINindex);
		writegridtofile(AVNMINTOTindex);
		writegridtofile(AVSONindex);
	}

	updateFieldCapacityGrids();
	updateSB1RELGrid();

	for (int z=1;z<=tracer_count;z++)
	{
		/* only one grid of the grid stack is considered for a zonal statistic.
The following loop copies the values of the level of interest into the grid */
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				int k=1; /* concentration of the first layer */
				//        k=(int)GWLEVEL[i][j]; /* concentration of the first partially saturated layer */
				/* level of interest: first level (for level of groundwater: ...= conc[z][i][j][(int)GWLEVEL[i][j]]) */
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
				//        if (zn[i][j]!=nodata) 
				{
					if (tracer_info[z].evaporating)
					{
						conc1[z][i][j]=conc[z][i][j][k];
					}
					else
					{
//						conc1[z][i][j] = theta[i][j][k] > 0.0001 ? conc[z][i][j][k] / (theta[i][j][k] * 1000.0*ST[(int)SOIL[i][j]].getDzui()[k]) : 99999.0;
	// ToDo: dicken should be initialized again here
						conc1[z][i][j] = theta[i][j][k] > 0.0001 ? conc[z][i][j][k] / (theta[i][j][k] * 1000.0*ST[(int)SOIL[i][j]].getDzui()[k]) : 99999.0;
						//            if ((conc1[z][i][j]>tracer_info[z].max_conc)&&(tracer_info[z].max_conc>0))
						//              conc1[z][i][j]=tracer_info[z].max_conc;
					}
				}
			}  
		}

		writestacktofile(tracer_module[z].storindex);
		/* mean zonal concentrations now have to be copied into the Qg-Array */
		for (int i=1;i<=zonmax;i++)
		{
			/* tracer concentrations in baseflow */
			if (tracer_module[z].zone[i].RelCellSum>0)
			{
				QBb[z][i]=tracer_module[z].zone[i].sx/tracer_module[z].zone[i].RelCellSum;
			}
			else
			{
				QBb[z][i]=0;
			}
			/* superposition weighed by the mean fractions of the runoff components for subbasins */
			if ((QBb[0][i]+QIi[0][i]+QDd[0][i])>0)
			{
				Qg[z][i] = (QBb[0][i]*QBb[z][i] + QIi[0][i]*QIi[z][i] + QDd[0][i]*QDd[z][i])/(QBb[0][i]+QIi[0][i]+QDd[0][i]);
			}
			else
			{
				Qg[z][i] = 0; /* if non discharge at all: all konzentrations are set to zero (if no or only negative baseflow 
					flows [and no other components flow]) */
			}
			tracer_module[z].zone[i].sx=Qg[z][i];
			tracer_module[z].zone[i].n=1;
			tracer_module[z].zone[i].RelCellSum=1;
		}
		/* concentration is the mixed concentration of direct flow, interflow and baseflow (exfiltration or deltagw) */
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde, minute);
	}

}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::ClearMemory(bool appendmode)
/* this function frees memeory for all the allocated variables */
{
	int j,z,k;
	unsigned long i;
	double *dicken;
	if (do_it > 0)
	{
	if (MAXPOND) writegridtofile(pondindex); // since POND could have been changed in routing model (Lake model changes content of reservoirs), the POND grid must be written again

		if (RecvBuffer)
		{
			for (int i = 0; i < mpi_size; i++)
			{
				if ((*RivLnkCellsToReceive)[i].size()>0)
				{
					if (RecvBuffer[i]) delete[] RecvBuffer[i];
				}
			}
			if (RecvBuffer)		delete[] RecvBuffer;
			if (RecvBufferSize) delete[] RecvBufferSize;
			if (RecvBufferPos)	delete[] RecvBufferPos;
		}
		
		if (SendBuffer)
		{
			for (int i = 0; i < RivLnkCellsToSend->size(); i++)
			{
				if ((*RivLnkCellsToSend)[i].size()>0)
				{
					if (SendBuffer[i]) delete[] SendBuffer[i];
				}
			}
			if (SendBuffer) delete[] SendBuffer;
		}

		if (tempGridRows) delete[] tempGridRows;
		if (tempStackRows) delete[] tempStackRows;

		if (initialized) WriteSlrContents(appendmode);
		if (ICE) { if (glacierlist) free(glacierlist);}
		// if (initialized) levels=ST[0].levels; else levels=0;
		levels=ST[0].getLevels();
		//    if (balance) free(balance);
		//    if (maxbalance) free(maxbalance);
		double deltasb=0;
		double deltapond=0;
		if (initialized) 
		{
			if (control_row>0)
			{
				if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
				{
					dicken=ST[(int)SOIL[control_row][control_col]].getDzui();
					for (int k=1;k<=levels;k++)
					deltasb+=(theta[control_row][control_col][k]*dicken[k]*1000.0);
					char logfilename[1024];
					sprintf(logfilename,"modell_%d.log",mpi_rank);
					logdatei=fopen(logfilename,"a+");
					fprintf(logdatei,"(mpi_rank %d) soil water content (control point r=%d c=%d) at model-end: %lf mm\n",mpi_rank,control_row,control_col,deltasb);
					fclose(logdatei);
				}
			}
			double sum_sb1=0.0;
			double anz_sb1=0.0;
			double sum_sb2=0.0;
			double anz_sb2=0.0;
			double *thetasat;
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
					if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
					else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
					if ((zn[i][j]>0)&&(MPIREG[i][j]==mpi_rank))
					{
						int s=(int)SOIL[i][j];
						dicken=ST[s].getDzui();
						thetasat = ST[s].getThetas();
						for (int k=1;k<=levels;k++)
						{
							sum_sb1+=(theta[i][j][k]*dicken[k])*(1.0-icevalue);
							sum_sb2+=(thetasat[k]*dicken[k])*icevalue;
						}
						anz_sb1+=RelCellArea!=NULL?RelCellArea[i][j]*(1.0-icevalue):(1.0-icevalue);
						anz_sb1+=icevalue;
						anz_sb2+=RelCellArea!=NULL?RelCellArea[i][j]*icevalue:icevalue;
						anz_sb2+=(1.0-icevalue);
					}
				}
			}
			double sumtemp1=sum_sb1;
			double anztemp1=anz_sb1;
			double sumtemp2=sum_sb2;
			double anztemp2=anz_sb2;
#ifdef _MPI
			// reduce anz_sb and sum_sb from all instances
			MPI_Reduce(&anz_sb1, &anztemp1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&sum_sb1, &sumtemp1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&anz_sb2, &anztemp2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&sum_sb2, &sumtemp2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
			if (mpi_rank==0)
			{
				char logfilename[1024];
				sprintf(logfilename,"modell_%d.log",mpi_rank);
				logdatei=fopen(logfilename,"a+");
				fprintf(logdatei,"mean soil water content for all cells at model-end: %lf m^3/m^2\n",sumtemp1/anztemp1+sumtemp2/anztemp2);
				fprintf(logdatei,"mean soil water content for unglacierized or only partly unglacierized cells at model-end: %lf m^3/m^2\n",sumtemp1/anztemp1);
				fprintf(logdatei,"mean soil water content for partly or completely glacierized cells at model-end: %lf m^3/m^2\n",sumtemp2/anztemp2);
				fclose(logdatei);
			}
		}
		substance_transport_delete();
		if (additionals) free(additionals);

		/* if there are  unsaved data in the buffers, they have to be saved now */
		if (control_row>0)
		{
			if ((MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
			{
				if (initialized) flushbuffer(   qudatei,   qudateiname,   &qubufcnt,   qubuffer);
				if (initialized) flushbuffer( thetadatei,thetadateiname,&thetabufcnt,thetabuffer);
				if (initialized) flushbuffer( hhyddatei, hhyddateiname,& hhydbufcnt, hhydbuffer);
				if (initialized) flushbuffer( otherdatadatei, otherdatadateiname,&otherdatabufcnt,otherdatabuffer,1,16+2+this->iNumVegetationLayers*6);
				if (initialized) flushbuffer( etrdetaildatei, etrdetaildateiname,& etrdetailbufcnt, etrdetailbuffer,0);
				if (initialized) flushbuffer( interflowdatei, interflowdateiname,& interflowbufcnt, interflowbuffer);

				for (int i=1;i<=qubufmax;i++) {if (qubuffer[i].wert) free(qubuffer[i].wert);}
				if (qubuffer) free(qubuffer);

				for (int i=1;i<=thetabufmax;i++) {if (thetabuffer[i].wert) free(thetabuffer[i].wert);}
				if (thetabuffer) free(thetabuffer);

				for (int i=1;i<=hhydbufmax;i++) {if (hhydbuffer[i].wert) free(hhydbuffer[i].wert);}
				if (hhydbuffer) free(hhydbuffer);

				for (int i=1;i<=otherdatabufmax;i++) {if (otherdatabuffer[i].wert) free(otherdatabuffer[i].wert);}
				if (otherdatabuffer) free(otherdatabuffer);
				
				for (int i=1;i<=etrdetailbufmax;i++) {if (etrdetailbuffer[i].wert) free(etrdetailbuffer[i].wert);}
				if (etrdetailbuffer) free(etrdetailbuffer);

				for (int i=1;i<=interflowbufmax;i++) {if (interflowbuffer[i].wert) free(interflowbuffer[i].wert);}
				if (interflowbuffer) free(interflowbuffer);

				if (initialized && do_TMod>0 && iTModMethod>0)
				{
					flushbuffer(temperfile,temperfilename,&temperbufcnt,temperbuffer);
		
					for (int i=1;i<=temperbufmax;i++) {if (temperbuffer[i].wert) free(temperbuffer[i].wert);}
					if (temperbuffer) free(temperbuffer);

				}
			}
		}

		INFEXstat->ClearMemory();   if (INFEXstat)   delete(INFEXstat   ); INFEXstat   = NULL;
		QDstat->ClearMemory();      if (QDstat)      delete(QDstat      ); QDstat      = NULL;
		QIstat->ClearMemory();      if (QIstat)      delete(QIstat      ); QIstat      = NULL;
		QBstat->ClearMemory();      if (QBstat)      delete(QBstat      ); QBstat      = NULL;
		if (gwflow_active) 
		{ GWINstat->ClearMemory();    if (GWINstat)  delete(GWINstat    ); GWINstat    = NULL;}
		if (gwflow_active) 
		{ GWEXstat->ClearMemory();  if (GWEXstat)    delete(GWEXstat    ); GWEXstat    = NULL;}
		Qgstat->ClearMemory();      if (Qgstat)      delete(Qgstat      ); Qgstat      = NULL;
		SB_1_stat->ClearMemory();   if (SB_1_stat)   delete(SB_1_stat   ); SB_1_stat   = NULL; 
		SB_2_stat->ClearMemory();   if (SB_2_stat)   delete(SB_2_stat   ); SB_2_stat   = NULL;
		if (ICEOLD)	SB_2_IMMO_stat->ClearMemory();   if (SB_2_IMMO_stat)   delete(SB_2_IMMO_stat   ); SB_2_IMMO_stat   = NULL;
		ROOTstat->ClearMemory();    if (ROOTstat)    delete(ROOTstat    ); ROOTstat    = NULL;
		GWSTANDstat->ClearMemory(); if (GWSTANDstat) delete(GWSTANDstat ); GWSTANDstat = NULL;
		GWNstat->ClearMemory();     if (GWNstat)     delete(GWNstat     ); GWNstat     = NULL;
		INFILTstat->ClearMemory();  if (INFILTstat)  delete(INFILTstat  ); INFILTstat  = NULL;
		MACROstat->ClearMemory();   if (MACROstat)   delete(MACROstat   ); MACROstat   = NULL;
		if (methode>2)
		{
			SUBSTEPSstat->ClearMemory();   if (SUBSTEPSstat)   delete(SUBSTEPSstat   ); SUBSTEPSstat   = NULL;
		}
		if (HDRAIN)  {QDRAINstat->ClearMemory(); if (QDRAINstat) delete(QDRAINstat  ); QDRAINstat  = NULL;}
		if (MAXPOND) {PONDstat->ClearMemory();   if (PONDstat)   delete(PONDstat    ); PONDstat    = NULL;}
		
		if (do_TMod>0 && iTModMethod>0)
		{
			if (thawdepth_stat)
			{
				thawdepth_stat->ClearMemory(); 
				delete(thawdepth_stat);
				thawdepth_stat=NULL;
			}
		}

#ifdef _WASIMDEBUG_
		iteration_stat->ClearMemory(); if (iteration_stat) delete(iteration_stat); iteration_stat=NULL;
		fclose(tempdatei);
#endif
		if (QDstatfilename     ) free(QDstatfilename     ); QDstatfilename       = NULL;
		if (QIstatfilename     ) free(QIstatfilename     ); QIstatfilename       = NULL;
		if (QBstatfilename     ) free(QBstatfilename     ); QBstatfilename       = NULL;
		if (GWINstatfilename   ) free(GWINstatfilename   ); GWINstatfilename     = NULL;
		if (GWEXstatfilename   ) free(GWEXstatfilename   ); GWEXstatfilename     = NULL;
		if (Qgstatfilename     ) free(Qgstatfilename     ); Qgstatfilename       = NULL;
		if (SB_1_statfilename  ) free(SB_1_statfilename  ); SB_1_statfilename    = NULL;
		if (SB_2_statfilename  ) free(SB_2_statfilename  ); SB_2_statfilename    = NULL;
		if (ICEOLD)	if (SB_2_IMMO_statfilename  ) free(SB_2_IMMO_statfilename  ); SB_2_IMMO_statfilename    = NULL;
		if (ROOTstatfilename   ) free(ROOTstatfilename   ); ROOTstatfilename     = NULL;
		if (QDRAINstatfilename ) free(QDRAINstatfilename ); QDRAINstatfilename   = NULL;
		if (GWSTANDstatfilename) free(GWSTANDstatfilename); GWSTANDstatfilename  = NULL;
		if (GWNstatfilename    ) free(GWNstatfilename    ); GWNstatfilename      = NULL;
		if (INFILTstatfilename ) free(INFILTstatfilename ); INFILTstatfilename   = NULL;
		if (MACROstatfilename  ) free(MACROstatfilename  ); MACROstatfilename    = NULL;
		if (methode>2)
		{
			if (SUBSTEPSstatfilename  ) free(SUBSTEPSstatfilename  ); SUBSTEPSstatfilename    = NULL;
		}
		if (PONDstatfilename   ) free(PONDstatfilename   ); PONDstatfilename     = NULL;
		if (INFEXstatfilename  ) free(INFEXstatfilename  ); INFEXstatfilename    = NULL;
		
		if (do_TMod>0 && this->iTModMethod>0)
		{
			if (thawdepth_statfilename) free(thawdepth_statfilename); thawdepth_statfilename  = NULL;
		}

#ifdef _WASIMDEBUG_
		if (iteration_statfilename) free(iteration_statfilename); iteration_statfilename= NULL;
#endif

		delete_grid( theta_index);
		delete_grid(    hh_index);
		delete_grid(  hgeo_index);
		delete_grid(    qu_index);
		delete_grid(gwstandindex);
                delete_grid(throotindex);
		delete_grid(sb_1_index);
		delete_grid(sb_2_index);
		delete_grid(rootindex);
		if (irrigationtable)
		{
			delete_grid(irrigindex);
			delete_grid(maxpumpindex);
			delete_grid(sblowindex);
			delete_grid(sbhighindex);
			if (etr2etpindex>0)   
			delete_grid(etr2etpindex);
		}
		delete_grid(balanceindex);
		delete_grid(gwtableexternindex);
		delete_grid(gwthetaindex);
		delete_grid(    gwnindex);
		delete_grid( percolindex);
		delete_grid( macroindex);
		if (methode>2)
		{
			delete_grid( substepsindex);
		}
		delete_grid( upriseindex);
		delete_grid(tanbetaindex);
		delete_grid(gwlevelindex);
		delete_grid(     qdindex);
		delete_grid(     qiindex);
		delete_grid(deltagwindex);
//		delete_grid(olddeltagwindex);
		if (gwflow_active) 
		{
			delete_grid(   gwinindex);
			delete_grid(   gwexindex);
		}
		delete_grid(   sattindex);
		delete_grid( threduindex);
		delete_grid(  infexindex);
		delete_grid(  infiltindex);
		delete_grid(qu1extindex);

		if (MAXPOND) delete_grid(pondindex);
		if (HDRAIN) { delete_grid(qdrainindex); delete_grid(drainlevelindex); }
		if (do_TMod>0 && iTModMethod>0)
		{
			delete_grid(temperindex);
			if (tlowbndindex>0) delete_grid(tlowbndindex);
			delete_grid(thawtmodindex);
		}
		if (do_TMod>0)
		{
			if (FROZENDEPTH) 
			{
				delete_grid(frozendepthingwindex);
				free_matrix(FROZENDEPTH, gridinfo[frozendepthingwindex].nrl, gridinfo[frozendepthingwindex].nrh, 1, ncols,nrows);
				FROZENDEPTH=NULL;
			}
		}

		free_drei_d_matrix(theta, gridinfo[theta_index].nrl, gridinfo[theta_index].nrh, 1, ncols, levels, nrows);      theta = NULL;
		free_drei_d_matrix(hh, gridinfo[hh_index].nrl, gridinfo[hh_index].nrh, 1, ncols, levels,nrows);         hh = NULL;
		free_drei_d_matrix(hgeo, gridinfo[hgeo_index].nrl, gridinfo[hgeo_index].nrh, 1, ncols, levels,nrows);       hgeo = NULL;
		free_drei_d_matrix(qu, gridinfo[qu_index].nrl, gridinfo[qu_index].nrh, 1, ncols, levels,nrows);         qu = NULL;

		free_matrix(theta1, gridinfo[theta_index].nrl, gridinfo[theta_index].nrh, 1, ncols,nrows); theta1 = NULL;
		free_matrix(hh1, gridinfo[hh_index].nrl, gridinfo[hh_index].nrh, 1, ncols, nrows);    hh1 = NULL;
		free_matrix(hgeo1, gridinfo[hgeo_index].nrl, gridinfo[hgeo_index].nrh, 1, ncols, nrows);  hgeo1 = NULL;
		free_matrix(qu1, gridinfo[qu_index].nrl, gridinfo[qu_index].nrh, 1, ncols, nrows);    qu1 = NULL;
		free_matrix(qu1_ext, gridinfo[qu1extindex].nrl, gridinfo[qu1extindex].nrh, 1, ncols, nrows);    qu1_ext = NULL;

		free_matrix(GWSTAND, gridinfo[gwstandindex].nrl, gridinfo[gwstandindex].nrh, 1, ncols, nrows);    GWSTAND = NULL;
		free_matrix(GWtheta, gridinfo[gwthetaindex].nrl, gridinfo[gwthetaindex].nrh, 1, ncols, nrows);    GWtheta = NULL;
		free_matrix(GWN, gridinfo[gwnindex].nrl, gridinfo[gwnindex].nrh, 1, ncols, nrows);        GWN = NULL;
		free_matrix(MACRO, gridinfo[macroindex].nrl, gridinfo[macroindex].nrh, 1, ncols, nrows);      MACRO = NULL;
		if (methode>2)
		{
			free_matrix(SUBSTEPS, gridinfo[substepsindex].nrl, gridinfo[substepsindex].nrh, 1, ncols, nrows);      SUBSTEPS = NULL;
		}
		
		free_matrix(UPRISE, gridinfo[upriseindex].nrl, gridinfo[upriseindex].nrh, 1, ncols, nrows);     UPRISE = NULL;
		free_matrix(PERCOL, gridinfo[percolindex].nrl, gridinfo[percolindex].nrh, 1, ncols, nrows);     PERCOL = NULL;
		free_matrix(GWLEVEL, gridinfo[gwlevelindex].nrl, gridinfo[gwlevelindex].nrh, 1, ncols, nrows);    GWLEVEL = NULL;
		free_matrix(TANBETA, gridinfo[tanbetaindex].nrl, gridinfo[tanbetaindex].nrh, 1, ncols, nrows);    TANBETA = NULL;
		//    free_matrix(THREDU    ,gridinfo[threduindex].nrl, gridinfo[threduindex].nrh, 1, ncols,nrows);     THREDU = NULL;
		free_matrix(GWTableExtern, gridinfo[gwtableexternindex].nrl, gridinfo[gwtableexternindex].nrh, 1, ncols, nrows); GWTableExtern = NULL;
		free_matrix(EXTBALANCE, gridinfo[balanceindex ].nrl, gridinfo[balanceindex ].nrh, 1, ncols, nrows); EXTBALANCE = NULL;
		free_matrix(SB_1_, gridinfo[sb_1_index].nrl, gridinfo[sb_1_index].nrh, 1, ncols, nrows);    SB_1_ = NULL;
		free_matrix(SB_2_, gridinfo[sb_2_index].nrl, gridinfo[sb_2_index].nrh, 1, ncols, nrows);    SB_2_ = NULL;
		free_matrix(ROOT, gridinfo[rootindex].nrl, gridinfo[rootindex].nrh, 1, ncols, nrows);    ROOT = NULL;
		free_matrix(THROOT, gridinfo[throotindex].nrl, gridinfo[throotindex].nrh, 1, ncols, nrows);    THROOT = NULL;

		if (irrigationtable) 
		{
			free_matrix(IRRIG, gridinfo[irrigindex].nrl, gridinfo[irrigindex].nrh, 1, ncols, nrows);    IRRIG = NULL;
			free_matrix(MAXPUMP, gridinfo[maxpumpindex].nrl, gridinfo[maxpumpindex].nrh, 1, ncols, nrows);    MAXPUMP = NULL;
			free_matrix(SBLOW, gridinfo[sblowindex].nrl, gridinfo[sblowindex].nrh, 1, ncols, nrows);    SBLOW = NULL;
			free_matrix(SBHIGH, gridinfo[sbhighindex].nrl, gridinfo[sbhighindex].nrh, 1, ncols, nrows);    SBHIGH = NULL;
			if (etr2etpindex>0)
			{
				free_matrix(ETR2ETP, gridinfo[etr2etpindex].nrl, gridinfo[etr2etpindex].nrh, 1, ncols, nrows);
				ETR2ETP=NULL;
			}
		}

		free_matrix(QDgrid, gridinfo[qdindex].nrl, gridinfo[qdindex].nrh, 1, ncols, nrows);     QDgrid = NULL;
		free_matrix(QIgrid, gridinfo[qiindex].nrl, gridinfo[qiindex].nrh, 1, ncols, nrows);     QIgrid = NULL;
		free_matrix(DELTAGW, gridinfo[deltagwindex].nrl, gridinfo[deltagwindex].nrh, 1, ncols, nrows);     DELTAGW = NULL;
		if (gwflow_active) 
		{
			free_matrix(GWIN, gridinfo[gwinindex].nrl, gridinfo[gwinindex].nrh, 1, ncols, nrows);       GWIN = NULL;
			free_matrix(GWEX, gridinfo[gwexindex].nrl, gridinfo[gwexindex].nrh, 1, ncols, nrows);       GWEX = NULL;
		}
//		free_matrix(oldDELTAGW ,mpistartrow, mpiendrow, 1, ncols,nrows);  oldDELTAGW = NULL;
		free_matrix(SATT, gridinfo[sattindex].nrl, gridinfo[sattindex].nrh, 1, ncols, nrows);       SATT = NULL;
		free_matrix(INFEX, gridinfo[infexindex].nrl, gridinfo[infexindex].nrh, 1, ncols, nrows);      INFEX = NULL;
		free_matrix(INFILT, gridinfo[infiltindex].nrl, gridinfo[infiltindex].nrh, 1, ncols, nrows);      INFILT = NULL;
		if (MAXPOND) { free_matrix(POND, gridinfo[pondindex].nrl, gridinfo[pondindex].nrh, 1, ncols, nrows);       POND = NULL; }
		if (HDRAIN)  {
			free_matrix(DRAINLEVEL, gridinfo[drainlevelindex].nrl, gridinfo[drainlevelindex].nrh, 1, ncols, nrows); DRAINLEVEL = NULL;
			free_matrix(QDRAIN, gridinfo[qdrainindex].nrl, gridinfo[qdrainindex].nrh, 1, ncols, nrows);     QDRAIN = NULL;
		}

		if (do_TMod>0 && iTModMethod>0)
		{
			free_drei_d_matrix(temper, gridinfo[temperindex].nrl, gridinfo[temperindex].nrh, 1, ncols, nrows, levels);      temper = NULL;
			free_matrix(temper1, gridinfo[temperindex].nrl, gridinfo[temperindex].nrh, 1, ncols, nrows);                   temper1 = NULL;

			free_drei_d_matrix(SEstack, gridinfo[seindex].nrl, gridinfo[seindex].nrh, 1, ncols, nrows, levels);            SEstack = NULL;
			free_matrix(SE1, gridinfo[seindex].nrl, gridinfo[seindex].nrh, 1, ncols, nrows);							       SE1 = NULL;

			free_matrix(THAWTMOD, gridinfo[thawtmodindex].nrl, gridinfo[thawtmodindex].nrh, 1, ncols, nrows);             THAWTMOD = NULL;
			if (tlowbndindex > 0)
			{
				free_matrix(TLOWBND, gridinfo[tlowbndindex].nrl, gridinfo[tlowbndindex].nrh, 1, ncols, nrows);             
				TLOWBND = NULL;
			}
		}



		if (Q0) free(Q0);
		//    if (QBb) free(QBb);
		//    if (QB) free(QB);
		if (flussdichte) free(flussdichte);
		if (code) free(code); code=NULL;
		if (kelsqd) free(kelsqd);  kelsqd = NULL;
		if (kelsqi) free(kelsqi);  kelsqi = NULL;
		if (kelsqlat) free(kelsqlat);  kelsqlat = NULL;
		if (SR_Anteil)      free(SR_Anteil);
		if (thetadateiname) free(thetadateiname);
		if (hhyddateiname) free(hhyddateiname);
		if (otherdatadateiname) free(otherdatadateiname);
		if (etrdetaildateiname) free(etrdetaildateiname);
		if (interflowdateiname) free(interflowdateiname);
		if (qudateiname) free(qudateiname);
		if (storagefilename) free(storagefilename);

		if (do_TMod>0 && iTModMethod>0)
		{
			if (temperfilename) free(temperfilename);
		}


		if (maxzeit) free(maxzeit);
		if ((zoneFZT_QD)&&(zoneFZT_QI))
		{
			for (z=0;z<=tracer_count;z++)
			{
				for (int j=1;j<=zonmax;j++)
				{
					if ((statrecord*)zoneFZT_QD[z][j]) free((statrecord*)zoneFZT_QD[z][j]);
					if ((statrecord*)zoneFZT_QI[z][j]) free((statrecord*)zoneFZT_QI[z][j]);
					zoneFZT_QD[z][j]=NULL;
					zoneFZT_QI[z][j]=NULL;
				}
				if ((statrecord**)zoneFZT_QD[z]) free((statrecord**)zoneFZT_QD[z]);
				if ((statrecord**)zoneFZT_QI[z]) free((statrecord**)zoneFZT_QI[z]);
				zoneFZT_QD[z]=NULL;
				zoneFZT_QI[z]=NULL;
			}
		}
		if (zoneFZT_QD) free(zoneFZT_QD);
		if (zoneFZT_QI) free(zoneFZT_QI);
		zoneFZT_QD=NULL;
		zoneFZT_QI=NULL;

		if (zarray) free(zarray); zarray=NULL; /* local newly allocated! */
		if (zone) free(zone);   zone = NULL; /* local newly allocated! */

		if (zoneINFEX) free(zoneINFEX);   zoneINFEX   = NULL;
		if (zoneQD) free(zoneQD);      zoneQD      = NULL;
		if (zoneQI) free(zoneQI);      zoneQI      = NULL;
		if (zoneQB) free(zoneQB);      zoneQB      = NULL;
		if (zoneQg) free(zoneQg);      zoneQg      = NULL;
		if (gwflow_active) 
		{
			if (zoneGWIN) free(zoneGWIN);    zoneGWIN    = NULL;
			if (zoneGWEX) free(zoneGWEX);    zoneGWEX    = NULL;
		}
		if (zoneSB_1_)   free(zoneSB_1_);    zoneSB_1_   = NULL;
		if (zoneSB_2_)   free(zoneSB_2_);    zoneSB_2_   = NULL;
		if (ICEOLD)	if (zoneSB_2_IMMO)   free(zoneSB_2_IMMO);    zoneSB_2_IMMO   = NULL;
		if (zoneROOT)    free(zoneROOT);     zoneROOT    = NULL;
		if (zoneGWSTAND) free(zoneGWSTAND);  zoneGWSTAND = NULL;
		if (zoneGWN)     free(zoneGWN);      zoneGWN     = NULL;
		if (zoneINFILT)  free(zoneINFILT);   zoneINFILT   = NULL;
		if (zoneMACRO)   free(zoneMACRO);    zoneMACRO   = NULL;
		if (methode>2)
		{
			if (zoneSUBSTEPS)   free(zoneSUBSTEPS);    zoneSUBSTEPS   = NULL;
		}
		if (HDRAIN)  {if (zoneQDRAIN) free(zoneQDRAIN);  zoneQDRAIN  = NULL;}
		if (MAXPOND) {if (zonePOND)   free(zonePOND);    zonePOND    = NULL; }
		if (do_TMod>0 && this->iTModMethod>0)
		{
			if (zonethawdepth) free(zonethawdepth); zonethawdepth = NULL;
		}
		
		if (LUSE) delete[] LUSE; LUSE = NULL;
		if (ETP ) delete[] ETP;  ETP  = NULL;
		if (ETR ) delete[] ETR;  ETR  = NULL;
		if (WURZ) delete[] WURZ; WURZ = NULL;
		if (VBG ) delete[] VBG;  VBG  = NULL;
		if (RSI ) delete[] RSI;  RSI = NULL;
		if (RSC ) delete[] RSC;  RSC = NULL;
		if (LAI ) delete[] LAI;  LAI = NULL;
		if (Z0  ) delete[] Z0 ;  Z0  = NULL;

		if (do_LayeredSnow)
		{

			free_drei_d_matrix(thetaSnow, gridinfo[thetaSnowindex].nrl, gridinfo[thetaSnowindex].nrh, 1, ncols, levelsSnow, nrows);      thetaSnow = NULL;
			free_drei_d_matrix(hhydSnow, gridinfo[hhydSnowindex].nrl, gridinfo[hhydSnowindex].nrh, 1, ncols, levelsSnow, nrows);         hhydSnow  = NULL;
			free_drei_d_matrix(hgeoSnow, gridinfo[hgeoSnowindex].nrl, gridinfo[hgeoSnowindex].nrh, 1, ncols, levelsSnow, nrows);         hgeoSnow  = NULL;
			free_drei_d_matrix(quSnow, gridinfo[quSnowindex].nrl, gridinfo[quSnowindex].nrh, 1, ncols, levelsSnow, nrows);               quSnow    = NULL;
			free_drei_d_matrix(swes, gridinfo[swesindex].nrl, gridinfo[swesindex].nrh, 1, ncols, levelsSnow, nrows);                     swes      = NULL;
			free_drei_d_matrix(swel, gridinfo[swelindex].nrl, gridinfo[swelindex].nrh, 1, ncols, levelsSnow, nrows);                     swel      = NULL;
			free_drei_d_matrix(thickSnow, gridinfo[thickSnowindex].nrl, gridinfo[thickSnowindex].nrh, 1, ncols, levelsSnow, nrows);      thickSnow = NULL;
			free_drei_d_matrix(densSnow, gridinfo[densSnowindex].nrl, gridinfo[densSnowindex].nrh, 1, ncols, levelsSnow, nrows);         densSnow  = NULL;

			free_matrix(thetaSnow1, gridinfo[thetaSnowindex].nrl, gridinfo[thetaSnowindex].nrh, 1, ncols, nrows); thetaSnow1 = NULL;
			free_matrix(hhydSnow1, gridinfo[hhydSnowindex].nrl, gridinfo[hhydSnowindex].nrh, 1, ncols, nrows);    hhydSnow1 = NULL;
			free_matrix(hgeoSnow1, gridinfo[hgeoSnowindex].nrl, gridinfo[hgeoSnowindex].nrh, 1, ncols, nrows);  hgeoSnow1 = NULL;
			free_matrix(quSnow1, gridinfo[quSnowindex].nrl, gridinfo[quSnowindex].nrh, 1, ncols, nrows);    quSnow1 = NULL;
			free_matrix(swes1, gridinfo[swesindex].nrl, gridinfo[swesindex].nrh, 1, ncols, nrows);    swes1 = NULL;
			free_matrix(swel1, gridinfo[swelindex].nrl, gridinfo[swelindex].nrh, 1, ncols, nrows);    swel1 = NULL;
			free_matrix(thickSnow1, gridinfo[thickSnowindex].nrl, gridinfo[thickSnowindex].nrh, 1, ncols, nrows);    thickSnow1 = NULL;
			free_matrix(densSnow1, gridinfo[densSnowindex].nrl, gridinfo[densSnowindex].nrh, 1, ncols, nrows);    densSnow1 = NULL;
			free_matrix(upperIndexSnow, gridinfo[upperIndexSnowindex].nrl, gridinfo[upperIndexSnowindex].nrh, 1, ncols, nrows);    upperIndexSnow = NULL;



			if (swes_stackname) free(swes_stackname); swes_stackname = NULL;
			if (swel_stackname) free(swel_stackname); swel_stackname = NULL;
			if (temperSnow_stackname) free(temperSnow_stackname); temperSnow_stackname = NULL;
			if (thetaSnow_stackname) free(thetaSnow_stackname); thetaSnow_stackname = NULL;
			if (hhydSnow_stackname) free(hhydSnow_stackname); hhydSnow_stackname = NULL;
			if (hgeoSnow_stackname) free(hgeoSnow_stackname); hgeoSnow_stackname = NULL;
			if (quSnow_stackname) free(quSnow_stackname); quSnow_stackname = NULL;
			if (thickSnow_stackname) free(thickSnow_stackname); thickSnow_stackname = NULL;
			if (densSnow_stackname) free(densSnow_stackname); densSnow_stackname = NULL;
			if (upperIndexSnow_name) free(upperIndexSnow_name); upperIndexSnow_name = NULL;
			
			delete_grid(thetaSnowindex);
			delete_grid(hhydSnowindex);
			delete_grid(hgeoSnowindex);
			delete_grid(quSnowindex);
			delete_grid(swesindex);
			delete_grid(swelindex);
			delete_grid(thickSnowindex);
			delete_grid(densSnowindex);
			delete_grid(upperIndexSnowindex);

			if (do_TMod)
			{
				free_drei_d_matrix(temperSnow, gridinfo[temperSnowindex].nrl, gridinfo[temperSnowindex].nrh, 1, ncols, levelsSnow, nrows);      temperSnow = NULL;
				free_matrix(temperSnow1, gridinfo[temperSnowindex].nrl, gridinfo[temperSnowindex].nrh, 1, ncols, nrows);      temperSnow1 = NULL;
				delete_grid(temperSnowindex);
			}

			if (LayeredSnowProperties)
			{
				delete[] LayeredSnowProperties;
				LayeredSnowProperties = NULL;
			}
		}

		if (do_ArcVeg)
		{
			free_drei_d_matrix(ExcessIceFactor, gridinfo[ExcessIceFactorindex].nrl, gridinfo[ExcessIceFactorindex].nrh, 1, ncols, levels, nrows);	ExcessIceFactor		= NULL;
			free_matrix(ExcessIceFactor1, gridinfo[ExcessIceFactorindex].nrl, gridinfo[ExcessIceFactorindex].nrh, 1, ncols, nrows);					ExcessIceFactor1	= NULL;

			free_drei_d_matrix(AVINHLP, gridinfo[AVINHLPindex].nrl, gridinfo[AVINHLPindex].nrh, 1, ncols, 8, nrows); 			AVINHLP			= NULL;
			free_matrix(AVINHLP1, gridinfo[AVINHLPindex].nrl, gridinfo[AVINHLPindex].nrh, 1, ncols, nrows);						AVINHLP1		= NULL;

			int nuse = MultiLayerLUTable[1].getNumFractions();
			free_drei_d_matrix(AVINSCAL, gridinfo[AVINSCALindex].nrl, gridinfo[AVINSCALindex].nrh, 1, ncols, nuse, nrows); 		AVINSCAL		= NULL;
			free_matrix(AVINSCAL1, gridinfo[AVINSCALindex].nrl, gridinfo[AVINSCALindex].nrh, 1, ncols, nrows);					AVINSCAL1		= NULL;

			int lv = MultiLayerLUTable[1].getNumLayers();
			free_drei_d_matrix(AVNMIN, gridinfo[AVNMINindex].nrl, gridinfo[AVNMINindex].nrh, 1, ncols, lv, nrows); 			    AVNMIN = NULL;
			free_matrix(AVNMIN1, gridinfo[AVNMINindex].nrl, gridinfo[AVNMINindex].nrh, 1, ncols, nrows);						AVNMIN1 = NULL;

			free_matrix(AVGDD,			gridinfo[AVGDDindex].nrl,			gridinfo[AVGDDindex].nrh,		1, ncols, nrows);   AVGDD			= NULL;
			free_matrix(AVINMAX7DAVG,	gridinfo[AVINMAX7DAVGindex].nrl,	gridinfo[AVINMAX7DAVGindex].nrh,1, ncols, nrows);   AVINMAX7DAVG	= NULL;
			free_matrix(AVTDAVG,		gridinfo[AVTDAVGindex].nrl,			gridinfo[AVTDAVGindex].nrh,		1, ncols, nrows);   AVTDAVG			= NULL;
			free_matrix(AVINDAVG,		gridinfo[AVINDAVGindex].nrl,		gridinfo[AVINDAVGindex].nrh,	1, ncols, nrows);   AVINDAVG		= NULL;
			free_matrix(AVSNOWFREE,		gridinfo[AVSNOWFREEindex].nrl,		gridinfo[AVSNOWFREEindex].nrh,	1, ncols, nrows);   AVSNOWFREE		= NULL;
			free_matrix(AVPGP,          gridinfo[AVPGPindex].nrl,           gridinfo[AVPGPindex].nrh, 1, ncols, nrows);   AVPGP = NULL;
			free_matrix(AVNMINTOT,      gridinfo[AVNMINTOTindex].nrl,       gridinfo[AVNMINTOTindex].nrh, 1, ncols, nrows);   AVNMINTOT = NULL;
			free_matrix(AVSON,          gridinfo[AVSONindex].nrl,           gridinfo[AVSONindex].nrh, 1, ncols, nrows);   AVSON = NULL;

			if (ExcessIceFactor_stackname)	free(ExcessIceFactor_stackname);	ExcessIceFactor_stackname	= NULL;
			if (AVINHLP_stackname)			free(AVINHLP_stackname);			AVINHLP_stackname			= NULL;
			if (AVGDD_name)					free(AVGDD_name);					AVGDD_name					= NULL;
			if (AVINMAX7DAVG_name)			free(AVINMAX7DAVG_name);			AVINMAX7DAVG_name			= NULL;
			if (AVTDAVG_name)				free(AVTDAVG_name);					AVTDAVG_name				= NULL;
			if (AVINDAVG_name)				free(AVINDAVG_name);				AVINDAVG_name				= NULL;
			if (AVSNOWFREE_name)			free(AVSNOWFREE_name);				AVSNOWFREE_name				= NULL;
			if (AVPGP_name)					free(AVPGP_name);					AVPGP_name					= NULL;
			if (AVINSCAL_stackname)			free(AVINSCAL_stackname);			AVINSCAL_stackname			= NULL;
			if (AVNMINTOT_name)				free(AVNMINTOT_name);				AVNMINTOT_name              = NULL;
			if (AVNMIN_stackname)			free(AVNMIN_stackname);				AVNMIN_stackname			= NULL;
			if (AVSON_name)					free(AVSON_name);					AVSON_name					= NULL;

			delete_grid(ExcessIceFactorindex);
			delete_grid(AVINHLPindex);
			delete_grid(AVGDDindex);
			delete_grid(AVINMAX7DAVGindex);
			delete_grid(AVTDAVGindex);
			delete_grid(AVINDAVGindex);
			delete_grid(AVSNOWFREEindex);
			delete_grid(AVPGPindex);
			delete_grid(AVINSCALindex);
			delete_grid(AVNMINindex);
			delete_grid(AVNMINTOTindex);
			delete_grid(AVSONindex);

			if (ARCVEGFRACTIONS) delete[] ARCVEGFRACTIONS; ARCVEGFRACTIONS = NULL;
			if (ARCVEGBIOMASS) delete[] ARCVEGBIOMASS; ARCVEGBIOMASS = NULL;

			if (AVNMinParam)
			{
				for (int i = 0; i <= AVNumClimZon; i++)
				{
					delete[](AVNMinParam[i]);
					AVNMinParam[i] = NULL;
				}
				delete[](AVNMinParam);
				AVNMinParam = NULL;
			}

			double** tempgrid;
			int gather = 1;
			if ((gather > 0) && (mpi_size > 1)) tempgrid = MPIGatherIntoNewGrid(2, gridinfo[2].pointer); else tempgrid = gridinfo[2].pointer;
			if (mpi_rank == 0)
			{
				std::string outname = std::string(gridinfo[2].name) + std::string("_out");
				char *coutname = new char[1000];
				sprintf(coutname, "%s", outname.c_str());
				printgridfile(coutname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids); // write back digital elevation model
				delete[] coutname;
			}
			if ((gather > 0) && (mpi_size > 1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
		}


	}
	modelobject::ClearMemory();
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info,TLastResults *LastResults,int maxLastResults)
{
	/* ---------------------------------------------------------------------------------------------------------*/
	/*                        initialization of unsaturated zone tracer module is different                     */
	/*                        from the other modules because of the stacks that are needed here                 */
	/* ---------------------------------------------------------------------------------------------------------*/
	try
	{
		retvalue=0;
		int l,lv_temp;
		unsigned long bufrow;
		char *identify,*tempstr1,*tempstr2,*tempstr3,*tempstr4;
		if (do_it==0) return;
		identify=newstring(300,1);
		tracer_count = p_tracer_count;

		/* calling the calculation of flow time zones */
		fzt_grid_auswerten();
		//    qd=(double*)calloc(tracer_count+1,sizeof(double));
		//    conc_ifl=(double*)calloc(tracer_count+1,sizeof(double));
		/*
conc_lat=(double*)calloc(tracer_count+1,sizeof(double));
conc_gwin=(double*)calloc(tracer_count+1,sizeof(double));
conc_gwex=(double*)calloc(tracer_count+1,sizeof(double));
conc_sl_old=(double**)calloc(tracer_count+1,sizeof(double*));
for (z=0;z<=tracer_count;z++)
conc_sl_old[z]=(double*)calloc(levels+1,sizeof(double));
*/
		tempstr1=newstring(200,sizeof(char));
		tempstr2=newstring(200,sizeof(char));
		tempstr3=newstring(200,sizeof(char));
		tempstr4=newstring(200,sizeof(char));
		if (tracer_count>0)
		{
			concbuffer    = (buffertype**)calloc(tracer_count+1,sizeof(buffertype*));
			concbufmax    = (unsigned long*)calloc(tracer_count+1,sizeof(unsigned long));
			concbufcnt    = (unsigned long*)calloc(tracer_count+1,sizeof(unsigned long));
			concdatei     = (FILE**)calloc(tracer_count+1,sizeof(FILE*));
			concdateiname = (char**)calloc(tracer_count+1,sizeof(char*));
			conc          = (double****)calloc(tracer_count+1,sizeof(double***));
			conc1         = (double***)calloc(tracer_count+1,sizeof(double**));
			//      conc_sl       = (double**)calloc(tracer_count+1,sizeof(double*));
			
			CONC_QD       = (double***)calloc(tracer_count+1,sizeof(double**));
			CONC_QI       = (double***)calloc(tracer_count+1,sizeof(double**));

			conc_qdindex  = new int[tracer_count+1];
			conc_qiindex  = new int[tracer_count+1];

			tracer_info = p_tracer_info;
			tracer_module = (tracer_module_type*)calloc(tracer_count+1,sizeof(tracer_module_type));
			for (int z=1;z<=tracer_count;z++)
			{

				if (tracer_info[z].radioactive==1)
				/* this is the factor, with which the actual concentration of a radioactive
		tracer has to be multiplied to get the radioactive decay for one time step */
				tracer_info[z].decfact = exp(-(InV/substepcnt/60/24)/(tracer_info[z].half_time*_1_ln2_));

				/* initializing a statistc register set for this tracer */
				tracer_module[z].zone = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
				if (!tracer_module[z].zone) speichererror(1);
				for (int j=0;j<=zonmax;j++)
				{ 
					tracer_module[z].zone[j].n         = tracer_module[z].zone[j].RelCellSum = tracer_module[z].zone[j].sx        = 0;
					tracer_module[z].zone[j].zonenwert = zone[j].zonenwert;
					tracer_module[z].zone[j].rel       = zone[j].rel;
				}

				/* initializing the output stream for the tracer statistics */
				tracer_module[z].statfilename=newstring(_stringsize_+1,1);
				sprintf(tracer_module[z].statfilename,"%s%.3s%.5s.%s",tracer_info[z].path,tracer_info[z].name,instanzname,tracer_info[z].statext);
				tracer_module[z].stat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&tracer_module[z].zone,(int)zonmax,(int)tracer_info[z].statcode,(int)99,(char*)tracer_module[z].statfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));         if (retvalue!=0) throw retvalue;
				sprintf(identify,"tracer concentration for tracer %d (%s) ",z,tracer_info[z].name);
				tracer_module[z].stat->schreibzonenwerte(identify,instanzname,LastResults,maxLastResults);
				/* searching input concentrations:
	if the tracer is evaporating, the input concentration is always the same as in the last storage, 
	but if it is not evaporating (salts), the input-concentration may not rise above the maximum
	concentration, therefore another grid containing the concentration in the flux is allocated and 
	consequently searched for here in such cases */
				sprintf(tempstr4,"concflux_tracer_%1d_%s",z,"input"); 
				sprintf(tempstr3,"concflux_tracer_%1d_%s",z,_irrigation_parameters_); 
				sprintf(tempstr2,"concflux_tracer_%1d_%s",z,_Schneemodell_);
				sprintf(tempstr1,"concflux_tracer_%1d_%s",z,_Interzeptionsmodell_);
				tracer_module[z].cinput=gridinfo[getindex_first_of_4(tempstr1,tempstr2,tempstr3,tempstr4)].pointer;
				// Note: the same tracer input is used in silting up model. The fraction of direct runoff, generated in 
				// the silting up model will have the same tracer concentration as the input here (only fluxes without storage
				// effects are regarded there), so we can safely use the old algorithms here. Later, when mixing tracers etc. the 
				// amount of direct runoff is added up to the runoff generated in the silting up model and the concentration 
				// of any tracer is taken from the normal input, since the already calculated from silting up runoff is a fraction 
				// of that input, too.

				/* initializing buffers for point output files */
				if ((control_row>0) && (MPIREG[control_row]) && (MPIREG[control_row][control_col]==mpi_rank)) /* writing the files only, if control_row is greater than zero */
				{
					/* file buffer for salt concentrations */
					concbufmax[z]      = (unsigned long)((long)buffersize/(4*sizeof(long)+(levels+1)*sizeof(double)));
					if (concbufmax[z]<1) concbufmax[z]=1;
					concbufcnt[z]      = 0; /* initialising the buffer index counter with 0 */
					concbuffer[z]      = (buffertype*)calloc(concbufmax[z]+1,sizeof(buffertype));
					for (bufrow = 1; bufrow <= concbufmax[z]; bufrow++)
					{
						concbuffer[z][bufrow].wert = (double*)calloc(levels + 1, sizeof(double));
					}
					/* allocating memory for each level and the mean + the zero-index-element */
					concdateiname[z]=newstring(_stringsize_+1,1);
					sprintf(concdateiname[z],"%s%.3s%.5s.%s",tracer_info[z].path,tracer_info[z].name,"point",tracer_info[z].statext);
					if ((concdatei[z]=fopen(concdateiname[z],"wt"))==NULL)
					{
						(*WSstream) << "FATAL: error opening tracer concentration file "<<concdateiname[z]<<"[end-entry]";
						throw(65);
					}
					fprintf(concdatei[z],"YY\tMM\tDD\tHH\t%s[%d][%d]",identify,control_row,control_col);
					for (int i = 1; i <= levels; i++)
					{
						fprintf(concdatei[z], "[%d]\t", i);
					}
					fprintf(concdatei[z],"\n");
					fclose(concdatei[z]);
				}

				/* allocating memory for the concentration in a storage */
				int gridindex=1;
				while (gridinfo[gridindex].identification != NULL)
				{ 
					gridindex++;
					if (gridindex == gridmaxanzahl)
					{ 
						(*WSstream) << "FATAL: no free grid left. please contact the author to set the number of grids to a larger value "<<gridmaxanzahl<<"[end-entry]";
						throw(12);
					}
				}
				sprintf(identify,"concstor_tracer_%d_%s",z,instanzname);
				tracer_module[z].gridname=newstring(_stringsize_+1,1);
				sprintf(tracer_module[z].gridname,"%s%.3s%.5s.%s",tracer_info[z].path,tracer_info[z].name,instanzname,tracer_info[z].gridext);
				if (gridslesen!=1)
				{ /*  generate internally */
					conc[z] = drei_d_matrix(mpistartrow,mpiendrow,1,ncols,levels,nrows);
					conc1[z] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
					initialize_grid(gridindex, tracer_module[z].gridname, conc1[z], conc[z], levels, identify, tracer_info[z].writegrid, mpistartrow, mpiendrow);
					gridinfo[gridindex].getfirstlevel = 1;
					for (int k = mpistartrow; k <= mpiendrow; k++)
					{
						for (int l=FirstLastCol[k].First;l<=FirstLastCol[k].Last;l++)
						{
							if ((zn[k][l] != nodata) && (MPIREG[k][l]==mpi_rank))
							{
								for (int lv=1;lv<=(int)GWLEVEL[k][l]-1;lv++)
								{ 
									conc[z][k][l][lv]=tracer_info[z].evaporating? 
									tracer_info[z].soilini : 
										// ToDo: replace ST[(int)SOIL[k][l]].getDzui() by dicken with excessIce consideration
									theta[k][l][lv]*1000.0*ST[(int)SOIL[k][l]].getDzui()[lv]*tracer_info[z].soilini; /* all concentrations(masses are set to soilini */
								}
								for (int lv=(int)GWLEVEL[k][l];lv<=levels;lv++) 
								{
									conc[z][k][l][lv]=tracer_info[z].evaporating? 
									tracer_info[z].gwini : 
										// ToDo: replace ST[(int)SOIL[k][l]].getDzui() by dicken with excessIce consideration
									theta[k][l][lv]*1000.0*ST[(int)SOIL[k][l]].getDzui()[lv]*tracer_info[z].gwini; /* all concentrations(masses are set to soilini */
								}
							}
							else
							{
								for (int lv=1;lv<=levels;lv++) 
								conc[z][k][l][lv]=nodata;
							}
						}
					}
					tracer_module[z].storindex = gridindex;
				}
				else
				{ /* read from disk */
					lv_temp=readstackfromfile(tracer_module[z].gridname,&(conc[z]),_MPIFullGridOff_);
					if (lv_temp!=levels)
					{
						(*WSstream) << "FATAL: wrong number of levels in stack "<<tracer_module[z].gridname
						<<" (there are "<<lv_temp<<" levels, it must be "<<levels<<"[end-entry]";
						throw(13);
					}
					conc1[z] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
					initialize_grid(gridindex, tracer_module[z].gridname, conc1[z], conc[z], levels, identify, tracer_info[z].writegrid, mpistartrow, mpiendrow);
					gridinfo[gridindex].getfirstlevel = 1;
					tracer_module[z].storindex = gridindex;
				}


				// if the silting up model is used, the tracer concentrations are available in another grid
				// we will get the grids address here and put it into tracer_module[z].concflux, a here unused pointer,
				// to have access when mixing the direct runoff fractions from siltingup-model with unsatzonmodel.
				tracer_module[z].concflux=NULL;
				char tracername[200];
				sprintf(tracername,"concflux_tracer_%d_SiltingUpModel",z);
				if (ist_vorhanden(tracername))
				{
					tracer_module[z].concflux=gridinfo[getindex(tracername)].pointer;
				}
				else
				tracer_module[z].concflux=NULL;

				char* lpszInternalName;
				lpszInternalName=new char[1024];

				gridindex++; // Index for the next grid in gridinfo-array
				CONC_QD[z] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				sprintf(lpszInternalName,"_CONC_QD_%i_",z);
				initialize_grid(gridindex, lpszInternalName, CONC_QD[z], NULL, 0, lpszInternalName, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int l=FirstLastCol[k].First;l<=FirstLastCol[k].Last;l++)  
					{
						if ((zn[k][l] != nodata) && (MPIREG[k][l] == mpi_rank))
						{
							CONC_QD[z][k][l] = 0;
						}
						else
						{
							CONC_QD[z][k][l] = nodata;
						}
					}
				}
				conc_qdindex[z] = gridindex;

				gridindex++; // Index for the next grid in gridinfo-array
				CONC_QI[z] = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
				sprintf(lpszInternalName,"_CONC_QI_%i_",z);
				initialize_grid(gridindex, lpszInternalName, CONC_QI[z], NULL, 0, lpszInternalName, 0, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int l = FirstLastCol[k].First; l <= FirstLastCol[k].Last; l++)
					{
						if ((zn[k][l] != nodata) && (MPIREG[k][l] == mpi_rank))
						{
							CONC_QD[z][k][l] = 0;
						}
						else
						{
							CONC_QD[z][k][l] = nodata;
						}
					}
				}
				conc_qiindex[z] = gridindex;

				delete[] lpszInternalName;
			}
		}

		/* ---------------------------------------------------------------------------------------------------------*/
		/*                        end of initialization of unsaturated zone tracer module                           */
		/*            the following sequences are done in every case (even if no tracers are modelled)              */
		/* ---------------------------------------------------------------------------------------------------------*/

		if (gwflow_active)
		{
			if (surface_interaction)
			{
				RIVDEP=gridinfo[getindex(_river_depth_)].pointer;
				RIVWIT=gridinfo[getindex(_river_width_)].pointer;
				RIVLNK=gridinfo[getindex(_river_links_)].pointer;
				KOLM  =gridinfo[getindex(_kolmation_)].pointer;
			}
			else
			{
				RIVDEP=NULL;
				RIVWIT=NULL;
				RIVLNK=NULL;
				KOLM=NULL;
			}

			zoneGWIN = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneGWIN) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{ 
				zoneGWIN[i].n         = zoneGWIN[i].RelCellSum = zoneGWIN[i].sx        = 0;
				zoneGWIN[i].zonenwert = zone[i].zonenwert;
				zoneGWIN[i].rel       = zone[i].rel;
			}

			zoneGWEX = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneGWEX) speichererror(1);
			for (int i=0;i<=zonmax;i++)
			{ 
				zoneGWEX[i].n         = zoneGWEX[i].RelCellSum =  zoneGWEX[i].sx        = 0;
				zoneGWEX[i].zonenwert = zone[i].zonenwert;
				zoneGWEX[i].rel       = zone[i].rel;
			}

			GWINstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneGWIN,(int)zonmax,(int)GWINstatcode,(int)7,(char*)GWINstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));       if (retvalue!=0) throw retvalue;
			GWINstat->schreibzonenwerte("re-infiltration from rivers [mm per Zone]",instanzname,LastResults,maxLastResults);
			GWEXstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&zoneGWEX,(int)zonmax,(int)GWEXstatcode,(int)7,(char*)GWEXstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));       if (retvalue!=0) throw retvalue;
			GWEXstat->schreibzonenwerte("base flow (exfiltration into rivers, only with gw-module active) [mm per Zone]",instanzname,LastResults,maxLastResults);
			/* first: looking for the first empty entry in the gridinfo-array */
			int gridindex = 1;
			while (gridinfo[gridindex].identification != NULL)
			{
				gridindex++;
				if (gridindex == gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no more free grids in memeory, (max. "<<gridmaxanzahl<<")!"<<"[end-entry]";
					throw(12);
				}
			}

			GWIN = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
			initialize_grid(gridindex, GWINgridname, GWIN, NULL, 0, _river_infiltration_, GWINwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						GWIN[k][j] = 0;
					}
					else
					{
						GWIN[k][j] = nodata;
					}
				}
			}
			gwinindex = gridindex;
			gridindex++;

			GWEX = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(gridindex, GWEXgridname, GWEX, NULL, 0, _river_exfiltration_, GWEXwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						GWEX[k][j] = 0;
					}
					else
					{
						GWEX[k][j] = nodata;
					}
				}
			}
			gwexindex = gridindex;
			gridindex++;

			/* allocating memeory for a pointer array for the concentrations within groundwater */
			gw_conc=(double***)calloc(tracer_count+1,sizeof(double**));
			for (int z=1;z<=tracer_count;z++)
			{
				sprintf(identify,"concstor_tracer_%d_layer_1_%s",z,_gwflow_);
				gw_conc[z]=gridinfo[getindex(identify)].pointer;
			}
			gw_thick=gridinfo[getindex(_aquifer_thickness_"1")].pointer;
			riverflow_table=(riverflow_table_type*)calloc(zonmax+1,sizeof(riverflow_table_type));

			/* initializing the storage coefficient for the first aquifer with the drainable porosity of the actual soil 
			-> rapid free drainage should take place up to a suction of 3.45 m*/
			S01=gridinfo[getindex(_gw_storage_coeff_"1")].pointer;

		}
		else
		{
			riverflow_table=NULL;
			gw_conc=NULL;
			gw_thick=NULL;
			GWIN=NULL;
			GWEX=NULL;
			S01=NULL;
		}
		if (GWINgridname    ) free(GWINgridname    );    GWINgridname = NULL;
		if (GWEXgridname    ) free(GWEXgridname    );    GWEXgridname = NULL;

// first old value of lateral groundwater inflow is set to the same value like
// the real first lateral inflow -> no storage change 
//		for (int i = 1; i <= nrows; i++)
//		{
//			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
//			{
//				if (zn[i][j] > nodata)
//				{
//					oldDELTAGW[i][j]=DELTAGW[i][j]; 
//				}
//			}
//		}

			
		Qg = (double**)calloc(tracer_count + 1, sizeof(double*));
			
		for (int z = 0; z <= tracer_count; z++)
		{
			Qg[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}

		QB=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QB[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}
		QBb=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QBb[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}

		QD=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QD[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}
		QDd=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QDd[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}

		QI=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QI[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}
		QIi=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z = 0; z <= tracer_count; z++)
		{
			QIi[z] = (double*)calloc(zonmax + 1, sizeof(double));
		}

		/* estimating first runoff and tracer concentrations for routing initialization (only baseflow is assumed) */
		for (int i = 1; i <= zonmax; i++)
		{
			/* setting zonal statistic registers to zero */
			zone[i].sx=zone[i].n=0.0;
			for (int z = 1; z <= tracer_count; z++)
			{
				tracer_module[z].zone[i].n=0.0;
				tracer_module[z].zone[i].RelCellSum=0.0;
				tracer_module[z].zone[i].sx=0.0;
			}
		}
		/* initializing conc1 - grid (tracer concentrations for baseflow/irrigation) */
		for (int z = 1; z <= tracer_count; z++)
		{
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					int k=(int)GWLEVEL[i][j]; /* concentration of the first partially saturated layer */
					/* level of interest: first level (for level of groundwater: ...= conc[z][i][j][(int)GWLEVEL[i][j]]) */
					if (zn[i][j]!=nodata) 
					{
						if (tracer_info[z].evaporating)
						{
							conc1[z][i][j]=conc[z][i][j][k];
						}
						else
						{
							// ToDo: replace ST[(int)SOIL[k][l]].getDzui() by dicken with excessIce consideration
							conc1[z][i][j]=conc[z][i][j][k]/(theta[i][j][k]*1000.0*ST[(int)SOIL[i][j]].getDzui()[k]);
							if ((conc1[z][i][j]>tracer_info[z].max_conc)&&(tracer_info[z].max_conc>0))
							conc1[z][i][j]=tracer_info[z].max_conc;
						}
						/* groundwater-concentration is set to the value from the concentration matrix conc1 (first saturated layer in conc[tracer][row][col][layer] */
						if (gwflow_active) gw_conc[z][i][j]=conc1[z][i][j];
					}
				}  
			}
		}
		int TG;
		double DT = InV/substepcnt*60.0;           /* interval in secounds (InV: interval in minutes) */
		double timestep = 6.0e4*InV/substepcnt;    /* factor for changing the measure of fluxes from m/s to mm/timestep */
		for (int i=mpistartrow;i<=mpiendrow;i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j]==mpi_rank))
				{
					TG = zarray[(int)zn[i][j]];
					/* estimating baseflow for each grid cell if this is not done in the groundwater flow module */
					if (!gwflow_active)
					{
						//15.01.2008 oldDELTAGW[i][j]=DELTAGW[i][j]=Q0[TG]*ST[(int)SOIL[i][j]].getK_sat()[1]*exp(GWSTAND[i][j]*kelsqlat[TG]);
//						oldDELTAGW[i][j]=
						DELTAGW[i][j]=Q0[TG]*exp((GWSTAND[i][j]*kelsqlat[TG]))/timestep;
						DELTAGW[i][j]*=timestep;
					}
					/* updating statistic registers */
					if (RelCellArea!=NULL)
					{
						zone[TG].sx+=RelCellArea[i][j]*DELTAGW[i][j];
						zone[TG].n++;
						zone[TG].RelCellSum+=RelCellArea[i][j];
					}
					else
					{
						zone[TG].sx+=DELTAGW[i][j];
						zone[TG].n++;
						zone[TG].RelCellSum+=1.0;
					}
					for (int z=1;z<=tracer_count;z++)
					{
						if (RelCellArea!=NULL)
						{
							/* initializing tracer concentrations for baseflow with Groundwater tracer concentrations */
							tracer_module[z].zone[TG].sx+=RelCellArea[i][j]*DELTAGW[i][j]*conc1[z][i][j];
							/* weighted averaging weighed by absolute flow values */
							tracer_module[z].zone[TG].n+=DELTAGW[i][j];
							tracer_module[z].zone[TG].RelCellSum+=RelCellArea[i][j]*DELTAGW[i][j];
						}
						else
						{
							/* initializing tracer concentrations for baseflow with Groundwater tracer concentrations */
							tracer_module[z].zone[TG].sx+=DELTAGW[i][j]*conc1[z][i][j];
							/* weighted averaging weighed by absolute flow values */
							tracer_module[z].zone[TG].n+=DELTAGW[i][j];
							tracer_module[z].zone[TG].RelCellSum+=DELTAGW[i][j];
						}
					}
				}
			}
		}
		for (int i=1;i<=zonmax;i++)
		{
			/* mean values for first baseflow are stored in first total flow for routing initialization */
			Qg[0][i]=zone[i].sx/zone[i].RelCellSum;
			/* restoring zero into the stat. registers */
			zone[i].sx=zone[i].n=zone[i].RelCellSum=0.0;
			for (int z=1;z<=tracer_count;z++)
			{
				Qg[z][i]=tracer_module[z].zone[i].RelCellSum>0?tracer_module[z].zone[i].sx/tracer_module[z].zone[i].RelCellSum:0.0;
				tracer_module[z].zone[i].n=0.0;
				tracer_module[z].zone[i].RelCellSum=0.0;
				tracer_module[z].zone[i].sx=0.0;
			}
		}
		/* reading storage content file in order to initialize the flow travel time zones for interflow and surface flow */
		if (gridslesen==1)
		{
			ReadSlrContents();
		}
		if (tempstr1) free(tempstr1);
		if (tempstr2) free(tempstr2);
		if (tempstr3) free(tempstr3);
		if (tempstr3) free(tempstr4);
		if (identify) free(identify);
		return;
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
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::substance_transport_delete()
{
	int z;
	unsigned int j;
	if (do_it==0) return;
	if (tracer_count>0)
	{
		for (z=1;z<=tracer_count;z++)
		{
			/* free memory for grids */
			free_drei_d_matrix(conc[z], mpistartrow, mpiendrow, 1, ncols, levels,nrows);
			conc[z] = NULL;
			free_matrix(conc1[z], mpistartrow, mpiendrow, 1, ncols, nrows);
			tracer_module[z].concstor = NULL;
			conc1[z] = NULL;
			delete_grid(tracer_module[z].storindex);
			if (tracer_module[z].gridname) free(tracer_module[z].gridname);
			/* free memory for zonal statistics */
			tracer_module[z].stat->ClearMemory();
			if (tracer_module[z].stat) delete(tracer_module[z].stat);
			if (tracer_module[z].statfilename) free(tracer_module[z].statfilename);
			if (tracer_module[z].zone) free(tracer_module[z].zone);
			/* writing the files only, if control_row is greater than zero */
			if ((control_row>0) && (MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
			{
				if (initialized) 
				flushbuffer( concdatei[z], concdateiname[z], &(concbufcnt[z]), concbuffer[z]);
				for (int j = 1; j <= concbufmax[z]; j++)
				{
					if (concbuffer[z][j].wert)
					{
						free(concbuffer[z][j].wert);
					}
				}
				if (concbuffer[z]) free(concbuffer[z]);
				if (concdateiname[z]) free(concdateiname[z]);
			}

			delete_grid(conc_qdindex[z]);
			free_matrix(CONC_QD[z], mpistartrow, mpiendrow, 1, ncols, nrows);    CONC_QD[z] = NULL;
			delete_grid(conc_qiindex[z]);
			free_matrix(CONC_QI[z], mpistartrow, mpiendrow, 1, ncols, nrows);    CONC_QI[z] = NULL;
		}

		if (concbuffer    ) free(concbuffer    ); concbuffer    =NULL;
		if (concbufmax    ) free(concbufmax    ); concbufmax    =NULL;
		if (concbufcnt    ) free(concbufcnt    ); concbufcnt    =NULL;
		if (concdatei     ) free(concdatei     ); concdatei     =NULL;
		if (concdateiname ) free(concdateiname ); concdateiname =NULL;
		if (conc          ) free(conc          ); conc          =NULL;
		if (conc1         ) free(conc1         ); conc1         =NULL;
		//    if (conc_sl       ) free(conc_sl       ); conc_sl       =NULL;

		if (conc_qdindex) delete[] conc_qdindex; conc_qdindex=NULL;
		if (conc_qiindex) delete[] conc_qiindex; conc_qiindex=NULL;

		if (CONC_QD) free(CONC_QD); CONC_QD=NULL;
		if (CONC_QI) free(CONC_QI); CONC_QI=NULL;

	}
	if (gw_conc) { free(gw_conc); }
	gw_conc=NULL;

	//  if (qd) free(qd);
	/*  if (conc_sl_old)
for (z=0;z<=tracer_count;z++) 
{
if (conc_sl_old[z]) 
	free(conc_sl_old[z]);
conc_sl_old[z]=NULL;
}
if (conc_sl_old) free(conc_sl_old); 
conc_sl_old=NULL;

//  if (conc_ifl) free(conc_ifl); conc_ifl=NULL;
if (conc_lat) free(conc_lat); conc_lat=NULL;
if (conc_gwin) free(conc_gwin); conc_gwin=NULL;
if (conc_gwex) free(conc_gwex); conc_gwex=NULL;
*/  
	if (Qg)
	for (z=0;z<=tracer_count;z++) 
	{
		if (Qg[z]) free(Qg[z]);
		Qg[z]=NULL;
	}
	if (Qg) free(Qg);Qg = NULL;

	if (QD)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QD[z]) free(QD[z]);
		QD[z]=NULL;
	}
	if (QD) free(QD);QD = NULL;

	if (QDd)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QDd[z]) free(QDd[z]);
		QDd[z]=NULL;
	}
	if (QDd) free(QDd);QDd = NULL;

	if (QB)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QB[z]) free(QB[z]);
		QB[z]=NULL;
	}
	if (QB) free(QB);QB = NULL;

	if (QBb)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QBb[z]) free(QBb[z]);
		QBb[z]=NULL;
	}
	if (QBb) free(QBb);QBb = NULL;

	if (QI)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QI[z]) free(QI[z]);
		QI[z]=NULL;
	}
	if (QI) free(QI);QI = NULL;

	if (QIi)
	for (z=0;z<=tracer_count;z++) 
	{
		if (QIi[z]) free(QIi[z]);
		QIi[z]=NULL;
	}
	if (QIi) free(QIi);QIi = NULL;

	if (riverflow_table) free(riverflow_table); riverflow_table=NULL;
	return;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void unsatzonmodellklasse::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
		if ((control_row>0) && (MPIREG[control_row]) && (MPIREG[control_row][control_col] == mpi_rank))
		{
			flushbuffer(   qudatei,   qudateiname,   &qubufcnt,   qubuffer);
			flushbuffer(thetadatei,thetadateiname,&thetabufcnt,thetabuffer);
			flushbuffer( hhyddatei, hhyddateiname, &hhydbufcnt, hhydbuffer);
			flushbuffer(otherdatadatei,otherdatadateiname,&otherdatabufcnt,otherdatabuffer,1,16+2+this->iNumVegetationLayers*6);
			flushbuffer( etrdetaildatei, etrdetaildateiname, &etrdetailbufcnt, etrdetailbuffer,0);
			flushbuffer( interflowdatei, interflowdateiname, &interflowbufcnt, interflowbuffer);
			if (do_TMod>0)
			{
				flushbuffer( temperfile, temperfilename, &temperbufcnt, temperbuffer);
			}
			for (z = 1; z <= tracer_count; z++)
			{
				flushbuffer(concdatei[z], concdateiname[z], &(concbufcnt[z]), concbuffer[z]);
			}
		}

		if (QDstat)   QDstat->flushbuffer();
		if (QIstat)   QIstat->flushbuffer();
		if (QBstat)   QBstat->flushbuffer();
		if (gwflow_active) 
		{
			if (GWINstat) GWINstat->flushbuffer();
			if (GWEXstat) GWEXstat->flushbuffer();
		}
		if (Qgstat)   Qgstat->flushbuffer();
#ifdef _WASIMDEBUG_
		if (iteration_stat) iteration_stat->flushbuffer();
#endif
		if (SB_1_stat)   SB_1_stat->flushbuffer();
		if (SB_2_stat)   SB_2_stat->flushbuffer();
		if (ICEOLD)	if (SB_2_IMMO_stat)   SB_2_IMMO_stat->flushbuffer();
		if (ROOTstat)    ROOTstat->flushbuffer();
		if (GWSTANDstat) GWSTANDstat->flushbuffer();
		if (GWNstat)     GWNstat->flushbuffer();
		if (INFILTstat)  INFILTstat->flushbuffer();
		if (MACROstat)   MACROstat->flushbuffer();
		if (methode>2)
		{
			if (SUBSTEPSstat)   SUBSTEPSstat->flushbuffer();
		}
		if (INFEXstat)   INFEXstat->flushbuffer();

		if (MAXPOND) if (PONDstat) PONDstat->flushbuffer();
		if (HDRAIN)  if (QDRAINstat) QDRAINstat->flushbuffer();
		if (do_TMod>0)
		{
			if (this->iTModMethod>0)
			{
				if (thawdepth_stat) thawdepth_stat->flushbuffer();
			}
		}
		
		for (z=1;z<=tracer_count;z++)
		tracer_module[z].stat->flushbuffer();
	}
}

int unsatzonmodellklasse::GetReadGrid()
{
	return gridslesen;
}

void unsatzonmodellklasse::GetExternalGWTable(int iExternalTimeStep)
{
	if (do_it==0) return;
	if (m_bCoupled==false)
	{
		// when coupling the model the first time, try to find the fillable porosity for a correct balance
		char* searchtext=new char[_stringsize_];
		sprintf(searchtext,"%s%d",_gw_storage_coeff_,1);
		int index=ist_vorhanden(searchtext);
		if (index>0)
		{
			FILLPOR = gridinfo[index].pointer;
		}
		else
		{
			FILLPOR=NULL;
		}
		m_bCoupled=true;
		delete[] searchtext;
	}
	// now we know, that we have to initialize the GW table by the values given in the external GWTable file 
	// the HEAD from WaSiM's GW-Model is initialized from the unsatzon-model's gwstand-grid, so we are on the save side
	// to initialize the GWSTAND in unsatzonmodel (gwmodel will not run without unsatzonmodel, except we deliver all 
	// required grid externally!)
	// GWSTAND = DHM - GWTableExterm	
	// When WaSiM runs with it's own GW-Model, then the grid DELTAGW contains the net change in water content for each zell as velocity
	// so we have simply to calculate this entity here.
	// the only problem is, that we only know water tables without information about the drainable porosities
	// so we have to estimate: fillable porosity is given by the difference between the saturation porosity and the actual theta value
	// and the drainable porosity must be estimated by the difference between theta at saturation and the theta above the groundwater in the last time step
	int i;
	int j;
	double drain_theta;
	drain_theta=0;
	double amount;
	double d;
	double sum_amount;
	int startschicht,endschicht;
	for (int i = mpistartrow; i <= mpiendrow; i++)  /* in each layer for all rows */
	{
		for (int j=1;j<=ncols;j++) /* and within the rows for all valid cells (columns) between the first an the last entry */
		{	
			if ((zn[i][j]>nodata)&&(GWTableExtern[i][j]>nodata))
			{
				if ((GWSTAND[i][j]+hgeo[i][j][0])>GWTableExtern[i][j])
				{ // falling gw
					// 1. Startschicht festlegen
					startschicht = (int)GWLEVEL[i][j];
					// 2. Endschicht festlegen
					endschicht = startschicht;
					while (hgeo[i][j][endschicht]>GWTableExtern[i][j]) 
					{
						endschicht++;
						if (endschicht>levels)
						{
							endschicht=levels;
							GWTableExtern[i][j]=hgeo[i][j][endschicht];
							break;
						}
					}
					// 3. Wenn beide Schichten gleich: normale Berechnung wie früher
					if (startschicht == endschicht)
					{
						// use a suction which should result in similar drainable porosities as may be used in a groundwater model						
						//d = ST[(int)SOIL[i][j]].getDzui()[startschicht]; 
						d=3.45;
						drain_theta = ST[(int)SOIL[i][j]].getThetas(startschicht)-ST[(int)SOIL[i][j]].get_theta(startschicht,d);
						amount = ((GWSTAND[i][j]+hgeo[i][j][0])-GWTableExtern[i][j])*drain_theta;
						// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
						theta[i][j][startschicht]-=amount/ST[(int)SOIL[i][j]].getDzui()[startschicht];
						sum_amount=amount;
					}
					// 4. wenn ungleiche Schichten: Schleife über alle Schichten: jeweils Boden bis auf 3.45m Saugspannung entwässern
					// 5. Wassermenge zusammenzählen -> dabei nach Methode 1 oder 2 vorgehen (alle Schritte iinerhalb der Schleife)
					else
					{
						int k;
						// erste Schicht: bis Unterkante dieser Schicht entleeren
						k=startschicht;
						// d = ST[(int)SOIL[i][j]].getDzui()[k]; 
						d=3.45;
						drain_theta = ST[(int)SOIL[i][j]].getThetas(k)-ST[(int)SOIL[i][j]].get_theta(k,d);
						amount = ((GWSTAND[i][j]+hgeo[i][j][0])-hgeo[i][j][k])*drain_theta;
						// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
						theta[i][j][k]-=amount/ST[(int)SOIL[i][j]].getDzui()[k];
						sum_amount=amount;
						// nun Schleife bis vor endschicht
						for (int k=startschicht+1;k<endschicht;k++)
						{
							//d = ST[(int)SOIL[i][j]].getDzui()[k]; 
							d=3.45;
							drain_theta = ST[(int)SOIL[i][j]].getThetas(k)-ST[(int)SOIL[i][j]].get_theta(k,d);
							// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
							amount = ST[(int)SOIL[i][j]].getDzui()[k];
							amount *= drain_theta;
							// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
							theta[i][j][k]-=amount/ST[(int)SOIL[i][j]].getDzui()[k];
							sum_amount+=amount;
						}
						// und nun dasselbe für die letzte Schicht
						k=endschicht;
						//d = ST[(int)SOIL[i][j]].getDzui()[k]; 
						d=3.45;
						drain_theta = ST[(int)SOIL[i][j]].getThetas(k)-ST[(int)SOIL[i][j]].get_theta(k,d);
						amount = (hgeo[i][j][k-1]-GWTableExtern[i][j]);
						amount *= drain_theta;
						// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
						theta[i][j][k]-=amount/ST[(int)SOIL[i][j]].getDzui()[k];
						sum_amount+=amount;
					}
					EXTBALANCE[i][j]=-sum_amount*1000.0; // in mm, falling GW means that water will be missing for the balance
					GWLEVEL[i][j]=endschicht;
				}
				else
				{ //rising GW
					// 1. Startschicht festlegen
					startschicht = (int)GWLEVEL[i][j];
					// 2. Endschicht festlegen
					endschicht = startschicht;
					if (endschicht>1) 
					{
						while ((hgeo[i][j][endschicht-1]<GWTableExtern[i][j]) && (endschicht>1))
						{
							endschicht--;
						}
					}
					// 3. Wenn beide Schichten gleich: normale Berechnung wie früher
					if (startschicht == endschicht)
					{
						// use a suction which should result in similar drainable porosities as may be used in a groundwater model						
						drain_theta = (ST[(int)SOIL[i][j]].getThetas(startschicht)-theta[i][j][startschicht]);
						amount = (GWTableExtern[i][j]-(GWSTAND[i][j]+hgeo[i][j][0]))*drain_theta;
						// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
						theta[i][j][startschicht]+=amount/ST[(int)SOIL[i][j]].getDzui()[startschicht];
						//						theta[i][j][startschicht]=ST[(int)SOIL[i][j]].getThetas(startschicht);
						sum_amount=amount;
					}
					// 4. wenn ungleiche Schichten: Schleife über alle Schichten: jeweils Boden bis auf 3.45m Saugspannung entwässern
					// 5. Wassermenge zusammenzählen -> dabei nach Methode 1 oder 2 vorgehen (alle Schritte iinerhalb der Schleife)
					else
					{
						int k;
						// erste Schicht: bis Unterkante dieser Schicht entleeren
						k=startschicht;
						drain_theta = (ST[(int)SOIL[i][j]].getThetas(k)-theta[i][j][k]);
						amount = (GWTableExtern[i][j]-(GWSTAND[i][j]+hgeo[i][j][0]))*drain_theta;
						//						theta[i][j][k]+=amount/ST[(int)SOIL[i][j]].getDzui()[k];
						theta[i][j][k]=ST[(int)SOIL[i][j]].getThetas(k);
						sum_amount=amount;
						// nun Schleife bis vor endschicht
						for (int k=startschicht-1;k>endschicht;k--)
						{
							drain_theta = (ST[(int)SOIL[i][j]].getThetas(k)-theta[i][j][k]);
							// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
							amount = ST[(int)SOIL[i][j]].getDzui()[k];
							amount *= drain_theta;
							//							theta[i][j][k]+=amount/ST[(int)SOIL[i][j]].getDzui()[k];
							theta[i][j][k]=ST[(int)SOIL[i][j]].getThetas(k);
							sum_amount+=amount;
						}
						// und nun dasselbe für die letzte Schicht
						k=endschicht;
						drain_theta = (ST[(int)SOIL[i][j]].getThetas(k)-theta[i][j][k]);
						amount = (GWTableExtern[i][j]-hgeo[i][j][k])*drain_theta;
						// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
						theta[i][j][k]+=amount/ST[(int)SOIL[i][j]].getDzui()[k];
						if (theta[i][j][k]>ST[(int)SOIL[i][j]].getThetas(k))
						{
							theta[i][j][k] = ST[(int)SOIL[i][j]].getThetas(k);
							// ToDo: replace ST[(int)SOIL[i][j]].getDzui() by dicken with excessIce consideration
							amount = ST[(int)SOIL[i][j]].getDzui()[k]*drain_theta;
						}
						//						theta[i][j][k]=ST[(int)SOIL[i][j]].getThetas(k);
						sum_amount+=amount;
					}
					EXTBALANCE[i][j]=sum_amount*1000.0; // in mm
					GWLEVEL[i][j]=endschicht;
				}
				//				DELTAGW[i][j]=EXTBALANCE[i][j]/(60.0*iExternalTimeStep);
				//				DELTAGW is zero because GWSTAND is directly taken from the external grid now -> this avoids even small errors due to rounding 
				//				when converting absolute amounts to fluxes and vice versa
				//				The EXTBALANCE-grid holds the balances of additional or missing water.
				DELTAGW[i][j]=0;
				GWSTAND[i][j]=GWTableExtern[i][j]-hgeo[i][j][0];

			}
			else
			{
				if (zn[i][j]>nodata)
				{   // valid cell in WaSiM but not in the external model
					EXTBALANCE[i][j]=0;
					// no change to the other values because WaSiM will use the internal values
					// DELTAGW[i][j]=0;
					// GWSTAND[i][j]=0;
				}
				else
				{   // nodata in WaSiM
					EXTBALANCE[i][j]=nodata;
					DELTAGW[i][j]=nodata;
					GWSTAND[i][j]=nodata;
				}
			}
		}
	}
}

void unsatzonmodellklasse::temperaturePreSpinUp(long me,long se,long tg,long mt,long jr)
{

	if ((!T) || (do_TMod==0))
	{
		(*WSstream) << "FATAL: heat transfer model not active in unsatzon_model. This is required for pre-SpinUp!" << "[end-entry]";
		throw 1040;
	}

	if (!FDD)
	{
		// looking for FDD grid, create if not found
		if (ist_vorhanden(_FreezingDegreeDays_))
		{
			FDD = gridinfo[getindex(_FreezingDegreeDays_)].pointer;
		}
		else
		{ // not found, create it
			fddindex=1;
			while (gridinfo[fddindex].identification!=NULL)
			{ 
				fddindex++;
				if (fddindex==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no more free grids in memory, (max. "<<gridmaxanzahl<<")!"<<"[end-entry]";
					throw(12);
				}
			}
			tddindex=fddindex+1;
			tdayindex=tddindex+1;
			tyearindex=tdayindex+1;

			FDD = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
			initialize_grid(fddindex, _FreezingDegreeDays_, FDD, NULL, 0, _FreezingDegreeDays_, 0, mpistartrow,mpiendrow);
			TDD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(tddindex, _ThawingDegreeDays_, TDD, NULL, 0, _ThawingDegreeDays_, 0, mpistartrow, mpiendrow);
			TDAY = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(tdayindex, _DailyTemperature_, TDD, NULL, 0, _DailyTemperature_, 0, mpistartrow, mpiendrow);
			TYEAR = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(tyearindex, _AnnualTemperature_, TDD, NULL, 0, _AnnualTemperature_, 0, mpistartrow, mpiendrow);

			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j=FirstLastCol[k].First;j<=FirstLastCol[k].Last;j++)  
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						FDD[k][j] = zn[k][j] == nodata ? nodata : 0;
						TDD[k][j] = zn[k][j] == nodata ? nodata : 0;
						TDAY[k][j] = zn[k][j] == nodata ? nodata : 0;
						TYEAR[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
			}

		}
	}

	numIntPreSpinUp++;
	numIntPreSpinUpAll++;

#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{																				 
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				// update for daily average temperature
				TDAY[i][j]+=T[i][j];
				// for annual average temperature 
				TYEAR[i][j]+=T[i][j];
				// if a day is complete, compute updates fr FDD and TDD
				if (60 * (se - 1 + int(hour_offset)) + me + InV == 1440)
//				if ((60*se + int(me+60*hour_offset)) == 1440)
				{
					// compute daily average temperature
					TDAY[i][j]/=numIntPreSpinUp;
					// adding up FDD and TDD here
					if (TDAY[i][j]>=0)
					{
						TDD[i][j]+=TDAY[i][j];
					} 
					else 
					{ 
						FDD[i][j]-=TDAY[i][j]; 
					} 
					// reset daily average temperature
					TDAY[i][j]=0;
					// and set counter to 0
				}
			}
		}
	}
	// if ((60 * se + int(me + 60 * hour_offset)) == 1440)
	if (60 * (se - 1 + int(hour_offset)) + me + InV == 1440)
	{
		numIntPreSpinUp=0;
	}
}

void unsatzonmodellklasse::initializeSoilTemperatures()
{
	double numYears=((double)numIntPreSpinUpAll*InV)/(double(1440.0)*double(365.25));
#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{																				 
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
		{
//			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank) && !(doSpinUp && HRUMaster[i][j]==nodata))
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				// annual average temperature
				TYEAR[i][j]/=numIntPreSpinUpAll;
				// Setting soil temperatures here
				int s=(int)SOIL[i][j];
				int use=(int)LUSE[1][i][j];
				int PFTableLayer=0; // numerical layer of the PermaFrostTable in the Richards an heat transfer model 
				double depth=0;
				double PFTDepth=1.0; // Permafrost Table Depth, can also be computed by a regresion of some type...
				while ((depth<PFTDepth)&&(PFTableLayer<levels))
				{
					PFTableLayer++;
					depth+=ST[s].getDzui()[PFTableLayer];
				}
				// the thermal parameters will be taken from near the surface, i.e. a layer that is much more affected by the 
				// atmospheric temperatur changes and can thaw and freeze etc.
				// if the permafrost table is at 1m depth and in layer 10, the reference layer for the parameters is 10/3 = 3 (at least 1)
				int representativeLayer=PFTableLayer/3;
				if (representativeLayer==0) representativeLayer=1;

				double nFactorFreezing = LUtable[use].GetNFactorFreezing(WSstream);
				double nFactorThawing  = LUtable[use].GetNFactorThawing(WSstream);
				if (nFactorFreezing<-9999)
				{
					nFactorFreezing=dTModNFactorFreezing;
				}
				if (nFactorThawing<-9999)
				{
					nFactorThawing=dTModNFactorThawing;
				}

				double LambdaWater = this->dTModLamWat;
				double LambdaIce = this->dTModLamIce;
				double lambda_air = 0.0262;

				double Kdry = ST[s].getDryThermalConduct(representativeLayer);
				double porosity = ST[s].getThetas()[representativeLayer];
				double acttheta = theta[i][j][representativeLayer];
				double soilfrac = 1-ST[s].getThetas()[representativeLayer]; // volume fraction for solid soil material
				double SoilAndWater = acttheta + soilfrac;

				double Kthawed = (porosity-acttheta)*lambda_air + SoilAndWater * pow(Kdry,soilfrac/SoilAndWater) * pow(LambdaWater,acttheta/SoilAndWater);
				double Kfrozen = (porosity-acttheta)*lambda_air + SoilAndWater * pow(Kdry,soilfrac/SoilAndWater) * pow(LambdaIce,acttheta/SoilAndWater);
		
				double Tpft=0; // temperature of the permafrost table
				// approach after Kudryavtsev et al. (1974)
				Tpft = (Kthawed * nFactorThawing * TDD[i][j] - Kfrozen * nFactorFreezing * FDD[i][j])/(Kfrozen*numYears*365.0);
				// setting all layers below the permafrost table to a temperature according to the geothermal heat flux
				if (Tpft<0)
				{
					// now setting the temperatures for the layer with the permafrost table to Tpft
					temper[i][j][PFTableLayer]=Tpft;
					for (int k=PFTableLayer;k<levels;k++)
					{
						Kdry = ST[s].getDryThermalConduct(k+1);
						porosity = ST[s].getThetas()[k+1];
						acttheta = theta[i][j][k+1];
						soilfrac = 1-ST[s].getThetas()[k+1]; // volume fraction for solid soil material
						SoilAndWater = acttheta + soilfrac;
						double lambda=0;
						if (temper[i][j][k]>0)
						{
							lambda = (porosity-acttheta)*lambda_air + SoilAndWater * pow(Kdry,soilfrac/SoilAndWater) * pow(LambdaWater,acttheta/SoilAndWater);
						}
						else
						{
							lambda = (porosity-acttheta)*lambda_air + SoilAndWater * pow(Kdry,soilfrac/SoilAndWater) * pow(LambdaIce,acttheta/SoilAndWater);
						}

						temper[i][j][k+1] = temper[i][j][k]+ST[s].getDzui()[k+1] *(useHeatFluxAsLowerBC>0?(QLOWBND[i][j]):0.06)/lambda; 
					}
					// setting all layers above the permafrost table to another temperature, e.g. a linear interpolation bteween Tpft and TYEAR
					for (int k=PFTableLayer;k>0;k--)
					{
						temper[i][j][k-1] = interpolate(hgeo[i][j][PFTableLayer],Tpft,hgeo[i][j][0],TYEAR[i][j],hgeo[i][j][k-1]);
					}
				}
				else
				{
					// using TYEAR and QLOWBND or a constant geothermal heatr flux to initialize temperatures
					temper[i][j][1]=TYEAR[i][j];
					for (int k=1;k<levels;k++)
					{
						Kdry = ST[s].getDryThermalConduct(k+1);
						porosity = ST[s].getThetas()[k+1];
						acttheta = theta[i][j][k+1];
						soilfrac = 1-ST[s].getThetas()[k+1]; // volume fraction for solid soil material
						SoilAndWater = acttheta + soilfrac;
						double lambda = (porosity-acttheta)*lambda_air + SoilAndWater * pow(Kdry,soilfrac/SoilAndWater) * pow(LambdaWater,acttheta/SoilAndWater);
						temper[i][j][k+1] = temper[i][j][k]+ST[s].getDzui()[k+1] *(useHeatFluxAsLowerBC>0?(QLOWBND[i][j]):0.06)/lambda; 
					}
				}
			}
		}
	}
	printgridfile((char *)_FreezingDegreeDays_, FDD, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	printgridfile((char *)_ThawingDegreeDays_, TDD, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	printgridfile((char *)_AnnualTemperature_, TYEAR, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	printstackfile((char *)_Soil_Temperature_Stack_, temper, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
	// free allocated memory
	free_matrix(FDD, mpistartrow, mpiendrow, 1, ncols,nrows);       FDD = NULL;
	free_matrix(TDD, mpistartrow, mpiendrow, 1, ncols, nrows);       TDD = NULL;
	free_matrix(TDAY, mpistartrow, mpiendrow, 1, ncols, nrows);      TDAY = NULL;
	free_matrix(TYEAR, mpistartrow, mpiendrow, 1, ncols, nrows);     TYEAR = NULL;
	delete_grid(fddindex);
	delete_grid(tddindex);
	delete_grid(tdayindex);
	delete_grid(tyearindex);
	do_initializeTemper=false; // otherwise temper-stack would be initialized during first time step from air temperature (as is done without pre-spin-up)
}

int unsatzonmodellklasse::readLayeredSnowModelParameters()
{
	do_LayeredSnow = 0;

	// now reading the layered snow model parameters if this section exists and do_LayeredSnow = 1 
	if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
	{
		(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
		throw(28);
	}
	if (findeabschnitt(_LayeredSnowModel_, false))
	{
		fscanfln(steuerdatei, "%d", &do_LayeredSnow, "", tstr);
		if (do_LayeredSnow>0)
		{
			if (this->do_TMod == 0)
			{
				(*WSstream) << "ERROR:\n-----------------\n Error in model configuration: when using [LayeredSnowModel], [heat_transfer] must be used as well.\n"
					<< "It is currently not implemented to rout water through the snow without accounting\n"
					<< "for the temperature of the snow pack. So you need to either switch heat_transfer on or LayeredSnowModel off\n-----------------\n" << "[end-entry]";
				throw(1041);
			}

			int dummyi;
			fscanfln(steuerdatei, "%d", &iLayeredSnowMethod, " method of layered snow model ", tstr);
			fscanfln(steuerdatei, "%lf", &dLayeredSnowMinDensity, " dLayeredSnowMaxDensity", tstr);
			fscanfln(steuerdatei, "%lf", &dLayeredSnowMaxDensity, " dLayeredSnowMaxDensity", tstr);
			fscanfln(steuerdatei, "%lf", &dTOffsetFreshSnow, " dTOffsetFreshSnow", tstr);
			fscanfln(steuerdatei, "%lf", &dDrainageDensitiySnow, " dDrainageDensitiySnow", tstr);
			fscanfln(steuerdatei, "%lf", &dWaterHoldingCapacity, " dWaterHoldingCapacity", tstr);
			fscanfln(steuerdatei, "%lf", &dSEeffectiveness, " dSEeffectiveness", tstr);
			fscanfln(steuerdatei, "%lf", &dSEeffectivenessSoil, " dSEeffectivenessSoil", tstr);

			swes_stackname = newstring(_stringsize_ + 1, 1);
			swel_stackname = newstring(_stringsize_ + 1, 1);
			temperSnow_stackname = newstring(_stringsize_ + 1, 1);
			thetaSnow_stackname = newstring(_stringsize_ + 1, 1);
			hhydSnow_stackname = newstring(_stringsize_ + 1, 1);
			hgeoSnow_stackname = newstring(_stringsize_ + 1, 1);
			quSnow_stackname = newstring(_stringsize_ + 1, 1);
			thickSnow_stackname = newstring(_stringsize_ + 1, 1);
			densSnow_stackname = newstring(_stringsize_ + 1, 1);
			upperIndexSnow_name = newstring(_stringsize_ + 1, 1);

			fscanfln(steuerdatei, "%s", swes_stackname, "name of the SWES-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &swes_writestack, "Writecode swes-stack in ", tstr);

			fscanfln(steuerdatei, "%s", swel_stackname, "name of the SWEL-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &swel_writestack, "Writecode swel-stack in ", tstr);

			fscanfln(steuerdatei, "%s", temperSnow_stackname, "name of the temperSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &temperSnow_writestack, "Writecode temperSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", thetaSnow_stackname, "name of the thetaSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &thetaSnow_writestack, "Writecode thetaSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", hhydSnow_stackname, "name of the hhydSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &hhydSnow_writestack, "Writecode hhydSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", hgeoSnow_stackname, "name of the hgeoSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &hgeoSnow_writestack, "Writecode hgeoSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", quSnow_stackname, "name of the quSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &quSnow_writestack, "Writecode quSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", thickSnow_stackname, "name of the thickSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &thickSnow_writestack, "Writecode thickSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", densSnow_stackname, "name of the densSnow-stack in ", tstr);
			fscanfln(steuerdatei, "%d", &densSnow_writestack, "Writecode densSnow-stack in ", tstr);

			fscanfln(steuerdatei, "%s", upperIndexSnow_name, "name of the upperIndexSnow_grid in ", tstr);
			fscanfln(steuerdatei, "%d", &upperIndexSnow_writegrid, "Writecode upperIndexSnow_writegrid in ", tstr);

		}
	}

	if (fclose(steuerdatei) == EOF)
	{
		(*WSstream) << "FATAL: control file \"" << steuerdateiname << "\" cannot be closed" << "[end-entry]";
		throw(11);
	}

	return do_LayeredSnow;
}

// reads the snow layers properties. This table is of the same type es the soil table, only that it is not included in the soil table in order to allow separate configurations for soil and nsow
// especially since snow and soil may use a different num,ber of layers
void unsatzonmodellklasse::readLayeredSnowProperties()
{
	int i, usenr;

	// create array and read entries from table LayeredSnowProperties

	if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
	{
		(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
		throw(28);
	}
	findeabschnitt(_LayeredSnowProperties_);
	fscanf(steuerdatei, "%d", &usenr);
	tstr = fgets(tstr, _stringsize_, steuerdatei);

	/* Erzeugen und initialisieren der Bodenartentabelle. Die Werte werden in
	defaultini aus der Steuerdatei eingelesen */
	LayeredSnowProperties = new SOILtabletype[usenr + 1];
	for (i = 1; i <= usenr; i++)
	{
		LayeredSnowProperties[i].setCode(i);
	}

	/* im nicht genutzten ersten Element des SOILtable-Arrays wird die gesamtzahl an Bodenarten
	festgehalten um beim spaeteren Zugriff diesen WErt noch zu kennen und um beim Freigeben
	keinen Fehler zu provozieren */
	LayeredSnowProperties[0].setCode(usenr);

	for (i = 1; i <= usenr; i++)
	{
		tstr = fgets(tstr, _stringsize_, steuerdatei);
		std::string record(tstr);
		if ((record.find("MultipleHorizons") != std::string::npos) || (record.find("method") != std::string::npos) || (record.find("{") != std::string::npos))
		{
			while (record.find("}") == std::string::npos)
			{
				tstr = fgets(tstr, _stringsize_, steuerdatei);
				record.append(std::string(tstr));
			}
			LayeredSnowProperties[i].initialize(WSstream, record, "MultipleHorizons");
		}
		else
		{
			(*WSstream) << "FATAL: entries in [LayeredSnowProperties] must use the MultipleHorizons syntax, i.e. no very-old-style soil table entries are allowed here\n\n" << "[end-entry]";
			throw(30);
		}
		if ((i>1) && (LayeredSnowProperties[i].getLevels() != LayeredSnowProperties[i - 1].getLevels()))
		{
			(*WSstream) << "FATAL: the numbers of snow layers has to be equal for all snow types\n"
				"(here: snow type " << LayeredSnowProperties[i].getCode() << ": " << LayeredSnowProperties[i].getLevels()
				<< " layers, but snow type " << LayeredSnowProperties[i - 1].getCode() << ": "
				<< LayeredSnowProperties[i - 1].getLevels() << " layers)\n\n" << "[end-entry]";
			throw(30);
		}
	}

	if (fclose(steuerdatei) == EOF)
	{
		(*WSstream) << "FATAL: control file \"" << steuerdateiname << "\" cannot be closed" << "[end-entry]";
		throw(11);
	}

	for (int i = 1; i <= LayeredSnowProperties[0].getCode(); i++) /* LayeredSnowProperties is of type soil table. In LayeredSnowProperties[0].code the number of "soil types" (one snow type for now) is stored */
	{
		for (int k = 1; k <= LayeredSnowProperties[i].getLevels(); k++)
		{
			if (LayeredSnowProperties[i].getDryDensity(k)<0)
			{
				LayeredSnowProperties[i].setDryDensity(k, dTModRhoDrySoil);
			}

			if (LayeredSnowProperties[i].getDryHeatCapacity(k)<0)
			{
				LayeredSnowProperties[i].setDryHeatCapacity(k, dTModCDrySoil);
			}

			if (LayeredSnowProperties[i].getDryThermalConduct(k)<0)
			{
				LayeredSnowProperties[i].setDryThermalConduct(k, dTModLamDrySoil);
			}
		}

		if (LayeredSnowProperties[i].getKMinFrozenSoil()<0)
		{
			LayeredSnowProperties[i].setKMinFrozenSoil(dTModKmin);
		}

		for (int k = 1; k <= LayeredSnowProperties[i].getLevels(); k++)
		{
			LayeredSnowProperties[i].LookUpTableGenerationTMod(k, this->dTModFactor);
		}
	}
}

// reading in or creating the grids and stacks required for layered snow model
void unsatzonmodellklasse::ReadOrCreateLayeredSnowGridsAndStacks()
{
	// looking for the surface temperature and energy balance grids that will be used as input/boundary condition grids for the heat transfer model
	// if they are not found, the heat transfer model cannot be used for snow.
	if (ist_vorhanden(_surftemp_))				SURFTEMP	= gridinfo[getindex(_surftemp_)].pointer;
	if (ist_vorhanden(_surftemp_onGlac_))		SURFTEMPGL	= gridinfo[getindex(_surftemp_onGlac_)].pointer;
	if (ist_vorhanden(_EnergyBalance_))			EBAL		= gridinfo[getindex(_EnergyBalance_)].pointer;
	if (ist_vorhanden(_EnergyBalance_onGlac_))	EBALGL		= gridinfo[getindex(_EnergyBalance_onGlac_)].pointer;
	if (ist_vorhanden(_Regenanteil_))			P_RAIN	    = gridinfo[getindex(_Regenanteil_)].pointer;
	if (ist_vorhanden(_snowstorage_solid_))		SSNO	    = gridinfo[getindex(_snowstorage_solid_)].pointer;
	if (ist_vorhanden(_snowstorage_liquid_))	SLIQ        = gridinfo[getindex(_snowstorage_liquid_)].pointer;
	if (ist_vorhanden(_Schneespeicher_Total_))	SSTO       = gridinfo[getindex(_Schneespeicher_Total_)].pointer;
	if (ist_vorhanden(_Snow_out_))				QSNOW_OUT  = gridinfo[getindex(_Snow_out_)].pointer;

	int lv_temp=0;

	/* first: looking for the first empty entry in the gridinfo-array */
	int gridindex = 1;
	while (gridinfo[gridindex].identification != NULL)
	{
		gridindex++;
		if (gridindex == gridmaxanzahl)
		{
			(*WSstream) << "FATAL: no more free grids in memeory, (max. " << gridmaxanzahl << ")!" << "[end-entry]";
			throw(12);
		}
	}

	// how many levels does the snow layers have (all table entries should have the same count, which was already checked when reading in the snowLayerProperties
	levelsSnow = 0;
	for (int i = 1; i <= LayeredSnowProperties[0].getCode(); i++)
	{
		if (LayeredSnowProperties[i].getLevels()>levelsSnow) levelsSnow = LayeredSnowProperties[i].getLevels(); 
		LayeredSnowProperties[0].setLevels(levelsSnow);
	}
	
	// now constructing new grids/stacks or rtead them from disk
	if ((gridslesen == 1) || (gridslesen == 3))
	{ // reading stacks from disk
		// the theta-stack
		lv_temp = readstackfromfile(thetaSnow_stackname, &thetaSnow, _MPIFullGridOff_); /* lv_temp is the number of levels from the stack, read in */
		if (lv_temp != levelsSnow) /* comparing the read number of levels with the maximum number of levels from the soil table */
		{
			(*WSstream) << "FATAL: wrong number of snow levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		thetaSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, thetaSnow_stackname, thetaSnow1, thetaSnow, lv_temp, _Theta_Snow_, thetaSnow_writestack, mpistartrow, mpiendrow); /* because theta is not a grid but a stack, the adress-parameter for the grid is NULL */
		thetaSnowindex = gridindex;
		gridindex++;

		/* hgeoSnow: the stack for topographic altitudes --> may change when snow gets more compact */
		lv_temp = readstackfromfile(hgeoSnow_stackname, &hgeoSnow, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		hgeoSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, hgeoSnow_stackname, hgeoSnow1, hgeoSnow, lv_temp, _geodetic_height_Snow_, hgeoSnow_writestack, mpistartrow, mpiendrow);
		hgeoSnowindex = gridindex;
		gridindex++;

		/* hhydSnow: the stack for hydraulic heads */
		lv_temp = readstackfromfile(hhydSnow_stackname, &hhydSnow, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		hhydSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, hhydSnow_stackname, hhydSnow1, hhydSnow, lv_temp, _hydraulic_head_Snow_, hhydSnow_writestack, mpistartrow, mpiendrow);
		hhydSnowindex = gridindex;
		gridindex++;

		// SWES Snow water equivalent solid, i.e. "real" snow, without water in the snow pack
		lv_temp = readstackfromfile(swes_stackname, &swes, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		swes1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, swes_stackname, swes1, swes, lv_temp, _SnowWaterEquivalentSolid_, swes_writestack, mpistartrow, mpiendrow);
		swesindex = gridindex;
		gridindex++;

		// SWEL Snow water equivalent liquid, i.e. the water in the snow pack
		lv_temp = readstackfromfile(swel_stackname, &swel, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		swel1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, swel_stackname, swel1, swel, lv_temp, _SnowWaterEquivalentLiquid_, swel_writestack, mpistartrow, mpiendrow);
		swelindex = gridindex;
		gridindex++;

		// SWEL Snow water equivalent liquid, i.e. the water in the snow pack
		lv_temp = readstackfromfile(thickSnow_stackname, &thickSnow, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		thickSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, thickSnow_stackname, thickSnow1, thickSnow, lv_temp, _SnowLayerThickness_, thickSnow_writestack, mpistartrow, mpiendrow);
		thickSnowindex = gridindex;
		gridindex++;

		// SWEL Snow water equivalent liquid, i.e. the water in the snow pack
		lv_temp = readstackfromfile(densSnow_stackname, &densSnow, _MPIFullGridOff_);
		if (lv_temp != levelsSnow)
		{
			(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
				<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
			throw(13);
		}
		densSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, densSnow_stackname, densSnow1, densSnow, lv_temp, _SnowDensity_, densSnow_writestack, mpistartrow, mpiendrow);
		densSnowindex = gridindex;
		gridindex++;

		readgridfromfile(upperIndexSnow_name, &upperIndexSnow, 1, _upperIndexSnow_, _ShowMessagesOn_, levelsSnow, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
		initialize_grid(gridindex, upperIndexSnow_name, upperIndexSnow, NULL, 0, _upperIndexSnow_, upperIndexSnow_writegrid, mpistartrow, mpiendrow);
		upperIndexSnowindex = gridindex;
		gridindex++;

		if (do_TMod>0)
		{
			/* temper: the stack for soil temperatures */
			lv_temp = readstackfromfile(temperSnow_stackname, &temperSnow, _MPIFullGridOff_);
			if (lv_temp != levelsSnow)
			{
				(*WSstream) << "FATAL: wrong number of levels in stack " << gridinfo[gridindex].name
					<< " (there are " << lv_temp << " level, it must be " << levelsSnow << "[end-entry]";
				throw(13);
			}

			temperSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, temperSnow_stackname, temperSnow1, temperSnow, lv_temp, _Soil_Temperature_Stack_Snow_, temperSnow_writestack, mpistartrow, mpiendrow);
			temperSnowindex = gridindex;
			gridindex++;
		}
		do_initializeTemperSnow = false;
	}
	else /* stacks and grids are not read from disk but generated internally */
	{
		/* 1. stacks */
		/* geographical altitude [m a.s.l.] */
		hgeoSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		hgeoSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, hgeoSnow_stackname, hgeoSnow1, hgeoSnow, levelsSnow, _geodetic_height_Snow_, hgeoSnow_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
//					hgeoSnow[k][j][levelsSnow] = gridinfo[2].pointer[k][j]; /* the zero-level contains then the elevation on the surface from the digital elevation model*/
					hgeoSnow[k][j][levelsSnow] = hgeo[k][j][0]; /* the zero-level contains then the elevation on the surface from the digital elevation model*/
					for (int lv = levelsSnow - 1; lv >= 0; lv--)
					{
						hgeoSnow[k][j][lv] = hgeoSnow[k][j][lv + 1] + LayeredSnowProperties[1].getDzui()[lv];
					}
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						hgeoSnow[k][j][lv] = nodata;
					}
				}
			}
		}
		hgeoSnowindex = gridindex;
		gridindex++;

		/* a stack for water content for all levels (theta-values) */
		thetaSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		thetaSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, thetaSnow_stackname, thetaSnow1, thetaSnow, levelsSnow, _Theta_Snow_, thetaSnow_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv = 1; lv <= levelsSnow; lv++)
					{
						thetaSnow[k][j][lv] = LayeredSnowProperties[1].getThetar()[lv];
					}
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						thetaSnow[k][j][lv] = nodata;
					}
				}
			}
		}
		thetaSnowindex = gridindex;
		gridindex++;

		/* hydraulic head */
		hhydSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		hhydSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, hhydSnow_stackname, hhydSnow1, hhydSnow, levelsSnow, _hydraulic_head_Snow_, hhydSnow_writestack, mpistartrow, mpiendrow);
		double h_temp;
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv = levelsSnow; lv >= 1; lv--)
					{
						h_temp = LayeredSnowProperties[1].h_calc(lv, thetaSnow[k][j][lv]);
						hhydSnow[k][j][lv] = (hgeoSnow[k][j][lv - 1] + hgeoSnow[k][j][lv])*0.5 - h_temp;
						/* hydraulic head is calculated from mean geograohical elevation minus suction */
					}
					hhydSnow[k][j][0] = hhydSnow[k][j][1] + 0.5*LayeredSnowProperties[1].getDzui()[1]; // hydraulic head is always valid for the center of a alyer, so we use here the value of the uppermost layer plus one half of the layer thickness
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						hhydSnow[k][j][lv] = nodata;
					}
				}
			}
		}
		hhydSnowindex = gridindex;
		gridindex++;

		/* swes head */
		swes = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		swes1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, swes_stackname, swes1, swes, levelsSnow, _SnowWaterEquivalentSolid_, swes_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					swes[k][j][0] = 0;
					for (int lv = 1; lv <= levelsSnow; lv++)
					{
						swes[k][j][lv] = 0;
					}
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						swes[k][j][lv] = nodata;
					}
				}
			}
		}
		swesindex = gridindex;
		gridindex++;

		/* swel head */
		swel = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		swel1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, swel_stackname, swel1, swel, levelsSnow, _SnowWaterEquivalentLiquid_, swel_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					swel[k][j][0] = 0;
					for (int lv = 1; lv <= levelsSnow; lv++)
					{
						swel[k][j][lv] = 0;
					}
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						swel[k][j][lv] = nodata;
					}
				}
			}
		}
		swelindex = gridindex;
		gridindex++;

		/* thickSnow */
		thickSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		thickSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, thickSnow_stackname, thickSnow1, thickSnow, levelsSnow, _SnowLayerThickness_, thickSnow_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						thickSnow[k][j][lv] = nodata;
					}
//					thickSnow[k][j][0] = 0;
//					for (int lv = 1; lv <= levelsSnow; lv++)
//					{
//						thickSnow[k][j][lv] = LayeredSnowProperties[1].getDzui()[lv];
//					}
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						thickSnow[k][j][lv] = nodata;
					}
				}
			}
		}
		thickSnowindex = gridindex;
		gridindex++;

		/* thickSnow */
		densSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
		densSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, densSnow_stackname, densSnow1, densSnow, levelsSnow, _SnowDensity_, densSnow_writestack, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						densSnow[k][j][lv] = nodata;
					}
//					for (int lv = 1; lv <= levelsSnow; lv++)
//					{
//						densSnow[k][j][lv] = LayeredSnowProperties[1].getDryDensity()[lv]/1000.0;
//					}
//					densSnow[k][j][0] = densSnow[k][j][0];
				}
				else
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						densSnow[k][j][lv] = nodata;
					}
				}
			}
		}
		densSnowindex = gridindex;
		gridindex++;

		upperIndexSnow = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
		initialize_grid(gridindex, upperIndexSnow_name, upperIndexSnow, NULL, 0, _upperIndexSnow_, upperIndexSnow_writegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					upperIndexSnow[k][j] = levelsSnow;
				}
				else
				{
					upperIndexSnow[k][j] = nodata;
				}
			}
		}
		upperIndexSnowindex = gridindex;
		gridindex++;


		if (do_TMod > 0)
		{
			/* soil temperatures */
			temperSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
			temperSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(gridindex, temperSnow_stackname, temperSnow1, temperSnow, levelsSnow, _Soil_Temperature_Stack_Snow_, temperSnow_writestack, mpistartrow, mpiendrow);
//			double temperSnow_temp, actdiff, maxdiff;
//			do_initializeTemperSnow = true;
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					for (int lv = 0; lv <= levelsSnow; lv++)
					{
						temperSnow[k][j][lv] = nodata;
					}
				}
			}
			do_initializeTemperSnow = true;
			temperSnowindex = gridindex;
			gridindex++;
		}
	}

	/* outflux from the layers (Q-bottom) */
	quSnow = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levelsSnow, nrows);
	quSnow1 = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
	initialize_grid(gridindex, quSnow_stackname, quSnow1, quSnow, levelsSnow, _soil_bottom_fluxes_Snow_, quSnow_writestack, mpistartrow, mpiendrow);
	for (int k = mpistartrow; k <= mpiendrow; k++)
	{
		for (int j = 1; j <= ncols; j++)
		{
			if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
			{
				for (int lv = 1; lv <= levelsSnow; lv++)
				{
					quSnow[k][j][lv] = 0;
				}
				quSnow1[k][j] = 0;
			}
			else
			{
				for (int lv = 0; lv <= levelsSnow; lv++)
				{
					quSnow[k][j][lv] = nodata;
				}
				quSnow1[k][j] = nodata;
			}
		}
	}
	quSnowindex = gridindex;
	gridindex++;

	
}

int unsatzonmodellklasse::allocateGridElementsInRow0(double **grid, unsigned long cnt)
{
	if (grid)
	{
		grid[0] = (double *)calloc((long)(cnt), sizeof(double));
// this shouldnot be done permanently but only temporarily in unsatzonmodell::run(), but not without saving the original value of grid[mpistartrow - 1] in another location and getting this value back at the end of the run() method in order to allow other modules to use the expected values in halos
//		if (mpistartrow > 1) 
//		{	// alias for better organizing the loop in ::run() later
//			grid[mpistartrow - 1] = grid[0]; 
//		} 
		// so for now we only count the number of grids we need to reserve space for the row pointers
		numTempGridRows++;
		return 1;
	}
	return 0;
}

int unsatzonmodellklasse::allocateStackElementsInRow0(double ***stack, unsigned long cnt, int level)
{
	if (stack)
	{
		stack[0] = (double **)calloc((size_t)(cnt), sizeof(double*));
		if (!stack[0])
		{
			(*WSstream) << "FATAL: memory allocation error (2) in allocateStackElementsInRow0()" << "[end-entry]";
			throw(1);
		}
		// assigning the memory linearly to a single temporary variable (will later be used for deleting the entire memory of a row)
		stack[0][0] = (double*)calloc((size_t)(level + 1)*cnt, sizeof(double));

		if (!stack[0][0])
		{
			(*WSstream) << "FATAL: memory allocation error (3) in allocateStackElementsInRow0()" << "[end-entry]";
			throw(1);
		}

		for (int j = 0; j < cnt; j++)
		{
			// instead of allocating memory for every stack, only the pointers are set now
			stack[0][j] = stack[0][0] + j*(level + 1);
		}

// this shouldnot be done permanently but only temporarily in unsatzonmodell::run(), but not without saving the original value of grid[mpistartrow - 1] in another location and getting this value back at the end of the run() method in order to allow other modules to use the expected values in halos
//		if (mpistartrow > 1)
//		{   // alias for better organizing the loop in ::run() later
//			// attention: when these grid/stack rows will be used in other ranks, it can cause errors
//			stack[mpistartrow - 1] = stack[0]; 
//		} 
		// so for now we only count the number of stacks we need to reserve space for the row pointers
		numTempStackRows++;

		return level;
	}
	else
	{
		return 0;
	}
}

// when WaSiM runs as MPI version, Re-infiltration is possible by gathering all river cells in the first few ranks and processing the cells in flow order (from highest to lowest cell, per river link).
// if rank 0 is not big enough (cell count should still not exceed a normal MPI ranks cell count), the river link processing can be split to more than one rank 
// here, the addionally required memory is allocated in row 0 of each grid and stack that is used in unsatzonmodel.
// In addition, in another method the actual values of all of these grids and stacks will be gathered in a buffer and send/received via MPI calls or copied directly if the processing takes place in the same rank
unsigned long unsatzonmodellklasse::prepareMPIReinfiltration()
{
	if (surface_interaction < 2) return 0;
	int cnt = 0;
	for (std::vector< std::vector<RiverLinkCell> >::iterator iter = RivLnkCellsToReceive->begin(); iter < RivLnkCellsToReceive->end();iter++) 
	{
		cnt += iter->size();
	}

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	numAllocatedDoubles = 0;
	if (cnt > 0)
	{
		FirstLastCol[0].First = 0;
		FirstLastCol[0].Last = cnt - 1;

		numAllocatedDoubles += allocateGridElementsInRow0(gridinfo[2].pointer, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ALBEDO            , cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(CLAY				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(DELTAGW			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(DRAINDIST			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(DRAINLEVEL		, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(EBAL				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ETR2ETP			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(EVAP				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(EVAR				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(FDD				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(FROZENDEPTH		, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(FZT				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GW1COLMATION		, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GW1COLMFROZEN		, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GW1COLMUNFROZEN	, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWEX				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWIN				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWLEVEL			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWN				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWSTAND			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(GWtheta			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(HRU				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(HRUMaster			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ICE				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ICEOLD			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(INFEX				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(INFILT			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(IRRIG				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(IRRIGCODE			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(KOLM				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(LAKE				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(MACRO				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(MAXPOND			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(MAXPUMP			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(MPIREG			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(P					, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(PERCOL			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(POND				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QDRAIN			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(P_RAIN			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QDSU				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QDgrid			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QIgrid			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QLOWBND			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(QSNOW_OUT			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(RIVDEP			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(RIVLNK			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(RIVWIT			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ROOT				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(RSE				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(RelCellArea		, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(S0				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SATT				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SBHIGH			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SB_1_				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SB_2_				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SLIQ				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SLP				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SNOW				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SOIL				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SSNO				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SSTO				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SUBSTEPS			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(SURFTEMP			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(T					, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(TANBETA			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(TDAY				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(TDD				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(THAWTMOD			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(THROOT			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(TLOWBND			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(TYEAR				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(UPRISE			, cnt);
//		numAllocatedDoubles += allocateGridElementsInRow0(gridinfo[2].pointer, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(gw_thick			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(qu1				, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(qu1_ext			, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(upperIndexSnow	, cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(zn				, cnt);

		numAllocatedDoubles += allocateGridElementsInRow0(ETP[0], cnt);
		numAllocatedDoubles += allocateGridElementsInRow0(ETR[0], cnt);

		for (int i = 1; i <= iNumVegetationLayers; i++)
		{
			numAllocatedDoubles += allocateGridElementsInRow0(LUSE[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(ETP[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(ETR[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(VBG[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(WURZ[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(LAI[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(RSC[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(RSI[i]	, cnt);
			numAllocatedDoubles += allocateGridElementsInRow0(Z0[i]		, cnt);
		}

		// stacks for soil
		numAllocatedDoubles += allocateStackElementsInRow0(theta	, cnt, levels+1);
		numAllocatedDoubles += allocateStackElementsInRow0(hh		, cnt, levels+1);
		numAllocatedDoubles += allocateStackElementsInRow0(hgeo		, cnt, levels+1);
		numAllocatedDoubles += allocateStackElementsInRow0(qu		, cnt, levels+1);
		numAllocatedDoubles += allocateStackElementsInRow0(temper, cnt, levels + 1);
		numAllocatedDoubles += allocateStackElementsInRow0(SEstack, cnt, levels + 1);

		// stacks for snow
		numAllocatedDoubles += allocateStackElementsInRow0(thetaSnow	, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(hhydSnow		, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(hgeoSnow		, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(quSnow		, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(temperSnow	, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(swes			, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(densSnow		, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(thickSnow	, cnt, levelsSnow+1);
		numAllocatedDoubles += allocateStackElementsInRow0(swel			, cnt, levelsSnow+1);

		numAllocatedDoubles += allocateStackElementsInRow0(ExcessIceFactor, cnt, levels + 1);


		// tracer grids and stacks
		if (tracer_count > 0)
		{
			(*WSstream) << "FATAL: prepareMPIReinfiltration() doesn't allow tracers yet" << "[end-entry]";
			throw(1150);
		}

		// allocating send and receive buffers
		// receive buffer: reserve a buffer for each instance for each instance
		// size is numAllocatedDoubles (which is per cell) times cells from this rank
		RecvBuffer = new double*[mpi_size];
		RecvBufferSize = new unsigned long[mpi_size];
		RecvBufferPos = new unsigned long[mpi_size];
		for (int i = 0; i < mpi_size; i++)
		{
			if ((*RivLnkCellsToReceive)[i].size()>0)
			{
				unsigned long bufsize = numAllocatedDoubles * (*RivLnkCellsToReceive)[i].size();
				RecvBuffer[i] = new double[bufsize];
				RecvBufferSize[i] = bufsize;
			}
			else
			{
				RecvBuffer[i] = NULL;
				RecvBufferSize[i] = 0;
			}
		}

		// allocate memory for the original addresses of the rows that will be used for MPI-Re-Infiltration computation.
		// the rows allocated in index 0 must be indirected to by the row before mpistartrow for using the for loop over all rows of the actual rank
		tempGridRows = new double*[numTempGridRows + 1];
		tempStackRows = new double**[numTempStackRows + 1];
	}

#ifdef _MPI
	// broadcasting numAllocatedDoubles in order to tell all ranks how much buffer needs to be allocated (so sending and receiving ranks are synchronous with respect to memory sent and received)
	MPI_Bcast(&numAllocatedDoubles, 1, MPI_INT, 0,MPI_COMM_WORLD);
#endif

	if (RivLnkCellsToSend->size() > 0)
	{
		// creating send buffer for sending link cell values later to the processing ranks and receiving values back from there
		unsigned long sz = RivLnkCellsToSend->size();
		SendBuffer = new double*[sz];
		SendBufferSize = new unsigned long[sz];
		for (int i = 0; i < sz; i++)
		{
			if ((*RivLnkCellsToSend)[i].size()>0)
			{
				unsigned long bufsize = numAllocatedDoubles * (*RivLnkCellsToSend)[i].size();
				SendBuffer[i] = new double[bufsize];
				SendBufferSize[i] = bufsize;
			}
			else
			{
				SendBuffer[i] = NULL;
				SendBufferSize[i] = 0;
			}
		}
	}

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// put cell values to the processing cells
	MPISyncRivLnkCells1();

//#if defined _MPIDEBUG && defined _MPI
#if defined _MPI

	// delete (or set to 0) all sendbuffers and receive buffers and then populate the recvBuffer again from the processing cells 
	// and broadcast ort send the RecvBuffers to the other ranks again, check sendBuffer then against the original values in the local river link cells
	// if every value is identical, the test is successfully finished
	if (RivLnkCellsToSend->size() > 0)
	{
		// creating send buffer for sending link cell values later to the processing ranks and receiving values back from there
		unsigned long sz = RivLnkCellsToSend->size();
		for (int i = 0; i < sz; i++)
		{
			for (int j = 0; j < SendBufferSize[i]; j++)
			{
				SendBuffer[i][j] = -987654321;
			}
		}
	}

	for (int i = 0; i < mpi_size; i++)
	{
		if ((*RivLnkCellsToReceive)[i].size()>0)
		{
			for (int j = 0; j < RecvBufferSize[i]; j++)
			{
				RecvBuffer[i][j] = 123456789;
			}
		}
	}

	MPISyncRivLnkCells2();

#endif

	return numAllocatedDoubles;
}

// puts cell values from a single cell into the send buffer
unsigned long unsatzonmodellklasse::putCellDataIntoBuffer(int i, int j, double *buffer, unsigned long startpos)
{
	int nextpos = startpos;

	if (gridinfo[2].pointer){ buffer[nextpos] = gridinfo[2].pointer[i][j];	nextpos++; }
	if (ALBEDO)             { buffer[nextpos] = ALBEDO[i][j];				nextpos++; }
	if (CLAY)				{ buffer[nextpos] = CLAY[i][j];					nextpos++; }
	if (DELTAGW)			{ buffer[nextpos] = DELTAGW[i][j];				nextpos++; }
	if (DRAINDIST)			{ buffer[nextpos] = DRAINDIST[i][j];			nextpos++; }
	if (DRAINLEVEL)			{ buffer[nextpos] = DRAINLEVEL[i][j];			nextpos++; }
	if (EBAL)				{ buffer[nextpos] = EBAL[i][j];					nextpos++; }
	if (ETR2ETP)			{ buffer[nextpos] = ETR2ETP[i][j];				nextpos++; }
	if (EVAP)				{ buffer[nextpos] = EVAP[i][j];					nextpos++; }
	if (EVAR)				{ buffer[nextpos] = EVAR[i][j];					nextpos++; }
	if (FDD)				{ buffer[nextpos] = FDD[i][j];					nextpos++; }
	if (FROZENDEPTH)		{ buffer[nextpos] = FROZENDEPTH[i][j];			nextpos++; }
	if (FZT)				{ buffer[nextpos] = FZT[i][j];					nextpos++; }
	if (GW1COLMATION)		{ buffer[nextpos] = GW1COLMATION[i][j];			nextpos++; }
	if (GW1COLMFROZEN)		{ buffer[nextpos] = GW1COLMFROZEN[i][j];		nextpos++; }
	if (GW1COLMUNFROZEN)	{ buffer[nextpos] = GW1COLMUNFROZEN[i][j];		nextpos++; }
	if (GWEX)				{ buffer[nextpos] = GWEX[i][j];					nextpos++; }
	if (GWIN)				{ buffer[nextpos] = GWIN[i][j];					nextpos++; }
	if (GWLEVEL)			{ buffer[nextpos] = GWLEVEL[i][j];				nextpos++; }
	if (GWN)				{ buffer[nextpos] = GWN[i][j];					nextpos++; }
	if (GWSTAND)			{ buffer[nextpos] = GWSTAND[i][j];				nextpos++; }
	if (GWtheta)			{ buffer[nextpos] = GWtheta[i][j];				nextpos++; }
	if (HRU)				{ buffer[nextpos] = HRU[i][j];					nextpos++; }
	if (HRUMaster)			{ buffer[nextpos] = HRUMaster[i][j];			nextpos++; }
	if (ICE)				{ buffer[nextpos] = ICE[i][j];					nextpos++; }
	if (ICEOLD)				{ buffer[nextpos] = ICEOLD[i][j];				nextpos++; }
	if (INFEX)				{ buffer[nextpos] = INFEX[i][j];				nextpos++; }
	if (INFILT)				{ buffer[nextpos] = INFILT[i][j];				nextpos++; }
	if (IRRIG)				{ buffer[nextpos] = IRRIG[i][j];				nextpos++; }
	if (IRRIGCODE)			{ buffer[nextpos] = IRRIGCODE[i][j];			nextpos++; }
	if (KOLM)				{ buffer[nextpos] = KOLM[i][j];					nextpos++; }
	if (LAKE)				{ buffer[nextpos] = LAKE[i][j];					nextpos++; }
	if (MACRO)				{ buffer[nextpos] = MACRO[i][j];				nextpos++; }
	if (MAXPOND)			{ buffer[nextpos] = MAXPOND[i][j];				nextpos++; }
	if (MAXPUMP)			{ buffer[nextpos] = MAXPUMP[i][j];				nextpos++; }
	if (MPIREG)				{ buffer[nextpos] = MPIREG[i][j];				nextpos++; }
	if (P)					{ buffer[nextpos] = P[i][j];					nextpos++; }
	if (PERCOL)				{ buffer[nextpos] = PERCOL[i][j];				nextpos++; }
	if (POND)				{ buffer[nextpos] = POND[i][j];					nextpos++; }
	if (P_RAIN)				{ buffer[nextpos] = P_RAIN[i][j];				nextpos++; }
	if (QDRAIN)				{ buffer[nextpos] = QDRAIN[i][j];				nextpos++; }
	if (QDSU)				{ buffer[nextpos] = QDSU[i][j];					nextpos++; }
	if (QDgrid)				{ buffer[nextpos] = QDgrid[i][j];				nextpos++; }
	if (QIgrid)				{ buffer[nextpos] = QIgrid[i][j];				nextpos++; }
	if (QLOWBND)			{ buffer[nextpos] = QLOWBND[i][j];				nextpos++; }
	if (QSNOW_OUT)			{ buffer[nextpos] = QSNOW_OUT[i][j];			nextpos++; }
	if (RIVDEP)				{ buffer[nextpos] = RIVDEP[i][j];				nextpos++; }
	if (RIVLNK)				{ buffer[nextpos] = RIVLNK[i][j];				nextpos++; }
	if (RIVWIT)				{ buffer[nextpos] = RIVWIT[i][j];				nextpos++; }
	if (ROOT)				{ buffer[nextpos] = ROOT[i][j];					nextpos++; }
	if (RSE)				{ buffer[nextpos] = RSE[i][j];					nextpos++; }
	if (RelCellArea)		{ buffer[nextpos] = RelCellArea[i][j];			nextpos++; }
	if (S0)					{ buffer[nextpos] = S0[i][j];					nextpos++; }
	if (SATT)				{ buffer[nextpos] = SATT[i][j];					nextpos++; }
	if (SBHIGH)				{ buffer[nextpos] = SBHIGH[i][j];				nextpos++; }
	if (SB_1_)				{ buffer[nextpos] = SB_1_[i][j];				nextpos++; }
	if (SB_2_)				{ buffer[nextpos] = SB_2_[i][j];				nextpos++; }
	if (SLIQ)				{ buffer[nextpos] = SLIQ[i][j];					nextpos++; }
	if (SLP)				{ buffer[nextpos] = SLP[i][j];					nextpos++; }
	if (SNOW)				{ buffer[nextpos] = SNOW[i][j];					nextpos++; }
	if (SOIL)				{ buffer[nextpos] = SOIL[i][j];					nextpos++; }
	if (SSNO)				{ buffer[nextpos] = SSNO[i][j];					nextpos++; }
	if (SSTO)				{ buffer[nextpos] = SSTO[i][j];					nextpos++; }
	if (SUBSTEPS)			{ buffer[nextpos] = SUBSTEPS[i][j];				nextpos++; }
	if (SURFTEMP)			{ buffer[nextpos] = SURFTEMP[i][j];				nextpos++; }
	if (T)					{ buffer[nextpos] = T[i][j];					nextpos++; }
	if (TANBETA)			{ buffer[nextpos] = TANBETA[i][j];				nextpos++; }
	if (TDAY)				{ buffer[nextpos] = TDAY[i][j];					nextpos++; }
	if (TDD)				{ buffer[nextpos] = TDD[i][j];					nextpos++; }
	if (THAWTMOD)			{ buffer[nextpos] = THAWTMOD[i][j];				nextpos++; }
	if (THROOT)				{ buffer[nextpos] = THROOT[i][j];				nextpos++; }
	if (TLOWBND)			{ buffer[nextpos] = TLOWBND[i][j];				nextpos++; }
	if (TYEAR)				{ buffer[nextpos] = TYEAR[i][j];				nextpos++; }
	if (UPRISE)				{ buffer[nextpos] = UPRISE[i][j];				nextpos++; }
//	if (gridinfo[2].pointer) { buffer[nextpos] = gridinfo[2].pointer[i][j];	nextpos++; }
	if (gw_thick)			{ buffer[nextpos] = gw_thick[i][j];				nextpos++; }
	if (qu1)				{ buffer[nextpos] = qu1[i][j];					nextpos++; }
	if (qu1_ext)			{ buffer[nextpos] = qu1_ext[i][j];				nextpos++; }
	if (upperIndexSnow)		{ buffer[nextpos] = upperIndexSnow[i][j];		nextpos++; }
	if (zn)					{ buffer[nextpos] = zn[i][j];					nextpos++; }

	if (ETP[0])				{ buffer[nextpos] = ETP[0][i][j];				nextpos++; }
	if (ETR[0])				{ buffer[nextpos] = ETR[0][i][j];				nextpos++; }

	for (int z = 1; z <= iNumVegetationLayers; z++)
	{
		if (LUSE[z]) { buffer[nextpos] = LUSE[z][i][j]; nextpos++; }
		if (ETP[z])  { buffer[nextpos] = ETP[z][i][j];  nextpos++; }
		if (ETR[z])  { buffer[nextpos] = ETR[z][i][j];  nextpos++; }
		if (VBG[z])  { buffer[nextpos] = VBG[z][i][j];  nextpos++; }
		if (WURZ[z]) { buffer[nextpos] = WURZ[z][i][j]; nextpos++; }
		if (LAI[z])  { buffer[nextpos] = LAI[z][i][j];  nextpos++; }
		if (RSC[z])  { buffer[nextpos] = RSC[z][i][j];  nextpos++; }
		if (RSI[z])  { buffer[nextpos] = RSI[z][i][j];  nextpos++; }
		if (Z0[z])   { buffer[nextpos] = Z0[z][i][j];   nextpos++; }
	}

	// stacks for soil
	if (theta)  { memmove((void*)&(buffer[nextpos]), &(theta[i][j][0]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (hh)		{ memmove((void*)&(buffer[nextpos]), &(hh[i][j][0]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (hgeo)	{ memmove((void*)&(buffer[nextpos]), &(hgeo[i][j][0]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (qu)		{ memmove((void*)&(buffer[nextpos]), &(qu[i][j][0]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (temper) { memmove((void*)&(buffer[nextpos]), &(temper[i][j][0]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }
	if (SEstack) { memmove((void*)&(buffer[nextpos]), &(SEstack[i][j][0]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }

	// stacks for snow
	if (thetaSnow)  { memmove((void*)&(buffer[nextpos]), &(thetaSnow[i][j][0]), (levelsSnow + 1)*sizeof(double));  nextpos += (levelsSnow + 1); }
	if (hhydSnow)	{ memmove((void*)&(buffer[nextpos]), &(hhydSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (hgeoSnow)	{ memmove((void*)&(buffer[nextpos]), &(hgeoSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (quSnow)		{ memmove((void*)&(buffer[nextpos]), &(quSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (temperSnow)	{ memmove((void*)&(buffer[nextpos]), &(temperSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (swes)		{ memmove((void*)&(buffer[nextpos]), &(swes[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (densSnow)	{ memmove((void*)&(buffer[nextpos]), &(densSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (thickSnow)	{ memmove((void*)&(buffer[nextpos]), &(thickSnow[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (swel)		{ memmove((void*)&(buffer[nextpos]), &(swel[i][j][0]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }

	if (ExcessIceFactor) { memmove((void*)&(buffer[nextpos]), &(ExcessIceFactor[i][j][0]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }

	// todo: tracer grids and stacks

	return nextpos;
}

unsigned long unsatzonmodellklasse::getCellDataFromBuffer(int i,int j,double *buffer, unsigned long nextpos)
{

	if (gridinfo[2].pointer){ gridinfo[2].pointer[i][j] = buffer[nextpos];	nextpos++; }
	if (ALBEDO)             { ALBEDO[i][j]              = buffer[nextpos];	nextpos++; }
	if (CLAY)				{ CLAY[i][j]				= buffer[nextpos];	nextpos++; }
	if (DELTAGW)			{ DELTAGW[i][j]				= buffer[nextpos];	nextpos++; }
	if (DRAINDIST)			{ DRAINDIST[i][j]			= buffer[nextpos];	nextpos++; }
	if (DRAINLEVEL)			{ DRAINLEVEL[i][j]			= buffer[nextpos];	nextpos++; }
	if (EBAL)				{ EBAL[i][j]				= buffer[nextpos];	nextpos++; }
	if (ETR2ETP)			{ ETR2ETP[i][j]				= buffer[nextpos];	nextpos++; }
	if (EVAP)				{ EVAP[i][j]				= buffer[nextpos];	nextpos++; }
	if (EVAR)				{ EVAR[i][j]				= buffer[nextpos];	nextpos++; }
	if (FDD)				{ FDD[i][j]					= buffer[nextpos];	nextpos++; }
	if (FROZENDEPTH)		{ FROZENDEPTH[i][j]			= buffer[nextpos];	nextpos++; }
	if (FZT)				{ FZT[i][j]					= buffer[nextpos];	nextpos++; }
	if (GW1COLMATION)		{ GW1COLMATION[i][j]		= buffer[nextpos];	nextpos++; }
	if (GW1COLMFROZEN)		{ GW1COLMFROZEN[i][j]		= buffer[nextpos];	nextpos++; }
	if (GW1COLMUNFROZEN)	{ GW1COLMUNFROZEN[i][j]		= buffer[nextpos];	nextpos++; }
	if (GWEX)				{ GWEX[i][j]				= buffer[nextpos];	nextpos++; }
	if (GWIN)				{ GWIN[i][j]				= buffer[nextpos];	nextpos++; }
	if (GWLEVEL)			{ GWLEVEL[i][j]				= buffer[nextpos];	nextpos++; }
	if (GWN)				{ GWN[i][j]					= buffer[nextpos];	nextpos++; }
	if (GWSTAND)			{ GWSTAND[i][j]				= buffer[nextpos];	nextpos++; }
	if (GWtheta)			{ GWtheta[i][j]				= buffer[nextpos];	nextpos++; }
	if (HRU)				{ HRU[i][j]					= buffer[nextpos];	nextpos++; }
	if (HRUMaster)			{ HRUMaster[i][j]			= buffer[nextpos];	nextpos++; }
	if (ICE)				{ ICE[i][j]					= buffer[nextpos];	nextpos++; }
	if (ICEOLD)				{ ICEOLD[i][j]				= buffer[nextpos];	nextpos++; }
	if (INFEX)				{ INFEX[i][j]				= buffer[nextpos];	nextpos++; }
	if (INFILT)				{ INFILT[i][j]				= buffer[nextpos];	nextpos++; }
	if (IRRIG)				{ IRRIG[i][j]				= buffer[nextpos];	nextpos++; }
	if (IRRIGCODE)			{ IRRIGCODE[i][j]			= buffer[nextpos];	nextpos++; }
	if (KOLM)				{ KOLM[i][j]				= buffer[nextpos];	nextpos++; }
	if (LAKE)				{ LAKE[i][j]				= buffer[nextpos];	nextpos++; }
	if (MACRO)				{ MACRO[i][j]				= buffer[nextpos];	nextpos++; }
	if (MAXPOND)			{ MAXPOND[i][j]				= buffer[nextpos];	nextpos++; }
	if (MAXPUMP)			{ MAXPUMP[i][j]				= buffer[nextpos];	nextpos++; }
	if (MPIREG)				{ MPIREG[i][j]				= buffer[nextpos];	nextpos++; }
	if (P)					{ P[i][j]					= buffer[nextpos];	nextpos++; }
	if (PERCOL)				{ PERCOL[i][j]				= buffer[nextpos];	nextpos++; }
	if (POND)				{ POND[i][j]				= buffer[nextpos];	nextpos++; }
	if (P_RAIN)				{ P_RAIN[i][j]				= buffer[nextpos];	nextpos++; }
	if (QDRAIN)				{ QDRAIN[i][j]				= buffer[nextpos];	nextpos++; }
	if (QDSU)				{ QDSU[i][j]				= buffer[nextpos];	nextpos++; }
	if (QDgrid)				{ QDgrid[i][j]				= buffer[nextpos];	nextpos++; }
	if (QIgrid)				{ QIgrid[i][j]				= buffer[nextpos];	nextpos++; }
	if (QLOWBND)			{ QLOWBND[i][j]				= buffer[nextpos];	nextpos++; }
	if (QSNOW_OUT)			{ QSNOW_OUT[i][j]			= buffer[nextpos];	nextpos++; }
	if (RIVDEP)				{ RIVDEP[i][j]				= buffer[nextpos];	nextpos++; }
	if (RIVLNK)				{ RIVLNK[i][j]				= buffer[nextpos];	nextpos++; }
	if (RIVWIT)				{ RIVWIT[i][j]				= buffer[nextpos];	nextpos++; }
	if (ROOT)				{ ROOT[i][j]				= buffer[nextpos];	nextpos++; }
	if (RSE)				{ RSE[i][j]					= buffer[nextpos];	nextpos++; }
	if (RelCellArea)		{ RelCellArea[i][j]			= buffer[nextpos];	nextpos++; }
	if (S0)					{ S0[i][j]					= buffer[nextpos];	nextpos++; }
	if (SATT)				{ SATT[i][j]				= buffer[nextpos];	nextpos++; }
	if (SBHIGH)				{ SBHIGH[i][j]				= buffer[nextpos];	nextpos++; }
	if (SB_1_)				{ SB_1_[i][j]				= buffer[nextpos];	nextpos++; }
	if (SB_2_)				{ SB_2_[i][j]				= buffer[nextpos];	nextpos++; }
	if (SLIQ)				{ SLIQ[i][j]				= buffer[nextpos];	nextpos++; }
	if (SLP)				{ SLP[i][j]					= buffer[nextpos];	nextpos++; }
	if (SNOW)				{ SNOW[i][j]				= buffer[nextpos];	nextpos++; }
	if (SOIL)				{ SOIL[i][j]				= buffer[nextpos];	nextpos++; }
	if (SSNO)				{ SSNO[i][j]				= buffer[nextpos];	nextpos++; }
	if (SSTO)				{ SSTO[i][j]				= buffer[nextpos];	nextpos++; }
	if (SUBSTEPS)			{ SUBSTEPS[i][j]			= buffer[nextpos];	nextpos++; }
	if (SURFTEMP)			{ SURFTEMP[i][j]			= buffer[nextpos];	nextpos++; }
	if (T)					{ T[i][j]					= buffer[nextpos];	nextpos++; }
	if (TANBETA)			{ TANBETA[i][j]				= buffer[nextpos];	nextpos++; }
	if (TDAY)				{ TDAY[i][j]				= buffer[nextpos];	nextpos++; }
	if (TDD)				{ TDD[i][j]					= buffer[nextpos];	nextpos++; }
	if (THAWTMOD)			{ THAWTMOD[i][j]			= buffer[nextpos];	nextpos++; }
	if (THROOT)				{ THROOT[i][j]				= buffer[nextpos];	nextpos++; }
	if (TLOWBND)			{ TLOWBND[i][j]				= buffer[nextpos];	nextpos++; }
	if (TYEAR)				{ TYEAR[i][j]				= buffer[nextpos];	nextpos++; }
	if (UPRISE)				{ UPRISE[i][j]				= buffer[nextpos];	nextpos++; }
//	if (gridinfo[2].pointer){ gridinfo[2].pointer[i][j]	= buffer[nextpos];	nextpos++; }
	if (gw_thick)			{ gw_thick[i][j]			= buffer[nextpos];	nextpos++; }
	if (qu1)				{ qu1[i][j]					= buffer[nextpos];	nextpos++; }
	if (qu1_ext)			{ qu1_ext[i][j]				= buffer[nextpos];	nextpos++; }
	if (upperIndexSnow)		{ upperIndexSnow[i][j]		= buffer[nextpos];	nextpos++; }
	if (zn)					{ zn[i][j]					= buffer[nextpos];	nextpos++; }

	if (ETP[0])				{ ETP[0][i][j]				= buffer[nextpos];	nextpos++; }
	if (ETR[0])				{ ETR[0][i][j]				= buffer[nextpos];	nextpos++; }

	for (int z = 1; z <= iNumVegetationLayers; z++)
	{
		if (LUSE[z]) { LUSE[z][i][j]= buffer[nextpos]; nextpos++; }
		if (ETP[z])  { ETP[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (ETR[z])  { ETR[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (VBG[z])  { VBG[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (WURZ[z]) { WURZ[z][i][j]= buffer[nextpos]; nextpos++; }
		if (LAI[z])  { LAI[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (RSC[z])  { RSC[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (RSI[z])  { RSI[z][i][j]	= buffer[nextpos]; nextpos++; }
		if (Z0[z])   { Z0[z][i][j]	= buffer[nextpos]; nextpos++; }
	}

	// stacks for soil
	if (theta)  { memmove(&(theta[i][j][0])	, (void*)&(buffer[nextpos]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (hh)		{ memmove(&(hh[i][j][0])	, (void*)&(buffer[nextpos]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (hgeo)	{ memmove(&(hgeo[i][j][0])	, (void*)&(buffer[nextpos]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (qu)		{ memmove(&(qu[i][j][0])	, (void*)&(buffer[nextpos]), (levels + 1)*sizeof(double)); nextpos += (levels + 1); }
	if (temper) { memmove(&(temper[i][j][0]), (void*)&(buffer[nextpos]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }
	if (SEstack) { memmove(&(SEstack[i][j][0]), (void*)&(buffer[nextpos]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }

	// stacks for snow
	if (thetaSnow)  { memmove(&(thetaSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (hhydSnow)	{ memmove(&(hhydSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (hgeoSnow)	{ memmove(&(hgeoSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (quSnow)		{ memmove(&(quSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (temperSnow)	{ memmove(&(temperSnow[i][j][0]), (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (swes)		{ memmove(&(swes[i][j][0])		, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (densSnow)	{ memmove(&(densSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (thickSnow)	{ memmove(&(thickSnow[i][j][0])	, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	if (swel)		{ memmove(&(swel[i][j][0])		, (void*)&(buffer[nextpos]), (levelsSnow + 1)*sizeof(double)); nextpos += (levelsSnow + 1); }
	
	if (ExcessIceFactor) { memmove(&(ExcessIceFactor[i][j][0]), (void*)&(buffer[nextpos]), (levels + 1) * sizeof(double)); nextpos += (levels + 1); }

	// todo: tracer grids and stacks

	return nextpos;
}


// function could be used in MPISyncRivLnkCells1, but is now implemented inline there
unsigned long unsatzonmodellklasse::fillSendBuffer(std::vector<RiverLinkCell> &RivLnkVector, double *buffer, unsigned long startpos)
{
	unsigned long nextpos = startpos;
	for (std::vector<RiverLinkCell>::iterator iter = RivLnkVector.begin(); iter != RivLnkVector.end(); iter++)
	{
		nextpos = putCellDataIntoBuffer(iter->i, iter->j, buffer, nextpos);
	}
	return nextpos;
}

void unsatzonmodellklasse::MPISyncRivLnkCells1()
{
	// filling send buffers now for the first time
	for (int i = 0; i < (*RivLnkCellsToSend).size(); i++)
	{  // RivLnkCellsToSend contains one vector per processing instance, i.e. usually only one or two entries at all.
		// those vectors, however, contain ordered cells for several routing links
		// that should be put together in the correct order in the processing rank
		unsigned long nextpos = 0;
		// fillSendBuffer((*RivLnkCellsToSend)[i], SendBuffer[i], pos);
		for (std::vector<RiverLinkCell>::iterator iter = (*RivLnkCellsToSend)[i].begin(); iter != (*RivLnkCellsToSend)[i].end(); iter++)
		{
			nextpos = putCellDataIntoBuffer(iter->i, iter->j, SendBuffer[i], nextpos);
		}
	}

#ifdef _MPI

	MPI_Request *requestVector;
	requestVector = (MPI_Request*)calloc((*RivLnkCellsToSend).size() + (*RivLnkCellsToReceive).size(), sizeof(MPI_Request));
	int request_cnt = 0;

	// if the send rank is the sme as the receive rank, the memory can simply be copied as well, but this is not mandatory, since a good MPI implementation will do this anyway and the code becomes more unreadable
	if (SendBufferSize != NULL)
	{
		for (int i = 0; i < (*RivLnkCellsToSend).size(); i++)
		{
			// sending 
			if (SendBufferSize[i] > 0)
			{
				int destrank = (*RivLnkCellsToSend)[i].begin()->processingrank;
				MPI_Isend(&(SendBuffer[i][0]), SendBufferSize[i], MPI_DOUBLE, destrank, 0,MPI_COMM_WORLD,&(requestVector[request_cnt]));
				request_cnt++;
			}
		}
	}

	if (RecvBufferSize != NULL)
	{
		for (int i = 0; i < (*RivLnkCellsToReceive).size(); i++)
		{
			// sending 
			if (RecvBufferSize[i] > 0)
			{
				int sourcerank = (*RivLnkCellsToReceive)[i].begin()->mpirank;
				MPI_Irecv(&(RecvBuffer[i][0]), RecvBufferSize[i], MPI_DOUBLE, sourcerank, 0,MPI_COMM_WORLD,&(requestVector[request_cnt]));
				request_cnt++;
			}
		}
	}

	if (request_cnt > 0)
	{
		MPI_Waitall(request_cnt, requestVector, statusVector);
	}
	free(requestVector);

#else
	// if not in MPI mode, the SendBuffer can simply be copied to the recvBuffer
	memcpy((void*)&(RecvBuffer[0][0]), (void*)&(SendBuffer[0][0]), (unsigned long)(SendBufferSize[0]*sizeof(double)));
#endif

	// When copying the buffer data after MPI transfer into the processing cells, the buffers are processed out-of-order (as defined in the RiverLinkVector)
	// that means, tht the actual position for each buffer must somewhere stored, which is the RecvBufferPos-array. This is now initialized to the first element of each buffer
	if (RecvBufferPos)
	{
		for (int i = 0; i < mpi_size; i++)
		{
			RecvBufferPos[i] = 0;
		}
	}

	// now extracting cell values from the receive buffer and putting them into the processing cells in row 0 of each grid
	int j = 0;
	for (int r = 0; r < (*RiverLinkVector).size(); r++)
	{
		if ((*RiverLinkVector)[r]->linkProcessingrank == mpi_rank) // only cells for the actual rank will be moved to the processing cells
		{
			// processing each cell: each riverLinkCell contains a member that tells us the source rank, so we know where to get the data from
			// for (int c = 0; c < (*RiverLinkVector)[r]->linkCellVector->size(); c++)
			for (std::vector<RiverLinkCell>::iterator iter = (*RiverLinkVector)[r]->linkCellVector->begin(); iter != (*RiverLinkVector)[r]->linkCellVector->end(); iter++)
			{
				int sourcerank = iter->mpirank;
				RecvBufferPos[sourcerank] = getCellDataFromBuffer(0, j, RecvBuffer[sourcerank], RecvBufferPos[sourcerank]);
				j++;
			}
		}
	}

	SetTemporaryRowPtrs();
// must be called before the sub steps are cicled through, i.e. in pre_run, not here 
//	MPIBCastRiverFlowTable(); // all ranks need to know how much water they can use for their infiltration links
}


void unsatzonmodellklasse::MPISyncRivLnkCells2() 
//todo: hier noch die Reihenfolge der Buffer umkehren bzw. die Werte aus den processing cells nun in den RecvBuffer schreiben , dann MPI exchange durchführen und dann die Werte in den einzelnen Ranks aus dem Sendbuffer in die lokalen Zellen kopieren bzw. Vergleichen mit dem Originalwert
{
	// filling RevcBuffers with cell values for back trensfer
	// When copying the buffer data from the processing cells, the buffers are processed out-of-order (as defined in the RiverLinkVector)
	// that means, tht the actual position for each buffer must be stored somewhere, which is the RecvBufferPos-array. This is now initialized to the first element of each buffer
	if (RecvBufferPos)
	{
		for (int i = 0; i < mpi_size; i++)
		{
			RecvBufferPos[i] = 0;
		}
	}

	// now getting values from the processing cells in row 0 of each grid and putting them into the receive buffer 
	int j = 0;
	for (int r = 0; r < (*RiverLinkVector).size(); r++)
	{
		if ((*RiverLinkVector)[r]->linkProcessingrank == mpi_rank) // only cells for the actual rank will be moved to the processing cells
		{
			// processing each cell: each riverLinkCell contains a member that tells us the source rank, so we know where to get the data from
			// for (int c = 0; c < (*RiverLinkVector)[r]->linkCellVector->size(); c++)
			for (std::vector<RiverLinkCell>::iterator iter = (*RiverLinkVector)[r]->linkCellVector->begin(); iter != (*RiverLinkVector)[r]->linkCellVector->end(); iter++)
			{
				int targetrank = iter->mpirank;
				RecvBufferPos[targetrank] = putCellDataIntoBuffer(0, j, RecvBuffer[targetrank], RecvBufferPos[targetrank]);
				j++;
			}
		}
	}

#ifdef _MPI
	MPI_Request *requestVector;
	requestVector = (MPI_Request*)calloc((*RivLnkCellsToSend).size() + (*RivLnkCellsToReceive).size(), sizeof(MPI_Request));
	int request_cnt = 0;

	if (RecvBufferSize != NULL)
	{
		for (int i = 0; i < (*RivLnkCellsToReceive).size(); i++)
		{
			// sending 
			if (RecvBufferSize[i] > 0)
			{
				int targetrank = (*RivLnkCellsToReceive)[i].begin()->mpirank;
				MPI_Isend(&(RecvBuffer[i][0]), RecvBufferSize[i], MPI_DOUBLE, targetrank, 0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
				request_cnt++;
			}
		}
	}

	// if the send rank is the sme as the receive rank, the memory can simply be copied as well, but this is not mandatory, since a good MPI implementation will do this anyway and the code becomes more unreadable
	if (SendBufferSize != NULL)
	{
		for (int i = 0; i < (*RivLnkCellsToSend).size(); i++)
		{
			// sending 
			if (SendBufferSize[i] > 0)
			{
				int sourcerank = (*RivLnkCellsToSend)[i].begin()->processingrank;
				MPI_Irecv(&(SendBuffer[i][0]), SendBufferSize[i], MPI_DOUBLE, sourcerank, 0, MPI_COMM_WORLD, &(requestVector[request_cnt]));
				request_cnt++;
			}
		}
	}

	if (request_cnt > 0)
	{
		MPI_Waitall(request_cnt, requestVector, statusVector);
	}
	free(requestVector);

#else
	// if not in MPI mode, the SendBuffer can simply be copied to the recvBuffer
	memcpy((void*)&(SendBuffer[0][0]), (void*)&(RecvBuffer[0][0]), (unsigned long)(SendBufferSize[0] * sizeof(double)));
#endif


	// extracting the balues from the sendBuffers into the local rank's cells
	for (int i = 0; i < (*RivLnkCellsToSend).size(); i++)
	{  // RivLnkCellsToSend contains one vector per processing instance, i.e. usually only one or two entries at all.
		// those vectors, however, contain ordered cells for several routing links
		// that should be put together in the correct order in the processing rank

//		fillSendBuffer((*RivLnkCellsToSend)[i], SendBuffer[i], pos);
		unsigned long nextpos = 0;
		for (std::vector<RiverLinkCell>::iterator iter = (*RivLnkCellsToSend)[i].begin(); iter != (*RivLnkCellsToSend)[i].end(); iter++)
		{
			nextpos = getCellDataFromBuffer(iter->i, iter->j, SendBuffer[i], nextpos);
		}
	}

	RecoverOrigRowPtrs();
// must be called in post_run in order to get the results only once gathered when multiple sub steps are computed
//	MPIGatherRiverFlowTable(); // all re-infiltrated water needs to be reported back to the routing module in rank 0 in order to let that instance know how much water is available as input into other river links in the next time step
	MPIPrepareRiverFlowTable();
}

void unsatzonmodellklasse::SetTemporaryRowPtrs()
{
	if (tempGridRows)
	{
		int nextpos = 0;
		if (gridinfo[2].pointer){ tempGridRows[nextpos] = gridinfo[2].pointer[mpistartrow - 1];	gridinfo[2].pointer[mpistartrow - 1] = gridinfo[2].pointer[0];	nextpos++; }
		if (ALBEDO)             { tempGridRows[nextpos] = ALBEDO[mpistartrow - 1];				ALBEDO[mpistartrow - 1]             = ALBEDO[0];				nextpos++; }
		if (CLAY)				{ tempGridRows[nextpos] = CLAY[mpistartrow - 1];				CLAY[mpistartrow - 1]				= CLAY[0];					nextpos++; }
		if (DELTAGW)			{ tempGridRows[nextpos] = DELTAGW[mpistartrow - 1];				DELTAGW[mpistartrow - 1]			= DELTAGW[0];				nextpos++; }
		if (DRAINDIST)			{ tempGridRows[nextpos] = DRAINDIST[mpistartrow - 1];			DRAINDIST[mpistartrow - 1]			= DRAINDIST[0];				nextpos++; }
		if (DRAINLEVEL)			{ tempGridRows[nextpos] = DRAINLEVEL[mpistartrow - 1];			DRAINLEVEL[mpistartrow - 1]			= DRAINLEVEL[0];			nextpos++; }
		if (EBAL)				{ tempGridRows[nextpos] = EBAL[mpistartrow - 1];				EBAL[mpistartrow - 1]				= EBAL[0];					nextpos++; }
		if (ETR2ETP)			{ tempGridRows[nextpos] = ETR2ETP[mpistartrow - 1];				ETR2ETP[mpistartrow - 1]			= ETR2ETP[0];				nextpos++; }
		if (EVAP)				{ tempGridRows[nextpos] = EVAP[mpistartrow - 1];				EVAP[mpistartrow - 1]				= EVAP[0];					nextpos++; }
		if (EVAR)				{ tempGridRows[nextpos] = EVAR[mpistartrow - 1];				EVAR[mpistartrow - 1]				= EVAR[0];					nextpos++; }
		if (FDD)				{ tempGridRows[nextpos] = FDD[mpistartrow - 1];					FDD[mpistartrow - 1]				= FDD[0];					nextpos++; }
		if (FROZENDEPTH)		{ tempGridRows[nextpos] = FROZENDEPTH[mpistartrow - 1];			FROZENDEPTH[mpistartrow - 1]		= FROZENDEPTH[0];			nextpos++; }
		if (FZT)				{ tempGridRows[nextpos] = FZT[mpistartrow - 1];					FZT[mpistartrow - 1]				= FZT[0];					nextpos++; }
		if (GW1COLMATION)		{ tempGridRows[nextpos] = GW1COLMATION[mpistartrow - 1];		GW1COLMATION[mpistartrow - 1]		= GW1COLMATION[0];			nextpos++; }
		if (GW1COLMFROZEN)		{ tempGridRows[nextpos] = GW1COLMFROZEN[mpistartrow - 1];		GW1COLMFROZEN[mpistartrow - 1]		= GW1COLMFROZEN[0];			nextpos++; }
		if (GW1COLMUNFROZEN)	{ tempGridRows[nextpos] = GW1COLMUNFROZEN[mpistartrow - 1];		GW1COLMUNFROZEN[mpistartrow - 1]	= GW1COLMUNFROZEN[0];		nextpos++; }
		if (GWEX)				{ tempGridRows[nextpos] = GWEX[mpistartrow - 1];				GWEX[mpistartrow - 1]				= GWEX[0];					nextpos++; }
		if (GWIN)				{ tempGridRows[nextpos] = GWIN[mpistartrow - 1];				GWIN[mpistartrow - 1]				= GWIN[0];					nextpos++; }
		if (GWLEVEL)			{ tempGridRows[nextpos] = GWLEVEL[mpistartrow - 1];				GWLEVEL[mpistartrow - 1]			= GWLEVEL[0];				nextpos++; }
		if (GWN)				{ tempGridRows[nextpos] = GWN[mpistartrow - 1];					GWN[mpistartrow - 1]				= GWN[0];					nextpos++; }
		if (GWSTAND)			{ tempGridRows[nextpos] = GWSTAND[mpistartrow - 1];				GWSTAND[mpistartrow - 1]			= GWSTAND[0];				nextpos++; }
		if (GWtheta)			{ tempGridRows[nextpos] = GWtheta[mpistartrow - 1];				GWtheta[mpistartrow - 1]			= GWtheta[0];				nextpos++; }
		if (HRU)				{ tempGridRows[nextpos] = HRU[mpistartrow - 1];					HRU[mpistartrow - 1]				= HRU[0];					nextpos++; }
		if (HRUMaster)			{ tempGridRows[nextpos] = HRUMaster[mpistartrow - 1];			HRUMaster[mpistartrow - 1]			= HRUMaster[0];				nextpos++; }
		if (ICE)				{ tempGridRows[nextpos] = ICE[mpistartrow - 1];					ICE[mpistartrow - 1]				= ICE[0];					nextpos++; }
		if (ICEOLD)				{ tempGridRows[nextpos] = ICEOLD[mpistartrow - 1];				ICEOLD[mpistartrow - 1]				= ICEOLD[0];				nextpos++; }
		if (INFEX)				{ tempGridRows[nextpos] = INFEX[mpistartrow - 1];				INFEX[mpistartrow - 1]				= INFEX[0];					nextpos++; }
		if (INFILT)				{ tempGridRows[nextpos] = INFILT[mpistartrow - 1];				INFILT[mpistartrow - 1]				= INFILT[0];				nextpos++; }
		if (IRRIG)				{ tempGridRows[nextpos] = IRRIG[mpistartrow - 1];				IRRIG[mpistartrow - 1]				= IRRIG[0];					nextpos++; }
		if (IRRIGCODE)			{ tempGridRows[nextpos] = IRRIGCODE[mpistartrow - 1];			IRRIGCODE[mpistartrow - 1]			= IRRIGCODE[0];				nextpos++; }
		if (KOLM)				{ tempGridRows[nextpos] = KOLM[mpistartrow - 1];				KOLM[mpistartrow - 1]				= KOLM[0];					nextpos++; }
		if (LAKE)				{ tempGridRows[nextpos] = LAKE[mpistartrow - 1];				LAKE[mpistartrow - 1]				= LAKE[0];					nextpos++; }
		if (MACRO)				{ tempGridRows[nextpos] = MACRO[mpistartrow - 1];				MACRO[mpistartrow - 1]				= MACRO[0];					nextpos++; }
		if (MAXPOND)			{ tempGridRows[nextpos] = MAXPOND[mpistartrow - 1];				MAXPOND[mpistartrow - 1]			= MAXPOND[0];				nextpos++; }
		if (MAXPUMP)			{ tempGridRows[nextpos] = MAXPUMP[mpistartrow - 1];				MAXPUMP[mpistartrow - 1]			= MAXPUMP[0];				nextpos++; }
		if (MPIREG)				{ tempGridRows[nextpos] = MPIREG[mpistartrow - 1];				MPIREG[mpistartrow - 1]				= MPIREG[0];				nextpos++; }
		if (P)					{ tempGridRows[nextpos] = P[mpistartrow - 1];					P[mpistartrow - 1]					= P[0];						nextpos++; }
		if (PERCOL)				{ tempGridRows[nextpos] = PERCOL[mpistartrow - 1];				PERCOL[mpistartrow - 1]				= PERCOL[0];				nextpos++; }
		if (POND)				{ tempGridRows[nextpos] = POND[mpistartrow - 1];				POND[mpistartrow - 1]				= POND[0];					nextpos++; }
		if (P_RAIN)				{ tempGridRows[nextpos] = P_RAIN[mpistartrow - 1];				P_RAIN[mpistartrow - 1]				= P_RAIN[0];				nextpos++; }
		if (QDRAIN)				{ tempGridRows[nextpos] = QDRAIN[mpistartrow - 1];				QDRAIN[mpistartrow - 1]				= QDRAIN[0];				nextpos++; }
		if (QDSU)				{ tempGridRows[nextpos] = QDSU[mpistartrow - 1];				QDSU[mpistartrow - 1]				= QDSU[0];					nextpos++; }
		if (QDgrid)				{ tempGridRows[nextpos] = QDgrid[mpistartrow - 1];				QDgrid[mpistartrow - 1]				= QDgrid[0];				nextpos++; }
		if (QIgrid)				{ tempGridRows[nextpos] = QIgrid[mpistartrow - 1];				QIgrid[mpistartrow - 1]				= QIgrid[0];				nextpos++; }
		if (QLOWBND)			{ tempGridRows[nextpos] = QLOWBND[mpistartrow - 1];				QLOWBND[mpistartrow - 1]			= QLOWBND[0];				nextpos++; }
		if (QSNOW_OUT)			{ tempGridRows[nextpos] = QSNOW_OUT[mpistartrow - 1];			QSNOW_OUT[mpistartrow - 1]			= QSNOW_OUT[0];				nextpos++; }
		if (RIVDEP)				{ tempGridRows[nextpos] = RIVDEP[mpistartrow - 1];				RIVDEP[mpistartrow - 1]				= RIVDEP[0];				nextpos++; }
		if (RIVLNK)				{ tempGridRows[nextpos] = RIVLNK[mpistartrow - 1];				RIVLNK[mpistartrow - 1]				= RIVLNK[0];				nextpos++; }
		if (RIVWIT)				{ tempGridRows[nextpos] = RIVWIT[mpistartrow - 1];				RIVWIT[mpistartrow - 1]				= RIVWIT[0];				nextpos++; }
		if (ROOT)				{ tempGridRows[nextpos] = ROOT[mpistartrow - 1];				ROOT[mpistartrow - 1]				= ROOT[0];					nextpos++; }
		if (RSE)				{ tempGridRows[nextpos] = RSE[mpistartrow - 1];					RSE[mpistartrow - 1]				= RSE[0];					nextpos++; }
		if (RelCellArea)		{ tempGridRows[nextpos] = RelCellArea[mpistartrow - 1];			RelCellArea[mpistartrow - 1]		= RelCellArea[0];			nextpos++; }
		if (S0)					{ tempGridRows[nextpos] = S0[mpistartrow - 1];					S0[mpistartrow - 1]					= S0[0];					nextpos++; }
		if (SATT)				{ tempGridRows[nextpos] = SATT[mpistartrow - 1];				SATT[mpistartrow - 1]				= SATT[0];					nextpos++; }
		if (SBHIGH)				{ tempGridRows[nextpos] = SBHIGH[mpistartrow - 1];				SBHIGH[mpistartrow - 1]				= SBHIGH[0];				nextpos++; }
		if (SB_1_)				{ tempGridRows[nextpos] = SB_1_[mpistartrow - 1];				SB_1_[mpistartrow - 1]				= SB_1_[0];					nextpos++; }
		if (SB_2_)				{ tempGridRows[nextpos] = SB_2_[mpistartrow - 1];				SB_2_[mpistartrow - 1]				= SB_2_[0];					nextpos++; }
		if (SLIQ)				{ tempGridRows[nextpos] = SLIQ[mpistartrow - 1];				SLIQ[mpistartrow - 1]				= SLIQ[0];					nextpos++; }
		if (SLP)				{ tempGridRows[nextpos] = SLP[mpistartrow - 1];					SLP[mpistartrow - 1]				= SLP[0];					nextpos++; }
		if (SNOW)				{ tempGridRows[nextpos] = SNOW[mpistartrow - 1];				SNOW[mpistartrow - 1]				= SNOW[0];					nextpos++; }
		if (SOIL)				{ tempGridRows[nextpos] = SOIL[mpistartrow - 1];				SOIL[mpistartrow - 1]				= SOIL[0];					nextpos++; }
		if (SSNO)				{ tempGridRows[nextpos] = SSNO[mpistartrow - 1];				SSNO[mpistartrow - 1]				= SSNO[0];					nextpos++; }
		if (SSTO)				{ tempGridRows[nextpos] = SSTO[mpistartrow - 1];				SSTO[mpistartrow - 1]				= SSTO[0];					nextpos++; }
		if (SUBSTEPS)			{ tempGridRows[nextpos] = SUBSTEPS[mpistartrow - 1];			SUBSTEPS[mpistartrow - 1]			= SUBSTEPS[0];				nextpos++; }
		if (SURFTEMP)			{ tempGridRows[nextpos] = SURFTEMP[mpistartrow - 1];			SURFTEMP[mpistartrow - 1]			= SURFTEMP[0];				nextpos++; }
		if (T)					{ tempGridRows[nextpos] = T[mpistartrow - 1];					T[mpistartrow - 1]					= T[0];						nextpos++; }
		if (TANBETA)			{ tempGridRows[nextpos] = TANBETA[mpistartrow - 1];				TANBETA[mpistartrow - 1]			= TANBETA[0];				nextpos++; }
		if (TDAY)				{ tempGridRows[nextpos] = TDAY[mpistartrow - 1];				TDAY[mpistartrow - 1]				= TDAY[0];					nextpos++; }
		if (TDD)				{ tempGridRows[nextpos] = TDD[mpistartrow - 1];					TDD[mpistartrow - 1]				= TDD[0];					nextpos++; }
		if (THAWTMOD)			{ tempGridRows[nextpos] = THAWTMOD[mpistartrow - 1];			THAWTMOD[mpistartrow - 1]			= THAWTMOD[0];				nextpos++; }
		if (THROOT)				{ tempGridRows[nextpos] = THROOT[mpistartrow - 1];				THROOT[mpistartrow - 1]				= THROOT[0];				nextpos++; }
		if (TLOWBND)			{ tempGridRows[nextpos] = TLOWBND[mpistartrow - 1];				TLOWBND[mpistartrow - 1]			= TLOWBND[0];				nextpos++; }
		if (TYEAR)				{ tempGridRows[nextpos] = TYEAR[mpistartrow - 1];				TYEAR[mpistartrow - 1]				= TYEAR[0];					nextpos++; }
		if (UPRISE)				{ tempGridRows[nextpos] = UPRISE[mpistartrow - 1];				UPRISE[mpistartrow - 1]				= UPRISE[0];				nextpos++; }
//		if (gridinfo[2].pointer){ tempGridRows[nextpos] = gridinfo[2].pointer[mpistartrow - 1];	gridinfo[2].pointer[mpistartrow - 1]= gridinfo[2].pointer[0];	nextpos++; }
		if (gw_thick)			{ tempGridRows[nextpos] = gw_thick[mpistartrow - 1];			gw_thick[mpistartrow - 1]			= gw_thick[0];				nextpos++; }
		if (qu1)				{ tempGridRows[nextpos] = qu1[mpistartrow - 1];					qu1[mpistartrow - 1]				= qu1[0];					nextpos++; }
		if (qu1_ext)			{ tempGridRows[nextpos] = qu1_ext[mpistartrow - 1];				qu1_ext[mpistartrow - 1]			= qu1_ext[0];				nextpos++; }
		if (upperIndexSnow)		{ tempGridRows[nextpos] = upperIndexSnow[mpistartrow - 1];		upperIndexSnow[mpistartrow - 1]		= upperIndexSnow[0];		nextpos++; }
		if (zn)					{ tempGridRows[nextpos] = zn[mpistartrow - 1];					zn[mpistartrow - 1]					= zn[0];					nextpos++; }
		
		if (ETP[0])				{ tempGridRows[nextpos] = ETP[0][mpistartrow - 1];				ETP[0][mpistartrow - 1]				= ETP[0][0];				nextpos++; }
		if (ETR[0])				{ tempGridRows[nextpos] = ETR[0][mpistartrow - 1];				ETR[0][mpistartrow - 1]				= ETR[0][0];				nextpos++; }

		for (int z = 1; z <= iNumVegetationLayers; z++)
		{
			if (LUSE[z]) { tempGridRows[nextpos] = LUSE[z][mpistartrow - 1]; LUSE[z][mpistartrow - 1]= LUSE[z][0]; nextpos++; }
			if (ETP[z])  { tempGridRows[nextpos] = ETP[z][mpistartrow - 1] ; ETP[z][mpistartrow - 1] = ETP[z][0];  nextpos++; }
			if (ETR[z])  { tempGridRows[nextpos] = ETR[z][mpistartrow - 1] ; ETR[z][mpistartrow - 1] = ETR[z][0];  nextpos++; }
			if (VBG[z])  { tempGridRows[nextpos] = VBG[z][mpistartrow - 1] ; VBG[z][mpistartrow - 1] = VBG[z][0];  nextpos++; }
			if (WURZ[z]) { tempGridRows[nextpos] = WURZ[z][mpistartrow - 1]; WURZ[z][mpistartrow - 1]= WURZ[z][0]; nextpos++; }
			if (LAI[z])  { tempGridRows[nextpos] = LAI[z][mpistartrow - 1] ; LAI[z][mpistartrow - 1] = LAI[z][0];  nextpos++; }
			if (RSC[z])  { tempGridRows[nextpos] = RSC[z][mpistartrow - 1] ; RSC[z][mpistartrow - 1] = RSC[z][0];  nextpos++; }
			if (RSI[z])  { tempGridRows[nextpos] = RSI[z][mpistartrow - 1] ; RSI[z][mpistartrow - 1] = RSI[z][0];  nextpos++; }
			if (Z0[z])   { tempGridRows[nextpos] = Z0[z][mpistartrow - 1]  ; Z0[z][mpistartrow - 1]  = Z0[z][0];   nextpos++; }
		}
		nextpos = 0;
		// stacks for soil
		if (theta)  { tempStackRows[nextpos] = theta[mpistartrow - 1]	;	theta[mpistartrow - 1]	= theta[0]	;	nextpos++; }
		if (hh)		{ tempStackRows[nextpos] = hh[mpistartrow - 1]		;	hh[mpistartrow - 1]		= hh[0]		;	nextpos++; }
		if (hgeo)	{ tempStackRows[nextpos] = hgeo[mpistartrow - 1]	;	hgeo[mpistartrow - 1]	= hgeo[0]	;	nextpos++; }
		if (qu)		{ tempStackRows[nextpos] = qu[mpistartrow - 1]		;	qu[mpistartrow - 1]		= qu[0]		;	nextpos++; }
		if (temper) { tempStackRows[nextpos] = temper[mpistartrow - 1]	;	temper[mpistartrow - 1] = temper[0]	;	nextpos++; }
		if (SEstack) { tempStackRows[nextpos] = SEstack[mpistartrow - 1];	SEstack[mpistartrow - 1] = SEstack[0];	nextpos++; }

		// stacks for snow
		if (thetaSnow)  { tempStackRows[nextpos] = thetaSnow[mpistartrow - 1]	;	thetaSnow[mpistartrow - 1]	= thetaSnow[0]	;nextpos++; }
		if (hhydSnow)	{ tempStackRows[nextpos] = hhydSnow[mpistartrow - 1]	;	hhydSnow[mpistartrow - 1]	= hhydSnow[0]	;nextpos++; }
		if (hgeoSnow)	{ tempStackRows[nextpos] = hgeoSnow[mpistartrow - 1]	;	hgeoSnow[mpistartrow - 1]	= hgeoSnow[0]	;nextpos++; }
		if (quSnow)		{ tempStackRows[nextpos] = quSnow[mpistartrow - 1]		;	quSnow[mpistartrow - 1]		= quSnow[0]		;nextpos++; }
		if (temperSnow)	{ tempStackRows[nextpos] = temperSnow[mpistartrow - 1]	;	temperSnow[mpistartrow - 1]	= temperSnow[0]	;nextpos++; }
		if (swes)		{ tempStackRows[nextpos] = swes[mpistartrow - 1]		;	swes[mpistartrow - 1]		= swes[0]		;nextpos++; }
		if (densSnow)	{ tempStackRows[nextpos] = densSnow[mpistartrow - 1]	;	densSnow[mpistartrow - 1]	= densSnow[0]	;nextpos++; }
		if (thickSnow)	{ tempStackRows[nextpos] = thickSnow[mpistartrow - 1]	;	thickSnow[mpistartrow - 1]	= thickSnow[0]	;nextpos++; }
		if (swel)		{ tempStackRows[nextpos] = swel[mpistartrow - 1]		;	swel[mpistartrow - 1]		= swel[0]		;nextpos++; }
		
		if (ExcessIceFactor) { tempStackRows[nextpos] = ExcessIceFactor[mpistartrow - 1];	ExcessIceFactor[mpistartrow - 1] = ExcessIceFactor[0]; nextpos++; }
		// todo: tracer grids and stacks

		// adapt FirstLastCol now temporary
		tempFirstLastCol.First = FirstLastCol[mpistartrow - 1].First; 
		FirstLastCol[mpistartrow - 1].First = FirstLastCol[0].First;
		tempFirstLastCol.Last = FirstLastCol[mpistartrow - 1].Last;
		FirstLastCol[mpistartrow - 1].Last = FirstLastCol[0].Last;

	}
}

void unsatzonmodellklasse::RecoverOrigRowPtrs()
{
	if (tempGridRows)
	{
		int nextpos = 0;
		if (gridinfo[2].pointer) { gridinfo[2].pointer[mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
		if (ALBEDO)				{ ALBEDO[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (CLAY)				{ CLAY[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (DELTAGW)			{ DELTAGW[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (DRAINDIST)			{ DRAINDIST[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (DRAINLEVEL)			{ DRAINLEVEL[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (EBAL)				{ EBAL[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (ETR2ETP)			{ ETR2ETP[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (EVAP)				{ EVAP[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (EVAR)				{ EVAR[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (FDD)				{ FDD[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (FROZENDEPTH)		{ FROZENDEPTH[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (FZT)				{ FZT[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (GW1COLMATION)		{ GW1COLMATION[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (GW1COLMFROZEN)		{ GW1COLMFROZEN[mpistartrow - 1]	= tempGridRows[nextpos]; nextpos++; }
		if (GW1COLMUNFROZEN)	{ GW1COLMUNFROZEN[mpistartrow - 1]	= tempGridRows[nextpos]; nextpos++; }
		if (GWEX)				{ GWEX[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (GWIN)				{ GWIN[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (GWLEVEL)			{ GWLEVEL[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (GWN)				{ GWN[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (GWSTAND)			{ GWSTAND[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (GWtheta)			{ GWtheta[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (HRU)				{ HRU[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (HRUMaster)			{ HRUMaster[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (ICE)				{ ICE[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (ICEOLD)				{ ICEOLD[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (INFEX)				{ INFEX[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (INFILT)				{ INFILT[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (IRRIG)				{ IRRIG[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (IRRIGCODE)			{ IRRIGCODE[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (KOLM)				{ KOLM[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (LAKE)				{ LAKE[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (MACRO)				{ MACRO[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (MAXPOND)			{ MAXPOND[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (MAXPUMP)			{ MAXPUMP[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (MPIREG)				{ MPIREG[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (P)					{ P[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (PERCOL)				{ PERCOL[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (POND)				{ POND[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (P_RAIN)				{ P_RAIN[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (QDRAIN)				{ QDRAIN[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (QDSU)				{ QDSU[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (QDgrid)				{ QDgrid[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (QIgrid)				{ QIgrid[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (QLOWBND)			{ QLOWBND[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (QSNOW_OUT)			{ QSNOW_OUT[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (RIVDEP)				{ RIVDEP[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (RIVLNK)				{ RIVLNK[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (RIVWIT)				{ RIVWIT[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (ROOT)				{ ROOT[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (RSE)				{ RSE[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (RelCellArea)		{ RelCellArea[mpistartrow - 1]		= tempGridRows[nextpos]; nextpos++; }
		if (S0)					{ S0[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SATT)				{ SATT[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SBHIGH)				{ SBHIGH[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (SB_1_)				{ SB_1_[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (SB_2_)				{ SB_2_[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (SLIQ)				{ SLIQ[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SLP)				{ SLP[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SNOW)				{ SNOW[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SOIL)				{ SOIL[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SSNO)				{ SSNO[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SSTO)				{ SSTO[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (SUBSTEPS)			{ SUBSTEPS[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (SURFTEMP)			{ SURFTEMP[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (T)					{ T[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (TANBETA)			{ TANBETA[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (TDAY)				{ TDAY[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (TDD)				{ TDD[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (THAWTMOD)			{ THAWTMOD[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (THROOT)				{ THROOT[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (TLOWBND)			{ TLOWBND[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (TYEAR)				{ TYEAR[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (UPRISE)				{ UPRISE[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
//		if (gridinfo[2].pointer){ gridinfo[2].pointer[mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
		if (gw_thick)			{ gw_thick[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (qu1)				{ qu1[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }
		if (qu1_ext)			{ qu1_ext[mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (upperIndexSnow)		{ upperIndexSnow[mpistartrow - 1]	= tempGridRows[nextpos]; nextpos++; }
		if (zn)					{ zn[mpistartrow - 1]				= tempGridRows[nextpos]; nextpos++; }

		if (ETP[0])				{ ETP[0][mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }
		if (ETR[0])				{ ETR[0][mpistartrow - 1]			= tempGridRows[nextpos]; nextpos++; }

		for (int z = 1; z <= iNumVegetationLayers; z++)
		{
			if (LUSE[z]) { LUSE[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (ETP[z])  { ETP[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (ETR[z])  { ETR[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (VBG[z])  { VBG[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (WURZ[z]) { WURZ[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (LAI[z])  { LAI[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (RSC[z])  { RSC[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (RSI[z])  { RSI[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
			if (Z0[z])   { Z0[z][mpistartrow - 1] = tempGridRows[nextpos]; nextpos++; }
		}
		nextpos = 0;
		// stacks for soil
		if (theta)  { theta[mpistartrow - 1]	= tempStackRows[nextpos];	nextpos++; }
		if (hh)		{ hh[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (hgeo)	{ hgeo[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (qu)		{ qu[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (temper) { temper[mpistartrow - 1]   = tempStackRows[nextpos];	nextpos++; }
		if (SEstack) { SEstack[mpistartrow - 1] = tempStackRows[nextpos];	nextpos++; }

		// stacks for snow
		if (thetaSnow)  { thetaSnow[mpistartrow - 1]	= tempStackRows[nextpos];	nextpos++; }
		if (hhydSnow)	{ hhydSnow[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (hgeoSnow)	{ hgeoSnow[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (quSnow)		{ quSnow[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (temperSnow)	{ temperSnow[mpistartrow - 1]	= tempStackRows[nextpos];	nextpos++; }
		if (swes)		{ swes[mpistartrow - 1]			= tempStackRows[nextpos];	nextpos++; }
		if (densSnow)	{ densSnow[mpistartrow - 1]		= tempStackRows[nextpos];	nextpos++; }
		if (thickSnow)	{ thickSnow[mpistartrow - 1]	= tempStackRows[nextpos];	nextpos++; }
		if (swel)		{ swel[mpistartrow - 1]			= tempStackRows[nextpos];	nextpos++; }
		
		if (ExcessIceFactor) { ExcessIceFactor[mpistartrow - 1] = tempStackRows[nextpos];	nextpos++; }

		// todo: tracer grids and stacks

		// adapt FirstLastCol now temporary
		FirstLastCol[mpistartrow - 1].First = tempFirstLastCol.First;
		FirstLastCol[mpistartrow - 1].Last = tempFirstLastCol.Last;

	}
}


// sends the river flow table (which in turn was filled by runoff routing) from rank 0 to all ranks in order to provide the actual flow in all link channels
void unsatzonmodellklasse::MPIBCastRiverFlowTable()
{
#ifdef _MPI

	if (riverflow_table)
	{
		//	struct riverflow_table_type
//	{
//		double flow_mm[10];      /* <- from routing  : the flow, entering the actual routing channel from upstream in mm */
//		double flowm3s;          /* <- from routing  : same variable but in m3/s */
//		double speed;            /* <- from routing  : velocity of the flow, to get depth: depth =  */
//		double upstream_area;    /* <- from routing  : total sum of the upstream tributaries */
//		double infilt_amount;    /* -> from unsatzon : total infitration amount in m3 */
//		double infilt_fraction;  /* <-> from both    : fraction of infiltrating water on river flow in a channel */
//		int    ROUT_index;       /* <- from routing  : in wich routing description is this subbasin included as target subbasin (e.g. as x in "TG x (..) from ..." */
//		int    OL_index;         /* <- from routing  : what OL-index has this subbasin within the routing description, important for SUMTRIBs, since they have in around 50% of all cases no own channel */
//	};

		const int nitems = 8;
		int blocklengths[8] = { 10, 1, 1, 1, 1, 1, 1, 1 };
		MPI_Datatype types[8] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT };
		MPI_Aint offsets[8] = { offsetof(riverflow_table_type, flow_mm),
								 offsetof(riverflow_table_type, flowm3s),
								 offsetof(riverflow_table_type, speed),
								 offsetof(riverflow_table_type, upstream_area),
								 offsetof(riverflow_table_type, infilt_amount),
								 offsetof(riverflow_table_type, infilt_fraction),
								 offsetof(riverflow_table_type, ROUT_index),
								 offsetof(riverflow_table_type, OL_index) };

		MPI_Datatype mpi_riverflow_table_type;
		MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_riverflow_table_type);
		MPI_Type_commit(&mpi_riverflow_table_type);

		MPI_Bcast(&(riverflow_table[0]), zonmax + 1, mpi_riverflow_table_type, 0, MPI_COMM_WORLD);
		MPI_Type_free(&mpi_riverflow_table_type);
	}
#endif

}


#ifdef _MPI
void infilt_amount_sum(riverflow_table_type const *in, riverflow_table_type *inout, int *len, const MPI_Datatype *datatype)
{
	for (int i = 0; i<*len; i++)
	{
		inout[i].infilt_amount += in[i].infilt_amount;
	}
}
#endif

// gather the results of infiltrated water in rank 0 in order to allow runoff routing in the next time step to use the reduced outflwo as inflow into other river links
void unsatzonmodellklasse::MPIGatherRiverFlowTable()
{
	// moving the weighted avarage from infilt_fraction to infilt_amount
	if (riverflow_table)
	{
		for (int i = 0; i <= zonmax; i++)
		{
			riverflow_table[i].infilt_amount = riverflow_table[i].infilt_fraction;
			riverflow_table[i].infilt_fraction = 0;
		}

#ifdef _MPI
		const int nitems = 8;
		int blocklengths[8] = { 10, 1, 1, 1, 1, 1, 1, 1 };
		MPI_Datatype types[8] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT };
		MPI_Aint offsets[8] = { offsetof(riverflow_table_type, flow_mm),
			offsetof(riverflow_table_type, flowm3s),
			offsetof(riverflow_table_type, speed),
			offsetof(riverflow_table_type, upstream_area),
			offsetof(riverflow_table_type, infilt_amount),
			offsetof(riverflow_table_type, infilt_fraction),
			offsetof(riverflow_table_type, ROUT_index),
			offsetof(riverflow_table_type, OL_index) };

		MPI_Datatype mpi_riverflow_table_type;
		MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_riverflow_table_type);
		MPI_Type_commit(&mpi_riverflow_table_type);

		void(*infilt_amount_sum_function)(riverflow_table_type const *in, riverflow_table_type *inout, int *len, const MPI_Datatype *datatype) = &infilt_amount_sum;

		MPI_Op mySumOp;
		MPI_Op_create((MPI_User_function*)infilt_amount_sum_function, true, &mySumOp);

		MPI_Barrier(MPI_COMM_WORLD);

		if (mpi_rank == 0)
		{
			MPI_Reduce(MPI_IN_PLACE, riverflow_table, zonmax + 1, mpi_riverflow_table_type, mySumOp, 0,MPI_COMM_WORLD);
		}
		else
		{
			MPI_Reduce(riverflow_table, NULL, zonmax + 1, mpi_riverflow_table_type, mySumOp, 0, MPI_COMM_WORLD);
		}


		// #ifdef _MPIDEBUG
		MPI_Barrier(MPI_COMM_WORLD);
		// #endif

		MPI_Op_free(&mySumOp);
		MPI_Type_free(&mpi_riverflow_table_type);
#endif
	}
}


void unsatzonmodellklasse::MPIPrepareRiverFlowTable()
{
	// since unsatzonmodell::run can be called multiple times, but the riverflow_table is taken from routing model only once (and will be sent only once, as well), the 
	// re-infiltrated amount of water for each substep needs to be stored in another field of that record, 
	// so each sub time step can work with an empty infilt_amount again (but flowm3s will be reduced in each sub step automatically according to the number of sub steps)
	if (riverflow_table)
	{
		for (int i = 0; i <= zonmax; i++)
		{
			riverflow_table[i].infilt_fraction += riverflow_table[i].infilt_amount / substepcnt;
			riverflow_table[i].infilt_amount = 0;
		}
	}
}
