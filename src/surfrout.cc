/*   $Revision: 30 $ 
	$Date: 8.06.10 0:05 $  */

/*
(c) 2007 by Joerg Schulla

*/

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <cmath>
#include <limits>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "surfrout.h"
#include "wasim.h"
#ifdef _MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


surfaceroutingclass::surfaceroutingclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults, int maxLastResults,
		double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int p_numsubb,int *p_listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
		int p_zonmax,int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll,double p_yll,
		double p_cellsize,double p_nodata,
		long mine,long ste , long tge , long mte , long jre, double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU,int substepcnt):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
		p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
		p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr,
		p_LUtable,p_SOILtable,p_MultiLayerLUTable,
		WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
		p_doSpinUp, p_routingmode,p_vHRU,substepcnt),
// initializing member variables here to 0/NULL or whatever...
method(0),storagefile(0), storagefilename(0),
DHM(0),T(0),LUSE(0),VCF(0),SOIL(0),SLOPE(0),SSSU(0),FLD(0),QDgrid(0),QIgrid(0),CONC_QD(0),CONC_QI(0),CONC_SSSU(0),RIVDEP(0),RIVWIT(0),FZT(0),
QB(0),QI(0),QIi(0),QD(0),QDd(0),Qg(0),
SRFCVEL(0),SRFCVEL1(0),
FLUX(0), FLUX1(0), WEIGHT(0), WEIGHT1(0), HEAD(0), HEADold(0), ROUGHX(0), ROUGHY(0), CellSizeX(0),CellSizeY(0),PLANCURV(0), SRQIN(0), SRQOUT(0), TAUCRIT(0),
SRFCVELwritegrid(0),FLUXwritegrid(0),
SRFCVELgridname(0),FLUXgridname(0),
srfcvelindex(0), fluxindex(0), 
headindex(0), headoldindex(0), roughxindex(0), roughyindex(0), manningindex(0),
plancurvindex(0), weightindex(0), srqinindex(0), srqoutindex(0), taucritindex(0),
QDSRstat(0),zoneQDSR(0),QDSRstatcode(0),QDSRstatfilename(0),
QISRstat(0),zoneQISR(0),QISRstatcode(0),QISRstatfilename(0),
QBSRstat(0),zoneQBSR(0),QBSRstatcode(0),QBSRstatfilename(0),
QGSRstat(0),zoneQGSR(0),QGSRstatcode(0),QGSRstatfilename(0),
SFSTSRstat(0),zoneSFSTSR(0),SFSTSRstatcode(0),SFSTSRstatfilename(0),
iNumVegetationLayers(0),
iMaxIterANL(0),iMaxIterVelocity(0),
dMaxDeltaANL(0),dMaxDeltaVelocity(0),
dMinSubTimeStep(0),dMaxSubTimeStep(0),
zoneFZT_QD(0),zoneFZT_QI(0),
dConcFact(1.0),dMinFlowDepth(0),
old_year(0), old_month(0), old_day(0), old_hour(0), old_minute(0),moduleSubstepCounter(0),
zoneHEADStart(0), zoneHEADEnd(0), zoneQDFromFlux(0), QGScalings(0), numCellsQD(0),
dumpmode(false), maxDumpFiles(0), actDumpFiles(0), closeAfterDump(false), repairAfterDump(false)
{
	try
	{
		retvalue=0;
		int i,j,k;

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

		/* ACHTUNG !!!!!: die Steuerdatei wird im Constructor der MutterKlasse (teilmodell)
	geoeffnet und bis zum Erreichen des

	ignoranzparameters 
	+ Intervallaenge

	gelesen.
	ist der Ignoranzparameter = 0, dann wird die Datei nach diesem Parameter (InV) nicht
	weitergelesen und gleich geschlossen.
	So kann der gesamte Ablauf der ini-Prozedur abgebrochen werden, da der
	Ignoranzparameter schon bekannt ist!
	*/

		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;

#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
#endif
		if (doSpinUp!=0 && HRU==NULL)
		{
			do_it=0;
			return;
		}

		if (ist_vorhanden(_GlacierizedCells_)||ist_vorhanden(_GlacierCodes_))
		{
			(*WSstream) << "FATAL: error: Surface routing module cannot be used when glacier model is active. Remove grids for glacier codes and glaciated cell fractions from the grid list"<<"[end-entry]";
			throw(2);
		}
		// getting addresses of existing grids (standardgrids or interpolation results etc.)

		// elevation model
		DHM = gridinfo[2].pointer; /* Hoehenmodell von jetzt an unter DHM ansprechbar */

		// Multilayer-Land use or landuse grid
		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}
		
		char name[500];
		int lay;

		LUSE   = new double**[iNumVegetationLayers+1];
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;
		}

		VCF   = new double**[iNumVegetationLayers+1];
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Vegetationsbedeckungsgrad_,lay);
			VCF[lay]   = gridinfo[getindex(name)].pointer;
		}

		// soil types
		SOIL  =  gridinfo[getindex(_Bodenarten_)].pointer; // soil types (texture)

		// searching surface storage of silting up module which is used in surface routing, too
		// to be available in surface routing AND silting up module, this grid should be read in 
		// as a variable grid
		SSSU  =  gridinfo[getindex(_SURFSTORSILTINGUP_)].pointer;

		// reserving memeory for the filenames (will be released in this->FreeMemeory like all other allocated memory)
		QDSRstatfilename	= newstring(_stringsize_+1,1);
		QISRstatfilename	= newstring(_stringsize_+1,1);
		QBSRstatfilename	= newstring(_stringsize_+1,1);
		QGSRstatfilename	= newstring(_stringsize_+1,1);
		SFSTSRstatfilename	= newstring(_stringsize_+1,1);
		SRFCVELgridname     = newstring(_stringsize_+1,1);
		FLUXgridname        = newstring(_stringsize_+1,1);

		storagefilename = (char*)calloc(_stringsize_, sizeof(char));

		/* lesen der Werte aus der Steuerdatei */
		fscanfln(steuerdatei,"%d",&method,"method for flowpath-calculation ",tstr);
		if ((method<1)||(method>6)) {
			(*WSstream) << "FATAL: \n-----------------error------------\n"
			"invalid method ("<<method<<") read for "<<instanzname<<"\n"<<"[end-entry]";
			throw(10);
		}

		fscanfln(steuerdatei,"%s%d",QDSRstatfilename,&QDSRstatcode ,"name of the statistics for direct runoff from surface routing (QDSR) ",tstr);
		fscanfln(steuerdatei,"%s%d",QISRstatfilename,&QISRstatcode ,"name of the statistics for interflow from surface routing (QISR) ",tstr);
		fscanfln(steuerdatei,"%s%d",QBSRstatfilename,&QBSRstatcode ,"name of the statistics for baseflow from surface routing (QBSR) ",tstr);
		fscanfln(steuerdatei,"%s%d",QGSRstatfilename,&QGSRstatcode ,"name of the statistics for total runoff from surface routing (QGSR) ",tstr);
		fscanfln(steuerdatei,"%s",SRFCVELgridname,"name of the surface-routing-velocity-stack, "    ,tstr);
		fscanfln(steuerdatei,"%d",&SRFCVELwritegrid,"Writegrid for surface-routing-velocity-stack, ",tstr);
		fscanfln(steuerdatei,"%s",FLUXgridname,"name of the surface-routing-fluxes-stack, "    ,tstr);
		fscanfln(steuerdatei,"%d",&FLUXwritegrid,"Writegrid for surface-routing-fluxes-stack, ",tstr);

		// reading other parameters

		fscanfln(steuerdatei,"%lf",&dMaxDeltaANL,"maximum iteration error for A_NL [m] in ",tstr);
		fscanfln(steuerdatei,"%d",&iMaxIterANL,"maximum iteration steps for A_NL in ",tstr);
		fscanfln(steuerdatei,"%lf",&dMaxDeltaVelocity,"maximum iteration error [m/s] for velocity in ",tstr);
		fscanfln(steuerdatei,"%d",&iMaxIterVelocity,"maximum iteration steps for velocity in ",tstr);
		fscanfln(steuerdatei,"%lf",&dMinSubTimeStep,"minimum sub time step [s] in ",tstr);
		fscanfln(steuerdatei,"%lf",&dMaxSubTimeStep,"maximum sub time step [s] in ",tstr);
		fscanfln(steuerdatei,"%lf",&dMinFlowDepth,"minimum flow depth to take into account the roughness of the crops [m] in ",tstr);
		fscanfln(steuerdatei,"%lf",&dConcFact,"Concentration factor [1/1] in ",tstr);

		// making sure that InV/dMaxSubTimeStep is integer
		double temp=dMaxSubTimeStep;
		if ((int)((InV/substepcnt*60)/dMaxSubTimeStep)>0)
			dMaxSubTimeStep = (InV/substepcnt*60.0)/((int)((InV/substepcnt*60)/dMaxSubTimeStep));
		else
			dMaxSubTimeStep=(InV/substepcnt*60.0);
		if (temp!=dMaxSubTimeStep) 
		{
			(*WSstream) << "INFO: longest allowed subtimestep changed from "<<temp<<" to "<<dMaxSubTimeStep<<". Must be integer fraction of time step in seconds."<<"[end-entry]";
		}

		// making sure that InV/dMinSubTimeStep is integer
		temp=dMinSubTimeStep;
		if ((int)((InV/substepcnt*60)/dMinSubTimeStep)>0)
			dMinSubTimeStep = (InV/substepcnt*60.0)/((int)((InV/substepcnt*60)/dMinSubTimeStep));
		else
			dMinSubTimeStep=(InV/substepcnt*60.0);
		if (temp!=dMinSubTimeStep) 
		{
			(*WSstream) << "WARNING: shortest allowed subtimestep changed from "<<temp<<" to "<<dMinSubTimeStep<<". Must be integer fraction of time step in seconds."<<"[end-entry]";
		}

		// reading the value of readgrids from control file
		if (fscanfln(steuerdatei,"%d",&gridslesen   ,"read control parameter for storage grids for ",tstr)!=1)
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
			"error reading readgrid-Parameter for "<<instanzname<<"\n"
			"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
			throw(10);
		}

		// reading the optional name of the thaw depth statistics grid
		if (fscanf(steuerdatei,"%s%d",SFSTSRstatfilename,&SFSTSRstatcode)!=2)
		{
			SFSTSRstatcode=0;
			SFSTSRstatfilename=NULL;
		}

//		// file name of the storage file containing all content of the single linear storages for all flow time zones 
//		if ((fscanfln(steuerdatei, "%s", storagefilename, "name of the storage content file ", tstr) != 1) || (strcmp(storagefilename,"[unsatzon_model]")==0))
//		{
//			(*WSstream) << "WARNING: \n-----------------error------------\n"
//				"error reading storage file name for " << instanzname << "\n"
//				"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block" << "[end-entry]";
//			// throw(10);
//		}
//		char tsn[1000];
//		sprintf(tsn, "%s_%d", storagefilename, mpi_rank);
//		sprintf(storagefilename, "%s", tsn);


		// now we can close the control file, since the last parameter was read successfully
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing control file \""<<steuerdateiname<<"\""<<"[end-entry]";
			throw(11);
		}

		// looking for the first free entry in the gridinfo array
		int sti=1;
		while (gridinfo[sti].identification!=NULL)
		{
			sti ++;
			if (sti==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid anymore, allowed number of grids is "<<gridmaxanzahl<<"!"<<"[end-entry]";
				throw(12);
			}
		}

		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;

		// create grids and stacks
		SRFCVEL = drei_d_matrix(sr,er,1,ncols,10,nrows);
		SRFCVEL1 = matrix(sr, er, 1, ncols,nrows);
		initialize_grid(sti, SRFCVELgridname, SRFCVEL1, SRFCVEL, 10, _SURFACEROUTING_VELOCITY_, SRFCVELwritegrid, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j=1;j<=ncols;j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j]==mpi_rank))
				{
					for (int lv=0;lv<=9;lv++)
					{
						SRFCVEL[k][j][lv] = 0;
					}
				}
				else
				{
					for (int lv=0;lv<=9;lv++)
					{
						SRFCVEL[k][j][lv] = nodata;
					}
				}
			}
		}
		srfcvelindex = sti;
		sti++;

		// create grids and stacks
		FLUX = drei_d_matrix(sr,er,1,ncols,10,nrows);
		FLUX1 = matrix(sr,er,1,ncols,nrows);
		initialize_grid(sti, FLUXgridname, FLUX1, FLUX, 10, _SURFACEROUTING_FLUXES_, FLUXwritegrid, sr,er);
		for (int k=mpistartrow;k<=mpiendrow;k++) 
		{
			for (int j=1;j<=ncols;j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j]==mpi_rank))
				{
					for (int lv=0;lv<=9;lv++)
					{
						FLUX[k][j][lv] = 0;
					}
				}
				else
				{
					for (int lv=0;lv<=9;lv++)
					{
						FLUX[k][j][lv] = nodata;
					}
				}
			}
		}
		fluxindex = sti;
		MPISynchronizeStackHalo(fluxindex, 10, FLUX);
		sti++;

		WEIGHT = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, 10,nrows);
		WEIGHT1 = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(sti, "", WEIGHT1, WEIGHT, 10, _SURFACEROUTING_WEIGHTS_, 0, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j=1;j<=ncols;j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					for (int lv=0;lv<8;lv++)
					{
						WEIGHT[k][j][lv] = 0;
					}
				}
				else
				{
					for (int lv=0;lv<8;lv++)
					{
						WEIGHT[k][j][lv] = nodata;
					}
				}
			}
		}
		weightindex = sti;
		sti++;

		HEAD = matrix(sr, er, 1, ncols, nrows);
		initialize_grid(sti, "", HEAD, NULL, 0, _SURFACEROUTING_HEAD_, 0, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			//for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					HEAD[k][j] = 0;
				}
				else
				{
					HEAD[k][j] = nodata;
				}
			}
		}
		headindex = sti;
		MPISynchronizeGridHalo(headindex, HEAD);
		//MPISynchronizeGridHaloJustRecv(headindex, HEAD); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex, HEAD); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex, HEAD); //LAURA
		sti++;

		HEADold = matrix(sr, er, 1, ncols, nrows);
		initialize_grid(sti, "", HEADold, NULL, 0, _SURFACEROUTING_HEAD_OLD_, 0, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			//if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					HEADold[k][j] = 0;
				}
				else
				{
					HEADold[k][j] = nodata;
				}
			}
		}
		headoldindex = sti;
		MPISynchronizeGridHalo(headoldindex, HEADold);
		//MPISynchronizeGridHaloJustRecv(headoldindex, HEADold); //LAURA
		//MPISynchronizeGridHaloJustSend(headoldindex, HEADold); //LAURA
		//MPISynchronizeGridHaloJustWait(headoldindex, HEADold); //LAURA
		sti++;

		ROUGHX = matrix(sr, er, 1, ncols, nrows);
		initialize_grid(sti, "", ROUGHX, NULL, 0, _SURFACEROUTING_ROUGHNESS_X_, 0, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			//if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					ROUGHX[k][j] = 10;
				}
				else
				{
					ROUGHX[k][j] = 0;
				}
			}
		}
		roughxindex = sti;
		MPISynchronizeGridHalo(roughxindex, ROUGHX);
		//MPISynchronizeGridHaloJustRecv(roughxindex, ROUGHX); //LAURA
		//MPISynchronizeGridHaloJustSend(roughxindex, ROUGHX); //LAURA
		//MPISynchronizeGridHaloJustWait(roughxindex, ROUGHX); //LAURA
		sti++;

		ROUGHY = matrix(sr, er, 1, ncols, nrows);
		initialize_grid(sti, "", ROUGHY, NULL, 0, _SURFACEROUTING_ROUGHNESS_Y_, 0, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			//if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					ROUGHY[k][j] = 10;
				}
				else
				{
					ROUGHY[k][j] = 0;
				}
			}
		}
		roughyindex = sti;
		MPISynchronizeGridHalo(roughyindex, ROUGHY);
		//MPISynchronizeGridHaloJustRecv(roughyindex, ROUGHY); //LAURA
		//MPISynchronizeGridHaloJustSend(roughyindex, ROUGHY); //LAURA
		//MPISynchronizeGridHaloJustWait(roughyindex, ROUGHY); //LAURA
		sti++;

		MANNING = matrix(sr, er, 1, ncols, nrows);
		initialize_grid(sti, "", MANNING, NULL, 0, _SURFACEROUTING_MANNING_COEFF_, 0, sr, er);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = 1; j <= ncols; j++)
			//if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					int luse = -9999;
					int reduction = 0;
					while ((luse == -9999) && (reduction<iNumVegetationLayers))
					{
						luse = LUSE[iNumVegetationLayers - reduction][k][j];
						reduction++;
					}
					MANNING[k][j] = this->LUtable[luse].GetSR_K_ST();
				}
				else
				{
					MANNING[k][j] = 0;
				}
			}
		}
		manningindex = sti;
		MPISynchronizeGridHalo(manningindex, MANNING);
		//MPISynchronizeGridHaloJustRecv(manningindex, MANNING); //LAURA
		//MPISynchronizeGridHaloJustSend(manningindex, MANNING); //LAURA
		//MPISynchronizeGridHaloJustWait(manningindex, MANNING); //LAURA
		sti++;

		if (ist_vorhanden(_PLANCURVATURE_))
		{
			plancurvindex=getindex(_PLANCURVATURE_);
			PLANCURV=gridinfo[plancurvindex].pointer;
			plancurvindex=-plancurvindex;
		}
		else
		{
			PLANCURV = matrix(sr,er,1,ncols,nrows);
			initialize_grid(sti, "", PLANCURV, NULL, 0, _PLANCURVATURE_, 0, sr,er);
			GetPlanCurvature();
			plancurvindex = sti;
			MPISynchronizeGridHalo(plancurvindex, PLANCURV);
			sti++;
		}

		if (ist_vorhanden(_SURFACEROUTING_QIN_))
		{
			srqinindex=getindex(_SURFACEROUTING_QIN_);
			SRQIN=gridinfo[srqinindex].pointer;
			srqinindex=-srqinindex;
		}
		else
		{
			SRQIN = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(sti, "", SRQIN, NULL, 0, _SURFACEROUTING_QIN_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j=1;j<=ncols;j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						SRQIN[k][j] = 0;
					} else	
					{
						SRQIN[k][j] = nodata;
					}
				}
			}
			srqinindex = sti;
			sti++;
		}

		if (ist_vorhanden(_SURFACEROUTING_QOUT_))
		{
			srqoutindex=getindex(_SURFACEROUTING_QOUT_);
			SRQOUT=gridinfo[srqoutindex].pointer;
			srqoutindex=-srqoutindex;
		}
		else
		{
			SRQOUT = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(sti, "", SRQOUT, NULL, 0, _SURFACEROUTING_QOUT_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j=1;j<=ncols;j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						SRQOUT[k][j] = 0;
					} else	
					{
						SRQOUT[k][j] = nodata;
					}
				}
			}
			srqoutindex = sti;
			sti++;
		}

		if (ist_vorhanden(_CRITICAL_SHEAR_STRESS_))
		{
			taucritindex=getindex(_CRITICAL_SHEAR_STRESS_);
			TAUCRIT=gridinfo[taucritindex].pointer;
			taucritindex=-taucritindex;
		}
		else
		{
			TAUCRIT = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
			initialize_grid(sti, "", TAUCRIT, NULL, 0, _CRITICAL_SHEAR_STRESS_, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j=1;j<=ncols;j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						TAUCRIT[k][j] = 0.9;
					} else	
					{
						TAUCRIT[k][j] = nodata;
					}
				}
			}
			taucritindex = sti;
			sti++;
		}

		// initializing the statistics registers for all statistics	tables we use in this module
		// statistics registers for direct discharge from surface routing module
		zoneQDSR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQDSR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneQDSR[i].n           = zoneQDSR[i].sx        = 0;
			zoneQDSR[i].zonenwert   = zone[i].zonenwert;
			zoneQDSR[i].rel         = zone[i].rel;
			zoneQDSR[i].RelCellSum  = zone[i].RelCellSum;
		}

		// statistics registers for interflow from surface routing module
		zoneQISR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQISR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneQISR[i].n           = zoneQISR[i].sx        = 0;
			zoneQISR[i].zonenwert   = zone[i].zonenwert;
			zoneQISR[i].rel         = zone[i].rel;
			zoneQISR[i].RelCellSum  = zone[i].RelCellSum;
		}

		// statistics registers for baseflow from surface routing module
		zoneQBSR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQBSR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneQBSR[i].n           = zoneQBSR[i].sx        = 0;
			zoneQBSR[i].zonenwert   = zone[i].zonenwert;
			zoneQBSR[i].rel         = zone[i].rel;
			zoneQBSR[i].RelCellSum  = zone[i].RelCellSum;
		}

		// statistics registers for direct discharge from surface routing module
		zoneQGSR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQGSR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneQGSR[i].n           = zoneQGSR[i].sx        = 0;
			zoneQGSR[i].zonenwert   = zone[i].zonenwert;
			zoneQGSR[i].rel         = zone[i].rel;
			zoneQGSR[i].RelCellSum  = zone[i].RelCellSum;
		}

		// statistics registers for surface storage in surface routing module
		zoneSFSTSR = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSFSTSR) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneSFSTSR[i].n           = zoneSFSTSR[i].sx        = 0;
			zoneSFSTSR[i].zonenwert   = zone[i].zonenwert;
			zoneSFSTSR[i].rel         = zone[i].rel;
			zoneSFSTSR[i].RelCellSum  = zone[i].RelCellSum;
		}

		// initializing the statistics class for each register set we use in this module
		// direct discharge from surface routing
		QDSRstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneQDSR,(int)zonmax,(int)QDSRstatcode,(int)5,(char*)QDSRstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		QDSRstat->schreibzonenwerte("Direct_Discharge_Surface_Routing",instanzname,LastResults,maxLastResults);

		// INTERFLOW from surface routing
		QISRstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneQISR,(int)zonmax,(int)QISRstatcode,(int)5,(char*)QISRstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		QISRstat->schreibzonenwerte("Interflow_Surface_Routing",instanzname,LastResults,maxLastResults);

		// BASEFLOW from surface routing
		QBSRstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneQBSR,(int)zonmax,(int)QBSRstatcode,(int)5,(char*)QBSRstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		QBSRstat->schreibzonenwerte("Baseflow_Surface_Routing",instanzname,LastResults,maxLastResults);

		// total runoff from surface routing
		QGSRstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneQGSR,(int)zonmax,(int)QGSRstatcode,(int)6,(char*)QGSRstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		QGSRstat->schreibzonenwerte("Total_Runoff_Surface_Routing",instanzname,LastResults,maxLastResults);

		// surface storage from surface routing
		SFSTSRstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneSFSTSR,(int)zonmax,(int)SFSTSRstatcode,(int)6,(char*)SFSTSRstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		SFSTSRstat->schreibzonenwerte("Surface_Storage_in_SurfaceRouting_Module",instanzname,LastResults,maxLastResults);

		// the following sequence is for initializing the sample code which works identically to the unsatzonmodel (and to be raplaced by specfic surface routing code)

		// getting pointers to other grids which must exist already
		//	T      = gridinfo[getindex(_Temperatur_)].pointer;
		//  SLOPE  = gridinfo[getindex(_Gefaelle_)].pointer;
		//	FLD    = gridinfo[getindex(_FLOWDIRECTION_)].pointer;
		QDgrid = gridinfo[getindex(_Direktabfluss_)].pointer;
		QIgrid = gridinfo[getindex(_Interflow_)].pointer;
		RIVDEP = gridinfo[getindex(_river_depth_)].pointer;
		RIVWIT = gridinfo[getindex(_river_width_)].pointer;
		RIVLNK = gridinfo[getindex(_river_links_)].pointer;
		if (ist_vorhanden(_Fliesszeiten_))
		{
			FZT = gridinfo[getindex(_Fliesszeiten_)].pointer;
		}
		else
		{
			FZT = NULL;
		}
		if (RelCellArea != NULL)
		{
			CellSizeX = gridinfo[getindex(_CellSizeX_)].pointer;
			CellSizeY = gridinfo[getindex(_CellSizeY_)].pointer;
		}


		// two statistics registers that will be used to scale the QDgrid values when using method 6 since the absolute fluxes are not easy to determine when using an implicite method like Gauss-Seidel
		zoneHEADStart = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
		if (!zoneHEADStart) speichererror(1);
		for (i = 0; i <= zonmax; i++)
		{
			zoneHEADStart[i].n = 0;
			zoneHEADStart[i].sx = 0;
			zoneHEADStart[i].zonenwert = zone[i].zonenwert;
			zoneHEADStart[i].rel = zone[i].rel;
			zoneHEADStart[i].RelCellSum = zone[i].RelCellSum;
		}

		
		zoneHEADEnd = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
		if (!zoneHEADEnd) speichererror(1);
		for (i = 0; i <= zonmax; i++)
		{
			zoneHEADEnd[i].n = 0;
			zoneHEADEnd[i].sx = 0;
			zoneHEADEnd[i].zonenwert = zone[i].zonenwert;
			zoneHEADEnd[i].rel = zone[i].rel;
			zoneHEADEnd[i].RelCellSum = zone[i].RelCellSum;
		}

		zoneQDFromFlux = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
		if (!zoneQDFromFlux) speichererror(1);
		for (i = 0; i <= zonmax; i++)
		{
			zoneQDFromFlux[i].n = 0;
			zoneQDFromFlux[i].sx = 0;
			zoneQDFromFlux[i].zonenwert = zone[i].zonenwert;
			zoneQDFromFlux[i].rel = zone[i].rel;
			zoneQDFromFlux[i].RelCellSum = zone[i].RelCellSum;
		}

		QGScalings = new double[zonmax + 1];
		numCellsQD = new double[zonmax + 1];

		// searching all CONC_QD and CONC_QI-grids from unsatzonmodel -> is be done in tracer_ini-method

		// now delete unneccessary memory
		if (SRFCVELgridname) free(SRFCVELgridname); SRFCVELgridname=NULL;
		if (FLUXgridname) free(FLUXgridname); FLUXgridname=NULL;
		
		initialized=1; 

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

/* a function wich analyses the flow travel time grid (the time, which
surface ruoff from a grid cell takes to reach the next subbasin outlet)
as a result there will allocated a set of statistic registers (for each
subbasin and each flow travel time zone one statistic register).
The distribution of travel zones is printed to the log file */
//#pragma optimize ("p",on)

void surfaceroutingclass::setSubStepCnt(int step)
{
	substepcnt=step;
}

#ifdef _MPI
void zonesumsurfrout(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype)
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

void surfaceroutingclass::fzt_grid_auswerten()
{
	int logdateiopen;
	int i,j,id1,id2,z;
	FILE *logdatei;
	double cellcount;

	/* searching in the flow time grid for the longest flow time in order to reserve enough memory (in steps of time intervals) */
	(*WSstream) << "looking for longest flow travel time" << "[end-entry]";
	maxzeit = (double*)calloc(zonmax + 1, sizeof(double));
	for (int i = 1; i <= zonmax; i++) maxzeit[i] = 0.0;
	if (FZT)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					// was already changed in unsatzon FZT[i][j] = ceil(FZT[i][j] / 60.0 / InV);
					if (FZT[i][j] > maxzeit[zarray[(int)zn[i][j]]])
					{
						maxzeit[zarray[(int)zn[i][j]]] = FZT[i][j];
					}
				}
			}
		}
	}
	else
	{
		for (int i = 1; i <= zonmax; i++)
		{
			maxzeit[i] = 1.0;
		}
	}

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	double *maxzeit2 = new double[zonmax + 1]; // 	double *maxzeit2 = (double*)calloc(zonmax + 1, sizeof(double));
	MPI_Allreduce(&(maxzeit[1]), &(maxzeit2[1]), zonmax, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
	for (int i = 1; i <= zonmax; i++)
	{
		maxzeit[i] = maxzeit2[i];
	}
	delete[] maxzeit2;
	MPI_Barrier(MPI_COMM_WORLD);
	(*WSstream) << "done and reduced to rank 0" << "[end-entry]";
#else
	(*WSstream) << "done" << "[end-entry]";
#endif


	/* per subbasin: allocating an array of statistk-registers as many as specified in "maxzeit[i]" plus one as a buffer,
	because the actual surface flow is routed through these registers by pointer rotating */

	zoneFZT_QD = (statrecord***)calloc(tracer_count + 1, sizeof(statrecord**));
	zoneFZT_QI = (statrecord***)calloc(tracer_count + 1, sizeof(statrecord**));
	for (int z = 0; z <= tracer_count; z++)
	{
		zoneFZT_QD[z] = (statrecord**)calloc(zonmax + 1, sizeof(statrecord*));
		zoneFZT_QI[z] = (statrecord**)calloc(zonmax + 1, sizeof(statrecord*));
		for (int i = 1; i <= zonmax; i++)
		{
			zoneFZT_QD[z][i] = (statrecord*)calloc((int)maxzeit[i] + 1, sizeof(statrecord));
			zoneFZT_QI[z][i] = (statrecord*)calloc((int)maxzeit[i] + 1, sizeof(statrecord));
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
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					id1 = zarray[(int)zn[i][j]];
					id2 = (int)FZT[i][j];
					if (RelCellArea != NULL)
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
		MPI_Type_contiguous(5,MPI_DOUBLE,&zonevectortype);
		MPI_Type_commit(&zonevectortype);

		void(*zonesumfunction)(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype) = &zonesumsurfrout;
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

			MPI_Allreduce(&(zoneFZT_QD[0][i][0]), tempzone, (int)(maxzeit[i] + 1), zonevectortype, mySumOp,MPI_COMM_WORLD);

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
				(*WSstream) << "WARNING:warning: error opening log-file "<<logfilename<<"[end-entry]";
				logdateiopen=false;
			}
			else logdateiopen=true;
			fprintf(logdatei,"Flow time distribution:\n");
		}

		for (int i = 1; i <= zonmax; i++)
		{
			if (mpi_rank == 0)
			{
				if (logdateiopen) fprintf(logdatei, "TG\t%4d\t(%3d\tZonen):", (int)zone[i].zonenwert, (int)maxzeit[i]);
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
			for (int z = 0; z <= tracer_count; z++)
			{
				for (int j = 1; j <= (int)maxzeit[i]; j++)
				{
					zoneFZT_QI[z][i][j].rel = zoneFZT_QD[z][i][j].rel = zoneFZT_QD[0][i][j].n / cellcount;
					zoneFZT_QI[z][i][j].zonenwert = zoneFZT_QD[z][i][j].zonenwert = j;
					if (z == 0)
					{
						if (mpi_rank == 0)
						{
							if (logdateiopen) fprintf(logdatei, "\t%5.3f", zoneFZT_QD[0][i][j].rel);
						}
						if (zoneFZT_QD[0][i][1].rel <= 0)
						{
							(*WSstream) << "FATAL: in flow-travel-time grid: the first time zone of subbasin " << zarray[i]
								<< " has no cells!" << "[end-entry]";
							throw(61);
						}
						zoneFZT_QD[0][i][j].n = 0;
					}
				}
			}
			if (mpi_rank == 0)
			{
				if (logdateiopen) fprintf(logdatei, "\n");
			}
		}
		if (mpi_rank == 0)
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
				(*WSstream) << "WARNING: log-file "<<logfilename<<" could not be opened"<<"[end-entry]";
			}
			else
			{
				fprintf(logdatei,"flow time distribution will not be calculated (no flow time grid provided)\n");
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
//#pragma optimize ("p",off)

/* fuer Zellen, die mindestens einen nodata-Wert als Rand haben kann keine
3x3-Matrix zur Berechnugn von slope, aspect, curvature aufgestellt werden.
Diese Kennwerte werden mit der folgenden Funktion aus den Werten der
direkten Nachbarn mit gueltigen Werten berechnet */
void surfaceroutingclass::getnachbar(int i,int j,double **Grid)
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


void surfaceroutingclass::GetPlanCurvature()
{
	if (mpistartrow == 1)
	{
		for (int j = 1; j <= ncols; j++)
		{
			PLANCURV[1][j] = nodata;
		}
	}
	if (mpiendrow == nrows)
	{
		for (int j = 1; j <= ncols; j++)
		{
			PLANCURV[nrows][j] = nodata;
		}
	}

	int sr = (mpistartrow == 1) ? mpistartrow + 1 : mpistartrow;
	int er = (mpiendrow == nrows) ? mpiendrow - 1 : mpiendrow;

	for (int i=sr;i<er;i++)
	{
		for (int j=2;j<ncols;j++)
		{
			if ( (DHM[i-1][j-1]!=nodata)&&
					(DHM[i-1][j  ]!=nodata)&&
					(DHM[i-1][j+1]!=nodata)&&
					(DHM[i  ][j-1]!=nodata)&&
					(DHM[i  ][j  ]!=nodata)&&
					(DHM[i  ][j+1]!=nodata)&&
					(DHM[i+1][j-1]!=nodata)&&
					(DHM[i+1][j  ]!=nodata)&&
					(DHM[i+1][j+1]!=nodata))
			{
				/* Ermitteln der Koefizienten des Gleichungssystems */
				double H = (DHM[i][j+1]-DHM[i][j-1])/(cellsize+cellsize);
				double G = (DHM[i-1][j]-DHM[i+1][j])/(cellsize+cellsize);
				double F = (-DHM[i-1][j-1]+DHM[i+1][j-1]+DHM[i-1][j+1]-DHM[i+1][j+1])/(4.0*cellsize*cellsize);
				double E = ((DHM[i][j-1]+DHM[i][j+1])/2-DHM[i][j])/(cellsize*cellsize);
				double D = ((DHM[i-1][j]+DHM[i+1][j])/2-DHM[i][j])/(cellsize*cellsize);
				if ((G!=0)||(H!=0)) PLANCURV[i][j] = 2*(D*H*H + E*G*G - F*G*H) / (G*G + H*H);   else PLANCURV[i][j]=0; 

			}
			else
			{
				PLANCURV[i][j]=nodata;
			}
			PLANCURV[i][1] = nodata; 
			PLANCURV[i][ncols] = nodata;
		}
	}
	/* Fehlwerte dort auffuellen, wo im Hoehenmodell Daten stehen -> naechster Nachbar oder Mittelwert der Nachbarn */
	for (int i=mpistartrow;i<=mpiendrow;i++)
	{
		for (int j=1;j<=ncols;j++)
		{
			if ((DHM[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(PLANCURV[i][j]==nodata))  
			{
				getnachbar(i,j,PLANCURV);
			}
		}
	}
}

void surfaceroutingclass::GetSlopeAspectAndCelltype(int i,int j,double &slope,double &aspect,int &celltype,int step)
{
	celltype = 0; // initializing to 0 which has no meaning at all.
	slope = 0; // default return value e.g. for river cells
	aspect = 0; // default return value e.g. for river cells
	double h[9]; // local elevations, taken from DEM with additional water amounts on top of it
	double head[9];
	double sl[9];
	
	if ((method==3)||(method==4)||(method==5)) // 3=multiple flow paths dynamically based on water level, 4 = same with single flow path, 5 = simple v calculation without iteration (suited for back water conditions)
	{
		head[0] = HEAD[i][j];

		if ((i>1)&&(zn[i-1][j]>0))		{
			head[1] = HEAD[i - 1][j];
		} else { 
			head[1] = 0;
		}

		if ((i>1)&&(j<ncols)&&(zn[i-1][j+1]>0))		{
			head[2] = HEAD[i - 1][j + 1];
		} else { 
			head[2] = 0;
		}

		if ((j<ncols)&&(zn[i][j+1]>0))		{
			head[3] = HEAD[i][j + 1];
		} else { 
			head[3] = 0;
		}

		if ((i<nrows)&&(j<ncols)&&(zn[i+1][j+1]>0))		{
			head[4] = HEAD[i + 1][j + 1];
		} else { 
			head[4] = 0;
		}

		if ((i<nrows)&&(zn[i+1][j]>0))		{
			head[5] = HEAD[i + 1][j];
		} else { 
			head[5] = 0;
		}

		if ((i<nrows)&&(j>1)&&(zn[i+1][j-1]>0))		{
			head[6] = HEAD[i + 1][j - 1];
		} else { 
			head[6] = 0;
		}

		if ((j>1)&&(zn[i][j-1]>0))		{
			head[7] = HEAD[i][j - 1];
		} else { 
			head[7] = 0;
		}

		if ((i>1)&&(j>1)&&(zn[i-1][j-1]>0))		{
			head[8] = HEAD[i - 1][j - 1];
		} else { 
			head[8] = 0;
		}

	}
	else
	{
		// no water level regarded
		head[0] = head[1] = head[2] = head[3] = head[4] = head[5] = head[6] = head[7] = head[8] = 0.0;
	}

	// analyse location of the cell relative to their neihgbors
	// get thy type of the cell
	// 1 = on a river cell -> immediate runoff
	// 2 = in the middle of nowhere = inmidst of the catchment but not on a river, multiple flow path possible (up to three)
	// 3 = a cell with at least one nodata-cell as neighbor or at the border of the model grid
	if ((RIVWIT[i][j]>0)&&(RIVDEP[i][j]>0)) {
		celltype = 1; 
	} else 	if ((i>1)&&(j>1)&&(i<nrows)&&(j<ncols)
			&& (zn[i-1][j-1]>0) && (zn[i-1][j]>0) && (zn[i-1][j+1]>0)
			&& (zn[i][j-1]>0)                     && (zn[i][j+1]>0)
			&& (zn[i+1][j-1]>0) && (zn[i+1][j]>0) && (zn[i+1][j+1]>0))	{
		celltype = 2; // this is preliminary, since the slopes are not checked yet.
	} else {
		celltype = 3;
	}

	if (celltype == 2)
	{
		// this is a cell inmidst the valid model domain
		h[0]=DHM[i  ][j  ]+head[0];
		h[1]=DHM[i-1][j  ]+head[1];
		h[2]=DHM[i-1][j+1]+head[2];
		h[3]=DHM[i  ][j+1]+head[3];
		h[4]=DHM[i+1][j+1]+head[4];
		h[5]=DHM[i+1][j  ]+head[5];
		h[6]=DHM[i+1][j-1]+head[6];
		h[7]=DHM[i  ][j-1]+head[7];
		h[8]=DHM[i-1][j-1]+head[8];

		// looking for the steepest slope and its exposition
		slope=0;

		sl[1]=(h[0]-h[1])/cellsize; // 1 is to the north, then the other indices are following clockwise
		sl[2]=(h[0]-h[2])/(cellsize*1.414);
		sl[3]=(h[0]-h[3])/cellsize;
		sl[4]=(h[0]-h[4])/(cellsize*1.414);
		sl[5]=(h[0]-h[5])/cellsize;
		sl[6]=(h[0]-h[6])/(cellsize*1.414);
		sl[7]=(h[0]-h[7])/cellsize;
		sl[8]=(h[0]-h[8])/(cellsize*1.414);

		double maxslope=sl[1];
		int maxslope_index=1;

		for (int cnt=2;cnt<=8;cnt++)
		{
			if (sl[cnt]>maxslope)
			{
				maxslope=sl[cnt];
				maxslope_index=cnt;
			}
		}
		slope = maxslope;
		aspect = (maxslope_index-1.0)*45.0;

		int maxcnt=0;
		for (int cnt=2;cnt<=8;cnt++)
		{
			if (sl[cnt]==maxslope) maxcnt++;
		}
		if (maxcnt>1) // more than one possible directions, so the aspect must be estimated ba another approach
		{
			double H = (h[3]-h[7])/(cellsize+cellsize);
			double G = (h[1]-h[5])/(cellsize+cellsize);
			/*        if ((G!=0)||(H!=0)) 
			{
				slope = 180.0/Pi*atan(sqrt(G*G + H*H));           
			}
			else 
			{
				slope=0;
			}
	*/
			if ( H!=0)          
			{
//				aspect = 180.0 - 180.0 / Pi*atan(G / H) + 90.0*(H / fabs(H));
				aspect = 180.0 - 180.0 / Pi*atan(G / H) + 90.0*(H / std::abs(H));
			}
			else 
			{
				aspect = G>0 ? 180.0 : 0;
			}
		}
	} 

	if (celltype==3) 
	{
		h[0]=DHM[i  ][j  ]+head[0];

		if ((i>1)&&(zn[i-1][j]>0))		{
			h[1]=DHM[i-1][j]+head[1];
			sl[1]=(h[0]-h[1])/cellsize; // 1 is to the north, then the other indices are following clockwise
		} else { 
			sl[1]=-9999;
		}

		if ((i>1)&&(j<ncols)&&(zn[i-1][j+1]>0))		{
			h[2]=DHM[i-1][j+1]+head[2];
			sl[2]=(h[0]-h[2])/(cellsize*1.414);
		} else { 
			sl[2]=-9999;
		}

		if ((j<ncols)&&(zn[i][j+1]>0))		{
			h[3]=DHM[i][j+1]+head[3];
			sl[3]=(h[0]-h[3])/cellsize;
		} else { 
			sl[3]=-9999;
		}

		if ((i<nrows)&&(j<ncols)&&(zn[i+1][j+1]>0))		{
			h[4]=DHM[i+1][j+1]+head[4];
			sl[4]=(h[0]-h[4])/(cellsize*1.414);
		} else { 
			sl[4]=-9999;
		}

		if ((i<nrows)&&(zn[i+1][j]>0))		{
			h[5]=DHM[i+1][j]+head[5];
			sl[5]=(h[0]-h[5])/cellsize;
		} else { 
			sl[5]=-9999;
		}

		if ((i<nrows)&&(j>1)&&(zn[i+1][j-1]>0))		{
			h[6]=DHM[i+1][j-1]+head[6];
			sl[6]=(h[0]-h[6])/(cellsize*1.414);
		} else { 
			sl[6]=-9999;
		}

		if ((j>1)&&(zn[i][j-1]>0))		{
			h[7]=DHM[i][j-1]+head[7];
			sl[7]=(h[0]-h[7])/cellsize;
		} else { 
			sl[7]=-9999;
		}

		if ((i>1)&&(j>1)&&(zn[i-1][j-1]>0))		{
			h[8]=DHM[i-1][j-1]+head[8];
			sl[8]=(h[0]-h[8])/(cellsize*1.414);
		} else { 
			sl[8]=-9999;
		}

		// looking for the steepest slope and its exposition
		slope=0;

		double maxslope=sl[1];
		int maxslope_index=1;

		for (int cnt=2;cnt<=8;cnt++)
		{
			if (sl[cnt]>maxslope)
			{
				maxslope=sl[cnt];
				maxslope_index=cnt;
			}
		}
		slope = maxslope;
		aspect = (maxslope_index-1.0)*45.0;

		int maxcnt=0;
		for (int cnt=2;cnt<=8;cnt++)
		{
			if (sl[cnt]==maxslope) maxcnt++;
		}
		if (maxcnt>1)
		{
			// replace missing values by the neighbors to the other side
			h[0]=DHM[i][j]+head[0];

			if ((i>1)&&(zn[i-1][j]>0))		{
				h[1]=DHM[i-1][j]+head[1];
			} else { // h1 ha s no valid data, so get it from h5 or from the actual cell
				h[1]=((i<nrows)&&(zn[i+1][j]>0))? (h[0]+h[0]-DHM[i+1][j]) : h[0];
			}

			if ((j<ncols)&&(zn[i][j+1]>0))		{
				h[3]=DHM[i][j+1]+head[3];
			} else { // h3 has no valid data, so get it from h7 or from the actual cell
				h[3]=((j>1)&&(zn[i][j-1]>0))? (h[0]+h[0]-DHM[i][j-1]) : h[0];
			}

			if ((i<nrows)&&(zn[i+1][j]>0))		{
				h[5]=DHM[i+1][j]+head[5];
			} else { // h5 has no valid data, so get it from h1 or from the actual cell
				h[5]=((i>1)&&(zn[i-1][j]>0))? (h[0]+h[0]-DHM[i-1][j]) : h[0];
			}

			if ((j>1)&&(zn[i][j-1]>0))		{
				h[7]=DHM[i][j-1]+head[7];
			} else { // h7 has no valid data, so get it from h3 or from the actual cell
				h[7]=((j<ncols)&&(zn[i][j+1]>0))? (h[0]+h[0]-DHM[i][j+1]) : h[0];
			}

			double H = (h[3]-h[7])/(cellsize+cellsize);
			double G = (h[1]-h[5])/(cellsize+cellsize);
			if ((G!=0)||(H!=0)) 
			{
				slope = 180.0/Pi*atan(sqrt(G*G + H*H));           
			}
			else 
			{
				slope=0;
			}
			if ( H!=0)          
			{
//				aspect = 180.0 - 180.0 / Pi*atan(G / H) + 90.0*(H / fabs(H));
				aspect = 180.0 - 180.0 / Pi*atan(G / H) + 90.0*(H / std::abs(H));
			}
			else 
			{
				aspect = G>0 ? 180.0 : 0;
			}
		}
	}
}

void surfaceroutingclass::GetMultipleFlowPath(int i,int j,double slope,double aspect,int &celltype,int step)
{
	double sum;
	double w1,w2,w3,alpha,beta;
	int MainDir=0;
	double Old_Weight[9];
	double totalin[9];
	// we have to check now if there is at least one positive slope to any cell within the same catchment or to a river cell outside the catchment
	// Remark: if no positive slope to a cell inside the same catchment can be found but to a cell neighboring the actual cell in another subcatchment
	// which is a river cell (in LNK-Grid only, not in RIVDEP or RIVWIT), then the actual cell is assumed to be the outlet of the actual subcatchment and will get type 1
	if (celltype==1) return; // river cells are not taken into account for surface routing
//	if ((celltype==2)||(celltype==3)) // normal cells will have opportunities to discharge into three neighboring cells
//	{  
	for (int cnt=1;cnt<=8;cnt++) 
	{
		Old_Weight[cnt]=WEIGHT[i][j][cnt];
		WEIGHT[i][j][cnt]=0.0;
	}

	if ((method==3)||(method==4)||(method==5)) // 3=multiple flow paths dynamic based on water level, 4 = same with single flow path, 5 = simple v calculation, else like 4
	{
		totalin[0]=HEAD[i  ][j  ];

		if (celltype == 2) // normal cell with 8 valid neighbors
		{ // no checks required, makes the code faster
			totalin[1] = HEAD[i - 1][j    ];
			totalin[2] = HEAD[i - 1][j + 1];
			totalin[3] = HEAD[i    ][j + 1];
			totalin[4] = HEAD[i + 1][j + 1];
			totalin[5] = HEAD[i + 1][j    ];
			totalin[6] = HEAD[i + 1][j - 1];
			totalin[7] = HEAD[i    ][j - 1];
			totalin[8] = HEAD[i - 1][j - 1];
		}
		else
		{
			if ((i > 1)     &&                (zn[i - 1][j    ] > 0))		{	totalin[1] = HEAD[i - 1][j    ];		} else {	totalin[1] = 0;		}
			if ((i > 1)     && (j < ncols) && (zn[i - 1][j + 1] > 0))		{	totalin[2] = HEAD[i - 1][j + 1];		} else {	totalin[2] = 0;		}
			if ((j < ncols) &&                (zn[i    ][j + 1] > 0))		{	totalin[3] = HEAD[i    ][j + 1];		} else {	totalin[3] = 0;		}
			if ((i < nrows) && (j < ncols) && (zn[i + 1][j + 1] > 0))		{	totalin[4] = HEAD[i + 1][j + 1];		} else {	totalin[4] = 0;		}
			if ((i < nrows) &&                (zn[i + 1][j    ] > 0))		{	totalin[5] = HEAD[i + 1][j    ];		} else {	totalin[5] = 0;		}
			if ((i < nrows) && (j > 1)     && (zn[i + 1][j - 1] > 0))		{	totalin[6] = HEAD[i + 1][j - 1];		} else {	totalin[6] = 0;		}
			if ((j > 1)     &&                (zn[i    ][j - 1] > 0))		{	totalin[7] = HEAD[i    ][j - 1];		} else {	totalin[7] = 0;		}
			if ((i > 1)     && (j > 1)     && (zn[i - 1][j - 1] > 0))		{	totalin[8] = HEAD[i - 1][j - 1];		} else {	totalin[8] = 0;		}
		}
	}
	else
	{
		// no water level regarded
		totalin[0]=totalin[1]=totalin[2]=totalin[3]=totalin[4]=totalin[5]=totalin[6]=totalin[7]=totalin[8]=0.0;
	}

	if ((slope!=0)&&(celltype==2))
	{ // there is a main slope with an aspect, so we can calculate 3 main runoff directions
		// N =1 (NW, NE as minor flow directions)
		// NE=2 (N and E as minor flow directions)
		// E =3 (NW and SE as minor flow directions)
		// SE=4 (E and S as minor flow directions)
		// S =5 (SE and SW as minor flow directions)
		// SW=6 (S and W as minor flow directions)
		// W =7 (SW and NW as minor flow directions)
		// NW=8 (W and N as minor flow directions)
		// aspect is counting clockwise from north from 0...360
		if (aspect > 337.5) aspect -= 360.0;
		MainDir = int(floor((aspect + 22.5) / 45.0)) + 1;
		// rotating coordinate system for main direction showing to N
		aspect -= (MainDir - 1)*45.0;
// that's the old version, probably a bit slower
//		if ((aspect<22.5)||(aspect>=337.5)) { // N (NW, NE as minor flow directions)
//			MainDir=1; // N
//		} else if ((aspect>=22.5)&&(aspect<67.5)) { // NE (N and E as minor flow directions)
//			MainDir=2; // NE
//			aspect-=45.0; // rotating coordinate system for main direction showing to N
//		} else if ((aspect>=67.5)&&(aspect<112.5)) { // E (NW and SE as minor flow directions)
//			MainDir=3; // E
//			aspect-=90.0; // rotating coordinate system for main direction showing to N
//		} else if ((aspect>=112.5)&&(aspect<157.5)) { // SE (E and S as minor flow directions)
//			MainDir=4; // SE
//			aspect-=135.0;// rotating coordinate system for main direction showing to N
//		} else if ((aspect>=157.5)&&(aspect<202.5)) { // S (SE and SW as minor flow directions)
//			MainDir=5; // S
//			aspect-=180.0;// rotating coordinate system for main direction showing to N
//		} else if ((aspect>=202.5)&&(aspect<247.5)) { // SW (S and W as minor flow directions)
//			MainDir=6.0; // SW
//			aspect-=225.0;// rotating coordinate system for main direction showing to N
//		} else if ((aspect>=247.5)&&(aspect<292.5)) { // W (SW and NW as minor flow directions)
//			MainDir=7; // W
//			aspect-=270.0;// rotating coordinate system for main direction showing to N
//		} else { // if ((aspect>=22.5)&&(aspect<67.5)) { // NW (W and N as minor flow directions)
//			MainDir=8; // NW
//			aspect-=315.0;// rotating coordinate system for main direction showing to N
//		}
		w1 = w3 = 0;
		w2 = 1 - 0;
		WEIGHT[i][j][MainDir] = w2; // thats all we need for single flow paths
		SRFCVEL1[i][j] = SRFCVEL[i][j][0]=MainDir;

		if ((method == 1) || (method == 3))
		{
			if (aspect > 22.5)
			{
				alpha = aspect - 360.0;
			}
			else
			{
				alpha = aspect;
			}
			if (aspect >= 0)
			{
				beta = 45.0 - alpha;
			}
			else
			{
				beta = 45.0 + alpha;
			}
			alpha = alpha*0.017453; // 0.017453 = Pi/180.0 -> conversion grad in rad
			beta = beta *0.017453;
//			w1 = (alpha > 0) ? (fabs(cos(alpha)*cos(beta))) : (fabs(cos(alpha)*sin(beta))); // minor direction to the left
//			w2 = cos(alpha); // main direction
//			w3 = (alpha > 0) ? (fabs(cos(alpha)*sin(beta))) : (fabs(cos(alpha)*cos(beta))); // minor direction to the left
			w1 = (alpha > 0) ? (std::abs(cos(alpha)*cos(beta))) : (std::abs(cos(alpha)*sin(beta))); // minor direction to the left
			w2 = cos(alpha); // main direction
			w3 = (alpha > 0) ? (std::abs(cos(alpha)*sin(beta))) : (std::abs(cos(alpha)*cos(beta))); // minor direction to the left
			if (w1 < 1e-4) w1 = 0;
			if (w2 < 1e-4) w2 = 0;
			if (w3 < 1e-4) w3 = 0;
			// single flow path by setting w1 and w3 to 0
			//		if ((method==2)||(method==4)||(method==5))
			//		{
			//			w1=w3=0.0;
			//		}
			sum = w1 + w2 + w3;
			// normalizing weights to 1.0
			w1 = w1 / sum;
			w2 = w2 / sum;
			w3 = w3 / sum;
			if (MainDir == 1) { // N
				WEIGHT[i][j][1] = w2;
				if (PLANCURV[i][j] <= 0)	{ // multiple flow path' are calculated for divergence areas only
					WEIGHT[i][j][8] = w1; WEIGHT[i][j][2] = w3;
				}
			}
			else if (MainDir == 2) { // NE
				WEIGHT[i][j][2] = w2;
				if (PLANCURV[i][j] <= 0)	{ // multiple flow path' are calculated for divergence areas only
					WEIGHT[i][j][1] = w1; WEIGHT[i][j][3] = w3;
				}
			}
			else if (MainDir == 3) { // E
				WEIGHT[i][j][3] = w2;
				if (PLANCURV[i][j] <= 0)	{ // see above...
					WEIGHT[i][j][2] = w1; WEIGHT[i][j][4] = w3;
				}
			}
			else if (MainDir == 4) { // SE
				WEIGHT[i][j][4] = w2;
				if (PLANCURV[i][j] <= 0)	{
					WEIGHT[i][j][3] = w1; WEIGHT[i][j][5] = w3;
				}
			}
			else if (MainDir == 5) { // S
				WEIGHT[i][j][5] = w2;
				if (PLANCURV[i][j] <= 0)	{
					WEIGHT[i][j][4] = w1; WEIGHT[i][j][6] = w3;
				}
			}
			else if (MainDir == 6) { // SW
				WEIGHT[i][j][6] = w2;
				if (PLANCURV[i][j] <= 0)	{
					WEIGHT[i][j][5] = w1; WEIGHT[i][j][7] = w3;
				}
			}
			else if (MainDir == 7) { // W
				WEIGHT[i][j][7] = w2;
				if (PLANCURV[i][j] <= 0)	{
					WEIGHT[i][j][6] = w1; WEIGHT[i][j][8] = w3;
				}
			}
			else { // NW
				WEIGHT[i][j][8] = w2;
				if (PLANCURV[i][j] <= 0)	{
					WEIGHT[i][j][7] = w1; WEIGHT[i][j][1] = w3;
				}
			}
		}

		double dhm=totalin[0]+DHM[i][j];

		if ((WEIGHT[i][j][1]>0)&&((dhm-DHM[i-1][j  ] - totalin[1])<0)) {WEIGHT[i][j][1]=0;}
		if ((WEIGHT[i][j][2]>0)&&((dhm-DHM[i-1][j+1] - totalin[2])<0)) {WEIGHT[i][j][2]=0;}
		if ((WEIGHT[i][j][3]>0)&&((dhm-DHM[i  ][j+1] - totalin[3])<0)) {WEIGHT[i][j][3]=0;}
		if ((WEIGHT[i][j][4]>0)&&((dhm-DHM[i+1][j+1] - totalin[4])<0)) {WEIGHT[i][j][4]=0;}
		if ((WEIGHT[i][j][5]>0)&&((dhm-DHM[i+1][j  ] - totalin[5])<0)) {WEIGHT[i][j][5]=0;}
		if ((WEIGHT[i][j][6]>0)&&((dhm-DHM[i+1][j-1] - totalin[6])<0)) {WEIGHT[i][j][6]=0;}
		if ((WEIGHT[i][j][7]>0)&&((dhm-DHM[i  ][j-1] - totalin[7])<0)) {WEIGHT[i][j][7]=0;}
		if ((WEIGHT[i][j][8]>0)&&((dhm-DHM[i-1][j-1] - totalin[8])<0)) {WEIGHT[i][j][8]=0;}
			
		// normalizing weights again to 1
		sum=0;
		for (int d=1;d<=8;d++) sum=sum+WEIGHT[i][j][d];
		if (sum==0) 
		{
			celltype=3;
		}
		else
		{
			for (int d=1;d<=8;d++)
			{
				WEIGHT[i][j][d]=WEIGHT[i][j][d]/sum;
//				fpe_check(i,j,celltype,d,WEIGHT[i][j][d],sum,w1,w2,w3);
			}
		}
	}
	// also a changed value of celltype is handled here (e.g. if slope was calculated downwards but it is upwards, in reality)
	if ((slope==0)||(celltype==3))
	{ // also cells with slope = 0 must be investigated because they may have an even slope to all directions!
		double dhm=totalin[0]+DHM[i][j];
		if ((i>1    )           && (zn[i-1][j  ])!=nodata) 
		{ w1=dhm-DHM[i-1][j]  - totalin[1];  } else w1=0; if (w1>1.0e-5) WEIGHT[i][j][1]=w1;

		if ((i>1    )&&(j<ncols)&& (zn[i-1][j+1])!=nodata) 
		{ w1=dhm-DHM[i-1][j+1]- totalin[2];} else w1=0; if (w1>1.0e-5) WEIGHT[i][j][2]=w1/1.4142135;

		if (           (j<ncols)&& (zn[i  ][j+1])!=nodata) 
		{ w1=dhm-DHM[i][j+1]  - totalin[3];  } else w1=0; if (w1>1.0e-5) WEIGHT[i][j][3]=w1;

		if ((i<nrows)&&(j<ncols)&& (zn[i+1][j+1])!=nodata) 
		{ w1=dhm-DHM[i+1][j+1]- totalin[4];} else w1=0; if (w1>1.0e-5) WEIGHT[i][j][4]=w1/1.4142135;

		if ((i<nrows)           && (zn[i+1][j  ])!=nodata) 
		{ w1=dhm-DHM[i+1][j]  - totalin[5];  } else w1=0; if (w1>1.0e-5) WEIGHT[i][j][5]=w1;

		if ((i<nrows)&&(j>1    )&& (zn[i+1][j-1])!=nodata) 
		{ w1=dhm-DHM[i+1][j-1]- totalin[6];} else w1=0; if (w1>1.0e-5) WEIGHT[i][j][6]=w1/1.4142135;

		if (           (j>1    )&& (zn[i  ][j-1])!=nodata) 
		{ w1=dhm-DHM[i][j-1]  - totalin[7];  } else w1=0; if (w1>1.0e-5) WEIGHT[i][j][7]=w1;

		if ((i>1    )&&(j>1    )&& (zn[i-1][j-1])!=nodata) 
		{ w1=dhm-DHM[i-1][j-1]- totalin[8];} else w1=0; if (w1>1.0e-5) WEIGHT[i][j][8]=w1/1.4142135;

		sum=0;
		double maxw=0;
		int maxdir=0;
		if ((method==2)||(method==4)||(method==5)) // single flow path in every case
		{
			for (int cnt=1;cnt<=8;cnt++) 
			{
				if (WEIGHT[i][j][cnt]>maxw)
				{
					maxw=WEIGHT[i][j][cnt];
					maxdir=cnt;
				}
				WEIGHT[i][j][cnt]=0;
			}
			if (maxdir>0)
			{
				WEIGHT[i][j][maxdir]=1;
			}
		}
		else
		{ // multiple flow paths
			for (int cnt=1;cnt<=8;cnt++) 
			{
				if (WEIGHT[i][j][cnt]>maxw)
				{
					maxw=WEIGHT[i][j][cnt];
					maxdir=cnt;
				}
				sum+=WEIGHT[i][j][cnt];
			}
			if (sum>0)
			{
				for (int cnt=1;cnt<=8;cnt++) 
				{
					WEIGHT[i][j][cnt]=WEIGHT[i][j][cnt]/sum;
//					fpe_check(i,j,celltype,cnt,WEIGHT[i][j][cnt],sum,0,0,0);
				}
			}
		}
		SRFCVEL1[i][j] = SRFCVEL[i][j][0] = maxdir;
	}
	SRFCVEL[i][j][9]=DHM[i][j]+SSSU[i][j]/1000.0;
	// for each additional sub time step we check for each cell, if the original weights have changed due to 
	// water level changes. If not, then the original SRFCVEL-Values may be used without doing iterations again and again
	WEIGHT[i][j][0]=0;
	if (step>0)
	{
		// changing WEIGHTS are indicated by a 0 in the zero-element of the WEIGHT-vector
		int changed=0;
		for (int cnt=1;cnt<=8;cnt++)
		{
//			if (fabs(Old_Weight[cnt] - WEIGHT[i][j][cnt]) <= 0.01)
			if (std::abs(Old_Weight[cnt] - WEIGHT[i][j][cnt]) <= 0.01)
			{
				// changes in WEIGHTs are below one per cent, so we assume the velocities to remain constant
				WEIGHT[i][j][0]+=1.0;
			}
		}

	}
//	}
}

// calculating hydraulic roughness k_st or M in m^(1/3)/s --> the so called Manning-Strickler-Coefficient
// i,j : row and column of the actual cell
// v: flow velocity (1st iteration value which may be needed to further enhance M)
// luse: landuse code for the actual cell
// SR_cultivation: type of cultivation: 1=conserving;2=non-conserving,3=other (like forest or grassland)
// RH: hydraulic radius in m -> flow depth as first iteration
// VCF: vegetation covered fraction
// Ie: energy slope(?) (Energieliniengefälle)
double surfaceroutingclass::GetRoughness(int i,int j,double v,int luse,int SR_cultivation,double RH,double VCFlocal,double Ie)
{
	double kst=this->LUtable[luse].GetSR_K_ST();
	
	double ReynoldsNumber=0;
	double M=kst; // the return value defaults to the K_st-Value in the control file for other cultivation techniques
	double a,b,c; // coefficients for calculation of n = 1/k for conserving cultivation
	double n;

	if (SR_cultivation==1) { // conserving
		ReynoldsNumber = 4.0e6*v*RH; // 4*v*RH/1e-6
		if (ReynoldsNumber<20000) {
			LUtable[luse].GetSR_rough_coeff1(a,b,c);
		}	else	{
			LUtable[luse].GetSR_rough_coeff2(a,b,c);
		}
		n=(a*pow((double)(VCFlocal*100.0),(double)b))/pow((double)ReynoldsNumber,(double)c);
		M=1.0/n;
	} else if (SR_cultivation==2) { // nonconserving
		double LambdaBed=0,LambdaCrop=0,LambdaTotal=0;
		const double g=9.8159;
		double StemDiameter=LUtable[luse].GetSR_StemDiameter(WSstream,tag,monat,jahr,DHM[i][j],i,j);
		double CropDistX=LUtable[luse].GetSR_CropDistX();
		double CropDistY=LUtable[luse].GetSR_CropDistY();
		double Cw=0;
		if ((StemDiameter>=0.001)&&(RH>=dMinFlowDepth))
		{
			// change LambdaCrop according to crop properties
			LambdaCrop=0;
			ReynoldsNumber=v*StemDiameter*1.0e6;
			if (ReynoldsNumber<800) {
				Cw=3.07*pow((double)ReynoldsNumber,(double)-0.168);
			} else if (ReynoldsNumber<8000) {
				Cw=1.0;
			} else if (ReynoldsNumber<10000) {
				Cw=1.2;
			} else {
				Cw=LUtable[luse].GetSR_C_WR();
			}
			// now iterating A_nl and A_nb
			double ANL=50.0*StemDiameter;
			double ANLold=0;
			double deltaANL=ANL;
			int iter=0;
			// check for convergence of the following iteration (will be usccessfull for combinations of low slope and high velocities only!
			if (Cw*Ie/(v*v)>0.04)
			{ // will converge only for iter > 20
				Ie=(v*v)/Cw/25.0;
			}
			while ((deltaANL>this->dMaxDeltaANL)&&(iter < this->iMaxIterANL))
			{
				ANL=128.9*Cw*StemDiameter*pow((double)(1.0+(g*ANL*Ie)/(v*v/2.0)),(double)-2.14);
				//				ANL=255.7*Cw*StemDiameter*(1.0+(2.0*g*ANL*Ie)/(v*v));
//				deltaANL = fabs(ANLold - ANL);
				deltaANL = std::abs(ANLold - ANL);
				ANLold = ANL;
				iter++;
			}
			double ANB = 0.24*pow((double)ANL,(double)0.59)*pow((double)(Cw*StemDiameter),(double)0.41);
			double Vn_Vvor2 = 1.15*pow((double)(ANL/CropDistX),(double)-0.48)+
			0.5*pow((double)(ANB/CropDistY),(double)1.1);
			// solving the 3rd-order equation Fr=(y*(y*y-1)))/(2*(y-CropDistY/(CropDistY-Stemdiameter)))

			/*			double Fr = v/sqrt(g*RH);
			double coeffp  = 3*(-1-2*Fr);
			double coeffq  = 27*(2*CropDistY*Fr)/(CropDistY-StemDiameter);
			double coeffu  = 0.5*pow((double)(-4.0*coeffq + 4*sqrt(coeffq*coeffq+4*coeffp*coeffp*coeffp)),(double)0.33333);
			double coeffv  = 0.5*pow((double)(-4.0*coeffq - 4*sqrt(coeffq*coeffq+4*coeffp*coeffp*coeffp)),(double)0.33333);
			double coeffz  = coeffu+coeffv;
			double coeffy  = coeffz/3.0;
			double dCw = 2.0/(Fr*Fr)*(1.0-coeffy);
*/
			double coeffy = 0.6;
			double Fr = v/sqrt(g*RH);
			double y_old=0;
			iter=0;
			// coeffy will converge only if a) distance between corps is at least 5 times the diameter of the stems
			// and b) Froude number times (cropdisty/(cropdisty-stemdiameter))^3 is less than 2/3
			if (CropDistY<=5*StemDiameter)
			{
				CropDistY=StemDiameter*5.0;
			}
			if (Fr*pow((double)(CropDistY/(CropDistY-StemDiameter)),(double)3.0)>0.666)
			{
				Fr = 2.0/(3.0*pow((double)(CropDistY/(CropDistY-StemDiameter)),(double)3.0));
			}
			if (Fr<0.05) 
			Fr=0.05;
//			while ((fabs(y_old-coeffy)>0.0001)&&(iter < this->iMaxIterANL))
			while ((std::abs(y_old - coeffy)>0.0001) && (iter < this->iMaxIterANL))
			{
				y_old=coeffy;
				//				coeffy=CropDistY/(CropDistY-StemDiameter)+(coeffy*(coeffy*coeffy-1.0))/(2.0*Fr);
				coeffy = (2 * Fr * (coeffy - CropDistY/(CropDistY-StemDiameter)))/(coeffy*coeffy-1.0);
				iter++;
			}
			double dCw = 2.0/(Fr*Fr)*(1.0-coeffy);
			double Cwr = 1.31 * Cw * Vn_Vvor2 + dCw;
			LambdaCrop = 4.0*RH*StemDiameter/(CropDistX*CropDistY)*Cwr;
		}
		else
		{
			LambdaCrop=0;
		}
		LambdaBed = (8.0 * g)/(pow((double)RH,(double)0.3333)*kst*kst);
		LambdaTotal=LambdaBed+LambdaCrop;

		M = (8.0*g)/(LambdaTotal*pow((double)RH,(double)0.3333));
		if (M>0)
		{
			M=sqrt(M);
		}
		else
		{	
			M=0.001; // should be obvious in model results as extremely low flow rates
		}
	} else {
		M=kst;
	}
	return M;
}

void surfaceroutingclass::GetFlowVelocities(int i,int j,int step,int NextMajorStep,int celltype)
{
//	if ((WEIGHT[i][j][0]<8)||(step==0)||(step==NextMajorStep))
	if ((WEIGHT[i][j][0]<8)||(step==0)||(step==NextMajorStep)||(method==5))
	{
		int luse=-9999;
		int reduction=0;
		while ((luse==-9999)&&(reduction<iNumVegetationLayers))
		{
			luse=LUSE[iNumVegetationLayers-reduction][i][j];
			reduction++;
		}
		int SR_cultivation = LUtable[luse].getSR_cultivation();
		//		double VCFlocal=VCF[iNumVegetationLayers][i][j];
		double VCFlocal = 0.5;
		if ((SR_cultivation == 1) || (SR_cultivation == 2))
		{
			// time consuming function, so call it only when really needed
			VCFlocal = LUtable[luse].GetSR_ResidCovrg(WSstream, tag, monat, jahr, DHM[i][j], iNumVegetationLayers, i, j);
		}
		FLUX[i][j][0]=FLUX1[i][j]=0;
		double totalin[9];
		
		if ((method==3)||(method==4)||(method==5)) // 3=multiple flow paths dynamic based on water level, 4 = same with single flow path, 5 = like 4 but simple v calculation
		{
			totalin[0]=HEAD[i  ][j  ];
			
			if ((i>1)    &&           (zn[i-1][j  ]>0))		{totalin[1]=HEAD[i-1][j  ];} else {totalin[1]=0;}
			if ((i>1)    &&(j<ncols)&&(zn[i-1][j+1]>0))		{totalin[2]=HEAD[i-1][j+1];} else {totalin[2]=0;}
			if ((j<ncols)&&           (zn[i  ][j+1]>0))		{totalin[3]=HEAD[i  ][j+1];} else {totalin[3]=0;}
			if ((i<nrows)&&(j<ncols)&&(zn[i+1][j+1]>0))		{totalin[4]=HEAD[i+1][j+1];} else {totalin[4]=0;}
			if ((i<nrows)&&           (zn[i+1][j  ]>0))		{totalin[5]=HEAD[i+1][j  ];} else {totalin[5]=0;}
			if ((i<nrows)&&(j>1)    &&(zn[i+1][j-1]>0))		{totalin[6]=HEAD[i+1][j-1];} else {totalin[6]=0;}
			if ((j>1)    &&           (zn[i  ][j-1]>0))		{totalin[7]=HEAD[i  ][j-1];} else {totalin[7]=0;}
			if ((i>1)    &&(j>1)    &&(zn[i-1][j-1]>0))		{totalin[8]=HEAD[i-1][j-1];} else {totalin[8]=0;}
		}
		else
		{
			// no water level regarded
			totalin[0]=totalin[1]=totalin[2]=totalin[3]=totalin[4]=totalin[5]=totalin[6]=totalin[7]=totalin[8]=0.0;
		}
		double Qmax = 0;
		for (int cnt=1;cnt<=8;cnt++)
		{
			if ((WEIGHT[i][j][cnt]>0)&&(HEAD[i][j]>0))
			{
				// first guess for velocity: old value or 0.1m/s
				if (SRFCVEL[i][j][cnt]<=0)
				{
					SRFCVEL[i][j][cnt]=0.1;
				}
				FLUX[i][j][cnt]=0.1;

				double S; // slope
				double dist = cellsize;
				double elevdiff = 0; // will be sued for computing Qmax. Is always the volume that can be filled up according to the actual water levels and elevations of the two cells where waterflux is calculated
				double elevdiff2; // maximum elevation difference to the neighboring cell if that cell is a river cell (i.e. a sink)
				// getting slope from flow direction
				if (cnt==1) { // N
					if (celltype == 2)
					{
						S = (DHM[i][j] <= (DHM[i - 1][j] + totalin[cnt])) ? // cell is of type 2 means it has all valid neighbors. second criterion is to make sure there is a gradient from the source to the target cell
							((DHM[i + 1][j] + totalin[5]) - (DHM[i - 1][j] + totalin[cnt])) / (cellsize * 2) : // if both conditions are true: take the gradient second order from the cell to the north top the cell to the south
							((DHM[i][j] + totalin[0]) - (DHM[i - 1][j] + totalin[cnt])) / cellsize; // otherwise take the noprmal gradient from the cell to the north to the actual cell
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i - 1][j] + totalin[cnt])) / cellsize; // otherwise take the noprmal gradient from the cell to the north to the actual cell
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j]);
					}
				} else if (cnt==2) { // NE				   
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i - 1][j + 1] + totalin[cnt]))) ?
							((DHM[i + 1][j - 1] + totalin[6]) - (DHM[i - 1][j + 1] + totalin[cnt])) / (cellsize * 2 * 1.4142135) :
							((DHM[i][j] + totalin[0]) - (DHM[i - 1][j + 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j + 1] + totalin[cnt]);
					}
					else
					{
						S =  ((DHM[i][j] + totalin[0]) - (DHM[i - 1][j + 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j + 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j + 1]);
					}
					dist = cellsize*1.4142135;
				}
				else if (cnt == 3) { // E		    	   
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i][j + 1] + totalin[cnt]))) ?
							((DHM[i][j - 1] + totalin[7]) - (DHM[i][j + 1] + totalin[cnt])) / (cellsize * 2) :
							((DHM[i][j] + totalin[0]) - (DHM[i][j + 1] + totalin[cnt])) / cellsize;
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i][j + 1] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i][j + 1] + totalin[cnt])) / cellsize;
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i][j + 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i][j + 1]);
					}
				}
				else if (cnt == 4) { // SE
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i + 1][j + 1] + totalin[cnt]))) ?
							((DHM[i - 1][j - 1] + totalin[8]) - (DHM[i + 1][j + 1] + totalin[cnt])) / (cellsize * 2 * 1.4142135) :
							((DHM[i][j] + totalin[0]) - (DHM[i + 1][j + 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j + 1] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i + 1][j + 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j + 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j + 1]);
					}
					dist = cellsize*1.4142135;
				}
				else if (cnt == 5) { // S		     	    
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i + 1][j] + totalin[cnt]))) ?
							((DHM[i - 1][j] + totalin[1]) - (DHM[i + 1][j] + totalin[cnt])) / (cellsize * 2) :
							((DHM[i][j] + totalin[0]) - (DHM[i + 1][j] + totalin[cnt])) / cellsize;
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i + 1][j] + totalin[cnt])) / cellsize;
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j]);
					}
				}
				else if (cnt == 6) { // SW		     	    
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i + 1][j - 1] + totalin[cnt]))) ?
							((DHM[i - 1][j + 1] + totalin[2]) - (DHM[i + 1][j - 1] + totalin[cnt])) / (cellsize * 2 * 1.4142135) :
							((DHM[i][j] + totalin[0]) - (DHM[i + 1][j - 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j - 1] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i + 1][j - 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j - 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i + 1][j - 1]);
					}
					dist = cellsize*1.4142135;
				}
				else if (cnt == 7) { // W		     	    
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i][j - 1] + totalin[cnt]))) ?
							((DHM[i][j + 1] + totalin[3]) - (DHM[i][j - 1] + totalin[cnt])) / (cellsize * 2) :
							((DHM[i][j] + totalin[0]) - (DHM[i][j - 1] + totalin[cnt])) / cellsize;
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i][j - 1] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i][j - 1] + totalin[cnt])) / cellsize;
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i][j - 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i][j - 1]);
					}
				}
				else { // NW, cnt==8					     	    
					if (celltype == 2)
					{
						S = ((DHM[i][j] <= (DHM[i - 1][j - 1] + totalin[cnt]))) ?
							((DHM[i + 1][j + 1] + totalin[4]) - (DHM[i - 1][j - 1] + totalin[cnt])) / (cellsize * 2 * 1.4142135) :
							((DHM[i][j] + totalin[0]) - (DHM[i - 1][j - 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff2 = elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j - 1] + totalin[cnt]);
					}
					else
					{
						S = ((DHM[i][j] + totalin[0]) - (DHM[i - 1][j - 1] + totalin[cnt])) / (cellsize*1.4142135);
						elevdiff = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j - 1] + totalin[cnt]);
						elevdiff2 = (DHM[i][j] + totalin[0]) - (DHM[i - 1][j - 1]);
					}
					dist = cellsize*1.4142135;
				}
				if (S>0)
				{
					S=sqrt(S);
				}
				else
				{
					S=0;
				}

				double v     = SRFCVEL[i][j][cnt];
				double v_old = v; // default
				//double Q   = cellsize*cellsize*WEIGHT[i][j][cnt]*HEAD[i][j]/(InV*60.0); // in m3/s
				double Q = cellsize*cellsize*WEIGHT[i][j][cnt] * HEAD[i][j] / dMaxSubTimeStep; // in m3/s

				Qmax = 2 * cellsize*cellsize*WEIGHT[i][j][cnt] * elevdiff / dMaxSubTimeStep;
//				Qmax = cellsize*cellsize*WEIGHT[i][j][cnt] * elevdiff2 / dMaxSubTimeStep;

				if (Q > Qmax) { Q = Qmax; }

//				double Q = cellsize*cellsize*WEIGHT[i][j][cnt] * (HEAD[i][j]-totalin[cnt]) / dMaxSubTimeStep; // in m3/s
				double RH = 0;
				double vtemp;
				double M;
				if (method!=5)
				{
					int iter=0;
					double dActConcFact=1; // will be set to dConcFact when shear stres tau is larger than tau_crit
//					while (((fabs(v_old - v)>dMaxDeltaVelocity) || (!iter)) && (iter < this->iMaxIterVelocity)) // 1 mm/s is iteration tolerance
					while (((std::abs(v_old - v)>dMaxDeltaVelocity) || (!iter)) && (iter < this->iMaxIterVelocity)) // 1 mm/s is iteration tolerance
					{
						// 1. store old v in vold for getting the difference later to stop the iteration
						vtemp=v_old;
						v_old=v;
						v=(vtemp+v)/2.0;
						// 1. calculate hydraulic radius as Q/v/cs=RH
						RH = Q/v/dist;
						// 2. calculate v with this new hydraulic radius
						M=GetRoughness(i,j,v,luse,SR_cultivation,RH,VCFlocal,S);
						v = M*pow((double)RH,(double)(0.66667))*S;
						iter++;
					}

					// calulating tau = rho * g * RH * S * (n_g/n)^(2/3) 
					// = 9815 Kg/m^3*m/s^-2 * RH [m] * S * (M/Kst)^(2/3) --> is already sqrt(S), so we have to write S*S
					double tau; // shear stress
					tau = 9815.0 * RH * S * S * pow(double(M/this->LUtable[luse].GetSR_K_ST()),double(1.5));
					// if tau > TAUCRIT then apply ConcentrationFactor and redo v-estimation
					if ((tau >= TAUCRIT[i][j])&&(dConcFact>1.0))
					{
						iter=0;
						dActConcFact=dConcFact; // will be set to dConcFact when shear stres tau is larger than tau_crit
						v_old = v; // default
						//while (((fabs(v_old - v)>dMaxDeltaVelocity) || (!iter)) && (iter < this->iMaxIterVelocity)) // 1 mm/s is iteration tolerance
						while (((std::abs(v_old - v)>dMaxDeltaVelocity) || (!iter)) && (iter < this->iMaxIterVelocity)) // 1 mm/s is iteration tolerance
						{
							// 1. store old v in vold for getting the difference later to stop the iteration
							vtemp=v_old;
							v_old=v;
							v=(vtemp+v)/2.0;
							// 1. calculate hydraulic radius as Q/v/cs=RH
							RH = Q/v/dist*dActConcFact;
							// 2. calculate v with this new hydraulic radius
							M=GetRoughness(i,j,v,luse,SR_cultivation,RH,VCFlocal,S);
							v = M*pow((double)RH,(double)(0.66667))*S;
							iter++;
						}
					}

					// limiting RH to the water level above ground -> could not be greater than precipitation plus inflow plus storage
//					if (RH / dActConcFact>(HEAD[i][j] - totalin[cnt]))
					if (RH / dActConcFact>(HEAD[i][j]))
					{
//						RH = (HEAD[i][j] - totalin[cnt])*dActConcFact;
						RH = HEAD[i][j]*dActConcFact;
						v_old = v;
						M=GetRoughness(i,j,v,luse,SR_cultivation,RH,VCFlocal,S);
						v = M*pow((double)RH,(double)0.66667)*S;
						int iter=0;
						// iteration because of roughness and thus v may vary when a fix RH is given
						//while (((fabs(v_old-v)>dMaxDeltaVelocity)||(!iter))&&(iter < this->iMaxIterVelocity))
						while (((std::abs(v_old - v)>dMaxDeltaVelocity) || (!iter)) && (iter < this->iMaxIterVelocity))
						{
							v_old=v;
							M=GetRoughness(i,j,v,luse,SR_cultivation,RH,VCFlocal,S);
							v = M*pow((double)RH,(double)(0.66667))*S;
							iter++;
						}
					}
					// make sure, the program is numerically stable, so runoff of less than 0.1mm will remain on the cell
					if (RH<0.0001)
					{
						v=0;
						RH=0;
					}
//					else if (v*global_mintime > cellsize/2.0)
//					{
//						v = cellsize / global_mintime / 2.0;
//					}
					SRFCVEL[i][j][cnt]=(double)v;
					FLUX[i][j][cnt]=(double)(v*dist*RH/dActConcFact);
					FLUX[i][j][0]+=FLUX[i][j][cnt];
				}
				else
				{
					// simple approach for method 5 (back water conditions with very short time steps)
					RH = HEAD[i][j];
					M=10;
					v = M*pow(RH,0.66667)*S;
					double maxv = ((cnt==2)||(cnt==4)||(cnt==6)||(cnt==8))?((cellsize*1.4142135)/(20*dMinSubTimeStep)):(cellsize/(20*dMinSubTimeStep));
					if (v>maxv) v=maxv;
					SRFCVEL[i][j][cnt]=v;
					FLUX[i][j][cnt]=(v*dist*RH);
					FLUX[i][j][0]+=FLUX[i][j][cnt];
				}
			}
			else
			{
				SRFCVEL[i][j][cnt]=0;
				FLUX[i][j][cnt]=0;
			}
			fpe_check(i,j,99,cnt,WEIGHT[i][j][cnt],SRFCVEL[i][j][cnt],FLUX[i][j][cnt],0,0);
		}
		FLUX1[i][j] = FLUX[i][j][0];
		// due to iteration inaccuracies we have to adjust fluxes to the maximum, which is the total input in HEAD
//		if (FLUX[i][j][0]>0)
		if ((method!=5)&&(FLUX[i][j][0]>0))
		{
			//			double Q       = cellsize*cellsize*HEAD[i][j]/(InV*60.0); // in m3/s
			double Q = cellsize*cellsize*HEAD[i][j]/dMaxSubTimeStep; // in m3/s
			if (Q > Qmax) { Q = Qmax; }

			double redfact = Q / FLUX[i][j][0];

			if (redfact<1.0)
			{
				FLUX[i][j][0] = 0;
				for (int cnt=1;cnt<=8;cnt++)
				{
					if (FLUX[i][j][cnt]>0)
					{
						FLUX[i][j][cnt]*=redfact;
						FLUX[i][j][0] += FLUX[i][j][cnt];
					}
				}
				FLUX1[i][j] = FLUX[i][j][0];
//				redfact = Q / FLUX[i][j][0];
			}
		}
	}
	else
	{
		if (WEIGHT[i][j][0]<8)
		{
			for (int cnt = 0; cnt <= 8; cnt++)
			{
				SRFCVEL[i][j][cnt]=0;
				FLUX[i][j][cnt]=0;
			}
		}
	}
	FLUX1[i][j] = FLUX[i][j][0] = FLUX[i][j][(int)SRFCVEL[i][j][0]];

}

double surfaceroutingclass::MinimumFlowTime(int step)
{
	double mintime=1e36;
	double flowtime;
	for (int i=mpistartrow;i<=mpiendrow;i++) /* Reihen des Grids */
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				for (int cnt=1;cnt<=8;cnt++)
				{
					if (method<6)
					{
						if ((WEIGHT[i][j][cnt]>1e-4) && (SRFCVEL[i][j][cnt] > 1e-4))
						{
							flowtime = cellsize / SRFCVEL[i][j][cnt];
							if (flowtime < mintime) mintime = flowtime;
						}
					}
					else
					{
						if (fabs(SRFCVEL[i][j][cnt]) > 1e-6)
						{
							flowtime = cellsize / SRFCVEL[i][j][cnt];
							if (fabs(flowtime) < mintime) mintime = fabs(flowtime);
						}
					}
				}
			}
		}
	}
	return mintime;
}

// this function balances all inflows into a cell in SRQIN[i][j] and all outflows out of the cell 
// in m^3 for the actual sub time step
void surfaceroutingclass::CalculateWaterFlows(int i,int j,double mintime)
{
	SRQOUT[i][j]=0.0;
	SRQIN[i][j]=0.0;
	// get inflows from all eight neighbors or totalizing outflow to all possible eight neighbors
	// result for SRQIN is in m^3
	for (int z=1;z<=tracer_count;z++)
	{
		tracer_module[z].concflux[i][j]=0; // here, we are summing up the tracer masses for later mixing with actual content and inflowing/outflowing amounts
	}
	if  (i>1) 
	//		if (zn[i-1][j  ]==zn[i][j])  
	if (zn[i-1][j  ]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i-1][j  ][5];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i-1][j  ][5]*tracer_module[z].concstor[i-1][j];
		}
	}
	if ((i>1)    &&(j<ncols)) 
	//		if (zn[i-1][j+1]==zn[i][j])  
	if (zn[i-1][j+1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i-1][j+1][6];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i-1][j+1][6]*tracer_module[z].concstor[i-1][j+1];
		}
	}
	if             (j<ncols)  
	//		if (zn[i  ][j+1]==zn[i][j])  
	if (zn[i  ][j+1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i  ][j+1][7];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i][j+1][7]*tracer_module[z].concstor[i][j+1];
		}
	}
	if ((i<nrows)&&(j<ncols)) 
	//		if (zn[i+1][j+1]==zn[i][j])  
	if (zn[i+1][j+1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i+1][j+1][8];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i+1][j+1][8]*tracer_module[z].concstor[i+1][j+1];
		}
	}
	if  (i<nrows)             
	//		if (zn[i+1][j  ]==zn[i][j])  
	if (zn[i+1][j  ]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i+1][j  ][1];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i+1][j  ][1]*tracer_module[z].concstor[i+1][j];
		}
	}
	if ((i<nrows)&&(j>1))     
	//		if (zn[i+1][j-1]==zn[i][j])  
	if (zn[i+1][j-1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i+1][j-1][2];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i+1][j-1][2]*tracer_module[z].concstor[i+1][j-1];
		}
	}
	if                 (j>1)  
	//		if (zn[i  ][j-1]==zn[i][j])  
	if (zn[i  ][j-1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i  ][j-1][3];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i][j-1][3]*tracer_module[z].concstor[i][j-1];
		}
	}
	if ((i>1)    &&(j>1))     
	//		if (zn[i-1][j-1]==zn[i][j])  
	if (zn[i-1][j-1]!=nodata)  
	{
		SRQIN[i][j]+=FLUX[i-1][j-1][4];
		for (int z=1;z<=tracer_count;z++)
		{
			tracer_module[z].concflux[i][j]+=FLUX[i-1][j-1][4]*tracer_module[z].concstor[i-1][j-1];
		}
	}
	// tracer_module[z].concflux[i][j] contains mass fluxes now (since we added flux in m3/s times Kg/Kg = 1000*Kg/s
	for (int z=1;z<=tracer_count;z++)
	{
		tracer_module[z].concflux[i][j]*=mintime*1000.0; // now in kg as absolute value
	}
	SRQIN[i][j]*=mintime; // conversion from m^3/s into m^3
	// SRQIN[i][j] contains the total input from all neighboring cells in m^3 for the actual sub time step

	// now SRQOUT --> sum of outflows into all possible cells in m^3 for the actual sub time step
	for (int cnt=1;cnt<=8;cnt++)
	{
		if (FLUX[i][j][cnt]>0)
		{
			SRQOUT[i][j]+=FLUX[i][j][cnt]; // result is in m^3/s
		}
	}
	SRQOUT[i][j]*=mintime; // conversion from m^3/s into m^3
}

void surfaceroutingclass::effective_roughness()
{
#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			// Roughness values for actual cell R0, cell to the right R1 and South R2 and temporary sum TS for geometric average
			double R0 = 0;
			double R1 = 0;
			double R2 = 0;
			double RS = 0; 
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
			{
				
				double const minManningScaling = 0.02; // water depth below this value: no reduction for transmissivity. Set to 0.01 m my help

				if (RIVDEP[i][j] > 0 && RIVWIT[i][j] > 0) 
					R0 = MANNING[i][j]; 
				else 
//					R0 = MANNING[i][j];
					R0 = (HEAD[i][j] - DHM[i][j])> minManningScaling ? MANNING[i][j] * (HEAD[i][j] - DHM[i][j]) : MANNING[i][j] * minManningScaling;

				if ((j < ncols) && (zn[i][j + 1] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i][j + 1])))
				{
					if (RIVDEP[i][j+1] > 0 && RIVWIT[i][j+1] > 0) 
						R1 = MANNING[i][j+1];
					else 
//						R1 = MANNING[i][j + 1];
//						R1 = MANNING[i][j + 1] * (HEAD[i][j + 1] - DHM[i][j + 1]);
						R1 = (HEAD[i][j+1] - DHM[i][j+1]) > minManningScaling ? MANNING[i][j+1] * (HEAD[i][j+1] - DHM[i][j+1]) : MANNING[i][j+1] * minManningScaling;
				}
				else 
				{
					R1 = 0;
				}

				if ((i<nrows) && (zn[i + 1][j] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i + 1][j])))
				{
					if (RIVDEP[i + 1][j] > 0 && RIVWIT[i + 1][j] > 0)
						R1 = MANNING[i + 1][j];
					else
//						R2 = MANNING[i + 1][j];
//						R2 = MANNING[i + 1][j] * (HEAD[i + 1][j] - DHM[i + 1][j]);
	    				R2 = (HEAD[i+1][j] - DHM[i+1][j]) > minManningScaling ? MANNING[i+1][j] * (HEAD[i+1][j] - DHM[i+1][j]) : MANNING[i+1][j] * minManningScaling;
				}
				else
				{
					R2 = 0;
				}

				/* if the cell is not a no-flow-boundary (marked by MANNING==0),
				it should be possible for water to enter the cell even if the cell is dry.
				in such a case, the transmissivities are build only from the neigboring
				cell from which water com3es into the actual dry cell.
				This should be true also for groundwater tables below the bottom of an aquifer
				(then is T<0 because of THICK+GWSTAND <0) */
				// making sure that a dry cell can be filled up but an empty cell cannot drain more water: actual cell
				if (R0 <= 0.0)
				{
					if ((R1>0) && (HEAD[i][j + 1]>HEAD[i][j])) // but only if the gradient points towards the actual cell; river cells are no source cells, so that is not important here
					{
						if (MANNING[i][j]>0.0)
						{
							R0 = R1;
						}
						else
						{
							R0 = 0.0;
						}
					}
					else
					{
						R0 = 0.0;
					}
				}
				
				// making sure that a dry cell can be filled up but an empty cell cannot drain more water: eastern cell
				// for river cells in the neighbourhood, R1 is always > 0 (if MANNING was set to > 0)
				if ((R1 <= 0.0) && (j<ncols) && (zn[i][j + 1] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i][j + 1])))
				{
					if ((R0>0) && (HEAD[i][j]>HEAD[i][j + 1]))
					{
						if (MANNING[i][j + 1]>0.0)
						{
							R1 = R0;
						}
						else
						{
							R1 = 0.0;
						}
					}
					else
					{
						R1 = 0.0;
					}
				}

				// making sure that a dry cell can be filled up but an empty cell cannot drain more water: southern cell
				// for river cells in the neighbourhood, R1 is always > 0 (if MANNING was set to > 0)
				if ((R2 <= 0.0) && (i<nrows) && (zn[i + 1][j] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i + 1][j])))
				{
					if ((R0>0) && (HEAD[i][j]>HEAD[i + 1][j]))
					{
						if (MANNING[i + 1][j]>0.0)
						{
							R2 = R0;
						}
						else
						{
							R2 = 0.0;
						}
					}
					else
					{
						R2 = 0.0;
					}
				}

				RS = R0 + R1;
				ROUGHX[i][j] = RS>0.0 ? 2 * R0*R1 / RS : 0.0;

				RS = R0 + R2;
				ROUGHY[i][j] = RS>0.0 ? 2 * R0*R2 / RS : 0.0;

//				if (R0 < R1) ROUGHX[i][j] = R0; else ROUGHX[i][j] = R1;
//				if (R0 < R2) ROUGHY[i][j] = R0; else ROUGHY[i][j] = R2;

			}
			else
			{
				// glacierized cells should be handled as having frozen ground, so the transmissivity is 0, as should be set for nodata cells in zn, since the algorithm is taking into account cells from the nodata neighborhood as well
				// do NOT call this line for MPIREG[i][j]==mpi_rank only, since in gauss_seidel_iteration the values are assumed to have valid values even for cells outside the actual zones (shoudl then be 0)
				ROUGHX[i][j] = ROUGHY[i][j] = 0.0;
			}
		}
	}

}

double surfaceroutingclass::gauss_seidel_iteration(double DT,int iter)
{
	double maxerror = 0;
	//double DT = 60.0*InV / substepcnt; /* timestep in seconds */

#ifdef _OPENMP
#ifndef _OMP3
	int nthreads = omp_get_max_threads();
	double *maxerrorarray = new double[nthreads];
	for (int n = 0; n < nthreads; n++) maxerrorarray[n] = 0;

#pragma omp parallel for schedule(runtime) shared(maxerrorarray)
#else
#pragma omp parallel for schedule(runtime) reduction(max:maxerror)
#endif
#endif
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
#ifdef _OPENMP
		int th_id = omp_get_thread_num();
#endif
		double err, headold, const_boundary_q; /* for calculating iteration error */
		double h1, h2, h3, h4, /* heads  (1) back, (2) right, (3) front and (4) left cell if looking on the cell from a higer row number */
			T1, T2, T3, T4, /* Transmissivities for the flow from (1) back, (2) right, (3) front and (4) left cell to the centered cell */
			f1, f2, f3, f4; /* flows from (1) back to center, (2) center to right, (3) center to front, (4) left to center*/
	//		double f1a, f2a, f3a, f4a, f5, f6, f7, f8;
		double fsum, Tsum;
		double atemp, btemp, ctemp, storage_coefficient, dx, dy;
//		double alpha = 1.0; // fully implicit
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j]>nodata) && (MPIREG[i][j] == mpi_rank)) // glacierized cells should not be taken into account (even partly glacierized)
			{
				if((RIVWIT[i][j]==nodata) || (RIVDEP[i][j]==nodata))
				{
					bool noflux = false;
					if (((j > 1     ) ? (ROUGHX[i][j - 1] <= 0) : true) &&
						((j < ncols ) ? (ROUGHX[i][j    ] <= 0) : true) &&
						((i>1       ) ? (ROUGHY[i - 1][j] <= 0) : true) &&
						((i < nrows ) ? (ROUGHY[i    ][j] <= 0) : true))
					{
						noflux = true;
					}

					if (!noflux)
					{
						if (RelCellArea != NULL)
						{
							dx = CellSizeX[i][j];
							dy = CellSizeY[i][j];
						}
						else
						{
							dx = dy = cellsize;
						}
						headold = HEAD[i][j]; // old head is used to calculate changes (errors) between iteration steps 
						// heads for t' in centerd coordinates 
						// note that h1 and h4 are newly calculated values, already, if the row index is 2 or larger. This may be a problem when using many cores in openMP or when using many ranks in MPI
						// since a growing number of rows has then no precalculated values as starting points and get more and more instable
						// depending on alpha, it is an explicite (alpha = 0), an fully implicite (alpha = 1) or an inbetween method (0 < alpha < 1) 

						if ((i > 1) && (zn[i - 1][j] != nodata))
//						if (i > 1)
						{
							h1 = HEAD[i - 1][j];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i - 1][j] > 0 && RIVWIT[i - 1][j] > 0)
							{
								h1 = DHM[i][j];
							}
							else if (DHM[i - 1][j] > HEAD[i][j])
							{
								h1 = HEAD[i][j] + (HEAD[i - 1][j] - DHM[i - 1][j]);
							}
							else if (DHM[i][j] > HEAD[i - 1][j])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h1 = DHM[i][j];
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h1 = HEAD[i][j];
						}

						if ((j < ncols) && (zn[i][j + 1] != nodata))
//						if (j < ncols)
						{
							h2 = HEAD[i][j + 1];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i][j+1] > 0 && RIVWIT[i][j+1] > 0)
							{
								h2 = DHM[i][j];
							}
							else if (DHM[i][j + 1] > HEAD[i][j])
							{
								h2 = HEAD[i][j] + (HEAD[i][j + 1] - DHM[i][j + 1]);
							}
							else if (DHM[i][j] > HEAD[i][j + 1])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h2 = DHM[i][j];
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h2 = HEAD[i][j];
						}

						if ((i < nrows) && (zn[i + 1][j] != nodata))
//						if (i < nrows)
						{
							h3 = HEAD[i + 1][j];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i + 1][j] > 0 && RIVWIT[i + 1][j] > 0)
							{
								h3 = DHM[i][j];
							}
							else if (DHM[i + 1][j] > HEAD[i][j])
							{
								h3 = HEAD[i][j] + (HEAD[i + 1][j] - DHM[i + 1][j]);
							}
							else if (DHM[i][j] > HEAD[i + 1][j])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h3 = DHM[i][j];
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h3 = HEAD[i][j];
						}

						if ((j > 1) && (zn[i][j-1] != nodata))
//						if (j > 1)
						{
							h4 = HEAD[i][j-1];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i][j-1] > 0 && RIVWIT[i][j-1] > 0)
							{
								h4 = DHM[i][j];
							}
							else if (DHM[i][j - 1] > HEAD[i][j])
							{
								h4 = HEAD[i][j] + (HEAD[i][j - 1] - DHM[i][j - 1]);
							}
							else if (DHM[i][j] > HEAD[i][j-1])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h4 = DHM[i][j];
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h4 = HEAD[i][j];
						}
						
						/* transmissivities in centered coordinates */
						T1 = i > 1 ?      dx / dy*ROUGHY[i - 1][j] : ROUGHY[i][j];
						T2 = j <= ncols ? dy / dx*ROUGHX[i    ][j] : ROUGHX[i][j];
						T3 = i <= nrows ? dx / dy*ROUGHY[i    ][j] : ROUGHY[i][j];
						T4 = j > 1 ?      dy / dx*ROUGHX[i][j - 1] : ROUGHX[i][j];
						/* lokal fluxes in centered coordinates */

//						h1 += -1000.0;
//						h2 += -1000.0;
//						h3 += -1000.0;
//						h4 += -1000.0;


						f1 = T1*h1;
						f2 = T2*h2;
						f3 = T3*h3;
						f4 = T4*h4;

//						f5 = T1 * h3;
//						f6 = T2 * h4;
//						f7 = T3 * h1;
//						f8 = T4 * h2;

						Tsum = T1 + T3 + T2 + T4;
						
						fsum = f1 + f2 + f3 + f4;
						// this does not work - it was an experiment to use the partial derivative of dT(x)/dx*dh/dx and dT(y)/dy*dh/dy as additional terms. But it becomes highly unstable after a few iterations already
//						fsum = 2*f1 + 2*f2 + 2*f3 + 2*f4 - f5 - f6 - f7 - f8;

						/* if any layer is unconfined, the storage coefficient of the
						uppermost aquifer is used which is assumed to be the one for unconfined conditions */
						storage_coefficient = 1.0; // could later be extended by snow porosity -> will then be like a groundwater flux in snow
						//          storage_coefficient=S0[lr][i][j];
						const_boundary_q = 0; // no constant flux boundary conditions used here

//						atemp = storage_coefficient / DT*HEADold[i][j] + const_boundary_q;
//						atemp *= (dx*dy);
//						btemp = fsum;
//						ctemp = Tsum + dx*dy*storage_coefficient / DT;
//
//						HEAD[i][j] = (ctemp != 0) ? (atemp + btemp) / ctemp : DHM[i][j]; // if for any reason ctemp were 0, the division would result in a NaN or -1.#inf etc. We limit it in this case to the actual DHM value (assuming that no water is any longer on the surface available for routing)

						double hsum = h1 + h2 + h3 + h4;
						ctemp = (dx*dy*storage_coefficient) / (Tsum*DT);
						btemp = hsum/4.0;
						atemp = ctemp * HEADold[i][j] + const_boundary_q * dx*dy / Tsum;
//						new head is calculated after dumping
						HEAD[i][j] = (ctemp != -1.0) ? (atemp + btemp) / (ctemp + 1.0) : DHM[i][j]; // if for any reason ctemp were 0, the division would result in a NaN or -1.#inf etc. We limit it in this case to the actual DHM value (assuming that no water is any longer on the surface available for routing)



						DumptypeSR dump;
						if (this->dumpmode == true)
						{
#ifdef _READDUMPFILESR
							bool readdump = true;
#else
							bool readdump = false;
#endif

							if (readdump)
							{
								std::ifstream dumpstream;
								dumpstream.open("dump.dmp", std::ios::binary);
								dumpstream.read((char *)(&DT), sizeof(double));
								double temp[5];
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // dhm
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // roughx
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // roughy
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // rivdep
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // rivwit
								dumpstream.read((char *)(&temp[0]), sizeof(double) * 5); // head

								dumpstream.read((char *)(&dump.h[0]), sizeof(double) * 5);
								dumpstream.read((char *)(&dump.T[0]), sizeof(double) * 5);
								dumpstream.read((char *)(&dump.f[0]), sizeof(double) * 5);
								dumpstream.read((char *)(&dump.dx), sizeof(double));
								dumpstream.read((char *)(&dump.dy), sizeof(double));
								dumpstream.read((char *)(&dump.headold), sizeof(double));
								dumpstream.read((char *)(&dump.Tsum), sizeof(double));
								dumpstream.read((char *)(&dump.fsum), sizeof(double));
								dumpstream.read((char *)(&dump.atemp), sizeof(double));
								dumpstream.read((char *)(&dump.btemp), sizeof(double));
								dumpstream.read((char *)(&dump.ctemp), sizeof(double));
							}
							else
							{
								dump.DT = DT;
								dump.dhm[0] = DHM[i][j];
								dump.dhm[1] = (i > 1) ? DHM[i - 1][j] : -9999;
								dump.dhm[2] = (j < ncols) ? DHM[i][j + 1] : -9999;
								dump.dhm[3] = (i < nrows) ? DHM[i + 1][j] : -9999;
								dump.dhm[4] = (j > 1) ? DHM[i][j - 1] : -9999;

								dump.roughx[0] = ROUGHX[i][j];
								dump.roughx[1] = (i > 1) ? ROUGHX[i - 1][j] : -9999;
								dump.roughx[2] = (j < ncols) ? ROUGHX[i][j + 1] : -9999;
								dump.roughx[3] = (i < nrows) ? ROUGHX[i + 1][j] : -9999;
								dump.roughx[4] = (j > 1) ? ROUGHX[i][j - 1] : -9999;

								dump.roughx[0] = ROUGHY[i][j];
								dump.roughx[1] = (i > 1) ? ROUGHY[i - 1][j] : -9999;
								dump.roughx[2] = (j < ncols) ? ROUGHY[i][j + 1] : -9999;
								dump.roughx[3] = (i < nrows) ? ROUGHY[i + 1][j] : -9999;
								dump.roughx[4] = (j > 1) ? ROUGHY[i][j - 1] : -9999;

								dump.rivdep[0] = RIVDEP[i][j];
								dump.rivdep[1] = (i > 1) ? RIVDEP[i - 1][j] : -9999;
								dump.rivdep[2] = (j < ncols) ? RIVDEP[i][j + 1] : -9999;
								dump.rivdep[3] = (i < nrows) ? RIVDEP[i + 1][j] : -9999;
								dump.rivdep[4] = (j > 1) ? RIVDEP[i][j - 1] : -9999;

								dump.rivwit[0] = RIVWIT[i][j];
								dump.rivwit[1] = (i > 1) ? RIVWIT[i - 1][j] : -9999;
								dump.rivwit[2] = (j < ncols) ? RIVWIT[i][j + 1] : -9999;
								dump.rivwit[3] = (i < nrows) ? RIVWIT[i + 1][j] : -9999;
								dump.rivwit[4] = (j > 1) ? RIVWIT[i][j - 1] : -9999;

								dump.head[0] = HEAD[i][j];
								dump.head[1] = (i > 1) ? HEAD[i - 1][j] : -9999;
								dump.head[2] = (j < ncols) ? HEAD[i][j + 1] : -9999;
								dump.head[3] = (i < nrows) ? HEAD[i + 1][j] : -9999;
								dump.head[4] = (j > 1) ? HEAD[i][j - 1] : -9999;

								dump.h[1] = h1;
								dump.h[2] = h2;
								dump.h[3] = h3;
								dump.h[4] = h4;

								dump.T[1] = T1;
								dump.T[2] = T2;
								dump.T[3] = T3;
								dump.T[4] = T4;

								dump.f[1] = f1;
								dump.f[2] = f2;
								dump.f[3] = f3;
								dump.f[4] = f4;

								dump.dx = dx;
								dump.dy = dy;
								dump.headold = headold;
								dump.Tsum = Tsum;
								dump.fsum = fsum;
								dump.atemp = atemp;
								dump.atemp = btemp;
								dump.atemp = ctemp;
							}
						}

//						HEAD[i][j] = (ctemp != 0) ? (atemp + btemp) / ctemp : DHM[i][j]; // if for any reason ctemp were 0, the division would result in a NaN or -1.#inf etc. We limit it in this case to the actual DHM value (assuming that no water is any longer on the surface available for routing)
//						HEAD[i][j] = (ctemp != 0) ? (atemp + btemp) / (ctemp + 1.0) : DHM[i][j]; // if for any reason ctemp were 0, the division would result in a NaN or -1.#inf etc. We limit it in this case to the actual DHM value (assuming that no water is any longer on the surface available for routing)

						if (this->dumpmode == true)
						{
							// HEAD[i][j] = log(-1.0); // for testig if it reacts to nan which is generated by log(-1.0)
							// HEAD[i][j] /= 0.0; // for testig if it reacts to inf which is generated by a floating point division by zero
							if (std::isnan(HEAD[i][j]) || std::isinf(HEAD[i][j]))
							{
								if ((actDumpFiles < maxDumpFiles) || maxDumpFiles == 0)
								{
									cout << "NaN or Inf found for HEAD[" << i << "][" << j << "] (value: " << HEAD[i][j] << ") in rank " << mpi_rank << endl;
									cout << "writing dump file..." << endl;
									std::ofstream dumpstream;
									std::stringstream fname;
									fname << "surfrout_";
									fname << "dump_row_" << i << "_col_" << j << "_iteration_" << iter << "_" << jahr << setw(2) << setfill('0') << monat << setw(2) << setfill('0') << tag << setw(2) << setfill('0') << stunde << ".dmp";
									std::string sname(fname.str());

									dump.newhead = HEAD[i][j];
									dumpstream.open(sname.c_str(), std::ios::binary);
									dumpstream.write((char *)(&dump.DT),        sizeof(double) * 5);
									dumpstream.write((char *)(&dump.dhm[0]),    sizeof(double) * 5);
									dumpstream.write((char *)(&dump.roughx[0]), sizeof(double) * 5);
									dumpstream.write((char *)(&dump.roughy[0]), sizeof(double) * 5);
									dumpstream.write((char *)(&dump.rivdep[0]), sizeof(double) * 5);
									dumpstream.write((char *)(&dump.rivwit[0]), sizeof(double) * 5);
									dumpstream.write((char *)(&dump.head[0]),   sizeof(double) * 5);
									dumpstream.write((char *)(&dump.newhead),   sizeof(double));
									dumpstream.write((char *)(&dump.h[0]),      sizeof(double) * 5);
									dumpstream.write((char *)(&dump.T[0]),      sizeof(double) * 5);
									dumpstream.write((char *)(&dump.f[0]),      sizeof(double) * 5);
									dumpstream.write((char *)(&dump.dx),        sizeof(double));
									dumpstream.write((char *)(&dump.dy),        sizeof(double));
									dumpstream.write((char *)(&dump.headold),   sizeof(double));
									dumpstream.write((char *)(&dump.Tsum),      sizeof(double));
									dumpstream.write((char *)(&dump.fsum),      sizeof(double));
									dumpstream.write((char *)(&dump.atemp),     sizeof(double));
									dumpstream.write((char *)(&dump.btemp),     sizeof(double));
									dumpstream.write((char *)(&dump.ctemp),     sizeof(double));
									dumpstream.close();
								}
								actDumpFiles++;


								if (closeAfterDump) 	{
									(*WSstream) << "FATAL: \n\n -------- command line option -closeAfterDump active! ---------------------------------------------\n\n"
										<< " a dumpfile has been written and WaSiM aborts now.\n"
										<< "-------------------------------------------------------------\n\n" << "[end-entry]";

									throw(199);
								}
							}
						}

						err = fabs(HEAD[i][j] - headold);
						/* due to the much smaller storage coefficient in confined conditions,
						it is possible to use larger errors (since the same error of heads between
						succesive iteration steps means much less difference in flux-balance than
						for unconfined conditions) */
#ifndef _OMP3
						//#pragma omp critical
#ifdef _OPENMP
						if (err > maxerrorarray[th_id])
						{
							maxerrorarray[th_id] = err;
						}
#else
						if (err > maxerror)
						{
							maxerror = err;
						}
#endif
#else
						if (err > maxerror)
						{
							maxerror = err;
						}
#endif

//#ifndef _OMP3
//#pragma omp critical
//						if (err>maxerror)
//						{
//							maxerror = err;
//						}
//#else
//						if (err>maxerror)
//						{
//							maxerror = err;
//						}
//#endif
					}
					else
					{
						// no change here for now
						int k = 1;
					}
				}
				else
				{
					// if a constant head consition is true for this cell, then all HEADs should be set to this constant head
					// in order to avoid oscillation indications
					if ((RIVWIT[i][j] > 0) && (RIVDEP[i][j] > 0))
					{
						HEAD[i][j] = HEADold[i][j] = DHM[i][j];
					}
				}
			}
		}
	}

#ifdef _OPENMP
#ifndef _OMP3
	for (int n = 0; n < nthreads; n++)
	{
		if (maxerrorarray[n] > maxerror)
		{
			maxerror = maxerrorarray[n];
		}
	}
	delete[] maxerrorarray;
#endif
#endif

	return maxerror;
}

void surfaceroutingclass::GetFlowVelocities_AfterGaussSeidel()
{

#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		double h, h1, h2, h3, h4, b0, b1, b2, b3, b4; // heads  (1) back, (2) right, (3) front and (4) left cell if looking on the cell from a higer row number 
		double T1, T2, T3, T4; // Transmissivities for the flow from (1) back, (2) right, (3) front and (4) left cell to the centered cell 
		double f1, f2, f3, f4; // fluxes from the surrounding cells from (1) back, (2) right, (3) front and (4) left cell to the centered cell
		double storage_coefficient, dx, dy;

		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j]>nodata) && (MPIREG[i][j] == mpi_rank)) // glacierized cells should not be taken into account (even partly glacierized)
			{
				if ((RIVWIT[i][j] == nodata) || (RIVDEP[i][j] == nodata))
				{
					bool noflux = false;
					if (((j > 1) ? (ROUGHX[i][j - 1] <= 0) : true) &&
						((j < ncols) ? (ROUGHX[i][j] <= 0) : true) &&
						((i>1) ? (ROUGHY[i - 1][j] <= 0) : true) &&
						((i < nrows) ? (ROUGHY[i][j] <= 0) : true))
					{
						noflux = true;
					}

					if (!noflux)
					{
						if (RelCellArea != NULL)
						{
							dx = CellSizeX[i][j];
							dy = CellSizeY[i][j];
						}
						else
						{
							dx = dy = cellsize;
						}
						h = HEAD[i][j];
						b0 = h-DHM[i][j];

						if (i > 1)
						{
							h1 = HEAD[i - 1][j];
							b1 = HEAD[i - 1][j]-DHM[i - 1][j];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i - 1][j] > 0 && RIVWIT[i - 1][j] > 0)
							{
								h1 = DHM[i][j];
								b1 = 0;
							}
							else if (DHM[i - 1][j] > HEAD[i][j])
							{
								h1 = HEAD[i][j] + (HEAD[i - 1][j] - DHM[i - 1][j]);
							}
							else if (DHM[i][j] > HEAD[i - 1][j])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h1 = DHM[i][j];
								b1 = 0;
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h1 = HEAD[i][j];
							b1 = b0;
						}

						if (j < ncols)
						{
							h2 = HEAD[i][j + 1];
							b2 = HEAD[i][j + 1] - DHM[i][j+1];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i][j + 1] > 0 && RIVWIT[i][j + 1] > 0)
							{
								h2 = DHM[i][j];
								b2 = 0;
							}
							else if (DHM[i][j + 1] > HEAD[i][j])
							{
								h2 = HEAD[i][j] + (HEAD[i][j + 1] - DHM[i][j + 1]);
								b2 = HEAD[i][j + 1] - DHM[i][j + 1];
							}
							else if (DHM[i][j] > HEAD[i][j + 1])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h2 = DHM[i][j];
								b2 = 0;
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h2 = HEAD[i][j];
							b2 = b0;
						}

						if (i < nrows)
						{
							h3 = HEAD[i + 1][j];
							b3 = HEAD[i + 1][j] - DHM[i + 1][j];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i + 1][j] > 0 && RIVWIT[i + 1][j] > 0)
							{
								h3 = DHM[i][j];
								b3 = 0;
							}
							else if (DHM[i + 1][j] > HEAD[i][j])
							{
								h3 = HEAD[i][j] + (HEAD[i + 1][j] - DHM[i + 1][j]);
								b3 = HEAD[i + 1][j] - DHM[i + 1][j];
							}
							else if (DHM[i][j] > HEAD[i + 1][j])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h3 = DHM[i][j];
								b3 = 0;
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h3 = HEAD[i][j];
							b3 = b0;
						}

						if (j > 1)
						{
							h4 = HEAD[i][j - 1];
							b4 = HEAD[i][j - 1] - DHM[i][j - 1];
							// check if there is a hydraulic connection to the neighboring cell at all, i.e. if the DHM of the neighbouring cell is already above the water table ofthe actual cell, the head of the neighboring cell 
							// must be de-coupled from the actual cell by reducing it for this cell to a no flow boundary + the water table that is on the neighbouring cell
							if (RIVDEP[i][j - 1] > 0 && RIVWIT[i][j - 1] > 0)
							{
								h4 = DHM[i][j];
								b4 = 0;
							}
							else if (DHM[i][j - 1] > HEAD[i][j])
							{
								h4 = HEAD[i][j] + (HEAD[i][j - 1] - DHM[i][j - 1]);
								b4 = HEAD[i][j - 1] - DHM[i][j - 1];
							}
							else if (DHM[i][j] > HEAD[i][j - 1])
							{ // other way round: when the actual cells DHM is already above the head of the other cell, decouple that cell by setting h1 to the actual cells DHM, making sure that the actual cells head can't drop below it's DHM value
								h4 = DHM[i][j];
								b4 = 0;
							}
						}
						else
						{
							// no-flow boundary is achieved by setting the head to the same value as for the actual cell
							h4 = HEAD[i][j];
							b4 = b0;
						}

						// transmissivities in centered coordinates, measure is m^2/s because the roughness in m/s and the flow thickness in m is included already
						T1 = i > 1 ? dx / dy*ROUGHY[i - 1][j] : ROUGHY[i][j];
						T2 = j <= ncols ? dy / dx*ROUGHX[i][j] : ROUGHX[i][j];
						T3 = i <= nrows ? dx / dy*ROUGHY[i][j] : ROUGHY[i][j];
						T4 = j > 1 ? dy / dx*ROUGHX[i][j - 1] : ROUGHX[i][j];

						// all fluxes in positive formulation into the direction of flow
						// measures in m^3/s
						f1 = T1 * dx * (h1 - h) / dy; // flows from above/back to the central cell dx is the width of the flow to convert the flow to m3/s
						f2 = T2 * dy * (h - h2) / dx; // flows from the central cell to the right cell
						f3 = T3 * dx * (h - h3) / dy; // flows from the central cell to the front/below cell
						f4 = T4 * dy * (h4 - h) / dx; // flows from the left cell to the center cell
						
						// all other fluxes will be given in positive formulation regarding to their flow direction
						FLUX[i][j][1] = f1;
						FLUX[i][j][2] = 0;
						FLUX[i][j][3] = f2;
						FLUX[i][j][4] = 0;
						FLUX[i][j][5] = f3;
						FLUX[i][j][6] = 0;
						FLUX[i][j][7] = f4;
						FLUX[i][j][8] = 0;

						if ((T1 != 0) && (i>1))
						{
							double b = (b0 + b1)/2.0;
							if (fabs(b) > 1e-6)
							{
								SRFCVEL[i][j][1] = FLUX[i][j][1] / (dx*b);
							}
							else
							{
								SRFCVEL[i][j][1] = 0;
							}
						}
						else
						{
							SRFCVEL[i][j][1] = 0;
						}

						SRFCVEL[i][j][2] = 0;

						if ((T2 != 0) && (j<ncols))
						{
							double b = (b0 + b2) / 2.0;
							if (fabs(b) > 1e-6)
							{
								SRFCVEL[i][j][3] = FLUX[i][j][3] / (dy*b);
							}
							else
							{
								SRFCVEL[i][j][3] = 0;
							}
						}
						else
						{
							SRFCVEL[i][j][3] = 0;
						}

						SRFCVEL[i][j][4] = 0;

						if ((T3 != 0) && (i<nrows))
						{
							double b = (b0 + b3) / 2.0;
							if (fabs(b) > 1e-6)
							{
								SRFCVEL[i][j][5] = FLUX[i][j][5] / (dx*b);
							}
							else
							{
								SRFCVEL[i][j][5] = 0;
							}
						}
						else
						{
							SRFCVEL[i][j][5] = 0;
						}

						SRFCVEL[i][j][6] = 0;

						if ((T4 != 0) && (j>1))
						{
							double b = (b0 + b4) / 2.0;
							if (fabs(b) > 1e-6)
							{
								SRFCVEL[i][j][7] = FLUX[i][j][7] / (dy*b);
							}
							else
							{
								SRFCVEL[i][j][7] = 0;
							}
						}
						else
						{
							SRFCVEL[i][j][7] = 0;
						}

						SRFCVEL[i][j][8] = 0;

						// here, for the balance we need to count outflows negative, i.e. fluxes to the right und bottom cell f2 and f3. 
						// The other fluxes are inflows from above (f1) and ritgh (f4) and count positive as usual
						// measure in m3/s
//						FLUX[i][j][0] = FLUX[i][j][1] - FLUX[i][j][3] - FLUX[i][j][5] + FLUX[i][j][7];
						double fx = (fabs(FLUX[i][j][1]) + fabs(FLUX[i][j][5])) / 2.0;
						double fy = (fabs(FLUX[i][j][3]) + fabs(FLUX[i][j][7])) / 2.0;
						double temp = fx*fx + fy*fy;
						FLUX1[i][j] = FLUX[i][j][0] = (temp > 0) ? sqrt(temp) : 0;

						fx = (fabs(SRFCVEL[i][j][1]) + fabs(SRFCVEL[i][j][5])) / 2.0;
						fy = (fabs(SRFCVEL[i][j][3]) + fabs(SRFCVEL[i][j][7])) / 2.0;
						temp = fx*fx + fy*fy;
						SRFCVEL1[i][j] = SRFCVEL[i][j][0] = (temp > 0) ? sqrt(temp) : 0;


					}
					else
					{
						for (int cnt = 0; cnt <= 8; cnt++)
						{
							SRFCVEL[i][j][cnt]=FLUX[i][j][cnt] = 0;
						}
						// no change here for now
					}
				}
			}
		}
	}

// synchronizing fluxes here!
#ifdef _MPI
	MPISynchronizeStackHalo(fluxindex, 10);
	MPISynchronizeStackHalo(srfcvelindex, 10);
#endif

// now filling the fluxes for river cells which were not known above and must no be constructed from the fluxes of the neighbouring cells
// the reason for this to be an extra loop is that the neighbouring fluxes were not all known before finishing the above loop
#pragma omp parallel for schedule(runtime) 
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j] > nodata) && (MPIREG[i][j] == mpi_rank)) // glacierized cells should not be taken into account (even partly glacierized)
			{
				if ((RIVWIT[i][j] != nodata) && (RIVDEP[i][j] != nodata))
				{
					if ((i > 1) && (zn[i-1][j] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i-1][j])))
					{
						FLUX[i][j][1] = (FLUX[i - 1][j][5] != nodata) ? FLUX[i - 1][j][5] : 0;
						SRFCVEL[i][j][1] = SRFCVEL[i - 1][j][5];
					}
					else
					{
						FLUX[i][j][1] = 0;
						SRFCVEL[i][j][1] = 0;
					}
					
					FLUX[i][j][2] = 0;
					SRFCVEL[i][j][2] = 0;

					if ((j < ncols) && (zn[i][j+1] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i][j+1])))
					{
						FLUX[i][j][3] = (FLUX[i][j + 1][7] != nodata) ? FLUX[i][j + 1][7]:0;
						SRFCVEL[i][j][3] = SRFCVEL[i][j + 1][7];
					}
					else
					{
						FLUX[i][j][3] = 0;
						SRFCVEL[i][j][3] = 0;
					}

					FLUX[i][j][4] = 0;
					SRFCVEL[i][j][4] = 0;

					if ((i<nrows) && (zn[i+1][j] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i+1][j])))
					{
						FLUX[i][j][5] = (FLUX[i + 1][j][1] != nodata) ? FLUX[i + 1][j][1]:0;
						SRFCVEL[i][j][5] = SRFCVEL[i + 1][j][1];
					}
					else
					{
						FLUX[i][j][5] = 0;
						SRFCVEL[i][j][5] = 0;
					}

					FLUX[i][j][6] = 0;
					SRFCVEL[i][j][6] = 0;

					if ((j>1) && (zn[i][j-1] != nodata) && !((this->do_it == 3) && (zn[i][j] != zn[i][j-1])))
					{
						FLUX[i][j][7] = (FLUX[i][j - 1][3] != nodata) ? FLUX[i][j - 1][3]:0;
						SRFCVEL[i][j][7] = SRFCVEL[i][j - 1][3];
					}
					else
					{
						FLUX[i][j][7] = 0;
						SRFCVEL[i][j][7] = 0;
					}

					FLUX[i][j][8] = 0;
					SRFCVEL[i][j][8] = 0;

					// sum of ll fluxes: here, outflows count negative, i.e. fluxes to the right (index 3) and down (index 5)
					double fx = (fabs(FLUX[i][j][1]) + fabs(FLUX[i][j][5])) ;
					double fy = (fabs(FLUX[i][j][3]) + fabs(FLUX[i][j][7])) ;
					double temp = fx*fx + fy*fy;
					FLUX1[i][j] = FLUX[i][j][0] = (temp > 0) ? sqrt(temp) : 0;

					fx = (fabs(SRFCVEL[i][j][1]) + fabs(SRFCVEL[i][j][5])) ;
					fy = (fabs(SRFCVEL[i][j][3]) + fabs(SRFCVEL[i][j][7])) ;
					temp = fx*fx + fy*fy;
					SRFCVEL1[i][j] = SRFCVEL[i][j][0] = (temp > 0) ? sqrt(temp) : 0;


				}
			}
		}
	}

}

void surfaceroutingclass::pre_run(long me,long se, long tg, long mt, long jr, long atcsubstep)
{
	// print a message to the output stream and register the actual date and time
	teilmodellklasse::run(me,se, tg, mt, jr,NULL,NULL);

	if (do_it == 0) return;

	//	if (mpi_rank>0) return; // this module is currently not able to run on a cluster. Only the instance running mpi rank 0 will do surface routing in order to keep the model running at all

	// initializing statistics register
	QDSRstat->init();
	QISRstat->init();
	QBSRstat->init();
	QGSRstat->init();
	SFSTSRstat->init();

	for (int z = 1; z <= tracer_count; z++)
	{
		tracer_module[z].stat->init();
	}

	tempsum = 0;
	// shifting registers for flow travel time zones -> this is part of the sample code (copy from unsatzonmodellklasse)

	for (int z = 0; z <= tracer_count; z++)
	{
		for (int i = 1; i <= zonmax; i++)
		{
			for (int j = 2; j <= (int)maxzeit[i]; j++)
			{
				zoneFZT_QD[z][i][j - 1].n = zoneFZT_QD[z][i][j].n;
				zoneFZT_QD[z][i][j - 1].RelCellSum = zoneFZT_QD[z][i][j].RelCellSum;
				zoneFZT_QD[z][i][j - 1].sx = zoneFZT_QD[z][i][j].sx;
				zoneFZT_QI[z][i][j - 1].n = zoneFZT_QI[z][i][j].n;
				zoneFZT_QI[z][i][j - 1].RelCellSum = zoneFZT_QI[z][i][j].RelCellSum;
				zoneFZT_QI[z][i][j - 1].sx = zoneFZT_QI[z][i][j].sx;
			}
			/* only the most distant zone is set to 0 (the other zones contain the surface runoff from the last intervals */
			zoneFZT_QD[z][i][(int)maxzeit[i]].sx = 0;
			zoneFZT_QD[z][i][(int)maxzeit[i]].n = 0;
			zoneFZT_QD[z][i][(int)maxzeit[i]].RelCellSum = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].sx = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].n = 0;
			zoneFZT_QI[z][i][(int)maxzeit[i]].RelCellSum = 0;
		}
	}
}

void surfaceroutingclass::run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
	
	if (old_year == jr && old_month == mt && old_day == tg && old_hour == se && old_minute == me)
	{
		moduleSubstepCounter++;
	}
	else
	{
		moduleSubstepCounter = 1;
		old_year = jr;
		old_month = mt;
		old_day = tg;
		old_hour = se;
		old_minute = me;
	}

	for (int i = 0; i <= zonmax; i++)
	{
		zoneHEADStart[i].n = 0;
		zoneHEADStart[i].sx = 0;
		zoneHEADEnd[i].n = 0;
		zoneHEADEnd[i].sx = 0;
		zoneQDFromFlux[i].n = 0;
		zoneQDFromFlux[i].sx = 0;
	}
	
	//	if (mpi_rank>0) return; // this module is currently not able to run on a cluster. Only the instance running mpi rank 0 will do surface routing in order to keep the model running at all
#pragma omp parallel for schedule(runtime) 
	// first estimate for HEAD: QDgrid (in unsatzonmodel, SSSU was already regarded as input, so QDgrid will contain all of the SSSU-Water again...
	for (int i=mpistartrow;i<=mpiendrow;i++) /* Reihen des Grids */
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				HEAD[i][j]=(QDgrid[i][j])/1000.0; // input in m, because all slope and velocity calculations are done with SI measures
				if (method == 6)
				{
					HEAD[i][j] += DHM[i][j]; // method 6 uses the absolute head while in methods 1 to 5 the local slopes are computed seperately
					HEADold[i][j] = HEAD[i][j];
				}
				SSSU[i][j]=QDgrid[i][j];
				for (int z=1;z<=tracer_count;z++)
				{
					// since we share the storage with silting up module, no decay must be calculated here. Exception: if silting up is switched off
					tracer_module[z].concstor[i][j] = CONC_QD[z][i][j]; // only concentrations are stored here, conversion into masses will follow when computing the fluxes
					if (SiltingUp==false)
					{
						if (tracer_info[z].radioactive)
						tracer_module[z].concstor[i][j]*=tracer_info[z].decfact;
					}
					CONC_SSSU[z][i][j]=0;
					CONC_QD[z][i][j]=0;
					// initialize tracer concentration in this module -> it will be given to silting up in CONC_SSSU later
				}
				QDgrid[i][j]=0; // will be filled again for river cells only within surface routing
			}
		}
	}

	int step = 0;
	int NextMajorStep = 0;
	int numsteps = 0;
	double subtimestepsum = 0;
	double mintime = 0;
	global_mintime = dMinSubTimeStep;
	double actualStepLength = global_mintime;

  //#pragma omp parallel for schedule(runtime) shared(step)
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
			{
				// for scaling of the QDgrid values with the change in overall HEAD (which should be identical, since suring routing there are no ther sources or sinks)
				zoneHEADStart[zarray[(int)zn[i][j]]].n++;
				zoneHEADStart[zarray[(int)zn[i][j]]].sx += (HEAD[i][j] - DHM[i][j]);
			}
		}
	}

	while (subtimestepsum < ((InV / substepcnt*60.0) - (dMinSubTimeStep / 2.0)))
	{

		int iter = 0;
		int maxiter = this->iMaxIterVelocity;
		double maxerror = 0;
		double maxerrorf = this->dMaxDeltaVelocity;
		double global_maxerror = 0;

		if (method < 6)
		{
#ifdef _MPI
			MPISynchronizeGridHalo(headindex); // all instances need the heighboring informations about HEAD
			//MPISynchronizeGridHaloJustRecv(headindex); //LAURA
			//MPISynchronizeGridHaloJustSend(headindex); //LAURA
			//MPISynchronizeGridHaloJustWait(headindex); //LAURA

			MPISynchronizeGridHalo(headoldindex); // all instances need the heighboring informations about HEAD
			//MPISynchronizeGridHaloJustRecv(headoldindex); //LAURA
			//MPISynchronizeGridHaloJustSend(headoldindex); //LAURA
			//MPISynchronizeGridHaloJustWait(headoldindex); //LAURA

			MPISynchronizeStackHalo(fluxindex, 10);

#endif
			// do not use any common variable in the openmp-loop for writing! multiple threads are accessing the same variable at the same time!
			// all variables with write access must be listed in the private-part of the openmp-clause. If the variables are even of interest 
			// for the results outside the loop, then the accumulation must be done using the suitable openmp-directives (or do it manually as 
			// we do it with the zonal statistics by doing such accumulation in a seperate loop after the central processing.
			// note: the input for routing is found in QDgrid, but for a substep it will be in HEAD including inflows from other cells
#pragma omp parallel for schedule(runtime) shared(step,NextMajorStep)
			// iterate for velocities (including inflows!)
			for (int i = mpistartrow; i <= mpiendrow; i++)  /* Reihen des Grids */
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) /* Spalten des Grids */
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
					{

						// do routing here in any way you want...
						// routing of water into other cells: use SSSU as storage grid and QDgrid as input grid
						// routing of the respective tracers: using CONC_QD with QDgrid as volumes as iput and CONC_SSSU[z] as storage

						// get thy type of the cell
						// 1 = at the border of the catchment, no river on the cell, multiple flow only allowed to cells within the subcatchment (up to three)
						// 2 = on a river cell -> immediate runoff
						// 3 = in the middle of nowhere = inmidst of the catchment but not on a river, multiple flow path possible (up to three)
						// 4 = not on a river cell but at the border of the catchment without slope into the catchment, the neighbor cell outside the catchment is a river -> will be handled like runoff into a river
						int celltype;
						double slope, aspect;

						GetSlopeAspectAndCelltype(i, j, slope, aspect, celltype, step); // slope and 

						//fpe_check(i,j,celltype,step,slope,aspect,HEAD[i][j],QDgrid[i][j],SSSU[i][j]);

						GetMultipleFlowPath(i, j, slope, aspect, celltype, step);
						// now, in WEIGHT[i][j][1..8] the weightings for the surface runoff are stored.
						// calculate velocities in each direction -> iteration!
						GetFlowVelocities(i, j, step, NextMajorStep, celltype); // and now, the matching flow velocities are stored in SRFCVEL[i][j][1..8]  as well as the fluxes are stored in FLUX[i][j][1..8] in m/s

					}
				}
			}
		}
		else
		{ // Gauss Seidel fully implicit iteration (stable, but many iterations needed); May be extended later to Crank-Nicholson, if that turns out to be  equally stable but faster 
			effective_roughness();
#ifdef _MPI
			MPISynchronizeGridHalo(roughxindex);
			//MPISynchronizeGridHaloJustRecv(roughxindex); //LAURA
			//MPISynchronizeGridHaloJustSend(roughxindex); //LAURA
			//MPISynchronizeGridHaloJustWait(roughxindex); //LAURA

			MPISynchronizeGridHalo(roughyindex);
			//MPISynchronizeGridHaloJustRecv(roughyindex); //LAURA
			//MPISynchronizeGridHaloJustSend(roughyindex); //LAURA
			//MPISynchronizeGridHaloJustWait(roughyindex); //LAURA
#endif
			do
			{
				iter++;
#ifdef _MPI
				MPISynchronizeGridHalo(headindex);
				//MPISynchronizeGridHaloJustRecv(headindex); //LAURA
				//MPISynchronizeGridHaloJustSend(headindex); //LAURA
				//MPISynchronizeGridHaloJustWait(headindex); //LAURA

				MPISynchronizeGridHalo(headoldindex);
				//MPISynchronizeGridHaloJustRecv(headoldindex); //LAURA
				//MPISynchronizeGridHaloJustSend(headoldindex); //LAURA
				//MPISynchronizeGridHaloJustWait(headoldindex); //LAURA

#endif

//roughness can also be estimated per iteration, but then the performance is not as good, although the results are a bit smoother, usually there is more retention effect because R tends to lower values
//				effective_roughness();
//#ifdef _MPI
//				MPISynchronizeGridHalo(roughxindex);
//				MPISynchronizeGridHalo(roughyindex);
//#endif

				maxerror = gauss_seidel_iteration(global_mintime,iter);

#ifdef _MPI
				// and also mpi-reduce maximum of maxerror
				MPI_Reduce(&maxerror, &global_maxerror, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
				MPI_Bcast(&global_maxerror, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
				global_maxerror = maxerror;
#endif
			} while (((iter<maxiter) && (global_maxerror>maxerrorf)) || (iter < 3));

		}

		// FLUX-stacks need to be synchronized for all MPI ranks
#ifdef _MPI
		MPISynchronizeGridHalo(headindex); // all instances need the heighboring informations about HEAD
		//MPISynchronizeGridHaloJustRecv(headindex); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex); //LAURA

		MPISynchronizeGridHalo(headoldindex); // all instances need the heighboring informations about HEAD
		//MPISynchronizeGridHaloJustRecv(headoldindex); //LAURA
		//MPISynchronizeGridHaloJustSend(headoldindex); //LAURA
		//MPISynchronizeGridHaloJustWait(headoldindex); //LAURA
		
		MPISynchronizeStackHalo(fluxindex, 10);
#endif
		if (method > 5)
		{
			GetFlowVelocities_AfterGaussSeidel();
		}
		// now we know the flow times for the input conditions (SSSU counts for slope calculation but not QDgrid (resp. only when SSSU > QDgrid)
		// so we have to get the minimum flow time -> must not be in a parallel section
		if (step == 0)
		{
			mintime = MinimumFlowTime(step);
#ifdef _MPI
			// reducing mintime and bcast it to every instance
			MPI_Reduce(&mintime, &global_mintime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			MPI_Bcast(&global_mintime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
			global_mintime=mintime;
#endif
			// round to the next smaller integer fraction of InV, so that InV is an integer multiple of mintime
			numsteps = ((int)((InV / substepcnt*60.0) / global_mintime) + 1);
			global_mintime = (InV / substepcnt*60.0) / numsteps;
			if (global_mintime<dMinSubTimeStep)
			{
				global_mintime = dMinSubTimeStep;
			}
			else if (global_mintime>dMaxSubTimeStep)
			{
				global_mintime = dMaxSubTimeStep;
			}
			numsteps = ((int)((InV / substepcnt*60.0) / global_mintime) + 1);
		}
		if (step == NextMajorStep)
		{
			NextMajorStep += numsteps / ((InV / substepcnt*60.0) / dMaxSubTimeStep);
		}

#ifdef _MPI
		//		MPI_Barrier(MPI_COMM_WORLD); // wait here for all the instances to resume the core loop 
#endif

//		char filename[1000];
//		sprintf(filename, "debug_surfacerouting_head_%.04d%.02d%.02d%.2d_loop_%.02d_step_%.03d.grd", jr, mt, tg, se, moduleSubstepCounter,step);
//		printgridfile(filename, HEAD, true, true);
	
		// now we can run the rest of the schema in a parallelized section again
#pragma omp parallel for schedule(runtime) shared(step)

		// again: do not use any common variable in the openmp-loop for writing! 

		// calculate the amount of water flowing into each cell (filling SSSU) from QDgrid and all inflows 
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					if (method < 6)
					{
						// here, we are not source but target oriented, meaning that not flows INTO other cells are 
						// calculated but flows coming FROM other cells into the actual cell. Thus, we can 
						// still stay in a parallel section without writing into other cells data (only reading the flow from these cells).
						CalculateWaterFlows(i, j, global_mintime);

						// do one really need to separate these two loops? there is no access to any other cell apart from the actual i,j indices, so it's probably ok to merge both loops
						//						}
						//					}
						//				}
						//		
						//		#pragma omp parallel for schedule(runtime) 
						//				// balancing these flows
						//				for (int i=1;i<=nrows;i++) /* Reihen des Grids */
						//				{
						//					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
						//					{
						//						if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
						//						{

						SRQIN[i][j] = 1000.0*SRQIN[i][j] / (cellsize*cellsize); // conversion into mm
						SRQOUT[i][j] = 1000.0*SRQOUT[i][j] / (cellsize*cellsize); // conversion into mm

						double TracerMassIn;
						double TracerMassOut;

						// balancing tracer concentrations according to the fluxes between the cells
						// this must be done as an integral part of the surfave routing for each sub interval
						for (int z = 1; z <= tracer_count; z++)
						{
							TracerMassOut = SRQOUT[i][j] * tracer_module[z].concstor[i][j]; // from Kg/l/m2 into Kg/m2
							TracerMassIn = tracer_module[z].concflux[i][j] / (cellsize*cellsize); // from Kg in Kg/m2
							tracer_module[z].concstor[i][j] = SSSU[i][j] * tracer_module[z].concstor[i][j] + (TracerMassIn - TracerMassOut); // in kg/m2
						}

						// balancing the water table in each cell from old storage, inflow and outflow
						SSSU[i][j] = SSSU[i][j] + SRQIN[i][j] - SRQOUT[i][j];
						if (SSSU[i][j] < 0)
						{
							SRQOUT[i][j] += SSSU[i][j];
							SSSU[i][j] = 0.0;
						}

						for (int z = 1; z <= tracer_count; z++)
						{
							if (SSSU[i][j] > 0)
							{
								tracer_module[z].concstor[i][j] /= (SSSU[i][j]); // now as concentration;
							}
							else
							{
								tracer_module[z].concstor[i][j] = 0;
							}
						}
						// we will now have to set SSSU to 0 in order to avoid 
						// the same quantity of water and tracers to be regarded as input into the silting up model in the next time step, again
						// SSSU is the storage of routed but not infiltrated water.
						if ((RIVDEP[i][j] > 0) && (RIVWIT[i][j] > 0))
						{
							QDgrid[i][j] += SSSU[i][j]; // in mm for statistics file
							for (int z = 1; z <= tracer_count; z++)
							{
								CONC_QD[z][i][j] += tracer_module[z].concstor[i][j] * SSSU[i][j];
								tracer_module[z].concstor[i][j] = 0;
							}
							SSSU[i][j] = 0;
						}
						else
						{
							QDgrid[i][j] = 0;
							for (int z = 1; z <= tracer_count; z++)
							{
								CONC_QD[z][i][j] = 0;
							}
						}
						// recalculation of total inflow for the next sup time step
						HEAD[i][j] = SSSU[i][j] / 1000.0; // in m for next sub time step
					}
					else
					{ // Gauss-Seidel method without tracers yet!
//						if ((RIVDEP[i][j] > 0) && (RIVWIT[i][j] > 0))
//						{
//							// FLUX[i][j][0] was already computed in GetFlowVelocities_AfterGaussSeidel() for each cell, so we only need to convert that value from m3/s into mm
//							double SumFromFlux = zoneQDFromFlux[zarray[(int)zn[i][j]]].sx;
//							double actFromFlux = (fabs(FLUX[i][j][0]) * 1000 * actualStepLength) / (CellSizeX[i][j] * CellSizeY[i][j]);
//							double scaling = (SumFromFlux>0)?QGScalings[zarray[(int)zn[i][j]]] / SumFromFlux:0;
//							QDgrid[i][j] += scaling*actFromFlux; // converting FLUX[i][j][0] from m3/s in mm for statistics of QD
//							//QDgrid[i][j] += QGScalings[zarray[(int)zn[i][j]]];
//						}
//						else
//						{
//							QDgrid[i][j] = 0;
//						}
						// HEAD[i][j] already contains the new values, so we only need to define the effective runoff (influx) to river cells
						if (HEAD[i][j] < DHM[i][j])
						{
							HEAD[i][j] = DHM[i][j];
						}
						SSSU[i][j] = (HEAD[i][j] - DHM[i][j])*1000.0;
						HEADold[i][j] = HEAD[i][j];
						if (SSSU[i][j] < 0)
						{
							SSSU[i][j] = 0.0;
						}
					}
				}
			}
		}

#ifdef _MPI
		MPISynchronizeGridHalo(headindex); // all instances need the heighboring informations about HEAD
		//MPISynchronizeGridHaloJustRecv(headindex); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex); //LAURA

		MPISynchronizeGridHalo(headoldindex); // all instances need the heighboring informations about HEAD
		//MPISynchronizeGridHaloJustRecv(headoldindex); //LAURA
		//MPISynchronizeGridHaloJustSend(headoldindex); //LAURA
		//MPISynchronizeGridHaloJustWait(headoldindex); //LAURA

		MPISynchronizeStackHalo(fluxindex, 10);
#endif

		step++;
		subtimestepsum += global_mintime;
		actualStepLength = global_mintime;
		
	}
	
//#pragma omp parallel for schedule(runtime) shared(step)
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
			{
				// for scaling of the QDgrid values with the change in overall HEAD (which should be identical, since suring routing there are no ther sources or sinks)
				zoneHEADEnd[zarray[(int)zn[i][j]]].n++;
				zoneHEADEnd[zarray[(int)zn[i][j]]].sx += (HEAD[i][j] - DHM[i][j])>0 ? (HEAD[i][j] - DHM[i][j]):0;
			}
		}
	}

	MPIGatherZone(zoneHEADStart, zonmax);
	MPIScatterZone(zoneHEADStart, zonmax);

	MPIGatherZone(zoneHEADEnd, zonmax);
	MPIScatterZone(zoneHEADEnd, zonmax);

	MPIGatherZone(zoneQDFromFlux, zonmax);
	MPIScatterZone(zoneQDFromFlux, zonmax);

	for (int i = 1; i <= zonmax; i++)
	{
		//double s = (zoneHEADStart[i].n > 0) ? zoneHEADStart[i].sx / zoneHEADStart[i].n : 0;
		//double e = (zoneHEADEnd[i].n > 0) ? zoneHEADEnd[i].sx / zoneHEADEnd[i].n : 0;
		double s = zoneHEADStart[i].sx;
		double e = zoneHEADEnd[i].sx;

		QGScalings[i] = 1000*(s - e); // positive formulation for comparison with QDgrid and conversion from m to mm
		if (QGScalings[i] < 0) QGScalings[i] = 0;

		numCellsQD[i] = 0.0;

	}
	
	if (method==6)
	{
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
				{
					if ((RIVDEP[i][j] > 0) && (RIVWIT[i][j] > 0))
					{
						numCellsQD[zarray[(int)zn[i][j]]]+=1.0;
					}
				}
			}
		}
	}

#ifdef _MPI
	double *numCellsQDglobal = new double[zonmax + 1];
	MPI_Allreduce(&(numCellsQD[1]), &(numCellsQDglobal[1]), zonmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (int i = 1; i <= zonmax; i++)
	{
		numCellsQD[i] = numCellsQDglobal[i];
	}
	delete[] numCellsQDglobal;
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	
	// statistics must be done outside the parallel section because common variables will be accessed by wertzufuegen
	// as well as by tracer_module[z].stat->zone[].sx and other operators
	// this is sample code, copied from unsatzonmodellklasse.
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) 
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				int TG = zarray[(int)zn[i][j]]; /* which subbasin */
				for (int z=1;z<=tracer_count;z++)
				{
					CONC_SSSU[z][i][j]=tracer_module[z].concstor[i][j];
					tracer_module[z].concstor[i][j]=0;
				}

				// normalizing concentrations (where added as multiple of conc times amount)
				if ((RIVDEP[i][j]>0)&&(RIVWIT[i][j]>0))
				{
					for (int z=1;z<=tracer_count;z++)
					{
						if (QDgrid[i][j]>0)
						{
							CONC_QD[z][i][j]/=QDgrid[i][j];
						}
						else
						{
							CONC_QD[z][i][j]=0.0;
						}
					}
				}

				/* updating the statistics */
				SFSTSRstat->wertzufuegen(i,j,TG,SSSU[i][j]);

				/* getting quick surface flow into the statistic register */
				int fzt_class;
				if (FZT)
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
					// only water reaching a river will be taken into account
					if ((RIVDEP[i][j]>0)&&(RIVWIT[i][j]>0)&& numCellsQD[zarray[(int)zn[i][j]]]!=0)
					{
						QDgrid[i][j] = QGScalings[zarray[(int)zn[i][j]]] / numCellsQD[zarray[(int)zn[i][j]]];
						zoneFZT_QD[0][TG][fzt_class].sx += RelCellArea[i][j]*QDgrid[i][j]*substepcnt;
					}
					else
					{
						QDgrid[i][j] = 0;
					}
					zoneFZT_QD[0][TG][fzt_class].n  += 1.0;
					zoneFZT_QD[0][TG][fzt_class].RelCellSum  += RelCellArea[i][j];

					zoneFZT_QI[0][TG][fzt_class].sx += RelCellArea[i][j]*QIgrid[i][j];
					zoneFZT_QI[0][TG][fzt_class].n  +=1.0;
					zoneFZT_QI[0][TG][fzt_class].RelCellSum  += RelCellArea[i][j];
				}
				else
				{
					// only water reaching a river will be taken into account
					if ((RIVDEP[i][j]>0)&&(RIVWIT[i][j]>0) && numCellsQD[zarray[(int)zn[i][j]]] != 0)
					{
						QDgrid[i][j] = QGScalings[zarray[(int)zn[i][j]]]/numCellsQD[zarray[(int)zn[i][j]]];
						zoneFZT_QD[0][TG][fzt_class].sx += QDgrid[i][j];
					}
					else
					{
						QDgrid[i][j] = 0;
					}
					zoneFZT_QD[0][TG][fzt_class].n  +=1.0;
					zoneFZT_QD[0][TG][fzt_class].RelCellSum  += 1.0;

					zoneFZT_QI[0][TG][fzt_class].sx += QIgrid[i][j];
					zoneFZT_QI[0][TG][fzt_class].n  +=1.0;
					zoneFZT_QI[0][TG][fzt_class].RelCellSum  += 1.0;
				}
				for (int z=1;z<=tracer_count;z++)
				{
					// for tracers the averaging has to be carried out weigted by absolute values 
					if (RelCellArea!=NULL)
					{
						// only water reaching a river will be taken into account
						if ((RIVDEP[i][j]>0)&&(RIVWIT[i][j]>0))
						{
							zoneFZT_QD[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
								RelCellArea[i][j]*QDgrid[i][j]*CONC_QD[z][i][j]*substepcnt : RelCellArea[i][j]*CONC_QD[z][i][j]*substepcnt; // if evaporating: handling masses, now konverted to concentrations 
						}
						zoneFZT_QD[z][TG][fzt_class].n  += QDgrid[i][j];
						zoneFZT_QD[z][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*QDgrid[i][j];
						
						zoneFZT_QI[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
							RelCellArea[i][j]*CONC_QI[z][i][j]*QIgrid[i][j] : RelCellArea[i][j]*CONC_QI[z][i][j];
						zoneFZT_QI[z][TG][fzt_class].n  += QIgrid[i][j];
						zoneFZT_QI[z][TG][fzt_class].RelCellSum  += RelCellArea[i][j]*QIgrid[i][j];
					}
					else
					{
						// only water reaching a river will be taken into account
						if ((RIVDEP[i][j]>0)&&(RIVWIT[i][j]>0))
						{
							zoneFZT_QD[z][TG][fzt_class].sx += tracer_info[z].evaporating ? 
								QDgrid[i][j]*CONC_QD[z][i][j]*substepcnt : CONC_QD[z][i][j]*substepcnt; // if evaporating: handling masses, now konverted to concentrations 
						}
						zoneFZT_QD[z][TG][fzt_class].n  += QDgrid[i][j];
						zoneFZT_QD[z][TG][fzt_class].RelCellSum  += QDgrid[i][j];

						zoneFZT_QI[z][TG][fzt_class].sx += tracer_info[z].evaporating ? CONC_QI[z][i][j]*QIgrid[i][j] : CONC_QI[z][i][j];
						zoneFZT_QI[z][TG][fzt_class].n  += QIgrid[i][j];
						zoneFZT_QI[z][TG][fzt_class].RelCellSum  += QIgrid[i][j];
					}
				}

				// tracer statistics for the output will contain the tracer concentrations within 
				// the surface storage after surface routing (which should be stored preveously in SSSU and concstor resp.)
				for (int z=1;z<=tracer_count;z++)
				{
					if (RelCellArea!=NULL)
					{
						tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concstor[i][j]*RelCellArea[i][j]*(SSSU[i][j]);
						tracer_module[z].stat->zone[TG].n+=(SSSU[i][j]);
						tracer_module[z].stat->zone[TG].RelCellSum+=(SSSU[i][j])*RelCellArea[i][j];
					}
					else
					{
						tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concstor[i][j]*(SSSU[i][j]);
						tracer_module[z].stat->zone[TG].n+=(SSSU[i][j]);
						tracer_module[z].stat->zone[TG].RelCellSum+=(SSSU[i][j]);
					}
				}
				//		QDgrid[i][j]=0;
			}
		}
	}
// statistics will be done in the post_run method in order to make the model ready for multiple loops within one regular time step
// for unsatzonmodel, surface routing and gwflow (and etr)
}

void surfaceroutingclass::post_run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
//	if (mpi_rank>0) return; // this module is currently not able to run on a cluster. Only the instance running mpi rank 0 will do surface routing in order to keep the model running at all
	// discharge components for subbasins --> copy from unsatzonmodell with a few changes (ELS uses ParA here as kelsqd, QI is already routed through the SLR-cascade QB is already known a.s.o)
	double t1,t2,t3,els_fact;

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
	MPIGatherZone(zoneQBSR,zonmax,true);
//	MPIGatherZone(zoneQI,zonmax,true); not neccessary, since zoneQD and zoneQI will be initialized hereafter from zoneFZT_QD and zoneFZT_QI, resp.
//	MPIGatherZone(zoneQD,zonmax,true);


	if (mpi_rank==0)
	{
		for (int i=1;i<=zonmax;i++)
		{
			QD[0][i]  = zoneFZT_QD[0][i][1].sx/zoneFZT_QD[0][i][1].RelCellSum; // FZT[i][1] is defined, always 
			for (int z=1;z<=tracer_count;z++)
			QD[z][i]  = zoneFZT_QD[z][i][1].RelCellSum>0 ? zoneFZT_QD[z][i][1].sx/zoneFZT_QD[z][i][1].RelCellSum : 0;

			//    QDSRstat->wertzufuegen(0,0,i,QD[0][i]);
			els_fact=exp(-InV/(kelsqd[i]));
			t3=QD[0][i]*(1-els_fact) + QDd[0][i]*els_fact;
			for (int z=1;z<=tracer_count;z++)
			{
				QDd[z][i] = QD[z][i]*QD[0][i]*(1-els_fact) + QDd[z][i]*QDd[0][i]*els_fact; // adding tracer masses for influx and storage outflow 
				QDd[z][i] = t3>0 ? QDd[z][i]/t3:0.0; // converting masses to concentrations 
			}
			QDd[0][i] = t3; // single linear storage for surface runoff
	//		QDSRstat->wertzufuegen(0,0,i,QDd[0][i]);
			zoneQDSR[i].sx			= QDd[0][i]*zoneQDSR[i].rel;
			zoneQDSR[i].n			= zoneQDSR[i].rel;
			zoneQDSR[i].RelCellSum	= zoneQDSR[i].rel;

			QI[0][i]  = zoneFZT_QI[0][i][1].sx/zoneFZT_QI[0][i][1].RelCellSum;    // average subbasin value of interflow is calculated to add it to the baseflow for getting later total flow
			for (int z=1;z<=tracer_count;z++)
			QI[z][i]  = zoneFZT_QI[z][i][1].RelCellSum>0 ? zoneFZT_QI[z][i][1].sx/zoneFZT_QI[z][i][1].RelCellSum : 0.0;
			els_fact=exp(-InV/(kelsqi[i]));

			t2 = QI[0][i]*(1-els_fact) + QIi[0][i]*els_fact; /* single linear storage effect for interflow */
			for (int z=1;z<=tracer_count;z++)
			{
				QIi[z][i] = QI[z][i]*QI[0][i]*(1-els_fact) + QIi[z][i]*QIi[0][i]*els_fact; /* adding tracer masses for influx and storage outflow */
				QIi[z][i] = t2>0 ? QIi[z][i]/t2:0.0; /* konverting masses to concentrations */
			}
			QIi[0][i] = t2; /* single linear storage effect for interflow */
	//		QISRstat->wertzufuegen(0,0,i,QIi[0][i]);
			zoneQISR[i].sx			= QIi[0][i]*zoneQISR[i].rel;
			zoneQISR[i].n			= zoneQISR[i].rel;
			zoneQISR[i].RelCellSum	= zoneQISR[i].rel;

			Qg[0][i] = QDd[0][i] + QIi[0][i] + QB[0][i];
			//    Qg[0][i] = QD[0][i] + QIi[0][i] + QB[0][i];
	//		QBSRstat->wertzufuegen(0,0,i,QB[0][i]);
			zoneQBSR[i].sx			= QB[0][i]*zoneQBSR[i].rel;
			zoneQBSR[i].n			= zoneQBSR[i].rel;
			zoneQBSR[i].RelCellSum	= zoneQBSR[i].rel;

	//		QGSRstat->wertzufuegen(0,0,i,Qg[0][i]);
			zoneQGSR[i].sx			= Qg[0][i]*zoneQGSR[i].rel;
			zoneQGSR[i].n			= zoneQGSR[i].rel;
			zoneQGSR[i].RelCellSum	= zoneQGSR[i].rel;
		}
	}

	QDSRstat->calc_and_done(jahr,monat,tag,stunde, minute);
	QISRstat->calc_and_done(jahr,monat,tag,stunde,minute);
	QGSRstat->calc_and_done(jahr,monat,tag,stunde,minute);
	QBSRstat->calc_and_done(jahr,monat,tag,stunde,minute);
	SFSTSRstat->calc_and_done(jahr,monat,tag,stunde, minute);

	// mean zonal concentrations now have to be copied into the Qg-Array -> taken from unsatzon with a few variations
	for (int z=1;z<=tracer_count;z++)
	{
		for (int i=1;i<=zonmax;i++)
		{
			/* superposition weighed by the mean fractions of the runoff components for subbasins */
			//		  if ((QB[0][i]+QIi[0][i]+QD[0][i])>0)
			if ((QB[0][i]+QIi[0][i]+QDd[0][i])>0)
			{
				Qg[z][i] = (QB[0][i]*QB[z][i] + QIi[0][i]*QIi[z][i] + QDd[0][i]*QDd[z][i])/(QB[0][i]+QIi[0][i]+QDd[0][i]);
				//			Qg[z][i] = (QB[0][i]*QB[z][i] + QIi[0][i]*QIi[z][i] + QD[0][i]*QD[z][i])/(QB[0][i]+QIi[0][i]+QD[0][i]);
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
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde, minute);
		//	  writegridtofile(tracer_module[z].fluxindex);
		writegridtofile(concsssuindex[z]);
	}

	writegridtofile(srfcvelindex);
	writestacktofile(srfcvelindex);

	writegridtofile(fluxindex);
	writestacktofile(fluxindex);
}

/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void surfaceroutingclass::ClearMemory(bool appendmode)
{
	if (do_it > 0)
	{
		if (initialized) WriteSlrContents(appendmode);
		substance_transport_delete();

		// free allocated memeory
		QDSRstat->ClearMemory();
		QISRstat->ClearMemory();
		QBSRstat->ClearMemory();
		QGSRstat->ClearMemory();
		SFSTSRstat->ClearMemory();
		
		if (QDSRstat) delete(QDSRstat);
		if (QISRstat) delete(QISRstat);
		if (QBSRstat) delete(QBSRstat);
		if (QGSRstat) delete(QGSRstat);
		if (SFSTSRstat) delete(SFSTSRstat);

		if (zoneQDSR)  free(zoneQDSR);  zoneQDSR=NULL;
		if (zoneQISR)  free(zoneQISR);  zoneQISR=NULL;
		if (zoneQBSR)  free(zoneQBSR);  zoneQBSR=NULL;
		if (zoneQGSR)  free(zoneQGSR);  zoneQGSR=NULL;
		if (zoneSFSTSR)  free(zoneSFSTSR);  zoneSFSTSR=NULL;

		if (zoneHEADStart) free(zoneHEADStart); zoneHEADStart = NULL;
		if (zoneHEADEnd) free(zoneHEADEnd); zoneHEADEnd = NULL;
		if (zoneQDFromFlux) free(zoneQDFromFlux); zoneQDFromFlux = NULL;
		delete[] QGScalings;

		if (QDSRstatfilename) free(QDSRstatfilename); QDSRstatfilename=NULL;
		if (QISRstatfilename) free(QISRstatfilename); QISRstatfilename=NULL;
		if (QBSRstatfilename) free(QBSRstatfilename); QBSRstatfilename=NULL;
		if (QGSRstatfilename) free(QGSRstatfilename); QGSRstatfilename=NULL;
		if (SFSTSRstatfilename) free(SFSTSRstatfilename); SFSTSRstatfilename=NULL;

		if (storagefilename) free(storagefilename);

		//    free_matrix(gridinfo[grid1index].pointer,1,nrows,1,ncols);
		//    delete_grid(grid1index);
		if (maxzeit) free(maxzeit);
		if ((zoneFZT_QD)&&(zoneFZT_QI))
		{
			for (int z=0;z<=tracer_count;z++)
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

		free_drei_d_matrix(gridinfo[srfcvelindex].stackpointer, gridinfo[srfcvelindex].nrl, gridinfo[srfcvelindex].nrh, 1, ncols, 10, nrows);
		free_matrix(gridinfo[srfcvelindex].pointer, gridinfo[srfcvelindex].nrl, gridinfo[srfcvelindex].nrh, 1, ncols, nrows);
		
		free_drei_d_matrix(gridinfo[fluxindex].stackpointer, gridinfo[fluxindex].nrl, gridinfo[fluxindex].nrh, 1, ncols, 10, nrows);
		free_matrix(gridinfo[fluxindex].pointer, gridinfo[fluxindex].nrl, gridinfo[fluxindex].nrh, 1, ncols, nrows);
		
		free_drei_d_matrix(gridinfo[weightindex].stackpointer, gridinfo[weightindex].nrl, gridinfo[weightindex].nrh, 1, ncols, 10, nrows);
		free_matrix(gridinfo[weightindex].pointer, gridinfo[weightindex].nrl, gridinfo[weightindex].nrh, 1, ncols, nrows);

		free_matrix(gridinfo[headindex].pointer, gridinfo[headindex].nrl, gridinfo[headindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[headoldindex].pointer, gridinfo[headoldindex].nrl, gridinfo[headoldindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[roughxindex].pointer, gridinfo[roughxindex].nrl, gridinfo[roughxindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[roughyindex].pointer, gridinfo[roughyindex].nrl, gridinfo[roughyindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[manningindex].pointer, gridinfo[manningindex].nrl, gridinfo[manningindex].nrh, 1, ncols, nrows);

		if (plancurvindex>0) free_matrix(gridinfo[plancurvindex].pointer, gridinfo[plancurvindex].nrl, gridinfo[plancurvindex].nrh, 1, ncols, nrows);
		if (srqinindex>0)    free_matrix(gridinfo[srqinindex].pointer, gridinfo[srqinindex].nrl, gridinfo[srqinindex].nrh, 1, ncols, nrows);
		if (srqoutindex>0)   free_matrix(gridinfo[srqoutindex].pointer, gridinfo[srqoutindex].nrl, gridinfo[srqoutindex].nrh, 1, ncols, nrows);
		if (taucritindex>0)   free_matrix(gridinfo[taucritindex].pointer, gridinfo[taucritindex].nrl, gridinfo[taucritindex].nrh, 1, ncols, nrows);

		delete_grid(srfcvelindex);
		delete_grid(fluxindex);
		delete_grid(weightindex);
		delete_grid(headindex);
		delete_grid(headoldindex);
		delete_grid(roughxindex);
		delete_grid(roughyindex);
		delete_grid(manningindex);
		if (plancurvindex>0)  delete_grid(plancurvindex);
		if (srqinindex>0)     delete_grid(srqinindex);
		if (srqoutindex>0)    delete_grid(srqoutindex);
		if (taucritindex>0)   delete_grid(taucritindex);
		
		if (LUSE) delete[] LUSE; LUSE=NULL;
		if (VCF) delete[] VCF; VCF=NULL;
	}
	modelobject::ClearMemory();
}

void surfaceroutingclass::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	int z;
	char *tempstr1,*tempstr2,*tempstr3,*tempstr4;
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (retvalue!=0) throw retvalue;

		fzt_grid_auswerten(); // identical to unsatzonmodellklasse

		if (tracer_count>0)
		{

			for (z=1;z<=tracer_count;z++)
			{
				// instead of using the cinput-grid, we use two grids for each tracer:
				// CONC_QD and CONC_QI. They both will be initialized later in this method
				tracer_module[z].cinput=NULL;
			}

			CONC_QD       = (double***)calloc(tracer_count+1,sizeof(double**));
			CONC_QI       = (double***)calloc(tracer_count+1,sizeof(double**));
			// CONC_SSSU will store the tracer concentrations for silting up module, if SSSU will contain a value > 0
			CONC_SSSU     = (double***)calloc(tracer_count+1,sizeof(double**));

			// looking for the first free entry in the gridinfo array
			int sti=1;
			while (gridinfo[sti].identification!=NULL)
			{
				sti ++;
				if (sti==gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no free grid anymore, allowed number of grids is "<<gridmaxanzahl<<"!"<<"[end-entry]";
					throw(12);
				}
			}

			int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
			int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;

			char lpszInternalName[200];
			SiltingUp=true;
			for (int z=1;z<=tracer_count;z++)
			{
				sprintf(lpszInternalName,"_CONC_QD_%i_",z);
				CONC_QD[z]=gridinfo[getindex(lpszInternalName)].pointer;
				sprintf(lpszInternalName,"_CONC_QI_%i_",z);
				CONC_QI[z]=gridinfo[getindex(lpszInternalName)].pointer;
				sprintf(lpszInternalName,"concstor_tracer_%d_%s",z,"SiltingUpModel");
				if (ist_vorhanden(lpszInternalName))
				{
					CONC_SSSU[z]=gridinfo[getindex(lpszInternalName)].pointer;
					concsssuindex[z]=-1;
				}
				else
				{
					SiltingUp=false;
					// create grids and stacks
					CONC_SSSU[z] = matrix(sr,er,1,ncols,nrows);
					initialize_grid(sti, tracer_module[z].gridname, CONC_SSSU[z], NULL, 0, lpszInternalName, tracer_info[z].writegrid, sr,er);
					for (int k=mpistartrow;k<=mpiendrow;k++) 
					{
						for (int j=1;j<=ncols;j++)
						{
							if ((zn[k][j] != nodata) && (MPIREG[k][j]=mpi_rank))
							{
								CONC_SSSU[z][k][j] = 0;
							} 
							else	
							{
								CONC_SSSU[z][k][j] = nodata;
							}
						}
					}
					concsssuindex[z] = sti;
					MPISynchronizeGridHalo(concsssuindex[z], CONC_SSSU[z]);
					//MPISynchronizeGridHaloJustRecv(concsssuindex[z], CONC_SSSU[z]); //LAURA
					//MPISynchronizeGridHaloJustSend(concsssuindex[z], CONC_SSSU[z]); //LAURA
					//MPISynchronizeGridHaloJustWait(concsssuindex[z], CONC_SSSU[z]); //LAURA
					sti++;
				}
			}
		}
		
		// do this is any case, wether tracers are used or not.
		// we initialize the QB, QI, QD and Qg_arrays here and not in the constructor, because the number of tracers was not known there.

		// allocating memeroy for zonal statistics as they are given to routing module (Qg) or used to calculate the zonal totals
		// index [0] in the 2-d-array means the water fluxes, indices 1..z mean tracers, the second index is the subbasin
		Qg=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z=0;z<=tracer_count;z++) Qg[z]=(double*)calloc(zonmax+1,sizeof(double));

		QD=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z=0;z<=tracer_count;z++) QD[z]=(double*)calloc(zonmax+1,sizeof(double));
		// surface runoff with ELS-effect applied
		QDd=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z=0;z<=tracer_count;z++) QDd[z]=(double*)calloc(zonmax+1,sizeof(double));
		
		QI=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z=0;z<=tracer_count;z++) QI[z]=(double*)calloc(zonmax+1,sizeof(double));
		// interflow with ELS-effect applied
		QIi=(double**)calloc(tracer_count+1,sizeof(double*));
		for (int z=0;z<=tracer_count;z++) QIi[z]=(double*)calloc(zonmax+1,sizeof(double));
		
		/* reading storage content file in order to initialize the flow travel time zones for interflow and surface flow */
		// is done in modapuzr now
//		if (gridslesen == 1)
//		{
//			ReadSlrContents();
//		}

		return;
	}
}

void surfaceroutingclass::substance_transport_delete()
{
	if (do_it==0) return;
	else
	{
		if (CONC_QD) free(CONC_QD); CONC_QD=NULL;
		if (CONC_QI) free(CONC_QI); CONC_QI=NULL;
		if (CONC_SSSU) free(CONC_SSSU); CONC_SSSU=NULL;

		if (Qg)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				if (Qg[z]) free(Qg[z]);
				Qg[z]=NULL;
			}
		}
		if (Qg) free(Qg);Qg = NULL;

		if (QD)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				if (QD[z]) free(QD[z]);
				QD[z]=NULL;
			}
		}
		if (QD) free(QD);QD = NULL;

		if (QDd)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				if (QDd[z]) free(QDd[z]);
				QDd[z]=NULL;
			}
		}
		if (QDd) free(QDd);QDd = NULL;

		if (QI)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				if (QI[z]) free(QI[z]);
				QI[z]=NULL;
			}
		}
		if (QI) free(QI);QI = NULL;

		if (QIi)
		{
			for (int z=0;z<=tracer_count;z++) 
			{
				if (QIi[z]) free(QIi[z]);
				QIi[z]=NULL;
			}
		}
		if (QIi) free(QIi);QIi = NULL;
		
		teilmodellklasse::substance_transport_delete();
	}
}

void surfaceroutingclass::StorageFileError(const char *message, double c1, double c2)
{
	(*WSstream) << "FATAL: error in storage file content reading " << message << ":\n"
		"value in actual control is:   " << c1 << "\n"
		"but value in storage file is: " << c2 << "\n"
		"WaSiM throws..." << "[end-entry]";
	throw(71);
}

void surfaceroutingclass::ReadFloatFromStoragefile(double *FValue, const char *ErrorMessage)
{
	double st_float;
	fscanf(storagefile, "%lf", &st_float);
	if (st_float != 0)
	{
		if (fabs((st_float - *FValue) / st_float)>1e-5) StorageFileError(ErrorMessage, (double)(*FValue), (double)st_float);
	}
	else
	{
		if (*FValue != 0) StorageFileError(ErrorMessage, (double)(*FValue), (double)st_float);
	}
}

void surfaceroutingclass::ReadSlrContents(char *_storagefilename, int flag)
{
	int i, j, z, k;
	int st_int, st_int2;
	char *satz;
	if (_storagefilename != NULL)
	{
		sprintf(this->storagefilename, "%s", _storagefilename);
	}
	// if readgrid in Unsatzonmodel is 2 or 4, then only the groundwater table is read in, but no storage file.
	if ((_storagefilename != NULL) && (flag != 2) && (flag != 4))
	{
		if (flag != 0) // read the storage file for unsatzonmodel.readgrid==1 or 3
		{

			satz = (char*)calloc(_stringsize_, sizeof(char));

			char *completedname = new char[_stringsize_ + 1];
			sprintf(completedname, "%s", storagefilename);
			ExtendFileName(DefaultInputDirectory, completedname);
			if ((storagefile = fopen(completedname, "rt")) == NULL)
			{
				(*WSstream) << "WARNING: error opening storage content file " << completedname <<
					". This means that the actual model run will not initialize the flow travel time zones for direct runoff and interflow, as it was before version 9.08.06." << "[end-entry]";
				//		throw(64);
				return;
			}

			do
			{
				satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (feof(storagefile))
				{
					(*WSstream) << "WARNING: error locking for the section [" << _SurfaceRoutingModell_ << "] in " << completedname <<
						". This means that the actual model run will not initialize the flow travel time zones for direct runoff and interflow, as it was before version 9.08.06." << "[end-entry]";
					// throw(71);
					return;
				}
			} while (!strstr(satz, _SurfaceRoutingModell_));
			delete[] completedname;

			/* TG_number */
			fscanf(storagefile, "%d", &st_int); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
			if (st_int != zonmax) StorageFileError("number of subbasins", (double)zonmax, (double)st_int);

			/* number of tracers */
			fscanf(storagefile, "%d", &st_int); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
			if (st_int != tracer_count) StorageFileError("number of tracers", (double)tracer_count, (double)st_int);

			for (z = 0; z <= tracer_count; z++)
			{
				/* actual tracer (0=runoff) */
				fscanf(storagefile, "%d", &st_int); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
				if (st_int != z) StorageFileError("actual tracers(0=runoff)", (double)z, (double)st_int);
				for (i = 1; i <= zonmax; i++)
				{
					/* actual tracer (0=runoff) */
					fscanf(storagefile, "%d%d", &st_int, &st_int2); satz = fgets(satz, _stringsize_, storagefile); // read the comment into a dummy string
					if (st_int != i) StorageFileError("actual zone", (double)i, (double)st_int);
					if (st_int2 != (int)maxzeit[i]) StorageFileError("maxzeit of actual subbasin", (double)((int)maxzeit[i]), (double)st_int);

					for (j = 1; j <= (int)maxzeit[i]; j++)
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
			fclose(storagefile);
			if (satz) free(satz);
		}
	}
	/* -----------ENDE STORAGEFILE -------------- */
}

void surfaceroutingclass::WriteSlrContents(bool appendmode)
{  /* writing contents of overland flow and interflow single linear reservoirs to a file named storage.ftz */
	int i, j, z;
//	if (mpi_rank==0)
//	{

		char *completedname = new char[_stringsize_ + 1];
		sprintf(completedname, "%s", storagefilename);
		ExtendFileName(DefaultOutputDirectory, completedname);
		if (appendmode == true) // if routing model was active
			storagefile = fopen(completedname, "a");
		else
			storagefile = fopen(completedname, "w");

		delete[] completedname;

		fprintf(storagefile, "[" _SurfaceRoutingModell_ "]\n");
		fprintf(storagefile, "%d (number of subasins)\n", zonmax);
		fprintf(storagefile, "%d (number of tracers)\n", tracer_count);
		for (z = 0; z <= tracer_count; z++)
		{
			fprintf(storagefile, "%d (0=runoff, 1...9 = tracers)\n", z);
			for (i = 1; i <= zonmax; i++)
			{
				fprintf(storagefile, "%d %d (zone %d and maxzeit)\n", i, (int)maxzeit[i], i);
				for (j = 1; j <= (int)maxzeit[i]; j++)
				{
					fprintf(storagefile, "%lg %lg %lg %lg %lg (zoneFZT_QD.n/.sx/.zonenwert/.rel/.RelCellSum)\n",
						zoneFZT_QD[z][i][j].n,
						zoneFZT_QD[z][i][j].sx,
						zoneFZT_QD[z][i][j].zonenwert,
						zoneFZT_QD[z][i][j].rel,
						zoneFZT_QD[z][i][j].RelCellSum);
					fprintf(storagefile, "%lg %lg %lg %lg %lg (zoneFZT_QI.n/.sx/.zonenwert/.rel/.RelCellSum)\n",
						zoneFZT_QI[z][i][j].n,
						zoneFZT_QI[z][i][j].sx,
						zoneFZT_QI[z][i][j].zonenwert,
						zoneFZT_QI[z][i][j].rel,
						zoneFZT_QI[z][i][j].RelCellSum);
				}
				fprintf(storagefile, "%g (QDd[%d][%d])\n", QDd[z][i], z, i);
				fprintf(storagefile, "%g (QIi[%d][%d])\n", QIi[z][i], z, i);
			}
		}
		fprintf(storagefile, "\n");
		fclose(storagefile);
//	}
	/* -----------ENDE STORAGEFILE -------------- */
}

void surfaceroutingclass::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
		if (QDSRstat) QDSRstat->flushbuffer();
		if (QISRstat) QISRstat->flushbuffer();
		if (QBSRstat) QBSRstat->flushbuffer();
		if (QGSRstat) QGSRstat->flushbuffer();
		if (SFSTSRstat) SFSTSRstat->flushbuffer();

		for (z=1;z<=tracer_count;z++)
		tracer_module[z].stat->flushbuffer();
	}
}

// this method gets the address of the 2D-array for all total flows incl. tracers 
// (first dimension are water fluxes and tracers, the 2nd dimension are the subbasins)
// the address is given to routing module as input.
double  **surfaceroutingclass::GetQgesAdresse()
{
	if (do_it!=0)
	{
		if (Qg!=NULL)
		{
			if (Qg[0]!=NULL)
			return (double**)Qg;
			else
			{
				(*WSstream) << "FATAL: Qg[0] (surface routing) not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
				throw(57);
				return NULL;
			}
		}
		else
		{
			(*WSstream) << "FATAL: Qg (surface routing) not initialized (may be, tracer module not linked? -> programmers error!)"<<"[end-entry]";
			throw(57);
			return NULL;
		}
	}
	else
	return NULL;
}

// this is used in the sample code in order to simulate unsatzonmodellklasse
void	 surfaceroutingclass::SetQBAddress(double **address)
{
	if (do_it==0) return;
	QB = address;
}

// if QI would be only a reference to QI in an other modul (e.g. unsatzonmodel)
// than use this function to set QI and do NOT allocate and deallocate memory for QI in this module
void	 surfaceroutingclass::SetQIAddress(double **address)
{
	if (do_it==0) return;
	QI = address;
}

// if QD would be only a reference to QD in an other modul, 
// than use this function to set QD and do NOT allocate and deallocate memory for QD in this module
void surfaceroutingclass::SetQDAddress(double **address)
{
	if (do_it==0) return;
	QD = address;
}

void surfaceroutingclass::SetKelsQIAddress(double *address)
{
	if (do_it!=0)
	{
		kelsqi=address;
	}
}

void surfaceroutingclass::SetKelsQDAddress(double *address)
{
	if (do_it!=0)
	{
		kelsqd=address;
	}
}



/* ++++++++++++++++++++ end of surface routing routines ++++++++++++++++++++++++++++++ */
