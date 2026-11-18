/*   $Revision: 48 $ 
$Date: 15.09.10 1:07 $  */
/*
(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

Utilities fuer WaSiM-ETH 4.1 und WaSiM-

wie Speicherverwaltung, Fehlermeldungen,
Grids Lesen, Initialisieren, schreiben, aufsummieren, mitteln usw.
Statistikdateien verwalten, schreiben,lesen
Boden- und Pflanzeneigenschaften aus Tabellen lesen
sowie im Speicher diverse Grids erzeugen (Aspect, Slope, k-Wert, nFK usw.)

Stand : Juli 97
Stand : September 97 Einfuegen einiger neuer functions zur Verarbeitung von Gridstapeln (matrix_drei_d)
						readstackfromfile (ohne Pruefung auf vollstaendige Werte bzw. Fehlwerte)
						printstackfile (ohne Auffuellen der Fehlwerte (im Zonengrid) mit nodata)
						writestacktofile (ohne Moeglichkeit der Bildung von Summen- oder Mittelwertstacks)

Stand : 12.12.1997
Version : 4          Pufferung von Dateiausgaben (in zonalstatistik und routstatistik) zur
						Beschleunigung der Modellierung durch effektiveren Plattenzugriff

Stand : 25.1.1998
Version : 4.5        Ausgabe von Minima/Maxima/Summe/Quadratsumme von Grids -> speicherintensiv!

Stand : 12.2.1998
Version : 4.6        Implmentierung Stofftransport (Einlesen derTracer_info tabelle)

Stand : 14.2.1998    Aenderung der Speicherverwaltung in matrix und d_matrix, damit es bei optimierenden
						Compilern keine Adressbereichsueberschreitungen mehr geben kann

date  : 20.2.1998     implementation of substance-transport -> reading the input table for tracers

date  : 16.6.1998     surface resitances (rsc) for Landuse-table is now interpolated dayly, so the array for the 12 monthly values 
						is extended by a 0-th and a 13-th value for interpolation for the first 15 and the last 15 days of the year

date  : 8.Juli 1998   reading a list of subbasins, for which outputs should be written (this is an extension to the concept of 
						output-codes for each statistic file, because if there a too much subbasins, output for all subbasins produces
						too much output data. Now: if code 5 is given, output is written only for subbasins indicated in the output list.
						function readsteuerdatei()

date    : march, 9 1999
version : 4.6a hour-entry can be given for start or for end of time step via min_to_int_end

date    : march, 19-23 1999
version : anazlyzing a region grid which is used for multiple regression regions


date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.
*/


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <signal.h>
#include <ctime>
#include <sstream>
#include <stddef.h>
#include "stringp.h"
#include "modconst.h"
#include "modef.h"
#include "wasim.h"
#include "gifreader.h"
#ifdef _NETCDF
#include "netcdfhlp.h"
#endif
#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;


/* ----------------------- Zonale Statistik ------------------------------------- */

void zonalstatistic::fileopen(const char *mode)
{
	if ((zonestat=fopen(dateiname,mode))==NULL)
	{
		(*WSstream) << "FATAL: error opening \n" << dateiname << "[end-entry]";
		throw 15;
	}
}

void zonalstatistic::fileclose()
{
	if ((fclose(zonestat))==EOF)
	{
		(*WSstream) << "FATAL: cannot close " << dateiname << "[end-entry]";
		throw 16;
	}
}

void zonalstatistic::flushbuffer()
/* for faster disk access the output is buffered in a way, that the output files are
opened only, when the buffer is flushed. With "setbuf" or "setvbuf" the buffering is
working for opened files only, and on some systems like Windows 95 it is not possible
to have more than 15 files open at the same time. So the files in WaSiM are opened
each time immediatly before a write access and are closed after. To reduce the number
of this actions, a internal buffer is used, using the binary format for the buffered data*/
{
	unsigned long i,j;
	if ((code>0)&&(bufcnt>0))
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		(*WSstream) << "write " << dateiname;
		fileopen("a+");
		for (i=1;i<=bufcnt;i++)
		{
			fprintf(zonestat,"%2ld\t%2ld\t%2ld\t%2ld",buffer[i].j,buffer[i].m,buffer[i].t,buffer[i].s);
			if (useMinutes) fprintf(zonestat, "\t%2ld", buffer[i].min);
			for (j=1;j<=wertcnt;j++)
			{
				if (prec<99)
				fprintf(zonestat,"\t%.*f",prec,buffer[i].wert[j]);
				else
				fprintf(zonestat,"\t%.6g",buffer[i].wert[j]);
			}
			fprintf(zonestat,"\n");
		}
		fileclose();
		(*WSstream) << "("<< (long)(bufcnt*((long)4*sizeof(long)+wertcnt*sizeof(double))) <<" of "<< buffersize << " buffer-bytes)"<<"[end-entry]";
		SetSignalAbfangfunktion();
		bufcnt=0;
		wertcnt=0;
	}
}

void zonalstatistic::wertzufuegen(int i,int j,int index,double x,double weight)
{
	// weight is a new parameter allowing also fractions of a cell to be avaluated. Value and area must be weighed
	if (code>0)
	{
		zone[index].n += weight;
		if ((i>0)&&(j>0)&&(i<=nrows)&&(j<=ncols)&&RelCellArea!=NULL)
		{
			if (RelCellArea[i][j]>nodata)
			{
				zone[index].RelCellSum+=RelCellArea[i][j]*weight;
				zone[index].sx += RelCellArea[i][j]*x*weight;
			}
			else
			{
				zone[index].RelCellSum+=weight;
				zone[index].sx += x*weight;
			}
		}
		else
		{
			zone[index].RelCellSum+=weight;
			zone[index].sx += x*weight;
		}
	}
}

void zonalstatistic::init()
{
	int i;
	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	if (code>0)
	{
		cnt++; /* Statistik fuer cnt intervalle addiert, wenn 2. bis 4 Stelle von code > 001 */
		for (i=0;i<=zonmax;i++)
		{
			zone[i].n  = 0.0;
			zone[i].sx = 0.0;
			zone[i].RelCellSum=0.0;
		}
	}
}

#ifdef _MPI
void zonesum(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype)
{
	for (int i=0;i<*len;i++)
	{
		inout[i].n+=in[i].n;
		inout[i].sx+=in[i].sx;
		//		inout[i].rel+=in[i].rel; --> do not add the already calculated relative values which must remain constant
		//		inout[i].zonenwert+=in[i].zonenwert; --> also the zonenwert is a constant
		inout[i].RelCellSum+=in[i].RelCellSum;
	}
}
#endif

void zonalstatistic::calc_and_done(long jahr,long monat,long tag,long stunde,long minute,int gather)
/* Schreibt die zonalen Ergebnisse in die Statistikdateien, Dezimalstellen p */
{
	int    i;
	double gesamtmittel;
	double gesamtsumme;
	double gesamtRelCelSum;
	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	if (code>0)
	{
		int mpi_rank=0;
		int mpi_size=1;
#ifdef _MPI

		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

		if (gather>0)
		{
			// method 1 is the recommended method, 2 to 4 are own implementations minimizing network traffic by using as few MPI calls as possible -> with not much effect as experiments show.
			int gathermethod=2;


			if (gathermethod==1)
			{
				// Version mit eigener Reduce Funktion
				MPI_Datatype zonevectortype;
				MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
				MPI_Type_commit(&zonevectortype);

				void(*zonesumfunction)(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype) = &zonesum;
				MPI_Op mySumOp;
				MPI_Op_create((MPI_User_function*)zonesumfunction, true, &mySumOp);

				statrecord *tempzone;
				tempzone=new statrecord[zonmax+1];

				for (i=0;i<=zonmax;i++)
				{
					tempzone[i].n=tempzone[i].sx=tempzone[i].RelCellSum=0;
				}

				MPI_Reduce(zone, tempzone, zonmax + 1, zonevectortype, mySumOp, 0,MPI_COMM_WORLD);

				for (i=0;i<=zonmax;i++)
				{
					zone[i].n=tempzone[i].n;
					zone[i].sx=tempzone[i].sx;
					zone[i].RelCellSum=tempzone[i].RelCellSum;
				}

				MPI_Op_free(&mySumOp);
				MPI_Type_free(&zonevectortype);
				delete[] tempzone;
			} else if (gathermethod == 2) {
				// JOERG this seems not to reduce all members of the struct or all vector elements correctly
				MPI_Datatype zonevectortype;
				MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
				MPI_Type_commit(&zonevectortype);
	
				void(*zonesumfunction)(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype) = &zonesum;
				MPI_Op mySumOp;
				MPI_Op_create((MPI_User_function*)zonesumfunction, true, &mySumOp);

				if (mpi_rank == 0) //LAURA
				{ //LAURA
					MPI_Reduce(MPI_IN_PLACE, zone, zonmax + 1, zonevectortype, mySumOp, 0, MPI_COMM_WORLD); //LAURA 
				}
				else{ //LAURA
					MPI_Reduce(zone, NULL, zonmax + 1, zonevectortype, mySumOp, 0, MPI_COMM_WORLD); //LAURA

				} //LAURA

				MPI_Op_free(&mySumOp);
				MPI_Type_free(&zonevectortype);

			} else if (gathermethod==3) {
				// Version ohne Reduce, aber mit send/Recv Funktion für Summenbildung, alle Zonen zusammen in langen Speicherblöcken
				double *n1,*sx1,*relcellsum1;
				double *ntot,*sxtot,*relcellsumtot;

				n1 = new double[zonmax+1];
				sx1 = new double[zonmax+1];
				relcellsum1 = new double[zonmax+1];
				ntot = new double[zonmax+1];
				sxtot = new double[zonmax+1];
				relcellsumtot = new double[zonmax+1];

				for (i=0;i<=zonmax;i++)
				{
					n1[i]=zone[i].n;
					sx1[i]=zone[i].sx;
					relcellsum1[i]=zone[i].RelCellSum;
				}

				if (mpi_rank>0)
				{
//				(*WSstream) << "sending statistic zones from rank "<<mpi_rank<<" for "<<this->dateiname<<"[end-entry]";
					MPI_Send(&(n1[1]), zonmax, MPI_DOUBLE, 0, 10000 + mpi_rank,MPI_COMM_WORLD);
					MPI_Send(&(sx1[1]), zonmax, MPI_DOUBLE, 0, 20000 + mpi_rank, MPI_COMM_WORLD);
					MPI_Send(&(relcellsum1[1]), zonmax, MPI_DOUBLE, 0, 30000 + mpi_rank, MPI_COMM_WORLD);
				}
				else
				{
					for (int j=1;j<=zonmax;j++)
					{
						ntot[j]  = n1[j];
						sxtot[j] = sx1[j];
						relcellsumtot[j] = relcellsum1[j];
					}
					for (int k=1;k<mpi_size;k++)
					{
//					(*WSstream) << "receiving statistic zones from rank "<<k<<" for "<<this->dateiname<<"[end-entry]";
						MPI_Status s;
						MPI_Recv(&(n1[1]), zonmax, MPI_DOUBLE, k, 10000 + k,MPI_COMM_WORLD,&s);
						MPI_Recv(&(sx1[1]), zonmax, MPI_DOUBLE, k, 20000 + k, MPI_COMM_WORLD, &s);
						MPI_Recv(&(relcellsum1[1]), zonmax, MPI_DOUBLE, k, 30000 + k, MPI_COMM_WORLD, &s);
						for (int j=1;j<=zonmax;j++)
						{
							ntot[j]  += n1[j];
							sxtot[j] += sx1[j];
							relcellsumtot[j] += relcellsum1[j];
						}
					}
				}

				if (mpi_rank==0)
				{
					for (i=0;i<=zonmax;i++)
					{
						zone[i].n=ntot[i];
						zone[i].sx=sxtot[i];
						zone[i].RelCellSum=relcellsumtot[i];
					}
				}

				delete[] n1;
				delete[] sx1;
				delete[] relcellsum1;
				delete[] ntot;
				delete[] sxtot;
				delete[] relcellsumtot;
			} else if (gathermethod==4) {
				// Version ohne Reduce, aber mit send/Recv Funktion für Summenbildung, alle Zonen zusammen in langen Speicherblöcken
				double *n1,*sx1,*relcellsum1;
				double *ntot,*sxtot,*relcellsumtot;

				n1 = new double[zonmax+1];
				sx1 = new double[zonmax+1];
				relcellsum1 = new double[zonmax+1];
				ntot = new double[zonmax+1];
				sxtot = new double[zonmax+1];
				relcellsumtot = new double[zonmax+1];

				double *vtemp = new double[3*(zonmax+1)];
				for (i=0;i<=zonmax;i++)
				{
					vtemp[i]=zone[i].n;
					vtemp[i+zonmax]=zone[i].sx;
					vtemp[i+2*zonmax]=zone[i].RelCellSum;
				}

				if (mpi_rank>0)
				{
	//				(*WSstream) << "sending statistic zones from rank "<<mpi_rank<<" for "<<this->dateiname<<"[end-entry]";
					MPI_Send(&(vtemp[1]), 3 * zonmax, MPI_DOUBLE, 0, 3000 + mpi_rank,MPI_COMM_WORLD);
				}
				else
				{
					for (int j=1;j<=zonmax;j++)
					{
						ntot[j]  = n1[j];
						sxtot[j] = sx1[j];
						relcellsumtot[j] = relcellsum1[j];
					}
					for (int k=1;k<mpi_size;k++)
					{
	//					(*WSstream) << "receiving statistic zones from rank "<<k<<" for "<<this->dateiname<<"[end-entry]";
						MPI_Status s;
						MPI_Recv(&(vtemp[1]), 3 * zonmax, MPI_DOUBLE, k, 3000 + k, MPI_COMM_WORLD, &s);
						for (int j=1;j<=zonmax;j++)
						{
							ntot[j]  += vtemp[j];
							sxtot[j] += vtemp[j+zonmax];
							relcellsumtot[j] += vtemp[j+2*zonmax];
						}
					}
				}

			
				if (mpi_rank==0)
				{
					for (i=0;i<=zonmax;i++)
					{
						zone[i].n=ntot[i];
						zone[i].sx=sxtot[i];
						zone[i].RelCellSum=relcellsumtot[i];
					}
				}

				delete[] n1;
				delete[] sx1;
				delete[] relcellsum1;
				delete[] ntot;
				delete[] sxtot;
				delete[] relcellsumtot;
			}

		}
#endif

		if (mpi_rank==0)
		{
			if (intervalle>1)
			for (i=0;i<=zonmax;i++)
			{
				zone[i].sx+=sumsum[i]; /* updateing fuer laufende Summen ueber mehrere Intervalle */
				if ((code==1)||(code==2)||(code==5))
				{
					zone[i].n+=sumn[i]; /* wenn Mittelwerte gebildet werden und nicht Summen ueber mehrere Intervalle, muss n von 1 an gezaehlt werden */
					zone[i].RelCellSum+=sumRelCellSum[i];
				}
			}
			gesamtmittel=0; /* Gebietsmittelwert fuer ein einzelnes Intervall */
			if (cnt==intervalle)
			{
				wertcnt=0;
				bufcnt++;
				buffer[bufcnt].j=jahr;
				buffer[bufcnt].m=monat;
				buffer[bufcnt].t=tag;
				buffer[bufcnt].s=stunde;
				buffer[bufcnt].min = minute;
				lastresults.j=jahr;
				lastresults.m=monat;
				lastresults.t=tag;
				lastresults.s=stunde;
				lastresults.min = minute;
			}
			gesamtsumme = 0;
			gesamtRelCelSum = 0;
			for (i=1;i<=zonmax;i++)
			{
				/* sollen die Zonenmittelwerte jeweils ueber mehrere Intervalle adiert werden,
			so wird das durch die jeweils am Beginn eines Intervalles erfolgende Nullsetzen
			des n-Zaehlers im Zonenrecord beruecksichtigt. Bei Mittelwertbildung ueber
			mehrere Intervalle wird zone[i].n nicht auf 0 gesetzt, das erfolgt erst nach
			dem Schreiben nach der entsprechenden Anzahl Intervallen */

				//      if (zone[i].n>0) mittelwert=zone[i].sx/zone[i].n; else mittelwert = 0;
				if (zone[i].n>0) mittelwert = zone[i].sx/zone[i].RelCellSum; else mittelwert = 0;
//				gesamtmittel += mittelwert * zone[i].rel;
				// neu 16.01.2016
				gesamtsumme += zone[i].sx;
				gesamtRelCelSum += zone[i].RelCellSum;

				if (cnt==intervalle)
				{
					if ((code==2)||(code==4))
					{
						buffer[bufcnt].wert[i]=mittelwert;
						lastresults.wert[i]=mittelwert;
						wertcnt++;
					}
					else
					if ((code==5)||(code==6))
					{
						if (wertcnt<numsubb)
						{
							if ((int)(zone[i].zonenwert)==listsubb[wertcnt+1])
							{
								wertcnt++;
								buffer[bufcnt].wert[wertcnt]=mittelwert;
								lastresults.wert[wertcnt]=mittelwert;
							}
						}
					}

					//			if (wertcnt<numsubb)
					//              if (i==listsubb[wertcnt+1])
					//              {
					//                wertcnt++;
					//                buffer[bufcnt].wert[wertcnt]=mittelwert;
					//                lastresults.wert[wertcnt]=mittelwert;
					//              }
				}
			}
			if (gesamtRelCelSum>0)
			{
				gesamtmittel=gesamtsumme/gesamtRelCelSum;
			}
			else
			{
				gesamtmittel=0;
			}

			if (intervalle>1)
			for (i=0;i<=zonmax;i++)
			{
				sumsum[i]        = zone[i].sx;
				sumn[i]          = zone[i].n;
				sumRelCellSum[i] = zone[i].RelCellSum;
			}

			/* nach festgelegter Zanhl von Intervallen werden Summen oder Mittelwerte rausgeschrieben.
		Wenn der Code = 1 oder 2 war, dann ist in der Initialisierung pro Zeitintervall zone[].n
		auf 0 gesetzt worden, so dass zone[].sx/zone[].n immer den Mittelwert fuer das gesamte Teil-
		gebiet der Summen der Werte bis zum aktuellen Intervall enthaelt.
		Ist Code = oder 4, so enthaelt zone[].n immer die Anzahl der Zellen des Teil-Gebietes mal Anzahl
		bisher durchlaufener Intervalle, so dass zone[].sx/zone[].n immer den fortlaufenden Mittelwert
		fuer das Teilgebiet enthaelt, gesamtmittel also den fortlaufenden Gesamtmittelwert fuer das
		Gesamtgebiet */
			if (cnt==intervalle)
			{
				wertcnt++;
				buffer[bufcnt].wert[wertcnt] = gesamtmittel;
				lastresults.wert[wertcnt]    = gesamtmittel;
				/* jetzt koennen alle Statistikregister auf 0 gesetzt werden */
				if (intervalle>1)
				for (i=0;i<=zonmax;i++)
				{
					sumsum[i]        = 0.0;
					sumn[i]          = 0.0;
					sumRelCellSum[i] = 0.0;
				}
				cnt=0;
				if (bufcnt>=bufmax) flushbuffer();
			}
		}
	}
}


int zonalstatistic::schreibzonenwerte(const char *titel,const char *instanzname,TLastResults *LastResults,int maxLastResults)
{
	/* schreibt in die zonestat-Datei den Dateikopf (Reihenfolge der Zonen anhand ihrer Werte) */
	int i;
	int lastresultsindex;

	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	int mpi_rank=0;

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#endif

	if ((code>0)&&(mpi_rank==0))  /* code==0: keine Ausgabedatei */
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		//    (*WSstream) << "Ctrl-C disabled (writing file headers..) ");
		/* Schreiben der Zonenfolge fuer die Statistikdatei */
		fileopen("wt");
		for (i=0;((LastResults[i].wert)&&(i<maxLastResults));i++); // get the first free index of the results-array
		if ((i<=maxLastResults-1)&&(!LastResults[i].wert)) 
		{
			if (lastresults.wert) // the internal storage array for the last results is valid (wert is initialized)
			{
				LastResults[i] = lastresults;
				lastresultsindex=i;
			}
		}

#ifdef ASGI
		fprintf(zonestat,"%s_submodel:%s\n",titel,instanzname);
		sprintf(lastresults.Name,"%s_submodel:%s",titel,instanzname); // print the name also into the LastResults-array
		fprintf(zonestat,"YY\tMM\tDD\tHH");
		if (useMinutes) fprintf(zonestat, "\tmm");
		if ((code==2)||(code==4))  /* alle Teilgebiete werden beruecksichtigt */
		{
			for (i=1;i<=zonmax;i++)
			{
				fprintf(zonestat,"\t%.0lf",zone[i].zonenwert);
				lastresults.wert[i]=0;
				lastresults.refvalue[i]=0;
				lastresults.ZoneCodes[i]=(double)zone[i].zonenwert;
			}
		}
		if ((code==5)||(code==6))  /* nur eine definierte Auswahl */
		{
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++)  
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[zarray[i]].zonenwert);
				lastresults.wert[i]=0;
				lastresults.refvalue[i]=0;
				lastresults.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
			}
		}
#if defined _ENGLISH_
		fprintf(zonestat,"\ttot_average\n");
#else
		fprintf(zonestat,"\tg.Mittel\n");
#endif

#else

		fprintf(zonestat,"YY\tMM\tDD\tHH");
		if (useMinutes) fprintf(zonestat, "\tmm");
		if ((code==2)||(code==4))  /* alle Teilgebiete werden beruecksichtigt */
		{
			for (i=1;i<=zonmax;i++)  
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[i].zonenwert);
				lastresults.wert[i]=0;
				lastresults.refvalue[i]=0;
				lastresults.ZoneCodes[i]=(double)zone[i].zonenwert;
			}
		}
		if ((code==5)||(code==6))  /* nur eine definierte Auswahl */
		{
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++)  
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[zarray[i]].zonenwert);
				lastresults.wert[i]=0;
				lastresults.refvalue[i]=0;
				lastresults.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
			}
		}
#if defined _ENGLISH_
		fprintf(zonestat,"\ttot_average\n");
#else
		fprintf(zonestat,"\tg.Mittel\n");
#endif
		fprintf(zonestat,"%s_submodel:%s\n",titel,instanzname);
		sprintf(lastresults.Name,"%s_submodel:%s",titel,instanzname); // print the name also into the LastResults-array
#endif

		fprintf(zonestat,"--\t--\t--\t--");
		if (useMinutes) fprintf(zonestat, "\t--");
		if ((code==2)||(code==4))	{
			for (i=1;i<=zonmax;i++)  
			{
				fprintf(zonestat,"\t%.7lf",zone[i].rel);
				lastresults.rel[i]=zone[i].rel;
			}
			lastresults.wert[zonmax+1]=0;
			lastresults.refvalue[zonmax+1]=0;
			lastresults.ZoneCodes[zonmax+1]=0;
			lastresults.rel[zonmax+1]=(double)1.0;
		}	else if ((code==5)||(code==6)) {
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++) 
			{
				fprintf(zonestat,"\t%.7lf",zone[zarray[i]].rel);
				lastresults.rel[i]=zone[zarray[i]].rel;
			}
			lastresults.wert[numsubb+1]=0;
			lastresults.refvalue[numsubb+1]=0;
			lastresults.ZoneCodes[numsubb+1]=0;
			lastresults.rel[numsubb+1]=(double)1.0;
		} else {
			lastresults.wert[1]=0;
			lastresults.refvalue[1]=0;
			lastresults.ZoneCodes[1]=0;
			lastresults.rel[1]=(double)1.0;
		}
		fprintf(zonestat,"\t1.0000\n");
		// zonecode and initial value for the entire basin
		fileclose();
		SetSignalDefaults();
		return lastresultsindex;
	}
	else
	return -1;
}

zonalstatistic::zonalstatistic(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double **p_RelCellArea,int p_nrows,int p_ncols,statrecord **p_zone,
int p_zonmax,int p_code,int p_prec,char *p_dateiname,char *p_DefaultOutputDirectory,int p_numsubb,int *p_listsubb, bool useMinutes):
initialized(0),WSstream(0), zone(0),buffer(0),bufmax(0),bufcnt(0),
wertcnt(0),RelCellArea(0),nrows(0),ncols(0),nodata(0),zonmax(0),
prec(0),code(0),intervalle(0),numsubb(0),listsubb(0),dateiname(0),DefaultOutputDirectory(0),
zonestat(0),mittelwert(0),cnt(0),gesamtsumme(0),zarray(0),useMinutes(false)
{
	try
	{
		retvalue=0;
		int i;
		unsigned long li;
		WSstream=Stream;
		if (p_code>1000)
		{
			RelCellArea = p_RelCellArea;
			nrows       = p_nrows;
			ncols       = p_ncols;
			nodata      = -9999;
			code        = p_code / 1000; /* von 1nnn ... 4nnn bleibt 1..4 */
			intervalle  = p_code % 1000; /* von 1nnn ... 4nnn bleibt nnn */
			zonmax      = p_zonmax;
			zone        = *p_zone;
			prec        = p_prec;
			dateiname   = p_dateiname;
			DefaultOutputDirectory=p_DefaultOutputDirectory;
			numsubb     = p_numsubb;
			if (numsubb > zonmax) numsubb=zonmax; /* no more zones are allowed for the output-list than are defined in the zonegrid */
			listsubb    = p_listsubb;
			gesamtsumme = 0.0;  /* Summe der pro Intervall errechneten Gebietsmittelwerte fuer cnt Intervalle */
			cnt         = 0;    /* haelt die Anzahl der abgelaufenen Intervalle zur Mittelbildung fest */
			this->useMinutes = useMinutes;
			// construct filename when a default output dîr was given
			if (DefaultOutputDirectory)
			{
				bool bIsComplete=false;
				std::string temp(dateiname);

#if defined(WIN32)
				if ((temp.find(":\\") == 1) ||  // like c:\... in traditional notation
					(temp.find(":/") == 1) ||  // like c:/... in alternative notation
					(temp.find("\\\\")==0)||  // like \\server\... in traditional notation
					(temp.find("\\\\")==0))   // like //server/... in alternative notation
					bIsComplete = true;
#elif defined(WIN64)
				if ((temp.find(":\\") == 1) ||  // like c:\... in traditional notation
					(temp.find(":/") == 1) ||  // like c:/... in alternative notation
					(temp.find("\\\\") == 0) ||  // like \\server\... in traditional notation
					(temp.find("\\\\") == 0))   // like //server/... in alternative notation
					bIsComplete = true;
#else
				if (temp.find_first_of("/~") == 0) bIsComplete = true; // starts with a / or a ~ meaning there are absolute paths defined under UNIX/LINUX
#endif
				
				if (!bIsComplete)
				{
					char *completename=new char[_stringsize_+1];
					sprintf(completename,"%s%s",DefaultOutputDirectory,dateiname);
					sprintf(dateiname,"%s",completename);
					delete[] completename;
				}
			}

			if (numsubb>0)
			{
				zarray = new int[numsubb+1];
				for (int i=0;i<=numsubb;i++)
				{
					zarray[i]=0;
				}
				for (int j=1;j<=numsubb;j++)
				{
					// suchen des zone[x].zonenwertes, der dem listsubb-EIntrag entspricht
					for (int k=j;k<=zonmax;k++)
					{
						if ((int)zone[k].zonenwert==listsubb[j])
						{
							zarray[j]=k;
							break;
						}
						if (k>=zonmax)
						{
							(*WSstream) << "\n\n ------ error ------\n\n there is a code in the output list which no zone exists for (" << listsubb[j] <<")"<<"[end-entry]";
							throw 1;
						}
					}
				}
			}
			else
			zarray=NULL;

			if ((code==1)||(code==3))
			bufmax = (unsigned)(buffersize/(4*sizeof(long)+2*sizeof(double)));
			else
			{
				if ((code==2)||(code==4))
				bufmax = (unsigned)(buffersize/(4*sizeof(long)+(zonmax+1)*sizeof(double)));
				else
				bufmax = (unsigned)(buffersize/(4*sizeof(long)+(numsubb+1)*sizeof(double)));
			}
			if (bufmax<1) bufmax=1;
			bufcnt      = 0; /* initialising the buffer index counter with 0 */
			buffer      = (buffertype*)calloc((size_t)bufmax+1,sizeof(buffertype));
			titel       = (char*)calloc(_stringsize_,sizeof(char));
			if ((code==1)||(code==3))
			{
				for (li=1;li<=bufmax;li++) buffer[li].wert=(double*)calloc(2,sizeof(double)); /* allocating memory for the mean + the zero-index-element */
				lastresults.wert = (double*)calloc(2,sizeof(double));
				lastresults.refvalue = (double*)calloc(2,sizeof(double));
				lastresults.rel  = (double*)calloc(2,sizeof(double));
				lastresults.ZoneCodes = (unsigned long*)calloc(2,sizeof(unsigned long));
				lastresults.NumZones  = 2;
				lastresults.Name = titel;
				lastresults.j = 0;
				lastresults.m = 0;
				lastresults.t = 0;
				lastresults.s = 0;
				lastresults.min = 0;
			}
			else
			{
				if ((code==2)||(code==4))
				{
					for (li=1;li<=bufmax;li++) buffer[li].wert=(double*)calloc(zonmax+2,sizeof(double)); /* allocating memory for each zone and the mean + the zero-index-element */
					lastresults.wert = (double*)calloc(zonmax+2,sizeof(double));
					lastresults.refvalue = (double*)calloc(zonmax+2,sizeof(double));
					lastresults.rel  = (double*)calloc(zonmax+2,sizeof(double));
					lastresults.ZoneCodes = (unsigned long*)calloc(zonmax+2,sizeof(unsigned long));
					lastresults.NumZones  = zonmax+1;
					lastresults.Name = titel;
					lastresults.j = 0;
					lastresults.m = 0;
					lastresults.t = 0;
					lastresults.s = 0;
					lastresults.min = 0;
				}
				else
				{
					for (li=1;li<=bufmax;li++) buffer[li].wert=(double*)calloc(numsubb+2,sizeof(double)); /* allocating memory for each zone and the mean + the zero-index-element */
					lastresults.wert = (double*)calloc(numsubb+2,sizeof(double));
					lastresults.refvalue = (double*)calloc(numsubb+2,sizeof(double));
					lastresults.rel  = (double*)calloc(numsubb+2,sizeof(double));
					lastresults.ZoneCodes = (unsigned long*)calloc(numsubb+2,sizeof(unsigned long));
					lastresults.NumZones  = numsubb+1;
					lastresults.Name = titel;
					lastresults.j = 0;
					lastresults.m = 0;
					lastresults.t = 0;
					lastresults.s = 0;
					lastresults.min = 0;
				}
			}
			for (i=0;i<=zonmax;i++)
			{
				zone[i].n          = 0.0;
				zone[i].sx         = 0.0;
				zone[i].RelCellSum = 0.0;
			}
			if (intervalle>1) /* Reservieren von Speicher zum Zwischenspeichern der Ergebnisse der STatistik,
						wenn Ergebnisse ueber mehrere Intervalle gemittelt oder summiert werden sollen */
			{
				if ((sumsum = (double*)calloc(zonmax+1,sizeof(double)))==NULL)
				{
					(*WSstream) << "FATAL: no free memory in zonalstatistik::zonalstatistik "<<"[end-entry]";
					throw 1;
				}
				if ((sumn   = (double*)calloc(zonmax+1,sizeof(double)))==NULL)
				{
					(*WSstream) << "FATAL: no free memory in zonalstatistik::zonalstatistik "<<"[end-entry]";
					throw 1;
				}
				if ((sumRelCellSum = (double*)calloc(zonmax+1,sizeof(double)))==NULL)
				{
					(*WSstream) << "FATAL: no free memory in zonalstatistik::zonalstatistik "<<"[end-entry]";
					throw 1;
				}
				for (i=0;i<=zonmax;i++)
				{
					sumsum[i]        = 0.0;
					sumn[i]          = 0.0;
					sumRelCellSum[i] = 0.0;
				}
			}
		}
		else code=-1;
		initialized=1; // this is a flag wich is initially set to 0 within the baseclass constructor (memory_handling) per default.
		// it is used for determining the kind of destruction of the members (resp. the kind of ClearMemory), 
		// because the "destructors" are also called from the DLL-Framework, if an error during initialization occurs
		// its neccesary to distinct between normal destruction and destruction which should only free memory (no buffer writing and so on)
	}
	catch (...)
	{ // all exceptions should be handled within this constructor in order to allow the complete construction
		// of the object. Zhis is required in order to call the destructor correctly (especially, if errors during
		// construction occure)
		retvalue=1;
		initialized=0;
	}
}

void zonalstatistic::ClearMemory()
{
	unsigned long i;
	//  if (!initialized) return; // no buffers should be written if an error during initialization occurs
	if (code>0)
	{
		if (bufcnt>0) 
		if (initialized) flushbuffer(); /* if some data is still in the internal buffer */
		for (i=1;i<=bufmax;i++) 
		{
			if (buffer[i].wert) free(buffer[i].wert);
			buffer[i].wert=NULL;
		}
		if (buffer) free(buffer); buffer=NULL;
		if (lastresults.wert) free(lastresults.wert); lastresults.wert = NULL;
		if (lastresults.refvalue) free(lastresults.refvalue); lastresults.refvalue = NULL;
		if (lastresults.rel ) free(lastresults.rel ); lastresults.rel  = NULL;
		if (lastresults.ZoneCodes) free(lastresults.ZoneCodes); lastresults.ZoneCodes = NULL;
		if (titel) free(titel); titel=NULL;
		if (zarray) delete[] zarray; zarray=NULL;
		if ((code>0)&&(intervalle>1))
		{
			if (sumsum) free(sumsum);sumsum = NULL;
			if (sumn) free(sumn);sumn=NULL;
			if (sumRelCellSum) free(sumRelCellSum);sumRelCellSum=NULL;
		}
	}
}

/* ~~~~~~~~~~~~~~~~~~ Ende der Routinen fuer Statistikfiles ~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~ Routinen fuer Abflussrouting : abgeleitet von normaler zonalstatistik ~~~~~~~~~~~~~~~ */

routstatistic::routstatistic(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,
int maxLastResults,double **p_RelCellArea,int p_nrows,int p_ncols,
statrecord **p_zone,double *p_abfrel,int p_zonmax,int p_code,int p_prec,
char *p_dateiname,char *p_DefaultOutputDirectory,int p_numsubb,int *p_listsubb, bool useMinutes):
zonalstatistic(retvalue,Stream,LastResults,maxLastResults,p_RelCellArea,p_nrows,p_ncols,p_zone,
p_zonmax,p_code,p_prec,p_dateiname,p_DefaultOutputDirectory,p_numsubb,p_listsubb,useMinutes),
abfrel(0)
{

	if (p_code>1000)
	{
		if ((code==1)||(code==3))
		{
			lastresults.NumZones--;
			observedRunoff.wert = (double*)calloc(2,sizeof(double));
			observedRunoff.refvalue = (double*)calloc(2,sizeof(double));
			observedRunoff.rel  = (double*)calloc(2,sizeof(double));
			observedRunoff.ZoneCodes = (unsigned long*)calloc(2,sizeof(unsigned long));
			observedRunoff.NumZones  = 1;
			observedRunoff.Name = (char*)calloc(500,sizeof(char));
			sprintf(observedRunoff.Name,"observed Runoff");
			observedRunoff.j = 0;
			observedRunoff.m = 0;
			observedRunoff.t = 0;
			observedRunoff.s = 0;
			observedRunoff.min = 0;
		}
		else
		{
			if ((code==2)||(code==4))
			{
				lastresults.NumZones--;
				observedRunoff.wert = (double*)calloc(zonmax+2,sizeof(double));
				observedRunoff.refvalue = (double*)calloc(zonmax+2,sizeof(double));
				observedRunoff.rel  = (double*)calloc(zonmax+2,sizeof(double));
				observedRunoff.ZoneCodes = (unsigned long*)calloc(zonmax+2,sizeof(unsigned long));
				observedRunoff.NumZones  = zonmax;
				observedRunoff.Name = (char*)calloc(500,sizeof(char));
				sprintf(observedRunoff.Name,"observed Runoff");
				observedRunoff.j = 0;
				observedRunoff.m = 0;
				observedRunoff.t = 0;
				observedRunoff.s = 0;
				observedRunoff.min = 0;
			}
			else
			{
				lastresults.NumZones--;
				observedRunoff.wert = (double*)calloc(numsubb+2,sizeof(double));
				observedRunoff.refvalue = (double*)calloc(numsubb+2,sizeof(double));
				observedRunoff.rel  = (double*)calloc(numsubb+2,sizeof(double));
				observedRunoff.ZoneCodes = (unsigned long*)calloc(numsubb+2,sizeof(unsigned long));
				observedRunoff.NumZones  = numsubb;
				observedRunoff.Name = (char*)calloc(500,sizeof(char));
				sprintf(observedRunoff.Name,"observed Runoff");
				observedRunoff.j = 0;
				observedRunoff.m = 0;
				observedRunoff.t = 0;
				observedRunoff.s = 0;
				observedRunoff.min = 0;
			}
		}
	}
	abfrel = p_abfrel;
}

void routstatistic::calc_and_done(long jahr,long monat,long tag,long stunde,long minute)
/* Schreibt die zonalen Ergebnisse in die Statistikdateien, Dezimalstellen p */
{
	int i;
	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	if (code>0)
	{
		if (intervalle>1)
		for (i=1;i<=zonmax;i++)
		{
			zone[i].sx+=sumsum[i]; /* updatting fuer laufende Summen ueber mehrere Intervalle */
			if ((code==1)||(code==2)||(code==5))
			{
				zone[i].n+=sumn[i]; /* wenn Mittelwerte gebildet werden und nicht Summen ueber mehrere Intervalle, muss n von 1 an gezaehlt werden */
				zone[i].RelCellSum+=sumRelCellSum[i];
			}
		}

		if (cnt==intervalle)
		{

			wertcnt=0;
			bufcnt++;
			//      buffer[bufcnt].j=jahr%100;
			buffer[bufcnt].j = jahr;
			buffer[bufcnt].m = monat;
			buffer[bufcnt].t = tag;
			buffer[bufcnt].s = stunde;
			buffer[bufcnt].min = minute;
			observedRunoff.j = lastresults.j = jahr;
			observedRunoff.m = lastresults.m = monat;
			observedRunoff.t = lastresults.t = tag;
			observedRunoff.s = lastresults.s = stunde;
			observedRunoff.min = lastresults.s = minute;

			for (i=1;i<=zonmax;i++)
			{
				/* sollen die Zonenmittelwerte jeweils ueber mehrere Intervalle addiert werden,
		so wird das durch die jeweils am Beginn eines Intervalles erfolgende Nullsetzen
		des n-Zaehlers im Zonenrecord beruecksichtigt (in zonalstatistic::init()).
		Bei Mittelwertbildung ueber mehrere Intervalle wird zone[i].n nicht auf 0 gesetzt,
		das erfolgt erst nach dem Schreiben nach der entsprechenden Anzahl Intervallen */

				//        if (zone[i].n>0) mittelwert=(zone[i].sx/zone[i].n)*(zone[i].n/zone[i].RelCellSum); else mittelwert = 0.0;
				if (zone[i].n>0) mittelwert=(zone[i].sx/zone[i].n)*(zone[i].n/zone[i].RelCellSum); else mittelwert = 0.0;
				if ((code==2)||(code==4))
				{
					buffer[bufcnt].wert[i]=mittelwert;
					lastresults.wert[i]=mittelwert;
					wertcnt++;
				}
				else
				if ((code==5)||(code==6))
				if (wertcnt<numsubb)
				{
					if ((int)(zone[i].zonenwert)==listsubb[wertcnt+1])
					{
						wertcnt++;
						buffer[bufcnt].wert[wertcnt]=mittelwert;
						lastresults.wert[wertcnt]=mittelwert;
					}
				}
				//            if (wertcnt<numsubb)
				//              if (i==listsubb[wertcnt+1])
				//              {
				//                wertcnt++;
				//                buffer[bufcnt].wert[wertcnt] = mittelwert;
				//                lastresults.wert[wertcnt]    = mittelwert;
				//              }
			}

			/* jetzt koennen alle Statistikregister auf 0 gesetzt werden */
			if (intervalle>1)
			for (i=0;i<=zonmax;i++)
			{
				sumsum[i]        = 0.0;
				sumn[i]          = 0.0;
				sumRelCellSum[i] = 0.0;
			}
			cnt=0;
			if (bufcnt>=bufmax) flushbuffer();
		}
		else
		{
			if (intervalle>1)
			for (i=0;i<=zonmax;i++)
			{
				sumsum[i]        = zone[i].sx;
				sumn[i]          = zone[i].n;
				sumRelCellSum[i] = zone[i].RelCellSum;
			}
		}
	}
}

int routstatistic::schreibzonenwerte(const char *titel,const char *instanzname,TLastResults *LastResults,int maxLastResults)
{
	/* schreibt in die zonestat-Datei den Dateikopf (Reihenfolge der Zonen anhand ihrer Werte) */
	int i;
	int lastresultsindex;
	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	/* oeffnen der Statistikdatei und Schreiben des Dateikopfes */
	if (code>0)  /* code==0: keine Ausgabedatei */
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		//    (*WSstream) << "Ctrl-C disabled (writing file headers..) ");
		/* Schreiben der Zonenfolge fuer die Statistikdatei */
		fileopen("wt");
		for (i=0;((LastResults[i].wert)&&(i<maxLastResults));i++); // get the first free index of the results-array
		if ((i<=maxLastResults-1)&&(!LastResults[i].wert)) 
		{
			if (lastresults.wert) // the internal storage array for the last results is valid (wert is initialized)
			{
				LastResults[i] = lastresults;
			}
		}
		i++; // now the observed Runoff
		if ((i<=maxLastResults-1)&&(!LastResults[i].wert)) 
		{
			if (observedRunoff.wert) // the internal storage array for the last results is valid (wert is initialized)
			{
				LastResults[i] = observedRunoff;
				lastresultsindex=i;
			}
		}
#ifdef ASGI

		/* Unter ASGI wird zuerst der Titel, dann die TEZ ohne "" erwartet */
		fprintf(zonestat,"%s_submodel:%s\n",titel,instanzname);
		sprintf(lastresults.Name,"%s_submodel:%s",titel,instanzname); // print the name also into the LastResults-array
		fprintf(zonestat,"YY\tMM\tDD\tHH");
		if (useMinutes) fprintf(zonestat, "\tmm");
		if ((code==2)||(code==4))
		{
			for (i=1;i<=zonmax;i++)  
			{
				fprintf(zonestat,"\t%.0lf",zone[i].zonenwert);
				lastresults.wert[i]=0;
				lastresults.ZoneCodes[i]=(double)zone[i].zonenwert;
				observedRunoff.wert[i]=0;
				observedRunoff.ZoneCodes[i]=(double)zone[i].zonenwert;
			}
		}
		if ((code==5)||(code==6))
		{
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++) 
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[zarray[i]].zonenwert);
				lastresults.wert[i]=0;
				lastresults.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
				observedRunoff.wert[i]=0;
				observedRunoff.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
			}
		}
		fprintf(zonestat,"\n");

#else

		fprintf(zonestat,"YY\tMM\tDD\tHH");
		if (useMinutes) fprintf(zonestat, "\tmm");
		if ((code==2)||(code==4))
		{
			for (i=1;i<=zonmax;i++)  
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[i].zonenwert);
				lastresults.wert[i]=0;
				lastresults.ZoneCodes[i]=(double)zone[i].zonenwert;
				observedRunoff.wert[i]=0;
				observedRunoff.ZoneCodes[i]=(double)zone[i].zonenwert;
			}
		}
		if ((code==5)||(code==6))
		{
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++) 
			{
				fprintf(zonestat,"\t\"%.0lf\"",zone[zarray[i]].zonenwert);
				lastresults.wert[i]=0;
				lastresults.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
				observedRunoff.wert[i]=0;
				observedRunoff.ZoneCodes[i]=(double)(zone[zarray[i]].zonenwert);
			}
		}
		fprintf(zonestat,"\n");
		fprintf(zonestat,"%s_submodel:%s\n",titel,instanzname);
		sprintf(lastresults.Name, "%s_submodel:%s", titel, instanzname); // print the name also into the LastResults-array
#endif

		fprintf(zonestat,"--\t--\t--\t--");
		if (useMinutes) fprintf(zonestat, "\t--");
		if ((code==2)||(code==4))
		{
			for (i=1;i<=zonmax;i++)  
			{
				fprintf(zonestat,"\t%.7lf",abfrel[i]);
				lastresults.rel[i]=abfrel[i];
				observedRunoff.rel[i]=abfrel[i];
			}
		}
		if ((code==5)||(code==6))
		{
			for (i=1;(i<=numsubb)&&(i<=zonmax);i++) 
			{
				fprintf(zonestat,"\t%.7lf",abfrel[zarray[i]]);
				lastresults.rel[i]=abfrel[zarray[i]];
				observedRunoff.rel[i]=abfrel[zarray[i]];
			}
		}
		fprintf(zonestat,"\n");
		fileclose();
		SetSignalDefaults();
		return lastresultsindex;
	}
	else
	return 0;
}

void routstatistic::Rsquare(Tstatist **p_eff,Tstatist **p_leff)
{
	int    i,j,wertcnt;
	double oristreu;
	/* code==0:    keine Ausgabe
	code==1nnn: pro nnn Intervalle Mittel nur fuer Gesamtgebiet
	code==2nnn: pro nnn Intervalle Mittel     fuer alle Zonen
	code==3nnn: pro nnn Intervalle Summe  nur fuer Gesamtgebiet
	code==4nnn: pro nnn Intervalle Summe      fuer alle Zonen
*/
	if (code>0)
	{
		flushbuffer(); /* if there are still any data in the buffer which must be written into the file NOW */
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		//    (*WSstream) << "Ctrl-C disabled (writing r-squares..) ");
		fileopen("a+");
		fprintf(zonestat,"----------------------\n");
		/* zunaechst lineares rsquare */
		for (i=0;i<=20;i++)
		{
			if (i!=10) fprintf(zonestat,"lin. R-square\t%d",i-10); else fprintf(zonestat,"LIN. R-SQUARE\t%d",i-10);
			wertcnt=0;
			for (j=1;j<=zonmax;j++)
			{
				/* Streuung der gemessenen Daten */
				if (p_eff[i][j].ngem > 2)
				oristreu =  (p_eff[i][j].sx2gem  - (p_eff[i][j].sxgem  * p_eff[i][j].sxgem  ) / p_eff[i][j].ngem ) / (p_eff[i][j].ngem-1 );
				else
				oristreu=0.0;
				/* Streuung der berechneten Residuen */
				if (p_eff[i][j].ndiff > 2)
				mittelwert=(p_eff[i][j].sx2diff) / (p_eff[i][j].ndiff-1);
				else
				mittelwert=0.0;
				if ((oristreu==0)||(mittelwert==0)) mittelwert=1.0; else mittelwert = 1 - mittelwert/oristreu;
				if ((code==2)||(code==4))
				{
					fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
				}
				else
				{
					if ((code==5)||(code==6))
					{
						if (wertcnt<numsubb)
						{
							if (j==zarray[wertcnt+1])
							{
								wertcnt++;
								fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
							}
						}
					}
				}
			}
			fprintf(zonestat,"\n");
		}
		fprintf(zonestat,"----------------------\n");
		/* nun logarithmisches rsquare */
		for (i=0;i<=20;i++)
		{
			if (i!=10) fprintf(zonestat,"log. R-square\t%d",i-10); else fprintf(zonestat,"LOG. R-SQUARE\t%d",i-10);
			wertcnt=0;
			for (j=1;j<=zonmax;j++)
			{
				/* Streuung der gemessenen Daten */
				if (p_leff[i][j].ngem > 2)
				oristreu =  (p_leff[i][j].sx2gem  - (p_leff[i][j].sxgem  * p_leff[i][j].sxgem  ) / p_leff[i][j].ngem ) / (p_leff[i][j].ngem-1 );
				else
				oristreu=0.0;
				/* Streuung der berechneten Residuen */
				if (p_leff[i][j].ndiff > 2)
				mittelwert=(p_leff[i][j].sx2diff) / (p_leff[i][j].ndiff-1);
				else
				mittelwert=0.0;
				if ((oristreu==0)||(mittelwert==0)) mittelwert=1.0; else mittelwert = 1 - mittelwert/oristreu;
				if ((code==2)||(code==4))
				fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
				else
				{
					if ((code==5)||(code==6))
					{
						if (wertcnt<numsubb)
						{
							if (j==zarray[wertcnt+1])
							{
								wertcnt++;
								fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
							}
						}
					}
				}
			}
			fprintf(zonestat,"\n");
		}
		fprintf(zonestat,"----------------------\n");
		/* jetzt lineare EV (Erklaerte Varianz) */
		for (i=0;i<=20;i++)
		{
#if defined _ENGLISH_
			if (i!=10) fprintf(zonestat,"coeff.var.lin\t%d",i-10); else fprintf(zonestat,"COEFF.VAR.LIN\t%d",i-10);
#else
			if (i!=10) fprintf(zonestat,"erkl. Var.lin\t%d",i-10); else fprintf(zonestat,"ERKL. VAR.LIN\t%d",i-10);
#endif
			wertcnt=0;
			for (j=1;j<=zonmax;j++)
			{
				/* Streuung der gemessenen Daten */
				if (p_eff[i][j].ngem > 2)
				oristreu =  (p_eff[i][j].sx2gem  - (p_eff[i][j].sxgem  * p_eff[i][j].sxgem  ) / p_eff[i][j].ngem ) / (p_eff[i][j].ngem-1 );
				else
				oristreu=0.0;
				/* Streuung der berechneten Residuen jeweils um den mittleren Fehler bereinigt (Verschiebung) */
				if (p_eff[i][j].ndiff > 2)
				{
					mittelwert =(p_eff[i][j].sx2diff - (p_eff[i][j].sxdiff * p_eff[i][j].sxdiff ) / p_eff[i][j].ndiff) / (p_eff[i][j].ndiff-1);
				}
				else
				mittelwert=0.0;
				if ((oristreu==0)||(mittelwert==0)) mittelwert=1.0; else mittelwert = 1 - mittelwert/oristreu;
				if ((code==2)||(code==4))
				fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
				else
				{
					if ((code==5)||(code==6))
					{
						if (wertcnt<numsubb)
						{
							if (j==zarray[wertcnt+1])
							{
								wertcnt++;
								fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
							}
						}
					}
				}
				if (i==10)
				{
					zone[j].sx=mittelwert;
					zone[j].n=1;
				}
			}
			fprintf(zonestat,"\n");
		}
		fprintf(zonestat,"----------------------\n");
		/* nun logarithmisches EV */
		for (i=0;i<=20;i++)
		{
#if defined _ENGLISH_
			if (i!=10) fprintf(zonestat,"coeff.var.log\t%d",i-10); else fprintf(zonestat,"COEFF.VAR.LOG\t%d",i-10);
#else
			if (i!=10) fprintf(zonestat,"erkl. Var.log\t%d",i-10); else fprintf(zonestat,"ERKL. VAR.LOG\t%d",i-10);
#endif
			wertcnt=0;
			for (j=1;j<=zonmax;j++)
			{
				/* Streuung der gemessenen Daten */
				if (p_leff[i][j].ngem > 2)
				oristreu =  (p_leff[i][j].sx2gem  - (p_leff[i][j].sxgem  * p_leff[i][j].sxgem  ) / p_leff[i][j].ngem ) / (p_leff[i][j].ngem-1 );
				else
				oristreu=0.0;
				/* Streuung der berechneten Residuen jeweils um den mittleren Fehler vermindert*/
				if (p_leff[i][j].ndiff > 2)
				mittelwert=(p_leff[i][j].sx2diff - (p_leff[i][j].sxdiff * p_leff[i][j].sxdiff ) / p_leff[i][j].ndiff) / (p_leff[i][j].ndiff-1);
				else
				mittelwert=0.0;
				if ((oristreu==0)||(mittelwert==0)) mittelwert=1.0; else mittelwert = 1 - mittelwert/oristreu;
				if ((code==2)||(code==4))
				fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
				else
				{
					if ((code==5)||(code==6))
					{
						if (wertcnt<numsubb)
						{
							if (j==zarray[wertcnt+1])
							{
								wertcnt++;
								fprintf(zonestat,"\t%.*lf",prec-1,mittelwert);
							}
						}
					}
				}
			}
			fprintf(zonestat,"\n");
		}
		fprintf(zonestat,"----------------------\n");
		fileclose();
		SetSignalDefaults();
	}
}

void routstatistic::ClearMemory()
{
	//  if (!initialized) return; // no buffers should be written if an error during initialization occurs
	if (code>0)
	{
		if (observedRunoff.wert) free(observedRunoff.wert); observedRunoff.wert = NULL;
		if (observedRunoff.refvalue) free(observedRunoff.refvalue); observedRunoff.refvalue = NULL;
		if (observedRunoff.rel ) free(observedRunoff.rel ); observedRunoff.rel  = NULL;
		if (observedRunoff.ZoneCodes) free(observedRunoff.ZoneCodes); observedRunoff.ZoneCodes = NULL;
		if (observedRunoff.Name) free(observedRunoff.Name); observedRunoff.Name = NULL;
	}
	zonalstatistic::ClearMemory();
}

/* ~~~~~~~~~~~~~~~~~~ Ende Routinen fuer Abflussrouting ~~~~~~~~~~~~~~~ */

void modelobject::steuererror(const char *text,const char *message)
{
	(*WSstream) << "FATAL: error in control file, "<<text<<"; "<<message<<"not found!"<<"[end-entry]";
	throw 17;
}

void modelobject::writeerror(char *dateiname,int fehler)
{
	(*WSstream) << "FATAL: write error "<<dateiname<<" ! ("<<fehler<<")"<<"[end-entry]";
	throw 18;
}

void modelobject::readerror(char *dateiname,int fehler)
{
	(*WSstream) << "FATAL: read error "<<dateiname<<" ! ("<<fehler<<")"<<"[end-entry]";
	throw 19;
}


int  modelobject::fscanfln(FILE *datei,const char *format,void *var,const char *errorstring,char *string)
{     
	int i;
	if ((i=fscanf(datei,format,var)) != 1)
	{          
		// ARNOLD 
		// --> START CHANGE  
		//       (make output, and read one additional line to indicate where error occurred!)
		(*WSstream) << "--> In ." << "\n"
		//					<< " Status of read statement: " << i << "\n"
		//					<< "     Format of input read: " << format << "\n"
		//					<< "                   string: " << *string << "\n"
		//					<< "                     &val: " << &var << " or direct as var: " << var<< "\n\n"
		<< "  Next line in Wasim control file: " << "[end-entry]";
		char c[1000];
		fscanf(datei,"%s",c);
		(*WSstream) << c << "[end-entry]";
		// <-- END CHANGE
		steuererror(errorstring,instanzname);
	}

	fgets(string,_stringsize_,datei);
	return i;
}

/*
int  modelobject::fscanfln(FILE *datei,const char *format,void *var,const char *errorstring,char *string)
{
int i;
if ((i=fscanf(datei,format,var)) != 1) steuererror(errorstring,instanzname);
fgets(string,_stringsize_,datei);
return i;
}
*/

/* Ueberladen der function fscanfln mit einer gleichnamigen function mit 2 einzulesenden Variablen */
int  modelobject::fscanfln(FILE *datei,const char *format,void *var1,void *var2,const char *errorstring,char *string)
{
	int i;
	if ((i=fscanf(datei,format,var1,var2)) != 2) steuererror(errorstring,instanzname);
	fgets(string,_stringsize_,datei);
	return i;
}

void modelobject::SetInitialized(int p_initialized)
{
	initialized=p_initialized;
}

void memory_handling::speichererror(int fehler)
{
	(*WSstream) << "FATAL: memory allocation error in "<<instanzname<<" ("<<fehler<<"[end-entry]";
	throw 81;
}

/*
double **memory_handling::d_matrix(long nrl,long nrh,long ncl,long nch)
{
long i, nrow=nrh+1, ncol=nch+1;
double **m;

m = (double **)calloc((size_t)(nrow),sizeof(double*));
if (!m) { (*WSstream) << "Speicherbelegungsfehler 1 in matrix()\n"); throw(81); }
m[0]=(double *)calloc((size_t)(nrow*ncol),sizeof(double));
if (!m[0]) { (*WSstream) << "Speicherbelegungsfehler 2 in matrix()\n"); throw(82); }
for (i=1;i<=nrh;i++) m[(int)i]=m[(int)(i-1)]+(int)ncol;
return m;
}


void memory_handling::free_d_matrix(double **m,long nrl,long nrh,long ncl,long nch)
{
long i;
free((char *)(m[0]);
free((char *)(m+(int)nrl));
}
*/


//double **memory_handling::d_matrix(long nrl, long nrh, long ncl, long nch, long nrows)
//{
//	long i,j,nrow,ncol,temprow,colcount,overhead;
//	double size,tempsize;
//	int exponent;
//	double **m;
//
//	nrow = nrh - nrl + 2;    /* because of zero-index there should be one row more. normally, the number of rows to reserve would be like nrh-nrl+1 */
//	ncol = nch + 1;    /* since we always reserve at least an entire row, the number of columns to allocate will be the same for all grids, wether in MPI mode or in normal mode  */
//	temprow = nrl;        /* the row, up to which memory is allocated */
//	colcount = ncol*sizeof(double);                                                     /* numberof bytesfor one row */
//	if ((unsigned long)colcount+24>(unsigned long)_maxlinmem_) 
//	{
//		(*WSstream) << "FATAL: too long grid rows"<<"[end-entry]"; 
//		throw(1);
//	}
//	size = nrow*ncol*sizeof(double);                                                   /* how much memory is needed? */
//	overhead = (unsigned long)(size * _maxmemoverhead_);                /* whatis the overhead we are willing to waste? (say, 5%) */
//	if (nrows>0)
//		m = (double **)calloc((unsigned long)(nrows+1), sizeof(double*));       /* allocating memory for the pointersto the rows (row indizes) */
//	else
//		m = (double **)calloc((unsigned long)(nrh + 1), sizeof(double*));       /* allocating memory for the pointersto the rows (row indizes) */
//    	
//	if (!m) 
//	{ 
//		(*WSstream) << "FATAL: memory allocation error 1 in matrix()"<<"[end-entry]"; 
//		throw(1); 
//	}
//	for (int i = 0; i < temprow; i++)
//	{
//		m[(int)i] = NULL;
//	}
//	do  /* subdividing the size into powers of 2 */
//	{
//		exponent=(int)ceil(log_2(size+24.0));  /* exponent for the largest power of 2 that is smaller than size */
//		if ( ((size<=overhead)&&((unsigned long)(pow(2.0,(double)exponent)<_maxlinmem_)))
//				|| (pow(2.0,(double)exponent-1.0)<(colcount+24)))
//		{
//			/* this is that part of the loop, which is done at last (if not so much memory was allocated or if
//		the rest of memory to allocate is smaller than overhead) */
//			m[(int)temprow]=(double *)calloc((size_t)((nrow-temprow+nrl-1)*ncol),sizeof(double));
//			if (!m[(int)temprow]) 
//			{ 
//				(*WSstream) << "FATAL: memory allocation error 2 in matrix()"<<"[end-entry]"; 
//				throw(1); 
//			}
//			for (i=temprow+1;i<=nrh;i++) m[(int)i]=m[(int)(i-1)]+(int)ncol;
//			size=-1.0;
//		}
//		else
//		{
//			tempsize = pow(2.0,(double)(exponent-1))-24.0; /* how much memory can be allocated by this exponent */
//			if ((unsigned long)tempsize>_maxlinmem_) tempsize=(double)_maxlinmem_; /* if this is more than _maxlinmem_, cut it to that value */
//			j = (unsigned long)floor(tempsize/colcount);                           /* largest integer number of rows, allocatable with tempsize bytes*/
//			size = (double)(nrow-j-temprow+nrl-1)*sizeof(double)*ncol;          /* remaining size of memory to allocate */
//			m[(int)temprow]=(double *)calloc((size_t)(j*ncol),sizeof(double));  /* now allocating the memory for <j> rows */
//			if (!m[(int)temprow]) 
//			{ 
//				(*WSstream) << "FATAL: memory allocation error 2 in matrix()"<<"[end-entry]"; 
//				throw(1); 
//			}
//			for (i=temprow+1;i<=temprow+j-1;i++) m[(int)i]=m[(int)(i-1)]+(int)ncol; /* setting the pointers for row indizes */
//			temprow+=j; /* in the next lap: this is the row, for wich is now to allocate memory */
//		}
//	}
//	while (size > 0);
//	for (int i = nrh + 1; i <= nrows; i++)
//	{
//		m[(int)i] = NULL;
//	}
//	return m;
//}
//
//void memory_handling::free_d_matrix(double **m, long nrl, long nrh, long ncl, long nch, long nrows)
//{
//	long j,nrow,ncol,temprow,colcount,overhead;
//	double size,tempsize;
//	int exponent;
//
//	if (!m) return;
//
//	nrow = nrh - nrl + 2;    /* because of zero-index there should be one row more. normally, the number of rows to reserve would be like nrh-nrl+1 */
//	ncol = nch + 1;    /* since we always reserve at least an entire row, the number of columns to allocate will be the same for all grids, wether in MPI mode or in normal mode  */
//	temprow = nrl;        /* the row, up to which memory is allocated */
//	colcount = ncol*sizeof(double);
//	if ((unsigned long)colcount+24>(unsigned long)_maxlinmem_) 
//	{
//		(*WSstream) << "FATAL: too long matrix rows"<<"[end-entry]"; 
//		throw(1);
//	}
//	size = nrow*ncol*sizeof(double);
//	overhead = (unsigned long)(size * _maxmemoverhead_);
//	do
//	{
//		exponent=(int)ceil(log_2(size+24));
//		if ( ((size<=overhead)&&((unsigned long)(pow(2.0,(double)exponent)<_maxlinmem_)))
//				|| (pow(2.0,(double)exponent-1.0)<(colcount+24)))
//		{
//			if (m[(int)temprow]) free(m[(int)temprow]);
//			m[(int)temprow]=NULL;
//			size=-1.0;
//		}
//		else
//		{
//			tempsize = pow(2.0,(double)(exponent-1.0))-24.0;
//			if ((unsigned long)tempsize>_maxlinmem_) tempsize=(double)_maxlinmem_;
//			j = (unsigned long)floor(tempsize/colcount);
//			size = (double)(nrow-j-temprow+nrl-1)*sizeof(double)*ncol;
//			if (m[(int)temprow]) free(m[(int)temprow]);
//			m[(int)temprow]=NULL;
//			temprow+=j;
//		}
//	}
//	while (size > 0.0);
//	if (m) free(m); m=NULL;
//}


double **memory_handling::matrix(long nrl, long nrh, long ncl, long nch, long nrows)
{
	long i,j,nrow,ncol,temprow,colcount,overhead;
	long size,tempsize;
	int exponent;
	double **m;
	int mpi_rank = 0;
#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	nrow    = nrh-nrl+2;    /* because of zero-index there should be one row more. normally, the number of rows to reserve would be like nrh-nrl+1 */
	ncol    = nch+1;    /* since we always reserve at least an entire row, the number of columns to allocate will be the same for all grids, wether in MPI mode or in normal mode  */
	temprow = nrl;        /* the row, up to which memory is allocated */
	colcount= ncol*sizeof(double);                                                     /* numberof bytesfor one row */
	/* rows are allowed only up to a lenght of _maxlinmem_ bytes */
	if ((unsigned long)colcount+24>(unsigned long)_maxlinmem_) 
	{
		(*WSstream) << "FATAL: too long matrix rows"<<"[end-entry]"; 
		throw(1);
	}
	size = nrow*ncol*sizeof(double);                                                   /* how much memory is needed? */
	overhead = (unsigned long)(size * (double)_maxmemoverhead_);                /* whatis the overhead we are willing to waste? (say, 5%) */
	if (nrows > 0)
	{
//cout << "in matrix: allocating for " << nrows + 1 << " rows" << " in rank " << mpi_rank << endl;
		m = (double **)calloc((long)(nrows + 1), sizeof(double*));       /* allocating memory for the pointersto the rows (row indizes) */
//		m = new double*[nrows+1];
//cout << "in matrix: after allocating for " << nrows + 1 << " rows" << " in rank " << mpi_rank << endl;
	}
	else
	{
//cout << "in matrix: allocating for " << nrh + 1 << " rows" << " in rank " << mpi_rank << endl;
		m = (double **)calloc((long)(nrh + 1), sizeof(double*));       /* allocating memory for the pointersto the rows (row indizes) */
//		m = new double*[nrh + 1];
	}

	if (!m) 
	{ 
		(*WSstream) << "FATAL: memory allocation error 1 in matrix()"<<"[end-entry]"; 
		throw(1); 
	}
//cout << "in matrix: setting rows to NULL for rows 0 to " << temprow-1 << " in rank " << mpi_rank << endl;
	for (int i = 0; i < temprow; i++)
	{
		m[(int)i] = NULL;
	}

//	for (int i = nrl; i <= nrh; i++)
//	{
//		m[(int)i] = new double[(size_t)ncol + 1];
//	}
	do  /* subdividing the size into powers of 2 */
	{
		exponent=(int)ceil(log_2(size+24.0));  /* exponent for the largest power of 2 that is smaller than size */
		if ( ((size<=overhead)&&((unsigned long)(pow(2.0,(double)exponent)<_maxlinmem_)))
				|| (pow(2.0,(double)exponent-1.0)<(colcount+24)))
		{
			/* this is that part of the loop, which is done at last (if not so much memory was allocated or if
		the rest of memory to allocate is smaller than overhead) */
			m[(int)temprow]=(double *)calloc((size_t)((nrow-temprow+nrl-1)*ncol),sizeof(double));
			if (!m[(int)temprow]) 
			{ 
				(*WSstream) << "FATAL: memory allocation error 2 in matrix()"<<"[end-entry]"; 
				throw(1); 
			}
			for (int i=temprow+1;i<=nrh;i++) m[(int)i]=m[(int)(i-1)]+(int)ncol;
			size=-1.0;
		}
		else
		{
			tempsize = pow(2.0,(double)(exponent-1))-24.0; /* how much memory can be allocated by this exponent */
			if ((unsigned long)tempsize>_maxlinmem_) tempsize=(double)_maxlinmem_; /* if this is more than _maxlinmem_, cut it to that value */
			j = (unsigned long)floor((double)(tempsize/colcount));                           /* largest integer number of rows, allocatable with tempsize bytes*/
			size = (double)(nrow-j-temprow+nrl-1)*sizeof(double)*ncol;          /* remaining size of memory to allocate */
//cout << "in matrix: setting temprow "<<temprow<<" in rank "<< mpi_rank << endl;
			m[(int)temprow] = (double *)calloc((size_t)(j*ncol), sizeof(double));  /* now allocating the memory for <j> rows */
			if (!m[(int)temprow]) 
			{ 
				(*WSstream) << "FATAL: memory allocation error 2 in matrix()"<<"[end-entry]"; 
				throw(1); 
			}
			for (int i = temprow + 1; i <= temprow + j - 1; i++)
			{
				m[(int)i] = m[(int)(i - 1)] + (int)ncol; /* setting the pointers for row indizes */
			}
			temprow+=j; /* in the next lap: this is the row, for wich is now to allocate memory */
		}
	}
	while (size > 0);

//	cout << "in matrix: setting trailing rows " << nrh + 1 << " to "<<nrows<<" in rank " << mpi_rank << endl;
	for (int i = nrh + 1; i <= nrows; i++)
	{
		m[(int)i] = NULL;
	}
	return m;
}

void memory_handling::free_matrix(double **m, long nrl, long nrh, long ncl, long nch, long nrows)
{
	long j,nrow,ncol,temprow,colcount,overhead;
	double size,tempsize;
	int exponent;

	if (!m) return; // return, if the adress is invalid

	nrow = nrh - nrl + 2;    /* because of zero-index there should be one row more. normally, the number of rows to reserve would be like nrh-nrl+1 */
	ncol = nch + 1;    /* since we always reserve at least an entire row, the number of columns to allocate will be the same for all grids, wether in MPI mode or in normal mode  */
	temprow = nrl;        /* the row, up to which memory is allocated */
	colcount = ncol*sizeof(double);
	if ((unsigned long)colcount+24>(unsigned long)_maxlinmem_) 
	{
		(*WSstream) << "FATAL: too long matrix rows"<<"[end-entry]"; 
		throw(1);
	}

//	for (int i = nrl; i <= nrh; i++)
//	{
//		delete[] m[i];
//	}
//	delete[] m;
	size = nrow*ncol*sizeof(double);
	overhead = (unsigned long)(size * _maxmemoverhead_);
	do
	{
		exponent=(int)ceil(log_2(size+24));
		if ( ((size<=overhead)&&((unsigned long)(pow(2.0,(double)exponent)<_maxlinmem_)))
				|| (pow(2.0,(double)exponent-1.0)<(colcount+24)))
		{
			if (m[(int)temprow]) free(m[(int)temprow]);
			m[(int)temprow]=NULL;
			size=-1.0;
		}
		else
		{
			tempsize = pow(2.0,(double)(exponent-1.0))-24.0;
			if ((unsigned long)tempsize>_maxlinmem_) tempsize=(double)_maxlinmem_;
			j = (unsigned long)floor(tempsize/colcount);
			size = (double)(nrow-j-temprow+nrl-1)*sizeof(double)*ncol;
			if (m[(int)temprow]) free(m[(int)temprow]);
			m[(int)temprow]=NULL;
			temprow+=j;
		}
	}
	while (size > 0.0);

	// free memory for row 0 if MPI reinfiltration was active
	if (m[0] != NULL) free(m[0]);

	if (m) free(m); m=NULL;
	
}

double ***memory_handling::drei_d_matrix(long nrl, long nrh, long ncl, long nch, long levels, long nrows)
{
	long nrow,ncol;
	double ***m;
	double *temp;

	nrow = nrh - nrl + 2;    /* because of zero-index there should be one row more. normally, the number of rows to reserve would be like nrh-nrl+1 */
	ncol = nch + 1;    /* since we always reserve at least an entire row, the number of columns to allocate will be the same for all grids, wether in MPI mode or in normal mode  */
	if (nrows>0)
		m = (double ***)calloc((size_t)(nrows + 1), sizeof(double**));
	else
		m = (double ***)calloc((size_t)(nrh + 1), sizeof(double**));
	if (!m)
	{ 
		(*WSstream) << "FATAL: memory allocation error (1) in 3_d_matrix()"<<"[end-entry]"; 
		throw(1); 
	}
	for (int i = 0; i < nrl; i++)
	{
		m[(int)i] = NULL;
	}
	for (int i = nrl; i <= nrh; i++)
	{
		m[i] = (double **)calloc((size_t)(ncol),sizeof(double*));
		if (!m[i])
		{
			(*WSstream) << "FATAL: memory allocation error (2) in 3_d_matrix()"<<"[end-entry]";
			throw(1);
		}
		// assigning the memory linearly to a single temporary variable (will later be used for deleting the entire memory of a row)
		m[i][0] = (double*)calloc((size_t)(levels + 1)*ncol, sizeof(double));
//		temp = (double*)calloc((size_t)(levels+1)*ncol,sizeof(double));

		if (!m[i][0])
		{
			(*WSstream) << "FATAL: memory allocation error (3) in 3_d_matrix()"<<"[end-entry]";
			throw(1);
		}

		for (int j=ncl;j<=nch;j++)
		{
			// instead of allocating memory for every stack, only the pointers are set now
			m[i][j] = m[i][0] + j*(levels + 1);
		}
//#endif	
	}
	for (int i = nrh+1; i <= nrows; i++)
	{
		m[i] = NULL;
	}
	return m;
}

void memory_handling::free_drei_d_matrix(double ***m, long nrl, long nrh, long ncl, long nch, long levels, long nrows)
{
	if (!m) return;

	for (int i=nrl;i<=nrh;i++)
	{
		if (m[i])
		{
			free(m[i][0]); /* Freigeben des Speichers auf die einzelnen Saeulen in den Gridzeilen */
		}
		free(m[i]); /* Freigeben des Speichers auf die Zeilen */
	}
	// free memory for row 0 if MPI reinfiltration was active
	if (m[0] != NULL) free(m[0]);

	/* Freigeben des Speichers auf das Pointerarray */
	if (m) free(m); m=NULL;
}


double memory_handling::log_2(double x)
{
	int i=0;
	while (x>1) {x/=2.0; i++;}
	return i;
	/*  return log(x)*_1_ln2_;             */
}


char *memory_handling::newstring(int count,int errorcode)
{
	char *string;
	if ((string=(char *)calloc(count,sizeof(char)))==NULL) speichererror(errorcode);
	return string;
}

void modelobject::clearstat(double nodata)
{
	stat_nodata = nodata;
	statanz=statmini=statmaxi=sx=sx2=0.0;
}

void modelobject::LogEntry(int iLevel,std::string sMessage)
{
	if (mpi_rank==0)
	{
		time_t *ltime;
		ltime = new(time_t);
		std::ofstream logfile;
		char logfilename[1024];
		sprintf(logfilename,"modell_%d.log",mpi_rank);
		logfile.open(logfilename,std::ios::app);
		time(ltime);
		logfile << ctime(ltime) << ": " << sMessage.c_str() << endl;
		logfile.close();
		delete ltime;
	}
}

void modelobject::addstat(double x,double weight)
{
	if (x>(stat_nodata+1))
	{
		statanz=statanz+weight;
		sx+=x*weight;
		sx2=sx2+x*x*weight;
		if (statanz==1)
		{
			statmini=x;
			statmaxi=x;
		}
		else
		if (x<statmini) statmini=x;
		else if (x>statmaxi) statmaxi=x;
	}
}

void modelobject::statistics(double *anzahl,double *mini,double *maxi,double *summe,double *mittel,double *stddev)
{
	*anzahl=(double)statanz;
	*mini=(double)statmini;
	*maxi=(double)statmaxi;
	*summe=(double)sx;

	if (*anzahl>0) *mittel=(double)( (double)sx / (double)*anzahl); else *mittel=-9999.0;
	if (*anzahl>2)
	{
		*stddev   =(double)( ((double)sx2-(double)sx*(double)sx / (double)*anzahl)/((double)*anzahl-1));
		/* wenn alle Werte gleich sind, ist *stddev == 0, sqrt(0) gibt einen Floating-Point-Fehler */
		if (*stddev>0) *stddev = (double)sqrt(*stddev);
	}
	else
	*stddev=(double)-9999.0;
}

void modelobject::ExtendFileName(char *directory,char *dateiname)
{
	bool bIsComplete=false;
	std::string temp(dateiname);
	// check if the dateiname contains a drive letter (Windows) or starts with a / or ~ (Linux)
#if defined(WIN32)
	if ((temp.find(":\\")==1) ||  // like c:\... in traditional notation
			(temp.find(":/")==1)  ||  // like c:/... in alternative notation
			(temp.find("\\\\")==0)||  // like \\server\... in traditional notation
			(temp.find("\\\\")==0))   // like //server/... in alternative notation
	bIsComplete=true;
#elif defined(WIN64)
	if ((temp.find(":\\")==1) ||  // like c:\... in traditional notation
			(temp.find(":/")==1)  ||  // like c:/... in alternative notation
			(temp.find("\\\\")==0)||  // like \\server\... in traditional notation
			(temp.find("\\\\")==0))   // like //server/... in alternative notation
	bIsComplete=true;
#else
	if (temp.find_first_of("/~")==0) bIsComplete=true; // starts with a / or a ~ meaning there are absolute paths defined under UNIX/LINUX
#endif
	
	if ((!bIsComplete)&&(directory!=NULL))
	{
		char *completename=new char[_stringsize_+1];
		sprintf(completename,"%s%s",directory,dateiname);
		sprintf(dateiname,"%s",completename);
		delete[] completename;
	}
}

void modelobject::SetExternalGridProperties(double xll,double yll,double cs,double **maskgrid)
{
	xllcorner_ext = xll;
	yllcorner_ext = yll;
	cellsize_ext = cs;
	Mask = maskgrid;
}

/* liest ein Grid von der Disk in einen Speicherbereich ein, dessen Adresse durch ***grid
an den Aufrufer zurueckgegebeb wird */
void modelobject::readgridfromfile(char *dateiname, double ***p_grid, int fillmiss, bool bAllocateGridMemory, 
									bool showmessages, double defaultValue, bool limitAutoFill, bool MPIFullGrid, bool MPIRank0Only, bool MPISyncHalo,bool isTimeSeries)
{
	int i,j,k,l,offset;
	long nodatafound;
	FILE *datei,*log;
	double anz,**zn;
	double **buffer;
	long maxdiff,format;
	char *FormatTest=NULL;
	double **GIFbuffer; // temporary grid holding the values from the GIF before interpolating the values to the model grid
	bool bMPINativeReading=false;
	bool doInterpolationAfterReading = false;

#ifdef _MPI
	MPI_File MPIfile;
	MPI_Status MPIstatus;
	MPI_Offset MPIoffset;
#endif

	signal(SIGINT,SIG_IGN); 
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif

	int maxerr=5;
	int retriesleft=maxerr;
	char *completedname = new char[_stringsize_ + 1];

	if (mpi_rank == 0)
	{
		// =============== step 1: extend filename ================================
		sprintf(completedname, "%s", dateiname);
#ifdef _NETCDF
        char *completednamenetcdf = (char *)calloc(_stringsize_+1,sizeof(char));
        int layerindex = index(completedname, "<");
        if(layerindex != 0) {
            // We have a < in the gridname, means the file is a netcdf.
            // We now cut off the layer index. We'll get it back later,
            // so we can access it in the netcdf routines.
            strcpy(completednamenetcdf, completedname);
            substr(completedname, 0, layerindex-1);
        }
#endif
		ExtendFileName(DefaultInputDirectory, completedname);
		//  delete[] completedname;

		// =============== step 2: open file, in case of an error retry up to 5 times ================================
		datei = NULL;
		while (datei == NULL)
		{
			if ((datei = fopen(completedname, "rb")) == NULL)
			{
				if (retriesleft > 0)
				{
					int waitinterval = 200;
					(*WSstream) << "FATAL: error opening " << completedname << ". Will try again in " << waitinterval << " ms... (retriesleft = " << retriesleft << ")" << "[end-entry]";
#ifdef WIN32
					_sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
					// however, _sleep is a deprecated C++ feature, so it must be replaced somewhen in the future
#elif WIN64
					_sleep(waitinterval);	// with Visual V++ 2005, the Sleep() call is not working. There, a .NET-Class Thread::Sleep() should be used, which would make WaSiM a .Net application. 
					// however, _sleep is a deprecated C++ feature, so it must be replaced somewhen in the future
#else
					clock_t goal;
					goal =  waitinterval + clock();
					while( goal > clock() );
					//		      m_sleep(waitinterval);
#endif
				}
				else
				{
					(*WSstream) << "FATAL: still error opening " << completedname << " - no retry." << "[end-entry]";
					throw(21);
				}
				retriesleft = retriesleft - 1; // retry until retriesleft==0
			}
		}
		retriesleft = maxerr;

		// =============== step 3: analyse format (text, binary single, or binary double ================================
		FormatTest = (char*)calloc(100, sizeof(char));
		fread(FormatTest, sizeof(char), 4, datei);
		FormatTest[4] = '\0';

#ifdef _NETCDF
        // Taking advantage of short circuit evaluation to ensure a legal substr() operation
		std::string fourletters = ((string)FormatTest).substr(0, 4);
		if (strlen(FormatTest) > 2 && ((fourletters.find("HDF") == 1) || (fourletters.find("CDF") == 0)))
		{
            format = 10;
        }
        else
#endif
        if ((strcmp(FormatTest, "ncol") == 0) || (strcmp(FormatTest, "NCOL") == 0))
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
					format = 0; // double
				}
				else
				{
					if (strcmp(FormatTest, "GIF8") == 0)
					{
						// call GIF-Decoder and interpolate GIF to model grid (assuming that at least some parts of the model domain are coveerd by the GIF (e.g. RADAR data)
						format = 3; // 0=float, 1=ascii, 2=double, 3=GIF with xll/yll-coordinates and cell siez from header of list file
					}
					else
					{
						(*WSstream) << "FATAL: could not recognise file format (neither ascii, nor single, double, or GIF) for " << completedname << "[end-entry]";
						throw(22);
					}
				}
			}
			// make sure the reading will start from the beginning
		}

		fseek(datei, 0, 0);
		// rewind(datei); not sure if this is really ANSII?

		// =============== step 4: reading header data (if the input grid is a GIF: read the complete grid) and compare with base data from zone grid ================================
		if (showmessages)
		{
			(*WSstream) << "reading " << completedname << " for " << instanzname << "..." << "[end-entry]";
		}
		if (format == 1)	{
			for (i = 1; i <= 6; i++)
			{
				char temp[100];
				if (fscanf(datei, "%s%lf", temp, &(zz[i])) != 2)
				{
					(*WSstream) << "FATAL: error reading ascii value #" << i << " from " << completedname << "[end-entry]";
					throw(22);
				}
			}
		}
		else if (format == 2)	{ // double
			if ((fread(&zz[1], sizeof(double), 12, datei)) != 12)
			{
				(*WSstream) << "FATAL: error reading double value from the header of " << completedname << "[end-entry]";
				throw(22);
			}
		}
		else if (format == 0) { // float (compatibility mode for old grids)
			float fzz[13];
			for (i = 1; i <= 12; i++)
			{
				if ((fread(&(fzz[i]), sizeof(float), 1, datei)) != 1)
				{
					(*WSstream) << "FATAL: error reading double value #" << i << " from the header of " << completedname << "[end-entry]";
					throw(22);
				}
				zz[i] = (double)fzz[i];
			}
		}
#ifdef _NETCDF
        else if (format == 10)  { // netcdf
            // TODO: Is this correct, or should we read stuff from the netcdf?
            zz[1] = ncols;
            zz[2] = nrows;
            zz[3] = xll;
            zz[4] = yll;
            zz[5] = cellsize;
            zz[6] = nodata;
        }
#endif
		else {
			// GIF is read later
			zz[1] = ncols;
			zz[2] = nrows;
			zz[3] = xll;
			zz[4] = yll;
			zz[5] = cellsize;
			zz[6] = nodata;
		}

		/* Bevor diese Routine aufgerufen werden kann, muessen die Werte
	der Gridbeschreibenden Daten bekannt sein (werden als Parameter an den
	Constructor uebergeben */

		//ARNOLD2006: Round xll, yll to next integer
		if (floor(zz[3]) != zz[3])
		{
			zz[3] = floor(zz[3]);
		}
		if (floor(zz[4]) != zz[4])
		{
			zz[4] = floor(zz[4]);
		}

		if ((ncols != (int)floor(zz[1])) ||
			(nrows != (int)floor(zz[2])) ||
			((int)floor(xll) != (int)floor(zz[3])) ||
			((int)floor(yll) != (int)floor(zz[4])) ||
			(cellsize != zz[5]) ||
			(nodata < zz[6]) // nodata value in inputgrid may be less than, but not greater than nodata in WaSiM
			)
		{
			// no error any more, ASCII, float and double grids will now be interpolated automatically. It is the responsibility of the user to guarantee, that category information are available in the same 
			// than the sub basin or zone grid, so e.g. land use and soil type information will not be subject to interpolation (or if, then only nearest neighbor)
			// (*WSstream) << "FATAL: " << completedname << " has not the specified grid measures (in " << instanzname << ")!\n\n"
			// 	<< "file: "<< gridinfo[1].name << " (how it should be) :\n"
			// 	<< "ncols:\t"	  << ncols << "\n" 
			// 	<< "nrows:\t"     << nrows << "\n" 
			// 	<< "xllcorner:\t" << xll << "\n" 
			// 	<< "yllcorner:\t" << yll << "\n" 
			// 	<< "cellsize:\t"  << cellsize << "\n"
			// 	<< "nodata_value:\t" << nodata << "\n" 
			// 	<< "----------------------------------------------------\n\n"
			// 	<< "file: " << completedname << " (how it is):\n"
			// 	<< "ncols:\t"         << zz[1] << "\n"
			// 	<< "nrows:\t"         << zz[2] << "\n"
			// 	<< "xllcorner:\t"     << zz[3] << "\n"
			// 	<< "yllcorner:\t"     << zz[4] << "\n"
			// 	<< "cellsize:\t"      << zz[5] << "\n"
			// 	<< "nodata_value:\t"  << zz[6] << "\n" << "[end-entry]";
			// 
			// if (fclose(datei) == EOF)
			// {
			// 	(*WSstream) << "FATAL: error closing grid " << completedname << "[end-entry]";
			// 	throw(23);
			// }
			// throw(24);
			doInterpolationAfterReading = true;
		}
#ifdef _NETCDF
        // Get the full filename back, including varname and layer
		if (strlen(completednamenetcdf) > 0)
		{
			strcpy(completedname, completednamenetcdf);
		}
        free(completednamenetcdf);
#endif
	}
#ifdef _MPI
	// all the other ranks doesn't know yet about the file name, the format and the header infos, so we need to broadcast these data now
	if (!MPIRank0Only)
	{
		MPI_Bcast(&(completedname[0]), _stringsize_ + 1, MPI_CHAR, 0,MPI_COMM_WORLD);
		MPI_Bcast(&format, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(zz[1]), 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

	// =============== step 5: reserve memory if requested by parameter bAllocateGridMemory ================================
	if ((format == 1) || (format == 3) || (format == 10) || (MPIFullGrid) || (fillmiss!=0) || (doInterpolationAfterReading==true))
	{
		// if the grid is in ascii format or a GIF or netCDF, or if the entire grid is needed because of filling in missing values is requested or 
		// because the entire grid is even needed in this instance (only rank 0 can use complete grids), the Scatter mode is set to false
		bMPINativeReading = false;
	}
	else
	{
		// for common binary files (singe or double) that doesnt need to be filled and that are no GIF, we can read optimized by using 
		// MPI collective I/O and so we do not need to scatter the read in data
		bMPINativeReading = true;
	}
	
	int sr = 1;
	int er = nrows;
	if (bMPINativeReading || (!bMPINativeReading && mpi_rank>0))
	{
		// only a stripe that is defined for the actual rank needs to be allocated
		if (MPISyncHalo == _MPISyncHaloOn_)
		{
			sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
			er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
		}
		else
		{
			sr = mpistartrow;
			er = mpiendrow;
		}
	}

#ifdef _MPI
	if (!MPIRank0Only)
	{
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

	if (bAllocateGridMemory == true)
	{
		// default sr=1 and er=nrows means that the entire grid is allocated
		*p_grid = matrix(sr, er, 1, ncols, nrows);
		for (int i = sr; i <= er; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				(*p_grid)[i][j] = nodata;
			}
		}
	}
	else
	{
		if (*p_grid==NULL)
		{
			(*WSstream) << "FATAL: error reading grid "<<completedname
			<<" no memory allocated (pointer given to function is NULL)"<<"[end-entry]";
			throw(22);
		}
	}

	// =============== step 6: if in MPI scatter mode, close the file in instance 0 and open it for all instances with MPI collective I/O routines ================================
#ifdef _MPI
	if (bMPINativeReading)
	{
		if (mpi_rank == 0)
		{
			if (fclose(datei) == EOF)
			{
				(*WSstream) << "FATAL: error closing grid " << completedname << "[end-entry]";
				throw(23);
			}
		}
		// open grid for MPI collective I/O
		MPI_File_open(MPI_COMM_WORLD, completedname, MPI_MODE_RDONLY, MPI_INFO_NULL,&MPIfile);
		// set Offset where to start reading from
		if (format == 0)
		{ // float: 
			MPIoffset = (12 + ncols*(sr - 1))*sizeof(float);
		}
		else
		{ // double
			MPIoffset = (12 + ncols*(sr - 1))*sizeof(double);
		}
		MPI_File_seek(MPIfile,MPIoffset, MPI_SEEK_SET);
	}
#endif


	// =============== step 7: finally reading the grid ================================
	float *fbuffer=NULL;
	if (format == 0)
	{
		if (doInterpolationAfterReading)
		{
			fbuffer = new float[(int)zz[1] + 1];
		}
		else
		{
			fbuffer = new float[ncols + 1];
		}
	}

#ifdef _NETCDF
    if(format==10) {
		if (mpi_rank == 0)
		{

			double *stackbuf = NULL; // will be allocated in method netcdfIn->readTimestep() and must be deleted by the caller!
			int nc_rows, nc_cols = 0;
			double nc_xll, nc_yll, nc_cellsize,nc_nodata;
			int nc_levels;
			if (isTimeSeries)
			{
				//int htemp = (int)((60.0 * stunde + hour_offset*60 - timestep) / 60.0); // 0...23 --> if hour = 24, then (60*24 + 0 - 1440)/60 = 0, should usually be 0 for daily time steps.
				int htemp = stunde;
				netcdfIn->readTimestep(true, stackbuf, completedname, NULL, nc_rows, nc_cols, nc_xll, nc_yll, nc_cellsize, nc_nodata, nc_levels,-1,jahr,monat,tag,htemp,minute);
			}
			else
			{
				netcdfIn->readTimestep(true, stackbuf, completedname, NULL, nc_rows, nc_cols, nc_xll, nc_yll, nc_cellsize, nc_nodata, nc_levels);
			}
			//  void readTimestep(bool doublemode, double **dataD, float **dataF, 
			//                          char* filename, int layer, char* varname, int nrows, 
			//                        int ncols, double nodata);

			// Wasim doesn't use column and row 0, so we have to copy the data
			// each process starts reading at the first row of its subdomain
			int count = sr - 1;

			// temporary allocate memory for a real 2D-matrix that can be used for interpolating nc-values to the wasim grid more effectively
			// the copy process may be done faster by moving memory with memmove or even by allocating a one dimensional array of double* as nc_grid* 
			// and assigning the adresses of the first elements of each row - but thats more or less what matrix() does anyway, so we don't bother 
			// with this little extra amount of time compared to the reading process itself and the following interpolation
			double **nc_grid = matrix(1, nc_rows, 1, nc_cols);
			for (j = 1; j <= nc_cols; j++) {
				for (i = 1; i <= nc_rows; i++) {
					nc_grid[i][j] = (stackbuf[count] == nc_nodata) ? nodata : stackbuf[count];
					count++;
				}
				count += nrows - (er - sr) - 1;
			}
			delete[] stackbuf; // now deleting the temporary buffer that was allocated in netcdfIn->readTimestep()

			//            if (mpi_rank) {
			//                for (j=1; j<=ncols; j++) {
			//                    for (i=sr; i<=er; i++) {
			//                        (*p_grid)[i][j] = stackbuf[count];
			//                        count++;
			//                    }
			//                    count += nrows - (er-sr) - 1;
			//                }
			//            } 
			//			else 
			//			{
			//				// this works by accident only: for grids that are not read in native mode, sr and er are equal to 1 and nrows, so the entire grid is read in in instance 0
			//				// that means, we could also read in the grid in instance 0 and scatter it in MPI mode.
			//				// also, this means that thsi special loop is not required for rank 0 but cell values for rank 0 can be assigned as for any other rank above
			//                for (j=1; j<=ncols; j++)
			//                    for (i=sr; i<=er; i++) {
			//                        (*p_grid)[i][j] = stackbuf[count];
			//                        count++;
			//                    }
			//            }
			// delete[] stackbuf;

			// now interpolate the nc-values to the wasim grid
			// interpolate only, if the header data from the 
			if (!(this->xllcorner_ext == 9999 && this->yllcorner_ext == 9999 && this->cellsize_ext == 9999))
			{
				GIFReader *gr = new GIFReader(completedname);
				gr->GIFReader::CutModelGrid(nc_grid, 1, nc_cols, 1, nc_rows, nc_xll - 0.5*nc_cellsize, nc_yll - 0.5*nc_cellsize, nc_cellsize,
					*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
//				gr->GIFReader::CutModelGrid(nc_grid, 1, nc_cols, 1, nc_rows, nc_xll, nc_yll, nc_cellsize,
//					*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
				delete gr;
			}
			else
			{
				// here we could check for the ccordinates in the netcdf as well. WaSiM already has two internal grids xcoord and ycoord, holdiung the coordinates of each cell.
				// however, for now we assume that a netCDF file containing the same dimensioned grids in x and y direction also has the same coordinates, so we simply map the values 1:1 here
				double **zn = gridinfo[1].pointer;
				if (nc_rows==nrows && nc_cols==ncols)
				{
					for (int i = 1; i <= nrows; i++)
					{
						for (j = 1; j <= ncols; j++)
						{
							if (modelobject::Mask != NULL)
							{
								if (modelobject::Mask[i][j] != nodata)
									(*p_grid)[i][j] = nc_grid[i][j];
								else
									(*p_grid)[i][j] = nodata;
							}
							else
							{
								if (zn)
								{
									if (zn[i][j] != nodata)
										(*p_grid)[i][j] = nc_grid[i][j];
									else
										(*p_grid)[i][j] = nodata;
								}
								else
								{
									(*p_grid)[i][j] = nc_grid[i][j];
								}
							}
						}
					}
				}
				else
				{
					(*WSstream) << "FATAL: error assigning cell values from a netCDF file to WaSiM grids from file " << completedname
						<< ". number of columns ("<<nc_cols<<") and/or rows ("<<nc_rows<<") in nc-file not matching those in WaSiM ("<<ncols<<" and "<<nrows<<" resp.).\n"
						<<" If the error occurs in interpolation method 9, please define xllcorner, yllcorner and cellsize in the x-, y- and z header rows of the input file for interpolation" << "[end-entry]";
					throw(22);
				}
			}
			free_matrix(nc_grid, 1, nc_rows, 1, nc_cols);
		}
    } else {
#endif
    if (format==0) 
	{ /* single floating point precision values (old binary format) */
		if (doInterpolationAfterReading == false)
		{
			for (i = sr; i <= er; i++)
			{
				if (bMPINativeReading)
				{
#ifdef _MPI
					MPI_Status s;
					MPI_File_read(MPIfile, &(fbuffer[1]), ncols, MPI_FLOAT,&s);
					for (int j = 1; j <= ncols; j++)
					{
						(*p_grid)[i][j] = (double)fbuffer[j];
					}
#endif
				}
				else
				{
					if (mpi_rank == 0)
					{
						if ((int)fread(&(fbuffer[1]), sizeof(float), ncols, datei) != ncols)
						{
							(*WSstream) << "FATAL: error reading double grid " << completedname << "[end-entry]";
							throw(22);
						}
						for (int j = 1; j <= ncols; j++)
						{
							(*p_grid)[i][j] = (double)fbuffer[j];
						}
					}
				}
			}
		}
		else
		{ // read in float grid and do interpolation afterwards
   		  // bMPINativeReading will always be false and se and er are not the ones we expect, so we have to read the entire grid in rank 0 as if it were an ASCII file and interpolate it later to the model grid extend and resolution
			if (mpi_rank == 0)
			{
				double **FullGrid;
				FullGrid = matrix(1, zz[2], 1, zz[1]);
				for (i = 1; i <= (int)zz[2]; i++)
				{
					if ((int)fread(&(fbuffer[1]), sizeof(float), (int)zz[1], datei) != (int)zz[1])
					{
						(*WSstream) << "FATAL: error reading (irregular) float  grid " << completedname << "[end-entry]";
						throw(22);
					}
					for (int j = 1; j <= (int)zz[1]; j++)
					{
						FullGrid[i][j] = (double)fbuffer[j];
					}
				}
				// the grid content is now in FullGrid		
				GIFReader *gr = new GIFReader(completedname);
				gr->GIFReader::CutModelGrid(FullGrid, 1, zz[1], 1, zz[2], zz[3], zz[4], zz[5],
					*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
				free_matrix(FullGrid, 1, zz[2], 1, zz[1]);
				delete gr;
//#ifdef _DEBUG				
//				char tempname[255];
//				sprintf(tempname, "%s_interpolated_float", completedname);
//				printgridfile(tempname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//#endif
			}
		}
	}
	else if (format == 1) 
	{ /* read in in ASCII-Format */
		if (mpi_rank == 0)
		{
//cout << "reading row " << i << " in rank " << mpi_rank << endl;
			if (doInterpolationAfterReading == false)
			{
				for (i = sr; i <= er; i++)
				{
					for (j = 1; j <= ncols; j++)
					{
						if (fscanf(datei, "%lf", &((*p_grid)[i][j])) != 1)
						{
							(*WSstream) << "FATAL: error reading ascii grid " << completedname
								<< " in ASCII-Format in row " << i << ", col " << j << "[end-entry]";
							throw(22);
						}
					}
				}
			}
			else
			{ // read with individual grid geomtry parameters into a temporary grid and interpolate later
				double **FullGrid;
				FullGrid = matrix(1, zz[2], 1, zz[1]);
				for (i = 1; i <= (int)zz[2]; i++)
				{
					for (j = 1; j <= (int)zz[1]; j++)
					{
						if (fscanf(datei, "%lf", &(FullGrid[i][j])) != 1)
						{
							(*WSstream) << "FATAL: error reading ascii grid " << completedname
								<< " in ASCII-Format in row " << i << ", col " << j << "[end-entry]";
							throw(22);
						}
					}
				}
				// the grid content is now in FullGrid		
				GIFReader *gr = new GIFReader(completedname);
				gr->GIFReader::CutModelGrid(FullGrid, 1, zz[1], 1, zz[2], zz[3], zz[4], zz[5],
					*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
				free_matrix(FullGrid, 1, zz[2], 1, zz[1]);
				delete gr;
//#ifdef _DEBUG				
//				char tempname[255];
//				sprintf(tempname, "%s_interpolated_ascii", completedname);
//				printgridfile(tempname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//#endif
			}
		}
	} 
	else if (format==2) 
	{
		if (doInterpolationAfterReading == false)
		{
			for (i = sr; i <= er; i++)
			{
				if (bMPINativeReading)
				{
#ifdef _MPI
					MPI_Status s;
					MPI_File_read(MPIfile,&(*p_grid)[i][1], ncols, MPI_DOUBLE,&s);
#endif
				}
				else
				{
					if (mpi_rank == 0)
					{
						if ((int)fread(&(*p_grid)[i][1], sizeof(double), ncols, datei) != ncols)
						{
							(*WSstream) << "FATAL: error reading double grid " << completedname << "[end-entry]";
							throw(22);
						}
					}
				}
			}
		}
		else
		{ // read with individual grid geomtry parameters into a temporary grid and interpolate later
			if (mpi_rank == 0)
			{
				double **FullGrid;
				FullGrid = matrix(1, zz[2], 1, zz[1]);
				for (i = 1; i <= (int)zz[2]; i++)
				{
					if ((int)fread(&(FullGrid[i][1]), sizeof(double), (int)zz[1], datei) != (int)zz[1])
					{
						(*WSstream) << "FATAL: error reading (irregular) float  grid " << completedname << "[end-entry]";
						throw(22);
					}
				}
				// the grid content is now in FullGrid		
				GIFReader *gr = new GIFReader(completedname);
				gr->GIFReader::CutModelGrid(FullGrid, 1, zz[1], 1, zz[2], zz[3], zz[4], zz[5],
					*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
				free_matrix(FullGrid, 1, zz[2], 1, zz[1]);
				delete gr;
//#ifdef _DEBUG				
//				char tempname[255];
//				sprintf(tempname, "%s_interpolated_double", completedname);
//				printgridfile(tempname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//#endif
			}

		}
	} 
	else 
	{
		if (mpi_rank == 0)
		{
			// read the GIF into a temporary grid which will be allocated within the GIF reading routine and must be deleted after reading and interpolating values into the model grid
			// first setting the grid dimensions to the required output (GIF doesn't contain thise values anyway)
			// setting header values is done when the header infos for other formats ar read
			// zz[1] = ncols;
			// zz[2] = nrows;
			// zz[3] = xll;
			// zz[4] = yll;
			// zz[5] = cellsize;
			// zz[6] = nodata;
			//  call GIF reader
			// Input: 
			// - file handle (already opened file!
			// output: 
			// - pointer to a temporary grid holding the rows and columns of the GIF (Values ranging from 0 to 255)
			// - column count of GIF
			// - row count of GIF
			// note: xllcorner and yllcorner of GIF content as well as cellsize are already read from the list file itself. The GIF does not contain such information
			GIFReader *gr = new GIFReader(completedname);

			int imgwidth;
			int imgheight;
			double **FullGrid;
			gr->GetDimensions(imgwidth, imgheight);
			FullGrid = matrix(1, imgheight, 1, imgwidth);
			try
			{
				gr->ReadGIF(FullGrid);
			}
			catch (std::string ex)
			{
				(*WSstream) << "FATAL: error reading GIF " << completedname << " (" << ex << ")" << "[end-entry]";
				throw(22);
			}
			// the GIF content is now in FullGrid		
			// if the gif was n ot read in as a precip grid, the values for xllcorner, yllcorner and cellsize are assumed to be Swiss coordinates
			if (xllcorner_ext == 0) xllcorner_ext = 255000;
			if (yllcorner_ext == 0) yllcorner_ext = -160000;
			if (cellsize_ext == 0) cellsize_ext = 1000;
			gr->GIFReader::CutModelGrid(FullGrid, 1, imgwidth, 1, imgheight, xllcorner_ext, yllcorner_ext, cellsize_ext,
				*p_grid, 1, ncols, 1, nrows, xll, yll, cellsize, modelobject::Mask);
			// mapping only if the parameters were set
			gr->MapValues(*p_grid, 1, ncols, 1, nrows, this->gmf1, this->gmx1, this->gmcut, this->gmf2, this->gmx2);
			free_matrix(FullGrid, 1, imgheight, 1, imgwidth);
			delete gr;
//#ifdef _DEBUG				
//			char tempname[255];
//			sprintf(tempname, "temp_interpolated_nc");
//			printgridfile(tempname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//#endif
		}
	}

	if ((format == 0) && (fbuffer))	delete[] fbuffer;

#ifdef _NETCDF
    }
#endif


	if (bMPINativeReading)
	{
#ifdef _MPI
		MPI_File_close(&MPIfile);
#endif
	}
	else
	{
		if (mpi_rank == 0)
		{
			if (fclose(datei) == EOF)
			{
				(*WSstream) << "FATAL: error closing grid " << completedname << "[end-entry]";
				throw(23);
			}
		}
	}

//>>> TEST-OUTPUT
/*char fname[512];
FILE *myFile;
sprintf(fname, "data.bin.%d\0", mpi_rank);
std::cout << "TEST OUT FILE: " << fname << std::endl;
std::cout << "nr, nc = " << nrows << ", " << ncols << std::endl;

myFile = fopen(fname, "wb");
for (int n=1; n<=nrows; n++) {
    std::cout << "n = " << n << std::endl;
    fwrite(&((*p_grid)[n][1]), sizeof(double), ncols, myFile);
}
fclose(myFile);

sleep(2);
MPI_Abort(MPI_COMM_WORLD,0);
exit(0);*/
//<<<


	if (FormatTest) free(FormatTest); FormatTest = NULL;

//cout << " before filling missing values " << completedname << " in rank " << mpi_rank << endl;

	// =============== step 8: fill missing values according to parameter fillmiss ================================
	if (((fillmiss == 1) || (fillmiss == 2)) && (mpi_rank==0))
	{
		if (showmessages)
		{
			(*WSstream) << "... " << completedname << " checking nodata values and replace" << "[end-entry]";
		}
		char logfilename[1024];
		sprintf(logfilename, "modell_%d.log", mpi_rank);
		if ((log = fopen(logfilename, "a")) == NULL)
		{
			(*WSstream) << "FATAL: error opening log file " << logfilename << "[end-entry]";
			throw(5);
		}
		nodatafound = 0;
		zn = gridinfo[1].pointer;
		maxdiff = 10; // 10 cells in each direction is the maximum search radius. This should become a dynamic parameter with a value of INF or -1 indicating no limits to the search radius
		buffer = matrix(1, nrows, 1, ncols);
		/* Belegen des Buffers mit Gridwerten */
		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
			{
				// synchronizing nodata with WaSiMs nodata pf -9999
				if ((*p_grid)[i][j] == zz[6])
				{
					(*p_grid)[i][j] = nodata;
				}
				buffer[i][j] = (*p_grid)[i][j];
			}
		}

		for (i = 1; i <= nrows; i++)
		{
			for (j = 1; j <= ncols; j++)
			{
				/* es interessieren die zwei Faelle: 1) ist grid kein Fehlwert, aber zone
												     2) ist grid Fehlwert, aber zone nicht */
				if ((((*p_grid)[i][j]) != zz[6]) && (zn[i][j] == nodata)) buffer[i][j] = nodata;
				if ((((*p_grid)[i][j]) == zz[6]) && (zn[i][j] != nodata)) /* Fehlwert gefunden */
				{
					nodatafound++;
					/* im Datengrid ist ein Fehlwert an einer vom Zonengrid
						bedeckten Stelle aufgetreten: diesen durch das Mittel der um
						gebenden Daten ersetzen */
					if (fillmiss == 1)
					{
						anz = 0;
						offset = 1;
						while ((anz == 0) && (offset <= maxdiff))
						{
							for (k = i - offset; k <= i + offset; k++)
							{
								for (l = j - offset + (i - k); l <= j + offset - (i - k); l++)
								{
									if ((k > 0) && (k <= nrows) && (l > 0) && (l <= ncols))
									{
										if ((*p_grid)[k][l] != nodata)
										{
											anz++;
											/* ersetzungen werden erst in den Buffer geschrieben, um ein gegenseitiges Beeinflussen
							von Fehlwerten zu vermeiden buffer wird dann wieder in Grid kopiert */
											buffer[i][j] = (*p_grid)[k][l];
											break;
										}
									}
								}
								if (anz > 0) break;
							}
							offset++;
						}
						if ((anz == 0) && (limitAutoFill))
						{
							(*WSstream) << "FATAL: In " << completedname << " is a too big hole at x/y = " << i << "/" << j << "!" << "[end-entry]";
							if (mpi_rank == 0) fprintf(log, "In %s is a too big hole at x/y = %d/%d\n!\n", completedname, j, i);
							throw(25);
						}
						if (!limitAutoFill)
						{
							buffer[i][j] = defaultValue;
						}
						else
						{
							(*WSstream) << "DEBUG: " << "y: " << i << " x: " << j << "in " << completedname << ": " << (*p_grid)[i][j] << " replaced by " << buffer[i][j] << "[end-entry]";
							if (mpi_rank == 0) fprintf(log, "y: %d x: %d  in %s %lf replaced by %lf!\n", i, j, completedname, (*p_grid)[i][j], buffer[i][j]);
						}
					}
					else if (fillmiss == 2) {
						buffer[i][j] = defaultValue;
					}
				}
			}

		}
		if (nodatafound > 0)
		{
			(*WSstream) << "DEBUG: In " << completedname << ": " << nodatafound << " nodata values replaced!" << "[end-entry]";
			if (mpi_rank == 0) fprintf(log, "in %s: %ld nodata values replaced!\n\n", completedname, nodatafound);
			for (i = 1; i <= nrows; i++)
			{
				for (j = 1; j <= ncols; j++)
				{
					(*p_grid)[i][j] = buffer[i][j];
				}
			}
		}
		if (mpi_rank == 0) fclose(log);
		free_matrix(buffer, 1, nrows, 1, ncols);
		buffer = NULL;
	}


#ifdef _MPI
	// =============== step 9: if only rank 0 read the grid, it must be scattered now to the other instances ================================
	if ((!bMPINativeReading) && (!MPIRank0Only))
	{
//char *aktname = newstring(_stringsize_, 1);
//sprintf(aktname, "zonegrid_%d_aa", mpi_rank);
//printgridfile(aktname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//cout << "before entering scatter " << completedname << " in rank " << mpi_rank << endl;
#ifdef _MPI
		MPIScatterGrid(0, *p_grid);
#endif
//cout << "after leaving scatter " << completedname<<" in rank " << mpi_rank << endl;
//sprintf(aktname, "zonegrid_%d_ab", mpi_rank);
//printgridfile(aktname, (*p_grid), WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
//delete[] aktname;
	}
	
	// =============== step 10: if bMPINativeReading was false and the entire grid was read in but MPIFullGrid was false, we can now release the allocated grid memory and re-allocate for rank 0 only the required amount of memory ================================
	if ((mpi_rank == 0) && (bMPINativeReading == false) && (MPIFullGrid==false))
	{
//cout << "re-allocating memory for the required fraction of " << completedname << " in rank " << mpi_rank << endl;
		// re-allocate memory for the actual ranks required grid extend only
		// but only for grids that are not used in their entire extend in rank 0
		double **newgrid;
		if (MPISyncHalo == _MPISyncHaloOn_)
		{
			sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow;
			er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow;
		}
		else
		{
			sr = mpistartrow;
			er = mpiendrow;
		}
		newgrid = matrix(sr,er, 1, ncols,nrows);
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = 1; j <= ncols; j++)
			{
				newgrid[i][j] = (*p_grid)[i][j];
			}
		}
		double **tempgrid = (*p_grid);
		(*p_grid) = newgrid;
		free_matrix(tempgrid, 1, nrows, 1, ncols);
	}

	// =============== step 11: snychronize halos (may be we will do this later only if a new parameter requests halos - saves some memory and time, but it's pure optimization) ================================
	if ((!MPIRank0Only) && (MPISyncHalo == _MPISyncHaloOn_))
	{
//cout << "start SyncHalo for " << completedname << " in rank " << mpi_rank << endl;
		MPISynchronizeGridHalo(0, *p_grid);
//cout << "finished SyncHalo for " << completedname << " in rank " << mpi_rank << endl;
	}
#endif
	// cleaning up 
	delete[] completedname;
//	SetSignalDefaults();
	SetSignalAbfangfunktion();
	
}


void modelobject::MPIScatterGrid(int idx,double **grid)
{
#ifdef _MPI

//	cout << "rank: " << mpi_rank << ", MPIREG: " << MPIREG << ", MPIHLP" << MPIHLP << endl;

	// real scattering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1)&&(MPIREG)&&(MPIHLP))
	{
//cout << "MPIHLP " << MPIHLP << "\nMPIREG " << MPIREG << endl;
		if (grid == NULL) grid = gridinfo[idx].pointer;

		if (mpi_rank==0)
		{
			for (int i=1;i<=nrows;i++)
			{
//cout << "FirstLastCol["<<i<<"].First=" << FirstLastCol[i].First << " in rank " << mpi_rank << endl;
				int n = 1;
				if (FirstLastCol[i].First>0)
				{
					int pos=FirstLastCol[i].First;
					int cnt=MPIHLP[i][2*n];
					int rank=MPIHLP[i][2*n-1];
					while (cnt>0)
					{
						if (rank>0)
						{
//cout << "scatter: sending " << cnt << " values in row " << i << " to rank " << rank << endl;
							MPI_Send(&(grid[i][pos]),cnt,MPI_DOUBLE,rank,i,MPI_COMM_WORLD);
						}
						pos+=cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2*n<=ncols)
						{
							cnt=MPIHLP[i][2*n];
							rank=MPIHLP[i][2*n-1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
//cout << "FirstLastCol["<<i<<"].First=" << FirstLastCol[i].First << " in rank " << mpi_rank << endl;
				int n = 1;
				int pos=MPIHLP[i][2*n-1];
				int cnt=MPIHLP[i][2*n];
				while (pos>0)
				{
//cout << "scatter: receiving " << cnt << " values in row " << i << " from rank 0 in rank " << mpi_rank << endl;
					MPI_Status s;
					MPI_Recv(&(grid[i][pos]), cnt, MPI_DOUBLE, 0, i,MPI_COMM_WORLD, &s);
					n++;
					if (2*n<=ncols)
					{
						pos=MPIHLP[i][2*n-1];
						cnt=MPIHLP[i][2*n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
						throw 120;
					}
				}
			}
		}

	}
//	else
//	{
//		cout << "either MPIHLP(" << MPIHLP << ") or MPIREG(" << MPIREG<< ") are NULL" << endl;
//	}
#endif
}



void modelobject::MPIGatherGrid(int idx,double **grid)
{
#ifdef _MPI

	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1)&&(MPIREG)&&(MPIHLP))
	{
		if (grid==NULL) grid=gridinfo[idx].pointer;

		if (mpi_rank==0)
		{
			for (int i=1;i<=nrows;i++)
			{
				int n=1;
				if (FirstLastCol[i].First>0)
				{
					int pos=FirstLastCol[i].First;
					int cnt=MPIHLP[i][2*n];
					int rank=MPIHLP[i][2*n-1];
					while (cnt>0)
					{
						if (rank>0)
						{
							MPI_Status s;
							MPI_Recv(&(grid[i][pos]), cnt, MPI_DOUBLE, rank, i,MPI_COMM_WORLD,&s);
						}
						pos+=cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2*n<=ncols)
						{
							cnt=MPIHLP[i][2*n];
							rank=MPIHLP[i][2*n-1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				int n=1;
				int pos=MPIHLP[i][2*n-1];
				int cnt=MPIHLP[i][2*n];
				while (pos>0)
				{
					MPI_Send(&(grid[i][pos]), cnt, MPI_DOUBLE, 0, i,MPI_COMM_WORLD);
					n++;
					if (2*n<=ncols)
					{
						pos=MPIHLP[i][2*n-1];
						cnt=MPIHLP[i][2*n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
						throw 120;
					}
				}
			}
		}

	}
#endif
}

double **modelobject::MPIGatherIntoNewGrid(int idx, double **grid)
{
	double **newgrid=NULL;
#ifdef _MPI
	
	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1) && (MPIREG) && (MPIHLP))
	{
		if (grid == NULL) grid = gridinfo[idx].pointer;
		if (mpi_rank == 0)
		{
			newgrid = matrix(1, nrows, 1, ncols);
			for (int i = 1; i <= nrows; i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					newgrid[i][j] = (MPIREG[i][j]==mpi_rank)?grid[i][j]:nodata;
				}
				int n = 1;
				if (FirstLastCol[i].First>0)
				{
					int pos = FirstLastCol[i].First;
					int cnt = MPIHLP[i][2 * n];
					int rank = MPIHLP[i][2 * n - 1];
					while (cnt>0)
					{
						if (rank>0)
						{
							MPI_Status s;
							MPI_Recv(&(newgrid[i][pos]), cnt, MPI_DOUBLE, rank, i, MPI_COMM_WORLD,&s);

						}
						pos += cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2 * n <= ncols)
						{
							cnt = MPIHLP[i][2 * n];
							rank = MPIHLP[i][2 * n - 1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row " << i << "[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				int n = 1;
				int pos = MPIHLP[i][2 * n - 1];
				int cnt = MPIHLP[i][2 * n];
				while (pos>0)
				{
					MPI_Send(&(grid[i][pos]), cnt, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
					n++;
					if (2 * n <= ncols)
					{
						pos = MPIHLP[i][2 * n - 1];
						cnt = MPIHLP[i][2 * n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row " << i << "[end-entry]";
						throw 120;
					}
				}
			}
		}

	}
#endif
	return newgrid;
}


void modelobject::MPIGatherStack(int idx,double ***stack)
{
#ifdef _MPI

	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1)&&(MPIREG)&&(MPIHLP))
	{
		if (stack==NULL) stack=gridinfo[idx].stackpointer;

		int levels=gridinfo[idx].levels;
		(*WSstream) << "synchronising stack " << gridinfo[idx].identification << " for mpi-rank " << mpi_rank <<"[end-entry]";
		if (mpi_rank==0)
		{
			for (int i=1;i<=nrows;i++)
			{
				int n=1;
				if (FirstLastCol[i].First>0)
				{
					int pos=FirstLastCol[i].First;
					int cnt=MPIHLP[i][2*n];
//int cnt=FirstLastCol[i].Last-pos+1;
					int rank=MPIHLP[i][2*n-1];
					while (cnt>0)
					{
						if (rank>0)
						{
							//							if ((i==245)||(i==246))
							//							(*WSstream) << "row " << i << ", col " << pos << ", cnt: " << cnt << "(" << cnt*(levels+1) << ")" << "[end-entry]";
//		(*WSstream) << "receiving stack from mpi-rank " << mpi_rank <<" tag: "<<((i-1)*ncols+pos)<<"[end-entry]";
							MPI_Status s;
							MPI_Recv(&(stack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, rank, i, MPI_COMM_WORLD,&s);
						}
						pos+=cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2*n<=ncols)
						{
							cnt=MPIHLP[i][2*n];
							rank=MPIHLP[i][2*n-1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row (in row "<<i<<")"<<"[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				int n=1;
				int pos=MPIHLP[i][2*n-1];
				int cnt=MPIHLP[i][2*n];
				while (pos>0)
				{
					MPI_Send(&(stack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
					
					n++;
					if (2*n<=ncols)
					{
						pos=MPIHLP[i][2*n-1];
						cnt=MPIHLP[i][2*n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row (in row "<<i<<")"<<"[end-entry]";
						throw 120;
					}
				}
			}
		}
	/*	*/
	}
#endif
}

double*** modelobject::MPIGatherIntoNewStack(int idx, double ***stack)
{
	double ***newstack=NULL;
#ifdef _MPI
	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1) && (MPIREG) && (MPIHLP))
	{
		if (stack == NULL) stack = gridinfo[idx].stackpointer;

		int levels = gridinfo[idx].levels;
		(*WSstream) << "synchronising stack " << gridinfo[idx].identification << " for mpi-rank " << mpi_rank << "[end-entry]";
		if (mpi_rank == 0)
		{
			newstack = drei_d_matrix(1, nrows, 1, ncols, levels, nrows);
			for (int i = 1; i <= nrows; i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if (MPIREG[i][j] == mpi_rank)
					{
						for (int k = 0; k <= levels; k++)
						{
							newstack[i][j][k] = stack[i][j][k];
						}
					}
					else
					{
						for (int k = 0; k <= levels; k++)
						{
							newstack[i][j][k] = nodata;
						}
					}

				}

				int n = 1;
				if (FirstLastCol[i].First>0)
				{
					int pos = FirstLastCol[i].First;
					int cnt = MPIHLP[i][2 * n];
					//int cnt=FirstLastCol[i].Last-pos+1;
					int rank = MPIHLP[i][2 * n - 1];
					while (cnt>0)
					{
						if (rank>0)
						{
							MPI_Status s;
							MPI_Recv(&(newstack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, rank, i, MPI_COMM_WORLD, &s);
						}
						pos += cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2 * n <= ncols)
						{
							cnt = MPIHLP[i][2 * n];
							rank = MPIHLP[i][2 * n - 1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row (in row " << i << ")" << "[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				int n = 1;
				int pos = MPIHLP[i][2 * n - 1];
				int cnt = MPIHLP[i][2 * n];
				while (pos>0)
				{
					MPI_Send(&(stack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, 0, i, MPI_COMM_WORLD);

					n++;
					if (2 * n <= ncols)
					{
						pos = MPIHLP[i][2 * n - 1];
						cnt = MPIHLP[i][2 * n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row (in row " << i << ")" << "[end-entry]";
						throw 120;
					}
				}
			}
		}
		/*	*/
	}
#endif
	return newstack;
}


void modelobject::MPIScatterStack(int idx,double ***stack,int levels)
{
#ifdef _MPI

	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)

	if ((mpi_size>1)&&(MPIREG)&&(MPIHLP))
	{
		if (stack == NULL)
		{
			stack = gridinfo[idx].stackpointer;
		}
		if (idx > 0)
		{
			levels = gridinfo[idx].levels;
		}
		else if (levels == 0)
		{
			(*WSstream) << "FATAL: no level-count given as (optional) parameter to MPIScatterStack. Please report to the author."<< "[end-entry]";
			throw 1010;
		}
		if (idx > 0)
		{
			(*WSstream) << "re-distributing stack " << gridinfo[idx].identification << " for mpi-rank " << mpi_rank << "[end-entry]";
		}
		if (mpi_rank==0)
		{
			for (int i=1;i<=nrows;i++)
			{
				int n=1;
				if (FirstLastCol[i].First>0)
				{
					int pos=FirstLastCol[i].First;
					int cnt=MPIHLP[i][2*n];
					int rank=MPIHLP[i][2*n-1];
					while (cnt>0)
					{
						if (rank>0)
						{
							MPI_Send(&(stack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, rank, i, MPI_COMM_WORLD);
						}
						pos+=cnt;
						n++;
						// get the count and the rank of the next strip - if there is nodata, then rank will be nodata and no receive call is made
						if (2*n<=ncols)
						{
							cnt=MPIHLP[i][2*n];
							rank=MPIHLP[i][2*n-1];
						}
						else
						{
							(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
							throw 120;
						}
					}
				}
			}
		}
		else
		{
			for (int i=mpistartrow;i<=mpiendrow;i++)
			{
				int n=1;
				int pos=MPIHLP[i][2*n-1];
				int cnt=MPIHLP[i][2*n];
				while (pos>0)
				{
					MPI_Status s;
					MPI_Recv(&(stack[i][pos][0]), cnt*(levels + 1), MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &s);
					n++;
					if (2*n<=ncols)
					{
						pos=MPIHLP[i][2*n-1];
						cnt=MPIHLP[i][2*n];
					}
					else
					{
						(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
						throw 120;
					}
				}
			}
		}
	}

#endif
}



//LAURA *************** created JustRecv  JustSend and JustWait
void modelobject::MPISynchronizeGridHaloJustRecv(int idx, double **grid)
{
#ifdef _MPI
	if (idx<0) return;
	if (grid == NULL) grid = gridinfo[idx].pointer;
	request_cnt_LAURA = 0;
	for (std::vector<mpihalo>::iterator iter = recvHalo->begin(); iter<recvHalo->end(); iter++)
	{
		int i = iter->row;
		int j = iter->col;
		int cnt = iter->cnt;
		int remote_rank = iter->rank;
		//cout << "LAURA : halo receive " << cnt << "cells from col " << j << " in row " << i << " from rank " << remote_rank <<" in rank "<<mpi_rank<< "  idx=" << idx << endl;
		MPI_Irecv(&(grid[i][j]), cnt, MPI_DOUBLE, remote_rank, i,MPI_COMM_WORLD, &(requestVectorLAURA[request_cnt_LAURA]));
		request_cnt_LAURA++;
	}

#endif

}


void modelobject::MPISynchronizeGridHaloJustSend(int idx, double **grid)
{
#ifdef _MPI
	if (idx<0) return;
	if (grid == NULL) grid = gridinfo[idx].pointer;

	for (std::vector<mpihalo>::iterator iter = sendHalo->begin(); iter<sendHalo->end(); iter++)
	{
		int i = iter->row;
		int j = iter->col;
		int cnt = iter->cnt;
		int remote_rank = iter->rank;
		//cout << "LAURA : halo send " << cnt << "cells from col " << j << " in row " << i << " to rank " << remote_rank <<" in rank "<<mpi_rank<< "  idx=" << idx << endl;
		MPI_Isend(&(grid[i][j]), cnt, MPI_DOUBLE, remote_rank, i,MPI_COMM_WORLD, &(requestVectorLAURA[request_cnt_LAURA]));
		request_cnt_LAURA++;
	}


#endif

}

void modelobject::MPISynchronizeGridHaloJustWait(int idx, double **grid)
{
#ifdef _MPI
	if (idx<0) return;

	MPI_Waitall(request_cnt_LAURA, requestVectorLAURA,statusVector);
	request_cnt_LAURA = 0;

#endif
}



void modelobject::MPISynchronizeGridHalo(int idx,double **grid)
{
#ifdef _MPI

	if (idx<0) return;
	if (grid == NULL) grid = gridinfo[idx].pointer;
	MPI_Request *requestVector;
	requestVector = (MPI_Request*)calloc(sendHalo->size() + recvHalo->size(), sizeof(MPI_Request));
	int request_cnt=0;

	for (std::vector<mpihalo>::iterator iter=recvHalo->begin();iter<recvHalo->end();iter++)
	{
		int i=iter->row;
		int j=iter->col;
		int cnt=iter->cnt;
		int remote_rank=iter->rank;
		MPI_Irecv(&(grid[i][j]), cnt, MPI_DOUBLE, remote_rank, i, MPI_COMM_WORLD,&(requestVector[request_cnt]));
		request_cnt++;
	}

	for (std::vector<mpihalo>::iterator iter = sendHalo->begin(); iter<sendHalo->end(); iter++)
	{
		int i = iter->row;
		int j = iter->col;
		int cnt = iter->cnt;
		int remote_rank = iter->rank;
		MPI_Isend(&(grid[i][j]), cnt, MPI_DOUBLE, remote_rank, i,MPI_COMM_WORLD,&(requestVector[request_cnt]));
		request_cnt++;
	}

	MPI_Waitall(request_cnt, requestVector, statusVector);
	free(requestVector);

#endif
}


void modelobject::MPISynchronizeStackHalo(int idx,int levels,double ***grid)
{
#ifdef _MPI


	if (grid==NULL) grid=gridinfo[idx].stackpointer;
	levels=gridinfo[idx].levels;

	MPI_Request *requestVector;
	requestVector = (MPI_Request*)calloc(sendHalo->size() + recvHalo->size(), sizeof(MPI_Request));
	int request_cnt=0;

	for (std::vector<mpihalo>::iterator iter=sendHalo->begin();iter<sendHalo->end();iter++)
	{
		int i=iter->row;
		int j=iter->col;
		int cnt=iter->cnt;
		int remote_rank=iter->rank;
		MPI_Isend(&(grid[i][j][0]), cnt*(levels + 1), MPI_DOUBLE, remote_rank, 0,MPI_COMM_WORLD,&(requestVector[request_cnt]));
		request_cnt++;
	}

	for (std::vector<mpihalo>::iterator iter=recvHalo->begin();iter<recvHalo->end();iter++)
	{
		int i=iter->row;
		int j=iter->col;
		int cnt=iter->cnt;
		int remote_rank=iter->rank;
		MPI_Irecv(&(grid[i][j][0]), cnt*(levels + 1), MPI_DOUBLE, remote_rank, 0,MPI_COMM_WORLD,&(requestVector[request_cnt]));
		request_cnt++;
	}


	MPI_Waitall(request_cnt,requestVector,statusVector);

	free(requestVector);

#endif
}

void modelobject::MPIGatherZone(statrecord *localzone, int zonecnt, bool gatherrelcellsum)
{
#ifdef _MPI

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	MPI_Datatype zonevectortype;
	MPI_Type_contiguous(5, MPI_DOUBLE, &zonevectortype);
	MPI_Type_commit(&zonevectortype);

	void(*zonesumfunction)(statrecord const *in, statrecord *inout, int *len, const MPI_Datatype *datatype) = &zonesum;
	MPI_Op mySumOp;
	MPI_Op_create((MPI_User_function*)zonesumfunction, true, &mySumOp);


	//LAURA opt 8

	statrecord *tempzone;
	tempzone = new statrecord[zonecnt + 1];

	for (int i = 0; i <= zonecnt; i++)
	{
		tempzone[i].RelCellSum = localzone[i].RelCellSum;  //create a copy
	}


	if (mpi_rank == 0) 												//LAURA
	{ 														//LAURA
		MPI_Reduce(MPI_IN_PLACE, localzone, zonecnt + 1, zonevectortype, mySumOp, 0,MPI_COMM_WORLD);
	}
	else {													//LAURA
		MPI_Reduce(localzone, NULL, zonecnt + 1, zonevectortype, mySumOp, 0,MPI_COMM_WORLD);
	} 	                                                                                                	//LAURA  

	//gatherrelcellsum set same for all ranks
	if (!gatherrelcellsum)												//LAURA
	{														//LAURA
		if (mpi_rank == 0) {
			for (int i = 0; i <= zonecnt; i++)
			{
				localzone[i].RelCellSum = tempzone[i].RelCellSum;  //restore ORIG
			}
		}
	}														//LAURA

//	mySumOp.Free();
//	zonevectortype.Free();
	MPI_Op_free(&mySumOp);
	MPI_Type_free(&zonevectortype);
	delete[] tempzone;
#endif
}



void modelobject::MPIScatterZone(statrecord *localzone, int zonecnt, bool scatterrelcellsum)
{
#ifdef _MPI

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	statrecord *tempzone;
	tempzone = new statrecord[zonecnt + 1];

	const int nitems = 5;
	int blocklengths[5] = { 1, 1, 1, 1, 1 };

	MPI_Datatype types[5] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	MPI_Aint offsets[5] = { offsetof(statrecord, n),
							offsetof(statrecord, sx),
							offsetof(statrecord, zonenwert),
							offsetof(statrecord, rel),
							offsetof(statrecord, RelCellSum) };

	MPI_Datatype zonevectortype;
	MPI_Type_create_struct(nitems, blocklengths, offsets, types,&zonevectortype);
	MPI_Type_commit(&zonevectortype);

	for (int i = 0; i <= zonecnt; i++)
	{
		tempzone[i].n = localzone[i].n;
		tempzone[i].sx = localzone[i].sx;
		tempzone[i].RelCellSum = localzone[i].RelCellSum;
	}

	MPI_Bcast(&(tempzone[0]), zonecnt + 1, zonevectortype, 0, MPI_COMM_WORLD);

	for (int i = 0; i <= zonecnt; i++)
	{
		localzone[i].n				= tempzone[i].n;
		localzone[i].sx				= tempzone[i].sx;
		localzone[i].RelCellSum		= tempzone[i].RelCellSum;
	}
	MPI_Type_free(&zonevectortype);
	delete[] tempzone;
#endif
}

void modelobject::MPIGatherGridWithMask(int idx,int maskidx,long *CellsPerRank,double **OptionalTargetGrid)
{
#ifdef _MPI
	// this version of the GatherFloatGrid method uses another grid as mask (e.g. the lake grid)
	//
	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1)&&(MPIREG))
	{
		double **GRID=gridinfo[idx].pointer;
		double **MASK=gridinfo[maskidx].pointer;
		if (mpi_rank==0)
		{
//			(*WSstream) << "DEBUG: Gather with mask grid "<<idx<<" for rank "<<mpi_rank<<"[end-entry]";
			// create an array of buffers, one for each MPI-Region. Buffer size is the count of MPI cells for that region
			double **buffer = new double*[mpi_size];
			long *counter = new long[mpi_size];
			for (int r=1;r<mpi_size;r++)
			{
				buffer[r] = new double[CellsPerRank[r]];
				MPI_Status s;
				MPI_Recv(&(buffer[r][0]), CellsPerRank[r], MPI_DOUBLE, r, r,MPI_COMM_WORLD,&s);
				counter[r]=0;
			}

			for (int i=1;i<=nrows;i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					int r=MPIREG[i][j];
					if (MASK[i][j]>nodata)
					{
						if (r > 0)
						{
							if (OptionalTargetGrid)
							{
								OptionalTargetGrid[i][j] = buffer[r][counter[r]];
							}
							else
							{
								GRID[i][j] = buffer[r][counter[r]];
							}
							counter[r]++;
						}
						else
						{
							if (r == 0)
							{
								if (OptionalTargetGrid)
								{
									OptionalTargetGrid[i][j] = GRID[i][j];
								}
							}
						}
					}
				}
			}

			delete[] counter;
			for (int r=1;r<mpi_size;r++)
			{
				delete[] buffer[r];
			}
			delete[] buffer;
		}
		else
		{
//			(*WSstream) << "DEBUG: Gather with mask grid "<<idx<<" for rank "<<mpi_rank<<"[end-entry]";
//			double *buffer = new double[v_mpiregioncnt[mpi_rank]];
			double *buffer = new double[CellsPerRank[mpi_rank]];
			long counter = 0;

			for (int i=1;i<=nrows;i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					int r=MPIREG[i][j];
					if ((MASK[i][j]>nodata)&&(r==mpi_rank))
					{
						buffer[counter]=GRID[i][j];
						counter++;
					}
				}
			}

			MPI_Send(&(buffer[0]), CellsPerRank[mpi_rank], MPI_DOUBLE, 0, mpi_rank,MPI_COMM_WORLD);
			delete[] buffer;

//			(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
		}

	}
#else // non-MPI versions should simply copy the non-masked values if the target grid is different from the source grid
	if (OptionalTargetGrid)
	{
		double **GRID=gridinfo[idx].pointer;
		double **MASK=gridinfo[maskidx].pointer;
		for (int i=1;i<=nrows;i++)
		{
			for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if (MASK[i][j]>nodata)
				{
					OptionalTargetGrid[i][j]=GRID[i][j];
				}
			}
		}
	}
#endif
}

void modelobject::MPIScatterGridWithMask(int idx,int maskidx,long *CellsPerRank)
{
#ifdef _MPI

//cout << "rank: " << mpi_rank << endl;

	// this version of the GatherFloatGrid method uses another grid as mask (e.g. the lake grid)
	//
	// real gathering must be done for more than one instances only 
	// (otherwise the only running instance will contain all valid values already)
	if ((mpi_size>1)&&(MPIREG))
	{
		double **GRID=gridinfo[idx].pointer;
		double **MASK=gridinfo[maskidx].pointer;
		if (mpi_rank==0)
		{
//			(*WSstream) << "DEBUG: Scatter grid "<<idx<<" with mask for rank "<<mpi_rank<<"[end-entry]";
			// create an array of buffers, one for each MPI-Region. Buffer size is the count of MPI cells for that region
			double **buffer = new double*[mpi_size];
			long *counter = new long[mpi_size];
			for (int r=1;r<mpi_size;r++)
			{
//				buffer[r] = new double[this->v_mpiregioncnt[r]];
				buffer[r] = new double[CellsPerRank[r]];
				counter[r]=0;
			}

			// filling the buffers for each mpi region with the values from the grid where mask is valid
			for (int i=1;i<=nrows;i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					int r=MPIREG[i][j];
					if ((MASK[i][j]>nodata)&&(r>0))
					{
						buffer[r][counter[r]]=GRID[i][j];
						counter[r]++;
					}
				}
			}

			for (int r=1;r<mpi_size;r++)
			{
				MPI_Send(&(buffer[r][0]), CellsPerRank[r], MPI_DOUBLE, r, r,MPI_COMM_WORLD);
			}

			delete[] counter;
			for (int r=1;r<mpi_size;r++)
			{
				delete[] buffer[r];
			}
			delete[] buffer;
		}
		else
		{
//			(*WSstream) << "DEBUG: Scatter grid "<<idx<<" with mask for rank "<<mpi_rank<<"[end-entry]";
			double *buffer = new double[v_mpiregioncnt[mpi_rank]];
			long counter = 0;

			MPI_Status s;
			MPI_Recv(&(buffer[0]), CellsPerRank[mpi_rank], MPI_DOUBLE, 0, mpi_rank,MPI_COMM_WORLD,&s);
			for (int i=1;i<=nrows;i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					int r=MPIREG[i][j];
					if ((MASK[i][j]>nodata)&&(r==mpi_rank))
					{
						GRID[i][j]=buffer[counter];
						counter++;
					}
				}
			}

			delete[] buffer;

//			(*WSstream) << "FATAL: MPI-Regions: to much strips to be stored in a grid row in row "<<i<<"[end-entry]";
		}

	}
#endif
}


void modelobject::printgridfile(char *dateiname, double **grid, bool bTextMode, bool doublemode, bool WriteNetcdfGrids, int year, int month, int day, int hour, int minute, double hour_offset)
{
	FILE *datei;
	int i, j;

	signal(SIGINT, SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif
	//  (*WSstream) << "Ctrl-C disabled (writing grid files..) ");

	char *completedname = new char[_stringsize_ + 1];
	sprintf(completedname, "%s", dateiname);
	ExtendFileName(DefaultOutputDirectory, completedname);
	//  delete[] completedname;

#ifdef _NETCDF

	if (!WriteNetcdfGrids) 
	{
#endif
		if (!bTextMode)
		{
			if ((datei = fopen(completedname, "wb")) == NULL)
			{
				(*WSstream) << "FATAL: error creating " << completedname << " for writing in binary mode." << "[end-entry]";
				throw(26);
			}
			zz[1] = ncols; zz[2] = nrows; zz[3] = xll; zz[4] = yll; zz[5] = cellsize; zz[6] = nodata;
			zz[7] = zz[8] = zz[9] = zz[10] = zz[11] = zz[12] = 0.0;

			if (!doublemode)
			{
				for (int i = 1; i <= 12; i++)
				{
					float v = zz[i];
					if ((fwrite(&v, sizeof(float), 1, datei)) != 1) writeerror(completedname, 18);
				}
			}
			else
			{
				if ((fwrite(&zz[1], sizeof(double), 12, datei)) != 12) writeerror(completedname, 18);
			}

			clearstat(nodata);
			float *row;
			double *drow;
			if (!doublemode)
			{
				row = new float[ncols + 1];
			}
			else
			{
				drow = new double[ncols + 1];
			}
			for (i = 1; i <= nrows; i++)
			{
				if (grid[i]) // when writing MPI grids per instance, statistic can be evaluated for valid rows only
				{
					//    for (j=1;j<FirstLastCol[i].First;j++) grid[i][j]=nodata;
					if (FirstLastCol[i].Last > 0)
					{
						//      for (j=FirstLastCol[i].Last+1;j<=ncols;j++) grid[i][j]=nodata;
						if (RelCellArea != NULL)
							for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) addstat(grid[i][j], RelCellArea[i][j]);
						else
							for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) addstat(grid[i][j], 1.0);
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							addstat(grid[i][j], 1.0);
						}
					}
				}
				if (!doublemode)
				{
					if (grid[i])
					{
						for (j = 1; j <= ncols; j++)
						{
							row[j] = (float)grid[i][j];
						}
						if (((int)fwrite(&(row[1]), sizeof(float), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							row[j] = nodata;
						}
						if (((int)fwrite(&(row[1]), sizeof(float), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
				}
				else
				{
					if (grid[i])
					{
						if (((int)fwrite(&grid[i][1], sizeof(double), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							drow[j] = nodata;
						}
						if (((int)fwrite(&(drow[1]), sizeof(double), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
				}

			}
			if (!doublemode)
			{
				delete[] row;
			}
			else
			{
				delete[] drow;
			}

			statistics(&(zz[7]), &(zz[8]), &(zz[9]), &(zz[10]), &(zz[11]), &(zz[12]));

			if (!doublemode)
			{
				fseek(datei, 6 * sizeof(float), 0);
				for (int i = 7; i <= 12; i++)
				{
					float v = zz[i];
					if ((fwrite(&v, sizeof(float), 1, datei)) != 1) writeerror(completedname, 18);
				}
			}
			else
			{
				fseek(datei, 6 * sizeof(double), 0);
				fwrite(&(zz[7]), sizeof(double), 6, datei);
			}

			if ((fclose(datei)) != 0)
			{
				(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
				throw(27);
			}
		}
		else
		{
			if ((datei = fopen(completedname, "wt")) == NULL)
			{
				(*WSstream) << "FATAL: error creating " << completedname << " for writing in text mode." << "[end-entry]";
				throw(26);
			}
			if (fprintf(datei,
				"ncols\t%d\n"
				"nrows\t%d\n"
				"xllcorner\t%lf\n"
				"yllcorner\t%lf\n"
				"cellsize\t%lf\n"
				"nodata_value\t%lf\n", ncols, nrows, xll, yll, cellsize, nodata) < 0) writeerror(completedname, 18);

			for (i = 1; i <= nrows; i++)
			{
				if (grid[i])
				{
					//if (mpi_rank==0)
					//cout << "writing " << completedname<< "row " << i << "in rank " << mpi_rank << "from " << grid[i] << endl;
					for (j = 1; j <= ncols; j++)
					{
						if (fprintf(datei, "%g", grid[i][j]) < 0)  writeerror(completedname, 18);
						if (j < ncols) fprintf(datei, "\t");
					}
				}
				else
				{
					//if (mpi_rank == 0)
					//cout << "writing " << completedname << " with nodata row " << i << "in rank " << mpi_rank << "from " << grid[i] << endl;
					for (j = 1; j <= ncols; j++)
					{
						if (fprintf(datei, "-9999") < 0)  writeerror(completedname, 18);
						if (j < ncols) fprintf(datei, "\t");
					}
				}
				fprintf(datei, "\n");
			}
			if ((fclose(datei)) != 0)
			{
				(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
				throw(27);
			}
		}

#ifdef _NETCDF
	}
    if(WriteNetcdfGrids) {
        
        netcdfOut->setGeometry(xll, yll, cellsize);
        
        // Determine the file name of the netcdf file and the name of
        // the variable that is supposed to be written. We extract this
        // from the file name of the grid file.
        char *completednamenetcdf = (char *)calloc(_stringsize_+1,sizeof(char));
        netcdfhlp::findFileNameFromGridName(completednamenetcdf, completedname, false);
        
        char *varname = (char *)calloc(_stringsize_+1,sizeof(char));
        netcdfhlp::findVarNameFromGridName(varname, completedname);
        
        double *latbuf, *lonbuf;
        
        if(gridinfo[1].has_latgrid_) {
            latbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
            for (j=1; j<=ncols; j++) {
				for (i = 1; i <= nrows; i++) {
						latbuf[count] = (gridinfo[gridinfo[1].latgrid_index_].pointer[i]) ? gridinfo[gridinfo[1].latgrid_index_].pointer[i][j] : nodata;
                    count++;
                }
            }
            netcdfOut->setLatGrid(latbuf);
        }
        
        if(gridinfo[1].has_longrid_) {
            lonbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					lonbuf[count] = (gridinfo[gridinfo[1].longrid_index_].pointer[i])?gridinfo[gridinfo[1].longrid_index_].pointer[i][j]:nodata;
                    count++;
                }
            }
            netcdfOut->setLonGrid(lonbuf);
        }
        
        // Since Wasim doesn't use row and column 0 of the array, we have to copy the data for netcdf
        if (doublemode) {

            double *gridbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
                    gridbuf[count] = (grid[i])?grid[i][j]:nodata;
                    count++;
                }
            }
			// int htemp = (int)((60.0 * hour + hour_offset * 60 - timestep) / 60.0); // 0...23 --> if hour = 24, then (60*24 + 0 - 1440)/60 = 0, should usually be 0 for daily time steps.
			int htemp = hour;
			netcdfOut->writeTimestep(true, gridbuf, NULL, completednamenetcdf, varname, nrows, ncols, nodata,cellsize,year,month,day,htemp,minute);
            delete[] gridbuf;
            
        } else {
            
            float *gridbuf = new float[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					gridbuf[count] = (grid[i]) ? grid[i][j] : nodata;
                    count++;
                }
            }
			// int htemp = (int)((60.0 * hour + hour_offset*60 - timestep) / 60.0); // 0...23 --> if hour = 24, then (60*24 + 0 - 1440)/60 = 0, should usually be 0 for daily time steps.
			int htemp = hour;
			netcdfOut->writeTimestep(false, NULL, gridbuf, completednamenetcdf, varname, nrows, ncols, nodata,cellsize, year, month, day, htemp, minute);
            delete[] gridbuf;
        }
        
        if(gridinfo[1].has_latgrid_) {
            delete[] latbuf;
        }

        if(gridinfo[1].has_longrid_) {
            delete[] lonbuf;
        }
        
        free(completednamenetcdf);
        free(varname);
    }
    
#else
    if(WriteNetcdfGrids) {
        (*WSstream) << "WARNING: Binary compiled without NetCDF support, skipping output."<<"[end-entry]";
    }
#endif
	
	delete[] completedname;
	SetSignalAbfangfunktion();
}

/*ueberladene Form fuer double **grids  */
void modelobject::printsummengridfile(char *dateiname, double **grid, bool bTextMode, bool bDoubleMode, bool WriteNetcdfGrids)
{
	FILE *datei;
	int i, j;
	double floatwert;
	signal(SIGINT, SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif
	//  (*WSstream) << "Ctrl-C disabled (writing sum grids..) ");
	char *completedname = new char[_stringsize_ + 1];
	sprintf(completedname, "%s", dateiname);
	ExtendFileName(DefaultOutputDirectory, completedname);
	//  delete[] completedname;

#ifdef _NETCDF

	if (!WriteNetcdfGrids) 
	{
#endif
		if (!bTextMode)
		{
			if ((datei = fopen(completedname, "wb")) == NULL)
			{
				(*WSstream) << "FATAL: error creating " << completedname << " for writing in binary mode." << "[end-entry]";
				throw(26);
			}
			zz[1] = ncols; zz[2] = nrows; zz[3] = xll; zz[4] = yll; zz[5] = cellsize; zz[6] = nodata;
			zz[7] = zz[8] = zz[9] = zz[10] = zz[11] = zz[12] = 0.0;

			if (!bDoubleMode)
			{
				for (int i = 1; i <= 12; i++)
				{
					float v = zz[i];
					if ((fwrite(&v, sizeof(float), 1, datei)) != 1) writeerror(completedname, 18);
				}
			}
			else
			{
				if ((fwrite(&zz[1], sizeof(double), 12, datei)) != 12) writeerror(completedname, 18);
			}

			clearstat(nodata);
			float *row;
			double *drow;
			if (!bDoubleMode)
			{
				row = new float[ncols + 1];
			}
			else
			{
				drow = new double[ncols + 1];
			}
			for (int i = 1; i <= nrows; i++)
			{
				if (grid[i])
				{
					//    (*WSstream) << "firstlastcol: %d ... %d\n",FirstLastCol[i].First,FirstLastCol[i].Last);
					for (j = 1; j < FirstLastCol[i].First; j++)
					{
						grid[i][j] = (double)nodata;
					}
					if (FirstLastCol[i].Last>0)
					{
						for (j = FirstLastCol[i].Last + 1; j <= ncols; j++) grid[i][j] = (double)nodata;
						if (RelCellArea != NULL)
							for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) addstat((double)grid[i][j], RelCellArea[i][j]);
						else
							for (j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++) addstat((double)grid[i][j], 1.0);
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							grid[i][j] = (double)nodata;
							addstat(grid[i][j], 1.0);
						}
					}
					if (!bDoubleMode)
					{
						for (j = 1; j <= ncols; j++)
						{
							row[j] = (float)grid[i][j];
						}
						if (((int)fwrite(&(row[1]), sizeof(float), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
					else
					{
						if (((int)fwrite(&grid[i][1], sizeof(double), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
				}
				else
				{
					if (!bDoubleMode)
					{
						for (j = 1; j <= ncols; j++)
						{
							row[j] = nodata;
						}
						if (((int)fwrite(&(row[1]), sizeof(float), ncols, datei)) != ncols) writeerror(completedname, 18);
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							drow[j] = nodata;
						}
						if (((int)fwrite(&(drow[1]), sizeof(double), ncols, datei)) != ncols) writeerror(completedname, 18);
					}

				}

			}
			if (!bDoubleMode)
			{
				delete[] row;
			}
			else
			{
				delete[] drow;
			}
			statistics(&(zz[7]), &(zz[8]), &(zz[9]), &(zz[10]), &(zz[11]), &(zz[12]));

			if (!bDoubleMode)
			{
				fseek(datei, 6 * sizeof(float), 0);
				for (int i = 7; i <= 12; i++)
				{
					float v = zz[i];
					if ((fwrite(&v, sizeof(float), 1, datei)) != 1) writeerror(completedname, 18);
				}
			}
			else
			{
				fseek(datei, 6 * sizeof(double), 0);
				fwrite(&(zz[7]), sizeof(double), 6, datei);
			}

			if ((fclose(datei)) != 0)
			{
				(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
				throw(27);
			}
		}
		else
		{
			if ((datei = fopen(completedname, "wt")) == NULL)
			{
				(*WSstream) << "FATAL: error creating " << completedname << " for writing in text mode." << "[end-entry]";
				throw(26);
			}

			if (fprintf(datei,
				"ncols\t%d\n"
				"nrows\t%d\n"
				"xllcorner\t%lf\n"
				"yllcorner\t%lf\n"
				"cellsize\t%lf\n"
				"nodata_value\t%lf\n", ncols, nrows, xll, yll, cellsize, nodata) < 0) writeerror(completedname, 18);

			for (i = 1; i <= nrows; i++)
			{
				if (grid[i])
				{
					for (j = 1; j < FirstLastCol[i].First; j++) grid[i][j] = (double)nodata;
					if (FirstLastCol[i].Last>0)
					{
						for (j = FirstLastCol[i].Last + 1; j <= ncols; j++) grid[i][j] = (double)nodata;
					}
					else
					{
						for (j = 1; j <= ncols; j++)
						{
							grid[i][j] = (double)nodata;
						}
					}
					for (j = 1; j <= ncols; j++)
					{
						if (fprintf(datei, "%lg", grid[i][j]) < 0)  writeerror(completedname, 18);
						if (j < ncols) fprintf(datei, "\t");
					}
					fprintf(datei, "\n");
				}
				else
				{
					for (j = 1; j <= ncols; j++)
					{
						if (fprintf(datei, "%lg", nodata) < 0)  writeerror(completedname, 18);
						if (j < ncols) fprintf(datei, "\t");
					}
					fprintf(datei, "\n");
				}
			}
			if ((fclose(datei)) != 0)
			{
				(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
				throw(27);
			}

		}

#ifdef _NETCDF
	}
    if(WriteNetcdfGrids) {

        netcdfOut->setGeometry(xll, yll, cellsize);

        char *completednamenetcdf = (char *)calloc(_stringsize_+1,sizeof(char));
        netcdfhlp::findFileNameFromGridName(completednamenetcdf, completedname, true);

        char *varname = (char *)calloc(_stringsize_+1,sizeof(char));
        netcdfhlp::findVarNameFromGridName(varname, completedname);
        
        double *latbuf, *lonbuf;
        
        if(gridinfo[1].has_latgrid_) {
            latbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					latbuf[count] = (gridinfo[gridinfo[1].latgrid_index_].pointer[i])?gridinfo[gridinfo[1].latgrid_index_].pointer[i][j]:nodata;
                    count++;
                }
            }
            netcdfOut->setLatGrid(latbuf);
        }
        
        if(gridinfo[1].has_longrid_) {
            lonbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					lonbuf[count] = (gridinfo[gridinfo[1].longrid_index_].pointer[i])?gridinfo[gridinfo[1].longrid_index_].pointer[i][j]:nodata;
                    count++;
                }
            }
            netcdfOut->setLonGrid(lonbuf);
        }

        // Since Wasim doesn't use row and column 0 of the array, we have to copy the data for netcdf
        if (bDoubleMode) {
            double *gridbuf = new double[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					gridbuf[count] = (grid[i])?grid[i][j]:nodata;
                    count++;
                }
            }
            netcdfOut->writeTimestep(true, gridbuf, NULL, completednamenetcdf, varname, nrows, ncols, nodata,cellsize);
            delete[] gridbuf;
        } else {
            float *gridbuf = new float[(size_t)ncols * nrows];
            size_t count = 0;
			for (j = 1; j <= ncols; j++) {
				for (i = 1; i <= nrows; i++) {
					gridbuf[count] = (grid[i])?grid[i][j]:nodata;
                    count++;
                }
            }
            netcdfOut->writeTimestep(false, NULL, gridbuf, completednamenetcdf, varname, nrows, ncols, nodata,cellsize);
            delete[] gridbuf;
        }

        free(completednamenetcdf);
        free(varname);
    }

#else
    if(WriteNetcdfGrids) {
        (*WSstream) << "WARNING: Binary compiled without NetCDF support, skipping output."<<"[end-entry]";
    }
#endif
	
	delete[] completedname;
	SetSignalAbfangfunktion();
}


void modelobject::writegridtofile(int idx, int gather, bool doNotPrintAverageYet)
/* schreibt ein durch **gridinfo[idx].pointer adressiertes Grid aus dem Speicher auf die Disk
Codes :
1000er Stelle: 1 = ein Summengrid über die Dauer der Vegetationsperiode wird gechrieben
				2 = ein Mittelwertgrid über die Dauer der Vegetationsperiode wird gechrieben
	100er Stelle: 1 = ein Minimumgrid wird geschrieben
				2 = ein Maximumgrid wird geschrieben
				3 = Minumum- und Maximumgrids werden geschrieben
	10er Stelle: 0 = kein Summengrid wird geschrieben
				1 = es wird nur am Modellende ein gesamtes Summengrid geschrieben
				2 = es wird pro Jahr ein Summengrid geschrieben
				3 = es wird pro Monat ein Summengrid geschrieben
				4 = es wird pro Tag ein Summengrid geschrieben (nur, wenn Intervalldauer < 1 Tag ist!)
				5 = es wird am Modellende ein Durchschnitssgrid geschrieben
				6 = es wird jedes Jahr ein Durchschnittsgrid geschrieben
				7 = es wird jeden Monat ein Durchschnitssgrid geschrieben
				8 = es wird fuer jeden Tag ein Durchschnitssgrid geschrieben
	1er Stelle:  (gilt für die aktuellen Werte der Groessen, also keine Summen oder Durchschnittswerte!)
				1 = schreiben nach jedem Intervall in dasselbe Grid
				2 = immer ein neues Grid pro Intervall (Name: "<%.4s>%02d%02d.%02d",name,monat,tag,stunde)
				3 = nur das letzte Grid des Modellaufes wird geschrieben
				4 = immer zu 24 Uhr wird ein Grid geschrieben (immer dergleiche Name, also Ueberschreiben des alten Grids)
				5 = immer zu 24 Uhr wird ein neues Grid geschrieben         
				6 = am Monatsende
				7 = am Jahresende
				8 = dieses Grid wird nach dem Einlesen sofort rausgeschrieben, so dass man ein gefülltest Grid ausgegeben erhält und nicht die Lücken manuell füllen muss.
				*/
{
	if ((idx>gridmaxanzahl)||(idx<1)) return;
	int      writegrid       = gridinfo[idx].writegrid;
	int      summenwritegrid = gridinfo[idx].summenwritegrid;
	int      minwritegrid    = gridinfo[idx].minwritegrid;
	int      maxwritegrid    = gridinfo[idx].maxwritegrid;
	int      vegWritegrid    = gridinfo[idx].vegWritegrid;

	if ((writegrid==0)&&(summenwritegrid==0)&&(maxwritegrid==0)&&(minwritegrid==0)&&(vegWritegrid==0)) return;

	char     *identification = gridinfo[idx].identification;
	char     *name        = gridinfo[idx].outname;
	double    **grid       = gridinfo[idx].pointer;
	double    ***stack     = gridinfo[idx].stackpointer;
	char     *summenname  = gridinfo[idx].summenname;
	double   **summengrid = gridinfo[idx].summenpointer;
	char     *minname     = gridinfo[idx].minname;
	double    **mingrid    = gridinfo[idx].minpointer;
	char     *maxname     = gridinfo[idx].maxname;
	double    **maxgrid    = gridinfo[idx].maxpointer;
	double   *anzahl      = &(gridinfo[idx].anzahl);       /* for averaging the number of summed intervals has to known */
	int      *minmaxflag  = &(gridinfo[idx].minmaxflag);   /* the first time this funtion is called, minmax is still 0 and is then set to 1 (0=initialization with grid values) */
	double    **zn         = gridinfo[1].pointer;
	int      i,j,p;
	char     *dateiname,*tempstr;

	int      *vegWriteFlag    = &(gridinfo[idx].vegWriteFlag);
	int      *vegResetFlag    = &(gridinfo[idx].vegResetFlag);
	double   **vegSumPointer = gridinfo[idx].vegSumPointer;
	char     *vegName        = gridinfo[idx].vegName;
	double    **VEGSTART      = gridinfo[idx].VEGSTART;
	double    **VEGSTOP       = gridinfo[idx].VEGSTOP;
	double    **VEGDURAT      = gridinfo[idx].VEGDURAT;

	if (grid)
	{
		if ((MPIREG==NULL)&&(ist_vorhanden(_MPIRegions_)))
		{
			MPIREG=gridinfo[getindex(_MPIRegions_)].pointer;
		}
		if ((MPIHLP==NULL)&&(ist_vorhanden(_MPIHelper_)))
		{
			MPIHLP=gridinfo[getindex(_MPIHelper_)].pointer;
		}
		/* wenn durchschnitsswerte berechnet werden sollen, muss die Anzahl der Zeitintervalle bekannt sein, fuer welche
	aus dem Summengrid nach Erreichen des jeweils festgelegten Zeitpunktes der Durchschnitsswert berechnet werden soll.
	Die ANzahl der Zeitintervalle wird fuer jedes Grid fuer alle Zellen gemeisam konstant in gridinfo[].anzahl gehalten.
	pro Zeitintervall wird der Wert um 1 erhoeht. Da das Gridinfo-array durch calloc eingerichtet wurde, ist anzahl am Beginn
	immer 0 und braucht nicht initialisiert zu werden */
		gridinfo[idx].anzahl++;

		/* if stack is valid AND there should be a sum-, mean-,min- or maxgrid, the values from the uppermost level are stored into the grid */
		if (stack)
		{
			/* gridinfo[idx].getfirstlevel is always zero, except, if it was set to 1 explicitly, e.g. to prevent the pointer[i][j]-field from overwriting like it is used in the concentration stack */
			if ((gridinfo[idx].getfirstlevel == 0) && ((summenwritegrid > 0) || (minwritegrid > 0) || (maxwritegrid > 0)))
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j] > nodata) && ((MPIREG[i][j] == mpi_rank) || (gather == 0)))
						{
							grid[i][j] = stack[i][j][1];
						}
					}
				}
			}
		}

		/* egal, welchen Code Summenwritegrid enthaelt,
	wenn er > 1 ist wird auf jeden Fall aufsummiert. Die Mittelwertbildung erfolgt erst am Ende des Aufsummierungsintervalles */
		if (summenwritegrid!=0)
		{
			for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
					{
						summengrid[i][j]+=(double)grid[i][j]; 
					}
					else 
					{
						summengrid[i][j] = nodata;
					}
				}
			}
		}

		if (*minmaxflag==0)
		{
			if (minwritegrid==1)
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							mingrid[i][j]=grid[i][j];
						}
						else
						{
							mingrid[i][j]=nodata;
						}
					}
				}
			}
			if (maxwritegrid==1)
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							maxgrid[i][j]=grid[i][j];
						}
						else
						{
							maxgrid[i][j]=nodata;
						}
					}
				}
			}
			*minmaxflag=1;
		}
		else
		{
			if (minwritegrid==1)
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							if (grid[i][j]<mingrid[i][j]) mingrid[i][j]=grid[i][j];
						}
					}
				}
			}
			if (maxwritegrid==1)
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							if (grid[i][j]>maxgrid[i][j]) maxgrid[i][j]=grid[i][j];
						}
					}
				}
			}
		}

		// updating vegSumPointer
		if (vegWritegrid!=0)
		{
			// check, if VEGSTART, VEGSTOP, VEGDURAT are already known; if not: find them!
			if (!VEGSTART) 
			{
				char *tempname;
				int pos=0;
				tempname=new char[_stringsize_+1];
				std::string s(name);
				if (s.find("_Layer_")!=std::string::npos)
				{   // try to find the Layer number at the end of the Name.
					// per default, the VEG* grids with Number 1 will be looked for
					s=s.substr(s.find("_Layer_")+7);
					s=s.substr(0,s.find("."));
					sprintf(tempname,"%s%s",_VegetationStart_,s.c_str());	
				}
				else
				{
					sprintf(tempname,"%s1",_VegetationStart_);	
				}
				int vegstartidx=getindex(tempname);
				VEGSTART=gridinfo[vegstartidx].pointer;
				gridinfo[idx].VEGSTART=VEGSTART;
				delete[] tempname;
			}
			if (!VEGSTOP) 
			{
				char *tempname;
				int pos=0;
				tempname=new char[_stringsize_+1];
				std::string s(name);
				if (s.find("_Layer_")!=std::string::npos)
				{   // try to find the Layer number at the end of the Name.
					// per default, the VEG* grids with Number 1 will be looked for
					s=s.substr(s.find("_Layer_")+7);
					s=s.substr(0,s.find("."));
					sprintf(tempname,"%s%s",_VegetationStop_,s.c_str());	
				}
				else
				{
					sprintf(tempname,"%s1",_VegetationStop_);	
				}
				int vegstopidx=getindex(tempname);
				VEGSTOP=gridinfo[vegstopidx].pointer;
				gridinfo[idx].VEGSTOP=VEGSTOP;
				delete[] tempname;
			}
			if (!VEGDURAT) 
			{
				char *tempname;
				int pos=0;
				tempname=new char[_stringsize_+1];
				std::string s(name);
				if (s.find("_Layer_")!=std::string::npos)
				{   // try to find the Layer number at the end of the Name.
					// per default, the VEG* grids with Number 1 will be looked for
					s=s.substr(s.find("_Layer_")+7);
					s=s.substr(0,s.find("."));
					sprintf(tempname,"%s%s",_VegetationDuration_,s.c_str());	
				}
				else
				{
					sprintf(tempname,"%s1",_VegetationDuration_);	
				}
				int vegduratidx=getindex(tempname);
				VEGDURAT=gridinfo[vegduratidx].pointer;
				gridinfo[idx].VEGDURAT=VEGDURAT;
				delete[] tempname;
			}

			// now Updating depending on content of VEGSTART and VEGSTOP 
			for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
					{
						if ((VEGSTART[i][j]>=0)&&(VEGSTOP[i][j]<=0))
						{
							if (grid[i][j]!=nodata)
							{
								vegSumPointer[i][j]+=grid[i][j];
							}
						}
					}
				}
			}

			if (*vegWriteFlag==1)
			{
				// write the grid now (either at the end of the model period or at the vegWrite-Day defined in the LUtable (usually at both times)
				// ... to be implemented
				*vegWriteFlag = 0;
				if (vegWritegrid == 2) // Mittelwert
				{
					for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
					{
						for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						{
							if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
							{
								vegSumPointer[i][j] = vegSumPointer[i][j]/((double)VEGDURAT[i][j]);
							}
						}
					}
				}

				/* Namen des Jahresgrids zusammensetzen */
				tempstr=new char[_stringsize_];

				/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
				if (strrchr(vegName,dirzeichen)!=NULL)
				strcpy(tempstr,strrchr(vegName,dirzeichen)+1); else strcpy(tempstr,vegName);
				dateiname = new char[_stringsize_];
				if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
				p = index(tempstr,".");
				//        if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;
				if (vegWritegrid == 2)
				sprintf(dateiname,"%.*s_%04ld_vegper.mit",p,tempstr,jahr); /* Summengrids fuers Jahr enden auf "xxxxjj00.sum" */
				else
				sprintf(dateiname,"%.*s_%04ld_vegper.sum",p,tempstr,jahr); /* Mittelwertgrids fuers Jahr enden auf "xxxxjj00.mit" */
				/* Pfad wird wieder angehaengt */
				if (strrchr(vegName,dirzeichen)!=NULL)
				{
					strcpy(tempstr,dateiname);
					p=(int)(strrchr(vegName,dirzeichen)-vegName)+1;
					sprintf(dateiname,"%.*s%s",p,vegName,tempstr);
				}
				/* das im Speicher an der Adresse **summengrid gehaltene Array wird in ein Grid geschrieben */
				double **tempgrid;
				if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, vegSumPointer); else tempgrid = vegSumPointer;
				if (mpi_rank == 0)
				{
					printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
				if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);

				/* nun muessen die im Summengrid gehaltenen Werte auf 0 zurueckgesetzt werden  */
				// will be done vegResetFlag-Handler above
				//			for (i=1;i<=nrows;i++)
				//			  for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				//				if (zn[i][j]>nodata) summengrid[i][j] = 0.0;

				if (dateiname) delete[] dateiname; dateiname=NULL;
				if (tempstr) delete[] tempstr;tempstr=NULL;
			}
			if (*vegResetFlag==1)
			{
				// reset the values in vegSumPointer and vegCount
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							vegSumPointer[i][j]=0;
						}
						else
						{
							vegSumPointer[i][j]=nodata;
						}
					}
				}
				*vegResetFlag = 0;
			}
		}

		if ((writegrid==3)&&((minf*minute+hf*stunde+df*tag+mf*monat+jf*jahr)<endtime))
		{  }
		else if (
				((writegrid==1)||(writegrid==2)||(writegrid==3)) || // end of model period
				( ((writegrid==4)||(writegrid==5))&&((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) ) || // // end of day
				( (writegrid==6)&&((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) && (tag==mon[monat]) ) || // end of month
				( (writegrid==7)&&((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) && (tag==31) && (monat==12) ) || // end of year
				(writegrid==8) // writing the grid after reading and filling it
				)
		{
			tempstr=newstring(strlen(name)+20,1);
			/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
			if (strrchr(name,dirzeichen)!=NULL)
			strcpy(tempstr,strrchr(name,dirzeichen)+1); else strcpy(tempstr,name);
			dateiname = newstring(strlen(name)+20,1);
			if ((writegrid==1)||(writegrid==3)||(writegrid==4)) 
			{
				strcpy(dateiname,name);
			} else if ((writegrid==2)||(writegrid==5)||(writegrid==6)||(writegrid==7))	{
				if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), "%s", ".");
				p = index(tempstr,".");
				//        if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;

				// the first call will always jump over this if statement since lastYear, lastMonth, lastDay and lastHour are always set to 0 during initialization of the gridinfo array
				// thus we make sure that the first call will always be compatible with the old file extension naming schema (only 3 digits for the hour)
				if ((jahr == gridinfo[idx].lastYear) && (monat == gridinfo[idx].lastMonth) && (tag == gridinfo[idx].lastDay) && (stunde == gridinfo[idx].lastHour) && (minute == gridinfo[idx].lastMinute))
				{
					gridinfo[idx].subIntervalCnt++; // increment sub interval counter and use it for creating the grid name
					if (WriteNetcdfGrids)
					{
//						sprintf(dateiname, "%.*s_%04ld", p, tempstr, jahr);
						sprintf(dateiname, "%.*s_%s", p, tempstr, "timeseries");
					}
					else
					{
						sprintf(dateiname, "%.*s_%04ld%02ld%02ld.%02ld%02ld_%.03ld", p, tempstr, jahr, monat, tag, stunde, minute, gridinfo[idx].subIntervalCnt);
					}
				}
				else
				{
					// initializing date time counters when this is the first time the writegridtofiel function is called for a grid
					if (gridinfo[idx].lastYear == 0)
					{
						gridinfo[idx].lastYear = jahr;
						gridinfo[idx].lastMonth = monat;
						gridinfo[idx].lastDay = tag;
						gridinfo[idx].lastHour = stunde;
						gridinfo[idx].lastMinute = minute;
					}
					gridinfo[idx].subIntervalCnt = 0;
					if (WriteNetcdfGrids)
					{
//						sprintf(dateiname, "%.*s_%04ld", p, tempstr, jahr);
						sprintf(dateiname, "%.*s_%s", p, tempstr, "timeseries");
					}
					else
					{
						sprintf(dateiname, "%.*s_%04ld%02ld%02ld.%02ld%02ld", p, tempstr, jahr, monat, tag, stunde, minute);
					}
				}
				gridinfo[idx].lastYear = jahr;
				gridinfo[idx].lastMonth = monat;
				gridinfo[idx].lastDay = tag;
				gridinfo[idx].lastHour = stunde;
				gridinfo[idx].lastMinute = minute;

				/* Pfad wird wieder angehaengt */
				if (strrchr(name,dirzeichen)!=NULL)
				{
					strcpy(tempstr,dateiname);
					p=(int)(strrchr(name,dirzeichen)-name)+1;
					sprintf(dateiname,"%.*s%s",p,name,tempstr);
				}
			} else if (writegrid==8) {
				sprintf(dateiname,"%s.filled",name);
				gridinfo[idx].writegrid = 0; // reset writegrid in order to prevent writing after each time step. code 8 is for initial writing only (after reading, missing values are filled and then the grid is written bak
			}
			/* das im Speicher an der Adresse **grid gehaltene Array wird in ein Grid geschrieben */
			double** tempgrid;
			if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, grid); else tempgrid = grid;
			if (mpi_rank==0)
			{
				printgridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids,jahr,monat,tag,stunde,minute,hour_offset);
				if (!((writegrid == 1) || (writegrid == 3) || (writegrid == 4)) && ((minf*minute+hf*stunde + df*tag + mf*monat + jf*jahr)>=endtime))
				{ // file extension is different from what is defined in the control file, so we will write that file in any case
					strcpy(dateiname, name);
					printgridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
			}
			if (dateiname) free(dateiname); dateiname=NULL;
			if (tempstr) free(tempstr); tempstr=NULL;
			if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
		}
		/* writing minimum grid */
		if ((minwritegrid == 1) && (!doNotPrintAverageYet) && ((minf*minute+hf*stunde + df*tag + mf*monat + jf*jahr) >= endtime))
		{
			/* Namen des minimum-grids zusammensetzen */
			tempstr=newstring(strlen(minname)+20,1);
			/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
			if (strrchr(minname,dirzeichen)!=NULL)
			strcpy(tempstr,strrchr(minname,dirzeichen)+1); else strcpy(tempstr,minname);
			dateiname = newstring(strlen(minname)+20,1);
			if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
			sprintf(dateiname,"%.*s.min",index(tempstr,"."),tempstr); /* Summengrids fuer ges.Zeit enden auf "xxxxyyyy.sum" */
			/* Pfad wird wieder angehaengt */
			if (strrchr(minname,dirzeichen)!=NULL)
			{
				strcpy(tempstr,dateiname);
				p=(int)(strrchr(minname,dirzeichen)-minname)+1;
				sprintf(dateiname,"%.*s%s",p,minname,tempstr);
			}
			/* das im Speicher an der Adresse **mingrid gehaltene Array wird in ein Grid geschrieben */
			double **tempgrid;
			if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, mingrid); else tempgrid = mingrid;
			if (mpi_rank == 0)
			{
				printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
			if (dateiname) free(dateiname);dateiname=NULL;
			if (tempstr) free(tempstr);tempstr=NULL;
		}

		/* writing maximum grid */
		if ((maxwritegrid == 1) && (!doNotPrintAverageYet) && ((minf*minute+hf*stunde + df*tag + mf*monat + jf*jahr) >= endtime))
		{
			/* Namen des maximum-grids zusammensetzen */
			tempstr=newstring(strlen(maxname)+20,1);
			/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
			if (strrchr(maxname,dirzeichen)!=NULL)
			strcpy(tempstr,strrchr(maxname,dirzeichen)+1); else strcpy(tempstr,maxname);
			dateiname = newstring(strlen(maxname)+20,1);
			if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
			sprintf(dateiname,"%.*s.max",index(tempstr,"."),tempstr); /* Summengrids fuer ges.Zeit enden auf "xxxxyyyy.sum" */
			/* Pfad wird wieder angehaengt */
			if (strrchr(maxname,dirzeichen)!=NULL)
			{
				strcpy(tempstr,dateiname);
				p=(int)(strrchr(maxname,dirzeichen)-maxname)+1;
				sprintf(dateiname,"%.*s%s",p,maxname,tempstr);
			}
			/* das im Speicher an der Adresse **maxgrid gehaltene Array wird in ein Grid geschrieben */
			double **tempgrid;
			if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, maxgrid); else tempgrid = maxgrid;
			if (mpi_rank == 0)
			{
				printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
			if (dateiname) free(dateiname);dateiname=NULL;
			if (tempstr) free(tempstr);tempstr=NULL;
		}

		/* Schreiben der Summengrids */
		if (((summenwritegrid == 1) || (summenwritegrid == 5)) && (!doNotPrintAverageYet) && ((minf*minute + hf*stunde + df*tag + mf*monat + jf*jahr) >= endtime))
		{
			/* Schreiben eines Grids am Ende des Modellaufes */
			if (summenwritegrid==5)
			{
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							summengrid[i][j] /= (double)*anzahl;
						}
					}
				}
			}
			/* Namen des Jahresgrids zusammensetzen */
			tempstr=newstring(strlen(summenname)+20,1);
			/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
			if (strrchr(summenname,dirzeichen)!=NULL)
			strcpy(tempstr,strrchr(summenname,dirzeichen)+1); else strcpy(tempstr,summenname);
			dateiname = newstring(strlen(summenname)+20,1);
			if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
			p = index(tempstr,".");
			if (summenwritegrid == 1)
			sprintf(dateiname,"%.*s.sum",p,tempstr); /* Summengrids fuer ges.Zeit enden auf "xxxxyyyy.sum" */
			else
			sprintf(dateiname,"%.*s.mit",p,tempstr); /* Mittelwertgrids fuer ges.Zeit enden auf "xxxxyyyy.mit" */
			/* Pfad wird wieder angehaengt */
			if (strrchr(summenname,dirzeichen)!=NULL)
			{
				strcpy(tempstr,dateiname);
				p=(int)(strrchr(summenname,dirzeichen)-summenname)+1;
				sprintf(dateiname,"%.*s%s",p,summenname,tempstr);
			}
			/* das im Speicher an der Adresse **summengrid gehaltene Array wird in ein Grid geschrieben */
			double **tempgrid;
			if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, summengrid); else tempgrid = summengrid;
			if (mpi_rank == 0)
			{
				printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
			/* nun muessen die im Summengrid gehaltenen Werte auf 0 zurueckgesetzt werden  */
			*anzahl = 0;
			for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
			{
				for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
					{
						summengrid[i][j] = 0.0;
					}
				}
			}
			if (dateiname) free(dateiname);dateiname=NULL;
			if (tempstr) free(tempstr);tempstr=NULL;
		}
		else
		{
			if (((summenwritegrid == 2) || (summenwritegrid == 6)) && (!doNotPrintAverageYet) && ((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) && (tag == 31) && (monat == 12))
			{ /* Schreiben eines Summen- oder Mittelwertgrids am Jahresende */
				if (summenwritegrid == 6)
				{
					for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
					{
						for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						{
							if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
							{
								summengrid[i][j] /= (double)*anzahl;
							}
						}
					}
				}
				/* Namen des Jahresgrids zusammensetzen */
				tempstr=newstring(strlen(summenname)+20,1);
				/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
				if (strrchr(summenname,dirzeichen)!=NULL)
				strcpy(tempstr,strrchr(summenname,dirzeichen)+1); else strcpy(tempstr,summenname);
				dateiname = newstring(strlen(summenname)+20,1);
				if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
				p = index(tempstr,".");
				//        if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;
				if (summenwritegrid == 2)
				sprintf(dateiname,"%.*s_%04ld00.sum",p,tempstr,jahr); /* Summengrids fuers Jahr enden auf "xxxxjj00.sum" */
				else
				sprintf(dateiname,"%.*s_%04ld00.mit",p,tempstr,jahr); /* Mittelwertgrids fuers Jahr enden auf "xxxxjj00.mit" */
				/* Pfad wird wieder angehaengt */
				if (strrchr(summenname,dirzeichen)!=NULL)
				{
					strcpy(tempstr,dateiname);
					p=(int)(strrchr(summenname,dirzeichen)-summenname)+1;
					sprintf(dateiname,"%.*s%s",p,summenname,tempstr);
				}
				/* das im Speicher an der Adresse **summengrid gehaltene Array wird in ein Grid geschrieben */
				double **tempgrid;
				if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, summengrid); else tempgrid = summengrid;
				if (mpi_rank == 0)
				{
					printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
				}
				if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
				/* nun muessen die im Summengrid gehaltenen Werte auf 0 zurueckgesetzt werden  */
				*anzahl = 0;
				for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
				{
					for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
					{
						if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
						{
							summengrid[i][j] = 0.0;
						}
					}
				}
				if (dateiname) free(dateiname);dateiname=NULL;
				if (tempstr) free(tempstr);tempstr=NULL;
			}
			else
			{
				if (((summenwritegrid == 3) || (summenwritegrid == 7)) && (!doNotPrintAverageYet) && ((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) && (tag == mon[(int)monat]))
				{
					/* Schreiben eines SUmmen- oder Mittelwertgrids am Monatsende */
					if (summenwritegrid == 7)
					{
						for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
						{
							for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
							{
								if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
								{
									summengrid[i][j] /= (double)*anzahl;
								}
							}
						}
					}
					/* Namen des Monatsgrids zusammensetzen */
					tempstr=newstring(strlen(summenname)+20,1);
					/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
					if (strrchr(summenname,dirzeichen)!=NULL)
						strcpy(tempstr,strrchr(summenname,dirzeichen)+1); else strcpy(tempstr,summenname);
					dateiname = newstring(strlen(summenname)+20,1);
					if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
					p = index(tempstr,".");
					//          if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;
					if (summenwritegrid == 3)
						// $kf 2017-02-08: 4 digit year
						// sprintf(dateiname,"%.*s_%02ld%02ld.sum",p,tempstr,jahr%100,monat); /* Monatssummen als xxxxjjmm.sum */
						sprintf(dateiname,"%.*s_%04ld%02ld.sum",p,tempstr,jahr,monat); /* Monatssummen als xxxxjjmm.sum */
					else
						// $kf 2017-02-08: 4 digit year
						// sprintf(dateiname,"%.*s_%02ld%02ld.mit",p,tempstr,jahr%100,monat); /* Monatsmittelwerte als xxxxjjmm.mit */
						sprintf(dateiname,"%.*s_%04ld%02ld.mit",p,tempstr,jahr,monat); /* Monatsmittelwerte als xxxxjjmm.mit */
					/* Pfad wird wieder angehaengt */
					if (strrchr(summenname,dirzeichen)!=NULL)
					{
						strcpy(tempstr,dateiname);
						p=(int)(strrchr(summenname,dirzeichen)-summenname)+1;
						sprintf(dateiname,"%.*s%s",p,summenname,tempstr);
					}
					/* das im Speicher an der Adresse **summengrid gehaltene Array wird in ein Grid geschrieben */
					double **tempgrid;
					if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, summengrid); else tempgrid = summengrid;
					if (mpi_rank == 0)
					{
						printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
					}
					if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
					/* nun muessen die im Summengrid gehaltenen Werte auf 0 zurueckgesetzt werden  */
					*anzahl = 0;
					for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
					{
						for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
						{
							if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
							{
								summengrid[i][j] = 0.0;
							}
						}
					}
					if (dateiname) free(dateiname);dateiname=NULL;
					if (tempstr) free(tempstr);tempstr=NULL;
				}
				else
				{
					if (((summenwritegrid == 4) || (summenwritegrid == 8)) && (!doNotPrintAverageYet) && ((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440))
					{
						/* Schreiben eines Summen- oder Mittelwertgrids am Tagesende (wenn InV < 1440) */
						if (summenwritegrid==8)
						{
							for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
							{
								for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
								{
									if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
									{
										summengrid[i][j] /= (double)*anzahl;
									}
								}
							}
						}
						/* Namen des tagesgrids zusammensetzen */
						tempstr=newstring(strlen(summenname)+20,1);
						/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
						if (strrchr(summenname,dirzeichen)!=NULL)
							strcpy(tempstr,strrchr(summenname,dirzeichen)+1); else strcpy(tempstr,summenname);
						dateiname = newstring(strlen(summenname)+20,1);
						if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), ".");
						p = index(tempstr,".");
						//            if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;
						if (summenwritegrid == 4)
							sprintf(dateiname,"%.*s_%04ld%02ld.%02lds",p,tempstr,jahr%100,monat,tag); /* Tagessummen als       xxxxjjMM.dds */
						else
							sprintf(dateiname,"%.*s_%04ld%02ld.%02ldm",p,tempstr,jahr%100,monat,tag); /* Tages-Mittelwerte als xxxxjjMM.ddm */
						/* Pfad wird wieder angehaengt */
						if (strrchr(summenname,dirzeichen)!=NULL)
						{
							strcpy(tempstr,dateiname);
							p=(int)(strrchr(summenname,dirzeichen)-summenname)+1;
							sprintf(dateiname,"%.*s%s",p,summenname,tempstr);
						}
						/* das im Speicher an der Adresse **summengrid gehaltene Array wird in ein Grid geschrieben */
						double **tempgrid;
						if ((gather > 0) && (mpi_size>1)) tempgrid = MPIGatherIntoNewGrid(idx, summengrid); else tempgrid = summengrid;
						if (mpi_rank==0)
						{
							printsummengridfile(dateiname, tempgrid, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
						}
						if ((gather > 0) && (mpi_size>1)) free_matrix(tempgrid, 1, nrows, 1, ncols);
						/* nun muessen die im Summengrid gehaltenen Werte auf 0 zurueckgesetzt werden  */
						*anzahl = 0;
						for (int i = gridinfo[idx].nrl; i <= gridinfo[idx].nrh; i++)
						{
							for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
							{
								if ((zn[i][j]>nodata)&&((MPIREG[i][j]==mpi_rank)||(gather==0)))
								{
									summengrid[i][j] = 0.0;
								}
							}
						}
						if (dateiname) free(dateiname);dateiname=NULL;
						if (tempstr) free(tempstr);tempstr=NULL;
					}
				}
			}
		}
	}
}

int modelobject::readstackfromfile(char *dateiname, double ****p_stack, bool MPIFullGrid)
/* !!! es findet keine Ueberpruefung auf nodata-Werte statt !!! */
{
	int i, j, levels;
	FILE *datei;
	int format;

#ifdef _MPI
	MPI_File MPIfile;
	MPI_Status MPIstatus;
	MPI_Offset MPIoffset;
#endif

	signal(SIGINT, SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK, SIG_IGN);
#endif

	char *completedname = new char[_stringsize_ + 1];
	sprintf(completedname, "%s", dateiname);
	ExtendFileName(DefaultInputDirectory, completedname);

	if (mpi_rank == 0)
	{

		/* Oeffnen des grids und lesen der Parameter */
		if ((datei = fopen(completedname, "rb")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << completedname << "[end-entry]";
			throw(21);
		}

		(*WSstream) << "reading " << completedname << " for " << instanzname << "..." << "[end-entry]";

		char FormatTest[100];
		fread(FormatTest, sizeof(char), 5, datei);

		FormatTest[5] = '\0';

#ifdef _NETCDF
		std::string fourletters = ((string)FormatTest).substr(0, 4);
		if (strlen(FormatTest) > 2 && ((fourletters.find("HDF") == 1) || (fourletters.find("CDF") == 0)))
		{
			format = 10; // netCDF
			fclose(datei); // will be opened again in netCDF reader internally
			// TODO: Is this correct, or should we read stuff from the netcdf?

			zz[1] = ncols;
			zz[2] = nrows;
			zz[3] = xll;
			zz[4] = yll;
			zz[5] = cellsize;
			zz[6] = nodata;
			levels = -1; // must be taken from variable t's dimension in netCDF
		}
		else
#endif
			if (strcmp(FormatTest, "stack") == 0)
			{
				format = 1;  //ascii
				fclose(datei);
				// reopen in Textmode
				datei = fopen(completedname, "rt");
				char text[30];
				// read header
				fscanf(datei, "%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf", text, &(zz[1]), text, &(zz[2]), text, &(zz[3]), text, &(zz[4]), text, &(zz[5]), text, &(zz[6]), text, &(zz[7]), text, &(zz[8]), text, &(zz[9]), text, &(zz[10]), text, &(zz[11]), text, &(zz[12]));
				levels = (int)floor(zz[12]);
			}
			else
			{
				fseek(datei, 0, 0);
				if ((fread(&zz[1], sizeof(double), 12, datei)) != 12)
				{
					(*WSstream) << "FATAL: error reading from " << completedname << "[end-entry]";
					throw(22);
				}
				// analyse binary format: is it a double or a double grid?
				if (zz[6] == -9999) {
					format = 2; // double
				}
				else {
					fseek(datei, 0, 0);
					float fzz[13];
					fread(&fzz[1], sizeof(float), 12, datei);
					if (fzz[6] == -9999)
					{
						format = 0; // float
						for (i = 1; i <= 12; i++)
						{
							zz[i] = (double)fzz[i];
						}
					}
					else
					{
						(*WSstream) << "FATAL: Error: could not recognize file format (ascii, single float or double float) for stack " << completedname << "[end-entry]";
						throw(22);
					}
				}
				levels = (int)floor(zz[12]);
			}

		/* Bevor diese Routine aufgerufen werden kann, muessen die Werte
	der Gridbeschreibenden Daten bekannt sein (werden als Parameter an den
	Constructor von modelobject bzw. an dessen Erben uebergeben */

		if ((ncols != (int)floor(zz[1])) ||
			(nrows != (int)floor(zz[2])) ||
			((int)floor(xll) != (int)floor(zz[3])) ||
			((int)floor(yll) != (int)floor(zz[4])) ||
			(cellsize != zz[5]) || (nodata != zz[6]))
		{
			(*WSstream) << "FATAL: " << completedname << " has not the specified grid measures (in " << instanzname << ")\n"
				<< (strlen(gridinfo[1].name) > 35 ? &(gridinfo[1].name[strlen(gridinfo[1].name) - 35]) : gridinfo[1].name)
				<< " " << (strlen(completedname) > 35 ? &(completedname[strlen(completedname) - 35]) : completedname)
				<< "\n" << ncols << " " << zz[1] << "\n" << nrows << " " << zz[2] << "\n" << xll
				<< " " << zz[3] << "\n" << yll << " " << zz[4] << "\n" << cellsize << " " << zz[5]
				<< "\n" << nodata << " " << zz[6] << "[end-entry]";
			if (fclose(datei) == EOF)
			{
				(*WSstream) << "FATAL: error closing stack (1) " << completedname << "[end-entry]";
				throw 23;
			}
			throw 24;
		}
	}
#ifdef _MPI
	// all the other ranks doesn't know yet about the file name, the format and the header infos, so we need to broadcast these data now
	MPI_Bcast(&(completedname[0]), _stringsize_ + 1, MPI_CHAR, 0,MPI_COMM_WORLD);
	MPI_Bcast(&format, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(zz[1]), 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&levels, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	int sr = mpistartrow;
	int er = mpiendrow;
	// only a stripe that is defined for the actual rank needs to be allocated
	//	sr = (mpistartrow>1) ? mpistartrow - 1 : mpistartrow;
	//	er = (mpiendrow<nrows) ? mpiendrow + 1 : mpiendrow;
	if ((mpi_rank == 0) && ((format == 1) || (format == 10) || (MPIFullGrid))) // ascii and netCDF or special Stacks that need to be present in rank 0: read the complete stack in instance 0 and scatter it to the other instances
	{
		sr = 1;
		er = nrows;
	}

#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if (!MPIFullGrid && ((format == 0) || (format == 2))) // single or double formats can be read with MPI functions in parallel if not the entire grid is required in rank 0. Otherwise the grid is read in rank 0 only and scattered to the other ranks
	{
		if (mpi_rank == 0)
		{
			if (fclose(datei) == EOF)
			{
				(*WSstream) << "FATAL: error closing grid " << completedname << "[end-entry]";
				throw(23);
			}
		}
		// open grid for MPI collective I/O
		MPI_File_open(MPI_COMM_WORLD, completedname, MPI_MODE_RDWR, MPI_INFO_NULL, &MPIfile);
		// set Offset where to start reading from
		if (format == 0)
		{ // float: 
			MPIoffset = (12 + ncols * (sr - 1)*(levels + 1)) * sizeof(float);
		}
		else
		{ // double
			MPIoffset = (12 + ncols * (sr - 1)*(levels + 1)) * sizeof(double);
		}
		MPI_File_seek(MPIfile,MPIoffset, MPI_SEEK_SET);
	}
#endif

	// now allocating the correct amount of memory for the stack
	if (format != 10)
	{ // for netCDF, memory can be allocated after reading the first layer only, since number of levels is returned as one of the parameters by reference
		*p_stack = drei_d_matrix(sr, er, 1, ncols, levels, nrows);
	}
	//cout << "mpi rank: " << mpi_rank << " reserved stack from row " << sr << " to " << er << ". This are approx. bytes: "<< (er-sr+1)*(ncols+1)*8*(levels+1) <<endl;
	/* Einlesen des Grids */
	double *dBuffer;
	float *fBuffer;
	unsigned long ulBufPos = 0;

	// in order to reduce the number of write accesses, a buffer is reserved and filled and then written as a single big chunk instead of potentially thousands of small chunks
	if (format == 2)
	{
		dBuffer = new double[(er - sr + 1)*ncols*(levels + 1)];
	}

	if (format == 0)
	{
		fBuffer = new float[(er - sr + 1)*ncols*(levels + 1)];
	}

	if (!MPIFullGrid && ((format == 0) || (format == 2)))
	{
#ifdef _MPI
		if (format == 0)
		{
			MPI_Status s;
			MPI_File_read(MPIfile,&(fBuffer[0]), (er - sr + 1)*ncols*(levels + 1), MPI_FLOAT,&s);
			for (i = sr; i <= er; i++)
			{
				ulBufPos = 0;
				for (j = 1; j <= ncols; j++)
				{
					for (int k = 0; k <= levels; k++)
					{
						(*p_stack)[i][j][k] = (double)fBuffer[(levels + 1)*(j - 1 + (i - sr)*ncols) + k];
						//						(*p_stack)[i][j][k] = (double)fBuffer[ulBufPos];
						//						ulBufPos++;
					}
				}
			}
		}
		else
		{
			MPI_Status s;
			MPI_File_read(MPIfile,&(dBuffer[0]), (er - sr + 1)*ncols*(levels + 1), MPI_DOUBLE,&s);
			for (i = sr; i <= er; i++)
			{
				memmove((*p_stack)[i][1], &(dBuffer[ulBufPos]), ncols*(levels + 1) * sizeof(double));
				ulBufPos += ncols * (levels + 1);
				//				for (j = 1; j <= ncols; j++)
				//				{
				//					for (int k = 0; k <= levels; k++)
				//					{
				//						(*p_stack)[i][j][k] = dBuffer[ulBufPos];
				//						ulBufPos++;
				//					}
				//				}
			}
		}
#else
		// normal C file I/O, 
		if (format == 0) // single precision format
		{
			for (i = sr; i <= er; i++)
			{
				ulBufPos = 0;
				if ((int)fread(&(fBuffer[0]), sizeof(float), (levels + 1)*ncols, datei) != (levels + 1)*ncols)
				{
					(*WSstream) << "FATAL: error reading float stack " << completedname << "at row " << i << ", column " << j << "[end-entry]";
					throw 22;
				}
				for (j = 1; j <= ncols; j++)
				{
					for (int k = 0; k <= levels; k++)
					{
						(*p_stack)[i][j][k] = (double)fBuffer[(levels + 1)*(j - 1) + k];
						//						(*p_stack)[i][j][k] = (double)fBuffer[ulBufPos];
						//						ulBufPos++;
					}
				}
			}
		}
		else // double format
		{
			for (i = sr; i <= er; i++)
			{
				if ((int)fread(&(*p_stack)[i][1][0], sizeof(double), (levels + 1)*ncols, datei) != (levels + 1)*ncols)
				{
					(*WSstream) << "FATAL: error reading double stack " << completedname << "at row " << i << "[end-entry]";
					throw 22;
				}
			}
		}
#endif
	}
	else // MPIFullGrid or ascii format or netCDF: read in rank 0 only
	{
		if (mpi_rank == 0)
		{
			if (format <= 2)
			{
				for (i = 1; i <= nrows; i++)
				{
					if (format == 0) // format=0: single
					{
						if ((int)fread(&(fBuffer[0]), sizeof(float), (levels + 1)*ncols, datei) != (levels + 1)*ncols)
						{
							(*WSstream) << "FATAL: error reading float stack " << completedname << "at row " << i << "[end-entry]";
							throw 22;
						}
						for (j = 1; j <= ncols; j++)
						{
							for (int k = 0; k <= levels; k++)
							{
								(*p_stack)[i][j][k] = (double)fBuffer[(levels + 1)*(j - 1) + k];
							}
						}
					}
					else if (format == 2) // format 2: double
					{
						if ((int)fread(&(*p_stack)[i][1][0], sizeof(double), (levels + 1)*ncols, datei) != (levels + 1)*ncols)
						{
							(*WSstream) << "FATAL: error reading double stack " << completedname << "at row " << i << "[end-entry]";
							throw 22;
						}
					}
					else if (format == 1)
					{
						for (j = 1; j <= ncols; j++)
						{
							for (int k = 0; k <= levels; k++)
							{
								if (fscanf(datei, "%lf", &(*p_stack)[i][j][k]) != 1)
								{
									(*WSstream) << "FATAL: error reading ascii stack " << completedname << "at row " << i << ", column " << j << ", level " << k << "[end-entry]";
									throw 22;
								}
							}
						}
					}
				}
			}
			else // format==10, i.e. netCDF
			{
				// Todo read netCDF here timestep nby timestep
				double *stackbuf = NULL; // will be allocated in method netcdfIn->readTimestep() and must be deleted by the caller!
				int nc_rows, nc_cols = 0;
				double nc_xll, nc_yll, nc_cellsize, nc_nodata;
				int nc_levels;
				int k = 0;
				// read first layer. Return values will give us also the number of levels
				netcdfIn->readTimestep(true, stackbuf, completedname, NULL, nc_rows, nc_cols, nc_xll, nc_yll, nc_cellsize, nc_nodata, nc_levels);
				levels = nc_levels - 1;

				if ((ncols != (int)floor(nc_cols)) ||
					(nrows != (int)floor(nc_rows)) ||
					((int)floor(xll) != (int)floor(nc_xll)) ||
					((int)floor(yll) != (int)floor(nc_yll)) ||
					(cellsize != nc_cellsize))
				{
					(*WSstream) << "FATAL: \nError: " << completedname << " has not the specified stack measures (in " << instanzname << ")\n"
						<< (strlen(gridinfo[1].name) > 125 ? &(gridinfo[1].name[strlen(gridinfo[1].name) - 125]) : gridinfo[1].name)
						<< " <--vs.--> " << (strlen(completedname) > 125 ? &(completedname[strlen(completedname) - 125]) : completedname)
						<< "\n expected ncols: " << ncols << " <-> read from stack: " << nc_cols
						<< "\n expected nrows: " << nrows << " <-> read from stack: " << nc_rows
						<< "\n expected xll:   " << xll << " <-> read from stack: " << nc_xll
						<< "\n expected yll:   " << yll << " <-> read from stack: " << nc_yll
						<< "\n expected cs:" << cellsize << " <-> read from stack: " << nc_cellsize
						<< "\n expected levels:" << levels << " <-> read from stack: " << nc_levels
						<< "\n\n netCDF stacks that are not meteorologic input data must still match the model domain structure as given by the first grid and the soil table" << "[end-entry]";
					if (fclose(datei) == EOF)
					{
						(*WSstream) << "FATAL: error closing stack (1) " << completedname << "[end-entry]";
						throw 23;
					}
					throw 24;
				}

				*p_stack = drei_d_matrix(sr, er, 1, ncols, levels, nrows);

				//now map the buffer data of the first layer into the stack 
				int count = 0; // start copying the buffer data into the structured 3D-matrix
				for (j = 1; j <= ncols; j++) {
					for (i = 1; i <= nrows; i++) {
						(*p_stack)[i][j][k] = (stackbuf[count] == nc_nodata) ? nodata : stackbuf[count];
						count++;
					}
					//count += nrows - (er - sr) - 1;
				}
				k++;

				// read all the remainig levels and map them to the respective layer in the result stack
				while (k <= levels)
				{
					// read next timestep
					netcdfIn->readTimestep(true, stackbuf, completedname, NULL, nc_rows, nc_cols, nc_xll, nc_yll, nc_cellsize, nc_nodata, nc_levels, k);
					// int count = sr - 1;
					int count = 0; // start copying the buffer data into the structured 3D-matrix
					for (j = 1; j <= ncols; j++) {
						for (i = 1; i <= nrows; i++) {
							(*p_stack)[i][j][k] = (stackbuf[count] == nc_nodata) ? nodata : stackbuf[count];
							count++;
						}
						//count += nrows - (er - sr) - 1;
					}
					k++;
				}
				delete[] stackbuf; // now deleting the temporary buffer that was allocated in netcdfIn->readTimestep()
			}
		}
	}

#ifdef _MPI

#ifdef _NETCDF
	if (format == 10)
	{
		MPI_Bcast(&levels, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
		if (mpi_rank > 0) // for netCDF memory was allocated later because level count was not available earlier. This was done in rank 0 only, so we need to allocate memeory for other ranks now
		{
			*p_stack = drei_d_matrix(sr, er, 1, ncols, levels, nrows);
		}
	}
#endif

	if (!MPIFullGrid && ((format == 0) || (format == 2)))
	{
		MPI_File_close(&MPIfile);
	}
	else
	{ // was only read in rank 0 because ascii format and/or MPIFullGrid
		if (mpi_rank == 0)
		{
			if (format != 10)
			{
				if (fclose(datei) == EOF)
				{
					(*WSstream) << "FATAL: error closing stack (2) " << completedname << "[end-entry]";
					throw(23);
				}
			}
		}
		MPIScatterStack(0, *p_stack,levels);
		if (mpi_rank > 0) // setting all nodata cells to nodata since these values will not be scattered
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if (MPIREG[i][j] != mpi_rank)
					{
						for (int k = 1; k <= levels; k++)
						{
							(*p_stack)[i][j][k]=nodata;
						}
					}
				}
			}
		}

		if ((!MPIFullGrid)&&(mpi_rank==0)) // re-allocating stack for rank 0
		{
			double ***newstack;
//			(*WSstream) << "DEBUG: calling  drei_d_matrix for " << completedname << " with mpistartrow="<<mpistartrow<<" and mpiendrow="<<mpiendrow<<"[end-entry]";
			newstack = drei_d_matrix(mpistartrow, mpiendrow, 1, ncols, levels, nrows);
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j=1;j<=ncols;j++)
				{
					for (int k = 0; k <= levels; k++)
					{
						newstack[i][j][k] = (*p_stack)[i][j][k];
					}
				}
			}
			double ***tempstack = (*p_stack);
			(*p_stack) = newstack;
//			(*WSstream) << "DEBUG: calling  free_drei_d_matrix for " << completedname << "[end-entry]";
			free_drei_d_matrix(tempstack, sr, er, 1, ncols, levels, nrows);
//cout << "mpi rank: " << mpi_rank << " releasing stack from row " << sr << " to " << er << " (" << (er - sr + 1)*(ncols + 1) * 8 * (levels + 1) << "bytes) and allocating new from " << mpistartrow << " to " << mpiendrow << ". This are approx. bytes: " << (mpiendrow - mpistartrow + 1)*(ncols + 1) * 8 * (levels + 1) << endl;
		}
	}
#else
	if (format != 10) // netCDF files will not be opened directly but within the netCDF routines, where they are closed as well, so we do not need a closing command here
	{
		if (fclose(datei) == EOF)
		{
			(*WSstream) << "FATAL: error closing stack (3) " << completedname << "[end-entry]";
			throw(23);
		}
	}
#endif		

	if (format == 0)	delete[] fBuffer;
	if (format == 2)	delete[] dBuffer;

	// debug
//	char dateinametemp[500];
//	sprintf(dateinametemp, "%s_reprinted", dateiname);
//	printstackfile(dateinametemp, *p_stack, levels, false, true);

	delete[] completedname;
	SetSignalDefaults();
	return levels;
}

void modelobject::printstackfile(char *dateiname, double ***stack, int levels, bool WriteASCIIGrids, bool bDoubleMode, bool WriteNetcdfGrids)
{
	// if called from writestacktofile, then following rules apply for MPI:
	// WriteASCIIGrids==true: only rank 0 writes the file - all ranks data were gathered before calling this routine (i.e. if calling this routine separately, gather must be called before when in MPI mode and format is ascii)
	// all binary formats will be written in parallel by MPI functions.
	FILE *datei;
	int i,j;
#ifdef _MPI
	MPI_File MPIfile;
	MPI_Status MPIstatus;
	MPI_Offset MPIoffset;
#endif


	signal(SIGINT,SIG_IGN);
#if defined _BCC_
	signal(SIGBREAK,SIG_IGN);
#endif
	char *completedname = new char[_stringsize_+1];
	sprintf(completedname,"%s",dateiname);
	ExtendFileName(DefaultOutputDirectory,completedname);
	//  delete[] completedname;

#ifdef _NETCDF

	if (!WriteNetcdfGrids)
	{
#endif
		//  (*WSstream) << "Ctrl-C disabled (writing stack file..) ");
		if (WriteASCIIGrids)
		{
			if (mpi_rank == 0) // will be true for non-MPI version as well
			{
				if ((datei = fopen(completedname, "wt")) == NULL)
				{
					(*WSstream) << "FATAL: error opening " << completedname << "[end-entry]";
					throw(26);
				}
			}
		}
		else
		{
#ifdef _MPI
			// open grid for MPI collective I/O
			MPI_File_open(MPI_COMM_WORLD, completedname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&MPIfile);
			// set Offset where to start write to (rank 0: header must be written first
			if (mpi_rank > 0)
			{
				// setting the base offset for the beginning of the actual row only (not yet for the column. This will be done in the write loop itself
				if (bDoubleMode == false)
				{ // float: 
					MPIoffset = (12 + ncols*(levels + 1)*(mpistartrow - 1))*sizeof(float);
				}
				else
				{ // double
					MPIoffset = (12 + ncols*(levels + 1)*(mpistartrow - 1))*sizeof(double);
				}
				MPI_File_seek(MPIfile,MPIoffset, MPI_SEEK_SET);
			}
			else
			{
				MPI_File_seek(MPIfile,0, MPI_SEEK_SET);
			}

#else 
			// here, no mpi_rank is required since this branch will be compiled for non-mpi-version only
			if ((datei = fopen(completedname, "wb")) == NULL)
			{
				(*WSstream) << "FATAL: error opening " << completedname << "[end-entry]";
				throw(26);
			}
#endif
		}

		zz[1] = ncols; zz[2] = nrows; zz[3] = xll; zz[4] = yll; zz[5] = cellsize; zz[6] = nodata;
		zz[12] = levels;
		zz[7] = zz[8] = zz[9] = zz[10] = zz[11] = 0.0;

		if (WriteASCIIGrids && (mpi_rank == 0)) {
			fprintf(datei, "stack_ncols\t%lf\nnrows\t%lf\nxll\t%lf\nyll\t%lf\ncs\t%lf\nnodata\t%lf\ncnt\t%lf\nmin\t%lf\nmax\t%lf\nsum\t%lf\naverage\t%lf\nnumlevels\t%lf\n",
				zz[1], zz[2], zz[3], zz[4], zz[5], zz[6], zz[7], zz[8], zz[9], zz[10], zz[11], zz[12]);
		}
		else if (bDoubleMode)
		{
#ifdef _MPI
			if (mpi_rank == 0)
			{
				MPI_File_seek(MPIfile,0, MPI_SEEK_SET);
				MPI_Status s;
				MPI_File_write(MPIfile, &(zz[1]), 12, MPI_DOUBLE,&s);
				MPI_File_seek(MPIfile, MPIoffset, MPI_SEEK_SET);
			}
#else
			// this branch requires no mpi_rank checking because it is compiled when _MPI is not deefined
			if ((fwrite(&zz[1], sizeof(double), 12, datei)) != 12) writeerror(completedname, 18);
#endif
		}
		else
		{
#ifdef _MPI
			if (mpi_rank == 0)
			{
				MPI_File_seek(MPIfile,0, MPI_SEEK_SET);
				for (int k = 1; k <= 12; k++)
				{
					float v = (float)zz[k];
					MPI_Status s;
					MPI_File_write(MPIfile,&v, 1, MPI_FLOAT,&s);
				}
				MPI_File_seek(MPIfile,MPIoffset, MPI_SEEK_SET);
			}
#else
			// this branch requires no mpi_rank checking because it is compiled when _MPI is not deefined
			for (int k = 1; k <= 12; k++)
			{
				float v = (float)zz[k];
				if ((fwrite(&v, sizeof(float), 1, datei)) != 1) writeerror(completedname, 18);
			}
#endif
		}

		float *scol = NULL;
		double *dcol = NULL;
		if ((!bDoubleMode) && (!WriteASCIIGrids))
		{
			scol = new float[levels + 2];
		}
		dcol = new double[levels + 1];
		if (WriteASCIIGrids)
		{
			if (mpi_rank == 0)
			{
				for (i = 1; i <= nrows; i++)
				{
					if (stack[i])
					{
						for (int j = 1; j <= ncols; j++)
						{
							for (int k = 0; k <= levels; k++)
							{
								fprintf(datei, "%lg\n", stack[i][j][k]);
							}
						}
					}
					else
					{
						for (int j = 1; j <= ncols; j++)
						{
							for (int k = 0; k <= levels; k++)
							{
								fprintf(datei, "%lg\n", nodata);
							}
						}
					}
				}
			}
		}
		else
		{
#ifdef _MPI
			int sr = mpistartrow;
			int er = mpiendrow;
			bool bInRegion = false;
			if (mpi_rank == 0)
			{
				er = nrows;
			}
			if (mpi_rank == mpi_size - 1)
			{
				sr = 1;
				bInRegion = true;
			}

			double *dBuffer;
			float *fBuffer;
			unsigned long ulBufPos=0;

			// in order to reduce the number of write accesses, a buffer is reserved and filled and then written as a single big chunk instead of potentially thousands of small chunks
			if (bDoubleMode)
			{
				dBuffer = new double[(er - sr + 1)*ncols*(levels + 1)];
			}
			else
			{
				fBuffer = new float[(er - sr + 1)*ncols*(levels + 1)];
			}

			MPIoffset = 0;

			for (i = sr; i <= er; i++)
			{
				if (stack[i])
				{
					for (j = 1; j <= ncols; j++)
					{
						if ((MPIREG[i][j] != mpi_rank) && (MPIREG[i][j] != nodata)) // first cell of next region found
						{
							bInRegion = false;
						}
						else if ((MPIREG[i][j] == mpi_rank) || (bInRegion))
						{
							bInRegion = true;
							if (bDoubleMode)
							{
								if (MPIoffset == 0)
								{
									MPIoffset = (12 + ncols*(levels + 1)*(i - 1) + (j - 1)*(levels + 1))*sizeof(double);
								}
								for (int k = 0; k <= levels; k++)
								{
									dBuffer[ulBufPos] = stack[i][j][k];
									ulBufPos++;
								}
							}
							else
							{
								if (MPIoffset == 0)
								{
									MPIoffset = (12 + ncols*(levels + 1)*(i - 1) + (j - 1)*(levels + 1))*sizeof(float);
								}
								for (int k = 0; k <= levels; k++)
								{
									//	scol[k] = stack[i][j][k];
									fBuffer[ulBufPos] = stack[i][j][k];
									ulBufPos++;
								}
							}
						}
						else
						{
							// leave bInRegion at the same state
						}
					}
				}
				else
				{
					// this branch is used for the first rows of the highest rank (northernmost rows that may be empty -> we start with bInRegion==true!
					// and it is also used for the last rows that may also have no allocated memory
					if (bInRegion)
					{
						for (j = 1; j <= ncols; j++)
						{
							if (bDoubleMode)
							{
								if (MPIoffset == 0)
								{
									MPIoffset = (12 + ncols*(levels + 1)*(i - 1) + (j - 1)*(levels + 1))*sizeof(double);
								}
								for (int k = 0; k <= levels; k++)
								{
									dBuffer[ulBufPos] = nodata;
									ulBufPos++;
								}
							}
							else
							{
								if (MPIoffset == 0)
								{
									MPIoffset = (12 + ncols*(levels + 1)*(i - 1) + (j - 1)*(levels + 1))*sizeof(float);
								}
								for (int k = 0; k <= levels; k++)
								{
									fBuffer[ulBufPos] = nodata;
									ulBufPos++;
								}
							}
						}
					}
				}
			}

			if (bDoubleMode)
			{
				MPI_Status s;
				MPI_File_write_at(MPIfile,MPIoffset, &(dBuffer[0]), ulBufPos, MPI_DOUBLE,&s);
				MPI_File_set_size(MPIfile,(12 + ncols * nrows*(levels + 1)) * sizeof(double));
				delete[] dBuffer;
			}
			else
			{
				MPI_Status s;
				MPI_File_write_at(MPIfile,MPIoffset, &(fBuffer[0]), ulBufPos, MPI_FLOAT,&s);
				MPI_File_set_size(MPIfile,(12 + ncols * nrows*(levels + 1)) * sizeof(float));
				delete[] fBuffer;
			}

#else
			for (i = 1; i <= nrows; i++)
			{
				for (j = 1; j <= ncols; j++)
				{
					if (bDoubleMode)
					{
						if (((int)fwrite(&(stack[i][j][0]), sizeof(double), levels + 1, datei)) != levels + 1) writeerror(completedname, 18);
					}
					else
					{
						for (int k = 0; k <= levels; k++)
						{
							scol[k] = stack[i][j][k];
						}
						if (((int)fwrite(&(scol[0]), sizeof(float), levels + 1, datei)) != levels + 1) writeerror(completedname, 18);
					}
				}
			}
#endif
		}
		delete[] dcol;

		if ((!bDoubleMode) && (!WriteASCIIGrids)) delete[] scol;


		if (WriteASCIIGrids)
		{
			if (mpi_rank == 0)
			{
				if ((fclose(datei)) != 0)
				{
					(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
					throw(27);
				}
			}
		}
		else
		{
#ifdef _MPI
			MPI_File_close(&MPIfile);
#else
			// this branch requires no mpi_rank checking because it is compiled when _MPI is not deefined
			if ((fclose(datei)) != 0)
			{
				(*WSstream) << "FATAL: error closing " << completedname << "[end-entry]";
				throw(27);
			}
#endif
		}

#ifdef _NETCDF
	}
    if(WriteNetcdfGrids) {
		if (mpi_rank == 0)
		{

			netcdfOut->setGeometry(xll, yll, cellsize);

			char *completednamenetcdf = (char *)calloc(_stringsize_ + 1, sizeof(char));
			netcdfhlp::findFileNameFromGridName(completednamenetcdf, completedname, false);

			char *varname = (char *)calloc(_stringsize_ + 1, sizeof(char));
			netcdfhlp::findVarNameFromGridName(varname, completedname);

			char *varname_and_layer = (char *)calloc(_stringsize_ + 1, sizeof(char));

			double *latbuf, *lonbuf;

			if (gridinfo[1].has_latgrid_) {
				latbuf = new double[(size_t)ncols * nrows];
				size_t count = 0;
				for (j = 1; j <= ncols; j++) {
					for (i = 1; i <= nrows; i++) {
						latbuf[count] = (gridinfo[gridinfo[1].latgrid_index_].pointer[i]) ? gridinfo[gridinfo[1].latgrid_index_].pointer[i][j] : nodata;
						count++;
					}
				}
				netcdfOut->setLatGrid(latbuf);
			}

			if (gridinfo[1].has_longrid_) {
				lonbuf = new double[(size_t)ncols * nrows];
				size_t count = 0;
				for (j = 1; j <= ncols; j++) {
					for (i = 1; i <= nrows; i++) {
						lonbuf[count] = (gridinfo[gridinfo[1].longrid_index_].pointer[i]) ? gridinfo[gridinfo[1].longrid_index_].pointer[i][j] : nodata;
						count++;
					}
				}
				netcdfOut->setLonGrid(lonbuf);
			}

			// Since Wasim doesn't use row and column 0 of the array, we have to copy the data for netcdf
			if (bDoubleMode) {
				double *stackbuf = new double[(size_t)ncols * nrows];
				for (int k = 0; k <= levels; k++) {
					size_t count = 0;
					for (j = 1; j <= ncols; j++) {
						for (i = 1; i <= nrows; i++) {
							stackbuf[count] = (stack[i]) ? stack[i][j][k] : nodata;
							count++;
						}
					}
					//                strcpy(varname_and_layer, varname);
					//                sprintf(varname_and_layer,"%s%u",varname_and_layer,k);
					//				netcdfOut->writeTimestep(true, stackbuf, NULL, completednamenetcdf, varname_and_layer, nrows, ncols, nodata);
					netcdfOut->writeTimestep(true, stackbuf, NULL, completednamenetcdf, varname, nrows, ncols, nodata,cellsize);
				}
				delete[] stackbuf;
			}
			else {
				float *stackbuf = new float[(size_t)ncols * nrows];
				for (int k = 0; k <= levels; k++) {
					size_t count = 0;
					for (j = 1; j <= ncols; j++) {
						for (i = 1; i <= nrows; i++) {
							stackbuf[count] = (stack[i]) ? stack[i][j][k] : 0;
							count++;
						}
					}
					//                strcpy(varname_and_layer, varname);
					//                sprintf(varname_and_layer,"%s%u",varname_and_layer,k);
					//				netcdfOut->writeTimestep(false, NULL, stackbuf, completednamenetcdf, varname_and_layer, nrows, ncols, nodata);
					netcdfOut->writeTimestep(false, NULL, stackbuf, completednamenetcdf, varname, nrows, ncols, nodata,cellsize);
				}
				delete[] stackbuf;
			}

			free(completednamenetcdf);
			free(varname);
			free(varname_and_layer);
		}
    }

#else
    if(WriteNetcdfGrids) {
        (*WSstream) << "WARNING: Binary compiled without NetCDF support, skipping output."<<"[end-entry]";
    }
#endif

	delete[] completedname;
	SetSignalAbfangfunktion();
}


void modelobject::writestacktofile(int idx)
{

	//this->MPIGatherGrid(idx);
	//if (mpi_rank==0)
	//{
	//	writegridtofile(idx);
	//}

	if ((idx>gridmaxanzahl)||(idx<1)) return;
	int      writestack       = gridinfo[idx].writegrid;
	int      summenwritestack = gridinfo[idx].summenwritegrid;
	int      levels           = gridinfo[idx].levels;

	if ((writestack==0)&&(summenwritestack==0)) return;

	char     *name        = gridinfo[idx].name;
	double    ***stack     = gridinfo[idx].stackpointer;
	int      p;
	char     *dateiname,*tempstr;

	if ((writestack==3)&&((minf*minute + hf*stunde+df*tag+mf*monat+jf*jahr)<endtime))
	{  }
	else if (((writestack==1)||(writestack==2)||(writestack==3)) || ( ((writestack==4)||(writestack==5))&&((60 * stunde + (int)((minute + 60.0*hour_offset))) == 1440) ) )
	{
		tempstr=newstring(strlen(name)+7,1);
		/* Pfad vom Namen abschneiden, wenn ein dirzeichen auftaucht (\ oder / je nach System) */
		if (strrchr(name,dirzeichen)!=NULL)
		strcpy(tempstr,strrchr(name,dirzeichen)+1); else strcpy(tempstr,name);
		dateiname = newstring(strlen(name)+7,1);
		if ((writestack == 1) || (writestack == 3) || (writestack == 4))
		{
			strcpy(dateiname, name);
		}
		else if ((writestack==2)||(writestack==5))
		{
			if ( index(tempstr,".")==0) sprintf(tempstr + strlen(tempstr), "%s", ".");
			if ( index(tempstr,".") < 4 ) p = index(tempstr,"."); else p=4;

// New:		// the first call will always jump over this if statement since lastYear, lastMonth, lastDay and lastHour are always set to 0 during initialization of the gridinfo array
			// thus we make sure that the first call will always be compatible with the old file extension naming schema (only 3 digits for the hour)
			if ((jahr == gridinfo[idx].lastYear) && (monat == gridinfo[idx].lastMonth) && (tag == gridinfo[idx].lastDay) && (stunde == gridinfo[idx].lastHour) && (minute == gridinfo[idx].lastMinute))
			{
				gridinfo[idx].subIntervalCnt++; // increment sub interval counter and use it for creating the grid name
				sprintf(dateiname, "%.*s_%04ld%02ld%02ld.%03ld_%.03ld", p, tempstr, jahr, monat, tag, stunde, gridinfo[idx].subIntervalCnt);
			}
			else
			{
				// initializing date time counters when this is the first time the writegridtofiel function is called for a grid
				if (gridinfo[idx].lastYear == 0)
				{
					gridinfo[idx].lastYear = jahr;
					gridinfo[idx].lastMonth = monat;
					gridinfo[idx].lastDay = tag;
					gridinfo[idx].lastHour = stunde;
					gridinfo[idx].lastMinute = minute;
				}
				gridinfo[idx].subIntervalCnt = 0;
				sprintf(dateiname, "%.*s_%04ld%02ld%02ld.%03ld", p, tempstr, jahr, monat, tag, stunde);
			}
			gridinfo[idx].lastYear = jahr;
			gridinfo[idx].lastMonth = monat;
			gridinfo[idx].lastDay = tag;
			gridinfo[idx].lastHour = stunde;
			gridinfo[idx].lastMinute = minute;
			// old:			sprintf(dateiname, "%.*s%04ld%02ld%02ld.%03ld", p, tempstr, jahr, monat, tag, stunde);
			/* Pfad wird wieder angehaengt */
			if (strrchr(name,dirzeichen)!=NULL)
			{
				strcpy(tempstr,dateiname);
				p=(int)(strrchr(name,dirzeichen)-name)+1;
				sprintf(dateiname,"%.*s%s",p,name,tempstr);
			}
		}
		/* das im Speicher an der Adresse **grid gehaltene Array wird in ein Grid geschrieben */

#ifdef _MPI
//		MPI_Barrier(MPI_COMM_WORLD);
#endif
		// gather stack and write in rank 0 only
		double ***newstack;
		if (mpi_size > 1 && (WriteASCIIGrids || WriteNetcdfGrids))		{
			newstack = MPIGatherIntoNewStack(idx, stack);
		}	else	{
			newstack = stack;
		}
		if (mpi_rank == 0)
		{
			printstackfile(dateiname, newstack, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			if (!((writestack == 1) || (writestack == 3) || (writestack == 4)) && ((minf*minute + hf*stunde + df*tag + mf*monat + jf*jahr) >= endtime))
			{ // file extension is different from what is defined in the control file, so we will write that file in any case
				strcpy(dateiname, name);
				printstackfile(dateiname, newstack, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
			if (mpi_size > 1 && (WriteASCIIGrids || WriteNetcdfGrids))
			{
				free_drei_d_matrix(newstack, 1, nrows, 1, ncols, levels, nrows);
			}
		}
		else
		{
			printstackfile(dateiname, stack, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			if (!((writestack == 1) || (writestack == 3) || (writestack == 4)) && ((minf*minute + hf*stunde + df*tag + mf*monat + jf*jahr) >= endtime))
			{ // file extension is different from what is defined in the control file, so we will write that file in any case
				strcpy(dateiname, name);
				printstackfile(dateiname, stack, levels, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids);
			}
		}
		if (dateiname) free(dateiname); dateiname = NULL;
		if (tempstr) free(tempstr); tempstr=NULL;
	}
}



/* eine utility, die am Index i im Gridinfo-array Speicherplatz fuer
den Dateinamen des Grids sowie fuer die Bezeichnung reserviert und diese beiden
Strings in die Felder des gridinfo-arrays kopiert. Ausserdem wird,
writegrid gesplittet um die jeweiligen Codes fuer das staendige und fuer das
Summengrid zu ermitteln und zu speichern. Wenn ein Summengrid gespeichert werden soll,
dann wird dessen Name erzeugt und gespeichert (wie der normale Gridname,
jedoch mit einem .sum am Ende). Das Feld anzahl wird explizit auf 0 gesetzt (wenn auch nicht unbedingt noetig) */
void modelobject::initialize_grid(int i,
const char *dateiname,
double **grid,
double ***stack,
int Levels,
const char *Bezeichnung,
int writegrid,
int pnrl,
int pnrh)
{
	if (i==gridmaxanzahl)
	{ 
		(*WSstream) << "FATAL: no free grid in memory,(max. "<<gridmaxanzahl<<")!"<<"[end-entry]"; 
		throw(12); 
	}
	if (gridinfo[i].name==NULL) // for a number of grids, the name is allocated already externally (this array gridinfo is no class but only a data structure!)
	{
		gridinfo[i].name           = newstring(_stringsize_+1,1);
	}
	
	char* tempdateiname = new char[_stringsize_ + 1];
	strcpy(tempdateiname, dateiname); // this is done to avoid overlapping memory areas in case dateiname is gridinfo[i].name
	strcpy(gridinfo[i].name,tempdateiname);
	delete[] tempdateiname;

	// for all grids initialized within some sub model, the outname is probably not reserved. Only standard grids may have such a value
	// but the writegridtofile-function nowadays relies on this value, so we create it here
	if (!gridinfo[i].outname)
	{
		gridinfo[i].outname=newstring(_stringsize_+1,1);
		strcpy(gridinfo[i].outname,gridinfo[i].name);
                gridinfo[i].outnameinternally=1;
	}
	if (gridinfo[i].identification==NULL)
	{
		gridinfo[i].identification = newstring(101,1);                
	}

	char* tempBezeichnung = new char[101];
	strcpy(tempBezeichnung, Bezeichnung); // this is done to avoid overlapping memory areas in case dateiname is gridinfo[i].name
	strcpy(gridinfo[i].identification,tempBezeichnung);
	delete[] tempBezeichnung;

	gridinfo[i].pointer         = grid;
	gridinfo[i].stackpointer    = stack;
	gridinfo[i].levels          = Levels;
	gridinfo[i].writegrid       = (writegrid % 10)    / 1;     /* lowest significant digit of writecode is for actual values */
	gridinfo[i].summenwritegrid = (writegrid % 100 )  / 10;    /* second lowest digit is for mean or sum values */
	gridinfo[i].minwritegrid    = (writegrid % 1000)  / 100;   /* second most significant digit is the writecode for minimum or maximum grids */
	gridinfo[i].vegWritegrid    = (writegrid % 10000) / 1000;  /* largest significant digit is the writecode for vegetation period dependent grids */

	gridinfo[i].nrh = pnrh;
	gridinfo[i].nrl = pnrl;

	gridinfo[i].anzahl=0;
	gridinfo[i].minmaxflag=0;
	if (gridinfo[i].minwritegrid==1)
	{    gridinfo[i].minwritegrid=1; gridinfo[i].maxwritegrid=0;  }
	else
	{
		if (gridinfo[i].minwritegrid==2)
		{       gridinfo[i].minwritegrid=0; gridinfo[i].maxwritegrid=1;    }
		else
		{
			if (gridinfo[i].minwritegrid==3)
			{        gridinfo[i].minwritegrid=1; gridinfo[i].maxwritegrid=1;      }
			else
			{        gridinfo[i].minwritegrid=0; gridinfo[i].maxwritegrid=0;      }
		}
	}
	if (gridinfo[i].summenwritegrid!=0)
	{
		/* es soll ein Grid mit der Gesamtsumme des Modellaufes eingerichtet werden */
		gridinfo[i].summenpointer = matrix(gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
		gridinfo[i].summenname = newstring(_stringsize_+1,1);
		strcpy(gridinfo[i].summenname,gridinfo[i].name);
	}
	if (gridinfo[i].minwritegrid!=0)
	{
		gridinfo[i].minpointer = matrix(gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
		gridinfo[i].minname = newstring(_stringsize_,1);
		strcpy(gridinfo[i].minname,gridinfo[i].name);
	}
	if (gridinfo[i].maxwritegrid!=0)
	{
		gridinfo[i].maxpointer = matrix(gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
		gridinfo[i].maxname = newstring(_stringsize_,1);
		strcpy(gridinfo[i].maxname,gridinfo[i].name);
	}

	if (gridinfo[i].vegWritegrid !=0)
	{
		gridinfo[i].vegSumPointer = matrix(gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
		gridinfo[i].vegName = newstring(_stringsize_,1);
		strcpy(gridinfo[i].vegName,gridinfo[i].name);
	}
//	if (gridinfo[i].writegrid == 8)
//	{
//		this->writegridtofile(i);
//	}
}

void modelobject::delete_grid(int i)
{
	if (i<=gridmaxanzahl)
	{
		if (gridinfo[i].name) free(gridinfo[i].name);
		gridinfo[i].name=NULL;
		if (gridinfo[i].identification) free(gridinfo[i].identification);
		gridinfo[i].identification=NULL;
		gridinfo[i].writegrid=0;
		gridinfo[i].levels=0;
		gridinfo[i].pointer=NULL;
		gridinfo[i].stackpointer=NULL;
		gridinfo[i].anzahl=0;
		if (gridinfo[i].summenwritegrid!=0)
		{
			free_matrix(gridinfo[i].summenpointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
			gridinfo[i].summenpointer=NULL;
			if (gridinfo[i].summenname) free(gridinfo[i].summenname);
			gridinfo[i].summenname=NULL;
		}
		if (gridinfo[i].minwritegrid!=0)
		{
			free_matrix(gridinfo[i].minpointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
			gridinfo[i].minpointer=NULL;
			if (gridinfo[i].minname) free(gridinfo[i].minname);
			gridinfo[i].minname=NULL;
		}
		if (gridinfo[i].maxwritegrid!=0)
		{
			free_matrix(gridinfo[i].maxpointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
			gridinfo[i].maxpointer=NULL;
			if (gridinfo[i].maxname) free(gridinfo[i].maxname);
			gridinfo[i].maxname=NULL;
		}
		if (gridinfo[i].vegWritegrid !=0)
		{
			free_matrix(gridinfo[i].vegSumPointer, gridinfo[i].nrl, gridinfo[i].nrh, 1, ncols, nrows);
			gridinfo[i].vegSumPointer=NULL;
			if (gridinfo[i].vegName) free(gridinfo[i].vegName );
			gridinfo[i].vegName =NULL;
		}
                if (gridinfo[i].outnameinternally==1)
                {
                   free(gridinfo[i].outname);
                   gridinfo[i].outname=NULL;
	        }
        }
}

void modelobject::putnachrichtcode(int code)
{
	nachricht=code;
	return;
}

/* sucht in der geoeffneten Steuerdatei nach der Position des
als Parameter zu uebergebenden Abschnittes, der Dateizeiger
wird auf den ersten Datensatz des Abschnittes gesetzt */

bool modelobject::findeabschnitt(const char *Abschnittsname,bool bErrorOnNoFound)
{
	char *tempstr,*satz,*ttempstr,*tsatz;
	bool retval=false;
	tempstr=(char*)calloc(_stringsize_,sizeof(char)); ttempstr=tempstr;
	satz=(char*)calloc(_stringsize_,sizeof(char)); tsatz=satz;

	fseek(steuerdatei,0,0);
	sprintf(satz,"");
	int iReadCnt=1;
	while ((strcmp(satz,Abschnittsname)!=0)&&(!feof(steuerdatei))&&(!(iReadCnt<=0)))
	{
		iReadCnt=fscanf(steuerdatei,"%s",satz);
		fgets(tempstr,_stringsize_,steuerdatei);
	}
	if (feof(steuerdatei))  
	{
		retval=false;
		if (bErrorOnNoFound)
		{
			steuererror(Abschnittsname,"(EOF)");
		}
	}
	else   retval=true;
	tempstr=ttempstr; if (tempstr) free(tempstr); tempstr=NULL;
	satz=tsatz;       if (satz) free(satz);    satz=NULL;
	return retval;
}

/* liefert die indexnummer eines mit der Identification als Suchkriterium im
gridinfo-array gesuchten Grids, wenn dieses vorhanden ist, ansonsten 0 */

int modelobject::ist_vorhanden(const char *name,bool ignorepointer)
{
	int i;
	i=1;
	while (strcmp(gridinfo[i].identification,name)!=0)
	{
		i++;
		if ((i==gridmaxanzahl)||((gridinfo[i].pointer==NULL)&&(ignorepointer==false))||(gridinfo[i].name==NULL)) return 0;
	}
	return i;
}

int modelobject::getindex(const char *name)
{
	int  i=1;
	while (strcmp(gridinfo[i].identification,name)!=0)
	{
		i++;
		if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
		{
			(*WSstream) << "FATAL: \n---------------- !!! ERROR !!!  --------------------------\n"
			<<name<<" not found in memory,\n"
			"These data are needed within "<<instanzname<<" \n"
			"if they are interpolated data, please specify an interpolator for\n"
			<<name<<" in the control file or read it in as a standardfile,\n"
			"-----------------------------------------------------------\n"<<"[end-entry]";
			throw(33);
		}
	}
	return (int) i;
}

int modelobject::getindex_first_of_3(const char *name1,const char *name2,const char *name3)
/* Liefert den Index im gridinfo-Array zurueck, an dem der name1 auftaucht (im Feld identifîcation), wenn es ihn gibt;
ansonsten wird, wenn name2 != "", das erste Auftreten von name2, so es existiert, geliefert, wenn name2 auch
nicht im gridinfo[x].identification auftaucht, dann wird dasselbe mit name3 versucht. Gibt es keinen der drei
Namen im gridinfo-array, dann gibt es eine Fehlermeldung. */
{
	int i=1;
	while (strcmp(gridinfo[i].identification,name1)!=0)
	{
		i++;
		if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
		{
			i=1;
			while (strcmp(gridinfo[i].identification,name2)!=0)
			{
				i++;
				if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
				{
					i=1;
					while (strcmp(gridinfo[i].identification,name3)!=0)
					{
						i++;
						if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
						{
							(*WSstream) << "FATAL: \n---------------- !!! ERROR !!!  ------------------------------------\n"
							<< "Wether "<<name1<<", "<<name2<<" or "<<name3<<" \nwhere found in memory,\n"
							<< "One of these data sets has to be present within "<<instanzname<<" \n"
							<< "If the data are to be interpolated, please specify an interpolatior for\n"
							<< name1 << ", " << name1 << " or " <<name1<< " in the control file,\n"
							<< "or read it as a standard grid \n"
							<< "-----------------------------------------------------------------------\n"<<"[end-entry]";
							throw(33);
						}
					}
					break;
				}
			}
			break;
		}
	}
	return (int) i;
}

int modelobject::getindex_first_of_4(const char *name1,const char *name2,const char *name3,const char *name4)
/* Liefert den Index im gridinfo-Array zurueck, an dem der name1 auftaucht (im Feld identifîcation), wenn es ihn gibt;
ansonsten wird, wenn name2 != "", das erste Auftreten von name2, so es existiert, geliefert, wenn name2 auch
nicht im gridinfo[x].identification auftaucht, dann wird dasselbe mit name3 versucht. Gibt es keinen der drei
Namen im gridinfo-array, dann gibt es eine Fehlermeldung. */
{
	int i=1;
	while (strcmp(gridinfo[i].identification,name1)!=0)
	{
		i++;
		if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
		{
			i=1;
			while (strcmp(gridinfo[i].identification,name2)!=0)
			{
				i++;
				if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
				{
					i=1;
					while (strcmp(gridinfo[i].identification,name3)!=0)
					{
						i++;
						if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
						{
							i=1;
							while (strcmp(gridinfo[i].identification,name4)!=0)
							{
								i++;
								if ((i>=gridmaxanzahl)||(gridinfo[i].identification==NULL))
								{
									(*WSstream) << "FATAL: \n---------------- !!! ERROR !!!  ------------------------------------\n"
									<< "Wether "<<name1<<", "<<name2<<", "<<name3<<" or "<<name4<<" \nwhere found in memory,\n"
									<< "One of these data sets has to be present within "<<instanzname<<" \n"
									<< "If the data are to be interpolated, please specify an interpolatior for\n"
									<< name1<<", "<<name2<<", "<<name3<<" or "<<name4<<" in the control file,\n"
									<< "or read it as a standard grid \n"
									<< "-----------------------------------------------------------------------\n"<<"[end-entry]";
									throw(33);
								}
							}
							break;
						}
					}
					break;
				}
			}
			break;
		}
	}
	return (int) i;
}

double modelobject::interpolate(double x1,double y1,double x2,double y2,double x)
{
	if ((x2-x1)==0.0) return y1; else return (x-x1)/(x2-x1)*y2 + (x2-x)/(x2-x1)*y1;
}

void modelobject::set_end(long y,long m,long d,long h,long min)
{
	minutee = min;
	stundee = h;
	tage    = d;
	monate  = m;
	jahre   = y;
	endtime = jf*jahre+mf*monate+df*tage+hf*stundee+ minf * minutee;
}

modelobject::modelobject() :memory_handling()
{

}

modelobject::modelobject(wasim::mystream *Stream, char *steuername, int LogLevel, bool p_WriteASCIIGrids, bool p_WriteDoubleGrids, bool p_WriteNetcdfGrids, char *p_DefaultInputDirectory, char *p_DefaultOutputDirectory,
	long *pv_mpiregioncnt,std::vector<mpihalo> *pv_sendHalo,std::vector<mpihalo> *pv_recvHalo,int p_doSpinUp,std::vector<HRUMapElement*> *p_vHRU):
memory_handling(),
numsubb(0),listsubb(0),min_to_int_end(0),hour_offset(0),
satz(0),tsatz(0),tstr(0),ttstr(0),steuerdateiname(0),
steuerdatei(0),statanz(0),statmini(0),statmaxi(0),sx(0),sx2(0),stat_nodata(0),
ncols(0),nrows(0),xll(0),yll(0),yol(0),cellsize(0),nodata(0),fehlwert(0),
RelCellArea(0),CellSizeX(0),CellSizeY(0),CoordX(0),CoordY(0),
cellsizexindex(0),cellsizeyindex(0),relcellareaindex(0),coordxindex(0),coordyindex(0),
zone(0),zarray(0),rarray(0),FirstLastCol(0),zonmax(0),zonmaxwert(0),transdist(0),
regiomax(0),regiomaxwert(0),center_of_mass(0),maxregionenzahl(0),minute(0),stunde(0),
tag(0),monat(0),jahr(0),minutea(0),stundea(0),taga(0),monata(0),jahra(0),minutee(0),stundee(0),tage(0),monate(0),jahre(0),timestep(0),
jf(0),mf(0),df(0),hf(0),endtime(0),LUtable(0),SOILtable(0),MultiLayerLUTable(0),irrigationtable(0),
writegrid(0),gridinfo(0),result(0),resultnr(0),nachricht(0),iLogLevel(0),
MPIREG(0), MPIHLP(0), mpi_rank(0), mpi_size(1), WriteASCIIGrids(false), WriteDoubleGrids(false), WriteNetcdfGrids(false), HRU(0), HRUMaster(0),
xllcorner_ext(9999),yllcorner_ext(9999),cellsize_ext(9999),Mask(NULL),gmcut(127),gmf1(8.916787e-5),gmx1(0.0767528404),gmf2(8.916787e-5),gmx2(0.0767528404)
{
	WSstream=Stream;
	if ((instanzname = (char*)calloc(101,sizeof(char)))==NULL) speichererror(1);
	steuerdateiname  = newstring(_stringsize_+1,1);
	satz             = newstring(_stringsize_+1,1); tsatz=satz;
	tstr             = newstring(_stringsize_+1,1); ttstr=tstr;
	WriteASCIIGrids  = p_WriteASCIIGrids;
	WriteDoubleGrids = p_WriteDoubleGrids;
	WriteNetcdfGrids = p_WriteNetcdfGrids;
	DefaultInputDirectory = p_DefaultInputDirectory;
	DefaultOutputDirectory = p_DefaultOutputDirectory;
	this->v_mpiregioncnt = pv_mpiregioncnt;
	this->sendHalo = pv_sendHalo;
	this->recvHalo = pv_recvHalo;

	this->doSpinUp = p_doSpinUp;
	this->vHRU     = p_vHRU;

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
#endif

	strcpy(steuerdateiname,steuername);

	/* initialisieren des Monatsfeldes */
	mon[1]  = 31;
	mon[2]  = 28;
	mon[3]  = 31;
	mon[4]  = 30;
	mon[5]  = 31;
	mon[6]  = 30;
	mon[7]  = 31;
	mon[8]  = 31;
	mon[9]  = 30;
	mon[10] = 31;
	mon[11] = 30;
	mon[12] = 31;

	jf = 535680; /* factors for comparing dates assembled from year, month, day, hour and minute*/
	mf = 44640;
	df = 1440;
	hf = 60;
	minf = 1; // factor per minute, if used at all

#ifdef _NETCDF
    netcdfOut = new netcdfWriter(WSstream);
    netcdfIn = new netcdfReader(WSstream);
    //gridinfo[1].has_latgrid_ = false;
    //gridinfo[1].has_longrid_ = false;
#endif

}


void modelobject::ClearMemory()
{
	(*WSstream) << "DEBUG: deleting module "<<instanzname<<"[end-entry]";
	if (instanzname    ) free(instanzname    );  instanzname=NULL;
	if (steuerdateiname) free(steuerdateiname);  steuerdateiname=NULL;
	satz=tsatz;  if (satz) free(satz);  satz=NULL;
	tstr=ttstr;  if (tstr) free(tstr);  tstr=NULL;
    
#ifdef _NETCDF
    delete netcdfOut;
    delete netcdfIn;
#endif

}
/* =================== Ende der modelobject-routinen =====================  */

// Filebuffer-Routinen

FileBuffer::FileBuffer() :
	//m_file(0),
	m_filename(""),
	stringbuffer1(0),
	stringbuffer2(0),
	m_buffer(0),
	m_memsize(0),
	m_bufpos(0),
	m_bufmax(0),
	m_valcnt(0),
	m_filepos(0),
	m_mode(0),
	m_status(0), // status 0 means: still uninitialized
	m_useMinutes(false)
{
}

FileBuffer::FileBuffer(std::string filename, int mode, long memorysize) :
	//m_file(0),
	m_filename(filename),
	stringbuffer1(0),
	stringbuffer2(0),
	m_buffer(0),
	m_memsize(memorysize),
	m_bufpos(0),
	m_bufmax(0),
	m_valcnt(0),
	m_filepos(0),
	m_mode(mode),
	m_status(0),
	m_useMinutes(false)
{
	int i=this->initialize(m_filename,m_mode,m_memsize);
}

FileBuffer::~FileBuffer()
{
	if (stringbuffer1) delete[] stringbuffer1;
	if (stringbuffer1) delete[] stringbuffer2;
	if (m_buffer)
	{
		for (int i=0;i<m_bufmax;i++)
		{
			delete[] m_buffer[i].wert;
		}
		delete[] m_buffer;
	}
}

void FileBuffer::ClearMemory()
{
	if (stringbuffer1) 
	{
		delete[] stringbuffer1;
		stringbuffer1=NULL;
	}
	if (stringbuffer2)
	{
		delete[] stringbuffer2;
		stringbuffer2=NULL;
	}
	if (m_buffer)
	{
		for (int i=0;i<m_bufmax;i++)
		{
			delete[] m_buffer[i].wert;
		}
		delete[] m_buffer;
		m_buffer=NULL;
	}
}

// clear all buffer elements, reset filepos, reset bufpos,clear filename
void FileBuffer::reset()
{
}

int FileBuffer::initialize(std::string filename,int mode,long memorysize)
{
	int Header_Rows;
	if (m_status==0)
	{
		std::stringstream sstr;
		int i,stringbuffersize;
		char c,nextchar;
		double value;
		unsigned int y,m,d,h,min;
		m_filename=filename,
		m_memsize=memorysize,

		stringbuffersize=4096;
		stringbuffer1=new char[stringbuffersize];
		stringbuffer2=new char[stringbuffersize];

		// check, if file can be opened
		m_file.open(m_filename.c_str(),std::ios_base::in|std::ios::binary);
		if (!m_file.good())
		{
			m_status=3;
			return m_status;
		}
		
		// read away the header: header rows must not start with a number,
		// but the first valid row must start with a number!
		Header_Rows=0;
		nextchar = m_file.peek();
		while ((nextchar<0x30 || nextchar>0x39) && (m_file.good())&& (!m_file.eof()))
		{
		
//			do 
//			{
//				m_file.get(c);
//			}
//			while ((c!='\n') && (m_file.good())&& (!m_file.eof()));
			m_file.getline(stringbuffer1, stringbuffersize);
			Header_Rows++;
			std::string row(stringbuffer1);
			int idx = 0;
			int len = 0;
			if ((row.find("YY") != std::string::npos) && (row.find("MM") != std::string::npos) && (row.find("DD") != std::string::npos) && (row.find("HH") != std::string::npos))
			{
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
					m_useMinutes = true;
				}
			}

			nextchar = m_file.peek();
		}

		if (!m_file.good()) {
			m_status=FB_STATE_FILENOTGOOD;	
			return m_status;
		}

		m_filepos=m_file.tellg();

		// now analyzing the first valid row for the number of values
		if (m_file.good() && (!m_file.eof()))
		{
			m_valcnt = 0;
			// read date
			if (!m_useMinutes)
			{
				m_file >> y >> m >> d >> h;
			}
			else
			{
				m_file >> y >> m >> d >> h >> min;
			}

			// read the values into a stringstream
			do
			{
				m_file.get(c);
				sstr << c;
			}
			while (c!='\n');

			if (sstr.good())
			{   
				do
				{
					sstr >> value;
					if (sstr.good()) m_valcnt++;
				}   
				while (sstr.good());
			} 
			else
			{
				throw 201;
			}
			m_file.seekg(0,ios_base::end);
			m_filesize = m_file.tellg();
			m_file.close();
		}
		else
		{
			throw 202;
		}

		// now allocating the buffer
		m_bufmax=m_memsize / (m_valcnt*sizeof(double)+4*sizeof(int));
		m_buffer=new buffertype[m_bufmax+1];
		for (int i=0;i<m_bufmax;i++)
		{
			m_buffer[i].wert=new double[m_valcnt];
		}
		m_bufpos=m_bufmax+1; // in order to be ready to read the file when the buffer is accessed for the first time
		m_status=1; // initialized
		ReadBufferFromFile();

	}
	return(Header_Rows);
}

// for free access to the buffer, if required
void FileBuffer::SetBufPos(long pos)
{
	if ((pos>=-1)&&(pos<=this->m_bufmax))
	{
		this->m_bufpos = pos;
	}
}

// for free access to the buffer, if required; 
// relpos is relative to the actual position, can be positive or negative
void FileBuffer::SetBufPosRel(long relpos)
{
	if (((m_bufpos+relpos)>=-1)&&((m_bufpos+relpos)<=this->m_bufmax))
	{
		this->m_bufpos += relpos;
	}
}

void FileBuffer::ReadBufferFromFile()
{
	if (m_status==1)
	{
		m_bufpos=0;
		m_file.open(m_filename.c_str(),std::ios_base::in|std::ios::binary);
		m_file.seekg(m_filepos);
//		char *record = new char[2000];
		while (m_file.good() && !m_file.eof() && m_bufpos<m_bufmax)
		{

			long actpos = m_file.tellg();
			if (m_filesize-actpos<10)
			{
				m_file.setstate(std::ios_base::eofbit);
				break; // end of file expected
			}
			

//			m_file.getline(record, 2000);
//			stringstream ss(record);
//			ss >> m_buffer[m_bufpos].j >> m_buffer[m_bufpos].m >> m_buffer[m_bufpos].t >> m_buffer[m_bufpos].s;

			if (!m_useMinutes)
			{
				m_file >> m_buffer[m_bufpos].j >> m_buffer[m_bufpos].m >> m_buffer[m_bufpos].t >> m_buffer[m_bufpos].s;
				m_buffer[m_bufpos].min = 0;
			}
			else
			{
				m_file >> m_buffer[m_bufpos].j >> m_buffer[m_bufpos].m >> m_buffer[m_bufpos].t >> m_buffer[m_bufpos].s >> m_buffer[m_bufpos].min;
			}

			if (!m_file.good()) {
				m_status=FB_STATE_FILENOTGOOD;	
				return;
			}
			
			for (int i=0;i<m_valcnt;i++)
			{
//				ss >> m_buffer[m_bufpos].wert[i];
				m_file >> m_buffer[m_bufpos].wert[i];
				if (!m_file.good()) {
					m_status=FB_STATE_FILENOTGOOD;	
					return;
				}
			}

			if (!m_file.eof())
			{
				m_bufpos++;
			}
			else
			{
				break;
			}
		}
		if (m_file.eof())
		{
			m_status=2; // file end reached
		}
		else 
		{
			m_status=1;
		}
//		delete[] record;
		m_bufmax=m_bufpos; 
		m_bufpos=0;
		m_filepos=m_file.tellg();
		m_file.close();
	}
}

void FileBuffer::WriteBufferToFile()
{
}

void FileBuffer::NextRecord()
{
	if (m_bufpos<m_bufmax-1)
	{
		m_bufpos++;
	}
	else
	{
		if (m_status==1)
		{
			ReadBufferFromFile();
		}
		else 
		{
			m_status=4; // last record already read or file read error
		}
	}
}

void FileBuffer::GetActdate(int &year,int &month,int &day,int &hour,int &minute)
{
	year  = m_buffer[m_bufpos].j;
	month = m_buffer[m_bufpos].m;
	day   = m_buffer[m_bufpos].t;
	hour = m_buffer[m_bufpos].s;
	minute = m_buffer[m_bufpos].min;
}

void FileBuffer::GetActValues(double *value)
{
	value=m_buffer[m_bufpos].wert;
}

double FileBuffer::GetActValue(int i)
{
	if ((i>=0)&&(i<m_valcnt))
	return m_buffer[m_bufpos].wert[i];
	else
	return -9999;
}

size_t StringUtils::StringSplit (const std::string &text, const std::string &separators, std::vector< std::string > &words, bool bAddEmptyStrings )
{
	size_t InsideString;
	size_t n = text.length();
	size_t start = text.find_first_not_of(separators);
	size_t stop;

	// if the string constists of separators only, 
	// it should be split into as many substrings 
	// as the separator characters indicate

	// one new if-statement at 07.07.2002, SUL
	if (start==std::string::npos)
	{
		if (bAddEmptyStrings)
		start = text.find_first_of(separators);
	}  

	while (start < n)
	{
		stop=start;

		if ( (text[start]=='\"') && (separators.find("\"")==std::string::npos) )
		{
			do
			{
				InsideString=0;
				stop = text.find_first_of("\"",stop+1);
				if ( (stop>0)&&(stop<n)&&(text[stop-1]=='\\'))
				{
					if ( !( (stop>1)&&(text[stop-2]=='\\') ) )
					InsideString=1;
				}
			} while (InsideString);
			stop++; // the \" itself belongs still to the string
		}
		else
		{
			stop = text.find_first_of(separators, start);
		}
		if (stop > n) stop = n;
		words.push_back(text.substr(start, stop - start));
		start=stop+1;
		start = text.find_first_not_of(separators, stop+1); 
		// would be used if multiple separator characters 
		// should be treated as a single one
	}
	// new if-statement added at 07.07.2002, SUL
	// adding an empty string if the last character 
	// of the input string was a separator character
	if (text.find_last_of(separators) == n-1)
	{
		if (bAddEmptyStrings)
		words.push_back(std::string(""));  
	}
	return words.size();
}
