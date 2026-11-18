/*   $Revision: 3 $:  $
	$Date: 16.05.10 0:09 $  */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include "cspline.h"

#include "modconst.h"
using namespace std;

// Constructor
TSplines::TSplines(double **tCoordXe,double **tCoordYe,double **tCoordXa,double **tCoordYa,long tnrowsein,long tncolsein,double txllein,double tyllein,double tcseinx,double tcseiny,double tnodataein,
long tnrowsaus,long tncolsaus,double txllaus,double tyllaus,double tcsaus,double tnodataaus,double tbeta,double tfinefact,long tequalarea,double **tzn):
memory_handling(),
yNRAbyNCE(0),ycNRAbyNCE(0),xcNRAbyNCE(0),tempgrid(0),
CoordXe(0),CoordYe(0),CoordXa(0),CoordYa(0),zn(0),
x(0),y2(0),xcol(0),y2col(0),finefact(0),beta(0),
y(0),ycol(0),
xtempmin(0),xtempmax(0),ytempmin(0),ytempmax(0),xlltemp(0),ylltemp(0),xlltemprot(0),ylltemprot(0),
nrtempmin(0),nrtempmax(0),nctempmin(0),nctempmax(0),nrtemp(0),nctemp(0),
equalarea(0),nrowsein(0),ncolsein(0),nrowsaus(0),ncolsaus(0),
xllein(0),yllein(0),cseinx(0),cseiny(0),nodataein(0),xllaus(0),yllaus(0),csaus(0),nodataaus(0),
i(0),j(0),k(0),m(0),klo(0),khi(0),anf(0),end(0),r(0),r1(0),r2(0),c(0),rowindex(0),colindex(0),
rf(0),tx1(0),ty1(0),shy(0),shx(0),
CYa(0),CXe(0),CXa(0),cosbeta(0),sinbeta(0),yoffset(0),xlleinrot(0),ylleinrot(0),
dx(0),dy(0),entf(0),alpha(0),gamma(0),tst(0),LookUpTable(0)
{
	int i,j;
	CoordXe    = tCoordXe;
	CoordYe    = tCoordYe;
	CoordXa    = tCoordXa;
	CoordYa    = tCoordYa;
	nrowsein   = tnrowsein;
	ncolsein   = tncolsein;   
	xllein     = txllein;     
	yllein     = tyllein;     
	cseinx     = tcseinx;      
	cseiny     = tcseiny;      
	nodataein  = tnodataein;  
	nrowsaus   = tnrowsaus;   
	ncolsaus   = tncolsaus;   
	xllaus     = txllaus;     
	yllaus     = tyllaus;     
	csaus      = tcsaus;      
	nodataaus  = tnodataaus;  
	beta       = tbeta; 
	finefact   = tfinefact;
	equalarea  = tequalarea;
	zn         = tzn;

	/* getting the coordinates range for the rotated input grid (considering skewness) */
	xtempmin=xtempmax=xllein;
	ytempmin=ytempmax=yllein;
	for (i=1;i<=nrowsein;i++)
	{
		if (xtempmin>CoordXe[i][1]) 
		xtempmin=CoordXe[i][1];
		if (xtempmax<CoordXe[i][ncolsein]) 
		xtempmax=CoordXe[i][ncolsein];
	}
	for (i=1;i<=ncolsein;i++)
	{
		if (ytempmin>CoordYe[nrowsein][i]) 
		ytempmin=CoordYe[nrowsein][i];
		if (ytempmax<CoordYe[1][i]) 
		ytempmax=CoordYe[1][i];
	}
	csaus/=finefact;
	nrtemp=(long)((ytempmax-ytempmin)/csaus+0.49999)+1;
	nctemp=(long)((xtempmax-xtempmin)/csaus+0.49999)+1;
	xlltemp=xtempmin;
	ylltemp=ytempmin;

	GetValidSplineArea();

	tempgrid   = matrix(1,nrtempmax-nrtempmin+3,1,nctempmax-nctempmin+3);
	LookUpTable=new unsigned short * [nrtempmax-nrtempmin+3];
	for (i=1;i<nrtempmax-nrtempmin+3;i++)
	{
		LookUpTable[i] = new unsigned short [nctempmax-nctempmin+3];
		for (j=1;j<nctempmax-nctempmin+3;j++)
		LookUpTable[i][j]=0;
	}
	LookUpTableInit();

	//  tst=fopen("lookup.asc","wt");
	//  fprintf(tst,"ncols       %d\n"
	//              "nrows       %d\n"
	//              "xll         %g\n"
	//              "yll         %g\n"
	//              "cs          %g\n"
	//              "nodatavalue %g\n",nctempmax-nctempmin+2,
	//              nrtempmax-nrtempmin+2,
	//              xlltemp+csaus*(nctempmin-1),
	//              ylltemp+csaus*(nrtemp-nrtempmax),csaus,nodataaus);
	//  for (i=1;i<nrtempmax-nrtempmin+3;i++)
	//  {
	//    for (j=1;j<nctempmax-nctempmin+3;j++)
	//      fprintf(tst,"%d\t",LookUpTable[i][j]);
	//    fprintf(tst,"\n");
	//  }
	//  fclose(tst);

	yNRAbyNCE  = matrix(1,nrtemp+2,1,ncolsein+2);
	ycNRAbyNCE = matrix(1,nrtemp+2,1,ncolsein+2);
	xcNRAbyNCE = matrix(1,nrtemp+2,1,ncolsein+2);
	xcol  = new double[ncolsein+2];
	ycol  = new double[ncolsein+2];
	y2col = new double[ncolsein+2];

	x=new double[nrowsein+2];
	y=new double[nrowsein+2];
	y2=new double[nrowsein+2];
	cosbeta   = cos(beta);
	sinbeta   = sin(beta);

	/* berechnen der rotierten Koordinaten der linken unteren Ecke des Stationengrids
	(das Ausgabegrid behält seine Eckkordinaten, wiel die Rotation um diese Ecke erfolgt) */
	xlleinrot=xllein;
	ylleinrot=yllein;
	RotateCoordinates(&xlleinrot,&ylleinrot,beta);
}

//Destructor
TSplines::~TSplines()
{
	int i;
	if (xcol) delete xcol;
	if (ycol) delete ycol;
	if (y2col) delete y2col;
	if (x) delete x;
	if (y) delete y;
	if (y2) delete y2;

	//  free_matrix(yNRAbyNCE,1,nrowsaus+2,1,ncolsein+2);
	//  free_matrix(ycNRAbyNCE,1,nrowsaus+2,1,ncolsein+2);
	//  free_matrix(xcNRAbyNCE,1,nrowsaus+2,1,ncolsein+2);

	free_matrix(yNRAbyNCE,1,nrtemp+2,1,ncolsein+2);
	free_matrix(ycNRAbyNCE,1,nrtemp+2,1,ncolsein+2);
	free_matrix(xcNRAbyNCE,1,nrtemp+2,1,ncolsein+2);
	free_matrix(tempgrid,1,nrtempmax-nrtempmin+3,1,nctempmax-nctempmin+3);
	for (i=1;i<nrtempmax-nrtempmin+3;i++)
	if (LookUpTable[i]) delete LookUpTable[i];
	if (LookUpTable) delete LookUpTable;
}

void TSplines::LookUpTableInit()
{
	int i,j;
	if (zn)
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if ((zn[i][j]>nodataaus+0.01)&&(CoordXa[i][j]>nodataaus+0.01)&&(CoordYa[i][j]>nodataaus+0.01))
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					if ((r>=nrtempmin)&&(c>=nctempmin)&&(r<=nrtempmax)&&(c<=nctempmax))
					{
						LookUpTable[r-nrtempmin+1][c-nctempmin+1]=1;
					}
				}
			}
		}
	}
	else
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if ((CoordXa[i][j]>nodataaus+0.01)&&(CoordYa[i][j]>nodataaus+0.01))
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					if ((r>=nrtempmin)&&(c>=nctempmin)&&(r<=nrtempmax)&&(c<=nctempmax))
					{
						LookUpTable[r-nrtempmin+1][c-nctempmin+1]=1;
					}
				}
			}
		}
	}
}


void TSplines::GetValidSplineArea()
{
	int i,j;
	nrtempmin=nrtemp;
	nctempmin=nctemp;
	nrtempmax=1;
	nctempmax=1;
	if (zn)
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if ((zn[i][j]>nodataaus+0.01)&&(CoordXa[i][j]>nodataaus+0.01)&&(CoordYa[i][j]>nodataaus+0.01))
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					if ((r>0)&&(c>0)&&(r<=nrtemp)&&(c<=nctemp))
					{
						if (r<nrtempmin)
						nrtempmin=r;
						if (r>nrtempmax)
						nrtempmax=r;
						if (c<nctempmin)
						nctempmin=c;
						if (c>nctempmax)
						nctempmax=c;
					}
				}
			}
		}
	}
	else
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if ((CoordXa[i][j]>nodataaus+0.01)&&(CoordYa[i][j]>nodataaus+0.01))
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					//          r=nrtemp-(CoordYa[i][j]-ylltemp)/csaus;
					//          c=(CoordXa[i][j]-xlltemp)/csaus+1;
					if ((r>0)&&(c>0)&&(r<=nrtemp)&&(c<=nctemp))
					{
						if (r<nrtempmin)
						nrtempmin=r;
						if (r>nrtempmax)
						nrtempmax=r;
						if (c<nctempmin)
						nctempmin=c;
						if (c>nctempmax)
						nctempmax=c;
					}
				}
			}
		}
	}
	if (nctempmin>1)      nctempmin--;
	if (nctempmax<nctemp) nctempmax++;
	if (nrtempmin>1)      nrtempmin--;
	if (nrtempmax<nrtemp) nrtempmax++;
}

void TSplines::ReRotateResults(double **za)
{
	int i,j;
	/* rotating the coordinates of the intermediate result into the real result grid */
	if (zn)
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if (zn[i][j]>nodataaus+0.01)
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					//          r=nrtemp-(CoordYa[i][j]-ylltemp)/csaus;
					//          c=(CoordXa[i][j]-xlltemp)/csaus+1;
					if ((r>=nrtempmin)&&(c>=nctempmin)&&(r<=nrtempmax)&&(c<=nctempmax))
					za[i][j]=tempgrid[r-nrtempmin+1][c-nctempmin+1];
					else
					za[i][j]=nodataaus;
				}
				else
				za[i][j]=nodataaus;
			}
		}
	}
	else
	{
		for (i=1;i<=nrowsaus;i++)
		{
			for (j=1;j<=ncolsaus;j++)
			{
				if ((CoordXa[i][j]>nodataaus+0.01)&&(CoordYa[i][j]>nodataaus+0.01))
				{
					r=nrtemp-(int)((CoordYa[i][j]-ylltemp)/csaus+0.4999);
					c=(int)((CoordXa[i][j]-xlltemp)/csaus)+1;
					//          r=nrtemp-(CoordYa[i][j]-ylltemp)/csaus;
					//          c=(CoordXa[i][j]-xlltemp)/csaus+1;
					if ((r>=nrtempmin)&&(c>=nctempmin)&&(r<=nrtempmax)&&(c<=nctempmax))
					//          if ((r>0)&&(c>0)&&(r<=nrtemp)&&(c<=nctemp))
					za[i][j]=tempgrid[r-nrtempmin+1][c-nctempmin+1];
					else
					za[i][j]=nodataaus;
				}
				else
				za[i][j]=nodataaus;
			}
		}
	}
}


void TSplines::RotateCoordinates(double *x,double *y,double beta)
{
	double dx,dy,entf,alpha,gamma;
	dx = (*x-xllaus);
	dy = (*y-yllaus);
	entf = sqrt(dx*dx+dy*dy);
	/* actual angle of the points vector wiht respect to xll and yll */
	if (dx==0.0)
	alpha = dy > 0.0 ? Pi/2.0 : -Pi/2.0;
	else alpha = atan2(dy,dx);
	/* Winkel nach rechts drehen um Aniotropiewinkel */
	gamma = alpha - beta;
	/* neue Koordinaten */
	*x = xllaus+cos(gamma)*entf;
	*y = yllaus+sin(gamma)*entf;
}

// memberfunction ConstructCubicSpline
// to construct the parameters for a 1-D-Cubic Spline (natural cubic spline if yp1 and ypn are >=1e30)
// version 1, 27.7.1999, Joerg Schulla
// parameters: 
// double *x         :  input x values as a vector [0..n+1]
// double *y         :  input y values as a vector [0..n+1]
// long anf          :  first index of tabulated values
// long end          :  last index of tabulated values
// double yp1 and ypn:  initial boundary conditions for the second derivatives, if >=1e30
//                     then the boundary conditions are set to 0 (natural spline)
// double *y2        :  result of the constructed spline [0..n+1](second derivatives at the given coordinates in *x) 
void TSplines::ConstructCubicSpline(double *x,double *y,long anf,long end,double yp1,double ypn,double *y2)
{
	long i,k;
	double p,qn,sig,un;
	double *u;
	u=new double[end+1];

	/* lower boundary condition */
	if (yp1>0.99e30)
	{  y2[anf]=u[anf]=0.0; }
	else
	{
		y2[anf]=-0.5;
		u[anf]=(double)((3.0/(x[anf+1]-x[anf]))*((y[anf+1]-y[anf])/(x[anf+1]-x[anf])-yp1));
	}

	/* decomposition loop forthe tridiagonal matrix, y2 and u are used for temp. storage */
	for (i=anf+1;i<end;i++)
	{
		sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p     = (double)(sig*y2[i-1]+2.0);
		y2[i] = (double)((sig-1.0)/p);
		u[i]  = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]  = (double)((6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p);
	}

	/* upper boundary condition */
	if (ypn>0.99e30)
	{  qn=un=0.0; }
	else
	{
		qn=-0.5;
		un=(double)((3.0/(x[end]-x[end-1]))*(ypn-(y[end]-y[end-1])/(x[end]-x[end-1])));
	}
	y2[end]=(double)((un-qn*u[end-1])/(qn*y2[end-1]+1.0));

	/* decomposition loop */
	for (k=end-1;k>=anf;k--)
	y2[k]=y2[k]*y2[k+1]+u[k];

	if (u) delete u;
}

// function EvalCubicSpline
// to evaluate a 1-D-Spline (get the real interpolation results)
// version 1, 28.7.1999, Joerg Schulla
// parameters: 
// double *x         :  input x values as a vector [0..n_ein+1], arranged in ascending order!!!
// double *y         :  input y values as a vector [0..n_ein+1]  (must match the same order as the xa's!!!)
// double *y2        :  input y2 (2nd derivatives) as a vector [0..n_ein+1] as result of ConstructCubicSpline
// long n            :  number of tabulated input values
// double xa         :  output-grid x value 
// double ya         :  result of the spline evaluation (interpolation results for the coordinates given in xa parameter) 
// long equalarea    :  if set to 1, then the last value in a row in repeated in ordwer to allow interpolation within single gridcells
void TSplines::EvalCubicSpline(double *x,double *y,long *klo,long *khi,double *y2,long n,double xa,double *ya,double nodataein,double nodataaus)
{
	long changeflag;
	double a,b,h;

	/* looking for klo and khi starting from the given values */
	while ((x[*klo]>xa)&&(*klo>0))
	{ /* xa is outside the lower boundary of the given interval */
		(*klo) -=1;
		if ((*klo)<=0) break;
	}
	while ((x[*klo]<=xa)&&(*klo<n))
	{ 
		if (x[(*klo)+1]>xa) break;
		(*klo)+=1;
		if ((*klo)>n) break;
	}
	*khi = (*klo)+1;

	while ((x[*khi]<xa)&&(*khi<=n))
	{ /* xa is outside the lower boundary of the given interval */
		(*khi)++;
		if ((*khi)>n) break;
	}
	/* in order to be able to interpolate in a grid with valid values for entire grid cells (not just station data which are 
	valid at a grid corner or at the center!) it is necessary for grid interpolation, to set the values of the outputgrid 
	at least to the old input grid value for all grid cells lying within the area of the input grid cell. Thus, even 
	a single grid cell may be interpolated (e.g. if one interpolates a single cell to a cellsize half as large as the original
	cell size one will have only one point, but this point has 4 corners which will be substituted by the originalvalue of 
	the only valid point - the result will be, that all cells have the same output value, but the same area will be covered!
	-> give the nodata values at the upper boundary a temporary value
*/
	changeflag=0;
	if (((y[*khi]<nodataein+0.01)||(y2[*khi]<nodataein+.01))&&(equalarea>0))
	{
		y[*khi]=y[*klo];
		y2[*khi]=y2[*klo];
		changeflag=1;
	}

	if (((*klo)<1)||((*khi)>n+1)||((*klo)>n)||((*khi)<2))
	{
		*ya=nodataaus;
	}
	else
	{
		if ((y[*khi]<nodataein+.01)||(y[*klo]<nodataein+.01)||(y2[*khi]<nodataein+.01)||(y2[*klo]<nodataein+.01))
		{
			*ya=nodataaus;
		}
		else
		{
			h=x[*khi]-x[*klo];
			if (h==0.0)
			{
				cout << "The two input coordinates "<<x[*klo]<<" (index "<<*klo<<") \n"
				<< "and "<<x[*khi]<<" (index "<<*khi<<") for the 1-D-Splineevaluation are the same!"<<endl;
				throw 82;
			}
			a=(x[*khi]-xa)/h;
			b=(xa-x[*klo])/h;
			*ya=(double)(a*y[*klo]+b*y[*khi]+
			((a*a*a-a)*y2[*klo]+(b*b*b-b)*y2[*khi])*(h*h)/(double)6.0);
		}
	}
	if (changeflag)
	{
		y[*khi]=nodataein;
		y2[*khi]=nodataein;
	}
	return;
}

void TSplines::BicubicColumnSpline(double **ze,double **y2a)
{
	long i,k,anf,end;

	/* do the spline construction for all line segments within a column 
	single pixels are handled in the way that their 2nd derivativs are set to 0 */
	/* 1. step: construction of a 1-D-cubic spline for all input columns */
	//  cout << "contructing the bicubic spline";
	//  printf("column splines...");
	for (i=1;i<=ncolsein;i++)
	{
		//    printf(".");
		anf=1; end=0;
		while (anf<nrowsein) 
		{
			anf=end+1;
			/* suchen Anfang ab letzter end-Position */
			while ((anf<=nrowsein)&&(ze[anf][i]<nodataein+.01)) 
			{
				y2a[anf][i]=nodataaus;
				anf++;
			}
			/* suchen Ende ab aktueller anf-Position */
			end=anf;
			while ((end<nrowsein)&&(ze[end+1][i]>nodataein+.01)) 
			end++;

			for (k=anf;k<=end && k<=nrowsein;k++)
			{
				//        x[k]=yllein+(nrowsein-k)*csein;
				x[k]=CoordYe[k][i];
				y[k]=ze[k][i];
			}
			/* the last edge of a grid sqare sequence is taken as the same value as the forelast, 
		in order to be able to interpolate also within single grid cells (and not to have to write nodata values there) */
			if (anf<=nrowsein) /* at least 2 zells are needed */
			{
				if (equalarea)
				{
					//          x[end+1]=yllein+(nrowsein-k)*csein;
					x[end+1]=x[end]-cseiny;
					y[end+1]=y[end];
				}
				ConstructCubicSpline(x,y,anf,end+equalarea,(double)1e30,(double)1e30,y2);
				for (k=anf;k<=end+equalarea;k++)
				y2a[k][i]=y2[k];
			}
			else
			{
				if (anf<=nrowsein)
				{
					if (ze[anf][i]>nodataein+.01)
					y2a[anf][i]=y2[k];
					else
					y2a[anf][i]=nodataaus;
				}
			}
		}
	}
	//  printf("\n");

}

void TSplines::BicubicRowSpline(double **ze,double **y2a,double **za)
{
	/* Es wird nun zunaechst ein Array belegt, dass interpolierte Werte fuer alle Spalten
	der Inputdaten enthält, aber darin jeweils fuer alle Reihen des Ausgabegrids.
	Diese Operation muss mit rotierten Koordinaten durchgefuehrt werden,
	es reicht, wenn sie einmal fuer das gesamte Grid ausgefuehrt wird! */
	for (j=1;j<=ncolsein;j++)
	{
		/* copying the y coordinates and the function values and derivatives into the temporary arrays */
		for (k=nrowsein;k>=1;k--)
		{
			x[nrowsein-k+1]  = CoordYe[k][j]>nodataein?CoordYe[k][j]:
			ylleinrot+cseiny*(nrowsein-k);
			y[nrowsein-k+1]  = ze[k][j];
			y2[nrowsein-k+1] = y2a[k][j];
		}
		klo=1;khi=1;
		if ((equalarea)&&(y[nrowsein]>nodataein+.01))
		{
			/* weil das Stations-Netz das Bezugskoordinatensystem enthält, kann man hier auf sinbeta und cobeta verzichten */
			x[nrowsein+1]  = x[nrowsein]+cseiny;
			y[nrowsein+1]  = ze[1][j];
			y2[nrowsein+1] = y2a[1][j];
		}
		else
		{
			x[nrowsein+1]  = nodataaus;
			y[nrowsein+1]  = nodataaus;
			y2[nrowsein+1] = nodataaus;
		}
		CYa=ylltemp+(nrtemp-(equalarea==0?1:0)+1)*csaus;
		rowindex=2;
		for (m=1;m<=nrtemp;m++)
		{
			CYa-=csaus;
			ycNRAbyNCE[m][j]=CYa;
			if (rowindex<nrowsein)
			{
				if (CYa<CoordYe[rowindex][j])
				{
					rowindex++;
				}
			}
			if (rowindex<=nrowsein)
			{
				tx1=(CYa-CoordYe[rowindex-1][j])/(CoordYe[rowindex][j]-CoordYe[rowindex-1][j]);
				CXa=(1.0-tx1)*CoordXe[rowindex-1][j]+tx1*CoordXe[rowindex][j];
			}
			else
			{
				CXa=CoordXe[rowindex][j];
			}
			xcNRAbyNCE[m][j]=CXa;
			EvalCubicSpline(x,y,&klo,&khi,y2,nrowsein,CYa,&(yNRAbyNCE[m][j]),nodataein,nodataaus);
		}
	}

	//  tst=fopen("tst1.asc","wt");
	//  fprintf(tst,"ncols       %d\n"
	//              "nrows       %d\n"
	//              "xll         %g\n"
	//              "yll         %g\n"
	//              "cs          %g\n"
	//              "nodatavalue %g\n",ncolsein,nrtemp,xllaus,yllaus,csaus,nodataaus);
	//  for (i=1;i<=nrtemp;i++)
	//  {
	//    for (j=1;j<=ncolsein;j++)
	//      fprintf(tst,"%.4g\t",yNRAbyNCE[i][j]);
	//    fprintf(tst,"\n");
	//  }
	//  fclose(tst);
	//  
	//  tst=fopen("tst1x.asc","wt");
	//  fprintf(tst,"ncols       %d\n"
	//              "nrows       %d\n"
	//              "xll         %g\n"
	//              "yll         %g\n"
	//              "cs          %g\n"
	//              "nodatavalue %g\n",ncolsein,nrtemp,xllaus,yllaus,csaus,nodataaus);
	//  for (i=1;i<=nrtemp;i++)
	//  {
	//    for (j=1;j<=ncolsein;j++)
	//      fprintf(tst,"%.6g\t",xcNRAbyNCE[i][j]);
	//    fprintf(tst,"\n");
	//  }
	//  fclose(tst);
	//
	//  tst=fopen("tst1y.asc","wt");
	//  fprintf(tst,"ncols       %d\n"
	//              "nrows       %d\n"
	//              "xll         %g\n"
	//              "yll         %g\n"
	//              "cs          %g\n"
	//              "nodatavalue %g\n",ncolsein,nrtemp,xllaus,yllaus,csaus,nodataaus);
	//  for (i=1;i<=nrtemp;i++)
	//  {
	//    for (j=1;j<=ncolsein;j++)
	//      fprintf(tst,"%.6g\t",ycNRAbyNCE[i][j]);
	//    fprintf(tst,"\n");
	//  }
	//  fclose(tst);

	for (i=nrtempmin;i<=nrtempmax;i++)
	{
		anf=1; end=0;
		while (anf<ncolsein) 
		{
			anf=end+1;
			/* suchen Anfang ab letzter end-Position */
			while ((anf<=ncolsein)&&(yNRAbyNCE[i][anf]<nodataein+.01)) 
			{
				y2col[anf]=nodataaus;
				anf++;
			}
			/* suchen Ende ab aktueller anf-Position */
			end=anf;
			while ((end<ncolsein)&&(yNRAbyNCE[i][end+1]>nodataein+.01)) 
			end++;

			for (k=anf;k<=end && k<=ncolsein;k++)
			{
				xcol[k]=xcNRAbyNCE[i][k];
				ycol[k]=yNRAbyNCE[i][k];       
			}
			/* the last edge of a grid sqare sequence is taken as the same value as the forelast, 
		in order to be able to interpolate also within single grid cells (and not to have to write nodata values there) */
			if (anf<=ncolsein) /* at least 2 zells are needed */
			{
				if (equalarea) 
				{
					xcol[end+1]=xlleinrot+cseinx*end;
					ycol[end+1]=ycol[end];
				}
				ConstructCubicSpline(xcol,ycol,anf,end+equalarea,(double)1e30,(double)1e30,y2col);
			}
		}

		/* now evaluate the row splines */
		klo=1;khi=1;
		if ((equalarea)&&(yNRAbyNCE[i][ncolsein]>nodataein+.01)) 
		{
			yNRAbyNCE[i][ncolsein+1]=yNRAbyNCE[i][ncolsein];
		}
		else
		{
			yNRAbyNCE[i][ncolsein+1]=nodataaus;
		}
		CXa=xlltemp+csaus*(nctempmin);
		colindex=2;
		for (m=nctempmin;m<=nctempmax;m++)
		{
			CXa+=csaus;
			if (LookUpTable[i-nrtempmin+1][m-nctempmin+1]>0) 
			EvalCubicSpline(xcol,yNRAbyNCE[i],&klo,&khi,y2col,ncolsein,CXa,&(tempgrid[i-nrtempmin+1][m-nctempmin+1]),nodataaus,nodataein);
		}
	}

	//  tst=fopen("tst2.asc","wt");
	//  fprintf(tst,"ncols       %d\n"
	//              "nrows       %d\n"
	//              "xll         %g\n"
	//              "yll         %g\n"
	//              "cs          %g\n"
	//              "nodatavalue %g\n",nctempmax-nctempmin+1,nrtempmax-nrtempmin+1,xlltempmin+csaus*(nctempmin-1),yllein+csaus*(nrtemp-nrtempmax),csaus,nodataaus);
	//  for (i=nrtempmin;i<=nrtempmax;i++)
	//  {
	//    for (j=1;j<=nctemp;j++)
	//      fprintf(tst,"%.4g\t",tempgrid[i][j]);
	//    fprintf(tst,"\n");
	//  }
	//  fclose(tst);

	ReRotateResults(za);

}


