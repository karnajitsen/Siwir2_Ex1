#include<iostream>
#include "Grid.h"
//#include "Stencil.h"
#include "stdio.h"
#include <ctime>


#define M_PI 3.14
#include <cmath>

#define XDOMLOW 0
#define XDOMHIGH 1
#define YDOMLOW 0
#define YDOMHIGH 1
#define TOLERR 0.01
#define V1 2
#define V2 1

using namespace std;

inline void calculatefx(double * fvec)
{

}

inline double gxy(const double x, const double y)
{
	return sin(M_PI * x) * sinh(M_PI * y);
}

Grid** initialize( const double hsize, const int level)
{
	Grid** arrGrid = NULL;
	int je = level;
	
	for (int i = 0; i < level; i++)
	{
		int gdim = pow(2, je--) + 1;
		arrGrid[i] = new Grid(gdim, gdim, hsize, hsize);
		Grid grd(gdim, gdim, hsize, hsize);
		for (int j = 0; j < gdim; j++)
		{
			double k = j*hsize;
			grd(0, j) = gxy(0,k);
			grd(j, 0) = gxy(k,0);
			grd(j, gdim - 1) = gxy(k, gdim - 1);
			grd(gdim - 1, j) = gxy(gdim - 1, k);
		}
		arrGrid[i] = &grd;
	}
	return arrGrid;

}

inline void calresidual(const Grid * grd, const double * fvec, double * rvec)
{
	int xlen = (*grd).getXsize()-1, m=0;
	for (int i = 1; i < xlen; i++)
	{
		for (int j = 1; j < xlen; j++,m++)
		{
			rvec[m] = fvec[m] - 4.0*(*grd)(i, j) - (*grd)(i, j-1) - (*grd)(i, j+1) - (*grd)(i-1, j) - (*grd)(i+1, j);
		}
	}
}

inline void rstrToCoarse( double * rvec, int const totdim)
{
	for (int i = 1, j=0; i < totdim; i += 2,j++)
	{
		rvec[j] = 0.5 * rvec[i] + 0.25 * (rvec[i - 1] + rvec[i + 1]);
	}
}

inline void interpolate(const Grid * grd, Grid * tmpgrd)
{
	int len = (*grd).getXsize();
	int hx = (*grd).getHx()/2.0;
	int nlen = len * 2 - 1;
	tmpgrd = new Grid(nlen,nlen,hx,hx);
	for (int i = 0; i < len; i++)
	{
		int k = 2 * i;
		for (int j = 0; j < len; j++)
		{
			int l = 2 * j;
			(*tmpgrd)(k,l) = (*grd)(i,j);
			(*tmpgrd)(k + 1, l) = 0.5*((*grd)(i, j) + (*grd)(i+1, j));
			(*tmpgrd)(k, l + 1) = 0.5*((*grd)(i, j) + (*grd)(i , j+1));
			(*tmpgrd)(k + 1, l + 1) = 0.25*((*grd)(i, j) + (*grd)(i + 1, j) + (*grd)(i, j + 1) + +(*grd)(i + 1, j + 1));
		}
	}
}

inline void rbgs(Grid* grd, const double * fvec, const int iter)
{
	int dimX = (*grd).getXsize();
	//double temp = 0.0;

	for (int i = 0; i < iter; i++)
	{
		for (int j = 1; j < dimX; j++)
		{
			for (int k = (j & 1) + 1; k < dimX; k+=2)
			{
				(*grd)(j, k) = 0.25*((*grd)(j + 1, k) + (*grd)(j - 1, k) + (*grd)(j, k + 1) + (*grd)(j, k - 1));
			}
			
		}

		for (int j = 1; j < dimX; j++)
		{
			for (int k = j & 1; k < dimX; k += 2)
			{
				(*grd)(j, k) = 0.25*((*grd)(j + 1, k) + (*grd)(j - 1, k) + (*grd)(j, k + 1) + (*grd)(j, k - 1));
			}

		}
	}

}

inline void calNorm(Grid* grd, const double * fvec, double* norm)
{
	int dimX = (*grd).getXsize();
	double r = 0.0;
	
	for (int j = 1; j < dimX; j++)
	{
		for (int k = (j & 1) + 1; k < dimX; k += 2)
		{
			r = 0.25*((*grd)(j + 1, k) + (*grd)(j - 1, k) + (*grd)(j, k + 1) + (*grd)(j, k - 1));
			*norm += r*r;
		}

	}

	for (int j = 1; j < dimX; j++)
	{
		for (int k = j & 1; k < dimX; k += 2)
		{
			r = 0.25*((*grd)(j + 1, k) + (*grd)(j - 1, k) + (*grd)(j, k + 1) + (*grd)(j, k - 1));
			*norm += r*r;
		}

	}

	*norm = sqrt(*norm / dimX / dimX);
}

int main(int argc, char** argv)
{
  
	std::cout << "1";
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(EXIT_FAILURE);
	}

	std::cout << "1";
	clock_t tim;
	int level = atoi(argv[0]);
	int vcycle = atoi(argv[1]);
	size_t gdim = pow(2,level) + 1;
	size_t vdim = gdim - 2;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1);
	size_t totdim = vdim*vdim;
	Grid ** grids = initialize(hsize,level);
	double * fvec = NULL, *rvec = NULL;
	fvec = new double[totdim]();

	tim = clock();
	std::cout << "2";
	for (int i = 0; i < vcycle; i++)
	{		
		rbgs(grids[0],fvec,V1);
		calresidual(grids[0], fvec, rvec);
		rstrToCoarse(rvec, totdim);
		for (int j = 1; j < level; j++)
		{
			vdim = vdim - 2;
			totdim = vdim*vdim;
			rbgs(grids[j], rvec, V1);
			calresidual(grids[j], rvec, rvec);
			rstrToCoarse(rvec, totdim);
		}

		rbgs(grids[level+1], rvec, V1);
		Grid * tmpgrid = NULL;
		interpolate(grids[level + 1], tmpgrid);
		//free(grids[level + 1]);
		for (int j = level; j > 0; j--)
		{  
			tmpgrid = NULL;
			rbgs(grids[j], rvec, V2);
			interpolate(grids[j], tmpgrid);
			//free(grids[j]);
		}
		oldnorm = newnorm;
		calNorm(grids[0], fvec, &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;
		//if (norm < TOLERR)
			//break;

		std::cout << "Residual after " << i + 1 << "V-Cycle = " << newnorm;
		std::cout << "Covergence rate after " << i + 1 << "V-Cycle = " << convrate;
	}
	
	tim = clock() - tim;

	std::cout << "Time spend for two V - cycles= " << tim;
	return 0;
}
