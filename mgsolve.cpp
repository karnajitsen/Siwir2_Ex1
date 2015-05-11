#include<iostream>
#include "Grid.h"
#include <stdlib.h>
//#include "Stencil.h"
#include "stdio.h"
#include <ctime>


#define M_PI 3.14
#include <cmath>

#define XDOMLOW 0.0
#define XDOMHIGH 1.0
#define YDOMLOW 0.0
#define YDOMHIGH 1.0
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

Grid** initialize(double hsize, const int level)
{
	int je = level;
	int gdim = pow(2, je) + 1;
	Grid** arrGrid = (Grid**)_aligned_malloc(level*sizeof(Grid*), 1024);
	for (int i = 0; i < level; i++)
	{
		std::cout << &arrGrid[i] << "\n";
		arrGrid[i] = new Grid(gdim, gdim, hsize, hsize);
		for (int j = 0.0; j < gdim; j++)
		{
			double k = j*hsize;
			double l = (gdim - 1.0)*hsize;
			(*arrGrid[i])(j, 0) = gxy(k, 0.0);
			(*arrGrid[i])(0, j) = gxy(0.0, k);
			(*arrGrid[i])(j, gdim - 1) = gxy(k, l);
			(*arrGrid[i])(gdim - 1, j) = gxy(l, k);

			//cout << j << " " << gdim - 1 << " " << grd(j, 0) << " " << grd(0, j) << " " << grd(j, gdim - 1) << " " << grd(gdim - 1, j) << "\n\n";
		}
		
		gdim = pow(2, --je) + 1;
		hsize *= 2.0;
	}
	return arrGrid;

}

inline void restriction(const Grid * grd, const Grid * fgrd, Grid* rgrid)
{
	int xlen = (*grd).getXsize() - 1, m = 0;
	Grid tmpgrd(xlen + 1, xlen + 1, 0.0, 0.0);
	for (int i = 1; i < xlen; i++)
	{
		for (int j = 1; j < xlen; j++, m++)
		{
			tmpgrd(i,j) = (*fgrd)(i,j) - 4.0*(*grd)(i, j) - (*grd)(i, j - 1) - (*grd)(i, j + 1) - (*grd)(i - 1, j) - (*grd)(i + 1, j);
		}
	}

	int rlen = (*rgrid).getXsize() - 2;
	for (int i = 1; i < rlen; i++)
	{
		for (int j = 1; j < rlen; j++)
		{
			(*rgrid)(i, j) = (tmpgrd(2 * i - 1, 2 * j - 1) + tmpgrd(2 * i - 1, 2 * j + 1) +
				tmpgrd(2 * i + 1, 2 * j - 1) + tmpgrd(2 * i + 1, 2 * j + 1) +
				2.0*(tmpgrd(2 * i, 2 * j - 1) + tmpgrd(2 * i, 2 * j + 1) +
				tmpgrd(2 * i - 1, 2 * j) + tmpgrd(2 * i + 1, 2 * j)) + 4 * tmpgrd(2 * i, 2 * j)) / 16.0;
		}
	}
}

inline void restriction(double * rvec, int const totdim)
{
	for (int i = 1, j = 0; i < totdim; i += 2, j++)
	{
		rvec[j] = 0.5 * rvec[i] + 0.25 * (rvec[i - 1] + rvec[i + 1]);
	}
}

inline void interpolate(const Grid * grd, Grid * tmpgrd)
{
	int len = (*grd).getXsize();
	int hx = (*grd).getHx() / 2.0;
	int nlen = len * 2 - 1;
	tmpgrd = new Grid(nlen, nlen, hx, hx);
	for (int i = 0; i < len; i++)
	{
		int k = 2 * i;
		for (int j = 0; j < len; j++)
		{
			int l = 2 * j;
			(*tmpgrd)(k, l) = (*grd)(i, j);
			(*tmpgrd)(k + 1, l) = 0.5*((*grd)(i, j) + (*grd)(i + 1, j));
			(*tmpgrd)(k, l + 1) = 0.5*((*grd)(i, j) + (*grd)(i, j + 1));
			(*tmpgrd)(k + 1, l + 1) = 0.25*((*grd)(i, j) + (*grd)(i + 1, j) + (*grd)(i, j + 1) + +(*grd)(i + 1, j + 1));
		}
	}
}

inline void rbgs(Grid* xgrd, const Grid* fgrd, const int iter)
{
	int dimX = (*xgrd).getXsize();
	//double temp = 0.0;

	for (int i = 0; i < iter; i++)
	{
		for (int j = 1; j < dimX-1; j++)
		{
			for (int k = ((j + 1) & 0x1) + 1; k < dimX-1; k += 2)
			{
				std::cout << "888888\n ";
				(*xgrd)(j, k) = 0.25*((*xgrd)(j + 1, k) + (*xgrd)(j - 1, k) + (*xgrd)(j, k + 1) + (*xgrd)(j, k - 1));
			}

		}

		for (int j = 1; j < dimX-1; j++)
		{
			for (int k = (j & 0x1)+1; k < dimX-1; k += 2)
			{
				std::cout << "99999";
				(*xgrd)(j, k) = 0.25*((*xgrd)(j + 1, k) + (*xgrd)(j - 1, k) + (*xgrd)(j, k + 1) + (*xgrd)(j, k - 1));
			}

		}
	}

}

inline void calNorm(Grid* xgrd, const Grid * fgrd, double* norm)
{
	int dimX = (*xgrd).getXsize();
	double r = 0.0;

	for (int j = 1; j < dimX; j++)
	{
		for (int k = (j & 1) + 1; k < dimX; k += 2)
		{
			r = 0.25*((*xgrd)(j + 1, k) + (*xgrd)(j - 1, k) + (*xgrd)(j, k + 1) + (*xgrd)(j, k - 1));
			*norm += r*r;
		}

	}

	for (int j = 1; j < dimX; j++)
	{
		for (int k = j & 1; k < dimX; k += 2)
		{
			r = 0.25*((*xgrd)(j + 1, k) + (*xgrd)(j - 1, k) + (*xgrd)(j, k + 1) + (*xgrd)(j, k - 1));
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
		exit(0);
	}

	//std::cout << argv[1] << " " << argv[2];
	clock_t tim;
	int level = atoi(argv[1]);
	int vcycle = atoi(argv[2]);
	size_t gdim = pow(2, level) + 1;
	size_t vdim = gdim - 2;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	size_t totdim = vdim*vdim;

	Grid ** xGrids = initialize(hsize, level);
	Grid ** fGrids = initialize(hsize, level);
	//fvec = new double[totdim]();

	tim = clock();
	

	for (int i = 0; i < vcycle; i++)
	{

		std::cout << "3";
		rbgs(xGrids[0], fGrids[0], V1);
		std::cout << "3";
		restriction(xGrids[0], fGrids[0],xGrids[1]);
		//rstrToCoarse(rvec, totdim);
		for (int j = 1; j < level; j++)
		{
			vdim = vdim - 2;
			totdim = vdim*vdim;
			rbgs(xGrids[j], fGrids[j], V1);
			restriction(xGrids[j], fGrids[0], xGrids[j+1]);
			//rstrToCoarse(rvec, totdim);
		}

		rbgs(xGrids[level + 1], fGrids[level + 1], V1);
		Grid * tmpgrid = NULL;
		interpolate(xGrids[level + 1], tmpgrid);
		for (int j = level; j > 0; j--)
		{
			tmpgrid = NULL;
			rbgs(xGrids[j], fGrids[j], V2);
			interpolate(xGrids[j], tmpgrid);
		}
		oldnorm = newnorm;
		calNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;
		
		std::cout << "Residual after " << i + 1 << "V-Cycle = " << newnorm;
		std::cout << "Covergence rate after " << i + 1 << "V-Cycle = " << convrate;
	}

	tim = clock() - tim;

	std::cout << "Time spend for two V - cycles= " << tim;
	return 0;
}
