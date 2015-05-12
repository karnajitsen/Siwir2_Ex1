#include<iostream>
#include "Grid.h"
#include <stdlib.h>
//#include "Stencil.h"
#include "stdio.h"
#include <ctime>
#include <memory>
#include <malloc.h>


//#define M_PI 3.14
#include <cmath>
#include <stdlib.h>
#define XDOMLOW 0.0
#define XDOMHIGH 1.0
#define YDOMLOW 0.0
#define YDOMHIGH 1.0
#define TOLERR 0.01
#define V1 2
#define V2 1

using namespace std;

Grid** initialize(double hsize, const size_t level)
{
	size_t je = level;
	size_t gdim = pow(2, je) + 1;
	Grid** arrGrid = NULL;
	posix_memalign((void **)arrGrid, ALLIGNMENT, level*sizeof(Grid*));
	for (size_t i = 0; i < level; i++)
	{
		//std::cout << &arrGrid[i] << "\n";
		arrGrid[i] = new Grid(gdim, gdim, hsize, hsize);
		gdim = pow(2, --je) + 1;
		hsize *= 2.0;
	}
	return arrGrid;

}

inline void restriction(const Grid * grd, const Grid * fgrd, Grid* rgrid)
{
	size_t xlen = (*grd).getXsize() - 1, m = 0;
	Grid tmpgrd(xlen + 1, xlen + 1, (*grd).getHx(), (*grd).getHx());
	for (size_t i = 1; i < xlen; i++)
	{
		for (size_t j = 1; j < xlen; j++, m++)
		{	
			//std::cout << "j= " << j << " " <<"888888\n ";
			tmpgrd(j,i) = (*fgrd)(j,i) - 4.0*((*grd)(j, i) - (*grd)(j, i - 1) - (*grd)(j, i + 1) 
					      - (*grd)(j - 1, i) - (*grd)(j + 1, i));
			
		}
	}
	
	size_t rlen = (*rgrid).getXsize() - 1;
	//std::cout << "j= " << (*rgrid).getXsize() << " " << "666666\n ";
	for (size_t i = 1; i < rlen; i++)
	{
		for (size_t j = 1; j < rlen; j++)
		{
			//std::cout << "j= " << j << " " << "999999\n ";
			(*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
				tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1) +
				2.0*(tmpgrd(2 * j, 2 * i - 1) + tmpgrd(2 * j, 2 * i + 1) +
				tmpgrd(2 * j - 1, 2 * i) + tmpgrd(2 * j + 1, 2 * i)) + 4 * tmpgrd(2 * j, 2 * i)) / 16.0;
		}
	}
}

inline void interpolate(const Grid * srcgrd, const Grid * tgtgrd)
{
	size_t len = (*srcgrd).getXsize() - 1;
	size_t hx = (*srcgrd).getHx() / 2.0;
	size_t nlen = len * 2 - 1;
	//Grid * tmpgrd = new Grid(nlen, nlen, hx, hx);
	for (size_t j = 0; j < len; j++)
	{
		size_t k = 2 * j;
		for (size_t i = 0; i < len; i++)
		{
			size_t l = 2 * i;
			(*tgtgrd)(k, l) += (*srcgrd)(i, j);
			(*tgtgrd)(k + 1, l) += 0.5*((*srcgrd)(i, j) + (*srcgrd)(i + 1, j));
			(*tgtgrd)(k, l + 1) += 0.5*((*srcgrd)(i, j) + (*srcgrd)(i, j + 1));
			(*tgtgrd)(k + 1, l + 1) += 0.25*((*srcgrd)(i, j) + (*srcgrd)(i + 1, j) + (*srcgrd)(i, j + 1)
											+(*srcgrd)(i + 1, j + 1));
			
		}
	}
}

inline void rbgs(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize();
	//double temp = 0.0;
	//std::cout << dimX << "#####";

	for (size_t i = 0; i < iter; i++)
	{
		for (size_t j = 1; j < dimX - 1; j++)
		{
			for (size_t k = ((j + 1) & 0x1) + 1; k < dimX - 1; k += 2)
			{
				//std::cout << "k= " << k << " " <<"888888\n ";
				(*xgrd)(k,j) = 0.25*((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1) + (*xgrd)(k, j - 1));
			}

		}

		for (size_t j = 1; j < dimX - 1; j++)
		{
			for (size_t k = (j & 0x1) + 1; k < dimX - 1; k += 2)
			{
				//std::cout << "k= " << k << " " << "999999\n ";
				(*xgrd)(k,j) = 0.25*((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1) + (*xgrd)(k, j - 1));
			}

		}
	}

}

void calNorm(Grid* xgrd, const Grid * fgrd, double* norm)
{
	size_t dimX = (*xgrd).getXsize() - 1;
	double hX = (*xgrd).getHx() ;
	double r = 0.0;

	*norm = 0.0;

	for (size_t j = 1; j < dimX; j++)
	{
		for (size_t k = 1; k < dimX; k++)
		{
			r = hX*(*fgrd)(j, k) - 4.0*(*xgrd)(j,k) + (*xgrd)(j + 1, k) + (*xgrd)(j - 1, k) + (*xgrd)(j, k + 1) 
				+ (*xgrd)(j, k - 1);
			*norm += r*r;
		}

	}
	
	*norm = sqrt(*norm / dimX / dimX);
}

int main(int argc, char** argv)
{

	//std::cout << "1";
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	//std::cout << argv[1] << " " << argv[2];
	clock_t tim;
	size_t level = atoi(argv[1]);
	size_t vcycle = atoi(argv[2]);
	size_t gdim = pow(2, level) + 1;
	size_t vdim = gdim - 2;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	//size_t totdim = vdim*vdim;

	Grid ** xGrids = initialize(hsize, level);
	Grid ** fGrids = initialize(hsize, level);
	//fvec = new double[totdim]();

	tim = clock();

	
	for (size_t i = 0; i < vcycle; i++)
	{
		rbgs(xGrids[0], fGrids[0], V1);
		//std::cout << "3";
				
		restriction(xGrids[0], fGrids[0],xGrids[1]);
		size_t jl = 0;
		//rstrToCoarse(rvec, totdim);
		for (jl = 1; jl < level-1; jl++)
		{
			vdim = vdim - 2;
			totdim = vdim*vdim;
			rbgs(xGrids[jl], fGrids[jl], V1); 
			restriction(xGrids[jl], fGrids[0], xGrids[jl+1]);
			//rstrToCoarse(rvec, totdim);
		}
		
		for (size_t j = level - 1; j > 0; j--)
		{	
			
			rbgs(xGrids[j], fGrids[j], V2);
			//std::cout << "66" << "\n";
			interpolate(xGrids[j], xGrids[j-1]);
			
		}
		//std::cout << "55" << "\n";
		oldnorm = newnorm;
		calNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;
		
		std::cout << "Residual after " << i + 1 << "V-Cycle = " << newnorm << '\n';
		std::cout << "Covergence rate after " << i + 1 << "V-Cycle = " << convrate << '\n';
	}

	tim = clock() - tim;

	std::cout << "Time spend for two V - cycles= " << ((float)tim) / CLOCKS_PER_SEC << '\n';

	return 0;
}
