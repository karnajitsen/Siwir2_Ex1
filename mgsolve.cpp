#include<iostream>
#include<fstream>
#include "Grid.h"
#include <stdlib.h>
#include <ctime>
#include <string>
#include <malloc.h>
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

Grid** initialize(double hsize, const size_t level, bool flag)
{
	size_t je = level;
	size_t gdim = pow(2, je) + 1;
	//bool flag = true;

	Grid** arrGrid = (Grid**)memalign(ALLIGNMENT, level*sizeof(Grid*));
	for (size_t i = 0; i < level; i++)
	{
		arrGrid[i] = new Grid(gdim, gdim, hsize, hsize,flag);
		gdim = pow(2, --je) + 1;
		hsize *= 2.0;
		flag = false;
	}
	return arrGrid;

}

inline void restriction(const Grid * xgrd, const Grid * fgrd, Grid* rgrid)
{
	size_t xlen = (*xgrd).getXsize() - 1;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
    double	alpha = 1.0/hx/hx;
    double	beta = 1.0/hx/hx;
	double	center =  (2.0 * alpha) + (2.0 * beta);
	
    Grid tmpgrd(xlen + 1, xlen + 1,hx,hx,false);
    for (size_t i = 1; i < xlen; i++)
    {
        for (size_t j = 1; j < xlen; j++)
        {
            tmpgrd(j, i) = (*fgrd)(j, i) + alpha*((*xgrd)(j + 1, i) + (*xgrd)(j - 1, i)) + beta * ((*xgrd)(j, i + 1)
                + (*xgrd)(j, i - 1)) - (*xgrd)(j, i) * center;
        }
    }
	
	size_t rlen = (*rgrid).getXsize() - 1;
	
	for (size_t i = 1; i < rlen; i++)
	{
		for (size_t j = 1; j < rlen; j++)
		{
			(*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
				tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1) +
				2.0*(tmpgrd(2 * j, 2* i - 1) + tmpgrd(2 * j,  2*i + 1) +
				tmpgrd(2 * j - 1, 2* i) + tmpgrd(2 * j + 1,  2*i)) + 4.0 * tmpgrd(2 * j,  2*i)) / 16.0;
		}
	}

}

inline void interpolate(Grid * srcgrd, Grid * tgtgrd)
{
	size_t len = (*srcgrd).getXsize() -1;
	size_t txlen = (*tgtgrd).getXsize();
	double hx = (*tgtgrd).getHx();
	Grid tmpgrd(txlen, txlen, hx, hx,false);
		
	for (size_t j = 0; j < len; j++)
	{
		size_t k = 2 * j;
		for (size_t i = 0; i < len; i++)
		{
			size_t l = 2 * i;
			tmpgrd(l, k) = (*srcgrd)(i, j);
			tmpgrd(l, k + 1) = 0.5*((*srcgrd)(i, j) + (*srcgrd)(i, j + 1));
			tmpgrd(l + 1, k) = 0.5*((*srcgrd)(i, j) + (*srcgrd)(i + 1, j));
			tmpgrd(l + 1, k + 1) = 0.25*((*srcgrd)(i, j) + (*srcgrd)(i + 1, j) + (*srcgrd)(i, j + 1)
											+(*srcgrd)(i + 1, j + 1));			
		}
	}

	for (size_t i = 1; i < txlen - 1; i++)
	{
		for (size_t j = 1; j < txlen - 1; j++)
		{
			(*tgtgrd)(j,i) += tmpgrd(j, i);	

		}
	}

}

inline void smooth(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize();
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
    double	center = (2.0 * alpha + 2.0 * beta);
	
	for (size_t i = 0; i < iter; i++)
	{
		for (size_t j = 1; j < dimX - 1; j++)
		{
			for (size_t k = ((j + 1) & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hx*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
                    +(*xgrd)(k, j - 1)))/center;

			}

		}

		for (size_t j = 1; j < dimX - 1; j++)
		{
			for (size_t k = (j & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hx*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) / center;

				
			}

		}
	}
}

inline void resdualNorm(Grid* xgrd, const Grid * fgrd, double* norm)
{

	size_t dimX = (*xgrd).getXsize() - 1;
	//double hx = (*xgrd).getHx() ;
	double r = 0.0;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
    double	alpha = 1.0;
    double	beta = 1.0;
	double	center = (2.0 * alpha + 2.0 * beta);

	*norm = 0.0;

	for (size_t j = 1; j < dimX; j++)
	{
		for (size_t k = 1; k < dimX; k++)
		{
            r = hx*hx*(*fgrd)(k,j) + alpha*((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
				+ (*xgrd)(k, j - 1)) - (*xgrd)(k,j) * center;
          
			*norm += r*r;
		}

	}
	
	*norm = sqrt(*norm / (dimX+1) / (dimX+1));
}

inline void errorNorm(Grid* xgrd, const Grid * sgrd, double* norm)
{

	size_t dimX = (*xgrd).getXsize() - 1;
	//double hx = (*xgrd).getHx() ;
	double r = 0.0;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
	double	center = (2.0 * alpha + 2.0 * beta);

	*norm = 0.0;

	for (size_t j = 1; j < dimX; j++)
	{
		for (size_t k = 1; k < dimX; k++)
		{
			r = (*sgrd)(k + 1, j) - (*xgrd)(k, j);

		*norm += r*r;
		}

	}

	*norm = sqrt(*norm / (dimX + 1) / (dimX + 1));
}

inline double gxy(const double x, const double y)
{
	return sin(M_PI * x) * sinh(M_PI * y);
}

int main(int argc, char** argv)
{

	//std::cout << "1";
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	clock_t tim;
	size_t level = atoi(argv[1]);
	size_t vcycle = atoi(argv[2]);
	size_t gdim = pow(2, level) + 1;
	size_t vdim = gdim - 2;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	
	Grid ** xGrids = initialize(hsize, level, true);
	Grid ** fGrids = initialize(hsize, level,false);
	Grid sGrid(gdim, gdim, hsize, hsize,true);

	for (size_t i = 0; i < gdim; i++)
	{
		for (size_t j = 0; j < gdim; j++)
		{
			sGrid(j, i) = sGrid.gxy(j*hsize, i*hsize);
		}
	}
	
   	tim = clock();

	smooth(xGrids[0], fGrids[0], V1);
	for (size_t i = 0; i < vcycle; i++)
    {
		//size_t jl = 0;
		//restriction(xGrids[0], fGrids[0], fGrids[1]);
		
		for (size_t jl = 0; jl < level - 1; jl++)
		{
			smooth(xGrids[jl], fGrids[jl], V1);
			restriction(xGrids[jl], fGrids[jl], fGrids[jl + 1]);
		}
		
		for (size_t j = level - 1; j > 0; j--)
		{	
            smooth(xGrids[j], fGrids[j], V2);
			interpolate(xGrids[j], xGrids[j - 1]);
			(*xGrids[j]).reset();
			(*fGrids[j]).reset();					
		}
		
		//smooth(xGrids[0], fGrids[0], V1);
        oldnorm = newnorm;
		resdualNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
            convrate = newnorm / oldnorm;
		
		std::cout << "Residual Norm after " << i + 1 << " V-Cycle = " << newnorm << '\n';
		std::cout << "Covergence rate after " << i + 1 << " V-Cycle = " << convrate << '\n';
     }

	tim = clock() - tim;

	std::cout << "Time spend for all V - cycles= " << ((float)tim) / CLOCKS_PER_SEC << '\n';

	errorNorm(xGrids[0], &sGrid, &newnorm);
	std::cout << "Error Norm for h as 1/" << gdim-1 <<" = " << newnorm << '\n';
	std::string fname = std::string("data/solution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOut(fname);
	std::string fnames = std::string("data/exactsolution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOutsolt(fnames);
	for (int y = 0; y < gdim; ++y) {
		for (int x = 0; x < gdim; ++x) {

			fOut << x*hsize << "\t" << y*hsize << "\t" << (*xGrids[0])(x, y) << std::endl;
			fOutsolt << x*hsize << "\t" << y*hsize << "\t" << sGrid(x, y) << std::endl;
		}
		fOut << std::endl;
		fOutsolt << std::endl;
	}
	fOut.close();
	fOutsolt.close();

	return 0;
}
