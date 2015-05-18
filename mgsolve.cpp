#include "Grid.h"
#include<iostream>
#include<fstream>
#include <sys/time.h>
#define XDOMLOW 0.0
#define XDOMHIGH 1.0
#define YDOMLOW 0.0
#define YDOMHIGH 1.0
#define TOLERR 0.00001
#define V1 2
#define V2 1

Grid ** xGrids = nullptr;
Grid ** fGrids = nullptr;
Grid *sGrid = nullptr;
size_t *ndflag = 1;

void init(double hsize, const size_t level, bool dirflg)
{
	size_t je = level;
	size_t ydim = pow(2, je) + 1;
	size_t xdim = ydim;
	bool flag = true;
	xGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
	fGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
	for (size_t i = 0; i < level; i++)
	{
		xGrids[i] = new Grid(xdim, ydim, hsize, hsize, flag, dirflg);
		fGrids[i] = new Grid(xdim, ydim, hsize, hsize, false, dirflg);
		ydim = pow(2, --je) + 1;
		xdim = ydim;
		hsize *= 2.0;
		flag = false;
	}


	for (size_t i = 1; i < (*fGrids[0]).getYsize()-1; i++)
	{
		//(*fGrids[0])(0, i) = 1.0;
		//(*fGrids[0])((*fGrids[0]).getXsize() - 1, i) = 1.0;
		for (size_t j = 1; j < (*fGrids[0]).getXsize(); j++)
		{
			(*fGrids[0])(j, i) = 2.0;
		}
	}
}


inline void smooth(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize();
	size_t dimY = (*xgrd).getYsize();
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
	double	center = 1.0 / (2.0 * alpha + 2.0 * beta);

	for (size_t i = 0; i < iter; i++)
	{
		for (size_t j = 1; j < dimY - 1; j++)
		{
			for (size_t k = ((j + 1) & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;

			}

			if (*ndflag == 0 && l == 1 && xgrd == xGrids[0])
			{
				(*xgrd)(0, j) = -hx + (*xgrd)(1, j);
				(*xgrd)(dimX - 1, j) = -hx + (*xgrd)(dimX - 2, j);

			}

		}

		for (size_t j = 1; j < dimY - 1; j++)
		{
			for (size_t k = (j & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;


			}

			if (*ndflag == 0 && l == 1 && xgrd == xGrids[0])
			{
				(*xgrd)(0, j) = -hx + (*xgrd)(1, j);
				(*xgrd)(dimX - 1, j) = -hx + (*xgrd)(dimX - 2, j);

			}

		}
	}
}

void restriction(const Grid * xgrd, const Grid * fgrd, Grid* rgrid)
{
	size_t xlen = (*xgrd).getXsize() - 1;
	size_t ylen = (*xgrd).getYsize() - 1;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0 / hx / hx;
	double	beta = 1.0 / hy / hy;
	double	center = (2.0 * alpha) + (2.0 * beta);


	/*cout << "====B4 restriction residual=== \n\n";
	for (size_t j = 0; j <= xlen; j++)
	{
		for (size_t k = 0; k <= xlen; k++)
		{
			cout << (*xgrd)(k, j) << " ";
		}

		cout << '\n';
	}*/


	Grid tmpgrd(xlen + 1, ylen + 1, hx, hy, false,true);
	for (size_t i = 1; i < ylen; i++)
	{
		for (size_t j = 1; j < xlen; j++)
		{
			tmpgrd(j, i) = (*fgrd)(j, i) + alpha*((*xgrd)(j + 1, i) + (*xgrd)(j - 1, i)) + beta * ((*xgrd)(j, i + 1)
				+ (*xgrd)(j, i - 1)) - (*xgrd)(j, i) * center;
		}

		/*if (*ndflag == 0)
		{
			tmpgrd(0, i) = (*fgrd)(0, i) - (2.0 / hx) + 2.0 * alpha*((*xgrd)(1, i)) + beta * ((*xgrd)(0, i + 1)
				+ (*xgrd)(0, i - 1)) 
				- (*xgrd)(0, i) * center ;

			tmpgrd(xlen, i) = (*fgrd)(xlen, i) - (2.0 / hx) + 2.0 * alpha*((*xgrd)(xlen - 1, i)) + beta * ((*xgrd)(xlen, i + 1)
				+ (*xgrd)(xlen, i - 1)) 
				- (*xgrd)(xlen, i) * center ;
		}*/
	}

	//cout << "====After restriction residual=== \n\n";
	//for (size_t j = 0; j <= xlen; j++)
	//{
	//	for (size_t k = 0; k <= xlen; k++)
	//	{
	//		cout << tmpgrd(k, j) << " ";
	//	}

	//	cout << '\n';
	//}

	size_t rxlen = (*rgrid).getXsize() - 1;
	size_t rylen = (*rgrid).getYsize() - 1;

	

	for (size_t i = 1; i < rylen; i++)
	{

		for (size_t j = 1; j < rxlen; j++)
		{
			(*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
				tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1) +
				2.0*(tmpgrd(2 * j, 2 * i - 1) + tmpgrd(2 * j, 2 * i + 1) +
				tmpgrd(2 * j - 1, 2 * i) + tmpgrd(2 * j + 1, 2 * i)) + 4.0 * tmpgrd(2 * j, 2 * i)) / 16.0;
		}

		//if (*ndflag == 0)
		//{
		//	(*rgrid)(0, i) = 0.5 * tmpgrd(0, 2 * i) + 0.25 * tmpgrd(1, 2 * i) + 0.125*(tmpgrd(1, 2 * i - 1) + tmpgrd(1, 2 * i + 1));
		//		
		//		/*(2.0*(tmpgrd(1, 2 * i - 1) + tmpgrd(1, 2 * i + 1) + 2.0 * hx) +
		//		2.0*(tmpgrd(0, 2 * i - 1) + tmpgrd(0, 2 * i + 1) +
		//		2.0 * (hx + tmpgrd(1, 2 * i))) + 4.0 * tmpgrd(0, 2 * i)) / 16.0;*/

		//	(*rgrid)(rxlen, i) = 0.5 * tmpgrd(xlen, 2 * i) + 0.25 * tmpgrd(xlen - 1, 2 * i) + 0.125*(tmpgrd(xlen - 1, 2 * i - 1) + tmpgrd(xlen - 1, 2 * i + 1));
		//		
		//	/*(2.0*(tmpgrd(xlen - 1, 2 * i - 1) + tmpgrd(xlen - 1, 2 * i + 1) - 2.0 * hx) +
		//		2.0*(tmpgrd(xlen, 2 * i - 1) + tmpgrd(xlen, 2 * i + 1) +
		//		2.0 * (hx + tmpgrd(xlen - 1, 2 * i))) + 4.0 * tmpgrd(xlen, 2 * i)) / 16.0;*/
		//}
	}

	//cout << "====After restriction === \n\n";
	//for (size_t j = 0; j <= rxlen; j++)
	//{
	//	for (size_t k = 0; k <= rxlen; k++)
	//	{
	//		cout << (*rgrid)(k, j) << " ";
	//	}

	//	cout << '\n';
	//}

}

inline void interpolate(Grid * srcgrd, Grid * tgtgrd)
{
	size_t len = (*srcgrd).getXsize() - 1;
	size_t txlen = (*tgtgrd).getXsize();
	double hx = (*tgtgrd).getHx();
	Grid tmpgrd(txlen, txlen, hx, hx, false,true);

	//cout << "====b4 Interpolate === \n\n";
	//for (size_t j = 0; j <= len; j++)
	//{
	//	for (size_t k = 0; k <= len; k++)
	//	{
	//		cout << (*srcgrd)(k, j) << " ";
	//	}

	//	cout << '\n';
	//}

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
				+ (*srcgrd)(i + 1, j + 1));
		}

		/*if (*ndflag == 0)
		{
			tmpgrd(2 * len, 2*j + 1) = 0.5*((*srcgrd)(len, j) + (*srcgrd)(len, j + 1));
			tmpgrd(2 * len, 2*j) = (*srcgrd)(len, j);
		}*/
	}

	//cout << "====b4 Interpolate Add === \n\n";
	//for (size_t j = 0; j < txlen; j++)
	//{
	//	for (size_t k = 0; k < txlen; k++)
	//	{
	//		cout << tmpgrd(k, j) << " ";
	//	}

	//	cout << '\n';
	//}

	for (size_t i = 1; i < txlen - 1; i++)
	{
		for (size_t j = 1; j < txlen - 1; j++)
		{
			(*tgtgrd)(j, i) += tmpgrd(j, i);

		}
	}


	//cout << "====After Interpolate Add === \n\n";
	//for (size_t j = 0; j < txlen; j++)
	//{
	//	for (size_t k = 0; k < txlen; k++)
	//	{
	//		cout << (*tgtgrd)(k, j) << " ";
	//	}

	//	cout << '\n';
	//}

}

inline void resdualNorm(const Grid* xgrd, const Grid * fgrd, double* norm)
{

	size_t dimX = (*xgrd).getXsize() - 1;
	size_t dimY = (*xgrd).getYsize() - 1;
	double r = 0.0;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
	double	center = (2.0 * alpha + 2.0 * beta);

	*norm = 0.0;

	for (size_t j = 1; j < dimY; j++)
	{
		for (size_t k = 1; k < dimX; k++)
		{
			r = hx*hy*(*fgrd)(k, j) + alpha*((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
				+ (*xgrd)(k, j - 1)) - (*xgrd)(k, j) * center;

			*norm += r*r;
		}

		*norm = sqrt(*norm / (dimX - 1) / (dimY - 1));
}



inline void errorNorm(const Grid* xgrd, const Grid * sgrd, double* norm)
{

	size_t dimX = (*xgrd).getXsize();
	size_t dimY = (*xgrd).getYsize();
	double r = 0.0;
	*norm = 0.0;

	for (size_t j = 0; j < dimY; j++)
	{
		for (size_t k = 0; k < dimX; k++)
		{
			r = (*sgrd)(k , j) - (*xgrd)(k, j);

			*norm += r*r;
		}

	}

	*norm = sqrt(*norm / dimX / dimY);
}

inline void mgsolve(size_t level, size_t vcycle)
{
	size_t gdim = pow(2, level) + 1;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);

	init(hsize, level, false);
	sGrid = new Grid(gdim, gdim, hsize, hsize, true, false);

	for (size_t i = 0; i < gdim; i++)
	{
		for (size_t j = 0; j < gdim; j++)
		{
			(*sGrid)(j, i) = (*sGrid).gxy2(j*hsize);
		}
	}

	for (size_t i = 0; i < vcycle; i++)
	{
		//orthogonalize(xGrids[0]);

		for (size_t jl = 0; jl < level - 1; jl++)
		{
			//orthogonalize(fGrids[jl]);
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

		oldnorm = newnorm;
		resdualNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;

		if (*ndflag == 0)
		{
			std::cout << "Neumann:: Residual Norm after " << i + 1 << " V-Cycle = " << newnorm << "\n\n";
			std::cout << "Neumann:: Covergence rate after " << i + 1 << " V-Cycle = " << convrate << "\n\n";
		}
		else
		{
			std::cout << "Dirichlet:: Residual Norm after " << i + 1 << " V-Cycle = " << newnorm << "\n\n";
			std::cout << "Dirichlet:: Covergence rate after " << i + 1 << " V-Cycle = " << convrate << "\n\n";
		}

	}
	//orthogonalize(xGrids[0]);
	errorNorm(xGrids[0], sGrid, &newnorm);
	if (*ndflag == 0)
	std::cout << "Neumann:: Error Norm for h as 1/" << gdim - 1 << " = " << newnorm << "\n\n";
	else
	std::cout << "Dirichlet:: Error Norm for h as 1/" << gdim - 1 << " = " << newnorm << "\n\n";
}

int main(int argc, char** argv)
{

	//std::cout << "1";
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	size_t level = atoi(argv[1]);
	size_t vcycle = atoi(argv[2]);

	timeval start, end;
	
	
	std::cout << "\n\n =============== Output for Dirichlet Boundary Value Problem 1 ===================\n\n";
	gettimeofday(&start, 0);
	*ndflag = 1;
	mgsolve(level, vcycle);
	
	gettimeofday(&end, 0);
	double elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
	std::cout << "Dirichlet:: Time spend for Multigrid Solver = " << elapsed << '\n';

	double hsize = (*xGrids[0]).getHx();
	double gdim = (*xGrids[0]).getXsize();

	std::string fname1 = std::string("data/Dirichlet/solution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOut1(fname1);
	std::string fnames1 = std::string("data/Dirichlet/exactsolution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOutsolt1(fnames1);
	for (size_t y = 0; y < gdim; ++y) {
	for (size_t x = 0; x < gdim; ++x) {

	fOut1 << x*hsize << "\t" << y*hsize << "\t" << (*xGrids[0])(x, y) << std::endl;
	fOutsolt1 << x*hsize << "\t" << y*hsize << "\t" << (*sGrid)(x, y) << std::endl;
	}
	fOut1 << std::endl;
	fOutsolt1 << std::endl;
	}
	fOut1.close();
	fOutsolt1.close();
	std::cout << "\n\n =============== Dirichlet Boundary Value Problem 1 ends here ===================\n\n";

	std::cout << "\n\n =============== Output for Neumann Boundary Value Problem 2 ===================\n\n";
	delete xGrids;
	delete fGrids;
	delete sGrid;
	*ndflag = 0;
	gettimeofday(&start, 0);
	mgsolve(level, vcycle);
	gettimeofday(&end, 0);


	elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
	std::cout << "Neumann:: Time spend for Multigrid Solver = " << elapsed << "\n";

	hsize = (*xGrids[0]).getHx();
	gdim = (*xGrids[0]).getXsize();

	std::string fname2 = std::string("data/Neumann/solution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOut2(fname2);
	std::string fnames2 = std::string("data/Neumann/exactsolution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
	std::ofstream	fOutsolt2(fnames2);
	for (size_t y = 0.0; y < gdim; ++y) {
	for (size_t x = 0.0; x < gdim; ++x) {

	fOut2 << x*hsize << "\t" << y*hsize << "\t" << (*xGrids[0])(x, y) << std::endl;
	fOutsolt2 << x*hsize << "\t" << y*hsize << "\t" << (*sGrid)(x, y) << std::endl;
	}
	fOut2 << std::endl;
	fOutsolt2 << std::endl;
	}
	fOut2.close();
	fOutsolt2.close();

	std::cout << "\n\n =============== Neumann Bounday Value Problem 2 ends here ===================\n\n";
	return 0;
}
