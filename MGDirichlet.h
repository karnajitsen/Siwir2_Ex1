void init(double, const size_t, bool);
void restriction(const Grid *, const Grid *, Grid*);
void interpolate(Grid *, Grid *);
void resdualNorm(const Grid*, const Grid *, double*);
void errorNorm(const Grid*, const Grid *, double*);

using namespace std;

inline void dirichsmooth(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize();
	size_t dimY = (*xgrd).getYsize();
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
	double	center = 1.0/(2.0 * alpha + 2.0 * beta);

	for (size_t i = 0; i < iter; i++)
	{
		for (size_t j = 1; j < dimY - 1; j++)
		{
			for (size_t k = ((j + 1) & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;

			}

		}

		for (size_t j = 1; j < dimY - 1; j++)
		{
			for (size_t k = (j & 0x1) + 1; k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;


			}

		}
	}
}


inline void mgDirichlet(size_t level, size_t vcycle)
{
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	init(hsize, level, true);
	/*size_t gdim = pow(2, level) + 1;
	double oldnorm = 0.0, newnorm = 0.0, convrate = 0.0;
	double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	std::cout << "00000";
	init(hsize, level, true);
	std::cout << "1111131131131";
	sGrid = new Grid(gdim, gdim, hsize, hsize, true, true);

	for (size_t i = 0; i < gdim; i++)
	{
		for (size_t j = 0; j < gdim; j++)
		{
			(*sGrid)(j, i) = (*sGrid).gxy1(j*hsize, i*hsize);
		}
	}
	std::cout << "11111";

	for (size_t i = 0; i < vcycle; i++)
	{

		for (size_t jl = 0; jl < level - 1; jl++)
		{
			dirichsmooth(xGrids[jl], fGrids[jl], V1);
			cout << "2222";
			restriction(xGrids[jl], fGrids[jl], fGrids[jl + 1]);
		}

		for (size_t j = level - 1; j > 0; j--)
		{
			dirichsmooth(xGrids[j], fGrids[j], V2);
			interpolate(xGrids[j], xGrids[j - 1]);
			(*xGrids[j]).reset();
			(*fGrids[j]).reset();
		}

		oldnorm = newnorm;
		resdualNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;
		
		std::cout << "Dirichlet:: Residual Norm after " << i + 1 << " V-Cycle = " << newnorm << "\n\n";
		std::cout << "Dirichlet:: Covergence rate after " << i + 1 << " V-Cycle = " << convrate << "\n\n";
	}
	errorNorm(xGrids[0], sGrid, &newnorm);
	std::cout << "Dirichlet:: Error Norm for h as 1/" << gdim - 1 << " = " << newnorm << "\n\n";*/

}

