inline void neumannsmooth(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize();
	size_t dimY = (*xgrd).getYsize();
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	double	alpha = 1.0;
	double	beta = 1.0;
	double	center = 1.0/(2.0 * alpha + 2.0 * beta);

	cout << "====B4 smoothing=== \n\n";
	for (size_t j = 0; j < dimY; j++)
	{
		for (size_t k = 0; k < dimY; k++)
		{
			cout << (*xgrd)(k, j) << " ";
		}

		cout << '\n';
	}

	for (size_t i = 0; i < iter; i++)
	{
		for (size_t j = 1; j < dimY - 1; j++)
		{
			size_t k = ((j + 1) & 0x1);
			if (k == 0)
			{
				(*xgrd)(k, j) = center*(hx*hy*(*fgrd)(k, j) + 2.0*hx + 2.0 * alpha * (*xgrd)(k + 1, j) + beta*(*xgrd)(k , j+1) + (*xgrd)(k , j-1));
				(*xgrd)(dimX - 1, j) = center*(hx*hy*(*fgrd)(k, j) - 2.0*hx + 2.0 * alpha * (*xgrd)(dimX - 2, j)
										+ beta * ((*xgrd)(dimX - 2, j + 1) + (*xgrd)(dimX - 2, j - 1)));
				k += 2;
			}
			for (; k < dimX - 1; k += 2)
			{
					(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;
			}

		}


		for (size_t j = 1; j < dimY - 1; j++)
		{
			size_t k = j & 0x1;
			if (k == 0)
			{
				(*xgrd)(k, j) = center*(hx*hy*(*fgrd)(k, j) + 2.0*hx + 2.0 * alpha * (*xgrd)(k + 1, j) + beta*(*xgrd)(k, j + 1) + (*xgrd)(k, j - 1));
				(*xgrd)(dimX - 1, j) = center*(hx*hy*(*fgrd)(k, j) - 2.0*hx + 2.0 * alpha * (*xgrd)(dimX - 2, j)
					+ beta * ((*xgrd)(dimX - 2, j + 1) + (*xgrd)(dimX - 2, j - 1)));
				k += 2;
			}
			for (;k < dimX - 1; k += 2)
			{
				(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + alpha * ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + beta * ((*xgrd)(k, j + 1)
					+ (*xgrd)(k, j - 1))) * center;
			}

		}

	}


	cout << "====After smoothing=== \n\n";
	for (size_t j = 0; j < dimY; j++)
	{
		for (size_t k = 0; k < dimY; k++)
		{
			cout << (*xgrd)(k, j) << " ";
		}

		cout << '\n';
	}

}


inline void orthogonalize(Grid* grd)
{
	size_t dimX = (*grd).getXsize();
	size_t dimY = (*grd).getYsize();
	double sum = 0.0;

	for (size_t y = 1; y < dimY - 1; y++)
	{
		for (size_t x = 0; x < dimX; x++)
		{
			sum += (*grd)(x, y);
		}
	}

	if (sum != 0.0)
	{
		for (size_t y = 1; y < dimY - 1; y++)
		{
			for (size_t x = 0; x < dimX; x++)
			{
				(*grd)(x, y) -= sum / dimX / (dimY-2.0);
			}
		}
	}

}

inline void MGNeumann(size_t level, size_t vcycle)
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

		for (size_t jl = 0; jl < level - 1; jl++)
		{
			orthogonalize(fGrids[jl]);
			neumannsmooth(xGrids[jl], fGrids[jl], V1);
			restriction(xGrids[jl], fGrids[jl], fGrids[jl + 1]);
		}

		for (size_t j = level - 1; j > 0; j--)
		{
			neumannsmooth(xGrids[j], fGrids[j], V2);
			interpolate(xGrids[j], xGrids[j - 1]);
			(*xGrids[j]).reset();
			(*fGrids[j]).reset();
		}

		oldnorm = newnorm;
		resdualNorm(xGrids[0], fGrids[0], &newnorm);
		if (oldnorm != 0.0)
			convrate = newnorm / oldnorm;

		std::cout << "Neumann:: Residual Norm after " << i + 1 << " V-Cycle = " << newnorm << "\n\n";
		std::cout << "Neumann:: Covergence rate after " << i + 1 << " V-Cycle = " << convrate << "\n\n";

		orthogonalize(xGrids[0]);
	}

	errorNorm(xGrids[0], sGrid, &newnorm);
	std::cout << "Neumann:: Error Norm for h as 1/" << gdim - 1 << " = " << newnorm << "\n\n";

}

