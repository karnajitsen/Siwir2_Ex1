#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <malloc.h>
#include <string>
#define LD 16
#define ALLIGNMENT 32
using namespace std;
//#define M_PI 3.14
class Grid
{

	//__declspec(align(128)) 
	double * __restrict data = NULL;
	size_t sizeX, sizeY, ld, totLength;
	double hx, hy;
	
public:
	explicit Grid()
	{
		data = (double*)memalign(ALLIGNMENT, 0);
		sizeX = 0;
		sizeY = 0;
	}

	explicit Grid(const size_t x, const size_t y, const double& _hx, const double& _hy , bool bndrYN)
	{
		sizeX = x;
		sizeY = y;
		hx = _hx;
		hy = _hy;
		ld = x + LD;
		totLength = (x - 2)*(y - 2);
		data = (double*) memalign(ALLIGNMENT, ld*y*sizeof(double));
		//data = (double*) _aligned_malloc(ld*y*sizeof(double), ALLIGNMENT);
		if (bndrYN)
		{
			double l = (sizeX - 1.0)*hx;
			for (int j = 0.0; (size_t)j < sizeX; j++)
			{
				double k = j*hx;

				data[j] = gxy(k, 0.0);
				data[j*ld] = gxy(0.0, k);
				data[j + ld * (sizeX - 1)] = gxy(k, l);
				data[j * ld + (sizeX - 1)] = gxy(l, k);
			}
		}
		//data++;
	}
	~Grid()
	{
		//--data;
		free(data);
	}

	inline long double gxy(const double x, const double y)
	{
		return (long double) sin(M_PI * x) * sinh(M_PI * y);
	}

	inline void reset()
	{
		for (size_t i = 1; i < sizeX - 1; i++)
		{
			for (size_t j = 1; j < sizeX - 1; j++)
			{
				data[i*ld + j] = 0.0;
			}
		}
	}

	inline void print(std::string text)
	{
		cout << '\n\n';
		cout << " ======== " << text << " ========" << '\n\n';
		for (size_t j = sizeX - 1; j >= 0; j--)
		{
			for (size_t k = 0; k < sizeX; k++)
			{
				cout << data[j*ld + k] << " ";
			}

			cout << '\n';
		}
	}

	inline double& operator()(const size_t x, const size_t y)
	{
		assert(x < sizeX);
		assert(y < sizeY);
		return data[y*ld + x];
	}

	inline double& operator()(const size_t x, const size_t y) const
	{
		assert(x < sizeX);
		assert(y < sizeY);
		return data[y*ld + x];
	}

	/*inline Grid * operator+=(const Grid * rhs)
	{
		return this;
	}

	inline Grid* operator-=(const Grid* rhs)
	{
		return this;
	}
	*/	
	inline size_t getXsize() const
	{
		return sizeX;
	}

	inline size_t getYsize() const
	{
		return sizeY;
	}

	inline size_t getTotLength() const
	{
		return totLength;
	}

	inline double getHx() const
	{
		return hx;
	}

	inline double getHy() const
	{
		return hy;
	}

};

