#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <memory>
#include <malloc.h>
#define LD 16
#define ALLIGNMENT 1024
//#define M_PI 3.14
class Grid
{

	//__declspec(align(128)) 
	double * __restrict data = NULL;
	size_t sizeX, sizeY, ld, totLength;
	double hx, hy;
	
public:
	explicit Grid
	{
		data = (double*)_aligned_malloc(0, ALLIGNMENT);
		sizeX = 0;
		sizeY = 0;
	}

	explicit Grid(const size_t x, const size_t y, const double& _hx, const double& _hy )
	{
		sizeX = x;
		sizeY = y;
		hx = _hx;
		hy = _hy;
		ld = x + LD;
		totLength = (x - 2)*(y - 2);
		data = (double*) _aligned_malloc(ld*y*sizeof(double), ALLIGNMENT);

		for (int j = 0.0; j < sizeX; j++)
		{
			double k = j*hx;
			double l = (sizeX - 1.0)*hx;
			data[j] = gxy(k, 0.0);
			data[j*ld] = gxy(0.0, k);
			data[j + ld * (sizeX - 1)] = gxy(k, l);
			data[j * ld + (sizeX - 1)] = gxy(l, k);
		}
		//data++;
	}
	~Grid()
	{
		//--data;
		free(data);
	}

	inline double gxy(const double x, const double y)
	{
		return sin(M_PI * x) * sinh(M_PI * y);
	}


	inline double& operator()(const int x, const int y)
	{
		assert((unsigned int)x < sizeX);
		assert((unsigned int)y < sizeY);
		return data[y*ld + x];
	}

	inline double& operator()(const int x, const int y) const
	{
		assert((unsigned int)x < sizeX);
		assert((unsigned int)y < sizeY);
		return data[y*ld + x];
	}

	inline Grid * operator+=(const Grid * rhs)
	{
		return this;
	}

	inline Grid* operator-=(const Grid* rhs)
	{

	}
		
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

	inline size_t getHx() const
	{
		return hx;
	}

	inline size_t getHy() const
	{
		return hy;
	}

};

