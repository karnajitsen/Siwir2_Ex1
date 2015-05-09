#pragma once
#include<iostream>
#include <assert.h>
#define LD 16
class Grid
{

	double * __restrict data = NULL;
	size_t sizeX, sizeY, ld, totLength;
	double hx, hy;
	
public:
	explicit Grid()
	{
		data = new double[0]();
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
		data = new double[ld*y]();
		data++;
	}
	~Grid()
	{
		--data;
		free(data);
	}

	inline double& operator()(const int x, const int y)
	{
		assert((unsigned) x < sizeX);
		assert((unsigned) y < sizeY);
		return data[y*ld + x];
	}

	inline double& operator()(const int x, const int y) const
	{
		assert((unsigned)x < sizeX);
		assert((unsigned)y < sizeY);
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

