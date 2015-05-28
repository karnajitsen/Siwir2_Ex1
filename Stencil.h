struct Stencil
{
private:
	double C, E, W, N, S, NE, NW, SE, SW;
	
public:
	explicit inline Stencil(const double _C, const double _E, const double _W, const double _N,
		const double _S, const double _NE, const double _NW, const double _SE, const double _SW) :C(_C), E(_E), W(_W), S(_S), N(_N),
		NE(_NE), NW(_NW), SE(_SE), SW(_SW){}

	explicit inline Stencil(const double _C, const double _E, const double _W, const double _N,
		const double _S) : C(_C), E(_E), W(_W), S(_S), N(_N){}	
};