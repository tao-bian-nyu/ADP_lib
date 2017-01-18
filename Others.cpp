#include <vector>
#include <math.h>
//#include <iostream>
#include "Controllers.h"
#include "Matrix.h"
#include "MatrixCalc.h"
#include "AlgorithmADP.h"
#include "Step.h"
#include "Others.h"
//#include "AlgorithmVI.h"
//#include "AlgorithmPI.h"

namespace ADP{

	const std::vector<double> sinusoidal(const unsigned int m, const double t)
	{

		std::vector<double>out(m,0);
		for (unsigned int k=0;k<m;++k)
			for(int j=1;j<=500;++j)
				out[k] += 0.05*sin(j*t+2*M_PI*k/m);
		return out; 
	}

	const std::vector<double> linInput(const std::vector<double>& x, const int m, const double t){
		Matrix K(x.size(),m,1);
		std::vector<double> u = vec(K * x);
		return u;
	}
}
