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

	const std::vector<double> sinusoidal(const double magnitude, const unsigned int m, const double t)
	{

		std::vector<double>out(m,0);
		//for (unsigned int k=0;k<m;++k)
		unsigned int k= 0;
		for (auto it=out.begin(); it!=out.end(); ++it)
			for(unsigned int j=1;j<=500;++j)
				*it += sin(j*t+2*M_PI*k++/m);
				//out[k] += sin(j*t+2*M_PI*k/m);
		return magnitude*out; 
	}

	const std::vector<double> linInput(const std::vector<double>& x, const int m, const double t){
		Matrix K(x.size(),m,1);
		std::vector<double> u = vec(K * x);
		return u;
	}
}
