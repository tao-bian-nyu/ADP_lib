//#include <algorithm>
#include <iostream>
//#include <vector>
//#include "Matrix.h"
#include "AlgorithmRLS.h"
#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
//#include "Diagonal.h"
#include "MatrixCalc.h"


namespace ADP
{
	AlgorithmRLS::AlgorithmRLS(const Matrix& Phi, const std::vector<double>& d, const long double epsilon)
		//:mW(Phi.size()[1],0),mP(Phi.size()[1],Phi.size()[1],0),mEpsilon(epsilon)
		:mW(Phi.size()[1],0),mP(Phi.size()[1]),mEpsilon(epsilon)
	{
		const unsigned int n = Phi.size()[0];
		const unsigned int m = Phi.size()[1];
		//std::vector<double> g(m,0);
		//Diagonal P0(m,1/mEpsilon);
		//SquareMatrix P(mP+1/mEpsilon); 
		mP = mP+1/mEpsilon; 

		for(unsigned int i=1;i<=n;++i)
		{
			mP = mP - 1 / (1 + double(t(Phi.col(i))*mP*Phi.col(i))) * mP * Phi.col(i) * t(mP * Phi.col(i));
			//g = vec(mP * Phi.col(i));
			double alpha = d[i-1] - mW * Phi.col(i);
			//mW = mW + g * alpha;
			mW = mW + alpha * vec(mP * Phi.col(i));
		}

	}


	const std::vector<double>& AlgorithmRLS::disp() const
	{
		return mW;
	}

	void AlgorithmRLS::print() const
	{
		std::cout << "w is " << std::endl;
		ADP::disp(mW);
		std::cout << "P is " << std::endl;
		ADP::disp(mP);
	}
}


