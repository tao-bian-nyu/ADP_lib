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
		:mW(Phi.size()[1],0),mP(Phi.size()[1]),mEpsilon(epsilon)
	{
		const unsigned int n = Phi.size()[0];
		const unsigned int m = Phi.size()[1];
		mP = mP+1/mEpsilon; 

		for(unsigned int i=1;i<=n;++i)
		{
			mP = mP - 1 / (1 + double(T(Phi.col(i))*mP*Phi.col(i))) * mP * Phi.col(i) * T(mP * Phi.col(i));
			mW = mW + (d[i-1] - mW * Phi.col(i)) * vec(mP * Phi.col(i));
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


