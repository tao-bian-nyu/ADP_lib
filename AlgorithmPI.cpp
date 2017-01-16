#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "matrix.h"
#include "AlgorithmPI.h"
#include "AlgorithmVI.h"
#include "AlgorithmADP.h"
#include "MatrixCalc.h"

namespace ADP
{
	std::vector<Matrix> AlgorithmPI::offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N, double eps)
	{
		const unsigned int n = sysA.size()[0];
		mP = SymmetricMatrix(n,0);
		std::vector<double> vecP(n*n,0);
		Matrix mKold(mK0);
		//error = P; 
		//Matrix BP(n,n,0.0); 
		for(mk = 0; mk<N; mk++)
		{
			//std::vector<Matrix> result = myalg.offline(sysA-sysB*R*BP, sysB*0, Q+BP.t()*R*BP, R);
			vecP = vec(inv(kSum(t(sysA-sysB*mK),t(sysA-sysB*mK))) * vec(-mQ-t(mK)*mR*mK));
			mP = Matrix(vecP,n);
			mKold = mK;
			mK = inv(mR)*t(sysB)*mP;
			if((mK-mKold).F()<eps){ 
				std::cout << std::endl; 
				mResult = std::vector<Matrix>({mK-mKold, mP, mK}); 
				return  mResult;
			}

			std::cout << "PI trial " << mk << ','; 
			std::cout.flush();

		} 
		std::cout << std::endl; 
		mResult = std::vector<Matrix>({mK-mK0, mP, mK}); 
		return  mResult;
	}

	AlgorithmPI::~AlgorithmPI(){};
}

