//#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
//#include "matrix.h"
#include "AlgorithmVI.h"
//#include "AlgorithmADP.h"

namespace ADP
{


	AlgorithmVI::AlgorithmVI(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, Step* stepf, const double bound)
		: mQ(Q), mR(R), mP0(P0), mP(P0), mbound(bound), mk(1),mStep(stepf){}


	std::vector<Matrix> AlgorithmVI::offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N, double eps)
	{
		//const int n = sysA.size()[0];
		//Matrix P(n,n,0.0); 
		//SymmetricMatrix P(mP0); 
		SymmetricMatrix error(mP);
		for(mk = 1; mk<N; mk++)
		{
			//const double step = 1/(mk+100.0);
			const double step = mStep->stepOut(mk);
			error = t(sysA) * mP + mP * sysA + mQ - mP * sysB * inv(mR) * t(sysB) * mP;
			//Matrix tempP = P + step * error;
			mP = mP + step * error;
			if(mP.F() > mbound || !(mP.F()>0))
			{
				mbound += 100;
				mP=mP0;
				continue;
			}else if(error.F()<eps){
				std::cout << "VI trial: " << mk << std::endl; 
				mResult = std::vector<Matrix>({error, mP, inv(mR)*t(sysB)*mP});
				//std::vector<Matrix> vecOut({error, mP, inv(R)*t(sysB)*mP});
				//mResult = vecOut; 
				return  mResult;
			}
			//else P = mP;
			//P.disp();

		}
		std::cout << "reach maximum loop number" << std::endl;
		mResult = std::vector<Matrix>({error, mP, inv(mR)*t(sysB)*mP});
		//std::vector<Matrix> vecOut({error, mP, inv(R)*sysB.t()*mP});
		return  mResult;


	}

	std::vector<Matrix> AlgorithmVI::onlineI(const std::vector<double>&  vec)
	{
		std::cout << "VI loop " << mk << std::endl;
		//const double step = 1/(mk+100.0);
		//const double step = Step(mk);
		const double step = mStep->stepOut(mk++);
		const unsigned int n = mQ.size()[0];
		//const int m = R.size()[0];
		auto first = vec.begin();
		const auto last = vec.begin()+n*(n+1)/2;
		const std::vector<double> vecH(first,last);
		//first = last;
		first = vec.end();
		const std::vector<double> vecK(last,first);
		//SymmetricMatrix H(vecH);
		//H.disp();
		const Matrix K(vecK,n);
		//K.disp();
		//mResult[1] = mResult[1] + eps*(H+Q-t(K)*R*K);
		SymmetricMatrix error(SymmetricMatrix(vecH)+mQ-t(K)*mR*K);
		//std::cout << step << std::endl;
		mP = mP + step * error;
		//mP = mP + step * error;

		if(mP.F() > mbound || !(mP.F()>0))
		{
			mbound+=100;
			mP=mP0;
		}

		//std::vector<Matrix> vecOut({error, mP, K});
		//std::vector<Matrix> mResult({error, P, K});
		mResult = std::vector<Matrix>({error, mP, K}); 
		return mResult;
	}


	std::vector<Matrix> AlgorithmVI::onlineB(const std::vector<double>&  state, const std::vector<double>& input)
	{
		return mResult;
	}


	void AlgorithmVI::resetStep()
	{
		mk=1;
	}

	//void AlgorithmVI::offline(const SquareMatrix& sysA, const Matrix& sysB, const SquareMatrix& Q, const SquareMatrix& R, double eps)
	//{
	//mResult = offline(const SquareMatrix& sysA, const Matrix& sysB, const SquareMatrix& Q, const SquareMatrix& R, double eps);
	//}

}
