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


	AlgorithmPI::AlgorithmPI(const SymmetricMatrix& Q, const SymmetricMatrix& R, const Matrix& K0)
		: mQ(Q), mR(R), mP(Q), mK0(K0), mK(K0), mk(1), mResult({mK-mK0, Q, mK}){}
	std::shared_ptr<AlgorithmADP> AlgorithmPI::Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound)
	{
		return std::shared_ptr<AlgorithmADP>(new AlgorithmPI(Q,R,K0));
	}

	std::vector<Matrix> AlgorithmPI::onlineI(const std::vector<double>&  vec)
	{
		//std::cout << "PI loop " << mk++ << std::endl;
		//const double step = 1/(mk+100.0);
		//const double step = Step(mk);
		const unsigned int n = mQ.size()[0];
		//const int m = R.size()[0];
		auto first = vec.begin();
		const auto last = vec.begin()+n*(n+1)/2;
		const std::vector<double> vecP(first,last);
		//first = last;
		first = vec.end();
		const std::vector<double> vecK(last,first);
		//SymmetricMatrix H(vecH);
		//H.disp();
		const Matrix K(vecK,n);
		//const Matrix err(mK-Matrix(vecK,n));
		//mResult[1] = mResult[1] + eps*(H+Q-t(K)*R*K);
		mP=SymmetricMatrix(vecP);
		//std::cout << step << std::endl;
		//mP = mP + step * error;

		//std::vector<Matrix> vecOut({error, mP, K});
		//std::vector<Matrix> mResult({error, P, K});
		mResult = std::vector<Matrix>({K-mK, mP, K}); 
		//mResult[0] = K-mK; 
		mK = K;
		return mResult;

	}



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
			disp(mP);
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
		mResult = std::vector<Matrix>({mK-mKold, mP, mK}); 
		return  mResult;
	}


	std::vector<Matrix> AlgorithmPI::online(const std::vector<double>& vec0, const std::vector<double>& vec1, const std::vector<double>& vec2, std::shared_ptr<Matrix> mBigr, SymmetricMatrix& mThetaInv, std::vector<double>& mBigV)
	{
		std::vector<double> phi;
		phi.reserve(vec1.size() + vec0.size());
		phi = vec0;
		//std::vector<double> vec2 = vec(2*kProd(Diagonal(mQ.size()[0]),mR*mKadp)*(mIxx.back() - *itxx)) + mIxu.back() - *itxu;
		std::vector<double> vec3 = vec(2*kProd(Diagonal(mQ.size()[0]),mR*mK)*vec1) + vec2;
		phi.insert(phi.end(), vec3.begin(), vec3.end());

		mThetaInv = mThetaInv - 1 / (1 + double(t(phi)*mThetaInv*phi)) * mThetaInv * phi * t(mThetaInv * phi);
		std::vector<double> err((vec1*vec(mQ+mK.t()*mR*mK) - mBigV * phi)* (vec(mThetaInv * phi)));
		mBigV = mBigV + err;

		if((err*err)<1e-10)
		{
			onlineI(mBigV);
			//disp(mResult[1]);
			long double eps = 1e-10;
			mThetaInv=mThetaInv*0+1/eps; 
			mBigV = mBigV * 0;

			return mResult;
		}
		return mResult;
	}


	AlgorithmPI::~AlgorithmPI(){};
}

