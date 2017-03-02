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
		: mQ(&Q), mR(&R), mP0(&P0), mP(P0), mK(), mbound(bound), mk(1),mStep(stepf){}


	std::shared_ptr<AlgorithmADP> AlgorithmVI::Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound) const
	{
		return std::shared_ptr<AlgorithmADP>(new AlgorithmVI(Q,R,P0,stepf,bound));
	}


	const std::vector<Matrix>& AlgorithmVI::offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N, double eps)
	{
		SymmetricMatrix error(mP);
		for(unsigned mk = 1; mk<N; mk++)
		{
			const double step = mStep->stepOut(mk);
			error = T(sysA) * mP + mP * sysA + *mQ - mP * sysB * inv(*mR) * T(sysB) * mP;
			mP = mP + step * error;
			if(norm(mP) > mbound || !(mP>0))
			{
				std::cout << "VI trial: " << mk << std::endl; 
				mbound += 100;
				mP=*mP0;
				continue;
			}else if(norm(error)<eps){
				std::cout << "VI trial: " << mk << std::endl; 
				mResult = std::vector<Matrix>({error, mP, inv(*mR)*T(sysB)*mP});
				return  mResult;
			}

		}
		std::cout << "reach maximum loop number" << std::endl;
		mResult = std::vector<Matrix>({error, mP, inv(*mR)*T(sysB)*mP});
		return  mResult;


	}

	void AlgorithmVI::onlineI(const std::vector<double>&  vec)
	{
		//std::cout << "VI loop " << mk << std::endl;
		//ADP::disp(vec);
		double step = 0;
		if (mStep==nullptr)
			step = 1.0/(100 + mk++);
		else
			step = mStep->stepOut(mk++);
		//
		////
		////
		//
		// add new term. if mStep is nullptr, use default step function; otherwise use mStep;
		//
		//
		//
		//
		std::cout << "step is " << step << std::endl;
		const unsigned int n = mQ->size()[0];
		//auto first = vec.begin();
		//const auto last = vec.begin()+n*(n+1)/2;
		const std::vector<double> vecH(vec.begin(),vec.begin()+n*(n+1)/2);
		//first = vec.end();
		const std::vector<double> vecK(vec.begin()+n*(n+1)/2,vec.end());
		mK = Matrix(vecK,n);
		SymmetricMatrix error(SymmetricMatrix(vecH)+*mQ-T(mK)* *mR *mK);
		mP = mP + step * error;
		ADP::disp(error);

		if(norm(mP) > mbound || !(mP>0))
		{
			mbound+=100;
			mP= *mP0;
		}

		mResult = std::vector<Matrix>({error, mP, mK}); 
	}


	void AlgorithmVI::resetStep()
	{
		mk=1;
	}

	const std::vector<Matrix>& AlgorithmVI::online(const std::vector<double>& vec0, const std::vector<double>& vec1, const std::vector<double>& vec2, std::shared_ptr<Matrix> mBigr, SymmetricMatrix& mThetaInv, std::vector<double>& mBigV)
	{
		//ADP::disp(vec1);
		//ADP::disp(vec2);
		//ADP::disp(vec0);
		std::vector<double> phi;
		phi.reserve(vec1.size() + vec2.size());
		phi = vec1;
		phi.insert(phi.end(), vec2.begin(), vec2.end());

		//ADP::disp(phi);
		//*mBigr= *mBigr + prod(phi,vec0);// half online case
		mThetaInv = mThetaInv - 1 / (1 + double(T(phi)*mThetaInv*phi)) * mThetaInv * phi * T(mThetaInv * phi);
		//mThetaInv =  1/0.9* ( mThetaInv - 1 / (0.9 + double(T(phi)*mThetaInv*phi)) * mThetaInv * phi * T(mThetaInv * phi));   // discounted
		*mBigr= *mBigr + mThetaInv*phi*T(vec0) - mThetaInv*phi*T(phi)* *mBigr;  // online case;
		mBigV = vec(*mBigr * vec(mP));
		//mBigV = mBigV + (vec0*vec(mP) - mBigV * phi)* vec(mThetaInv * phi);// half online case

		//ADP::disp(mBigV);

		onlineI(mBigV);
		return mResult;
	}	

	void AlgorithmVI::disp() const
	{
		mResult[0].disp();
		mResult[1].disp();
		mResult[2].disp();

	}
}
