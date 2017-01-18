#ifndef CONTROLLERADP_H
#define CONTROLLERADP_H
#include <vector>
#include <list>
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
	template <typename T>
	class ControllerADP: public Controllers
	{
		public: 
			//typedef const std::vector<double> (*noise)(const unsigned int m, const double t);
			ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, Step* stepf=nullptr);
			ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, Step* stepf=nullptr);
			ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const Matrix& K, Step* stepf=nullptr);
			ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, const Matrix& K, Step* stepf=nullptr);
			virtual const std::vector<double> input(const std::vector<double>& x, const double dt, const double t=0, noise noif=&sinusoidal);
			virtual ~ControllerADP(){};
			void dispAll();
		private:
			std::list<std::vector<double>> mxx;
			//std::list<std::vector<double>> mu;
			std::list<std::vector<double>> mIxx;
			std::list<std::vector<double>> mIxu;
			//std::list<std::vector<double>> mTheta;
			//std::list<std::vector<double>> mXi;
			//int mn;
			//int mm;
			Matrix mK0;
			Matrix mKadp;
			SymmetricMatrix mQ;
			SymmetricMatrix mR;
			SymmetricMatrix mP;
			double mdelta;
			std::list<std::vector<double>>::iterator itx;
			std::list<std::vector<double>>::iterator itxx;
			std::list<std::vector<double>>::iterator itxu;
			Matrix mThetaInv;
			std::vector<double> mBigV;
			std::shared_ptr<Matrix> mBigTheta;
			std::shared_ptr<Matrix> mBigr;
			std::shared_ptr<AlgorithmADP> mADPalg;
			void LS(const std::vector<double>& phi, const double d);
			//Matrix mxT;
			//Matrix muT;
			//Matrix mIxT;
			//Matrix mTheta;
	};
#include "ControllerADP.cpp"
}

#endif

