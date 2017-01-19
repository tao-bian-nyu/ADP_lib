#ifndef DYNAMICAL_H
#define DYNAMICAL_H
#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>
//#include "Matrix.h"
#include "SquareMatrix.h"
#include "Controllers.h"
#include "Others.h"

namespace ADP
{

	class Dynamical
	{
		public:
			typedef const std::vector<double> (*inputfun) (const std::vector<double>&, const int, const double);
			typedef const std::vector<double> (Dynamical::*simulate) (const double, const std::vector<double>&, const double);
			Dynamical(const SquareMatrix A, const Matrix B, inputfun input = &linInput, const long double dt=1e-5);
			//Dynamical(const SquareMatrix A, const Matrix B, const std::vector<double> x0, Controllers* controller, const long double dt=1e-5);
			Dynamical(const SquareMatrix A, const Matrix B = Matrix(), Controllers* controller = nullptr, const long double dt=1e-5);
			//const std::vector<double>& x(const double t, simulate = &EU);
			template <typename E>
			const std::vector<double>& x(const double t, const std::vector<double> x0); 
			//const std::vector<double>& x(const double t, const std::vector<double> x0, simulate simFun= &Dynamical::EU); 
			void printAll() const;
			//const std::vector<double> EU(const double t0, const std::vector<double>& x, const double t);
		private:
			SquareMatrix mA;
			Matrix mB;
			std::unordered_map<int,std::vector<double>> mStateAll;
			std::unordered_map<int,std::vector<double>> mInputAll;
			std::vector<double> mx;
			std::vector<double> mu;
			//std::shared_ptr<Controllers> mCtrl;                  // pointer to controller object
			Controllers* mCtrl;                  // pointer to controller object
			inputfun mInFun;                  // pointer to controller object
			double mT;
			double mdt;
			Dynamical(const Dynamical&);     
			const Dynamical& operator=(const Dynamical&);
			template <typename E>
			const std::vector<double> Run(const double t0, const std::vector<double>& x, const double t);
			//const std::vector<double> RK(const double t0, const std::vector<double>& x, const double t);
			//const std::vector<double> EU(const std::vector<double>&, const double t){return mx;};
	};
#include "Dynamical.cpp"

}
#endif
