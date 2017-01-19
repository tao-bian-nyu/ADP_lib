#ifndef DYNAMICAL_H
#define DYNAMICAL_H
#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>
//#include "Matrix.h"
#include "SquareMatrix.h"
#include "Controllers.h"

namespace ADP
{

	class Dynamical
	{
		public:
			typedef const std::vector<double> (*inputfun) (const std::vector<double>&, const int m, const double t);
			Dynamical(const SquareMatrix A, const Matrix B, const std::vector<double> x0, inputfun input, const long double dt=1e-5);
			Dynamical(const SquareMatrix A, const Matrix B, const std::vector<double> x0, Controllers* controller, const long double dt=1e-5);
			Dynamical(const SquareMatrix A, const std::vector<double> x0, const long double dt=1e-5);
			const std::vector<double>& x(const double t);
			const std::vector<double>& x(const double t, inputfun input); 
			void printAll() const;
		private:
			SquareMatrix mA;
			Matrix mB;
			std::unordered_map<int,std::vector<double>> mStateAll;
			std::unordered_map<int,std::vector<double>> mInputAll;
			std::vector<double> mx;
			std::vector<double> mu;
			//std::shared_ptr<Controllers> mCtrl;                  // pointer to controller object
			Controllers* mCtrl;                  // pointer to controller object
			double mT;
			double mdt;
			Dynamical(const Dynamical&);     
			const Dynamical& operator=(const Dynamical&);
	};

}
#endif
