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
#include "EU.h"
#include "RK.h"

namespace ADP
{

	class Dynamical
	{
		public:
			typedef const std::vector<double> (*inputfun) (const std::vector<double>&, const int, const double);
			//typedef const std::vector<double> (Dynamical::*simulate) (const double, const std::vector<double>&, const double);
			Dynamical(const SquareMatrix& A, const Matrix& B, inputfun input, const long double dt=1e-5);
			Dynamical(const SquareMatrix& A, const Matrix& B = Matrix(), Controllers* controller = nullptr, const long double dt=1e-5);
			template <typename E=EU>
			const std::vector<double>& x(const double t, const std::vector<double> x0); 
			void printAll() const;
			std::unordered_map<int,std::vector<double>>::const_iterator expState() const;
		private:
			const SquareMatrix* mA;
			const Matrix* mB;
			std::unordered_map<int,std::vector<double>> mStateAll;                  // all the state are in a Hashtable
			std::unordered_map<int,std::vector<double>> mInputAll;                  // all the input are in a Hashtable
			std::vector<double> mx;                                        // current state
			std::vector<double> mu;                                     // current input
			Controllers* mCtrl;                  // pointer to controller object
			inputfun mInFun;                     // pointer to controller function 
			double mT;                            // current time
			double mdt;                      // simulation time
			Dynamical(const Dynamical&);    // a dynamical system is unique, cannot be copied  
			const Dynamical& operator=(const Dynamical&);   // a dynamical system is unique, cannot be assigned 
			template <typename E>
			const std::vector<double> Run(const double t0, const std::vector<double>& x, const double t);
	};
#include "Dynamical.cpp"

}
#endif
