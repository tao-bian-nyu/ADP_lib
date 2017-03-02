#include <iostream>
#include <vector>
#include <memory>
#include <exception>
#include <iterator> 
#include <algorithm>
#include "../Matrix.h"
#include "../MatrixCalc.h"
#include "../SquareMatrix.h"
#include "../SymmetricMatrix.h"
#include "../Diagonal.h"
#include "../AlgorithmADP.h"
#include "../AlgorithmRLS.h"
#include "../AlgorithmVI.h"
#include "../AlgorithmPI.h"
#include "../Dynamical.h"
#include "../ControllerADP.h"
#include "../Step.h"
#include "../RK.h"
#include "../EU.h"
#include "../ControllerADP.h"

namespace ADP{
	extern "C" {
		typedef   ControllerADP<AlgorithmVI> ControllerVI;
		//ControllerADP<AlgorithmVI>* ControllerVI_new(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, Step* stepf=nullptr){
		//return new ControllerADP<AlgorithmVI>(Q, R, delta, stepf);
		//}
		ControllerVI* ControllerVI_new(){
			std::cout << "this is a test" << std::endl;
			SymmetricMatrix Q({1,0,0,1,0,1});
			SymmetricMatrix P(Q*0+0.01);
			SymmetricMatrix R(1,1);
			//std::shared_ptr<Step> mystep(new Step(1,10,1));
			//Step* mystep = new Step(1,10,1);
			return new ControllerVI(Q, R, 0.1, P);
			//return new ControllerVI(Q, R, 0.1, P, mystep);
		}
		//ControllerADP<AlgorithmVI>* ControllerVI_new(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const Matrix& K, Step* stepf=nullptr){
		//return new ControllerADP<AlgorithmVI>(Q, R, delta, K, stepf);
		//}
		//ControllerADP<AlgorithmVI>* ControllerVI_new(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, const Matrix& K, Step* stepf=nullptr){
		//return new ControllerADP<AlgorithmVI>(Q, R, delta, P, K, stepf);
		//}

		//const std::vector<double> ControllerVI_input(ControllerVI* self, const std::vector<double>& x, const double dt, const double t=0, Controllers::noise noif=&sinusoidal){
		//return self->input(x, dt, t, noif); 
		//}


		//const Matrix ControllerVI_learner(ControllerVI* self, const std::vector<double>& x, const std::vector<double>& u, const double dt, const double t=0){
		//return self->learner(x, u, dt, t);
		//}


		// This is the real one


		//std::vector LuaScript::getVector(const std::string& name) {
		//std::vector<double> v;
		//lua_getglobal(L, name.c_str());
		//if(lua_isnil(L, -1)) {
		//return std::vector();
		//}
		//lua_pushnil(L);
		//while(lua_next(L, -2)) {
		//v.push_back((double)lua_tonumber(L, -1));
		//lua_pop(L, 1);
		//}

		//double n = lua_gettop(L);
		//lua_pop(L, n);
		//return v;
		//}

		//ControllerVI* ControllerVI_new(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, Step* stepf=nullptr){
		//return new ControllerVI(Q, R, delta, P, stepf);
		//}
double* ControllerVI_learner(ControllerVI* self, const double* x, const double* u, const int n, const int m, const double dt, const double t=0){ std::vector<double> xVec(x, x + n); std::vector<double> uVec(u, u + m);
			//disp(xVec);
			//disp(uVec);
			//std::cout << dt << ',' << t << std::endl;
			//std::vector<double> xVec({1,2,3});
			//std::vector<double> uVec({1});
			//std::cout << uVec.size()  << std::endl;
			//disp(xVec);
			Matrix K = self->learner(xVec, uVec, dt, t);
			//disp(K);
			std::vector<double> vecK = vec(K);
			//disp(vecK)
			double* vecIt = new double[vecK.size()];
			std::copy(vecK.begin(), vecK.end(), vecIt);
			//std::cout << vecIt << std::endl;
			return vecIt;
		}

		void ControllerVI_dispAll(ControllerVI* self){
			self->dispAll();
		}	
		void ControllerVI_delete(ControllerVI* self){

			delete self; 
		}

	}
}
