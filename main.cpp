#include <iostream>
#include <vector>
#include <memory>
#include <exception>
#include "Matrix.h"
#include "MatrixCalc.h"
#include "SquareMatrix.h"
#include "SymmetricMatrix.h"
#include "Diagonal.h"
#include "AlgorithmADP.h"
#include "AlgorithmRLS.h"
#include "AlgorithmVI.h"
#include "AlgorithmPI.h"
#include "Dynamical.h"
#include "ControllerADP.h"
#include "Step.h"
#include "RK.h"
#include "EU.h"

int main()
{
	using namespace ADP;
	using VI = AlgorithmVI; 
	using PI = AlgorithmPI; 

	std::cout<< "This is the test of ADPlib v1.0"<<std::endl;

	SquareMatrix sysA({-1,0.1,0,2,-1,0.1,1, 2, 1});
	Matrix sysB({0,0,1},1);
	SymmetricMatrix Q({1,0,0,1,0,1});
	SymmetricMatrix P(Q*0+0.01);
	SymmetricMatrix R(1,1.0);
	Step mystep(1,10,1);
	std::vector<double> x0({2,0,0});
	disp(sysA);
	disp(sysB);
	disp(Q);
	disp(R);

	std::cout << norm<'F'>(sysA) << std::endl;
	std::cout << norm<'E'>(sysA) << std::endl;
	std::cout << norm<2>(x0) << std::endl;

	// simulate ADP + dynamical system
	int n = sysA.size()[0];
	int m = sysB.size()[0];
	Matrix K0({0,0.0,0},n);
	ControllerADP<VI> myADP(Q,R,0.1, P, &mystep);
	ControllerADP<PI> myADP2(Q,R,0.1, K0);
	double t = 0;
	double dt = 0.0001;

	Controllers* myCtrl = &myADP;

	Dynamical myADPsys(sysA, sysB, myCtrl, dt);
	//disp(myADPsys.x(50,x0,&Dynamical::RK));
	disp(myADPsys.x(50,x0));

	myCtrl = &myADP2;

	//Dynamical myADPsys2(sysA, sysB, myCtrl, dt);
	//disp(myADPsys2.x<RK>(50,x0));


	 //offline VI
	std::shared_ptr<AlgorithmADP> myalg(new AlgorithmVI(Q,R,P,&mystep));
	std::vector<Matrix> result = myalg->offline(sysA,sysB);

	result[0].disp();
	result[1].disp();
	result[2].disp();

	
	 //offline PI
	std::shared_ptr<AlgorithmADP> myalg2(new AlgorithmPI(Q,R,K0));
	std::vector<Matrix> result2 = myalg2->offline(sysA,sysB);


	std::cout << "test" << std::endl;
	result2[0].disp();
	result2[1].disp();
	result2[2].disp();

	return 0;
}
