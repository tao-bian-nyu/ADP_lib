#include <iostream>
#include <vector>
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
#include <memory>

int main()
{
	using namespace ADP;
	std::cout<< "This is the test of ADPlib v1.0"<<std::endl;
	int ncol=5;
	int nrow=3;
	Matrix mats(ncol,nrow);
	SquareMatrix squmat(ncol);

	std::vector<double> vect01({2,3.1});

	Matrix sysB0v(vect01);
	sysB0v.disp();

	Matrix wd = t(sysB0v) * sysB0v;
	std::cout<< "This is wd"<<std::endl;
	wd.disp();
	mats.disp();


	std::vector<double> vect(5,2.3);

	mats(1,1) = 112.0;
	Matrix newmat =  mats.t() * mats;
	newmat.disp();
	Matrix mm = mats * vect;
	mm.disp();
	std::cout << sizeof(double) << std::endl;

	double x = mats(1,2)+1;
	x++;
	Matrix mat2(ncol,nrow);
	Matrix mat3(mats);
	mat2=mats+mat3;

	std::vector<double> vec(9,2.3);
	mat2.col(1,vec);
	Matrix mat4 = prod(vec,vec);
	mat4.disp();

	std::cout <<"ithis " <<mat2.col(1)[2] << ',' << std::endl;
	std::cout << mat2.size()[0] << ',' << mat2.size()[1] << std::endl;
	SquareMatrix addMat(4);
	addMat.add(1.0);
	Matrix temp(3.0 + addMat);

	Matrix madmat(4,4,0.0);

	//SquareMatrix* poi = (SquareMatrix*)&madmat;
	SquareMatrix* sMat3 = dynamic_cast<SquareMatrix*>(&madmat);
	//SquareMatrix sMat2 = *poi;
	//matDisplay(*sMat3);
	if(sMat3)
		sMat3->disp();

	temp.disp();
	kProd(temp,temp).disp();
	//matDisplay(kProd(temp,temp));

	vec[2]=1.0;
	Diagonal sMat(vec);
	// SquareMatrix sMat(mat2);
	std::cout << sMat.size()[0] << ',' << sMat.size()[1] << std::endl;

	std::cout << std::endl;
	sMat.disp();


	Matrix mat5(3,10, 100);
	mat2 = mat5.t();
	std::cout << mat2.size()[0] << ',' << mat2.size()[1] << std::endl;
	mat2.disp();


	//my RLS algorithm
	SquareMatrix phi({-1,0,2,-1});
	std::vector<double> d({1,1});
	AlgorithmRLS myRLS(phi, d, 0.0001);

	myRLS.print();
	phi.inv().disp();

	// std::cout << sMat.size()[0] << ',' << sMat.size()[1] << std::endl;
	// std::cout << sMat(1,1) << ','  << std::endl;

	SquareMatrix sysA({-1,-0.1,2,-1});
	Matrix sysB({0,1},1);
	SymmetricMatrix Q({1,0,1});
	sysA.disp();
	sysB.disp();
	Q.disp();
	//SymmetricMatrix P({0.01,0,0.01});
	SymmetricMatrix P({0.0197,0.000188,0.0197});
	P.disp();
	SquareMatrix R(1,1.0);
	R.disp();
	std::vector<double> x0(2,2);

	//Dynamical mysys(sysA, x0, 10);
	//sysA.disp();
	//Matrix(mysys.disp(10)).disp();
	//Matrix K0(t(sysB));
	Matrix K0({0,0.1},2);
	//mysys.printAll();


	ControllerADP myADP(sysA.size()[0], sysB.size()[0], Q, R, P, 0.1, K0);
	double t = 0;
	double dt = 0.0001;

	Controllers* myCtrl = &myADP;

	Dynamical myADPsys(sysA, sysB, x0, 100, myCtrl, dt);

	//for(int j=0;j<10;++j)
	//{
		//x0 = x0;
		//t += dt;
		//myADP.input(x0, dt,t);
	//}

	//myADP.dispAll();

	std::shared_ptr<AlgorithmADP> myalg(new AlgorithmVI(P));
	std::vector<Matrix> result = myalg->offline(sysA,sysB,Q,R);

	result[0].disp();
	result[1].disp();
	result[2].disp();

	result = myalg->offline(SquareMatrix({-1,-0.0,3,-1}),sysB,Q,R);

	result[0].disp();
	result[1].disp();
	result[2].disp();

	//std::shared_ptr<AlgorithmADP> myalg2(new AlgorithmPI);
	//std::vector<Matrix> result2 = myalg2->offline(sysA,sysB,Q,R);


	//result2[0].disp();
	//result2[1].disp();
	//result2[2].disp();

	return 0;
}
