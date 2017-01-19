//#include "Matrix.h"
//#include "Dynamical.h"
//#include <iostream>
//#include <math.h>
//#include "MatrixCalc.h"
//#include "Others.h"

//namespace ADP
//{

	Dynamical::Dynamical(const SquareMatrix A, const Matrix B, inputfun input, const long double dt)
:mA(A),mB(B),mT(0),mdt(dt),mCtrl(nullptr), mInFun(input)
	//,mStateAll(floor(T/dt),A.size()[0],0)
{
	//mStateAll.insert({0,mx});
	//mInputAll.insert({0,mu});
	if (mA.size()[0]!=mB.size()[1])
	{
		std::cout << "unmatched system matrices" << std::endl;
		throw;
	}
	//double t = 0;
	//mStateAll.col(1,mx);
	//for(int i=0;i<=floor(T/mdt);++i){
	//t += mdt;
	//mu = input(mx,B.size()[0],t);
	//mx = mx + mdt*vec(mA * mx + mB * mu);
	//mStateAll.col(i+1,mx);
	//}
}

//Dynamical::Dynamical(const SquareMatrix A, const std::vector<double> x0,  const long double dt)
//:mA(A),mB(1,A.size()[0],0),mx(x0),mu(1,0),mT(0),mdt(dt),mCtrl(nullptr)
////,mStateAll(floor(T/dt),A.size()[0],0)
//{
//mStateAll.insert({0,mx});
//mInputAll.insert({0,mu});
//if (mA.size()[0]!=mB.size()[1])
//{
//std::cout << "unmatched system matrices" << std::endl;
//throw;
//}
////mStateAll.col(1,mx);
////for(int i=0;i<=floor(T/mdt);++i){
////mx = mx + mdt*vec(mA * mx);
////mStateAll.col(i+1,mx);
////}
//}

	Dynamical::Dynamical(const SquareMatrix A, const Matrix B, Controllers* controller, const long double dt)
:mA(A),mB(B),mT(0),mdt(dt),mCtrl(controller), mInFun(nullptr)
	//,mStateAll(floor(T/dt),A.size()[0],0)
{
	//mStateAll.insert({0,mx});
	//mInputAll.insert({0,mu});
	if (mA.size()[0]!=mB.size()[1])
	{
		std::cout << "unmatched system matrices" << std::endl;
		throw;
	}
	//double t = 0;
	//mStateAll.col(1,mx);
	//for(int i=0;i<=floor(T/mdt);++i){
	//t += mdt;
	////if(t==50)
	////{
	////mA=SquareMatrix({-1,-0.0,3,-1});
	////}
	//mu = mCtrl->input(mx,mdt,t);
	//mx = mx + mdt*vec(mA * mx + mB * mu);
	////disp(mx);
	//mStateAll.col(i+1,mx);
	//}
}

//const std::vector<double>& Dynamical::x(const double t, simulate = &EU) 
//{
//if (t<mT) return mStateAll[floor(t/mdt)];

////.col(1,mx);
//for(unsigned int i=floor(mT/mdt); i<=floor(t/mdt);++i){
//mT += mdt;
////if(t==50)
////{
////mA=SquareMatrix({-1,-0.0,3,-1});
////}
////if (mCtrl!=nullptr)
//mu = *mInFun(mx,mB.size()[0],mT);
////mu = mCtrl->input(mx,mdt,mT);
////else (mu=0);
//mx = mx + mdt*vec(mA * mx + mB * mu);
////disp(mx);
////mStateAll.col(i+1,mx);
//mStateAll.insert({i+1,mx});
//mInputAll.insert({i+1,mu});
//}
//return mx;

//}
	template <typename E>
const std::vector<double>& Dynamical::x(const double t, const std::vector<double> x0) 
{
	if (mT==0 || x0!=mStateAll[0]) 
	{
		//(this->*simFun)(0, x0,t);
		Run<E>(0, x0,t);
		return mx;

	}
	else if (t<mT) 
		return mStateAll[floor(t/mdt)];
	else{

		//(this->*simFun)(mT, mx, t);
		Run<E>(mT, mx, t);
		return mx;
	}


	//return .col(floor(t/mdt));
}


void Dynamical::printAll() const
{
	std::cout<< "A is:" << std::endl;
	disp(mA);
	std::cout<< "B is:" << std::endl;
	disp(mB);
	std::cout<< "current state is:" << std::endl;
	disp(mx);
}


	template <typename E>
const std::vector<double> Dynamical::Run(const double t0, const std::vector<double>& x, const double t)
{
	mT = t0;
	mx = x;

	E sim(mA,mB,mdt);
	for(unsigned int i=floor(mT/mdt);i<=floor(t/mdt);++i){

		mx = sim.linear(mx, mCtrl, mInFun, mT);
		mT += mdt;
		//if (mCtrl!=nullptr)
		//mu = mCtrl->input(mx,mdt,mT);
		//else if (mInFun!=nullptr)
		//mu = mInFun(mx,mB.size()[0],mT);

		mStateAll[i] = mx;
		mInputAll[i] = mu;

		//mx = mx + mdt*vec(mA * mx + mB * mu);
	}
	return mx;
}	



//const std::vector<double> Dynamical::RK(const double t0, const std::vector<double>& x, const double t)
//{
//mT = t0;
//mx = x;

//for(unsigned int i=floor(mT/mdt);i<=floor(t/mdt);++i){

////if(t==50)
////{
////mA=SquareMatrix({-1,-0.0,3,-1});
////}
//if (mCtrl!=nullptr)
//mu = mCtrl->input(mx,mdt,mT);
//else if (mInFun!=nullptr)
//mu = mInFun(mx,mB.size()[0],mT);

//mStateAll[i] = mx;
//mInputAll[i] = mu;

//mT += mdt;
//mx = mx + mdt*vec(mA * mx + mB * mu);
////disp(mx);
////mStateAll.col(i+1,mx);
//}
//return mx;

//}
//}
