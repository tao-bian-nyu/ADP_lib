#include "Matrix.h"
#include "Dynamical.h"
#include <iostream>
#include <math.h>
#include "MatrixCalc.h"
#include "Others.h"

namespace ADP
{

	Dynamical::Dynamical(const SquareMatrix A, const Matrix B, const std::vector<double> x0, inputfun input, const long double dt)
		:mA(A),mB(B),mx(x0),mu(B.size()[0],0),mT(0),mdt(dt),mCtrl(nullptr)
		 //,mxT(floor(T/dt),A.size()[0],0)
	{
		mxT.insert({0,mx});
		muT.insert({0,mu});
		//double t = 0;
		//mxT.col(1,mx);
		//for(int i=0;i<=floor(T/mdt);++i){
			//t += mdt;
			//mu = input(mx,B.size()[0],t);
			//mx = mx + mdt*vec(mA * mx + mB * mu);
			//mxT.col(i+1,mx);
		//}
	}
	 
	Dynamical::Dynamical(const SquareMatrix A, const std::vector<double> x0,  const long double dt)
		:mA(A),mB(1,A.size()[0],0),mx(x0),mu(1,0),mT(0),mdt(dt),mCtrl(nullptr)
		 //,mxT(floor(T/dt),A.size()[0],0)
	{
		mxT.insert({0,mx});
		muT.insert({0,mu});
		//mxT.col(1,mx);
		//for(int i=0;i<=floor(T/mdt);++i){
			//mx = mx + mdt*vec(mA * mx);
			//mxT.col(i+1,mx);
		//}
	}

	Dynamical::Dynamical(const SquareMatrix A, const Matrix B, const std::vector<double> x0, Controllers* controller, const long double dt)
		:mA(A),mB(B),mx(x0),mu(B.size()[0],0),mT(0),mdt(dt),mCtrl(controller)
		 //,mxT(floor(T/dt),A.size()[0],0)
	{
		mxT.insert({0,mx});
		muT.insert({0,mu});
		//double t = 0;
		//mxT.col(1,mx);
		//for(int i=0;i<=floor(T/mdt);++i){
			//t += mdt;
			////if(t==50)
			////{
				////mA=SquareMatrix({-1,-0.0,3,-1});
			////}
			//mu = mCtrl->input(mx,mdt,t);
			//mx = mx + mdt*vec(mA * mx + mB * mu);
			////disp(mx);
			//mxT.col(i+1,mx);
		//}
	}

	const std::vector<double>& Dynamical::x(const double t, inputfun input) 
	{
		if (t<mT) return mxT[floor(t/mdt)];

		//mxT.col(1,mx);
		for(unsigned int i=floor(mT/mdt); i<=floor(t/mdt);++i){
			mT += mdt;
			//if(t==50)
			//{
				//mA=SquareMatrix({-1,-0.0,3,-1});
			//}
			//if (mCtrl!=nullptr)
			mu = input(mx,mB.size()[0],mT);
				//mu = mCtrl->input(mx,mdt,mT);
			//else (mu=0);
			mx = mx + mdt*vec(mA * mx + mB * mu);
			//disp(mx);
			//mxT.col(i+1,mx);
			mxT.insert({i+1,mx});
			muT.insert({i+1,mu});
		}
		return mx;
	
	}
	const std::vector<double>& Dynamical::x(const double t) 
	{
		if (t<mT) return mxT[floor(t/mdt)];

		//mxT.col(1,mx);
		for(unsigned int i=floor(mT/mdt);i<=floor(t/mdt);++i){
			mT += mdt;
			//if(t==50)
			//{
				//mA=SquareMatrix({-1,-0.0,3,-1});
			//}
			if (mCtrl!=nullptr)
				mu = mCtrl->input(mx,mdt,mT);
			//else (mu=0);
			mx = mx + mdt*vec(mA * mx + mB * mu);
			//disp(mx);
			//mxT.col(i+1,mx);
			mxT.insert({i+1,mx});
			muT.insert({i+1,mu});
		}
		return mx;
		//return mxT.col(floor(t/mdt));
	}

	
	void Dynamical::printAll() const
	{
		std::cout<< mT << std::endl;
		disp(mx);
	}

}
