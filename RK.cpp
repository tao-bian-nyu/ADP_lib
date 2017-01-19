#include "RK.h"
#include "MatrixCalc.h"

namespace ADP{

	const std::vector<double> RK::linear(const std::vector<double>& x, Controllers* controller, inputfun input, const double t)
	{
			if (controller!=nullptr)
			{
				std::vector<double> u = controller->input(x,mdt,t);
				std::vector<double> k1 (vec(mA * x + mB * u));
				u = controller->input(x+k1*mdt/2,0,t+mdt/2);
				std::vector<double> k2 (vec(mA * (x+k1*mdt/2) + mB * u));
				u = controller->input(x+k2*mdt/2,0,t+mdt/2);
				std::vector<double> k3 (vec(mA * (x+k2*mdt/2) + mB * u));
				u = controller->input(x+k3*mdt,0,t+mdt);
				std::vector<double> k4 (vec(mA * (x+k3*mdt) + mB * u));
				return std::vector<double>(x+mdt/6*(k1+2*k2+2*k3+k4));

			}
			else if (input!=nullptr)
			{
				std::vector<double> u = input(x,mB.size()[0],t);
				std::vector<double> k1 (vec(mA * x + mB * u));
				u = input(x+k1*mdt/2,mB.size()[0],t+mdt/2);
				std::vector<double> k2 (vec(mA * (x+k1*mdt/2) + mB * u));
				u = input(x+k2*mdt/2,mB.size()[0],t+mdt/2);
				std::vector<double> k3 (vec(mA * (x+k2*mdt/2) + mB * u));
				u = input(x+k3*mdt,mB.size()[0],t+mdt);
				std::vector<double> k4 (vec(mA * (x+k3*mdt) + mB * u));


				return std::vector<double>(x+mdt/6*(k1+2*k2+2*k3+k4));
			}
			else
			{

				std::vector<double> k1 (vec(mA * x));
				std::vector<double> k2 (vec(mA * (x+k1*mdt/2)));
				std::vector<double> k3 (vec(mA * (x+k2*mdt/2)));
				std::vector<double> k4 (vec(mA * (x+k3*mdt)));

				return std::vector<double>(x+mdt/6*(k1+2*k2+2*k3+k4));

			}


	}
}
