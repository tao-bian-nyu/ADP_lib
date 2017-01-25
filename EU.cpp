#include "EU.h"
#include "MatrixCalc.h"

namespace ADP{

	//const std::vector<double> EU::linear(const std::vector<double>& x, Controllers* controller, inputfun input, const double t)
	void EU::linear(const std::vector<double>& x, Controllers* controller, inputfun input, const double t)
	{
			if (controller!=nullptr)
			{
				std::vector<double> u = controller->input(x,mdt,t);
				mu = u;
				mx = (x+mdt*vec(mA * x + mB * u));

			}
			else if (input!=nullptr)
			{
				std::vector<double> u = input(x,mB.size()[0],t);
				mu = u;
				mx = (x+mdt*vec(mA * x + mB * u));
			}
			else
				mx = (x+mdt*vec(mA * x));

	}


	void EU::dispAll() const
	{
		std::cout << "mA is";
		disp(mA);
		std::cout << "mB is";
		disp(mB);
		std::cout << "mx is";
		disp(mx);
		std::cout << "mu is";
		disp(mu);
		std::cout << "mdt is " << mdt << std::endl;
	}
}
