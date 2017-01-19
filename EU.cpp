#include "EU.h"
#include "MatrixCalc.h"

namespace ADP{

	const std::vector<double> EU::linear(const std::vector<double>& x, Controllers* controller, inputfun input, const double t)
	{
			if (controller!=nullptr)
			{
				std::vector<double> u = controller->input(x,mdt,t);
				return (x+mdt*vec(mA * x + mB * u));

			}
			else if (input!=nullptr)
				
			{
				std::vector<double> u = input(x,mB.size()[0],t);
				return (x+mdt*vec(mA * x + mB * u));
			}
			else
				return (x+mdt*vec(mA * x));
			

	}
}
