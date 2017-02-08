#ifndef STEP_H 
#define STEP_H 

#include <math.h>

namespace ADP{

	class Step
	{
		public:
			Step(const double a, const double b, const double c):ma(a),mb(b),mc(c){};
			const double stepOut(const unsigned int n){return ma/(mb+pow(n,mc));};
			void change(const double a, const double b, const double c){ma=a;mb=b;mc=c;};
		private:
			double ma;
			double mb;
			double mc;
	};


}

#endif
