#ifndef CONTROLLERS_H
#define CONTROLLERS_H
#include <vector>
#include "Others.h" 

namespace ADP{
	class Controllers{
		public: 
			typedef const std::vector<double> (*noise)(const unsigned int, const double);
			virtual const std::vector<double> input(const std::vector<double>& x, const double dt, const double t=0, noise noif=&sinusoidal)=0;
			virtual void dispAll() const =0;
			virtual ~Controllers(){};
	};
}

#endif
