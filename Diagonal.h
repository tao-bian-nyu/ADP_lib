#ifndef DIAGONAL_H
#define DIAGONAL_H
//#include <vector>
#include "SymmetricMatrix.h"

namespace ADP
{

  class Diagonal:public SymmetricMatrix
  {
  public:
    explicit Diagonal(const int ncol = 1, const double val=1.0);
    Diagonal(const Diagonal& mat):SymmetricMatrix(mat){};
    Diagonal(const Matrix& mat);
    explicit Diagonal(const std::vector<double> & input);
    virtual ~Diagonal(){};
  };

}
#endif
