#ifndef _SLOWNESSVECTORMATRIX_H_
#define _SLOWNESSVECTORMATRIX_H_
#include <sstream>
#include <vector>
#include "slowness.h"
#include "SeisppError.h"
namespace PWMIG
{
using namespace PWMIG;
using namespace SEISPP;
class SlownessVectorMatrix
{
public:
    /*! Default constructor - 0x0 */
    SlownessVectorMatrix();
    /*! \brief Initializing constructor fills with 0 vectors. 

      This is the main constructor for this object.   Because 
      of the somewhat odd use it is totally nonstandard in its
      behavior.   It does not build a complete object, but initializes
      a matrix of slowness vectors of size nrowxncol to all zeros.
      The operator()method should be used to fill the matrix with 
      actual data.   This was done because the entire purpose of this
      object is to divorce pwmigcore routines from global travel
      time calculators.   This was a maintenance issue choice 
      created by the lack of a clean solution for C++ to the dsap
      version of the taup library and not common alternative.  

      \param nrow - number of rows grid (matrix) of SlownessVectors
      \param ncol - number of columns grid (matrix) of SlownessVectors
      */
    SlownessVectorMatrix(int nrow, int ncol);
    /*! Copy constructor */
    SlownessVectorMatrix(const SlownessVectorMatrix& parent);
    /*! \brief return number of rows in this matrix */
    int rows()
    {
        return nrow;
    };
    /*! return number of columns in this matrix */
    int columns()
    {
        return ncol;
    };
    /*! \brief overloaded subscript operator.

      \param i is row index (C order)
      \param j is col index (C order)

      \exception - Will throw a SeisppError if the index is outside
         row and column range.
      */
    SlownessVector& operator()(int i, int j);
    SlownessVectorMatrix& operator=(const SlownessVectorMatrix& parent);
    friend ostream& operator<<(ostream& ostr, SlownessVectorMatrix& svm);
private:
    int nrow;
    int ncol;
    vector<SlownessVector> uarray;
};
} // End namespace encapsulation 
#endif
