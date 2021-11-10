/*!
 * @file Laplacian.cpp
 *
 * @brief Defines a class representing the discrete laplacian matrix.
 *
 * @authors Gabriel Suau, Lucas Trautmann, Geoffrey Lebaud
 *
 * @version 0.1.0
 *
 * @copyright © 2021 Gabriel Suau
 * @copyright © 2021 Lucas Trautmann
 * @copyright © 2021 Geoffrey Lebaud
 *
 * @copyright This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include "Laplacian.h"
#include "MPIUtils.h"
#include "Function.h"
#include "DataFile.h"
#include "Vector.h"

#include <cmath>
#include <mpi.h>


Laplacian::Laplacian(DataFile* DF, Function* function):
  _DF(DF), _function(function), _N(localSize)
{
}


void Laplacian::buildMat() {
  
}


void Laplacian::updateMat() {
  
}


/*!
 * @details The discrete laplacian matrix is very sparse (block tridiagonal), so
 * only the strictly necessary operations are performed for this matrix-vector product.
 *
 * @param x DVector multiplied by the matrix.
 *
 * @return The result of the matvec product (a DVector object).
 */
DVector Laplacian::matVecProd(const DVector& x)
{
  // Vecteur resultat
  DVector result;
  int size(_N);
  int i, k, k1, k2;
  result.resize(size, 0.);

  for (i = 0 ; i < _N ; i++) {
    k1 = _IA[i];
    k2 = _IA[i+1];
    for (k = k1 ; k < k2 ; k++) {
      result[i] += _AA[k] * x[_JA[k]];
    }
  }

  // // MPI Communications
  // // Each proc has to communicate with proc - 1 and proc + 1
  // if (MPI_Rank + 1 < MPI_Size)
  //   {
  //     MPI_Sendrecv(&x[size - _Nx], _Nx, MPI_DOUBLE, MPI_Rank + 1, 0,
  //                  &next[0], _Nx, MPI_DOUBLE, MPI_Rank + 1, 1,
  //                  MPI_COMM_WORLD, &status);
  //   }
  // if (MPI_Rank - 1 >= 0)
  //   {
  //     MPI_Sendrecv(&x[0], _Nx, MPI_DOUBLE, MPI_Rank - 1, 1,
  //                  &prev[0], _Nx, MPI_DOUBLE, MPI_Rank - 1, 0,
  //                  MPI_COMM_WORLD, &status);
  //   }

  return result;
}


/*!
 * @param b Right hand side vector.
 * @param x0 Initial guess for the iterative solver.
 * @param tolerance Tolerance for the residuals of the conjugate gradient method.
 * @param maxIterations Maximum number of iterations. When we reach this number of iterations, the linear solver stops and returns the current result, even if the result has not converged enough.
 *
 * @return void
 */
void Laplacian::solveConjGrad(const DVector& b, DVector& x, double tolerance, int maxIterations)
{
  // Variables intermédiaires
  DVector res(b - this->matVecProd(x));
  DVector p(res);
  // Compute the initial global residual
  double resDotRes(res.dot(res));
  MPI_Allreduce(MPI_IN_PLACE, &resDotRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double beta(sqrt(resDotRes));

  // Iterations of the method
  int k(0);
  while ((beta > tolerance) && (k < maxIterations))
    {
      DVector z(this->matVecProd(p));
      double resDotP(res.dot(p)), zDotP(z.dot(p));
      double alpha(resDotP/zDotP);
      x = x + alpha * p;
      res = res - alpha * z;
      double resDotRes(res.dot(res));
      double gamma(resDotRes/pow(beta,2));
      p = res + gamma * p;
      beta = sqrt(resDotRes);
      ++k;
    }

  // Logs
#if VERBOSITY>1
  if (MPI_Rank == 0)
    {
      if ((k == maxIterations) && (beta > tolerance))
        {
          std::cout << termcolor::yellow << "SOLVER::GC::WARNING : The GC method did not converge. Residual L2 norm = " << beta << " (" << maxIterations << " iterations)" << std::endl;
          std::cout << termcolor::reset;
        }
      else
        {
          std::cout << termcolor::green << "SOLVER::GC::SUCCESS : The GC method converged in " << k << " iterations ! Residual L2 norm = " << beta << std::endl;
          std::cout << termcolor::reset;
        }
    }
#endif

}
