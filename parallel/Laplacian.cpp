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
 * @param resFile Name of the file in which to save the residuals L2 norm.
 *
 * @return The solution of the linear system (a DVector object).
 */
DVector Laplacian::solveConjGrad(const DVector& b, const DVector& x0, double tolerance, int maxIterations, std::ofstream& resFile)
{
  // Variables intermédiaires
  DVector x(x0);
  DVector res(b - this->matVecProd(x0));
  DVector p(res);
  // Compute the initial global residual
  double resDotRes(res.dot(res));
  MPI_Allreduce(MPI_IN_PLACE, &resDotRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double beta(sqrt(resDotRes));
  if (MPI_Rank == 0)
    resFile << beta << std::endl;

  // Iterations of the method
  int k(0);
  while ((beta > tolerance) && (k < maxIterations))
    {
      // Matvec product
      DVector z(this->matVecProd(p));
      // Dot products
      double resDotP(res.dot(p)), zDotP(z.dot(p));
      MPI_Allreduce(MPI_IN_PLACE, &resDotP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &zDotP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      // Compute alpha
      double alpha(resDotP/zDotP);
      // Update the solution
      x = x + alpha * p;
      // Update the residual
      res = res - alpha * z;
      // Dot product
      double resDotRes(res.dot(res));
      MPI_Allreduce(MPI_IN_PLACE, &resDotRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      // Compute gamma
      double gamma(resDotRes/pow(beta,2));
      // Update p
      p = res + gamma * p;
      beta = sqrt(resDotRes);
      ++k;
      if (MPI_Rank == 0)
        resFile << beta << std::endl;
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

  return x;
}
