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
#include <cassert>

Laplacian::Laplacian(DataFile* DF, Function* function):
  _DF(DF), _function(function), _N(localSize)
{
}


void Laplacian::setFromTriplets(const std::vector<triplet_t>& triplets)
{
  triplet_t triplet;
  int i;

  assert (triplets.size() == _NNZ);

  _AA.resize(_NNZ);
  _JA.resize(_NNZ);
  _IA.resize(_N + 1, 0);

  for (i = 0 ; i < _NNZ ; i++) {
    triplet = triplets[i];
    _AA[i] = triplet.coef;
    _JA[i] = triplet.j;
    _IA[triplet.i + 1]++;
  }

  for (i = 0 ; i < _N ; i++) {
    _IA[i+1] += _IA[i];
  }

}


void Laplacian::buildMat() {
  int i, Nx, Ny, count;
  double dx, dy, dt, one_over_dx2, one_over_dy2;
  double alpha, beta, gamma, D;
  std::vector<triplet_t> triplets;

  Nx = _DF->getNx();
  Ny = _DF->getNy();
  dx = _DF->getDx();
  dy = _DF->getDy();
  dt = _DF->getTimeStep();
  D = _DF->getDiffCoeff();

  one_over_dx2 = 1.0 / (dx * dx);
  one_over_dy2 = 1.0 / (dy * dy);

  if (_DF->getTimeScheme() == "ExplicitEuler") {
    alpha = - 2.0 * D * dt * (one_over_dx2 + one_over_dy2);
    beta = dt * D * one_over_dx2;
    gamma = dt * D * one_over_dy2;
  }
  else if (_DF->getTimeScheme() == "ImplicitEuler") {
    alpha = 1.0 + 2.0 * D * dt * (one_over_dx2 + one_over_dy2);
    beta = - dt * D * one_over_dx2;
    gamma = - dt * D * one_over_dy2;
  }

  // _NNZ = _N + (_N - Nx) + (_N - _Nx) + (_N - Ny) + (_N - Ny);
  _NNZ = 5 * _N - 2 * Nx - 2 * Ny;
  count = 0;
  triplets.resize(_NNZ);

  for (i = 0 ; i < _N ; i++) {
    if (i >= Nx) triplets[count++] = {i, i-Nx, gamma};
    if (i % Nx != 0)  triplets[count++] = {i, i-1, beta};
    triplets[count++] = {i, i, alpha};
    if (i % Nx != Nx - 1) triplets[count++] = {i, i+1, beta};
    if (i < _N - Nx) triplets[count++] = {i, i+Nx, gamma};
  }

  assert (count == _NNZ);

  this->setFromTriplets(triplets);

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
  int i, k, k1, k2;
  result.resize(_N, 0.);

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
 * @param x Initial guess for the iterative solver (and solution at the end of the routine).
 * @param tolerance Tolerance for the residuals of the method.
 * @param maxIterations Maximum number of iterations. When we reach this number of iterations, the linear solver stops and returns the current result, even if the result has not converged enough.
 *
 * @return void
 */
void Laplacian::CG(const DVector& b, DVector& x, double tolerance, int maxIterations)
{
  // Variables
  DVector r, p, z;
  double rDotR, rDotP, zDotP, alpha, beta, gamma;
  int k;
  // Compute the initial ridual
  r = b - this->matVecProd(x);
  rDotR = r.dot(r);
  beta = sqrt(rDotR);

  p = r;

  // Iterations of the method
  k = 0;
  while ((beta > tolerance) && (k < maxIterations)) {
    z = this->matVecProd(p);
    rDotP = r.dot(p);
    zDotP = z.dot(p);
    alpha = rDotP/zDotP;
    x = x + alpha * p;
    r = r - alpha * z;
    rDotR = r.dot(r);
    gamma = rDotR/pow(beta,2);
    p = r + gamma * p;
    beta = sqrt(rDotR);
    ++k;
  }

  // Logs
#if VERBOSITY>1
  if (MPI_Rank == 0) {
    if ((k == maxIterations) && (beta > tolerance)) {
      std::cout << termcolor::yellow << "SOLVER::GC::WARNING : The GC method did not converge. Residual L2 norm = " << beta << " (" << maxIterations << " iterations)" << std::endl;
      std::cout << termcolor::reset;
    } else {
      std::cout << termcolor::green << "SOLVER::GC::SUCCESS : The GC method converged in " << k << " iterations ! Residual L2 norm = " << beta << std::endl;
      std::cout << termcolor::reset;
    }
  }
#endif

}

/*!
 * @param b Right hand side vector.
 * @param x Initial guess for the iterative solver (and solution at the end of the routine).
 * @param tolerance Tolerance for the residuals of the method.
 * @param maxIterations Maximum number of iterations. When we reach this number of iterations, the linear solver stops and returns the current result, even if the result has not converged enough.
 *
 * @return void
 */
void Laplacian::BICGSTAB(const DVector& b, DVector& x, double tolerance, int maxIterations)
{
  // Variables
  DVector r, rs, p, s, Ap, As;
  double rDotR, alpha, beta, gamma, omega;
  int k;
  // Compute the initial residual
  r = b - this->matVecProd(x);
  rDotR = r.dot(r);
  beta = sqrt(rDotR);
  p = r;
  rs = r;

  // Iterations of the method
  k = 0;
  while ((beta > tolerance) && (k < maxIterations)) {
    Ap = this->matVecProd(p);
    gamma = r.dot(rs);
    alpha = gamma / Ap.dot(rs);
    s = r - alpha * Ap;
    As = this->matVecProd(s);
    omega = As.dot(s) / As.dot(As);
    x = x + alpha * p + omega * s;
    r = s - omega * As;
    beta = r.dot(rs) * alpha / (omega * gamma);
    p = r + beta * (p - omega * Ap);
    ++k;
  }

  // Logs
#if VERBOSITY>1
  if (MPI_Rank == 0) {
    if ((k == maxIterations) && (beta > tolerance)) {
      std::cout << termcolor::yellow << "SOLVER::BICGSTAB::WARNING : The BICGSTAB method did not converge. Residual L2 norm = " << beta << " (" << maxIterations << " iterations)" << std::endl;
      std::cout << termcolor::reset;
    } else {
      std::cout << termcolor::green << "SOLVER::BICGSTAB::SUCCESS : The BICGSTAB method converged in " << k << " iterations ! Residual L2 norm = " << beta << std::endl;
      std::cout << termcolor::reset;
    }
  }
#endif

}
