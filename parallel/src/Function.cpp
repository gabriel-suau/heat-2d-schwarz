/*!
 * @file Function.cpp
 *
 * @brief Defines classes for the source terms and the boundary conditions.
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


#include "Function.h"
#include "MPIUtils.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Vector.h"

#include <cmath>
#include <fstream>


//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//

Function::Function(DataFile* DF):
  _DF(DF), _xmin(DF->getxMin()), _ymin(DF->getyMin()), _xmax(DF->getxMax()), _ymax(DF->getyMax()), _Lx(DF->getLx()), _Ly(DF->getLy()), _dx(_DF->getDx()), _dy(_DF->getDy()), _Nx(DF->getNx()), _Ny(DF->getNy()), _N(localSize), _Sol0(localSize), _sourceTerm(localSize), _exactSol(localSize)
{
}


/*!
 * @details Resizes the initial condition, the source term and the exact solution. Sets
 * the initial condition to 1 in the whole domain.
 */
void Function::Initialize()
{
  int i;
  // Logs de début
#if VERBOSITY>0
  if (MPI_Rank == 0) {
    std::cout << "====================================================================================================" << std::endl;
    std::cout << "Building initial condition..." << std::endl;
  }
#endif

  // Resize la CI, le terme source et la solution exacte
  _Sol0.resize(_N);
  _sourceTerm.resize(_N);
  _exactSol.resize(_N);
  for (i = 0 ; i < _N ; ++i) {
    _Sol0[i] = 1.0;
  }

  // Logs de fin
#if VERBOSITY>0
  if (MPI_Rank == 0) {
    std::cout << termcolor::green << "SUCCESS::FUNCTION : Initial Condition was successfully built." << std::endl;
    std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
  }
#endif
}


/*!
 * @details Builds the source term vector while taking into account the boundary conditions.
 *
 * @param t Current time (used when the source term is time dependent, i.e.in scenario 3).
 */
void Function::buildSourceTerm(double t)
{
  int i, j, k, kloc;
  double x, y, D, one_over_dx2, one_over_dy2;

  D = _DF->getDiffCoeff();
  one_over_dx2 = 1.0 / (_dx * _dx);
  one_over_dy2 = 1.0 / (_dy * _dy);

  for (k = kBegin ; k <= kEnd ; k++) {
    // Intérieur du domaine
    i = k % _Nx;
    j = k / _Nx;
    x = _xmin + (i+1) * _dx;
    y = _ymin + (j+1) * _dy;
    kloc = k - kBegin;
    _sourceTerm[kloc] = f(x, y, t);
    // Conditions aux limites
    // Bord bas
    if (j == 0) {
      _sourceTerm[kloc] += D * g(x, _ymin, t) * one_over_dy2;
    }
    // Bord haut
    else if (j == _Ny - 1) {
      _sourceTerm[kloc] += D * g(x, _ymax, t) * one_over_dy2;
    }
    // Bord gauche
    if (i == 0) {
      _sourceTerm[kloc] += D * h(_xmin, y, t) * one_over_dx2;
    }
    // Bord droit
    else if (i == _Nx - 1) {
      _sourceTerm[kloc] += D * h(_xmax, y, t) * one_over_dx2;
    }
  }

}

/*!
 * @param fileName Name of the file in which to write the exact solution.
 */
void Function::saveCurrentExactSolution(std::string &fileName) const
{
  int i, j, k;
  double x, y;
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(10);

  for (k = kBegin ; k <= kEnd ; ++k) {
    j = k / _Nx;
    i = k % _Nx;
    x = _xmin + (i+1) * _dx;
    y = _ymin + (j+1) * _dy;
    outputFile << x << " " << y << " " << _exactSol[k - kBegin] << std::endl;
  }
}

/*!
 * @param fileName Name of the file in which to write the source term.
 *
 * @deprecated This function is only here for debugging purposes. It's never used in the final code.
 */
void Function::saveSourceTerm(std::string& fileName) const
{
  int i, j, k;
  double x, y;
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(10);

  for (k = kBegin ; k <= kEnd ; ++k) {
    j = k / _Nx;
    i = k % _Nx;
    x = _xmin + (i+1) * _dx;
    y = _ymin + (j+1) * _dy;
    outputFile << x << " " << y << " " << _sourceTerm[k - kBegin] << std::endl;
  }
}


//--------------------------------------------------//
//--------------------Function 1--------------------//
//--------------------------------------------------//

Function1::Function1(DataFile* DF):
  Function(DF)
{
}


double Function1::f(const double x, const double y, const double t)
{
  return 2 * (y - y * y + x - x * x);
}


double Function1::g(const double x, const double y, const double t)
{
  return 0.;
}


double Function1::h(const double x, const double y, const double t)
{
  return 0.;
}


void Function1::buildExactSolution(double t)
{
  int i, j, k;
  double x, y;

  for (k = kBegin ; k <= kEnd ; ++k) {
    j = k / _Nx;
    i = k % _Nx;
    x = _xmin + (i+1) * _dx;
    y = _ymin + (j+1) * _dy;
    _exactSol[k - kBegin] = x * (1 - x) * y * (1 - y);
  }
}


//--------------------------------------------------//
//--------------------Function 2--------------------//
//--------------------------------------------------//

Function2::Function2(DataFile* DF):
  Function(DF)
{
}


double Function2::f(const double x, const double y, const double t)
{
  return sin(x) + cos(y);
}


double Function2::g(const double x, const double y, const double t)
{
  return sin(x) + cos(y);
}


double Function2::h(const double x, const double y, const double t)
{
  return sin(x) + cos(y);
}


void Function2::buildExactSolution(double t)
{
  int i, j, k;
  double x, y;
  
  for (k = kBegin ; k <= kEnd ; ++k) {
    j = k / _Nx;
    i = k % _Nx;
    x = _xmin + (i+1) * _dx;
    y = _ymin + (j+1) * _dy;
    _exactSol[k - kBegin] = sin(x) + cos(y);
  }
}


//--------------------------------------------------//
//--------------------Function 3--------------------//
//--------------------------------------------------//

Function3::Function3(DataFile* DF):
  Function(DF)
{
}


double Function3::f(const double x, const double y, const double t)
{
  return exp(-pow(x - 0.5 * _Lx, 2)) * exp(-pow(y - 0.5 * _Ly, 2)) * cos(0.5 * M_PI * t);
}


double Function3::g(const double x, const double y, const double t)
{
  return 0.;
}


double Function3::h(const double x, const double y, const double t)
{
  return 1.;
}


void Function3::buildExactSolution(double t)
{
  _exactSol.resize(localSize, 0.);
}
