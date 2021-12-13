/*!
 * @file TimeScheme.cpp
 *
 * @brief Defines classes for time integration.
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


#include "TimeScheme.h"
#include "Vector.h"
#include "MPIUtils.h"
#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//

/*!
 * @param DF DataFile object.
 * @param function Function object.
 * @param laplacian Laplacian object.
 */
TimeScheme::TimeScheme(DataFile* DF, Function* function, Laplacian* laplacian):
  _DF(DF), _function(function), _laplacian(laplacian), _dt(DF->getTimeStep()), _t0(DF->getInitialTime()), _tf(DF->getFinalTime()), _t(_t0), _resultsDir(DF->getResultsDirectory())
{
  _pSol = new DVector(_function->getInitialCondition());
  _ptmp = new DVector(localSize);
}


/*!
 * @param fileName File in which the solution must be saved.
 */
void TimeScheme::saveCurrentSolution(std::string &fileName) const
{
  int nx, i, j, k;
  double xmin, ymin, dx, dy, x, y;

  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(10);

  // Récupération des variables utiles
  nx = _DF->getNx();
  xmin = _DF->getxMin();
  ymin = _DF->getyMin();
  dx = _DF->getDx();
  dy = _DF->getDy();

  for (k = kBegin ; k <= kEnd ; ++k) {
    j = k / nx;
    i = k % nx;
    x = xmin + (i+1) * dx;
    y = ymin + (j+1) * dy;
    outputFile << x << " " << y << " " << (*_pSol)[k - kBegin] << std::endl;
  }
}


/*!
 * @details Main function of the TimeScheme class. It performs the time iterations of the selected time integration methods until the final time is reached. It also saves the solution, and the residuals when they need to be saved.
*/
void TimeScheme::solve()
{
  int n, scenario;
  std::string solFileName;

  // Logs de début
#if VERBOSITY>0
  if (MPI_Rank == 0) {
    std::cout << "====================================================================================================" << std::endl;
    std::cout << "Time loop..." << std::endl;
  }
#endif

  // Variables pratiques
  n = 0;
  scenario = _DF->getScenario();

  // Sauvegarde la condition initiale
  solFileName = _resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n) + ".dat";
  std::cout << solFileName << std::endl;
  saveCurrentSolution(solFileName);

  // Démarrage du chrono
  auto start = MPI_Wtime();

  // Construction de la matrice
  _laplacian->buildMat();

  // Boucle en temps
  while (_t < _tf) {
    oneStep();
    ++n;
    _t += _dt;

    // Save the numerical solution
    if (!_DF->isSaveFinalResultOnly() && n % _DF->getSaveFrequency() == 0) {
#if VERBOSITY>0
      if (MPI_Rank == 0)
        std::cout << "Saving solution at t = " << _t << std::endl;
#endif
      solFileName = _resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".dat";
      saveCurrentSolution(solFileName);
    }
  }

  // Calcul et sauvegarde du temps CPU
  auto finish = MPI_Wtime();
  auto cpuTime = finish - start;
  if (MPI_Rank == 0) {
    std::string cputimeFileName(_DF->getErrorAndCPUTimeDir() + "/cputime.dat");
    std::ofstream cputimeFile(cputimeFileName, std::ios::app);
    cputimeFile.precision(10);
    cputimeFile << MPI_Size << " " << _DF->getNx() << " " << _DF->getNy() << " " << _DF->getNx() * _DF->getNy() << " " << _DF->getDx() << " " << _DF->getDy() << " " <<_DF->getDx() * _DF->getDy() << " " << cpuTime << std::endl;
  }


  // Save the final solution
  if (_DF->isSaveFinalResultOnly()) {
#if VERBOSITY>0
    if (MPI_Rank == 0)
      std::cout << "Saving solution at t = " << _t << std::endl;
#endif
    std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".dat");
    saveCurrentSolution(solFileName);
  }

  // Save the exact solution
  if (_DF->getScenario() == 1 || _DF->getScenario() == 2) {
    _function->buildExactSolution(_t);
#if VERBOSITY>0
    if (MPI_Rank == 0)
      std::cout << "Saving exact solution at t = " << _t << std::endl;
#endif
    std::string exactSolFileName(_resultsDir + "/solution_exacte_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + ".dat");
    _function->saveCurrentExactSolution(exactSolFileName);
  }

  // Calcul et sauvegarde de l'erreur pour les scenario 1 et 2
  if (_DF->getScenario() == 1 || _DF->getScenario() == 2) {
    double L2error(computeCurrentL2Error());
    double L1error(computeCurrentL1Error());
    if (MPI_Rank == 0) {
      std::string errorFileName(_DF->getErrorAndCPUTimeDir() + "/error.dat");
      std::ofstream errorFile(errorFileName, std::ios::app);
      errorFile.precision(10);
      errorFile << _DF->getNx() << " " << _DF->getNy() << " " << _DF->getNx() * _DF->getNy() << " " << _DF->getDx() << " " << _DF->getDy() << " " <<_DF->getDx() * _DF->getDy() << " " << L2error << " " << L1error << std::endl;
      std::cout << "Error L2 = " << L2error << " at t = " << _t << " for Nx = " << _DF->getNx() << ", Ny = " << _DF->getNy() << std::endl << std::endl;
      std::cout << "Error L1 = " << L1error << " at t = " << _t << " for Nx = " << _DF->getNx() << ", Ny = " << _DF->getNy() << std::endl << std::endl;
    }
  }

  // Logs de fin
#if VERBOSITY>0
  if (MPI_Rank == 0) {
    std::cout << termcolor::green << "SUCCESS::TIMESCHEME : Time loop completed successfully in " << cpuTime << " seconds !" << std::endl;
    std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
  }
#endif
}

double TimeScheme::computeCurrentL2Error()
{
  DVector errorVec;
  double error;

  errorVec = (*_pSol) - _function->getExactSolution();
  error = errorVec.dot(errorVec);
  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = sqrt(_DF->getDx() * _DF->getDy() * error);

  return error;
}

double TimeScheme::computeCurrentL1Error()
{
  DVector errorVec;
  double error;
  int i;

  errorVec = (*_pSol) - _function->getExactSolution();
  error = 0.0;

  for (i = 0 ; i < localSize ; ++i) {
    error += std::abs(errorVec[i]);
  }

  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = _DF->getDx() * _DF->getDy() * error;

  return error;
}


//-------------------------------------------------------------//
//--------------------Explicit Euler scheme--------------------//
//-------------------------------------------------------------//

ExplicitEuler::ExplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian):
  TimeScheme(DF, function, laplacian)
{
}


void ExplicitEuler::oneStep()
{
  int i, maxIt;
  double tol, error;

  i = 0;
  maxIt = _DF->getSchwarzMaxIterations();
  tol = _DF->getSchwarzTolerance();
  error = computeCurrentL2Error();

  // Calcul du terme source
  _function->buildSourceTerm(_t);

  // Schwarz
  while (i < maxIt && error > tol) {
    _function->updateSourceTerm(*_pSol);
    _laplacian->updateMat();
    (*_pSol) = (*_pSol) + _laplacian->matVecProd(*_pSol) + _dt * _function->getSourceTerm();
    error = computeCurrentL2Error();
    i++;
  }

}


//-------------------------------------------------------------//
//--------------------Implicit Euler scheme--------------------//
//-------------------------------------------------------------//

ImplicitEuler::ImplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian):
  TimeScheme(DF, function, laplacian)
{
}


void ImplicitEuler::oneStep()
{
  int i, maxIt;
  double tol, error;

  i = 0;
  maxIt = _DF->getSchwarzMaxIterations();
  tol = _DF->getSchwarzTolerance();
  error = computeCurrentL2Error();

  // Calcul du terme source
  _function->buildSourceTerm(_t + _dt);

  // Schwarz
  while (i < maxIt && error > tol) {
    _function->updateSourceTerm(*_pSol);
    _laplacian->updateMat();
    if (_DF->getLinearSolver() == "CG") {
      _laplacian->CG(*_pSol + _dt * _function->getSourceTerm(), *_pSol, _DF->getTolerance(), _DF->getMaxIterations());
    }
    else if (_DF->getLinearSolver() == "BICGSTAB") {
      _laplacian->BICGSTAB(*_pSol + _dt * _function->getSourceTerm(), *_pSol, _DF->getTolerance(), _DF->getMaxIterations());
    }
    error = computeCurrentL2Error();
    i++;
  }

}
