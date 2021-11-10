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
TimeScheme::TimeScheme()
{
}


/*!
 * @param DF DataFile object.
 * @param function Function object.
 * @param laplacian Laplacian object.
 */
TimeScheme::TimeScheme(DataFile* DF, Function* function, Laplacian* laplacian):
  _DF(DF), _function(function), _laplacian(laplacian), _Sol(_function->getInitialCondition()), _timeStep(DF->getTimeStep()), _initialTime(DF->getInitialTime()), _finalTime(DF->getFinalTime()), _currentTime(_initialTime), _resultsDir(DF->getResultsDirectory()), _resFileName(_resultsDir + "/" + DF->getResFile())
{
}


/*!
 * @param DF DataFile object.
 * @param function Function object.
 * @param laplacian Laplacian object.
 *
 * @deprecated This method could be useful if we decided to construct an empty TimeScheme object.
 * But we never use it in this code...
 */
void TimeScheme::Initialize(DataFile* DF, Function* function, Laplacian* laplacian)
{
  _DF = DF;
  _function = function;
  _laplacian = laplacian;
  _Sol = _function->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
  _resultsDir = DF->getResultsDirectory();
  _resFileName = DF->getResFile();
}


/*!
 * @param fileName File in which the solution must be saved.
 */
void TimeScheme::saveCurrentSolution(std::string &fileName) const
{
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(10);

  // Récupération des variables utiles
  int Nx(_DF->getNx());
  double xmin(_DF->getxMin()), ymin(_DF->getyMin());
  double dx(_DF->getDx()), dy(_DF->getDy());

  for (int k(kBegin) ; k <= kEnd ; ++k)
    {
      int j(k/Nx), i(k%Nx);
      double x(xmin + (i+1) * dx), y(ymin + (j+1) * dy);
      outputFile << x << " " << y << " " << _Sol[k - kBegin] << std::endl;
    }
}


/*!
 * @details Main function of the TimeScheme class. It performs the time iterations of the selected time integration methods until the final time is reached. It also saves the solution, and the residuals when they need to be saved.
*/
void TimeScheme::solve()
{
  // Logs de début
#if VERBOSITY>0
  if (MPI_Rank == 0)
    {
      std::cout << "====================================================================================================" << std::endl;
      std::cout << "Time loop..." << std::endl;
    }
#endif

  // Variables pratiques
  int n(0);
  int scenario(_DF->getScenario());

  // Sauvegarde la condition initiale
  std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n) + ".dat");
  saveCurrentSolution(solFileName);

  // Démarrage du chrono
  auto start = MPI_Wtime();

  // Construction de la matrice
  _laplacian->buildMat();

  // Boucle en temps
  while (_currentTime < _finalTime)
    {
      oneStep();
      ++n;
      _currentTime += _timeStep;
      if (!_DF->isSaveFinalResultOnly() && n % _DF->getSaveFrequency() == 0)
        {
          // Save the numerical solution
#if VERBOSITY>0
          if (MPI_Rank == 0)
            std::cout << "Saving solution at t = " << _currentTime << std::endl;
#endif
          std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".dat");
          saveCurrentSolution(solFileName);
        }
    }

  // Calcul et sauvegarde du temps CPU
  auto finish = MPI_Wtime();
  double cpuTime = finish - start;
  if (MPI_Rank == 0)
    {
      std::string cputimeFileName(_DF->getErrorAndCPUTimeDir() + "/cputime.dat");
      std::ofstream cputimeFile(cputimeFileName, std::ios::app);
      cputimeFile.precision(10);
      cputimeFile << MPI_Size << " " << _DF->getNx() << " " << _DF->getNy() << " " << _DF->getNx() * _DF->getNy() << " " << _DF->getDx() << " " << _DF->getDy() << " " <<_DF->getDx() * _DF->getDy() << " " << cpuTime << std::endl;
    }


  // Save the final solution
  if (_DF->isSaveFinalResultOnly())
    {
#if VERBOSITY>0
      if (MPI_Rank == 0)
        std::cout << "Saving solution at t = " << _currentTime << std::endl;
#endif
      std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".dat");
      saveCurrentSolution(solFileName);
    }

  // Save the exact solution
  if (_DF->getScenario() == 1 || _DF->getScenario() == 2)
    {
      _function->buildExactSolution(_currentTime);
#if VERBOSITY>0
      if (MPI_Rank == 0)
        std::cout << "Saving exact solution at t = " << _currentTime << std::endl;
#endif
      std::string exactSolFileName(_resultsDir + "/solution_exacte_scenario_" + std::to_string(scenario) + "_" + std::to_string(MPI_Rank) + ".dat");
      _function->saveCurrentExactSolution(exactSolFileName);
    }

  // Calcul et sauvegarde de l'erreur pour les scenario 1 et 2
  if (_DF->getScenario() == 1 || _DF->getScenario() == 2)
    {
      double L2error(computeCurrentL2Error());
      double L1error(computeCurrentL1Error());
      if (MPI_Rank == 0)
        {
          std::string errorFileName(_DF->getErrorAndCPUTimeDir() + "/error.dat");
          std::ofstream errorFile(errorFileName, std::ios::app);
          errorFile.precision(10);
          errorFile << _DF->getNx() << " " << _DF->getNy() << " " << _DF->getNx() * _DF->getNy() << " " << _DF->getDx() << " " << _DF->getDy() << " " <<_DF->getDx() * _DF->getDy() << " " << L2error << " " << L1error << std::endl;
        }
      if (MPI_Rank == 0)
        {
          std::cout << "Error L2 = " << L2error << " at t = " << _currentTime << " for Nx = " << _DF->getNx() << ", Ny = " << _DF->getNy() << std::endl << std::endl;
          std::cout << "Error L1 = " << L1error << " at t = " << _currentTime << " for Nx = " << _DF->getNx() << ", Ny = " << _DF->getNy() << std::endl << std::endl;
        }
    }

  // Logs de fin
#if VERBOSITY>0
  if (MPI_Rank == 0)
    {
      std::cout << termcolor::green << "SUCCESS::TIMESCHEME : Time loop completed successfully in " << cpuTime << " seconds !" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
    }
#endif
}

double TimeScheme::computeCurrentL2Error()
{
  double error(0.);
  DVector errorVec(_Sol - _function->getExactSolution());
  error = errorVec.dot(errorVec);
  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = sqrt(_DF->getDx() * _DF->getDy() * error);
  return error;
}

double TimeScheme::computeCurrentL1Error()
{
  double error(0.);
  DVector errorVec(_Sol - _function->getExactSolution());
  for (int i(0) ; i < localSize ; ++i)
    {
      error += std::abs(errorVec[i]);
    }
  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = _DF->getDx() * _DF->getDy() * error;
  return error;
}


//-------------------------------------------------------------//
//--------------------Explicit Euler scheme--------------------//
//-------------------------------------------------------------//
ExplicitEuler::ExplicitEuler():
  TimeScheme()
{
}


ExplicitEuler::ExplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian):
  TimeScheme(DF, function, laplacian)
{
}


void ExplicitEuler::Initialize(DataFile* DF, Function* function, Laplacian* laplacian)
{
  _DF = DF;
  _function = function;
  _laplacian = laplacian;
  _Sol = _function->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}


void ExplicitEuler::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);

  // Calcul du terme source
  _function->buildSourceTerm(_currentTime + dt);

  // Update la matrice du laplacien si besoin (conditions mixtes)
  _laplacian->updateMat();

  // Calcul de la solution
  _Sol = _Sol + _laplacian->matVecProd(_Sol) + dt * _function->getSourceTerm();
}


//-------------------------------------------------------------//
//--------------------Implicit Euler scheme--------------------//
//-------------------------------------------------------------//
ImplicitEuler::ImplicitEuler():
  TimeScheme()
{
}


ImplicitEuler::ImplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian):
  TimeScheme(DF, function, laplacian)
{
}


void ImplicitEuler::Initialize(DataFile* DF, Function* function, Laplacian* laplacian)
{
  _DF = DF;
  _function = function;
  _laplacian = laplacian;
  _Sol = _function->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}


void ImplicitEuler::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);
  double tolerance(_DF->getTolerance());
  int maxIt(_DF->getMaxIterations());
  // Ouverture du fichier pour les résidus
  std::ofstream resFile(_resFileName, std::ios::app);
  resFile.precision(10);

  // Calcul du terme source
  _function->buildSourceTerm(_currentTime + dt);

  // Update la matrice du laplacien si besoin (conditions mixtes)
  _laplacian->updateMat();

  // Calcul de la solution
  _Sol = _laplacian->solveConjGrad(_Sol + dt * _function->getSourceTerm(), _Sol, tolerance, maxIt, resFile);
}
