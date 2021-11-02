#include "TimeScheme.h"
#include "Vector.h"
#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>

//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//
TimeScheme::TimeScheme()
{
}

TimeScheme::TimeScheme(DataFile* DF, Function* function, Laplacian* laplacian):
  _DF(DF), _function(function), _laplacian(laplacian), _Sol(_function->getInitialCondition()), _timeStep(DF->getTimeStep()), _initialTime(DF->getInitialTime()), _finalTime(DF->getFinalTime()), _currentTime(_initialTime), _resultsDir(DF->getResultsDirectory()), _resFileName(_resultsDir + "/" + DF->getResFile())
{
}

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

void TimeScheme::saveCurrentSolution(std::string &fileName) const
{
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(10);

  // Récupération des variables utiles
  int Nx(_DF->getNx()), Ny(_DF->getNy());
  double xmin(_DF->getxMin()), ymin(_DF->getyMin());
  double dx(_DF->getDx()), dy(_DF->getDy());

  outputFile << "# vtk DataFile Version 3.0" << std::endl;
  outputFile << "sol" << std::endl;
  outputFile << "ASCII" << std::endl;
  outputFile << "DATASET STRUCTURED_POINTS" << std::endl;
  outputFile << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl;
  outputFile << "ORIGIN " << xmin << " " << ymin << " " << 0 << std::endl;
  outputFile << "SPACING " << dx << " " << dy << " " << 1 << std::endl;;
  outputFile << "POINT_DATA " << Nx*Ny << std::endl;
  outputFile << "SCALARS sol float" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;

  for(int j=0; j<Ny; ++j)
    {
      for(int i=0; i<Nx; ++i)
        {
          outputFile << _Sol[i+j*Nx] << " ";
        }
      outputFile << std::endl;
    }

  // for(int j=0; j<Ny; ++j)
  //   {
  //     for(int i=0; i<Nx; ++i)
  //       {
  //         double x(xmin + (i+1) * dx), y(ymin + (j+1) * dy);
  //         outputFile << x << " " << y << " " << _Sol[i+j*Nx] << std::endl;
  //       }
  //   }
}

void TimeScheme::solve()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Time loop..." << std::endl;
  
  // Variables pratiques
  int n(0);
  int scenario(_DF->getScenario());

  // Sauvegarde la condition initiale
  std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(n) + ".vtk");
  std::string exactSolFileName(_resultsDir + "/solution_exacte_scenario_" + std::to_string(scenario) + "_" + std::to_string(n) + ".vtk");
  saveCurrentSolution(solFileName);
  _function->buildExactSolution(_currentTime);
  _function->saveCurrentExactSolution(exactSolFileName);
             
  // Démarrage du chrono
  auto start = std::chrono::high_resolution_clock::now();

  // Boucle en temps
  while (_currentTime < _finalTime)
    {
      oneStep();
      ++n;
      _currentTime += _timeStep;
      _function->buildExactSolution(_currentTime);
      if (n % _DF->getSaveFrequency() == 0)
        {
          // Save numerical solution
          std::cout << "Saving solution at t = " << _currentTime << std::endl;
          std::string solFileName(_resultsDir + "/solution_scenario_" + std::to_string(scenario) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".vtk");
          saveCurrentSolution(solFileName);
          // Save exact solution
          if (_DF->getScenario() == 1 || _DF->getScenario() == 2)
            {
              std::cout << "Saving exact solution at t = " << _currentTime << std::endl;
              std::string exactSolFileName(_resultsDir + "/solution_exacte_scenario_" + std::to_string(scenario) + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".vtk");
              _function->saveCurrentExactSolution(exactSolFileName);
            }
        }
    }

  // Fin du chrono
  auto finish = std::chrono::high_resolution_clock::now();
  double duration = std::chrono::duration<double, std::milli>(finish-start).count() * 1e-3;
  
  // Calcul de l'erreur pour les scenario 1 et 2
  if (_DF->getScenario() == 1 || _DF->getScenario() == 2)
    {
      double error(computeCurrentError());
      std::cout << "L2 error = " << error << " at t = " << _currentTime << " for Nx = " << _DF->getNx() << ", Ny = " << _DF->getNy() << std::endl;
    }
  // Logs de fin
  std::cout << termcolor::green << "SUCCESS::TIMESCHEME : Time loop completed successfully in " << duration << " seconds !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

double TimeScheme::computeCurrentError()
{
  double error(0.);
  DVector errorVec(_Sol - _function->getExactSolution());
  error = errorVec.dot(errorVec);
  error = sqrt(_DF->getDx() * _DF->getDy() * error);
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
  // Calcul de la solution
  _Sol = _Sol + (_laplacian->matVecProd(_Sol) + dt * _function->getSourceTerm());
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
  resFile.precision(7);
  // Calcul du terme source
  _function->buildSourceTerm(_currentTime + dt);
  // Calcul de la solution
  _Sol = _laplacian->solveConjGrad(_Sol + dt * _function->getSourceTerm(), _Sol, tolerance, maxIt, resFile);
}
