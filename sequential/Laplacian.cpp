#include "Laplacian.h"
#include "Function.h"
#include "DataFile.h"
#include "Vector.h"

#include <cmath>

Laplacian::Laplacian()
{
}


Laplacian::Laplacian(DataFile* DF, Function* function):
  _DF(DF), _function(function)
{
}


void Laplacian::Initialize(DataFile* DF, Function* function)
{
  _DF = DF;
  _function = function;
  this->Initialize();
}


void Laplacian::Initialize()
{
  // On récupère les paramètres necessaires pour construire la matrice
  double dx(_DF->getDx()), dy(_DF->getDy());
  double D(_DF->getDiffCoeff());
  double dt(_DF->getTimeStep());
  _Nx = _DF->getNx();
  _Ny = _DF->getNy();

  // Calcul des coefficients de la matrice
  if (_DF->getTimeScheme() == "ExplicitEuler")
    {
      _alpha = dt * D / pow(dy,2);
      _beta = dt * D / pow(dx,2);
      _gamma = - 2. * dt * D * (1./pow(dx,2) + 1./pow(dy,2));
    }
  else if (_DF->getTimeScheme() == "ImplicitEuler")
    {
      _alpha = - dt * D / pow(dy,2);
      _beta = - dt * D / pow(dx,2);
      _gamma = 1 + 2. * dt * D * (1./pow(dx,2) + 1./pow(dy,2));
    }
}


DVector Laplacian::matVecProd(const DVector& x)
{
  // Vecteur resultat
  DVector result;
  int size(x.size());
  result.resize(size, 0.);

  for (int i(0) ; i < size ; ++i)
    {
      result[i] = _gamma * x[i];
      if (i % _Nx != 0)
        {
          result[i] += _beta * x[i-1];
        }
      if (i % _Nx != _Nx - 1)
        {
          result[i] += _beta * x[i+1];
        }
      if (i < _Nx * _Ny - _Nx)
        {
          result[i] += _alpha * x[i+_Nx];
        }
      if (i >= _Nx)
        {
          result[i] += _alpha*x[i-_Nx];
        }
    }
  return result;
}


DVector Laplacian::solveConjGrad(const DVector& b, const DVector& x0, double tolerance, int maxIterations, std::ofstream& resFile)
{
  // Variables intermédiaires
  DVector x(x0);
  DVector res(b - this->matVecProd(x0));
  DVector p(res);
  double beta(sqrt(res.dot(res)));
  resFile << beta << std::endl;
  
  // Itérations de la méthode
  int k(0);
  while ((beta > tolerance) && (k < maxIterations))
    {
      DVector z(this->matVecProd(p));
      double alpha(res.dot(p)/z.dot(p));
      x = x + alpha * p;
      res = res - alpha * z;
      double gamma(res.dot(res)/pow(beta,2));
      p = res + gamma * p;
      beta = sqrt(res.dot(res));
      ++k;
      resFile << beta << std::endl;
    }
  // Logs
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
  return x;
}
