#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "DataFile.h"
#include "Function.h"
#include "Vector.h"
#include <fstream>

class Laplacian
{
private:
  // Pointeurs vers les trucs importants
  DataFile* _DF;
  Function* _function;
  
  // Coefficients de la matrice du laplacien 2D
  double _alpha, _beta, _gamma;
  int _Nx, _Ny;
  
public:
  // Constructeurs
  Laplacian();
  Laplacian(DataFile* DF, Function* function);

  // Destructeur
  ~Laplacian() = default;

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Function* function);
  
  // Getters
  double getAlpha() const {return _alpha;};
  double getBeta() const {return _beta;};
  double getGamma() const {return _gamma;};
  int getNx() const {return _Nx;}
  int getNy() const {return _Ny;};
  
  // Produit matVec et GC
  DVector matVecProd(const DVector& x);
  DVector solveConjGrad(const DVector& b, const DVector& x0, double tolerance, int maxIterations, std::ofstream& resFile);
  
  // Printer (pour debugger)
  void print() const;
};

#endif // LAPLACIAN_H
