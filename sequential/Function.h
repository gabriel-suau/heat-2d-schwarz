#ifndef FUNCTION_H
#define FUNCTION_H

#include "DataFile.h"
#include "termcolor.h"
#include "Vector.h"

class Function
{
protected:
  // Pointeur vers le fichier de paramètres
  DataFile* _DF;
  
  // Variables utiles
  double _xmin, _ymin, _xmax, _ymax, _Lx, _Ly, _dx, _dy;
  int _Nx, _Ny;

  // Vecteur solution initiale
  DVector _Sol0;
  // Terme source
  DVector _sourceTerm;
  // Solution exacte
  DVector _exactSol;
  
public:
  // Constructeur
  Function();
  Function(DataFile* DF);

  // Destructeur
  virtual ~Function() = default;
  
  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF);

  // Construction du terme source
  void buildSourceTerm(double t);

  // Construction de la solution exacte
  virtual void buildExactSolution(double t) = 0;
  
  // Sauvegarde de la solution exacte
  void saveCurrentExactSolution(std::string& fileName) const;

  // Getters
  const DVector& getInitialCondition() const {return _Sol0;};
  const DVector& getSourceTerm() const {return _sourceTerm;};
  const DVector& getExactSolution() const {return _exactSol;};
  
  // Pour construire les CI/CL en fonction du système
  virtual double f(const double x, const double y, const double t) = 0;
  virtual double g(const double x, const double y, const double t) = 0;
  virtual double h(const double x, const double y, const double t) = 0;
};

// Classe fille scenario 1
class Function1 : public Function
{
public:
  // Constructeur
  Function1();
  Function1(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

// Classe fille scenario 2
class Function2 : public Function
{
public:
  // Constructeur
  Function2();
  Function2(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

// Classe fille scenario 3
class Function3 : public Function
{
public:
  // Constructeur
  Function3();
  Function3(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

#endif // FUNCTION_H
