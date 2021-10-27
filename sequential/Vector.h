#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

// Classe de vecteur pratique pour les opérations élémentaires
class DVector: public std::vector<double>
{
public:
  // Constructors
  DVector();
  DVector(size_type count);

  // Additional methods
  DVector add(const DVector& vec);
  DVector sub(const DVector& vec);
  double dot(const DVector& vec);
};

std::ostream& operator<< (std::ostream &os, const DVector& v);
DVector operator+ (const DVector& u, const DVector& v);
DVector operator- (const DVector& u, const DVector& v);
DVector operator* (double alpha, const DVector& u);
DVector operator* (const DVector& u, double alpha);

#endif //VECTOR_H
