#include "Vector.h"
#include <vector>
#include <iostream>

// Constructors
DVector::DVector():
  std::vector<double>()
{
}


// Constructors
DVector::DVector(size_type count):
  std::vector<double>(count)
{
}

// Additional methods
DVector DVector::add(const DVector& vec)
{
  int N(vec.size());
  for (int i(0) ; i < N ; ++i)
    {
      this->operator[](i) += vec[i];
    }
  return *this;
}


DVector DVector::sub(const DVector& vec)
{
  int N(vec.size());
  for (int i(0) ; i < N ; ++i)
    {
      this->operator[](i) -= vec[i];
    }
  return *this;
}


double DVector::dot(const DVector& vec)
{
  int N(vec.size());
  double dot(0.);
  for (int i(0) ; i < N ; ++i)
    {
      dot += this->operator[](i) * vec[i];
    }
  return dot;
}


// Operators
std::ostream& operator<< (std::ostream &os, const DVector& v)
{
  int N(v.size());
  for (int i(0) ; i < N - 1 ; ++i)
    {
      os << v[i] << " ";
    }
  os << v[N-1] << std::endl;
  return os;
}

DVector operator+ (const DVector& u, const DVector& v)
{
  if (u.size() == v.size())
    {
      DVector w;
      w.resize(u.size());
      for (int i(0) ; i < u.size() ; ++i)
        {
          w[i] = u[i] + v[i];
        }
      return w;
    }
  else
    {
      std::cout << "ERROR : DVector sizes do not match (" << u.size() << " and " << v.size() << ")" << std::endl;
      exit(1);
    }
}

DVector operator- (const DVector& u, const DVector& v)
{
  if (u.size() == v.size())
    {
      DVector w;
      w.resize(u.size());
      for (int i(0) ; i < u.size() ; ++i)
        {
          w[i] = u[i] - v[i];
        }
      return w;
    }
  else
    {
      std::cout << "ERROR : DVector sizes do not match (" << u.size() << " and " << v.size() << ")" << std::endl;
      exit(1);
    }
}


DVector operator* (double alpha, const DVector& u)
{
  DVector w;
  w.resize(u.size());
  for (int i(0) ; i < u.size() ; ++i)
    {
      w[i] = alpha * u[i];
    }
  return w;
}


DVector operator* (const DVector& u, double alpha)
{
  DVector w;
  w.resize(u.size());
  for (int i(0) ; i < u.size() ; ++i)
    {
      w[i] = alpha * u[i];
    }
  return w;
}
