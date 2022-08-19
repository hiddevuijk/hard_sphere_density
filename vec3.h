/*

  3 dimensional vector class.
  operators (+,-, +=, -=, *, / ...) are defined element wise

*/



#ifndef GUARD_VEC3_H
#define GUARD_VEC3_H

#include <math.h>
#include <iostream>

class Vec3 {
 public:
  Vec3(double x, double y, double z)
    : x(x), y(y), z(z) {}

  Vec3()
    : x(0.0), y(0.0), z(0.0) {}
  
  double x, y, z;

  double Length() const {
    return sqrt(x * x + y * y + z * z);
  }
  double LengthSquared() const {
    return x * x + y * y + z * z;
  }

  // set length to d
  void Normalize(double d=1.0) {
    double l = Length();
    x *= d / l;
    y *= d / l;
    z *= d / l;
  }

  // apply periodic boundary conditions for a square box with size L 
  // and the origin at the bottom left corner
  void PeriodicBoundaryCondition(double L) {
    x -= L * floor(x / L);
    y -= L * floor(y / L);
    z -= L * floor(z / L);
  }
    

  // apply periodic boundary condition in one direction
  void PeriodicBoundaryX(double L) {
    x -= L * floor(x / L);
  }

  void PeriodicBoundaryY(double L) {
    y -= L * floor(y / L);
  }

  void PeriodicBoundaryZ(double L) {
    z -= L * floor(z / L);
  }

  // operators
  Vec3 operator+=(const Vec3& r) {
    x += r.x; y += r.y; z += r.z;
    return *this;
  }

  Vec3 operator+=(const double& add) {
    x += add; y += add; z += add;
    return *this;
  }

  Vec3 operator-=(const Vec3& r) {
    x -= r.x; y -= r.y; z -= r.z;
    return *this;
  }  

  Vec3 operator-=(const double& minus) {
    x -= minus; y -= minus; z -= minus;
    return *this;
  }

  Vec3 operator*=(const Vec3& r) {
    x *= r.x; y *= r.y; z *= r.z;
    return *this;
  }

  Vec3 operator*=(const double& mult) {
    x *= mult; y *= mult; z *= mult;
    return *this;
  }

  Vec3 operator/=(const Vec3& r) {
    x /= r.x; y /= r.y; z /= r.z;
    return *this;
  }

  Vec3 operator/=(const double& div) {
    x /= div; y /= div; z /= div;
    return *this;
  }

};


////////////////////////////////////////
//
// Nonmember functions
//
////////////////////////////////////////

#endif
