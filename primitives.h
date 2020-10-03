// primitives.h

#pragma once

#include <vector>
#include <math.h>

using namespace std;

const double PI = 3.14159265; // Pi constant

class Point { // Point Class
public:
    double x;
    double y;
    double z;
	Point();
    Point(double newX, double newY, double newZ);
};

class Color { // Color Class
public:
    double r;
    double g;
    double b;
	Color();
    Color(double newR, double newG, double newB);
};

class Material {
public:
	Color* Kdiff;		// Diffuse reflection coefficient
	Color* Kspec;		// Specular reflection coefficient
	double SpecExp;		// Specular exponent
	double Kreflect;	// Reflected light coefficient
	double Krefract;	// Refracted light coefficient
	double IndRefract;  // Index of refraction
	
	Material();
	Material(Color* diffuse, Color* specular, double exponent, double reflect, double refract, double index);
};

class Sphere{ // Sphere Class

private:
	Material* material;
	Point* center;
	double radius;

public:
    Sphere();
	Sphere(Point* ctr, double rad, Material* mat, vector<Sphere*>* sphereList);

	Point* getCenter();
	Material* getMaterial();
	vector<double> findIntersect(Point* rayStart, vector<double> rayVect);
};

class Plane{ // Plane Class

private:
	Material* material;
	double A;
	double B;
	double C;
	double D;

public:
	Plane();
	Plane(double a, double b, double c, double d, Material* mat, vector<Plane*>* planeList);

	vector<double> getNormal();
	Material* getMaterial();
	double findIntersect(Point* rayStart, vector<double> rayVect);
};
