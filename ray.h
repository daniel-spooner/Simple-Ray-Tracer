// ray.h

#pragma once

#include <GL/glut.h>
#include <vector>
#include <math.h>
#include <iostream>

#include "primitives.h"

using namespace std;

// Ray Class
class Ray{

private:
    Point* startPoint; // Starting point vector
    vector<double> rayVector; // Ray directional vector
	double currRefractedIndex;

public:
    Ray();
    Ray(double x, double y, int windowSize, double fov);
	Ray(vector<double> vec, Point* start);
    Ray(vector<double> vec, Point* start, double refIndex);

	Point* getStart();
	vector<double> getVec();
	double getIndex();

};
