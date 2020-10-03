// ray.cpp

#include "ray.h"
#include <iostream>

Ray::Ray(){ // Default constructor
    startPoint = new Point(0.0, 0.0, 0.0);
    rayVector = {0.0f, 0.0f, -1.0f};
	currRefractedIndex = 1.0;
}

Ray::Ray(double x, double y, int windowSize, double fov) { // Starting Ray Constructor

    startPoint = new Point(0.0, 0.0, 0.0);
	currRefractedIndex = 1.0;
    double Vx = (x / (windowSize/2.0)) * tan((fov / 2.0) * PI / 180.0);
    double Vy = (y / (windowSize / 2.0)) * tan((fov / 2.0) * PI / 180.0);
    
	rayVector = {Vx, Vy, -1.0f};
}

Ray::Ray(vector<double> vec, Point* start) { // Intermediate Ray Constructor
	rayVector = vec;
	startPoint = start;
	currRefractedIndex = 1.0;
}

Ray::Ray(vector<double> vec, Point* start, double refIndex) { // refracted Ray Constructor
	rayVector = vec;
	startPoint = start;
	currRefractedIndex = refIndex;
}

Point* Ray::getStart() { // start point getter
	return startPoint;
}

vector<double> Ray::getVec(){ // ray direction getter
    return rayVector;
}

double Ray::getIndex() { // refractive index getter
	return currRefractedIndex;
}

