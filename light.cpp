// light.cpp

#include "light.h"

Light::Light() { // Default Constructor
	point = new Point();
	intensity = new Color();
}

Light::Light(Point* pnt, Color* clr, vector<Light*>* lightList) { // Constructor
	point = pnt;
	intensity = clr;
	lightList->push_back(this);
}

Point* Light::getPoint() { // Light Point getter
	return point;
}

Color* Light::getIntensity() { // Light intensity getter
	return intensity;
}