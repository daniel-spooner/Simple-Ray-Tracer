// light.h

#include <vector>
#include "primitives.h"

class Light { // Light Class
private:
	Point* point;
	Color* intensity;

public:
	Light();
	Light(Point* pnt, Color* clr, vector<Light*>* lightList);

	Point* getPoint();
	Color* getIntensity();
};
