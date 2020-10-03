// primitives.cpp

#include "primitives.h"

// Point
Point::Point() { // Default Constructor
	x = 0.0; y = 0.0, z = 0.0;
}

Point::Point(double newX, double newY, double newZ){ // Constructor
    x = newX; y = newY, z = newZ;
}

// Color
Color::Color() { // Default Constructor
	r = 0.0; g = 0.0; b = 0.0;
}

Color::Color(double newR, double newG, double newB){ // Constructor
    r = newR; g = newG; b = newB;
}

// Material
Material::Material() { // Default Constructor
	Kdiff = new Color();
	Kspec = new Color();
	SpecExp = 0.0; Kreflect = 0.0; Krefract = 0.0; IndRefract = 0.0;
}

Material::Material(Color* diffuse, Color* specular, double exponent, double reflect, double refract, double index) { // Constructor
	Kdiff = diffuse;
	Kspec = specular;
	SpecExp = exponent;
	Kreflect = reflect;
	Krefract = refract;
	IndRefract = index;
}

// Sphere
Sphere::Sphere(){ // Default Constructor
	center = NULL; material = NULL; radius = 0.0;
}

Sphere::Sphere(Point* ctr, double rad, Material* mat, vector<Sphere*>* sphereList) { // Constructor
	center = ctr;
	radius = rad;
	material = mat;
	sphereList->push_back(this);
}

Point* Sphere::getCenter(){ // sphere center point getter
	return center;
}

Material* Sphere::getMaterial(){ // sphere material getter
	return material;
}

vector<double> Sphere::findIntersect(Point* rayStart, vector<double> rayVect) { // Vector sphere intersect
	
	// finding roots from the equation:
	// t^2(Vx^2 + Vy^2 + Vz^2) + t(2QxVx + 2QyVy + 2QzVz) + (Qx^2 + Qy^2 + Qz^2) - r^2 = 0 

	double Qx = rayStart->x - center->x;	double Qy = rayStart->y - center->y;	double Qz = rayStart->z - center->z;
	double Vx = rayVect[0];	double Vy = rayVect[1];	double Vz = rayVect[2];
	
	double t2Constant = (Vx * Vx) + (Vy * Vy) + (Vz * Vz);
	double t1Constant = (2 * Qx * Vx) + (2 * Qy * Vy) + (2 * Qz * Vz); 
	double t0Constant = ((Qx * Qx) + (Qy * Qy) + (Qz * Qz)) - (radius * radius);

	if (((t1Constant * t1Constant) - (4 * t2Constant * t0Constant)) < 0) {
		return {-1.0, -1.0};
	}

	double root1 = (-1 * t1Constant + sqrt((t1Constant * t1Constant) - (4 * t2Constant * t0Constant))) / (2.0 * t2Constant);
	double root2 = (-1 * t1Constant - sqrt((t1Constant * t1Constant) - (4 * t2Constant * t0Constant))) / (2.0 * t2Constant);

	return {root1, root2};
}

// Plane
Plane::Plane() { // Default Constructor
	A = 0.0; B = 0.0; C = 0.0; D = 0.0;
	material = new Material();
}

Plane::Plane(double a, double b, double c, double d, Material* mat, vector<Plane*>* planeList){ // Constructor
	A = a;
	B = b;
	C = c;
	D = d;
	material = mat;
	planeList->push_back(this);
}

vector<double> Plane::getNormal(){ // Plane normal getter
	return {A, B, C};
}

Material* Plane::getMaterial(){ // Plane Material Getter
	return material;
}

double Plane::findIntersect(Point* rayStart, vector<double> rayVect) { // finding intersect of a vector and a plane

	if (rayVect.size() != 3)
		return 0.0;

	double P0x = rayStart->x;	double P0y = rayStart->y;	double P0z = rayStart->z;
	double Vx = rayVect[0];	double Vy = rayVect[1];	double Vz = rayVect[2];

	if (((A * Vx) + (B * Vy) + (C * Vz)) == 0.0) {
		return -1.0;
	}

	double tIntersect = -1 * ((A * P0x) + (B * P0y) + (C * P0z) + D) / ((A * Vx) + (B * Vy) + (C * Vz));

	if (tIntersect < 0.0)
		return -1.0;

	return tIntersect;

}