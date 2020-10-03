/*
Ray Tracer

Author: Daniel Spooner

*/
// main.cpp

#include <GL/glut.h>
#include <iostream>
#include <vector>

#include "ray.h"
#include "primitives.h"
#include "light.h"


// Global Variables
const unsigned int windowSize = 400;
const unsigned int depth = 4;
const double FOVy = 60.0;
const double FOVx = FOVy;
int sceneNumber;
Color* ambient = new Color(0.05, 0.05, 0.05); // Ambient render Color

// Global Object Lists
vector<Light*>* lightList = new vector<Light*>();
vector<Sphere*>* sphereList = new vector<Sphere*>();
vector<Plane*>* planeList = new vector<Plane*>();

void printVector(vector<double> vec) { // print all items in a vector
	cout << "[ ";
	for (int i = 0; i < vec.size(); i++) {
		cout << vec[i] << " ";
	}
	cout << "]" << endl;
}

double dotProd(vector<double> v1, vector<double> v2){ // taking dot product of two vectors
	if(v1.size() != 3 && v2.size() != 3)
		return 0.0;
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

vector<double> normalize(vector<double> vec){ // normalizing a vector
	vector<double> newVec = {0.0, 0.0, 0.0};
	
	if(vec.size() != 3)
		return vec;
	
	double size = sqrt(dotProd(vec, vec));
	newVec[0] = vec[0] / size;
	newVec[1] = vec[1] / size;
	newVec[2] = vec[2] / size;
	
	return newVec;
}

vector<double> subtractVectors(vector<double> v1, vector<double> v2) { // subtracting two vectors
	if (v1.size() != v2.size())
		return v1;

	for (int i = 0; i < v1.size(); i++) {
		v1[i] = v1[i] - v2[i];
	}
	return v1;
}

vector<double> addVectors(vector<double> v1, vector<double> v2) { // adding two vectors
	if (v1.size() != v2.size())
		return v1;

	for (int i = 0; i < v1.size(); i++) {
		v1[i] = v1[i] + v2[i];
	}
	return v1;
}

vector<double> mulScalarVector(double scalar, vector<double> vec) { // multiply scalar and vector
	for (int i = 0; i < vec.size(); i++)
		vec[i] *= scalar;
	return vec;
}

double properRoot(vector<double> roots) { // finding correct root of sphere intersection
	double smallestRoot = roots[0];

	if(roots[0] < roots[1] && roots[0] >= 0.01){
		return roots[0];
	}else if(roots[1] < roots[0] && roots[1] >= 0.01){
		return roots[1];
	}else {
		return -1.0;
	}

	return -1.0;
}

// Trace function: traces rays into scene, recursively calls reflected and refracted rays
Color* trace(Ray* ray, unsigned int rayDepth) {

	if (rayDepth == 0) { // Returning ambient light if ray depth is maxed 
		return ambient;
	}

	int i;
	double minT = -1.0;

	double sphereT;
	double minSphereT = -1.0;
	Sphere* minSphere = NULL;
	for (i = 0; i < sphereList->size(); i++) { // finding smallest t-value of all sphere intersects
		sphereT = properRoot((*sphereList)[i]->findIntersect(ray->getStart(), ray->getVec()));
		if ((minSphereT < 0.01 && sphereT > 0.01) || (sphereT < minSphereT && sphereT > 0.01)) {
			minSphereT = sphereT;
			minSphere = (*sphereList)[i];
		}
	}

	double planeT;
	double minPlaneT = -1.0;
	Plane* minPlane = NULL;
	for (i = 0; i < planeList->size(); i++) { // finding smallest t-value of all plane intersects
		planeT = (*planeList)[i]->findIntersect(ray->getStart(), ray->getVec());
		if (planeT > 0.01) {
			if (minPlaneT < 0.0 || planeT < minPlaneT) {
				minPlaneT = planeT;
				minPlane = (*planeList)[i];
			}
		}
	}

	// finding minimum t-value from sphere min and plane min
	if (minSphereT < 0.01 && minPlaneT < 0.01) {
		return ambient;
	}
	else if (minSphereT < 0.01) {
		minT = minPlaneT;
		minSphere = NULL;
	}
	else if (minPlaneT < 0.01) {
		minT = minSphereT;
		minPlane = NULL;
	}
	else if (minSphereT < minPlaneT) {
		minT = minSphereT;
		minPlane = NULL;
	}
	else {
		minT = minPlaneT;
		minSphere = NULL;
	}

	vector<double> rayVec = ray->getVec();
	Point* rayStart = ray->getStart();
	Point* intersect = new Point(rayStart->x + minT * rayVec[0],
		rayStart->y + minT * rayVec[1],
		rayStart->z + minT * rayVec[2]);

	// VIEW VECTOR
	vector<double> viewVec = normalize(rayVec);
	viewVec[0] *= -1; viewVec[1] *= -1; viewVec[2] *= -1;

	// NORMAL VECTOR
	vector<double> normalVec;
	if (minSphere == NULL) {
		normalVec = normalize(minPlane->getNormal());
	}
	else {
		Point* sphrCenter = minSphere->getCenter();
		normalVec = normalize({ intersect->x - sphrCenter->x, intersect->y - sphrCenter->y, intersect->z - sphrCenter->z });
	}

	Material* mat; // getting material
	if (minSphere == NULL)
		mat = minPlane->getMaterial();
	else
		mat = minSphere->getMaterial();

	vector<double> lightVec;
	vector<double> reflectVec;
	vector<double> shadowVec;
	Color* finalIntensity = new Color(0.0, 0.0, 0.0);
	Color* summationIntensity = new Color(0.0, 0.0, 0.0);
	for (i = 0; i < lightList->size(); i++) { // summation in phong lighting equation
		// LIGHT VECTOR
		lightVec = { 0.0, 0.0, 0.0 };
		Point* light = (*lightList)[i]->getPoint();
		if (lightList->size() > 0) {
			lightVec = normalize({ light->x - intersect->x, light->y - intersect->y, light->z - intersect->z });
		}

		// SHADOW VECTOR
		shadowVec = { intersect->x - light->x, intersect->y - light->y, intersect->z - light->z };
		double Si = 1.0;
		double currLightT;
		for(int j = 0; j < sphereList->size(); j ++){
			currLightT = properRoot((*sphereList)[j]->findIntersect((*lightList)[i]->getPoint(), shadowVec));
			if (currLightT > 0.0 && currLightT < 0.95) {
				Si = 0.0 + (*sphereList)[j]->getMaterial()->Krefract;
			}
		}
		
		// REFLECTED LIGHT VECTOR
		reflectVec = normalize(subtractVectors(mulScalarVector(2.0 * dotProd(normalVec, lightVec), normalVec), lightVec));

		double RiV = dotProd(reflectVec, viewVec);
		if (RiV < 0.0 || mat->SpecExp == 0.0)
			RiV = 0.0;
		double RiVq = pow(RiV, mat->SpecExp);
		Color* KsRiVq = new Color(mat->Kspec->r * RiVq, mat->Kspec->g * RiVq, mat->Kspec->b * RiVq);

		double NLi = dotProd(normalVec, lightVec);
		if (NLi < 0.0)
			NLi = 0.0;
		Color* KdNLi = new Color(mat->Kdiff->r * NLi, mat->Kdiff->g * NLi, mat->Kdiff->b * NLi);

		Color* lightClr = (*lightList)[i]->getIntensity();
		// summing light intensity
		summationIntensity->r += Si * lightClr->r * (KsRiVq->r + KdNLi->r);
		summationIntensity->g += Si * lightClr->g * (KsRiVq->g + KdNLi->g);
		summationIntensity->b += Si * lightClr->b * (KsRiVq->b + KdNLi->b);
	}

	Color* IaKd = new Color(ambient->r * mat->Kdiff->r, ambient->g * mat->Kdiff->g, ambient->b * mat->Kdiff->b);
	
	finalIntensity->r = IaKd->r + summationIntensity->r;
	finalIntensity->g = IaKd->g + summationIntensity->g;
	finalIntensity->b = IaKd->b + summationIntensity->b;


	// if there is reflection, trace reflected ray to find reflected color
	if (mat->Kreflect > 0.0) {
		vector<double> newRayVec = normalize(subtractVectors(mulScalarVector(2.0 * dotProd(normalVec, viewVec), normalVec), viewVec));
		Color* reflectedColor;
		reflectedColor = trace(new Ray(newRayVec, intersect), rayDepth - 1);

		finalIntensity->r += (mat->Kreflect * reflectedColor->r);
		finalIntensity->g += (mat->Kreflect * reflectedColor->g);
		finalIntensity->b += (mat->Kreflect * reflectedColor->b);
	}
	else {
		finalIntensity->r = IaKd->r + summationIntensity->r; finalIntensity->g = IaKd->g + summationIntensity->g; finalIntensity->b = IaKd->b + summationIntensity->b;
	}

	// if there is refraction, trace refracted ray to find refracted color
	if(mat->Krefract > 0.0 && minSphere != NULL){	
		
		vector<double> newRayVec;
		Color* refractedColor;
		double n1 = ray->getIndex();
		double n2 = mat->IndRefract;
			
		double n = n1 / n2;

		double cos1 = -1 * dotProd(normalVec, normalize(ray->getVec()));
		double sin2 = n * n * (1.0 - cos1 * cos1);

		if (sin2 > 1.0){
			vector<double> reflect = normalize(subtractVectors(mulScalarVector(2.0 * dotProd(normalVec, viewVec), normalVec), viewVec));
			refractedColor = trace (new Ray(reflect, intersect, n1), rayDepth - 1);
		}else{
			double cost = sqrt(1.0 - sin2);
			newRayVec = addVectors(mulScalarVector(n, ray->getVec()), mulScalarVector(n * cos1 - cost, normalVec));
			refractedColor = trace (new Ray(newRayVec, intersect, n2), rayDepth - 1);
		}

		if(refractedColor != NULL){
			finalIntensity->r += (mat->Krefract * refractedColor->r);
			finalIntensity->g += (mat->Krefract * refractedColor->g);
			finalIntensity->b += (mat->Krefract * refractedColor->b);
		}
	}
	
	// returning final intensity
	return finalIntensity;
}


void render() { // render - casts rays into scene
	double lowerBound = -1.0 * (double)windowSize / 2.0 + 0.5;
	double upperBound = (double)windowSize / 2.0;

	Ray* r = NULL;

	double currT;
	Color* clr = NULL;

	// Nested loop casting a ray into each pixel
	for (double x = lowerBound; x < upperBound; x += 1.0) {
		for (double y = lowerBound; y < upperBound; y += 1.0) {
			// tracing a ray and setting retrieved pixel color
			clr = trace(new Ray(x, y, windowSize, FOVy), depth);
			glColor3f((GLfloat)clr->r, (GLfloat)clr->g, (GLfloat)clr->b);
			glBegin(GL_POINTS);
				glVertex2f((int)(x - 0.5) + (windowSize / 2.0), (int)(y - 0.5) + (windowSize / 2.0));
			glEnd();
		}
	}
}

void scene1() { // scene 1
	
	Material* planeColor = new Material(new Color(0.05, 0.05, 0.05), new Color(0.1, 0.1, 0.1), 1.0, 0.0, 0.0, 0.0);

	Material* diffuse1 = new Material(new Color(0.8, 0.8, 0.2), new Color(0.0, 0.0, 0.0), 1.0, 0.0, 0.0, 1.0);
	Material* diffuse2 = new Material(new Color(0.2, 0.8, 0.8), new Color(0.0, 0.0, 0.0), 1.0, 0.0, 0.0, 1.0);
	Material* diffuse3 = new Material(new Color(0.8, 0.2, 0.8), new Color(0.0, 0.0, 0.0), 1.0, 0.0, 0.0, 1.0);
	
	Material* specular1 = new Material(new Color(0.5, 0.5, 0.8), new Color(0.8, 0.8, 0.8), 2.0, 0.0, 0.0, 1.0);
	Material* specular2 = new Material(new Color(0.5, 0.5, 0.8), new Color(0.8, 0.8, 0.8), 8.0, 0.0, 0.0, 1.0);
	Material* specular3 = new Material(new Color(0.5, 0.5, 0.8), new Color(0.8, 0.8, 0.8), 64.0, 0.0, 0.0, 1.0);

	new Plane(0.0, 0.0, 1.0, 15.0, planeColor, planeList);

	new Light(new Point(0.0, 0.0, -5.0), new Color(1.0, 1.0, 1.0), lightList);

	new Sphere(new Point(-3.5, 2.5, -10.0), 1.5, diffuse1, sphereList);
	new Sphere(new Point(0, 2.5, -10.0), 1.5, diffuse2, sphereList);
	new Sphere(new Point(3.5, 2.5, -10.0), 1.5, diffuse3, sphereList);

	new Sphere(new Point(-3.5, -2.5, -10.0), 1.5, specular1, sphereList);
	new Sphere(new Point(0, -2.5, -10.0), 1.5, specular2, sphereList);
	new Sphere(new Point(3.5, -2.5, -10.0), 1.5, specular3, sphereList);

	render();
	glFlush();
}

void scene2() { // scene 2

	
	Material* planeColor = new Material(new Color(0.15, 0.15, 0.15), new Color(0.1, 0.1, 0.1), 1.0, 0.5, 0.0, 1.0);
	Material* sphereColor1 = new Material(new Color(0.2, 0.2, 0.8), new Color(0.5, 0.5, 0.5), 64.0, 1.0, 0.0, 1.0);
	Material* sphereColor2 = new Material(new Color(0.8, 0.2, 0.2), new Color(0.5, 0.5, 0.5), 64.0, 1.0, 0.0, 1.0);
	Material* sphereColor3 = new Material(new Color(0.8, 0.4, 0.1), new Color(0.5, 0.5, 0.5), 64.0, 1.0, 0.0, 1.0);

	new Plane(0.0, 0.0, 1.0, 8.0, planeColor, planeList);

	new Sphere(new Point(0.0, 1.5, -5.0), 1.0, sphereColor1, sphereList);
	new Sphere(new Point(0.0, -1.5, -5.0), 1.0, sphereColor2, sphereList);
	new Sphere(new Point(-2.0, 0, -5.0), 0.75, sphereColor3, sphereList);

	new Light(new Point(-2.0, -2.0, -3.0), new Color(0.4, 0.4, 0.9), lightList);
	new Light(new Point(-2.0, 2.0, -3.0), new Color(0.9, 0.4, 0.4), lightList);
	

	render();
	glFlush();
}

void scene3() { // scene 3
	
	Material* planeColor = new Material(new Color(0.65, 0.65, 0.75), new Color(0.1, 0.1, 0.1), 1.0, 0.0, 0.0, 1.0);
	Material* sphereColor1 = new Material(new Color(0.0, 0.3, 1.0), new Color(1.0, 1.0, 1.0), 64.0, 0.1, 0.0, 1.0);
	Material* sphereColor2 = new Material(new Color(0.8, 0.8, 0.2), new Color(0.0, 0.0, 0.0), 1.0, 0.1, 0.0, 1.0);
	
	Material* refractR = new Material(new Color(0.2, 0.0, 0.0), new Color(1.0, 1.0, 1.0), 64.0, 0.0, 0.85, 0.95);
	Material* refractG = new Material(new Color(0.0, 0.2, 0.0), new Color(1.0, 1.0, 1.0), 64.0, 0.0, 0.85, 1.05);
	Material* refractY = new Material(new Color(0.2, 0.2, 0.0), new Color(1.0, 1.0, 1.0), 64.0, 0.0, 0.85, 0.95);

	new Plane(0.0, 1.0, 0.0, 8.0, planeColor, planeList);

	new Sphere(new Point(-1.25, -0.75, -5.0), 1.5, refractR, sphereList);
	new Sphere(new Point(1.25, -0.4, -4.0), 0.75, refractG, sphereList);
	new Sphere(new Point(1.25, 2.0, -6.0), 1.0, refractY, sphereList);
	new Sphere(new Point(2.5, -2.0, -7.0), 1.0, sphereColor2, sphereList);
	new Sphere(new Point(0.0, 6.0, -20.0), 10.0, sphereColor1, sphereList);

	new Light(new Point(10.0, 0.0, 0.0), new Color(0.8, 0.8, 0.8), lightList);
	new Light(new Point(0.0, 0.0, 0.0), new Color(0.8, 0.8, 0.8), lightList);

	render();
	glFlush();
}

void scene4() { // scene 4

	Material* planeColor1 = new Material(new Color(0.3, 0.3, 0.1), new Color(0.8, 0.8, 0.8), 32.0, 0.9, 0.0, 0.0);

	Material* diffuseColor = new Material(new Color(0.3, 0.9, 0.3), new Color(0.0, 0.0, 0.0), 1.0, 0.0, 0.0, 0.0);
	Material* specularColor = new Material(new Color(0.8, 0.8, 0.2), new Color(1.0, 1.0, 1.0), 64.0, 0.0, 0.0, 0.0);

	Material* reflect = new Material(new Color(1.0, 1.0, 1.0), new Color(0.8, 0.8, 0.8), 64.0, 0.9, 0.0, 0.0);
	Material* refract = new Material(new Color(0.1, 0.0, 0.0), new Color(0.8, 0.8, 0.8), 64.0, 0.0, 0.90, 0.95);

	new Plane(0.0, 1.0, 0.0, 2.0, planeColor1, planeList);
	new Plane(-1.0, 0.0, 1.0, 50.0, planeColor1, planeList);

	new Sphere(new Point(-3.5, 0.0, -10.0), 1.0, reflect, sphereList);
	new Sphere(new Point(-2.0, 2.0, -10.0), 1.0, reflect, sphereList);
	new Sphere(new Point(-0.5, 0.0, -10.0), 1.0, reflect, sphereList);

	new Sphere(new Point(2.5, -0.75, -6.0), 1.0, diffuseColor, sphereList);
	new Sphere(new Point(-7, 7, -16.0), 1.0, specularColor, sphereList);

	new Sphere(new Point(1.0, 0.5, -3.5), 0.8, refract, sphereList);
	new Sphere(new Point(-0.59, -0.45, -3.0), 0.85, refract, sphereList);

	new Light(new Point(-2.0, 1.0, -10.0), new Color(1.0, 1.0, 1.0), lightList);
	new Light(new Point(-2.0, 20.0, -10.0), new Color(0.4, 0.4, 0.4), lightList);
	
	new Light(new Point(-10.0, 1.5, 5.0), new Color(1.0, 0.0, 0.0), lightList);
	new Light(new Point(10.0, 0.0, 5.0), new Color(0.0, 0.0, 1.0), lightList);

	render();
	glFlush();
}

void init(char* scene) { // glut initialization
	glClearColor(GLfloat(ambient->r), GLfloat(ambient->g), GLfloat(ambient->b), 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, windowSize, 0.0, windowSize);

	if (scene != NULL) {
		switch (atoi(scene)) {
			case(1): scene1();
				break;
			case(2): scene2();
				break;
			case(3): scene3();
				break;
			case(4): scene4();
				break;
			default: cout << "No valid scene argument passed. Loading scene1 by default" << endl;
				scene1();
				break;
		}
	}
	else {
		cout << "No valid scene argument passed. Loading scene1 by default" << endl;
		scene1();
	}
}

int main(int argc, char** argv){ // main

	// GLUT setup
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(500, 100);
	glutInitWindowSize(windowSize, windowSize);
	glutCreateWindow("Ray Tracer");
	
	if (argc > 1)
		init(argv[1]);
	else
		init(NULL);

	glutMainLoop();

	return 0;
}
