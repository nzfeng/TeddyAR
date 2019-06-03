#pragma once

#include <string>
#include <vector>
#include <math.h>
#include "../Eigen/Dense"

struct Vec3f
{
	float x, y, z;
};

struct Vertex
{
    float x, y, z;
    bool isBoundary;
};

/* A transformation vector struct that represents either a translation, 
 * scaling, or rotation vector.
 */
struct transvector {
	std::string type; // 'r', 's', or 't'
	float x, y, z, theta;
};


struct camera {
	// a struct representing initial camera translation vector
	transvector position;
	// a four-entry array (x, y, z, theta) representing camera rotation
	transvector orientation;
};

struct frustum {
	// perspective parameters
	float near;
	float far;
	float left;
	float right;
	float top;
	float bottom;
};

struct point_light {
	// homogeneous coordinates of point light in world space
	float position[4];
	// r,g,b values of point light
	float rgb[3];
	// attentuation parameter of light
	float k;
};