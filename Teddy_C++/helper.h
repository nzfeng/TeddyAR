#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <unistd.h>

#include "structs.h"
#include "quaternion.h"
#include "halfedge.h"
#include "poly2tri.h"
#include "main.h"

#include "../Eigen/Dense"

enum TriangleType { TERMINAL, SLEEVE, JUNCTION };

bool isSamePoint(p2t::Point *p1, p2t::Point *p2);

int findIndexOfVertex(const std::vector<p2t::Point*> &polyline, p2t::Point *p);

void constructHalfedge(const std::vector<p2t::Point*> &polyline, const std::vector<p2t::Triangle*> &triangles, 
	std::vector<HEV*> *hevs, std::vector<HEF*> *hefs);

void setTriangleType(HEF *hef);

void labelTriangles(std::vector<HEF*> *hefs);

void prune_branches(std::vector<HEV*> *hevs, std::vector<HEF*> *hef, std::vector<HEF*> *new_hefs);

Eigen::Vector4f getBBoxDimensions(const std::vector<Vertex> &vertices);

bool isOutsideCircle(const float x, const float y, const float r, const float x_c, const float y_c);

Eigen::Vector3f centroidOfTriangle(HEF *hef);


void halfedgeToBuffers(std::vector<HEV*> *hevs, std::vector<HEF*> *hefs, 
	std::vector<Vertex> &vertex_buffer, std::vector<Vec3f> &normal_buffer);

void inflate(const std::vector<p2t::Point*> &polyline, const std::vector<p2t::Triangle*> &triangles,
	std::vector<Vertex> &vertex_buffer, std::vector<Vec3f> &normal_buffer);