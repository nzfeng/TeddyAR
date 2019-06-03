#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <unistd.h>

#include "structs.h"
#include "quaternion.h"
#include "halfedge.h"
#include "poly2tri.h"
#include "main.h"
#include "helper.h"

#include "../Eigen/Dense"

using namespace std;

const float tol = 1e-10;



bool isSamePoint(p2t::Point *p1, p2t::Point *p2) {

	return (abs(p1->x - p2->x) < tol && abs(p1->y - p2->y) < tol);
}

int findIndexOfVertex(const vector<p2t::Point*> &polyline, p2t::Point *p) {

	for (int i=0; i < polyline.size(); i++) {

		if (isSamePoint(polyline[i], p)) {
			return i;
		}
	}
	return -1;
}

int findIndexOfHEV(vector<HEV*> *hevs, const Eigen::Vector3f &v) {

	for (int i=0; i < hevs->size(); i++) {
		float x = hevs->at(i)->x;
		float y = hevs->at(i)->y;
		if (abs(x - v(0)) < tol && abs(y - v(1))) {
			return i;
		}
	}
	return -1;
}

void constructHalfedge(const vector<p2t::Point*> &polyline, const vector<p2t::Triangle*> &triangles, 
	vector<HEV*> *hevs, vector<HEF*> *hefs) {

	map<pair<int, int>, HE*> edge_hash;

	for (int i=0; i < polyline.size(); i++) {
		p2t::Point *p = polyline[i]; 
		add_vertex(p->x, p->y, 0.0, true, hevs);
	}

	for (int i=0; i < triangles.size(); i++) {
		int indices[3];
		for (int j=0; j < 3; j++) {
			p2t::Point *p = triangles[i]->GetPoint(j);
			int idx = findIndexOfVertex(polyline, p);
			assert(idx != -1);
			indices[j] = idx;
		}

		add_face(indices[0], indices[1], indices[2], hevs, hefs, edge_hash);
	}

	// for(map<pair<int,int>, HE*>::iterator it = edge_hash.begin(); it != edge_hash.end(); ++it) {
	//   	cerr << it->first.first << " " << it->first.second << "\n";
	// }
}

/* Determine if triangle is TERMINAL, SLEEVE, or JUNCTION. */
void setTriangleType(HEF *hef) {

	int numBoundaryEdges = 0;
	HE *he = hef->he;
	// Go around all the edges of the polygon
	do {
		if (he->flip == NULL) {
			numBoundaryEdges++;
		}
		he = he->next;
	} while (he != hef->he);

	// assert(numBoundaryEdges < 3);

	if (numBoundaryEdges == 0) {
		hef->type = JUNCTION;
	}
	else if (numBoundaryEdges == 1) {
		hef->type = SLEEVE;
	}
	else if (numBoundaryEdges == 2) {
		hef->type = TERMINAL;
	}
}

/* Label each triangle as TERMINAL, SLEEVE, or JUNCTION. */
void labelTriangles(vector<HEF*> *hefs) {

	for (int i=0; i < hefs->size(); i++) {

		setTriangleType(hefs->at(i));
	}
}



/* Visit each terminal triangle, and perform fanning operation if necessary. */
void prune_branches(vector<HEV*> *hevs, vector<HEF*> *hefs, vector<HEF*> *new_hefs) {

	map<pair<int, int>, HE*> edge_hash;

	// Vertices to be used to form new triangles, and faces to be replaced.
	vector<int> fanpoints;

	for (int i=0; i < hefs->size(); i++) {
		HEF *hef = hefs->at(i);

		HE *intHE;

		if (hef->type == TERMINAL) {
			intHE = hef->he;
			Eigen::Vector3f fan_center;
			int center_idx;
			fanpoints.clear();

			while (true) {
				assert(intHE != NULL);
				
				// go around polygon edges
				HE *prev = intHE;
				do {
					intHE = intHE->next;
					if (intHE->flip != NULL) {
						break;
					}
				} while (intHE != prev);

				hef = intHE->face; hef->visited = true;
				// If junction triangle, break
				if (hef->type == JUNCTION) {
					fan_center = centroidOfTriangle(hef);
					// only add to HEVs if it doesn't exist already
					int exists = findIndexOfHEV(hevs, fan_center);
					if (exists == -1) {
						add_vertex(fan_center(0), fan_center(1), 0.0, false, hevs); // add fan center
						center_idx = hevs->size()-1;
					}
					else {
						center_idx = exists;
					}
					
					break;
				}

				// add current triangle's vertices to fanpoints
				
				// the point not part of the interior edge
				int idx1 = intHE->next->next->vertex->index;
				if (count(fanpoints.begin(), fanpoints.end(), idx1) == 0) {
					fanpoints.push_back(idx1);
				}

				// intHE should now be the interior halfedge. Do hemisphere test.
				float x1 = intHE->vertex->x;
				float y1 = intHE->vertex->y;
				float x2 = intHE->next->vertex->x;
				float y2 = intHE->next->vertex->y;
				bool isOutside = false;
				for (int k=0; k < fanpoints.size(); k++) {
					if (isOutsideCircle(hevs->at(fanpoints[k])->x, hevs->at(fanpoints[k])->y, sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/2.0, (x1+x2)/2.0, (y1+y2)/2.0)) {
						fan_center << (x1+x2)/2.0, (y1+y2)/2.0, 0.0;
						add_vertex(fan_center(0), fan_center(1), 0.0, false, hevs); // add fan center 
						center_idx = hevs->size()-1;
						isOutside = true;
						break;
					} 
				}
				// add rest of triangle points (if not junction triangle)
				if (count(fanpoints.begin(), fanpoints.end(), intHE->vertex->index) == 0) {
					fanpoints.push_back(intHE->vertex->index);
				}
				if (count(fanpoints.begin(), fanpoints.end(), intHE->next->vertex->index) == 0) {
					fanpoints.push_back(intHE->next->vertex->index);
				}

				if (isOutside) {
					break;
				}
				
				assert(intHE->flip != NULL);
				intHE = intHE->flip;

			}
			// Perform fanning.	
			// Sort point indices from least to greatest (modulo number of boundary points)
			sort(fanpoints.begin(), fanpoints.end()); 
			for (int k=0; k < fanpoints.size()-1; k++) {
				int diff = fanpoints[k+1] - fanpoints[k];
				if (diff > 1) {
					rotate(fanpoints.begin(), fanpoints.begin() + k + 1, fanpoints.end());
				}
			}
			for (int k=0; k < fanpoints.size()-1; k++) {
				add_face(center_idx, fanpoints[k], fanpoints[(k+1)], hevs, new_hefs, edge_hash);
			}
		}
	}
	// Now add all the other untouched triangles.
	for (int i=0; i < hefs->size(); i++) {
		HEF *hef = hefs->at(i);
		if (hef->visited == false) {
			if (hef->type == SLEEVE) {
				// determine the boundary edge
				HE *he = hef->he;
				do {
					he = he->next;
					if (he->flip == NULL) {
						break;
					}
				} while (he != hef->he);
				add_vertex((he->next->vertex->x + he->next->next->vertex->x)/2.0, (he->next->vertex->y + he->next->next->vertex->y)/2.0, 0.0, false, hevs); 
				add_vertex((he->vertex->x + he->next->next->next->vertex->x)/2.0, (he->vertex->y + he->next->next->next->vertex->y)/2.0, 0.0, false, hevs);
				add_face(he->vertex->index, he->next->vertex->index, hevs->size()-2, 
					hevs, new_hefs, edge_hash);
				add_face(he->vertex->index, hevs->size()-2, hevs->size()-1,
					hevs, new_hefs, edge_hash);
				add_face(hevs->size()-2, he->next->next->vertex->index, hevs->size()-1, 
					hevs, new_hefs, edge_hash);
			}
			// add_face(hef->he->vertex->index, hef->he->next->vertex->index, hef->he->next->next->vertex->index, 
			// 	hevs, new_hefs, edge_hash);
		}
	}

	delete_HEF(hefs);
}

// Determine the final spine (the final 2D triangulation.)
void final_triangulation(vector<HEV*> *hevs, vector<HEF*> *hefs, vector<HEF*> *new_hefs) {

	map<pair<int, int>, HE*> edge_hash;
}


// Return x_min, y_min, x_max, y_max.
Eigen::Vector4f getBBoxDimensions(const vector<Vertex> &vertices) {

	float x_min = numeric_limits<float>::infinity();
	float y_min = numeric_limits<float>::infinity();
	float x_max = 0;
	float y_max = 0;

	for (int i=0; i < vertices.size(); i++) {
		float x = vertices[i].x;
		float y = vertices[i].y;
		if (x < x_min) {
			x_min = x;
		}
		if (x > x_max) {
			x_max = x;
		}
		if (y < y_min) {
			y_min = y;
		}
		if (y > y_max) {
			y_max = y;
		}
	}

	Eigen::Vector4f result;
	result << x_min, y_min, x_max, y_max;
	return result;
}

/* Determine if a point given by (x, y) is outside a circle defined by radius r
 * and center at (x_c, y_c).
 */
bool isOutsideCircle(const float x, const float y, const float r, const float x_c, const float y_c) {

	return ((x-x_c)*(x-x_c) + (y-y_c)*(y-y_c) > r*r);
}


/* Return the centroid of a triangle given by its HEF. */
Eigen::Vector3f centroidOfTriangle(HEF *hef) {
	HEV *a = hef->he->vertex;
	HEV *b = hef->he->next->vertex;
	HEV *c = hef->he->next->next->vertex;

	float x = (1/3.0) * (a->x + b->x + c->x);
	float y = (1/3.0) * (a->y + b->y + c->y);
	
	Eigen::Vector3f v; v << x, y, 0.0;
	return v;
}

void halfedgeToBuffers(vector<HEV*> *hevs, vector<HEF*> *hefs, 
	vector<Vertex> &vertex_buffer, vector<Vec3f> &normal_buffer) {

	// not the actual function, using for testing right now
	for (int i=0; i < hefs->size(); i++) {
		Vertex v1; 
		v1.x = hefs->at(i)->he->vertex->x;
		v1.y = hefs->at(i)->he->vertex->y;
		v1.z = hefs->at(i)->he->vertex->z;

		Vertex v2; 
		v2.x = hefs->at(i)->he->next->vertex->x;
		v2.y = hefs->at(i)->he->next->vertex->y;
		v2.z = hefs->at(i)->he->next->vertex->z;

		Vertex v3; 
		v3.x = hefs->at(i)->he->next->next->vertex->x;
		v3.y = hefs->at(i)->he->next->next->vertex->y;
		v3.z = hefs->at(i)->he->next->next->vertex->z;

		vertex_buffer.push_back(v1);
		vertex_buffer.push_back(v2);
		vertex_buffer.push_back(v3);
	}
}

void inflate(const vector<p2t::Point*> &polyline, const vector<p2t::Triangle*> &triangles,
	vector<Vertex> &vertex_buffer, vector<Vec3f> &normal_buffer) {


	vector<HEV*> *hevs = new vector<HEV*>();
    vector<HEF*> *hefs = new vector<HEF*>();
	constructHalfedge(polyline, triangles, hevs, hefs);
	cerr << "Halfedges constructed..." << endl;

	labelTriangles(hefs);
	cerr << "Triangles labeled..." << endl;

	vector<HEF*> *new_hefs = new vector<HEF*>();
	prune_branches(hevs, hefs, new_hefs);
	cerr << "Branches pruned..." << endl;

	halfedgeToBuffers(hevs, new_hefs, vertex_buffer, normal_buffer);

	delete_HEV(hevs);
	delete_HEF(new_hefs);
	
}