#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <map>

#include "structs.h"
#include "quaternion.h"
#include "halfedge.h"
#include "poly2tri.h"
#include "main.h"
#include "helper.h"

#include "../Eigen/Dense"

using namespace std;

const float tol = 1e-5;



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

int findIndexOfHEV(vector<HEV*> *hevs, float x_t, float y_t) {

	for (int i=0; i < hevs->size(); i++) {
		float x = hevs->at(i)->x;
		float y = hevs->at(i)->y;
		if (abs(x - x_t) < tol && abs(y - y_t) < tol) {
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

void moduloSort(vector<int> &fanpoints) {

	sort(fanpoints.begin(), fanpoints.end()); 
	for (int k=0; k < fanpoints.size()-1; k++) {
		int diff = fanpoints[k+1] - fanpoints[k];
		if (diff > 1) {
			rotate(fanpoints.begin(), fanpoints.begin() + k + 1, fanpoints.end());
		}
	}
}

/* Visit each terminal triangle, and perform fanning operation if necessary. 
 * 
 * Also performs re-triangulation; new_hefs should contain all necessary triangles 
 * for sewing up after calling this function.
 */
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
					int exists = findIndexOfHEV(hevs, fan_center(0), fan_center(1));
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
				moduloSort(fanpoints);
				for (int k=1; k < fanpoints.size()-1; k++) // exclude the vertices that make up the circle diameter
				{
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
			moduloSort(fanpoints);
			for (int k=0; k < fanpoints.size()-1; k++) {
				add_face(center_idx, fanpoints[k], fanpoints[(k+1)], hevs, new_hefs, edge_hash);
			}
		}
	}
	// Now add all the other untouched triangles.
	for (int i=0; i < hefs->size(); i++) {
		HEF *hef = hefs->at(i);
		if (hef->visited == false) {
			// unvisited triangles should always be SLEEVE.
			// assert(hef->type == SLEEVE);
			// determine the boundary edge
			HE *he = hef->he;
			do {
				if (he->flip == NULL) {
					break;
				}
				he = he->next;
			} while (he != hef->he);

			float x_t = (he->next->vertex->x + he->next->next->vertex->x)/2.0;
			float y_t = (he->next->vertex->y + he->next->next->vertex->y)/2.0;

			int idx1 = findIndexOfHEV(hevs, x_t, y_t);
			if (idx1 == -1) {
				add_vertex(x_t, y_t, 0.0, false, hevs);
				idx1 = hevs->size()-1;
			}

			x_t = (he->vertex->x + he->next->next->vertex->x)/2.0;
			y_t = (he->vertex->y + he->next->next->vertex->y)/2.0;
			int idx2 = findIndexOfHEV(hevs, x_t, y_t);
			if (idx2 == -1) {
				add_vertex(x_t, y_t, 0.0, false, hevs);
				idx2 = hevs->size()-1;
			}

			add_face(he->vertex->index, he->next->vertex->index, idx1, 
				hevs, new_hefs, edge_hash);
			add_face(he->vertex->index, idx1, idx2,
				hevs, new_hefs, edge_hash);
			add_face(idx1, he->next->next->vertex->index, idx2, 
				hevs, new_hefs, edge_hash);
		}
		// add_face(hef->he->vertex->index, hef->he->next->vertex->index, hef->he->next->next->vertex->index, 
		// 	hevs, new_hefs, edge_hash);

		else if (hef->type == JUNCTION) {
			// check if there's an edge that is either shared w/ an unvisited triangle or a JUNCTION, or its midpoint already exists as an HEV
			HE *he = hef->he;
			do {
				if (he->flip != NULL && (he->flip->face->visited == false || he->flip->face->type == JUNCTION || findIndexOfHEV(hevs, (he->vertex->x+he->next->vertex->x)/2.0, (he->vertex->y+he->next->vertex->y)/2.0) != -1)) {
					Eigen::Vector3f centroid = centroidOfTriangle(hef);
					int centroid_idx = findIndexOfHEV(hevs, centroid(0), centroid(1));
					assert(centroid_idx != -1);
					float x = (he->vertex->x + he->next->vertex->x)/2.0;
					float y = (he->vertex->y + he->next->vertex->y)/2.0;
					int midpoint_idx = findIndexOfHEV(hevs, x, y);
					if (midpoint_idx == -1) {
						add_vertex(x, y, 0.0, false, hevs);
						midpoint_idx = hevs->size()-1;
					}
					add_face(he->vertex->index, midpoint_idx, centroid_idx, 
						hevs, new_hefs, edge_hash);
					add_face(midpoint_idx, he->next->vertex->index, centroid_idx, 
						hevs, new_hefs, edge_hash);
				}
				he = he->next;
			} while (he != hef->he);
		}
	}

	delete_HEF(hefs);
}

void addTriangleToBuffer(vector<Vertex> &vertex_buffer, vector<int> &triangle_buffer, 
	int idx1, int idx2, int idx3) {

	// Make sure the normal is facing the camera.
	Eigen::Vector3f v1; v1 << vertex_buffer[idx2].x - vertex_buffer[idx1].x, vertex_buffer[idx2].y - vertex_buffer[idx1].y, 0.0;
    Eigen::Vector3f v2; v2 << vertex_buffer[idx3].x - vertex_buffer[idx1].x, vertex_buffer[idx3].y - vertex_buffer[idx1].y, 0.0;
    Eigen::Vector3f cross = v1.cross(v2);
    if (cross(2) < 0.0) {
        int tmp = idx1;
        idx1 = idx3;
        idx3 = tmp;
    }

    triangle_buffer.push_back(idx1);
    triangle_buffer.push_back(idx2);
    triangle_buffer.push_back(idx3);
}

/* Elevate spine, and sew up triangles. 
 * 
 * Results in a vector of vertices, and a vector of triangle vertex indices (which
 * will be 3x the length of the vertex vector.)
 */
void create_mesh(vector<HEV*> *hevs, vector<HEF*> *hefs, 
	vector<Vertex> &vertex_buffer, vector<int> &triangle_buffer) {

	vertex_buffer.clear();
	triangle_buffer.clear();

	// number of points to add to each edge; ad hoc length for now
	int numpts;
	// Elevate the spine.
	for (int i=0; i < hevs->size(); i++) {
		// For all vertices on the spine...
		if (hevs->at(i)->isBoundary == false) {
			float length = 0.0;
			int num = 0;
			HE *he = hevs->at(i)->out;
			do {
				if (he->next->vertex->isBoundary == true) {
					float x = (he->vertex->x - he->next->vertex->x);
					float y = (he->vertex->y - he->next->vertex->y);
					length += sqrt(x*x + y*y);
					num++;
				}

        		he = he->flip->next;
			} while (he != hevs->at(i)->out && he->flip != NULL);

			hevs->at(i)->z = length / num;

			numpts = int(length / num / 0.7);
		}
	}

	// Sew up the triangles!

	// Re-index vertices.
	for (int i=0; i < hevs->size(); i++) {
		Vertex v; v.x = hevs->at(i)->x; v.y = hevs->at(i)->y; v.z = hevs->at(i)->z;
		vertex_buffer.push_back(v);
		hevs->at(i)->index = i;
	}

	// Create new points.
	map<pair<int, int>, vector<int>> v_per_he;

	// iterate over all edges
	for (int i=0; i < hefs->size(); i++) {

		HEF *hef = hefs->at(i);
		HE *he = hef->he;
		do {
			if ((he->vertex->isBoundary && !he->next->vertex->isBoundary)
				|| (!he->vertex->isBoundary && he->next->vertex->isBoundary)) {
				
				Eigen::Vector3f e;
				HEV *src = he->vertex;
				HEV *dst = he->next->vertex;

				// if (v_per_he.count(get_edge_key(src->index, dst->index)) > 0) {
				// 	break;
				// }

				// edge is from chordal axis to boundary
				if (src->isBoundary) {
					src = he->next->vertex;
					dst = he->vertex;
				}
				vector<int> new_indices;

				e << dst->x - src->x, dst->y - src->y, 0.0;

				float b = e.norm();
				float a = src->z;
				float scale = 1.0/(numpts+1);

				for (int k=0; k < numpts; k++) {
					Vertex v;
					v.x = src->x + (k+1)*scale*e(0);
					assert((k+1)*scale < 1.0);
					v.y = src->y + (k+1)*scale*e(1);
					float x = (k+1)*scale*b;
					v.z = sqrt((a*a*b*b - a*a*x*x)/(b*b));
					vertex_buffer.push_back(v);
					new_indices.push_back(vertex_buffer.size()-1);
				}

				v_per_he[get_edge_key(src->index, dst->index)] = new_indices;
			}
			he = he->next;
		} while (he != hef->he);
	}


	// Create face array
	for (int i=0; i < hefs->size(); i++) {

		HEF *hef = hefs->at(i);
		// Find the edge that is either on the boundary or on the chordal axis 
		// (between either two interior points or exterior points)
		bool boundary = false;
		HE *he = hef->he;
		do {
			if (he->vertex->isBoundary && he->next->vertex->isBoundary) {
				boundary = true;
				break;
			}
			if (!he->vertex->isBoundary && !he->next->vertex->isBoundary) {
				boundary = false;
				break;
			}
			he = he->next;
		} while (he != hef->he);

		vector<int> v_per_e1 = v_per_he[get_edge_key(he->next->vertex->index, he->next->next->vertex->index)];
		vector<int> v_per_e2 = v_per_he[get_edge_key(he->next->next->vertex->index, he->next->next->next->vertex->index)];

		if (boundary) {
			reverse(v_per_e1.begin(), v_per_e1.end());
			reverse(v_per_e2.begin(), v_per_e2.end());
		}

		assert(v_per_e1.size() == numpts);
		assert(v_per_e2.size() == numpts);

		/* Add triangles! */
		
		// first triangles
		addTriangleToBuffer(vertex_buffer, triangle_buffer, 
			he->vertex->index, he->next->vertex->index, v_per_e2[0]);

		addTriangleToBuffer(vertex_buffer, triangle_buffer, 
			v_per_e1[0], v_per_e2[0], he->next->vertex->index);

		for (int k=0; k < numpts-1; k++) {

			addTriangleToBuffer(vertex_buffer, triangle_buffer, 
				v_per_e1[k], v_per_e2[k], v_per_e1[k+1]);

			addTriangleToBuffer(vertex_buffer, triangle_buffer, 
				v_per_e2[k], v_per_e2[k+1], v_per_e1[k+1]);
		}
		// end triangle
		addTriangleToBuffer(vertex_buffer, triangle_buffer, 
			v_per_e1[numpts-1], he->next->next->vertex->index, v_per_e2[numpts-1]);
	}
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

	vertex_buffer.clear();
	normal_buffer.clear();

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
	vector<Vertex> &vertex_buffer, vector<Vec3f> &normal_buffer, vector<int> &triangle_buffer) {


	vector<HEV*> *hevs = new vector<HEV*>();
    vector<HEF*> *hefs = new vector<HEF*>();
	constructHalfedge(polyline, triangles, hevs, hefs);
	cerr << "Halfedges constructed..." << endl;

	labelTriangles(hefs);
	cerr << "Triangles labeled..." << endl;

	vector<HEF*> *new_hefs = new vector<HEF*>();
	prune_branches(hevs, hefs, new_hefs);
	cerr << "Branches pruned..." << endl;

	create_mesh(hevs, new_hefs, vertex_buffer, triangle_buffer);
	cerr << "Mesh created..." << endl;

	// halfedgeToBuffers(hevs, new_hefs, vertex_buffer, normal_buffer);

	delete_HEV(hevs);
	delete_HEF(new_hefs);
	
}