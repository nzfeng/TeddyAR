#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include <iterator>

#include "structs.h"
#include "quaternion.h"
#include "poly2tri.h"
#include "halfedge.h"
#include "main.h"

#include "../Eigen/Dense"

using namespace std;




// /* Given a terminal triangle, determine the index of the point that does not 
//  * lie on the interior edge. 
//  */
// int getExteriorPoint(Triangle &t, const vector<Point> &points) {

// }


// /* Constrained Delaunay triangulation. */

// /* Obtain the chordal axis by connecting the midpoints of the internal edges. */

// /* Prune insignificant branches. */
// vector<Point> removedPoints;
// for (int i=0; i < triangles.size(); i++) {
// 	// for all terminal triangles...
// 	if (getTriangleType(triangles[i]) == TERMINAL) {
// 		Triangle t = triangles[i];
// 		while (true) {
// 			int extPoint = getExteriorPoint(t);
// 			removedPoints.push_back(getExteriorPoint(t), points);

// 		}
// 		// There should only be one adjacent triangle, since terminal triangles have 1 interior edge.

// 	}
// }

/* Re-triangulate the mesh. Remove short edges and small triangles. */

/* Elevate each vertex on the spine. */

/* Elevate each non-spinal edges. */

/* Sew together the neighboring elevated edges. */
