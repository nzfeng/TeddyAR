using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Poly2Tri;

public class Teddy
{

    // The final method combining everything. 
    public static TeddyResult inflate(List<Vector3> points)
    {
        Debug.Log("Inside Teddy: " + points.Count);
        TeddyResult res = new TeddyResult(new List<Vector3>(), new List<Vector3>(), new List<int>());

        if (isSelfIntersecting(points) || points.Count < 3)
        {
            Debug.Log("Invalid input");
            return res;
        }

        points = resample(points);
        Debug.Log("Points resampled");
        Debug.Log("Teddy # points: " + points.Count);

        if (!isPolygonClockwise(points))
        {
            Debug.Log("Points reversed");
            points.Reverse();
        }

        List<Vertex> hevs = new List<Vertex>();
        List<Triangle> hefs = p2t_CDT(points, hevs);

        //List<Triangle> hefs = ConstrainedDelaunay.GenerateTriangulation(points);

        for (int i=0; i < hefs.Count; i++)
        {
            hefs[i].visited = false;
        }
        Debug.Log("CDT generated...");

        List<HalfEdge> hes = TransformRepresentation.TransformFromTriangleToHalfEdge(hefs);
        Debug.Log("Halfedges generated...");

        //// Get vertices from halfedges
        //List<Vertex> hevs = new List<Vertex>(hes.Count);
        //for (int i=0; i < hes.Count; i++)
        //{
        //    hevs[hes[i].v.index] = hes[i].v;
        //}
        //Debug.Log("Halfedges labeled...");

        // Label triangles.
        labelTriangles(hefs);
        Debug.Log("Triangles labeled...");

        // Prune branches.
        List<Triangle> new_hefs = prune_branches(hevs, hefs);
        hes = TransformRepresentation.TransformFromTriangleToHalfEdge(new_hefs);
        Debug.Log("Branches pruned...");

        Debug.Log("points " + points.Count);
        Debug.Log("hevs " + hevs.Count);
        Debug.Log("new_hefs " + hevs.Count);

        //List<Vector3> pts = new List<Vector3>();
        //for (int i=0; i < new_hefs.Count; i++)
        //{
        //    pts.Add(new Vector3(new_hefs[i].v1.position.x, new_hefs[i].v1.position.y, new_hefs[i].v1.position.z));
        //    pts.Add(new Vector3(new_hefs[i].v2.position.x, new_hefs[i].v2.position.y, new_hefs[i].v2.position.z));
        //    pts.Add(new Vector3(new_hefs[i].v3.position.x, new_hefs[i].v3.position.y, new_hefs[i].v3.position.z));
        //    pts.Add(new Vector3(new_hefs[i].v1.position.x, new_hefs[i].v1.position.y, new_hefs[i].v1.position.z));
        //}
        //res.points = pts;
        //return res;

        // Construct mesh.
        res = create_mesh(hevs, new_hefs);
        Debug.Log("Mesh generated...");
        return res;
    }


    public static List<Triangle> p2t_CDT(List<Vector3> points, List<Vertex> hevs)
    {
        List<PolygonPoint> pts = new List<PolygonPoint>();
        for (int i=0; i < points.Count; i++)
        {
            PolygonPoint p = new PolygonPoint(points[i].x, points[i].y);
            pts.Add(p);
        }
        Polygon poly = new Polygon(pts);
        P2T.Triangulate(poly);
        Debug.Log("p2T: " + poly.Triangles.Count);

        // Generate Triangles from p2t triangles. First create Vertex list.
        for (int i=0; i < points.Count; i++)
        {
            hevs.Add(new Vertex(points[i]));
            hevs[i].index = i;
            hevs[i].isBoundary = true;
        }
        List<Triangle> triangles = new List<Triangle>();
        for (int i=0; i < poly.Triangles.Count; i++)
        {
            int idx1 = FindIndexOfVertex(hevs, poly.Triangles[i].Points[0]);
            int idx2 = FindIndexOfVertex(hevs, poly.Triangles[i].Points[1]);
            int idx3 = FindIndexOfVertex(hevs, poly.Triangles[i].Points[2]);
            triangles.Add(new Triangle(hevs[idx1], hevs[idx2], hevs[idx3]));
        }

        return triangles;
    }

    public static bool IsSamePoint(TriangulationPoint p1, TriangulationPoint p2)
    {
        double tol = 1e-10;
        return (Math.Abs(p1.X - p2.X) < tol && Math.Abs(p1.Y - p2.Y) < tol);
    }

    public static int FindIndexOfVertex(List<Vertex> vertices, TriangulationPoint p)
    {
        double tol = 1e-10;
        for (int i=0; i < vertices.Count; i++)
        {
            if (Math.Abs(vertices[i].position.x - p.X) < tol && Math.Abs(vertices[i].position.y - p.Y) < tol)
            {
                return i;
            }
        }
        return -1;
    }

    /* Check if the input polyline intersects itself. This is the most naive implementation, O(n^2). */
    public static bool isSelfIntersecting(List<Vector3> points)
    {
        int n = points.Count;

        for (int i=0; i < n; i++)
        {

            for (int j=i+1; j < n; j++)
            {

                if (Intersections.AreLinesIntersecting(points[i].XY(), points[(i+1)%n].XY(), points[j].XY(), points[(j+1)%n].XY(), true))
                {
                    return true;
                }
            }
        }
        return false;
    }


    /* Re-sample points to remove redundant points and short edges. */
    public static List<Vector3> resample(List<Vector3> points)
    {
        List<Vector3> newPoints = new List<Vector3>();
        int n = points.Count;
        float length = 0.05f; // change this if necessary
        newPoints.Add(new Vector3(points[0].x, points[0].y, 0.0f));
        int i = 0;

        while (i < n-1)
        {
            Vector3 p2 = points[(i + 1) % n];
            float x = newPoints[newPoints.Count - 1].x;
            float y = newPoints[newPoints.Count - 1].y;
            float norm = (float) Math.Sqrt((x - p2.x) * (x - p2.x) + (y - p2.y) * (y - p2.y));

            if (norm > length)
            {
                float a = length / norm;
                newPoints.Add(new Vector3(x + a * (p2.x - x), y + a * (p2.y - y), 0.0f));
            }
            else
            {
                i += 1;
            }
        }

        return newPoints;
    }

    /* Return true if the input points are already ordered clockwise. 
     * (The CDT assumes constraints are ordered clockwise.) */
    public static bool isPolygonClockwise(List<Vector3> points)
    {
        float tol = 1e-10f;
        // Find lowermost vertex
        int k = 0;
        int n = points.Count;
        for (int i=1; i < n; i++)
        {
            if ((points[i].y < points[k].y) || ((float) Math.Abs(points[i].y - points[k].y) < tol && points[i].x > points[k].x))
            {
                k = i;
            }
        }
        // deal with negative mod
        int before = (k - 1)%n;
        if (before < 0)
        {
            before += n;
        }
        // return Geometry.IsTriangleOrientedClockwise(points[before].XY(), points[k].XY(), points[(k + 1) % n].XY());

        Vector3 e2 = points[(k + 1) % n] - points[k];
        Vector3 e1 = points[before] - points[k];
        Vector3 r = Vector3.Cross(e1, e2);
        if (r.z < 0.0f)
        {
            return false;
        }
        return true;
    }


    /* Visit each terminal triangle, and perform fanning operation if necessary.
     * 
     * Also performs re-triangulation.
     */
    public static List<Triangle> prune_branches(List<Vertex> hevs, List<Triangle> hefs)
    {
        List<Triangle> new_hefs = new List<Triangle>();

        List<int> fanpoints = new List<int>();

        for (int i=0; i < hefs.Count; i++)
        {
            Triangle hef = hefs[i];
            HalfEdge intHE;

            if (hef.type == 0)
            {
                intHE = hef.halfEdge;
                Vector3 fan_center = new Vector3(0.0f, 0.0f, 0.0f);
                int center_idx = 0;
                fanpoints.Clear();

                while (true)
                {
                    // go around polygon edges
                    HalfEdge prev = intHE;
                    do
                    {
                        intHE = intHE.nextEdge;
                        if (intHE.oppositeEdge != null)
                        {
                            break;
                        }
                    } while (intHE != prev);

                    hef = intHE.t; hef.visited = true;
                    // If junction triangle, break.
                    if (hef.type == 2)
                    {
                        fan_center = centroidOfTriangle(hef);
                        // only add to HEVs if it doesn't exist already
                        int exists = findIndexOfHEV(hevs, fan_center[0], fan_center[1]);
                        if (exists == -1)
                        {
                            add_vertex(new Vertex(new Vector3(fan_center[0], fan_center[1], 0.0f)), false, hevs);
                            center_idx = hevs.Count - 1;
                        }
                        else
                        {
                            center_idx = exists;
                        }
                        break;
                    }
                    // the point not part of the interior edge
                    int idx1 = intHE.nextEdge.nextEdge.src.index;
                    if (fanpoints.IndexOf(idx1) == -1)
                    {
                        fanpoints.Add(idx1);
                    }
                    // intHE should now be the interior halfedge. Do hemisphere test.
                    float x1 = intHE.src.position.x;
                    float y1 = intHE.src.position.y;
                    float x2 = intHE.nextEdge.src.position.x;
                    float y2 = intHE.nextEdge.src.position.y;
                    bool isOutside = false;
                    fanpoints = moduloSort(fanpoints);
                    for (int k=1; k < fanpoints.Count-1; k++) // exclude the vertices that make up the circle diameter
                    {
                        if (isOutsideCircle(hevs[fanpoints[k]].position.x, hevs[fanpoints[k]].position.y,
                            (float)Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2.0f, (x1 + x2) / 2.0f, (y1+y2)/2.0f))
                        {
                            fan_center.x = (x1 + x2) / 2.0f;
                            fan_center.y = (y1 + y2) / 2.0f;
                            add_vertex(new Vertex(fan_center), false, hevs); // add fan center
                            center_idx = hevs.Count - 1;
                            isOutside = true;
                            break;
                        }
                    }
                    // add rest of triangle points (if not junction triangle)
                    if (fanpoints.IndexOf(intHE.src.index) == -1)
                    {
                        fanpoints.Add(intHE.src.index);
                    }
                    if (fanpoints.IndexOf(intHE.nextEdge.src.index) == -1)
                    {
                        fanpoints.Add(intHE.nextEdge.src.index);
                    }
                    if (isOutside)
                    {
                        break;
                    }

                    intHE = intHE.oppositeEdge;
                }
                // Perform fanning.
                fanpoints = moduloSort(fanpoints);
                for (int k=0; k < fanpoints.Count-1; k++)
                {
                    add_face(center_idx, fanpoints[k], fanpoints[k + 1], hevs, new_hefs);
                }
            }
        }
        // Now add all the other untouched triangles.
        for (int i=0; i < hefs.Count; i++)
        {
            Triangle hef = hefs[i];
            if (hef.visited == false)
            {
                HalfEdge he = hef.halfEdge;
                do
                {
                    if (he.oppositeEdge == null)
                    {
                        break;
                    }
                    he = he.nextEdge;
                } while (he != hef.halfEdge);

                float x_t = (he.nextEdge.src.position.x + he.nextEdge.nextEdge.src.position.x) / 2.0f;
                float y_t = (he.nextEdge.src.position.y + he.nextEdge.nextEdge.src.position.y) / 2.0f;

                int idx1 = findIndexOfHEV(hevs, x_t, y_t);
                if (idx1 == -1)
                {
                    add_vertex(new Vertex(new Vector3(x_t, y_t, 0.0f)), false, hevs);
                    idx1 = hevs.Count - 1;
                }

                x_t = (he.src.position.x + he.nextEdge.nextEdge.src.position.x) / 2.0f;
                y_t = (he.src.position.y + he.nextEdge.nextEdge.src.position.y) / 2.0f;
                int idx2 = findIndexOfHEV(hevs, x_t, y_t);
                if (idx2 == -1)
                {
                    add_vertex(new Vertex(new Vector3(x_t, y_t, 0.0f)), false, hevs);
                    idx2 = hevs.Count - 1;
                }

                add_face(he.src.index, he.nextEdge.src.index, idx1, hevs, new_hefs);
                add_face(he.src.index, idx1, idx2, hevs, new_hefs);
                add_face(idx1, he.nextEdge.nextEdge.src.index, idx2, hevs, new_hefs);
            }
            // if JUNCTION
            else if (hef.type == 2)
            {
                // check if there's an edge that is either shared w/ an unvisited triangle or a JUNCTION, or its midpoint already exists as an HEV
                HalfEdge he = hef.halfEdge;
                do
                {
                    if (he.oppositeEdge != null && (he.oppositeEdge.t.visited == false || he.oppositeEdge.t.type == 2 || findIndexOfHEV(hevs, (he.src.position.x + he.nextEdge.src.position.x) / 2.0f, (he.src.position.y + he.nextEdge.src.position.y) / 2.0f) != -1))
                    {
                        Vector3 centroid = centroidOfTriangle(hef);
                        int centroid_idx = findIndexOfHEV(hevs, centroid.x, centroid.y);
                        float x = (he.src.position.x + he.nextEdge.src.position.x) / 2.0f;
                        float y = (he.src.position.y + he.nextEdge.src.position.y) / 2.0f;
                        int midpoint_idx = findIndexOfHEV(hevs, x, y);
                        if (midpoint_idx == -1)
                        {
                            add_vertex(new Vertex(new Vector3(x, y, 0.0f)), false, hevs);
                            midpoint_idx = hevs.Count - 1;
                        }
                        add_face(he.src.index, midpoint_idx, centroid_idx, hevs, new_hefs);
                        add_face(midpoint_idx, he.nextEdge.src.index, centroid_idx, hevs, new_hefs);
                    }
                    he = he.nextEdge;
                } while (he != hef.halfEdge);
            }
        }
        return new_hefs;
    }

    /* Elevate spine, and sew up triangles. 
     * 
     * Results in a vector of vertices, a vector of vertex normals, and a vector 
     * of triangle vertex indices (which will be 3x the length of the vertex vector.)
     */    
    public static TeddyResult create_mesh(List<Vertex> hevs, List<Triangle> hefs)
    {
        List<Vertex> vertices_lst = new List<Vertex>();
        List<Triangle> triangles_lst = new List<Triangle>();

        Debug.Log("CM hevs " + hevs.Count);
        Debug.Log("CM new_hefs " + hevs.Count);


        // number of points to add to each edge; this number is ad hoc.
        // Initialized to random number (5) due to C# idiosyncrasy
        int numpts = 5;
        // Elevate the spine.
        for (int i=0; i < hevs.Count; i++)
        {
            // For all vertices on the spine...
            if (hevs[i].isBoundary == false)
            {
                int num = 0;
                float length = 0.0f;
                HalfEdge he = hevs[i].halfEdge;
                do
                {
                    if (he.nextEdge.src.isBoundary == true)
                    {
                        float x = (he.src.position.x - he.nextEdge.src.position.x);
                        float y = (he.src.position.y - he.nextEdge.src.position.y);
                        length += (float)Math.Sqrt(x * x + y * y);
                        num++;
                    }
                    if (he.oppositeEdge == null)
                    {
                        break;
                    }
                    he = he.oppositeEdge.nextEdge;

                } while (he != hevs[i].halfEdge);
                // Debug.Log("Elevate: " + i);
                hevs[i].position = new Vector3(hevs[i].position.x, hevs[i].position.y, length / num);

                numpts = (int) (length / num / 0.01f); // this float is ad hoc
            }
        }
        Debug.Assert(numpts > 0);


        // Sew up the triangles!

        // Re-index the vertices.
        for (int i=0; i < hevs.Count; i++)
        {
            hevs[i].index = i;
            Vertex v = new Vertex(new Vector3(hevs[i].position.x, hevs[i].position.y, hevs[i].position.z));
            add_vertex(v, true, vertices_lst);
        }

        // Create new points.
        Dictionary<Pair, List<int>> v_per_he = new Dictionary<Pair, List<int>>(); // Pair is <int, int>

        // Iterate over all the edges.
        for (int i=0; i < hefs.Count; i++)
        {
            Triangle hef = hefs[i];
            HalfEdge he = hef.halfEdge;
            do
            {
                if ((he.src.isBoundary && !he.nextEdge.src.isBoundary)
                    || (!he.src.isBoundary && he.nextEdge.src.isBoundary))
                {
                    Vertex src = he.src;
                    Vertex dst = he.nextEdge.src;

                    if (src.isBoundary)
                    {
                        src = he.nextEdge.src;
                        dst = he.src;
                    }
                    List<int> new_indices = new List<int>();

                    Vector3 e = new Vector3(dst.position.x - src.position.x, dst.position.y - src.position.y, 0.0f);

                    float b = e.magnitude;
                    float a = src.position.z;
                    float scale = 1.0f / (numpts + 1);

                    for (int k = 0; k < numpts; k++)
                    {
                        float x = (k + 1) * scale * b;
                        Vertex v = new Vertex(new Vector3(src.position.x + (k + 1) * scale * e.x, src.position.y + (k + 1) * scale * e.y, (float)Math.Sqrt((a * a * b * b - a * a * x * x) / (b * b))));
                        add_vertex(v, false, vertices_lst);
                        new_indices.Add(vertices_lst.Count - 1);
                    }
                    if (!IsKeyInDict(v_per_he, get_edge_key(src.index, dst.index))) 
                    {
                        v_per_he.Add(get_edge_key(src.index, dst.index), new_indices);
                    }

                }
                he = he.nextEdge;
            } while (he != hef.halfEdge);
        }


        // Create face array.
        for (int i=0; i < hefs.Count; i++)
        {
            // Debug.Log("Create face array: " + i);
            Triangle hef = hefs[i];
            // Find the edge that is either on the boundary or on the chordal axis 
            // (between either two interior points or exterior points)
            bool boundary = false;
            HalfEdge he = hef.halfEdge;
            do
            {
                Debug.Assert(he.nextEdge != null); 
                if (he.src.isBoundary && he.nextEdge.src.isBoundary)
                {
                    boundary = true;
                    break;
                }
                if (!he.src.isBoundary && !he.nextEdge.src.isBoundary)
                {
                    boundary = false;
                    break;
                }

                he = he.nextEdge;
            } while (he != hef.halfEdge);

            List<int> v_per_e1 = GetListFromDict(v_per_he, get_edge_key(he.nextEdge.src.index, he.nextEdge.nextEdge.src.index));
            List<int> v_per_e2 = GetListFromDict(v_per_he, get_edge_key(he.nextEdge.nextEdge.src.index, he.src.index));

            Debug.Assert(v_per_e1.Count == numpts);
            Debug.Assert(v_per_e2.Count == numpts);

            if (boundary)
            {
                Debug.Log("Reversed");
                v_per_e1.Reverse();
                v_per_e2.Reverse();
            }

            // Add triangles! 

            // first triangles
            add_face(he.src.index, he.nextEdge.src.index, v_per_e2[0], vertices_lst, triangles_lst);

            add_face(v_per_e1[0], v_per_e2[0], he.nextEdge.src.index, vertices_lst, triangles_lst);

            for (int k=0; k < numpts-1; k++)
            {
                add_face(v_per_e2[k], v_per_e1[k], v_per_e2[k + 1], vertices_lst, triangles_lst);

                add_face(v_per_e1[k], v_per_e1[k + 1], v_per_e2[k + 1], vertices_lst, triangles_lst);

            }
            // end triangle
            add_face(v_per_e1[numpts - 1], he.nextEdge.nextEdge.src.index, v_per_e2[numpts - 1], vertices_lst, triangles_lst);
        }

        // Build arrays.

        List<HalfEdge> hes = TransformRepresentation.TransformFromTriangleToHalfEdge(triangles_lst); // construct halfedge representation

        List<Vector3> vertices_vec = new List<Vector3>();
        List<Vector3> normals_vec = new List<Vector3>();
        List<int> triangles_vec = new List<int>();

        for (int i=0; i < vertices_lst.Count; i++)
        {
            vertices_vec.Add(new Vector3(vertices_lst[i].position.x*Screen.width/Screen.height, vertices_lst[i].position.y, vertices_lst[i].position.z));
            // normals_vec.Add(compute_vnormal(vertices_lst[i]));
        }
        for (int i=0; i < triangles_lst.Count; i++)
        {
            addTriangleToBuffer(vertices_lst, triangles_lst[i], triangles_vec);
        }

        // Add reflection
        int n = vertices_lst.Count;
        for (int i = 0; i < n; i++)
        {
            Vertex v = new Vertex(new Vector3(vertices_lst[i].position.x * Screen.width / Screen.height, vertices_lst[i].position.y, -vertices_lst[i].position.z));
            add_vertex(v, false, vertices_lst);
            vertices_vec.Add(v.position);
        }
        for (int i = 0; i < triangles_lst.Count; i++)
        {
            int idx1 = triangles_lst[i].v1.index + n;
            int idx2 = triangles_lst[i].v2.index + n;
            int idx3 = triangles_lst[i].v3.index + n;

            if (!IsTriangleClockwise(idx1, idx2, idx3, vertices_lst))
            {
                int tmp = idx1;
                idx1 = idx3;
                idx3 = tmp;
            }

            triangles_vec.Add(idx1);
            triangles_vec.Add(idx2);
            triangles_vec.Add(idx3);
        }

        TeddyResult res = new TeddyResult(vertices_vec, vertices_vec, triangles_vec);
        return res;
    }




    public static List<int> GetListFromDict(Dictionary<Pair, List<int>> v_per_he, Pair pair)
    {
        List<int> newList = new List<int>(); // deep copy

        foreach (KeyValuePair<Pair, List<int>> item in v_per_he)
        {
            if (item.Key.first == pair.first && item.Key.second == pair.second)
            {
                for (int i=0; i < item.Value.Count; i++)
                {
                    newList.Add(item.Value[i]);
                }
            }
        }
        return newList; 
    }

    public static bool IsKeyInDict(Dictionary<Pair, List<int>> v_per_he, Pair pair)
    {
        foreach (KeyValuePair<Pair, List<int>> item in v_per_he)
        {
            if (item.Key.first == pair.first && item.Key.second == pair.second)
            {
                return true;
            }
        }
        return false;
    }

    public static Vector3 compute_vnormal(Vertex v)
    {
        Vector3 n = new Vector3(0.0f, 0.0f, 0.0f);
        HalfEdge he = v.halfEdge;

        do
        {
            // compute normal of the plane of the face
            Triangle f = he.t;
            Vector3 face_normal = compute_fnormal(f);
            // compute area of the triangular face
            float face_area = 0.5f * face_normal.magnitude;
            face_normal = face_normal.normalized;
            // accumulate onto our normal vector
            n += face_normal * face_area;
            if (he.oppositeEdge == null)
            {
                break;
            }
            he = he.oppositeEdge.nextEdge;
        } while (he != v.halfEdge);

        n = n.normalized;
        return n;
    }

    public static Vector3 compute_fnormal(Triangle f)
    {
        // Unity has a clockwise winding order
        Vertex v1 = f.halfEdge.src;
        Vertex v2 = f.halfEdge.nextEdge.src;
        Vertex v3 = f.halfEdge.nextEdge.nextEdge.src;

        Vector3 vec1 = v2.position - v1.position;
        Vector3 vec2 = v3.position - v1.position;
        Vector3 face_normal = Vector3.Cross(vec1, vec2);
        return face_normal;
    }

    public static Pair get_edge_key(int x, int y)
    {
        Pair p = new Pair(Math.Min(x, y), Math.Max(x, y));
        return p;
    }

    public static void addTriangleToBuffer(List<Vertex> vertices, Triangle t, List<int> triangles)
    {
        int idx1 = t.v1.index;
        int idx2 = t.v2.index;
        int idx3 = t.v3.index;

        if (!IsTriangleClockwise(idx1, idx2, idx3, vertices))
        {
            int tmp = idx1;
            idx1 = idx3;
            idx3 = tmp;
        }

        triangles.Add(idx1);
        triangles.Add(idx2);
        triangles.Add(idx3);
    }

    public static bool IsTriangleClockwise(int idx1, int idx2, int idx3, List<Vertex> hevs)
    {
        Vector3 e2 = hevs[idx3].position - hevs[idx2].position;
        Vector3 e1 = hevs[idx1].position - hevs[idx2].position;
        Vector3 r = Vector3.Cross(e1, e2);
        if (r.z < 0.0f)
        {
            return false;
        }
        return true;

    }

    public static void add_vertex(Vertex v, bool isBoundary, List<Vertex> hevs)
    {
        hevs.Add(v);
        hevs[hevs.Count - 1].isBoundary = isBoundary;
        hevs[hevs.Count - 1].index = hevs.Count - 1;
    }

    public static void add_face(int idx1, int idx2, int idx3, List<Vertex> hevs, List<Triangle> new_hefs)
    {
        if (!IsTriangleClockwise(idx1, idx2, idx3, hevs))
        {
            int tmp = idx1;
            idx1 = idx3;
            idx3 = tmp;
        }
        new_hefs.Add(new Triangle(hevs[idx1], hevs[idx2], hevs[idx3]));
        new_hefs[new_hefs.Count - 1].visited = false;
    }

    public static Vector3 centroidOfTriangle(Triangle hef)
    {
        Vertex a = hef.halfEdge.src;
        Vertex b = hef.halfEdge.nextEdge.src;
        Vertex c = hef.halfEdge.nextEdge.nextEdge.src;

        float x = (1.0f / 3.0f) * (a.position.x + b.position.x + c.position.x);
        float y = (1.0f / 3.0f) * (a.position.y + b.position.y + c.position.y);

        Vector3 v = new Vector3(x, y, 0.0f);
        return v;
    }

    public static int findIndexOfHEV(List<Vertex> hevs, float x_t, float y_t)
    {
        float tol = 0.00001f;
        for (int i=0; i < hevs.Count; i++)
        {
            float x = hevs[i].position.x;
            float y = hevs[i].position.y;
            if ((float) Math.Abs(x-x_t) < tol && (float) Math.Abs(y-y_t) < tol)
            {
                return i;
            }
        }
        return -1;
    }

    public static List<int> moduloSort(List<int> fanpoints)
    {
        List<int> newList = new List<int>();

        fanpoints.Sort();
        int rotate_point = 0;

        for (int k=0; k < fanpoints.Count-1; k++)
        {
            int diff = fanpoints[k + 1] - fanpoints[k];
            if (diff > 1)
            {
                rotate_point = k + 1;
                break;
            }
        }

        for (int i=rotate_point; i < fanpoints.Count; i++)
        {
            newList.Add(fanpoints[i]);
        }
        for (int i=0; i < rotate_point; i++)
        {
            newList.Add(fanpoints[i]);
        }

        return newList;
    }


    // Determine if triangle is TERMINAL, SLEEVE, or JUNCTION.
    public static void setTriangleType(Triangle triangle)
    {
        int numBoundaryEdges = 0;
        HalfEdge he = triangle.halfEdge;
        do
        {
            if (he.oppositeEdge == null)
            {
                numBoundaryEdges++;
            }
            he = he.nextEdge;
        } while (he != triangle.halfEdge);

        if (numBoundaryEdges == 0)
        {
            triangle.type = 2; // JUNCTION
        }
        else if (numBoundaryEdges == 1)
        {
            triangle.type = 1; // SLEEVE
        }
        else if (numBoundaryEdges == 2)
        {
            triangle.type = 0; // TERMINAL
        }

    }

    public static void labelTriangles(List<Triangle> triangles)
    {
        for (int i=0; i < triangles.Count; i++)
        {
            setTriangleType(triangles[i]);
        }
    } 


    /* Return true if (x,y) is outside the circle at (x_c, y_c) with radius r. */
    public static bool isOutsideCircle(float x, float y, float r, float x_c, float y_c)
    {
        if ((x-x_c)*(x-x_c) + (y-y_c)*(y-y_c) > r*r)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

}