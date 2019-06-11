/* Using Habrador https://www.habrador.com/tutorials/math/. */
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

//Sort the points along one axis. The first 3 points form a triangle. Consider the next point and connect it with all
//previously connected points which are visible to the point. An edge is visible if the center of the edge is visible to the point.
public class IncrementalTriangulation
{
    public static List<Triangle> TriangulatePoints(List<Vertex> points)
    {
        List<Triangle> triangles = new List<Triangle>();

        //Sort the points along x-axis
        //OrderBy is always soring in ascending order - use OrderByDescending to get in the other order
        points = points.OrderBy(n => n.position.x).ToList();
        Debug.Log("IT: order points");

        //The first 3 vertices are always forming a triangle
        Debug.Log("IT # of vertices: " + points.Count);
        Triangle newTriangle = new Triangle(points[0].position, points[1].position, points[2].position);
        Debug.Log("IT: newTriangle");

        triangles.Add(newTriangle);
        Debug.Log("IT: add new triangle");

        //All edges that form the triangles, so we have something to test against
        List<Edge> edges = new List<Edge>();
        Debug.Log("IT: instantiate edges");

        edges.Add(new Edge(newTriangle.v1, newTriangle.v2));
        edges.Add(new Edge(newTriangle.v2, newTriangle.v3));
        edges.Add(new Edge(newTriangle.v3, newTriangle.v1));
        Debug.Log("IT: add first triangle");

        //Add the other triangles one by one
        //Starts at 3 because we have already added 0,1,2
        for (int i = 3; i < points.Count; i++)
        {
            Vector3 currentPoint = points[i].position;

            //The edges we add this loop or we will get stuck in an endless loop
            List<Edge> newEdges = new List<Edge>();

            //Is this edge visible? We only need to check if the midpoint of the edge is visible 
            for (int j = 0; j < edges.Count; j++)
            {
                Edge currentEdge = edges[j];

                Vector3 midPoint = (currentEdge.v1.position + currentEdge.v2.position) / 2f;

                Edge edgeToMidpoint = new Edge(currentPoint, midPoint);

                //Check if this line is intersecting
                bool canSeeEdge = true;

                for (int k = 0; k < edges.Count; k++)
                {
                    //Dont compare the edge with itself
                    if (k == j)
                    {
                        continue;
                    }

                    if (AreEdgesIntersecting(edgeToMidpoint, edges[k]))
                    {
                        canSeeEdge = false;

                        break;
                    }
                }

                //This is a valid triangle
                if (canSeeEdge)
                {
                    Edge edgeToPoint1 = new Edge(currentEdge.v1, new Vertex(currentPoint));
                    Edge edgeToPoint2 = new Edge(currentEdge.v2, new Vertex(currentPoint));

                    newEdges.Add(edgeToPoint1);
                    newEdges.Add(edgeToPoint2);

                    Triangle newTri = new Triangle(edgeToPoint1.v1, edgeToPoint1.v2, edgeToPoint2.v1);

                    triangles.Add(newTri);
                }
            }


            for (int j = 0; j < newEdges.Count; j++)
            {
                edges.Add(newEdges[j]);
            }
        }
        Debug.Log("IT: triangles added");


        return triangles;
    }



    private static bool AreEdgesIntersecting(Edge edge1, Edge edge2)
    {
        Vector2 l1_p1 = new Vector2(edge1.v1.position.x, edge1.v1.position.z);
        Vector2 l1_p2 = new Vector2(edge1.v2.position.x, edge1.v2.position.z);
        
        Vector2 l2_p1 = new Vector2(edge2.v1.position.x, edge2.v1.position.z);
        Vector2 l2_p2 = new Vector2(edge2.v2.position.x, edge2.v2.position.z);

        bool isIntersecting = Intersections.AreLinesIntersecting(l1_p1, l1_p2, l2_p1, l2_p2, true);

        return isIntersecting;
    }
}