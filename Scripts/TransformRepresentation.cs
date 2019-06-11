/* Using Habrador https://www.habrador.com/tutorials/math/. */
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TransformRepresentation
{
    //From triangle where each triangle has one vertex to half edge
    public static List<HalfEdge> TransformFromTriangleToHalfEdge(List<Triangle> triangles)
    {
        //Make sure the triangles have the same orientation
        MeshOperations.OrientTrianglesClockwise(triangles);

        //First create a list with all possible half-edges
        List<HalfEdge> halfEdges = new List<HalfEdge>(triangles.Count * 3);

        for (int i = 0; i < triangles.Count; i++)
        {
            Triangle t = triangles[i];

            HalfEdge he1 = new HalfEdge(t.v1);
            HalfEdge he2 = new HalfEdge(t.v2);
            HalfEdge he3 = new HalfEdge(t.v3);

            he1.nextEdge = he2;
            he2.nextEdge = he3;
            he3.nextEdge = he1;

            he1.prevEdge = he3;
            he2.prevEdge = he1;
            he3.prevEdge = he2;

            //The vertex needs to know of an edge going from it
            he1.v.halfEdge = he2;
            he2.v.halfEdge = he3;
            he3.v.halfEdge = he1;

            he1.src = t.v3;
            he2.src = t.v1;
            he3.src = t.v2;

            //The face the half-edge is connected to
            t.halfEdge = he1;

            he1.t = t;
            he2.t = t;
            he3.t = t;

            //Add the half-edges to the list
            halfEdges.Add(he1);
            halfEdges.Add(he2);
            halfEdges.Add(he3);
        }

        //Find the half-edges going in the opposite direction
        for (int i = 0; i < halfEdges.Count; i++)
        {
            HalfEdge he = halfEdges[i];

            Vertex goingToVertex = he.v;
            Vertex goingFromVertex = he.prevEdge.v;

            for (int j = 0; j < halfEdges.Count; j++)
            {
                //Dont compare with itself
                if (i == j)
                {
                    continue;
                }

                HalfEdge heOpposite = halfEdges[j];

                //Is this edge going between the vertices in the opposite direction
                if (goingFromVertex.position == heOpposite.v.position && goingToVertex.position == heOpposite.prevEdge.v.position)
                {
                    he.oppositeEdge = heOpposite;

                    break;
                }
            }
        }

        return halfEdges;
    }
}

