/* Using Habrador https://www.habrador.com/tutorials/math/. */
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshOperations
{
    //Orient triangles so they have the correct orientation
    public static void OrientTrianglesClockwise(List<Triangle> triangles)
    {
        for (int i = 0; i < triangles.Count; i++)
        {
            Triangle tri = triangles[i];

            Vector2 v1 = new Vector2(tri.v1.position.x, tri.v1.position.z);
            Vector2 v2 = new Vector2(tri.v2.position.x, tri.v2.position.z);
            Vector2 v3 = new Vector2(tri.v3.position.x, tri.v3.position.z);

            if (!Geometry.IsTriangleOrientedClockwise(v1, v2, v3))
            {
                tri.ChangeOrientation();
            }
        }
    }
}