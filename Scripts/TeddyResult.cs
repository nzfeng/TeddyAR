using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TeddyResult
{
    public List<Vector3> points;
    public List<Vector3> vertices;
    public List<int> triangles;

    public TeddyResult(List<Vector3> points, List<Vector3> vertices, List<int> triangles)
    {
        this.points = points;
        this.vertices = vertices;
        this.triangles = triangles;
    }
}