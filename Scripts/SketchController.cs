using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using GoogleARCore;
using GoogleARCoreInternal;

public class SketchController : MonoBehaviour {

    static Material lineMaterial;
    public Material material;
    public Vector3 spawnPos;
    public List<Vector3> points = new List<Vector3>();
    public Mesh mesh;
    public bool isSketching;

    private Material createMaterial(Color color)
    {
        return new Material("Shader \"Lines/Background\" { Properties { _Color (\"Main Color\", Color) = (" + color.r + "," + color.g + "," + color.b + "," + color.a + ") } SubShader { Pass { ZWrite on  Blend SrcAlpha OneMinusSrcAlpha Colormask RGBA Lighting Off Offset 1, 1 Color[_Color] }}}");
    }

    // Use this for initialization
    void Start () {

        //Mesh mesh = new Mesh();
        //if (!GetComponent<MeshFilter>() || !GetComponent<MeshRenderer>())
        //{
        //    gameObject.AddComponent<MeshFilter>();
        //    // gameObject.AddComponent<MeshRenderer>();
        //}
        //GetComponent<MeshFilter>().mesh = mesh;
        mesh = new Mesh();
        mesh.Clear();
        // GetComponent<Renderer>().material = material;

        // material = createMaterial(Color.red);
    }
	
	// Update is called once per frame
	void Update () {

        if (isSketching == false)
        {
            Debug.Log("false");
            isSketching = true;
            //Mesh mesh = GetComponent<MeshFilter>().mesh;
            mesh.Clear();
            //Shader shader = Shader.Find("Custom/mesh_shader");
            //material = new Material(shader);

            TeddyResult res = Teddy.inflate(points);
            mesh.vertices = res.vertices.ToArray();
            mesh.triangles = res.triangles.ToArray();
            // points = res.points;
            points.Clear();
        }
        DrawMesh();
    }

	static void CreateLineMaterial()
    {
        if (!lineMaterial)
        {
            // Unity has a built-in shader that is useful for drawing
            // simple colored things.
            Shader shader = Shader.Find("Hidden/Internal-Colored");
            lineMaterial = new Material(shader);
            lineMaterial.hideFlags = HideFlags.HideAndDontSave;
            // Turn on alpha blending
            lineMaterial.SetInt("_SrcBlend", (int)UnityEngine.Rendering.BlendMode.SrcAlpha);
            lineMaterial.SetInt("_DstBlend", (int)UnityEngine.Rendering.BlendMode.OneMinusSrcAlpha);
            // Turn backface culling off
            lineMaterial.SetInt("_Cull", (int)UnityEngine.Rendering.CullMode.Off);
            // Turn off depth writes
            lineMaterial.SetInt("_ZWrite", 0);
        }
    }

    // Draw sketched lines
    public void DrawLines()
    {

        CreateLineMaterial();
        // Apply the line material
        lineMaterial.SetPass(0);

        GL.PushMatrix();
        GL.LoadOrtho(); // for drawing in 2D

        // Draw lines
        int n = points.Count;
        // Debug.Log(n); // to confirm that OnPostRender is being called every frame

        for (int i = 0; i < n - 1; i++)
        {

            GL.Begin(GL.LINES);
            GL.Vertex3(points[i].x+0.5f, points[i].y+0.5f, points[i].z);
            GL.Vertex3(points[(i + 1)].x+0.5f, points[(i + 1)].y+0.5f, points[(i + 1)].z);
            GL.End();
        }

        GL.PopMatrix();
    }

    public void DrawMidLines()
    {
        CreateLineMaterial();
        // Apply the line material
        lineMaterial.SetPass(0);

        GL.PushMatrix();
        // GL.MultMatrix(transform.localToWorldMatrix);
        GL.LoadOrtho();

        // Draw lines
        int n = points.Count;

        for (int i = 0; i < n - 2; i += 3)
        {

            GL.Begin(GL.LINES);
            GL.Vertex3(points[i].x + 0.5f, points[i].y + 0.5f, points[i].z);
            GL.Vertex3(points[(i + 1)].x + 0.5f, points[(i + 1)].y + 0.5f, points[(i + 1)].z);
            GL.End();

            GL.Begin(GL.LINES);
            GL.Vertex3(points[i+1].x + 0.5f, points[i+1].y + 0.5f, points[i+1].z);
            GL.Vertex3(points[(i + 2)].x + 0.5f, points[(i + 2)].y + 0.5f, points[(i + 2)].z);
            GL.End();

            GL.Begin(GL.LINES);
            GL.Vertex3(points[i+2].x + 0.5f, points[i+2].y + 0.5f, points[i+2].z);
            GL.Vertex3(points[i].x + 0.5f, points[(i)].y + 0.5f, points[(i)].z);
            GL.End();
        }

        GL.PopMatrix();
    }

    public void DrawMesh()
    {

        //Mesh mesh = GetComponent<MeshFilter>().mesh;
        for (int i=0; i < mesh.vertices.Length; i++)
        {
            mesh.vertices[i] = transform.InverseTransformPoint(mesh.vertices[i]);
        }

        mesh.RecalculateNormals();
        Graphics.DrawMesh(mesh, new Vector3(spawnPos.x, spawnPos.y, spawnPos.z), Quaternion.LookRotation(transform.forward, transform.up), material, 0);
    }

    public void DrawTriangles()
    {
        //Mesh mesh = GetComponent<MeshFilter>().mesh;

        material.SetPass(0);

        GL.PushMatrix();
        // GL.LoadOrtho();
        GL.MultMatrix(transform.localToWorldMatrix);

        for (int i=0; i < mesh.triangles.Length; i += 3)
        {
            GL.Begin(GL.TRIANGLES);
            GL.Vertex3(mesh.vertices[mesh.triangles[i]].x, mesh.vertices[mesh.triangles[i]].y, mesh.vertices[mesh.triangles[i]].z+spawnPos.z);
            GL.Vertex3(mesh.vertices[mesh.triangles[i+1]].x, mesh.vertices[mesh.triangles[i+1]].y, mesh.vertices[mesh.triangles[i+1]].z +spawnPos.z);
            GL.Vertex3(mesh.vertices[mesh.triangles[i+2]].x, mesh.vertices[mesh.triangles[i+2]].y, mesh.vertices[mesh.triangles[i+2]].z + spawnPos.z);
            // GL.Vertex3(mesh.vertices[mesh.triangles[i]].x +  0.5f, mesh.vertices[mesh.triangles[i]].y + 0.5f, 0.0f);
            // GL.Vertex3(mesh.vertices[mesh.triangles[i + 1]].x +0.5f, mesh.vertices[mesh.triangles[i + 1]].y + 0.5f, 0.5f);
            // GL.Vertex3(mesh.vertices[mesh.triangles[i + 2]].x + 0.5f, mesh.vertices[mesh.triangles[i + 2]].y + 0.5f, 0.5f);
            GL.End();
        }
        GL.PopMatrix();

        // Orthographic
        //CreateLineMaterial();
        //lineMaterial.SetPass(0);
        //Mesh mesh = GetComponent<MeshFilter>().mesh;
        //GL.PushMatrix();
        //GL.LoadOrtho();

        //int n = points.Count;

        //for (int i = 0; i < mesh.triangles.Length-2; i += 3)
        //{

        //    GL.Begin(GL.LINES);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i]].x + 0.5f, mesh.vertices[mesh.triangles[i]].y + 0.5f, 0.0f);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i+1]].x + 0.5f, mesh.vertices[mesh.triangles[i+1]].y + 0.5f, 0.0f);
        //    GL.End();

        //    GL.Begin(GL.LINES);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i + 1]].x + 0.5f, mesh.vertices[mesh.triangles[i + 1]].y + 0.5f, 0.0f);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i + 2]].x + 0.5f, mesh.vertices[mesh.triangles[i + 2]].y + 0.5f, 0.0f);
        //    GL.End();

        //    GL.Begin(GL.LINES);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i + 2]].x + 0.5f, mesh.vertices[mesh.triangles[i + 2]].y + 0.5f, 0.0f);
        //    GL.Vertex3(mesh.vertices[mesh.triangles[i]].x + 0.5f, mesh.vertices[mesh.triangles[i]].y + 0.5f, 0.0f);
        //    GL.End();
        //}

        //GL.PopMatrix();
    }

	public void OnPostRender()
    {
        DrawLines();

        // DrawMidLines(); // draw the result of CDT + re-triangulating

        // DrawMesh();

        // DrawTriangles();

    }
}
