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
    public bool isSketching;

    // Use this for initialization
    void Start () {

        Mesh mesh = new Mesh();
        if (!GetComponent<MeshFilter>() || !GetComponent<MeshRenderer>())
        {
            gameObject.AddComponent<MeshFilter>();
            gameObject.AddComponent<MeshRenderer>();
        }
        GetComponent<MeshFilter>().mesh = mesh;
        mesh.Clear();
        GetComponent<Renderer>().material = material;
    }
	
	// Update is called once per frame
	void Update () {

        if (isSketching == false)
        {
            Debug.Log("false");
            isSketching = true;
            Mesh mesh = GetComponent<MeshFilter>().mesh;
            Teddy.inflate(points, mesh.vertices, mesh.normals, mesh.triangles);
            points.Clear();
        }
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

    public void DrawLines()
    {

        CreateLineMaterial();
        // Apply the line material
        lineMaterial.SetPass(0);

        GL.PushMatrix();
        GL.LoadOrtho(); // for drawing in 2D

        // Draw lines
        int n = points.Count;
        // Debug.Log(n); to confirm that OnPostRender is being called every frame

        for (int i = 0; i < n - 1; i++)
        {

            GL.Begin(GL.LINES);
            GL.Vertex3(points[i].x+Screen.width/2.0f/Screen.width, points[i].y+Screen.height/2.0f/Screen.height, points[i].z);
            GL.Vertex3(points[(i + 1)].x+Screen.width / 2.0f / Screen.width, points[(i + 1)].y + Screen.height / 2.0f / Screen.height, points[(i + 1)].z);
            GL.End();
        }

        GL.PopMatrix();
    }

    public void DrawMesh(Vector3 spawnPos)
    { 
    
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        // Graphics.DrawMesh(mesh, spawnPos, Quaternion.identity, material, 0);
    }

	public void OnPostRender()
    {
        DrawLines();

        DrawMesh(spawnPos);

    }
}
