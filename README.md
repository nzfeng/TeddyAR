# CS174

This is an augmented reality (AR) app for Android implementing *Teddy* ([Igarashi 1999](https://www.cs.toronto.edu/~jacobson/seminar/igarashi-et-al-1999.pdf)), an inflation-based sketch-to-mesh algorithm.

Constrained Delaunay Triangulation (CDT) is currently implemented using [poly2tri](https://github.com/MaulingMonkey/poly2tri-cs), since it was fast enough to run in real-time. Specifically, the approach of [Sloan 1992](https://www.newcastle.edu.au/__data/assets/pdf_file/0019/22519/23_A-fast-algortithm-for-generating-constrained-Delaunay-triangulations.pdf) was way too slow.

"Teddy_C++" contains a C++ implementation which uses OpenGL.

"Teddy" contains the materials used for the Android app, which was developed in Unity (with scripts in C#.)

There are currently some issues with speed (takes ~10 - 15 seconds to do sketch-to-mesh conversion on a Samsung Galaxy S9, mesh manipulations like translation and enlarge/shrink cause freezing), and a minor issue with mesh appearance. This project was intended as a proof of concept and as an exercise in AR, so the user interface is only intuitive to the creator.

This project may be useful as an example of:
* using Google ARCore + Unity
* drawing in 2D screen coordinates
* creating, drawing, and manipulating meshes
* using plane detection
* C# programming
* an implementation of the Teddy drawing system
