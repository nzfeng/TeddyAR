# TeddyAR

This is an augmented reality (AR) app for Android implementing *Teddy* ([Igarashi 1999](https://www.cs.toronto.edu/~jacobson/seminar/igarashi-et-al-1999.pdf)), an inflation-based sketch-to-mesh algorithm.

![](heart_demo.gif)

![](swan_demo.gif)

Constrained Delaunay Triangulation (CDT) is currently implemented using [poly2tri](https://github.com/MaulingMonkey/poly2tri-cs), since it was fast enough to run in real-time. Specifically, the approach of [Sloan 1992](https://www.newcastle.edu.au/__data/assets/pdf_file/0019/22519/23_A-fast-algortithm-for-generating-constrained-Delaunay-triangulations.pdf) (implemented here using tutorials from habrador.com) was much too slow.

"Teddy_C++" contains a C++ implementation which uses OpenGL. This is useful if you want to see Teddy in action with minimal setup. Just run "make" in the directory of "Teddy_C++", and run "./main". Draw in 2D using the cursor. The mesh will be generated automatically, after which you can drag using the cursor to rotate the camera and inspect the mesh. (However, the C++ program does not check for bad input and will segfault if the sketch is self-intersecting. The C# version does check for bad input, and simply has the user re-draw.) You can also hit "c" to toggle the results of the CDT, or "p" to toggle the results of the re-triangulation procedure.

"Teddy" contains the materials used for the Android app, which was developed in Unity (with scripts in C#.)

There are currently some issues with speed (takes ~10 - 15 seconds to do sketch-to-mesh conversion on a Samsung Galaxy S9, mesh manipulations like translation and enlarge/shrink cause freezing), and a minor issue with mesh appearance. This project was intended as a proof of concept and as an exercise in AR, so the user interface is really only intuitive to the creator.

This project may be useful as an example of:
* using Google ARCore + Unity
* drawing in 2D screen coordinates
* creating, drawing, and manipulating meshes
* using plane detection
* C# programming
* an implementation of the Teddy drawing system
