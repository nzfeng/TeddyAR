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
#include <cstdlib>

#include "structs.h"
#include "quaternion.h"
#include "poly2tri.h"
#include "main.h"
#include "halfedge.h"
#include "helper.h"

#include "../Eigen/Dense"

#include <OpenGl/glu.h>

#ifdef __APPLE__
	#include <OpenGL/gl.h>
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
	#include <GL/gl.h>
#endif

using namespace std;


int xres = 600;
int yres = 600;
int idx = 0;
float delta_t = 0.05;
int mode = 0; // 0 = drawing, 1 = orbit
float drawPlaneZ = 0.0;
const float tol = 1e-10;

// the camera
camera cam;
// perspective parameters
frustum frust;
// viector of lights in this scene
vector<point_light> lights;

/* Parameters for ArcBall UI. */
Quaternion last_rotation, current_rotation; // should be initialized to 1
Eigen::Vector3f p_start, p_curr;

bool is_pressed = false;

vector<Vertex> vertex_buffer;
vector<Vec3f> normal_buffer;
vector<Vec3f> color_buffer;

vector<Vertex> drawnPoints;
vector<p2t::Point*> polyline;
vector<p2t::Triangle*> triangles;
p2t::CDT *cdt;

vector<Vertex> test;

float deg2rad(float angle) {
    return angle * M_PI / 180.0;
}

float rad2deg(float angle) {
    return angle * 180.0 / M_PI;
}

void setup_scene(void) {

	// camera position
	transvector position;
	position.x = 0;
	position.y = 0;
	position.z = 300;
	cam.position = position;
	// camera orientation
	transvector orientation;
	float norm = 1.0;
	orientation.x = 0 / norm;
	orientation.y = 1 / norm;
	orientation.z = 0 / norm;
	orientation.theta = 0.0;
	cam.orientation = orientation;
	// frustum
	frust.far = 800;
	frust.near = 200;
	frust.left = -200;
	frust.right = 200;
	frust.bottom = -200;
	frust.top = 200;

	// a single point light
	point_light light;
	light.position[0] = 0;
	light.position[1] = 0;
	light.position[2] = 50;
	light.position[3] = 1.0;
	light.rgb[0] = 1.0;
	light.rgb[1] = 1.0;
	light.rgb[2] = 1.0;
	light.k = 0.2;
	lights.push_back(light);
}

/* Initialize OpenGL to the states we want it to be in. */
void init(void) {

	setup_scene();

	/* Use Gouraud shading. */
	glShadeModel(GL_SMOOTH);
	/* Enable backface cullling. */
	glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    /* Use depth buffering when rendering. */
    glEnable(GL_DEPTH_TEST);
    /* Automatically normalize normal vectors before passing them into 
     * normal arrays. 
     */
    glEnable(GL_NORMALIZE);
    /* Enable vertex and normal arrays. */
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    /* Set up projection and modelview matrices (akin to our perspective 
     * projection and world-space-to-camera-space matrices.)
     */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // glFrustum(frust.left, frust.right, frust.bottom, frust.top, 
    // 	frust.near, frust.far);
    // gluPerspective(50.0, 1.0, 3.0, 500);
    glOrtho(-50, 50, -50, 50, 3.0, 500);

    glMatrixMode(GL_MODELVIEW);

    /* Set up lights. */
    // init_lights();

}


/* Called when the window gets resized. */
void reshape(int width, int height) {

	/* Prevent the window from shrinking to 0. */
	height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;

    /* 'glViewport' determines how to convert from NDC to screen coordinates 
     * given the dimensions of the window. We let (0,0) represent the lower
     * left corner of the window.
     */
    glViewport(0, 0, width, height);

    /* Re-render. */
    glutPostRedisplay();

}

/* Handle processing of points in world and camera space. */
void display(void) {

	/* Reset the color and depth buffers. */
	glClearColor(1.0, 1.0, 1.0, 0.0); // set background color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* Reset modelview matrix. */
	glLoadIdentity();
    /* Create matrix representing the inverse rotation of the camera. */
    glRotatef(-rad2deg(cam.orientation.theta), cam.orientation.x, cam.orientation.y, cam.orientation.z);
    /* Create matrix representing the inverse translation of the camera. */
    glTranslatef(-cam.position.x, -cam.position.y, -cam.position.z);
    /* Create rotation matrices to represent rotations caused by the mouse. */
	glMultMatrixf(get_current_rotation().toMatrix().data());
    /* Set up lights in their specified positions. */
    set_lights();
    /* Specify the points and faces that we want drawn. */
    // draw_axes();
    /* draw sketched lines */
    // draw_lines();
    /* draw mesh */
    // draw_mesh();
    draw_triangulation();
    /* Swap active and hidden buffers. */
    glutSwapBuffers();
}

/* Enable lights and set their colors. */
void init_lights(void) {

	/* Enable lighting calculations during rendering. */
	glEnable(GL_LIGHTING);

	/* Associate each point light with one of OpenGL's built-in lights. */
	int num_lights = lights.size();
	for(int i = 0; i < num_lights; i++) {

		int light_id = GL_LIGHT0 + i;
		glEnable(light_id);
		/* Specify colors for the ambient, diffuse, and specular components 
		 * of the light, and the attenutation. 
		 */
        glLightfv(light_id, GL_AMBIENT, lights[i].rgb);
        glLightfv(light_id, GL_DIFFUSE, lights[i].rgb);
        glLightfv(light_id, GL_SPECULAR, lights[i].rgb);
        glLightf(light_id, GL_QUADRATIC_ATTENUATION, lights[i].k);
	}

}

/* Position the lights. */
void set_lights(void) {

	int num_lights = lights.size();
    
    for(int i = 0; i < num_lights; i++)
    {
        int light_id = GL_LIGHT0 + i;
        /* Apply the current Modelview Matrix to the given light position. */
        glLightfv(light_id, GL_POSITION, lights[i].position);
    }
}

void draw_lines() {

	if (drawnPoints.size() < 2) {
		return;
	}

	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);

    int n = drawnPoints.size();
    
    for (int i=0; i < n; i++) {
    	Vertex p1 = drawnPoints[i];
    	Vertex p2 = drawnPoints[(i+1)%n];
    	glBegin(GL_LINES);
	    glColor4f(0.0, 0.0, 0.0, 1.0);
	    glVertex3f(p1.x, p1.y, drawPlaneZ);
	    glVertex3f(p2.x, p2.y, drawPlaneZ);
	    glEnd();
    }
}


void draw_axes() {

	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);
    // positive x-axis
    glBegin(GL_LINES);
    glColor4f(1.0, 0.0, 0.0, 1.0);
    glVertex3f(0, 0, 0);
    glVertex3f(1000, 0, 0);
    glEnd();
    // y-axis
    glBegin(GL_LINES);
    glColor4f(0.0, 1.0, 0.0, 1.0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1000, 0);
    glEnd();
    // z-axis
    glBegin(GL_LINES);
    glColor4f(0.0, 0.0, 1.0, 1.0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1000);
    glEnd();
}

void draw_triangulation(void) {

	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);

	// for (int i=0; i < triangles.size(); i++) {
	// 	for (int j=0; j < 3; j++) {
	// 		glBegin(GL_LINES);
	// 	    glColor4f(1.0, 0.0, 0.0, 1.0);
	// 	    glVertex3f(triangles[i]->GetPoint(j)->x, triangles[i]->GetPoint(j)->y, drawPlaneZ);
	// 	    glVertex3f(triangles[i]->GetPoint((j+1)%3)->x, triangles[i]->GetPoint((j+1)%3)->y, drawPlaneZ);
	// 	    glEnd();
	// 	}
	// }

	glLineWidth(3.0);
	for (int i=0; i < test.size(); i += 3) {
		for (int j=0; j < 3; j++) {
			glBegin(GL_LINES);
		    glColor4f(0.0, 0.0, 1.0, 1.0);
		    glVertex3f(test[i+j].x, test[i+j].y, test[i+j].z);
		    glVertex3f(test[i+(j+1)%3].x, test[i+(j+1)%3].y, test[i+(j+1)%3].z);
		    glEnd();
		}
	}
}


/* Render mesh. */
void draw_mesh(void) {

	float ambient[3] = {0.5, 0.0, 0.0};
	float diffuse[3] = {0.5, 0.0, 0.0};
	float specular[3] = {0.5, 0.0, 0.0};
	float shininess = 0.2;

	// Generate color buffer.
	for (int i=0; i < vertex_buffer.size(); i++) {
		Vec3f c; 
		c.x = 0.0; c.y = 0.7; c.z = 0.7;
		color_buffer.push_back(c);
	}

	glPushMatrix();

	/* Specify object material properties. */
	glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);

    /* Specify object vertex and normals. */
    glVertexPointer(3, GL_FLOAT, 0, &vertex_buffer[0]);
    glNormalPointer(GL_FLOAT, 0, &normal_buffer[0]);
    glColorPointer(3, GL_FLOAT, 0, &color_buffer[0]);

    int buffer_size = vertex_buffer.size();
    glDrawArrays(GL_TRIANGLES, 0, buffer_size);

    /* Pop this transformation off the stack to make way for the next. */
	glPopMatrix();
}

/* Resample drawn points to remove redundant points and short edges. */
void resample() {

	int n = drawnPoints.size();
	float length = 5.0;
	polyline.push_back(new p2t::Point(drawnPoints[0].x, drawnPoints[0].y));
	int i=0;

	while (i < n-1) {

		Vertex p2 = drawnPoints[(i+1)%n];
		float x = polyline[polyline.size()-1]->x;
		float y = polyline[polyline.size()-1]->y;
		float norm = sqrt((x-p2.x)*(x-p2.x) + (y-p2.y)*(y-p2.y));

		if (norm > length) {
			float a = length / norm;
			polyline.push_back(new p2t::Point(x + a*(p2.x-x), y + a*(p2.y-y)));
		}
		else {
			i++;
		}
	}
}

// Save current screenshot.
void save_image(int i) {

	const char filenamepattern[] = "frames/frame%s%d.ppm";
	char filename[1024];
	FILE *f;
	string s = "";
	if (i < 10) {s = "0";}
	sprintf(filename, filenamepattern, s.c_str(), i);
	f = fopen(filename, "wb");
	//fflush(stdout);

	// Read pixels
	GLint V[4];
    /* V has 4 values; the x and y window coords of the viewport, followed
     * by its width and height. x,y represents the lower left corner of of the viewport rectangle, in pixels.
     */
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];

	char *buf = new char[width*height*3];
    /* GL_PACK_ALIGNMENT specifies the alignment requirements for the start of each pixel
     * row in memory. Here, 1 signifies byte-alignment.
     */
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
    /* glReadPixels(x, y, width, height, format, type, data) reads a block of pixels from the framebuffer.
     * x, y specify the window coords of the first pixel that is read. This location is the lower left corner of a rectangular block of pixels.
     * width, height specify the dimensions of the pixel rectangle.
     * 'format' specifies the format of the pixel data. Options: GL_ALPHA, GL_RGB, GL_RGBA
     * 'type' specifies the data type of the pixel data. Options: GL_UNSIGNED_BYTE, GL_UNSIGNED_SHORT_5_6_5, GL_UNSIGNED_SHORT_4_4_4_4, GL_UNSIGNED_SHORT_5_5_5_1
     * 'data' returns the pixel data.
     */
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;
}

Vertex screen_to_world(int x, int y) {

	GLint viewport[4];  // corner, width, height
	glGetIntegerv(GL_VIEWPORT, viewport); 
	GLdouble modelview[16]; 
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	GLdouble projection[16]; 
	glGetDoublev(GL_PROJECTION_MATRIX, projection);  
	float winX = (float) x;
	float winY = (float) y;
	float winZ = 0.0;
	winY = (float)viewport[3] - winY; 
	//glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	GLdouble posX, posY, posZ;    

	gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
	//Vertex v; v.x = posX; v.y = posY; v.z = 0.0;
	Vertex v; v.x = (winX-viewport[2]/2)*(100.0/viewport[2]); v.y = (winY-(viewport[3]/2))*(100.0/viewport[3]); v.z = 0.0;

	return v;
}

/* Respond to mouse clicks and releases. */
void mouse_pressed(int button, int state, int x, int y) {

	/* if the left-mouse button was clicked down */
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    	is_pressed = true;
    	if (mode == 1) {
	        // Get viewport width and height and convert from screen to NDC.
	        GLfloat V[4];
	        glGetFloatv(GL_VIEWPORT, V);
	        float x_NDC = (float)x / (V[2] - 1.0) * 2.0 - 1.0; 
	        float y_NDC = (float)y / (V[3] - 1.0) * 2.0 - 1.0;
	        p_start[0] = x_NDC; 
	        p_start[1] = -y_NDC; 
	        p_start[2] = sqrt(max(0.0, 1.0 - x_NDC*x_NDC - y_NDC*y_NDC));
    	}
    	// drawing mode
        else if (mode == 0) {

    		delete cdt;
    		for (int i=0; i < polyline.size(); i++) {
    			cerr << "clear: " << i << endl;
    			delete polyline[i];
    		}
    		polyline.clear();

        	// push the initial point
        	Vertex p = screen_to_world(x, y);
        	drawnPoints.push_back(p);
        	
        }
    }
    /* if the left-mouse button was released */
    else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
    	is_pressed = false;
    	if (mode == 1) {
    		last_rotation = current_rotation * last_rotation;
    		current_rotation = Quaternion();
    	}
    	// drawing mode
    	else if (mode == 0) {
    		// Perform CDT.
    		resample();
    		cdt = new p2t::CDT(polyline);
    		cdt->Triangulate();
    		triangles = cdt->GetTriangles();

    		inflate(polyline, triangles, test, normal_buffer);

  			glutPostRedisplay();
    	}
    }
}

/* Respond to mouse movement (for ArcBall UI). */
void mouse_moved(int x, int y) {
	// orbit mode
	if (is_pressed && mode == 1) {
		// Get viewport width and height and convert from screen to NDC.
        GLfloat V[4];
        glGetFloatv(GL_VIEWPORT, V);
        float x_NDC = (float)x / (V[2] - 1.0) * 2.0 - 1.0; 
        float y_NDC = (float)y / (V[3] - 1.0) * 2.0 - 1.0;
		p_curr[0] = x_NDC;
		p_curr[1] = -y_NDC;
		p_curr[2] = sqrt(max(0.0, 1.0 - x_NDC*x_NDC - y_NDC*y_NDC));
		current_rotation = compute_rotation_quaternion(p_curr, p_start);
        glutPostRedisplay();
	}
	// drawing mode
	else if (is_pressed && mode == 0) {
		Vertex p = screen_to_world(x, y);
		drawnPoints.push_back(p);

		glutPostRedisplay();
	}
}

Quaternion get_current_rotation(void) {
	Quaternion q = current_rotation * last_rotation;
	q = q / q.norm();
	return q;
}

/* Return the quaternion that represents the given ArcBall rotation. */
Quaternion compute_rotation_quaternion(const Eigen::Vector3f &p_curr, 
	const Eigen::Vector3f &p_start) {

	Quaternion q;
	float theta = acos(min(1.0f, p_start.dot(p_curr) / 
		(p_start.norm()*p_curr.norm())));

	Eigen::Vector3f u = p_start.cross(p_curr);
	u = u / u.norm();

	q.s = cos(theta/2.0);
	if (theta == 0.0) { q.s = 1.0; }
	q.v << u[0] * sin(theta/2.0), u[1] * sin(theta/2.0), u[2] * sin(theta/2.0);

	return q;
}

/* Respond to keyboard input. */
void key_pressed(unsigned char key, int x, int y) {

	switch(key) {
		/* If 'q' is pressed, quit the program. */
		case 'q':
			exit(0);
			break;
        /* Switch between drawing and panning modes. */
        case 'f':
        	mode = (mode + 1) % 2;
        	cerr << "Mode: " << mode << endl;
            break;
        /* Clear drawing data. */
        case 'c':
        	break;
		/* If spacebar is pressed, reset to original view. */
		case ' ':
			/* Reset modelview matrix. */
			last_rotation = Quaternion();
	    	glutPostRedisplay();
			break;
		default:
			break;
	}
}


int main (int argc, char **argv) {
	if (argc != 1) {
		return 1;
	}
	else {

		/* Intialize the GLUT (Graphics Library Utility Toolkit) library. In 
		 * this case, we don't need to have GLUT process any command-line 
		 * arugments, so what we pass in doesn't really matter.
		 */
		glutInit(&argc, argv);
		/* Tell OpenGL that we need a RGB pixel buffer and a depth buffer,
		 * and enable double-buffering.
	     */
	    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	    /* Create an xres-by-yres window. */
	    glutInitWindowSize(xres, yres);
	    /* Tell OpenGL to set the program window at (0,0), the top-left corner 
	     * of the computer screen.
	     */
	    glutInitWindowPosition(0, 0);
	    /* Name the window. */
	    glutCreateWindow("Teddy");

	    /* Set up the program. */
	    init();
	    /* Specify to OpenGL our display function. */
	    glutDisplayFunc(display);
	    /* Specify to OpenGL our reshape function. */
	    glutReshapeFunc(reshape);
	    /* Specify to OpenGL our function for handling mouse presses. */
	    glutMouseFunc(mouse_pressed);
	    /* Specify to OpenGL our function for handling mouse movement. */
	    glutMotionFunc(mouse_moved);
	    /* Specify to OpenGL our function for handling key presses. */
	    glutKeyboardFunc(key_pressed);
	    /* Start the "event processing loop". */
	    glutMainLoop();

	}

	return 0;
}
