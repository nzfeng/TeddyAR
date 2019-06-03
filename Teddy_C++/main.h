

#include <algorithm>
#include <iterator>
#include "structs.h"
#include "quaternion.h"
#include "halfedge.h"
#include "poly2tri.h"
#include "../Eigen/Dense"

void setup_scene(void);

void reshape(int width, int height);

void display(void);

void init_lights(void);

void set_lights(void);

void draw_objects(void);

void draw_lines(void);

void draw_triangulation(void);

void draw_mesh(void);

void draw_axes(void);

void save_image(int i);

void resample(void);

Vertex screen_to_world(float x, float y);

void mouse_pressed(int button, int state, int x, int y);

void mouse_moved(int x, int y);

Quaternion get_current_rotation(void);

Quaternion compute_rotation_quaternion(const Eigen::Vector3f &p_curr, 
	const Eigen::Vector3f &p_start);

void key_pressed(unsigned char key, int x, int y);