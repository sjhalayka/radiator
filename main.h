#ifndef main_H
#define main_H

#include "uv_camera.h"
#include "custom_math.h"

using custom_math::vector_3;
using custom_math::vector_4;
using custom_math::line_segment_3;



#include <cstdlib>
#include <GL/glut.h>       //GLUT Library

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;







size_t n = 1000;// 100000;
double emitter_distance = 10;
    
double dimension = 2.1;

const double G = 6.6743e-11;
const double pi = 4.0 * atan(1.0);
const double c_meters = 299792458;

vector<vector_3> ellipse_points;



void idle_func(void);
void init_opengl(const int& width, const int& height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);

void render_string(int x, const int y, void* font, const string& text);
void draw_objects(void);

custom_math::vector_3 background_colour(0.0f, 0.0f, 0.0f);
custom_math::vector_3 control_list_colour(1.0f, 1.0f, 1.0f);

bool draw_axis = true;
bool draw_control_list = true;

uv_camera main_camera;


GLint win_id = 0;
GLint win_x = 800, win_y = 600;
double camera_w = 10.0;

double camera_fov = 45.0f;
double camera_x_transform = 0;
double camera_y_transform = 0;
double u_spacer = 0.01f;
double v_spacer = 0.5f*u_spacer;
double w_spacer = camera_w*0.01f;
double camera_near = 1.0f;
double camera_far = 10000.0;

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;

//
//vector_4 iEllipsoid(vector_3 ro, vector_3 rd, vector_3 r)
//{
//    vector_3 r2 = r * r;
//    double a = rd.dot(rd / r2);
//    double b = ro.dot(rd / r2);
//    double c = ro.dot(ro / r2);
//    double h = b * b - a * (c - 1.0);
//    if (h < 0.0) return vector_4(-1.0, -1.0, -1.0, -1.0);
//    double t = (-b - sqrt(h)) / a;
//    vector_3 n = (ro + rd*t) / r2;
//    n.normalize();
//
//    return vector_4(t, n.x, n.y, n.z);
//}


//
//
//// https://paulbourke.net/geometry/circlesphere/raysphere.c
//int RaySphere(vector_3 p1, vector_3 p2, vector_3 sc, double r, double* mu1, double* mu2)
//{
//    double a, b, c;
//    double bb4ac;
//    vector_3 dp;
//
//    dp.x = p2.x - p1.x;
//    dp.y = p2.y - p1.y;
//    dp.z = p2.z - p1.z;
//    a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
//    b = 2 * (dp.x * (p1.x - sc.x) + dp.y * (p1.y - sc.y) + dp.z * (p1.z - sc.z));
//    c = sc.x * sc.x + sc.y * sc.y + sc.z * sc.z;
//    c += p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;
//    c -= 2 * (sc.x * p1.x + sc.y * p1.y + sc.z * p1.z);
//    c -= r * r;
//    bb4ac = b * b - 4 * a * c;
//    if (fabs(a) < 1e-5 || bb4ac < 0) {
//        *mu1 = 0;
//        *mu2 = 0;
//        return(FALSE);
//    }
//
//    *mu1 = (-b + sqrt(bb4ac)) / (2 * a);
//    *mu2 = (-b - sqrt(bb4ac)) / (2 * a);
//
//    return(TRUE);
//}






//
//void get_line_segments(const vector_3 sphere_location,
//	const double sphere_radius,
//	const double dimension)
//{
//	srand(0);
//
//	const double disk_like = 3 - dimension;
//
//	if (1)//true == redo_line_segments)
//	{
//		threeD_line_segments.clear();
//
//		for (size_t i = 0; i < n; i++)
//		{
//			line_segment_3 ls;
//			ls.start = threeD_oscillators[i];
//
//			for (size_t j = 0; j < n; j++)
//			{
//				if (i == j)
//					continue;
//
//				ls.end = threeD_oscillators[j];
//
//				line_segment_3 ls_;
//
//				ls_.start = ls.start;
//				ls_.end = ls.start + (ls.start - ls.end).normalize() * 1e30;// end;// +(ls.end - ls.start).normalize() * 10.0f;
//
//				threeD_line_segments.push_back(ls_);
//			}
//		}
//	}
//
//	threeD_line_segments_intersected.clear();
//
//	for (size_t i = 0; i < threeD_line_segments.size(); i++)
//	{
//		const vector_3 dir = (threeD_line_segments[i].end - threeD_line_segments[i].start).normalize();
//
//		if (dir.dot(sphere_location) > 0)
//		{
//			double mu1 = 0, mu2 = 0;
//
//			if (RaySphere(threeD_line_segments[i].start, threeD_line_segments[i].end, sphere_location, sphere_radius, &mu1, &mu2))
//			{
//				line_segment_3 ls_;
//				ls_.start = threeD_line_segments[i].start;
//				ls_.end = threeD_line_segments[i].start + threeD_line_segments[i].end * mu2;
//
//				threeD_line_segments_intersected.push_back(ls_);
//			}
//		}
//	}
//
//	//	cout << static_cast<float>(threeD_line_segments_intersected.size()) / static_cast<float>(threeD_line_segments.size()) << ", " << endl;
//
//	vector<vector_3> vectors;
//
//	for (size_t i = 0; i < threeD_line_segments_intersected.size(); i++)
//	{
//		vector_3 v = (threeD_line_segments_intersected[i].end - threeD_line_segments_intersected[i].start);
//		v.normalize();
//		vectors.push_back(v);
//	}
//
//	double parallelity = 0;
//	size_t count = 0;
//
//	for (size_t i = 0; i < vectors.size() - 1; i++)
//	{
//		for (size_t j = (i + 1); j < vectors.size(); j++)
//		{
//			const double d = vectors[i].dot(vectors[j]);
//
//			parallelity += d;
//			count++;
//		}
//	}
//
//	//cout << count << " " << vectors.size() * (vectors.size() - 1) / 2.0 << endl;
//
//	parallelity /= count;
//
//	parallelity = abs(parallelity);
//
//	cout << "parallelity: " << parallelity << " " << parallelity / sphere_location.length() << endl;
//
//
//	double avg_strength = pow(c_meters, disk_like);
//
//	double a = 1.0 / pow(sphere_location.x, 2.0);
//	// static_cast<double>(threeD_line_segments_intersected.size()) / static_cast<double>(threeD_line_segments.size()) 
//	//
//	double g = (1.0 - parallelity);// / pow(sphere_location.x, 2.0);
//
//	//cout << g << endl;
//	//cout << a << endl;
//	//cout << g / a << endl;
//	//cout << endl;
//}
//


#endif
