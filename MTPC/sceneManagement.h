#pragma once
#include <iostream>
#include <vector>
//#include <GL/glew.h>
//#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <string>
#include <fstream>
#include <map>
#include "morton code.h"
#include <opencv2/opencv.hpp>
//#include <atlimage.h>


double dmin(double p1, double p2, double p3);
double dmax(double p1, double p2, double p3);

class Vertex { //顶点坐标
public:
	double x, y, z;
	Vertex() {
		x = y = z = 0;
	}
	Vertex(double mx, double my, double mz):x(mx),y(my),z(mz){}
	Vertex operator+(Vertex v2) {
		Vertex ret;
		ret.x = x + v2.x;
		ret.y = y + v2.y;
		ret.z = z + v2.z;
		return ret;
	}
	Vertex operator/(double m) {
		Vertex ret;
		ret.x = x / m;
		ret.y = y / m;
		ret.z = z / m;
		return ret;
	}
	Vertex operator-(Vertex v2) {
		Vertex ret;
		ret.x = x - v2.x;
		ret.y = y - v2.y;
		ret.z = z - v2.z;
		return ret;
	}
	Vertex operator*(double t) {
		Vertex ret;
		ret.x = x * t;
		ret.y = y * t;
		ret.z = z * t;
		return ret;
	}
	Vertex negative()
	{
		Vertex ret;
		ret.x = -x, ret.y = -y, ret.z = -z;
		return ret;
	}
	//double operator*(glm::fvec3 v2) {
	//	double ret;
	//	ret = x * v2.x + y * v2.y + z * v2.z;
	//	return ret;
	//}
	double operator*(Vertex v2) {
		double ret;
		ret = x * v2.x + y * v2.y + z * v2.z;
		return ret;
	}
	Vertex cross(Vertex v2) {
		Vertex ret;
		ret.x = y * v2.z - v2.y * z;
		ret.y = v2.x * z - x * v2.z;
		ret.z = x * v2.y - v2.x * y;
		return ret;
	}
	double norm() {
		return sqrt(x*x + y * y + z * z);
	}
	void print() {
		std::cout << x << " " << y << " " << z << std::endl;
	}
	void normalize() {
		double dis = norm();
		x = x / dis, y = y / dis, z = z / dis;
	}
};

double distance(Vertex v1, Vertex v2);

//class VertexN {//法向量
//public:
//	glm::fvec3 vn;
//	double norm() {
//		return sqrt(pow(vn.x, 2) + pow(vn.y, 2) + pow(vn.z, 2));
//	}
//	void print() {
//		std::cout << "vn: " << vn.x << " " << vn.y << " " << vn.z << std::endl;
//	}
//};

class VertexT {//贴图坐标
public:
	double vtx, vty;
	void print() {
		std:cout << "vt: " << vtx << " " << vty << std::endl;
	}
};


class Face {//面对应的顶点索引
public:
	Vertex v1, v2, v3;
	Vertex vn1, vn2, vn3;
	VertexT vt1, vt2, vt3;
	Vertex norm;
	std::string material;
	unsigned int morton_code;
	void print();
	double calAera();
	void calNorm();
};

class Material {
public:
	std::string name;
	Vertex kd, ks;
	double Ns;
	double Ni;
	vector <Face> f; //存放所有对应该材质的三角形面片
	std::string map_Kd_filename;
	bool mapping_flag;
	cv::Mat img;
	int map_height, map_width;
	void readinMap() {
		std::cout << "Reading file " << map_Kd_filename << std::endl;
		
		img = cv::imread(map_Kd_filename); //read in file

		if (img.empty()) {
			std::cout << "Cannot read file: " << map_Kd_filename << std::endl;
		}
		map_height = img.rows, map_width = img.cols;
	}
	Material() {
		mapping_flag = false;
	}
	void print();
};

class Camera {
public:
	Vertex eye, look_at, up;
	double fovy;
	int width, height;
	void print();
};

class Light {
public:
	std::string name;
	Vertex radiance;
	void print();
};

class boundingBox {
public:
	double max_x, max_y, max_z, min_x, min_y, min_z; //boundingBox
	void print() {
		cout << min_x << " " << min_y << " " << min_z << " " << max_x << " " << max_y << " " << max_z << endl;
	}
};

class scene_data {
public:
	std::vector <Vertex> v;
	std::vector <Vertex> vn;
	std::vector <VertexT> vt;
	std::vector <Face> f;
	std::vector <Light> l;
	Camera camera;
	std::map<std::string, Material*> material_map;
	std::map<std::string, Light*> light_map;
	bool read_scene(std::string filename);
	bool read_obj(std::string filename);
	bool read_mtl(std::string filename);
	bool read_xml(std::string filename);
	~scene_data() {
		map <string, Material*>::iterator t = material_map.begin();
		while (t != material_map.end()) {
			delete t->second;
			t++;
		}
		map<string, Light*>::iterator t2 = light_map.begin();
		while (t2 != light_map.end()) {
			delete t2->second;
			t2++;
		}
	}
};

bool compare(Face t1, Face t2);

enum type {
	DIFFUSE, SPECULAR, TRANSMISSION
};

class Ray {
public:
	Vertex start;
	Vertex direction;
	enum type ray_type;
	Ray(Vertex s, Vertex d, type t) :start(s), direction(d), ray_type(t){}
	Ray(){}
};

class image {
public:
	double *img;
	double *avg_img;
	int width, height;
	image(int width, int height) {
		this->width = width, this->height = height;
		img = new double[width * height * 3];
		avg_img = new double[width * height * 3];
		memset(img, 0.0, sizeof(double)*width*height * 3);//initialize img
		memset(avg_img, 0.0, sizeof(double)*width*height * 3);//initialize img
	}
	~image() {
		delete[]img;
		delete[]avg_img;
	}
	int getIndex(int x, int y, int rgb) { //x colunmn y row
		return y * width * 3 + 3 * x + rgb;
	}
	void average(int core) {
		for (int i = core; i + core < width; i++) {
			for (int j = core; j + core < height; j++) {
				double current_r = 0, current_g = 0, current_b = 0;
				for (int k = -core; k <= core; k++) {
					for (int m = -core; m <= core; m++) {
						current_r += img[getIndex(i + k, j + m, 0)];
						current_g += img[getIndex(i + k, j + m, 1)];
						current_b += img[getIndex(i + k, j + m, 2)];
					}
				}
				current_r /= pow(2 * core + 1, 2), current_g /= pow(2 * core + 1, 2), current_b /= pow(2 * core + 1, 2);
				avg_img[getIndex(i, j, 0)] = current_r;
				avg_img[getIndex(i, j, 1)] = current_g;
				avg_img[getIndex(i, j, 2)] = current_b;
			}
		}
	}
};

//record the intersection information
class intersection {
public:
	Ray ray;  //ray = start_point + t * direction
	Vertex p; //intersection point
	Vertex pn;//intersection point normal
	double t; //ray parameters t
	Face f;   //intersection face
};

//judge if a ray intersects with a triangle, ret is the intersection point
bool intersect(Ray &ray, Face &triangle, Vertex &ret);

//judge if a ray intersects with a boungingbox
bool intersect(Ray &ray, boundingBox &b);
