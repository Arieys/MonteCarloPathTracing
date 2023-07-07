#pragma once
#include <string>
#include <vector>

#include "material.h"

#include <glm/glm.hpp>

#include <opencv2/opencv.hpp>
//顶点
struct Vertex {
	// 位置
	glm::vec3 Position;
	// 法向量
	glm::vec3 Normal;
	// 纹理坐标
	glm::vec2 TexCoords;
	// u向量
	glm::vec3 Tangent;
	// v向量
	glm::vec3 Bitangent;
};
//纹理
struct Texture {
	cv::Mat img;
	std::string type;
	std::string path;
};

class Face {//面对应的顶点索引
public:
	std::vector <Vertex> v;
	glm::vec3 normal;
	Material material;
	unsigned int morton_code;
	double calAera();
	void calNormal();
};

bool compare(Face t1, Face t2);
double dmin(double p1, double p2, double p3);
double dmax(double p1, double p2, double p3);
